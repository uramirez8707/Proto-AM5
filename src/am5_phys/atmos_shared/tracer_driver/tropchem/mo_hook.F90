      module MOZ_HOOK_MOD

      use diag_manager_mod, only : register_diag_field, send_data
      use atmos_cmip_diag_mod, only : register_cmip_diag_field_2d
      use time_manager_mod, only : time_type
      use constants_mod,    only : PI
      use fms_mod, only : mpp_root_pe, mpp_pe

      implicit none

      private
      public  :: moz_hook_init, moz_hook

!     save

!----------------------------------------------------------------------
!        Set global lightning NOx scaling factor
!----------------------------------------------------------------------
      real :: factor = 1.                    ! user-controlled scaling factor to achieve arbitrary NO prod.
      logical :: normalize_by_area = .false. ! normalize lightning NOx production by grid cell area
      logical :: allow_small_storms = .false. ! modify area normalization for high-res grids
      real :: min_land_frac = -999.          ! minimum land fraction for flash frequency calculation (default=-999)
      real, parameter :: AREA_PER_STORM = 1.e10 ! m2 (100km x 100km)
      real :: vdist(16,3)                    ! vertical distribution of lightning
      integer :: id_prod_no_col, id_flash_freq, id_prod_no_col_lght, id_eminox_lght
      real :: lat25
      real, parameter :: MW_N = 14.00674, &                 ! molecular weight of nirogen (AMU)
                         ONE_OVER_AVO = 1.65979e-24         ! reciprocal of Avogadros number (mole)
      
character(len=128), parameter :: version     = '$Id$'
character(len=128), parameter :: tagname     = '$Name$'
logical                       :: module_is_initialized = .false.

      CONTAINS

      subroutine moz_hook_init( lght_no_prd_factor, normalize_lght_no_prd_area, &
                                allow_small_storms_lght_no_prd, min_land_frac_lght, &
                                Time, axes, verbose )
!----------------------------------------------------------------------
!       ... Initialize the chemistry "hook" routine
!----------------------------------------------------------------------

      implicit none

!----------------------------------------------------------------------
!        ... Dummy args
!----------------------------------------------------------------------
      type(time_type), intent(in) :: Time
      integer,         intent(in) :: axes(4)
      real,            intent(in) :: lght_no_prd_factor         ! lightning no production factor
      logical,         intent(in) :: normalize_lght_no_prd_area ! normalize lightning NOx production by grid cell area
      logical,         intent(in) :: allow_small_storms_lght_no_prd ! modify area normalization for high-res grids
      real,            intent(in) :: min_land_frac_lght         ! minimum land fraction for flash frequency calculation (default=-999)
      integer,         intent(in) :: verbose

!----------------------------------------------------------------------
!        ... local variables
!----------------------------------------------------------------------

      if (module_is_initialized) return

      factor = lght_no_prd_factor
      normalize_by_area = normalize_lght_no_prd_area
      allow_small_storms = allow_small_storms_lght_no_prd
      min_land_frac = min_land_frac_lght
      if (verbose >= 2) then
         if (mpp_root_pe().eq.mpp_pe()) then
            write(*,*) 'MOZ_HOOK_INIT: Lightning NO production scaling factor = ',factor
            if (normalize_lght_no_prd_area) then
               write(*,*) 'MOZ_HOOK_INIT: Normalize lightning NO production by grid cell area'
               if (allow_small_storms) &
                  write(*,*) 'MOZ_HOOK_INIT: ... if grid cell area .GT. AREA_PER_STORM'
            else
               write(*,*) 'MOZ_HOOK_INIT: Normalize lightning NO production by grid cell (not area)'
            end if
         end if
      end if

      lat25 = 25. * PI/180.

!----------------------------------------------------------------------
!       ... vdist(kk,itype) = % of lightning NOx between (kk-1) and (kk)
!           km for profile itype
!----------------------------------------------------------------------
      vdist(:,1) = (/ 20.1, 2.3, 0.8, 1.5, 3.4, 5.3, 3.6, 3.8, &       ! midlat cont
                       5.4, 6.6, 8.3, 9.6,12.8,10.0, 6.2, 0.3 /)
      vdist(:,2) = (/  5.8, 2.9, 2.6, 2.4, 2.2, 2.1, 2.3, 6.1, &       ! trop marine
                      16.5,14.1,13.7,12.8,12.5, 2.8, 0.9, 0.3 /)
      vdist(:,3) = (/  8.2, 1.9, 2.1, 1.6, 1.1, 1.6, 3.0, 5.8, &       ! trop cont
                       7.6, 9.6,10.5,12.3,11.8,12.5, 8.1, 2.3 /)
      id_prod_no_col = register_diag_field('tracers','prod_no_col',axes(1:2),Time, &
                                           'prod_no_col','TgN/y')
      id_prod_no_col_lght = register_diag_field('tracers','prod_no_col_lght',axes(1:2),Time, &
                                           'prod_no_col','molec cm-2 s-1')
      id_flash_freq  = register_diag_field('tracers','flash_freq',axes(1:2),Time, &
                                           'flash_freq','cm-2 s-1')
      id_eminox_lght = register_cmip_diag_field_2d ( 'tracers', 'eminox_lght', Time, &
                              'Total Emission Rate of NOx from lightning', 'kg m-2 s-1', &
                standard_name='tendency_of_atmosphere_mass_content_of_nox_expressed_as_nitrogen_due_to_emission')
      module_is_initialized = .true.
      
      end subroutine MOZ_HOOK_INIT


      subroutine moz_hook( cldtop, cldbot, oro, zm, zint, t, &
                           prod_no, area, lat, &
                           Time,is,js )
!----------------------------------------------------------------------
!        ... General purpose chemistry "hook" routine.
!           Update deposition velocity and sulfur input fields,
!           and calculate lightning NOx source & Rn emissions.
!----------------------------------------------------------------------


      implicit none

!----------------------------------------------------------------------
!        ... Dummy args
!----------------------------------------------------------------------
      integer, intent(in) :: cldtop(:,:)             ! cloud top level index
      integer, intent(in) :: cldbot(:,:)             ! cloud bottom level index
      real, intent(in) :: oro(:,:)                   ! orography "flag"
      real, intent(in) :: zm(:,:,:)                  ! geopot height above surface at midpoints (m)
      real, intent(in) :: zint(:,:,:)                ! geopot height above surface at interfaces (m)
      real, intent(in) :: t(:,:,:)                   ! temperature
      
      real, intent(out) :: prod_no(:,:,:)            ! production of NOx (molec cm^-3 s^-1)
      real, intent(in)  :: area(:,:)                 ! area (m^2)
      real, intent(in) :: lat(:,:)                   ! latitude
      type(time_type), intent(in) :: Time            ! time
      integer, intent(in) :: is, js                  ! lon,lat indices

!----------------------------------------------------------------------
!        ... Local variables
!----------------------------------------------------------------------
      integer, parameter :: LAND = 1, OCEAN = 0
      real, parameter    :: dayspy = 365.
      real, parameter    :: secpyr = dayspy * 8.64e4

      integer :: i, j, &
                 cldtind, &         ! level index for cloud top
                 cldbind            ! level index for cloud base > 273K
      integer :: k, kk, zlow_ind, zhigh_ind, itype
      real    :: frac_sum
      real       :: zlow, zhigh, zlow_scal, zhigh_scal, fraction
      real, dimension( size(prod_no,1),size(prod_no,2) ) :: &
                 dchgzone, &        ! depth of discharge zone [km]
                 cldhgt, &          ! cloud top height [km]
                 cgic, &            ! Cloud-Ground/Intracloud discharge ratio
                 flash_energy, &    ! Energy of flashes per second
                 glob_prod_no_col   ! Global NO production rate for diagnostics
      real :: prod_no_col(size(prod_no,1),size(prod_no,2)) ! production of NOx (molec cm^-2 s^-1)
      real :: flash_freq(size(prod_no,1),size(prod_no,2))  
      real :: local_area(size(area,1),size(area,2)) ! storm area (cm^2)
      logical :: used
      integer :: platl, plonl, plev
!----------------------------------------------------------------------
!         ... Parameters to determine CG/IC ratio [Price and Rind, 1993]
!----------------------------------------------------------------------
      real, parameter  :: ca = .021, cb = -.648, cc = 7.49, cd = -36.54, ce = 64.09


      plonl = size(prod_no,1)
      platl = size(prod_no,2)
      plev  = size(prod_no,3)

!----------------------------------------------------------------------
!        Lightning NO production : Initialize ...
!----------------------------------------------------------------------
      flash_freq(:,:) = 0.
      cldhgt(:,:)     = 0.
      dchgzone(:,:)   = 0.
      cgic(:,:)       = 0.
      flash_energy(:,:) = 0.
      prod_no(:,:,:)    = 0.
      prod_no_col(:,:)  = 0.
      glob_prod_no_col(:,:) = 0.
      if (normalize_by_area) then 
         if (allow_small_storms) then
            local_area(:,:) = MIN(area,AREA_PER_STORM)*1.e4
         else
            local_area(:,:) = AREA_PER_STORM*1.e4
         end if
      else
         local_area(:,:) = area(:,:)*1.e4
      end if

!----------------------------------------------------------------------
!        Check whether tropchem is active
!----------------------------------------------------------------------
      if (.not. module_is_initialized) return

!--------------------------------------------------------------------------------
!        ... Estimate flash frequency and resulting NO emissions
!           [Price, Penner, Prather, 1997 (JGR)]
!    Lightning only occurs in convective clouds with a discharge zone, i.e.
!    an altitude range where liquid water, ice crystals, and graupel coexist.
!    We test this by examining the temperature at the cloud base.
!    It is assumed that only one thunderstorm occurs per grid box, and its
!    flash frequency is determined by the maximum cloud top height (not the
!    depth of the discharge zone). This is somewhat speculative but yields
!    reasonable results.
!
!       The CG/IC ratio is determined by an empirical formula from Price and
!    Rind [1993]. The average energy of a CG flash is estimated as 6.7e9 J,
!    and the average energy of a IC flash is assumed to be 1/10 of that value.
!       The NO production rate is assumed proportional to the discharge energy
!    with 1e17 N atoms per J. The total number of N atoms is then distributed
!    over the complete column of grid boxes.
!--------------------------------------------------------------------------------
         do j = 1,platl
            do i = 1,plonl
!--------------------------------------------------------------------------------
!         ... Find cloud top and bottom level above 273K
!--------------------------------------------------------------------------------
               cldtind =  cldtop(i,j) 
               cldbind =  cldbot(i,j) 
               do
                  if( cldbind <= cldtind ) then
                    exit
                  else if ( t(i,j,cldbind) < 273. ) then
                    exit
                  end if
                  cldbind = cldbind - 1
               end do
               if( cldtind < plev .and. cldtind > 0 .and. cldtind < cldbind ) then
!--------------------------------------------------------------------------------
!       ... Compute cloud top height and depth of charging zone
!--------------------------------------------------------------------------------
                  cldhgt(i,j) = 1.e-3*MAX( 0.,zint(i,j,cldtind) )
                  dchgzone(i,j) = cldhgt(i,j)-1.e-3*zm(i,j,cldbind)
!--------------------------------------------------------------------------------
!       ... Compute flash frequency for given cloud top height
!           (flashes storm^-1 min^-1)
!--------------------------------------------------------------------------------
                  if( ( NINT(oro(i,j))==LAND .and. min_land_frac<0. ) .or. &
                      ( oro(i,j) > min_land_frac .and. min_land_frac>=0. ) ) then
                     flash_freq(i,j) = 3.44e-5 * cldhgt(i,j)**4.9 
                  else
                     flash_freq(i,j) = 6.40e-4 * cldhgt(i,j)**1.7
                  end if
!--------------------------------------------------------------------------------
!       ... Compute CG/IC ratio
!           cgic = proportion of CG flashes (=PG from PPP paper)
!--------------------------------------------------------------------------------
                  cgic(i,j) = 1./((((ca*dchgzone(i,j) + cb)*dchgzone(i,j) + cc) &
                                      *dchgzone(i,j) + cd)*dchgzone(i,j) + ce)
                  if( dchgzone(i,j) < 5.5 ) then
                     cgic(i,j) = 0.
                  end if
                  if( dchgzone(i,j) > 14. ) then
                     cgic(i,j) = .02
                  end if
!--------------------------------------------------------------------------------
!       ... Compute flash energy (CG*6.7e9 + IC*6.7e8)
!           and convert to total energy per second
!--------------------------------------------------------------------------------
                  flash_energy(i,j) = cgic(i,j)*6.7e9 + (1. - cgic(i,j))*6.7e8
                  flash_energy(i,j) = flash_energy(i,j)*flash_freq(i,j)/60.

!--------------------------------------------------------------------------------
!         ... Compute number of N atoms produced per second
!           and convert to N atoms per second per cm2 and apply fudge factor
!         ... If (normalize_by_area), then assume storm is fixed area, and scale
!           flashes to grid cell. Otherwise, assume one storm per grid cell.
!         ... If (allow_small_storms), then assume storm occupies full grid cell, if
!           grid cell area .LT. AREA_PER_STORM
!--------------------------------------------------------------------------------
                  prod_no_col(i,j) = 1.e17*flash_energy(i,j) / local_area(i,j) * factor
!--------------------------------------------------------------------------------
!         ... Compute global NO production rate in TgN/yr:
!           TgN per second: * MW_N * 1.e-12 / AVO
!           TgN per year: * secpyr
!--------------------------------------------------------------------------------
                  glob_prod_no_col(i,j) = 1.e17*flash_energy(i,j) &
                                        * MW_N * ONE_OVER_AVO * 1.e-12 * secpyr * factor
               end if
            end do
         end do
         if(id_prod_no_col >0) &
            used=send_data(id_prod_no_col,glob_prod_no_col,Time,is_in=is,js_in=js)
         if(id_prod_no_col_lght >0) &
            used=send_data(id_prod_no_col_lght,prod_no_col,Time,is_in=is,js_in=js)
!--------------------------------------------------------------------------------
!         ... Convert from molec cm-2 s-1 to kgN m-2 s-1
!           * (MW_N*1e-3) * (1/AVO) * 1e4
!--------------------------------------------------------------------------------
         if(id_eminox_lght >0) &
            used=send_data(id_eminox_lght,prod_no_col*MW_N*ONE_OVER_AVO*10.,Time,is_in=is,js_in=js)
         if(id_flash_freq >0) &
            used=send_data(id_flash_freq,flash_freq/60./local_area,Time,is_in=is,js_in=js)


!--------------------------------------------------------------------------------
!        ... distribute production up to cloud top [Pickering et al., 1998 (JGR)]
!--------------------------------------------------------------------------------
            do j = 1,platl
               do i = 1,plonl
                  cldtind =  cldtop(i,j) 
                   if( prod_no_col(i,j) > 0. ) then
                     if( cldhgt(i,j) > 0. ) then
                        if(  ABS(lat(i,j)) > lat25 ) then
                           itype = 1                              ! midlatitude continental
                        else if ( ( NINT(oro(i,j))==LAND .and. min_land_frac<0. ) .or. &
                                  ( oro(i,j) > min_land_frac .and. min_land_frac>=0. ) ) then
                           itype = 3                              ! tropical continental
                        else
                           itype = 2                              ! topical marine
                        end if
                        frac_sum = 0.
                        
                        do k = cldtind,plev
                           zlow       = zint(i,j,k+1) * 1.e-3   ! lower interface height (km)
                           zlow_scal  = zlow * 16./cldhgt(i,j)  ! scale to 16 km convection height
                           zlow_ind   = MAX( 1,INT(zlow_scal)+1 )  ! lowest vdist index to include in layer
                           zhigh      = zint(i,j,k) * 1.e-3     ! upper interface height (km)
                           zhigh_scal = zhigh * 16./cldhgt(i,j) ! height (km) scaled to 16km convection height
                           zhigh_ind  = MAX( 1,MIN( 16,INT(zhigh_scal)+1 ) )  ! highest vdist index to include in layer
                           do kk = zlow_ind,zhigh_ind
                              fraction = MIN( zhigh_scal,REAL(kk) ) &         ! fraction of vdist in this model layer
                                         - MAX( zlow_scal,REAL(kk-1) )
                              fraction = MAX( 0., MIN( 1.,fraction ) )
                              frac_sum = frac_sum + fraction*vdist(kk,itype)
                              prod_no(i,j,k) = prod_no(i,j,k) &         ! sum the fraction of column NOx in layer k
                                             + fraction*vdist(kk,itype)*.01
                           end do
                           prod_no(i,j,k) = prod_no_col(i,j) * prod_no(i,j,k) & ! multiply fraction by column amount
                                               / (1.e5*(zhigh - zlow))                   ! and convert to atom N cm^-3 s^-1
                        end do
                     end if
                  end if
               end do
            end do

      end subroutine MOZ_HOOK

      end module MOZ_HOOK_MOD
