module gw_beres_mod
  use fms, only: FATAL, fms_mpp_npes, fms_mpp_pe, fms_mpp_root_pe, &
                 fms_mpp_stdlog, check_nml_error, fms_mpp_input_nml_file, &
                 fms2_io_open_file, FmsNetcdfFile_t, fms_mpp_get_current_pelist, &
                 fms2_io_get_dimension_size, fms2_io_read_data, &
                 fms2_io_close_file
  use mpp_mod, only: mpp_error
  use fms_mod, only: write_version_number
  use FMSconstants,     only: PI, RAD_TO_DEG, SECONDS_PER_DAY
  use gw_common_utils,        only: GWBand, energy_momentum_adjust, &
                gw_drag_prof, gw_prof, gw_common_init, gw_common_end
  use diag_manager_mod,       only: register_diag_field, send_data
  use time_manager_mod,       only:  time_type
  implicit none
  private

!-----------------------------------------------------
! Parameterization for convective gravity waves based on
! Beres et al. 2004, A method of specifying the gravity wave spectrum above convection
! based on latent heating properties and background wind, JAS, 61, 324-337. 
! Source codes are modified from the implementation in the GEOS model in Sep 2023.
! The GEOS implementation including additional heating rate-based GW source over the
! extratropics, representing GW from frontogenesis.
! Dec 2023: include waves corresponds to the deep plume only.
!-----------------------------------------------------

!---------------------------------------------------------------------
!----------- ****** VERSION NUMBER ******* ---------------------------
character(len=128)  :: version =  '$Id$'
character(len=128)  :: tagname =  '$Name$'

!---------------------------------------------------------------------
!-------  interfaces --------
  public :: gw_beres_init, gw_beres_ifc, gw_beres_end

  private :: gw_beres_src, gw_beres_readfile

!------Declar BeresSourceDesc type--------
  type :: BeresSourceDesc
    logical           :: active
    logical           :: storm_shift     !< Whether wind speeds are shifted to be relative to storm cells.
    real              :: min_hdepth      !< Heating depths below this value [m] will be ignored.
    real              :: spectrum_source !< Source for wave spectrum
    real, allocatable :: k(:)            !< Index for level where wind speed is used as the source speed.
    real              :: tndmax          !< tendency limiter
    integer           :: maxh            !< Table bounds, for convenience. (Could be inferred from shape(mfcc).)
    integer           :: maxuh           !< Table bounds, for convenience. (Could be inferred from shape(mfcc).)
    real, allocatable :: hd(:)           !< Heating depths [m].
    real, allocatable :: mfcc(:,:,:)     !< Table of source spectra.
    real, allocatable :: taubck(:,:)     !< Forced background for extratropics
  end type BeresSourceDesc

  ! module variables here
  logical :: module_is_initialized = .false.
  type(GWBand) :: beres_band
  type(BeresSourceDesc) :: beres_dc_desc ! Beres DeepCu
!  type(BeresSourceDesc) :: beres_sc_desc ! Beres ShallowCu
  real                                 :: min_hdepth=1000.0 ! heating depths below this value will be ignored in unit of meters.
  real                                 :: fcrit2=1.0 ! Critical Froude number square
  real                                 :: bkg_wavelength=1.0E5 ! Horizontal wavelength in meters.
  logical                             :: do_beres_dc=.true. ! Activate gravity waves generated from the parameterized deep convection.

  !-----------------------------------------
  !----- namelist ------
  character(len=255) :: BERES_FILE_NAME = 'newmfspectra40_dc25.nc' ! Beres Scheme File containing the lookup table for deep convection
  real               :: BERES_EFFGW=1.0 ! Efficiency of gravity wave momentum transfer, tuning parameter
  real               :: BERES_DC_SRC_LEVEL = 70000.0 ! pressure level of the gravity wave source associated with deep convection in Pa
!  real               :: BERES_EXTROP_SCALING =  50.0 ! Extratropical background frontalgenesis GW stress source strength
  real               :: BERES_NH = 50.0 ! NH extratropical background GW stress source amplitude
  real               :: BERES_SH = 50.0 ! SH extratropical background GW stress source amplitude
  real               :: BERES_BKG_TNDMAX = 100.0 ! Maximum wind tendency from stress divergence in m/s/day
  real               :: HRCF = 20 ! If >0, the inverse of the fractional area covered by deep convection in the grid cell. If <0, use progonstic convection area, and the (-HRCF) is timed to HRCF.
  real               :: QBO_HDEPTH_SCALING=0.25 ! tuning factor to reduce GW phase speeds for better QBO
  logical           :: DO_STORM_SHIFT = .true. ! Whether wind speeds are shifted to be relative to storm cells.
  integer          :: BERES_PGWV = 32 ! Maximum wave number and width of spectrum bins
!  logical           :: DO_BERES_SC = .true. ! Whether to active shallow convective GW

  namelist / gw_beres_nml / BERES_FILE_NAME, BERES_EFFGW, BERES_DC_SRC_LEVEL, &
          BERES_NH, BERES_SH, BERES_BKG_TNDMAX, HRCF, QBO_HDEPTH_SCALING, &
          DO_STORM_SHIFT, BERES_PGWV

!----------------------------------------
!   variables for netcdf diagnostic fields.
integer            :: id_udt_beres, id_vdt_beres, id_tdt_beres, id_hdepth_beres, id_maxq_beres
real                 :: missing_value = -999.
character(len=8) :: mod_name = 'gw_beres'

  contains

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!                      PUBLIC SUBROUTINES
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!####################################################################

  subroutine gw_beres_init(lonb, latb, pref_in, Time, axes)
!--------------------------------------------------------------
! Read namelist, set up phase speed coordinate (band), set up beres_dc_desc
! Read in the look-up table for source spectra
! This subroutine cooresponds to gw_common_init and gw_beres_init in GEOS
!-------------------------------------------------------------------------------
    real, dimension(:,:),       intent(in) :: lonb, latb
    real, dimension(:),          intent(in) :: pref_in
    type(time_type),         intent(in) :: Time
    integer, dimension(4), intent(in) :: axes
!-------------------------------------------------------------------
!   intent(in) variables:
!
!       lonb      2d array of model longitudes on cell corners [radians]
!       latb      2d array of model latitudes at cell corners [radians]
!       pref      array of reference pressures at full levels (plus
!                 surface value at nlev+1), based on 1013.25hPa pstar
!                 [ Pa ]
!       Time      current time (time_type)
!       axes      data axes for diagnostics
!
!------------------------------------------------------------------

!   local variables
    integer                         :: io, ierr !< Error codes when reading the namelist
    integer                         :: logunit
    type(FmsNetcdfFile_t)           :: fileobj
    integer,            allocatable :: pes(:) !< Array of the pes in the current pelist

!    real                            :: min_hdepth
!    logical                         :: storm_shift
!    real                             :: dc !< phase speed coordinate interval from the input file
    integer                      :: nlon, nlat, ncol, nlev, j, i
    real                           :: lats((size(lonb,1) -1)*(size(latb,2)-1))
    real                           :: pref_half(size(pref_in,1))

!    if routine has already been executed, return.
    if (module_is_initialized) return

!    read namelist.
    read (fms_mpp_input_nml_file, nml=gw_beres_nml, iostat=io)
    ierr = check_nml_error(io,"gw_beres_nml")

!    write version number and namelist to logfile.
    call write_version_number( version, tagname)
    logunit = fms_mpp_stdlog()
    if (fms_mpp_pe() == fms_mpp_root_pe()) write (logunit, nml=gw_beres_nml)

!  define the grid dimensions. idf and jdf are the (i,j) dimensions of
!  domain on this processor. Calculate latitude at the center of each grid.
!  Conver latitude array to latitude column.
      nlat  = size(latb,2) - 1
      nlon  = size(lonb,1) - 1
      ncol = nlat * nlon
      do j=1,nlat
        do i=1,nlon
          lats(i+(j-1)*nlon)=  0.5*( latb(i,j+1)+latb(i,j) )
        end do
     end do

  ! reference phalf level
   nlev = size(pref_in) -1
   pref_half(nlev+1)=pref_in(nlev+1)
   pref_half(2:nlev)=sqrt(pref_in(1:nlev-1)*pref_in(2:nlev))
   pref_half(1)=pref_in(1)**2/pref_half(2)

    !< Read the beres file, and pass relevent info to beres_dc_desc
    if (do_beres_dc) then
        call gw_beres_readfile(beres_dc_desc, beres_band, ncol, nlev, lats, pref_half, BERES_FILE_NAME, BERES_DC_SRC_LEVEL, min_hdepth, BERES_PGWV)
    end if
    ! No shallow convection is currently used to generate GW. If shallow convection GW is used, a separate input file will be
    ! be provided. Following commands are to read the ShallowCu look-up table and set up beres_sc_desc.
    ! if (do_beres_sc) then
    !  call gw_beres_readfile(beres_sc_desc, beres_band, beres_sc_file, beres_sc_src_level, 0.0, ncl, nlats)
    ! end if

    call gw_common_init (pref_half)
    !< Set up diagnostics
    id_udt_beres = register_diag_field(mod_name, 'udt_beres', axes(1:3), Time,&
            'zonal wind tendency from Beres scheme', 'm/s^2', missing_value=missing_value)
    id_vdt_beres = register_diag_field(mod_name, 'vdt_beres', axes(1:3), Time,&
            'meridional wind tendency from Beres scheme', 'm/s^2', missing_value=missing_value)
    id_tdt_beres = register_diag_field(mod_name, 'tdt_beres', axes(1:3), Time,&
            'temperature tendency from Beres scheme', 'K/s', missing_value=missing_value)
    id_hdepth_beres = register_diag_field (mod_name, 'hdepth_beres', axes(1:2), Time, &
            'depth of the convective heating used in Beres scheme', 'm', missing_value=missing_value)
    id_maxq_beres = register_diag_field (mod_name, 'maxq_beres', axes(1:2), Time, &
            'maxium convective heating rate used in Beres scheme', 'K/s', missing_value=missing_value)

    module_is_initialized = .true.

  end subroutine gw_beres_init

  subroutine gw_beres_ifc(is, js, lat_in, Time, dt, pfull, phalf, zfull, &
  u_in, v_in, temp_in, dtconv_in, cqa_in, utgw_out, vtgw_out, ttgw_out, &
  dqcdt_in )
    integer,                          intent(in) :: is, js
    real,  dimension(:,:),       intent(in) :: lat_in      ! latitudes
    type(time_type),            intent(in) :: Time       ! Time
    real,                                intent(in) :: dt           ! Time step.
    real, dimension(:,:,:),      intent(in) :: pfull, phalf, zfull, u_in, v_in, temp_in
    real, dimension(:,:,:),      intent(in) :: dtconv_in, cqa_in ! grid-mean convective heating rate and updraft area from parameterized convection

    real, dimension(:,:,:),      intent(out):: utgw_out, vtgw_out, ttgw_out
    real, optional, intent(in), dimension(:,:,:) :: dqcdt_in
!---------------------------Local storage-------------------------------
!-------------------- interface variables reshaped into 2D-------------
    integer      :: ncol, nx, ny       ! number of atmospheric columns
    integer      :: pver         ! number of vertical layers

    real          :: u(size(u_in,1)*size(u_in,2),size(u_in,3))      ! Midpoint zonal winds. ( m s-1)
    real          :: v(size(u_in,1)*size(u_in,2),size(u_in,3))      ! Midpoint meridional winds. ( m s-1)
    real          :: t(size(u_in,1)*size(u_in,2),size(u_in,3))      ! Midpoint temperatures. (K)
    real          :: netdt(size(u_in,1)*size(u_in,2),size(u_in,3))  ! Local (instead of grid-mean) convective heating rate (K s-1)
    real          :: pint(size(u_in,1)*size(u_in,2),size(u_in,3)+1) ! Interface pressures (phalf). (Pa)
    real          :: pmid(size(u_in,1)*size(u_in,2),size(u_in,3))  ! Midpoint pressure (pfull) (Pa)
    real          :: zm(size(u_in,1)*size(u_in,2),size(u_in,3))     ! Midpoint altitudes above ground (m).

    real          :: lats(size(u_in,1)*size(u_in,2))      ! latitudes (radius)

    real          :: utgw(size(u_in,1)*size(u_in,2),size(u_in,3))       ! zonal wind tendency (m/s/s)
    real          :: vtgw(size(u_in,1)*size(u_in,2),size(u_in,3))       ! meridional wind tendency (m/s/s)
    real          :: ttgw(size(u_in,1)*size(u_in,2),size(u_in,3))       ! temperature tendency (K/s)
 
    real, allocatable, dimension(:,:) :: dqcdt  ! Condensate tendency due to large-scale (kg kg-1 s-1)
 
    !---------------------------Other local storage-------------------------------
 
    integer :: k, m, nn

   ! reference phalf level (Pa)
    real :: prefint(size(u_in,3)+1)
    ! air density at phalf (kg m-3)
    real :: rhoi(size(u_in,1)*size(u_in,2), size(u_in,3)+1)
    ! Brunt-Vaisalla frequencies at phalf (s-1)
    real :: ni(size(u_in,1)*size(u_in,2), size(u_in,3)+1)

    real, allocatable :: tau(:,:,:)  ! wave Reynolds stress
    ! gravity wave wind tendency for each wave
    real, allocatable :: gwut(:,:,:)
   ! ! Wave phase speeds for each column
   ! real, allocatable :: c(:,:)
 
    ! Efficiency for a gravity wave source.
    real :: effgw(size(u_in,1)*size(u_in,2))
 
    ! Momentum fluxes used by fixer.
    real :: um_flux(size(u_in,1)*size(u_in,2)), vm_flux(size(u_in,1)*size(u_in,2))
 
    ! Energy change used by fixer.
    real :: de(size(u_in,1)*size(u_in,2))
 
    ! Reynolds stress for waves propagating in each cardinal direction.
    real :: taucd(size(u_in,1)*size(u_in,2),size(u_in,3)+1,4)
 
    ! Indices of top gravity wave source level and lowest level where wind
    ! tendencies are allowed.
    integer :: src_level(size(u_in,1)*size(u_in,2))
    integer :: tend_level(size(u_in,1)*size(u_in,2))
 
    ! Projection of wind at midpoints and interfaces.
    real :: ubm(size(u_in,1)*size(u_in,2),size(u_in,3))
    real :: ubi(size(u_in,1)*size(u_in,2),size(u_in,3)+1)
 
    ! Unit vectors of source wind (zonal and meridional components).
    real :: xv(size(u_in,1)*size(u_in,2))
    real :: yv(size(u_in,1)*size(u_in,2))
 
    ! Heating depth [m] and maximum heating in each column.
    real :: hdepth(size(u_in,1)*size(u_in,2)), maxq0(size(u_in,1)*size(u_in,2))
 
    real :: pint_adj(size(u_in,1)*size(u_in,2),size(u_in,3)+1)
    real :: zfac_layer
 
    character(len=1) :: cn
    character(len=9) :: fname(4)
 
    integer :: i,j,l
    ! for diagnostics output
    logical  :: used

    !----------------------------------------------------------------------------
   ! !debug
   ! if (fms_mpp_pe()==fms_mpp_root_pe()) then
   !    print*, 'maximum heating rates in Beres: ', maxval(dtconv_in)
   !    print*, 'before reshape u:', maxval(u_in)
   !    print*, 'before reshape t:', maxval(temp_in)
   ! endif

   ! convert 3D inputs into 2D
   nx = size(u_in, 1)
   ny = size(u_in, 2)
   ncol=nx*ny
   pver=size(u_in,3)
   u=reshape(u_in,(/ncol,pver/))
   v=reshape(v_in,(/ncol,pver/))
   t=reshape(temp_in,(/ncol,pver/))
   if (HRCF<0) then
! Use progonostic updraft area.
! (-HRCF) now corresponds to a constant factor adjusting the amplitude of heating rate
       do i=1,nx
         do j = 1,ny
           do k = 1, pver
             if ( cqa_in(i,j,k)>0 ) then
               netdt(i*j,k)=dtconv_in(i,j,k)/cqa_in(i,j,k)*(-HRCF)
             else
               netdt(i*j,k)=0.
             endif
           enddo
        enddo
     enddo
   else
! Assume a constant convetive area of 1/HRCF
      netdt=reshape(dtconv_in*HRCF, (/ncol, pver/))
   endif
   pint = reshape(phalf, (/ncol, pver/))
   pmid = reshape(pfull, (/ncol, pver/))
   zm = reshape(zfull, (/ncol, pver/))
   lats = reshape(lat_in, (/ncol/))
   if (present(dqcdt_in)) then
       allocate (dqcdt(ncol,pver))
       dqcdt = reshape(dqcdt_in, (/ncol, pver/))
   endif

 !  if (fms_mpp_pe()==fms_mpp_root_pe()) then
 !     print*, 'after reshape u:', maxval(u)
 !     print*, 'after reshape t:', maxval(t)
 !     print*, 'after reshape tdt', maxval(netdt)
 !  endif

    call gw_prof (ncol , pver, pint , pmid , t , rhoi, ni )

    ! Allocate wavenumber fields.
    allocate(tau(ncol,-beres_band%ngwv:beres_band%ngwv,pver+1))
    allocate(gwut(ncol,pver,-beres_band%ngwv:beres_band%ngwv))
   ! allocate(c(ncol,-beres_band%ngwv:beres_band%ngwv))
 
    ! Efficiency of gravity wave momentum transfer.
    ! This is really only to remove the pole points.
    where (pi/2.0 - abs(lats(:ncol)) >= 1.e-4 )  !-4*epsilon(1.0))
      effgw = BERES_EFFGW
    elsewhere
      effgw = 0.0
    end where
 
!   !debug
!    if (fms_mpp_pe()==fms_mpp_root_pe()) then
!       print*, 'before Beres source: ', maxval(dqcdt)
!    endif
    ! Determine wave sources for Beres deep scheme
    if (present(dqcdt_in)) then
      call gw_beres_src(ncol, pver, beres_band, beres_dc_desc, &
        u, v, netdt, zm, src_level, tend_level, tau, &
        ubm, ubi, xv, yv, hdepth, maxq0, dqcdt=dqcdt )
    else
      call gw_beres_src(ncol, pver, beres_band, beres_dc_desc, &
        u, v, netdt, zm, src_level, tend_level, tau, &
        ubm, ubi, xv, yv, hdepth, maxq0 )
    endif
    !WMP pressure scaling near model top
    !!!  pint_adj = 1.0
    zfac_layer = 300.0 ! 3mb
    do k=1,pver+1
      do i=1,ncol
        pint_adj(i,k) = MIN(1.0,MAX(0.0,(pint(i,k)/zfac_layer)**3))
      enddo
    enddo
 
    ! Solve for the drag profile with orographic sources.
    call gw_drag_prof(ncol, pver, beres_band, pint, &
      src_level, tend_level, dt, t,    &
      rhoi, ni, ubm, ubi, xv, yv, tau, utgw, vtgw, &
      ttgw, gwut, tau_adjust=pint_adj)

    ! Apply efficiency and limiters
    call energy_momentum_adjust(ncol, pver, beres_band, pint,  u, v, dt, tau, &
                                effgw, t, ubm, ubi, xv, yv, utgw, vtgw, ttgw, &
                                tend_level, tndmax_in=beres_dc_desc%tndmax)
  
    deallocate(tau, gwut)

  ! convert 2D outputs to 3D
  utgw_out = reshape(utgw, (/nx, ny, pver/))
  vtgw_out = reshape(vtgw, (/nx, ny, pver/))
  ttgw_out = reshape(ttgw, (/nx, ny, pver/))

 ! send out diagnostics
    if ( id_udt_beres > 0 ) then
        used = send_data (id_udt_beres, utgw_out, Time, is, js)
    endif
    if ( id_vdt_beres > 0 ) then
        used = send_data (id_vdt_beres, vtgw_out, Time, is, js)
    endif
    if ( id_tdt_beres > 0 ) then
        used = send_data (id_tdt_beres, ttgw_out, Time, is, js)
    endif

    if (id_hdepth_beres > 0) then
       used = send_data (id_hdepth_beres, reshape(hdepth, (/nx, ny/)), &
                                       Time, is, js)
    endif
   if (id_maxq_beres > 0) then
       used = send_data (id_maxq_beres, reshape(maxq0, (/nx, ny/)), &
                                       Time, is, js)
    endif
  end subroutine gw_beres_ifc

!---------------------------------------------------------------------
  subroutine gw_beres_end
!---------------------------------------------------------------------
!    deallocate the module variables.
!---------------------------------------------------------------------
     beres_band%is_set = .false.
     deallocate(beres_band%cref)
     deallocate(beres_dc_desc%k, beres_dc_desc%hd, beres_dc_desc%mfcc, &
               beres_dc_desc%taubck)
     call gw_common_end
!---------------------------------------------------------------------
!    mark the module as uninitialized.
!---------------------------------------------------------------------
      module_is_initialized = .false.
  end subroutine gw_beres_end

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!                      PRIVATE SUBROUTINES
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!####################################################################
  !------------------------------------
  subroutine gw_beres_readfile(desc, band, ncol, pver, lats, prefint, file_name, spectrum_source_in, min_hdepth_in, pgwv)
  !----------------------------------------
  ! Read input file, set up the phase speed band structure and pass the relevant information to "desc".
  !------------------------------Arguments--------------------------------
      type(BeresSourceDesc), intent(inout) :: desc
      type(GWBand),           intent(inout)    :: band
      character(len=255),    intent(in)  :: file_name
      real,                   intent(in)    :: spectrum_source_in
      real,                   intent(in)    :: min_hdepth_in
      integer,              intent(in)    :: ncol
      integer,              intent(in)    :: pver
      real,                   intent(in)    :: lats(ncol)
      real,                   intent(in)    :: prefint(pver+1) !< reference phalf
      integer,              intent(in)    :: pgwv
!-----------------local variables----------------------------
    type(FmsNetcdfFile_t)           :: fileobj
    integer,            allocatable :: pes(:) !< Array of the pes in the current pelist
    real, allocatable :: mfcc(:,:,:) !< loop-up table for GW source as a function of heating depth, mean wind and phase speed
    real, allocatable :: hdcc(:) !<heating depth coordinate of the look-up table
    real, allocatable :: pscc(:) !< phase speed coordinate of the look-up table
    integer                         :: nhd_file !<dimension of the heating depth coordinate
    integer                         :: nmw_file !<dimension of the mean wind coordinate
    integer                         :: nps_file !< dimension of the phase speed coordinate
    integer                         :: nps_file_half !< number of positive phase speed
    integer                         :: nps_start !<index of the first phase speed coordinate to use
    integer                         :: nmw_file_half !< number of positive mean wind
    real                              :: dc !< phase speed coordinate interval
      ! For forced background extratropical wave speed
      real    :: c4, latdeg, flat_gw
      real, allocatable :: c0(:), cw4(:)
      integer :: i, kc

      ! read input file
    allocate(pes(fms_mpp_npes()))
    call fms_mpp_get_current_pelist(pes)
    if (.not. fms2_io_open_file(fileobj, trim(file_name), "read", pelist=pes)) &
      call mpp_error(FATAL, "gw_beres_readfile:: Unable to open the file "//trim(file_name))

    call fms2_io_get_dimension_size(fileobj, "PS", nps_file)
    call fms2_io_get_dimension_size(fileobj, "MW", nmw_file)
    call fms2_io_get_dimension_size(fileobj, "HD", nhd_file)

    allocate( mfcc(nhd_file, nmw_file, nps_file) )
    allocate( hdcc(nhd_file) )
    allocate(pscc(nps_file))

    call fms2_io_read_data(fileobj, "HD", hdcc)
    call fms2_io_read_data(fileobj, "mfcc", mfcc)
    call fms2_io_read_data(fileobj, "PS", pscc)
    call fms2_io_close_file(fileobj)

    dc=pscc(2)-pscc(1)
    nps_file_half= (nps_file-1)/2
    nps_start = nps_file_half - pgwv + 1
    nmw_file_half = (nmw_file - 1)/2

    if (beres_band%is_set) then ! If the band is already set up, check if it is consistent with the input file
       if (dc .ne. beres_band%dc ) &
          call mpp_error(FATAL, "gw_beres_readfile:: dc is inconsistent with input file "//trim(file_name))
       if (nps_file_half .lt. beres_band%ngwv ) &
          call mpp_error(FATAL, "gw_beres_readfile:: Not enough ngwv in input file "//trim(file_name))
    else ! band has not been set up. Create a new "band" following the input file.
       call beres_band%create_new_GWband(pgwv, dc, fcrit2, bkg_wavelength )
    endif

      ! Get HD (heating depth) dimension.
      desc%maxh = nhd_file

      ! Get MW (mean wind) dimension.
      desc%maxuh = nmw_file_half

      allocate(desc%hd(nhd_file))
      allocate(desc%mfcc(nhd_file,-nmw_file_half:nmw_file_half,-pgwv:pgwv))

      desc%mfcc=mfcc(:,:,nps_start:(nps_start+2*pgwv))

      ! While not currently documented in the file, it uses kilometers. Convert
      ! to meters.
      desc%hd = hdcc * 1000.0

      ! Source level index allocated, filled later
      desc%spectrum_source = spectrum_source_in
      allocate(desc%k(ncol))

      do i = 0, pver
      ! spectrum source index
         if (prefint(i+1) < beres_dc_desc%spectrum_source) beres_dc_desc%k(:) = i+1
      end do

      desc%min_hdepth = min_hdepth_in

      desc%storm_shift = DO_STORM_SHIFT

      desc%tndmax = BERES_BKG_TNDMAX

      ! Intialize forced background wave speeds
      ! Represent extratropical frontalgenesis GW
      allocate(desc%taubck(ncol,-band%ngwv:band%ngwv))
      allocate(c0(-band%ngwv:band%ngwv))
      allocate(cw4(-band%ngwv:band%ngwv))
      desc%taubck = 0.0
      c0  = 0.0
      cw4 = 0.0
      do kc = -4,4
         c4 =  10.0*kc
         cw4(kc) =  exp(-(c4/30.)**2)
      enddo
      do kc = -band%ngwv,band%ngwv
         c0(kc) =  10.0*(4.0/real(band%ngwv))*kc
         desc%taubck(:,kc) =  exp(-(c0(kc)/30.)**2)
      enddo
      do i=1,ncol
         ! include forced background stress in extra tropics
         ! Determine the background stress at c=0
         ! Include dependence on latitude:
         latdeg = lats(i)*RAD_TO_DEG
!         if (-15.3 < latdeg .and. latdeg < 15.3) then
!            flat_gw =  0.10
!         else if (latdeg > -31. .and. latdeg <= -15.3) then
!            flat_gw =  0.10
!         else if (latdeg <  31. .and. latdeg >=  15.3) then
!            flat_gw =  0.10
!         else if (latdeg > -60. .and. latdeg <= -31.) then
!            flat_gw =  0.50*exp(-((abs(latdeg)-60.)/23.)**2)
!         else if (latdeg <  60. .and. latdeg >=  31.) then
!            flat_gw =  0.50*exp(-((abs(latdeg)-60.)/23.)**2)
!         else if (latdeg <= -60.) then
!            flat_gw =  0.50*exp(-((abs(latdeg)-60.)/70.)**2)
!         else if (latdeg >=  60.) then
!            flat_gw =  0.50*exp(-((abs(latdeg)-60.)/70.)**2)
!         end if

! Modified profile, NH and SH amplitude are controlled separately, no flat
          if ( latdeg >= 0. .and. latdeg < 60.) then
             flat_gw = BERES_NH*0.50*exp(-((latdeg-60.)/23.)**2)
          else if ( latdeg >= 60.) then
             flat_gw = BERES_NH* 0.50*exp(-((latdeg-60.)/70.)**2)
          else if ( latdeg < 0. .and. latdeg > -60. ) then
             flat_gw = BERES_SH*0.50*exp(-((latdeg+60.)/23.)**2)
          else if ( latdeg <= -60. ) then
             flat_gw = BERES_SH*0.50*exp(-((latdeg+60.)/70.)**2)
          end if

! use the latitudinal profile of source_amp from AD99
!         flat_gw=BERES_NH*0.2*(1+tanh((latdeg-30.)/5.))+BERES_SH*0.2*(1+tanh(-(latdeg+30.)/5.))
         desc%taubck(i,:) = 0.001*flat_gw*desc%taubck(i,:)*(sum(cw4)/sum(desc%taubck(i,:)))
      enddo
      deallocate( c0, cw4 )
      deallocate( mfcc, hdcc, pscc )
      deallocate(pes)
    end subroutine gw_beres_readfile


  !------------------------------------
  subroutine gw_beres_src(ncol, pver, band, desc, u, v, &
     netdt, zm, src_level, tend_level, tau, ubm, ubi, xv, yv, &
     hdepth, maxq0, dqcdt)
!-----------------------------------------------------------------------
! Driver for multiple gravity wave drag parameterization.
!
! The parameterization is assumed to operate only where water vapor
! concentrations are negligible in determining the density.
!
! Beres, J.H., M.J. Alexander, and J.R. Holton, 2004: "A method of
! specifying the gravity wave spectrum above convection based on latent
! heating properties and background wind". J. Atmos. Sci., Vol 61, No. 3,
! pp. 324-337.
!
!-----------------------------------------------------------------------
  !use these subroutines/functions: get_unit_vector, dot_2d, GWBand

!------------------------------Arguments--------------------------------
  ! Column and vertical dimensions.
  integer, intent(in) :: ncol, pver

  ! Wavelengths triggered by convection.
  type(GWBand), intent(in) :: band

  ! Settings for convection type (e.g. deep vs shallow).
  type(BeresSourceDesc), intent(in) :: desc

  ! Midpoint zonal/meridional winds.
  real, intent(in) :: u(ncol,pver), v(ncol,pver)
  ! Local heating rate due to convection (K/s)
  ! Convection area is already counted
  real, intent(in) :: netdt(:,:)
  ! Midpoint altitudes.
  real, intent(in) :: zm(ncol,pver)

  ! Indices of top gravity wave source level and lowest level where wind
  ! tendencies are allowed.
  integer, intent(out) :: src_level(ncol)
  integer, intent(out) :: tend_level(ncol)

  ! Wave Reynolds stress.
  real, intent(out) :: tau(ncol,-band%ngwv:band%ngwv,pver+1)
  ! Projection of wind at midpoints and interfaces.
  real, intent(out) :: ubm(ncol,pver)
  real, intent(out) :: ubi(ncol,pver+1)
  ! Unit vectors of source wind (zonal and meridional components).
  real, intent(out) :: xv(ncol), yv(ncol)
  !! Phase speeds.
  !real(kind=r8_kind), intent(out) :: c(ncol,-band%ngwv:band%ngwv)

  ! Heating depth [m] and maximum heating in each column.
  real, intent(out) :: hdepth(ncol), maxq0(ncol)

  ! Condensate tendency due to large-scale (kg kg-1 s-1)
  real, optional, intent(in) :: dqcdt(ncol,pver)  ! Condensate tendency due to large-scale (kg kg-1 s-1)

!---------------------------Local Storage-------------------------------
  ! Column and level indices.
  integer :: i, k

  ! Zonal/meridional wind at roughly the level where the convection occurs.
  real :: uconv(ncol), vconv(ncol), ubi1d(ncol)

  ! Maximum heating rate.
  real :: q0(ncol)

  ! Bottom/top heating range index.
  integer  :: boti(ncol), topi(ncol)
  ! Index for looking up heating depth dimension in the table.
  integer  :: hd_idx(ncol)
  ! Mean wind in heating region.
  real :: uh(ncol)
  ! Min/max wavenumber for critical level filtering.
  integer :: Umini(ncol), Umaxi(ncol)
  ! Source level tau for a column.
  real :: tau0(-band%ngwv:band%ngwv)
  ! Speed of convective cells relative to storm.
  real :: CS(ncol)
  ! Index to shift spectra relative to ground.
  integer :: shift

  ! Averaging length.
  real, parameter :: AL = 1.0e5
  integer :: thread

  !----------------------------------------------------------------------
  ! Initialize tau array
  !----------------------------------------------------------------------
  tau = 0.0
  hdepth = 0.0
  q0 = 0.0
  tau0 = 0.0
  ubi = 0.0

  !------------------------------------------------------------------------
  ! Determine wind and unit vectors approximately at the source level, then
  ! project winds.
  !------------------------------------------------------------------------

  ! Source wind speed and direction.
  do i=1,ncol
   uconv(i) = u(i,desc%k(i))
   vconv(i) = v(i,desc%k(i))
  enddo

  ! Get the unit vector components and magnitude at the source level.
  ubi1d = 0.0
  call get_unit_vector(uconv, vconv, xv, yv, ubi1d)
  do i=1,ncol
   ubi(i,desc%k(i)+1) = ubi1d(i)
  enddo

  ! Project the local wind at midpoints onto the source wind.
  do k = 1, pver
     ubm(:,k) = dot_2d(u(:,k), v(:,k), xv, yv)
  end do

  ! Compute the interface wind projection by averaging the midpoint winds.
  ! Use the top level wind at the top interface.
  ubi(:,1) = ubm(:,1)

  ubi(:,2:pver) = (ubm(:,1:pver-1)+ubm(:,2:pver))/2.0

  !-----------------------------------------------------------------------
  ! Calculate heating depth.
  !
  ! Heating depth is defined as the first height range from the bottom in
  ! which heating rate is continuously positive.
  !-----------------------------------------------------------------------

  ! First find the indices for the top and bottom of the heating range.
  boti = 0
  topi = 0
  do k = pver, 1, -1
     do i = 1, ncol
        if (boti(i) == 0) then
           ! Detect if we are outside the maximum range (where z = 20 km).
           if (zm(i,k) >= 20000.0) then
              boti(i) = k
              topi(i) = k
           else
              ! First spot where heating rate is positive.
              if (netdt(i,k) > 0.0) boti(i) = k
           end if
        else if (topi(i) == 0) then
           ! Detect if we are outside the maximum range (z = 20 km).
           if (zm(i,k) >= 20000.0) then
              topi(i) = k
           else
              ! First spot where heating rate is no longer positive.
              if (netdt(i,k) <= 0.0) topi(i) = k
           end if
        end if
     end do
     ! When all done, exit.
     if (all(topi /= 0)) exit
  end do

  ! Heating depth in m.
  hdepth = [ ( (zm(i,topi(i))-zm(i,boti(i))), i = 1, ncol ) ]

  ! !debug
  !  if (fms_mpp_pe()==fms_mpp_root_pe()) then
  !     print*, 'maximum heating rates in Beres source: ', maxval(netdt)
  !     print*, 'max hdepth before: ', maxval(hdepth)
  !  endif

  ! J. Richter: this is an effective reduction of the GW phase speeds (needed to drive the QBO)
  hdepth = hdepth*QBO_HDEPTH_SCALING

  hd_idx = index_of_nearest(hdepth, desc%hd)

  ! hd_idx=0 signals that a heating depth is too shallow, i.e. that it is
  ! either not big enough for the lowest table entry, or it is below the
  ! minimum allowed for this convection type.
  ! Values above the max in the table still get the highest value, though.
  where (hdepth < max(desc%min_hdepth, desc%hd(1))) hd_idx = 0.0

  ! Maximum heating rate.
  do k = minval(topi), maxval(boti)
     where (k >= topi .and. k <= boti)
        q0 = max(q0, netdt(:,k))
     end where
  end do

  !output max heating rate in K/day
  maxq0 = q0*SECONDS_PER_DAY

  if (desc%storm_shift) then

     ! Find the cell speed where the storm speed is > 10 m/s.
     ! Storm speed is taken to be the source wind speed.
     do i=1,ncol
       CS(i) = sign(max(abs(ubm(i,desc%k(i)))-10.0, 0.0), ubm(i,desc%k(i)))
     enddo

     ! Average wind in heating region, relative to storm cells.
     uh = 0.0
     do k = minval(topi), maxval(boti)
        where (k >= topi .and. k <= boti)
           uh = uh + ubm(:,k)/(boti-topi+1)
        end where
     end do

     uh = uh - CS

  else

     ! For shallow convection, wind is relative to ground, and "heating
     ! region" wind is just the source level wind.
     do i=1,ncol
       uh(i) = ubm(i,desc%k(i))
     enddo

  end if

  ! Limit uh to table range.
  uh = min(uh,  real(desc%maxuh))
  uh = max(uh, -real(desc%maxuh))

  ! Speeds for critical level filtering.
  Umini =  band%ngwv
  Umaxi = -band%ngwv
  do k = minval(topi), maxval(boti)
     where (k >= topi .and. k <= boti)
        Umini = min(Umini, nint(ubm(:,k)/band%dc))
        Umaxi = max(Umaxi, nint(ubm(:,k)/band%dc))
     end where
  end do

  Umini = max(Umini, -band%ngwv)
  Umaxi = min(Umaxi, band%ngwv)

  !-----------------------------------------------------------------------
  ! Gravity wave sources
  !-----------------------------------------------------------------------
  ! Start loop over all columns.
  !-----------------------------------------------------------------------
  do i=1,ncol

     !---------------------------------------------------------------------
     ! Look up spectrum only if the heating depth is large enough, else set
     ! tau0 = 0.
     !---------------------------------------------------------------------

     if (hd_idx(i) > 0) then

        !------------------------------------------------------------------
        ! Look up the spectrum using depth and uh.
        !------------------------------------------------------------------

        tau0 = desc%mfcc(hd_idx(i),nint(uh(i)),:)

        if (desc%storm_shift) then
           ! For deep convection, the wind was relative to storm cells, so
           ! shift the spectrum so that it is now relative to the ground.
           shift = -nint(CS(i)/band%dc)
           tau0 = eoshift(tau0, shift)
        end if

        ! Adjust magnitude.
        tau0 = tau0*(q0(i)**2)/AL

        ! Adjust for critical level filtering.
        tau0(Umini(i):Umaxi(i)) = 0.0
 
        tau(i,:,topi(i)+1) = tau0

     else
      if (present(dqcdt)) then
        if (dqcdt(i,desc%k(i)) > 1.e-8) then ! frontal region (large-scale forcing)
        ! include forced background stress in extra tropical large-scale systems
        ! Set the phase speeds and wave numbers in the direction of the source wind.
        ! Set the source stress magnitude (positive only, note that the sign of the 
        ! stress is the same as (c-u).
         tau(i,:,desc%k(i)+1) = desc%taubck(i,:)
         topi(i) = desc%k(i)
        endif
      endif
     endif
  enddo
  !-----------------------------------------------------------------------
  ! End loop over all columns.
  !-----------------------------------------------------------------------

  ! Output the source level.
  src_level = topi
  tend_level = topi

  !! Set phase speeds; just use reference speeds.
 ! c = spread(band%cref, 1, ncol)

    !debug
 !   if (fms_mpp_pe()==fms_mpp_root_pe()) then
 !      print*, 'maximum heating rates in Beres source after: ', maxval(maxq0)
 !      print*, 'max hdepth after: ', maxval(hdepth)
 !   endif

end subroutine gw_beres_src

!==========================================
! Short routine to get the indices of a set of values rounded to their
! nearest points on a grid.
function index_of_nearest(x, grid) result(idx)
  real,     intent(in) :: x(:)
  real, intent(in) :: grid(:)

  integer :: idx(size(x))

  real :: interfaces(size(grid)-1)
  integer :: i, n

  n = size(grid)
  interfaces = (grid(:n-1) + grid(2:))/2.d0

  idx = 1
  do i = 1, n-1
     where (x > interfaces(i)) idx = i + 1
  end do

end function index_of_nearest

!==========================================
! Take two components of a vector, and find the unit vector components and
! total magnitude.
subroutine get_unit_vector(u, v, u_n, v_n, mag)
  real, intent(in) :: u(:)
  real, intent(in) :: v(:)
  real, intent(out) :: u_n(:)
  real, intent(out) :: v_n(:)
  real, intent(out) :: mag(:)

  integer :: i

  mag = sqrt(u*u + v*v)

  ! Has to be a loop/if instead of a where, because floating point
  ! exceptions can trigger even on a masked divide-by-zero operation
  ! (especially on Intel).
  do i = 1, size(mag)
     if (mag(i) > tiny(u)) then
        u_n(i) = u(i)/mag(i)
        v_n(i) = v(i)/mag(i)
     else
        u_n(i) = 0.0
        v_n(i) = 0.0
     end if
  end do

end subroutine get_unit_vector

!================================================
! Vectorized version of a 2D dot product (since the intrinsic dot_product
! is more suitable for arrays of contiguous vectors).
function dot_2d(u1, v1, u2, v2) result(dot_2d_out)
  real, intent(in) :: u1(:), v1(:)
  real, intent(in) :: u2(:), v2(:)

  real :: dot_2d_out(size(u1))

  dot_2d_out = u1*u2 + v1*v2
end function dot_2d

end module gw_beres_mod
