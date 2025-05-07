module cg_drag_mod

use mpp_mod,                only:  input_nml_file, mpp_get_current_pelist
use mpp_domains_mod,        only:  domain2D, mpp_get_ntile_count
use fms_mod,                only:  fms_init, mpp_pe, mpp_root_pe, mpp_npes, &
                                   check_nml_error,  &
                                   error_mesg,  FATAL, NOTE, &
                                   stdlog, write_version_number
use fms2_io_mod,            only:  FmsNetcdfFile_t, FmsNetcdfDomainFile_t, &
                                   register_restart_field, register_axis, unlimited, &
                                   open_file, read_restart, write_restart, close_file, &
                                   register_field, write_data, get_global_io_domain_indices, &
                                   register_variable_attribute
use time_manager_mod,       only:  time_manager_init, time_type, day_of_year
use diag_manager_mod,       only:  diag_manager_init,   &
                                   register_diag_field, send_data
use constants_mod,          only:  constants_init, PI, RDGAS, GRAV, CP_AIR, &
                                   SECONDS_PER_DAY
use FMSconstants,           only:  RAD_TO_DEG

!-------------------------------------------------------------------

implicit none
private

!---------------------------------------------------------------------
!    cg_drag_mod computes the convective gravity wave forcing on 
!    the zonal flow. the parameterization is described in Alexander and 
!    Dunkerton [JAS, 15 December 1999]. 
!--------------------------------------------------------------------
  

!---------------------------------------------------------------------
!----------- ****** VERSION NUMBER ******* ---------------------------


character(len=128)  :: version =  '$Id$'
character(len=128)  :: tagname =  '$Name$'



!---------------------------------------------------------------------
!-------  interfaces --------

public    cg_drag_init, cg_drag_calc, cg_drag_end, cg_drag_restart, &
          cg_drag_time_vary, cg_drag_endts


private   read_nc_restart_file, gwfc

!--- for netcdf restart
type (domain2D), pointer               :: cg_domain !< Atmosphere domain
integer                                :: old_time_step

!wfc++ Addition for regular use
      integer, allocatable, dimension(:,:)     ::  source_level

      real,     allocatable, dimension(:,:)     ::  source_amp
      real,     allocatable, dimension(:,:,:)   ::  gwd_u, gwd_v
!wfc--


!--------------------------------------------------------------------
!---- namelist -----

integer     :: cg_drag_freq=0     ! calculation frequency [ s ]
integer     :: cg_drag_offset=0   ! offset of calculation from 00Z [ s ]
                                  ! only has use if restarts are written
                                  ! at 00Z and calculations are not done
                                  ! every time step

real        :: source_level_pressure= 315.e+02    
                                  ! highest model level with  pressure 
                                  ! greater than this value (or sigma
                                  ! greater than this value normalized
                                  ! by 1013.25 hPa) will be the gravity
                                  ! wave source level at the equator 
                                  ! [ Pa ]
integer     :: nk=1               ! number of wavelengths contained in 
                                  ! the gravity wave spectrum
real        :: cmax=99.6          ! maximum phase speed in gravity wave
                                  ! spectrum [ m/s ]
real        :: dc=2.4             ! gravity wave spectral resolution 
                                  ! [ m/s ]
                                  ! previous values: 0.6
real        :: Bt_0=.003          ! sum across the wave spectrum of 
                                  ! the magnitude of momentum flux, 
                                  ! divided by density [ m^2/s^2 ]
            
real        :: Bt_aug=.000        ! magnitude of momentum flux divided by density 

real        :: Bt_nh=.003         ! magnitude of momentum flux divided by density   (SH limit )

real        :: Bt_sh=.003         ! magnitude of momentum flux divided by density  (SH limit )

real        :: phi0n = 30., phi0s = -30., dphin = 5., dphis = -5.

logical     :: dump_flux=.true.
                                  ! deposit remaining flux at the model top ?
logical     :: do_conserve_energy=.true.
                                  ! conserve total energy?
real         :: nh_fixer_amp=0. ! ad-hoc fixer on the seasonal cycle of the NH source amplitude

integer      :: nh_fixer_type = 0 ! ad-hoc fixer: 0=not active; 1=dipole; 2=apply to early winter

namelist / cg_drag_nml /         &
                          cg_drag_freq, cg_drag_offset, &
                          source_level_pressure,   &
                          nk, cmax, dc, Bt_0, Bt_aug,  &
                          Bt_sh, Bt_nh,  &
                          phi0n,phi0s,dphin,dphis,      &
                          dump_flux, do_conserve_energy, nh_fixer_amp, &
                          nh_fixer_type

!--------------------------------------------------------------------
!-------- public data  -----


!--------------------------------------------------------------------
!------ private data ------

!--------------------------------------------------------------------
!   these are the arrays which define the gravity wave source spectrum:
!
!   c0       gravity wave phase speeds [ m/s ]
!   kwv      horizontal wavenumbers of gravity waves  [  /m ]
!   k2       squares of wavenumbers [ /(m^2) ]
!
!-------------------------------------------------------------------
real,    dimension(:),     allocatable   :: c0, kwv, k2


!---------------------------------------------------------------------
!   wave spectrum parameters.
!---------------------------------------------------------------------
integer    :: nc        ! number of wave speeds in spectrum
                        ! (symmetric around c = 0)
integer    :: flag = 1  ! flag = 1  for peak flux at  c    = 0
                        ! flag = 0  for peak flux at (c-u) = 0
real       :: Bw = 0.4  ! amplitude for the wide spectrum [ m^2/s^2 ]  
                        ! ~ u'w'
real       :: Bn = 0.0  ! amplitude for the narrow spectrum [ m^2/s^2 ] 
                        ! ~ u'w';  previous values: 5.4
real       :: cw = 40.0 ! half-width for the wide c spectrum [ m/s ]
                        ! previous values: 50.0, 25.0 
real       :: cn =  2.0 ! half-width for the narrow c spectrum  [ m/s ]
integer    :: klevel_of_source
                        ! k index of the gravity wave source level at
                        ! the equator in a standard atmosphere

!---------------------------------------------------------------------
!   variables which control module calculations:
!   
!   cgdrag_alarm time remaining until next cg_drag calculation  [ s ]
!
!---------------------------------------------------------------------
integer          :: cgdrag_alarm

!---------------------------------------------------------------------
!   variables for netcdf diagnostic fields.
!---------------------------------------------------------------------
integer          :: id_kedx_cgwd, id_kedy_cgwd, id_bf_cgwd, &
                    id_gwfx_cgwd, id_gwfy_cgwd
real             :: missing_value = -999.
character(len=7) :: mod_name = 'cg_drag'


logical          :: module_is_initialized=.false.

!-------------------------------------------------------------------
!-------------------------------------------------------------------



                        contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!                      PUBLIC SUBROUTINES
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


!####################################################################

subroutine cg_drag_init (domain, lonb, latb, pref, Time, axes)

!-------------------------------------------------------------------
!   cg_drag_init is the constructor for cg_drag_mod.
!-------------------------------------------------------------------

!-------------------------------------------------------------------
type(domain2D), target,  intent(in)      :: domain !< Atmosphere domain
real,    dimension(:,:), intent(in)      :: lonb, latb
real,    dimension(:),   intent(in)      :: pref
integer, dimension(4),   intent(in)      :: axes
type(time_type),         intent(in)      :: Time
!-------------------------------------------------------------------

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

!-------------------------------------------------------------------
!   local variables: 

      integer                 :: unit, ierr, io, logunit
      integer                 :: n, i, j, k
      integer                 :: idf, jdf, kmax
!      real                    :: pif = 3.14159265358979/180.
!      real                    :: pif = PI/180.

!      real, allocatable       :: lat(:,:)
      real                    :: lat(size(lonb,1) - 1, size(latb,2) - 1)
!-------------------------------------------------------------------
!   local variables: 
!   
!       unit           unit number for nml file 
!       ierr           error return flag 
!       io             error return code 
!       n              loop index
!       k              loop index
!       idf            number of i points on this processor
!       jdf            number of j points on this processor
!       kmax           number of k points on this processor
!
!---------------------------------------------------------------------

!---------------------------------------------------------------------
!    if routine has already been executed, return.
!---------------------------------------------------------------------
      if (module_is_initialized) return

!---------------------------------------------------------------------
!    verify that all modules used by this module have been initialized.
!---------------------------------------------------------------------
      call fms_init
      call time_manager_init
      call diag_manager_init
      call constants_init
!---------------------------------------------------------------------
!    read namelist.
!---------------------------------------------------------------------
      read (input_nml_file, nml=cg_drag_nml, iostat=io)
      ierr = check_nml_error(io,"cg_drag_nml")

!---------------------------------------------------------------------
!    write version number and namelist to logfile.
!---------------------------------------------------------------------
      call write_version_number (version, tagname)
      logunit = stdlog()
      if (mpp_pe() == mpp_root_pe()) write (logunit, nml=cg_drag_nml)

!-------------------------------------------------------------------
!  define the grid dimensions. idf and jdf are the (i,j) dimensions of 
!  domain on this processor, kmax is the number of model layers.
!-------------------------------------------------------------------
      kmax = size(pref(:)) - 1
      jdf  = size(latb,2) - 1
      idf  = size(lonb,1) - 1

      allocate(  source_level(idf,jdf)  )
      allocate(  source_amp(idf,jdf)  )
!      allocate(  lat(idf,jdf)  )

!--------------------------------------------------------------------
!    define the k level which will serve as source level for the grav-
!    ity waves. it is that model level just below the pressure specif-
!    ied as the source location via namelist input.
!--------------------------------------------------------------------
      do k=1,kmax
        if (pref(k) > source_level_pressure) then
          klevel_of_source = k
          exit
        endif
      end do

      do j=1,jdf
        do i=1,idf
          lat(i,j)=  0.5*( latb(i,j+1)+latb(i,j) )
          source_level(i,j) = (kmax + 1) - ((kmax + 1 -    &
                              klevel_of_source)*cos(lat(i,j)) + 0.5)
          source_amp(i,j) = Bt_0 +                         &
                      Bt_nh*0.5*(1.+tanh((lat(i,j)*RAD_TO_DEG-phi0n)/dphin)) + &
                      Bt_sh*0.5*(1.+tanh((lat(i,j)*RAD_TO_DEG-phi0s)/dphis))
        end do
      end do
      source_level = MIN (source_level, kmax-1)

!      deallocate( lat )

!---------------------------------------------------------------------
!    define the number of waves in the gravity wave spectrum, and define
!    an array of their speeds. They are defined symmetrically around
!    c = 0.0 m/s.
!---------------------------------------------------------------------
      nc = 2.0*cmax/dc + 1
      allocate ( c0(nc) )
      do n=1,nc
        c0(n) = (n-1)*dc - cmax
      end do
 
!--------------------------------------------------------------------
!    define the wavenumber kwv and its square k2 for the gravity waves 
!    contained in the spectrum. currently nk = 1, which means that the 
!    wavelength of all gravity waves considered is 300 km. 
!--------------------------------------------------------------------
      allocate ( kwv(nk) )
      allocate ( k2 (nk) )
      do n=1,nk
        kwv(n) = 2.*PI/((30.*(10.**n))*1.e3)
        k2(n) = kwv(n)*kwv(n)
      end do

!--------------------------------------------------------------------
!    initialize netcdf diagnostic fields.
!-------------------------------------------------------------------
      id_bf_cgwd =  &
         register_diag_field (mod_name, 'bf_cgwd', axes(1:3), Time, &
              'buoyancy frequency from cg_drag', ' /s',   &
              missing_value=missing_value)
      id_gwfx_cgwd =  &
         register_diag_field (mod_name, 'gwfx_cgwd', axes(1:3), Time, &
              'gravity wave forcing on mean flow', &
              'm/s^2',  missing_value=missing_value)
      id_gwfy_cgwd =  &
         register_diag_field (mod_name, 'gwfy_cgwd', axes(1:3), Time, &
              'gravity wave forcing on mean flow', &
              'm/s^2',  missing_value=missing_value)
      id_kedx_cgwd =  &
         register_diag_field (mod_name, 'kedx_cgwd', axes(1:3), Time, &
               'effective eddy viscosity from cg_drag', 'm^2/s',   &
               missing_value=missing_value)
      id_kedy_cgwd =  &
         register_diag_field (mod_name, 'kedy_cgwd', axes(1:3), Time, &
               'effective eddy viscosity from cg_drag', 'm^2/s',   &
               missing_value=missing_value)

!--------------------------------------------------------------------
!    allocate and define module variables to hold values across 
!    timesteps, in the event that cg_drag is not called on every step.
!--------------------------------------------------------------------
     allocate ( gwd_u(idf,jdf,kmax) )
     allocate ( gwd_v(idf,jdf,kmax) )

!--------------------------------------------------------------------
!    if present, read the restart data file.
!---------------------------------------------------------------------

      cg_domain => domain
      call read_nc_restart_file

      old_time_step = cgdrag_alarm 
!---------------------------------------------------------------------
!    mark the module as initialized.
!---------------------------------------------------------------------
      module_is_initialized = .true.

!---------------------------------------------------------------------



end subroutine cg_drag_init


!####################################################################
 
subroutine cg_drag_time_vary (delt)

real           ,        intent(in)      :: delt

!---------------------------------------------------------------------
!    decrement the time remaining until the next cg_drag calculation.
!---------------------------------------------------------------------
      cgdrag_alarm = cgdrag_alarm - delt

!---------------------------------------------------------------------
 
end subroutine cg_drag_time_vary


!####################################################################
 
subroutine cg_drag_endts
 
!--------------------------------------------------------------------
!    if this was a calculation step, reset cgdrag_alarm to indicate 
!    the time remaining before the next calculation of gravity wave 
!    forcing.
!--------------------------------------------------------------------
      if (cgdrag_alarm <= 0 ) then
        cgdrag_alarm = cgdrag_alarm + cg_drag_freq
      endif

end subroutine cg_drag_endts


!####################################################################

subroutine cg_drag_calc (is, js, lat, pfull, zfull, temp, uuu, vvv,  &
                         Time, delt, gwfcng_x, gwfcng_y, dtemp)
!--------------------------------------------------------------------  
!    cg_drag_calc defines the arrays needed to calculate the convective
!    gravity wave forcing, calls gwfc to calculate the forcing, returns 
!    the desired output fields, and saves the values for later retrieval
!    if they are not calculated on every timestep.
!
!---------------------------------------------------------------------

!---------------------------------------------------------------------
integer,                intent(in)      :: is, js
real, dimension(:,:),   intent(in)      :: lat
real, dimension(:,:,:), intent(in)      :: pfull, zfull, temp, uuu, vvv
type(time_type),        intent(in)      :: Time
real           ,        intent(in)      :: delt
real, dimension(:,:,:), intent(out)     :: gwfcng_x, gwfcng_y
real, dimension(:,:,:), intent(out)     :: dtemp 

!-------------------------------------------------------------------
!    intent(in) variables:
!
!       is,js    starting subdomain i,j indices of data in 
!                the physics_window being integrated
!       lat      array of model latitudes at cell boundaries [radians]
!       pfull    pressure at model full levels [ Pa ]
!       zfull    height at model full levels [ m ]
!       temp     temperature at model levels [ deg K ]
!       uuu      zonal wind  [ m/s ]
!       vvv      meridional wind  [ m/s ]
!       Time     current time, needed for diagnostics [ time_type ]
!       delt     physics time step [ s ]
!
!    intent(out) variables:
!
!       gwfcng_x time tendency for u eqn due to gravity-wave forcing
!                [ m/s^2 ]
!       gwfcng_y time tendency for v eqn due to gravity-wave forcing
!                [ m/s^2 ]
!       dtemp    time tendency of the temperature in K/s due to gravity-wave forcing
!                [ K/s ]
!
!-------------------------------------------------------------------

!-------------------------------------------------------------------
!    local variables:

      real,    dimension (size(uuu,1), size(uuu,2), size(uuu,3))  ::  &
                                         dtdz, ked_gwfc_x, ked_gwfc_y

      real,    dimension (size(uuu,1),size(uuu,2), 0:size(uuu,3)) ::  &
                                         zzchm, zu, zv, zden, zbf,    &
                                         gwd_xtnd, ked_xtnd,          &
                                         gwd_ytnd, ked_ytnd

      integer           :: iz0
      logical           :: used
      real              :: bflim = 2.5E-5
      integer           :: ie, je
      integer           :: imax, jmax, kmax
      integer           :: i, j, k
      integer           :: nday ! the number of day in year used for ad-hoc fixer
      real              :: fixer_amp ! amplitude of the ad-hoc fixer
      real, dimension(size(uuu,1),size(uuu,2)) :: fa_lat ! ad-hoc fixer amplitude 
                                 ! with latitudinal variance

!-------------------------------------------------------------------
!    local variables:
!
!       dtdz          temperature lapse rate [ deg K/m ]
!       ked_gwfc      effective diffusion coefficient from cg_drag_mod 
!                     [ m^2/s ]
!       zzchm         heights at model levels [ m ]
!       zu            zonal velocity [ m/s ]
!       zden          atmospheric density [ kg/m^3 ]
!       zbf           buoyancy frequency [ /s ]
!       gwd_xtnd      zonal wind tendency resulting from cg_drag_mod 
!                     [ m/s^2 ]
!       ked_xtnd      effective diffusion coefficient from cg_drag_mod 
!                     [ m^2/s ]
!       source_level  k index of gravity wave source level ((i,j) array)
!       iz0           k index of gravity wave source level in a column
!       used          return code for netcdf diagnostics
!       bflim         minimum allowable value of squared buoyancy 
!                     frequency [ /s^2 ]
!       ie, je        ending subdomain indices of data in the current 
!                     physics window being integrated
!       imax, jmax, kmax 
!                     physics window dimensions
!       i, j, k, nn   do loop indices
!
!---------------------------------------------------------------------

!---------------------------------------------------------------------
!    define processor extents and loop limits.
!---------------------------------------------------------------------
      imax = size(uuu,1)
      jmax = size(uuu,2)
      kmax = size(uuu,3)
      ie = is + imax - 1
      je = js + jmax - 1

!---------------------------------------------------------------------
!    if the convective gravity wave forcing should be calculated on 
!    this timestep (i.e., the alarm has gone off), proceed with the
!    calculation.
!---------------------------------------------------------------------


      if (cgdrag_alarm <= 0) then

!-----------------------------------------------------------------------
!    calculate temperature lapse rate. do one-sided differences over 
!    delta z at upper boundary and centered differences over 2 delta z 
!    in the interior.  dtdz is not needed at the lower boundary, since
!    the source level is constrained to be above level kmax.
!----------------------------------------------------------------------
        do j=1,jmax
          do i=1,imax
! The following index-offsets are needed in case a physics_window is being used.
            iz0 = source_level(i +is-1,j+js-1)
            dtdz(i,j,1) = (temp  (i,j,1) - temp  (i,j,2))/    &
                          (zfull(i,j,1) - zfull(i,j,2))
            do k=2,iz0
              dtdz(i,j,k) = (temp  (i,j,k-1) - temp  (i,j,k+1))/   &
                            (zfull(i,j,k-1) - zfull(i,j,k+1))
            end do

!--------------------------------------------------------------------
!    calculate air density.
!--------------------------------------------------------------------
            do k=1,iz0+1
              zden(i,j,k  ) = pfull(i,j,k)/(temp(i,j,k)*RDGAS)
            end do

!----------------------------------------------------------------------
!    calculate buoyancy frequency. restrict the squared buoyancy 
!    frequency to be no smaller than bflim.
!----------------------------------------------------------------------
            do k=1,iz0 
              zbf(i,j,k) = (GRAV/temp(i,j,k))*(dtdz(i,j,k) + GRAV/CP_AIR)
              if (zbf(i,j,k) < bflim) then
                zbf(i,j,k) = sqrt(bflim)
              else 
                zbf(i,j,k) = sqrt(zbf(i,j,k))
              endif
            end do

!----------------------------------------------------------------------
!    if zbf is to be saved for netcdf output, the remaining vertical
!    levels must be initialized.
!----------------------------------------------------------------------
            if (id_bf_cgwd > 0) then
              zbf(i,j,iz0+1:) = 0.0
            endif

!----------------------------------------------------------------------
!    define an array of heights at model levels and an array containing
!    the zonal wind component.
!----------------------------------------------------------------------
            do k=1,iz0+1
              zzchm(i,j,k) = zfull(i,j,k)
            end do
            do k=1,iz0   
              zu(i,j,k) = uuu(i,j,k)
              zv(i,j,k) = vvv(i,j,k)
            end do

!----------------------------------------------------------------------
!    add an extra level above model top so that the gravity wave forcing
!    occurring between the topmost model level and the upper boundary
!    may be calculated. define variable values at the new top level as
!    follows: z - use delta z of layer just below; u - extend vertical 
!    gradient occurring just below; density - geometric mean; buoyancy 
!    frequency - constant across model top.
!----------------------------------------------------------------------
            zzchm(i,j,0) = zzchm(i,j,1) + zzchm(i,j,1) - zzchm(i,j,2)
            zu(i,j,0)    = 2.*zu(i,j,1) - zu(i,j,2)
            zv(i,j,0)    = 2.*zv(i,j,1) - zv(i,j,2)
            zden(i,j,0)  = zden(i,j,1)*zden(i,j,1)/zden(i,j,2)
            zbf(i,j,0)   = zbf(i,j,1)
          end do
        end do

 !----------------------------------------------------------------------
 !   Ad-hoc fixer applied to NH wave source.
 !   The seasonal dependence of the wave source amplitude is controlled by nh_fixer type.
 !   Only apply to NH
 !----------------------------------------------------------------------
     if (nh_fixer_type>0 ) then
     nday=day_of_year(Time)

     ! nh_fixer_tyep=1: 
     ! reduce wave source amplitude in Oct/Nov and increase in Dec/Jan
     if (nh_fixer_type==1) then 
     if (nday>152) then
         fixer_amp= nh_fixer_amp*((nday-335.)/40.)*exp(-((nday-335.)/40.)**2)
     else
         fixer_amp= nh_fixer_amp*((nday+30.)/40.)*exp(-((nday+30.)/40.)**2)
     endif
     ! nh_fixer_type=2:
     ! reduce wave amplitude in Oct/Nov following a Gaussian function
     elseif (nh_fixer_type==2) then
         fixer_amp = - nh_fixer_amp*exp(-((nday-300.)/20.)**2)
     endif ! nh_fixer_type

     do i = 1,imax
         do j = 1,jmax
             fa_lat(i,j)=(1.+fixer_amp* &
             0.5*(1.+tanh((lat(i,j)*RAD_TO_DEG-phi0n)/dphin)))
         end do
     end do
     endif

!---------------------------------------------------------------------
!    pass the vertically-extended input arrays to gwfc. gwfc will cal-
!    culate the gravity-wave forcing and, if desired, an effective eddy 
!    diffusion coefficient at each level above the source level. output
!    is returned in the vertically-extended arrays gwfcng and ked_gwfc.
!    upon return move the output fields into model-sized arrays. 
!---------------------------------------------------------------------
     if (nh_fixer_type>0) then ! use ad-hoc fixer
       call gwfc (is, ie, js, je, source_level, source_amp,    &
              zden, zu, zbf, zzchm, gwd_xtnd, ked_xtnd, fixer=fa_lat)
     else ! ad-hoc fixer not activated
       call gwfc (is, ie, js, je, source_level, source_amp,    &
                     zden, zu, zbf,zzchm, gwd_xtnd, ked_xtnd)
     endif
          gwfcng_x  (:,:,1:kmax) = gwd_xtnd(:,:,1:kmax  )
          ked_gwfc_x(:,:,1:kmax) = ked_xtnd(:,:,1:kmax  )

     if (abs(nh_fixer_amp)>0.01) then ! use ad-hoc fixer
       call gwfc (is, ie, js, je, source_level, source_amp,    &
              zden, zv, zbf, zzchm, gwd_ytnd, ked_ytnd, fixer=fa_lat)
     else ! ad-hoc fixer not activated
       call gwfc (is, ie, js, je, source_level, source_amp,    &
                     zden, zv, zbf,zzchm, gwd_ytnd, ked_ytnd)
     endif
          gwfcng_y  (:,:,1:kmax) = gwd_ytnd(:,:,1:kmax  )
          ked_gwfc_y(:,:,1:kmax) = ked_ytnd(:,:,1:kmax  )

!--------------------------------------------------------------------
!    store the gravity wave forcing into a processor-global array.
!-------------------------------------------------------------------
          gwd_u(is:ie,js:je,:) = gwfcng_x(:,:,:)
          gwd_v(is:ie,js:je,:) = gwfcng_y(:,:,:)


!--------------------------------------------------------------------
!    if activated, store the effective eddy diffusivity into a 
!    processor-global array, and if desired as a netcdf diagnostic, 
!    send the data to diag_manager_mod.
!-------------------------------------------------------------------

          if (id_kedx_cgwd > 0) then
            used = send_data (id_kedx_cgwd, ked_gwfc_x, Time, is, js, 1)
          endif

          if (id_kedy_cgwd > 0) then
            used = send_data (id_kedy_cgwd, ked_gwfc_y, Time, is, js, 1)
          endif



!--------------------------------------------------------------------
!    save any other netcdf file diagnostics that are desired.
!--------------------------------------------------------------------
        if (id_bf_cgwd > 0) then
          used = send_data (id_bf_cgwd,  zbf(:,:,1:), Time, is, js )
        endif

        if (id_gwfx_cgwd > 0) then
          used = send_data (id_gwfx_cgwd, gwfcng_x, Time, is, js, 1)
        endif
        if (id_gwfy_cgwd > 0) then
          used = send_data (id_gwfy_cgwd, gwfcng_y, Time, is, js, 1)
        endif



!--------------------------------------------------------------------
!    if this is not a timestep on which gravity wave forcing is to be 
!    calculated, retrieve the values calculated previously from storage
!    and return to the calling subroutine.
!--------------------------------------------------------------------
      else   ! (cgdrag_alarm <= 0)
        gwfcng_x(:,:,:) = gwd_u(is:ie,js:je,:)
        gwfcng_y(:,:,:) = gwd_v(is:ie,js:je,:)
     endif  ! (cgdrag_alarm <= 0)

!--------------------------------------------------------------------



! CALCULATE HEATING TO CONSERVE TOTAL ENERGY

  if (do_conserve_energy) then
     dtemp = -((uuu + 0.5*delt*gwfcng_x)*gwfcng_x                           &
             + (vvv + 0.5*delt*gwfcng_y)*gwfcng_y)/Cp_Air
  else
     dtemp = 0.0
  endif

end subroutine cg_drag_calc



!###################################################################

subroutine cg_drag_end

!--------------------------------------------------------------------
!    cg_drag_end is the destructor for cg_drag_mod.
!--------------------------------------------------------------------

!--------------------------------------------------------------------
!    local variables

!For version 3 and after, use NetCDF restarts.
      if (mpp_pe() == mpp_root_pe() ) &
            call error_mesg ('cg_drag_mod', 'write_restart_nc: &
              &Writing netCDF formatted restart file as &
                &requested. ', NOTE)
      call cg_drag_restart


!---------------------------------------------------------------------
!    mark the module as uninitialized.
!---------------------------------------------------------------------
      module_is_initialized = .false.

!---------------------------------------------------------------------


end subroutine cg_drag_end



!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!                     PRIVATE SUBROUTINES
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


subroutine read_nc_restart_file
!-----------------------------------------------------------------------
!    subroutine read_restart_nc reads a netcdf restart file to obtain 
!    the variables needed upon experiment restart. 
!-----------------------------------------------------------------------

!---------------------------------------------------------------------
!   local variables:

      character(len=64)     :: fname='INPUT/cg_drag.res.nc'
      type(FmsNetcdfFile_t)       :: Cg_restart !< Fms2io fileobj
      type(FmsNetcdfDomainFile_t) :: Til_restart !< Fms2io domain fileobj
      integer, allocatable, dimension(:) :: pes !< Array of the pes in the current pelist

!---------------------------------------------------------------------
!   local variables:
!
!        fname            restart file name
!
!----------------------------------------------------------------------

!--------------------------------------------------------------------
!    output a message indicating entrance into this routine.
!--------------------------------------------------------------------
      if (mpp_pe() == mpp_root_pe() ) then
        call error_mesg ('cg_drag_mod',  'read_restart_nc:&
             &Reading netCDF formatted restart file:'//trim(fname), NOTE)
      endif

      !< Get the current pelist
      allocate(pes(mpp_npes()))
      call mpp_get_current_pelist(pes)

      !< Open the scalar file with the current pelist, so that only the root pe opens and reads the file and
      !! distributes the data to the other pes
      !> read the values of gwd_u and gwd_v
      if (open_file(Cg_restart, fname, "read", is_restart = .true., pelist=pes)) then
        call cg_drag_register_restart(Cg_restart)
        call read_restart(Cg_restart)
        call close_file(Cg_restart)

        !> if current cg_drag calling frequency differs from that previously
        !! used, adjust the time remaining before the next calculation.
        if (cg_drag_freq /= old_time_step) then
          cgdrag_alarm = cgdrag_alarm - old_time_step + cg_drag_freq
          if (mpp_pe() == mpp_root_pe() ) then
            call error_mesg ('cg_drag_mod',   &
                  'cgdrag time step has changed, &
                  &next cgdrag time also changed', NOTE)
          endif
          old_time_step = cg_drag_freq
        endif

        !> if cg_drag_offset is specified and is smaller than the time remaining
        !! until the next calculation, modify the time remaining to be
        !! that offset time. the assumption is made that the restart was
        !! written at 00Z.
        if (cg_drag_offset /= 0) then
          if (cgdrag_alarm > cg_drag_offset) then
            cgdrag_alarm = cg_drag_offset
          endif
        endif
      else

        !> if no restart file is present, initialize the gwd field to zero.
        !! define the time remaining until the next cg_drag calculation from
        !! the namelist inputs.
        if (cg_drag_offset > 0) then
          cgdrag_alarm = cg_drag_offset
        else
          cgdrag_alarm = cg_drag_freq
        endif
      endif
      deallocate(pes)

      if (open_file(Til_restart, fname, "read", cg_domain, is_restart = .true.)) then
        call cg_drag_register_tile_restart(Til_restart)
        call read_restart(Til_restart)
        call close_file(Til_restart)
      else
        !> if no restart file is present, initialize the gwd field to zero.
        !! define the time remaining until the next cg_drag calculation from
        !! the namelist inputs.
        gwd_u(:,:,:) = 0.0
        gwd_v(:,:,:) = 0.0
      endif

!---------------------------------------------------------------------
end subroutine read_nc_restart_file

!####################################################################
! register restart field to be read and written through save_restart and restore_state.
subroutine cg_drag_register_restart(Cg_restart)

  type(FmsNetcdfFile_t), intent(inout) :: Cg_restart !< Fms2io file obj
  character(len=8), dimension(1)       :: dim_names !< Array of dimension names

  dim_names(1) = "Time"
  call register_axis(Cg_restart, dim_names(1), unlimited)
  call register_restart_field(Cg_restart, "cgdrag_alarm", cgdrag_alarm, dim_names)
  call register_restart_field(Cg_restart, "cg_drag_freq", old_time_step, dim_names)

end subroutine cg_drag_register_restart

!> \brief register restart field to be read and written through save_restart and restore_state.
subroutine cg_drag_register_tile_restart (Til_restart)

  type(FmsNetcdfDomainFile_t), intent(inout) :: Til_restart !< Fms2io domain file obj
  character(len=8), dimension(4)             :: dim_names !< Array of dimension names

  dim_names(1) = "xaxis_1"
  dim_names(2) = "yaxis_1"
  dim_names(3) = "zaxis_1"
  dim_names(4) = "Time"

  call register_axis(Til_restart, dim_names(1), "x")
  call register_axis(Til_restart, dim_names(2), "y")
  call register_axis(Til_restart, dim_names(3), size(gwd_u, 3))
  if (.not. Til_restart%mode_is_append)  call register_axis(Til_restart, dim_names(4), unlimited)

  !< Register the domain decomposed dimensions as variables so that the combiner can work
  !! correctly
  call register_field(Til_restart, dim_names(1), "double", (/dim_names(1)/))
  call register_field(Til_restart, dim_names(2), "double", (/dim_names(2)/))

  call register_restart_field(Til_restart, "gwd_u", gwd_u, dim_names)
  call register_restart_field(Til_restart, "gwd_v", gwd_v, dim_names)

end subroutine cg_drag_register_tile_restart

!####################################################################
subroutine cg_drag_restart(timestamp)

! write out restart file.

  character(len=*), intent(in), optional :: timestamp !< A character string that represents the model time,
                                                      !! used for writing restart. timestamp will append to
                                                      !! the any restart file name as a prefix.

  character(len=128)           :: fname
  type(FmsNetcdfFile_t)       :: Cg_restart
  type(FmsNetcdfDomainFile_t) :: Til_restart
  logical :: tile_file_exist
  integer, allocatable, dimension(:) :: pes

  if (present(timestamp)) then
    fname='RESTART/'//trim(timestamp)//'.cg_drag.res.nc'
  else
    fname='RESTART/cg_drag.res.nc'
  endif

  !< Get the current pelist
  allocate(pes(mpp_npes()))
  call mpp_get_current_pelist(pes)

  !< Open the scalar file with the current pelist, so that only the root pe opens and writes the file
  if (open_file(Cg_restart, fname, "overwrite", is_restart = .true., pelist=pes)) then
     call cg_drag_register_restart(Cg_restart)
     call write_restart(Cg_restart)
     call close_file(Cg_restart)
  endif
  deallocate(pes)

  if (mpp_get_ntile_count(cg_domain) == 1) then
     tile_file_exist = open_file(Til_restart, fname, "append", cg_domain, is_restart = .true.)
  else
     tile_file_exist = open_file(Til_restart, fname, "overwrite", cg_domain, is_restart = .true.)
  endif

  if (tile_file_exist) then
     call cg_drag_register_tile_restart(Til_restart)
     call write_restart(Til_restart)
     call add_domain_dimension_data(Til_restart)
     call close_file(Til_restart)
  endif

end subroutine cg_drag_restart

!< Add_dimension_data: Adds dummy data for the domain decomposed axis
subroutine add_domain_dimension_data(fileobj)
  type(FmsNetcdfDomainFile_t) :: fileobj !< Fms2io domain decomposed fileobj
  integer, dimension(:), allocatable :: buffer !< Buffer with axis data
  integer :: is, ie !< Starting and Ending indices for data

    call get_global_io_domain_indices(fileobj, "xaxis_1", is, ie, indices=buffer)
    call write_data(fileobj, "xaxis_1", buffer)
    deallocate(buffer)

    call get_global_io_domain_indices(fileobj, "yaxis_1", is, ie, indices=buffer)
    call write_data(fileobj, "yaxis_1", buffer)
    deallocate(buffer)

end subroutine add_domain_dimension_data

!####################################################################

subroutine gwfc (is, ie, js, je, source_level, source_amp, rho, u,    &
                 bf, z, gwf, ked, fixer)

!-------------------------------------------------------------------
!    subroutine gwfc computes the gravity wave-driven-forcing on the
!    zonal wind given vertical profiles of wind, density, and buoyancy 
!    frequency. 
!    Based on version implemented in SKYHI -- 27 Oct 1998 by M.J. 
!    Alexander and L. Bruhwiler.
!-------------------------------------------------------------------

!-------------------------------------------------------------------
integer,                     intent(in)             :: is, ie, js, je
integer, dimension(:,:),     intent(in)             :: source_level
real,    dimension(:,:),     intent(in)             :: source_amp
real,    dimension(:,:,0:),  intent(in)             :: rho, u, bf, z
real,    dimension(:,:,0:),  intent(out)            :: gwf
real,    dimension(:,:,0:),  intent(out)            :: ked
real,    dimension(:,:),     intent(in), optional   :: fixer

!-------------------------------------------------------------------
!  intent(in) variables:
!
!      is, ie, js, je   starting/ending subdomain i,j indices of data
!                       in the physics_window being integrated
!      source_level     k index of model level serving as gravity wave
!                       source
!      source_amp     amplitude of  gravity wave source
!                       
!      rho              atmospheric density [ kg/m^3 ] 
!      u                zonal wind component [ m/s ]
!      bf               buoyancy frequency [ /s ]
!      z                height of model levels  [ m ]
!
!  intent(out) variables:
!
!      gwf              gravity wave forcing in u equation  [ m/s^2 ]
!
!  intent(out), optional variables:
!
!      ked              eddy diffusion coefficient from gravity wave 
!                       forcing [ m^2/s ]
!
!------------------------------------------------------------------

!------------------------------------------------------------------
!  local variables

      real,    dimension (0:size(u,3)-1 ) :: wv_frcng, diff_coeff
      logical, dimension (nc) ::   msk
      real   , dimension (nc) ::   c0mu0, B0
      integer                 ::   iz0, i, j, k, ink, n
      real                    ::   fm, fe, Hb, alp2, Foc, c, test, rbh,&
                                   c0mu, eps, Bsum, Bexp, ampl,        &
                                   size_u1, size_u2, dz, fac, omc

      size_u1 = size(u,1)
      size_u2 = size(u,2)

!------------------------------------------------------------------
!  local variables:
! 
!      wv_frcng    gravity wave forcing tendency [ m/s^2 ]
!      diff_coeff  eddy diffusion coefficient [ m2/s ]
!      c0mu        difference between phase speed of wave n and u 
!                  [ m/s ]
!      dz          delta z between model levels [ m ]
!      fac         factor used in determining if wave is breaking 
!                  [ s/m ]
!      omc         critical frequency that marks total internal 
!                  reflection  [ /s ]
!      msk         indicator as to whether wave n is still propagating 
!                  upwards (msk=1), or has been removed from the 
!                  spectrum because of breaking or reflection (msk=0)
!      c0mu0       difference between phase speed of wave n and u at the
!                  source level [ m/s ]
!      B0          wave momentum flux amplitude for wave n [ (m/s)^2 ]
!      fm          used to sum up momentum flux from all waves n 
!                  deposited at a level [ (m/s)^2 ]
!      fe          used to sum up contributions to diffusion coefficient
!                  from all waves n at a level [ (m/s)^3 ]
!      Hb          density scale height [ m ]
!      alp2        scale height factor: 1/(2*Hb)**2  [ /m^2 ]
!      Foc         wave breaking threshold [ s/m ]
!      c           wave phase speed used in defining wave momentum flux
!                  amplitude [ m/s ]
!      test        condition defining internal reflection [ /s ]
!      rbh         atmospheric density at half-level (geometric mean)
!                  [ kg/m^3 ]
!      eps         intermittency factor
!      Bsum        total mag of gravity wave momentum flux at source 
!                  level, divided by the density  [ m^2/s^2 ]
!      iz0         source level vertical index for the given column
!      i,j,k       spatial do loop indices
!      ink         wavenumber loop index 
!      n           phase speed loop index 
!      ampl        phase speed loop index 
!
!--------------------------------------------------------------------

!-------------------------------------------------------------------
!    initialize the output arrays. these will hold values at each 
!    (i,j,k) point, summed over the wavelengths and phase speeds
!    defining the gravity wave spectrum.
!-------------------------------------------------------------------
      gwf = 0.0
      ked = 0.0

      do j=1,size_u2
        do i=1,size_u1
! The following index-offsets are needed in case a physics_window is being used.
          iz0 = source_level(i+is-1,j+js-1)
          if (present(fixer)) then
             ampl=source_amp(i+is-1,j+js-1)*fixer(i,j)
          else
          ampl= source_amp(i+is-1,j+js-1)
          endif
!--------------------------------------------------------------------
!    define wave momentum flux (B0) at source level for each phase 
!    speed n, and the sum over all phase speeds (Bsum), which is needed 
!    to calculate the intermittency. 
!-------------------------------------------------------------------
          Bsum = 0.
          do n=1,nc
            c0mu0(n) = c0(n) - u(i,j,iz0)   

!---------------------------------------------------------------------
!    when the wave phase speed is same as wind speed, there is no
!    momentum flux.
!---------------------------------------------------------------------
            if (c0mu0(n) == 0.0)  then
              B0(n) = 0.0
            else 

!---------------------------------------------------------------------
!    define wave momentum flux at source level for phase speed n. Add
!    the contribution from this phase speed to the previous sum.
!---------------------------------------------------------------------
              c = c0(n)*flag + c0mu0(n)*(1 - flag)
              Bexp = exp(-alog(2.0)*(c/cw)**2)
              if (c0mu0(n) < 0.0) then
                B0(n) = -(Bw*Bexp + Bn*Bexp)
              else 
                B0(n) = (Bw*Bexp + Bn*Bexp)
              endif
              Bsum = Bsum + abs(B0(n))
            endif
          end do

!---------------------------------------------------------------------
!    define the intermittency factor eps. the factor of 1.5 is currently
!    unexplained.
!---------------------------------------------------------------------
          if (Bsum == 0.0) then
            call error_mesg ('cg_drag_mod', &
               ' zero flux input at source level', FATAL)
          endif
          eps = (ampl*1.5/nk)/Bsum

!--------------------------------------------------------------------
!    loop over the nk different wavelengths in the spectrum.
!--------------------------------------------------------------------
          do ink=1,nk   ! wavelength loop

!---------------------------------------------------------------------
!    initialize a flag which will indicate which waves are still 
!    propagating upwards.
!---------------------------------------------------------------------
            msk = .true.

!----------------------------------------------------------------------
!    integrate upwards from the source level.  define variables over 
!    which to sum the deposited flux and effective eddy diffusivity 
!    from all waves breaking at a given level.
!----------------------------------------------------------------------
            do k=iz0, 0, -1
!----------------------------------------------------------------------
!    define variables needed at levels above the source level.
!---------------------------------------------------------------------
              fac = 0.5*(rho(i,j,k)/rho(i,j,iz0))*kwv(ink)/bf(i,j,k)
              dz = z(i,j,k) - z(i,j,k+1)
              Hb = -dz/alog(rho(i,j,k)/rho(i,j,k+1))
              alp2 = 0.25/(Hb*Hb)
              omc = sqrt((bf(i,j,k)*bf(i,j,k)*k2(ink))/    &
                            (k2(ink) + alp2))

              fm = 0.
              fe = 0.
              do n=1,nc     ! phase speed loop

!----------------------------------------------------------------------
!    check only those waves which are still propagating, i.e., msk = .true.
!----------------------------------------------------------------------
                if (msk(n)) then
                  c0mu = c0(n) - u(i,j,k)

!----------------------------------------------------------------------
!    if phase speed matches the wind speed, remove c0(n) from the 
!    set of propagating waves.
!----------------------------------------------------------------------
                  if (c0mu == 0.) then
                    msk(n) = .false.
                  else

!---------------------------------------------------------------------
!    define the criterion which determines if wave is reflected at this 
!    level (test).
!---------------------------------------------------------------------
                    test = abs(c0mu)*kwv(ink) - omc
                    if (test >= 0.0) then

!---------------------------------------------------------------------
!    wave has undergone total internal reflection. remove it from the
!    propagating set.
!---------------------------------------------------------------------
                      msk(n) = .false.
                    else 

!---------------------------------------------------------------------
!    if wave is  not reflected at this level, determine if it is 
!    breaking at this level (Foc >= 0),  or if wave speed relative to 
!    windspeed has changed sign from its value at the source level 
!    (c0mu0(n)*c0mu <= 0). if it is above the source level and is
!    breaking, then add its momentum flux to the accumulated sum at 
!    this level, and increase the effective diffusivity accordingly. 
!    set flag to remove phase speed c0(n) from the set of active waves
!    moving upwards to the next level.
!---------------------------------------------------------------------
                      if(c0mu0(n)*c0mu <= 0.0 ) then
                        msk(n) = .false.
                        if (k  < iz0) then
                          fm = fm + B0(n)
                          fe = fe + c0mu*B0(n)
                        endif
                      else
                        Foc = B0(n)/(c0mu)**3 - fac
                        if (Foc >= 0.0) then
                          msk(n) = .false.
                          if (k  < iz0) then
                            fm = fm + B0(n)
                            fe = fe + c0mu*B0(n)
                          endif
                        endif
                      endif
                    endif   ! (test >= 0.0)
                  endif ! (c0mu == 0.0)
                endif   ! (msk == .true.)
              end do  ! phase speed loop


    if( dump_flux ) then 
!   optional upper boundary
!           Dump remaining flux in the top model level. 
             if ( k == 0  )  then
                do n= 1, nc
                    if ( msk(n))  then
                          fm= fm + B0(n)
                          fe = fe + c0mu*B0(n)
                    endif
                enddo 
             endif 
    endif 



!----------------------------------------------------------------------
!    compute the gravity wave momentum flux forcing and eddy 
!    diffusion coefficient obtained across the entire wave spectrum
!    at this level.
!----------------------------------------------------------------------
              if ( k < iz0) then
                rbh = sqrt(rho(i,j,k)*rho(i,j,k+1))
                wv_frcng(k) = ( rho(i,j,iz0)/rbh)*fm*eps/dz
                wv_frcng(k+1) =  0.5*(wv_frcng(k+1) + wv_frcng(k))
                diff_coeff(k) = (rho(i,j,iz0)/rbh)*fe*eps/(dz*   &
                                 bf(i,j,k)*bf(i,j,k))
                diff_coeff(k+1) = 0.5*(diff_coeff(k+1) +    &
                                       diff_coeff(k))
              else 
                wv_frcng(iz0) = 0.0
                diff_coeff(iz0) = 0.0
              endif
            end do  ! (k loop)   
       

!---------------------------------------------------------------------
!    increment the total forcing at each point with that obtained from
!    the set of waves with the current wavenumber.
!---------------------------------------------------------------------
            do k=0,iz0      
              gwf(i,j,k) = gwf(i,j,k) + wv_frcng(k)
              ked(i,j,k) = ked(i,j,k) + diff_coeff(k)
            end do  



          end do   ! wavelength loop
        end do  ! i loop                      
      end do   ! j loop                 

!--------------------------------------------------------------------



end subroutine gwfc



!####################################################################


end module cg_drag_mod
