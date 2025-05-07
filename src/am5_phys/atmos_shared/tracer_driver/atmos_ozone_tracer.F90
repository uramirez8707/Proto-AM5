module atmos_ozone_tracer_mod
! <CONTACT EMAIL="Pu.Lin@noaa.gov">
!   Pu Lin
! </CONTACT>

! <REVIEWER EMAIL="">
!  
! </REVIEWER>

! <HISTORY>

! <OVERVIEW>
!     This code implements an ozone tracer. The chemical tendency of ozone is 
!     replaced by the difference between production and loss. The production 
!     rate is prescribed, the loss rate is calculated as concentration devided
!     by life time. Life time (or -1/relaxation coefficient) is prescribed. 
!     Use implicit scheme for the differentiation
! </OVERVIEW>

! <DESCRIPTION>
!     This code implements an ozone tracer, no ozone chemistry.
! </DESCRIPTION>

!-----------------------------------------------------------------------

use mpp_mod,             only: input_nml_file
use              fms_mod, only : write_version_number, &
                                 mpp_pe,               &
                                 mpp_root_pe,          &
                                 lowercase,   &
                                 check_nml_error, error_mesg, &
                                 FATAL,WARNING, NOTE,  &
                                 stdlog
use          fms2_io_mod, only : file_exists,          &
                                 open_file,            &
                                 close_file,           &
                                 FmsNetcdfDomainFile_t,&
                                 variable_exists
use      mpp_domains_mod, only : domain2D
use     time_manager_mod, only : time_type, &
                                 print_date, time_manager_init
use     time_interp_mod,  only : time_interp_init
use     diag_manager_mod, only : send_data,            &
                                 register_static_field,&
                                 diag_manager_init, &
				 get_base_time
use   tracer_manager_mod, only : get_tracer_index,     &
                                 query_method
use    field_manager_mod, only : MODEL_ATMOS,          &
                                 parse
use     interpolator_mod, only : interpolate_type,     &
                                 interpolator_init,    &
                                 interpolator_end,     &
                                 interpolator,         &
				 obtain_interpolator_time_slices, &
				 unset_interpolator_time_flag, &
                                 CONSTANT,             &
                                 INTERP_WEIGHTED_P  



implicit none
private
!-----------------------------------------------------------------------
!----- interfaces -------

public  atmos_ozone_tracer, atmos_ozone_tracer_init,  &
	atmos_ozone_tracer_time_vary, atmos_ozone_tracer_endts, &
        atmos_ozone_tracer_end


!-------namelist --------------------------
character(len=64) :: o3_relax_file = 'ko3.COPCAT.nc' ! file for chemical relaxation coefficients
character(len=32) :: o3_relax_name = 'ko3'
character(len=64) :: o3_prod_file = 'O3_Prod.nc' ! file for net chemical production rate
character(len=32) :: o3_prod_name = 'prod' 
namelist /ozone_tracer_nml/ o3_relax_file, o3_prod_file, &
                            o3_relax_name, o3_prod_name

!--- Arrays to help calculate tracer sources/sinks ---

character(len=6), parameter :: module_name = 'tracer'

!--- identification numbers for  diagnostic fields and axes ----


logical :: module_is_initialized=.FALSE.

!type(interpolate_type), save :: o3conc  ! used to read in the ozone climatology
type(interpolate_type), save :: k_relax_in ! used to read in the ozone relaxation coefficient
type(interpolate_type), save :: o3p_in ! used to read in the climatology of ozone net production

!character(len=32) :: specname='ozone'
!character(len=32) :: o3units=''
!---- version number -----
character(len=128) :: version = '$Id: atmos_ozone_tracer.F90, 2018/05/17 $'
character(len=128) :: tagname = '$Name: $'
!-----------------------------------------------------------------------


contains


!#######################################################################
!<SUBROUTINE NAME="atmos_ozone_tracer">
!<OVERVIEW>
! The routine that calculate the sources and sinks of ozone tracer.
!</OVERVIEW>
!<DESCRIPTION>
!</DESCRIPTION>
!<TEMPLATE>
!call atmos_ozone_tracer (lon, lat, pwt, ozone, ozone_dt, Time, kbot)
!</TEMPLATE>
!   <IN NAME="lon" TYPE="real" DIM="(:,:)">
!     Longitude of the centre of the model gridcells
!   </IN>
!   <IN NAME="lat" TYPE="real" DIM="(:,:)">
!     Latitude of the centre of the model gridcells
!   </IN>
!   <IN NAME="pwt" TYPE="real" DIM="(:,:,:)">
!     The pressure weighting array. = dP/grav
!   </IN>
!   <IN NAME="ozone" TYPE="real" DIM="(:,:,:)">
!     The array of the ozone tracer
!   <IN NAME="ozone_clim" TYPE="real" DIM="(:,:,:)">
!     The array of the ozone tracer climatology
!   </IN>
!   <IN NAME="Time" TYPE="type(time_type)">
!     Model time.
!   </IN>
!   <IN NAME="dt" TYPE="real">
!     Model physics timestep (s)
!   </IN>
!   <IN NAME="kbot" TYPE="integer, optional" DIM="(:,:)">
!     Integer array describing which model layer intercepts the surface.
!   </IN>
!   <IN NAME="is, js" TYPE="integer">
!     Local domain start indices
!   </IN>

!   <OUT NAME="ozone_dt" TYPE="real" DIM="(:,:,:)">
!     The array of the tendency of the ozone tracer
!   </OUT>
 subroutine atmos_ozone_tracer (lon, lat, phalf, ozone, ozone_dt, Time, dt, is, js, kbot)

!-----------------------------------------------------------------------
   real, intent(in),  dimension(:,:)   :: lon, lat
   real, intent(in),  dimension(:,:,:) :: ozone
   real, intent(out), dimension(:,:,:) :: ozone_dt
   type(time_type), intent(in)         :: Time     
   real, intent(in)                    :: dt     !timesteps in seconds
   integer, intent(in)        :: is, js
   integer, intent(in),  dimension(:,:), optional :: kbot
   real, intent(in),  dimension(:,:,:) :: phalf

!-----------------------------------------------------------------------

   integer :: i,j,k,id,jd,kd
!   real, dimension(size(ozone,1),size(ozone,2),size(ozone,3)) :: ozone_clim
   real, dimension(size(ozone,1),size(ozone,2),size(ozone,3)) :: k_relax, o3_p, tau_relax

!-----------------------------------------------------------------------

   if (.not. module_is_initialized)  &
      call error_mesg ('atmos_ozone_tracer','atmos_ozone_tracer_init must be called first.', FATAL)

   id=size(ozone,1); jd=size(ozone,2); kd=size(ozone,3)

!!----------- read in ozone climatology------------
!  call interpolator(o3conc, Time, phalf,ozone_clim,trim(specname),is,js)
!  if (trim(o3units)=='kg/kg' .or. trim(o3units)=='mmr')  then
!  ozone_clim = ozone_clim * 29./48.  ! convert from kg/kg to mol/mol
!  endif

!------------read in relaxation coefficient------------------------------
  call interpolator(k_relax_in, Time, phalf, k_relax, trim(o3_relax_name),is,js)  ! in units of 1/s
  tau_relax=-1.0/k_relax

!------------read climatolgy of net ozone production rate----------------
  call interpolator(o3p_in, Time, phalf, o3_p , trim(o3_prod_name), is, js)  ! in units of VMR/s

!  calculate the tendency 
!  o3_p and k_relax are assumed to have no longitudinal variations
   do i = 1,id
   do k = 1,kd
   do j = 1,jd
      ozone_dt(i,j,k)= o3_p(i,j,k)/(1.0-k_relax(i,j,k)*dt)-(1.0/(tau_relax(i,j,k)+dt)) *ozone(i,j,k)
   end do
   end do
   end do

 end subroutine atmos_ozone_tracer
!</SUBROUTINE>

!#######################################################################

!<SUBROUTINE NAME="atmos_ozone_tracer_init">
!<OVERVIEW>
! The constructor routine for the ozone tracer module.
!</OVERVIEW>
!<DESCRIPTION>
! A routine to initialize the ozone tracer module. Also read in data for dry deposition
!</DESCRIPTION>
!<TEMPLATE>
!call atmos_ozone_tracer_init (r, axes, Time, nozone, lonb_mod, latb_mod, phalf, mask)
!</TEMPLATE>
!   <INOUT NAME="r" TYPE="real" DIM="(:,:,:,:)">
!     Tracer fields dimensioned as (nlon,nlat,nlev,ntrace). 
!   </INOUT>
!   <IN NAME="mask" TYPE="real, optional" DIM="(:,:,:)">
!      optional mask (0. or 1.) that designates which grid points
!           are above (=1.) or below (=0.) the ground dimensioned as
!           (nlon,nlat,nlev).
!   </IN>
!   <IN NAME="Time" TYPE="type(time_type)">
!     Model time.
!   </IN>
!   <IN NAME="axes" TYPE="integer" DIM="(4)">
!     The axes relating to the tracer array dimensioned as
!      (nlon, nlat, nlev, ntime)
!   </IN>
!   <IN NAME="lonb_mod" TYPE="real" DIM="(:,:)">
!     The longitude corners for the local domain.
!   </IN>
!   <IN NAME="latb_mod" TYPE="real" DIM="(:,:)">
!     The latitude corners for the local domain.
!   </IN>
!   <IN NAME="phalf" TYPE="real" DIM="(:,:,:)">
!     Pressure on the model half levels (Pa)
!   </IN>
!   <OUT NAME="o3ddep_data" TYPE="interpolate_type" DIM="(1)">
!     file contains ozone dry deposition velocity
!   </OUT>
subroutine atmos_ozone_tracer_init( domain, r, axes, Time, nozone, &
                                   lonb_mod, latb_mod, phalf, mask,o3ddep_data)

!-----------------------------------------------------------------------
!
!   r    = tracer fields dimensioned as (nlon,nlat,nlev,ntrace)
!   mask = optional mask (0. or 1.) that designates which grid points
!          are above (=1.) or below (=0.) the ground dimensioned as
!          (nlon,nlat,nlev).
!
!-----------------------------------------------------------------------
type(domain2D),target,intent(in) :: domain
real,             intent(inout), dimension(:,:,:,:) :: r
type(time_type),  intent(in)                        :: Time
integer,          intent(in)                        :: axes(4)
real,             intent(in),    dimension(:,:)     :: lonb_mod
real,             intent(in),    dimension(:,:)     :: latb_mod
real,             intent(in),    dimension(:,:,:)   :: phalf
integer,          intent(out)                       :: nozone
real, intent(in), dimension(:,:,:), optional        :: mask
type(interpolate_type), intent(out) :: o3ddep_data
!
!-----------------------------------------------------------------------
!
      integer :: n, unit, ierr, io,logunit
      integer :: flag_file, flag_spec, flag_file2
      character(len=64) :: control='', name='', filename='', specname=''
      character(len=64) :: control2='',name2='',filename2=''
      logical :: tracer_initialized
      real, dimension(size(r,1),size(r,2),size(r,3)) :: ozonex
      type(time_type) :: Model_init_time
      character(len=64) :: file_o3dry = 'depvel.nc'  ! default dry deposition velocity file
      type(interpolate_type) :: o3conc  ! used to read in the ozone initial condition
      character(len=32) :: o3units=''

      if (module_is_initialized) return

!---- write namelist ------------------

      call write_version_number (version, tagname)
      if ( mpp_pe() == mpp_root_pe() ) &
       write ( stdlog(), nml=ozone_tracer_nml )

!----------------------------------------------------------------------
!    verify that modules used by this module that are not called later
!    have already been initialized.
!-------------------------------------------------------------------
      call time_interp_init

!-----------------------------------------------------------------------
!    read namelist.
!-----------------------------------------------------------------------
      if ( file_exists('input.nml')) then
        read (input_nml_file,nml=ozone_tracer_nml, iostat=io)
        ierr = check_nml_error(io, 'ozone_tracer_nml')
      endif

!---------------------------------------------------------------------
!    write version number and namelist to logfile.
!---------------------------------------------------------------------
      call write_version_number (version, tagname)
      logunit = stdlog()
      if (mpp_pe() == mpp_root_pe() ) &
                        write (logunit, nml=ozone_tracer_nml)

!----------set initial value-----------------------------------------
      nozone = -1
      n = get_tracer_index(MODEL_ATMOS,'O3' )
!     if (n>0) nozone = n
      if (n>0) then
        nozone = n
        if (mpp_pe() == mpp_root_pe()) then
           write(*,*) 'Ozone was initiated as tracer number ',n
           write(logunit,*) 'Ozone was initiated as tracer number ', n
         end if


!----------------------------------------------------------------------
      
      Model_init_time = get_base_time()
      call print_date (Model_init_time , str='Ozone data is mapped to &
                                              &model time:')

!-----------------------------------------------------------------------
!     ... set initial condition for ozone 
!-----------------------------------------------------------------------
      tracer_initialized = .false.
!      if ( field_exist('INPUT/atmos_tracers.res.nc', 'o3') .or. &
!           field_exist('INPUT/fv_tracer.res.nc', 'o3') .or. &
!           field_exist('INPUT/tracer_o3.res', 'o3') ) then
!         tracer_initialized = .true.
!      end if
      tracer_initialized = check_if_tracer_initialized("o3",domain)

      if (.not. tracer_initialized ) then
      if (query_method('init_conc',MODEL_ATMOS,n,name,control) ) then
	 if (trim(name) == 'file') then
	     flag_file = parse(control, 'file',filename)
	     flag_spec = parse(control,'name',specname)

	     if (flag_file>0 ) then
		call interpolator_init( o3conc,trim(filename),lonb_mod,latb_mod,&
                                        data_out_of_bounds=(/CONSTANT/), &
                                        vert_interp=(/INTERP_WEIGHTED_P/) )
                if( flag_spec > 0 ) then
                   specname = lowercase(specname)
                else
                   specname = 'ozone'
                end if
                call interpolator(o3conc, Time, phalf,ozonex,trim(specname), &
                           clim_units=o3units)
                o3units=lowercase(o3units)
                if (trim(o3units)=='kg/kg' .or. trim(o3units)=='mmr' ) then
                    r(:,:,:,n)=ozonex * 29./48.   ! convert from kg/kg to mol/mol
                else
                    r(:,:,:,n)=ozonex
                end if
             end if
         end if
      end if
      end if

!---------------initiate "k_relax_in" and "o3p_in"-------------------
      call interpolator_init (k_relax_in, trim(o3_relax_file),lonb_mod, &
                              latb_mod, data_out_of_bounds=(/CONSTANT/), &
                              vert_interp=(/INTERP_WEIGHTED_P/) )

      call interpolator_init (o3p_in, trim(o3_prod_file),lonb_mod, &
                              latb_mod, data_out_of_bounds=(/CONSTANT/), &
                              vert_interp=(/INTERP_WEIGHTED_P/) )

!------------------initiate the dry deposition velocity file ----------------
      if (query_method('dry_deposition',MODEL_ATMOS,n,name2,control2) ) then
	 if( trim(name2) == 'file' ) then
            flag_file2 = parse(control2, 'file',filename2)
            if(flag_file2 <= 0 .or. len(trim(filename2))== 0) then
               filename2=file_o3dry
            end if
            call interpolator_init( o3ddep_data, trim(filename2), lonb_mod, &
                                     latb_mod,data_out_of_bounds=(/CONSTANT/), &
                                       vert_interp=(/INTERP_WEIGHTED_P/))
         end if
      end if

!-----------------print out settings for ozone tracer--------------------------
      if (mpp_pe() == mpp_root_pe() ) then
         write(*,*)'Initial concentration from file: ',trim(filename), &
                           ', with the name of ', trim(specname)
         write(*,*)'Dry deposition velocity from file: ',trim(filename2)
         write(logunit,*)'Initial concentration from file: ',trim(filename), &
                           ', with the name of ', trim(specname)
         write(logunit,*)'Dry deposition velocity from file: ',trim(filename2)
      endif

!----------------------------------------------------------------------------- 
      else !if n<=0
!<ERROR MSG="Ozone tracer not found in field table" STATUS="WARNING">
!   Ozone tracer was not included in the field table
!</ERROR>
         call error_mesg ('atmos_ozone_tracer_init', 'ozone' // ' is not found', WARNING)
      end if !if n>0

      module_is_initialized = .TRUE.

!-----------------------------------------------------------------------

 end subroutine atmos_ozone_tracer_init
!</SUBROUTINE>

!> @brief This function just checks if a tracer initialized in a set of files
!! @return flag indicating if a tracer initialized
function check_if_tracer_initialized(tracername, domain) result (tracer_initialized)
  character(len=*), intent(in), optional :: tracername
  type(domain2D),target,intent(in) :: domain !< Atmosphere domain

  logical :: tracer_initialized

  type(FmsNetcdfDomainFile_t) :: fileobj !< fms2io fileobj for domain decomposed

  tracer_initialized = .false.

  if (open_file(fileobj, 'INPUT/atmos_tracers.res.nc', "read", domain)) then
     tracer_initialized = variable_exists(fileobj, tracername)
     call close_file(fileobj)
  elseif (open_file(fileobj, 'INPUT/fv_tracer.res.nc', "read", domain)) then
     tracer_initialized = variable_exists(fileobj, tracername)
     call close_file(fileobj)
  elseif (open_file(fileobj, 'INPUT/tracer_'//trim(lowercase(tracername))//'.res.nc', "read", domain)) then
     tracer_initialized = variable_exists(fileobj, tracername)
     call close_file(fileobj)
  endif

end function check_if_tracer_initialized

!#################################################################### 
 
subroutine atmos_ozone_tracer_time_vary (model_time)
 
!----------------------------------------------------------------------
!     subroutine ozone_time_vary calculates time-dependent, 
!     space-independent variables needed by this module.
!---------------------------------------------------------------------

type(time_type),    intent(in)   :: model_time
 
!----------------------------------------------------------------------

     call obtain_interpolator_time_slices (k_relax_in, model_time)
     call obtain_interpolator_time_slices (o3p_in, model_time)
 
end subroutine atmos_ozone_tracer_time_vary

!</SUBROUTINE>

!#######################################################################

 subroutine atmos_ozone_tracer_endts

      call unset_interpolator_time_flag(k_relax_in)
      call unset_interpolator_time_flag(o3p_in)

 end subroutine atmos_ozone_tracer_endts


!</SUBROUTINE>

!#######################################################################

!<SUBROUTINE NAME="atmos_ozone_tracer_end">
!<OVERVIEW>
!  The destructor routine for the ozone tracer module.
!</OVERVIEW>
! <DESCRIPTION>
! This subroutine writes the version name to logfile and exits.
! </DESCRIPTION>
!<TEMPLATE>
! call atmos_ozone_tracer_end
!</TEMPLATE>
 subroutine atmos_ozone_tracer_end

       call interpolator_end (k_relax_in)
       call interpolator_end (o3p_in)
       module_is_initialized = .FALSE.
 end subroutine atmos_ozone_tracer_end

end module atmos_ozone_tracer_mod
