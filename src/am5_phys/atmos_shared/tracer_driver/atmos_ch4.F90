module atmos_ch4_mod
! <CONTACT EMAIL="Vaishali.Naik@noaa.gov">
!   Vaishali Naik
! </CONTACT>

! <REVIEWER EMAIL="none@noaa.gov">
!   none yet
! </REVIEWER>

! <HISTORY SRC="http://www.gfdl.noaa.gov/fms-cgi-bin/cvsweb.cgi/FMS/"/>

! <OVERVIEW>
! This module initializes and calls ch4 radiation override if an external 
! field is used for radiation calculations. Not used at the moment
! </OVERVIEW>

! <DESCRIPTION>
! </DESCRIPTION>


use mpp_mod, only: input_nml_file 
use              fms_mod, only : write_version_number,    &
                                 mpp_pe, mpp_root_pe,                 &
                                 stdlog, stdout,          &
                                 check_nml_error, error_mesg,         &
                                 FATAL, NOTE, WARNING
use          fms2_io_mod, only : file_exists
use   tracer_manager_mod, only : get_tracer_index, tracer_manager_init
use    field_manager_mod, only : MODEL_ATMOS
use     time_manager_mod, only : time_type
use    data_override_mod, only : data_override
use              mpp_mod, only : mpp_pe, mpp_root_pe


implicit none

private

!-----------------------------------------------------------------------
!----- interfaces -------

public  atmos_ch4_rad
public  atmos_ch4_rad_init

!-----------------------------------------------------------------------
!----------- module data -------------------
!-----------------------------------------------------------------------

character(len=48), parameter    :: mod_name = 'atmos_ch4_mod'

integer, save   :: ind_ch4  = 0
real            :: radiation_ch4_dvmr = -1


!---------------------------------------------------------------------
!-------- namelist  ---------
!-----------------------------------------------------------------------

logical  :: ch4_radiation_override = .false.

namelist /atmos_ch4_nml/  &
          ch4_radiation_override

!PUBLIC VARIABLES
public :: ch4_radiation_override

logical :: module_is_initialized = .FALSE.
integer :: logunit

!---- version number -----
character(len=128) :: version = '$$'
character(len=128) :: tagname = '$$'

!-----------------------------------------------------------------------

contains

!#######################################################################

!<SUBROUTINE NAME ="atmos_ch4_rad">

!<OVERVIEW>
! Subroutine to get global avg ch4 to be used in radiation.
! input ch4 field is from data override 
!</OVERVIEW>

 subroutine atmos_ch4_rad(Time, radiation_ch4_dvmr)

!
!-----------------------------------------------------------------------
!     arguments
!-----------------------------------------------------------------------
!
   type (time_type),      intent(in)    :: Time
   real,                  intent(inout) :: radiation_ch4_dvmr
!
!-----------------------------------------------------------------------
!     local variables
!-----------------------------------------------------------------------
!

!
!-----------------------------------------------------------------------
!     local parameters
!-----------------------------------------------------------------------
!

character(len=64), parameter    :: sub_name = 'atmos_ch4_rad'
character(len=256), parameter   :: error_header =                               &
     '==>Error from ' // trim(mod_name) // '(' // trim(sub_name) // '):'
character(len=256), parameter   :: warn_header =                                &
     '==>Warning from ' // trim(mod_name) // '(' // trim(sub_name) // '):'
character(len=256), parameter   :: note_header =                                &
     '==>Note from ' // trim(mod_name) // '(' // trim(sub_name) // '):'
logical   :: used
!-----------------------------------------------------------------------

if (ind_ch4 > 0 .and. ch4_radiation_override) then

! input is in moist vmr (mol/mol) that is converted to dry in vmr in radiation_types.F90

  call data_override('ATM', 'ch4_dvmr_rad', radiation_ch4_dvmr, Time, override=used)
  if (.not. used) then
    call error_mesg (trim(error_header), ' data override needed for ch4_dvmr_rad ', FATAL)
  endif
!  if (mpp_pe() == mpp_root_pe() ) &
!      write (logunit,*)' atmos_ch4_rad       : mean radiation ch4_dvmr = ',radiation_ch4_dvmr

!else
!  if (mpp_pe() == mpp_root_pe() ) &
!      write (logunit,*)' atmos_ch4_rad: ch4 radiation override not active: ',ch4_radiation_override
endif


!-----------------------------------------------------------------------

end subroutine atmos_ch4_rad
!</SUBROUTINE>

!#######################################################################

!<SUBROUTINE NAME ="atmos_ch4_rad_init">

!<OVERVIEW>
! Subroutine to initialize the CH4 module for radiation calculations.
!</OVERVIEW>

 subroutine atmos_ch4_rad_init

!
!-----------------------------------------------------------------------
!     arguments
!-----------------------------------------------------------------------
!


!-----------------------------------------------------------------------
!     local variables
!         unit       io unit number used to read namelist file
!         ierr       error code
!         io         error status returned from io operation
!-----------------------------------------------------------------------
!
integer :: ierr, io
integer :: n
integer :: outunit
!
!-----------------------------------------------------------------------
!     local parameters
!-----------------------------------------------------------------------
!

character(len=64), parameter    :: sub_name = 'atmos_ch4_init'
character(len=256), parameter   :: error_header =                               &
     '==>Error from ' // trim(mod_name) // '(' // trim(sub_name) // '):'
character(len=256), parameter   :: warn_header =                                &
     '==>Warning from ' // trim(mod_name) // '(' // trim(sub_name) // '):'
character(len=256), parameter   :: note_header =                                &
     '==>Note from ' // trim(mod_name) // '(' // trim(sub_name) // '):'


     if (module_is_initialized) return

     call write_version_number (version, tagname)

!-----------------------------------------------------------------------
!    read namelist.
!-----------------------------------------------------------------------
      if ( file_exists('input.nml')) then
        read (input_nml_file, nml=atmos_ch4_nml, iostat=io)
        ierr = check_nml_error(io,'atmos_ch4_nml')
      endif

!---------------------------------------------------------------------
!    write namelist to logfile.
!---------------------------------------------------------------------
      logunit=stdlog()
      if (mpp_pe() == mpp_root_pe() ) &
                          write (logunit, nml=atmos_ch4_nml)

!----- set initial value of methane ------------

n = get_tracer_index(MODEL_ATMOS,'ch4')
if (n > 0) then
  ind_ch4 = n
    outunit=stdout()
    write (outunit,*) trim(note_header), ' ch4 was initialized as tracer number ', ind_ch4
endif

if (.not.(ind_ch4 > 0 .and. ch4_radiation_override)) then
   if (mpp_pe() == mpp_root_pe() ) &
     write (logunit,*)' ch4 radiation override not active:ch4_radiation_override = ',ch4_radiation_override
endif


call write_version_number (version, tagname)
module_is_initialized = .TRUE.


!-----------------------------------------------------------------------

end subroutine atmos_ch4_rad_init
!</SUBROUTINE>

end module atmos_ch4_mod
