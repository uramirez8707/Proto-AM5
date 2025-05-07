                    module ls_cloud_macrophysics_mod

!-----------------------------------------------------------------------
!
!         interface module for grid-scale moisture processes when
!                 prognostic clouds are active
!
!                ---------------------------------------
!                      OPTIONS AVAILABLE:
!               a) Tiedtke stratiform prognostic cloud scheme
!
!-----------------------------------------------------------------------

! fms modules
use time_manager_mod,      only: time_type
use mpp_mod,               only: input_nml_file
use fms_mod,               only: error_mesg, FATAL, NOTE,        &
                                 check_nml_error,    &
                                 write_version_number,           &
                                 mpp_pe, mpp_root_pe, stdlog,    &
                                 mpp_clock_id, mpp_clock_begin,  &
                                 mpp_clock_end, CLOCK_MODULE,    &
                                 CLOCK_MODULE_DRIVER, &
                                 MPP_CLOCK_SYNC
use fms2_io_mod,           only: file_exists
use physics_types_mod,     only: physics_control_type

! atmos_param modules
use tiedtke_macro_mod,     only: tiedtke_macro, tiedtke_macro_init, &
                                 tiedtke_macro_end,  &
                                 tiedtke_macro_time_vary, &
                                 tiedtke_macro_diagnostics
use aerosol_cloud_mod,     only: determine_activated_aerosol
use lscloud_types_mod,     only: atmos_state_type, &
                                 particles_type, cloud_state_type, &
                                 lsc_constants_type, lscloud_nml_type, &
                                 cloud_processes_type, precip_state_type

use aerosol_types_mod,     only: aerosol_type
use moist_proc_utils_mod,  only: mp_input_type, mp_output_type, &
                                 mp_lsdiag_type, mp_nml_type,   &
                                 mp_lsdiag_control_type, &
                                 mp_conv2ls_type, mp_tendency_type
use physics_radiation_exch_mod,            &
                          only : exchange_control_type

implicit none
private

!-----------------------------------------------------------------------
!-------------------- public data/interfaces ---------------------------

public   ls_cloud_macrophysics_init, ls_cloud_macrophysics, &
         ls_cloud_macrophysics_end, ls_cloud_macrophysics_time_vary


!-----------------------------------------------------------------------
!-------------------- private data -------------------------------------

!--------------------- version number ----------------------------------
character(len=128) :: version = '$Id: $'
character(len=128) :: tagname = '$Name: $'


!-------------------- namelist data (private) --------------------------


!------ variables obtained during initialization from other nmls -----

real       :: qmin
logical    :: do_liq_num, do_ice_num
logical    :: do_pdf_clouds
logical    :: tiedtke_macrophysics
integer    :: nsphum, nql, nqi, nqa, nqn, nqni, nqr, nqs, nqg

!-------------------- clock definitions --------------------------------

integer   ::  tiedtke_clock
integer   :: active_clock, main_clock, diag_clock

logical    :: module_is_initialized = .false.




                             contains


!#######################################################################

subroutine ls_cloud_macrophysics_init     &
                    ( id, jd, kd, lon, lat, axes, Time, phalf, Nml_mp, &
                     Constants_lsc, Physics_control, Nml_lsc, Exch_ctrl)

!-------------------------------------------------------------------------
!   initialize the active large-scale cloud scheme. currently 
!   tiedtke is available.
!-------------------------------------------------------------------------

integer,                     intent(in)     :: id, jd, kd
integer,                     intent(in)     :: axes(4)
type(time_type),             intent(in)     :: Time
real,dimension(:,:),         intent(in)     :: lon,  lat    ! h1g
real,dimension(:,:,:),       intent(in)     :: phalf        ! h1g
type(lsc_constants_type),    intent (inout) :: Constants_lsc
type(physics_control_type),  intent(in)     :: Physics_control
type(lscloud_nml_type),      intent (in)    :: Nml_lsc
type(mp_nml_type),           intent (in)    :: Nml_mp
type(exchange_control_type), intent(in)     :: Exch_ctrl

!------------------------------------------------------------------------
!  local variables:

      integer :: ierr, io, logunit
      integer :: tiedtke_init_clock

!-----------------------------------------------------------------------

      if (module_is_initialized) return

!-------------------------------------------------------------------------
!    save control variables obtained from derived type inputs which will
!    be needed in this module.
!-------------------------------------------------------------------------
      do_liq_num = Exch_ctrl%do_liq_num
      do_ice_num = Exch_ctrl%do_ice_num
      do_pdf_clouds = Nml_lsc%do_pdf_clouds
      tiedtke_macrophysics = Constants_lsc%tiedtke_macrophysics
      qmin = Exch_ctrl%qmin

      nsphum = Physics_control%nsphum
      nql    = Physics_control%nql
      nqi    = Physics_control%nqi
      nqa    = Physics_control%nqa
      nqn    = Physics_control%nqn
      nqni   = Physics_control%nqni
      nqr    = Physics_control%nqr
      nqs    = Physics_control%nqs
      nqg    = Physics_control%nqg

!------------------------------------------------------------------------
!    if doing prognostic clouds, then tiedtke microphysics must be activated.
!------------------------------------------------------------------------
      if (.not. tiedtke_macrophysics ) then
        call error_mesg ('ls_cloud_macrophysics_mod', &
              'with prognostic clouds tiedtke must be active &
              to include large-scale condensation', FATAL)
      endif

!------------------------------------------------------------------------
!    define clocks for large-scale cloud schemes initialization (local
!    variables) and for the prognostic loop clocks (module variables).
!------------------------------------------------------------------------
      tiedtke_init_clock = mpp_clock_id(     &
                     '   Lscloud_driver: tiedtke:Initialization' , &
                                                grain=CLOCK_MODULE_DRIVER )
      tiedtke_clock = mpp_clock_id(     &
                     '   Lscloud_driver: tiedtke' , &
                                                grain=CLOCK_MODULE_DRIVER )
      active_clock = mpp_clock_id(    &
                     '   Ls_macrophysics: active' , &
                                                grain=CLOCK_MODULE_DRIVER )
      main_clock = mpp_clock_id(    &
                     '   Ls_macrophysics: main' , &
                                                grain=CLOCK_MODULE_DRIVER )
      diag_clock = mpp_clock_id(    &
                     '   Ls_macrophysics: diag' , &
                                                grain=CLOCK_MODULE_DRIVER )

!------------------------------------------------------------------------
!    if tiedtke macrophysics is active, initialize it here.
!------------------------------------------------------------------------
      if (tiedtke_macrophysics) then
        call mpp_clock_begin (tiedtke_init_clock)
        call tiedtke_macro_init (Nml_mp, Nml_lsc, Constants_lsc, Exch_ctrl)
        call mpp_clock_end   (tiedtke_init_clock)
      endif

!----------------------------------------------------------------------

end subroutine ls_cloud_macrophysics_init



!#######################################################################

subroutine ls_cloud_macrophysics_time_vary (dtcloud_in)

real, intent(in) :: dtcloud_in

!-------------------------------------------------------------------------

      if (tiedtke_macrophysics) then
        call tiedtke_macro_time_vary (dtcloud_in)
      endif

!-------------------------------------------------------------------------

end subroutine ls_cloud_macrophysics_time_vary

!#######################################################################

subroutine ls_cloud_macrophysics (   &
             is, ie, js, je, Time, dt, rdiag, Input_mp, Output_mp, Tend_mp,&
             C2ls_mp, Lsdiag_mp_control, Lsdiag_mp, Atmos_state,   &
             Cloud_state, Particles, Precip_state, Cloud_processes, Aerosol)

!-------------------------------------------------------------------------

integer,                    intent(in)           :: is,ie,js,je
type(time_type),            intent(in)           :: Time
real,                       intent(in)           :: dt
type(mp_input_type),        intent(inout)        :: Input_mp
type(mp_output_type),       intent(inout)        :: Output_mp
type(mp_tendency_type),     intent(inout)        :: Tend_mp
type(mp_conv2ls_type),      intent(inout)        :: C2ls_mp
type(mp_lsdiag_type),       intent(inout)        :: Lsdiag_mp
type(mp_lsdiag_control_type),intent(inout)       :: Lsdiag_mp_control
type(atmos_state_type),     intent(inout)        :: Atmos_state
type(cloud_state_type),     intent(inout)        :: Cloud_state
type(particles_type),       intent(inout)        :: Particles
type(precip_state_type),    intent(inout)        :: Precip_state
type(cloud_processes_type), intent(inout)        :: Cloud_processes
type(aerosol_type),         intent(in)           :: Aerosol
real, dimension(:,:,:,size(Output_mp%rdt,4)+1:),     &
                            intent(inout)        ::  rdiag


!-----------------------------------------------------------------------
!   local variables:

      integer :: idim, jdim, kdim

!--------------------------------------------------------------------------
!    define window dimensions.
!--------------------------------------------------------------------------
      idim = size(Input_mp%t,1)
      jdim = size(Input_mp%t,2)
      kdim = size(Input_mp%t,3)

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!       A. TIEDTKE PROGNOSTIC CLOUD SCHEME
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

      if (tiedtke_macrophysics) then

        call mpp_clock_begin (tiedtke_clock)
!--------------------------------------------------------------------------
!   if the tiedtke scheme is active, determine aerosol available for use as
!   condensation nuclei, then call tiedtke_macro to calculate the changes
!   in the large-scale cloud amount and area, and then call
!   tiedtke_macro_diagnostics to save relevant diagnostics.
!--------------------------------------------------------------------------
        call mpp_clock_begin (active_clock)
        call determine_activated_aerosol (   &
             idim, jdim, kdim, Lsdiag_mp_control%n_diag_4d, C2ls_mp,   &
             Input_mp, Atmos_state, Particles, Cloud_state%qa_upd,   &
             Lsdiag_mp%diag_4d, Lsdiag_mp_control%diag_id,    &
                                                 Lsdiag_mp_control%diag_pt )
        call mpp_clock_end (active_clock)

        call mpp_clock_begin (main_clock)
        CALL tiedtke_macro (   &
             idim, jdim, kdim, C2ls_mp, Input_mp, Atmos_state,   &
             Cloud_state, Tend_mp%ttnd, Tend_mp%qtnd, Cloud_processes,  &
                                  Particles, Lsdiag_mp, Lsdiag_mp_control )
        call mpp_clock_end (main_clock)

        call mpp_clock_begin (diag_clock)
        call tiedtke_macro_diagnostics (    &
             idim, jdim, kdim, Lsdiag_mp_control%n_diag_4d,   &
             Lsdiag_mp_control%diag_id,  &
             Lsdiag_mp%diag_4d, Lsdiag_mp_control%diag_pt,    &
             Cloud_processes, Cloud_state)
        call mpp_clock_end (diag_clock)

        call mpp_clock_end (tiedtke_clock)
      endif

!------------------------------------------------------------------------


end subroutine ls_cloud_macrophysics


!######################################################################

subroutine ls_cloud_macrophysics_end

!------------------------------------------------------------------------
      integer :: tiedtke_term_clock

!------------------------------------------------------------------------
!   define clocks for termination.
!------------------------------------------------------------------------
      tiedtke_term_clock = mpp_clock_id(     &
                    '   Lscloud_driver: tiedtke:Termination' , &
                                                grain=CLOCK_MODULE_DRIVER )
      if (tiedtke_macrophysics) then
        call mpp_clock_begin ( tiedtke_term_clock )
        call tiedtke_macro_end
        call mpp_clock_end   ( tiedtke_term_clock )
      endif

!-----------------------------------------------------------------------

end subroutine ls_cloud_macrophysics_end


!######################################################################




                 end module ls_cloud_macrophysics_mod
