
                    module convection_driver_mod

!-----------------------------------------------------------------------
!
!         This module accesses the following convective parameterizations.
!         Access is controlled by logical variables found in 
!         moist_processes_nml.
!
!         ---------------------------------------
!         1)  uw convection
!
!         The following convective implementations are available:
!         -------------------------------------------------
!         1)  uw convection  alone
!             
!------------------------------------------------------------------------
!
!       The module also calculates cumulus momentum transport (if desired,
!    using either a diffusive or a non-local scheme), and produces needed
!    fields for input to the large-scale cloud scheme and COSP, updates 
!    the time tendencies of model variables as a result of convection, and
!    outputs convective diagnostics.
!    
!-----------------------------------------------------------------------



! infrastructure modules
use sat_vapor_pres_mod,     only: compute_qs
use time_manager_mod,       only: time_type, get_time, assignment(=)
use diag_manager_mod,       only: register_diag_field, send_data, &
                                  get_diag_field_id, DIAG_FIELD_NOT_FOUND
use diag_data_mod,          only: CMOR_MISSING_VALUE
use mpp_domains_mod,        only: domain2D
use mpp_mod,                only: input_nml_file
use fms2_io_mod,            only: file_exists
use fms_mod,                only: error_mesg, FATAL, WARNING,NOTE,&
                                  check_nml_error,    &
                                  write_version_number,           &
                                  mpp_pe, mpp_root_pe, stdlog,    &
                                  mpp_clock_id, mpp_clock_begin,  &
                                  mpp_clock_end, CLOCK_MODULE,    &
                                  CLOCK_MODULE_DRIVER, &
                                  MPP_CLOCK_SYNC
use field_manager_mod,      only: MODEL_ATMOS
use tracer_manager_mod,     only: get_tracer_index,&
                                  get_tracer_names, &
                                  NO_TRACER
use constants_mod,          only: CP_AIR, HLV, HLS, HLF, RDGAS, RVGAS, &
                                  SECONDS_PER_DAY, KAPPA, AVOGNO
use atmos_global_diag_mod,  only: register_global_diag_field, &
                                  buffer_global_diag
use atmos_cmip_diag_mod,    only: register_cmip_diag_field_2d, &
                                  register_cmip_diag_field_3d, &
                                  send_cmip_data_3d, &
                                  query_cmip_diag_id, &
                                  cmip_diag_id_type

! atmos modules
use physics_types_mod,      only: physics_control_type, phys_mp_exch_type
use vert_diff_driver_mod,   only: surf_diff_type
use physics_radiation_exch_mod,        &
                            only: clouds_from_moist_block_type, &
                                  exchange_control_type,  &
                                  cloud_scheme_data_type
use uw_conv_mod,            only: uw_conv, uw_conv_end, uw_conv_init
use detr_ice_num_mod,       only: detr_ice_num, detr_ice_num_init,   &
                                  detr_ice_num_end
use cu_mo_trans_mod,        only: cu_mo_trans_init, cu_mo_trans,   &
                                  cu_mo_trans_end
use moz_hook_mod,           only: moz_hook
use aerosol_types_mod,      only: aerosol_type
use moist_proc_utils_mod,   only: capecalcnew, column_diag, rh_calc, &
                                  mp_nml_type, mp_input_type,   &
                                  mp_tendency_type, mp_removal_type, &
                                  mp_removal_control_type, &
                                  mp_conv2ls_type, mp_output_type
use convection_utilities_mod,         &
                            only: conv_tendency_type, &
                                  conv_output_type,   &
                                  conv_results_type
use atmos_tracer_utilities_mod,         &
                            only: wet_deposition

implicit none
private

!-----------------------------------------------------------------------
!-------------------- public data/interfaces ---------------------------

public     convection_driver_init, convection_driver_time_vary,   &
           convection_driver, convection_driver_endts,  & 
           convection_driver_end,    &
           convection_driver_restart, cape_cin_diagnostics
  
private                             &
!   initialization:
           diag_field_init,    &
!   called in prognostic loop:
           convection_driver_alloc, define_total_convective_output,     &
           convective_diagnostics,  define_convective_area, &
           compute_convective_area,  & 
           define_inputs_for_cosp, prevent_neg_precip_fluxes, &
           convection_driver_dealloc, &
!   associated with uw convection:
           uw_conv_driver, uw_conv_driver_part,  uw_alloc,  &
           finalize_uw_outputs, update_inputs, uw_diagnostics, &
           uw_dealloc, &
!   routines accessed by multiple convection implementations:
           update_outputs, define_and_apply_scale

!-----------------------------------------------------------------------
!-------------------- private data -------------------------------------

!--------------------- version number ----------------------------------
character(len=128) ::  version = '$Id: $'
character(len=128) :: tagname = '$Name: $'


!-------------------- namelist data (private) --------------------------

!---------------- namelist variable definitions -----------------------

!
!   do_cmt     = compute cumulus momentum transport (default = T).
!   cmt_mass_flux_source = parameterization(s) being used to supply the 
!                mass flux profiles seen by the cumulus momentum transport
!                module; currently 'uw', 'all' (default = 'all').
 
!   do_gust_cv = switch to use convective gustiness (default = F).
!   do_gust_cv_new = switch to use newer convective gustiness expression   
!                                                    (default = F).
!   gustmax    = maximum convective gustiness (m/s); default = 3.
!   gustconst  = precip rate which defines precip rate which begins to
!                matter for convective gustiness (kg/m2/sec)
!                default = 1. cm/day = 10. mm/da

!   do_limit_uw = limit UW shallow tendencies to prevent the formation
!                of grid points with negative total water specific 
!                humidities. This situation can occur because both
!                shallow and deep convection operate on the same
!                soundings without knowledge of what the other is doing
!                (default = T).
!   detrain_ice_num = if true, convective ice particles may be detrained
!                into the large-scale clouds (default = F).
!   reproduce_AM4 = setting this to T will reproduce legacy 
!                (warsaw) results because the model will use temperature 
!                and tracer fields that have not been updated with the uw 
!                convective tendency for the calculation of cu_mo_trans 
!                and lightning (updated values should have been used for 
!                consistency). ANSWERS THUS WILL CHANGE from warsaw results
!                for any runs using both uw convection and cu_mo_trans (or
!                tracer "no") when this variable is set to F, since in that
!                case all fields are updated with uw tendencies before the 
!                cumulus momentum transport and lightning calculations 
!                are done. (default = F).
!----------------------------------------------------------------------

logical            :: do_cmt=.true.
character(len=64)  :: cmt_mass_flux_source = 'uw'
logical            :: do_gust_cv = .false.
logical            :: do_gust_cv_new = .false.
real               :: gustmax = 3.     
real               :: gustconst = 10./SECONDS_PER_DAY 
logical            :: do_limit_uw = .false.
logical            :: detrain_ice_num = .false.
logical            :: reproduce_AM4 = .true.


namelist /convection_driver_nml/    &
              do_cmt, cmt_mass_flux_source,   &
              do_gust_cv,  do_gust_cv_new, gustmax, gustconst, &
              do_limit_uw, &
              detrain_ice_num, &
              reproduce_AM4


!------------------- other global variables  ------------------------


!  variables imported during initialization and saved as module variables:
real    :: qmin                    ! minimum value for condensate specific
                                   ! humidity
logical :: do_liq_num              ! prognostic cloud droplet number ?
logical :: do_ice_num              ! prognostic ice particle number ?
logical :: do_cosp                 ! call COSP diagnostic package ?
logical :: do_lsc                  ! using bulk large scale condensation ?
logical :: do_uw_conv              ! uw convection scheme is active ?
integer :: num_uw_tracers          ! number of tracers to be transported 
                                   ! by the uw scheme
logical :: doing_prog_clouds       ! prognostic clouds are active ?
integer :: nsphum,   &             ! tracer index for specific humidity
           nql,      &             ! tracer index for cloud water
           nqi,      &             ! tracer index for cloud ice
           nqa,      &             ! tracer index for cloud area
           nqn,      &             ! tracer index for cloud droplet number
           nqni,     &             ! tracer index for cloud ice particle
                                   !                                number
           nqr,      &             ! tracer index for rain water
           nqs,      &             ! tracer index for snow
           nqg                     ! tracer index for graupel
integer :: num_prog_tracers        ! total number of prognostic tracers
logical, dimension(:), allocatable ::   &
           cloud_tracer            ! logical array indicating which tracers
                                   ! are cloud tracers 
logical, dimension(:), allocatable ::   &
           tracers_in_uw           ! logical array indicating which tracers
                                   ! are transported by uw convection

!  variables imported during  the time_vary code segment and saved 
!  as module variables:
real    :: dt                      ! model timestep [s]
real    :: dtinv                   ! inverse of model timestep
type(time_type) :: Time            ! Time at end of current step, used for
                                   ! diagnostics
integer :: i_shallow               ! index of uw clouds in cloud array

!  variables used to define the active convective implementation: 
logical :: luwconv = .false.       ! uw only


!-------------------- clock definitions --------------------------------

integer :: convection_clock,  &    ! clock to time total convection 
           uw_clock,     &         ! clock to time uw parameterization
           cmt_clock               ! clock to time cumulus momentum
                                   !                  transport calculation

!-------------------- diagnostics variables ----------------------------

!CMIP diagnostics:
integer, public          :: id_pr_g, id_prc_g, id_prsn_g, id_prsnc, id_prrc
integer                  :: id_prc, id_ci, id_ccb, id_cct
type(cmip_diag_id_type)  :: ID_tntc, ID_tnhusc, ID_mc, ID_emilnox_area

! uw diagnostics:
integer :: id_tdt_uw, id_qdt_uw, &
           id_qadt_uw, id_qldt_uw, id_qidt_uw, id_qndt_uw, id_qnidt_uw, &
           id_scale_uw, &
           id_uw_precip, id_uw_snow, id_uw_freq

! total convection diagnostics:
integer :: id_tdt_conv, id_qdt_conv, id_prec_conv, id_snow_conv, &
           id_conv_freq, id_gust_conv, id_mc_full, id_mc_half,   &
           id_mc_conv_up, id_conv_cld_base, id_conv_cld_top, &
           id_qadt_conv, id_qldt_conv, id_qndt_conv, id_qidt_conv, &
           id_qnidt_conv, &
           id_qa_conv_col, id_ql_conv_col, id_qn_conv_col,  &
           id_qni_conv_col, id_qi_conv_col, &
           id_q_conv_col, id_t_conv_col,              &
           id_enth_conv_col, id_wat_conv_col, &
           id_enth_uw_col, id_wat_uw_col, &
           id_prod_no, id_conv_rain3d, id_conv_snow3d
integer, dimension(:), allocatable ::    &
                                      id_tracerdt_conv,  &
                                      id_tracerdt_conv_col, &
                                      id_conv_tracer,  &
                                      id_conv_tracer_col

! cape-cin diagnostics:
integer :: id_cape, id_cin, id_tp, id_rp, id_lcl, id_lfc, id_lzb


character(len=5) :: mod_name = 'moist'
character(len=8) :: mod_name2 = 'moist_tr'
real             :: missing_value = -999.


logical :: module_is_initialized = .false.



                          contains


!*******************************************************************
!
!                     PUBLIC SUBROUTINES
!
!*******************************************************************



!#######################################################################

subroutine convection_driver_init       &
             (domain, id, jd, kd, axes, Time, Physics_control, Exch_ctrl,    &
                                      Nml_mp, Control, lonb, latb, pref )
 
!---------------------------------------------------------------------
!    subroutine convection_driver_init processes convection_driver_nml,
!    saves input variables that are needed as module variables, checks for
!    consistency between specified integration choices, sets up logical
!    variables defining the integration path for convection, initializes
!    the active convective parameterization(s), and initializes the 
!    convective clocks and convective diagnostics.
!---------------------------------------------------------------------

!---------------------------------------------------------------------
type(domain2D), target,        intent(in)    :: domain !< Atmosphere domain
integer,                       intent(in)    :: id, jd, kd
integer,                       intent(in)    :: axes(4)
type(time_type),               intent(in)    :: Time
type(physics_control_type),    intent(in)    :: Physics_control
type(exchange_control_type),   intent(in)    :: Exch_ctrl
type(mp_nml_type),             intent(in)    :: Nml_mp
type(mp_removal_control_type), intent(in)    :: Control   
real, dimension(:,:),          intent(in)    :: lonb, latb
real, dimension(:),            intent(in)    :: pref
!---------------------------------------------------------------------
 
!------------------------------------------------------------------------
!      id, jd, kd  dimensions of processor window
!      axes        data axis indices, (x,y,pf,ph) for diagnostics 
!      Time        time used for diagnostics [time_type]
!      Physics_control 
!                  derived type containing control variables needed by
!                  multiple physics modules
!      Exch_ctrl   derived type variable containing control variables
!                  needed by both physics and radiation modules
!      Nml_mp      derived type variable containing the moist_processes_nml
!                  variables
!      Control     derived type variable containing control variables
!                  associated with tracer removal and transport by
!                  available convective schemes
!      lonb        longitude of grid box corners [ radians ]
!      latb        latitude of grid box corners [ radians ]
!      pref        array of reference pressures at full levels (plus 
!                  surface value at nlev+1), based on 1013.25 hPa pstar
!                  [ Pa ]
!------------------------------------------------------------------------

      integer :: secs, days  
      integer :: io, ierr, logunit

!------------------------------------------------------------------------
!      secs    seconds component of time_type variable Time
!      days    days component of time_type variable Time
!      unit    unit number used to read nml file
!      io      error return code
!      ierr    error return flag
!      logunit unit number used for stdlog file
!------------------------------------------------------------------------

!-----------------------------------------------------------------------
      if ( module_is_initialized ) return

!-----------------------------------------------------------------------
!    process the convection_driver_nml.
!-----------------------------------------------------------------------
      if ( file_exists('input.nml')) then
        read (input_nml_file, nml=convection_driver_nml, iostat=io)
        ierr = check_nml_error(io,'convection_driver_nml')
!----------------------------------------------------------------------
!    write version and namelist to standard logfile.
!----------------------------------------------------------------------
        call write_version_number ( version, tagname )
        logunit = stdlog()
        if ( mpp_pe() == mpp_root_pe() ) &
                        write ( logunit, nml=convection_driver_nml )
      endif

!----------------------------------------------------------------------
!    define needed module variables supplied by the derived types input 
!    to this subroutine.
!----------------------------------------------------------------------
      qmin = Exch_ctrl%qmin
      do_liq_num = Exch_ctrl%do_liq_num
      do_ice_num = Exch_ctrl%do_ice_num
      do_lsc = Nml_mp%do_lsc
      do_cosp = Exch_ctrl%do_cosp
      do_uw_conv = Nml_mp%do_uw_conv
      doing_prog_clouds = Exch_ctrl%doing_prog_clouds
      nsphum = Physics_control%nsphum
      nql = Physics_control%nql
      nqi = Physics_control%nqi
      nqa = Physics_control%nqa
      nqn = Physics_control%nqn
      nqni = Physics_control%nqni
      nqr = Physics_control%nqr
      nqs = Physics_control%nqs
      nqg = Physics_control%nqg
      num_prog_tracers = Physics_control%num_prog_tracers

      num_uw_tracers = Control%num_uw_tracers
 
      allocate (tracers_in_uw (num_prog_tracers))

      tracers_in_uw = Control%tracers_in_uw    

      allocate (cloud_tracer(size(Physics_control%cloud_tracer)))
      cloud_tracer = Physics_control%cloud_tracer

!----------------------------------------------------------------------
!    define logical controls indicating the status of the available
!    convective implementations in this experiment.
!----------------------------------------------------------------------
      if (do_uw_conv) luwconv = .true.

      if (do_cmt) then
        if ( .not.  do_uw_conv) then
          call error_mesg ( 'convection_driver_init', &
                'do_cmt is active but no cumulus schemes activated', &
                                                              FATAL)
        endif
      endif

!-------------------------------------------------------------------------
!    initialize the cumulus momentum transport module, defining logicals 
!    indicating which convective schemes are to be seen by that module.
!-------------------------------------------------------------------------
      if (do_cmt) then
        call cu_mo_trans_init (axes, Time, Nml_mp, cmt_mass_flux_source)
      endif 

!--------------------------------------------------------------------
!    continue the initialization of the convection scheme modules.
!--------------------------------------------------------------------
      if (do_uw_conv) call uw_conv_init (doing_prog_clouds, axes, Time,   &
                                         kd, Nml_mp, Control%tracers_in_uw)

!-----------------------------------------------------------------------
!    initialize clocks.
!-----------------------------------------------------------------------
      convection_clock = mpp_clock_id( '   Physics_up: Moist Proc: Conv' ,&
                                             grain=CLOCK_MODULE_DRIVER )
      uw_clock         = mpp_clock_id( '   Moist Processes: UW'  ,&
                                             grain=CLOCK_MODULE_DRIVER )
      cmt_clock        = mpp_clock_id( '   Moist Processes: CMT'         ,&
                                             grain=CLOCK_MODULE_DRIVER )
!------------------------------------------------------------------------
!    call diag_field_init to register the netcdf diagnostic fields.
!------------------------------------------------------------------------
      call diag_field_init (axes, Time, Control)

!------------------------------------------------------------------------
!    initialize the ice number detrain module.
!------------------------------------------------------------------------
      call detr_ice_num_init

!---------------------------------------------------------------------
      module_is_initialized = .true.


end subroutine convection_driver_init


!#######################################################################

subroutine convection_driver_time_vary(Time_in, dt_in, i_shallow_in)

!------------------------------------------------------------------------
!    subroutine convection_driver_time_vary saves needed input arguments
!    as module variables for use within the prognostic spatial loops 
!    and verifies that needed cloud inputs will be available.
!------------------------------------------------------------------------

!-------------------------------------------------------------------------
type(time_type), intent(in) :: Time_in 
real,            intent(in) :: dt_in  
integer,         intent(in) :: i_shallow_in
!-------------------------------------------------------------------------

!-------------------------------------------------------------------------
!    Time_in       ! time used for diagnostics [ time_type ]
!    dt_in         ! time step [ seconds ]
!    i_shallow_in  ! index in cloud arrays for uw clouds
!-------------------------------------------------------------------------

!---------------------------------------------------------------------
!    verify that the module has been initialized.
!---------------------------------------------------------------------
      if (.not. module_is_initialized) then
        call error_mesg (   &
                 'convection_driver_mod/convection_driver_time_vary',  &
                     'convection_driver_init has not been called.', FATAL)
      endif

!------------------------------------------------------------------------
!    save the input arguments as module variables of this module. Pass them
!    to dependent modules as required.
!------------------------------------------------------------------------
      dt = dt_in
      dtinv = 1./dt

      i_shallow = i_shallow_in

      Time = Time_in
      
!---------------------------------------------------------------------
!    if uw parameterization is active, be sure the uw cloud field argument 
!    is valid.
!---------------------------------------------------------------------
      if (do_uw_conv) then
        if (i_shallow /= 0) then
        else
          call error_mesg ('convection_driver_mod',   &
                  'input args for uw shallow clouds not correct', FATAL)
        endif
      endif

!----------------------------------------------------------------------


end subroutine convection_driver_time_vary


!########################################################################

subroutine convection_driver   &
                   (is, ie, js, je, Surf_diff, Phys_mp_exch, &
                       Moist_clouds_block, Input_mp, Tend_mp, C2ls_mp, &
                                          Output_mp, Removal_mp, Aerosol)

!------------------------------------------------------------------------
!    subroutine convection_driver saves needed variables on input for
!    later use, calls the driver routines for the activated convective
!    implementation to compute convective effects on the model variables,
!    computes cumulus momentum transport and other effects of convection,
!    and outputs various convective diagnostics.
!------------------------------------------------------------------------

integer,                intent(in)           :: is, ie, js, je
type(surf_diff_type),   intent(in)           :: Surf_diff
type(phys_mp_exch_type),intent(inout)        :: Phys_mp_exch
type(clouds_from_moist_block_type), &
                        intent(inout)        :: Moist_clouds_block
type(mp_input_type),    intent(inout)        :: Input_mp
type(mp_tendency_type), intent(inout)        :: Tend_mp
type(mp_conv2ls_type),  intent(inout)        :: C2ls_mp
type(mp_output_type),   intent(inout)        :: Output_mp
type(mp_removal_type),  intent(inout)        :: Removal_mp
type(aerosol_type),     intent(in), optional :: Aerosol

!-----------------------------------------------------------------------
!    is,ie      starting and ending i indices for window
!    js,je      starting and ending j indices for window
!    Surf_diff  derived type which will transport surface fluxes from
!               atmos_model to convection_driver  via physics_driver
!               and moist_processes (not yet implemented)
!    Phys_mp_exch 
!               derived type used to transfer data between physics_driver
!               and convection_driver via moist_processes
!    Moist_clouds_block 
!               derived type used to transfer cloud data between 
!               atmos_model and convection_driver via physics_driver and
!               moist_processes
!    Input_mp   derived type used to transfer needed input data between
!               moist_processes and convection_driver
!    Tend_mp    derived type used to transfer calculated tendency data
!               between convection_driver and moist_processes
!    C2ls_mp    derived type used to transfer data from convection_driver
!               to lscloud_driver via moist_processes.
!    Output_mp  derived type used to transfer output fields between
!               convection_driver and moist_processes
!    Removal_mp derived type used to transfer precipitation and tracer
!               removal fields between convection_driver and 
!               moist_processes
!    Aerosol    derived type containing model aerosol fields to be input
!               to model convective schemes
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------

      type(conv_results_type) :: Conv_results
      integer                 :: ix, jx, kx, nt

!-----------------------------------------------------------------------
!    Conv_results     derived type variable containing convective results
!                     that are collected during integration and passed to
!                     the diagnostics output subroutine.                  
!    ix, jx, kx       i, j, and k dimensions of the current physics_window
!    nt               number of activated prognostic tracers
!-----------------------------------------------------------------------

!---------------------------------------------------------------------
!    turn on the convective parameterizations clock.
!---------------------------------------------------------------------
      call mpp_clock_begin (convection_clock)

!---------------------------------------------------------------------
!    verify that the module has been initialized.
!---------------------------------------------------------------------
      if (.not. module_is_initialized) then
        call error_mesg ('convection_driver_mod',  &
                 'convection_driver_init has not been called.', FATAL)
      endif

!-----------------------------------------------------------------------
!    be sure optional argument is present if needed.
!-----------------------------------------------------------------------
      if (present(Aerosol)) then
      else
        if (do_uw_conv) then
          call error_mesg ('convection_driver_mod',  &
            'Aerosol argument required when either do_uw_conv with&
               &  prognostic droplet number is activated.', FATAL)
        endif
      endif

!-----------------------------------------------------------------------
!    save the value of rdt upon entry to the routine so that the total
!    convective contribution to the tendency (as defined by this module)
!    may later be defined.
!    save input temperature, specific humidity and tracer fields. the 
!    original values may be needed by multiple parameterizations called 
!    from within this module; however, these values may be updated within
!    this module, so those values would not be available.
!-----------------------------------------------------------------------
      Output_mp%rdt_init = Output_mp%rdt
      Input_mp%tin_orig = Input_mp%tin
      Input_mp%qin_orig = Input_mp%qin
      Input_mp%tracer_orig = Input_mp%tracer

!-----------------------------------------------------------------------
!    define array dimensions.
!-----------------------------------------------------------------------
      ix = size(Input_mp%t,1) 
      jx = size(Input_mp%t,2) 
      kx = size(Input_mp%t,3) 
      nt = size(Output_mp%rdt,4)

!--------------------------------------------------------------------
!    call convection_driver-alloc to allocate and initialize the components
!    of the conv_results_type variable Conv_results.
!--------------------------------------------------------------------
      call convection_driver_alloc (ix, jx, kx, Conv_results)


!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!&                                                                      &
!&               AVAILABLE CONVECTIVE IMPLEMENTATIONS                   &
!&                                                                      &
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

!------------------------------------------------------------------------
!    integrate the active convective implementation.
!------------------------------------------------------------------------     

!        @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!        @                                            @
!        @         UW CONVECTION SCHEME ONLY          @
!        @                                            @
!        @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!---------------------------------------------------------------------
!    if only the uw convection scheme is desired call subroutine 
!    uw_conv_driver.
!---------------------------------------------------------------------
      if (luwconv) then
        call uw_conv_driver    &
                ( is, ie, js, je, Input_mp, Aerosol, Phys_mp_exch,   &
                  Output_mp, Tend_mp, Conv_results, Removal_mp,      &
                                 Moist_clouds_block%cloud_data(i_shallow))
      end if
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!&                                                                      &
!&             END OF INDIVIDUAL CONVECTIVE SCHEMES                     &
!&                                                                      &
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&


!---------------------------------------------------------------------
!    now that all potential convection schemes have been processed, 
!    calculate cumulus momentum transport, if desired. 
!---------------------------------------------------------------------
      if (do_cmt) then

!----------------------------------------------------------------------
!    activate the cmt clock, call cu_mo_trans, and then deactivate 
!    the clock.
!----------------------------------------------------------------------
        call mpp_clock_begin (cmt_clock)
        call cu_mo_trans (is, js, Time, dt, num_prog_tracers, Input_mp, &
                         Conv_results, Output_mp, Tend_mp%ttnd_conv)
        call mpp_clock_end (cmt_clock)
      endif 

!------------------------------------------------------------------------
!    call define_total_convective_output to  a) define total cloud mass
!    flux, b) define cloud base and cloud top, c) define lightning
!    source of nox, d) define convective gustiness, e) set flag indicating
!    columns with convection, f) update the Input_mp%tracer field with
!    the changes produced in this module, and g) update the tin and
!    rdt fields with the uw tendencies if they have not already been so
!    updated. Updating fields at this point only occurs when legacy warsaw 
!    results are desired; otherwise the updates would have been done at
!    the more appropriate time.
!------------------------------------------------------------------------
      call define_total_convective_output (is, js, nt, C2ls_mp,  &
                          Conv_results, Input_mp, Output_mp, Phys_mp_exch) 

!------------------------------------------------------------------------
!    call convective_diagnostics to produce and output desired diagnostics
!    reflecting the model's total convection, summed over the active
!    convective parameterization(s).      
!------------------------------------------------------------------------
      call convective_diagnostics (is, js, C2ls_mp, Conv_results,  &
                                           Input_mp, Tend_mp, Output_mp)

!-----------------------------------------------------------------------
!    call define_convective_area to compute the grid box area taken up by 
!    convective clouds, so that this information may be supplied to the 
!    large-scale cloud module.  
!-----------------------------------------------------------------------
      call define_convective_area    &
                      (C2ls_mp, Moist_clouds_block, Input_mp)

!----------------------------------------------------------------------
!    define the interface-level precip fluxes needed for input to the 
!    COSP simulator package.
!---------------------------------------------------------------------
      if (do_cosp) then
        call define_inputs_for_cosp (Removal_mp)
      endif
!------------------------------------------------------------------------
!    deallocate the components of the Conv_results derived type variable.
!------------------------------------------------------------------------
      call convection_driver_dealloc (Conv_results)

!---------------------------------------------------------------------
!    end the timing of the convection code section.
!---------------------------------------------------------------------
      call mpp_clock_end (convection_clock)


!------------------------------------------------------------------------


end subroutine convection_driver


!######################################################################

subroutine convection_driver_endts 

!-----------------------------------------------------------------------
!    subroutine convection_driver_endts performs needed calculations 
!    upon exiting the prognostic loop.
!-----------------------------------------------------------------------

end subroutine convection_driver_endts

!########################################################################

subroutine convection_driver_end 

!--------------------------------------------------------------------- 
!    subroutine convection_driver_end calls destructor routines for the
!    modules initialized within this module, and deallocates module
!    variables.
!--------------------------------------------------------------------- 

!--------------------------------------------------------------------- 
!    call the destructor routines for the active convection modules.
!--------------------------------------------------------------------- 
      if (do_uw_conv    ) call uw_conv_end 
      if (do_cmt        ) call cu_mo_trans_end 
      call detr_ice_num_end 

!----------------------------------------------------------------------
!    deallocate module variables.
!----------------------------------------------------------------------
      deallocate (tracers_in_uw    )

      deallocate (cloud_tracer)

      deallocate(id_tracerdt_conv)        ! h1g, 2017-02-02
      deallocate (id_tracerdt_conv_col)   ! h1g, 2017-02-02
      deallocate (id_conv_tracer)         ! h1g, 2017-02-02
      deallocate (id_conv_tracer_col)     ! h1g, 2017-02-02

!--------------------------------------------------------------------

   
 end subroutine convection_driver_end



!######################################################################

subroutine convection_driver_restart (timestamp)

!------------------------------------------------------------------------
!    subroutine convection_driver_restart controls the writing of any 
!    restart files associated with activated convection schemes.
!------------------------------------------------------------------------

character(len=*), intent(in), optional :: timestamp

!-----------------------------------------------------------------------
!       timestamp   character string that represents the model time, 
!                   used for writing restart. timestamp will append to
!                   the any restart file name as a prefix.
!-----------------------------------------------------------------------

end subroutine convection_driver_restart


!#######################################################################

subroutine cape_cin_diagnostics (is, ie, js, je, Input_mp, Time)

!-----------------------------------------------------------------------
!    subroutine cape_cin_diagnostics calls subroutine capecalcnew to 
!    compute a parcel's ascent, computing cape and cin of the 
!    environment as it does, if diagnostics for cape or cin are requested, 
!-----------------------------------------------------------------------

integer,             intent(in) :: is, ie, js, je
type(mp_input_type), intent(in) :: Input_mp
type(time_type),     intent(in) :: Time

!-----------------------------------------------------------------------
!    is,ie      starting and ending i indices for window
!    js,je      starting and ending j indices for window
!    Input_mp   derived type variable containing model profiles
!    Time       variable containing current diagnostic time [ time_type ]
!-----------------------------------------------------------------------

      real, dimension(size(Input_mp%tin,1), size(Input_mp%tin,2),  &
                                            size(Input_mp%tin,3)) ::  &
                                                             rin, rp, tp
      real, dimension(size(Input_mp%tin,1), size(Input_mp%tin,2)) ::  &
                                                                cape, cin
      integer, dimension(size(Input_mp%tin,1), size(Input_mp%tin,2)) ::  &
                                                       klcl, klfc, klzb

      logical :: avgbl, used
      integer :: i, j, ix, jx, kx
 
!--------------------------------------------------------------------
!    rin          model h2o mixing ratio [ kg [h2o] / kg [dry air ] 
!    rp           rising parcel mixing ratio profile
!    tp           rising parcel temperature profile
!    cape         convective available potential energy
!    cin          convective inhibition
!    klcl         model lifting condensation level for rising parcel
!    klfc         model level of free convection for rising parcel
!    klzb         model level of zero buoyancy for rising parcel
!    avgbl        outdated variable no longer used in capecalcnew
!    used         logical used to indicate data has been received by
!                 diag_manager_mod
!    i, j         do loop indices
!    ix, jx, kx   array  spatial dimensions; size of physics_window
!------------------------------------------------------

!-------------------------------------------------------------------
!    proceed with computation if diagnostics for cape or cin are requested,
!------------------------------------------------------------------
      if ( id_cape > 0 .or. id_cin > 0) then

!---------------------------------------------------------------------
!    define physics window dimensions.
!--------------------------------------------------------------------
        kx = size(Input_mp%tin,3)
        ix = size(Input_mp%tin,1)
        jx = size(Input_mp%tin,2)

!----------------------------------------------
!    calculate mixing ratio.
!----------------------------------------------
        rin = Input_mp%qin/(1.0 - Input_mp%qin) 

!-----------------------------------------------------------------------
!    call routine to calculate cape and cin based on parcel rise.
!-----------------------------------------------------------------------
        avgbl = .false.
        do j = 1,jx 
          do i = 1,ix 
            call capecalcnew   &
                ( kx, Input_mp%pfull(i,j,:), Input_mp%phalf(i,j,:),   &
                  CP_AIR, RDGAS, RVGAS, HLV, KAPPA, Input_mp%tin(i,j,:), &
                  rin(i,j,:), avgbl, cape(i,j), cin(i,j), tp(i,j,:), &
                  rp(i,j,:), klcl(i,j), klfc(i,j), klzb(i,j))
          end do
        end do

!-------------------------------------------------------------------------
!    output any requested diagnostics.
!-------------------------------------------------------------------------
        if (id_cape > 0) used = send_data ( id_cape, cape, Time, is, js )
        if ( id_cin > 0 ) used = send_data ( id_cin, cin, Time, is, js )
        if ( id_tp  > 0 ) used = send_data ( id_tp,  tp, Time, is, js )
        if ( id_rp  > 0 ) used = send_data ( id_rp,  rp, Time, is, js )
        if ( id_lcl > 0 ) used = send_data ( id_lcl, 1.0*klcl, Time,  &
                                                                 is, js )
        if ( id_lfc > 0 ) used = send_data ( id_lfc, 1.0*klfc, Time,  &
                                                                 is, js )
        if ( id_lzb > 0 ) used = send_data ( id_lzb, 1.0*klzb, Time,  &
                                                                 is, js )
      end if

!-----------------------------------------------------------------------



end subroutine cape_cin_diagnostics



!*******************************************************************
!
!                     PRIVATE INITIALIZATION SUBROUTINES
!
!*******************************************************************



!#########################################################################

subroutine diag_field_init ( axes, Time, Control)

!-----------------------------------------------------------------------
!    this subroutine initializes diagnostic fields from this module. it
!    also initializes global integrals for netCDF output.
!-----------------------------------------------------------------------

integer,                       intent(in) :: axes(4)
type(time_type),               intent(in) :: Time
type(mp_removal_control_type), intent(in) :: Control

!------------------------------------------------------------------------

!------------------------------------------------------------------------
!      axes        data axis indices, (x,y,pf,ph) for diagnostics 
!      Time        time used for diagnostics [time_type]
!      Control     derived type variable containing control variables
!                  associated with tracer removal and transport by
!                  available convective schemes
!------------------------------------------------------------------------

      character(len=32)     :: tracer_units, tracer_name
      character(len=128)    :: diaglname
      integer, dimension(3) :: half = (/1,2,4/)
      integer               :: n, nn

!-----------------------------------------------------------------------
!      tracer_units  units assigned to each tracer field
!      tracer_name   name assigned to each tracer field
!      diaglname     long name associated with each tracer diagnostic field
!      half          axis indices for x, y and model half-levels
!      n             do loop index
!      nn            counter for subset of tracers meeting a particular
!                    condition
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!    initialize global integrals for netCDF output.
!-----------------------------------------------------------------------
      id_pr_g = register_global_diag_field ('pr', Time, 'Precipitation', &
                  'kg m-2 s-1', standard_name='precipitation_flux',   &
                                                           buffer=.true. )
      id_prc_g = register_global_diag_field ('prc', Time,    &
                  'Convective Precipitation', 'kg m-2 s-1',   &
                        standard_name='convective_precipitation_flux', &
                                                           buffer=.true. )
      id_prsn_g = register_global_diag_field ('prsn', Time,    &
                  'Snowfall Flux', 'kg m-2 s-1', &
                            standard_name='snowfall_flux', buffer=.true. )

!-------------------------------------------------------------------------
!    diagnostics related to total convective tendencies of temperature,
!    vapor and precipitation.
!-------------------------------------------------------------------------
      id_tdt_conv = register_diag_field ( mod_name, &
                   'tdt_conv', axes(1:3), Time, &
                   'Temperature tendency from convection ',  'deg_K/s',  &
                                           missing_value=missing_value)

      ID_tntc = register_cmip_diag_field_3d ( mod_name, 'tntc', Time, &
                  'Tendency of Air Temperature Due to Convection ',   &
                      'K s-1', standard_name=      &
                          'tendency_of_air_temperature_due_to_convection' )

      id_qdt_conv = register_diag_field ( mod_name, &
                  'qdt_conv', axes(1:3), Time, &
                  'Spec humidity tendency from convection ',  'kg/kg/s',  &
                                            missing_value=missing_value)

      ID_tnhusc = register_cmip_diag_field_3d ( mod_name, 'tnhusc', Time, &
                  'Tendency of Specific Humidity Due to Convection ',   &
                     's-1', standard_name=   &
                        'tendency_of_specific_humidity_due_to_convection' )

      id_q_conv_col = register_diag_field ( mod_name, &
                      'q_conv_col', axes(1:2), Time, &
                      'Water vapor path tendency from convection ',  &
                                                              'kg/m2/s' )
   
      id_t_conv_col = register_diag_field ( mod_name, &
                      't_conv_col', axes(1:2), Time, &
                      'Column static energy tendency from convection ', &
                                                                'W/m2' )
   
      id_enth_conv_col = register_diag_field ( mod_name, &
                         'enth_conv_col', axes(1:2), Time, &
                         'Column enthalpy tendency from convection',  &
                                                                'W/m2' )
 
      id_wat_conv_col = register_diag_field ( mod_name, &
                        'wat_conv_col', axes(1:2), Time, &
                        'Column total water tendency from convection', &
                                                        'kg(h2o)/m2/s' )

      id_prec_conv = register_diag_field ( mod_name, &
                     'prec_conv', axes(1:2), Time, &
                     'Precipitation rate from convection ',    &
                     'kg(h2o)/m2/s', interp_method = "conserve_order1" )

      id_prc = register_cmip_diag_field_2d ( mod_name, 'prc', Time, &
                        'Convective Precipitation',   'kg m-2 s-1', &
                    standard_name = 'convective_precipitation_flux', &
                                  interp_method = "conserve_order1" ) 
 
      id_prrc = register_cmip_diag_field_2d ( mod_name, 'prrc', Time, &
                             'Convective Rainfall Rate', 'kg m-2 s-1', &
                             standard_name='convective_rainfall_flux', &
                                       interp_method="conserve_order1" )

      id_snow_conv = register_diag_field ( mod_name, &
                     'snow_conv', axes(1:2), Time, &
                     'Frozen precip rate from convection ',  &
                     'kg(h2o)/m2/s', interp_method = "conserve_order1" )

      id_conv_freq = register_diag_field ( mod_name, &
                     'conv_freq', axes(1:2), Time, &
                     'frequency of convection ',       '1', &
                     missing_value = missing_value, &
                      interp_method = "conserve_order1"      )

      id_prsnc = register_cmip_diag_field_2d ( mod_name, 'prsnc', Time, &
                              'Convective Snowfall Flux', 'kg m-2 s-1', &
                               standard_name='convective_snowfall_flux', &
                                         interp_method="conserve_order1" )

      id_ci = register_cmip_diag_field_2d ( mod_name, 'ci', Time, &
                 'Fraction of Time Convection Occurs in Cell',  '1.0',  &
                             standard_name='convection_time_fraction', &
                             interp_method = "conserve_order1" )

      id_gust_conv = register_diag_field ( mod_name, &
                     'gust_conv', axes(1:2), Time, &
                          'Gustiness resulting from convection ', 'm/s' )

      id_conv_rain3d= register_diag_field ( mod_name, &
                   'conv_rain3d', axes(half), Time, &
                      'Rain fall rate from convection -3D ',    &
                        'kg(h2o)/m2/s', interp_method = "conserve_order1" )

      id_conv_snow3d= register_diag_field ( mod_name, &
                    'conv_snow3d', axes(half), Time, &
                     'Snow fall rate from convection -3D',   &
                      'kg(h2o)/m2/s', interp_method = "conserve_order1" )

!----------------------------------------------------------------------
!    tendencies of cloud tracers resulting from convection.
!----------------------------------------------------------------------
      if (doing_prog_clouds ) then

        id_qldt_conv = register_diag_field ( mod_name, &
                       'qldt_conv', axes(1:3), Time, &
                           'Liquid water tendency from convection',  &
                               'kg/kg/s', missing_value=missing_value  )

        if (do_liq_num) then
          id_qndt_conv = register_diag_field ( mod_name, &
                          'qndt_conv', axes(1:3), Time, &
                            'Liquid drop tendency from convection',  &
                                 '#/kg/s', missing_value=missing_value)
        endif

        id_qidt_conv = register_diag_field ( mod_name, &
                         'qidt_conv', axes(1:3), Time, &
                            'Ice water tendency from convection',   &
                                'kg/kg/s', missing_value=missing_value )

        id_qadt_conv = register_diag_field ( mod_name, &
                         'qadt_conv', axes(1:3), Time, &
                           'Cloud fraction tendency from convection',  &
                              '1/sec', missing_value=missing_value )

        id_ql_conv_col = register_diag_field ( mod_name, &
                          'ql_conv_col', axes(1:2), Time, &
                           'Liquid water path tendency from convection',  &
                                                              'kg/m2/s' )
   
        if (do_liq_num) then
          id_qn_conv_col = register_diag_field ( mod_name, &
                            'qn_conv_col', axes(1:2), Time, &
                               'Liquid drp tendency from convection',  &
                                                               'kg/m2/s' )
        endif
 
        id_qi_conv_col = register_diag_field ( mod_name, &
                           'qi_conv_col', axes(1:2), Time, &
                            'Ice water path tendency from convection',  &
                                                              'kg/m2/s' )
   
        id_qa_conv_col = register_diag_field ( mod_name, &
                           'qa_conv_col', axes(1:2), Time, &
                             'Cloud mass tendency from convection',   &
                                                               'kg/m2/s' )
      
        if (do_ice_num) then
          id_qnidt_conv = register_diag_field ( mod_name, &
                            'qnidt_conv', axes(1:3), Time, &
                              'Ice number tendency from convection',   &
                                  '#/kg/s', missing_value=missing_value )

          id_qni_conv_col = register_diag_field ( mod_name, &
                              'qni_conv_col', axes(1:2), Time, &
                               'Ice number tendency from convection',   &
                                                               'kg/m2/s' )
        endif
      endif ! (doing_prog_clouds)

!-----------------------------------------------------------------------
!    diagnostics for cloud base and cloud top.
!-----------------------------------------------------------------------
      id_conv_cld_base = register_diag_field ( mod_name, &
                            'conv_cld_base', axes(1:2), Time, &
                               'pressure at convective cloud base',  &
                                  'Pa', mask_variant = .true., &
                                        missing_value=missing_value )

      id_ccb = register_cmip_diag_field_2d ( mod_name, 'ccb', Time, &
                  'Air Pressure at Convective Cloud Base', 'Pa', &
                        standard_name =    &
                               'air_pressure_at_convective_cloud_base', &
                                                   mask_variant = .true. )

      id_conv_cld_top = register_diag_field ( mod_name, &
                          'conv_cld_top', axes(1:2), Time, &
                           'pressure at convective cloud top',   'Pa', &
                               mask_variant = .true., &
                                     missing_value=missing_value )

      id_cct = register_cmip_diag_field_2d ( mod_name, 'cct', Time, &
                 'Air Pressure at Convective Cloud Top', 'Pa', &
                      standard_name =   &
                              'air_pressure_at_convective_cloud_top', &
                                                  mask_variant = .true. )

!-----------------------------------------------------------------------
!    convective mass flux diagnostics.
!-----------------------------------------------------------------------
      id_mc_full = register_diag_field ( mod_name, &
                      'mc_full', axes(1:3), Time, &
                         'Net Mass Flux from convection',   'kg/m2/s', &
                                           missing_value=missing_value )

      id_mc_half = register_diag_field ( mod_name, &
                     'mc_half', axes(half), Time, &
                       'Net Mass Flux from convection on half levs',   &
                              'kg/m2/s', missing_value=missing_value )

      ID_mc = register_cmip_diag_field_3d ( mod_name, 'mc', Time, &
                    'Convective Mass Flux',   'kg m-2 s-1', &
                        standard_name=  &
                         'atmosphere_net_upward_convective_mass_flux', &
                           interp_method = "conserve_order1", axis="half" )
   
      id_mc_conv_up = register_diag_field ( mod_name, &
                        'mc_conv_up', axes(1:3), Time, &
                           'Upward Mass Flux from convection', 'kg/m2/s', &
                                            missing_value=missing_value )

!---------------------------------------------------------------------
!    register diagnostics for lightning NOx.
!---------------------------------------------------------------------
      if (get_tracer_index(MODEL_ATMOS,'no') > 0) then
        id_prod_no = register_diag_field ( 'tracers', &
                                         'hook_no', axes(1:3), Time, &
                                               'hook_no',   'molec/cm3/s')
        ID_emilnox_area = register_cmip_diag_field_3d ( mod_name,  &
                        'emilnox_area', Time, &
           'Layer-integrated Lightning Production of NOx', 'mol m-2 s-1', &
             standard_name=   &
              'tendency_of_atmosphere_moles_of_nox_expressed_as_nitrogen')
      end if

!-------------------------------------------------------------------------
!    register diagnostics associated with CAPE / CIN calculations.
!-------------------------------------------------------------------------
      id_cape = register_diag_field ( mod_name, &
                    'cape', axes(1:2), Time, &
                      'Convectively available potential energy', 'J/Kg')
      
      id_cin = register_diag_field ( mod_name, &
                    'cin', axes(1:2), Time, 'Convective inhibition','J/Kg')

      id_tp = register_diag_field ( mod_name, &
              'tp', axes(1:3), Time, 'Temperature of lifted parcel', 'K')

      id_rp = register_diag_field ( mod_name, &
              'rp', axes(1:3), Time, 'Humidity of lifted parcel', 'kg/kg')

      id_lcl = register_diag_field ( mod_name, &
               'klcl', axes(1:2), Time, 'Index of LCL', 'none')

      id_lfc = register_diag_field ( mod_name, &
               'klfc', axes(1:2), Time, 'Index of LFC', 'none')

      id_lzb = register_diag_field ( mod_name, &
               'klzb', axes(1:2), Time, 'Index of LZB', 'none') 

!------------------------------------------------------------------------
!    register diagnostics specific to the uw parameterization.
!------------------------------------------------------------------------
      if (do_uw_conv) then
        id_uw_precip = register_diag_field ( mod_name, &
                       'uw_precip', axes(1:2), Time, &
                       'Precipitation rate from uw shallow',  'kg/m2/s', &
                       interp_method = "conserve_order1" )

        id_uw_snow = register_diag_field ( mod_name, &
                     'uw_snow', axes(1:2), Time, &
                     'Snow rate from uw shallow',       'kg/m2/s' , &
                     interp_method = "conserve_order1" )

        id_uw_freq = register_diag_field ( mod_name, &
                     'uw_freq', axes(1:2), Time, &
                     'frequency of precip from uw shallow ',  'number' , &
                           missing_value = missing_value,   &
                                   interp_method = "conserve_order1"  )

        id_enth_uw_col = register_diag_field ( mod_name, &
                         'enth_uw_col', axes(1:2), Time, &
                         'Column enthalpy tendency from UW convection', &
                         'W/m2' )
 
        id_wat_uw_col = register_diag_field ( mod_name, &
                        'wat_uw_col', axes(1:2), Time, &
                        'Column total water tendency from UW convection',&
                        'kg(h2o)/m2/s' )

        id_scale_uw = register_diag_field ( mod_name, &
                      'scale_uw', axes(1:2), Time, &
                      'Scaling factor applied to UW convection&
                      & tendencies','1' )
          
        id_tdt_uw = register_diag_field ( mod_name, &
                    'tdt_uw', axes(1:3), Time, &
                    'UW convection heating rate', 'deg K/s', &
                    missing_value=missing_value               )

        id_qdt_uw = register_diag_field ( mod_name, &
                    'qdt_uw', axes(1:3), Time, &
                    'UW convection moistening rate', 'kg/kg/s', &
                    missing_value=missing_value               )

        id_qadt_uw = register_diag_field ( mod_name, &
                     'qadt_uw', axes(1:3), Time, &
                     'UW convection cloud amount tendency', '1/s', &
                     missing_value=missing_value               )

        id_qldt_uw = register_diag_field ( mod_name, &
                     'qldt_uw', axes(1:3), Time, &
                     'UW convection cloud liquid tendency', 'kg/kg/s', &
                     missing_value=missing_value               )

        id_qidt_uw = register_diag_field ( mod_name, &
                     'qidt_uw', axes(1:3), Time, &
                     'UW convection ice water tendency', 'kg/kg/s', &
                     missing_value=missing_value               )

        if (do_liq_num) &
          id_qndt_uw = register_diag_field ( mod_name, &
                       'qndt_uw', axes(1:3), Time, &
                       'UW convection cloud drop tendency', '#/kg/s', &
                       missing_value=missing_value               )

        if (do_ice_num) &
          id_qnidt_uw = register_diag_field ( mod_name, &
                        'qnidt_uw', axes(1:3), Time, &
                        'UW convection ice number tendency', '#/kg/s', &
                        missing_value=missing_value               )

      endif

!---------------------------------------------------------------------
!    allocate and initialize arrays to hold the diagnostic ids for each 
!    active tracer. diagnostics for tendency due to convection, 
!    column tendency due to convection, the tracer amount and tracer 
!    column amount are available.
!---------------------------------------------------------------------
      allocate (id_tracerdt_conv     (num_prog_tracers))
      allocate (id_tracerdt_conv_col (num_prog_tracers))
      allocate (id_conv_tracer       (num_prog_tracers))
      allocate (id_conv_tracer_col   (num_prog_tracers))

      id_tracerdt_conv = NO_TRACER
      id_tracerdt_conv_col = NO_TRACER
      id_conv_tracer = NO_TRACER
      id_conv_tracer_col = NO_TRACER
 
!------------------------------------------------------------------------
!    define the diagnostics names that are requested and register the
!    diagnostics for those tracers that were specified to be affected
!    by a convection scheme.
!------------------------------------------------------------------------
      do n = 1,num_prog_tracers
        call get_tracer_names (MODEL_ATMOS, n, name = tracer_name,  &
                                                  units = tracer_units)
        if (Control%tracers_in_uw(n)) then
          diaglname = trim(tracer_name)//  &
                        ' total tendency from moist convection'
          id_tracerdt_conv(n) =    &
                         register_diag_field ( mod_name, &
                         TRIM(tracer_name)//'dt_conv',  &
                         axes(1:3), Time, trim(diaglname), &
                         TRIM(tracer_units)//'/s',  &
                         missing_value=missing_value)

          diaglname = trim(tracer_name)//  &
                       ' total path tendency from moist convection'
          id_tracerdt_conv_col(n) =  &
                         register_diag_field ( mod_name, &
                         TRIM(tracer_name)//'dt_conv_col', &
                         axes(1:2), Time, trim(diaglname), &
                         TRIM(tracer_units)//'*(kg/m2)/s',   &
                         missing_value=missing_value)
        endif
 
!----------------------------------------------------------------------
!    output the distribution and column values of any tracer for which
!    they are requested, even if not transported by convection.
!----------------------------------------------------------------------
        diaglname = trim(tracer_name)

!---------------------------------------------------------------------
!    RSH:
!    temporary get-around for the fact that 'cl' may be both a tracer 
!    variable (full chemistry) and a CMIP6 cloud diagnostic variable 
!    in module 'moist', and so they need to be registered differently. 
!    In the future, all the tracers should be registered under mod_name2,
!    after existing scripts / experiments are replaced. 
!---------------------------------------------------------------------
        if (trim(diaglname) == 'cl') then
          id_conv_tracer(n) =    &
                        register_diag_field ( mod_name2,    &
                        TRIM(tracer_name),  &
                        axes(1:3), Time, trim(diaglname), &
                        TRIM(tracer_units)      ,  &
                        missing_value=missing_value)
        else
          id_conv_tracer(n) =    &
                        register_diag_field ( mod_name,    &
                        TRIM(tracer_name),  &
                        axes(1:3), Time, trim(diaglname), &
                        TRIM(tracer_units)      ,  &
                        missing_value=missing_value)
        endif
        diaglname =  ' column integrated' // trim(tracer_name)
        id_conv_tracer_col(n) =  &
                        register_diag_field ( mod_name, &
                        TRIM(tracer_name)//'_col', &
                        axes(1:2), Time, trim(diaglname), &
                        TRIM(tracer_units)      ,   &
                        missing_value=missing_value)
      end do

!---------------------------------------------------------------------


end subroutine diag_field_init




!*******************************************************************
!
!                     PRIVATE DRIVER-CALLED SUBROUTINES,
!                       NON-CONVECTION-SCHEME SPECIFIC
!
!*******************************************************************


!#######################################################################

subroutine convection_driver_alloc (ix, jx, kx, Conv_results)

!---------------------------------------------------------------------
!    subroutine convection_driver_alloc allocates and initializes
!    variables that are calculated in this module and transferred
!    between convective parameterizations and / or to the diagnostics 
!    routine for output.
!---------------------------------------------------------------------

integer,                 intent(in)      :: ix, jx, kx
type(conv_results_type), intent(inout)   :: Conv_results

!--------------------------------------------------------------------
!      ix, jx, kx      dimensions of physics window
!      Conv_results    conv_results_type variable containing local 
!                      variables used in multiple convective 
!                      parameterizations and for diagnostic output
!--------------------------------------------------------------------

!--------------------------------------------------------------------
!    allocate and initialize arrays which transfer data between individual
!    parameterizations within an implementation.
!--------------------------------------------------------------------

      allocate (Conv_results%conv_calc_completed  (ix, jx))
      allocate (Conv_results%available_cf_for_uw  (ix, jx, kx))
      Conv_results%conv_calc_completed = .false.
      Conv_results%available_cf_for_uw = 1.0

!------------------------------------------------------------------------
!    allocate and initialize the massflux-related components of the 
!    conv_results_type variable Conv_results.
!------------------------------------------------------------------------
      allocate (Conv_results%uw_mflux(ix, jx, kx+1))  ! cmf
      Conv_results%uw_mflux = 0.

!------------------------------------------------------------------------
!    allocate the components of the conv_results_type variable Conv_results
!    which depend on the results of the convective implementation, not each
!    convective parameterization included within it.
!------------------------------------------------------------------------
      allocate (Conv_results%cldtop        (ix, jx)) 
      allocate (Conv_results%cldbot        (ix, jx))
      allocate (Conv_results%prod_no       (ix, jx, kx))
      Conv_results%cldtop = 0
      Conv_results%cldbot = 0
      Conv_results%prod_no = 0.

!----------------------------------------------------------------------


end subroutine convection_driver_alloc



!######################################################################

subroutine define_total_convective_output    &
              (is, js, nt, C2ls_mp, Conv_results, Input_mp, Output_mp,   &
                                                             Phys_mp_exch)

!----------------------------------------------------------------------
!    subroutine define_total_convective_output: a) defines total cloud 
!    mass flux, b) defines cloud base and cloud top, c) defines lightning
!    source of nox, d) defines convective gustiness, e) sets flag 
!    indicating columns with convection, f) updates the Input_mp%tracer 
!    field with the changes produced in this module, and g) updates the 
!    tin and rdt fields with the uw tendencies if they have not already 
!    been so updated. Updating fields at this point only occurs when 
!    legacy warsaw results are desired; otherwise the updates would have 
!    been done at the more appropriate earlier time in the integration.
!----------------------------------------------------------------------

!-----------------------------------------------------------------------
integer,                  intent(in)    :: is, js, nt
type(mp_conv2ls_type),    intent(inout) :: C2ls_mp
type(conv_results_type),  intent(inout) :: Conv_results
type(mp_input_type),      intent(inout) :: Input_mp
type(mp_output_type),     intent(inout) :: Output_mp
type(phys_mp_exch_type),  intent(inout) :: Phys_mp_exch
!-----------------------------------------------------------------------

!----------------------------------------------------------------------
!    is, js     starting i,j physics window indices
!    nt         number of prognostic tracers
!    C2ls_mp    derived type used to transfer data from 
!               convection_driver to lscloud_driver via moist_processes.
!    Conv_results  
!               conv_results_type variable containing local variables 
!               used in multiple convective parameterizations and for
!               diagnostic output
!    Input_mp   derived type used to transfer needed input data between
!               moist_processes and convection_driver
!    Output_mp  derived type used to transfer output fields between
!               convection_driver and moist_processes
!    Phys_mp_exch
!               derived type used to transfer data between physics_driver
!               and convection_driver via moist_processes
!----------------------------------------------------------------------

      real, parameter :: boltz = 1.38044e-16
      integer :: i, j, k, n
      integer :: ix, jx, kx

!---------------------------------------------------------------------
!       boltz           boltzmann constant
!       i, j, k, n      do-loop indices
!       ix, jx, kx      physics window dimensions
!---------------------------------------------------------------------

!------------------------------------------------------------------------
!    define array dimensions.
!------------------------------------------------------------------------
      ix = size(Input_mp%tin,1)
      jx = size(Input_mp%tin,2)
      kx = size (Input_mp%t, 3)

!------------------------------------------------------------------------
!    define total convective mass flux from all sources, at both full
!    levels and at half levels.
!------------------------------------------------------------------------
      C2ls_mp%mc_full(:,:,1)=0.; 
      C2ls_mp%mc_half(:,:,1)=0.; 
      do k=2,kx   
        C2ls_mp%mc_full(:,:,k) = 0.5*(Conv_results%uw_mflux(:,:,k)+   &
                                      Conv_results%uw_mflux(:,:,k-1)) 
      end do
      do k=2,kx+1   
        C2ls_mp%mc_half(:,:,k) = Conv_results%uw_mflux(:,:,k-1)
      end do

!------------------------------------------------------------------------ 
!    define convective cloud base and cloud top. these are needed if 
!    diagnostics defining the temporal and spatial location of convection 
!    are desired or if tracer "no" is present, so that the nox tendency 
!    due to lightning may be calculated.
!------------------------------------------------------------------------ 
      if (get_tracer_index(MODEL_ATMOS,'no') .ne. NO_TRACER &
           .or. id_conv_freq > 0 &
           .or. id_ci > 0 &
           .or. id_conv_cld_base > 0 &
           .or. id_ccb           > 0 &
           .or. id_cct           > 0 &
           .or. id_conv_cld_top > 0 ) then

        do j=1,jx
          do i=1,ix
            do k=1,kx
              if (C2ls_mp%mc_full(i,j,k) /= 0.0 ) then
                Conv_results%cldtop(i,j) = k
                exit
              endif
            enddo
            do k = size(Input_mp%r,3),1,-1
              if (C2ls_mp%mc_full(i,j,k) /= 0.0 ) then
                Conv_results%cldbot(i,j) = k
                exit
              endif
            enddo
          enddo
        enddo
      end if

!-----------------------------------------------------------------------
!    calculate NOx tendency from lightning and add it to the tendency
!    field.  
!-----------------------------------------------------------------------
      if ( get_tracer_index(MODEL_ATMOS,'no') .ne. NO_TRACER ) then
        call moz_hook       &
              (Conv_results%cldtop, Conv_results%cldbot, Input_mp%land, &
               Input_mp%zfull, Input_mp%zhalf, Input_mp%t,   &
               Conv_results%prod_no, Input_mp%area, Input_mp%lat,   &
               Time, is, js)

        Output_mp%rdt(:,:,:,get_tracer_index(MODEL_ATMOS,'no')) =  &
              Output_mp%rdt(:,:,:,get_tracer_index(MODEL_ATMOS,'no')) + &
                  Conv_results%prod_no*((boltz*Input_mp%t)/      &
                                                  (10. * Input_mp%pfull)) 
      endif

!-----------------------------------------------------------------------
!    calculate convective gustiness, if desired. two forms are available.
!-----------------------------------------------------------------------
      if (do_gust_cv) then
        where((Output_mp%precip) > 0.0)
          Output_mp%gust_cv = gustmax*   &
                    sqrt(Output_mp%precip/(gustconst + Output_mp%precip))
        end where
      end if

      if (do_gust_cv_new) then
        Output_mp%gust_cv = sqrt(Phys_mp_exch%cgust)
      end if

!---------------------------------------------------------------------
!    save a field indicating whether or not convection has occurred
!    within the column.
!---------------------------------------------------------------------
      where (Output_mp%precip > 0.) Output_mp%convect = .true.

!------------------------------------------------------------------------
!    update the current temp and rdt tendencies with the contributions 
!    obtained from uw transport (tin_tentative, rdt_tentative), if  
!    reproduce_AM4 was set true. otherwise these fields have 
!    already been updated, and the _tentative fields contain 0.0.
!------------------------------------------------------------------------
      Input_mp%tin = Input_mp%tin + Input_mp%tin_tentative
      Output_mp%rdt = Output_mp%rdt + Output_mp%rdt_tentative

!---------------------------------------------------------------------
!    update the input tracer arrays with the tendencies obtained in this
!    module.
!---------------------------------------------------------------------
      do n=1,nt                       
        if (.not. cloud_tracer(n)) then
          Input_mp%tracer(:,:,:,n) = Input_mp%tracer_orig(:,:,:,n) +   &
                 (Output_mp%rdt(:,:,:,n) - Output_mp%rdt_init(:,:,:,n))*dt
        endif
      end do

!------------------------------------------------------------------------


end subroutine define_total_convective_output



!######################################################################

subroutine convective_diagnostics    &
             (is, js, C2ls_mp, Conv_results, Input_mp, Tend_mp, Output_mp)

!-----------------------------------------------------------------------
!    subroutine convective_diagnostics outputs requested convective
!    diagnostics associated with the convective implementation rather than
!    any specific convective parameterization.
!-----------------------------------------------------------------------

integer,                     intent(in)    :: is, js
type(mp_conv2ls_type),       intent(inout) :: C2ls_mp
type(conv_results_type),     intent(inout) :: Conv_results
type(mp_input_type),         intent(inout) :: Input_mp
type(mp_tendency_type),      intent(inout) :: Tend_mp
type(mp_output_type),        intent(in)    :: Output_mp

!--------------------------------------------------------------------
!    is, js     starting i,j physics window indices
!    C2ls_mp    derived type used to transfer data from convection_driver
!               to lscloud_driver via moist_processes.
!    Conv_results
!               conv_results_type variable containing local variables 
!               used in multiple convective parameterizations and for
!               diagnostic output
!    Input_mp   derived type used to transfer needed input data between
!               moist_processes and convection_driver
!    Tend_mp    derived type used to transfer calculated tendency data
!               between convection_driver and moist_processes
!    Output_mp  derived type used to transfer output fields between
!               convection_driver and moist_processes
!--------------------------------------------------------------------

      real, dimension(size(Input_mp%t,1),size(Input_mp%t,2)) ::      &
                                                  temp_2d, freq_count
      logical, dimension(size(Input_mp%t,1),size(Input_mp%t,2)) ::      &
                                                               ltemp
      real, dimension(size(Input_mp%t,1),size(Input_mp%t,2), &
                           size(Input_mp%t,3)) ::      &
                                  temp_3d1, temp_3d2, uw_massflx_full
      integer :: i,j,k, n
      integer :: ix, jx, kx
      logical :: used

!--------------------------------------------------------------------
!      temp_2d            array used to hold diagnostic fields when
!                         they are sent for output
!      freq_count         array used to hold diagnostic field when sent
!                         for output
!      ltemp              logical used to define condition needed for
!                         several diagnostics
!      temp_3d1           temporary array used in computing diagnostics
!      temp_3d2           temporary array used in computing diagnostics
!      uw_massflx_full    uw massflux on full levels  [ kg/m2/s ]
!      i, j, k, n         do loop indices
!      ix, jx, kx         array spatial dimensions, size of physics window
!      used               logical used to indicate data has been received 
!                         by diag_manager_mod
!--------------------------------------------------------------------

!-------------------------------------------------------------------
!    define array dimensions.
!-------------------------------------------------------------------
      ix = size(Input_mp%r,1)
      jx = size(Input_mp%r,2)
      kx = size(Input_mp%r,3)

!-----------------------------------------------------------------------
!    output the NOx tendency from lightning.
!-----------------------------------------------------------------------
      if ( get_tracer_index(MODEL_ATMOS,'no') .ne. NO_TRACER ) then
        used = send_data(id_prod_no, Conv_results%prod_no, Time, is, js)
        used = send_cmip_data_3d (ID_emilnox_area, & ! convert molec/cm3/s to mol/m2/s
                                  Conv_results%prod_no*1.e6*(Input_mp%zhalf(:,:,1:kx)-Input_mp%zhalf(:,:,2:kx+1))/AVOGNO, &
                                  Time, is, js, 1, phalf=log(Input_mp%phalf))
      endif

!----------------------------------------------------------------------
!    output any requested convectively-transported tracer fields 
!    and / or their column sums before convective transport.
!----------------------------------------------------------------------
      do n=1, num_prog_tracers
        used = send_data (id_conv_tracer(n),   &
                           Input_mp%tracer_orig(:,:,:,n), Time, is, js, 1)
        if (id_conv_tracer_col(n) > 0)  &
          call column_diag(id_conv_tracer_col(n), is, js, Time, &
                       Input_mp%tracer_orig(:,:,:,n), 1.0, Input_mp%pmass) 
      end do

!---------------------------------------------------------------------
!    output diagnostics:
!    total cumulus mass flux on full levels,
!    total cumulus mass flux on half levels,
!    total cumulus mass flux on half levels (CMOR standard).
!---------------------------------------------------------------------
      used = send_data (id_mc_full, C2ls_mp%mc_full, Time, is, js, 1)
      used = send_data (id_mc_half, C2ls_mp%mc_half, Time, is, js, 1)
      used = send_cmip_data_3d (ID_mc, C2ls_mp%mc_half, Time, is, js, 1)

!---------------------------------------------------------------------
!    total convective updraft mass flux on full levels.
!---------------------------------------------------------------------
      if (id_mc_conv_up > 0 ) then
        do k=1,kx
          uw_massflx_full(:,:,k) = 0.5*(Conv_results%uw_mflux(:,:,k) + &
                                         Conv_results%uw_mflux(:,:,k+1))
        end do
        used = send_data (id_mc_conv_up, uw_massflx_full(:,:,:), Time, is, js, 1 )
      endif

!------------------------------------------------------------------------
!    output diagnostics related to convective cloud base and cloud top.
!    both FMS-standard and CMOR-standard output variables are currently
!    present.
!------------------------------------------------------------------------
      if ( id_conv_cld_base > 0 ) then
        temp_2d = missing_value
        do j = 1,jx
          do i = 1,ix
            if ( Conv_results%cldbot(i,j) > 0 ) temp_2d(i,j) =    &
                              Input_mp%pfull(i,j,Conv_results%cldbot(i,j))
          end do
        end do
        used = send_data(id_conv_cld_base, temp_2d, Time, is_in=is,   &
                               js_in=js,  mask = Conv_results%cldbot > 0)
      end if

      if ( id_ccb > 0 ) then
        temp_2d = CMOR_MISSING_VALUE
        do j = 1,jx  
          do i = 1,ix
            if ( Conv_results%cldbot(i,j) > 0 ) temp_2d(i,j) =    &
                             Input_mp%pfull(i,j,Conv_results%cldbot(i,j))
          end do
        end do
        used = send_data(id_ccb, temp_2d, Time, is_in=is,   &
                                js_in=js,  mask = Conv_results%cldbot > 0)
      end if

      if ( id_conv_cld_top > 0 ) then
        temp_2d = missing_value
        do j = 1,jx
          do i = 1,ix
            if ( Conv_results%cldtop(i,j) > 0 ) temp_2d(i,j) =   &
                               Input_mp%pfull(i,j,Conv_results%cldtop(i,j))
          end do
        end do
        used = send_data(id_conv_cld_top, temp_2d, Time, is_in=is, &
                                js_in=js,  mask = Conv_results%cldtop > 0)
      end if

      if ( id_cct > 0 ) then
        temp_2d = CMOR_MISSING_VALUE
        do j = 1,jx
          do i = 1,ix
            if ( Conv_results%cldtop(i,j) > 0 ) temp_2d(i,j) =    &
                              Input_mp%pfull(i,j,Conv_results%cldtop(i,j))
          end do
        end do
        used = send_data(id_cct, temp_2d, Time, is_in=is, &
                               js_in=js,  mask = Conv_results%cldtop > 0)
      end if


!---------------------------------------------------------------------
!    temperature change due to dry and moist convection:
!---------------------------------------------------------------------
      used = send_data (id_tdt_conv, Tend_mp%ttnd_conv, Time, is, js, 1)
      if (query_cmip_diag_id(ID_tntc)) then
        used = send_cmip_data_3d (ID_tntc, Tend_mp%ttnd_conv, Time, is,  &
                                        js,1, phalf=log(Input_mp%phalf))
      endif

!---------------------------------------------------------------------
!    vapor specific humidity change due to convection:
!---------------------------------------------------------------------
      used = send_data (id_qdt_conv, Tend_mp%qtnd_conv, Time, is, js, 1)
      if (query_cmip_diag_id(ID_tnhusc)) then
        used = send_cmip_data_3d (ID_tnhusc, Tend_mp%qtnd_conv, Time,  &
                                   is, js, 1, phalf=log(Input_mp%phalf))
      endif

!---------------------------------------------------------------------
!    total precipitation due to convection (both FMS and CMOR standards):
!---------------------------------------------------------------------
      used = send_data (id_prec_conv, Output_mp%precip, Time, is, js)
      used = send_data (id_prc, Output_mp%precip, Time, is, js)
      if (id_prc_g > 0)  call buffer_global_diag     &
                          (id_prc_g, Output_mp%precip(:,:), Time, is, js)

!---------------------------------------------------------------------
!    frozen precipitation (snow) due to convection:
!---------------------------------------------------------------------
      used = send_data (id_snow_conv, Output_mp%fprec, Time, is, js)
      used = send_data (id_prsnc, Output_mp%fprec, Time, is, js)
!---------------------------------------------------------------------
!    liquid precipitation (rain) due to convection:
!---------------------------------------------------------------------
      used = send_data (id_prrc, Output_mp%lprec, Time, is, js)

!---------------------------------------------------------------------
!    convective frequency (both FMS and CMOR standards).
!---------------------------------------------------------------------
      if (id_conv_freq > 0 .or. id_ci > 0) then
        ltemp = Output_mp%precip > 0. .or. Conv_results%cldtop > 0
        where (ltemp)
          freq_count = 1.
        elsewhere
          freq_count = 0.
        end where
        if (id_conv_freq > 0) &
            used = send_data (id_conv_freq, freq_count, Time, is, js )
        if (id_ci > 0)     &
            used = send_data (id_ci, freq_count, Time, is, js )
      endif

!---------------------------------------------------------------------
!    surface wind gustiness due to convection:
!---------------------------------------------------------------------
      used = send_data (id_gust_conv, Output_mp%gust_cv, Time, is, js)

!---------------------------------------------------------------------
!    water vapor path tendency due to convection:
!---------------------------------------------------------------------
      if (id_q_conv_col > 0)   &
           call column_diag (id_q_conv_col, is, js, Time, &
                              Tend_mp%qtnd_conv, 1.0, Input_mp%pmass)
  
!---------------------------------------------------------------------
!    dry static energy tendency due to dry and moist convection:
!---------------------------------------------------------------------
      if (id_t_conv_col > 0)   &
            call column_diag (id_t_conv_col, is, js, Time, &
                               Tend_mp%ttnd_conv, CP_AIR, Input_mp%pmass)
   
!---------------------------------------------------------------------
!    define the total prognostic cloud liquid, ice, drop number, 
!    ice number and area tendencies due to convection.
!---------------------------------------------------------------------
      if (doing_prog_clouds) then
        Tend_mp%qldt_conv = Output_mp%rdt(:,:,:,nql) -  &
                                          Output_mp%rdt_init(:,:,:,nql)
        Tend_mp%qidt_conv = Output_mp%rdt(:,:,:,nqi) -   &
                                         Output_mp%rdt_init(:,:,:,nqi)
        if (do_liq_num) Tend_mp%qndt_conv =    &
                           Output_mp%rdt(:,:,:,nqn) -   &
                                       Output_mp%rdt_init(:,:,:,nqn)
        if (do_ice_num) Tend_mp%qnidt_conv =    &
                              Output_mp%rdt(:,:,:,nqni) -   &
                                      Output_mp%rdt_init(:,:,:,nqni)
        Tend_mp%qadt_conv = Output_mp%rdt(:,:,:,nqa) -    &
                                       Output_mp%rdt_init(:,:,:,nqa)

!---------------------------------------------------------------------
!    output diagnostics for cloud liquid tendency and liquid water path 
!    tendency due to convection.
!---------------------------------------------------------------------
        if (id_qldt_conv > 0 .or. id_ql_conv_col > 0) then
          used = send_data (id_qldt_conv, Tend_mp%qldt_conv,    &
                                                       Time, is, js, 1)
          if (id_ql_conv_col > 0)    &
                call column_diag (id_ql_conv_col, is, js, Time,   &
                                  Tend_mp%qldt_conv, 1.0, Input_mp%pmass)

        endif

!---------------------------------------------------------------------
!    output diagnostics for cloud drop number tendency and cloud drop 
!    number path tendency due to convection.
!---------------------------------------------------------------------
        if (do_liq_num) then
          if (id_qndt_conv > 0 .or. id_qn_conv_col > 0) then
            used = send_data (id_qndt_conv, Tend_mp%qndt_conv,    &
                                                       Time, is, js, 1)
            if (id_qn_conv_col > 0)     &
                call column_diag (id_qn_conv_col, is, js, Time,   &
                                  Tend_mp%qndt_conv, 1.0, Input_mp%pmass)
          endif
        endif

!---------------------------------------------------------------------
!    output diagnostics for cloud ice tendency and cloud ice water path 
!    tendency due to convection.
!---------------------------------------------------------------------
        if (id_qidt_conv > 0 .or. id_qi_conv_col > 0) then
          used = send_data (id_qidt_conv, Tend_mp%qidt_conv, Time,   &
                                                              is, js, 1)
          if (id_qi_conv_col > 0)    &
                call column_diag (id_qi_conv_col, is, js, Time,    &
                                  Tend_mp%qidt_conv, 1.0, Input_mp%pmass)
        endif        


!---------------------------------------------------------------------
!    output diagnostics for cloud ice number tendency and cloud ice number
!    path tendency due to convection.
!---------------------------------------------------------------------
        if (do_ice_num) then
          if (id_qnidt_conv > 0 .or. id_qni_conv_col > 0) then
            used = send_data (id_qnidt_conv, Tend_mp%qnidt_conv,    &
                                                         Time, is, js, 1)
            if (id_qni_conv_col > 0)   &
                 call column_diag (id_qni_conv_col, is, js, Time,   &
                                Tend_mp%qnidt_conv, 1.0, Input_mp%pmass)
          endif
        endif

!---------------------------------------------------------------------
!    output diagnostics for cloud area tendency and column integrated 
!    cloud mass tendency due to convection.
!---------------------------------------------------------------------
        if (id_qadt_conv > 0 .or.  id_qa_conv_col > 0 ) then
          used = send_data (id_qadt_conv, Tend_mp%qadt_conv,    &
                                                        Time, is, js, 1)
          if (id_qa_conv_col > 0)     &
               call column_diag (id_qa_conv_col, is, js, Time,    &
                                Tend_mp%qadt_conv, 1.0, Input_mp%pmass)
        endif
      endif !(doing_prog_clouds)
         
!---------------------------------------------------------------------
!    compute the column integrated enthalpy and total water tendencies 
!    due to convection parameterizations, if those diagnostics are desired.
!---------------------------------------------------------------------
      if (id_enth_conv_col > 0 .or. id_wat_conv_col > 0) then
        temp_3d1 = Output_mp%rdt(:,:,:,nql) - Output_mp%rdt_init(:,:,:,nql)
        temp_3d2 = Output_mp%rdt(:,:,:,nqi) - Output_mp%rdt_init(:,:,:,nqi)

        if (id_enth_conv_col > 0) then
          temp_2d = -HLV*Output_mp%precip -HLF*Output_mp%fprec
          call column_diag    &
             (id_enth_conv_col, is, js, Time, Tend_mp%ttnd_conv, CP_AIR, &
                  temp_3d1, -HLV, temp_3d2, -HLS, Input_mp%pmass, temp_2d)
        endif

        if (id_wat_conv_col > 0) then
          temp_2d = Output_mp%precip
          call column_diag   &
             (id_wat_conv_col, is, js, Time, Tend_mp%qtnd_conv, 1.0,   &
                   temp_3d1, 1.0, temp_3d2, 1.0, Input_mp%pmass, temp_2d)
        endif
      endif


!---------------------------------------------------------------------
!    compute the tracer tendencies due to convection for any tracers that
!    are to be transported by any convective parameterization.
!---------------------------------------------------------------------
      do n=1,size(Output_mp%rdt,4)
        if ( tracers_in_uw(n))    then
 
!---------------------------------------------------------------------
!    output diagnostics for tracer tendency and column integrated 
!    tracer tendency due to convection.
!---------------------------------------------------------------------
          if (id_tracerdt_conv(n) > 0 .or.    &
                                 id_tracerdt_conv_col(n) > 0) then
            temp_3d1 = Output_mp%rdt(:,:,:,n) - Output_mp%rdt_init(:,:,:,n)
            used = send_data (id_tracerdt_conv(n), temp_3d1,    &
                                                       Time, is, js, 1 )

            if (id_tracerdt_conv_col(n) > 0) &
              call column_diag    &
                  (id_tracerdt_conv_col(n), is, js, Time, temp_3d1,   &
                                                      1.0, Input_mp%pmass)
          endif        
        endif
      end do

!------------------------------------------------------------------


end subroutine convective_diagnostics



!######################################################################

subroutine define_convective_area (C2ls_mp, Moist_clouds_block, Input_mp)

!----------------------------------------------------------------------
!    subroutine define_convective_area computes the grid box area taken up
!    by convective clouds and thus unavailable to the large-scale cloud
!    scheme, and the ratio of gridbox relative humidity to that in the
!    convective cloud environment.
!----------------------------------------------------------------------

!---------------------------------------------------------------------
type(mp_conv2ls_type),              intent(inout) :: C2ls_mp
type(clouds_from_moist_block_type), intent(in)    :: Moist_clouds_block
type(mp_input_type),                intent(in)    :: Input_mp

!-----------------------------------------------------------------------
!    C2ls_mp    derived type used to transfer data from convection_driver
!               to lscloud_driver via moist_processes.
!    Moist_clouds_block 
!               derived type used to transfer cloud data between 
!               atmos_model and convection_driver via physics_driver and
!               moist_processes
!    Input_mp   derived type used to transfer needed input data between
!               moist_processes and convection_driver
!
!----------------------------------------------------------------------

      real, dimension(size(Input_mp%t,1),size(Input_mp%t,2),   &
                                         size(Input_mp%t,3)) ::  &
                                             conv_area_input,  &
                                             rh_wtd_conv_area

!---------------------------------------------------------------------
!      conv_area_input  area taken up by convective clouds, summed over 
!                       the active convective schemes which predict clouds.
!                       this is the area unavailable for large-scale
!                       clouds.
!      rh_wtd_conv_area sum of convective area times relative humidity,
!                       summed over active convective cloud schemes 
!                       (relative-humidity-weighted convective area). 
!                       uw, updraft level each contribute  CF*1.0 (their 
!                       cloud areas are assumed saturated)
!---------------------------------------------------------------------

!-----------------------------------------------------------------------
!    define the total convective area and relative-humidity-weighted 
!    convective area for the current convective implementation.
!    if no convective scheme which produces convective clouds is active,
!    set these fields to 0.0.
!-----------------------------------------------------------------------
      if (do_uw_conv) then
        conv_area_input =   &
                    Moist_clouds_block%cloud_data(i_shallow)%cloud_area
        rh_wtd_conv_area = &
                    Moist_clouds_block%cloud_data(i_shallow)%cloud_area
      else
        conv_area_input = 0.
        rh_wtd_conv_area = 0.
      endif

!-----------------------------------------------------------------------
!    call compute_convective_area to compute the grid box area taken up 
!    by convective clouds, and the ratio of gridbox relative humidity to
!    that in the cloud environment.  
!-----------------------------------------------------------------------
      if (.not. do_lsc) then
        call compute_convective_area     &
            (Input_mp%tin, Input_mp%pfull, Input_mp%qin, conv_area_input, &
             rh_wtd_conv_area, 1.0, C2ls_mp%convective_humidity_ratio, &
                                         C2ls_mp%convective_humidity_area)
      endif

!---------------------------------------------------------------------


end subroutine define_convective_area   



!#######################################################################

subroutine compute_convective_area     &
                 (t, pfull, q, conv_area_input, rh_wtd_conv_area,   &
                           max_cnv_frac, humidity_ratio, convective_area)

!-------------------------------------------------------------------------
!    subroutine compute_convective_area defines the grid box area affected
!    by the convective clouds (convective_area) and the ratio of the 
!    grid-box relative humidity to the humidity in the environment of the 
!    convective clouds (humidity_ratio). 
!-------------------------------------------------------------------------

real, dimension(:,:,:),  intent(in)   :: t, pfull, q, &
                                         conv_area_input, rh_wtd_conv_area
real,                    intent(in)   :: max_cnv_frac
real, dimension(:,:,:),  intent(out)  :: humidity_ratio, convective_area

!------------------------------------------------------------------------

!----------------------------------------------------------------------
!         t        temperature            [ K ]
!         pfull    pressure on full levels [ Pa ]
!         q        specific humidity [ kg h2o / kg moist air ]
!         conv_area_input
!                  convective cloud area as determined by active convective
!                  schemes
!         rh_wtd_conv_area
!                  sum of products of convective area * rh over all active
!                  convective cloud schemes
!         max_cnv_frac
!                 largest area in a gridbox which may be taken up by 
!                 convective clouds
!         humidity_ratio
!                 ratio of the grid-box relative humidity to the humidity 
!                 in the environment of the convective clouds  
!         convective_area 
!                 grid box area affected by the convective clouds 
!         
!----------------------------------------------------------------------

!----------------------------------------------------------------------
      real, dimension(size(t,1), size(t,2),   &
                             size(t,3)) :: qs, qrf, env_fraction, env_qv
      integer :: i, j ,k
      integer :: ix, jx, kx

!-----------------------------------------------------------------------
!      qs              saturation specific humidity
!      qrf             model specific humidity, forced to be realizable
!      env_fraction    portion of gridbox free of convective cloud 
!                      influence
!      env_qv          temporary variable used in calculation of
!                      humidity_ratio
!      i, j, k         do loop indices
!      ix, jx, kx      dimensions of physics window
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!    define array dimensions.
!-----------------------------------------------------------------------
      ix = size(t,1)
      jx = size(t,2)
      kx = size(t,3)

!-----------------------------------------------------------------------
!    define the grid box area whose humidity is affected by the 
!    convective clouds (convective_area) initially as that obtained from 
!    the active cloud schemes (conv_area_input) and passed into this 
!    subroutine. here it may be limited by max_cnv_frac, which may be set
!    either arbitrarily or by an nml variable. the convective environment 
!    fraction (env_fraction) is  defined as the remainder of the box.
!-----------------------------------------------------------------------
      do k=1, kx
        do j=1,jx   
          do i=1,ix  
            convective_area(i,j,k) = min (conv_area_input(i,j,k), &
                                                             max_cnv_frac)
            env_fraction(i,j,k) = 1.0 - convective_area(i,j,k)
          end do
        end do
      end do

!-----------------------------------------------------------------------
!    calculate the ratio of the gridbox relative humidity to that
!    in the cloud environment (humidity_ratio). to do so, define a 
!    realizable grid box specific humidity (qrf) and the saturation 
!    specific humidity (qs).
!------------------------------------------------------------------
      qrf = MAX(q, 0.0)
      call compute_qs (t, pfull, qs)

!----------------------------------------------------------------------
!    given the gridbox specific humidity (qrf) and the convective area
!    specific humidity (based on qs), the environmental specific humidity
!    must be given by
!
!      q(ENVIRONMENT) =  (q(GRIDBOX) -q(ConvectiveArea)*ConvectiveArea)/ &
!                                                     (1 - ConvectiveArea).
!
!    the convective cloud area is assumed saturated for the uw clouds
!      
!    variable env_qv is defined as the numerator in the above expression.
!    where the ConvectiveArea has been passed in as rh_wtd_conv_area, 
!    taking account of the different treatment of qs in the cloud area by 
!    uw clouds.
!-------------------------------------------------------------------
      env_qv = qrf - qs*rh_wtd_conv_area
      do k=1,kx
        do j=1,jx   
          do i=1,ix  

!---------------------------------------------------------------------
!    one can define the ratio of the grid-box relative humidity to the 
!    humidity in the environment of the convective clouds only if the 
!    grid box contains vapor (qrf > 0.0) and there is some vapor
!    outside of the convective clouds (env_qv > 0.).
!----------------------------------------------------------------------
            if (qrf(i,j,k) /= 0.0 .and. env_qv(i,j,k) > 0.0) then
 
!--------------------------------------------------------------------
!    there must also be some grid box area that does not contain convective
!    clouds (env_fraction > 0.).
!--------------------------------------------------------------------  
              if (env_fraction(i,j,k) > 0.0) then
                humidity_ratio(i,j,k) =    &
                   MAX (qrf(i,j,k)*env_fraction(i,j,k)/env_qv(i,j,k), 1.0)
 
!---------------------------------------------------------------------
!    if the grid box is filled with convective clouds, set humidity ratio 
!    to a flag value. this will not happen if max_cnv_frac is < 1.0.
!----------------------------------------------------------------------
              else
                humidity_ratio(i,j,k) = -10.0
              endif

!--------------------------------------------------------------------
!    if there either is no vapor in the gridbox or all the vapor has been 
!    taken up by the convective clouds, set the humidity_ratio to 1.0.
!---------------------------------------------------------------------
            else
              humidity_ratio(i,j,k) = 1.0
            endif
          end do
        end do
      end do

!------------------------------------------------------------------------


end subroutine compute_convective_area



!######################################################################

subroutine define_inputs_for_cosp (Removal_mp)

!---------------------------------------------------------------------
!    subroutine define_inputs_for_cosp  converts the precip fluxes from
!    =mid-layer values to interface values.
!---------------------------------------------------------------------

type(mp_removal_type),  intent(inout) :: Removal_mp

!---------------------------------------------------------------------
!    Removal_mp derived type used to transfer precipitation and tracer
!               removal fields between convection_driver and 
!               moist_processes
!---------------------------------------------------------------------

      integer :: k

!--------------------------------------------------------------------
!    k     do loop index
!--------------------------------------------------------------------

!---------------------------------------------------------------------
!    define precip fluxes (index 1 is model lid)
!---------------------------------------------------------------------
      do k=2, size(Removal_mp%ice_precflxh,3)
        Removal_mp%ice_precflxh(:,:,k) =                  &
                                      Removal_mp%ice_precflxh(:,:,k-1) +  &
                                      Removal_mp%ice_precflx(:,:,k-1)
        Removal_mp%liq_precflxh(:,:,k) =        &
                                      Removal_mp%liq_precflxh(:,:,k-1) +  &
                                      Removal_mp%liq_precflx(:,:,k-1)
      end do

!--------------------------------------------------------------------
!    adjust precip fluxes to remove any negative values that were produced.
!    precip contribution is determined as the negative of the total 
!    moisture tendency, so at top of clouds a positive moisture tendency 
!    sometimes results in a negative precipitation contribution. 
!----------------------------------------------------------------------
      call prevent_neg_precip_fluxes (Removal_mp%ice_precflxh)
      call prevent_neg_precip_fluxes (Removal_mp%liq_precflxh)
!-----------------------------------------------------------------------


end subroutine define_inputs_for_cosp 



!#######################################################################

subroutine prevent_neg_precip_fluxes (fluxh)

!----------------------------------------------------------------------
!    subroutine prevent_neg_precip_fluxes checks for negative precip
!    fluxes (implying precip is moving upwards in the atmosphere) at each 
!    level, and if encountered the flux is eliminated.
!----------------------------------------------------------------------

real, dimension(:,:,:), intent(inout) :: fluxh

!---------------------------------------------------------------------
!  fluxh  precip flux at model half-level
!---------------------------------------------------------------------

      real, dimension(size(fluxh,1), size(fluxh,2)) :: sumneg
      integer :: i,j,k

!---------------------------------------------------------------------
!  sumneg      the accumulated vertical sum of unbalanced negative 
!              precip fluxes in the column (beginning at the top)
!  i, j, k     do loop indices
!---------------------------------------------------------------------

!-----------------------------------------------------------------------
!    move down each column looking for negative precip fluxes at each
!    level. if found, the negative flux is eliminated at the level, and 
!    positive fluxes lower down will be reduced until the negative flux 
!    is balanced.
!-----------------------------------------------------------------------
      sumneg(:,:) = 0.
      do k=2, size(fluxh,3)
        do j=1,size(fluxh,2)
          do i=1,size(fluxh,1)
            if (fluxh(i,j,k) > 0.0) then
              if (fluxh(i,j,k) > ABS(sumneg(i,j))) then
                fluxh(i,j,k) = fluxh(i,j,k) + sumneg(i,j)
                sumneg(i,j) = 0.
              else
                sumneg(i,j) = sumneg(i,j) + fluxh(i,j,k)
                fluxh(i,j,k) = 0.
              endif
            else
              sumneg(i,j) = sumneg(i,j) + fluxh(i,j,k)
              fluxh(i,j,k) = 0.
            endif
          end do
        end do
      end do

!----------------------------------------------------------------------


end subroutine prevent_neg_precip_fluxes



!#######################################################################

subroutine convection_driver_dealloc (Conv_results)

!--------------------------------------------------------------------
!    subroutine convection_driver_dealloc deallocates the components of
!    the conv_results_type variable Conv_results.
!--------------------------------------------------------------------

type(conv_results_type), intent(inout)   :: Conv_results

!---------------------------------------------------------------------
!    Conv_results
!               conv_results_type variable containing local variables 
!               used in multiple convective parameterizations and for
!               diagnostic output
!---------------------------------------------------------------------

!------------------------------------------------------------------------
!    deallocate the components of the conv_results_type variable
!    Conv_results.
!------------------------------------------------------------------------
      deallocate (Conv_results%uw_mflux)  ! cmf

      deallocate(Conv_results%available_cf_for_uw)
      deallocate(Conv_results%conv_calc_completed)

      deallocate (Conv_results%cldtop) 
      deallocate (Conv_results%cldbot) 
      deallocate (Conv_results%prod_no) 

!----------------------------------------------------------------------


end subroutine convection_driver_dealloc


!#######################################################################

!*******************************************************************
!
!                     PRIVATE, UW-RELATED SUBROUTINES
!
!*******************************************************************

subroutine uw_conv_driver    &
              (is, ie, js, je, Input_mp, Aerosol, Phys_mp_exch, &
              Output_mp, Tend_mp, Conv_results, Removal_mp, Cld_props)

!----------------------------------------------------------------------
!   subroutine uw_conv_driver prepares for and executes the uw convection
!   parameterization, and then processes its output appropriately. both
!   parts of uw_conv_driver_part are executed when called from this 
!   subroutine.
!----------------------------------------------------------------------

integer,                      intent(in)     :: is, ie, js, je
type(mp_input_type),          intent(inout)  :: Input_mp
type(aerosol_type),           intent(in)     :: Aerosol
type(phys_mp_exch_type),      intent(inout)  :: Phys_mp_exch
type(mp_output_type),         intent(inout)  :: Output_mp
type(mp_tendency_type),       intent(inout)  :: Tend_mp
type(conv_results_type),      intent(inout)  :: Conv_results
type(mp_removal_type),        intent(inout)  :: Removal_mp
type(cloud_scheme_data_type), intent(inout)  :: Cld_props

!---------------------------------------------------------------------
!    is,js      starting i and j indices for window
!    ie,je      ending i and j indices for window
!    Input_mp   derived type used to transfer needed input data between
!               moist_processes and convection_driver
!    Aerosol    derived type containing model aerosol fields to be input
!               to model convective schemes
!    Phys_mp_exch
!               derived type used to transfer data between physics_driver
!               and convection_driver via moist_processes
!    Output_mp  derived type used to transfer output fields between
!               convection_driver and moist_processes
!    Tend_mp    derived type used to transfer calculated tendency data
!               between convection_driver and moist_processes
!    Conv_results
!               conv_results_type variable containing variables 
!               used in multiple convective parameterizations and for
!               diagnostic output
!    Removal_mp derived type used to transfer precipitation and tracer
!               removal fields between convection_driver and 
!               moist_processes
!    Cld_props  derived type used to transfer uw cloud data between 
!               atmos_model and convection_driver via physics_driver and
!               moist_processes
!-----------------------------------------------------------------------

      type(conv_tendency_type) :: Uw_tend
      type(conv_output_type)   :: Output_uw

!------------------------------------------------------------------------
!    Uw_tend      conv_tendency_type variable containing tendency
!                 output from uw convection
!    Output_uw    conv_output_type variable containing output
!                 fields from uw convection
!------------------------------------------------------------------------

!-----------------------------------------------------------------------
!    activate the uw clock.
!-----------------------------------------------------------------------
      call mpp_clock_begin (uw_clock)

!-----------------------------------------------------------------------
!    call uw_conv_driver_part, executing both parts of the driver.
!-----------------------------------------------------------------------
      call uw_conv_driver_part  &
              (is, ie, js, je, Input_mp, Aerosol, Phys_mp_exch, &
              Output_mp, Tend_mp, Conv_results, Removal_mp, Cld_props, &
                                    Uw_tend, Output_uw, .true., .true.)

!----------------------------------------------------------------------
!    turn off the uw clock.
!----------------------------------------------------------------------
      call mpp_clock_end (uw_clock)

!------------------------------------------------------------------------


end subroutine uw_conv_driver



!########################################################################

subroutine uw_conv_driver_part    &
              (is, ie, js, je, Input_mp, Aerosol, Phys_mp_exch, &
              Output_mp, Tend_mp, Conv_results, Removal_mp, Cld_props, &
                           Uw_tend, Output_uw, do_segment1, do_segment2)

!-----------------------------------------------------------------------
integer,                        intent(in)    :: is, ie, js, je
type(mp_input_type),            intent(inout) :: Input_mp
type(aerosol_type),             intent(in)    :: Aerosol
type(phys_mp_exch_type),        intent(inout) :: Phys_mp_exch
type(mp_output_type),           intent(inout) :: Output_mp
type(mp_tendency_type),         intent(inout) :: Tend_mp
type(conv_results_type),        intent(inout) :: Conv_results
type(mp_removal_type),          intent(inout) :: Removal_mp
type(cloud_scheme_data_type),   intent(inout) :: Cld_props
type(conv_tendency_type),       intent(inout) :: Uw_tend
type(conv_output_type),         intent(inout) :: Output_uw
logical,                        intent(in)    :: do_segment1, do_segment2
                                           
!-----------------------------------------------------------------------
!    is,js      starting i and j indices for window
!    ie,je      ending i and j indices for window
!    Input_mp   derived type used to transfer needed input data between
!               moist_processes and convection_driver
!    Aerosol    derived type containing model aerosol fields to be input
!               to model convective schemes
!    Phys_mp_exch
!               derived type used to transfer data between physics_driver
!               and convection_driver via moist_processes
!    Output_mp  derived type used to transfer output fields between
!               convection_driver and moist_processes
!    Tend_mp    derived type used to transfer calculated tendency data
!               between convection_driver and moist_processes
!    Conv_results
!               conv_results_type variable containing variables 
!               used in multiple convective parameterizations and for
!               diagnostic output
!    Removal_mp derived type used to transfer precipitation and tracer
!               removal fields between convection_driver and 
!               moist_processes
!    Cld_props  derived type used to transfer uw cloud data between 
!               atmos_model and convection_driver via physics_driver and
!               moist_processes
!    Uw_tend    conv_tendency_type variable containing tendency
!               output from uw convection
!    Output_uw  conv_output_type variable containing output
!               fields from uw convection
!    do_segment1
!               logical indicating if first part of uw_conv_driver is to
!               be executed
!    do_segment2
!               logical indicating if second part of uw_conv_driver is to
!               be executed
!-----------------------------------------------------------------------

      real, dimension(size(Input_mp%t,1), size(Input_mp%t,2),      &
                                         size(Input_mp%t,3))  :: targ, qarg
      real, dimension(size(Input_mp%t,1), size(Input_mp%t,2),      &
                         size(Input_mp%t,3), num_uw_tracers)  :: trcr
      real, dimension(size(Input_mp%t,1), size(Input_mp%t,2),      &
                         size(Input_mp%t,3),       &
                                   num_prog_tracers)          :: tracerarg
      integer :: ix, jx, kx
      logical :: used
      integer :: nt                           
      integer :: n                           
      integer :: nn

!---------------------------------------------------------------------
!   targ           temperature field passed to uw_conv
!   qarg           specific humidity field passed to uw_conv
!   trcr           set of tracers actually transported by uw_conv
!   tracerarg      tracer fields passed to uw_conv
!   ix, jx, kx     physics window dimensions
!   nt             number of prognostic tracers
!   n              do loop index
!   nn             counter
!----------------------------------------------------------------------

!-----------------------------------------------------------------------
!    define array dimensions and number of prognostic tracers.
!-----------------------------------------------------------------------
      ix = size(Input_mp%t,1) 
      jx = size(Input_mp%t,2) 
      kx = size(Input_mp%t,3) 
      nt = size(Output_mp%rdt,4)

!-----------------------------------------------------------------------
!    this first part is executed when do_segment1 = .true.. to reproduce
!    the warsaw results it is necessary to execute the first part
!-----------------------------------------------------------------------
      if (do_segment1) then

!----------------------------------------------------------------------
!    call uw_alloc to allocate the needed components of the local derived
!    type varaibles resident in this module.
!----------------------------------------------------------------------
        call uw_alloc (ix, jx, kx, Uw_tend, Output_uw)

        targ = Input_mp%tin_orig
        qarg = Input_mp%qin_orig
        tracerarg = Input_mp%tracer_orig

!----------------------------------------------------------------------
!    if any tracers are to be transported by UW convection, check each
!    active tracer to find those to be transported and fill the 
!    trcr array with these fields.
!---------------------------------------------------------------------
        nn = 1
        do n=1, num_prog_tracers
          if (tracers_in_uw(n)) then
            trcr(:,:,:,nn) = tracerarg(:,:,:,n)
            nn = nn + 1
          endif
        end do

!-------------------------------------------------------------------------
!     call uw_conv to calculate the effects of shallow convection.
!-------------------------------------------------------------------------
        call uw_conv (is, js, Time, targ, qarg, Input_mp%uin,   &
                      Input_mp%vin, Input_mp%pfull, Input_mp%phalf,&
                      Input_mp% zfull, Input_mp%zhalf, tracerarg,   &
                      Input_mp%omega, dt, Input_mp%pblht, &
                      Input_mp%ustar, Input_mp%bstar, Input_mp%qstar,  &
                      Input_mp%land, Input_mp%coldT, Aerosol,   &
                      Input_mp%lat, Input_mp%lon, Input_mp%cush, &
                      doing_prog_clouds,  &
                      Conv_results%conv_calc_completed,   &
                      Conv_results%available_cf_for_uw, Uw_tend%ttnd,    &
                      Uw_tend%qtnd, Uw_tend%qltnd, Uw_tend%qitnd, &
                      Uw_tend%qatnd, Uw_tend%qntnd,      &
                      Uw_tend%utnd, Uw_tend%vtnd, Uw_tend%rain, &
                      Uw_tend%snow, Conv_results%uw_mflux,   &
                      Removal_mp%liq_precflx, &
                      Removal_mp%ice_precflx, Cld_props%liquid_amt,   &
                      Cld_props%ice_amt, Cld_props%cloud_area,   &
                      Cld_props%droplet_number, trcr, Uw_tend%qtr,   &
                      Removal_mp%uw_wetdep, Input_mp%cbmf,    &
                      Phys_mp_exch%cgust, Phys_mp_exch%tten_gw, Phys_mp_exch%cqa_gw)

!-------------------------------------------------------------------------
!    call detr_ice_num to calculate the ice number tendency due to 
!    detrainment, which is proportional to the ice mass.
!-------------------------------------------------------------------------
        if (do_ice_num .and. detrain_ice_num) then
          CALL detr_ice_num (targ, Uw_tend%qitnd(:,:,:),    &
                                                 Uw_tend%qnitnd(:,:,:) )
        end if
      endif ! (do_segment1)

!----------------------------------------------------------------------
!    the second segment is executed whe do_segment2 is .true..
!----------------------------------------------------------------------
      if (do_segment2) then

!---------------------------------------------------------------------
!    if desired, call define_and_apply_scale to compute any adjustment 
!    needed in order to preserve realizability for the water variables.
!---------------------------------------------------------------------
        if (do_limit_uw) then
          call define_and_apply_scale   &
              (Input_mp, Uw_tend, Output_uw, .true., Uw_tend%qtr)
        else  
          Output_uw%scale = 1.0
        endif 

!-----------------------------------------------------------------------
!    call finalize_uw_outputs to define output fields, update input
!    fields as needed, and output uw-related diagnostics.
!-----------------------------------------------------------------------
        call finalize_uw_outputs (is, js, Input_mp, Uw_tend, Output_uw, &
                                                              Output_mp)

!-----------------------------------------------------------------------
!    call update_outputs to update the arrays which will return the
!    convective tendencies to moist_processes.
!-----------------------------------------------------------------------
        call update_outputs (Uw_tend, Output_mp, Tend_mp)

!----------------------------------------------------------------------
!    call uw_dealloc to deallocate the components of the derived type
!    variables Uw_tend and Output_uw.
!----------------------------------------------------------------------
        call uw_dealloc (Uw_tend, Output_uw)

!----------------------------------------------------------------------
!    end of segment2.
!----------------------------------------------------------------------
      endif ! (do_segment2)

!------------------------------------------------------------------------


end subroutine uw_conv_driver_part



!#####################################################################

subroutine uw_alloc (ix, jx, kx, Uw_tend, Output_uw)

!----------------------------------------------------------------------
!    subroutine uw_alloc allocates the needed components of the
!    conv_tendency_type and conv_output_type variables.
!----------------------------------------------------------------------
integer,                  intent(in)    :: ix, jx, kx
type(conv_tendency_type), intent(inout) :: Uw_tend
type(conv_output_type),   intent(inout) :: Output_uw

!------------------------------------------------------------------------
!    ix, jx, kx      physics window dimensions
!    Uw_tend         conv_tendency_type variable containing tendency
!                    output from uw convection
!    Output_uw       conv_output_type variable containing output
!                    fields from uw convection
!------------------------------------------------------------------------

!----------------------------------------------------------------------
!    allocate and initialize the needed components of Uw_tend.
!----------------------------------------------------------------------
      allocate (Uw_tend%delta_q   (ix, jx,kx))
      allocate (Uw_tend%rain      (ix, jx))
      allocate (Uw_tend%snow      (ix, jx))
      allocate (Uw_tend%ttnd      (ix, jx, kx))
      allocate (Uw_tend%qtnd      (ix, jx, kx))
      allocate (Uw_tend%utnd      (ix, jx, kx))
      allocate (Uw_tend%vtnd      (ix, jx, kx))
      allocate (Uw_tend%qltnd     (ix, jx, kx))
      allocate (Uw_tend%qitnd     (ix, jx, kx))
      allocate (Uw_tend%qatnd     (ix, jx, kx))
      allocate (Uw_tend%qntnd     (ix, jx, kx))
      allocate (Uw_tend%qnitnd    (ix, jx, kx))
      allocate (Uw_tend%qtr       (ix, jx, kx, num_uw_tracers))

      Uw_tend%delta_q  = 0.
      Uw_tend%rain     = 0.
      Uw_tend%snow     = 0.
      Uw_tend%ttnd     = 0.
      Uw_tend%qtnd     = 0.
      Uw_tend%utnd     = 0.
      Uw_tend%vtnd     = 0.
      Uw_tend%qltnd    = 0.
      Uw_tend%qitnd    = 0.
      Uw_tend%qatnd    = 0.
      Uw_tend%qntnd    = 0.
      Uw_tend%qnitnd   = 0.
      Uw_tend%qtr      = 0.

!----------------------------------------------------------------------
!    allocate and initialize the needed components of Output_uw.
!----------------------------------------------------------------------
      allocate (Output_uw%liquid_precip  (ix, jx, kx))
      allocate (Output_uw%frozen_precip  (ix, jx, kx))
      allocate (Output_uw%total_precip   (ix, jx))
      allocate (Output_uw%scale          (ix, jx))
      allocate (Output_uw%scale_REV      (ix, jx))

      Output_uw%liquid_precip  = 0.
      Output_uw%frozen_precip  = 0.
      Output_uw%total_precip   = 0.
      Output_uw%scale          = 1.0
      Output_uw%scale_REV      = 1.0

!----------------------------------------------------------------------


end subroutine uw_alloc



!#######################################################################

subroutine finalize_uw_outputs (is, js, Input_mp, Uw_tend, Output_uw,  &
                                                             Output_mp)

!---------------------------------------------------------------------
!    subroutine finalize_uw_outputs finishes processing the uw convection
!    results, updates the needed output variables, and produces the
!    diagnostics related to uw convection.
!---------------------------------------------------------------------

integer,                  intent(in)    :: is, js
type(mp_input_type),      intent(inout) :: Input_mp
type(conv_tendency_type), intent(inout) :: Uw_tend
type(conv_output_type),   intent(inout) :: Output_uw
type(mp_output_type),     intent(inout) :: Output_mp

!---------------------------------------------------------------------
!    is,js      starting i and j indices for window
!    Input_mp   derived type used to transfer needed input data between
!               moist_processes and convection_driver
!    Uw_tend    conv_tendency_type variable containing tendency
!               output from uw convection
!    Output_uw  conv_output_type variable containing output
!               fields from uw convection
!    Output_mp  derived type used to transfer output fields between
!               convection_driver and moist_processes
!-----------------------------------------------------------------------

      logical  :: used
      integer  :: n
      integer  :: nn

!---------------------------------------------------------------------
!   used          logical used to indicate data has been received by
!                 diag_manager_mod
!   n             do loop index
!   nn            counter
!---------------------------------------------------------------------

!---------------------------------------------------------------------
!    call update_inputs to update the components of Input_mp that have
!    been modified by uw convection. 
!---------------------------------------------------------------------
      call update_inputs (Uw_tend, Input_mp)

!---------------------------------------------------------------------
!    call uw_diagnostics to output desired diagnostics related to the
!    uw convection scheme.
!---------------------------------------------------------------------
      call uw_diagnostics (is, js, Input_mp, Uw_tend, Output_uw) 

!-----------------------------------------------------------------------
!    if the warsaw order of calculation (inconsistent) is to be 
!    retained, save the changes to temperature and to the tracers in
!    the variables xxx_tentative, so that they may be applied at a
!    later time. 
!-----------------------------------------------------------------------
      if (reproduce_AM4) then
        Input_mp%tin_tentative = Uw_tend%ttnd*dt 

!------------------------------------------------------------------------
!    save the current tracer tendencies obtained from uw transport.
!------------------------------------------------------------------------
        if (do_limit_uw) then
          nn = 1
          do n=1, num_prog_tracers
            if (tracers_in_uw(n)) then
              Output_mp%rdt_tentative(:,:,:,n) = Uw_tend%qtr(:,:,:,nn)
              nn = nn + 1
            else
              Output_mp%rdt_tentative(:,:,:,n) = 0.
            endif
          end do
        else

!------------------------------------------------------------------------
!    update the current tracer tendencies with the contributions 
!    obtained from uw transport.
!------------------------------------------------------------------------
          nn = 1
          do n=1, num_prog_tracers
            if (tracers_in_uw(n)) then
              Output_mp%rdt(:,:,:,n) = Output_mp%rdt(:,:,:,n) +   &
                                                   Uw_tend%qtr(:,:,:,nn)
              nn = nn + 1
            endif
          end do
        endif

!-----------------------------------------------------------------------
!    if warsaw results are to be corrected, the Input_mp%tin
!    and Output_mp%rdt fields are updated with the uw tendencies at this
!    point.
!-----------------------------------------------------------------------
      else
        Input_mp%tin = Input_mp%tin + Uw_tend%ttnd*dt

!------------------------------------------------------------------------
!    update the current tracer tendencies with the contributions 
!    obtained from uw transport.
!------------------------------------------------------------------------
        nn = 1
        do n=1, num_prog_tracers
          if (tracers_in_uw(n)) then
            Output_mp%rdt(:,:,:,n) = Output_mp%rdt(:,:,:,n) +   &
                                                   Uw_tend%qtr(:,:,:,nn)
            nn = nn + 1
          endif
        end do
      endif

!---------------------------------------------------------------------


end subroutine finalize_uw_outputs



!#######################################################################

subroutine update_inputs (Uw_tend, Input_mp)

!---------------------------------------------------------------------
!    subroutine update_inputs updates the fields contained in Input_mp
!    that have been modified by the uw convection calculation.
!---------------------------------------------------------------------

type(conv_tendency_type), intent(in)    :: Uw_tend
type(mp_input_type),      intent(inout) :: Input_mp

!----------------------------------------------------------------------
!    Uw_tend    conv_tendency_type variable containing tendency
!               output from uw convection
!    Input_mp   derived type used to transfer needed input data between
!               moist_processes and convection_driver
!--------------------------------------------------------------------

!-------------------------------------------------------------------------
!    update Input_mp fields with changes from uw_conv.
!-------------------------------------------------------------------------
      Input_mp%qin = Input_mp%qin + Uw_tend%qtnd*dt
      Input_mp%uin = Input_mp%uin + Uw_tend%utnd*dt
      Input_mp%vin = Input_mp%vin + Uw_tend%vtnd*dt
      Input_mp%tracer(:,:,:,nql) = Input_mp%tracer(:,:,:,nql) +    &
                                                       Uw_tend%qltnd*dt
      Input_mp%tracer(:,:,:,nqi) = Input_mp%tracer(:,:,:,nqi) +     &
                                                       Uw_tend%qitnd*dt
      Input_mp%tracer(:,:,:,nqa) = Input_mp%tracer(:,:,:,nqa) +    &
                                                       Uw_tend%qatnd*dt
      if (do_liq_num) then
        Input_mp%tracer(:,:,:,nqn) = Input_mp%tracer(:,:,:,nqn) +   &
                                                       Uw_tend%qntnd*dt
      endif
      if (do_ice_num) then
        Input_mp%tracer(:,:,:,nqni) = Input_mp%tracer(:,:,:,nqni) +   &
                                                       Uw_tend%qnitnd*dt
      endif

!------------------------------------------------------------------


end subroutine update_inputs



!#########################################################################

subroutine uw_diagnostics (is, js, Input_mp, Uw_tend, Output_uw) 

!--------------------------------------------------------------------
!    subroutine uw_diagnostics outputs uw-related diagnostics. 
!--------------------------------------------------------------------

integer,                   intent(in) :: is, js
type(mp_input_type),       intent(in) :: Input_mp
type(conv_tendency_type),  intent(in) :: Uw_tend
type(conv_output_type),    intent(in) :: Output_uw

!-------------------------------------------------------------------
!    is,js      starting i and j indices for window
!    Input_mp   derived type used to transfer needed input data between
!               moist_processes and convection_driver
!    Uw_tend    conv_tendency_type variable containing tendency
!               output from uw convection
!    Output_uw  conv_output_type variable containing output
!               fields from uw convection
!-------------------------------------------------------------------

      real, dimension(size(Input_mp%tin,1),   &
                                        size(Input_mp%tin,2)) :: temp_2d
      logical, dimension(size(Input_mp%tin,1),     &
                                        size(Input_mp%tin,2)) :: ltemp
      logical :: used

!----------------------------------------------------------------------
!   temp_2d       temporary real array
!   ltemp         temporary logical array
!   used          logical used to indicate data has been received by
!                 diag_manager_mod
!----------------------------------------------------------------------

!-----------------------------------------------------------------------
!    output the scaling factor.
!-----------------------------------------------------------------------
      used = send_data (id_scale_uw, Output_uw%scale, Time, is, js )

!-----------------------------------------------------------------------
!    output total precip and snow from the uw convection scheme.
!-----------------------------------------------------------------------
      used = send_data (id_uw_precip, Uw_tend%rain + Uw_tend%snow,   &
                                                           Time, is, js)
      used = send_data (id_uw_snow, Uw_tend%snow, Time, is, js)

!-----------------------------------------------------------------------
!    prognostic variable tendencies from uw convection.
!-----------------------------------------------------------------------
      used = send_data (id_tdt_uw, Uw_tend%ttnd, Time, is, js, 1)
      used = send_data (id_qdt_uw, Uw_tend%qtnd, Time, is, js, 1)
      used = send_data (id_qadt_uw, Uw_tend%qatnd, Time, is, js, 1)
      used = send_data (id_qldt_uw, Uw_tend%qltnd, Time, is, js, 1)
      used = send_data (id_qidt_uw, Uw_tend%qitnd, Time, is, js, 1)
      if (do_liq_num) then
        used = send_data (id_qndt_uw, Uw_tend%qntnd, Time, is, js, 1)
      endif
      if (do_ice_num) then
        used = send_data (id_qnidt_uw, Uw_tend%qnitnd, Time, is, js, 1)
      end if

!-------------------------------------------------------------------
!    enthalpy and water column tendencies from uw.
!-------------------------------------------------------------------
      if (id_enth_uw_col > 0) then
        temp_2d = -HLV*Uw_tend%rain -HLS*Uw_tend%snow
        call column_diag (id_enth_uw_col, is, js, Time, Uw_tend%ttnd,  &
                          CP_AIR, Uw_tend%qltnd, -HLV, Uw_tend%qitnd,   &
                                             -HLS, Input_mp%pmass, temp_2d)
      endif

      if (id_wat_uw_col > 0) then
        temp_2d = Uw_tend%rain + Uw_tend%snow
        call column_diag(id_wat_uw_col, is, js, Time, Uw_tend%qtnd, 1.0,  &
                          Uw_tend%qltnd, 1.0, Uw_tend%qitnd, 1.0, &
                                                  Input_mp%pmass, temp_2d)
      endif
        
!----------------------------------------------------------------------
!    uw convection scheme frequency diagnostics.
!----------------------------------------------------------------------
      if (id_uw_freq > 0) then
        ltemp = Uw_tend%rain > 0. .or. Uw_tend%snow > 0.0
        where (ltemp) 
          temp_2d = 1.
        elsewhere
          temp_2d = 0.
        end where
        used = send_data (id_uw_freq, temp_2d, Time, is, js)
      endif

!----------------------------------------------------------------------


end subroutine uw_diagnostics 



!#######################################################################

subroutine uw_dealloc (Uw_tend, Output_uw)

!-----------------------------------------------------------------------
!    subroutine uw_dealloc deallocates the components of the derived type
!    variables Uw_tend and Output_uw.
!-----------------------------------------------------------------------

type(conv_tendency_type), intent(inout) :: Uw_tend
type(conv_output_type),   intent(inout) :: Output_uw

!------------------------------------------------------------------------
!    Uw_tend      conv_tendency_type variable containing tendency
!                 output from uw convection
!    Output_uw    conv_output_type variable containing output
!                 fields from uw convection
!------------------------------------------------------------------------
!---------------------------------------------------------------------
!    deallocate the Uw_tend components.
!---------------------------------------------------------------------
      deallocate (Uw_tend%delta_q)
      deallocate (Uw_tend%rain)
      deallocate (Uw_tend%snow)
      deallocate (Uw_tend%ttnd)
      deallocate (Uw_tend%qtnd)
      deallocate (Uw_tend%utnd)
      deallocate (Uw_tend%vtnd)
      deallocate (Uw_tend%qltnd)
      deallocate (Uw_tend%qitnd)
      deallocate (Uw_tend%qatnd)
      deallocate (Uw_tend%qntnd)
      deallocate (Uw_tend%qnitnd)
      deallocate (Uw_tend%qtr    )

!---------------------------------------------------------------------
!    deallocate the Output_uw components.
!---------------------------------------------------------------------
      deallocate (Output_uw%liquid_precip  )
      deallocate (Output_uw%frozen_precip  )
      deallocate (Output_uw%total_precip   )
      deallocate (Output_uw%scale          )
      deallocate (Output_uw%scale_REV   )

!---------------------------------------------------------------------


end subroutine uw_dealloc 



!#######################################################################




!*******************************************************************
!
!      PRIVATE SUBROUTINES USED BY MULTIPLE CONVECTION SCHEMES
!
!*******************************************************************


!#######################################################################

subroutine update_outputs (Conv_tend, Output_mp, Tend_mp)
                      
!-----------------------------------------------------------------------
!    subroutine update_outputs updates the physics tendency, convection 
!    tendency and moist_processes tendency arrays with the contributions 
!    from the current convection parameterization.
!-----------------------------------------------------------------------

type(conv_tendency_type),  intent(in)    :: Conv_tend
type(mp_output_type),      intent(inout) :: Output_mp
type(mp_tendency_type),    intent(inout) :: Tend_mp

!----------------------------------------------------------------------
!    Conv_tend  conv_tendency_type variable containing tendency
!               output from the current convective parameterization
!    Output_mp  derived type used to transfer output fields between
!               convection_driver and moist_processes
!    Tend_mp    derived type used to transfer calculated tendency data
!               between convection_driver and moist_processes
!----------------------------------------------------------------------

!-----------------------------------------------------------------------
!    update the physics tendency, convection tendency and moist_processes
!    tendency arrays with the contributions from the current convection
!    scheme. dependent on scheme, some of these fields may not be relevant
!    and so are not allocated.
!-------------------------------------------------------------------
      if (allocated ( Conv_tend%ttnd)) then 
        Output_mp%tdt = Output_mp%tdt + Conv_tend%ttnd
        Tend_mp%ttnd_conv = Tend_mp%ttnd_conv + Conv_tend%ttnd
        Tend_mp%ttnd = Tend_mp%ttnd + Conv_tend%ttnd
      endif
      if (allocated ( Conv_tend%qtnd)) then 
        Output_mp%rdt(:,:,:,1) = Output_mp%rdt(:,:,:,1) + Conv_tend%qtnd
        Tend_mp%qtnd_conv = Tend_mp%qtnd_conv + Conv_tend%qtnd
        Tend_mp%qtnd = Tend_mp%qtnd + Conv_tend%qtnd
      endif
      if (allocated ( Conv_tend%utnd)) &
        Output_mp%udt = Output_mp%udt + Conv_tend%utnd
      if (allocated ( Conv_tend%vtnd)) &
        Output_mp%vdt = Output_mp%vdt + Conv_tend%vtnd
      if (allocated ( Conv_tend%rain)) &
        Output_mp%lprec = Output_mp%lprec + Conv_tend%rain
      if (allocated ( Conv_tend%snow)) &
        Output_mp%fprec = Output_mp%fprec + Conv_tend%snow

!-----------------------------------------------------------------------
!    define the total precipitation rate (precip) for all schemes.
!    the different definition here is done to preserve order of
!    operations with the warsaw code release and avoid answer change
!    with this revised code.
!-----------------------------------------------------------------------
      Output_mp%precip = Output_mp%lprec + Output_mp%fprec     

!-----------------------------------------------------------------------
!    define tendencies for prognostic clloud fields, if that option is
!    active.
!-----------------------------------------------------------------------
      if (doing_prog_clouds) then
        if (allocated ( Conv_tend%qltnd)) &
          Output_mp%rdt(:,:,:,nql) = Output_mp%rdt(:,:,:,nql) +   &
                                                       Conv_tend%qltnd
        if (allocated ( Conv_tend%qitnd)) &
          Output_mp%rdt(:,:,:,nqi) = Output_mp%rdt(:,:,:,nqi) +  &
                                                       Conv_tend%qitnd
        if (allocated ( Conv_tend%qatnd)) &
          Output_mp%rdt(:,:,:,nqa) = Output_mp%rdt(:,:,:,nqa) +     &
                                                       Conv_tend%qatnd
        if (allocated ( Conv_tend%qntnd)) then 
          if (do_liq_num) Output_mp%rdt(:,:,:,nqn) =    &
                              Output_mp%rdt(:,:,:,nqn) + Conv_tend%qntnd
        endif
        if (allocated ( Conv_tend%qnitnd))  then 
          if (do_ice_num) Output_mp%rdt(:,:,:,nqni) =    &
                            Output_mp%rdt(:,:,:,nqni) + Conv_tend%qnitnd
        endif
      endif

!-----------------------------------------------------------------------


end subroutine update_outputs



!########################################################################

subroutine define_and_apply_scale (Input_mp, Conv_tend, Output_conv,&
                                          uw_scheme, qtr)

!-----------------------------------------------------------------------
!    subroutine define_and_apply_scale defines a factor to modify
!    predicted tendencies so that negative values of water and water
!    phases are not produced by the convection scheme, and values lower
!    than a specified minimum are not retained.
!    it is called by uw parameterizations, but in
!    slightly different ways in earlier code versions (warsaw and earlier).
!    those differences are preserved here to avoid changing answers; 
!    ultimately it is desirable to treat the functionality of this
!    subroutine in the same way whenever it is employed.
!-----------------------------------------------------------------------

type(mp_input_type),       intent(in)     :: Input_mp
type(conv_tendency_type),  intent(inout)  :: Conv_tend
type(conv_output_type),    intent(inout)  :: Output_conv
logical,                   intent(in)     :: uw_scheme
real, dimension(:,:,:,:),  intent(inout)  :: qtr

!-----------------------------------------------------------------------
!    Input_mp   derived type used to transfer needed input data between
!               moist_processes and convection_driver
!    Conv_tend  conv_tendency_type variable containing tendency
!               output from the current convective parameterization
!    Output_conv 
!               conv_output_type variable containing output
!               fields from the convection parameterization being processed
!    uw_scheme  logical indicating if the routine is to be handled as the
!               original uw convection code did
!    qtr        set of tracers being transported by the current convective
!               parameterization
!------------------------------------------------------------------------

      real, dimension(size(Input_mp%qin,1), size(Input_mp%qin,2),   &
                                           size(Input_mp%qin,3)) :: temp
      real          :: posdef, delta_posdef
      integer       :: ix, jx, kx                          
      integer       :: i, j, k, n 
      integer       :: nn                          

!------------------------------------------------------------------------
!    temp            temporary array
!    posdef          value of field that is to be kept non-negative before
!                    convection was calculated
!    delta_posdef    change in non-negative field due to convective 
!                    parameterization
!    ix, jx, kx      physics window dimensions
!    i, j, k, n      do loop indices
!    nn              counter
!------------------------------------------------------------------------

!-----------------------------------------------------------------------
!    define array dimensions.
!-----------------------------------------------------------------------
      ix = size(Input_mp%qin,1)
      jx = size(Input_mp%qin,2)
      kx = size(Input_mp%qin,3)

      if (uw_scheme) then
!------------------------------------------------------------------------
!    prevent the formation of negative liquid and ice, following the 
!    method employed in the warsaw code for the uw parameterization.
!------------------------------------------------------------------------
        temp = Input_mp%tracer(:,:,:,nql)/dt + Conv_tend%qltnd
        where (temp(:,:,:) .lt. 0.)
          Conv_tend%ttnd  = Conv_tend%ttnd  - temp*HLV/CP_AIR
          Conv_tend%qtnd  = Conv_tend%qtnd  + temp
          Conv_tend%qltnd = Conv_tend%qltnd - temp
        end where

        temp = Input_mp%tracer(:,:,:,nqi)/dt + Conv_tend%qitnd
        where (temp .lt. 0.)
          Conv_tend%ttnd  = Conv_tend%ttnd  - temp*HLS/CP_AIR
          Conv_tend%qtnd  = Conv_tend%qtnd  + temp
          Conv_tend%qitnd = Conv_tend%qitnd - temp
        end where

      end if
        
!-----------------------------------------------------------------------
!    compute a scaling factor for each grid point.  when this factor is
!    multiplied by the predicted tendencies, they will be reduced in
!    mgnitude to prevent the creation of negative water.
!-----------------------------------------------------------------------
      do k=1,kx
        do j=1,jx
          do i=1,ix

!-----------------------------------------------------------------------
!    the uw scheme used the total water substance as the positive definite
!    quantity that was to be preserved.
!-----------------------------------------------------------------------
            if (uw_scheme) then
              posdef = Input_mp%qin(i,j,k) + Input_mp%tracer(i,j,k,nql) + &
                                                Input_mp%tracer(i,j,k,nqi)
              delta_posdef  = ( Conv_tend%qtnd(i,j,k) +    &
                                      Conv_tend%qltnd(i,j,k) +   &
                                              Conv_tend%qitnd(i,j,k) )*dt
            endif

!-------------------------------------------------------------------------
!    if the positive definite quantity is being reduced on this step and
!    the value of that quantity after the timestep will be lower than
!    the minimum value specified, then the tendencies must be reduced to
!    preserve realizability. the percentage of predicted tendency that
!    will not cause negative values is given by the negative of the ratio 
!    of the initial field value to the predicted change (temp). 
!    temp will be a nonnegative number since posdef is positive.
!------------------------------------------------------------------------
            if (delta_posdef .lt.0 .and.    &
                            posdef + delta_posdef .lt. qmin   ) then
              temp(i,j,k) = max( 0.0, -(posdef - qmin)/delta_posdef )
            else
              temp(i,j,k) = 1.0
            endif
          end do
        end do
      end do

!-----------------------------------------------------------------------
!    the scaling factor for each column is defined as the minimum value 
!    of that ratio that is found within that column.  that value is applied
!    at each point in the column. this assures both that fields will not 
!    become negative, and that column conservation of enthalpy and water 
!    substance will not be affected by this adjustment.
!-----------------------------------------------------------------------
      Output_conv%scale = minval( temp, dim=3 )

!-----------------------------------------------------------------------
!    now apply the scaling factor to the water tracer, momentum, 
!    temperature, precipitation  and transported tracer tendencies 
!    returned from the convection scheme. NOte again that uw 
!    were originally treated differently, and that different treatment
!    is retained.
!    NOTE THAT THE TRANSPORTED TRACERS WERE NOT SCALED IN THE WARSAW
!    AND EARLIER CODE VERSIONS. THIS INCLUSION HERE MAY CHANGE ANSWERS.
!-----------------------------------------------------------------------
      if (uw_scheme) then
        do k=1,kx
          Conv_tend%utnd(:,:,k)  = Output_conv%scale*Conv_tend%utnd(:,:,k)
          Conv_tend%vtnd(:,:,k)  = Output_conv%scale*Conv_tend%vtnd(:,:,k)
          Conv_tend%ttnd(:,:,k)  = Output_conv%scale*Conv_tend%ttnd(:,:,k)
          Conv_tend%qtnd(:,:,k)  = Output_conv%scale*Conv_tend%qtnd(:,:,k)
          Conv_tend%qltnd(:,:,k) = Output_conv%scale*Conv_tend%qltnd(:,:,k)
          Conv_tend%qitnd(:,:,k) = Output_conv%scale*Conv_tend%qitnd(:,:,k)
          Conv_tend%qatnd(:,:,k) = Output_conv%scale*Conv_tend%qatnd(:,:,k)
        end do

        if (do_liq_num) then
          do k=1,kx
            Conv_tend%qntnd(:,:,k) = Output_conv%scale*  &
                                                  Conv_tend%qntnd(:,:,k)
          end do
        end if

        if (do_ice_num) then
          do k=1,kx
            Conv_tend%qnitnd(:,:,k) = Output_conv%scale*   &
                                                 Conv_tend%qnitnd(:,:,k)
          end do
        end if

        Conv_tend%rain(:,:) = Output_conv%scale*Conv_tend%rain(:,:)
        Conv_tend%snow(:,:) = Output_conv%scale*Conv_tend%snow(:,:)

      endif

!-------------------------------------------------------------------------


end subroutine define_and_apply_scale



!#####################################################################



                  end module convection_driver_mod
