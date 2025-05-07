module physics_driver_mod
! contact: Fei.Liu@noaa.gov - fil
!
!     Provides high level interfaces for calling the entire
!     FMS atmospheric physics package.
!
!    physics_driver_mod accesses the model's physics modules and
!    obtains tendencies and boundary fluxes due to the physical
!    processes that drive atmospheric time tendencies and supply 
!    boundary forcing to the surface models.
!
!     This version of physics_driver_mod has been designed around the implicit
!     version diffusion scheme of the GCM. It requires two routines to advance
!     the model one time step into the future. These two routines
!     correspond to the down and up sweeps of the standard tridiagonal solver.
!     Radiation, Rayleigh damping, gravity wave drag, vertical diffusion of
!     momentum and tracers, and the downward pass of vertical diffusion for
!     temperature and specific humidity are performed in the down routine.
!     The up routine finishes the vertical diffusion and computes moisture
!     related terms (convection,large-scale condensation, and precipitation).
!
! native format restart file
! FUTURE Deal with conservation of total energy?

use time_manager_mod,        only: time_type, get_time, operator (-), &
                                   time_manager_init, operator(*)
use field_manager_mod,       only: field_manager_init, MODEL_ATMOS
use tracer_manager_mod,      only: tracer_manager_init, &
                                   get_number_tracers, &
                                   get_tracer_names, &
                                   get_tracer_index, NO_TRACER
use block_control_mod,       only: block_control_type
use atmos_tracer_driver_mod, only: atmos_tracer_driver_init,    &
                                   atmos_tracer_driver_time_vary, &
                                   atmos_tracer_driver_endts, &
                                   atmos_tracer_driver,  &
                                   atmos_tracer_driver_end
use mpp_mod,                 only: input_nml_file, mpp_get_current_pelist
use mpp_domains_mod,         only: domain2D, mpp_get_ntile_count
use fms_mod,                 only: mpp_clock_id, mpp_clock_begin,   &
                                   mpp_clock_end, CLOCK_MODULE_DRIVER, &
                                   fms_init, stdlog, stdout,  &
                                   write_version_number, &
                                   error_mesg, FATAL,   &
                                   NOTE, check_nml_error, mpp_pe, &
                                   mpp_root_pe, mpp_chksum, string, mpp_npes
use fms2_io_mod,             only:  FmsNetcdfFile_t, FmsNetcdfDomainFile_t, &
                                   register_restart_field, register_axis, unlimited, &
                                   open_file, read_restart, write_restart, close_file, &
                                   register_field, write_data, get_global_io_domain_indices, &
                                   register_variable_attribute, variable_exists
use diag_manager_mod,        only: register_diag_field, send_data
!use fms,                     only: fms_mpp_pe, fms_mpp_root_pe


! shared atmospheric package modules:

use atmos_cmip_diag_mod,     only: register_cmip_diag_field_3d, &
                                   send_cmip_data_3d, &
                                   cmip_diag_id_type, &
                                   query_cmip_diag_id

!    shared radiation package modules:

use aerosol_types_mod,       only: aerosol_type, aerosol_time_vary_type

use physics_radiation_exch_mod, only: exchange_control_type, &
                                      clouds_from_moist_type, &
                                      clouds_from_moist_block_type, &
                                      cosp_from_rad_type, &
                                      cosp_from_rad_control_type, &
                                      cosp_from_rad_block_type, &
                                      radiation_flux_control_type, & 
                                      radiation_flux_block_type, & 
                                      alloc_clouds_from_moist_type, &
                                      alloc_cloud_scheme_data_type

use physics_types_mod,       only: alloc_physics_tendency_type, &
                                   physics_tendency_type, & 
                                   phys_mp_exch_type, &
                                   phys2cosp_type, precip_flux_type, &
                                   physics_tendency_block_type, &
                                   physics_type, & 
                                   physics_control_type, & 
                                   physics_input_block_type, &
                                   dealloc_physics_tendency_type

use moist_proc_utils_mod, only:    mp_removal_type
  
use aerosol_mod,             only: aerosol_init, aerosol_driver, &
                                   aerosol_time_vary, &
                                   aerosol_endts, &
                                   aerosol_dealloc, aerosol_end
!    component modules:

use cosp_driver_mod,         only: cosp_driver_init, cosp_driver, &
                                   cosp_driver_end, cosp_driver_time_vary, &
                                   cosp_driver_endts
use  moist_processes_mod,    only: moist_processes,    &
                                   moist_processes_init,  &
                                   set_cosp_precip_sources, &
                                   define_cosp_precip_fluxes, &
                                   moist_processes_time_vary, &
                                   moist_processes_endts, &
                                   moist_processes_restart, &
                                   moist_processes_end

use vert_turb_driver_mod,    only: vert_turb_driver,  &
                                   vert_turb_driver_init,  &
                                   vert_turb_driver_end

use vert_diff_driver_mod,    only: vert_diff_driver_down,  &
                                   vert_diff_driver_up,    &
                                   vert_diff_driver_init,  &
                                   vert_diff_driver_end,   &
                                   surf_diff_type
 
use damping_driver_mod,      only: damping_driver,      &
                                   damping_driver_init, &
                                   damping_driver_time_vary,  &
                                   damping_driver_endts, &
                                   damping_driver_end,  &
                                   damping_driver_restart

use monin_obukhov_mod,        only: monin_obukhov_init

#ifdef SCM
! Option to add SCM radiative tendencies from forcing to lw_tendency
! and radturbten

use scm_forc_mod,            only: use_scm_rad, add_scm_tdtlw, add_scm_tdtsw

#endif

!-----------------------------------------------------------------

implicit none
private

!---------------------------------------------------------------------
!    physics_driver_mod accesses the model's physics modules and
!    obtains tendencies and boundary fluxes due to the physical
!    processes that drive atmospheric time tendencies and supply 
!    boundary forcing to the surface models.
!---------------------------------------------------------------------

!---------------------------------------------------------------------
!----------- version number for this module -------------------

character(len=128) :: version = '$Id$'
character(len=128) :: tagname = '$Name$'


!---------------------------------------------------------------------
!-------  interfaces --------

public  physics_driver_init, physics_driver_down,   &
        physics_driver_down_time_vary, physics_driver_up_time_vary, &
        physics_driver_down_endts, physics_driver_up_endts, &
        physics_driver_up, physics_driver_end, &
        do_moist_in_phys_up, get_diff_t, &
        get_radturbten, zero_radturbten, physics_driver_restart, &
        cosp_driver_init, set_cosp_precip_sources

private          &

!  called from physics_driver_down:
         check_args, &

!  called from physics_driver_init:
         physics_driver_register_restart_domain,  physics_driver_register_restart_scalars, &

!  called from physics_driver_restart:
         physics_driver_netcdf, &

!  called from check_args:
         check_dim

interface check_dim
     module procedure check_dim_2d, check_dim_3d, check_dim_4d
end interface


!---------------------------------------------------------------------
!------- namelist ------
 
! <NAMELIST NAME="physics_driver_nml">
!  do_radiation  calculating radiative fluxes and  heating rates?
!  do_cosp       activate COSP simulator ?
!  do_modis_yim  activate simple modis simulator ?
!  do_moist_processes call moist_processes routines ?
!  tau_diff           time scale for smoothing diffusion coefficients
!  diff_min           minimum value of a diffusion coefficient beneath which the 
!                     coefficient is reset to zero
!  diffusion_smooth diffusion coefficients should be smoothed in time?
!  override_aerosols_cloud    use offline aerosols for cloud calculation
!                             (via data_override in aerosol_driver)?
!  qmin  UNITS="kg h2o/kg air"  minimum permissible value of cloud liquid, cloud ice, 
!                               saturated volume fraction, or rain and snow areas.
!               NOTE: qmin should be chosen such that the range of {qmin, max(qa,ql,qi)}
!                     is resolved by the precision of the numbers used.
!  N_land UNITS="1/(m*m*m)  assumed number of cloud drops per unit volume in liquid clouds
!                            over land when droplet number is not prdicted.
!  N_ocean UNITS="1/(m*m*m) assumed number of cloud drops per unit volume in liquid clouds
!                            over ocean when droplet number is not predicted.
!  do_liq_num  the prognostic droplet number option is activated ?
!  do_ice_num  the prognostic ice particle number option is activated ?
!  qcvar  1 / relative variance of sub-grid cloud water distribution
!         see morrison and gettelman, 2007, J. Climate for details
!  overlap
!      overlap        integer variable indicating which overlap 
!                     assumption to use:
!                     overlap = 1. means condensate in adjacent levels 
!                                  is treated as part of the same cloud
!                                  i.e. maximum-random overlap
!                     overlap = 2. means condensate in adjacent levels 
!                                  is treated as different clouds
!                                  i.e. random overlap
!  N_min         UNITS="1/(m*m*m) " minimum number of droplets allowed in a grid box when predicted
!                                   droplet number code is activated
!  min_diam_ice  UNITS="m" minimum size of ice crystals allowed in microphysics and radiation
!                          calculations
!  dcs            UNITS="m" ice crystal size at which autoconversion to falling ice occurs;
!                 used in radiation and microphysics calculations
!  min_diam_drop  UNITS="m"  minimum size of droplets allowed in microphysics and radiation
!                            calculations
!  max_diam_drop  UNITS="m"  maximum size of droplets allowed in microphysics and radiation
!                            calculations
!  use_tau   switch to determine whether current time level (tau)
!            will be used or else future time level (tau+1).
!            if use_tau = true then the input values for t,q, and r
!            are used; if use_tau = false then input values
!            tm+tdt*dt, etc. are used in moist_processes and vert_turb_driver.
!  cosp_frequency  UNITS="sec"  frequency at which the COSP simulator is to be called

logical :: do_radiation = .true.
logical :: do_cosp = .false.   
logical :: do_modis_yim = .true.
logical :: do_moist_processes = .true.
real    :: tau_diff = 3600.    
real    :: diff_min = 1.e-3   
logical :: diffusion_smooth = .true.
logical :: override_aerosols_cloud = .false.
real    :: qmin = 1.0e-10
real    :: N_land = 3.e8
real    :: N_ocean = 1.e8
logical :: do_liq_num = .true.
logical :: do_ice_num = .false.
real    :: qcvar = 1.
integer :: overlap = 2
real    :: N_min = 1.E6
real    :: min_diam_ice    = 10.e-6
real    :: min_diam_drop    = 2.e-6
real    :: max_diam_drop    = 50.e-6
real    :: dcs    =  200.e-6
logical :: use_tau = .false.
logical :: aer_in_phys_dn = .false.  ! ZNT 06/23/2023: Do aerosol in phys_dn instead of phys_up
real    :: cosp_frequency = 10800.


namelist / physics_driver_nml / do_radiation, do_cosp, &
                                do_modis_yim, &
                                do_moist_processes, tau_diff,      &
                                diff_min, diffusion_smooth, &
                                override_aerosols_cloud,    &
                                qmin, N_land, N_ocean, do_liq_num,  &
                                do_ice_num, qcvar, overlap, N_min, &
                                min_diam_ice, dcs, min_diam_drop, &
                                max_diam_drop, use_tau, cosp_frequency, &
                                aer_in_phys_dn


!---------------------------------------------------------------------
!------- public data ------
! 
public  surf_diff_type   ! defined in  vert_diff_driver_mod, republished
                         ! here

!---------------------------------------------------------------------
!------- private data ------
!--------------------------------------------------------------------
!    the following allocatable arrays are either used to hold physics 
!    data between timesteps when required, or hold physics data between
!    physics_down and physics_up.
!  
!    diff_cu_mo     contains contribution to difusion coefficient
!                   coming from cu_mo_trans_mod (called from 
!                   moist_processes in physics_driver_up) and then used 
!                   as input on the next time step to vert_diff_down 
!                   called in physics_driver_down.
!    diff_t         vertical diffusion coefficient for temperature
!                   which optionally may be time smoothed, meaning
!                   values must be saved between steps
!    diff_m         vertical diffusion coefficient for momentum
!                   which optionally may be time smoothed, meaning
!                   values must be saved between steps
!    radturbten     the sum of the radiational and turbulent heating,
!                   generated in both physics_driver_down (radiation)
!                   and physics_driver_up (turbulence) and then used
!                   in moist_processes
!    pbltop         top of boundary layer obtained from vert_turb_driver
!                   and then used on the next timestep in topo_drag_mod
!                   called from damping_driver_down        
!    cush
!    cbmf
!    convect        flag indicating whether convection is occurring in
!                   a grid column. generated in physics_driver_up and
!                   then used in vert_turb_driver called from 
!                   physics_driver_down on the next step.
!    temp_last
!    q_last
!----------------------------------------------------------------------
real,    dimension(:,:,:), allocatable,target :: diff_cu_mo, diff_t, diff_m
real,    dimension(:,:,:), allocatable,target :: radturbten, qturbten   ! ZNT: 06/08/2021
real,    dimension(:,:,:), allocatable,target :: tten_gw, cqa_gw ! P1L: 09/12/2024: grid-mean convective heating rate and updraft area to be used in gw_beres
real,    dimension(:,:,:), allocatable,target :: dqcdt_gw ! P1L: 10/21/2024: total large-scale condensation to be used in gw_beres
real,    dimension(:,:)  , allocatable,target :: pbltop, cush, cbmf
real,    dimension(:,:)  , allocatable,target :: hmint, cgust
real,    dimension(:,:)  , allocatable,target :: pblhto, rkmo, taudpo
logical, dimension(:,:)  , allocatable,target :: convect
integer, dimension(:,:,:), allocatable,target :: exist_shconv, exist_dpconv
real,    dimension(:,:,:), allocatable,target :: pblht_prev, hlsrc_prev, &
                                                 qtsrc_prev, cape_prev,  &
                                                 cin_prev

real,    dimension(:,:,:), allocatable        :: temp_last, q_last

integer                                :: vers
integer                                :: now_doing_strat = 0
integer                                :: now_doing_entrain = 0
real, allocatable                      :: r_convect(:,:)

type(aerosol_time_vary_type)           :: Aerosol_cld

!---------------------------------------------------------------------
!    internal timing clock variables:
!---------------------------------------------------------------------
integer :: damping_clock, turb_clock,   &
           tracer_clock, diff_up_clock, diff_down_clock, &
           moist_processes_clock, cosp_clock

!--------------------------------------------------------------------
!    miscellaneous control variables:
!---------------------------------------------------------------------
logical   :: do_check_args = .true.   ! argument dimensions should 
                                      ! be checked ?
logical   :: module_is_initialized = .false.
                                      ! module has been initialized ?
logical   :: doing_entrain            ! entrain_mod has been activated ?
logical   :: doing_uw_conv            ! uw_conv shallow cu mod has been 
                                      ! activated ?
logical   :: doing_liq_num = .false.  ! Prognostic cloud droplet number has 
                                      ! been activated?
integer   :: nt                       ! total no. of tracers
integer   :: ntp                      ! total no. of prognostic tracers
!integer   :: ncol                     ! number of stochastic columns
integer   ::  nsphum                  ! index for specific humidity tracer
 

logical   :: step_to_call_cosp
logical   :: doing_prog_clouds
real      :: rad_time_step

!---------------------------------------------------------------------
!---------------------------------------------------------------------

character(len=4)     :: mod_name = 'phys'
character(len=32)    :: tracer_units, tracer_name
character(len=128)   :: diaglname
real                 :: missing_value = -999.

integer                            :: id_tdt_phys,         &
                                      id_tdt_phys_vdif_dn, &
                                      id_tdt_phys_vdif_up, &
                                      id_tdt_phys_turb,    &
                                      id_tdt_phys_moist

integer, dimension(:), allocatable :: id_tracer_phys,         &
                                      id_tracer_phys_vdif_dn, &
                                      id_tracer_phys_vdif_up, &
                                      id_tracer_phys_turb,    &
                                      id_tracer_phys_moist

type(cmip_diag_id_type) :: ID_tntmp, ID_tnhusmp, &
                           ID_pfull, ID_phalf, ID_zg !, ID_zfull

type (clouds_from_moist_block_type) :: Restart

type(precip_flux_type)              :: Precip_flux

integer :: i_shallow
type (domain2D)               :: physics_domain !< Atmosphere domain

                            contains



!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!                     PUBLIC SUBROUTINES
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!#####################################################################
! <SUBROUTINE NAME="physics_driver_init">
!  <OVERVIEW>
!    physics_driver_init is the constructor for physics_driver_mod.
!  </OVERVIEW>
!  <DESCRIPTION>
!    physics_driver_init is the constructor for physics_driver_mod.
!  </DESCRIPTION>
!  <TEMPLATE>
!   call physics_driver_init (Time, lonb, latb, axes, pref, &
!                             trs, Surf_diff, phalf )
!  </TEMPLATE>
!  <IN NAME="Time" TYPE="time_type">
!   current time
!  </IN>
!  <IN NAME="pref" TYPE="real">
!   reference prssure profiles
!  </IN>
!  <IN NAME="latb" TYPE="real">
!   array of model latitudes at cell corners [radians]
!  </IN>
!  <IN NAME="lonb" TYPE="real">
!   array of model longitudes at cell corners [radians]
!  </IN>
!  <IN NAME="axes" TYPE="integer">
!   axis indices, (/x,y,pf,ph/)
!                (returned from diag axis manager)
!  </IN>
!  <INOUT NAME="trs" TYPE="real">
!   atmospheric tracer fields
!  </INOUT>
!  <INOUT NAME="Surf_diff" TYPE="surf_diff_type">
!   surface diffusion derived type
!  </INOUT>
!  <IN NAME="phalf" TYPE="real">
!   pressure at model interface levels
!  </IN>
! <ERROR MSG="physics_driver_init must be called first" STATUS="FATAL">
! </ERROR>
! </SUBROUTINE>
!
subroutine physics_driver_init (Time, lonb, latb, lon, lat, axes, &
                                Surf_diff, Exch_ctrl, Atm_block,   &
                                Moist_clouds, Physics, Physics_tendency, &
                                diffm, difft)

!---------------------------------------------------------------------
!    physics_driver_init is the constructor for physics_driver_mod.
!---------------------------------------------------------------------

type(time_type),              intent(in)    :: Time
real,    dimension(:,:),      intent(in)    :: lonb, latb
real,    dimension(:,:),      intent(in)    :: lon, lat
integer, dimension(4),        intent(in)    :: axes
type(surf_diff_type),         intent(inout) :: Surf_diff
type (exchange_control_type), intent(inout) :: Exch_ctrl
type (block_control_type),    intent(in)    :: Atm_block
type(clouds_from_moist_type), intent(inout) :: Moist_clouds(:)
type(physics_type),           intent(inout) :: Physics
type(physics_tendency_type),  intent(inout) :: Physics_tendency
real,    dimension(:,:,:),    intent(out),  optional :: diffm, difft

!---------------------------------------------------------------------
!  intent(in) variables:
!     Time       current time (time_type)
!     lonb       longitude of the grid box corners [ radians ]
!     latb       latitude of the grid box corners [ radians ]
!     lon
!     lat
!     axes       axis indices, (/x,y,pf,ph/)
!                (returned from diag axis manager)
!
!   intent(inout) variables:
!     Surf_diff  surface diffusion derived type variable
!
!   intent(out), optional variables:
!
!     diffm
!     difft
!
!---------------------------------------------------------------------

!---------------------------------------------------------------------
!  local variables:

      integer :: nb, ibs, ibe, jbs, jbe
      real, dimension (size(lonb,1)-1, size(latb,2)-1) :: sgsmtn
      integer          ::  id, jd, kd, n, k, nc
      integer          ::  ierr, io, unit, logunit, outunit
      integer          ::  ndum

      integer          ::  moist_processes_init_clock, damping_init_clock,&
                           turb_init_clock, diff_init_clock, &
                           aerosol_init_clock, &
                           tracer_init_clock
      real, dimension(:,:,:),   allocatable :: phalf
      real, dimension(:,:,:,:), allocatable :: trs
      type(FmsNetcdfFile_t)       ::  Phy_restart !< Fms2io fileobj
      type(FmsNetcdfDomainFile_t) ::  Til_restart !< Fms2io domain decomposed fileobj
      integer, dimension(:),    allocatable :: pes !< Array of pes in the current pelist
!---------------------------------------------------------------------
!  local variables:
!
!       sgsmtn        sgs orography obtained from mg_drag_mod;
!                     appears to not be currently used
!       id,jd,kd      model dimensions on the processor  
!       n             loop index
!       ierr          error code
!       io            io status returned from an io call
!       unit          unit number used for an i/ operation
!       logunit
!       outunit
!       ndum          dummy argument
!       x_clock_init  clock for timing the initialization of process x
!                     where x is moist_processes, damping, turb, diff,
!                     aerosol, tracer 
!-----------------------------------------------------------------------

!---------------------------------------------------------------------
!    if routine has already been executed, return.
!---------------------------------------------------------------------
      if (module_is_initialized) return

!---------------------------------------------------------------------
!    verify that the modules used by this module that are not called 
!    later in this subroutine have already been initialized.
!---------------------------------------------------------------------
      call fms_init
      call time_manager_init
      call tracer_manager_init
      call field_manager_init (ndum)
 
!--------------------------------------------------------------------
!    read namelist.
!--------------------------------------------------------------------
      read (input_nml_file, nml=physics_driver_nml, iostat=io)
      ierr = check_nml_error(io,"physics_driver_nml")

!--------------------------------------------------------------------
!    consistency checks for namelist options
!--------------------------------------------------------------------
      if (do_cosp .and. .not. do_radiation) &
        call error_mesg('physics_driver_init',  &
            'do_radiation must be .true. if do_cosp is .true.',FATAL)
      if (do_cosp .and. .not. do_moist_processes) &
         call error_mesg('physics_driver_init',  &
         'do_moist_processes must be .true. if do_cosp is .true.', FATAL)

!--------------------------------------------------------------------
!    write version number and namelist to log file.
!--------------------------------------------------------------------
      call write_version_number (version, tagname)
      logunit = stdlog()
      if (mpp_pe() == mpp_root_pe() ) &
               write(logunit, nml=physics_driver_nml)
 
!---------------------------------------------------------------------
!    define the model dimensions on the local processor (id, jd, kd). 
!    retrieve the total number of tracers (nt) and prognostic 
!    tracers (ntp). Save the number of prognostic tracers in a 
!    physics_control_type for use in other modules.
!---------------------------------------------------------------------
      id = size(lonb,1)-1 
      jd = size(latb,2)-1 
      kd = Atm_block%npz
      call get_number_tracers (MODEL_ATMOS, num_tracers=nt, &
                               num_prog=ntp)
      Physics%control%num_prog_tracers = ntp

!---------------------------------------------------------------------
      cosp_clock       =       &
                mpp_clock_id( '   Physics_down: COSP',    &
                   grain=CLOCK_MODULE_DRIVER )
      damping_clock         =     &
                mpp_clock_id( '   Physics_down: Damping',    &
                  grain=CLOCK_MODULE_DRIVER )
      turb_clock            =      &
                mpp_clock_id( '   Physics_down: Vert. Turb.', &
                  grain=CLOCK_MODULE_DRIVER )
      tracer_clock          =      &
                mpp_clock_id( '   Physics_down: Tracer',    &
                 grain=CLOCK_MODULE_DRIVER )
      diff_down_clock       =     &
                mpp_clock_id( '   Physics_down: Vert. Diff.',   &
                 grain=CLOCK_MODULE_DRIVER )
      diff_up_clock         =     &
                mpp_clock_id( '   Physics_up: Vert. Diff.',     &
                grain=CLOCK_MODULE_DRIVER )
      moist_processes_clock =      &
                mpp_clock_id( '   Physics_up: Moist Processes', &
                grain=CLOCK_MODULE_DRIVER )

      moist_processes_init_clock =      &
        mpp_clock_id( '   Physics_driver_init: Moist Processes: Initialization', &
                grain=CLOCK_MODULE_DRIVER )
      damping_init_clock         =     &
        mpp_clock_id( '   Physics_driver_init: Damping: Initialization',    &
                  grain=CLOCK_MODULE_DRIVER )
      turb_init_clock            =      &
        mpp_clock_id( '   Physics_driver_init: Vert. Turb.: Initialization', &
                  grain=CLOCK_MODULE_DRIVER )
      diff_init_clock       =     &
        mpp_clock_id( '   Physics_driver_init: Vert. Diff.: Initialization',   &
                 grain=CLOCK_MODULE_DRIVER )
      aerosol_init_clock       =       &
        mpp_clock_id( '   Physics_driver_init: Aerosol: Initialization', &
                       grain=CLOCK_MODULE_DRIVER )
      tracer_init_clock          =      &
        mpp_clock_id( '   Physics_driver_init: Tracer: Initialization',    &
                 grain=CLOCK_MODULE_DRIVER )


!-----------------------------------------------------------------------
!     dummy checks
!-----------------------------------------------------------------------
      if (do_ice_num .and. .not. do_liq_num) then
        call error_mesg ('physics_driver_mod',  &
           'do_ice_num can only be .true. if do_liq_num is .true.', FATAL)
      endif

!------------------------------------------------------------------------
!   place some control variables needed in multiple physics modules
!   into the Physics%control derived type for easy movement.
!------------------------------------------------------------------------
      Physics%control%use_tau = use_tau
      Physics%control%nsphum = get_tracer_index ( MODEL_ATMOS, 'sphum' )
      Physics%control%nql = get_tracer_index ( MODEL_ATMOS, 'liq_wat' )
      Physics%control%nqi = get_tracer_index ( MODEL_ATMOS, 'ice_wat' )
      Physics%control%nqa = get_tracer_index ( MODEL_ATMOS, 'cld_amt' )
      Physics%control%nqn = get_tracer_index ( MODEL_ATMOS, 'liq_drp' )
      Physics%control%nqni = get_tracer_index ( MODEL_ATMOS, 'ice_num' )
      Physics%control%nqr = get_tracer_index (MODEL_ATMOS, 'rainwat')
      Physics%control%nqs = get_tracer_index (MODEL_ATMOS, 'snowwat')
      Physics%control%nqg = get_tracer_index (MODEL_ATMOS, 'graupel')

      Physics%control%nqnr = get_tracer_index (MODEL_ATMOS, 'rain_num')
      Physics%control%nqns = get_tracer_index (MODEL_ATMOS, 'snow_num')

      physics_domain = Physics%control%domain
!-----------------------------------------------------------------------
!   allocate a logical array to define whether a tracer is a cloud tracer
!   (one of those defined above), or not.
!-----------------------------------------------------------------------
      allocate (Physics%control%cloud_tracer   &
                                   (Physics%control%num_prog_tracers))
      Physics%control%cloud_tracer = .FALSE. 

      if (Physics%control%nsphum /= NO_TRACER) then
        Physics%control%cloud_tracer(Physics%control%nsphum) = .TRUE.
      endif
      if (Physics%control%nql    /= NO_TRACER) then
        Physics%control%cloud_tracer(Physics%control%nql   ) = .TRUE.
      endif
      if (Physics%control%nqi    /= NO_TRACER) then
        Physics%control%cloud_tracer(Physics%control%nqi   ) = .TRUE.
      endif
      if (Physics%control%nqa    /= NO_TRACER) then
        Physics%control%cloud_tracer(Physics%control%nqa   ) = .TRUE.
      endif
      if (Physics%control%nqn    /= NO_TRACER) then
        Physics%control%cloud_tracer(Physics%control%nqn   ) = .TRUE.
      endif
      if (Physics%control%nqni   /= NO_TRACER) then
        Physics%control%cloud_tracer(Physics%control%nqni  ) = .TRUE.
      endif
      if (Physics%control%nqr    /= NO_TRACER) then
        Physics%control%cloud_tracer(Physics%control%nqr   ) = .TRUE.
      endif
      if (Physics%control%nqs    /= NO_TRACER) then
        Physics%control%cloud_tracer(Physics%control%nqs   ) = .TRUE.
      endif
      if (Physics%control%nqg    /= NO_TRACER) then
        Physics%control%cloud_tracer(Physics%control%nqg   ) = .TRUE.
      endif

      if (Physics%control%nqnr    /= NO_TRACER) then
        Physics%control%cloud_tracer(Physics%control%nqnr   ) = .TRUE.
      endif
      if (Physics%control%nqns    /= NO_TRACER) then
        Physics%control%cloud_tracer(Physics%control%nqns   ) = .TRUE.
      endif
      
!----------------------------------------------------------------------
!   define logical variable indicating whether prognostic clouds (using 
!   tracer fields) are active. 
!----------------------------------------------------------------------
      if (min(Physics%control%nql, Physics%control%nqi,   &
                                            Physics%control%nqa) > 0) then
        doing_prog_clouds = .true.
      else
        doing_prog_clouds = .false.
      endif
      
!----------------------------------------------------------------------
!do some dummy checks on the tracer indices.
!----------------------------------------------------------------------
      if (doing_prog_clouds) then
        if (min(Physics%control%nql, Physics%control%nqi,    &
                                           Physics%control%nqa) <= 0) &
          call error_mesg ('physics_driver_init', &
                        'stratiform cloud tracer(s) not found', FATAL)
        if (Physics%control%nql == Physics%control%nqi .or.   &
            Physics%control%nqa == Physics%control%nqi .or.  &
            Physics%control%nql == Physics%control%nqa)    & 
          call error_mesg ('physics_driver_init',  &
          'tracers indices cannot be the same (i.e., nql=nqi=nqa).', FATAL)
        if (mpp_pe() == mpp_root_pe()) &
            write (logunit,'(a,3i4)')   &
              'Stratiform cloud tracer indices: nql,nqi,nqa =',  &
                  Physics%control%nql, Physics%control%nqi,   &
                                                    Physics%control%nqa
                                  
        if (Physics%control%nqn == NO_TRACER .and.   &
                              Exch_ctrl%do_liq_num ) &
           call error_mesg ('physics_driver_init', &
                     'prognostic droplet number scheme requested but&
                                           &  tracer not found', FATAL)
        if (Physics%control%nqni == NO_TRACER .and.   &
                                  Exch_ctrl%do_ice_num ) &
                call error_mesg ('physics_driver_init', &
                     'prognostic ice number scheme requested but &
                                              &tracer not found', FATAL)
      endif

!----------------------------------------------------------------------
!    place the physics_driver_nml variables that are needed by both the 
!    moist_processes and radiation codes into Exch_ctrl.
!----------------------------------------------------------------------
      Exch_ctrl%cosp_frequency = cosp_frequency
      Exch_ctrl%N_min = N_min
      Exch_ctrl%qmin = qmin
      Exch_ctrl%qcvar   = qcvar
      Exch_ctrl%N_land = N_land
      Exch_ctrl%N_ocean = N_ocean
      if (overlap.ne.1 .and. overlap.ne.2) &
          call error_mesg  ('physics_driver_init',&
                               'overlap must be either 1 or 2 ', FATAL)
      Exch_ctrl%overlap = overlap
      Exch_ctrl%do_liq_num = do_liq_num
      Exch_ctrl%do_ice_num = do_ice_num

      Exch_ctrl%min_diam_ice = min_diam_ice
      Exch_ctrl%dcs          = dcs            
      Exch_ctrl%min_diam_drop = min_diam_drop
      Exch_ctrl%max_diam_drop = max_diam_drop

      Exch_ctrl%do_cosp = do_cosp
      Exch_ctrl%do_modis_yim = do_modis_yim
      Exch_ctrl%doing_prog_clouds = doing_prog_clouds

!-----------------------------------------------------------------------
!--- allocate Physics_tendency to hold the physics tendencies
!-----------------------------------------------------------------------
      call alloc_physics_tendency_type (Physics_tendency, Atm_block)

!--- define trs and p_half on the full domain 
      allocate (trs(id,jd,kd,nt), phalf(id,jd,kd+1))
      do nb = 1, Atm_block%nblks
        ibs = Atm_block%ibs(nb)-Atm_block%isc+1
        ibe = Atm_block%ibe(nb)-Atm_block%isc+1
        jbs = Atm_block%jbs(nb)-Atm_block%jsc+1
        jbe = Atm_block%jbe(nb)-Atm_block%jsc+1
        trs(ibs:ibe,jbs:jbe,:,1:ntp)    = Physics%block(nb)%q
        trs(ibs:ibe,jbs:jbe,:,ntp+1:nt) = Physics%block(nb)%tmp_4d
        phalf(ibs:ibe,jbs:jbe,:)        = Physics%block(nb)%p_half
!--- the 'temp' variable inside of Physics is no longer needed - deallocate it
        deallocate(Physics%block(nb)%tmp_4d)
      enddo

!-----------------------------------------------------------------------
!---------- initialize physics -------

      if (do_moist_processes) then
        call mpp_clock_begin ( moist_processes_init_clock )
        call moist_processes_init (physics_domain, id, jd, kd, lonb, latb, lon, lat,  &
                                   phalf, Physics%glbl_qty%pref(:,1),&
                                   axes, Time, Physics%control, Exch_ctrl) 

        call mpp_clock_end ( moist_processes_init_clock )
      else
        diff_cu_mo = 0.0
        convect = .false.
      endif
     
!-----------------------------------------------------------------------
!    initialize damping_driver_mod.
!-----------------------------------------------------------------------
      call mpp_clock_begin ( damping_init_clock )
      call damping_driver_init (physics_domain, lonb, latb, Physics%glbl_qty%pref(:,1), &
                                axes, Time, sgsmtn)
      call mpp_clock_end ( damping_init_clock )

!-----------------------------------------------------------------------
!    initialize vert_turb_driver_mod.
!-----------------------------------------------------------------------
      call mpp_clock_begin ( turb_init_clock )
      call vert_turb_driver_init (physics_domain, lonb, latb, id, jd, kd, axes, Time, &
                                  Exch_ctrl, Physics%control,  &
                                  doing_entrain)
      call mpp_clock_end ( turb_init_clock )

!-----------------------------------------------------------------------
!    initialize vert_diff_driver_mod.
!-----------------------------------------------------------------------
      call mpp_clock_begin ( diff_init_clock )
      call vert_diff_driver_init (Surf_diff, id, jd, kd, axes, Time )
      call mpp_clock_end ( diff_init_clock )

      ! if (do_moist_processes) then   ! ZNT 06/23/2023: unconditionally initialize aerosol because EDMF needs it
!-----------------------------------------------------------------------
!    initialize aerosol_mod     
!-----------------------------------------------------------------------
      call mpp_clock_begin ( aerosol_init_clock )
      call aerosol_init (lonb, latb, Aerosol_cld)
      call mpp_clock_end ( aerosol_init_clock )
      ! endif ! do_moist_processes

!-----------------------------------------------------------------------
!    initialize atmos_tracer_driver_mod.
!-----------------------------------------------------------------------
      call mpp_clock_begin ( tracer_init_clock )
      call atmos_tracer_driver_init (physics_domain, lonb, latb, trs, axes, Time, phalf)
      call mpp_clock_end ( tracer_init_clock )

!---------------------------------------------------------------------
!    allocate space for the module variables and initialize them.
!---------------------------------------------------------------------
      allocate ( diff_t     (id, jd, kd) ) ; diff_t = 0.0
      allocate ( diff_m     (id, jd, kd) ) ; diff_m = 0.0
      allocate ( diff_cu_mo (id, jd, kd) ) ; diff_cu_mo = 0.0
      allocate ( pbltop     (id, jd) )     ; pbltop     = -999.0
      allocate ( cush       (id, jd) )     ; cush=-1. !miz
      allocate ( cbmf       (id, jd) )     ; cbmf=0.0 !miz
      allocate ( hmint      (id, jd) )     ; hmint=0. !miz
      allocate ( cgust      (id, jd) )     ; cgust=0.0 !miz
      allocate ( pblhto     (id, jd) )     ; pblhto=0.0 !miz
      allocate ( rkmo       (id, jd) )     ; rkmo=15.0 !miz
      allocate ( taudpo     (id, jd) )     ; taudpo=28800.    !miz
      allocate ( exist_shconv(id, jd,48) ) ; exist_shconv = 0 !miz
      allocate ( exist_dpconv(id, jd,48) ) ; exist_dpconv = 0 !miz
      allocate ( pblht_prev  (id, jd,48) ) ; pblht_prev   = 0.!miz
      allocate ( hlsrc_prev  (id, jd,48) ) ; hlsrc_prev   = 0.!miz
      allocate ( qtsrc_prev  (id, jd,48) ) ; qtsrc_prev   = 0.!miz
      allocate ( cape_prev   (id, jd,48) ) ; cape_prev    = 0.!miz
      allocate ( cin_prev    (id, jd,48) ) ; cin_prev     = 0.!miz

      allocate ( convect    (id, jd) )     ; convect = .false.
      allocate ( radturbten (id, jd, kd))  ; radturbten = 0.0
      allocate ( qturbten   (id, jd, kd))  ; qturbten = 0.0   ! ZNT: 06/08/2021
      allocate ( r_convect  (id, jd) )     ; r_convect   = 0.0
      allocate ( tten_gw (id, jd, kd) ) ; tten_gw = 0.0 ! P1L: 03/12/2024
      allocate ( cqa_gw (id, jd, kd) ) ; cqa_gw = 0.0 ! P1L: 09/12/2024
      allocate ( dqcdt_gw (id, jd, kd) ) ; dqcdt_gw = 0.0 ! P1L: 10/21/2024

      if (do_cosp) then
!--------------------------------------------------------------------
!    these variables are needed to preserve values of rain fluxes, q and T
!    from the step preceding the COSP call for use in the COSP simulator 
!    on the next step.
!--------------------------------------------------------------------
        allocate ( Precip_flux%fl_lsrain  (id, jd, kd))
        allocate ( Precip_flux%fl_lssnow  (id, jd, kd))
        allocate ( Precip_flux%fl_lsgrpl  (id, jd, kd))
        allocate ( Precip_flux%fl_ccrain  (id, jd, kd))
        allocate ( Precip_flux%fl_ccsnow  (id, jd, kd))
        allocate ( temp_last (id, jd, kd))
        allocate ( q_last    (id, jd, kd))
        Precip_flux%fl_lsrain = 0.
        Precip_flux%fl_lssnow = 0.
        Precip_flux%fl_lsgrpl = 0.
        Precip_flux%fl_ccrain = 0.
        Precip_flux%fl_ccsnow = 0.
        temp_last = 0.
        q_last    = 0.
      endif

!-----------------------------------------------------------------------
!    allocate derived-type that stores cloud properties
!    return from moist processes
!-----------------------------------------------------------------------
     call error_mesg('physics_driver_mod', 'number of cloud schemes found = '//trim(string(Exch_ctrl%ncld)), NOTE)

     call alloc_clouds_from_moist_type(Moist_clouds, Exch_ctrl, Atm_block)

!------------------------------------------------------------------------
!    save convective cloud indices to be passed to convection_driver_mod.
!------------------------------------------------------------------------
     i_shallow = Moist_clouds(1)%block(1)%index_uw_conv

!--------------------------------------------------------------------
!    call physics_driver_read_restart to obtain initial values for the module
!    variables. Also register restart fields to be ready for intermediate 
!    restart.
!--------------------------------------------------------------------
      allocate(Restart%Cloud_data(Exch_ctrl%ncld))

      do nc = 1, Exch_ctrl%ncld
        ! restart values allocated on the full domain
        call alloc_cloud_scheme_data_type( Moist_clouds(1)%block(1)%Cloud_data(nc)%scheme_name, &
                                           id, jd, kd, Restart%Cloud_data(nc))
      enddo

      !< Get the current pelist
      allocate(pes(mpp_npes()))
      call mpp_get_current_pelist(pes)

      !< Open the scalar file with the current pelist, so that only the root pe opens and reads the file and
      !! distributes the data to the other pes
      if (open_file(Phy_restart,"INPUT/physics_driver.res.nc","read", is_restart=.true., pelist=pes)) then !scalar file
         call physics_driver_register_restart_scalars(Restart, Phy_restart)
         call read_restart(Phy_restart)
         call close_file(Phy_restart)
      endif
      deallocate(pes)

      if (open_file(Til_restart,"INPUT/physics_driver.res.nc","read", physics_domain, is_restart=.true.)) then !domain file
         call physics_driver_register_restart_domain(Restart, Til_restart)
         call read_restart(Til_restart)
         call close_file(Til_restart)
      endif

!---------------------------------------------------------------------
!    convert the real variable (r_convect) indicating columns with 
!    convection to a logical variable (convect). this will be used in 
!    vert_turb_driver_mod.
!---------------------------------------------------------------------
      convect = .false.
      where(r_convect .GT. 0.) 
         convect = .true.
      end where
         
100 FORMAT("CHECKSUM::",A32," = ",Z20)
      outunit = stdout()
      write(outunit,*) 'BEGIN CHECKSUM(physics_driver_init):: '
      write(outunit,100) 'diff_cu_mo             ', mpp_chksum(diff_cu_mo            )
      write(outunit,100) 'pbltop                 ', mpp_chksum(pbltop                )
      write(outunit,100) 'cush                   ', mpp_chksum(cush                  )
      write(outunit,100) 'cbmf                   ', mpp_chksum(cbmf                  )
      write(outunit,100) 'hmint                  ', mpp_chksum(hmint                 )
      write(outunit,100) 'cgust                  ', mpp_chksum(cgust                 )
      write(outunit,100) 'pblhto                 ', mpp_chksum(pblhto                )
      write(outunit,100) 'rkmo                   ', mpp_chksum(rkmo                  )
      write(outunit,100) 'taudpo                 ', mpp_chksum(taudpo                )
      write(outunit,100) 'exist_shconv           ', mpp_chksum(exist_shconv          )
      write(outunit,100) 'exist_dpconv           ', mpp_chksum(exist_dpconv          )
      write(outunit,100) 'pblht_prev             ', mpp_chksum(pblht_prev            )
      write(outunit,100) 'hlsrc_prev             ', mpp_chksum(hlsrc_prev            )
      write(outunit,100) 'qtsrc_prev             ', mpp_chksum(qtsrc_prev            )
      write(outunit,100) 'cape_prev              ', mpp_chksum(cape_prev             )
      write(outunit,100) 'cin_prev               ', mpp_chksum(cin_prev              )
      write(outunit,100) 'diff_t                 ', mpp_chksum(diff_t                )
      write(outunit,100) 'diff_m                 ', mpp_chksum(diff_m                )
      write(outunit,100) 'r_convect              ', mpp_chksum(r_convect             )
      write(outunit,100) 'tten_gw              ', mpp_chksum(tten_gw               ) ! P1L
      write(outunit,100) 'cqa_gw               ', mpp_chksum(cqa_gw                ) ! P1L
      write(outunit,100) 'dqcdt_gw             ', mpp_chksum(dqcdt_gw              ) ! P1L
   if (doing_prog_clouds) then
      write(outunit,100) 'radturbten             ', mpp_chksum(radturbten            )
      write(outunit,100) 'qturbten               ', mpp_chksum(qturbten              )  ! ZNT
   endif
      do nc = 1, size(Restart%Cloud_data,1)
        ! NOTE: the order of the checksums in stdout will be different
        if ( trim(Restart%Cloud_data(nc)%scheme_name).eq.'uw_conv' ) then
          write(outunit,100) 'shallow_cloud_area     ', mpp_chksum(Restart%Cloud_data(nc)%cloud_area    )
          write(outunit,100) 'shallow_liquid         ', mpp_chksum(Restart%Cloud_data(nc)%liquid_amt    )
          write(outunit,100) 'shallow_ice            ', mpp_chksum(Restart%Cloud_data(nc)%ice_amt       )
          write(outunit,100) 'shallow_droplet_number ', mpp_chksum(Restart%Cloud_data(nc)%droplet_number)
          write(outunit,100) 'shallow_ice_number     ', mpp_chksum(Restart%Cloud_data(nc)%ice_number    )
        endif
        if ( trim(Restart%Cloud_data(nc)%scheme_name).eq.'strat_cloud' ) then
          write(outunit,100) 'lsc_cloud_area         ', mpp_chksum(Restart%Cloud_data(nc)%cloud_area    )
          write(outunit,100) 'lsc_liquid             ', mpp_chksum(Restart%Cloud_data(nc)%liquid_amt    )
          write(outunit,100) 'lsc_ice                ', mpp_chksum(Restart%Cloud_data(nc)%ice_amt       )
          write(outunit,100) 'lsc_droplet_number     ', mpp_chksum(Restart%Cloud_data(nc)%droplet_number)
          write(outunit,100) 'lsc_ice_number         ', mpp_chksum(Restart%Cloud_data(nc)%ice_number    )
          write(outunit,100) 'lsc_snow               ', mpp_chksum(Restart%Cloud_data(nc)%snow          )
          write(outunit,100) 'lsc_rain               ', mpp_chksum(Restart%Cloud_data(nc)%rain          )
          write(outunit,100) 'lsc_snow_size          ', mpp_chksum(Restart%Cloud_data(nc)%snow_size     )
          write(outunit,100) 'lsc_rain_size          ', mpp_chksum(Restart%Cloud_data(nc)%rain_size     )
        endif
      enddo ! nc

      do nb = 1, Atm_block%nblks
        ibs = Atm_block%ibs(nb)-Atm_block%isc+1
        ibe = Atm_block%ibe(nb)-Atm_block%isc+1
        jbs = Atm_block%jbs(nb)-Atm_block%jsc+1
        jbe = Atm_block%jbe(nb)-Atm_block%jsc+1

        !-- copy cloud data from restart
        do nc = 1, size(Restart%Cloud_data,1)

          ! common to all cloud schemes
          Moist_clouds(1)%block(nb)%Cloud_data(nc)%cloud_area     = Restart%Cloud_data(nc)%cloud_area    (ibs:ibe,jbs:jbe,:)
          Moist_clouds(1)%block(nb)%Cloud_data(nc)%liquid_amt     = Restart%Cloud_data(nc)%liquid_amt    (ibs:ibe,jbs:jbe,:)
          Moist_clouds(1)%block(nb)%Cloud_data(nc)%ice_amt        = Restart%Cloud_data(nc)%ice_amt       (ibs:ibe,jbs:jbe,:)
          Moist_clouds(1)%block(nb)%Cloud_data(nc)%droplet_number = Restart%Cloud_data(nc)%droplet_number(ibs:ibe,jbs:jbe,:)
          Moist_clouds(1)%block(nb)%Cloud_data(nc)%scheme_name    = Restart%Cloud_data(nc)%scheme_name

          ! properties specific to large-scale/stratiform clouds
          if (trim(Restart%Cloud_data(nc)%scheme_name) .eq. 'strat_cloud') then
            Moist_clouds(1)%block(nb)%Cloud_data(nc)%ice_number = Restart%Cloud_data(nc)%ice_number (ibs:ibe,jbs:jbe,:)
            Moist_clouds(1)%block(nb)%Cloud_data(nc)%rain       = Restart%Cloud_data(nc)%rain       (ibs:ibe,jbs:jbe,:)
            Moist_clouds(1)%block(nb)%Cloud_data(nc)%snow       = Restart%Cloud_data(nc)%snow       (ibs:ibe,jbs:jbe,:)
            Moist_clouds(1)%block(nb)%Cloud_data(nc)%rain_size  = Restart%Cloud_data(nc)%rain_size  (ibs:ibe,jbs:jbe,:)
            Moist_clouds(1)%block(nb)%Cloud_data(nc)%snow_size  = Restart%Cloud_data(nc)%snow_size  (ibs:ibe,jbs:jbe,:)
          endif
  
          ! properties specific to uw shallow convective clouds
          if (trim(Restart%Cloud_data(nc)%scheme_name) .eq. 'uw_conv') then
            Moist_clouds(1)%block(nb)%Cloud_data(nc)%ice_number = Restart%Cloud_data(nc)%ice_number (ibs:ibe,jbs:jbe,:)
          endif

        enddo

        !--- return trs to the blocked data structure
        Physics%block(nb)%q = trs(ibs:ibe,jbs:jbe,:,1:ntp)
        Physics_tendency%block(nb)%qdiag = trs(ibs:ibe,jbs:jbe,:,ntp+1:nt)
      enddo
      deallocate (trs, phalf)

!---------------------------------------------------------------------
!    if desired, define variables to return  the restart fields of 
!    diff_m and diff_t.
!---------------------------------------------------------------------
      if (present(difft)) then
        difft = diff_t
      endif
      if (present(diffm)) then
        diffm = diff_m
      endif

!---------------------------------------------------------------------
!    register module diagnostics
!---------------------------------------------------------------------

      id_tdt_phys_vdif_dn = register_diag_field ( mod_name,    &
         'tdt_phys_vdif_dn', axes(1:3), Time,                  &
         'temperature tendency from physics driver vdif down', &
         'K/s', missing_value=missing_value)

      id_tdt_phys_vdif_up = register_diag_field ( mod_name,    &
         'tdt_phys_vdif_up', axes(1:3), Time,                  &
         'temperature tendency from physics driver vdif up',   &
         'K/s', missing_value=missing_value)

      id_tdt_phys_turb = register_diag_field ( mod_name,       &
         'tdt_phys_turb', axes(1:3), Time,                     &
         'temperature tendency from physics driver vdif turb', &
         'K/s', missing_value=missing_value)

      id_tdt_phys_moist = register_diag_field ( mod_name,            &
         'tdt_phys_moist', axes(1:3), Time,                          &
         'temperature tendency from physics driver moist processes', &
         'K/s', missing_value=missing_value)

      id_tdt_phys = register_diag_field ( mod_name,            &
         'tdt_phys', axes(1:3), Time,                          &
         'temperature tendency from physics ', &
         'K/s', missing_value=missing_value)

     !-------- CMIP diagnostics --------
      ID_pfull = register_cmip_diag_field_3d ( mod_name, 'pfull', Time, &
                                     'Pressure on Model Levels', 'Pa', &
                                      standard_name = 'air_pressure')

      ID_phalf = register_cmip_diag_field_3d ( mod_name, 'phalf', Time, &
                                'Pressure on Model Half-Levels', 'Pa', &
                              standard_name='air_pressure', axis='half')

      ID_zg = register_cmip_diag_field_3d ( mod_name, 'zg', Time, &
                                      'Geopotential Height', 'm', &
                             standard_name = 'geopotential_height')

     !ID_zfull = register_cmip_diag_field_3d ( mod_name, 'zfull', Time, &
     !                            'Altitude of Model Full-Levels', 'm', &
     !                standard_name = 'height_above_reference_ellipsoid')

     !-------- CMIP diagnostics (tendencies due to physics) --------
      ID_tntmp = register_cmip_diag_field_3d ( mod_name, 'tntmp', Time, &
                  'Tendency of Air Temperature due to Model Physics', 'K s-1', & 
             standard_name='tendency_of_air_temperature_due_to_model_physics' )

      nsphum = get_tracer_index(MODEL_ATMOS,'sphum')
      if (nsphum /= NO_TRACER) then
        ID_tnhusmp = register_cmip_diag_field_3d ( mod_name, 'tnhusmp', Time, &
                    'Tendency of Specific Humidity due to Model Physics', 's-1', &
               standard_name='tendency_of_specific_humidity_due_to_model_physics' )
      endif

      allocate (id_tracer_phys(ntp))
      allocate (id_tracer_phys_vdif_dn(ntp))
      allocate (id_tracer_phys_vdif_up(ntp))
      allocate (id_tracer_phys_turb(ntp))
      allocate (id_tracer_phys_moist(ntp))

      do n = 1,ntp

        call get_tracer_names (MODEL_ATMOS, n, name = tracer_name,  &
                               units = tracer_units)
        
        diaglname = trim(tracer_name)//  &
                    ' tendency from physics'
        id_tracer_phys(n) =    &
                         register_diag_field ( mod_name, &
                         TRIM(tracer_name)//'_phys',  &
                         axes(1:3), Time, trim(diaglname), &
                         TRIM(tracer_units)//'/s',  &
                         missing_value=missing_value)

        diaglname = trim(tracer_name)//  &
                    ' tendency from physics driver vdif down'
        id_tracer_phys_vdif_dn(n) =    &
                         register_diag_field ( mod_name, &
                         TRIM(tracer_name)//'_phys_vdif_dn',  &
                         axes(1:3), Time, trim(diaglname), &
                         TRIM(tracer_units)//'/s',  &
                         missing_value=missing_value)

        diaglname = trim(tracer_name)//  &
                    ' tendency from physics driver vdif up'
        id_tracer_phys_vdif_up(n) =    &
                         register_diag_field ( mod_name, &
                         TRIM(tracer_name)//'_phys_vdif_up',  &
                         axes(1:3), Time, trim(diaglname), &
                         TRIM(tracer_units)//'/s',  &
                         missing_value=missing_value)

        diaglname = trim(tracer_name)//  &
                    ' tendency from physics driver vert turb'
        id_tracer_phys_turb(n) =    &
                         register_diag_field ( mod_name, &
                         TRIM(tracer_name)//'_phys_turb',  &
                         axes(1:3), Time, trim(diaglname), &
                         TRIM(tracer_units)//'/s',  &
                         missing_value=missing_value)

        diaglname = trim(tracer_name)//  &
                    ' tendency from physics driver moist processes'
        id_tracer_phys_moist(n) =    &
                         register_diag_field ( mod_name, &
                         TRIM(tracer_name)//'_phys_moist',  &
                         axes(1:3), Time, trim(diaglname), &
                         TRIM(tracer_units)//'/s',  &
                         missing_value=missing_value)

      end do

      call monin_obukhov_init
!---------------------------------------------------------------------
!    mark the module as initialized.
!---------------------------------------------------------------------
      module_is_initialized = .true.

!-----------------------------------------------------------------------

 end subroutine physics_driver_init


!######################################################################
! <SUBROUTINE NAME="physics_driver_down_time_vary">
!  <OVERVIEW>
!    physics_driver_time_vary makes sure that all time-dependent, spacially-
!    independent calculations are completed before entering window or thread
!    loops. Resultant fields are usually saved as module variables in the
!    module where needed.
!  </OVERVIEW>
!  <DESCRIPTION>
!    physics_driver_time_vary makes sure that all time-dependent, spacially-
!    independent calculations are completed before entering window or thread
!    loops. Resultant fields are usually saved as module variables in the
!    module where needed.
!  </DESCRIPTION>
!  <TEMPLATE>
!   call physics_driver_down_time_vary (Time, Time_next)
!
!  </TEMPLATE>
!  <IN NAME="Time" TYPE="time_type">
!   current time
!  </IN>
!  <IN NAME="Time_next" TYPE="time_type">
!   time of next time step
!  </IN>
! </SUBROUTINE>
!

subroutine physics_driver_down_time_vary (Time, Time_next, dt)

!---------------------------------------------------------------------
!    physics_driver_down_time_vary makes sure that all time-dependent, 
!    spacially-independent calculations are completed before entering window 
!    or thread loops. Resultant fields are usually saved as module variables in 
!    the module where needed.
!-----------------------------------------------------------------------

type(time_type),         intent(in)             :: Time, Time_next
real,                    intent(in)             :: dt

type(time_type) :: Time_last
!---------------------------------------------------------------------      
!------------------------------------------------------------------------
!    call damping_driver_time_vary to update the counter determining when
!    the convective drag module will be again called.
!------------------------------------------------------------------------
      call damping_driver_time_vary (dt)

!------------------------------------------------------------------------
!    call atmos_tracer_driver_time_vary to obtain values from the tracer 
!    climatology and emission data at the appropriate time, if needed. 
!------------------------------------------------------------------------
      call atmos_tracer_driver_time_vary (Time)

!----------------------------------------------------------------------
!    call aerosol_time_vary to retrieve appropriate aerosol fields from
!    the climatology, if that source of aerosol is being used.
!----------------------------------------------------------------------
      ! ZNT 06/23/2023: do aerosol in physics down
      if (aer_in_phys_dn) then 
        call aerosol_time_vary (Time, Aerosol_cld)
      endif

!-------------------------------------------------------------------------      

end subroutine physics_driver_down_time_vary



!######################################################################

subroutine physics_driver_down_endts(is,js)

integer, intent(in)  :: is,js

!-----------------------------------------------------------------------
!    call the component xxx_endts routines to perform needed updates to
!    primarily flag and counter variables at the end of the time step, 
!    after all spacial-dependent calculations are completed.
!-----------------------------------------------------------------------
      call damping_driver_endts
      call atmos_tracer_driver_endts

!--------------------------------------------------------------------
!    set a flag to indicate that this check was done and need not be
!    done again.
!--------------------------------------------------------------------
      do_check_args = .false.

end subroutine physics_driver_down_endts

!###################################################################

subroutine physics_driver_up_time_vary (Time, Time_next, dt, &
                                        step_to_call_cosp_in)

!---------------------------------------------------------------------
!    physics_driver_up_time_vary makes sure that all time-dependent, 
!    spacially-independent calculations are completed before entering 
!    window or thread loops. Resultant fields are usually saved as 
!    module variables in the module where needed.
!-----------------------------------------------------------------------

type(time_type),         intent(in)             :: Time
type(time_type),         intent(in)             :: Time_next
real,                    intent(in)             :: dt
logical,                 intent(in)             :: step_to_call_cosp_in

   
!----------------------------------------------------------------------
!    save the flag indicating if this is step to call cosp.
!----------------------------------------------------------------------
    step_to_call_cosp = step_to_call_cosp_in

!----------------------------------------------------------------------
!    call aerosol_time_vary to retrieve appropriate aerosol fields from
!    the climatology, if that source of aerosol is being used.
!----------------------------------------------------------------------
    ! ZNT 06/23/2023: unconditionally call aerosol, and do aerosol in physics down
    if (.not. aer_in_phys_dn) then
      call aerosol_time_vary (Time, Aerosol_cld)
    endif

    if (do_moist_processes) then
!----------------------------------------------------------------------
!    call moist_processes_time_vary to pass needed time-dependent fields 
!    to subordinate modules.
!----------------------------------------------------------------------
      call moist_processes_time_vary (Time_next, dt, i_shallow)
    endif
!----------------------------------------------------------------------
!    call cosp_driver_time_vary to obtain satellite location at current
!    time if orbital data is being collected.
!----------------------------------------------------------------------
    if (do_cosp) call cosp_driver_time_vary (Time_next)

!----------------------------------------------------------------------      

end subroutine physics_driver_up_time_vary


!######################################################################

subroutine physics_driver_up_endts 

!-----------------------------------------------------------------------
!    call the component xxx_endts routines to perform needed updates to
!    primarily flag and counter variables at the end of the time step, 
!    after all spacial-dependent calculations are completed.
!-----------------------------------------------------------------------
    if (do_cosp) call cosp_driver_endts
    
    if (do_moist_processes) then
      call moist_processes_endts 
    endif
    
    ! ZNT 06/23/2023: unconditionally call aerosol, and do aerosol in physics down
    !     06/27/2023: always deallocate after the upward sweep
    call aerosol_endts (Aerosol_cld)

end subroutine physics_driver_up_endts


!######################################################################


!######################################################################
subroutine physics_driver_down (is, ie, js, je, npz,              &
                                Time_prev, Time, Time_next,       &
                                lat, lon, area,                   &
                                Physics_input_block,              &
                                frac_land, rough_mom,             &
                                frac_open_sea,                    &
                                albedo,                           &
                                t_surf_rad, u_ref, v_ref,         & !bqx+ u_ref, v_ref
                                t_ref, q_ref,                     &
                                u_star,    b_star, q_star,        &
                                dtau_du, dtau_dv,  tau_x,  tau_y, &
                                Physics_tendency_block,           &
                                Surf_diff,                        &
                                gust,                             &
                                Rad_flux_control,                 &
                                Rad_flux_block,                   &
                                shflx, lhflx,                     & ! optional input, required by NCEP TKE-based EDMF
                                wind, thv_atm, thv_surf,          & ! optional input, required by NCEP TKE-based EDMF
                                diffm, difft  )

!---------------------------------------------------------------------
!    physics_driver_down calculates "first pass" physics tendencies,
!    associated with radiation, damping and turbulence, and obtains
!    the vertical diffusion tendencies to be passed to the surface and
!    used in the semi-implicit vertical diffusion calculation.
!-----------------------------------------------------------------------

integer,                 intent(in)             :: is, ie, js, je, npz
type(time_type),         intent(in)             :: Time_prev, Time, Time_next
real,dimension(:,:),     intent(in)             :: lat, lon, area
type(physics_input_block_type), intent(in)      :: Physics_input_block
real,dimension(:,:),     intent(in)             :: frac_land,   &
                                                   rough_mom, &
                                                   albedo, t_surf_rad, &
                                                   u_ref, v_ref, & !bqx+
                                                   t_ref, q_ref, &  ! cjg: PBL depth mods
                                                   u_star, b_star, &
                                                   q_star, dtau_du,   &
                                                   dtau_dv, frac_open_sea
real,dimension(:,:),     intent(inout)          :: tau_x,  tau_y
type(physics_tendency_block_type), intent(inout):: Physics_tendency_block
real,dimension(:,:),     intent(out)            :: gust
type(surf_diff_type),    intent(inout)          :: Surf_diff
type(radiation_flux_control_type),  intent(in)  :: Rad_flux_control
type(radiation_flux_block_type),    intent(in)  :: Rad_flux_block
real,  dimension(:,:),   intent(in), optional   :: shflx, lhflx            ! optional input, required by NCEP TKE-based EDMF
real,  dimension(:,:),   intent(in), optional   :: wind, thv_atm, thv_surf ! optional input, required by NCEP TKE-based EDMF
real,  dimension(:,:,:), intent(out)  ,optional :: diffm, difft 

!-----------------------------------------------------------------------
!   intent(in) variables:
!
!      is,ie,js,je    starting/ending subdomain i,j indices of data in 
!                     the physics_window being integrated
!      npz            number of model levels
!      Time_prev      previous time, for variables um,vm,tm,qm,rm 
!                     (time_type)
!      Time           current time, for variables u,v,t,q,r  (time_type)
!      Time_next      next time, used for diagnostics   (time_type)
!      lat            latitude of model points [ radians ]
!      lon            longitude of model points [ radians ]
!      area           grid box area - currently not used [ m**2 ]
!      Physics_input_block  derived type variable containing: 
!         1) p_half         pressure at half levels (offset from t,q,u,v,r)
!                          [ Pa ]
!         2) p_full         pressure at full levels [ Pa }
!         3) z_half         height at half levels [ m ]
!         4) z_full         height at full levels [ m ]
!         5) u              zonal wind at current time step [ m / s ]
!         6) v              meridional wind at current time step [ m / s ]
!         7) t              temperature at current time step [ deg k ]
!         9) q              multiple 3d tracer fields at current time step
!        10) um,vm          zonal and meridional wind at previous time step
!        11) tm             temperature at previous time step
!        12) qm             multiple 3d tracer fields at previous time step
!      Rad_flux_control
!      Rad_flux_block
!      frac_land      fraction of land coverage in a model grid point
!      rough_mom      boundary layer roughness
!      frac_open_sea
!      albedo         surface albedo
!      t_surf_rad     surface radiative temperature
!      t_ref 
!      q_ref
!      u_ref   10m zonal wind
!      v_ref   10m meridional wind 
!      u_star  boundary layer wind speed (frictional speed) 
!      b_star  
!      q_star  boundary layer specific humidity
!      dtau_du derivative of zonal surface stress w.r.t zonal wind speed
!      dtau_dv derivative of meridional surface stress w.r.t meridional wind speed 
!
!  intent(inout) variables:
!
!      tau_x
!      tau_y
!      Physics_tendency_block derived type variable containing:
!          1) u_dt           zonal wind tendency [ m / s**2 ]
!          2) v_dt           meridional wind tendency [ m / s**2 ]
!          3) t_dt           temperature tendency [ deg k / sec ]
!          4) q_dt           multiple tracer tendencies 
!                            (index 1 = specific humidity) 
!                            [ unit / unit / sec ]
!          5) qdiag          multiple 3d diagnostic tracer fields 
!                            [ unit / unit ]
!      Surf_diff      surface_diffusion_type variable
!
!   intent(out) variables:
!
!      gust
!
!   intent(out), optional variables:
!
!      diffm
!      difft
!
!-----------------------------------------------------------------------

!---------------------------------------------------------------------
!    local variables:

      real, dimension(ie-is+1,je-js+1,npz)   :: diff_t_vert, diff_m_vert
      real, dimension(ie-is+1,je-js+1,npz)   :: tdt_rad, tdt_lw
      real, dimension(ie-is+1,je-js+1,npz+1) :: pflux, lphalf   ! ZNT 06/23/2023: pflux
      real, dimension(ie-is+1,je-js+1)       :: z_pbl
      integer                                :: sec, day, n, nextinct
      real                                   :: dt, alpha, dt2
      logical                                :: used

!---------------------------------------------------------------------
!   local variables:
!
!      diff_t_vert     vertical diffusion coefficient for temperature
!                      calculated on the current step
!      diff_m_vert     vertical diffusion coefficient for momentum   
!                      calculated on the current step
!      z_pbl           height of planetary boundary layer
!      sec, day        second and day components of the time_type 
!                      variable
!      dt              model physics time step [ seconds ]
!      alpha           ratio of physics time step to diffusion-smoothing
!                      time scale
!
!---------------------------------------------------------------------
      real, dimension(:,:,:,:), pointer :: r, rm
      real, dimension(:,:,:), pointer :: p_full, p_half, z_full, z_half
      real, dimension(:,:,:), pointer :: udt, vdt, tdt
      real, dimension(:,:,:,:), pointer :: rdt, rdiag
      real, dimension(:,:,:), pointer :: u, v, t, um, vm, tm 

      integer :: kk, k

      type(aerosol_type)               :: Aerosol   ! ZNT 06/23/2023

!---------------------------------------------------------------------
!    set up local pointers into the physics input and physics tendency
!    blocks.
!---------------------------------------------------------------------
      u => Physics_input_block%u
      v => Physics_input_block%v
      t => Physics_input_block%t
      r => Physics_input_block%q
      if (associated(Physics_input_block%um)) then
        um => Physics_input_block%um
        vm => Physics_input_block%vm
        tm => Physics_input_block%tm
        rm => Physics_input_block%qm
      else
        um => Physics_input_block%u
        vm => Physics_input_block%v
        tm => Physics_input_block%t
        rm => Physics_input_block%q
      endif
      p_full => Physics_input_block%p_full
      p_half => Physics_input_block%p_half
      z_full => Physics_input_block%z_full
      z_half => Physics_input_block%z_half
      udt => Physics_tendency_block%u_dt
      vdt => Physics_tendency_block%v_dt
      tdt => Physics_tendency_block%t_dt
      rdt => Physics_tendency_block%q_dt
      rdiag => Physics_tendency_block%qdiag

!---------------------------------------------------------------------
!    verify that the module is initialized.
!---------------------------------------------------------------------
      if ( .not. module_is_initialized) then
        call error_mesg ('physics_driver_mod',  &
                         'module has not been initialized', FATAL)
      endif

!---------------------------------------------------------------------
!    check the size of the input arguments. this is only done on the
!    first call to physics_driver_down.
!---------------------------------------------------------------------
      if (do_check_args) call check_args  &
                   (lat, lon, area, p_half, p_full, z_half, z_full, &
                    u, v, t, r(:,:,:,1), r, um, vm, tm, r(:,:,:,1), rm, &
                    udt, vdt, tdt, rdt(:,:,:,1), rdt, rdiag)

!---------------------------------------------------------------------
!    compute the physics time step (from tau-1 to tau+1).
!---------------------------------------------------------------------
      call get_time (Time_next - Time_prev, sec, day)
      dt = real(sec + day*86400)

!---------------------------------------------------------------------
!     save cmip diagnostics
!---------------------------------------------------------------------
      if (query_cmip_diag_id(ID_pfull)) then
         used = send_cmip_data_3d (ID_pfull, p_full, Time_next, is, js, 1)
      endif

      if (query_cmip_diag_id(ID_phalf)) then
         used = send_cmip_data_3d (ID_phalf, p_half, Time_next, is, js, 1)
      endif

     ! only allow output of geop height on model levels
     ! output from dycore for pressure levels
      if (query_cmip_diag_id(ID_zg)) then
         used = send_cmip_data_3d (ID_zg, z_full, Time_next, is, js, 1)
      endif

!---------------------------------------------------------------------

      if (do_radiation) then
        radturbten(is:ie,js:je,:) = radturbten(is:ie,js:je,:) + Rad_flux_block%tdt_rad(:,:,:)
        surf_diff%tdt_rad(is:ie,js:je,:)=Rad_flux_block%tdt_rad(:,:,:) !miz
      endif
      
#ifdef SCM
! Option to add SCM radiative tendencies from forcing to Rad_flux_block%tdt_lw
! and radturbten

      if (use_scm_rad) then
        call add_scm_tdtlw( Rad_flux_block%tdt_lw )
        call add_scm_tdtlw( radturbten (is:ie,js:je,:) )
        call add_scm_tdtsw( radturbten (is:ie,js:je,:) )
      endif

#endif

!----------------------------------------------------------------------
!    call damping_driver to calculate the various model dampings that
!    are desired. 
!----------------------------------------------------------------------
      z_pbl(:,:) = pbltop(is:ie,js:je) 
      call mpp_clock_begin ( damping_clock )
!!debug
!      if (fms_mpp_pe()==fms_mpp_root_pe()) then
!          print*, 'maximum condensation in physics driver: ',maxval(dqcdt_gw)
!          print*, 'rdt before calling damping_driver: ',maxval(rdt(:,:,:,1))
!      endif
      call damping_driver (is, js, lat, Time_next, dt, area,        &
                           p_full, p_half, z_full, z_half,          &
                           um, vm, tm, rm(:,:,:,1), rm(:,:,:,1:ntp),&
                           u_ref, v_ref, z_pbl,                     &  !bqx+
                           udt, vdt, tdt, dqcdt_gw(is:ie, js:je, :), &
                           tten_gw(is:ie, js:je, :), cqa_gw(is:ie, js:je, :)) !p1l
     call mpp_clock_end ( damping_clock )

!---------------------------------------------------------------------
!    call vert_turb_driver to calculate diffusion coefficients. save
!    the planetary boundary layer height on return.
!---------------------------------------------------------------------

      if (id_tdt_phys_turb > 0) then
        used = send_data ( id_tdt_phys_turb, -2.0*tdt(:,:,:), &
                           Time_next, is, js, 1)
      endif

      do n=1,ntp
        if (id_tracer_phys_turb(n) > 0) then
          used = send_data ( id_tracer_phys_turb(n), -2.0*rdt(:,:,:,n), &
                             Time_next, is, js, 1)
        endif
      end do

      call mpp_clock_begin ( turb_clock )
!-----------------------------------------------------------------------------
!    call aerosol driver to obtain aerosol data for condensation calculations. 
!-----------------------------------------------------------------------------
      if (aer_in_phys_dn) then
          pflux(:,:,1) = 0.0E+00
          do k=2,size(p_full,3)
            pflux(:,:,k) = 0.5E+00*(p_full(:,:,k-1) + p_full(:,:,k))
          end do
          pflux(:,:,size(p_full,3)+1) = p_full(:,:,size(p_full,3))

          ! ZNT 06/27/2023: aerosol is computed from tracer amounts before the tracer driver
          call aerosol_driver (is, js, Time, r, p_half, pflux, &
                             Aerosol_cld,Aerosol, override_aerosols_cloud)

          call vert_turb_driver (is, js, Time, Time_next, dt,            &
                             Rad_flux_block%tdt_lw, frac_land,  &
                             p_half, p_full, z_half, z_full,         &
                             u_ref, v_ref, tau_x, tau_y, t_surf_rad, &  ! ZNT
                             shflx, lhflx,                           &  ! ZNT
                             wind, thv_atm, thv_surf,                &  ! ZNT
                             t_ref, q_ref,                           &  ! cjg: PBL depth mods
                             u_star, b_star, q_star, rough_mom,      &
                             lat, convect(is:ie,js:je),              &
                             u, v, t, r(:,:,:,1), r, um, vm,                  &
                             tm, rm(:,:,:,1), rm, rdiag,                      &
                             udt, vdt, tdt, rdt(:,:,:,1), rdt,                &
                             diff_t_vert, diff_m_vert, gust, z_pbl,           &
                             Aerosol = Aerosol)
           call aerosol_dealloc (Aerosol)
      else
           call vert_turb_driver (is, js, Time, Time_next, dt,            &
                             Rad_flux_block%tdt_lw, frac_land,  &
                             p_half, p_full, z_half, z_full,         &
                             u_ref, v_ref, tau_x, tau_y, t_surf_rad, &  ! ZNT
                             shflx, lhflx,                           &  ! ZNT
                             wind, thv_atm, thv_surf,                &  ! ZNT
                             t_ref, q_ref,                           &  ! cjg: PBL depth mods
                             u_star, b_star, q_star, rough_mom,      &
                             lat, convect(is:ie,js:je),              &
                             u, v, t, r(:,:,:,1), r, um, vm,                  &
                             tm, rm(:,:,:,1), rm, rdiag,                      &
                             udt, vdt, tdt, rdt(:,:,:,1), rdt,                &
                             diff_t_vert, diff_m_vert, gust, z_pbl)
      endif ! aer_in_phys_dn
     call mpp_clock_end ( turb_clock )
     pbltop(is:ie,js:je) = z_pbl(:,:)

      if (id_tdt_phys_turb > 0) then
        used = send_data ( id_tdt_phys_turb, +2.0*tdt(:,:,:), &
                           Time_next, is, js, 1)
      endif

      do n=1,ntp
        if (id_tracer_phys_turb(n) > 0) then
          used = send_data ( id_tracer_phys_turb(n), +2.0*rdt(:,:,:,n), &
                             Time_next, is, js, 1)
        endif
      end do

!-----------------------------------------------------------------------
!    process any tracer fields.
!-----------------------------------------------------------------------

      nextinct = get_tracer_index(MODEL_ATMOS,'Extinction')
      if (Rad_flux_control%do_rad .and. nextinct /= NO_TRACER) then
        rdiag(:,:,:,nextinct) = Rad_flux_block%extinction(:,:,:)
      endif

      call mpp_clock_begin ( tracer_clock )
      call atmos_tracer_driver (is, ie, js, je, Time, lon, lat,  &
                                area, z_pbl, rough_mom,         &
                                frac_open_sea, frac_land, &
                                p_half, p_full,  &
                                u, v, t, r(:,:,:,1), r, rm, rdt, rdiag, dt, &
                                u_star, b_star, q_star, z_half, z_full, &
                                t_surf_rad, albedo, Time_next, &
                                Rad_flux_block%flux_sw_down_vis_dir, &
                                Rad_flux_block%flux_sw_down_vis_dif)
      call mpp_clock_end ( tracer_clock )

!-----------------------------------------------------------------------
!    optionally use an implicit calculation of the vertical diffusion 
!    coefficients.
!
!    the vertical diffusion coefficients are solved using an implicit
!    solution to the following equation:
!
!    dK/dt   = - ( K - K_cur) / tau_diff
!
!    where K         = diffusion coefficient
!          K_cur     = diffusion coefficient diagnosed from current 
!                      time steps' state
!          tau_diff  = time scale for adjustment
!
!    in the code below alpha = dt / tau_diff
!---------------------------------------------------------------------
      if (diffusion_smooth) then
        call get_time (Time_next - Time, sec, day)
        dt2 = real(sec + day*86400)
        alpha = dt2/tau_diff
        diff_m(is:ie,js:je,:) = (diff_m(is:ie,js:je,:) +       &
                                 alpha*(diff_m_vert(:,:,:) +  &
                                 diff_cu_mo(is:ie,js:je,:)) )/&
                                 (1. + alpha)
        where (diff_m(is:ie,js:je,:) < diff_min)
          diff_m(is:ie,js:je,:) = 0.0
        end where
        diff_t(is:ie,js:je,:) = (diff_t(is:ie,js:je,:) +      &
                                 alpha*diff_t_vert(:,:,:) )/  &
                                 (1. + alpha)
        where (diff_t(is:ie,js:je,:) < diff_min)
          diff_t(is:ie,js:je,:) = 0.0
        end where
      else
        diff_t(is:ie,js:je,:) = diff_t_vert
        diff_m(is:ie,js:je,:) = diff_m_vert + diff_cu_mo(is:ie, js:je,:)
      end if

!-----------------------------------------------------------------------
!    call vert_diff_driver_down to calculate the first pass atmos-
!    pheric vertical diffusion.
!-----------------------------------------------------------------------

      if (id_tdt_phys_vdif_dn > 0) then
        used = send_data ( id_tdt_phys_vdif_dn, -2.0*tdt(:,:,:), &
                           Time_next, is, js, 1)
      endif

      do n=1,ntp
        if (id_tracer_phys_vdif_dn(n) > 0) then
          used = send_data ( id_tracer_phys_vdif_dn(n), -2.0*rdt(:,:,:,n), &
                             Time_next, is, js, 1)
        endif
      end do

      call mpp_clock_begin ( diff_down_clock )
      radturbten(is:ie,js:je,:) = radturbten(is:ie,js:je,:) - tdt(:,:,:)
        qturbten(is:ie,js:je,:) =   qturbten(is:ie,js:je,:) - rdt(:,:,:,1)   ! ZNT
      call vert_diff_driver_down (is, js, Time_next, dt, p_half,   &
          p_full, z_full,   &
          diff_m(is:ie,js:je,:),         &
          diff_t(is:ie,js:je,:),         &
          u ,v ,t ,r(:,:,:,1) ,r(:,:,:,1:ntp), &
          dtau_du, dtau_dv, tau_x, tau_y,  &
          udt, vdt, tdt, rdt(:,:,:,1), rdt,        &
          Surf_diff)

      if (id_tdt_phys_vdif_dn > 0) then
        used = send_data ( id_tdt_phys_vdif_dn, +2.0*tdt(:,:,:), &
                           Time_next, is, js, 1)
      endif

      do n=1,ntp
        if (id_tracer_phys_vdif_dn(n) > 0) then
          used = send_data ( id_tracer_phys_vdif_dn(n), +2.0*rdt(:,:,:,n), &
                             Time_next, is, js, 1)
        endif
      end do

!---------------------------------------------------------------------
!    if desired, return diff_m and diff_t to calling routine.
!-----------------------------------------------------------------------
      if (present(difft)) then
        difft = diff_t(is:ie,js:je,:)
      endif
      if (present(diffm)) then
        diffm = diff_m(is:ie,js:je,:)
      endif

     call mpp_clock_end ( diff_down_clock )

      u => null()
      v => null()
      t => null()
      r => null()
      um => null()
      vm => null()
      tm => null()
      rm => null()
      p_full => null()
      p_half => null()
      z_full => null()
      z_half => null()
      udt => null()
      vdt => null()
      tdt => null()
      rdt => null()
      rdiag => null()

 end subroutine physics_driver_down



!#######################################################################
 subroutine physics_driver_up (is, ie, js, je, npz,        &
                               Time_prev, Time, Time_next, &
                               lat, lon, area,             &
                               Physics_input_block,        &
                               frac_land,                  &
                               u_star, b_star, q_star,     &
                               shflx, lhflx,               &!miz
                               Physics_tendency_block,     &
                               Moist_clouds_block,         &
                               Cosp_block, Surf_diff,      &
                               lprec, fprec, gust)

!----------------------------------------------------------------------
!    physics_driver_up completes the calculation of vertical diffusion 
!    and also handles moist physical processes.
!---------------------------------------------------------------------

integer,                intent(in)                :: is, ie, js, je, npz
type(time_type),        intent(in)                :: Time_prev, Time, Time_next
real,dimension(:,:),    intent(in)                :: lat, lon, area
type(physics_input_block_type), intent(inout)     :: Physics_input_block
real,dimension(:,:),    intent(in)                :: frac_land
real,dimension(:,:),    intent(in)                :: u_star, b_star, q_star, shflx, lhflx!miz
type(physics_tendency_block_type), intent(inout)  :: Physics_tendency_block
type(clouds_from_moist_block_type), intent(inout) :: Moist_clouds_block
type(cosp_from_rad_block_type),     intent(inout) :: Cosp_block
type(surf_diff_type),   intent(inout)             :: Surf_diff
real,dimension(:,:),    intent(out)               :: lprec, fprec
real,dimension(:,:),    intent(inout)             :: gust

!-----------------------------------------------------------------------
!   intent(in) variables:
!
!      is,ie,js,je    starting/ending subdomain i,j indices of data in 
!                     the physics_window being integrated
!      npz            number of vertical levels
!      Time_prev      previous time, for variables um,vm,tm,qm,rm 
!                     (time_type)
!      Time           current time, for variables u,v,t,q,r  (time_type)
!      Time_next      next time, used for diagnostics   (time_type)
!      lat            latitude of model points [ radians ]
!      lon            longitude of model points [ radians ]
!      area           grid box area - currently not used [ m**2 ]
!      frac_land      fraction of land coverage in a model grid point
!      u_star
!      b_star
!      q_star
!
!  intent(inout) variables:
!
!      Physics_input_block  derived type variable containing: 
!         1) p_half         pressure at half levels (offset from t,q,u,v,r)
!                          [ Pa ]
!         2) p_full         pressure at full levels [ Pa }
!         3) z_half         height at half levels [ m ]
!         4) z_full         height at full levels [ m ]
!         5) u              zonal wind at current time step [ m / s ]
!         6) v              meridional wind at current time step [ m / s ]
!         7) t              temperature at current time step [ deg k ]
!         9) q              multiple 3d tracer fields at current time step
!        10) um,vm          zonal and meridional wind at previous time step
!        11) tm             temperature at previous time step
!        12) qm             multiple 3d tracer fields at previous time step
!        13) omega          Veritical pressure tendency 
!      Physics_tendency_block derived type variable containing:
!          1) u_dt           zonal wind tendency [ m / s**2 ]
!          2) v_dt           meridional wind tendency [ m / s**2 ]
!          3) t_dt           temperature tendency [ deg k / sec ]
!          4) q_dt           multiple tracer tendencies 
!                            (index 1 = specific humidity) 
!                            [ unit / unit / sec ]
!          5) qdiag          multiple 3d diagnostic tracer fields 
!                            [ unit / unit ]
!      Moist_clouds_block
!      Cosp_block
!      Surf_diff      surface_diffusion_type variable
!      gust
!
!   intent(out) variables:
!
!      lprec     
!      fprec       
!
!   intent(in), optional variables:
!
!--------------------------------------------------------------------
 
!--------------------------------------------------------------------
!   local variables:

      type (precip_flux_type)          :: MP2cosp
      type (Phys2cosp_type)            :: Phys2cosp
      type (phys_mp_exch_type)         :: Phys_mp_exch
      type(aerosol_type)               :: Aerosol
      real, dimension(ie-is+1, je-js+1)          :: gust_cv
      real, dimension(ie-is+1, je-js+1, npz+1)   :: pflux, lphalf
      real, dimension(ie-is+1, je-js+1), target  :: tdt_shf,  qdt_lhf
      integer :: sec, day
      real    :: dt
      integer :: i, j , k, n
      real    :: alphb
      integer :: imax, jmax, kmax
      logical :: used

      type(MP_removal_type) :: Removal_mp
   
!---------------------------------------------------------------------
!   local variables:
!
!        MP2cosp
!        Phys2cosp
!        Phys_mp_exch
!        Aerosol          aerosol_type variable describing the aerosol
!                         fields to be seen by the moist processes routines
!        gust_cv
!        pflux
!        tdt_shf          temperature tendency from sensible heat flux
!        qdt_lhf          moisture tendency from latent heat flux
!        sec, day         second and day components of the time_type 
!                         variable
!        dt               physics time step [ seconds ]
!        i,j,k,n
!        alphb
!        imax,jmax,kmax
!        used
!
!---------------------------------------------------------------------

!----------------------------------------------------------------------
!   local pointers to derived type components
!----------------------------------------------------------------------
      real, dimension(:,:,:),   pointer :: t                    
      real, dimension(:,:,:,:), pointer :: r
      real, dimension(:,:,:),   pointer :: p_full, p_half                 
      real, dimension(:,:,:),   pointer :: tdt
      real, dimension(:,:,:,:), pointer :: rdt          

      t => Physics_input_block%t
      r => Physics_input_block%q
      p_full => Physics_input_block%p_full
      p_half => Physics_input_block%p_half
      tdt => Physics_tendency_block%t_dt
      rdt => Physics_tendency_block%q_dt


!---------------------------------------------------------------------
!    verify that the module is initialized.
!---------------------------------------------------------------------
      if ( .not. module_is_initialized) then
        call error_mesg ('physics_driver_mod',  &
             'module has not been initialized', FATAL)
      endif

!----------------------------------------------------------------------
!    define model spatial dimensions.
!----------------------------------------------------------------------
      imax = ie -is + 1
      jmax = je- js + 1
      kmax = npz        

!-------------------------------------------------------------------------
!    if cosp is activated and this is a step on which cosp input data is
!    to be collected, set up pointers or allocate the necessary derived 
!    type variable components.    
!-------------------------------------------------------------------------
      if (do_cosp) then
        if (step_to_call_cosp) then
          allocate (MP2cosp%fl_lsrain(imax, jmax, kmax))
          allocate (MP2cosp%fl_lssnow(imax, jmax, kmax))
          allocate (MP2cosp%fl_lsgrpl(imax, jmax, kmax))
          allocate (MP2cosp%fl_ccrain(imax, jmax, kmax))
          allocate (MP2cosp%fl_ccsnow(imax, jmax, kmax))

          allocate (Phys2cosp%temp_last(imax, jmax, kmax))
          allocate (Phys2cosp%q_last(imax, jmax, kmax))
          Phys2cosp%p_full => Physics_input_block%p_full
          Phys2cosp%z_full => Physics_input_block%z_full
          Phys2cosp%p_half => Physics_input_block%p_half
          Phys2cosp%z_half => Physics_input_block%z_half
          Phys2cosp%u  => Physics_input_block%u
          Phys2cosp%v  => Physics_input_block%v
          allocate (Phys2cosp%frac_land(imax, jmax      ))
          allocate (Phys2cosp%lat   (imax, jmax      ))
         allocate (Phys2cosp%lon   (imax, jmax      ))
        endif
      endif

!----------------------------------------------------------------------
!    compute the physics time step (from tau-1 to tau+1).
!---------------------------------------------------------------------
      call get_time (Time_next-Time_prev, sec, day)
      dt = real(sec+day*86400)

!-------------------------------------------------------------------------
!    save temp and moisture tendencies before calculating vertical
!    diffusion. 
!-------------------------------------------------------------------------
      if (id_tdt_phys_vdif_up > 0) then
        used = send_data ( id_tdt_phys_vdif_up, -2.0*tdt(:,:,:), &
                           Time_next, is, js, 1)
      endif

      do n=1,ntp
        if (id_tracer_phys_vdif_up(n) > 0) then
          used = send_data ( id_tracer_phys_vdif_up(n), -2.0*rdt(:,:,:,n), &
                             Time_next, is, js, 1)
        endif
      end do

      call mpp_clock_begin ( diff_up_clock )
!------------------------------------------------------------------
!    call vert_diff_driver_up to complete the vertical diffusion
!    calculation.
!------------------------------------------------------------------
      call vert_diff_driver_up (is, js, Time_next, dt, p_half,   &
                                Surf_diff, tdt, rdt(:,:,:,1), rdt )

!-----------------------------------------------------------------------
!    add the temperature tendency due to vertical  diffusion to radturbten.
!-----------------------------------------------------------------------
      radturbten(is:ie,js:je,:) = radturbten(is:ie,js:je,:) + tdt(:,:,:)
        qturbten(is:ie,js:je,:) =   qturbten(is:ie,js:je,:) + rdt(:,:,:,1)   ! ZNT
      call mpp_clock_end ( diff_up_clock )

!-------------------------------------------------------------------------
!    complete calculation of vertical diffusion tendency diagnostics.
!-------------------------------------------------------------------------
      if (id_tdt_phys_vdif_up > 0) then
        used = send_data ( id_tdt_phys_vdif_up, +2.0*tdt(:,:,:), &
                           Time_next, is, js, 1)
      endif

      do n=1,ntp
        if (id_tracer_phys_vdif_up(n) > 0) then
          used = send_data ( id_tracer_phys_vdif_up(n), +2.0*rdt(:,:,:,n), &
                             Time_next, is, js, 1)
        endif
      end do

!-----------------------------------------------------------------------
!    prepare to call moist_processes, which calculates moist physics terms,
!    including convection and processes involving condensation, if 
!    desired.
!-----------------------------------------------------------------------
      if (do_moist_processes) then

!-----------------------------------------------------------------------
!    set up diagnostics to capture tendencies due to moist_processes.
!-----------------------------------------------------------------------
        if (id_tdt_phys_moist > 0) then
          used = send_data ( id_tdt_phys_moist, -2.0*tdt(:,:,:), &
                             Time_next, is, js, 1)
        endif

        do n=1,ntp
          if (id_tracer_phys_moist(n) > 0) then
            used = send_data ( id_tracer_phys_moist(n), -2.0*rdt(:,:,:,n), &
                               Time_next, is, js, 1)
          endif
        end do

        call mpp_clock_begin ( moist_processes_clock )

!-----------------------------------------------------------------------
!    call aerosol driver to obtain aerosol data needed in condensation 
!    calculations.
!-----------------------------------------------------------------------
      ! ZNT 06/27/2023 
      ! Note: aerosol driver may be called twice (down and up), and here
      !       aerosol is computed from tracer amounts AFTER the tracer driver

        pflux(:,:,1) = 0.0E+00
        do k=2,size(p_full,3)
          pflux(:,:,k) = 0.5E+00*(p_full(:,:,k-1) + p_full(:,:,k))
        end do
        pflux(:,:,size(p_full,3)+1) = p_full(:,:,size(p_full,3))
        call aerosol_driver (is, js, Time, r, p_half, pflux, &
            Aerosol_cld,Aerosol, override_aerosols_cloud)

!------------------------------------------------------------------------
!   set up pointers to the module variables that are transferred between
!   physics_driver and moist_processes.
!------------------------------------------------------------------------
        Phys_mp_exch%diff_t => diff_t(is:ie,js:je,:)
        Phys_mp_exch%radturbten => radturbten(is:ie,js:je,:)
        Phys_mp_exch%qturbten   => qturbten  (is:ie,js:je,:)  ! ZNT
        Phys_mp_exch%cush       => cush      (is:ie,js:je  )
        Phys_mp_exch%cbmf       => cbmf      (is:ie,js:je  )
        Phys_mp_exch%pbltop     => pbltop    (is:ie,js:je  )
        Phys_mp_exch%diff_cu_mo => diff_cu_mo(is:ie,js:je,:)
        Phys_mp_exch%convect    => convect   (is:ie,js:je  )
        Phys_mp_exch%tdt_shf    => tdt_shf 
        Phys_mp_exch%qdt_lhf    => qdt_lhf 
        Phys_mp_exch%hmint      => hmint     (is:ie,js:je  )
        Phys_mp_exch%cgust      => cgust    (is:ie,js:je  )
        Phys_mp_exch%pblhto     => pblhto    (is:ie,js:je  )
        Phys_mp_exch%rkmo       => rkmo      (is:ie,js:je  )
        Phys_mp_exch%taudpo     => taudpo    (is:ie,js:je  )
        Phys_mp_exch%exist_shconv  => exist_shconv (is:ie,js:je,:)
        Phys_mp_exch%exist_dpconv  => exist_dpconv (is:ie,js:je,:)
        Phys_mp_exch%pblht_prev    => pblht_prev   (is:ie,js:je,:)
        Phys_mp_exch%hlsrc_prev    => pblht_prev   (is:ie,js:je,:)
        Phys_mp_exch%qtsrc_prev    => pblht_prev   (is:ie,js:je,:)
        Phys_mp_exch%cape_prev     => pblht_prev   (is:ie,js:je,:)
        Phys_mp_exch%cin_prev      => pblht_prev   (is:ie,js:je,:)
        Phys_mp_exch%tten_gw        => tten_gw (is:ie,js:je,:)  ! P1L
        Phys_mp_exch%cqa_gw       => cqa_gw (is:ie, js:je, :) ! P1L

!-----------------------------------------------------------------------
!    call moist processes to compute moist physics, including convection 
!    and processes involving condenstion.
!-----------------------------------------------------------------------
        call moist_processes (    &
              is, ie, js, je, npz, Time_next,     frac_land, u_star,  &
              b_star, q_star, area, lon, lat, Physics_input_block,   &
              Moist_clouds_block, Physics_tendency_block, Phys_mp_exch, &
              Surf_diff, Removal_mp, shflx, lhflx,  &
              lprec, fprec, gust_cv, Aerosol=Aerosol)
        call mpp_clock_end ( moist_processes_clock )

!-------------------------------------------------------------------------
!    save total condensation to be used in gw_beres
!-------------------------------------------------------------------------
     dqcdt_gw(is:ie,js:je,:)=-rdt(:,:,:,1)

!-------------------------------------------------------------------------
!    save the cumulus momentum output field in a module variable for use
!    in vert diff calculation done in physics_driver_down. reinitialize
!    radturbten for use on next time step.
!-------------------------------------------------------------------------
        radturbten(is:ie,js:je,:) = 0.0
          qturbten(is:ie,js:je,:) = 0.0   ! ZNT

!---------------------------------------------------------------------
!    add the convective gustiness effect to that previously obtained 
!    from non-convective parameterizations.
!---------------------------------------------------------------------
        gust = sqrt( gust*gust + gust_cv*gust_cv)

!------------------------------------------------------------------------
!    complete calculation of moist processes tendency diagnostics.
!------------------------------------------------------------------------
        if (id_tdt_phys_moist > 0) then
          used = send_data ( id_tdt_phys_moist, +2.0*tdt(:,:,:), &
                             Time_next, is, js, 1)
        endif
        do n=1,ntp
          if (id_tracer_phys_moist(n) > 0) then
            used = send_data ( id_tracer_phys_moist(n), +2.0*rdt(:,:,:,n), &
                               Time_next, is, js, 1)
          endif
        end do

!------------------------------------------------------------------------
!    calculate temperature and tracer tendencies due to model physics.
!------------------------------------------------------------------------
        if (id_tdt_phys > 0) then
           used = send_data ( id_tdt_phys, tdt(:,:,:), &
                              Time_next, is, js, 1)
        endif
        do n=1,ntp
          if (id_tracer_phys(n) > 0) then
            used = send_data ( id_tracer_phys(n), rdt(:,:,:,n), &
                               Time_next, is, js, 1)
          endif
        end do

!----------------------------------------------------------------------
!    if the Aerosol derived type variable component arrays were allocated, 
!    call aerosol_dealloc to deallocate them.
!----------------------------------------------------------------------
        call aerosol_dealloc (Aerosol)

      !------ CMIP diagnostics (tendencies due to physics) ------
      if (query_cmip_diag_id(ID_tntmp) .or. query_cmip_diag_id(ID_tnhusmp)) then
        lphalf = log(p_half)
      endif
      if (query_cmip_diag_id(ID_tntmp)) then 
         used = send_cmip_data_3d (ID_tntmp, tdt(:,:,:), Time_next, is, js,1, phalf=lphalf)
      endif
      if (query_cmip_diag_id(ID_tnhusmp)) then
        used = send_cmip_data_3d (ID_tnhusmp, rdt(:,:,:,nsphum), Time_next, is, js, 1, phalf=lphalf)
      endif


!-----------------------------------------------------------------------
!    code needed to execute COSP
!-----------------------------------------------------------------------
        if (do_cosp) then
          call mpp_clock_begin ( cosp_clock )
          alphb = SUM(temp_last(is:ie,js:je,:))

!---------------------------------------------------------------------
!    on the first step of a job segment, the values of t,q and precip 
!    flux will not be available at the proper time level. in this case
!    denoted by temp-_last = 0.0, use values from the current step for 
!    t, q and precip flux.
!---------------------------------------------------------------------
          if (alphb == 0.) then
            call define_cosp_precip_fluxes (is, js, Precip_flux,   &
                                                               Removal_mp)
            if (step_to_call_cosp) then
              Phys2cosp%temp_last(:,:,:) = t(:,:,:) + dt*tdt(:,:,:)
              Phys2cosp%q_last(:,:,:) = r(:,:,:,1) + dt*rdt(:,:,:,1)
              MP2cosp%fl_lsrain(:,:,:) =  &
                                       Precip_flux%fl_lsrain(is:ie,js:je,:)
              MP2cosp%fl_lssnow(:,:,:) =   &
                                       Precip_flux%fl_lssnow(is:ie,js:je,:)
              MP2cosp%fl_lsgrpl(:,:,:) =    &
                                       Precip_flux%fl_lsgrpl(is:ie,js:je,:)
              MP2cosp%fl_ccrain(:,:,:) =    &
                                       Precip_flux%fl_ccrain(is:ie,js:je,:)
              MP2cosp%fl_ccsnow(:,:,:) =     &
                                       Precip_flux%fl_ccsnow(is:ie,js:je,:)
            endif
          else

!--------------------------------------------------------------------
!    on all other steps of the job on which the cosp simulator is 
!    called, define input variables needed by COSP from values computed on
!    the last step that are currently available, before calculating new
!    values for the current step.
!--------------------------------------------------------------------
            if (step_to_call_cosp) then
              MP2cosp%fl_lsrain(:,:,:) =    &
                                    Precip_flux%fl_lsrain(is:ie,js:je,:)
              MP2cosp%fl_lssnow(:,:,:) =    &
                                     Precip_flux%fl_lssnow(is:ie,js:je,:)
              MP2cosp%fl_lsgrpl(:,:,:) =    &
                                     Precip_flux%fl_lsgrpl(is:ie,js:je,:)
              MP2cosp%fl_ccrain(:,:,:) =    &
                                     Precip_flux%fl_ccrain(is:ie,js:je,:)
              MP2cosp%fl_ccsnow(:,:,:) =    &
                                     Precip_flux%fl_ccsnow(is:ie,js:je,:)
              Phys2cosp%temp_last(:,:,:) = temp_last(is:ie,js:je,:)
              Phys2cosp%q_last(:,:,:)    = q_last(is:ie,js:je,:)
            endif
            call define_cosp_precip_fluxes (is, js, Precip_flux,   &
                                                               Removal_mp)
          endif

          if (step_to_call_cosp) then
!----------------------------------------------------------------------
!    define the remaining input fields needed by cosp_driver.
!----------------------------------------------------------------------
            Phys2cosp%lat    = lat*180./ACOS(-1.0)
            Phys2cosp%lon    = lon*180./ACOS(-1.0)
            Phys2cosp%frac_land = frac_land

            call cosp_driver (   &
             is, ie, js, je, Time_next, MP2cosp, Phys2cosp, Cosp_block)

!-----------------------------------------------------------------------
!    deallocate arrays used to hold cosp input fields.
!-----------------------------------------------------------------------
            deallocate (MP2cosp%fl_lsrain)
            deallocate (MP2cosp%fl_lssnow)
            deallocate (MP2cosp%fl_lsgrpl)
            deallocate (MP2cosp%fl_ccrain)
            deallocate (MP2cosp%fl_ccsnow)
 
            deallocate (Phys2cosp%temp_last)
            deallocate (Phys2cosp%q_last)
            Phys2cosp%p_full => null()
            Phys2cosp%z_full => null()
            Phys2cosp%p_half => null()
            Phys2cosp%z_half => null()
            Phys2cosp%u => null()
            Phys2cosp%v => null()
            deallocate (Phys2cosp%frac_land)
            deallocate (Phys2cosp%lat   )
            deallocate (Phys2cosp%lon   )
          endif ! (step_to_call_cosp)

!--------------------------------------------------------------------
!    save t and q from end of step for use with next call to COSP
!--------------------------------------------------------------------
          temp_last(is:ie,js:je,:) = t(:,:,:) + tdt(:,:,:)*dt
          q_last   (is:ie,js:je,:) = r(:,:,:,1) + rdt(:,:,:,1)*dt
          call mpp_clock_end ( cosp_clock )
        else
            if (allocated(Removal_mp%ice_precflxh)) then
                deallocate(Removal_mp%ice_precflxh)
            endif
            if (allocated(Removal_mp%liq_precflxh)) then
                deallocate(Removal_mp%liq_precflxh)
            endif
            if (allocated(Removal_mp%rain3d)) then
                deallocate(Removal_mp%rain3d)
            endif
            if (allocated(Removal_mp%snowclr3d)) then
                deallocate(Removal_mp%snowclr3d)
            endif
        endif ! (do_cosp)
      endif ! do_moist_processes
       
!-----------------------------------------------------------------------
!    nullify all local pointers.
!-----------------------------------------------------------------------

      t => null()
      r => null()
      p_full => null()
      p_half => null()
      tdt => null()
      rdt => null()

      Phys_mp_exch%diff_t => null()
      Phys_mp_exch%radturbten => null()
      Phys_mp_exch%qturbten => null()        ! ZNT
      Phys_mp_exch%diff_cu_mo => null()
      Phys_mp_exch%cush   => null()
      Phys_mp_exch%cbmf   => null()
      Phys_mp_exch%pbltop => null()
      Phys_mp_exch%convect => null()
      Phys_mp_exch%tdt_shf => null()
      Phys_mp_exch%qdt_lhf => null()
      Phys_mp_exch%hmint   => null()
      Phys_mp_exch%cgust   => null()
      Phys_mp_exch%pblhto  => null()
      Phys_mp_exch%rkmo    => null()
      Phys_mp_exch%taudpo  => null()
      Phys_mp_exch%exist_shconv => null()
      Phys_mp_exch%exist_dpconv => null()
      Phys_mp_exch%pblht_prev   => null()
      Phys_mp_exch%hlsrc_prev   => null()
      Phys_mp_exch%qtsrc_prev   => null()
      Phys_mp_exch%cape_prev   => null()
      Phys_mp_exch%cin_prev   => null()
      Phys_mp_exch%tten_gw => null()
      Phys_mp_exch%cqa_gw => null()

!-----------------------------------------------------------------------



 end subroutine physics_driver_up


!#######################################################################
subroutine physics_driver_end (Time, Physics, Moist_clouds,  &
                               Physics_tendency, Atm_block)

!---------------------------------------------------------------------
!    physics_driver_end is the destructor for physics_driver_mod.
!---------------------------------------------------------------------

type(time_type), intent(in) :: Time 
type(physics_type), intent(in) :: Physics
type(clouds_from_moist_type), intent(inout) :: Moist_clouds(:)
type(physics_tendency_type),  intent(inout) :: Physics_tendency
type(block_control_type), intent(in) :: Atm_block

!--------------------------------------------------------------------
!   intent(in) variables:
! 
!      Time      current time [ time_type(days, seconds) ]
!
!--------------------------------------------------------------------
integer :: n, nb, nc, ibs, ibe, jbs, jbe
integer :: moist_processes_term_clock, damping_term_clock, turb_term_clock, &
           diff_term_clock, aerosol_term_clock,  &
           tracer_term_clock, cosp_term_clock

!---------------------------------------------------------------------
!    verify that the module is initialized.
!---------------------------------------------------------------------
      if ( .not. module_is_initialized) then
        call error_mesg ('physics_driver_mod',  &
              'module has not been initialized', FATAL)
      endif

      moist_processes_term_clock =      &
        mpp_clock_id( '   Phys_driver_term: MP: Termination', &
                grain=CLOCK_MODULE_DRIVER )
      damping_term_clock         =     &
        mpp_clock_id( '   Phys_driver_term: Damping: Termination',    &
                  grain=CLOCK_MODULE_DRIVER )
      turb_term_clock            =      &
        mpp_clock_id( '   Phys_driver_term: Vert. Turb.: Termination', &
                  grain=CLOCK_MODULE_DRIVER )
      diff_term_clock       =     &
        mpp_clock_id( '   Phys_driver_term: Vert. Diff.: Termination',   &
                 grain=CLOCK_MODULE_DRIVER )
      cosp_term_clock       =       &
        mpp_clock_id( '   Phys_driver_term: COSP: Termination', &
                       grain=CLOCK_MODULE_DRIVER )
      ! if (do_moist_processes) &   ! ZNT 06/23/2023: unconditionally call aerosol_end
      aerosol_term_clock       =       &
        mpp_clock_id( '   Phys_driver_term: Aerosol: Termination', &
                       grain=CLOCK_MODULE_DRIVER )
      tracer_term_clock          =      &
        mpp_clock_id( '   Phys_driver_term: Tracer: Termination',    &
                 grain=CLOCK_MODULE_DRIVER )

      do nb = 1, Atm_block%nblks
        ibs = Atm_block%ibs(nb)-Atm_block%isc+1
        ibe = Atm_block%ibe(nb)-Atm_block%isc+1
        jbs = Atm_block%jbs(nb)-Atm_block%jsc+1
        jbe = Atm_block%jbe(nb)-Atm_block%jsc+1

        do nc = 1, size(Moist_clouds(1)%block(nb)%Cloud_data,1)

          ! common to all cloud schemes
          Restart%Cloud_data(nc)%cloud_area    (ibs:ibe,jbs:jbe,:) = Moist_clouds(1)%block(nb)%Cloud_data(nc)%cloud_area
          Restart%Cloud_data(nc)%liquid_amt    (ibs:ibe,jbs:jbe,:) = Moist_clouds(1)%block(nb)%Cloud_data(nc)%liquid_amt
          Restart%Cloud_data(nc)%ice_amt       (ibs:ibe,jbs:jbe,:) = Moist_clouds(1)%block(nb)%Cloud_data(nc)%ice_amt
          Restart%Cloud_data(nc)%droplet_number(ibs:ibe,jbs:jbe,:) = Moist_clouds(1)%block(nb)%Cloud_data(nc)%droplet_number
          Restart%Cloud_data(nc)%scheme_name                       = Moist_clouds(1)%block(nb)%Cloud_data(nc)%scheme_name

          ! properties specific to large-scale/stratiform clouds
          if (trim(Moist_clouds(1)%block(nb)%Cloud_data(nc)%scheme_name) .eq. 'strat_cloud') then
            Restart%Cloud_data(nc)%ice_number(ibs:ibe,jbs:jbe,:) = Moist_clouds(1)%block(nb)%Cloud_data(nc)%ice_number
            Restart%Cloud_data(nc)%rain      (ibs:ibe,jbs:jbe,:) = Moist_clouds(1)%block(nb)%Cloud_data(nc)%rain
            Restart%Cloud_data(nc)%snow      (ibs:ibe,jbs:jbe,:) = Moist_clouds(1)%block(nb)%Cloud_data(nc)%snow
            Restart%Cloud_data(nc)%rain_size (ibs:ibe,jbs:jbe,:) = Moist_clouds(1)%block(nb)%Cloud_data(nc)%rain_size
            Restart%Cloud_data(nc)%snow_size (ibs:ibe,jbs:jbe,:) = Moist_clouds(1)%block(nb)%Cloud_data(nc)%snow_size
          endif
 
          ! properties specific to uw shallow convective clouds
          if (trim(Moist_clouds(1)%block(nb)%Cloud_data(nc)%scheme_name) .eq. 'uw_conv') then
            Restart%Cloud_data(nc)%ice_number(ibs:ibe,jbs:jbe,:) = Moist_clouds(1)%block(nb)%Cloud_data(nc)%ice_number
          endif

        enddo
      enddo

!----------------------------------------------------------------------
!    call physics_driver_netcdf to process physics_driver restart file.
!----------------------------------------------------------------------
      call physics_driver_netcdf

!--------------------------------------------------------------------
!    call the destructor routines for those modules who were initial-
!    ized from this module.
!--------------------------------------------------------------------
      call mpp_clock_begin ( turb_term_clock )
      call vert_turb_driver_end
      call mpp_clock_end ( turb_term_clock )
      call mpp_clock_begin ( diff_term_clock )
      call vert_diff_driver_end
      call mpp_clock_end ( diff_term_clock )


!--------------------------------------------------------------------
!    terminate radiation routines and data
!--------------------------------------------------------------------
      ! if (do_moist_processes) then   ! ZNT 06/23/2023: unconditionally call aerosol_end
      call mpp_clock_begin ( aerosol_term_clock )
      call aerosol_end (Aerosol_cld)
      call mpp_clock_end ( aerosol_term_clock )
      ! endif


      if (do_moist_processes) then  
        call mpp_clock_begin ( moist_processes_term_clock )
        call moist_processes_end ()
        call mpp_clock_end ( moist_processes_term_clock )
      endif

      call mpp_clock_begin ( tracer_term_clock )
      call atmos_tracer_driver_end
      call mpp_clock_end ( tracer_term_clock )
      call mpp_clock_begin ( damping_term_clock )
      call damping_driver_end
      call mpp_clock_end ( damping_term_clock )
      if (do_cosp) then
        call mpp_clock_begin ( cosp_term_clock )
        call cosp_driver_end
        call mpp_clock_end ( cosp_term_clock )
      endif

!---------------------------------------------------------------------
!    deallocate the module variables.
!---------------------------------------------------------------------
      deallocate (diff_cu_mo, diff_t, diff_m, pbltop, cush, cbmf,  &
                  hmint, cgust, pblhto, rkmo, taudpo, exist_shconv, &  ! h1g, 2017-01-31
                  exist_dpconv, & 
                  pblht_prev, hlsrc_prev, qtsrc_prev, cape_prev, cin_prev, & !h1g, 2017-01-31
                  convect, radturbten, qturbten, r_convect, &   ! ZNT: qturbten, 06/08/2021
                  tten_gw, cqa_gw, dqcdt_gw ) ! P1L 09/12/2024

      if (do_cosp) then
        deallocate ( temp_last, q_last)
        deallocate (&
           Precip_flux%fl_lsrain, Precip_flux%fl_lssnow,   &
           Precip_flux%fl_lsgrpl, Precip_flux%fl_ccrain,   &
           Precip_flux%fl_ccsnow)
      endif
 
      deallocate (id_tracer_phys_vdif_dn)
      deallocate (id_tracer_phys_vdif_up)
      deallocate (id_tracer_phys_turb)
      deallocate (id_tracer_phys_moist)

      call dealloc_physics_tendency_type (Physics_tendency)

      deallocate ( Restart%Cloud_data           )   ! h1g, 2017-02-02

!---------------------------------------------------------------------
!    mark the module as uninitialized.
!---------------------------------------------------------------------
      module_is_initialized = .false.


!-----------------------------------------------------------------------

 end subroutine physics_driver_end

!#######################################################################
subroutine physics_driver_restart(timestamp)

! write out restart file.
! Arguments: 
!   timestamp (optional, intent(in)) : A character string that represents the model time, 
!                                      used for writing restart. timestamp will append to
!                                      the any restart file name as a prefix. 

  character(len=*), intent(in), optional :: timestamp


  if (mpp_pe() == mpp_root_pe() ) then
     call error_mesg('physics_driver_mod', 'Writing netCDF formatted restart file: RESTART/physics_driver.res.nc', NOTE)
  endif
  call physics_driver_netcdf(timestamp)

  call moist_processes_restart(timestamp)
  call damping_driver_restart(timestamp)

end subroutine physics_driver_restart

subroutine physics_driver_netcdf(timestamp)

! Write out restart file for physics driver.
! This routine is needed so that physics_driver_restart and physics_driver_end
! can call a routine which will not result in multiple copies of restart files 
! being written by the destructor routines.
! Arguments: 
!   timestamp (optional, intent(in)) : A character string that represents the model time, 
!                                      used for writing restart. timestamp will append to
!                                      the any restart file name as a prefix. 


  character(len=*), intent(in), optional :: timestamp
  type(FmsNetcdfFile_t)       ::  Phy_restart !< Fms2io fileobj
  type(FmsNetcdfDomainFile_t) ::  Til_restart !< Fms2io domain decomposed fileobj
  logical :: tile_file_exist !< Flag indicating if the file was opened
  character(len=128) :: filename !< String of filename
  integer, allocatable, dimension(:)   :: pes !< Array of pes in the current pelist

  r_convect = 0.
  where(convect)
     r_convect = 1.0
  end where

  if (present(timestamp)) then
     filename = "RESTART/"//trim(timestamp)//".physics_driver.res.nc"
  else
     filename = "RESTART/physics_driver.res.nc"
  endif

  !< Get the current pelist
  allocate(pes(mpp_npes()))
  call mpp_get_current_pelist(pes)

  !< Open the scalar file with the current pelist, so that only the root pe opens and writes the file
  if (open_file(Phy_restart, trim(filename),"overwrite", is_restart=.true., pelist=pes)) then !scalar file
     call physics_driver_register_restart_scalars(Restart, Phy_restart)
     call write_restart(Phy_restart)
     call close_file(Phy_restart)
  endif
  deallocate(pes)

  if (mpp_get_ntile_count(physics_domain) == 1) then
     tile_file_exist = open_file(Til_restart, trim(filename),  "append", physics_domain, is_restart=.true.)
  else
     tile_file_exist = open_file(Til_restart, trim(filename),  "overwrite", physics_domain, is_restart=.true.)
  endif

  if (tile_file_exist) then !domain file
      call physics_driver_register_restart_domain(Restart, Til_restart)
      call write_restart(Til_restart)
      call add_domain_dimension_data(Til_restart)
      call close_file(Til_restart)
  endif

end subroutine physics_driver_netcdf

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

!#######################################################################
function do_moist_in_phys_up()

!--------------------------------------------------------------------
!    do_moist_in_phys_up returns the value of do_moist_processes
!----------------------------------------------------------------------

logical :: do_moist_in_phys_up

!---------------------------------------------------------------------
!    verify that the module is initialized.
!---------------------------------------------------------------------
      if ( .not. module_is_initialized) then
        call error_mesg ('do_moist_in_phys_up',  &
              'module has not been initialized', FATAL)
      endif
 
!-------------------------------------------------------------------
!    define output variable.
!-------------------------------------------------------------------
      do_moist_in_phys_up = do_moist_processes

 
end function do_moist_in_phys_up

!#####################################################################
! <FUNCTION NAME="get_diff_t">
!  <OVERVIEW>
!    returns the values of array diff_t
!  </OVERVIEW>
!  <DESCRIPTION>
!    returns the values of array diff_t
!  </DESCRIPTION>
!  <TEMPLATE>
!   diff_t(:,:,:) = get_diff_t()
!  </TEMPLATE>
! </FUNCTION>
!
!#####################################################################
function get_diff_t() result(diff_t_out)
real, dimension(size(diff_t,1),size(diff_t,2),size(diff_t,3)) :: diff_t_out

  if ( .not. module_is_initialized) then
    call error_mesg ('get_diff_t','module has not been initialized', FATAL)
  endif

  diff_t_out = diff_t

end function get_diff_t

!#####################################################################
function get_radturbten() result(radturbten_out)

!    returns the values of array radturbten

real, dimension(size(radturbten,1),size(radturbten,2),size(radturbten,3)) :: radturbten_out

  if ( .not. module_is_initialized) then
    call error_mesg ('get_radturbten','module has not been initialized', FATAL)
  endif

  radturbten_out = radturbten

end function get_radturbten
!#####################################################################
subroutine zero_radturbten()

!    sets all values of array radturbten to zero

  if ( .not. module_is_initialized) then
    call error_mesg ('zero_radturbten','module has not been initialized', FATAL)
  endif

  radturbten = 0.0

end subroutine zero_radturbten



!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!                    PRIVATE SUBROUTINES
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
          
               
     
!#####################################################################
subroutine physics_driver_register_restart_scalars (Restart, Phy_restart)

!    physics_driver_register_restart_scalars will register restart field when do_netcdf file
!    is true.

  type(clouds_from_moist_block_type), intent(inout), target :: Restart
  type(FmsNetcdfFile_t), intent(inout)       ::  Phy_restart !< Fms2io fileobj

  character(len=8), dimension(1)       :: dim_names !< Array of dimension names

  if (do_moist_processes) then
    if(doing_prog_clouds) then
       now_doing_strat = 1
    else
       now_doing_strat = 0
    endif

    if(doing_entrain) then
       now_doing_entrain = 1
    else
       now_doing_entrain = 0
    endif
  endif

  dim_names(1) = "Time"
  call register_axis(Phy_restart, dim_names(1), unlimited)

  call register_restart_field(Phy_restart, 'doing_strat',   now_doing_strat, dim_names)
  call register_restart_field(Phy_restart, 'doing_entrain', now_doing_entrain, dim_names)

  if (.not. Phy_restart%is_readonly) then !If not reading the file,
    call register_variable_attribute(Phy_restart, "doing_strat", "long_name", "doing_strat", str_len=len_trim("doing_strat"))
    call register_variable_attribute(Phy_restart, "doing_entrain", "long_name", "doing_entrain", str_len=len_trim("doing_entrain"))

    call register_variable_attribute(Phy_restart, "doing_strat", "units", "none", str_len=4)
    call register_variable_attribute(Phy_restart, "doing_entrain", "units", "none", str_len=4)
  endif

end subroutine physics_driver_register_restart_scalars

!#####################################################################
subroutine physics_driver_register_restart_domain (Restart, Til_restart)

!    physics_driver_register_restart_domain will register restart field when do_netcdf file
!    is true.

  type(clouds_from_moist_block_type), intent(inout), target :: Restart
  type(FmsNetcdfDomainFile_t), intent(inout) ::  Til_restart !< Fms2io domain decomposed fileobj

  integer           :: nc
  logical           :: reproduce_ulm_restart = .true.
  integer           :: index_strat
  character(len=8), dimension(4)             :: dim_names_4d, dim_names_4d2 !< Array of dimension names
  character(len=8), dimension(3)             :: dim_names_3d !< Array of dimension names

  dim_names_4d =  (/"xaxis_1", "yaxis_1", "zaxis_1", "Time   "/)
  dim_names_4d2 =  (/"xaxis_1", "yaxis_1", "zaxis_2", "Time   "/)
  dim_names_3d =  (/"xaxis_1", "yaxis_1", "Time   "/)

  call register_axis(Til_restart, "xaxis_1", "x")
  call register_axis(Til_restart, "yaxis_1", "y")
  call register_axis(Til_restart, "zaxis_1",  size(diff_cu_mo, 3))
  call register_axis(Til_restart, "zaxis_2", size(exist_shconv, 3))
  if (.not. Til_restart%mode_is_append) call register_axis(Til_restart, "Time", unlimited)

  !< Register the domain decomposed dimensions as variables so that the combiner can work
  !! correctly
  call register_field(Til_restart, "xaxis_1", "double", (/"xaxis_1"/))
  call register_field(Til_restart, "yaxis_1", "double", (/"yaxis_1"/))

  call register_restart_field(Til_restart, 'diff_cu_mo', diff_cu_mo, dim_names_4d)
  call register_restart_field(Til_restart, 'pbltop',     pbltop, dim_names_3d)
  call register_restart_field(Til_restart, 'cush',       cush, dim_names_3d, is_optional = .true.)
  call register_restart_field(Til_restart, 'cbmf',       cbmf, dim_names_3d, is_optional = .true.)
  call register_restart_field(Til_restart, 'hmint',      hmint, dim_names_3d, is_optional = .true.)
  call register_restart_field(Til_restart, 'cgust',      cgust, dim_names_3d, is_optional = .true.)
  call register_restart_field(Til_restart, 'pblhto',     pblhto, dim_names_3d, is_optional = .true.)
  call register_restart_field(Til_restart, 'rkmo',       rkmo, dim_names_3d, is_optional = .true.)
  call register_restart_field(Til_restart, 'taudpo',     taudpo, dim_names_3d, is_optional = .true.)
  call register_restart_field(Til_restart, 'exist_shconv', exist_shconv, dim_names_4d2, is_optional = .true.)
  call register_restart_field(Til_restart, 'exist_dpconv', exist_dpconv, dim_names_4d2, is_optional = .true.)
  call register_restart_field(Til_restart, 'pblht_prev',   pblht_prev, dim_names_4d2, is_optional = .true.)
  call register_restart_field(Til_restart, 'hlsrc_prev',   hlsrc_prev, dim_names_4d2, is_optional = .true.)
  call register_restart_field(Til_restart, 'qtsrc_prev',   qtsrc_prev, dim_names_4d2, is_optional = .true.)
  call register_restart_field(Til_restart, 'cape_prev',    cape_prev, dim_names_4d2, is_optional = .true.)
  call register_restart_field(Til_restart, 'cin_prev',     cin_prev, dim_names_4d2, is_optional = .true.)
  call register_restart_field(Til_restart, 'diff_t',     diff_t, dim_names_4d)
  call register_restart_field(Til_restart, 'diff_m',     diff_m, dim_names_4d)
  call register_restart_field(Til_restart, 'convect',    r_convect, dim_names_3d)
  call register_restart_field(Til_restart, 'tten_gw',  tten_gw, dim_names_4d, is_optional = .true. ) ! P1L
  call register_restart_field(Til_restart, 'cqa_gw',  cqa_gw, dim_names_4d, is_optional = .true. ) ! P1L
  call register_restart_field(Til_restart, 'dqcdt_gw',dqcdt_gw, dim_names_4d, is_optional = .true. ) ! P1L

  if (doing_prog_clouds) then
    call register_restart_field(Til_restart, 'radturbten',       radturbten, dim_names_4d)
    call register_restart_field(Til_restart,   'qturbten',         qturbten, dim_names_4d, is_optional = .true.)   ! ZNT
  endif

  index_strat = 0
  do nc = 1, size(Restart%Cloud_data,1)
    if (trim(Restart%Cloud_data(nc)%scheme_name).eq.'strat_cloud' .and. .not. reproduce_ulm_restart) then
      call register_restart_field(Til_restart, 'lsc_cloud_area',     Restart%Cloud_data(nc)%cloud_area, dim_names_4d,     is_optional = .true.)
      call register_restart_field(Til_restart, 'lsc_liquid',         Restart%Cloud_data(nc)%liquid_amt, dim_names_4d,     is_optional = .true.)
      call register_restart_field(Til_restart, 'lsc_ice',            Restart%Cloud_data(nc)%ice_amt, dim_names_4d,        is_optional = .true.)
      call register_restart_field(Til_restart, 'lsc_droplet_number', Restart%Cloud_data(nc)%droplet_number, dim_names_4d, is_optional = .true.)
      call register_restart_field(Til_restart, 'lsc_ice_number',     Restart%Cloud_data(nc)%ice_number, dim_names_4d,     is_optional = .true.)
      call register_restart_field(Til_restart, 'lsc_snow',           Restart%Cloud_data(nc)%snow, dim_names_4d,           is_optional = .true.)
      call register_restart_field(Til_restart, 'lsc_rain',           Restart%Cloud_data(nc)%rain, dim_names_4d,           is_optional = .true.)
      call register_restart_field(Til_restart, 'lsc_snow_size',      Restart%Cloud_data(nc)%snow_size, dim_names_4d,      is_optional = .true.)
      call register_restart_field(Til_restart, 'lsc_rain_size',      Restart%Cloud_data(nc)%rain_size, dim_names_4d,      is_optional = .true.)
    endif
    if (trim(Restart%Cloud_data(nc)%scheme_name).eq.'strat_cloud' .and. reproduce_ulm_restart) index_strat = nc

    if (trim(Restart%Cloud_data(nc)%scheme_name).eq.'uw_conv') then
      call register_restart_field(Til_restart, 'shallow_cloud_area',     Restart%Cloud_data(nc)%cloud_area, dim_names_4d,     is_optional = .true.)
      call register_restart_field(Til_restart, 'shallow_liquid',         Restart%Cloud_data(nc)%liquid_amt, dim_names_4d,     is_optional = .true.)
      call register_restart_field(Til_restart, 'shallow_ice',            Restart%Cloud_data(nc)%ice_amt, dim_names_4d,        is_optional = .true.)
      call register_restart_field(Til_restart, 'shallow_droplet_number', Restart%Cloud_data(nc)%droplet_number, dim_names_4d, is_optional = .true.)
      call register_restart_field(Til_restart, 'shallow_ice_number',     Restart%Cloud_data(nc)%ice_number, dim_names_4d,     is_optional = .true.)
    endif
  enddo

    ! save large-scale clouds last to reproduce ulm code
    if (index_strat > 0) then
      nc = index_strat
      call register_restart_field(Til_restart, 'lsc_cloud_area',     Restart%Cloud_data(nc)%cloud_area, dim_names_4d,     is_optional = .true.)
      call register_restart_field(Til_restart, 'lsc_liquid',         Restart%Cloud_data(nc)%liquid_amt, dim_names_4d,     is_optional = .true.)
      call register_restart_field(Til_restart, 'lsc_ice',            Restart%Cloud_data(nc)%ice_amt, dim_names_4d,        is_optional = .true.)
      call register_restart_field(Til_restart, 'lsc_droplet_number', Restart%Cloud_data(nc)%droplet_number, dim_names_4d, is_optional = .true.)
      call register_restart_field(Til_restart, 'lsc_ice_number',     Restart%Cloud_data(nc)%ice_number, dim_names_4d,     is_optional = .true.)
      call register_restart_field(Til_restart, 'lsc_snow',           Restart%Cloud_data(nc)%snow, dim_names_4d,           is_optional = .true.)
      call register_restart_field(Til_restart, 'lsc_rain',           Restart%Cloud_data(nc)%rain, dim_names_4d,           is_optional = .true.)
      call register_restart_field(Til_restart, 'lsc_snow_size',      Restart%Cloud_data(nc)%snow_size, dim_names_4d,      is_optional = .true.)
      call register_restart_field(Til_restart, 'lsc_rain_size',      Restart%Cloud_data(nc)%rain_size, dim_names_4d,      is_optional = .true.)
    endif

end subroutine physics_driver_register_restart_domain

! </SUBROUTINE>    
!#####################################################################
subroutine check_args (lat, lon, area, p_half, p_full, z_half, z_full,&
                        u, v, t, q, r, um, vm, tm, qm, rm,             &
                        udt, vdt, tdt, qdt, rdt, rdiag)

!----------------------------------------------------------------------
!    check_args determines if the input arrays to physics_driver_down
!    are of a consistent size.
!-----------------------------------------------------------------------

real,    dimension(:,:),    intent(in)           :: lat, lon, area
real,    dimension(:,:,:),  intent(in)           :: p_half, p_full,   &
                                                    z_half, z_full,   &
                                                    u, v, t, q, um, vm, &
                                                    tm, qm
real,    dimension(:,:,:,:),intent(in)           :: r, rm
real,    dimension(:,:,:),  intent(in)           :: udt, vdt, tdt, qdt
real,    dimension(:,:,:,:),intent(in)           :: rdt
real,    dimension(:,:,:,ntp+1:),intent(in)      :: rdiag

!-----------------------------------------------------------------------
!   intent(in) variables:
!
!      lat            latitude of model points [ radians ]
!      lon            longitude of model points [ radians ]
!      area           grid box area - currently not used [ m**2 ]
!      p_half         pressure at half levels (offset from t,q,u,v,r)
!                     [ Pa ]
!      p_full         pressure at full levels [ Pa }
!      z_half         height at half levels [ m ]
!      z_full         height at full levels [ m ]
!      u              zonal wind at current time step [ m / s ]
!      v              meridional wind at current time step [ m / s ]
!      t              temperature at current time step [ deg k ]
!      q              specific humidity at current time step  kg / kg ]
!      r              multiple 3d tracer fields at current time step
!      um,vm          zonal and meridional wind at previous time step
!      tm,qm          temperature and specific humidity at previous 
!                     time step
!      rm             multiple 3d tracer fields at previous time step
!      udt            zonal wind tendency [ m / s**2 ]
!      vdt            meridional wind tendency [ m / s**2 ]
!      tdt            temperature tendency [ deg k / sec ]
!      qdt            specific humidity tendency 
!                     [  kg vapor / kg air / sec ]
!      rdt            multiple tracer tendencies [ unit / unit / sec ]
!
!   intent(in), optional:
!
!---------------------------------------------------------------------

!----------------------------------------------------------------------
!   local variables:

      integer ::  id, jd, kd  ! model dimensions on the processor  
      integer ::  ierr        ! error flag

!--------------------------------------------------------------------
!    define the sizes that the arrays should be.
!--------------------------------------------------------------------
      id = size(u,1) 
      jd = size(u,2) 
      kd = size(u,3) 

!--------------------------------------------------------------------
!    check the dimensions of each input array. if they are incompat-
!    ible in size with the standard, the error flag is set to so
!    indicate.
!--------------------------------------------------------------------
      ierr = 0
      ierr = ierr + check_dim (lat, 'lat',  id,jd)
      ierr = ierr + check_dim (lon, 'lon',  id,jd)
      ierr = ierr + check_dim (area,'area', id,jd)

      ierr = ierr + check_dim (p_half,'p_half', id,jd,kd+1)
      ierr = ierr + check_dim (p_full,'p_full', id,jd,kd)
      ierr = ierr + check_dim (z_half,'z_half', id,jd,kd+1)
      ierr = ierr + check_dim (z_full,'z_full', id,jd,kd)

      ierr = ierr + check_dim (u, 'u',  id,jd,kd)
      ierr = ierr + check_dim (v, 'v',  id,jd,kd)
      ierr = ierr + check_dim (t, 't',  id,jd,kd)
      ierr = ierr + check_dim (q, 'q',  id,jd,kd)
      ierr = ierr + check_dim (um,'um', id,jd,kd)
      ierr = ierr + check_dim (vm,'vm', id,jd,kd)
      ierr = ierr + check_dim (tm,'tm', id,jd,kd)
      ierr = ierr + check_dim (qm,'qm', id,jd,kd)

      ierr = ierr + check_dim (udt,'udt', id,jd,kd)
      ierr = ierr + check_dim (vdt,'vdt', id,jd,kd)
      ierr = ierr + check_dim (tdt,'tdt', id,jd,kd)
      ierr = ierr + check_dim (qdt,'qdt', id,jd,kd)

      if (ntp > 0) then
        ierr = ierr + check_dim (r,  'r',   id,jd,kd,ntp)
        ierr = ierr + check_dim (rm, 'rm',  id,jd,kd,ntp)
      endif
      if (ntp > 0) then
        ierr = ierr + check_dim (rdt,'rdt', id,jd,kd,ntp)
      endif
      if (nt > ntp) then
        ierr = ierr + check_dim (rdiag,'rdiag', id,jd,kd,nt-ntp)
      endif

!--------------------------------------------------------------------
!    if any problems were detected, exit with an error message.
!--------------------------------------------------------------------
      if (ierr > 0) then
        call error_mesg ('physics_driver_mod', 'bad dimensions', FATAL)
      endif

!-----------------------------------------------------------------------


      end subroutine check_args


!#######################################################################
function check_dim_2d (data,name,id,jd) result (ierr)

!--------------------------------------------------------------------
!    check_dim_2d compares the size of two-dimensional input arrays
!    with supplied expected dimensions and returns an error if any
!    inconsistency is found.
!--------------------------------------------------------------------

real,    intent(in), dimension(:,:) :: data
character(len=*), intent(in)        :: name
integer, intent(in)                 :: id, jd
integer                             :: ierr

!---------------------------------------------------------------------
!  intent(in) variables:
!
!     data        array to be checked
!     name        name associated with array to be checked
!     id, jd      expected i and j dimensions
!     
!  result variable:
!
!     ierr        set to 0 if ok, otherwise is a count of the number
!                 of incompatible dimensions
!
!--------------------------------------------------------------------

      ierr = 0
      if (size(data,1) /= id) then
        call error_mesg ('physics_driver_mod',  &
             'dimension 1 of argument ' //  &
              name(1:len_trim(name)) // ' has wrong size.', NOTE)
        ierr = ierr + 1
      endif
      if (size(data,2) /= jd) then
           call error_mesg ('physics_driver_mod',  &
                'dimension 2 of argument ' //  &
                name(1:len_trim(name)) // ' has wrong size.', NOTE)
           ierr = ierr + 1
      endif

!----------------------------------------------------------------------

      end function check_dim_2d

!#######################################################################
function check_dim_3d (data,name,id,jd,kd) result (ierr)

!--------------------------------------------------------------------
!    check_dim_3d compares the size of thr1eedimensional input arrays
!    with supplied expected dimensions and returns an error if any
!    inconsistency is found.
!--------------------------------------------------------------------

real,    intent(in), dimension(:,:,:) :: data
character(len=*), intent(in)          :: name
integer, intent(in)                   :: id, jd, kd
integer  ierr

!---------------------------------------------------------------------
!  intent(in) variables:
!
!     data        array to be checked
!     name        name associated with array to be checked
!     id, jd,kd   expected i, j and k dimensions
!     
!  result variable:
!
!     ierr        set to 0 if ok, otherwise is a count of the number
!                 of incompatible dimensions
!
!--------------------------------------------------------------------

      ierr = 0
      if (size(data,1) /= id) then
        call error_mesg ('physics_driver_mod',  &
                'dimension 1 of argument ' //  &
                name(1:len_trim(name)) // ' has wrong size.', NOTE)
        ierr = ierr + 1
      endif
      if (size(data,2) /= jd) then
        call error_mesg ('physics_driver_mod',  &
              'dimension 2 of argument ' //  &
                name(1:len_trim(name)) // ' has wrong size.', NOTE)
        ierr = ierr + 1
      endif
      if (size(data,3) /= kd) then
        call error_mesg ('physics_driver_mod',  &
                'dimension 3 of argument ' //  &
                name(1:len_trim(name)) // ' has wrong size.', NOTE)
        ierr = ierr + 1
      endif

!---------------------------------------------------------------------


      end function check_dim_3d


!#######################################################################
function check_dim_4d (data,name,id,jd,kd,nt) result (ierr)

!--------------------------------------------------------------------
!    check_dim_4d compares the size of four dimensional input arrays
!    with supplied expected dimensions and returns an error if any
!    inconsistency is found.
!--------------------------------------------------------------------
real,    intent(in), dimension(:,:,:,:) :: data
character(len=*), intent(in)            :: name
integer, intent(in)                     :: id, jd, kd, nt
integer                                 :: ierr

!---------------------------------------------------------------------
!  intent(in) variables:
!
!     data          array to be checked
!     name          name associated with array to be checked
!     id,jd,kd,nt   expected i, j and k dimensions
!     
!  result variable:
!
!     ierr          set to 0 if ok, otherwise is a count of the number
!                   of incompatible dimensions
!
!--------------------------------------------------------------------

      ierr = 0
      if (size(data,1) /= id) then
        call error_mesg ('physics_driver_mod',  &
                'dimension 1 of argument ' //  &
                name(1:len_trim(name)) // ' has wrong size.', NOTE)
        ierr = ierr + 1
      endif
      if (size(data,2) /= jd) then
        call error_mesg ('physics_driver_mod',  &
                'dimension 2 of argument ' //  &
                name(1:len_trim(name)) // ' has wrong size.', NOTE)
        ierr = ierr + 1
      endif
      if (size(data,3) /= kd) then
        call error_mesg ('physics_driver_mod',  &
                'dimension 3 of argument ' //  &
                name(1:len_trim(name)) // ' has wrong size.', NOTE)
        ierr = ierr + 1
      endif
      if (size(data,4) /= nt) then
        call error_mesg ('physics_driver_mod',  &
                'dimension 4 of argument ' //  &
                name(1:len_trim(name)) // ' has wrong size.', NOTE)
        ierr = ierr + 1
      endif

!---------------------------------------------------------------------


      end function check_dim_4d


 
end module physics_driver_mod
