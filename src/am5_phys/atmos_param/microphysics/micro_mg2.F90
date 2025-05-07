module micro_mg2_mod

! this is ncar routine micro_mg2

!---------------------------------------------------------------------------------
! Purpose:
!   MG microphysics version 2, prognostic precipitation.
!       point for the development of MG2
!
! Author: Andrew Gettelman, Hugh Morrison.
! Contributions from: Peter Caldwell, Xiaohong Liu and Steve Ghan
!!
! NOTE: If do_cldice is false, then MG microphysics should not update CLDICE
! or NUMICE; however, it is assumed that the other microphysics scheme will have
! updated CLDICE and NUMICE. The other microphysics should handle the following
! processes that would have been done by MG:
!   - Detrainment (liquid and ice)
!   - Homogeneous ice nucleation
!   - Heterogeneous ice nucleation
!   - Bergeron process
!   - Melting of ice
!   - Freezing of cloud drops
!   - Autoconversion (ice -> snow)
!   - Growth/Sublimation of ice
!   - Sedimentation of ice
!---------------------------------------------------------------------------------
! Based on micro_mg (restructuring of former cldwat2m_micro)
! Author: Andrew Gettelman, Hugh Morrison.
! Contributions from: Xiaohong Liu and Steve Ghan
! December 2005-May 2010
! Description in: Morrison and Gettelman, 2008. J. Climate (MG2008)
!                 Gettelman et al., 2010 J. Geophys. Res. - Atmospheres (G2010)
! for questions contact Hugh Morrison, Andrew Gettelman
! e-mail: morrison@ucar.edu, andrew@ucar.edu
!---------------------------------------------------------------------------------

!--------------------------------------------------------------------------
! Interfaces, diagnostics, used modules, constants modified for use in GFDL
! based models
! Huan Guo
!--------------------------------------------------------------------------

! Code comments added by HM, 093011
! General code structure:
!
! Code is divided into two main subroutines:
!   subroutine micro_mg_init --> initializes microphysics routine, should be called
!                                  once at start of simulation
!   subroutine micro_mg_tend --> main microphysics routine to be called each time step
!
! List of external functions:
!   qsat_water --> for calculating saturation vapor pressure with respect to liquid water
!   qsat_ice --> for calculating saturation vapor pressure with respect to ice
!   gamma   --> standard mathematical gamma function
! .........................................................................
! List of inputs through use statement in fortran90:
! Variable Name                      Description                Units
! .........................................................................
! gravit          acceleration due to gravity                    m s-2
! rair            dry air gas constant for air                  J kg-1 K-1
! tmelt           temperature of melting point for water          K
! cpair           specific heat at constant pressure for dry air J kg-1 K-1
! rh2o            gas constant for water vapor                  J kg-1 K-1
! latvap          latent heat of vaporization                   J kg-1
! latice          latent heat of fusion                         J kg-1
! qsat_water      external function for calculating liquid water
!                 saturation vapor pressure/humidity              -
! qsat_ice        external function for calculating ice
!                 saturation vapor pressure/humidity              pa
! rhmini          relative humidity threshold parameter for
!                 nucleating ice                                  -
! .........................................................................
! NOTE: List of all inputs/outputs passed through the call/subroutine statement
!       for micro_mg_tend is given below at the start of subroutine micro_mg_tend.
!---------------------------------------------------------------------------------

! Procedures required:
! 1) An implementation of the gamma function (if not intrinsic).
! 2) saturation vapor pressure and specific humidity over water
! 3) svp over ice

use gamma_mg_mod,              only: gamma =>gamma_mg
use lscloud_types_mod,         only: diag_id_type, diag_pt_type

use mpp_mod,                   only: input_nml_file
use fms_mod,                   only: mpp_pe, error_mesg,  &
                                     FATAL, &
                                     stdlog, write_version_number, &
                                     check_nml_error, &
                                     mpp_root_pe,  mpp_chksum
use sat_vapor_pres_mod,         only: lookup_es2, lookup_es3, compute_qs
use physics_radiation_exch_mod, only : exchange_control_type

! Parameters from the utilities module.
use micro_mg2_utils, only: &
     r8, &
     pi, &
     omsm, &
     qsmall, &
     mincld, &
     rhosn, &
     rhoi, &
     rhow, &
     rhows, &
     ac, bc, &
     ai, bi, &
     aj, bj, &
     ar, br, &
     as, bs, &
     mi0, &
     rising_factorial

! Constituent properties.
use micro_mg2_utils, only: &
       mg_liq_props, &
       mg_ice_props, &
       mg_rain_props, &
       mg_snow_props

! Size calculation functions.
use micro_mg2_utils, only: &
       size_dist_param_liq, &
       size_dist_param_basic, &
       avg_diameter


use micro_mg2_utils, only: &
     micro_mg_utils_init, &
     size_dist_param_liq, &
     size_dist_param_basic, &
     avg_diameter, &
     rising_factorial, &
     ice_deposition_sublimation, &
     sb2001v2_liq_autoconversion,&
     sb2001v2_accre_cld_water_rain,&
     kk2000_liq_autoconversion, &
     ice_autoconversion, &
     immersion_freezing, &
     contact_freezing, &
     snow_self_aggregation, &
     accrete_cloud_water_snow, &
     secondary_ice_production, &
     accrete_rain_snow, &
     heterogeneous_rain_freezing, &
     accrete_cloud_water_rain, &
     self_collection_rain, &
     accrete_cloud_ice_snow, &
     evaporate_sublimate_precip, &
     bergeron_process_snow,      & ! h1g, 2020-03-23
     cotton_liq_autoconversion     ! h1g, 2020-03-23


implicit none
private
save

public :: &
     micro_mg2_init, &
     micro_mg2_get_cols, &
     micro_mg2_tend

!------------------------------------------------------------------------
!--version number--------------------------------------------------------

character(len=128) :: Version = '$Id: micro_mg2.F90,v 1.1.2.2.2.1 2017/01/23 14:50:42 Huan.Guo Exp $'
character(len=128) :: Tagname = '$Name:  $'

!=========================================================
! Constants set in initialization
!=========================================================

! Set using arguments to micro_mg_init
real(r8) :: g           ! gravity
real(r8) :: r           ! dry air gas constant
real(r8) :: rv          ! water vapor gas constant
real(r8) :: cpp         ! specific heat of dry air
real(r8) :: tmelt       ! freezing point of water (K)

! latent heats of:
real(r8) :: xxlv        ! vaporization
real(r8) :: xlf         ! freezing
real(r8) :: xxls        ! sublimation

real(r8) :: rhosu       ! typical 850mn air density

real(r8) :: icenuct     ! ice nucleation temperature: currently -5 degrees C

real(r8) :: snowmelt    ! what temp to melt all snow: currently 2 degrees C
real(r8) :: rainfrze    ! what temp to freeze all rain: currently -5 degrees C

real       :: dcs       !autoconversion size threshold
logical    :: rho_factor_in_max_vt = .true.
real       :: max_rho_factor_in_vt = 1.0

! switch for specification rather than prediction of droplet and crystal number
! note: number will be adjusted as needed to keep mean size within bounds,
! even when specified droplet or ice number is used

! ***note: Even if constant cloud ice number is set, ice number is allowed
! to evolve based on process rates. This is needed in order to calculate
! the change in mass due to ice nucleation. All other ice microphysical
! processes are consistent with the specified constant ice number if
! this switch is turned on.

! nccons = .true. to specify constant cloud droplet number
! nicons = .true. to specify constant cloud ice number

logical :: nccons = .false.
logical :: nicons = .false.

! parameters for specified ice and droplet number concentration
! note: these are local in-cloud values, not grid-mean
real(r8) :: ncnst = 100.e6_r8    ! droplet num concentration when
                                 ! nccons=.true. (m-3)
real(r8) :: ninst = 0.1e6_r8     ! ice num concentration when
                                 ! nicons=.true. (m-3)
! <---h1g, 2012-06-12
logical           :: liu_in = .false.
                             ! True = Liu et al 2007 Ice nucleation
                             ! False = cooper fixed ice nucleation (MG2008)

logical           :: use_Meyers = .false.
! Ni (/m3) = 1000 * exp( (12.96* [(esl-esi)/esi]) - 0.639 )
!  Figure 9.3 of Rogers and Yau (1998) shows the nearly linear
!       variation of [(esl-esi)/esi] from 0. at 273.16K to 0.5 at
!       233.16K.  Analytically this is parameterized as (tfreeze-T)/80.
!
!  Ni (/m3) = 1000 * exp( 12.96* (tfreeze-T)/80 - 0.639 )

logical           :: use_Fan2019 = .false.
real(r8)          :: Nice_max_Fan = 40.0

!---> h1g, 2014-05-19
real              :: tc_cooper   = -35.0
real              :: SmallIceFallFac  = 1.0  ! h1g, 2023-06-12
real              :: IceFallFac  = 1.0
real              :: SnowFallFac = 1.0
real              :: RainFallFac = 1.0       ! h1g, 2023-06-12

logical           :: include_contact_freeze_in_berg = .false.

character(len=16)  :: micro_mg_precip_frac_method = "max_overlap"  ! type of precipitation fraction method
real(r8)           :: micro_mg_bergs_eff_factor = 1.0_r8     ! bergs efficiency factor (liquid to snow)
real(r8)           :: micro_mg_berg_eff_factor  = 1.0_r8      ! berg efficiency factor (liquid to ice)

logical            :: allow_sed_supersat = .true. ! Allow supersaturated conditions after sedimentation loop

logical            :: do_sb_physics = .false.
logical            :: use_hetfrz_classnuc = .false.
real               :: tau_act_liq = 1800.0
real               :: tau_act_ice = 900.0
real               :: tc_act  = -35.0

logical            :: allow_rain_num_evap        = .false.
logical            :: allow_snow_num_sublimation = .false.

logical            :: no_evap_in_sedimentation = .false.
real(r8)           :: vfactor = 1.0
real(r8)           :: vfac_drop = 1.0   ! h1g, 2020-06-18
real(r8)           :: vfac_ice  = 1.0   ! h1g, 2020-06-18

real(r8)           :: icld_cri = -0.2    ! h1g, 2020-07-02
real(r8)           :: evap_subl_fac = 1.0  ! h1g, 2020-07-06

real(r8)           :: evap_RH_cri = 1.0  ! h1g, 2023-06-26
real(r8)           :: subl_RH_cri = 1.0  ! h1g, 2023-06-26


real(r8)           :: ice_sublim_factor = 1.0 ! h1g, 2020-07-16

! minimum mass of new crystal due to freezing of cloud droplets done
! externally (kg)
real(r8), parameter :: mi0l_min = 4._r8/3._r8*pi*rhow*(4.e-6_r8)**3

real              ::  rhmini=0.80     ! minimum rh for ice cld fraction > 0
logical           ::  microp_uniform = .false.
                                      ! .true. = configure uniform for
                                      ! sub-columns
                                      ! .false. = use w/o sub-columns
                                      ! (default)
logical           ::  do_cldice = .true.
                                      ! .true. = do all processes (default)
                                      ! .false. = skip all processes
                                      ! affecting cloud ice
logical           ::  do_ice_nucl_wpdf
logical           ::  do_Ni_linear_interp = .false.
logical           ::  do_implicit_fall    = .false.


logical           ::  do_qc_implicit_fall    = .true.
logical           ::  do_qi_implicit_fall    = .true.
logical           ::  do_qr_implicit_fall    = .true.
logical           ::  do_qs_implicit_fall    = .true.

logical           ::  include_ice_in_snowflx = .false.   ! add ice flux into snow flux for aerosol wet scavenging, h1g, 2024-01-31

! additional constants to help speed up code
real(r8) :: gamma_br_plus1
real(r8) :: gamma_br_plus4
real(r8) :: gamma_bs_plus1
real(r8) :: gamma_bs_plus4
real(r8) :: gamma_bi_plus1
real(r8) :: gamma_bi_plus4
real(r8) :: gamma_bj_plus1
real(r8) :: gamma_bj_plus4
real(r8) :: xxlv_squared
real(r8) :: xxls_squared
!<--- h1g, 2014-05-19

real(r8)           ::  dum_5, dum_30, dum_tmp
integer            ::  iter = 1
logical            ::  do_liq_num_adjust = .true.
logical            ::  do_liq_num_riming = .true.
logical            ::  do_liq_num_ihom   = .true.  ! h1g, 2020-06-22
logical            ::  do_ice_num_adjust = .true.  ! h1g, 2020-07-01

logical            ::  do_cotton_auto    = .false.
real(r8)           ::  rthresh = 8.6
logical            ::  do_HM_splinter    = .true.
logical            ::  remove_super_RK   = .false.
logical            ::  use_const_ELI     = .false.    ! --> h1g, 2020-04-16
real(r8)           ::  ELI_RK = 0.7                   ! --> h1g, 2020-04-16
logical            ::  use_FanAndCooper  = .false.    ! --> h1g, 2020-04-18
real(r8)           ::  sublim_factor = 0.0
real(r8)           ::  ice_nucl_factor = 1.0

logical            ::  include_homogeneous_for_wetdep = .true.  ! --> h1g, 2024-01-31, will change to false in NML
namelist / micro_mg2_nml /   &
                 max_rho_factor_in_vt, &
                 rho_factor_in_max_vt, &
                 nccons, ncnst,                &  ! cjg
                 nicons, ninst,                &  ! cjg
                 liu_in,                       &  ! h1g
                 use_Meyers, tc_cooper,      & !h1g
                 SmallIceFallFac,            & !h1g, 2023-06-12
                 IceFallFac, SnowFallFac,    & !h1g
                 RainFallFac,            & !h1g, 2023-06-12
                 include_contact_freeze_in_berg, & !h1g
                 do_sb_physics,  allow_sed_supersat, & ! h1g
                 use_hetfrz_classnuc, micro_mg_precip_frac_method, micro_mg_berg_eff_factor,  &
                 tau_act_liq, tau_act_ice, tc_act, allow_rain_num_evap, allow_snow_num_sublimation, &
                 rhmini, microp_uniform,  do_cldice, do_Ni_linear_interp, &
                 do_implicit_fall, vfactor, no_evap_in_sedimentation, iter, &
                 use_Fan2019, Nice_max_Fan, &
                 do_qc_implicit_fall,  do_qi_implicit_fall, do_qr_implicit_fall, do_qs_implicit_fall, &
                 do_liq_num_adjust, do_liq_num_riming, do_cotton_auto, rthresh, do_HM_splinter,       & ! h1g, 2020-03-06
                 remove_super_RK, use_const_ELI, ELI_RK, use_FanAndCooper, sublim_factor, &  ! h1g, 2020-04-18
                 ice_nucl_factor, vfac_drop, vfac_ice, do_liq_num_ihom, micro_mg_bergs_eff_factor, &  ! h1g, 2020-06-22
                 do_ice_num_adjust, icld_cri, evap_subl_fac, ice_sublim_factor,  &     ! h1g, 2020-07-06
                 evap_RH_cri,  subl_RH_cri,                                      &     ! h1g, 2023-06-26
                 include_homogeneous_for_wetdep, include_ice_in_snowflx                ! h1g, 2024-01-31
                 

!===============================================================================
contains
!===============================================================================

subroutine micro_mg2_init( &
     kind, gravit, rair, rh2o, cpair,    &
     tmelt_in, latvap, latice,           &
    ! rhmini_in, microp_uniform_in, do_cldice_in, &
     do_ice_nucl_wpdf_in, errstring, Exch_ctrl)

  !-----------------------------------------------------------------------
  !
  ! Purpose:
  ! initialize constants for MG microphysics
  !
  ! Author: Andrew Gettelman Dec 2005
  !
  !-----------------------------------------------------------------------

  integer,  intent(in)  :: kind         ! Kind used for reals
  real(r8), intent(in)  :: gravit
  real(r8), intent(in)  :: rair
  real(r8), intent(in)  :: rh2o
  real(r8), intent(in)  :: cpair
  real(r8), intent(in)  :: tmelt_in     ! Freezing point of water (K)
  real(r8), intent(in)  :: latvap
  real(r8), intent(in)  :: latice

  logical, intent(in) :: do_ice_nucl_wpdf_in
  type(exchange_control_type), intent(in) :: Exch_ctrl

  character(128), intent(out) :: errstring    ! Output status (non-blank for error return)

  INTEGER   :: io, ierr, logunit

  !-----------------------------------------------------------------------

  dcs = Exch_ctrl%dcs
  call micro_mg_utils_init(kind, rh2o, cpair, tmelt_in, latvap, latice, &
       dcs, errstring)

  errstring = ' '

  if( kind .ne. r8 ) then
     errstring = 'micro_mg2_init: KIND of reals does not match'
     return
  endif

  ! declarations for MG code (transforms variable names)

  g= gravit                 ! gravity
  r= rair                   ! dry air gas constant: note units(phys_constants are in J/K/kmol)
  rv= rh2o                  ! water vapor gas constant
  cpp = cpair               ! specific heat of dry air
  tmelt = tmelt_in

! latent heats
  xxlv = latvap         ! latent heat vaporization
  xlf  = latice         ! latent heat freezing
  xxls = xxlv + xlf     ! latent heat of sublimation

!---------------------------------------------------------------
!     process namelist
!---------------------------------------------------------------
  read (input_nml_file, nml=micro_mg2_nml, iostat=io)
  ierr = check_nml_error(io,'micro_mg2_nml')

!-----------------------------------------------------------------------
!    write version and namelist to stdlog.
!-----------------------------------------------------------------------
  call write_version_number (version, tagname)
  logunit = stdlog()
  if (mpp_pe() == mpp_root_pe()) &
        write (logunit, nml=micro_mg2_nml)

  do_ice_nucl_wpdf = do_ice_nucl_wpdf_in

  ! typical air density at 850 mb
  rhosu = 85000._r8/(rair * tmelt)

  ! Maximum temperature at which snow is allowed to exist
  snowmelt = tmelt + 2._r8
  ! Minimum temperature at which rain is allowed to exist
  rainfrze = tmelt - 40._r8

  ! Ice nucleation temperature
  icenuct  = tmelt - 5._r8

  ! Define constants to help speed up code (this limits calls to gamma function)
  gamma_br_plus1=gamma(1._r8+br)
  gamma_br_plus4=gamma(4._r8+br)
  gamma_bs_plus1=gamma(1._r8+bs)
  gamma_bs_plus4=gamma(4._r8+bs)
  gamma_bi_plus1=gamma(1._r8+bi)
  gamma_bi_plus4=gamma(4._r8+bi)
  gamma_bj_plus1=gamma(1._r8+bj)
  gamma_bj_plus4=gamma(4._r8+bj)

  xxlv_squared=xxlv**2
  xxls_squared=xxls**2

end subroutine micro_mg2_init

!===============================================================================
!microphysics routine for each timestep goes here...
subroutine micro_mg2_tend (  lon, lat, &
      dqa_activation, total_activation, tiedtke_macrophysics, &
      j, jdim, &
      mgncol,             nlev,               deltatin,           &
      concen_dust_sub,                                            &
      tn,                           qn,                           &
      qcn,                          qin,                          &
      ncn,                          nin,                          &
      qrn,                          qsn,                          &
      nrn,                          nsn,                          &
      relvar,              accre_enhan,                           &
      p,                   pdel,       zhalf,                     &
      cldn,               liqcldf,            icecldf,            &
      delta_cf, D_eros_l, nerosc, D_eros_i, nerosi, dqcdt, dqidt, &
      naai,               npccn,                                  &
      rndst,              nacon,  &
      tlat,               qvlat,  &
      qctend,             qitend, &
      nctend,             nitend, &
      qrtend,             qstend, &
      nrtend,             nstend, &
      prect,              preci,  &
      qsout,             rflx,    sflx,    &
      qrout,             reff_rain,         reff_snow,   &
      errstring,                                         &
      f_snow_berg,       ssat_disposal,                  &
      n_diag_4d, diag_4l, diag_id, diag_pt)


  ! input arguments
  real(r8), intent(in) :: lon(mgncol), lat(mgncol)

  logical, intent (in) :: dqa_activation
  logical, intent (in) :: total_activation
  logical, intent (in) :: tiedtke_macrophysics
  integer,  intent(in) :: j, jdim
  integer,  intent(in) :: mgncol                ! number of microphysics columns
  integer,  intent(in) :: nlev                  ! number of layers
  real(r8), intent(in) :: deltatin              ! time step (s)

  real(r8), intent(in) :: concen_dust_sub(mgncol,nlev)  ! sub-micro dust aerosol concentration (ug/m3)

  real(r8), intent(in) :: tn(mgncol,nlev)               ! input temperature (K)
  real(r8), intent(in) :: qn(mgncol,nlev)               ! input h20 vapor mixing ratio (kg/kg)
  real(r8), intent(in) :: relvar(mgncol,nlev)          ! relative variance of cloud water (-)
  real(r8), intent(in) :: accre_enhan(mgncol,nlev)     ! optional accretion enhancement factor (-)

  ! note: all input cloud variables are grid-averaged
  real(r8), intent(in) :: qcn(mgncol,nlev)       ! cloud water mixing ratio (kg/kg)
  real(r8), intent(in) :: qin(mgncol,nlev)       ! cloud ice mixing ratio (kg/kg)
  real(r8), intent(in) :: ncn(mgncol,nlev)       ! cloud water number conc (1/kg)
  real(r8), intent(in) :: nin(mgncol,nlev)       ! cloud ice number conc (1/kg)

  real(r8), intent(in) :: qrn(mgncol,nlev)       ! rain mixing ratio (kg/kg)
  real(r8), intent(in) :: qsn(mgncol,nlev)       ! snow mixing ratio (kg/kg)
  real(r8), intent(in) :: nrn(mgncol,nlev)       ! rain number conc (1/kg)
  real(r8), intent(in) :: nsn(mgncol,nlev)       ! snow number conc (1/kg)

  real(r8), intent(in) :: p(mgncol,nlev)         ! air pressure (pa)
  real(r8), intent(in) :: pdel(mgncol,nlev)      ! pressure difference across level (pa)

  real(r8), intent(in) :: zhalf(mgncol,nlev+1)   ! half-pressure level height (m)
  ! hm add 11-16-11, interface pressure

  real(r8), intent(in) :: cldn(mgncol,nlev)      ! cloud fraction (no units)
  real(r8), intent(in) :: liqcldf(mgncol,nlev)   ! liquid cloud fraction (no units)
  real(r8), intent(in) :: icecldf(mgncol,nlev)   ! ice cloud fraction (no units)


  real(r8), intent(in) :: delta_cf(mgncol,nlev)
  real(r8), intent(inout) :: D_eros_l(mgncol,nlev)
  real(r8), intent(inout) :: nerosc(mgncol,nlev)
  real(r8), intent(inout) :: D_eros_i(mgncol,nlev)
  real(r8), intent(inout) :: nerosi(mgncol,nlev)
  real(r8), intent(inout) :: dqcdt(mgncol,nlev)
  real(r8), intent(inout) :: dqidt(mgncol,nlev)

  ! used for scavenging
  ! Inputs for aerosol activation
  real(r8), intent(inout) :: naai(mgncol,nlev)     ! ice nucleation number (from microp_aero_ts) (1/kg)
  real(r8), intent(in)    :: npccn(mgncol,nlev)    ! ccn activated number tendency (from microp_aero_ts) (1/kg*s)

  ! Note that for these variables, the dust bin is assumed to be the last index.
  ! (For example, in CAM, the last dimension is always size 4.)
  real(r8), intent(in) :: rndst(:,:,:)  ! radius of each dust bin, for contact freezing (from microp_aero_ts) (m)
  real(r8), intent(in) :: nacon(:,:,:) ! number in each dust bin, for contact freezing  (from microp_aero_ts) (1/m^3)

  ! output arguments
  ! direct cw to precip conversion
  real(r8), intent(out) :: tlat(mgncol,nlev)         ! latent heating rate       (W/kg)
  real(r8), intent(out) :: qvlat(mgncol,nlev)        ! microphysical tendency qv (1/s)
  real(r8), intent(out) :: qctend(mgncol,nlev)       ! microphysical tendency qc (1/s)
  real(r8), intent(out) :: qitend(mgncol,nlev)       ! microphysical tendency qi (1/s)
  real(r8), intent(out) :: nctend(mgncol,nlev)       ! microphysical tendency nc (1/(kg*s))
  real(r8), intent(out) :: nitend(mgncol,nlev)       ! microphysical tendency ni (1/(kg*s))

  real(r8), intent(out) :: qrtend(mgncol,nlev)       ! microphysical tendency qr (1/s)
  real(r8), intent(out) :: qstend(mgncol,nlev)       ! microphysical tendency qs (1/s)
  real(r8), intent(out) :: nrtend(mgncol,nlev)       ! microphysical tendency nr (1/(kg*s))
  real(r8), intent(out) :: nstend(mgncol,nlev)       ! microphysical tendency ns (1/(kg*s))


  real(r8), intent(out) :: prect(mgncol)          ! surface precip rate (m/s)
  real(r8), intent(out) :: preci(mgncol)          ! cloud ice/snow precip rate (m/s)


  real(r8), intent(out) :: qsout(mgncol,nlev)        ! snow mixing ratio (kg/kg)
  real(r8), intent(out) :: rflx(mgncol,nlev+1)         ! grid-box average rain flux (kg m^-2 s^-1)
  real(r8), intent(out) :: sflx(mgncol,nlev+1)         ! grid-box average snow flux (kg m^-2 s^-1)
  real(r8), intent(out) :: qrout(mgncol,nlev)        ! grid-box average rain mixing ratio (kg/kg)
  real(r8), intent(out) :: reff_rain(mgncol,nlev)    ! rain effective radius (micron)
  real(r8), intent(out) :: reff_snow(mgncol,nlev)    ! snow effective radius (micron)

  character(128),   intent(out) :: errstring        ! output status (non-blank for error return)
  real(r8), intent(out) :: f_snow_berg  (mgncol,nlev) ! ratio of bergeron
                                                     ! production of qi to
                                                     ! sum of bergeron,
                                                     ! riming and freezing
  real(r8), intent(out) :: ssat_disposal(mgncol,nlev)
                                 ! disposition of supersaturation at end
                                 ! of step; 0.= no ssat, 1.= liq, 2.=ice)
  INTEGER,INTENT(IN) :: n_diag_4d
  REAL, dimension( mgncol,jdim, nlev, 0:n_diag_4d ), INTENT(INOUT) ::  diag_4l
  TYPE(diag_id_type),INTENT(IN) :: diag_id
  TYPE(diag_pt_type),INTENT(INout) :: diag_pt


!--> h1g, 2019-12-05
  ! temporary variables for sub-stepping
  real(r8) :: tlat1(mgncol,nlev)
  real(r8) :: qvlat1(mgncol,nlev)
  real(r8) :: qctend1(mgncol,nlev)
  real(r8) :: qitend1(mgncol,nlev)
  real(r8) :: nctend1(mgncol,nlev)
  real(r8) :: nitend1(mgncol,nlev)
  real(r8) :: qrtend1(mgncol,nlev)
  real(r8) :: qstend1(mgncol,nlev)
  real(r8) :: nrtend1(mgncol,nlev)
  real(r8) :: nstend1(mgncol,nlev)

  real(r8) :: prect1(mgncol)
  real(r8) :: preci1(mgncol)
  
!<-- h1g, 2019-12-05

  ! local workspace
  ! all units mks unless otherwise stated

  ! local copies of input variables
  real(r8) :: q(mgncol,nlev)           ! water vapor mixing ratio (kg/kg)
  real(r8) :: t(mgncol,nlev)           ! temperature (K)
  real(r8) :: qc(mgncol,nlev)      ! cloud liquid mixing ratio (kg/kg)
  real(r8) :: qi(mgncol,nlev)      ! cloud ice mixing ratio (kg/kg)
  real(r8) :: nc(mgncol,nlev)      ! cloud liquid number concentration (1/kg)
  real(r8) :: ni(mgncol,nlev)      ! cloud liquid number concentration (1/kg)
  real(r8) :: qr(mgncol,nlev)      ! rain mixing ratio (kg/kg)
  real(r8) :: qs(mgncol,nlev)      ! snow mixing ratio (kg/kg)
  real(r8) :: nr(mgncol,nlev)      ! rain number concentration (1/kg)
  real(r8) :: ns(mgncol,nlev)      ! snow number concentration (1/kg)

!--> h1g, 2019-12-06
! temporary copies of snow and rain variables
  real(r8) :: qrtmp(mgncol,nlev)      ! rain mixing ratio (kg/kg)
  real(r8) :: qstmp(mgncol,nlev)      ! snow mixing ratio (kg/kg)
  real(r8) :: nrtmp(mgncol,nlev)      ! rain number concentration (1/kg)
  real(r8) :: nstmp(mgncol,nlev)      ! snow number concentration (1/kg)
!<-- h1g, 2019-12-06

  ! general purpose variables
  real(r8) :: deltat            ! sub-time step (s)
  real(r8) :: mtime             ! the assumed ice nucleation timescale

  ! physical properties of the air at a given point
  real(r8) :: rho(mgncol,nlev)    ! density (kg m-3)
  real(r8) :: dv(mgncol,nlev)     ! diffusivity of water vapor
  real(r8) :: mu(mgncol,nlev)     ! viscosity
  real(r8) :: sc(mgncol,nlev)     ! schmidt number
  real(r8) :: rhof(mgncol,nlev)   ! density correction factor for fallspeed


  ! cloud fractions
  real(r8) :: precip_frac(mgncol,nlev) ! precip fraction assuming maximum overlap
  real(r8) :: cldm(mgncol,nlev)   ! cloud fraction
  real(r8) :: icldm(mgncol,nlev)  ! ice cloud fraction
  real(r8) :: lcldm(mgncol,nlev)  ! liq cloud fraction
  real(r8) :: qsfm(mgncol,nlev)   ! subgrid cloud water saturation scaling factor

  ! mass mixing ratios
  real(r8) :: qcic(mgncol,nlev)   ! in-cloud cloud liquid
  real(r8) :: qiic(mgncol,nlev)   ! in-cloud cloud ice
  real(r8) :: qsic(mgncol,nlev)   ! in-precip snow
  real(r8) :: qric(mgncol,nlev)   ! in-precip rain

  ! number concentrations
  real(r8) :: ncic(mgncol,nlev)   ! in-cloud droplet
  real(r8) :: niic(mgncol,nlev)   ! in-cloud cloud ice
  real(r8) :: nsic(mgncol,nlev)   ! in-precip snow
  real(r8) :: nric(mgncol,nlev)   ! in-precip rain
  ! maximum allowed ni value
  real(r8) :: nimax(mgncol,nlev)

  ! Size distribution parameters for:
  ! cloud ice
  real(r8) :: lami(mgncol,nlev)   ! slope
  real(r8) :: n0i(mgncol,nlev)    ! intercept
  ! cloud liquid
  real(r8) :: lamc(mgncol,nlev)   ! slope
  real(r8) :: pgam(mgncol,nlev)   ! spectral width parameter
  ! snow
  real(r8) :: lams(mgncol,nlev)   ! slope
  real(r8) :: n0s(mgncol,nlev)    ! intercept
  ! rain
  real(r8) :: lamr(mgncol,nlev)   ! slope
  real(r8) :: n0r(mgncol,nlev)    ! intercept

  ! Rates/tendencies due to:

  ! Instantaneous snow melting
  real(r8) :: minstsm(mgncol,nlev)    ! mass mixing ratio
  real(r8) :: ninstsm(mgncol,nlev)    ! number concentration
  ! Instantaneous rain freezing
  real(r8) :: minstrf(mgncol,nlev)    ! mass mixing ratio
  real(r8) :: ninstrf(mgncol,nlev)    ! number concentration

  ! deposition of cloud ice
  real(r8) :: vap_dep(mgncol,nlev)    ! deposition from vapor to ice PMC 12/3/12
  ! sublimation of cloud ice
  real(r8) :: ice_sublim(mgncol,nlev) ! sublimation from ice to vapor PMC 12/3/12
  ! ice nucleation
  real(r8) :: nnuccd(mgncol,nlev) ! number rate from deposition/cond.-freezing
  real(r8) :: mnuccd(mgncol,nlev) ! mass mixing ratio
  ! freezing of cloud water
  real(r8) :: mnuccc(mgncol,nlev) ! mass mixing ratio
  real(r8) :: nnuccc(mgncol,nlev) ! number concentration
  ! contact freezing of cloud water
  real(r8) :: mnucct(mgncol,nlev) ! mass mixing ratio
  real(r8) :: nnucct(mgncol,nlev) ! number concentration
  ! deposition nucleation in mixed-phase clouds (from external scheme)
  real(r8) :: mnudep(mgncol,nlev) ! mass mixing ratio
  real(r8) :: nnudep(mgncol,nlev) ! number concentration
  ! ice multiplication
  real(r8) :: msacwi(mgncol,nlev) ! mass mixing ratio
  real(r8) :: nsacwi(mgncol,nlev) ! number concentration
  ! autoconversion of cloud droplets
  real(r8) :: prc(mgncol,nlev)    ! mass mixing ratio
  real(r8) :: nprc(mgncol,nlev)   ! number concentration (rain)
  real(r8) :: nprc1(mgncol,nlev)  ! number concentration (cloud droplets)
  ! self-aggregation of snow
  real(r8) :: nsagg(mgncol,nlev)  ! number concentration
  ! self-collection of rain
  real(r8) :: nragg(mgncol,nlev)  ! number concentration
  ! collection of droplets by snow
  real(r8) :: psacws(mgncol,nlev)     ! mass mixing ratio
  real(r8) :: npsacws(mgncol,nlev)    ! number concentration
  ! collection of rain by snow
  real(r8) :: pracs(mgncol,nlev)  ! mass mixing ratio
  real(r8) :: npracs(mgncol,nlev) ! number concentration
  ! freezing of rain
  real(r8) :: mnuccr(mgncol,nlev) ! mass mixing ratio
  real(r8) :: nnuccr(mgncol,nlev) ! number concentration
  ! freezing of rain to form ice (mg add 4/26/13)
  real(r8) :: mnuccri(mgncol,nlev)    ! mass mixing ratio
  real(r8) :: nnuccri(mgncol,nlev)    ! number concentration
  ! accretion of droplets by rain
  real(r8) :: pra(mgncol,nlev)    ! mass mixing ratio
  real(r8) :: npra(mgncol,nlev)   ! number concentration
  ! autoconversion of cloud ice to snow
  real(r8) :: prci(mgncol,nlev)   ! mass mixing ratio
  real(r8) :: nprci(mgncol,nlev)  ! number concentration
  ! accretion of cloud ice by snow
  real(r8) :: prai(mgncol,nlev)   ! mass mixing ratio
  real(r8) :: nprai(mgncol,nlev)  ! number concentration
  ! evaporation of rain
  real(r8) :: pre(mgncol,nlev)    ! mass mixing ratio
  ! sublimation of snow
  real(r8) :: prds(mgncol,nlev)   ! mass mixing ratio
  ! number evaporation
  real(r8) :: nsubi(mgncol,nlev)  ! cloud ice
  real(r8) :: nsubc(mgncol,nlev)  ! droplet
  real(r8) :: nsubs(mgncol,nlev)  ! snow
  real(r8) :: nsubr(mgncol,nlev)  ! rain
  ! bergeron process
  real(r8) :: berg(mgncol,nlev)   ! mass mixing ratio (cloud ice)
  real(r8) :: bergs(mgncol,nlev)  ! mass mixing ratio (snow)


  ! fallspeeds
  ! number-weighted
  real(r8) :: uns(mgncol,nlev)    ! snow
  real(r8) :: unr(mgncol,nlev)    ! rain
  ! air density corrected fallspeed parameters
  real(r8) :: arn(mgncol,nlev)    ! rain
  real(r8) :: asn(mgncol,nlev)    ! snow
  real(r8) :: acn(mgncol,nlev)    ! cloud droplet
  real(r8) :: ain(mgncol,nlev)    ! cloud ice
  real(r8) :: ajn(mgncol,nlev)    ! cloud small ice

  ! Mass of liquid droplets used with external heterogeneous freezing.
  real(r8) :: mi0l(mgncol)

  ! saturation vapor pressures
  real(r8) :: esl(mgncol,nlev)    ! liquid
  real(r8) :: esi(mgncol,nlev)    ! ice
  real(r8) :: esn                 ! checking for RH after rain evap

  ! saturation vapor mixing ratios
  real(r8) :: qvl(mgncol,nlev)    ! liquid
  real(r8) :: qvi(mgncol,nlev)    ! ice
  real(r8) :: qvn                 ! checking for RH after rain evap

  ! relative humidity
  real(r8) :: relhum(mgncol,nlev)

  ! parameters for cloud water and cloud ice sedimentation calculations
  real(r8) :: fc(mgncol,nlev)
  real(r8) :: fnc(mgncol,nlev)
  real(r8) :: fi(mgncol,nlev)
  real(r8) :: fni(mgncol,nlev)

  real(r8) :: fr(mgncol,nlev)
  real(r8) :: fnr(mgncol,nlev)
  real(r8) :: fs(mgncol,nlev)
  real(r8) :: fns(mgncol,nlev)

  real(r8) :: faloutc(nlev)
  real(r8) :: faloutnc(nlev)
  real(r8) :: falouti(nlev)
  real(r8) :: faloutni(nlev)

  real(r8) :: faloutr(nlev)
  real(r8) :: faloutnr(nlev)
  real(r8) :: falouts(nlev)
  real(r8) :: faloutns(nlev)

  real(r8) :: faltndc
  real(r8) :: faltndnc
  real(r8) :: faltndi
  real(r8) :: faltndni
  real(r8) :: faltndqie
  real(r8) :: faltndqce

  real(r8) :: faltndr
  real(r8) :: faltndnr
  real(r8) :: faltnds
  real(r8) :: faltndns

  real(r8) :: rainrt(mgncol,nlev)     ! rain rate for reflectivity calculation

  ! dummy variables
  real(r8) :: dum
  real(r8) :: dum1
  real(r8) :: dum2
  real(r8) :: dumni0
  real(r8) :: dumns0
  ! dummies for checking RH
  real(r8) :: qtmp
  real(r8) :: ttmp
  ! dummies for conservation check
  real(r8) :: ratio
  real(r8) :: tmpfrz
  ! dummies for in-cloud variables
  real(r8) :: dumc(mgncol,nlev)   ! qc
  real(r8) :: dumnc(mgncol,nlev)  ! nc
  real(r8) :: dumi(mgncol,nlev)   ! qi
  real(r8) :: dumni(mgncol,nlev)  ! ni
  real(r8) :: dumr(mgncol,nlev)   ! rain mixing ratio
  real(r8) :: dumnr(mgncol,nlev)  ! rain number concentration
  real(r8) :: dums(mgncol,nlev)   ! snow mixing ratio
  real(r8) :: dumns(mgncol,nlev)  ! snow number concentration
  ! Array dummy variable
  real(r8) :: dum_2D(mgncol,nlev)
  real(r8) :: pdel_inv(mgncol,nlev)

  ! loop array variables
  ! "i" and "k" are column/level iterators for internal (MG) variables
  ! "n" is used for other looping (currently just sedimentation)
  integer i, k, n, it

  ! number of sub-steps for loops over "n" (for sedimentation)
  integer nstep
  integer mdust

  ! Varaibles to scale fall velocity between small and regular ice regimes.
  real(r8) :: irad
  real(r8) :: ifrac

  real(r8) :: tnd_qsnow(mgncol,nlev)    ! snow mass tendency (kg/kg/s)
  real(r8) :: tnd_nsnow(mgncol,nlev)    ! snow number tendency (#/kg/s)
  real(r8) :: re_ice(mgncol,nlev)       ! ice effective radius (m)

  ! From external ice nucleation.
  real(r8) :: frzimm(mgncol,nlev)       ! Number tendency due to immersion freezing (1/cm3)
  real(r8) :: frzcnt(mgncol,nlev)       ! Number tendency due to contact freezing (1/cm3)
  real(r8) :: frzdep(mgncol,nlev)       ! Number tendency due to deposition nucleation (1/cm3)

  real(r8) :: lflx(mgncol,nlev+1)       ! grid-box average liquid condensate flux (kg m^-2 s^-1)
  real(r8) :: iflx(mgncol,nlev+1)       ! grid-box average ice condensate flux (kg m^-2 s^-1)

  real(r8) :: qcsinksum_rate1ord(mgncol,nlev) ! 1st order rate for direct cw to precip conversion

  real(r8) :: effc(mgncol,nlev)         ! droplet effective radius (micron)
  real(r8) :: effc_fn(mgncol,nlev)      ! droplet effective radius, assuming nc = 1.e8 kg-1
  real(r8) :: effi(mgncol,nlev)         ! cloud ice effective radius (micron)
  real(r8) :: sadice(mgncol,nlev)       ! cloud ice surface area density (cm2/cm3)
  real(r8) :: sadsnow(mgncol,nlev)      ! cloud snow surface area density (cm2/cm3)

  real(r8) :: nevapr(mgncol,nlev)       ! evaporation rate of rain + snow (1/s)
  real(r8) :: evapsnow(mgncol,nlev)     ! sublimation rate of snow (1/s)
  real(r8) :: am_evp_st(mgncol,nlev)    ! stratiform evaporation area (frac)
  real(r8) :: prain(mgncol,nlev)        ! production of rain + snow (1/s)
  real(r8) :: prodsnow(mgncol,nlev)     ! production of snow (1/s)
  real(r8) :: cmeout(mgncol,nlev)       ! evap/sub of cloud (1/s)
  real(r8) :: deffi(mgncol,nlev)        ! ice effective diameter for optics (radiation) (micron)
  real(r8) :: pgamrad(mgncol,nlev)      ! ice gamma parameter for optics (radiation) (no units)
  real(r8) :: lamcrad(mgncol,nlev)      ! slope of droplet distribution for optics (radiation) (1/m)

  real(r8) :: dsout(mgncol,nlev)        ! snow diameter (m)

  real(r8) :: qcsevap(mgncol,nlev)      ! cloud water evaporation due to sedimentation (1/s)
  real(r8) :: qisevap(mgncol,nlev)      ! cloud ice sublimation due to sublimation (1/s)
  real(r8) :: qvres(mgncol,nlev)        ! residual condensation term to ensure RH < 100% (1/s)

  real(r8) :: cmeitot(mgncol,nlev)     ! grid-mean cloud ice sub/dep (1/s)
  real(r8) :: vtrmc(mgncol,nlev)        ! mass-weighted cloud water fallspeed (m/s)
  real(r8) :: vtrmi(mgncol,nlev)        ! mass-weighted cloud ice fallspeed (m/s)
  real(r8) :: umr(mgncol,nlev)          ! mass weighted rain fallspeed (m/s)
  real(r8) :: ums(mgncol,nlev)          ! mass weighted snow fallspeed (m/s)
  real(r8) :: qcsedten(mgncol,nlev)     ! qc sedimentation tendency (1/s)
  real(r8) :: qisedten(mgncol,nlev)     ! qi sedimentation tendency (1/s)
  real(r8) :: qrsedten(mgncol,nlev)     ! qr sedimentation tendency (1/s)
  real(r8) :: qssedten(mgncol,nlev)     ! qs sedimentation tendency (1/s)

  real(r8) :: pratot(mgncol,nlev)          ! accretion of cloud by rain
  real(r8) :: prctot(mgncol,nlev)          ! autoconversion of cloud to rain
  real(r8) :: mnuccctot(mgncol,nlev)       ! mixing ratio tend due to immersion freezing
  real(r8) :: mnuccttot(mgncol,nlev)       ! mixing ratio tend due to contact freezing
  real(r8) :: msacwitot(mgncol,nlev)       ! mixing ratio tend due to H-M splintering
  real(r8) :: psacwstot(mgncol,nlev)       ! collection of cloud water by snow
  real(r8) :: bergstot(mgncol,nlev)        ! bergeron process on snow
  real(r8) :: bergtot(mgncol,nlev)         ! bergeron process on cloud ice
  real(r8) :: melttot(mgncol,nlev)         ! melting of cloud ice
  real(r8) :: homotot(mgncol,nlev)         ! homogeneous freezing cloud water

  real(r8) :: qcrestot(mgncol,nlev)        ! residual cloud condensation due to removal of excess supersat
  real(r8) :: prcitot(mgncol,nlev)         ! autoconversion of cloud ice to snow
  real(r8) :: praitot(mgncol,nlev)         ! accretion of cloud ice by snow
  real(r8) :: qirestot(mgncol,nlev)        ! residual ice deposition due to removal of excess supersat
  real(r8) :: mnuccrtot(mgncol,nlev)       ! mixing ratio tendency due to heterogeneous freezing of rain to snow (1/s1)
  real(r8) :: mnuccritot(mgncol,nlev)      ! mixing ratio tendency due to heterogeneous freezing of rain to ice for small rain drops (1/s1)
  real(r8) :: pracstot(mgncol,nlev)        ! mixing ratio tendency due to accretion of rain by snow (1/s)
  real(r8) :: mnuccdtot(mgncol,nlev)      ! mass tendency from ice nucleation
  real(r8) :: meltsdttot(mgncol,nlev)      ! latent heating rate due to melting of snow  (W/kg)
  real(r8) :: frzrdttot(mgncol,nlev)       ! latent heating rate due to homogeneous freezing of rain (W/kg)
  real(r8) :: nrout(mgncol,nlev)        ! rain number concentration (1/m3)
  real(r8) :: nsout(mgncol,nlev)        ! snow number concentration (1/m3)
  real(r8) :: refl(mgncol,nlev)         ! analytic radar reflectivity
  real(r8) :: arefl(mgncol,nlev)        ! average reflectivity will zero points outside valid range
  real(r8) :: areflz(mgncol,nlev)       ! average reflectivity in z.
  real(r8) :: frefl(mgncol,nlev)        ! fractional occurrence of radar reflectivity
  real(r8) :: csrfl(mgncol,nlev)        ! cloudsat reflectivity
  real(r8) :: acsrfl(mgncol,nlev)       ! cloudsat average
  real(r8) :: fcsrfl(mgncol,nlev)       ! cloudsat fractional occurrence of radar reflectivity
  real(r8) :: rercld(mgncol,nlev)       ! effective radius calculation for rain + cloud
  real(r8) :: ncai(mgncol,nlev)         ! output number conc of ice nuclei available (1/m3)
  real(r8) :: ncal(mgncol,nlev)         ! output number conc of CCN (1/m3)
  real(r8) :: qrout2(mgncol,nlev)       ! copy of qrout as used to compute drout2
  real(r8) :: qsout2(mgncol,nlev)       ! copy of qsout as used to compute dsout2
  real(r8) :: nrout2(mgncol,nlev)       ! copy of nrout as used to compute drout2
  real(r8) :: nsout2(mgncol,nlev)       ! copy of nsout as used to compute dsout2
  real(r8) :: drout2(mgncol,nlev)       ! mean rain particle diameter (m)
  real(r8) :: dsout2(mgncol,nlev)       ! mean snow particle diameter (m)
  real(r8) :: freqs(mgncol,nlev)        ! fractional occurrence of snow
  real(r8) :: freqr(mgncol,nlev)        ! fractional occurrence of rain
  real(r8) :: nfice(mgncol,nlev)        ! fractional occurrence of ice
  real(r8) :: qcrat(mgncol,nlev)        ! limiter for qc process rates (1=no limit --> 0. no qc)

  real(r8) :: sum_freeze(mgncol,nlev)
  real(r8) :: sum_freeze2(mgncol,nlev)
  real(r8) :: sum_rime(mgncol,nlev)
  real(r8) :: sum_splinter(mgncol,nlev)
  real(r8) :: sum_bergs(mgncol,nlev)
  real(r8) :: sum_cond(mgncol,nlev)
  real(r8) :: sum_berg(mgncol,nlev)
  real(r8) :: sum_ice_adj(mgncol,nlev)
  real(r8) :: qldt_sum

!  these variables are only used in the GFDL implementation

  real(r8) :: cmelo(mgncol,nlev)   ! liquid condensation
  real(r8) :: eroslo(mgncol,nlev)  ! liquid erosion
  real(r8) :: erosio(mgncol,nlev)  ! ice erosion
  real(r8) :: preo(mgncol,nlev)    ! rain evaporation
  real(r8) :: prdso(mgncol,nlev)   ! snow sublimation
  real(r8) :: npccn2(mgncol,nlev)   ! ccn activated number tendency (from microp_aero_ts) (1/kg*s)
  logical  :: do_berg1
  logical  :: limit_berg = .false.
  real(r8) :: berg_lim = 1.0e-6_r8
  real(r8) :: dum3                ! temporary dummy variable


! droplet number
   real(r8) :: nucclim(mgncol,nlev)
   real(r8) :: nucclimo(mgncol,nlev)
  real(r8) :: npccno(mgncol,nlev)
   real(r8) :: nnuccco(mgncol,nlev)
  real(r8) :: nnuccto(mgncol,nlev)
  real(r8) :: npsacwso(mgncol,nlev)
   real(r8) :: nsubco(mgncol,nlev)
   real(r8) :: nerosco(mgncol,nlev)
   real(r8) :: nprao(mgncol,nlev)
   real(r8) :: nprc1o(mgncol,nlev)

! cloud ice number
  real(r8) :: nucclim1i(mgncol,nlev)
  real(r8) :: nucclim1io(mgncol,nlev)
  real(r8) :: nnuccdo(mgncol,nlev)
  real(r8) :: nsacwio(mgncol,nlev)
  real(r8) :: nsubio(mgncol,nlev)
  real(r8) :: nerosio(mgncol,nlev)
  real(r8) :: nprcio(mgncol,nlev)
  real(r8) :: npraio(mgncol,nlev)
  real(r8) :: nnuccrio(mgncol,nlev)

  real(r8) :: dum2i(mgncol,nlev)   ! used with ice nuleation
  real(r8) :: dum2l(mgncol,nlev)   ! used with drop nuleation
  real(r8) :: dum2a(mgncol,nlev)   ! used with ice nuleation

  real(r8) ::  flx(nlev), precip, dum_1D(nlev)

  !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

  ! default return error message
  errstring = ' '

  ! Process inputs

  ! assign variable deltat for sub-stepping...
  deltat  = deltatin


  ! Copies of input concentrations that may be changed internally.

  t  = tn
  q  = qn
  qc = qcn
  nc = ncn
  qi = qin
  ni = nin
  qr = qrn
  nr = nrn
  qs = qsn
  ns = nsn

  frzimm = 0.0_r8
  frzcnt = 0.0_r8
  frzdep = 0.0_r8

  tnd_qsnow = 0.0_r8
  tnd_nsnow = 0.0_r8
  re_ice    = 20e-6

  ssat_disposal = 0.0_r8
  sum_freeze    = 0.0_r8
  sum_freeze2   = 0.0_r8

  ! cldn: used to set cldm, unused for subcolumns
  ! liqcldf: used to set lcldm, unused for subcolumns
  ! icecldf: used to set icldm, unused for subcolumns

  if (microp_uniform) then
     ! subcolumns, set cloud fraction variables to one
     ! if cloud water or ice is present, if not present
     ! set to mincld (mincld used instead of zero, to prevent
     ! possible division by zero errors).

     where (qc >= qsmall)
        lcldm = 1._r8
     elsewhere
        lcldm = mincld
     end where

     where (qi >= qsmall)
        icldm = 1._r8
     elsewhere
        icldm = mincld
     end where

     cldm = max(icldm, lcldm)
     qsfm = 1._r8

  else
    ! get cloud fraction, check for minimum
     cldm = max(cldn,mincld)
     lcldm = max(liqcldf,mincld)
     icldm = max(icecldf,mincld)
     qsfm = 1.0_r8
  end if

  ! Initialize local variables

  ! local physical properties
  rho = p/(r*t)
  dv = 8.794E-5_r8 * t**1.81_r8 / p
  mu = 1.496E-6_r8 * t**1.5_r8 / (t + 120._r8)
  sc = mu/(rho*dv)


  ! air density adjustment for fallspeed parameters
  ! includes air density correction factor to the
  ! power of 0.54 following Heymsfield and Bansemer 2007

  rhof=(rhosu/rho)**0.54_r8
  if (.not. rho_factor_in_max_vt) rhof = 1.0
  rhof = MIN (rhof, max_rho_factor_in_vt)

! --->h1g, add namelist variables, 2014-07-01
! Zhao et al., ACP 2013, Table 1,
! ai: 350-1400 (s^-1);     as: 5.86-23.44 (m^0.59 s^-1)
! IceFallFac: 0.5 -- 2;    SnowFallFac: 0.5 -- 2
  arn=ar* RainFallFac * rhof                          !h1g, 2023-06-12
  asn=as* SnowFallFac * rhof
  acn=g*rhow/(18._r8*mu)
  ain=ai * IceFallFac * (rhosu/rho)**0.35_r8
  ajn=aj * SmallIceFallFac * (rhosu/rho)**0.35_r8     !h1g, 2023-06-12

  diag_4l(:,j,:,diag_pt%qidt_tiny)  = 0.0
  diag_4l(:,j,:,diag_pt%qnidt_tiny) = 0.0
  diag_4l(:,j,:,diag_pt%qrdt_tiny)  = 0.0
  diag_4l(:,j,:,diag_pt%qnrdt_tiny) = 0.0

  diag_4l(:,j,:,diag_pt%qldt_tiny)  = 0.0
  diag_4l(:,j,:,diag_pt%qndt_tiny)  = 0.0
  diag_4l(:,j,:,diag_pt%qsdt_tiny)  = 0.0
  diag_4l(:,j,:,diag_pt%qnsdt_tiny) = 0.0


  !INITIALIZE STUFF FOR SUBSTEPPING
  !===============================================

!--> h1g, 2019-12-05
  ! get sub-step time step
  deltat=deltat/real(iter)
!<-- h1g, 2019-12-05

  ! set mtime here to avoid answer-changing
  mtime=deltat

  ! initialize tendencies to zero
!--> h1g, 2019-12-05
  tlat1 = 0._r8
  qvlat1 = 0._r8
  qctend1 = 0._r8
  qitend1 = 0._r8
  nctend1 = 0._r8
  nitend1 = 0._r8
  qrtend1 = 0._r8
  qstend1 = 0._r8
  nrtend1 = 0._r8
  nstend1 = 0._r8

  prect1 = 0._r8
  preci1 = 0._r8 
!<-- h1g, 2019-12-05

  ! initialize microphysics output
  qcsevap=0._r8
  qisevap=0._r8
  qvres  =0._r8
  cmeitot =0._r8
  vtrmc =0._r8
  vtrmi =0._r8
  qcsedten =0._r8
  qisedten =0._r8
  qrsedten =0._r8
  qssedten =0._r8

  diag_4l(:,j,:,diag_pt%qndt_sedi) = 0.0
  diag_4l(:,j,:,diag_pt%qnidt_sedi) = 0.0
  diag_4l(:,j,:,diag_pt%rain_num_sedi) = 0.0
  diag_4l(:,j,:,diag_pt%snow_num_sedi) = 0.0

!initialize gfdl-only arrays-- this may not be needed and will be reviewed later.
    preo =0._r8
    prdso=0._r8
    cmelo =0._r8
    eroslo =0._r8
    erosio =0._r8
!droplet number
    nucclimo   = 0._r8
    npccno     = 0._r8
    nnuccco    = 0._r8
    nnuccto    = 0._r8
    npsacwso   = 0._r8
    nsubco     = 0._r8
    nerosco    = 0._r8
    nprao      = 0._r8
    nprc1o     = 0._r8
!ice number
    nucclim1io = 0._r8
    nnuccdo    = 0._r8
    nsacwio    = 0._r8
    nsubio     = 0._r8
    nerosio    = 0._r8
    nprcio     = 0._r8
    npraio     = 0._r8
    nnuccrio   = 0._r8

  pratot=0._r8
  prctot=0._r8
  mnuccctot=0._r8
  mnuccttot=0._r8
  msacwitot=0._r8
  psacwstot=0._r8
  bergstot=0._r8
  bergtot=0._r8
  melttot=0._r8
  homotot=0._r8
  qcrestot=0._r8
  prcitot=0._r8
  praitot=0._r8
  qirestot=0._r8
  mnuccrtot=0._r8
  mnuccritot=0._r8
  pracstot=0._r8
  meltsdttot=0._r8
  frzrdttot=0._r8
  mnuccdtot=0._r8

  rflx=0._r8
  sflx=0._r8
  lflx=0._r8
  iflx=0._r8

  ! initialize precip output

  qrout=0._r8
  qsout=0._r8
  nrout=0._r8
  nsout=0._r8

  ! for refl calc
  rainrt = 0._r8

  ! initialize rain size
  rercld=0._r8

  qcsinksum_rate1ord = 0._r8

  ! initialize variables for trop_mozart
  nevapr = 0._r8
  evapsnow = 0._r8
  am_evp_st = 0._r8
  prain = 0._r8
  prodsnow = 0._r8
  cmeout = 0._r8

  precip_frac = mincld

  lamc=0._r8

  !*********DO SUBSTEPPING!***************
  !============================================
  substepping: do it=1,iter

     ! initialize sub-step microphysical tendencies

  tlat=0._r8
  qvlat=0._r8
  qctend=0._r8
  qitend=0._r8
  qstend = 0._r8
  qrtend = 0._r8
  nctend=0._r8
  nitend=0._r8
  nrtend = 0._r8
  nstend = 0._r8

  ! initialize in-cloud and in-precip quantities to zero
  qcic  = 0._r8
  qiic  = 0._r8
  qsic  = 0._r8
  qric  = 0._r8

  ncic  = 0._r8
  niic  = 0._r8
  nsic  = 0._r8
  nric  = 0._r8

  ! initialize precip at surface

  prect = 0._r8
  preci = 0._r8

  ! initialize precip fallspeeds to zero
  ums = 0._r8
  uns = 0._r8
  umr = 0._r8
  unr = 0._r8

  ! initialize limiter for output
  qcrat = 1._r8

  ! Many outputs have to be initialized here at the top to work around
  ! ifort problems, even if they are always overwritten later.
  effc = 10._r8
  lamcrad = 0._r8
  pgamrad = 0._r8
  effc_fn = 10._r8
  effi = 25._r8
  sadice = 0._r8
  sadsnow = 0._r8
  deffi = 50._r8

  qrout2 = 0._r8
  nrout2 = 0._r8
  drout2 = 0._r8
  qsout2 = 0._r8
  nsout2 = 0._r8
  dsout = 0._r8
  dsout2 = 0._r8

  freqr = 0._r8
  freqs = 0._r8

  reff_rain = 0._r8
  reff_snow = 0._r8

  refl = -9999._r8
  arefl = 0._r8
  areflz = 0._r8
  frefl = 0._r8
  csrfl = 0._r8
  acsrfl = 0._r8
  fcsrfl = 0._r8

  ncal = 0._r8
  ncai = 0._r8

  nfice = 0._r8

!--> h1g, 2019-12-06
! re-calculate saturation vapor pressure for liquid and ice
  do k = 1, nlev
    do i = 1, mgncol
      call compute_qs (t(i,k), p(i,k), qvl(i,k), q = q(i,k),   &
                       esat = esl(i,k), es_over_liq = .true.)
      call compute_qs (t(i,k), p(i,k), qvi(i,k), q = q(i,k),   &
                       esat = esi(i,k), es_over_liq_and_ice = .true.)

      if ( t(i,k) < tmelt) then
           ! Scale the water saturation values to reflect subgrid scale
           ! ice cloud fraction, where ice clouds begin forming at a
           ! gridbox average relative humidity of rhmini (not 1).
           !
           ! NOTE: For subcolumns and other non-subgrid clouds, qsfm willi
           ! be 1.
           qvi(i,k) = qsfm(i,k) * qvi(i,k)
           esi(i,k) = qsfm(i,k) * esi(i,k)
           qvl(i,k) = qsfm(i,k) * qvl(i,k)
           esl(i,k) = qsfm(i,k) * esl(i,k)
      endif
    end do
  end do

  relhum = q / max(qvl, qsmall)
!<-- h1g, 2019-12-06

! --->h1g, 2019-10-25
              dum_30 = 5.0_r8*exp(0.304_r8*( 30.0))
              dum_5  = 5.0_r8*exp(0.304_r8*( 5.0))
! <---h1g, 2019-10-25

! calculate ice nucleation dum2i
  do k=1,nlev
   do i=1,mgncol
          if (t(i,k).lt. icenuct       ) then

            if ( liu_in .or. use_Fan2019 ) then
              dum2i(i,k) = naai(i,k)
              if ( lat(i)*180.0/3.14159 > 60.0 ) dum2i(i,k) = dum2i(i,k) * ice_nucl_factor  ! h1g,  2019-12-17

! --->h1g, 2014-05-30 add Meyers ice nucleation formula (Only temperature dependent)
            elseif ( use_Meyers ) then
               dum2i(i,k) =  (exp(12.96* (tmelt -t(i,k))/80 - 0.639)) *1000._r8
               dum2i(i,k)= ( dum2i(i,k) )/rho(i,k) ! convert from m-3 to kg-1
! <--- h1g,  2014-05-30

! <--  h1g,  2020-04-18
            elseif ( use_FanAndCooper ) then
                dum2i(i,k) = 2.74 * concen_dust_sub(i,k) * exp(0.412_r8*(tmelt-t(i,k)))
                dum2i(i,k) = min( dum2i(i,k), Nice_max_Fan*1000._r8 )
                dum2i(i,k) = dum2i(i,k) * p(i,k)/ 95000.
                dum2i(i,k)= dum2i(i,k)/rho(i,k) ! convert from m-3 to kg-1

                dum_tmp = 0.005_r8*exp(0.304_r8*(tmelt-t(i,k)))*1000._r8
                dum_tmp = min( dum_tmp, 5.0_r8*exp(0.304_r8*(-tc_cooper)))
                dum_tmp = dum_tmp/rho(i,k) ! convert from m-3 to kg-1
                dum2i(i,k)= dum2i(i,k) + dum_tmp
! -->  h1g,  2020-04-18

            else
! cooper curve (factor of 1000 is to convert from L-1 to m-3)
              dum2i(i,k)=0.005_r8*exp(0.304_r8*(tmelt-t(i,k)))*1000._r8
! put limit on number of nucleated crystals, set to number at T=-30 C
! cooper (limit to value at -35 C)
              dum2i(i,k)= min(dum2i(i,k),5.0_r8*exp(0.304_r8*(-tc_cooper)))

! --->h1g, 2019-10-25
              dum_tmp = dum_5 + (dum_30-dum_5)/25.0*( tmelt-t(i,k) - 5.0)
              if ( do_Ni_linear_interp ) &
              dum2i(i,k)= max(dum2i(i,k), dum_tmp )
! <---h1g, 2019-10-25

              dum2i(i,k)=dum2i(i,k)/rho(i,k) ! convert from m-3 to kg-1
            endif
          else
            dum2i(i,k)=0._r8
          end if  ! t(i,k).lt. icenuct
         ! naai(i,k) = dum2i(i,k)
     end do
     end do

  dum2l = 0.

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! droplet activation
   do k=1,nlev
    do i=1,mgncol
     if ( qc(i,k).ge.qsmall ) then
       dum2l(i,k) = max(0._r8, npccn  (i,k))

!RSH npccn2 is the change in droplet number on this step. In the
! non-GFDL_COMPATIBLE_MICROP, this is not calculated since the input
! droplet number already includes the newly activated droplets.
        npccn2(i,k) = ((dum2l(i,k) - nc(i,k)/cldm(i,k))/tau_act_liq )*cldm(i,k)
        npccn2(i,k) = max(0._r8,npccn2(i,k))
      else
       npccn2(i,k)=0._r8
     end if
  !   nc(i,k)    = nc(i,k)+npccn2(i,k)*deltat   ! from MG2, additional changes from h1g
    end do
   end do

! ice activation
  if (do_cldice) then
! --->h1g, 2014-05-30 add Meyers ice nucleation option
      where (dum2i > 0._r8 .and. t < icenuct .and. &
        relhum*esl/esi > rhmini+0.05_r8 .and. icldm > icld_cri)

        !if NAAI > 0. then set numice = naai (as before)
        !note: this is gridbox averaged
        nnuccd = (dum2i-ni/icldm)/tau_act_ice * icldm
        nnuccd = max(nnuccd,0._r8)
        nimax = dum2i*icldm
     !Calc mass of new particles using new crystal mass...
     !also this will be multiplied by mtime as nnuccd is...
        mnuccd = nnuccd * mi0

      elsewhere
        nnuccd = 0._r8
        nimax = 0._r8
        mnuccd = 0._r8
      end where
    !  ni    = ni + nnuccd*deltat    ! h1g, 2019-12-10
  end if  ! do_cldice

  !=============================================================================
!--> h1g, 2019-12-06
! temporary copy of qs, ns, qr, nr for calculating qs, ns, qr, nr tendencies
  qstmp = qs
  nstmp = ns
  qrtmp = qr
  nrtmp = nr
!<-- h1g, 2019-12-06

  do k=1,nlev
     do i=1,mgncol

        ! calculate instantaneous precip processes (melting and homogeneous freezing)

        ! melting of snow at +2 C

        if (t(i,k) > snowmelt) then
           if (qs(i,k) > 0._r8) then

              ! make sure melting snow doesn't reduce temperature below threshold
              dum = -xlf/cpp*qs(i,k)
              if (t(i,k)+dum < snowmelt) then
                 dum = (t(i,k)-snowmelt)*cpp/xlf
                 dum = dum/qs(i,k)
                 dum = max(0._r8,dum)
                 dum = min(1._r8,dum)
              else
                 dum = 1._r8
              end if

              minstsm(i,k) = dum*qs(i,k)
              ninstsm(i,k) = dum*ns(i,k)

              dum1=-xlf*minstsm(i,k)/deltat
              tlat(i,k)=tlat(i,k)+dum1
              meltsdttot(i,k)=meltsdttot(i,k) + dum1

              qs(i,k) = max(qs(i,k) - minstsm(i,k), 0._r8)
              ns(i,k) = max(ns(i,k) - ninstsm(i,k), 0._r8)
              qr(i,k) = max(qr(i,k) + minstsm(i,k), 0._r8)
              nr(i,k) = max(nr(i,k) + ninstsm(i,k), 0._r8)
           end if
        end if

     end do
  end do


  do k=1,nlev
    do i=1,mgncol
        ! freezing of rain at -5 C

        if (t(i,k) < rainfrze) then

           if (qr(i,k) > 0._r8) then

              ! make sure freezing rain doesn't increase temperature above threshold
              dum = xlf/cpp*qr(i,k)
              if (t(i,k)+dum > rainfrze) then
                 dum = -(t(i,k)-rainfrze)*cpp/xlf
                 dum = dum/qr(i,k)
                 dum = max(0._r8,dum)
                 dum = min(1._r8,dum)
              else
                 dum = 1._r8
              end if

              minstrf(i,k) = dum*qr(i,k)
              ninstrf(i,k) = dum*nr(i,k)

              ! heating tendency
              dum1 = xlf*minstrf(i,k)/deltat
              tlat(i,k)=tlat(i,k)+dum1
              frzrdttot(i,k)=frzrdttot(i,k) + dum1

              qr(i,k) = max(qr(i,k) - minstrf(i,k), 0._r8)
              nr(i,k) = max(nr(i,k) - ninstrf(i,k), 0._r8)
              qs(i,k) = max(qs(i,k) + minstrf(i,k), 0._r8)
              ns(i,k) = max(ns(i,k) + ninstrf(i,k), 0._r8)

           end if
        end if
     end do
  end do


  do k=1,nlev
    do i=1,mgncol
        ! obtain in-cloud values of cloud water/ice mixing ratios and number concentrations
        !-------------------------------------------------------
        ! for microphysical process calculations
        ! units are kg/kg for mixing ratio, 1/kg for number conc

        if (qc(i,k).ge.qsmall) then
           ! limit in-cloud values to 0.005 kg/kg
           qcic(i,k)=min(qc(i,k)/lcldm(i,k),5.e-3_r8)
           ncic(i,k)=max(nc(i,k)/lcldm(i,k),0._r8)

           ! specify droplet concentration
           if (nccons) then
              ncic(i,k)=ncnst/rho(i,k)
           end if
        else
           qcic(i,k)=0._r8
           ncic(i,k)=0._r8
        end if

        if (qi(i,k).ge.qsmall) then
           ! limit in-cloud values to 0.005 kg/kg
           qiic(i,k)=min(qi(i,k)/icldm(i,k),5.e-3_r8)
           niic(i,k)=max( ( ni(i,k) +nnuccd(i,k)*deltat )/icldm(i,k),0._r8)  ! h1g, 2019-12-10
        !   niic(i,k)=max( ni(i,k)/icldm(i,k),0._r8)  ! h1g, 2019-12-10

           ! switch for specification of cloud ice number
           if (nicons) then
              niic(i,k)=ninst/rho(i,k)
           end if
        else
           qiic(i,k)=0._r8
           niic(i,k)=0._r8
        end if

     end do
  end do

  !========================================================================

  ! for sub-columns cldm has already been set to 1 if cloud
  ! water or ice is present, so precip_frac will be correctly set below
  ! and nothing extra needs to be done here

  precip_frac = cldm

  micro_vert_loop: do k=1,nlev

    if (trim(micro_mg_precip_frac_method) == 'in_cloud') then

      if (k /= 1) then
        where (qc(:,k) < qsmall .and. qi(:,k) < qsmall)
              precip_frac(:,k) = precip_frac(:,k-1)
        end where
      endif

    else if (trim(micro_mg_precip_frac_method) == 'max_overlap') then

    ! calculate precip fraction based on maximum overlap assumption

    ! if rain or snow mix ratios are smaller than threshold,
    ! then leave precip_frac as cloud fraction at current level
      if (k /= 1) then
        where (qr(:,k-1) >= qsmall .or. qs(:,k-1) >= qsmall)
              precip_frac(:,k)=max(precip_frac(:,k-1),precip_frac(:,k))
        end where
      end if

    endif


  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  ! get size distribution parameters based on in-cloud cloud water
  ! these calculations also ensure consistency between number and mixing ratio
  !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

  ! cloud liquid
  !-------------------------------------------
     call size_dist_param_liq(mg_liq_props, qcic(:,k), ncic(:,k), rho(:,k), &
          pgam(:,k), lamc(:,k), mgncol)

        !========================================================================
        ! autoconversion of cloud liquid water to rain
        ! formula from Khrouditnov and Kogan (2000), modified for sub-grid distribution of qc
        ! minimum qc of 1 x 10^-8 prevents floating point error

     if (.not. do_sb_physics .and.  .not.do_cotton_auto) then
       call kk2000_liq_autoconversion(microp_uniform, qcic(:,k), &
          ncic(:,k), rho(:,k), relvar(:,k), prc(:,k), nprc(:,k), nprc1(:,k), mgncol)
     endif

     ! assign qric based on prognostic qr, using assumed precip fraction
     ! note: this could be moved above for consistency with qcic and qiic calculations
     qric(:,k) = qr(:,k)/precip_frac(:,k)
     nric(:,k) = nr(:,k)/precip_frac(:,k)

     ! limit in-precip mixing ratios to 10 g/kg
     qric(:,k)=min(qric(:,k),0.01_r8)

     ! add autoconversion to precip from above to get provisional rain mixing ratio
     ! and number concentration (qric and nric)

     where (qric(:,k).lt.qsmall)
        qric(:,k)=0._r8
        nric(:,k)=0._r8
     end where

    ! make sure number concentration is a positive number to avoid
    ! taking root of negative later

    nric(:,k)=max(nric(:,k),0._r8)

    ! Get size distribution parameters for cloud ice

    call size_dist_param_basic(mg_ice_props, qiic(:,k), niic(:,k), &
          lami(:,k), mgncol, n0=n0i(:,k))

     ! Alternative autoconversion
     if (do_cotton_auto ) then
       call cotton_liq_autoconversion(.false., qcic(:,k), &
          ncic(:,k), rho(:,k), relvar(:,k), rthresh, prc(:,k), nprc(:,k), nprc1(:,k), mgncol)
     elseif (do_sb_physics) then
       call sb2001v2_liq_autoconversion(pgam(:,k),qcic(:,k),ncic(:,k), &
            qric(:,k),rho(:,k),relvar(:,k),prc(:,k),nprc(:,k),nprc1(:,k), mgncol)
     endif


    !.......................................................................
    ! Autoconversion of cloud ice to snow
    ! similar to Ferrier (1994)

    if (do_cldice) then
       call ice_autoconversion(t(:,k), qiic(:,k), lami(:,k), n0i(:,k), &
            dcs, prci(:,k), nprci(:,k), mgncol)

    else
       ! Add in the particles that we have already converted to snow, and
       ! don't do any further autoconversion of ice.
       prci(:,k)  = tnd_qsnow(:,k) / cldm(:,k)
       nprci(:,k) = tnd_nsnow(:,k) / cldm(:,k)
    end if

     ! note, currently we don't have this
     ! inside the do_cldice block, should be changed later
     ! assign qsic based on prognostic qs, using assumed precip fraction
     qsic(:,k) = qs(:,k)/precip_frac(:,k)
     nsic(:,k) = ns(:,k)/precip_frac(:,k)

     ! limit in-precip mixing ratios to 10 g/kg
     qsic(:,k)=min(qsic(:,k),0.01_r8)

     ! if precip mix ratio is zero so should number concentration

     where (qsic(:,k) < qsmall)
        qsic(:,k)=0._r8
        nsic(:,k)=0._r8
     end where

     ! make sure number concentration is a positive number to avoid
     ! taking root of negative later

     nsic(:,k)=max(nsic(:,k),0._r8)
     !.......................................................................
     ! get size distribution parameters for precip
     !......................................................................
     ! rain

     call size_dist_param_basic(mg_rain_props, qric(:,k), nric(:,k), &
          lamr(:,k), mgncol, n0=n0r(:,k))

     where (lamr(:,k) >= qsmall)

        ! provisional rain number and mass weighted mean fallspeed (m/s)

        unr(:,k) = min(arn(:,k)*gamma_br_plus1/lamr(:,k)**br,9.1_r8*rhof(:,k))
        umr(:,k) = min(arn(:,k)*gamma_br_plus4/(6._r8*lamr(:,k)**br),9.1_r8*rhof(:,k))

     elsewhere
        umr(:,k) = 0._r8
        unr(:,k) = 0._r8
     end where

     !......................................................................
     ! snow

     call size_dist_param_basic(mg_snow_props, qsic(:,k), nsic(:,k), &
          lams(:,k), mgncol, n0=n0s(:,k))

     where (lams(:,k) > 0._r8)
        ! provisional snow number and mass weighted mean fallspeed (m/s)
        ums(:,k) = min(asn(:,k)*gamma_bs_plus4/(6._r8*lams(:,k)**bs),1.2_r8*rhof(:,k))
        uns(:,k) = min(asn(:,k)*gamma_bs_plus1/lams(:,k)**bs,1.2_r8*rhof(:,k))

     elsewhere
        ums(:,k) = 0._r8
        uns(:,k) = 0._r8
     end where

     if (do_cldice) then
        if (.not. use_hetfrz_classnuc) then

           ! heterogeneous freezing of cloud water
           !----------------------------------------------

           call immersion_freezing(microp_uniform, t(:,k), pgam(:,k), lamc(:,k), &
                qcic(:,k), ncic(:,k), relvar(:,k), mnuccc(:,k), nnuccc(:,k), mgncol)


           ! make sure number of droplets frozen does not exceed available ice nuclei concentration
           ! this prevents 'runaway' droplet freezing

           where (qcic(:,k).ge.qsmall .and. t(:,k).lt.269.15_r8)
              where (nnuccc(:,k)*lcldm(:,k).gt.nnuccd(:,k))
                 ! scale mixing ratio of droplet freezing with limit
                 mnuccc(:,k)=mnuccc(:,k)*(nnuccd(:,k)/(nnuccc(:,k)*lcldm(:,k)))
                 nnuccc(:,k)=nnuccd(:,k)/lcldm(:,k)
              end where
           end where

           mdust = size(rndst,3)
           call contact_freezing(microp_uniform, t(:,k), p(:,k), rndst(:,k,:), &
                nacon(:,k,:), pgam(:,k), lamc(:,k), qcic(:,k), ncic(:,k), &
                relvar(:,k), mnucct(:,k), nnucct(:,k), mgncol, mdust)

           mnudep(:,k)=0._r8
           nnudep(:,k)=0._r8
        else

           ! Mass of droplets frozen is the average droplet mass, except
           ! with two limiters: concentration must be at least 1/cm^3, and
           ! mass must be at least the minimum defined above.
           mi0l = qcic(:,k)/max(ncic(:,k), 1.0e6_r8/rho(:,k))
           mi0l = max(mi0l_min, mi0l)

           where (qcic(:,k) >= qsmall)
              nnuccc(:,k) = frzimm(:,k)*1.0e6_r8/rho(:,k)
              mnuccc(:,k) = nnuccc(:,k)*mi0l

              nnucct(:,k) = frzcnt(:,k)*1.0e6_r8/rho(:,k)
              mnucct(:,k) = nnucct(:,k)*mi0l

              nnudep(:,k) = frzdep(:,k)*1.0e6_r8/rho(:,k)
              mnudep(:,k) = nnudep(:,k)*mi0
           elsewhere
              nnuccc(:,k) = 0._r8
              mnuccc(:,k) = 0._r8

              nnucct(:,k) = 0._r8
              mnucct(:,k) = 0._r8

              nnudep(:,k) = 0._r8
              mnudep(:,k) = 0._r8
           end where

        end if

     else
        mnuccc(:,k)=0._r8
        nnuccc(:,k)=0._r8
        mnucct(:,k)=0._r8
        nnucct(:,k)=0._r8
        mnudep(:,k)=0._r8
        nnudep(:,k)=0._r8
     end if

     call snow_self_aggregation(t(:,k), rho(:,k), asn(:,k), rhosn, qsic(:,k), nsic(:,k), &
          nsagg(:,k), mgncol)

     call accrete_cloud_water_snow(t(:,k), rho(:,k), asn(:,k), uns(:,k), mu(:,k), &
          qcic(:,k), ncic(:,k), qsic(:,k), pgam(:,k), lamc(:,k), lams(:,k), n0s(:,k), &
          psacws(:,k), npsacws(:,k), use_const_ELI, ELI_RK, mgncol)   ! h1g, 2020-04-16
      if ( .not. do_liq_num_riming ) npsacws(:,k) = 0.0    ! h1g, 2020-03-10


     if (do_cldice .and. do_HM_splinter ) then
        call secondary_ice_production(t(:,k), psacws(:,k), msacwi(:,k), nsacwi(:,k), mgncol)
     else
        nsacwi(:,k) = 0.0_r8
        msacwi(:,k) = 0.0_r8
     end if

     call accrete_rain_snow(t(:,k), rho(:,k), umr(:,k), ums(:,k), unr(:,k), uns(:,k), &
          qric(:,k), qsic(:,k), lamr(:,k), n0r(:,k), lams(:,k), n0s(:,k), &
          pracs(:,k), npracs(:,k), mgncol)

     call heterogeneous_rain_freezing(t(:,k), qric(:,k), nric(:,k), lamr(:,k), &
          mnuccr(:,k), nnuccr(:,k), mgncol)

     if (do_sb_physics) then
       call sb2001v2_accre_cld_water_rain(qcic(:,k), ncic(:,k), qric(:,k), &
            rho(:,k), relvar(:,k), pra(:,k), npra(:,k), mgncol)
     else
       call accrete_cloud_water_rain(microp_uniform, qric(:,k), qcic(:,k), &
            ncic(:,k), relvar(:,k), pra(:,k), npra(:,k), mgncol)
     endif

     if (.not. microp_uniform) then
       pra(:,k)  = pra(:,k)  * accre_enhan(:,k)
       npra(:,k) = npra(:,k) * accre_enhan(:,k)
     endif


     call self_collection_rain(rho(:,k), qric(:,k), nric(:,k), nragg(:,k), mgncol)

     if (do_cldice) then
        call accrete_cloud_ice_snow(t(:,k), rho(:,k), asn(:,k), qiic(:,k), niic(:,k), &
             qsic(:,k), lams(:,k), n0s(:,k), prai(:,k), nprai(:,k), mgncol)
     else
           prai(:,k) = 0._r8
           nprai(:,k) = 0._r8
     end if

     call evaporate_sublimate_precip(t(:,k), rho(:,k), &
          dv(:,k), mu(:,k), sc(:,k), q(:,k), qvl(:,k), qvi(:,k), &
          lcldm(:,k), precip_frac(:,k), arn(:,k), asn(:,k), qcic(:,k), qiic(:,k), &
          qric(:,k), qsic(:,k), lamr(:,k), n0r(:,k), lams(:,k), n0s(:,k), &
          pre(:,k), prds(:,k), am_evp_st(:,k), evap_RH_cri, subl_RH_cri, mgncol)
       pre(:,k)  = evap_subl_fac * pre(:,k)
       prds(:,k) = evap_subl_fac * prds(:,k)

!       where (relhum(:,k) > evap_RH_cri ) 
!         pre(:,k)  = 0.0
!      endwhere
!       where (relhum(:,k) > subl_RH_cri ) 
!         prds(:,k)  = 0.0
!      endwhere

     call bergeron_process_snow(t(:,k), rho(:,k), dv(:,k), mu(:,k), sc(:,k), &
          qvl(:,k), qvi(:,k), asn(:,k), qcic(:,k), qsic(:,k), lams(:,k), n0s(:,k), &
          bergs(:,k), mgncol)

     bergs(:,k)=bergs(:,k)*micro_mg_bergs_eff_factor

     !+++PMC 12/3/12 - NEW VAPOR DEP/SUBLIMATION GOES HERE!!!
     if (do_cldice) then

        call ice_deposition_sublimation(t(:,k), q(:,k), qc(:,k), qi(:,k), ni(:,k), &
             icldm(:,k), rho(:,k), dv(:,k), qvl(:,k), qvi(:,k), &
             berg(:,k), vap_dep(:,k), ice_sublim(:,k), mgncol)
        vap_dep(:,k)    = vap_dep(:,k)    * ice_sublim_factor  !h1g, 2017-07-17
        ice_sublim(:,k) = ice_sublim(:,k) * ice_sublim_factor  !h1g, 2017-07-16
        berg(:,k)=berg(:,k)*micro_mg_berg_eff_factor
      !  vap_dep(:,k) = vap_dep(:,k) + dqidt(:,k)    ! h1g, 2016-12-19

        where ( ice_sublim(:,k) < 0._r8 .and. qi(:,k) > qsmall .and. icldm(:,k) > mincld)
           nsubi(:,k) = sublim_factor*ice_sublim(:,k) / qi(:,k) * ni(:,k) / icldm(:,k)
        elsewhere
           nsubi(:,k) = 0._r8
        end where

        ! bergeron process should not reduce nc unless
        ! all ql is removed (which is handled elsewhere)
        !in fact, nothing in this entire file makes nsubc nonzero.
        nsubc(:,k) = 0._r8

     end if !do_cldice
     !---PMC 12/3/12

     do i=1,mgncol

        ! conservation to ensure no negative values of cloud water/precipitation
        ! in case microphysical process rates are large
        !===================================================================

        ! note: for check on conservation, processes are multiplied by omsm
        ! to prevent problems due to round off error

        ! conservation of qc
        !-------------------------------------------------------------------

        dum = ((prc(i,k)+pra(i,k)+mnuccc(i,k)+mnucct(i,k)+msacwi(i,k)+ &
             psacws(i,k)+bergs(i,k))*lcldm(i,k)+berg(i,k))*deltat

        if (dum.gt.qc(i,k)+dqcdt(i,k)*deltat+D_eros_l(i,k)*deltat) then

           ratio = (qc(i,k)/deltat+ max(dqcdt(i,k),0.0))/((prc(i,k)+pra(i,k)+mnuccc(i,k)+mnucct(i,k)+ &
                msacwi(i,k)+psacws(i,k)+bergs(i,k))*lcldm(i,k)+berg(i,k)-D_eros_l(i,k)+max(-dqcdt(i,k),0.0))*omsm
           prc(i,k) = prc(i,k)*ratio
           pra(i,k) = pra(i,k)*ratio
           mnuccc(i,k) = mnuccc(i,k)*ratio
           mnucct(i,k) = mnucct(i,k)*ratio
           msacwi(i,k) = msacwi(i,k)*ratio
           psacws(i,k) = psacws(i,k)*ratio
           bergs(i,k) = bergs(i,k)*ratio
           berg(i,k) = berg(i,k)*ratio
           D_eros_l(i,k) = D_eros_l(i,k)*ratio
           if( dqcdt(i,k) < 0.0 ) dqcdt(i,k)=dqcdt(i,k)*ratio
           qcrat(i,k) = ratio
           if( ratio > 1) print*,'error ratio>1 in conservation of qc', i,k,ratio
           if( ratio < 0) print*,'error ratio<0 in conservation of qc', i,k,ratio
         else
           qcrat(i,k) = 1._r8
        end if


        !PMC 12/3/12: ratio is also frac of step w/ liquid.
        !thus we apply berg for "ratio" of timestep and vapor
        !deposition for the remaining frac of the timestep.
        if (qc(i,k) >= qsmall) then
           vap_dep(i,k) = vap_dep(i,k)*(1._r8-qcrat(i,k))
        end if

     end do

     do i=1,mgncol

        !=================================================================
        ! apply limiter to ensure that ice/snow sublimation and rain evap
        ! don't push conditions into supersaturation, and ice deposition/nucleation don't
        ! push conditions into sub-saturation
        ! note this is done after qc conservation since we don't know how large
        ! vap_dep is before then
        ! estimates are only approximate since other process terms haven't been limited
        ! for conservation yet

        ! first limit ice deposition/nucleation vap_dep + mnuccd
        dum1 = vap_dep(i,k) + mnuccd(i,k)
        if (dum1 > 1.e-20_r8) then
           dum = (q(i,k)-qvi(i,k))/(1._r8 + xxls_squared*qvi(i,k)/(cpp*rv*t(i,k)**2))/deltat
           dum = max(dum,0._r8)
           if (dum1 > dum) then
              ! Allocate the limited "dum" tendency to mnuccd and vap_dep
              ! processes. Don't divide by cloud fraction; these are grid-
              ! mean rates.
              dum1 = mnuccd(i,k) / (vap_dep(i,k)+mnuccd(i,k))
              mnuccd(i,k) = dum*dum1
              vap_dep(i,k) = dum - mnuccd(i,k)
           end if
        end if

     end do

     do i=1,mgncol

        !===================================================================
        ! conservation of nc
        !-------------------------------------------------------------------
        dum = ((nprc1(i,k)+npra(i,k)+nnuccc(i,k)+nnucct(i,k)+ &
             npsacws(i,k)-nsubc(i,k))*lcldm(i,k) - npccn2(i,k) )*deltat

        if (dum.gt.nc(i,k)+nerosc(i,k)*lcldm(i,k)*deltat) then
           ratio = ( nc(i,k)/deltat + npccn2(i,k) )/( (nprc1(i,k)+npra(i,k)+nnuccc(i,k)+nnucct(i,k)+&
                npsacws(i,k)-nsubc(i,k)-nerosc(i,k) )*lcldm(i,k) )*omsm

           nprc1(i,k)   = nprc1(i,k)  * ratio
           npra(i,k)    = npra(i,k)   * ratio
           nnuccc(i,k)  = nnuccc(i,k) * ratio
           nnucct(i,k)  = nnucct(i,k) * ratio
           npsacws(i,k) = npsacws(i,k)* ratio
           nsubc(i,k)   = nsubc(i,k)  * ratio
           nerosc(i,k)  = nerosc(i,k) * ratio
        end if

        mnuccri(i,k)=0._r8
        nnuccri(i,k)=0._r8

        if (do_cldice) then

           ! freezing of rain to produce ice if mean rain size is smaller than Dcs
           if (lamr(i,k) > qsmall .and. 1._r8/lamr(i,k) < Dcs) then
              mnuccri(i,k)=mnuccr(i,k)
              nnuccri(i,k)=nnuccr(i,k)
              mnuccr(i,k)=0._r8
              nnuccr(i,k)=0._r8
           end if
        end if

     end do

     do i=1,mgncol
        ! conservation of rain mixing ratio
        !-------------------------------------------------------------------
        dum = ((-pre(i,k)+pracs(i,k)+mnuccr(i,k)+mnuccri(i,k))*precip_frac(i,k)- &
             (pra(i,k)+prc(i,k))*lcldm(i,k))*deltat

        ! note that qrtend is included below because of instantaneous freezing/melt
        if (dum.gt.qr(i,k).and. &
             (-pre(i,k)+pracs(i,k)+mnuccr(i,k)+mnuccri(i,k)).ge.qsmall) then
           ratio = (qr(i,k)/deltat+(pra(i,k)+prc(i,k))*lcldm(i,k))/   &
                precip_frac(i,k)/(-pre(i,k)+pracs(i,k)+mnuccr(i,k)+mnuccri(i,k))*omsm
           pre(i,k)=pre(i,k)*ratio
           pracs(i,k)=pracs(i,k)*ratio
           mnuccr(i,k)=mnuccr(i,k)*ratio
           mnuccri(i,k)=mnuccri(i,k)*ratio
        end if
     end do

     do i=1,mgncol
        ! conservation of rain number
        !-------------------------------------------------------------------
        ! Add evaporation of rain number.
        if (allow_rain_num_evap)  then
          if (pre(i,k) < 0._r8 .and. qr(i,k) > qsmall) then
            nsubr(i,k) = max(-nr(i,k)/deltat, pre(i,k)*nr(i,k)/qr(i,k))
          else
            nsubr(i,k) = 0._r8
          endif
        else
        ! neglect evaporation of nr, allow_rain_num_evap = .false. h1g, 2017-02-27
           nsubr(i,k) = 0._r8
        end if
     end do

     do i=1,mgncol
        dum = ((-nsubr(i,k)+npracs(i,k)+nnuccr(i,k)+nnuccri(i,k)-nragg(i,k))*precip_frac(i,k)- &
             nprc(i,k)*lcldm(i,k))*deltat

        if (dum.gt.nr(i,k)) then
           ratio = (nr(i,k)/deltat+nprc(i,k)*lcldm(i,k))/precip_frac(i,k)/ &
                (-nsubr(i,k)+npracs(i,k)+nnuccr(i,k)+nnuccri(i,k)-nragg(i,k))*omsm

           nragg(i,k)=nragg(i,k)*ratio
           npracs(i,k)=npracs(i,k)*ratio
           nnuccr(i,k)=nnuccr(i,k)*ratio
           nsubr(i,k)=nsubr(i,k)*ratio
           nnuccri(i,k)=nnuccri(i,k)*ratio
        end if
     end do

     if (do_cldice) then

        do i=1,mgncol

              ! conservation of qi
              !-------------------------------------------------------------------

           dum = ((-mnuccc(i,k)-mnucct(i,k)-mnudep(i,k)-msacwi(i,k))*lcldm(i,k)+(prci(i,k)+ &
                prai(i,k))*icldm(i,k)-mnuccri(i,k)*precip_frac(i,k) &
                -ice_sublim(i,k)-vap_dep(i,k)-berg(i,k)-mnuccd(i,k))*deltat

           if (dum.gt.qi(i,k)+dqidt(i,k)*deltat+D_eros_i(i,k)*deltat) then
              ratio = (qi(i,k)/deltat+max(dqidt(i,k),0.0)+vap_dep(i,k)+berg(i,k)+mnuccd(i,k)+ &
                   (mnuccc(i,k)+mnucct(i,k)+mnudep(i,k)+msacwi(i,k))*lcldm(i,k)+ &
                   mnuccri(i,k)*precip_frac(i,k))/ &
                   ((prci(i,k)+prai(i,k))*icldm(i,k)-ice_sublim(i,k)-D_eros_i(i,k)+max(-dqidt(i,k),0.0))*omsm

              if( dqidt(i,k) < 0.0 ) dqidt(i,k)=dqidt(i,k)*ratio
              if( ratio > 1) print*,"error ratio>1 in conservation of qi",i,k,ratio
              if( ratio < 0) then
                print*,"error ratio<0 in conservation of qi",i,k, lon(i), lat(i), it
                print*,"dum=", dum, qi(i,k), dqidt(i,k), D_eros_i(i,k), deltat
                print*,"conv of qi",  -mnuccc(i,k), -mnucct(i,k), -mnudep(i,k), -msacwi(i,k), lcldm(i,k)
                print*, "==", prci(i,k), prai(i,k), icldm(i,k), -mnuccri(i,k), precip_frac(i,k)
                print*, "+++", -ice_sublim(i,k), -vap_dep(i,k), -berg(i,k), -mnuccd(i,k)
               endif

              prci(i,k) = prci(i,k)*ratio
              prai(i,k) = prai(i,k)*ratio
              ice_sublim(i,k) = ice_sublim(i,k)*ratio
              D_eros_i(i,k)   = D_eros_i(i,k)*ratio
           end if
        end do

     end if

     if (do_cldice) then

        do i=1,mgncol

              ! conservation of ni
              !-------------------------------------------------------------------
           if (use_hetfrz_classnuc) then
              tmpfrz = nnuccc(i,k)
           else
              tmpfrz = 0._r8
           end if
           dum = ((-nnucct(i,k)-tmpfrz-nnudep(i,k)-nsacwi(i,k))*lcldm(i,k)+(nprci(i,k)+ &
                nprai(i,k)-nsubi(i,k))*icldm(i,k)-nnuccri(i,k)*precip_frac(i,k) - nnuccd(i,k) )*deltat

           if (dum.gt.ni(i,k)+nerosi(i,k)*icldm(i,k)*deltat) then
              ratio = (ni(i,k)/deltat+ nnuccd(i,k)+&
                   (nnucct(i,k)+tmpfrz+nnudep(i,k)+nsacwi(i,k))*lcldm(i,k)+ &
                   nnuccri(i,k)*precip_frac(i,k))/ &
                   ((nprci(i,k)+nprai(i,k)-nsubi(i,k)-nerosi(i,k))*icldm(i,k))*omsm
              nprci(i,k) = nprci(i,k)*ratio
              nprai(i,k) = nprai(i,k)*ratio
              nsubi(i,k) = nsubi(i,k)*ratio
              nerosi(i,k)= nerosi(i,k)*ratio
           end if
        end do
     end if

     do i=1,mgncol

        ! conservation of snow mixing ratio
        !-------------------------------------------------------------------
        dum = (-(prds(i,k)+pracs(i,k)+mnuccr(i,k))*precip_frac(i,k)-(prai(i,k)+prci(i,k))*icldm(i,k) &
             -(bergs(i,k)+psacws(i,k))*lcldm(i,k))*deltat

        if (dum.gt.qs(i,k).and.-prds(i,k).ge.qsmall) then
           ratio = (qs(i,k)/deltat+(prai(i,k)+prci(i,k))*icldm(i,k)+ &
                (bergs(i,k)+psacws(i,k))*lcldm(i,k)+(pracs(i,k)+mnuccr(i,k))*precip_frac(i,k))/ &
                precip_frac(i,k)/(-prds(i,k))*omsm
           prds(i,k)=prds(i,k)*ratio
        end if

     end do

     do i=1,mgncol

        ! conservation of snow number
        !-------------------------------------------------------------------
        ! calculate loss of number due to sublimation
        ! for now neglect sublimation of ns

!--> h1g, 2019-12-03
        if (allow_snow_num_sublimation)  then
          if (prds(i,k) < 0._r8 .and. qs(i,k) > qsmall ) then
            dum = prds(i,k)*deltat/qs(i,k)
            dum = max(-1._r8,dum)
            nsubs(i,k) = dum*ns(i,k)/deltat
          else
            nsubs(i,k)=0._r8
          endif
        else
        ! neglect sublimation of ns, allow_snow_num_sublimation = .false. h1g, 2019-12-03
          nsubs(i,k)=0._r8
        endif
!<-- h1g, 2019-12-03


        dum = ((-nsagg(i,k)-nsubs(i,k)-nnuccr(i,k))*precip_frac(i,k)-nprci(i,k)*icldm(i,k))*deltat

        if (dum.gt.ns(i,k)) then
           ratio = (ns(i,k)/deltat+nnuccr(i,k)* &
                precip_frac(i,k)+nprci(i,k)*icldm(i,k))/precip_frac(i,k)/ &
                (-nsubs(i,k)-nsagg(i,k))*omsm
           nsubs(i,k)=nsubs(i,k)*ratio
           nsagg(i,k)=nsagg(i,k)*ratio
        end if

     end do

     do i=1,mgncol

        ! next limit ice and snow sublimation and rain evaporation
        ! get estimate of q and t at end of time step
        ! don't include other microphysical processes since they haven't
        ! been limited via conservation checks yet

        if ((pre(i,k)+prds(i,k))*precip_frac(i,k)+ice_sublim(i,k) < -1.e-20_r8) then

           qtmp=q(i,k)-(ice_sublim(i,k)+vap_dep(i,k)+mnuccd(i,k)+ &
                (pre(i,k)+prds(i,k))*precip_frac(i,k))*deltat
           ttmp=t(i,k)+((pre(i,k)*precip_frac(i,k))*xxlv+ &
                (prds(i,k)*precip_frac(i,k)+vap_dep(i,k)+ice_sublim(i,k)+mnuccd(i,k))*xxls)*deltat/cpp

           if ( ttmp .lt.-150.0+273.15 .or. ttmp .gt.90+273.15)  &
 write(*,'(a, i4, 2f9.4, 15e12.3)') 'MG2: bad temperature@2117', k, lon(i), lat(i),   &
            ttmp, t(i,k), pre(i,k), precip_frac(i,k), prds(i,k),vap_dep(i,k), ice_sublim(i,k), mnuccd(i,k), &
            qtmp, q(i,k), p(i,k)

           ! use rhw to allow ice supersaturation
           call compute_qs(ttmp, p(i,k), qvn, q = q(i,k),   &
                             esat = esn, es_over_liq = .true.)

           ! modify ice/precip evaporation rate if q > qsat
           if (qtmp > qvn) then

              dum1=pre(i,k)*precip_frac(i,k)/((pre(i,k)+prds(i,k))*precip_frac(i,k)+ice_sublim(i,k))
              dum2=prds(i,k)*precip_frac(i,k)/((pre(i,k)+prds(i,k))*precip_frac(i,k)+ice_sublim(i,k))
              ! recalculate q and t after vap_dep and mnuccd but without evap or sublim
              qtmp=q(i,k)-(vap_dep(i,k)+mnuccd(i,k))*deltat
              ttmp=t(i,k)+((vap_dep(i,k)+mnuccd(i,k))*xxls)*deltat/cpp

              ! use rhw to allow ice supersaturation
              call compute_qs(ttmp, p(i,k), qvn, q = q(i,k),   &
                                esat = esn, es_over_liq = .true.)

              dum=(qtmp-qvn)/(1._r8 + xxlv_squared*qvn/(cpp*rv*ttmp**2))
              dum=min(dum,0._r8)

              ! modify rates if needed, divide by precip_frac to get local (in-precip) value
              pre(i,k)=dum*dum1/deltat/precip_frac(i,k)

              ! do separately using RHI for prds and ice_sublim
               call compute_qs(ttmp, p(i,k), qvn, q = q(i,k),   &
                                esat = esn, es_over_liq_and_ice = .true. )

              dum=(qtmp-qvn)/(1._r8 + xxls_squared*qvn/(cpp*rv*ttmp**2))
              dum=min(dum,0._r8)

              ! modify rates if needed, divide by precip_frac to get local (in-precip) value
              prds(i,k) = dum*dum2/deltat/precip_frac(i,k)

              ! don't divide ice_sublim by cloud fraction since it is grid-averaged
              dum1 = (1._r8-dum1-dum2)
              ice_sublim(i,k) = dum*dum1/deltat
           end if
        end if

     end do

     ! Big "administration" loop enforces conservation, updates variables
     ! that accumulate over substeps, and sets output variables.

     do i=1,mgncol

        ! get tendencies due to microphysical conversion processes
        !==========================================================
        ! note: tendencies are multiplied by appropriate cloud/precip
        ! fraction to get grid-scale values
        ! note: vap_dep is already grid-average values

        ! The net tendencies need to be added to rather than overwritten,
        ! because they may have a value already set for instantaneous
        ! melting/freezing.

        qvlat(i,k) = qvlat(i,k)-dqcdt(i,k)-dqidt(i,k)-(pre(i,k)+prds(i,k))*precip_frac(i,k)-&
             vap_dep(i,k)-ice_sublim(i,k)-mnuccd(i,k)-mnudep(i,k)*lcldm(i,k) - D_eros_l(i,k) - D_eros_i(i,k)

      !  if ( i==7 .and. k==31 ) print*,'2001 tlat= ', i,k, 'dqcdt',dqcdt(i,k), 'D_eros_l',D_eros_l(i,k), &
      !                                                 'D_eros_i', D_eros_i(i,k), &
      !                                                 'pre',pre(i,k), 'precip_frac', precip_frac(i,k),  &
      !                                                 'prds',  prds(i,k), 'vap_dep', vap_dep(i,k),  &
      !                                                 'ice_sublim', ice_sublim(i,k), 'mnuccd', mnuccd(i,k),&
      !                                                 'mnudep', mnudep(i,k), 'bergs',bergs(i,k),         &
      !                                                 'psacws',psacws(i,k), 'mnuccc', mnuccc(i,k),    &
      !                                                 'mnucct', mnucct(i,k), 'msacwi',msacwi(i,k),     &
      !                                                 'mnuccr',mnuccr(i,k), 'pracs',pracs(i,k),   &
      !                                                 'mnuccri',mnuccri(i,k),  'berg', berg(i,k)

        tlat(i,k) = tlat(i,k)+dqcdt(i,k)*xxlv+D_eros_l(i,k)*xxlv+dqidt(i,k)*xxls+D_eros_i(i,k)*xxls+ &
                   ((pre(i,k)*precip_frac(i,k))*xxlv &
               +(prds(i,k)*precip_frac(i,k)+vap_dep(i,k)+ice_sublim(i,k)+mnuccd(i,k)+mnudep(i,k)*lcldm(i,k))*xxls+ &
             ((bergs(i,k)+psacws(i,k)+mnuccc(i,k)+mnucct(i,k)+msacwi(i,k))*lcldm(i,k)+(mnuccr(i,k)+ &
             pracs(i,k)+mnuccri(i,k))*precip_frac(i,k)+berg(i,k))*xlf)

        qctend(i,k) = qctend(i,k)+dqcdt(i,k)+D_eros_l(i,k)+ &
             (-pra(i,k)-prc(i,k)-mnuccc(i,k)-mnucct(i,k)-msacwi(i,k)- &
             psacws(i,k)-bergs(i,k))*lcldm(i,k)-berg(i,k)

        if (do_cldice) then
! Note by h1g 2017-02-24, mnudep = 0.0 be default ( use_hetfrz_classnuc = .false.  )
           qitend(i,k) = qitend(i,k)+dqidt(i,k)+D_eros_i(i,k)+ &
                (mnuccc(i,k)+mnucct(i,k)+mnudep(i,k)+msacwi(i,k))*lcldm(i,k)+(-prci(i,k)- &
                prai(i,k))*icldm(i,k)+vap_dep(i,k)+berg(i,k)+ice_sublim(i,k)+ &
                mnuccd(i,k)+mnuccri(i,k)*precip_frac(i,k)
        end if

        qrtend(i,k) = qrtend(i,k)+ &
             (pra(i,k)+prc(i,k))*lcldm(i,k)+(pre(i,k)-pracs(i,k)- &
             mnuccr(i,k)-mnuccri(i,k))*precip_frac(i,k)

        qstend(i,k) = qstend(i,k)+ &
             (prai(i,k)+prci(i,k))*icldm(i,k)+(psacws(i,k)+bergs(i,k))*lcldm(i,k)+(prds(i,k)+ &
             pracs(i,k)+mnuccr(i,k))*precip_frac(i,k)

        cmeout(i,k) = cmeout(i,k) + vap_dep(i,k) + ice_sublim(i,k) + mnuccd(i,k)

        ! add output for cmei (accumulate)
        cmeitot(i,k) = cmeitot(i,k) + vap_dep(i,k) + ice_sublim(i,k) + mnuccd(i,k) + dqidt(i,k)

        ! assign variables for trop_mozart, these are grid-average
        !-------------------------------------------------------------------
        ! evaporation/sublimation is stored here as positive term

!--> h1g, 2019-12-06
        evapsnow(i,k)  = evapsnow(i,k) - prds(i,k)*precip_frac(i,k)
        nevapr(i,k)    = nevapr(i,k)   - pre(i,k)*precip_frac(i,k)
!<-- h1g, 2019-12-06

        ! change to make sure prain is positive: do not remove snow from
        ! prain used for wet deposition
        prain(i,k) = prain(i,k)+(pra(i,k)+prc(i,k))*lcldm(i,k)+(-pracs(i,k)- &
             mnuccr(i,k)-mnuccri(i,k))*precip_frac(i,k)
        prodsnow(i,k) = prodsnow(i,k)+(prai(i,k)+prci(i,k))*icldm(i,k)+(psacws(i,k)+bergs(i,k))*lcldm(i,k)+(&
             pracs(i,k)+mnuccr(i,k))*precip_frac(i,k)

        ! following are used to calculate 1st order conversion rate of cloud water
        !    to rain and snow (1/s), for later use in aerosol wet removal routine
        ! previously, wetdepa used (prain/qc) for this, and the qc in wetdepa may be smaller than the qc
        !    used to calculate pra, prc, ... in this routine
        ! qcsinksum_rate1ord = { rate of direct transfer of cloud water to rain & snow }
        !                      (no cloud ice or bergeron terms)
!--> h1g, 2019-12-06
        qcsinksum_rate1ord(i,k) = qcsinksum_rate1ord(i,k) + (pra(i,k)+prc(i,k)+psacws(i,k))*lcldm(i,k)
        ! Avoid zero/near-zero division.
      !  qcsinksum_rate1ord(i,k) = qcsinksum_rate1ord(i,k) / &
      !       max(qc(i,k),1.0e-30_r8)
!<-- h1g, 2019-12-06

        preo(i,k) = preo(i,k) + pre(i,k)*precip_frac(i,k)
        prdso(i,k) = prdso(i,k) + prds(i,k)*precip_frac(i,k)

        eroslo(i,k) = eroslo(i,k) + D_eros_l(i,k)
        erosio(i,k) = erosio(i,k) + D_eros_i(i,k)

        cmelo(i,k)  = cmelo(i,k) + dqcdt(i,k)

        ! microphysics output, note this is grid-averaged
        pratot(i,k) = pratot(i,k)+pra(i,k)*lcldm(i,k)
        prctot(i,k) = prctot(i,k)+prc(i,k)*lcldm(i,k)
        mnuccctot(i,k) = mnuccctot(i,k)+mnuccc(i,k)*lcldm(i,k)
        mnuccttot(i,k) = mnuccttot(i,k)+mnucct(i,k)*lcldm(i,k)
        mnuccdtot(i,k) = mnuccdtot(i,k)+mnuccd(i,k)

        msacwitot(i,k) = msacwitot(i,k)+msacwi(i,k)*lcldm(i,k)
        psacwstot(i,k) = psacwstot(i,k)+psacws(i,k)*lcldm(i,k)
        bergstot(i,k) = bergstot(i,k)+bergs(i,k)*lcldm(i,k)
        bergtot(i,k) = bergtot(i,k)+berg(i,k)
        prcitot(i,k) = prcitot(i,k)+prci(i,k)*icldm(i,k)
        praitot(i,k) = praitot(i,k)+prai(i,k)*icldm(i,k)

        pracstot(i,k) = pracstot(i,k)+pracs(i,k)*precip_frac(i,k)
        mnuccrtot(i,k) = mnuccrtot(i,k)+mnuccr(i,k)*precip_frac(i,k)
        mnuccritot(i,k) = mnuccritot(i,k)+mnuccri(i,k)*precip_frac(i,k)  ! h1g, 2020-02-11

        nctend(i,k) = nctend(i,k)+nerosc(i,k)*lcldm(i,k) + npccn2(i,k) + &
             (-nnuccc(i,k)-nnucct(i,k)-npsacws(i,k)+nsubc(i,k) &
             -npra(i,k)-nprc1(i,k))*lcldm(i,k)

        if (do_cldice) then
           if (use_hetfrz_classnuc) then
              tmpfrz = nnuccc(i,k)
           else
              tmpfrz = 0._r8
           end if
           nitend(i,k) = nitend(i,k)+nerosi(i,k)*icldm(i,k)+ nnuccd(i,k)+&
                (nnucct(i,k)+tmpfrz+nnudep(i,k)+nsacwi(i,k))*lcldm(i,k)+(nsubi(i,k)-nprci(i,k)- &
                nprai(i,k))*icldm(i,k)+nnuccri(i,k)*precip_frac(i,k)
        end if

        nstend(i,k) = nstend(i,k)+(nsubs(i,k)+ &
             nsagg(i,k)+nnuccr(i,k))*precip_frac(i,k)+nprci(i,k)*icldm(i,k)

        nrtend(i,k) = nrtend(i,k)+ &
             nprc(i,k)*lcldm(i,k)+(nsubr(i,k)-npracs(i,k)-nnuccr(i,k) &
             -nnuccri(i,k)+nragg(i,k))*precip_frac(i,k)

        ! make sure that ni at advanced time step does not exceed
        ! maximum (existing N + source terms*dt), which is possible if mtime < deltat
        ! note that currently mtime = deltat
        !================================================================

! ---> h1g, 2017-03-03
        IF (diag_id%qnidt_nucclim1 + diag_id%qni_nucclim1_col > 0) &
          diag_4l(i,j,k,diag_pt%qnidt_nucclim1)  = nitend(i,k)

        if (do_cldice .and. nitend(i,k).gt.0._r8.and.ni(i,k)+nitend(i,k)*deltat.gt.nimax(i,k)) then
         !  nitend(i,k)=max(0._r8,(nimax(i,k)-ni(i,k))/deltat)
        end if

        IF (diag_id%qnidt_nucclim1 + diag_id%qni_nucclim1_col > 0) &
          diag_4l(i,j,k,diag_pt%qnidt_nucclim1)  = nitend(i,k) - diag_4l(i,j,k,diag_pt%qnidt_nucclim1)
! <--- h1g, 2017-03-03

     end do

     ! End of "administration" loop

  end do micro_vert_loop ! end k loop

  !-----------------------------------------------------
  ! convert rain/snow q and N for output to history, note,
  ! output is for gridbox average

  qrout = qrout + qr
  nrout = nrout + nr * rho
  qsout = qsout + qs
  nsout = nsout + ns * rho

  ! calculate n0r and lamr from rain mass and number
  ! divide by precip fraction to get in-precip (local) values of
  ! rain mass and number, divide by rhow to get rain number in kg^-1

  do k=1,nlev

     call size_dist_param_basic(mg_rain_props, qric(:,k), nric(:,k), lamr(:,k), mgncol, n0=n0r(:,k))

     ! Calculate rercld

     ! calculate mean size of combined rain and cloud water

     call calc_rercld(lamr(:,k), n0r(:,k), lamc(:,k), pgam(:,k), qric(:,k), qcic(:,k), ncic(:,k), &
          rercld(:,k), mgncol)

  enddo

  ! Assign variables back to start-of-timestep values
  ! Some state variables are changed before the main microphysics loop
  ! to make "instantaneous" adjustments. Afterward, we must move those changes
  ! back into the tendencies.
  ! These processes:
  !  - Droplet activation (npccn, impacts nc)
  !  - Instantaneous snow melting  (minstsm/ninstsm, impacts qr/qs/nr/ns)
  !  - Instantaneous rain freezing (minstfr/ninstrf, impacts qr/qs/nr/ns)
  !================================================================================

  ! Re-apply droplet activation tendency
 ! nc = ncn
 ! nctend = nctend + npccn2

 ! ni = nin
 ! nitend = nitend + nnuccd


  ! Re-apply rain freezing and snow melting.
  dum_2D = qs
  qs = qstmp
  qstend = qstend + (dum_2D-qs)/deltat
  if (diag_id%snow_inst + diag_id%snow_inst_col > 0)    &
             diag_4l(:,j,:,diag_pt%snow_inst )  = diag_4l(:,j,:,diag_pt%snow_inst )+(dum_2D-qs)/deltat

  dum_2D = ns
  ns = nstmp
  nstend = nstend + (dum_2D-ns)/deltat
  if (diag_id%snow_num_inst + diag_id%snow_num_inst_col > 0)    &
             diag_4l(:,j,:,diag_pt%snow_num_inst )  = diag_4l(:,j,:,diag_pt%snow_num_inst)+(dum_2D-ns)/deltat

  dum_2D = qr
  qr = qrtmp
  qrtend = qrtend + (dum_2D-qr)/deltat
  if (diag_id%rain_inst + diag_id%rain_inst_col > 0)    &
             diag_4l(:,j,:,diag_pt%rain_inst )  = diag_4l(:,j,:,diag_pt%rain_inst )+(dum_2D-qr)/deltat


  dum_2D = nr
  nr = nrtmp
  nrtend = nrtend + (dum_2D-nr)/deltat
  if (diag_id%rain_num_inst + diag_id%rain_num_inst_col > 0)    &
             diag_4l(:,j,:,diag_pt%rain_num_inst )  = diag_4l(:,j,:,diag_pt%rain_num_inst )+(dum_2D-nr)/deltat

  !.............................................................................
  !================================================================================
  ! modify to include snow. in prain & evap (diagnostic here: for wet dep)
  nevapr = nevapr + evapsnow
  prain = prain + prodsnow



  do k=1,nlev
     do i=1,mgncol

        ! calculate sedimentation for cloud water and ice
        !================================================================================

        ! update in-cloud cloud mixing ratio and number concentration
        ! with microphysical tendencies to calculate sedimentation, assign to dummy vars
        ! note: these are in-cloud values***, hence we divide by cloud fraction

        dumc(i,k) = (qc(i,k)+qctend(i,k)*deltat)/lcldm(i,k)
        dumi(i,k) = (qi(i,k)+qitend(i,k)*deltat)/icldm(i,k)
        dumnc(i,k) = max((nc(i,k)+nctend(i,k)*deltat)/lcldm(i,k),0._r8)
        dumni(i,k) = max((ni(i,k)+nitend(i,k)*deltat)/icldm(i,k),0._r8)

        dumr(i,k) = (qr(i,k)+qrtend(i,k)*deltat)/precip_frac(i,k)
        dumnr(i,k) = max((nr(i,k)+nrtend(i,k)*deltat)/precip_frac(i,k),0._r8)
        dums(i,k) = (qs(i,k)+qstend(i,k)*deltat)/precip_frac(i,k)
        dumns(i,k) = max((ns(i,k)+nstend(i,k)*deltat)/precip_frac(i,k),0._r8)


        ! switch for specification of droplet and crystal number
        if (nccons) then
           dumnc(i,k)=ncnst/rho(i,k)
        end if

        ! switch for specification of cloud ice number
        if (nicons) then
           dumni(i,k)=ninst/rho(i,k)
        end if
     enddo
  enddo

  do k=1,nlev

     ! obtain new slope parameter to avoid possible singularity

     call size_dist_param_basic(mg_ice_props, dumi(:,k), dumni(:,k), &
          lami(:,k), mgncol)

     call size_dist_param_liq(mg_liq_props, dumc(:,k), dumnc(:,k), rho(:,k), &
          pgam(:,k), lamc(:,k), mgncol)

  enddo

  do k=1,nlev
     do i=1,mgncol

        ! calculate number and mass weighted fall velocity for droplets and cloud ice
        !-------------------------------------------------------------------


        if (dumc(i,k).ge.qsmall) then

           vtrmc(i,k)=acn(i,k)*gamma(4._r8+bc+pgam(i,k))/ &
                (lamc(i,k)**bc*gamma(pgam(i,k)+4._r8))

           fc(i,k) = g*rho(i,k)*vtrmc(i,k)

           fnc(i,k) = g*rho(i,k)* &
                acn(i,k)*gamma(1._r8+bc+pgam(i,k))/ &
                (lamc(i,k)**bc*gamma(pgam(i,k)+1._r8))
        else
           fc(i,k) = 0._r8
           fnc(i,k)= 0._r8
        end if

        ! calculate number and mass weighted fall velocity for cloud ice

        if (dumi(i,k).ge.qsmall) then

           vtrmi(i,k)=min(ain(i,k)*gamma_bi_plus4/(6._r8*lami(i,k)**bi), &
                1.2_r8*rhof(i,k))

           fi(i,k) = g*rho(i,k)*vtrmi(i,k)
           fni(i,k) = g*rho(i,k)* &
                min(ain(i,k)*gamma_bi_plus1/lami(i,k)**bi,1.2_r8*rhof(i,k))

           ! adjust the ice fall velocity for smaller (r < 20 um) ice
           ! particles (blend over 18-20 um)
           irad = 1.5_r8 / lami(i,k) * 1e6_r8
           ifrac = min(1._r8, max(0._r8, (irad - 18._r8) / 2._r8))

           if (ifrac .lt. 1._r8) then
              vtrmi(i,k) = ifrac * vtrmi(i,k) + &
                 (1._r8 - ifrac) * &
                 min(ajn(i,k)*gamma_bj_plus4/(6._r8*lami(i,k)**bj), &
                 1.2_r8*rhof(i,k))

              fi(i,k) = g*rho(i,k)*vtrmi(i,k)
              fni(i,k) = ifrac * fni(i,k) + &
                 (1._r8 - ifrac) * &
                 g*rho(i,k)* &
                 min(ajn(i,k)*gamma_bj_plus1/lami(i,k)**bj,1.2_r8*rhof(i,k))
           end if
        else
           fi(i,k) = 0._r8
           fni(i,k)= 0._r8
        end if

     enddo

  enddo


  do k=1,nlev

        ! fallspeed for rain

        call size_dist_param_basic(mg_rain_props, dumr(:,k), dumnr(:,k), &
             lamr(:,k), mgncol)
  enddo

  do k=1,nlev

     do i=1,mgncol
        if (lamr(i,k).ge.qsmall) then

           ! 'final' values of number and mass weighted mean fallspeed for rain (m/s)

           unr(i,k) = min(arn(i,k)*gamma_br_plus1/lamr(i,k)**br,9.1_r8*rhof(i,k))
           umr(i,k) = min(arn(i,k)*gamma_br_plus4/(6._r8*lamr(i,k)**br),9.1_r8*rhof(i,k))

           fr(i,k) = g*rho(i,k)*umr(i,k)
           fnr(i,k) = g*rho(i,k)*unr(i,k)

        else
           fr(i,k)=0._r8
           fnr(i,k)=0._r8
        end if

        ! fallspeed for snow

        call size_dist_param_basic(mg_snow_props, dums(i,k), dumns(i,k), &
             lams(i,k))

        if (lams(i,k).ge.qsmall) then

           ! 'final' values of number and mass weighted mean fallspeed for snow (m/s)
           ums(i,k) = min(asn(i,k)*gamma_bs_plus4/(6._r8*lams(i,k)**bs),1.2_r8*rhof(i,k))
           uns(i,k) = min(asn(i,k)*gamma_bs_plus1/lams(i,k)**bs,1.2_r8*rhof(i,k))

           fs(i,k) = g*rho(i,k)*ums(i,k)
           fns(i,k) = g*rho(i,k)*uns(i,k)

        else
           fs(i,k)=0._r8
           fns(i,k)=0._r8
        end if

        ! redefine dummy variables - sedimentation is calculated over grid-scale
        ! quantities to ensure conservation

        dumc(i,k) = (qc(i,k)+qctend(i,k)*deltat)
        dumnc(i,k) = max((nc(i,k)+nctend(i,k)*deltat),0._r8)
        dumi(i,k) = (qi(i,k)+qitend(i,k)*deltat)
        dumni(i,k) = max((ni(i,k)+nitend(i,k)*deltat),0._r8)
        dumr(i,k) = (qr(i,k)+qrtend(i,k)*deltat)
        dumnr(i,k) = max((nr(i,k)+nrtend(i,k)*deltat),0._r8)
        dums(i,k) = (qs(i,k)+qstend(i,k)*deltat)
        dumns(i,k) = max((ns(i,k)+nstend(i,k)*deltat),0._r8)

        if (dumc(i,k).lt.qsmall) dumnc(i,k)=0._r8
        if (dumi(i,k).lt.qsmall) dumni(i,k)=0._r8
        if (dumr(i,k).lt.qsmall) dumnr(i,k)=0._r8
        if (dums(i,k).lt.qsmall) dumns(i,k)=0._r8

     enddo
  end do       !!! vertical loop

if ( do_implicit_fall ) then
  fc  = vfac_drop * fc/g/rho
  fnc = vfac_drop * fnc/g/rho
  fi  = vfac_ice  * fi/g/rho
  fni = vfac_ice  * fni/g/rho

  fr  = vfactor * fr/g/rho
  fnr = vfactor * fnr/g/rho
  fs  = vfactor * fs/g/rho
  fns = vfactor * fns/g/rho

 ! cloud water (mass) sedimentation
  do i=1,mgncol
    dum_1D(:) = dumc(i,:)
    call implicit_fall ( deltat, 1, nlev, zhalf(i,:) , fc(i,:), pdel(i,:), dum_1D, precip, flx)
    do k=1,nlev
      lflx(i,k+1) = lflx(i,k+1) + flx(k)/g/deltat
      qcsedten(i,k)= qcsedten(i,k) + (dum_1D(k) - dumc(i,k))/deltat
      qctend(i,k)  = qctend(i,k)   + (dum_1D(k) - dumc(i,k))/deltat
    enddo
    if ( precip .ge. 0.0_r8 ) then  !h1g, 2019-11-26, ensure numerical stability
      prect(i) = prect(i)+precip/g/deltat/1000._r8
    else 
      qvlat(i,nlev) = qvlat(i,nlev) + precip/deltat/pdel(i,nlev)
      tlat(i,nlev)  = tlat(i,nlev) - precip/deltat/pdel(i,nlev) * xxlv
    endif
  enddo

 ! cloud water (number) sedimentation
  do i=1,mgncol
    dum_1D(:) = dumnc(i,:)
    call implicit_fall ( deltat, 1, nlev, zhalf(i,:) , fnc(i,:), pdel(i,:), dum_1D, precip, flx)
    do k=1,nlev
      nctend(i,k)  = nctend(i,k)   + (dum_1D(k) - dumnc(i,k))/deltat
      IF ( diag_id%qndt_sedi + diag_id%qn_sedi_col > 0 ) &
           diag_4l(i,j,k,diag_pt%qndt_sedi) =     &
           diag_4l(i,j,k,diag_pt%qndt_sedi) + (dum_1D(k) - dumnc(i,k))/deltat
    enddo
  enddo

 ! cloud ice (mass) sedimentation
  do i=1,mgncol
    dum_1D(:) = dumi(i,:)
    call implicit_fall ( deltat, 1, nlev, zhalf(i,:) , fi(i,:), pdel(i,:), dum_1D, precip, flx)
    do k=1,nlev
      iflx(i,k+1) = iflx(i,k+1) + flx(k)/g/deltat
      qisedten(i,k)= qisedten(i,k) + (dum_1D(k) - dumi(i,k))/deltat
      qitend(i,k)  = qitend(i,k)   + (dum_1D(k) - dumi(i,k))/deltat
    enddo
    if ( precip .ge. 0.0_r8 ) then !h1g, 2019-11-26, ensure numerical stability
      prect(i) = prect(i) + precip/g/deltat/1000._r8
      preci(i) = preci(i) + precip/g/deltat/1000._r8
    else 
      qvlat(i,nlev) = qvlat(i,nlev) + precip/deltat/pdel(i,nlev)
      tlat(i,nlev)  = tlat(i,nlev) - precip/deltat/pdel(i,nlev) * xxls
    endif
  enddo

 ! cloud ice (number) sedimentation
  do i=1,mgncol
    dum_1D(:) = dumni(i,:)
    call implicit_fall ( deltat, 1, nlev, zhalf(i,:) , fni(i,:), pdel(i,:), dum_1D, precip, flx)
    do k=1,nlev
      nitend(i,k)  = nitend(i,k)   + (dum_1D(k) - dumni(i,k))/deltat
      IF ( diag_id%qnidt_sedi +  diag_id%qni_sedi_col > 0 ) &
           diag_4l(i,j,k,diag_pt%qnidt_sedi) =    &
           diag_4l(i,j,k,diag_pt%qnidt_sedi) + (dum_1D(k) - dumni(i,k))/deltat
    enddo
  enddo

 ! rain water (mass) sedimentation
  do i=1,mgncol
     dum_1D(:) = dumr(i,:)
     call implicit_fall ( deltat, 1, nlev, zhalf(i,:) , fr(i,:), pdel(i,:), dum_1D, precip, flx)
     do k=1,nlev
       if ( flx(k) .ge. qsmall ) rflx(i,k+1) = rflx(i,k+1) + flx(k)/g/deltat !h1g, 2019-11-26, ensure numerical stability
       qrsedten(i,k)= qrsedten(i,k) + (dum_1D(k) - dumr(i,k))/deltat
       qrtend (i,k) = qrtend(i,k)   + (dum_1D(k) - dumr(i,k))/deltat
     enddo
     if ( precip .ge. 0.0_r8 ) then !h1g, 2019-11-26, ensure numerical stability
       prect(i) = prect(i)+precip/g/deltat/1000._r8
     else 
      qvlat(i,nlev) = qvlat(i,nlev) + precip/deltat/pdel(i,nlev)
      tlat(i,nlev)  = tlat(i,nlev) - precip/deltat/pdel(i,nlev) * xxlv
     endif
  enddo

 ! rain water (number) sedimentation
  do i=1,mgncol
    dum_1D(:) = dumnr(i,:)
    call implicit_fall ( deltat, 1, nlev, zhalf(i,:) , fnr(i,:), pdel(i,:), dum_1D, precip, flx)
    do k=1,nlev
      nrtend(i,k)  = nrtend(i,k)   + (dum_1D(k) - dumnr(i,k))/deltat
      IF ( diag_id%rain_num_sedi +  diag_id%rain_num_sedi_col > 0 ) &
           diag_4l(i,j,k,diag_pt%rain_num_sedi) =    &
           diag_4l(i,j,k,diag_pt%rain_num_sedi) + (dum_1D(k) - dumnr(i,k))/deltat
    enddo
  enddo


 ! snow water (mass) sedimentation
  do i=1,mgncol
    dum_1D(:) = dums(i,:)
    call implicit_fall ( deltat, 1, nlev, zhalf(i,:) , fs(i,:), pdel(i,:), dum_1D, precip, flx)
    do k=1,nlev
       if ( flx(k) .ge. qsmall ) sflx(i,k+1) = sflx(i,k+1) + flx(k)/g/deltat !h1g, 2019-11-26, ensure numerical stability
       qssedten(i,k)= qssedten(i,k) + (dum_1D(k) - dums(i,k))/deltat
       qstend(i,k)  = qstend(i,k)   + (dum_1D(k) - dums(i,k))/deltat
    enddo
    if ( precip .ge. 0.0_r8 ) then !h1g, 2019-11-26, ensure numerical stability
      prect(i) = prect(i)+precip/g/deltat/1000._r8
      preci(i) = preci(i)+precip/g/deltat/1000._r8
    else 
      qvlat(i,nlev) = qvlat(i,nlev) + precip/deltat/pdel(i,nlev)
      tlat(i,nlev)  = tlat(i,nlev) - precip/deltat/pdel(i,nlev) * xxls
    endif
  enddo
 ! snow water (number) sedimentation
  do i=1,mgncol
     dum_1D(:) = dumns(i,:)
     call implicit_fall ( deltat, 1, nlev, zhalf(i,:) , fns(i,:), pdel(i,:), dum_1D, precip, flx)
     do k=1,nlev
       nstend(i,k)  = nstend(i,k)  + (dum_1D(k) - dumns(i,k))/deltat
        IF ( diag_id%snow_num_sedi +  diag_id%snow_num_sedi_col > 0 ) &
             diag_4l(i,j,k,diag_pt%snow_num_sedi) =    &
             diag_4l(i,j,k,diag_pt%snow_num_sedi) + (dum_1D(k) - dumns(i,k))/deltat
     enddo
  enddo

else

  do k=1,nlev
     do i=1,mgncol
       pdel_inv(i,k) = 1._r8/pdel(i,k)
     enddo
  enddo

  ! initialize nstep for sedimentation sub-steps

  ! calculate number of split time steps to ensure courant stability criteria
  ! for sedimentation calculations
  !-------------------------------------------------------------------

  do i=1,mgncol

     nstep = 1 + int(max( &
          maxval( fi(i,:)*pdel_inv(i,:)), &
          maxval(fni(i,:)*pdel_inv(i,:))) &
          * deltat)

!     if ( mpp_pe() == mpp_root_pe() ) &
!     write(*, *)  "nstep = ", nstep


     ! loop over sedimentation sub-time step to ensure stability
     !==============================================================
     do n = 1,nstep

        if (do_cldice) then
           falouti  = fi(i,:)  * dumi(i,:)
           faloutni = fni(i,:) * dumni(i,:)
        else
           falouti  = 0._r8
           faloutni = 0._r8
        end if

        ! top of model
        k = 1

        ! add fallout terms to microphysical tendencies
        faltndi = falouti(k)/pdel(i,k)
        faltndni = faloutni(k)/pdel(i,k)
        qitend(i,k) = qitend(i,k)-faltndi/nstep
        nitend(i,k) = nitend(i,k)-faltndni/nstep

        ! sedimentation tendency for output
        qisedten(i,k)=qisedten(i,k)-faltndi/nstep
        IF ( diag_id%qnidt_sedi +  diag_id%qni_sedi_col > 0 ) &
             diag_4l(i,j,k,diag_pt%qnidt_sedi) =    &
                   diag_4l(i,j,k,diag_pt%qnidt_sedi) - faltndni/nstep

        dumi(i,k) = dumi(i,k)-faltndi*deltat/nstep
        dumni(i,k) = dumni(i,k)-faltndni*deltat/nstep

        do k = 2,nlev

           ! for cloud liquid and ice, if cloud fraction increases with height
           ! then add flux from above to both vapor and cloud water of current level
           ! this means that flux entering clear portion of cell from above evaporates
           ! instantly

           ! note: this is not an issue with precip, since we assume max overlap
           dum1=icldm(i,k)/icldm(i,k-1)
           dum1=min(dum1,1._r8)

!--> h1g, 2019-12-18
           if ( no_evap_in_sedimentation ) dum1 = 1.0
!<-- h1g, 2019-12-18

           faltndqie=(falouti(k)-falouti(k-1))/pdel(i,k)

           faltndi=(falouti(k)-dum1*falouti(k-1))/pdel(i,k)
           faltndni=(faloutni(k)-dum1*faloutni(k-1))/pdel(i,k)

           ! add fallout terms to eulerian tendencies

           qitend(i,k) = qitend(i,k)-faltndi/nstep
           nitend(i,k) = nitend(i,k)-faltndni/nstep

           ! sedimentation tendency for output
           qisedten(i,k)=qisedten(i,k)-faltndi/nstep
           IF ( diag_id%qnidt_sedi +  diag_id%qni_sedi_col > 0 ) &
                diag_4l(i,j,k,diag_pt%qnidt_sedi) =    &
                      diag_4l(i,j,k,diag_pt%qnidt_sedi) - faltndni/nstep

           ! add terms to to evap/sub of ice water

           qvlat(i,k)=qvlat(i,k)-(faltndqie-faltndi)/nstep
           ! for output
           qisevap(i,k)=qisevap(i,k)-(faltndqie-faltndi)/nstep
           tlat(i,k)=tlat(i,k)+(faltndqie-faltndi)*xxls/nstep

           dumi(i,k) = dumi(i,k)-faltndi*deltat/nstep
           dumni(i,k) = dumni(i,k)-faltndni*deltat/nstep

        end do

        ! Ice flux
        do k = 1,nlev
          iflx(i,k+1) = iflx(i,k+1) + falouti(k) / g / real(nstep)
        end do

        ! units below are m/s
        ! sedimentation flux at surface is added to precip flux at surface
        ! to get total precip (cloud + precip water) rate

        prect(i) = prect(i)+falouti(nlev)/g/real(nstep)/1000._r8
        preci(i) = preci(i)+falouti(nlev)/g/real(nstep)/1000._r8
     end do


     ! calculate number of split time steps to ensure courant stability criteria
     ! for sedimentation calculations
     !-------------------------------------------------------------------
     nstep = 1 + int(max( &
          maxval( fc(i,:)*pdel_inv(i,:)), &
          maxval(fnc(i,:)*pdel_inv(i,:))) &
          * deltat)

     ! loop over sedimentation sub-time step to ensure stability
     !==============================================================
     do n = 1,nstep

        faloutc  = fc(i,:)  * dumc(i,:)
        faloutnc = fnc(i,:) * dumnc(i,:)

        ! top of model
        k = 1

        ! add fallout terms to microphysical tendencies
        faltndc = faloutc(k)/pdel(i,k)
        faltndnc = faloutnc(k)/pdel(i,k)
        qctend(i,k) = qctend(i,k)-faltndc/nstep
        nctend(i,k) = nctend(i,k)-faltndnc/nstep

        ! sedimentation tendency for output
        qcsedten(i,k)=qcsedten(i,k)-faltndc/nstep

        IF ( diag_id%qndt_sedi + diag_id%qn_sedi_col > 0 ) &
               diag_4l(i,j,k,diag_pt%qndt_sedi) =     &
               diag_4l(i,j,k,diag_pt%qndt_sedi) - faltndnc/nstep


        dumc(i,k) = dumc(i,k)-faltndc*deltat/nstep
        dumnc(i,k) = dumnc(i,k)-faltndnc*deltat/nstep

        do k = 2,nlev

           dum=lcldm(i,k)/lcldm(i,k-1)
           dum=min(dum,1._r8)

!--> h1g, 2019-12-18
           if ( no_evap_in_sedimentation ) dum = 1.0
!<-- h1g, 2019-12-18

           faltndqce=(faloutc(k)-faloutc(k-1))/pdel(i,k)
           faltndc=(faloutc(k)-dum*faloutc(k-1))/pdel(i,k)
           faltndnc=(faloutnc(k)-dum*faloutnc(k-1))/pdel(i,k)

           ! add fallout terms to eulerian tendencies
           qctend(i,k) = qctend(i,k)-faltndc/nstep
           nctend(i,k) = nctend(i,k)-faltndnc/nstep

           ! sedimentation tendency for output
           qcsedten(i,k)=qcsedten(i,k)-faltndc/nstep
           IF ( diag_id%qndt_sedi + diag_id%qn_sedi_col > 0 ) &
                diag_4l(i,j,k,diag_pt%qndt_sedi) =     &
                diag_4l(i,j,k,diag_pt%qndt_sedi) - faltndnc/nstep


           ! add terms to to evap/sub of cloud water
           qvlat(i,k)=qvlat(i,k)-(faltndqce-faltndc)/nstep
           ! for output
           qcsevap(i,k)=qcsevap(i,k)-(faltndqce-faltndc)/nstep
           tlat(i,k)=tlat(i,k)+(faltndqce-faltndc)*xxlv/nstep

           dumc(i,k) = dumc(i,k)-faltndc*deltat/nstep
           dumnc(i,k) = dumnc(i,k)-faltndnc*deltat/nstep

        end do

        !Liquid condensate flux here
        do k = 1,nlev
           lflx(i,k+1) = lflx(i,k+1) + faloutc(k) / g / real(nstep)
        end do

        prect(i) = prect(i)+faloutc(nlev)/g/real(nstep)/1000._r8
     end do

     ! calculate number of split time steps to ensure courant stability criteria
     ! for sedimentation calculations
     !-------------------------------------------------------------------
     nstep = 1 + int(max( &
          maxval( fr(i,:)*pdel_inv(i,:)), &
          maxval(fnr(i,:)*pdel_inv(i,:))) &
          * deltat)


    ! loop over sedimentation sub-time step to ensure stability
     !==============================================================
     do n = 1,nstep

        faloutr  = fr(i,:)  * dumr(i,:)
        faloutnr = fnr(i,:) * dumnr(i,:)

        ! top of model
        k = 1
        ! add fallout terms to microphysical tendencies
        faltndr = faloutr(k)/pdel(i,k)
        faltndnr = faloutnr(k)/pdel(i,k)
        qrtend(i,k) = qrtend(i,k)-faltndr/nstep
        nrtend(i,k) = nrtend(i,k)-faltndnr/nstep

        ! sedimentation tendency for output
        qrsedten(i,k)=qrsedten(i,k)-faltndr/nstep
        IF ( diag_id%rain_num_sedi +  diag_id%rain_num_sedi_col > 0 ) &
             diag_4l(i,j,k,diag_pt%rain_num_sedi) =    &
                   diag_4l(i,j,k,diag_pt%rain_num_sedi) - faltndnr/nstep

        dumr(i,k) = dumr(i,k)-faltndr*deltat/real(nstep)
        dumnr(i,k) = dumnr(i,k)-faltndnr*deltat/real(nstep)

        do k = 2,nlev

           faltndr=(faloutr(k)-faloutr(k-1))/pdel(i,k)
           faltndnr=(faloutnr(k)-faloutnr(k-1))/pdel(i,k)

           ! add fallout terms to eulerian tendencies
           qrtend(i,k) = qrtend(i,k)-faltndr/nstep
           nrtend(i,k) = nrtend(i,k)-faltndnr/nstep

           ! sedimentation tendency for output
           qrsedten(i,k)=qrsedten(i,k)-faltndr/nstep
        IF ( diag_id%rain_num_sedi +  diag_id%rain_num_sedi_col > 0 ) &
             diag_4l(i,j,k,diag_pt%rain_num_sedi) =    &
                   diag_4l(i,j,k,diag_pt%rain_num_sedi) - faltndnr/nstep

           dumr(i,k) = dumr(i,k)-faltndr*deltat/real(nstep)
           dumnr(i,k) = dumnr(i,k)-faltndnr*deltat/real(nstep)
        end do

        ! Rain Flux
        do k = 1,nlev
           rflx(i,k+1) = rflx(i,k+1) + faloutr(k) / g / real(nstep)
        end do
        prect(i) = prect(i)+faloutr(nlev)/g/real(nstep)/1000._r8

     end do

     ! calculate number of split time steps to ensure courant stability criteria
     ! for sedimentation calculations
     !-------------------------------------------------------------------
     nstep = 1 + int(max( &
          maxval( fs(i,:)*pdel_inv(i,:)), &
          maxval(fns(i,:)*pdel_inv(i,:))) &
          * deltat)

     ! loop over sedimentation sub-time step to ensure stability
     !==============================================================
     do n = 1,nstep

        falouts  = fs(i,:)  * dums(i,:)
        faloutns = fns(i,:) * dumns(i,:)

        ! top of model
        k = 1

        ! add fallout terms to microphysical tendencies
        faltnds = falouts(k)/pdel(i,k)
        faltndns = faloutns(k)/pdel(i,k)
        qstend(i,k) = qstend(i,k)-faltnds/nstep
        nstend(i,k) = nstend(i,k)-faltndns/nstep

        ! sedimentation tendency for output
        qssedten(i,k)=qssedten(i,k)-faltnds/nstep
        IF ( diag_id%snow_num_sedi +  diag_id%snow_num_sedi_col > 0 ) &
             diag_4l(i,j,k,diag_pt%snow_num_sedi) =    &
                   diag_4l(i,j,k,diag_pt%snow_num_sedi) - faltndns/nstep

        dums(i,k) = dums(i,k)-faltnds*deltat/real(nstep)
        dumns(i,k) = dumns(i,k)-faltndns*deltat/real(nstep)

        do k = 2,nlev

           faltnds=(falouts(k)-falouts(k-1))/pdel(i,k)
           faltndns=(faloutns(k)-faloutns(k-1))/pdel(i,k)

           ! add fallout terms to eulerian tendencies
           qstend(i,k) = qstend(i,k)-faltnds/nstep
           nstend(i,k) = nstend(i,k)-faltndns/nstep

           ! sedimentation tendency for output
           qssedten(i,k)=qssedten(i,k)-faltnds/nstep

           dums(i,k) = dums(i,k)-faltnds*deltat/real(nstep)
           dumns(i,k) = dumns(i,k)-faltndns*deltat/real(nstep)

        end do   !! k loop

        ! Snow Flux
        do k = 1,nlev
           sflx(i,k+1) = sflx(i,k+1) + falouts(k) / g / real(nstep)
        end do
        prect(i) = prect(i)+falouts(nlev)/g/real(nstep)/1000._r8
        preci(i) = preci(i)+falouts(nlev)/g/real(nstep)/1000._r8

     end do   !! nstep loop

  enddo !! i loop

endif
! end sedimentation

  tlat1 = tlat1 + tlat
  t = t + tlat*deltat/cpp

  qvlat1 = qvlat1 + qvlat
  q = q + qvlat*deltat

  qctend1 = qctend1 + qctend
  qc = qc + qctend*deltat

  qitend1 = qitend1 + qitend
  qi = qi + qitend*deltat

  nctend1 = nctend1 + nctend
  nc = nc + nctend*deltat

  nitend1 = nitend1 + nitend
  ni = ni + nitend*deltat

  qrtend1 = qrtend1 + qrtend
  qr = qr + qrtend*deltat

  qstend1 = qstend1 + qstend
  qs = qs + qstend*deltat

  nrtend1 = nrtend1 + nrtend
  nr = nr + nrtend*deltat

  nstend1 = nstend1 + nstend
  ns = ns + nstend*deltat

  prect1 = prect1 + prect
  preci1 = preci1 + preci

!--> h1g, 2019-12-12, remove tiny or negative hydrometeor mass or number
!--> in order to avoid numerical instability
  do k=1, nlev
    do i=1, mgncol
      if ( qi(i,k).lt. qsmall ) then
        if (diag_id%qidt_tiny + diag_id%qi_tiny_col > 0) &
            diag_4l(i,j,k,diag_pt%qidt_tiny) =           &
            diag_4l(i,j,k,diag_pt%qidt_tiny) - qi(i,k)/deltat
        qitend1(i,k)     = qitend1(i,k) - qi(i,k)/deltat
        qvlat1(i,k)      = qvlat1(i,k) + qi(i,k)/deltat
        q(i,k)           = q(i,k)      + qi(i,k)
        tlat1(i,k)       = tlat1(i,k)  - qi(i,k)/deltat*xxls
        t(i,k)           = t(i,k)      - qi(i,k)*xxls/cpp
        qi(i,k)          = 0.0
      endif
    enddo
  enddo

  do k=1, nlev
    do i=1, mgncol
      if ( ni(i,k).lt. qsmall ) then
        if (diag_id%qnidt_tiny + diag_id%qni_tiny_col > 0) &
            diag_4l(i,j,k,diag_pt%qnidt_tiny) =           &
            diag_4l(i,j,k,diag_pt%qnidt_tiny) - ni(i,k)/deltat
        nitend1(i,k)     = nitend1(i,k) - ni(i,k)/deltat
        ni(i,k)          = 0.0
      endif
    enddo
  enddo

  do k=1, nlev
    do i=1, mgncol
      if ( qs(i,k).lt. qsmall ) then
        if (diag_id%qsdt_tiny + diag_id%qs_tiny_col > 0) &
            diag_4l(i,j,k,diag_pt%qsdt_tiny) =           &
            diag_4l(i,j,k,diag_pt%qsdt_tiny) - qs(i,k)/deltat
        qstend1(i,k)     = qstend1(i,k) - qs(i,k)/deltat
        qvlat1(i,k)      = qvlat1(i,k) + qs(i,k)/deltat
        q(i,k)           = q(i,k)      + qs(i,k)
        tlat1(i,k)       = tlat1(i,k)  - qs(i,k)/deltat*xxls
        t(i,k)           = t(i,k)      - qs(i,k)*xxls/cpp
        qs(i,k)          = 0.0
      endif
    enddo
  enddo

  do k=1, nlev
    do i=1, mgncol
      if ( ns(i,k).lt. qsmall ) then
        if (diag_id%qnsdt_tiny + diag_id%qns_tiny_col > 0) &
            diag_4l(i,j,k,diag_pt%qnsdt_tiny) =           &
            diag_4l(i,j,k,diag_pt%qnsdt_tiny) - ns(i,k)/deltat
        nstend1(i,k)     = nstend1(i,k) - ns(i,k)/deltat
        ns(i,k)          = 0.0
      endif
    enddo
  enddo

  do k=1, nlev
    do i=1, mgncol
      if ( qc(i,k).lt. qsmall ) then
        if (diag_id%qldt_tiny + diag_id%ql_tiny_col > 0) &
            diag_4l(i,j,k,diag_pt%qldt_tiny) =           &
            diag_4l(i,j,k,diag_pt%qldt_tiny) - qc(i,k)/deltat
        qctend1(i,k)     = qctend1(i,k) - qc(i,k)/deltat
        qvlat1(i,k)      = qvlat1(i,k) + qc(i,k)/deltat
        q(i,k)           = q(i,k)      + qc(i,k)
        tlat1(i,k)       = tlat1(i,k)  - qc(i,k)/deltat*xxlv
        t(i,k)           = t(i,k)      - qc(i,k)*xxlv/cpp
        qc(i,k)          = 0.0
      endif
    enddo
  enddo

  do k=1, nlev
    do i=1, mgncol
      if ( nc(i,k).lt.qsmall ) then
        if (diag_id%qndt_tiny + diag_id%qn_tiny_col > 0) &
            diag_4l(i,j,k,diag_pt%qndt_tiny) =           &
            diag_4l(i,j,k,diag_pt%qndt_tiny) - nc(i,k)/deltat
        nctend1(i,k)     = nctend1(i,k) - nc(i,k)/deltat
        nc(i,k)          = 0.0
      endif
    enddo
  enddo

  do k=1, nlev
    do i=1, mgncol
      if ( qr(i,k).lt. qsmall ) then
        if (diag_id%qrdt_tiny + diag_id%qr_tiny_col > 0) &
            diag_4l(i,j,k,diag_pt%qrdt_tiny) =           &
            diag_4l(i,j,k,diag_pt%qrdt_tiny) - qr(i,k)/deltat
        qrtend1(i,k)     = qrtend1(i,k) - qr(i,k)/deltat
        qvlat1(i,k)      = qvlat1(i,k) + qr(i,k)/deltat
        q(i,k)           = q(i,k)      + qr(i,k)
        tlat1(i,k)       = tlat1(i,k)  - qr(i,k)/deltat*xxlv
        t(i,k)           = t(i,k)      - qr(i,k)*xxlv/cpp
        qr(i,k)          = 0.0
      endif
    enddo
  enddo

  do k=1, nlev
    do i=1, mgncol
      if ( nr(i,k).lt. qsmall ) then
        if (diag_id%qnrdt_tiny + diag_id%qnr_tiny_col > 0) &
            diag_4l(i,j,k,diag_pt%qnrdt_tiny) =           &
            diag_4l(i,j,k,diag_pt%qnrdt_tiny) - nr(i,k)/deltat
        nrtend1(i,k)     = nrtend1(i,k) - nr(i,k)/deltat
        nr(i,k)          = 0.0
      endif
    enddo
  enddo

!<-- h1g, 2019-12-12


  npccno     = npccno   + npccn2
  nprao      = nprao    - npra*lcldm
  nprc1o     = nprc1o   - nprc1*lcldm
  nerosco    = nerosco  + nerosc*lcldm

  nnuccco    = nnuccco - nnuccc*lcldm
  nnuccto    = nnuccto  - nnucct*lcldm
  npsacwso   = npsacwso - npsacws*lcldm
  nsubco     = nsubco   + nsubc*lcldm
  nucclimo   = nucclimo + nucclim

  nnuccdo    = nnuccdo + nnuccd
  nerosio    = nerosio + nerosi*icldm
  nsacwio    = nsacwio + nsacwi*lcldm
  nsubio     = nsubio  + nsubi*icldm
  nprcio     = nprcio  - nprci*icldm
  npraio     = npraio  - nprai*icldm
  nnuccrio   = nnuccrio+ nnuccri*precip_frac
  nucclim1io = nucclim1io + nucclim1i

end do substepping ! iter loop, sub-step
  deltat = deltatin

prect  = prect1/real(iter)
preci  = preci1/real(iter)

lflx   = lflx /real(iter)
iflx   = iflx /real(iter)
rflx   = rflx /real(iter)
sflx   = sflx /real(iter)

if ( include_ice_in_snowflx ) sflx   = sflx + iflx  ! h1g, 2024-01-31

qcsedten = qcsedten/real(iter)
qisedten = qisedten/real(iter)
qrsedten = qrsedten/real(iter)
qssedten = qssedten/real(iter)

diag_4l(:,j,:,diag_pt%qndt_sedi)     = diag_4l(:,j,:,diag_pt%qndt_sedi)  /real(iter)
diag_4l(:,j,:,diag_pt%qnidt_sedi)    = diag_4l(:,j,:,diag_pt%qnidt_sedi) /real(iter)
diag_4l(:,j,:,diag_pt%rain_num_sedi) = diag_4l(:,j,:,diag_pt%rain_num_sedi)/real(iter)
diag_4l(:,j,:,diag_pt%snow_num_sedi) = diag_4l(:,j,:,diag_pt%snow_num_sedi)/real(iter)

  ! assign variables back to start-of-timestep values before updating after sub-steps
  !================================================================================

  t  = tn
  q  = qn
  qc = qcn
  nc = ncn
  qi = qin
  ni = nin
  qr = qrn
  nr = nrn
  qs = qsn
  ns = nsn

  tlat = tlat1/real(iter)
  qvlat = qvlat1/real(iter)
  qctend = qctend1/real(iter)
  qitend = qitend1/real(iter)
  nctend = nctend1/real(iter)
  nitend = nitend1/real(iter)

  qrtend = qrtend1/real(iter)
  qstend = qstend1/real(iter)
  nrtend = nrtend1/real(iter)
  nstend = nstend1/real(iter)

  diag_4l(:,j,:,diag_pt%qidt_tiny)  = diag_4l(:,j,:,diag_pt%qidt_tiny)/real(iter)
  diag_4l(:,j,:,diag_pt%qnidt_tiny) = diag_4l(:,j,:,diag_pt%qnidt_tiny)/real(iter)

  diag_4l(:,j,:,diag_pt%qsdt_tiny)  = diag_4l(:,j,:,diag_pt%qsdt_tiny)/real(iter)
  diag_4l(:,j,:,diag_pt%qnsdt_tiny) = diag_4l(:,j,:,diag_pt%qnsdt_tiny)/real(iter)

  diag_4l(:,j,:,diag_pt%qldt_tiny)  = diag_4l(:,j,:,diag_pt%qldt_tiny)/real(iter)
  diag_4l(:,j,:,diag_pt%qndt_tiny)  = diag_4l(:,j,:,diag_pt%qndt_tiny)/real(iter)

  diag_4l(:,j,:,diag_pt%qrdt_tiny)  = diag_4l(:,j,:,diag_pt%qrdt_tiny)/real(iter)
  diag_4l(:,j,:,diag_pt%qnrdt_tiny) = diag_4l(:,j,:,diag_pt%qnrdt_tiny)/real(iter)



  ! divide output precip q and N by number of sub-steps to get average over time step
  !================================================================================

  qrout = qrout/real(iter)
  qsout = qsout/real(iter)
  nrout = nrout/real(iter)
  nsout = nsout/real(iter)

  ! divide trop_mozart variables by number of sub-steps to get average over time step
  !================================================================================

  nevapr = nevapr/real(iter)
  evapsnow = evapsnow/real(iter)
  prain = prain/real(iter)
  prodsnow = prodsnow/real(iter)

  cmeout = cmeout/real(iter)

  cmeitot = cmeitot/real(iter)
  meltsdttot = meltsdttot/real(iter)
  frzrdttot  = frzrdttot /real(iter)

  where ( qc .gt. 0.0)
    qcsinksum_rate1ord = qcsinksum_rate1ord/qc/real(iter)
  end where

  preo = preo/real(iter)
  prdso = prdso/real(iter)

  cmelo =cmelo/real(iter)

  eroslo=eroslo/real(iter)
  erosio=erosio/real(iter)

  pratot=pratot/real(iter)
  prctot=prctot/real(iter)
  mnuccctot=mnuccctot/real(iter)
  mnuccttot=mnuccttot/real(iter)
  mnuccdtot=mnuccdtot/real(iter)

  msacwitot=msacwitot/real(iter)
  psacwstot=psacwstot/real(iter)
  bergstot=bergstot/real(iter)
  bergtot=bergtot/real(iter)
  prcitot=prcitot/real(iter)
  praitot=praitot/real(iter)

  mnuccrtot =mnuccrtot/real(iter)
  mnuccritot=mnuccritot/real(iter)   ! h1g, 2020-02-11
  pracstot =pracstot /real(iter)

  npccno = npccno/real(iter)     ! h1g, 2020-03-09
  nprc1o = nprc1o/real(iter)     ! h1g, 2020-03-09
  nprao  = nprao/real(iter)      ! h1g, 2020-03-09
  nerosco= nerosco/real(iter)    ! h1g, 2020-03-09
  nnuccco= nnuccco/real(iter)    ! h1g, 2020-03-09
  nnuccto= nnuccto/real(iter)    ! h1g, 2020-03-09
  npsacwso= npsacwso/real(iter)  ! h1g, 2020-03-09
  nsubco  = nsubco/real(iter)    ! h1g, 2020-03-09

  nnuccdo    = nnuccdo/real(iter)    ! h1g, 2020-06-29
  nerosio    = nerosio/real(iter)    ! h1g, 2020-06-30
  nsacwio    = nsacwio/real(iter)    ! h1g, 2020-06-29
  nsubio     = nsubio/real(iter)     ! h1g, 2020-06-29
  nprcio     = nprcio/real(iter)     ! h1g, 2020-06-29
  npraio     = npraio/real(iter)     ! h1g, 2020-06-29
  nnuccrio   = nnuccrio/real(iter)   ! h1g, 2020-06-30
  nucclim1io = nucclim1io/real(iter) ! h1g, 2020-06-29



  if (diag_id%snow_inst + diag_id%snow_inst_col > 0)    &
       diag_4l(:,j,:,diag_pt%snow_inst )  = diag_4l(:,j,:,diag_pt%snow_inst )/real(iter)
  if (diag_id%snow_num_inst + diag_id%snow_num_inst_col > 0)    &
             diag_4l(:,j,:,diag_pt%snow_num_inst )  = diag_4l(:,j,:,diag_pt%snow_num_inst)/real(iter)
  if (diag_id%rain_inst + diag_id%rain_inst_col > 0)    &
             diag_4l(:,j,:,diag_pt%rain_inst )  = diag_4l(:,j,:,diag_pt%rain_inst )/real(iter)
  if (diag_id%rain_num_inst + diag_id%rain_num_inst_col > 0)    &
             diag_4l(:,j,:,diag_pt%rain_num_inst )  = diag_4l(:,j,:,diag_pt%rain_num_inst )/real(iter)





  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

  ! get new update for variables that includes sedimentation tendency
  ! note : here dum variables are grid-average, NOT in-cloud

  do k=1,nlev
     do i=1,mgncol
        dumc(i,k) = max(qc(i,k)+qctend(i,k)*deltat,0._r8)
        dumi(i,k) = max(qi(i,k)+qitend(i,k)*deltat,0._r8)
        dumnc(i,k) = max(nc(i,k)+nctend(i,k)*deltat,0._r8)
        dumni(i,k) = max(ni(i,k)+nitend(i,k)*deltat,0._r8)

        dumr(i,k) = max(qr(i,k)+qrtend(i,k)*deltat,0._r8)
        dumnr(i,k) = max(nr(i,k)+nrtend(i,k)*deltat,0._r8)
        dums(i,k) = max(qs(i,k)+qstend(i,k)*deltat,0._r8)
        dumns(i,k) = max(ns(i,k)+nstend(i,k)*deltat,0._r8)

        ! switch for specification of droplet and crystal number
        if (nccons) then
           dumnc(i,k)=ncnst/rho(i,k)*lcldm(i,k)
        end if

        ! switch for specification of cloud ice number
        if (nicons) then
           dumni(i,k)=ninst/rho(i,k)*icldm(i,k)
        end if

        if (dumc(i,k).lt.qsmall) dumnc(i,k)=0._r8
        if (dumi(i,k).lt.qsmall) dumni(i,k)=0._r8
        if (dumr(i,k).lt.qsmall) dumnr(i,k)=0._r8
        if (dums(i,k).lt.qsmall) dumns(i,k)=0._r8

     enddo
  enddo

  ! calculate instantaneous processes (melting, homogeneous freezing)
  !====================================================================

  ! melting of snow at +2 C
  do k=1,nlev

     do i=1,mgncol

        if (t(i,k)+tlat(i,k)/cpp*deltat > snowmelt) then
           if (dums(i,k) > 0._r8) then

              ! make sure melting snow doesn't reduce temperature below threshold
              dum = -xlf/cpp*dums(i,k)
              if (t(i,k)+tlat(i,k)/cpp*deltat+dum.lt. snowmelt) then
                 dum = (t(i,k)+tlat(i,k)/cpp*deltat-snowmelt)*cpp/xlf
                 dum = dum/dums(i,k)
                 dum = max(0._r8,dum)
                 dum = min(1._r8,dum)
              else
                 dum = 1._r8
              end if

              qstend(i,k)=qstend(i,k)-dum*dums(i,k)/deltat
              nstend(i,k)=nstend(i,k)-dum*dumns(i,k)/deltat
              qrtend(i,k)=qrtend(i,k)+dum*dums(i,k)/deltat
              nrtend(i,k)=nrtend(i,k)+dum*dumns(i,k)/deltat

              if (diag_id%snow_melt  + diag_id%snow_melt_col > 0) &
                  diag_4l(i,j,k, diag_pt%snow_melt) = diag_4l(i,j,k, diag_pt%snow_melt)-dum*dums(i,k)/deltat
              if (diag_id%snow_num_melt  + diag_id%snow_num_melt_col > 0) &
                  diag_4l(i,j,k, diag_pt%snow_num_melt) = diag_4l(i,j,k, diag_pt%snow_num_melt)-dum*dumns(i,k)/deltat

              dum1=-xlf*dum*dums(i,k)/deltat
              tlat(i,k)=tlat(i,k)+dum1

              meltsdttot(i,k)=meltsdttot(i,k) + dum1
           end if
        end if
     enddo
  enddo
  do k=1,nlev
      do i=1,mgncol

        ! freezing of rain at -5 C

        if (t(i,k)+tlat(i,k)/cpp*deltat < rainfrze) then

           if (dumr(i,k) > 0._r8) then

              ! make sure freezing rain doesn't increase temperature above threshold
              dum = xlf/cpp*dumr(i,k)
              if (t(i,k)+tlat(i,k)/cpp*deltat+dum.gt.rainfrze) then
                 dum = -(t(i,k)+tlat(i,k)/cpp*deltat-rainfrze)*cpp/xlf
                 dum = dum/dumr(i,k)
                 dum = max(0._r8,dum)
                 dum = min(1._r8,dum)
              else
                 dum = 1._r8
              end if

              qrtend(i,k)=qrtend(i,k)-dum*dumr(i,k)/deltat
              nrtend(i,k)=nrtend(i,k)-dum*dumnr(i,k)/deltat

              ! get mean size of rain = 1/lamr, add frozen rain to either snow or cloud ice
              ! depending on mean rain size

              call size_dist_param_basic(mg_rain_props, dumr(i,k), dumnr(i,k), &
                   lamr(i,k))

              if (lamr(i,k) < 1._r8/Dcs) then
                 qstend(i,k)=qstend(i,k)+dum*dumr(i,k)/deltat
                 nstend(i,k)=nstend(i,k)+dum*dumnr(i,k)/deltat

                 if (diag_id%srfrain_accrs + diag_id%srfrain_accrs_col > 0)    &
                     diag_4l(i,j,k,diag_pt%srfrain_accrs) = diag_4l(i,j,k,diag_pt%srfrain_accrs)-dum*dumr(i,k)/deltat

                 if (diag_id%rain_num2snow + diag_id%rain_num2snow_col > 0)    &
                     diag_4l(i,j,k,diag_pt%rain_num2snow ) = diag_4l(i,j,k,diag_pt%rain_num2snow )-dum*dumnr(i,k)/deltat

              else
                 qitend(i,k)=qitend(i,k)+dum*dumr(i,k)/deltat
                 nitend(i,k)=nitend(i,k)+dum*dumnr(i,k)/deltat
                 if (diag_id%qidt_rain2ice  + diag_id%qi_rain2ice_col > 0) &
                     diag_4l(i,j,k, diag_pt%qidt_rain2ice) = diag_4l(i,j,k, diag_pt%qidt_rain2ice)+dum*dumr(i,k)/deltat
                 if (diag_id%qnidt_rain2ice  + diag_id%qni_rain2ice_col > 0) &
                     diag_4l(i,j,k, diag_pt%qnidt_rain2ice) = diag_4l(i,j,k, diag_pt%qnidt_rain2ice)+dum*dumnr(i,k)/deltat
              end if

              ! heating tendency
              dum1 = xlf*dum*dumr(i,k)/deltat
              frzrdttot(i,k)=frzrdttot(i,k) + dum1
              tlat(i,k)=tlat(i,k)+dum1
           end if
        end if

      enddo
  enddo
   if (do_cldice) then
      do k=1,nlev
        do i=1,mgncol
           if (t(i,k)+tlat(i,k)/cpp*deltat > tmelt) then
              if (dumi(i,k) > 0._r8) then

                 ! limit so that melting does not push temperature below freezing
                 !-----------------------------------------------------------------
                 dum = -dumi(i,k)*xlf/cpp
                 if (t(i,k)+tlat(i,k)/cpp*deltat+dum.lt.tmelt) then
                    dum = (t(i,k)+tlat(i,k)/cpp*deltat-tmelt)*cpp/xlf
                    dum = dum/dumi(i,k)
                    dum = max(0._r8,dum)
                    dum = min(1._r8,dum)
                 else
                    dum = 1._r8
                 end if

                 qctend(i,k)=qctend(i,k)+dum*dumi(i,k)/deltat

                 ! for output
                 melttot(i,k)= dum*dumi(i,k)/deltat

                 ! assume melting ice produces droplet
                 ! mean volume radius of 8 micron

              IF (diag_id%qndt_melt + diag_id%qn_melt_col > 0) &
                   diag_4l(i,j,k,diag_pt%qndt_melt) = diag_4l(i,j,k,diag_pt%qndt_melt) &
                    + 3._r8*dum*dumi(i,k)/deltat/ &
                      (4._r8*pi*5.12e-16_r8*rhow)
              IF (diag_id%qidt_melt2  + diag_id%qi_melt2_col  > 0) &
                   diag_4l(i,j,k,diag_pt%qidt_melt2) =  qitend(i,k)
              IF (diag_id%qnidt_melt +  diag_id%qni_melt_col > 0) &
                   diag_4l(i,j,k,diag_pt%qnidt_melt) =  nitend(i,k)

                 nctend(i,k)=nctend(i,k)+3._r8*dum*dumi(i,k)/deltat/ &
                      (4._r8*pi*5.12e-16_r8*rhow)

                 qitend(i,k)=((1._r8-dum)*dumi(i,k)-qi(i,k))/deltat
                 nitend(i,k)=((1._r8-dum)*dumni(i,k)-ni(i,k))/deltat
                 tlat(i,k)=tlat(i,k)-xlf*dum*dumi(i,k)/deltat

              IF (diag_id%qidt_melt2  + diag_id%qi_melt2_col  > 0) &
                   diag_4l(i,j,k,diag_pt%qidt_melt2) =    &
                          qitend(i,k) - diag_4l(i,j,k,diag_pt%qidt_melt2)
              IF (diag_id%qnidt_melt +  diag_id%qni_melt_col > 0) &
                   diag_4l(i,j,k,diag_pt%qnidt_melt) =    &
                          nitend(i,k) - diag_4l(i,j,k,diag_pt%qnidt_melt)
              end if
           end if
        enddo
     enddo

     ! homogeneously freeze droplets at -40 C
     !-----------------------------------------------------------------

     do k=1,nlev
        do i=1,mgncol
           if (t(i,k)+tlat(i,k)/cpp*deltat < 233.15_r8) then
              if (dumc(i,k) > 0._r8) then

                 ! limit so that freezing does not push temperature above threshold
                 dum = dumc(i,k)*xlf/cpp
                 if (t(i,k)+tlat(i,k)/cpp*deltat+dum.gt.233.15_r8) then
                    dum = -(t(i,k)+tlat(i,k)/cpp*deltat-233.15_r8)*cpp/xlf
                    dum = dum/dumc(i,k)
                    dum = max(0._r8,dum)
                    dum = min(1._r8,dum)
                 else
                    dum = 1._r8
                 end if

                 qitend(i,k)=qitend(i,k)+dum*dumc(i,k)/deltat
                 ! for output
                 homotot(i,k)=dum*dumc(i,k)/deltat

              IF (diag_id%qldt_freez + diag_id%ql_freez_col > 0) &
                    diag_4l(i,j,k,diag_pt%qldt_freez) = qctend(i,k)
               sum_freeze(i,k) = qctend(i,k)
               IF ( diag_id%qndt_ihom  + diag_id%qn_ihom_col > 0) &
                    diag_4l(i,j,k,diag_pt%qndt_ihom) =  nctend(i,k)

                 ! assume 25 micron mean volume radius of homogeneously frozen droplets
                 ! consistent with size of detrained ice in stratiform.F90
                 nitend(i,k)=nitend(i,k)+dum*3._r8*dumc(i,k)/(4._r8*3.14_r8*1.563e-14_r8* &
                      500._r8)/deltat
                 qctend(i,k)=((1._r8-dum)*dumc(i,k)-qc(i,k))/deltat

                 if ( do_liq_num_ihom ) &   ! h1g, 2020-06-22
                 nctend(i,k)=((1._r8-dum)*dumnc(i,k)-nc(i,k))/deltat

                 tlat(i,k)=tlat(i,k)+xlf*dum*dumc(i,k)/deltat

              if( isnan( tlat(i,k) ) ) &
             write(*,'(a, 2i5, 5f8.2, i5 )') 'NaN@2907',  i, k, tlat(i,k), xlf, dum, dumc(i,k), deltat


               IF (diag_id%qldt_freez + diag_id%ql_freez_col > 0) &
                    diag_4l(i,j,k,diag_pt%qldt_freez) =    &
                          qctend(i,k) - diag_4l(i,j,k,diag_pt%qldt_freez)
               sum_freeze(i,k) = -(qctend(i,k) - sum_freeze(i,k))
               IF (diag_id%qndt_ihom  + diag_id%qn_ihom_col > 0) &
                    diag_4l(i,j,k,diag_pt%qndt_ihom) =    &
                            nctend(i,k) - diag_4l(i,j,k,diag_pt%qndt_ihom)
              IF ( diag_id%qnidt_ihom +  diag_id%qni_ihom_col > 0 ) &
                    diag_4l(i,j,k,diag_pt%qnidt_ihom) =    &
                           dum*3._r8*dumc(i,k)/   &
                                (4._r8*3.14_r8*1.563e-14_r8*500._r8)/deltat
              end if
           end if
        enddo
     enddo
     ! remove any excess over-saturation, which is possible due to non-linearity when adding
     ! together all microphysical processes
     !-----------------------------------------------------------------
     ! follow code similar to old CAM scheme

     do k=1,nlev
        do i=1,mgncol

           qtmp=q(i,k)+qvlat(i,k)*deltat
           ttmp=t(i,k)+tlat(i,k)/cpp*deltat

           if ( ttmp .lt.-150.0+273.15 .or. ttmp .gt.90+273.15)  &
 write(*,'(a,6f8.2)') 'MG2: bad temperature@2930',   &
            ttmp, t(i,k), tlat(i,k)/cpp*deltat, tlat(i,k), cpp, deltat

          ! use rhw to allow ice supersaturation
           call compute_qs(ttmp, p(i,k), qvn, q = q(i,k),   &
                             esat = esn, es_over_liq = .true.)

           if (qtmp > qvn .and. qvn > 0 .and. allow_sed_supersat) then
              ! expression below is approximate since there may be ice deposition
              dum = (qtmp-qvn)/(1._r8+xxlv_squared*qvn/(cpp*rv*ttmp**2))/deltat
              ! add to output cme
              cmeout(i,k) = cmeout(i,k)+dum
              ! now add to tendencies, partition between liquid and ice based on temperature

              if( remove_super_RK ) then    ! h1g, 2020-03-30
                if (ttmp > 233.15_r8) then
                   dum1=0.0_r8
                   if ( tiedtke_macrophysics ) ssat_disposal(i,k) = 1._r8
                else
                   dum1=1.0_r8
                   if ( tiedtke_macrophysics ) ssat_disposal(i,k) = 2._r8
                endif
              else  ! h1g, 2020-03-30

                if (ttmp > 268.15_r8) then
                   dum1=0.0_r8
                   ! now add to tendencies, partition between liquid and ice based on temperature
                   if ( tiedtke_macrophysics ) ssat_disposal(i,k) = 1._r8
                   !-------------------------------------------------------
                else if (ttmp < 238.15_r8) then
                   dum1=1.0_r8
                   if ( tiedtke_macrophysics ) ssat_disposal(i,k) = 2._r8
                else
                   dum1=(268.15_r8-ttmp)/30._r8
                   if ( tiedtke_macrophysics ) ssat_disposal(i,k) = 2._r8
                end if
              endif  ! h1g, 2020-03-30

              dum = (qtmp-qvn)/(1._r8+(xxls*dum1+xxlv*(1._r8-dum1))**2 &
                   *qvn/(cpp*rv*ttmp**2))/deltat
              qctend(i,k)=qctend(i,k)+dum*(1._r8-dum1)
              ! for output
              qcrestot(i,k)=dum*(1._r8-dum1)
              qitend(i,k)=qitend(i,k)+dum*dum1
              qirestot(i,k)=dum*dum1
              qvlat(i,k)=qvlat(i,k)-dum
              ! for output
              qvres(i,k)=-dum
              tlat(i,k)=tlat(i,k)+dum*(1._r8-dum1)*xxlv+dum*dum1*xxls
           end if
        enddo
     enddo
  end if

  ! calculate effective radius for pass to radiation code
  !=========================================================
  ! if no cloud water, default value is 10 micron for droplets,
  ! 25 micron for cloud ice

  ! update cloud variables after instantaneous processes to get effective radius
  ! variables are in-cloud to calculate size dist parameters
  do k=1,nlev
     do i=1,mgncol
        dumc(i,k) = max(qc(i,k)+qctend(i,k)*deltat,0._r8)/lcldm(i,k)
        dumi(i,k) = max(qi(i,k)+qitend(i,k)*deltat,0._r8)/icldm(i,k)
        dumnc(i,k) = max(nc(i,k)+nctend(i,k)*deltat,0._r8)/lcldm(i,k)
        dumni(i,k) = max(ni(i,k)+nitend(i,k)*deltat,0._r8)/icldm(i,k)

        dumr(i,k) = max(qr(i,k)+qrtend(i,k)*deltat,0._r8)/precip_frac(i,k)
        dumnr(i,k) = max(nr(i,k)+nrtend(i,k)*deltat,0._r8)/precip_frac(i,k)
        dums(i,k) = max(qs(i,k)+qstend(i,k)*deltat,0._r8)/precip_frac(i,k)
        dumns(i,k) = max(ns(i,k)+nstend(i,k)*deltat,0._r8)/precip_frac(i,k)

        ! switch for specification of droplet and crystal number
        if (nccons) then
           dumnc(i,k)=ncnst/rho(i,k)
        end if

        ! switch for specification of cloud ice number
        if (nicons) then
           dumni(i,k)=ninst/rho(i,k)
        end if

        ! limit in-cloud mixing ratio to reasonable value of 5 g kg-1
        dumc(i,k)=min(dumc(i,k),5.e-3_r8)
        dumi(i,k)=min(dumi(i,k),5.e-3_r8)
        ! limit in-precip mixing ratios
        dumr(i,k)=min(dumr(i,k),10.e-3_r8)
        dums(i,k)=min(dums(i,k),10.e-3_r8)
     enddo
  enddo
  ! cloud ice effective radius
  !-----------------------------------------------------------------

  if (do_cldice) then
     do k=1,nlev
        do i=1,mgncol
           if (dumi(i,k).ge.qsmall) then

              dum_2D(i,k) = dumni(i,k)
              call size_dist_param_basic(mg_ice_props, dumi(i,k), dumni(i,k), &
                   lami(i,k), dumni0)

              if (dumni(i,k) /=dum_2D(i,k)) then
                 ! adjust number conc if needed to keep mean size in reasonable range
                  if (diag_id%qnidt_size_adj + diag_id%qni_size_adj_col  > 0)    &
                      diag_4l(i,j,k,diag_pt%qnidt_size_adj)  = nitend(i,k)

                 if ( do_ice_num_adjust ) &  ! h1g, 2020-07-01
                 nitend(i,k)=(dumni(i,k)*icldm(i,k)-ni(i,k))/deltat

                  if (diag_id%qnidt_size_adj + diag_id%qni_size_adj_col  > 0)    &
                      diag_4l(i,j,k,diag_pt%qnidt_size_adj)  = nitend(i,k)   &
                                                             - diag_4l(i,j,k,diag_pt%qnidt_size_adj)
              end if

              effi(i,k) = 1.5_r8/lami(i,k)*1.e6_r8
              sadice(i,k) = 2._r8*pi*(lami(i,k)**(-3))*dumni0*rho(i,k)*1.e-2_r8  ! m2/m3 -> cm2/cm3

           else
              effi(i,k) = 25._r8
              sadice(i,k) = 0._r8
           end if

           ! ice effective diameter for david mitchell's optics
           deffi(i,k)=effi(i,k)*rhoi/rhows*2._r8
        enddo
     enddo
  else
     do k=1,nlev
        do i=1,mgncol
           ! NOTE: If CARMA is doing the ice microphysics, then the ice effective
           ! radius has already been determined from the size distribution.
           effi(i,k) = re_ice(i,k) * 1.e6_r8      ! m -> um
           deffi(i,k)=effi(i,k) * 2._r8
           sadice(i,k) = 4._r8*pi*(effi(i,k)**2)*ni(i,k)*rho(i,k)*1e-2_r8
        enddo
     enddo
  end if

  ! cloud droplet effective radius
  !-----------------------------------------------------------------
  do k=1,nlev
     do i=1,mgncol
        if (dumc(i,k).ge.qsmall  .and. do_liq_num_adjust) then

           ! switch for specification of droplet and crystal number
           if (nccons) then
              ! make sure nc is consistence with the constant N by adjusting tendency, need
              ! to multiply by cloud fraction
              ! note that nctend may be further adjusted below if mean droplet size is
              ! out of bounds
              nctend(i,k)=(ncnst/rho(i,k)*lcldm(i,k)-nc(i,k))/deltat
           end if
           IF (diag_id%qndt_size_adj + diag_id%qn_size_adj_col  > 0) &
                 diag_4l(i,j,k,diag_pt%qndt_size_adj ) =  nctend(i,k)

           dum = dumnc(i,k)
           call size_dist_param_liq(mg_liq_props, dumc(i,k), dumnc(i,k), rho(i,k), &
                pgam(i,k), lamc(i,k))

           if (dum /= dumnc(i,k)) then
              ! adjust number conc if needed to keep mean size in reasonable range
              nctend(i,k)=(dumnc(i,k)*lcldm(i,k)-nc(i,k))/deltat
           end if
            IF (diag_id%qndt_size_adj + diag_id%qn_size_adj_col  > 0) &
                 diag_4l(i,j,k,diag_pt%qndt_size_adj ) =     &
                        nctend(i,k) - diag_4l(i,j,k,diag_pt%qndt_size_adj )

           effc(i,k) = (pgam(i,k)+3._r8)/lamc(i,k)/2._r8*1.e6_r8
           !assign output fields for shape here
           lamcrad(i,k)=lamc(i,k)
           pgamrad(i,k)=pgam(i,k)


           ! recalculate effective radius for constant number, in order to separate
           ! first and second indirect effects
           !======================================
           ! assume constant number of 10^8 kg-1

           dumnc(i,k)=1.e8_r8

           ! Pass in "false" adjust flag to prevent number from being changed within
           ! size distribution subroutine.
           call size_dist_param_liq(mg_liq_props, dumc(i,k), dumnc(i,k), rho(i,k), &
                pgam(i,k), lamc(i,k))

           effc_fn(i,k) = (pgam(i,k)+3._r8)/lamc(i,k)/2._r8*1.e6_r8

        else
           effc(i,k) = 10._r8
           lamcrad(i,k)=0._r8
           pgamrad(i,k)=0._r8
           effc_fn(i,k) = 10._r8
        end if
     enddo
  enddo
  ! recalculate 'final' rain size distribution parameters
  ! to ensure that rain size is in bounds, adjust rain number if needed
  do k=1,nlev
     do i=1,mgncol

        if (dumr(i,k).ge.qsmall) then

           dum = dumnr(i,k)

           call size_dist_param_basic(mg_rain_props, dumr(i,k), dumnr(i,k), &
                lamr(i,k))

           if (dum /= dumnr(i,k)) then
             if (diag_id%rain_num_adj + diag_id%rain_num_adj_col > 0)    &
              diag_4l(:,j,:,diag_pt%rain_num_adj)  = nrtend(i,k)

              ! adjust number conc if needed to keep mean size in reasonable range
              nrtend(i,k)=(dumnr(i,k)*precip_frac(i,k)-nr(i,k))/deltat

             if (diag_id%rain_num_adj + diag_id%rain_num_adj_col > 0)    &
              diag_4l(:,j,:,diag_pt%rain_num_adj)  = nrtend(i,k)-diag_4l(:,j,:,diag_pt%rain_num_adj)

           end if

        end if
     enddo
  enddo
  ! recalculate 'final' snow size distribution parameters
  ! to ensure that snow size is in bounds, adjust snow number if needed
  do k=1,nlev
     do i=1,mgncol
        if (dums(i,k).ge.qsmall) then

           dum = dumns(i,k)

           call size_dist_param_basic(mg_snow_props, dums(i,k), dumns(i,k), &
                lams(i,k), n0=dumns0)

           if (dum /= dumns(i,k)) then
             if (diag_id%snow_num_adj + diag_id%snow_num_adj_col > 0)    &
              diag_4l(:,j,:,diag_pt%snow_num_adj)  = nstend(i,k)

              ! adjust number conc if needed to keep mean size in reasonable range
              nstend(i,k)=(dumns(i,k)*precip_frac(i,k)-ns(i,k))/deltat

             if (diag_id%snow_num_adj + diag_id%snow_num_adj_col > 0)    &
              diag_4l(:,j,:,diag_pt%snow_num_adj)  = nstend(i,k)-diag_4l(:,j,:,diag_pt%snow_num_adj)
           end if

           sadsnow(i,k) = 2._r8*pi*(lams(i,k)**(-3))*dumns0*rho(i,k)*1.e-2_r8  ! m2/m3 -> cm2/cm3

        end if


     end do ! vertical k loop
  enddo

  ! DO STUFF FOR OUTPUT:
  !==================================================
  ! averaging for snow and rain number and diameter
  !--------------------------------------------------

  ! drout2/dsout2:
  ! diameter of rain and snow
  ! dsout:
  ! scaled diameter of snow (passed to radiation in CAM)
  ! reff_rain/reff_snow:
  ! calculate effective radius of rain and snow in microns for COSP using Eq. 9 of COSP v1.3 manual

  where (qrout .gt. 1.e-7_r8 &
       .and. nrout.gt.0._r8)
     qrout2 = qrout * precip_frac
     nrout2 = nrout * precip_frac
     ! The avg_diameter call does the actual calculation; other diameter
     ! outputs are just drout2 times constants.
     drout2 = avg_diameter(qrout, nrout, rho, rhow)
     freqr = precip_frac

     reff_rain=1.5_r8*drout2*1.e6_r8
  elsewhere
     qrout2 = 0._r8
     nrout2 = 0._r8
     drout2 = 0._r8
     freqr = 0._r8
     reff_rain = 0._r8
  end where

  where (qsout .gt. 1.e-7_r8 &
       .and. nsout.gt.0._r8)
     qsout2 = qsout * precip_frac
     nsout2 = nsout * precip_frac
     ! The avg_diameter call does the actual calculation; other diameter
     ! outputs are just dsout2 times constants.
     dsout2 = avg_diameter(qsout, nsout, rho, rhosn)
     freqs = precip_frac

     dsout=3._r8*rhosn/rhows*dsout2

     reff_snow=1.5_r8*dsout2*1.e6_r8
  elsewhere
     dsout  = 0._r8
     qsout2 = 0._r8
     nsout2 = 0._r8
     dsout2 = 0._r8
     freqs  = 0._r8
     reff_snow=0._r8
  end where

!--> h1g, 2010-01-15, add limits for rain drop radius
  reff_rain          = max(  30.0_r8, reff_rain          )
  reff_rain          = min( 750.0_r8, reff_rain          )
!<-- h1g, 2010-01-15, add limits for rain drop radius



! diagnostics for water tendencies
! water  vapor specific humicity

 !     if (diag_id%qdt_cond   > 0) &
 !             diag_4l(:,j,:,diag_pt%qdt_cond)  = -cmelo(:,:)

      if  (diag_id%qdt_deposition > 0)  &
              diag_4l(:,j,:,diag_pt%qdt_deposition )  = -cmeitot( : , : )
      if  (diag_id%qdt_sedi_ice2vapor> 0)  &
              diag_4l(:,j,:,diag_pt%qdt_sedi_ice2vapor) = qisevap( : , : )
      if  (diag_id%qdt_sedi_liquid2vapor> 0)  &
            diag_4l(:,j,:,diag_pt%qdt_sedi_liquid2vapor) = qcsevap( : , : )
      if  (diag_id%qdt_super_sat_rm > 0)  &
              diag_4l(:,j,:,diag_pt%qdt_super_sat_rm) = qvres( : , : )

! cloud liquid water
      if (diag_id%qldt_accr  + diag_id%ql_accr_col > 0) &
              diag_4l(:,j,:,diag_pt%qldt_accr)  = - pratot(:,:)
      if (diag_id%qldt_auto  + diag_id%ql_auto_col > 0)&
              diag_4l(:,j,:,diag_pt%qldt_auto)  = -prctot(:,:)
      if (diag_id%qldt_freez2 + diag_id%ql_freez2_col > 0) &
              diag_4l(:,j,:,diag_pt%qldt_freez2) =   &
                                          -(mnuccctot(:,:) + mnuccttot(:,:) )
      sum_freeze2(:,:) =  mnuccctot(:,:) + mnuccttot(:,:)
      if (diag_id%qldt_accrs  + diag_id%ql_accrs_col > 0) &
              diag_4l(:,j,:,diag_pt%qldt_accrs)  = -psacwstot(:,:)
      sum_rime(:,:) =  psacwstot(:,:)
      if (diag_id%qldt_HM_splinter + diag_id%ql_HM_splinter_col > 0)&
              diag_4l(:,j,:,diag_pt%qldt_HM_splinter)  = -msacwitot(:,:)
      sum_splinter(:,:) =  msacwitot(:,:)
      if (diag_id%qldt_bergs + diag_id%ql_bergs_col > 0) &
              diag_4l(:,j,:,diag_pt%qldt_bergs)  = -bergstot(:,:)
      sum_bergs(:,:) =  bergstot (:,:)


      if (diag_id%qidt_dep + diag_id%qi_dep_col > 0)    &
              diag_4l(:,j,:,diag_pt%qidt_dep)  = max(cmeitot(:,:),0._r8)
      if (diag_id%qidt_subl + diag_id%qi_subl_col > 0)   &
              diag_4l(:,j,:,diag_pt%qidt_subl) = - max(-1._r8*cmeitot(:,:),0._r8)

      sum_cond(:,:) = max(cmeitot(:,:),0._r8)

      if (diag_id%qldt_cond + diag_id%ql_cond_col > 0) &
              diag_4l(:,j,:,diag_pt%qldt_cond)  =  max(cmelo(:,:), 0._r8)
      if (diag_id%qldt_evap  + diag_id%ql_evap_col > 0) &
              diag_4l(:,j,:,diag_pt%qldt_evap)  =   &
                                           - max(-1._r8*cmelo(:,:),0._r8)


      if (diag_id%qldt_eros + diag_id%ql_eros_col > 0) &
              diag_4l(:,j,:,diag_pt%qldt_eros)  =   eroslo(:,:)
      if (diag_id%qdt_eros_l                       > 0) &
              diag_4l(:,j,:,diag_pt%qdt_eros_l)  = -eroslo(:,:)

      if (diag_id%qidt_eros + diag_id%qi_eros_col > 0) &
              diag_4l(:,j,:,diag_pt%qidt_eros)  =   erosio(:,:)
      if (diag_id%qdt_eros_i                       > 0) &
              diag_4l(:,j,:,diag_pt%qdt_eros_i)  = -erosio(:,:)


      if (diag_id%qldt_berg + diag_id%ql_berg_col > 0) &
              diag_4l(:,j,:,diag_pt%qldt_berg)  =  -bergtot(:,:)
      sum_berg(:,:) =  bergtot(:,:)
      IF ( diag_id%qldt_sedi  + diag_id%ql_sedi_col > 0 ) &
              diag_4l(:,j,1:nlev,diag_pt%qldt_sedi) = qcsedten(:,:)
      IF ( diag_id%liq_adj  + diag_id%liq_adj_col > 0 ) &
              diag_4l(:,j,:,diag_pt%liq_adj) = qcrestot(:,:)

! cloud ice water
      if (diag_id%qidt_auto + diag_id%qi_auto_col > 0) &
             diag_4l(:,j,:,diag_pt%qidt_auto) = -prcitot(:,:)
      if (diag_id%qidt_accr  + diag_id%qi_accr_col > 0) &
             diag_4l(:,j,:,diag_pt%qidt_accr) = -praitot(:,:)
      if (diag_id%qidt_rain2ice  + diag_id%qi_rain2ice_col > 0) &
             diag_4l(:,j,:,diag_pt%qidt_rain2ice) =  diag_4l(:,j,:,diag_pt%qidt_rain2ice) &
                                                    + mnuccritot(:,:)
      IF ( diag_id%qidt_fall  + diag_id%qi_fall_col > 0 ) &
              diag_4l(:,j,1:nlev,diag_pt%qidt_fall) = qisedten(:,1:nlev)
      IF ( diag_id%ice_adj  +  diag_id%ice_adj_col > 0 ) &
              diag_4l(:,j,:,diag_pt%ice_adj) = qirestot(:,:)
      sum_ice_adj(:,:) = qirestot(:,:)

! ---> rain water mixing ratio
      if (diag_id%srfrain_accrs + diag_id%srfrain_accrs_col > 0)    &
             diag_4l(:,j,:,diag_pt%srfrain_accrs)  =  diag_4l(:,j,:,diag_pt%srfrain_accrs) &
                                                     -pracstot(:,:)
      if (diag_id%rain_freeze + diag_id%rain_freeze_col > 0)    &
             diag_4l(:,j,:,diag_pt%rain_freeze)  = -mnuccrtot(:,:)
     if  (diag_id%rain_evap + diag_id%rain_evap_col > 0)  &
             diag_4l(:,j,:,diag_pt%rain_evap)  =  pre( : , : )*precip_frac(:,:)
     if  (diag_id%rain_sedi + diag_id%rain_sedi_col > 0)  &
             diag_4l(:,j,:,diag_pt%rain_sedi)  = qrsedten(:,:)

! ---> snow mixing ratio
     if  (diag_id%qdt_snow_sublim + diag_id%q_snow_sublim_col > 0 )  &
             diag_4l(:,j,:,diag_pt%qdt_snow_sublim )  =  -prds( : , : )*precip_frac(:,:)
     if  (diag_id%snow_sedi + diag_id%snow_sedi_col > 0)  &
             diag_4l(:,j,:,diag_pt%snow_sedi)  = qssedten(:,:)

! ---> rain number mixing ratio
      if ( diag_id%rain_num2snow + diag_id%rain_num2snow_col > 0 )    &
             diag_4l(:,j,:,diag_pt%rain_num2snow)  =  diag_4l(:,j,:,diag_pt%rain_num2snow) &
                                                     -npracs(:,:)*precip_frac(:,:)
      if ( diag_id%rain_num_evap + diag_id%rain_num_evap_col > 0 )    &
             diag_4l(:,j,:,diag_pt%rain_num_evap)  =  nsubr( : , : )*precip_frac(:,:)

      if ( diag_id%rain_num_freez + diag_id%rain_num_freez_col > 0 )    &
             diag_4l(:,j,:,diag_pt%rain_num_freez)  =  -nnuccr( : , : )*precip_frac(:,:)


! --->liquid droplet number
!      if ( diag_id%qndt_cond  + diag_id%qn_cond_col  > 0 ) &
!            diag_4l(:,j,:,diag_pt%qndt_cond)  = npccn2(:,:)
      if (diag_id%qndt_cond + diag_id%qn_cond_col > 0)    &
             diag_4l(:,j,:,diag_pt%qndt_cond)  = npccno(:,:)
      if (diag_id%qndt_eros + diag_id%qn_eros_col > 0)    &
             diag_4l(:,j,:,diag_pt%qndt_eros)  = nerosco(:,:)
      if (diag_id%qndt_pra + diag_id%qn_pra_col > 0)    &
             diag_4l(:,j,:,diag_pt%qndt_pra)  = nprao(:,:)
      if (diag_id%qndt_auto + diag_id%qn_auto_col > 0)    &
              diag_4l(:,j,:,diag_pt%qndt_auto)  = nprc1o(:,:)
      if (diag_id%qndt_freez + diag_id%qn_freez_col > 0)    &
             diag_4l(:,j,:,diag_pt%qndt_freez)  = nnuccco(:,:)
      if (diag_id%qndt_contact_frz + diag_id%qn_contact_frz_col > 0)    &
             diag_4l(:,j,:,diag_pt%qndt_contact_frz)  = nnuccto(:,:)
      if (diag_id%qndt_sacws + diag_id%qn_sacws_col > 0)    &
             diag_4l(:,j,:,diag_pt%qndt_sacws)  = npsacwso(:,:)
      if (diag_id%qndt_evap + diag_id%qn_evap_col > 0)    &
             diag_4l(:,j,:,diag_pt%qndt_evap)  = nsubco(:,:)

! ---> ice number
      if (diag_id%qnidt_nnuccd +  diag_id%qni_nnuccd_col > 0)    &
             diag_4l(:,j,:,diag_pt%qnidt_nnuccd)  =  nnuccdo(:,:)
      if (diag_id%qnidt_nsacwi> 0)    &
             diag_4l(:,j,:,diag_pt%qnidt_nsacwi)  =  nsacwio

      if (diag_id%qnidt_nsubi  + diag_id%qni_nsubi_col  > 0)    &
             diag_4l(:,j,:,diag_pt%qnidt_nsubi)  = nsubio(:,:)
     if (diag_id%qnidt_nerosi  + diag_id%qni_nerosi_col  > 0)    &
             diag_4l(:,j,:,diag_pt%qnidt_nerosi)  = nerosio
      if (diag_id%qnidt_auto  + diag_id%qni_auto_col  > 0)    &
             diag_4l(:,j,:,diag_pt%qnidt_auto)  = nprcio
      if (diag_id%qnidt_accr  + diag_id%qni_accr_col  > 0)    &
             diag_4l(:,j,:,diag_pt%qnidt_accr)  = npraio
      if (diag_id%qnidt_rain2ice  + diag_id%qni_rain2ice_col > 0) &
             diag_4l(:,j,:,diag_pt%qnidt_rain2ice) =  diag_4l(:,j,:,diag_pt%qnidt_rain2ice) &
                                                    +nnuccrio

!10/23/13: NOTE: STILL NEEDS CONVERSION !!!!
!RSH:
!   calculate fraction of total ice / snow creation that requires
!   ice-forming nuclei
      do k=1,nlev
        do i=1,mgncol
          qldt_sum = sum_cond(i,k) + sum_rime(i,k) + sum_berg(i,k) + &
                     sum_ice_adj(i,k) + MAX(sum_bergs(i,k), 0.0) + &
                     sum_freeze(i,k) + sum_freeze2(i,k) + sum_splinter(i,k)
          if ( ABS(qldt_sum) > 0.0            ) then
! ---> h1g, 2014-07-18, add option of including contact freeze in bergeron
               if (include_homogeneous_for_wetdep) then
           
                  if( include_contact_freeze_in_berg ) then
                    f_snow_berg(i,k) = (sum_berg(i,k) + sum_cond(i,k) +   &
                                  sum_ice_adj(i,k) +    &
                                  MAX( sum_bergs(i,k), 0.0) +     &
                                  sum_freeze (i,k) + sum_freeze2(i,k) )/qldt_sum
                  else
                    f_snow_berg(i,k) = (sum_berg(i,k) + sum_cond(i,k) +   &
                                  sum_ice_adj(i,k) +    &
                                  MAX( sum_bergs(i,k), 0.0) +     &
                                  sum_freeze (i,k) + mnuccctot(i,k) )/qldt_sum  ! h1g 2015-06-05

                  endif
               else
                   f_snow_berg(i,k) = (sum_berg(i,k) + MAX( sum_bergs(i,k), 0.0))/qldt_sum  ! h1g 2024-01-31
                 
               endif
! <--- h1g, 2014-07-18
          else
            f_snow_berg(i,k) = 0._r8
          endif
        end do
      end do

end subroutine micro_mg2_tend


!========================================================================
!OUTPUT CALCULATIONS
!========================================================================

subroutine calc_rercld(lamr, n0r, lamc, pgam, qric, qcic, ncic, rercld, mgncol)
  integer, intent(in) :: mgncol
  real(r8), dimension(mgncol), intent(in) :: lamr          ! rain size parameter (slope)
  real(r8), dimension(mgncol), intent(in) :: n0r           ! rain size parameter (intercept)
  real(r8), dimension(mgncol), intent(in) :: lamc          ! size distribution parameter (slope)
  real(r8), dimension(mgncol), intent(in) :: pgam          ! droplet size parameter
  real(r8), dimension(mgncol), intent(in) :: qric          ! in-cloud rain mass mixing ratio
  real(r8), dimension(mgncol), intent(in) :: qcic          ! in-cloud cloud liquid
  real(r8), dimension(mgncol), intent(in) :: ncic          ! in-cloud droplet number concentration

  real(r8), dimension(mgncol), intent(inout) :: rercld     ! effective radius calculation for rain + cloud

  ! combined size of precip & cloud drops
  real(r8) :: Atmp

  integer :: i

  do i=1,mgncol
     ! Rain drops
     if (lamr(i) > 0._r8) then
        Atmp = n0r(i) * pi / (2._r8 * lamr(i)**3._r8)
     else
        Atmp = 0._r8
     end if

     ! Add cloud drops
     if (lamc(i) > 0._r8) then
        Atmp = Atmp + &
             ncic(i) * pi * rising_factorial(pgam(i)+1._r8, 2)/(4._r8 * lamc(i)**2._r8)
     end if

     if (Atmp > 0._r8) then
        rercld(i) = rercld(i) + 3._r8 *(qric(i) + qcic(i)) / (4._r8 * rhow * Atmp)
     end if
  enddo
end subroutine calc_rercld


!========================================================================
!UTILITIES
!========================================================================

  subroutine micro_mg2_get_cols(ncol, nlev, top_lev, qcn, qin, &
     qrn, qsn, mgncol, mgcols)

  ! Determines which columns microphysics should operate over by
  ! checking for non-zero cloud water/ice.

  integer, intent(in) :: ncol      ! Number of columns with meaningful data
  integer, intent(in) :: nlev      ! Number of levels to use
  integer, intent(in) :: top_lev   ! Top level for microphysics

  real(r8), intent(in) :: qcn(:,:) ! cloud water mixing ratio (kg/kg)
  real(r8), intent(in) :: qin(:,:) ! cloud ice mixing ratio (kg/kg)
  real(r8), intent(in) :: qrn(:,:) ! rain mixing ratio (kg/kg)
  real(r8), intent(in) :: qsn(:,:) ! snow mixing ratio (kg/kg)

  integer, intent(out) :: mgncol   ! Number of columns MG will use
  integer, allocatable, intent(out) :: mgcols(:) ! column indices

  integer :: lev_offset  ! top_lev - 1 (defined here for consistency)
  logical :: ltrue(ncol) ! store tests for each column

  integer :: i, ii ! column indices

  if (allocated(mgcols)) deallocate(mgcols)

  lev_offset = top_lev - 1

  ! Using "any" along dimension 2 collapses across levels, but
  ! not columns, so we know if water is present at any level
  ! in each column.

  ltrue = any(qcn(:ncol,top_lev:(nlev+lev_offset)) >= qsmall, 2)
  ltrue = ltrue .or. any(qin(:ncol,top_lev:(nlev+lev_offset)) >= qsmall, 2)
  ltrue = ltrue .or. any(qrn(:ncol,top_lev:(nlev+lev_offset)) >= qsmall, 2)
  ltrue = ltrue .or. any(qsn(:ncol,top_lev:(nlev+lev_offset)) >= qsmall, 2)

  ! Scan for true values to get a usable list of indices.

  mgncol = count(ltrue)
  allocate(mgcols(mgncol))
  i = 0
  do ii = 1,ncol
     if (ltrue(ii)) then
        i = i + 1
        mgcols(i) = ii
     end if
  end do

end subroutine micro_mg2_get_cols




! =======================================================================
! time - implicit monotonic scheme
! developed by sj lin, 2016
! =======================================================================

subroutine implicit_fall (dt, ktop, kbot, ze, vt, dp, q, precip, m1)

    implicit none

    integer, intent (in) :: ktop, kbot

    real(r8), intent (in) :: dt

    real(r8), intent (in), dimension (ktop:kbot + 1) :: ze

    real(r8), intent (in), dimension (ktop:kbot) :: vt, dp

    real(r8), intent (inout), dimension (ktop:kbot) :: q

    real(r8), intent (out), dimension (ktop:kbot) :: m1

    real(r8), intent (out) :: precip

    real(r8), dimension (ktop:kbot) :: dz, qm, dd

    integer :: k

    do k = ktop, kbot
        dz (k) = ze (k) - ze (k + 1)
        dd (k) = dt * vt (k)
        q (k) = q (k) * dp (k)
    enddo

    ! -----------------------------------------------------------------------
    ! sedimentation: non - vectorizable loop
    ! -----------------------------------------------------------------------

    qm (ktop) = q (ktop) / (dz (ktop) + dd (ktop))
    do k = ktop + 1, kbot
        qm (k) = (q (k) + dd (k - 1) * qm (k - 1)) / (dz (k) + dd (k))
    enddo

    ! -----------------------------------------------------------------------
    ! qm is density at this stage
    ! -----------------------------------------------------------------------

    do k = ktop, kbot
        qm (k) = qm (k) * dz (k)
    enddo

    ! -----------------------------------------------------------------------
    ! output mass fluxes: non - vectorizable loop
    ! -----------------------------------------------------------------------

    m1 (ktop) = q (ktop) - qm (ktop)
    do k = ktop + 1, kbot
        m1 (k) = m1 (k - 1) + q (k) - qm (k)
    enddo
    precip = m1 (kbot)

    ! -----------------------------------------------------------------------
    ! update:
    ! -----------------------------------------------------------------------

    do k = ktop, kbot
        q (k) = qm (k) / dp (k)
    enddo

end subroutine implicit_fall




end module micro_mg2_mod
