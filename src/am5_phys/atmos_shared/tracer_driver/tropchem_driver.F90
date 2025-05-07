module tropchem_driver_mod
!
! <CONTACT EMAIL="Larry.Horowitz@noaa.gov">
!   Larry W. Horowitz
! </CONTACT>

! <OVERVIEW>
!     This code calculates tracer tendencies due to tropospheric chemistry
! </OVERVIEW>

! <DESCRIPTION>
!
! This code calculates chemical production and loss of tracers due
! to tropospheric chemistry. It also includes dry deposition, upper
! boundary conditions, emissions. Off-line sulfate concentrations are
! read in for use in calculating heterogeneous reaction rates (if SO4
! is not included as a tracer).
!
! This module is only activated if do_tropchem=T in tropchem_driver_nml
!
! </DESCRIPTION>


!-----------------------------------------------------------------------

use                    mpp_mod, only : input_nml_file
use                fms2_io_mod, only : file_exists, open_file, close_file, FmsNetcdfDomainFile_t, variable_exists
use            mpp_domains_mod, only : domain2D
use                    fms_mod, only : write_version_number, &
                                       mpp_pe,  &
                                       mpp_root_pe, &
                                       lowercase,   &
                                       uppercase, &
                                       stdlog, &
                                       check_nml_error, &
                                       error_mesg, &
                                       FATAL, &
                                       WARNING, &
                                       NOTE
use         tropchem_types_mod, only : tropchem_opt, tropchem_diag, tropchem_types_init, missing_value
use           time_manager_mod, only : time_type,print_date, &
                                       get_date, &
                                       set_date, &
                                       set_time, &
                                       days_in_year, &
                                       days_in_month, &
                                       real_to_time_type, &
                                       time_type_to_real, &
                                       operator(+), operator(-)
use           diag_manager_mod, only : send_data,            &
                                       register_diag_field,  &
                                       register_static_field, &
                                       get_base_time
use        atmos_cmip_diag_mod, only : register_cmip_diag_field_2d
use         tracer_manager_mod, only : get_tracer_index,     &
                                       get_tracer_names,     &
                                       get_number_tracers,   &
                                       query_method,         &
                                       check_if_prognostic,  &
                                       NO_TRACER
use          field_manager_mod, only : MODEL_ATMOS, MODEL_LAND, parse
use atmos_tracer_utilities_mod, only : dry_deposition, get_chem_param
use              constants_mod, only : grav, rdgas, WTMAIR, WTMH2O, AVOGNO, &
                                       PI, DEG_TO_RAD, SECONDS_PER_DAY
use                    mpp_mod, only : mpp_clock_id,         &
                                       mpp_clock_begin,      &
                                       mpp_clock_end
use           interpolator_mod, only : interpolate_type,     &
                                       interpolate_type_eq,  &
                                       interpolator_init,    &
                                       obtain_interpolator_time_slices, &
                                       unset_interpolator_time_flag, &
                                       interpolator_end,     &
                                       interpolator,         &
                                       query_interpolator,   &
                                       init_clim_diag,       &
                                       CONSTANT,             &
                                       INTERP_WEIGHTED_P
use            time_interp_mod, only : time_interp_init, time_interp
use              mo_chemdr_mod, only : chemdr, chemdr_init
use              mo_setsox_mod, only : setsox_init
use             mo_chemini_mod, only : chemini
use             M_TRACNAME_MOD, only : tracnam
use                MO_GRID_MOD, only : pcnstm1
use              CHEM_MODS_MOD, only : phtcnt, gascnt
use               MOZ_HOOK_MOD, only : moz_hook_init
use   strat_chem_utilities_mod, only : strat_chem_utilities_init, &
                                       strat_chem_dcly_dt, &
                                       strat_chem_dcly_dt_time_vary, &
                                       strat_chem_dcly_dt_endts, &
                                       strat_chem_get_aerosol, &
                                       psc_type, &
                                       strat_chem_get_h2so4, &
                                       strat_chem_get_psc, &
                                       strat_chem_destroy_psc, &
                                       strat_chem_psc_sediment, &
                                       strat_chem_get_extra_h2o
use           mo_chem_utls_mod, only : get_spc_ndx
use          atmos_sulfate_mod, only : atmos_sulfate_init, &
                                       atmos_sulfate_time_vary, &
                                       atmos_DMS_emission
use       shortwave_driver_mod, only : shortwave_number_of_bands, &
                                       get_solar_flux_by_band
use astronomy_mod,         only : diurnal_solar, universal_time
use horiz_interp_mod, only: horiz_interp_type, horiz_interp_init, &
                            horiz_interp_new, horiz_interp

use cloud_chem, only: CLOUD_CHEM_PH_LEGACY, CLOUD_CHEM_PH_BISECTION, &
                      CLOUD_CHEM_PH_CUBIC, CLOUD_CHEM_F1P
use aerosol_thermodynamics, only: AERO_ISORROPIA, AERO_LEGACY, NO_AERO
use mo_usrrxt_mod, only: HET_CHEM_LEGACY, HET_CHEM_J1M, &
                         GSO2_WANG2014, GSO2_ZHENG2015, &
                         GSO2_ZHENG2015_LOW
use mo_chem_utls_mod, only : get_rxt_ndx

use atmos_cmip_diag_mod,   only : register_cmip_diag_field_3d, &
                                  register_cmip_diag_field_2d, &
                                  send_cmip_data_3d, &
                                  cmip_diag_id_type, &
                                  query_cmip_diag_id

use atmos_dust_mod, only: n_dust_tracers, dust_tracers

implicit none

private

!-----------------------------------------------------------------------
!     ... interfaces
!-----------------------------------------------------------------------
public  tropchem_driver, tropchem_driver_init,  &
        tropchem_driver_time_vary, tropchem_driver_endts

!-----------------------------------------------------------------------
!     ...  declare type that will store the field infomation for the
!          emission file
!-----------------------------------------------------------------------
type,public :: field_init_type
   character(len=64), pointer :: field_names(:)
   real, pointer :: scale_emis(:)
end type field_init_type


!-----------------------------------------------------------------------
!     ... namelist
!-----------------------------------------------------------------------
integer, parameter :: maxinv = 100
real               :: relaxed_dt = SECONDS_PER_DAY*10.,     & ! relaxation timescale (sec) for the upper boundary values
                      relaxed_dt_lbc = SECONDS_PER_DAY*10., & ! relaxation timescale (sec) for the lower boundary values
                      ub_pres = 100.e2,               & ! pressure (Pa) above which to apply chemical upper boundary conditions
                      lb_pres = 950.e2                  ! pressure (Pa) below which to apply chemical lower boundary conditions
character(len=64)  :: file_sulfate = 'sulfate.nc',    & ! NetCDF file for sulfate concentrations
                      file_conc = 'conc_all.nc',      & ! NetCDF file for tracer concentrations (initial and fixed)
                      file_emis_1 = 'emissions.',     & ! NetCDF file name (beginning) for emissions
                      file_emis_2 = '.nc',            & ! NetCDF file name (end) for emissions
                      file_emis3d_1 = 'emissions3d.', & ! NetCDF file name (beginning) for 3-D emissions
                      file_emis3d_2 = '.nc',          & ! NetCDF file name (end) for 3-D emissions
                      file_ub = 'ub_vals.nc'            ! NetCDF file for chemical upper boundary conditions
character(len=64)  :: file_aircraft = 'aircraft.nc',  & ! NetCDF file for aircraft emissions
                      file_jval_lut = 'jvals.v5',     & ! ascii file for photolysis rate lookup table
                      file_jval_lut_min = ''            ! ascii file for photolysis rate LUT (for solar min)
character(len=10), dimension(maxinv) :: inv_list =''    ! list of invariant (fixed) tracers
real               :: aircraft_scale_factor = -999.     ! aircraft emissions scale factor
real               :: lght_no_prd_factor = 1.           ! lightning NOx scale factor
logical            :: normalize_lght_no_prd_area = .false. ! normalize lightning NOx production by grid cell area
logical            :: allow_small_storms_lght_no_prd = .false. ! modify area normalization for high-res grids
real               :: min_land_frac_lght = -999.        ! minimum land fraction for lightning NOx calculation
real               :: strat_chem_age_factor = 1.        ! scale factor for age of air
real               :: strat_chem_dclydt_factor = 1.     ! scale factor for dcly/dt
logical            :: do_tropchem = .false.             ! Do tropospheric chemistry?
logical            :: use_tdep_jvals = .false.          ! Use explicit temperature dependence for photolysis rates
real               :: o3_column_top = 10.               ! O3 column above model top (DU)
real               :: jno_scale_factor = 1.             ! scale factor for NO photolysis rate (jNO)
logical            :: repartition_water_tracers = .false. ! Allow PSC scheme to act on total water (vapor+condensed)
logical            :: allow_negative_cosz = .false.     ! Allow negative values for cosine of solar zenith angle
logical            :: allow_psc_settling_type1 = .false.! Allow Type-I (NAT) PSCs to settle
logical            :: allow_psc_settling_type2 = .false.! Allow Type-II (ice) PSCs to settle
logical            :: force_cly_conservation = .false.  ! Force chemical conservation of Cly
logical            :: rescale_cly_components = .false.  ! Rescale individual Cly components to total Cly VMR
logical            :: set_min_h2o_strat = .false.       ! Don't allow total water concentration in the stratosphere to fall below 2*CH4_trop
character(len=64)  :: ch4_filename = 'ch4_gblannualdata'! Methane timeseries filename
real               :: ch4_scale_factor = 1.             ! Methane scale factor to convert to VMR (mol/mol)
character(len=64)  :: cfc_lbc_filename = 'chemlbf'      ! Input file for CFC lower boundary conditions
logical            :: time_varying_cfc_lbc = .true.     ! Allow time variation of CFC lower boundary conditions
integer, dimension(6) :: cfc_lbc_dataset_entry = (/ 1, 1, 1, 0, 0, 0 /) ! Entry date for CFC lower boundary condition file
integer            :: verbose = 3                       ! level of diagnostic output
logical            :: do_fastjx_photo = .false.         ! use fastjx routine ?
character(len=32)  :: clouds_in_fastjx = 'lsc_only'     ! nature of clouds seen in fastjx calculation; may currently be 'none' or 'lsc_only' (default)
logical            :: check_convergence = .false.       ! if T, non-converged chem tendencies will not be used
real               :: e90_tropopause_vmr = 9.e-8        ! e90 tropopause concentration
logical            :: time_varying_solarflux = .false.  ! allow solar cycle on fastjx v7.1

!co2
real*8             :: co2_fixed_value   = 330e-6
!character(len=64)  :: co2_filename = 'co2_gblannualdata'
character(len=64)  :: co2_filename = 'no_file'
real               :: co2_scale_factor = 1.e-6
real               :: co2_fixed_year   = -999

character(len=64)  :: cloud_chem_pH_solver     = 'bisection'
character(len=64)  :: cloud_chem_type          = 'f1p'
logical            :: het_chem_fine_aerosol_only = .false.
real               :: min_lwc_for_cloud_chem     = 1.e-8
real               :: frac_dust_incloud          = 0
real               :: frac_aerosol_incloud       = 1
real               :: cloud_pH                   = -999  !<0 do not force
real               :: max_rh_aerosol             = 9999      !max rh used for aerosol thermo (to make sure no filter)
logical            :: limit_no3                  = .true.   !for isorropia/stratosphere
logical            :: scale_dust_uptake      = .true.  !adjust reaction rate by available alkalinity

character(len=64)  :: aerosol_thermo_method = 'legacy'               ! other choice isorropia
character(len=64)  :: het_chem_type         = 'legacy'
real               :: gN2O5                 = 0.1
real               :: gNO2                  = 1e-4
real               :: gSO2                  = 0.
real               :: gSO2_dust             = 0.
real               :: gNH3                  = 0.05
real               :: gHNO3_dust            = 0.
real               :: gNO3_dust             = -999.
real               :: gN2O5_dust            = -999.
real               :: gH2SO4_dust           = 0.
logical            :: do_h2so4_nucleation   = .false.
logical            :: cloud_ho2_h2o2        = .true.
real               :: gNO3                  = 0.1
real               :: gHO2                  = 1.

real               :: rh_het_max            = 9999. !maximum rh used to calculate surface area.

character(len=128) :: sim_data_filename = 'sim.dat'      ! Input file for chemistry pre-processor

character(len=64)  :: gso2_dynamic          = 'none'
character(len=64)  :: ghno3_dust_dynamic    = 'none'

real               :: min_t_sfc_cld_chem    = -999   !by default, this filter is turned off. A more reasonable choice would be 273.15K

logical            :: modulate_frac_ic = .false. !modulate the fraction of acids and aerosol in clouds based on frac_liq

real               :: NO2_SO2_max = -999 !NO2 max concentration below which gso2 is scaled by NO2/NO2_SO2_max

integer, parameter :: max_scale_emis_fields = 10
real               :: scale_emis_field_values(max_scale_emis_fields)
character(len=64)  :: scale_emis_field_names(max_scale_emis_fields)

logical            :: do_terpene_emis_bug   = .false. !error double counting biogenic terpene emissions


type(tropchem_diag),  save :: trop_diag
type(tropchem_opt),   save :: trop_option
type (domain2D), pointer :: tropchem_domain !< Atmosphere domain




integer :: n_hno3d, n_so4d

namelist /tropchem_driver_nml/    &
                               relaxed_dt, &
                               relaxed_dt_lbc, &
                               ub_pres, &
                               lb_pres, &
                               file_sulfate, &
                               file_conc, &
                               file_emis_1, &
                               file_emis_2, &
                               file_emis3d_1, &
                               file_emis3d_2, &
                               file_ub, &
                               inv_list, &
                               file_aircraft,&
                               aircraft_scale_factor, &
                               lght_no_prd_factor, &
                               normalize_lght_no_prd_area, &
                               allow_small_storms_lght_no_prd, &
                               min_land_frac_lght, &
                               strat_chem_age_factor, &
                               strat_chem_dclydt_factor, &
                               do_tropchem, &
                               use_tdep_jvals, &
                               file_jval_lut, &
                               file_jval_lut_min, &
                               o3_column_top, &
                               jno_scale_factor, &
                               repartition_water_tracers, &
                               allow_negative_cosz, &
                               allow_psc_settling_type1, &
                               allow_psc_settling_type2, &
                               force_cly_conservation, &
                               rescale_cly_components, &
                               set_min_h2o_strat, &
                               ch4_filename, &
                               ch4_scale_factor, &
                               co2_fixed_value, &
                               co2_fixed_year, &
                               co2_filename, &
                               co2_scale_factor, &
                               cfc_lbc_filename, &
                               time_varying_cfc_lbc, &
                               cfc_lbc_dataset_entry, &
                               verbose,   &
                               do_fastjx_photo, &
                               clouds_in_fastjx, &
                               check_convergence, &
                               e90_tropopause_vmr, &
                               aerosol_thermo_method, &
                               het_chem_type, &
                               gn2o5,gno2,gno3,gso2,gnh3,ghno3_dust,gh2so4_dust,gho2,ghno3_dust_dynamic,gso2_dust,gn2o5_dust,gno3_dust, &
                               do_h2so4_nucleation, &
                               check_convergence, &
                               cloud_chem_pH_solver, &
                               cloud_chem_type, &
                               min_lwc_for_cloud_chem, &
                               het_chem_fine_aerosol_only, &
                               cloud_pH, &
                               frac_dust_incloud, frac_aerosol_incloud, &
                               max_rh_aerosol, limit_no3, cloud_ho2_h2o2, &
                               sim_data_filename,time_varying_solarflux, gso2_dynamic, &
                               rh_het_max, &
                               modulate_frac_ic, &
                               min_t_sfc_cld_chem, &
                               NO2_SO2_max, &
                               scale_emis_field_names, scale_emis_field_values, &
                               do_terpene_emis_bug


integer                     :: nco2 = 0
character(len=7), parameter :: module_name = 'tracers'
real, parameter :: g_to_kg    = 1.e-3,    & !conversion factor (kg/g)
                   m2_to_cm2  = 1.e4,     & !conversion factor (cm2/m2)
                   twopi      = 2.*PI

real, parameter :: mw_so4     = 96e-3 !kg/mol

real, parameter :: emis_cons = WTMAIR * g_to_kg * m2_to_cm2 / AVOGNO
logical, dimension(pcnstm1) :: has_emis = .false., &      ! does tracer have surface emissions?
                               has_emis3d = .false., &    ! does tracer have 3-D emissions?
                               land_does_emission = .false., &    ! surface emission in land
                               has_xactive_emis = .false., & ! does tracer have interactive emissions?
                               diurnal_emis = .false., &   ! diurnally varying emissions?
                               diurnal_emis3d = .false.    ! diurnally varying 3-D emissions?

type(interpolate_type),dimension(pcnstm1), save :: inter_emis, &
                                                   inter_emis3d, &
                                                   inter_aircraft_emis
type(interpolate_type), save :: airc_default
type(field_init_type),dimension(pcnstm1) :: emis_field_names, &
                                            emis3d_field_names
logical, dimension(pcnstm1) :: has_ubc = .false., &
     has_lbc    = .false., &
     has_lbc_2d = .false., &
     fixed_lbc_time = .false.
type(time_type), dimension(pcnstm1) :: lbc_entry
logical, dimension(pcnstm1) :: has_airc = .false., lbc_dry=.false.
character(len=64),dimension(pcnstm1) :: ub_names, airc_names,lbc_names
real, dimension(pcnstm1) :: lbc_factor
real, parameter :: small = 1.e-50
integer :: sphum_ndx=0, cl_ndx=0, clo_ndx=0, hcl_ndx=0, hocl_ndx=0, clono2_ndx=0, &
           cl2o2_ndx=0, cl2_ndx=0, clno2_ndx=0, br_ndx=0, bro_ndx=0, hbr_ndx=0, &
           hobr_ndx=0, brono2_ndx=0, brcl_ndx=0, &
           hno3_ndx=0, o3_ndx=0, &
           no_ndx=0, no2_ndx=0, no3_ndx=0, n_ndx=0, n2o5_ndx=0, ho2no2_ndx=0, &
           pan_ndx=0, onit_ndx=0, mpan_ndx=0, isopno3_ndx=0, onitr_ndx=0, &
           extinct_ndx=0, noy_ndx=0, cly_ndx=0, bry_ndx=0, ch4_ndx=0, &
           dms_ndx=0, so4_ndx(6)=0, co_ndx=0, n2o_ndx=0

integer :: o3s_ndx=0
integer :: o3s_e90_ndx=0
integer :: e90_ndx=0

logical :: do_interactive_h2o = .false.         ! Include chemical sources/sinks of water vapor?
real, parameter :: solarflux_min = 1.09082, &   ! solar minimum flux (band 18) [W/m2]
                   solarflux_max = 1.14694      ! solar maximum flux (band 18) [W/m2]
integer :: num_solar_bands = 0

!-----------------------------------------------------------------------
!     ... identification numbers for diagnostic fields
!-----------------------------------------------------------------------
integer :: id_sul, id_temp, id_dclydt, id_dbrydt, id_dclydt_chem, &
           id_psc_sat, id_psc_nat, id_psc_ice, id_volc_aer, &
           id_imp_slv_nonconv, id_srf_o3, id_coszen, id_h2o_chem
integer :: inqa, inql, inqi !index of the three water species(nqa, nql, nqi)
integer :: age_ndx ! index of age tracer
logical :: module_is_initialized=.false.
logical :: use_lsc_in_fastjx

!cmip6 diagnostics
type(cmip_diag_id_type) :: ID_pso4_aq_kg_m2_s, ID_pso4_gas_kg_m2_s, &
                           ID_jno2, ID_jo1d, &
			   ID_lossch4, ID_lossco, ID_lossn2o, ID_o3loss, ID_o3prod
integer :: jno2_ndx, jo1d_ndx

integer, dimension(pcnstm1) :: indices, id_prod, id_loss, id_chem_tend, &
                               id_emis, id_emis3d, id_xactive_emis, &
                               id_ub, id_lb, id_airc
!new diagnostics (f1p)
integer, parameter :: max_dust = 5
integer, dimension(pcnstm1) :: id_prod_mol, id_loss_mol
integer :: id_pso4_h2o2,id_pso4_o3,id_ghno3_d,id_phno3_d(max_dust), id_phno3_g_d, id_pso4_d(max_dust), &
           id_pso4_g_d, id_gso2, id_aerosol_pH, id_cloud_pH, id_cloud_pHw, id_cloud_pHwl, id_cld_amt_chem, id_cld_liq_chem, id_sa_aerosol, id_sa_so4, id_sa_bc, id_sa_oa, id_sa_ss, id_sa_dust

integer :: usr_so2_dust(max_dust), usr_so4_dust(max_dust), usr_hno3_dust(max_dust), usr_no3_dust(max_dust), usr_n2o5_dust(max_dust)
integer :: id_rx_so2_dust(max_dust), id_rx_so4_dust(max_dust), id_rx_hno3_dust(max_dust), id_rx_no3_dust(max_dust), id_rx_n2o5_dust(max_dust)

integer :: id_so2_emis_cmip, id_nh3_emis_cmip
integer :: id_co_emis_cmip, id_no_emis_cmip
integer :: id_co_emis_cmip2, id_no_emis_cmip2
integer :: id_so2_emis_cmip2, id_nh3_emis_cmip2
integer :: id_emico, id_emiso2, id_eminh3
integer :: id_eminox_woL, id_emiisop_woB
logical :: has_ts_avg = .true.   ! currently reading in from monthly mean files.
integer, dimension(phtcnt)  :: id_jval
integer, dimension(gascnt)  :: id_rate_const
integer :: id_prodox, id_lossox ! for production and loss of ox.(jmao,1/1/2011)

type(interpolate_type), save :: conc       ! used to read in the concentration of OH and CH4
type(interpolate_type), save :: sulfate    ! used to read in the data for sulate
type(interpolate_type), save :: ub_default ! used for the upper bound data
type(interpolate_type),dimension(pcnstm1), save :: ub

type :: lb_type
   real, dimension(:), pointer :: gas_value
   type(time_type), dimension(:), pointer :: gas_time
end type lb_type
type(lb_type), dimension(pcnstm1) :: lb
type(interpolate_type),dimension(pcnstm1), save :: lbc_interp

type :: co2_type
   logical                                :: use_fix_value
   real                                   :: fixed_value
   real,dimension(:), pointer             :: gas_value
   type(time_type), dimension(:), pointer :: gas_time
   logical                                :: use_fix_time
   type(time_type)                        :: fixed_entry
end type co2_type
type(co2_type) :: co2_t
integer :: clock_id,ndiag

!++van
real, allocatable    :: nb_N_Ox(:)
!--van

type (horiz_interp_type), save :: Interp


!---- version number ---------------------------------------------------
character(len=128), parameter :: version     = '$Id$'
character(len=128), parameter :: tagname     = '$Name$'
!-----------------------------------------------------------------------

contains


!#######################################################################

! <SUBROUTINE NAME="tropchem_driver">
!   <OVERVIEW>
!     Tropospheric chemistry driver.
!   </OVERVIEW>
!   <DESCRIPTION>
!     This subroutine calculates the sources and sinks of tracers
!     due to tropospheric chemistry. It is called from atmos_tracer_driver.
!   </DESCRIPTION>
!   <TEMPLATE>
!     call tropchem_driver (lon, lat, land, ocn_flx_fraction, pwt, r, chem_dt,           &
!                           Time, phalf, pfull, t, is, ie, js, je, dt, &
!                           z_half, z_full, q, tsurf, albedo, coszen,  &
!                           area, w10m, flux_sw_down_vis_dir, flux_sw_down_vis_dif, &
!                           half_day, &
!                           Time_next, rdiag,  do_nh3_atm_ocean_exchange, kbot)
!   </TEMPLATE>
!   <IN NAME="lon" TYPE="real" DIM="(:,:)">
!     The longitudes for the local domain.
!   </IN>
!   <IN NAME="lat" TYPE="real" DIM="(:,:)">
!     The latitudes for the local domain.
!   </IN>
!   <IN NAME="land" TYPE="real" DIM="(:,:)">
!     Land fraction
!   </IN>
!   <IN NAME="ocn_flx_fraction" TYPE="real" DIM="(:,:)">
!     Fraction of cell through which ocean flux is allowed
!   </IN>
!   <IN NAME="pwt" TYPE="real" DIM="(:,:,:)">
!     Pressure weighting (air mass) for each layer (kg/m2)
!   </IN>
!   <IN NAME="r" TYPE="real" DIM="(:,:,:,:)">
!     Tracer mixing ratios (tropchem tracers in VMR)
!   </IN>
!   <IN NAME="Time, Time_next" TYPE="time_type">
!     Model time
!   </IN>
!   <IN NAME="phalf" TYPE="real" DIM="(:,:,:)">
!     Pressure on the model half levels (Pa)
!   </IN>
!   <IN NAME="pfull" TYPE="real" DIM="(:,:,:)">
!     Pressure on the model full levels (Pa)
!   </IN>
!   <IN NAME="t" TYPE="real" DIM="(:,:,:)">
!     Temperature.
!   </IN>
!   <IN NAME="is, js" TYPE="integer">
!     Local domain start indices
!   </IN>
!   <IN NAME="ie, je" TYPE="integer">
!     Local domain end indices
!   </IN>
!   <IN NAME="dt" TYPE="real">
!     Model physics timestep (s)
!   </IN>
!   <IN NAME="z_half" TYPE="real" DIM="(:,:,:)">
!     Height at model half levels (m)
!   </IN>
!   <IN NAME="z_full" TYPE="real" DIM="(:,:,:)">
!     Height at model full levels (m)
!   </IN>
!   <IN NAME="q" TYPE="real" DIM="(:,:,:)">
!     Specific humidity (kg/kg)
!   </IN>
!   <IN NAME="tsurf" TYPE="real" DIM="(:,:)">
!     Surface temperature (K)
!   </IN>
!   <IN NAME="albedo" TYPE="real" DIM="(:,:)">
!     Surface albedo
!   </IN>
!   <IN NAME="coszen" TYPE="real" DIM="(:,:)">
!     Cosine of the solar zenith angle
!   </IN>
!   <IN NAME="area" TYPE="real" DIM="(:,:)">
!     Grid box area (m^2)
!   </IN>
!   <IN NAME="w10m" TYPE="real" DIM="(:,:)">
!     Windspeed at 10m (m/s)
!   </IN>
!   <IN NAME="half_day" TYPE="real" DIM="(:,:)">
!     Half-day length  (dimensionless; 0 to pi)
!   </IN>
!   <IN NAME="do_nh3_atm_ocean_exchange," TYPE="logical">
!     Allow interactive atm-ocn exchange of NH3?
!   </IN>
!   <OUT NAME="chem_dt" TYPE="real" DIM="(:,:,:,:)">
!     Tracer tendencies from tropospheric chemistry (VMR/s)
!   </OUT>
!   <INOUT NAME="rdiag" TYPE="real" DIM="(:,:,:,:)">
!     Diagnostic tracer mixing ratios (tropchem tracers in VMR),
!     updated on output
!   </INOUT>
!   <IN NAME="kbot" TYPE="integer, optional" DIM="(:,:)">
!     Integer array describing which model layer intercepts the surface.
!   </IN>

subroutine tropchem_driver( lon, lat, land, ocn_flx_fraction, pwt, r, chem_dt, &
                            Time, phalf, pfull, t, is, ie, js, je, dt,         &
                            z_half, z_full, q, tsurf, albedo, coszen, rrsun,   &
                            area, w10m, half_day,                              &
                            Time_next, rdiag, do_nh3_atm_ocean_exchange, kbot )

!-----------------------------------------------------------------------
   real, intent(in),    dimension(:,:)            :: lon, lat
   real, intent(in),    dimension(:,:)            :: land    ! land fraction
   real, intent(in),    dimension(:,:)            :: ocn_flx_fraction ! grid box fraction over which DMS flux from ocean occurs
   real, intent(in),    dimension(:,:,:)          :: pwt
   real, intent(in),    dimension(:,:,:,:)        :: r
   real, intent(out),   dimension(:,:,:,:)        :: chem_dt
   type(time_type), intent(in)                    :: Time, Time_next
   integer, intent(in)                            :: is, ie, js, je
   real, intent(in),    dimension(:,:,:)          :: phalf,pfull,t
   real, intent(in)                               :: dt      ! timestep (s)
   real, intent(in),    dimension(:,:,:)          :: z_half  ! height in meters at half levels
   real, intent(in),    dimension(:,:,:)          :: z_full  ! height in meters at full levels
   real, intent(in),    dimension(:,:,:)          :: q       ! specific humidity at current time step (kg/kg)
   real, intent(in),    dimension(:,:)            :: tsurf   ! surface temperature (K)
   real, intent(in),    dimension(:,:)            :: albedo  ! surface albedo
   real, intent(in),    dimension(:,:)            :: coszen  ! cosine of solar zenith angle
   real, intent(in)                               :: rrsun   ! earth-sun distance factor (r_avg/r)^2
   real, intent(in),    dimension(:,:)            :: area    ! grid box area (m^2)
   real, intent(in),    dimension(:,:)            :: w10m    ! wind speed at 10m (m/s)
   real, intent(in), dimension(:,:)               :: half_day! half-day length (0 to pi)
   real, intent(inout), dimension(:,:,:,:)        :: rdiag   ! diagnostic tracer concentrations
   logical, intent(in)                            :: do_nh3_atm_ocean_exchange
   integer, intent(in),  dimension(:,:), optional :: kbot
!-----------------------------------------------------------------------
   real, dimension(size(r,1),size(r,2),size(r,3)) :: sulfate_data
!  real, dimension(size(r,1),size(r,2),size(r,3)) :: ub_temp,rno
   real, dimension(size(r,1),size(r,2),size(r,3),maxinv) :: inv_data
   real, dimension(size(r,1),size(r,2)) :: emis, r_lb_2d
   real, dimension(size(r,1),size(r,2), pcnstm1) :: emisz
   real, dimension(size(r,1),size(r,2),size(r,3)) :: emis3d, xactive_emis
   real, dimension(size(r,1),size(r,2),size(r,3)) :: age, cly0, cly, cly_ratio, &
                                                     bry, dclydt, dbrydt, noy, &
                                                     extinct, strat_aerosol
   real, dimension(size(r,1),size(r,2),size(r,3),3) :: psc_vmr_save, dpsc_vmr
   real, dimension(size(r,1),size(r,2)) :: tsfcair, pwtsfc, flux_sw_down_vis
   integer :: i,j,k,n,kb,id,jd,kd,ninv,nt,ntp, idx, index1, index2
!  integer :: nno,nno2
   integer :: inv_index
   integer :: plonl
   logical :: used
   real :: scale_factor, frac, ico2
   real,  dimension(size(r,1),size(r,3)) :: pdel, h2so4, h2o_temp, qlocal, cloud_water, co2_2d, frac_liq, s_down
   real, dimension(size(r,1),size(r,2),size(r,3),pcnstm1)  :: r_temp, r_in, emis_source, r_ub, airc_emis
   real, dimension(size(r,1),size(r,2),size(r,3)) :: tend_tmp, extra_h2o
   real, dimension(pcnstm1) :: r_lb
   real, dimension(size(land,1), size(land,2)) :: oro ! 0 and 1 rep. of land
   real, dimension(size(r,1),size(r,2)) :: coszen_local, fracday_local, scale_dry_lbc
   real :: rrsun_local
   real, dimension(size(r,1),size(r,2),size(r,3),pcnstm1) :: prod, loss
   real, dimension(size(r,1),size(r,2),size(r,3)):: prodox, lossox
   real, dimension(size(r,1),size(r,2),size(r,3),phtcnt) :: jvals
   real, dimension(size(r,1),size(r,2),size(r,3),gascnt) :: rate_constants
   real, dimension(size(r,1),size(r,2),size(r,3)) :: imp_slv_nonconv
   real, dimension(size(r,1),size(r,2),size(r,3)):: e90_vmr
   real, dimension(size(r,1),size(r,2),size(r,3)):: dz
   real, dimension(size(r,1),size(r,2),size(r,3)):: temp_so4
   real :: solar_phase
   real :: solflxband(num_solar_bands)
   type(psc_type) :: psc
   !f1p
   !trop diag arrays
   real, dimension(size(r,1),size(r,2),size(r,3),trop_diag%nb_diag) :: trop_diag_array
   type(time_type) :: co2_time
   character(len=32) :: tracer_name, noytracer
   type(time_type) :: lbc_time
!-----------------------------------------------------------------------

!<ERROR MSG="tropchem_driver_init must be called first." STATUS="FATAL">
!   Tropchem_driver_init needs to be called before tracer_driver.
!</ERROR>
   if (.not. module_is_initialized)  &
      call error_mesg ('Tropchem_driver','tropchem_driver_init must be called first.', FATAL)

   ntp = size(r,4)
   plonl = size(r,1)

!initialize diagnostic array
   trop_diag_array(:,:,:,:) = 0.
   where(land(:,:) >= 0.5)
      oro(:,:) = 1.
   elsewhere
      oro(:,:) = 0.
   endwhere

   id=size(r,1); jd=size(r,2); kd=size(r,3)

   ninv=0
   do n = 1, size(inv_list)
      if(inv_list(n) /= '') then
         ninv = ninv + 1
      else
         exit
      end if
   end do

   emis_source(:,:,:,:) = 0.0
   airc_emis(:,:,:,:) = 0.0
   emisz(:,:,:) = 0.0

   tsfcair(:,:) = t(:,:,kd)
   pwtsfc(:,:) = t(:,:,kd)
   dz(:,:,:) = z_half(:,:,1:kd) - z_half(:,:,2:(kd+1))

   do n = 1, pcnstm1
!-----------------------------------------------------------------------
!     ... read in the surface emissions, using interpolator
!-----------------------------------------------------------------------
      if (has_emis(n)) then
         if (do_nh3_atm_ocean_exchange .and. tracnam(n)(1:3).eq.'NH3') then
            call read_2D_emis_data( inter_emis(n), emis, Time, Time_next, &
                 emis_field_names(n)%field_names, &
                 diurnal_emis(n), coszen, half_day, lon, &
                 is, js, has_xactive_emis(n),emis_field_names(n)%scale_emis,'ocean')
         else
            call read_2D_emis_data( inter_emis(n), emis, Time, Time_next, &
                 emis_field_names(n)%field_names, &
                 diurnal_emis(n), coszen, half_day, lon, &
                 is, js, has_xactive_emis(n),emis_field_names(n)%scale_emis)
         end if
         if ( land_does_emission(n) ) then
            emis = emis * ( 1. - land )
         end if

         if (id_emis(n) > 0) then
            used = send_data(id_emis(n),emis,Time_next,is_in=is,js_in=js)
         end if

         emisz(:,:,n) = emis(:,:)

         if (tracnam(n) == 'NO' .and. id_no_emis_cmip > 0) then
           used = send_data(id_no_emis_cmip,emis*1.0e04*0.030/AVOGNO, &
                              Time_next, &
                              is_in=is,js_in=js)
         endif
         if (tracnam(n) == 'CO' .and. id_co_emis_cmip > 0) then
           used = send_data(id_co_emis_cmip,emis*1.0e04*0.028/AVOGNO,Time_next, &
                              is_in=is,js_in=js)
         endif
         if (tracnam(n) == 'SO2' .and. id_so2_emis_cmip > 0) then
           used = send_data(id_so2_emis_cmip,emis*1.0e04*0.064/AVOGNO,Time_next, &
                              is_in=is,js_in=js)
         endif
         if (tracnam(n) == 'NH3' .and. id_nh3_emis_cmip > 0) then
           used = send_data(id_nh3_emis_cmip,emis*1.0e04*0.017/AVOGNO,Time_next, &
                              is_in=is,js_in=js)
         endif

         if (present(kbot)) then
            do j=1,jd
               do i=1,id
                  kb=kbot(i,j)
                  emis_source(i,j,kb,n) = emis(i,j)/pwt(i,j,kb) * emis_cons
               end do
            end do
         else
            emis_source(:,:,kd,n) = emis(:,:)/pwt(:,:,kd) * emis_cons
         end if
      end if

!-----------------------------------------------------------------------
!     ... read in the 3-D emissions, using interpolator
!-----------------------------------------------------------------------
      if (has_emis3d(n)) then
         call read_3D_emis_data( inter_emis3d(n), emis3d, Time, Time_next,phalf, &
                                 emis3d_field_names(n)%field_names, &
                                 diurnal_emis3d(n), coszen, half_day, lon, &
                                 is, js, id_emis3d(n), emis3d_field_names(n)%scale_emis )

         emis_source(:,:,:,n) = emis_source(:,:,:,n) &
                              + emis3d(:,:,:)/pwt(:,:,:) * emis_cons
         do k=1, size(emis3d,3)
           emisz(:,:,n) = emisz(:,:,n) + emis3d(:,:,k)
         end do
      end if

!-----------------------------------------------------------------------
!     ... calculate interactive (DMS only) emissions
!-----------------------------------------------------------------------
      if ( has_xactive_emis(n) .or. id_xactive_emis(n)>0 ) then
         if (trim(tracnam(n)) .eq. "DMS") then
            call calc_xactive_emis( n, Time, Time_next,lon, lat, pwt, is, ie, js, je, &
                 area, land, ocn_flx_fraction,tsurf, w10m, xactive_emis, &
                 kbot=kbot, id_emis_diag=id_xactive_emis(n) )
            if (has_xactive_emis(n)) then
               do k=1, size(emis3d,3)
                  emisz(:,:,n) = emisz(:,:,n) + xactive_emis(:,:,k) * pwt(:,:,k) / emis_cons
               end do
               emis_source(:,:,:,n) = emis_source(:,:,:,n) + xactive_emis(:,:,:)
            end if
         endif
      end if

!-----------------------------------------------------------------------
!     ... read in the aircraft emissions
!-----------------------------------------------------------------------
      if(has_airc(n)) then
         call interpolator( inter_aircraft_emis(n), Time, phalf, &
                            airc_emis(:,:,:,n), trim(airc_names(n)),is,js)
         if (aircraft_scale_factor >= 0.) airc_emis(:,:,:,n) = airc_emis(:,:,:,n) * aircraft_scale_factor
         if(id_airc(n) > 0) &
              used = send_data(id_airc(n),airc_emis(:,:,:,n),Time_next, is_in=is, js_in=js)

         do k=1, size(emis3d,3)
           emisz(:,:,n) = emisz(:,:,n) + airc_emis(:,:,k,n)
         end do
         airc_emis(:,:,:,n) = airc_emis(:,:,:,n)/pwt(:,:,:)*emis_cons
      end if
         if (tracnam(n) == 'NO') then
           if (id_no_emis_cmip2 > 0) then
             used = send_data(id_no_emis_cmip2,emisz(:,:,n)*1.0e04*0.030/AVOGNO,Time_next, &
                                                 is_in=is,js_in=js)
           endif
           if (id_eminox_woL > 0) then
             used = send_data(id_eminox_woL, emisz(:,:,n)*1.0e04*0.014/AVOGNO, Time_next, is_in=is,js_in=js)
           endif
         endif
         if (tracnam(n) == 'CO') then
           if (id_co_emis_cmip2 > 0) then
             used = send_data(id_co_emis_cmip2,emisz(:,:,n)*1.0e04*0.028/AVOGNO,Time_next, &
                                                 is_in=is,js_in=js)
           endif
           if (id_emico > 0) then
             used = send_data(id_emico, emisz(:,:,n)*1.0e04*0.028/AVOGNO, Time_next, is_in=is,js_in=js)
           endif
         endif
         if (tracnam(n) == 'SO2') then
           if (id_so2_emis_cmip2 > 0) then
             used = send_data(id_so2_emis_cmip2,emisz(:,:,n)*1.0e04*0.064/AVOGNO,Time_next, &
                                                 is_in=is,js_in=js)
           endif
           if (id_emiso2 > 0) then
             used = send_data(id_emiso2, emisz(:,:,n)*1.0e04*0.064/AVOGNO, Time_next, is_in=is,js_in=js)
           endif
         endif
         if (tracnam(n) == 'NH3') then
           if (id_nh3_emis_cmip2 > 0) then
             used = send_data(id_nh3_emis_cmip2,emisz(:,:,n)*1.0e04*0.017/AVOGNO,Time_next, &
                                                  is_in=is,js_in=js)
           endif
           if (id_eminh3 > 0) then
             used = send_data(id_eminh3, emisz(:,:,n)*1.0e04*0.017/AVOGNO, Time_next, is_in=is,js_in=js)
           endif
         endif
         if (tracnam(n) == 'ISOP') then
           if (id_emiisop_woB > 0) then
             used = send_data(id_emiisop_woB, emisz(:,:,n)*1.0e04*0.068/AVOGNO, Time_next, is_in=is,js_in=js)
           endif
         endif
   end do

!-----------------------------------------------------------------------
!     ... read in the concentrations of "invariant" (i.e., prescribed)
!         species
!-----------------------------------------------------------------------
   do n = 1,ninv
      call interpolator( conc, Time, phalf, inv_data(:,:,:,n), &
                         trim(inv_list(n)), is, js)
      inv_index = get_tracer_index( MODEL_ATMOS, trim(inv_list(n)) ) - ntp
      rdiag(:,:,:,inv_index) = inv_data(:,:,:,n)
   end do

!-----------------------------------------------------------------------
!     ... read in the sulfate aerosol concentrations
!-----------------------------------------------------------------------
   call interpolator(sulfate, Time, phalf, sulfate_data, 'sulfate', is,js)
   used = send_data(id_sul, sulfate_data, Time_next, is_in=is, js_in=js)

   call mpp_clock_begin(clock_id)

   chem_dt(:,:,:,:) =0.

!-----------------------------------------------------------------------
!     ... assign concentrations of prognostic (r) and diagnostic (rdiag)
!         species to r_temp
!-----------------------------------------------------------------------
   do n = 1,pcnstm1
      if(indices(n) <= ntp) then
         r_temp(:,:,:,n) = r(:,:,:,indices(n))
      else
         r_temp(:,:,:,n) = rdiag(:,:,:,indices(n)-ntp)
      end if
   end do

!-----------------------------------------------------------------------
!     ... convert to H2O VMR
!-----------------------------------------------------------------------
   if (sphum_ndx > 0) then
      r_temp(:,:,:,sphum_ndx) = r_temp(:,:,:,sphum_ndx) * WTMAIR / WTMH2O
   end if

!-----------------------------------------------------------------------
!     ... convert volcanic aerosol extinction into aerosol surface area
!-----------------------------------------------------------------------
   if (extinct_ndx > 0 .and. extinct_ndx <= ntp) then
      extinct(:,:,:) = r(:,:,:,extinct_ndx)
   else if (extinct_ndx > ntp) then
      extinct(:,:,:) = rdiag(:,:,:,extinct_ndx-ntp)
   else
      extinct(:,:,:) = 0.
   end if
   call strat_chem_get_aerosol( extinct, strat_aerosol )

!-----------------------------------------------------------------------
!     ... get age of air
!-----------------------------------------------------------------------
   if(age_ndx > 0 .and. age_ndx <= ntp) then
      age(:,:,:) = r(:,:,:,age_ndx)
   else
      age(:,:,:) = 0.
   end if

!-----------------------------------------------------------------------
!     ... Chemical families
!-----------------------------------------------------------------------
   cly0(:,:,:) = 0.
   if (cl_ndx>0) then
      cly0(:,:,:) = cly0(:,:,:) + r_temp(:,:,:,cl_ndx)
   end if
   if (clo_ndx>0) then
      cly0(:,:,:) = cly0(:,:,:) + r_temp(:,:,:,clo_ndx)
   end if
   if (hcl_ndx>0) then
      cly0(:,:,:) = cly0(:,:,:) + r_temp(:,:,:,hcl_ndx)
   end if
   if (hocl_ndx>0) then
      cly0(:,:,:) = cly0(:,:,:) + r_temp(:,:,:,hocl_ndx)
   end if
   if (clono2_ndx>0) then
      cly0(:,:,:) = cly0(:,:,:) + r_temp(:,:,:,clono2_ndx)
   end if
   if (cl2o2_ndx>0) then
      cly0(:,:,:) = cly0(:,:,:) + r_temp(:,:,:,cl2o2_ndx)*2
   end if
   if (cl2_ndx>0) then
      cly0(:,:,:) = cly0(:,:,:) + r_temp(:,:,:,cl2_ndx)*2
   end if
   if (clno2_ndx>0) then
      cly0(:,:,:) = cly0(:,:,:) + r_temp(:,:,:,clno2_ndx)
   end if
   if (brcl_ndx>0) then
      cly0(:,:,:) = cly0(:,:,:) + r_temp(:,:,:,brcl_ndx)
   end if

!-----------------------------------------------------------------------
!     ... cosine of solar zenith angle
!-----------------------------------------------------------------------
   if (allow_negative_cosz) then
      call diurnal_solar( lat, lon, Time, coszen_local, fracday_local, &
                          rrsun_local, dt_time=real_to_time_type(dt), &
                          allow_negative_cosz=.true. )
   else
      coszen_local(:,:) = coszen(:,:)
      rrsun_local = rrsun
   end if

   r_temp(:,:,:,:) = MAX(r_temp(:,:,:,:),small)

!set CO2
   if (nco2 == NO_TRACER) then
      if (co2_t%use_fix_value) then
         co2_2d(:,:) = co2_t%fixed_value
      else
         if (co2_t%use_fix_time) then
            co2_time = co2_t%fixed_entry
         else
            co2_time = Time
         end if
         call time_interp( co2_time, co2_t%gas_time(:), frac, index1, index2 )
         ico2 = co2_t%gas_value(index1) + &
              frac*( co2_t%gas_value(index2) - co2_t%gas_value(index1) )
         co2_2d(:,:) = ico2
      end if
   end if
   do j = 1,jd
      do k = 1,kd
         pdel(:,k) = phalf(:,j,k+1) - phalf(:,j,k)
      end do
      qlocal(:,:) = q(:,j,:)
      if (nco2>0) then
         co2_2d(:,:) = r(:,j,:,nco2)
      end if

!-----------------------------------------------------------------------
!     ... get stratospheric h2so4
!-----------------------------------------------------------------------
      call strat_chem_get_h2so4( pfull(:,j,:), age(:,j,:), h2so4 )

!-----------------------------------------------------------------------
!     ... compute PSC amounts
!-----------------------------------------------------------------------
      if (sphum_ndx>0) then
         h2o_temp(:,:) = r_temp(:,j,:,sphum_ndx)
      else
         h2o_temp(:,:) = qlocal(:,:) * WTMAIR/WTMH2O
      end if
      cloud_water(:,:) = MAX(r(:,j,:,inql)+r(:,j,:,inqi),0.)
      where ( cloud_water(:,:) .gt. 0 )
          frac_liq(:,:) =  MAX( r(:,j,:,inql) , 0.)/cloud_water
      elsewhere
      frac_liq(:,:) = 1.
      endwhere

      !only do cld chemistry at the surface if t is greater than freezing temperature. This is because water saturation vapor pressure over ice is much higher than the saturation vapor pressure of water over liquid water (f1p,sm)
      where (tsurf(:,j) .lt.min_t_sfc_cld_chem)
         frac_liq(:,size(r,3)) = 0.
      endwhere


      if (repartition_water_tracers) then
         h2o_temp(:,:) = h2o_temp(:,:) + cloud_water(:,:) * WTMAIR/WTMH2O
      end if
      if (set_min_h2o_strat) then
         call strat_chem_get_extra_h2o( h2o_temp, age(:,j,:), r_temp(:,j,:,ch4_ndx), Time, extra_h2o(:,j,:) )
         h2o_temp(:,:) = h2o_temp(:,:) + extra_h2o(:,j,:)
      end if

      call strat_chem_get_psc( t(:,j,:), pfull(:,j,:), &
                               r_temp(:,j,:,hno3_ndx), h2o_temp(:,:), &
                               h2so4, strat_aerosol(:,j,:), psc, psc_vmr_out=psc_vmr_save(:,j,:,:) )

      if (repartition_water_tracers) then
         cloud_water(:,:) = MAX(0.,cloud_water(:,:) - psc_vmr_save(:,j,:,3)*WTMH2O/WTMAIR) ! reduce cloud_water by amount of type-II PSC
         h2o_temp(:,:) = h2o_temp(:,:) - cloud_water(:,:) * WTMAIR/WTMH2O                  ! remaining water is present as vapor
      end if
      if (sphum_ndx>0) then
         r_temp(:,j,:,sphum_ndx) = h2o_temp(:,:)
      end if
      qlocal(:,:) = h2o_temp(:,:) * WTMH2O/WTMAIR
      r_in(:,j,:,:) = r_temp(:,j,:,:)

!NEED TO BE MERGED
!!!!!!!
!!!!!!!
!      where( cloud_liq(:,:)+cloud_ice(:,:) .gt. small )
!         s_down(:,:)    = cloud_water(:,:)/(cloud_liq(:,:)+cloud_ice(:,:))
!         cloud_liq(:,:) = cloud_liq(:,:)*s_down(:,:)
!         cloud_ice(:,:) = cloud_ice(:,:)*s_down(:,:)
!      end where
!!!!!!
!!!!!!

!-----------------------------------------------------------------------
!     ... get solar cycle phase (use radiation band #18)
!-----------------------------------------------------------------------
      call get_solar_flux_by_band(solflxband)
      solar_phase = solflxband(num_solar_bands)
      solar_phase = (solar_phase-solarflux_min)/(solarflux_max-solarflux_min)

!-----------------------------------------------------------------------
!     ... get e90 concentrations
!-----------------------------------------------------------------------

   e90_ndx = get_tracer_index( MODEL_ATMOS,'e90' )
   if (e90_ndx > 0) then
      e90_vmr(:,j,:) = r(:,j,:,e90_ndx)
   else
      e90_vmr(:,j,:) = 0.
   end if

!-----------------------------------------------------------------------
!     ... call chemistry driver
!-----------------------------------------------------------------------
      call chemdr(r_temp(:,j,:,:),             & ! species volume mixing ratios (VMR)
                  r(:,j,:,:),                  &
                  phalf(:,j,:),                & ! pressure at boundaries (Pa)
                  pwt(:,j,:) ,                 & ! column air density (Kg/m2)
                  j,                           & ! j
                  Time_next,                   & ! time
                  lat(:,j),                    & ! latitude
                  lon(:,j),                    & ! longitude
                  dt,                          & ! timestep in seconds
                  phalf(:,j,SIZE(phalf,3)),    & ! surface press ( pascals )
                  phalf(:,j,1),                & ! model top pressure (pascals)
                  pfull(:,j,:),                & ! midpoint press ( pascals )
                  pdel,                        & ! delta press across midpoints
                  z_full(:,j,:),               & ! height at midpoints ( m )
                  z_half(:,j,:),               & ! height at interfaces ( m )
                  MAX(r(:,j,:,inqa),0.),       & ! cloud fraction
                  cloud_water(:,:),            & ! total cloud water (kg/kg)
                  frac_liq(:,:),               & ! fraction of liquid water
                  t(:,j,:),                    & ! temperature
                  inv_data(:,j,:,:),           & ! invariant species
                  qlocal(:,:),                 & ! specific humidity ( kg/kg )
                  albedo(:,j),                 & ! surface albedo
                  coszen_local(:,j),           & ! cosine of solar zenith angle
                  rrsun_local,                 & ! earth-sun distance factor
                  prod(:,j,:,:),               & ! chemical production rate
                  loss(:,j,:,:),               & ! chemical loss rate
                  jvals(:,j,:,:),              & ! photolysis rates (s^-1)
                  rate_constants(:,j,:,:),     & ! kinetic rxn rate constants (cm^3 molec^-1 s^-1 for 2nd order)
                  sulfate_data(:,j,:),         & ! sulfate aerosol
                  psc,                         & ! polar stratospheric clouds (PSCs)
                  do_interactive_h2o,          & ! include h2o sources/sinks?
                  solar_phase,                 & ! solar cycle phase (1=max, 0=min)
                  imp_slv_nonconv(:,j,:),      & ! flag for non-convergence of implicit solver
                  plonl,                       & ! number of longitudes
                  prodox(:,j,:),               & ! production of ox(jmao,1/1/2011)
                  lossox(:,j,:),               & ! loss of ox(jmao,1/1/2011)
                  e90_vmr(:,j,:),              & ! e90 concentrations
                  e90_tropopause_vmr,          & ! e90 tropopause threshold
                  co2_2d,                      &
                  trop_diag_array(:,j,:,:),    &
                  trop_option,                 &
                  trop_diag)

      call strat_chem_destroy_psc( psc )

   end do

   r_temp(:,:,:,:) = MAX( r_temp(:,:,:,:), small )
   if (allow_psc_settling_type1 .or. allow_psc_settling_type2) then
      call strat_chem_psc_sediment( psc_vmr_save, pfull, dt, dpsc_vmr )
      if (.not. allow_psc_settling_type1) dpsc_vmr(:,:,:,2) = 0.
      if (.not. allow_psc_settling_type2) dpsc_vmr(:,:,:,3) = 0.
   end if

!-----------------------------------------------------------------------
!     ... output diagnostics
!-----------------------------------------------------------------------
   do n = 1,pcnstm1
      if(id_prod(n)>0) then
         used = send_data(id_prod(n),prod(:,:,:,n),Time_next,is_in=is,js_in=js)
      end if
      if(id_loss(n)>0) then
         used = send_data(id_loss(n),loss(:,:,:,n),Time_next,is_in=is,js_in=js)
      end if
      if(id_prod_mol(n)>0) then
         used = send_data(id_prod_mol(n),prod(:,:,:,n)*pwt(:,:,:)*1.e3/WTMAIR, &
              Time_next,is_in=is,js_in=js)
      end if
      if(id_loss_mol(n)>0) then
         used = send_data(id_loss_mol(n),loss(:,:,:,n)*pwt(:,:,:)*1.e3/WTMAIR, &
              Time_next,is_in=is,js_in=js)
      end if

      if (n == sphum_ndx) then
         scale_factor = WTMAIR/WTMH2O
! add PSC ice back to H2O
!         if (repartition_water_tracers) then
!           r_temp(:,:,:,n) = r_temp(:,:,:,n) + &
!                             MAX( 0.,psc_vmr_save(:,:,:,3) - MAX(r(:,:,:,inql)+r(:,:,:,inqi),0.)*scale_factor ) + &
!                             dpsc_vmr(:,:,:,3)
!        else
!           r_temp(:,:,:,n) = r_temp(:,:,:,n) + psc_vmr_save(:,:,:,3)+dpsc_vmr(:,:,:,3)
!        end if
!        if (set_min_h2o_strat) then
!           r_temp(:,:,:,n) = MAX( r_temp(:,:,:,n) - extra_h2o(:,:,:), small )
!        end if
      else
         scale_factor = 1.
      end if

!     if (n == hno3_ndx) then
!        r_temp(:,:,:,n) = r_temp(:,:,:,n) + psc_vmr_save(:,:,:,2)+dpsc_vmr(:,:,:,2) ! add PSC NAT back to gas-phase HNO3
!     end if

!-----------------------------------------------------------------------
!     ... compute tendency
!-----------------------------------------------------------------------
      tend_tmp(:,:,:) = ( r_temp(:,:,:,n) - r_in(:,:,:,n) )/dt
      if(indices(n) <= ntp) then
!-----------------------------------------------------------------------
!     ... prognostic species
!-----------------------------------------------------------------------
!        tend_tmp(:,:,:) = ( r_temp(:,:,:,n) - MAX(r(:,:,:,indices(n))*scale_factor,small) )/dt
         chem_dt(:,:,:,indices(n)) = airc_emis(:,:,:,n) + emis_source(:,:,:,n) + tend_tmp(:,:,:)
      else
!-----------------------------------------------------------------------
!     ... diagnostic species
!-----------------------------------------------------------------------
!        tend_tmp(:,:,:) = ( r_temp(:,:,:,n) - MAX(rdiag(:,:,:,indices(n)-ntp)*scale_factor,small) )/dt
         rdiag(:,:,:,indices(n)-ntp) = r_temp(:,:,:,n)
      end if
!-----------------------------------------------------------------------
!     ... output diagnostic tendency
!-----------------------------------------------------------------------
      if(id_chem_tend(n)>0) then
         used = send_data( id_chem_tend(n), tend_tmp(:,:,:), Time_next, is_in=is,js_in=js)
      end if

!-----------------------------------------------------------------------
!     ... apply upper boundary condition
!-----------------------------------------------------------------------
      if(has_ubc(n)) then
         call interpolator(ub(n), Time, phalf, r_ub(:,:,:,n), trim(ub_names(n)), is, js)
         if(id_ub(n)>0) then
            used = send_data(id_ub(n), r_ub(:,:,:,n), Time_next, is_in=is, js_in=js)
         end if
         where (pfull(:,:,:) < ub_pres)
            chem_dt(:,:,:,indices(n)) = (r_ub(:,:,:,n) - r(:,:,:,indices(n))) / relaxed_dt
         endwhere
      end if

!-----------------------------------------------------------------------
!     ... apply lower boundary condition
!-----------------------------------------------------------------------
      if(has_lbc(n)) then
         if (.not.has_lbc_2d(n)) then
            if (fixed_lbc_time(n)) then
               lbc_Time = lbc_entry(n)
            else
               lbc_Time = Time
            end if
            call time_interp( lbc_Time, lb(n)%gas_time(:), frac, index1, index2 )
            r_lb(n) = lb(n)%gas_value(index1) + frac*( lb(n)%gas_value(index2) - lb(n)%gas_value(index1))
            if(id_lb(n)>0) then
               used = send_data(id_lb(n), r_lb(n), Time_next)
            end if
            scale_dry_lbc = 1.            
            do k=1,size(chem_dt,3)
               if (lbc_dry(n)) then
                  scale_dry_lbc = (1.-r_temp(:,:,k,sphum_ndx))
               end if
               !r_temp is in vmr
               where (pfull(:,:,k) > lb_pres)
                  chem_dt(:,:,k,indices(n)) = (r_lb(n)*scale_dry_lbc - r(:,:,k,indices(n))) / relaxed_dt_lbc
               end where
            end do                        
         else
            if (fixed_lbc_time(n)) then
               lbc_Time = lbc_entry(n)
            else
               lbc_Time = Time
            end if
                        
            call interpolator(lbc_interp(n), lbc_time, r_lb_2d, trim(lbc_names(n)), is,js)
            r_lb_2d = r_lb_2d*lbc_factor(n)
            if (id_lb(n)>0) then
               used = send_data(id_lb(n), r_lb_2d, Time_next, is_in=is, js_in=js)   
            end if
            scale_dry_lbc = 1.
            do k=1,size(chem_dt,3)
               if (lbc_dry(n)) then
                  scale_dry_lbc = (1.-r_temp(:,:,k,sphum_ndx))
               end if
               !r_temp is in vmr
               where (pfull(:,:,k) > lb_pres)
                  chem_dt(:,:,k,indices(n)) = (r_lb_2d*scale_dry_lbc - r(:,:,k,indices(n))) / relaxed_dt_lbc
               end where
            end do
            
         end if
         
      end if

   end do

!-----------------------------------------------------------------------
!     ... ox budget diagnostics
!-----------------------------------------------------------------------
   if(id_prodox>0) then
      used = send_data(id_prodox, prodox(:,:,:), Time_next, is_in=is, js_in=js)
   end if
   if(id_lossox>0) then
      used = send_data(id_lossox, lossox(:,:,:), Time_next, is_in=is, js_in=js)
   end if

!-----------------------------------------------------------------------
!     ... production/loss diagnostics
!-----------------------------------------------------------------------
   if (query_cmip_diag_id(ID_lossch4) .and. ch4_ndx>0) then
      used = send_cmip_data_3d (ID_lossch4, loss(:,:,:,ch4_ndx)*pwt(:,:,:)*1.e3/(WTMAIR*dz(:,:,:)), &
           Time_next, is_in=is, js_in=js, ks_in=1)
   end if
   if (query_cmip_diag_id(ID_lossco) .and. co_ndx>0) then
      used = send_cmip_data_3d (ID_lossco, loss(:,:,:,co_ndx)*pwt(:,:,:)*1.e3/(WTMAIR*dz(:,:,:)), &
           Time_next, is_in=is, js_in=js, ks_in=1)
   end if
   if (query_cmip_diag_id(ID_lossn2o) .and. n2o_ndx>0) then
      used = send_cmip_data_3d (ID_lossn2o, loss(:,:,:,n2o_ndx)*pwt(:,:,:)*1.e3/(WTMAIR*dz(:,:,:)), &
           Time_next, is_in=is, js_in=js, ks_in=1)
   end if
   if (query_cmip_diag_id(ID_o3loss)) then
      used = send_cmip_data_3d (ID_o3loss, lossox(:,:,:)*pwt(:,:,:)*1.e3/(WTMAIR*dz(:,:,:)), &
           Time_next, is_in=is, js_in=js, ks_in=1)
   end if
   if (query_cmip_diag_id(ID_o3prod)) then
      used = send_cmip_data_3d (ID_o3prod, prodox(:,:,:)*pwt(:,:,:)*1.e3/(WTMAIR*dz(:,:,:)), &
           Time_next, is_in=is, js_in=js, ks_in=1)
   end if

!-----------------------------------------------------------------------
!     ... surface concentration diagnostics
!-----------------------------------------------------------------------
      if ( o3_ndx>0 ) then
         used = send_data(id_srf_o3, r_temp(:,:,size(r_temp,3),o3_ndx), Time_next, is_in=is, js_in=js)
      end if

!-----------------------------------------------------------------------
!     ... aerosol chemistry diagnostics
!-----------------------------------------------------------------------
   if (id_pso4_h2o2>0 .and. trop_diag%ind_pso4_h2o2>0) then
      used = send_data(id_pso4_h2o2,trop_diag_array(:,:,:,trop_diag%ind_pso4_h2o2)*pwt(:,:,:)/(WTMAIR*1e-3), &
           Time_next,is_in=is,js_in=js)
   end if
   if (id_pso4_o3>0 .and. trop_diag%ind_pso4_o3>0) then
      used = send_data(id_pso4_o3,trop_diag_array(:,:,:,trop_diag%ind_pso4_o3)*pwt(:,:,:)/(WTMAIR*1e-3), &
           Time_next,is_in=is,js_in=js)
   end if

   if (query_cmip_diag_id(ID_pso4_gas_kg_m2_s)) then
      temp_so4 = 0.
      do n=1,size(so4_ndx)
         if (so4_ndx(n)>0) then
            temp_so4 = temp_so4 + mw_so4 * prod(:,:,:,so4_ndx(n))*pwt(:,:,:)/(WTMAIR*1e-3)
         end if
      end do
      used = send_cmip_data_3d (ID_pso4_gas_kg_m2_s,  &
           temp_so4, &
           Time_next, is_in=is, js_in=js, ks_in=1)
   end if
   if (query_cmip_diag_id(ID_pso4_aq_kg_m2_s)) then
      used = send_cmip_data_3d (ID_pso4_aq_kg_m2_s,  &
           mw_so4 * (trop_diag_array(:,:,:,trop_diag%ind_pso4_o3)+ trop_diag_array(:,:,:,trop_diag%ind_pso4_h2o2))*pwt(:,:,:)/(WTMAIR*1e-3), &
           Time_next, is_in=is, js_in=js, ks_in=1)
   end if
   do n=1,5
      if (id_phno3_d(n) >0) then
         !phno3_d is in VMR/s convert to mole/m2/s
         used = send_data(id_phno3_d(n),trop_diag_array(:,:,:,trop_diag%ind_phno3_d(n))*pwt(:,:,:)*1.e3/WTMAIR,Time_next,is_in=is,js_in=js)
      end if
   end do
   do n=1,5
      if (id_pso4_d(n) >0) then
         !phno3_d is in VMR/s convert to mole/m2/s
         used = send_data(id_pso4_d(n),trop_diag_array(:,:,:,trop_diag%ind_pso4_d(n))*pwt(:,:,:)*1.e3/WTMAIR,Time_next,is_in=is,js_in=js)
      end if
   end do
   if (id_ghno3_d>0) then
      used = send_data(id_ghno3_d,trop_diag_array(:,:,:,trop_diag%ind_ghno3_d),Time_next,is_in=is,js_in=js)
   end if
   if (id_gso2>0) then
      used = send_data(id_gso2,trop_diag_array(:,:,:,trop_diag%ind_gso2),Time_next,is_in=is,js_in=js)
   end if

   if (id_sa_aerosol>0) then
      used = send_data(id_sa_aerosol,trop_diag_array(:,:,:,trop_diag%ind_sa_aerosol),Time_next,is_in=is,js_in=js)
   end if
   if (id_sa_so4>0) then
      used = send_data(id_sa_so4,trop_diag_array(:,:,:,trop_diag%ind_sa_so4),Time_next,is_in=is,js_in=js)
   end if
   if (id_sa_bc>0) then
      used = send_data(id_sa_bc,trop_diag_array(:,:,:,trop_diag%ind_sa_bc),Time_next,is_in=is,js_in=js)
   end if
   if (id_sa_oa>0) then
      used = send_data(id_sa_oa,trop_diag_array(:,:,:,trop_diag%ind_sa_oa),Time_next,is_in=is,js_in=js)
   end if
   if (id_sa_dust>0) then
      used = send_data(id_sa_dust,trop_diag_array(:,:,:,trop_diag%ind_sa_dust),Time_next,is_in=is,js_in=js)
   end if
   if (id_sa_ss>0) then
      used = send_data(id_sa_ss,trop_diag_array(:,:,:,trop_diag%ind_sa_ss),Time_next,is_in=is,js_in=js)
   end if


   if (id_aerosol_pH>0) then
      used = send_data(id_aerosol_pH,trop_diag_array(:,:,:,trop_diag%ind_aerosol_pH),Time_next,is_in=is,js_in=js, mask = & 
           ( trop_diag_array(:,:,:,trop_diag%ind_aerosol_pH) .gt. (missing_value + tiny(missing_value))))
   end if
   if (id_cloud_pH>0) then
      used = send_data(id_cloud_pH,trop_diag_array(:,:,:,trop_diag%ind_cloud_pH),Time_next,is_in=is,js_in=js, mask = &
           (trop_diag_array(:,:,:,trop_diag%ind_cloud_pH).gt. (missing_value + tiny(missing_value))))
   end if
   if (id_cloud_pHw>0) then
      used = send_data(id_cloud_pHw,trop_diag_array(:,:,:,trop_diag%ind_cloud_pH)*max(r(:,:,:,inqa),0.),Time_next,is_in=is,js_in=js, &
                       mask = (trop_diag_array(:,:,:,trop_diag%ind_cloud_pH).gt. (missing_value + tiny(missing_value))))
   end if
   if (id_cld_amt_chem>0) then
      used = send_data(id_cld_amt_chem,max(r(:,:,:,inqa),0.),Time_next,is_in=is,js_in=js, &
                       mask = (trop_diag_array(:,:,:,trop_diag%ind_cloud_pH).gt. (missing_value + tiny(missing_value))))
   end if

   if (id_cld_liq_chem>0) then
      used = send_data(id_cld_liq_chem,max(r(:,:,:,inql),0.),Time_next,is_in=is,js_in=js, &
                       mask = (trop_diag_array(:,:,:,trop_diag%ind_cloud_pH).gt. (missing_value + tiny(missing_value))))
   end if

   if (id_cloud_pHwl>0) then
      used = send_data(id_cloud_pHwl,trop_diag_array(:,:,:,trop_diag%ind_cloud_pH)*max(r(:,:,:,inql),0.),Time_next,is_in=is,js_in=js, &
                       mask = (trop_diag_array(:,:,:,trop_diag%ind_cloud_pH).gt. (missing_value + tiny(missing_value))))
   end if

   if (id_phno3_g_d>0) then
      used = send_data(id_phno3_g_d,trop_diag_array(:,:,:,trop_diag%ind_phno3_g_d)*pwt(:,:,:)*1.e3/WTMAIR,Time_next,is_in=is,js_in=js)
   end if

!-----------------------------------------------------------------------
!     ... special case(nox = no + no2)
!-----------------------------------------------------------------------
!  nno = get_tracer_index(MODEL_ATMOS,'no')
!  nno2 = get_tracer_index(MODEL_ATMOS,'no2')
!  if((nno /= 0) .and. (nno2 /= 0)) then
!     rno(:,:,:) = r(:,:,:,nno)/ MAX((r(:,:,:,nno) + r(:,:,:,nno2)),1.e-30)

!     call interpolator(ub, Time,phalf,ub_temp,'nox',is,js)

!     where(pfull(:,:,:) < ub_pres)
!        chem_dt(:,:,:,nno) =((rno(:,:,:)*ub_temp(:,:,:))-r(:,:,:,nno)) / relaxed_dt
!        chem_dt(:,:,:,nno2) = (((1.-rno(:,:,:))*ub_temp(:,:,:))-r(:,:,:,nno2)) / &
!             relaxed_dt
!     endwhere
!  end if

!-----------------------------------------------------------------------
!     ... Chemical families (Cly)
!-----------------------------------------------------------------------
   cly(:,:,:) = 0.
   if (cl_ndx>0) then
      cly(:,:,:) = cly(:,:,:) + r_temp(:,:,:,cl_ndx)
   end if
   if (clo_ndx>0) then
      cly(:,:,:) = cly(:,:,:) + r_temp(:,:,:,clo_ndx)
   end if
   if (hcl_ndx>0) then
      cly(:,:,:) = cly(:,:,:) + r_temp(:,:,:,hcl_ndx)
   end if
   if (hocl_ndx>0) then
      cly(:,:,:) = cly(:,:,:) + r_temp(:,:,:,hocl_ndx)
   end if
   if (clono2_ndx>0) then
      cly(:,:,:) = cly(:,:,:) + r_temp(:,:,:,clono2_ndx)
   end if
   if (cl2o2_ndx>0) then
      cly(:,:,:) = cly(:,:,:) + r_temp(:,:,:,cl2o2_ndx)*2
   end if
   if (cl2_ndx>0) then
      cly(:,:,:) = cly(:,:,:) + r_temp(:,:,:,cl2_ndx)*2
   end if
   if (clno2_ndx>0) then
      cly(:,:,:) = cly(:,:,:) + r_temp(:,:,:,clno2_ndx)
   end if
   if (brcl_ndx>0) then
      cly(:,:,:) = cly(:,:,:) + r_temp(:,:,:,brcl_ndx)
   end if

!-----------------------------------------------------------------------
!     ... Cly chemical tendency diagnostic
!-----------------------------------------------------------------------
   if (id_dclydt_chem>0) then
      used = send_data(id_dclydt_chem, (cly(:,:,:)-cly0(:,:,:))/dt, Time_next, is_in=is, js_in=js)
   end if

!-----------------------------------------------------------------------
!     ... Cly conservation
!-----------------------------------------------------------------------
   if (force_cly_conservation .or. rescale_cly_components) then
      if (rescale_cly_components) then
         cly_ratio(:,:,:) = r(:,:,:,cly_ndx) / MAX( cly(:,:,:), small )
         cly(:,:,:) = r(:,:,:,cly_ndx)
      else if (force_cly_conservation) then
         cly_ratio(:,:,:) = cly0(:,:,:) / MAX( cly(:,:,:), small )
         cly(:,:,:) = cly0(:,:,:)
      end if
      if (cl_ndx>0) then
         r_temp(:,:,:,cl_ndx) = r_temp(:,:,:,cl_ndx) * cly_ratio(:,:,:)
      end if
      if (clo_ndx>0) then
         r_temp(:,:,:,clo_ndx) = r_temp(:,:,:,clo_ndx) * cly_ratio(:,:,:)
      end if
      if (hcl_ndx>0) then
         r_temp(:,:,:,hcl_ndx) = r_temp(:,:,:,hcl_ndx) * cly_ratio(:,:,:)
      end if
      if (hocl_ndx>0) then
         r_temp(:,:,:,hocl_ndx) = r_temp(:,:,:,hocl_ndx) * cly_ratio(:,:,:)
      end if
      if (clono2_ndx>0) then
         r_temp(:,:,:,clono2_ndx) = r_temp(:,:,:,clono2_ndx) * cly_ratio(:,:,:)
      end if
      if (cl2o2_ndx>0) then
         r_temp(:,:,:,cl2o2_ndx) = r_temp(:,:,:,cl2o2_ndx) * cly_ratio(:,:,:)
      end if
      if (cl2_ndx>0) then
         r_temp(:,:,:,cl2_ndx) = r_temp(:,:,:,cl2_ndx) * cly_ratio(:,:,:)
      end if
      if (clno2_ndx>0) then
         r_temp(:,:,:,clno2_ndx) = r_temp(:,:,:,clno2_ndx) * cly_ratio(:,:,:)
      end if
      if (brcl_ndx>0) then
         r_temp(:,:,:,brcl_ndx) = r_temp(:,:,:,brcl_ndx) * cly_ratio(:,:,:)
      end if
   end if

!-----------------------------------------------------------------------
!     ... Chemical families (Bry, NOy)
!-----------------------------------------------------------------------
   bry(:,:,:) = 0.
   if (br_ndx>0) then
      bry(:,:,:) = bry(:,:,:) + r_temp(:,:,:,br_ndx)
   end if
   if (bro_ndx>0) then
      bry(:,:,:) = bry(:,:,:) + r_temp(:,:,:,bro_ndx)
   end if
   if (hbr_ndx>0) then
      bry(:,:,:) = bry(:,:,:) + r_temp(:,:,:,hbr_ndx)
   end if
   if (hobr_ndx>0) then
      bry(:,:,:) = bry(:,:,:) + r_temp(:,:,:,hobr_ndx)
   end if
   if (brono2_ndx>0) then
      bry(:,:,:) = bry(:,:,:) + r_temp(:,:,:,brono2_ndx)
   end if
   if (brcl_ndx>0) then
      bry(:,:,:) = bry(:,:,:) + r_temp(:,:,:,brcl_ndx)
   end if
   
!++van
   noy(:,:,:) = 0.
! Loop over total number of atmospheric tracers (nt), not just solver tracers (pcnstm1)
   call get_number_tracers(MODEL_ATMOS, num_tracers=nt)
   do n = 1,nt
      if ( nb_N_Ox(n) .gt. 0.) then
        call get_tracer_names (MODEL_ATMOS, n, tracer_name)
        if (tracer_name .eq. 'brono2') then
            noytracer = 'BrONO2' 
        else if (tracer_name .eq. 'clono2') then
            noytracer = 'ClONO2'
        else
            noytracer =  uppercase(tracer_name)
        end if
        idx = get_spc_ndx(noytracer)
        noy(:,:,:) = noy(:,:,:) + r_temp(:,:,:,idx)*nb_N_ox(n)
      end if
   end do
!--van


!-----------------------------------------------------------------------
!     ... stratospheric Cly and Bry source
!-----------------------------------------------------------------------
   if(age_ndx > 0 .and. age_ndx <= ntp) then
      call strat_chem_dcly_dt(Time, phalf, is, js, age, cly, bry, dclydt, dbrydt)
      do k = 1,kd
         where( coszen(:,:) > 0. )
            dclydt(:,:,k) = 2*dclydt(:,:,k)
            dbrydt(:,:,k) = 2*dbrydt(:,:,k)
         elsewhere
            dclydt(:,:,k) = 0.
            dbrydt(:,:,k) = 0.
         end where
      end do
   else
      dclydt(:,:,:) = 0.
      dbrydt(:,:,:) = 0.
   end if
   if (cl_ndx>0) then
      chem_dt(:,:,:,indices(cl_ndx)) = chem_dt(:,:,:,indices(cl_ndx)) + dclydt(:,:,:)
      used = send_data(id_dclydt, dclydt, Time_next, is_in=is, js_in=js)
   end if
   if (br_ndx>0) then
      chem_dt(:,:,:,indices(br_ndx)) = chem_dt(:,:,:,indices(br_ndx)) + dbrydt(:,:,:)
      used = send_data(id_dbrydt, dbrydt, Time_next, is_in=is, js_in=js)
   end if

!-----------------------------------------------------------------------
!     ... Set diagnostic tracers for chemical families
!-----------------------------------------------------------------------
   if (noy_ndx > ntp) then
      rdiag(:,:,:,noy_ndx-ntp) = noy(:,:,:)
   end if
   if (cly_ndx > ntp) then
      rdiag(:,:,:,cly_ndx-ntp) = cly(:,:,:) + dclydt(:,:,:)*dt
   else if (cly_ndx > 0) then
      chem_dt(:,:,:,cly_ndx) = dclydt(:,:,:)
   end if
   if (bry_ndx > ntp) then
      rdiag(:,:,:,bry_ndx-ntp) = bry(:,:,:) + dbrydt(:,:,:)*dt
   end if

!-----------------------------------------------------------------------
!     ... Photolysis rates
!-----------------------------------------------------------------------
   do n = 1,phtcnt
      if(id_jval(n)>0) then
         used = send_data(id_jval(n),jvals(:,:,:,n),Time_next,is_in=is,js_in=js)
      end if
   end do
   if (query_cmip_diag_id(ID_jno2) .and. jno2_ndx>0) then
      used = send_cmip_data_3d (ID_jno2, jvals(:,:,:,jno2_ndx), &
           Time_next, is_in=is, js_in=js, ks_in=1)
   end if
   if (query_cmip_diag_id(ID_jo1d) .and. jo1d_ndx>0) then
      used = send_cmip_data_3d (ID_jo1d, jvals(:,:,:,jo1d_ndx), &
           Time_next, is_in=is, js_in=js, ks_in=1)
   end if

!-----------------------------------------------------------------------
!     ... Kinetic reaction rates
!-----------------------------------------------------------------------
   do n = 1,gascnt
      if(id_rate_const(n)>0) then
         used = send_data(id_rate_const(n),rate_constants(:,:,:,n),Time_next,is_in=is,js_in=js)
      end if
   end do


!reaction on dust
   do n=1,max_dust
      if (usr_hno3_dust(n)>0 .and. id_rx_hno3_dust(n)>0) then
         used = send_data(id_rx_hno3_dust(n),rate_constants(:,:,:,usr_hno3_dust(n)),Time_next,is_in=is,js_in=js)
      end if
      if (usr_n2o5_dust(n)>0 .and. id_rx_n2o5_dust(n)>0) then
         used = send_data(id_rx_n2o5_dust(n),rate_constants(:,:,:,usr_n2o5_dust(n)),Time_next,is_in=is,js_in=js)
      end if
      if (usr_no3_dust(n)>0 .and. id_rx_no3_dust(n)>0) then
         used = send_data(id_rx_no3_dust(n),rate_constants(:,:,:,usr_no3_dust(n)),Time_next,is_in=is,js_in=js)
      end if
      if (usr_so2_dust(n)>0 .and. id_rx_so2_dust(n)>0) then
         used = send_data(id_rx_so2_dust(n),rate_constants(:,:,:,usr_so2_dust(n)),Time_next,is_in=is,js_in=js)
      end if
      if (usr_so4_dust(n)>0 .and. id_rx_so4_dust(n)>0) then
         used = send_data(id_rx_so4_dust(n),rate_constants(:,:,:,usr_so4_dust(n)),Time_next,is_in=is,js_in=js)
      end if
   end do


!-----------------------------------------------------------------------
!     ... Output diagnostics
!-----------------------------------------------------------------------
   used = send_data(id_volc_aer, strat_aerosol, Time_next, is_in=is, js_in=js)
   used = send_data(id_psc_sat, psc_vmr_save(:,:,:,1), Time_next, is_in=is, js_in=js)
   used = send_data(id_psc_nat, psc_vmr_save(:,:,:,2), Time_next, is_in=is, js_in=js)
   used = send_data(id_psc_ice, psc_vmr_save(:,:,:,3), Time_next, is_in=is, js_in=js)
   if (id_h2o_chem>0) then
      if (sphum_ndx>0) then
         used = send_data(id_h2o_chem, r_temp(:,:,:,sphum_ndx), Time_next, is_in=is, js_in=js)
      else
         used = send_data(id_h2o_chem, q(:,:,:)*WTMAIR/WTMH2O, Time_next, is_in=is, js_in=js)
      end if
   end if
   used = send_data(id_coszen, coszen_local(:,:), Time_next, is_in=is, js_in=js)
   used = send_data(id_imp_slv_nonconv,imp_slv_nonconv(:,:,:),Time_next,is_in=is,js_in=js)

!-----------------------------------------------------------------------
!     ... convert H2O VMR tendency to specific humidity tendency
!-----------------------------------------------------------------------
   if (sphum_ndx > 0) then
      n = indices(sphum_ndx)
      chem_dt(:,:,:,n) = chem_dt(:,:,:,n) * WTMH2O / WTMAIR
!     chem_dt(:,:,:,n) = 0.
   end if

   call mpp_clock_end(clock_id)

!-----------------------------------------------------------------------

end subroutine tropchem_driver
!</SUBROUTINE>

!#######################################################################

! <FUNCTION NAME="tropchem_driver_init">
!   <OVERVIEW>
!     Initializes the tropospheric chemistry driver.
!   </OVERVIEW>
!   <DESCRIPTION>
!     This subroutine initializes the tropospheric chemistry module.
!     It is called from atmos_tracer_driver_init.
!     Data sets are read in for dry deposition, upper boundary conditions,
!     and emissions. Off-line sulfate concentrations are also read in for
!     use in calculating heterogeneous reaction rates (if SO4 is not
!     included as a tracer).
!   </DESCRIPTION>
!   <TEMPLATE>
!     Ltropchem = tropchem_driver_init( r, mask, axes, Time, &
!                                       lonb_mod, latb_mod, phalf)
!   </TEMPLATE>
!   <IN NAME="mask" TYPE="real, optional" DIM="(:,:,:)">
!      optional mask that designates which grid points
!      are above (1) or below (0) the ground
!   </IN>
!   <IN NAME="axes" TYPE="integer" DIM="(4)">
!     The axes relating to the tracer array
!   </IN>
!   <IN NAME="Time" TYPE="type(time_type)">
!     Model time.
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
!   <OUT NAME="drydep_data" TYPE="interpolate_type" DIM="(:)">
!     Tracer dry deposition velocities
!   </OUT>
!   <OUT NAME="Ltropchem" TYPE="logical">
!     Do tropospheric chemistry? (Output as function value)
!   </OUT>
!   <INOUT NAME="r" TYPE="real" DIM="(:,:,:,:)">
!     Tracer mixing ratios (tropchem tracers in VMR)
!   </INOUT>

function tropchem_driver_init( domain, r, mask, axes, Time, &
                               lonb_mod, latb_mod, phalf) result(Ltropchem)

!-----------------------------------------------------------------------
!
!   r    = tracer fields dimensioned as (nlon,nlat,nlev,ntrace)
!   mask = optional mask (0. or 1.) that designates which grid points
!          are above (=1.) or below (=0.) the ground dimensioned as
!          (nlon,nlat,nlev).
!
!-----------------------------------------------------------------------
   type(domain2D),target,intent(in) :: domain !< Atmosphere domain
   real, intent(inout), dimension(:,:,:,:) :: r
   real, intent(in),    dimension(:,:,:), optional :: mask
   type(time_type), intent(in) :: Time
   integer        , intent(in) :: axes(4)
   real, intent(in), dimension(:,:) :: lonb_mod
   real, intent(in), dimension(:,:) :: latb_mod
   real, intent(in),dimension(:,:,:) :: phalf

   real    :: small_value

   logical :: Ltropchem
   integer :: flag_file, flag_spec, flag_fixed
   integer :: n, i, nt
   integer :: ierr, io, logunit
   character(len=64) :: nc_file,filename,specname
   character(len=256) :: control=''
   character(len=64) :: name=''
   type(interpolate_type) :: init_conc
   character(len=64),dimension(pcnstm1) :: emis_files = '', &
                                           emis3d_files = '', &
                                           conc_files = '', &
                                           ub_files = '', &
                                           lb_files = '', &
                                           dry_files, &
                                           wet_ind, &
                                           conc_names, &
                                           dry_names, &
                                           airc_files
   logical :: tracer_initialized

   character(len=16) :: fld
   character(len=32) :: tracer_name

   integer :: flb, series_length, year, diy
   real :: input_time
   real :: scale_factor, extra_seconds, fixed_year
   type(time_type) :: Year_t
!
!-----------------------------------------------------------------------
!
   if (module_is_initialized) return

!-----------------------------------------------------------------------
!     ... write version number
!-----------------------------------------------------------------------
   call write_version_number(version, tagname)

!-----------------------------------------------------------------------
!     ... read namelist
!-----------------------------------------------------------------------
   if(file_exists('input.nml')) then
      read (input_nml_file, nml=tropchem_driver_nml, iostat=io)
      ierr = check_nml_error(io,'tropchem_driver_nml')
   end if

   logunit = stdlog()
   if(mpp_pe() == mpp_root_pe()) then
      write(logunit, nml=tropchem_driver_nml)
      verbose = verbose + 1
   end if

   Ltropchem = do_tropchem
   if (.not. Ltropchem) then
      return
   end if

  tropchem_domain => domain

!-------------------------------------------------------------------------
!     ... Make sure input value for clouds_in_fastjx is a valid option.
!-------------------------------------------------------------------------
   if (trim(clouds_in_fastjx) == 'none') then
     use_lsc_in_fastjx = .false.
   else if (trim(clouds_in_fastjx) == 'lsc_only') then
     use_lsc_in_fastjx = .true.
   else
     call error_mesg ('tropchem_driver_init', &
                     ' invalid string for clouds_in_fastjx', FATAL)
   endif


if ( (gNO3_dust .gt. 0. .or. gN2O5_dust .gt. 0.) .and. .not. het_chem_fine_aerosol_only ) then
   call error_mesg ('tropchem_driver_init', 'uptake on dust + all aerosol => double counting', FATAL )
end if

if ( trim(cloud_chem_type) == 'f1p' ) then
  trop_option%cloud_chem = CLOUD_CHEM_F1P
else
  call error_mesg ('trop_chem_driver', 'unknown cloud chem type', FATAL)
end if

trop_option%verbose = verbose

!cloud chem pH solver


if ( trim(cloud_chem_pH_solver) == 'am3' ) then
   trop_option%cloud_chem_pH_solver  = CLOUD_CHEM_PH_LEGACY
   if(mpp_pe() == mpp_root_pe()) write(*,*) 'legacypH'
elseif ( trim(cloud_chem_pH_solver) == 'bisection' ) then
   trop_option%cloud_chem_pH_solver   = CLOUD_CHEM_PH_BISECTION
elseif ( trim(cloud_chem_pH_solver) == 'cubic' ) then
   trop_option%cloud_chem_pH_solver   = CLOUD_CHEM_PH_CUBIC
else
   call error_mesg ('tropchem_driver_init', 'undefined cloud chem', FATAL )
end if

if ( cloud_pH .lt. 0 ) then
   trop_option%cloud_H = -999
else
   trop_option%cloud_H = 10**(-cloud_pH)
end if

   if(mpp_pe() == mpp_root_pe()) write(*,*) 'cloud_H',trop_option%cloud_H


if ( trim(het_chem_type) == 'legacy' ) then
   trop_option%het_chem = HET_CHEM_LEGACY
   if(mpp_pe() == mpp_root_pe()) write(*,*) 'Using legacy heterogeneous chemistry'
elseif ( trim(het_chem_type) == 'j1m' ) then
   trop_option%het_chem = HET_CHEM_J1M
   if(mpp_pe() == mpp_root_pe()) write(*,*) 'Using new (J1M) heterogeneous chemistry'
end if

! Chemical pre-processor input filename
trop_option%sim_data_flsp = sim_data_filename
if(mpp_pe() == mpp_root_pe()) write(*,*) 'Chemical pre-processor input filename: ', TRIM(sim_data_filename)

!gammas to be added when het chem is working
trop_option%gN2O5                    = gN2O5
if(mpp_pe() == mpp_root_pe())    write(*,*)     "gN2O5:",trop_option%gN2O5
trop_option%gNO3                     = gNO3
if(mpp_pe() == mpp_root_pe())    write(*,*)     "gNO3:",trop_option%gNO3
trop_option%gNO2                     = gNO2
if(mpp_pe() == mpp_root_pe())    write(*,*)     "gNO2:",trop_option%gNO2
trop_option%gHO2                     = gHO2
if(mpp_pe() == mpp_root_pe())    write(*,*)     "gHO2:",trop_option%gHO2
trop_option%gNH3                     = gNH3
if(mpp_pe() == mpp_root_pe())    write(*,*)     "gNH3:",trop_option%gNH3


trop_option%gHNO3_dust               = gHNO3_dust
trop_option%gN2O5_dust               = gN2O5_dust
trop_option%gSO2_dust                = gSO2_dust
trop_option%gNO3_dust                = gNO3_dust
trop_option%scale_dust_uptake        = scale_dust_uptake

trop_option%do_fastjx_photo = do_fastjx_photo
trop_option%min_lwc_for_cloud_chem = min_lwc_for_cloud_chem
trop_optioN%check_convergence = check_convergence
trop_optioN%use_lsc_in_fastjx = use_lsc_in_fastjx
trop_option%het_chem_fine_aerosol_only = het_chem_fine_aerosol_only
trop_option%cloud_ho2_h2o2 = cloud_ho2_h2o2
trop_option%max_rh_aerosol = max_rh_aerosol
trop_option%limit_no3      = limit_no3
trop_option%frac_aerosol_incloud = frac_aerosol_incloud
trop_option%time_varying_solarflux = time_varying_solarflux

trop_option%gSO2                     = gSO2
if(mpp_pe() == mpp_root_pe()) write(*,*) 'gSO2: ',trop_option%gSO2
if (trim(gso2_dynamic).eq.'none') then
   trop_option%gSO2_dynamic             = -1
else if (trim(gso2_dynamic).eq.'wang2014') then
   trop_option%gSO2_dynamic             = GSO2_WANG2014   
   !http://onlinelibrary.wiley.com/doi/10.1002/2013JD021426/full  
else if (trim(gso2_dynamic).eq.'zheng2015') then
   trop_option%gSO2_dynamic             = GSO2_ZHENG2015
!   http://www.atmos-chem-phys.net/15/2031/2015/
else if (trim(gso2_dynamic).eq.'zheng2015_low') then
   trop_option%gSO2_dynamic             = GSO2_ZHENG2015_LOW
!   http://www.atmos-chem-phys.net/15/2031/2015/
end if

trop_option%NO2_SO2_max = NO2_SO2_max

if(mpp_pe() == mpp_root_pe()) write(*,*) 'gso2_dynamic case:',trop_option%gSO2_dynamic

if (trim(ghno3_dust_dynamic).eq.'none') then
   trop_option%ghno3_dust_dynamic = -1
elseif (trim(ghno3_dust_dynamic).eq.'dynamic') then
   trop_option%ghno3_dust_dynamic = 1
end if


trop_option%modulate_frac_ic = modulate_frac_ic
if (mpp_pe()==mpp_root_pe()) write(*,*) 'modulate_frac_ic',modulate_frac_ic

!aerosol thermo
if    ( trim(aerosol_thermo_method)   == 'legacy' ) then
   trop_option%aerosol_thermo = AERO_LEGACY
   if(mpp_pe() == mpp_root_pe()) write(*,*) 'legacy no3'
elseif ( trim(aerosol_thermo_method)   == 'isorropia' ) then
   trop_option%aerosol_thermo = AERO_ISORROPIA
elseif ( trim(aerosol_thermo_method)   == 'no_thermo' ) then
   trop_option%aerosol_thermo = NO_AERO
else
   call error_mesg ('tropchem_driver_init', 'undefined aerosol thermo', FATAL )
end if

trop_option%rh_het_max    = rh_het_max
if (mpp_pe()==mpp_root_pe()) then
   write(*,*) 'rh_het_max',trop_option%rh_het_max
end if

!-----------------------------------------------------------------------
!     ... Setup sulfate input/interpolation
!-----------------------------------------------------------------------
   call interpolator_init( sulfate, trim(file_sulfate), lonb_mod, latb_mod, &
                           data_out_of_bounds=(/CONSTANT/),      &
                           vert_interp=(/INTERP_WEIGHTED_P/) )

!-----------------------------------------------------------------------
!     ... Initialize chemistry driver
!-----------------------------------------------------------------------
   call chemini( file_jval_lut, file_jval_lut_min, use_tdep_jvals, &
                 o3_column_top, jno_scale_factor, verbose,   &
                 do_fastjx_photo, trop_option)

!-----------------------------------------------------------------------
!     ... set initial value of indices
!-----------------------------------------------------------------------
   indices(:) = 0
   do i=1,pcnstm1
      n = get_tracer_index(MODEL_ATMOS, tracnam(i))
      if (trim(tracnam(i)) == 'H2O') then
         if (n <= 0) then
            n = get_tracer_index(MODEL_ATMOS, 'sphum')
         end if
         sphum_ndx = i
         do_interactive_h2o = .true.
      end if
      if (n >0) then
         indices(i) = n
         if (indices(i) > 0 .and. mpp_pe() == mpp_root_pe()) then
            write(*,30) tracnam(i),indices(i)
            write(logunit,30) trim(tracnam(i)),indices(i)
         end if
      else
!<ERROR MSG="Tropospheric chemistry tracer not found in field table" STATUS="WARNING">
!   A tropospheric chemistry tracer was not included in the field table
!</ERROR>
         call error_mesg ('tropchem_driver_init', trim(tracnam(i)) // ' is not found', WARNING)
      end if
   end do
30 format (A,' was initialized as tracer number ',i3)

   nco2 = get_tracer_index(MODEL_ATMOS, 'co2' )
   if (nco2 > 0 .and. mpp_pe() == mpp_root_pe()) then
      call error_mesg ('tropchem_driver', 'CO2 is active',NOTE)
   endif
   cl_ndx     = get_spc_ndx('Cl')
   clo_ndx    = get_spc_ndx('ClO')
   hcl_ndx    = get_spc_ndx('HCl')
   hocl_ndx   = get_spc_ndx('HOCl')
   clono2_ndx = get_spc_ndx('ClONO2')
   cl2o2_ndx  = get_spc_ndx('Cl2O2')
   cl2_ndx    = get_spc_ndx('Cl2')
   clno2_ndx  = get_spc_ndx('ClNO2')
   br_ndx     = get_spc_ndx('Br')
   bro_ndx    = get_spc_ndx('BrO')
   hbr_ndx    = get_spc_ndx('HBr')
   hobr_ndx   = get_spc_ndx('HOBr')
   brono2_ndx = get_spc_ndx('BrONO2')
   brcl_ndx   = get_spc_ndx('BrCl')
   hno3_ndx   = get_spc_ndx('HNO3')
   no_ndx     = get_spc_ndx('NO')
   no2_ndx    = get_spc_ndx('NO2')
   no3_ndx    = get_spc_ndx('NO3')
   n_ndx      = get_spc_ndx('N')
   n2o5_ndx   = get_spc_ndx('N2O5')
   ho2no2_ndx = get_spc_ndx('HO2NO2')
   pan_ndx    = get_spc_ndx('PAN')
   onit_ndx   = get_spc_ndx('ONIT')
   mpan_ndx   = get_spc_ndx('MPAN')
   isopno3_ndx= get_spc_ndx('ISOPNO3')
   onitr_ndx  = get_spc_ndx('ONITR')
   o3_ndx     = get_spc_ndx('O3')
   ch4_ndx    = get_spc_ndx('CH4')
   dms_ndx    = get_spc_ndx('DMS')
   so4_ndx(1) = get_spc_ndx('SO4')
   so4_ndx(2) = get_spc_ndx('SO4_D1')
   so4_ndx(3) = get_spc_ndx('SO4_D2')
   so4_ndx(4) = get_spc_ndx('SO4_D3')
   so4_ndx(5) = get_spc_ndx('SO4_D4')
   so4_ndx(6) = get_spc_ndx('SO4_D5')
   co_ndx     = get_spc_ndx('CO')
   n2o_ndx    = get_spc_ndx('N2O')

   o3s_ndx     = get_spc_ndx('O3S')
   o3s_e90_ndx = get_spc_ndx('O3S_E90')
   e90_ndx     = get_tracer_index(MODEL_ATMOS,'e90')
   extinct_ndx = get_tracer_index(MODEL_ATMOS, 'Extinction')
   noy_ndx     = get_tracer_index(MODEL_ATMOS, 'NOy')
   cly_ndx     = get_tracer_index(MODEL_ATMOS, 'Cly')
   bry_ndx     = get_tracer_index(MODEL_ATMOS, 'Bry')

   jno2_ndx    = get_rxt_ndx( 'jno2' )
   jo1d_ndx    = get_rxt_ndx( 'jo1d' )

!-----------------------------------------------------------------------
!     ... Check Cly settings
!-----------------------------------------------------------------------
   if (rescale_cly_components) then
      if (cly_ndx == NO_TRACER .or. .not. check_if_prognostic(MODEL_ATMOS,cly_ndx)) then
         call error_mesg ('tropchem_driver_init', &
                          'rescale_cly_components=T requires Cly to be registered as a prognostic tracer', FATAL)
      end if
      if (force_cly_conservation) then
         call error_mesg ('tropchem_driver_init', &
                          'rescale_cly_components=T incompatible with force_cly_conservation=T setting', FATAL)
      end if
   end if

!++van
!----------------------------------------
!     ... For calculating NOy
!----------------------------------------

    call get_number_tracers(MODEL_ATMOS, num_tracers=nt)
    allocate(nb_N_Ox(nt)) 
    if(mpp_pe() == mpp_root_pe()) then
       write (*,*) 'NOTE: tropchem_driver_init, nt = ', nt
       write (*,*) 'NOy is composed of :'
    end if
    do n = 1,nt
      call  get_chem_param (n, nb_N_Ox=nb_N_Ox(n))
      if ( nb_N_Ox(n) .gt. 0.) then
        call get_tracer_names (MODEL_ATMOS, n, tracer_name)
        if(mpp_pe() == mpp_root_pe()) write (*,'(2a,g14.6)') trim(tracer_name),', nb_N_ox=',nb_N_ox(n)
      end if
    end do
!--van    

!-----------------------------------------------------------------------
!     ... Setup upper boundary condition data
!-----------------------------------------------------------------------
   call interpolator_init( ub_default, trim(file_ub), lonb_mod, latb_mod, &
                           data_out_of_bounds=(/CONSTANT/),          &
                           vert_interp=(/INTERP_WEIGHTED_P/) )

!-----------------------------------------------------------------------
!     ... Set up concentration input/interpolation
!-----------------------------------------------------------------------
   call interpolator_init( conc, trim(file_conc), lonb_mod, latb_mod, &
                           data_out_of_bounds=(/CONSTANT/),&
                           vert_interp=(/INTERP_WEIGHTED_P/) )

!-----------------------------------------------------------------------
!     ... Set up aircraft emissions interpolation
!-----------------------------------------------------------------------
   call interpolator_init( airc_default, trim(file_aircraft), lonb_mod, latb_mod,&
                           data_out_of_bounds=(/CONSTANT/), &
                           vert_interp=(/INTERP_WEIGHTED_P/))

!-----------------------------------------------------------------------
!     ... Setup emissions input/interpolation
!-----------------------------------------------------------------------
   do i = 1,pcnstm1
      nc_file = trim(file_emis_1)//lowercase(trim(tracnam(i)))//trim(file_emis_2)
      call init_emis_data( inter_emis(i), MODEL_ATMOS, 'emissions', indices(i), nc_file, &
                           lonb_mod, latb_mod, emis_field_names(i), &
                           has_emis(i), diurnal_emis(i), axes, Time, land_does_emission(i) )
      if( has_emis(i) ) emis_files(i) = trim(nc_file)

!-----------------------------------------------------------------------
!     ... Vertically-distributed emissions
!-----------------------------------------------------------------------
      nc_file = trim(file_emis3d_1)//lowercase(trim(tracnam(i)))//trim(file_emis3d_2)
      call init_emis_data( inter_emis3d(i), MODEL_ATMOS, 'emissions3d', indices(i), nc_file, &
                           lonb_mod, latb_mod, emis3d_field_names(i), &
                           has_emis3d(i), diurnal_emis3d(i), axes, Time )
      if( has_emis3d(i) ) emis3d_files(i) = trim(nc_file)

!-----------------------------------------------------------------------
!     ... Interactive emissions
!-----------------------------------------------------------------------
      if ( .not. do_terpene_emis_bug .or. trim(tracnam(i))=="DMS" ) then
         call init_xactive_emis( MODEL_ATMOS, 'xactive_emissions', indices(i), tracnam(i), &
                                 axes, Time, lonb_mod, latb_mod, phalf, &
                                 has_xactive_emis(i), id_xactive_emis(i), mask )
      endif
!-----------------------------------------------------------------------
!     ... Upper boundary condition
!-----------------------------------------------------------------------
      if( query_method('upper_bound', MODEL_ATMOS,indices(i),name,control) ) then
         if( trim(name)=='file' ) then
            flag_file = parse(control, 'file',filename)
            flag_spec = parse(control, 'name',specname)

            if( flag_file > 0 .and. trim(filename) /= trim(file_ub) ) then
               ub_files(i) = trim(filename)
               call interpolator_init(ub(i), trim(filename), lonb_mod, latb_mod, &
                       data_out_of_bounds=(/CONSTANT/),          &
                       vert_interp=(/INTERP_WEIGHTED_P/))
            else
               ub_files(i) = trim(file_ub)
               ub(i) = ub_default
            end if
            if(flag_spec > 0) then
               ub_names(i) = trim(specname)
            else
               ub_names(i) = trim(lowercase(tracnam(i)))
            end if

            has_ubc(i) = .true.

         end if
      end if

!-----------------------------------------------------------------------
!     ... Lower boundary condition
!-----------------------------------------------------------------------
      lbc_entry(i) = get_base_time()
      if( query_method('lower_bound', MODEL_ATMOS,indices(i),name,control) ) then
         if( trim(name)=='file' ) then

            flag_spec = parse(control,'dry',specname)
            if (flag_spec>0) lbc_dry(i) = .true.
            
            flag_file = parse(control, 'file', filename)
            flag_spec = parse(control, 'factor', scale_factor)
            flag_fixed = parse(control, 'fixed_year', fixed_year)
            
            if( flag_file > 0 ) then
               lb_files(i) = 'INPUT/' // trim(filename)
               if( file_exists(lb_files(i)) ) then
                  flb = tropchem_open_file( lb_files(i) )
                  read(flb, FMT='(i12)') series_length
                  allocate( lb(i)%gas_value(series_length), &
                            lb(i)%gas_time(series_length) )
!---------------------------------------------------------------------
!    convert the time stamps of the series to time_type variables.
!---------------------------------------------------------------------
                  do n = 1,series_length
                     read (flb, FMT = '(2f12.4)') input_time, lb(i)%gas_value(n)
                     year = INT(input_time)
                     Year_t = set_date(year,1,1,0,0,0)
                     diy = days_in_year (Year_t)
                     extra_seconds = (input_time - year)*diy*SECONDS_PER_DAY
                     lb(i)%gas_time(n) = Year_t + set_time(NINT(extra_seconds), 0)
                  end do
                  if (flag_spec > 0) then
                     lb(i)%gas_value(:) = lb(i)%gas_value(:) * scale_factor
                  end if
                  close ( flb )
                  if( flag_fixed > 0 ) then
                     fixed_lbc_time(i) = .true.
                     year = INT(fixed_year)
                     Year_t = set_date(year,1,1,0,0,0)
                     diy = days_in_year (Year_t)
                     extra_seconds = (fixed_year - year)*diy*SECONDS_PER_DAY
                     lbc_entry(i) = Year_t + set_time(NINT(extra_seconds), 0)
                  end if                  
               else
                  call error_mesg ('tropchem_driver_init', &
                                   'Failed to find input file '//trim(lb_files(i)), FATAL)
               end if
            else
               flag_file = parse(control, 'file_nc', filename)
               lb_files(i) =  trim(filename)
               if (flag_file>0) then
                  if( file_exists('INPUT/' // lb_files(i)) ) then
                     call interpolator_init( lbc_interp(i), trim(lb_files(i)), lonb_mod, latb_mod, &
                          data_out_of_bounds=(/CONSTANT/),      &
                          vert_interp=(/INTERP_WEIGHTED_P/) )
                     has_lbc_2d(i) = .True.

                     flag_spec = parse(control, 'factor',scale_factor)
                     
                     if(flag_spec > 0) then
                        lbc_factor(i) = scale_factor
                     else
                        lbc_factor(i) = 1.
                     end if
                     
                     flag_spec = parse(control, 'name',specname)
                     
                     if(flag_spec > 0) then
                        lbc_names(i) = trim(specname)
                     else
                        lbc_names(i) = trim(lowercase(tracnam(i)))
                     end if

                     if (mpp_root_pe().eq.mpp_pe()) write(*,*) 'fixed_year',flag_fixed,fixed_year
                     
                     if( flag_fixed > 0 ) then                        
                        fixed_lbc_time(i) = .true.
                        year = INT(fixed_year)
                        Year_t = set_date(year,1,1,0,0,0)
                        diy = days_in_year (Year_t)
                        extra_seconds = (fixed_year - year)*diy*SECONDS_PER_DAY
                        lbc_entry(i) = Year_t + set_time(NINT(extra_seconds), 0)
                     end if
                     
                  else
                     call error_mesg ('tropchem_driver_init', &
                          'Failed to find input file '//trim(lb_files(i))//' '//trim(control)//' '//tracnam(i), FATAL)
                  end if
               else
                  call error_mesg ('tropchem_driver_init', 'Tracer '//trim(lowercase(tracnam(i)))// &
                       ' has lower_bound specified without a filename', FATAL)
               end if
            end if
            has_lbc(i) = .true.
         end if
      end if

!-----------------------------------------------------------------------
!     ... Initial conditions
!-----------------------------------------------------------------------
      tracer_initialized = check_if_tracer_initialized(lowercase(tracnam(i)))

      if(.not. tracer_initialized) then
         if( query_method('init_conc',MODEL_ATMOS,indices(i),name,control) ) then
            if( trim(name) == 'file' ) then
               flag_file = parse(control, 'file',filename)
               flag_spec = parse(control, 'name',specname)

               if( flag_file>0 .and. trim(filename) /= trim(file_conc) ) then
                  conc_files(i) = trim(filename)
                  call interpolator_init( init_conc,trim(filename),lonb_mod,latb_mod,&
                                          data_out_of_bounds=(/CONSTANT/), &
                                          vert_interp=(/INTERP_WEIGHTED_P/) )
               else
                  conc_files(i) = trim(file_conc)
                  init_conc = conc
               end if

               if( flag_spec > 0 ) then
                  conc_names(i) = trim(lowercase(specname))
                  specname = lowercase(specname)
               else
                  conc_names(i) = trim(lowercase(tracnam(i)))
                  specname = trim(lowercase(tracnam(i)))
               end if

               call interpolator(init_conc, Time, phalf,r(:,:,:,indices(i)),trim(specname))
            end if
         end if
      end if

!-----------------------------------------------------------------------
!     ... Aircraft emissions
!-----------------------------------------------------------------------
      if( query_method('aircraft_emis',MODEL_ATMOS,indices(i),name,control) ) then
         has_airc(i) = .true.
         if( trim(name) == 'file' ) then
            flag_file = parse(control,'file',filename)
            flag_spec = parse(control,'name',specname)

            if( flag_file >0 .and. trim(filename) /= trim(lowercase(file_aircraft)) ) then
               airc_files(i) = trim(filename)
               call interpolator_init( inter_aircraft_emis(i),trim(filename), lonb_mod, latb_mod, &
                                       data_out_of_bounds=(/CONSTANT/), &
                                       vert_interp=(/INTERP_WEIGHTED_P/) )
            else
               airc_files(i) = trim(file_aircraft)
               inter_aircraft_emis(i) = airc_default
            end if

            if( flag_spec >0 ) then
               airc_names(i) = trim(specname)
            else
               airc_names(i) = trim(lowercase(tracnam(i)))
            end if
         end if
      end if


!-----------------------------------------------------------------------
!     ... Wet deposition
!-----------------------------------------------------------------------
      if( query_method('wet_deposition',MODEL_ATMOS,indices(i),name,control) ) then
         wet_ind(i) = 'This species has wet deposition'
      else
         wet_ind(i) = ''
      end if

   end do

!move CO2 input out of the loop of "do i = 1,pcnstm1", 2016-07-25
!fp
!CO2
      if ( file_exists('INPUT/' // trim(co2_filename) ) ) then
         co2_t%use_fix_value  = .false.
         !read from file
         flb = tropchem_open_file( 'INPUT/' // trim(co2_filename) )
         read(flb,FMT='(i12)') series_length
         allocate( co2_t%gas_value(series_length), co2_t%gas_time(series_length) )
         do n = 1,series_length
            read (flb, FMT = '(2f12.4)') input_time, co2_t%gas_value(n)
            year = INT(input_time)
            Year_t = set_date(year,1,1,0,0,0)
            diy = days_in_year (Year_t)
            extra_seconds = (input_time - year)*diy*SECONDS_PER_DAY
            co2_t%gas_time(n) = Year_t + set_time(NINT(extra_seconds), 0)
         end do
         close(flb)         
         if (co2_scale_factor .gt. 0) then
            co2_t%gas_value = co2_t%gas_value * co2_scale_factor
         end if
         if (co2_fixed_year .gt. 0) then
            co2_t%use_fix_time = .true.
            year = INT(co2_fixed_year)
            Year_t = set_date(year,1,1,0,0,0)
            diy = days_in_year (Year_t)
            extra_seconds = (co2_fixed_year - year)*diy*SECONDS_PER_DAY
            co2_t%fixed_entry = Year_t + set_time(NINT(extra_seconds), 0)
         end if
      else
         co2_t%use_fix_value  = .true.
         co2_t%fixed_value    = co2_fixed_value
      end if

!-----------------------------------------------------------------------
!     ... Print out settings for tracer
!-----------------------------------------------------------------------
   if( mpp_pe() == mpp_root_pe() ) then
      write(logunit,*) '---------------------------------------------------------------------------------------'
      do i = 1,pcnstm1
         write(logunit,*) 'The tracname index is ',i
         write(logunit,*) 'The tracname is ',tracnam(i)
         if(check_if_prognostic(MODEL_ATMOS,indices(i))) then
            write(logunit,*) 'This is a prognostic tracer.'
         else
            write(logunit,*) 'This is a diagnostic tracer.'
         end if
         if(has_emis(i)) then
            write(logunit,*)'Emissions from file: ',trim(emis_files(i))
            if ( land_does_emission(i) ) then
               if (get_tracer_index(MODEL_LAND,trim(tracnam(i)))<=0) then
                  call error_mesg('atmos_tracer_utilities_init', &
                       'Emission of atmospheric tracer //"'//trim(tracnam(i))//&
                       '" is done on land side, but corresponding land tracer is not defined in the field table.', FATAL)
                  write(logunit,*) 'Emissions done in land'
               endif
            end if
         end if
         if(has_emis3d(i)) then
            write(logunit,*)'3-D Emissions from file: ',trim(emis3d_files(i))
         end if
         if(has_ubc(i)) then
            write(logunit,*)'Upper BC from file: ',trim(ub_files(i)), &
                             ', with the name of ',trim(ub_names(i))
         end if
         if(has_lbc(i)) then
            write(logunit,*)'Lower BC from file: ',trim(lb_files(i))
            if (fixed_lbc_time(i)) then
               write(logunit,*) '... with fixed year'
            end if
            if (has_lbc_2d(i)) then
               write(logunit,*) '... with 2d file'
            end if
         end if
         if(conc_files(i) /= '') then
            write(logunit,*)'Concentration from file: ',trim(conc_files(i)), &
                             ', with the name of ',trim(conc_names(i))
         end if
         if(dry_files(i) /= '') then
            write(logunit,*)'Dry deposition velocity from file: ',trim(dry_files(i)), &
                             ' with the name of '//trim(dry_names(i))
         end if
         if(wet_ind(i) /= '') then
            write(logunit,*) wet_ind(i)
         end if
         if(has_airc(i)) then
            write(logunit,*)'Aircraft emissions from file: ',trim(airc_files(i)), &
                             ' with the name of '//trim(airc_names(i))
         end if
         write(logunit,*) '---------------------------------------------------------------------------------------'
      end do
   end if


!-----------------------------------------------------------------------
!     ... Get the index number for the cloud variables
!-----------------------------------------------------------------------
   inqa = get_tracer_index(MODEL_ATMOS,'cld_amt') ! cloud fraction
   inql = get_tracer_index(MODEL_ATMOS,'liq_wat') ! cloud liquid specific humidity
   inqi = get_tracer_index(MODEL_ATMOS,'ice_wat') ! cloud ice water specific humidity

   age_ndx = get_tracer_index(MODEL_ATMOS,'age')  ! age tracer

!-----------------------------------------------------------------------
!     ... Call the chemistry hook init routine
!-----------------------------------------------------------------------
   call moz_hook_init( lght_no_prd_factor, normalize_lght_no_prd_area, &
                       allow_small_storms_lght_no_prd, min_land_frac_lght, &
                       Time, axes, verbose )

!-----------------------------------------------------------------------
!     ... Initializations for stratospheric chemistry
!-----------------------------------------------------------------------
!++lwh
   if (set_min_h2o_strat) then
      if (ch4_ndx>0) then
         if (.not. has_lbc(ch4_ndx)) then
            call error_mesg ('Tropchem_driver','set_min_h2o_strat=T, but LBC not set for CH4', FATAL)
         end if
      else
         call error_mesg ('Tropchem_driver','set_min_h2o_strat=T, but CH4 not included in chemistry solver', FATAL)
      end if
   end if
   call strat_chem_utilities_init( lonb_mod, latb_mod, &
                                   strat_chem_age_factor, strat_chem_dclydt_factor, &
                                   set_min_h2o_strat, ch4_filename, ch4_scale_factor, &
                                   fixed_lbc_time(ch4_ndx), lbc_entry(ch4_ndx), &
                                   cfc_lbc_filename, time_varying_cfc_lbc, cfc_lbc_dataset_entry)
!--lwh
   id_dclydt      = register_diag_field( module_name, 'cly_chem_dt', axes(1:3), Time, 'cly_chem_dt', 'VMR/s' )
   id_dclydt_chem = register_diag_field( module_name, 'cly_chem_dt_diag', axes(1:3), Time, 'cly_chem_dt_diag', 'VMR/s' )
   id_dbrydt      = register_diag_field( module_name, 'bry_chem_dt', axes(1:3), Time, 'bry_chem_dt', 'VMR/s' )

   id_volc_aer = register_diag_field( module_name, 'volc_aer_SA', axes(1:3), Time, 'volcanic_aerosol_surface_area', 'cm2/cm3' )
   id_psc_sat  = register_diag_field( module_name, 'psc_sat', axes(1:3), Time, 'psc_sat', 'VMR' )
   id_psc_nat  = register_diag_field( module_name, 'psc_nat', axes(1:3), Time, 'psc_nat', 'VMR' )
   id_psc_ice  = register_diag_field( module_name, 'psc_ice', axes(1:3), Time, 'psc_ice', 'VMR' )
   id_h2o_chem = register_diag_field( module_name, 'h2o_chem', axes(1:3), Time, 'h2o_chem', 'VMR' )

!-----------------------------------------------------------------------
!     ... Initialize additional diagnostics
!-----------------------------------------------------------------------
   id_sul = register_diag_field( module_name, 'sulfate', axes(1:3), Time, 'sulfate', 'VMR' )
   id_coszen = register_diag_field( module_name, 'coszen_tropchem', axes(1:2), Time, &
                                             'cosine_sza_tropchem', 'none' )
   id_imp_slv_nonconv = register_diag_field( module_name, 'imp_slv_nonconv', axes(1:3), Time, &
                                             'tropchem_implicit_solver_not_converged', 'VMR' )
   id_srf_o3 = register_diag_field( module_name, 'o3_srf', axes(1:2), Time, 'o3_srf', 'VMR' )

!-----------------------------------------------------------------------
!     ... Register diagnostic fields for species tendencies
!-----------------------------------------------------------------------
   id_co_emis_cmip =     &
        register_diag_field( module_name, 'co_emis_cmip', axes(1:2), &
                             Time, 'co_emis_cmip', 'kg/m2/s')
   id_no_emis_cmip =     &
        register_diag_field( module_name, 'no_emis_cmip', axes(1:2), &
                            Time, 'no_emis_cmip', 'kg/m2/s')
   id_so2_emis_cmip =     &
        register_diag_field( module_name, 'so2_emis_cmip', axes(1:2), &
                             Time, 'so2_emis_cmip', 'kg/m2/s')
   id_nh3_emis_cmip =     &
        register_diag_field( module_name, 'nh3_emis_cmip', axes(1:2), &
                            Time, 'nh3_emis_cmip', 'kg/m2/s')

   id_co_emis_cmip2 =     &
        register_diag_field( module_name, 'co_emis_cmip2', axes(1:2), &
                             Time, 'co_emis_cmip2', 'kg/m2/s')
   id_no_emis_cmip2 =     &
        register_diag_field( module_name, 'no_emis_cmip2', axes(1:2), &
                            Time, 'no_emis_cmip2', 'kg/m2/s')
   id_so2_emis_cmip2 =     &
        register_diag_field( module_name, 'so2_emis_cmip2', axes(1:2), &
                             Time, 'so2_emis_cmip2', 'kg/m2/s')
   id_nh3_emis_cmip2 =     &
        register_diag_field( module_name, 'nh3_emis_cmip2', axes(1:2), &
                            Time, 'nh3_emis_cmip2', 'kg/m2/s')
   !---- register cmip-named variables ----
   id_emico = register_cmip_diag_field_2d ( module_name, 'emico', Time, &
                             'Total Emission Rate of CO', 'kg m-2 s-1', &
               standard_name='tendency_of_atmosphere_mass_content_of_carbon_monoxide_due_to_emission')
   id_emiso2 = register_cmip_diag_field_2d ( module_name, 'emiso2', Time, &
                              'Total Emission Rate of SO2', 'kg m-2 s-1', &
                standard_name='tendency_of_atmosphere_mass_content_of_sulfur_dioxide_due_to_emission')
   id_eminh3 = register_cmip_diag_field_2d ( module_name, 'eminh3', Time, &
                              'Total Emission Rate of NH3', 'kg m-2 s-1', &
                standard_name='tendency_of_atmosphere_mass_content_of_ammonia_due_to_emission')
   id_eminox_woL = register_cmip_diag_field_2d ( module_name, 'eminox_woL', Time, &
                              'Total Emission Rate of NOx without lightning NOx', 'kg m-2 s-1', &
                standard_name='tendency_of_atmosphere_mass_content_of_nox_expressed_as_nitrogen_due_to_emission')
   id_emiisop_woB = register_cmip_diag_field_2d ( module_name, 'emiisop_woB', Time, &
                              'Total Emission Rate of Isoprene without biogenic', 'kg m-2 s-1', &
                standard_name='tendency_of_atmosphere_mass_content_of_isoprene_due_to_emission')

!--for Ox(jmao,1/1/2011)
   id_prodox = register_diag_field( module_name, 'Ox_prod', axes(1:3), &
        Time, 'Ox_prod','VMR/s')
   id_lossox = register_diag_field( module_name, 'Ox_loss', axes(1:3), &
        Time, 'Ox_loss','VMR/s')
!f1p
   id_pso4_h2o2   = register_diag_field( module_name, 'PSO4_H2O2',axes(1:3), Time, 'PSO4_H2O2','mole/m2/s')
   id_pso4_o3     = register_diag_field( module_name, 'PSO4_O3',axes(1:3), Time, 'PSO4_O3','mole/m2/s')

   id_gso2        = register_diag_field( module_name, 'gamma_so2',axes(1:3), Time, 'gamma_so2','unitless')
   id_ghno3_d     = register_diag_field( module_name, 'gamma_hno3_d',axes(1:3), Time, 'gamma_hno3_d','unitless')



   n_hno3d = 0
   n_so4d  = 0

   do i=1,n_dust_tracers
      if (dust_tracers(i)%is_hno3d) then
         n_hno3d = n_hno3d+1
         id_phno3_d(n_hno3d)    = register_diag_field( module_name, 'P'//trim(dust_tracers(i)%name),axes(1:3), &
              Time,  'P'//trim(dust_tracers(i)%name),'mole/m2/s')         

         write(fld,'(A6,I1.1,9X)') 'hno3_d',n_hno3d
         usr_hno3_dust(n_hno3d) = get_rxt_ndx(trim(fld)) - phtcnt
         id_rx_hno3_dust(n_hno3d) = register_diag_field( module_name, 'rx_'//TRIM(fld), axes(1:3), Time, 'rx_'//TRIM(fld),'1/s')

         write(fld,'(A6,I1.1,9X)') 'n2o5_d',n_hno3d 
         usr_n2o5_dust(n_hno3d) = get_rxt_ndx(trim(fld)) - phtcnt
         id_rx_n2o5_dust(n_hno3d) = register_diag_field( module_name, 'rx_'//TRIM(fld), axes(1:3), Time, 'rx_'//TRIM(fld),'1/s')

         write(fld,'(A5,I1.1,10X)') 'no3_d',n_hno3d
         usr_no3_dust(n_hno3d) = get_rxt_ndx(trim(fld)) - phtcnt
         id_rx_no3_dust(n_hno3d) = register_diag_field( module_name, 'rx_'//TRIM(fld), axes(1:3), Time, 'rx_'//TRIM(fld),'1/s')

      end if
      if (dust_tracers(i)%is_so4d) then
         n_so4d = n_so4d+1
         id_pso4_d(n_so4d)    = register_diag_field( module_name, 'P'//trim(dust_tracers(i)%name),axes(1:3), &
              Time,  'P'//trim(dust_tracers(i)%name),'mole/m2/s')         

         write(fld,'(A5,I1.1,10X)') 'so4_d',n_so4d
         usr_so4_dust(n_so4d) = get_rxt_ndx(trim(fld)) - phtcnt
         id_rx_so4_dust(n_so4d) = register_diag_field( module_name, 'rx_'//TRIM(fld), axes(1:3), Time, 'rx_'//TRIM(fld),'1/s')

         write(fld,'(A5,I1.1,10X)') 'so2_d',n_so4d
         usr_so2_dust(n_so4d) = get_rxt_ndx(trim(fld)) - phtcnt
         id_rx_so2_dust(n_so4d) = register_diag_field( module_name, 'rx_'//TRIM(fld), axes(1:3), Time, 'rx_'//TRIM(fld),'1/s')

      end if
   end do

   if (mpp_root_pe().eq.mpp_pe()) then
      write(*,*) 'usr_n2o5_dust', usr_n2o5_dust,id_rx_n2o5_dust
      write(*,*) 'usr_hno3_dust', usr_hno3_dust,id_rx_hno3_dust
      write(*,*) 'usr_so4_dust', usr_so4_dust,id_rx_so4_dust
      write(*,*) 'usr_so2_dust', usr_so2_dust,id_rx_so2_dust
      write(*,*) 'usr_no3_dust', usr_no3_dust,id_rx_no3_dust
   end if

   

   id_phno3_g_d    = register_diag_field( module_name, 'PHNO3_G_D',axes(1:3), &
        Time, 'PHNO3_G_D','mole/m2/s')

   if (max_dust<n_hno3d .or. max_dust<n_so4d) then
      call error_mesg ('tropchem_driver_init', &
           'max_dust<n_hno3d .or. max_dust<n_so4d', FATAL)
   end if


   id_sa_aerosol  = register_diag_field( module_name, 'sa_aerosol',axes(1:3), Time, 'sa_aerosol','cm2/cm3')
   id_sa_so4  = register_diag_field( module_name, 'sa_so4',axes(1:3), Time, 'sa_so4','cm2/cm3')
   id_sa_bc  = register_diag_field( module_name, 'sa_bc',axes(1:3), Time, 'sa_bc','cm2/cm3')
   id_sa_oa  = register_diag_field( module_name, 'sa_oa',axes(1:3), Time, 'sa_oa','cm2/cm3')
   id_sa_ss  = register_diag_field( module_name, 'sa_ss',axes(1:3), Time, 'sa_ss','cm2/cm3')
   id_sa_dust  = register_diag_field( module_name, 'sa_dust',axes(1:3), Time, 'sa_dust','cm2/cm3')




   id_aerosol_pH  = register_diag_field( module_name, 'aerosol_pH',axes(1:3), Time, 'aerosol_ph','unitless', mask_variant = .true.,missing_value=missing_value)
   id_cloud_pH  = register_diag_field( module_name, 'cloud_pH',axes(1:3), Time, 'cloud_ph','unitless', mask_variant = .true.,missing_value=missing_value)
   id_cloud_pHw  = register_diag_field( module_name, 'cloud_pHw',axes(1:3), Time, 'cloud_ph weighted by cloud fraction','unitless', mask_variant = .true.,missing_value=missing_value)
   id_cloud_pHwl  = register_diag_field( module_name, 'cloud_pHwl',axes(1:3), Time, 'cloud_ph weighted by cloud water','unitless', mask_variant = .true.,missing_value=missing_value)
   id_cld_amt_chem  = register_diag_field( module_name, 'cld_amt_chem',axes(1:3), Time, 'cloud fraction for chemistry','unitless', mask_variant = .true.,missing_value=missing_value)
   id_cld_liq_chem  = register_diag_field( module_name, 'cld_liq_chem',axes(1:3), Time, 'cloud liquid water for chemistry','unitless', mask_variant = .true.,missing_value=missing_value)

   if (id_cloud_pHw .gt. 0 .and. id_cld_amt_chem.le.0) then
      call error_mesg ('tropchem_driver_init', &
           'cld_amt_chem needs to be archived if cloud_pHw is requested', FATAL)
   end if
   if (id_cloud_pHwl .gt. 0 .and. id_cld_liq_chem.le.0) then
      call error_mesg ('tropchem_driver_init', &
           'cld_liq_chem needs to be archived if cloud_pHwl is requested', FATAL)
   end if


!for cmip6
   ID_pso4_aq_kg_m2_s      = register_cmip_diag_field_3d (  module_name,'pso4_aq_kg_m2_s', Time, &
                             'Aqueous-phase Production Rate of Sulfate Aerosol', 'kg m-2 s-1',  &
      standard_name='tendency_of_atmosphere_mass_content_of_sulfate_dry_aerosol_particles_due_to_aqueous_phase_net_chemical_production')
   ID_pso4_gas_kg_m2_s     = register_cmip_diag_field_3d (  module_name,'pso4_gas_kg_m2_s', Time, &
                            'Gas-phase Production Rate of Sulfate Aerosol', 'kg m-2 s-1',  &
      standard_name='tendency_of_atmosphere_mass_content_of_sulfate_dry_aerosol_particles_due_to_gaseous_phase_net_chemical_production')


   ID_lossch4      = register_cmip_diag_field_3d (  module_name,'lossch4', Time, &
                     'Monthly Loss of atmospheric Methane', 'mol m-3 s-1',  &
      standard_name='tendency_of_atmosphere_mole_concentration_of_methane_due_to_chemical_destruction')
   ID_lossco       = register_cmip_diag_field_3d (  module_name,'lossco', Time, &
                     'Monthly Loss of atmospheric Carbon Monoxide', 'mol m-3 s-1',  &
      standard_name='tendency_of_atmosphere_mole_concentration_of_carbon_monoxide_due_to_chemical_destruction')
   ID_lossn2o      = register_cmip_diag_field_3d (  module_name,'lossn2o', Time, &
                     'Monthly Loss of atmospheric Nitrous Oxide', 'mol m-3 s-1',  &
      standard_name='tendency_of_atmosphere_mole_concentration_of_nitrous_oxide_due_to_chemical_destruction')
   ID_o3loss       = register_cmip_diag_field_3d (  module_name,'o3loss', Time, &
                     'O3 destruction rate', 'mol m-3 s-1',  &
      standard_name='tendency_of_atmosphere_mole_concentration_of_ozone_due_to_chemical_destruction')
   ID_o3prod       = register_cmip_diag_field_3d (  module_name,'o3prod', Time, &
                     'O3 production rate', 'mol m-3 s-1',  &
      standard_name='tendency_of_atmosphere_mole_concentration_of_ozone_due_to_chemical_production')

   do i=1,pcnstm1
      id_chem_tend(i) = register_diag_field( module_name, trim(tracnam(i))//'_chem_dt', axes(1:3), &
                                             Time, trim(tracnam(i))//'_chem_dt','VMR/s' )
      id_prod(i) = register_diag_field( module_name, trim(tracnam(i))//'_prod', axes(1:3), &
                                        Time, trim(tracnam(i))//'_prod','VMR/s')
      id_loss(i) = register_diag_field( module_name, trim(tracnam(i))//'_loss', axes(1:3), &
                                        Time, trim(tracnam(i))//'_loss','VMR/s')
      id_prod_mol(i) = register_diag_field( module_name, trim(tracnam(i))//'_prodm', axes(1:3), &
                                        Time, trim(tracnam(i))//'_prodm','mole/m2/s')
      id_loss_mol(i) = register_diag_field( module_name, trim(tracnam(i))//'_lossm', axes(1:3), &
                                        Time, trim(tracnam(i))//'_lossm','mole/m2/s')
      if( has_emis(i) ) then
         id_emis(i) = register_diag_field( module_name, trim(tracnam(i))//'_emis', axes(1:2), &
                                           Time, trim(tracnam(i))//'_emis', 'molec/cm2/s')
      else
         id_emis(i) = 0
      end if
      if( has_emis3d(i) ) then
         id_emis3d(i) = register_diag_field( module_name, trim(tracnam(i))//'_emis3d', axes(1:3), &
                                             Time, trim(tracnam(i))//'_emis3d', 'molec/cm2/s')
      else
         id_emis3d(i) = 0
      end if

      if( has_ubc(i) ) then
         id_ub(i) = register_diag_field( module_name, trim(tracnam(i))//'_up', axes(1:3), &
                                         Time, trim(tracnam(i))//'_up','VMR' )
      else
         id_ub(i) = 0
      end if
      if( has_lbc(i) ) then
         if ( lbc_dry(i) ) then
            fld=' dry'
         else
            fld=''
         end if
         if (has_lbc_2d(i)) then
            id_lb(i) = register_diag_field( module_name, trim(tracnam(i))//'_lbc', axes(1:2), &
                 Time, trim(tracnam(i))//'_lbc','VMR'//trim(fld) )            
         else
            id_lb(i) = register_diag_field( module_name, trim(tracnam(i))//'_lbc', &
                 Time, trim(tracnam(i))//'_lbc','VMR'//trim(fld) )
         end if
      else
         id_lb(i) = 0
      end if
      if( has_airc(i) ) then
         id_airc(i) = register_diag_field( module_name, trim(tracnam(i))//'_airc_emis', axes(1:3), &
                                           Time, trim(tracnam(i))//'_airc_emis','molec/cm2/s' )
      else
         id_airc(i) = 0
      end if
   end do

!-----------------------------------------------------------------------
!     ... Register diagnostic fields for photolysis rates
!-----------------------------------------------------------------------
   do i=1,phtcnt
      write(fld,'(''jval_'',I3.3,8x)') i
      id_jval(i) = register_diag_field( module_name, TRIM(fld), axes(1:3), Time, TRIM(fld),'1/s')
   end do
   ID_jno2 = register_cmip_diag_field_3d (  module_name,'jno2', Time, &
                'Photolysis Rate of NO2', 's-1',  &
                standard_name='photolysis_rate_of_nitrogen_dioxide')
   ID_jo1d = register_cmip_diag_field_3d (  module_name,'photo1d', Time, &
                'Photolysis Rate of Ozone (O3) to Excited Atomic Oxygen (singlet d: O1d)', 's-1',  &
                standard_name='photolysis_rate_of_ozone_to_1D_oxygen_atom')

!-----------------------------------------------------------------------
!     ... Register diagnostic fields for kinetic rate constants
!-----------------------------------------------------------------------
   do i=1,gascnt
      write(fld,'(''k_rxn'',I3.3,8x)') i
      id_rate_const(i) = register_diag_field( module_name, TRIM(fld), axes(1:3), Time, TRIM(fld),'cm3/molec/s')
   end do

!-----------------------------------------------------------------------
!     ... initialize time_interp
!-----------------------------------------------------------------------
   call time_interp_init


!-----------------------------------------------------------------------
!     ... initialize mpp clock id
!-----------------------------------------------------------------------
   clock_id = mpp_clock_id('Chemistry')
   call setsox_init(trop_option)
   call chemdr_init(trop_option)

!initialize diag array
   if ( trop_option%aerosol_thermo == AERO_ISORROPIA ) then
      small_value = 1.e-25
   else
      small_value = 1.e-20
   end if

   call tropchem_types_init(trop_diag,small_value)
   if ( query_cmip_diag_id(ID_pso4_aq_kg_m2_s) .or. id_pso4_h2o2 > 0 ) then
      trop_diag%nb_diag       = trop_diag%nb_diag + 1
      trop_diag%ind_pso4_h2o2 = trop_diag%nb_diag
   end if
   if ( query_cmip_diag_id(ID_pso4_aq_kg_m2_s) .or. id_pso4_o3 > 0 ) then
      trop_diag%nb_diag       = trop_diag%nb_diag + 1
      trop_diag%ind_pso4_o3   = trop_diag%nb_diag
   end if
   if ( id_phno3_g_d > 0 ) then
      trop_diag%nb_diag          = trop_diag%nb_diag + 1
      trop_diag%ind_phno3_g_d    = trop_diag%nb_diag
   end if
   if ( id_pso4_g_d > 0 ) then
      trop_diag%nb_diag          = trop_diag%nb_diag + 1
      trop_diag%ind_pso4_g_d     = trop_diag%nb_diag
   end if
   if ( id_ghno3_d > 0 ) then
      trop_diag%nb_diag          = trop_diag%nb_diag + 1
      trop_diag%ind_ghno3_d      = trop_diag%nb_diag
   end if
   do i=1,5
      if ( id_phno3_d(i) > 0 ) then
         trop_diag%nb_diag          = trop_diag%nb_diag + 1
         trop_diag%ind_phno3_d(i)   = trop_diag%nb_diag
      end if
   end do
   do i=1,5
      if ( id_pso4_d(i) > 0 ) then
         trop_diag%nb_diag          = trop_diag%nb_diag + 1
         trop_diag%ind_pso4_d(i)   = trop_diag%nb_diag
      end if
   end do
   if ( id_gso2 > 0 ) then
      trop_diag%nb_diag       = trop_diag%nb_diag + 1
      trop_diag%ind_gso2      = trop_diag%nb_diag
   end if
   if ( id_sa_aerosol > 0 ) then
      trop_diag%nb_diag       = trop_diag%nb_diag + 1
      trop_diag%ind_SA_aerosol= trop_diag%nb_diag
   end if
   if ( id_sa_so4 > 0 ) then
      trop_diag%nb_diag       = trop_diag%nb_diag + 1
      trop_diag%ind_SA_SO4    = trop_diag%nb_diag
   end if
   if ( id_sa_bc > 0 ) then
      trop_diag%nb_diag       = trop_diag%nb_diag + 1
      trop_diag%ind_SA_BC     = trop_diag%nb_diag
   end if
   if ( id_sa_oa > 0 ) then
      trop_diag%nb_diag       = trop_diag%nb_diag + 1
      trop_diag%ind_SA_OA     = trop_diag%nb_diag
   end if
   if ( id_sa_ss > 0 ) then
      trop_diag%nb_diag       = trop_diag%nb_diag + 1
      trop_diag%ind_SA_SS     = trop_diag%nb_diag
   end if
   if ( id_sa_dust > 0 ) then
      trop_diag%nb_diag       = trop_diag%nb_diag + 1
      trop_diag%ind_SA_DUST   = trop_diag%nb_diag
   end if
   if ( id_aerosol_pH > 0 ) then
      trop_diag%nb_diag           = trop_diag%nb_diag + 1
      trop_diag%ind_aerosol_pH    = trop_diag%nb_diag
   end if
   if ( id_cloud_pH > 0 .or. id_cloud_pHw > 0) then
      trop_diag%nb_diag           = trop_diag%nb_diag + 1
      trop_diag%ind_cloud_pH      = trop_diag%nb_diag
   end if

   module_is_initialized = .true.


!-----------------------------------------------------------------------

end function tropchem_driver_init
!</FUNCTION>

!#####################################################################

subroutine tropchem_driver_time_vary (Time)

      type(time_type), intent(in) :: Time
      type(time_type) :: lbc_time

      integer :: yr, mo,day, hr,min, sec, mo_yr, dum, dayspmn
      integer :: n


!-----------------------------------------------------------------------
!     ... initialize number of shortwave bands
!-----------------------------------------------------------------------
      if (num_solar_bands == 0) then
         call shortwave_number_of_bands (num_solar_bands)
      endif


      do n=1, size(inter_emis,1)
        if (has_emis(n)) then
          call obtain_interpolator_time_slices (inter_emis(n), Time)
        endif
      end do

      do n=1, size(inter_emis3d,1)
        if (has_emis3d(n)) then
          call obtain_interpolator_time_slices (inter_emis3d(n), Time)
        endif
      end do

      do n=1, size(inter_aircraft_emis,1)
        if (has_airc(n)) then
          call obtain_interpolator_time_slices   &
                                         (inter_aircraft_emis(n), Time)
        endif
      end do

      call obtain_interpolator_time_slices (conc, Time)

      call obtain_interpolator_time_slices (sulfate, Time)

      do n=1, size(ub,1)
        if (has_ubc(n)) then
          call obtain_interpolator_time_slices (ub(n), Time)
        endif
      end do

      do n=1,size(has_lbc,1)
         if (has_lbc_2d(n)) then
            if (fixed_lbc_time(n)) then
               lbc_Time = lbc_entry(n)
            else
               lbc_Time = Time
            end if            
            call obtain_interpolator_time_slices (lbc_interp(n), lbc_time)
         end if
      end do
     

      call strat_chem_dcly_dt_time_vary (Time)

!----------------------------------------------------------------------
!    determine if this time step starts a new month; if so, then
!    new interactive isoprene emission data is needed, and the
!    necessary flag is set.
!----------------------------------------------------------------------
      do n=1,pcnstm1
        if ( has_xactive_emis(n) .or. id_xactive_emis(n)>0 ) then
          if (tracnam(n) .eq. "DMS") then
              call atmos_sulfate_time_vary (Time)
              exit
          endif
        endif
      end do

end subroutine tropchem_driver_time_vary




!#####################################################################

subroutine tropchem_driver_endts


      integer :: n

      do n=1, size(inter_emis,1)
        if (has_emis(n)) then
          call unset_interpolator_time_flag(inter_emis(n))
         endif
      end do

      do n=1, size(inter_emis3d,1)
        if (has_emis3d(n)) then
         call unset_interpolator_time_flag(inter_emis3d(n))
        endif
      end do

      do n=1, size(inter_aircraft_emis,1)
        if (has_airc(n)) then
          call unset_interpolator_time_flag(inter_aircraft_emis(n))
        endif
      end do

      call unset_interpolator_time_flag(conc)
      call unset_interpolator_time_flag(sulfate)

      do n=1, size(ub,1)
        if (has_ubc(n)) then
          call unset_interpolator_time_flag(ub(n))
        endif
      end do

      do n=1,size(has_lbc_2d,1)
         if (has_lbc_2d(n)) then
            call unset_interpolator_time_flag(lbc_interp(n))
         end if
      end do

      call strat_chem_dcly_dt_endts



end subroutine tropchem_driver_endts


!######################################################################

subroutine tropchem_driver_end

!-----------------------------------------------------------------------
!     ... initialize mpp clock id
!-----------------------------------------------------------------------
   
   deallocate(nb_N_Ox)
   module_is_initialized = .false.
   
    
!-----------------------------------------------------------------------

end subroutine tropchem_driver_end

!#######################################################################

! <SUBROUTINE NAME="read_2D_emis_data">
!   <OVERVIEW>
!     Read emissions file
!   </OVERVIEW>
!   <DESCRIPTION>
!     Reads tracer surface emissions from a NetCDF file
!   </DESCRIPTION>
!   <TEMPLATE>
!     call read_2D_emis_data( emis_type, emis, Time, &
!                             field_names, &
!                             Ldiurnal, coszen, half_day, lon, &
!                             is, js, id_emis_diag )
!   </TEMPLATE>

subroutine read_2D_emis_data( emis_type, emis, Time, Time_next, &
                              field_names, &
                              Ldiurnal, coszen, half_day, lon, &
                              is, js, skip_biogenic_emis, scale_emis, skip_field )

   type(interpolate_type),intent(inout) :: emis_type
   real, dimension(:,:),intent(out) :: emis
   type(time_type),intent(in) :: Time, Time_next
   character(len=*),dimension(:), intent(in) :: field_names
   real,dimension(:), intent(in) :: scale_emis
   character(len=*), intent(in), optional :: skip_field
   logical, intent(in) :: Ldiurnal
   real, dimension(:,:), intent(in) :: coszen, half_day, lon
   integer, intent(in) :: is, js
   logical, intent(in) :: skip_biogenic_emis

   integer :: i, j, k
   logical :: used
   real, dimension(size(emis,1),size(emis,2)) :: temp_data
   real :: diurnal_scale_factor, gmt, iso_on, iso_off, dayfrac
   real :: local_angle, factor_tmp
   integer :: n
   
   emis(:,:) = 0.
   temp_data(:,:) = 0.
   do k = 1,size(field_names)
      if (skip_biogenic_emis .and. field_names(k) .eq. "biogenic") then
          temp_data(:,:) = 0.
      else
         if (present(skip_field) .and. trim(field_names(k)).eq.trim(skip_field)) then
            temp_data(:,:) = 0.                      
         else
            call interpolator(emis_type,Time,temp_data,field_names(k),is,js)
         end if
         temp_data = temp_data*scale_emis(k)
      end if
      emis(:,:) = emis(:,:) + temp_data(:,:)
   end do

   if (Ldiurnal) then
      do j=1,size(emis,2)
      do i=1,size(emis,1)
         if( coszen(i,j) < 0. ) then
            diurnal_scale_factor = 0.
         else
            iso_off = .8 * half_day(i,j)
            iso_on  = -iso_off
            dayfrac = iso_off/PI
            gmt = universal_time(Time)
            local_angle = gmt + lon(i,j) - PI
            if (local_angle >= PI) local_angle = local_angle - twopi
            if (local_angle < -PI) local_angle = local_angle + twopi
            if( local_angle >= iso_off .or. local_angle <= iso_on ) then
               diurnal_scale_factor = 0.
            else
               factor_tmp = local_angle - iso_on
               factor_tmp = factor_tmp / MAX(2.*iso_off,1.e-6)
               diurnal_scale_factor = 2. / dayfrac * (sin(PI*factor_tmp))**2
            end if
         end if
         emis(i,j) = emis(i,j) * diurnal_scale_factor
      end do
      end do
   end if

end subroutine read_2D_emis_data
!</SUBROUTINE>

!#######################################################################

! <SUBROUTINE NAME="read_3D_emis_data">
!   <OVERVIEW>
!     Read emissions file
!   </OVERVIEW>
!   <DESCRIPTION>
!     Reads tracer 3-D emissions from a NetCDF file
!   </DESCRIPTION>
!   <TEMPLATE>
!     call read_3D_emis_data( emis_type, emis, Time, phalf, &
!                             field_names, &
!                             Ldiurnal, coszen, half_day, lon, &
!                             is, js, id_emis_diag )
!   </TEMPLATE>

subroutine read_3D_emis_data( emis_type, emis, Time, Time_next, phalf, &
                              field_names, &
                              Ldiurnal, coszen, half_day, lon, &
                              is, js, id_emis_diag, scale_emis )

   type(interpolate_type),intent(inout) :: emis_type
   real, dimension(:,:,:),intent(in) :: phalf
   real, dimension(:,:,:),intent(out) :: emis
   type(time_type),intent(in) :: Time, Time_next
   character(len=*),dimension(:), intent(in) :: field_names
   real,dimension(:), intent(in) :: scale_emis
   logical, intent(in) :: Ldiurnal
   real, dimension(:,:), intent(in) :: coszen, half_day, lon
   integer, intent(in) :: is, js
   integer, intent(in),optional :: id_emis_diag ! id for diagnostic


   integer :: i, j, k
   logical :: used
   real, dimension(size(emis,1),size(emis,2),size(emis,3)) :: temp_data
   real :: diurnal_scale_factor, gmt, iso_on, iso_off, dayfrac
   real :: local_angle, factor_tmp

   emis(:,:,:) = 0.
   temp_data(:,:,:) = 0.
   do k = 1,size(field_names)
      call interpolator(emis_type,Time,phalf,temp_data,field_names(k),is,js)
      emis(:,:,:) = emis(:,:,:) + temp_data(:,:,:)*scale_emis(k)
   end do
   if (Ldiurnal) then
      do j=1,size(emis,2)
      do i=1,size(emis,1)
         if( coszen(i,j) < 0. ) then
            diurnal_scale_factor = 0.
         else
            iso_off = .8 * half_day(i,j)
            iso_on  = -iso_off
            dayfrac = iso_off/PI
            gmt = universal_time(Time)
            local_angle = gmt + lon(i,j) + PI
            if (local_angle >= PI) local_angle = local_angle - twopi
            if (local_angle < -PI) local_angle = local_angle + twopi
            if( local_angle >= iso_off .or. local_angle <= iso_on ) then
               diurnal_scale_factor = 0.
            else
               factor_tmp = local_angle - iso_on
               factor_tmp = factor_tmp / MAX(2.*iso_off,1.e-6)
               diurnal_scale_factor = 2. / dayfrac * (sin(PI*factor_tmp))**2
            end if
         end if
         emis(i,j,:) = emis(i,j,:) * diurnal_scale_factor
      end do
      end do
   end if

   if (present(id_emis_diag)) then
      if (id_emis_diag > 0) then
         used = send_data(id_emis_diag,emis,Time_next,is_in=is,js_in=js)
      end if
   end if
end subroutine read_3D_emis_data
!</SUBROUTINE>

!#######################################################################

! <SUBROUTINE NAME="calc_xactive_emis">
!   <OVERVIEW>
!     Calculates interactive emissions
!   </OVERVIEW>
!   <DESCRIPTION>
!     Calculates interactive emissions
!   </DESCRIPTION>
!   <TEMPLATE>
!     call calc_xactive_emis( index, emis, Time, is, js, id_emis_diag )
!   </TEMPLATE>

subroutine calc_xactive_emis( index, Time, Time_next, lon, lat, pwt, is, ie, js, je, &
                              area, land, ocn_flx_fraction, tsurf, w10m, emis, &
                              kbot, id_emis_diag )

   integer,intent(in) :: index
   type(time_type),intent(in) :: Time, Time_next
   real, intent(in), dimension(:,:) :: lon, lat
   real, intent(in), dimension(:,:,:) :: pwt
   integer, intent(in) :: is, ie, js, je
   real, intent(in), dimension(:,:) :: area    ! grid box area (m^2)
   real, intent(in), dimension(:,:) :: land    ! land fraction
   real, intent(in), dimension(:,:) :: ocn_flx_fraction
   real, intent(in), dimension(:,:) :: tsurf   ! surface temperature (K)
   real, intent(in), dimension(:,:) :: w10m    ! wind speed at 10m (m/s)
   real, dimension(:,:,:),intent(out) :: emis  ! VMR/s
   integer, intent(in), dimension(:,:), optional :: kbot
   integer, intent(in),optional :: id_emis_diag ! id for diagnostic

   logical :: used


   if (index == dms_ndx) then
      call atmos_DMS_emission( lon, lat, area, ocn_flx_fraction, tsurf, w10m, pwt, &
                               emis, Time, Time_next, is, ie, js, je, kbot )
   else
      call error_mesg ('calc_xactive_emis', &
                       'Interactive emissions not defined for species: '//trim(tracnam(index)), FATAL)
   end if

   if (present(id_emis_diag)) then
      if (id_emis_diag > 0) then
         used = send_data( id_emis_diag, emis, Time_next, is_in=is, js_in=js)
      end if
   end if
end subroutine calc_xactive_emis
!</SUBROUTINE>

!#######################################################################

! <SUBROUTINE NAME="init_emis_data">
!   <OVERVIEW>
!     Open emissions file
!   </OVERVIEW>
!   <DESCRIPTION>
!     Opens NetCDF file of tracer surface emissions for reading,
!     and set up interpolation to model grid/time
!   </DESCRIPTION>
!   <TEMPLATE>
!     call init_emis_data( emis_type, model, method_type, index, file_name, &
!                          lonb_mod, latb_mod, field_type, flag, diurnal )
!   </TEMPLATE>

subroutine init_emis_data( emis_type, model, method_type, pos, file_name, &
                           lonb_mod, latb_mod, field_type, flag, diurnal, &
                           axes, Time, land_does_emis )

   type(interpolate_type),intent(inout) :: emis_type
   integer, intent(in) :: model,pos
   character(len=*),intent(in) :: method_type
   character(len=*),intent(inout) ::file_name
   real,intent(in),dimension(:,:) :: lonb_mod,latb_mod
   type(field_init_type),intent(out) :: field_type
   logical, intent(out) :: flag, diurnal
   logical, intent(out), optional :: land_does_emis
   integer        , intent(in)  :: axes(4)
   type(time_type), intent(in)  :: Time

   character(len=128) :: name, control
   integer :: nfields
   integer :: flag_name, flag_file, flag_diurnal, flag_scale
   character(len=64) :: emis_name, emis_file, control_diurnal

   integer :: n,n2
   real    :: scale_emis

   flag = .false.
   diurnal = .false.
   control = ''
   if( query_method(trim(method_type),model,pos,name,control) ) then
      if( trim(name(1:4)) == 'file' ) then
         flag = .true.
         flag_file = parse(control, 'file', emis_file)
         flag_name = parse(control, 'name', emis_name)
         flag_diurnal = parse(control, 'diurnal', control_diurnal)
         if(flag_file > 0) then
            file_name = emis_file
         else if (flag_name > 0) then
            select case (trim(method_type))
               case ('emissions3d')
                  file_name  = trim(file_emis3d_1)//trim(emis_name)//trim(file_emis3d_2)
               case default
                  file_name  = trim(file_emis_1)//trim(emis_name)//trim(file_emis_2)
            end select
         end if
         diurnal = (flag_diurnal > 0)

         call interpolator_init( emis_type, trim(file_name), &
                                 lonb_mod, latb_mod,  &
                                 data_out_of_bounds=(/CONSTANT/), &
                                 vert_interp=(/INTERP_WEIGHTED_P/) )
         call query_interpolator(emis_type,nfields=nfields)
         allocate(field_type%field_names(nfields))
         allocate(field_type%scale_emis(nfields))         
         call query_interpolator(emis_type,field_names=field_type%field_names)
         do n=1,nfields
            field_type%scale_emis(n) = 1.

            flag_scale = parse(control, field_type%field_names(n), scale_emis)
            if (flag_scale > 0) then
               field_type%scale_emis(n) = scale_emis
            else
               do n2=1,max_scale_emis_fields
                  if (trim(field_type%field_names(n)).eq.trim(scale_emis_field_names(n2))) then
                     field_type%scale_emis(n) = scale_emis_field_values(n2)
                  end if
               end do
            end if

            if (mpp_root_pe().eq.mpp_pe()) write(*,40) field_type%field_names(n), &
                 field_type%scale_emis(n)
            
         end do         
      end if
      if ( present(land_does_emis) )  land_does_emis  = (index(lowercase(name),'land:lm3')>0)
   end if

40 format (A,F5.3)
end subroutine init_emis_data
!</SUBROUTINE>

!#######################################################################

! <SUBROUTINE NAME="init_xactive_emis">
!   <OVERVIEW>
!     Set up interactive emission calculations
!   </OVERVIEW>
!   <DESCRIPTION>
!     Set up interactive emission calculations
!   </DESCRIPTION>
!   <TEMPLATE>
!     call init_xactive_emis( model, method_type, index, species, &
!                             axes, Time, lonb_mod, latb_mod, phalf, &
!                             flag, mask )
!   </TEMPLATE>

subroutine init_xactive_emis( model, method_type, index, species, &
                              axes, Time, lonb_mod, latb_mod, phalf, &
                              flag, id_xemis, mask )

   integer,         intent(in)  :: model, index
   character(len=*),intent(in)  :: method_type, species
   integer        , intent(in)  :: axes(4)
   type(time_type), intent(in)  :: Time
   real,            intent(in), dimension(:,:)   :: lonb_mod,latb_mod
   real,            intent(in), dimension(:,:,:) :: phalf
   real,            intent(in), dimension(:,:,:), optional :: mask
   logical,         intent(out) :: flag
   integer,         intent(out) :: id_xemis

   character(len=64) :: name, control
   integer :: nhalf, nfull

   flag = .false.
   control = ''
   nhalf = SIZE(phalf,3)
   nfull = nhalf - 1

   flag = query_method(trim(method_type),model,index,name,control)

   if (flag) then
   if (trim(species) .eq. "DMS") then
            id_xemis = &
               register_diag_field( module_name, trim(species)//'_xactive_emis', axes(1:3), &
                                    Time, trim(species)//'_xactive_emis', 'VMR/s')
            if (flag .or. id_xemis>0) then
               call atmos_sulfate_init( lonb_mod, latb_mod, nfull, axes, Time, mask )
            end if
   endif
   end if

end subroutine init_xactive_emis
!</SUBROUTINE>

!> @brief This function is just a wrapper for Fortran's `open`
!! @return Unique unit number
function tropchem_open_file (filename) result (funit)
  character(len=*), intent(in), optional :: filename
  integer :: funit

  open(file=filename, form='formatted',action='read', newunit=funit)

end function tropchem_open_file

!> @brief This function just checks if a tracer initialized in a set of files
!! @return flag indicating if a tracer initialized
function check_if_tracer_initialized(tracername) result (tracer_initialized)
  character(len=*), intent(in), optional :: tracername
  logical :: tracer_initialized

  type(FmsNetcdfDomainFile_t) :: fileobj !< fms2io fileobj for domain decomposed

  tracer_initialized = .false.

  if (open_file(fileobj, 'INPUT/atmos_tracers.res.nc', "read", tropchem_domain)) then
     tracer_initialized = variable_exists(fileobj, tracername)
     call close_file(fileobj)
  elseif (open_file(fileobj, 'INPUT/fv_tracer.res.nc', "read", tropchem_domain)) then
     tracer_initialized = variable_exists(fileobj, tracername)
     call close_file(fileobj)
  elseif (open_file(fileobj, 'INPUT/tracer_'//trim(lowercase(tracername))//'.res', "read", tropchem_domain)) then
     tracer_initialized = variable_exists(fileobj, tracername)
     call close_file(fileobj)
  endif

end function check_if_tracer_initialized

!############################################################################
end module tropchem_driver_mod
