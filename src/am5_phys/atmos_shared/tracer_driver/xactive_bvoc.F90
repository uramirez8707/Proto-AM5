module xactive_bvoc_mod
!
! <CONTACT EMAIL="Jordan.Schnell@noaa.gov">
!   Jordan L. Schnell
! </CONTACT>
!
! <OVERVIEW>
!   This code calculates interactive biogenic emissions
!   due to variations in environmental condtions
! </OVERVIEW
!
! <DESCRIPTION>
!
! This code calculates interactive biogenic emisisons
! (VOCs + NO + CO) due to variations in environmental conditions
! such as temperature and light. There are currently 2 options
! to choose how these are calculated. Their descriptions and
! the required input datasets are listed below.
!
!   *** All default input datasets for AM3 isoprene, MEGANv2 and
!       MEGANv3 are in the tarball:
!       /lustre/f1/unswept/Jordan.Schnell/input/xactive_emissions/megan.xactive.bvoc.tar
!
!
!   (0) All version are backward compatiable such that they can reproduce
!       AM3 interactive isoprene emissions. This is controlled by the namelist
!       variable: "do_AM3_ISOP". Version specific functions/subroutines have
!       the format *_AM3
!
!   (1) xactive_algorithm = 'MEGAN2'
!
!        a) All species (see below) emisions are calculated following the
!           MEGAN PECEEA algorithm of Guenther et al., 2006
!        b) Emission factor files must be included.
!        c) Version specific functions/subroutines have the format '***_megan2'
!
!
!        REQUIRED DATASETS
!        --------------------
!        LAI file         = mksrf_lai.060929.nc
!        PFT file         = mksrf_pft.060929.nc
!        Surf T file      = tas_monthly_clim_1980-2000.nc
!        sw down file     = dswrf_monthly_clim_1980-2000.nc
!        ISOP EM capacity = megan2.ISOP.nc
!        Other EM files   = megan2.$(SPECIES).nc (in megan2.xactive.bvoc.tar)
!
!       Updates from the original xactive algorithms which were
!       implemented in GFDL-AM3 by Arlene M. Fiore and Vaishali A. Naik,
!       include:
!
!       (A) Compounds beyond isoprene including:
!         **(1) Monoterpenes--------- C10H16 (and sequisterpenes, see below)
!           (2) Methanol------------- CH3OH
!           (3) Acetone-------------- CH3COCH3
!           (4) Carbon Monoxide------ CO
!           (5) Ethanol-------------- C2H5OH
!           (6) Acetylaldehyde------- CH3CHO
!           (7) Formaldehyde--------- CH2O
!           (8) Acetic Acid---------- CH3COOH
!           (9) Ethene--------------- C2H4
!          (10) Ethane--------------- C2H6
!          (11) Propene-------------- C3H6
!          (12) Propane-------------- C3H8
!
!        (B)  **Monoterpenes vs. Sequisterpenes
!          - Two namelist variables control how these are handled
!            (1) "do_SESQTERP"
!                   AND
!            (2) "do_PARSED_TERP"
!
!               (a) IF do_SESQTERP==TRUE && do_PARSED_TERP==TRUE THEN
!                   Emissions are calculated for each monoterpene AND
!                   sequisterpene species using their species specific
!                   MEGAN parameters. Currently, however they are all
!                   lumped as monoterpenes in the tracer array since
!                   sequisterpenes are not explicitly modeled (yet).
!
!               (b) IF do_SESQTERP==TRUE && do_PARSED_TERP==FALSE THEN
!                   Emisson factors for sequsiterpenes are included in the
!                   monoterpene emissions, but all monoterpenes and
!                   sequisterpenes use the alpha/beta pinene MEGAN parameters.
!
!               (c) IF do_SESQTERP==FALSE && do_PARSED_TERP==TRUE
!                   Same as (i) except sequisterpene emission factors are
!                   NOT included
!
!               (d) IF do_SESQTERP==FALSE && do_PARSED_TERP==TRUE
!                   Emission factors are summed over only MONOterpene species
!                   and the MEGAN parameters all follow alpha/beta pinene

!
!---------------------------------------------------------------------------------------------
!
!   (2) xactive_algorithm  = 'MEGAN3'
!
!       All species follow MEGAN v3.0 from Guenther et al., 2018
!       The code is in large part a duplication of the MEGAN3
!       source code located at https://sites.google.com/uci.edu/bai/megan/versions
!
!       At minimum, this option requires
!          (1) emisions capacities in main tarball
!          (2) the monthly mean surface temperature file from above
!          (3) LAIv (i.e., LAI divided by fractional vegetation cover...in tarball)
!
!       *** The MEGAN3 Emission factors are combined into one and are not pft specific
!       *** MEGAN3 uses LAIv, which is sum(LAI of all types) / vegetation cover fraction

!       (A) Depending on which online fields are used, and which gammas are specified,
!           additional input data may be required
!           (1) LAIv
!           (2) FCOVER (if do_ONLINE_LAI)
!           (3) W126 ozone
!           (4) CO2
!           (5) Soil moisture
!
!       (B) MEGAN3 has additional (Optional) gammas (i.e., envrionmental responses) that can be
!           applied that were not part of MEGAN2, these include:
!           (1) High temperatures
!           (2) Low temperatures
!           (3) High winds
!           (4) Air quality (W126 ozone)
!           (5) CO2 (isoprene only)
!           (6) Bi-directional LAI (only ethanol and acetylaldehyde)
!           (7) Soil mositure (not yet included)
!
!-----------------------------------------------------------------------------------------------
!
!   (3) xactive_algorithm = 'BEIS3'
!
!       NOT YET IMPLEMENTED
!
! </DESCRIPTION>


!------------------------------------------------------------------------------

use                mpp_mod,  only : input_nml_file, mpp_get_current_pelist
use                fms_mod,  only : write_version_number,  &
                                    mpp_pe,                &
                                    mpp_root_pe,           &
                                    mpp_npes,              &
                                    stdlog,                &
                                    check_nml_error,       &
                                    error_mesg,            &
                                    FATAL,                 &
                                    WARNING,               &
                                    NOTE
use        mpp_domains_mod,  only : domain2D, mpp_get_ntile_count
use            MO_GRID_MOD,  only : pcnstm1
use            fms2_io_mod,  only : FmsNetcdfFile_t, FmsNetcdfDomainFile_t, &
                                    register_restart_field, register_axis, unlimited, &
                                    open_file, read_restart, write_restart, close_file, &
                                    register_field, write_data, get_global_io_domain_indices, &
                                    register_variable_attribute, read_data, file_exists
use         M_TRACNAME_MOD,  only : tracnam
use      tracer_manager_mod, only : get_tracer_index,      &
                                    query_method
use      field_manager_mod,  only : MODEL_ATMOS
use       time_manager_mod,  only : time_type,             &
                                    get_date,              &
                                    set_date,              &
                                    time_type_to_real,     &
                                    operator(+),           &
                                    operator(-)
use           constants_mod, only : SECONDS_PER_DAY,       &
                                    DEG_TO_RAD,            &
                                    WTMAIR, rdgas,         &
                                    AVOGNO, PI, EPSLN
use         horiz_interp_mod, only: horiz_interp_type,     &
                                    horiz_interp_init,     &
                                    horiz_interp_new,      &
                                    horiz_interp,          &
                                    horiz_interp_del
use       diag_manager_mod, only  : send_data,             &
                                    register_diag_field,   &
                                    register_static_field, &
                                    get_base_time

!use something_in_the_land_mod, only : get_LAI,            &
!                                      get_SOIL_PARAM



implicit none

private

!------------------------------------------------------------------------------
!     ... INTERFACES ...
!------------------------------------------------------------------------------

public xactive_bvoc, xactive_bvoc_init, xactive_bvoc_end, &
       ind_xbvoc_ISOP, ind_xbvoc_TERP

!------------------------------------------------------------------------------
!      ... NAMELIST ...
!------------------------------------------------------------------------------


character(len=64)   :: xactive_algorithm = 'MEGAN2',                           & ! Algorithm used to calculate emisisons
                       file_LAI    = 'INPUT/mksrf_lai.060929.nc',              & ! filename: leaf area index (LAI)
                       file_PFT    = 'INPUT/mksrf_pft.060929.nc',              & ! filename: plant functional types (PFT)
                       file_PPFD   = 'INPUT/dswrf_monthly_clim_1980-2000.nc',  & ! filename: monthly avg sw down
                       file_TEMP   = 'INPUT/tas_monthly_clim_1980-2000.nc',    & ! filename: monthly avg surface T
                       file_FCOVER = 'INPUT/fcover_2008.nc',                   & ! filename: fractional vegetation cover
                       file_O3     = 'surface_ozone.nc',                       & ! filename: W2126 ozone
                       file_SOILM  = 'soil_moisture.nc',                       & ! filename: soil moisture
                       file_WS     = 'wind_speed.nc',                          & ! filename: for wind speed
                       file_CO2    = 'co2_conc.nc'                               ! filename: CO2 conc.

logical             :: do_AM3_ISOP    = .false., &             ! flag: Reproduce AM3 isoprene emissions?
                       do_SESQTERP    = .false., &             ! flag: compute sequisterpenes?
                       do_PARSED_TERP = .false., &             ! flag: parse terpenes (i.e., each species individually)?
                       do_GAMMA_BDLAI = .false., &             ! flag: Gamma for Bi-directional LAI?
                       do_GAMMA_CO2   = .false., &             ! flag: Gamma for CO2 effect?
                       do_GAMMA_AQ    = .false., &             ! flag: Gamma for Air Quality effects (i.e. ozone)?
                       do_GAMMA_SM    = .false., &             ! flag: Gamma for soil moisture?
                       do_GAMMA_HT    = .false., &             ! flag: Gamma for high temperature?
                       do_GAMMA_LT    = .false., &             ! flag: Gamma for low temperature?
                       do_GAMMA_HW    = .false.                ! flag: Gamma for wind storms?

logical             :: do_ONLINE_LAI     = .false.,  &         ! flag: online leaf area index (LAI)?
                       do_ONLINE_PFT     = .false.,  &         ! flag: online plant functional type (PFT)?
                       do_ONLINE_TEMP    = .TRUE.,   &         ! flag: online temperatures?
                       do_ONLINE_PPFD    = .TRUE.,   &         ! flag: online radiation?
                       do_ONLINE_O3      = .false.,  &         ! flag: online AQI/ W126 O3
                       do_ONLINE_WIND    = .false.,  &         ! flag: online wind speed
                       do_ONLINE_SM      = .false.,  &         ! flag: online soil mositure
                       do_ONLINE_CO2     = .false.             ! flag: online CO2 concentrations

logical             :: use_isop_shrub_crop_bug = .false.       ! if T, swap shrub and crop emission potentials and distributions
logical             :: fix_megan2_isop   = .TRUE.              ! if T, increases isop emis factors by 50% to achieve ~500 Tg/yr global

real                :: T_s = 297.                              ! Temperature that represents standard conditions [K]

real                :: scale_isoprene_emissions = 1.           ! Global scale for isoprene emission
real                :: scale_terpene_emissions = 1.            ! Global scale for terpene emission
real                :: min_land_frac = 0.01                    ! Fraction of land required to calculate emissions
integer             :: verbose = 3                             ! level of diagnostic output

integer, parameter  :: nPFT = 17, nVEG = 5                     ! Number of plant funct types, vegetation types (MEGAN2)
integer, parameter  :: nMOS = 12, nHOUR = 24                   ! Number of months,hours for STORE Arrays
real, parameter     :: twopi = 2.*PI
character(len=7), parameter :: module_name = 'tracers'
real                :: RHO_CANOPY = 1.                         ! Emisions lost in the canopy = (1 - RHO)

integer, parameter  :: ind_xbvoc_ISOP = 1, &                   ! Index for isoprene emissions in xbvoc4soa
                       ind_xbvoc_TERP = 2                      ! Index for terpene emissions in xbvoc4soa

namelist /xactive_bvoc_nml/                     &
                             xactive_algorithm, &
                             do_AM3_ISOP,       &
                             do_SESQTERP,       &
                             do_PARSED_TERP,    &
                             file_LAI,          &
                             file_PFT,          &
                             file_FCOVER,       &
                             file_PPFD,         &
                             file_TEMP,         &
                             file_O3,           &
                             file_CO2,          &
                             file_SOILM,        &
                             file_WS,           &
                             do_GAMMA_BDLAI,    &
                             do_GAMMA_CO2,      &
                             do_GAMMA_AQ,       &
                             do_GAMMA_SM,       &
                             do_GAMMA_HT,       &
                             do_GAMMA_LT,       &
                             do_GAMMA_HW,       &
                             do_ONLINE_LAI,     &
                             do_ONLINE_PFT,     &
                             do_ONLINE_TEMP,    &
                             do_ONLINE_PPFD,    &
                             do_ONLINE_O3,      &
                             do_ONLINE_CO2,     &
                             do_ONLINE_WIND,    &
                             do_ONLINE_SM,      &
                             RHO_CANOPY,        &
                             min_land_frac,     &
                             T_s,               &
                             scale_isoprene_emissions, &
                             scale_terpene_emissions,  &
                             use_isop_shrub_crop_bug,  &
                             fix_megan2_isop,   &
                             verbose

logical                     :: Ldebug = .false.
logical                     :: module_is_initialized = .false.
logical, dimension(pcnstm1) :: has_xactive_emis = .false.  ! Does the tracer have xactive emissions?

integer, dimension(pcnstm1) :: indices,     &
                               id_EMIS,     &
                               id_G_TEMP,   &
                               id_G_PAR,    &
                               id_G_AGE,    &
                               id_G_LAI,    &
                               id_G_BDLAI,  &
                               id_G_CO2,    &
                               id_G_AQ,     &
                               id_G_SM,     &
                               id_G_HT,     &
                               id_G_LT,     &
                               id_G_HW


real, allocatable, dimension(:,:)         :: MEGAN_PARAM      ! MEGAN MODEL PARAMETERS
real, allocatable, dimension(:,:)         :: TERP_PARAM       ! Parameters for parsed terpenes

real, allocatable, dimension(:,:,:)       :: ECISOP_AM3       ! Emisison capapcites for AM3 isoprene (MEGAN2)

real, allocatable, dimension(:,:,:,:)     :: ECBVOC           ! Emisison capacities (lat x lon x PFT/VEG x SPEC) MEGAN2
real, allocatable, dimension(:,:,:)       :: ECBVOC_MEGAN3    ! " " for MEGAN3 (combined for vegetation types)

real, allocatable, dimension(:,:,:,:)     :: ECTERP           ! Emissions capacities for parsed terpenes
real, allocatable, dimension(:,:,:)       :: ECTERP_MEGAN3    ! "" MEGAN3

real, allocatable, dimension(:,:,:)       :: LDFg             ! Gridded light dependent fractions (MEGAN3)
real, allocatable, dimension(:,:,:)       :: LDFg_TERP        ! same but for parsed terpenes (MEGAN3)


real, allocatable, dimension(:,:,:,:)     :: MLAI             ! monthly lai for each pft (could eventually tie to LM3)
real, allocatable, dimension(:,:,:)       :: MLAI_MEGAN3      ! LAIv data for MEGAN3 which isn't pft specific

real, allocatable, dimension(:,:,:)       :: PCTPFT           ! Percent coverage of each PFT
real, allocatable, dimension(:,:,:)       :: FCOVER           ! fractional cover of vegetation, used in MEGAN3 to
                                                              ! calc LAIv = LAI/FCOVER


real, allocatable, dimension(:,:,:)       :: T24_STORE,    &  ! Array to hold hourly temperature ? IS 'SAVE' necessary?
                                             P24_STORE,    &  ! '' PAR
                                             WS_STORE,     &  ! '' wind speed
                                             O3_STORE         ! '' surf O3

! Monthly mean sfc air T and sw down at surface, from C. Wiedinmyer 2/18/09 - Sheffield inputs (Princeton)
! CW provided 1948-2000; current input files take 1980-2000 average
real, allocatable, dimension(:,:,:)       :: Tmo, Pmo        ! Monthly mean temp and par

real, allocatable, dimension(:,:)         :: CO2_STORE,   &
                                             SOILM, WILT
type (domain2D), pointer                  :: xactive_domain !< Atmosphere domain

logical                                   :: in_different_file = .false.
integer                                   :: vers = 1


! Diagnostics
real, allocatable, dimension(:,:)         :: diag_gamma_temp, &
                                             diag_gamma_par,  &
                                             diag_gamma_ht,   &
                                             diag_gamma_lt,   &
                                             diag_gamma_hw,   &
                                             diag_gamma_aq,   &
                                             diag_gamma_co2,  &
                                             diag_gamma_sm
real, allocatable, dimension(:,:,:)       :: diag_gamma_age,  &
                                             diag_gamma_lai

real, allocatable, dimension(:,:)         :: diag_gamma_age_megan3, &
                                             diag_gamma_lai_megan3, &
                                             diag_gamma_bdlai_megan3

!---- version number ---------------------------------------------------
character(len=128), parameter :: version     = '$Id$'
character(len=128), parameter :: tagname     = '$Name$'
!-----------------------------------------------------------------------

contains

!#################################################################################
!
! <SUBROUTINE NAME="xactive_bvoc">
!   <OVERVIEW>
!     Calculates interactive bigenic VOC and CO emisisons
!   </OVERVIEW>
!   <DESCRIPTION>
!     This subroutine calculates interactive biogenic emissions from
!     vegetation in response to envrionmental conditions. The routine is
!     called from atmos_tracer_driver. First, the arrays that hold values
!     to calculate daily averages are updated. Then the species are looped
!     over, and depending on the namelist option, the appropriate subroutine
!     is called (e.g., reproduce AM3 isoprene, parsed terpenes, etc.). The
!     emissions are passed back as a tendency in the lowest model layer.
!   </DESCRIPTION>
!   <TEMPLATE>
!      call xactive_bvoc(lon, lat, land, is, ie, js je, Time, Time_next, coszen, &
!                        pwt, T1, P1, WS1, CO2, O3, rtnd_xactive)
!   </TEMPLATE>
!   <IN NAME="lon" TYPE="real" DIM="(:,:)">
!     Longitude of the center of the model gridcells
!   </IN>
!   <IN NAME="lat" TYPE="real" DIM="(:,:)">
!     Latitude of the center of the model gridcells
!   </IN>
!   <IN NAME="land" TYPE="real" DIM="(:,:)">
!     Fraction of land in grid cell
!   </IN>
!   <IN NAME="is,ie,js,je" TYPE="integer" DIM="(:,:)">
!     Local domain start and end indices
!   </IN>
!   <IN NAME="Time, Time_next" TYPE="type(time_type)">
!     Model time
!   </IN>
!   <IN NAME="coszen" TYPE="real" DIM="(:,:)">
!     Cosine of the solar zenith angle
!   </IN>
!   <IN NAME="pwt" TYPE="real" DIM="(:,:)">
!     Air mass in the bottom model layer
!   </IN>
!   <IN NAME="T1" TYPE="real" DIM="(:,:)">
!     Temperature in the bottom model layer
!   </IN>
!   <IN NAME="P1" TYPE="real" DIM="(:,:)">
!     Flux of photosynthetically active radiation
!   </IN>
!   <IN NAME="WS1" TYPE="real" DIM="(:,:)">
!     10-m wind speed
!   </IN>
!   <IN NAME="CO2" TYPE="real" DIM="(:,:)">
!     Surface CO2 concentration
!   </IN>
!   <IN NAME="O3" TYPE="real" DIM="(:,:)">
!     Surface ozone concentration
!   </IN>
!   <OUT NAME="rtnd_xactive" TYPE="real" DIM="(:,:,:)">
!     xactive tracer tendencies
!   </OUT>
!
subroutine xactive_bvoc( lon, lat, land, is, ie, js, je, Time, Time_next, coszen, &
                         pwtsfc, T1, P1, WS1, CO2, O3, rtnd_xactive, xbvoc4soa   )

   real, intent(in), dimension(:,:)            :: lon, lat        ! Longitude, latitude []
   real, intent(in), dimension(:,:)            :: land            ! Land fraction []
   integer, intent(in)                         :: is, ie, js, je  ! Local domain indices
   type(time_type), intent(in)                 :: Time, Time_next ! Time now, time next []
   real, intent(in), dimension(:,:)            :: coszen          ! cosine solar zenith angle []
   real, intent(in), dimension(:,:)            :: pwtsfc          ! Air mass in lowest layer [kg/m2]
   real, intent(in), dimension(:,:)            :: T1              ! surface temp [K]
   real, intent(in), dimension(:,:)            :: P1              ! PAR/PPFD [umoles/m2/s]
   real, intent(in), dimension(:,:)            :: WS1             ! 10m wind speed [m/s]
   real, intent(in), dimension(:,:)            :: CO2             ! surface CO2 conc. [VMR]
   real, intent(in), dimension(:,:)            :: O3              ! surface O3 conc.  [VMR]
   real, intent(out), dimension(:,:,:)         :: rtnd_xactive    ! xactive tracer tendencies [VMR/s]
   real, intent(out), dimension(:,:,:)         :: xbvoc4soa       ! biogenic emissions (for SOA) [molec/cm2/s]

!-------------------------------------------------------------------------------------------------
!-------------------------------------  Local Variables  -----------------------------------------

   real, dimension(size(T1,1),size(T1,2))      :: T24, P24, TMAX, TMIN
   real, dimension(size(T1,1),size(T1,2))      :: WSMAX, AQI
   real, dimension(size(T1,1),size(T1,2))      :: EMIS, EMIS_TERP
   real, dimension(size(T1,1),size(T1,2),nPFT) :: LAIp, LAIc   ! Used in MEGAN2
   real, dimension(size(T1,1),size(T1,2))      :: LAIp3, LAIc3      ! Used in MEGAN3
   integer, dimension(size(T1,1),size(T1,2))   :: DAY_BEGIN
   integer                                     :: yr, month, day, hr, minute, sec, month_p
   integer                                     :: nlon, nlat, i, j
   integer                                     :: xactive_knt, nTERP
   logical                                     :: used

   IF ( .NOT. module_is_initialized )    &
      call error_mesg ('xactive_bvoc',   &
           'xactive_bvoc_init must be called first.', FATAL)

! Get model time (GMT)
   call get_date(Time,yr,month,day,hr,minute,sec)
! Covert hour to an index (i.e., so it is non-zero)
   hr = hr + 1
! Set the previous month
   IF ( month == 1 ) THEN
      month_p = 12
   ELSE
      month_p = month - 1
   ENDIF

   nlon = size(lon,1)
   nlat = size(lat,2)
   if (Ldebug .and. mpp_pe()==mpp_root_pe()) then
      write(*,*) 'xactive_bvoc: nlat,nlon=', nlat,nlon
      write(*,*) 'xactive_bvoc: is,ie,js,je=', is,ie,js,je
      write(*,*) 'xactive_bvoc: lon1,lat1=', lon(1,1)*180./PI,lat(1,1)*180./PI
      write(*,*) 'xactive_bvoc: SIZE(T1)=', SIZE(T1)
      write(*,*) 'xactive_bvoc: SIZE(T24_STORE)=', SIZE(T24_STORE)
   endif

! Index for Local 8AM, the starting index for summing over O3 values
! to calculate the W126 Air Quality Index (8AM - 8PM)
! ..................................................................
   DO j = 1,nlat
   DO i = 1, nlon
      DAY_BEGIN(i,j) = int(ceiling((112.5 - lon(i,j)/15.)))
      IF ( DAY_BEGIN(i,j) .lt. 0 ) THEN
         DAY_BEGIN(i,j) = DAY_BEGIN(i,j) + 24
      ENDIF
      DAY_BEGIN(i,j) = DAY_BEGIN(i,j) + 1
   ENDDO
   ENDDO

! Update the Daily Average/Max Arrays, i.e, "___STORE"
! ......................................................
   IF ( do_ONLINE_TEMP ) THEN
      T24_STORE(is:ie,js:je,hr) = T1
   ENDIF
   IF ( do_ONLINE_PPFD ) THEN
      P24_STORE(is:ie,js:je,hr) = P1
   ENDIF
   IF ( do_ONLINE_WIND .AND. do_GAMMA_HW ) THEN
      WS_STORE(is:ie,js:je,hr)  = WS1
   ENDIF
   IF ( do_ONLINE_CO2 .AND. do_GAMMA_CO2) THEN
      CO2_STORE(is:ie,js:je) = CO2
   ENDIF
   IF ( do_ONLINE_O3 .AND. do_GAMMA_AQ ) THEN
! Convert to ppm and calculate W126
      O3_STORE(is:ie,js:je,hr) = (O3 * 1.e6) *  &
                         ( 1. / (1. + (4403. *exp(-126. * (O3*1.e6)))))
   ENDIF

! Calculate daily fields if necessary
! ...................................
   T24(:,:) = 0.
   P24(:,:) = 0.
   TMAX(:,:) = 0.
   TMIN(:,:) = 0.
   WSMAX(:,:) = 0.
   AQI(:,:) = 0.
   DO j = 1,nlat
   DO i = 1,nlon
! Only do for land fraction > 0.01
   IF ( land(i,j) > min_land_frac ) THEN
! Temperatures
      IF ( do_ONLINE_TEMP ) THEN
         T24(i,j)   = SUM(T24_STORE(i+is-1,j+js-1,:))/24.
         IF ( do_GAMMA_HT ) THEN
            TMAX(i,j)  = MAXVAL(T24_STORE(i+is-1,j+js-1,:))
         ENDIF
         IF ( do_GAMMA_LT ) THEN
             TMIN(i,j)  = MINVAL(T24_STORE(i+is-1,j+js-1,:))
         ENDIF
      ELSE
         T24(i,j) = T1(i,j)
      ENDIF
! PAR/PPFD
      IF ( do_ONLINE_PPFD ) THEN
         P24(i,j)   = SUM(P24_STORE(i+is-1,j+js-1,:))/24.
      ELSE
         P24(i,j)   = P1(i,j)
      ENDIF
! Wind Speed
      IF ( do_ONLINE_WIND .AND. do_GAMMA_HW ) THEN
         WSMAX(i,j) = MAXVAL(WS_STORE(i+is-1,j+js-1,:))
      ENDIF
! Air quality
      IF ( do_ONLINE_O3 .AND. do_GAMMA_AQ ) THEN
         IF ( DAY_BEGIN(i,j) > 13 ) THEN
           AQI(i,j) = SUM(O3_STORE(i+is-1,j+js-1,DAY_BEGIN(i,j):24)) +             &
                      SUM(O3_STORE(i+is-1,j+js-1,1:(12 - (25 - DAY_BEGIN(i,j)))))
         ELSE
           AQI(i,j) = SUM(O3_STORE(i+is-1,j+js-1,DAY_BEGIN(i,j):DAY_BEGIN(i,j)+11))
         ENDIF
      ENDIF
   ENDIF ! land fraction
   ENDDO ! lon/i
   ENDDO ! lat/j

! Update the LAI data, using either online or from read in data
! .............................................................
   IF ( do_ONLINE_LAI ) THEN
   ! Do something that brings in LAI from the land model
   ! Then either format for MEGAN2 LAI = (lon x lat x pft(17))
   ! Or for MEGAN3, LAIv = LAI/FCOVER = (lon x lat)
   ELSE
      IF ( xactive_algorithm == 'MEGAN2' .or. do_AM3_ISOP ) THEN
         LAIc = MLAI(is:ie,js:je,:,month)
         LAIp = MLAI(is:ie,js:je,:,month_p)
      ENDIF
      IF (xactive_algorithm == 'MEGAN3' ) THEN
         LAIc3 = MLAI_MEGAN3(is:ie,js:je,month)
         LAIp3 = MLAI_MEGAN3(is:ie,js:je,month_p)
      ENDIF
   ENDIF

   IF ( do_ONLINE_PFT ) THEN
      ! call send_PFT ( PCTPFT )
   ENDIF

! Intialized the xactive counter
! ..............................
   xactive_knt = 0
! Initialize terpene emissions (in the case that they are parsed)
! ................................................................
   IF ( ALLOCATED(ECTERP) ) THEN
      nTERP = size(ECTERP,4)
   ENDIF

!--------------------------------------------------------------------
!  MAIN LOOP
!  Determine if it has xactive emissions, if so increase the counter
!  and store the index to be sent back to the tracer driver
!  Call the appropriate routine (megan version, terpenes) and return
!  emissions, convert to tendency to be sent back to tracer driver
!....................................................................
   rtnd_xactive(:,:,:) = 0.
   xbvoc4soa(:,:,:) = 0.
   DO i = 1, pcnstm1

      IF ( has_xactive_emis(i) ) THEN

         EMIS(:,:) = 0.
         xactive_knt = xactive_knt + 1

         IF ( trim(tracnam(i))=='DMS' ) THEN
            ! SKIP - calculated in tropchem driver
         ELSEIF ( do_AM3_ISOP .AND. trim(tracnam(i))=='ISOP' ) THEN
! Reproduces AM3 isoprene emissions
            call calc_xactive_bvoc_AM3 ( Time, Time_next, is, js,         &
                                         lon, lat, land, coszen,          &
                                         P1, T1, LAIp, LAIc,              &
                                         Pmo(is:ie,js:je,:),              &
                                         Tmo(is:ie,js:je,:),              &
                                         ECISOP_AM3(is:ie,js:je,:), month, EMIS, &
                                         id_GAMMA_TEMP=id_G_TEMP(i),      &
                                         id_GAMMA_PAR=id_G_PAR(i),        &
                                         id_GAMMA_LAI=id_G_LAI(i),        &
                                         id_GAMMA_AGE=id_G_AGE(i))
         ELSEIF ( do_PARSED_TERP .AND. trim(tracnam(i))=='C10H16' ) THEN
! Calculates the emissions of each terpene separately
            DO j = 1, nTERP
               EMIS_TERP(:,:) = 0.
               IF ( xactive_algorithm == 'MEGAN2' ) THEN
                  call calc_xactive_bvoc_megan2 ( Time, Time_next, is, js,      &
                                                  lon, lat, land, coszen,       &
                                                  P1, P24, T1, T24, LAIp, LAIc, &
                                                  Tmo(is:ie,js:je,:),           &
                                                  TERP_PARAM(:,j),              &
                                                  ECTERP(is:ie,js:je,:,j),      &
                                                  month, tracnam(i), EMIS_TERP, &
                                                  id_GAMMA_TEMP=id_G_TEMP(i),   &
                                                  id_GAMMA_PAR=id_G_PAR(i),     &
                                                  id_GAMMA_LAI=id_G_LAI(i),     &
                                                  id_GAMMA_AGE=id_G_AGE(i),     &
                                                  id_GAMMA_CO2=id_G_CO2(i),     &
                                                  id_GAMMA_SM=id_G_SM(i))
               ELSEIF ( xactive_algorithm == 'MEGAN3' ) THEN
                  call calc_xactive_bvoc_megan3 ( Time, Time_next, is, js,      &
                                                  lon, lat, land, coszen,       &
                                                  P1, P24, T1, T24,             &
                                                  LAIp3, LAIc3,                 &
                                                  TMAX, TMIN,                   &
                                                  Tmo(is:ie,js:je,:), WSMAX, AQI, &
                                                  TERP_PARAM(:,j),              &
                                                  LDFg_TERP(is:ie,js:je,j),     &
                                                  ECTERP_MEGAN3(is:ie,js:je,j), &
                                                  month, tracnam(i), EMIS_TERP, &
                                                  id_GAMMA_TEMP=id_G_TEMP(i),   &
                                                  id_GAMMA_PAR=id_G_PAR(i),     &
                                                  id_GAMMA_LAI=id_G_LAI(i),     &
                                                  id_GAMMA_AGE=id_G_AGE(i),     &
                                                  id_GAMMA_BDLAI=id_G_BDLAI(i), &
                                                  id_GAMMA_CO2=id_G_CO2(i),     &
                                                  id_GAMMA_AQ=id_G_AQ(i),       &
                                                  id_GAMMA_SM=id_G_SM(i),       &
                                                  id_GAMMA_HT=id_G_HT(i),       &
                                                  id_GAMMA_LT=id_G_LT(i),       &
                                                  id_GAMMA_HW=id_G_HW(i))
                ENDIF ! /megan version
                EMIS(:,:) = EMIS(:,:) + EMIS_TERP(:,:)
             ENDDO !/nTERP
         ELSE
            IF ( xactive_algorithm == 'MEGAN2' ) THEN
            call calc_xactive_bvoc_megan2 ( Time, Time_next, is, js,         &
                                            lon, lat, land, coszen,          &
                                            P1, P24, T1, T24, LAIp, LAIc,    &
                                            Tmo(is:ie,js:je,:),              &
                                            MEGAN_PARAM(:,xactive_knt),      &
                                            ECBVOC(is:ie,js:je,:,xactive_knt), &
                                            month, tracnam(i), EMIS,         &
                                            id_GAMMA_TEMP=id_G_TEMP(i),      &
                                            id_GAMMA_PAR=id_G_PAR(i),        &
                                            id_GAMMA_LAI=id_G_LAI(i),        &
                                            id_GAMMA_AGE=id_G_AGE(i),        &
                                            id_GAMMA_CO2=id_G_CO2(i),        &
                                            id_GAMMA_SM=id_G_SM(i))
            ELSEIF ( xactive_algorithm == 'MEGAN3' ) THEN
            call calc_xactive_bvoc_megan3 ( Time, Time_next, is, js,         &
                                            lon, lat, land, coszen,          &
                                            P1, P24, T1, T24,                &
                                            LAIp3, LAIc3,                    &
                                            TMAX, TMIN,                      &
                                            Tmo(is:ie,js:je,:), WSMAX, AQI,  &
                                            MEGAN_PARAM(:,xactive_knt),      &
                                            LDFg(is:ie,js:je,xactive_knt),   &
                                            ECBVOC_MEGAN3(ie:ie,js:je,xactive_knt), &
                                            month, tracnam(i), EMIS,         &
                                            id_GAMMA_TEMP=id_G_TEMP(i),      &
                                            id_GAMMA_PAR=id_G_PAR(i),        &
                                            id_GAMMA_LAI=id_G_LAI(i),        &
                                            id_GAMMA_AGE=id_G_AGE(i),        &
                                            id_GAMMA_BDLAI=id_G_BDLAI(i),    &
                                            id_GAMMA_CO2=id_G_CO2(i),        &
                                            id_GAMMA_AQ=id_G_AQ(i),          &
                                            id_GAMMA_SM=id_G_SM(i),          &
                                            id_GAMMA_HT=id_G_HT(i),          &
                                            id_GAMMA_LT=id_G_LT(i),          &
                                            id_GAMMA_HW=id_G_HW(i))
            ENDIF !/megan version
         ENDIF !/species
         ! Send emissions diagnostics
         if ( trim(tracnam(i))=='ISOP' ) then
            if (abs(scale_isoprene_emissions - 1.).gt.epsln) then
               EMIS = EMIS*scale_isoprene_emissions
            end if
         end if
         if ( trim(tracnam(i))=='C10H16' ) then
            if (abs(scale_terpene_emissions - 1.).gt.epsln) then
               EMIS = EMIS*scale_terpene_emissions
            end if
         end if

         
         IF ( id_EMIS(i) > 0 ) THEN
            used = send_data ( id_EMIS(i), EMIS, Time_next, is_in=is, js_in=js)
         ENDIF
! Convert Emisisons (molecules/cm2/s) to VMR/s
!  rdt(VMR/s) = EMIS(molec/cm2/s) * 1e4(cm2/m2) / &
!                 ( pwt(kg/m2) / WTMAIR(g/mol) * 1e3(g/kg) * AVOGNO(molec/mole) )
           rtnd_xactive(:,:,xactive_knt) = EMIS * 10. / pwtsfc * WTMAIR / AVOGNO
           IF ( trim(tracnam(i))=='ISOP' ) THEN
              xbvoc4soa(:,:,ind_xbvoc_ISOP) = EMIS
           ELSEIF ( trim(tracnam(i))=='C10H16' ) THEN
              xbvoc4soa(:,:,ind_xbvoc_TERP) = EMIS
           ENDIF

      ENDIF !has_xactive_emis
   ENDDO !pctnstm1

end subroutine xactive_bvoc
!</SUBROUTINE>


!########################################################################################

! <SUBROUTINE NAME="xactive_bvoc_init">
!   <OVERVIEW>
!     Initializes the interactive biogenic VOC emission module
!   </OVERVIEW>
!   <DESCRIPTION>
!     This subroutine initializes the BVOC emission module.
!     It is called from atmos_tracer_driver_init.
!     1) Namelist variables are parsed
!     2) MEGAN Model paramters and emission factors are allocated, read in, & stored
!     3) Diagnostics IDs are assigned
!     4) Datasets of monthly average temprerature, PAR, LAI, and fractional grid cell
!        coverage of plant functional types are read in (if required) by seperate
!        initialization subroutines.
!   </DESCRIPTION>
!   <TEMPLATE>
!      call xactive_bvoc_init ( lonb, latb, Time, axes, nxactive )
!   </TEMPLATE>
!   <IN NAME="lonb" TYPE="real" DIM="(:,:)">
!     Longitude corners for the local domain
!   </IN>
!   <IN NAME="latb" TYPE="real" DIM="(:,:)">
!     Latitude corners for the local domain
!   </IN>
!   <IN NAME="Time" TYPE="type(time_type)">
!     Model time
!   </IN>
!   <IN NAME="axes" TYPE="integer" DIM="(4)">
!     The axes relating to the tracer array
!   </IN>
!   <OUT NAME="xactive_ndx" TYPE="integer" DIM="(:)">
!     Index/Location of each xactive species in
!     the tracer array
!   </OUT>

subroutine xactive_bvoc_init(domain, lonb, latb, Time, axes, xactive_ndx)

   type(domain2D),target,intent(in)    :: domain !< Atmosphere domain
   real, intent(in), dimension(:,:)    :: lonb, latb     ! Lat/Lon corners
   type(time_type), intent(in)         :: Time           ! Model time
   integer, intent(in)                 :: axes(4)        ! Diagnostics axes
   integer, intent(out), dimension(:)  :: xactive_ndx    ! index into tracer array

!----------------Local Variables---------------------------------------------------------

   character(len=5) :: pftnames(nPFT)         =  (/ 'pft01','pft02','pft03','pft04',             &
                                                    'pft05','pft06','pft07','pft08',             &
                                                    'pft09','pft10','pft11','pft12',             &
                                                     'pft13','pft14','pft15','pft16', 'pft17'/)
!  character(len=3) :: vegnames(nVEG)         =  (/ 'ntr', 'btr', 'crp', 'grs', 'shr' /)
   character(len=3) :: vegnames(nVEG)         =  (/ 'ntr', 'btr', 'shr', 'grs', 'crp' /)

   character(len=7) :: terpnames_megan3(8)    =  (/'MT_PINE', 'MT_ACYC', 'MT_CAMP',         &
                                                   'MT_SABI', 'MT_AROM', 'MT_OXY ',         &
                                                   'SQT_HR ', 'SQT_LR '/)

   character(len=4) :: terpnames_megan2(11)   =  (/'MYRC','SABI','LIMO','CARE','OCIM', &
                                                    'BPIN','APIN','OMTP','FARN','CARY','OSQT'/)


   character(len=5) :: paramnames_megan3(20)  =  (/ 'BETA ', 'C_t1 ', 'C_eo ', 'A_new',        &
                                                    'A_gro', 'A_mat', 'A_old', 'C_AQ ',     &
                                                    'C_HW ', 'C_HT ', 'C_LT ', 'T_AQ ',        &
                                                    'T_HW ', 'T_HT ', 'T_LT ', 'DT_AQ',       &
                                                    'DT_HW', 'DT_HT', 'DT_LT', 'mw   '  /)
   character(len=5) :: paramnames_megan2(9)   =   (/'BETA ', 'LDF  ', 'C_t1 ', 'C_eo ',        &
                                                    'A_new', 'A_gro', 'A_mat', 'A_old',    &
                                                    'mw   '/)

   integer          :: nlon, nlat, i, j, k, n, xknt, nTERP, nxactive
   integer          :: ierr, io, logunit, nPARAMS

   integer, parameter             :: nlonin = 720, nlatin = 360
   real, dimension(nlonin)        :: inlon
   real, dimension(nlatin)        :: inlat
   real, dimension(nlonin+1)      :: inlone
   real, dimension(nlatin+1)      :: inlate
   real, dimension(nlonin,nlatin) :: AM3_ISOP_DATAIN

! Higher resolution input data can be created upon request (jschnell)
   integer, parameter                 :: m3nlonin = 720, m3nlatin = 360
   real, dimension(m3nlonin)          :: m3inlon
   real, dimension(m3nlatin)          :: m3inlat
   real, dimension(m3nlonin+1)        :: m3inlone
   real, dimension(m3nlatin+1)        :: m3inlate
   real, dimension(m3nlonin,m3nlatin) :: MEGAN3_DATAIN

   real               :: dlon, dlat
   real               :: megan2_isop_sf
   real               :: toss
   character(len=64)  :: ecfile
   character(len=256) :: control=''
   character(len=64)  :: name=''
   integer, parameter :: nPARAMS_megan2 = 9
   integer, parameter :: nPARAMS_megan3 = 20
   type(FmsNetcdfFile_t)       ::  Xbvoc_restart !< Netcdf fileobj
   type(FmsNetcdfFile_t)       ::  megan3_isop !< Netcdf fileobj
   type(FmsNetcdfFile_t)       ::  ecfile_obj !< Netcdf fileobj
   type(FmsNetcdfDomainFile_t) ::  Til_restart !< Domain decomposed fileobj
   integer, allocatable, dimension(:) :: pes !< Array of pes in the current pelist
   type (horiz_interp_type)    :: Interp

   nlon = size(lonb,1) - 1
   nlat = size(latb,2) - 1
   if (Ldebug .and. mpp_pe()==mpp_root_pe()) then
      write(*,*) 'xactive_bvoc_init: nlat,nlon=', nlat,nlon
      write(*,*) 'xactive_bvoc: lonb1,latb1=', lonb(1,1)*180./PI,latb(1,1)*180./PI
   endif

   IF ( module_is_initialized ) RETURN

!-----------------------------------------------------------------------
!     ... write version number
!-----------------------------------------------------------------------
   call write_version_number(version, tagname)

!-----------------------------------------------------------------------
!     ... read namelist
!-----------------------------------------------------------------------
   IF (file_exists('input.nml')) THEN
      read (input_nml_file, nml=xactive_bvoc_nml, iostat=io)
      ierr = check_nml_error(io,'xactive_bvoc_nml')
   ENDIF

   logunit = stdlog()
   IF(mpp_pe() == mpp_root_pe()) THEN
      write(logunit, nml=xactive_bvoc_nml)
      verbose = verbose + 1
   ENDIF

! Overwrite ( do_ONLINE_TEMP/PPFD ) since currently this is needed
! to get daily averages as there are no files available
   IF ( xactive_algorithm == 'MEGAN3' ) THEN
      do_ONLINE_TEMP = .TRUE.
      do_ONLINE_PPFD = .TRUE.
   ENDIF
! Set a few version specific parameters
   IF ( xactive_algorithm == 'MEGAN2' ) THEN
      nPARAMS = nPARAMS_megan2
      IF ( fix_megan2_isop ) THEN
         megan2_isop_sf = 1.5
      ELSE
         megan2_isop_sf = 1.0
      ENDIF
   ELSE IF ( xactive_algorithm == 'MEGAN3' ) THEN
      nPARAMS = nPARAMS_megan3
   ENDIF

!-----------------------------------------------------------------------
!     ... Set up xactive emissions, diagnostics, etc.
!-----------------------------------------------------------------------
   nxactive = SIZE(xactive_ndx)
   ALLOCATE ( MEGAN_PARAM(nPARAMS,nxactive) )
   IF ( xactive_algorithm == 'MEGAN2' ) THEN
      ALLOCATE ( ECBVOC(nlon,nlat,nPFT,nxactive))
   ELSEIF ( xactive_algorithm == 'MEGAN3' ) THEN
      ALLOCATE ( ECBVOC_MEGAN3(nlon,nlat,nxactive) )
      ALLOCATE ( LDFg(nlon,nlat,nxactive) )
! Populate these so it doesn't need to be done each time a new species is read in
      if (open_file(megan3_isop,"INPUT/megan3.xactive.ISOP.nc","read")) then
         call read_data (megan3_isop, 'lon', m3inlon)
         call read_data (megan3_isop, 'lat', m3inlat)
         call close_file(megan3_isop)
      endif
      m3inlon = m3inlon*DEG_TO_RAD
      m3inlat = m3inlat*DEG_TO_RAD
      dlat = m3inlat(2)-m3inlat(1)
      dlon = m3inlon(2)-m3inlon(1)
      m3inlone(1:m3nlonin) = m3inlon-(dlon/2.)
      m3inlone(m3nlonin+1) = m3inlon(m3nlonin)+(dlon/2.)
      m3inlate(1:m3nlatin) = m3inlat-(dlat/2.)
      m3inlate(m3nlatin+1) = m3inlat(m3nlatin)+(dlat/2.)
   END IF
   IF ( do_AM3_ISOP ) THEN
      ALLOCATE( ECISOP_AM3(nlon,nlat,nVEG) )
   ENDIF
!----------------------------------------------------------------------
!     ... Set up the required arrays if parsing terpenes
!---------------------------------------------------------------------
   IF ( do_PARSED_TERP ) THEN
      IF ( do_SESQTERP ) THEN
         IF ( xactive_algorithm == 'MEGAN2' ) THEN
             nTERP = 11
         ELSE
             nTERP = 8
         ENDIF
      ELSE
         IF ( xactive_algorithm == 'MEGAN2' ) THEN
            nTERP = 8
         ELSE
            nTERP = 6
         ENDIF
      ENDIF
      ALLOCATE( TERP_PARAM(nPARAMS,nTERP) )
      IF ( xactive_algorithm == 'MEGAN2' ) THEN
         ALLOCATE( ECTERP(nlon,nlat,nPFT,nTERP) )
      ELSE IF ( xactive_algorithm == 'MEGAN3' ) THEN
         ALLOCATE( ECTERP_MEGAN3(nlon,nlat,nTERP) )
         ALLOCATE( LDFg_TERP (nlon,nlat,nTERP) )
      ENDIF
   ENDIF


   if (use_isop_shrub_crop_bug) vegnames(:) =  (/ 'ntr', 'btr', 'crp', 'grs', 'shr' /)

   indices(:) = 0
   xknt = 0
   DO i = 1, pcnstm1
      n = get_tracer_index(MODEL_ATMOS, tracnam(i))
      if (Ldebug .and. mpp_pe()==mpp_root_pe()) &
         write(*,*) 'xactive_bvoc_init:', TRIM(tracnam(i)),i,n
      IF ( n .le. 0 ) THEN
         IF ( mpp_pe()==mpp_root_pe()) call error_mesg('xactive_bvoc_init',       &
              trim(tracnam(i)) // ' is not found', WARNING)
         cycle
      ENDIF
      indices(i) = n
      has_xactive_emis(i) = query_method('xactive_emissions',MODEL_ATMOS,         &
                                         indices(i),name,control)
      IF ( has_xactive_emis(i) ) THEN
         xknt = xknt + 1
         xactive_ndx(xknt) = get_tracer_index(MODEL_ATMOS,trim(tracnam(i)))
      ENDIF

      IF ( trim(tracnam(i))=='DMS' ) THEN
         IF ( mpp_pe()==mpp_root_pe()) call error_mesg('xactive_bvoc_init',       &
              'skipping set up for non-BVOC tracer '//trim(tracnam(i)),NOTE)
         cycle
      ENDIF

! Register the diagnostics for emissions and all possible gammas
      IF ( has_xactive_emis(i) ) THEN
         IF ( mpp_pe()==mpp_root_pe()) call error_mesg('xactive_bvoc_init',       &
              'Initializing xactive emissions for '//trim(tracnam(i)),NOTE)

! Emissions and standard gamma diagnostics for all species
         id_EMIS(i)         = register_diag_field(module_name,                    &
                              trim(tracnam(i))//'_xactive_emis', axes(1:2),       &
                              Time, trim(tracnam(i))//'_xactive_emis',            &
                              'molecules/cm2/s')
         id_G_TEMP(i)       = register_diag_field(module_name,                    &
                              trim(tracnam(i))//'_gamma_temp',                    &
                              axes(1:2), Time, trim(tracnam(i))//'_gamma_temp',   &
                              'unitless')
         id_G_PAR(i)        = register_diag_field(module_name,                    &
                              trim(tracnam(i))//'_gamma_light',                   &
                              axes(1:2), Time, trim(tracnam(i))//'_gamma_light',  &
                              'unitless')
         id_G_LAI(i)        = register_diag_field(module_name,                    &
                              trim(tracnam(i))//'_gamma_lai',                     &
                              axes(1:2), Time, trim(tracnam(i))//'_gamma_lai',    &
                              'unitless')
         id_G_AGE(i)        = register_diag_field(module_name,                    &
                              trim(tracnam(i))//'_gamma_age',                     &
                              axes(1:2), Time, trim(tracnam(i))//'_gamma_age',    &
                              'unitless')
! Optional gamma diagnostics (only for MEGAN3)
         IF ( do_GAMMA_HT ) THEN
            id_G_HT(i)      = register_diag_field(module_name,                    &
                              trim(tracnam(i))//'_gamma_high_temp', axes(1:2),    &
                              Time, trim(tracnam(i))//'_gamma_high_temp',         &
                              'unitless')
         ENDIF
         IF ( do_GAMMA_LT ) THEN
            id_G_LT(i)      = register_diag_field(module_name,                    &
                              trim(tracnam(i))//'_gamma_low_temp', axes(1:2),     &
                              Time, trim(tracnam(i))//'_gamma_low_temp',          &
                              'unitless')
         ENDIF
         IF ( do_GAMMA_HW ) THEN
            id_G_HW(i)      = register_diag_field(module_name,                    &
                              trim(tracnam(i))//'_gamma_high_wind', axes(1:2),    &
                              Time, trim(tracnam(i))//'_gamma_high_wind',         &
                              'unitless')
         ENDIF
         IF ( do_GAMMA_SM ) THEN
            id_G_SM(i)      = register_diag_field(module_name,                    &
                              trim(tracnam(i))//'_gamma_soil', axes(1:2),         &
                              Time, trim(tracnam(i))//'_gamma_soil',              &
                              'unitless')
         ENDIF
         IF ( do_GAMMA_AQ ) THEN
            id_G_AQ(i)      = register_diag_field(module_name,                    &
                              trim(tracnam(i))//'_gamma_AQI', axes(1:2),          &
                              Time, trim(tracnam(i))//'_gamma_AQI',               &
                              'unitless')
         ENDIF
         IF ( do_GAMMA_CO2 ) THEN
            id_G_CO2(i)     = register_diag_field(module_name,                    &
                              trim(tracnam(i))//'_gamma_high_temp', axes(1:2),    &
                              Time, trim(tracnam(i))//'_gamma_high_temp',         &
                              'unitless')
         ENDIF
         IF ( do_GAMMA_BDLAI) THEN
            id_G_BDLAI(i)     = register_diag_field(module_name,                  &
                                trim(tracnam(i))//'_gamma_BDLAI', axes(1:2),      &
                                Time, trim(tracnam(i))//'_gamma_BDLAI',           &
                                'unitless')
         ENDIF

!--------------------------------------------------------------------------------------
!  ... Read in the MEGAN model paramters and emission capacities
!  ... NOTE: Reproducing AM3 isoprene is a special case
!  ... NOTE: terpenes only vs. monoterpenes + sequisterpenes
!  ... >>>>>>> parsed vs. lumped terpenes
!--------------------------------------------------------------------------------------
         IF ( trim(tracnam(i))=='ISOP' .AND. do_AM3_ISOP ) THEN
            ecfile = 'INPUT/megan.ISOP.nc'
            if (open_file(ecfile_obj,ecfile,"read")) then
               call read_data (ecfile_obj, 'lon', inlon)
               call read_data (ecfile_obj, 'lat', inlat)
               inlon = inlon*DEG_TO_RAD
               inlat = inlat*DEG_TO_RAD
               dlat = inlat(2)-inlat(1)
               dlon = inlon(2)-inlon(1)
               inlone(1:nlonin) = inlon-(dlon/2.)
               inlone(nlonin+1) = inlon(nlonin)+(dlon/2.)
               inlate(1:nlatin) = inlat-(dlat/2.)
               inlate(nlatin+1) = inlat(nlatin)+(dlat/2.)
               call horiz_interp_init
               call horiz_interp_new ( Interp, inlone, inlate, lonb, latb )
               DO j = 1, nVEG
                  call read_data (ecfile_obj,vegnames(j),AM3_ISOP_DATAIN)
                  call horiz_interp (Interp,AM3_ISOP_DATAIN,ECISOP_AM3(:,:,j), verbose=verbose)
               ENDDO
               call horiz_interp_del( Interp )
               call close_file(ecfile_obj)
            ELSE
                  call error_mesg ('xactive_bvoc_init',  &
                     ' AM3 isoprene emission capacity file does not exist', FATAL)
            ENDIF
         ELSE IF ( trim(tracnam(i))=='C10H16') THEN
              IF ( xactive_algorithm == 'MEGAN2' ) THEN
                 IF ( do_PARSED_TERP ) THEN
! Both mono- and sesq- terpenes are included in this file,
! but sesq may not be used (i.e., if do_SESQTERP = .false.')
                       ecfile = 'INPUT/megan2.xactive.parsed_terpenes.nc'
                       if (open_file(ecfile_obj,ecfile,"read")) then
                        DO k = 1, nTERP
                           DO j = 1, nPFT
                              call read_data(ecfile_obj,terpnames_megan2(k)//'_'//pftnames(j), &
                                             toss)
                              ECTERP(:,:,j,k) = toss
                           ENDDO
                        ENDDO
                        call close_file(ecfile_obj)
                       else
                        call error_mesg ('xactive_bvoc_init',  &
                        'MEGAN file (megan2.xactive.parsed_terpenes.nc) for '//trim(tracnam(i))//' does not exist', FATAL)
                       endif
                 ELSE
                    IF ( do_SESQTERP ) THEN
                       ecfile = 'INPUT/megan2.xactive.lumped_terpenes_sesq.nc'
                    ELSE
                       ecfile = 'INPUT/megan2.xactive.lumped_terpenes_mono.nc'
                    ENDIF
                    if (open_file(ecfile_obj,ecfile,"read")) then
                     DO j = 1, nPFT
                        call read_data(ecfile_obj,pftnames(j),toss)
                        ECBVOC(:,:,j,xknt) = toss
                     ENDDO
                     call close_file(ecfile_obj)
                    else
                     call error_mesg ('xactive_bvoc_init',  &
                        'MEGAN file '//trim(ecfile)//'for '//trim(tracnam(i))//' does not exist', FATAL)
                    endif
                 ENDIF
              ELSE IF ( xactive_algorithm == 'MEGAN3' ) THEN
                 IF ( do_PARSED_TERP ) THEN
! Both mono- and sesq- terpenes are included in this file,
! but sesq may not be used (i.e., if do_SESQTERP = .false.')
                    ecfile = 'INPUT/megan3.xactive.parsed_terpenes.nc'
                    if (open_file(ecfile_obj,ecfile,"read")) then
                     call horiz_interp_init
                     call horiz_interp_new ( Interp, m3inlone, m3inlate, lonb, latb )
                     DO k = 1, nTERP
                        call read_data (ecfile_obj,trim(terpnames_megan3(k))//'_EF',        &
                                        MEGAN3_DATAIN)
                        call horiz_interp (Interp,MEGAN3_DATAIN,                  &
                                           ECTERP_MEGAN3(:,:,k), verbose=verbose)
                        call read_data (ecfile_obj,trim(terpnames_megan3(k))//'_LDF',       &
                                        MEGAN3_DATAIN)
                        call horiz_interp (Interp,MEGAN3_DATAIN,                  &
                                           LDFg_TERP(:,:,k), verbose=verbose)
                     ENDDO!nterp
                     call horiz_interp_del( Interp )
                     call close_file(ecfile_obj)
                    else
                     call error_mesg ('xactive_bvoc_init',  &
                        'MEGAN file '//trim(ecfile)//'for '//trim(tracnam(i))//' does not exist', FATAL)
                    endif
                 ELSE
                    IF ( do_SESQTERP ) THEN
                       ecfile = 'INPUT/megan3.xactive.lumped_terpenes_sesq.nc'
                    ELSE
                       ecfile = 'INPUT/megan3.xactive.lumped_terpenes_mono.nc'
                    ENDIF
                    if (open_file(ecfile_obj,ecfile,"read")) then
                     call horiz_interp_init
                     call horiz_interp_new ( Interp, m3inlone, m3inlate, lonb, latb )
                     call read_data (ecfile_obj,'EF',MEGAN3_DATAIN)
                     call horiz_interp (Interp, MEGAN3_DATAIN,ECBVOC_MEGAN3(:,:,xknt))
                     call read_data (ecfile_obj,'LDF',MEGAN3_DATAIN)
                     call horiz_interp (Interp, MEGAN3_DATAIN, LDFg(:,:,xknt ))
                     call horiz_interp_del( Interp )
                     call close_file(ecfile_obj)
                    else
                     call error_mesg ('xactive_bvoc_init',  &
                        'MEGAN file '//trim(ecfile)//'for '//trim(tracnam(i))//' does not exist', FATAL)
                    endif
                 ENDIF
              ENDIF
         ELSE
            IF ( xactive_algorithm == 'MEGAN2') THEN
               ecfile = 'INPUT/megan2.xactive.'//trim(tracnam(i))//'.nc'
               if (open_file(ecfile_obj,ecfile,"read")) then
                  IF (mpp_pe() == mpp_root_pe()) call error_mesg ( 'xactive_bvoc_init', &
                      'Reading EF from file ' //ecfile, NOTE)
                  DO j = 1, nPFT
                     call read_data(ecfile_obj,pftnames(j),toss)
                     IF ( trim(tracnam(i)) == 'ISOP' ) THEN
                        toss = megan2_isop_sf * toss
                     ENDIF
                     ECBVOC(:,:,j,xknt) = toss
                  ENDDO
                  call close_file(ecfile_obj)
               ELSE
                  call error_mesg ('xactive_bvoc_init',  &
                     'MEGAN file for '//trim(tracnam(i))//' does not exist', FATAL)
               ENDIF
            ELSE IF ( xactive_algorithm == 'MEGAN3') THEN
               ecfile = 'INPUT/megan3.xactive.'//trim(tracnam(i))//'.nc'
               IF (open_file(ecfile_obj,ecfile,"read")) then
                  IF (mpp_pe() == mpp_root_pe()) call error_mesg ( 'xactive_bvoc_init', &
                      'Reading EF from file ' //ecfile, NOTE)
                     call horiz_interp_init
                     call horiz_interp_new ( Interp, m3inlone, m3inlate, lonb, latb )
                     call read_data (ecfile_obj,'EF',MEGAN3_DATAIN)
                     call horiz_interp (Interp,MEGAN3_DATAIN,ECBVOC_MEGAN3(:,:,xknt), verbose=verbose)
                     call read_data (ecfile_obj,'LDF',MEGAN3_DATAIN)
                     call horiz_interp (Interp,MEGAN3_DATAIN,LDFg(:,:,xknt), verbose=verbose)
                     call horiz_interp_del( Interp )
                     call close_file(ecfile_obj)
               ENDIF
            ENDIF ! xactive_algorithm
         ENDIF ! AM3 isop, terpene, other

! Read in all of the megan model parameters
!----------------------------------------------------------------------
!----------------------------------------------------------------------
!                ... MEGANv2.1 (MEGAN2)
!----------------------------------------------------------------------
! BETA:  temperature coefficient (emission type 2: light independent)
! LDF:   light dependent fraction (scalar in MEGAN2, gridded in MEGAN3)
! C_t1:  temperature coefficient (emission type 1: light dependent)
! C_eo:  temperature coefficient (emission type 1: light dependent)
! A_new: coefficient for 'new' leaves
! A_gro: coefficient for 'growing' leaves
! A_mat: coefficient for 'mature' leaves
! A_old: coefficient for old/senescing leaves
!-----------------------------------------------------------------------
!               ... NEW to MEGANv3 (MEGAN3)
!-----------------------------------------------------------------------
! C_AQ:  coefficient for poor air quality stress
! C_HW:  coefficient for high wind speed stress
! C_HT:  coefficient for high temperature stress
! C_LT:  coefficient for low temperature stress
! T_AQ:  threshold for poor air quality stress        (units = ppm-hours)
! T_HW:  threshold for high wind speed stress         (units = m/s)
! T_HT:  threshold for high temperature stress        (units = Celsius)
! T_LT:  threshold for low temperature stress         (units = Celsius)
! DT_AQ: delta threshold for poor air quality stress  (units = ppm-hours)
! DT_HW: delta threshold for high wind speed stress   (units = m/s)
! DT_HT: delta threshold for high temperature stress  (units = Celsius)
! DT_LT: delta threshold for low temperature stress   (units = Celsius)
! mw: Molecular weigth (units = g/mole)
!------------------------------------------------------------------------
         IF ( trim(tracnam(i))=='ISOP' .and. do_AM3_ISOP ) THEN
            IF ( mpp_pe() == mpp_root_pe()) call error_mesg ('xactive_bvoc_init', &
                 'MEGAN Parameters for AM3 ISOP hardcoded in subroutine, skipping',NOTE)
         ELSE
            IF (.not. open_file(ecfile_obj,ecfile,"read")) then
               call error_mesg ('xactive_bvoc_init',  &
                        'File '//trim(ecfile)//'for '//trim(tracnam(i))//' does not exist', FATAL)
            ENDIF
            DO j = 1, nPARAMS
               IF ( trim(tracnam(i))=='C10H16') THEN
                  IF ( do_PARSED_TERP ) THEN
                     DO k = 1, nTERP
                        IF ( xactive_algorithm == 'MEGAN2' ) THEN
                           call read_data(ecfile_obj, &
                                trim(paramnames_megan2(j))//'_'//terpnames_megan2(k), &
                                          TERP_PARAM(j,k))
                        ELSE IF ( xactive_algorithm == 'MEGAN3' ) THEN
                           call read_data(ecfile_obj, &
                                trim(paramnames_megan3(j))//'_'//trim(terpnames_megan3(k)), &
                                       TERP_PARAM(j,k))
                        ENDIF
                     ENDDO
                  ELSE
                     IF ( xactive_algorithm == 'MEGAN2' ) THEN
                        call read_data(ecfile_obj, trim(paramnames_megan2(j)), &
                                       MEGAN_PARAM(j,xknt))
                     ELSE IF ( xactive_algorithm == 'MEGAN3' ) THEN
                        call read_data(ecfile_obj, trim(paramnames_megan3(j)), &
                                       MEGAN_PARAM(j,xknt))
                     ENDIF
                  ENDIF
               ELSE
                     IF ( xactive_algorithm == 'MEGAN2' ) THEN
                        call read_data(ecfile_obj, trim(paramnames_megan2(j)), &
                                       MEGAN_PARAM(j,xknt))
                     ELSE IF (xactive_algorithm == 'MEGAN3' ) THEN
                         call read_data(ecfile_obj, trim(paramnames_megan3(j)), &
                                       MEGAN_PARAM(j,xknt))
                     ENDIF
               ENDIF
            ENDDO ! j/ nparams
            call close_file(ecfile_obj)
         ENDIF !/if do_AM3_ISOP
      ENDIF ! has_xactive
   ENDDO ! i/ species

!----------------------------------------------------------------------------
!  ... Set up the data for Temperature, Downward Shortwave Radiation,
!      LAI, PFT, and optionally - winds, soil moisture, ozone, and CO2.
!      If reproducing AM3 ISOP emis, the files used cannot be modified.
!----------------------------------------------------------------------------

!----------------------------------------------------------------------------
!  ...  Surface temperature setup
!----------------------------------------------------------------------------
   ALLOCATE( diag_gamma_temp(nlon,nlat) )
   diag_gamma_temp(:,:) = 0.
   IF ( do_ONLINE_TEMP ) THEN
      ALLOCATE( T24_STORE(nlon,nlat,24) )
      T24_STORE(:,:,:) = 0.
! Currently, the monthly average surface temperature must still be read in
! since it is still required by "fGAMMA_AGE", and until it is coded to be
! read from a restart file, the array is too large to be justifiably
! created for it's minimum impact.
      !IF ( do_AM3_ISOP ) THEN
         ALLOCATE( Tmo(nlon,nlat,nMOS) )
         call temp_init_AM3( lonb, latb, axes)
      !ENDIF
   ELSE
      IF ( do_GAMMA_HT .or. do_GAMMA_LT ) THEN
         call error_mesg( 'xactive_bvoc_init', &
          'High/Low temperature gammas must use online temperatures', WARNING)
         do_ONLINE_TEMP = .TRUE.
         ALLOCATE( T24_STORE(nlon,nlat,24) )
         T24_STORE(:,:,:) = 0.
      ENDIF
      ALLOCATE( Tmo(nlon,nlat,nMOS) )
      call temp_init_AM3( lonb, latb, axes)
   ENDIF

!---------------------------------------------------------------
!   ...  PAR/PPFD setup
!---------------------------------------------------------------
   ALLOCATE( diag_gamma_par(nlon,nlat) )
   diag_gamma_par(:,:) = 0.
   ALLOCATE( Pmo(nlon,nlat,nMOS ) )
   IF ( do_ONLINE_PPFD ) THEN
      ALLOCATE( P24_STORE(nlon,nlat,24) )
      P24_STORE(:,:,:) = 0.
      IF ( do_AM3_ISOP ) THEN
         call ppfd_init_AM3( lonb, latb, axes )
      ENDIF
   ELSE
      call ppfd_init_AM3(lonb, latb, axes)
   ENDIF

!-----------------------------------
!  ... Plant functional type setup
!-----------------------------------
   IF ( xactive_algorithm == 'MEGAN2' .OR. do_AM3_ISOP )  THEN
      ALLOCATE( diag_gamma_age(nlon,nlat,nPFT) )
      diag_gamma_age(:,:,:) = 0.
   ENDIF
   IF ( xactive_algorithm == 'MEGAN3' ) THEN
      ALLOCATE(diag_gamma_age_megan3(nlon,nlat) )
      diag_gamma_age_megan3(:,:) = 0.
   ENDIF
   IF ( do_ONLINE_PFT ) THEN
      IF ( xactive_algorithm == 'MEGAN2' ) THEN
         ALLOCATE( PCTPFT(nlon,nlat,nPFT) )
      ENDIF
      IF ( do_AM3_ISOP ) THEN
         call pft_init_AM3( lonb, latb, axes )
      ENDIF
   ELSE
      IF ( xactive_algorithm == 'MEGAN2' .OR. do_AM3_ISOP ) THEN
         ALLOCATE( PCTPFT(nlon,nlat,nPFT) )
         call pft_init_AM3 (lonb, latb, axes)
      ELSE IF ( xactive_algorithm == 'MEGAN3' ) THEN
         ! MEGAN3 does not require pctpft data
      ENDIF
   ENDIF

!---------------------
!  ... LAI setup
!---------------------
   IF ( xactive_algorithm == 'MEGAN2' .OR. do_AM3_ISOP ) THEN
      ALLOCATE( diag_gamma_lai(nlon,nlat,nPFT) )
      diag_gamma_lai(:,:,:) = 0.
   ENDIF
   IF ( xactive_algorithm == 'MEGAN3' ) THEN
      ALLOCATE( diag_gamma_lai_megan3(nlon,nlat) )
      diag_gamma_lai_megan3(:,:) = 0.
   ENDIF

   IF ( do_ONLINE_LAI ) THEN
      IF ( xactive_algorithm == 'MEGAN3' ) THEN
         ALLOCATE( MLAI_MEGAN3(nlon,nlat,2) )
         ALLOCATE( FCOVER(nlon,nlat,2) )
         call fcover_init_megan3 (lonb,latb,axes)
      ELSE IF ( xactive_algorithm == 'MEGAN2' ) THEN
         ALLOCATE( MLAI(nlon,nlat,nPFT,2) )
      ENDIF
      IF ( do_AM3_ISOP ) THEN
         call lai_init_AM3 (lonb, latb, axes)
      ENDIF
   ELSE
      IF ( xactive_algorithm == 'MEGAN2' .OR. do_AM3_ISOP ) THEN
         ALLOCATE( MLAI(nlon,nlat,nPFT,nMOS) )
         call lai_init_AM3( lonb, latb, axes)
      ENDIF
      IF ( xactive_algorithm == 'MEGAN3' ) THEN
         ALLOCATE ( MLAI_MEGAN3(nlon,nlat,nMOS) )
         MLAI_MEGAN3(:,:,:) = 0.
         call lai_init_megan3 (lonb, latb, axes )
      ENDIF
   ENDIF

!-----------------------
!  ... BDLAI DIAG
!----------------------
   IF ( do_GAMMA_BDLAI ) THEN
      ALLOCATE(diag_gamma_bdlai_megan3(nlon,nlat))
      diag_gamma_bdlai_megan3(:,:) = 0.
   ENDIF

!--------------------
!  ...  High winds
!--------------------
   IF ( do_GAMMA_HW ) THEN
      ALLOCATE( diag_gamma_hw(nlon,nlat ))
      diag_gamma_hw(:,:) = 0.
      ALLOCATE( WS_STORE(nlon,nlat,24) )
      WS_STORE(:,:,:) = 0.
      IF ( .not. do_ONLINE_WIND ) THEN
         call error_mesg( 'xactive_bvoc_init', &
         'No wind speed file available, must be calculated online', WARNING)
         do_ONLINE_WIND = .TRUE.
      ENDIF
   ENDIF

!-----------------------
!  ... Soil moisture
!-----------------------
   IF ( do_GAMMA_SM ) THEN
      ALLOCATE(diag_gamma_sm(nlon,nlat))
      diag_gamma_sm(:,:) = 0.
      IF (.not. do_ONLINE_SM ) THEN
         call error_mesg( 'xactive_bvoc_init', &
         'Offline soil moisture gamma N/A, must be online',WARNING)
         do_ONLINE_SM = .TRUE.
      ENDIF
   ENDIF


!-----------------------
!  ... CO2
!-----------------------
   IF ( do_GAMMA_CO2 ) THEN
      ALLOCATE(diag_gamma_co2(nlon,nlat) )
      diag_gamma_co2(:,:) = 0.
      ALLOCATE(CO2_STORE(nlon,nlat))
      CO2_STORE(:,:) = 0.
      IF (.not. do_ONLINE_CO2 ) THEN
         call error_mesg('xactive_bvoc_init', &
         'Offline CO2 N/A, must be online',WARNING)
         do_ONLINE_CO2 = .TRUE.
      ENDIF
   ENDIF

!-----------------------
!  ... Ozone
!-----------------------
   IF ( do_GAMMA_AQ ) THEN
      ALLOCATE(O3_STORE(nlon,nlat,24))
      O3_STORE(:,:,:) = 0.
      ALLOCATE(diag_gamma_aq(nlon,nlat))
      diag_gamma_aq(:,:) = 0.
      IF (.not. do_ONLINE_O3 ) THEN
         call error_mesg('xactive_bvoc_init', &
         'Offline O3 n/a, must be online', WARNING)
         do_ONLINE_O3 = .TRUE.
      ENDIF
   ENDIF

   xactive_domain => domain
   if (Ldebug .and. mpp_pe()==mpp_root_pe()) &
      write(*,*) 'xactive_bvoc_init: calling xactive_bvoc_register_restart'

   !< Get the current pelist
   allocate(pes(mpp_npes()))
   call mpp_get_current_pelist(pes)

   !< Open the scalar file with the current pelist, so that only the root pe opens and reads the file and
   !! distributes the data to the other pes
   if (open_file(Xbvoc_restart,"INPUT/xactive_bvoc.res.nc","read", is_restart=.true., pelist=pes)) then
      if (mpp_pe() == mpp_root_pe() ) &
      call error_mesg ('xactive_bvoc_mod',  'xactive_bvoc_init:&
         &Reading netCDF formatted restart file: xactive_bvoc.res.nc', NOTE)
      call xactive_bvoc_register_restart_scalars(Xbvoc_restart)
      call read_restart(Xbvoc_restart)
      call close_file(Xbvoc_restart)
   endif
   deallocate(pes)

   if (open_file(Til_restart,"INPUT/xactive_bvoc.res.nc","read", xactive_domain, is_restart=.true.)) then
      call xactive_bvoc_register_restart_domains(Til_restart)
      call read_restart(Til_restart)
      call close_file(Til_restart)
   endif

   module_is_initialized = .TRUE.

   IF (mpp_pe() == mpp_root_pe()) THEN
      write(*,*) ' xactive BVOC module is initialized...'
      IF ( xactive_algorithm == 'MEGAN2' ) THEN
         write(*,*) 'Using the algorithms and emission capacities from Megan v2.1'
      ELSE IF ( xactive_algorithm == 'MEGAN3') THEN
         write(*,*) 'Using the algorithms and emission capacities from Megan v3.0'
      ENDIF
      IF ( do_AM3_ISOP ) write(*,*) 'Reproducing AM3 legacy isoprene emissions'
   ENDIF

end subroutine xactive_bvoc_init

!----------------------
!-------------------


!########################################################################################
!
!<SUBROUTINE NAME="calc_xactive_bvoc_AM3">
!  <OVERVIEW>
!    Calculates interactive biogenic isoprene emissions following the algorithms used in AM3
!    as implemented by Arlene M. Fiore and Vaishali A. Naik.
!  </OVERVIEW
!  <DESCRIPTION>
!     Calculates interactive BVOC emissions using algorithms from
!     PCEEA MEGAN model in Guenther, ACP, 2006.
!     Note - gamma soil moisture is assumed constant (at one)
!  </DESCRIPTION>
!  <TEMPLATE>
!    call calc_xactive_bvoc_AM3 ( Time, Time_next, is, js, lon, lat, land, coszen,  &
!                                 PPFD1, T1, LAIp, LAIc, Pclim, Tclim,              &
!                                 ECBVOC_S, month, EMIS,                            &
!                                 id_GAMMA_TEMP, id_GAMMA_PAR, id_GAMMA_LAI,        &
!                                 id_GAMMA_AGE )
!  </TEMPLATE>
!  <IN NAME="Time, Time_next" TYPE="type(time_type)">
!    Model time
!  </IN>
!  <IN NAME="is, js" TYPE="integer">
!    Local domain start indices
!  </IN>
!  <IN NAME="lon,lat" TYPE="real" DIM="(:,:)">
!    Longitude/Latitude centers of the local domain
!  </IN>
!  <IN NAME="land" TYPE="real" DIM="(:,:)">
!    Land fraction
!  </IN>
!  <IN NAME="coszen" TYPE="real" DIM="(:,:)">
!    Cosine of the solar zenith angle
!  </IN>
!  <IN NAME="PPFD1" TYPE="real" DIM="(:,:)">
!    Current photosynthetic photon flux density
!  </IN>
!  <IN NAME="T1" TYPE="real" DIM="(:,:)">
!    Current temperature
!  </IN>
!  <IN NAME="LAIp" TYPE="real" DIM="(:,:)">
!    Previous month's LAI
!  </IN>
!  <IN NAME="LAIc" TYPE="real" DIM="(:,:)">
!    Current month's LAI
!  </IN>
!  <IN NAME="Tclim" TYPE="real" DIM="(:,:)">
!    Climatological temperature
!  </IN>
!  <IN NAME="Pclim" TYPE="real" DIM="(:,:)">
!    Climatological PPFD
!  </IN>
!  <IN NAME="ECBVOC_S" TYPE="real" DIM="(:,:,:)">
!    Isoprene Emission capacities for each vegetation type
!  </IN>
!  <IN NAME="month" TYPE="integer">
!    Current month
!  </IN>
!  <IN NAME="species" TYPE="character">
!    Species name
!  </IN>
!  <OUT NAME="EMIS" TYPE="real" DIM="(:,:)">
!    Output emissions for this timestep
!  </OUT>
!  <IN NAME="id_*" TYPE="integer, optional">
!    IDs for diagnotiscs (gammas and emissions)
!  </IN>
!
subroutine calc_xactive_bvoc_AM3( Time, Time_next, is, js, lon, lat, land, coszen,    &
                                  PPFD1, T1, LAIp, LAIc, Pclim, Tclim,                &
                                  ECBVOC_S, month, EMIS,       &
                                  id_GAMMA_TEMP, id_GAMMA_PAR, id_GAMMA_LAI,          &
                                  id_GAMMA_AGE )

   type(time_type), intent(in)            :: Time, Time_next
   integer, intent(in)                    :: is, js
   real, intent(in), dimension(:,:)       :: lon, lat
   real, intent(in), dimension(:,:)       :: land
   real, intent(in), dimension(:,:)       :: coszen
   real, intent(in), dimension(:,:)       :: PPFD1
   real, intent(in), dimension(:,:)       :: T1
   real, intent(in), dimension(:,:,:)     :: LAIp, LAIc
   real, intent(in), dimension(:,:,:)     :: Pclim, Tclim
   real, intent(in), dimension(:,:,:)     :: ECBVOC_S
   integer, intent(in)                    :: month

   real, intent(out), dimension(:,:)      :: EMIS

   integer, intent(in), optional          :: id_GAMMA_TEMP, id_GAMMA_PAR
   integer, intent(in), optional          :: id_GAMMA_LAI, id_GAMMA_AGE

!------------ Local Variables ------------------------------------------------------------
   type(time_type) :: Year_t
   real     :: BETA, LDF, C_t1, C_eo, A_new, A_gro, A_mat, A_old, MW_sp
   real     :: GAMMA_TLD, GAMMA_PAR,  calday
   real     :: GAMMA_AGE(nPFT), GAMMA_LAI(nPFT), GAMMA_LAIAGE(nPFT)
   logical  :: age_work(nPFT)
   logical  :: do_age(nPFT)
   logical  :: used
   real     :: work_emis(nVEG)
   integer  :: pft_li(nVEG)
   integer  :: pft_lu(nVEG)
   integer  :: i, j, n, nlat, nlon, nl, nu, ie, je
   integer  :: yr, mo, day, hr, minute, sec

   ie = is + size(lon,1) -1
   je = js + size(lon,2) -1

   nlon = size(lon,1)
   nlat = size(lon,2)

   BETA  = 1.                 ! See Table 4, Guenther et al., 2012
   LDF   = 1.                 ! ------------------------
   C_t1  = 80.                ! for description of parameters
   C_eo  = 1.75
   A_new = 0.05
   A_gro = 0.6
   A_mat = 1.125
   A_old = 1.0
   MW_sp = 68.           !Molecular weight for conversion from ug to molecules

   do_age(:) = .TRUE.               !Calculate gamma age for this PFTs?
   do_age((/2,3,5,6,10/)) = .FALSE. !gamma LAI = 1 for evergreen PFTs

   pft_li(:) = (/ 2,5,10,13,16 /)
   pft_lu(1:nVEG-1) = (/ 4,9,12,15 /)
   pft_lu(nVEG) = nPFT

   call get_date(Time,yr,mo,day,hr,minute,sec)  !model GMT
   !Get Julian date (fraction) = calday
   Year_t = set_date(yr,1,1,0,0,0)
   calday = time_type_to_real( Time-Year_t) / SECONDS_PER_DAY

   EMIS(:,:) = 0.
   DO j = 1, nlat
   DO i = 1, nlon
      IF ( land(i,j) > min_land_frac ) THEN

!------------------------------------------------------------------------------------
!                     ... GAMMA LIGHT CALCULATION
!------------------------------------------------------------------------------------
        IF ( coszen(i,j) <= 0. ) THEN
           GAMMA_PAR  = 0.
        ELSE
           GAMMA_PAR = fGAMMA_PAR_AM3(PPFD1(i,j), coszen(i,j), Pclim(i,j,month), calday)
        ENDIF
        GAMMA_PAR = max(GAMMA_PAR, 0.)
!------------------------------------------------------------------------------------
!                     ... GAMMA TEMP CALCULATION
!------------------------------------------------------------------------------------
        GAMMA_TLD = fGAMMA_TLD_AM3(T1(i,j), Tclim(i,j,month), C_t1, C_eo)
        GAMMA_TLD = max(GAMMA_TLD, 0.)
!------------------------------------------------------------------------------------
!                     ... GAMMA AGE CALCULATION
!------------------------------------------------------------------------------------
        GAMMA_AGE(:) = fGAMMA_AGE_MEGAN2(LAIp(i,j,:), LAIc(i,j,:), Tclim(i,j,month), &
                                     month, A_new, A_gro, A_mat, A_old, do_age)
!------------------------------------------------------------------------------------
!                     ... GAMMA LAI CALCULATION
!------------------------------------------------------------------------------------
        GAMMA_LAI(:) = 0.49 * LAIc(i,j,:) /  &
                       sqrt( 1. + 0.2 * LAIc(i,j,:)*LAIc(i,j,:))
!------------------------------------------------------------------------------------
!                     ... Combine LAI and AGE
!------------------------------------------------------------------------------------
        GAMMA_LAIAGE(:) = GAMMA_LAI * GAMMA_AGE

!------------------------------------------------------------------------------------
!                     ... Sum over the PFTs part of each VEG
!------------------------------------------------------------------------------------
        DO n = 1, nVEG
           nl = pft_li(n)
           nu = pft_lu(n)
           work_emis(n) = dot_product( GAMMA_LAIAGE(nl:nu), PCTPFT(i+is-1,j+js-1,nl:nu)) &
                          * ECBVOC_S(i,j,n)
        ENDDO

           EMIS(i,j) = sum(work_emis) * GAMMA_TLD * GAMMA_PAR
! Update diagnostics
           diag_gamma_temp(i+is-1,j+js-1)  = GAMMA_TLD
           diag_gamma_par(i+is-1,j+js-1)   = GAMMA_PAR
           diag_gamma_lai(i+is-1,j+js-1,:) = GAMMA_LAI(:)
           diag_gamma_age(i+is-1,j+js-1,:) = GAMMA_AGE(:)

      ENDIF!land
   ENDDO
   ENDDO

! Apply the canopy loss factor and convert from ug/m2/hr to molecules/m2/s
     EMIS(:,:) = EMIS(:,:) * RHO_CANOPY  *  (1.67e10 / MW_sp)

!---------------------------------------------------------------------------------------
!              ... Send diagnostics
!---------------------------------------------------------------------------------------

! Gamma temperature
      IF (present(id_GAMMA_TEMP) .AND. id_GAMMA_TEMP > 0) THEN
         used = send_data( id_GAMMA_TEMP, diag_gamma_temp(is:ie,js:je), &
                           Time_next, is_in=is, js_in=js )
      ENDIF
! Gamma light
      IF (present(id_GAMMA_PAR) .AND. id_GAMMA_PAR > 0) THEN
         used = send_data( id_GAMMA_PAR, diag_gamma_par(is:ie,js:je), &
                           Time_next, is_in=is, js_in=js )
      ENDIF
! Gamma LAI
      IF (present(id_GAMMA_LAI) .AND. id_GAMMA_LAI > 0) THEN
         used = send_data( id_GAMMA_LAI, diag_gamma_lai(is:ie,js:je,:), &
                           Time_next, is_in=is, js_in=js )
      ENDIF
! Gamma Age
      IF (present(id_GAMMA_AGE) .AND. id_GAMMA_AGE > 0) THEN
         used = send_data( id_GAMMA_AGE, diag_gamma_age(is:ie,js:je,:), &
                           Time_next, is_in=is, js_in=js )
      ENDIF

end subroutine calc_xactive_bvoc_AM3
!</SUBROUTINE>


!########################################################################################
!<SUBROUTINE NAME="calc_xactive_bvoc_megan2">
!  <OVERVIEW>
!    Calculates interactive biogenic emissions
!  </OVERVIEW
!  <DESCRIPTION>
!    Calculates interactive BVOC (and CO) using algorithms from MEGAN v3
!     documented in Guenther et al. (2018).  Each gamma has a
!     seperate function and only called if the gamma is required (e.g., PAR) or
!    - if it is an optional gamma - the namelist value is set, otherwise the gamma
!    is set =1. The gammas are applied to the emissions factors for each PFT/VEG type.
!    Diagnostics are accumulated and sent.
!  </DESCRIPTION>
!  <TEMPLATE>
!    call calc_xactive_bvoc_megan2 ( Time, Time_next, is, js, ,lon, lat, land,  &
!                                 coszen,  PPFD1, PPFD24, T1, T24, LAIp, LAIc,  &
!                                 Tclim,                                        &
!                                 MEGAN_PARAM, ECBVOC_S, month, species, EMIS,  &
!                                 id_GAMMA_TEMP, id_GAMMA_PAR, id_GAMMA_LAI,    &
!                                 id_GAMMA_AGE, id_GAMMA_CO2,  id_GAMMA_SM )
!  </TEMPLATE>
!  <IN NAME="Time, Time_next" TYPE="type(time_type)">
!    Model time
!  </IN>
!  <IN NAME="is, js" TYPE="integer">
!    Local domain start indices
!  </IN>
!  <IN NAME="lon,lat" TYPE="real" DIM="(:,:)">
!    Longitude/Latitude centers of the local domain
!  </IN>
!  <IN NAME="land" TYPE="real" DIM="(:,:)">
!    Land fraction
!  </IN>
!  <IN NAME="coszen" TYPE="real" DIM="(:,:)">
!    Cosine of the solar zenith angle
!  </IN>
!  <IN NAME="PPFD1" TYPE="real" DIM="(:,:)">
!    Instantaneous photosynthetic photon flux density
!  </IN>
!  <IN NAME="PPFD24" TYPE="real" DIM="(:,:)">
!    Previous 24h average PPFD
!  </IN>
!  <IN NAME="T1" TYPE="real" DIM="(:,:)">
!    Current temperature
!  </IN>
!  <IN NAME="T24" TYPE="real" DIM="(:,:)">
!    Previous 24h avg temperature
!  </IN>
!  <IN NAME="LAIp" TYPE="real" DIM="(:,:,:)">
!    Previous month's LAI
!  </IN>
!  <IN NAME="LAIc" TYPE="real" DIM="(:,:,:)">
!    Current month's LAI
!  </IN>
!  <IN NAME="Tclim" TYPE="real" DIM="(:,:)">
!    Climatological temperature
!  </IN>
!  <IN NAME="MEGAN_PARAM" TYPE="real" DIM="(:,:)">
!    Megan model parameters for this species
!  </IN>
!  <IN NAME="ECBVOC_S" TYPE="real" DIM="(:,:)">
!    Emission capacities for each pft/veg type for this species
!  </IN>
!  <IN NAME="month" TYPE="integer">
!    Current month
!  </IN>
!  <IN NAME="species" TYPE="character">
!    Species name
!  </IN>
!  <OUT NAME="EMIS" TYPE="real" DIM="(:,:)">
!    Output emissions for this timestep
!  </OUT>
!  <IN NAME="id_*" TYPE="integer, optional">
!    IDs for diagnotiscs (gammas and emissions)
!  </IN>
!
subroutine calc_xactive_bvoc_megan2 ( Time, Time_next, is, js, lon, lat, land,     &
                                      coszen,  PPFD1, PPFD24, T1, T24, LAIp, LAIc, &
                                      Tclim,                                       &
                                      MEGAN_PARAM, ECBVOC_S, month, species, EMIS, &
                                      id_GAMMA_TEMP, id_GAMMA_PAR, id_GAMMA_LAI,   &
                                      id_GAMMA_AGE, id_GAMMA_CO2, id_GAMMA_SM )


      type(time_type), intent(in)        :: Time, Time_next
      integer, intent(in)                :: is, js          ! Local domain start indices
      real, intent(in), dimension(:,:)   :: lon, lat        ! Lat/lon centers
      real, intent(in), dimension(:,:)   :: land            ! Land Fraction
      real, intent(in), dimension(:,:)   :: coszen          ! cosine zenith angle
      real, intent(in), dimension(:,:)   :: PPFD1           ! PPPFD, this timestep [umoles/m2/s]
      real, intent(in), dimension(:,:)   :: PPFD24          ! PPFD, daily avg.    [umoles/m2/s]
      real, intent(in), dimension(:,:)   :: T1              ! Surf T, this timestep, [K]
      real, intent(in), dimension(:,:)   :: T24             ! Surf T (daily avg.)    [K]
      real, intent(in), dimension(:,:,:) :: LAIp, LAIc      ! Previous & current LAI [m2/m2]
      real, intent(in), dimension(:,:,:) :: Tclim           ! Climatological temperature [K]
      real, intent(in), dimension(:)     :: MEGAN_PARAM     ! MEGAN model parameters
      real, intent(in), dimension(:,:,:) :: ECBVOC_S        ! Emission capacities  [ug/m2/h]
      integer, intent(in)                :: month           ! Current month index
      character(len=*), intent(in)       :: species         ! species name

      real, intent(out), dimension(:,:)  :: EMIS            ! Emissions for this timestep [molec/m2/s]

! Optional diagnostic IDs for emissions and gammas
      integer, intent(in), optional      :: id_GAMMA_PAR, id_GAMMA_TEMP
      integer, intent(in), optional      :: id_GAMMA_LAI, id_GAMMA_AGE
      integer, intent(in), optional      :: id_GAMMA_SM, id_GAMMA_CO2
!------------------------------------------------------------------------------------
!                 ... local variables
!------------------------------------------------------------------------------------

      real     :: BETA, LDF, C_t1, C_eo           ! temperature emission factors
      real     :: A_new, A_gro, A_mat, A_old      ! Leaf age factors
      real     :: MW_sp                           ! molecular weight
      real     :: GAMMA_TLD, GAMMA_TLI, GAMMA_TMP
      real     :: GAMMA_PAR
      real     :: GAMMA_LAI(nPFT), GAMMA_AGE(nPFT)
      real     :: GAMMA_SM, GAMMA_CO2
      real     :: dummy
      logical  :: do_age(nPFT)                    ! flag: calc gamma age
      logical  :: used
      integer  :: i, j, n, nlat, nlon, ie, je

      dummy = 1.

      nlon = size(lon,1)
      nlat = size(lon,2)

      ie = is + size(lon,1) -1
      je = js + size(lon,2) -1

      BETA  = MEGAN_PARAM(1)
      LDF   = MEGAN_PARAM(2)
      C_t1  = MEGAN_PARAM(3)
      C_eo  = MEGAN_PARAM(4)
      A_new = MEGAN_PARAM(5)
      A_gro = MEGAN_PARAM(6)
      A_mat = MEGAN_PARAM(7)
      A_old = MEGAN_PARAM(8)
      MW_sp = MEGAN_PARAM(9)

      do_age(:) = .TRUE.               !Calculate gamma age for this PFTs?
      do_age((/2,3,5,6,10/)) = .FALSE. !gamma LAI = 1 for evergreen PFTs

      EMIS(:,:) = 0.
      DO i = 1,nlon
      DO j = 1,nlat
! only do land
      IF ( land(i,j) > min_land_frac ) THEN
!------------------------------------------------------------------------------------
!                     ... GAMMA LIGHT CALCULATION
!------------------------------------------------------------------------------------
         IF ( coszen(i,j) <= 0. ) THEN
            GAMMA_PAR  = 0.
         ELSE
            GAMMA_PAR = fGAMMA_PAR_AM4(PPFD1(i,j), PPFD24(i,j))
         ENDIF
! Apply light dependent factor
        GAMMA_PAR = LDF*GAMMA_PAR + (1. - LDF)
! Prevent negative values
        GAMMA_PAR = max(GAMMA_PAR,0.)
!------------------------------------------------------------------------------------
!                     ...  GAMMA TEMPERATURE CALCULATION
!------------------------------------------------------------------------------------
! Light dependent
         GAMMA_TLD = fGAMMA_TLD_AM4(T1(i,j), T24(i,j), C_t1, C_eo)
! Light independent
         GAMMA_TLI = exp(BETA * (T1(i,j) - T_s))
! Combined
         GAMMA_TMP = LDF*GAMMA_TLD + (1. - LDF)*GAMMA_TLI
! Prevent negative
         GAMMA_TMP = max(GAMMA_TMP,0.)

! Soil moisture (According to Alex Guenther via Collete Heald, thoughts are
! mixed on whether soil moisture gamma should be used,)
         IF ( do_GAMMA_SM ) THEN
            GAMMA_SM = fGAMMA_SM(dummy)
         ELSE
            GAMMA_SM = 1.
         ENDIF

! CO2 (only for isoprene)
         IF ( do_GAMMA_CO2 .AND. trim(species) == 'ISOP' ) THEN
            GAMMA_CO2 = fGAMMA_CO2(CO2_STORE(i+is-1,j+js-1))
         ELSE
            GAMMA_CO2 = 1.
         ENDIF

! Gamma age
         GAMMA_AGE(:) = fGAMMA_AGE_MEGAN2(LAIp(i,j,:), LAIc(i,j,:), Tclim(i,j,month), &
                                     month, A_new, A_gro, A_mat, A_old, do_age)
! Gamma LAI
         GAMMA_LAI(:) = 0.49 * LAIc(i,j,:) / &
                          sqrt(1. + 0.2 * LAIc(i,j,:)*LAIc(i,j,:))
! APPLY THE GAMMAS !
         DO n = 1, nPFT
            EMIS(i,j) = EMIS(i,j)      + (ECBVOC_S(i,j,n)  * PCTPFT(i+is-1,j+js-1,n) *   &
                        GAMMA_TMP      * GAMMA_PAR         * GAMMA_AGE(n)  *   &
                        GAMMA_LAI(n)   * GAMMA_SM          * GAMMA_CO2        )
         ENDDO !PFTs
!--------------------------------------------------------------------------------------------

! Update diagnostics
         diag_gamma_temp(i+is-1,j+js-1)  = GAMMA_TMP
         diag_gamma_par(i+is-1,j+js-1)   = GAMMA_PAR
         diag_gamma_lai(i+is-1,j+js-1,:) = GAMMA_LAI(:)
         diag_gamma_age(i+is-1,j+js-1,:) = GAMMA_AGE(:)
! Update (optional) diagnostics
         IF (present(id_GAMMA_SM) .AND. id_GAMMA_SM > 0) THEN
            diag_gamma_sm(i+is-1,j+js-1) = GAMMA_SM
         ENDIF
         IF (present(id_GAMMA_CO2) .AND. id_GAMMA_CO2 > 0) THEN
            diag_gamma_co2(i+is-1,j+js-1) = GAMMA_CO2
         ENDIF
      ENDIF !land
      ENDDO !lon
      ENDDO !lat

! Apply the canopy loss factor and convert from ug/m2/hr to molecules/m2/s
      EMIS(:,:) = EMIS(:,:) * RHO_CANOPY  *  (1.67e10 / MW_sp)
!---------------------------------------------------------------------------------------
!              ... Accumulate diagnostics
!---------------------------------------------------------------------------------------

! Gamma temperature
      IF (present(id_GAMMA_TEMP) .AND. id_GAMMA_TEMP > 0) THEN
         used = send_data( id_GAMMA_TEMP, diag_gamma_temp(is:ie,js:je), &
                           Time_next, is_in=is, js_in=js )
      ENDIF
! Gamma light
      IF (present(id_GAMMA_PAR) .AND. id_GAMMA_PAR > 0) THEN
         used = send_data( id_GAMMA_PAR, diag_gamma_par(is:ie,js:je), &
                           Time_next, is_in=is, js_in=js )
      ENDIF
! Gamma LAI
      IF (present(id_GAMMA_LAI) .AND. id_GAMMA_LAI > 0) THEN
         used = send_data( id_GAMMA_LAI, diag_gamma_lai(is:ie,js:je,:), &
                           Time_next, is_in=is, js_in=js )
      ENDIF
! Gamma Age
      IF (present(id_GAMMA_AGE) .AND. id_GAMMA_AGE > 0) THEN
         used = send_data( id_GAMMA_AGE, diag_gamma_age(is:ie,js:je,:), &
                           Time_next, is_in=is, js_in=js )
      ENDIF
! Gamma soil
      IF (present(id_GAMMA_SM) .AND. id_GAMMA_SM > 0) THEN
         used = send_data( id_GAMMA_SM, diag_gamma_sm(is:ie,js:je), &
                           Time_next, is_in=is, js_in=js)
      ENDIF
! Gamma co2
      IF (present(id_GAMMA_CO2) .AND. id_GAMMA_CO2 > 0) THEN
         used = send_data(id_GAMMA_CO2, diag_gamma_co2(is:ie,js:je), &
                          Time_next, is_in=is, js_in=js)
      ENDIF

end subroutine calc_xactive_bvoc_megan2
!</SUBROUTINE



!########################################################################################
!<SUBROUTINE NAME="calc_xactive_bvoc_megan3">
!  <OVERVIEW>
!    Calculates interactive biogenic emissions
!  </OVERVIEW
!  <DESCRIPTION>
!    Calculates interactive BVOC (and CO) using algorithms, emission capacities,
!    and light dependent fractions (LDF) from MEGANv3 documented in
!    Guenther et al. (2018).  Each gamma has a seperate function and only called if
!    the gamma is required (e.g., PAR) or - if it is an optional gamma - the namelist
!    value is set, otherwise the gamma is set =1. The gammas are applied to the
!    emissions factors for each VEG type. Diagnostics are accumulated and sent.
!  </DESCRIPTION>
!  <TEMPLATE>
!    call calc_xactive_bvoc_megan3 ( Time, Time_next, is, js, ,lon, lat, land,   &
!                                 coszen,  PPFD1, PPFD24, T1, T24, LAIp, LAIc,   &
!                                 TMAX, TMIN, Tclim, WSMAX, AQI,                 &
!                                 MEGAN_PARAM, LDFg_S, ECBVOC_S, month, species, &
!                                 EMIS, id_GAMMA_TEMP, id_GAMMA_PAR,             &
!                                 id_GAMMA_LAI, id_GAMMA_AGE, id_GAMMA_BDLAI,    &
!                                 id_GAMMA_CO2, id_GAMMA_AQ, id_GAMMA_SM,        &
!                                 id_GAMMA_HT, id_GAMMA_LT, id_GAMMA_HW)
!  </TEMPLATE>
!  <IN NAME="Time, Time_next" TYPE="type(time_type)">
!    Model time
!  </IN>
!  <IN NAME="is, js" TYPE="integer">
!    Local domain start indices
!  </IN>
!  <IN NAME="lon,lat" TYPE="real" DIM="(:,:)">
!    Longitude/Latitude centers of the local domain
!  </IN>
!  <IN NAME="land" TYPE="real" DIM="(:,:)">
!    Land fraction
!  </IN>
!  <IN NAME="coszen" TYPE="real" DIM="(:,:)">
!    Cosine of the solar zenith angle
!  </IN>
!  <IN NAME="PPFD1" TYPE="real" DIM="(:,:)">
!    Instantaneous photosynthetic photon flux density
!  </IN>
!  <IN NAME="PPFD24" TYPE="real" DIM="(:,:)">
!    Previous 24h average PPFD
!  </IN>
!  <IN NAME="T1" TYPE="real" DIM="(:,:)">
!    Current temperature
!  </IN>
!  <IN NAME="T24" TYPE="real" DIM="(:,:)">
!    Previous 24h avg temperature
!  </IN>
!  <IN NAME="LAIp" TYPE="real" DIM="(:,:)">
!    Previous month's LAI
!  </IN>
!  <IN NAME="LAIc" TYPE="real" DIM="(:,:)">
!    Current month's LAI
!  </IN>
!  <IN NAME="TMAX" TYPE="real" DIM="(:,:)">
!    Daily max temperature
!  </IN>
!  <IN NAME="TMIN" TYPE="real" DIM="(:,:)">
!    Daily min temperature
!  </IN>
!  <IN NAME="Tclim" TYPE="real" DIM="(:,:)">
!    Climatological temperature
!  </IN>
!  <IN NAME="WSMAX" TYPE="real" DIM="(:,:)">
!    Daily max 10-m wind speed
!  </IN>
!  <IN NAME="AQI" TYPE="real" DIM="(:,:)">
!    Aiq quality index
!  </IN>
!  <IN NAME="MEGAN_PARAM" TYPE="real" DIM="(:,:)">
!    Megan model parameters for this species
!  </IN>
!  <IN NAME="LDFg_S" TYPE="real" DIM="(:,:)">
!    Gridded LDF fractions
!  </IN>
!  <IN NAME="ECBVOC_S" TYPE="real" DIM="(:,:)">
!    Emission capacities for each pft/veg type for this species
!  </IN>
!  <IN NAME="month" TYPE="integer">
!    Current month
!  </IN>
!  <IN NAME="species" TYPE="character">
!    Species name
!  </IN>
!  <OUT NAME="EMIS" TYPE="real" DIM="(:,:)">
!    Output emissions for this timestep
!  </OUT>
!  <IN NAME="id_*" TYPE="integer, optional">
!    IDs for diagnotiscs (gammas)
!  </IN>
!
subroutine calc_xactive_bvoc_megan3 ( Time, Time_next, is, js, lon, lat, land,     &
                                      coszen,  PPFD1, PPFD24, T1, T24, LAIp, LAIc, &
                                      TMAX, TMIN, Tclim, WSMAX, AQI,               &
                                      MEGAN_PARAM, LDFg_S, ECBVOC_S, month,        &
                                      species, EMIS,                               &
                                      id_GAMMA_TEMP, id_GAMMA_PAR, id_GAMMA_LAI,   &
                                      id_GAMMA_AGE, id_GAMMA_BDLAI, id_GAMMA_CO2,  &
                                      id_GAMMA_AQ, id_GAMMA_SM, id_GAMMA_HT,       &
                                      id_GAMMA_LT, id_GAMMA_HW )


      type(time_type), intent(in)        :: Time, Time_next
      integer, intent(in)                :: is, js          ! Local domain start indices
      real, intent(in), dimension(:,:)   :: lon, lat        ! Lat/lon centers
      real, intent(in), dimension(:,:)   :: land
      real, intent(in), dimension(:,:)   :: coszen          ! cosine zenith angle
      real, intent(in), dimension(:,:)   :: PPFD1           ! PPPFD, this timestep [umoles/m2/s]
      real, intent(in), dimension(:,:)   :: PPFD24          ! PPFD, daily avg.    [umoles/m2/s]
      real, intent(in), dimension(:,:)   :: T1              ! Surf T, this timestep, [K]
      real, intent(in), dimension(:,:)   :: T24             ! Surf T (daily avg.)    [K]
      real, intent(in), dimension(:,:)   :: LAIp, LAIc      ! Previous & current LAIv [m2/m2]
      real, intent(in), dimension(:,:)   :: TMAX, TMIN      ! Daily max/min temperature [K]
      real, intent(in), dimension(:,:,:) :: Tclim           ! Climatological temperature [K]
      real, intent(in), dimension(:,:)   :: WSMAX           ! Daily max 10-m wind speed [m/s]
      real, intent(in), dimension(:,:)   :: AQI             ! Air quality index
      real, intent(in), dimension(:)     :: MEGAN_PARAM     ! MEGAN model parameters
      real, intent(in), dimension(:,:)   :: ECBVOC_S        ! Emission capacities  [ug/m2/h]
      real, intent(in), dimension(:,:)   :: LDFg_S          ! Light dependent fraction
      integer, intent(in)                :: month           ! Current month index
      character(len=*), intent(in)       :: species         ! species name

      real, intent(out), dimension(:,:)  :: EMIS            ! Emissions for this timestep [molec/m2/s]

! Optional diagnostic IDs for emissions and gammas
      integer, intent(in), optional      :: id_GAMMA_PAR, id_GAMMA_TEMP
      integer, intent(in), optional      :: id_GAMMA_LAI, id_GAMMA_AGE
      integer, intent(in), optional      :: id_GAMMA_BDLAI, id_GAMMA_CO2
      integer, intent(in), optional      :: id_GAMMA_AQ, id_GAMMA_SM, id_GAMMA_HT
      integer, intent(in), optional      :: id_GAMMA_LT, id_GAMMA_HW

!------------------------------------------------------------------------------------
!                 ... local variables
!------------------------------------------------------------------------------------

      real     :: BETA, C_t1, C_eo                ! temperature emission factors
      real     :: A_new, A_gro, A_mat, A_old      ! Leaf age factors
      real     :: T_HT, DT_HT, C_HT               ! High temperature parameters
      real     :: T_LT, DT_LT, C_LT               ! Low temperature ...
      real     :: T_HW, DT_HW, C_HW               ! High wind ...
      real     :: T_AQ, DT_AQ, C_AQ               ! Air quality ...
      real     :: MW_sp                           ! molecular weight
      real     :: GAMMA_TLD, GAMMA_TLI, GAMMA_TMP
      real     :: GAMMA_PAR
      real     :: GAMMA_LAI, GAMMA_AGE
      real     :: GAMMA_HT, GAMMA_LT, GAMMA_SM
      real     :: GAMMA_AQ, GAMMA_HW, GAMMA_CO2
      real     :: GAMMA_BDLAI
      real     :: dummy
      logical  :: used
      integer  :: i, j, n, nlat, nlon, ie, je

      dummy = 1. ! Something to send to soil mosisture gamma

      nlon = size(lon,1)
      nlat = size(lon,2)

      ie = is + size(lon,1) -1
      je = js + size(lon,2) -1

      BETA  = MEGAN_PARAM(1)
      C_t1  = MEGAN_PARAM(2)
      C_eo  = MEGAN_PARAM(3)
      A_new = MEGAN_PARAM(4)
      A_gro = MEGAN_PARAM(5)
      A_mat = MEGAN_PARAM(6)
      A_old = MEGAN_PARAM(7)
      C_AQ  = MEGAN_PARAM(8)
      C_HW  = MEGAN_PARAM(9)
      C_HT  = MEGAN_PARAM(10)
      C_LT  = MEGAN_PARAM(11)
      T_AQ  = MEGAN_PARAM(12)
      T_HW  = MEGAN_PARAM(13)
      T_HT  = MEGAN_PARAM(14)
      T_LT  = MEGAN_PARAM(15)
      DT_AQ = MEGAN_PARAM(16)
      DT_HW = MEGAN_PARAM(17)
      DT_HT = MEGAN_PARAM(18)
      DT_LT = MEGAN_PARAM(19)
      MW_sp = MEGAN_PARAM(20)

      EMIS(:,:) = 0.0
      DO i = 1,nlon
      DO j = 1,nlat
      IF ( land(i,j) > min_land_frac ) THEN
!------------------------------------------------------------------------------------
!                     ... GAMMA LIGHT CALCULATION
!------------------------------------------------------------------------------------
         IF ( coszen(i,j) <= 0. ) THEN
            GAMMA_PAR  = 0.
         ELSE
            GAMMA_PAR = fGAMMA_PAR_AM4(PPFD1(i,j), PPFD24(i,j))
         ENDIF
! Apply light dependent factor
        GAMMA_PAR = LDFg_S(i,j)*GAMMA_PAR + (1. - LDFg_S(i,j))
! Prevent negative value
        GAMMA_PAR = MAX(GAMMA_PAR, 0.)
!------------------------------------------------------------------------------------
!                     ...  GAMMA TEMPERATURE CALCULATION
!------------------------------------------------------------------------------------
! Light dependent
        GAMMA_TLD = fGAMMA_TLD_AM4(T1(i,j), T24(i,j), C_t1, C_eo)
! Light independent
        GAMMA_TLI = exp(BETA * (T1(i,j) - T_s))
! Combined
        GAMMA_TMP = LDFg_S(i,j)*GAMMA_TLD + (1. - LDFg_S(i,j))*GAMMA_TLI
! Prvent negative
        GAMMA_TMP = MAX(GAMMA_TMP, 0.)
!------------------------------------------------------------------------------------
!                    ...  Optional Gammas
!------------------------------------------------------------------------------------

! High temperature
         IF ( do_GAMMA_HT ) THEN
             GAMMA_HT = fGAMMA_HT(TMAX(i,j), T_HT, DT_HT, C_HT)
         ELSE
             GAMMA_HT = 1.
         ENDIF

! Low Temperature
         IF ( do_GAMMA_LT ) THEN
            GAMMA_LT = fGAMMA_LT(TMIN(i,j), T_LT, DT_LT, C_LT)
         ELSE
            GAMMA_LT = 1.
         ENDIF

! High Winds
         IF ( do_GAMMA_HW ) THEN
            GAMMA_HW = fGAMMA_HW(WSMAX(i,j), T_HW, DT_HW, C_HW)
         ELSE
            GAMMA_HW = 1.
         ENDIF

! Air Quality
         IF ( do_GAMMA_AQ ) THEN
            GAMMA_AQ = fGAMMA_AQ(AQI(i,j), T_AQ, DT_AQ, C_AQ)
         ELSE
            GAMMA_AQ = 1.
         ENDIF

! Soil moisture (According to Alex Guenther via Collete Heald, thoughts are
! mixed on whether soil moisture gamma should be used,)
         IF ( do_GAMMA_SM ) THEN
            GAMMA_SM = fGAMMA_SM(dummy)
         ELSE
            GAMMA_SM = 1.
         ENDIF

! CO2 (only for isoprene)
         IF ( do_GAMMA_CO2 .AND. trim(species) == 'ISOP' ) THEN
            GAMMA_CO2 = fGAMMA_CO2(CO2_STORE(i+is-1,j+js-1))
         ELSE
            GAMMA_CO2 = 1.0
         ENDIF

!------------------------------------------------------------------------------------
!             ... PFT SPECIFIC GAMMA (AGE, LAI, bidirectional LAI) CALCULATIONS
!------------------------------------------------------------------------------------

! Bi-directional LAI (only for ethanol and acetylaldehyde in MEGAN3)
         IF ( do_GAMMA_BDLAI .AND. (trim(species) == 'C2H5OH'  &
              .OR. trim(species) == 'CH3CHO' )) THEN
            GAMMA_BDLAI = fGAMMA_BDLAI(LAIc(i,j))
         ELSE
            GAMMA_BDLAI = 1.
         ENDIF

! Gamma age
         GAMMA_AGE = fGAMMA_AGE_MEGAN3(LAIp(i,j), LAIc(i,j), Tclim(i,j,month), &
                                       month, A_new, A_gro, A_mat, A_old  )
! Gamma LAI
         !GAMMA_LAI = 0.49 * LAIc(i,j) / sqrt(1. + 0.2 * LAIc(i,j)*LAIc(i,j))
         GAMMA_LAI = LAIc(i,j)

!-------------------------------------------------------------------------------------------
! APPLY THE GAMMAS !
         EMIS(i,j) =  LDFg_S(i,j)     * ECBVOC_S(i,j)    * GAMMA_TMP    * GAMMA_PAR    *   &
                      GAMMA_AGE       * GAMMA_LAI        * GAMMA_BDLAI  *   &
                      GAMMA_SM        * GAMMA_AQ         * GAMMA_CO2    *   &
                      GAMMA_HT        * GAMMA_LT         * GAMMA_HW

!--------------------------------------------------------------------------------------------

! Update diagnostics
         diag_gamma_temp(i+is-1,j+js-1)  = GAMMA_TMP
         diag_gamma_par(i+is-1,j+js-1)   = GAMMA_PAR
         diag_gamma_lai_megan3(i+is-1,j+js-1) = GAMMA_LAI
         diag_gamma_age_megan3(i+is-1,j+js-1) = GAMMA_AGE
! Update (optional) diagnostics
         IF (present(id_GAMMA_BDLAI) .AND. id_GAMMA_BDLAI > 0 ) THEN
            diag_gamma_bdlai_megan3(i+is-1,j+js-1) = GAMMA_BDLAI
         ENDIF
         IF (present(id_GAMMA_HT) .AND. id_GAMMA_HT > 0) THEN
            diag_gamma_ht(i+is-1,j+js-1) = GAMMA_HT
         ENDIF
         IF (present(id_GAMMA_LT) .AND. id_GAMMA_LT > 0) THEN
            diag_gamma_lt(i+is-1,j+js-1) = GAMMA_LT
         ENDIF
         IF (present(id_GAMMA_HW) .AND. id_GAMMA_HW > 0) THEN
            diag_gamma_hw(i+is-1,j+js-1) = GAMMA_HW
         ENDIF
         IF (present(id_GAMMA_SM) .AND. id_GAMMA_SM > 0) THEN
            diag_gamma_sm(i+is-1,j+js-1) = GAMMA_SM
         ENDIF
         IF (present(id_GAMMA_AQ) .AND. id_GAMMA_AQ > 0) THEN
            diag_gamma_aq(i+is-1,j+js-1) = GAMMA_AQ
         ENDIF
         IF (present(id_GAMMA_CO2) .AND. id_GAMMA_CO2 > 0) THEN
            diag_gamma_co2(i+is-1,j+js-1) = GAMMA_CO2
         ENDIF
      ENDIF !land
      ENDDO !lon
      ENDDO !lat

! Apply the canopy loss factor and convert from ug/m2/hr to molecules/m2/s
      EMIS(:,:) = EMIS(:,:) * RHO_CANOPY  *  (1.67e10 / MW_sp)

!---------------------------------------------------------------------------------------
!              ... Accumulate diagnostics
!---------------------------------------------------------------------------------------

! Gamma temperature
      IF (present(id_GAMMA_TEMP) .AND. id_GAMMA_TEMP > 0) THEN
         used = send_data( id_GAMMA_TEMP, diag_gamma_temp(is:ie,js:je),  &
                           Time_next, is_in=is, js_in=js )
      ENDIF
! Gamma light
      IF (present(id_GAMMA_PAR) .AND. id_GAMMA_PAR > 0) THEN
         used = send_data( id_GAMMA_PAR, diag_gamma_par(is:ie,js:je), &
                           Time_next, is_in=is, js_in=js )
      ENDIF
! Gamma LAI
      IF (present(id_GAMMA_LAI) .AND. id_GAMMA_LAI > 0) THEN
         used = send_data( id_GAMMA_LAI, diag_gamma_lai_megan3(is:ie,js:je), &
                           Time_next, is_in=is, js_in=js )
      ENDIF
! Gamma bi-directional LAI
      IF (present(id_GAMMA_BDLAI) .AND. id_GAMMA_BDLAI > 0) THEN
         used = send_data( id_GAMMA_BDLAI, diag_gamma_bdlai_megan3(is:ie,js:je), &
                           Time_next, is_in=is, js_in=js )
      ENDIF
! Gamma Age
      IF (present(id_GAMMA_AGE) .AND. id_GAMMA_AGE > 0) THEN
         used = send_data( id_GAMMA_AGE, diag_gamma_age_megan3(is:ie,js:je), &
                           Time_next, is_in=is, js_in=js )
      ENDIF
! Gamma high temperature
      IF (present(id_GAMMA_HT) .AND. id_GAMMA_HT > 0) THEN
         used = send_data( id_GAMMA_HT, diag_gamma_ht(is:ie,js:je), &
                           Time_next, is_in=is, js_in=js)
      ENDIF
! Gamma low temperature
      IF (present(id_GAMMA_LT) .AND. id_GAMMA_LT > 0) THEN
         used = send_data( id_GAMMA_LT, diag_gamma_lt(is:ie,js:je), &
                           Time_next, is_in=is, js_in=js)
      ENDIF
! Gamma high winds
      IF (present(id_GAMMA_HW) .AND. id_GAMMA_HW > 0) THEN
         used = send_data( id_GAMMA_HW, diag_gamma_hw(is:ie,js:je), &
                           Time_next, is_in=is, js_in=js)
      ENDIF
! Gamma soil
      IF (present(id_GAMMA_SM) .AND. id_GAMMA_SM > 0) THEN
         used = send_data( id_GAMMA_SM, diag_gamma_sm(is:ie,js:je), &
                           Time_next, is_in=is, js_in=js)
      ENDIF
! Gamma air quality
      IF (present(id_GAMMA_AQ) .AND. id_GAMMA_AQ > 0) THEN
         used = send_data(id_GAMMA_AQ, diag_gamma_aq(is:ie,js:je), &
                           Time_next, is_in=is, js_in=js)
      ENDIF
! Gamma co2
      IF (present(id_GAMMA_CO2) .AND. id_GAMMA_CO2 > 0) THEN
         used = send_data(id_GAMMA_CO2, diag_gamma_co2(is:ie,js:je), &
                          Time_next, is_in=is, js_in=js)
      ENDIF

end subroutine calc_xactive_bvoc_megan3
!</SUBROUTINE


!...................................................................................
! ..................................................................................
!...................................................................................




!###################################################################################
!<FUNCTION NAME="fGAMMA_TLD_AM3">
!    ! Light dependendent gamma temperature calculation (AM3)
!
function fGAMMA_TLD_AM3(T1, T_wrk, C_t1, C_eo)

   implicit none

   real, intent(in)   :: T1, T_wrk, C_t1, C_eo
   real               :: T_opt, X, E_opt
   real, parameter    :: C_t2 = 200.

   real               :: fGAMMA_TLD_AM3


      T_opt = 313. + 0.6 * (T_wrk - T_s)
      X     = ((1. / T_opt) - (1. / T1)) / 0.00831
      E_opt = C_eo * exp(0.08 * (T_wrk - T_s))

      fGAMMA_TLD_AM3 = E_opt * C_t2 * exp(  C_t1 * X) /     &
                      (C_t2 - C_t1 * (1. - exp(C_t2 * X)))

end function fGAMMA_TLD_AM3
!-------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------




!###################################################################################
!        FUNCTION GAMMA_TLD
!                       ...  GAMMA TEMPERATURE CALCULATION (AM4)
!------------------------------------------------------------------------------------
function fGAMMA_TLD_AM4(T1, T24, C_t1, C_eo)

   implicit none

   real, intent(in)    :: T1, T24, C_t1, C_eo
   real                :: T_opt, X, E_opt, T240
   real, parameter     :: C_t2 = 230.

   real                :: fGAMMA_TLD_AM4

   T240 = T24       ! Set 10 day average = 24 hr average (Guenther et al., 2018)

   IF ( T1 < 260. ) THEN
      fGAMMA_TLD_AM4 = 0.
   ELSE
      T_opt = 312.5 + 0.6 * (T240 - T_s)
      X     = ((1. / T_opt) - (1. / T1)) / 0.00831
      E_opt = C_eo * exp(0.05 * (T24 - T_s)) * &
                exp(0.05 * (T240 - T_s))

      fGAMMA_TLD_AM4 = E_opt * C_t2 * exp(  C_t1 * X) / &
                       (C_t2 - C_t1 * (1. - exp(C_t2 * X)))
   ENDIF

end function fGAMMA_TLD_AM4
!--------------------------------------------------------------------
!--------------------------------------------------------------------




!###################################################################################
!        FUNCTION GAMMA_PAR
!                       ...  GAMMA LIGHT CALCULATION (AM3)

function fGAMMA_PAR_AM3(PPFD, coszen, P_wrk, calday)

   implicit none

   real, intent(in)   :: PPFD, coszen, P_wrk, calday
   real               :: P_toa, PHI
   real               :: calc_GAMMA_LHT

   real               :: fGAMMA_PAR_AM3

   P_toa              = 3000. + 99. * cos( twopi * (calday - 10.) / 365. ) !G06, Eq. 13
   PHI                = MIN(PPFD / (coszen * P_toa), 1.)                   !G06, Eq. 12

   fGAMMA_PAR_AM3     = coszen * (2.46 * (1. + 0.0005 * (P_wrk - 400.))      &
                        *  PHI - 0.9*PHI*PHI)

end function fGAMMA_PAR_AM3
!------------------------------------------------------------------------
!------------------------------------------------------------------------


!###################################################################################
!        FUNCTION GAMMA_PAR
!                       ...  GAMMA LIGHT CALCULATION (AM4)
!
function fGAMMA_PAR_AM4(PPFD1, PPFD24)

   implicit none

   real, intent(in)   :: PPFD1, PPFD24
   real, parameter    :: ALPHA = 0.004
   real, parameter    :: C1    = 1.03

   real               :: fGAMMA_PAR_AM4

   IF ( PPFD24 < 0.01 ) THEN
      fGAMMA_PAR_AM4 = 0.
   ELSE
      fGAMMA_PAR_AM4 = (ALPHA * C1 * PPFD1) /       &
                       ((1 + ALPHA**2. * PPFD1**2.)**0.5)
   ENDIF

end function fGAMMA_PAR_AM4
!------------------------------------------------------------------------
!------------------------------------------------------------------------



!##########################################################################
!
!        FUNCTION GAMMA_AGE
!                       ...  Cacluate gamma age

!    Modfied from code written by Arlene M. Fiore and Vaishali A. Naik
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

   function fGAMMA_AGE_MEGAN2(LAIp, LAIc, T_wrk, month, &
                       A_new, A_gro, A_mat, A_old, doage)

   implicit none

   real, intent(in)    :: LAIp(:), LAIc(:)
   real, intent(in)    :: T_wrk
   real, intent(in)    :: A_new, A_gro, A_mat, A_old
   integer, intent(in) :: month
   logical, intent(in) :: doage(:)

   integer, parameter :: t(12) = (/ 31,31,28,31,30,31,30,31,31,30,31,30 /)

   integer :: n
   real    :: x, wrk, GAMMA_AGE
   real    :: F_new, F_gro, F_mat, F_old
   real    :: ti, tm                      ! ti = # days btw budbreak and
                                          ! induction of isoprene emissions
                                          ! tm = # days btw budbreak + initiation
                                          ! of peak isop emission rates
   real :: fGAMMA_AGE_MEGAN2(nPFT)

!-------------------------------------------------------------------------------------
!       ... calculations following equations 17&18 in Guenther et al. [2006]
!           -- getting terms needed for gamma age
!-------------------------------------------------------------------------------------
       IF( T_wrk <= 303. ) THEN
          ti = 5. + 0.7*(300. - T_wrk)                     !eqn 18a (see corrigendum)
       ELSE
          ti = 2.9                                          ! eqn 18b (corrigendum)
       ENDIF
       tm = 2.3*ti                                   ! Eq 19
!++amf/van
      fGAMMA_AGE_MEGAN2(1) = 0.
      DO n = 2, nPFT
         IF (doage(n)) THEN
            IF( LAIc(n) == LAIp(n) ) THEN                  !previous LAI = current LAI - p.392 G06
              F_new = 0.0
              F_gro = 0.1
              F_mat = 0.8
              F_old = 0.1
            ELSE IF( LAIp(n) > LAIc(n) ) THEN
              F_new = 0.0
              F_gro = 0.0
              F_old = (LAIp(n) - LAIc(n)) / LAIp(n)
              F_mat = 1. - F_old
            ELSE IF( LAIp(n) < LAIc(n) ) THEN              !LAIp < LAIc
               F_old = 0.
               x    = LAIp(n)/LAIc(n)                           !--amf/van
               wrk  = 1. - x                                        ! = Fnew
               IF( t(month) <= tm ) THEN
                  F_mat =  x                                         ! Eqn 17c
               ELSE
                  F_mat = x + (((t(month) - tm)/t(month) ) * wrk)    ! Eqn 17d
               ENDIF
               IF( t(month) <= ti ) THEN
                  F_new = wrk                                        ! Eqn 17a
                  F_gro = 1. - (F_new + F_mat)                         ! Eqn 17e
               ELSE
                  F_new = (ti/t(month)) * wrk                        ! Eqn 17b
                  F_gro = 1. - (F_new + F_mat)                         ! Eqn 17e
               ENDIF
            end IF
         GAMMA_AGE   = A_new*F_new + A_gro*F_gro + A_mat*F_mat + A_old*F_old   !!  Eq 16
      ELSE
        GAMMA_AGE   = 1.
      ENDIF
      fGAMMA_AGE_MEGAN2(n) = GAMMA_AGE
      ENDDO
   end function fGAMMA_AGE_MEGAN2




!##########################################################################
!
!        FUNCTION GAMMA_AGE
!                       ...  Cacluate gamma age

!    Modfied from code written by Arlene M. Fiore and Vaishali A. Naik
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

   function fGAMMA_AGE_MEGAN3(LAIp, LAIc, T_wrk, month, &
                       A_new, A_gro, A_mat, A_old)

   implicit none

   real, intent(in)    :: LAIp, LAIc
   real, intent(in)    :: T_wrk
   real, intent(in)    :: A_new, A_gro, A_mat, A_old
   integer, intent(in) :: month

   integer, parameter :: t(12) = (/ 31,31,28,31,30,31,30,31,31,30,31,30 /)

   integer :: n
   real    :: x, wrk, GAMMA_AGE
   real    :: F_new, F_gro, F_mat, F_old
   real    :: ti, tm                      ! ti = # days btw budbreak and
                                          ! induction of isoprene emissions
                                          ! tm = # days btw budbreak + initiation
                                          ! of peak isop emission rates
   real :: fGAMMA_AGE_MEGAN3

!-------------------------------------------------------------------------------------
!       ... calculations following equations 17&18 in Guenther et al. [2006]
!           -- getting terms needed for gamma age
!-------------------------------------------------------------------------------------
   IF ( T_wrk <= 303. ) THEN
      ti = 5. + 0.7*(300. - T_wrk)                     !eqn 18a (see corrigendum)
   ELSE
      ti = 2.9                                          ! eqn 18b (corrigendum)
   ENDIF
   tm = 2.3*ti                                   ! Eq 19
   IF ( LAIc == LAIp ) THEN                  !previous LAI = current LAI - p.392 G06
      F_new = 0.0
      F_gro = 0.1
      F_mat = 0.8
      F_old = 0.1
   ELSE IF ( LAIp > LAIc ) THEN
      F_new = 0.0
      F_gro = 0.0
      F_old = (LAIp - LAIc) / LAIp
      F_mat = 1. - F_old
   ELSE IF( LAIp < LAIc ) THEN              !LAIp < LAIc
      F_old = 0.
      x     = LAIp/LAIc                           !--amf/van
      wrk   = 1. - x                                        ! = Fnew
      IF ( t(month) <= tm ) THEN
         F_mat =  x                                         ! Eqn 17c
      ELSE
         F_mat = x + (((t(month) - tm)/t(month) ) * wrk)    ! Eqn 17d
      ENDIF
      IF ( t(month) <= ti ) THEN
         F_new = wrk                                        ! Eqn 17a
         F_gro = 1. - (F_new + F_mat)                         ! Eqn 17e
      ELSE
         F_new = (ti/t(month)) * wrk                        ! Eqn 17b
         F_gro = 1. - (F_new + F_mat)                         ! Eqn 17e
      ENDIF
   ENDIF

   GAMMA_AGE   = A_new*F_new + A_gro*F_gro + A_mat*F_mat + A_old*F_old   !!  Eq 16

   fGAMMA_AGE_MEGAN3 = GAMMA_AGE

   end function fGAMMA_AGE_MEGAN3


!##########################################################################
!   FUNCTION fGAMMA_BDLAI
!                      ... bi-directional gamma lai calculation
!
!  !!ONLY FOR ETHANOL AND ACETALALDEHYDE!!
!!##########################################################################

   function fGAMMA_BDLAI(LAI)

   implicit none

   real, intent(in)   :: LAI

   real               :: fGAMMA_BDLAI

   IF ( LAI < 2. ) THEN
      fGAMMA_BDLAI = 0.5 * LAI
   ELSE IF ( LAI <= 6. .AND. LAI >= 2. ) THEN
      fGAMMA_BDLAI = 1. - 0.0625 * (LAI - 2.)
   ELSE
      fGAMMA_BDLAI = 0.75
   ENDIF


   end function fGAMMA_BDLAI

!!##########################################################################
!     FUNCTION fGAMMA_SM
!              soil moisture gamma calculation
!!##########################################################################

   function fGAMMA_SM(dummy)

   implicit none

   real, intent(in) :: dummy
   real             :: fGAMMA_SM


   !   wilt = WWLT(SLTYP(I,J))
   !   t1 = wilt + d1
   !   IF ( SOILM(I,J) < wilt ) THEN
   !      GAMSM(I,J) = 0
   !   ELSE IF ( SOILM(I,J) >= wilt .AND. SOILM(I,J) < t1 ) THEN
   !       GAMSM(I,J) = (SOILM(I,J) - wilt)/d1
   !   ELSE
   !       GAMSM(I,J) = 1
   !   END IF

      fGAMMA_SM = 1.0
   end function fGAMMA_SM


!##########################################################################
!   FUNCTION fGAMMA_HW
!                      ... Response to high temperatures
!##########################################################################

   function fGAMMA_HT(TMAX1, T_HT, DT_HT, C_HT)

   implicit none

   real, intent(in) :: TMAX1
   real, intent(in) :: T_HT
   real, intent(in) :: DT_HT
   real, intent(in) :: C_HT
!  Local variables
   real             :: THTK, t1
!  Function Declaration
   real             :: fGAMMA_HT

   THTK = 273.15 + T_HT
   t1   = THTK + DT_HT
   IF ( TMAX1 <= THTK ) THEN
      fGAMMA_HT = 1.0
   ELSE IF ( TMAX1 > THTK .AND. TMAX1 < t1 ) THEN
      fGAMMA_HT = 1.0 + (C_HT - 1.0) * (TMAX1 - THTK) / DT_HT
   ELSE
      fGAMMA_HT = C_HT
   ENDIF

   end function fGAMMA_HT


!##########################################################################
!   FUNCTION fGAMMA_LT
!                      ... Response to low temperatures
!##########################################################################

   function fGAMMA_LT(TMIN1, T_LT, DT_LT, C_LT)

   implicit none

   real, intent(in)   :: TMIN1
   real, intent(in)   :: T_LT
   real, intent(in)   :: DT_LT
   real, intent(in)   :: C_LT
   ! Local variables
   real               :: TLTK, t1
   ! Function declaration
   real               :: fGAMMA_LT

   TLTK = 273.15 + T_LT
   t1   = TLTK - DT_LT
   IF ( TMIN1 >= TLTK ) THEN
      fGAMMA_LT = 1.0
   ELSE IF ( TMIN1 < TLTK .AND. TMIN1 > t1 ) THEN
      fGAMMA_LT = 1.0 + (C_LT - 1.0) * (TLTK - TMIN1) / DT_LT
   ELSE
      fGAMMA_LT = C_LT
   ENDIF

   end function fGAMMA_LT



!##########################################################################
!   FUNCTION fGAMMA_HW
!                      ... Response to high winds/storms
!##########################################################################

   function fGAMMA_HW(WSmax1, T_HW, DT_HW, C_HW)

   implicit none

   real, intent(in)   :: WSmax1
   real, intent(in)   :: T_HW
   real, intent(in)   :: DT_HW
   real, intent(in)   :: C_HW
   ! Local Variables
   real               :: t1
   ! Function declaration
   real               :: fGAMMA_HW

   t1 = T_HW + DT_HW
   IF ( WSmax1 <= T_HW ) THEN
      fGAMMA_HW = 1.0
   ELSE IF ( WSmax1 > T_HW .AND. WSmax1 < t1 ) THEN
      fGAMMA_HW = 1.0 + (C_HW - 1.0) * (WSmax1 - T_HW) / DT_HW
   ELSE
      fGAMMA_HW = C_HW
   ENDIF

   end function fGAMMA_HW
!--------------------------------------------------------------------


!##########################################################################
!   FUNCTION fGAMMA_AQ
!                      ... Response to air quality
!##########################################################################

   function fGAMMA_AQ(AQI1, T_AQ, DT_AQ, C_AQ)

   implicit none

   real, intent(in)   :: AQI1      ! Air quality index, W126 of O3
   real, intent(in)   :: T_AQ
   real, intent(in)   :: DT_AQ
   real, intent(in)   :: C_AQ

   real               :: t1

   real               :: fGAMMA_AQ

   t1 = T_AQ + DT_AQ
   IF ( AQI1 <= T_AQ ) THEN
      fGAMMA_AQ = 1.0
   ELSE IF ( AQI1 > T_AQ .AND. AQI1 < t1 ) THEN
      fGAMMA_AQ = 1.0 + (C_AQ - 1.0) * (AQI1 - T_AQ) / DT_AQ
   ELSE
      fGAMMA_AQ = C_AQ
   ENDIF

   end function fGAMMA_AQ
!--------------------------------------------------------------------



!##########################################################################
!   FUNCTION fGAMMA_CO2
!                       ... Response to CO2
!##########################################################################
   function fGAMMA_CO2(CO2)

   implicit none

   real, intent(in)  :: CO2     !CO2 concentration

   real, parameter   :: ISmax = 1.344
   real, parameter   :: hCO2  = 1.4614
   real, parameter   :: Cstar = 585.

   real              :: Ci

   real              :: fGAMMA_CO2

   Ci = 0.7 * CO2

   fGAMMA_CO2 = ISmax - ((ISmax*Ci**hCO2) / (Cstar**hCO2 + Ci**hCO2))

   end function fGAMMA_CO2
!--------------------------------------------------------------------



!##########################################################################
!<SUBROUTINE NAME="temp_init_AM3">
!
!Read in the monthly average temperature data
!##########################################################################
   subroutine temp_init_AM3(lonb,latb,axes)

   real, intent(in), dimension(:,:)      :: lonb, latb
   integer, intent(in)                   :: axes(4)

   integer                               :: nlon, nlat, i, m

   integer, parameter                    :: nlonin = 720
   integer, parameter                    :: nlatin = 360
   integer, parameter                    :: metlonin = 360
   integer, parameter                    :: metlatin = 180

   character(len=35), parameter  :: tasfile = 'INPUT/tas_monthly_clim_1980-2000.nc'
   character(len=5)              :: tasnames(12) = (/'tas01', 'tas02', 'tas03', &
                                                     'tas04', 'tas05', 'tas06', &
                                                     'tas07', 'tas08', 'tas09', &
                                                     'tas10', 'tas11', 'tas12'/)

   integer                               :: id_tas(12)
   real, dimension(nlonin,nlatin)        :: datain
   real, dimension(metlonin,metlatin,12) :: tas
   real, dimension(metlonin)             :: metlon
   real, dimension(metlatin)             :: metlat
   real, dimension(metlonin+1)           :: metlone
   real, dimension(metlatin+1)           :: metlate
   integer, dimension(12)                :: mos
   logical                               :: used
   real                                  :: dlat, dlon
   type(FmsNetcdfFile_t)                 :: tasfile_obj !< Fms2io fileobj
   type (horiz_interp_type)              :: Interp

   nlon = size(lonb,1) - 1
   nlat = size(latb,2) - 1

   IF (open_file(tasfile_obj,tasfile,"read")) then
      IF (mpp_pe() == mpp_root_pe()) call error_mesg ('temp_init_AM3', &
         'Reading NetCDF formatted input file: tas_monthly_clim_1980-2000.nc',NOTE)

!read in lat & lon from input file, get boundaries and convert to radians
      call read_data (tasfile_obj, 'lon', metlon)
      call read_data (tasfile_obj, 'lat', metlat)

      dlon = 0.5*(metlon(1)-metlon(2))
      dlat = 0.5*(metlat(2)-metlat(1))

      DO i = 1, metlatin
         metlate(i) = metlat(i)-dlat
      ENDDO

      metlate(metlatin+1) = metlat(metlatin)+dlat

      DO i = 1, metlonin
         metlone(i) = metlon(i)-dlon
      ENDDO
      metlone(metlonin+1) = metlon(metlonin)+dlon

      metlone = metlone*DEG_TO_RAD
      metlate = metlate*DEG_TO_RAD

      call horiz_interp_init
      call horiz_interp_new ( Interp, metlone, metlate, lonb, latb )

      call read_data (tasfile_obj, 'time', mos)
      call read_data (tasfile_obj,'tas_clim', tas(:,:,:))

      DO m = 1, 12
         call horiz_interp (Interp, tas(:,:,m), Tmo(:,:,m), verbose=verbose)
  !register diagnostic field
         id_tas(m) = register_static_field( 'tracers', tasnames(m), axes(1:2), &
              tasnames(m), 'unitless')
  !send data to diagnostic
         IF (id_tas(m) > 0) THEN
            used = send_data(id_tas(m),Tmo(:,:,m))
         ENDIF
      ENDDO
      call horiz_interp_del( Interp )
      call close_file(tasfile_obj)
   ELSE
      call error_mesg ('temp_init_AM3',  &
          'tasfile :'//tasfile//' does not exist', FATAL)
   ENDIF

end subroutine temp_init_AM3



!##########################################################################
!<SUBROUTINE NAME="ppfd_init_AM3">
!   <OVERVIEW>
!
! Read in the the monthly average PAR data
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
subroutine ppfd_init_AM3 (lonb, latb, axes)

   real, intent(in), dimension(:,:)        :: lonb, latb
   integer, intent(in)                     :: axes(4)

   integer                                 :: nlon, nlat, i, m
   integer, parameter                      :: nlonin = 720, nlatin = 360
   integer, parameter                      :: metlonin = 360, metlatin= 180
   character(len=37), parameter :: dswfile = 'INPUT/dswrf_monthly_clim_1980-2000.nc'

   character(len=5) :: dswnames(12)        = (/'dsw01','dsw02','dsw03','dsw04', &
                                               'dsw05','dsw06','dsw07','dsw08', &
                                               'dsw09','dsw10','dsw11','dsw12' /)
   integer                                 :: id_dsw(12)
   real, dimension(nlonin, nlatin)         :: datain
   real, dimension(metlonin,metlatin,12)   :: dswrf
   real, dimension(metlonin)               :: metlon
   real, dimension(metlatin)               :: metlat
   real, dimension(metlonin+1)             :: metlone
   real, dimension(metlatin+1)             :: metlate
   integer, dimension(12)                  :: mos
   logical                                 :: used
   real                                    :: dlat, dlon
   real, parameter                         :: const0 = 4.766
   type(FmsNetcdfFile_t)                   :: dswfile_obj !< Fms2io fileobj
   type (horiz_interp_type)                :: Interp

   nlon = size(lonb,1) - 1
   nlat = size(latb,2) - 1

!  --- check existence of input file containing climatological (1980-2000)
!  monthly surface down SW radiation --------
  IF (open_file(dswfile_obj,dswfile,"read")) then

!set up for input grid
     IF (mpp_pe() == mpp_root_pe()) call error_mesg ('ppfd_init_AM3',  &
          'Reading NetCDF formatted input file: dswrf_monthly_clim_1980-2000.nc', NOTE)

!read in lat & lon from input file, get boundaries and convert to radians
     call read_data (dswfile_obj, 'lon', metlon)
     call read_data (dswfile_obj, 'lat', metlat)

     dlon = 0.5*(metlon(1)-metlon(2))
     dlat = 0.5*(metlat(2)-metlat(1))

     DO i = 1, metlatin
        metlate(i) = metlat(i)-dlat
     ENDDO

     metlate(metlatin+1) = metlat(metlatin)+dlat

     DO i = 1, metlonin
        metlone(i) = metlon(i)-dlon
     ENDDO
     metlone(metlonin+1) = metlon(metlonin)+dlon

     metlone = metlone*DEG_TO_RAD
     metlate = metlate*DEG_TO_RAD

     call horiz_interp_init
     call horiz_interp_new ( Interp, metlone, metlate, lonb, latb )

     call read_data (dswfile_obj, 'time', mos)
     call read_data (dswfile_obj,'dswrf_clim', dswrf(:,:,:))

     DO m = 1, 12
        call horiz_interp (Interp, dswrf(:,:,m), Pmo(:,:,m),verbose=verbose)
!!! Convert total shortwave to PPFD !!!!!
        Pmo(:,:,m) = Pmo(:,:,m) * const0 * 0.5
!register diagnostic field
        id_dsw(m) = register_static_field( 'tracers', dswnames(m), axes(1:2), &
             dswnames(m), 'unitless')
!send data to diagnostic
        IF (id_dsw(m) > 0) THEN
           used = send_data(id_dsw(m),Pmo(:,:,m))
        ENDIF
     ENDDO
     call horiz_interp_del( Interp )
     call close_file(dswfile_obj)
  ELSE
     call error_mesg ('ppfd_init_AM3',  &
          'dswfile does not exist', FATAL)
  ENDIF

end subroutine ppfd_init_AM3
!</SUBROUTINE>



!#############################################################################################
!
!<SUBROUTINE NAME="pft_init_AM3">
!   <OVERVIEW>
!     Reads in and stores the monthly average LAI data
!   </OVERVIEW>
!   <DESCRIPTION
!     Reads in the file "mksrf_pft
!   </DESCRIPTION>
!   <TEMPLATE>
!
!   </TEMPLATE>
!
subroutine pft_init_AM3( lonb, latb, axes )

   real, intent(in), dimension(:,:)   :: lonb, latb
   integer, intent(in)                :: axes(4)

   integer                            :: nlon, nlat, i
   integer, parameter                 :: nlonin = 720, nlatin = 360
   integer, dimension(nPFT)           :: pft
   real, dimension(nlonin)            :: inlon, lonpft
   real, dimension(nlatin)            :: inlat, latpft
   real, dimension(nlonin+1)          :: lonpfte
   real, dimension(nlatin+1)          :: latpfte
   real                               :: edgen, edgee, edges, edgew, dlat, dlon
   logical                            :: used

   character(len=5) :: pftnames(nPFT) = (/'pft01','pft02','pft03','pft04',  &
                                          'pft05','pft06','pft07','pft08',  &
                                          'pft09','pft10','pft11','pft12',  &
                                          'pft13','pft14','pft15','pft16',  &
                                          'pft17'/)

   integer                                  :: id_pft(nPFT)
   real, dimension(nlonin,nlatin,nPFT)      :: datapft
   type(FmsNetcdfFile_t)                    :: file_PFT_obj !< Fms2io fileobj
   type (horiz_interp_type)                 :: Interp

   nlon = size(lonb,1) - 1
   nlat = size(latb,1) - 1

   IF ( do_AM3_ISOP .AND. file_PFT/='INPUT/mksrf_pft.060929.nc' ) THEN
      call error_mesg ('pft_init', 'incorrect file to reproduce AM3'//  &
                       'isoprene emissions--correcting',WARNING)
      file_PFT = 'INPUT/mksrf_pft.060929.nc'
   ENDIF

   IF ( open_file(file_PFT_obj, file_PFT, "read") ) THEN
      IF ( mpp_pe() == mpp_root_pe() ) call error_mesg ( 'pft_init_AM3', &
           'Reading NetCDF formatted input file: mksrf_pft.060929.nc', NOTE)
! Read in lat & lon from input file, get boundaries and convert to radians
      call read_data (file_PFT_obj, 'lon', lonpft)
      call read_data (file_PFT_obj, 'lat', latpft)
      call read_data (file_PFT_obj, 'EDGEW', edgew)
      call read_data (file_PFT_obj, 'EDGES', edges)
      call read_data (file_PFT_obj, 'EDGEE', edgee)
      call read_data (file_PFT_obj, 'EDGEN', edgen)

      lonpfte(1) = edgew
      latpfte(1) = edges
      latpfte(nlatin+1) = edgen
      lonpfte(nlonin+1) = edgee

      dlon = 2. * (lonpft(1) - edgew)
      dlat = 2. * (latpft(1) - edges)

      DO i = 2, nlatin
         latpfte(i) = latpfte(i-1) + dlat
      ENDDO
      DO i = 2, nlonin
         lonpfte(i) = lonpfte(i-1) + dlon
      ENDDO

      lonpfte = lonpfte*DEG_TO_RAD
      latpfte = latpfte*DEG_TO_RAD

      call horiz_interp_init
      call horiz_interp_new ( Interp, lonpfte, latpfte, lonb, latb )

      call read_data (file_PFT_obj, 'pft', pft)

! Read pct_pft field
      call read_data (file_PFT_obj, 'PCT_PFT', datapft)

! Loop over pftnames
      DO i = 1, nPFT
! Register diagnostic field
         id_pft(i) = register_static_field( 'tracers', pftnames(i), axes(1:2), &
                     pftnames(i), 'unitless')
         call horiz_interp (Interp, datapft(:,:,i), PCTPFT(:,:,i), verbose=verbose)
! Send data to diagnostic
         IF ( id_pft(i) > 0 ) THEN
            used = send_data(id_pft(i), PCTPFT(:,:,i))
         ENDIF
      ENDDO
! Scale the percentages to a fraction
      PCTPFT(:,:,:) = 0.01 * PCTPFT(:,:,:)
      call horiz_interp_del( Interp )
      call close_file(file_PFT_obj)
   ELSE
      call error_mesg ('lai_pft_init', &
           'PFT file: '//file_PFT//' does not exist.', FATAL )
   ENDIF

end subroutine pft_init_AM3


!#############################################################################################
!
! <SUBROUTINE NAME="lai_init_AM3">
!   <OVERVIEW>
!     Reads in monthly LAI data for each PFT and places into the array 'MLAI'
!   </OVERVIEW>
!   <DESCRIPTION
!
!   </DESCRIPTION>
!   <TEMPLATE>
!
!   </TEMPLATE>

!#############################################################################################!
subroutine lai_init_AM3( lonb,latb, axes )

   real, intent(in), dimension(:,:)          :: lonb, latb
   integer, intent(in)                       :: axes(4)
   integer                                   :: nlon, nlat, i, m
   integer, parameter                        :: nlonin = 720, nlatin = 360
   real, dimension(nlonin+1)                 :: lonlaie
   real, dimension(nlatin+1)                 :: latlaie
   real, dimension(nlonin)                   :: inlon, lonlai
   real, dimension(nlatin)                   :: inlat, latlai
   real                                      :: edgen, edgee, edges, edgew
   real                                      :: dlat, dlon
   integer                                   :: id_lai(nPFT)
   integer, dimension(nMOS)                  :: mos
   logical                                   :: used
   real, dimension(nlonin,nlatin,nPFT,nMOS)  :: datalai
   type(FmsNetcdfFile_t)                     :: file_LAI_obj !< Fms2io fileobj
   character(len=5)  :: lainames(nPFT) =  (/'lai01','lai02','lai03','lai04', &
                                            'lai05','lai06','lai07','lai08', &
                                            'lai09','lai10','lai11','lai12', &
                                            'lai13','lai14','lai15','lai16', &
                                            'lai17'/)
   type (horiz_interp_type)                  :: Interp

   nlon = size(lonb,1) - 1
   nlat = size(latb,1) - 1

   IF ( do_AM3_ISOP .AND. file_LAI/='INPUT/mksrf_lai.060929.nc' ) THEN
      call error_mesg ('lai_init, incorrect file to reproduce AM3', &
                       'isoprene emissions--correcting',WARNING)
      file_LAI = 'INPUT/mksrf_lai.060929.nc'
   ENDIF

!  --- check existence of input file containing monthly lai, for each pft
!  --------
   IF (open_file(file_LAI_obj, file_LAI, "read")) THEN
! Set up for input grid
      IF(mpp_pe() == mpp_root_pe()) call error_mesg ('lai_pft_init',  &
           'Reading NetCDF formatted input file: mksrf_lai.060929.nc', NOTE)
! Read in lat & lon from input file, get boundaries and convert to radians
         call read_data (file_LAI_obj, 'lon', lonlai)
         call read_data (file_LAI_obj, 'lat', latlai)
         call read_data (file_LAI_obj, 'EDGEW', edgew)
         call read_data (file_LAI_obj, 'EDGES', edges)
         call read_data (file_LAI_obj, 'EDGEE', edgee)
         call read_data (file_LAI_obj, 'EDGEN', edgen)
! Get lat/lon edges and spacing
         lonlaie(1) = edgew
         latlaie(1) = edges
         latlaie(nlatin+1) = edgen
         lonlaie(nlonin+1) = edgee
         dlon = 2. * (lonlai(1) - edgew)
         dlat = 2. * (latlai(1) - edges)
         DO i = 2, nlatin
            latlaie(i) = latlaie(i-1) + dlat
         ENDDO
         DO i = 2, nlonin
            lonlaie(i) = lonlaie(i-1) + dlon
         ENDDO
         lonlaie = lonlaie*DEG_TO_RAD
         latlaie = latlaie*DEG_TO_RAD

         call horiz_interp_init
         call horiz_interp_new ( Interp, lonlaie, latlaie, lonb, latb )
! Read in pft and time dimensions from lai file
         call read_data (file_LAI_obj, 'time', mos)
! Loop over pftnames
         DO i = 1, nPFT
            call read_data (file_LAI_obj,lainames(i),datalai(:,:,i,:))
            DO m = 1, nMOS
               call horiz_interp (Interp, datalai(:,:,i,m),        &
                                  MLAI(:,:,i,m),verbose=verbose)
! Store diagnostics for one month only - choose July for now
            IF (m .eq. 7) THEN
! Register diagnostic field
               id_lai(i) = register_static_field( 'tracers', lainames(i), &
                           axes(1:2), lainames(i), 'unitless')
! Send data to diagnostic
               IF (id_lai(i) > 0) THEN
                  used = send_data(id_lai(i),MLAI(:,:,i,m))
               ENDIF
            ENDIF
         ENDDO
      ENDDO
      call horiz_interp_del( Interp )
      call close_file(file_LAI_obj)
   ELSE
      call error_mesg ('lai_init_AM3',  &
           'laifile: '//file_LAI//' does not exist', FATAL)
   ENDIF
end subroutine lai_init_AM3
!</SUBROUTINE>



!#############################################################################################
!
! <SUBROUTINE NAME="lai_init_megan3">
!   <OVERVIEW>
!     Reads in monthly LAI data and places into the array 'MLAI_MEGAN3'
!   </OVERVIEW>
!   <DESCRIPTION
!
!   </DESCRIPTION>
!   <TEMPLATE>
!
!   </TEMPLATE>

!#############################################################################################!

subroutine lai_init_megan3( lonb,latb, axes )

   real, intent(in), dimension(:,:)          :: lonb, latb
   integer, intent(in)                       :: axes(4)
   integer                                   :: m
   integer, parameter                        :: nlonin = 720, nlatin = 360
   real, dimension(nlonin+1)                 :: inlone
   real, dimension(nlatin+1)                 :: inlate
   real, dimension(nlonin)                   :: inlon
   real, dimension(nlatin)                   :: inlat
   real                                      :: dlat, dlon
   integer                                   :: id_lai
   logical                                   :: used
   real, dimension(nlonin,nlatin,nMOS)       :: datalai
   type(FmsNetcdfFile_t)                     :: file_LAI_obj !< Fms2io fileobj
   type (horiz_interp_type)                  :: Interp

   IF ( file_LAI =='INPUT/mksrf_lai.060929.nc' ) THEN
      call error_mesg ('lai_init_megan3, incorrect file for MEGAN3', &
                       'needs to be combined LAI--correcting',WARNING)
      file_LAI = 'INPUT/mksrf_lai.060929.combined_pft.nc'
   ENDIF

   IF (open_file(file_LAI_obj, file_LAI, "read")) THEN
! Set up for input grid
      IF (mpp_pe() == mpp_root_pe()) call error_mesg ('lai_init_megan3',  &
           'Reading NetCDF formatted input file'//file_LAI, NOTE)

      call read_data (file_LAI_obj, 'lon', inlon)
      call read_data (file_LAI_obj, 'lat', inlat)
      inlon = inlon*DEG_TO_RAD
      inlat = inlat*DEG_TO_RAD
      dlat  = inlat(2)-inlat(1)
      dlon  = inlon(2)-inlon(1)
      inlone(1:nlonin) = inlon-(dlon/2.)
      inlone(nlonin+1) = inlon(nlonin)+(dlon/2.)
      inlate(1:nlatin) = inlat-(dlat/2.)
      inlate(nlatin+1) = inlat(nlatin)+(dlat/2.)

      call horiz_interp_init
      call horiz_interp_new ( Interp, inlone, inlate, lonb, latb )
      call read_data (file_LAI_obj,'LAI',datalai)
      DO m = 1, nMOS
         call horiz_interp (Interp, datalai(:,:,m), MLAI_MEGAN3(:,:,m),verbose=verbose)
! Store diagnostics for one month only - choose July for now
         IF (m .eq. 7) THEN
! Register diagnostic field
            id_lai = register_static_field( 'tracers', 'LAIv', &
                        axes(1:2), 'LAIv', 'unitless')
! Send data to diagnostic
            IF (id_lai > 0) THEN
                  used = send_data(id_lai,MLAI_MEGAN3(:,:,m))
            ENDIF
         ENDIF
      ENDDO
      call horiz_interp_del( Interp )
      call close_file(file_LAI_obj)
   ELSE
      call error_mesg ('lai_init_megan3',  &
           'laifile: '//file_LAI//' does not exist', FATAL)
   ENDIF
end subroutine lai_init_megan3
!</SUBROUTINE>


!#############################################################################################
!
! <SUBROUTINE NAME="fcover_init_megan3">
!   <OVERVIEW>
!     Reads in monthly LAI data and places into the array 'MLAI_MEGAN3'
!   </OVERVIEW>
!   <DESCRIPTION
!
!   </DESCRIPTION>
!   <TEMPLATE>
!
!   </TEMPLATE>

!#############################################################################################!

subroutine fcover_init_megan3( lonb,latb, axes )

   real, intent(in), dimension(:,:)          :: lonb, latb
   integer, intent(in)                       :: axes(4)
   integer                                   :: m
   integer, parameter                        :: nlonin = 720, nlatin = 360
   real, dimension(nlonin+1)                 :: inlone
   real, dimension(nlatin+1)                 :: inlate
   real, dimension(nlonin)                   :: inlon
   real, dimension(nlatin)                   :: inlat
   real                                      :: dlat, dlon
   integer                                   :: id_fcover
   logical                                   :: used
   real, dimension(nlonin,nlatin,nMOS)       :: datain
   type(FmsNetcdfFile_t)                     :: file_FCOVER_obj !< Fms2io fileobj
   type (horiz_interp_type)                  :: Interp

   IF (open_file(file_FCOVER_obj, file_FCOVER, "read")) THEN
! Set up for input grid
      IF (mpp_pe() == mpp_root_pe()) call error_mesg ('fcover_init_megan3',  &
           'Reading NetCDF formatted input file'//file_FCOVER, NOTE)

      call read_data (file_FCOVER_obj, 'lon', inlon)
      call read_data (file_FCOVER_obj, 'lat', inlat)
      inlon = inlon*DEG_TO_RAD
      inlat = inlat*DEG_TO_RAD
      dlat  = inlat(2)-inlat(1)
      dlon  = inlon(2)-inlon(1)
      inlone(1:nlonin) = inlon-(dlon/2.)
      inlone(nlonin+1) = inlon(nlonin)+(dlon/2.)
      inlate(1:nlatin) = inlat-(dlat/2.)
      inlate(nlatin+1) = inlat(nlatin)+(dlat/2.)

      call horiz_interp_init
      call horiz_interp_new ( Interp, inlone, inlate, lonb, latb )
      call read_data (file_FCOVER_obj,'FCOVER',datain(:,:,:))
      DO m = 1, nMOS
         call horiz_interp (Interp, datain(:,:,m), FCOVER(:,:,m),verbose=verbose)
! Store diagnostics for one month only - choose July for now
         IF (m .eq. 7) THEN
! Register diagnostic field
            id_fcover = register_static_field( 'tracers', 'FCOVER', &
                        axes(1:2), 'FCOVER', 'unitless')
! Send data to diagnostic
            IF (id_fcover > 0) THEN
                  used = send_data(id_fcover,FCOVER(:,:,m))
            ENDIF
         ENDIF
      ENDDO
      call horiz_interp_del( Interp )
      call close_file(file_FCOVER_obj)
   ELSE
      call error_mesg ('fcover_init_megan3',  &
           'fcover file: '//file_FCOVER//' does not exist', FATAL)
   ENDIF
end subroutine fcover_init_megan3
!</SUBROUTINE>


!#####################################################################
! <SUBROUTINE NAME="xactive_bvoc_register_restart">
!  <OVERVIEW>
!    xactive_bvoc_register_restart registers restart fields
!  </OVERVIEW>
subroutine xactive_bvoc_register_restart_scalars(Xbvoc_restart)
   type(FmsNetcdfFile_t), intent(inout)       ::  Xbvoc_restart !< Fms2io fileobj
   character(len=8), dimension(1)             :: dim_names !< Array of dimension names

   dim_names =  (/"Time"/)
   call register_axis(Xbvoc_restart, "Time", unlimited)

   call register_restart_field(Xbvoc_restart, 'version', vers)

end subroutine xactive_bvoc_register_restart_scalars

!< xactive_bvoc_register_restart_domains: register netcdf restart variables
subroutine xactive_bvoc_register_restart_domains(Til_restart)
   type(FmsNetcdfDomainFile_t), intent(inout)       ::  Til_restart !< Fms2io domain decomposed fileobj
   character(len=2)  :: mon_string
   character(len=8), dimension(3)             :: dim_names !< Array of dimension names
   integer                                    :: ihour

   if (Ldebug .and. mpp_pe()==mpp_root_pe()) &
      write(*,*) 'xactive_bvoc_register_restart: ', &
      'T24_STORE,P24_STORE,WS_STORE,O3_STORE=', ALLOCATED(T24_STORE), &
      ALLOCATED(P24_STORE), ALLOCATED(WS_STORE), ALLOCATED(O3_STORE)

  dim_names =  (/"xaxis_1", "yaxis_1", "Time   "/)

  call register_axis(Til_restart, "xaxis_1", "x")
  call register_axis(Til_restart, "yaxis_1", "y")
  call register_axis(Til_restart, "Time", unlimited)

  !< Register the domain decomposed dimensions as variables so that the combiner can work
  !! correctly
  call register_field(Til_restart, dim_names(1), "double", (/dim_names(1)/))
  call register_field(Til_restart, dim_names(2), "double", (/dim_names(2)/))

  do ihour = 1,24
     write(mon_string,'(i2.2)') ihour
     if (Ldebug .and. mpp_pe()==mpp_root_pe()) &
        write(*,*) 'xactive_bvoc_register_restart: register field T24_STORE_'//mon_string
     if (ALLOCATED(T24_STORE)) call &
        register_restart_field(Til_restart, 'T24_STORE_'//mon_string, T24_STORE(:,:,ihour), dim_names, is_optional = .true.)
     if (ALLOCATED(P24_STORE)) call &
        register_restart_field(Til_restart, 'P24_STORE_'//mon_string, P24_STORE(:,:,ihour), dim_names, is_optional = .true.)
     if (ALLOCATED(WS_STORE)) call &
        register_restart_field(Til_restart, 'WS_STORE_'//mon_string,  WS_STORE(:,:,ihour), dim_names, is_optional = .true.)
     if (ALLOCATED(O3_STORE)) call &
        register_restart_field(Til_restart, 'O3_STORE_'//mon_string,  O3_STORE(:,:,ihour), dim_names, is_optional = .true.)
   end do
end subroutine xactive_bvoc_register_restart_domains
! </SUBROUTINE>




!##########################################################################
!
! <SUBROUTINE NAME="xactive_bvoc_end">
!   <OVERVIEW>
!     Ends the xactive bvoc module by deallocating any used arrays
!   </OVERVIEW>
!
subroutine xactive_bvoc_end
   type(FmsNetcdfFile_t)       ::  Xbvoc_restart !< Fms2io fileobj
   type(FmsNetcdfDomainFile_t) ::  Til_restart !< Fms2io domain decomposed fileobj
   logical :: tile_file_open !< Flag indicated whether the restart file was opened sucessfully
   integer, allocatable, dimension(:) :: pes !< Array of pes in the current pelist

   if (Ldebug .and. mpp_pe()==mpp_root_pe()) write(*,*) 'xactive_bvoc_end: calling save_restart'

   !< Get the current pelist
   allocate(pes(mpp_npes()))
   call mpp_get_current_pelist(pes)

   !< Open the scalar file with the current pelist, so that only the root pe opens and writes the file
   if (open_file(Xbvoc_restart,"RESTART/xactive_bvoc.res.nc","overwrite", is_restart=.true., pelist=pes)) then
      call xactive_bvoc_register_restart_scalars(Xbvoc_restart)
      call write_restart(Xbvoc_restart)
      call close_file(Xbvoc_restart)
   endif
   deallocate(pes)

   if (Ldebug .and. mpp_pe()==mpp_root_pe()) write(*,*) 'xactive_bvoc_end: calling save_restart_Til'
   if (mpp_get_ntile_count(xactive_domain) == 1) then
      tile_file_open = open_file(Til_restart,"RESTART/xactive_bvoc.res.nc","append", xactive_domain, is_restart=.true.)
   else
      tile_file_open = open_file(Til_restart,"RESTART/xactive_bvoc.res.nc","overwrite", xactive_domain, is_restart=.true.)
   endif

   if (tile_file_open) then
      call xactive_bvoc_register_restart_domains(Til_restart)
      call write_restart(Til_restart)
      call add_domain_dimension_data(Til_restart)
      call close_file(Til_restart)
   endif

   if (Ldebug .and. mpp_pe()==mpp_root_pe()) write(*,*) 'xactive_bvoc_end: back from save_restart_Til'

   IF (mpp_pe() == mpp_root_pe()) THEN
      write(*,*) 'xactive_bvoc_end: Deallocating xactive arrays'
   ENDIF

   IF ( ALLOCATED(ECISOP_AM3) )              DEALLOCATE(ECISOP_AM3)
   IF ( ALLOCATED(ECBVOC) )                  DEALLOCATE(ECBVOC)
   IF ( ALLOCATED(ECBVOC_MEGAN3) )           DEALLOCATE(ECBVOC_MEGAN3)
   IF ( ALLOCATED(ECTERP) )                  DEALLOCATE(ECTERP)
   IF ( ALLOCATED(ECTERP_MEGAN3) )           DEALLOCATE(ECTERP_MEGAN3)
   IF ( ALLOCATED(MEGAN_PARAM) )             DEALLOCATE(MEGAN_PARAM)
   IF ( ALLOCATED(TERP_PARAM) )              DEALLOCATE(TERP_PARAM)
   IF ( ALLOCATED(LDFg) )                    DEALLOCATE(LDFg)
   IF ( ALLOCATED(LDFg_TERP) )               DEALLOCATE(LDFg_TERP)
   IF ( ALLOCATED(MLAI) )                    DEALLOCATE(MLAI)
   IF ( ALLOCATED(MLAI_MEGAN3) )             DEALLOCATE(MLAI_MEGAN3)
   IF ( ALLOCATED(PCTPFT) )                  DEALLOCATE(PCTPFT)
   IF ( ALLOCATED(FCOVER) )                  DEALLOCATE(FCOVER)
   IF ( ALLOCATED(T24_STORE) )               DEALLOCATE(T24_STORE)
   IF ( ALLOCATED(P24_STORE) )               DEALLOCATE(P24_STORE)
   IF ( ALLOCATED(WS_STORE) )                DEALLOCATE(WS_STORE)
   IF ( ALLOCATED(O3_STORE) )                DEALLOCATE(O3_STORE)
   IF ( ALLOCATED(Tmo) )                     DEALLOCATE(Tmo)
   IF ( ALLOCATED(Pmo) )                     DEALLOCATE(Pmo)
   IF ( ALLOCATED(CO2_STORE) )               DEALLOCATE(CO2_STORE)
   IF ( ALLOCATED(SOILM) )                   DEALLOCATE(SOILM)
   IF ( ALLOCATED(WILT) )                    DEALLOCATE(WILT)
   IF ( ALLOCATED(diag_gamma_temp) )         DEALLOCATE(diag_gamma_temp)
   IF ( ALLOCATED(diag_gamma_par) )          DEALLOCATE(diag_gamma_par)
   IF ( ALLOCATED(diag_gamma_ht) )           DEALLOCATE(diag_gamma_ht)
   IF ( ALLOCATED(diag_gamma_lt) )           DEALLOCATE(diag_gamma_lt)
   IF ( ALLOCATED(diag_gamma_hw) )           DEALLOCATE(diag_gamma_hw)
   IF ( ALLOCATED(diag_gamma_aq) )           DEALLOCATE(diag_gamma_aq)
   IF ( ALLOCATED(diag_gamma_co2) )          DEALLOCATE(diag_gamma_co2)
   IF ( ALLOCATED(diag_gamma_sm) )           DEALLOCATE(diag_gamma_sm)
   IF ( ALLOCATED(diag_gamma_age) )          DEALLOCATE(diag_gamma_age)
   IF ( ALLOCATED(diag_gamma_lai) )          DEALLOCATE(diag_gamma_lai)
   IF ( ALLOCATED(diag_gamma_age_megan3) )   DEALLOCATE(diag_gamma_age_megan3)
   IF ( ALLOCATED(diag_gamma_lai_megan3) )   DEALLOCATE(diag_gamma_lai_megan3)
   IF ( ALLOCATED(diag_gamma_bdlai_megan3) ) DEALLOCATE(diag_gamma_bdlai_megan3)

   IF (mpp_pe() == mpp_root_pe()) THEN
       write(*,*) 'Finished deallocating xactive arrays'
   ENDIF

   module_is_initialized = .FALSE.


end subroutine xactive_bvoc_end
!</SUBROUTINE>

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


!############################################################################
end module xactive_bvoc_mod
