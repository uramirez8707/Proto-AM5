                 module cloud_spec_mod
! <CONTACT EMAIL="Fei.Liu@noaa.gov">
!  fil
! </CONTACT>
! <REVIEWER EMAIL="Dan.Schwarzkopf@noaa.gov">
!  ds
! </REVIEWER>
! <HISTORY SRC="http://www.gfdl.noaa.gov/fms-cgi-bin/cvsweb.cgi/FMS/"/>
! <OVERVIEW>
!    cloud_spec_mod defines the variables that are used in a partic-
!    ular cloud parameterization to specify the cloud location, cloud
!    type and cloud magnitude for the active cloud parameterization(s).
! </OVERVIEW>
! <DESCRIPTION>
!    if microphysically-based radiative properties are desired, then
!    cloud_spec_mod also provides the microphysical parameters used in
!    determining the radiative properties, either from the cloud scheme
!    itself if they are present.
! </DESCRIPTION>

!   shared modules:

use time_manager_mod,         only: time_type, time_manager_init, &
                                    set_time, operator (+)
use mpp_mod,                  only: input_nml_file
use fms_mod,                  only: mpp_pe, &
                                    mpp_root_pe, stdlog,  fms_init, &
                                    write_version_number, &
                                    check_nml_error, error_mesg,   &
                                    FATAL, NOTE, stdout
use fms2_io_mod,              only: file_exists
use tracer_manager_mod,       only:         &
!                                   tracer_manager_init,  &
                                    get_tracer_index, NO_TRACER
use field_manager_mod,        only:       &
                                    field_manager_init, &
                                    MODEL_ATMOS
use data_override_mod,        only: data_override
use random_number_streams_mod, only: random_number_streams_init, &
                                     get_random_number_streams, &
                                     random_number_streams_end
use random_numbers_mod,    only:  randomNumberStream,   &
                                  getRandomNumbers
use constants_mod,         only : radian, RDGAS

use aerosol_types_mod,        only: aerosol_type

! atmos param modules:

use physics_radiation_exch_mod, only: clouds_from_moist_block_type, &
                                      exchange_control_type

! cloud radiation modules

use cloudrad_types_mod,       only: cld_specification_type, &
                                    microphysics_type,  &
                                    cloudrad_control_type

use strat_clouds_W_mod,       only: strat_clouds_W_init,   &
                                    strat_clouds_amt, strat_clouds_W_end

use uw_clouds_W_mod,          only: uw_clouds_W_init, &
                                    uw_clouds_amt, &
                                    uw_clouds_W_end

!BW use rh_based_clouds_mod,      only: rh_based_clouds_init,  &
!BW                                     rh_clouds_amt, &
!BW                                     rh_based_clouds_end

!--------------------------------------------------------------------

implicit none
private

!--------------------------------------------------------------------
!    cloud_spec_mod defines the variables that are used in a partic-
!    ular cloud parameterization to specify the cloud location, cloud
!    type and cloud magnitude for the active cloud parameterization(s).
!    if microphysically-based radiative properties are desired, then
!    cloud_spec_mod also provides the microphysical parameters used in
!    determining the radiative properties, either from the cloud scheme
!    itself if they are present.
!--------------------------------------------------------------------


!---------------------------------------------------------------------
!----------- version number for this module --------------------------

character(len=128)  :: version =  '$Id$'
character(len=128)  :: tagname =  '$Name$'


!---------------------------------------------------------------------
!-------  interfaces --------

public          &
         cloud_spec_init, cloud_spec, cloud_spec_end

private    &

!  called from cloud_spec:
         combine_cloud_properties


!---------------------------------------------------------------------
!-------- namelist  ---------

character(len=16)  ::      &
              cloud_type_form = 'stratuw' ! cloud parameterization being
                                          ! used; either 'strat', 'rh',
                                          ! 'stratuw', or 'none'
real :: wtr_cld_reff=10.                ! assumed cloud drop efective
                                        ! radius [ microns ]
real :: ice_cld_reff=50.                ! assumed ice cloud effective
                                        ! size [ microns ]
real :: rain_reff=250.                  ! assumed rain drop effective
                                        ! radius [ microns ]
character(len=16) :: overlap_type = 'random'
                                        ! cloud overlap assumption;
                                        ! allowable values are 'random'
                                        ! or 'max-random'
logical :: doing_data_override=.false.
logical :: do_fu2007 = .false.
logical :: do_rain   = .false. !sjl
logical :: do_snow   = .false. !miz
logical :: do_graupel  = .false. !sjl

logical   :: do_stochastic_clouds = .true.

logical :: use_cloud_tracers_in_radiation = .false.
                               ! if true, use lsc cloud tracer fields
                               ! in radiation (these transported on
                               ! current step, will have non-realizable
                               ! total cloud areas at some points); if
                               ! false, then use balanced (realizable)
                               ! fields saved at end of last step
                               ! only an issue when both lsc and conv
                               ! clouds are active (AM3)


namelist /cloud_spec_nml / cloud_type_form, wtr_cld_reff,   &
                           ice_cld_reff, rain_reff, overlap_type, &
                           doing_data_override, do_fu2007,    &
                           do_rain, do_snow, do_graupel, &
                           do_stochastic_clouds, &
                           use_cloud_tracers_in_radiation

!----------------------------------------------------------------------
!----  public data -------


!----------------------------------------------------------------------
!----  private data -------

!--------------------------------------------------------------------
!    assumed water paths.
!--------------------------------------------------------------------
real   ::  lwpath_hi  = 6.313929   ! assumed water path for high clouds
                                   ! [ grams / m**2 ]
real   ::  lwpath_mid = 18.94179   ! assumed water path for middle
                                   ! clouds [ grams / m**2 ]
real   ::  lwpath_low = 75.76714   ! assumed water path for low clouds
                                   ! [ grams / m**2 ]

!---------------------------------------------------------------------
!    logical  flags.

logical :: module_is_initialized = .false.   ! module initialized ?

!---------------------------------------------------------------------
!    time-step related constants.

integer :: num_pts       !  number of grid columns processed so far that
                         !  have cloud data present (used to identify
                         !  module coldstart condition)
integer :: tot_pts       !  total number of grid columns in the
                         !  processor's domain

!---------------------------------------------------------------------
!     indices for cloud tracers

integer :: nql           ! tracer index for liquid water
integer :: nqi           ! tracer index for ice water
integer :: nqa           ! tracer index for cloud area
integer :: nqn           ! tracer index for cloud droplet number
integer :: nqni          ! tracer index for ice crystal number
integer :: nqr, nqs, nqg ! tracer index for rainwat, snowwat and graupel


!----------------------------------------------------------------------
!     miscellaneous variables:

!BW integer :: num_slingo_bands  ! number of radiative bands over which
!BW                              ! cloud optical depth is calculated in the
!BW                              ! gordon diag_cloud parameterization

integer :: id, jd, kmax

type(time_type) :: Radiation_time_step  ! saved for data override

logical :: doing_prog_clouds

!----------------------------------------------------------------------
!----------------------------------------------------------------------



                         contains



!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!                     PUBLIC SUBROUTINES
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!#####################################################################
! <SUBROUTINE NAME="cloud_spec_init">
!  <OVERVIEW>
!   Contructor of cloud_spec_package module
!  </OVERVIEW>
!  <DESCRIPTION>
!   Contructor of cloud_spec_package module
!  </DESCRIPTION>
!  <TEMPLATE>
!   call cloud_spec_init ( pref, lonb, latb, axes, Time)
!  </TEMPLATE>
!  <IN NAME="pref" TYPE="real">
!   reference pressure levels containing two reference pressure profiles
!                 for use in defining transmission functions [ Pa ]
!  </IN>
!  <IN NAME="lonb" TYPE="real">
!   the longitude array of the model grid box corners
!  </IN>
!  <IN NAME="latb" TYPE="real">
!   the latitude array of the model grid box corners
!  </IN>
!  <IN NAME="axes" TYPE="real">
!   diagnostic variable axes for netcdf files
!  </IN>
!  <IN NAME="Time" TYPE="time_type">
!   current time [ time_type(days, seconds) ]
!  </IN>
! </SUBROUTINE>
!
subroutine cloud_spec_init (Exch_ctrl, pref, lonb, latb, axes, Time,   &
                            rad_time_step, Cldrad_control)

!---------------------------------------------------------------------
!    cloud_spec_init is the constructor for cloud_spec_mod.
!---------------------------------------------------------------------

type(exchange_control_type), intent(inout) :: Exch_ctrl
real, dimension(:,:),        intent(in)    ::  pref
real, dimension(:,:),        intent(in)    ::  lonb, latb
integer, dimension(4),       intent(in)    ::  axes
type(time_type),             intent(in)    ::  Time
integer,                     intent(in)    ::  rad_time_step
type(cloudrad_control_type), intent(inout) ::  Cldrad_control

!-------------------------------------------------------------------
!    intent(in) variables:
!
!       pref      array containing two reference pressure profiles
!                 for use in defining transmission functions [ Pa ]
!       lonb      array of model longitudes at cell corners [ radians ]
!       latb      array of model latitudes at cell corners [radians]
!       axes      diagnostic variable axes
!       Time      current time [time_type(days, seconds)]
!
!----------------------------------------------------------------------

!----------------------------------------------------------------------
!   local variables:

      integer   ::   ierr, io, logunit
      integer   ::   ndum, i, j, ii, jj


!--------------------------------------------------------------------
!   local variables:
!
!      ierr     error code
!      io       error status returned from io operation
!      ndum     dummy argument needed for call to field_manager_init
!
!--------------------------------------------------------------------

!---------------------------------------------------------------------
!    if routine has already been executed, exit.
!---------------------------------------------------------------------
      if (module_is_initialized) return

!---------------------------------------------------------------------
!    verify that modules used by this module that are not called later
!    have already been initialized.
!---------------------------------------------------------------------
      call fms_init
      call time_manager_init
      call field_manager_init (ndum)
!  not yet compliant:
!     call tracer_manager_init  ! not public

!---------------------------------------------------------------------
!    read namelist.
      read (input_nml_file, nml=cloud_spec_nml, iostat=io)
      ierr = check_nml_error(io,"cloud_spec_nml")

!----------------------------------------------------------------------
!    write version number and namelist to logfile.
!---------------------------------------------------------------------
      call write_version_number (version, tagname)
      logunit = stdlog()
      if (mpp_pe() == mpp_root_pe() ) &
           write (logunit, nml=cloud_spec_nml)

      id = size(lonb,1) - 1
      jd = size(latb,2) - 1
      kmax = size(pref,1) - 1

!-----------------------------------------------------------------------
!    define output field.
!-----------------------------------------------------------------------
      Exch_ctrl%cloud_type_form = cloud_type_form

!--------------------------------------------------------------------
!    verify a valid type of cloud overlap. set logical variables
!    based on the namelist value.
!--------------------------------------------------------------------
      if (trim(overlap_type) == 'random') then
        Cldrad_control%do_random_overlap = .true.
      else if (trim(overlap_type) == 'max-random') then
        Cldrad_control%do_max_random_overlap = .true.
      else
        call error_mesg ('cloud_spec_mod',  &
         ' invalid specification of overlap_type', FATAL)
      endif

      doing_prog_clouds = Exch_ctrl%doing_prog_clouds

!-------------------------------------------------------------------
!    set the variables indicating that the above control variables have
!    been set.
!--------------------------------------------------------------------
!BW   Cldrad_control%do_random_overlap_iz = .true.
!BW   Cldrad_control%do_max_random_overlap_iz = .true.

!--------------------------------------------------------------------
!    save the flags indicating whether stochastic clouds are to be
!    used.
!--------------------------------------------------------------------
      Cldrad_control%do_stochastic_clouds = do_stochastic_clouds
!BW   Cldrad_control%do_stochastic_clouds_iz = .true.

!--------------------------------------------------------------------
!    if stochastic clouds is active, be sure that the
!    cloud_generator module has been initialized.
!--------------------------------------------------------------------
      if (Cldrad_control%do_stochastic_clouds) then
          call random_number_streams_init ( lonb, latb, Cldrad_control )
      endif

!-------------------------------------------------------------------
!    verify that the nml variable cloud_type_form specifies a valid
!    cloud parameterization. set the appropriate logical control
!    variable(s) to .true.. call the constructor modules for the
!    specific cloud scheme(s) requested.
!-------------------------------------------------------------------
      if (trim(cloud_type_form) == 'strat')  then

!-------------------------------------------------------------------
!    cloud fractions, heights are predicted by the model based on klein
!    parameterization.
!-------------------------------------------------------------------
         Cldrad_control%do_strat_clouds = .true.

!-------------------------------------------------------------------
!    cloud fractions, heights are diagnosed based on model relative
!    humidity.
!-------------------------------------------------------------------
!BW   else if (trim(cloud_type_form)  == 'rh')   then
!BW      Cldrad_control%do_rh_clouds = .true.
!BW      call rh_based_clouds_init

!------------------------------------------------------------------
!    cloud fractions, heights are provided by the uw_conv shallow
!    convection scheme
!-------------------------------------------------------------------
      else if (trim(cloud_type_form) == 'uw')  then
         Cldrad_control%do_uw_clouds = .true.

!-------------------------------------------------------------------
!    cloud fractions, heights are provided by the klein large-scale
!    and uw_conv shallow convection cloud parameterizations.
!-------------------------------------------------------------------
      else if (trim(cloud_type_form) == 'stratuw')  then
         Cldrad_control%do_strat_clouds = .true.
         Cldrad_control%do_uw_clouds = .true.

!---------------------------------------------------------------
!    model is run without clouds.
!-------------------------------------------------------------------
      else if (trim(cloud_type_form) == 'none')  then
         Cldrad_control%do_no_clouds = .true.

!-------------------------------------------------------------------
!    failure message if none of the above options was chosen.
!-------------------------------------------------------------------
      else
         call error_mesg ('cloud_spec_mod',  &
              'invalid cloud_type_form specified', FATAL)
      endif  ! (strat)

!--------------------------------------------------------------------
!    initialize cloud schemes that are used
!--------------------------------------------------------------------
      if (Cldrad_control%do_strat_clouds) then
         call strat_clouds_W_init(latb, lonb, Cldrad_control, Exch_ctrl)
      endif
      if (Cldrad_control%do_uw_clouds) then
         call uw_clouds_W_init (Exch_ctrl)
      endif

!--------------------------------------------------------------------
!    define the dimensions of the model subdomain assigned to the
!    processor.
!--------------------------------------------------------------------
      tot_pts = (size(latb,2)-1)*(size(lonb,1)-1)

!--------------------------------------------------------------------
!    determine if the current run is cold-starting this module. if a
!    restart file is present, then this is not a coldstart. in that case
!    set num_pts to tot_pts so that if cloud data is not available an
!    error message can be generated. if this is a coldstart, cloud data
!    will not be available until num_pts equals or exceeds tot_pts, so
!    continue processing without issuing an error message.
!--------------------------------------------------------------------
      if (file_exists ('INPUT/tracer_cld_amt.res') .or.  &
          file_exists ('INPUT/strat_cloud.res') ) then
        num_pts = tot_pts
      else
        num_pts = 0
      endif

!---------------------------------------------------------------------
!    obtain the tracer indices for the strat_cloud variables when
!    running gcm.
!---------------------------------------------------------------------
      if (Cldrad_control%do_strat_clouds) then
          nql = get_tracer_index ( MODEL_ATMOS, 'liq_wat' )
          nqi = get_tracer_index ( MODEL_ATMOS, 'ice_wat' )
          nqa = get_tracer_index ( MODEL_ATMOS, 'cld_amt' )

          if (do_rain) then !sjl
             nqr = get_tracer_index ( MODEL_ATMOS, 'rainwat' )
             if (nqr < 0 ) call error_mesg ('cloud_spec_mod', &
                'rainwat tracer not found, but do_rain is true', FATAL)
          end if
          if (do_snow) then !miz
             nqs = get_tracer_index ( MODEL_ATMOS, 'snowwat' )
             if (nqs < 0 ) call error_mesg ('cloud_spec_mod', &
                'snowwat tracer not found, but do_snow is true', FATAL)
          end if
          if (do_graupel) then !sjl
             nqg = get_tracer_index ( MODEL_ATMOS, 'graupel' )
             if (nqg < 0 ) call error_mesg ('cloud_spec_mod', &
                'graupel tracer not found, but do_graupel is true', FATAL)
          end if

          if (mpp_pe() == mpp_root_pe()) &
            write (logunit,'(a,3i4)') 'Stratiform cloud tracer ind&
                &ices: nql,nqi,nqa =',nql,nqi,nqa
          if (min(nql,nqi,nqa) <= 0)   &
             call error_mesg ('cloud_spec_mod', &
             'stratiform cloud tracer(s) not found', FATAL)
          if (nql == nqi .or. nqa == nqi .or. nql == nqa)   &
              call error_mesg ('cloud_spec_mod',  &
            'tracers indices cannot be the same (i.e., nql=nqi=nqa).', &
                                                              FATAL)
          nqn = get_tracer_index ( MODEL_ATMOS, 'liq_drp' )
          if (nqn /= NO_TRACER)  then
            Cldrad_control%do_liq_num = .true.
          else
            Cldrad_control%do_liq_num = .false.
          endif
          nqni = get_tracer_index ( MODEL_ATMOS, 'ice_num' )
          if (nqni /= NO_TRACER)  then
            Cldrad_control%do_ice_num = .true.
          else
            Cldrad_control%do_ice_num = .false.
          endif
      else
          Cldrad_control%do_liq_num = .false.
          Cldrad_control%do_ice_num = .false.
      endif

!BW   Cldrad_control%do_liq_num_iz = .true.
!BW   Cldrad_control%do_ice_num_iz = .true.

!---------------------------------------------------------------------
!    define the variables indicating that the cloud parameterization
!    control variables have been defined.
!---------------------------------------------------------------------
!BW   Cldrad_control%do_rh_clouds_iz = .true.
!BW   Cldrad_control%do_strat_clouds_iz = .true.
!BW   Cldrad_control%do_no_clouds_iz = .true.
!BW   Cldrad_control%do_uw_clouds_iz = .true.

!--------------------------------------------------------------------
!    include do_fu2007 in the cloudrad_control_type variable for use
!    in other modules.
!--------------------------------------------------------------------
      Cldrad_control%using_fu2007 = do_fu2007
!BW   Cldrad_control%using_fu2007_iz = .true.

!--------------------------------------------------------------------
!    save the radiative time step (needed when data override is active)
!--------------------------------------------------------------------
     !if (doing_data_override) then
          Radiation_time_step = set_time (rad_time_step, 0)
     !endif

!---------------------------------------------------------------------
!    mark the module initialized.
!---------------------------------------------------------------------
      module_is_initialized = .true.

!-------------------------------------------------------------------

end subroutine cloud_spec_init

!######################################################################
! <SUBROUTINE NAME="cloud_spec">
!  <OVERVIEW>
!    cloud_radiative_properties defines the cloud radiative properties
!    appropriate for the radiation options that are active.
!  </OVERVIEW>
!  <DESCRIPTION>
!    cloud_radiative_properties defines the cloud radiative properties
!    appropriate for the radiation options that are active.
!  </DESCRIPTION>
!  <TEMPLATE>
!   call cloud_spec (is, ie, js, je, lat, z_half, z_full, Rad_time,
!                       Atmos_input, &
!                       Surface, Cld_spec, Cloud_microphys,  &
!                       Moist_clouds, r)
!  </TEMPLATE>
!  <IN NAME="is,ie,js,je" TYPE="integer">
!   starting/ending subdomain i,j indices of data in
!                   the physics_window being integrated
!  </IN>
!  <IN NAME="Rad_time" TYPE="time_type">
!   time at which radiation calculation is to apply
!  </IN>
!  <INOUT NAME="Atmos_input" TYPE="atmos_input_type">
!    atmospheric input fields on model grid,
!  </INOUT>
!  <INOUT NAME="Cld_spec" TYPE="cld_specification_type">
!   cloud specification properties on model grid,
!  </INOUT>
!  <INOUT NAME="Lsc_microphys" TYPE="microphysics_type">
!   microphysical specification for large-scale
!                        clouds
!  </INOUT>
!  <INOUT NAME="Shallow_microphys" TYPE="microphysics_type">
!   microphysical specification for
!                        clouds associated with uw shallow convection
!  </INOUT>
!  <INOUT NAME="Surface" TYPE="Surface">
!   Surface boundary condition to radiation package
!  </INOUT>
!  <IN NAME="lsc_liquid_in" TYPE="real">
!   OPTIONAL: lsc cloud water mixing ratio
!  </IN>
!  <IN NAME="lsc_ice_in" TYPE="real">
!   OPTIONAL: cloud ice mixing ratio
!  </IN>
!  <IN NAME="lsc_area_in" TYPE="real">
!   OPTIONAL: fractional cloud area
!  </IN>
!  <IN NAME="r" TYPE="real">
!   OPTIONAL: model tracer fields on the current time step
!  </IN>
! </SUBROUTINE>
!
subroutine cloud_spec (is, ie, js, je, lat, land, z_half, z_full, Rad_time, &
                       press, pflux, temp, cloudtemp, cloudvapor, clouddeltaz, &
                       r, Cldrad_control, Cld_spec, Cloud_microphys, Aerosol, Moist_clouds_block)

!----------------------------------------------------------------------
!    cloud_spec specifies the cloud field seen by the radiation package.
!----------------------------------------------------------------------

integer,                      intent(in)             :: is, ie, js, je
real, dimension(:,:),         intent(in)             :: lat
real, dimension(:,:),         intent(in)             :: land
real, dimension(:,:,:),       intent(in)             :: z_half, z_full
type(time_type),              intent(in)             :: Rad_time
real, dimension(:,:,:),       intent(in)             :: press, pflux, &
                                                        temp, cloudtemp, &
                                                        cloudvapor, clouddeltaz
real, dimension(:,:,:,:),     intent(in)             :: r
type(cloudrad_control_type),  intent(in)             :: Cldrad_control
type(cld_specification_type), intent(inout)          :: Cld_spec
type(microphysics_type),      intent(inout), dimension(:), allocatable :: Cloud_microphys
type(aerosol_type),           intent(in)             :: Aerosol
type(clouds_from_moist_block_type), intent(in)       :: Moist_clouds_block

!-------------------------------------------------------------------

!---------------------------------------------------------------------
!   intent(in) variables:
!
!      is,ie,js,je       starting/ending subdomain i,j indices of data
!                        in the physics_window being integrated
!      lat               latitude of model points  [ radians ]
!      z_half            height asl at half levels [ m ]
!      z_full            height asl at full levels [ m ]
!      Rad_time          time at which radiation calculation is to apply
!                        [ time_type (days, seconds) ]
!
!   intent(inout) variables:
!
!      Atmos_input       atmospheric input fields on model grid,
!                        [ atmos_input_type ]
!      Surface           variables defining the surface albedo and land
!                        fraction
!                        [ surface_type ]
!      Cld_spec          variables on the model grid which define all or
!                        some of the following, dependent on the
!                        specific cloud parameterization: cloud optical
!                        paths, particle sizes, cloud fractions, cloud
!                        thickness, number of clouds in a column,
!                        and /or cloud type (high/mid/low, ice/liq or
!                        random/max overlap)
!                        [ cld_specification_type ]
!      Lsc_microphys     variables describing the microphysical proper-
!                        ties of the large-scale clouds
!                        [ microphysics_type ]
!   intent(in), optional variables:
!
!      lsc_liquid_in     cloud water mixing ratio (or specific humidity
!                        ????) [ non-dimensional ]
!      lsc_ice_in        cloud ice mixing ratio (or specific humidity
!                         ????) [ non-dimensional ]
!      lsc_area_in       fractional cloud area [ non-dimensional ]
!      r                 model tracer fields on the current time step
!
!----------------------------------------------------------------------

!---------------------------------------------------------------------
!   local variables:

      integer   :: ix, jx, kx
      logical   :: override
      logical   :: strat_data_found
      type(time_type) :: Data_time
      real, dimension (size (press,1), &
                       size (press,2), &
                       size (press,3)) :: rho
      integer   :: ncld, n
      integer   :: ncld_used
      character(len=64), dimension(size(Moist_clouds_block%Cloud_data,1)) :: scheme_names_used
      type(microphysics_type) :: Lsc_microphys

! locally define indices to make code more readable
      integer :: index_strat, index_shallow

!     indices for cloud schemes
!     Microphysics index (clouds actually used)

integer :: istrat, ishallow

!---------------------------------------------------------------------
!   local variables:
!
!        ix      number of grid points in x direction (on processor)
!        jx      number of grid points in y direction (on processor)
!        kx      number of model layers
!        rho     atmospheric density [ kg / m**3 ]
!        ncld    number of cloud schemes
!
!---------------------------------------------------------------------

      ncld = size(Moist_clouds_block%Cloud_data,1)

!---------------------------------------------------------------------
!    check for the presence of known cloud schemes
!---------------------------------------------------------------------
      ncld_used = 0

     !-------------------
     ! stratiform clouds
     !-------------------
      istrat = 0
      strat_data_found = .false.
      if (Cldrad_control%do_strat_clouds) then
       ! maybe this condition can be allowed if tracers are to be used for large-scale
        if (Moist_clouds_block%index_strat == 0) call error_mesg ('cloud_spec_mod', &
             'stratiform cloud properties not found when &
             &stratiform clouds requested', FATAL)
        ncld_used = 1
        scheme_names_used(ncld_used) = 'strat_cloud'
        istrat = ncld_used
        strat_data_found = .true.
      endif

     !-------------------------------------
     ! check for shallow cloud input data
     !-------------------------------------
      ishallow = 0
      if (Cldrad_control%do_uw_clouds) then
        if (Moist_clouds_block%index_uw_conv == 0) call error_mesg ('cloud_spec_mod',  &
                 'shallow cloud properties not found when &
                        &uw shallow clouds requested', FATAL)
        ncld_used = ncld_used+1
        scheme_names_used(ncld_used) = 'uw_conv'
        ishallow = ncld_used
      endif

      !! DEBUGGING !!
!----------------------------------------------------------------------
!    define model dimensions.
!----------------------------------------------------------------------
      ix = size(press,1)
      jx = size(press,2)
      kx = size(press,3)

!----------------------------------------------------------------------
!    allocate and initialize the arrays contained in the structures
!    used to specify the cloud amounts, types and locations and
!    the microphysical parameters.
!----------------------------------------------------------------------
      allocate(Cloud_microphys(ncld_used))
      do n = 1, ncld_used
        call Cloud_microphys(n)%alloc ( ix, jx, kx, &
                                scheme_names_used(n), Cldrad_control)
      enddo

      call Cld_spec%alloc (ix, jx, kx, Cldrad_control)

!---------------------------------------------------------------------
!    define the cloud_water, cloud_ice and cloud_area components of
!    Cld_spec.
!---------------------------------------------------------------------
      if (Moist_clouds_block%index_strat > 0 .and. .not. use_cloud_tracers_in_radiation) then
        index_strat = Moist_clouds_block%index_strat
        Cld_spec%cloud_ice   = Moist_clouds_block%Cloud_data(index_strat)%ice_amt
        Cld_spec%cloud_water = Moist_clouds_block%Cloud_data(index_strat)%liquid_amt
        Cld_spec%cloud_area  = Moist_clouds_block%Cloud_data(index_strat)%cloud_area
        if (Cldrad_control%do_liq_num) then
            Cld_spec%cloud_droplet = Moist_clouds_block%Cloud_data(index_strat)%droplet_number
        endif
        if (Cldrad_control%do_ice_num) then
            Cld_spec%cloud_ice_num = Moist_clouds_block%Cloud_data(index_strat)%ice_number
        endif
        Cld_spec%snow       = Moist_clouds_block%Cloud_data(index_strat)%snow
        Cld_spec%rain       = Moist_clouds_block%Cloud_data(index_strat)%rain
        Cld_spec%snow_size  = Moist_clouds_block%Cloud_data(index_strat)%snow_size
        Cld_spec%rain_size  = Moist_clouds_block%Cloud_data(index_strat)%rain_size
      endif

!----------------------------------------------------------------------
!    if a cloud scheme is activated (in contrast to running without any
!    clouds), call the appropriate subroutine to define the cloud
!    location, type, amount or whatever other arrays the particular
!    parameterization uses to specify its clouds. if the model is being
!    run with do_no_clouds = .true., exit from this routine, leaving
!    the cloud specification variables as they were initialized (to a
!    condition of no clouds).
!---------------------------------------------------------------------
      if (.not. Cldrad_control%do_no_clouds) then

!---------------------------------------------------------------------
!    if the rh diagnostic cloud scheme is active, call rh_clouds_amt
!    to define the needed cloud specification variables.
!---------------------------------------------------------------------
!BW     if (Cldrad_control%do_rh_clouds) then
!BW       call rh_clouds_amt (is, ie, js, je, press, lat,  &
!BW                           Cld_spec)
!BW     endif ! (do_rh_clouds)

!--------------------------------------------------------------------
!    if klein prognostic clouds are active, call strat_clouds_amt to
!    obtain the needed cloud specification variables.
!--------------------------------------------------------------------
        if (Cldrad_control%do_strat_clouds) then

!---------------------------------------------------------------------
!    if the gcm is being executed, call strat_cloud_avg to obtain the
!    appropriate (either instantaneous or time-averaged) values of
!    cloud water, cloud ice and cloud fraction. if the sa_gcm or the
!    standalone columns mode is being executed with the strat cloud
!    option, then values for the cloud water, cloud ice and when needed
!    cloud area have been input as optional arguments to this sub-
!    routine.
!---------------------------------------------------------------------
          if (Moist_clouds_block%index_strat > 0 .and. .not. use_cloud_tracers_in_radiation) then

            if (Cld_spec%cloud_area(1,1,1) == -99.) then
                Cld_spec%cloud_water(:,:,:) = r(:,:,:,nql)
                Cld_spec%cloud_ice  (:,:,:) = r(:,:,:,nqi)
                Cld_spec%cloud_area (:,:,:) = r(:,:,:,nqa)
                if (Cldrad_control%do_liq_num) then
                  Cld_spec%cloud_droplet (:,:,:) = r(:,:,:,nqn)
                endif
            endif

            if (Cld_spec%cloud_ice_num(1,1,1) == -99.) then
                if (Cldrad_control%do_ice_num) then
                  Cld_spec%cloud_ice_num (:,:,:) = r(:,:,:,nqni)
                endif
            endif

          else  ! (present (lsc_liquid_in))

              Cld_spec%cloud_water(:,:,:) = r(:,:,:,nql)
              Cld_spec%cloud_ice  (:,:,:) = r(:,:,:,nqi)
              Cld_spec%cloud_area (:,:,:) = r(:,:,:,nqa)
              if (Cldrad_control%do_liq_num) then
                Cld_spec%cloud_droplet (:,:,:) = r(:,:,:,nqn)
              endif
              if (Cldrad_control%do_ice_num) then
                Cld_spec%cloud_ice_num (:,:,:) = r(:,:,:,nqni)
              endif
          endif ! (present(lsc_liquid_in))


            if (do_rain) then !sjl
              Cld_spec%cloud_water(:,:,:) = Cld_spec%cloud_water(:,:,:)+r(:,:,:,nqr)
            end if
            if (do_snow) then !miz
              Cld_spec%cloud_ice(:,:,:) = Cld_spec%cloud_ice(:,:,:)+r(:,:,:,nqs)
            end if
            if (do_graupel) then !SJL
              Cld_spec%cloud_ice(:,:,:) = Cld_spec%cloud_ice(:,:,:)+r(:,:,:,nqg)
            end if

!---------------------------------------------------------------------
!    if the cloud input data is to be overriden, define the time slice
!    of data which is to be used. allocate storage for the cloud data.
!---------------------------------------------------------------------
          if (doing_data_override) then
            Data_time = Rad_time + Radiation_time_step

!---------------------------------------------------------------------
!    call data_override to retrieve the processor subdomain's cloud
!    water data from the override file. if the process fails, write
!    an error message; if it succeeds move the data for the current
!    physics window, into the appropriate Cld_spec% array.
!---------------------------------------------------------------------
            call data_override ('ATM', 'qlnew', Cld_spec%cloud_water,   &
                                Data_time, override=override,           &
                                is_in=is, ie_in=ie, js_in=js, je_in=je)
            if ( .not. override) then
              call error_mesg ('radiation_driver_mod', &
                'ql => cloud_water not overridden successfully', FATAL)
            endif

!---------------------------------------------------------------------
!    call data_override to retrieve the processor subdomain's cloud
!    ice data from the override file. if the process fails, write
!    an error message; if it succeeds move the data for the current
!    physics window, into the appropriate Cld_spec% array.
!---------------------------------------------------------------------
            call data_override ('ATM', 'qinew', Cld_spec%cloud_ice,   &
                                Data_time, override=override,         &
                                is_in=is, ie_in=ie, js_in=js, je_in=je)
            if ( .not. override) then
              call error_mesg ('radiation_driver_mod', &
                'qi => cloud_ice   not overridden successfully', FATAL)
            endif

!---------------------------------------------------------------------
!    call data_override to retrieve the processor subdomain's cloud
!    fraction data from the override file. if the process fails, write
!    an error message; if it succeeds move the data for the current
!    physics window, into the appropriate Cld_spec% array.
!---------------------------------------------------------------------
            call data_override ('ATM', 'qanew', Cld_spec%cloud_area,   &
                                Data_time, override=override,         &
                                is_in=is, ie_in=ie, js_in=js, je_in=je)
            if ( .not. override) then
              call error_mesg ('radiation_driver_mod', &
               'qa => cloud_area not overridden successfully', FATAL)
            endif
            strat_data_found = .true.
          endif ! (doing_override)

!---------------------------------------------------------------------
!    if values for the cloud variables have been successfully obtained,
!    call strat_clouds_amt to define the appropriate cloud specification
!    variables.
!---------------------------------------------------------------------
          if (strat_data_found) then
            call strat_clouds_amt (is, ie, js, je, Rad_time, &
                                   pflux, press, cloudtemp, &
                                   cloudvapor(:,:,:)/(1.0+cloudvapor(:,:,:)), &
                                   land, Cldrad_control, &
                                   Cld_spec, Cloud_microphys(istrat), Aerosol)

!----------------------------------------------------------------------
!    cloud data was not successfully obtained.
!    if this is not the coldstart step, write an error message and
!    stop execution.
!----------------------------------------------------------------------
          else
            if (num_pts >= tot_pts) then
              call error_mesg ('cloud_spec_mod',  &
                     'no strat cloud data available', FATAL)

!----------------------------------------------------------------------
!    if this is the coldstart step, retain the input values corres-
!    ponding to no clouds, increment the points counter, and continue.
!----------------------------------------------------------------------
            else
!$OMP ATOMIC UPDATE
              num_pts = num_pts + size(press,1)*size(press,2)
            endif
          endif

        else ! (do_strat_clouds)
!----------------------------------------------------------------------
!    only allocate part of the variables in this derived type
!    when strat_clouds is not active
!       (WHERE DOES THIS GET DEALLOCATED?)
!----------------------------------------------------------------------
          allocate(Lsc_microphys%conc_drop(id,jd,kmax))
          allocate(Lsc_microphys%conc_ice (id,jd,kmax))
          allocate(Lsc_microphys%size_drop(id,jd,kmax))
          allocate(Lsc_microphys%size_ice (id,jd,kmax))
          allocate(Lsc_microphys%size_rain(id,jd,kmax))
          Lsc_microphys%conc_drop = 0.0
          Lsc_microphys%conc_ice  = 0.0
          Lsc_microphys%size_drop = 0.0
          Lsc_microphys%size_ice  = 0.0
          Lsc_microphys%size_rain = 0.0

        endif ! (do_strat_clouds)

!--------------------------------------------------------------------
!    since uw_clouds may be active along with strat clouds
!    the associated properties are determined
!    outside of the above loop. these properties are placed in
!    Shallow_microphys.
!----------------------------------------------------------------------
        if (Cldrad_control%do_uw_clouds) then
          index_shallow = Moist_clouds_block%index_uw_conv
          call uw_clouds_amt (is, ie, js, je,  &
                           Moist_clouds_block%Cloud_data(index_shallow)%cloud_area,     &
                           Moist_clouds_block%Cloud_data(index_shallow)%liquid_amt,     &
                           Moist_clouds_block%Cloud_data(index_shallow)%ice_amt,        &
                           Moist_clouds_block%Cloud_data(index_shallow)%droplet_number, &
                           Moist_clouds_block%Cloud_data(index_shallow)%ice_number,     &
                           land, press, cloudtemp, Cldrad_control, Cloud_microphys(ishallow) )
        endif

!---------------------------------------------------------------------
!    call combine_cloud_properties to combine (if necessary) the cloud
!    properties from multiple cloud types (large-scale,
!    uw shallow) into a single set for use by the radiation package.
!    this is only needed when microphysically-based properties are
!    present, and when either strat clouds and / or uw
!    shallow clouds is activated.
!---------------------------------------------------------------------
!BW     if (Cldrad_control%do_sw_micro .or. Cldrad_control%do_lw_micro) then
          if (ncld_used > 0) then
              call combine_cloud_properties ( is, js,  &
                                             temp(:,:,1), Rad_time, &
                                             Cldrad_control, Cloud_microphys, Cld_spec)
          endif
!BW     endif


!--------------------------------------------------------------------
!    if microphysics is active and strat_clouds is not, define the water
!    paths (in units of kg / m**2).  if strat_clouds is active, these
!    values will have already been defined. when microphysics is active,
!    define the effective sizes for the liquid and ice particles.
!--------------------------------------------------------------------
!BW   if (Cldrad_control%do_lw_micro .or.    &
!BW       Cldrad_control%do_sw_micro)  then
    !BW if (.not. Cldrad_control%do_strat_clouds) then

        if (Cldrad_control%do_strat_clouds) then
          Cld_spec%reff_liq_micro = Cloud_microphys(istrat)%size_drop
          Cld_spec%reff_ice_micro = Cloud_microphys(istrat)%size_ice
        endif

      endif  !  (.not. do_no_clouds)
!---------------------------------------------------------------------


end subroutine cloud_spec


!######################################################################

subroutine cloud_spec_end (Cldrad_control)

type(cloudrad_control_type), intent(in) :: Cldrad_control

!---------------------------------------------------------------------
!    cloud_spec_end is the destructor for cloud_spec_mod.
!---------------------------------------------------------------------

!--------------------------------------------------------------------
!    close the modules that were initialized by this module.
!--------------------------------------------------------------------
      if (.not. Cldrad_control%do_no_clouds) then

!-------------------------------------------------------------------
!    rh-based diagnostic clouds.
!-------------------------------------------------------------------
!BW     if (Cldrad_control%do_rh_clouds) then
!BW       call rh_based_clouds_end

!------------------------------------------------------------------
!    cloud types which may coexist must be processed outside of if loop
!------------------------------------------------------------------
!BW     else

          if (Cldrad_control%do_strat_clouds) then
            call strat_clouds_W_end (Cldrad_control)
          endif
          if (Cldrad_control%do_uw_clouds) then
            call uw_clouds_W_end
          endif
!BW     endif

        if (Cldrad_control%do_stochastic_clouds) then
          call random_number_streams_end
        endif

      endif  ! (not do_no_clouds)

!--------------------------------------------------------------------
!    mark the module as no longer initialized.
!--------------------------------------------------------------------
      module_is_initialized = .false.

!--------------------------------------------------------------------

end subroutine cloud_spec_end

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!                    PRIVATE SUBROUTINES
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!###################################################################
! <SUBROUTINE NAME="combine_cloud_properties">
!  <OVERVIEW>
!    combine_cloud_properties produces cloud specification property
!    arrays for the total cloud field in each grid box, using as input
!    the specification of the component cloud types that may be present
!    (large-scale, mesoscale and cell-scale).
!  </OVERVIEW>
!  <DESCRIPTION>
!    combine_cloud_properties produces cloud specification property
!    arrays for the total cloud field in each grid box, using as input
!    the specification of the component cloud types that may be present
!    (large-scale).
!  </DESCRIPTION>
!  <TEMPLATE>
!   call combine_cloud_properties (Lsc_microphys, Cld_spec)
!  </TEMPLATE>
!  <INOUT NAME="Cld_spec" TYPE="cld_specification_type">
!   cloud specification properties on model grid,
!  </INOUT>
!  <IN NAME="Lsc_microphys" TYPE="microphysics_type">
!   microphysical specification for large-scale
!                        clouds
!  </IN>
!  </IN>
!  <IN NAME="Shallow_microphys" TYPE="microphysics_type">
!   microphysical specification for
!                        clouds associated with uw shallow convection
!  </IN>
! </SUBROUTINE>
!
subroutine combine_cloud_properties (is, js, temp, Rad_time, &
                                     Cldrad_control, Cloud_microphys, Cld_spec)

!----------------------------------------------------------------------
!    combine_cloud_properties produces cloud specification property
!    arrays for the total cloud field in each grid box, using as input
!    the specification of the component cloud types that may be present
!    (large-scale, uw shallow).
!----------------------------------------------------------------------

integer, intent(in)  :: is, js
real, dimension(:,:), intent(in) :: temp
type(time_type), intent(in) :: Rad_time
type(cloudrad_control_type),            intent(in)    :: Cldrad_control
type(microphysics_type),  dimension(:), intent(in)    :: Cloud_microphys
type(cld_specification_type),           intent(inout) :: Cld_spec

!----------------------------------------------------------------------
!   intent(in) variables:
!
!       Lsc_microphys  microphysical specification for large-scale
!                      clouds
!                      [ microphysics_type ]
!       Shallow_microphys
!                      microphysical specification for
!                      clouds associated with uw shallow convection
!                      [ microphysics_type ]
!
!    intent(inout) variables:
!
!       Cld_spec       variables on the model grid which define all or
!                      some of the following, dependent on the specific
!                      cloud parameterization: cloud optical paths,
!                      particle sizes, cloud fractions, cloud thickness,
!                      number of clouds in a column, and /or cloud type
!                      (high/mid/low, ice/liq or random/max overlap)
!                      [ cld_specification_type ]
!
!---------------------------------------------------------------------

      type(randomNumberStream),   &
                    dimension(size(Cld_spec%camtsw,1),   &
                              size(Cld_spec%camtsw,2)) :: streams
      real, &
                    dimension(size(Cld_spec%camtsw,1),   &
                              size(Cld_spec%camtsw,2),   &
                              size(Cld_spec%camtsw,3),   &
                              Cldrad_control%num_lw_cloud_bands+ &
                              Cldrad_control%num_sw_cloud_bands) :: &
                                                     randomNumbers
      integer :: nn, nsubcols

      integer :: i, j, k, n, ncld
      integer :: istrat, ishallow

!---------------------------------------------------------------------
!    total-cloud specification properties need be defined only when
!    strat_cloud and/or uw shallow clouds are active.
!---------------------------------------------------------------------

      ncld = size(Cloud_microphys,1)

  ! indices for cloud types in microphysics type
      istrat=0
      ishallow=0
      do n = 1, ncld
         if (trim(Cloud_microphys(n)%scheme_name) == 'strat_cloud') istrat = n
         if (trim(Cloud_microphys(n)%scheme_name) == 'uw_conv')     ishallow = n
      enddo

!----------------------------------------------------------------------
!    define the random overlap cloud fraction as the sum of the
!    fractions of all cloud schemes
!---------------------------------------------------------------------

      ! general case: does not reproduce ulm version
      Cld_spec%crndlw = 0.0
      do n = 1, ncld
        Cld_spec%crndlw = Cld_spec%crndlw + Cloud_microphys(n)%cldamt
      enddo

!---------------------------------------------------------------------
!    randomly-overlapped clouds are being assumed for 
!    strat cloud module clouds. set the max overlap cloud fraction to
!    zero, be certain that the random overlap fraction is .le. 1. after
!    the summing of the component cloud fractions, and define the total
!    cloud fraction to be used by the sw code.
!---------------------------------------------------------------------
      Cld_spec%cmxolw = 0.0
      Cld_spec%crndlw = MIN (Cld_spec%crndlw, 1.00)
      Cld_spec%camtsw = Cld_spec%crndlw

!--------------------------------------------------------------------
!    if stochastic clouds are being used, define the cloud type to be
!    seen by the radiation code in each stochastic subcolumn.
!--------------------------------------------------------------------
      if (Cldrad_control%do_stochastic_clouds) then
        nsubcols = Cldrad_control%num_sw_cloud_bands + &
                   Cldrad_control%num_lw_cloud_bands

!--------------------------------------------------------------------
!   assign either a 1 or a 0 to each subcolumn indicating whether
!   lsc cloud is present or not.
!--------------------------------------------------------------------
        if (istrat > 0) then
          do n=1,nsubcols
            if ( n > Cldrad_control%num_sw_cloud_bands) then
              nn = n - Cldrad_control%num_sw_cloud_bands
            else
              nn = n + Cldrad_control%num_lw_cloud_bands
            endif
            do k=1,size(Cld_spec%camtsw,3) ! Levels
               do j=1,size(Cld_spec%camtsw,2) ! Lons
                do i=1,size(Cld_spec%camtsw,1) ! Lats
                  if (Cloud_microphys(istrat)%stoch_cldamt(i,j,k,nn) > 0.) then
                    !----------------------------------------------
                    ! fill it in with the large-scale cloud values
                    !-----------------------------------------------
                     Cld_spec%stoch_cloud_type(i,j,k,n) = istrat
                    !Cld_spec%stoch_cloud_type(i,j,k,n) = 1
                  else
                     Cld_spec%stoch_cloud_type(i,j,k,n) = 0
                  endif
                enddo
              enddo
            enddo
          enddo
        endif ! strat_cloud

!----------------------------------------------------------------------
!    compare the uw shallow cloud amount to a random number, and replace
!    the large-scale cloud or clear sky previously
!    assigned in each subcolumn with an assignment of uw shallow cloud
!    when the number is less than the cloud fraction. use the maximum
!    overlap assumption. treat the random number as the location with
!    the PDF of total water. uw shallow clouds are at the top of this
!    PDF, then large-scale clouds and clear sky.
!------------------------------------------------------------
        if (Cldrad_control%do_uw_clouds) then
          call get_random_number_streams (is, js, Rad_time, temp, streams, perm=2)

!----------------------------------------------------------------------
!    get the random numbers to do both sw and lw at oncer.
!----------------------------------------------------------------------
          do j=1,size(Cld_spec%camtsw,2) ! Lons
            do i=1,size(Cld_spec%camtsw,1) ! Lats
              call getRandomNumbers (streams(i,j), randomNumbers(i,j,1,:))
            end do
          end do

!----------------------------------------------------------------------
!    here is maximum overlap. we use a 3D arrary for the random numbers
!    for flexibility.
!----------------------------------------------------------------------
          do k=2,size(Cld_spec%camtsw,3)
            randomNumbers(:,:,k,:) = randomNumbers(:,:,1,:)
          end do

!----------------------------------------------------------------------
!    assign cloud type, band by band
!----------------------------------------------------------------------
          do n=1,nsubcols
            where( randomNumbers(:,:,:,n) > &
                   (1. - Cloud_microphys(ishallow)%cldamt)) &
                 ! assign uw shallow cloud
                   Cld_spec%stoch_cloud_type(:,:,:,n) = ishallow
                 ! Cld_spec%stoch_cloud_type(:,:,:,n) = 4
          enddo
        endif

!---------------------------------------------------------------------
!     define the cloud amount in each stochastic subcolumn to be either
!     1.0 if cloud is present, or 0.0 if no cloud exists.
!---------------------------------------------------------------------
        do n=1,Cldrad_control%num_sw_cloud_bands
          do k=1,size(Cld_spec%camtsw,3) ! Levels
            do j=1,size(Cld_spec%camtsw,2) ! Lons
              do i=1,size(Cld_spec%camtsw,1) ! Lats
                if (Cld_spec%stoch_cloud_type(i,j,k,n) /= 0) then
                  Cld_spec%camtsw_band(i,j,k,n) = 1.0
                else
                  Cld_spec%camtsw_band(i,j,k,n) = 0.0
                endif
              end do
            end do
          end do
        end do

        do n=1,Cldrad_control%num_lw_cloud_bands
          nn = Cldrad_control%num_sw_cloud_bands + n
          do k=1,size(Cld_spec%camtsw,3) ! Levels
            do j=1,size(Cld_spec%camtsw,2) ! Lons
              do i=1,size(Cld_spec%camtsw,1) ! Lats
                if (Cld_spec%stoch_cloud_type(i,j,k,nn) /= 0) then
                  Cld_spec%crndlw_band(i,j,k,n) = 1.0
                else
                  Cld_spec%crndlw_band(i,j,k,n) = 0.0
                endif
              end do
            end do
          end do
        end do
      endif  ! (do_stochastic)


!-------------------------------------------------------------------


end subroutine combine_cloud_properties


!#################################################################


                end module cloud_spec_mod
