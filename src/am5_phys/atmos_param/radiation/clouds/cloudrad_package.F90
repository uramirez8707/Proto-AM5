                 module cloudrad_package_mod
! Contact:  Fei.Liu@noaa.gov - fil
!  Module that supplies cloud radiative properties

! shared modules:

use mpp_mod,                  only: input_nml_file
use fms_mod,                  only: fms_init, &
                                    write_version_number, mpp_pe, &
                                    mpp_root_pe, stdlog, &
                                    check_nml_error, error_mesg,   &
                                    FATAL
use time_manager_mod,         only: time_type, time_manager_init

! cloud radiation modules:

use cloudrad_types_mod,       only: cld_specification_type, &
                                    cldrad_properties_type, &
                                    microrad_properties_type, &
                                    microphysics_type, &
                                    cloudrad_control_type, &
                                    microrad_properties_alloc

use cloudrad_diagnostics_mod, only: cloudrad_diagnostics_init, &
                                    cloudrad_netcdf, &
                                    cloudrad_diagnostics_end

use microphys_rad_mod,        only: lwemiss_calc, comb_cldprops_calc, &
                                    microphys_rad_init, &
                                    microphys_rad_end, &
                                    microphys_lw_driver, &
                                    microphys_sw_driver

!--------------------------------------------------------------------

implicit none
private

!--------------------------------------------------------------------
!    cloudrad_package_mod computes cloud radiative properties consistent
!    with the activated radiation package options and returns them to
!    radiation_driver_mod.
!--------------------------------------------------------------------


!---------------------------------------------------------------------
!----------- version number for this module --------------------------

character(len=128)  :: version =  '$Id$'
character(len=128)  :: tagname =  '$Name$'


!---------------------------------------------------------------------
!-------  interfaces --------

public          &
         cloudrad_package_init, cloud_radiative_properties, &
         cloudrad_package_end


private          &
!  called from cloud_radiative_properties:
         combine_cloud_properties

!---------------------------------------------------------------------
!-------- namelist  ---------

character(len=16)  :: microphys_form =  'predicted' ! level of microphysics
                                                    ! being used; either
                                                    ! 'none', or 'predicted'
real               :: min_cld_drop_rad = 4.2    ! smallest allowable
                                                ! cloud drop radius
                                                ! (microns) allowed in
                                                ! slingo scheme
real               :: max_cld_drop_rad = 16.6   ! largest allowable
                                                ! cloud drop radius
                                                ! (microns) allowed in
                                                ! slingo scheme
real               :: min_cld_ice_size = 18.6   ! smallest allowable
                                                ! cloud ice size
                                                ! (microns) allowed in
                                                ! fu-liou scheme
real               :: max_cld_ice_size = 130.2  ! largest allowable
                                                ! cloud ice size
                                                ! (microns) allowed in
                                                ! fu-liou scheme

 logical     ::  do_ica_calcs=.false.           ! do independent column
                                                ! calculations when sto-
                                                ! chastic clouds are
                                                ! active ?


namelist /cloudrad_package_nml /     &
                               min_cld_drop_rad, max_cld_drop_rad, &
                               min_cld_ice_size, max_cld_ice_size, &
                               microphys_form, do_ica_calcs

!----------------------------------------------------------------------
!----  public data -------


!----------------------------------------------------------------------
!----  private data -------

logical :: module_is_initialized = .false.  ! module initialized?



!----------------------------------------------------------------------
!----------------------------------------------------------------------



                          contains



!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!                     PUBLIC SUBROUTINES
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!#####################################################################
subroutine cloudrad_package_init (pref, lonb, latb, axes, Time, &
                                  Cldrad_control)

!---------------------------------------------------------------------
!    cloudrad_package_init is the constructor for cloudrad_package_mod.
!----------------------------------------------------------------------

real,    dimension(:,:),     intent(in)    ::  pref
real,    dimension(:,:),     intent(in)    ::  lonb, latb
integer, dimension(4),       intent(in)    ::  axes
type(time_type),             intent(in)    ::  Time
type(cloudrad_control_type), intent(inout) ::  Cldrad_control

!---------------------------------------------------------------------
!   intent(in) variables:
!
!       pref      array containing two reference pressure profiles
!                 for use in defining transmission functions [ Pa ]
!       lonb      2d array of model longitudes on cell corners[ radians ]
!       latb      2d array of model latitudes at cell corners [radians]
!       axes      diagnostic variable axes
!       Time      current time [time_type(days, seconds)]
!
!----------------------------------------------------------------------

!---------------------------------------------------------------------
!   local variables:

      integer         :: io, ierr, logunit

!---------------------------------------------------------------------
!   local variables:
!
!      io       error status returned from io operation
!      ierr     error code
!
!--------------------------------------------------------------------

!-------------------------------------------------------------------
!    if routine has already been executed, exit.
!---------------------------------------------------------------------
      if (module_is_initialized) return

!---------------------------------------------------------------------
!    verify that modules used by this module that are not called later
!    have already been initialized.
!---------------------------------------------------------------------
      call fms_init
      call time_manager_init

!---------------------------------------------------------------------
!    read namelist.
      read (input_nml_file, nml=cloudrad_package_nml, iostat=io)
      ierr = check_nml_error(io,"cloudrad_package_nml")
!---------------------------------------------------------------------
!    write version number and namelist to logfile.
!---------------------------------------------------------------------
      call write_version_number (version, tagname)
      logunit = stdlog()
      if (mpp_pe() == mpp_root_pe() )    &
                       write (logunit, nml=cloudrad_package_nml)

!---------------------------------------------------------------------
!    define flag indicating if ICA calculations being done.
!    THIS MAY NEED TO BE DEFINED ELSEWHERE
!---------------------------------------------------------------------
      Cldrad_control%do_ica_calcs = do_ica_calcs
!BW   Cldrad_control%do_ica_calcs_iz = .true.

!-------------------------------------------------------------------
!    verify that the component cloud scheme elements of Cldrad_control
!    have been defined.
!-------------------------------------------------------------------
!BW   if (Cldrad_control%do_rh_clouds_iz .and.  &
!BW   if (Cldrad_control%do_strat_clouds_iz .and.  &
!BW       Cldrad_control%do_no_clouds_iz .and.  &
!BW       Cldrad_control%do_uw_clouds_iz .and.  &
!BW       Cldrad_control%do_donner_deep_clouds_iz ) then
!BW   else
!BW     call error_mesg ('cloudrad_package_mod', &
!BW      'Cldrad_control%do_{some cloud type} not yet been defined', &
!BW                                                              FATAL)
!BW   endif

!-------------------------------------------------------------------
!    define variables which denote whether microphysically-based cloud
!    radiative properties will be defined for the longwave and shortwave
!    radiation calculations. if esf sw is being used, then do_sw_micro
!    will be .true.; if number of lw bands > 1 or if esf sw is being
!    used along with strat clouds or diag clouds or specified clouds,
!    then do_lw_micro will be .true..
!----------------------------------------------------------------------
     !if (Lw_control%do_lwcldemiss) then
      if (Cldrad_control%num_lw_cloud_bands > 1) then
        Cldrad_control%do_lw_micro = .true.
      endif
!BW   if (Sw_control%do_esfsw) then
        Cldrad_control%do_sw_micro = .true.   ! always true
        if (Cldrad_control%do_strat_clouds) then
          Cldrad_control%do_lw_micro = .true.
        endif
!BW   endif

!--------------------------------------------------------------------
!    mark the logical controls as initialized.
!--------------------------------------------------------------------
!BW   Cldrad_control%do_lw_micro_iz = .true.
!BW   Cldrad_control%do_sw_micro_iz = .true.

!----------------------------------------------------------------------
!    define the microphysical use desired for this experiment. different
!    levels may be used for the lw and sw parameterizations.
!----------------------------------------------------------------------
      if (trim(microphys_form) == 'predicted') then

!---------------------------------------------------------------------
!    if microphys_form asks for predicted microphysics, then either
!    strat or donner deep clouds must be activated, and either one or
!    both of do_lw_micro and do_sw_micro must be .true..  if these
!    conditions are met, set do_pred_cld_microphys to .true..
!    if  neither do_lw_micro or do_sw_micro are .true.
!    or if a different cloud scheme has been activated, stop execution
!    with an error message.
!---------------------------------------------------------------------
        if (Cldrad_control%do_strat_clouds .or. &
            Cldrad_control%do_uw_clouds) then
          if (Cldrad_control%do_sw_micro .and. &
              Cldrad_control%do_lw_micro) then
            Cldrad_control%do_pred_cld_microphys = .true.
          else
            call error_mesg( 'cloudrad_package_mod',  &
             ' not using microphysics -- set microphys_form '//&
                                                    'to none.', FATAL)
          endif
        else
          call error_mesg( 'cloudrad_package_mod',  &
                    ' predicted microphys not available with this '//&
                                                 'cloud scheme.', FATAL)
        endif

!---------------------------------------------------------------------
!    if no microphysics is requested, make sure that donner_deep clouds
!    has not been requested (must use predicted cloud microphysics
!    for that scheme -- all others can be run without microphysics).
!    also verify that the lw and sw schemes requested are not micro-
!    physically_based.  if not ok, write an error message and stop.
!---------------------------------------------------------------------
      else if (trim(microphys_form) == 'none') then
        if (Cldrad_control%do_uw_clouds) then
          call error_mesg( 'cloudrad_package_mod',  &
            ' use predicted microphys with donner or uw clouds.', FATAL)
        else
          if (Cldrad_control%do_sw_micro .or.    &
              Cldrad_control%do_lw_micro) then
            call error_mesg ('cloudrad_package_mod', &
               'must specify microphys_form when using microphysica'//&
                'lly-based cld rad scheme', FATAL)
          endif
        endif

!----------------------------------------------------------------------
!    error condition.
!----------------------------------------------------------------------
      else
        call error_mesg( 'cloudrad_package_mod',  &
           ' microphys_form is not an acceptable value.', FATAL)
      endif

!---------------------------------------------------------------------
!    define variables indicating that the desired cloud microphysics
!    type control variables have been defined.
!---------------------------------------------------------------------
!BW   Cldrad_control%do_bulk_microphys_iz = .true.
!BW   Cldrad_control%do_pred_cld_microphys_iz = .true.
!BW   Cldrad_control%do_presc_cld_microphys_iz = .true.

!---------------------------------------------------------------------
!    if a microphysically-based scheme is being used, initialize the
!    microphys_rad module.
!---------------------------------------------------------------------
      if (Cldrad_control%do_pred_cld_microphys) then
        call microphys_rad_init ( min_cld_drop_rad, max_cld_drop_rad, &
                                  min_cld_ice_size, max_cld_ice_size, &
                                  axes, Time, lonb, latb, Cldrad_control )
      endif

!-------------------------------------------------------------------
!    when running the gcm or the standalone model with a cloud scheme,
!    call cloudrad_diagnostics_init to initialize the netcdf diagnostics
!    associated with the cloudrad package.
!-------------------------------------------------------------------
      if (.not. Cldrad_control%do_no_clouds) then
        call cloudrad_diagnostics_init (min_cld_drop_rad,   &
                                        max_cld_drop_rad, &
                                  min_cld_ice_size, max_cld_ice_size, &
                                        axes, Time, &
                                  Cldrad_control)
      endif

!--------------------------------------------------------------------
!    mark the module initialized.
!--------------------------------------------------------------------
      module_is_initialized= .true.

!--------------------------------------------------------------------


end subroutine cloudrad_package_init


!####################################################################
subroutine cloud_radiative_properties (is, ie, js, je, Rad_time,   &
                                       Time, Time_next,  &
                                       cosz, tsfc, press, pflux, temp, &
                                       deltaz, cloudtemp, cloudvapor, &
                                       clouddeltaz, Cldrad_control, Cld_spec, &
                                       Cloud_microphys, Model_microphys, &
                                       crndlw, cmxolw, emrndlw, emmxolw, &
                                       camtsw, cldsct, cldext, cldasymm)

!----------------------------------------------------------------------
!    cloud_radiative_properties defines the cloud radiative properties
!    appropriate for the radiation options that are active.
!----------------------------------------------------------------------

!----------------------------------------------------------------------
integer,                      intent(in)               :: is, ie, js, je
type(time_type),              intent(in)               :: Rad_time, Time, &
                                                          Time_next
real, dimension(:,:),         intent(in)               :: cosz, tsfc
real, dimension(:,:,:),       intent(in)               :: press, pflux, &
                                                          temp, deltaz, &
                                                          cloudtemp,  &
                                                          cloudvapor, &
                                                          clouddeltaz
type(cloudrad_control_type),  intent(in)               :: Cldrad_control
type(cld_specification_type), intent(inout)            :: Cld_spec
type(microphysics_type),      intent(in), dimension(:) :: Cloud_microphys
type(microphysics_type),      intent(inout)            :: Model_microphys
real, dimension(:,:,:,:),     intent(out)              :: crndlw, camtsw
real, dimension(:,:,:),       intent(out)              :: cmxolw
real, dimension(:,:,:,:,:),   intent(out)              :: emrndlw, emmxolw
real, dimension(:,:,:,:,:),   intent(out)              :: cldsct, cldext, cldasymm
!-------------------------------------------------------------------

!---------------------------------------------------------------------
!   intent(in) variables:
!
!      is,ie,js,je       starting/ending subdomain i,j indices of data
!                        in the physics_window being integrated
!      Time_next         time on next timestep, used as stamp for
!                        diagnostic output [ time_type (days, seconds) ]
!      Time              The current time.  [time_type]
!      Astro             astronomical properties needed by radiation
!                        package
!                        [ astronomy_type ]
!      Atmos_input       atmospheric input fields on model grid,
!                        [ atmos_input_type ]
!      Cld_spec          cloud specification properties on model grid,
!                        [ cld_specification_type ]
!      Lsc_microphys     microphysical specification for large-scale
!                        clouds
!                        [ microphysics_type ]
!      Shallow_microphys
!                        microphysical specification for
!                        clouds associated with uw shallow convection
!                        [ microphysics_type ]
!
!   intent(out) variables:
!
!      emrndlw           longwave cloud emissivity for random overlap clouds by band
!      emmxolw           longwave cloud emissivity for maximum overlap clouds by band
!      cldsct, cldext, cldasymm
!                        shortwave cloud radiative properties
!
!----------------------------------------------------------------------

!---------------------------------------------------------------------
!   local variables:

      type(cldrad_properties_type)   :: Cldrad_props
      type(microrad_properties_type), &
                      dimension(size(Cloud_microphys(:))) :: &
                                        Microrad_props
      integer  ::   ix, jx, kx, nb, np, n
      logical  ::   donner_flag_uw = .false.
      integer  ::   strat_index, shallow_index


!---------------------------------------------------------------------
!   local variables:
!
!       Lscrad_props   cloud radiative properties for the large-scale
!                      clouds
!                      [ microrad_properties_type ]
!       Shallowrad_props
!                      cloud radiative properties for
!                      clouds associated with uw shallow convection
!                      [ microrad_properties_type ]
!       ix             x dimension of current physics window
!       jx             y dimension of current physics window
!       kx             z dimension of current physics window
!       donner_flag    optional argument to microphys_rad module
!                      indicating that the meso or cell cloud compon-
!                      ent is currently being processed. needed because
!                      a different fu ice parameterization is used in
!                      these cases than is used for large-scale clouds.
!
!---------------------------------------------------------------------

!--------------------------------------------------------------------
!    define the dimensions of the current physics window.
!---------------------------------------------------------------------
      ix = size(Cld_spec%camtsw,1)
      jx = size(Cld_spec%camtsw,2)
      kx = size(Cld_spec%camtsw,3)

!--------------------------------------------------------------------
!    allocate and initialize the cloud radiative property arrays
!---------------------------------------------------------------------
      call Cldrad_props%alloc (ix, jx, kx, Cldrad_control)

      do n = 1, size(Microrad_props,1)
         call microrad_properties_alloc (Microrad_props(n), ix, jx, kx, Cldrad_control)
        !Microrad_props(n)%alloc (ix, jx, kx, Cldrad_control)
      enddo

!----------------------------------------------------------------------
!    if a cloud scheme is activated (in contrast to running without any
!    clouds), call the microphysically-based
!    modules to define the cloud lw and sw radiative properties. if the
!    model is being run with do_no_clouds = .true., exit from this
!    routine, leaving the cloud radiative property variables as they
!    were initialized (to values compatible to the non-existence of
!    cloudiness).
!---------------------------------------------------------------------
      if (.not. Cldrad_control%do_no_clouds) then

        ! define indexing of derived-type cloud properties
        strat_index = 0
        shallow_index = 0
        do n = 1, size(Cloud_microphys(:))
           if (trim(Cloud_microphys(n)%scheme_name) == 'strat_cloud') strat_index = n
           if (trim(Cloud_microphys(n)%scheme_name) == 'uw_conv')     shallow_index = n
           ! copy cloud scheme names
           Microrad_props(n)%scheme_name = trim(Cloud_microphys(n)%scheme_name)
        enddo

!RSH : the if loop is needed in case radiation is being run without seeing
!      strat_clouds
if (Cldrad_control%do_strat_clouds) then

        if (Cldrad_control%do_lw_micro) then
          if (Cldrad_control%do_ica_calcs) then
            call microphys_lw_driver (is, ie, js, je, Cldrad_control, &
                                      Cloud_microphys(strat_index),  &
                                      Cloud_rad_props=Cldrad_props)
          else
            call microphys_lw_driver (is, ie, js, je, Cldrad_control, &
                                      Cloud_microphys(strat_index),  &
                                      Micro_rad_props=Microrad_props(strat_index))
          endif
        endif
!BW     if (Cldrad_control%do_sw_micro) then
          if (Cldrad_control%do_ica_calcs) then
            call microphys_sw_driver (is, ie, js, je, Cldrad_control, &
                                      Cloud_microphys(strat_index),  &
                                      Cloud_rad_props=Cldrad_props)
          else
            call microphys_sw_driver (is, ie, js, je, Cldrad_control, &
                                      Cloud_microphys(strat_index),  &
                                      Micro_rad_props=Microrad_props(strat_index))
          endif

!RSH adds:
endif
!--------------------------------------------------------------------
!    if the uw shallow convection scheme is active, obtain the cloud
!    radiative properties associated with its clouds. only micro-
!    physically-based properties are available.
!    NOTE FOR NOW:
!   the optional argument  donner_flag is set to .false. when processing
!   shallow clouds. the ice cloud radiative properties are obtained from
!    the parameterization used by strat_cloud (effective size), rather
!    than that used by donner_deep (generalized effective size).
!----------------------------------------------------------------------
       if (Cldrad_control%do_uw_clouds) then
         if (shallow_index == 0) call error_mesg( 'cloudrad_package_mod',  &
           'shallow cloud scheme not found when do_uw_clouds = true', FATAL)
         donner_flag_uw = .false.
         call microphys_lw_driver (is, ie, js, je, Cldrad_control,  &
                                   Cloud_microphys(shallow_index),  &
                                   Micro_rad_props=Microrad_props(shallow_index), &
                                   donner_flag=donner_flag_uw)
         call microphys_sw_driver (is, ie, js, je, Cldrad_control,  &
                                   Cloud_microphys(shallow_index),  &
                                   Micro_rad_props=Microrad_props(shallow_index), &
                                   donner_flag=donner_flag_uw)
        endif
!RSH  endif ! ( .not. do_no_clouds)

!---------------------------------------------------------------------
!    call combine_cloud_properties to define a set of total-cloud cloud
!    radiative properties. if donner deep and / or uw shallow clouds
!    are active, this requires the combination of the cloud properties
!    associated with the different types of cloud present (large-scale,
!    donner meso and  cell, uw shallow). for other schemes the
!    total-cloud values are simply the large-scale cloud values.
!    this procedure is only needed when microphysically-based properties
!    are being used.
!---------------------------------------------------------------------
      if (.not. Cldrad_control%do_ica_calcs) then
        if (Cldrad_control%do_sw_micro  .or. &
            Cldrad_control%do_lw_micro) then
          call combine_cloud_properties (is, js, Rad_time, Time_next, &
                                         deltaz, Cldrad_control, Cld_spec, &
                                         Cloud_microphys, Microrad_props, &
                                         Cldrad_props)
        endif
      endif

   endif ! ( .not. do_no_clouds)

!----------------------------------------------------------------------
!    call lwemiss_calc to compute lw emissivity from the absorption
!    coefficient when a microphysically-based lw emissivity scheme
!    is being used.
!----------------------------------------------------------------------
      if (Cldrad_control%do_lw_micro) then
        call lwemiss_calc (clouddeltaz,   &
                           Cldrad_props%abscoeff, Cldrad_props%cldemiss)
        Cldrad_props%emmxolw = Cldrad_props%cldemiss
        Cldrad_props%emrndlw = Cldrad_props%cldemiss
      endif

!-------------------------------------------------------------------
!    if running the gcm or feedback program with a cloud scheme active,
!    call cloudrad_netcdf to generate netcdf output fields.
!-------------------------------------------------------------------
      if (.not. Cldrad_control%do_no_clouds) then
          call cloudrad_netcdf (is, js, Time, Time_next, cosz, &
                                tsfc, press, pflux, temp, deltaz, &
                                cloudtemp, cloudvapor, &
                                Cldrad_control, Cloud_microphys, &
                                Microrad_props,  Cldrad_props, &
                                Cld_spec, Model_microphys)
      endif   ! (do_no_clouds)

!--------------------------------------------------------------------
!    copy output cloud properties
!--------------------------------------------------------------------
      cmxolw   = Cld_spec%cmxolw
      if (Cldrad_control%do_stochastic_clouds) then
         crndlw = Cld_spec%crndlw_band
         camtsw = Cld_spec%camtsw_band
      else
         crndlw(:,:,:,1) = Cld_spec%crndlw
         camtsw(:,:,:,1) = Cld_spec%camtsw
      endif

      emrndlw  = Cldrad_props%emrndlw
      emmxolw  = Cldrad_props%emmxolw
      cldasymm = Cldrad_props%cldasymm
      ! single scattering and extinction are weighted with the
      ! cloud deltaz (this weighting is omitted in the shortwave code)
      cldsct = Cldrad_props%cldsct
      cldext = Cldrad_props%cldext
  !BW do np = 1, size(cldsct,5)
  !BW do nb = 1, size(cldsct,4)
  !BW    cldsct(:,:,:,nb,np)   = 1.0e-3*Cldrad_props%cldsct(:,:,:,nb,np)*clouddeltaz
  !BW    cldext(:,:,:,nb,np)   = 1.0e-3*Cldrad_props%cldext(:,:,:,nb,np)*clouddeltaz
  !BW enddo
  !BW enddo

!--------------------------------------------------------------------
!    call cloudrad_package_dealloc to deallocate the local derived type
!    variable arrays.
!--------------------------------------------------------------------
     !call microrad_properties_dealloc (Microrad_props)
     !call cldrad_props_dealloc (Cldrad_props)
      do n = 1, size(Microrad_props(:))
        call Microrad_props(n)%dealloc
      enddo
      call Cldrad_props%dealloc

!---------------------------------------------------------------------


end subroutine cloud_radiative_properties


!####################################################################
subroutine cloudrad_package_end (Cldrad_control)

type(cloudrad_control_type), intent(in) :: Cldrad_control

!--------------------------------------------------------------------
!    cloudrad_package_end is the destructor for cloudrad_package_mod.
!--------------------------------------------------------------------

!-------------------------------------------------------------------
!    deactivate the modules which are component modules of
!    cloudrad_package_mod.
!-------------------------------------------------------------------
      if (.not. Cldrad_control%do_no_clouds) then
        call cloudrad_diagnostics_end (Cldrad_control)
      endif
      if (Cldrad_control%do_pred_cld_microphys) then
        call microphys_rad_end
      endif
!---------------------------------------------------------------------
!    mark the module as not initialized.
!---------------------------------------------------------------------
      module_is_initialized = .false.

!---------------------------------------------------------------------


end subroutine cloudrad_package_end


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!                    PRIVATE SUBROUTINES
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!#####################################################################
subroutine combine_cloud_properties (is, js, Rad_time, Time_next,  &
                                     deltaz,  Cldrad_control, Cld_spec, &
                                     Cloud_microphys, Microrad_props, &
                                     Cldrad_props)

!----------------------------------------------------------------------
!    combine_cloud_properties produces cloud-radiative properties fields
!    for the total-cloud field in each grid box, using as input the
!    properties and characteristics of the various cloud types that may
!    be present (large-scale, uw shallow).
!----------------------------------------------------------------------

integer,                        intent(in)    :: is, js
type(time_type),                intent(in)    :: Rad_time, Time_next
real, dimension(:,:,:),         intent(in)    :: deltaz
type(cloudrad_control_type),    intent(in)    :: Cldrad_control
type(cld_specification_type),   intent(inout) :: Cld_spec
type(microphysics_type),        intent(in)    :: Cloud_microphys(:)
type(microrad_properties_type), intent(in)    :: Microrad_props(:)
type(cldrad_properties_type),   intent(inout) :: Cldrad_props

!----------------------------------------------------------------------
!   intent(in) variables:
!
!       Cloud_microphys  microphysical specification for
!                        all cloud schemes [ microphysics_type ]
!       Microrad_props   cloud radiative properties for all
!                        cloud schemes [ microrad_properties_type ]
!
!    intent(inout) variables:
!
!      Cldrad_props    cloud radiative properties on model grid,
!                      [ cldrad_properties_type ]
!
!---------------------------------------------------------------------

!---------------------------------------------------------------------
!    call comb_cldprops_calc with the appropriate arguments dependent
!    upon the cloud / convection scheme which is active.
!---------------------------------------------------------------------
!     if (Cldrad_control%do_strat_clouds .and.    &
!         Cldrad_control%do_uw_clouds .and. &
!         Cldrad_control%do_donner_deep_clouds) then

!----------------------------------------------------------------------
!    if strat_cloud, donner_deep and uw shallow are all active, then
!    lw and sw cloud radiative properties are microphysically based,
!    and large-scale, donner mesoscale and cell-scale, and uw shallow
!    cloud properties are available. call comb_cldprops_calc to combine
!    these into a single set of cloud radiative properties to be used
!    by the radiation package.
!---------------------------------------------------------------------

      if (size(Microrad_props(:)) > 1) then
        call comb_cldprops_calc (deltaz, Cld_spec%stoch_cloud_type, &
                                 Cldrad_props%cldext,   &
                                 Cldrad_props%cldsct,   &
                                 Cldrad_props%cldasymm,  &
                                 Cldrad_props%abscoeff,   &
                                 Cldrad_control, Cloud_microphys, Microrad_props)

!---------------------------------------------------------------------
!    if one cloud scheme alone is active, then total cloud values are
!    defined as just those values
!----------------------------------------------------------------------
      else if (size(Microrad_props(:)) == 1) then
!BW     if (Cldrad_control%do_sw_micro) then
          Cldrad_props%cldsct(:,:,:,:,1) = Microrad_props(1)%cldsct(:,:,:,:)
          Cldrad_props%cldext(:,:,:,:,1) = Microrad_props(1)%cldext(:,:,:,:)
          Cldrad_props%cldasymm(:,:,:,:,1) = Microrad_props(1)%cldasymm(:,:,:,:)
!BW     endif
        if (Cldrad_control%do_lw_micro) then
          Cldrad_props%abscoeff(:,:,:,:,1) = Microrad_props(1)%abscoeff(:,:,:,:)
        endif


      endif

!--------------------------------------------------------------------


end subroutine  combine_cloud_properties


!###################################################################

                   end module cloudrad_package_mod
