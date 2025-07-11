                 module strat_clouds_W_mod
! <CONTACT:  Fei.Liu@noaa.gov - fil
!    strat_clouds_W_mod obtains the cloud specification variables
!    for the klein strat cloud parameterization from cloud_rad_mod
!    and makes them available to the radiation package.

!   shared modules:

use constants_mod,          only: radian
use time_manager_mod,       only: time_type, time_manager_init
use mpp_mod,                only: input_nml_file
use fms_mod,                only: mpp_pe, &
                                  mpp_root_pe, stdlog,  fms_init, &
                                  write_version_number, &
                                  check_nml_error, error_mesg,   &
                                  FATAL

!   atmos shared modules:

use aerosol_types_mod,      only: aerosol_type
use physics_radiation_exch_mod, only : exchange_control_type

!   atmos param module

use cloud_rad_mod,       only: cloud_rad_init, cloud_summary3, &
                               snow_and_rain

!   cloud radiation shared module

use cloudrad_types_mod,  only: cld_specification_type, &
                               cldrad_properties_type,  &
                               microphysics_type, &
                               cloudrad_control_type

!    stochastic cloud generator modules

use random_number_streams_mod, only: get_random_number_streams
use random_numbers_mod,        only: randomNumberStream
use cloud_generator_mod,       only: cloud_generator_init, &
                                     generate_stochastic_clouds,&
                                     cloud_generator_end
!--------------------------------------------------------------------

implicit none
private

!--------------------------------------------------------------------
!    strat_clouds_W_mod obtains the cloud specification variables
!    for the klein strat cloud parameterization from cloud_rad_mod
!    and makes them available to the radiation package.
!--------------------------------------------------------------------


!---------------------------------------------------------------------
!----------- version number for this module --------------------------

character(len=128)  :: version =  '$Id$'
character(len=128)  :: tagname =  '$Name$'


!---------------------------------------------------------------------
!-------  interfaces --------

public          &
          strat_clouds_W_init, strat_clouds_amt, &
          strat_clouds_W_end

!---------------------------------------------------------------------
!-------- namelist  ---------

integer   :: seedperm = 0
logical   :: one_generator_call = .false.


namelist /strat_clouds_W_nml /                      &
                                seedperm, &
                                one_generator_call


!----------------------------------------------------------------------
!----  public data -------


!----------------------------------------------------------------------
!----  private data -------


logical                               :: module_is_initialized = .false.  ! module is initialized ?
integer          :: num_sw_bands, num_lw_bands
!----------------------------------------------------------------------
!----------------------------------------------------------------------



                              contains



!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!                     PUBLIC SUBROUTINES
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!#####################################################################
subroutine strat_clouds_W_init(latb, lonb, Cldrad_control, Exch_ctrl)

  real, dimension(:,:),        intent(in) :: latb, lonb
  type(cloudrad_control_type), intent(in) :: Cldrad_control
  type(exchange_control_type), intent(in), optional :: Exch_ctrl

!---------------------------------------------------------------------
!    strat_clouds_W_init is the constructor for strat_clouds_W_mod.
!---------------------------------------------------------------------
!       lonb      2d array of model longitudes on cell corners [ radians ]
!       latb      2d array of model latitudes at cell corners [radians]


!----------------------------------------------------------------------
!   local variables:

      integer   ::   ierr, io, logunit

!--------------------------------------------------------------------
!   local variables:
!
!      ierr     error code
!      io       error status returned from io operation
!
!--------------------------------------------------------------------

!---------------------------------------------------------------------
!    if routine has already been executed, exit.
!---------------------------------------------------------------------
      if (module_is_initialized) return

!---------------------------------------------------------------------
!    Save copy of number of sw and lw bands
!---------------------------------------------------------------------
      num_sw_bands = Cldrad_control%num_sw_cloud_bands
      num_lw_bands = Cldrad_control%num_lw_cloud_bands
!---------------------------------------------------------------------
!    verify that modules used by this module that are not called later
!    have already been initialized.
!---------------------------------------------------------------------
      call fms_init
      call time_manager_init
!BW   call radiation_clouds_util_init
      if (present(Exch_ctrl)) then
        call cloud_rad_init (Exch_ctrl)
      else
        call cloud_rad_init
      endif

!---------------------------------------------------------------------
!    read namelist.
!---------------------------------------------------------------------
   read (input_nml_file, nml=strat_clouds_W_nml, iostat=io)
   ierr = check_nml_error(io,'strat_clouds_W_nml')
!----------------------------------------------------------------------
!    write version number and namelist to logfile.
!---------------------------------------------------------------------
      call write_version_number (version, tagname)
      logunit = stdlog()
      if (mpp_pe() == mpp_root_pe() ) &
                 write (logunit, nml=strat_clouds_W_nml)

!--------------------------------------------------------------------
!    if stochastic clouds is active, be sure that the
!    cloud_generator module has been initialized.
!--------------------------------------------------------------------
      if (Cldrad_control%do_stochastic_clouds) then
          call cloud_generator_init
      endif

!--------------------------------------------------------------------
!    mark the module as initialized.
!---------------------------------------------------------------------
      module_is_initialized = .true.

!----------------------------------------------------------------------

end subroutine strat_clouds_W_init


!######################################################################
subroutine strat_clouds_amt (is, ie, js, je, Rad_time, pflux, &
                             press, temp, qv, land, &
                             Cldrad_control, Cld_spec, Lsc_microphys, Aerosol)

!---------------------------------------------------------------------
!    strat_clouds_amt defines the location, amount (cloud fraction),
!    and number of clouds present on the model grid, in addition to
!    liquid and ice-water paths, cloud thickness, and effective drop
!    and crystal sizes. if a microphysically-based cloud parameter-
!    ization is being used, particle sizes and concentrations are also
!    provided.
!----------------------------------------------------------------------

integer,                      intent(in)        :: is, ie, js, je
type(time_type),              intent(in)        :: Rad_time
real,    dimension(:,:,:),    intent(in)        :: pflux, press, temp, qv
real,    dimension(:,:),      intent(in)        :: land
type(cloudrad_control_type),  intent(in)        :: Cldrad_control
type(cld_specification_type), intent(inout)     :: Cld_spec
type(microphysics_type),      intent(inout)     :: Lsc_microphys
type(aerosol_type),           intent(in)        :: Aerosol

!----------------------------------------------------------------------
!   intent(in) variables:
!
!      is,ie,js,je  starting/ending subdomain i,j indices of data in
!                   the physics_window being integrated
!      Rad_time     time type variable containing radiation time
!      pflux        average of pressure at adjacent model levels
!                   [ (kg /( m s^2) ]
!      press        pressure at model levels (1:nlev), surface
!                   pressure is stored at index value nlev+1
!                   [ (kg /( m s^2) ]
!      temp         temperature at model levels (1:nlev), to be used
!                   in cloud calculations
!                   [ deg K ]
!      qv           water vapor specific humidity at model levels
!                   (1:nlev), to be used in cloud calculations
!      land         fraction of grid box covered by land
!                   [ non-dimensional ]
!
!   intent(inout), optional variables:
!
!      Cld_spec     cld_specification_type variable containing the
!                   cloud specification input fields needed by the
!                   radiation package
!
!               the following elements of Cld_spec are defined here:
!
!                  %cmxolw  fraction of maximally overlapped clouds
!                           seen by the longwave radiation
!                           [ dimensionless ]
!                  %crndlw  fraction of randomly overlapped clouds
!                           seen by the longwave radiation
!                           [ dimensionless ]
!                  %camtsw  cloud fraction seen by the shortwave
!                           radiation; the sum of the maximally
!                           overlapped and randomly overlapped
!                           longwave cloud fractions  [ dimensionless ]
!                  %nmxolw  number of maximally overlapped longwave
!                           clouds in each grid column.
!                  %nrndlw  number of randomly overlapped longwave
!                           clouds in each grid column.
!                  %ncldsw  number of clouds seen by the shortwave
!                           radiation in each grid column.
!                  %cloud_thickness
!                           number of model layers over which the cloud
!                           in this grid box extends
!                  %lwp     liquid water path
!                           [ kg / m^2 ]
!                  %iwp     ice water path
!                           [ kg / m^2 ]
!                  %reff_liq
!                           effective drop radius [ microns ]
!                  %reff_ice
!                           effective ice particle size [ microns ]
!
!      Lsc_microphys
!                   microphysics_type variable containing the size,
!                   concentration and fraction of the four condensate
!                   types (cloud drop, cloud ice, rain, snow) in the
!                   grid box, present when microphysically-based
!                   cloud radiation properties are desired.
!
!               the following components of this variable are output
!               from this routine when microphysically-based properties
!               are desired:
!
!                  %conc_ice  ice particle concentration [ g /m^3 ]
!                  %conc_drop cloud droplet concentration [ g /m^3 ]
!                  %size_ice  ice particle effective diameter
!                  [ microns ]
!                  %size_drop cloud droplet effective diameter
!                  [ microns ]
!
!----------------------------------------------------------------------

!-------------------------------------------------------------------
!    local variables

      real, dimension (size(pflux,1), size(pflux,2),  &
                       size(pflux,3)-1) ::      cldamt

      real, dimension (size(pflux,1), size(pflux,2),  &
                       size(pflux,3)-1, num_lw_bands) :: &
                         ql_stoch_lw2, qi_stoch_lw2, qa_stoch_lw2, &
                         qn_stoch_lw2, qni_stoch_lw2

      real, dimension (size(pflux,1), size(pflux,2),  &
                       size(pflux,3)-1, num_sw_bands) :: &
                         ql_stoch_sw2, qi_stoch_sw2, qa_stoch_sw2, &
                         qn_stoch_sw2, qni_stoch_sw2

      real, dimension (size(pflux,1), size(pflux,2),                 &
                       size(pflux,3)-1,                              &
                       num_lw_bands + num_sw_bands), &
             target :: ql_stoch, qi_stoch, qa_stoch, qn_stoch, qni_stoch

      real, dimension(:, :, :, :), pointer :: &
                  ql_stoch_lw, qi_stoch_lw, qa_stoch_lw, qn_stoch_lw, &
                  ql_stoch_sw, qi_stoch_sw, qa_stoch_sw, qn_stoch_sw, &
                  qni_stoch_lw, qni_stoch_sw

!      integer, dimension(size(Cld_spec%cld_thickness_lw_band, 1), &
!                         size(Cld_spec%cld_thickness_lw_band, 2), &
!                         size(Cld_spec%cld_thickness_lw_band, 3), &
      integer, dimension(size(temp, 1), size(temp, 2), size(temp, 3), &
          num_lw_bands + num_sw_bands) ::         &
                        cld_thickness

      integer, dimension (size(pflux,1), size(pflux,2), &
                          size(pflux,3)-1) ::   ktop, kbtm

      integer, dimension (size(pflux,1), size(pflux,2)) :: &
                                                 ncldlvls

      type(randomNumberStream), &
                dimension(size(pflux,1), size(pflux,2)) :: streams
      integer     ::    kx
      integer     ::    i, j, k, kc, nb
      real        ::   seedwts(8) = (/3000.,1000.,300.,100.,30.,10.,3.,1./)

!-------------------------------------------------------------------
!    local variables:
!
!       cldamt          cloud fraction, in cloud-space when microphysics
!                       not being used, in model-space when microphysics
!                       is active
!                       [ dimensionless ]
!       lwp             cloud liquid water path in cloud-space
!                       [ kg condensate / m^2 ]
!       iwp             cloud ice path, in cloud-space
!                       [ kg condensate / m^2 ]
!       reff_liq        effective radius for liquid clouds,
!                       in cloud-space [ microns ]
!       reff_ice        effective particle size for ice clouds
!                       in cloud-space [ microns ]
!       ktop            index of the model level which is cloud top
!       kbtm            index of the model level which is cloud base
!       ncldlvls        number of layers with cloud in a column
!       kx              number of model layers
!       i,j,k,kc        do-loop indices
!
!--------------------------------------------------------------------

!---------------------------------------------------------------------
!    be sure module has been initialized.
!---------------------------------------------------------------------
      if (.not. module_is_initialized ) then
        call error_mesg ('strat_clouds_W_mod',   &
             'module has not been initialized', FATAL )
      endif

!---------------------------------------------------------------------
!    define the number of model layers.
!---------------------------------------------------------------------
      kx = size (press,3)

!----------------------------------------------------------------------
!    compute the cloud specification properties under the assumption
!    of random cloud overlap.
!----------------------------------------------------------------------
      if (Cldrad_control%do_random_overlap) then

!----------------------------------------------------------------------
!    if microphysically-based radiative properties are needed, call
!    cloud_summary3 with the Lsc_microphys% optional arguments.
!----------------------------------------------------------------------
        if (Cldrad_control%do_pred_cld_microphys) then

!--------------------------------------------------------------------
!    call cloud_summary3 with the full cloud field, regardless of
!    whether or not stochastic clouds are active.
!    the full cloud field is assumed random overlap when stochastic
!    clouds are activated.
!--------------------------------------------------------------------
            where (Cld_spec%cloud_area(:,:,:) > 0.0)
              Cld_spec%cld_thickness(:,:,:) = 1
            end where
          call cloud_summary3 (is, js, land,   &
                               Cldrad_control%using_fu2007, &
                               Cld_spec%cloud_water, &
                               Cld_spec%cloud_ice, Cld_spec%cloud_area,&
                               Cld_spec%cloud_droplet, &
                               Cld_spec%cloud_ice_num, &
                               press(:,:,1:kx), pflux, temp, Aerosol, ncldlvls, &
                               cldamt, Cld_spec%lwp, Cld_spec%iwp,   &
                               Cld_spec%reff_liq, Cld_spec%reff_ice, &
                               conc_drop= Lsc_microphys%conc_drop, &
                               conc_ice = Lsc_microphys%conc_ice, &
                               size_drop =Lsc_microphys%size_drop,   &
                               size_ice = Lsc_microphys%size_ice, &
                         droplet_number = Lsc_microphys%droplet_number, &
                            ice_number = Lsc_microphys%ice_number )


          call snow_and_rain(Cld_spec%cloud_area, press(:,:,1:kx),  &
                             pflux, temp,  cldamt, Cld_spec%snow,  &
                             Cld_spec%rain, Cld_spec%snow_size,  &
                             Cld_spec%rain_size,   &
                             Lsc_microphys%conc_rain,   &
                             Lsc_microphys%conc_snow,   &
                             Lsc_microphys%size_rain,   &
                             Lsc_microphys%size_snow )


          cldamt = MIN (cldamt, 1.0)
          Cld_spec%ncldsw        = ncldlvls
          Cld_spec%nrndlw        = ncldlvls
          Cld_spec%camtsw        = cldamt
          Cld_spec%crndlw        = cldamt
          Lsc_microphys%cldamt   = cldamt

!---------------------------------------------------------------------
!    if using stochastic clouds for either sw or lw, Initialize the random number streams,
!       one per grid cell, with unique and replicable integer based
!       on grid location and model date/time
!---------------------------------------------------------------------
          if (Cldrad_control%do_stochastic_clouds) then
            call get_random_number_streams (is, js, Rad_time, temp(:,:,1), streams, perm=seedperm)

            if (one_generator_call) then
!---------------------------------------------------------------------
!    then generate all the subcolumns at once and divide them into
!    those needed  for the sw and lw bands.
!    call routine to obtain  band-dependent values of ql, qi and qa.
!---------------------------------------------------------------------
              call generate_stochastic_clouds (        &
                      streams,                &
                      Cld_spec%cloud_water,   &
                      Cld_spec%cloud_ice,     &
                      Cld_spec%cloud_area,    &
                      Cld_spec%cloud_droplet, &
                      Cld_spec%cloud_ice_num, &
                      pFull = press(:, :, :kx),&
                      pHalf = pflux, &
                      temperature = temp(:, :, :kx),        &
                      qv= qv(:, :, :kx), &
                      cld_thickness = cld_thickness, &
                      ql_stoch = ql_stoch, &
                      qi_stoch = qi_stoch, &
                      qa_stoch = qa_stoch, &
                      qn_stoch = qn_stoch, &
                      qni_stoch = qni_stoch )

          ql_stoch_lw => ql_stoch(:, :, :, 1:num_lw_bands)
          qi_stoch_lw => qi_stoch(:, :, :, 1:num_lw_bands)
          qa_stoch_lw => qa_stoch(:, :, :, 1:num_lw_bands)
          qn_stoch_lw => qn_stoch(:, :, :, 1:num_lw_bands)
          qni_stoch_lw => qni_stoch(:, :, :, 1:num_lw_bands)
          Cld_spec%cld_thickness_lw_band = &
                       cld_thickness(:, :, :, 1:num_lw_bands)

          ql_stoch_sw => ql_stoch(:, :, :, num_lw_bands+1:)
          qi_stoch_sw => qi_stoch(:, :, :, num_lw_bands+1:)
          qa_stoch_sw => qa_stoch(:, :, :, num_lw_bands+1:)
          qn_stoch_sw => qn_stoch(:, :, :, num_lw_bands+1:)
          qni_stoch_sw => qni_stoch(:, :, :, num_lw_bands+1:)
          Cld_spec%cld_thickness_sw_band = &
                     cld_thickness(:, :, :, num_lw_bands+1:)

!---------------------------------------------------------------------
!    call cloud_summary3 for each lw band, using the band-dependent
!    cloud inputs, to obtain band-dependent values of liquid and ice
!    size and concentration.
!---------------------------------------------------------------------
            do nb=1,num_lw_bands
              call cloud_summary3 (          &
                is, js, land, &
                Cldrad_control%using_fu2007, &
                ql_stoch_lw(:,:,:,nb),&
                qi_stoch_lw(:,:,:,nb), qa_stoch_lw(:,:,:,nb),&
                qn_stoch_lw(:,:,:,nb), &
                qni_stoch_lw(:,:,:,nb),&
                press(:,:,1:kx), pflux, temp, Aerosol, ncldlvls, &
                cldamt, Cld_spec%lwp_lw_band(:,:,:,nb),&
                Cld_spec%iwp_lw_band(:,:,:,nb),   &
                Cld_spec%reff_liq_lw_band(:,:,:,nb), &
                Cld_spec%reff_ice_lw_band(:,:,:,nb), &
                conc_drop= Lsc_microphys%lw_stoch_conc_drop(:,:,:,nb), &
                conc_ice = Lsc_microphys%lw_stoch_conc_ice(:,:,:,nb), &
                size_drop =Lsc_microphys%lw_stoch_size_drop(:,:,:,nb), &
                size_ice = Lsc_microphys%lw_stoch_size_ice(:,:,:,nb), &
       droplet_number = Lsc_microphys%lw_stoch_droplet_number(:,:,:,nb), &
          ice_number =  Lsc_microphys%lw_stoch_ice_number(:,:,:,nb) )

              !now that the vertical cloud fraction has been used to
              !properly calculate the in-cloud particle size, rescale
              !the concentrations and cloud amounts to that the cloud
              !amount is unity in any partially cloudy sub-column.
              !
              !This is necessary so that the radiation code will not
              !do cloud fraction weights of cloudy and clear sky fluxes.
              !
              !The rescaling of the concentrations is necessary so that the
              !total optical depth of the layer is constant.  Note that this
              !works because cloud extinction is linear in the concentration
              Lsc_microphys%lw_stoch_conc_drop(:,:,:,nb) = &
                   Lsc_microphys%lw_stoch_conc_drop(:,:,:,nb) * cldamt(:,:,:)
              Lsc_microphys%lw_stoch_conc_ice(:,:,:,nb) = &
                   Lsc_microphys%lw_stoch_conc_ice(:,:,:,nb) * cldamt(:,:,:)
              where (cldamt .gt. 0.) cldamt = 1.

              Cld_spec%nrndlw_band(:,:,nb) = ncldlvls(:,:)
              Lsc_microphys%lw_stoch_cldamt(:,:,:,nb) = cldamt
            end do

!---------------------------------------------------------------------
!    call cloud_summary3 for each sw band, using the band-dependent
!    cloud inputs, to obtain band-dependent values of liquid and ice
!    size and concentration.
!---------------------------------------------------------------------
          do nb=1,size(Lsc_microphys%sw_stoch_conc_ice,4)
            call cloud_summary3 (                            &
              is, js, land,  &
              Cldrad_control%using_fu2007, &
              ql_stoch_sw(:,:,:,nb), &
              qi_stoch_sw(:,:,:,nb), qa_stoch_sw(:,:,:,nb),&
              qn_stoch_sw(:,:,:,nb), &
              qni_stoch_sw(:,:,:,nb), &
              press(:,:,1:kx), pflux, temp, Aerosol, ncldlvls, &
              cldamt, Cld_spec%lwp_sw_band(:,:,:,nb), &
              Cld_spec%iwp_sw_band(:,:,:,nb),   &
              Cld_spec%reff_liq_sw_band(:,:,:,nb), &
              Cld_spec%reff_ice_sw_band(:,:,:,nb), &
              conc_drop= Lsc_microphys%sw_stoch_conc_drop(:,:,:,nb), &
              conc_ice = Lsc_microphys%sw_stoch_conc_ice(:,:,:,nb), &
              size_drop =Lsc_microphys%sw_stoch_size_drop(:,:,:,nb), &
              size_ice = Lsc_microphys%sw_stoch_size_ice(:,:,:,nb), &
       droplet_number = Lsc_microphys%sw_stoch_droplet_number(:,:,:,nb), &
         ice_number = Lsc_microphys%sw_stoch_ice_number(:,:,:,nb) )

         !now that the vertical cloud fraction has been used to
         !properly calculate the in-cloud particle size, rescale
         !the concentrations and cloud amounts to that the cloud
         !amount is unity in any partially cloudy sub-column.
         !
         !This is necessary so that the radiation code will not
         !do cloud fraction weights of cloudy and clear sky fluxes.
         !
         !The rescaling of the concentrations is necessary so that         the
         !total optical depth of the layer is constant.  Note that         this
         !works because cloud extinction is linear in the concentra        tion
             Lsc_microphys%sw_stoch_conc_drop(:,:,:,nb) = &
             Lsc_microphys%sw_stoch_conc_drop(:,:,:,nb) * cldamt(:,:,:)
             Lsc_microphys%sw_stoch_conc_ice(:,:,:,nb) = &
             Lsc_microphys%sw_stoch_conc_ice(:,:,:,nb) * cldamt(:,:,:)
             where (cldamt .gt. 0.) cldamt = 1.

             Cld_spec%ncldsw_band(:,:,nb) = ncldlvls(:,:)
             Lsc_microphys%sw_stoch_cldamt(:,:,:,nb) = cldamt
           end do

        else  ! (one_call)
!---------------------------------------------------------------------
!    call routine to obtain lw band-dependent values of ql, qi and qa.
!---------------------------------------------------------------------
            call generate_stochastic_clouds (        &
                     streams,                &
                     Cld_spec%cloud_water,   &
                     Cld_spec%cloud_ice,     &
                     Cld_spec%cloud_area,    &
                     Cld_spec%cloud_droplet, &
                     Cld_spec%cloud_ice_num, &
                     pFull    = press(:, :, :kx),     &
                     pHalf    = pflux,&
                     temperature = temp(:, :, :kx),   &
                     qv = qv(:,:, :kx),   &
                     cld_thickness = Cld_spec%cld_thickness_lw_band, &
                     ql_stoch = ql_stoch_lw2, &
                     qi_stoch = qi_stoch_lw2, &
                     qa_stoch = qa_stoch_lw2, &
                     qn_stoch = qn_stoch_lw2, &
                     qni_stoch = qni_stoch_lw2 )

!---------------------------------------------------------------------
!    call routine to obtain sw band-dependent values of ql, qi and qa.
!---------------------------------------------------------------------
            call generate_stochastic_clouds (        &
                     streams,                &
                     Cld_spec%cloud_water,   &
                     Cld_spec%cloud_ice,     &
                     Cld_spec%cloud_area,    &
                     Cld_spec%cloud_droplet, &
                     Cld_spec%cloud_ice_num, &
                     pFull    = press(:, :, :kx),     &
                     pHalf    = pflux,&
                     temperature = temp(:, :, :kx),   &
                     qv = qv(:,:, :kx),   &
                     cld_thickness = Cld_spec%cld_thickness_sw_band, &
                     ql_stoch = ql_stoch_sw2, &
                     qi_stoch = qi_stoch_sw2, &
                     qa_stoch = qa_stoch_sw2, &
                     qn_stoch = qn_stoch_sw2, &
                     qni_stoch = qni_stoch_sw2 )

!---------------------------------------------------------------------
!    call cloud_summary3 for each lw band, using the band-dependent
!    cloud inputs, to obtain band-dependent values of liquid and ice
!    size and concentration.
!---------------------------------------------------------------------
            do nb=1,num_lw_bands
              call cloud_summary3 (          &
                is, js, land,  &
                Cldrad_control%using_fu2007, &
                ql_stoch_lw2(:,:,:,nb),&
                qi_stoch_lw2(:,:,:,nb), qa_stoch_lw2(:,:,:,nb),&
                qn_stoch_lw2(:,:,:,nb), &
                qni_stoch_lw2(:,:,:,nb), &
                press(:,:,1:kx), pflux, temp, Aerosol, ncldlvls, &
                cldamt, Cld_spec%lwp_lw_band(:,:,:,nb),&
                Cld_spec%iwp_lw_band(:,:,:,nb),   &
                Cld_spec%reff_liq_lw_band(:,:,:,nb), &
                Cld_spec%reff_ice_lw_band(:,:,:,nb), &
                conc_drop= Lsc_microphys%lw_stoch_conc_drop(:,:,:,nb), &
                conc_ice = Lsc_microphys%lw_stoch_conc_ice(:,:,:,nb), &
                size_drop =Lsc_microphys%lw_stoch_size_drop(:,:,:,nb), &
                size_ice = Lsc_microphys%lw_stoch_size_ice(:,:,:,nb), &
       droplet_number = Lsc_microphys%lw_stoch_droplet_number(:,:,:,nb), &
         ice_number =   Lsc_microphys%lw_stoch_ice_number(:,:,:,nb) )

              !now that the vertical cloud fraction has been used to
              !properly calculate the in-cloud particle size, rescale
              !the concentrations and cloud amounts to that the cloud
              !amount is unity in any partially cloudy sub-column.
              !
              !This is necessary so that the radiation code will not
              !do cloud fraction weights of cloudy and clear sky fluxes.
              !
              !The rescaling of the concentrations is necessary so that the
              !total optical depth of the layer is constant.  Note that this
              !works because cloud extinction is linear in the concentration
              Lsc_microphys%lw_stoch_conc_drop(:,:,:,nb) = &
                   Lsc_microphys%lw_stoch_conc_drop(:,:,:,nb) * cldamt(:,:,:)
              Lsc_microphys%lw_stoch_conc_ice(:,:,:,nb) = &
                   Lsc_microphys%lw_stoch_conc_ice(:,:,:,nb) * cldamt(:,:,:)
              where (cldamt .gt. 0.) cldamt = 1.

              Cld_spec%nrndlw_band(:,:,nb) = ncldlvls(:,:)
              Lsc_microphys%lw_stoch_cldamt(:,:,:,nb) = cldamt
            end do

!---------------------------------------------------------------------
!    call cloud_summary3 for each sw band, using the band-dependent
!    cloud inputs, to obtain band-dependent values of liquid and ice
!    size and concentration.
!---------------------------------------------------------------------
            do nb=1,size(Lsc_microphys%sw_stoch_conc_ice,4)
              call cloud_summary3 (                            &
                is, js, land, &
                Cldrad_control%using_fu2007, &
                ql_stoch_sw2(:,:,:,nb), &
                qi_stoch_sw2(:,:,:,nb), qa_stoch_sw2(:,:,:,nb),&
                qn_stoch_sw2(:,:,:,nb), &
                qni_stoch_sw2(:,:,:,nb), &
                press(:,:,1:kx), pflux, temp, Aerosol, ncldlvls, &
                cldamt, Cld_spec%lwp_sw_band(:,:,:,nb), &
                Cld_spec%iwp_sw_band(:,:,:,nb),   &
                Cld_spec%reff_liq_sw_band(:,:,:,nb), &
                Cld_spec%reff_ice_sw_band(:,:,:,nb), &
                conc_drop= Lsc_microphys%sw_stoch_conc_drop(:,:,:,nb), &
                conc_ice = Lsc_microphys%sw_stoch_conc_ice(:,:,:,nb), &
                size_drop =Lsc_microphys%sw_stoch_size_drop(:,:,:,nb), &
                size_ice = Lsc_microphys%sw_stoch_size_ice(:,:,:,nb), &
       droplet_number = Lsc_microphys%sw_stoch_droplet_number(:,:,:,nb), &
         ice_number =  Lsc_microphys%sw_stoch_ice_number(:,:,:,nb))

              !now that the vertical cloud fraction has been used to
              !properly calculate the in-cloud particle size, rescale
              !the concentrations and cloud amounts to that the cloud
              !amount is unity in any partially cloudy sub-column.
              !
              !This is necessary so that the radiation code will not
              !do cloud fraction weights of cloudy and clear sky fluxes.
              !
              !The rescaling of the concentrations is necessary so that the
              !total optical depth of the layer is constant.  Note that this
              !works because cloud extinction is linear in the concentration
              Lsc_microphys%sw_stoch_conc_drop(:,:,:,nb) = &
                   Lsc_microphys%sw_stoch_conc_drop(:,:,:,nb) * cldamt(:,:,:)
              Lsc_microphys%sw_stoch_conc_ice(:,:,:,nb) = &
                   Lsc_microphys%sw_stoch_conc_ice(:,:,:,nb) * cldamt(:,:,:)
              where (cldamt .gt. 0.) cldamt = 1.

              Cld_spec%ncldsw_band(:,:,nb) = ncldlvls(:,:)
              Lsc_microphys%sw_stoch_cldamt(:,:,:,nb) = cldamt
            end do
         endif ! (one_generator_call)
       endif  ! (do_stochastic_clouds)

!---------------------------------------------------------------------
!    if microphysically-based radiative properties are not needed, call
!    cloud_summary3 without the Lsc_microphys% optional arguments.
!----------------------------------------------------------------------
        else  ! (not micro)

!--------------------------------------------------------------------
!    define the cloud thickness to be 1 at those points with cloud
!    present.
!---------------------------------------------------------------------
          where (Cld_spec%cloud_area(:,:,:) > 0.0)
            Cld_spec%cld_thickness(:,:,:) = 1
          end where
          call cloud_summary3 (is, js, land,  &
                               Cldrad_control%using_fu2007, &
                               Cld_spec%cloud_water, &
                               Cld_spec%cloud_ice, Cld_spec%cloud_area,&
                               Cld_spec%cloud_droplet, &
                               Cld_spec%cloud_ice_num, &
                               press(:,:,1:kx), pflux, temp, Aerosol, ncldlvls, &
                               cldamt, Cld_spec%lwp,   &
                               Cld_spec%iwp, Cld_spec%reff_liq,   &
                               Cld_spec%reff_ice)
          cldamt = MIN (cldamt, 1.0)
          Cld_spec%ncldsw        = ncldlvls
          Cld_spec%nrndlw        = ncldlvls
          Cld_spec%camtsw        = cldamt
          Cld_spec%crndlw        = cldamt
          Lsc_microphys%cldamt   = cldamt
        endif

!----------------------------------------------------------------------
!    all clouds are assumed to be randomly overlapped
!----------------------------------------------------------------------
        Cld_spec%nmxolw        = 0
        Cld_spec%cmxolw        = 0.0E+00

!---------------------------------------------------------------------
!    define cloud specification properties when max-random overlap is
!    assumed. in this case cloud in adjacent layers is assumed to be
!    part of the same cloud.
!---------------------------------------------------------------------
      else if (Cldrad_control%do_max_random_overlap) then

!----------------------------------------------------------------------
!    microphysically-based radiative properties are not implemented
!    with the max-random overlap assumption.
!----------------------------------------------------------------------
        if (Cldrad_control%do_pred_cld_microphys) then
          call error_mesg ('strat_clouds_W_mod', &
               'must use random overlap cloud assumption with strat '//&
              'clouds when microphysics are desired', FATAL)

!---------------------------------------------------------------------
!    if microphysically-based radiative properties are not needed, call
!    cloud_summary3 without the Lsc_microphys% optional arguments.
!----------------------------------------------------------------------
        else
          call cloud_summary3 (is, js, land,  &
                               Cldrad_control%using_fu2007, &
                               Cld_spec%cloud_water, &
                               Cld_spec%cloud_ice, Cld_spec%cloud_area,&
                               Cld_spec%cloud_droplet, &
                               Cld_spec%cloud_ice_num, &
                               press(:,:,1:kx), pflux, temp, Aerosol, ncldlvls, &
                               Cld_spec%camtsw, Cld_spec%lwp,   &
                               Cld_spec%iwp, Cld_spec%reff_liq,   &
                               Cld_spec%reff_ice,  &
                               ktop=ktop, kbot=kbtm)

!---------------------------------------------------------------------
!    when only bulk properties are returned, they are in cloud space,
!    and must be converted to physical space before being stored in
!    Cld_spec. random overlap and max overlap properties are assigned
!    according to the cloud thickness - multi layer clouds are assumed
!    to be max overlap.
!-------------------------------------------------------------------
          do j=1, size(press,2)
            do i=1, size(press,1)
              Cld_spec%ncldsw(i,j) = ncldlvls(i,j)
              do kc=1, Cld_spec%ncldsw(i,j)
                do k=ktop(i,j,kc), kbtm(i,j,kc)
                  if (ktop(i,j,kc) == kbtm(i,j,kc)) then
                    Cld_spec%crndlw(i,j,k) = Cld_spec%camtsw(i,j,k)
                    Cld_spec%cmxolw(i,j,k) = 0.0
                    Cld_spec%cld_thickness(i,j,k) = 1
                  else
                    Cld_spec%cmxolw(i,j,k) = Cld_spec%camtsw(i,j,k)
                    Cld_spec%crndlw(i,j,k) = 0.0
                    Cld_spec%cld_thickness(i,j,k) = kbtm(i,j,kc) -    &
                                                    ktop(i,j,kc) + 1
                  endif
                end do
                if (ktop(i,j,kc) == kbtm(i,j,kc)) then
                  Cld_spec%nrndlw(i,j) = Cld_spec%nrndlw(i,j) + 1
                else
                  Cld_spec%nmxolw(i,j) = Cld_spec%nmxolw(i,j) + 1
                endif
              end do
            end do
          end do
        endif ! (do_pred_micro)
      endif ! (do_random_overlap)

!---------------------------------------------------------------------


end subroutine strat_clouds_amt


!####################################################################
subroutine strat_clouds_W_end (Cldrad_control)

type(cloudrad_control_type), intent(in) :: Cldrad_control

!----------------------------------------------------------------------
!    strat_clouds_W_end is the destructor for strat_clouds_W_mod.
!----------------------------------------------------------------------

!---------------------------------------------------------------------
!    be sure module has been initialized.
!---------------------------------------------------------------------
      if (.not. module_is_initialized ) then
        call error_mesg ('strat_clouds_W_mod',   &
             'module has not been initialized', FATAL )
      endif

!---------------------------------------------------------------------
! close cloud_generator
!---------------------------------------------------------------------
      if (Cldrad_control%do_stochastic_clouds) then
          call cloud_generator_end()
      endif

!---------------------------------------------------------------------
!    mark the module as not initialized.
!---------------------------------------------------------------------
      module_is_initialized = .false.

!--------------------------------------------------------------------

end subroutine strat_clouds_W_end

!#################################################################

                    end module strat_clouds_W_mod
