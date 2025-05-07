                 module cloud_rad_mod
! Contact - Steve Klein
!        The cloud radiation module uses the stored values of the
!     prognostic cloud variables, and computes the cloud albedo and
!     absorption for the two shortwave bands (ultra-violet/visible and
!     near-infrared), the longwave cloud emissivity, and the 
!     fractional areas covered by clouds.
!
!      The cloud radiation module condenses the cloud information 
!     provided by the stratiform cloud scheme and converts it into
!     the areas covered by, the water paths and the effective particle 
!     sizes of liquid and ice. 
!
!  ************************Note***********************
!   This part of the documentation needs to be updated
!  ***************************************************
!
!  Diagnostic fields may be output to a netcdf file by specifying the
!  module name cloud_rad and the desired field names (given below)
!  in file diag_table. See the documentation for diag_manager.
!  
!  Diagnostic fields for module name: cloud_rad
!  
!     nisccp       frequency of sunlit times at the times of the radiation
!                  calculation at each point {fraction} [real,dimension(:,:)]
!  
!     pc#tau%      where # is a number from 1 to 7
!                  and   % is a number from 1 to 7
!                  {fraction} [real,dimension(:,:)]
!  
!                  Thus there are 49 diagnostic fields of this type.  All
!                  of them are necessary to receive the complete decomposition
!                  of clouds visible from space into the ISCCP categories.
!  
!                  The 7 cloud top pressure ("pc") categories and 7 optical
!                  depth ("tau") categories are defined as:
!  
!                  pc #      pc range (mb)    tau %        tau range
!                  ----    ----------------   -----    ---------------------
!  
!                   1              pc < 180     1     0.0    < tau < taumin 
!                   2        180 < pc < 310     2     taumin < tau < 1.3
!                   3        310 < pc < 440     3     1.3    < tau < 3.6
!                   4        440 < pc < 560     4     3.6    < tau < 9.4
!                   5        560 < pc < 680     5     9.4    < tau < 23
!                   6        680 < pc < 800     6     23     < tau < 60
!                   7        800 < pc                 60     < tau
!  
!                  What is saved in these diagnostics is the time mean of
!                  the area covered by clouds of this type when the sun is
!                  above the horizon. This is done so that the calculation 
!                  will mimic the ISCCP product which is broken down into 
!                  these categories only for sunlit places.
!  
!                  NOTE:  TO DETERMINE THE MEAN AREA COVERED BY A CLOUD TYPE 
!                         WHEN THE SUN IS ABOVE THE HORIZON YOU MUST DIVIDE
!                         BY NISCCP:
!  
!                         area of cloud type pc#tau% =   pc#tau% / nisccp
!  
!     aice         fractional area of sunlit clouds seen from space whose cloud 
!                  top contains ice. {fraction} [real,dimension(:,:)]
!  
!     reffice      time mean ice effective radius of cloud tops visible from
!                  space including areas where there is no such cloud {microns} 
!                  [real,dimension(:,:)]
!  
!                  NOTE:  THUS THE TIME MEAN CLOUD TOP EFFECTIVE RADIUS OF CLOUD 
!                         TOPS WITH ICE VISIBLE FROM SPACE IS:
!  
!                         mean reffice  =    reffice /  aice
!        
!     aliq         fractional area of sunlit clouds seen from space whose cloud 
!                  top contains liquid. {fraction} [real,dimension(:,:)]
!  
!     reffliq      time mean cloud droplet effective radius of cloud tops 
!                  visible from space including areas where there is no such 
!                  cloud {microns} [real,dimension(:,:)]
!     
!                  NOTE:  mean reffliq  =    reffliq / aliq
!  
!     alow         fractional area of sunlit clouds seen from space whose cloud 
!                  tops are low (pc > 680 mb). {fraction} [real,dimension(:,:)]
!  
!     tauicelow    time mean optical depth of ice for cloud tops visible from 
!                  space including areas where there is no such cloud {microns} 
!                  [real,dimension(:,:)]
!    
!     tauliqlow    time mean optical depth of liquid for cloud tops visible from 
!                  space including areas where there is no such cloud {microns} 
!                  [real,dimension(:,:)]
!     
!     tlaylow      time mean of the low level mean temperature (pc > 680 mb) 
!                  when low cloud tops are visible from space including times 
!                  where there is no such cloud {microns} 
!                  [real,dimension(:,:)]
!        
!     tcldlow      time mean of the cloud top temperature for cloud tops visible 
!                  from space including times where there is no such cloud 
!                  {microns}  [real,dimension(:,:)]
!        
!                  NOTE:  mean tauicelow  =    tauicelow / alow
!                         mean tauliqlow  =    tauliqlow / alow
!                         mean tlaylow    =    tlaylow   / alow
!                         mean tcldlow    =    tcldlow   / alow

! REFERENCES
! The shortwave properties of liquid clouds come from:
! 
 !     Slingo, A., 1989: A GCM parameterization for the shortwave 
 !     radiative properties of water clouds. J. Atmos. Sci., vol. 46, 
 !     pp. 1419-1427.
!
! The shortwave and longwave properties of ice clouds come from:
! 
 !     Ebert, E. E. and J. A. Curry, 1992: A parameterization of ice cloud
 !     optical properties for climate models. J. Geophys. Res., vol. 97,
 !     D1, pp. 3831-3836.
!
! The longwave emissivity parameterization of liquid clouds comes from:
! 
 !     Stephens, G. L., 1978: Radiation profiles in extended water clouds.
 !     II: Parameterization schemes. J. Atmos. Sci., vol. 35, 
 !     pp. 2123-2132.
!
! The parameterization of liquid cloud effective radius comes from:
! 
 !     Martin, G. M., D. W. Johnson, and A. Spice, 1994: The measurement 
 !     and parameterization of effective radius of droplets in warm stratocumulus
 !     clouds. J. Atmos. Sci, vol 51, pp. 1823-1842.
!
! The parameterization of ice cloud effective radius comes from:
! 
 !     Donner, L. J., C. J. Seman, B. J. Soden, R. S. Hemler, J. C. Warren,
 !     J. Strom, and K.-N. Liou, 1997: Large-scale ice clouds in the GFDL
 !     SKYHI general circulation model. J. Geophys. Res., vol. 102, D18,
 !     pp. 21,745-21,768.
!
! The algorithm to reproduce the ISCCP satellite view of clouds comes from:
! 
 !     Klein, S. A., and C. Jakob, 1999: Validation and sensitivities of 
 !     frontal clouds simulated by the ECMWF model. Monthly Weather Review,
 !     127(10),  2514-2531.
!

! FUTURE
! The optical depth and particle size for every model level will
!     become a diagnostic output field.

!   shared modules:

use mpp_mod,                    only: input_nml_file
use fms_mod,                    only: fms_init,       &
                                      stdlog, mpp_pe, mpp_root_pe, &
                                      write_version_number,  &
                                      error_mesg, FATAL,  &
                                      check_nml_error
use constants_mod,              only: RDGAS, GRAV, TFREEZE, DENS_H2O, &
                                      constants_init, pi
use gamma_mg_mod,               ONLY: gamma_mg, gamma_mg_init, gamma_mg_end
use lscloud_constants_mod,      ONLY: lscloud_constants_init, rhow, &
                                      di_mg, ci_mg
use time_manager_mod,           only: time_type, time_manager_init
use aerosol_types_mod,          only: aerosol_type
use physics_radiation_exch_mod, only: exchange_control_type

implicit none
private

!---------------------------------------------------------------------
!    cloud_rad_mod does the following:           
!     
!    (a)  subroutine cloud_summary3 returns cloud specification var-
!         iables that are used in calculating the cloud radiative prop-
!         erties. these include cloud locations, water paths and effect-
!         ive particle sizes for use in determining bulk properties and
!         concentrations and drop sizes if microphysically-based prop-
!         erties are desired.
!---------------------------------------------------------------------

!---------------------------------------------------------------------
!------------ version number for this module -------------------------
        
character(len=128) :: version = '$Id$'
character(len=128) :: tagname = '$Name$'


!---------------------------------------------------------------------- 
!--------- interfaces --------

public     &
         cloud_rad_init, cloud_rad_end, cloud_summary3, snow_and_rain
!---------------------------------------------------------------------
!    public subroutines:
!
!      cloud_rad_init
!                        Initializes values of qmin, N_land, and 
!                        N_ocean using values from strat_cloud namelist
!                        as well as reads its own namelist variables. 
!                        In addition, it registed diagnostic fields
!                        if needed, and returns the value of the
!                        cloud overlap to strat_cloud.

!----------------------------------------------------------------------

private     &
         max_rnd_overlap, rnd_overlap

!---------------------------------------------------------------------
!    private subroutines:
!
!       max_rnd_overlap
!       rnd_overlap
!---------------------------------------------------------------------

!---------------------------------------------------------------------
!-------- namelist  ---------

real         :: taucrit = 1.
logical      :: adjust_top = .true.
real         :: scale_factor = 1
real         :: qamin = 1.E-2
logical      :: do_brenguier = .false.
logical      :: snow_in_cloudrad = .false.
logical      :: rain_in_cloudrad = .false.

!  taucrit       critical optical depth for switching direct beam 
!                     to diffuse beam for use in Delta-Eddington 
!                     solution [ dimensionless] 
!  adjust_top    logical variable indicating whether or not to use 
!                     the code which places the top and bottom of the 
!                     cloud at the faces which are most in view from
!                     the top and bottom of the cloud block. this is 
!                     done to avoid undue influence of very small cloud
!                     fractions. if true this adjustment of tops is 
!                     performed; if false this is not performed.
!  scale_factor  factor which multiplies actual cloud optical 
!                     depths to account for the plane-parallel homo-
!                     genous cloud bias  (e.g. Cahalan effect).
!                     [ dimensionless] 
!  qamin         minimum permissible cloud fraction 
!                     [ dimensionless] 
!  do_brenguier  should drops at top of stratocumulus clouds be
!                     scaled?

namelist /cloud_rad_nml/                                       &
                         taucrit, adjust_top, scale_factor,    &
                         qamin, do_brenguier, snow_in_cloudrad, &
                         rain_in_cloudrad

!------------------------------------------------------------------
!---- public data ------


!-------------------------------------------------------------------
!---- private data ------

!----------------------------------------------------------------------
!   module variables obtained form other modules during initialization
!----------------------------------------------------------------------
integer :: overlap      !integer variable indicating which overlap 
                        !assumption to use:
                        !overlap = 1. means condensate in adjacent levels 
                        !          is treated as part of the same cloud
                        !          i.e. maximum-random overlap
                        !overlap = 2. means condensate in adjacent levels 
                        !          is treated as different clouds
                        !          i.e. random overlap

real    :: N_land       ! number of cloud droplets in liquid
                        ! clouds over land  [ m**(-3) ]
real    :: N_ocean      ! number of cloud droplets in liquid
                        ! clouds over ocean [ m**(-3) ]
real    :: N_min
real    :: qmin         ! minimum permissible cloud 
                        ! condensate [ kg condensate / kg air ]       
real    :: qcvar
real    :: dcs
real    :: min_diam_ice
real    :: min_diam_drop
real    :: max_diam_drop
logical :: do_liq_num   ! use prog. droplet number ?
logical :: do_ice_num   ! use prog ice crystal number?

!-------------------------------------------------------------------
!   needed physical parameters:
!-------------------------------------------------------------------
real, parameter :: taumin = 1.E-06  ! minimum permissible tau  
                                    ! [ dimensionless ]
real, parameter :: k_land = 1.143   ! ratio of effective radius to 
                                    ! volume radius for continental
                                    ! air masses  [ dimensionless ]
real, parameter :: k_ocean = 1.077  ! ratio of effective radius to 
                                    ! volume radius for continental
                                    ! air masses  [ dimensionless ]
!----------------------------------------------------------------------
!    diagnostics variables.        
!----------------------------------------------------------------------
character(len=8)    :: mod_name = 'cloud_rad'
real                :: missing_value = -999.

!----------------------------------------------------------------------
!    misc. variables.        
!----------------------------------------------------------------------
logical   :: module_is_initialized = .false.  ! is module initialized ?

!---------------------------------------------------------------------




                         contains



!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!                     PUBLIC SUBROUTINES
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!#######################################################################
subroutine cloud_rad_init (Exch_ctrl)
  
!--------------------------------------------------------------------
!    cloud_rad_init is the constructor for cloud_rad_mod.
!  
!   Called once to initialize cloud_rad module.   This routine reads the
!   namelist and obtains module variables this module needs that are 
!   defined in other modules. 
!
!   Fatal crashes occur in initialization of the module if:
!   1. overlap does not equal 1 or 2
!   2. taucrit < 0.
!   3. scale_factor < 0.
!   4. qamin outside of the range of 0 to 1.
!--------------------------------------------------------------------

type(exchange_control_type), intent(in), optional :: Exch_ctrl     

!----------------------------------------------------------------------
!  Internal variables

      integer  :: io, ierr, logunit

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
      call constants_init
      call lscloud_constants_init
      call gamma_mg_init

!--------------------------------------------------------------------
!    read namelist.
!--------------------------------------------------------------------
      read (input_nml_file, nml=cloud_rad_nml, iostat=io)
      ierr = check_nml_error(io,'cloud_rad_nml')

!---------------------------------------------------------------------
!    write version number and namelist to logfile.
!---------------------------------------------------------------------
      call write_version_number (version, tagname)
      logunit = stdlog()
      if (mpp_pe() == mpp_root_pe() ) &
                        write (logunit, nml=cloud_rad_nml)

!-----------------------------------------------------------------------
!    prevent unreasonable values.
!-----------------------------------------------------------------------
      if (taucrit .lt. 0.) &
                call error_mesg  ('cloud_rad_init in cloud_rad module',&
                  'taucrit must be greater than or equal to 0. ', FATAL)
      if (scale_factor .lt. 0.) &
                call error_mesg  ('cloud_rad_init in cloud_rad module',&
                         'scale_factor must be greater than 0. ', FATAL)
      if (qamin .le. 0. .or. qamin .ge. 1.) &
                call error_mesg  ('cloud_rad_init in cloud_rad module',&
                               'qamin must be between 0. and 1.', FATAL)
        
!-----------------------------------------------------------------------
!    define needed module variables obtained from other modules.
!-----------------------------------------------------------------------
      if (present (Exch_ctrl)) then
        qmin = Exch_ctrl%qmin
        N_min = Exch_ctrl%N_min
        N_land = Exch_ctrl%N_land
        N_ocean = Exch_ctrl%N_ocean
        qcvar = Exch_ctrl%qcvar
        overlap = Exch_ctrl%overlap
        do_liq_num = Exch_ctrl%do_liq_num
        do_ice_num = Exch_ctrl%do_ice_num 
        min_diam_ice = Exch_ctrl%min_diam_ice
        min_diam_drop = Exch_ctrl%min_diam_drop
        max_diam_drop = Exch_ctrl%max_diam_drop
        dcs = Exch_ctrl%dcs
      else
        call error_mesg ('cloud_rad_init', ' first call to &
                   &cloud_rad_init does not have Exch_ctrl as arg',&
                                                                    FATAL)
      endif
 
!-----------------------------------------------------------------------
!    mark this module as initialized.
!-----------------------------------------------------------------------
      module_is_initialized = .true.


!-----------------------------------------------------------------------


end subroutine cloud_rad_init


!######################################################################

subroutine cloud_summary3 (is, js, land,  use_fu2007, ql, qi, qa, qn, &
                           qni, pfull, phalf, tkel, Aerosol, nclds,   &
                           cldamt, lwp, iwp, reff_liq, reff_ice, &
                           ktop, kbot, conc_drop, conc_ice, size_drop,  &
                           size_ice, droplet_number, ice_number)
   
!---------------------------------------------------------------------
!   cloud_summary3 returns the specification properties of the clouds
!   present in the strat_cloud_mod.
!
!   cloud_summary3 returns the cloud properties that are needed by the
!   radiation package. If microphysically-based properties are desired,
!   concentrations (g/m**3) and particle size (microns) are returned for
!   liquid and ice cloud particles. If bulk physics are being used, then
!   cloud thickness and the number and area of random and max-overlapped 
!   clouds are returned.
!---------------------------------------------------------------------
 
integer,                   intent(in)            :: is,js
real, dimension(:,:),      intent(in)            :: land
logical,                   intent(in)            :: use_fu2007
real, dimension(:,:,:),    intent(in)            :: ql, qi, qa, qn, pfull,&
                                                    phalf, tkel, qni
type(aerosol_type),        intent(in)            :: Aerosol
integer, dimension(:,:),   intent(out)           :: nclds          
real, dimension(:,:,:),    intent(out)           :: cldamt, lwp, iwp, &
                                                    reff_liq, reff_ice
integer, dimension(:,:,:), intent(out), optional :: ktop, kbot 
real,    dimension(:,:,:), intent(out), optional :: conc_drop,conc_ice,&
                                                    size_drop,size_ice,&
                                                    droplet_number, &
                                                    ice_number

!---------------------------------------------------------------------
!    intent(in) variables:
!
!       is,js        Indices for model slab
!       land         Fraction of the grid box covered by land
!                    [ dimensionless ]
!       ql           Cloud liquid condensate [ kg condensate/kg air ]
!       qi           Cloud ice condensate [ kg condensate/kg air ]
!       qa           Cloud volume fraction [ fraction ]
!       qn           Cloud droplet number [ #/kg air]
!       qni          Ice particle number [ #/kg air]
!       pfull        Pressure at full levels [ Pascals ]
!       phalf        Pressure at half levels [ Pascals ]
!                    NOTE: it is assumed that phalf(j+1) > phalf(j)
!       tkel         Temperature [ deg. Kelvin ] 
!       Aerosol
!
!    intent(out) variables:
!
!       nclds        Number of random-overlap clouds in a column
!       cldamt       Cloud amount of condensed cloud
!       lwp          Liquid water path 
!       iwp          Ice water path
!       reff_liq     Effective radius of cloud drops
!       reff_ice     Effective radius of ice crystals
!
!   intent(out), optional variables:
! 
!       ktop         Integer level for top of cloud, present when 
!                    max-random overlap assumption made
!       kbot         Integer level for bottom of cloud, present when
!                    max-random overlap assumption made
!       conc_drop    Liquid cloud droplet mass concentration, present 
!                    when microphysically-based cloud radiative
!                    properties are desired
!       conc_ice     Ice cloud mass concentration, present when
!                    microphysically-based cloud radiative
!                    properties are desired
!       size_drop    Effective diameter of liquid cloud droplets, 
!                    present when microphysically-based cloud radiative
!                    properties are desired
!       size_ice     Effective diameter of ice cloud, present when 
!                    microphysically-based cloud radiative
!                    properties are desired
!       droplet_number
!                    number of cloud droplets [ # / kg(air) ]
!       ice_number
!                    number of cloud ice particles [ # / kg(air) ]
!
!--------------------------------------------------------------------
 
!--------------------------------------------------------------------
!    local variables:

      real,dimension (size(ql,1),size(ql,2),   &
                                 size(ql,3)) :: qa_local, ql_local, &
                                                qi_local, N_drop

      real,dimension (size(ql,1),size(ql,2)) :: k_ratio
      integer  :: i, j, k
      real     :: depth, aerosols_concen

!--------------------------------------------------------------------
!    local variables:
!
!       qa_local     local value of qa (fraction)
!       ql_local     local value of ql (kg condensate / kg air)
!       qi_local     local value of qi (kg condensate / kg air)
!       N_drop       number of cloud droplets per cubic meter (or per kg)
!       k_ratio      ratio of effective radius to mean volume radius
!       depth        depth of model layer
!       aerosols_concen
!                    aerosol concentration used to determine activated
!                    cloud droplet number (lin microphysics)
!
!---------------------------------------------------------------------

!--------------------------------------------------------------------
!    create local values of ql and qi. this step is necessary to remove 
!    the values of (qi,ql) which are 0 < (qi,ql) < qmin   or
!    (qi,ql) > qmin and qa <= qamin.
!--------------------------------------------------------------------

      do k=1, size(ql,3)
        do j=1, size(ql,2)
          do i=1, size(ql,1)
            qa_local(i,j,k) = 0.
            if ((qa(i,j,k) > qamin) .and. (ql(i,j,k) > qmin) ) then
              ql_local(i,j,k) = ql(i,j,k)
              qa_local(i,j,k) = qa(i,j,k)
            else
              ql_local(i,j,k) = 0.
            endif      
            if ((qa(i,j,k) > qamin) .and. (qi(i,j,k) > qmin) ) then
              qi_local(i,j,k) = qi(i,j,k)
              qa_local(i,j,k) = qa(i,j,k)
            else
              qi_local(i,j,k) = 0.
            endif       
          end do
        end do
      end do

!--------------------------------------------------------------------
!    define the cloud droplet concentration. 
      if (do_liq_num) then
        N_drop = qn    
      else
        do k=1, size(ql,3)
          N_drop(:,:,k)  = N_land*land(:,:) + N_ocean*(1. - land(:,:))
        end do
      endif

!--------------------------------------------------------------------
!    define the ratio of the effective drop radius to the mean volume 
!    radius (k_ratio). 
!--------------------------------------------------------------------
      k_ratio(:,:) = k_land*land(:,:) + k_ocean*(1. - land(:,:))

!----------------------------------------------------------------------
!    define optional output argument droplet_number (in units of #/kg). 
!----------------------------------------------------------------------
      if (present(droplet_number)) then
        if (do_liq_num) then
          droplet_number = qn
        else 
          do k=1, size(ql,3)
            do j=1, size(ql,2)
              do i=1, size(ql,1)
                droplet_number(i,j,k) = N_drop(i,j,k)/(pfull(i,j,k)/  &
                                           (RDGAS*tkel(i,j,k)))
              end do
            end do
          end do
        endif    
      endif

!----------------------------------------------------------------------
!    define optional output argument ice_number (in units of #/kg).
!----------------------------------------------------------------------
      if (present(ice_number)) then
        if ( do_ice_num ) then
          ice_number = qni   ! #/kg
        else
          ice_number = -999.
        end if
      endif

!--------------------------------------------------------------------
!    execute the following when  the max-random overlap assumption 
!    is being made. 
!--------------------------------------------------------------------
      if (present(ktop) .and. present(kbot)) then    ! max-rnd

!--------------------------------------------------------------------
!    if microphysics output is required, only the random overlap assump-
!    tion is allowed; if max-random overlap is requested, an error
!    message will be issued. if random overlap is requested, call
!    subroutine rnd_overlap to obtain the cloud specification proper-
!    ties, including the microphysical parameters.
!--------------------------------------------------------------------
        if (present (conc_drop) .and.  present (conc_ice ) .and. &
            present (size_ice ) .and.  present (size_drop)) then      
          call error_mesg ( 'cloud_rad_mod', &
           ' max-random overlap not currently available for radiation '//&
              'scheme requiring microphysically-based outputs', FATAL)
     
!----------------------------------------------------------------------
!    if some but not all of the microphysics variables are present,
!    stop execution.
!---------------------------------------------------------------------
        else if (present (conc_drop) .or.  present (conc_ice ) .or. &
                 present (size_ice ) .or.  present (size_drop)) then
          call error_mesg ('cloud_rad_mod', &
                ' if any microphysical args present, all must be '//&
                                                    'present', FATAL)

        else
          call  max_rnd_overlap (ql_local, qi_local, qa_local, pfull,  &
                                 phalf, tkel, N_drop, k_ratio, nclds,  &
                                 ktop, kbot, cldamt, lwp, iwp,   &
                                 reff_liq, reff_ice, qni  )
        endif
     
!---------------------------------------------------------------------
!    if only one of ktop and kbot is present, stop execution; both are 
!    needed for max-random overlap and neither are permitted when the 
!    random overlap assumption is made.
!---------------------------------------------------------------------
      else if (present(ktop) .or. present(kbot)) then ! error
        call error_mesg ('cloud_rad_mod',  &
                  'kbot and ktop must either both be absent or both '//&
                    'be present', FATAL)

!---------------------------------------------------------------------
!    if neither are present, then random overlap is assumed.
!---------------------------------------------------------------------
      else                 

!---------------------------------------------------------------------
!    if microphysical properties are desired, call subroutine 
!    rnd_overlap to obtain the cloud specification properties, including
!    the microphysical parameters.
!--------------------------------------------------------------------
        if (present (conc_drop) .and.  present (conc_ice ) .and. &
            present (size_ice ) .and.  present (size_drop)) then      
          call rnd_overlap (ql_local, qi_local, qa_local,  &
                            use_fu2007, pfull, phalf, tkel,  &
                            N_drop, qni,   k_ratio, nclds, &
                            lwp, iwp, reff_liq, reff_ice,   &
                            conc_drop_org=conc_drop,&
                            conc_ice_org =conc_ice,&
                            size_drop_org=size_drop,&
                            size_ice_org =size_ice)
          cldamt = qa_local

!--------------------------------------------------------------------
!    account for the plane-parallel homogeneous cloud bias.
!--------------------------------------------------------------------
          conc_drop = scale_factor*conc_drop
          conc_ice  = scale_factor*conc_ice 

!----------------------------------------------------------------------
!    if some but not all of the microphysics variables are present,
!    stop execution.
!---------------------------------------------------------------------
        else if (present (conc_drop) .or.  present (conc_ice ) .or. &
                 present (size_ice ) .or.  present (size_drop)) then   
          call error_mesg ('cloud_rad_mod', &
                ' if any microphysical args present, all must '//&
                                                'be present', FATAL)

!----------------------------------------------------------------------
!    if microphysics terms are not required, call rnd_overlap to obtain
!    the cloud specification variables from a bulk physics formulation.
!----------------------------------------------------------------------
        else
           call  rnd_overlap (ql_local, qi_local, qa_local,  &
                              use_fu2007, pfull, phalf, tkel,  &
                              N_drop, qni, k_ratio, nclds,  &
                              lwp, iwp, reff_liq, reff_ice)
           cldamt = qa_local
        endif
      endif ! (present(ktop and kbot))

!---------------------------------------------------------------------
    


end subroutine cloud_summary3



!###################################################################

subroutine cloud_rad_end

!------------------------------------------------------------------------
!    mark the module as uninitialized.
!    A destructor routine for the cloud_rad module.                                                                                                        !------------------------------------------------------------------------
      call gamma_mg_end

      module_is_initialized = .false.

!------------------------------------------------------------------------


end subroutine cloud_rad_end




!######################################################################
subroutine max_rnd_overlap (ql, qi, qa, pfull, phalf, tkel, N_drop,  &
                           k_ratio, nclds, ktop, kbot, cldamt, lwp,  &
                           iwp, reff_liq, reff_ice, qni)

!----------------------------------------------------------------------
!    max_rnd_overlap returns various cloud specification properties
!    obtained with the maximum-random overlap assumption.
!----------------------------------------------------------------------
 
real,    dimension(:,:,:), intent(in)             :: ql, qi, qa,  &
                                                     pfull, phalf, tkel, &
                                                     N_drop, qni  
real,    dimension(:,:),   intent(in)             :: k_ratio
integer, dimension(:,:),   intent(out)            :: nclds
integer, dimension(:,:,:), intent(out)            :: ktop, kbot
real,    dimension(:,:,:), intent(out)            :: cldamt, lwp, iwp, &
                                                     reff_liq, reff_ice

!---------------------------------------------------------------------
!   intent(in) variables:
!
!       ql           Cloud liquid condensate [ kg condensate/kg air ]
!       qi           Cloud ice condensate [ kg condensate/kg air ]
!       qa           Cloud volume fraction [ fraction ]
!       pfull        Pressure at full levels [ Pascals ]
!       phalf        Pressure at half levels, index 1 at model top 
!                    [ Pascals ]
!       tkel         Temperature [ deg Kelvin ]
!       N_drop       Number of cloud droplets per cubic meter
!       qni          Number of ice particles per kg
!       k_ratio      Ratio of effective radius to mean volume radius
!
!   intent(out) variables:
!
!       nclds        Number of (random overlapping) clouds in column 
!       ktop         Level of the top of the cloud
!       kbot         Level of the bottom of the cloud
!       cldamt       Cloud amount of condensed cloud [ dimensionless ]
!       lwp          Cloud liquid water path [ kg condensate / m **2 ]
!       iwp          Cloud ice path [ kg condensate / m **2 ]
!       reff_liq     Effective radius for liquid clouds [ microns ]
!       reff_ice     Effective particle size for ice clouds [ microns ]
!
!---------------------------------------------------------------------

!--------------------------------------------------------------------
!   local variables:

      real, dimension (size(ql,1), size(ql,2), size(ql,3))  :: &
                    cldamt_cs, lwp_cs, iwp_cs, reff_liq_cs, reff_ice_cs 

      integer    :: kdim
      integer    :: top_t, bot_t
      integer    :: tmp_top, tmp_bot, nlev
      logical    :: already_in_cloud, cloud_bottom_reached
      real       :: sum_liq, sum_ice, maxcldfrac
      real       :: totcld_bot, max_bot
      real       :: totcld_top, max_top, tmp_val
      real       :: reff_liq_local, sum_reff_liq
      real       :: reff_ice_local, sum_reff_ice
      integer    :: i, j, k, kc, t

      real       :: dumc, dumnc, rho, pgam, lamc, lammax, lammin, &
                    dumi, dumni, lami
!--------------------------------------------------------------------
!   local variables:
!
!       kdim              number of model layers
!       top_t             used temporarily as tag for cloud top index
!       bot_t             used temporarily as tag for cloud bottom index
!       tmp_top           used temporarily as tag for cloud top index
!       tmp_bot           used temporarily as tag for cloud bottom index
!       nlev              number of levels in the cloud
!       already_in_cloud  if true, previous layer contained cloud
!       cloud_bottom_reached
!                         if true, the cloud-free layer beneath a cloud
!                         has been reached
!       sum_liq           sum of liquid in cloud 
!                         [ kg condensate / m**2 ]
!       sum_ice           sum of ice in cloud 
!                         [ kg condensate / m**2 ]
!       maxcldfrac        maximum cloud fraction in any layer of cloud
!                         [ fraction ]
!       totcld_bot        total cloud fraction from bottom view
!       max_bot           largest cloud fraction face from bottom view
!       totcld_top        total cloud fraction from top view
!       max_top           largest cloud fraction face from top view
!       tmp_val           temporary number used in the assigning of top 
!                         and bottom
!       reff_liq_local    gridpoint value of reff of liquid clouds 
!                         [ microns ]
!       sum_reff_liq      condensate-weighted sum over cloud of 
!                         reff_liq_local  
!                         [ (kg condensate / m**2) * microns ]
!       reff_ice_local    gridpoint value ofreff of ice clouds  
!                         [ microns ]
!       sum_reff_ice      condensate-weighted sum over cloud of 
!                         reff_ice_local 
!                         [ (kg condensate / m**2) * microns ]
!       i,j,k,kc,t        do-loop indices
!
!----------------------------------------------------------------------

!---------------------------------------------------------------------
!    define the number of vertical layers in the model. initialize the
!    output fields to correspond to the absence of clouds.
!---------------------------------------------------------------------
      kdim     = size(ql,3)
      nclds    = 0
      ktop     = 1
      kbot     = 0
      cldamt   = 0.
      lwp      = 0.
      iwp      = 0.
      reff_liq = 10.
      reff_ice = 30.

!--------------------------------------------------------------------
!    find the levels with cloud in each column. determine the vertical
!    extent of each individual cloud, treating cloud in adjacent layers
!    as components of a multi-layer cloud, and then calculate appropr-
!    iate values of water paths and effective particle size.
!--------------------------------------------------------------------

      do j=1,size(ql,2)
        do i=1,size(ql,1)

!--------------------------------------------------------------------
!    set a flag indicating that we are searching for the next cloud top.
!--------------------------------------------------------------------
          already_in_cloud  = .false.
          cloud_bottom_reached = .false.

!--------------------------------------------------------------------
!    march down the column.
!--------------------------------------------------------------------
          do k=1,kdim      

!--------------------------------------------------------------------
!    find a layer containing cloud in the column. 
!--------------------------------------------------------------------
            if ( (ql(i,j,k) .gt. qmin) .or. &
                 (qi(i,j,k) .gt. qmin) ) then      

!--------------------------------------------------------------------
!    if the previous layer was not cloudy, then a new cloud has been
!    found. increment the cloud counter, set the flag to indicate the 
!    layer is in a cloud, save its cloud top level, initialize the 
!    values of its ice and liquid contents and fractional area and 
!    effective crystal and drop sizes. 
!--------------------------------------------------------------------
              if (.not. already_in_cloud)  then
                nclds(i,j) = nclds(i,j) + 1
                already_in_cloud = .true.
                cloud_bottom_reached = .false.
                ktop(i,j,nclds(i,j)) = k
                sum_liq          = 0.
                sum_ice          = 0.
                maxcldfrac       = 0.
                sum_reff_liq     = 0.
                sum_reff_ice     = 0.        
              endif

              if (ql(i,j,k) .gt. qmin) then
                call define_liquid_particle_size   & 
                          (k, ql(i,j,k), pfull(i,j,k), k_ratio(i,j),    &
                               qa(i,j,:), tkel(i,j,k), N_drop(i,j,k),    &
                                                        reff_liq(i,j,k))
              endif

              if (qi(i,j,k) .gt. qmin) then
                call define_ice_particle_size    &
                          (.FALSE., .FALSE., qi(i,j,k), qa(i,j,k),    &
                                qni(i,j,k), tkel(i,j,k), reff_ice(i,j,k))
              endif

!---------------------------------------------------------------------
!    add this layer's contributions to the current cloud. total liquid
!    content, ice content, largest cloud fraction and condensate-
!    weighted effective droplet and crystal radii are accumulated over 
!    the cloud.
!---------------------------------------------------------------------
              sum_liq = sum_liq + ql(i,j,k)*  &
                        (phalf(i,j,k+1) - phalf(i,j,k))/GRAV
              sum_ice = sum_ice + qi(i,j,k)* &
                        (phalf(i,j,k+1) - phalf(i,j,k))/GRAV
              maxcldfrac = MAX(maxcldfrac,qa(i,j,k))
              sum_reff_liq  = sum_reff_liq + (reff_liq(i,j,k)*ql(i,j,k)*&
                              (phalf(i,j,k+1) - phalf(i,j,k))/GRAV)
              sum_reff_ice  = sum_reff_ice + &
                              (reff_ice(i,j,k) * qi(i,j,k) * &
                              (phalf(i,j,k+1) - phalf(i,j,k))/GRAV)
            endif ! (ql > qmin or qi > qmin)

!--------------------------------------------------------------------
!    when the cloud-free layer below a cloud is reached, or if the
!    bottom model level is reached, define the cloud bottom level and
!    set a flag indicating that mean values for the cloud may now be
!    calculated.
!--------------------------------------------------------------------
            if (ql(i,j,k) <= qmin .and. qi(i,j,k) <= qmin .and. &
                already_in_cloud) then                 
              cloud_bottom_reached = .true.
              kbot(i,j,nclds(i,j)) = k - 1
            else if (already_in_cloud .and. k == kdim) then
              cloud_bottom_reached = .true.
              kbot(i,j,nclds(i,j)) = kdim
            endif

!--------------------------------------------------------------------
!    define the cloud fraction as the largest value of any layer in the
!    cloud. define the water paths as the total liquid normalized by the
!    fractional area of the cloud. define the condensate-weighted 
!    effective water and ice radii. 
!--------------------------------------------------------------------
            if (cloud_bottom_reached) then
              cldamt_cs(i,j,nclds(i,j)) = maxcldfrac
              lwp_cs(i,j,nclds(i,j)) = sum_liq/cldamt_cs(i,j,nclds(i,j))
              iwp_cs(i,j,nclds(i,j)) = sum_ice/cldamt_cs(i,j,nclds(i,j))
              if (sum_liq > 0.) then
                reff_liq_cs(i,j,nclds(i,j)) = sum_reff_liq/sum_liq
              else
                reff_liq_cs(i,j,nclds(i,j)) = 10.0
              end if
              if (sum_ice > 0.) then
                reff_ice_cs(i,j,nclds(i,j)) = sum_reff_ice/sum_ice
              else
                reff_ice_cs(i,j,nclds(i,j)) = 30.0
              end if

!----------------------------------------------------------------------
!    if adjust_top is true, the top and bottom indices of multi-layer
!    clouds are adjusted to be those that are the most exposed to top 
!    and bottom view.
!----------------------------------------------------------------------
              if (adjust_top) then
    
!---------------------------------------------------------------------
!    define the cloud thickness.
!---------------------------------------------------------------------
                nlev = kbot(i,j,nclds(i,j)) - ktop(i,j,nclds(i,j)) + 1
                if (nlev > 1) then

!---------------------------------------------------------------------
!    use the current top and bottom as the first guess for the new 
!    values.
!---------------------------------------------------------------------
                  tmp_top = ktop(i,j,nclds(i,j))
                  tmp_bot = kbot(i,j,nclds(i,j))

!--------------------------------------------------------------------
!    initialize local search variables.
!--------------------------------------------------------------------
                  totcld_bot = 0.
                  totcld_top = 0.
                  max_bot    = 0.
                  max_top    = 0.
          
!--------------------------------------------------------------------
!    to find the adjusted cloud top, begin at current top and work 
!    downward. find the layer which is most exposed when viewed from
!    the top; i.e., the cloud fraction increase is largest for that
!    layer. the adjusted cloud base is found equivalently, starting
!    from the actual cloud base and working upwards.
!--------------------------------------------------------------------
                  do t=1,nlev

!--------------------------------------------------------------------
!    find adjusted cloud top.
!--------------------------------------------------------------------
                    top_t   = ktop(i,j,nclds(i,j)) + t - 1
                    tmp_val = MAX(0., qa(i,j,top_t) - totcld_top)
                    if (tmp_val > max_top) then
                      max_top = tmp_val
                      tmp_top = top_t
                    end if
                    totcld_top = totcld_top + tmp_val         
                              
!--------------------------------------------------------------------
!    find adjusted cloud base.
!--------------------------------------------------------------------
                    bot_t   = kbot(i,j,nclds(i,j)) - t + 1
                    tmp_val = MAX(0., qa(i,j,bot_t) - totcld_bot)
                    if (tmp_val > max_bot) then
                      max_bot = tmp_val
                      tmp_bot = bot_t
                    end if
                    totcld_bot = totcld_bot + tmp_val         
                  end do
                       
!--------------------------------------------------------------------
!    assign tmp_top and tmp_bot as the new ktop and kbot.
!--------------------------------------------------------------------
                  ktop(i,j,nclds(i,j)) = tmp_top
                  kbot(i,j,nclds(i,j)) = tmp_bot
                endif  !(nlev > 1)  
              endif  ! (adjust_top)

!---------------------------------------------------------------------
!    reset already_in_cloud and cloud_bottom_reached to indicate that
!    the current cloud has been exited.
!---------------------------------------------------------------------
              already_in_cloud     = .false.
              cloud_bottom_reached = .false.
            endif   ! (cloud_bottom_reached)
          end do
        end do
      end do

!---------------------------------------------------------------------
!    place cloud properties into physical-space arrays for return to
!    calling routine. NOTE THAT ALL LEVELS IN A GIVEN CLOUD ARE
!    ASSIGNED THE SAME PROPERTIES.
!---------------------------------------------------------------------
      do j=1,size(ql,2)
        do i=1,size(ql,1)
          do kc=1, nclds(i,j)
            do k= ktop(i,j,kc), kbot(i,j,kc)
              cldamt(i,j,k)   = cldamt_cs(i,j,kc)
              lwp(i,j,k)      = lwp_cs(i,j,kc)
              iwp(i,j,k)      = iwp_cs(i,j,kc)
              reff_liq(i,j,k) = reff_liq_cs(i,j,kc)
              reff_ice(i,j,k) = reff_ice_cs(i,j,kc)
            end do
          end do
        end do
      end do
     
!---------------------------------------------------------------------

end subroutine max_rnd_overlap




!#####################################################################
subroutine rnd_overlap    (ql, qi, qa, use_fu2007, pfull, phalf,   &
                           tkel, N_drop, qni, k_ratio, nclds, lwp, iwp,  &
                           reff_liq, reff_ice, conc_drop_org,  &
                           conc_ice_org, size_drop_org, size_ice_org)

!----------------------------------------------------------------------
!    rnd_overlap returns various cloud specification properties, 
!    obtained with the random-overlap assumption. implicit in this
!    asusmption is that all clouds are only a single layer thick; i.e.,
!    clouds at adjacent levels in the same column are independent of
!    one another.
!----------------------------------------------------------------------
 
real,    dimension(:,:,:), intent(in)             :: ql, qi, qa,  &
                                                     pfull, phalf, tkel, &
                                                     N_drop, qni  
logical,                   intent(in)             :: use_fu2007
real,    dimension(:,:),   intent(in)             :: k_ratio
integer, dimension(:,:),   intent(out)            :: nclds
real,    dimension(:,:,:), intent(out)            :: lwp, iwp, &
                                                     reff_liq, reff_ice
real,    dimension(:,:,:), intent(out), optional  :: conc_drop_org,  &
                                                     conc_ice_org,  &
                                                     size_drop_org,  &
                                                     size_ice_org

!---------------------------------------------------------------------
!   intent(in) variables:
!
!       ql           Cloud liquid condensate [ kg condensate/kg air ]
!       qi           Cloud ice condensate [ kg condensate/kg air ]
!       qa           Cloud volume fraction [ fraction ]
!       pfull        Pressure at full levels [ Pascals ]
!       phalf        Pressure at half levels, index 1 at model top 
!                    [ Pascals ]
!       tkel         Temperature [ deg Kelvin ]
!       N_drop       Number of cloud droplets per cubic meter
!       qni          Number of ice particles per kg
!       k_ratio      Ratio of effective radius to mean volume radius
!
!   intent(out) variables:
!
!       nclds        Number of (random overlapping) clouds in column 
!       lwp          Cloud liquid water path [ kg condensate / m **2 ]
!       iwp          Cloud ice path [ kg condensate / m **2 ]
!       reff_liq     Effective radius for liquid clouds [ microns ]
!       reff_ice     Effective particle size for ice clouds [ microns ]
!
!    intent(out), optional variables:
!
!       conc_drop_org Liquid cloud droplet mass concentration 
!                     [ g / m**3 ]
!       conc_ice_org  Ice cloud mass concentration [ g / m**3 ]
!       size_drop_org Effective diameter of liquid cloud droplets 
!                     [ microns ]
!       size_ice_org  Effective diameter of ice clouds { microns ]
!
!---------------------------------------------------------------------

!--------------------------------------------------------------------
!   local variables:

      logical    ::  want_microphysics
      integer    ::  i, j, k

!--------------------------------------------------------------------
!   local variables:
!
!       want_microphysics   logical indicating if microphysical 
!                           parameters are to be calculated
!       i,j,k               do-loop indices
!
!----------------------------------------------------------------------

!--------------------------------------------------------------------
!    define a logical, indicating the presence (or not) of the optional
!    microphysics arguments.
!--------------------------------------------------------------------
      if (present(conc_drop_org) .and. present(conc_ice_org ) .and. &
          present(size_ice_org ) .and. present(size_drop_org)) then  
        want_microphysics = .true.

!----------------------------------------------------------------------
!    if some but not all of the microphysics variables are present,
!    stop execution.
!---------------------------------------------------------------------
      else if (present(conc_drop_org) .or. present(conc_ice_org ) .or. &
               present(size_ice_org ) .or. present(size_drop_org)) then 
        call error_mesg ('cloud_rad_mod', &
            ' if any microphysical args present, all must be present',&
                                                                FATAL)

!----------------------------------------------------------------------
!    if the optional arguments are not present, set the appropriate
!    flag to indicate that only bulk properties will be calculated.
!---------------------------------------------------------------------
      else
         want_microphysics = .false.
      end if

!--------------------------------------------------------------------
!    count the layers with cloud in each column, starting at model top. 
!--------------------------------------------------------------------
      nclds = 0
      do k=1,size(ql,3)
        do j=1,size(ql,2)
          do i=1,size(ql,1)
            if (qa(i,j,k) .gt. qmin ) then
               
!---------------------------------------------------------------------
!    when cloud is found, increment the cloud column counter.
!---------------------------------------------------------------------
              nclds(i,j) = nclds(i,j) + 1
            endif
          end do
        end do
      end do

!---------------------------------------------------------------------
!    if cloud water is present (> qmin) compute cloud water path and 
!    if a microphysics-based scheme is active, the droplet concentration.
!    if cloud ice is present (> qmin) compute cloud ice path and 
!    if a microphysics-based scheme is active, the ice concentration.
!---------------------------------------------------------------------
      if(want_microphysics) then
        do k=1,size(ql,3)
          do j=1,size(ql,2)
            do i=1,size(ql,1)
!---------------------------------------------------------------------
!    if liquid water is present, compute the liquid water path. 
!---------------------------------------------------------------------
              if (ql(i,j,k) .gt. qmin) then
                lwp(i,j,k) = ql(i,j,k)*(phalf(i,j,k+1) - phalf(i,j,k))/ &
                                                              GRAV/qa(i,j,k)
!----------------------------------------------------------------------
!    if microphysical properties are desired, calculate the droplet
!    concentration. units of concentration are g / m**3.
!----------------------------------------------------------------------
                conc_drop_org(i,j,k) =     &
                      1000.*ql(i,j,k)*(phalf(i,j,k+1) - phalf(i,j,k))/&
                      RDGAS/tkel(i,j,k)/log(phalf(i,j,k+1)/  &
                      MAX(phalf(i,j,k), pfull(i,j,1)))/qa(i,j,k)
              else
                lwp(i,j,k) = 0.
                conc_drop_org(i,j,k) = 0.
              endif
!---------------------------------------------------------------------
!    if ice is present, compute the ice water path.
!---------------------------------------------------------------------
              if (qi(i,j,k) .gt. qmin) then
                iwp(i,j,k) = qi(i,j,k)*(phalf(i,j,k+1) - phalf(i,j,k))/ &
                                                             GRAV/qa(i,j,k)
!----------------------------------------------------------------------
!    if microphysical properties are desired, calculate the ice
!    concentration. units of concentration are in g / m**3.
!----------------------------------------------------------------------
                conc_ice_org (i,j,k) =     &
                      1000.*qi(i,j,k)*(phalf(i,j,k+1) - phalf(i,j,k))/ &
                      RDGAS/tkel(i,j,k)/log(phalf(i,j,k+1)/   &
                      MAX(phalf(i,j,k), pfull(i,j,1)))/ qa(i,j,k)
              else
                iwp(i,j,k) = 0.
                conc_ice_org(i,j,k) = 0.
              endif            
            end do
          end do
        end do

!------------------------------------------------------------------------
!    call define_liquid_particle_size to compute the effective radius 
!    and /or the effective drop size diameter to be used by the radiation 
!    code.
!------------------------------------------------------------------------
        do k=1,size(ql,3)
          do j=1,size(ql,2)
            do i=1,size(ql,1)
              if (ql(i,j,k) .gt. qmin) then
                call define_liquid_particle_size (k,  ql(i,j,k),   &
                     pfull(i,j,k), k_ratio(i,j), qa(i,j,:), tkel(i,j,k), &
                                            N_drop(i,j,k),  reff_liq(i,j,k))
!               size_drop_org(i,j,k) = 2.*reff_liq(i,j,k)
!-------------------------------------------------------------------------
!    if ql is not present (or below qmin), set the output fields to
!    appropriate values.
!-------------------------------------------------------------------------
              else
                reff_liq(i,j,k) = 10.
!               size_drop_org(i,j,k) = 20.
              endif  ! (ql > qmin)
            end do
          end do
        end do

        size_drop_org = 2.*reff_liq

!------------------------------------------------------------------------
!    call define_ice_particle_size to compute the effective radius 
!    and /or the effective crystal size to be used by the radiation 
!    code.
!------------------------------------------------------------------------
        do k=1,size(ql,3)
          do j=1,size(ql,2)
            do i=1,size(ql,1)
              if (qi(i,j,k) .gt. qmin) then
                call define_ice_particle_size (want_microphysics,   &
                    use_fu2007, qi(i,j,k), qa(i,j,k), qni(i,j,k),   &
                    tkel(i,j,k), reff_ice(i,j,k), size_ice_org(i,j,k))

!-----------------------------------------------------------------------
!    if insufficient ice is present, set output fields to appropriate
!    values.
!-----------------------------------------------------------------------
              else  ! (qi > qmin)
                reff_ice(i,j,k) = 30.
                size_ice_org(i,j,k) = 60.
              end if ! (qi > qmin)                    
            end do
          end do
        end do
      else ! want_microphysics
        do k=1,size(ql,3)
          do j=1,size(ql,2)
            do i=1,size(ql,1)
              if (ql(i,j,k) .gt. qmin) then
                lwp(i,j,k) = ql(i,j,k)*(phalf(i,j,k+1) - phalf(i,j,k))/ &
                                                              GRAV/qa(i,j,k) 
              else
                lwp(i,j,k) = 0.
              endif
              if (qi(i,j,k) .gt. qmin) then
                iwp(i,j,k) = qi(i,j,k)*(phalf(i,j,k+1) - phalf(i,j,k))/ &
                                                             GRAV/qa(i,j,k)
              else
                iwp(i,j,k) = 0.
              endif            
            end do
          end do
        end do

        do k=1,size(ql,3)
          do j=1,size(ql,2)
            do i=1,size(ql,1)
              if (ql(i,j,k) .gt. qmin) then
                call define_liquid_particle_size (k,  ql(i,j,k),   &
                     pfull(i,j,k), k_ratio(i,j), qa(i,j,:), tkel(i,j,k), &
                                            N_drop(i,j,k),  reff_liq(i,j,k))
              else
                reff_liq(i,j,k) = 10.
              endif  ! (ql > qmin)
            end do
          end do
        end do

      do k=1,size(ql,3)
        do j=1,size(ql,2)
          do i=1,size(ql,1)
            if (qi(i,j,k) .gt. qmin) then
              call define_ice_particle_size (want_microphysics,   &
                  use_fu2007, qi(i,j,k), qa(i,j,k), qni(i,j,k),   &
                  tkel(i,j,k), reff_ice(i,j,k), size_ice_org(i,j,k))
            else  ! (qi > qmin)
              reff_ice(i,j,k) = 30.
            end if ! (qi > qmin)
          end do
        end do
      end do
    endif ! want_microphysics

!-------------------------------------------------------------------

end subroutine rnd_overlap   



!########################################################################

SUBROUTINE snow_and_rain (qa, pfull, phalf, tkel, cldamt, snow, rain,  &
                          size_snow_in, size_rain_in, conc_rain,  &
                          conc_snow,  size_rain, size_snow )

!----------------------------------------------------------------------
real, dimension(:,:,:), intent(in)     :: qa, pfull, phalf, tkel
real, dimension(:,:,:), intent(in)     :: snow, rain, size_snow_in,  &
                                          size_rain_in
real, dimension(:,:,:), intent(inout)  :: cldamt
real, dimension(:,:,:), intent(out)    :: conc_rain, conc_snow, &
                                          size_rain, size_snow

! NOTE: size_snow currently not used 
!--------------------------------------------------------------------


      REAL, DIMENSION( size(qa,1), size(qa,2), size(qa,3))  :: cldmax_loc

      REAL    :: dum
      integer :: i,j,k


      conc_snow = 0.
      conc_rain = 0.
      size_snow = 1.e-20
      size_rain = 1.e-20

 

      IF (snow_in_cloudrad .OR. rain_in_cloudrad ) THEN

!--------------------------------------------------------------------
!    define the locally seen cloud max above each level.
!--------------------------------------------------------------------
        cldmax_loc = 0.
        do k=1,size(qa,3)
          do j=1,size(qa,2)
            do i=1,size(qa,1)
! this might be o.k. as long as stochastic clouds are used  
              if (k.eq.1) then
                cldmax_loc(i,j,k)  = qa(i,j,k) !max overlap for precip
              else
                cldmax_loc(i,j,k) = max(cldmax_loc(i,j,k-1), qa(i,j,k))
              end if
            end do
          end do
        end do

!-----------------------------------------------------------------------
!    define the snow field to be seen by the radiation code.
!-----------------------------------------------------------------------
        IF (snow_in_cloudrad) THEN
          do k=1,size(qa,3) 
            do j=1,size(qa,2)
              do i=1,size(qa,1)
                IF (cldmax_loc(i,j,k) .GE. qamin) THEN
                  dum =  snow(i,j,k)/cldmax_loc(i,j,k)
                  dum = MIN(dum, 60.e-3)
                  IF (dum .GE. qmin) THEN
                    conc_snow (i,j,k) =     &
                       1000.*dum*(phalf(i,j,k+1) - phalf(i,j,k))/ &
                        RDGAS/tkel(i,j,k)/log(phalf(i,j,k+1)/   &
                        MAX(phalf(i,j,k), pfull(i,j,1)))
!size_snow currently not used
!!                  size_snow (i,j,k) =  MAX( 1.e-20, size_snow_in( i,j,k))
!!                  cldamt(i,j,k) = cldmax_loc(i,j,k)
                  END IF
                END IF 
              end do
            end do
          end do
        END IF  

!-----------------------------------------------------------------------
!    define the rain field to be seen by the radiation code.
!-----------------------------------------------------------------------
        IF (rain_in_cloudrad) THEN
          do k=1,size(qa,3)
            do j=1,size(qa,2)
              do i=1,size(qa,1)
                IF (cldmax_loc(i,j,k) .GE. qamin) THEN
                  dum = rain(i,j,k)/cldmax_loc(i,j,k)
                  dum = MIN(dum, 60.e-3)
                  IF (dum .GE. qmin) THEN
                    conc_rain (i,j,k) =     &
                       1000.*dum*(phalf(i,j,k+1) - phalf(i,j,k))/ &
                       RDGAS/tkel(i,j,k)/log(phalf(i,j,k+1)/   &
                       MAX(phalf(i,j,k), pfull(i,j,1)))
                    size_rain (i,j,k) = MAX( 1.e-20, size_rain_in(i,j,k))
!!                  cldamt(i,j,k) = cldmax_loc(i,j,k)
                  END IF
                END IF 
              end do
            end do
          end do
        END IF 
      END IF !(snow_in_cloudrad .OR. rain_in_cloudrad ) 

!-----------------------------------------------------------------------

END SUBROUTINE snow_and_rain


!########################################################################

subroutine define_liquid_particle_size (k, ql, pfull, k_ratio, qa, tkel, &
                                                  N_drop,  reff_liq_local)

!------------------------------------------------------------------------
integer,            intent(in)  :: k
real,               intent(in)  :: ql, pfull, tkel, N_drop
real, dimension(:), intent(in)  :: qa
real,               intent(in)  :: k_ratio
real,               intent(out) :: reff_liq_local
!------------------------------------------------------------------------

      integer :: kdim
      real :: dumc, dumnc, rho, pgam, lamc, lammax, lammin

!-----------------------------------------------------------------------

      kdim = size(qa,1)

!---------------------------------------------------------------------
!    compute the effective cloud droplet radius. 
!    when cloud droplets are not being predicted, the following formula
!    for liquid clouds is used, as recommended by 
!    Martin et al., J. Atmos. Sci, vol 51, pp. 1823-1842:
!
!    reff (in microns) =  k * 1.E+06 *
!                    (3*airdens*(ql/qa)/(4*pi*Dens_h2o*N_liq))**(1/3)
!
!    where airdens = density of air in kg air/m3
!               ql = liquid condensate in kg cond/kg air
!               qa = cloud fraction
!               pi = 3.14159
!         Dens_h2o = density of pure liquid water (kg liq/m3) 
!            N_liq = density of cloud droplets (number per cubic meter)
!                k = factor to account for difference between 
!                    mean volume radius and effective radius
!---------------------------------------------------------------------
      if (ql > qmin) then
        if (.not. do_liq_num) then
          reff_liq_local = k_ratio* 620350.49 *  &
                  (pfull*ql/qa(k)/RDGAS/tkel/DENS_H2O/N_drop  )**(1./3.)
        else 
          IF ( .NOT. do_ice_num ) THEN

!--------------------------------------------------------------------
!    use the following formula for  cloud drop effective radius for the 
!    case when cloud droplet number is being predicted, but ice particle 
!    number is not being predicted.
!
! yim: a variant for prognostic droplet number
!    reff (in microns) =  k * 1.E+06 *
!                    (3*ql/(4*pi*Dens_h2o*N_liq))**(1/3)
!
!    where airdens = density of air in kg air/m3
!               ql = liquid condensate in kg cond/kg air
!               qa = cloud fraction
!               pi = 3.14159
!         Dens_h2o = density of pure liquid water (kg liq/m3) 
!            N_liq = mixing ratio of cloud droplets (number/kg air)
!                k = factor to account for difference between 
!                    mean volume radius and effective radius
!--------------------------------------------------------------------
            if (ql > qmin) then
              reff_liq_local = k_ratio*620350.49*    &
                         (ql/DENS_H2O/max(N_drop,N_min*max(qa(k),qmin)/  &
                                            (pfull/RDGAS/tkel)))**(1./3.)
            else
              reff_liq_local = 0.0
            endif
!--------------------------------------------------------------------
!    use the following formula for cloud droplet effective radius in the 
!    case when both cloud droplet number and ice particle number are 
!    being predicted.
!-----------------------------------------------------------------------
          ELSE 

!-----------------------------------------------------------------------
!    calculate in-cloud mixing ratio and number concentration.
!    limit the in-cloud mixing ratio to reasonable value of 5 g kg-1.
!    put a lower limit on the droplet number concentration.
!-----------------------------------------------------------------------
            dumc = ql/qa(k)
            dumnc = N_drop/qa(k)
            dumc = min(dumc, 5.e-3)
            dumnc = max(dumnc, N_min)

!-----------------------------------------------------------------------
!    if sufficient liquid is present, compute the effective droplet
!    radius using the gamma function formulation.
!-----------------------------------------------------------------------
            if ( dumc > qmin ) then

!-----------------------------------------------------------------------
!    add upper limit to in-cloud number concentration to prevent 
!    numerical error.
!-----------------------------------------------------------------------
              dumnc = min(dumnc,dumc*1.e20)
              rho = pfull/(RDGAS*tkel)
              pgam = 0.0005714*(dumnc/1.e6/rho) + 0.2714
              pgam = 1./(pgam**2) - 1.
              pgam = max(pgam,2.)
              pgam = min(pgam,15.)

              lamc = (pi/6.*rhow*dumnc*gamma_mg(pgam + 4.)/ &
                              (dumc*gamma_mg(pgam + 1.)))**(1./3.)
              lammin = (pgam + 1.)/max_diam_drop
              lammax = (pgam + 1.)/min_diam_drop
              if (lamc.lt.lammin) then
                lamc = lammin
              else if (lamc.gt.lammax) then
                lamc = lammax
              end if
              reff_liq_local = gamma_mg(qcvar+1./3.)/   &
                                (gamma_mg(qcvar)*qcvar**(1./3.))* &
                                    gamma_mg(pgam + 4.)/ &
                                        gamma_mg(pgam + 3.)/lamc/2.*1.e6
            else
              reff_liq_local = 10.
            end if ! (dumc > qmin)
          END IF  ! (.not. do_ice_num)
        end if  ! (.not. do_liq_num) 

!----------------------------------------------------------------------
!    for single layer liquid or mixed phase clouds it is assumed that
!    cloud liquid is vertically stratified within the cloud.  under
!    such situations for observed stratocumulus clouds it is found
!    that the cloud mean effective radius is between 80 and 100% of
!    the cloud top effective radius. (Brenguier et al., Journal of
!    Atmospheric Sciences, vol. 57, pp. 803-821 (2000))  for linearly 
!    stratified cloud in liquid specific humidity, the cloud top 
!    effective radius is greater than the effective radius of the 
!    cloud mean specific humidity by a factor of 2**(1./3.).
!    this correction, 0.9*(2**(1./3.)) = 1.134, is applied only to 
!    single layer liquid or mixed phase clouds.
!    (as coded this effect is not applied if there is any chance of
!    cloud overlap in adjacent layers. Only if there is no cloud in
!    both adjacent layers is the effect applied.)
!----------------------------------------------------------------------
        if (do_brenguier) then
          if ( k == 1 ) then
            if (qa(2) < qamin) then
              reff_liq_local = 1.134*reff_liq_local
            endif
          else if (k == kdim ) then
            if ( qa(kdim-1) < qamin) then
              reff_liq_local = 1.134*reff_liq_local
            endif
          else if (qa(k-1) .lt. qamin .and. & 
                   qa(k+1) .lt. qamin)  then
            reff_liq_local = 1.134*reff_liq_local
          end if
        end if

!-------------------------------------------------------------------------
!    if ql is not present (or below qmin), set the output fields to
!    appropriate values.
!-------------------------------------------------------------------------
      else
        reff_liq_local = 0.
      endif  ! (ql > qmin)

!----------------------------------------------------------------------


end subroutine define_liquid_particle_size



!########################################################################

subroutine define_ice_particle_size ( want_microphysics, use_fu2007,   &
                          qi, qa, qni, tkel, reff_ice_local, size_ice_org)

!------------------------------------------------------------------------
logical, intent(in)            :: use_fu2007, want_microphysics
real,    intent(in)            :: qi, qa, qni, tkel
real,    intent(out)           :: reff_ice_local
real,    intent(out), optional :: size_ice_org
!-----------------------------------------------------------------------

      real :: dumi, dumni, lami, lammax, lammin


!---------------------------------------------------------------------
!    compute the effective ice crystal size. 
!    1) for bulk physics cases, the effective radius is taken from the 
!       formulation in Donner (1997, J. Geophys. Res., 102, 
!       pp. 21745-21768) which is based on Heymsfield and Platt (1984) 
!       with enhancement for particles smaller than 20 microns.  
!    2) if microphysical properties are requested, then the size of the
!       ice crystals comes from the Deff column [ reference ?? ].     
!    3) alternatively, if use_fu2007 is .true, Fu's parameterization of 
!       dge is used --yim
!    4) if  the ice particle number is being predicted the appropriate 
!       gamma function expression is used.

!              T Range (K)               Reff (microns)   Deff (microns)
!     -------------------------------    --------------   --------------
!
!     Tfreeze-25. < T                       92.46298         100.6
!     Tfreeze-30. < T <= Tfreeze-25.        72.35392          80.8
!     Tfreeze-35. < T <= Tfreeze-30.        85.19071          93.5
!     Tfreeze-40. < T <= Tfreeze-35.        55.65818          63.9
!     Tfreeze-45. < T <= Tfreeze-40.        35.29989          42.5
!     Tfreeze-50. < T <= Tfreeze-45.        32.89967          39.9
!     Tfreeze-55  < T <= Tfreeze-50         16.60895          21.6
!                   T <= Tfreeze-55.        15.41627          20.2
!---------------------------------------------------------------------

!---------------------------------------------------------------------
!    calculate the effective ice crystal size using the appropriate data.
!    if microphysics being used, define size_ice_org.
!---------------------------------------------------------------------
      if (qi .gt. qmin) then
        IF ( .not. do_ice_num ) THEN

          if (use_fu2007) then
            reff_ice_local = 47.05 +   &
                                0.6624*(tkel - TFREEZE) +&
                                0.001741*(tkel-TFREEZE)**2
            size_ice_org = reff_ice_local
          else ! (use_fu2007)
            if (want_microphysics) then
              if (tkel > TFREEZE - 25.) then
                reff_ice_local = 100.6      
              else if (tkel >  TFREEZE - 30. .and. &
                       tkel <= TFREEZE - 25.) then
                reff_ice_local = 80.8        
              else if (tkel >  TFREEZE - 35. .and. &
                       tkel <= TFREEZE - 30.) then
                reff_ice_local = 93.5       
              else if (tkel >  TFREEZE - 40. .and. &
                       tkel <= TFREEZE - 35.) then
                reff_ice_local = 63.9         
              else if (tkel >  TFREEZE - 45. .and. &
                       tkel <= TFREEZE - 40.) then
                reff_ice_local = 42.5       
              else if (tkel >  TFREEZE - 50. .and. &
                       tkel <= TFREEZE - 45.) then
                reff_ice_local = 39.9           
              else if (tkel >  TFREEZE - 55. .and. &
                       tkel <= TFREEZE - 50.) then
                reff_ice_local = 21.6         
              else
                reff_ice_local = 20.2             
              endif

              size_ice_org = reff_ice_local
            endif
          endif ! (use_fu2007)

!---------------------------------------------------------------------
!    calculate reff_ice using the bulk physics data when bulk physics
!    is being used.
!---------------------------------------------------------------------
          if (tkel > TFREEZE - 25.) then
            reff_ice_local = 92.46298
          else if (tkel >  TFREEZE - 30. .and. &
                   tkel <= TFREEZE - 25.) then
            reff_ice_local = 72.35392
          else if (tkel >  TFREEZE - 35. .and. &
                   tkel <= TFREEZE - 30.) then
            reff_ice_local = 85.19071 
          else if (tkel >  TFREEZE - 40. .and. &
                   tkel <= TFREEZE - 35.) then
            reff_ice_local = 55.65818
          else if (tkel >  TFREEZE - 45. .and. &
                   tkel <= TFREEZE - 40.) then
            reff_ice_local = 35.29989
          else if (tkel >  TFREEZE - 50. .and. &
                   tkel <= TFREEZE - 45.) then
            reff_ice_local = 32.89967
          else if (tkel >  TFREEZE - 55. .and. &
                   tkel <= TFREEZE - 50.) then
            reff_ice_local = 16.60895
          else
            reff_ice_local = 15.41627
          endif


!----------------------------------------------------------------------
!    when ice particle number is predicted use the gamma function 
!    expression from the Morrison Gettelman code.
!----------------------------------------------------------------------
        ELSE   ! (.not. do_ice_num)

!-----------------------------------------------------------------------
!    calculate in-cloud mixing ratio and number concentration.
!    limit the in-cloud mixing ratio to reasonable value of 5 g kg-1.
!    put a lower limit on the droplet number concentration.
!-----------------------------------------------------------------------
          dumi = qi/qa
          dumni = qni/qa
          dumi = min(dumi ,5.e-3)

!-----------------------------------------------------------------------
!    if sufficient ice is present, compute the cloud ice effective
!    radius using the gamma function formulation.
!-----------------------------------------------------------------------
          if (dumi > qmin) then

!-----------------------------------------------------------------------
!    add upper limit to in-cloud number concentration to prevent 
!    numerical error.
!-----------------------------------------------------------------------
            dumni = min(dumni, dumi*1.e20)

!-----------------------------------------------------------------------
!    define the effective radius as 3rd moment/2nd moment as in Fu 
!    and Liou paper.  This assumes spherical ice particles.
!-----------------------------------------------------------------------
            lami = (gamma_mg(1. + di_mg)*ci_mg*dumni/dumi)**(1./di_mg)
            lammax = 1./ min_diam_ice
            lammin = 1./(2.*dcs)
            if (lami.lt.lammin) then
              lami = lammin
            else if (lami.gt.lammax) then
              lami = lammax
            end if
            reff_ice_local = 1.5/lami*1.e6
          else
            reff_ice_local = 25.
          end if

          reff_ice_local = 2.*reff_ice_local 
          size_ice_org = reff_ice_local
        END IF  

!-----------------------------------------------------------------------
!    if insufficient ice is present, set output fields to appropriate
!    values.
!-----------------------------------------------------------------------
      else  ! (qi > qmin)
        reff_ice_local = 30.
        if (present(size_ice_org)) size_ice_org = 60.
      end if ! (qi > qmin)                    

!----------------------------------------------------------------------


end subroutine define_ice_particle_size


!#####################################################################

                  end module cloud_rad_mod

