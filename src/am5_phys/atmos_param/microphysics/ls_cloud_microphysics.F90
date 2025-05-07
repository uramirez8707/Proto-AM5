                    module ls_cloud_microphysics_mod

!-----------------------------------------------------------------------
!
!         interface module for cloud microphysics
!         ---------------------------------------
!         OPTIONS AVAILABLE:
!             Rotstayn-Klein microphysics (this used in old strat_cloud)
!             NCAR microphysics version 2.0 (became available? )
!
!-----------------------------------------------------------------------

! fms modules
use time_manager_mod,      only: time_type, get_time, set_date
use mpp_mod,               only: input_nml_file
use fms_mod,               only: error_mesg, FATAL, NOTE,        &
                                 check_nml_error,    &
                                 write_version_number, stdlog,   &
                                 mpp_pe, mpp_root_pe, stdlog,    &
                                 mpp_clock_id, mpp_clock_begin,  &
                                 mpp_clock_end, CLOCK_MODULE,    &
                                 CLOCK_MODULE_DRIVER, &
                                 MPP_CLOCK_SYNC
use field_manager_mod,     only: MODEL_ATMOS
use tracer_manager_mod,    only: get_tracer_index,&
                                 get_number_tracers, &
                                 get_tracer_names, &
                                 query_method, &
                                 NO_TRACER
use constants_mod,         only: CP_AIR, GRAV, HLV, HLS, HLF, &
                                 RDGAS, RVGAS, TFREEZE, WTMAIR, &
                                 SECONDS_PER_DAY, KAPPA

!  shared physics and physics utilities modules

use physics_types_mod,     only: physics_control_type
use lscloud_types_mod,     only: lscloud_types_init, atmos_state_type, &
                                 diag_id_type, diag_pt_type, &
                                 lsc_constants_type, lscloud_nml_type, &
                                 lscloud_debug_type, particles_type,  &
                                 cloud_state_type, cloud_processes_type,  &
                                 precip_state_type
use aerosol_types_mod,     only: aerosol_type
use physics_radiation_exch_mod,      &
                           only: exchange_control_type
use moist_proc_utils_mod,  only: mp_input_type, mp_output_type,  &
                                 mp_lsdiag_type, mp_nml_type,  &
                                 mp_lsdiag_control_type,  &
                                 mp_conv2ls_type, mp_tendency_type,  &
                                 mp_removal_type

! physics modules

use lscloud_debug_mod,     only: write_debug_output
use rotstayn_klein_mp_mod, only: rotstayn_klein_microp, &
                                 rotstayn_klein_microp_init,  &
                                 rotstayn_klein_microp_end
use micro_mg2_mod,         only: micro_mg2_init, micro_mg2_get_cols,&
                                 micro_mg2_tend

implicit none
private

!-----------------------------------------------------------------------
!-------------------- public data/interfaces ---------------------------

public   ls_cloud_microphysics_init, ls_cloud_microphysics, &
         ls_cloud_microphysics_end, ls_cloud_microphysics_time_vary

!-----------------------------------------------------------------------
!-------------------- private data -------------------------------------

private  adjust_precip_fields, adjust_for_supersaturation_removal,  &
         destroy_tiny_clouds


!--------------------- version number ----------------------------------
character(len=128) :: version = '$Id: $'
character(len=128) :: tagname = '$Name: $'

!--------------------------------------------------------------------------
!---namelist---------------------------------------------------------------

logical :: mass_cons = .true.         ! should we ensure water mass
                                      ! conservation by adjusting precip
                                      ! to balance column water
                                      ! mass change ?
integer :: override_liq_num = 0       ! override model predicted droplet
                                      ! number ? 0 = no, otherwise = y
integer :: override_ice_num = 0       ! override model predicted ice
                                      ! number ? 0 = no, otherwise = y
integer, dimension(6) :: init_date = (/ 1, 1, 1, 0, 0, 0 /)
                                      ! date to use as base for
                                      ! defining microphysics start time
real    :: micro_begin_sec  = 0.0     ! begin microphysics this many
                                      ! seconds after init_date
integer :: top_lev = 1                ! topmost level for ncar microphysics
real    :: min_precip_needing_adjustment     = 0.0      
real    :: lowest_allowed_precip = 0.0
logical :: use_ndust = .false.
real    :: accretion_scale = 1.0
real    :: liq_num_eros_fac = 1.0
real    :: ice_num_eros_fac = 1.0
real    :: ls_cond_max      = 0.99   ! Tiedtke sometimes converts all vapor to condensate and leads to negative vapor.
                                     ! ls_cond_max determins the maximum vapor for large-scale condensation/deposition.
logical :: limit_ls_cond = .false.   ! limit large-scale condensation/deposition from Tiedtke to avoid negative vapor. 

logical :: do_cleanup = .false.
logical :: debug_cld_microphysics = .false.
namelist / ls_cloud_microphysics_nml /   &
                               mass_cons, &
                               override_liq_num, override_ice_num, &
                               init_date, &
                               micro_begin_sec, &
                               min_precip_needing_adjustment, &
                               lowest_allowed_precip, use_ndust, accretion_scale, &
                               do_cleanup, liq_num_eros_fac, ice_num_eros_fac,    &
                               debug_cld_microphysics, ls_cond_max, limit_ls_cond   !h1g, 2020-06-22

!-------------------- clock definitions --------------------------------

integer  :: rk_micro_clock, ncar_micro_clock

!----------------------------------------------------------------------
!    module variables retrieved from other modules
!----------------------------------------------------------------------
real    :: qmin
logical :: limit_conv_cloud_frac
integer :: super_ice_opt
logical :: do_pdf_clouds
logical :: doing_prog_clouds
real    :: dtcloud, inv_dtcloud
logical :: do_rk_microphys, do_ncar_MG2
logical :: tiedtke_macrophysics
logical :: dqa_activation, total_activation
integer :: nsphum, nql, nqi, nqa, nqn, nqni, nqr, nqs, nqg, nqnr, nqns

!--------------------------------------------------------------------
!    other module variables
!--------------------------------------------------------------------

integer, parameter          :: r8 = selected_real_kind(12)
                                  ! 8 byte real
integer                     :: current_days0, current_sec0
                                  ! variables related to delayed initiation
                                  ! of microphysics



logical            :: module_is_initialized = .false.




                             contains



!#######################################################################

subroutine ls_cloud_microphysics_init  (    &
                     Nml_mp, Constants_lsc, Physics_control, &
                     id, jd, kd, Time, axes, Nml_lsc, Exch_ctrl)

!------------------------------------------------------------------------

type(mp_nml_type),           intent(in)    :: Nml_mp
type(lsc_constants_type),    intent(in)    :: Constants_lsc
type(physics_control_type),  intent(in)    :: Physics_control
integer,                     intent(in)    :: id, jd, kd
integer,                     intent(in)    :: axes(4)
type(time_type),             intent(in)    :: Time
type(lscloud_nml_type),      intent(in)    :: Nml_lsc
type(exchange_control_type), intent(in)    :: Exch_ctrl

!------------------------------------------------------------------------
! local variables:
      integer              :: logunit, io, ierr
      type(time_type)      :: Time_init
      character(len=128)   :: errstring ! Output status: non-blank for
                                        ! error return
      integer              :: rk_micro_init_clock, &
                              ncar_micro_init_clock

!-----------------------------------------------------------------------

      if (module_is_initialized) return

!-----------------------------------------------------------------------
!    save variables needed from other modules as module variables
!-----------------------------------------------------------------------
      qmin = Exch_ctrl%qmin
      limit_conv_cloud_frac = Nml_mp%limit_conv_cloud_frac
      super_ice_opt = Nml_lsc%super_ice_opt
      do_pdf_clouds = Nml_lsc%do_pdf_clouds
      doing_prog_clouds = Exch_ctrl%doing_prog_clouds
      do_rk_microphys = Constants_lsc%do_rk_microphys
      do_ncar_MG2       = Constants_lsc%do_ncar_MG2

      tiedtke_macrophysics = Constants_lsc%tiedtke_macrophysics
      dqa_activation = Constants_lsc%dqa_activation
      total_activation = Constants_lsc%total_activation
      nql = Physics_control%nql
      nqi = Physics_control%nqi
      nqa = Physics_control%nqa
      nqn = Physics_control%nqn
      nqni = Physics_control%nqni
      nqg = Physics_control%nqg
      nqr = Physics_control%nqr
      nqs = Physics_control%nqs

      nqnr = Physics_control%nqnr
      nqns = Physics_control%nqns

!------------------------------------------------------------------------
!    define clocks for large-scale cloud schemes initialization (local
!    variables) and for the prognostic loop clocks (module variables).
!------------------------------------------------------------------------
      if (do_rk_microphys) then
        rk_micro_init_clock = mpp_clock_id(   &
               '   Ls_cld_micro: rk_micro:Initialization' , &
                                                grain=CLOCK_MODULE_DRIVER )
      else if (do_ncar_MG2) then
        ncar_micro_init_clock = mpp_clock_id(    &
               '   Ls_cld_micro: ncar_MG2:Initialization' , &
                                                grain=CLOCK_MODULE_DRIVER)
      endif

      if (do_rk_microphys) then
        rk_micro_clock = mpp_clock_id(   &
               '   Ls_cld_micro: rk_micro' , &
                                                grain=CLOCK_MODULE_DRIVER )
      else if (do_ncar_MG2) then
        ncar_micro_clock = mpp_clock_id(    &
               '   Ls_cld_micro: ncar_MG2' , &
                                                grain=CLOCK_MODULE_DRIVER )
      endif

!-------------------------------------------------------------------------
!    process namelist.
!-------------------------------------------------------------------------
      read (input_nml_file, nml=ls_cloud_microphysics_nml, iostat=io)
      ierr = check_nml_error(io,'ls_cloud_microphysics_nml')

!-------------------------------------------------------------------------
!    write version number and namelist to standard log.
!-------------------------------------------------------------------------
      call write_version_number (version, tagname)
      logunit = stdlog()
      if (mpp_pe() == mpp_root_pe()) &
                          write (logunit, nml=ls_cloud_microphysics_nml)

!-------------------------------------------------------------------------
!    make sure needed modules have been initialized.
!-------------------------------------------------------------------------
      call lscloud_types_init

!-----------------------------------------------------------------------
!    initialize the active microphysics scheme module.
!-----------------------------------------------------------------------
      if (doing_prog_clouds) then

!-----------------------------------------------------------------------
!  rotstayn-klein microphysics
!-----------------------------------------------------------------------
        if (do_rk_microphys) then
          call mpp_clock_begin (rk_micro_init_clock)
          call rotstayn_klein_microp_init (Nml_lsc, Exch_ctrl)
          call mpp_clock_end   (rk_micro_init_clock)

!-----------------------------------------------------------------------
!  ncar microphysics  (ncar v2.0)
!-----------------------------------------------------------------------
        else if (do_ncar_MG2) then
          call mpp_clock_begin (ncar_micro_init_clock)
          call micro_mg2_init (r8, GRAV,  RDGAS, RVGAS, CP_AIR, TFREEZE, &
                               HLV, HLF, Nml_lsc%do_ice_nucl_wpdf,   &
                              errstring, Exch_ctrl)
          if (trim(errstring) /= '') then
            call error_mesg ('ls_cloud_microphysics/micro_mg2_init', &
                                                         errstring, FATAL)
          endif
          call mpp_clock_end (ncar_micro_init_clock)

!-----------------------------------------------------------------------
!  no valid microphys scheme chosen
!-----------------------------------------------------------------------
        else
           call error_mesg   &
              ('ls_cloud_microphysics/ls_cloud_microphysics_init', &
               'invalid microphys_scheme option in lscloud_driver nml',  &
                                                                     FATAL)
        endif
      endif  ! (doing_prog_clouds)

!-------------------------------------------------------------------------
!    get namelist initial time from namelist to determine whether
!    it is time for microphysics to be active.
!-------------------------------------------------------------------------
      Time_init = set_date( init_date(1), init_date(2), init_date(3),  &
                            init_date(4), init_date(5), init_date(6) )
      call get_time( Time_init, current_sec0, current_days0)

!------------------------------------------------------------------------
      module_is_initialized = .true.

!------------------------------------------------------------------------


end subroutine ls_cloud_microphysics_init





!########################################################################

subroutine ls_cloud_microphysics_time_vary (dtcloud_in)

real, intent(in) :: dtcloud_in

!-----------------------------------------------------------------------
!    define current timestep and its inverse.
!-----------------------------------------------------------------------
      dtcloud = dtcloud_in
      inv_dtcloud = 1.0/dtcloud

!----------------------------------------------------------------------

end subroutine ls_cloud_microphysics_time_vary


!########################################################################

subroutine ls_cloud_microphysics (    &
                 is, ie, js, je, Time, dt, lon, lat, Input_mp, Output_mp, C2ls_mp,&
                 Tend_mp, Lsdiag_mp, Lsdiag_mp_control, Atmos_state,   &
                 Cloud_state, Particles, Precip_state, Cloud_processes, &
                                                      Removal_mp, Aerosol)

!-----------------------------------------------------------------------

type(mp_input_type),        intent(inout)        :: Input_mp
type(mp_output_type),       intent(inout)        :: Output_mp
type(mp_removal_type),      intent(inout)        :: Removal_mp
type(mp_conv2ls_type),      intent(inout)        :: C2ls_mp
type(mp_tendency_type),     intent(inout)        :: Tend_mp
type(mp_lsdiag_type),       intent(inout)        :: Lsdiag_mp
type(mp_lsdiag_control_type), intent(inout)      :: Lsdiag_mp_control
type(atmos_state_type),     intent(inout)        :: Atmos_state
type(cloud_state_type),     intent(inout)        :: Cloud_state
type(particles_type),       intent(inout)        :: Particles
type(precip_state_type),    intent(inout)        :: Precip_state
type(cloud_processes_type), intent(inout)        :: Cloud_processes
type(time_type),            intent(in)           :: Time
integer,                    intent(in)           :: is, ie, js, je
real,                       intent(in)           :: dt
type(aerosol_type),         intent(in), optional :: Aerosol
real,                        intent(in), dimension(:,:) :: lon, lat


!------------------------------------------------------------------------
!   local variables:

      real, dimension(size(Input_mp%tin,1),    &
                            size(Input_mp%tin,2)) ::    &
                               ice_lin, graupel_lin,   &
                               enth_micro_col,  wat_micro_col

      real, dimension( size(Input_mp%tin,1), size(Input_mp%tin,2),   &
                            size(Input_mp%tin,3))   ::   &
                               delp, delz, &
                               ST_micro, SQ_micro, SL_micro, SI_micro, &
                               SN_micro, SNI_micro,                    &
                               SR_micro, SS_micro, SNR_micro,SNS_micro,&
                               D_eros_l, D_eros_i,  &
                               nerosc, nerosi, dqcdt, dqidt, qa_new, &
                               ssat_disposal, ql_new,  qi_new,           &
                               nctend, nitend, qn_new, qni_new, &
                               rho, liqcldf, icecldf, tmp2s,  &
                               accre_enhann, tnd_qsnown, &
                               tnd_nsnown, re_icen, relvarn, &
                               crystal1, rho_air,  &
                               aerosols_concen, droplets_concen, test_bqx, dte3d

      real, dimension( size(Input_mp%tin,1), size(Input_mp%tin,2),   &
                            size(Input_mp%tin,3),4) ::  &
                               rbar_dust_4bin, ndust_4bin

      integer,dimension(:),allocatable  :: mgcols

      integer               :: mgncol
      integer               :: i, j, k, n
      integer               :: ix, jx, kx
      integer               :: nlev
      character(len=128)    :: errstring
      real                  :: current_total_sec
      integer               :: current_sec, current_days
      real                  :: depth

!-------------------------------------------------------------------------
!   define array dimensions
!-------------------------------------------------------------------------
      ix = size(Input_mp%tin,1)
      jx = size(Input_mp%tin,2)
      kx = size(Input_mp%tin,3)

!------------------------------------------------------------------------
!   call selected microphysics scheme.
!------------------------------------------------------------------------
!------------------------------------------------------------------------
!    Rotstayn-Klein microphysics
!------------------------------------------------------------------------
      if (do_rk_microphys) then
        call mpp_clock_begin (rk_micro_clock)
        call rotstayn_klein_microp ( &
                      ix, jx, kx, Particles%N3D, total_activation,  &
                      dtcloud, inv_dtcloud, Input_mp%pfull,&
                      Input_mp%pmass, Atmos_state%airdens,     &
                      Atmos_state%esat0, Cloud_state%ql_in,  &
                      Cloud_state%qi_in, Cloud_state%qa_in,   &
                      Cloud_state%ql_mean, Cloud_state%qa_mean, &
                      Cloud_state%qn_mean, Input_mp%omega,  &
                      Input_mp%tin, Atmos_state%U_ca, &
                      Input_mp%qin, Atmos_state%qs,  &
                      Cloud_processes%D_eros, Cloud_processes%dcond_ls, &
                      Cloud_processes% dcond_ls_ice,         &
                      Cloud_processes%qvg, Atmos_state%gamma,   &
                      Cloud_processes%delta_cf, Particles%drop1,    &
                      Particles%concen_dust_sub, Cloud_state%ql_upd,   &
                      Cloud_state%qi_upd, Cloud_state%qn_upd,       &
                      Cloud_state%qi_mean, Cloud_state%qa_upd,   &
                      C2ls_mp%convective_humidity_area,   &
                      Lsdiag_mp_control%n_diag_4d, Lsdiag_mp%diag_4d,   &
                      Lsdiag_mp_control%diag_id,     &
                      Lsdiag_mp_control%diag_pt,  &
                      Lsdiag_mp_control%n_diag_4d_kp1,   &
                      Lsdiag_mp%diag_4d_kp1,  &
                      limit_conv_cloud_frac, &
                      Cloud_state%SA_out, Cloud_state%SN_out,        &
                      Tend_mp%ttnd, Tend_mp%qtnd, Cloud_state%SL_out,  &
                      Cloud_state%SI_out, Removal_mp%rain3d,  &
                      Removal_mp%snow3d, Removal_mp%snowclr3d,    &
                      Precip_state%surfrain, Precip_state%surfsnow,  &
                      Cloud_processes%f_snow_berg )

!-----------------------------------------------------------------------
!   define the output tendency fields.
!-----------------------------------------------------------------------
        Tend_mp%q_tnd(:,:,:,nql) = Cloud_state%SL_out
        Tend_mp%q_tnd(:,:,:,nqi) = Cloud_state%SI_out
        Tend_mp%q_tnd(:,:,:,nqa) = Cloud_state%SA_out
        if (nqn /= NO_TRACER) &
          Tend_mp%q_tnd(:,:,:,nqn) = Cloud_state%SN_out(:,:,:)
        if (nqni /= NO_TRACER) &
          Tend_mp%q_tnd(:,:,:,nqni) = Cloud_state%SNI_out(:,:,:)
        call mpp_clock_end   (rk_micro_clock)

!-----------------------------------------------------------------------
!  NCAR microphysics
!-----------------------------------------------------------------------
      else if ( do_ncar_MG2 ) then
        call mpp_clock_begin (ncar_micro_clock)

        ST_micro(:,:,:)  = 0.0
        SQ_micro(:,:,:)  = 0.0
        SL_micro(:,:,:)  = 0.0
        SI_micro(:,:,:)  = 0.0
        SN_micro(:,:,:)  = 0.0
        SNI_micro(:,:,:) = 0.0

        SR_micro(:,:,:)  = 0.0
        SNR_micro(:,:,:) = 0.0
        SS_micro(:,:,:)  = 0.0
        SNS_micro(:,:,:) = 0.0

!-----------------------------------------------------------------------
!   determine whether NCAR microphysics are active at the current time
!-----------------------------------------------------------------------
          Atmos_state%tn = Input_mp%tin + Tend_mp%ttnd
          Atmos_state%qvn = Input_mp%qin + Tend_mp%qtnd

!--------------------------------------------------------------------------
!     define some input fields related to ls condensation and the  cloud
!     erosion process that are needed when tiedtke macrophysics are active,
!     since the magnitude of these processes is still subject to change
!     based on what the microphysics does. for the non-tiedtke case, the
!     magnitude of these processes have been locked in before the micro-
!     physics tendencies are calculated, and so these input fields are
!     set to 0.0.
!--------------------------------------------------------------------------
          do k=1,kx
            do j=1,jx
              do i=1,ix
                Cloud_processes%dcond_ls_tot(i,j,k) =   &
                               Cloud_processes%dcond_ls(i,j,k) +   &
                                      Cloud_processes%dcond_ls_ice(i,j,k)
                if (tiedtke_macrophysics) then
                  D_eros_i(i,j,k) = -Cloud_state%qi_upd(i,j,k)* &
                                        Cloud_processes%D_eros(i,j,k)/ &
                                                                  dtcloud
                  D_eros_l(i,j,k) = -Cloud_state%ql_upd(i,j,k)* &
                                        Cloud_processes%D_eros(i,j,k)/ &
                                                                  dtcloud
                  if (Cloud_state%ql_upd(i,j,k) >= qmin) then
                    nerosc(i,j,k) = liq_num_eros_fac * D_eros_l(i,j,k)/  &
                                      Cloud_state%ql_upd(i,j,k)* &
                               Cloud_state%qn_upd(i,j,k)/MAX(0.0001, &
                                               Cloud_state%qa_upd(i,j,k))
                  else
                    nerosc(i,j,k) = 0.
                  endif
                  if (Cloud_state%qi_upd(i,j,k) >= qmin) then
                    nerosi(i,j,k) = ice_num_eros_fac * D_eros_i(i,j,k)/   &
                                          Cloud_state%qi_upd(i,j,k)* &
                               Cloud_state%qni_upd(i,j,k)/MAX(0.0001, &
                                              Cloud_state%qa_upd(i,j,k))
                  else
                    nerosi(i,j,k) = 0.
                  endif
                  if (Cloud_processes%dcond_ls_tot(i,j,k) > 0.) then
                    if (Atmos_state%tn(i,j,k) <= (tfreeze - 40.) ) then
                      dqcdt (i,j,k) = 0.
                      if ( limit_ls_cond ) then
                        dqidt(i,j,k) = min( ls_cond_max*Atmos_state%qvn(i,j,k), &
                              Cloud_processes%dcond_ls_tot(i,j,k) )*inv_dtcloud
                      else
                        dqidt(i,j,k) = Cloud_processes%dcond_ls_tot(i,j,k)*inv_dtcloud
                      endif
                    else
                      dqidt (i,j,k) = 0.
                      if ( limit_ls_cond ) then
                        dqcdt(i,j,k) = min( ls_cond_max*Atmos_state%qvn(i,j,k), &
                              Cloud_processes%dcond_ls_tot(i,j,k) )*inv_dtcloud
                      else
                        dqcdt(i,j,k) = Cloud_processes%dcond_ls_tot(i,j,k)*inv_dtcloud
                      endif
                    endif
                  else
                    if (Atmos_state%tn(i,j,k) <= tfreeze) then
                      dqcdt(i,j,k) = MAX(  &
                                     Cloud_processes%dcond_ls_tot(i,j,k),&
                                          -Cloud_state%ql_upd(i,j,k))
                      dqidt(i,j,k) = MAX(   &
                                     Cloud_processes%dcond_ls_tot(i,j,k) -&
                                                        dqcdt(i,j,k),   &
                                                -Cloud_state%qi_upd(i,j,k))
                      dqcdt(i,j,k) = dqcdt(i,j,k)*inv_dtcloud
                      dqidt(i,j,k) = dqidt(i,j,k)*inv_dtcloud
                    else
                      dqidt(i,j,k) = 0.
                      dqcdt(i,j,k) = MAX(   &
                                      Cloud_processes%dcond_ls_tot(i,j,k),&
                                           -Cloud_state%ql_upd(i,j,k))* &
                                                              inv_dtcloud
                    endif
                  endif
                else ! (tiedtke)
                  dqidt(i,j,k) = 0.
                  dqcdt(i,j,k) = 0.
                  nerosi(i,j,k) = 0.
                  nerosc(i,j,k) = 0.
                  D_eros_l(i,j,k) = 0.
                  D_eros_i(i,j,k) = 0.
                endif  ! (tiedtke)
              end do
            end do
          end do

!-------------------------------------------------------------------------
!    executed ncar microphysics:
!-------------------------------------------------------------------------
            rho = Input_mp%pfull/(RDGAS*Atmos_state%tn)

!------------------------------------------------------------------------
!   define amount of activated aerosol to be passed to microphysics.
!------------------------------------------------------------------------
!-----------------------------------------------------------------------
!   use the values previously calculated and
!   input to this routine. convert to units of #/kg.
!-----------------------------------------------------------------------
            crystal1 = Particles%crystal1/rho

!-----------------------------------------------------------------------
!   define activated droplets in units of #/kg (drop1 is #/cc).
!-----------------------------------------------------------------------
            Particles%drop2 = Particles%drop1*1.e6/Atmos_state%airdens

!------------------------------------------------------------------------
!   set liquid and ice cloud fraction to be the same as total large-scale
!   cloud fraction.
!------------------------------------------------------------------------
            liqcldf  = Cloud_state%qa_upd
            icecldf  = Cloud_state%qa_upd

!------------------------------------------------------------------------
!    are these the actual bin centers, or arbitrary values ??
!    radius of 4 dust bins for contact freezing (in the unit of m)
!------------------------------------------------------------------------
            do k = 1,kx
              do j = 1,jx
                do i = 1,ix
                  rbar_dust_4bin(i,j,k,1) = 5.e-6
                  rbar_dust_4bin(i,j,k,2) = 10.e-6
                  rbar_dust_4bin(i,j,k,3) = 15.e-6
                  rbar_dust_4bin(i,j,k,4) = 20.e-6

!------------------------------------------------------------------------
!    define the number of particles in each of the 4 dust bins, if
!    contact freezing is to be done. the active code below assigns
!    the total number to each size bin, as is done with CLUBB.
!    is this OK, or should the total number be distributed across all the
!    bins? (commented code distributes total equally across bins)
!------------------------------------------------------------------------
                  if ( use_ndust ) then
                    ndust_4bin(i,j,k,1)     = Particles%ndust(i, j, k)
                    ndust_4bin(i,j,k,2)     = Particles%ndust(i, j, k)
                    ndust_4bin(i,j,k,3)     = Particles%ndust(i, j, k)
                    ndust_4bin(i,j,k,4)     = Particles%ndust(i, j, k)
!                   ndust_4bin(i,j,k,1) = 0.25*Particles%ndust(i,j,k)
!                   ndust_4bin(i,j,k,2) = 0.25*Particles%ndust(i,j,k)
!                   ndust_4bin(i,j,k,3) = 0.25*Particles%ndust(i,j,k)
!                   ndust_4bin(i,j,k,4) = 0.25*Particles%ndust(i,j,k)

!------------------------------------------------------------------------
!    if contact freezing not desired, set ndust = 0. in each bin.
!------------------------------------------------------------------------
                  else
                    ndust_4bin(i,j,k,1)     = 0.0
                    ndust_4bin(i,j,k,2)     = 0.0
                    ndust_4bin(i,j,k,3)     = 0.0
                    ndust_4bin(i,j,k,4)     = 0.0
                  endif
                enddo
              enddo
            enddo

!------------------------------------------------------------------------
!    define the relative variance of the cloud water within each gridbox.
!    with Tiedtke macrophysics up to this time only a constant
!    value has been used, though spatial dependence could be introduced.
!------------------------------------------------------------------------
              relvarn(:,:,:) = Cloud_state%relvarn(:,:,:)

              nlev = kx
              mgncol = ix
              accre_enhann(:,:,:) = accretion_scale  ! accretion enhancement factor

              call get_time( time, current_sec, current_days)
            !  if ( mpp_pe() == mpp_root_pe() ) &
            !   write(*,*)  'current_sec =',   current_sec

            if ( debug_cld_microphysics ) then
              do k=1,kx
                do j=1,jx
                  do i=1,ix
                    if ( Atmos_state%tn(i,j,k) .lt.-150.0+273.15 .or. &
                         Atmos_state%tn(i,j,k) .gt.90+273.15)         &
            write(*,'(a,3i5, 5f12.5, 15e12.3)') 'before MG2: bad temperature@1213',   &
              i,j,k, current_sec/3600.0, Atmos_state%tn(i,j,k), Atmos_state%qvn(i,j,k), Input_mp%pfull(i,j,k), dtcloud
                  enddo
                enddo
              enddo
             endif

              do j=1,jx
               call  micro_mg2_tend (  lon(:,j), lat(:,j), &
                   dqa_activation, total_activation, &
                   tiedtke_macrophysics, j, jx,      &
                   mgncol,     nlev,     dtcloud,  &
                   Particles%concen_dust_sub(:,j,:),  &
                   Atmos_state%tn(:,j,:),     &
                   Atmos_state%qvn(:,j,:),   &
                   Cloud_state%ql_upd(:,j,:), Cloud_state%qi_upd(:,j,:), &
                   Cloud_state%qn_upd(:,j,:), Cloud_state%qni_upd(:,j,:), &
                   Cloud_state%qr_upd(:,j,:), Cloud_state%qs_upd(:,j,:), &
                   Cloud_state%qnr_upd(:,j,:), Cloud_state%qns_upd(:,j,:), &
                   relvarn(:,j,:), accre_enhann(:,j,:),                   &
                   Input_mp%pfull(:,j,:),  &
                   Atmos_state%delp(:,j,:),   &
                   Input_mp%zhalf(:,j,:),  &
                   Cloud_state%qa_upd(:,j,:),  &
                   liqcldf(:,j,:)       , icecldf(:,j,:),   &
                   Cloud_processes%delta_cf(:,j,:), &
                   D_eros_l(:,j,:), nerosc(:,j,:), &
                   D_eros_i(:,j,:), nerosi(:,j,:), &
                   dqcdt(:,j,:),    dqidt(:,j,:), &
                   crystal1(:,j,:), &
                   Particles%drop2(:,j,:),   &
                   rbar_dust_4bin(:,j,:,:), ndust_4bin(:,j,:,:), &
                   ST_micro(:,j,:), SQ_micro(:,j,:), SL_micro(:,j,:), &
                   SI_micro(:,j,:), SN_micro(:,j,:), SNI_micro(:,j,:), &
                   SR_micro(:,j,:), SS_micro(:,j,:), SNR_micro(:,j,:), SNS_micro(:,j,:),&
                   Precip_state%surfrain(:,j),   &
                   Precip_state%surfsnow(:,j),   &
                   Precip_state%lsc_snow(:,j,:), &
                   Removal_mp%rain3d(:,j,:),   &
                   Removal_mp%snow3d(:,j,:),   &
                   Precip_state%lsc_rain(:,j,:),   &
                   Precip_state%lsc_rain_size(:,j,:),  &
                   Precip_state%lsc_snow_size(:,j,:),   &
                   errstring, Cloud_processes%f_snow_berg(:,j,:), &
                   ssat_disposal (:,j,:), &
                   Lsdiag_mp_control%n_diag_4d, Lsdiag_mp%diag_4d,  &
                       Lsdiag_mp_control%diag_id,    &
                                                 Lsdiag_mp_control%diag_pt)
!Convert from effective radius to diameter for use in radiation.
! in old NCAR, diameter was returned from mmicro_pcond routine.
                   Precip_state%lsc_rain_size(:,j,:) =     &
                                   2.0*Precip_state%lsc_rain_size(:,j,:)
                   Precip_state%lsc_snow_size(:,j,:) =      &
                                  2.0*Precip_state%lsc_snow_size(:,j,:)

               if ( debug_cld_microphysics ) then
                   if(maxval(Removal_mp%snow3d(:,j,:)) > 1.e-2) write(*,*) 'max snow3d',maxval(Removal_mp%snow3d(:,j,:))
                   if(minval(Removal_mp%snow3d(:,j,:)) <-1.e-2) write(*,*) 'min snow3d',minval(Removal_mp%snow3d(:,j,:))

                do i=1,ix
                  do k=1,kx
                    if( Cloud_state%qn_upd(i,j,k) + dtcloud * SN_micro(i,j,k) < -1.e-3 ) then
                      print*, 'negative drop number @1281', lon(i,j), lat(i,j), k, Cloud_state%qn_upd(i,j,k),  &
                               SN_micro(i,j,k), Cloud_state%qn_upd(i,j,k) + dtcloud * SN_micro(i,j,k)
                    endif

                    if( Cloud_state%qr_upd(i,j,k) + dtcloud * SR_micro(i,j,k) < -1.e-3 ) then
                      print*, 'negative rain mass @1281', lon(i,j), lat(i,j), k, Cloud_state%qr_upd(i,j,k),  &
                               SR_micro(i,j,k), Cloud_state%qr_upd(i,j,k) + dtcloud * SR_micro(i,j,k)
                    endif

                    if( Cloud_state%qnr_upd(i,j,k) + dtcloud * SNR_micro(i,j,k) < -1.e-3 ) then
                      print*, 'negative rain number @1281', lon(i,j), lat(i,j), k, Cloud_state%qnr_upd(i,j,k),  &
                               SNR_micro(i,j,k), Cloud_state%qnr_upd(i,j,k) + dtcloud * SNR_micro(i,j,k)
                    endif

                    if( Cloud_state%qs_upd(i,j,k) + dtcloud * SS_micro(i,j,k) < -1.e-3 ) then
                      print*, 'negative snow mass @1281', lon(i,j), lat(i,j), k, Cloud_state%qs_upd(i,j,k),  &
                               SS_micro(i,j,k), Cloud_state%qs_upd(i,j,k) + dtcloud * SS_micro(i,j,k)
                    endif

                    if( Cloud_state%qns_upd(i,j,k) + dtcloud * SNS_micro(i,j,k) < -1.e-3 ) then
                      print*, 'negative snow number @1281', lon(i,j), lat(i,j), k, Cloud_state%qns_upd(i,j,k),  &
                               SNS_micro(i,j,k), Cloud_state%qns_upd(i,j,k) + dtcloud * SNS_micro(i,j,k)
                    endif
                  enddo
                enddo
               endif

              enddo  ! end of j loop

!------------------------------------------------------------------------
!    calculate column enthalpy and total water changes
!    Note: in MG2, temperature tendency is multiplied by Cp_air.
!------------------------------------------------------------------------
              enth_micro_col(:,:) = 0.0
              wat_micro_col(:,:)  = 0.0
              do j=1,jx
                do i=1,ix
                  do k=1,kx
                    enth_micro_col(i,j) = enth_micro_col(i,j)   +         &
                        ( ST_micro(i,j,k) - HLV*SL_micro(i,j,k) - HLV*SR_micro(i,j,k)  &
                                          - HLS*SI_micro(i,j,k) - HLS*SS_micro(i,j,k) )*    &
                                             Atmos_state%delp(i,j,k)/grav

                    wat_micro_col(i,j) = wat_micro_col(i,j)  +            &
                         ( SQ_micro(i,j,k) + SL_micro(i,j,k) +  SR_micro(i,j,k) &
                                           + SI_micro(i,j,k) +  SS_micro(i,j,k) )*   &
                                             Atmos_state%delp(i,j,k)/grav
                  enddo

                  enth_micro_col(i,j) = enth_micro_col(i,j) +   &
                                                 (-HLV*1000.0* &
                                        (Precip_state%surfrain(i,j) -  &
                                         Precip_state%surfsnow(i,j)) -  &
                                 HLS*1000.0 * Precip_state%surfsnow(i,j) )

                  wat_micro_col(i,j) = wat_micro_col(i,j) +   &
                                       Precip_state%surfrain(i,j) *1000.0

                   if ( debug_cld_microphysics ) then
                    if(  abs(enth_micro_col(i,j)) > 1.e-8 ) then
                      print*, 'enth-imb after MG2@1236', lon(i,j), lat(i,j), enth_micro_col(i,j)  
                    endif
                    if(  abs(wat_micro_col(i,j))  > 1.e-10  ) then
                      print*, 'wat-imb after MG2@1239', lon(i,j), lat(i,j), wat_micro_col(i,j)  
                    endif
                   endif
                              
                enddo
              enddo

!------------------------------------------------------------------------
!    adjust precip fields to assure mass conservation and realizable
!    values.
!------------------------------------------------------------------------
          call adjust_precip_fields (   &
                              ix, jx, kx, SQ_micro, SL_micro, SI_micro,  SR_micro, SS_micro, &
                                  Atmos_state, Precip_state, Lsdiag_mp, &
                                                       Lsdiag_mp_control )

!-----------------------------------------------------------------------
!    update prognostic tendencies due to microphysics terms.
!-----------------------------------------------------------------------
          Tend_mp%qtnd = Tend_mp%qtnd + SQ_micro*dtcloud
          Tend_mp%ttnd = Tend_mp%ttnd + ST_micro/cp_air*dtcloud
          Cloud_state%SL_out = Cloud_state%SL_out + SL_micro*dtcloud
          Cloud_state%SI_out = Cloud_state%SI_out + SI_micro*dtcloud
          Cloud_state%SN_out = Cloud_state%SN_out + SN_micro*dtcloud
          Cloud_state%SNI_out = Cloud_state%SNI_out + SNI_micro*dtcloud

          if (nqr /= NO_TRACER) then
            Cloud_state%SR_out = Cloud_state%SR_out + SR_micro*dtcloud
          endif
          if (nqs /= NO_TRACER) then
            Cloud_state%SS_out = Cloud_state%SS_out + SS_micro*dtcloud
          endif
          if (nqnr /= NO_TRACER) then
            Cloud_state%SNR_out = Cloud_state%SNR_out + SNR_micro*dtcloud
          endif
          if (nqns /= NO_TRACER) then
            Cloud_state%SNS_out = Cloud_state%SNS_out + SNS_micro*dtcloud
          endif

!------------------------------------------------------------------------
!    adjustment to fields needed after removing supersaturation
!------------------------------------------------------------------------
          call adjust_for_supersaturation_removal (  &
              ix, jx, kx, C2ls_mp, Input_mp, Atmos_state,  &
              ssat_disposal, Particles, Cloud_state, Lsdiag_mp,&
              Lsdiag_mp_control )

!------------------------------------------------------------------------
!    process fields after microphysics is completed 
!------------------------------------------------------------------------

!-------------------------------------------------------------------------
!    remove any clouds with less condensate present than the specified
!    allowable minimum.
!-------------------------------------------------------------------------
            call destroy_tiny_clouds (    &
                   ix, jx, kx, Cloud_state, Tend_mp, Lsdiag_mp, &
                      Lsdiag_mp_control, C2ls_mp, Input_mp, Atmos_state)

!-----------------------------------------------------------------------
!    define output fields.
!-----------------------------------------------------------------------
            Tend_mp%q_tnd(:,:,:,nql) = Cloud_state%SL_out
            Tend_mp%q_tnd(:,:,:,nqi) = Cloud_state%SI_out
            Tend_mp%q_tnd(:,:,:,nqa) = Cloud_state%SA_out
            if (nqn /= NO_TRACER) &
                     Tend_mp%q_tnd(:,:,:,nqn) = Cloud_state%SN_out(:,:,:)
            if (nqni /= NO_TRACER) &
                     Tend_mp%q_tnd(:,:,:,nqni) = Cloud_state%SNI_out(:,:,:)
            if (nqr /= NO_TRACER) &
                     Tend_mp%q_tnd(:,:,:,nqr)  = Cloud_state%SR_out(:,:,:)
            if (nqs /= NO_TRACER) &
                     Tend_mp%q_tnd(:,:,:,nqs)  = Cloud_state%SS_out(:,:,:)
            if (nqnr /= NO_TRACER) &
                     Tend_mp%q_tnd(:,:,:,nqnr) = Cloud_state%SNR_out(:,:,:)
            if (nqns /= NO_TRACER) &
                     Tend_mp%q_tnd(:,:,:,nqns) = Cloud_state%SNS_out(:,:,:)

!-----------------------------------------------------------------------
!    define the total precipitating ice field for use in COSP (stored in
!    Removal_mp%snowclr3d).it is moved to this field so that total cloud
!    ice plus snow is contained in (Removal_mp%snowclr3d + the cloud ice
!    field), as in the R-K microphysics case.
!-----------------------------------------------------------------------
            Removal_mp%snowclr3d = Removal_mp%snow3d

        call mpp_clock_end   (ncar_micro_clock)
!-------------------------------------------------------------------------
!    exit with error if no valid microphysics scheme was specified.
!-------------------------------------------------------------------------
      else    ! do rk
        call error_mesg ('ls_cloud_microphysics/ls_cloud_microphysics', &
              'invalid lscloud_driver_nml microphys_scheme option', FATAL)
      endif    ! (do_rk)

!---------------------------------------------------------------------


end subroutine ls_cloud_microphysics



!#######################################################################

subroutine ls_cloud_microphysics_end

!------------------------------------------------------------------------

      integer   :: rk_micro_term_clock, &
                   ncar_micro_term_clock

!------------------------------------------------------------------------

      if (.not. module_is_initialized) return

!------------------------------------------------------------------------
!    define clocks for each microphysics scheme.
!------------------------------------------------------------------------
      if (do_rk_microphys) then
        rk_micro_term_clock = mpp_clock_id(   &
               '   Ls_cld_micro: rk_micro:Termination' , &
                                                grain=CLOCK_MODULE_DRIVER )
      endif

!-------------------------------------------------------------------------
!    call the termination routines for each of the available
!    microphysics schemes which has one.
!-------------------------------------------------------------------------
      if (do_rk_microphys ) then
        call mpp_clock_begin (rk_micro_term_clock)
        call rotstayn_klein_microp_end
        call mpp_clock_end   (rk_micro_term_clock)
      endif

      module_is_initialized = .false.

!----------------------------------------------------------------------

end subroutine ls_cloud_microphysics_end


!########################################################################

subroutine adjust_precip_fields (    &
              ix, jx, kx, SQ_micro, SL_micro, SI_micro,  SR_micro, SS_micro,   &
                                   Atmos_state, Precip_state, Lsdiag_mp, &
                                                 Lsdiag_mp_control )

!------------------------------------------------------------------------
!    subroutine adjust_precip_fields modifies the surface precipitation to
!    balance the atmospheric tendencies of water, and thus conserve water
!    mass. Any needed adjustments are available for examination as
!    model netcdf diagnostics.
!------------------------------------------------------------------------

integer,                    intent(in)    :: ix, jx, kx
type(mp_lsdiag_type),       intent(inout) :: Lsdiag_mp
type(mp_lsdiag_control_type), intent(inout) :: Lsdiag_mp_control
type(atmos_state_type),     intent(inout) :: Atmos_state
type(precip_state_type),    intent(inout) :: Precip_state
real, dimension (:,:,:),    intent(in)    :: SL_micro, SI_micro, SQ_micro, SR_micro, SS_micro


!----------------------------------------------------------------------
!   local variables:

      real, dimension (ix, jx)     :: m1, m2, scalef
      integer                      :: i, j, k

!----------------------------------------------------------------------
!    if enforcement of water mass conservation is desired, compute
!    the precip reaching the ground and the net change in vapor,
!    cloud water and cloud ice.
!----------------------------------------------------------------------
      if (mass_cons) then
        do j=1,jx
          do i=1,ix
            m1(i,j) = 0.
            do k=1,kx
              m1(i,j) = m1(i,j) +   &
                   (SQ_micro(i,j,k) + SL_micro(i,j,k) + SI_micro(i,j,k) + SR_micro(i,j,k) + SS_micro(i,j,k) )* &
                                    dtcloud*Atmos_state%delp(i,j,k)/grav
            end do
            m2(i,j) = 1.e3*Precip_state%surfrain(i,j)*dtcloud
!------------------------------------------------------------------------
!    for small precip, adjustment for conservation may be ignored. other-
!    wise, compute the ratio of condensate loss to precip at the surface.
!------------------------------------------------------------------------
            if ( m2(i,j) .ne. 0.0 ) THEN
              scalef(i,j) = -m1(i,j)/m2(i,j)
!-----------------------------------------------------------------------
!   define diagnostics capturing the rate (kg/m2/s) that the precip
!   field is adjusted to balance the loss of atmospheric water mass.
!-----------------------------------------------------------------------
              if (Lsdiag_mp_control%diag_id%rain_mass_conv > 0   ) &
              Lsdiag_mp%diag_4d(i,j,1,   &
                            Lsdiag_mp_control%diag_pt%rain_mass_conv) = &
                           (scalef(i,j)*Precip_state%surfrain(i,j) -    &
                                        Precip_state%surfrain(i,j))*1.0e3
              if (Lsdiag_mp_control%diag_id%snow_mass_conv > 0   ) &
              Lsdiag_mp%diag_4d(i,j,1,   &
                           Lsdiag_mp_control%diag_pt%snow_mass_conv) = &
                              (scalef(i,j)*Precip_state%surfsnow(i,j) -  &
                                         Precip_state%surfsnow(i,j))*1.0e3
!------------------------------------------------------------------------
!    modify the output rain and snow precip fields.
!------------------------------------------------------------------------
              Precip_state%surfrain(i,j) =    &
                                   scalef(i,j)*Precip_state%surfrain(i,j)
              Precip_state%surfsnow(i,j) =    &
                                   scalef(i,j)*Precip_state%surfsnow(i,j)
            end if
          end do
        end do
      end if
      
!------------------------------------------------------------------------
!    save the rain and snow precipitation fields before any lower limit
!    is imposed (usually 0.0).
!------------------------------------------------------------------------
      if (Lsdiag_mp_control%diag_id%neg_rain > 0) &
        Lsdiag_mp%diag_4d(:,:,1,    &
                      Lsdiag_mp_control%diag_pt%neg_rain) = 1.0e3*    &
          (Precip_state%surfrain(:,:) - Precip_state%surfsnow(:,:))* &
                                                                   dtcloud
      if (Lsdiag_mp_control%diag_id%neg_snow > 0) &
        Lsdiag_mp%diag_4d(:,:,1,    &
                        Lsdiag_mp_control%diag_pt%neg_snow) = 1.0e3*    &
          (Precip_state%surfsnow(:,:))*dtcloud

!-----------------------------------------------------------------------
!    impose lower limit.
!-----------------------------------------------------------------------
      Precip_state%surfrain = max(     &
             1.e3*(Precip_state%surfrain - Precip_state%surfsnow)*   &
                                        dtcloud , lowest_allowed_precip)
      Precip_state%surfsnow = max(    &
             1.e3*Precip_state%surfsnow*dtcloud,   &
                                                  lowest_allowed_precip)

!-----------------------------------------------------------------------
!    compute amount of precip which has been eliminated by this
!    adjustment.
!-----------------------------------------------------------------------
      if (Lsdiag_mp_control%diag_id%neg_rain > 0) &
          Lsdiag_mp%diag_4d(:,:,1,Lsdiag_mp_control%diag_pt%neg_rain) =   &
          -1.0*( (Precip_state%surfrain(:,:))  -   &
                  Lsdiag_mp%diag_4d(:,:,1,   &
                              Lsdiag_mp_control%diag_pt%neg_rain))/dtcloud
      if (Lsdiag_mp_control%diag_id%neg_snow > 0) &
          Lsdiag_mp%diag_4d(:,:,1,Lsdiag_mp_control%diag_pt%neg_snow) =   &
          -1.0*( (Precip_state%surfsnow(:,:))  -   &
                  Lsdiag_mp%diag_4d(:,:,1,    &
                             Lsdiag_mp_control%diag_pt%neg_snow))/dtcloud

!-------------------------------------------------------------------------

end subroutine adjust_precip_fields

!########################################################################

subroutine adjust_for_supersaturation_removal (  &
                      ix, jx, kx, C2ls_mp, Input_mp, Atmos_state, &
                       ssat_disposal, Particles, Cloud_state, Lsdiag_mp, &
                                                       lsdiag_mp_control )

!-----------------------------------------------------------------------
!    with tiedtke macrophysics, supersaturation removal results in an
!    increase in cloudiness to the max allowable cloudiness in the grid
!    box and a consequent increase in activated aerosols due to this
!    increase in coverage when the Ming dqa activation is being used.
!-----------------------------------------------------------------------

integer,                    intent(in)    :: ix, jx, kx
type(mp_lsdiag_type),       intent(inout) :: Lsdiag_mp
type(mp_lsdiag_control_type), intent(inout) :: Lsdiag_mp_control
type(atmos_state_type),     intent(inout) :: Atmos_state
type(mp_input_type),        intent(inout) :: Input_mp
type(mp_conv2ls_type),      intent(inout) :: C2ls_mp
type(cloud_state_type),     intent(inout) :: Cloud_state
type(particles_type),       intent(inout) :: Particles
real, dimension(:,:,:),     intent(in)    :: ssat_disposal


!-----------------------------------------------------------------------
!   local variables:

      real, dimension(ix, jx, kx)  :: rho, tmp2s
      integer                      :: i, j, k

!----------------------------------------------------------------------
!    process only occurs if using tiedtke macrophysics and not doing
!    pdf clouds.
!----------------------------------------------------------------------
      if (tiedtke_macrophysics .and.  .not. do_pdf_clouds) then

!-----------------------------------------------------------------------
!    where supersaturation is present, define the effects of removing it
!    on the cloud area and cloud particle / ice crystal number.
!-----------------------------------------------------------------------
        do k=1,kx
          do j=1,jx
            do i=1,ix
              if (ssat_disposal(i,j,k) > 0.0) then

!-----------------------------------------------------------------------
!    define the density (rho).
!-----------------------------------------------------------------------
                rho(i,j,k) = Input_mp%pfull(i,j,k)/   &
                                        (RDGAS*Atmos_state%tn(i,j,k))

!-----------------------------------------------------------------------
!    define the area unavailable for large-scale clouds due to it
!    containing convective cloud (tmp2s).
!-----------------------------------------------------------------------
                if (limit_conv_cloud_frac) then
                  tmp2s(i,j,k) = C2ls_mp%convective_humidity_area(i,j,k)
                else
                  tmp2s(i,j,k) = 0.
                endif

!-----------------------------------------------------------------------
!    when dqa activation is being used, the increase in cloud area results
!    in an increase in activated ice particles and cloud nuclei,
!    proportional to the cloud area increase. save the incremental
!    increase due to removing superstauration as diagnostics.
!-----------------------------------------------------------------------
             !   if (dqa_activation) then    ! h1g, 2020-03-19
                  if (ssat_disposal(i,j,k) == 2.) then
                    Cloud_state%SNi_out(i,j,k) =     &
                        Cloud_state%SNi_out(i,j,k) + &
                          (Particles%crystal1(i,j,k)/rho(i,j,k)*  &
                               (1. - Cloud_state%qa_upd(i,j,k) -  &
                                   tmp2s(i,j,k))/dtcloud)*dtcloud
                    if (Lsdiag_mp_control%diag_id%qnidt_super +   &
                            Lsdiag_mp_control%diag_id%qni_super_col > 0 ) &
                        Lsdiag_mp%diag_4d(i,j,k,  &
                              Lsdiag_mp_control%diag_pt%qnidt_super) =    &
                        Particles%crystal1(i,j,k)/rho(i,j,k)*  &
                     (1. - Cloud_state%qa_upd(i,j,k) - tmp2s(i,j,k))/  &
                                                                  dtcloud

                  else if (ssat_disposal(i,j,k) == 1.) then
                    Cloud_state%SN_out(i,j,k) =    &
                        Cloud_state%SN_out(i,j,k) +  &
                          (Particles%drop2(i,j,k)*    &
                               (1. - Cloud_state%qa_upd(i,j,k) -   &
                                  tmp2s(i,j,k))/dtcloud)*dtcloud
                    if (Lsdiag_mp_control%diag_id%qndt_super +   &
                            Lsdiag_mp_control%diag_id%qn_super_col > 0 ) &
                         Lsdiag_mp%diag_4d(i,j,k,   &
                             Lsdiag_mp_control%diag_pt%qndt_super ) =    &
                           Particles%drop2(i,j,k)*               &
                     (1. - Cloud_state%qa_upd(i,j,k) - tmp2s(i,j,k))/  &
                                                                  dtcloud
                  endif
             !   end if ! dqa_activation ! h1g, 2020-03-19
                if (max(Lsdiag_mp_control%diag_id%qadt_super,  &
                            Lsdiag_mp_control%diag_id%qa_super_col) > 0) then
                  Lsdiag_mp%diag_4d(i,j,k,  &
                        Lsdiag_mp_control%diag_pt%qadt_super ) = &
                     (1. - Cloud_state%qa_upd(i,j,k) - tmp2s(i,j,k))/  &
                                                                  dtcloud
                endif

!-------------------------------------------------------------------------
!    add the change to the cloud area increment resulting from this
!    process (SA_out), and update the model cloud area after this process
!    is completed (Cloud_state%qa_upd).
!-------------------------------------------------------------------------
                Cloud_state%SA_out(i,j,k) =   &
                      Cloud_state%SA_out(i,j,k) + &
                          (1. - Cloud_state%qa_upd(i,j,k) - tmp2s(i,j,k))
                Cloud_state%qa_upd(i,j,k) = 1. - tmp2s(i,j,k)
              endif ! ssat_disposal > 0.0
            end do
          end do
        end do
      end if

!-------------------------------------------------------------------------


end subroutine adjust_for_supersaturation_removal


!######################################################################

subroutine destroy_tiny_clouds (   &
                  ix, jx, kx, Cloud_state, Tend_mp, Lsdiag_mp,   &
                        Lsdiag_mp_control, C2ls_mp, Input_mp, Atmos_state)

!-----------------------------------------------------------------------
!    routine to conservatively remove unacceptably small clouds.
!-----------------------------------------------------------------------

integer,                    intent(in)    :: ix, jx, kx
type(mp_lsdiag_type),       intent(inout) :: Lsdiag_mp
type(mp_lsdiag_control_type), intent(inout) :: Lsdiag_mp_control
type(mp_tendency_type),     intent(inout) :: Tend_mp
type(atmos_state_type),     intent(inout) :: Atmos_state
type(mp_input_type),        intent(inout) :: Input_mp
type(mp_conv2ls_type),      intent(inout) :: C2ls_mp
type(cloud_state_type),     intent(inout) :: Cloud_state

!-----------------------------------------------------------------------
!   local variables:

      real, dimension (ix,jx,kx)   :: ql_new, qi_new, qn_new, qni_new, &
                                      qa_new
      integer :: i,j,k
      real, dimension (ix,jx,kx)   :: qr_new, qnr_new, qs_new, qns_new

!-----------------------------------------------------------------------
!    define current cloud and particle values.
!----------------------------------------------------------------------
      ql_new  = Cloud_state%ql_in  + Cloud_state%SL_out
      qi_new  = Cloud_state%qi_in  + Cloud_state%SI_out
      qn_new  = Cloud_state%qn_in  + Cloud_state%SN_out
      qni_new = Cloud_state%qni_in + Cloud_state%SNi_out

      qr_new  = Cloud_state%qr_in  + Cloud_state%SR_out
      qnr_new = Cloud_state%qnr_in + Cloud_state%SNR_out
      qs_new  = Cloud_state%qs_in  + Cloud_state%SS_out
      qns_new = Cloud_state%qns_in + Cloud_state%SNS_out
!-----------------------------------------------------------------------
!    if these values are lower than acceptable, or if the new cloud area
!    is lower than acceptable, set the tendency to balance the input value,
!    so that the field is 0. upon exiting this routine. include
!    adjustments to temp and vapor to conserve energy and water mass.
!-----------------------------------------------------------------------
      do k=1,kx
        do j=1,jx
          do i=1,ix
            if ((ql_new(i,j,k) <= qmin  .and.    &
                 qi_new(i,j,k) <= qmin)   .or.   &
                (Cloud_state%qa_upd(i,j,k) <=        qmin)) then
              Cloud_state%SL_out(i,j,k) = Cloud_state%SL_out(i,j,k) -  &
                                                              ql_new(i,j,k)
              Cloud_state%SI_out(i,j,k) = Cloud_state%SI_out(i,j,k) -  &
                                                              qi_new(i,j,k)
              Cloud_state%SA_out(i,j,k) = Cloud_state%SA_out(i,j,k) -  &
                                                  Cloud_state%qa_upd(i,j,k)
              Tend_mp%ttnd(i,j,k) = Tend_mp%ttnd(i,j,k) -    &
                            (hlv*ql_new(i,j,k) + hls*qi_new(i,j,k))/cp_air
              Tend_mp%qtnd(i,j,k) = Tend_mp%qtnd(i,j,k) +    &
                                           (ql_new(i,j,k) + qi_new(i,j,k))
              Cloud_state%SN_out(i,j,k) = Cloud_state%SN_out(i,j,k) -  &
                                                             qn_new(i,j,k)
              Cloud_state%SNi_out(i,j,k) = Cloud_state%SNi_out(i,j,k) -  &
                                                            qni_new(i,j,k)

!------------------------------------------------------------------------
!    save diagnostics defining the adjustments made here to destroy the
!    clouds.
!------------------------------------------------------------------------
              if (Lsdiag_mp_control%diag_id%qldt_destr > 0 .or.   &
                             Lsdiag_mp_control%diag_id%ql_destr_col > 0) &
                  Lsdiag_mp%diag_4d(i,j,k,  &
                             Lsdiag_mp_control%diag_pt%qldt_destr) =  &
                                      - ql_new(i,j,k)/dtcloud
              if (Lsdiag_mp_control%diag_id%qidt_destr > 0 .or.   &
                             Lsdiag_mp_control%diag_id%qi_destr_col > 0) &
                Lsdiag_mp%diag_4d(i,j,k,     &
                                Lsdiag_mp_control%diag_pt%qidt_destr) =  &
                                      - qi_new(i,j,k)/dtcloud
              if (Lsdiag_mp_control%diag_id%qadt_destr > 0 .or.   &
                             Lsdiag_mp_control%diag_id%qa_destr_col > 0) &
                   Lsdiag_mp%diag_4d(i,j,k,    &
                                Lsdiag_mp_control%diag_pt%qadt_destr) =  &
                          - Cloud_state%qa_upd(i,j,k)/dtcloud
              if (Lsdiag_mp_control%diag_id%qndt_destr > 0 .or.   &
                             Lsdiag_mp_control%diag_id%qn_destr_col > 0) &
                Lsdiag_mp%diag_4d(i,j,k,   &
                              Lsdiag_mp_control%diag_pt%qndt_destr) =  &
                                      - qn_new(i,j,k)/dtcloud
              if (Lsdiag_mp_control%diag_id%qnidt_destr +    &
                            Lsdiag_mp_control%diag_id%qni_destr_col > 0) &
                Lsdiag_mp%diag_4d(i,j,k,   &
                            Lsdiag_mp_control%diag_pt%qnidt_destr) = &
                                    - qni_new(i,j,k)/dtcloud
              if (Lsdiag_mp_control%diag_id%qdt_destr +    &
                              Lsdiag_mp_control%diag_id%q_destr_col > 0) &
                Lsdiag_mp%diag_4d(i,j,k,    &
                           Lsdiag_mp_control%diag_pt%qdt_destr) =   &
                      (ql_new(i,j,k) + qi_new(i,j,k))/dtcloud
            endif
          end do
        end do
      end do

!-----------------------------------------------------------------------
!    redefine the new cloud tracer values.
!-----------------------------------------------------------------------

    if ( do_cleanup ) then  ! --> h1g, 20200317
      ql_new  =  Cloud_state%ql_in  + Cloud_state%SL_out
      qi_new  =  Cloud_state%qi_in  + Cloud_state%SI_out
      qn_new  =  Cloud_state%qn_in  + Cloud_state%SN_out
      qni_new =  Cloud_state%qni_in + Cloud_state%SNI_out

!-----------------------------------------------------------------------
!    if the new value of cloud water is too small (including negative
!    roundoff values), and the vapor will remain positive when
!    conservatively adjusted, eliminate the cloudwater by adjusting the
!    vapor.
!-----------------------------------------------------------------------
      do k=1,kx
        do j=1,jx
          do i=1,ix
            if (abs(ql_new(i,j,k)) .le. qmin  .and.  &
                 Input_mp%qin(i,j,k) + Tend_mp%qtnd(i,j,k) +   &
                                                 ql_new(i,j,k) > 0.0) then
              Cloud_state%SL_out(i,j,k) =  - Cloud_state%ql_in(i,j,k)
              Tend_mp%qtnd(i,j,k) = Tend_mp%qtnd(i,j,k) + ql_new(i,j,k)
              Tend_mp%ttnd(i,j,k) = Tend_mp%ttnd(i,j,k) -    &
                                               (hlv*ql_new(i,j,k))/cp_air

!------------------------------------------------------------------------
!    compute diagnostic for this liquid loss due to this cleanup.
!------------------------------------------------------------------------
              if (Lsdiag_mp_control%diag_id%qdt_cleanup_liquid +    &
                    Lsdiag_mp_control%diag_id%q_cleanup_liquid_col > 0) &
                 Lsdiag_mp%diag_4d(i,j,k,   &
                       Lsdiag_mp_control%diag_pt%qdt_cleanup_liquid) =   &
                                                   ql_new(i,j,k)/dtcloud

!------------------------------------------------------------------------
!    with the removal of all liquid, the cloud droplet number must also
!    be set to 0.0. define diagnostic for droplet loss due to this cleanup.
!------------------------------------------------------------------------
              Cloud_state%SN_out(i,j,k) = Cloud_state%SN_out(i,j,k) -   &
                                                           qn_new(i,j,k)
              if (Lsdiag_mp_control%diag_id%qndt_cleanup +   &
                          Lsdiag_mp_control%diag_id%qn_cleanup_col > 0) &
                Lsdiag_mp%diag_4d(i,j,k,   &
                         Lsdiag_mp_control%diag_pt%qndt_cleanup) = &
                                              - qn_new(i,j,k)/dtcloud
            endif
          end do
        end do
      end do

!-----------------------------------------------------------------------
!    if the new value of cloud ice is too small (including negative
!    roundoff values), and the vapor will remain positive when
!    conservatively adjusted, eliminate the cloudice by adjusting the
!    vapor.
!-----------------------------------------------------------------------
      do k=1,kx
        do j=1,jx
          do i=1,ix
            if (abs(qi_new(i,j,k)) .le. qmin  .and.  &
                  Input_mp%qin(i,j,k) + Tend_mp%qtnd(i,j,k) +   &
                                                qi_new(i,j,k) > 0.0) then
              Cloud_state%SI_out(i,j,k) =  - Cloud_state%qi_in(i,j,k)
              Tend_mp%qtnd(i,j,k) = Tend_mp%qtnd(i,j,k) + qi_new(i,j,k)
              Tend_mp%ttnd(i,j,k) = Tend_mp%ttnd(i,j,k) -    &
                                                (hls*qi_new(i,j,k))/cp_air

!------------------------------------------------------------------------
!    compute diagnostic for this liquid loss due to this cleanup.
!------------------------------------------------------------------------
              if (Lsdiag_mp_control%diag_id%qdt_cleanup_ice +    &
                        Lsdiag_mp_control%diag_id%q_cleanup_ice_col > 0) &
                Lsdiag_mp%diag_4d(i,j,k,   &
                         Lsdiag_mp_control%diag_pt%qdt_cleanup_ice) =   &
                                       qi_new(i,j,k)/dtcloud

!------------------------------------------------------------------------
!    with the removal of all ice, the ice particle number must also
!    be set to 0.0. define diagnostic for crystal loss due to this cleanup.
!------------------------------------------------------------------------
              Cloud_state%SNI_out(i,j,k) = Cloud_state%SNI_out(i,j,k) -  &
                                                             qni_new(i,j,k)
              if (Lsdiag_mp_control%diag_id%qnidt_cleanup +    &
                    Lsdiag_mp_control%diag_id%qni_cleanup_col > 0) &
               Lsdiag_mp%diag_4d(i,j,k,     &
                         Lsdiag_mp_control%diag_pt%qnidt_cleanup) = &
                                  - qni_new(i,j,k)/dtcloud
            endif
          end do
        end do
      end do

!-----------------------------------------------------------------------
!    force the change in ice crystal number to not be so large as to
!    eliminate more crystals than were present initially. save a diagnostic
!    if desired.
!-----------------------------------------------------------------------
      do k=1,kx
        do j=1,jx
          do i=1,ix
            if (Lsdiag_mp_control%diag_id%qnidt_cleanup2 +    &
                   Lsdiag_mp_control%diag_id%qni_cleanup2_col > 0) &
              Lsdiag_mp%diag_4d(i,j,k,    &
                              Lsdiag_mp_control%diag_pt%qnidt_cleanup2) = &
                                           Cloud_state%SNi_out(i,j,k)

            Cloud_state%SNi_out(i,j,k) = MAX(Cloud_state%SNi_out(i,j,k), &
                                            - Cloud_state%qni_in(i,j,k))
            if (Lsdiag_mp_control%diag_id%qnidt_cleanup2 +    &
                        Lsdiag_mp_control%diag_id%qni_cleanup2_col > 0) &
                  Lsdiag_mp%diag_4d(i,j,k,   &
                          Lsdiag_mp_control%diag_pt%qnidt_cleanup2) =    &
            (Lsdiag_mp%diag_4d(i,j,k,    &
                    Lsdiag_mp_control%diag_pt%qnidt_cleanup2) - &
                     Cloud_state%SNi_out(i,j,k))*inv_dtcloud
          end do
        end do
      end do


!-----------------------------------------------------------------------
!    force the change in cloud droplet number to not be so large as to
!    eliminate more droplets than were present initially. save a diagnostic
!    if desired.
!-----------------------------------------------------------------------
      do k=1,kx
        do j=1,jx
          do i=1,ix
            if (Lsdiag_mp_control%diag_id%qndt_cleanup2 +     &
                     Lsdiag_mp_control%diag_id%qn_cleanup2_col > 0) &
                   Lsdiag_mp%diag_4d(i,j,k,    &
                          Lsdiag_mp_control%diag_pt%qndt_cleanup2) = &
                                            Cloud_state%SN_out(i,j,k)
            Cloud_state%SN_out(i,j,k) = MAX(Cloud_state%SN_out(i,j,k),  &
                                              - Cloud_state%qn_in(i,j,k))
            if (Lsdiag_mp_control%diag_id%qndt_cleanup2 +    &
                          Lsdiag_mp_control%diag_id%qn_cleanup2_col > 0) &
               Lsdiag_mp%diag_4d(i,j,k,    &
                           Lsdiag_mp_control%diag_pt%qndt_cleanup2) = &
                 (Lsdiag_mp%diag_4d(i,j,k,   &
                            Lsdiag_mp_control%diag_pt%qndt_cleanup2) -   &
                     Cloud_state%SN_out(i,j,k))*inv_dtcloud
          end do
        end do
      end do

     endif ! do_cleanup  --> h1g 20200317
!----------------------------------------------------------------------
!    make sure the new cloud area is not smaller than the minimum
!    allowable. if not set the tendency so that cloud area is reduced to
!    zero after the step. save a diagnostic if desired.
!----------------------------------------------------------------------
      if (Lsdiag_mp_control%diag_id%qadt_destr +    &
                          Lsdiag_mp_control%diag_id%qa_destr_col > 0)    &
      Lsdiag_mp%diag_4d(:,:,:,Lsdiag_mp_control%diag_pt%qadt_destr) =    &
         Lsdiag_mp%diag_4d(:,:,:,Lsdiag_mp_control%diag_pt%qadt_destr) +  &
                              Cloud_state%SA_out*inv_dtcloud

      qa_new = Cloud_state%qa_in + Cloud_state%SA_out

      where ( abs(qa_new) .le. qmin )
        Cloud_state%SA_out  = -Cloud_state%qa_in
      endwhere

      if (Lsdiag_mp_control%diag_id%qadt_destr +   &
                          Lsdiag_mp_control%diag_id%qa_destr_col > 0)    &
      Lsdiag_mp%diag_4d(:,:,:,Lsdiag_mp_control%diag_pt%qadt_destr) =    &
      Lsdiag_mp%diag_4d(:,:,:,Lsdiag_mp_control%diag_pt%qadt_destr) -     &
                              Cloud_state%SA_out*inv_dtcloud

!------------------------------------------------------------------------
!    enforce constraints on the cloud area. the change must not be so
!    large as to eliminate more cloud than is present. also it must not
!    be so large as to more than fill the available area in the grid box
!    (some area may have been taken up by the convective system, so the
!    max available area is (1 - conv area). Include a diagnostic if
!    desired. this constraint has already been imposed with r-k
!    microphysics, as part of the destruction diagnostic.
!------------------------------------------------------------------------
      if ( do_ncar_MG2 ) then
        if (Lsdiag_mp_control%diag_id%qadt_limits +    &
                         Lsdiag_mp_control%diag_id%qa_limits_col > 0)    &
       Lsdiag_mp%diag_4d(:,:,:,Lsdiag_mp_control%diag_pt%qadt_limits) =   &
                                  Cloud_state%SA_out(:,:,:)

        Cloud_state%SA_out = MAX(Cloud_state%SA_out,-Cloud_state%qa_in)
        Cloud_state%SA_out = MIN(Cloud_state%SA_out,   &
                             1. - C2ls_mp%convective_humidity_area - &
                                                       Cloud_state%qa_in)
        if (Lsdiag_mp_control%diag_id%qadt_limits +   &
                    Lsdiag_mp_control%diag_id%qa_limits_col      > 0)    &
         Lsdiag_mp%diag_4d(:,:,:,Lsdiag_mp_control%diag_pt%qadt_limits) = &
                      (Cloud_state%SA_out(:,:,:) - &
        Lsdiag_mp%diag_4d(:,:,:,Lsdiag_mp_control%diag_pt%qadt_limits))* &
                                                               inv_dtcloud
      endif

! rain destruction
      if ( do_ncar_MG2  ) then
        do k=1,kx
          do j=1,jx
            do i=1,ix
              if ( qr_new(i,j,k) <= qmin .or. qnr_new(i,j,k) <= qmin ) then
                Cloud_state%SR_out(i,j,k) = Cloud_state%SR_out(i,j,k) -  &
                                                              qr_new(i,j,k)
                Cloud_state%SNR_out(i,j,k) = Cloud_state%SNR_out(i,j,k) -  &
                                                              qnr_new(i,j,k)

                Tend_mp%qtnd(i,j,k) = Tend_mp%qtnd(i,j,k) + qr_new(i,j,k)
                Tend_mp%ttnd(i,j,k) = Tend_mp%ttnd(i,j,k) - (hlv*qr_new(i,j,k))/cp_air

                if (Lsdiag_mp_control%diag_id%qrdt_destr > 0 .or. &
                    Lsdiag_mp_control%diag_id%qr_destr_col > 0) &
                    Lsdiag_mp%diag_4d(i,j,k,Lsdiag_mp_control%diag_pt%qrdt_destr) =    &
                                          - qr_new(i,j,k)/dtcloud
                if (Lsdiag_mp_control%diag_id%qnrdt_destr > 0 .or. &
                    Lsdiag_mp_control%diag_id%qnr_destr_col > 0) &
                    Lsdiag_mp%diag_4d(i,j,k,Lsdiag_mp_control%diag_pt%qnrdt_destr) =    &
                                          - qnr_new(i,j,k)/dtcloud
                if (Lsdiag_mp_control%diag_id%qdt_destr +    &
                              Lsdiag_mp_control%diag_id%q_destr_col > 0) &
                     Lsdiag_mp%diag_4d(i,j,k,Lsdiag_mp_control%diag_pt%qdt_destr) = &
                     Lsdiag_mp%diag_4d(i,j,k,Lsdiag_mp_control%diag_pt%qdt_destr) + &
                        qr_new(i,j,k)/dtcloud
              endif
            enddo
          enddo
        enddo

! snow destruction
        do k=1,kx
          do j=1,jx
            do i=1,ix
              if ( qs_new(i,j,k) <= qmin .or. qns_new(i,j,k) <= qmin) then
                Cloud_state%SS_out(i,j,k) = Cloud_state%SS_out(i,j,k) -  &
                                                              qs_new(i,j,k)
                Cloud_state%SNS_out(i,j,k) = Cloud_state%SNS_out(i,j,k) -  &
                                                              qns_new(i,j,k)
                Tend_mp%qtnd(i,j,k) = Tend_mp%qtnd(i,j,k) + qs_new(i,j,k)
                Tend_mp%ttnd(i,j,k) = Tend_mp%ttnd(i,j,k) - (hls*qs_new(i,j,k))/cp_air

                if (Lsdiag_mp_control%diag_id%qsdt_destr > 0 .or. &
                    Lsdiag_mp_control%diag_id%qs_destr_col > 0) &
                    Lsdiag_mp%diag_4d(i,j,k,Lsdiag_mp_control%diag_pt%qsdt_destr) =    &
                                          - qs_new(i,j,k)/dtcloud
                if (Lsdiag_mp_control%diag_id%qnsdt_destr > 0 .or. &
                    Lsdiag_mp_control%diag_id%qns_destr_col > 0) &
                    Lsdiag_mp%diag_4d(i,j,k,Lsdiag_mp_control%diag_pt%qnsdt_destr) =    &
                                          - qns_new(i,j,k)/dtcloud
                if (Lsdiag_mp_control%diag_id%qdt_destr +      &
                    Lsdiag_mp_control%diag_id%q_destr_col > 0) &
                    Lsdiag_mp%diag_4d(i,j,k,Lsdiag_mp_control%diag_pt%qdt_destr) = &
                    Lsdiag_mp%diag_4d(i,j,k,Lsdiag_mp_control%diag_pt%qdt_destr) + &
                         qs_new(i,j,k)/dtcloud
              endif
            enddo
          enddo
        enddo
      endif   ! do_ncar_MG2
!-----------------------------------------------------------------------
end subroutine destroy_tiny_clouds


!########################################################################


          end module ls_cloud_microphysics_mod
