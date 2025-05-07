
module hanedmf_mod

!-----------------------------------------------------------------------
use            mpp_mod, only:  input_nml_file
use            fms_mod, only:  mpp_pe, mpp_root_pe, stdlog, &
                               error_mesg, &
                               check_nml_error, FATAL, &
                               write_version_number, &
                               stdout, mpp_chksum
use sat_vapor_pres_mod, only:  compute_qs
use      constants_mod, only:  HLv,HLf,Cp_Air,Grav,rdgas,rvgas, &
                               kappa,vonkarm,tfreeze,HLs
use   diag_manager_mod, only:  register_diag_field, send_data
use   time_manager_mod, only:  time_type, get_date
use  satmedmfvdifq_mod, only:  satmedmfvdifq_init, satmedmfvdifq_run
use  monin_obukhov_mod, only:  mo_diff
use  aerosol_types_mod, only : aerosol_type      ! ZNT 06/22/2023
use  aer_ccn_act_mod,   only : aer_ccn_act_init  ! ZNT 06/22/2023


implicit none
private
!-----------------------------------------------------------------------
!  ---- public interfaces ----

   public  hanedmf, hanedmf_init

!-----------------------------------------------------------------------
!   ---- version number ----

 character(len=128) :: version = '$Id$'
 character(len=128) :: tagname = '$Name$'
 logical            :: module_is_initialized=.false.

!-----------------------------------------------------------------------
!   ---- local/private data ----

    real, parameter :: fv  = rvgas/rdgas-1.
    real, parameter :: eps = rdgas/rvgas
    real, parameter :: epsm1 = rdgas/rvgas-1.


!-----------------------------------------------------------------------
!   --- namelist ----

integer :: edmf_level = 2  ! Set the EDMF running mode (replacing the logical flags in vert_turb_driver).
                           ! edmf_level = 0: equivalent to previous 'do_hanedmf_diag',
                           !                 Run EDMF for diagnostics only (not using EDMF for SGS transport)
                           ! edmf_level = 1: equivalent to previous 'do_hanedmf_mfonly'
                           !                 Use the MF tendencies only (not using the ED coefficients)
                           ! edmf_level = 2: equivalent to previous 'do_hanedmf_full',
                           !                 Use both MF tendencies and ED coefficients for SGS transport
                           !                 * This is the full implementation of NCEP TKE-based EDMF 

logical :: do_debug =.false.
logical :: do_tke_eq = .false.
integer :: tke_eq_iter = 1
real    :: cfl_crit = 1.0
real    :: ck0 = 0.4
real    :: ck1 = 0.15
real    :: ch0 = 0.4
real    :: ch1 = 0.15
real    :: ce0 = 0.4
integer :: nstep_edmf = 1
logical :: do_cmp_qs = .false.
integer :: use_hl = 0     ! ZNT 06/12/2023: use hl instead of thl as updraft/downdraft conserved variable
                          !                 to be more consistent with ED and subs-detr form for MF
                          !                 0=no (original, use thl); 
                          !                 1=Tu (hl only for computation of Tu for subs-detr form of grid-mean tendencies); 
                          !                 2=bu (hl only for computation of bu for wu/wd and MF calculation);
                          !                 3=both (hl for both).
integer :: tkeh_flg = 0   ! ZNT 10/13/2020: option for TKE interpolation 
                          ! to half level for diffusivity
                          ! 0 = arithmetic mean; 1 = geometric mean; 2 = harmonic mean
real    :: cldtime = 500. ! ZNT 11/11/2020: Scu top cooling timescale (for downdraft initial Tpert)
real    :: upd_fac = 0.13
real    :: dnd_fac = 0.10
real    :: qlia_tnd_fac = 1.0
real    :: qa_mf_fac  = 0.0 ! ZNT 06/01/2023: MF factor for cloud fraction (qa). Model before 202305 = 0.0
real    :: qn_mf_fac  = 0.0 ! ZNT 05/12/2023: MF factor for droplet number (qn). Model before 202304 = 0.0
real    :: qni_mf_fac = 0.0 ! ZNT 05/12/2023: MF factor for ice number (qni). Model before 202304 = 0.0
logical :: do_tqp_mf = .false. ! ZNT 05/19/2023: MF transports passive T,qv,ql,qi (instead of interactive with condensation)
logical :: do_subdtr = .false. ! ZNT 06/01/2023: MF transports interactive T,qv,ql,qi but using subsidence-detrainment form
real    :: shflx_max = 1000.0  ! ZNT 10/14/2020: capping sfcflx to prevent huge T changes during the split step.
                               ! This should only be active during model spin up where Sfc does not match Atm.
                               ! Default = 1000.0 W/m2 corresponds to deltaT=45K in the lowest layer over 1800s.
                               ! The cap of sfcflx only affects the diffusivity, not actual 'vdif' tendencies.
logical :: do_switch_ed_bug = .true. ! ZNT 10/15/2020: the previous version erronously mess up the ED for 
                                     ! for momentum and heat. Setting this to true retains the results pre-10/14/2020.
                                     ! The correct run should be with .false. for this flag. 
logical :: do_upd_ovrsht = .false.   ! ZNT 11/18/2020: .true. reproduces model behavior before mid-2020,
                                     ! i.e. updraft can overshoot the dry PBL height and dissipates Sc.
logical :: use_lock_zpbl = .false.   ! ZNT 12/08/2020: Use Lock scheme's ZPBL calculation
real    :: parcel_buoy = 0.25  ! ZNT 12/08/2020: From Lock scheme's ZPBL calculation.
                               ! scaling factor for surface parcel buoyancy. The value 0.25 is from AM4-orig.
real    :: critjump    = 0.1   ! ZNT 12/08/2020: From Lock scheme's ZPBL calculation.
                               ! critical jump for finding stable interfaces to bound convective layers, or to 
                               ! identify ambigous layers (K). The value 0.1 is from AM4-orig.

logical :: do_tkevdif = .true. ! ZNT 06/05/2021: True: tkevdif computed in EDMF; False: tkevdif computed in ver_diff. 
                               !                 Should set to false for AM4-EDMF?
logical :: do_tkediss_imp = .false.  ! ZNT 06/05/2021: implicit computation of TKE dissipation. Not implemented yet.     
real    :: mf_buop_fac = 1.0   ! ZNT 06/05/2021: Debugging: factor for buoyancy production by MF
real    :: mf_shrp_fac = 1.0   ! ZNT 06/05/2021: Debugging: factor for shear production by MF
logical :: do_pos_mf_buop = .false.  ! ZNT 06/05/2021: True: limit buoyancy production by MF to be positive only
logical :: do_cfl_col  = .false.     ! ZNT 06/06/2021: True: column-wise MF limiter; False: local MF limiter
real    :: mf_tkeflx_fac = 1.0 ! ZNT 06/06/2021: factor to limit the MF transport of TKE
real    :: c_entu = 0.4        ! ZNT 06/21/2021: entrainment parameter for updraft
real    :: c_entd = 0.4        ! ZNT 06/21/2021: entrainment parameter for downdraft
real    :: actei  = 0.7        ! ZNT 06/21/2021: cloud top entrainment instability parameter

logical :: use_z_above_sfc = .false. ! ZNT 03/08/2023: .false. reproduces results in 2021-2022; .true. is the correct implementation
logical :: do_km_L65 = .false. ! ZNT 10/02/2023: False for L33 where kmpbl and kmscu = km/2; 
                               !                  True for L65 where kmpbl and kmscu = km/3. (10/04/2023: Corrected from 2/3)
real    :: z_ed_top = 1.0e8    ! ZNT 12/16/2023: Top height for TKE-based ED. 
                               !                 Above this height and convective PBL top, the TKE-based ED is set to zero.
logical :: do_wush = .false.   ! ZNT 12/24/2023: (From GFS) Shear impact on updraft/downdraft velocity
logical :: do_tkemean = .false.! ZNT 12/24/2023: (From GFS) Subcloud mean TKE impact on updraft/downdraft entrainment


! ZNT 06/22/2023: Namelist parameters for aerosol activation settings using AM4 default parameters for now. 
!                 They should be changed to using directly the parameters in moist_processes_nml in the future.
logical :: use_online_aerosol = .true.
logical :: use_sub_seasalt = .false.
real    :: sea_salt_scale = 0.1
real    :: om_to_oc = 1.67


namelist /hanedmf_nml/ edmf_level, do_debug, do_tke_eq, tke_eq_iter, cfl_crit, &   ! ZNT 01/25/2024: edmf_level
        ck0,ck1,ch0,ch1,ce0, nstep_edmf, do_cmp_qs, use_hl, tkeh_flg, cldtime,  &  ! ZNT 06/12/2023: use_hl
        c_entu, c_entd, actei, &
        upd_fac, dnd_fac, qlia_tnd_fac, shflx_max, do_switch_ed_bug, do_upd_ovrsht, &
        use_lock_zpbl, parcel_buoy, critjump, do_tkevdif, do_tkediss_imp, &
        mf_buop_fac, mf_shrp_fac, do_pos_mf_buop, do_cfl_col, mf_tkeflx_fac, &
        use_z_above_sfc, &    ! ZNT 03/08/2023
        do_subdtr, qa_mf_fac, & ! ZNT 06/01/2023
        qn_mf_fac, qni_mf_fac, do_tqp_mf, & ! ZNT 05/19/2023
        do_km_L65, z_ed_top, do_wush, do_tkemean, &  ! ZNT 10/02/2023; 12/16/2023; 12/24/2023
        use_online_aerosol, use_sub_seasalt, sea_salt_scale, om_to_oc  ! ZNT 06/22/2023

!---------------------------------------------------------------------
! DIAGNOSTICS FIELDS 
!---------------------------------------------------------------------
integer :: id_zfull_han, id_zhalf_han, id_pfull_han, id_phalf_han   ! ZNT 06/20/2023

integer :: id_tdt_han, id_qdt_han, id_udt_han, id_vdt_han, &
           id_qldt_han, id_qidt_han, id_qadt_han, id_tke_han, &
           id_qndt_han, id_qnidt_han,                                   &         ! ZNT 05/03/2023
           id_km_han, id_kt_han, id_ke_han, id_mu_han, id_md_han, &
           id_shf_han, id_lhf_han, id_taux_han, id_tauy_han, & 
           id_risfc_han, id_fm_han, id_ft_han, id_zpbl_han
! ZNT 05/08/2020: MF diagnostics
integer :: id_tdt_mf_han, id_qdt_mf_han, id_udt_mf_han, id_vdt_mf_han, &
           id_shf_mf_han, id_lhf_mf_han, id_taux_mf_han, id_tauy_mf_han, &
           id_qldt_mf_han, id_qidt_mf_han, id_qadt_mf_han,                &
           id_qndt_mf_han, id_qnidt_mf_han,                               &       ! ZNT 05/03/2023
           id_tdt_mfd_han, id_qdt_mfd_han, id_udt_mfd_han, id_vdt_mfd_han, &
           id_shf_mfd_han, id_lhf_mfd_han, id_taux_mfd_han, id_tauy_mfd_han, &
           id_qldt_mfd_han, id_qidt_mfd_han, id_qadt_mfd_han,  &
           id_qndt_mfd_han, id_qnidt_mfd_han                                      ! ZNT 05/03/2023
! End ZNT 05/08/2020

! ZNT 10/03/2020: TKE diagnostics
integer :: id_de_buop_han, id_de_shrp_han, id_de_diss_han
! End ZNT 10/03/2020

integer :: id_tcko_han, id_qvcko_han, id_qlcko_han, id_qicko_han, id_buou_han, &  ! ZNT 06/19/2020
           id_tcdo_han, id_qvcdo_han, id_qlcdo_han, id_qicdo_han, id_buod_han     ! ZNT 06/19/2020
integer :: id_qacko_han, id_qncko_han, id_qnicko_han, &  ! ZNT 05/12/2023; 05/24/2023
           id_qacdo_han, id_qncdo_han, id_qnicdo_han     ! ZNT 05/12/2023; 05/24/2023


! ZNT 05/17/2023: passive MF diagnostics
integer :: id_tpcko_han, id_qvpcko_han, id_qlpcko_han, id_qipcko_han, &  ! ZNT 05/17/2023: passive T,q,ql,qi in updraft
           id_qapcko_han, id_qnpcko_han, id_qnipcko_han, &               ! ZNT 07/05/2023: passive qa,qn,qni in updraft
           id_tpcdo_han, id_qvpcdo_han, id_qlpcdo_han, id_qipcdo_han, &  ! ZNT 05/17/2023: passive T,q,ql,qi in downdraft
           id_qapcdo_han, id_qnpcdo_han, id_qnipcdo_han, &               ! ZNT 07/05/2023: passive qa,qn,qni in downdraft
           id_tpdt_han,  id_qpdt_han,   id_qlpdt_han,  id_qipdt_han, &   ! ZNT 05/17/2023: passive T,q,ql,qi tendencies by EDMF
           id_qapdt_han, id_qnpdt_han,  id_qnipdt_han, &                 ! ZNT 07/05/2023: passive qa,qn,qni tendencies by EDMF
           id_tpdt_mf_han, id_qpdt_mf_han, id_qlpdt_mf_han, id_qipdt_mf_han, &     ! ZNT 05/17/2023: passive tendencies by updraft
           id_qapdt_mf_han, id_qnpdt_mf_han, id_qnipdt_mf_han, &         ! ZNT 07/05/2023: passive qa,qn,qni tendencies by updraft
           id_tpdt_mfd_han, id_qpdt_mfd_han, id_qlpdt_mfd_han, id_qipdt_mfd_han, & ! ZNT 05/17/2023: passive tendencies by downdraft
           id_qapdt_mfd_han, id_qnpdt_mfd_han, id_qnipdt_mfd_han         ! ZNT 07/05/2023: passive qa,qn,qni tendencies by downdraft


! ZNT 05/24/2023: Subs-detr tendencies diagnostics (Only for T, q, ql, qa, qn, qni; without u, v)
integer :: id_tdt_sub_han,   id_qdt_sub_han,   id_qldt_sub_han,   id_qidt_sub_han,  &
           id_qadt_sub_han,  id_qndt_sub_han,  id_qnidt_sub_han,  & 
           id_tdt_subd_han,  id_qdt_subd_han,  id_qldt_subd_han,  id_qidt_subd_han, &
           id_qadt_subd_han, id_qndt_subd_han, id_qnidt_subd_han
integer :: id_tdt_dtr_han,   id_qdt_dtr_han,   id_qldt_dtr_han,   id_qidt_dtr_han,  &
           id_qadt_dtr_han,  id_qndt_dtr_han,  id_qnidt_dtr_han,  & 
           id_tdt_dtrd_han,  id_qdt_dtrd_han,  id_qldt_dtrd_han,  id_qidt_dtrd_han, &
           id_qadt_dtrd_han, id_qndt_dtrd_han, id_qnidt_dtrd_han
                     
! ZNT 05/23/2023: passive MF diagnostics for subs-detr decomposition
integer :: id_tpdt_sub_han, id_qpdt_sub_han, id_qlpdt_sub_han, id_qipdt_sub_han,  &  ! ZNT 05/23/2023: passive tendencies by updraft-subsidence
           id_qapdt_sub_han, id_qnpdt_sub_han, id_qnipdt_sub_han                     ! ZNT 07/05/2023
integer :: id_tpdt_subd_han,id_qpdt_subd_han,id_qlpdt_subd_han,id_qipdt_subd_han, &  ! ZNT 05/23/2023: passive tendencies by downdraft-ascent
           id_qapdt_subd_han, id_qnpdt_subd_han, id_qnipdt_subd_han                  ! ZNT 07/05/2023
integer :: id_tpdt_dtr_han, id_qpdt_dtr_han, id_qlpdt_dtr_han, id_qipdt_dtr_han,  &  ! ZNT 05/23/2023: passive tendencies by updraft-detrainment
           id_qapdt_dtr_han, id_qnpdt_dtr_han, id_qnipdt_dtr_han                     ! ZNT 07/05/2023
integer :: id_tpdt_dtrd_han,id_qpdt_dtrd_han,id_qlpdt_dtrd_han,id_qipdt_dtrd_han, &  ! ZNT 05/23/2023: passive tendencies by downdraft-detrainment
           id_qapdt_dtrd_han, id_qnpdt_dtrd_han, id_qnipdt_dtrd_han                  ! ZNT 07/05/2023


integer :: id_zsml_han ! ZNT 12/08/2020

character(len=14) :: mod_name = 'hanedmf'

real :: missing_value = -999.
integer :: icount = 0

integer :: do_print_aer = 0   ! ZNT 06/27/2023: debug output aerosol array

contains

!#######################################################################

   subroutine hanedmf (is,ie,js,je,dt,Time_next,tdtlw_in,rough,u_star,b_star, &
                       u_ref,v_ref,tau_x,tau_y,t_surf_rad, &
                       shflx,lhflx,wind,thv_atm,thv_surf,  &
                       t,qv,ql,qi,qa,qn,qni,  &                        ! ZNT 05/03/2023: adding qn, qni
                       u,v,zfull,pfull,zhalf,phalf,tke, &
                       udt_mfud,vdt_mfud,tdt_mfud,qdt_mfud, &          ! ZNT 06/05/2020
                       qldt_mfud,qidt_mfud,qadt_mfud, &                
                       qndt_mfud,qnidt_mfud, &                         ! ZNT 05/03/2023: adding qn, qni
                       diff_m,diff_t, &                                ! ZNT 10/11/2020 
                       hpbl_out, &
                       use_z_above_sfc_out, &                          ! ZNT 03/08/2023
                       kbot, asol)                                     ! ZNT 12/08/2020; 06/27/2023
!-----------------------------------------------------------------------
!
!                      driver for han's edmf
!
!-----------------------------------------------------------------------
!
!      is,ie,js,je  i,j indices marking the slab of model working on
!      time      variable needed for netcdf diagnostics
!
!      tdtlw_in  lw radiative tendency (K/s)
!      rough     surface roughness (m)
!      u_star    friction velocity (m/s)
!      b_star    buoyancy scale (m/s**2)
!
!      three dimensional fields on model full levels, reals dimensioned
!      (:,:,pressure), third index running from top of atmosphere to 
!      bottom
!          
!      t         temperature (K)
!      qv        water vapor specific humidity (kg vapor/kg air)
!      ql        liquid water specific humidity (kg cond/kg air)
!      qi        ice water specific humidity (kg cond/kg air)
!      qa        cloud fraction 
!      qn        cloud droplet number (1/kg air, ZNT added 05/2023)
!      qni       cloud ice number (1/kg air, ZNT added 05/2023)
!      zfull     height of full levels (m)
!      pfull     pressure (Pa)
!      u         zonal wind (m/s)
!      v         meridional wind (m/s)
!
!      the following two fields are on the model half levels, with
!      size(zhalf,3) = size(t,3) +1, zhalf(:,:,size(zhalf,3)) 
!      must be height of surface (if you are not using eta-model)
!
!      zhalf     height at half levels (m)
!      phalf     pressure at half levels (Pa)
!
!-----------------------------------------------------------------------
!--------------------- interface arguments -----------------------------
integer,         intent(in)                      :: is,ie,js,je
type(time_type), intent(in)                      :: Time_next
real,            intent(in)                      :: dt
real,            intent(in),    dimension(:,:,:) :: tdtlw_in
real,            intent(in),    dimension(:,:)   :: rough,u_star,b_star
real,            intent(in),    dimension(:,:)   :: u_ref, v_ref
real,            intent(in),    dimension(:,:)   :: tau_x, tau_y
real,            intent(in),    dimension(:,:)   :: t_surf_rad
real,            intent(in),    dimension(:,:)   :: shflx, lhflx
real,            intent(in),    dimension(:,:)   :: wind,thv_atm,thv_surf
real,            intent(in),    dimension(:,:,:) :: t,qv,ql,qi,qa
real,            intent(in),    dimension(:,:,:) :: qn,qni                ! ZNT 05/03/2023 
real,            intent(in),    dimension(:,:,:) :: u,v,zfull,pfull
real,            intent(in),    dimension(:,:,:) :: zhalf, phalf
real,            intent(inout), dimension(:,:,:) :: tke
real,            intent(out),   dimension(:,:,:) :: udt_mfud,  vdt_mfud,  &
                                                    tdt_mfud,  qdt_mfud   ! ZNT 06/05/2020 - dimension nx,ny,nz; MF only
real,            intent(out),   dimension(:,:,:) :: qldt_mfud, qidt_mfud, &
                                                    qadt_mfud             ! ZNT 10/11/2020 - dimension nx,ny,nz; MF only
real,            intent(out),   dimension(:,:,:) :: qndt_mfud, qnidt_mfud ! ZNT 05/03/2023
real,            intent(out),   dimension(:,:,:) :: diff_m, diff_t        ! ZNT 06/05/2020 - dimension nx,ny,nz+1
real,            intent(out),   dimension(:,:)   :: hpbl_out              ! ZNT 12/08/2020
logical,         intent(out)                     :: use_z_above_sfc_out   ! ZNT 03/09/2023 - to be used in vert_turb_driver
integer,  intent(in),   dimension(:,:), optional :: kbot
type(aerosol_type), intent(in),         optional :: asol     ! ZNT 06/27/2023

!-----------------------------------------------------------------------
!---------------------- local data -------------------------------------

   real,dimension(size(t,1),size(t,2),size(t,3)) :: &
         tdt_out, qvdt_out, qldt_out, qidt_out, qadt_out, udt_out, vdt_out,  &
         tcko_out, qvcko_out, qlcko_out, qicko_out, qacko_out, buou_out, &    ! ZNT 06/19/2020; 05/24/2023: adding qa
         tcdo_out, qvcdo_out, qlcdo_out, qicdo_out, qacdo_out, buod_out, &    ! ZNT 06/19/2020; 05/24/2023: adding qa
         qndt_out, qnidt_out, qncko_out, qnicko_out, qncdo_out, qnicdo_out, & ! ZNT 05/03/2023: adding qn, qni
         xmf_out, xmfd_out ! , tke
! ZNT 05/22/2020: TKE as a diagnostic tracer now
   real,dimension(size(t,1),size(t,2),size(t,3)+1) :: &
         dku_out, dkt_out, dkq_out
   real,dimension(size(t,1),size(t,2)) :: &
         dusfc_out, dvsfc_out, dtsfc_out, dqsfc_out, &   ! ZNT 12/08/2020: hpbl_out moved to output
         fm_out, ft_out, risfc_out
! ZNT 05/08/2020, 10/11/2020: Diagnostic MF output
   real,dimension(size(t,1),size(t,2),size(t,3)) :: &
         tdt_mf_out,  qvdt_mf_out,  udt_mf_out,  vdt_mf_out,  &
         qldt_mf_out, qidt_mf_out, qadt_mf_out,               &
         qndt_mf_out, qnidt_mf_out,                 &   ! ZNT 05/03/2023: adding qn, qni
         tdt_mfd_out, qvdt_mfd_out, udt_mfd_out, vdt_mfd_out, &
         qldt_mfd_out, qidt_mfd_out, qadt_mfd_out,            &
         qndt_mfd_out, qnidt_mfd_out                    ! ZNT 05/03/2023: adding qn, qni
   real,dimension(size(t,1),size(t,2)) :: &
         dusfc_mf_out,  dvsfc_mf_out,  dtsfc_mf_out,  dqsfc_mf_out, &
         dusfc_mfd_out, dvsfc_mfd_out, dtsfc_mfd_out, dqsfc_mfd_out
! End ZNT 05/08/2020, 10/11/2020

! ZNT 05/17/2023: Arrays for passive variables
   real,dimension(size(t,1),size(t,2),size(t,3),7) :: &    ! ZNT 07/05/2023: add qa,qn,qni (4 -> 7)
         tqpdt_out, tqpcko_out, tqpcdo_out, tqpdt_mf_out, tqpdt_mfd_out

! ZNT 05/24/2023: Arrays for Subs-detr tendencies  (Only for T, q, ql, qa, qn, qni; without u, v)
   real,dimension(size(t,1),size(t,2),size(t,3)) :: &
         tdt_sub_out,   qvdt_sub_out,  qldt_sub_out,  qidt_sub_out,  &
         qadt_sub_out,  qndt_sub_out,  qnidt_sub_out,  &
         tdt_subd_out,  qvdt_subd_out, qldt_subd_out, qidt_subd_out, &
         qadt_subd_out, qndt_subd_out, qnidt_subd_out, &
         tdt_dtr_out,   qvdt_dtr_out,  qldt_dtr_out,  qidt_dtr_out,  &
         qadt_dtr_out,  qndt_dtr_out,  qnidt_dtr_out,  &
         tdt_dtrd_out,  qvdt_dtrd_out, qldt_dtrd_out, qidt_dtrd_out, &
         qadt_dtrd_out, qndt_dtrd_out, qnidt_dtrd_out
         
! ZNT 05/23/2023: Arrays for Subs-detr tendencies for passive variables
   real,dimension(size(t,1),size(t,2),size(t,3),7) :: &   ! ZNT 07/05/2023: add qa,qn,qni (4 -> 7)
         tqpdt_sub_out, tqpdt_subd_out, tqpdt_dtr_out, tqpdt_dtrd_out
         
! ZNT 10/03/2020: Diagnostic TKE output
   real,dimension(size(t,1),size(t,2),size(t,3)) :: &
         de_buop_out, de_shrp_out, de_diss_out
! End ZNT 10/03/2020

! ZNT 12/08/2020: variables for zpbl_depth
   real,dimension(size(t,1),size(t,2),size(t,3)) :: zfull_ag, slv, hleff
   real,dimension(size(t,1),size(t,2),size(t,3)+1) :: zhalf_ag
   real,dimension(size(t,1),size(t,2)) :: zsml, zsurf
   integer :: ipbl, ibot, nlev
   real :: tmpjump
   real :: d608 = (rvgas-rdgas)/rdgas

! ZNT 06/22/2023: additional variables for aerosol activation
   real, dimension(size(t,1),size(t,2),size(t,3)) :: pmass,    &  ! layer mass (kg/m2)
          am1, am2, am3, am4, am5, amx1, amx2, amx3, amx4, amx5


!-------------------- Han EDMF arrays ----------------------------------
! ZNT 10/02/2020: Added '_1' variables for single-step values

integer, dimension(1)  :: kinver       ! (I) index of highest T-inversion
integer, dimension(1)  :: kpbl,kpbl_1  ! (O) index of PBL top

real :: delt,delt_1 ! (I) time step
real :: xkzm_m      ! (I) bkgnd momentum diff
real :: xkzm_h      ! (I) bkgnd heat diff
real :: xkzm_s      ! (I) sigma threshold for bkgnd diff
real :: dspfac      ! (I) TKE diss. heating factor
real :: bl_upfr     ! (I) updraft frac.
real :: bl_dnfr     ! (I) downdraft frac.

real, dimension(1,size(t,3))   :: dv,dv_1,   & ! (IO) dV/dt
                                  du,du_1,   & ! (IO) dU/dt
                                  tdt,tdt_1    ! (IO) dT/dt
! real, dimension(1,size(t,3),7) :: rtg,rtg_1    ! (IO) d(trac)/dt  <- ZNT 05/03/2023: adding qn, qni (5+2=7 tracers)
! real, dimension(1,size(t,3),11) :: rtg,rtg_1    ! (IO) d(trac)/dt  <- ZNT 05/17/2023: adding passive T,q,ql,qi (7+4=11 tracers)
real, dimension(1,size(t,3),14) :: rtg,rtg_1    ! (IO) d(trac)/dt  <- ZNT 07/05/2023: adding passive qa,qn,qni (11+3=14 tracers)

real, dimension(1,size(t,3))   :: u1,u1_1,   & ! (I) U-wind
                                  v1,v1_1,   & ! (I) V-wind
                                  t1,t1_1      ! (I) Temperature
! real, dimension(1,size(t,3),7) :: q1,q1_1      ! (I) Tracers      <- ZNT 05/03/2023: adding qn, qni (5+2=7 tracers)
! real, dimension(1,size(t,3),11) :: q1,q1_1      ! (I) Tracers      <- ZNT 05/17/2023: adding passive T,q,ql,qi (7+4=11 tracers)
real, dimension(1,size(t,3),14) :: q1,q1_1      ! (I) Tracers      <- ZNT 07/05/2023: adding passive qa,qn,qni (11+3=14 tracers)

real, dimension(1,size(t,3),4) :: aer1, aerx1  ! (I) Aerosol    <- ZNT 06/28/2023: adding aerosol arrays to pass into updraft

real, dimension(1,size(t,3))   :: swh,  & ! (I) SW heating rate
                                  hlw     ! (I) LW heating rate
real, dimension(1)             :: xmu,  & ! (I) Zenith angle factor
                                garea,  & ! (I) Area of grid cell
                                zvfun,  & ! (I) Vegetation parameter
                                  psk,  & ! (I) Exner func. at SFC
                               rbsoil,  & ! (I) Bulk Ri-number at SFC
                                 zorl,  & ! (I) SFC roughness length
                                 tsea,  & ! (I) SFC skin temperature
                                 u10m,  & ! (I) U-wind at 10m
                                 v10m,  & ! (I) V-wind at 10m
                                   fm,  & ! (I) M-O func for mom.
                                   fh,  & ! (I) M-O func for heat
                                 evap,  & ! (I) SFC moisture flux
                                 heat,  & ! (I) SFC temperature flux
                               stress,  & ! (I) SFC wind stress
                                 spd1     ! (I) lowest level wind speed

real, dimension(1,size(t,3)+1):: prsi,  & ! (I) pressure at half levels
                                 phii     ! (I) geopotential at half levels  <- ZNT 03/08/2023: Changed to relative phi above surface (instead of sea-level)
real, dimension(1,size(t,3))  ::  del,  & ! (I) pressure thickness
                                 prsl,  & ! (I) pressure at full levels
                                prslk,  & ! (I) Exner func. at full levels
                                 phil     ! (I) geopotential at full levels  <- ZNT 03/08/2023: Changed to relative phi above surface (instead of sea-level)

logical :: dspheat   ! (I) flag for using TKE diss. heating
integer :: do_mf_aer ! (I) flag for using aerosol activation in updraft <- ZNT 06/28/2023: 0 = skip; 1 = 'do_online_aerosol = .false.'; 2 = 'do_online_aerosol = .true.'


real, dimension(1)            :: dusfc,dusfc_1, & ! (O) diag. X-mom. flux
                                 dvsfc,dvsfc_1, & ! (O) diag. Y-mom. flux
                                 dtsfc,dtsfc_1, & ! (O) diag. SFC SHF
                                 dqsfc,dqsfc_1, & ! (O) diag. SFC LHF
                                  hpbl,hpbl_1     ! (O) height of PBL top

real, dimension(1,size(t,3)-1)::   dku,dku_1, & ! (O) momentum eddy diff.
                                   dkt,dkt_1, & ! (O) heat eddy diff.
                                   dkq,dkq_1    ! (O) tke eddy diff.
real, dimension(1,size(t,3))  ::   xmf,xmf_1, & ! (O) updraft mass flux
                                  xmfd,xmfd_1   ! (O) downdraft mass flux

! ZNT 06/19/2020: Plume variables
real, dimension(1,size(t,3))  ::  tcko, & ! (O) updraft temperature
                                 qvcko, & ! (O) updraft qv
                                 qlcko, & ! (O) updraft ql
                                 qicko, & ! (O) updraft qi
                                 qacko, & ! (O) updraft qa            <- ZNT 05/24/2023
                                 qncko, & ! (O) updraft qn            <- ZNT 05/03/2023
                                qnicko, & ! (O) updraft qni           <- ZNT 05/03/2023
                                  buou, & ! (O) updraft buoyancy
                                  tcdo, & ! (O) downdraft temperature
                                 qvcdo, & ! (O) downdraft qv
                                 qlcdo, & ! (O) downdraft ql
                                 qicdo, & ! (O) downdraft qi
                                 qacdo, & ! (O) downdraft qa          <- ZNT 05/24/2023
                                 qncdo, & ! (O) downdraft qn          <- ZNT 05/03/2023
                                qnicdo, & ! (O) downdraft qni         <- ZNT 05/03/2023 
                                  buod    ! (O) downdraft buoyancy

! ZNT 05/08/2020, 10/11/2020: Diagnostic MF output
real, dimension(1,size(t,3))   :: dv_mf,   dv_mf_1,  & ! (O) dV/dt by updraft MF
                                  du_mf,   du_mf_1,  & ! (O) dU/dt by updraft MF
                                 tdt_mf,  tdt_mf_1,  & ! (O) dT/dt by updraft MF
                                 qdt_mf,  qdt_mf_1,  & ! (O) dq/dt by updraft MF
                                qldt_mf, qldt_mf_1,  & ! (O) dql/dt by updraft MF
                                qidt_mf, qidt_mf_1,  & ! (O) dqi/dt by updraft MF
                                qadt_mf, qadt_mf_1,  & ! (O) dqa/dt by updraft MF
                                qndt_mf, qndt_mf_1,  & ! (O) dqn/dt by updraft MF    <- ZNT 05/03/2023
                               qnidt_mf, qnidt_mf_1, & ! (O) dqni/dt by updraft MF   <- ZNT 05/03/2023
                                  dv_mfd,  dv_mfd_1, & ! (O) dV/dt by downdraft MF
                                  du_mfd,  du_mfd_1, & ! (O) dU/dt by downdraft MF
                                 tdt_mfd, tdt_mfd_1, & ! (O) dT/dt by downdraft MF
                                 qdt_mfd, qdt_mfd_1, & ! (O) dq/dt by downdraft MF
                                qldt_mfd,qldt_mfd_1, & ! (O) dql/dt by downdraft MF
                                qidt_mfd,qidt_mfd_1, & ! (O) dql/dt by downdraft MF
                                qadt_mfd,qadt_mfd_1, & ! (O) dqa/dt by downdraft MF
                                qndt_mfd,qndt_mfd_1, & ! (O) dqn/dt by downdraft MF  <- ZNT 05/03/2023
                               qnidt_mfd,qnidt_mfd_1   ! (O) dqni/dt by downdraft MF <- ZNT 05/03/2023

! ZNT 05/17/2023: Arrays for passive variables; 07/05/2023: add qa,qn,qni (4 -> 7)
real, dimension(1,size(t,3),7) ::  tqpcko, & ! (O) updraft passive T,q,ql,qi...
                                   tqpcdo, & ! (O) downdraft passive T,q,ql,qi...
                                   tqpdt_mf, tqpdt_mf_1, & ! (O) d(T,q,ql,qi...)/dt by updraft MF
                                   tqpdt_mfd, tqpdt_mfd_1  ! (O) d(T,q,ql,qi...)/dt by downdraft MF

! ZNT 05/24/2023: Arrays for Subs-detr tendencies  (Only for T, q, ql, qa, qn, qni; without u, v)
real, dimension(1,size(t,3))  :: tdt_sub,  tdt_sub_1,  & ! (O) dT/dt by updraft MF - subsidence component
                                 qdt_sub,  qdt_sub_1,  & ! (O) dq/dt by updraft MF - subsidence component
                                qldt_sub, qldt_sub_1,  & ! (O) dql/dt by updraft MF - subsidence component
                                qidt_sub, qidt_sub_1,  & ! (O) dqi/dt by updraft MF - subsidence component
                                qadt_sub, qadt_sub_1,  & ! (O) dqa/dt by updraft MF - subsidence component
                                qndt_sub, qndt_sub_1,  & ! (O) dqn/dt by updraft MF - subsidence component
                               qnidt_sub, qnidt_sub_1, & ! (O) dqni/dt by updraft MF - subsidence component
                                 tdt_subd, tdt_subd_1, & ! (O) dT/dt by downdraft MF - ascent component
                                 qdt_subd, qdt_subd_1, & ! (O) dq/dt by downdraft MF - ascent component
                                qldt_subd,qldt_subd_1, & ! (O) dql/dt by downdraft MF - ascent component
                                qidt_subd,qidt_subd_1, & ! (O) dql/dt by downdraft MF - ascent component
                                qadt_subd,qadt_subd_1, & ! (O) dqa/dt by downdraft MF - ascent component
                                qndt_subd,qndt_subd_1, & ! (O) dqn/dt by downdraft MF - ascent component
                               qnidt_subd,qnidt_subd_1   ! (O) dqni/dt by downdraft MF - ascent component
real, dimension(1,size(t,3))  :: tdt_dtr,  tdt_dtr_1,  & ! (O) dT/dt by updraft MF - detrainment component
                                 qdt_dtr,  qdt_dtr_1,  & ! (O) dq/dt by updraft MF - detrainment component
                                qldt_dtr, qldt_dtr_1,  & ! (O) dql/dt by updraft MF - detrainment component
                                qidt_dtr, qidt_dtr_1,  & ! (O) dqi/dt by updraft MF - detrainment component
                                qadt_dtr, qadt_dtr_1,  & ! (O) dqa/dt by updraft MF - detrainment component
                                qndt_dtr, qndt_dtr_1,  & ! (O) dqn/dt by updraft MF - detrainment component
                               qnidt_dtr, qnidt_dtr_1, & ! (O) dqni/dt by updraft MF - detrainment component
                                 tdt_dtrd, tdt_dtrd_1, & ! (O) dT/dt by downdraft MF - detrainment component
                                 qdt_dtrd, qdt_dtrd_1, & ! (O) dq/dt by downdraft MF - detrainment component
                                qldt_dtrd,qldt_dtrd_1, & ! (O) dql/dt by downdraft MF - detrainment component
                                qidt_dtrd,qidt_dtrd_1, & ! (O) dql/dt by downdraft MF - detrainment component
                                qadt_dtrd,qadt_dtrd_1, & ! (O) dqa/dt by downdraft MF - detrainment component
                                qndt_dtrd,qndt_dtrd_1, & ! (O) dqn/dt by downdraft MF - detrainment component
                               qnidt_dtrd,qnidt_dtrd_1   ! (O) dqni/dt by downdraft MF - detrainment component

! ZNT 05/23/2023: Arrays for Subs-detr tendencies for passive variables
real, dimension(1,size(t,3),7) ::  tqpdt_sub,  tqpdt_sub_1, & ! (O) d(T,q,ql,qi...)/dt by updraft MF - subsidence component
                                   tqpdt_subd, tqpdt_subd_1   ! (O) d(T,q,ql,qi...)/dt by downdraft MF - ascent component     
real, dimension(1,size(t,3),7) ::  tqpdt_dtr,  tqpdt_dtr_1, & ! (O) d(T,q,ql,qi...)/dt by updraft MF - detrainment component
                                   tqpdt_dtrd, tqpdt_dtrd_1   ! (O) d(T,q,ql,qi...)/dt by downdraft MF - detrainment component     

real, dimension(1)            :: dusfc_mf, dusfc_mf_1, & ! (O) diag. X-mom. flux by Up MF (should be 0)
                                 dvsfc_mf, dvsfc_mf_1, & ! (O) diag. Y-mom. flux by Up MF (should be 0)
                                 dtsfc_mf, dtsfc_mf_1, & ! (O) diag. SFC SHF by Up MF (should be 0)
                                 dqsfc_mf, dqsfc_mf_1, & ! (O) diag. SFC LHF by Up MF (should be 0)
                                 dusfc_mfd,dusfc_mfd_1,& ! (O) diag. X-mom. flux by Dn MF (should be 0)
                                 dvsfc_mfd,dvsfc_mfd_1,& ! (O) diag. Y-mom. flux by Dn MF (should be 0)
                                 dtsfc_mfd,dtsfc_mfd_1,& ! (O) diag. SFC SHF by Dn MF (should be 0)
                                 dqsfc_mfd,dqsfc_mfd_1   ! (O) diag. SFC LHF by Dn MF (should be 0)
! End ZNT 05/08/2020, 10/11/2020

! ZNT 10/03/2020: Diagnostic TKE output
real, dimension(1,size(t,3))   :: de_buop, de_buop_1,  & ! (O) Buoyancy production of TKE
                                  de_shrp, de_shrp_1,  & ! (O) Shear production of TKE
                                  de_diss, de_diss_1     ! (O) Viscous dissipation of TKE
! End ZNT 10/03/2020

! Local variables
character(len=100):: errmsg
integer :: errflg
integer :: i,j,k, kx,kxp,istep_tke,nstep_tke,istep_edmf, itqp
real :: rho_a, delb
logical :: used

   if (.not. module_is_initialized) call error_mesg ('hanedmf',  &
                         'hanedmf_init has not been called.', FATAL)


   ! ZNT 03/08/2023: moved the calculation of zfull_ag and zhalf_ag to the beginning
   nlev = size(t,3)
   if (present(kbot)) then
     do i=1,size(t,1)
       do j=1,size(t,2)
          zsurf(i,j) = zhalf(i,j,kbot(i,j)+1)
       enddo
     enddo
   else
      zsurf(:,:) = zhalf(:,:,nlev+1)
   end if
   do k = 1, nlev
      zfull_ag(:,:,k) = zfull(:,:,k) - zsurf(:,:)
      zhalf_ag(:,:,k) = zhalf(:,:,k) - zsurf(:,:)
   enddo
   zhalf_ag(:,:,nlev+1) = zhalf(:,:,nlev+1) - zsurf(:,:)

   use_z_above_sfc_out = use_z_above_sfc

   ! Common HanEDMF parameters   
   kx=size(t,3); kxp = kx+1; delt = dt
   xkzm_m  = 1.0; xkzm_h  = 1.0; xkzm_s  = 1.0
   dspfac  = 1.0; dspheat = .true.
   ! bl_upfr = 0.13; bl_dnfr = 0.1
   bl_upfr = upd_fac; bl_dnfr = dnd_fac    ! ZNT 10/09/2020
   ! NOTE: Setting kinver = 0 disables background diffusivity.
   ! (diff. is still ~1.4e-4 m2/s, due to limiters of TKE and mixlen) 
   ! Setting kinver = kx enables backg. diff. for the entire column.
   kinver = 0 
   dku_out = 0.; dkt_out = 0.; dkq_out = 0.

   ! ZNT 06/22/2023: get aerosol following the algorithm in uw_conv.F90
   ! <TO DO #1: NEED TO PASS asol INTO HANEDMF > <- Done 06/27/2023
   do k = 1, nlev
      pmass(:,:,k) = (phalf(:,:,k+1) - phalf(:,:,k))/grav
   enddo

   if (present(asol)) then 
      if(use_online_aerosol) then
         do_mf_aer = 2
      else
         do_mf_aer = 1
      endif
      call get_aerosol (zhalf, pmass, asol,         &
                        am1, am2, am3, am4, am5 ,   &
                        amx1, amx2, amx3, amx4, amx5 )
   else 
      do_mf_aer = 0
      am1=0.; am2=0.; am3=0.; am4=0.; am5=0.;
      amx1=0.; amx2=0.; amx3=0.; amx4=0.; amx5=0.;
   endif
   
   ! ZNT 06/27/2023: debug output all aerosol values
   if (do_debug .and. do_print_aer < 2) then
      write(*,*) 'do_mf_aer', do_mf_aer
      write(*,*) 'am1', am1
      write(*,*) 'am2', am2
      write(*,*) 'am3', am3
      write(*,*) 'am4', am4
      write(*,*) 'am5', am5
      write(*,*) 'amx1', amx1
      write(*,*) 'amx2', amx2
      write(*,*) 'amx3', amx3
      write(*,*) 'amx4', amx4
      write(*,*) 'amx5', amx5
      do_print_aer = do_print_aer + 1
   endif
   ! <TO DO #2: NEED TO PASS am1-4 and amx1-4 INTO satmedmfvdifq_run and then draft solvers>
   ! <TO DO #3: NEED TO COMPUTE ACTIVATION (follow 'do_new_qnact' and 'do_2nd_act')>


   do i=1,size(t,1)
     do j=1,size(t,2)

       ! ZNT 04/27/2020: need to convert INPUT fields: 
       ! u1,v1,t1,q1,swh,hlw,xmu,garea,psk,rbsoil,zorl,u10m,v10m,fm,fh,
       ! tsea,heat,evap,stress,spd1,prsi,del,prsl,prslk,phii,phil,delt,
       ! dspheat,kinver,xkzm_m,xkzm_h,xkzm_s,dspfac,bl_upfr,bl_dnfr,
       !
       ! Note that several fields are missing! e.g., swh,garea
       ! Need to change physics driver to pass in: u_ref, v_ref, t_surf_rad

       if (use_z_above_sfc) then 
       ! ZNT 03/08/2023: Use geopot above surface. This should be the correct implementation
           phii (1,:) =  zhalf_ag(i,j,kxp:1:-1)*grav
           phil (1,:) =  zfull_ag(i,j,kx:1:-1)*grav
       else 
       ! ZNT 03/08/2023: Use geopot above sea level. Used in 2021-22 but incorrect.
           phii (1,:) =  zhalf(i,j,kxp:1:-1)*grav
           phil (1,:) =  zfull(i,j,kx:1:-1)*grav
       endif

       u1(1,:)   = u(i,j,kx:1:-1) 
       v1(1,:)   = v(i,j,kx:1:-1)
       t1(1,:)   = t(i,j,kx:1:-1)
       ! tracer order: sphum, cld_liq, cld_ice, cld_amt, tke
       q1(1,:,1) = qv(i,j,kx:1:-1)
       q1(1,:,2) = ql(i,j,kx:1:-1)
       q1(1,:,3) = qi(i,j,kx:1:-1)
       q1(1,:,4) = qa(i,j,kx:1:-1)
       q1(1,:,5) = qn(i,j,kx:1:-1)        ! ZNT 05/03/2023: adding qn
       q1(1,:,6) = qni(i,j,kx:1:-1)       ! ZNT 05/03/2023: adding qni
       q1(1,:,7) = t(i,j,kx:1:-1) + &
                    phil(1,:)/cp_air      ! ZNT 05/17/2023: adding passive T+gz/cp
       q1(1,:,8) = qv(i,j,kx:1:-1)        ! ZNT 05/17/2023: adding passive qv
       q1(1,:,9) = ql(i,j,kx:1:-1)        ! ZNT 05/17/2023: adding passive ql
       q1(1,:,10) = qi(i,j,kx:1:-1)       ! ZNT 05/17/2023: adding passive qi
       q1(1,:,11) = qa(i,j,kx:1:-1)       ! ZNT 07/05/2023: adding passive qa
       q1(1,:,12) = qn(i,j,kx:1:-1)       ! ZNT 07/05/2023: adding passive qn
       q1(1,:,13) = qni(i,j,kx:1:-1)      ! ZNT 07/05/2023: adding passive qni
       q1(1,:,14) = tke(i,j,kx:1:-1)      ! ZNT 05/22/2020: TKE as tracer now

       ! ZNT 06/28/2023: Pass am1-4 and amx1-4 to aer1 and aerx1 arrays as input for satmedmfvdifq_run
       aer1(1,:,1)  = am1(i,j,kx:1:-1)
       aer1(1,:,2)  = am2(i,j,kx:1:-1)
       aer1(1,:,3)  = am3(i,j,kx:1:-1)
       aer1(1,:,4)  = am4(i,j,kx:1:-1)
       aerx1(1,:,1) = amx1(i,j,kx:1:-1)
       aerx1(1,:,2) = amx2(i,j,kx:1:-1)
       aerx1(1,:,3) = amx3(i,j,kx:1:-1)
       aerx1(1,:,4) = amx4(i,j,kx:1:-1)
       

       swh  (1,:) = 0.0
       hlw  (1,:) = tdtlw_in(i,j,kx:1:-1)
       xmu  (1)   = 1.0
       ! For now, assume grid size = (100km)^2, which
       ! effectively disables resol. dependence
       garea(1)   = 100.e3*100.e3  
       ! For now, zero out vegetation
       zvfun(1)   = 0.0

       psk  (1)   = (phalf(i,j,kxp)/1.0e5)**kappa 
       prsi (1,:) =  phalf(i,j,kxp:1:-1)
       del  (1,:) =  phalf(i,j,kxp:2:-1) - phalf(i,j,kx:1:-1)
       prsl (1,:) =  pfull(i,j,kx:1:-1)
       prslk(1,:) = (pfull(i,j,kx:1:-1)/1.0e5)**kappa  

       u10m (1) = u_ref(i,j)
       v10m (1) = v_ref(i,j)
       zorl (1) = rough(i,j)*100.0
       tsea (1) = t_surf_rad(i,j)

       ! From Monin-Obukhov scheme: (need to add as an output from SFCFlx scheme)
       ! delb = grav*(pt0 - pt)/pt0; rbsoil = -z*delta_b/(speed*speed + small)
       rho_a = -1.0/grav*(phalf(i,j,kxp) - phalf(i,j,kx))/  &
                         (zhalf(i,j,kxp) - zhalf(i,j,kx))
       delb = grav*(thv_surf(i,j) - thv_atm(i,j))/thv_surf(i,j)
       ! heat (1) = shflx(i,j)/rho_a/Cp_Air   ! in K*m/s
       heat (1) = min(max(shflx(i,j),-shflx_max),shflx_max)/rho_a/Cp_Air  ! ZNT 10/14/2020: added limiter
       evap (1) = lhflx(i,j)/rho_a          ! in kg/kg*m/s
       stress(1) = u_star(i,j)**2.0         ! in m2/s2
       ! this differs from (tau_x^2+tau_y^2)^.5 by gust and density

       ! spd1(1) = wind(i,j) ! Note wind(i,j) is not equal to ((uatm)^2+(vatm)^2), need to correct
       spd1(1) = max(sqrt(u1(1,1)**2.0 + v1(1,1)**2.0), 1.e-02)
       fm(1) = vonkarm*wind(i,j)/u_star(i,j)
       fh(1) = vonkarm*delb/b_star(i,j)      ! ZNT 03/09/2023: Note - unintended results when b_star \approx 0?

       if (use_z_above_sfc) then 
       ! ZNT 03/08/2023: Use height above surface. This should be the correct implementation
           rbsoil(1) = -zfull_ag(i,j,kx)*delb/max(wind(i,j)**2.0,1.e-04)
       else 
       ! ZNT 03/08/2023: Use height above sea level. Used in 2021-22 but incorrect.
           rbsoil(1) = -zfull(i,j,kx)*delb/max(wind(i,j)**2.0,1.e-04)
       endif

       ! DEBUG - 05/03/2020

       if (do_debug) then
          icount = icount + 1
          write(*,*) ' '
          write(*,*) '******** ICOUNT = ', icount, ' ********'
          write(*,*) '********** INPUT **********'
          write(*,*) 'kx', kx
          write(*,*) 'grav', grav, 'rdgas', rdgas, 'cp_air', cp_air, &
                     'hlv', hlv, 'hlf', hlf, 'fv', fv, &
                     'eps', eps, 'epsm1', epsm1
          write(*,*) 'u1', u1
          write(*,*) 'v1', v1
          write(*,*) 't1', t1
          write(*,*) 'qv=q1(:,1)', q1(:,:,1)
          write(*,*) 'ql=q1(:,2)', q1(:,:,2)
          write(*,*) 'qi=q1(:,3)', q1(:,:,3)
          write(*,*) 'qa=q1(:,4)', q1(:,:,4)
          write(*,*) 'qn=q1(:,5)', q1(:,:,5)   ! ZNT 05/03/2023
          write(*,*) 'qni=q1(:,6)', q1(:,:,6)  ! ZNT 05/03/2023
          write(*,*) 'T+gz/cp=q1(:,7)', q1(:,:,7)  ! ZNT 05/17/2023
          write(*,*) 'qv=q1(:,8)', q1(:,:,8)   ! ZNT 05/17/2023
          write(*,*) 'ql=q1(:,9)', q1(:,:,9)   ! ZNT 05/17/2023
          write(*,*) 'qi=q1(:,10)', q1(:,:,10) ! ZNT 05/17/2023
          write(*,*) 'qa=q1(:,11)', q1(:,:,11)  ! ZNT 07/05/2023
          write(*,*) 'qn=q1(:,12)', q1(:,:,12)  ! ZNT 07/05/2023
          write(*,*) 'qni=q1(:,13)', q1(:,:,13) ! ZNT 07/05/2023
          write(*,*) 'tke=q1(:,14)', q1(:,:,14)
          ! ZNT 06/28/2023: aerosol input
          write(*,*) 'am1=aer1(:,1)', aer1(:,:,1)
          write(*,*) 'am2=aer1(:,2)', aer1(:,:,2)
          write(*,*) 'am3=aer1(:,3)', aer1(:,:,3)
          write(*,*) 'am4=aer1(:,4)', aer1(:,:,4)
          write(*,*) 'amx1=aerx1(:,1)', aerx1(:,:,1)
          write(*,*) 'amx2=aerx1(:,2)', aerx1(:,:,2)
          write(*,*) 'amx3=aerx1(:,3)', aerx1(:,:,3)
          write(*,*) 'amx4=aerx1(:,4)', aerx1(:,:,4)

          write(*,*) 'swh', swh
          write(*,*) 'hlw', hlw
          write(*,*) 'xmu', xmu
          write(*,*) 'garea', garea
          write(*,*) 'zvfun', zvfun
          write(*,*) 'psk', psk
          write(*,*) 'rbsoil', rbsoil
          write(*,*) 'zorl', zorl
          write(*,*) 'u10m', u10m
          write(*,*) 'v10m', v10m
          write(*,*) 'fm', fm
          write(*,*) 'fh', fh
          write(*,*) 'tsea', tsea
          write(*,*) 'heat', heat
          write(*,*) 'evap', evap
          write(*,*) 'stress', stress
          write(*,*) 'spd1', spd1
          write(*,*) 'prsi', prsi
          write(*,*) 'del', del
          write(*,*) 'prsl', prsl
          write(*,*) 'prslk', prslk
          write(*,*) 'phii', phii
          write(*,*) 'phil', phil
          write(*,*) 'delt', delt
          write(*,*) 'dspheat', dspheat
          write(*,*) 'do_mf_aer', do_mf_aer
          write(*,*) 'kinver', kinver
          write(*,*) 'xkzm_m', xkzm_m
          write(*,*) 'xkzm_h', xkzm_h
          write(*,*) 'xkzm_s', xkzm_s
          write(*,*) 'dspfac', dspfac
          write(*,*) 'bl_upfr', bl_upfr
          write(*,*) 'bl_dnfr', bl_dnfr
       endif

       if (do_tke_eq) then
          nstep_tke = tke_eq_iter   ! ZNT: do_tke_eq: zero out TKE and iterate
          ! q1(1,:,7) = 0.0           ! ZNT 05/03/2023: changing TKE index from 5 to 7
          ! q1(1,:,11) = 0.0           ! ZNT 05/17/2023: changing TKE index from 7 to 11
          q1(1,:,14) = 0.0           ! ZNT 07/05/2023: changing TKE index from 11 to 14
       else
          nstep_tke = 1             ! ZNT: not do_tke_eq: simple forward stepping
                                ! (Or do RK2?)
       endif

       do istep_tke = 1,nstep_tke
          ! Nullify total (large-step) output arrays
          dv = 0.0; du = 0.0; tdt = 0.0; rtg = 0.0
          dvsfc = 0.0; dusfc = 0.0
          dtsfc = 0.0; dqsfc = 0.0
          dku = 0.0; dkt = 0.0; dkq = 0.0
          xmf = 0.0; xmfd = 0.0
          hpbl = 0.0; kpbl = 0
          dv_mf = 0.0; du_mf = 0.0
          tdt_mf = 0.0; qdt_mf = 0.0
          qldt_mf = 0.0; qidt_mf = 0.0; qadt_mf = 0.0
          qndt_mf = 0.0; qnidt_mf = 0.0; tqpdt_mf = 0.0     ! ZNT 05/03/2023; 05/17/2023
          dv_mfd = 0.0; du_mfd = 0.0
          tdt_mfd = 0.0; qdt_mfd = 0.0
          qldt_mfd = 0.0; qidt_mfd = 0.0; qadt_mfd = 0.0
          qndt_mfd = 0.0; qnidt_mfd = 0.0; tqpdt_mfd = 0.0  ! ZNT 05/03/2023; 05/17/2023
          dusfc_mf = 0.0; dvsfc_mf = 0.0
          dtsfc_mf = 0.0; dqsfc_mf = 0.0
          dusfc_mfd = 0.0; dvsfc_mfd = 0.0
          dtsfc_mfd = 0.0; dqsfc_mfd = 0.0
          
          ! ZNT 05/24/2023: subs-detr tendencies
          tdt_sub = 0.0;   qdt_sub = 0.0;   qldt_sub = 0.0;   qidt_sub = 0.0
          qadt_sub = 0.0;  qndt_sub = 0.0;  qnidt_sub = 0.0
          tdt_subd = 0.0;  qdt_subd = 0.0;  qldt_subd = 0.0;  qidt_subd = 0.0
          qadt_subd = 0.0; qndt_subd = 0.0; qnidt_subd = 0.0
          tdt_dtr = 0.0;   qdt_dtr = 0.0;   qldt_dtr = 0.0;   qidt_dtr = 0.0
          qadt_dtr = 0.0;  qndt_dtr = 0.0;  qnidt_dtr = 0.0
          tdt_dtrd = 0.0;  qdt_dtrd = 0.0;  qldt_dtrd = 0.0;  qidt_dtrd = 0.0
          qadt_dtrd = 0.0; qndt_dtrd = 0.0; qnidt_dtrd = 0.0
          
          tqpdt_sub = 0.0; tqpdt_subd = 0.0   ! ZNT 05/23/2023
          tqpdt_dtr = 0.0; tqpdt_dtrd = 0.0   ! ZNT 05/23/2023
          
          de_buop = 0.0; de_shrp = 0.0; de_diss = 0.0          

          ! Initialize single-step input fields, ZNT: 10/02/2020
          u1_1 = u1; v1_1 = v1; t1_1 = t1; q1_1 = q1
          do istep_edmf = 1,nstep_edmf
             ! Nullify single-step output arrays, ZNT: 10/02/2020
             delt_1 = delt/nstep_edmf
             dv_1 = 0.0; du_1 = 0.0; tdt_1 = 0.0; rtg_1 = 0.0
             dvsfc_1 = 0.0; dusfc_1 = 0.0
             dtsfc_1 = 0.0; dqsfc_1 = 0.0
             dku_1 = 0.0; dkt_1 = 0.0; dkq_1 = 0.0
             xmf_1 = 0.0; xmfd_1 = 0.0
             hpbl_1 = 0.0; kpbl_1 = 0
             dv_mf_1 = 0.0; du_mf_1 = 0.0
             tdt_mf_1 = 0.0; qdt_mf_1 = 0.0
             qldt_mf_1 = 0.0; qidt_mf_1 = 0.0; qadt_mf_1 = 0.0
             qndt_mf_1 = 0.0; qnidt_mf_1 = 0.0; tqpdt_mf_1 = 0.0    ! ZNT 05/03/2023; 05/17/2023
             dv_mfd_1 = 0.0; du_mfd_1 = 0.0
             tdt_mfd_1 = 0.0; qdt_mfd_1 = 0.0
             qldt_mfd_1 = 0.0; qidt_mfd_1 = 0.0; qadt_mfd_1 = 0.0
             qndt_mfd_1 = 0.0; qnidt_mfd_1 = 0.0; tqpdt_mfd_1 = 0.0 ! ZNT 05/03/2023; 05/17/2023
             dusfc_mf_1 = 0.0; dvsfc_mf_1 = 0.0
             dtsfc_mf_1 = 0.0; dqsfc_mf_1 = 0.0
             dusfc_mfd_1 = 0.0; dvsfc_mfd_1 = 0.0
             dtsfc_mfd_1 = 0.0; dqsfc_mfd_1 = 0.0
             
             ! ZNT 05/24/2023: subs-detr tendencies
             tdt_sub_1 = 0.0;   qdt_sub_1 = 0.0;   qldt_sub_1 = 0.0;   qidt_sub_1 = 0.0
             qadt_sub_1 = 0.0;  qndt_sub_1 = 0.0;  qnidt_sub_1 = 0.0
             tdt_subd_1 = 0.0;  qdt_subd_1 = 0.0;  qldt_subd_1 = 0.0;  qidt_subd_1 = 0.0
             qadt_subd_1 = 0.0; qndt_subd_1 = 0.0; qnidt_subd_1 = 0.0
             tdt_dtr_1 = 0.0;   qdt_dtr_1= 0.0;    qldt_dtr_1 = 0.0;   qidt_dtr_1 = 0.0
             qadt_dtr_1 = 0.0;  qndt_dtr_1 = 0.0;  qnidt_dtr_1 = 0.0
             tdt_dtrd_1 = 0.0;  qdt_dtrd_1 = 0.0;  qldt_dtrd_1 = 0.0;  qidt_dtrd_1 = 0.0
             qadt_dtrd_1 = 0.0; qndt_dtrd_1 = 0.0; qnidt_dtrd_1 = 0.0
          
             tqpdt_sub_1 = 0.0; tqpdt_subd_1 = 0.0   ! ZNT 05/23/2023
             tqpdt_dtr_1 = 0.0; tqpdt_dtrd_1 = 0.0   ! ZNT 05/23/2023
             
             de_buop_1 = 0.0; de_shrp_1 = 0.0; de_diss_1 = 0.0

             ! call satmedmfvdifq_run(1,1,kx,7,2,3,4,5,6,7,             & ! ZNT: 05/12/2023, added indices for qn/qni
             ! call satmedmfvdifq_run(1,1,kx,11,2,3,4,5,6,7,11,         & ! ZNT: 05/17/2023, added indices for passive T,q,ql,qi
             call satmedmfvdifq_run(1,1,kx,14,2,3,4,5,6,7,14,         & ! ZNT: 07/05/2023, added indices for passive qa,qn,qni
                  grav,rdgas,Cp_Air,rvgas,HLv,HLf,fv,eps,epsm1,       &
                  dv_1,du_1,tdt_1,rtg_1,u1_1,v1_1,t1_1,q1_1,          & 
                  aer1,aerx1,swh,hlw,xmu,garea,zvfun,                 & ! ZNT: 06/28/2023, add aer1 and aerx1
                  cfl_crit,do_cfl_col,                                & ! ZNT: 09/24-28/2020; 06/06/2021
                  ck0,ck1,ch0,ch1,ce0,do_cmp_qs,use_hl,do_mf_aer,     & ! ZNT: 09/24-28/2020; 06/06/2021; 06/12/2023; 06/28/2023 do_mf_aer
                  tkeh_flg,cldtime, do_upd_ovrsht,c_entu,c_entd,actei,& ! ZNT: 11/11-18/2020; 06/21/2021
                  do_tkevdif, mf_tkeflx_fac, do_tkediss_imp,          & ! ZNT: 06/05-06/2021 
                  mf_buop_fac, mf_shrp_fac, do_pos_mf_buop,           & ! ZNT: 06/05/2021
                  do_km_L65, do_wush, do_tkemean,                     & ! ZNT: 10/03/2023; 12/24/2023
                  psk,rbsoil,zorl,u10m,v10m,fm,fh,                    &
                  tsea,heat,evap,stress,spd1,kpbl_1,                  &
                  prsi,del,prsl,prslk,phii,phil,delt_1,               &
                  dspheat,dusfc_1,dvsfc_1,dtsfc_1,dqsfc_1,hpbl_1,     &
                  kinver,xkzm_m,xkzm_h,xkzm_s,dspfac,bl_upfr,bl_dnfr, &
                  dku_1, dkt_1, dkq_1, xmf_1, xmfd_1,                 &
                  tcko,qvcko,qlcko,qicko,qacko,qncko,qnicko,          & ! ZNT: 06/19/2020; qn/qni added 05/03/2023; qa added 05/24/2023
                  tqpcko,buou,                                        & ! ZNT: 05/17/2023: added passive T,q,ql,qi 
                  tcdo,qvcdo,qlcdo,qicdo,qacdo,qncdo,qnicdo,          & ! ZNT: 06/19/2020; qn/qni added 05/03/2023; qa added 05/24/2023
                  tqpcdo,buod,                                        & ! ZNT: 05/17/2023: added passive T,q,ql,qi 
                  dv_mf_1,du_mf_1,tdt_mf_1,qdt_mf_1,                  & ! ZNT: 05/08/2020, MF diag
                  qldt_mf_1,qidt_mf_1,qadt_mf_1,                      & ! ZNT: 10/11/2020, MF diag
                  qndt_mf_1,qnidt_mf_1,tqpdt_mf_1,                    & ! ZNT: 05/03/2023; 05/17/2023
                  dv_mfd_1,du_mfd_1,tdt_mfd_1,qdt_mfd_1,              & ! ZNT: 05/08/2020, MF diag
                  qldt_mfd_1,qidt_mfd_1,qadt_mfd_1,                   & ! ZNT: 10/11/2020, MF diag
                  qndt_mfd_1,qnidt_mfd_1,tqpdt_mfd_1,                 & ! ZNT: 05/03/2023; 05/17/2023
                  dusfc_mf_1,dvsfc_mf_1,dtsfc_mf_1,dqsfc_mf_1,        & ! ZNT: 05/08/2020, MF diag
                  dusfc_mfd_1,dvsfc_mfd_1,dtsfc_mfd_1,dqsfc_mfd_1,    & ! ZNT: 05/08/2020, MF diag
                  tdt_sub_1,qdt_sub_1,qldt_sub_1,qidt_sub_1,          & ! ZNT: 05/24/2023, subs-detr
                  qadt_sub_1,qndt_sub_1,qnidt_sub_1,                  & ! ZNT: 05/24/2023
                  tdt_subd_1,qdt_subd_1,qldt_subd_1,qidt_subd_1,      & ! ZNT: 05/24/2023
                  qadt_subd_1,qndt_subd_1,qnidt_subd_1,               & ! ZNT: 05/24/2023
                  tdt_dtr_1,qdt_dtr_1,qldt_dtr_1,qidt_dtr_1,          & ! ZNT: 05/24/2023, subs-detr
                  qadt_dtr_1,qndt_dtr_1,qnidt_dtr_1,                  & ! ZNT: 05/24/2023
                  tdt_dtrd_1,qdt_dtrd_1,qldt_dtrd_1,qidt_dtrd_1,      & ! ZNT: 05/24/2023
                  qadt_dtrd_1,qndt_dtrd_1,qnidt_dtrd_1,               & ! ZNT: 05/24/2023
                  tqpdt_sub_1,tqpdt_subd_1,tqpdt_dtr_1,tqpdt_dtrd_1,  & ! ZNT: 05/23/2023
                  de_buop_1,de_shrp_1,de_diss_1,                      & ! ZNT: 10/03/2020, TKE diag
                  errmsg,errflg)

              ! Update single-step input fields, ZNT: 10/02/2020
              u1_1 = u1_1 + du_1*delt_1;  v1_1 = v1_1 + du_1*delt_1
              t1_1 = t1_1 + tdt_1*delt_1; q1_1 = q1_1 + rtg_1*delt_1

              ! Accumulate total (large-step) output arrays, ZNT: 10/02/2020 
              dv = dv + dv_1/nstep_edmf
              du = du + du_1/nstep_edmf
              tdt = tdt + tdt_1/nstep_edmf
              rtg = rtg + rtg_1/nstep_edmf
              dvsfc = dvsfc + dvsfc_1/nstep_edmf
              dusfc = dusfc + dusfc_1/nstep_edmf
              dtsfc = dtsfc + dtsfc_1/nstep_edmf
              dqsfc = dqsfc + dqsfc_1/nstep_edmf
              dku = dku + dku_1/nstep_edmf
              dkt = dkt + dkt_1/nstep_edmf
              dkq = dkq + dkq_1/nstep_edmf
              xmf = xmf + xmf_1/nstep_edmf
              xmfd = xmfd + xmfd_1/nstep_edmf

              dv_mf = dv_mf + dv_mf_1/nstep_edmf
              du_mf = du_mf + du_mf_1/nstep_edmf
              tdt_mf = tdt_mf + tdt_mf_1/nstep_edmf
              qdt_mf = qdt_mf + qdt_mf_1/nstep_edmf
              qldt_mf = qldt_mf + qldt_mf_1/nstep_edmf
              qidt_mf = qidt_mf + qidt_mf_1/nstep_edmf
              qadt_mf = qadt_mf + qadt_mf_1/nstep_edmf
              qndt_mf = qndt_mf + qndt_mf_1/nstep_edmf     ! ZNT 05/03/2023
              qnidt_mf = qnidt_mf + qnidt_mf_1/nstep_edmf  ! ZNT 05/03/2023
              tqpdt_mf = tqpdt_mf + tqpdt_mf_1/nstep_edmf  ! ZNT 05/17/2023
              dvsfc_mf = dvsfc_mf + dvsfc_mf_1/nstep_edmf
              dusfc_mf = dusfc_mf + dusfc_mf_1/nstep_edmf
              dtsfc_mf = dtsfc_mf + dtsfc_mf_1/nstep_edmf
              dqsfc_mf = dqsfc_mf + dqsfc_mf_1/nstep_edmf

              dv_mfd = dv_mfd + dv_mfd_1/nstep_edmf
              du_mfd = du_mfd + du_mfd_1/nstep_edmf
              tdt_mfd = tdt_mfd + tdt_mfd_1/nstep_edmf
              qdt_mfd = qdt_mfd + qdt_mfd_1/nstep_edmf
              qldt_mfd = qldt_mfd + qldt_mfd_1/nstep_edmf
              qidt_mfd = qidt_mfd + qidt_mfd_1/nstep_edmf
              qadt_mfd = qadt_mfd + qadt_mfd_1/nstep_edmf
              qndt_mfd = qndt_mfd + qndt_mfd_1/nstep_edmf     ! ZNT 05/03/2023
              qnidt_mfd = qnidt_mfd + qnidt_mfd_1/nstep_edmf  ! ZNT 05/03/2023
              tqpdt_mfd = tqpdt_mfd + tqpdt_mfd_1/nstep_edmf  ! ZNT 05/17/2023
              dvsfc_mfd = dvsfc_mfd + dvsfc_mfd_1/nstep_edmf
              dusfc_mfd = dusfc_mfd + dusfc_mfd_1/nstep_edmf
              dtsfc_mfd = dtsfc_mfd + dtsfc_mfd_1/nstep_edmf
              dqsfc_mfd = dqsfc_mfd + dqsfc_mfd_1/nstep_edmf

              ! ZNT 05/24/2023: subs-detr tendencies
              tdt_sub = tdt_sub + tdt_sub_1/nstep_edmf
              qdt_sub = qdt_sub + qdt_sub_1/nstep_edmf
              qldt_sub = qldt_sub + qldt_sub_1/nstep_edmf
              qidt_sub = qidt_sub + qidt_sub_1/nstep_edmf
              qadt_sub = qadt_sub + qadt_sub_1/nstep_edmf
              qndt_sub = qndt_sub + qndt_sub_1/nstep_edmf
              qnidt_sub = qnidt_sub + qnidt_sub_1/nstep_edmf
              tdt_subd = tdt_subd + tdt_subd_1/nstep_edmf
              qdt_subd = qdt_subd + qdt_subd_1/nstep_edmf
              qldt_subd = qldt_subd + qldt_subd_1/nstep_edmf
              qidt_subd = qidt_subd + qidt_subd_1/nstep_edmf
              qadt_subd = qadt_subd + qadt_subd_1/nstep_edmf
              qndt_subd = qndt_subd + qndt_subd_1/nstep_edmf
              qnidt_subd = qnidt_subd + qnidt_subd_1/nstep_edmf
              
              tdt_dtr = tdt_dtr + tdt_dtr_1/nstep_edmf
              qdt_dtr = qdt_dtr + qdt_dtr_1/nstep_edmf
              qldt_dtr = qldt_dtr + qldt_dtr_1/nstep_edmf
              qidt_dtr = qidt_dtr + qidt_dtr_1/nstep_edmf
              qadt_dtr = qadt_dtr + qadt_dtr_1/nstep_edmf
              qndt_dtr = qndt_dtr + qndt_dtr_1/nstep_edmf
              qnidt_dtr = qnidt_dtr + qnidt_dtr_1/nstep_edmf
              tdt_dtrd = tdt_dtrd + tdt_dtrd_1/nstep_edmf
              qdt_dtrd = qdt_dtrd + qdt_dtrd_1/nstep_edmf
              qldt_dtrd = qldt_dtrd + qldt_dtrd_1/nstep_edmf
              qidt_dtrd = qidt_dtrd + qidt_dtrd_1/nstep_edmf
              qadt_dtrd = qadt_dtrd + qadt_dtrd_1/nstep_edmf
              qndt_dtrd = qndt_dtrd + qndt_dtrd_1/nstep_edmf
              qnidt_dtrd = qnidt_dtrd + qnidt_dtrd_1/nstep_edmf
              
              tqpdt_sub  = tqpdt_sub  + tqpdt_sub_1/nstep_edmf  ! ZNT 05/23/2023
              tqpdt_subd = tqpdt_subd + tqpdt_subd_1/nstep_edmf ! ZNT 05/23/2023
              tqpdt_dtr  = tqpdt_dtr  + tqpdt_dtr_1/nstep_edmf  ! ZNT 05/23/2023
              tqpdt_dtrd = tqpdt_dtrd + tqpdt_dtrd_1/nstep_edmf ! ZNT 05/23/2023


              de_buop = de_buop + de_buop_1/nstep_edmf
              de_shrp = de_shrp + de_shrp_1/nstep_edmf
              de_diss = de_diss + de_diss_1/nstep_edmf

              ! Note: not averaging the updraft variables ('*ko','*do','buo*') yet
              ! Because they need mass-flux weighting (ZNT 10/02/2020)

              ! For PBL height, need to take the highest/largest index
              hpbl = max(hpbl, hpbl_1)
              kpbl = max(kpbl, kpbl_1)

          enddo
          ! q1(1,:,7) = q1(1,:,7) + rtg(1,:,7)*delt  ! ZNT 05/03/2023: changing TKE index from 5 to 7
          ! q1(1,:,11) = q1(1,:,11) + rtg(1,:,11)*delt  ! ZNT 05/17/2023: changing TKE index from 7 to 11
          q1(1,:,14) = q1(1,:,14) + rtg(1,:,14)*delt  ! ZNT 07/05/2023: changing TKE index from 11 to 14
       enddo

       if (do_debug) then
          write(*,*) '*********** OUTPUT **********'
          write(*,*) 'du/dt', du
          write(*,*) 'dv/dt', dv
          write(*,*) 'dT/dt', tdt
          write(*,*) 'dqv/dt', rtg(:,:,1)
          write(*,*) 'dql/dt', rtg(:,:,2)
          write(*,*) 'dqi/dt', rtg(:,:,3)
          write(*,*) 'dqa/dt', rtg(:,:,4)
          write(*,*) 'dqn/dt', rtg(:,:,5)  ! ZNT 05/03/2023
          write(*,*) 'dqni/dt', rtg(:,:,6) ! ZNT 05/03/2023
          write(*,*) 'd(T+gz/cp)/dt,passive', rtg(:,:,7)  ! ZNT 05/17/2023
          write(*,*) 'dqv/dt,passive', rtg(:,:,8)  ! ZNT 05/17/2023
          write(*,*) 'dql/dt,passive', rtg(:,:,9)  ! ZNT 05/17/2023
          write(*,*) 'dqi/dt,passive', rtg(:,:,10) ! ZNT 05/17/2023
          write(*,*) 'dqa/dt,passive', rtg(:,:,11)  ! ZNT 07/05/2023
          write(*,*) 'dqn/dt,passive', rtg(:,:,12)  ! ZNT 07/05/2023
          write(*,*) 'dqni/dt,passive', rtg(:,:,13) ! ZNT 07/05/2023
          write(*,*) 'dtke/dt', rtg(:,:,14)
          write(*,*) 'SHF', dtsfc
          write(*,*) 'LHF', dqsfc
          write(*,*) 'dragU', dusfc
          write(*,*) 'dragV', dvsfc
          write(*,*) 'hpbl', hpbl
          write(*,*) 'kpbl', kpbl
          write(*,*) 'Km', dku
          write(*,*) 'Kt', dkt
          write(*,*) 'Ke', dkq
          write(*,*) 'Mu', xmf
          write(*,*) 'Md', xmfd 
          write(*,*) 'buou', buou
          write(*,*) 'buod', buod      
       endif

       ! ZNT 04/27/2020: need to convert OUTPUT fields: 
       ! dv,du,tdt,rtg,dusfc,dvsfc,dtsfc,dqsfc,hpbl,dku,dkt,dkq,xmf,xmfd

       udt_out(i,j,kx:1:-1)  = du(1,:)
       vdt_out(i,j,kx:1:-1)  = dv(1,:)
       tdt_out(i,j,kx:1:-1)  = tdt(1,:)
       ! tracer order: sphum, cld_liq, cld_ice, cld_amt, tke
       qvdt_out(i,j,kx:1:-1) = rtg(1,:,1)
       qldt_out(i,j,kx:1:-1) = rtg(1,:,2)
       qidt_out(i,j,kx:1:-1) = rtg(1,:,3)
       qadt_out(i,j,kx:1:-1) = rtg(1,:,4)
       qndt_out(i,j,kx:1:-1)  = rtg(1,:,5)  ! ZNT 05/03/2023: adding qn
       qnidt_out(i,j,kx:1:-1) = rtg(1,:,6)  ! ZNT 05/03/2023: adding qni
       tke(i,j,kx:1:-1) = q1(1,:,14)      ! ZNT 05/22/2020: TKE is tracer now
       dku_out(i,j,kx:2:-1)  = dku(1,:)  ! ZNT 05/25/2020: Dimension corrected
       dkt_out(i,j,kx:2:-1)  = dkt(1,:)  ! ZNT (as above)
       dkq_out(i,j,kx:2:-1)  = dkq(1,:)  ! ZNT (as above)
       xmf_out(i,j,kx:1:-1)  = xmf(1,:)
       xmfd_out(i,j,kx:1:-1) = xmfd(1,:)
       tcko_out(i,j,kx:1:-1)  = tcko(1,:)   ! ZNT: 06/19/2020
       qvcko_out(i,j,kx:1:-1) = qvcko(1,:)  ! ZNT: 06/19/2020
       qlcko_out(i,j,kx:1:-1) = qlcko(1,:)  ! ZNT: 06/19/2020
       qicko_out(i,j,kx:1:-1) = qicko(1,:)  ! ZNT: 06/19/2020
       qacko_out(i,j,kx:1:-1) = qacko(1,:)   ! ZNT: 05/24/2023: adding qa
       qncko_out(i,j,kx:1:-1) = qncko(1,:)   ! ZNT: 05/03/2023: adding qn
       qnicko_out(i,j,kx:1:-1) = qnicko(1,:) ! ZNT: 05/03/2023: adding qni
       buou_out(i,j,kx:1:-1)  = buou(1,:)   ! ZNT: 06/19/2020
       tcdo_out(i,j,kx:1:-1)  = tcdo(1,:)   ! ZNT: 06/19/2020
       qvcdo_out(i,j,kx:1:-1) = qvcdo(1,:)  ! ZNT: 06/19/2020
       qlcdo_out(i,j,kx:1:-1) = qlcdo(1,:)  ! ZNT: 06/19/2020
       qicdo_out(i,j,kx:1:-1) = qicdo(1,:)  ! ZNT: 06/19/2020
       qacdo_out(i,j,kx:1:-1) = qacdo(1,:)   ! ZNT: 05/24/2023: adding qa
       qncdo_out(i,j,kx:1:-1) = qncdo(1,:)   ! ZNT: 05/03/2023: adding qn
       qnicdo_out(i,j,kx:1:-1) = qnicdo(1,:) ! ZNT: 05/03/2023: adding qni
       buod_out(i,j,kx:1:-1)  = buod(1,:)   ! ZNT: 06/19/2020

       dusfc_out(i,j) = dusfc(1)
       dvsfc_out(i,j) = dvsfc(1)
       dtsfc_out(i,j) = dtsfc(1)
       dqsfc_out(i,j) = dqsfc(1)
       hpbl_out(i,j) = hpbl(1)
       fm_out(i,j) = fm(1)
       ft_out(i,j) = fh(1)
       risfc_out(i,j) = rbsoil(1)
! ZNT 05/08/2020, 10/11/2020: Diagnostic MF output
       udt_mf_out(i,j,kx:1:-1)   = du_mf(1,:)
       vdt_mf_out(i,j,kx:1:-1)   = dv_mf(1,:)
       tdt_mf_out(i,j,kx:1:-1)   = tdt_mf(1,:)
       qvdt_mf_out(i,j,kx:1:-1)  = qdt_mf(1,:)
       qldt_mf_out(i,j,kx:1:-1)  = qldt_mf(1,:)
       qidt_mf_out(i,j,kx:1:-1)  = qidt_mf(1,:)
       qadt_mf_out(i,j,kx:1:-1)  = qadt_mf(1,:)
       qndt_mf_out(i,j,kx:1:-1)  = qndt_mf(1,:)    ! ZNT 05/03/2023: adding qn
       qnidt_mf_out(i,j,kx:1:-1)  = qnidt_mf(1,:)  ! ZNT 05/03/2023: adding qni

       udt_mfd_out(i,j,kx:1:-1)  = du_mfd(1,:)
       vdt_mfd_out(i,j,kx:1:-1)  = dv_mfd(1,:)
       tdt_mfd_out(i,j,kx:1:-1)  = tdt_mfd(1,:)
       qvdt_mfd_out(i,j,kx:1:-1) = qdt_mfd(1,:)
       qldt_mfd_out(i,j,kx:1:-1) = qldt_mfd(1,:)
       qidt_mfd_out(i,j,kx:1:-1) = qidt_mfd(1,:)
       qadt_mfd_out(i,j,kx:1:-1) = qadt_mfd(1,:)
       qndt_mfd_out(i,j,kx:1:-1) = qndt_mfd(1,:)   ! ZNT 05/03/2023: adding qn
       qnidt_mfd_out(i,j,kx:1:-1) = qnidt_mfd(1,:) ! ZNT 05/03/2023: adding qni

       ! ZNT 05/24/2023: subs-detr tendencies
       tdt_sub_out(i,j,kx:1:-1)   = tdt_sub(1,:)
       qvdt_sub_out(i,j,kx:1:-1)  = qdt_sub(1,:)
       qldt_sub_out(i,j,kx:1:-1)  = qldt_sub(1,:)
       qidt_sub_out(i,j,kx:1:-1)  = qidt_sub(1,:)
       qadt_sub_out(i,j,kx:1:-1)  = qadt_sub(1,:)
       qndt_sub_out(i,j,kx:1:-1)  = qndt_sub(1,:)
       qnidt_sub_out(i,j,kx:1:-1)  = qnidt_sub(1,:)

       tdt_subd_out(i,j,kx:1:-1)  = tdt_subd(1,:)
       qvdt_subd_out(i,j,kx:1:-1) = qdt_subd(1,:)
       qldt_subd_out(i,j,kx:1:-1) = qldt_subd(1,:)
       qidt_subd_out(i,j,kx:1:-1) = qidt_subd(1,:)
       qadt_subd_out(i,j,kx:1:-1) = qadt_subd(1,:)
       qndt_subd_out(i,j,kx:1:-1) = qndt_subd(1,:)
       qnidt_subd_out(i,j,kx:1:-1) = qnidt_subd(1,:)
       
       tdt_dtr_out(i,j,kx:1:-1)   = tdt_dtr(1,:)
       qvdt_dtr_out(i,j,kx:1:-1)  = qdt_dtr(1,:)
       qldt_dtr_out(i,j,kx:1:-1)  = qldt_dtr(1,:)
       qidt_dtr_out(i,j,kx:1:-1)  = qidt_dtr(1,:)
       qadt_dtr_out(i,j,kx:1:-1)  = qadt_dtr(1,:)
       qndt_dtr_out(i,j,kx:1:-1)  = qndt_dtr(1,:)
       qnidt_dtr_out(i,j,kx:1:-1)  = qnidt_dtr(1,:)

       tdt_dtrd_out(i,j,kx:1:-1)  = tdt_dtrd(1,:)
       qvdt_dtrd_out(i,j,kx:1:-1) = qdt_dtrd(1,:)
       qldt_dtrd_out(i,j,kx:1:-1) = qldt_dtrd(1,:)
       qidt_dtrd_out(i,j,kx:1:-1) = qidt_dtrd(1,:)
       qadt_dtrd_out(i,j,kx:1:-1) = qadt_dtrd(1,:)
       qndt_dtrd_out(i,j,kx:1:-1) = qndt_dtrd(1,:)
       qnidt_dtrd_out(i,j,kx:1:-1) = qnidt_dtrd(1,:)
              
       do itqp = 1,7  ! ZNT 05/17/2023: adding passive T,q,ql,qi; 07/05/2023: adding passive qa,qn,qni
          tqpdt_out(i,j,kx:1:-1,itqp) = rtg(1,:,6+itqp)
          tqpcko_out(i,j,kx:1:-1,itqp) = tqpcko(1,:,itqp)
          tqpcdo_out(i,j,kx:1:-1,itqp) = tqpcdo(1,:,itqp)
          tqpdt_mf_out(i,j,kx:1:-1,itqp) = tqpdt_mf(1,:,itqp)
          tqpdt_mfd_out(i,j,kx:1:-1,itqp) = tqpdt_mfd(1,:,itqp)
          ! ZNT 05/23/2023: Subs-detr decomposition of MF tendencies
          tqpdt_sub_out(i,j,kx:1:-1,itqp) = tqpdt_sub(1,:,itqp)
          tqpdt_subd_out(i,j,kx:1:-1,itqp) = tqpdt_subd(1,:,itqp)
          tqpdt_dtr_out(i,j,kx:1:-1,itqp) = tqpdt_dtr(1,:,itqp)
          tqpdt_dtrd_out(i,j,kx:1:-1,itqp) = tqpdt_dtrd(1,:,itqp)
       end do

       dusfc_mf_out(i,j)  = dusfc_mf(1)
       dvsfc_mf_out(i,j)  = dvsfc_mf(1)
       dtsfc_mf_out(i,j)  = dtsfc_mf(1)
       dqsfc_mf_out(i,j)  = dqsfc_mf(1)
       dusfc_mfd_out(i,j) = dusfc_mfd(1)
       dvsfc_mfd_out(i,j) = dvsfc_mfd(1)
       dtsfc_mfd_out(i,j) = dtsfc_mfd(1)
       dqsfc_mfd_out(i,j) = dqsfc_mfd(1)
! End ZNT 05/08/2020, 10/11/2020
! ZNT 10/03/2020: Diagnostic TKE output
       de_buop_out(i,j,kx:1:-1) = de_buop(1,:)
       de_shrp_out(i,j,kx:1:-1) = de_shrp(1,:)
       de_diss_out(i,j,kx:1:-1) = de_diss(1,:)
! End ZNT 10/03/2020
     enddo
   enddo

   ! Generate Output to vert_turb_driver
   ! ZNT 06/05/2020, 10/11/2020: MF only
   ! qlia_tnd_fac = 0 to reproduce model behavior before 10/11/2020
   ! ZNT 05/17/2023: add option for passive T,q,ql,qi
   udt_mfud = udt_mf_out + udt_mfd_out
   vdt_mfud = vdt_mf_out + vdt_mfd_out
   if (do_tqp_mf) then ! ZNT 05/17/2023: use passive T,q,ql,qi tendencies
      tdt_mfud = tqpdt_mf_out(:,:,:,1) + tqpdt_mfd_out(:,:,:,1)
      qdt_mfud = tqpdt_mf_out(:,:,:,2) + tqpdt_mfd_out(:,:,:,2)
      qldt_mfud = (tqpdt_mf_out(:,:,:,3) + tqpdt_mfd_out(:,:,:,3))*qlia_tnd_fac
      qidt_mfud = (tqpdt_mf_out(:,:,:,4) + tqpdt_mfd_out(:,:,:,4))*qlia_tnd_fac
      ! ZNT 07/05/2023: Passive qa, qn, qni tendencies
      qadt_mfud = (tqpdt_mf_out(:,:,:,5) + tqpdt_mfd_out(:,:,:,5))*qlia_tnd_fac * qa_mf_fac
      qndt_mfud = (tqpdt_mf_out(:,:,:,6) + tqpdt_mfd_out(:,:,:,6))*qlia_tnd_fac * qn_mf_fac
      qnidt_mfud = (tqpdt_mf_out(:,:,:,7) + tqpdt_mfd_out(:,:,:,7))*qlia_tnd_fac * qni_mf_fac
      ! ZNT 07/05/2023: The formulas below reproduce model behavior before 05/2023
      ! qadt_mfud  = (qadt_mf_out + qadt_mfd_out)  *qlia_tnd_fac * qa_mf_fac
      ! qndt_mfud  = (qndt_mf_out + qndt_mfd_out)  *qlia_tnd_fac * qn_mf_fac
      ! qnidt_mfud = (qnidt_mf_out + qnidt_mfd_out)*qlia_tnd_fac * qni_mf_fac
      
   elseif (do_subdtr) then  ! ZNT 06/01/2023: use subsidence-detrainment form
      tdt_mfud = tdt_sub_out + tdt_subd_out + &
                 tdt_dtr_out + tdt_dtrd_out
      qdt_mfud = qvdt_sub_out + qvdt_subd_out + &
                 qvdt_dtr_out + qvdt_dtrd_out
      qldt_mfud = (qldt_sub_out + qldt_subd_out + &
                   qldt_dtr_out + qldt_dtrd_out)*qlia_tnd_fac
      qidt_mfud = (qidt_sub_out + qidt_subd_out + &
                   qidt_dtr_out + qidt_dtrd_out)*qlia_tnd_fac
      ! ZNT 05/12/2023: Output of qa, qn and qni tendencies due to MF
      qadt_mfud  = (qadt_sub_out + qadt_subd_out + &
                    qadt_dtr_out + qadt_dtrd_out)  *qlia_tnd_fac * qa_mf_fac
      qndt_mfud  = (qndt_sub_out + qndt_subd_out + &
                    qndt_dtr_out + qndt_dtrd_out)  *qlia_tnd_fac * qn_mf_fac
      qnidt_mfud = (qnidt_sub_out + qnidt_subd_out + &
                    qnidt_dtr_out + qnidt_dtrd_out)*qlia_tnd_fac * qni_mf_fac
                    
   else   ! ZNT 06/01/2023: old MF-divergence form (missing the source terms)
      tdt_mfud = tdt_mf_out + tdt_mfd_out
      qdt_mfud = qvdt_mf_out + qvdt_mfd_out
      qldt_mfud = (qldt_mf_out + qldt_mfd_out)*qlia_tnd_fac
      qidt_mfud = (qidt_mf_out + qidt_mfd_out)*qlia_tnd_fac
      ! ZNT 05/12/2023: Output of qa, qn and qni tendencies due to MF
      qadt_mfud  = (qadt_mf_out + qadt_mfd_out)  *qlia_tnd_fac * qa_mf_fac
      qndt_mfud  = (qndt_mf_out + qndt_mfd_out)  *qlia_tnd_fac * qn_mf_fac
      qnidt_mfud = (qnidt_mf_out + qnidt_mfd_out)*qlia_tnd_fac * qni_mf_fac
   endif



   ! ZNT 06/05/2020: Pass ED directly for use in vert_diff_driver
   if (do_switch_ed_bug) then ! ZNT 10/15/2020: buggy code
       diff_t = dku_out
       diff_m = dkt_out
   else                       ! This is the correct pass
       diff_t = dkt_out
       diff_m = dku_out
   endif


   ! ZNT 12/16/2023: Zero out the TKE-based ED above z_ed_top AND above convective PBL-top
   if (z_ed_top < 1.0e5) then     ! AM5 L65 model top is at 1 Pa ~ 80 km, thus no need to zero out ED
                                  ! if z_ed_top is higher than 100km = 1e5 m

       if (use_z_above_sfc) then  ! ZNT 03/08/2023: Use height above surface (correct implementation)
          do k = 1, nlev+1
             where (zhalf_ag(:,:,k) > z_ed_top .and. zhalf_ag(:,:,k) > hpbl_out(:,:) )
                diff_t(:,:,k) = 0.0
                diff_m(:,:,k) = 0.0
             endwhere
          enddo

       else  ! ZNT 03/08/2023: Use height above sea level (used in 2021-22 but incorrect)
          do k = 1, nlev+1
             where (zhalf_ag(:,:,k) > z_ed_top .and. zhalf(:,:,k) > hpbl_out(:,:) ) 
                diff_t(:,:,k) = 0.0
                diff_m(:,:,k) = 0.0
             endwhere
          enddo
       endif ! use_z_above_sfc
   endif ! z_ed_top


   ! Generate Diagnostic Output
   if ( id_zfull_han > 0 ) then  ! ZNT 05/20/2023
      used = send_data ( id_zfull_han, zfull, Time_next, is, js, 1)
   endif  
   if ( id_zhalf_han > 0 ) then  ! ZNT 05/20/2023
      used = send_data ( id_zhalf_han, zhalf, Time_next, is, js, 1)
   endif  
   if ( id_pfull_han > 0 ) then  ! ZNT 05/20/2023
      used = send_data ( id_pfull_han, pfull, Time_next, is, js, 1)
   endif  
   if ( id_phalf_han > 0 ) then  ! ZNT 05/20/2023
      used = send_data ( id_phalf_han, phalf, Time_next, is, js, 1)
   endif  
   
   if ( id_tdt_han > 0 ) then
      used = send_data ( id_tdt_han, tdt_out, Time_next, is, js, 1)
   endif

   if ( id_qdt_han > 0 ) then
      used = send_data ( id_qdt_han, qvdt_out, Time_next, is, js, 1)
   endif

   if ( id_qldt_han > 0 ) then
      used = send_data ( id_qldt_han, qldt_out, Time_next, is, js, 1)
   endif

   if ( id_qidt_han > 0 ) then
      used = send_data ( id_qidt_han, qidt_out, Time_next, is, js, 1)
   endif

   if ( id_qndt_han > 0 ) then  ! ZNT 05/12/2023
      used = send_data ( id_qndt_han, qndt_out, Time_next, is, js, 1)
   endif

   if ( id_qnidt_han > 0 ) then  ! ZNT 05/12/2023
      used = send_data ( id_qnidt_han, qnidt_out, Time_next, is, js, 1)
   endif

   if ( id_qadt_han > 0 ) then
      used = send_data ( id_qadt_han, qadt_out, Time_next, is, js, 1)
   endif

   if ( id_udt_han > 0 ) then
      used = send_data ( id_udt_han, udt_out, Time_next, is, js, 1)
   endif

   if ( id_vdt_han > 0 ) then
      used = send_data ( id_vdt_han, vdt_out, Time_next, is, js, 1)
   endif

   if ( id_tke_han > 0 ) then
      used = send_data ( id_tke_han, tke, Time_next, is, js, 1)
   endif

   if ( id_km_han > 0 ) then
      used = send_data ( id_km_han, dku_out, Time_next, is, js, 1)
   endif

   if ( id_kt_han > 0 ) then
      used = send_data ( id_kt_han, dkt_out, Time_next, is, js, 1)
   endif

   if ( id_ke_han > 0 ) then
      used = send_data ( id_ke_han, dkq_out, Time_next, is, js, 1)
   endif

   if ( id_mu_han > 0 ) then
      used = send_data ( id_mu_han, xmf_out, Time_next, is, js, 1)
   endif

   if ( id_md_han > 0 ) then
      used = send_data ( id_md_han, xmfd_out, Time_next, is, js, 1)
   endif

   ! ZNT 06/19/2020
   if ( id_tcko_han > 0 ) then
      used = send_data ( id_tcko_han, tcko_out, Time_next, is, js, 1)
   endif
   if ( id_qvcko_han > 0 ) then
      used = send_data ( id_qvcko_han, qvcko_out, Time_next, is, js, 1)
   endif
   if ( id_qlcko_han > 0 ) then
      used = send_data ( id_qlcko_han, qlcko_out, Time_next, is, js, 1)
   endif
   if ( id_qicko_han > 0 ) then
      used = send_data ( id_qicko_han, qicko_out, Time_next, is, js, 1)
   endif
   if ( id_qacko_han > 0 ) then  ! ZNT 05/24/2023
      used = send_data ( id_qacko_han, qacko_out, Time_next, is, js, 1)
   endif
   if ( id_qncko_han > 0 ) then  ! ZNT 05/12/2023
      used = send_data ( id_qncko_han, qncko_out, Time_next, is, js, 1)
   endif
   if ( id_qnicko_han > 0 ) then ! ZNT 05/12/2023
      used = send_data ( id_qnicko_han, qnicko_out, Time_next, is, js, 1)
   endif
   if ( id_buou_han > 0 ) then
      used = send_data ( id_buou_han, buou_out, Time_next, is, js, 1)
   endif

   if ( id_tcdo_han > 0 ) then
      used = send_data ( id_tcdo_han, tcdo_out, Time_next, is, js, 1)
   endif
   if ( id_qvcdo_han > 0 ) then
      used = send_data ( id_qvcdo_han, qvcdo_out, Time_next, is, js, 1)
   endif
   if ( id_qlcdo_han > 0 ) then
      used = send_data ( id_qlcdo_han, qlcdo_out, Time_next, is, js, 1)
   endif
   if ( id_qicdo_han > 0 ) then
      used = send_data ( id_qicdo_han, qicdo_out, Time_next, is, js, 1)
   endif
   if ( id_qacdo_han > 0 ) then  ! ZNT 05/24/2023
      used = send_data ( id_qacdo_han, qacdo_out, Time_next, is, js, 1)
   endif
   if ( id_qncdo_han > 0 ) then  ! ZNT 05/12/2023
      used = send_data ( id_qncdo_han, qncdo_out, Time_next, is, js, 1)
   endif
   if ( id_qnicdo_han > 0 ) then ! ZNT 05/12/2023
      used = send_data ( id_qnicdo_han, qnicdo_out, Time_next, is, js, 1)
   endif
   if ( id_buod_han > 0 ) then
      used = send_data ( id_buod_han, buod_out, Time_next, is, js, 1)
   endif
   ! End ZNT 06/19/2020

   if ( id_taux_han > 0 ) then
      used = send_data ( id_taux_han, dusfc_out, Time_next, is, js)
   endif

   if ( id_tauy_han > 0 ) then
      used = send_data ( id_tauy_han, dvsfc_out, Time_next, is, js)
   endif

   if ( id_shf_han > 0 ) then
      used = send_data ( id_shf_han, dtsfc_out, Time_next, is, js)
   endif

   if ( id_lhf_han > 0 ) then
      used = send_data ( id_lhf_han, dqsfc_out, Time_next, is, js)
   endif

   if ( id_fm_han > 0 ) then
      used = send_data ( id_fm_han, fm_out, Time_next, is, js)
   endif

   if ( id_ft_han > 0 ) then
      used = send_data ( id_ft_han, ft_out, Time_next, is, js)
   endif

   if ( id_risfc_han > 0 ) then
      used = send_data ( id_risfc_han, risfc_out, Time_next, is, js)
   endif

   if ( id_zpbl_han > 0 ) then
      used = send_data ( id_zpbl_han, hpbl_out, Time_next, is, js)
   endif

! ZNT 05/08/2020, 10/11/2020: MF diagnostics
   if ( id_tdt_mf_han > 0 ) then
      used = send_data ( id_tdt_mf_han, tdt_mf_out, Time_next, is, js, 1)
   endif

   if ( id_qdt_mf_han > 0 ) then
      used = send_data ( id_qdt_mf_han, qvdt_mf_out, Time_next, is, js, 1)
   endif

   if ( id_qldt_mf_han > 0 ) then
      used = send_data ( id_qldt_mf_han, qldt_mf_out, Time_next, is, js, 1)
   endif

   if ( id_qidt_mf_han > 0 ) then
      used = send_data ( id_qidt_mf_han, qidt_mf_out, Time_next, is, js, 1)
   endif

   if ( id_qndt_mf_han > 0 ) then  ! ZNT 05/12/2023
      used = send_data ( id_qndt_mf_han, qndt_mf_out, Time_next, is, js, 1)
   endif

   if ( id_qnidt_mf_han > 0 ) then ! ZNT 05/12/2023
      used = send_data ( id_qnidt_mf_han, qnidt_mf_out, Time_next, is, js, 1)
   endif

   if ( id_qadt_mf_han > 0 ) then
      used = send_data ( id_qadt_mf_han, qadt_mf_out, Time_next, is, js, 1)
   endif

   if ( id_udt_mf_han > 0 ) then
      used = send_data ( id_udt_mf_han, udt_mf_out, Time_next, is, js, 1)
   endif

   if ( id_vdt_mf_han > 0 ) then
      used = send_data ( id_vdt_mf_han, vdt_mf_out, Time_next, is, js, 1)
   endif

   if ( id_tdt_mfd_han > 0 ) then
      used = send_data ( id_tdt_mfd_han, tdt_mfd_out, Time_next, is, js, 1)
   endif

   if ( id_qdt_mfd_han > 0 ) then
      used = send_data ( id_qdt_mfd_han, qvdt_mfd_out, Time_next, is, js, 1)
   endif

   if ( id_qldt_mfd_han > 0 ) then
      used = send_data ( id_qldt_mfd_han, qldt_mfd_out, Time_next, is, js, 1)
   endif

   if ( id_qidt_mfd_han > 0 ) then
      used = send_data ( id_qidt_mfd_han, qidt_mfd_out, Time_next, is, js, 1)
   endif

   if ( id_qndt_mfd_han > 0 ) then  ! ZNT 05/12/2023
      used = send_data ( id_qndt_mfd_han, qndt_mfd_out, Time_next, is, js, 1)
   endif

   if ( id_qnidt_mfd_han > 0 ) then ! ZNT 05/12/2023
      used = send_data ( id_qnidt_mfd_han, qnidt_mfd_out, Time_next, is, js, 1)
   endif

   if ( id_qadt_mfd_han > 0 ) then
      used = send_data ( id_qadt_mfd_han, qadt_mfd_out, Time_next, is, js, 1)
   endif

   if ( id_udt_mfd_han > 0 ) then
      used = send_data ( id_udt_mfd_han, udt_mfd_out, Time_next, is, js, 1)
   endif

   if ( id_vdt_mfd_han > 0 ) then
      used = send_data ( id_vdt_mfd_han, vdt_mfd_out, Time_next, is, js, 1)
   endif

   if ( id_taux_mf_han > 0 ) then
      used = send_data ( id_taux_mf_han, dusfc_mf_out, Time_next, is, js)
   endif

   if ( id_tauy_mf_han > 0 ) then
      used = send_data ( id_tauy_mf_han, dvsfc_mf_out, Time_next, is, js)
   endif

   if ( id_shf_mf_han > 0 ) then
      used = send_data ( id_shf_mf_han, dtsfc_mf_out, Time_next, is, js)
   endif

   if ( id_lhf_mf_han > 0 ) then
      used = send_data ( id_lhf_mf_han, dqsfc_mf_out, Time_next, is, js)
   endif

   if ( id_taux_mfd_han > 0 ) then
      used = send_data ( id_taux_mfd_han, dusfc_mfd_out, Time_next, is, js)
   endif

   if ( id_tauy_mfd_han > 0 ) then
      used = send_data ( id_tauy_mfd_han, dvsfc_mfd_out, Time_next, is, js)
   endif

   if ( id_shf_mfd_han > 0 ) then
      used = send_data ( id_shf_mfd_han, dtsfc_mfd_out, Time_next, is, js)
   endif

   if ( id_lhf_mfd_han > 0 ) then
      used = send_data ( id_lhf_mfd_han, dqsfc_mfd_out, Time_next, is, js)
   endif

! End ZNT 05/08/2020

! ZNT 05/24/2023: subs-detr diagnostics
   if ( id_tdt_sub_han > 0 ) then
      used = send_data ( id_tdt_sub_han, tdt_sub_out, Time_next, is, js, 1)
   endif

   if ( id_qdt_sub_han > 0 ) then
      used = send_data ( id_qdt_sub_han, qvdt_sub_out, Time_next, is, js, 1)
   endif

   if ( id_qldt_sub_han > 0 ) then
      used = send_data ( id_qldt_sub_han, qldt_sub_out, Time_next, is, js, 1)
   endif

   if ( id_qidt_sub_han > 0 ) then
      used = send_data ( id_qidt_sub_han, qidt_sub_out, Time_next, is, js, 1)
   endif

   if ( id_qndt_sub_han > 0 ) then
      used = send_data ( id_qndt_sub_han, qndt_sub_out, Time_next, is, js, 1)
   endif

   if ( id_qnidt_sub_han > 0 ) then
      used = send_data ( id_qnidt_sub_han, qnidt_sub_out, Time_next, is, js, 1)
   endif

   if ( id_qadt_sub_han > 0 ) then
      used = send_data ( id_qadt_sub_han, qadt_sub_out, Time_next, is, js, 1)
   endif

   if ( id_tdt_subd_han > 0 ) then
      used = send_data ( id_tdt_subd_han, tdt_subd_out, Time_next, is, js, 1)
   endif

   if ( id_qdt_subd_han > 0 ) then
      used = send_data ( id_qdt_subd_han, qvdt_subd_out, Time_next, is, js, 1)
   endif

   if ( id_qldt_subd_han > 0 ) then
      used = send_data ( id_qldt_subd_han, qldt_subd_out, Time_next, is, js, 1)
   endif

   if ( id_qidt_subd_han > 0 ) then
      used = send_data ( id_qidt_subd_han, qidt_subd_out, Time_next, is, js, 1)
   endif

   if ( id_qndt_subd_han > 0 ) then
      used = send_data ( id_qndt_subd_han, qndt_subd_out, Time_next, is, js, 1)
   endif

   if ( id_qnidt_subd_han > 0 ) then
      used = send_data ( id_qnidt_subd_han, qnidt_subd_out, Time_next, is, js, 1)
   endif

   if ( id_qadt_subd_han > 0 ) then
      used = send_data ( id_qadt_subd_han, qadt_subd_out, Time_next, is, js, 1)
   endif

   if ( id_tdt_dtr_han > 0 ) then
      used = send_data ( id_tdt_dtr_han, tdt_dtr_out, Time_next, is, js, 1)
   endif

   if ( id_qdt_dtr_han > 0 ) then
      used = send_data ( id_qdt_dtr_han, qvdt_dtr_out, Time_next, is, js, 1)
   endif

   if ( id_qldt_dtr_han > 0 ) then
      used = send_data ( id_qldt_dtr_han, qldt_dtr_out, Time_next, is, js, 1)
   endif

   if ( id_qidt_dtr_han > 0 ) then
      used = send_data ( id_qidt_dtr_han, qidt_dtr_out, Time_next, is, js, 1)
   endif

   if ( id_qndt_dtr_han > 0 ) then
      used = send_data ( id_qndt_dtr_han, qndt_dtr_out, Time_next, is, js, 1)
   endif

   if ( id_qnidt_dtr_han > 0 ) then
      used = send_data ( id_qnidt_dtr_han, qnidt_dtr_out, Time_next, is, js, 1)
   endif

   if ( id_qadt_dtr_han > 0 ) then
      used = send_data ( id_qadt_dtr_han, qadt_dtr_out, Time_next, is, js, 1)
   endif

   if ( id_tdt_dtrd_han > 0 ) then
      used = send_data ( id_tdt_dtrd_han, tdt_dtrd_out, Time_next, is, js, 1)
   endif

   if ( id_qdt_dtrd_han > 0 ) then
      used = send_data ( id_qdt_dtrd_han, qvdt_dtrd_out, Time_next, is, js, 1)
   endif

   if ( id_qldt_dtrd_han > 0 ) then
      used = send_data ( id_qldt_dtrd_han, qldt_dtrd_out, Time_next, is, js, 1)
   endif

   if ( id_qidt_dtrd_han > 0 ) then
      used = send_data ( id_qidt_dtrd_han, qidt_dtrd_out, Time_next, is, js, 1)
   endif

   if ( id_qndt_dtrd_han > 0 ) then
      used = send_data ( id_qndt_dtrd_han, qndt_dtrd_out, Time_next, is, js, 1)
   endif

   if ( id_qnidt_dtrd_han > 0 ) then
      used = send_data ( id_qnidt_dtrd_han, qnidt_dtrd_out, Time_next, is, js, 1)
   endif

   if ( id_qadt_dtrd_han > 0 ) then
      used = send_data ( id_qadt_dtrd_han, qadt_dtrd_out, Time_next, is, js, 1)
   endif


! ZNT 10/03/2020
   if ( id_de_buop_han > 0 ) then
      used = send_data ( id_de_buop_han, de_buop_out, Time_next, is, js, 1)
   endif
   if ( id_de_shrp_han > 0 ) then
      used = send_data ( id_de_shrp_han, de_shrp_out, Time_next, is, js, 1)
   endif
   if ( id_de_diss_han > 0 ) then
      used = send_data ( id_de_diss_han, de_diss_out, Time_next, is, js, 1)
   endif
! End ZNT 10/03/2020


! ZNT 05/17/2023: Diagnostic Outputs for passive T,q,ql,qi; 07/05/2023: passive qa,qn,qni
! T
   if ( id_tpcko_han > 0 ) then
      used = send_data ( id_tpcko_han, tqpcko_out(:,:,:,1), Time_next,is,js,1)
   endif
   if ( id_tpcdo_han > 0 ) then
      used = send_data ( id_tpcdo_han, tqpcdo_out(:,:,:,1), Time_next,is,js,1)
   endif
   if ( id_tpdt_han > 0 ) then
      used = send_data ( id_tpdt_han, tqpdt_out(:,:,:,1), Time_next,is,js,1)
   endif
   if ( id_tpdt_mf_han > 0 ) then ! ZNT 05/12/2023
      used = send_data ( id_tpdt_mf_han, tqpdt_mf_out(:,:,:,1), &
                         Time_next, is, js, 1)
   endif
   if ( id_tpdt_mfd_han > 0 ) then ! ZNT 05/12/2023
      used = send_data ( id_tpdt_mfd_han, tqpdt_mfd_out(:,:,:,1), &
                         Time_next, is, js, 1)
   endif
   if ( id_tpdt_sub_han > 0 ) then ! ZNT 05/23/2023
      used = send_data ( id_tpdt_sub_han, tqpdt_sub_out(:,:,:,1), &
                         Time_next, is, js, 1)
   endif
   if ( id_tpdt_subd_han > 0 ) then ! ZNT 05/23/2023
      used = send_data ( id_tpdt_subd_han, tqpdt_subd_out(:,:,:,1), &
                         Time_next, is, js, 1)
   endif
   if ( id_tpdt_dtr_han > 0 ) then ! ZNT 05/23/2023
      used = send_data ( id_tpdt_dtr_han, tqpdt_dtr_out(:,:,:,1), &
                         Time_next, is, js, 1)
   endif
   if ( id_tpdt_dtrd_han > 0 ) then ! ZNT 05/23/2023
      used = send_data ( id_tpdt_dtrd_han, tqpdt_dtrd_out(:,:,:,1), &
                         Time_next, is, js, 1)
   endif
! q (=qv)
   if ( id_qvpcko_han > 0 ) then
      used = send_data ( id_qvpcko_han, tqpcko_out(:,:,:,2), Time_next,is,js,1)
   endif
   if ( id_qvpcdo_han > 0 ) then
      used = send_data ( id_qvpcdo_han, tqpcdo_out(:,:,:,2), Time_next,is,js,1)
   endif
   if ( id_qpdt_han > 0 ) then
      used = send_data ( id_qpdt_han, tqpdt_out(:,:,:,2), Time_next,is,js,1)
   endif
   if ( id_qpdt_mf_han > 0 ) then ! ZNT 05/12/2023
      used = send_data ( id_qpdt_mf_han, tqpdt_mf_out(:,:,:,2), &
                         Time_next, is, js, 1)
   endif
   if ( id_qpdt_mfd_han > 0 ) then ! ZNT 05/12/2023
      used = send_data ( id_qpdt_mfd_han, tqpdt_mfd_out(:,:,:,2), &
                         Time_next, is, js, 1)
   endif
   if ( id_qpdt_sub_han > 0 ) then ! ZNT 05/23/2023
      used = send_data ( id_qpdt_sub_han, tqpdt_sub_out(:,:,:,2), &
                         Time_next, is, js, 1)
   endif
   if ( id_qpdt_subd_han > 0 ) then ! ZNT 05/23/2023
      used = send_data ( id_qpdt_subd_han, tqpdt_subd_out(:,:,:,2), &
                         Time_next, is, js, 1)
   endif
   if ( id_qpdt_dtr_han > 0 ) then ! ZNT 05/23/2023
      used = send_data ( id_qpdt_dtr_han, tqpdt_dtr_out(:,:,:,2), &
                         Time_next, is, js, 1)
   endif
   if ( id_qpdt_dtrd_han > 0 ) then ! ZNT 05/23/2023
      used = send_data ( id_qpdt_dtrd_han, tqpdt_dtrd_out(:,:,:,2), &
                         Time_next, is, js, 1)
   endif
! ql
   if ( id_qlpcko_han > 0 ) then
      used = send_data ( id_qlpcko_han, tqpcko_out(:,:,:,3), Time_next,is,js,1)
   endif
   if ( id_qlpcdo_han > 0 ) then
      used = send_data ( id_qlpcdo_han, tqpcdo_out(:,:,:,3), Time_next,is,js,1)
   endif
   if ( id_qlpdt_han > 0 ) then
      used = send_data ( id_qlpdt_han, tqpdt_out(:,:,:,3), Time_next,is,js,1)
   endif
   if ( id_qlpdt_mf_han > 0 ) then ! ZNT 05/12/2023
      used = send_data ( id_qlpdt_mf_han, tqpdt_mf_out(:,:,:,3), &
                         Time_next, is, js, 1)
   endif
   if ( id_qlpdt_mfd_han > 0 ) then ! ZNT 05/12/2023
      used = send_data ( id_qlpdt_mfd_han, tqpdt_mfd_out(:,:,:,3), &
                         Time_next, is, js, 1)
   endif
   if ( id_qlpdt_sub_han > 0 ) then ! ZNT 05/23/2023
      used = send_data ( id_qlpdt_sub_han, tqpdt_sub_out(:,:,:,3), &
                         Time_next, is, js, 1)
   endif
   if ( id_qlpdt_subd_han > 0 ) then ! ZNT 05/23/2023
      used = send_data ( id_qlpdt_subd_han, tqpdt_subd_out(:,:,:,3), &
                         Time_next, is, js, 1)
   endif
   if ( id_qlpdt_dtr_han > 0 ) then ! ZNT 05/23/2023
      used = send_data ( id_qlpdt_dtr_han, tqpdt_dtr_out(:,:,:,3), &
                         Time_next, is, js, 1)
   endif
   if ( id_qlpdt_dtrd_han > 0 ) then ! ZNT 05/23/2023
      used = send_data ( id_qlpdt_dtrd_han, tqpdt_dtrd_out(:,:,:,3), &
                         Time_next, is, js, 1)
   endif
! qi
   if ( id_qipcko_han > 0 ) then
      used = send_data ( id_qipcko_han, tqpcko_out(:,:,:,4), Time_next,is,js,1)
   endif
   if ( id_qipcdo_han > 0 ) then
      used = send_data ( id_qipcdo_han, tqpcdo_out(:,:,:,4), Time_next,is,js,1)
   endif
   if ( id_qipdt_han > 0 ) then
      used = send_data ( id_qipdt_han, tqpdt_out(:,:,:,4), Time_next,is,js,1)
   endif
   if ( id_qipdt_mf_han > 0 ) then ! ZNT 05/12/2023
      used = send_data ( id_qipdt_mf_han, tqpdt_mf_out(:,:,:,4), &
                         Time_next, is, js, 1)
   endif
   if ( id_qipdt_mfd_han > 0 ) then ! ZNT 05/12/2023
      used = send_data ( id_qipdt_mfd_han, tqpdt_mfd_out(:,:,:,4), &
                         Time_next, is, js, 1)
   endif
   if ( id_qipdt_sub_han > 0 ) then ! ZNT 05/23/2023
      used = send_data ( id_qipdt_sub_han, tqpdt_sub_out(:,:,:,4), &
                         Time_next, is, js, 1)
   endif
   if ( id_qipdt_subd_han > 0 ) then ! ZNT 05/23/2023
      used = send_data ( id_qipdt_subd_han, tqpdt_subd_out(:,:,:,4), &
                         Time_next, is, js, 1)
   endif
   if ( id_qipdt_dtr_han > 0 ) then ! ZNT 05/23/2023
      used = send_data ( id_qipdt_dtr_han, tqpdt_dtr_out(:,:,:,4), &
                         Time_next, is, js, 1)
   endif
   if ( id_qipdt_dtrd_han > 0 ) then ! ZNT 05/23/2023
      used = send_data ( id_qipdt_dtrd_han, tqpdt_dtrd_out(:,:,:,4), &
                         Time_next, is, js, 1)
   endif
! qa
   if ( id_qapcko_han > 0 ) then
      used = send_data ( id_qapcko_han, tqpcko_out(:,:,:,5), Time_next,is,js,1)
   endif
   if ( id_qapcdo_han > 0 ) then
      used = send_data ( id_qapcdo_han, tqpcdo_out(:,:,:,5), Time_next,is,js,1)
   endif
   if ( id_qapdt_han > 0 ) then
      used = send_data ( id_qapdt_han, tqpdt_out(:,:,:,5), Time_next,is,js,1)
   endif
   if ( id_qapdt_mf_han > 0 ) then ! ZNT 05/12/2023
      used = send_data ( id_qapdt_mf_han, tqpdt_mf_out(:,:,:,5), &
                         Time_next, is, js, 1)
   endif
   if ( id_qapdt_mfd_han > 0 ) then ! ZNT 05/12/2023
      used = send_data ( id_qapdt_mfd_han, tqpdt_mfd_out(:,:,:,5), &
                         Time_next, is, js, 1)
   endif
   if ( id_qapdt_sub_han > 0 ) then ! ZNT 05/23/2023
      used = send_data ( id_qapdt_sub_han, tqpdt_sub_out(:,:,:,5), &
                         Time_next, is, js, 1)
   endif
   if ( id_qapdt_subd_han > 0 ) then ! ZNT 05/23/2023
      used = send_data ( id_qapdt_subd_han, tqpdt_subd_out(:,:,:,5), &
                         Time_next, is, js, 1)
   endif
   if ( id_qapdt_dtr_han > 0 ) then ! ZNT 05/23/2023
      used = send_data ( id_qapdt_dtr_han, tqpdt_dtr_out(:,:,:,5), &
                         Time_next, is, js, 1)
   endif
   if ( id_qapdt_dtrd_han > 0 ) then ! ZNT 05/23/2023
      used = send_data ( id_qapdt_dtrd_han, tqpdt_dtrd_out(:,:,:,5), &
                         Time_next, is, js, 1)
   endif
! qn
   if ( id_qnpcko_han > 0 ) then
      used = send_data ( id_qnpcko_han, tqpcko_out(:,:,:,6), Time_next,is,js,1)
   endif
   if ( id_qnpcdo_han > 0 ) then
      used = send_data ( id_qnpcdo_han, tqpcdo_out(:,:,:,6), Time_next,is,js,1)
   endif
   if ( id_qnpdt_han > 0 ) then
      used = send_data ( id_qnpdt_han, tqpdt_out(:,:,:,6), Time_next,is,js,1)
   endif
   if ( id_qnpdt_mf_han > 0 ) then ! ZNT 05/12/2023
      used = send_data ( id_qnpdt_mf_han, tqpdt_mf_out(:,:,:,6), &
                         Time_next, is, js, 1)
   endif
   if ( id_qnpdt_mfd_han > 0 ) then ! ZNT 05/12/2023
      used = send_data ( id_qnpdt_mfd_han, tqpdt_mfd_out(:,:,:,6), &
                         Time_next, is, js, 1)
   endif
   if ( id_qnpdt_sub_han > 0 ) then ! ZNT 05/23/2023
      used = send_data ( id_qnpdt_sub_han, tqpdt_sub_out(:,:,:,6), &
                         Time_next, is, js, 1)
   endif
   if ( id_qnpdt_subd_han > 0 ) then ! ZNT 05/23/2023
      used = send_data ( id_qnpdt_subd_han, tqpdt_subd_out(:,:,:,6), &
                         Time_next, is, js, 1)
   endif
   if ( id_qnpdt_dtr_han > 0 ) then ! ZNT 05/23/2023
      used = send_data ( id_qnpdt_dtr_han, tqpdt_dtr_out(:,:,:,6), &
                         Time_next, is, js, 1)
   endif
   if ( id_qnpdt_dtrd_han > 0 ) then ! ZNT 05/23/2023
      used = send_data ( id_qnpdt_dtrd_han, tqpdt_dtrd_out(:,:,:,6), &
                         Time_next, is, js, 1)
   endif
! qni
   if ( id_qnipcko_han > 0 ) then
      used = send_data ( id_qnipcko_han, tqpcko_out(:,:,:,7), Time_next,is,js,1)
   endif
   if ( id_qnipcdo_han > 0 ) then
      used = send_data ( id_qnipcdo_han, tqpcdo_out(:,:,:,7), Time_next,is,js,1)
   endif
   if ( id_qnipdt_han > 0 ) then
      used = send_data ( id_qnipdt_han, tqpdt_out(:,:,:,7), Time_next,is,js,1)
   endif
   if ( id_qnipdt_mf_han > 0 ) then ! ZNT 05/12/2023
      used = send_data ( id_qnipdt_mf_han, tqpdt_mf_out(:,:,:,7), &
                         Time_next, is, js, 1)
   endif
   if ( id_qnipdt_mfd_han > 0 ) then ! ZNT 05/12/2023
      used = send_data ( id_qnipdt_mfd_han, tqpdt_mfd_out(:,:,:,7), &
                         Time_next, is, js, 1)
   endif
   if ( id_qnipdt_sub_han > 0 ) then ! ZNT 05/23/2023
      used = send_data ( id_qnipdt_sub_han, tqpdt_sub_out(:,:,:,7), &
                         Time_next, is, js, 1)
   endif
   if ( id_qnipdt_subd_han > 0 ) then ! ZNT 05/23/2023
      used = send_data ( id_qnipdt_subd_han, tqpdt_subd_out(:,:,:,7), &
                         Time_next, is, js, 1)
   endif
   if ( id_qnipdt_dtr_han > 0 ) then ! ZNT 05/23/2023
      used = send_data ( id_qnipdt_dtr_han, tqpdt_dtr_out(:,:,:,7), &
                         Time_next, is, js, 1)
   endif
   if ( id_qnipdt_dtrd_han > 0 ) then ! ZNT 05/23/2023
      used = send_data ( id_qnipdt_dtrd_han, tqpdt_dtrd_out(:,:,:,7), &
                         Time_next, is, js, 1)
   endif
! End ZNT 05/17/2023; 07/05/2023


! ZNT 12/08/2020: Lock scheme PBL depth algorithm
! Need to allocate: zsml(nlon,nlat), zfull_ag(nlon,nlat,nlev), slv(nlon,nlat,nlev), 
!                   zsurf(nlon,nlat), ipbl, ibot, tmpjump, critjump
! Need to pass kbot as input to hanedmf
! Need to use: mo_diff
   
   hleff   = (min(1.,max(0.,0.05*(t       -tfreeze+20.)))*hlv + &
              min(1.,max(0.,0.05*(tfreeze -t          )))*hls)
   slv     = cp_air*t + grav*zfull_ag - hleff*(ql + qi)
   slv     = slv*(1+d608*(qv+ql+qi))

   ibot = nlev
   zsml = 0.0   ! note that this must be zero as this is 
                ! indicates stable surface layer and this
                ! value is output for use in gravity
                ! wave drag scheme
   do i=1,size(t,1)
     do j=1,size(t,2)
       ipbl    = -1
       if (present(kbot)) ibot = kbot(i,j)
       if (b_star(i,j) .gt. 0.) then
         call pbl_depth(slv(i,j,1:ibot)/cp_air, &
              zfull_ag(i,j,1:ibot),u_star(i,j),b_star(i,j),    &
              ipbl,zsml(i,j))
         ! if ((i==1).and.(j==1)) write (*,*) zsml(i,j)
         !------------------------------------------------------
         ! Following Lock et al. 2000, limit height of surface
         ! well mixed layer if interior stable interface is
         ! found.  An interior stable interface is diagnosed if
         ! the slope between 2 full levels is greater than
         ! the namelist parameter critjump
         if (ipbl .lt. ibot) then 
              do k = ibot, ipbl+1, -1
                   tmpjump =(slv(i,j,k-1)-slv(i,j,k))/cp_air 
                   if (tmpjump .gt. critjump) then
                         ipbl = k
                         zsml(i,j) = zhalf_ag(i,j,ipbl)
                         exit
                   end if
             enddo
         end if ! (ipbl .lt. ibot)
       end if ! (b_star(i,j) .gt. 0.)
     enddo
   enddo
   ! write (*,*) zsml(1,1)

   if (use_lock_zpbl) then
       if (use_z_above_sfc) then 
       ! ZNT 03/08/2023: Use height above surface. This should be the correct implementation
          hpbl_out = zsml(:,:)
       else 
       ! ZNT 03/08/2023: Use height above sea level. Used in 2021-22 but incorrect.
          hpbl_out = zsml(:,:) + zsurf(:,:)
       endif
   endif

   if ( id_zsml_han > 0 ) then
     used = send_data ( id_zsml_han, zsml, Time_next, is, js)
   endif
! End ZNT 12/08/2020

! ZNT 01/25/2024: zero out the MF tendencies or ED coefficients if edmf_level is 0 or 1.
   if (edmf_level == 0) then                       ! zero out MF tends (for diag run)
      udt_mfud   = 0.0
      vdt_mfud   = 0.0
      tdt_mfud   = 0.0
      qdt_mfud   = 0.0
      qldt_mfud  = 0.0
      qidt_mfud  = 0.0
      qadt_mfud  = 0.0
      qndt_mfud  = 0.0
      qnidt_mfud = 0.0
   endif

   if (edmf_level == 0 .or. edmf_level == 1) then  ! zero out ED coeff (for diag or mfonly run)
      diff_m = 0.0
      diff_t = 0.0
      if (use_z_above_sfc) then
      ! ZNT 03/08/2023: Use height above surface. This should be the correct implementation
         hpbl_out = 0.0
      else
      ! ZNT 03/08/2023: Use height above sea level. Used in 2021-22 but incorrect.
         hpbl_out = 0.0 + zsurf(:,:)
      endif 
   endif

   end subroutine hanedmf
!#######################################################################

subroutine pbl_depth(t, z, u_star, b_star, ipbl, h)  ! From Lock scheme

!
!  -----
!  INPUT
!  -----
!
!  t (= slv/cp)  liquid water virtual static energy divided by cp (K)
!  u_star        friction velocity (m/s)
!  b_star        buoyancy scale (m/s**2)
!       
!  ------
!  OUTPUT
!  ------
!
!  ipbl          half level containing pbl height
!  h             pbl height (m)

real,    intent(in) ,  dimension(:) :: t, z
real,    intent(in)                 :: u_star, b_star
integer, intent(out)                :: ipbl
real,    intent(out)                :: h

real     :: svp,h1,h2,t1,t2
real     :: ws,k_t_ref
integer  :: k,nlev
real     :: small  = 1.e-4

!initialize zsml
h = 0.

!compute # of levels
nlev = size(t,1)

!calculate surface parcel properties
svp  = t(nlev)
h1   = z(nlev)
call mo_diff(h1, u_star, b_star, ws, k_t_ref)
ws = max(small,ws/vonkarm/h1)
svp  = svp*(1.+(parcel_buoy*u_star*b_star/grav/ws) )

!search for level where this is exceeded              
h    = h1
t1   = t(nlev)
do k = nlev-1 , 2, -1
     h2 = z(k)
     t2 = t(k)
     if (t2.gt.svp) then
          h = h2 + (h1 - h2)*(t2 - svp)/(t2 - t1 )
          ipbl = k+1
          return
     end if
     h1 = h2
     t1 = t2
enddo

!one shouldn't end up here but nonetheless for safety this is put here
h = h2
ipbl = k+1

return

end subroutine pbl_depth


!#######################################################################

   subroutine hanedmf_init (axes, Time)

!-----------------------------------------------------------------------
!
!        initialization for large scale condensation
!
!-----------------------------------------------------------------------
  integer,         intent(in) :: axes(4)
  type(time_type), intent(in) :: Time

  integer  unit,io,ierr, logunit
  character(len=100):: errmsg
  integer :: errflg
  integer, dimension(3) :: full = (/1,2,3/), half = (/1,2,4/)


      if (module_is_initialized)  &
          call error_mesg  &
                   ('hanedmf_init in hanedmf_mod',  &
                    'attempting to call initialization twice', FATAL)

!-----------------------------------------------------------------------
!--------------- read namelist ------------------

      read (input_nml_file, nml=hanedmf_nml, iostat=io)
      ierr = check_nml_error(io,"hanedmf_nml")

!---------- output namelist --------------------------------------------

      logunit = stdlog()
      if ( mpp_pe() == mpp_root_pe() ) then
           call write_version_number(version, tagname)
           write (logunit,nml=hanedmf_nml)
      endif


      if ( use_online_aerosol ) call aer_ccn_act_init ! ZNT 06/22/2023

!---------- do native initialization -----------------------------------
call satmedmfvdifq_init (1,1,errmsg,errflg)

!---------- initialize diag fields -------------------------------------

   id_zfull_han = register_diag_field ( mod_name, &   ! ZNT 05/20/2023
     'zfull_han', axes(full), Time, &
     'zfull seen by Han EDMF', 'm' )
   id_zhalf_han = register_diag_field ( mod_name, &   ! ZNT 05/20/2023
     'zhalf_han', axes(half), Time, &
     'zhalf seen by Han EDMF', 'm' )
   id_pfull_han = register_diag_field ( mod_name, &   ! ZNT 05/20/2023
     'pfull_han', axes(full), Time, &
     'pfull seen by Han EDMF', 'Pa' )
   id_phalf_han = register_diag_field ( mod_name, &   ! ZNT 05/20/2023
     'phalf_han', axes(half), Time, &
     'phalf seen by Han EDMF', 'Pa' )

   id_tdt_han = register_diag_field ( mod_name, &
     'tdt_han', axes(full), Time, &
     'Temperature tendency by Han EDMF', 'K/s' )

   id_qdt_han = register_diag_field ( mod_name, &
     'qdt_han', axes(full), Time, &
     'Specific humidity tendency by Han EDMF', 'kg/kg/s' )

   id_udt_han = register_diag_field ( mod_name, &
     'udt_han', axes(full), Time, &
     'U wind tendency by Han EDMF', 'm/s^2' )

   id_vdt_han = register_diag_field ( mod_name, &
     'vdt_han', axes(full), Time, &
     'V wind tendency by Han EDMF', 'm/s^2' )

   id_qldt_han = register_diag_field ( mod_name, &
     'qldt_han', axes(full), Time, &
     'Cloud liquid tendency by Han EDMF', 'kg/kg/s' )

   id_qidt_han = register_diag_field ( mod_name, &
     'qidt_han', axes(full), Time, &
     'Cloud ice tendency by Han EDMF', 'kg/kg/s' )

   id_qndt_han = register_diag_field ( mod_name, &   ! ZNT 05/12/2023
     'qndt_han', axes(full), Time, &
     'Cloud droplet number tendency by Han EDMF', '1/kg/s' )

   id_qnidt_han = register_diag_field ( mod_name, &  ! ZNT 05/12/2023
     'qnidt_han', axes(full), Time, &
     'Cloud ice number tendency by Han EDMF', '1/kg/s' )

   id_qadt_han = register_diag_field ( mod_name, &
     'qadt_han', axes(full), Time, &
     'Cloud amount tendency by Han EDMF', 'kg/kg/s' )

   id_tke_han = register_diag_field ( mod_name, &
     'tke_han', axes(full), Time, &
     'TKE (diagnostic) by Han EDMF', 'm^2/s^2' )

   id_km_han = register_diag_field ( mod_name, &
     'Km_han', axes(half), Time, &
     'Momentum eddy diffusivity by Han EDMF', 'm/s^2' )

   id_kt_han = register_diag_field ( mod_name, &
     'Kt_han', axes(half), Time, &
     'T and tracer eddy diffusivity by Han EDMF', 'm/s^2' )

   id_ke_han = register_diag_field ( mod_name, &
     'Ke_han', axes(half), Time, &
     'TKE eddy diffusivity by Han EDMF', 'm/s^2' )

   id_mu_han = register_diag_field ( mod_name, &
     'Mu_han', axes(full), Time, &
     'Updraft mass flux by Han EDMF', 'm/s' )

   id_md_han = register_diag_field ( mod_name, &
     'Md_han', axes(full), Time, &
     'Stratocumulus downdraft mass flux by Han EDMF', 'm/s' )

! ZNT 06/19/2020
   id_tcko_han = register_diag_field ( mod_name, &
     'tcko_han', axes(full), Time, &
     'Updraft temperature', 'K' )
   id_qvcko_han = register_diag_field ( mod_name, &
     'qvcko_han', axes(full), Time, &
     'Updraft qv', 'kg/kg' )
   id_qlcko_han = register_diag_field ( mod_name, &
     'qlcko_han', axes(full), Time, &
     'Updraft ql', 'kg/kg' )
   id_qicko_han = register_diag_field ( mod_name, &
     'qicko_han', axes(full), Time, &
     'Updraft qi', 'kg/kg' )
   id_qacko_han = register_diag_field ( mod_name, &  ! ZNT 05/24/2023
     'qacko_han', axes(full), Time, &
     'Updraft qa', '1' )
   id_qncko_han = register_diag_field ( mod_name, &  ! ZNT 05/12/2023
     'qncko_han', axes(full), Time, &
     'Updraft qn', '1/kg' )
   id_qnicko_han = register_diag_field ( mod_name, & ! ZNT 05/12/2023
     'qnicko_han', axes(full), Time, &
     'Updraft qni', '1/kg' )
   id_buou_han = register_diag_field ( mod_name, &
     'buou_han', axes(full), Time, &
     'Updraft buoyancy', 'm/s2' )
   id_tcdo_han = register_diag_field ( mod_name, &
     'tcdo_han', axes(full), Time, &
     'Downdraft temperature', 'K' )
   id_qvcdo_han = register_diag_field ( mod_name, &
     'qvcdo_han', axes(full), Time, &
     'Downdraft qv', 'kg/kg' )
   id_qlcdo_han = register_diag_field ( mod_name, &
     'qlcdo_han', axes(full), Time, &
     'Downdraft ql', 'kg/kg' )
   id_qicdo_han = register_diag_field ( mod_name, &
     'qicdo_han', axes(full), Time, &
     'Downdraft qi', 'kg/kg' )
   id_qacdo_han = register_diag_field ( mod_name, &  ! ZNT 05/24/2023
     'qacdo_han', axes(full), Time, &
     'Downdraft qa', '1' )
   id_qncdo_han = register_diag_field ( mod_name, &  ! ZNT 05/12/2023
     'qncdo_han', axes(full), Time, &
     'Downdraft qn', '1/kg' )
   id_qnicdo_han = register_diag_field ( mod_name, & ! ZNT 05/12/2023
     'qnicdo_han', axes(full), Time, &
     'Downdraft qni', '1/kg' )
   id_buod_han = register_diag_field ( mod_name, &
     'buod_han', axes(full), Time, &
     'Downdraft buoyancy', 'm/s2' )
! End of ZNT 06/19/2020

   id_shf_han = register_diag_field ( mod_name, &
     'shf_han', axes(1:2), Time, &
     'Surface sensible heat flux by Han EDMF', 'W/m2' )

   id_lhf_han = register_diag_field ( mod_name, &
     'lhf_han', axes(1:2), Time, &
     'Surface latent heat flux by Han EDMF', 'W/m2' )

   id_taux_han = register_diag_field ( mod_name, &
     'taux_han', axes(1:2), Time, &
     'Surface X-drag by Han EDMF', 'Pa' )

   id_tauy_han = register_diag_field ( mod_name, &
     'tauy_han', axes(1:2), Time, &
     'Surface Y-drag by Han EDMF', 'Pa' )

   id_risfc_han = register_diag_field ( mod_name, &
     'risfc_han', axes(1:2), Time, &
     'Surface Richardson number by Han EDMF', 'none' )

   id_fm_han = register_diag_field ( mod_name, &
     'fm_han', axes(1:2), Time, &
     'M-O function for momentum by Han EDMF', 'none' )

   id_ft_han = register_diag_field ( mod_name, &
     'ft_han', axes(1:2), Time, &
     'M-O function for heat by Han EDMF', 'none' )


   id_zpbl_han = register_diag_field ( mod_name, &
     'zpbl_han', axes(1:2), Time, &
     'Depth of Han boundary layer', 'm' )

! ZNT 05/08/2020: MF diagnostics
   id_tdt_mf_han = register_diag_field ( mod_name, &
     'tdt_mf_han', axes(full), Time, &
     'Temperature tendency by Han EDMF - Updraft', 'K/s' )

   id_qdt_mf_han = register_diag_field ( mod_name, &
     'qdt_mf_han', axes(full), Time, &
     'Specific humidity tendency by Han EDMF - Updraft', 'kg/kg/s' )

   id_qldt_mf_han = register_diag_field ( mod_name, &
     'qldt_mf_han', axes(full), Time, &
     'Cloud liquid tendency by Han EDMF - Updraft', 'kg/kg/s' )

   id_qidt_mf_han = register_diag_field ( mod_name, &
     'qidt_mf_han', axes(full), Time, &
     'Cloud ice tendency by Han EDMF - Updraft', 'kg/kg/s' )

   id_qndt_mf_han = register_diag_field ( mod_name, &  ! ZNT 05/12/2023
     'qndt_mf_han', axes(full), Time, &
     'Cloud droplet number  tendency by Han EDMF - Updraft', '1/kg/s' )

   id_qnidt_mf_han = register_diag_field ( mod_name, & ! ZNT 05/12/2023
     'qnidt_mf_han', axes(full), Time, &
     'Cloud ice number  tendency by Han EDMF - Updraft', '1/kg/s' )

   id_qadt_mf_han = register_diag_field ( mod_name, &
     'qadt_mf_han', axes(full), Time, &
     'Cloud amount tendency by Han EDMF - Updraft', 'kg/kg/s' )

   id_udt_mf_han = register_diag_field ( mod_name, &
     'udt_mf_han', axes(full), Time, &
     'U wind tendency by Han EDMF - Updraft', 'm/s^2' )

   id_vdt_mf_han = register_diag_field ( mod_name, &
     'vdt_mf_han', axes(full), Time, &
     'V wind tendency by Han EDMF - Updraft', 'm/s^2' )

   id_tdt_mfd_han = register_diag_field ( mod_name, &
     'tdt_mfd_han', axes(full), Time, &
     'Temperature tendency by Han EDMF - Downdraft', 'K/s' )

   id_qdt_mfd_han = register_diag_field ( mod_name, &
     'qdt_mfd_han', axes(full), Time, &
     'Specific humidity tendency by Han EDMF - Downdraft', 'kg/kg/s' )

   id_qldt_mfd_han = register_diag_field ( mod_name, &
     'qldt_mfd_han', axes(full), Time, &
     'Cloud liquid tendency by Han EDMF - Downdraft', 'kg/kg/s' )

   id_qidt_mfd_han = register_diag_field ( mod_name, &
     'qidt_mfd_han', axes(full), Time, &
     'Cloud ice tendency by Han EDMF - Downdraft', 'kg/kg/s' )

   id_qndt_mfd_han = register_diag_field ( mod_name, &  ! ZNT 05/12/2023
     'qndt_mfd_han', axes(full), Time, &
     'Cloud droplet number  tendency by Han EDMF - Downdraft', '1/kg/s' )

   id_qnidt_mfd_han = register_diag_field ( mod_name, & ! ZNT 05/12/2023
     'qnidt_mfd_han', axes(full), Time, &
     'Cloud ice number  tendency by Han EDMF - Downdraft', '1/kg/s' )

   id_qadt_mfd_han = register_diag_field ( mod_name, &
     'qadt_mfd_han', axes(full), Time, &
     'Cloud amount tendency by Han EDMF - Downdraft', 'kg/kg/s' )

   id_udt_mfd_han = register_diag_field ( mod_name, &
     'udt_mfd_han', axes(full), Time, &
     'U wind tendency by Han EDMF - Downdraft', 'm/s^2' )

   id_vdt_mfd_han = register_diag_field ( mod_name, &
     'vdt_mfd_han', axes(full), Time, &
     'V wind tendency by Han EDMF - Downdraft', 'm/s^2' )

   id_shf_mf_han = register_diag_field ( mod_name, &
     'shf_mf_han', axes(1:2), Time, &
     'Surface sensible heat flux by Han EDMF - Updraft', 'W/m2' )

   id_lhf_mf_han = register_diag_field ( mod_name, &
     'lhf_mf_han', axes(1:2), Time, &
     'Surface latent heat flux by Han EDMF - Updraft', 'W/m2' )

   id_taux_mf_han = register_diag_field ( mod_name, &
     'taux_mf_han', axes(1:2), Time, &
     'Surface X-drag by Han EDMF - Updraft', 'Pa' )

   id_tauy_mf_han = register_diag_field ( mod_name, &
     'tauy_mf_han', axes(1:2), Time, &
     'Surface Y-drag by Han EDMF - Updraft', 'Pa' )

   id_shf_mfd_han = register_diag_field ( mod_name, &
     'shf_mfd_han', axes(1:2), Time, &
     'Surface sensible heat flux by Han EDMF - Downdraft', 'W/m2' )

   id_lhf_mfd_han = register_diag_field ( mod_name, &
     'lhf_mfd_han', axes(1:2), Time, &
     'Surface latent heat flux by Han EDMF - Downdraft', 'W/m2' )

   id_taux_mfd_han = register_diag_field ( mod_name, &
     'taux_mfd_han', axes(1:2), Time, &
     'Surface X-drag by Han EDMF - Downdraft', 'Pa' )

   id_tauy_mfd_han = register_diag_field ( mod_name, &
     'tauy_mfd_han', axes(1:2), Time, &
     'Surface Y-drag by Han EDMF - Downdraft', 'Pa' )

! End ZNT 05/08/2020


! ZNT 05/23/2023: subs-detr diagnostics
   id_tdt_sub_han = register_diag_field ( mod_name, &
     'tdt_sub_han', axes(full), Time, &
     'Temperature tendency by Han EDMF - Updraft-subsidence', 'K/s' )

   id_qdt_sub_han = register_diag_field ( mod_name, &
     'qdt_sub_han', axes(full), Time, &
     'Specific humidity tendency by Han EDMF - Updraft-subsidence', 'kg/kg/s' )

   id_qldt_sub_han = register_diag_field ( mod_name, &
     'qldt_sub_han', axes(full), Time, &
     'Cloud liquid tendency by Han EDMF - Updraft-subsidence', 'kg/kg/s' )

   id_qidt_sub_han = register_diag_field ( mod_name, &
     'qidt_sub_han', axes(full), Time, &
     'Cloud ice tendency by Han EDMF - Updraft-subsidence', 'kg/kg/s' )

   id_qndt_sub_han = register_diag_field ( mod_name, &  ! ZNT 05/12/2023
     'qndt_sub_han', axes(full), Time, &
     'Cloud droplet number  tendency by Han EDMF - Updraft-subsidence', '1/kg/s' )

   id_qnidt_sub_han = register_diag_field ( mod_name, & ! ZNT 05/12/2023
     'qnidt_sub_han', axes(full), Time, &
     'Cloud ice number  tendency by Han EDMF - Updraft-subsidence', '1/kg/s' )

   id_qadt_sub_han = register_diag_field ( mod_name, &
     'qadt_sub_han', axes(full), Time, &
     'Cloud amount tendency by Han EDMF - Updraft-subsidence', 'kg/kg/s' )

   id_tdt_subd_han = register_diag_field ( mod_name, &
     'tdt_subd_han', axes(full), Time, &
     'Temperature tendency by Han EDMF - Downdraft-ascent', 'K/s' )

   id_qdt_subd_han = register_diag_field ( mod_name, &
     'qdt_subd_han', axes(full), Time, &
     'Specific humidity tendency by Han EDMF - Downdraft-ascent', 'kg/kg/s' )

   id_qldt_subd_han = register_diag_field ( mod_name, &
     'qldt_subd_han', axes(full), Time, &
     'Cloud liquid tendency by Han EDMF - Downdraft-ascent', 'kg/kg/s' )

   id_qidt_subd_han = register_diag_field ( mod_name, &
     'qidt_subd_han', axes(full), Time, &
     'Cloud ice tendency by Han EDMF - Downdraft-ascent', 'kg/kg/s' )

   id_qndt_subd_han = register_diag_field ( mod_name, &  ! ZNT 05/12/2023
     'qndt_subd_han', axes(full), Time, &
     'Cloud droplet number  tendency by Han EDMF - Downdraft-ascent', '1/kg/s' )

   id_qnidt_subd_han = register_diag_field ( mod_name, & ! ZNT 05/12/2023
     'qnidt_subd_han', axes(full), Time, &
     'Cloud ice number  tendency by Han EDMF - Downdraft-ascent', '1/kg/s' )

   id_qadt_subd_han = register_diag_field ( mod_name, &
     'qadt_subd_han', axes(full), Time, &
     'Cloud amount tendency by Han EDMF - Downdraft-ascent', 'kg/kg/s' )

   id_tdt_dtr_han = register_diag_field ( mod_name, &
     'tdt_dtr_han', axes(full), Time, &
     'Temperature tendency by Han EDMF - Updraft-detrainment', 'K/s' )

   id_qdt_dtr_han = register_diag_field ( mod_name, &
     'qdt_dtr_han', axes(full), Time, &
     'Specific humidity tendency by Han EDMF - Updraft-detrainment', 'kg/kg/s' )

   id_qldt_dtr_han = register_diag_field ( mod_name, &
     'qldt_dtr_han', axes(full), Time, &
     'Cloud liquid tendency by Han EDMF - Updraft-detrainment', 'kg/kg/s' )

   id_qidt_dtr_han = register_diag_field ( mod_name, &
     'qidt_dtr_han', axes(full), Time, &
     'Cloud ice tendency by Han EDMF - Updraft-detrainment', 'kg/kg/s' )

   id_qndt_dtr_han = register_diag_field ( mod_name, &  ! ZNT 05/12/2023
     'qndt_dtr_han', axes(full), Time, &
     'Cloud droplet number  tendency by Han EDMF - Updraft-detrainment', '1/kg/s' )

   id_qnidt_dtr_han = register_diag_field ( mod_name, & ! ZNT 05/12/2023
     'qnidt_dtr_han', axes(full), Time, &
     'Cloud ice number  tendency by Han EDMF - Updraft-detrainment', '1/kg/s' )

   id_qadt_dtr_han = register_diag_field ( mod_name, &
     'qadt_dtr_han', axes(full), Time, &
     'Cloud amount tendency by Han EDMF - Updraft-detrainment', 'kg/kg/s' )

   id_tdt_dtrd_han = register_diag_field ( mod_name, &
     'tdt_dtrd_han', axes(full), Time, &
     'Temperature tendency by Han EDMF - Downdraft-detrainment', 'K/s' )

   id_qdt_dtrd_han = register_diag_field ( mod_name, &
     'qdt_dtrd_han', axes(full), Time, &
     'Specific humidity tendency by Han EDMF - Downdraft-detrainment', 'kg/kg/s' )

   id_qldt_dtrd_han = register_diag_field ( mod_name, &
     'qldt_dtrd_han', axes(full), Time, &
     'Cloud liquid tendency by Han EDMF - Downdraft-detrainment', 'kg/kg/s' )

   id_qidt_dtrd_han = register_diag_field ( mod_name, &
     'qidt_dtrd_han', axes(full), Time, &
     'Cloud ice tendency by Han EDMF - Downdraft-detrainment', 'kg/kg/s' )

   id_qndt_dtrd_han = register_diag_field ( mod_name, &  ! ZNT 05/12/2023
     'qndt_dtrd_han', axes(full), Time, &
     'Cloud droplet number  tendency by Han EDMF - Downdraft-detrainment', '1/kg/s' )

   id_qnidt_dtrd_han = register_diag_field ( mod_name, & ! ZNT 05/12/2023
     'qnidt_dtrd_han', axes(full), Time, &
     'Cloud ice number  tendency by Han EDMF - Downdraft-detrainment', '1/kg/s' )

   id_qadt_dtrd_han = register_diag_field ( mod_name, &
     'qadt_dtrd_han', axes(full), Time, &
     'Cloud amount tendency by Han EDMF - Downdraft-detrainment', 'kg/kg/s' )
! End ZNT 05/23/2023

   
! ZNT 10/03/2020: TKE diagnostics
   id_de_buop_han = register_diag_field ( mod_name, &
     'de_buop_han', axes(full), Time, &
     'TKE tendency by Han EDMF - Buoyancy Production', 'm^2/s^3' )
   id_de_shrp_han  = register_diag_field ( mod_name, &
     'de_shrp_han', axes(full), Time, &
     'TKE tendency by Han EDMF - Shear Production', 'm^2/s^3' )
   id_de_diss_han = register_diag_field ( mod_name, &
     'de_diss_han', axes(full), Time, &
     'TKE tendency by Han EDMF - Viscous Dissipation', 'm^2/s^3' )
! End ZNT 10/03/2020

   id_zsml_han = register_diag_field ( mod_name, &
     'zsml_han', axes(1:2), Time, &
     'Depth of Lock boundary layer diagnosed from Han scheme', 'm' )


! ZNT 05/17/2023: Diagnostic Outputs for passive T,q,ql,qi; 07/05/2023: passive qa,qn,qni
   id_tpcko_han = register_diag_field ( mod_name, &
     'tpcko_han', axes(full), Time, &
     'Updraft temperature (passive)', 'K' )
   id_qvpcko_han = register_diag_field ( mod_name, &
     'qvpcko_han', axes(full), Time, &
     'Updraft qv (passive)', 'kg/kg' )
   id_qlpcko_han = register_diag_field ( mod_name, &
     'qlpcko_han', axes(full), Time, &
     'Updraft ql (passive)', 'kg/kg' )
   id_qipcko_han = register_diag_field ( mod_name, &
     'qipcko_han', axes(full), Time, &
     'Updraft qi (passive)', 'kg/kg' )
   id_qapcko_han = register_diag_field ( mod_name, &
     'qapcko_han', axes(full), Time, &
     'Updraft qa (passive)', '1' )
   id_qnpcko_han = register_diag_field ( mod_name, &
     'qnpcko_han', axes(full), Time, &
     'Updraft qn (passive)', '1/kg' )
   id_qnipcko_han = register_diag_field ( mod_name, &
     'qnipcko_han', axes(full), Time, &
     'Updraft qni (passive)', '1/kg' )
!
   id_tpcdo_han = register_diag_field ( mod_name, &
     'tpcdo_han', axes(full), Time, &
     'Downdraft temperature (passive)', 'K' )
   id_qvpcdo_han = register_diag_field ( mod_name, &
     'qvpcdo_han', axes(full), Time, &
     'Downdraft qv (passive)', 'kg/kg' )
   id_qlpcdo_han = register_diag_field ( mod_name, &
     'qlpcdo_han', axes(full), Time, &
     'Downdraft ql (passive)', 'kg/kg' )
   id_qipcdo_han = register_diag_field ( mod_name, &
     'qipcdo_han', axes(full), Time, &
     'Downdraft qi (passive)', 'kg/kg' )
   id_qapcdo_han = register_diag_field ( mod_name, &
     'qapcdo_han', axes(full), Time, &
     'Downdraft qa (passive)', '1' )
   id_qnpcdo_han = register_diag_field ( mod_name, &
     'qnpcdo_han', axes(full), Time, &
     'Downdraft qn (passive)', '1/kg' )
   id_qnipcdo_han = register_diag_field ( mod_name, &
     'qnipcdo_han', axes(full), Time, &
     'Downdraft qni (passive)', '1/kg' )
!
   id_tpdt_han = register_diag_field ( mod_name, &
     'tpdt_han', axes(full), Time, &
     'Temperature (passive) tendency by Han EDMF', 'K/s' )
   id_qpdt_han = register_diag_field ( mod_name, &
     'qpdt_han', axes(full), Time, &
     'Specific humidity (passive) tendency by Han EDMF', 'kg/kg/s' )
   id_qlpdt_han = register_diag_field ( mod_name, &
     'qlpdt_han', axes(full), Time, &
     'Cloud liquid (passive) tendency by Han EDMF', 'kg/kg/s' )
   id_qipdt_han = register_diag_field ( mod_name, &
     'qipdt_han', axes(full), Time, &
     'Cloud ice (passive) tendency by Han EDMF', 'kg/kg/s' )
   id_qapdt_han = register_diag_field ( mod_name, &
     'qapdt_han', axes(full), Time, &
     'Cloud amount (passive) tendency by Han EDMF', '1/s' )
   id_qnpdt_han = register_diag_field ( mod_name, &
     'qnpdt_han', axes(full), Time, &
     'Cloud droplet number (passive) tendency by Han EDMF', '1/kg/s' )
   id_qnipdt_han = register_diag_field ( mod_name, &
     'qnipdt_han', axes(full), Time, &
     'Cloud ice number (passive) tendency by Han EDMF', '1/kg/s' )
!
   id_tpdt_mf_han = register_diag_field ( mod_name, &
     'tpdt_mf_han', axes(full), Time, &
     'Temperature (passive) tendency by Han EDMF - Updraft', 'K/s' )
   id_qpdt_mf_han = register_diag_field ( mod_name, &
     'qpdt_mf_han', axes(full), Time, &
     'Specific humidity (passive) tendency by Han EDMF - Updraft', 'kg/kg/s' )
   id_qlpdt_mf_han = register_diag_field ( mod_name, &
     'qlpdt_mf_han', axes(full), Time, &
     'Cloud liquid (passive) tendency by Han EDMF - Updraft', 'kg/kg/s' )
   id_qipdt_mf_han = register_diag_field ( mod_name, &
     'qipdt_mf_han', axes(full), Time, &
     'Cloud ice (passive) tendency by Han EDMF - Updraft', 'kg/kg/s' )
   id_qapdt_mf_han = register_diag_field ( mod_name, &
     'qapdt_mf_han', axes(full), Time, &
     'Cloud amount (passive) tendency by Han EDMF - Updraft', '1/s' )
   id_qnpdt_mf_han = register_diag_field ( mod_name, &
     'qnpdt_mf_han', axes(full), Time, &
     'Cloud liquid number (passive) tendency by Han EDMF - Updraft', '1/kg/s' )
   id_qnipdt_mf_han = register_diag_field ( mod_name, &
     'qnipdt_mf_han', axes(full), Time, &
     'Cloud ice number (passive) tendency by Han EDMF - Updraft', '1/kg/s' )
!
   id_tpdt_mfd_han = register_diag_field ( mod_name, &
     'tpdt_mfd_han', axes(full), Time, &
     'Temperature (passive) tendency by Han EDMF - Downdraft', 'K/s' )
   id_qpdt_mfd_han = register_diag_field ( mod_name, &
     'qpdt_mfd_han', axes(full), Time, &
     'Specific humidity (passive) tendency by Han EDMF - Downdraft', 'kg/kg/s' )
   id_qlpdt_mfd_han = register_diag_field ( mod_name, &
     'qlpdt_mfd_han', axes(full), Time, &
     'Cloud liquid (passive) tendency by Han EDMF - Downdraft', 'kg/kg/s' )
   id_qipdt_mfd_han = register_diag_field ( mod_name, &
     'qipdt_mfd_han', axes(full), Time, &
     'Cloud ice (passive) tendency by Han EDMF - Downdraft', 'kg/kg/s' )
   id_qapdt_mfd_han = register_diag_field ( mod_name, &
     'qapdt_mfd_han', axes(full), Time, &
     'Cloud amount (passive) tendency by Han EDMF - Downdraft', '1/s' )
   id_qnpdt_mfd_han = register_diag_field ( mod_name, &
     'qnpdt_mfd_han', axes(full), Time, &
     'Cloud liquid number (passive) tendency by Han EDMF - Downdraft', '1/kg/s' )
   id_qnipdt_mfd_han = register_diag_field ( mod_name, &
     'qnipdt_mfd_han', axes(full), Time, &
     'Cloud ice number (passive) tendency by Han EDMF - Downdraft', '1/kg/s' )

! ZNT 05/23/2023: Diagnostic Outputs for subs-detr decomposition of MF tendencies for passive T,q,ql,qi
   id_tpdt_sub_han = register_diag_field ( mod_name, &
     'tpdt_sub_han', axes(full), Time, &
     'Temperature (passive) tendency by Han EDMF - Updraft-subsidence', 'K/s' )
   id_qpdt_sub_han = register_diag_field ( mod_name, &
     'qpdt_sub_han', axes(full), Time, &
     'Specific humidity (passive) tendency by Han EDMF - Updraft-subsidence', 'kg/kg/s' )
   id_qlpdt_sub_han = register_diag_field ( mod_name, &
     'qlpdt_sub_han', axes(full), Time, &
     'Cloud liquid (passive) tendency by Han EDMF - Updraft-subsidence', 'kg/kg/s' )
   id_qipdt_sub_han = register_diag_field ( mod_name, &
     'qipdt_sub_han', axes(full), Time, &
     'Cloud ice (passive) tendency by Han EDMF - Updraft-subsidence', 'kg/kg/s' )
   id_qapdt_sub_han = register_diag_field ( mod_name, &
     'qapdt_sub_han', axes(full), Time, &
     'Cloud amount (passive) tendency by Han EDMF - Updraft-subsidence', '1/s' )
   id_qnpdt_sub_han = register_diag_field ( mod_name, &
     'qnpdt_sub_han', axes(full), Time, &
     'Cloud liquid number (passive) tendency by Han EDMF - Updraft-subsidence', '1/kg/s' )
   id_qnipdt_sub_han = register_diag_field ( mod_name, &
     'qnipdt_sub_han', axes(full), Time, &
     'Cloud ice number (passive) tendency by Han EDMF - Updraft-subsidence', '1/kg/s' )
!
   id_tpdt_subd_han = register_diag_field ( mod_name, &
     'tpdt_subd_han', axes(full), Time, &
     'Temperature (passive) tendency by Han EDMF - Downdraft-ascent', 'K/s' )
   id_qpdt_subd_han = register_diag_field ( mod_name, &
     'qpdt_subd_han', axes(full), Time, &
     'Specific humidity (passive) tendency by Han EDMF - Downdraft-ascent', 'kg/kg/s' )
   id_qlpdt_subd_han = register_diag_field ( mod_name, &
     'qlpdt_subd_han', axes(full), Time, &
     'Cloud liquid (passive) tendency by Han EDMF - Downdraft-ascent', 'kg/kg/s' )
   id_qipdt_subd_han = register_diag_field ( mod_name, &
     'qipdt_subd_han', axes(full), Time, &
     'Cloud ice (passive) tendency by Han EDMF - Downdraft-ascent', 'kg/kg/s' )
   id_qapdt_subd_han = register_diag_field ( mod_name, &
     'qapdt_subd_han', axes(full), Time, &
     'Cloud amount (passive) tendency by Han EDMF - Downdraft-ascent', '1/s' )
   id_qnpdt_subd_han = register_diag_field ( mod_name, &
     'qnpdt_subd_han', axes(full), Time, &
     'Cloud liquid number (passive) tendency by Han EDMF - Downdraft-ascent', '1/kg/s' )
   id_qnipdt_subd_han = register_diag_field ( mod_name, &
     'qnipdt_subd_han', axes(full), Time, &
     'Cloud ice number (passive) tendency by Han EDMF - Downdraft-ascent', '1/kg/s' )
!
   id_tpdt_dtr_han = register_diag_field ( mod_name, &
     'tpdt_dtr_han', axes(full), Time, &
     'Temperature (passive) tendency by Han EDMF - Updraft-detrainment', 'K/s' )
   id_qpdt_dtr_han = register_diag_field ( mod_name, &
     'qpdt_dtr_han', axes(full), Time, &
     'Specific humidity (passive) tendency by Han EDMF - Updraft-detrainment', 'kg/kg/s' )
   id_qlpdt_dtr_han = register_diag_field ( mod_name, &
     'qlpdt_dtr_han', axes(full), Time, &
     'Cloud liquid (passive) tendency by Han EDMF - Updraft-detrainment', 'kg/kg/s' )
   id_qipdt_dtr_han = register_diag_field ( mod_name, &
     'qipdt_dtr_han', axes(full), Time, &
     'Cloud ice (passive) tendency by Han EDMF - Updraft-detrainment', 'kg/kg/s' )
   id_qapdt_dtr_han = register_diag_field ( mod_name, &
     'qapdt_dtr_han', axes(full), Time, &
     'Cloud amount (passive) tendency by Han EDMF - Updraft-detrainment', '1/s' )
   id_qnpdt_dtr_han = register_diag_field ( mod_name, &
     'qnpdt_dtr_han', axes(full), Time, &
     'Cloud liquid number (passive) tendency by Han EDMF - Updraft-detrainment', '1/kg/s' )
   id_qnipdt_dtr_han = register_diag_field ( mod_name, &
     'qnipdt_dtr_han', axes(full), Time, &
     'Cloud ice number (passive) tendency by Han EDMF - Updraft-detrainment', '1/kg/s' )
!
   id_tpdt_dtrd_han = register_diag_field ( mod_name, &
     'tpdt_dtrd_han', axes(full), Time, &
     'Temperature (passive) tendency by Han EDMF - Downdraft-detrainment', 'K/s' )
   id_qpdt_dtrd_han = register_diag_field ( mod_name, &
     'qpdt_dtrd_han', axes(full), Time, &
     'Specific humidity (passive) tendency by Han EDMF - Downdraft-detrainment', 'kg/kg/s' )
   id_qlpdt_dtrd_han = register_diag_field ( mod_name, &
     'qlpdt_dtrd_han', axes(full), Time, &
     'Cloud liquid (passive) tendency by Han EDMF - Downdraft-detrainment', 'kg/kg/s' )
   id_qipdt_dtrd_han = register_diag_field ( mod_name, &
     'qipdt_dtrd_han', axes(full), Time, &
     'Cloud ice (passive) tendency by Han EDMF - Downdraft-detrainment', 'kg/kg/s' )
   id_qapdt_dtrd_han = register_diag_field ( mod_name, &
     'qapdt_dtrd_han', axes(full), Time, &
     'Cloud amount (passive) tendency by Han EDMF - Downdraft-detrainment', '1/s' )
   id_qnpdt_dtrd_han = register_diag_field ( mod_name, &
     'qnpdt_dtrd_han', axes(full), Time, &
     'Cloud liquid number (passive) tendency by Han EDMF - Downdraft-detrainment', '1/kg/s' )
   id_qnipdt_dtrd_han = register_diag_field ( mod_name, &
     'qnipdt_dtrd_han', axes(full), Time, &
     'Cloud ice number (passive) tendency by Han EDMF - Downdraft-detrainment', '1/kg/s' )

   ! Add more fields later...

      module_is_initialized=.true.

   end subroutine hanedmf_init

!#######################################################################

   subroutine get_aerosol (zint, pmass, asol,          &
                           am1, am2, am3, am4, am5 ,   &
                           amx1, amx2, amx3, amx4, amx5 )
   
   type(aerosol_type),  intent (in)      :: asol
   
   real,  intent(in),   dimension(:,:,:) :: zint, pmass
   real,  intent(out),  dimension(:,:,:) :: am1, am2, am3, am4, am5
   real,  intent(out),  dimension(:,:,:) :: amx1, amx2, amx3, amx4, amx5
   
   real, dimension(size(am1,1), size(am1,2), size(am1,3))  :: am_tmp
   real, dimension(size(amx1,1),size(amx1,2),size(amx1,3)) :: amx_tmp
   
   integer :: i, j, k, imax, jmax, kmax
   integer :: na, naer
   
   imax  = size( am1, 1 )
   jmax  = size( am1, 2 )
   kmax  = size( am1, 3 )
   naer  = size( asol%aerosol, 4)

   am1=0.; am2=0.; am3=0.; am4=0.; am5=0.;
   amx1=0.; amx2=0.; amx3=0.; amx4=0.; amx5=0.;
   
   if(use_online_aerosol) then

      do k=1,kmax
        am_tmp(:,:,k)  = 1./(zint(:,:,k)-zint(:,:,k+1)) * 1.0e-3
        amx_tmp(:,:,k) = 1./pmass(:,:,k)
      end do

      do na = 1,naer
        if(asol%aerosol_names(na) == 'so4' .or. &
           asol%aerosol_names(na) == 'so4_anthro' .or. &
           asol%aerosol_names(na) == 'so4_natural') then     !aerosol unit: kg/m2
          do k=1,kmax
            do j = 1, jmax
              do i=1, imax
                am1(i,j,k)=am1(i,j,k)+asol%aerosol(i,j,k,na)*am_tmp(i,j,k)  !am1 unit: g/cm3
                amx1(i,j,k)=amx1(i,j,k)+asol%aerosol(i,j,k,na)*amx_tmp(i,j,k)  !amx unit: kg/kg
              end do
            end do
          end do
        else if(asol%aerosol_names(na) == 'omphilic' .or. &
                asol%aerosol_names(na) == 'omphobic') then
          do k=1,kmax
            do j = 1, jmax
              do i=1, imax
                am4(i,j,k)=am4(i,j,k)+asol%aerosol(i,j,k,na)*am_tmp(i,j,k)
                amx4(i,j,k)=amx4(i,j,k)+asol%aerosol(i,j,k,na)*amx_tmp(i,j,k)
              end do
            end do
          end do
        else if(asol%aerosol_names(na) == 'bcphilic' .or. &
                asol%aerosol_names(na) == 'bcphobic' .or. &
                asol%aerosol_names(na) == 'dust1' .or. &
                asol%aerosol_names(na) == 'dust2' .or. &
                asol%aerosol_names(na) == 'dust3' .or. &
                asol%aerosol_names(na) == 'dust_mode1_of_2') then   !h1g, 2015-09-19
          do k=1,kmax
            do j = 1, jmax
              do i=1, imax
                am2(i,j,k)=am2(i,j,k)+asol%aerosol(i,j,k,na)*am_tmp(i,j,k)
                amx2(i,j,k)=amx2(i,j,k)+asol%aerosol(i,j,k,na)*amx_tmp(i,j,k)
              end do
            end do
          end do
        else if(asol%aerosol_names(na) == 'seasalt1' .or. &
                asol%aerosol_names(na) == 'seasalt2' .or. &
                asol%aerosol_names(na) == 'seasalt_aitken' .or. &   !h1g, 2015-09-19
                asol%aerosol_names(na) == 'seasalt_fine'  ) then    !h1g, 2015-09-19
          do k=1,kmax
            do j = 1, jmax
              do i=1, imax
                am3(i,j,k)=am3(i,j,k)+asol%aerosol(i,j,k,na)*am_tmp(i,j,k)
                amx3(i,j,k)=amx3(i,j,k)+asol%aerosol(i,j,k,na)*amx_tmp(i,j,k)
              end do
            end do
          end do
        else if(.not. use_sub_seasalt) then
          if(asol%aerosol_names(na) == 'seasalt3' .or. &
                asol%aerosol_names(na) == 'seasalt4' .or. &
                asol%aerosol_names(na) == 'seasalt5' .or. &
                asol%aerosol_names(na) == 'seasalt_coarse') then    !h1g, 2015-09-19
            do k=1,kmax
              do j = 1, jmax
                do i=1, imax
                  am5(i,j,k)=am5(i,j,k)+asol%aerosol(i,j,k,na)*am_tmp(i,j,k)
                  amx5(i,j,k)=amx5(i,j,k)+asol%aerosol(i,j,k,na)*amx_tmp(i,j,k)
                end do
              end do
            end do
          end if
        end if
      end do

      am2(:,:,:) =am2(:,:,:) +am3(:,:,:) +am4(:,:,:)
      amx2(:,:,:)=amx2(:,:,:)+amx3(:,:,:)+amx4(:,:,:)
      if(.not. use_sub_seasalt) then
        am3(:,:,:)=am3(:,:,:)+am5(:,:,:)
        amx3(:,:,:)=amx3(:,:,:)+amx5(:,:,:)
      end if
    else ! use_online_aerosol
      do k=1,kmax
        do j = 1, jmax
          do i=1, imax
            am2(i,j,k) = asol%aerosol(i,j,k,1)*am_tmp(i,j,k)
            amx2(i,j,k)= asol%aerosol(i,j,k,1)*amx_tmp(i,j,k)
          end do
        end do
      end do
      do k=1,kmax
        do j = 1, jmax
          do i=1, imax
            am1(i,j,k) = asol%aerosol(i,j,k,2)*am_tmp(i,j,k)
            amx1(i,j,k)= asol%aerosol(i,j,k,2)*amx_tmp(i,j,k)
          end do
        end do
      end do
      do k=1,kmax
        do j = 1, jmax
          do i=1, imax
            am4(i,j,k) = om_to_oc*asol%aerosol(i,j,k,3)*am_tmp(i,j,k)
            amx4(i,j,k)= om_to_oc*asol%aerosol(i,j,k,3)*amx_tmp(i,j,k)
          end do
        end do
      end do
      do k=1,kmax
        do j = 1, jmax
          do i=1, imax
            am3(i,j,k) = sea_salt_scale*asol%aerosol(i,j,k,5)*am_tmp(i,j,k)
            amx3(i,j,k)= sea_salt_scale*asol%aerosol(i,j,k,5)*amx_tmp(i,j,k)
          end do
        end do
      end do
    endif ! use_online_aerosol
     
   end subroutine get_aerosol


end module hanedmf_mod

