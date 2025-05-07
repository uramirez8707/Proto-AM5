!> \file satmedmfvdifq.F
!! This file contains the CCPP-compliant SATMEDMF scheme (updated version) which
!! computes subgrid vertical turbulence mixing using scale-aware TKE-based moist
!! eddy-diffusion mass-flux (TKE-EDMF) parameterization (by Jongil Han).

    module satmedmfvdifq_mod

    contains

!> \defgroup satmedmfvdifq GFS Scale-aware TKE-based Moist Eddy-Diffusivity Mass-flux (TKE-EDMF, updated version) Scheme Module

!! \brief This subroutine contains all of the logic for the
!! scale-aware TKE-based moist eddy-diffusion mass-flux (TKE-EDMF, updated version) scheme.
!! For local turbulence mixing, a TKE closure model is used.
!! Updated version of satmedmfvdif.f (May 2019) to have better low level
!! inversion, to reduce the cold bias in lower troposphere,
!! and to reduce the negative wind speed bias in upper troposphere

!> \section arg_table_satmedmfvdifq_init Argument Table
!! \htmlinclude satmedmfvdifq_init.html
!!
    subroutine satmedmfvdifq_init (isatmedmf,isatmedmf_vdifq, &
                                   errmsg,errflg)

    integer, intent(in) :: isatmedmf,isatmedmf_vdifq
    character(len=*), intent(out) :: errmsg
    integer,          intent(out) :: errflg

! Initialize CCPP error handling variables
    errmsg = ''
    errflg = 0

    if ( .NOT. isatmedmf==isatmedmf_vdifq) then
        write(errmsg,fmt='(*(a))') 'Logic error: satmedmfvdif is ', &
        'called, but isatmedmf/=isatmedmf_vdifq.'
        errflg = 1
        return
    end if

    end subroutine satmedmfvdifq_init

    subroutine satmedmfvdifq_finalize ()
    end subroutine satmedmfvdifq_finalize

!> \section arg_table_satmedmfvdifq_run Argument Table
!! \htmlinclude satmedmfvdifq_run.html
!!
!!\section gen_satmedmfvdifq GFS satmedmfvdifq General Algorithm
!! satmedmfvdifq_run() computes subgrid vertical turbulence mixing
!! using the scale-aware TKE-based moist eddy-diffusion mass-flux (EDMF) parameterization of
!! Han and Bretherton (2019) \cite Han_2019 .
!! -# The local turbulent mixing is represented by an eddy-diffusivity scheme which
!! is a function of a prognostic TKE.
!! -# For the convective boundary layer, nonlocal transport by large eddies
!! (mfpbltq.f), is represented using a mass flux approach (Siebesma et al.(2007) \cite Siebesma_2007 ).
!! -# A mass-flux approach is also used to represent the stratocumulus-top-induced turbulence
!! (mfscuq.f).
!! \section detail_satmedmfvidfq GFS satmedmfvdifq Detailed Algorithm

    subroutine satmedmfvdifq_run(ix,im,km,ntrac,ntcw,ntiw,ntaw,  & ! &  ZNT: 10/11/2020
    ntnw,ntniw,ntqpw,ntke,                                       & ! &  ZNT: 05/12-19/2023, add qn/qni; passive T,q,ql,qi
    grav,rd,cp,rv,hvap,hfus,fv,eps,epsm1,                        & 
    dv,du,tdt,rtg,u1,v1,t1,q1,aer1,aerx1,swh,hlw,xmu,garea,zvfun,& ! &  ZNT: 06/28/2023, add aer1 and aerx1
    cfl_crit,do_cfl_col,ck0,ck1,ch0,ch1,ce0,do_cmp_qs,use_hl,    & ! &  ZNT: 09/24-28/2020; 06/06/2021; 06/12/2023
    do_mf_aer,tkeh_flg,cldtime,do_upd_ovrsht,c_entu,c_entd,actei,& ! &  ZNT: 11/11-18/2020; 06/21/2021; 06/28/2023 do_mf_aer
    do_tkevdif, mf_tkeflx_fac, do_tkediss_imp,                   & ! &  ZNT: 06/05-06/2021
    mf_buop_fac, mf_shrp_fac, do_pos_mf_buop,                    & ! &  ZNT: 06/05/2021
    do_km_L65, do_wush, do_tkemean,                              & ! &  ZNT: 10/02/2023; 12/24/2023
    psk,rbsoil,zorl,u10m,v10m,fm,fh,                             & 
    tsea,heat,evap,stress,spd1,kpbl,                             & 
    prsi,del,prsl,prslk,phii,phil,delt,                          & 
    dspheat,dusfc,dvsfc,dtsfc,dqsfc,hpbl,                        & 
    kinver,xkzm_m,xkzm_h,xkzm_s,dspfac,bl_upfr,bl_dnfr,          & 
    dku, dkt, dkq, xmf, xmfd,                                    & 
    tcko_o,qvcko_o,qlcko_o,qicko_o,qacko_o,qncko_o,qnicko_o,     & ! &  ZNT: 06/19/2020; qn/qni added 05/03/2023; qa added 05/24/2023
    tqpcko_o,buou_o,                                             & ! &  ZNT: 05/19/2023: passive T,q,ql,qi
    tcdo_o,qvcdo_o,qlcdo_o,qicdo_o,qacdo_o,qncdo_o,qnicdo_o,     & ! &  ZNT: 06/19/2020; qn/qni added 05/03/2023; qa added 05/24/2023
    tqpcdo_o,buod_o,                                             & ! &  ZNT: 05/19/2023: passive T,q,ql,qi
    dv_mf,du_mf,tdt_mf,qdt_mf,qldt_mf,qidt_mf,qadt_mf,           & 
    qndt_mf,qnidt_mf,tqpdt_mf,                                   & ! &  ZNT: 05/12-19/2023
    dv_mfd,du_mfd,tdt_mfd,qdt_mfd,qldt_mfd,qidt_mfd,qadt_mfd,    & 
    qndt_mfd,qnidt_mfd,tqpdt_mfd,                                & ! &  ZNT: 05/12-19/2023
    dusfc_mf,dvsfc_mf,dtsfc_mf,dqsfc_mf,                         & 
    dusfc_mfd,dvsfc_mfd,dtsfc_mfd,dqsfc_mfd,                     & 
    tdt_sub,qdt_sub,qldt_sub,qidt_sub,                           & ! &  ZNT: 05/24/2023, subs-detr
    qadt_sub,qndt_sub,qnidt_sub,                                 & ! &  ZNT: 05/24/2023
    tdt_subd,qdt_subd,qldt_subd,qidt_subd,                       & ! &  ZNT: 05/24/2023
    qadt_subd,qndt_subd,qnidt_subd,                              & ! &  ZNT: 05/24/2023
    tdt_dtr,qdt_dtr,qldt_dtr,qidt_dtr,                           & ! &  ZNT: 05/24/2023, subs-detr
    qadt_dtr,qndt_dtr,qnidt_dtr,                                 & ! &  ZNT: 05/24/2023
    tdt_dtrd,qdt_dtrd,qldt_dtrd,qidt_dtrd,                       & ! &  ZNT: 05/24/2023
    qadt_dtrd,qndt_dtrd,qnidt_dtrd,                              & ! &  ZNT: 05/24/2023
    tqpdt_sub, tqpdt_subd, tqpdt_dtr, tqpdt_dtrd,                & ! &  ZNT: 05/23/2023
    de_buop,de_shrp,de_diss,                                     & ! &  ZNT: 10/03/2020
    errmsg,errflg)

    use sat_vapor_pres_mod, only : lookup_es, compute_qs                ! ZNT 04/27/2020
!      use funcphys , only : fpvs
    use tridi_mod

    implicit none

!----------------------------------------------------------------------
    integer, intent(in)  :: ix, im, km, ntrac, ntcw, ntiw, ntaw, ntke   ! ZNT: 10/11/2020
    integer, intent(in)  :: ntnw, ntniw, ntqpw                          ! ZNT: 05/12-19/2023
    integer, intent(in)  :: kinver(im)
    integer, intent(out) :: kpbl(im)

    real   , intent(in)  :: grav,rd,cp,rv,hvap,hfus,fv,eps,epsm1
    real   , intent(in)  :: cfl_crit,ck0,ck1,ch0,ch1,ce0    ! ZNT 09/24-28/2020; 06/06/2021
    logical, intent(in)  :: do_cfl_col                                  ! ZNT 06/06/2021
    real   , intent(in)  :: cldtime                         ! ZNT 11/11/2020
    real   , intent(in)  :: c_entu, c_entd, actei           ! ZNT 06/21/2021
    real   , intent(in)  :: delt, xkzm_m, xkzm_h, xkzm_s
    real   , intent(in)  :: dspfac, bl_upfr, bl_dnfr
    real   , intent(inout) :: dv(im,km),  du(im,km),  & 
                              tdt(im,km), rtg(im,km,ntrac)
    real   , intent(in)  :: aer1(im,km,4), aerx1(im,km,4)   ! ZNT: 06/28/2023, add aer1 and aerx1
    real   , intent(in)  :: u1(ix,km),     v1(ix,km),        &
                            t1(ix,km),     q1(ix,km,ntrac),  &
                            swh(ix,km),    hlw(ix,km),       &
                            xmu(im),       garea(im),        &
                            zvfun(im),                       &
                            psk(ix),       rbsoil(im),       &
                            zorl(im),      tsea(im),         &
                            u10m(im),      v10m(im),         &
                            fm(im),        fh(im),           &
                            evap(im),      heat(im),         &
                            stress(im),    spd1(im),         &
                            prsi(ix,km+1), del(ix,km),       &
                            prsl(ix,km),   prslk(ix,km),     &
                            phii(ix,km+1), phil(ix,km)
    real   , intent(out) :: dusfc(im),     dvsfc(im),        &
                            dtsfc(im),     dqsfc(im),        &
                            hpbl(im)

! ZNT: MF tendencies
    real   , intent(out) :: dv_mf(im,km),   du_mf(im,km),    &
                            dv_mfd(im,km),  du_mfd(im,km),   &
                            tdt_mf(im,km),  qdt_mf(im,km),   &
                            tdt_mfd(im,km), qdt_mfd(im,km),  &
                            qldt_mf(im,km), qidt_mf(im,km),  &
                            qadt_mf(im,km), qndt_mf(im,km),  &
                            qnidt_mf(im,km),                 &
                            qldt_mfd(im,km),qidt_mfd(im,km), &
                            qadt_mfd(im,km),qndt_mfd(im,km), &
                            qnidt_mfd(im,km)
    real   , intent(out) :: dusfc_mf(im),   dvsfc_mf(im),    &
                            dtsfc_mf(im),   dqsfc_mf(im),    &
                            dusfc_mfd(im),  dvsfc_mfd(im),   &
                            dtsfc_mfd(im),  dqsfc_mfd(im)
! ZNT: End of MF tendencies

! ZNT: 05/24/2023, 07/05/2023 - subs-detr decomposition of MF tendencies
    real   , intent(out) :: tdt_sub(im,km),   qdt_sub(im,km),   &
                            tdt_subd(im,km),  qdt_subd(im,km),  &
                            qldt_sub(im,km),  qidt_sub(im,km),  &
                            qadt_sub(im,km),  qndt_sub(im,km),  &
                            qnidt_sub(im,km),                   &
                            qldt_subd(im,km), qidt_subd(im,km), &
                            qadt_subd(im,km), qndt_subd(im,km), &
                            qnidt_subd(im,km)
    real   , intent(out) :: tdt_dtr(im,km),   qdt_dtr(im,km),   &
                            tdt_dtrd(im,km),  qdt_dtrd(im,km),  &
                            qldt_dtr(im,km),  qidt_dtr(im,km),  &
                            qadt_dtr(im,km),  qndt_dtr(im,km),  &
                            qnidt_dtr(im,km),                   &
                            qldt_dtrd(im,km), qidt_dtrd(im,km), &
                            qadt_dtrd(im,km), qndt_dtrd(im,km), &
                            qnidt_dtrd(im,km)
    real   , intent(out) :: tqpdt_sub(im,km,7), tqpdt_subd(im,km,7), & 
                            tqpdt_dtr(im,km,7), tqpdt_dtrd(im,km,7)
! ZNT: End of subs-detr tendencies
         
! ZNT: TKE tendencies, 10/03/2020
    real   , intent(out) ::  de_buop(im,km),             & 
                             de_shrp(im,km), de_diss(im,km)


    real   , intent(out) ::  dku(im,km-1), dkt(im,km-1), &
                             dkq(im,km-1),               & 
                             xmf(im,km),   xmfd(im,km)

! ZNT: plume diagnostics for T/qv/ql/qi/b, 06/19/2020
    real   , intent(out) ::  tcko_o(im,km), qvcko_o(im,km),  &
                             qlcko_o(im,km), qicko_o(im,km), &
                             buou_o(im,km),  &
                             tcdo_o(im,km), qvcdo_o(im,km),  &
                             qlcdo_o(im,km), qicdo_o(im,km), &
                             buod_o(im,km)

! ZNT: plume diagnostics for qa/qn/qni, 05/12/2023; 05/24/2023
    real   , intent(out) ::  qacko_o(im,km), qncko_o(im,km), &
                             qnicko_o(im,km),  &
                             qacdo_o(im,km), qncdo_o(im,km), &
                             qnicdo_o(im,km)

! ZNT: plume diagnostics for passive variables, 05/19/2023; 07/05/2023
    real   , intent(out) ::  tqpcko_o(im,km,7), tqpcdo_o(im,km,7),   & 
                             tqpdt_mf(im,km,7), tqpdt_mfd(im,km,7)

    logical, intent(in)  :: dspheat
    integer, intent(in)  :: do_mf_aer                   ! ZNT: 06/30/2023, option for aerosol activation in updrafts
                                                        !      0 = skip; 1 = 'do_online_aerosol = .false.'
                                                        !      2 = 'do_online_aerosol = .true.'
    logical, intent(in)  :: do_cmp_qs                   ! ZNT: 10/09/2020
    integer, intent(in)  :: use_hl                      ! ZNT: 06/12/2023: 0=no; 1=Tu; 2=bu; 3=both
    integer, intent(in)  :: tkeh_flg                    ! ZNT: 10/13/2020
    logical, intent(in)  :: do_upd_ovrsht               ! ZNT: 11/18/2020
    logical, intent(in)  :: do_tkevdif                  ! ZNT: 06/05/2021
    real   , intent(in)  :: mf_tkeflx_fac               ! ZNT: 06/06/2021
    logical, intent(in)  :: do_tkediss_imp              ! ZNT: 06/05/2021
    real   , intent(in)  :: mf_buop_fac, mf_shrp_fac    ! ZNT: 06/05/2021
    logical, intent(in)  :: do_pos_mf_buop              ! ZNT: 06/05/2021
    logical, intent(in)  :: do_km_L65                   ! ZNT: 10/02/2023
    logical, intent(in)  :: do_wush                     ! ZNT: 12/24/2023
    logical, intent(in)  :: do_tkemean                  ! ZNT: 12/24/2023

    character(len=*), intent(out) :: errmsg
    integer,          intent(out) :: errflg

!          flag for tke dissipative heating

!----------------------------------------------------------------------
!***
!***  local variables
!***
    integer :: i,is,k,kk,ktqp,n,ndt,km1,kmpbl,kmscu,ntrac1
    integer :: kps                                      ! ZNT: 12/24/2023
    integer :: lcld(im),kcld(im),krad(im),mrad(im)
    integer :: kx1(im), kpblx(im)

    real    :: tke(im,km),  tkeh(im,km-1)

    real    :: theta(im,km),thvx(im,km),  thlvx(im,km),  &
               qlx(im,km),  thetae(im,km),thlx(im,km),   &
               slx(im,km),  svx(im,km),   qtx(im,km),    &
               tvx(im,km),  pix(im,km),   radx(im,km-1), &
               cku(im,km-1),ckt(im,km-1)
! ZNT: dku(im,km-1),dkt(im,km-1), dkq(im,km-1) are now output

    real    :: plyr(im,km), rhly(im,km),  cfly(im,km), &
               qstl(im,km)

    real    :: dtdz1(im), gdx(im),   &
               phih(im),  phim(im),  prn(im,km-1), &
               rbdn(im),  rbup(im),  thermal(im),  &
               ustar(im), wstar(im), hpblx(im),    &
               ust3(im),  wst3(im),  &
               z0(im),    crb(im),   tkemean(im),  &   ! ZNT 12/24/2023: tkemean
               hgamt(im), hgamq(im), &
               wscale(im),vpert(im), &
               zol(im),   sflux(im), &
               sumx(im),  tx1(im),   tx2(im)           ! ZNT 12/24/2023: sumx

    real    :: radmin(im)

    real    :: zi(im,km+1),  zl(im,km),   zm(im,km),   &
               xkzo(im,km),  xkzmo(im,km),   &
               xkzm_hx(im),  xkzm_mx(im),    &
               ri(im,km-1),  tkmnz(im,km-1), &
               rdzt(im,km-1),rlmnz(im,km),   &
               al(im,km-1),  ad(im,km),   au(im,km-1), &
               f1(im,km),    f2(im,km*(ntrac-1))

    real    :: elm(im,km),   ele(im,km), &
               ckz(im,km),   chz(im,km), &
               diss(im,km-1),prod(im,km-1), &
               bf(im,km-1),  shr2(im,km-1), wush(im,km), & ! ZNT 12/24/2023: wush
               xlamue(im,km-1), xlamde(im,km-1), &
               gotvx(im,km), rlam(im,km-1)

!   variables for updrafts (thermals)

    real    :: tcko(im,km),  qcko(im,km,ntrac), &
               ucko(im,km),  vcko(im,km), &
               buou(im,km)  ! xmf(im,km) is an output now

!   variables for stratocumulus-top induced downdrafts

    real    :: tcdo(im,km),  qcdo(im,km,ntrac), &
               ucdo(im,km),  vcdo(im,km), &
               buod(im,km) ! xmfd(im,km) is an output now

    logical :: pblflg(im), sfcflg(im), flg(im)
    logical :: scuflg(im), pcnvflg(im)
    logical :: mlenflg

!  pcnvflg: true for unstable pbl

    real    :: aphi16,  aphi5, &
               wfac,    cfac,  &
               gamcrt,  gamcrq, sfcfrac, &
               conq,    cont,   conw,    &
               dsdz2,   dsdzt,  dkmax,   &
               dsig,    dt2,    dtodsd,  &
               dtodsu,  g,      factor, dz,     &
               gocp,    gravi,  zol1,   zolcru, &
               buop,    shrp,   dtn,     &
               prnum,   prmax,  prmin,  prtke,  &
               prscu,   pr0,   &
               dw2,     dw2min, zk,      &
               elmfac,  elefac, dspmax,  &
               alp,     clwt,   cql,     &
               f0,      robn,   crbmin, crbmax, &
               es,      qs,     value,  onemrh, &
               cfh,     gamma,  elocp,  el2orc, &
               epsi,    beta,   chx,    cqx,    &
               rdt,     rdz,    qmin,   qlmin,  &
               rimin,   rbcr,   rbint,  tdzmin, &
               rlmn,    rlmn1,  rlmn2,   &
               rlmx,    elmx,  &
               ttend,   utend,  vtend,  qtend,  &
               zfac,    zfmin,  vk,     spdk2,  &
               tkmin,   tkminx, xkgdx,  xkinv,  &
               zlup,    zldn,   bsum, cs0,csmf, &  ! ZNT 12/24/2023: csmf
               tem,     tem1,   tem2,   tem3,   &
               ptem,    ptem0,  ptem1,  ptem2

!   real    :: ck0, ck1, ch0, ch1, ce0, rchck
    real    :: rchck

    real    :: qlcr, zstblmax, hcrinv

    real    :: h1
!!
    parameter(wfac=7.0,cfac=4.5)
    parameter(gamcrt=3.,gamcrq=0.,sfcfrac=0.1)
    parameter(vk=0.4,rimin=-100.)
    parameter(rbcr=0.25,zolcru=-0.02,tdzmin=1.e-3)
    parameter(rlmn=30.,rlmn1=5.,rlmn2=10.)
    parameter(rlmx=300.,elmx=300.)
    parameter(prmin=0.25,prmax=4.0)
    parameter(pr0=1.0,prtke=1.0,prscu=0.67)
    parameter(f0=1.e-4,crbmin=0.15,crbmax=0.35)
    parameter(tkmin=1.e-9,tkminx=0.2,dspmax=10.0)
    parameter(qmin=1.e-8,qlmin=1.e-12,zfmin=1.e-8)
    parameter(aphi5=5.,aphi16=16.)
    parameter(elmfac=1.0,elefac=1.0,cql=100.)
    parameter(dw2min=1.e-4,dkmax=1000.,xkgdx=5000.)
    parameter(qlcr=3.5e-5,zstblmax=2500.)
    parameter(xkinv=0.3)
    parameter(h1=0.33333333,hcrinv=250.)
!     parameter(ck0=0.4,ck1=0.15,ch0=0.4,ch1=0.15,ce0=0.4)
!     parameter(cs0=0.5)
    parameter(cs0=0.2,csmf=0.5)
    parameter(rchck=1.5,ndt=20)

    gravi=1.0/grav
    g=grav
    gocp=g/cp
    cont=cp/g
    conq=hvap/g
    conw=1.0/g  ! for del in pa
!     parameter(cont=1000.*cp/g,conq=1000.*hvap/g,conw=1000./g) !kpa
    elocp=hvap/cp
    el2orc=hvap*hvap/(rv*cp)

!************************************************************************
! Initialize CCPP error handling variables
    errmsg = ''
    errflg = 0

!> ## Compute preliminary variables from input arguments
    dt2  = delt
    rdt = 1. / dt2

! the code is written assuming ntke=ntrac
! if ntrac > ntke, the code needs to be modified

    ntrac1 = ntrac - 1
    km1 = km - 1
    if (do_km_L65) then
! ZNT 10/04/2023: corrected from 2/3 to 1/3.  
       kmpbl = km / 3
       kmscu = km / 3
    else
       kmpbl = km / 2
       kmscu = km / 2
    endif
!>  - Compute physical height of the layer centers and interfaces from
!! the geopotential height (\p zi and \p zl)
    do k=1,km
        do i=1,im
            zi(i,k) = phii(i,k) * gravi
            zl(i,k) = phil(i,k) * gravi
            xmf(i,k) = 0.
            xmfd(i,k) = 0.
            buou(i,k) = 0.
            buod(i,k) = 0.
            wush(i,k) = 0.
            ckz(i,k) = ck1
            chz(i,k) = ch1
            rlmnz(i,k) = rlmn
        enddo
    enddo
    do i=1,im
        zi(i,km+1) = phii(i,km+1) * gravi
    enddo
    do k=1,km
        do i=1,im
            zm(i,k) = zi(i,k+1)
        enddo
    enddo
!>  - Compute horizontal grid size (\p gdx)
    do i=1,im
        gdx(i) = sqrt(garea(i))
    enddo
!>  - Initialize tke value at vertical layer centers and interfaces
!! from tracer (\p tke and \p tkeh)
    do k=1,km
        do i=1,im
            tke(i,k) = max(q1(i,k,ntke), tkmin)
        enddo
    enddo

    if (tkeh_flg == 1) then      ! ZNT 10/13/2020, use geometric mean for TKEh
        do k=1,km1
            do i=1,im
                tkeh(i,k) = sqrt(tke(i,k)*tke(i,k+1))
            enddo
        enddo
    elseif (tkeh_flg == 2) then  ! ZNT 10/13/2020, use harmonic mean for TKEh
        do k=1,km1
            do i=1,im
                tkeh(i,k) = 2.0/(1./tke(i,k)+1./tke(i,k+1))
            enddo
        enddo
    else  ! ZNT 10/13/2020, use simple mean for TKEh by default
        do k=1,km1
            do i=1,im
                tkeh(i,k) = 0.5 * (tke(i,k) + tke(i,k+1))
            enddo
        enddo
    endif
!>  - Compute reciprocal of \f$ \Delta z \f$ (rdzt)
    do k = 1,km1
        do i=1,im
            rdzt(i,k) = 1.0 / (zl(i,k+1) - zl(i,k))
            prn(i,k)  = pr0
        enddo
    enddo

!>  - Compute reciprocal of pressure (tx1, tx2)

!>  - Compute minimum turbulent mixing length (rlmnz)

!>  - Compute background vertical diffusivities for scalars and momentum (xkzo and xkzmo)

!>  - set background diffusivities as a function of
!!    horizontal grid size with xkzm_h & xkzm_m for gdx >= 25km
!!    and 0.01 for gdx=5m, i.e.,
!!    \n  xkzm_hx = 0.01 + (xkzm_h - 0.01)/(xkgdx-5.) * (gdx-5.)
!!    \n  xkzm_mx = 0.01 + (xkzm_m - 0.01)/(xkgdx-5.) * (gdx-5.)

    do i=1,im
        kx1(i) = 1
        tx1(i) = 1.0 / prsi(i,1)
        tx2(i) = tx1(i)
        if(gdx(i) >= xkgdx) then
            xkzm_hx(i) = xkzm_h
            xkzm_mx(i) = xkzm_m
        else
            tem  = 1. / (xkgdx - 5.)
            tem1 = (xkzm_h - 0.01) * tem
            tem2 = (xkzm_m - 0.01) * tem
            ptem = gdx(i) - 5.
            xkzm_hx(i) = 0.01 + tem1 * ptem
            xkzm_mx(i) = 0.01 + tem2 * ptem
        endif
    enddo
    do k = 1,km
        do i=1,im
            xkzo(i,k)  = 0.0
            xkzmo(i,k) = 0.0
            if (k < kinver(i)) then
            !                               minimum turbulent mixing length
                ptem      = prsl(i,k) * tx1(i)
                tem1      = 1.0 - ptem
                tem2      = tem1 * tem1 * 2.5
                tem2      = min(1.0, exp(-tem2))
                rlmnz(i,k)= rlmn * tem2
                rlmnz(i,k)= max(rlmnz(i,k), rlmn1)
            !                               vertical background diffusivity
                tem2      = tem1 * tem1 * 10.0
                tem2      = min(1.0, exp(-tem2))
                xkzo(i,k) = xkzm_hx(i) * tem2
            !                               vertical background diffusivity for momentum
                if (ptem >= xkzm_s) then
                    xkzmo(i,k) = xkzm_mx(i)
                    kx1(i)     = k + 1
                else
                    if (k == kx1(i) .AND. k > 1) tx2(i) = 1.0 / prsi(i,k)
                    tem1 = 1.0 - prsl(i,k+1) * tx2(i)
                    tem1 = tem1 * tem1 * 5.0
                    xkzmo(i,k) = xkzm_mx(i) * min(1.0, exp(-tem1))
                endif
            endif
        enddo
    enddo

!>  - Some output variables and logical flags are initialized
    do i = 1,im
        z0(i)    = 0.01 * zorl(i)
        dusfc(i) = 0.
        dvsfc(i) = 0.
        dtsfc(i) = 0.
        dqsfc(i) = 0.
        kpbl(i) = 1
        hpbl(i) = 0.
        kpblx(i) = 1
        hpblx(i) = 0.
        pblflg(i)= .true.
        sfcflg(i)= .true.
        if(rbsoil(i) > 0.) sfcflg(i) = .FALSE. 
        pcnvflg(i)= .false.
        scuflg(i)= .true.
        if(scuflg(i)) then
            radmin(i)= 0.
            mrad(i)  = km1
            krad(i)  = 1
            lcld(i)  = km1
            kcld(i)  = km1
        endif
    enddo

!>  - Compute \f$\theta\f$(theta), and \f$q_l\f$(qlx), \f$\theta_e\f$(thetae),
!! \f$\theta_v\f$(thvx),\f$\theta_{l,v}\f$ (thlvx) including ice water
    do k=1,km
        do i=1,im
            pix(i,k)   = psk(i) / prslk(i,k)
            theta(i,k) = t1(i,k) * pix(i,k)
            if(ntiw > 0) then
                tem = max(q1(i,k,ntcw),qlmin)
                tem1 = max(q1(i,k,ntiw),qlmin)
                qlx(i,k) = tem + tem1
                ptem = hvap*tem + (hvap+hfus)*tem1
                slx(i,k)   = cp * t1(i,k) + phil(i,k) - ptem
            else
                qlx(i,k) = max(q1(i,k,ntcw),qlmin)
                slx(i,k)   = cp * t1(i,k) + phil(i,k) - hvap*qlx(i,k)
            endif
            tem2       = 1.+fv*max(q1(i,k,1),qmin)-qlx(i,k)
            thvx(i,k)  = theta(i,k) * tem2
            tvx(i,k)   = t1(i,k) * tem2
            qtx(i,k) = max(q1(i,k,1),qmin)+qlx(i,k)
            thlx(i,k)  = theta(i,k) - pix(i,k)*elocp*qlx(i,k)
            thlvx(i,k) = thlx(i,k) * (1. + fv * qtx(i,k))
            svx(i,k)   = cp * tvx(i,k)
            ptem1      = elocp * pix(i,k) * max(q1(i,k,1),qmin)
            thetae(i,k)= theta(i,k) +  ptem1
            gotvx(i,k) = g / tvx(i,k)
        enddo
    enddo

!>  - Compute an empirical cloud fraction based on
!! Xu and Randall (1996) \cite xu_and_randall_1996
    do k = 1, km
        do i = 1, im
            plyr(i,k)   = 0.01 * prsl(i,k)   ! pa to mb (hpa)
        !  --- ...  compute relative humidity
            if (do_cmp_qs) then  ! ZNT 10/09/2020: note es is not used
                call compute_qs(t1(i,k), prsl(i,k), qs)
                qs = max(qmin,qs)
            else
                call lookup_es(t1(i,k), es)
                es  = 0.01 * es             ! ZNT 04/27/2020
            ! es  = 0.01 * fpvs(t1(i,k))       ! fpvs in pa
                qs  = max(qmin, eps * es / (plyr(i,k) + epsm1*es))
            endif
            rhly(i,k) = max(0.0, min(1.0, max(qmin, q1(i,k,1))/qs))
            qstl(i,k) = qs
        enddo
    enddo

    do k = 1, km
        do i = 1, im
            cfly(i,k) = 0.
            clwt = 1.0e-6 * (plyr(i,k)*0.001)
            if (qlx(i,k) > clwt) then
                onemrh= max(1.e-10, 1.0-rhly(i,k))
                tem1  = min(max((onemrh*qstl(i,k))**0.49,0.0001),1.0)
                tem1  = cql / tem1
                value = max(min( tem1*qlx(i,k), 50.0), 0.0)
                tem2  = sqrt(sqrt(rhly(i,k)))
                cfly(i,k) = min(max(tem2*(1.0-exp(-value)), 0.0), 1.0)
            endif
        enddo
    enddo

!>  - Compute buoyancy modified by clouds

    do k = 1, km1
        do i = 1, im
            tem  = 0.5 * (svx(i,k) + svx(i,k+1))
            tem1 = 0.5 * (t1(i,k) + t1(i,k+1))
            tem2 = 0.5 * (qstl(i,k) + qstl(i,k+1))
            cfh  = min(cfly(i,k+1),0.5*(cfly(i,k)+cfly(i,k+1)))
            alp  = g / tem
            gamma = el2orc * tem2 / (tem1**2)
            epsi  = tem1 / elocp
            beta  = (1. + gamma*epsi*(1.+fv)) / (1. + gamma)
            chx   = cfh * alp * beta + (1. - cfh) * alp
            cqx   = cfh * alp * hvap * (beta - epsi)
            cqx   = cqx + (1. - cfh) * fv * g
            ptem1 = (slx(i,k+1)-slx(i,k))*rdzt(i,k)
            ptem2 = (qtx(i,k+1)-qtx(i,k))*rdzt(i,k)
            bf(i,k) = chx * ptem1 + cqx * ptem2
        enddo
    enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!>  - Initialize diffusion coefficients to 0 and calculate the total
!! radiative heating rate (dku, dkt, radx)
    do k=1,km1
        do i=1,im
            dku(i,k)  = 0.
            dkt(i,k)  = 0.
            dkq(i,k)  = 0.
            cku(i,k)  = 0.
            ckt(i,k)  = 0.
            tem       = zi(i,k+1)-zi(i,k)
            radx(i,k) = tem*(swh(i,k)*xmu(i)+hlw(i,k))
        enddo
    enddo
!>  - Compute stable/unstable PBL flag (pblflg) based on the total
!! surface energy flux (\e false if the total surface energy flux
!! is into the surface)
    do i = 1,im
        sflux(i)  = heat(i) + evap(i)*fv*theta(i,1)
        if( .NOT. sfcflg(i) .OR. sflux(i) <= 0.) pblflg(i)= .FALSE. 
    enddo

!>  ## Calculate the PBL height
!!  The calculation of the boundary layer height follows Troen and Mahrt (1986)
!!  \cite troen_and_mahrt_1986 section 3. The approach is to find the level in 
!!  the column where a modified bulk Richardson number exceeds a critical value.
!!  - Compute critical bulk Richardson number (\f$Rb_{cr}\f$) (crb)
!!  - For the unstable PBL, crb is a constant (0.25)
!!  - For the stable boundary layer (SBL), \f$Rb_{cr}\f$ varies
!!     with the surface Rossby number, \f$R_{0}\f$, as given by
!!     Vickers and Mahrt (2004) \cite Vickers_2004
!!     \f[
!!      Rb_{cr}=0.16(10^{-7}R_{0})^{-0.18}
!!     \f]
!!     \f[
!!      R_{0}=\frac{U_{10}}{f_{0}z_{0}}
!!     \f]
!!     where \f$U_{10}\f$ is the wind speed at 10m above the ground surface,
!!     \f$f_0\f$ the Coriolis parameter, and \f$z_{0}\f$ the surface roughness
!!     length. To avoid too much variation, we restrict \f$Rb_{cr}\f$ to vary
!!     within the range of 0.15~0.35
    do i = 1,im
        if(pblflg(i)) then
        !         thermal(i) = thvx(i,1)
            thermal(i) = thlvx(i,1)
            crb(i) = rbcr
        else
            thermal(i) = tsea(i)*(1.+fv*max(q1(i,1,1),qmin))
            tem = sqrt(u10m(i)**2+v10m(i)**2)
            tem = max(tem, 1.)
            robn = tem / (f0 * z0(i))
            tem1 = 1.e-7 * robn
            crb(i) = 0.16 * (tem1 ** (-0.18))
            crb(i) = max(min(crb(i), crbmax), crbmin)
        endif
    enddo
!>  - Compute \f$\frac{\Delta t}{\Delta z}\f$ , \f$u_*\f$
    do i=1,im
        dtdz1(i)  = dt2 / (zi(i,2)-zi(i,1))
    enddo

    do i=1,im
        ustar(i) = sqrt(stress(i))
    enddo

!>  - Compute buoyancy \f$\frac{\partial \theta_v}{\partial z}\f$ (bf)
!! and the wind shear squared (shr2)

    do k = 1, km1
        do i = 1, im
            rdz  = rdzt(i,k)
        !        bf(i,k) = gotvx(i,k)*(thvx(i,k+1)-thvx(i,k))*rdz
            dw2  = (u1(i,k)-u1(i,k+1))**2 &
            + (v1(i,k)-v1(i,k+1))**2
            shr2(i,k) = max(dw2,dw2min)*rdz*rdz
            ri(i,k) = max(bf(i,k)/shr2(i,k),rimin)
        enddo
    enddo

! Find pbl height based on bulk richardson number (mrf pbl scheme)
!   and also for diagnostic purpose

    do i=1,im
        flg(i) = .false.
        rbup(i) = rbsoil(i)
    enddo
!>  - Given the thermal's properties and the critical Richardson number,
!! a loop is executed to find the first level above the surface (kpblx) where
!! the modified Richardson number is greater than the critical Richardson
!! number, using equation 10a from Troen and Mahrt (1996) \cite troen_and_mahrt_1986
!! (also equation 8 from Hong and Pan (1996) \cite hong_and_pan_1996):
    do k = 1, kmpbl
        do i = 1, im
            if( .NOT. flg(i)) then
                rbdn(i) = rbup(i)
                spdk2   = max((u1(i,k)**2+v1(i,k)**2),1.)
            !         rbup(i) = (thvx(i,k)-thermal(i))*
            !    &              (g*zl(i,k)/thvx(i,1))/spdk2
                rbup(i) = (thlvx(i,k)-thermal(i))* &
                (g*zl(i,k)/thlvx(i,1))/spdk2
                kpblx(i) = k
                flg(i)  = rbup(i) > crb(i)
            endif
        enddo
    enddo
!>  - Once the level is found, some linear interpolation is performed to find
!! the exact height of the boundary layer top (where \f$R_{i} > Rb_{cr}\f$)
!! and the PBL height (hpbl and kpbl) and the PBL top index are saved.
    do i = 1,im
        if(kpblx(i) > 1) then
            k = kpblx(i)
            if(rbdn(i) >= crb(i)) then
                rbint = 0.
            elseif(rbup(i) <= crb(i)) then
                rbint = 1.
            else
                rbint = (crb(i)-rbdn(i))/(rbup(i)-rbdn(i))
            endif
            hpblx(i) = zl(i,k-1) + rbint*(zl(i,k)-zl(i,k-1))
            if(hpblx(i) < zi(i,kpblx(i))) kpblx(i)=kpblx(i)-1
        else
            hpblx(i) = zl(i,1)
            kpblx(i) = 1
        endif
        hpbl(i) = hpblx(i)
        kpbl(i) = kpblx(i)
        if(kpbl(i) <= 1) pblflg(i)= .FALSE. 
    enddo

!>  - Compute mean tke within pbl - ZNT 12/24/2023 from GFS
      do i = 1, im
          sumx(i) = 0.
          tkemean(i) = 0.
      enddo
      do k = 1, kmpbl
          do i = 1, im
              if(k < kpbl(i)) then
                  dz = zi(i,k+1) - zi(i,k)
                  tkemean(i) = tkemean(i) + tke(i,k) * dz
                  sumx(i) = sumx(i) + dz
              endif
          enddo
      enddo
      do i = 1, im
          if(tkemean(i) > 0. .and. sumx(i) > 0.) then
              tkemean(i) = tkemean(i) / sumx(i)
          endif
      enddo

!>  - Compute wind shear term as a sink term for updraft and downdraft
!!  velocity  - ZNT 12/24/2023 from GFS
      kps = max(kmpbl, kmscu)
      do k = 2, kps
          do i = 1, im
              dz = zi(i,k+1) - zi(i,k)
              tem = (0.5*(u1(i,k-1)-u1(i,k+1))/dz)**2
              tem1 = tem+(0.5*(v1(i,k-1)-v1(i,k+1))/dz)**2
              wush(i,k) = csmf * sqrt(tem1)
          enddo
      enddo

!> ## Compute Monin-Obukhov similarity parameters
!!  - Calculate the Monin-Obukhov nondimensional stability paramter, commonly
!!    referred to as \f$\zeta\f$ using the following equation from Businger et al.(1971) \cite businger_et_al_1971
!!    (eqn 28):
!!    \f[
!!    \zeta = Ri_{sfc}\frac{F_m^2}{F_h} = \frac{z}{L}
!!    \f]
!!    where \f$F_m\f$ and \f$F_h\f$ are surface Monin-Obukhov stability functions calculated in sfc_diff.f and
!!    \f$L\f$ is the Obukhov length.
    do i=1,im
        zol(i) = max(rbsoil(i)*fm(i)*fm(i)/fh(i),rimin)
        if(sfcflg(i)) then
            zol(i) = min(zol(i),-zfmin)
        else
            zol(i) = max(zol(i),zfmin)
        endif
    !>  - Calculate the nondimensional gradients of momentum and temperature (\f$\phi_m\f$ (phim) and \f$\phi_h\f$(phih)) are calculated using
    !! eqns 5 and 6 from Hong and Pan (1996) \cite hong_and_pan_1996 depending on the surface layer stability:
    !!   - For the unstable and neutral conditions:
    !!     \f[
    !!     \phi_m=(1-16\frac{0.1h}{L})^{-1/4}
    !!     \phi_h=(1-16\frac{0.1h}{L})^{-1/2}
    !!     \f]
    !!   - For the stable regime
    !!     \f[
    !!     \phi_m=\phi_t=(1+5\frac{0.1h}{L})
    !!     \f]
        zol1 = zol(i)*sfcfrac*hpbl(i)/zl(i,1)
        if(sfcflg(i)) then
            tem     = 1.0 / (1. - aphi16*zol1)
            phih(i) = sqrt(tem)
            phim(i) = sqrt(phih(i))
        else
            phim(i) = 1. + aphi5*zol1
            phih(i) = phim(i)
        endif
    enddo

!>  - The \f$z/L\f$ (zol) is used as the stability criterion for the PBL.Currently,
!!    strong unstable (convective) PBL for \f$z/L < -0.02\f$ and weakly and moderately
!!    unstable PBL for \f$0>z/L>-0.02\f$
!>  - Compute the velocity scale \f$w_s\f$ (wscale) (eqn 22 of Han et al. 2019). It
!!    is represented by the value scaled at the top of the surface layer:
!!    \f[
!!    w_s=(u_*^3+7\alpha\kappa w_*^3)^{1/3}
!!    \f]
!!    where \f$u_*\f$ (ustar) is the surface friction velocity,\f$\alpha\f$ is the ratio
!!    of the surface layer height to the PBL height (specified as sfcfrac =0.1),
!!    \f$\kappa =0.4\f$ is the von Karman constant, and \f$w_*\f$ is the convective velocity
!!    scale defined as eqn23 of Han et al.(2019):
!!    \f[
!!    w_{*}=[(g/T)\overline{(w'\theta_v^{'})}_0h]^{1/3}
!!    \f]
    do i=1,im
        if(pblflg(i)) then
            if(zol(i) < zolcru) then
                pcnvflg(i) = .true.
            endif
            wst3(i) = gotvx(i,1)*sflux(i)*hpbl(i)
            wstar(i)= wst3(i)**h1
            ust3(i) = ustar(i)**3.
            wscale(i)=(ust3(i)+wfac*vk*wst3(i)*sfcfrac)**h1
            ptem = ustar(i)/aphi5
            wscale(i) = max(wscale(i),ptem)
        endif
    enddo

!>  ## The counter-gradient terms for temperature and humidity are calculated.
!! -  Equation 4 of Hong and Pan (1996) \cite hong_and_pan_1996 and are used to calculate the "scaled virtual temperature excess near the surface" (equation 9 in Hong and Pan (1996) \cite hong_and_pan_1996) for use in the mass-flux algorithm.

    do i = 1,im
        if(pcnvflg(i)) then
            hgamt(i) = heat(i)/wscale(i)
            hgamq(i) = evap(i)/wscale(i)
            vpert(i) = hgamt(i) + hgamq(i)*fv*theta(i,1)
            vpert(i) = max(vpert(i),0.)
            tem = min(cfac*vpert(i),gamcrt)
            thermal(i)= thermal(i) + tem
        endif
    enddo

!  enhance the pbl height by considering the thermal excess
!     (overshoot pbl top)

    do i=1,im
        flg(i)  = .true.
        if(pcnvflg(i)) then
            flg(i)  = .false.
            rbup(i) = rbsoil(i)
        endif
    enddo
    do k = 2, kmpbl
        do i = 1, im
            if( .NOT. flg(i)) then
                rbdn(i) = rbup(i)
                spdk2   = max((u1(i,k)**2+v1(i,k)**2),1.)
                rbup(i) = (thlvx(i,k)-thermal(i))* &
                (g*zl(i,k)/thlvx(i,1))/spdk2
                kpbl(i) = k
                flg(i)  = rbup(i) > crb(i)
            endif
        enddo
    enddo
    do i = 1,im
        if(pcnvflg(i)) then
            k = kpbl(i)
            if(rbdn(i) >= crb(i)) then
                rbint = 0.
            elseif(rbup(i) <= crb(i)) then
                rbint = 1.
            else
                rbint = (crb(i)-rbdn(i))/(rbup(i)-rbdn(i))
            endif
            hpbl(i) = zl(i,k-1) + rbint*(zl(i,k)-zl(i,k-1))
            if(hpbl(i) < zi(i,kpbl(i))) then
                kpbl(i) = kpbl(i) - 1
            endif
            if(kpbl(i) <= 1) then
                pcnvflg(i) = .false.
                pblflg(i) = .false.
            endif
        endif
    enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  look for stratocumulus
!> ## Determine whether stratocumulus layers exist and compute quantities
!!  - Starting at the PBL top and going downward, if the level is less than 2.5 km
!!    and \f$q_l\geq q_{lcr}\f$ then set kcld = k (find the cloud top index in the PBL.
!!    If no cloud water above the threshold is hound, \e scuflg is set to F.
    do i=1,im
        flg(i)  = scuflg(i)
    enddo
    do k = 1, km1
        do i=1,im
            if(flg(i) .AND. zl(i,k) >= zstblmax) then
                lcld(i)=k
                flg(i)=.false.
            endif
        enddo
    enddo
    do i = 1, im
        flg(i)=scuflg(i)
    enddo
    do k = kmscu,1,-1
        do i = 1, im
            if(flg(i) .AND. k <= lcld(i)) then
                if(qlx(i,k) >= qlcr) then
                    kcld(i)=k
                    flg(i)=.false.
                endif
            endif
        enddo
    enddo
    do i = 1, im
        if(scuflg(i) .AND. kcld(i)==km1) scuflg(i)= .FALSE. 
    enddo
!>  - Starting at the PBL top and going downward, if the level is less
!!    than the cloud top, find the level of the minimum radiative heating
!!    rate wihin the cloud. If the level of the minimum is the lowest model
!!    level or the minimum radiative heating rate is positive, then set
!!    scuflg to F.
    do i = 1, im
        flg(i)=scuflg(i)
    enddo
    do k = kmscu,1,-1
        do i = 1, im
            if(flg(i) .AND. k <= kcld(i)) then
                if(qlx(i,k) >= qlcr) then
                    if(radx(i,k) < radmin(i)) then
                        radmin(i)=radx(i,k)
                        krad(i)=k
                    endif
                else
                    flg(i)=.false.
                endif
            endif
        enddo
    enddo
    do i = 1, im
        if(scuflg(i) .AND. krad(i) <= 1) scuflg(i)= .FALSE. 
        if(scuflg(i) .AND. radmin(i)>=0.) scuflg(i)= .FALSE. 
    enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> ## Compute components for mass flux mixing by large thermals
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!>  - If the PBL is convective, the updraft properties are initialized
!!    to be the same as the state variables.
    do k = 1, km
        do i = 1, im
            if(pcnvflg(i)) then
                tcko(i,k) = t1(i,k)
                ucko(i,k) = u1(i,k)
                vcko(i,k) = v1(i,k)
            endif
            if(scuflg(i)) then
                tcdo(i,k) = t1(i,k)
                ucdo(i,k) = u1(i,k)
                vcdo(i,k) = v1(i,k)
            endif
        enddo
    enddo
    do kk = 1, ntrac1
        do k = 1, km
            do i = 1, im
                if(pcnvflg(i)) then
                    qcko(i,k,kk) = q1(i,k,kk)
                else
                endif
                if(scuflg(i)) then
                    qcdo(i,k,kk) = q1(i,k,kk)
                endif
            enddo
        enddo
    enddo

! ZNT 06/28/2023: DEBUG FOR AEROSOL
!      write(*,*) 'aer1', aer1
!      write(*,*) 'aerx1', aerx1

!>  - Call mfpbltq(), which is an EDMF parameterization (Siebesma et al.(2007) \cite Siebesma_2007)
!!    to take into account nonlocal transport by large eddies. For details of the mfpbltq subroutine, step into its documentation ::mfpbltq
    call mfpbltq(im,ix,km,kmpbl,ntcw,ntaw,ntnw,ntrac1,dt2, &
    grav,cp,rd,rv,hvap,fv,eps,epsm1, &       ! ZNT 04/27/2020
    pcnvflg,zl,zm,q1,t1,u1,v1,plyr,pix,thlx,thvx, &
    do_mf_aer,aer1,aerx1, &                  ! ZNT 06/30/2023, add aerosol activation
    gdx,hpbl,kpbl,vpert,  &
    buou,wush,tkemean,xmf, &                 ! ZNT 12/24/2023: wush, tkemean
    tcko,qcko,ucko,vcko,xlamue,bl_upfr, &
    cfl_crit,do_cfl_col,do_cmp_qs,use_hl, &  ! ZNT 09/24/2020; 06/06/2021; 06/12/2023
    do_upd_ovrsht,c_entu, &                  ! ZNT 11/18/2020; 06/21/2021
    do_wush, do_tkemean)                     ! ZNT 12/24/2023

!>  - Call mfscuq(), which is a new mass-flux parameterization for
!!    stratocumulus-top-induced turbulence mixing. For details of the mfscuq subroutine, step into its documentation ::mfscuq
    call mfscuq(im,ix,km,kmscu,ntcw,ntaw,ntnw,ntrac1,dt2, &
    grav,cp,rv,hvap,fv,eps,epsm1, &          ! ZNT 04/27/2020
    scuflg,zl,zm,q1,t1,u1,v1,plyr,pix, &
    thlx,thvx,thlvx,gdx,thetae, &
    krad,mrad,radmin, &
    buod,wush,tkemean,xmfd, &                ! ZNT 12/24/2023: wush, tkemean
    tcdo,qcdo,ucdo,vcdo,xlamde,bl_dnfr, &
    cfl_crit,do_cfl_col,do_cmp_qs,use_hl, &  ! ZNT 09/24/2020; 06/06/2021; 06/12/2023
    cldtime,c_entd,actei, &                  ! ZNT 11/11/2020; 06/21/2021
    do_wush, do_tkemean)                     ! ZNT 12/24/2023

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! ZNT 06/29/2020: Overwriting pbl height
!      do i = 1, im
!         hpbl(i) = hpblx(i)
!         kpbl(i) = kpblx(i)
!      enddo

!> ## Compute Prandtl number \f$P_r\f$ (prn) and exchange coefficient varying with height
    do k = 1, kmpbl
        do i = 1, im
            if(k < kpbl(i)) then
                tem = phih(i)/phim(i)
                ptem = sfcfrac*hpbl(i)
                tem1 = max(zi(i,k+1)-ptem, 0.)
                tem2 = tem1 / (hpbl(i) - ptem)
                if(pcnvflg(i)) then
                    tem = min(tem, pr0)
                    prn(i,k) = tem + (pr0 - tem) * tem2
                else
                    tem = max(tem, pr0)
                    prn(i,k) = tem
                endif
                prn(i,k) = min(prn(i,k),prmax)
                prn(i,k) = max(prn(i,k),prmin)
            
                ckz(i,k) = ck0 + (ck1 - ck0) * tem2
                ckz(i,k) = max(min(ckz(i,k), ck0), ck1)
                chz(i,k) = ch0 + (ch1 - ch0) * tem2
                chz(i,k) = max(min(chz(i,k), ch0), ch1)
            
            endif
        enddo
    enddo

! Above a threshold height (hcrinv), the background vertical diffusivities & mixing length
!    in the inversion layers are set to much smaller values (xkinv & rlmn2)

! Below the threshold height (hcrinv), the background vertical diffusivities & mixing length
!    in the inversion layers are increased with increasing roughness length & vegetation fraction

    do k = 1,km1
        do i=1,im
            if(zi(i,k+1) > hcrinv) then
                tem1 = tvx(i,k+1)-tvx(i,k)
                if(tem1 >= 0.) then
                    xkzo(i,k)  = min(xkzo(i,k), xkinv)
                    xkzmo(i,k) = min(xkzmo(i,k), xkinv)
                    rlmnz(i,k) = min(rlmnz(i,k), rlmn2)
                endif
            else
                tem1 = tvx(i,k+1)-tvx(i,k)
                if(tem1 > 0.) then
                    ptem = xkzo(i,k) * zvfun(i)
                    xkzo(i,k) = min(max(ptem, xkinv), xkzo(i,k))
                    ptem = xkzmo(i,k) * zvfun(i)
                    xkzmo(i,k) = min(max(ptem, xkinv), xkzmo(i,k))
                    ptem = rlmnz(i,k) * zvfun(i)
                    rlmnz(i,k) = min(max(ptem, rlmn2), rlmnz(i,k))
                endif
            endif
        enddo
    enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> ## Compute an asymtotic mixing length

    do k = 1, km1
        do i = 1, im
            zlup = 0.0
            bsum = 0.0
            mlenflg = .true.
            do n = k, km1
                if(mlenflg) then
                    dz = zl(i,n+1) - zl(i,n)
                !             tem1 = 0.5 * (thvx(i,n) + thvx(i,n+1))
                !!            tem1 = 0.5 * (thlvx(i,n) + thlvx(i,n+1))
                    tem3=((u1(i,n+1)-u1(i,n))/dz)**2
                    tem3=tem3+((v1(i,n+1)-v1(i,n))/dz)**2
                    tem3=cs0*sqrt(tem3)*sqrt(tke(i,k))
                    ptem = (gotvx(i,n)*(thvx(i,n+1)-thvx(i,k))+tem3)*dz
                !             ptem = (gotvx(i,n)*(thlvx(i,n+1)-thlvx(i,k)+tem3)*dz
                !             ptem = (gotvx(i,n)*(tem1-thvx(i,k))+tem3)*dz
                !!            ptem = (gotvx(i,n)*(tem1-thlvx(i,k)+tem3)*dz
                    bsum = bsum + ptem
                    zlup = zlup + dz
                    if(bsum >= tke(i,k)) then
                        if(ptem >= 0.) then
                            tem2 = max(ptem, zfmin)
                        else
                            tem2 = min(ptem, -zfmin)
                        endif
                        ptem1 = (bsum - tke(i,k)) / tem2
                        zlup = zlup - ptem1 * dz
                        zlup = max(zlup, 0.)
                        mlenflg = .false.
                    endif
                endif
            enddo
            zldn = 0.0
            bsum = 0.0
            mlenflg = .true.
            do n = k, 1, -1
                if(mlenflg) then
                    if(n == 1) then
                        dz = zl(i,1)
                        tem1 = tsea(i)*(1.+fv*max(q1(i,1,1),qmin))
                    !               tem1 = 0.5 * (tem1 + thvx(i,n))
                    !!              tem1 = 0.5 * (tem1 + thlvx(i,n))
                        tem3 = (u1(i,1)/dz)**2
                        tem3 = tem3+(v1(i,1)/dz)**2
                        tem3 = cs0*sqrt(tem3)*sqrt(tke(i,1))
                    else
                        dz = zl(i,n) - zl(i,n-1)
                        tem1 = thvx(i,n-1)
                    !               tem1 = thlvx(i,n-1)
                    !               tem1 = 0.5 * (thvx(i,n-1) + thvx(i,n))
                    !!              tem1 = 0.5 * (thlvx(i,n-1) + thlvx(i,n))
                        tem3 = ((u1(i,n)-u1(i,n-1))/dz)**2
                        tem3 = tem3+((v1(i,n)-v1(i,n-1))/dz)**2
                        tem3 = cs0*sqrt(tem3)*sqrt(tke(i,k))
                    endif
                    ptem = (gotvx(i,n)*(thvx(i,k)-tem1)+tem3)*dz
                !             ptem = (gotvx(i,n)*(thlvx(i,k)-tem1)+tem3)*dz
                    bsum = bsum + ptem
                    zldn = zldn + dz
                    if(bsum >= tke(i,k)) then
                        if(ptem >= 0.) then
                            tem2 = max(ptem, zfmin)
                        else
                            tem2 = min(ptem, -zfmin)
                        endif
                        ptem1 = (bsum - tke(i,k)) / tem2
                        zldn = zldn - ptem1 * dz
                        zldn = max(zldn, 0.)
                        mlenflg = .false.
                    endif
                endif
            enddo
        
            tem = 0.5 * (zi(i,k+1)-zi(i,k))
            tem1 = min(tem, rlmnz(i,k))
        !>  - Following Bougeault and Lacarrere(1989), the characteristic length
        !! scale (\f$l_2\f$) (eqn 10 in Han et al.(2019) \cite Han_2019) is given by:
        !!\f[
        !!   l_2=min(l_{up},l_{down})
        !!\f]
        !! and dissipation length scale \f$l_d\f$ is given by:
        !!\f[
        !!   l_d=(l_{up}l_{down})^{1/2}
        !!\f]
        !!    where \f$l_{up}\f$ and \f$l_{down}\f$ are the distances that a parcel
        !! having an initial TKE can travel upward and downward before being stopped
        !! by buoyancy effects.
            ptem2 = min(zlup,zldn)
            rlam(i,k) = elmfac * ptem2
            rlam(i,k) = max(rlam(i,k), tem1)
            rlam(i,k) = min(rlam(i,k), rlmx)
        
            ptem2 = sqrt(zlup*zldn)
            ele(i,k) = elefac * ptem2
            ele(i,k) = max(ele(i,k), tem1)
            ele(i,k) = min(ele(i,k), elmx)
        
        enddo
    enddo
!>  - Compute the surface layer length scale (\f$l_1\f$) following
!! Nakanishi (2001) \cite Nakanish_2001 (eqn 9 of Han et al.(2019) \cite Han_2019)
    do k = 1, km1
        do i = 1, im
            tem = vk * zl(i,k)
            if (zol(i) < 0.) then
                ptem = 1. - 100. * zol(i)
                ptem1 = ptem**0.2
                zk = tem * ptem1
            elseif (zol(i) >= 1.) then
                zk = tem / 3.7
            else
                ptem = 1. + 2.7 * zol(i)
                zk = tem / ptem
            endif
            elm(i,k) = zk*rlam(i,k)/(rlam(i,k)+zk)
        
            dz = zi(i,k+1) - zi(i,k)
            tem = max(gdx(i),dz)
            elm(i,k) = min(elm(i,k), tem)
            ele(i,k) = min(ele(i,k), tem)
        
        enddo
    enddo
    do i = 1, im
        elm(i,km) = elm(i,km1)
        ele(i,km) = ele(i,km1)
    enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> ## Compute eddy diffusivities
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    do k = 1, km1
        do i = 1, im
            xkzo(i,k) = 0.5 * (xkzo(i,k) + xkzo(i,k+1))
            xkzmo(i,k) = 0.5 * (xkzmo(i,k) + xkzmo(i,k+1))
        enddo
    enddo
    do k = 1, km1
        do i = 1, im
            tem = 0.5 * (elm(i,k) + elm(i,k+1))
            tem = tem * sqrt(tkeh(i,k))
            if(k < kpbl(i)) then
                if(pcnvflg(i)) then
                    dku(i,k) = ckz(i,k) * tem
                    dkt(i,k) = dku(i,k) / prn(i,k)
                else
                    if(ri(i,k) < 0.) then ! unstable regime
                        dku(i,k) = ckz(i,k) * tem
                        dkt(i,k) = dku(i,k) / prn(i,k)
                    else             ! stable regime
                        dkt(i,k) = chz(i,k) * tem
                        dku(i,k) = dkt(i,k) * prn(i,k)
                    endif
                endif
            else
                if(ri(i,k) < 0.) then ! unstable regime
                    dku(i,k) = ck1 * tem
                    dkt(i,k) = rchck * dku(i,k)
                else             ! stable regime
                    dkt(i,k) = ch1 * tem
                    prnum = 1.0 + 2.1 * ri(i,k)
                    prnum = min(prnum,prmax)
                    dku(i,k) = dkt(i,k) * prnum
                endif
            endif
        
            if(scuflg(i)) then
                if(k >= mrad(i) .AND. k < krad(i)) then
                    tem1 = ckz(i,k) * tem
                    ptem1 = tem1 / prscu
                    dku(i,k) = max(dku(i,k), tem1)
                    dkt(i,k) = max(dkt(i,k), ptem1)
                endif
            endif
        
            dkq(i,k) = prtke * dkt(i,k)
        
            dkt(i,k) = min(dkt(i,k),dkmax)
            dkt(i,k) = max(dkt(i,k),xkzo(i,k))
            dkq(i,k) = min(dkq(i,k),dkmax)
            dkq(i,k) = max(dkq(i,k),xkzo(i,k))
            dku(i,k) = min(dku(i,k),dkmax)
            dku(i,k) = max(dku(i,k),xkzmo(i,k))
        
        enddo
    enddo
!> ## Compute TKE.
!!  - Compute a minimum TKE deduced from background diffusivity for momentum.

    do k = 1, km1
        do i = 1, im
            if(k == 1) then
                tem = ckz(i,1)
                tem1 = 0.5 * xkzmo(i,1)
            else
                tem = 0.5 * (ckz(i,k-1) + ckz(i,k))
                tem1 = 0.5 * (xkzmo(i,k-1) + xkzmo(i,k))
            endif
            ptem = tem1 / (tem * elm(i,k))
            tkmnz(i,k) = ptem * ptem
            tkmnz(i,k) = min(tkmnz(i,k), tkminx)
            tkmnz(i,k) = max(tkmnz(i,k), tkmin)
        enddo
    enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!>  - Compute buoyancy and shear productions of TKE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    do k = 1, km1
        do i = 1, im
            if (k == 1) then
                tem = -dkt(i,1) * bf(i,1)
            !           if(pcnvflg(i)) then
            !             ptem1 = xmf(i,1) * buou(i,1)
            !           else
                ptem1 = 0.
            !           endif
                if(scuflg(i) .AND. mrad(i) == 1) then
                    ptem2 = xmfd(i,1) * buod(i,1)
                else
                    ptem2 = 0.
                endif
                if (do_pos_mf_buop) then       ! ZNT 06/05/2021
                    ptem1 = max(ptem1, 0.0)
                    ptem2 = max(ptem2, 0.0)
                endif
                tem = tem + (ptem1 + ptem2)* mf_buop_fac       ! ZNT 06/05/2021
                buop = 0.5 * (gotvx(i,1) * sflux(i) + tem)
            
                tem1 = dku(i,1) * shr2(i,1)
            
                tem = (u1(i,2)-u1(i,1))*rdzt(i,1)
            !           if(pcnvflg(i)) then
            !             ptem = xmf(i,1) * tem
            !             ptem1 = 0.5 * ptem * (u1(i,2)-ucko(i,2))
            !           else
                ptem1 = 0.
            !           endif
                if(scuflg(i) .AND. mrad(i) == 1) then
                    ptem = ucdo(i,1)+ucdo(i,2)-u1(i,1)-u1(i,2)
                    ptem = 0.5 * tem * xmfd(i,1) * ptem
                else
                    ptem = 0.
                endif
                ptem1 = ptem1 + ptem
            
                tem = (v1(i,2)-v1(i,1))*rdzt(i,1)
            !           if(pcnvflg(i)) then
            !             ptem = xmf(i,1) * tem
            !             ptem2 = 0.5 * ptem * (v1(i,2)-vcko(i,2))
            !           else
                ptem2 = 0.
            !           endif
                if(scuflg(i) .AND. mrad(i) == 1) then
                    ptem = vcdo(i,1)+vcdo(i,2)-v1(i,1)-v1(i,2)
                    ptem = 0.5 * tem * xmfd(i,1) * ptem
                else
                    ptem = 0.
                endif
                ptem2 = ptem2 + ptem
            
            !           tem2 = stress(i)*spd1(i)/zl(i,1)
                tem2 = stress(i)*ustar(i)*phim(i)/(vk*zl(i,1))
                shrp = 0.5 * (tem1 + (ptem1 + ptem2)*mf_shrp_fac + tem2)   ! ZNT 06/05/2021
            else
                tem1 = -dkt(i,k-1) * bf(i,k-1)
                tem2 = -dkt(i,k) * bf(i,k)
                tem  = 0.5 * (tem1 + tem2)
                if(pcnvflg(i) .AND. k <= kpbl(i)) then
                    ptem = 0.5 * (xmf(i,k-1) + xmf(i,k))
                    ptem1 = ptem * buou(i,k)
                else
                    ptem1 = 0.
                endif
                if(scuflg(i)) then
                    if(k >= mrad(i) .AND. k < krad(i)) then
                        ptem0 = 0.5 * (xmfd(i,k-1) + xmfd(i,k))
                        ptem2 = ptem0 * buod(i,k)
                    else
                        ptem2 = 0.
                    endif
                else
                    ptem2 = 0.
                endif
                if (do_pos_mf_buop) then       ! ZNT 06/05/2021
                    ptem1 = max(ptem1, 0.0)
                    ptem2 = max(ptem2, 0.0)
                endif
                buop = tem + (ptem1 + ptem2)* mf_buop_fac       ! ZNT 06/05/2021
            
                tem1 = dku(i,k-1) * shr2(i,k-1)
                tem2 = dku(i,k) * shr2(i,k)
                tem  = 0.5 * (tem1 + tem2)
                tem1 = (u1(i,k+1)-u1(i,k))*rdzt(i,k)
                tem2 = (u1(i,k)-u1(i,k-1))*rdzt(i,k-1)
                if(pcnvflg(i) .AND. k <= kpbl(i)) then
                    ptem = xmf(i,k) * tem1 + xmf(i,k-1) * tem2
                    ptem1 = 0.5 * ptem * (u1(i,k)-ucko(i,k))
                else
                    ptem1 = 0.
                endif
                if(scuflg(i)) then
                    if(k >= mrad(i) .AND. k < krad(i)) then
                        ptem0 = xmfd(i,k) * tem1 + xmfd(i,k-1) * tem2
                        ptem2 = 0.5 * ptem0 * (ucdo(i,k)-u1(i,k))
                    else
                        ptem2 = 0.
                    endif
                else
                    ptem2 = 0.
                endif
                shrp = tem + ptem1 + ptem2
                tem1 = (v1(i,k+1)-v1(i,k))*rdzt(i,k)
                tem2 = (v1(i,k)-v1(i,k-1))*rdzt(i,k-1)
                if(pcnvflg(i) .AND. k <= kpbl(i)) then
                    ptem = xmf(i,k) * tem1 + xmf(i,k-1) * tem2
                    ptem1 = 0.5 * ptem * (v1(i,k)-vcko(i,k))
                else
                    ptem1 = 0.
                endif
                if(scuflg(i)) then
                    if(k >= mrad(i) .AND. k < krad(i)) then
                        ptem0 = xmfd(i,k) * tem1 + xmfd(i,k-1) * tem2
                        ptem2 = 0.5 * ptem0 * (vcdo(i,k)-v1(i,k))
                    else
                        ptem2 = 0.
                    endif
                else
                    ptem2 = 0.
                endif
                shrp = shrp + (ptem1 + ptem2)*mf_shrp_fac  ! ZNT 06/05/2021
            endif
            prod(i,k) = buop + shrp
            de_buop(i,k) = buop   ! ZNT 10/03/2020
            de_shrp(i,k) = shrp   ! ZNT 10/03/2020
        enddo
    enddo

!----------------------------------------------------------------------
!>  - First predict tke due to tke production & dissipation(diss)

    dtn = dt2 / float(ndt)
    do n = 1, ndt
        do k = 1,km1
            do i=1,im
                tem = sqrt(tke(i,k))
                ptem = ce0 / ele(i,k)
            ! if (do_tkediss_imp) then  ! ZNT 06/05/2021: new algorithm
            ! (implicit, not implemented yet)
            ! tke_n = tke_o + (prod - ptem*tke_n^1.5)*dtn
            ! tke_n + ptem*dtn*tke_n^1.5 = (tke_o + prod*dtn)
            ! --> tke_n =
            ! --> diss = ptem*tke_n^1.5 = prod - (tke_n-tke_o)/dtn
            
            
            ! diss(i,k) =
            ! else     ! ZNT 06/05/2021: old algorithm (explicit) with ndt=20 substeps
            ! tke_n = tke_o + (prod - ptem*tke_o^1.5)*dtn
                diss(i,k) = ptem * tke(i,k) * tem
                tem1 = prod(i,k) + tke(i,k) / dtn
                diss(i,k)=max(min(diss(i,k), tem1), 0.)
            ! endif
                de_diss(i,k) = diss(i,k)   ! ZNT 10/03/2020
                tke(i,k) = tke(i,k) + dtn * (prod(i,k)-diss(i,k))
            !          tke(i,k) = max(tke(i,k), tkmin)
                tke(i,k) = max(tke(i,k), tkmnz(i,k))
            enddo
        enddo
    enddo

!>  - Compute updraft & downdraft properties for TKE

    do k = 1, km
        do i = 1, im
            if(pcnvflg(i)) then
                qcko(i,k,ntke) = tke(i,k)
            endif
            if(scuflg(i)) then
                qcdo(i,k,ntke) = tke(i,k)
            endif
        enddo
    enddo
    do k = 2, kmpbl
        do i = 1, im
            if (pcnvflg(i) .AND. k <= kpbl(i)) then
                dz   = zl(i,k) - zl(i,k-1)
                tem  = 0.5 * xlamue(i,k-1) * dz
                factor = 1. + tem
                qcko(i,k,ntke)=((1.-tem)*qcko(i,k-1,ntke)+tem* &
                (tke(i,k)+tke(i,k-1)))/factor
            endif
        enddo
    enddo
    do k = kmscu, 1, -1
        do i = 1, im
            if (scuflg(i) .AND. k < krad(i)) then
                if(k >= mrad(i)) then
                    dz = zl(i,k+1) - zl(i,k)
                    tem  = 0.5 * xlamde(i,k) * dz
                    factor = 1. + tem
                    qcdo(i,k,ntke)=((1.-tem)*qcdo(i,k+1,ntke)+tem* &
                    (tke(i,k)+tke(i,k+1)))/factor
                endif
            endif
        enddo
    enddo

!----------------------------------------------------------------------
!>  - Compute tridiagonal matrix elements for turbulent kinetic energy

    do i=1,im
        ad(i,1) = 1.0
        f1(i,1) = tke(i,1)
    enddo

    if ( .NOT. do_tkevdif) dkq(:,:) = xkzo(:,:)
          
    do k = 1,km1
        do i=1,im
            dtodsd  = dt2/del(i,k)
            dtodsu  = dt2/del(i,k+1)
            dsig    = prsl(i,k)-prsl(i,k+1)
            rdz     = rdzt(i,k)
            tem1    = dsig * dkq(i,k) * rdz
            dsdz2   = tem1 * rdz
            au(i,k) = -dtodsd*dsdz2
            al(i,k) = -dtodsu*dsdz2
            ad(i,k) = ad(i,k)-au(i,k)
            ad(i,k+1)= 1.-al(i,k)
            tem2    = dsig * rdz
        
            if(pcnvflg(i) .AND. k < kpbl(i)) then
                ptem      = 0.5 * tem2 * xmf(i,k) * mf_tkeflx_fac    ! ZNT: 06/06/2021
                ptem1     = dtodsd * ptem
                ptem2     = dtodsu * ptem
                tem       = tke(i,k) + tke(i,k+1)
                ptem      = qcko(i,k,ntke) + qcko(i,k+1,ntke)
                f1(i,k)   = f1(i,k)-(ptem-tem)*ptem1
                f1(i,k+1) = tke(i,k+1)+(ptem-tem)*ptem2
            else
                f1(i,k+1) = tke(i,k+1)
            endif
        
            if(scuflg(i)) then
                if(k >= mrad(i) .AND. k < krad(i)) then
                    ptem      = 0.5 * tem2 * xmfd(i,k) * mf_tkeflx_fac  ! ZNT: 06/06/2021
                    ptem1     = dtodsd * ptem
                    ptem2     = dtodsu * ptem
                    tem       = tke(i,k) + tke(i,k+1)
                    ptem      = qcdo(i,k,ntke) + qcdo(i,k+1,ntke)
                    f1(i,k)   = f1(i,k) + (ptem - tem) * ptem1
                    f1(i,k+1) = f1(i,k+1) - (ptem - tem) * ptem2
                endif
            endif
        
        enddo
    enddo

!>  - Call tridit() to solve tridiagonal problem for TKE

    call tridit(im,km,1,al,ad,au,f1,au,f1)

!>  - Recover the tendency of tke

    do k = 1,km
        do i = 1,im
        !           f1(i,k) = max(f1(i,k), tkmin)
            qtend = (f1(i,k)-q1(i,k,ntke))*rdt
            rtg(i,k,ntke) = rtg(i,k,ntke)+qtend
        enddo
    enddo
! ZNT 05/08/2020: Nullify MF diagnostic fields
    dv_mf = 0.;  dv_mfd = 0.
    du_mf = 0.;  du_mfd = 0.
    tdt_mf = 0.; tdt_mfd = 0.
    qdt_mf = 0.; qdt_mfd = 0.
    qldt_mf = 0.; qldt_mfd = 0.
    qidt_mf = 0.; qidt_mfd = 0.
    qadt_mf = 0.; qadt_mfd = 0.
    qndt_mf = 0.; qndt_mfd = 0.   ! ZNT: 05/12/2023
    qnidt_mf = 0.; qnidt_mfd = 0. ! ZNT: 05/12/2023
    tqpdt_mf = 0.; tqpdt_mfd = 0. ! ZNT: 05/19/2023
    dvsfc_mf = 0.; dvsfc_mfd = 0.
    dusfc_mf = 0.; dusfc_mfd = 0.
    dtsfc_mf = 0.; dtsfc_mfd = 0.
    dqsfc_mf = 0.; dqsfc_mfd = 0.
! ZNT 05/24/2023: Nullify MF subs-detr decomposition fields
    tdt_sub = 0.; tdt_subd = 0.
    qdt_sub = 0.; qdt_subd = 0.
    qldt_sub = 0.; qldt_subd = 0.
    qidt_sub = 0.; qidt_subd = 0.
    qadt_sub = 0.; qadt_subd = 0.
    qndt_sub = 0.; qndt_subd = 0.
    qnidt_sub = 0.; qnidt_subd = 0.
    tdt_dtr = 0.; tdt_dtrd = 0.
    qdt_dtr = 0.; qdt_dtrd = 0.
    qldt_dtr = 0.; qldt_dtrd = 0.
    qidt_dtr = 0.; qidt_dtrd = 0.
    qadt_dtr = 0.; qadt_dtrd = 0.
    qndt_dtr = 0.; qndt_dtrd = 0.
    qnidt_dtr = 0.; qnidt_dtrd = 0.
    tqpdt_sub = 0.; tqpdt_subd = 0. ! ZNT: 05/23/2023
    tqpdt_dtr = 0.; tqpdt_dtrd = 0. ! ZNT: 05/23/2023

!> ## Compute tridiagonal matrix elements for heat and moisture

    do i=1,im
        ad(i,1) = 1.
        f1(i,1) = t1(i,1)   + dtdz1(i) * heat(i)
        f2(i,1) = q1(i,1,1) + dtdz1(i) * evap(i)
    enddo
    if(ntrac1 >= 2) then
        do kk = 2, ntrac1
            is = (kk-1) * km
            do i = 1, im
                f2(i,1+is) = q1(i,1,kk)
            enddo
        enddo
    endif

    do k = 1,km1
        do i = 1,im
            dtodsd  = dt2/del(i,k)
            dtodsu  = dt2/del(i,k+1)
            dsig    = prsl(i,k)-prsl(i,k+1)
            rdz     = rdzt(i,k)
            tem1    = dsig * dkt(i,k) * rdz
            dsdzt   = tem1 * gocp
            dsdz2   = tem1 * rdz
            au(i,k) = -dtodsd*dsdz2
            al(i,k) = -dtodsu*dsdz2
            ad(i,k) = ad(i,k)-au(i,k)
            ad(i,k+1)= 1.-al(i,k)
            tem2    = dsig * rdz
        
            if(pcnvflg(i) .AND. k < kpbl(i)) then
                ptem      = 0.5 * tem2 * xmf(i,k)
                ptem1     = dtodsd * ptem
                ptem2     = dtodsu * ptem
                tem       = t1(i,k) + t1(i,k+1)
                ptem      = tcko(i,k) + tcko(i,k+1)
                f1(i,k)   = f1(i,k)+dtodsd*dsdzt-(ptem-tem)*ptem1
                f1(i,k+1) = t1(i,k+1)-dtodsu*dsdzt+(ptem-tem)*ptem2
            ! ZNT 05/08/2020: Diagnose Ttend by updraft
                tdt_mf(i,k)   = tdt_mf(i,k)  -(ptem-tem)*ptem1
                tdt_mf(i,k+1) = tdt_mf(i,k+1)+(ptem-tem)*ptem2
            ! ZNT 05/24/2023: Add subs-detr decomp (using dry static energy)
                tdt_sub(i,k)   = tdt_sub(i,k)   + &
                (t1(i,k+1)-t1(i,k)+gocp/rdz)*ptem1
                tdt_sub(i,k+1) = tdt_sub(i,k+1) + &
                (t1(i,k+1)-t1(i,k)+gocp/rdz)*ptem2
                tdt_dtr(i,k)   = tdt_dtr(i,k)   + &
                xlamue(i,k)/rdz*(ptem-tem)/2.0*ptem1 - &
                (tcko(i,k)   - t1(i,k))  *2.0*ptem1
                tdt_dtr(i,k+1) = tdt_dtr(i,k+1) + &
                xlamue(i,k)/rdz*(ptem-tem)/2.0*ptem2 + &
                (tcko(i,k+1) - t1(i,k+1))*2.0*ptem2
            ! End ZNT 05/08/2020
                tem       = q1(i,k,1) + q1(i,k+1,1)
                ptem      = qcko(i,k,1) + qcko(i,k+1,1)
                f2(i,k)   = f2(i,k) - (ptem - tem) * ptem1
                f2(i,k+1) = q1(i,k+1,1) + (ptem - tem) * ptem2
            ! ZNT 05/08/2020: Diagnose qtend by updraft
                qdt_mf(i,k)   = qdt_mf(i,k)  -(ptem-tem)*ptem1
                qdt_mf(i,k+1) = qdt_mf(i,k+1)+(ptem-tem)*ptem2
            ! ZNT 05/24/2023: Add subs-detr decomp
                qdt_sub(i,k)   = qdt_sub(i,k)   + &
                (q1(i,k+1,1)-q1(i,k,1))*ptem1
                qdt_sub(i,k+1) = qdt_sub(i,k+1) + &
                (q1(i,k+1,1)-q1(i,k,1))*ptem2
                qdt_dtr(i,k)   = qdt_dtr(i,k)   + &
                xlamue(i,k)/rdz*(ptem-tem)/2.0*ptem1 - &
                (qcko(i,k,1)   - q1(i,k,1))  *2.0*ptem1
                qdt_dtr(i,k+1) = qdt_dtr(i,k+1) + &
                xlamue(i,k)/rdz*(ptem-tem)/2.0*ptem2 + &
                (qcko(i,k+1,1) - q1(i,k+1,1))*2.0*ptem2
            ! End ZNT 05/08/2020
            else
                f1(i,k)   = f1(i,k)+dtodsd*dsdzt
                f1(i,k+1) = t1(i,k+1)-dtodsu*dsdzt
                f2(i,k+1) = q1(i,k+1,1)
            endif
        
            if(scuflg(i)) then
                if(k >= mrad(i) .AND. k < krad(i)) then
                    ptem      = 0.5 * tem2 * xmfd(i,k)
                    ptem1     = dtodsd * ptem
                    ptem2     = dtodsu * ptem
                    ptem      = tcdo(i,k) + tcdo(i,k+1)
                    tem       = t1(i,k) + t1(i,k+1)
                    f1(i,k)   = f1(i,k) + (ptem - tem) * ptem1
                    f1(i,k+1) = f1(i,k+1) - (ptem - tem) * ptem2
                ! ZNT 05/08/2020: Diagnose Ttend by downdraft
                    tdt_mfd(i,k)   = tdt_mfd(i,k)  +(ptem-tem)*ptem1
                    tdt_mfd(i,k+1) = tdt_mfd(i,k+1)-(ptem-tem)*ptem2
                ! ZNT 05/24/2023: Add subs-detr decomp (using dry static energy)
                    tdt_subd(i,k)  = tdt_subd(i,k)   - &
                    (t1(i,k+1)-t1(i,k)+gocp/rdz)*ptem1
                    tdt_subd(i,k+1) = tdt_subd(i,k+1) - &
                    (t1(i,k+1)-t1(i,k)+gocp/rdz)*ptem2
                    tdt_dtrd(i,k)   = tdt_dtrd(i,k)   + &
                    xlamde(i,k)/rdz*(ptem-tem)/2.0*ptem1 + &
                    (tcdo(i,k)   - t1(i,k))  *2.0*ptem1
                    tdt_dtrd(i,k+1) = tdt_dtrd(i,k+1) + &
                    xlamde(i,k)/rdz*(ptem-tem)/2.0*ptem2 - &
                    (tcdo(i,k+1) - t1(i,k+1))*2.0*ptem2
                ! End ZNT 05/08/2020
                    tem       = q1(i,k,1) + q1(i,k+1,1)
                    ptem      = qcdo(i,k,1) + qcdo(i,k+1,1)
                    f2(i,k)   = f2(i,k) + (ptem - tem) * ptem1
                    f2(i,k+1) = f2(i,k+1) - (ptem - tem) * ptem2
                ! ZNT 05/08/2020: Diagnose qtend by downdraft
                    qdt_mfd(i,k)   = qdt_mfd(i,k)  +(ptem-tem)*ptem1
                    qdt_mfd(i,k+1) = qdt_mfd(i,k+1)-(ptem-tem)*ptem2
                ! ZNT 05/24/2023: Add subs-detr decomp
                    qdt_subd(i,k)  = qdt_subd(i,k)   - &
                    (q1(i,k+1,1)-q1(i,k,1))*ptem1
                    qdt_subd(i,k+1) = qdt_subd(i,k+1) - &
                    (q1(i,k+1,1)-q1(i,k,1))*ptem2
                    qdt_dtrd(i,k)   = qdt_dtrd(i,k)   + &
                    xlamde(i,k)/rdz*(ptem-tem)/2.0*ptem1 + &
                    (qcdo(i,k,1)   - q1(i,k,1))  *2.0*ptem1
                    qdt_dtrd(i,k+1) = qdt_dtrd(i,k+1) + &
                    xlamde(i,k)/rdz*(ptem-tem)/2.0*ptem2 - &
                    (qcdo(i,k+1,1) - q1(i,k+1,1))*2.0*ptem2
                ! End ZNT 05/08/2020
                endif
            endif
        enddo
    enddo

    if(ntrac1 >= 2) then
        do kk = 2, ntrac1
            is = (kk-1) * km
            do k = 1, km1
                do i = 1, im
                    if(pcnvflg(i) .AND. k < kpbl(i)) then
                        dtodsd = dt2/del(i,k)
                        dtodsu = dt2/del(i,k+1)
                        dsig  = prsl(i,k)-prsl(i,k+1)
                        tem   = dsig * rdzt(i,k)
                        ptem  = 0.5 * tem * xmf(i,k)
                        ptem1 = dtodsd * ptem
                        ptem2 = dtodsu * ptem
                        tem1  = qcko(i,k,kk) + qcko(i,k+1,kk)
                        tem2  = q1(i,k,kk) + q1(i,k+1,kk)
                    ! ZNT 10/11/2020: Diagnose ql and qa tend by updraft
                    ! ZNT 05/24/2023: Add qn, qni tend; add subs-detr decomp
                        if(kk == ntcw) then
                            qldt_mf(i,k)   = qldt_mf(i,k)  -(tem1-tem2)*ptem1
                            qldt_mf(i,k+1) = qldt_mf(i,k+1)+(tem1-tem2)*ptem2
                            qldt_sub(i,k)   = qldt_sub(i,k)   + &
                            (q1(i,k+1,kk)-q1(i,k,kk))*ptem1
                            qldt_sub(i,k+1) = qldt_sub(i,k+1) + &
                            (q1(i,k+1,kk)-q1(i,k,kk))*ptem2
                            qldt_dtr(i,k)   = qldt_dtr(i,k)   + &
                            xlamue(i,k)/rdzt(i,k)*(tem1-tem2)/2.0*ptem1 - &
                            (qcko(i,k,kk)   - q1(i,k,kk))  *2.0*ptem1
                            qldt_dtr(i,k+1) = qldt_dtr(i,k+1) + &
                            xlamue(i,k)/rdzt(i,k)*(tem1-tem2)/2.0*ptem2 + &
                            (qcko(i,k+1,kk) - q1(i,k+1,kk))*2.0*ptem2
                        elseif(kk == ntiw) then
                            qidt_mf(i,k)   = qidt_mf(i,k)  -(tem1-tem2)*ptem1
                            qidt_mf(i,k+1) = qidt_mf(i,k+1)+(tem1-tem2)*ptem2
                            qidt_sub(i,k)   = qidt_sub(i,k)   + &
                            (q1(i,k+1,kk)-q1(i,k,kk))*ptem1
                            qidt_sub(i,k+1) = qidt_sub(i,k+1) + &
                            (q1(i,k+1,kk)-q1(i,k,kk))*ptem2
                            qidt_dtr(i,k)   = qidt_dtr(i,k)   + &
                            xlamue(i,k)/rdzt(i,k)*(tem1-tem2)/2.0*ptem1 - &
                            (qcko(i,k,kk)   - q1(i,k,kk))  *2.0*ptem1
                            qidt_dtr(i,k+1) = qidt_dtr(i,k+1) + &
                            xlamue(i,k)/rdzt(i,k)*(tem1-tem2)/2.0*ptem2 + &
                            (qcko(i,k+1,kk) - q1(i,k+1,kk))*2.0*ptem2
                        elseif(kk == ntaw) then
                            qadt_mf(i,k)   = qadt_mf(i,k)  -(tem1-tem2)*ptem1
                            qadt_mf(i,k+1) = qadt_mf(i,k+1)+(tem1-tem2)*ptem2
                            qadt_sub(i,k)   = qadt_sub(i,k)   + &
                            (q1(i,k+1,kk)-q1(i,k,kk))*ptem1
                            qadt_sub(i,k+1) = qadt_sub(i,k+1) + &
                            (q1(i,k+1,kk)-q1(i,k,kk))*ptem2
                            qadt_dtr(i,k)   = qadt_dtr(i,k)   + &
                            xlamue(i,k)/rdzt(i,k)*(tem1-tem2)/2.0*ptem1 - &
                            (qcko(i,k,kk)   - q1(i,k,kk))  *2.0*ptem1
                            qadt_dtr(i,k+1) = qadt_dtr(i,k+1) + &
                            xlamue(i,k)/rdzt(i,k)*(tem1-tem2)/2.0*ptem2 + &
                            (qcko(i,k+1,kk) - q1(i,k+1,kk))*2.0*ptem2
                        elseif(kk == ntnw) then           ! ZNT: 05/12/2023
                            qndt_mf(i,k)   = qndt_mf(i,k)  -(tem1-tem2)*ptem1
                            qndt_mf(i,k+1) = qndt_mf(i,k+1)+(tem1-tem2)*ptem2
                            qndt_sub(i,k)   = qndt_sub(i,k)   + &
                            (q1(i,k+1,kk)-q1(i,k,kk))*ptem1
                            qndt_sub(i,k+1) = qndt_sub(i,k+1) + &
                            (q1(i,k+1,kk)-q1(i,k,kk))*ptem2
                            qndt_dtr(i,k)   = qndt_dtr(i,k)   + &
                            xlamue(i,k)/rdzt(i,k)*(tem1-tem2)/2.0*ptem1 - &
                            (qcko(i,k,kk)   - q1(i,k,kk))  *2.0*ptem1
                            qndt_dtr(i,k+1) = qndt_dtr(i,k+1) + &
                            xlamue(i,k)/rdzt(i,k)*(tem1-tem2)/2.0*ptem2 + &
                            (qcko(i,k+1,kk) - q1(i,k+1,kk))*2.0*ptem2
                        elseif(kk == ntniw) then          ! ZNT: 05/12/2023
                            qnidt_mf(i,k)   = qnidt_mf(i,k)  -(tem1-tem2)*ptem1
                            qnidt_mf(i,k+1) = qnidt_mf(i,k+1)+(tem1-tem2)*ptem2
                            qnidt_sub(i,k)   = qnidt_sub(i,k)   + &
                            (q1(i,k+1,kk)-q1(i,k,kk))*ptem1
                            qnidt_sub(i,k+1) = qnidt_sub(i,k+1) + &
                            (q1(i,k+1,kk)-q1(i,k,kk))*ptem2
                            qnidt_dtr(i,k)   = qnidt_dtr(i,k)   + &
                            xlamue(i,k)/rdzt(i,k)*(tem1-tem2)/2.0*ptem1 - &
                            (qcko(i,k,kk)   - q1(i,k,kk))  *2.0*ptem1
                            qnidt_dtr(i,k+1) = qnidt_dtr(i,k+1) + &
                            xlamue(i,k)/rdzt(i,k)*(tem1-tem2)/2.0*ptem2 + &
                            (qcko(i,k+1,kk) - q1(i,k+1,kk))*2.0*ptem2
                        elseif(kk >= ntqpw .AND. kk <= (ntqpw+6)) then ! ZNT: 05/19/2023; 07/05/2023
                            ktqp = kk-ntqpw+1
                            tqpdt_mf(i,k,ktqp)   = &
                            tqpdt_mf(i,k,ktqp)  -(tem1-tem2)*ptem1
                            tqpdt_mf(i,k+1,ktqp) = &
                            tqpdt_mf(i,k+1,ktqp)+(tem1-tem2)*ptem2
                        ! ZNT: 05/23/2023: diagnose subs-detr form of MF tendencies,
                        ! and verify they add up to MF tendency for conserved variables.
                        ! tqpdt_sub, tqpdt_subd, tqpdt_dtr, tqpdt_dtrd
                            tqpdt_sub(i,k,ktqp)   = tqpdt_sub(i,k,ktqp)   + &
                            (q1(i,k+1,kk)-q1(i,k,kk))*ptem1
                            tqpdt_sub(i,k+1,ktqp) = tqpdt_sub(i,k+1,ktqp) + &
                            (q1(i,k+1,kk)-q1(i,k,kk))*ptem2
                            tqpdt_dtr(i,k,ktqp)   = tqpdt_dtr(i,k,ktqp)   + &
                            xlamue(i,k)/rdzt(i,k)*(tem1-tem2)/2.0*ptem1 - &
                            (qcko(i,k,kk)   - q1(i,k,kk))  *2.0*ptem1
                            tqpdt_dtr(i,k+1,ktqp) = tqpdt_dtr(i,k+1,ktqp) + &
                            xlamue(i,k)/rdzt(i,k)*(tem1-tem2)/2.0*ptem2 + &
                            (qcko(i,k+1,kk) - q1(i,k+1,kk))*2.0*ptem2
                        endif
                    ! End ZNT 10/11/2020
                        f2(i,k+is) = f2(i,k+is) - (tem1 - tem2) * ptem1
                        f2(i,k+1+is)= q1(i,k+1,kk) + (tem1 - tem2) * ptem2
                    else
                        f2(i,k+1+is) = q1(i,k+1,kk)
                    endif
                
                    if(scuflg(i)) then
                        if(k >= mrad(i) .AND. k < krad(i)) then
                            dtodsd = dt2/del(i,k)
                            dtodsu = dt2/del(i,k+1)
                            dsig  = prsl(i,k)-prsl(i,k+1)
                            tem   = dsig * rdzt(i,k)
                            ptem  = 0.5 * tem * xmfd(i,k)
                            ptem1 = dtodsd * ptem
                            ptem2 = dtodsu * ptem
                            tem1  = qcdo(i,k,kk) + qcdo(i,k+1,kk)
                            tem2  = q1(i,k,kk) + q1(i,k+1,kk)
                        ! ZNT 10/11/2020: Diagnose ql and qa tend by updraft
                        ! ZNT 05/24/2023: Add qn, qni tend; add subs-detr decomp
                            if(kk == ntcw) then
                                qldt_mfd(i,k)   = qldt_mfd(i,k)  +(tem1-tem2)*ptem1
                                qldt_mfd(i,k+1) = qldt_mfd(i,k+1)-(tem1-tem2)*ptem2
                                qldt_subd(i,k)   = qldt_subd(i,k)   - &
                                (q1(i,k+1,kk)-q1(i,k,kk))*ptem1
                                qldt_subd(i,k+1) = qldt_subd(i,k+1) - &
                                (q1(i,k+1,kk)-q1(i,k,kk))*ptem2
                                qldt_dtrd(i,k)   = qldt_dtrd(i,k)   + &
                                xlamde(i,k)/rdzt(i,k)*(tem1-tem2)/2.0*ptem1 + &
                                (qcdo(i,k,kk)   - q1(i,k,kk))  *2.0*ptem1
                                qldt_dtrd(i,k+1) = qldt_dtrd(i,k+1) + &
                                xlamde(i,k)/rdzt(i,k)*(tem1-tem2)/2.0*ptem2 - &
                                (qcdo(i,k+1,kk) - q1(i,k+1,kk))*2.0*ptem2
                            elseif(kk == ntiw) then
                                qidt_mfd(i,k)   = qidt_mfd(i,k)  +(tem1-tem2)*ptem1
                                qidt_mfd(i,k+1) = qidt_mfd(i,k+1)-(tem1-tem2)*ptem2
                                qidt_subd(i,k)   = qidt_subd(i,k)   - &
                                (q1(i,k+1,kk)-q1(i,k,kk))*ptem1
                                qidt_subd(i,k+1) = qidt_subd(i,k+1) - &
                                (q1(i,k+1,kk)-q1(i,k,kk))*ptem2
                                qidt_dtrd(i,k)   = qidt_dtrd(i,k)   + &
                                xlamde(i,k)/rdzt(i,k)*(tem1-tem2)/2.0*ptem1 + &
                                (qcdo(i,k,kk)   - q1(i,k,kk))  *2.0*ptem1
                                qidt_dtrd(i,k+1) = qidt_dtrd(i,k+1) + &
                                xlamde(i,k)/rdzt(i,k)*(tem1-tem2)/2.0*ptem2 - &
                                (qcdo(i,k+1,kk) - q1(i,k+1,kk))*2.0*ptem2
                            elseif(kk == ntaw) then
                                qadt_mfd(i,k)   = qadt_mfd(i,k)  +(tem1-tem2)*ptem1
                                qadt_mfd(i,k+1) = qadt_mfd(i,k+1)-(tem1-tem2)*ptem2
                                qadt_subd(i,k)   = qadt_subd(i,k)   - &
                                (q1(i,k+1,kk)-q1(i,k,kk))*ptem1
                                qadt_subd(i,k+1) = qadt_subd(i,k+1) - &
                                (q1(i,k+1,kk)-q1(i,k,kk))*ptem2
                                qadt_dtrd(i,k)   = qadt_dtrd(i,k)   + &
                                xlamde(i,k)/rdzt(i,k)*(tem1-tem2)/2.0*ptem1 + &
                                (qcdo(i,k,kk)   - q1(i,k,kk))  *2.0*ptem1
                                qadt_dtrd(i,k+1) = qadt_dtrd(i,k+1) + &
                                xlamde(i,k)/rdzt(i,k)*(tem1-tem2)/2.0*ptem2 - &
                                (qcdo(i,k+1,kk) - q1(i,k+1,kk))*2.0*ptem2
                            elseif(kk == ntnw) then           ! ZNT: 05/12/2023
                                qndt_mfd(i,k)   = qndt_mfd(i,k)  +(tem1-tem2)*ptem1
                                qndt_mfd(i,k+1) = qndt_mfd(i,k+1)-(tem1-tem2)*ptem2
                                qndt_subd(i,k)   = qndt_subd(i,k)   - &
                                (q1(i,k+1,kk)-q1(i,k,kk))*ptem1
                                qndt_subd(i,k+1) = qndt_subd(i,k+1) - &
                                (q1(i,k+1,kk)-q1(i,k,kk))*ptem2
                                qndt_dtrd(i,k)   = qndt_dtrd(i,k)   + &
                                xlamde(i,k)/rdzt(i,k)*(tem1-tem2)/2.0*ptem1 + &
                                (qcdo(i,k,kk)   - q1(i,k,kk))  *2.0*ptem1
                                qndt_dtrd(i,k+1) = qndt_dtrd(i,k+1) + &
                                xlamde(i,k)/rdzt(i,k)*(tem1-tem2)/2.0*ptem2 - &
                                (qcdo(i,k+1,kk) - q1(i,k+1,kk))*2.0*ptem2
                            elseif(kk == ntniw) then          ! ZNT: 05/12/2023
                                qnidt_mfd(i,k)   = qnidt_mfd(i,k)  +(tem1-tem2)*ptem1
                                qnidt_mfd(i,k+1) = qnidt_mfd(i,k+1)-(tem1-tem2)*ptem2
                                qnidt_subd(i,k)   = qnidt_subd(i,k)   - &
                                (q1(i,k+1,kk)-q1(i,k,kk))*ptem1
                                qnidt_subd(i,k+1) = qnidt_subd(i,k+1) - &
                                (q1(i,k+1,kk)-q1(i,k,kk))*ptem2
                                qnidt_dtrd(i,k)   = qnidt_dtrd(i,k)   + &
                                xlamde(i,k)/rdzt(i,k)*(tem1-tem2)/2.0*ptem1 + &
                                (qcdo(i,k,kk)   - q1(i,k,kk))  *2.0*ptem1
                                qnidt_dtrd(i,k+1) = qnidt_dtrd(i,k+1) + &
                                xlamde(i,k)/rdzt(i,k)*(tem1-tem2)/2.0*ptem2 - &
                                (qcdo(i,k+1,kk) - q1(i,k+1,kk))*2.0*ptem2
                            elseif(kk >= ntqpw .AND. kk <= (ntqpw+6)) then ! ZNT: 05/19/2023; 07/05/2023
                                ktqp = kk-ntqpw+1
                                tqpdt_mfd(i,k,ktqp)   = &
                                tqpdt_mfd(i,k,ktqp)  +(tem1-tem2)*ptem1
                                tqpdt_mfd(i,k+1,ktqp) = &
                                tqpdt_mfd(i,k+1,ktqp)-(tem1-tem2)*ptem2
                            ! ZNT: 05/23/2023: diagnose subs-detr form of MF tendencies,
                            ! and verify they add up to MF tendency for conserved variables.
                            ! tqpdt_sub, tqpdt_subd, tqpdt_dtr, tqpdt_dtrd
                                tqpdt_subd(i,k,ktqp)   = tqpdt_subd(i,k,ktqp)   - &
                                (q1(i,k+1,kk)-q1(i,k,kk))*ptem1
                                tqpdt_subd(i,k+1,ktqp) = tqpdt_subd(i,k+1,ktqp) - &
                                (q1(i,k+1,kk)-q1(i,k,kk))*ptem2
                                tqpdt_dtrd(i,k,ktqp)   = tqpdt_dtrd(i,k,ktqp)   + &
                                xlamde(i,k)/rdzt(i,k)*(tem1-tem2)/2.0*ptem1 + &
                                (qcdo(i,k,kk)   - q1(i,k,kk))  *2.0*ptem1
                                tqpdt_dtrd(i,k+1,ktqp) = tqpdt_dtrd(i,k+1,ktqp) + &
                                xlamde(i,k)/rdzt(i,k)*(tem1-tem2)/2.0*ptem2 - &
                                (qcdo(i,k+1,kk) - q1(i,k+1,kk))*2.0*ptem2
                            endif
                        ! End ZNT 10/11/2020
                            f2(i,k+is)  = f2(i,k+is) + (tem1 - tem2) * ptem1
                            f2(i,k+1+is)= f2(i,k+1+is) - (tem1 - tem2) * ptem2
                        endif
                    endif
                
                enddo
            enddo
        enddo
    endif

!> - Call tridin() to solve tridiagonal problem for heat and moisture

    call tridin(im,km,ntrac1,al,ad,au,f1,f2,au,f1,f2)

!> - Recover the tendencies of heat and moisture

    do  k = 1,km
        do i = 1,im
            ttend      = (f1(i,k)-t1(i,k))*rdt
            qtend      = (f2(i,k)-q1(i,k,1))*rdt
            tdt(i,k)   = tdt(i,k)+ttend
            rtg(i,k,1) = rtg(i,k,1)+qtend
            dtsfc(i)   = dtsfc(i)+cont*del(i,k)*ttend
            dqsfc(i)   = dqsfc(i)+conq*del(i,k)*qtend
        ! ZNT 05/08/2020: Diagnose T and q tend by MF
            tdt_mf(i,k)  = tdt_mf(i,k)*rdt
            tdt_mfd(i,k) = tdt_mfd(i,k)*rdt
            qdt_mf(i,k)  = qdt_mf(i,k)*rdt
            qdt_mfd(i,k) = qdt_mfd(i,k)*rdt
            dtsfc_mf(i)  = dtsfc_mf(i)+cont*del(i,k)*tdt_mf(i,k)
            dtsfc_mfd(i) = dtsfc_mfd(i)+cont*del(i,k)*tdt_mfd(i,k)
            dqsfc_mf(i)  = dqsfc_mf(i)+conq*del(i,k)*qdt_mf(i,k)
            dqsfc_mfd(i) = dqsfc_mfd(i)+conq*del(i,k)*qdt_mfd(i,k)
        ! End ZNT 05/08/2020
        ! ZNT 10/11/2020: Diagnose ql and qa tend by MF
            qldt_mf(i,k)  = qldt_mf(i,k)*rdt
            qldt_mfd(i,k) = qldt_mfd(i,k)*rdt
            qidt_mf(i,k)  = qidt_mf(i,k)*rdt
            qidt_mfd(i,k) = qidt_mfd(i,k)*rdt
            qadt_mf(i,k)  = qadt_mf(i,k)*rdt
            qadt_mfd(i,k) = qadt_mfd(i,k)*rdt
            qndt_mf(i,k)  = qndt_mf(i,k)*rdt   ! ZNT: 05/12/2023
            qndt_mfd(i,k) = qndt_mfd(i,k)*rdt
            qnidt_mf(i,k)  = qnidt_mf(i,k)*rdt ! ZNT: 05/12/2023
            qnidt_mfd(i,k) = qnidt_mfd(i,k)*rdt
            tqpdt_mf(i,k,1:7)  = tqpdt_mf(i,k,1:7)*rdt ! ZNT: 05/19/2023; 07/05/2023
            tqpdt_mfd(i,k,1:7) = tqpdt_mfd(i,k,1:7)*rdt
        ! End ZNT 10/11/2020
        ! ZNT: 05/24/2023: Diagnose subs-detr form of MF tendencies
            tdt_sub(i,k)  = tdt_sub(i,k)*rdt
            tdt_subd(i,k) = tdt_subd(i,k)*rdt
            qdt_sub(i,k)  = qdt_sub(i,k)*rdt
            qdt_subd(i,k) = qdt_subd(i,k)*rdt
            qldt_sub(i,k)  = qldt_sub(i,k)*rdt
            qldt_subd(i,k) = qldt_subd(i,k)*rdt
            qidt_sub(i,k)  = qidt_sub(i,k)*rdt
            qidt_subd(i,k) = qidt_subd(i,k)*rdt
            qadt_sub(i,k)  = qadt_sub(i,k)*rdt
            qadt_subd(i,k) = qadt_subd(i,k)*rdt
            qndt_sub(i,k)  = qndt_sub(i,k)*rdt
            qndt_subd(i,k) = qndt_subd(i,k)*rdt
            qnidt_sub(i,k)  = qnidt_sub(i,k)*rdt
            qnidt_subd(i,k) = qnidt_subd(i,k)*rdt
        
            tdt_dtr(i,k)  = tdt_dtr(i,k)*rdt
            tdt_dtrd(i,k) = tdt_dtrd(i,k)*rdt
            qdt_dtr(i,k)  = qdt_dtr(i,k)*rdt
            qdt_dtrd(i,k) = qdt_dtrd(i,k)*rdt
            qldt_dtr(i,k)  = qldt_dtr(i,k)*rdt
            qldt_dtrd(i,k) = qldt_dtrd(i,k)*rdt
            qidt_dtr(i,k)  = qidt_dtr(i,k)*rdt
            qidt_dtrd(i,k) = qidt_dtrd(i,k)*rdt
            qadt_dtr(i,k)  = qadt_dtr(i,k)*rdt
            qadt_dtrd(i,k) = qadt_dtrd(i,k)*rdt
            qndt_dtr(i,k)  = qndt_dtr(i,k)*rdt
            qndt_dtrd(i,k) = qndt_dtrd(i,k)*rdt
            qnidt_dtr(i,k)  = qnidt_dtr(i,k)*rdt
            qnidt_dtrd(i,k) = qnidt_dtrd(i,k)*rdt
        
            tqpdt_sub(i,k,1:7)  = tqpdt_sub(i,k,1:7)*rdt  ! ZNT: 05/23/2023; 07/05/2023
            tqpdt_subd(i,k,1:7) = tqpdt_subd(i,k,1:7)*rdt ! ZNT: 05/23/2023; 07/05/2023
            tqpdt_dtr(i,k,1:7)  = tqpdt_dtr(i,k,1:7)*rdt  ! ZNT: 05/23/2023; 07/05/2023
            tqpdt_dtrd(i,k,1:7) = tqpdt_dtrd(i,k,1:7)*rdt ! ZNT: 05/23/2023; 07/05/2023
        ! End ZNT 05/24/2023
        enddo
    enddo

    if(ntrac1 >= 2) then
        do kk = 2, ntrac1
            is = (kk-1) * km
            do k = 1, km
                do i = 1, im
                    qtend = (f2(i,k+is)-q1(i,k,kk))*rdt
                    rtg(i,k,kk) = rtg(i,k,kk)+qtend
                enddo
            enddo
        enddo
    endif

!> ## Add TKE dissipative heating to temperature tendency

    if(dspheat) then
        do k = 1,km1
            do i = 1,im
            !         tem = min(diss(i,k), dspmax)
            !         ttend = tem / cp
                ttend = diss(i,k) / cp
                tdt(i,k) = tdt(i,k) + dspfac * ttend
            enddo
        enddo
    endif

!> ## Compute tridiagonal matrix elements for momentum

    do i=1,im
    ! Implicit
        ad(i,1) = 1.0 + dtdz1(i) * stress(i) / spd1(i)
        f1(i,1) = u1(i,1)
        f2(i,1) = v1(i,1)
    ! Explicit - ZNT 05/07/2020, test only
    !         ad(i,1) = 1.
    !         f1(i,1) = u1(i,1)-dtdz1(i)*stress(i)/spd1(i)*u1(i,1)
    !         f2(i,1) = v1(i,1)-dtdz1(i)*stress(i)/spd1(i)*v1(i,1)
    enddo

    do k = 1,km1
        do i=1,im
            dtodsd  = dt2/del(i,k)
            dtodsu  = dt2/del(i,k+1)
            dsig    = prsl(i,k)-prsl(i,k+1)
            rdz     = rdzt(i,k)
            tem1    = dsig * dku(i,k) * rdz
            dsdz2   = tem1*rdz
            au(i,k) = -dtodsd*dsdz2
            al(i,k) = -dtodsu*dsdz2
            ad(i,k) = ad(i,k)-au(i,k)
            ad(i,k+1)= 1.-al(i,k)
            tem2    = dsig * rdz
        
            if(pcnvflg(i) .AND. k < kpbl(i)) then
                ptem      = 0.5 * tem2 * xmf(i,k)
                ptem1     = dtodsd * ptem
                ptem2     = dtodsu * ptem
                tem       = u1(i,k) + u1(i,k+1)
                ptem      = ucko(i,k) + ucko(i,k+1)
                f1(i,k)   = f1(i,k) - (ptem - tem) * ptem1
                f1(i,k+1) = u1(i,k+1) + (ptem - tem) * ptem2
            ! ZNT 05/08/2020: Diagnose Utend by updraft
                du_mf(i,k)   = du_mf(i,k)  -(ptem-tem)*ptem1
                du_mf(i,k+1) = du_mf(i,k+1)+(ptem-tem)*ptem2
            ! End ZNT 05/08/2020
                tem       = v1(i,k) + v1(i,k+1)
                ptem      = vcko(i,k) + vcko(i,k+1)
                f2(i,k)   = f2(i,k) - (ptem - tem) * ptem1
                f2(i,k+1) = v1(i,k+1) + (ptem - tem) * ptem2
            ! ZNT 05/08/2020: Diagnose Vtend by updraft
                dv_mf(i,k)   = dv_mf(i,k)  -(ptem-tem)*ptem1
                dv_mf(i,k+1) = dv_mf(i,k+1)+(ptem-tem)*ptem2
            ! End ZNT 05/08/2020
            else
                f1(i,k+1) = u1(i,k+1)
                f2(i,k+1) = v1(i,k+1)
            endif
        
            if(scuflg(i)) then
                if(k >= mrad(i) .AND. k < krad(i)) then
                    ptem      = 0.5 * tem2 * xmfd(i,k)
                    ptem1     = dtodsd * ptem
                    ptem2     = dtodsu * ptem
                    tem       = u1(i,k) + u1(i,k+1)
                    ptem      = ucdo(i,k) + ucdo(i,k+1)
                    f1(i,k)   = f1(i,k) + (ptem - tem) *ptem1
                    f1(i,k+1) = f1(i,k+1) - (ptem - tem) *ptem2
                ! ZNT 05/08/2020: Diagnose Utend by downdraft
                    du_mfd(i,k)   = du_mfd(i,k)  +(ptem-tem)*ptem1
                    du_mfd(i,k+1) = du_mfd(i,k+1)-(ptem-tem)*ptem2
                ! End ZNT 05/08/2020
                    tem       = v1(i,k) + v1(i,k+1)
                    ptem      = vcdo(i,k) + vcdo(i,k+1)
                    f2(i,k)   = f2(i,k) + (ptem - tem) * ptem1
                    f2(i,k+1) = f2(i,k+1) - (ptem - tem) * ptem2
                ! ZNT 05/08/2020: Diagnose Vtend by downdraft
                    dv_mfd(i,k)   = dv_mfd(i,k)  +(ptem-tem)*ptem1
                    dv_mfd(i,k+1) = dv_mfd(i,k+1)-(ptem-tem)*ptem2
                ! End ZNT 05/08/2020
                endif
            endif
        
        enddo
    enddo

!> - Call tridi2() to solve tridiagonal problem for momentum

    call tridi2(im,km,al,ad,au,f1,f2,au,f1,f2)

!> - Recover the tendencies of momentum

    do k = 1,km
        do i = 1,im
            utend = (f1(i,k)-u1(i,k))*rdt
            vtend = (f2(i,k)-v1(i,k))*rdt
            du(i,k)  = du(i,k)+utend
            dv(i,k)  = dv(i,k)+vtend
            dusfc(i) = dusfc(i)+conw*del(i,k)*utend
            dvsfc(i) = dvsfc(i)+conw*del(i,k)*vtend
        ! ZNT 05/08/2020: Diagnose U and V tend by MF
            du_mf(i,k)  = du_mf(i,k)*rdt
            du_mfd(i,k) = du_mfd(i,k)*rdt
            dv_mf(i,k)  = dv_mf(i,k)*rdt
            dv_mfd(i,k) = dv_mfd(i,k)*rdt
            dusfc_mf(i)  = dusfc_mf(i)+conw*del(i,k)*du_mf(i,k)
            dusfc_mfd(i) = dusfc_mfd(i)+conw*del(i,k)*du_mfd(i,k)
            dvsfc_mf(i)  = dvsfc_mf(i)+conw*del(i,k)*dv_mf(i,k)
            dvsfc_mfd(i) = dvsfc_mfd(i)+conw*del(i,k)*dv_mfd(i,k)
        ! End ZNT 05/08/2020
        enddo
    enddo

! ZNT 06/19/2020
    do k = 1,km
        do i = 1,im
            if(pcnvflg(i)) then
                tcko_o(i,k)  = tcko(i,k)
                qvcko_o(i,k) = qcko(i,k,1)
                qlcko_o(i,k) = qcko(i,k,ntcw)
                qicko_o(i,k) = qcko(i,k,ntiw)
                qacko_o(i,k) = qcko(i,k,ntaw)   ! ZNT: 05/23/2023
                qncko_o(i,k) = qcko(i,k,ntnw)   ! ZNT: 05/12/2023
                qnicko_o(i,k) = qcko(i,k,ntniw) ! ZNT: 05/12/2023
                tqpcko_o(i,k,1:7) = qcko(i,k,ntqpw:(ntqpw+6)) ! ZNT: 05/19/2023; 07/05/2023
                buou_o(i,k)  = buou(i,k)
            else
                tcko_o(i,k)  = 0.0
                qvcko_o(i,k) = 0.0
                qlcko_o(i,k) = 0.0
                qicko_o(i,k) = 0.0
                qacko_o(i,k) = 0.0   ! ZNT: 05/23/2023
                qncko_o(i,k) = 0.0   ! ZNT: 05/12/2023
                qnicko_o(i,k) = 0.0  ! ZNT: 05/12/2023
                tqpcko_o(i,k,1:7) = 0.0 ! ZNT: 05/19/2023; 07/05/2023
                buou_o(i,k)  = 0.0
            endif
            if(scuflg(i)) then
                tcdo_o(i,k)  = tcdo(i,k)
                qvcdo_o(i,k) = qcdo(i,k,1)
                qlcdo_o(i,k) = qcdo(i,k,ntcw)
                qicdo_o(i,k) = qcdo(i,k,ntiw)
                qacdo_o(i,k) = qcdo(i,k,ntaw)   ! ZNT: 05/23/2023
                qncdo_o(i,k) = qcdo(i,k,ntnw)   ! ZNT: 05/12/2023
                qnicdo_o(i,k) = qcdo(i,k,ntniw) ! ZNT: 05/12/2023
                tqpcdo_o(i,k,1:7) = qcdo(i,k,ntqpw:(ntqpw+6)) ! ZNT: 05/19/2023; 07/05/2023
                buod_o(i,k)  = buod(i,k)
            else
                tcdo_o(i,k)  = 0.0
                qvcdo_o(i,k) = 0.0
                qlcdo_o(i,k) = 0.0
                qicdo_o(i,k) = 0.0
                qacdo_o(i,k) = 0.0   ! ZNT: 05/23/2023
                qncdo_o(i,k) = 0.0   ! ZNT: 05/12/2023
                qnicdo_o(i,k) = 0.0  ! ZNT: 05/12/2023
                tqpcdo_o(i,k,1:7) = 0.0 ! ZNT: 05/19/2023; 07/05/2023
                buod_o(i,k)  = 0.0
            endif
        enddo
    enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> ## Save PBL height for diagnostic purpose

    do i = 1, im
        hpbl(i) = hpblx(i)
        kpbl(i) = kpblx(i)
    enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    return
    end subroutine satmedmfvdifq_run

    end module satmedmfvdifq_mod
