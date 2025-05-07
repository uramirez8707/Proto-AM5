!>\file mfpbltq.f
!! This file contains the subroutine that calculates mass flux and
!! updraft parcel properties for thermals driven by surface heating
!! for use in the TKE-EDMF PBL scheme (updated version).

!>\ingroup satmedmfvdifq
!! This subroutine computes mass flux and updraft parcel properties for
!! thermals driven by surface heating.
!!\section mfpbltq_gen GFS mfpblt General Algorithm

    subroutine mfpbltq(im,ix,km,kmpbl, & 
    ntcw,ntaw,ntnw,ntrac1,delt, &                    ! ZNT 06/30/2023: ntaw, ntnw
    grav,cp,rd,rv,hvap,fv,eps,epsm1, &               ! ZNT 04/27/2020
    cnvflg,zl,zm,q1,t1,u1,v1,plyr,pix,thlx,thvx, &
    do_mf_aer,aer1,aerx1, &                          ! ZNT 06/30/2023, add aerosol activation
    gdx,hpbl,kpbl,vpert,buo,wush,tkemean,xmf, &      ! ZNT 12/24/2023: wush, tkemean
    tcko,qcko,ucko,vcko,xlamue,a1, &
    cfl_crit, do_cfl_col, do_cmp_qs, use_hl, &       ! ZNT 09/24/2020; 06/06/2021; 06/12/2023
    do_upd_ovrsht, c_entu, &                         ! ZNT 11/18/2020; 06/21/2021
    do_wush, do_tkemean)                             ! ZNT 12/24/2023

    use sat_vapor_pres_mod, only :lookup_es,compute_qs ! ZNT 04/27/2020
!      use funcphys , only : fpvs
    use aer_ccn_act_k_mod, only : aer_ccn_act_k
    use fms_mod, only : error_mesg, FATAL, NOTE      ! ZNT 06/30/2023

    implicit none

    integer :: im, ix, km, kmpbl, ntcw, ntaw, ntnw, ntrac1   ! ZNT 06/30/2023: ntaw,ntnw
    integer :: kpbl(im)
    logical :: cnvflg(im)
    real    :: grav,cp,rd,rv,hvap,fv,eps,epsm1  ! ZNT 04/27/2020
    real    :: cfl_crit                      ! ZNT 09/24/2020
    logical :: do_cfl_col                    ! ZNT 06/06/2021
    logical :: do_cmp_qs                     ! ZNT 10/09/2020
    integer :: use_hl                        ! ZNT 06/12/2023: 0=no; 1=Tu; 2=bu; 3=both
    logical :: do_upd_ovrsht                 ! ZNT 11/18/2020
    real    :: c_entu                        ! ZNT 06/21/2021
    logical :: do_wush                       ! ZNT 12/24/2023
    logical :: do_tkemean                    ! ZNT 12/24/2023

    real    :: delt
    real    :: q1(ix,km,ntrac1), &
               t1(ix,km),  u1(ix,km), v1(ix,km), &
               plyr(im,km),pix(im,km),thlx(im,km), &
               thvx(im,km),zl(im,km), zm(im,km), &
               gdx(im),    hpbl(im),  vpert(im), &
               buo(im,km), wush(im,km),  &
               tkemean(im),xmf(im,km), &
               tcko(im,km),qcko(im,km,ntrac1), &
               ucko(im,km),vcko(im,km), &
               xlamue(im,km-1)

    integer :: do_mf_aer                      ! ZNT: 06/30/2023, option for aerosol activation in updrafts
                                              !      0 = skip; 1 = 'do_online_aerosol = .false.'
                                              !      2 = 'do_online_aerosol = .true.'
    real    :: aer1(im,km,4), aerx1(im,km,4)  ! ZNT: 06/30/2023, add aer1 and aerx1
    

!  local variables and arrays

    integer :: i, j, k, n, ndc
    integer :: kpblx(im), kpbly(im)

    real    :: dt2,     dz,      ce0,     cm, &
               tkcrt,   cmxfac,               &  ! ZNT 12/24/2023
               factor,  gocp, &
               g,       b1,      f1, &
               bb1,     bb2, &
               alp,     vprtmax, a1,      pgcon, &
               qmin,    qlmin,   xmmx,    rbint, &
               tem,     tem1,    tem2, &
               ptem,    ptem1,   ptem2

    real    :: elocp,   el2orc,  qs,      es, &
               tlu,     gamma,   qlu, &
               thup,    dq

    real    :: rbdn(im), rbup(im), hpblx(im), &
               xlamuem(im,km-1)
    real    :: delz(im), xlamax(im), ce0t(im)   !  ZNT 12/24/2023 - ce0t

    real    :: wu2(im,km), thlu(im,km), thvu(im,km), &
               qtx(im,km), qtu(im,km), thvu0(im,km), &
               hlx(im,km), hlu(im,km), &       ! ZNT 06/12/2023
               qnx(im,km), qnu(im,km)          ! ZNT 06/30/2023

    real    :: xlamavg(im),   sigma(im), &
               scaldfunc(im), sumx(im)

    real    :: cfl_fac(im)                     ! ZNT 06/06/2021

    logical :: totflg, flg(im)
    logical :: adjzpbl

    ! ZNT 06/30/2023: variables for aerosol activation
    logical :: below_lcl
    real    :: totalmass(4)
    integer :: tym, ier
    real    :: drop, tvu, rhou, wu, qn_asc, qn_ent, qn_act, emass
    character(len=256)   :: ermesg


!  physical parameters
!     parameter(ce0=0.4,cm=1.0)
    parameter(cm=1.0)
    parameter(tkcrt=2.,cmxfac=5.)        !  ZNT 12/24/2023
    parameter(qmin=1.e-8,qlmin=1.e-12)
    parameter(alp=1.5,vprtmax=3.0,pgcon=0.55)
    parameter(b1=0.5,f1=0.15)

! ZNT 06/21/2021
    ce0 = c_entu
! ZNT 04/27/2020
    g=grav; gocp=g/cp; elocp=hvap/cp; el2orc=hvap*hvap/(rv*cp)
! ZNT 06/30/2023
    ier = 0; ermesg = ' '; tym = size(totalmass,1)


!************************************************************************
!!
    adjzpbl = .true.   ! ZNT 07/01/2020: false = kill updraft at diagnosed PBL top
!       (default) true will instead adjust PBL top
    totflg = .true.
    do i=1,im
        totflg = totflg .and. (.not. cnvflg(i))
    enddo
    if(totflg) return
!!

    dt2 = delt

    do k = 1, km
        do i=1,im
            if (cnvflg(i)) then
                buo(i,k) = 0.
                wu2(i,k) = 0.
                qtx(i,k) = q1(i,k,1) + q1(i,k,ntcw)
                hlx(i,k) = thlx(i,k) / pix(i,k) + zl(i,k)*gocp   ! ZNT 06/12/2023: hlx is h_l/c_p
                qnx(i,k) = q1(i,k,ntnw)   ! ZNT 06/30/2023
            endif
        enddo
    enddo

!> - Compute thermal excess

    do i=1,im
        if(cnvflg(i)) then
            ptem = alp * vpert(i)
            ptem = min(ptem, vprtmax)
            thlu(i,1)= thlx(i,1) + ptem
            hlu(i,1) = thlu(i,1) / pix(i,1) + zl(i,1)*gocp     ! ZNT 06/12/2023
            qtu(i,1) = qtx(i,1)
            qnu(i,1) = qnx(i,1)           ! ZNT 06/30/2023
            buo(i,1) = g * ptem / thvx(i,1)
            thvu(i,1)= thvx(i,1) + ptem   ! ZNT
            thvu0(i,1)=thvu(i,1)          ! ZNT
        endif
    enddo

!> - Compute entrainment rate

    ! ZNT 12/24/2023 - compute TKE impact on entrainment rate
    do i=1,im
       if(cnvflg(i)) then
          ce0t(i) = ce0
          if(do_tkemean .and. tkemean(i) > tkcrt) then
            tem = sqrt(tkemean(i)/tkcrt)
            tem1 = min(tem, cmxfac)
            tem2 = tem1 * ce0
            ce0t(i) = max(ce0t(i), tem2)
          endif
       endif
    enddo

    do i=1,im
        if(cnvflg(i)) then
            k = kpbl(i) / 2
            k = max(k, 1)
            delz(i) = zl(i,k+1) - zl(i,k)
            xlamax(i) = ce0t(i) / delz(i)     ! ZNT 12/24/2023 - ce0t
        endif
    enddo

    do k = 1, kmpbl
        do i=1,im
            if(cnvflg(i)) then
                if(k < kpbl(i)) then
                    ptem = 1./(zm(i,k)+delz(i))
                    tem = max((hpbl(i)-zm(i,k)+delz(i)) ,delz(i))
                    ptem1 = 1./tem
                    xlamue(i,k) = ce0t(i) * (ptem+ptem1)  ! ZNT 12/24/2023 - ce0t
                else
                    xlamue(i,k) = xlamax(i)
                endif
            
                xlamuem(i,k) = cm * xlamue(i,k)
            endif
        enddo
    enddo

!> - Compute buoyancy for updraft air parcel

    do k = 2, kmpbl
        do i=1,im
            if(cnvflg(i)) then
                dz   = zl(i,k) - zl(i,k-1)
                tem  = 0.5 * xlamue(i,k-1) * dz
                factor = 1. + tem
            
                thlu(i,k) = ((1.-tem)*thlu(i,k-1)+tem* &
                            (thlx(i,k-1)+thlx(i,k)))/factor
                hlu(i,k) = ((1.-tem)*hlu(i,k-1)+tem* &
                            (hlx(i,k-1)+hlx(i,k)))/factor
                qtu(i,k) = ((1.-tem)*qtu(i,k-1)+tem* &
                            (qtx(i,k-1)+qtx(i,k)))/factor
            ! ZNT 06/23/2020: Diagnose thvu just with mixing
                thvu0(i,k) = ((1.-tem)*thvu0(i,k-1)+tem* &
                            (thvx(i,k-1)+thvx(i,k)))/factor
            
                if(use_hl==2 .OR. use_hl==3) then      ! ZNT 06/12/2023: use h_l instead of th_l for better conservation
                    tlu = hlu(i,k) - zl(i,k)*gocp
                else
                    tlu = thlu(i,k) / pix(i,k)
                endif
                if (do_cmp_qs) then  ! ZNT 10/09/2020: note es is not used
                    call compute_qs(tlu, plyr(i,k)*100.0, qs)
                    qs = max(qmin,qs)
                else
                    call lookup_es(tlu, es)
                    es = 0.01*es             ! ZNT 04/27/2020
                ! es = 0.01 * fpvs(tlu)      ! fpvs in pa
                    qs = max(qmin, eps * es / (plyr(i,k)+epsm1*es))
                endif
                dq = qtu(i,k) - qs
            
                if (dq > 0.) then
                !! ZNT 06/23/2020: Saturation adjustment:
                !!     Tl = T-Lv/cp*(qt-qs(T)) = T-Lv/cp*(qt-qs(Tl)-dqs/dT*(T-Tl))
                !! ==> set G=Lv/cp*dqs/dT ==> (T-Tl)(1+G) = Lv/cp*(qt-qs(Tl))
                !! ==> T-Tl = (Lv/cp*dq)/(1+G); ql =(T-Tl)/(Lv/cp)=dq/(1+G)
                !!     qs(T)=qs(Tl)+dqs/dT*Lv/cp*dq/(1+G) = qs(Tl)+dq*G/(1+G)
                !!     With CC, dqs/dT~(Lv/Rv/T^2)*qs => G~(Lv^2/cp/Rv/T^2)*qs
                !!     Or change to more strict formula using lookup_des?
                    gamma = el2orc * qs / (tlu**2)
                    qlu = dq / (1. + gamma)
                !! ZNT 06/23/2020: Why isn't qs updated below? qtu is not conserved!
                !!     If qs=qs(old)+qlu*gamma ==> qs+qlu=qs(old)+dq=qtu, conserved!
                !!     Start output analysis for k==8 (just below cloud top of Rf01).
                !              qtu(i,k) = qs + qlu    ! Changed ZNT 06/23/2020
                    qs = qtu(i,k) - qlu
                    tem1 = 1. + fv * qs - qlu
                    if(use_hl==2 .OR. use_hl==3) then      ! ZNT 06/12/2023: use h_l instead of th_l for better conservation
                        thup = pix(i,k) * (tlu + elocp * qlu)
                    else
                        thup = thlu(i,k) + pix(i,k) * elocp * qlu
                    endif
                    thvu(i,k) = thup * tem1
                else
                    tem1 = 1. + fv * qtu(i,k)
                    if(use_hl==2 .OR. use_hl==3) then      ! ZNT 06/12/2023: use h_l instead of th_l for better conservation
                        thvu(i,k) = pix(i,k)*tlu * tem1
                    else
                        thvu(i,k) = thlu(i,k) * tem1
                    endif
                endif
                buo(i,k) = g * (thvu(i,k) / thvx(i,k) - 1.)
            
            endif
        enddo
    enddo

!> - Compute updraft velocity square(wu2, eqn 13 in
!! Han et al.(2019) \cite Han_2019)

!     tem = 1.-2.*f1
!     bb1 = 2. * b1 / tem
!     bb2 = 2. / tem
!  from Soares et al. (2004,QJRMS)
!     bb1 = 2.
!     bb2 = 4.

!  from Bretherton et al. (2004, MWR)
!     bb1 = 4.
!     bb2 = 2.

!  from our tuning
    bb1 = 2.0
    bb2 = 4.0

    do i = 1, im
        if(cnvflg(i)) then
            dz   = zm(i,1)
            tem  = 0.5*bb1*xlamue(i,1)*dz
            tem1 = bb2 * buo(i,1) * dz
            ptem1 = 1. + tem
            wu2(i,1) = tem1 / ptem1
        endif
    enddo
    do k = 2, kmpbl
        do i = 1, im
            if(cnvflg(i)) then
                dz    = zm(i,k) - zm(i,k-1)
                tem  = 0.25*bb1*(xlamue(i,k)+xlamue(i,k-1))*dz
                if (do_wush) then  ! ZNT 12/24/2023
                   tem1 = max(wu2(i,k-1), 0.)
                   tem1 = bb2 * buo(i,k) - wush(i,k) * sqrt(tem1)
                   tem2 = tem1 * dz
                   ptem = (1. - tem) * wu2(i,k-1)
                   ptem1 = 1. + tem
                   wu2(i,k) = (ptem + tem2) / ptem1
                else
                   tem1 = bb2 * buo(i,k) * dz
                   ptem = (1. - tem) * wu2(i,k-1)
                   ptem1 = 1. + tem
                   wu2(i,k) = (ptem + tem1) / ptem1
                endif ! do_wush

                if( .NOT. adjzpbl) then    ! ZNT: 07/01/2020
                    wu2(i,k) = max(wu2(i,k), 0.0)  ! wu2 must be non-negative
                    if(k >= kpbl(i)) then  ! artificially stop updraft
                        wu2(i,k) = 0.0       ! when it overshoots zpbl
                    endif
                endif
            endif
        enddo
    enddo

!> - Update pbl height as the height where updraft velocity vanishes

    if(adjzpbl) then
        if(do_upd_ovrsht) then    ! ZNT 11/18/2020: use codes before mid-2020
            do i=1,im
                flg(i)  = .true.
                if(cnvflg(i)) then
                    flg(i)  = .false.
                    rbup(i) = wu2(i,1)
                endif
            enddo
            do k = 2, kmpbl
                do i = 1, im
                    if( .NOT. flg(i)) then
                        rbdn(i) = rbup(i)
                        rbup(i) = wu2(i,k)
                        kpbl(i)= k
                        flg(i)  = rbup(i).le.0.
                    endif
                enddo
            enddo
            do i = 1,im
                if(cnvflg(i)) then
                    k = kpbl(i)
                    if(rbdn(i) <= 0.) then
                        rbint = 0.
                    elseif(rbup(i) >= 0.) then
                        rbint = 1.
                    else
                        rbint = rbdn(i)/(rbdn(i)-rbup(i))
                    endif
                    hpbl(i) = zm(i,k-1) + rbint*(zm(i,k)-zm(i,k-1))
                endif
            enddo

        else   ! (do_upd_ovrsht)
            do i=1,im
                flg(i)  = .true.
                kpblx(i) = 1
                kpbly(i) = kpbl(i)
                if(cnvflg(i)) then
                    flg(i)  = .false.
                    rbup(i) = wu2(i,1)
                endif
            enddo
            do k = 2, kmpbl
                do i = 1, im
                    if( .NOT. flg(i)) then
                        rbdn(i) = rbup(i)
                        rbup(i) = wu2(i,k)
                        kpblx(i)= k
                        flg(i)  = rbup(i).le.0.
                    endif
                enddo
            enddo
            do i = 1,im
                if(cnvflg(i)) then
                    k = kpblx(i)
                    if(rbdn(i) <= 0.) then
                        rbint = 0.
                    elseif(rbup(i) >= 0.) then
                        rbint = 1.
                    else
                        rbint = rbdn(i)/(rbdn(i)-rbup(i))
                    endif
                    hpblx(i) = zm(i,k-1) + rbint*(zm(i,k)-zm(i,k-1))
                endif
            enddo
        
            do i = 1,im
                if(cnvflg(i)) then
                    if(kpblx(i) < kpbl(i)) then
                        kpbl(i) = kpblx(i)
                        hpbl(i) = hpblx(i)
                    endif
                    if(kpbl(i) <= 1) cnvflg(i)= .FALSE. 
                endif
            enddo
        endif ! (do_upd_ovrsht)
    
    endif ! (adjzpbl)

!> - Update entrainment rate

    do i=1,im
        if(cnvflg(i)) then
            k = kpbl(i) / 2
            k = max(k, 1)
            delz(i) = zl(i,k+1) - zl(i,k)
            xlamax(i) = ce0t(i) / delz(i)  ! ZNT 12/24/2023 - ce0t
        endif
    enddo

    do k = 1, kmpbl
        do i=1,im
            if(do_upd_ovrsht) then       ! ZNT 11/18/2020: use codes before mid-2020
                if(cnvflg(i)) then
                    if(k < kpbl(i)) then
                        ptem = 1./(zm(i,k)+delz(i))
                        tem = max((hpbl(i)-zm(i,k)+delz(i)) ,delz(i))
                        ptem1 = 1./tem
                        xlamue(i,k) = ce0t(i) * (ptem+ptem1)   ! ZNT 12/24/2023 - ce0t
                    else
                        xlamue(i,k) = xlamax(i)
                    endif
                
                    xlamuem(i,k) = cm * xlamue(i,k)
                endif

            else ! (do_upd_ovrsht)
                if(cnvflg(i) .AND. kpblx(i) < kpbly(i)) then
                !         if(cnvflg(i)) then
                    if(k < kpbl(i)) then
                        ptem = 1./(zm(i,k)+delz(i))
                        tem = max((hpbl(i)-zm(i,k)+delz(i)) ,delz(i))
                        ptem1 = 1./tem
                        xlamue(i,k) = ce0 * (ptem+ptem1)
                    else
                        xlamue(i,k) = xlamax(i)
                    endif
                
                    xlamuem(i,k) = cm * xlamue(i,k)
                endif
            endif ! (do_upd_ovrsht)
        enddo
    enddo

!> - Compute entrainment rate averaged over the whole pbl

    do i = 1, im
        xlamavg(i) = 0.
        sumx(i) = 0.
    enddo
    do k = 1, kmpbl
        do i = 1, im
            if (cnvflg(i) .AND. k < kpbl(i)) then
                dz = zl(i,k+1) - zl(i,k)
                xlamavg(i) = xlamavg(i) + xlamue(i,k) * dz
                sumx(i) = sumx(i) + dz
            endif
        enddo
    enddo
    do i = 1, im
        if(cnvflg(i)) then
            xlamavg(i) = xlamavg(i) / sumx(i)
        endif
    enddo

!> - Updraft mass flux as a function of updraft velocity profile

    do k = 1, kmpbl
        do i = 1, im
            if (cnvflg(i) .AND. k < kpbl(i)) then
                xmf(i,k) = a1 * sqrt(wu2(i,k))
            endif
        enddo
    enddo

!> - Compute updraft fraction as a function of mean entrainment rate
!!(Grell and Freitas (2014) \cite grell_and_freitas_2014

    do i = 1, im
        if(cnvflg(i)) then
            tem = 0.2 / xlamavg(i)
            tem1 = 3.14 * tem * tem
            sigma(i) = tem1 / (gdx(i) * gdx(i))
            sigma(i) = max(sigma(i), 0.001)
            sigma(i) = min(sigma(i), 0.999)
        endif
    enddo

!> - Compute scale-aware function based on
!! Arakawa and Wu (2013) \cite arakawa_and_wu_2013

    do i = 1, im
        if(cnvflg(i)) then
            if (sigma(i) > a1) then
                scaldfunc(i) = (1.-sigma(i)) * (1.-sigma(i))
                scaldfunc(i) = max(min(scaldfunc(i), 1.0), 0.)
            else
                scaldfunc(i) = 1.0
            endif
        endif
    enddo

!> - Final scale-aware updraft mass flux

    cfl_fac(:) = 1.0     ! ZNT 06/06/2021: column CFL limiter
    do k = 1, kmpbl
        do i = 1, im
            if (cnvflg(i) .AND. k < kpbl(i)) then
                xmf(i,k) = scaldfunc(i) * xmf(i,k)
                dz   = zl(i,k+1) - zl(i,k)
                xmmx = dz / dt2 * cfl_crit      ! ZNT 09/24/2020: max MF allowed by CFL
                if (do_cfl_col) then            ! ZNT 06/06/2021: column CFL limiter
                    if (xmmx < xmf(i,k)) then
                        cfl_fac(i) = min(cfl_fac(i), xmmx/xmf(i,k))
                    endif
                else                            ! ZNT 06/06/2021: local CFL limiter
                    xmf(i,k) = min(xmf(i,k),xmmx)
                endif
            endif
        enddo
    enddo
    if (do_cfl_col) then ! ZNT 06/06/2021: column CFL limiter
        do k = 1, kmpbl
            do i = 1, im
                xmf(i,k) = xmf(i,k)*cfl_fac(i)
            enddo
        enddo
    endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> - Compute updraft property using updated entranment rate

    do i=1,im
        if(cnvflg(i)) then
            thlu(i,1)= thlx(i,1)
            hlu(i,1) = hlx(i,1)    ! ZNT 06/12/2023: use h_l instead of th_l for better conservation
        endif
    enddo

!     do i=1,im
!       if(cnvflg(i)) then
!         ptem1 = max(qcko(i,1,ntcw), 0.)
!         tlu = thlu(i,1) / pix(i,1)
!         tcko(i,1) = tlu +  elocp * ptem1
!       endif
!     enddo

    below_lcl = .true.   ! ZNT 06/30/2023: loop from the surface

    do k = 2, kmpbl
        do i=1,im
            if(cnvflg(i) .AND. k <= kpbl(i)) then
                dz   = zl(i,k) - zl(i,k-1)
                tem  = 0.5 * xlamue(i,k-1) * dz
                factor = 1. + tem
            
                thlu(i,k) = ((1.-tem)*thlu(i,k-1)+tem* &
                            (thlx(i,k-1)+thlx(i,k)))/factor
                hlu(i,k)  = ((1.-tem)*hlu(i,k-1)+tem* &
                            (hlx(i,k-1)+hlx(i,k)))/factor
                qtu(i,k)  = ((1.-tem)*qtu(i,k-1)+tem* &
                            (qtx(i,k-1)+qtx(i,k)))/factor
                            
                qn_asc = (1.-tem)*qnu(i,k-1)/factor        ! ZNT 06/30/2023: qn_u from level below
                qn_ent = tem*(qnx(i,k-1)+qnx(i,k))/factor  ! ZNT 06/30/2023: qn_u from passive entrainment
                qnu(i,k) = qn_asc + qn_ent                 ! ZNT 06/30/2023: passive transport of qn_u
                               
                if(use_hl==1 .OR. use_hl==3) then      ! ZNT 06/12/2023: use h_l instead of th_l for better conservation
                    tlu = hlu(i,k) - zl(i,k)*gocp
                else
                    tlu = thlu(i,k) / pix(i,k)
                endif
                if (do_cmp_qs) then  ! ZNT 10/09/2020: note es is not used
                    call compute_qs(tlu, plyr(i,k)*100.0, qs)
                    qs = max(qmin,qs)
                else
                    call lookup_es(tlu, es)
                    es = 0.01*es             ! ZNT 04/27/2020
                ! es = 0.01 * fpvs(tlu)      ! fpvs in pa
                    qs = max(qmin, eps * es / (plyr(i,k)+epsm1*es))
                endif
                dq = qtu(i,k) - qs
            
                if (dq > 0.) then
                    gamma = el2orc * qs / (tlu**2)
                    qlu = dq / (1. + gamma)
                !              qtu(i,k) = qs + qlu
                    qs = qtu(i,k) - qlu
                    qcko(i,k,1) = qs
                    qcko(i,k,ntcw) = qlu
                    qcko(i,k,ntaw) = 1.      ! ZNT 06/30/2023: adjust cloud amount
                    tcko(i,k) = tlu + elocp * qlu
                else
                    qnu(i,k) = 0.            ! ZNT 06/30/2023: set qn_u to zero when updraft is not saturated
                    qcko(i,k,1) = qtu(i,k)
                    qcko(i,k,ntcw) = 0.
                    qcko(i,k,ntaw) = 0.      ! ZNT 06/30/2023: adjust cloud amount
                    tcko(i,k) = tlu
                endif
                
                ! ZNT 06/30/2023: Add aerosol activation (and qa output) here
                !     NOTE: the current form is simply copied from uw_conv and all variables need to be modified...
                !           if do_mf_aer == .false., the following would be skipped and qn_u would be passively  
                !           transported, except that it would be set to zero when updraft is not saturated.
                !     Need to declare variables: tvu, rhou, wu, totalmass, tym, drop, qn_act, ier, ermesg
                !     Algorithm follows 'do_new_qnact' and 'do_2nd_act' in conv_plumes_k.F90 (UW scheme)
                
                if ((do_mf_aer .ge. 1) .and. (dq > 0.)) then
                   tvu = tcko(i,k) * (1.+fv*qcko(i,k,1)-qcko(i,k,ntcw))   ! Tvu = Tu*(1+0.608*q - ql)
                   rhou= plyr(i,k) * 100.0 /(rd * tvu)                    ! rhou = p/(Rd*Tvu) in kg/m3; plyr is in hPa
                   wu = sqrt(max(wu2(i,k), 0.0))
                   
                   ! Q: (1) grid staggering; (2) possible double counting if environment is cloudy
                   
                   if (below_lcl) then               
                      emass = 1.0            ! Cloud base: all updraft air mass is treated as newly entrained
                      below_lcl = .false.
                   else
                      emass= 2.0*tem/factor  ! Lateral: entrained ambient air mass (kg) per unit kg of updraft air mass (~ entr*dz)
                   endif  ! (below_lcl) 
                   
                   totalmass(1:4)=(aerx1(i,k-1,1:4) + aerx1(i,k,1:4))/2.0 * & 
                                   emass *     &     ! entrained aerosol mass per unit kg of updraft air mass
                                   rhou*1.0e-3       ! convert mixing ratio kg/kg to g/cm3 (with 1.0e-3 factor)
                   
                   ! ZNT 06/30/2023: DEBUG
                   ! write(*,*) 'lev,T,p,rho,w,tym,totalmass',  k,tcko(i,k),plyr(i,k)*100.0,rhou,wu,tym,totalmass
                   
                   if ((do_mf_aer .ge. 2) .and. SUM(totalmass(:))/=0. .and. (wu.gt.0.)) then
                      call aer_ccn_act_k(tcko(i,k), plyr(i,k)*100.0, wu, totalmass, &   ! Note: plyr is in hPa
                                         tym, drop, ier, ermesg)         
                      if (ier /= 0) then
                         call error_mesg ('aer_ccn_act_k in mfpbltq (Han EDMF scheme)', ermesg, FATAL)
                      endif
                      qn_act = drop*1.0e6 /rhou
              
                   else
                      drop = 0.
                      qn_act = 0.0
                   endif
                   
                   qnu(i,k) = qnu(i,k) + qn_act
                   ! Question: Should qn_ent be neglected, as activation is computed again for entrained air mass? 
                   ! qnu(i,k) = qn_asc + qn_act
                   ! qnu(i,k) = qn_asc + max(qn_ent,qn_act)
                   
                   ! ZNT 06/30/2023: DEBUG
                   ! write(*,*) 'do_mf_aer,drop,qnu,qn_asc,qn_ent,qn_act,qne', &
                   !     do_mf_aer,drop,qnu(i,k),qn_asc,qn_ent,qn_act,q1(i,k,ntnw)
                         
                endif   ! (do_mf_aer .and. (dq > 0.))
                
                qcko(i,k,ntnw) = qnu(i,k)              ! copy the final qn_u value to the tracer output array
                
            endif ! (cnvflg(i) .AND. k <= kpbl(i))
        enddo
    enddo

    do k = 2, kmpbl
        do i = 1, im
            if (cnvflg(i) .AND. k <= kpbl(i)) then
                dz   = zl(i,k) - zl(i,k-1)
                tem  = 0.5 * xlamuem(i,k-1) * dz
                factor = 1. + tem
                ptem = tem + pgcon
                ptem1= tem - pgcon
                ucko(i,k) = ((1.-tem)*ucko(i,k-1)+ptem*u1(i,k) &
                             +ptem1*u1(i,k-1))/factor
                vcko(i,k) = ((1.-tem)*vcko(i,k-1)+ptem*v1(i,k) &
                             +ptem1*v1(i,k-1))/factor
            endif
        enddo
    enddo

    ! ZNT 06/30/2023: Change to full loop over all tracers and skip qc, qa, qn
    do n = 2, ntrac1
       if (n .eq. ntcw) cycle
       if (n .eq. ntaw) cycle    ! ZNT 06/30/2023: skip the passive transport for cloud amount
       if (n .eq. ntnw) cycle    ! ZNT 06/30/2023: skip the passive transport for qn
       do k = 2, kmpbl
          do i = 1, im
             if (cnvflg(i) .AND. k <= kpbl(i)) then
                 dz   = zl(i,k) - zl(i,k-1)
                 tem  = 0.5 * xlamue(i,k-1) * dz
                 factor = 1. + tem
                    
                 qcko(i,k,n) = ((1.-tem)*qcko(i,k-1,n)+tem* &
                               (q1(i,k,n)+q1(i,k-1,n)))/factor
             endif
           enddo
       enddo
    enddo
    
    return
    end subroutine mfpbltq
