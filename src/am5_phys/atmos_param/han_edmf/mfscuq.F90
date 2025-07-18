!>\file mfscuq.f
!! This file contains the mass flux and downdraft parcel preperties
!! parameterization for stratocumulus-top-driven turbulence (updated version).

!>\ingroup satmedmfvdifq
!! This subroutine computes mass flux and downdraft parcel properties
!! for stratocumulus-top-driven turbulence.
!! \section mfscuq GFS mfscu General Algorithm

    subroutine mfscuq(im,ix,km,kmscu, & 
    ntcw,ntaw,ntnw,ntrac1,delt, &                   ! ZNT 06/30/2023: ntaw, ntnw
    grav,cp,rv,hvap,fv,eps,epsm1, &                 ! ZNT 04/27/2020
    cnvflg,zl,zm,q1,t1,u1,v1,plyr,pix, &
    thlx,thvx,thlvx,gdx,thetae, &
    krad,mrad,radmin,buo,wush,tkemean,xmfd, &       ! ZNT 12/24/2023: wush, tkemean
    tcdo,qcdo,ucdo,vcdo,xlamde,a1, &
    cfl_crit, do_cfl_col, do_cmp_qs, use_hl, &      ! ZNT 09/24/2020; 06/06/2021; 06/12/2023
    cldtime, c_entd, actei, &                       ! ZNT 11/11/2020; 06/21/2021
    do_wush, do_tkemean)                            ! ZNT 12/24/2023

    use sat_vapor_pres_mod, only :lookup_es,compute_qs ! ZNT 04/27/2020
!      use funcphys , only : fpvs


    implicit none

    integer :: im, ix,  km, kmscu, ntcw, ntaw, ntnw, ntrac1   ! ZNT 06/30/2023: ntaw,ntnw
    integer :: krad(im), mrad(im)

    logical :: cnvflg(im)
    real    :: delt
    real    :: grav,cp,rv,hvap,fv,eps,epsm1  ! ZNT 04/27/2020
    real    :: cfl_crit                      ! ZNT 09/24/2020
    logical :: do_cfl_col                    ! ZNT 06/06/2021
    logical :: do_cmp_qs                     ! ZNT 10/09/2020
    integer :: use_hl                        ! ZNT 06/12/2023: 0=no; 1=Tu; 2=bu; 3=both
    real    :: cldtime                       ! ZNT 11/11/2020
    real    :: c_entd                        ! ZNT 06/21/2021
    real    :: actei
    logical :: do_wush                       ! ZNT 12/24/2023
    logical :: do_tkemean                    ! ZNT 12/24/2023


    real    :: q1(ix,km,ntrac1),t1(ix,km),   &
               u1(ix,km),      v1(ix,km),    &
               plyr(im,km),    pix(im,km),   &
               thlx(im,km), &
               thvx(im,km),    thlvx(im,km), &
               gdx(im), &
               zl(im,km),      zm(im,km),    &
               thetae(im,km),  radmin(im),   &
               buo(im,km),     wush(im,km),  &
               tkemean(im),    xmfd(im,km),  &
               tcdo(im,km),    qcdo(im,km,ntrac1), &
               ucdo(im,km),    vcdo(im,km),  &
               xlamde(im,km-1)

!  local variables and arrays


    integer ::   i,j,indx, k, n, kk, ndc
    integer ::   krad1(im)

    real   :: dt2,     dz,      ce0,     cm,   &
              tkcrt,   cmxfac,                 &  ! ZNT 12/24/2023
              gocp,    factor,  g,       tau,  &
              b1,      f1,      bb1,     bb2,  &
              a1,      a2, &
              cteit,   pgcon, &
              qmin,    qlmin, &
              xmmx,    tem,     tem1,    tem2, &
              ptem,    ptem1,   ptem2

    real   :: elocp,   el2orc,  qs,      es,   &
              tld,     gamma,   qld,     thdn, &
              thvd,    dq

    real   :: wd2(im,km), thld(im,km), &
              qtx(im,km), qtd(im,km),  &
              hlx(im,km), hld(im,km),  &       !  ZNT 06/12/2023
              qnx(im,km), qnd(im,km),  &       !  ZNT 06/30/2023
              thlvd(im),  hrad(im),    &
              xlamdem(im,km-1), ra1(im)

    real   :: delz(im), xlamax(im), ce0t(im)   !  ZNT 12/24/2023 - ce0t

    real   :: xlamavg(im),   sigma(im), &
              scaldfunc(im), sumx(im)

    real   :: cfl_fac(im)                     ! ZNT 06/06/2021

    logical :: totflg, flg(im)

!  physical parameters
!     parameter(ce0=0.4,cm=1.0,pgcon=0.55)
    parameter(cm=1.0,pgcon=0.55)
    parameter(tkcrt=2.,cmxfac=5.)        !  ZNT 12/24/2023
    parameter(qmin=1.e-8,qlmin=1.e-12)
    parameter(b1=0.45,f1=0.15)
    parameter(a2=0.5)
!     parameter(cldtime=500.)
!     parameter(actei = 0.7)
!     parameter(actei = 0.23)

! ZNT 06/21/2021
    ce0 = c_entd
! ZNT 04/27/2020
    g=grav; gocp=g/cp; elocp=hvap/cp; el2orc=hvap*hvap/(rv*cp)

!************************************************************************
!!
    totflg = .true.
    do i=1,im
        totflg = totflg .and. (.not. cnvflg(i))
    enddo
    if(totflg) return
!!
    dt2 = delt

    do k = 1, km
        do i=1,im
            if(cnvflg(i)) then
                buo(i,k) = 0.
                wd2(i,k) = 0.
                qtx(i,k) = q1(i,k,1) + q1(i,k,ntcw)
                hlx(i,k) = thlx(i,k) / pix(i,k) + zl(i,k)*gocp   ! ZNT 06/12/2023: hlx is h_l/c_p
                qnx(i,k) = q1(i,k,ntnw)   ! ZNT 08/03/2023
            endif
        enddo
    enddo

    do i = 1, im
        if(cnvflg(i)) then
            hrad(i) = zm(i,krad(i))
            krad1(i) = krad(i)-1
        endif
    enddo

    do i = 1, im
        if(cnvflg(i)) then
            k    = krad(i)
            tem  = zm(i,k+1)-zm(i,k)
            tem1 = cldtime*radmin(i)/tem
            tem1 = max(tem1, -3.0)
            thld(i,k)= thlx(i,k) + tem1
            hld(i,k) = thld(i,k) / pix(i,k) + zl(i,k)*gocp     ! ZNT 06/12/2023
            qtd(i,k) = qtx(i,k)
            qnd(i,k) = qnx(i,k)           ! ZNT 06/30/2023
            thlvd(i) = thlvx(i,k) + tem1
            buo(i,k) = - g * tem1 / thvx(i,k)
        endif
    enddo

!> - Specify downdraft fraction

    do i=1,im
        if(cnvflg(i)) then
            ra1(i) = a1
        endif
    enddo

!> - If the condition for cloud-top instability is met,
!! increase downdraft fraction

    do i = 1, im
        if(cnvflg(i)) then
            k = krad(i)
            tem = thetae(i,k) - thetae(i,k+1)
            tem1 = qtx(i,k) - qtx(i,k+1)
            if (tem > 0. .AND. tem1 > 0.) then
                cteit= cp*tem/(hvap*tem1)
                if(cteit > actei) then
                    ra1(i) = a2
                endif
            endif
        endif
    enddo

!> - First-guess level of downdraft extension (mrad)

    do i = 1, im
        flg(i) = cnvflg(i)
        mrad(i) = krad(i)
    enddo
    do k = kmscu,1,-1
        do i = 1, im
            if(flg(i) .AND. k < krad(i)) then
                if(thlvd(i) <= thlvx(i,k)) then
                    mrad(i) = k
                else
                    flg(i)=.false.
                endif
            endif
        enddo
    enddo
    do i=1,im
        if (cnvflg(i)) then
            kk = krad(i)-mrad(i)
            if(kk < 1) cnvflg(i)= .FALSE. 
        endif
    enddo
!!
    totflg = .true.
    do i=1,im
        totflg = totflg .and. (.not. cnvflg(i))
    enddo
    if(totflg) return
!!

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
            k = mrad(i) + (krad(i)-mrad(i)) / 2
            k = max(k, mrad(i))
            delz(i) = zl(i,k+1) - zl(i,k)
            xlamax(i) = ce0t(i) / delz(i)     ! ZNT 12/24/2023 - ce0t
        endif
    enddo

    do k = 1, kmscu
        do i=1,im
            if(cnvflg(i)) then
                if(k >= mrad(i) .AND. k < krad(i)) then
                    if(mrad(i) == 1) then
                        ptem = 1./(zm(i,k)+delz(i))
                    else
                        ptem = 1./(zm(i,k)-zm(i,mrad(i)-1)+delz(i))
                    endif
                    tem = max((hrad(i)-zm(i,k)+delz(i)) ,delz(i))
                    ptem1 = 1./tem
                    xlamde(i,k) = ce0t(i) * (ptem+ptem1)  ! ZNT 12/24/2023 - ce0t
                else
                    xlamde(i,k) = xlamax(i)
                endif
            
                xlamdem(i,k) = cm * xlamde(i,k)
            endif
        enddo
    enddo

!> - Compute buoyancy for downdraft air parcel

    do k = kmscu,1,-1
        do i=1,im
            if(cnvflg(i) .AND. k < krad(i)) then
                dz = zl(i,k+1) - zl(i,k)
                tem  = 0.5 * xlamde(i,k) * dz
                factor = 1. + tem
            
                thld(i,k) = ((1.-tem)*thld(i,k+1)+tem* &
                            (thlx(i,k)+thlx(i,k+1)))/factor
                hld(i,k)  = ((1.-tem)*hld(i,k+1)+tem* &
                            (hlx(i,k)+hlx(i,k+1)))/factor
                qtd(i,k)  = ((1.-tem)*qtd(i,k+1)+tem* &
                            (qtx(i,k)+qtx(i,k+1)))/factor
            
                if(use_hl==2 .OR. use_hl==3) then      ! ZNT 06/12/2023: use h_l instead of th_l for better conservation
                    tld = hld(i,k) - zl(i,k)*gocp
                else
                    tld = thld(i,k) / pix(i,k)
                endif
                if (do_cmp_qs) then  ! ZNT 10/09/2020: note es is not used
                    call compute_qs(tld, plyr(i,k)*100.0, qs)
                    qs = max(qmin,qs)
                else
                    call lookup_es(tld, es)
                    es = 0.01*es             ! ZNT 04/27/2020
                ! es = 0.01 * fpvs(tlu)      ! fpvs in pa
                    qs = max(qmin, eps * es / (plyr(i,k)+epsm1*es))
                endif
                dq = qtd(i,k) - qs
            
                if (dq > 0.) then
                    gamma = el2orc * qs / (tld**2)
                    qld = dq / (1. + gamma)
                !              qtd(i,k) = qs + qld
                    qs = qtd(i,k) - qld
                    tem1 = 1. + fv * qs - qld
                    if(use_hl==2 .OR. use_hl==3) then      ! ZNT 06/12/2023: use h_l instead of th_l for better conservation
                        thdn = pix(i,k) * (tld + elocp * qld)
                    else
                        thdn = thld(i,k) + pix(i,k) * elocp * qld
                    endif
                    thvd = thdn * tem1
                else
                    tem1 = 1. + fv * qtd(i,k)
                    if(use_hl==2 .OR. use_hl==3) then      ! ZNT 06/12/2023: use h_l instead of th_l for better conservation
                        thvd = pix(i,k)*tld * tem1
                    else
                        thvd = thld(i,k) * tem1
                    endif
                endif
                buo(i,k) = g * (1. - thvd / thvx(i,k))
            
            endif
        enddo
    enddo

!> - Compute downdraft velocity square(wd2)

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
            k = krad1(i)
            dz = zm(i,k+1) - zm(i,k)
        !         tem = 0.25*bb1*(xlamde(i,k)+xlamde(i,k+1))*dz
            tem = 0.5*bb1*xlamde(i,k)*dz
            tem1 = bb2 * buo(i,k+1) * dz
            ptem1 = 1. + tem
            wd2(i,k) = tem1 / ptem1
        endif
    enddo
    do k = kmscu,1,-1
        do i = 1, im
            if(cnvflg(i) .AND. k < krad1(i)) then
                dz    = zm(i,k+1) - zm(i,k)
                tem  = 0.25*bb1*(xlamde(i,k)+xlamde(i,k+1))*dz
                if (do_wush) then  ! ZNT 12/24/2023
                   tem1 = max(wd2(i,k+1), 0.)
                   tem1 = bb2*buo(i,k+1) - wush(i,k+1)*sqrt(tem1)
                   tem2 = tem1 * dz
                   ptem = (1. - tem) * wd2(i,k+1)
                   ptem1 = 1. + tem
                   wd2(i,k) = (ptem + tem2) / ptem1
                else
                   tem1 = bb2 * buo(i,k+1) * dz
                   ptem = (1. - tem) * wd2(i,k+1)
                   ptem1 = 1. + tem
                   wd2(i,k) = (ptem + tem1) / ptem1
                endif ! do_wush
            endif
        enddo
    enddo

    do i = 1, im
        flg(i) = cnvflg(i)
        if(flg(i)) mrad(i) = krad(i)
    enddo
    do k = kmscu,1,-1
        do i = 1, im
            if(flg(i) .AND. k < krad(i)) then
                if(wd2(i,k) > 0.) then
                    mrad(i) = k
                else
                    flg(i)=.false.
                endif
            endif
        enddo
    enddo

    do i=1,im
        if (cnvflg(i)) then
            kk = krad(i)-mrad(i)
            if(kk < 1) cnvflg(i)= .FALSE. 
        endif
    enddo
!!
    totflg = .true.
    do i=1,im
        totflg = totflg .and. (.not. cnvflg(i))
    enddo
    if(totflg) return
!!

!> - Update entrainment rate

    do i=1,im
        if(cnvflg(i)) then
            k = mrad(i) + (krad(i)-mrad(i)) / 2
            k = max(k, mrad(i))
            delz(i) = zl(i,k+1) - zl(i,k)
            xlamax(i) = ce0t(i) / delz(i)  ! ZNT 12/24/2023 - ce0t
        endif
    enddo

    do k = 1, kmscu
        do i=1,im
            if(cnvflg(i)) then
                if(k >= mrad(i) .AND. k < krad(i)) then
                    if(mrad(i) == 1) then
                        ptem = 1./(zm(i,k)+delz(i))
                    else
                        ptem = 1./(zm(i,k)-zm(i,mrad(i)-1)+delz(i))
                    endif
                    tem = max((hrad(i)-zm(i,k)+delz(i)) ,delz(i))
                    ptem1 = 1./tem
                    xlamde(i,k) = ce0t(i) * (ptem+ptem1)   ! ZNT 12/24/2023 - ce0t
                else
                    xlamde(i,k) = xlamax(i)
                endif
            
                xlamdem(i,k) = cm * xlamde(i,k)
            endif
        enddo
    enddo

!> - Compute entrainment rate averaged over the whole downdraft layers

    do i = 1, im
        xlamavg(i) = 0.
        sumx(i) = 0.
    enddo
    do k = kmscu, 1, -1
        do i = 1, im
            if(cnvflg(i) .AND. &
            (k >= mrad(i) .AND. k < krad(i))) then
                dz = zl(i,k+1) - zl(i,k)
                xlamavg(i) = xlamavg(i) + xlamde(i,k) * dz
                sumx(i) = sumx(i) + dz
            endif
        enddo
    enddo
    do i = 1, im
        if(cnvflg(i)) then
            xlamavg(i) = xlamavg(i) / sumx(i)
        endif
    enddo

!> - Compute downdraft mass flux

    do k = kmscu, 1, -1
        do i = 1, im
            if(cnvflg(i) .AND. &
            (k >= mrad(i) .AND. k < krad(i))) then
                xmfd(i,k) = ra1(i) * sqrt(wd2(i,k))
            endif
        enddo
    enddo

!> - Compute downdraft fraction as a function of mean entrainment rate
!! (Grell and Freitas(2014) \cite grell_and_freitas_2014

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
            if (sigma(i) > ra1(i)) then
                scaldfunc(i) = (1.-sigma(i)) * (1.-sigma(i))
                scaldfunc(i) = max(min(scaldfunc(i), 1.0), 0.)
            else
                scaldfunc(i) = 1.0
            endif
        endif
    enddo

!> - Compute final scale-aware downdraft mass flux

    cfl_fac(:) = 1.0     ! ZNT 06/06/2021: column CFL limiter
    do k = kmscu, 1, -1
        do i = 1, im
            if(cnvflg(i) .AND. &
            (k >= mrad(i) .AND. k < krad(i))) then
                xmfd(i,k) = scaldfunc(i) * xmfd(i,k)
                dz   = zl(i,k+1) - zl(i,k)
                xmmx = dz / dt2 * cfl_crit      ! ZNT 09/24/2020: max MF allowed by CFL
                if (do_cfl_col) then            ! ZNT 06/06/2021: column CFL limiter
                    if (xmmx < xmfd(i,k)) then
                        cfl_fac(i) = min(cfl_fac(i), xmmx/xmfd(i,k))
                    endif
                else                            ! ZNT 06/06/2021: local CFL limiter
                    xmfd(i,k) = min(xmfd(i,k),xmmx)
                endif
            endif
        enddo
    enddo
    if (do_cfl_col) then ! ZNT 06/06/2021: column CFL limiter
        do k = kmscu, 1, -1
            do i = 1, im
                xmfd(i,k) = xmfd(i,k)*cfl_fac(i)
            enddo
        enddo
    endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> - Compute downdraft property using updated entranment rate

    do i = 1, im
        if(cnvflg(i)) then
            k = krad(i)
            thld(i,k)= thlx(i,k)
            hld(i,k) = hlx(i,k)    ! ZNT 06/12/2023: use h_l instead of th_l for better conservation
        endif
    enddo

!     do i = 1, im
!       if(cnvflg(i)) then
!         k = krad(i)
!         ptem1 = max(qcdo(i,k,ntcw), 0.)
!         tld = thld(i,k) / pix(i,k)
!         tcdo(i,k) = tld +  elocp * ptem1
!         qcdo(i,k,1) = qcdo(i,k,1)+0.2*qcdo(i,k,1)
!         qcdo(i,k,ntcw) = qcdo(i,k,ntcw)+0.2*qcdo(i,k,ntcw)
!       endif
!     enddo

    do k = kmscu,1,-1
        do i=1,im
            if(cnvflg(i) .AND. &
            (k >= mrad(i) .AND. k < krad(i))) then
                dz = zl(i,k+1) - zl(i,k)
                tem  = 0.5 * xlamde(i,k) * dz
                factor = 1. + tem
            
                thld(i,k) = ((1.-tem)*thld(i,k+1)+tem* &
                            (thlx(i,k)+thlx(i,k+1)))/factor
                hld(i,k)  = ((1.-tem)*hld(i,k+1)+tem* &
                            (hlx(i,k)+hlx(i,k+1)))/factor
                qtd(i,k)  = ((1.-tem)*qtd(i,k+1)+tem* &
                            (qtx(i,k)+qtx(i,k+1)))/factor
                qnd(i,k)  = ((1.-tem)*qnd(i,k+1)+tem* &
                            (qnx(i,k)+qnx(i,k+1)))/factor
            
                if(use_hl==1 .OR. use_hl==3) then      ! ZNT 06/12/2023: use h_l instead of th_l for better conservation
                    tld = hld(i,k) - zl(i,k)*gocp
                else
                    tld = thld(i,k) / pix(i,k)
                endif
                if (do_cmp_qs) then  ! ZNT 10/09/2020: note es is not used
                    call compute_qs(tld, plyr(i,k)*100.0, qs)
                    qs = max(qmin,qs)
                else
                    call lookup_es(tld, es)
                    es = 0.01*es             ! ZNT 04/27/2020
                ! es = 0.01 * fpvs(tlu)      ! fpvs in pa
                    qs = max(qmin, eps * es / (plyr(i,k)+epsm1*es))
                endif
                dq = qtd(i,k) - qs
            
                if (dq > 0.) then
                    gamma = el2orc * qs / (tld**2)
                    qld = dq / (1. + gamma)
                !              qtd(i,k) = qs + qld
                    qs = qtd(i,k) - qld
                    qcdo(i,k,1) = qs
                    qcdo(i,k,ntcw) = qld
                    qcdo(i,k,ntaw) = 1.      ! ZNT 06/30/2023: adjust cloud amount
                    tcdo(i,k) = tld + elocp * qld
                else
                    qnd(i,k) = 0. 
                    qcdo(i,k,1) = qtd(i,k)
                    qcdo(i,k,ntcw) = 0.
                    qcdo(i,k,ntaw) = 0.      ! ZNT 06/30/2023: adjust cloud amount
                    tcdo(i,k) = tld
                endif
            
                ! ZNT 06/30/2023: No aerosol activation for downdraft
            
                qcdo(i,k,ntnw) = qnd(i,k)    ! copy the final qn_u value to the tracer output array
            endif
        enddo
    enddo

    do k = kmscu, 1, -1
        do i = 1, im
            if (cnvflg(i) .AND. k < krad(i)) then
                if(k >= mrad(i)) then
                    dz = zl(i,k+1) - zl(i,k)
                    tem  = 0.5 * xlamdem(i,k) * dz
                    factor = 1. + tem
                    ptem = tem - pgcon
                    ptem1= tem + pgcon
                
                    ucdo(i,k) = ((1.-tem)*ucdo(i,k+1)+ptem*u1(i,k+1) &
                                +ptem1*u1(i,k))/factor
                    vcdo(i,k) = ((1.-tem)*vcdo(i,k+1)+ptem*v1(i,k+1) &
                                +ptem1*v1(i,k))/factor
                endif
            endif
        enddo
    enddo


    ! ZNT 06/30/2023: Change to full loop over all tracers and skip qc, qa, qn
    do n = 2, ntrac1
       if (n .eq. ntcw) cycle
       if (n .eq. ntaw) cycle    ! ZNT 06/30/2023: skip the passive transport for cloud amount
       if (n .eq. ntnw) cycle    ! ZNT 06/30/2023: skip the passive transport for qn
       do k = kmscu, 1, -1
          do i = 1, im
             if (cnvflg(i) .AND. k < krad(i)) then
                if(k >= mrad(i)) then
                   dz = zl(i,k+1) - zl(i,k)
                   tem  = 0.5 * xlamde(i,k) * dz
                   factor = 1. + tem
                        
                   qcdo(i,k,n) = ((1.-tem)*qcdo(i,k+1,n)+tem* &
                                 (q1(i,k,n)+q1(i,k+1,n)))/factor
                endif
             endif
          enddo
       enddo
    enddo

    return
    end subroutine mfscuq
