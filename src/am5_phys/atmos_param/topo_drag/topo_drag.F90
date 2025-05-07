module topo_drag_mod

!=======================================================================
! TOPOGRAPHIC DRAG CLOSURE -- Garner (2005)
!=======================================================================

!-----------------------------------------------------------------------
!  Calculates horizontal velocity tendency due to topographic drag
!-----------------------------------------------------------------------

use          mpp_mod, only: input_nml_file
use  mpp_domains_mod, only: domain2D
use          fms_mod, only: error_mesg, FATAL, NOTE, &
                            mpp_pe, mpp_root_pe, stdout, stdlog,       &
                            check_nml_error, write_version_number,lowercase
use      fms2_io_mod, only: read_data, get_variable_size, variable_exists, file_exists, &
                            FmsNetcdfDomainFile_t, register_variable_attribute, &
                            register_restart_field, register_axis, unlimited, &
                            open_file, read_restart, write_restart, close_file, &
                            register_field, write_data, get_global_io_domain_indices, &
                            FmsNetcdfFile_t
use    constants_mod, only: Grav, Cp_Air, Rdgas, Pi, Radian
use horiz_interp_mod, only: horiz_interp_type, horiz_interp_init, &
                            horiz_interp_new, horiz_interp, horiz_interp_del

implicit none

private

character(len=128) :: version = '$Id$'
character(len=128) :: tagname = '$Name$'

logical :: module_is_initialized = .false.

! horizontal array size

integer :: nlon, nlat
integer :: kd=0

! arrays defined by topo_drag_init:

real, allocatable, dimension(:,:) :: t11, t21, t12, t22    ! drag tensor
real, allocatable, dimension(:,:) :: hmin, hmax
real, allocatable, dimension(:,:) :: lon, lat

! parameters:

real, parameter :: u0=1.0       ! arbitrary velocity scale for diagnostics
real, parameter :: xl=80.0e3    ! arbitrary horiz length scale for diagnostics
real, parameter :: ro=1.2       ! arbitrary density scale for diagnostics
real, parameter :: lapse=Grav/Cp_Air ! adiabatic temperature lapse rate
real, parameter :: tiny=1.0e-20
real, parameter :: Tpi = 2.*Pi
real, parameter :: resolution=60.0 ! # of points per degree in topo datasets
real, parameter :: frint=0.5

integer, parameter :: ipts=360*resolution
integer, parameter :: jpts=180*resolution

type(domain2D), pointer     :: topo_domain

! parameters in namelist (topo_drag_nml):

real :: &
   frcrit=0.7   &      ! critical value of Froude number for nonlinear flow
  ,alin=0.9     &      ! amplitude of propagating drag
  ,anonlin=3    &      ! amplitude of nonpropagating drag
  ,gamma=0.4    &      ! exponent in aspect ratio power law
  ,epsi=0.0     &      ! exponent in distribution power law
  ,beta=0.5     &      ! bluntness of topographic features
  ,h_frac=0.1   &      ! ratio of min to max subgrid mountain height
  ,zref_fac=1.0 &      ! adjusts level separating breaking/laminar flow
  ,tboost=1.5   &      ! surface T boost to improve PBL height estimate
  ,pcut=0.0     &      ! high-level cutoff pressure for momentum forcing
  ,samp=2.0     &      ! correction for coarse sampling of d2v/dz2
  ,max_udt=3.e-3     & ! upper bound on acceleration [m/s2]
  ,no_drag_frac=0.10 & ! fraction of lower atmosphere with no breaking
  ,max_pbl_frac=0.50  ! max fraction of lower atmosphere in PBL
logical :: &
   do_conserve_energy=.true. & ! conserve total energy?
  ,keep_residual_flux=.true. & ! redistribute residual pseudomomentum?
  ,do_pbl_average=.false.    & ! average u,rho,N over PBL for baseflux?
  ,use_mg_scaling=.true.     & ! base flux saturates with value 'usat'?
  ,use_mask_for_pbl=.false.  & ! use bottom no_drag_layer as pbl?
  ,use_pbl_from_lock=.true.  & ! use pbl height from Lock boundary scheme? 
  ,use_uref_4stable=.true.   & ! uref instead of u(kdim) if stable PBL?
  ,write_restarts_and_stop=.false. ! if .true. the run stops after generating the topo_drag.res files 
  character(len=28)  :: interp_method = 'conservative'
NAMELIST /topo_drag_nml/                                               &
  frcrit, alin, anonlin, beta, gamma, epsi,                            &
  h_frac, zref_fac, tboost, pcut, samp, max_udt,                       &
  no_drag_frac, max_pbl_frac,                                          &
  do_conserve_energy, keep_residual_flux, do_pbl_average,              &
  use_mg_scaling, use_mask_for_pbl, use_pbl_from_lock,                 &
  use_uref_4stable, write_restarts_and_stop, interp_method

public topo_drag, topo_drag_init, topo_drag_end
public topo_drag_restart

contains

!#######################################################################

subroutine topo_drag (is, js, delt, uwnd, vwnd, atmp, &
                      pfull, phalf, zfull, zhalf, & 
                      lat, u_ref, v_ref, z_pbl, &
      dtaux, dtauy, dtaux_np, dtauy_np, dtemp, taux, tauy, taus, kbot )

integer, intent(in) :: is, js
real,    intent(in) :: delt
integer, intent(in), optional, dimension(:,:) :: kbot

! INPUT
! -----

! UWND     Zonal wind (dimensioned IDIM x JDIM x KDIM)
! VWND     Meridional wind (dimensioned IDIM x JDIM x KDIM)
! ATMP     Temperature at full levels (IDIM x JDIM x KDIM)
! PFULL    Pressure at full levels (IDIM x JDIM x KDIM)
! PHALF    Pressure at half levels (IDIM x JDIM x KDIM+1)
! ZFULL    Height at full levels (IDIM x JDIM x KDIM)
! ZHALF    Height at half levels (IDIM x JDIM x KDIM+1)

real, intent(in), dimension(:,:,:) :: uwnd, vwnd, atmp
real, intent(in), dimension(:,:)   :: lat, u_ref, v_ref, z_pbl
real, intent(in), dimension(:,:,:) :: pfull, phalf, zfull, zhalf

! OUTPUT
! ------

! DTAUX,DTAUY  Tendency of the vector wind in m/s^2 (IDIM x JDIM x KDIM)
! DTEMP        Tendency of the temperature in K/s (IDIM x JDIM x KDIM)
! TAUX,TAUY    Base momentum flux in kg/m/s^2 (IDIM x JDIM) for diagnostics
! TAUS         clipped saturation momentum flux (IDIM x JDIM x KDIM) for diagnostics

real, intent(out), dimension(:,:)   :: taux, tauy
real, intent(out), dimension(:,:,:) :: dtaux, dtauy, dtaux_np, dtauy_np, dtemp, taus

integer, dimension(size(zfull,1),size(zfull,2)) :: kpbl, knod, kcut
real,    dimension(size(zhalf,1),size(zhalf,2),size(zhalf,3)) :: tausat

! work arrays

real, dimension(size(zfull,1),size(zfull,2)) :: taub, taul, taup, taun
real, dimension(size(zfull,1),size(zfull,2)) :: frulo, fruhi, frunl, rnorm

integer :: idim
integer :: jdim
integer :: i,j, k, kdim, km
real    :: dz

  idim = size(uwnd,1)
  jdim = size(uwnd,2)
  kdim = size(uwnd,3)

! estimate height of pbl

  call get_pbl ( atmp, zfull, pfull, phalf, kpbl, knod, kcut )
!
  if (use_pbl_from_lock) then
     kpbl = kdim
     do k = kdim, 2, -1
       where ( zfull(:,:,k) < zhalf(:,:,kdim+1) + z_pbl(:,:) )
           kpbl(:,:) = k - 1           ! the first full model level above PBL
       endwhere
    enddo
  endif

! calculate base flux

  call base_flux (                                                     &
                                             is, js, uwnd, vwnd, atmp, &
                                             lat, u_ref, v_ref, z_pbl, &  !bqx+
                                             taux, tauy, dtaux, dtauy, &
                                               taub, taul, taup, taun, &
                                           frulo, fruhi, frunl, rnorm, &
                                     zfull, zhalf, pfull, phalf, kpbl )


! calculate saturation flux profile

  call satur_flux (                                                    &
                                                     uwnd, vwnd, atmp, &
                                                   taup, taub, tausat, &
                                                  frulo, fruhi, frunl, &
                        dtaux, dtauy, zfull, pfull, phalf, kpbl, kcut )

! calculate momentum tendency

  call topo_drag_tend (                                                &
                                               delt, uwnd, vwnd, atmp, &
                                       taux, tauy, taul, taun, tausat, &
  dtaux, dtauy, dtaux_np, dtauy_np, dtemp, zfull, zhalf, pfull, phalf, kpbl ) !bqx+ dtaux_np, dtauy_np

! put saturation flux profile into 'taus' for diagnostics

  do k=1,kdim
!     taus(:,:,k) = 0.5*rnorm(:,:)*(tausat(:,:,k) + tausat(:,:,k+1))
     taus(:,:,k) = tausat(:,:,k+1)*taub(:,:)/taul(:,:)
  enddo

! put total drag into 'taux,tauy' for diagnostics

  taup = taup - tausat(:,:,1)
  taub = (taup + taun)/taul
  taux = taux*taub
  tauy = tauy*taub

end subroutine topo_drag

!=======================================================================
                                  
subroutine base_flux (                                                 &
                                             is, js, uwnd, vwnd, atmp, &
                                             lat, u_ref, v_ref, z_pbl, & !bqx
                                             taux, tauy, dtaux, dtauy, &
                                               taub, taul, taup, taun, &
                                           frulo, fruhi, frunl, rnorm, &
                                     zfull, zhalf, pfull, phalf, kpbl )

integer, intent(in) :: is, js
real, intent(in),  dimension(:,:,:) :: uwnd, vwnd, atmp
real, intent(in),  dimension(:,:)   :: lat, u_ref, v_ref, z_pbl
real, intent(in),  dimension(:,:,:) :: zfull, zhalf, pfull, phalf
real, intent(out), dimension(:,:)   :: taux, tauy
real, intent(out), dimension(:,:,:) :: dtaux, dtauy
real, intent(out), dimension(:,:)   :: taub, taul, taup, taun
real, intent(out), dimension(:,:)   :: frulo, fruhi, frunl, rnorm
integer, intent(in), dimension(:,:) :: kpbl

real, dimension(size(uwnd,1),size(uwnd,2)) :: ubar, vbar

integer :: i, idim, id
integer :: j, jdim, jd
integer :: k, kdim, kb, kbp, kt, km

real :: usat, bfreq2, bfreq, dphdz, vtau, d2udz2, d2vdz2
real :: dzfull, dzhalf, dzhalf1, dzhalf2, density
real :: frmin, frmax, frmed, frumin, frumax, frumed, fruclp, fruclm
real :: rnormal, gterm, hterm, fru0, frusat
real :: usum, vsum, n2sum, delp

  idim = size(uwnd,1)
  jdim = size(uwnd,2)
  kdim = size(uwnd,3)

! compute base flux

  do j=1,jdim
     do i=1,idim
        usum = 0.
        vsum = 0.
        kt = kpbl(i,j)
        kb = max(kd,kt)
        do k=kt,kb
           delp = phalf(i,j,k+1) - phalf(i,j,k)
           usum = usum + uwnd(i,j,k)*delp
           vsum = vsum + vwnd(i,j,k)*delp
        enddo
        ubar(i,j) = usum/(phalf(i,j,kb+1) - phalf(i,j,kt))
        vbar(i,j) = vsum/(phalf(i,j,kb+1) - phalf(i,j,kt))
     enddo
  enddo

  if (use_uref_4stable) then
   where( z_pbl (:,:) == 0. )
    ubar(:,:) = u_ref(:,:)
    vbar(:,:) = v_ref(:,:)
   endwhere
  endif

  do j=1,jdim
     jd = js+j-1
     do i=1,idim
        id = is+i-1
        kt = kpbl(i,j)
        kb = max(kd,kt)
        kbp = min(kdim, kb+1) !bqx+
        dzfull = zhalf(i,j,kt) - zhalf(i,j,kb+1)
        density = (phalf(i,j,kb+1) - phalf(i,j,kt))/(Grav*dzfull)
        dzfull = zfull(i,j,kt-1) - zfull(i,j,kbp)
        bfreq2 = Grav*((atmp(i,j,kt-1) - atmp(i,j,kbp))/dzfull+lapse)/&
                  (0.5*(atmp(i,j,kt-1) + atmp(i,j,kbp)))
!
        bfreq = sqrt(max(tiny, bfreq2))

!       included 'alin' 4/2015

        taux(i,j) = (ubar(i,j)*t11(id,jd) + vbar(i,j)*t21(id,jd))      &
                                                   *bfreq*density
        tauy(i,j) = (ubar(i,j)*t12(id,jd) + vbar(i,j)*t22(id,jd))      &
                                                   *bfreq*density

        taub(i,j) = max(tiny, sqrt(taux(i,j)**2 + tauy(i,j)**2))

!       min/max Froude numbers based on low-level flow

        vtau = max(tiny, -(ubar(i,j)*taux(i,j)                         &
                         + vbar(i,j)*tauy(i,j))/taub(i,j))
        frmax = hmax(id,jd)*bfreq / vtau
        frmin = hmin(id,jd)*bfreq / vtau
        frmed = frcrit + frint

!       linear momentum flux associated with min/max Froude numbers

        dphdz = bfreq / vtau
        usat = sqrt(density/ro) * vtau / sqrt(dphdz*xl)
        frusat = frcrit*usat

        frumin = frmin*usat
        frumax = frmax*usat
        frumed = frmed*usat

        frumax = max(frumax,frumin + tiny)
        fruclp = min(frumax,max(frumin,frusat))
        fruclm = min(frumax,max(frumin,frumed))
        fru0 = (u0/vtau)*usat

!       total drag in linear limit

        rnormal =                                                      &
                   (frumax**(2.0*gamma - epsi)                         &
                  - frumin**(2.0*gamma - epsi))/(2.0*gamma - epsi)
        rnormal = fru0**gamma * ro/rnormal  

        taul(i,j) =                                                    &
                   (frumax**(2.0 + gamma - epsi)                       &
                  - frumin**(2.0 + gamma - epsi))/(2.0 + gamma - epsi)

!       separate propagating and nonpropagating parts of total drag

        gterm = frusat**(beta + 1.0)*                                  &
                (frumax**(gamma - epsi - beta)                         &
               - fruclp**(gamma - epsi - beta))/(gamma - epsi - beta)

        taup(i,j) =  alin *                                             &
                  ( (fruclp**(2.0 + gamma - epsi)                       &
                   - frumin**(2.0 + gamma - epsi))/(2.0 + gamma - epsi) &
                                                       + frusat*gterm )

        taun(i,j) = anonlin*usat/(1.0 + beta) *                        &
                 ( (frumax**(1.0 + gamma - epsi)                       &
                  - fruclp**(1.0 + gamma - epsi))/(1.0 + gamma - epsi) &
                                                      - gterm )

!       5/2015 mg option: depth of blocking ~ U/N, not h

        if (use_mg_scaling) taun(i,j) = taun(i,j)/max(frmax,frcrit)

        fruhi(i,j) = frumax
        frulo(i,j) = frumin
        frunl(i,j) = frusat
        rnorm(i,j) = rnormal

     enddo
  enddo

! wind component opposite the drag at full levels (stored as 'dtaux')

  do k=1,kdim
     do j=1,jdim
        do i=1,idim
           dtaux(i,j,k) =                                              &
             -(uwnd(i,j,k)*taux(i,j) + vwnd(i,j,k)*tauy(i,j))/taub(i,j)
        enddo
     enddo
  enddo

! curvature of wind at full levels (stored as 'dtauy')

  dtauy = 0.

  do j=1,jdim
     do i=1,idim
!        kt = kpbl(i,j)
        kt = min(kdim-1, kpbl(i,j)) !bqx+
        do k=2,kt
           dzfull = zhalf(i,j,k) - zhalf(i,j,k+1)
           dzhalf1 = zfull(i,j,k-1) - zfull(i,j,k)
           dzhalf2 = zfull(i,j,k) - zfull(i,j,k+1)
           d2udz2 = ((uwnd(i,j,k-1) - uwnd(i,j,k  ))/dzhalf1           &
                   - (uwnd(i,j,k  ) - uwnd(i,j,k+1))/dzhalf2)/dzfull
           d2vdz2 = ((vwnd(i,j,k-1) - vwnd(i,j,k  ))/dzhalf1           &
                   - (vwnd(i,j,k  ) - vwnd(i,j,k+1))/dzhalf2)/dzfull
           dtauy(i,j,k) = -(d2udz2*taux(i,j) + d2vdz2*tauy(i,j))/      &
                                                              taub(i,j)
        enddo

     enddo
  enddo

end subroutine base_flux

!=======================================================================

subroutine satur_flux (                                                &
                                                     uwnd, vwnd, atmp, &
                                                   taup, taub, tausat, &
                                                  frulo, fruhi, frunl, &
                        dtaux, dtauy, zfull, pfull, phalf, kpbl, kcut )

real, intent(in),  dimension (:,:,:) :: uwnd, vwnd, atmp
real, intent(in),  dimension (:,:,:) :: dtaux, dtauy
real, intent(in),  dimension (:,:,:) :: zfull, pfull, phalf
real, intent(in),  dimension (:,:)   :: taup
real, intent(out), dimension (:,:)   :: taub
real, intent(out), dimension (:,:,:) :: tausat
real, intent(in),  dimension (:,:)   :: frulo, fruhi, frunl
integer, intent(in), dimension (:,:) :: kpbl, kcut

real, dimension(size(zfull,1),size(zfull,2)) :: usat

real :: dzhalf, gterm, gterm0, density
real :: bfreq2, bfreq, vtau, d2vtau, dphdz, xl1
real :: frumin, frumax, fruclp, frusat, frusat0, fruclp0

integer :: i, idim
integer :: j, jdim
integer :: k, kdim, k1

  idim = size(uwnd,1)
  jdim = size(uwnd,2)
  kdim = size(uwnd,3)

! get vertical profile of propagating part of momentum flux

  usat = frunl/frcrit

  do k=kdim,2,-1
     do j=1,jdim
        do i=1,idim

!          buoyancy frequency, velocity and density at half levels

           dzhalf = zfull(i,j,k-1) - zfull(i,j,k)
           density = (pfull(i,j,k) - pfull(i,j,k-1))/(Grav*dzhalf)
           bfreq2 = Grav*                                              &
                       ((atmp(i,j,k-1) - atmp(i,j,k))/dzhalf + lapse)/ & 
                   (0.5*(atmp(i,j,k-1) + atmp(i,j,k)))
           bfreq = sqrt(max(tiny, bfreq2))
           
           vtau = max(tiny, 0.5*(dtaux(i,j,k-1) + dtaux(i,j,k)))

!          WKB correction of vertical wavelength

           d2vtau = 0.5*(dtauy(i,j,k-1) + dtauy(i,j,k))
           xl1 = xl*max(0.5, min(2.0, 1.0 - samp*vtau*d2vtau/(bfreq*bfreq)))

!          min/max and critical momentum flux values at half levels

           dphdz = bfreq / vtau
           usat(i,j) = min(usat(i,j),sqrt(density/ro) * vtau/sqrt(dphdz*xl1))
           frusat = frcrit*usat(i,j)

           frumin = frulo(i,j)
           frumax = fruhi(i,j)
           fruclp = min(frumax,max(frumin,frusat))
           frusat0 = frunl(i,j)
           fruclp0 = min(frumax,max(frumin,frusat0))

!          propagating part of momentum flux (from WKB or EP)

           gterm0 = (frumax**(gamma - epsi - beta)                     &
                - fruclp0**(gamma - epsi - beta))/(gamma - epsi - beta)
           gterm = (fruclp0**(gamma - epsi)                            &
                               - fruclp**(gamma - epsi))/(gamma - epsi)

           tausat(i,j,k) = alin *                                      &
                 ( (fruclp**(2.0 + gamma - epsi)                       &
                  - frumin**(2.0 + gamma - epsi))/(2.0 + gamma - epsi) &
                         + frusat**2.0*(gterm0*frusat0**beta + gterm) )
        enddo
     enddo
  enddo

! make propagating flux constant with height in zero-drag top layer
! changed 5/2014

  k1 = maxval(kcut)
  do k=k1,1,-1
     where (k <= kcut)
        tausat(:,:,k) = tausat(:,:,k+1)
     endwhere
  enddo

! make propagating flux constant with height in zero-drag surface layer

  k1 = minval(kpbl)
  do k=kdim+1,k1+1,-1
     where (k > kpbl)
        tausat(:,:,k) = taup
     endwhere
  enddo

! redistribute residual forcing

  if ( keep_residual_flux ) then
     taub(:,:) = tausat(:,:,1)/(phalf(:,:,kdim+1) - phalf(:,:,1))
     do k=1,kdim
        tausat(:,:,k) = tausat(:,:,k)                                  &
                        - taub(:,:)*(phalf(:,:,kdim+1) - phalf(:,:,k))
     enddo
  endif

endsubroutine satur_flux

!=======================================================================

subroutine topo_drag_tend (                                            &
                                               delt, uwnd, vwnd, atmp, &
                                       taux, tauy, taul, taun, tausat, &
                dtaux, dtauy, dtaux_np, dtauy_np, dtemp, zfull, zhalf, &
                                                    pfull, phalf, kpbl )

real, intent(in) :: delt
real, intent(in), dimension(:,:,:)    :: uwnd, vwnd, atmp
real, intent(in), dimension(:,:,:)    :: zfull, zhalf, pfull, phalf
real, intent(in), dimension(:,:)      :: taux, tauy, taul, taun
real, intent(inout), dimension(:,:,:) :: tausat
real, intent(inout), dimension(:,:,:) :: dtaux, dtauy, dtemp
real, intent(out), dimension(:,:,:)   :: dtaux_np, dtauy_np !bqx
integer, intent(in), dimension (:,:)  :: kpbl

real, parameter :: bfmin=0.7e-2, bfmax=1.7e-2  ! min/max buoyancy freq [1/s]
real, parameter :: vvmin=1.0                   ! minimum surface wind [m/s]

integer,dimension(size(zfull,1),size(zfull,2)) :: kref
real :: dzhalf, zlast, rscale, phase, bfreq, bfreq2, vtau
real :: gfac, gfac1, dp, weight, wtsum, taunon, taunon1

integer :: i, idim
integer :: j, jdim
integer :: k, kdim, kr, kt

real,dimension(size(zfull,1),size(zfull,2)) :: dx, dy

  idim = size(uwnd,1)
  jdim = size(uwnd,2)
  kdim = size(uwnd,3)

! find reference level for non-propagating drag (z ~ pi U/N)

!the following re-orients the drag to align with low-level wind
!do j=1,jdim
!do i=1,idim
!k = knod(i,j)
!gfac = sqrt( (taux(i,j)**2 + tauy(i,j)**2) / max (tiny, uwnd(i,j,k)**2 + vwnd(i,j,k)**2) )
!dx(i,j) = gfac * uwnd(i,j,k)
!dy(i,j) = gfac * vwnd(i,j,k)
!enddo
!enddo

  do j=1,jdim
     do i=1,idim
        k = kpbl(i,j)
!stg        k = kdim
        phase = 0.0
        zlast = zhalf(i,j,k)
        do while (phase <= Pi*zref_fac .and. k > 1)
           k = k-1
           vtau = 0.5*(dtaux(i,j,k-1) + dtaux(i,j,k))
           dzhalf = zfull(i,j,k-1) - zfull(i,j,k)
           bfreq2 = Grav*                                              &
                       ((atmp(i,j,k-1) - atmp(i,j,k))/dzhalf + lapse)/ &
                   (0.5*(atmp(i,j,k-1) + atmp(i,j,k)))
           bfreq = sqrt(max(tiny, bfreq2))
           rscale = max(bfmin, min(bfmax, bfreq))/max(vvmin, vtau)
           dzhalf = zfull(i,j,k-1) - zlast
           phase = phase + dzhalf*rscale
           zlast = zfull(i,j,k-1)
        enddo
        kref(i,j) = k
     enddo
  enddo

! CALCULATE DECELERATION DUE TO PROPAGATING DRAG (~-rho^-1 dtau/dz)

  do k=1,kdim
     do j=1,jdim
        do i=1,idim
          dp = phalf(i,j,k+1) - phalf(i,j,k)
          gfac = tausat(i,j,k+1) - tausat(i,j,k)
          gfac1 = gfac*Grav/(dp*taul(i,j))
          dtaux(i,j,k) = gfac1*taux(i,j)
          dtauy(i,j,k) = gfac1*tauy(i,j)
!dtaux(i,j,k) = -gfac1*dx(i,j)
!dtauy(i,j,k) = -gfac1*dy(i,j)
        enddo
     enddo
  enddo

! CALCULATE DECELERATION DUE TO NON-PROPAGATING DRAG
  dtaux_np = 0. !bqx
  dtauy_np = 0. !bqx

  do j=1,jdim
     do i=1,idim
        kr = kref(i,j)
        kt = kpbl(i,j)
!stg        kt = kdim
        gfac = taun(i,j)/taul(i,j) * Grav
        wtsum = 0.0
        do k=kr,kt
           dp = phalf(i,j,k+1) - phalf(i,j,k)
           weight = pfull(i,j,k) - phalf(i,j,kr)
           wtsum = wtsum + dp*weight
        enddo
        taunon = 0.0
        do k=kr,kt
           weight = pfull(i,j,k) - phalf(i,j,kr)
           gfac1 = gfac*weight/wtsum
           dtaux(i,j,k) = dtaux(i,j,k) + gfac1*taux(i,j)
           dtauy(i,j,k) = dtauy(i,j,k) + gfac1*tauy(i,j)
!bqx
           dtaux_np(i,j,k) = gfac1*taux(i,j)
           dtauy_np(i,j,k) = gfac1*tauy(i,j)
!dtaux(i,j,k) = dtaux(i,j,k) - gfac1*dx(i,j)
!dtauy(i,j,k) = dtauy(i,j,k) - gfac1*dy(i,j)
           dp = phalf(i,j,k+1) - phalf(i,j,k)
           taunon = taunon + gfac1*dp
           taunon1 = taunon*taul(i,j)/Grav
           tausat(i,j,k) = tausat(i,j,k) + taunon1
        enddo
        do k=kt+1,kdim
           tausat(i,j,k) = tausat(i,j,k) + taunon1
        enddo
     enddo
  enddo

  dtaux = max(-max_udt, min(max_udt, dtaux))   !stg
  dtauy = max(-max_udt, min(max_udt, dtauy))

! CALCULATE HEATING TO CONSERVE TOTAL ENERGY

  if (do_conserve_energy) then
     dtemp = -((uwnd + 0.5*delt*dtaux)*dtaux                           &
             + (vwnd + 0.5*delt*dtauy)*dtauy)/Cp_Air
  else
     dtemp = 0.0
  endif

end subroutine topo_drag_tend

!=======================================================================

subroutine get_pbl ( atmp, zfull, pfull, phalf, kpbl, knod, kcut )

integer, intent(out), dimension(:,:) :: kpbl, knod, kcut
real, intent(in), dimension(:,:,:)   :: atmp
real, intent(in), dimension(:,:,:)   :: zfull, pfull, phalf

real, dimension(size(pfull,1),size(pfull,2)) :: ppbl, pbot
real, dimension(size(pfull,1),size(pfull,2)) :: tbot, zbot

integer :: i, idim
integer :: j, jdim
integer :: k, kdim

  idim = size(atmp,1)
  jdim = size(atmp,2)
  kdim = size(atmp,3)

  do j=1,jdim
     do i=1,idim
        ppbl(i,j) = (1.0 - max_pbl_frac)*phalf(i,j,kdim+1)
        pbot(i,j) = (1.0 - no_drag_frac)*phalf(i,j,kdim+1)
        tbot(i,j) = atmp(i,j,kdim) + tboost
        zbot(i,j) = zfull(i,j,kdim)
     enddo
  enddo

! find highest model level in no-drag surface layer

  knod = kdim-1

  do k=kdim-2,2,-1
     where ( pfull(:,:,k) >= pbot(:,:) )
        knod = k
     endwhere
  enddo

! find lowest model level in no-drag top layer

  kcut = 1

  do k=2,kdim
     where ( pfull(:,:,k) <= pcut )
        kcut = k
     endwhere
  enddo

  if (kd == 0 .and. do_pbl_average) kd = kdim-1

  if ( use_mask_for_pbl ) then
     kpbl = knod
     return
  endif

! find the first layer above PBL

  kpbl = kdim-1

  do k=kdim-2,2,-1
     where ( pfull(:,:,k) >= ppbl(:,:) .and.                          &
        tbot(:,:) - atmp(:,:,k) > lapse*(zfull(:,:,k) - zbot(:,:)) )
        kpbl = k -1
     endwhere
  enddo

end subroutine get_pbl

subroutine map_horiz_grid (                                            & 
                                   lon_des, lat_des, lon_src, lat_src, &
                                                   ib2, ie2, jb2, je2, &
                                                                 fold )

real, intent(in), dimension(:,:) :: lon_des, lat_des
real, intent(in), dimension(:)   :: lon_src, lat_src
integer, intent(out) :: ib2, ie2, jb2, je2
logical, intent(out) :: fold

!-----------------------------------------------------------------------
! local allocations
!-----------------------------------------------------------------------

real, dimension(size(lon_des,1),size(lon_des,2)) :: lon
real :: slat, nlat, wlon, elon
integer :: idim, jdim
logical :: pole

  idim = size(lon_des,1)
  jdim = size(lon_des,2)

  slat = minval(lat_des)
  nlat = maxval(lat_des)

  pole = ((lat_des(2,2)-lat_des(1,1))*                                 &
          (lat_des(idim,jdim)-lat_des(idim-1,jdim-1)) < 0)

  fold = (abs(lon_des(idim,jdim) - lon_des(1,1)) > pi)

  idim = size(lon_src)
  jdim = size(lat_src)

  if (pole) then

     call get_indices ( jdim, jb2, je2, lat_src, slat, nlat )
     ib2 = 1 ; ie2 = idim
     je2 = jdim

  else

     lon = lon_des

     if (fold) where (lon < pi) lon = lon + Tpi

     wlon = minval(lon)
     elon = mod(maxval(lon), Tpi)

     call get_indices ( idim, ib2, ie2, lon_src, wlon, elon )
     call get_indices ( jdim, jb2, je2, lat_src, slat, nlat )

  endif

  if ( fold .and. .NOT. pole ) ie2 = ie2 + ipts

  return
end subroutine map_horiz_grid

!#######################################################################

subroutine get_indices (                                               &
                                                     ndim, nbeg, nend, &
                                 axis_src, axis_des_beg, axis_des_end )

integer,                    intent (in)  :: ndim
integer,                    intent (out) :: nbeg, nend
real, dimension(ndim),      intent (in)  :: axis_src
real,                       intent (in)  :: axis_des_beg, axis_des_end

!-----------------------------------------------------------------------
! local allocations
!-----------------------------------------------------------------------

  integer, parameter :: npad=3
  integer :: n

  do n=1,ndim
     nbeg = n-1
     if ( axis_src(n) > axis_des_beg ) exit
  enddo

  do n=ndim,1,-1
     nend = n+1
     if ( axis_src(n) < axis_des_end ) exit
  enddo

  if ( nbeg == 0 .or. nend == ndim+1 ) &
     call error_mesg('topo_drag_mod','destination grid is not inside source grid', FATAL)

  nbeg = max(1,nbeg-npad)
  nend = min(ndim,nend+npad)

  return
end subroutine get_indices

!=======================================================================
 subroutine topo_drag_register_tile_restart(restart)
     type(FmsNetcdfDomainFile_t), intent(inout) :: restart
     character(len=8), dimension(3)             :: dim_names

     dim_names(1) = "xaxis_1"
     dim_names(2) = "yaxis_1"
     dim_names(3) = "Time"
     call register_axis(restart, dim_names(1), "x")
     call register_axis(restart, dim_names(2), "y")
     call register_axis(restart, dim_names(3), unlimited)

     !< Register the domain decomposed dimensions as variables so that the combiner can work
     !! correctly
     call register_field(restart, dim_names(1), "double", (/dim_names(1)/))
     call register_field(restart, dim_names(2), "double", (/dim_names(2)/))

     call register_restart_field(restart, "t11", t11, dim_names)
     call register_restart_field(restart, "t12", t12, dim_names)
     call register_restart_field(restart, "t21", t21, dim_names)
     call register_restart_field(restart, "t22", t22, dim_names)
     call register_restart_field(restart, "hmin", hmin, dim_names)
     call register_restart_field(restart, "hmax", hmax, dim_names)

end subroutine

subroutine topo_drag_init (domain, lonb, latb)

type(domain2D), target, intent(in) :: domain
real, intent(in), dimension(:,:) :: lonb, latb

character(len=128) :: msg
character(len=64)  :: restart_fname='INPUT/topo_drag.res.nc'
character(len=64)  :: topography_file='INPUT/poztopog.nc'
character(len=64)  :: dragtensor_file='INPUT/dragelements.nc'
character(len=3)   :: tensornames(4) = (/ 't11', 't21', 't12', 't22' /)

logical :: found_field(4), fold

real, parameter :: bfscale=1.0e-2      ! buoyancy frequency scale [1/s]

real, parameter :: flonmin = 0.
real, parameter :: flonmax = 360.
real, parameter :: flatmin = -90.
real, parameter :: flatmax = 90.

real, allocatable, dimension(:)   :: flonb, flatb, glonb, glatb
real, allocatable, dimension(:,:) :: clonb, clatb
real, allocatable, dimension(:,:) :: zdat, zout, zoutb
type (horiz_interp_type) :: Interp
real :: exponent

integer :: n, nfold
integer :: io, ierr, unit_nml, logunit
integer :: i, j
integer :: ib, ie, jb, je, ib2, ie2, jb2, je2
integer :: siz(2)
type(FmsNetcdfDomainFile_t) :: Topo_restart !< Fms2io domain decomposed fileobj
type(FmsNetcdfFile_t) :: topography_fileobj, dragtensor_fileobj !< Fms2io fileobj

  if (module_is_initialized) return

  nlon = size(lonb,1)-1
  nlat = size(latb,2)-1

! read namelist

   read (input_nml_file, nml=topo_drag_nml, iostat=io)
   ierr = check_nml_error(io,'topo_drag_nml')

! write version number and namelist to logfile

  call write_version_number (version, tagname)
  logunit = stdlog()
  if (mpp_pe() == mpp_root_pe())                                       &
                                    write (logunit, nml=topo_drag_nml)

  topo_domain => domain
  allocate (t11(nlon,nlat))
  allocate (t21(nlon,nlat))
  allocate (t12(nlon,nlat))
  allocate (t22(nlon,nlat))
  allocate (hmin(nlon,nlat))
  allocate (hmax(nlon,nlat))

  if (gamma == beta + epsi) gamma = gamma + tiny

! read restart file

  if ( open_file(Topo_restart, restart_fname, "read", topo_domain, is_restart = .true.) ) then

     if (mpp_pe() == mpp_root_pe()) then
        write ( msg, '("Reading restart file: ",a40)' ) restart_fname
        call error_mesg('topo_drag_mod', msg, NOTE)
     endif

     call topo_drag_register_tile_restart(Topo_restart)
     call read_restart(Topo_restart)
     call close_file(Topo_restart)

  else if (file_exists(topography_file) .and.                           &
           file_exists(dragtensor_file)) then

!    read and interpolate topography datasets

     ! check for correct field size in topography file
     if (mpp_pe() == mpp_root_pe()) then
        write ( msg, '("Reading topography file: ",a)') trim(topography_file)
        call error_mesg('topo_drag_mod', msg, NOTE)
     endif

     if( .not. open_file(topography_fileobj, topography_file, "read")) then
         call error_mesg('topo_drag_mod', "Error opening topography file", FATAL)
     endif

     if( .not. open_file(dragtensor_fileobj, dragtensor_file, "read")) then
         call error_mesg('topo_drag_mod', "Error opening dragtensor file", FATAL)
     endif

     ! check for correct field size in topography
     call get_variable_size(topography_fileobj, 'hpoz', siz)
     if (siz(1) /= ipts .or. siz(2) /= jpts) then
         call error_mesg('topo_drag_mod', 'Field \"hpoz\" in file '//  &
                   trim(topography_file)//' has the wrong size', FATAL)
     endif

     ! check for correct field size in tensor file

     if (mpp_pe() == mpp_root_pe()) then
        write ( msg, '("Reading drag tensor file: ",a)')               &
                                                trim(dragtensor_file)
        call error_mesg('topo_drag_mod', msg, NOTE)
     endif
     call get_variable_size(dragtensor_fileobj, tensornames(1), siz)
     if (siz(1) /= ipts .or. siz(2) /= jpts) then
         call error_mesg('topo_drag_mod', 'Fields in file ' &
         //trim(dragtensor_file)//' have the wrong size', FATAL)
     endif

     allocate (zout(1:nlon,1:nlat))
     allocate (zoutb(1:nlon+1,1:nlat+1))
     allocate (flonb(ipts+1))
     allocate (flatb(jpts+1))

     do i=1,ipts+1
        flonb(i) = (i-1)/resolution / Radian
     enddo
     do j=1,jpts+1
        flatb(j) = (-90.0 + (j-1)/resolution) / Radian
     enddo

     ! initialize horizontal interpolation
     call horiz_interp_init

     ! new data input procedure
     call map_horiz_grid (lonb, latb, flonb, flatb, ib2, ie2, jb2, je2, fold )

        if (fold) then
           nfold = ipts
        else
           nfold = 0
        endif

        allocate (glonb(ib2:ie2))
        allocate (glatb(jb2:je2))

        do i=ib2,ie2
           glonb(i) = (i-1)/resolution / Radian
        enddo
        do j=jb2,je2
           glatb(j) = (-90.0 + (j-1)/resolution) / Radian
        enddo
 
        allocate (clonb(nlon+1,nlat+1))
        allocate (clatb(nlon+1,nlat+1))

        clonb = lonb
        clatb = latb
        if (fold) then
           where (clonb(:,:) < pi) clonb(:,:) = clonb(:,:) + Tpi
        endif

        ! initialize horizontal interpolation

        call horiz_interp_new ( Interp, glonb, glatb,                &
                         clonb, clatb, interp_method=interp_method)

        allocate (zdat(ib2:ie2-1,jb2:je2-1))
        ! read relief file
        call nc_read_util ( topography_file, 'hpoz', zdat,          &
                            nfold, ib2, ie2-1, jb2, je2-1 )

        exponent = 2. - gamma
        zdat = max(0., zdat)**exponent
        ! regrid hpoz
        if (lowercase(trim(interp_method)) == 'bilinear') then
          call horiz_interp( Interp, zdat, zoutb)
          zout(1:nlon,1:nlat)=zoutb(1:nlon,1:nlat)
        else
          call horiz_interp( Interp, zdat, zout)
        endif
        hmax = (abs(zout)*(gamma + 2.)/(2.*gamma) * &
          (1. - h_frac**(2.*gamma))/(1. - h_frac**(gamma + 2.)))**(1.0/exponent)
        hmin = hmax*h_frac
        if(write_restarts_and_stop) print*,'topo_drag maxvals of zdat,zout',maxval(zdat),maxval(zout)

        ! read tensor file

        do n=1,4
           found_field(n) = variable_exists(dragtensor_fileobj, tensornames(n))
           if (.not. found_field(n)) cycle

           call nc_read_util ( dragtensor_file, tensornames(n), zdat,  &
                               nfold, ib2, ie2-1, jb2, je2-1 )

           ! regrid tensor elements
           if (lowercase(trim(interp_method)) == 'bilinear') then
             call horiz_interp( Interp, zdat, zoutb)
             zout(1:nlon,1:nlat)=zoutb(1:nlon,1:nlat)
           else
             call horiz_interp( Interp, zdat, zout)
           endif

           if ( tensornames(n) == 't11' ) then
              t11 = zout/bfscale
           else if ( tensornames(n) == 't21' ) then
              t21 = zout/bfscale
           else if ( tensornames(n) == 't12' ) then
              t12 = zout/bfscale
           else if ( tensornames(n) == 't22' ) then
              t22 = zout/bfscale
           endif
        enddo

        deallocate (glonb, glatb)
        deallocate (clonb, clatb)

     if (.not. found_field(3)) t12 = t21

     deallocate (flonb, flatb)
     deallocate (zdat, zout, zoutb)
     call horiz_interp_del ( Interp )

     call close_file(topography_fileobj)
     call close_file(dragtensor_fileobj)

     if(write_restarts_and_stop) then
       call topo_drag_restart
       call error_mesg ('topo_drag_init','Write topo_drag restarts and exit!', NOTE)
       stop
     endif
  else
     call error_mesg ('topo_drag_init','No sub-grid orography available for topo_drag', FATAL)

  endif

  module_is_initialized = .true.

end subroutine topo_drag_init

!=======================================================================

subroutine topo_drag_end

! writes static arrays to restart file

  if (mpp_pe() == mpp_root_pe() ) then
     call error_mesg('topo_drag_mod', 'Writing netCDF formatted restart file: RESTART/topo_drag.res.nc', NOTE)
  endif

  call topo_drag_restart

  module_is_initialized = .false.

end subroutine topo_drag_end

!#######################################################################
subroutine topo_drag_restart(timestamp)
! write out restart file.
! Arguments: 
!   timestamp (optional, intent(in)) : A character string that represents the model time, 
!                                      used for writing restart. timestamp will append to
!                                      the any restart file name as a prefix. 

      character(len=*), intent(in), optional :: timestamp
      type(FmsNetcdfDomainFile_t) :: Topo_restart
      character(len=128)  :: restart_fname

      if (present(timestamp)) then
          restart_fname='RESTART/'//trim(timestamp)//'.topo_drag.res.nc'
      else
          restart_fname='RESTART/topo_drag.res.nc'
      endif

      if (.not. open_file(Topo_restart, restart_fname, "overwrite", topo_domain, is_restart = .true.)) then
         call error_mesg("topo_drag_mod", "The topo_drag tiled restart file does not exist", fatal)
      endif

      call topo_drag_register_tile_restart(Topo_restart)
      call write_restart(Topo_restart)
      call add_domain_dimension_data(Topo_restart)
      call close_file(Topo_restart)
end subroutine topo_drag_restart

!#######################################################################

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

!#######################################################################
!< nc_read_util: utility to read in hpoz data file
subroutine nc_read_util ( filename, name, var, lx, ib, ie, jb, je )
  use netcdf
  character(len=*),      intent (in)    :: filename, name
  real, dimension(:,:),  intent (out)   :: var
  integer,               intent (in)    :: lx, ib, ie, jb, je
  ! local allocations
  real, allocatable, dimension(:,:) :: var_west, var_east
  integer :: status, rstatus
  integer :: NCID, VARID
  integer :: start(3), count(3)
  integer :: lxwest, lxeast

  status = NF90_OPEN (filename, NF90_NOWRITE, NCID)

  if ( status /= 0 ) then
     call error_handler ('netcdf file not found', FATAL) 
  endif

  status = NF90_INQ_VARID (NCID, name, VARID)

  if ( status /= 0 ) then
     call error_handler ('netcdf variable not found', FATAL) 
  endif

  if ( lx == 0 ) then
     start(1) = ib
     start(2) = jb
     start(3) = 1
     count(1) = size(var,1)
     count(2) = size(var,2)
     count(3) = 1

     rstatus = NF90_GET_VAR ( NCID, VARID, var, start, count )
     if ( rstatus /= 0 ) then
        call error_handler ('read failed', FATAL)
     endif
  else
     lxeast = lx - ib + 1
     lxwest = ie - lx
     allocate ( var_west(lxwest, size(var,2)) )
     allocate ( var_east(lxeast, size(var,2)) )
     start(1) = 1
     start(2) = jb
     start(3) = 1
     count(1) = lxwest
     count(2) = size(var,2)  
     count(3) = 1
     rstatus = NF90_GET_VAR ( NCID, VARID, var_west, start, count )
     if ( rstatus /= 0 ) then
        call error_handler ('read failed east of fold', FATAL) 
     endif
     var(lxeast+1:lxeast+lxwest,:) = var_west
     start(1) = ib
     count(1) = lxeast
     rstatus = NF90_GET_VAR ( NCID, VARID, var_east, start, count )
     if ( rstatus /= 0 ) then
        call error_handler ('read failed west of fold', FATAL)
     endif
     var(1:lxeast,:) = var_east
     deallocate ( var_west, var_east )
  endif
  status = NF90_CLOSE (NCID)
  return
end subroutine nc_read_util

subroutine error_handler ( message, level )
  character(len=*), intent(in) :: message
  integer, intent(in) :: level
  call error_mesg ( 'topo_drag_mod', message, level )
  return
end subroutine error_handler

endmodule topo_drag_mod
