module gw_common_utils
  use FMSconstants, only: GRAV, CP_AIR, RDGAS, SECONDS_PER_DAY, PI
  use axis_utils2_mod, only: interp_1d

  implicit none
  public :: GWBand, gw_drag_prof, energy_momentum_adjust, gw_prof, &
              gw_common_init, gw_common_end
  private

!------Declare GWBand type----------------------------
  type :: GWBand
    integer           :: ngwv    !< Dimension of the spectrum.
    real              :: dc      !< Delta between nearest phase speeds [m/s].
    real, allocatable :: cref(:) !< Reference speeds [m/s].
    real              :: effkwv  !< Effective horizontal wave number [1/m] (fcrit2*kwv).
    logical         :: is_set !< Whether it has been setup
    contains
      procedure :: create_new_GWband
  end type GWBand

!--------------------------------------------------------------------
!------ private data ------
  integer               :: ktop = 1 !< index for model top layer
  real , parameter :: dback = 0.05   !< Background diffusivity.
  real , parameter :: taumin = 1.e-10 !< Minimum non-zero stress.
  real , parameter :: umcfac = 0.5    !< Maximum allowed change in u-c
                                                           !! (before efficiency applied).
  real , parameter :: ubmc2mn = 0.01  !< Minimum value of (u-c)**2.

  logical            :: tau_0_ubc = .true.        !< Whether or not to enforce an upper
                                                  !! boundary condition of tau = 0. ()

 ! Index the cardinal directions.
integer, parameter :: west = 1
integer, parameter :: east = 2
integer, parameter :: south = 3
integer, parameter :: north = 4

logical                    :: module_is_initialized = .false.
! Newtonian cooling coefficient [1/s]
real,  allocatable    :: alpha(:)
  contains

  !=================================================
   !> @brief  Constructor for a GWBand that calculates derived components.
  !! Used directly to set the type's components.
  !! @code{.F90}
  !! type(GW_band) :: beres_band
  !! call beres_band%create_new_GWband(ngwv, dc, fcrit2, wavelength)
  !! @endcode
  subroutine create_new_GWband(this, ngwv, dc, fcrit2, wavelength)
    class(GWBand), intent(inout) :: this
    integer,       intent(in)    :: ngwv
    real,          intent(in)    :: dc
    real,          intent(in)    :: fcrit2
    real,          intent(in)    :: wavelength !< Wavelength in meters.

    integer :: l !< Wavenumber index.
    real :: kwv !<  horizontal wave number [1/m]

    ! Simple assignments.
    this%ngwv = ngwv
    this%dc = dc
   ! this%fcrit2 = fcrit2

    ! Uniform phase speed reference grid.
    allocate(this%cref(-ngwv:ngwv))
    this%cref = [( dc * l, l = -ngwv, ngwv )]

    ! Wavenumber and effective wavenumber come from the wavelength.
    kwv = 2.0*PI / wavelength
    this%effkwv = fcrit2 * kwv
    this%is_set = .true.
  end subroutine create_new_GWband

  !=============================================
   subroutine gw_common_init(prefint)
     real, dimension(:), intent(in) :: prefint !< reference phalf [Pa]
     integer :: nlev
     if (module_is_initialized) return

     nlev = size(prefint(:))-1
     allocate( alpha(nlev) )
     ! interpolation to get profiles of Newtownian cooling coefficients
     call calculate_alpha( nlev, prefint, alpha )

     module_is_initialized = .true.

   end subroutine gw_common_init

  !=============================================
   subroutine gw_common_end
   ! deallocate module variables and reset model flag
      deallocate (alpha)
      module_is_initialized = .false.
   end subroutine gw_common_end

  !=============================================
    subroutine gw_drag_prof(ncol, pver, band, pint, &
     src_level, tend_level, dt, t,    &
     rhoi, ni, ubm,  ubi,  xv,    yv,   &
     tau,  utgw,  vtgw, &
     ttgw,  gwut, ro_adjust, tau_adjust, &
     kwvrdg, satfac_in )

  !-----------------------------------------------------------------------
  ! Solve for the drag profile given GW source tau
  ! Currently only used for gw_beres. Potentially can be used for other GWs.
  ! 1. scan up from the wave source to determine the stress profile
  ! 2. scan down the stress profile to determine the tendencies
  !     => apply bounds to the tendency
  !          a. from wkb solution
  !          b. from computational stability constraints
  !     => adjust stress on interface below to reflect actual bounded
  !        tendency
  !--------------------------------------------------------------------

  !------------------------------Arguments--------------------------------
  ! Column and vertical dimension.
  integer, intent(in) :: ncol,pver
  ! Wavelengths.
  type(GWBand), intent(in) :: band
  ! Pressure coordinates.
  ! Interface pressures.
  real, intent(in) :: pint(ncol,pver+1)
  ! Level from which gravity waves are propagated upward.
  integer, intent(in) :: src_level(ncol)
  ! Lowest level where wind tendencies are calculated.
  integer, intent(in) :: tend_level(ncol)
  ! Using tend_level > src_level allows the orographic waves to prescribe
  ! wave propagation up to a certain level, but then allow wind tendencies
  ! and adjustments to tau below that level.

  ! Time step.
  real, intent(in) :: dt

  ! Midpoint and interface temperatures.
  real, intent(in) :: t(ncol,pver)
  ! Interface densities.
  real, intent(in) :: rhoi(ncol,pver+1)
  ! Midpoint and interface Brunt-Vaisalla frequencies.
  real, intent(in) :: ni(ncol,pver+1)
  ! Projection of wind at midpoints and interfaces.
  real, intent(in) :: ubm(ncol,pver), ubi(ncol,pver+1)
  ! Unit vectors of source wind (zonal and meridional components).
  real, intent(in) :: xv(ncol), yv(ncol)

  ! Wave Reynolds stress.
  real , intent(inout) :: tau(ncol,-band%ngwv:band%ngwv,pver+1)
  ! Zonal/meridional wind tendencies.
  real, intent(out) :: utgw(ncol,pver), vtgw(ncol,pver)
  ! Gravity wave heating tendency.
  real, intent(out) :: ttgw(ncol,pver)

  ! Gravity wave wind tendency for each wave.
  real , intent(out) :: gwut(ncol,pver,-band%ngwv:band%ngwv)

  ! Adjustment parameter for IGWs.
  real, intent(in), optional :: &
       ro_adjust(ncol,-band%ngwv:band%ngwv,pver+1)

  ! Adjustment parameter for TAU.
  real, intent(in), optional :: &
       tau_adjust(ncol,pver+1)

  ! Diagnosed horizontal wavenumber for ridges.
  real, intent(in), optional :: &
       kwvrdg(ncol)

  ! Factor for saturation calculation. Here backwards
  ! compatibility. I believe it should be 1.0 (jtb).
  ! Looks like it has been 2.0 for a while in CAM.
  real, intent(in), optional :: &
       satfac_in

  !---------------------------Local storage-------------------------------

  ! Level, wavenumber, constituent and column loop indices.
  integer :: k, l, m, i, kmax

  ! Lowest tendency and source levels.
  integer :: kbot_tend, kbot_src

  ! Wave phase speeds for each column.
  real   :: cref(-band%ngwv:band%ngwv)

  ! Imaginary part of vertical wavenumber.
  real  :: mi(ncol)
  ! Stress after damping.
  real  :: taudmp(ncol)
  ! Saturation stress.
  real  :: tausat(ncol)
  ! (ub-c) and (ub-c)**2
  real  :: ubmc(ncol), ubmc2(ncol)
  ! Temporary ubar tendencies (overall, and at wave l).
  real  :: ubt(ncol,pver), ubtl(ncol)
  real  :: wrk(ncol)
  ! Temporary effkwv
  real  :: effkwv(ncol)

    ! Delta Interface pressures.
  real :: delp(ncol,pver)
    ! Log of interface pressures.
  real :: piln(ncol,pver+1)

  ! saturation factor. Defaults to 2.0
  ! unless overidden by satfac_in
  real  :: satfac

  real  :: near_zero = tiny(1.0)

! pressure level related variables
  delp=pint(:,2:pver+1)-pint(:,1:pver) !< pressure layer thickness
  piln=log(pint) !< log of phalf

  if (present(satfac_in)) then
     satfac = satfac_in
  else
     satfac = 2.0
  endif

  ! read phase speed from "band"
  cref=band%cref

  ! Lowest levels that loops need to iterate over.
  kbot_tend = maxval(tend_level)
  kbot_src = maxval(src_level)

  ! Initialize gravity wave drag tendencies to zero.

  utgw    = 0.0
  vtgw    = 0.0
  gwut    = 0.0
  ttgw    = 0.0

  ! Workaround floating point exception issues on Intel by initializing
  ! everything that's first set in a where block.
  mi     = 0.0
  taudmp = 0.0
  tausat = 0.0
  ubmc   = 0.0
  ubmc2  = 0.0
  wrk    = 0.0

  if (present(kwvrdg)) then
    effkwv = kwvrdg
  else
    effkwv = band%effkwv
  endif

  !------------------------------------------------------------------------
  ! Compute the stress profiles and diffusivities
  !------------------------------------------------------------------------

  ! Loop from bottom to top to get stress profiles.
! !$OMP parallel do default(none) shared(kbot_src,ktop,kvtt,band,ubi,c,effkwv,rhoi,ni,satfac, &
! !$OMP                                  ro_adjust,ncol,alpha,piln,t,rog,src_level,tau_adjust,tau) &
! !$OMP                          private(k,d,l,i,tausat,taudmp,ubmc,ubmc2,wrk,mi)
  do k = kbot_src, ktop, -1  !++ but this is in model now

     do l = -band%ngwv, band%ngwv

        do i=1,ncol

        tausat(i) = 0.0
        taudmp(i) = 0.0

        if (src_level(i) >= k) then

          ! Determine the absolute value of the saturation stress.
          ! Define critical levels where the sign of (u-c) changes between
          ! interfaces.
          ubmc(i) = ubi(i,k) - cref(l)

          ! Test to see if u-c has the same sign here as the level below.
          if (ubmc(i) > 0.0 .eqv. ubi(i,k+1) > cref(l)) then
             if (ni(i,k) /= 0.0) &
                 tausat(i) = abs( effkwv(i) * rhoi(i,k) * ubmc(i)**3 / &
                                  (satfac*ni(i,k)) )
             if (present(ro_adjust)) &
                 tausat(i) = tausat(i) * sqrt(ro_adjust(i,l,k))
             if (present(tau_adjust)) &
                 tausat(i) = tausat(i) * tau_adjust(i,k)
          endif

          ! Compute stress for each wave. The stress at this level is the
          ! min of the saturation stress and the stress at the level below
          ! reduced by damping. The sign of the stress must be the same as
          ! at the level below.
          ubmc2(i) = max(ubmc(i)**2, ubmc2mn)
          mi(i) = ni(i,k) / (2.0 * effkwv(i) * ubmc2(i)) * &  ! Is this 2.0 related to satfac?
                 (alpha(k) + ni(i,k)**2/ubmc2(i) * dback)
          wrk(i) = -2.0*mi(i)*t(i,k)*(piln(i,k+1) - piln(i,k))*RDGAS/GRAV
          wrk(i) = max( wrk(i), -200.0 ) * exp(wrk(i))
          taudmp(i) = tau(i,l,k+1) * exp(wrk(i))
          ! For some reason, PGI 14.1 loses bit-for-bit reproducibility if
          ! we limit tau, so instead limit the arrays used to set it.
          if (tausat(i) <= taumin) tausat(i) = 0.0
          if (taudmp(i) <= taumin) taudmp(i) = 0.0
          tau(i,l,k) = min(taudmp(i), tausat(i))

        endif
        end do
     end do
  end do

  ! Force tau at the top of the model to zero, if requested.
  if (tau_0_ubc) tau(:,:,ktop) = 0.0

  !------------------------------------------------------------------------
  ! Compute the tendencies from the stress divergence.
  !------------------------------------------------------------------------


  ! Loop over levels from top to bottom
! !$OMP parallel do default(none) shared(kbot_tend,ktop,band,ncol,tau,delp,rdelp,c,ubm,dt,gravit,utgw,vtgw, &
! !$OMP                                  gwut,ubt,xv,yv,tend_level,near_zero) &
! !$OMP                          private(k,l,i,ubtl)
  do k = ktop, kbot_tend

     ! Accumulate the mean wind tendency over wavenumber.
     ubt(:,k) = 0.0

     do l = -band%ngwv, band%ngwv    ! loop over wave

        do i=1,ncol

          ! Determine the wind tendency, including excess stress carried down
          ! from above.
          ubtl(i) = GRAV * (tau(i,l,k+1)-tau(i,l,k)) / delp(i,k)


          ! Apply first tendency limit to maintain numerical stability.
          ! Enforce du/dt < |c-u|/dt  so u-c cannot change sign
          !    (u^n+1 = u^n + du/dt * dt)
          ! The limiter is somewhat stricter, so that we don't come anywhere
          ! near reversing c-u.
          ubtl(i) = min(ubtl(i), umcfac * abs(cref(l)-ubm(i,k)) / dt)

          if (k <= tend_level(i)) then

           ! Save tendency for each wave (for later computation of kzz).
           ! sign function returns magnitude of ubtl with sign of c-ubm
           ! Renders ubt/ubm check for mountain waves unecessary
           gwut(i,k,l) = sign(ubtl(i), cref(l)-ubm(i,k))
           ubt(i,k) = ubt(i,k) + gwut(i,k,l)

          end if

        end do

     end do

     do l = -band%ngwv, band%ngwv
        ! Redetermine the effective stress on the interface below from the
        ! wind tendency. If the wind tendency was limited above, then the
        ! new stress will be smaller than the old stress, causing stress
        ! divergence in the next layer down. This smoothes large stress
        ! divergences downward while conserving total stress.
        ! Include a protection on SMALL gwut to prevent floating point
        ! issues.
        !--------------------------------------------------
        do i=1,ncol
         if ( abs(gwut(i,k,l)) < near_zero ) then
           gwut(i,k,l) = 0.0
         end if
         if (k <= tend_level(i)) then
           tau(i,l,k+1) = tau(i,l,k) + &
                abs(gwut(i,k,l)) * delp(i,k) / GRAV
                !!! abs(gwut(i,k,l)) * p%del(i,k) / GRAV
         end if
        end do
     end do

     ! Project the mean wind tendency onto the components.
     do i=1,ncol
       if (k <= tend_level(i)) then
        utgw(i,k) = ubt(i,k) * xv(i)
        vtgw(i,k) = ubt(i,k) * yv(i)
       end if
     end do

     ! End of level loop.
  end do

  ! Evaluate second temperature tendency term: Conversion of kinetic
  ! energy into thermal.
! !$OMP parallel do default(none) shared(kbot_tend,ktop,band,ttgw,ubm,c,gwut) &
! !$OMP                          private(k,l)
  do k = ktop, kbot_tend
     do l = -band%ngwv, band%ngwv
        ttgw(:,k) = ttgw(:,k) - (ubm(:,k) - cref(l)) * gwut(:,k,l)
     end do
  end do
  ttgw = ttgw / CP_AIR

  end subroutine gw_drag_prof

  !========================================================
  subroutine energy_momentum_adjust(ncol, pver, band, pint,  u, v, dt, tau, &
                               effgw, t, ubm, ubi, xv, yv, utgw, vtgw, ttgw, &
                               tend_level, tndmax_in)
!-------------------Arguments-------------------------------------
  integer, intent(in) :: ncol, pver
  ! Wavelengths.
  type(GWBand), intent(in) :: band
  ! Pressure interfaces.
  real, intent(in) :: pint(ncol,pver+1)
  ! Winds
  real, intent(in) :: u(ncol,pver)      ! Midpoint zonal winds. ( m s-1)
  real, intent(in) :: v(ncol,pver)      ! Midpoint meridional winds. ( m s-1)
  ! timestep
  real, intent(in) :: dt
  ! Wave Reynolds stress.
  real , intent(in) :: tau(ncol,-band%ngwv:band%ngwv,pver+1)
  ! Efficiency
  real, intent(in) :: effgw(ncol)
  ! Temperature and winds
  real, intent(in) :: t(ncol,pver), ubm(ncol,pver), ubi(ncol,pver)
  ! projected winds.
  real, intent(in) :: xv(ncol), yv(ncol)
  ! tendency level index index
  integer, intent(in) :: tend_level(ncol)
  ! Tendency limiter
  real, intent(in), optional :: tndmax_in
  ! Tendencies.
  real, intent(inout) :: utgw(ncol,pver)
  real, intent(inout) :: vtgw(ncol,pver)
  real, intent(inout) :: ttgw(ncol,pver)
!--------------------------local variables------------------------------
  real :: taucd(ncol,pver+1,4)
  real :: um_flux(ncol), vm_flux(ncol)
  real :: de(ncol)

    ! Pressure thickness.
  real :: delp(ncol,pver)

  real :: tndmax,utfac,uhtmax
  integer :: kbot(ncol)

  ! Level index.
  integer :: i,k,l

! Pressure thickness
  delp=pint(:,2:pver+1)-pint(:,1:pver)
! Maximum wind tendency from stress divergence (before efficiency applied).
  if (present(tndmax_in)) then
     tndmax = tndmax_in
  else
     tndmax = 400.0 / SECONDS_PER_DAY
  endif

  do i=1,ncol
   !---------------------------------------------------------------------------------------
   ! Apply efficiency factor and tendency limiter to prevent unrealistically strong forcing
   !---------------------------------------------------------------------------------------
    uhtmax = 0.0
    utfac  = 1.0
    do k = ktop, pver
       uhtmax = max(sqrt(utgw(i,k)**2 + vtgw(i,k)**2), uhtmax)
    end do
    if (uhtmax > tndmax) utfac = tndmax/uhtmax
    do k = ktop, pver
       utgw(i,k) = utgw(i,k)*effgw(i)*utfac
       vtgw(i,k) = vtgw(i,k)*effgw(i)*utfac
       ttgw(i,k) = ttgw(i,k)*effgw(i)*utfac
    end do
  end do  ! i=1,ncol

  if (band%ngwv /= 0) then
   ! compute momentum and energy flux changes
    taucd = calc_taucd(ncol, pver, band%ngwv, tend_level, tau, band%cref, xv, yv, ubi)
    call momentum_flux(tend_level, taucd, um_flux, vm_flux)
    call energy_change(ncol,pver, dt, delp, u, v, utgw, vtgw, ttgw, de)
   ! add sub-source fixers to tendencies
    call momentum_fixer(ncol, pver, tend_level, pint, um_flux, vm_flux, utgw, vtgw)
    call energy_fixer(ncol, pver, tend_level, pint, de, ttgw)
  endif
end subroutine energy_momentum_adjust

 !============================================================
    subroutine calculate_alpha( pver, phalf, alpha)
    integer, intent(in):: pver !< dimension of the vertical coordinate
    real, dimension(pver+1), intent(in) :: phalf !< reference pressure at half levels [Pa]
    real, dimension(pver+1), intent(out) :: alpha !<Newtonian cooling coefficient needed in gw_drag_prof
    ! parivate data
    integer :: k,i
    ! Levels of pre-calculated Newtonian cooling (1/day).
    ! The following profile is digitized from:
    ! Wehrbein and Leovy (JAS, 39, 1532-1544, 1982) figure 5
     integer, parameter :: nalph = 71
       real :: alpha0(nalph) = [ &
       0.1,         0.1,         0.1,         0.1,         &
       0.1,         0.1,         0.1,         0.1,         &
       0.1,         0.1,         0.10133333,  0.104,       &
       0.108,       0.112,       0.116,       0.12066667,  &
       0.126,       0.132,       0.138,       0.144,       &
       0.15133333,  0.16,        0.17,        0.18,        &
       0.19,        0.19933333,  0.208,       0.216,       &
       0.224,       0.232,       0.23466667,  0.232,       &
       0.224,       0.216,       0.208,       0.20133333,  &
       0.196,       0.192,       0.188,       0.184,       &
       0.18266667,  0.184,       0.188,       0.192,       &
       0.196,       0.19333333,  0.184,       0.168,       &
       0.152,       0.136,       0.12133333,  0.108,       &
       0.096,       0.084,       0.072,       0.061,       &
       0.051,       0.042,       0.033,       0.024,       &
       0.017666667, 0.014,       0.013,       0.012,       &
       0.011,       0.010333333, 0.01,        0.01,        &
       0.01,        0.01,        0.01                         &
       ]

   ! Pressure levels that were used to calculate alpha0 (hPa).
    real :: palph(nalph) = [ &
       2.06115E-06, 2.74280E-06, 3.64988E-06, 4.85694E-06, &
       6.46319E-06, 8.60065E-06, 1.14450E-05, 1.52300E-05, &
       2.02667E-05, 2.69692E-05, 3.58882E-05, 4.77568E-05, &
       6.35507E-05, 8.45676E-05, 0.000112535, 0.000149752, &
       0.000199277, 0.000265180, 0.000352878, 0.000469579, &
       0.000624875, 0.000831529, 0.00110653,  0.00147247,  &
       0.00195943,  0.00260744,  0.00346975,  0.00461724,  &
       0.00614421,  0.00817618,  0.0108801,   0.0144783,   &
       0.0192665,   0.0256382,   0.0341170,   0.0453999,   &
       0.0604142,   0.0803939,   0.106981,    0.142361,    &
       0.189442,    0.252093,    0.335463,    0.446404,    &
       0.594036,    0.790490,    1.05192,     1.39980,     &
       1.86273,     2.47875,     3.29851,     4.38936,     &
       5.84098,     7.77266,     10.3432,     13.7638,     &
       18.3156,     24.3728,     32.4332,     43.1593,     &
       57.4326,     76.4263,     101.701,     135.335,     &
       180.092,     239.651,     318.907,     424.373,     &
       564.718,     751.477,     1200.                        &
       ]

  ! pre-calculated newtonian damping:
  !     * convert to 1/s
  !     * ensure it is not smaller than 1e-6
  !     * convert palph from hpa to pa

  do k=1,nalph
     alpha0(k) = alpha0(k) / SECONDS_PER_DAY
     alpha0(k) = max(alpha0(k), 1.e-6)
     palph(k) = palph(k)*1.e2
  end do

    call interp_1d(palph, phalf, alpha0, alpha,"linear")
  end subroutine calculate_alpha

  !===================================================
  subroutine gw_prof (ncol, pver, pint, pmid , t, rhoi, ni)
  !-----------------------------------------------------------------------
  ! Compute profiles of background state quantities for the multiple
  ! gravity wave drag parameterization.
  !
  ! The parameterization is assumed to operate only where water vapor
  ! concentrations are negligible in determining the density.
  !-----------------------------------------------------------------------
  !------------------------------Arguments--------------------------------
  ! Column and vertical dimensions.
  integer, intent(in) :: ncol,pver
  ! Pressure coordinates.
  real, intent(in) :: pmid(ncol,pver)
  real, intent(in) :: pint(ncol,pver+1)

  ! Midpoint temperatures.
  real, intent(in) :: t(ncol,pver)

  ! Interface density.
  real, intent(out) :: rhoi(ncol,pver+1)
  ! Midpoint and interface Brunt-Vaisalla frequencies.
  real, intent(out) :: ni(ncol,pver+1)

  !---------------------------Local Storage-------------------------------
  ! Column and level indices.
  integer :: i,k

  ! dt/dp
  real :: dtdp
  ! Brunt-Vaisalla frequency squared.
  real :: n2

  ! Interface temperature.
  real :: ti(ncol,pver+1)

  ! Minimum value of Brunt-Vaisalla frequency squared.
  real, parameter :: n2min = 5.e-5

  !------------------------------------------------------------------------
  ! Determine the interface densities and Brunt-Vaisala frequencies.
  !------------------------------------------------------------------------

  ! The top interface values are calculated assuming an isothermal
  ! atmosphere above the top level.
  k = 1
  do i = 1, ncol
     ti(i,k) = t(i,k)
     rhoi(i,k) = pint(i,k) / (RDGAS*ti(i,k))
     ni(i,k) = sqrt(GRAV*GRAV / (CP_AIR*ti(i,k)))
  end do

  ! Interior points use centered differences.
  ti(:,2:pver) = (t(:,1:pver-1)+t(:,2:pver))/2.0
  do k = 2, pver
     do i = 1, ncol
        rhoi(i,k) = pint(i,k) / (RDGAS*ti(i,k))
        dtdp = (t(i,k)-t(i,k-1)) / ( pmid(i,k)-pmid(i,k-1) ) ! * p%rdst(i,k-1)
        n2 = GRAV*GRAV/ti(i,k) * (1.0/CP_AIR - rhoi(i,k)*dtdp)
        ni(i,k) = sqrt(max(n2min, n2))
     end do
  end do

  ! Bottom interface uses bottom level temperature, density; next interface
  ! B-V frequency.
  k = pver+1
  do i = 1, ncol
     ti(i,k) = t(i,k-1)
     rhoi(i,k) = pint(i,k) / (RDGAS*ti(i,k))
     ni(i,k) = ni(i,k-1)
  end do

  !------------------------------------------------------------------------
  ! Determine the midpoint Brunt-Vaisala frequencies.
  !------------------------------------------------------------------------
  !nm=(ni(:,1:pver)+ni(2:pver+1))/2.0

end subroutine gw_prof


!==========================================================================

! Calculate the amount of momentum conveyed from below the gravity wave
! region, to the region where gravity waves are calculated.
subroutine momentum_flux(tend_level, taucd, um_flux, vm_flux)

  ! Bottom stress level.
  integer, intent(in) :: tend_level(:)
  ! Projected stresses.
  real, intent(in) :: taucd(:,:,:)
  ! Components of momentum change sourced from the bottom.
  real, intent(out) :: um_flux(:), vm_flux(:)

  integer :: i

  ! Tendency for U & V below source level.
  do i = 1, size(tend_level)
     um_flux(i) = taucd(i,tend_level(i)+1, east) + &
                  taucd(i,tend_level(i)+1, west)
     vm_flux(i) = taucd(i,tend_level(i)+1,north) + &
                  taucd(i,tend_level(i)+1,south)
  end do

end subroutine momentum_flux


!==========================================================================

! Subtracts a change in momentum in the gravity wave levels from wind
! tendencies in lower levels, ensuring momentum conservation.
subroutine momentum_fixer(ncol, pver, tend_level, phalf, um_flux, vm_flux, utgw, vtgw)

  integer, intent(in) :: ncol, pver
  ! Bottom stress level.
  integer, intent(in) :: tend_level(ncol)
  ! Pressure coordinates.
  real, intent(in) :: phalf(ncol,pver+1)
  ! Components of momentum change sourced from the bottom.
  real, intent(in) :: um_flux(ncol), vm_flux(ncol)
  ! Wind tendencies.
  real, intent(inout) :: utgw(ncol,pver), vtgw(ncol,pver)

  ! Indices.
  integer :: i, k
  ! Reciprocal of total mass.
  real :: rdm(ncol)

  ! Total mass from ground to source level: rho*dz = dp/gravit
  do i = 1, ncol
     rdm(i) = GRAV/(phalf(i,pver+1)-phalf(i,tend_level(i)+1))
  end do

  do k = minval(tend_level)+1, pver
     where (k > tend_level)
        utgw(:,k) = utgw(:,k) + -um_flux*rdm
        vtgw(:,k) = vtgw(:,k) + -vm_flux*rdm
     end where
  end do

end subroutine momentum_fixer

!==========================================================================

! Calculate the change in total energy from tendencies up to this point.
subroutine energy_change(ncol,pver, dt, pdel, u, v, dudt, dvdt, dtdt, de)

  integer, intent(in) :: ncol, pver
  ! Time step.
  real, intent(in) :: dt
  ! Pressure coordinates interval.
  real, intent(in) :: pdel(ncol,pver)
  ! Winds at start of time step.
  real, intent(in) :: u(ncol,pver), v(ncol,pver)
  ! Wind tendencies.
  real, intent(in) :: dudt(ncol,pver), dvdt(ncol,pver)
  ! Heating tendency.
  real, intent(in) :: dtdt(ncol,pver)
  ! Change in energy.
  real, intent(out) :: de(ncol)

  ! Level index.
  integer :: k

  ! Net gain/loss of total energy in the column.
  de = 0.0
  do k = 1, pver
     de = de + pdel(:,k)/GRAV * (CP_AIR*dtdt(:,k) + &
          dudt(:,k)*(u(:,k)+dudt(:,k)*0.5*dt) + &
          dvdt(:,k)*(v(:,k)+dvdt(:,k)*0.5*dt) )
  end do

end subroutine energy_change

!==========================================================================

! Subtract change in energy from the heating tendency in the levels below
! the gravity wave region.
subroutine energy_fixer(ncol, pver, tend_level, pint, de, ttgw)

  integer, intent(in) :: ncol, pver
  ! Bottom stress level.
  integer, intent(in) :: tend_level(ncol)
  ! Pressure coordinates at half levels.
  real, intent(in) :: pint(ncol,pver+1)
  ! Change in energy.
  real, intent(in) :: de(ncol)
  ! Heating tendency.
  real, intent(inout) :: ttgw(ncol,pver)

  ! Column/level indices.
  integer :: i, k
  ! Energy change to apply divided by all the mass it is spread across.
  real :: de_dm(ncol)

  do i = 1, ncol
     de_dm(i) = -de(i)*GRAV/(pint(i,pver+1)-pint(i,tend_level(i)+1))
  end do

  ! Subtract net gain/loss of total energy below tend_level.
  do k = minval(tend_level)+1, pver
     where (k > tend_level)
        ttgw(:,k) = ttgw(:,k) + de_dm/CP_AIR
     end where
  end do

end subroutine energy_fixer

!==========================================================================

! Calculate Reynolds stress for waves propagating in each cardinal
! direction.

function calc_taucd(ncol, pver, ngwv, tend_level, tau, c, xv, yv, ubi) &
     result(taucd)

  ! Column and gravity wave wavenumber dimensions.
  integer, intent(in) :: ncol, pver, ngwv
  ! Lowest level where wind tendencies are calculated.
  integer, intent(in) :: tend_level(:)
  ! Wave Reynolds stress.
  real, intent(in) :: tau(:,-ngwv:,:)
  ! Wave phase speeds, same for all columns
  real, intent(in) :: c(-ngwv:)
  ! Unit vectors of source wind (zonal and meridional components).
  real, intent(in) :: xv(:), yv(:)
  ! Projection of wind at interfaces.
  real, intent(in) :: ubi(:,:)

  real :: taucd(ncol,pver+1,4)

  ! Indices.
  integer :: i, k, l

  ! ubi at tend_level.
  real :: ubi_tend(ncol)

  ! Signed wave Reynolds stress.
  real :: tausg(ncol)

  ! Reynolds stress for waves propagating behind and forward of the wind.
  real :: taub(ncol)
  real :: tauf(ncol)

  taucd = 0.
  tausg = 0.

  ubi_tend = (/ (ubi(i,tend_level(i)+1), i = 1, ncol) /)

  do k = ktop, maxval(tend_level)+1

     taub = 0.
     tauf = 0.

     do l = -ngwv, ngwv
        where (k-1 <= tend_level)

           tausg = sign(tau(:,l,k), c(l)-ubi(:,k))

           where ( c(l) < ubi_tend )
              taub = taub + tausg
           elsewhere
              tauf = tauf + tausg
           end where

        end where
     end do

     where (k-1 <= tend_level)
        where (xv > 0.)
           taucd(:,k,east) = tauf * xv
           taucd(:,k,west) = taub * xv
        elsewhere
           taucd(:,k,east) = taub * xv
           taucd(:,k,west) = tauf * xv
        end where

        where ( yv > 0.)
           taucd(:,k,north) = tauf * yv
           taucd(:,k,south) = taub * yv
        elsewhere
           taucd(:,k,north) = taub * yv
           taucd(:,k,south) = tauf * yv
        end where
     end where

  end do

end function calc_taucd
end module gw_common_utils
