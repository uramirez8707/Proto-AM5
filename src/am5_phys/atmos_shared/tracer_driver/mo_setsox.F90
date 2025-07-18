module MO_SETSOX_MOD

  use mo_chem_utls_mod, only : get_spc_ndx
  use tropchem_types_mod, only : tropchem_opt, tropchem_diag, small_value, missing_value
  use fms_mod, only :  mpp_pe, mpp_root_pe, FATAL, error_mesg
  use cloud_chem, only : cloud_pH, cloud_so2_chem, cloud_nb_diag 
  use field_manager_mod,  only : MODEL_ATMOS
  use tracer_manager_mod, only : get_tracer_index        
  use tracer_manager_mod,    only: get_tracer_index,  query_method
  use field_manager_mod,     only: parse

  use cloud_chem, only : CLOUD_CHEM_F1P
  use aerosol_thermodynamics,   only : aerosol_thermo, AERO_LEGACY, AERO_ISORROPIA, NO_AERO
  use mpp_mod,            only : mpp_clock_id,         &
                                 mpp_clock_begin,      &
                                 mpp_clock_end

  implicit none

  character(len=128), parameter :: version     = '$Id: mo_setsox.F90,v 19.0 2012/01/06 20:34:08 fms Exp $'
  character(len=128), parameter :: tagname     = '$Name: siena_201305 $'
  logical                       :: module_is_initialized = .false.

  !<f1p
  integer    ::      ox_ndx, hno3_ndx, h2o2_ndx, so2_ndx, so4_ndx, nh3_ndx, nh4no3_ndx
  integer    ::      ho2_ndx, nh4_ndx, co2_ndx, hcooh_ndx, ch3cooh_ndx, n2o5_ndx
  real       ::      frac_ic_so4, frac_ic_no3, frac_ic_nh4
  real       ::      frac_ic_so4_snow, frac_ic_no3_snow, frac_ic_nh4_snow

  logical    ::      nh4no3_is_no3
  integer    ::      isoropia_clock_id

  public     ::      setsox, setsox_init
  private
  !>

CONTAINS

  subroutine setsox( press, plonl, dtime, tfld, qfld, &
       cwat, frac_liq, cldfr,xhnm,                    &
       qin, xco2, trop_diag_array, trop_option, trop_diag)

    !-----------------------------------------------------------------------      
    !          ... Compute heterogeneous reactions of SOX
    !
    !       (0) using initial PH to calculate PH
    !           (a) HENRY's law constants
    !           (b) PARTIONING
    !           (c) PH values
    !
    !       (1) using new PH to repeat
    !           (a) HENRY's law constants
    !           (b) PARTIONING
    !           (c) REACTION rates
    !           (d) PREDICTION
    !-----------------------------------------------------------------------      


    implicit none
    !-----------------------------------------------------------------------      
    !      ... Dummy arguments
    !-----------------------------------------------------------------------      
    integer, intent(in)  ::    plonl               ! number of local longitude points
    real, intent(in)     ::    dtime               ! time step (sec)
    real, intent(inout)  ::    qin(:,:,:)          ! xported species ( vmr )
    real, intent(inout)  ::    trop_diag_array(:,:,:)
    real, intent(in)     ::    xco2(:,:), cldfr(:,:) !co2 (vmr), cloud fraction
    real, intent(in)     ::    xhnm(:,:)           ! total atms density ( /cm**3)
    real, dimension(:,:), intent(in) ::  &
         tfld, &               ! temperature
         qfld, &               ! specific humidity( kg/kg )
         cwat, frac_liq,  &               ! cloud liquid water content (kg/kg)
         press                 ! midpoint pressure ( Pa )
    type(tropchem_opt),  intent(in) :: trop_option
    type(tropchem_diag), intent(in) :: trop_diag


    !-----------------------------------------------------------------------      
    !      ... Local variables
    !
    !           xhno3 ... in mixing ratio
    !-----------------------------------------------------------------------      
    integer, parameter :: itermax = 20
    real, parameter ::  const0 = 1.e3/6.023e23
    real, parameter ::  kh0 = 9.e3, &           ! HO2(g)          -> Ho2(a)
         kh1 = 2.05e-5, &        ! HO2(a)          -> H+ + O2-
         kh2 = 8.6e5,   &        ! HO2(a) + ho2(a) -> h2o2(a) + o2
         kh3 = 1.e8,    &        ! HO2(a) + o2-    -> h2o2(a) + o2
         Ra = 8314./101325., &   ! universal constant   (atm)/(M-K)
         xkw = 1.e-14            ! water acidity

    integer    ::      k, i, iter
    integer    ::      plev
    real       ::      wrk, delta
    real       ::      xk, xe, x2
    real       ::      tz, xl, px, qz, pz, es, qs, patm
    real       ::      Eso2, Eso4, Ehno3, Eco2, Eh2o, Enh3
    real       ::      hno3g, nh3g, so2g, h2o2g, o3g
    real       ::      rah2o2, rao3, ccc
    real       ::      cnh3, chno3, com, com1, com2, xra
    real       ::      RH
    !<f1p
    real, parameter ::  xH0 = 1.e-5
    real, parameter ::  Mw_air_over_Na = 1.e6*1.38e-23/287.*1.e-3                     ! /m3(a) Kg(a)/m3(a)  Kg(a)/L(a)  !28.97/6.0221413e23  
    integer    ::      n
    real       ::      thno3, tnh3,  rso2_o3,rso2_h2o2, pso4_o3,pso4_h2o2, pso4
    real       ::      exp_factor, EF, so2_after, ratio
    real       ::      correction

    real       :: xx0,yy1,xkp

    !-----------------------------------------------------------------------      
    !            for Ho2(g) -> H2o2(a) formation 
    !            schwartz JGR, 1984, 11589
    !-----------------------------------------------------------------------      
    real       ::      kh4    ! kh2+kh3
    real       ::      xam    ! air density /cm3
    real       ::      ho2s   ! ho2s = ho2(a)+o2-
    real       ::      r1h2o2 ! prod(h2o2) by ho2 in mole/L(w)/s
    real       ::      r2h2o2 ! prod(h2o2) by ho2 in mix/s



    real, dimension(SIZE(tfld,1),SIZE(tfld,2))  :: &
         xhno3, xh2o2, xso2, xso4,&
         xnh3, xo3,               &
         xlwc, cfact,             &
         xant,                    &
         hehno3,                  &            ! henry law const for hno3
         heh2o2,                  &            ! henry law const for h2o2
         heso2,                   &            ! henry law const for so2
         henh3,                   &            ! henry law const for nh3
         xnh4,xalk,               &
         heo3, xho2,              &            ! henry law const for nh3
         xH, xhcooh, xch3cooh

    real, dimension(size(tfld,1),size(tfld,2))       :: frac_ic_no3_eff, frac_ic_nh4_eff, frac_ic_so4_eff
    real, dimension(plonl)  :: t_fac
    real    :: ediag(cloud_nb_diag)
    real    :: xso4_tmp, xnh4_tmp, xant_tmp, xnh3_tmp, xhno3_tmp
    real    :: xpH
    logical :: converged

    real :: frac_acid_ic !fraction of acid available for cloud chemistry

    if ( .not. module_is_initialized ) then
       call error_mesg ('setsox','setsox_init must be called first.', FATAL)
    end if

    plev = SIZE(tfld,2)
    xalk = 0.


    !-----------------------------------------------------------------
    !       ... NOTE: The press array is in pascals and must be
    !                 mutiplied by 10 to yield dynes/cm**2.
    !-----------------------------------------------------------------

    !==================================================================
    !       ... First set the PH
    !==================================================================
    !      ... Initial values
    !           The values of so2, so4 are after (1) SLT, and CHEM
    !-----------------------------------------------------------------

    !do k = 1,plev
    !        precip(:,k) = cmfdqr(:,k) + rain(:,k) - evapr(:,k)
    !end do

    do k = 1,plev
       cfact(1:,k) = xhnm(1:,k)             &
                            * 1.e6          &             ! /m3(a)
                            * 1.38e-23/287. &             ! Kg(a)/m3(a)
                            * 1.e-3                       ! Kg(a)/L(a)
    end do

    frac_ic_so4_eff   = 0.
    frac_ic_no3_eff   = 0.
    frac_ic_nh4_eff   = 0.

    do k = 1,plev
       xH(:,k) = xH0                                    ! initial PH value
       if ( trop_option%cloud_chem .eq. CLOUD_CHEM_F1P ) then
         xlwc(:,k) = frac_liq(:,k) * cwat(:,k) *cfact(:,k)           ! cloud water  L(water)/L(air)
       end if
       if( hno3_ndx > 0 ) then
          xhno3(:,k) = qin(:,k,hno3_ndx)               ! mixing ratio
       else
          xhno3(:,k) = 0.
       end if
       if( h2o2_ndx > 0 ) then
          xh2o2(:,k) = qin(:,k,h2o2_ndx)                  ! mixing ratio
       else
          xh2o2(:,k) = 0.
       end if
       if( so2_ndx > 0 ) then
          xso2(:,k) = qin(:,k,so2_ndx)                   ! mixing ratio
       else
          xso2(:,k) = 0.
       end if
       if( so4_ndx > 0 ) then
          xso4(:,k) = qin(:,k,so4_ndx)                   ! mixing ratio
       else
          xso4(:,k) = 0.
       end if
       if( nh3_ndx > 0 ) then
          xnh3(:,k) = qin(:,k,nh3_ndx)                   ! mixing ratio
       else
          xnh3(:,k) = 0.
       end if
       if( nh4no3_ndx > 0 ) then
          xant(:,k) = qin(:,k,nh4no3_ndx)                   ! mixing ratio
       else
          xant(:,k) = 0.
       end if
       if( nh4_ndx > 0 ) then
          xnh4(:,k) = qin(:,k,nh4_ndx)                   ! mixing ratio
       else
          xnh4(:,k) = 0.
       end if
       if( ox_ndx > 0 ) then
          xo3(:,k) = qin(:,k,ox_ndx)                    ! mixing ratio
       else
          xo3(:,k) = 0.
       end if
       if ( hcooh_ndx > 0 ) then
          xhcooh(:,k) = qin(:,k,hcooh_ndx)
       else
          xhcooh(:,k) = 0.
       end if
       if ( ch3cooh_ndx > 0 ) then
          xch3cooh(:,k) = qin(:,k,ch3cooh_ndx)
       else
          xch3cooh(:,k) = 0.
       end if
       if ( ho2_ndx > 0 ) then
          xho2(:,k) = qin(:,k,ho2_ndx)
       else
          xho2(:,k) = 0.
       end if

       do i = 1,plonl

          if (trop_option%modulate_frac_ic) then

          if ( tfld(i,k) .gt. 258. ) then
             !assume ice is from riming
             if ( frac_ic_so4 .gt. 0. ) frac_ic_so4_eff(i,k) = frac_ic_so4 * frac_liq(i,k)
             if ( frac_ic_no3 .gt. 0. ) frac_ic_no3_eff(i,k) = frac_ic_no3 * frac_liq(i,k)
             if ( frac_ic_nh4 .gt. 0. ) frac_ic_nh4_eff(i,k) = frac_ic_nh4 * frac_liq(i,k)                
          else
             !0 from Bergeron
             if ( frac_ic_so4 .gt. 0. ) frac_ic_so4_eff(i,k) = frac_ic_so4
             if ( frac_ic_no3 .gt. 0. ) frac_ic_no3_eff(i,k) = frac_ic_no3
             if ( frac_ic_nh4 .gt. 0. ) frac_ic_nh4_eff(i,k) = frac_ic_nh4
          end if

          else

             frac_ic_so4_eff(i,k) = frac_ic_so4
             frac_ic_no3_eff(i,k) = frac_ic_no3
             frac_ic_nh4_eff(i,k) = frac_ic_nh4

          end if

       end do

    end do

    if (trop_diag%ind_cloud_pH>0)  trop_diag_array(:,:,trop_diag%ind_cloud_pH)  = missing_value

    !-----------------------------------------------------------------
    !       ... Temperature dependent Henry constants
    !-----------------------------------------------------------------
    do k = 1,plev                                             !! plev loop for STEP 0
       do i = 1,plonl

          xl = xlwc(i,k) 

          qz = qfld(i,k)             ! H2O mass mxing ratio Kg/Kg
          pz = .01*press(i,k)        ! pressure in mb
          tz = tfld(i,k)
          wrk = tz - 273.
          es = 6.11*10.**(7.63*wrk/(241.9 + wrk))            ! Magnus EQ
          qs = .622*es/pz                                    ! sat mass mix (H2O)
          RH = 100.*qz/qs                                    ! relative huminity(%)
          RH = MIN( 100.,MAX( RH,0. ) )

          if( xl >= trop_option%min_lwc_for_cloud_chem ) then

             pso4_o3   = 0.
             pso4_h2o2 = 0.                
             
             if ( trop_option%cloud_chem .eq. CLOUD_CHEM_F1P ) then

                patm = pz/1013.25

                if (trop_option%modulate_frac_ic) then
                   frac_acid_ic = frac_liq(i,k)
                else
                   frac_acid_ic = 1.
                end if


                if ( cldfr(i,k) .gt. 1.e-10 .and. xso2(i,k) .gt. 0. ) then
                   xl = xl/cldfr(i,k)

                   if ( nh4no3_is_no3 ) then
                      tnh3 = xnh3(i,k) +  frac_ic_nh4_eff(i,k)*xnh4(i,k)
                   else
                      tnh3 = xnh3(i,k) &
                           + frac_ic_nh4_eff(i,k)*xnh4(i,k) &
                           + frac_ic_no3_eff(i,k)*xant(i,k)
                   end if

                   thno3= frac_acid_ic * xhno3(i,k) &                !for acids partioning between ice and liquid
                        + frac_ic_no3_eff(i,k)*xant(i,k)
                   

                   if (trop_option%cloud_H .lt. 0.) then                   
                      call cloud_pH(thno3, xso2(i,k), tnh3, xhcooh(i,k)*frac_acid_ic, xch3cooh(i,k)*frac_acid_ic, &
		                    xso4(i,k)*frac_ic_so4_eff(i,k), xco2(i,k), xalk(i,k), tfld(i,k), patm, xl, xhnm(i,k),&
				    trop_option%cloud_chem_ph_solver,xH(i,k), ediag)
                   else
                      xH(i,k)  = trop_option%cloud_H
                      ediag(:) = 0.
                   end if

                   call cloud_so2_chem(patm, xH(i,k), tfld(i,k), xl, rso2_h2o2, rso2_o3)

                   !amount of so4 formed
                   !SO2i*(1 - exp(kdt/a))/(1-SO2i/H2O2i exp(kdt/a)) with a = 1/SO20 - H2O20

                   exp_factor = rso2_h2o2 * (xso2(i,k) - xh2o2(i,k)) * dtime            

                   !this should have little impact
                   !this adds some safety when xso2~=xh2o2. There is also no reason to not calculate the exp when exp_factor<-600, although this will lead to the same results
                   if ( exp_factor .lt. 600. .and. abs(exp_factor) .gt. small_value ) then
                     EF          = exp( exp_factor )
                     pso4_h2o2   = max(xso2(i,k) * xh2o2(i,k) * ( 1. - EF ) / ( xh2o2(i,k) -  xso2(i,k) * EF ),0.)                        
                   elseif (abs(exp_factor) .le. small_value ) then
                     pso4_h2o2   = rso2_h2o2 * xso2(i,k)**2 * dtime / (1 + rso2_h2o2*dtime*xso2(i,k))
                   else
                     pso4_h2o2   = min(xso2(i,k),xh2o2(i,k))
                   end if
                   
                   !                   so2_after  = max(xso2(i,k) - pso4_h2o2,0.)
                   so2_after = xso2(i,k)

                   exp_factor = rso2_o3 * (so2_after - xo3(i,k)) * dtime            

                   if ( abs(exp_factor) .lt. 600. .and. abs(exp_factor) .gt. small_value ) then
                     EF          = exp( exp_factor )
                     pso4_o3     = max(so2_after * xo3(i,k) * ( 1. - EF ) / ( xo3(i,k) -  so2_after * EF ),0.)
                   elseif (abs(exp_factor) .le. small_value ) then
                     pso4_o3     = rso2_o3 * xso2(i,k)**2 * dtime / (1 + rso2_o3*dtime*xso2(i,k))
                   else
                     pso4_o3 = min(so2_after,xo3(i,k))
                   end if
                   !>        

                   if ( (pso4_o3 + pso4_h2o2) .gt. xso2(i,k) ) then
                      ratio     = max(min( xso2(i,k) / ( pso4_o3+pso4_h2o2) ,1.),0.)
                      pso4_o3   = pso4_o3   * ratio
                      pso4_h2o2 = pso4_h2o2 * ratio
                   end if

                   pso4_o3   = pso4_o3*cldfr(i,k)
                   pso4_h2o2 = pso4_h2o2*cldfr(i,k)

                   ccc = pso4_h2o2
                   ccc = MAX( MIN( ccc, xso2(i,k), xh2o2(i,k) ), 0. )

                   if (trop_diag%ind_pso4_h2o2>0) &
                        trop_diag_array(i,k,trop_diag%ind_pso4_h2o2)   = ccc/dtime

                   xso4(i,k)  = xso4(i,k)  + ccc
                   xh2o2(i,k) = MAX( xh2o2(i,k) - ccc, small_value )
                   xso2(i,k)  = MAX( xso2(i,k)  - ccc, small_value )

                   !O3
                   ccc = pso4_o3
                   ccc = MAX( MIN( ccc,xo3(i,k),xso2(i,k) ), 0. )

                   if (trop_diag%ind_pso4_o3>0) &
                        trop_diag_array(i,k,trop_diag%ind_pso4_o3)   = ccc/dtime

                   xso4(i,k)  = xso4(i,k)  + ccc
                   xso2(i,k)  = MAX( xso2(i,k)  - ccc, small_value )                 
                   xo3(i,k)   = MAX( xo3(i,k)   - ccc, small_value )                 

                   !<f1p diag
                   if (trop_diag%ind_cloud_pH>0)      trop_diag_array(i,k,trop_diag%ind_cloud_pH)  =-log10(max(xH(i,k),1.e-20))

                 end if

               end if


!              ccc = pso4_h2o2
!              ccc = MAX( MIN( ccc, xso2(i,k), xh2o2(i,k) ), 0. )

!              if (trop_diag%ind_pso4_h2o2>0) &
!                   trop_diag_array(i,k,trop_diag%ind_pso4_h2o2)   = ccc/dtime

!              xso4(i,k)  = xso4(i,k)  + ccc
!              xh2o2(i,k) = MAX( xh2o2(i,k) - ccc, small_value )
!              xso2(i,k)  = MAX( xso2(i,k)  - ccc, small_value )

!              !O3
!              ccc = pso4_o3
!              ccc = MAX( MIN( ccc,xo3(i,k),xso2(i,k) ), 0. )

!              if (trop_diag%ind_pso4_o3>0) &
!                   trop_diag_array(i,k,trop_diag%ind_pso4_o3)   = ccc/dtime

!              xso4(i,k)  = xso4(i,k)  + ccc
!              xso2(i,k)  = MAX( xso2(i,k)  - ccc, small_value )                 
!              xo3(i,k)   = MAX( xo3(i,k)   - ccc, small_value )                 

!              !<f1p diag
!              if (trop_diag%ind_pH>0)      trop_diag_array(i,k,trop_diag%ind_pH)     =-log10(xH(i,k))*xlwc(i,k)
!              if (trop_diag%ind_cH>0)      trop_diag_array(i,k,trop_diag%ind_cH)     = xH(i,k)*xlwc(i,k)
!              if (trop_diag%ind_lwc>0)     trop_diag_array(i,k,trop_diag%ind_lwc)    = xlwc(i,k)
!              if (trop_diag%ind_eno3>0)    trop_diag_array(i,k,trop_diag%ind_eno3)   = ediag(1) *xlwc(i,k)
!              if (trop_diag%ind_ehcoo>0)   trop_diag_array(i,k,trop_diag%ind_ehcoo)  = ediag(2) *xlwc(i,k)
!              if (trop_diag%ind_ech3coo>0) trop_diag_array(i,k,trop_diag%ind_ech3coo)= ediag(3) *xlwc(i,k)
!              if (trop_diag%ind_ehso3>0)   trop_diag_array(i,k,trop_diag%ind_ehso3)  = ediag(4) *xlwc(i,k)
!              if (trop_diag%ind_eso3>0)    trop_diag_array(i,k,trop_diag%ind_eso3)   = ediag(5) *xlwc(i,k)
!              if (trop_diag%ind_ehco3>0)   trop_diag_array(i,k,trop_diag%ind_ehco3)  = ediag(6) *xlwc(i,k)
!              if (trop_diag%ind_eco3>0)    trop_diag_array(i,k,trop_diag%ind_eco3)   = ediag(7) *xlwc(i,k)
!              if (trop_diag%ind_eso4>0)    trop_diag_array(i,k,trop_diag%ind_eso4)   = ediag(8) *xlwc(i,k)
!              if (trop_diag%ind_enh4>0)    trop_diag_array(i,k,trop_diag%ind_enh4)   = ediag(9) *xlwc(i,k)
!              if (trop_diag%ind_eoh>0)     trop_diag_array(i,k,trop_diag%ind_eoh)    = ediag(10)*xlwc(i,k)
!              if (trop_diag%ind_ealk>0)    trop_diag_array(i,k,trop_diag%ind_ealk)   = ediag(11)*xlwc(i,k)


             !>         
          end if

          call mpp_clock_begin(isoropia_clock_id)
          call aerosol_thermo( trop_option%aerosol_thermo, min(rh,trop_option%max_rh_aerosol), tz, press(i,k), xso4(i,k), xnh3(i,k), xnh4(i,k), xhno3(i,k), xant(i,k), xph)
          if (trop_diag%ind_aerosol_pH>0) &
               trop_diag_array(i,k,trop_diag%ind_aerosol_pH)   = xph
          call mpp_clock_end(isoropia_clock_id)

          if ( trop_option%limit_no3 .and. trop_option%aerosol_thermo .eq. AERO_ISORROPIA ) then
             if (  (xnh4(i,k)+xant(i,k) .gt. small_value) .and. xant(i,k) / ( xnh4(i,k) + xant(i,k)) .gt. 0.75 ) then
                !force xant to be no more than xnh4
                correction = max(xnh4(i,k) - 2*xso4(i,k),0.)
                xhno3(i,k) = xhno3(i,k) + xant(i,k) - correction
                xant(i,k)  = correction
             end if
          end if

          
          !-----------------------------------------------------------------
          !      ... Washout SO2, SO4 and NH3
          !-----------------------------------------------------------------
          xso4(i,k)  = MAX( xso4(i,k),  small_value )
          xant(i,k)  = MAX( xant(i,k),  small_value )
          !xnh4(i,k)  = MAX( xnh4(i,k),  small_value ) !<f1p>
          xnh3(i,k)  = MAX( xnh3(i,k),  small_value )
          xso2(i,k)  = MAX( xso2(i,k),  small_value )
          !xo3(i,k)   = MAX( xo3(i,k),   small_value ) !<f1p>
          !xhno3(i,k) = MAX( xhno3(i,k), small_value ) !<f1p>

       end do
    end do

    !==============================================================
    !       ... Update the mixing ratios
    !==============================================================
    do k = 1,plev
       if( so2_ndx > 0 ) then
          qin(:,k,so2_ndx) =  MAX( xso2(:,k), small_value )
       end if
       if( so4_ndx > 0 ) then
          qin(:,k,so4_ndx) =  MAX( xso4(:,k), small_value )
       end if
       if( h2o2_ndx > 0 ) then
          qin(:,k,h2o2_ndx) =  MAX( xh2o2(:,k), small_value ) 
       end if
       if( nh3_ndx > 0 ) then
          qin(:,k,nh3_ndx) =  MAX( xnh3(:,k), small_value )
       end if
       if( nh4no3_ndx > 0 ) then
          qin(:,k,nh4no3_ndx) =  MAX( xant(:,k), small_value )
       end if
       if (nh4_ndx > 0 .and. trop_option%aerosol_thermo .ne. AERO_LEGACY) then
          qin(:,k,nh4_ndx) =  MAX( xnh4(:,k), small_value ) 
       end if
       if( hno3_ndx > 0 ) then
          qin(:,k,hno3_ndx) =  MAX( xhno3(:,k), small_value )
       end if
       if( ox_ndx > 0 ) then
          qin(:,k,ox_ndx) =  MAX( xo3(:,k), small_value )
       end if
    end do

  end subroutine SETSOX

  !<f1p
  subroutine setsox_init(trop_option)

    integer :: iflag
    logical :: flag
    character (len=500) :: text_in_scheme, text_in_param
    integer :: n
    type(tropchem_opt), intent(in) :: trop_option
    real :: nit_scale_factor

    ox_ndx      = get_spc_ndx( 'O3' )
    hno3_ndx    = get_spc_ndx( 'HNO3' )
    h2o2_ndx    = get_spc_ndx( 'H2O2' )
    so2_ndx     = get_spc_ndx( 'SO2' )
    so4_ndx     = get_spc_ndx( 'SO4' )
    nh3_ndx     = get_spc_ndx( 'NH3' )
    hcooh_ndx   = get_spc_ndx( 'HCOOH' )
    ch3cooh_ndx = get_spc_ndx( 'CH3COOH' )
    ho2_ndx     = get_spc_ndx( 'HO2' )
    nh4_ndx     = get_spc_ndx( 'NH4' )
    n2o5_ndx    = get_spc_ndx( 'N2O5' )
    nh4no3_ndx  = get_spc_ndx( 'NH4NO3' )

    frac_ic_nh4     = 0.
    frac_ic_no3     = 0.
    frac_ic_so4     = 0.


    flag = query_method('wet_deposition',MODEL_ATMOS,&
         get_tracer_index(MODEL_ATMOS,'so4'), &
         text_in_scheme,text_in_param)
    iflag=parse(text_in_param,'frac_incloud',frac_ic_so4)

    iflag=parse(text_in_param,'frac_incloud_snow',frac_ic_so4_snow)
    if (iflag == 0) then
       frac_ic_so4_snow = frac_ic_so4
    end if

    flag = query_method ('wet_deposition',MODEL_ATMOS,&
         get_tracer_index(MODEL_ATMOS,'nh4'), &
         text_in_scheme,text_in_param)
    iflag=parse(text_in_param,'frac_incloud',frac_ic_nh4)
    iflag=parse(text_in_param,'frac_incloud_snow',frac_ic_nh4_snow)
    if (iflag == 0) then
       frac_ic_nh4_snow = frac_ic_nh4
    end if

    if ( nh4no3_ndx .gt. 0 ) then
       flag = query_method ('wet_deposition',MODEL_ATMOS,&
            get_tracer_index(MODEL_ATMOS,'nh4no3'), &
            text_in_scheme,text_in_param)
       iflag=parse(text_in_param,'frac_incloud',frac_ic_no3)
       iflag=parse(text_in_param,'frac_incloud_snow',frac_ic_no3_snow)
       if (iflag == 0) then
          frac_ic_no3_snow = frac_ic_no3
       end if
    else
       call error_mesg ('setsox','nh4no3 must be definedx.', FATAL)
    end if


    if ( trop_option%aerosol_thermo .eq. AERO_LEGACY .or. trop_option%aerosol_thermo .eq. NO_AERO ) then
       nh4no3_is_no3 = .false.
    elseif ( trop_option%aerosol_thermo .eq. AERO_ISORROPIA ) then
       nh4no3_is_no3 = .true.
    end if


    if ( trop_option%frac_aerosol_incloud .ge. 0. ) then
       !overwrite
       frac_ic_nh4 = trop_option%frac_aerosol_incloud
       frac_ic_no3 = trop_option%frac_aerosol_incloud
       frac_ic_so4 = trop_option%frac_aerosol_incloud
    end if


    if (mpp_root_pe() .eq. mpp_pe()) then
       write(*,'(a,2e18.3)') 'frac_ic_nh4',frac_ic_nh4,frac_ic_nh4_snow
       write(*,'(a,2e18.3)') 'frac_ic_no3',frac_ic_no3,frac_ic_no3_snow
       write(*,'(a,2e18.3)') 'frac_ic_so4',frac_ic_so4,frac_ic_so4_snow
    endif


 ! ELSEIF ( trop_option%aerosol_thermo .eq. AERO_ISORROPIA ) then
!        nh4no3_ndx  = get_spc_ndx( 'ANO3' )
!        if ( nh4no3_ndx .gt. 0 ) then
!           flag = query_method ('wet_deposition',MODEL_ATMOS,&
!                get_tracer_index(MODEL_ATMOS,'ano3'), &
!                text_in_scheme,text_in_param)
!           flag=parse(text_in_param,'frac_incloud',frac_ic_no3)
!           flag=parse(text_in_param,'frac_incloud_snow',frac_ic_no3_snow)
!           if (flag == 0) then
!              frac_ic_no3_snow = frac_ic_no3
!           end if
!        else
!           call error_mesg ('setsox','ano3 needs to be defined', FATAL)
!        end if
!     END IF
    isoropia_clock_id = mpp_clock_id('Chemistry: Cloud: Isoropia')

    module_is_initialized = .true.

  end subroutine setsox_init
  !>

end module MO_SETSOX_MOD
