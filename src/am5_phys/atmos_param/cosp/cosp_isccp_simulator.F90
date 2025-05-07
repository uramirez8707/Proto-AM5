#include "cosp_defs.H"
! version number = 1.4.3
! (c) British Crown Copyright 2008, the Met Office.
! All rights reserved.
! $Revision: 23 $, $Date: 2011-03-31 09:41:37 -0400 (Thu, 31 Mar 2011) $
! $URL: http://cfmip-obs-sim.googlecode.com/svn/stable/v1.4.0/cosp_isccp_simulator.F90 $
! 
! Redistribution and use in source and binary forms, with or without modification, are permitted 
! provided that the following conditions are met:
! 
!     * Redistributions of source code must retain the above copyright notice, this list 
!       of conditions and the following disclaimer.
!     * Redistributions in binary form must reproduce the above copyright notice, this list
!       of conditions and the following disclaimer in the documentation and/or other materials 
!       provided with the distribution.
!     * Neither the name of the Met Office nor the names of its contributors may be used 
!       to endorse or promote products derived from this software without specific prior written 
!       permission.
! 
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR 
! IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND 
! FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR 
! CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL 
! DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, 
! DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER 
! IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT 
! OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

MODULE MOD_COSP_ISCCP_SIMULATOR
  USE MOD_COSP_CONSTANTS
  USE MOD_COSP_TYPES
  IMPLICIT NONE

INTERFACE 
#ifdef COSP_GFDL
      SUBROUTINE ICARUS(          &
     &     debug,                 &
     &     debugcol,              &
     &     npoints,               &
     &     sunlit,                &
     &     nlev,                  &
     &     ncol,                  &
     &     pfull,                 &
     &     phalf,                 &
     &     qv,                    &
     &     cc,                    &
     &     conv,                  &
     &     dtau_s,                &
     &     dtau_c,                &
     &     top_height,            &
     &     top_height_direction,  &
     &     overlap,               &
     &     frac_out,              &
     &     skt,                   &
     &     emsfc_lw,              &
     &     at,                    &
     &     dem_s,                 &
     &     dem_c,                 &
     &     fq_isccp,              &
     &     totalcldarea,          &
     &     meanptop,              &
     &     meantaucld,            &
     &     meanalbedocld,         &
     &     meantb,                &
     &     meantbclr,             &
     &     boxtau,                &
     &     boxptop,               &
     &     dtau_col,              &
     &     dem_col                &
     &)
#else
      SUBROUTINE ICARUS(
     &     debug,
     &     debugcol,
     &     npoints,
     &     sunlit,
     &     nlev,
     &     ncol,
     &     pfull,
     &     phalf,
     &     qv,
     &     cc,
     &     conv,
     &     dtau_s,
     &     dtau_c,
     &     top_height,
     &     top_height_direction,
     &     overlap,
     &     frac_out,
     &     skt,
     &     emsfc_lw,
     &     at,
     &     dem_s,
     &     dem_c,
     &     fq_isccp,
     &     totalcldarea,
     &     meanptop,
     &     meantaucld,
     &     meanalbedocld,
     &     meantb,
     &     meantbclr,
     &     boxtau,
     &     boxptop
     &)
#endif


#ifdef COSP_GFDL
use mpp_mod,only: get_unit
use fms_mod,only: stdlog, error_mesg, FATAL
#endif
      implicit none

!     NOTE:   the maximum number of levels and columns is set by
!             the following parameter statement

      INTEGER ncolprint
      
!     -----
!     Input 
!     -----

      integer debug       ! set to non-zero value to print out inputs
                    ! with step debug
      integer debugcol    ! set to non-zero value to print out column
                 ! decomposition with step debugcol
      INTEGER npoints       !  number of model points in the horizontal
      INTEGER nlev          !  number of model levels in column
      INTEGER ncol          !  number of subcolumns

      INTEGER sunlit(npoints) !  1 for day points, 0 for night time

      REAL pfull(npoints,nlev)
                       !  pressure of full model levels (Pascals)
                  !  pfull(npoints,1) is top level of model
                  !  pfull(npoints,nlev) is bot of model

      REAL phalf(npoints,nlev+1)
                  !  pressure of half model levels (Pascals)
                  !  phalf(npoints,1) is top of model
                  !  phalf(npoints,nlev+1) is the surface pressure

      REAL qv(npoints,nlev)
                  !  water vapor specific humidity (kg vapor/ kg air)
                  !         on full model levels

      REAL cc(npoints,nlev)   
                  !  input cloud cover in each model level (fraction) 
                  !  NOTE:  This is the HORIZONTAL area of each
                  !         grid box covered by clouds

      REAL conv(npoints,nlev) 
                  !  input convective cloud cover in each model
                  !   level (fraction) 
                  !  NOTE:  This is the HORIZONTAL area of each
                  !         grid box covered by convective clouds

      REAL dtau_s(npoints,nlev) 
                  !  mean 0.67 micron optical depth of stratiform
                !  clouds in each model level
                  !  NOTE:  this the cloud optical depth of only the
                  !  cloudy part of the grid box, it is not weighted
                  !  with the 0 cloud optical depth of the clear
                  !         part of the grid box

      REAL dtau_c(npoints,nlev) 
                  !  mean 0.67 micron optical depth of convective
                !  clouds in each
                  !  model level.  Same note applies as in dtau_s.

      INTEGER overlap                   !  overlap type
                              !  1=max
                              !  2=rand
                              !  3=max/rand

      INTEGER top_height                !  1 = adjust top height using both a computed
                                        !  infrared brightness temperature and the visible
                              !  optical depth to adjust cloud top pressure. Note
                              !  that this calculation is most appropriate to compare
                              !  to ISCCP data during sunlit hours.
                                        !  2 = do not adjust top height, that is cloud top
                                        !  pressure is the actual cloud top pressure
                                        !  in the model
                              !  3 = adjust top height using only the computed
                              !  infrared brightness temperature. Note that this
                              !  calculation is most appropriate to compare to ISCCP
                              !  IR only algortihm (i.e. you can compare to nighttime
                              !  ISCCP data with this option)

      INTEGER top_height_direction ! direction for finding atmosphere pressure level
                                 ! with interpolated temperature equal to the radiance
				 ! determined cloud-top temperature
				 !
				 ! 1 = find the *lowest* altitude (highest pressure) level
				 ! with interpolated temperature equal to the radiance
				 ! determined cloud-top temperature
				 !
				 ! 2 = find the *highest* altitude (lowest pressure) level
				 ! with interpolated temperature equal to the radiance 
				 ! determined cloud-top temperature
				 ! 
				 ! ONLY APPLICABLE IF top_height EQUALS 1 or 3
				 !				 !
				 ! 1 = old setting: matches all versions of 
				 ! ISCCP simulator with versions numbers 3.5.1 and lower
				 !
				 ! 2 = default setting: for version numbers 4.0 and higher
!
!     The following input variables are used only if top_height = 1 or top_height = 3
!
      REAL skt(npoints)                 !  skin Temperature (K)
      REAL emsfc_lw                     !  10.5 micron emissivity of surface (fraction)                                            
      REAL at(npoints,nlev)                   !  temperature in each model level (K)
      REAL dem_s(npoints,nlev)                !  10.5 micron longwave emissivity of stratiform
                              !  clouds in each
                                        !  model level.  Same note applies as in dtau_s.
      REAL dem_c(npoints,nlev)                  !  10.5 micron longwave emissivity of convective
                              !  clouds in each
                                        !  model level.  Same note applies as in dtau_s.

      REAL frac_out(npoints,ncol,nlev) ! boxes gridbox divided up into
                              ! Equivalent of BOX in original version, but
                              ! indexed by column then row, rather than
                              ! by row then column

#ifdef COSP_GFDL
       REAL, optional :: dtau_col(npoints,ncol,nlev)
                              ! tau values obtained from model
                              ! stochastic columns

       REAL, optional :: dem_col(npoints,ncol,nlev)
                              ! lw emissivity values obtained
                              ! from model stochastic columns
 

#endif


!     ------
!     Output
!     ------

      REAL fq_isccp(npoints,7,7)        !  the fraction of the model grid box covered by
                                        !  each of the 49 ISCCP D level cloud types

      REAL totalcldarea(npoints)        !  the fraction of model grid box columns
                                        !  with cloud somewhere in them.  NOTE: This diagnostic
					! does not count model clouds with tau < isccp_taumin
                              ! Thus this diagnostic does not equal the sum over all entries of fq_isccp.
			      ! However, this diagnostic does equal the sum over entries of fq_isccp with
			      ! itau = 2:7 (omitting itau = 1)
      
      
      ! The following three means are averages only over the cloudy areas with tau > isccp_taumin.  
      ! If no clouds with tau > isccp_taumin are in grid box all three quantities should equal zero.      
                              
      REAL meanptop(npoints)            !  mean cloud top pressure (mb) - linear averaging
                                        !  in cloud top pressure.
                              
      REAL meantaucld(npoints)          !  mean optical thickness 
                                        !  linear averaging in albedo performed.
      
      real meanalbedocld(npoints)        ! mean cloud albedo
                                        ! linear averaging in albedo performed
					
      real meantb(npoints)              ! mean all-sky 10.5 micron brightness temperature
      
      real meantbclr(npoints)           ! mean clear-sky 10.5 micron brightness temperature
      
      REAL boxtau(npoints,ncol)         !  optical thickness in each column
      
      REAL boxptop(npoints,ncol)        !  cloud top pressure (mb) in each column

    end subroutine icarus

END INTERFACE 
CONTAINS


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!-------------- SUBROUTINE COSP_ISCCP_SIMULATOR -----------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SUBROUTINE COSP_ISCCP_SIMULATOR(gbx,sgx,y)
  
  ! Arguments
  type(cosp_gridbox),intent(in) :: gbx  ! Gridbox info
  type(cosp_subgrid),intent(in) :: sgx  ! Subgridbox info
  type(cosp_isccp),intent(inout) :: y   ! ISCCP simulator output
  
  ! Local variables 
  integer :: Nlevels,Npoints
  real :: pfull(gbx%Npoints, gbx%Nlevels)
  real :: phalf(gbx%Npoints, gbx%Nlevels + 1)
  real :: qv(gbx%Npoints, gbx%Nlevels)
  real :: cc(gbx%Npoints, gbx%Nlevels)
  real :: conv(gbx%Npoints, gbx%Nlevels)
  real :: dtau_s(gbx%Npoints, gbx%Nlevels)
  real :: dtau_c(gbx%Npoints, gbx%Nlevels)
  real :: at(gbx%Npoints, gbx%Nlevels)
  real :: dem_s(gbx%Npoints, gbx%Nlevels)
  real :: dem_c(gbx%Npoints, gbx%Nlevels)
  real :: frac_out(gbx%Npoints, gbx%Ncolumns, gbx%Nlevels)
  integer :: sunlit(gbx%Npoints)
  
  Nlevels = gbx%Nlevels
  Npoints = gbx%Npoints
  ! Flip inputs. Levels from TOA to surface
  pfull  = gbx%p(:,Nlevels:1:-1) 
  phalf(:,1)         = 0.0 ! Top level
  phalf(:,2:Nlevels+1) = gbx%ph(:,Nlevels:1:-1)
  qv     = gbx%sh(:,Nlevels:1:-1) 
  cc     = 0.999999*gbx%tca(:,Nlevels:1:-1) 
  conv   = 0.999999*gbx%cca(:,Nlevels:1:-1) 
  dtau_s = gbx%dtau_s(:,Nlevels:1:-1) 
  dtau_c = gbx%dtau_c(:,Nlevels:1:-1) 
  at     = gbx%T(:,Nlevels:1:-1) 
  dem_s  = gbx%dem_s(:,Nlevels:1:-1) 
  dem_c  = gbx%dem_c(:,Nlevels:1:-1) 
  frac_out(1:Npoints,:,1:Nlevels) = sgx%frac_out(1:Npoints,:,Nlevels:1:-1)
  sunlit = int(gbx%sunlit)
#ifdef COSP_GFDL
  if (sgx%cols_input_from_model) then
    call icarus(0,0,gbx%npoints,sunlit,gbx%nlevels,gbx%ncolumns, &
            pfull,phalf,qv,cc,conv,dtau_s,dtau_c, &
            gbx%isccp_top_height,gbx%isccp_top_height_direction, &
            gbx%isccp_overlap,frac_out, &
            gbx%skt,gbx%isccp_emsfc_lw,at,dem_s,dem_c,y%fq_isccp,y%totalcldarea, &
            y%meanptop,y%meantaucld,y%meanalbedocld, &
            y%meantb,y%meantbclr,y%boxtau,y%boxptop, &
            sgx%dtau_col, sgx%dem_col)
  else
    call icarus(0,0,gbx%npoints,sunlit,gbx%nlevels,gbx%ncolumns, &
            pfull,phalf,qv,cc,conv,dtau_s,dtau_c, &
            gbx%isccp_top_height,gbx%isccp_top_height_direction, &
            gbx%isccp_overlap,frac_out, &
            gbx%skt,gbx%isccp_emsfc_lw,at,dem_s,dem_c,y%fq_isccp,y%totalcldarea, &
            y%meanptop,y%meantaucld,y%meanalbedocld, &
            y%meantb,y%meantbclr,y%boxtau,y%boxptop)
   endif
#else
    call icarus(0,0,gbx%npoints,sunlit,gbx%nlevels,gbx%ncolumns, &
            pfull,phalf,qv,cc,conv,dtau_s,dtau_c, &
            gbx%isccp_top_height,gbx%isccp_top_height_direction, &
            gbx%isccp_overlap,frac_out, &
            gbx%skt,gbx%isccp_emsfc_lw,at,dem_s,dem_c,y%fq_isccp,y%totalcldarea, &
            y%meanptop,y%meantaucld,y%meanalbedocld, &
            y%meantb,y%meantbclr,y%boxtau,y%boxptop)
#endif

  ! Flip outputs. Levels from surface to TOA
  ! --- (npoints,tau=7,pressure=7)
  y%fq_isccp(:,:,:) = y%fq_isccp(:,:,7:1:-1)
     
 
  ! Check if there is any value slightly greater than 1
  where ((y%totalcldarea > 1.0-1.e-5) .and. (y%totalcldarea < 1.0+1.e-5))
    y%totalcldarea = 1.0
  endwhere
              
END SUBROUTINE COSP_ISCCP_SIMULATOR

END MODULE MOD_COSP_ISCCP_SIMULATOR
