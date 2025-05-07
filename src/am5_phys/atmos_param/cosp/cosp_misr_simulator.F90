#include "cosp_defs.H"
! version number = 1.4.3
! (c) British Crown Copyright 2008, the Met Office.
! All rights reserved.
! $Revision: 23 $, $Date: 2011-03-31 09:41:37 -0400 (Thu, 31 Mar 2011) $
! $URL: http://cfmip-obs-sim.googlecode.com/svn/stable/v1.4.0/cosp_misr_simulator.F90 $
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

!
! History:
! Nov 2008 - A. Bodas-Salcedo - Initial version
!
!

MODULE MOD_COSP_MISR_SIMULATOR
  USE MOD_COSP_CONSTANTS
  USE MOD_COSP_TYPES
  IMPLICIT NONE

INTERFACE

#ifdef COSP_GFDL
      SUBROUTINE MISR_simulator(   &
     &     npoints,   &
     &     nlev,    &
     &     ncol,   &
     &     sunlit,  &
     &     zfull,   &
     &     at,    &
     &     dtau_s,   &
     &     dtau_c,   &
     &     frac_out,     &
     &     missing_value,     &
     &     fq_MISR_TAU_v_CTH,   &
     &     dist_model_layertops,  &
     &     MISR_mean_ztop,  &
     &     MISR_cldarea,   &
     &     dtau_col   &
     & )
#else
      SUBROUTINE MISR_simulator(   &
     &     npoints,   &
     &     nlev,    &
     &     ncol,    &
     &     sunlit,   &
     &     zfull,    &
     &     at,     &
     &     dtau_s,  &
     &     dtau_c,   &
     &     frac_out,    &
     &     missing_value, & *
     &     fq_MISR_TAU_v_CTH,     &
     &     dist_model_layertops,   &
     &     MISR_mean_ztop,    &
     &     MISR_cldarea   &
     & )
#endif
    

      implicit none
      integer n_MISR_CTH
      parameter(n_MISR_CTH=16)
         
!     -----
!     Input 
!     -----

      INTEGER npoints                   !  if ncol ==1, the number of model points in the horizontal grid  
                            !   else    the number of GCM grid points
                            
      INTEGER nlev                      !  number of model vertical levels
      
      INTEGER ncol                      !  number of model sub columns 
                        !  (must already be generated in via scops and passed to this
                        !   routine via the variable frac_out )
  
      INTEGER sunlit(npoints)           !  1 for day points, 0 for night time

      REAL zfull(npoints,nlev)          !  height (in meters) of full model levels (i.e. midpoints)
                                        !  zfull(npoints,1)    is    top level of model
                                        !  zfull(npoints,nlev) is bottom level of model (closest point to surface)  

      REAL at(npoints,nlev)             !  temperature in each model level (K)
 
      REAL dtau_s(npoints,nlev)         !  visible wavelength cloud optical depth ... for "stratiform" condensate
                                        !  NOTE:  this the cloud optical depth of only the
                    !     the model cell (i,j)
                    
      REAL dtau_c(npoints,nlev)         !  visible wavelength cloud optical depth ... for "convective" condensate
                                        !  NOTE:  this the cloud optical depth of only the
                    !     the model cell (i,j)
                                     
      REAL frac_out(npoints,ncol,nlev)  !  NOTE: only need if columns>1 ... subgrid scheme in use.
      
      REAL missing_value
#ifdef COSP_GFDL
      REAL,optional ::          dtau_col(npoints,ncol,nlev)
                               ! tau values obtained from model
                               ! stochastic columns


#endif
                                 
!     ------
!     Outputs
!     ------
            
      REAL fq_MISR_TAU_v_CTH(npoints,7,n_MISR_CTH)      
      REAL dist_model_layertops(npoints,n_MISR_CTH)
      REAL MISR_cldarea(npoints)               ! fractional area coverged by clouds 
      REAL MISR_mean_ztop(npoints)             ! mean cloud top hieght(m) MISR would observe
                                   ! NOTE: == 0 if area ==0
                            

      END SUBROUTINE MISR_simulator

END INTERFACE

CONTAINS


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!-------------- SUBROUTINE COSP_MISR_SIMULATOR -----------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SUBROUTINE COSP_MISR_SIMULATOR(gbx,sgx,y)
  
  ! Arguments
  type(cosp_gridbox),intent(in) :: gbx  ! Gridbox info
  type(cosp_subgrid),intent(in) :: sgx  ! Subgridbox info
  type(cosp_misr),intent(inout) :: y    ! MISR simulator output
  
  ! Local variables 
  integer :: Nlevels,Npoints
  real :: dtau_s(gbx%Npoints, gbx%Nlevels)
  real :: dtau_c(gbx%Npoints, gbx%Nlevels)
  real :: at(gbx%Npoints, gbx%Nlevels)
  real :: frac_out(gbx%Npoints, gbx%Ncolumns, gbx%Nlevels)
  integer :: sunlit(gbx%Npoints)
  
  real :: zfull(gbx%Npoints, gbx%Nlevels) !  height (in meters) of full model levels (i.e. midpoints)
                                          !  zfull(npoints,1)    is    top level of model
                                          !  zfull(npoints,nlev) is bottom level of model
     
    
  Nlevels = gbx%Nlevels
  Npoints = gbx%Npoints
  ! Levels from TOA to surface
  zfull  = gbx%zlev(:,Nlevels:1:-1)
  at     = gbx%T(:,Nlevels:1:-1) 
  dtau_s = gbx%dtau_s(:,Nlevels:1:-1) 
  dtau_c = gbx%dtau_c(:,Nlevels:1:-1) 
  frac_out(1:Npoints,:,1:Nlevels) = sgx%frac_out(1:Npoints,:,Nlevels:1:-1)
  sunlit = int(gbx%sunlit)
 
#ifdef COSP_GFDL
 if (sgx%cols_input_from_model) then
  call MISR_simulator(gbx%npoints,gbx%nlevels,gbx%ncolumns,&
                     sunlit,zfull,at,dtau_s,dtau_c,frac_out, &
                     R_UNDEF, &
                     y%fq_MISR,y%MISR_dist_model_layertops,  &
                     y%MISR_meanztop,y%MISR_cldarea,   &
                     sgx%dtau_col)
 else
#endif
  call MISR_simulator(gbx%npoints,gbx%nlevels,gbx%ncolumns,&
                     sunlit,zfull,at,dtau_s,dtau_c,frac_out, R_UNDEF, &
                     y%fq_MISR,y%MISR_dist_model_layertops,y%MISR_meanztop,y%MISR_cldarea)
#ifdef COSP_GFDL
 endif
#endif
            
END SUBROUTINE COSP_MISR_SIMULATOR

END MODULE MOD_COSP_MISR_SIMULATOR
