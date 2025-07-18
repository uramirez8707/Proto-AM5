  MODULE STABLE_BL_TURB_MOD

!=======================================================================
 use           mpp_mod, only: input_nml_file
 use           fms_Mod, ONLY: ERROR_MESG, FATAL, mpp_pe, mpp_root_pe,  &
                              check_nml_error, write_version_number,   &
                              stdlog
 use  Diag_Manager_Mod, ONLY: register_diag_field, send_data
 use  Time_Manager_Mod, ONLY: time_type
 use     Constants_Mod, ONLY: cp_air, hlv, hls, grav, vonkarm, tfreeze,&
                              rdgas, rvgas, omega
 use Monin_Obukhov_Mod, ONLY: stable_mix

 implicit none
 private
 public :: STABLE_BL_TURB, STABLE_BL_TURB_INIT, STABLE_BL_TURB_END


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  character(len=128) :: version = '$Id$'
  character(len=128) :: tagname = '$Name$'
  logical            :: module_is_initialized = .false.
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!---------------------------------------------------------------------
! --- CONSTANTS
!---------------------------------------------------------------------

  real :: oalsm, oalsh
  real, parameter :: d608 = (rvgas-rdgas)/rdgas
  
!---------------------------------------------------------------------
! --- NAMELIST
!---------------------------------------------------------------------

  real    :: akmax        = 1.e4 ! maximum diffusion coefficient value 
                                 !  (m2/s)
  real    :: alpha        = 0.5
  real    :: alsm         = 500.0 !500m may be too large --
                                  !the original Louis (1979) paper used 100m,
                                  !and ECMWF-IFS uses 150m.
  real    :: alsh         = 500.0 !500m may be too large --
                                  !the original Louis (1979) paper used 100m,
                                  !and ECMWF-IFS uses 150m.
  real    :: fmin         = 5.0e-5
  real    :: hpbl_cap     = 1000.
  real    :: ri_crit      = 0.2
  real    :: diff_min     = 0.001
  real    :: winddifmin   = 0.01
  real    :: small        = 1.e-5
  real    :: b_louis      = 9.4
  real    :: cmstar_louis = 7.4
  real    :: chstar_louis = 5.3

! add by danli
  real    :: pr           = 1
  real    :: lRi          = 0.5
  real    :: alsm_HBC     = 7
  logical :: do_stable_HBC = .false.



  NAMELIST / stable_bl_turb_nml / akmax, alpha, alsm, alsh, fmin, &
                                  hpbl_cap, diff_min, ri_crit, &
                                  pr, lRi, alsm_HBC, do_stable_HBC

!---------------------------------------------------------------------
! DIAGNOSTICS FIELDS 
!---------------------------------------------------------------------

integer :: id_z_sbl, id_f_sbl

character(len=14) :: mod_name = 'stable_bl_turb'

real :: missing_value = -999.

!---------------------------------------------------------------------
 contains

!#######################################################################
 subroutine STABLE_BL_TURB( is, js, Time, temp, qv,  ql,  qi,  um,  vm,&
                            zhalf, zfull, u_star, b_star, lat,         &
                            akm, akh, vspblcap, kbot )

!=======================================================================
!---------------------------------------------------------------------
! Arguments (Intent in)
!
!       is, js   -  Starting indices for window
!       Time     -  Time used for diagnostics [time_type]
!       temp     -  temperature (K)
!       qv       -  specific humidity of water vapor (kg/kg)
!       ql       -  specific humidity of cloud liquid (kg/kg)
!       qi       -  specific humidity of cloud ice (kg/kg)
!       um, vm   -  Wind components (m/s)
!       u_star   -  surface friction velocity (m/s)
!       b_star   -  surface buoyancy (m/s2)
!       lat      -  latitude in radians
!       zhalf    -  Height at half levels (m)
!       zfull    -  Height at full levels (m)
!
!      --------------
!      optional input
!      --------------
!
!      vspblcap  cap to height of very stable PBL mixing, coming
!                from any other module (m) (in usual application
!                this might be entrain_mod or edt mod)
!
!      kbot      integer indicating the lowest true layer of atmosphere
!                this is used only for eta coordinate model
!
!---------------------------------------------------------------------
! Arguments (Intent out)
!       akm  -  mixing coefficient for momentum
!       akh  -  mixing coefficient for heat and moisture
!---------------------------------------------------------------------

  type(time_type), intent(in)                    :: Time
  integer,         intent(in)                    :: is, js
  real,            intent(in),  dimension(:,:)   :: u_star, b_star, lat
  real,            intent(in),  dimension(:,:,:) :: temp, qv, ql, qi
  real,            intent(in),  dimension(:,:,:) :: um, vm 
  real,            intent(in),  dimension(:,:,:) :: zhalf,  zfull
  real,            intent(out), dimension(:,:,:) :: akm,    akh

  real,     intent(in),   dimension(:,:), optional :: vspblcap
  integer,  intent(in),   dimension(:,:), optional :: kbot
  
!---------------------------------------------------------------------

  real, dimension(SIZE(um,1),SIZE(um,2),SIZE(um,3)-1) ::             &
        dsdzh, shear, buoync, Ri, Ritmp, fm, fh, lm, lh, xxm1, xxm2, &
        phi, zfunc, cmtmp, chtmp, fmtmp, fhtmp

  real, dimension(SIZE(um,1),SIZE(um,2),SIZE(um,3))  :: mask, hleff, akm_hbc
  real, dimension(SIZE(um,1),SIZE(um,2),SIZE(um,3))  :: zfull_ag, slv
  real, dimension(SIZE(um,1),SIZE(um,2),SIZE(um,3)+1):: zhalf_ag

  real, dimension(SIZE(um,1),SIZE(um,2)) ::  zsurf,  &
        fcor, hpbl, z_sbl, f_sbl

  integer :: ix, jx, kx, i, j, k,  kxm
  integer :: shape1(1), shape3(3)
  logical :: used

!=======================================================================
! --- Initalize
!=======================================================================

! --- Check to see if STABLE_BL_TURB has been initialized
  if( .not. module_is_initialized ) CALL ERROR_MESG( ' STABLE_BL_TURB',                          &
       ' STABLE_BL_TURB_INIT has not been called', FATAL)

! --- Zero out output arrays
    akm(:,:,:) = 0.0
    akh(:,:,:) = 0.0
    akm_hbc(:,:,:) = 0.0
  z_sbl(:,:)   = 0.0
  f_sbl(:,:)   = 0.0
  
! --- Set dimensions etc
  ix  = SIZE( um, 1 )
  jx  = SIZE( um, 2 )
  kx  = SIZE( um, 3 )
  kxm = kx - 1

  shape1 =    ix * jx * kxm
  shape3 = (/ ix,  jx,  kxm /)

!====================================================================
! --- COMPUTE HEIGHT ABOVE SURFACE            
!====================================================================


       mask = 1.0
                   
       if (present(kbot)) then
            do j=1,jx
            do i=1,ix
                 zsurf(i,j) = zhalf(i,j,kbot(i,j)+1)
                 if (kbot(i,j).lt.kx) then
                    do k = kbot(i,j)+1,kx
                       mask(i,j,k) = 0.0
                    enddo
                 end if      
            enddo
            enddo
       else
            zsurf(:,:) = zhalf(:,:,kx+1)
       end if

       do k = 1, kx
            zfull_ag(:,:,k) = zfull(:,:,k) - zsurf(:,:)
            zhalf_ag(:,:,k) = zhalf(:,:,k) - zsurf(:,:)
       end do
       zhalf_ag(:,:,kx+1) = zhalf(:,:,kx+1) - zsurf(:,:)
       
!====================================================================
! --- DYNAMIC HEIGHT - also height relative to the surface     
!====================================================================

  fcor(:,:) = 2.0 * omega * SIN( lat(:,:) )
  fcor(:,:) = ABS( fcor(:,:) )
  fcor(:,:) = MAX( fcor(:,:), fmin )
  hpbl(:,:) = alpha * u_star(:,:) / fcor(:,:)

! --- bound
  hpbl(:,:) = MIN( hpbl(:,:), hpbl_cap )
  
! --- cap from entrainment turbulence
  if (present(vspblcap)) hpbl(:,:) = MIN( hpbl(:,:), vspblcap )
  
! --- height relative to the surface
! --- zfunc = zhalf_ag / hpbl, where stable conditions exist
! ---       = 1.0 for unstable conditions

  zfunc = 1.0
  do k = 1, kxm
      where( b_star(:,:) < 0.0) 
           zfunc(:,:,k) = min(1.,max(0.,zhalf_ag(:,:,k+1)/             &
                                     max(0.1,hpbl(:,:))))
      endwhere
  enddo
    
!====================================================================
! --- COMPUTE LIQUID WATER VIRTUAL STATIC ENERGY             
!====================================================================

   hleff   = (min(1.,max(0.,0.05*(temp   -tfreeze+20.)))*hlv + &
              min(1.,max(0.,0.05*(tfreeze -temp      )))*hls)
     
   slv     = cp_air*temp + grav*zfull_ag - hleff*(ql + qi)
   slv     = slv*(1+d608*(qv+ql+qi))
       
!====================================================================
! --- COMPUTE RICHARDSON NUMBER                 
!====================================================================

! --- D( )/DZ OPERATOR  
  
  dsdzh(:,:,1:kxm) = 1.0 / (zfull_ag(:,:,1:kxm) - zfull_ag(:,:,2:kx))

! --- WIND SHEAR SQUARED

  xxm1(:,:,1:kxm) = dsdzh(:,:,1:kxm)*( um(:,:,1:kxm) - um(:,:,2:kx) )
  xxm2(:,:,1:kxm) = dsdzh(:,:,1:kxm)*( vm(:,:,1:kxm) - vm(:,:,2:kx) )



  shear(:,:,:) = xxm1(:,:,:)*xxm1(:,:,:) + xxm2(:,:,:)*xxm2(:,:,:)

  where (shear .lt. (dsdzh*winddifmin*dsdzh*winddifmin)) 
         shear = dsdzh*winddifmin*dsdzh*winddifmin
  end where         

! --- BUOYANCY 
  xxm1(:,:,1:kxm) =       slv(:,:,1:kxm) - slv(:,:,2:kx) 
  xxm2(:,:,1:kxm) = 0.5*( slv(:,:,1:kxm) + slv(:,:,2:kx) )
 
 
  buoync(:,:,:) = grav * dsdzh(:,:,:) * xxm1(:,:,:) / xxm2(:,:,:)

! --- RICHARDSON NUMBER

  Ri(:,:,:) = buoync(:,:,:) / shear(:,:,:)   

!====================================================================
! --- MASK OUT UNDERGROUND VALUES FOR ETA COORDINATE
!====================================================================

  if( PRESENT( kbot ) ) then
     shear(:,:,1:kxm) =  shear(:,:,1:kxm) * mask(:,:,2:kx) 
    buoync(:,:,1:kxm) = buoync(:,:,1:kxm) * mask(:,:,2:kx) 
        Ri(:,:,1:kxm) =     Ri(:,:,1:kxm) * mask(:,:,2:kx) 
  endif

!====================================================================
! --- MIXING LENGTHS                 
!====================================================================

 do k = 1,kxm
   xxm1(:,:,k) = 1.0 / (vonkarm*zhalf_ag(:,:,k+1))
 end do 

  lm(:,:,:) = 1.0 / ( xxm1(:,:,1:kxm) + oalsm )
  lh(:,:,:) = 1.0 / ( xxm1(:,:,1:kxm) + oalsh )

!====================================================================
! --- STABILITY FUNCTIONS : STABLE SIDE       
!
! Note the very stable form of stability function acquired from 
! monin obukhov is weighted with the traditional stable form
! (phi = 1 + zeta/zeta_crit  or fm = (1 - Ri/Ricrit)**2)
! For Ricrit   = 0.2, phi retains the usual 1 + 5*zeta form.
!
! For heights greater than hpbl, the usual form is used.  For heights
! less than hpbl, the weight of the traditional form is given by
! zfunc which is linear in z/hpbl (see code above).
!====================================================================

  Ritmp = Ri
  where (Ritmp .lt. small) Ritmp = small

  CALL STABLE_MIX( Ritmp, phi)

  phi = (1-zfunc)*phi + zfunc* ((1-min(1.,(Ritmp/ri_crit)))**2.)
  
  fm(:,:,:) = phi(:,:,:)
  fh(:,:,:) =  fm(:,:,:) 
  
!====================================================================
! --- STABILITY FUNCTIONS : UNSTABLE SIDE (Louis 1979)
!
! f = 1.  - b * Ri / (1 + c*sqrt(-Ri))
!
! where b = 9.4 and
!
!              l * l * b * ( (  (1+(dz/z))**(1/3) - 1 )**(3/2))
! c = C_star * --------------------------------------------------
!              sqrt(z) * (dz**3/2)  
!
! where C_star(momentum) = 7.4, and C_star(heat) = 5.3
!     
!====================================================================
 
  Ritmp = Ri
  where (Ri .gt. 0.) Ritmp = 0.  
  
  zfunc(:,:,1:kxm) = 1.+(1./(dsdzh(:,:,1:kxm)*zhalf_ag(:,:,2:(kxm+1))))
  zfunc = zfunc **(1./3.) - 1.
  zfunc = zfunc **1.5
  zfunc = zfunc /  sqrt(zhalf_ag(:,:,2:(kxm+1))) 
  zfunc = zfunc * ( dsdzh(:,:,1:kxm) ** 1.5 )
  
  cmtmp =  cmstar_louis*lm(:,:,:)*lm(:,:,:)*b_louis*zfunc(:,:,:)
  chtmp =  chstar_louis*lh(:,:,:)*lh(:,:,:)*b_louis*zfunc(:,:,:)
  fmtmp(:,:,:) = 1. - (b_louis*Ritmp/(1.+cmtmp*sqrt(-1.*Ritmp)))
  fhtmp(:,:,:) = 1. - (b_louis*Ritmp/(1.+chtmp*sqrt(-1.*Ritmp)))
  
  where (Ri .lt. small)
      fm = fmtmp
      fh = fhtmp
  end where      
 
!====================================================================
! --- MIXING COEFFICENTS                 
!====================================================================

  shear(:,:,:) = SQRT( shear(:,:,:) )


! --- Momentum
   xxm1(:,:,:)    = lm(:,:,:) * lm(:,:,:) * fm(:,:,:)
   akm(:,:,2:kx) = xxm1(:,:,1:kxm) * shear(:,:,1:kxm) 
! --- Heat and Moisture
  xxm1(:,:,:)    = lh(:,:,:) * lh(:,:,:) * fh(:,:,:)
   akh(:,:,2:kx) = xxm1(:,:,1:kxm) * shear(:,:,1:kxm)

! the following code is modified by danli
! note the modifications only occur in the atmospheric boundary layer
! above the atmopsheric boundary layer, no modifications 


  if (do_stable_HBC) then


	where ( Ri .gt. 0.)

        akm(:,:,2:kx) =  (1.0 / (  1.0 / (vonkarm*zhalf_ag(:,:,2:kx)) + 1.0/alsm_HBC + Ri(:,:,1:kxm) / lRi) ) **2 * shear(:,:,1:kxm)
	akh(:,:,2:kx) =  akm(:,:,2:kx)/pr 
	
        end where

  endif   


       




  where (akm .lt. diff_min) akm = 0.0
  where (akm .gt. akmax) akm = akmax
  where (akh .lt. diff_min) akh = 0.0
  where (akh .gt. akmax) akh = akmax

!====================================================================
! --- Extra diagnostics
!====================================================================

  where( b_star(:,:)  < 0.0 .and. hpbl (:,:) > 0.0 )
          z_sbl(:,:) = hpbl(:,:)
          f_sbl(:,:) = 1.0
  endwhere
  
  if ( id_z_sbl > 0 ) then
     used = send_data ( id_z_sbl, z_sbl, Time, is, js )
  endif
  if ( id_f_sbl > 0 ) then
     used = send_data ( id_f_sbl, f_sbl, Time, is, js )
  endif
  
!=======================================================================
  end subroutine STABLE_BL_TURB

subroutine STABLE_BL_TURB_INIT ( axes, Time )
!=======================================================================
!     Initializes stable_bl_turb_mod: Reads and records namelist, 
!     sets up netcdf output if desired.

   
 integer,         intent(in) :: axes(4) !< Vector of axes integers
 type(time_type), intent(in) :: Time    !< time variable

 integer :: logunit, io, ierr

!=======================================================================

!---------------------------------------------------------------------
! --- Read namelist
!---------------------------------------------------------------------

  read (input_nml_file, nml=stable_bl_turb_nml, iostat=io)
  ierr = check_nml_error(io,'stable_bl_turb_nml')

!---------------------------------------------------------------------
! --- Output version
!---------------------------------------------------------------------

  if ( mpp_pe() == mpp_root_pe() ) then
       call write_version_number(version, tagname)
       logunit = stdlog()
       WRITE( logunit, nml = stable_bl_turb_nml ) 
  endif

!---------------------------------------------------------------------
! --- CONSTANTS
!---------------------------------------------------------------------

  oalsm = 1.0 / alsm
  oalsh = 1.0 / alsh

!---------------------------------------------------------------------
! --- initialize quantities for diagnostics output
!---------------------------------------------------------------------

   id_z_sbl = register_diag_field ( mod_name, &
     'z_sbl', axes(1:2), Time, &
     'Depth of stable boundary layer',              'm' )

   id_f_sbl = register_diag_field ( mod_name, &
     'f_sbl', axes(1:2), Time, &
     'Frequency of stable boundary layer',          ' ' )

!---------------------------------------------------------------------
 module_is_initialized = .true.
!=======================================================================
 end subroutine STABLE_BL_TURB_INIT

!#######################################################################
subroutine STABLE_BL_TURB_END
!    Closes down stable_bl_turb.
!---------------------------------------------------------------------
 module_is_initialized = .false.
!=======================================================================
 end subroutine STABLE_BL_TURB_END

!#######################################################################
  end MODULE STABLE_BL_TURB_MOD        

