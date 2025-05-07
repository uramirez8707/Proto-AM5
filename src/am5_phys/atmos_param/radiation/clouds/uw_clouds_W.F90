!FDOC_TAG_GFDL

                 module uw_clouds_W_mod
! CONTACT:  fei.liu@noaa.gov - fil
!
!          uw shallow convection cloud radiative properties module
!

use   time_manager_mod, only: time_type
use      constants_mod, only: DENS_H2O, RDGAS, TFREEZE, pi
use physics_radiation_exch_mod, only:    &
                              exchange_control_type

use cloudrad_types_mod, only: cld_specification_type, &
                              microphysics_type,  &
                              microrad_properties_type, &
                              cloudrad_control_type

!--------------------------------------------------------------------

implicit none
private

!--------------------------------------------------------------------
!          uw shallow convection cloud radiative properties module
!
!--------------------------------------------------------------------



!---------------------------------------------------------------------
!----------- ****** VERSION NUMBER ******* ---------------------------

   character(len=128)  :: version =  '$Id$'
   character(len=128)  :: tagname =  '$Name$'



!---------------------------------------------------------------------
!-------  interfaces --------

public          &
          uw_clouds_W_init,   &
          uw_clouds_W_end , uw_clouds_amt


!----------------------------------------------------------------------
!----  public data -------


!----------------------------------------------------------------------
!----  private data -------

real, parameter :: K_LAND  = 1.143
real, parameter :: K_OCEAN = 1.077
real, parameter :: QAMIN   = 1.0e-2

real, parameter :: N_LAND_ORIG  = 250.E+06
real, parameter :: N_OCEAN_ORIG = 100.E+06
real, parameter :: N_MIN_ORIG   = 1.0e06
real, parameter :: QMIN_ORIG    = 1.0e-10


real :: N_LAND, N_ocean, N_min, qmin

logical :: module_is_initialized = .false.



!----------------------------------------------------------------------
!----------------------------------------------------------------------




contains 


subroutine uw_clouds_W_init  (Exch_ctrl)

type(exchange_control_type), intent(in) :: Exch_ctrl

      integer            :: ierr, io, logunit

      if (module_is_initialized) return

      N_land = Exch_ctrl%N_land
      N_ocean = Exch_ctrl%N_ocean
      N_min = Exch_ctrl%N_min
      qmin = Exch_ctrl%qmin

!---------------------------------------------------------------------

      module_is_initialized = .true.


end subroutine uw_clouds_W_init

subroutine uw_clouds_W_end
       
!----------------------------------------------------------------------
!    uw_clouds_W_end is the destructor for uw_clouds_W_mod.
!----------------------------------------------------------------------
       
!---------------------------------------------------------------------
!    mark the module as not initialized.
!---------------------------------------------------------------------
      module_is_initialized = .false.
       
!--------------------------------------------------------------------


end subroutine uw_clouds_W_end


!#################################################################


!---------------------------------------------------------------------
subroutine uw_clouds_amt (is, ie, js, je,   &
                   shallow_cloud_area, shallow_liquid, shallow_ice, &
                   shallow_droplet_number, shallow_ice_number, land,  &
                   pfull, tkel, Cldrad_control, Shallow_microphys)

!---------------------------------------------------------------------
!    uw_clouds_amt defines the distribution of cloud water and cloud ice 
!    amounts [ g / m**3 ] and liquid and ice particle sizes and total cloud 
!    fraction for the clouds associated with uw shallow convection. these 
!    values will later be combined with other cloud fields to produce the 
!    cloud radiative properties that will be seen by the radiation package.
!----------------------------------------------------------------------

integer,                     intent(in)    :: is,ie,js,je
real, dimension(:,:,:),      intent(in)    :: shallow_cloud_area,  &
                                              shallow_liquid, shallow_ice, &
                                              shallow_droplet_number, &
                                              shallow_ice_number
real, dimension(:,:),        intent(in)    :: land
real, dimension(:,:,:),      intent(in)    :: pfull, tkel
type(cloudrad_control_type), intent(in)    :: Cldrad_control
type(microphysics_type),     intent(inout) :: Shallow_microphys

!--------------------------------------------------------------------------
!   local variables:

      real, dimension (size(shallow_liquid,1),                             &
                       size(shallow_liquid,2),                             &
                       size(shallow_liquid,3)) ::  liq_local, area_local,  &
                                                   ice_local, rho
      real, dimension (size(shallow_liquid,1),                             &
                       size(shallow_liquid,2)) ::  k_ratio, ndrops
      integer     :: ix, jx, kx
      integer     :: i, j, k

!cms++
!should better also pre-calc some const. in initialization ...
!    REAL :: dumc, dumnc, rho1, pgam, lamc, lammax, lammin, &
!           dumi, dumni, lami
 
!cms--
      
!----------------------------------------------------------------------
!    define vertical index.
!----------------------------------------------------------------------
      ix = size(Shallow_microphys%size_drop,1)
      jx = size(Shallow_microphys%size_drop,2)
      kx = size(Shallow_microphys%size_drop,3)

!---------------------------------------------------------------------
!    define k_ratio as appropriate mix of land and ocean values.
!----------------------------------------------------------------------
      k_ratio(:,:) = K_LAND*land(:,:) + K_OCEAN*(1.0 - land(:,:))
      ndrops(:,:)  = N_LAND*land(:,:) + N_OCEAN*(1.0 - land(:,:))
      
!---------------------------------------------------------------------
!    limit shallow cloud area and cloud amount so that non-zero clouds 
!    exist only in boxes with non-zero cloud areas.
!---------------------------------------------------------------------
      do k=1,kx                                      
        do j=1,jx
          do i=1,ix
            area_local(i,j,k) = 0.
            if ( (shallow_cloud_area(i,j,k) > QAMIN) .and.  &
                 (shallow_liquid(i,j,k) > QMIN) ) then
              liq_local(i,j,k) = shallow_liquid(i,j,k)
              area_local(i,j,k) = shallow_cloud_area(i,j,k)
            else
              liq_local(i,j,k) = 0.0                     
            endif

            if ( (shallow_cloud_area(i,j,k) > QAMIN) .and.  &
                 (shallow_ice(i,j,k) > QMIN)) then
              ice_local(i,j,k) = shallow_ice(i,j,k)
              area_local(i,j,k) = shallow_cloud_area(i,j,k)
            else
              ice_local(i,j,k) = 0.0                     
            endif
          end do
        end do
      end do
          
!!$!cms++
!!$   ice_num_if: if ( .NOT. Cldrad_control%do_ice_num) then
!!$!cms--


!----------------------------------------------------------------------
!    define droplet diameter based on cloud amount and droplet number. 
!    use formula as in cloud_rad_mod.
!----------------------------------------------------------------------
      if (Cldrad_control%do_liq_num) then
        do k=1,kx
          do j=1,jx
            do i=1,ix
              if (liq_local(i,j,k) > QMIN) then
                Shallow_microphys%size_drop(i,j,k)  =  2.0* &
                       (k_ratio(i,j)*620350.49*(liq_local(i,j,k)/DENS_H2O/ &
                       MAX(shallow_droplet_number(i,j,k), N_MIN*   &
                               max(area_local(i,j,k), QAMIN)/ &
                       (pfull(i,j,k)/RDGAS/tkel(i,j,k))))**(1./3.) )
              else
                Shallow_microphys%size_drop(i,j,k)  =  0.0 
              endif
            end do
          end do
        end do
        Shallow_microphys%droplet_number = shallow_droplet_number
      else

!----------------------------------------------------------------------
!  case of non-prognostic droplet number
!----------------------------------------------------------------------
        do k=1,kx
          do j=1,jx
            do i=1,ix
              if (liq_local(i,j,k) > QMIN) then
                Shallow_microphys%size_drop(i,j,k)  =    &
                 2.0* k_ratio(i,j)*620350.49*(pfull(i,j,k)*  &
                 shallow_liquid(i,j,k)/area_local(i,j,k)/RDGAS/  &
                 tkel(i,j,k)/DENS_H2O/ndrops(i,j))**(1./3.)
              else
                Shallow_microphys%size_drop(i,j,k)  =  0.0 
              endif
              Shallow_microphys%droplet_number(i,j,k) =   &
                                        ndrops(i,j)/(pfull(i,j,k)/  &
                                                    (RDGAS*tkel(i,j,k)))
            end do
          end do
        end do
      endif

!--------------------------------------------------------------------
!    if ice crystals are present, define their effective size, which
!    is a function of temperature. for ice clouds the effective radius
!    is taken from the formulation in Donner (1997, J. Geophys. Res., 
!    102, pp. 21745-21768) which is based on Heymsfield and Platt (1984)
!    with enhancement for particles smaller than 20 microns.  
!
!              T Range (K)               Reff (microns) 
!     -------------------------------    --------------
!
!     tfreeze-25. < T                      100.6           
!     tfreeze-30. < T <= Tfreeze-25.        80.8         
!     tfreeze-35. < T <= Tfreeze-30.        93.5             
!     tfreeze-40. < T <= Tfreeze-35.        63.9            
!     tfreeze-45. < T <= Tfreeze-40.        42.5           
!     tfreeze-50. < T <= Tfreeze-45.        39.9         
!     Tfreeze-55  < T <= Tfreeze-50         21.6          
!                   T <= Tfreeze-55.        20.2        
!
!--------------------------------------------------------------------
      if (Cldrad_control%using_fu2007) then
        do k=1,kx
          do j=1,jx
            do i=1,ix
!+yim Fu's parameterization of dge
              Shallow_microphys%size_ice(i,j,k) =   &
                            47.05 + 0.6624*(tkel(i,j,k) - TFREEZE) +  &
                            0.001741*(tkel(i,j,k)-TFREEZE)**2 
            end do
          end do
        end do
      else ! (using_fu2007)
        do k=1,kx
          do j=1,jx
            do i=1,ix
              if (ice_local(i,j,k) > QMIN) then
                if (tkel(i,j,k) > TFREEZE - 25. ) then
                  Shallow_microphys%size_ice(i,j,k) = 100.6    
                else if (tkel(i,j,k) >  TFREEZE - 30. .and. &
                       tkel(i,j,k) <= TFREEZE - 25.) then
                  Shallow_microphys%size_ice(i,j,k) = 80.8       
                else if (tkel(i,j,k) >  TFREEZE - 35. .and. &
                         tkel(i,j,k) <= TFREEZE - 30.) then
                  Shallow_microphys%size_ice(i,j,k)  = 93.5     
                else if (tkel(i,j,k) >  TFREEZE - 40. .and. &
                         tkel(i,j,k) <= TFREEZE - 35.) then
                  Shallow_microphys%size_ice(i,j,k) = 63.9      
                else if (tkel(i,j,k) >  TFREEZE - 45. .and. &
                         tkel(i,j,k) <= TFREEZE - 40.) then
                  Shallow_microphys%size_ice(i,j,k) = 42.5    
                else if (tkel(i,j,k) >  TFREEZE - 50. .and. &
                         tkel(i,j,k) <= TFREEZE - 45.) then
                  Shallow_microphys%size_ice(i,j,k) = 39.9       
                else if (tkel(i,j,k) >  TFREEZE - 55. .and. &
                         tkel(i,j,k) <= TFREEZE - 50.) then
                  Shallow_microphys%size_ice(i,j,k) = 21.6        
                else
                  Shallow_microphys%size_ice(i,j,k) = 20.2        
                end if
              else
                Shallow_microphys%size_ice(i,j,k) = 0.0
              endif
            end do
          end do
        end do
      endif  ! (using_fu2007)


!------------------------------------------------------------------------
!!$else  ! ice_num_if
!!$
!!$!also cal. liq. as in Morrison scheme
!!$ 
!!$      do k=1,kx
!!$        do j=1,jx
!!$          do i=1,ix
!!$
!!$
!!$! calc. based on Morrison Gettelman code
!!$!     
!!$          qamin_if: IF ( shallow_cloud_area(i,j,k) .GT. qamin ) THEN 
!!$
!!$! in-cloud mixing ratio and number conc
!!$         dumi = ice_local(i,j,k)/area_local(i,j,k)
!!$         dumni =  shallow_ice_number(i,j,k)/area_local(i,j,k)
!!$
!!$
!!$! limit in-cloud mixing ratio to reasonable value of 5 g kg-1
!!$         dumi =min(dumi ,5.e-3)
!!$
!!$!...................
!!$! cloud ice effective radius
!!$      if ( dumi > qmin ) then
!!$
!!$! add upper limit to in-cloud number concentration to prevent numerical error
!!$     dumni=min(dumni,dumi*1.e20)
!!$     lami = (gamma_mg(1.+di_mg)*ci_mg* &
!!$              dumni/dumi)**(1./di_mg)
!!$     lammax = 1./10.e-6
!!$     lammin = 1./(2.*dcs_mg)
!!$
!!$     if (lami.lt.lammin) then
!!$     lami = lammin
!!$     else if (lami.gt.lammax) then
!!$     lami = lammax
!!$     end if
!!$        Shallow_microphys%size_ice(i,j,k)  = 1.5/lami*1.e6
!!$        else
!!$        Shallow_microphys%size_ice(i,j,k)  =  0.
!!$        end if
!!$
!!$
!!$   must mult by 2 
!!$
!!$!-----------------------------------------------------------
!!$!
!!$! cloud droplet effective radius
!!$!
!!$
!!$! in-cloud mixing ratio and number conc
!!$         dumc = liq_local(i,j,k)/area_local(i,j,k)
!!$         dumnc = shallow_droplet_number(i,j,k)/area_local(i,j,k)
!!$! limit in-cloud mixing ratio to reasonable value of 5 g kg-1
!!$
!!$      dumc =min(dumc,5.e-3)
!!$
!!$
!!$      if ( dumc > qmin ) then
!!$
!!$
!!$! add upper limit to in-cloud number concentration to prevent numerical error
!!$       dumnc = min(dumnc,dumc*1.e20)
!!$
!!$         rho1=pfull(i,j,k)/(RDGAS*tkel(i,j,k))
!!$
!!$         pgam=0.0005714*(dumnc/1.e6/rho1)+0.2714
!!$         pgam=1./(pgam**2)-1.
!!$         pgam=max(pgam,2.)
!!$         pgam=min(pgam,15.)
!!$
!!$       lamc = (pi/6.*rhow*dumnc*gamma_mg(pgam+4.)/ &
!!$                 (dumc*gamma_mg(pgam+1.)))**(1./3.)
!!$       lammin = (pgam+1.)/50.e-6
!!$       lammax = (pgam+1.)/2.e-6
!!$       if (lamc.lt.lammin) then
!!$         lamc = lammin
!!$        else if (lamc.gt.lammax) then
!!$       lamc = lammax
!!$       end if
!!$       Shallow_microphys%size_drop(i,j,k)  =  gamma_mg(qcvar+1./3.)/(gamma_mg(qcvar)*qcvar**(1./3.))* &
!!$             gamma_mg(pgam+4.)/ &
!!$             gamma_mg(pgam+3.)/lamc/2.*1.e6
!!$        else
!!$        Shallow_microphys%size_drop(i,j,k)  = 0.
!!$        end if
!!$
!!$
!!$
!!$
!!$           END IF qamin_if
!!$          end do
!!$        end do
!!$      end do
!!$
!!$
!!$
!!$
!!$      Shallow_microphys%ice_number = shallow_ice_number
!!$      Shallow_microphys%droplet_number = shallow_droplet_number
!!$
!!$  
!!$
!!$    end if ice_num_if

!---------------------------------------------------------------------
!    convert the cloud and ice amounts from kg(h2o) / kg(air) to 
!    g(h2o) / m**3, as required for use in the microphys_rad routines
!    which compute cloud radiative properties.
!---------------------------------------------------------------------
      rho(:,:,:) = pfull(:,:,1:kx)/(RDGAS*tkel(:,:,1:kx))
      Shallow_microphys%conc_drop = 1.0e03*rho*liq_local  
      Shallow_microphys%conc_ice  = 1.0e03*rho* ice_local  

!---------------------------------------------------------------------
!    define the cloud area to be that after adjustment for trivial cloud
!    amounts.
!---------------------------------------------------------------------
      Shallow_microphys%cldamt = area_local


!---------------------------------------------------------------------



end subroutine uw_clouds_amt  



!####################################################################


                     end module uw_clouds_W_mod

