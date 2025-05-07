                  module longwave_fluxes_mod
!
! CONTACT: Fei Liu
! REVIEWER: Dan.Schwarzkopf@noaa.gov - Dan Schwartzkopf
!  This code is a helper module that provides various operations on 
!  longwave flux variables.

!  shared modules:

use fms_mod,            only: fms_init, &
                              error_mesg, &
                              FATAL

!  shared radiation package modules:

use longwave_types_mod, only: lw_diagnostics_type

!---------------------------------------------------------------------

implicit none
private

!---------------------------------------------------------------------
!    longwave_fluxes calculates the longwave fluxes between each model
!    level and every other modle level for each of the longwave
!    spectral parameterization bands.
!---------------------------------------------------------------------


!---------------------------------------------------------------------
!----------- version number for this module -------------------

character(len=128)  :: version =  '$Id$'
character(len=128)  :: tagname =  '$Name$'


!---------------------------------------------------------------------
!-------  interfaces --------

public    &
       longwave_fluxes_init, &
       longwave_fluxes_ks, longwave_fluxes_k_down, &
       longwave_fluxes_KE_KEp1, longwave_fluxes_diag, &
       longwave_fluxes_sum, longwave_fluxes_end



!---------------------------------------------------------------------
!------- public data ------


!---------------------------------------------------------------------
!------- private data ------

logical :: module_is_initialized = .false. ! module is initialized ?


!---------------------------------------------------------------------
!---------------------------------------------------------------------



                         contains

subroutine longwave_fluxes_init 

!---------------------------------------------------------------------
!    longwave_fluxes_init is the constructor for longwave_fluxes_mod.
!---------------------------------------------------------------------

!---------------------------------------------------------------------
     integer    ::  ierr, io, logunit

!---------------------------------------------------------------------
!  local variables:
!
!        ierr            error code
!        io              error status returned from io operation
!
!---------------------------------------------------------------------

!---------------------------------------------------------------------
!    if routine has already been executed, exit.
!---------------------------------------------------------------------
      if (module_is_initialized) return
 
!---------------------------------------------------------------------
!    verify that modules used by this module that are not called later
!    have already been initialized.
!---------------------------------------------------------------------
      call fms_init

!---------------------------------------------------------------------
!    mark the module as initialized.
!---------------------------------------------------------------------
      module_is_initialized = .true.

!---------------------------------------------------------------------

end subroutine longwave_fluxes_init



!#####################################################################
subroutine longwave_fluxes_ks (source, trans, source2, trans2,  &
                               cld_trans, cld_ind, do_totcld_forcing, &
                               Lw_diagnostics)

!---------------------------------------------------------------------
! subroutine to calculate longwave diagnostic fluxes
!---------------------------------------------------------------------

!---------------------------------------------------------------------
integer, dimension(:),      intent(in)    :: cld_ind
real, dimension (:,:,:,:),  intent(in)    :: source
real, dimension (:,:,:,:),  intent(in)    :: source2
real, dimension (:,:,:,:),  intent(in)    :: trans2, trans
real, dimension (:,:,:,:),  intent(in)    :: cld_trans
logical,                    intent(in)    :: do_totcld_forcing
type(lw_diagnostics_type),  intent(inout) :: Lw_diagnostics
!---------------------------------------------------------------------

!-------------------------------------------------------------------
!  intent(in) variables:
!
!     cld_ind  lookup table to tranlate longwave band index to cloud index
!     source   longwave source function
!     source2  longwave source functioin
!     trans    longwave transmittance function
!     trans2   longwave transmittance function
!     cld_trans longwave cloud transmittance function
!
!  intent(inout) variables:
!
!     Lw_diagnostics  contains the longwave diagnostic flux values
!
!---------------------------------------------------------------------

!------------------------------------------------------------------
!  local variables:

      real, dimension (size(source2,1), &
                       size(source2,2), &
                       size(source2,3) ) ::    flux_tmp, flux_tmp2

      integer   ::   k, ks, ke, nbands, m

!---------------------------------------------------------------------
!  local variables:
!
!      flux_tmp
!      flux_tmp2
!      k
!      ks
!      ke
!      nbands
!      m
!
!--------------------------------------------------------------------

!---------------------------------------------------------------------
!    be sure module has been initialized.
!---------------------------------------------------------------------
      if (.not. module_is_initialized ) then
        call error_mesg ('longwave_fluxes_mod',   &
              'module has not been initialized', FATAL )
      endif

!---------------------------------------------------------------------
!
!---------------------------------------------------------------------
      ks =1
      ke = size(source2,3)-1
!DS   nbands = size(source,4)
      nbands = size(cld_trans,4)  ! DS change

!---------------------------------------------------------------------
!
!---------------------------------------------------------------------
      do m = 1, nbands
        do k=KS+1, KE+1
          flux_tmp(:,:,k) = source(:,:,KS,m)*trans(:,:,k,m    )
          flux_tmp2(:,:,k) = source2(:,:,k,m)*trans2(:,:,k ,m    )
        end do

!---------------------------------------------------------------------
!
!---------------------------------------------------------------------
        if ((m  == 1) .or. (m  >= 7)) then
          Lw_diagnostics%fluxn(:,:,KS,m) =    &
                                     Lw_diagnostics%fluxn(:,:,KS,m) + &
                                     source(:,:,KS,m)*trans(:,:,KS,m)
        else
          Lw_diagnostics%fluxn(:,:,KS,m) =    &
                                     Lw_diagnostics%fluxn(:,:,KS,m) + &
                                     source(:,:,KS,m)
        endif

!---------------------------------------------------------------------
!
!---------------------------------------------------------------------
        do k=KS+1,KE+1
          Lw_diagnostics%fluxn(:,:,k,m) =   &
                                     Lw_diagnostics%fluxn(:,:,k,m) +  &
                                     flux_tmp(:,:,k)* & 
                                     cld_trans(:,:,k,cld_ind(m))
          Lw_diagnostics%fluxn(:,:,KS,m) =  &
                                     Lw_diagnostics%fluxn(:,:,KS,m) + &
                                     flux_tmp2(:,:,k)*   &
                                     cld_trans(:,:,k,cld_ind(m))
        end do

!---------------------------------------------------------------------
!
!---------------------------------------------------------------------
        if (do_totcld_forcing) then
          if ((m  == 1) .or. (m  >= 7)) then
            Lw_diagnostics%fluxncf(:,:,KS,m) = source(:,:,KS,m)*   &
                                               trans(:,:,KS,m)
          else
            Lw_diagnostics%fluxncf(:,:,KS,m) =  source(:,:,KS,m)
          endif
          do k=KS+1,KE+1
            Lw_diagnostics%fluxncf(:,:,k,m) =   &
                                    Lw_diagnostics%fluxncf(:,:,k,m) + &
                                    flux_tmp(:,:,k)
            Lw_diagnostics%fluxncf(:,:,KS,m) =  &
                                    Lw_diagnostics%fluxncf(:,:,KS,m) + &
                                    flux_tmp2(:,:,k)
          end do
        endif
     end do   ! (m loop)

!---------------------------------------------------------------------


end subroutine longwave_fluxes_ks



!####################################################################
subroutine longwave_fluxes_k_down (klevel, source, trans, trans2,   &
                                   cld_trans, cld_ind, do_totcld_forcing, &
                                   Lw_diagnostics)

!---------------------------------------------------------------------
! subroutine to calculate longwave diagnostic fluxes
!---------------------------------------------------------------------

integer,                      intent(in)     ::  klevel
real,    dimension (:,:,:,:), intent(in)     ::  source, trans, &
                                                 trans2, cld_trans
type(lw_diagnostics_type),    intent(inout)  ::  Lw_diagnostics
integer, dimension(:),        intent(in)     ::  cld_ind
logical,                      intent(in)     :: do_totcld_forcing

!-------------------------------------------------------------------
!  intent(in) variables:
!
!     klevel  starting vertical level to calculate longwave fluxes
!     source  longwave flux source function
!     trans   longwave flux transmittance function
!     trans2  longwave flux transmittance function
!     cld_trans longwave cloud transmittance function
!     cld_ind   lookup table to translate longwave band index to cloud index
!
!  intent(inout) variables:
!
!     Lw_diagnostics contains the longwave diagnotics flux values
!
!---------------------------------------------------------------------

!------------------------------------------------------------------
!  local variables:

      real, dimension (size(source,1), size(source,2)) :: flux4, flux4a

      real    ::  flux_tmp, flux_tmp2
      integer ::  kp, i, j, israd, ierad, jsrad, jerad
      integer :: ke
      integer :: m, nbands

!---------------------------------------------------------------------
!  local variables:
!
!      flux4
!      flux4a
!      flux3a
!      flux_tmp
!      flux_tmp2
!      kp
!      i,j,k
!      nn
!      ntot
!      israd,ierad
!      jsrad,jerad
!      ke
!      m
!      nbands
!
!--------------------------------------------------------------------

!---------------------------------------------------------------------
!    be sure module has been initialized.
!---------------------------------------------------------------------
      if (.not. module_is_initialized ) then
        call error_mesg ('longwave_fluxes_mod',   &
              'module has not been initialized', FATAL )
      endif

!---------------------------------------------------------------------
!
!---------------------------------------------------------------------
      ierad  = size(source,1)
      jerad  = size(source,2)
      israd  = 1
      jsrad  = 1
      ke     = size(source,3)-1
!DS   nbands = size(trans,4)
      nbands = size(cld_trans,4)  ! DS change

!---------------------------------------------------------------------
!
!---------------------------------------------------------------------
      if (do_totcld_forcing) then
!---------------------------------------------------------------------
!
!---------------------------------------------------------------------
        do m=1, nbands
          do kp=klevel+1,KE+1
            do j=jsrad,jerad
              do i=israd,ierad   
                flux_tmp = source(i,j,klevel,m)*trans(i,j,kp,m)
                Lw_diagnostics%fluxn(i,j,kp,m) =    &
                            Lw_diagnostics%fluxn(i,j,kp,m) + flux_tmp*   &
                            cld_trans(i,j,kp, cld_ind(m))
                Lw_diagnostics%fluxncf(i,j,kp,m) =   &
                            Lw_diagnostics%fluxncf(i,j,kp,m) + flux_tmp
              end do
            end do
          end do

          flux4(:,:)  = 0.0
          flux4a(:,:) = 0.0
          do kp=klevel+1,KE+1
            do j=jsrad,jerad
              do i=israd,ierad   
                flux_tmp2 = source(i,j,kp,m)*trans2(i,j,kp,m)
                flux4(i,j) = flux4(i,j) + flux_tmp2*  &
                             cld_trans(i,j,kp, cld_ind(m))
                flux4a(i,j) = flux4a (i,j) + flux_tmp2
              end do
            end do
          end do

          do j=jsrad,jerad
            do i=israd,ierad   
              Lw_diagnostics%fluxn  (i,j,klevel,m) =   &
                                 Lw_diagnostics%fluxn  (i,j,klevel,m) +  &
                                 flux4(i,j)
              Lw_diagnostics%fluxncf(i,j,klevel,m) =  &
                                 Lw_diagnostics%fluxncf(i,j,klevel,m) +  &
                                 flux4a(i,j)
            end do
          end do
        end do  ! (nbands loop)
      else ! do_totcld_forcing
        do m=1, nbands
!---------------------------------------------------------------------
!
!---------------------------------------------------------------------
          do kp=klevel+1,KE+1
            do j=jsrad,jerad
              do i=israd,ierad   
                flux_tmp = source(i,j,klevel,m)*trans(i,j,kp,m)
                Lw_diagnostics%fluxn(i,j,kp,m) =    &
                            Lw_diagnostics%fluxn(i,j,kp,m) + flux_tmp*   &
                            cld_trans(i,j,kp, cld_ind(m))
              end do
            end do
          end do
!---------------------------------------------------------------------
!
!---------------------------------------------------------------------
          flux4(:,:)  = 0.0
          flux4a(:,:) = 0.0
          do kp=klevel+1,KE+1
            do j=jsrad,jerad
              do i=israd,ierad   
                flux_tmp2 = source(i,j,kp,m)*trans2(i,j,kp,m)
                flux4(i,j) = flux4(i,j) + flux_tmp2*  &
                             cld_trans(i,j,kp, cld_ind(m))
              end do
            end do
          end do
!---------------------------------------------------------------------
!
!---------------------------------------------------------------------
          do j=jsrad,jerad
            do i=israd,ierad   
              Lw_diagnostics%fluxn(i,j,klevel,m) =   &
                                 Lw_diagnostics%fluxn(i,j,klevel,m) +  &
                                 flux4(i,j)
            end do
          end do
        end do  ! (nbands loop)
      endif ! do_totcld_forcing

!---------------------------------------------------------------------

end subroutine longwave_fluxes_k_down


!####################################################################
subroutine longwave_fluxes_KE_KEp1 (source, trans, trans2, cld_trans,&
                                    cld_ind, do_totcld_forcing, Lw_diagnostics)

!---------------------------------------------------------------------
! subroutine to calculate longwave diagnostic fluxe
!---------------------------------------------------------------------

real,    dimension (:,:,:,:),   intent(in)    :: source, cld_trans
real,    dimension (:,:,:),     intent(in)    :: trans, trans2
integer, dimension(:),          intent(in)    :: cld_ind
logical,                        intent(in)    :: do_totcld_forcing
type(lw_diagnostics_type),      intent(inout) :: Lw_diagnostics

!-------------------------------------------------------------------
!  intent(in) variables:
!
!     source    longwave flux source function
!     cld_trans longwave flux transmittance function
!     trans     longwave flux transmittance function
!     trans2    longwave flux transmittance function
!     cld_ind   lookup table to translate longwave band index to cloud index
!
!  intent(inout) variables:
!
!     Lw_diagnostics contains the longwave diagnostic flux values
!
!---------------------------------------------------------------------

!--------------------------------------------------------------------
!  local variables:

      real, dimension (size(trans,1), size(trans,2)) ::  &
                                                   flux_tmp, flux_tmp2

      integer :: ke
      integer :: m, nbands

!---------------------------------------------------------------------
!  local variables:
!
!      flux_tmp
!      flux_tmp2
!      ke
!      m
!      nbands
!
!--------------------------------------------------------------------

!---------------------------------------------------------------------
!    be sure module has been initialized.
!---------------------------------------------------------------------
      if (.not. module_is_initialized ) then
        call error_mesg ('longwave_fluxes_mod',   &
              'module has not been initialized', FATAL )
      endif

!---------------------------------------------------------------------
!
!---------------------------------------------------------------------
      ke     = size(source,3) - 1
!DS   nbands = size(trans,3) 
      nbands = size(cld_trans,4)  ! DS change

!---------------------------------------------------------------------
!
!---------------------------------------------------------------------
      do m=1,nbands
        flux_tmp(:,:) = source(:,:,KE+1,m)*trans(:,:,m)
        flux_tmp2(:,:) = source(:,:,KE,m)*trans2(:,:,m)
        Lw_diagnostics%fluxn(:,:,KE,m) =    &
                       Lw_diagnostics%fluxn(:,:,KE,m) + &
                       flux_tmp(:,:)*cld_trans(:,:,KE+1,cld_ind(m))
        Lw_diagnostics%fluxn(:,:,KE+1,m) =  &
                       Lw_diagnostics%fluxn(:,:,KE+1,m) +   &
                       flux_tmp2(:,:)*cld_trans(:,:,KE+1,cld_ind(m))

!---------------------------------------------------------------------
!
!---------------------------------------------------------------------
        if (do_totcld_forcing) then
          Lw_diagnostics%fluxncf(:,:,KE,m) =   &
                       Lw_diagnostics%fluxncf(:,:,KE,m) + flux_tmp(:,:)
          Lw_diagnostics%fluxncf(:,:,KE+1,m) =  &
                       Lw_diagnostics%fluxncf(:,:,KE+1,m) +  &
                       flux_tmp2(:,:)
        endif
      end do  ! (nbands loop)

!---------------------------------------------------------------------

end subroutine longwave_fluxes_KE_KEp1



!####################################################################
subroutine longwave_fluxes_diag (source, trans, cld_trans, cld_ind, &
                                 do_totcld_forcing, Lw_diagnostics)

!---------------------------------------------------------------------
! subroutine to calculate longwave diagnotic fluxes
!---------------------------------------------------------------------
!---------------------------------------------------------------------
real, dimension (:,:,:,:), intent(in)    :: source, trans, cld_trans
integer, dimension(:),     intent(in)    :: cld_ind
logical,                   intent(in)    :: do_totcld_forcing
type(lw_diagnostics_type), intent(inout) :: Lw_diagnostics
!-------------------------------------------------------------------
!  intent(in) variables:
!
!     source     longwave flux source function
!     trans      longwave flux transmittance function
!     cld_trans  longwave cloud transmittance function
!     cld_ind    lookup table to translate longwave band index to cloud index
!
!  intent(inout) variables:
!
!     Lw_diagnostics contains the longwave diagnostics flux values
!
!---------------------------------------------------------------------

!---------------------------------------------------------------------
!   local variables

      real, dimension (size(trans,1), &
                       size(trans,2), &
                       size(trans,3)) ::   flux_tmp

      integer   :: k, ks, ke
      integer   :: m, nbands

!---------------------------------------------------------------------
!  local variables:
!
!      flux_tmp
!      k
!      ks
!      ke
!      m
!      nbands
!
!--------------------------------------------------------------------

!---------------------------------------------------------------------
!    be sure module has been initialized.
!---------------------------------------------------------------------
      if (.not. module_is_initialized ) then
        call error_mesg ('longwave_fluxes_mod',   &
              'module has not been initialized', FATAL )
      endif

!---------------------------------------------------------------------
!
!---------------------------------------------------------------------
      ks     = 1
      ke     = size(trans,3) - 1
!DS   nbands = size(trans,4)
      nbands = size(cld_trans,4)  ! DS change

!---------------------------------------------------------------------
!
!---------------------------------------------------------------------
      do m=1,nbands

!---------------------------------------------------------------------
!
!---------------------------------------------------------------------

        do k=KS+1, KE+1
          flux_tmp(:,:,k) = source(:,:,k,m)*trans(:,:,k,m)
        end do

!---------------------------------------------------------------------
!
!---------------------------------------------------------------------
        do k=KS+1,KE+1
          Lw_diagnostics%fluxn(:,:,k,m) =   &
                             Lw_diagnostics%fluxn(:,:,k,m) +  &
                             flux_tmp(:,:,k)*cld_trans(:,:,k,cld_ind(m))
        end do

!---------------------------------------------------------------------
!
!---------------------------------------------------------------------
        if (do_totcld_forcing) then
          do k=KS+1,KE+1
            Lw_diagnostics%fluxncf(:,:,k,m) =    &
                               Lw_diagnostics%fluxncf(:,:,k,m) +   &
                               flux_tmp(:,:,k)
          end do
        endif
      end do ! (m loop)

!---------------------------------------------------------------------



end subroutine longwave_fluxes_diag




!###################################################################
subroutine longwave_fluxes_sum (is, ie, js, je, flux, nbtrge,         &
                                Lw_diagnostics, fluxcf)

!---------------------------------------------------------------------
! subroutine to compute summation of diagnostic longwave fluxes over all bands
!---------------------------------------------------------------------

integer,                          intent(in)    :: is, ie, js, &
                                                   je, nbtrge
real, dimension(:,:,:),           intent(out)   :: flux
real, dimension(:,:,:), optional, intent(out)   :: fluxcf
type(lw_diagnostics_type),        intent(in)    :: Lw_diagnostics

!-------------------------------------------------------------------
!  intent(in) variables:
!
!     is,ie,js,je obsolete
!     nbtrge      number of longwave flux bands
!     Lw_diagnostics longwave flux diagnostics
!
!  intent(out) variables:
!
!     flux   all sky total longwave flux
!     fluxcf clear sky total longwave flux
!
!---------------------------------------------------------------------

!--------------------------------------------------------------------
!   local variables
!--------------------------------------------------------------------
    integer       ::   m
    logical       ::   do_totcld_forcing

!---------------------------------------------------------------------
!  local variables:
!
!      j,m        
!
!--------------------------------------------------------------------

!---------------------------------------------------------------------
!    be sure module has been initialized.
!---------------------------------------------------------------------
      if (.not. module_is_initialized ) then
        call error_mesg ('longwave_fluxes_mod',   &
              'module has not been initialized', FATAL )
      endif

!---------------------------------------------------------------------
!    are cloud forcing diagnostics needed?
!---------------------------------------------------------------------

      do_totcld_forcing = .false.
      if (present(fluxcf)) do_totcld_forcing = .true.

!---------------------------------------------------------------------
      flux = 0.
      do m= 1, 6+NBTRGE               
        flux(:,:,:) = flux(:,:,:) + Lw_diagnostics%fluxn(:,:,:,m)
      end do

!---------------------------------------------------------------------
!
!---------------------------------------------------------------------
      if (do_totcld_forcing) then 
        fluxcf = 0.
        do m= 1, 6+NBTRGE               
          fluxcf(:,:,:) = fluxcf(:,:,:) +    &
                          Lw_diagnostics%fluxncf(:,:,:,m)
        end do
      endif

!--------------------------------------------------------------------


end subroutine longwave_fluxes_sum


!#####################################################################

subroutine longwave_fluxes_end

!--------------------------------------------------------------------
!    longwave_fluxes_end is the destructor for the longwave_fluxes_mod.
!--------------------------------------------------------------------

!---------------------------------------------------------------------
!    be sure module has been initialized.
!---------------------------------------------------------------------
      if (.not. module_is_initialized ) then
        call error_mesg ('longwave_fluxes_mod',   &
              'module has not been initialized', FATAL )
      endif

!--------------------------------------------------------------------
!    mark the module as uninitialized.
!--------------------------------------------------------------------
      module_is_initialized = .false.

!--------------------------------------------------------------------


end subroutine longwave_fluxes_end


!#####################################################################


                end module longwave_fluxes_mod
