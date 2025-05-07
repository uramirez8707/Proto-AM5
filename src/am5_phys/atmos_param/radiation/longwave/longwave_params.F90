
                    module longwave_params_mod

! CONTACT: Fei.Liu@noaa.gov - fil
! REVIEWER: Dan.Schwarzkopf@noaa.gov - ds
!
!  This code has the number of bands defined for the longwave gases.

!   shared modules:

 use fms_mod, only: fms_init, &
                    error_mesg, &
                    FATAL

!--------------------------------------------------------------------

implicit none
private

!-------------------------------------------------------------------
!    longwave_params_mod defines basic parameters used by the
!    longwave radiation code.
!------------------------------------------------------------------


!---------------------------------------------------------------------
!----------- version number for this module -------------------

character(len=128)  :: version =  '$Id$'
character(len=128)  :: tagname =  '$Name$'


!--------------------------------------------------------------------
!----- interfaces ------

public     &
         longwave_params_init, &
         longwave_params_end

!private   &


!-------------------------------------------------------------------
!----- public data --------

!--------------------------------------------------------------------
!       NBCO215
!       NBLY_RSB
!       NBLY_CKD
!       NBLW
!       NBLX
!---------------------------------------------------------------------
integer, parameter, public   :: NBCO215     = 3
integer, parameter, public   :: NBCO210     = 1
integer, parameter, public   :: NBCO243     = 1
integer, parameter, public   :: NBLY_RSB    = 16
integer, parameter, public   :: NBLY_CKD    = 48
integer, parameter, public   :: NBLW        = 300
integer, parameter, public   :: NBLX        = 48




!-------------------------------------------------------------------
!----- private data --------

logical :: module_is_initialized = .false.  ! module is initialized ?


!---------------------------------------------------------------------
!---------------------------------------------------------------------



                          contains



!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!                     PUBLIC SUBROUTINES
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!####################################################################
subroutine longwave_params_init

!------------------------------------------------------------------
!    longwave_params_init is the constructor for longwave_params_mod.
!------------------------------------------------------------------

!------------------------------------------------------------------
!  local variables:

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

!----------------------------------------------------------------------
!    mark the module as initialized.
!----------------------------------------------------------------------
     module_is_initialized = .true.

!------------------------------------------------------------------
9000 format ( 'NBCO215=', i3,'  NBCO210=', i3,'  NBCO243=', i3,'  NBLY_RSB=', i4,   &
              '  NBLY_CKD=', i4, '  NBLW= ', i4, '  NBLX=', i4 )

!-------------------------------------------------------------------


end  subroutine longwave_params_init



!####################################################################
subroutine longwave_params_end

!-------------------------------------------------------------------
!    longwave_params_end is the destructor for longwave_params_mod.
!-------------------------------------------------------------------

!---------------------------------------------------------------------
!    be sure module has been initialized.
!---------------------------------------------------------------------
      if (.not. module_is_initialized ) then
        call error_mesg ('longwave_params_mod',   &
              'module has not been initialized', FATAL )
      endif

!---------------------------------------------------------------------
!    mark the module as uninitialized.
!---------------------------------------------------------------------
      module_is_initialized = .false.

!-------------------------------------------------------------------



end subroutine longwave_params_end


!####################################################################


                   end module longwave_params_mod
