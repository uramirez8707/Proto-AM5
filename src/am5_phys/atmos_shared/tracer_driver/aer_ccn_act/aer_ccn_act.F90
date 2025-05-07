        module aer_ccn_act_mod

use fms_mod,             only: error_mesg, FATAL, &
                               mpp_pe, mpp_root_pe, stdlog, &
                               write_version_number, &
                               check_nml_error
use fms2_io_mod,         only: file_exists, ascii_read
use mpp_mod,             only: input_nml_file, get_unit
use aer_ccn_act_k_mod,   only: aer_ccn_act_k, aer_ccn_act2_k, &
                               aer_ccn_act_wpdf_k, aer_ccn_act_k_init, &
                               aer_ccn_act_k_end, aer_ccn_act_wpdf_m_k

implicit none
private
    private Loading
      
    public aer_ccn_act, aer_ccn_act2, aer_ccn_act_wpdf, &
           aer_ccn_act_wpdf_m, aer_ccn_act_init, aer_ccn_act_end

!--------------------- version number ---------------------------------

character(len=128) :: version = '$Id$'
character(len=128) :: tagname = '$Name$'

!---------------- private data -------------------


!-----------------------------------------------------------------------
!-------------------- namelist -----------------------------------------

logical  :: nooc = .false.   ! include organic aerosols as ccns ?
real     :: sul_concen = 0.1
real     :: low_concen = 0.1
real     :: high_concen = 1.
 !Parameters for look-up tables
 
real ::  lowup=0.3 !m/s
real ::  highup=10.

! earlier values: lowup2 = 0.001, highmass2 = 1000., highmass3 = 1000.
!real ::  lowup2=0.0001 !m/s
real ::  lowup2=0.01   !m/s
real ::  highup2=0.3
real ::  lowmass2=0.01 !ug m-3
!real ::  highmass2=1000.
real ::  highmass2=100.
real ::  lowmass3=0.01 !ug m-3
!real ::  highmass3=1000.
real ::  highmass3=100.
real ::  lowmass4=0.01 !ug m-3
real ::  highmass4=100.
real ::  lowmass5=0.01 !ug m-3
real ::  highmass5=100.
real :: lowT2=243.15 !K
real :: highT2=308.15

namelist /aer_ccn_act_nml/ nooc, sul_concen, low_concen, high_concen, &
                           lowup, highup, lowup2, highup2, lowmass2, &
                           highmass2, lowmass3, highmass3,  &
                           lowmass4, highmass4, lowmass5, highmass5, &
                           lowT2, highT2


logical :: module_is_initialized  = .false.
 
contains

subroutine aer_ccn_act (T1, P1, Updraft1, TotalMass, Drop)
real, dimension(:), intent(inout) :: TotalMass
real, intent(in) :: T1, P1, Updraft1
real, intent(inout) :: Drop
    
  integer :: tym, ier
  character(len=256) :: ermesg

  if(.not. module_is_initialized) call aer_ccn_act_init()

  tym = size (totalmass,1)

  call aer_ccn_act_k (T1, P1, Updraft1, TotalMass, tym, Drop, ier,  &
                      ermesg)
  if (ier /= 0) call error_mesg ('aer_ccn_act', ermesg, FATAL)

  
end subroutine aer_ccn_act

subroutine aer_ccn_act2 (T1, P1, Updraft1, TotalMass, mu,airdens,Nc,qc,qt,qe,tc,te,Drop)

!T1 temperature (K)
!P1 pressure (Pa)
!Updraft1 updraft velocity (m/s)
!TotalMass aerosol mass ()
!mu entrainment coef. (/s)
!airdens air density (kg/m3 air)
!Nc droplet mixing ratio (#/kg air)
!qc in-cloud vapor mixing ratio (kg water/kg air)
!qt in-cloud total water mixing ratio qc + ql (kg water/kg air)
!qe environment vapor mixing ratio (kg water/kg air)
!tc in-cloud temperature (K)
!te environment temperature (K)
!Drop droplet number concentration (#/cc)

real, dimension(:), intent(in) :: TotalMass
real, intent(in) :: T1, P1, Updraft1, mu,airdens, Nc, qc, qt, qe, tc, te
real, intent(inout) :: Drop

  integer :: tym, ier
  character(len=256) :: ermesg
        
  if(.not. module_is_initialized) call aer_ccn_act_init()
  tym = size (totalmass,1)

  call aer_ccn_act2_k (T1, P1, Updraft1, TotalMass, tym, mu,  &
                       airdens,Nc,qc,qt,qe,tc,te,Drop, ier, ermesg)
  if (ier /= 0) call error_mesg ('aer_ccn_act2', ermesg, FATAL)

end subroutine aer_ccn_act2

!-->cjg: addition
!
! Additional subroutines to compute CCN activation by integrating
! over an assumed subgrid-scale PDF of w

subroutine aer_ccn_act_wpdf(T, p, wm, wp2, totalmass, drop)

! Compute CCN activation assuming a normal distribution of w
! given by its mean (wm) and second moment (wp2)

real, intent(in)    :: T, p, wm, wp2
real, intent(inout) :: totalmass(4)
real, intent(out)   :: drop

  integer :: tym, ier
  character(len=256) :: ermesg

  if(.not. module_is_initialized) call aer_ccn_act_init()
  tym = size (totalmass,1)

   call aer_ccn_act_wpdf_k (T, p, wm, wp2, totalmass, tym,           &
                            drop, ier, ermesg)
  if (ier /= 0) call error_mesg ('aer_ccn_act_wpdf', ermesg, FATAL)

end subroutine aer_ccn_act_wpdf

!----------------------------------------------------------------------

subroutine aer_ccn_act_wpdf_m(T, p, wm, wp2, offs, totalmass, drop)


! Compute CCN activation assuming a normal distribution of w
! given by its mean (wm) and second moment (wp2)

real, intent(in)    :: T, p, wm, wp2
integer, intent(in) :: offs
real, intent(inout) :: totalmass(4)
real, intent(out)   :: drop

  integer :: tym, ier
  character(len=256) :: ermesg

  if(.not. module_is_initialized) call aer_ccn_act_init()
  tym = size (totalmass,1)

  call aer_ccn_act_wpdf_m_k (T, p, wm, wp2, offs, totalmass, tym,       &
                             drop, ier, ermesg)
  if (ier /= 0) call error_mesg ('aer_ccn_act_wpdf_m', ermesg, FATAL)

end subroutine aer_ccn_act_wpdf_m

!------------------------------------------------------------------------

subroutine aer_ccn_act_init ()

!--------------------------------------------------------------------  
!  local variables:
      
      integer   ::   unit, ierr, io, logunit
      integer, parameter :: res = 20 !
      real, dimension(res,res,res,res,res) :: droplets

      integer, parameter :: res2 = 20 !
      real, dimension(res2,res2,res2,res2,res2) :: droplets2

      if (module_is_initialized) return

!--------------------------------------------------------------------- 
!    read namelist.
!--------------------------------------------------------------------
      if ( file_exists('input.nml')) then
        read (input_nml_file, nml=aer_ccn_act_nml, iostat=io)
        ierr = check_nml_error(io,'aer_ccn_act_nmliostat=io')
      endif
!---------------------------------------------------------------------
!    write version number and namelist to logfile.
!--------------------------------------------------------------------
       call write_version_number (version, tagname)
       logunit=stdlog()
       if (mpp_pe() == mpp_root_pe() ) &
                        write (logunit, nml=aer_ccn_act_nml)

       call Loading( droplets, droplets2)

       call aer_ccn_act_k_init (droplets,   &
                      droplets2, res, res2, nooc,  &
                       sul_concen, low_concen, high_concen, &
                       lowup, highup, lowup2, highup2, lowmass2, &
                       highmass2, lowmass3, highmass3,  &
                       lowmass4, highmass4, lowmass5, highmass5, &
                      lowT2, highT2  )
       module_is_initialized  = .true.

end subroutine aer_ccn_act_init


subroutine Loading(droplets, droplets2)

real, dimension(:,:,:,:,:), intent(out) :: droplets, droplets2

character(len=:), dimension(:), allocatable :: droplets_file !< Restart file saved as a string
integer :: ios=0

  call ascii_read('INPUT/droplets.dat', droplets_file)
  read(droplets_file,fmt=*,iostat=ios) droplets
  if (ios /= 0) call error_mesg ('aer_ccn_act_init', 'Read of INPUT/droplets.dat failed', FATAL)
  deallocate(droplets_file)

  call ascii_read('INPUT/droplets2.dat', droplets_file)
  read(droplets_file,fmt=*,iostat=ios) droplets2
  if (ios /= 0) call error_mesg ('aer_ccn_act_init', 'Read of INPUT/droplets2.dat failed', FATAL)
  deallocate(droplets_file)

end subroutine Loading


subroutine aer_ccn_act_end()

  call aer_ccn_act_k_end 
  module_is_initialized  = .false.

end subroutine aer_ccn_act_end

end module aer_ccn_act_mod
