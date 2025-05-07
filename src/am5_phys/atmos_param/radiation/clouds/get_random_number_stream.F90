
module random_number_streams_mod

!----------------------------------------------------------------------

use mpp_mod,            only:  input_nml_file
use fms_mod,            only:  mpp_pe, &
                               mpp_root_pe, stdlog, fms_init, &
                               write_version_number, &
                               check_nml_error, &
                               error_mesg, FATAL, NOTE
use time_manager_mod,   only:  time_type
use constants_mod,      only:  RADIAN

use random_numbers_mod, only:  randomNumberStream,   &
                               initializeRandomNumberStream, &
                               constructSeed

use cloudrad_types_mod, only:  cloudrad_control_type

!--------------------------------------------------------------------

implicit none 
private

!--------------------------------------------------------------------
!----------- version number for this module --------------------------

character(len=128)  :: version =  '$Id$'
character(len=128)  :: tagname =  '$Name$'
!---------------------------------------------------------------------
!-------  interfaces --------

public  random_number_streams_init, &
        random_number_streams_end, &
        get_random_number_streams

!---------------------------------------------------------------------
!-------- namelist  ---------

logical :: force_use_of_temp_for_seed = .true.
                                        ! if true, when using stochastic 
                                        ! clouds, force the seed to use 
                                        ! top-model-level temps as input to
                                        ! random number generator
namelist /random_number_streams_nml/ force_use_of_temp_for_seed

!---------------------------------------------------------------------
!----  private data -------

logical :: module_is_initialized = .false.

!----------------------------------------------------------------------
!   variables needed for legacy random number seed:
!----------------------------------------------------------------------
real, dimension(:,:), allocatable  :: lats, lons ! lat and lon of columns
                                                 ! in this processor's
                                                 ! domain [ degrees ]
logical :: use_temp_for_seed


CONTAINS

!######################################################################

subroutine random_number_streams_init ( lonb, latb, Cldrad_control )
real, dimension(:,:),        intent(in)    ::  lonb, latb
type(cloudrad_control_type), intent(inout) ::  Cldrad_control

!----------------------------------------------------------------------
!   local variables:

      integer  ::   ierr, io
      integer  ::   id, jd, i, j, ii, jj

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
!    read namelist.
      read (input_nml_file, nml=random_number_streams_nml, iostat=io)
      ierr = check_nml_error(io,"random_number_streams_nml")
!----------------------------------------------------------------------
!    write version number and namelist to logfile.
!---------------------------------------------------------------------
      call write_version_number (version, tagname)
      if (mpp_pe() == mpp_root_pe() ) &
           write (stdlog(), nml=random_number_streams_nml)

!---------------------------------------------------------------------
!     determine the source of the random number seed generator to be used
!     for the stochastic cloud generation.
!---------------------------------------------------------------------
      use_temp_for_seed = .true.

!---------------------------------------------------------------------
!    mark the module initialized.
!---------------------------------------------------------------------
      module_is_initialized = .true.

!--------------------------------------------------------------------

end subroutine random_number_streams_init
  
!######################################################################

subroutine get_random_number_streams ( is, js, Rad_time, temp, streams, perm )
integer,                  intent(in)                  :: is, js
type(time_type),          intent(in)                  :: Rad_time
real,                     intent(in),  dimension(:,:) :: temp
type(randomNumberStream), intent(out), dimension(:,:) :: streams
integer,                  intent(in),  optional       :: perm

integer :: i, j
real    :: seedwts(8) = (/3000.,1000.,300.,100.,30.,10.,3.,1./)

      if (use_temp_for_seed) then
        do j = 1, size(streams,2)
          do i = 1, size(streams,1)
            streams(i,j) = initializeRandomNumberStream(  &
                               ishftc(nint(temp(i,j)*seedwts),perm))

          enddo
        enddo
      else
        do j = 1, size(streams,2)
          do i = 1, size(streams,1)
            streams(i,j) = initializeRandomNumberStream(  &
                            constructSeed(nint(lons(is + i - 1, js + j - 1)), &
                                          nint(lats(is + i - 1, js + j - 1)), &
                                                    Rad_time, perm=perm))
          enddo
        enddo
      endif

!----------------------------------------------------------------------

end subroutine get_random_number_streams

!######################################################################

subroutine random_number_streams_end

!----------------------------------------------------------------------
!    skip if the module is not initialized.
!----------------------------------------------------------------------
      if (.not.module_is_initialized) return

!----------------------------------------------------------------------
!    deallocate module arrays
!----------------------------------------------------------------------
      if (allocated(lats)) deallocate(lats)
      if (allocated(lons)) deallocate(lons)

!----------------------------------------------------------------------
!    mark the module as no longer initialized.
!----------------------------------------------------------------------
      module_is_initialized = .false.

!----------------------------------------------------------------------

end subroutine random_number_streams_end

!######################################################################

end module random_number_streams_mod

