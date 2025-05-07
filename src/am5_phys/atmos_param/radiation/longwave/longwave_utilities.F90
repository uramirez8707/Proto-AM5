 
module longwave_utilities_mod


use fms_mod,      only: write_version_number, &
                        error_mesg, FATAL
!--------------------------------------------------------------------

implicit none
private

!---------------------------------------------------------------------
!-------  version number --------

character(len=128)  :: version =  '$Id$'
character(len=128)  :: tagname =  '$Name$'

!---------------------------------------------------------------------
!-------  interfaces --------

public :: longwave_utilities_init, &
          longwave_utilities_end, &
          locate_in_table,  &
          looktab, table_alloc

interface looktab
    module procedure  looktab_type1, looktab_type2, looktab_type3
end interface

interface table_alloc
   module procedure    table1_alloc, table2_alloc, table3_alloc
end interface

!---------------------------------------------------------------------
!------- public derived-types ------

public :: gas_tf_type

type gas_tf_type
     real, dimension(:,:,:),   allocatable :: tdav,    &
                                          tlsqu,   &
                                          tmpdiff, &
                                          tstdav,  &
                                          n2o9c,   &
                                          tn2o17
     real, dimension(:,:,:),   allocatable :: co2nbl
     real, dimension(:,:,:),   allocatable :: co2990nbl, &
                                          co2900nbl, &
                                          co21070nbl
     real, dimension(:,:,:,:), allocatable :: co2spnb
     real, dimension(:,:,:),   allocatable :: co2990spnb
     real, dimension(:,:,:),   allocatable :: co2900spnb
     real, dimension(:,:,:),   allocatable :: co21070spnb
     real, dimension(:,:),     allocatable :: a1,    &
                                          a2
end type gas_tf_type

!------------------------------------------------------------------

public longwave_tables1_type

type longwave_tables1_type
    real, dimension(:,:), allocatable  ::  vae,   &
                                       td, &
                                       md, &
                                       cd
end type longwave_tables1_type

!--------------------------------------------------------------------

public longwave_tables2_type

type longwave_tables2_type
    real, dimension(:,:,:), allocatable  ::  vae,  &
                                         td,  &
                                         md,   &
                                         cd
end type longwave_tables2_type

!---------------------------------------------------------------------

public longwave_tables3_type

type longwave_tables3_type
     real,  dimension(:,:), allocatable    ::  vae,   &
                                           td
end type longwave_tables3_type

!---------------------------------------------------------------------

public lw_clouds_type

type lw_clouds_type
     real, dimension(:,:,:,:),   allocatable :: taucld_rndlw, &
                                            taucld_mxolw, &
                                            taunbl_mxolw
end type lw_clouds_type

!------------------------------------------------------------------

public lw_table_type

type lw_table_type
     real, dimension(:),    allocatable :: bdlocm,   &
                                       bdhicm,  &
                                       bandlo,  &
                                       bandhi
     integer, dimension(:), allocatable :: iband
end type lw_table_type

!------------------------------------------------------------------

public optical_path_type

type optical_path_type
     real, dimension (:,:,:,:), allocatable :: empl1f,  &
                                           empl2f,  &
                                           vrpfh2o, &
                                           xch2obd,  &
                                           tphfh2o, &
                                           avephif, &
                                           totaerooptdep
     real, dimension (:,:,:),   allocatable :: empl1, &
                                           empl2,  &
                                           var1, &
                                           var2, &
                                           emx1f,   &
                                           emx2f,   &
                                           totvo2,  &
                                           avephi,&
                                           totch2obdwd, &
                                           xch2obdwd, &
                                           totphi,   &
                                           cntval, &
                                           toto3,   &
                                           tphio3,  &
                                           var3,  &
                                           var4,        &
                                           wk,         &
                                           rh2os,  &
                                           rfrgn,  &
                                           tfac, &
                                           totaerooptdep_15, &
                                           totf11,   &
                                           totf12,  &
                                           totf113,   &
                                           totf22
      real, dimension (:,:), allocatable    :: emx1,  &
                                           emx2,  &
                                           csfah2o, &
                                           aerooptdep_KE_15
end type optical_path_type

!------------------------------------------------------------------

public sealw99_control_type

type sealw99_control_type
    character(len=16) :: continuum_form
    character(len=16) :: linecatalog_form
    logical           :: do_ch4lbltmpint
    logical           :: do_n2olbltmpint
    logical           :: do_lwcldemiss
    logical           :: do_h2o
    logical           :: do_o3
    logical           :: do_ch4
    logical           :: do_n2o
    logical           :: do_co2
    logical           :: do_co2_10um
    logical           :: do_cfc
end type sealw99_control_type

!------------------------------------------------------------------

public table_axis_type

type table_axis_type
  integer :: first_col
  real    :: min_val
  real    :: max_val
  real    :: tab_inc
end type table_axis_type

!---------------------------------------------------------------------
!------- public data -------

type (table_axis_type),        public   ::    &
               temp_1 = table_axis_type( 1, 100.0, 370.0, 10.0), &
               mass_1 = table_axis_type( 1, -16.0,   1.9,  0.1)

type (sealw99_control_type),  public   :: &
      Sealw99_control = sealw99_control_type(  '    ', '    ', &
                                              .false., .false., .false., &
                                              .false., .false., .false., &
                                              .false., .false., .false., &
                                              .false. )

!---------------------------------------------------------------------
!------- private data -------

logical :: module_is_initialized=.false.   ! module is initialized ?

!---------------------------------------------------------------------

CONTAINS

!#####################################################################
!
!                     PUBLIC SUBROUTINES
!
!#####################################################################

subroutine longwave_utilities_init

!---------------------------------------------------------------------
!    if routine has already been executed, exit.
!---------------------------------------------------------------------
      if (module_is_initialized) return

!---------------------------------------------------------------------
!    write version number to logfile.
!---------------------------------------------------------------------
      call write_version_number (version, tagname)

!-------------------------------------------------------------------
!    mark the module as initialized.
!-------------------------------------------------------------------
      module_is_initialized = .true.

!------------------------------------------------------------------

end subroutine longwave_utilities_init

!#####################################################################

subroutine longwave_utilities_end

!--------------------------------------------------------------------
!    this is the destructor for longwave_utilities_mod
!---------------------------------------------------------------------

!---------------------------------------------------------------------
!    be sure module has been initialized.
!---------------------------------------------------------------------
      if (.not. module_is_initialized ) then
        call error_mesg ('longwave_utilites_mod',   &
             'module has not been initialized', FATAL )
      endif

!---------------------------------------------------------------------
!    mark the module as uninitialized.
!---------------------------------------------------------------------
      module_is_initialized = .false.

!-------------------------------------------------------------------

end subroutine longwave_utilities_end

!#####################################################################
subroutine locate_in_table (table_axis, x, dx, ix, k_min, k_max)

!---------------------------------------------------------------------
!    given array x and an arithmetic sequence of table column headings
!    tabxmin, tabxmin+tabdeltax, ..., corresponding to column ixlow, 
!    ixlow+1, ..., ixupp, locate_in_table returns the array ix of
!    column indices and the array dx of residuals.
!    author: c. h. goldberg
!    revised: 1/1/93
!    certified:  radiation version 1.0
!----------------------------------------------------------------------

type(table_axis_type),     intent(in)  :: table_axis
real,    dimension(:,:,:), intent(in)  :: x
integer,                   intent(in)  :: k_min, k_max
real,    dimension(:,:,:), intent(out) :: dx
integer, dimension(:,:,:), intent(out) :: ix

!--------------------------------------------------------------------
!  intent(in) variables:
!
!    table_axis  contains the axis information such a min, increment, and first column value
!    x  array from which data is to be searched
!    k_min minimum k value of the search domain
!    k_max maximum k value of the search domain
!
!  intent(out) variables:
!
!    dx residual between x and x(ix+first_column)
!    ix index value of the searched domain in the array
!
!---------------------------------------------------------------------

!---------------------------------------------------------------------
!  local variables:

      real, dimension (size(x,1), size(x,2), size(x,3))  ::  fx
      integer     ::  k

!---------------------------------------------------------------------
!  local variables:
!
!     fx
!     table_min
!     table_inc
!     k
!     table_col
!
!---------------------------------------------------------------------

!---------------------------------------------------------------------
!    be sure module has been initialized.
!---------------------------------------------------------------------
      if (.not. module_is_initialized ) then
        call error_mesg ('longwave_utilities_mod',   &
               'module has not been initialized', FATAL )
      endif

!---------------------------------------------------------------------
!
!---------------------------------------------------------------------
      do k=k_min,k_max
        fx (:,:,k) = AINT((x(:,:,k) - table_axis%min_val )/  &
                     table_axis%tab_inc)
        ix (:,:,k) = INT(fx(:,:,k)) + table_axis%first_col
        dx (:,:,k) = x(:,:,k) - fx(:,:,k)*table_axis%tab_inc - &
                     table_axis%min_val
      end do

!---------------------------------------------------------------------

end subroutine locate_in_table

!####################################################################
subroutine looktab_type1 (tab, ix, iy, dx, dy, answer, k_min, k_max)

!----------------------------------------------------------------------
!    given arrays ix(:,:,:) and iy(:,:,:) of integral subscripts and
!    arrays dx(:,:,:) and dy(:,:,:) of differences from x(:,:,:) and
!    y(:,:,:), calculate answer(:,:,:) = f(x(:,:,:), y(:,:,:))
!    from four tables of values, f, df/dx, df/dy, and d2f/dxdy.
!    author: c. h. goldberg
!    revised: 1/1/93
!    certified:  radiation version 1.0
!--------------------------------------------------------------------

type(longwave_tables1_type), intent(in)  :: tab
integer,dimension(:,:,:),    intent(in)  :: ix, iy
real,   dimension(:,:,:),    intent(in)  :: dx, dy
real,   dimension(:,:,:),    intent(out) :: answer
integer,                     intent(in)  :: k_min, k_max

!---------------------------------------------------------------------
!  intent(in) variables:
!
!    tab data array that contains function values and differentials
!    ix  x subscript of input data array
!    iy  y subscript of input data array
!    dx  x step in the x subscript space
!    dy  y step in the y subscript space
!    k_min the minimum k value of the domain
!    k_max the maximum k value of the domain
!
!  intent(out) variables:
!
!    answer  the answer to be calculated
!
!---------------------------------------------------------------------

!--------------------------------------------------------------------
!  local variables:

      integer    ::  i_min, i_max, j_min, j_max, i, j, k

!--------------------------------------------------------------------
!  local variables:
!
!    i_min
!    i_max
!    j_min
!    j_max
!    i,j,k
!
!--------------------------------------------------------------------

!---------------------------------------------------------------------
!    be sure module has been initialized.
!---------------------------------------------------------------------
      if (.not. module_is_initialized ) then
        call error_mesg ('longwave_utilities_mod',   &
               'module has not been initialized', FATAL )
      endif

      i_min = lbound(ix,1)
      i_max = ubound(ix,1)
      j_min = lbound(ix,2)
      j_max = ubound(ix,2)

      do k=k_min, k_max
        do j=j_min, j_max
          do i=i_min, i_max
            answer(i,j,k) =                                         &
                                      tab%vae (ix(i,j,k), iy(i,j,k)) + &
                            dx(i,j,k)*tab%td  (ix(i,j,k), iy(i,j,k)) + &
                            dy(i,j,k)*tab%md  (ix(i,j,k), iy(i,j,k)) + &
                  dx(i,j,k)*dy(i,j,k)*tab%cd(ix(i,j,k), iy(i,j,k))
          end do
        end do
      end do

!---------------------------------------------------------------------

end subroutine looktab_type1

!#####################################################################
subroutine looktab_type2 (tab, ix, iy, dx, dy, answer, k_min, k_max, m)

!-------------------------------------------------------------------
!    given arrays ix(:,:,:) and iy(:,:,:) of integral subscripts and
!    arrays dx(:,:,:) and dy(:,:,:) of differences from x(:,:,:) and
!    y(:,:,:), calculate answer(:,:,:) = f(x(:,:,:), y(:,:,:))
!    from four tables of values, f, df/dx, df/dy, and d2f/dxdy.
!    author: c. h. goldberg
!    revised: 1/1/93
!    certified:  radiation version 1.0
!--------------------------------------------------------------------

type(longwave_tables2_type), intent(in)   :: tab
integer, dimension (:,:,:),  intent(in)   :: ix, iy
integer,                     intent(in)   :: m
real, dimension (:,:,:),     intent(in)   :: dx, dy
real, dimension (:,:,:),     intent(out)  :: answer
integer,                     intent(in)   :: k_min, k_max

!---------------------------------------------------------------------
!  intent(in) variables:
!
!    tab  the data array that contains function values and differentials
!    ix  x subscript of input data array
!    iy  y subscript of input data array
!    m   z index of the differential arrays
!    dx  x step in the x subscript space
!    dy  y step in the y subscript space
!    k_min the minimum k value of the domain
!    k_max the maximum k value of the domain
!
!  intent(out) variables:
!
!    answer
!
!---------------------------------------------------------------------

!--------------------------------------------------------------------
!  local variables:

       integer    ::    i_min, i_max, j_min, j_max
       integer    ::    i, j, k

!--------------------------------------------------------------------
!  local variables:
!
!    i_min
!    i_max
!    j_min
!    j_max
!    i,j,k
!
!-------------------------------------------------------------------

!---------------------------------------------------------------------
!    be sure module has been initialized.
!---------------------------------------------------------------------
      if (.not. module_is_initialized ) then
        call error_mesg ('longwave_utilities_mod',   &
               'module has not been initialized', FATAL )
      endif

      i_min = lbound(ix,1)
      i_max = ubound(ix,1)
      j_min = lbound(ix,2)
      j_max = ubound(ix,2)

      do k=k_min, k_max
        do j=j_min, j_max
          do i=i_min, i_max
            answer(i,j,k) =                                           &
                                   tab%vae (ix(i,j,k), iy(i,j,k),m) + &
                         dx(i,j,k)*tab%td (ix(i,j,k), iy(i,j,k),m) + &
                         dy(i,j,k)*tab%md (ix(i,j,k), iy(i,j,k),m) + &
               dx(i,j,k)*dy(i,j,k)*tab%cd   (ix(i,j,k), iy(i,j,k),m)
           end do
        end do
      end do

!--------------------------------------------------------------------

end subroutine looktab_type2

!###################################################################
subroutine looktab_type3 (tab, ix, dx,  answer, k_min, k_max, n)

!----------------------------------------------------------------------
!
!    given arrays ix(:,:,:) and dx(:,:,:) of integer subscripts and!
!    differences from x(:,:,:) and constant column subscript iyconst, 
!    calculate answer(:,:,:) = f(x(:,:,:), y(:,:,:)) from four tables
!    of values f, df/dx, df/dy, and d2f/dxdy.
!    author: c. h. goldberg
!    revised: 1/1/93
!    certified:  radiation version 1.0
!-----------------------------------------------------------------------

type(longwave_tables3_type), intent(in)  :: tab
integer, dimension (:,:,:),  intent(in)  :: ix
integer,                     intent(in)  :: n
real,    dimension(:,:,:),   intent(in)  :: dx
real,    dimension(:,:,:),   intent(out) :: answer
integer,                     intent(in)  :: k_min, k_max

!---------------------------------------------------------------------
!  intent(in) variables:
!
!    tab  the data array that contains function values and differentials
!    ix  x subscript of input data array
!    iy  y subscript of input data array
!    m   z index of the differential arrays
!    dx  x step in the x subscript space
!    dy  y step in the y subscript space
!    k_min the minimum k value of the domain
!    k_max the maximum k value of the domain
!
!  intent(out) variables:
!
!    answer
!
!---------------------------------------------------------------------

!--------------------------------------------------------------------
!  local variables:

      integer    :: i_min, i_max, j_min, j_max
      integer    :: i, j, k

!--------------------------------------------------------------------
!  local variables:
!
!    i_min
!    i_max
!    j_min
!    j_max
!    i,j,k
!
!-------------------------------------------------------------------

!---------------------------------------------------------------------
!    be sure module has been initialized.
!---------------------------------------------------------------------
      if (.not. module_is_initialized ) then
        call error_mesg ('longwave_utilities_mod',   &
               'module has not been initialized', FATAL )
      endif

!-----------------------------------------------------------------
      i_min = lbound(ix,1)
      i_max = ubound(ix,1)
      j_min = lbound(ix,2)
      j_max = ubound(ix,2)

      do k=k_min, k_max
        do j=j_min, j_max
          do i=i_min, i_max
                answer(i,j,k) =                                 &
                                      tab%vae (ix(i,j,k),n) +   &
                            dx(i,j,k)*tab%td(ix(i,j,k),n)
          end do
        end do
      end do

!------------------------------------------------------------------

end subroutine  looktab_type3

!#####################################################################
subroutine table1_alloc (tab, dim1, dim2)

!------------------------------------------------------------------
!    table1_alloc allocates the arrays contained in a 
!    longwave_tables1_type variable.
!------------------------------------------------------------------

type(longwave_tables1_type), intent (inout) :: tab
integer,                     intent(in)     :: dim1, dim2

!-------------------------------------------------------------------
! intent(in) variables:
!
!     dim1 size of the x dimension
!     dim2 size of the y dimension
!
!  intent(inout) variables:
!
!     tab the longwave table
!
!---------------------------------------------------------------------

!---------------------------------------------------------------------
!    be sure module has been initialized.
!---------------------------------------------------------------------
      if (.not. module_is_initialized ) then
        call error_mesg ('longwave_utilities_mod',   &
               'module has not been initialized', FATAL )
      endif

!---------------------------------------------------------------------
!    allocate the component arrays of a longwave_tables1_type variable.
!---------------------------------------------------------------------
      allocate (tab%vae(dim1, dim2))
      allocate (tab%td (dim1, dim2))
      allocate (tab%md (dim1, dim2))
      allocate (tab%cd (dim1, dim2))

!---------------------------------------------------------------------

end subroutine table1_alloc

!####################################################################
subroutine table2_alloc (tab, dim1, dim2, dim3)

!------------------------------------------------------------------
!    table2_alloc allocates the arrays contained in a 
!    longwave_tables2_type variable.
!------------------------------------------------------------------

type(longwave_tables2_type), intent (inout) :: tab
integer,                     intent(in)     :: dim1, dim2, dim3

!-------------------------------------------------------------------
! intent(in) variables:
!
!     dim1  size of the x dimension
!     dim2  size of the y dimension
!     dim3  size of the z dimension
!
!  intent(inout) variables:
! 
!     tab the longwave tables
!
!---------------------------------------------------------------------

!---------------------------------------------------------------------
!    be sure module has been initialized.
!---------------------------------------------------------------------
      if (.not. module_is_initialized ) then
        call error_mesg ('longwave_utilities_mod',   &
               'module has not been initialized', FATAL )
      endif

!---------------------------------------------------------------------
!    allocate the component arrays of a longwave_tables2_type variable.
!---------------------------------------------------------------------
      allocate (tab%vae(dim1, dim2, dim3))
      allocate (tab%td (dim1, dim2, dim3))
      allocate (tab%md (dim1, dim2, dim3))
      allocate (tab%cd (dim1, dim2, dim3))

!--------------------------------------------------------------------

end subroutine table2_alloc


!#####################################################################
subroutine table3_alloc (tab, dim1, dim2)

type(longwave_tables3_type), intent (inout) :: tab
integer,                     intent(in)     :: dim1, dim2

!-------------------------------------------------------------------
! intent(in) variables:
!
!     dim1 size of the x dimension
!     dim2 size of the y dimension
!
!  intent(inout) variables:
! 
!     tab the longwave tables
!
!---------------------------------------------------------------------

!---------------------------------------------------------------------
!    be sure module has been initialized.
!---------------------------------------------------------------------
      if (.not. module_is_initialized ) then
        call error_mesg ('longwave_utilities_mod',   &
               'module has not been initialized', FATAL )
      endif

!---------------------------------------------------------------------
!    allocate the component arrays of a longwave_tables3_type variable.
!---------------------------------------------------------------------
      allocate (tab%vae(dim1, dim2))
      allocate (tab%td (dim1, dim2))

end subroutine table3_alloc

!##################################################################

end module longwave_utilities_mod

