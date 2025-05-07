#***********************************************************************
#*                   GNU Lesser General Public License
#*
#* This file is part of the GFDL Flexible Modeling System (FMS).
#*
#* FMS is free software: you can redistribute it and/or modify it under
#* the terms of the GNU Lesser General Public License as published by
#* the Free Software Foundation, either version 3 of the License, or (at
#* your option) any later version.
#*
#* FMS is distributed in the hope that it will be useful, but WITHOUT
#* ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
#* FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
#* for more details.
#*
#* You should have received a copy of the GNU Lesser General Public
#* License along with FMS.  If not, see <http://www.gnu.org/licenses/>.
#***********************************************************************

AC_PREREQ([2.59])

# Initialize with name, version, and support email address.
AC_INIT([GFDL atmos_phys Library],
  [main],
  [gfdl.climate.model.info@noaa.gov],
  [atmos_phys],
  [https://www.gfdl.noaa.gov/atmos_phys])

# Find out about the host we're building on.
AC_CANONICAL_HOST
AC_CANONICAL_BUILD

# Find out about the target we're building for.
AC_CANONICAL_TARGET

AM_INIT_AUTOMAKE([foreign dist-zip subdir-objects])
AC_REQUIRE_AUX_FILE([tap-driver.sh])

# Keep libtool macros in an m4 directory.
AC_CONFIG_MACRO_DIR([m4])

# Set up libtool.
LT_PREREQ([2.4])
LT_INIT()

AC_ARG_ENABLE([setting-flags],
  [AS_HELP_STRING([--enable-setting-flags],
    [Allow configure to set some compiler flags.  Disabling this will also disable any other --with or --enable options that set flags, and will only use user-provided falgs. (Default yes)])])
AS_IF([test ${enable_setting_flags:-yes} = yes],
  [enable_setting_flags=yes],
  [enable_setting_flags=no])

# Find the C compiler.
AC_PROG_CC
AM_PROG_CC_C_O
AC_C_CONST
AC_PROG_CPP

# Find the Fortran compiler.
AC_PROG_FC
AC_PROG_F77

# Find the install program.
AC_PROG_INSTALL

# Check to see if any macros must be set to enable large (>2GB) files.
AC_SYS_LARGEFILE

# Fortran specific checks
AC_LANG_PUSH(Fortran)
AC_OPENMP()
if test $enable_setting_flags = yes; then
  # Make sure the compiler is seeing this as free-formatted, not
  # fixed-formatted, fortran code.
  AC_FC_FREEFORM()

  # Check that long lines of Fortran code can be handled. This will add
  # necessary fortran flags.
  AC_FC_LINE_LENGTH([unlimited])
  GX_FC_DEFAULT_REAL_KIND8_FLAG([dnl
    FCFLAGS="$FCFLAGS $FC_DEFAULT_REAL_KIND8_FLAG"])
  FFLAGS="$FFLAGS $FC_DEFAULT_REAL_KIND8_FLAG"
  # Add OpenMP flags
  if test ! -z "$OPENMP_CFLAGS"; then
    CFLAGS="$CFLAGS $OPENMP_CFLAGS"
    CPPFLAGS="$CPPFLAGS $OPENMP_CFLAGS"
  fi
  if test ! -z "$OPENMP_FCFLAGS"; then
    FCFLAGS="$FCFLAGS $OPENMP_FCFLAGS"
    FFLAGS="$FFLAGS $OPENMP_FCFLAGS"
  fi
fi
GX_FC_CHECK_MOD([fms_mod], [], [], [AC_MSG_ERROR([Can't find the FMS module.  Set CPPFLAGS/FCFLAGS])])
GX_FORTRAN_SEARCH_LIBS([fms_init], [FMS], [use fms_mod], [call fms_init()], [],
  [AC_MSG_ERROR([Can't find the FMS library.  Set LDFLAGS/LIBS])])
AC_LANG_POP(Fortran)

# Define an AM_CONDITIONAL to determine if you are on a CRAY
AM_CONDITIONAL([CRAY], [test `env | grep CRAY | wc -l` -gt 0])

# Find which mpi launcher to use
AC_CHECK_PROGS([MPI_LAUNCHER],[srun aprun mpirun])

# Check if the launcher can oversubscribe the MPI processes
AS_IF([$MPI_LAUNCHER --oversubscribe hostname >/dev/null 2>&1], \
      [ AC_SUBST([OVERSUBSCRIBE], [--oversubscribe])])

# These files will be created when the configure script is run.
AC_CONFIG_FILES([Makefile])
AC_OUTPUT()
