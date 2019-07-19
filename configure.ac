## Process this file with autoconf to produce a configure script.

# The version set here will propagate to other files from here
AC_INIT([rLVIS], 0.0.1)

# Checks for common programs using default macros
AC_PROG_CC

## Use gsl-config to find arguments for compiler and linker flags
##
## Check for non-standard programs: gsl-config(1)
AC_PATH_PROG([GSL_CONFIG], [gsl-config])
## If gsl-config was found, let's use it
if test "${GSL_CONFIG}" != ""; then
# Use gsl-config for header and linker arguments
GSL_CFLAGS=`${GSL_CONFIG} --cflags`
GSL_LIBS=`${GSL_CONFIG} --libs`
else
  AC_MSG_ERROR([gsl-config not found, is GSL installed?])
fi

AC_CHECK_LIB(hdf5, H5Dcreate2, [],
      AC_MSG_ERROR([Static hdf5 library not found.
      Please install the hdf5 library with C++ interface enabled (--enable-cxx=yes).
  	The required HDF5 library files can be installed as follows:
  	    - Debian-based (e.g. Debian >= 8.0, Ubuntu >= 15.04): 'sudo apt-get install libhdf5-dev'
        - Old Debian-based (e.g Debian < 8.0, Ubuntu < 15.04): Install from source (see INSTALL)
        - OS X using Homebrew: 'brew install homebrew/science/hdf5 --enable-cxx'
        - RPM-based (e.g Fedora): 'sudo yum install hdf5-devel']))


# Now substitute these variables in src/Makevars.in to create src/Makevars
AC_SUBST(GSL_CFLAGS)
AC_SUBST(GSL_LIBS)

AC_OUTPUT(src/Makevars)
