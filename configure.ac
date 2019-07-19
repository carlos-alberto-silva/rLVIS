## Process this file with autoconf to produce a configure script.
LDFLAGS="$LDFLAGS -L /usr/lib/x86_64-linux-gnu/hdf5/serial

# The version set here will propagate to other files from here
AC_INIT([rLVIS], 0.0.1)

# Checks for common programs using default macros
AC_PROG_CC

##################
## GREP
##################
AC_CHECK_PROG(GREP,grep,`which grep`)
AC_MSG_CHECKING([whether ${GREP} accepts -o])
if `echo '' | ${GREP} -o '' >/dev/null 2>&1`
then
  AC_MSG_RESULT([yes])
else
  AC_MSG_RESULT([no])
  AC_MSG_ERROR([${GREP} does not support option -o])
fi


##################
## GSL
##################

AC_PATH_PROG([GSL_CONFIG], [gsl-config])

## If gsl-config was found, let's use it
if test "${GSL_CONFIG}" != ""; then

# Use gsl-config for header and linker arguments
GSL_CFLAGS=`${GSL_CONFIG} --cflags`
GSL_LIBS=`${GSL_CONFIG} --libs`
else
  AC_MSG_ERROR([gsl-config not found, is libgsl installed?
    Please install the libgsl.
  	The required gsl library files can be installed as follows:
  	    - Debian-based: 'sudo apt-get install libgsl-dev'
        - OS X using Homebrew: 'brew install gsl'
        - RPM-based (e.g Fedora): 'sudo yum install gsl-devel'])
fi


##################
## HDF5
##################
AC_CHECK_PROG(H5CC,h5cc,`which h5cc`)
if test -z "${H5CC}"; then
  AC_MSG_ERROR([Static hdf5 library not found.
    Please install the hdf5 library.
  	The required HDF5 library files can be installed as follows:
  	    - Debian-based (e.g. Debian >= 8.0, Ubuntu >= 15.04): 'sudo apt-get install libhdf5-dev'
        - Old Debian-based (e.g Debian < 8.0, Ubuntu < 15.04): Install from source (see INSTALL)
        - OS X using Homebrew: 'brew install homebrew/science/hdf5'
        - RPM-based (e.g Fedora): 'sudo yum install hdf5-devel'])
fi
HDF5_LIBS=`${H5CC} -show | ${GREP} -o "\-L[^[[:space:]]]\+\|\-l[a-zA-Z0-9]\+ " | tr "\n" " " | sed -e 's/  / /g'`
HDF5_LIBS+="-lhdf5 -lhdf5_hl"

echo 'int main() { return; }' > myconftest.c
HDF5_CFLAGS=`${H5CC} -show -c myconftest.c | ${GREP} -o "\-I[a-zA-Z0-9\_\-\/]\+" | tr "\n" " "`
rm -f myconftest.c



# Now substitute these variables in src/Makevars.in to create src/Makevars
AC_SUBST(GSL_CFLAGS)
AC_SUBST(GSL_LIBS)
AC_SUBST(HDF5_CFLAGS)
AC_SUBST(HDF5_LIBS)

AC_OUTPUT(src/Makevars)