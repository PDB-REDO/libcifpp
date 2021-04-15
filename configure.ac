AC_PREREQ([2.69])

m4_define([libcifpp_version_major],[1])
m4_define([libcifpp_version_minor],[1])
m4_define([libcifpp_version_micro],[0])
m4_define([libcifpp_version_extra],[])
m4_define([libcifpp_version],[libcifpp_version_major().libcifpp_version_minor().libcifpp_version_micro()libcifpp_version_extra])

AC_INIT([libcifpp], [libcifpp_version], [m.hekkelman@nki.nl])

dnl Switch to a C++ compiler, and check if it works.
AC_LANG(C++)
AX_CXX_COMPILE_STDCXX_17([noext])

AX_CHECK_COMPILE_FLAG([-fstandalone-debug],
	[
		CXXFLAGS="$CXXFLAGS -fstandalone-debug"
	] , , [-Werror])

AC_CONFIG_SRCDIR([src/Cif++.cpp])
AC_CONFIG_AUX_DIR(config)
AC_CONFIG_MACRO_DIR([config/m4])
AC_CONFIG_HEADERS([src/Config.hpp])

AC_PREFIX_DEFAULT(/usr/local)

AX_FILESYSTEM

AC_PROG_INSTALL
PKG_PROG_PKG_CONFIG

AX_PTHREAD

AC_ARG_VAR([DEBUG], [Build a debug version of the application])

dnl disable shared libraries, they break the build on MacOS... sigh
LT_INIT([disable-shared pic-only])

AC_SUBST(LIBTOOL_DEPS)

dnl versioning, first for libtool
LIBCIF_CURRENT=1
LIBCIF_REVISION=1
LIBCIF_AGE=1

LIBCIF_LT_CURRENT="${LIBCIF_CURRENT}"
LIBCIF_LT_VERSION="${LIBCIF_CURRENT}:${LIBCIF_REVISION}:${LIBCIF_AGE}"

AC_SUBST(LIBCIF_LT_CURRENT)
AC_SUBST(LIBCIF_LT_VERSION)

dnl and now for the semantic version
LIBCIF_SEMANTIC_VERSION=libcifpp_version
AC_SUBST(LIBCIF_SEMANTIC_VERSION)

AC_ARG_VAR([CCP4], [The location where CCP4 is installed. CCP4 is only required if you want to rebuild the symmetry operations table.])

AX_MRC

dnl revision numbering is something used internally at the NKI
AC_ARG_ENABLE(
	revision,
	[AS_HELP_STRING([--disable-revision], [Create a build number as revision])])

AS_IF([test "x$enable_revision" != "xno" ], [
	UPDATE_REVISION=1
])

AC_SUBST([UPDATE_REVISION], [$UPDATE_REVISION])

AX_FILESYSTEM

AC_CHECK_FUNCS([floor pow rint sqrt strchr strerror])
AC_CHECK_HEADERS([sys/ioctl.h])
AC_CHECK_HEADERS([termios.h])
AC_CHECK_HEADER_STDBOOL
AC_CHECK_TYPES([ptrdiff_t])
AC_PROG_MAKE_SET

AX_BOOST_BASE([1.71], [], [AC_MSG_ERROR([Sorry, your boost is not found or not up-to-date.])])
AX_BOOST_IOSTREAMS
AX_BOOST_DATE_TIME
AX_BOOST_REGEX

AX_IOSTREAMS_Z
AX_IOSTREAMS_BZ2

dnl AC_CHECK_LIB([atomic], [atomic_flag_clear])

dnl Set output variables for the various directories
AC_ARG_VAR([DATA_LIB_DIR], [The location where to store packaged dictionary and CCD files])
AC_SUBST([DATA_LIB_DIR])
DATA_LIB_DIR=$datadir/libcifpp

AC_ARG_VAR([DATA_CACHE_DIR], [The location where to store cached dictionary and CCD files])
AC_SUBST([DATA_CACHE_DIR])
DATA_CACHE_DIR=/var/cache/libcifpp

AC_OUTPUT([GNUmakefile
		   libcifpp.pc
		   tools/update-dictionary-script])
