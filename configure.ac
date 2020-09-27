AC_PREREQ([2.69])

AC_INIT([libcif++], 1.0, [m.hekkelman@nki.nl])

dnl Switch to a C++ compiler, and check if it works.
AC_LANG(C++)
AX_CXX_COMPILE_STDCXX_17([noext])

AX_CHECK_COMPILE_FLAG([-fstandalone-debug], [ CXXFLAGS="$CXXFLAGS -fstandalone-debug" ], , [-Werror])

AC_CONFIG_SRCDIR([src/Cif++.cpp])
AC_CONFIG_AUX_DIR(config)
AC_CONFIG_MACRO_DIR([config/m4])
AC_CONFIG_HEADERS([include/cif++/Config.hpp])

AC_PREFIX_DEFAULT(/usr/local)

dnl AC_DEFUN([read_test], [AC_LANG_SOURCE(
dnl 	esyscmd(tools/m4esc.sh config-tests/$1))])
dnl 
dnl AC_MSG_CHECKING([compiler standards compliance])
dnl AC_COMPILE_IFELSE(
dnl 	[read_test(cpp-17-test.cpp)], [],
dnl 	[AC_MSG_ERROR([Your c++ compiler is not capable of compiling libcifpp, please upgrade])])
dnl AC_MSG_RESULT(ok)

AC_PROG_INSTALL
AC_ENABLE_STATIC
AC_DISABLE_SHARED
AC_PROG_LIBTOOL
AC_SUBST(LIBTOOL_DEPS)

AC_ARG_VAR([DEBUG], [Build a debug version of the library])

AC_ARG_VAR([CCP4], [The location where CCP4 is installed])

AC_PATH_PROG([PKG_CONFIG], [pkg-config])

AC_CHECK_FUNCS([floor pow rint sqrt strchr strerror])
AC_CHECK_HEADERS([sys/ioctl.h])
AC_CHECK_HEADERS([termios.h])
AC_CHECK_HEADER_STDBOOL
AC_CHECK_TYPES([ptrdiff_t])
AC_PROG_MAKE_SET

AX_BOOST_BASE([1.71], [], [
	AC_MSG_ERROR([Sorry, your boost is not found or not up-to-date.])
])
AX_BOOST_REGEX

AX_CHECK_LIBRARY([LIBZ], [zlib.h], [z], [],
                 [AC_MSG_ERROR([libz not found - compressed files not supported])])
AX_CHECK_LIBRARY([LIBBZ2], [bzlib.h], [bz2], [],
                 [AC_MSG_ERROR([libbz2 not found - compressed files not supported])])

AC_OUTPUT([GNUmakefile
		   libcif++.pc])
