AC_PREREQ([2.69])

AC_INIT([libcifpp], 1.0, [m.hekkelman@nki.nl])

dnl Switch to a C++ compiler, and check if it works.
AC_LANG(C++)
AX_CXX_COMPILE_STDCXX_17([noext])

AX_CHECK_COMPILE_FLAG([-fstandalone-debug], [ CXXFLAGS="$CXXFLAGS -fstandalone-debug" ], , [-Werror])

AC_CONFIG_SRCDIR([src/Cif++.cpp])
AC_CONFIG_AUX_DIR(config)
AC_CONFIG_MACRO_DIR([config/m4])
AC_CONFIG_HEADERS([include/cif++/Config.hpp])

AC_PREFIX_DEFAULT(/usr/local)

AC_CHECK_HEADER([filesystem], [], [AC_MSG_ERROR([The file <filesystem> is missing, perhaps you should install a more recent libstdc++ implementation.])])

dnl check if we need stdc++fs as library
AC_MSG_CHECKING([if linking to stdc++fs is required])
AC_TRY_LINK(
	[#include <filesystem>],
	[(void)std::filesystem::current_path();],
	[
		AC_MSG_RESULT(no)
	],
	[
		LIBS="$LIBS -lstdc++fs"

		AC_TRY_LINK(
			[#include <filesystem>],
			[(void)std::filesystem::current_path();],
			[
				AC_MSG_RESULT(yes)
			],
			[
				AC_MSG_ERROR([Could not link filesystem])
			]
		)
	]
)

dnl AC_DEFUN([read_test], [AC_LANG_SOURCE(
dnl 	esyscmd(tools/m4esc.sh config-tests/$1))])
dnl 
dnl AC_MSG_CHECKING([compiler standards compliance])
dnl AC_COMPILE_IFELSE(
dnl 	[read_test(cpp-17-test.cpp)], [],
dnl 	[AC_MSG_ERROR([Your c++ compiler is not capable of compiling libcifpp, please upgrade])])
dnl AC_MSG_RESULT(ok)

AC_PROG_INSTALL

dnl Shared libraries are not convenient
dnl LT_INIT([disable-shared])
LT_INIT
AC_SUBST(LIBTOOL_DEPS)

dnl versioning, first for libtool
LIBCIF_CURRENT=1
LIBCIF_REVISION=0
LIBCIF_AGE=0

LIBCIF_LT_CURRENT="${LIBCIF_CURRENT}"
LIBCIF_LT_VERSION="${LIBCIF_CURRENT}:${LIBCIF_REVISION}:${LIBCIF_AGE}"

AC_SUBST(LIBCIF_LT_CURRENT)
AC_SUBST(LIBCIF_LT_VERSION)

dnl and now for the semantic version
LIBCIF_SEMANTIC_VERSION=1.0.0
AC_SUBST(LIBCIF_SEMANTIC_VERSION)

AC_ARG_VAR([DEBUG], [Build a debug version of the library])

dnl revision numbering is something used internally at the NKI
AC_ARG_ENABLE(
	revision,
	[AS_HELP_STRING([--enable-revision], [Create a build number as revision])])

AS_IF([test "x$enable_revision" != "xyes" ], [
	UPDATE_REVISION=1
])

dnl disable resources on macOS
AS_CASE([$host],
	[*-apple-*], [USE_RSRC=0],
	[USE_RSRC=1]
)
AC_SUBST([USE_RSRC])

AC_SUBST([UPDATE_REVISION], [$UPDATE_REVISION])

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
AX_BOOST_IOSTREAMS
AX_BOOST_DATE_TIME

AC_OUTPUT([GNUmakefile
		   libcifpp.pc])
