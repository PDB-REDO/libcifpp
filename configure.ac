AC_INIT([libcif++], 1.0, [m.hekkelman@nki.nl])

dnl Switch to a C++ compiler, and check if it works.
AC_LANG(C++)
AX_CXX_COMPILE_STDCXX_14([noext])

AC_CONFIG_SRCDIR([src/AtomShape.cpp])
AC_CONFIG_AUX_DIR(config)
AC_CONFIG_MACRO_DIR([config/m4])

AC_PREFIX_DEFAULT(/usr/local)

AC_DEFUN([read_test], [AC_LANG_SOURCE(
	esyscmd(tools/m4esc.sh config-tests/$1))])

AX_BOOST_BASE([1.73])

AC_PROG_INSTALL
dnl AC_DISABLE_SHARED
dnl AC_PROG_LIBTOOL
dnl AC_SUBST(LIBTOOL_DEPS)

AC_ARG_VAR([DEBUG], [Build a debug version of the library])

AC_PATH_PROG([PKG_CONFIG], [pkg-config])

AC_ARG_VAR([CCP4], [The location where CCP4 is installed])

AS_IF([test x"$CCP4" != x""],
	  [
		CPPFLAGS="$CPPFLAGS -I ${CCP4}/include"
		CXXFLAGS="$CXXFLAGS -I ${CCP4}/include"
	  ])

AC_ARG_WITH([clipper],
	AS_HELP_STRING([--with-clipper=@<:@location@:>@],
		[Use the clipper library as specified.
		@<:@location=$CCP4@:>@]),
		[
			CPPFLAGS="$CPPFLAGS -I ${withval}/include"
			CXXFLAGS="$CXXFLAGS -I ${withval}/include"
		])

AC_ARG_WITH([newuoa],
	AS_HELP_STRING([--with-newuoa=@<:@location@:>@],
		[Use the newuoa library as specified.]),
		[
			CPPFLAGS="$CPPFLAGS -I ${withval}/include"
			CXXFLAGS="$CXXFLAGS -I ${withval}/include"
		])

AC_CHECK_HEADER(
	[clipper/clipper.h],
	[],
	[AC_MSG_ERROR([
Can't find the main clipper include file clipper/clipper.h. Please
specify either the CCP4 location of the location of the installed
clipper using --with-clipper=<path>
])])

AC_MSG_CHECKING([clipper version])
AC_COMPILE_IFELSE(
	[read_test(clipper-test.cpp)],
	[],
	[AC_MSG_ERROR([The version of clipper is not up to date])])

AC_CHECK_HEADER(
	[newuoa.h],
	[],
	[AC_MSG_ERROR([
Can't find the newuoa include file newuoa.h. Please install this
library and specify its location with --with-newuoa=<path>.
You can find newuoa at https://github.com/elsid/newuoa-cpp.git
])])

# AC_ARG_ENABLE(
# 	debug,
# 	[AS_HELP_STRING([--enable-debug], [Generate a debug version of the library])],
# 	[],
# 	[enable_debug=auto])

dnl Process Makefile.in to create Makefile
AC_OUTPUT([GNUmakefile])
