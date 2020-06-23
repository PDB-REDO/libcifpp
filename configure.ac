AC_INIT([libcif++], 1.0, [m.hekkelman@nki.nl])

dnl Switch to a C++ compiler, and check if it works.
AC_LANG(C++)
AX_CXX_COMPILE_STDCXX_17([noext])

AC_CONFIG_SRCDIR([src/AtomShape.cpp])
AC_CONFIG_AUX_DIR(config)
AC_CONFIG_MACRO_DIR([config/m4])

AC_PREFIX_DEFAULT(/usr/local)

AC_DEFUN([read_test], [AC_LANG_SOURCE(
	esyscmd(tools/m4esc.sh config-tests/$1))])

AC_PROG_INSTALL
AC_ENABLE_STATIC
AC_DISABLE_SHARED
AC_PROG_LIBTOOL
AC_SUBST(LIBTOOL_DEPS)

AC_ARG_VAR([DEBUG], [Build a debug version of the library])

AC_PATH_PROG([PKG_CONFIG], [pkg-config])

AC_ARG_VAR([CCP4], [The location where CCP4 is installed])

AS_IF([test x"$CCP4" != x""],
	  [
		CPPFLAGS="$CPPFLAGS -I ${CCP4}/include"
		CXXFLAGS="$CXXFLAGS -I ${CCP4}/include"
		LDFLAGS="$LDFLAGS -L${CCP4}/lib"
	  ])

AX_BOOST_BASE([1.73])
AX_BOOST_REGEX

dnl The ccp4 distribution comes bundled with boost libraries that are
dnl most likely out of date. Since they also include .so files, this will
dnl break your software.

OLD_CXXFLAGS="$CXXFLAGS"
CXXFLAGS="$CXXFLAGS $BOOST_CPPFLAGS"
OLD_CPPFLAGS="$CPPFLAGS"
CPPFLAGS="$CPPFLAGS $BOOST_CPPFLAGS"
OLD_LDFLAGS="$LDFLAGS"
LDFLAGS="$LDFLAGS $BOOST_LDFLAGS"
OLD_LIBS="$LIBS"
LIBS="-lboost_regex"

AC_MSG_CHECKING([boost regex library version])
AC_LINK_IFELSE(
	[read_test(boost-regex-test.cpp)],
	[AC_MSG_RESULT([yes])],
	[AC_MSG_ERROR([You most likely have a conflict between the version
of the boost header files and the boost libraries located in ccp4.
This can be solved by deleting the libboost_*.so files in your CCP4
library directory.
Note: only the files with names ending in .so should be deleted.
])])

CXXFLAGS="$OLD_CXXFLAGS"
CPPFLAGS="$OLD_CPPFLAGS"
LDFLAGS="$OLD_LDFLAGS"
LIBS="$OLD_LIBS"

AX_BOOST_IOSTREAMS
dnl AX_BOOST_FILESYSTEM
AX_BOOST_THREAD

AC_ARG_WITH([clipper],
	AS_HELP_STRING([--with-clipper=@<:@location@:>@],
		[Use the clipper library as specified.
		@<:@location=$CCP4@:>@]),
		[
			CPPFLAGS="$CPPFLAGS -I ${withval}/include"
			CXXFLAGS="$CXXFLAGS -I ${withval}/include"
			LDFLAGS="$LDFLAGS -L${withval}/lib"
		])

AC_ARG_WITH([newuoa],
	AS_HELP_STRING([--with-newuoa=@<:@location@:>@],
		[Use the newuoa library as specified.]),
		[
			CPPFLAGS="$CPPFLAGS -I ${withval}/include"
			CXXFLAGS="$CXXFLAGS -I ${withval}/include"
			LDFLAGS="$LDFLAGS -L${withval}/lib"
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

AX_CHECK_LIBRARY([LIBZ], [zlib.h], [z], [],
                 [AC_MSG_ERROR([libz not found - compressed files not supported])])
AX_CHECK_LIBRARY([LIBBZ2], [bzlib.h], [bz2], [],
                 [AC_MSG_ERROR([libbz2 not found - compressed files not supported])])

# AC_ARG_ENABLE(
# 	debug,
# 	[AS_HELP_STRING([--enable-debug], [Generate a debug version of the library])],
# 	[],
# 	[enable_debug=auto])

dnl Process files
AC_OUTPUT([GNUmakefile
		   libcif++.pc])
