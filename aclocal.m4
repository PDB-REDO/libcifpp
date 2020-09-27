# generated automatically by aclocal 1.16.1 -*- Autoconf -*-

# Copyright (C) 1996-2018 Free Software Foundation, Inc.

# This file is free software; the Free Software Foundation
# gives unlimited permission to copy and/or distribute it,
# with or without modifications, as long as this notice is preserved.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY, to the extent permitted by law; without
# even the implied warranty of MERCHANTABILITY or FITNESS FOR A
# PARTICULAR PURPOSE.

m4_ifndef([AC_CONFIG_MACRO_DIRS], [m4_defun([_AM_CONFIG_MACRO_DIRS], [])m4_defun([AC_CONFIG_MACRO_DIRS], [_AM_CONFIG_MACRO_DIRS($@)])])
# ===========================================================================
#      https://www.gnu.org/software/autoconf-archive/ax_boost_base.html
# ===========================================================================
#
# SYNOPSIS
#
#   AX_BOOST_BASE([MINIMUM-VERSION], [ACTION-IF-FOUND], [ACTION-IF-NOT-FOUND])
#
# DESCRIPTION
#
#   Test for the Boost C++ libraries of a particular version (or newer)
#
#   If no path to the installed boost library is given the macro searchs
#   under /usr, /usr/local, /opt and /opt/local and evaluates the
#   $BOOST_ROOT environment variable. Further documentation is available at
#   <http://randspringer.de/boost/index.html>.
#
#   This macro calls:
#
#     AC_SUBST(BOOST_CPPFLAGS) / AC_SUBST(BOOST_LDFLAGS)
#
#   And sets:
#
#     HAVE_BOOST
#
# LICENSE
#
#   Copyright (c) 2008 Thomas Porschberg <thomas@randspringer.de>
#   Copyright (c) 2009 Peter Adolphs
#
#   Copying and distribution of this file, with or without modification, are
#   permitted in any medium without royalty provided the copyright notice
#   and this notice are preserved. This file is offered as-is, without any
#   warranty.

#serial 45

# example boost program (need to pass version)
m4_define([_AX_BOOST_BASE_PROGRAM],
          [AC_LANG_PROGRAM([[
#include <boost/version.hpp>
]],[[
(void) ((void)sizeof(char[1 - 2*!!((BOOST_VERSION) < ($1))]));
]])])

AC_DEFUN([AX_BOOST_BASE],
[
AC_ARG_WITH([boost],
  [AS_HELP_STRING([--with-boost@<:@=ARG@:>@],
    [use Boost library from a standard location (ARG=yes),
     from the specified location (ARG=<path>),
     or disable it (ARG=no)
     @<:@ARG=yes@:>@ ])],
    [
     AS_CASE([$withval],
       [no],[want_boost="no";_AX_BOOST_BASE_boost_path=""],
       [yes],[want_boost="yes";_AX_BOOST_BASE_boost_path=""],
       [want_boost="yes";_AX_BOOST_BASE_boost_path="$withval"])
    ],
    [want_boost="yes"])


AC_ARG_WITH([boost-libdir],
  [AS_HELP_STRING([--with-boost-libdir=LIB_DIR],
    [Force given directory for boost libraries.
     Note that this will override library path detection,
     so use this parameter only if default library detection fails
     and you know exactly where your boost libraries are located.])],
  [
   AS_IF([test -d "$withval"],
         [_AX_BOOST_BASE_boost_lib_path="$withval"],
    [AC_MSG_ERROR([--with-boost-libdir expected directory name])])
  ],
  [_AX_BOOST_BASE_boost_lib_path=""])

BOOST_LDFLAGS=""
BOOST_CPPFLAGS=""
AS_IF([test "x$want_boost" = "xyes"],
      [_AX_BOOST_BASE_RUNDETECT([$1],[$2],[$3])])
AC_SUBST(BOOST_CPPFLAGS)
AC_SUBST(BOOST_LDFLAGS)
])


# convert a version string in $2 to numeric and affect to polymorphic var $1
AC_DEFUN([_AX_BOOST_BASE_TONUMERICVERSION],[
  AS_IF([test "x$2" = "x"],[_AX_BOOST_BASE_TONUMERICVERSION_req="1.20.0"],[_AX_BOOST_BASE_TONUMERICVERSION_req="$2"])
  _AX_BOOST_BASE_TONUMERICVERSION_req_shorten=`expr $_AX_BOOST_BASE_TONUMERICVERSION_req : '\([[0-9]]*\.[[0-9]]*\)'`
  _AX_BOOST_BASE_TONUMERICVERSION_req_major=`expr $_AX_BOOST_BASE_TONUMERICVERSION_req : '\([[0-9]]*\)'`
  AS_IF([test "x$_AX_BOOST_BASE_TONUMERICVERSION_req_major" = "x"],
        [AC_MSG_ERROR([You should at least specify libboost major version])])
  _AX_BOOST_BASE_TONUMERICVERSION_req_minor=`expr $_AX_BOOST_BASE_TONUMERICVERSION_req : '[[0-9]]*\.\([[0-9]]*\)'`
  AS_IF([test "x$_AX_BOOST_BASE_TONUMERICVERSION_req_minor" = "x"],
        [_AX_BOOST_BASE_TONUMERICVERSION_req_minor="0"])
  _AX_BOOST_BASE_TONUMERICVERSION_req_sub_minor=`expr $_AX_BOOST_BASE_TONUMERICVERSION_req : '[[0-9]]*\.[[0-9]]*\.\([[0-9]]*\)'`
  AS_IF([test "X$_AX_BOOST_BASE_TONUMERICVERSION_req_sub_minor" = "X"],
        [_AX_BOOST_BASE_TONUMERICVERSION_req_sub_minor="0"])
  _AX_BOOST_BASE_TONUMERICVERSION_RET=`expr $_AX_BOOST_BASE_TONUMERICVERSION_req_major \* 100000 \+  $_AX_BOOST_BASE_TONUMERICVERSION_req_minor \* 100 \+ $_AX_BOOST_BASE_TONUMERICVERSION_req_sub_minor`
  AS_VAR_SET($1,$_AX_BOOST_BASE_TONUMERICVERSION_RET)
])

dnl Run the detection of boost should be run only if $want_boost
AC_DEFUN([_AX_BOOST_BASE_RUNDETECT],[
    _AX_BOOST_BASE_TONUMERICVERSION(WANT_BOOST_VERSION,[$1])
    succeeded=no


    AC_REQUIRE([AC_CANONICAL_HOST])
    dnl On 64-bit systems check for system libraries in both lib64 and lib.
    dnl The former is specified by FHS, but e.g. Debian does not adhere to
    dnl this (as it rises problems for generic multi-arch support).
    dnl The last entry in the list is chosen by default when no libraries
    dnl are found, e.g. when only header-only libraries are installed!
    AS_CASE([${host_cpu}],
      [x86_64],[libsubdirs="lib64 libx32 lib lib64"],
      [ppc64|powerpc64|s390x|sparc64|aarch64|ppc64le|powerpc64le|riscv64],[libsubdirs="lib64 lib lib64"],
      [libsubdirs="lib"]
    )

    dnl allow for real multi-arch paths e.g. /usr/lib/x86_64-linux-gnu. Give
    dnl them priority over the other paths since, if libs are found there, they
    dnl are almost assuredly the ones desired.
    AS_CASE([${host_cpu}],
      [i?86],[multiarch_libsubdir="lib/i386-${host_os}"],
      [multiarch_libsubdir="lib/${host_cpu}-${host_os}"]
    )

    dnl first we check the system location for boost libraries
    dnl this location ist chosen if boost libraries are installed with the --layout=system option
    dnl or if you install boost with RPM
    AS_IF([test "x$_AX_BOOST_BASE_boost_path" != "x"],[
        AC_MSG_CHECKING([for boostlib >= $1 ($WANT_BOOST_VERSION) includes in "$_AX_BOOST_BASE_boost_path/include"])
         AS_IF([test -d "$_AX_BOOST_BASE_boost_path/include" && test -r "$_AX_BOOST_BASE_boost_path/include"],[
           AC_MSG_RESULT([yes])
           BOOST_CPPFLAGS="-I$_AX_BOOST_BASE_boost_path/include"
           for _AX_BOOST_BASE_boost_path_tmp in $multiarch_libsubdir $libsubdirs; do
                AC_MSG_CHECKING([for boostlib >= $1 ($WANT_BOOST_VERSION) lib path in "$_AX_BOOST_BASE_boost_path/$_AX_BOOST_BASE_boost_path_tmp"])
                AS_IF([test -d "$_AX_BOOST_BASE_boost_path/$_AX_BOOST_BASE_boost_path_tmp" && test -r "$_AX_BOOST_BASE_boost_path/$_AX_BOOST_BASE_boost_path_tmp" ],[
                        AC_MSG_RESULT([yes])
                        BOOST_LDFLAGS="-L$_AX_BOOST_BASE_boost_path/$_AX_BOOST_BASE_boost_path_tmp";
                        break;
                ],
      [AC_MSG_RESULT([no])])
           done],[
      AC_MSG_RESULT([no])])
    ],[
        if test X"$cross_compiling" = Xyes; then
            search_libsubdirs=$multiarch_libsubdir
        else
            search_libsubdirs="$multiarch_libsubdir $libsubdirs"
        fi
        for _AX_BOOST_BASE_boost_path_tmp in /usr /usr/local /opt /opt/local ; do
            if test -d "$_AX_BOOST_BASE_boost_path_tmp/include/boost" && test -r "$_AX_BOOST_BASE_boost_path_tmp/include/boost" ; then
                for libsubdir in $search_libsubdirs ; do
                    if ls "$_AX_BOOST_BASE_boost_path_tmp/$libsubdir/libboost_"* >/dev/null 2>&1 ; then break; fi
                done
                BOOST_LDFLAGS="-L$_AX_BOOST_BASE_boost_path_tmp/$libsubdir"
                BOOST_CPPFLAGS="-I$_AX_BOOST_BASE_boost_path_tmp/include"
                break;
            fi
        done
    ])

    dnl overwrite ld flags if we have required special directory with
    dnl --with-boost-libdir parameter
    AS_IF([test "x$_AX_BOOST_BASE_boost_lib_path" != "x"],
          [BOOST_LDFLAGS="-L$_AX_BOOST_BASE_boost_lib_path"])

    AC_MSG_CHECKING([for boostlib >= $1 ($WANT_BOOST_VERSION)])
    CPPFLAGS_SAVED="$CPPFLAGS"
    CPPFLAGS="$CPPFLAGS $BOOST_CPPFLAGS"
    export CPPFLAGS

    LDFLAGS_SAVED="$LDFLAGS"
    LDFLAGS="$LDFLAGS $BOOST_LDFLAGS"
    export LDFLAGS

    AC_REQUIRE([AC_PROG_CXX])
    AC_LANG_PUSH(C++)
        AC_COMPILE_IFELSE([_AX_BOOST_BASE_PROGRAM($WANT_BOOST_VERSION)],[
        AC_MSG_RESULT(yes)
    succeeded=yes
    found_system=yes
        ],[
        ])
    AC_LANG_POP([C++])



    dnl if we found no boost with system layout we search for boost libraries
    dnl built and installed without the --layout=system option or for a staged(not installed) version
    if test "x$succeeded" != "xyes" ; then
        CPPFLAGS="$CPPFLAGS_SAVED"
        LDFLAGS="$LDFLAGS_SAVED"
        BOOST_CPPFLAGS=
        if test -z "$_AX_BOOST_BASE_boost_lib_path" ; then
            BOOST_LDFLAGS=
        fi
        _version=0
        if test -n "$_AX_BOOST_BASE_boost_path" ; then
            if test -d "$_AX_BOOST_BASE_boost_path" && test -r "$_AX_BOOST_BASE_boost_path"; then
                for i in `ls -d $_AX_BOOST_BASE_boost_path/include/boost-* 2>/dev/null`; do
                    _version_tmp=`echo $i | sed "s#$_AX_BOOST_BASE_boost_path##" | sed 's/\/include\/boost-//' | sed 's/_/./'`
                    V_CHECK=`expr $_version_tmp \> $_version`
                    if test "x$V_CHECK" = "x1" ; then
                        _version=$_version_tmp
                    fi
                    VERSION_UNDERSCORE=`echo $_version | sed 's/\./_/'`
                    BOOST_CPPFLAGS="-I$_AX_BOOST_BASE_boost_path/include/boost-$VERSION_UNDERSCORE"
                done
                dnl if nothing found search for layout used in Windows distributions
                if test -z "$BOOST_CPPFLAGS"; then
                    if test -d "$_AX_BOOST_BASE_boost_path/boost" && test -r "$_AX_BOOST_BASE_boost_path/boost"; then
                        BOOST_CPPFLAGS="-I$_AX_BOOST_BASE_boost_path"
                    fi
                fi
                dnl if we found something and BOOST_LDFLAGS was unset before
                dnl (because "$_AX_BOOST_BASE_boost_lib_path" = ""), set it here.
                if test -n "$BOOST_CPPFLAGS" && test -z "$BOOST_LDFLAGS"; then
                    for libsubdir in $libsubdirs ; do
                        if ls "$_AX_BOOST_BASE_boost_path/$libsubdir/libboost_"* >/dev/null 2>&1 ; then break; fi
                    done
                    BOOST_LDFLAGS="-L$_AX_BOOST_BASE_boost_path/$libsubdir"
                fi
            fi
        else
            if test "x$cross_compiling" != "xyes" ; then
                for _AX_BOOST_BASE_boost_path in /usr /usr/local /opt /opt/local ; do
                    if test -d "$_AX_BOOST_BASE_boost_path" && test -r "$_AX_BOOST_BASE_boost_path" ; then
                        for i in `ls -d $_AX_BOOST_BASE_boost_path/include/boost-* 2>/dev/null`; do
                            _version_tmp=`echo $i | sed "s#$_AX_BOOST_BASE_boost_path##" | sed 's/\/include\/boost-//' | sed 's/_/./'`
                            V_CHECK=`expr $_version_tmp \> $_version`
                            if test "x$V_CHECK" = "x1" ; then
                                _version=$_version_tmp
                                best_path=$_AX_BOOST_BASE_boost_path
                            fi
                        done
                    fi
                done

                VERSION_UNDERSCORE=`echo $_version | sed 's/\./_/'`
                BOOST_CPPFLAGS="-I$best_path/include/boost-$VERSION_UNDERSCORE"
                if test -z "$_AX_BOOST_BASE_boost_lib_path" ; then
                    for libsubdir in $libsubdirs ; do
                        if ls "$best_path/$libsubdir/libboost_"* >/dev/null 2>&1 ; then break; fi
                    done
                    BOOST_LDFLAGS="-L$best_path/$libsubdir"
                fi
            fi

            if test -n "$BOOST_ROOT" ; then
                for libsubdir in $libsubdirs ; do
                    if ls "$BOOST_ROOT/stage/$libsubdir/libboost_"* >/dev/null 2>&1 ; then break; fi
                done
                if test -d "$BOOST_ROOT" && test -r "$BOOST_ROOT" && test -d "$BOOST_ROOT/stage/$libsubdir" && test -r "$BOOST_ROOT/stage/$libsubdir"; then
                    version_dir=`expr //$BOOST_ROOT : '.*/\(.*\)'`
                    stage_version=`echo $version_dir | sed 's/boost_//' | sed 's/_/./g'`
                        stage_version_shorten=`expr $stage_version : '\([[0-9]]*\.[[0-9]]*\)'`
                    V_CHECK=`expr $stage_version_shorten \>\= $_version`
                    if test "x$V_CHECK" = "x1" && test -z "$_AX_BOOST_BASE_boost_lib_path" ; then
                        AC_MSG_NOTICE(We will use a staged boost library from $BOOST_ROOT)
                        BOOST_CPPFLAGS="-I$BOOST_ROOT"
                        BOOST_LDFLAGS="-L$BOOST_ROOT/stage/$libsubdir"
                    fi
                fi
            fi
        fi

        CPPFLAGS="$CPPFLAGS $BOOST_CPPFLAGS"
        export CPPFLAGS
        LDFLAGS="$LDFLAGS $BOOST_LDFLAGS"
        export LDFLAGS

        AC_LANG_PUSH(C++)
            AC_COMPILE_IFELSE([_AX_BOOST_BASE_PROGRAM($WANT_BOOST_VERSION)],[
            AC_MSG_RESULT(yes)
        succeeded=yes
        found_system=yes
            ],[
            ])
        AC_LANG_POP([C++])
    fi

    if test "x$succeeded" != "xyes" ; then
        if test "x$_version" = "x0" ; then
            AC_MSG_NOTICE([[We could not detect the boost libraries (version $1 or higher). If you have a staged boost library (still not installed) please specify \$BOOST_ROOT in your environment and do not give a PATH to --with-boost option.  If you are sure you have boost installed, then check your version number looking in <boost/version.hpp>. See http://randspringer.de/boost for more documentation.]])
        else
            AC_MSG_NOTICE([Your boost libraries seems to old (version $_version).])
        fi
        # execute ACTION-IF-NOT-FOUND (if present):
        ifelse([$3], , :, [$3])
    else
        AC_DEFINE(HAVE_BOOST,,[define if the Boost library is available])
        # execute ACTION-IF-FOUND (if present):
        ifelse([$2], , :, [$2])
    fi

    CPPFLAGS="$CPPFLAGS_SAVED"
    LDFLAGS="$LDFLAGS_SAVED"

])

# ===========================================================================
#    https://www.gnu.org/software/autoconf-archive/ax_boost_iostreams.html
# ===========================================================================
#
# SYNOPSIS
#
#   AX_BOOST_IOSTREAMS
#
# DESCRIPTION
#
#   Test for IOStreams library from the Boost C++ libraries. The macro
#   requires a preceding call to AX_BOOST_BASE. Further documentation is
#   available at <http://randspringer.de/boost/index.html>.
#
#   This macro calls:
#
#     AC_SUBST(BOOST_IOSTREAMS_LIB)
#
#   And sets:
#
#     HAVE_BOOST_IOSTREAMS
#
# LICENSE
#
#   Copyright (c) 2008 Thomas Porschberg <thomas@randspringer.de>
#
#   Copying and distribution of this file, with or without modification, are
#   permitted in any medium without royalty provided the copyright notice
#   and this notice are preserved. This file is offered as-is, without any
#   warranty.

#serial 21

AC_DEFUN([AX_BOOST_IOSTREAMS],
[
	AC_ARG_WITH([boost-iostreams],
	AS_HELP_STRING([--with-boost-iostreams@<:@=special-lib@:>@],
                   [use the IOStreams library from boost - it is possible to specify a certain library for the linker
                        e.g. --with-boost-iostreams=boost_iostreams-gcc-mt-d-1_33_1 ]),
        [
        if test "$withval" = "no"; then
			want_boost="no"
        elif test "$withval" = "yes"; then
            want_boost="yes"
            ax_boost_user_iostreams_lib=""
        else
		    want_boost="yes"
		ax_boost_user_iostreams_lib="$withval"
		fi
        ],
        [want_boost="yes"]
	)

	if test "x$want_boost" = "xyes"; then
        AC_REQUIRE([AC_PROG_CC])
		CPPFLAGS_SAVED="$CPPFLAGS"
		CPPFLAGS="$CPPFLAGS $BOOST_CPPFLAGS"
		export CPPFLAGS

		LDFLAGS_SAVED="$LDFLAGS"
		LDFLAGS="$LDFLAGS $BOOST_LDFLAGS"
		export LDFLAGS

        AC_CACHE_CHECK(whether the Boost::IOStreams library is available,
					   ax_cv_boost_iostreams,
        [AC_LANG_PUSH([C++])
		 AC_COMPILE_IFELSE([AC_LANG_PROGRAM([[@%:@include <boost/iostreams/filtering_stream.hpp>
											 @%:@include <boost/range/iterator_range.hpp>
											]],
                                  [[std::string  input = "Hello World!";
								 namespace io = boost::iostreams;
									 io::filtering_istream  in(boost::make_iterator_range(input));
									 return 0;
                                   ]])],
                             ax_cv_boost_iostreams=yes, ax_cv_boost_iostreams=no)
         AC_LANG_POP([C++])
		])
		if test "x$ax_cv_boost_iostreams" = "xyes"; then
			AC_DEFINE(HAVE_BOOST_IOSTREAMS,,[define if the Boost::IOStreams library is available])
            BOOSTLIBDIR=`echo $BOOST_LDFLAGS | sed -e 's/@<:@^\/@:>@*//'`
            if test "x$ax_boost_user_iostreams_lib" = "x"; then
                for libextension in `ls $BOOSTLIBDIR/libboost_iostreams*.so* $BOOSTLIBDIR/libboost_iostream*.dylib* $BOOSTLIBDIR/libboost_iostreams*.a* 2>/dev/null | sed 's,.*/,,' | sed -e 's;^lib\(boost_iostreams.*\)\.so.*$;\1;' -e 's;^lib\(boost_iostream.*\)\.dylib.*$;\1;' -e 's;^lib\(boost_iostreams.*\)\.a.*$;\1;'` ; do
                     ax_lib=${libextension}
				    AC_CHECK_LIB($ax_lib, exit,
                                 [BOOST_IOSTREAMS_LIB="-l$ax_lib"; AC_SUBST(BOOST_IOSTREAMS_LIB) link_iostreams="yes"; break],
                                 [link_iostreams="no"])
				done
                if test "x$link_iostreams" != "xyes"; then
                for libextension in `ls $BOOSTLIBDIR/boost_iostreams*.dll* $BOOSTLIBDIR/boost_iostreams*.a* 2>/dev/null | sed 's,.*/,,' | sed -e 's;^\(boost_iostreams.*\)\.dll.*$;\1;' -e 's;^\(boost_iostreams.*\)\.a.*$;\1;'` ; do
                     ax_lib=${libextension}
				    AC_CHECK_LIB($ax_lib, exit,
                                 [BOOST_IOSTREAMS_LIB="-l$ax_lib"; AC_SUBST(BOOST_IOSTREAMS_LIB) link_iostreams="yes"; break],
                                 [link_iostreams="no"])
				done
                fi

            else
               for ax_lib in $ax_boost_user_iostreams_lib boost_iostreams-$ax_boost_user_iostreams_lib; do
				      AC_CHECK_LIB($ax_lib, main,
                                   [BOOST_IOSTREAMS_LIB="-l$ax_lib"; AC_SUBST(BOOST_IOSTREAMS_LIB) link_iostreams="yes"; break],
                                   [link_iostreams="no"])
                  done

            fi
            if test "x$ax_lib" = "x"; then
                AC_MSG_ERROR(Could not find a version of the library!)
            fi
			if test "x$link_iostreams" != "xyes"; then
				AC_MSG_ERROR(Could not link against $ax_lib !)
			fi
		fi

		CPPFLAGS="$CPPFLAGS_SAVED"
	LDFLAGS="$LDFLAGS_SAVED"
	fi
])

# ===========================================================================
#      https://www.gnu.org/software/autoconf-archive/ax_boost_regex.html
# ===========================================================================
#
# SYNOPSIS
#
#   AX_BOOST_REGEX
#
# DESCRIPTION
#
#   Test for Regex library from the Boost C++ libraries. The macro requires
#   a preceding call to AX_BOOST_BASE. Further documentation is available at
#   <http://randspringer.de/boost/index.html>.
#
#   This macro calls:
#
#     AC_SUBST(BOOST_REGEX_LIB)
#
#   And sets:
#
#     HAVE_BOOST_REGEX
#
# LICENSE
#
#   Copyright (c) 2008 Thomas Porschberg <thomas@randspringer.de>
#   Copyright (c) 2008 Michael Tindal
#
#   Copying and distribution of this file, with or without modification, are
#   permitted in any medium without royalty provided the copyright notice
#   and this notice are preserved. This file is offered as-is, without any
#   warranty.

#serial 23

AC_DEFUN([AX_BOOST_REGEX],
[
	AC_ARG_WITH([boost-regex],
	AS_HELP_STRING([--with-boost-regex@<:@=special-lib@:>@],
                   [use the Regex library from boost - it is possible to specify a certain library for the linker
                        e.g. --with-boost-regex=boost_regex-gcc-mt-d-1_33_1 ]),
        [
        if test "$withval" = "no"; then
			want_boost="no"
        elif test "$withval" = "yes"; then
            want_boost="yes"
            ax_boost_user_regex_lib=""
        else
		    want_boost="yes"
		ax_boost_user_regex_lib="$withval"
		fi
        ],
        [want_boost="yes"]
	)

	if test "x$want_boost" = "xyes"; then
        AC_REQUIRE([AC_PROG_CC])
		CPPFLAGS_SAVED="$CPPFLAGS"
		CPPFLAGS="$CPPFLAGS $BOOST_CPPFLAGS"
		export CPPFLAGS

		LDFLAGS_SAVED="$LDFLAGS"
		LDFLAGS="$LDFLAGS $BOOST_LDFLAGS"
		export LDFLAGS

        AC_CACHE_CHECK(whether the Boost::Regex library is available,
					   ax_cv_boost_regex,
        [AC_LANG_PUSH([C++])
			 AC_COMPILE_IFELSE([AC_LANG_PROGRAM([[@%:@include <boost/regex.hpp>
												]],
                                   [[boost::regex r(); return 0;]])],
                   ax_cv_boost_regex=yes, ax_cv_boost_regex=no)
         AC_LANG_POP([C++])
		])
		if test "x$ax_cv_boost_regex" = "xyes"; then
			AC_DEFINE(HAVE_BOOST_REGEX,,[define if the Boost::Regex library is available])
            BOOSTLIBDIR=`echo $BOOST_LDFLAGS | sed -e 's/@<:@^\/@:>@*//'`
            if test "x$ax_boost_user_regex_lib" = "x"; then
                for libextension in `ls $BOOSTLIBDIR/libboost_regex*.so* $BOOSTLIBDIR/libboost_regex*.dylib* $BOOSTLIBDIR/libboost_regex*.a* 2>/dev/null | sed 's,.*/,,' | sed -e 's;^lib\(boost_regex.*\)\.so.*$;\1;' -e 's;^lib\(boost_regex.*\)\.dylib.*;\1;' -e 's;^lib\(boost_regex.*\)\.a.*$;\1;'` ; do
                     ax_lib=${libextension}
				    AC_CHECK_LIB($ax_lib, exit,
                                 [BOOST_REGEX_LIB="-l$ax_lib"; AC_SUBST(BOOST_REGEX_LIB) link_regex="yes"; break],
                                 [link_regex="no"])
				done
                if test "x$link_regex" != "xyes"; then
                for libextension in `ls $BOOSTLIBDIR/boost_regex*.dll* $BOOSTLIBDIR/boost_regex*.a* 2>/dev/null | sed 's,.*/,,' | sed -e 's;^\(boost_regex.*\)\.dll.*$;\1;' -e 's;^\(boost_regex.*\)\.a.*$;\1;'` ; do
                     ax_lib=${libextension}
				    AC_CHECK_LIB($ax_lib, exit,
                                 [BOOST_REGEX_LIB="-l$ax_lib"; AC_SUBST(BOOST_REGEX_LIB) link_regex="yes"; break],
                                 [link_regex="no"])
				done
                fi

            else
               for ax_lib in $ax_boost_user_regex_lib boost_regex-$ax_boost_user_regex_lib; do
				      AC_CHECK_LIB($ax_lib, main,
                                   [BOOST_REGEX_LIB="-l$ax_lib"; AC_SUBST(BOOST_REGEX_LIB) link_regex="yes"; break],
                                   [link_regex="no"])
               done
            fi
            if test "x$ax_lib" = "x"; then
                AC_MSG_ERROR(Could not find a version of the Boost::Regex library!)
            fi
			if test "x$link_regex" != "xyes"; then
				AC_MSG_ERROR(Could not link against $ax_lib !)
			fi
		fi

		CPPFLAGS="$CPPFLAGS_SAVED"
	LDFLAGS="$LDFLAGS_SAVED"
	fi
])

# ===========================================================================
#  https://www.gnu.org/software/autoconf-archive/ax_check_compile_flag.html
# ===========================================================================
#
# SYNOPSIS
#
#   AX_CHECK_COMPILE_FLAG(FLAG, [ACTION-SUCCESS], [ACTION-FAILURE], [EXTRA-FLAGS], [INPUT])
#
# DESCRIPTION
#
#   Check whether the given FLAG works with the current language's compiler
#   or gives an error.  (Warnings, however, are ignored)
#
#   ACTION-SUCCESS/ACTION-FAILURE are shell commands to execute on
#   success/failure.
#
#   If EXTRA-FLAGS is defined, it is added to the current language's default
#   flags (e.g. CFLAGS) when the check is done.  The check is thus made with
#   the flags: "CFLAGS EXTRA-FLAGS FLAG".  This can for example be used to
#   force the compiler to issue an error when a bad flag is given.
#
#   INPUT gives an alternative input source to AC_COMPILE_IFELSE.
#
#   NOTE: Implementation based on AX_CFLAGS_GCC_OPTION. Please keep this
#   macro in sync with AX_CHECK_{PREPROC,LINK}_FLAG.
#
# LICENSE
#
#   Copyright (c) 2008 Guido U. Draheim <guidod@gmx.de>
#   Copyright (c) 2011 Maarten Bosmans <mkbosmans@gmail.com>
#
#   Copying and distribution of this file, with or without modification, are
#   permitted in any medium without royalty provided the copyright notice
#   and this notice are preserved.  This file is offered as-is, without any
#   warranty.

#serial 6

AC_DEFUN([AX_CHECK_COMPILE_FLAG],
[AC_PREREQ(2.64)dnl for _AC_LANG_PREFIX and AS_VAR_IF
AS_VAR_PUSHDEF([CACHEVAR],[ax_cv_check_[]_AC_LANG_ABBREV[]flags_$4_$1])dnl
AC_CACHE_CHECK([whether _AC_LANG compiler accepts $1], CACHEVAR, [
  ax_check_save_flags=$[]_AC_LANG_PREFIX[]FLAGS
  _AC_LANG_PREFIX[]FLAGS="$[]_AC_LANG_PREFIX[]FLAGS $4 $1"
  AC_COMPILE_IFELSE([m4_default([$5],[AC_LANG_PROGRAM()])],
    [AS_VAR_SET(CACHEVAR,[yes])],
    [AS_VAR_SET(CACHEVAR,[no])])
  _AC_LANG_PREFIX[]FLAGS=$ax_check_save_flags])
AS_VAR_IF(CACHEVAR,yes,
  [m4_default([$2], :)],
  [m4_default([$3], :)])
AS_VAR_POPDEF([CACHEVAR])dnl
])dnl AX_CHECK_COMPILE_FLAGS

# ===========================================================================
#     https://www.gnu.org/software/autoconf-archive/ax_check_library.html
# ===========================================================================
#
# SYNOPSIS
#
#   AX_CHECK_LIBRARY(VARIABLE-PREFIX, HEADER-FILE, LIBRARY-FILE,
#                    [ACTION-IF-FOUND], [ACTION-IF-NOT_FOUND])
#
# DESCRIPTION
#
#   Provides a generic test for a given library, similar in concept to the
#   PKG_CHECK_MODULES macro used by pkg-config.
#
#   Most simplest libraries can be checked against simply through the
#   presence of a header file and a library to link to. This macro allows to
#   wrap around the test so that it doesn't have to be recreated each time.
#
#   Rather than define --with-$LIBRARY arguments, it uses variables in the
#   same way that PKG_CHECK_MODULES does. It doesn't, though, use the same
#   names, since you shouldn't provide a value for LIBS or CFLAGS but rather
#   for LDFLAGS and CPPFLAGS, to tell the linker and compiler where to find
#   libraries and headers respectively.
#
#   If the library is find, HAVE_PREFIX is defined, and in all cases
#   PREFIX_LDFLAGS and PREFIX_CPPFLAGS are substituted.
#
#   Example:
#
#     AX_CHECK_LIBRARY([LIBEVENT], [event.h], [event], [],
#                      [AC_MSG_ERROR([Unable to find libevent])])
#
# LICENSE
#
#   Copyright (c) 2010 Diego Elio Petteno` <flameeyes@gmail.com>
#
#   This program is free software: you can redistribute it and/or modify it
#   under the terms of the GNU General Public License as published by the
#   Free Software Foundation, either version 3 of the License, or (at your
#   option) any later version.
#
#   This program is distributed in the hope that it will be useful, but
#   WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
#   Public License for more details.
#
#   You should have received a copy of the GNU General Public License along
#   with this program. If not, see <https://www.gnu.org/licenses/>.
#
#   As a special exception, the respective Autoconf Macro's copyright owner
#   gives unlimited permission to copy, distribute and modify the configure
#   scripts that are the output of Autoconf when processing the Macro. You
#   need not follow the terms of the GNU General Public License when using
#   or distributing such scripts, even though portions of the text of the
#   Macro appear in them. The GNU General Public License (GPL) does govern
#   all other use of the material that constitutes the Autoconf Macro.
#
#   This special exception to the GPL applies to versions of the Autoconf
#   Macro released by the Autoconf Archive. When you make and distribute a
#   modified version of the Autoconf Macro, you may extend this special
#   exception to the GPL to apply to your modified version as well.

#serial 5

AC_DEFUN([AX_CHECK_LIBRARY], [
  AC_ARG_VAR($1[_CPPFLAGS], [C preprocessor flags for ]$1[ headers])
  AC_ARG_VAR($1[_LDFLAGS], [linker flags for ]$1[ libraries])

  AC_CACHE_VAL(AS_TR_SH([ax_cv_have_]$1),
    [save_CPPFLAGS="$CPPFLAGS"
     save_LDFLAGS="$LDFLAGS"
     save_LIBS="$LIBS"

     AS_IF([test "x$]$1[_CPPFLAGS" != "x"],
       [CPPFLAGS="$CPPFLAGS $]$1[_CPPFLAGS"])

     AS_IF([test "x$]$1[_LDFLAGS" != "x"],
       [LDFLAGS="$LDFLAGS $]$1[_LDFLAGS"])

     AC_CHECK_HEADER($2, [
       AC_CHECK_LIB($3, [main],
         [AS_TR_SH([ax_cv_have_]$1)=yes],
         [AS_TR_SH([ax_cv_have_]$1)=no])
     ], [AS_TR_SH([ax_cv_have_]$1)=no])

     CPPFLAGS="$save_CPPFLAGS"
     LDFLAGS="$save_LDFLAGS"
     LIBS="$save_LIBS"
    ])

  AS_IF([test "$]AS_TR_SH([ax_cv_have_]$1)[" = "yes"],
    AC_DEFINE([HAVE_]$1, [1], [Define to 1 if ]$1[ is found])
    [$4],
    [$5])
])

# ===========================================================================
#  https://www.gnu.org/software/autoconf-archive/ax_cxx_compile_stdcxx.html
# ===========================================================================
#
# SYNOPSIS
#
#   AX_CXX_COMPILE_STDCXX(VERSION, [ext|noext], [mandatory|optional])
#
# DESCRIPTION
#
#   Check for baseline language coverage in the compiler for the specified
#   version of the C++ standard.  If necessary, add switches to CXX and
#   CXXCPP to enable support.  VERSION may be '11' (for the C++11 standard)
#   or '14' (for the C++14 standard).
#
#   The second argument, if specified, indicates whether you insist on an
#   extended mode (e.g. -std=gnu++11) or a strict conformance mode (e.g.
#   -std=c++11).  If neither is specified, you get whatever works, with
#   preference for an extended mode.
#
#   The third argument, if specified 'mandatory' or if left unspecified,
#   indicates that baseline support for the specified C++ standard is
#   required and that the macro should error out if no mode with that
#   support is found.  If specified 'optional', then configuration proceeds
#   regardless, after defining HAVE_CXX${VERSION} if and only if a
#   supporting mode is found.
#
# LICENSE
#
#   Copyright (c) 2008 Benjamin Kosnik <bkoz@redhat.com>
#   Copyright (c) 2012 Zack Weinberg <zackw@panix.com>
#   Copyright (c) 2013 Roy Stogner <roystgnr@ices.utexas.edu>
#   Copyright (c) 2014, 2015 Google Inc.; contributed by Alexey Sokolov <sokolov@google.com>
#   Copyright (c) 2015 Paul Norman <penorman@mac.com>
#   Copyright (c) 2015 Moritz Klammler <moritz@klammler.eu>
#   Copyright (c) 2016, 2018 Krzesimir Nowak <qdlacz@gmail.com>
#
#   Copying and distribution of this file, with or without modification, are
#   permitted in any medium without royalty provided the copyright notice
#   and this notice are preserved.  This file is offered as-is, without any
#   warranty.

#serial 10

dnl  This macro is based on the code from the AX_CXX_COMPILE_STDCXX_11 macro
dnl  (serial version number 13).

AC_DEFUN([AX_CXX_COMPILE_STDCXX], [dnl
  m4_if([$1], [11], [ax_cxx_compile_alternatives="11 0x"],
        [$1], [14], [ax_cxx_compile_alternatives="14 1y"],
        [$1], [17], [ax_cxx_compile_alternatives="17 1z"],
        [m4_fatal([invalid first argument `$1' to AX_CXX_COMPILE_STDCXX])])dnl
  m4_if([$2], [], [],
        [$2], [ext], [],
        [$2], [noext], [],
        [m4_fatal([invalid second argument `$2' to AX_CXX_COMPILE_STDCXX])])dnl
  m4_if([$3], [], [ax_cxx_compile_cxx$1_required=true],
        [$3], [mandatory], [ax_cxx_compile_cxx$1_required=true],
        [$3], [optional], [ax_cxx_compile_cxx$1_required=false],
        [m4_fatal([invalid third argument `$3' to AX_CXX_COMPILE_STDCXX])])
  AC_LANG_PUSH([C++])dnl
  ac_success=no

  m4_if([$2], [noext], [], [dnl
  if test x$ac_success = xno; then
    for alternative in ${ax_cxx_compile_alternatives}; do
      switch="-std=gnu++${alternative}"
      cachevar=AS_TR_SH([ax_cv_cxx_compile_cxx$1_$switch])
      AC_CACHE_CHECK(whether $CXX supports C++$1 features with $switch,
                     $cachevar,
        [ac_save_CXX="$CXX"
         CXX="$CXX $switch"
         AC_COMPILE_IFELSE([AC_LANG_SOURCE([_AX_CXX_COMPILE_STDCXX_testbody_$1])],
          [eval $cachevar=yes],
          [eval $cachevar=no])
         CXX="$ac_save_CXX"])
      if eval test x\$$cachevar = xyes; then
        CXX="$CXX $switch"
        if test -n "$CXXCPP" ; then
          CXXCPP="$CXXCPP $switch"
        fi
        ac_success=yes
        break
      fi
    done
  fi])

  m4_if([$2], [ext], [], [dnl
  if test x$ac_success = xno; then
    dnl HP's aCC needs +std=c++11 according to:
    dnl http://h21007.www2.hp.com/portal/download/files/unprot/aCxx/PDF_Release_Notes/769149-001.pdf
    dnl Cray's crayCC needs "-h std=c++11"
    for alternative in ${ax_cxx_compile_alternatives}; do
      for switch in -std=c++${alternative} +std=c++${alternative} "-h std=c++${alternative}"; do
        cachevar=AS_TR_SH([ax_cv_cxx_compile_cxx$1_$switch])
        AC_CACHE_CHECK(whether $CXX supports C++$1 features with $switch,
                       $cachevar,
          [ac_save_CXX="$CXX"
           CXX="$CXX $switch"
           AC_COMPILE_IFELSE([AC_LANG_SOURCE([_AX_CXX_COMPILE_STDCXX_testbody_$1])],
            [eval $cachevar=yes],
            [eval $cachevar=no])
           CXX="$ac_save_CXX"])
        if eval test x\$$cachevar = xyes; then
          CXX="$CXX $switch"
          if test -n "$CXXCPP" ; then
            CXXCPP="$CXXCPP $switch"
          fi
          ac_success=yes
          break
        fi
      done
      if test x$ac_success = xyes; then
        break
      fi
    done
  fi])
  AC_LANG_POP([C++])
  if test x$ax_cxx_compile_cxx$1_required = xtrue; then
    if test x$ac_success = xno; then
      AC_MSG_ERROR([*** A compiler with support for C++$1 language features is required.])
    fi
  fi
  if test x$ac_success = xno; then
    HAVE_CXX$1=0
    AC_MSG_NOTICE([No compiler with C++$1 support was found])
  else
    HAVE_CXX$1=1
    AC_DEFINE(HAVE_CXX$1,1,
              [define if the compiler supports basic C++$1 syntax])
  fi
  AC_SUBST(HAVE_CXX$1)
])


dnl  Test body for checking C++11 support

m4_define([_AX_CXX_COMPILE_STDCXX_testbody_11],
  _AX_CXX_COMPILE_STDCXX_testbody_new_in_11
)


dnl  Test body for checking C++14 support

m4_define([_AX_CXX_COMPILE_STDCXX_testbody_14],
  _AX_CXX_COMPILE_STDCXX_testbody_new_in_11
  _AX_CXX_COMPILE_STDCXX_testbody_new_in_14
)

m4_define([_AX_CXX_COMPILE_STDCXX_testbody_17],
  _AX_CXX_COMPILE_STDCXX_testbody_new_in_11
  _AX_CXX_COMPILE_STDCXX_testbody_new_in_14
  _AX_CXX_COMPILE_STDCXX_testbody_new_in_17
)

dnl  Tests for new features in C++11

m4_define([_AX_CXX_COMPILE_STDCXX_testbody_new_in_11], [[

// If the compiler admits that it is not ready for C++11, why torture it?
// Hopefully, this will speed up the test.

#ifndef __cplusplus

#error "This is not a C++ compiler"

#elif __cplusplus < 201103L

#error "This is not a C++11 compiler"

#else

namespace cxx11
{

  namespace test_static_assert
  {

    template <typename T>
    struct check
    {
      static_assert(sizeof(int) <= sizeof(T), "not big enough");
    };

  }

  namespace test_final_override
  {

    struct Base
    {
      virtual void f() {}
    };

    struct Derived : public Base
    {
      virtual void f() override {}
    };

  }

  namespace test_double_right_angle_brackets
  {

    template < typename T >
    struct check {};

    typedef check<void> single_type;
    typedef check<check<void>> double_type;
    typedef check<check<check<void>>> triple_type;
    typedef check<check<check<check<void>>>> quadruple_type;

  }

  namespace test_decltype
  {

    int
    f()
    {
      int a = 1;
      decltype(a) b = 2;
      return a + b;
    }

  }

  namespace test_type_deduction
  {

    template < typename T1, typename T2 >
    struct is_same
    {
      static const bool value = false;
    };

    template < typename T >
    struct is_same<T, T>
    {
      static const bool value = true;
    };

    template < typename T1, typename T2 >
    auto
    add(T1 a1, T2 a2) -> decltype(a1 + a2)
    {
      return a1 + a2;
    }

    int
    test(const int c, volatile int v)
    {
      static_assert(is_same<int, decltype(0)>::value == true, "");
      static_assert(is_same<int, decltype(c)>::value == false, "");
      static_assert(is_same<int, decltype(v)>::value == false, "");
      auto ac = c;
      auto av = v;
      auto sumi = ac + av + 'x';
      auto sumf = ac + av + 1.0;
      static_assert(is_same<int, decltype(ac)>::value == true, "");
      static_assert(is_same<int, decltype(av)>::value == true, "");
      static_assert(is_same<int, decltype(sumi)>::value == true, "");
      static_assert(is_same<int, decltype(sumf)>::value == false, "");
      static_assert(is_same<int, decltype(add(c, v))>::value == true, "");
      return (sumf > 0.0) ? sumi : add(c, v);
    }

  }

  namespace test_noexcept
  {

    int f() { return 0; }
    int g() noexcept { return 0; }

    static_assert(noexcept(f()) == false, "");
    static_assert(noexcept(g()) == true, "");

  }

  namespace test_constexpr
  {

    template < typename CharT >
    unsigned long constexpr
    strlen_c_r(const CharT *const s, const unsigned long acc) noexcept
    {
      return *s ? strlen_c_r(s + 1, acc + 1) : acc;
    }

    template < typename CharT >
    unsigned long constexpr
    strlen_c(const CharT *const s) noexcept
    {
      return strlen_c_r(s, 0UL);
    }

    static_assert(strlen_c("") == 0UL, "");
    static_assert(strlen_c("1") == 1UL, "");
    static_assert(strlen_c("example") == 7UL, "");
    static_assert(strlen_c("another\0example") == 7UL, "");

  }

  namespace test_rvalue_references
  {

    template < int N >
    struct answer
    {
      static constexpr int value = N;
    };

    answer<1> f(int&)       { return answer<1>(); }
    answer<2> f(const int&) { return answer<2>(); }
    answer<3> f(int&&)      { return answer<3>(); }

    void
    test()
    {
      int i = 0;
      const int c = 0;
      static_assert(decltype(f(i))::value == 1, "");
      static_assert(decltype(f(c))::value == 2, "");
      static_assert(decltype(f(0))::value == 3, "");
    }

  }

  namespace test_uniform_initialization
  {

    struct test
    {
      static const int zero {};
      static const int one {1};
    };

    static_assert(test::zero == 0, "");
    static_assert(test::one == 1, "");

  }

  namespace test_lambdas
  {

    void
    test1()
    {
      auto lambda1 = [](){};
      auto lambda2 = lambda1;
      lambda1();
      lambda2();
    }

    int
    test2()
    {
      auto a = [](int i, int j){ return i + j; }(1, 2);
      auto b = []() -> int { return '0'; }();
      auto c = [=](){ return a + b; }();
      auto d = [&](){ return c; }();
      auto e = [a, &b](int x) mutable {
        const auto identity = [](int y){ return y; };
        for (auto i = 0; i < a; ++i)
          a += b--;
        return x + identity(a + b);
      }(0);
      return a + b + c + d + e;
    }

    int
    test3()
    {
      const auto nullary = [](){ return 0; };
      const auto unary = [](int x){ return x; };
      using nullary_t = decltype(nullary);
      using unary_t = decltype(unary);
      const auto higher1st = [](nullary_t f){ return f(); };
      const auto higher2nd = [unary](nullary_t f1){
        return [unary, f1](unary_t f2){ return f2(unary(f1())); };
      };
      return higher1st(nullary) + higher2nd(nullary)(unary);
    }

  }

  namespace test_variadic_templates
  {

    template <int...>
    struct sum;

    template <int N0, int... N1toN>
    struct sum<N0, N1toN...>
    {
      static constexpr auto value = N0 + sum<N1toN...>::value;
    };

    template <>
    struct sum<>
    {
      static constexpr auto value = 0;
    };

    static_assert(sum<>::value == 0, "");
    static_assert(sum<1>::value == 1, "");
    static_assert(sum<23>::value == 23, "");
    static_assert(sum<1, 2>::value == 3, "");
    static_assert(sum<5, 5, 11>::value == 21, "");
    static_assert(sum<2, 3, 5, 7, 11, 13>::value == 41, "");

  }

  // http://stackoverflow.com/questions/13728184/template-aliases-and-sfinae
  // Clang 3.1 fails with headers of libstd++ 4.8.3 when using std::function
  // because of this.
  namespace test_template_alias_sfinae
  {

    struct foo {};

    template<typename T>
    using member = typename T::member_type;

    template<typename T>
    void func(...) {}

    template<typename T>
    void func(member<T>*) {}

    void test();

    void test() { func<foo>(0); }

  }

}  // namespace cxx11

#endif  // __cplusplus >= 201103L

]])


dnl  Tests for new features in C++14

m4_define([_AX_CXX_COMPILE_STDCXX_testbody_new_in_14], [[

// If the compiler admits that it is not ready for C++14, why torture it?
// Hopefully, this will speed up the test.

#ifndef __cplusplus

#error "This is not a C++ compiler"

#elif __cplusplus < 201402L

#error "This is not a C++14 compiler"

#else

namespace cxx14
{

  namespace test_polymorphic_lambdas
  {

    int
    test()
    {
      const auto lambda = [](auto&&... args){
        const auto istiny = [](auto x){
          return (sizeof(x) == 1UL) ? 1 : 0;
        };
        const int aretiny[] = { istiny(args)... };
        return aretiny[0];
      };
      return lambda(1, 1L, 1.0f, '1');
    }

  }

  namespace test_binary_literals
  {

    constexpr auto ivii = 0b0000000000101010;
    static_assert(ivii == 42, "wrong value");

  }

  namespace test_generalized_constexpr
  {

    template < typename CharT >
    constexpr unsigned long
    strlen_c(const CharT *const s) noexcept
    {
      auto length = 0UL;
      for (auto p = s; *p; ++p)
        ++length;
      return length;
    }

    static_assert(strlen_c("") == 0UL, "");
    static_assert(strlen_c("x") == 1UL, "");
    static_assert(strlen_c("test") == 4UL, "");
    static_assert(strlen_c("another\0test") == 7UL, "");

  }

  namespace test_lambda_init_capture
  {

    int
    test()
    {
      auto x = 0;
      const auto lambda1 = [a = x](int b){ return a + b; };
      const auto lambda2 = [a = lambda1(x)](){ return a; };
      return lambda2();
    }

  }

  namespace test_digit_separators
  {

    constexpr auto ten_million = 100'000'000;
    static_assert(ten_million == 100000000, "");

  }

  namespace test_return_type_deduction
  {

    auto f(int& x) { return x; }
    decltype(auto) g(int& x) { return x; }

    template < typename T1, typename T2 >
    struct is_same
    {
      static constexpr auto value = false;
    };

    template < typename T >
    struct is_same<T, T>
    {
      static constexpr auto value = true;
    };

    int
    test()
    {
      auto x = 0;
      static_assert(is_same<int, decltype(f(x))>::value, "");
      static_assert(is_same<int&, decltype(g(x))>::value, "");
      return x;
    }

  }

}  // namespace cxx14

#endif  // __cplusplus >= 201402L

]])


dnl  Tests for new features in C++17

m4_define([_AX_CXX_COMPILE_STDCXX_testbody_new_in_17], [[

// If the compiler admits that it is not ready for C++17, why torture it?
// Hopefully, this will speed up the test.

#ifndef __cplusplus

#error "This is not a C++ compiler"

#elif __cplusplus < 201703L

#error "This is not a C++17 compiler"

#else

#include <initializer_list>
#include <utility>
#include <type_traits>

namespace cxx17
{

  namespace test_constexpr_lambdas
  {

    constexpr int foo = [](){return 42;}();

  }

  namespace test::nested_namespace::definitions
  {

  }

  namespace test_fold_expression
  {

    template<typename... Args>
    int multiply(Args... args)
    {
      return (args * ... * 1);
    }

    template<typename... Args>
    bool all(Args... args)
    {
      return (args && ...);
    }

  }

  namespace test_extended_static_assert
  {

    static_assert (true);

  }

  namespace test_auto_brace_init_list
  {

    auto foo = {5};
    auto bar {5};

    static_assert(std::is_same<std::initializer_list<int>, decltype(foo)>::value);
    static_assert(std::is_same<int, decltype(bar)>::value);
  }

  namespace test_typename_in_template_template_parameter
  {

    template<template<typename> typename X> struct D;

  }

  namespace test_fallthrough_nodiscard_maybe_unused_attributes
  {

    int f1()
    {
      return 42;
    }

    [[nodiscard]] int f2()
    {
      [[maybe_unused]] auto unused = f1();

      switch (f1())
      {
      case 17:
        f1();
        [[fallthrough]];
      case 42:
        f1();
      }
      return f1();
    }

  }

  namespace test_extended_aggregate_initialization
  {

    struct base1
    {
      int b1, b2 = 42;
    };

    struct base2
    {
      base2() {
        b3 = 42;
      }
      int b3;
    };

    struct derived : base1, base2
    {
        int d;
    };

    derived d1 {{1, 2}, {}, 4};  // full initialization
    derived d2 {{}, {}, 4};      // value-initialized bases

  }

  namespace test_general_range_based_for_loop
  {

    struct iter
    {
      int i;

      int& operator* ()
      {
        return i;
      }

      const int& operator* () const
      {
        return i;
      }

      iter& operator++()
      {
        ++i;
        return *this;
      }
    };

    struct sentinel
    {
      int i;
    };

    bool operator== (const iter& i, const sentinel& s)
    {
      return i.i == s.i;
    }

    bool operator!= (const iter& i, const sentinel& s)
    {
      return !(i == s);
    }

    struct range
    {
      iter begin() const
      {
        return {0};
      }

      sentinel end() const
      {
        return {5};
      }
    };

    void f()
    {
      range r {};

      for (auto i : r)
      {
        [[maybe_unused]] auto v = i;
      }
    }

  }

  namespace test_lambda_capture_asterisk_this_by_value
  {

    struct t
    {
      int i;
      int foo()
      {
        return [*this]()
        {
          return i;
        }();
      }
    };

  }

  namespace test_enum_class_construction
  {

    enum class byte : unsigned char
    {};

    byte foo {42};

  }

  namespace test_constexpr_if
  {

    template <bool cond>
    int f ()
    {
      if constexpr(cond)
      {
        return 13;
      }
      else
      {
        return 42;
      }
    }

  }

  namespace test_selection_statement_with_initializer
  {

    int f()
    {
      return 13;
    }

    int f2()
    {
      if (auto i = f(); i > 0)
      {
        return 3;
      }

      switch (auto i = f(); i + 4)
      {
      case 17:
        return 2;

      default:
        return 1;
      }
    }

  }

  namespace test_template_argument_deduction_for_class_templates
  {

    template <typename T1, typename T2>
    struct pair
    {
      pair (T1 p1, T2 p2)
        : m1 {p1},
          m2 {p2}
      {}

      T1 m1;
      T2 m2;
    };

    void f()
    {
      [[maybe_unused]] auto p = pair{13, 42u};
    }

  }

  namespace test_non_type_auto_template_parameters
  {

    template <auto n>
    struct B
    {};

    B<5> b1;
    B<'a'> b2;

  }

  namespace test_structured_bindings
  {

    int arr[2] = { 1, 2 };
    std::pair<int, int> pr = { 1, 2 };

    auto f1() -> int(&)[2]
    {
      return arr;
    }

    auto f2() -> std::pair<int, int>&
    {
      return pr;
    }

    struct S
    {
      int x1 : 2;
      volatile double y1;
    };

    S f3()
    {
      return {};
    }

    auto [ x1, y1 ] = f1();
    auto& [ xr1, yr1 ] = f1();
    auto [ x2, y2 ] = f2();
    auto& [ xr2, yr2 ] = f2();
    const auto [ x3, y3 ] = f3();

  }

  namespace test_exception_spec_type_system
  {

    struct Good {};
    struct Bad {};

    void g1() noexcept;
    void g2();

    template<typename T>
    Bad
    f(T*, T*);

    template<typename T1, typename T2>
    Good
    f(T1*, T2*);

    static_assert (std::is_same_v<Good, decltype(f(g1, g2))>);

  }

  namespace test_inline_variables
  {

    template<class T> void f(T)
    {}

    template<class T> inline T g(T)
    {
      return T{};
    }

    template<> inline void f<>(int)
    {}

    template<> int g<>(int)
    {
      return 5;
    }

  }

}  // namespace cxx17

#endif  // __cplusplus < 201703L

]])

# =============================================================================
#  https://www.gnu.org/software/autoconf-archive/ax_cxx_compile_stdcxx_17.html
# =============================================================================
#
# SYNOPSIS
#
#   AX_CXX_COMPILE_STDCXX_17([ext|noext], [mandatory|optional])
#
# DESCRIPTION
#
#   Check for baseline language coverage in the compiler for the C++17
#   standard; if necessary, add switches to CXX and CXXCPP to enable
#   support.
#
#   This macro is a convenience alias for calling the AX_CXX_COMPILE_STDCXX
#   macro with the version set to C++17.  The two optional arguments are
#   forwarded literally as the second and third argument respectively.
#   Please see the documentation for the AX_CXX_COMPILE_STDCXX macro for
#   more information.  If you want to use this macro, you also need to
#   download the ax_cxx_compile_stdcxx.m4 file.
#
# LICENSE
#
#   Copyright (c) 2015 Moritz Klammler <moritz@klammler.eu>
#   Copyright (c) 2016 Krzesimir Nowak <qdlacz@gmail.com>
#
#   Copying and distribution of this file, with or without modification, are
#   permitted in any medium without royalty provided the copyright notice
#   and this notice are preserved. This file is offered as-is, without any
#   warranty.

#serial 2

AX_REQUIRE_DEFINED([AX_CXX_COMPILE_STDCXX])
AC_DEFUN([AX_CXX_COMPILE_STDCXX_17], [AX_CXX_COMPILE_STDCXX([17], [$1], [$2])])

m4_include([config/m4/libtool.m4])
m4_include([config/m4/ltoptions.m4])
m4_include([config/m4/ltsugar.m4])
m4_include([config/m4/ltversion.m4])
m4_include([config/m4/lt~obsolete.m4])
