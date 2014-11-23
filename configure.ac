#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.65])
AC_INIT([RSTools-playground], [0.9.5], [andre.hoffmann@meduniwien.ac.at])

ACLOCAL_AMFLAGS=-I m4
AC_CONFIG_MACRO_DIR([m4])
AC_CONFIG_AUX_DIR([m4])

CFLAGS="-ansi -pedantic -Woverlength-strings -Wno-deprecated-register -g -std=c99 -fopenmp $CFLAGS"
CXXFLAGS="-ansi -pedantic -Woverlength-strings -Wno-deprecated-register -g -std=c++0x -fopenmp $CXXFLAGS"

# prepare for automake
AM_INIT_AUTOMAKE([foreign -Wall])

# Checks for programs.
AC_PROG_CC
AC_PROG_CXX
AM_PROG_AR
AC_PROG_RANLIB
#AC_LIBTOOL_DLOPEN
AC_ENABLE_SHARED
AC_DISABLE_STATIC
AC_PROG_LIBTOOL
AC_PROG_INSTALL
LT_INIT

# Check if on OS X
AC_CANONICAL_HOST
AM_CONDITIONAL(BUILD_OS_IS_DARWIN, [test x"$build_os" = darwin])

# determine git version information
RSTOOLS_VERSION_HASH=esyscmd([git --git-dir=./.git --work-tree=. rev-parse --short HEAD])
RSTOOLS_VERSION_DATE="esyscmd([git --git-dir=./.git --work-tree=. log -1 --pretty=format:'%ci'])"

AC_DEFINE_UNQUOTED(RSTOOLS_VERSION_HASH,
                   ["$RSTOOLS_VERSION_HASH"],
                   [A string containing the Git hash of the last commit]
)

AC_DEFINE_UNQUOTED(RSTOOLS_VERSION_DATE,
                   ["$RSTOOLS_VERSION_DATE"],
                   [A string containing the date of the last commit]
)

AC_DEFINE_UNQUOTED(GETTEXT_PACKAGE,
                   ["rstools"],
                   [The translation package for libgettext(unused)]
)

# Checks for header files.
AC_CHECK_HEADERS([float.h string.h strings.h])

# Check for RSTools which will already come with all the necessary dependencies
PKG_CHECK_MODULES([RSTOOLS], [rstools >= 0.9.5])
AC_SUBST([RSTOOLS_CFLAGS])
AC_SUBST([RSTOOLS_LIBS])
LIBS=["$RSTOOLS_LIBS $LIBS"]
CFLAGS=["$RSTOOLS_CFLAGS $CFLAGS"]
CXXFLAGS=["$RSTOOLS_CFLAGS $CXXFLAGS"]

# Checks for typedefs, structures, and compiler characteristics.
AC_C_INLINE
AC_TYPE_SIZE_T
AC_TYPE_SSIZE_T

# Checks for library functions.
#AC_FUNC_MALLOC
#AC_FUNC_REALLOC
AC_CHECK_FUNCS([floor pow sqrt])

# Define included makefiles
AC_CONFIG_SRCDIR([src])
AC_CONFIG_HEADERS([config.h])
AC_CONFIG_FILES([Makefile
		 src/Makefile
 		 src/playground/Makefile])

AC_OUTPUT
