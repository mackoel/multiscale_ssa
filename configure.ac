dnl Process this file with autoconf to produce a configure script.
dnl Created by Anjuta application wizard.

AC_INIT(multiscale_ssa, 0.1)

AM_INIT_AUTOMAKE(AC_PACKAGE_NAME, AC_PACKAGE_VERSION)
AC_CONFIG_HEADERS([config.h])
AM_MAINTAINER_MODE

dnl AC_CONFIG_HEADERS([config.h])

dnl AM_INIT_AUTOMAKE([1.11])

dnl AM_SILENT_RULES([yes])

dnl AC_PROG_CC

AC_ISC_POSIX
AC_PROG_CC
AM_PROG_CC_STDC
AC_HEADER_STDC


dnl ***************************************************************************
dnl Internationalization
dnl ***************************************************************************
IT_PROG_INTLTOOL([0.35.0])

GETTEXT_PACKAGE=multiscale_ssa
AC_SUBST(GETTEXT_PACKAGE)
AC_DEFINE_UNQUOTED(GETTEXT_PACKAGE,"$GETTEXT_PACKAGE", [GETTEXT package name])

AM_GLIB_GNU_GETTEXT



dnl ***************************************************************************
dnl Check for Windows
dnl ***************************************************************************
AC_CANONICAL_HOST

case $host_os in
  *mingw*)
    platform_win32=yes
    native_win32=yes
    ;;
  pw32* | *cygwin*)
    platform_win32=yes
    native_win32=no
    ;;
  *)
    platform_win32=no
    native_win32=no
    ;;
esac
AM_CONDITIONAL(PLATFORM_WIN32, test x"$platform_win32" = "xyes")
AM_CONDITIONAL(NATIVE_WIN32, test x"$native_win32" = "xyes")


LT_INIT([win32-dll])


PKG_CHECK_MODULES(MULTISCALE_SSA, [glib-2.0 ])




AC_OUTPUT([
Makefile
src/Makefile
po/Makefile.in
])
