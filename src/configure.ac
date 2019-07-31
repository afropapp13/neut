AC_PREREQ([2.68])

AC_INIT([neut], [5.4.1])

AC_SUBST(NEUT_VERSION, AC_PACKAGE_VERSION)
AC_SUBST(NEUT_VERSION_NUM, $(echo AC_PACKAGE_VERSION | tr -d '.' ))

AC_CONFIG_SRCDIR([neutsmpl/neutroot.F])

AC_ARG_VAR(CERN, Top level directory of build cernlib)
AC_ARG_VAR(CERN_LEVEL, [The CERN 'level', almost certainly 2005])

AC_CHECK_FILE("${CERN}/${CERN_LEVEL}/include/pdf804/pilot.h",
  [],
  [AC_MSG_ERROR("pdf804/pilot.h does not exist in ${CERN}/${CERN_LEVEL}/include/ (expecting CERN and CERN_LEVEL to be defined.)")]
)

AM_INIT_AUTOMAKE([foreign 1.13 subdir-objects])

AC_PROG_CC
AC_PROG_CXX
AC_PROG_FC
AC_PROG_F77
AC_PROG_RANLIB

AC_CHECK_PROG(FOUND_ROOTCONFIG,root-config,"FOUND","NOT-FOUND")

case "${FOUND_ROOTCONFIG}" in
  "FOUND") AC_MSG_RESULT(Found foot-config.);;
  *) AC_MSG_ERROR(Failed to find root-config. Is root set up?);;
esac

AC_CONFIG_FILES([
  Makefile
  neutcore/Makefile
  nuccorspl/Makefile
  nuceff/Makefile
  n1p1h/Makefile
  ht2p2h/Makefile
  radcorr/Makefile
])

AC_OUTPUT
