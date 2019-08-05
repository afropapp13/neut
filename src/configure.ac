AC_PREREQ([2.68])

AC_INIT([neut], [5.4.1])

AC_SUBST(NEUT_VERSION, AC_PACKAGE_VERSION)
AC_SUBST(NEUT_VERSION_NUM, $(echo AC_PACKAGE_VERSION | tr -d '.' ))
AC_SUBST(NEUT_VERSION_MAJ, $(echo AC_PACKAGE_VERSION | cut -d . -f 1 ))
AC_SUBST(NEUT_VERSION_MIN, $(echo AC_PACKAGE_VERSION | cut -d . -f 2 ))
AC_SUBST(NEUT_VERSION_REV, $(echo AC_PACKAGE_VERSION | cut -d . -f 3 ))

AC_CONFIG_SRCDIR([neutsmpl/neutroot.F])
AC_CONFIG_MACRO_DIRS([m4])

################ Check CERNLIB
AC_ARG_VAR(CERN, Top level directory of build cernlib)
AC_ARG_VAR(CERN_LEVEL, [The CERN 'level', almost certainly 2005])

if test x"${CERN}" == x; then
  AC_MSG_ERROR(CERN environment variable is not defined.)
else
  echo "CERN = ${CERN}"
fi

if test x"${CERN_LEVEL}" == x; then
  AC_MSG_ERROR(CERN_LEVEL environment variable is not defined (likely should be 2005).)
else
  echo "CERN_LEVEL = ${CERN_LEVEL}"
fi

AC_CHECK_FILE("${CERN}/${CERN_LEVEL}/bin/cernlib",
  [],
  [AC_MSG_ERROR("${CERN}/${CERN_LEVEL}/bin/cernlib does not exist. Is CERNLIB built correctly?")]
)

AC_CHECK_FILE("${CERN}/${CERN_LEVEL}/include/pdf804/pilot.h",
  [],
  [AC_MSG_ERROR("pdf804/pilot.h does not exist in ${CERN}/${CERN_LEVEL}/include/ (expecting CERN and CERN_LEVEL to be defined.)")]
)

AM_INIT_AUTOMAKE([foreign 1.13 subdir-objects])

AC_PROG_CC
AC_PROG_CXX
AC_PROG_FC
AC_PROG_F77
AC_PROG_LIBTOOL
AC_PROG_LN_S

LT_INIT

########## Check ROOT
AC_CHECK_PROG(FOUND_ROOTCONFIG,root-config,"FOUND","NOT-FOUND")

case "${FOUND_ROOTCONFIG}" in
  "FOUND") AC_MSG_RESULT(Found root-config.);;
  *) AC_MSG_ERROR(Failed to find root-config. Is root set up?);
esac

ROOT_MAJOR_VERSION=$(root-config --version  | cut -d . -f 1)
if test ${ROOT_MAJOR_VERSION} = 6; then
  echo "Using ROOT 6."
else
  echo "Assuming ROOT 5."
fi
AM_CONDITIONAL([ROOT_MAJOR_IS_SIX], [test ${ROOT_MAJOR_VERSION} = 6])

######## Check compat mode
AC_ARG_ENABLE([compatnames],
    AS_HELP_STRING([--enable-compatnames], [Enable 5.4.0 and previous library naming scheme]))

AM_CONDITIONAL([NEUT_LIBCOMPAT], [test "x$enable_compatnames" = "xyes"])

AS_IF([test "x$enable_compatnames" = "xyes"], [
  AC_SUBST(LIBPREFIX, [])
], [
  AC_SUBST(LIBPREFIX, [NEUT])
])

AC_CONFIG_FILES([
  Makefile
  neutcore/Makefile
  nuccorspl/Makefile
  nuceff/Makefile
  n1p1h/Makefile
  ht2p2h/Makefile
  radcorr/Makefile
  tauola/Makefile
  specfunc/Makefile
  skmcsvc/Makefile
  zbsfns/Makefile
  partnuck/Makefile
  neutclass/Makefile
  neutsmpl/Makefile
  nuclsmpl/Makefile
  reweight/Makefile
  tools/Makefile
  neututils/Makefile
  neututils/calc_ccqe_crs/Makefile
  neututils/calc_coh_crs.rs/Makefile
  neututils/calc_dis_crs/Makefile
  neututils/calc_nuc_dens/Makefile
  neututils/calc_res_crs/Makefile
])

AC_OUTPUT
