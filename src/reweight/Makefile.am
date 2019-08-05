include $(top_srcdir)/build.common

AM_CXXFLAGS = -DNREWEIGHT_VERSION=@NEUT_VERSION_NUM@ $(DEFAULT_CXXFLAGS) $(NEUT_CLASS_INCLUDES)

LIBNAME = lib@LIBPREFIX@NReWeight.@NEUT_VERSION_MAJ@.@NEUT_VERSION_MIN@.@NEUT_VERSION_REV@.a

lib_LIBRARIES = $(LIBNAME)

lib@LIBPREFIX@NReWeight_@NEUT_VERSION_MAJ@_@NEUT_VERSION_MIN@_@NEUT_VERSION_REV@_a_SOURCES = NFortFns.cc NModeDefn.cc NReWeightCasc.cc NReWeight.cc NReWeightINuke.cc NReWeightNuclPiless.cc NReWeightNuXSecCCQE.cc NReWeightNuXSecCCRES.cc NReWeightNuXSecCOH.cc NReWeightNuXSecDIS.cc NReWeightNuXSecNC.cc NReWeightNuXSecNCEL.cc NReWeightNuXSecNCRES.cc NReWeightNuXSecRES.cc NReWeightUtils.cc NSystSet.cc NSystUncertainty.cc NTotCrs.cc

include_HEADERS = Controls.h NFortFns.h NModeDefn.h NReWeightCasc.h NReWeight.h NReWeightI.h NReWeightINuke.h NReWeightNuclPiless.h NReWeightNuXSecCCQE.h NReWeightNuXSecCCRES.h NReWeightNuXSecCOH.h NReWeightNuXSecDIS.h NReWeightNuXSecNCEL.h NReWeightNuXSecNC.h NReWeightNuXSecNCRES.h NReWeightNuXSecRES.h NReWeightUtils.h NSyst.h NSystSet.h NSystUncertainty.h NTotCrs.h PDGCodes.h

# Unused sources that live in the directory
UNUSED_srcs =

install-exec-hook:
	-(cd ${libdir}; $(LN_S) $(LIBNAME) $(shell echo $(LIBNAME) | cut -d . -f 1).a)

uninstall-hook:
	-(cd ${libdir}; rm -f $(shell echo $(LIBNAME) | cut -d . -f 1).a)
