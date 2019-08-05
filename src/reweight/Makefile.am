include $(top_srcdir)/build.common

AM_CXXFLAGS = -DNREWEIGHT_VERSION=@NEUT_VERSION_NUM@ $(DEFAULT_CXXFLAGS) $(NEUT_CLASS_INCLUDES)

if NEUT_LIBCOMPAT
LIBNAME = libNReWeight.@NEUT_VERSION_NUM@.a
else
LIBNAME = libNEUTReWeight.@NEUT_VERSION_NUM@.a
endif

lib_LIBRARIES = $(LIBNAME)

# Don't seem to be able to have canonical names easily build from makefile variables, this is ugly, but works

SRCLIST = NFortFns.cc NModeDefn.cc NReWeightCasc.cc NReWeight.cc NReWeightINuke.cc NReWeightNuclPiless.cc NReWeightNuXSecCCQE.cc NReWeightNuXSecCCRES.cc NReWeightNuXSecCOH.cc NReWeightNuXSecDIS.cc NReWeightNuXSecNC.cc NReWeightNuXSecNCEL.cc NReWeightNuXSecNCRES.cc NReWeightNuXSecRES.cc NReWeightUtils.cc NSystSet.cc NSystUncertainty.cc NTotCrs.cc

if NEUT_LIBCOMPAT
libNReWeight_@NEUT_VERSION_NUM@_a_SOURCES = $(SRCLIST)
else
libNEUTReWeight_@NEUT_VERSION_NUM@_a_SOURCES = $(SRCLIST)
endif

include_HEADERS = NReWeightControls.h NFortFns.h NModeDefn.h NReWeightCasc.h NReWeight.h NReWeightI.h NReWeightINuke.h NReWeightNuclPiless.h NReWeightNuXSecCCQE.h NReWeightNuXSecCCRES.h NReWeightNuXSecCOH.h NReWeightNuXSecDIS.h NReWeightNuXSecNCEL.h NReWeightNuXSecNC.h NReWeightNuXSecNCRES.h NReWeightNuXSecRES.h NReWeightUtils.h NSyst.h NSystSet.h NSystUncertainty.h NTotCrs.h PDGCodes.h

# Unused sources that live in the directory
UNUSED_srcs =

install-exec-hook:
	-(cd ${libdir}; $(LN_S) $(LIBNAME) $(shell echo $(LIBNAME) | cut -d . -f 1).a)

uninstall-hook:
	-(cd ${libdir}; rm -f $(shell echo $(LIBNAME) | cut -d . -f 1).a)
