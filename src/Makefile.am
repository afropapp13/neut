%C.h: %.h
	$(top_srcdir)/tools/fh2h.pl $< > $@

n1p1h/nieves1p1hC.h: n1p1h/inc/nieves1p1h.h
	$(top_srcdir)/tools/fh2h.pl $< > $@

neutcore/neutmodelC.h: neutcore/neutmodel.h n1p1h/nieves1p1hC.h
	$(top_srcdir)/tools/fh2h.pl $< > $@

BUILT_SOURCES = neutcore/necardC.h neutcore/neworkC.h neutcore/neutfilepathC.h neutcore/posinnucC.h neutcore/neutparamsC.h neutcore/neutmodelC.h neutcore/nefillverC.h neutcore/neutcrsC.h neutcore/rsconsC.h skmcsvc/vcworkC.h skmcsvc/vcvrtxC.h n1p1h/nieves1p1hC.h nuccorspl/nrcardC.h nuceff/fsihistC.h nuccorspl/nucleonfsihistC.h
CLEANFILES = $(BUILT_SOURCES)

include_HEADERS = $(BUILT_SOURCES)

SUBDIRS = neutcore nuccorspl nuceff n1p1h ht2p2h radcorr tauola specfunc skmcsvc zbsfns partnuck neutclass neutsmpl reweight tools
