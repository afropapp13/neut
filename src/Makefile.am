include build.common

F77 = $(FC)
FFLAGS = $(FCFLAGS)

%C.h: %.h
	$(top_srcdir)/tools/fh2h.pl $< > $@

n1p1h/nieves1p1hC.h: n1p1h/inc/nieves1p1h.h
	$(top_srcdir)/tools/fh2h.pl $< > $@

neutcore/neutmodelC.h: neutcore/neutmodel.h n1p1h/nieves1p1hC.h
	$(top_srcdir)/tools/fh2h.pl $< > $@

BUILT_SOURCES = neutcore/cardnameC.h neutcore/necardC.h neutcore/neworkC.h neutcore/neutfilepathC.h neutcore/posinnucC.h neutcore/neutparamsC.h neutcore/neutmodelC.h neutcore/nefillverC.h neutcore/neutcrsC.h neutcore/rsconsC.h skmcsvc/necardevC.h skmcsvc/vcworkC.h skmcsvc/vcvrtxC.h n1p1h/nieves1p1hC.h nuccorspl/nrcardC.h nuceff/fsihistC.h nuceff/efpionC.h nuccorspl/nucleonfsihistC.h
CLEANFILES = $(BUILT_SOURCES)

include_HEADERS = $(BUILT_SOURCES)

if NEUT_LIBCOMPAT
# Order here matters, should build core libraries, combine them in this dir (.), then build reweight, programs, and tools.
# Cannot currently build neututils dirs in compat mode.
SUBDIRS = neutcore nuccorspl nuceff n1p1h ht2p2h radcorr skmcsvc partnuck neutclass specfunc . neutsmpl tools
else
# Order here matters, should build core libraries, combine them in this dir (.), then build reweight, programs, and tools.
SUBDIRS = neutcore nuccorspl nuceff n1p1h mec ht2p2h radcorr skmcsvc partnuck neutclass specfunc . neutsmpl tools neututils nuclsmpl pionsmpl neutgeom
endif

SUBDIRS += reweight

SUBDIRS += test

# When I can get a hold of the zebra code I will add an --enable flag for the ZEBRA_SUBDIRS, but for the moment, disable them.
ZEBRA_SUBDIRS = zbsfns

if !NEUT_LIBCOMPAT

# Combine all of the core static libraries into one!
# The lines below trick/override automake into building this library out of all
# of the object files it used to build the other libraries.
# This will cause a make warning about overriding a target, but you can safetly
# ignore the warning.
LIBNAME = libNEUT.@NEUT_VERSION_NUM@.a

lib_LIBRARIES = $(LIBNAME)

libNEUT_@NEUT_VERSION_NUM@_a_SOURCES =

$(LIBNAME): $(NEUT_STATIC_LIBS)
	$(LD) -r -o neut.o --whole-archive --start-group $(NEUT_STATIC_LIBS) --end-group
	$(AR) -rs $@ neut.o

install-exec-hook:
	-(cd ${libdir}; $(LN_S) $(LIBNAME) $(shell echo $(LIBNAME) | cut -d . -f 1).a)
	-(mkdir -p ${bindir}/../share/neut; cp -r $(top_srcdir)/neutsmpl/Cards ${bindir}/../share/neut/ )
	-(mkdir -p ${bindir}/../share/neut/crsdat; cp -r $(top_srcdir)/crsdat/* ${bindir}/../share/neut/crsdat/ )

uninstall-hook:
	-(cd ${libdir}; rm -f $(shell echo $(LIBNAME) | cut -d . -f 1).a)

endif
