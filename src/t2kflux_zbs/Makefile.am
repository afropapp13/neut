include $(top_srcdir)/build.common

AM_FFLAGS = -DT2KFLUX_ZBSVER=@NEUT_VERSION_NUM@ -DFLUX_11B $(DEFAULT_FFLAGS)

LIBNAME = libt2kflux_zbs.@NEUT_VERSION_NUM@.a

noinst_LIBRARIES = $(LIBNAME)

libt2kflux_zbs_@NEUT_VERSION_NUM@_a_SOURCES = necardbm.F 
# 						necpjndfxcm.F nemkjndfx.F nerdnufxv.F rnbmom2802.F \
# 						neopnufxv.F nerd1nufxv.F \
# 						neopskfxv.F rnpos.F reset_signal_.c

# Unused sources that live in the directory
UNUSED_srcs =
