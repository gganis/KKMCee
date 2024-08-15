AUTOMAKE_OPTIONS = foreign subdir-objects
ACLOCAL_AMFLAGS = -I m4

# all_includes := $(all_includes)
# AM_CPPFLAGS = $(all_includes) -fPIC -Wall -Wno-write-strings # used while compiling
AM_CXXFLAGS = $(ROOTAUXCFLAGS) # used while compiling and linking

#F77=$(FC)
#AM_FFLAGS=$(AM_FCFLAGS) $(FCFLAGS)

lib_LTLIBRARIES = libKKpol.la

AM_FFLAGS = @MY_F77_FLAGS@

KKPOL_SRC = \
    anomwt2-Z.f

libKKpol_la_SOURCES = \
	$(KKPOL_SRC)
# libTauolaPhotos_la_DEPENDENCIES = HepEvt.h

libKKpol_la_LIBADD = $(FLIBS)
##	$(top_builddir)/jetset/libJetset.a
libKKpol_la_LDFLAGS = #-Wl,-z,defs

DISTCLEANFILES = $(BUILT_SOURCES)
