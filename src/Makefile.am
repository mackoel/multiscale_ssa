## Process this file with automake to produce Makefile.in

## Created by Anjuta

AM_CPPFLAGS = \
	-DPACKAGE_LOCALE_DIR=\""$(localedir)"\" \
	-DPACKAGE_SRC_DIR=\""$(srcdir)"\" \
	-DPACKAGE_DATA_DIR=\""$(pkgdatadir)"\" \
	$(MULTISCALE_SSA_CFLAGS)

AM_CFLAGS =\
	 -Wall\
	 -g

bin_PROGRAMS = multiscale_ssa

multiscale_ssa_SOURCES = \
	main.c

multiscale_ssa_CFLAGS = -fopenmp

multiscale_ssa_LDFLAGS = -lm

multiscale_ssa_LDADD = $(MULTISCALE_SSA_LIBS)

