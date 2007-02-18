###
bin_PROGRAMS = pfshdrcreate
bin_SCRIPTS  = jpeg2hdrgen pfsinhdrgen

pfshdrcreate_SOURCES = pfshdrcreate.cpp \
			  responses.cpp responses.h \
			  robertson02.cpp robertson02.h \
			  icip06.cpp icip06.h

man_MANS =  jpeg2hdrgen.1 pfsinhdrgen.1 pfshdrcreate.1

LIBS = $(PFS_LIBS) -lm
INCLUDES = $(PFS_CFLAGS)
EXTRA_DIST = $(man_MANS) $(bin_SCRIPTS)
