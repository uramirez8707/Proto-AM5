# Makefile created by mkmf 2022.01.00

CPPDEFS = -DINTERNAL_FILE_NML -DUSE_FMS2_IO
OTHERFLAGS = -I$(BUILDROOT)/FMS 
include $(MK_TEMPLATE)


.DEFAULT:
	-echo $@ does not exist.
all: libicebergs.a
icebergs.o: $(SRCROOT)/icebergs/src/icebergs.F90 icebergs_framework.o icebergs_io.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)/icebergs/src/icebergs.F90
icebergs_fms2io.o: $(SRCROOT)/icebergs/src/icebergs_fms2io.F90 icebergs_framework.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)/icebergs/src/icebergs_fms2io.F90
icebergs_fmsio.o: $(SRCROOT)/icebergs/src/icebergs_fmsio.F90 icebergs_framework.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)/icebergs/src/icebergs_fmsio.F90
icebergs_framework.o: $(SRCROOT)/icebergs/src/icebergs_framework.F90
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)/icebergs/src/icebergs_framework.F90
icebergs_io.o: $(SRCROOT)/icebergs/src/icebergs_io.F90 icebergs_fmsio.o icebergs_fms2io.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)/icebergs/src/icebergs_io.F90
./icebergs_io.F90: $(SRCROOT)/icebergs/src/icebergs_io.F90
	cp $(SRCROOT)/icebergs/src/icebergs_io.F90 .
./icebergs_framework.F90: $(SRCROOT)/icebergs/src/icebergs_framework.F90
	cp $(SRCROOT)/icebergs/src/icebergs_framework.F90 .
./icebergs_fmsio.F90: $(SRCROOT)/icebergs/src/icebergs_fmsio.F90
	cp $(SRCROOT)/icebergs/src/icebergs_fmsio.F90 .
./icebergs_fms2io.F90: $(SRCROOT)/icebergs/src/icebergs_fms2io.F90
	cp $(SRCROOT)/icebergs/src/icebergs_fms2io.F90 .
./icebergs.F90: $(SRCROOT)/icebergs/src/icebergs.F90
	cp $(SRCROOT)/icebergs/src/icebergs.F90 .
SRC = $(SRCROOT)/icebergs/src/icebergs_fms2io.F90 $(SRCROOT)/icebergs/src/icebergs_fmsio.F90 $(SRCROOT)/icebergs/src/icebergs_io.F90 $(SRCROOT)/icebergs/src/icebergs_framework.F90 $(SRCROOT)/icebergs/src/icebergs.F90
OBJ = icebergs_fms2io.o icebergs_fmsio.o icebergs_io.o icebergs_framework.o icebergs.o
OFF = $(SRCROOT)/icebergs/src/icebergs_io.F90 $(SRCROOT)/icebergs/src/icebergs_framework.F90 $(SRCROOT)/icebergs/src/icebergs_fmsio.F90 $(SRCROOT)/icebergs/src/icebergs_fms2io.F90 $(SRCROOT)/icebergs/src/icebergs.F90
clean: neat
	-rm -f .libicebergs.a.cppdefs $(OBJ) libicebergs.a

neat:
	-rm -f $(TMPFILES)

localize: $(OFF)
	cp $(OFF) .

TAGS: $(SRC)
	etags $(SRC)

tags: $(SRC)
	ctags $(SRC)

libicebergs.a: $(OBJ)
	$(AR) $(ARFLAGS) libicebergs.a $(OBJ)

