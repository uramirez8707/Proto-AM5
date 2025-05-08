# Makefile created by mkmf 2022.01.00

CPPDEFS = -DINTERNAL_FILE_NML -DSPMD
OTHERFLAGS = -I$(BUILDROOT)/FMS -I$(BUILDROOT)/am5_phys -I$(BUILDROOT)/GFDL_atmos_cubed_sphere 
include $(MK_TEMPLATE)


.DEFAULT:
	-echo $@ does not exist.
all: libatmos_drivers.a
atmos_model.o: $(SRCROOT)/atmos_drivers/coupled/atmos_model.F90 $(SRCROOT)FMS/include/file_version.h
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)FMS/include	$(SRCROOT)/atmos_drivers/coupled/atmos_model.F90
./atmos_model.F90: $(SRCROOT)/atmos_drivers/coupled/atmos_model.F90
	cp $(SRCROOT)/atmos_drivers/coupled/atmos_model.F90 .
./file_version.h: $(SRCROOT)FMS/include/file_version.h
	cp $(SRCROOT)FMS/include/file_version.h .
SRC = $(SRCROOT)/atmos_drivers/coupled/atmos_model.F90 $(SRCROOT)FMS/include/file_version.h
OBJ = atmos_model.o
OFF = $(SRCROOT)/atmos_drivers/coupled/atmos_model.F90 $(SRCROOT)FMS/include/file_version.h
clean: neat
	-rm -f .libatmos_drivers.a.cppdefs $(OBJ) libatmos_drivers.a

neat:
	-rm -f $(TMPFILES)

localize: $(OFF)
	cp $(OFF) .

TAGS: $(SRC)
	etags $(SRC)

tags: $(SRC)
	ctags $(SRC)

libatmos_drivers.a: $(OBJ)
	$(AR) $(ARFLAGS) libatmos_drivers.a $(OBJ)

