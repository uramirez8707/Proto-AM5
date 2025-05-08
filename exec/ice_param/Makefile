# Makefile created by mkmf 2022.01.00

CPPDEFS = -DINTERNAL_FILE_NML
OTHERFLAGS = -I$(BUILDROOT)/FMS 
include $(MK_TEMPLATE)


.DEFAULT:
	-echo $@ does not exist.
all: libice_param.a
ice_albedo.o: $(SRCROOT)/ice_param/ice_albedo.F90
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)/ice_param/ice_albedo.F90
ocean_albedo.o: $(SRCROOT)/ice_param/ocean_albedo.F90
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)/ice_param/ocean_albedo.F90
ocean_rough.o: $(SRCROOT)/ice_param/ocean_rough.F90
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)/ice_param/ocean_rough.F90
./ice_albedo.F90: $(SRCROOT)/ice_param/ice_albedo.F90
	cp $(SRCROOT)/ice_param/ice_albedo.F90 .
./ocean_rough.F90: $(SRCROOT)/ice_param/ocean_rough.F90
	cp $(SRCROOT)/ice_param/ocean_rough.F90 .
./ocean_albedo.F90: $(SRCROOT)/ice_param/ocean_albedo.F90
	cp $(SRCROOT)/ice_param/ocean_albedo.F90 .
SRC = $(SRCROOT)/ice_param/ocean_albedo.F90 $(SRCROOT)/ice_param/ocean_rough.F90 $(SRCROOT)/ice_param/ice_albedo.F90
OBJ = ocean_albedo.o ocean_rough.o ice_albedo.o
OFF = $(SRCROOT)/ice_param/ice_albedo.F90 $(SRCROOT)/ice_param/ocean_rough.F90 $(SRCROOT)/ice_param/ocean_albedo.F90
clean: neat
	-rm -f .libice_param.a.cppdefs $(OBJ) libice_param.a

neat:
	-rm -f $(TMPFILES)

localize: $(OFF)
	cp $(OFF) .

TAGS: $(SRC)
	etags $(SRC)

tags: $(SRC)
	ctags $(SRC)

libice_param.a: $(OBJ)
	$(AR) $(ARFLAGS) libice_param.a $(OBJ)

