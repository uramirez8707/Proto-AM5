# Makefile created by mkmf 2022.01.00

OTHERFLAGS = -I$(BUILDROOT)/FMS -I$(BUILDROOT)/GFDL_atmos_cubed_sphere -I$(BUILDROOT)/atmos_drivers -I$(BUILDROOT)/am5_phys -I$(BUILDROOT)/land_lad2 -I$(BUILDROOT)/SIS2 -I$(BUILDROOT)/MOM6 -I$(BUILDROOT)/ice_param -I$(BUILDROOT)/icebergs 
include $(MK_TEMPLATE)

.DEFAULT:
	-echo $@ does not exist.
all: libFMScoupler.a
atm_land_ice_flux_exchange.o: $(SRCROOT)/FMScoupler/full/atm_land_ice_flux_exchange.F90 surface_flux.o atmos_ocean_fluxes_calc.o atmos_ocean_dep_fluxes_calc.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)/FMScoupler/full/atm_land_ice_flux_exchange.F90
atmos_ocean_dep_fluxes_calc.o: $(SRCROOT)/FMScoupler/full/atmos_ocean_dep_fluxes_calc.F90
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)/FMScoupler/full/atmos_ocean_dep_fluxes_calc.F90
atmos_ocean_fluxes_calc.o: $(SRCROOT)/FMScoupler/full/atmos_ocean_fluxes_calc.F90
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)/FMScoupler/full/atmos_ocean_fluxes_calc.F90
coupler_main.o: $(SRCROOT)/FMScoupler/full/coupler_main.F90 flux_exchange.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)/FMScoupler/full/coupler_main.F90
flux_exchange.o: $(SRCROOT)/FMScoupler/full/flux_exchange.F90 atm_land_ice_flux_exchange.o land_ice_flux_exchange.o ice_ocean_flux_exchange.o atmos_ocean_fluxes_calc.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)/FMScoupler/full/flux_exchange.F90
ice_ocean_flux_exchange.o: $(SRCROOT)/FMScoupler/full/ice_ocean_flux_exchange.F90
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)/FMScoupler/full/ice_ocean_flux_exchange.F90
land_ice_flux_exchange.o: $(SRCROOT)/FMScoupler/full/land_ice_flux_exchange.F90
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)/FMScoupler/full/land_ice_flux_exchange.F90
surface_flux.o: $(SRCROOT)/FMScoupler/shared/surface_flux.F90
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)/FMScoupler/shared/surface_flux.F90
./coupler_main.F90: $(SRCROOT)/FMScoupler/full/coupler_main.F90
	cp $(SRCROOT)/FMScoupler/full/coupler_main.F90 .
./surface_flux.F90: $(SRCROOT)/FMScoupler/shared/surface_flux.F90
	cp $(SRCROOT)/FMScoupler/shared/surface_flux.F90 .
./flux_exchange.F90: $(SRCROOT)/FMScoupler/full/flux_exchange.F90
	cp $(SRCROOT)/FMScoupler/full/flux_exchange.F90 .
./land_ice_flux_exchange.F90: $(SRCROOT)/FMScoupler/full/land_ice_flux_exchange.F90
	cp $(SRCROOT)/FMScoupler/full/land_ice_flux_exchange.F90 .
./atm_land_ice_flux_exchange.F90: $(SRCROOT)/FMScoupler/full/atm_land_ice_flux_exchange.F90
	cp $(SRCROOT)/FMScoupler/full/atm_land_ice_flux_exchange.F90 .
./atmos_ocean_dep_fluxes_calc.F90: $(SRCROOT)/FMScoupler/full/atmos_ocean_dep_fluxes_calc.F90
	cp $(SRCROOT)/FMScoupler/full/atmos_ocean_dep_fluxes_calc.F90 .
./atmos_ocean_fluxes_calc.F90: $(SRCROOT)/FMScoupler/full/atmos_ocean_fluxes_calc.F90
	cp $(SRCROOT)/FMScoupler/full/atmos_ocean_fluxes_calc.F90 .
./ice_ocean_flux_exchange.F90: $(SRCROOT)/FMScoupler/full/ice_ocean_flux_exchange.F90
	cp $(SRCROOT)/FMScoupler/full/ice_ocean_flux_exchange.F90 .
SRC = $(SRCROOT)/FMScoupler/full/flux_exchange.F90 $(SRCROOT)/FMScoupler/full/ice_ocean_flux_exchange.F90 $(SRCROOT)/FMScoupler/full/atm_land_ice_flux_exchange.F90 $(SRCROOT)/FMScoupler/shared/surface_flux.F90 $(SRCROOT)/FMScoupler/full/atmos_ocean_dep_fluxes_calc.F90 $(SRCROOT)/FMScoupler/full/coupler_main.F90 $(SRCROOT)/FMScoupler/full/land_ice_flux_exchange.F90 $(SRCROOT)/FMScoupler/full/atmos_ocean_fluxes_calc.F90
OBJ = flux_exchange.o ice_ocean_flux_exchange.o atm_land_ice_flux_exchange.o surface_flux.o atmos_ocean_dep_fluxes_calc.o coupler_main.o land_ice_flux_exchange.o atmos_ocean_fluxes_calc.o
OFF = $(SRCROOT)/FMScoupler/full/coupler_main.F90 $(SRCROOT)/FMScoupler/shared/surface_flux.F90 $(SRCROOT)/FMScoupler/full/flux_exchange.F90 $(SRCROOT)/FMScoupler/full/land_ice_flux_exchange.F90 $(SRCROOT)/FMScoupler/full/atm_land_ice_flux_exchange.F90 $(SRCROOT)/FMScoupler/full/atmos_ocean_dep_fluxes_calc.F90 $(SRCROOT)/FMScoupler/full/atmos_ocean_fluxes_calc.F90 $(SRCROOT)/FMScoupler/full/ice_ocean_flux_exchange.F90
clean: neat
	-rm -f .libFMScoupler.a.cppdefs $(OBJ) libFMScoupler.a

neat:
	-rm -f $(TMPFILES)

localize: $(OFF)
	cp $(OFF) .

TAGS: $(SRC)
	etags $(SRC)

tags: $(SRC)
	ctags $(SRC)

libFMScoupler.a: $(OBJ)
	$(AR) $(ARFLAGS) libFMScoupler.a $(OBJ)

