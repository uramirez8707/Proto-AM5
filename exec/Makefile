# Makefile for experiment_name
SRCROOT = ${BASEDIR}/src/
BUILDROOT = ${BASEDIR}/exec/
MK_TEMPLATE = ${BASEDIR}/exec/intel-classic.mk
include $(MK_TEMPLATE)

experiment_name.x:  FMScoupler/libFMScoupler.a SIS2/libSIS2.a atmos_drivers/libatmos_drivers.a GFDL_atmos_cubed_sphere/libGFDL_atmos_cubed_sphere.a am5_phys/libam5_phys.a land_lad2/libland_lad2.a MOM6/libMOM6.a ice_param/libice_param.a icebergs/libicebergs.a FMS/libFMS.a 
	$(LD) $^ $(LDFLAGS) -o $@ $(STATIC_LIBS)
FMScoupler/libFMScoupler.a:  FMS/libFMS.a GFDL_atmos_cubed_sphere/libGFDL_atmos_cubed_sphere.a atmos_drivers/libatmos_drivers.a am5_phys/libam5_phys.a land_lad2/libland_lad2.a SIS2/libSIS2.a MOM6/libMOM6.a ice_param/libice_param.a icebergs/libicebergs.a  FORCE
	$(MAKE) SRCROOT=$(SRCROOT) BUILDROOT=$(BUILDROOT) MK_TEMPLATE=$(MK_TEMPLATE) --directory=FMScoupler $(@F)
SIS2/libSIS2.a:  FMS/libFMS.a MOM6/libMOM6.a icebergs/libicebergs.a ice_param/libice_param.a  FORCE
	$(MAKE) SRCROOT=$(SRCROOT) BUILDROOT=$(BUILDROOT) MK_TEMPLATE=$(MK_TEMPLATE) --directory=SIS2 $(@F)
atmos_drivers/libatmos_drivers.a:  FMS/libFMS.a am5_phys/libam5_phys.a GFDL_atmos_cubed_sphere/libGFDL_atmos_cubed_sphere.a  FORCE
	$(MAKE) SRCROOT=$(SRCROOT) BUILDROOT=$(BUILDROOT) MK_TEMPLATE=$(MK_TEMPLATE) --directory=atmos_drivers $(@F)
GFDL_atmos_cubed_sphere/libGFDL_atmos_cubed_sphere.a:  FMS/libFMS.a am5_phys/libam5_phys.a  FORCE
	$(MAKE) SRCROOT=$(SRCROOT) BUILDROOT=$(BUILDROOT) MK_TEMPLATE=$(MK_TEMPLATE) --directory=GFDL_atmos_cubed_sphere $(@F)
am5_phys/libam5_phys.a:  FMS/libFMS.a  FORCE
	$(MAKE) SRCROOT=$(SRCROOT) BUILDROOT=$(BUILDROOT) MK_TEMPLATE=$(MK_TEMPLATE) --directory=am5_phys $(@F)
land_lad2/libland_lad2.a:  FMS/libFMS.a  FORCE
	$(MAKE) SRCROOT=$(SRCROOT) BUILDROOT=$(BUILDROOT) MK_TEMPLATE=$(MK_TEMPLATE) --directory=land_lad2 $(@F)
MOM6/libMOM6.a:  FMS/libFMS.a  FORCE
	$(MAKE) SRCROOT=$(SRCROOT) BUILDROOT=$(BUILDROOT) MK_TEMPLATE=$(MK_TEMPLATE) OPENMP="" --directory=MOM6 $(@F)
ice_param/libice_param.a:  FMS/libFMS.a  FORCE
	$(MAKE) SRCROOT=$(SRCROOT) BUILDROOT=$(BUILDROOT) MK_TEMPLATE=$(MK_TEMPLATE) --directory=ice_param $(@F)
icebergs/libicebergs.a:  FMS/libFMS.a  FORCE
	$(MAKE) SRCROOT=$(SRCROOT) BUILDROOT=$(BUILDROOT) MK_TEMPLATE=$(MK_TEMPLATE) --directory=icebergs $(@F)
FMS/libFMS.a:   FORCE
	$(MAKE) SRCROOT=$(SRCROOT) BUILDROOT=$(BUILDROOT) MK_TEMPLATE=$(MK_TEMPLATE) --directory=FMS $(@F)
FORCE:

clean:
	$(MAKE) --directory=FMS clean
	$(MAKE) --directory=am5_phys clean
	$(MAKE) --directory=GFDL_atmos_cubed_sphere clean
	$(MAKE) --directory=atmos_drivers clean
	$(MAKE) --directory=land_lad2 clean
	$(MAKE) --directory=MOM6 clean
	$(MAKE) --directory=ice_param clean
	$(MAKE) --directory=SIS2 clean
	$(MAKE) --directory=icebergs clean
	$(MAKE) --directory=FMScoupler clean
localize:
	$(MAKE) -f $(BUILDROOT)FMS localize
	$(MAKE) -f $(BUILDROOT)am5_phys localize
	$(MAKE) -f $(BUILDROOT)GFDL_atmos_cubed_sphere localize
	$(MAKE) -f $(BUILDROOT)atmos_drivers localize
	$(MAKE) -f $(BUILDROOT)land_lad2 localize
	$(MAKE) -f $(BUILDROOT)MOM6 localize
	$(MAKE) -f $(BUILDROOT)ice_param localize
	$(MAKE) -f $(BUILDROOT)SIS2 localize
	$(MAKE) -f $(BUILDROOT)icebergs localize
	$(MAKE) -f $(BUILDROOT)FMScoupler localize
distclean:
	$(RM) -r FMS
	$(RM) -r am5_phys
	$(RM) -r GFDL_atmos_cubed_sphere
	$(RM) -r atmos_drivers
	$(RM) -r land_lad2
	$(RM) -r MOM6
	$(RM) -r ice_param
	$(RM) -r SIS2
	$(RM) -r icebergs
	$(RM) -r FMScoupler
	$(RM) -r experiment_name
	$(RM) -r Makefile 
