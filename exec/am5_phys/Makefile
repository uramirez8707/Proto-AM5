# Makefile created by mkmf 2022.01.00

CPPDEFS = -DINTERNAL_FILE_NML
OTHERFLAGS = -I$(BUILDROOT)/FMS 
include $(MK_TEMPLATE)


.DEFAULT:
	-echo $@ does not exist.
all: libam5_phys.a
MISR_simulator.o: $(SRCROOT)/am5_phys/atmos_param/cosp/MISR_simulator/MISR_simulator.F $(SRCROOT)/am5_phys/atmos_param/cosp/cosp_defs.H
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)/am5_phys/atmos_param/cosp	$(SRCROOT)/am5_phys/atmos_param/cosp/MISR_simulator/MISR_simulator.F
aer_ccn_act.o: $(SRCROOT)/am5_phys/atmos_shared/tracer_driver/aer_ccn_act/aer_ccn_act.F90 aer_ccn_act_k.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)/am5_phys/atmos_shared/tracer_driver/aer_ccn_act/aer_ccn_act.F90
aer_ccn_act_k.o: $(SRCROOT)/am5_phys/atmos_shared/tracer_driver/aer_ccn_act/aer_ccn_act_k.F90
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)/am5_phys/atmos_shared/tracer_driver/aer_ccn_act/aer_ccn_act_k.F90
aer_in_act.o: $(SRCROOT)/am5_phys/atmos_shared/tracer_driver/aer_ccn_act/aer_in_act.F90
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)/am5_phys/atmos_shared/tracer_driver/aer_ccn_act/aer_in_act.F90
aerosol.o: $(SRCROOT)/am5_phys/atmos_shared/aerosol/aerosol.F90 aerosol_types.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)/am5_phys/atmos_shared/aerosol/aerosol.F90
aerosol_cloud.o: $(SRCROOT)/am5_phys/atmos_param/aerosol_cloud/aerosol_cloud.F90 aerosol_types.o polysvp.o aerosol_params.o aer_ccn_act.o ice_nucl.o lscloud_types.o lscloud_debug.o moist_processes_utils.o physics_radiation_exch.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)/am5_phys/atmos_param/aerosol_cloud/aerosol_cloud.F90
aerosol_params.o: $(SRCROOT)/am5_phys/atmos_shared/tracer_driver/aer_ccn_act/aerosol_params.F90
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)/am5_phys/atmos_shared/tracer_driver/aer_ccn_act/aerosol_params.F90
aerosol_thermodynamics.o: $(SRCROOT)/am5_phys/atmos_shared/tracer_driver/tropchem/aerosol_thermodynamics.F90 tropchem_types_mod.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)/am5_phys/atmos_shared/tracer_driver/tropchem/aerosol_thermodynamics.F90
aerosol_types.o: $(SRCROOT)/am5_phys/atmos_shared/aerosol/aerosol_types.F90
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)/am5_phys/atmos_shared/aerosol/aerosol_types.F90
aerosolrad_driver.o: $(SRCROOT)/am5_phys/atmos_param/radiation/aerosols/aerosolrad_driver.F90 aerosol_types.o aerosol.o aerosolrad_types.o aerosolrad_package.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)/am5_phys/atmos_param/radiation/aerosols/aerosolrad_driver.F90
aerosolrad_package.o: $(SRCROOT)/am5_phys/atmos_param/radiation/aerosols/aerosolrad_package.F90 aerosol_types.o aerosolrad_types.o sealw99.o esfsw_driver.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)/am5_phys/atmos_param/radiation/aerosols/aerosolrad_package.F90
aerosolrad_types.o: $(SRCROOT)/am5_phys/atmos_param/radiation/aerosols/aerosolrad_types.F90
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)/am5_phys/atmos_param/radiation/aerosols/aerosolrad_types.F90
array_lib.o: $(SRCROOT)/am5_phys/atmos_param/cosp/quickbeam/array_lib.F90 $(SRCROOT)/am5_phys/atmos_param/cosp/cosp_defs.H mrgrnk.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)/am5_phys/atmos_param/cosp	$(SRCROOT)/am5_phys/atmos_param/cosp/quickbeam/array_lib.F90
atmos_age_tracer.o: $(SRCROOT)/am5_phys/atmos_shared/tracer_driver/atmos_age_tracer.F90
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)/am5_phys/atmos_shared/tracer_driver/atmos_age_tracer.F90
atmos_carbon_aerosol.o: $(SRCROOT)/am5_phys/atmos_shared/tracer_driver/atmos_carbon_aerosol.F90 atmos_cmip_diag.o atmos_tracer_utilities.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)/am5_phys/atmos_shared/tracer_driver/atmos_carbon_aerosol.F90
atmos_ch3i.o: $(SRCROOT)/am5_phys/atmos_shared/tracer_driver/atmos_ch3i.F90
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)/am5_phys/atmos_shared/tracer_driver/atmos_ch3i.F90
atmos_ch4.o: $(SRCROOT)/am5_phys/atmos_shared/tracer_driver/atmos_ch4.F90
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)/am5_phys/atmos_shared/tracer_driver/atmos_ch4.F90
atmos_cmip_diag.o: $(SRCROOT)/am5_phys/atmos_shared/atmos_cmip_diag/atmos_cmip_diag.F90 $(SRCROOT)/am5_phys/atmos_shared/atmos_cmip_diag/atmos_cmip_interp.inc
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)/am5_phys/atmos_shared/atmos_cmip_diag	$(SRCROOT)/am5_phys/atmos_shared/atmos_cmip_diag/atmos_cmip_diag.F90
atmos_co2.o: $(SRCROOT)/am5_phys/atmos_shared/tracer_driver/atmos_co2.F90
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)/am5_phys/atmos_shared/tracer_driver/atmos_co2.F90
atmos_convection_tracer.o: $(SRCROOT)/am5_phys/atmos_shared/tracer_driver/atmos_convection_tracer.F90 atmos_tracer_utilities.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)/am5_phys/atmos_shared/tracer_driver/atmos_convection_tracer.F90
atmos_dust.o: $(SRCROOT)/am5_phys/atmos_shared/tracer_driver/atmos_dust.F90 atmos_cmip_diag.o atmos_tracer_utilities.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)/am5_phys/atmos_shared/tracer_driver/atmos_dust.F90
atmos_global_diag.o: $(SRCROOT)/am5_phys/atmos_shared/atmos_cmip_diag/atmos_global_diag.F90
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)/am5_phys/atmos_shared/atmos_cmip_diag/atmos_global_diag.F90
atmos_lib.o: $(SRCROOT)/am5_phys/atmos_param/cosp/quickbeam/atmos_lib.F90 $(SRCROOT)/am5_phys/atmos_param/cosp/cosp_defs.H
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)/am5_phys/atmos_param/cosp	$(SRCROOT)/am5_phys/atmos_param/cosp/quickbeam/atmos_lib.F90
atmos_nh3_tag.o: $(SRCROOT)/am5_phys/atmos_shared/tracer_driver/atmos_nh3_tag.F90 atmos_cmip_diag.o atmos_tracer_utilities.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)/am5_phys/atmos_shared/tracer_driver/atmos_nh3_tag.F90
atmos_nudge.o: $(SRCROOT)/am5_phys/atmos_shared/atmos_nudge/atmos_nudge.F90
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)/am5_phys/atmos_shared/atmos_nudge/atmos_nudge.F90
atmos_ozone_tracer.o: $(SRCROOT)/am5_phys/atmos_shared/tracer_driver/atmos_ozone_tracer.F90
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)/am5_phys/atmos_shared/tracer_driver/atmos_ozone_tracer.F90
atmos_radon.o: $(SRCROOT)/am5_phys/atmos_shared/tracer_driver/atmos_radon.F90 atmos_tracer_utilities.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)/am5_phys/atmos_shared/tracer_driver/atmos_radon.F90
atmos_regional_tracer_driver.o: $(SRCROOT)/am5_phys/atmos_shared/tracer_driver/atmos_regional_tracer_driver.F90
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)/am5_phys/atmos_shared/tracer_driver/atmos_regional_tracer_driver.F90
atmos_sea_salt.o: $(SRCROOT)/am5_phys/atmos_shared/tracer_driver/atmos_sea_salt.F90 atmos_cmip_diag.o atmos_tracer_utilities.o atmos_carbon_aerosol.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)/am5_phys/atmos_shared/tracer_driver/atmos_sea_salt.F90
atmos_soa.o: $(SRCROOT)/am5_phys/atmos_shared/tracer_driver/atmos_soa.F90 atmos_cmip_diag.o xactive_bvoc.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)/am5_phys/atmos_shared/tracer_driver/atmos_soa.F90
atmos_sulfate.o: $(SRCROOT)/am5_phys/atmos_shared/tracer_driver/atmos_sulfate.F90 atmos_cmip_diag.o cloud_chem.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)/am5_phys/atmos_shared/tracer_driver/atmos_sulfate.F90
atmos_sulfur_hex.o: $(SRCROOT)/am5_phys/atmos_shared/tracer_driver/atmos_sulfur_hex.F90 atmos_tracer_utilities.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)/am5_phys/atmos_shared/tracer_driver/atmos_sulfur_hex.F90
atmos_tracer_driver.o: $(SRCROOT)/am5_phys/atmos_shared/tracer_driver/atmos_tracer_driver.F90 atmos_cmip_diag.o atmos_tracer_utilities.o atmos_radon.o atmos_carbon_aerosol.o atmos_convection_tracer.o atmos_sulfur_hex.o atmos_ch3i.o atmos_sea_salt.o atmos_dust.o atmos_sulfate.o atmos_soa.o atmos_nh3_tag.o tropchem_driver.o atmos_regional_tracer_driver.o strat_chem_driver.o atmos_age_tracer.o atmos_co2.o atmos_ch4.o atmos_tropopause.o xactive_bvoc.o atmos_ozone_tracer.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)/am5_phys/atmos_shared/tracer_driver/atmos_tracer_driver.F90
atmos_tracer_utilities.o: $(SRCROOT)/am5_phys/atmos_shared/tracer_driver/atmos_tracer_utilities.F90 atmos_cmip_diag.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)/am5_phys/atmos_shared/tracer_driver/atmos_tracer_utilities.F90
atmos_tropopause.o: $(SRCROOT)/am5_phys/atmos_shared/tracer_driver/atmos_tropopause.F90 atmos_cmip_diag.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)/am5_phys/atmos_shared/tracer_driver/atmos_tropopause.F90
betaDistribution.o: $(SRCROOT)/am5_phys/atmos_param/cloud_generator/betaDistribution.F90
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)/am5_phys/atmos_param/cloud_generator/betaDistribution.F90
calc_Re.o: $(SRCROOT)/am5_phys/atmos_param/cosp/quickbeam/calc_Re.F90 $(SRCROOT)/am5_phys/atmos_param/cosp/cosp_defs.H math_lib.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)/am5_phys/atmos_param/cosp	$(SRCROOT)/am5_phys/atmos_param/cosp/quickbeam/calc_Re.F90
cg_drag.o: $(SRCROOT)/am5_phys/atmos_param/cg_drag/cg_drag.F90
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)/am5_phys/atmos_param/cg_drag/cg_drag.F90
check_nan.o: $(SRCROOT)/am5_phys/atmos_param/lscloud_driver/check_nan.F90
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)/am5_phys/atmos_param/lscloud_driver/check_nan.F90
cloud_chem.o: $(SRCROOT)/am5_phys/atmos_shared/tracer_driver/cloud_chem.F90
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)/am5_phys/atmos_shared/tracer_driver/cloud_chem.F90
cloud_generator.o: $(SRCROOT)/am5_phys/atmos_param/cloud_generator/cloud_generator.F90 betaDistribution.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)/am5_phys/atmos_param/cloud_generator/cloud_generator.F90
cloud_rad.o: $(SRCROOT)/am5_phys/atmos_param/cloud_rad/cloud_rad.F90 gamma_mg.o lscloud_constants.o aerosol_types.o physics_radiation_exch.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)/am5_phys/atmos_param/cloud_rad/cloud_rad.F90
cloud_spec.o: $(SRCROOT)/am5_phys/atmos_param/radiation/clouds/cloud_spec.F90 get_random_number_stream.o aerosol_types.o physics_radiation_exch.o cloudrad_types.o strat_clouds_W.o uw_clouds_W.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)/am5_phys/atmos_param/radiation/clouds/cloud_spec.F90
cloudrad_diagnostics.o: $(SRCROOT)/am5_phys/atmos_param/radiation/clouds/cloudrad_diagnostics.F90 atmos_cmip_diag.o cloudrad_types.o microphys_rad.o isccp_clouds.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)/am5_phys/atmos_param/radiation/clouds/cloudrad_diagnostics.F90
cloudrad_driver.o: $(SRCROOT)/am5_phys/atmos_param/radiation/clouds/cloudrad_driver.F90 aerosol_types.o physics_radiation_exch.o cloudrad_types.o cloud_spec.o cloudrad_package.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)/am5_phys/atmos_param/radiation/clouds/cloudrad_driver.F90
cloudrad_package.o: $(SRCROOT)/am5_phys/atmos_param/radiation/clouds/cloudrad_package.F90 cloudrad_types.o cloudrad_diagnostics.o microphys_rad.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)/am5_phys/atmos_param/radiation/clouds/cloudrad_package.F90
cloudrad_types.o: $(SRCROOT)/am5_phys/atmos_param/radiation/clouds/cloudrad_types.F90
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)/am5_phys/atmos_param/radiation/clouds/cloudrad_types.F90
conv_closures.o: $(SRCROOT)/am5_phys/atmos_param/shallow_cu/conv_closures.F90 conv_utilities_k.o conv_plumes_k.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)/am5_phys/atmos_param/shallow_cu/conv_closures.F90
conv_plumes.o: $(SRCROOT)/am5_phys/atmos_param/shallow_cu/conv_plumes.F90
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)/am5_phys/atmos_param/shallow_cu/conv_plumes.F90
conv_plumes_k.o: $(SRCROOT)/am5_phys/atmos_param/shallow_cu/conv_plumes_k.F90 $(SRCROOT)FMS/include/fms_platform.h aer_ccn_act_k.o conv_utilities_k.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)FMS/include	$(SRCROOT)/am5_phys/atmos_param/shallow_cu/conv_plumes_k.F90
conv_utilities.o: $(SRCROOT)/am5_phys/atmos_param/shallow_cu/conv_utilities.F90 conv_utilities_k.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)/am5_phys/atmos_param/shallow_cu/conv_utilities.F90
conv_utilities_k.o: $(SRCROOT)/am5_phys/atmos_param/shallow_cu/conv_utilities_k.F90 $(SRCROOT)FMS/include/fms_platform.h
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)FMS/include	$(SRCROOT)/am5_phys/atmos_param/shallow_cu/conv_utilities_k.F90
convection_driver.o: $(SRCROOT)/am5_phys/atmos_param/convection_driver/convection_driver.F90 atmos_global_diag.o atmos_cmip_diag.o physics_types.o vert_diff_driver.o physics_radiation_exch.o uw_conv.o detr_ice_num.o cu_mo_trans.o mo_hook.o aerosol_types.o moist_processes_utils.o convection_utilities.o atmos_tracer_utilities.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)/am5_phys/atmos_param/convection_driver/convection_driver.F90
convection_utilities.o: $(SRCROOT)/am5_phys/atmos_param/convection_driver/convection_utilities.F90
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)/am5_phys/atmos_param/convection_driver/convection_utilities.F90
cosp.o: $(SRCROOT)/am5_phys/atmos_param/cosp/cosp.F90 $(SRCROOT)/am5_phys/atmos_param/cosp/cosp_defs.H cosp_types.o cosp_simulator.o cosp_modis_simulator.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)/am5_phys/atmos_param/cosp	$(SRCROOT)/am5_phys/atmos_param/cosp/cosp.F90
cosp_constants.o: $(SRCROOT)/am5_phys/atmos_param/cosp/cosp_constants.F90 $(SRCROOT)/am5_phys/atmos_param/cosp/cosp_defs.H
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)/am5_phys/atmos_param/cosp	$(SRCROOT)/am5_phys/atmos_param/cosp/cosp_constants.F90
cosp_diagnostics.o: $(SRCROOT)/am5_phys/atmos_param/cosp/cosp_diagnostics.F90 $(SRCROOT)/am5_phys/atmos_param/cosp/cosp_defs.H cosp_types.o cosp_io.o cosp_constants.o lmd_ipsl_stats.o cosp_modis_simulator.o modis_simulator.o cosp_utils.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)/am5_phys/atmos_param/cosp	$(SRCROOT)/am5_phys/atmos_param/cosp/cosp_diagnostics.F90
cosp_driver.o: $(SRCROOT)/am5_phys/atmos_param/cosp/cosp_driver.F90 $(SRCROOT)/am5_phys/atmos_param/cosp/cosp_defs.H physics_radiation_exch.o physics_types.o cosp_types.o cosp.o cosp_io.o cosp_constants.o cosp_modis_simulator.o cosp_diagnostics.o cosp_utils.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)/am5_phys/atmos_param/cosp	$(SRCROOT)/am5_phys/atmos_param/cosp/cosp_driver.F90
cosp_io.o: $(SRCROOT)/am5_phys/atmos_param/cosp/cosp_io.F90 $(SRCROOT)/am5_phys/atmos_param/cosp/cosp_defs.H cosp_constants.o cosp_types.o cosp_modis_simulator.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)/am5_phys/atmos_param/cosp	$(SRCROOT)/am5_phys/atmos_param/cosp/cosp_io.F90
cosp_isccp_simulator.o: $(SRCROOT)/am5_phys/atmos_param/cosp/cosp_isccp_simulator.F90 $(SRCROOT)/am5_phys/atmos_param/cosp/cosp_defs.H cosp_constants.o cosp_types.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)/am5_phys/atmos_param/cosp	$(SRCROOT)/am5_phys/atmos_param/cosp/cosp_isccp_simulator.F90
cosp_lidar.o: $(SRCROOT)/am5_phys/atmos_param/cosp/cosp_lidar.F90 $(SRCROOT)/am5_phys/atmos_param/cosp/cosp_defs.H cosp_constants.o cosp_types.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)/am5_phys/atmos_param/cosp	$(SRCROOT)/am5_phys/atmos_param/cosp/cosp_lidar.F90
cosp_misr_simulator.o: $(SRCROOT)/am5_phys/atmos_param/cosp/cosp_misr_simulator.F90 $(SRCROOT)/am5_phys/atmos_param/cosp/cosp_defs.H cosp_constants.o cosp_types.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)/am5_phys/atmos_param/cosp	$(SRCROOT)/am5_phys/atmos_param/cosp/cosp_misr_simulator.F90
cosp_modis_simulator.o: $(SRCROOT)/am5_phys/atmos_param/cosp/cosp_modis_simulator.F90 $(SRCROOT)/am5_phys/atmos_param/cosp/cosp_defs.H cosp_constants.o cosp_types.o modis_simulator.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)/am5_phys/atmos_param/cosp	$(SRCROOT)/am5_phys/atmos_param/cosp/cosp_modis_simulator.F90
cosp_radar.o: $(SRCROOT)/am5_phys/atmos_param/cosp/llnl/cosp_radar.F90 $(SRCROOT)/am5_phys/atmos_param/cosp/cosp_defs.H cosp_constants.o cosp_types.o cosp_utils.o radar_simulator_types.o array_lib.o atmos_lib.o format_input.o mrgrnk.o math_lib.o optics_lib.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)/am5_phys/atmos_param/cosp	$(SRCROOT)/am5_phys/atmos_param/cosp/llnl/cosp_radar.F90
cosp_rttov_simulator.o: $(SRCROOT)/am5_phys/atmos_param/cosp/cosp_rttov_simulator.F90 $(SRCROOT)/am5_phys/atmos_param/cosp/cosp_defs.H cosp_constants.o cosp_types.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)/am5_phys/atmos_param/cosp	$(SRCROOT)/am5_phys/atmos_param/cosp/cosp_rttov_simulator.F90
cosp_simulator.o: $(SRCROOT)/am5_phys/atmos_param/cosp/cosp_simulator.F90 $(SRCROOT)/am5_phys/atmos_param/cosp/cosp_defs.H cosp_constants.o cosp_types.o cosp_radar.o cosp_lidar.o cosp_isccp_simulator.o cosp_modis_simulator.o cosp_misr_simulator.o cosp_rttov_simulator.o cosp_stats.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)/am5_phys/atmos_param/cosp	$(SRCROOT)/am5_phys/atmos_param/cosp/cosp_simulator.F90
cosp_stats.o: $(SRCROOT)/am5_phys/atmos_param/cosp/cosp_stats.F90 $(SRCROOT)/am5_phys/atmos_param/cosp/cosp_defs.H cosp_constants.o cosp_types.o llnl_stats.o lmd_ipsl_stats.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)/am5_phys/atmos_param/cosp	$(SRCROOT)/am5_phys/atmos_param/cosp/cosp_stats.F90
cosp_types.o: $(SRCROOT)/am5_phys/atmos_param/cosp/cosp_types.F90 $(SRCROOT)/am5_phys/atmos_param/cosp/cosp_defs.H cosp_constants.o cosp_utils.o radar_simulator_types.o scale_LUTs_io.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)/am5_phys/atmos_param/cosp	$(SRCROOT)/am5_phys/atmos_param/cosp/cosp_types.F90
cosp_utils.o: $(SRCROOT)/am5_phys/atmos_param/cosp/cosp_utils.F90 $(SRCROOT)/am5_phys/atmos_param/cosp/cosp_defs.H cosp_constants.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)/am5_phys/atmos_param/cosp	$(SRCROOT)/am5_phys/atmos_param/cosp/cosp_utils.F90
cu_mo_trans.o: $(SRCROOT)/am5_phys/atmos_param/cu_mo_trans/cu_mo_trans.F90 convection_utilities.o moist_processes_utils.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)/am5_phys/atmos_param/cu_mo_trans/cu_mo_trans.F90
damping_driver.o: $(SRCROOT)/am5_phys/atmos_param/damping_driver/damping_driver.F90 mg_drag.o cg_drag.o topo_drag.o gw_beres.o atmos_cmip_diag.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)/am5_phys/atmos_param/damping_driver/damping_driver.F90
deep_conv.o: $(SRCROOT)/am5_phys/atmos_param/shallow_cu/deep_conv.F90 conv_utilities.o conv_utilities_k.o conv_plumes_k.o conv_closures.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)/am5_phys/atmos_param/shallow_cu/deep_conv.F90
detr_ice_num.o: $(SRCROOT)/am5_phys/atmos_param/convection_driver/detr_ice_num.F90 lscloud_constants.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)/am5_phys/atmos_param/convection_driver/detr_ice_num.F90
dsd.o: $(SRCROOT)/am5_phys/atmos_param/cosp/quickbeam/dsd.F90 $(SRCROOT)/am5_phys/atmos_param/cosp/cosp_defs.H array_lib.o math_lib.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)/am5_phys/atmos_param/cosp	$(SRCROOT)/am5_phys/atmos_param/cosp/quickbeam/dsd.F90
entrain.o: $(SRCROOT)/am5_phys/atmos_param/entrain/entrain.F90
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)/am5_phys/atmos_param/entrain/entrain.F90
esfsw_bands.o: $(SRCROOT)/am5_phys/atmos_param/radiation/shortwave/esfsw_bands.F90 esfsw_utilities.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)/am5_phys/atmos_param/radiation/shortwave/esfsw_bands.F90
esfsw_driver.o: $(SRCROOT)/am5_phys/atmos_param/radiation/shortwave/esfsw_driver.F90 esfsw_parameters.o esfsw_bands.o esfsw_utilities.o shortwave_types.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)/am5_phys/atmos_param/radiation/shortwave/esfsw_driver.F90
esfsw_parameters.o: $(SRCROOT)/am5_phys/atmos_param/radiation/shortwave/esfsw_parameters.F90
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)/am5_phys/atmos_param/radiation/shortwave/esfsw_parameters.F90
esfsw_utilities.o: $(SRCROOT)/am5_phys/atmos_param/radiation/shortwave/esfsw_utilities.F90
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)/am5_phys/atmos_param/radiation/shortwave/esfsw_utilities.F90
format_input.o: $(SRCROOT)/am5_phys/atmos_param/cosp/quickbeam/format_input.F90 $(SRCROOT)/am5_phys/atmos_param/cosp/cosp_defs.H array_lib.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)/am5_phys/atmos_param/cosp	$(SRCROOT)/am5_phys/atmos_param/cosp/quickbeam/format_input.F90
gamma_mg.o: $(SRCROOT)/am5_phys/atmos_param/microphysics/gamma_mg.F90 lscloud_constants.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)/am5_phys/atmos_param/microphysics/gamma_mg.F90
gas_tf.o: $(SRCROOT)/am5_phys/atmos_param/radiation/longwave/gas_tf.F90 longwave_utilities.o longwave_params.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)/am5_phys/atmos_param/radiation/longwave/gas_tf.F90
gases.o: $(SRCROOT)/am5_phys/atmos_param/cosp/quickbeam/gases.F90 $(SRCROOT)/am5_phys/atmos_param/cosp/cosp_defs.H
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)/am5_phys/atmos_param/cosp	$(SRCROOT)/am5_phys/atmos_param/cosp/quickbeam/gases.F90
get_random_number_stream.o: $(SRCROOT)/am5_phys/atmos_param/radiation/clouds/get_random_number_stream.F90 cloudrad_types.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)/am5_phys/atmos_param/radiation/clouds/get_random_number_stream.F90
gw_beres.o: $(SRCROOT)/am5_phys/atmos_param/gw_beres/gw_beres.F90 gw_common_utils.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)/am5_phys/atmos_param/gw_beres/gw_beres.F90
gw_common_utils.o: $(SRCROOT)/am5_phys/atmos_param/gw_beres/gw_common_utils.F90
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)/am5_phys/atmos_param/gw_beres/gw_common_utils.F90
hanedmf.o: $(SRCROOT)/am5_phys/atmos_param/han_edmf/hanedmf.F90 satmedmfvdifq.o aerosol_types.o aer_ccn_act.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)/am5_phys/atmos_param/han_edmf/hanedmf.F90
icarus.o: $(SRCROOT)/am5_phys/atmos_param/cosp/icarus-scops-3.7/icarus.F90 $(SRCROOT)/am5_phys/atmos_param/cosp/cosp_defs.H
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)/am5_phys/atmos_param/cosp	$(SRCROOT)/am5_phys/atmos_param/cosp/icarus-scops-3.7/icarus.F90
ice_nucl.o: $(SRCROOT)/am5_phys/atmos_shared/tracer_driver/aer_ccn_act/ice_nucl.F90 aer_ccn_act_k.o aerosol_params.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)/am5_phys/atmos_shared/tracer_driver/aer_ccn_act/ice_nucl.F90
isccp_cloud_types.o: $(SRCROOT)/am5_phys/atmos_param/cosp/icarus-scops-3.7/isccp_cloud_types.F $(SRCROOT)/am5_phys/atmos_param/cosp/cosp_defs.H
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)/am5_phys/atmos_param/cosp	$(SRCROOT)/am5_phys/atmos_param/cosp/icarus-scops-3.7/isccp_cloud_types.F
isccp_clouds.o: $(SRCROOT)/am5_phys/atmos_param/radiation/clouds/isccp_clouds.F90
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)/am5_phys/atmos_param/radiation/clouds/isccp_clouds.F90
isoropiaIIcode.o: $(SRCROOT)/am5_phys/atmos_shared/tracer_driver/tropchem/isoropiaIIcode.F $(SRCROOT)/am5_phys/atmos_shared/tracer_driver/tropchem/isrpia.inc
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)/am5_phys/atmos_shared/tracer_driver/tropchem	$(SRCROOT)/am5_phys/atmos_shared/tracer_driver/tropchem/isoropiaIIcode.F
lidar_simulator.o: $(SRCROOT)/am5_phys/atmos_param/cosp/actsim/lidar_simulator.F90 $(SRCROOT)/am5_phys/atmos_param/cosp/cosp_defs.H
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)/am5_phys/atmos_param/cosp	$(SRCROOT)/am5_phys/atmos_param/cosp/actsim/lidar_simulator.F90
llnl_stats.o: $(SRCROOT)/am5_phys/atmos_param/cosp/llnl/llnl_stats.F90 $(SRCROOT)/am5_phys/atmos_param/cosp/cosp_defs.H cosp_constants.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)/am5_phys/atmos_param/cosp	$(SRCROOT)/am5_phys/atmos_param/cosp/llnl/llnl_stats.F90
lmd_ipsl_stats.o: $(SRCROOT)/am5_phys/atmos_param/cosp/actsim/lmd_ipsl_stats.F90 $(SRCROOT)/am5_phys/atmos_param/cosp/cosp_defs.H llnl_stats.o cosp_constants.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)/am5_phys/atmos_param/cosp	$(SRCROOT)/am5_phys/atmos_param/cosp/actsim/lmd_ipsl_stats.F90
longwave_clouds.o: $(SRCROOT)/am5_phys/atmos_param/radiation/longwave/longwave_clouds.F90 longwave_utilities.o longwave_types.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)/am5_phys/atmos_param/radiation/longwave/longwave_clouds.F90
longwave_driver.o: $(SRCROOT)/am5_phys/atmos_param/radiation/driver/longwave_driver.F90 radiation_driver_types.o aerosolrad_types.o sealw99.o longwave_types.o radiative_gases.o radiative_gases_types.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)/am5_phys/atmos_param/radiation/driver/longwave_driver.F90
longwave_fluxes.o: $(SRCROOT)/am5_phys/atmos_param/radiation/longwave/longwave_fluxes.F90 longwave_types.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)/am5_phys/atmos_param/radiation/longwave/longwave_fluxes.F90
longwave_params.o: $(SRCROOT)/am5_phys/atmos_param/radiation/longwave/longwave_params.F90
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)/am5_phys/atmos_param/radiation/longwave/longwave_params.F90
longwave_tables.o: $(SRCROOT)/am5_phys/atmos_param/radiation/longwave/longwave_tables.F90 longwave_utilities.o longwave_params.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)/am5_phys/atmos_param/radiation/longwave/longwave_tables.F90
longwave_types.o: $(SRCROOT)/am5_phys/atmos_param/radiation/longwave/longwave_types.F90
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)/am5_phys/atmos_param/radiation/longwave/longwave_types.F90
longwave_utilities.o: $(SRCROOT)/am5_phys/atmos_param/radiation/longwave/longwave_utilities.F90
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)/am5_phys/atmos_param/radiation/longwave/longwave_utilities.F90
ls_cloud_macrophysics.o: $(SRCROOT)/am5_phys/atmos_param/macrophysics/ls_cloud_macrophysics.F90 physics_types.o tiedtke_macro.o aerosol_cloud.o lscloud_types.o aerosol_types.o moist_processes_utils.o physics_radiation_exch.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)/am5_phys/atmos_param/macrophysics/ls_cloud_macrophysics.F90
ls_cloud_microphysics.o: $(SRCROOT)/am5_phys/atmos_param/microphysics/ls_cloud_microphysics.F90 physics_types.o lscloud_types.o aerosol_types.o physics_radiation_exch.o moist_processes_utils.o lscloud_debug.o rotstayn_klein_mp.o micro_mg2.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)/am5_phys/atmos_param/microphysics/ls_cloud_microphysics.F90
lscale_cond.o: $(SRCROOT)/am5_phys/atmos_param/lscale_cond/lscale_cond.F90
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)/am5_phys/atmos_param/lscale_cond/lscale_cond.F90
lscloud_constants.o: $(SRCROOT)/am5_phys/atmos_param/lscloud_driver/lscloud_constants.F90
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)/am5_phys/atmos_param/lscloud_driver/lscloud_constants.F90
lscloud_debug.o: $(SRCROOT)/am5_phys/atmos_param/lscloud_driver/lscloud_debug.F90 lscloud_types.o moist_processes_utils.o check_nan.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)/am5_phys/atmos_param/lscloud_driver/lscloud_debug.F90
lscloud_driver.o: $(SRCROOT)/am5_phys/atmos_param/lscloud_driver/lscloud_driver.F90 physics_types.o lscale_cond.o lscloud_debug.o lscloud_types.o polysvp.o lscloud_netcdf.o physics_radiation_exch.o ls_cloud_macrophysics.o ls_cloud_microphysics.o aerosol_cloud.o moist_processes_utils.o aerosol_types.o atmos_tracer_utilities.o atmos_global_diag.o atmos_cmip_diag.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)/am5_phys/atmos_param/lscloud_driver/lscloud_driver.F90
lscloud_netcdf.o: $(SRCROOT)/am5_phys/atmos_param/lscloud_driver/lscloud_netcdf.F90 lscloud_types.o atmos_cmip_diag.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)/am5_phys/atmos_param/lscloud_driver/lscloud_netcdf.F90
lscloud_types.o: $(SRCROOT)/am5_phys/atmos_param/lscloud_driver/lscloud_types.F90 atmos_cmip_diag.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)/am5_phys/atmos_param/lscloud_driver/lscloud_types.F90
lw_gases_stdtf.o: $(SRCROOT)/am5_phys/atmos_param/radiation/longwave/lw_gases_stdtf.F90 longwave_utilities.o gas_tf.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)/am5_phys/atmos_param/radiation/longwave/lw_gases_stdtf.F90
m_tracname.o: $(SRCROOT)/am5_phys/atmos_shared/tracer_driver/tropchem/m_tracname.F90 moz.mods.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)/am5_phys/atmos_shared/tracer_driver/tropchem/m_tracname.F90
math_lib.o: $(SRCROOT)/am5_phys/atmos_param/cosp/quickbeam/math_lib.F90 $(SRCROOT)/am5_phys/atmos_param/cosp/cosp_defs.H mrgrnk.o array_lib.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)/am5_phys/atmos_param/cosp	$(SRCROOT)/am5_phys/atmos_param/cosp/quickbeam/math_lib.F90
mfpbltq.o: $(SRCROOT)/am5_phys/atmos_param/han_edmf/mfpbltq.F90 aer_ccn_act_k.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)/am5_phys/atmos_param/han_edmf/mfpbltq.F90
mfscuq.o: $(SRCROOT)/am5_phys/atmos_param/han_edmf/mfscuq.F90
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)/am5_phys/atmos_param/han_edmf/mfscuq.F90
mg_drag.o: $(SRCROOT)/am5_phys/atmos_param/mg_drag/mg_drag.F90
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)/am5_phys/atmos_param/mg_drag/mg_drag.F90
micro_mg2.o: $(SRCROOT)/am5_phys/atmos_param/microphysics/micro_mg2.F90 gamma_mg.o lscloud_types.o physics_radiation_exch.o micro_mg2_utils.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)/am5_phys/atmos_param/microphysics/micro_mg2.F90
micro_mg2_utils.o: $(SRCROOT)/am5_phys/atmos_param/microphysics/micro_mg2_utils.F90 gamma_mg.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)/am5_phys/atmos_param/microphysics/micro_mg2_utils.F90
microphys_rad.o: $(SRCROOT)/am5_phys/atmos_param/radiation/clouds/microphys_rad.F90 cloudrad_types.o sealw99.o esfsw_driver.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)/am5_phys/atmos_param/radiation/clouds/microphys_rad.F90
mo_chem_utls.o: $(SRCROOT)/am5_phys/atmos_shared/tracer_driver/tropchem/mo_chem_utls.F90 moz.mods.o m_tracname.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)/am5_phys/atmos_shared/tracer_driver/tropchem/mo_chem_utls.F90
mo_chemdr.o: $(SRCROOT)/am5_phys/atmos_shared/tracer_driver/tropchem/mo_chemdr.F90 tropchem_types_mod.o mo_chem_utls.o atmos_dust.o moz.mods.o mo_photo.o mo_exp_slv.o mo_imp_slv.o mo_rodas_slv.o mo_usrrxt.o mo_setinv.o moz.subs.o mo_setsox.o mo_fphoto.o strat_chem_utilities.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)/am5_phys/atmos_shared/tracer_driver/tropchem/mo_chemdr.F90
mo_chemini.o: $(SRCROOT)/am5_phys/atmos_shared/tracer_driver/tropchem/mo_chemini.F90 mo_photo.o mo_chem_utls.o mo_usrrxt.o moz.mods.o mo_exp_slv.o mo_imp_slv.o mo_rodas_slv.o mo_read_sim_chm.o mo_fphoto.o tropchem_types_mod.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)/am5_phys/atmos_shared/tracer_driver/tropchem/mo_chemini.F90
mo_exp_slv.o: $(SRCROOT)/am5_phys/atmos_shared/tracer_driver/tropchem/mo_exp_slv.F90 moz.mods.o mo_chem_utls.o moz.mat.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)/am5_phys/atmos_shared/tracer_driver/tropchem/mo_exp_slv.F90
mo_fastjx.o: $(SRCROOT)/am5_phys/atmos_shared/tracer_driver/tropchem/AM4_chem/mo_fastjx.F90 tropchem_types_mod.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)/am5_phys/atmos_shared/tracer_driver/tropchem/AM4_chem/mo_fastjx.F90
mo_fphoto.o: $(SRCROOT)/am5_phys/atmos_shared/tracer_driver/tropchem/AM4_chem/mo_fphoto.F90 mo_fastjx.o tropchem_types_mod.o aerosol_thermodynamics.o mo_chem_utls.o moz.mods.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)/am5_phys/atmos_shared/tracer_driver/tropchem/AM4_chem/mo_fphoto.F90
mo_hook.o: $(SRCROOT)/am5_phys/atmos_shared/tracer_driver/tropchem/mo_hook.F90 atmos_cmip_diag.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)/am5_phys/atmos_shared/tracer_driver/tropchem/mo_hook.F90
mo_imp_slv.o: $(SRCROOT)/am5_phys/atmos_shared/tracer_driver/tropchem/mo_imp_slv.F90 moz.mods.o mo_chem_utls.o m_tracname.o moz.mat.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)/am5_phys/atmos_shared/tracer_driver/tropchem/mo_imp_slv.F90
mo_jpl.o: $(SRCROOT)/am5_phys/atmos_shared/tracer_driver/tropchem/mo_jpl.F90
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)/am5_phys/atmos_shared/tracer_driver/tropchem/mo_jpl.F90
mo_photo.o: $(SRCROOT)/am5_phys/atmos_shared/tracer_driver/tropchem/mo_photo.F90 mo_chem_utls.o moz.mods.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)/am5_phys/atmos_shared/tracer_driver/tropchem/mo_photo.F90
mo_read_sim_chm.o: $(SRCROOT)/am5_phys/atmos_shared/tracer_driver/tropchem/mo_read_sim_chm.F90 moz.mods.o m_tracname.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)/am5_phys/atmos_shared/tracer_driver/tropchem/mo_read_sim_chm.F90
mo_rodas_slv.o: $(SRCROOT)/am5_phys/atmos_shared/tracer_driver/tropchem/mo_rodas_slv.F90 moz.mods.o mo_chem_utls.o moz.mat.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)/am5_phys/atmos_shared/tracer_driver/tropchem/mo_rodas_slv.F90
mo_setinv.o: $(SRCROOT)/am5_phys/atmos_shared/tracer_driver/tropchem/mo_setinv.F90
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)/am5_phys/atmos_shared/tracer_driver/tropchem/mo_setinv.F90
mo_setsox.o: $(SRCROOT)/am5_phys/atmos_shared/tracer_driver/mo_setsox.F90 mo_chem_utls.o tropchem_types_mod.o cloud_chem.o aerosol_thermodynamics.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)/am5_phys/atmos_shared/tracer_driver/mo_setsox.F90
mo_usrrxt.o: $(SRCROOT)/am5_phys/atmos_shared/tracer_driver/tropchem/mo_usrrxt.F90 aerosol_thermodynamics.o strat_chem_utilities.o tropchem_types_mod.o atmos_dust.o mo_chem_utls.o moz.mods.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)/am5_phys/atmos_shared/tracer_driver/tropchem/mo_usrrxt.F90
modis_simulator.o: $(SRCROOT)/am5_phys/atmos_param/cosp/MODIS_simulator/modis_simulator.F90 $(SRCROOT)/am5_phys/atmos_param/cosp/cosp_defs.H cosp_types.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)/am5_phys/atmos_param/cosp	$(SRCROOT)/am5_phys/atmos_param/cosp/MODIS_simulator/modis_simulator.F90
moist_processes.o: $(SRCROOT)/am5_phys/atmos_param/moist_processes/moist_processes.F90 physics_types.o physics_radiation_exch.o lscloud_driver.o convection_driver.o convection_utilities.o atmos_global_diag.o vert_diff_driver.o aerosol_types.o atmos_tracer_utilities.o moist_processes_utils.o atmos_dust.o atmos_tracer_driver.o atmos_sea_salt.o atmos_cmip_diag.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)/am5_phys/atmos_param/moist_processes/moist_processes.F90
moist_processes_utils.o: $(SRCROOT)/am5_phys/atmos_param/moist_processes/moist_processes_utils.F90 lscloud_types.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)/am5_phys/atmos_param/moist_processes/moist_processes_utils.F90
moz.mat.o: $(SRCROOT)/am5_phys/atmos_shared/tracer_driver/tropchem/AM4_chem/moz.mat.F90 moz.mods.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)/am5_phys/atmos_shared/tracer_driver/tropchem/AM4_chem/moz.mat.F90
moz.mods.o: $(SRCROOT)/am5_phys/atmos_shared/tracer_driver/tropchem/AM4_chem/moz.mods.F90
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)/am5_phys/atmos_shared/tracer_driver/tropchem/AM4_chem/moz.mods.F90
moz.subs.o: $(SRCROOT)/am5_phys/atmos_shared/tracer_driver/tropchem/AM4_chem/moz.subs.F90 moz.mods.o mo_jpl.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)/am5_phys/atmos_shared/tracer_driver/tropchem/AM4_chem/moz.subs.F90
mrgrnk.o: $(SRCROOT)/am5_phys/atmos_param/cosp/quickbeam/mrgrnk.F90 $(SRCROOT)/am5_phys/atmos_param/cosp/cosp_defs.H
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)/am5_phys/atmos_param/cosp	$(SRCROOT)/am5_phys/atmos_param/cosp/quickbeam/mrgrnk.F90
optical_path.o: $(SRCROOT)/am5_phys/atmos_param/radiation/longwave/optical_path.F90 longwave_utilities.o longwave_params.o lw_gases_stdtf.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)/am5_phys/atmos_param/radiation/longwave/optical_path.F90
optics_lib.o: $(SRCROOT)/am5_phys/atmos_param/cosp/quickbeam/optics_lib.F90 $(SRCROOT)/am5_phys/atmos_param/cosp/cosp_defs.H
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)/am5_phys/atmos_param/cosp	$(SRCROOT)/am5_phys/atmos_param/cosp/quickbeam/optics_lib.F90
ozone.o: $(SRCROOT)/am5_phys/atmos_param/radiation/radiative_gases/ozone.F90 radiative_gases_types.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)/am5_phys/atmos_param/radiation/radiative_gases/ozone.F90
pf_to_mr.o: $(SRCROOT)/am5_phys/atmos_param/cosp/llnl/pf_to_mr.F $(SRCROOT)/am5_phys/atmos_param/cosp/cosp_defs.H
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)/am5_phys/atmos_param/cosp	$(SRCROOT)/am5_phys/atmos_param/cosp/llnl/pf_to_mr.F
physics_driver.o: $(SRCROOT)/am5_phys/atmos_param/physics_driver/physics_driver.F90 atmos_tracer_driver.o atmos_cmip_diag.o aerosol_types.o physics_radiation_exch.o physics_types.o moist_processes_utils.o aerosol.o cosp_driver.o moist_processes.o vert_turb_driver.o vert_diff_driver.o damping_driver.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)/am5_phys/atmos_param/physics_driver/physics_driver.F90
physics_radiation_exch.o: $(SRCROOT)/am5_phys/atmos_param/physics_radiation_exch/physics_radiation_exch.F90
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)/am5_phys/atmos_param/physics_radiation_exch/physics_radiation_exch.F90
physics_types.o: $(SRCROOT)/am5_phys/atmos_param/physics_driver/physics_types.F90
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)/am5_phys/atmos_param/physics_driver/physics_types.F90
polysvp.o: $(SRCROOT)/am5_phys/atmos_param/lscloud_driver/polysvp.F90 lscloud_constants.o lscloud_types.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)/am5_phys/atmos_param/lscloud_driver/polysvp.F90
prec_scops.o: $(SRCROOT)/am5_phys/atmos_param/cosp/llnl/prec_scops.F $(SRCROOT)/am5_phys/atmos_param/cosp/cosp_defs.H
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)/am5_phys/atmos_param/cosp	$(SRCROOT)/am5_phys/atmos_param/cosp/llnl/prec_scops.F
rad_output_file.o: $(SRCROOT)/am5_phys/atmos_param/radiation/driver/rad_output_file.F90 aerosol_types.o atmos_cmip_diag.o radiation_driver_types.o aerosolrad_types.o cloudrad_types.o esfsw_driver.o shortwave_types.o longwave_types.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)/am5_phys/atmos_param/radiation/driver/rad_output_file.F90
radar_simulator.o: $(SRCROOT)/am5_phys/atmos_param/cosp/quickbeam/radar_simulator.F90 $(SRCROOT)/am5_phys/atmos_param/cosp/cosp_defs.H mrgrnk.o array_lib.o math_lib.o optics_lib.o radar_simulator_types.o scale_LUTs_io.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)/am5_phys/atmos_param/cosp	$(SRCROOT)/am5_phys/atmos_param/cosp/quickbeam/radar_simulator.F90
radar_simulator_init.o: $(SRCROOT)/am5_phys/atmos_param/cosp/quickbeam/radar_simulator_init.F90 $(SRCROOT)/am5_phys/atmos_param/cosp/cosp_defs.H radar_simulator_types.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)/am5_phys/atmos_param/cosp	$(SRCROOT)/am5_phys/atmos_param/cosp/quickbeam/radar_simulator_init.F90
radar_simulator_types.o: $(SRCROOT)/am5_phys/atmos_param/cosp/quickbeam/radar_simulator_types.F90 $(SRCROOT)/am5_phys/atmos_param/cosp/cosp_defs.H
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)/am5_phys/atmos_param/cosp	$(SRCROOT)/am5_phys/atmos_param/cosp/quickbeam/radar_simulator_types.F90
radiation_diag.o: $(SRCROOT)/am5_phys/atmos_param/radiation/util/radiation_diag.F90 radiation_driver.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)/am5_phys/atmos_param/radiation/util/radiation_diag.F90
radiation_driver.o: $(SRCROOT)/am5_phys/atmos_param/radiation/driver/radiation_driver.F90 radiation_types.o physics_radiation_exch.o radiation_driver_types.o aerosol_types.o aerosolrad_types.o cloudrad_types.o longwave_driver.o longwave_types.o shortwave_driver.o shortwave_types.o rad_output_file.o aerosolrad_driver.o cloudrad_driver.o cloudrad_diagnostics.o radiative_gases.o radiative_gases_types.o radiation_driver_diag.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)/am5_phys/atmos_param/radiation/driver/radiation_driver.F90
radiation_driver_diag.o: $(SRCROOT)/am5_phys/atmos_param/radiation/driver/radiation_driver_diag.F90 atmos_cmip_diag.o atmos_global_diag.o radiation_driver_types.o aerosolrad_types.o shortwave_types.o longwave_types.o radiative_gases_types.o shortwave_driver.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)/am5_phys/atmos_param/radiation/driver/radiation_driver_diag.F90
radiation_driver_types.o: $(SRCROOT)/am5_phys/atmos_param/radiation/driver/radiation_driver_types.F90
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)/am5_phys/atmos_param/radiation/driver/radiation_driver_types.F90
radiation_types.o: $(SRCROOT)/am5_phys/atmos_param/radiation/driver/radiation_types.F90 atmos_co2.o atmos_ch4.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)/am5_phys/atmos_param/radiation/driver/radiation_types.F90
radiative_gases.o: $(SRCROOT)/am5_phys/atmos_param/radiation/radiative_gases/radiative_gases.F90 sealw99.o radiative_gases_types.o ozone.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)/am5_phys/atmos_param/radiation/radiative_gases/radiative_gases.F90
radiative_gases_types.o: $(SRCROOT)/am5_phys/atmos_param/radiation/radiative_gases/radiative_gases_types.F90
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)/am5_phys/atmos_param/radiation/radiative_gases/radiative_gases_types.F90
rotstayn_klein_mp.o: $(SRCROOT)/am5_phys/atmos_param/microphysics/rotstayn_klein_mp.F90 cloud_generator.o aer_in_act.o polysvp.o lscloud_types.o lscloud_debug.o physics_radiation_exch.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)/am5_phys/atmos_param/microphysics/rotstayn_klein_mp.F90
satmedmfvdifq.o: $(SRCROOT)/am5_phys/atmos_param/han_edmf/satmedmfvdifq.F90 tridi.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)/am5_phys/atmos_param/han_edmf/satmedmfvdifq.F90
scale_LUTs_io.o: $(SRCROOT)/am5_phys/atmos_param/cosp/quickbeam/scale_LUTs_io.F90 $(SRCROOT)/am5_phys/atmos_param/cosp/cosp_defs.H radar_simulator_types.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)/am5_phys/atmos_param/cosp	$(SRCROOT)/am5_phys/atmos_param/cosp/quickbeam/scale_LUTs_io.F90
scops.o: $(SRCROOT)/am5_phys/atmos_param/cosp/icarus-scops-3.7/scops.F $(SRCROOT)/am5_phys/atmos_param/cosp/cosp_defs.H $(SRCROOT)/am5_phys/atmos_param/cosp/icarus-scops-3.7/congvec.H
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)/am5_phys/atmos_param/cosp -I$(SRCROOT)/am5_phys/atmos_param/cosp/icarus-scops-3.7	$(SRCROOT)/am5_phys/atmos_param/cosp/icarus-scops-3.7/scops.F
sealw99.o: $(SRCROOT)/am5_phys/atmos_param/radiation/longwave/sealw99.F90 longwave_utilities.o longwave_types.o longwave_params.o longwave_clouds.o longwave_fluxes.o longwave_tables.o optical_path.o gas_tf.o lw_gases_stdtf.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)/am5_phys/atmos_param/radiation/longwave/sealw99.F90
shortwave_driver.o: $(SRCROOT)/am5_phys/atmos_param/radiation/driver/shortwave_driver.F90 radiation_driver_types.o esfsw_driver.o shortwave_types.o radiative_gases_types.o solar_data_driver.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)/am5_phys/atmos_param/radiation/driver/shortwave_driver.F90
shortwave_types.o: $(SRCROOT)/am5_phys/atmos_param/radiation/shortwave/shortwave_types.F90
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)/am5_phys/atmos_param/radiation/shortwave/shortwave_types.F90
solar_data_driver.o: $(SRCROOT)/am5_phys/atmos_param/radiation/driver/solar_data_driver.F90
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)/am5_phys/atmos_param/radiation/driver/solar_data_driver.F90
stable_bl_turb.o: $(SRCROOT)/am5_phys/atmos_param/stable_bl_turb/stable_bl_turb.F90
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)/am5_phys/atmos_param/stable_bl_turb/stable_bl_turb.F90
strat_chem_driver.o: $(SRCROOT)/am5_phys/atmos_shared/tracer_driver/stratchem/strat_chem_driver.F90 strat_chem_model.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)/am5_phys/atmos_shared/tracer_driver/stratchem/strat_chem_driver.F90
strat_chem_model.o: $(SRCROOT)/am5_phys/atmos_shared/tracer_driver/stratchem/strat_chem_model.F90
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)/am5_phys/atmos_shared/tracer_driver/stratchem/strat_chem_model.F90
strat_chem_utilities.o: $(SRCROOT)/am5_phys/atmos_shared/tracer_driver/tropchem/strat_chem_utilities.F90
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)/am5_phys/atmos_shared/tracer_driver/tropchem/strat_chem_utilities.F90
strat_clouds_W.o: $(SRCROOT)/am5_phys/atmos_param/radiation/clouds/strat_clouds_W.F90 aerosol_types.o physics_radiation_exch.o cloud_rad.o cloudrad_types.o get_random_number_stream.o cloud_generator.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)/am5_phys/atmos_param/radiation/clouds/strat_clouds_W.F90
tiedtke_macro.o: $(SRCROOT)/am5_phys/atmos_param/macrophysics/tiedtke_macro.F90 betaDistribution.o lscloud_types.o lscloud_debug.o polysvp.o moist_processes_utils.o physics_radiation_exch.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)/am5_phys/atmos_param/macrophysics/tiedtke_macro.F90
topo_drag.o: $(SRCROOT)/am5_phys/atmos_param/topo_drag/topo_drag.F90
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)/am5_phys/atmos_param/topo_drag/topo_drag.F90
tridi.o: $(SRCROOT)/am5_phys/atmos_param/han_edmf/tridi.F90
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)/am5_phys/atmos_param/han_edmf/tridi.F90
tropchem_driver.o: $(SRCROOT)/am5_phys/atmos_shared/tracer_driver/tropchem_driver.F90 tropchem_types_mod.o atmos_cmip_diag.o atmos_tracer_utilities.o mo_chemdr.o mo_setsox.o mo_chemini.o m_tracname.o moz.mods.o mo_hook.o strat_chem_utilities.o mo_chem_utls.o atmos_sulfate.o shortwave_driver.o cloud_chem.o aerosol_thermodynamics.o mo_usrrxt.o atmos_dust.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)/am5_phys/atmos_shared/tracer_driver/tropchem_driver.F90
tropchem_types_mod.o: $(SRCROOT)/am5_phys/atmos_shared/tracer_driver/tropchem/tropchem_types_mod.F90
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)/am5_phys/atmos_shared/tracer_driver/tropchem/tropchem_types_mod.F90
uw_clouds_W.o: $(SRCROOT)/am5_phys/atmos_param/radiation/clouds/uw_clouds_W.F90 physics_radiation_exch.o cloudrad_types.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)/am5_phys/atmos_param/radiation/clouds/uw_clouds_W.F90
uw_conv.o: $(SRCROOT)/am5_phys/atmos_param/shallow_cu/uw_conv.F90 $(SRCROOT)FMS/include/fms_platform.h atmos_cmip_diag.o atmos_tracer_utilities.o moist_processes_utils.o aerosol_types.o aer_ccn_act.o conv_utilities.o conv_utilities_k.o conv_plumes_k.o conv_closures.o deep_conv.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)FMS/include	$(SRCROOT)/am5_phys/atmos_param/shallow_cu/uw_conv.F90
vert_advection.o: $(SRCROOT)/am5_phys/atmos_shared/vert_advection/vert_advection.F90
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)/am5_phys/atmos_shared/vert_advection/vert_advection.F90
vert_diff.o: $(SRCROOT)/am5_phys/atmos_param/vert_diff/vert_diff.F90
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)/am5_phys/atmos_param/vert_diff/vert_diff.F90
vert_diff_driver.o: $(SRCROOT)/am5_phys/atmos_param/vert_diff_driver/vert_diff_driver.F90 vert_diff.o atmos_cmip_diag.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)/am5_phys/atmos_param/vert_diff_driver/vert_diff_driver.F90
vert_turb_driver.o: $(SRCROOT)/am5_phys/atmos_param/vert_turb_driver/vert_turb_driver.F90 physics_radiation_exch.o physics_types.o stable_bl_turb.o entrain.o hanedmf.o aerosol_types.o moist_processes_utils.o atmos_cmip_diag.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)/am5_phys/atmos_param/vert_turb_driver/vert_turb_driver.F90
wet_deposition_0D.o: $(SRCROOT)/am5_phys/atmos_param/shallow_cu/wet_deposition_0D.F90
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)/am5_phys/atmos_param/shallow_cu/wet_deposition_0D.F90
xactive_bvoc.o: $(SRCROOT)/am5_phys/atmos_shared/tracer_driver/xactive_bvoc.F90 moz.mods.o m_tracname.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)/am5_phys/atmos_shared/tracer_driver/xactive_bvoc.F90
zeff.o: $(SRCROOT)/am5_phys/atmos_param/cosp/quickbeam/zeff.F90 $(SRCROOT)/am5_phys/atmos_param/cosp/cosp_defs.H math_lib.o optics_lib.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)/am5_phys/atmos_param/cosp	$(SRCROOT)/am5_phys/atmos_param/cosp/quickbeam/zeff.F90
./atmos_dust.F90: $(SRCROOT)/am5_phys/atmos_shared/tracer_driver/atmos_dust.F90
	cp $(SRCROOT)/am5_phys/atmos_shared/tracer_driver/atmos_dust.F90 .
./cosp_stats.F90: $(SRCROOT)/am5_phys/atmos_param/cosp/cosp_stats.F90
	cp $(SRCROOT)/am5_phys/atmos_param/cosp/cosp_stats.F90 .
./get_random_number_stream.F90: $(SRCROOT)/am5_phys/atmos_param/radiation/clouds/get_random_number_stream.F90
	cp $(SRCROOT)/am5_phys/atmos_param/radiation/clouds/get_random_number_stream.F90 .
./cosp.F90: $(SRCROOT)/am5_phys/atmos_param/cosp/cosp.F90
	cp $(SRCROOT)/am5_phys/atmos_param/cosp/cosp.F90 .
./lscale_cond.F90: $(SRCROOT)/am5_phys/atmos_param/lscale_cond/lscale_cond.F90
	cp $(SRCROOT)/am5_phys/atmos_param/lscale_cond/lscale_cond.F90 .
./prec_scops.F: $(SRCROOT)/am5_phys/atmos_param/cosp/llnl/prec_scops.F
	cp $(SRCROOT)/am5_phys/atmos_param/cosp/llnl/prec_scops.F .
./radiation_types.F90: $(SRCROOT)/am5_phys/atmos_param/radiation/driver/radiation_types.F90
	cp $(SRCROOT)/am5_phys/atmos_param/radiation/driver/radiation_types.F90 .
./satmedmfvdifq.F90: $(SRCROOT)/am5_phys/atmos_param/han_edmf/satmedmfvdifq.F90
	cp $(SRCROOT)/am5_phys/atmos_param/han_edmf/satmedmfvdifq.F90 .
./hanedmf.F90: $(SRCROOT)/am5_phys/atmos_param/han_edmf/hanedmf.F90
	cp $(SRCROOT)/am5_phys/atmos_param/han_edmf/hanedmf.F90 .
./wet_deposition_0D.F90: $(SRCROOT)/am5_phys/atmos_param/shallow_cu/wet_deposition_0D.F90
	cp $(SRCROOT)/am5_phys/atmos_param/shallow_cu/wet_deposition_0D.F90 .
./conv_utilities_k.F90: $(SRCROOT)/am5_phys/atmos_param/shallow_cu/conv_utilities_k.F90
	cp $(SRCROOT)/am5_phys/atmos_param/shallow_cu/conv_utilities_k.F90 .
./optical_path.F90: $(SRCROOT)/am5_phys/atmos_param/radiation/longwave/optical_path.F90
	cp $(SRCROOT)/am5_phys/atmos_param/radiation/longwave/optical_path.F90 .
./cosp_io.F90: $(SRCROOT)/am5_phys/atmos_param/cosp/cosp_io.F90
	cp $(SRCROOT)/am5_phys/atmos_param/cosp/cosp_io.F90 .
./gw_common_utils.F90: $(SRCROOT)/am5_phys/atmos_param/gw_beres/gw_common_utils.F90
	cp $(SRCROOT)/am5_phys/atmos_param/gw_beres/gw_common_utils.F90 .
./ozone.F90: $(SRCROOT)/am5_phys/atmos_param/radiation/radiative_gases/ozone.F90
	cp $(SRCROOT)/am5_phys/atmos_param/radiation/radiative_gases/ozone.F90 .
./radiation_driver.F90: $(SRCROOT)/am5_phys/atmos_param/radiation/driver/radiation_driver.F90
	cp $(SRCROOT)/am5_phys/atmos_param/radiation/driver/radiation_driver.F90 .
./cloudrad_driver.F90: $(SRCROOT)/am5_phys/atmos_param/radiation/clouds/cloudrad_driver.F90
	cp $(SRCROOT)/am5_phys/atmos_param/radiation/clouds/cloudrad_driver.F90 .
./esfsw_bands.F90: $(SRCROOT)/am5_phys/atmos_param/radiation/shortwave/esfsw_bands.F90
	cp $(SRCROOT)/am5_phys/atmos_param/radiation/shortwave/esfsw_bands.F90 .
./atmos_ozone_tracer.F90: $(SRCROOT)/am5_phys/atmos_shared/tracer_driver/atmos_ozone_tracer.F90
	cp $(SRCROOT)/am5_phys/atmos_shared/tracer_driver/atmos_ozone_tracer.F90 .
./shortwave_driver.F90: $(SRCROOT)/am5_phys/atmos_param/radiation/driver/shortwave_driver.F90
	cp $(SRCROOT)/am5_phys/atmos_param/radiation/driver/shortwave_driver.F90 .
./rotstayn_klein_mp.F90: $(SRCROOT)/am5_phys/atmos_param/microphysics/rotstayn_klein_mp.F90
	cp $(SRCROOT)/am5_phys/atmos_param/microphysics/rotstayn_klein_mp.F90 .
./lscloud_netcdf.F90: $(SRCROOT)/am5_phys/atmos_param/lscloud_driver/lscloud_netcdf.F90
	cp $(SRCROOT)/am5_phys/atmos_param/lscloud_driver/lscloud_netcdf.F90 .
./vert_diff.F90: $(SRCROOT)/am5_phys/atmos_param/vert_diff/vert_diff.F90
	cp $(SRCROOT)/am5_phys/atmos_param/vert_diff/vert_diff.F90 .
./atmos_convection_tracer.F90: $(SRCROOT)/am5_phys/atmos_shared/tracer_driver/atmos_convection_tracer.F90
	cp $(SRCROOT)/am5_phys/atmos_shared/tracer_driver/atmos_convection_tracer.F90 .
./atmos_regional_tracer_driver.F90: $(SRCROOT)/am5_phys/atmos_shared/tracer_driver/atmos_regional_tracer_driver.F90
	cp $(SRCROOT)/am5_phys/atmos_shared/tracer_driver/atmos_regional_tracer_driver.F90 .
./aer_ccn_act_k.F90: $(SRCROOT)/am5_phys/atmos_shared/tracer_driver/aer_ccn_act/aer_ccn_act_k.F90
	cp $(SRCROOT)/am5_phys/atmos_shared/tracer_driver/aer_ccn_act/aer_ccn_act_k.F90 .
./lscloud_types.F90: $(SRCROOT)/am5_phys/atmos_param/lscloud_driver/lscloud_types.F90
	cp $(SRCROOT)/am5_phys/atmos_param/lscloud_driver/lscloud_types.F90 .
./check_nan.F90: $(SRCROOT)/am5_phys/atmos_param/lscloud_driver/check_nan.F90
	cp $(SRCROOT)/am5_phys/atmos_param/lscloud_driver/check_nan.F90 .
./isoropiaIIcode.F: $(SRCROOT)/am5_phys/atmos_shared/tracer_driver/tropchem/isoropiaIIcode.F
	cp $(SRCROOT)/am5_phys/atmos_shared/tracer_driver/tropchem/isoropiaIIcode.F .
./mo_chemdr.F90: $(SRCROOT)/am5_phys/atmos_shared/tracer_driver/tropchem/mo_chemdr.F90
	cp $(SRCROOT)/am5_phys/atmos_shared/tracer_driver/tropchem/mo_chemdr.F90 .
./moz.mat.F90: $(SRCROOT)/am5_phys/atmos_shared/tracer_driver/tropchem/AM4_chem/moz.mat.F90
	cp $(SRCROOT)/am5_phys/atmos_shared/tracer_driver/tropchem/AM4_chem/moz.mat.F90 .
./solar_data_driver.F90: $(SRCROOT)/am5_phys/atmos_param/radiation/driver/solar_data_driver.F90
	cp $(SRCROOT)/am5_phys/atmos_param/radiation/driver/solar_data_driver.F90 .
./atmos_sea_salt.F90: $(SRCROOT)/am5_phys/atmos_shared/tracer_driver/atmos_sea_salt.F90
	cp $(SRCROOT)/am5_phys/atmos_shared/tracer_driver/atmos_sea_salt.F90 .
./atmos_carbon_aerosol.F90: $(SRCROOT)/am5_phys/atmos_shared/tracer_driver/atmos_carbon_aerosol.F90
	cp $(SRCROOT)/am5_phys/atmos_shared/tracer_driver/atmos_carbon_aerosol.F90 .
./atmos_sulfur_hex.F90: $(SRCROOT)/am5_phys/atmos_shared/tracer_driver/atmos_sulfur_hex.F90
	cp $(SRCROOT)/am5_phys/atmos_shared/tracer_driver/atmos_sulfur_hex.F90 .
./mrgrnk.F90: $(SRCROOT)/am5_phys/atmos_param/cosp/quickbeam/mrgrnk.F90
	cp $(SRCROOT)/am5_phys/atmos_param/cosp/quickbeam/mrgrnk.F90 .
./scale_LUTs_io.F90: $(SRCROOT)/am5_phys/atmos_param/cosp/quickbeam/scale_LUTs_io.F90
	cp $(SRCROOT)/am5_phys/atmos_param/cosp/quickbeam/scale_LUTs_io.F90 .
./stable_bl_turb.F90: $(SRCROOT)/am5_phys/atmos_param/stable_bl_turb/stable_bl_turb.F90
	cp $(SRCROOT)/am5_phys/atmos_param/stable_bl_turb/stable_bl_turb.F90 .
./betaDistribution.F90: $(SRCROOT)/am5_phys/atmos_param/cloud_generator/betaDistribution.F90
	cp $(SRCROOT)/am5_phys/atmos_param/cloud_generator/betaDistribution.F90 .
./cosp_driver.F90: $(SRCROOT)/am5_phys/atmos_param/cosp/cosp_driver.F90
	cp $(SRCROOT)/am5_phys/atmos_param/cosp/cosp_driver.F90 .
./atmos_tracer_utilities.F90: $(SRCROOT)/am5_phys/atmos_shared/tracer_driver/atmos_tracer_utilities.F90
	cp $(SRCROOT)/am5_phys/atmos_shared/tracer_driver/atmos_tracer_utilities.F90 .
./mo_imp_slv.F90: $(SRCROOT)/am5_phys/atmos_shared/tracer_driver/tropchem/mo_imp_slv.F90
	cp $(SRCROOT)/am5_phys/atmos_shared/tracer_driver/tropchem/mo_imp_slv.F90 .
./tiedtke_macro.F90: $(SRCROOT)/am5_phys/atmos_param/macrophysics/tiedtke_macro.F90
	cp $(SRCROOT)/am5_phys/atmos_param/macrophysics/tiedtke_macro.F90 .
./atmos_cmip_interp.inc: $(SRCROOT)/am5_phys/atmos_shared/atmos_cmip_diag/atmos_cmip_interp.inc
	cp $(SRCROOT)/am5_phys/atmos_shared/atmos_cmip_diag/atmos_cmip_interp.inc .
./radiation_diag.F90: $(SRCROOT)/am5_phys/atmos_param/radiation/util/radiation_diag.F90
	cp $(SRCROOT)/am5_phys/atmos_param/radiation/util/radiation_diag.F90 .
./conv_utilities.F90: $(SRCROOT)/am5_phys/atmos_param/shallow_cu/conv_utilities.F90
	cp $(SRCROOT)/am5_phys/atmos_param/shallow_cu/conv_utilities.F90 .
./esfsw_utilities.F90: $(SRCROOT)/am5_phys/atmos_param/radiation/shortwave/esfsw_utilities.F90
	cp $(SRCROOT)/am5_phys/atmos_param/radiation/shortwave/esfsw_utilities.F90 .
./gases.F90: $(SRCROOT)/am5_phys/atmos_param/cosp/quickbeam/gases.F90
	cp $(SRCROOT)/am5_phys/atmos_param/cosp/quickbeam/gases.F90 .
./convection_driver.F90: $(SRCROOT)/am5_phys/atmos_param/convection_driver/convection_driver.F90
	cp $(SRCROOT)/am5_phys/atmos_param/convection_driver/convection_driver.F90 .
./isccp_clouds.F90: $(SRCROOT)/am5_phys/atmos_param/radiation/clouds/isccp_clouds.F90
	cp $(SRCROOT)/am5_phys/atmos_param/radiation/clouds/isccp_clouds.F90 .
./array_lib.F90: $(SRCROOT)/am5_phys/atmos_param/cosp/quickbeam/array_lib.F90
	cp $(SRCROOT)/am5_phys/atmos_param/cosp/quickbeam/array_lib.F90 .
./atmos_sulfate.F90: $(SRCROOT)/am5_phys/atmos_shared/tracer_driver/atmos_sulfate.F90
	cp $(SRCROOT)/am5_phys/atmos_shared/tracer_driver/atmos_sulfate.F90 .
./mo_hook.F90: $(SRCROOT)/am5_phys/atmos_shared/tracer_driver/tropchem/mo_hook.F90
	cp $(SRCROOT)/am5_phys/atmos_shared/tracer_driver/tropchem/mo_hook.F90 .
./conv_closures.F90: $(SRCROOT)/am5_phys/atmos_param/shallow_cu/conv_closures.F90
	cp $(SRCROOT)/am5_phys/atmos_param/shallow_cu/conv_closures.F90 .
./mo_jpl.F90: $(SRCROOT)/am5_phys/atmos_shared/tracer_driver/tropchem/mo_jpl.F90
	cp $(SRCROOT)/am5_phys/atmos_shared/tracer_driver/tropchem/mo_jpl.F90 .
./radiation_driver_types.F90: $(SRCROOT)/am5_phys/atmos_param/radiation/driver/radiation_driver_types.F90
	cp $(SRCROOT)/am5_phys/atmos_param/radiation/driver/radiation_driver_types.F90 .
./llnl_stats.F90: $(SRCROOT)/am5_phys/atmos_param/cosp/llnl/llnl_stats.F90
	cp $(SRCROOT)/am5_phys/atmos_param/cosp/llnl/llnl_stats.F90 .
./isccp_cloud_types.F: $(SRCROOT)/am5_phys/atmos_param/cosp/icarus-scops-3.7/isccp_cloud_types.F
	cp $(SRCROOT)/am5_phys/atmos_param/cosp/icarus-scops-3.7/isccp_cloud_types.F .
./optics_lib.F90: $(SRCROOT)/am5_phys/atmos_param/cosp/quickbeam/optics_lib.F90
	cp $(SRCROOT)/am5_phys/atmos_param/cosp/quickbeam/optics_lib.F90 .
./esfsw_parameters.F90: $(SRCROOT)/am5_phys/atmos_param/radiation/shortwave/esfsw_parameters.F90
	cp $(SRCROOT)/am5_phys/atmos_param/radiation/shortwave/esfsw_parameters.F90 .
./topo_drag.F90: $(SRCROOT)/am5_phys/atmos_param/topo_drag/topo_drag.F90
	cp $(SRCROOT)/am5_phys/atmos_param/topo_drag/topo_drag.F90 .
./micro_mg2.F90: $(SRCROOT)/am5_phys/atmos_param/microphysics/micro_mg2.F90
	cp $(SRCROOT)/am5_phys/atmos_param/microphysics/micro_mg2.F90 .
./m_tracname.F90: $(SRCROOT)/am5_phys/atmos_shared/tracer_driver/tropchem/m_tracname.F90
	cp $(SRCROOT)/am5_phys/atmos_shared/tracer_driver/tropchem/m_tracname.F90 .
./lscloud_driver.F90: $(SRCROOT)/am5_phys/atmos_param/lscloud_driver/lscloud_driver.F90
	cp $(SRCROOT)/am5_phys/atmos_param/lscloud_driver/lscloud_driver.F90 .
./cloud_rad.F90: $(SRCROOT)/am5_phys/atmos_param/cloud_rad/cloud_rad.F90
	cp $(SRCROOT)/am5_phys/atmos_param/cloud_rad/cloud_rad.F90 .
./aerosolrad_types.F90: $(SRCROOT)/am5_phys/atmos_param/radiation/aerosols/aerosolrad_types.F90
	cp $(SRCROOT)/am5_phys/atmos_param/radiation/aerosols/aerosolrad_types.F90 .
./cosp_simulator.F90: $(SRCROOT)/am5_phys/atmos_param/cosp/cosp_simulator.F90
	cp $(SRCROOT)/am5_phys/atmos_param/cosp/cosp_simulator.F90 .
./cosp_radar.F90: $(SRCROOT)/am5_phys/atmos_param/cosp/llnl/cosp_radar.F90
	cp $(SRCROOT)/am5_phys/atmos_param/cosp/llnl/cosp_radar.F90 .
./cosp_isccp_simulator.F90: $(SRCROOT)/am5_phys/atmos_param/cosp/cosp_isccp_simulator.F90
	cp $(SRCROOT)/am5_phys/atmos_param/cosp/cosp_isccp_simulator.F90 .
./mo_setinv.F90: $(SRCROOT)/am5_phys/atmos_shared/tracer_driver/tropchem/mo_setinv.F90
	cp $(SRCROOT)/am5_phys/atmos_shared/tracer_driver/tropchem/mo_setinv.F90 .
./deep_conv.F90: $(SRCROOT)/am5_phys/atmos_param/shallow_cu/deep_conv.F90
	cp $(SRCROOT)/am5_phys/atmos_param/shallow_cu/deep_conv.F90 .
./vert_advection.F90: $(SRCROOT)/am5_phys/atmos_shared/vert_advection/vert_advection.F90
	cp $(SRCROOT)/am5_phys/atmos_shared/vert_advection/vert_advection.F90 .
./cosp_modis_simulator.F90: $(SRCROOT)/am5_phys/atmos_param/cosp/cosp_modis_simulator.F90
	cp $(SRCROOT)/am5_phys/atmos_param/cosp/cosp_modis_simulator.F90 .
./moist_processes.F90: $(SRCROOT)/am5_phys/atmos_param/moist_processes/moist_processes.F90
	cp $(SRCROOT)/am5_phys/atmos_param/moist_processes/moist_processes.F90 .
./longwave_types.F90: $(SRCROOT)/am5_phys/atmos_param/radiation/longwave/longwave_types.F90
	cp $(SRCROOT)/am5_phys/atmos_param/radiation/longwave/longwave_types.F90 .
./tridi.F90: $(SRCROOT)/am5_phys/atmos_param/han_edmf/tridi.F90
	cp $(SRCROOT)/am5_phys/atmos_param/han_edmf/tridi.F90 .
./zeff.F90: $(SRCROOT)/am5_phys/atmos_param/cosp/quickbeam/zeff.F90
	cp $(SRCROOT)/am5_phys/atmos_param/cosp/quickbeam/zeff.F90 .
./cosp_rttov_simulator.F90: $(SRCROOT)/am5_phys/atmos_param/cosp/cosp_rttov_simulator.F90
	cp $(SRCROOT)/am5_phys/atmos_param/cosp/cosp_rttov_simulator.F90 .
./mo_setsox.F90: $(SRCROOT)/am5_phys/atmos_shared/tracer_driver/mo_setsox.F90
	cp $(SRCROOT)/am5_phys/atmos_shared/tracer_driver/mo_setsox.F90 .
./radiative_gases_types.F90: $(SRCROOT)/am5_phys/atmos_param/radiation/radiative_gases/radiative_gases_types.F90
	cp $(SRCROOT)/am5_phys/atmos_param/radiation/radiative_gases/radiative_gases_types.F90 .
./damping_driver.F90: $(SRCROOT)/am5_phys/atmos_param/damping_driver/damping_driver.F90
	cp $(SRCROOT)/am5_phys/atmos_param/damping_driver/damping_driver.F90 .
./uw_conv.F90: $(SRCROOT)/am5_phys/atmos_param/shallow_cu/uw_conv.F90
	cp $(SRCROOT)/am5_phys/atmos_param/shallow_cu/uw_conv.F90 .
./aerosol.F90: $(SRCROOT)/am5_phys/atmos_shared/aerosol/aerosol.F90
	cp $(SRCROOT)/am5_phys/atmos_shared/aerosol/aerosol.F90 .
./atmos_nh3_tag.F90: $(SRCROOT)/am5_phys/atmos_shared/tracer_driver/atmos_nh3_tag.F90
	cp $(SRCROOT)/am5_phys/atmos_shared/tracer_driver/atmos_nh3_tag.F90 .
./radar_simulator.F90: $(SRCROOT)/am5_phys/atmos_param/cosp/quickbeam/radar_simulator.F90
	cp $(SRCROOT)/am5_phys/atmos_param/cosp/quickbeam/radar_simulator.F90 .
./cloud_chem.F90: $(SRCROOT)/am5_phys/atmos_shared/tracer_driver/cloud_chem.F90
	cp $(SRCROOT)/am5_phys/atmos_shared/tracer_driver/cloud_chem.F90 .
./longwave_utilities.F90: $(SRCROOT)/am5_phys/atmos_param/radiation/longwave/longwave_utilities.F90
	cp $(SRCROOT)/am5_phys/atmos_param/radiation/longwave/longwave_utilities.F90 .
./mo_usrrxt.F90: $(SRCROOT)/am5_phys/atmos_shared/tracer_driver/tropchem/mo_usrrxt.F90
	cp $(SRCROOT)/am5_phys/atmos_shared/tracer_driver/tropchem/mo_usrrxt.F90 .
./strat_chem_driver.F90: $(SRCROOT)/am5_phys/atmos_shared/tracer_driver/stratchem/strat_chem_driver.F90
	cp $(SRCROOT)/am5_phys/atmos_shared/tracer_driver/stratchem/strat_chem_driver.F90 .
./vert_diff_driver.F90: $(SRCROOT)/am5_phys/atmos_param/vert_diff_driver/vert_diff_driver.F90
	cp $(SRCROOT)/am5_phys/atmos_param/vert_diff_driver/vert_diff_driver.F90 .
./cosp_utils.F90: $(SRCROOT)/am5_phys/atmos_param/cosp/cosp_utils.F90
	cp $(SRCROOT)/am5_phys/atmos_param/cosp/cosp_utils.F90 .
./moz.mods.F90: $(SRCROOT)/am5_phys/atmos_shared/tracer_driver/tropchem/AM4_chem/moz.mods.F90
	cp $(SRCROOT)/am5_phys/atmos_shared/tracer_driver/tropchem/AM4_chem/moz.mods.F90 .
./cloud_spec.F90: $(SRCROOT)/am5_phys/atmos_param/radiation/clouds/cloud_spec.F90
	cp $(SRCROOT)/am5_phys/atmos_param/radiation/clouds/cloud_spec.F90 .
./aerosol_cloud.F90: $(SRCROOT)/am5_phys/atmos_param/aerosol_cloud/aerosol_cloud.F90
	cp $(SRCROOT)/am5_phys/atmos_param/aerosol_cloud/aerosol_cloud.F90 .
./lscloud_constants.F90: $(SRCROOT)/am5_phys/atmos_param/lscloud_driver/lscloud_constants.F90
	cp $(SRCROOT)/am5_phys/atmos_param/lscloud_driver/lscloud_constants.F90 .
./mo_chem_utls.F90: $(SRCROOT)/am5_phys/atmos_shared/tracer_driver/tropchem/mo_chem_utls.F90
	cp $(SRCROOT)/am5_phys/atmos_shared/tracer_driver/tropchem/mo_chem_utls.F90 .
./cloudrad_types.F90: $(SRCROOT)/am5_phys/atmos_param/radiation/clouds/cloudrad_types.F90
	cp $(SRCROOT)/am5_phys/atmos_param/radiation/clouds/cloudrad_types.F90 .
./mo_read_sim_chm.F90: $(SRCROOT)/am5_phys/atmos_shared/tracer_driver/tropchem/mo_read_sim_chm.F90
	cp $(SRCROOT)/am5_phys/atmos_shared/tracer_driver/tropchem/mo_read_sim_chm.F90 .
./atmos_global_diag.F90: $(SRCROOT)/am5_phys/atmos_shared/atmos_cmip_diag/atmos_global_diag.F90
	cp $(SRCROOT)/am5_phys/atmos_shared/atmos_cmip_diag/atmos_global_diag.F90 .
./mo_exp_slv.F90: $(SRCROOT)/am5_phys/atmos_shared/tracer_driver/tropchem/mo_exp_slv.F90
	cp $(SRCROOT)/am5_phys/atmos_shared/tracer_driver/tropchem/mo_exp_slv.F90 .
./isrpia.inc: $(SRCROOT)/am5_phys/atmos_shared/tracer_driver/tropchem/isrpia.inc
	cp $(SRCROOT)/am5_phys/atmos_shared/tracer_driver/tropchem/isrpia.inc .
./longwave_params.F90: $(SRCROOT)/am5_phys/atmos_param/radiation/longwave/longwave_params.F90
	cp $(SRCROOT)/am5_phys/atmos_param/radiation/longwave/longwave_params.F90 .
./aerosol_types.F90: $(SRCROOT)/am5_phys/atmos_shared/aerosol/aerosol_types.F90
	cp $(SRCROOT)/am5_phys/atmos_shared/aerosol/aerosol_types.F90 .
./longwave_driver.F90: $(SRCROOT)/am5_phys/atmos_param/radiation/driver/longwave_driver.F90
	cp $(SRCROOT)/am5_phys/atmos_param/radiation/driver/longwave_driver.F90 .
./xactive_bvoc.F90: $(SRCROOT)/am5_phys/atmos_shared/tracer_driver/xactive_bvoc.F90
	cp $(SRCROOT)/am5_phys/atmos_shared/tracer_driver/xactive_bvoc.F90 .
./strat_chem_model.F90: $(SRCROOT)/am5_phys/atmos_shared/tracer_driver/stratchem/strat_chem_model.F90
	cp $(SRCROOT)/am5_phys/atmos_shared/tracer_driver/stratchem/strat_chem_model.F90 .
./aerosolrad_package.F90: $(SRCROOT)/am5_phys/atmos_param/radiation/aerosols/aerosolrad_package.F90
	cp $(SRCROOT)/am5_phys/atmos_param/radiation/aerosols/aerosolrad_package.F90 .
./polysvp.F90: $(SRCROOT)/am5_phys/atmos_param/lscloud_driver/polysvp.F90
	cp $(SRCROOT)/am5_phys/atmos_param/lscloud_driver/polysvp.F90 .
./atmos_age_tracer.F90: $(SRCROOT)/am5_phys/atmos_shared/tracer_driver/atmos_age_tracer.F90
	cp $(SRCROOT)/am5_phys/atmos_shared/tracer_driver/atmos_age_tracer.F90 .
./radar_simulator_types.F90: $(SRCROOT)/am5_phys/atmos_param/cosp/quickbeam/radar_simulator_types.F90
	cp $(SRCROOT)/am5_phys/atmos_param/cosp/quickbeam/radar_simulator_types.F90 .
./mg_drag.F90: $(SRCROOT)/am5_phys/atmos_param/mg_drag/mg_drag.F90
	cp $(SRCROOT)/am5_phys/atmos_param/mg_drag/mg_drag.F90 .
./moz.subs.F90: $(SRCROOT)/am5_phys/atmos_shared/tracer_driver/tropchem/AM4_chem/moz.subs.F90
	cp $(SRCROOT)/am5_phys/atmos_shared/tracer_driver/tropchem/AM4_chem/moz.subs.F90 .
./atmos_ch4.F90: $(SRCROOT)/am5_phys/atmos_shared/tracer_driver/atmos_ch4.F90
	cp $(SRCROOT)/am5_phys/atmos_shared/tracer_driver/atmos_ch4.F90 .
./longwave_clouds.F90: $(SRCROOT)/am5_phys/atmos_param/radiation/longwave/longwave_clouds.F90
	cp $(SRCROOT)/am5_phys/atmos_param/radiation/longwave/longwave_clouds.F90 .
./detr_ice_num.F90: $(SRCROOT)/am5_phys/atmos_param/convection_driver/detr_ice_num.F90
	cp $(SRCROOT)/am5_phys/atmos_param/convection_driver/detr_ice_num.F90 .
./icarus.F90: $(SRCROOT)/am5_phys/atmos_param/cosp/icarus-scops-3.7/icarus.F90
	cp $(SRCROOT)/am5_phys/atmos_param/cosp/icarus-scops-3.7/icarus.F90 .
./radiative_gases.F90: $(SRCROOT)/am5_phys/atmos_param/radiation/radiative_gases/radiative_gases.F90
	cp $(SRCROOT)/am5_phys/atmos_param/radiation/radiative_gases/radiative_gases.F90 .
./pf_to_mr.F: $(SRCROOT)/am5_phys/atmos_param/cosp/llnl/pf_to_mr.F
	cp $(SRCROOT)/am5_phys/atmos_param/cosp/llnl/pf_to_mr.F .
./aer_ccn_act.F90: $(SRCROOT)/am5_phys/atmos_shared/tracer_driver/aer_ccn_act/aer_ccn_act.F90
	cp $(SRCROOT)/am5_phys/atmos_shared/tracer_driver/aer_ccn_act/aer_ccn_act.F90 .
./sealw99.F90: $(SRCROOT)/am5_phys/atmos_param/radiation/longwave/sealw99.F90
	cp $(SRCROOT)/am5_phys/atmos_param/radiation/longwave/sealw99.F90 .
./mfscuq.F90: $(SRCROOT)/am5_phys/atmos_param/han_edmf/mfscuq.F90
	cp $(SRCROOT)/am5_phys/atmos_param/han_edmf/mfscuq.F90 .
./atmos_tracer_driver.F90: $(SRCROOT)/am5_phys/atmos_shared/tracer_driver/atmos_tracer_driver.F90
	cp $(SRCROOT)/am5_phys/atmos_shared/tracer_driver/atmos_tracer_driver.F90 .
./MISR_simulator.F: $(SRCROOT)/am5_phys/atmos_param/cosp/MISR_simulator/MISR_simulator.F
	cp $(SRCROOT)/am5_phys/atmos_param/cosp/MISR_simulator/MISR_simulator.F .
./aerosolrad_driver.F90: $(SRCROOT)/am5_phys/atmos_param/radiation/aerosols/aerosolrad_driver.F90
	cp $(SRCROOT)/am5_phys/atmos_param/radiation/aerosols/aerosolrad_driver.F90 .
./cloudrad_package.F90: $(SRCROOT)/am5_phys/atmos_param/radiation/clouds/cloudrad_package.F90
	cp $(SRCROOT)/am5_phys/atmos_param/radiation/clouds/cloudrad_package.F90 .
./format_input.F90: $(SRCROOT)/am5_phys/atmos_param/cosp/quickbeam/format_input.F90
	cp $(SRCROOT)/am5_phys/atmos_param/cosp/quickbeam/format_input.F90 .
./atmos_nudge.F90: $(SRCROOT)/am5_phys/atmos_shared/atmos_nudge/atmos_nudge.F90
	cp $(SRCROOT)/am5_phys/atmos_shared/atmos_nudge/atmos_nudge.F90 .
./rad_output_file.F90: $(SRCROOT)/am5_phys/atmos_param/radiation/driver/rad_output_file.F90
	cp $(SRCROOT)/am5_phys/atmos_param/radiation/driver/rad_output_file.F90 .
./lw_gases_stdtf.F90: $(SRCROOT)/am5_phys/atmos_param/radiation/longwave/lw_gases_stdtf.F90
	cp $(SRCROOT)/am5_phys/atmos_param/radiation/longwave/lw_gases_stdtf.F90 .
./cosp_constants.F90: $(SRCROOT)/am5_phys/atmos_param/cosp/cosp_constants.F90
	cp $(SRCROOT)/am5_phys/atmos_param/cosp/cosp_constants.F90 .
./ice_nucl.F90: $(SRCROOT)/am5_phys/atmos_shared/tracer_driver/aer_ccn_act/ice_nucl.F90
	cp $(SRCROOT)/am5_phys/atmos_shared/tracer_driver/aer_ccn_act/ice_nucl.F90 .
./aer_in_act.F90: $(SRCROOT)/am5_phys/atmos_shared/tracer_driver/aer_ccn_act/aer_in_act.F90
	cp $(SRCROOT)/am5_phys/atmos_shared/tracer_driver/aer_ccn_act/aer_in_act.F90 .
./microphys_rad.F90: $(SRCROOT)/am5_phys/atmos_param/radiation/clouds/microphys_rad.F90
	cp $(SRCROOT)/am5_phys/atmos_param/radiation/clouds/microphys_rad.F90 .
./cloud_generator.F90: $(SRCROOT)/am5_phys/atmos_param/cloud_generator/cloud_generator.F90
	cp $(SRCROOT)/am5_phys/atmos_param/cloud_generator/cloud_generator.F90 .
./mo_fastjx.F90: $(SRCROOT)/am5_phys/atmos_shared/tracer_driver/tropchem/AM4_chem/mo_fastjx.F90
	cp $(SRCROOT)/am5_phys/atmos_shared/tracer_driver/tropchem/AM4_chem/mo_fastjx.F90 .
./mo_chemini.F90: $(SRCROOT)/am5_phys/atmos_shared/tracer_driver/tropchem/mo_chemini.F90
	cp $(SRCROOT)/am5_phys/atmos_shared/tracer_driver/tropchem/mo_chemini.F90 .
./math_lib.F90: $(SRCROOT)/am5_phys/atmos_param/cosp/quickbeam/math_lib.F90
	cp $(SRCROOT)/am5_phys/atmos_param/cosp/quickbeam/math_lib.F90 .
./esfsw_driver.F90: $(SRCROOT)/am5_phys/atmos_param/radiation/shortwave/esfsw_driver.F90
	cp $(SRCROOT)/am5_phys/atmos_param/radiation/shortwave/esfsw_driver.F90 .
./mfpbltq.F90: $(SRCROOT)/am5_phys/atmos_param/han_edmf/mfpbltq.F90
	cp $(SRCROOT)/am5_phys/atmos_param/han_edmf/mfpbltq.F90 .
./cosp_lidar.F90: $(SRCROOT)/am5_phys/atmos_param/cosp/cosp_lidar.F90
	cp $(SRCROOT)/am5_phys/atmos_param/cosp/cosp_lidar.F90 .
./atmos_co2.F90: $(SRCROOT)/am5_phys/atmos_shared/tracer_driver/atmos_co2.F90
	cp $(SRCROOT)/am5_phys/atmos_shared/tracer_driver/atmos_co2.F90 .
./fms_platform.h: $(SRCROOT)FMS/include/fms_platform.h
	cp $(SRCROOT)FMS/include/fms_platform.h .
./vert_turb_driver.F90: $(SRCROOT)/am5_phys/atmos_param/vert_turb_driver/vert_turb_driver.F90
	cp $(SRCROOT)/am5_phys/atmos_param/vert_turb_driver/vert_turb_driver.F90 .
./cosp_defs.H: $(SRCROOT)/am5_phys/atmos_param/cosp/cosp_defs.H
	cp $(SRCROOT)/am5_phys/atmos_param/cosp/cosp_defs.H .
./strat_chem_utilities.F90: $(SRCROOT)/am5_phys/atmos_shared/tracer_driver/tropchem/strat_chem_utilities.F90
	cp $(SRCROOT)/am5_phys/atmos_shared/tracer_driver/tropchem/strat_chem_utilities.F90 .
./atmos_ch3i.F90: $(SRCROOT)/am5_phys/atmos_shared/tracer_driver/atmos_ch3i.F90
	cp $(SRCROOT)/am5_phys/atmos_shared/tracer_driver/atmos_ch3i.F90 .
./lidar_simulator.F90: $(SRCROOT)/am5_phys/atmos_param/cosp/actsim/lidar_simulator.F90
	cp $(SRCROOT)/am5_phys/atmos_param/cosp/actsim/lidar_simulator.F90 .
./lmd_ipsl_stats.F90: $(SRCROOT)/am5_phys/atmos_param/cosp/actsim/lmd_ipsl_stats.F90
	cp $(SRCROOT)/am5_phys/atmos_param/cosp/actsim/lmd_ipsl_stats.F90 .
./conv_plumes.F90: $(SRCROOT)/am5_phys/atmos_param/shallow_cu/conv_plumes.F90
	cp $(SRCROOT)/am5_phys/atmos_param/shallow_cu/conv_plumes.F90 .
./conv_plumes_k.F90: $(SRCROOT)/am5_phys/atmos_param/shallow_cu/conv_plumes_k.F90
	cp $(SRCROOT)/am5_phys/atmos_param/shallow_cu/conv_plumes_k.F90 .
./scops.F: $(SRCROOT)/am5_phys/atmos_param/cosp/icarus-scops-3.7/scops.F
	cp $(SRCROOT)/am5_phys/atmos_param/cosp/icarus-scops-3.7/scops.F .
./shortwave_types.F90: $(SRCROOT)/am5_phys/atmos_param/radiation/shortwave/shortwave_types.F90
	cp $(SRCROOT)/am5_phys/atmos_param/radiation/shortwave/shortwave_types.F90 .
./modis_simulator.F90: $(SRCROOT)/am5_phys/atmos_param/cosp/MODIS_simulator/modis_simulator.F90
	cp $(SRCROOT)/am5_phys/atmos_param/cosp/MODIS_simulator/modis_simulator.F90 .
./cu_mo_trans.F90: $(SRCROOT)/am5_phys/atmos_param/cu_mo_trans/cu_mo_trans.F90
	cp $(SRCROOT)/am5_phys/atmos_param/cu_mo_trans/cu_mo_trans.F90 .
./lscloud_debug.F90: $(SRCROOT)/am5_phys/atmos_param/lscloud_driver/lscloud_debug.F90
	cp $(SRCROOT)/am5_phys/atmos_param/lscloud_driver/lscloud_debug.F90 .
./strat_clouds_W.F90: $(SRCROOT)/am5_phys/atmos_param/radiation/clouds/strat_clouds_W.F90
	cp $(SRCROOT)/am5_phys/atmos_param/radiation/clouds/strat_clouds_W.F90 .
./uw_clouds_W.F90: $(SRCROOT)/am5_phys/atmos_param/radiation/clouds/uw_clouds_W.F90
	cp $(SRCROOT)/am5_phys/atmos_param/radiation/clouds/uw_clouds_W.F90 .
./tropchem_types_mod.F90: $(SRCROOT)/am5_phys/atmos_shared/tracer_driver/tropchem/tropchem_types_mod.F90
	cp $(SRCROOT)/am5_phys/atmos_shared/tracer_driver/tropchem/tropchem_types_mod.F90 .
./cg_drag.F90: $(SRCROOT)/am5_phys/atmos_param/cg_drag/cg_drag.F90
	cp $(SRCROOT)/am5_phys/atmos_param/cg_drag/cg_drag.F90 .
./atmos_lib.F90: $(SRCROOT)/am5_phys/atmos_param/cosp/quickbeam/atmos_lib.F90
	cp $(SRCROOT)/am5_phys/atmos_param/cosp/quickbeam/atmos_lib.F90 .
./atmos_radon.F90: $(SRCROOT)/am5_phys/atmos_shared/tracer_driver/atmos_radon.F90
	cp $(SRCROOT)/am5_phys/atmos_shared/tracer_driver/atmos_radon.F90 .
./micro_mg2_utils.F90: $(SRCROOT)/am5_phys/atmos_param/microphysics/micro_mg2_utils.F90
	cp $(SRCROOT)/am5_phys/atmos_param/microphysics/micro_mg2_utils.F90 .
./longwave_fluxes.F90: $(SRCROOT)/am5_phys/atmos_param/radiation/longwave/longwave_fluxes.F90
	cp $(SRCROOT)/am5_phys/atmos_param/radiation/longwave/longwave_fluxes.F90 .
./gw_beres.F90: $(SRCROOT)/am5_phys/atmos_param/gw_beres/gw_beres.F90
	cp $(SRCROOT)/am5_phys/atmos_param/gw_beres/gw_beres.F90 .
./moist_processes_utils.F90: $(SRCROOT)/am5_phys/atmos_param/moist_processes/moist_processes_utils.F90
	cp $(SRCROOT)/am5_phys/atmos_param/moist_processes/moist_processes_utils.F90 .
./ls_cloud_macrophysics.F90: $(SRCROOT)/am5_phys/atmos_param/macrophysics/ls_cloud_macrophysics.F90
	cp $(SRCROOT)/am5_phys/atmos_param/macrophysics/ls_cloud_macrophysics.F90 .
./ls_cloud_microphysics.F90: $(SRCROOT)/am5_phys/atmos_param/microphysics/ls_cloud_microphysics.F90
	cp $(SRCROOT)/am5_phys/atmos_param/microphysics/ls_cloud_microphysics.F90 .
./entrain.F90: $(SRCROOT)/am5_phys/atmos_param/entrain/entrain.F90
	cp $(SRCROOT)/am5_phys/atmos_param/entrain/entrain.F90 .
./dsd.F90: $(SRCROOT)/am5_phys/atmos_param/cosp/quickbeam/dsd.F90
	cp $(SRCROOT)/am5_phys/atmos_param/cosp/quickbeam/dsd.F90 .
./cosp_misr_simulator.F90: $(SRCROOT)/am5_phys/atmos_param/cosp/cosp_misr_simulator.F90
	cp $(SRCROOT)/am5_phys/atmos_param/cosp/cosp_misr_simulator.F90 .
./physics_types.F90: $(SRCROOT)/am5_phys/atmos_param/physics_driver/physics_types.F90
	cp $(SRCROOT)/am5_phys/atmos_param/physics_driver/physics_types.F90 .
./mo_fphoto.F90: $(SRCROOT)/am5_phys/atmos_shared/tracer_driver/tropchem/AM4_chem/mo_fphoto.F90
	cp $(SRCROOT)/am5_phys/atmos_shared/tracer_driver/tropchem/AM4_chem/mo_fphoto.F90 .
./gas_tf.F90: $(SRCROOT)/am5_phys/atmos_param/radiation/longwave/gas_tf.F90
	cp $(SRCROOT)/am5_phys/atmos_param/radiation/longwave/gas_tf.F90 .
./atmos_soa.F90: $(SRCROOT)/am5_phys/atmos_shared/tracer_driver/atmos_soa.F90
	cp $(SRCROOT)/am5_phys/atmos_shared/tracer_driver/atmos_soa.F90 .
./aerosol_thermodynamics.F90: $(SRCROOT)/am5_phys/atmos_shared/tracer_driver/tropchem/aerosol_thermodynamics.F90
	cp $(SRCROOT)/am5_phys/atmos_shared/tracer_driver/tropchem/aerosol_thermodynamics.F90 .
./congvec.H: $(SRCROOT)/am5_phys/atmos_param/cosp/icarus-scops-3.7/congvec.H
	cp $(SRCROOT)/am5_phys/atmos_param/cosp/icarus-scops-3.7/congvec.H .
./physics_radiation_exch.F90: $(SRCROOT)/am5_phys/atmos_param/physics_radiation_exch/physics_radiation_exch.F90
	cp $(SRCROOT)/am5_phys/atmos_param/physics_radiation_exch/physics_radiation_exch.F90 .
./aerosol_params.F90: $(SRCROOT)/am5_phys/atmos_shared/tracer_driver/aer_ccn_act/aerosol_params.F90
	cp $(SRCROOT)/am5_phys/atmos_shared/tracer_driver/aer_ccn_act/aerosol_params.F90 .
./calc_Re.F90: $(SRCROOT)/am5_phys/atmos_param/cosp/quickbeam/calc_Re.F90
	cp $(SRCROOT)/am5_phys/atmos_param/cosp/quickbeam/calc_Re.F90 .
./radar_simulator_init.F90: $(SRCROOT)/am5_phys/atmos_param/cosp/quickbeam/radar_simulator_init.F90
	cp $(SRCROOT)/am5_phys/atmos_param/cosp/quickbeam/radar_simulator_init.F90 .
./physics_driver.F90: $(SRCROOT)/am5_phys/atmos_param/physics_driver/physics_driver.F90
	cp $(SRCROOT)/am5_phys/atmos_param/physics_driver/physics_driver.F90 .
./cosp_diagnostics.F90: $(SRCROOT)/am5_phys/atmos_param/cosp/cosp_diagnostics.F90
	cp $(SRCROOT)/am5_phys/atmos_param/cosp/cosp_diagnostics.F90 .
./cosp_types.F90: $(SRCROOT)/am5_phys/atmos_param/cosp/cosp_types.F90
	cp $(SRCROOT)/am5_phys/atmos_param/cosp/cosp_types.F90 .
./atmos_cmip_diag.F90: $(SRCROOT)/am5_phys/atmos_shared/atmos_cmip_diag/atmos_cmip_diag.F90
	cp $(SRCROOT)/am5_phys/atmos_shared/atmos_cmip_diag/atmos_cmip_diag.F90 .
./radiation_driver_diag.F90: $(SRCROOT)/am5_phys/atmos_param/radiation/driver/radiation_driver_diag.F90
	cp $(SRCROOT)/am5_phys/atmos_param/radiation/driver/radiation_driver_diag.F90 .
./tropchem_driver.F90: $(SRCROOT)/am5_phys/atmos_shared/tracer_driver/tropchem_driver.F90
	cp $(SRCROOT)/am5_phys/atmos_shared/tracer_driver/tropchem_driver.F90 .
./cloudrad_diagnostics.F90: $(SRCROOT)/am5_phys/atmos_param/radiation/clouds/cloudrad_diagnostics.F90
	cp $(SRCROOT)/am5_phys/atmos_param/radiation/clouds/cloudrad_diagnostics.F90 .
./mo_photo.F90: $(SRCROOT)/am5_phys/atmos_shared/tracer_driver/tropchem/mo_photo.F90
	cp $(SRCROOT)/am5_phys/atmos_shared/tracer_driver/tropchem/mo_photo.F90 .
./mo_rodas_slv.F90: $(SRCROOT)/am5_phys/atmos_shared/tracer_driver/tropchem/mo_rodas_slv.F90
	cp $(SRCROOT)/am5_phys/atmos_shared/tracer_driver/tropchem/mo_rodas_slv.F90 .
./gamma_mg.F90: $(SRCROOT)/am5_phys/atmos_param/microphysics/gamma_mg.F90
	cp $(SRCROOT)/am5_phys/atmos_param/microphysics/gamma_mg.F90 .
./longwave_tables.F90: $(SRCROOT)/am5_phys/atmos_param/radiation/longwave/longwave_tables.F90
	cp $(SRCROOT)/am5_phys/atmos_param/radiation/longwave/longwave_tables.F90 .
./atmos_tropopause.F90: $(SRCROOT)/am5_phys/atmos_shared/tracer_driver/atmos_tropopause.F90
	cp $(SRCROOT)/am5_phys/atmos_shared/tracer_driver/atmos_tropopause.F90 .
./convection_utilities.F90: $(SRCROOT)/am5_phys/atmos_param/convection_driver/convection_utilities.F90
	cp $(SRCROOT)/am5_phys/atmos_param/convection_driver/convection_utilities.F90 .
SRC = $(SRCROOT)/am5_phys/atmos_param/radiation/radiative_gases/radiative_gases_types.F90 $(SRCROOT)/am5_phys/atmos_param/radiation/driver/solar_data_driver.F90 $(SRCROOT)/am5_phys/atmos_shared/tracer_driver/atmos_tracer_utilities.F90 $(SRCROOT)/am5_phys/atmos_param/radiation/longwave/optical_path.F90 $(SRCROOT)/am5_phys/atmos_shared/tracer_driver/tropchem/m_tracname.F90 $(SRCROOT)/am5_phys/atmos_param/radiation/longwave/longwave_tables.F90 $(SRCROOT)/am5_phys/atmos_param/cosp/quickbeam/radar_simulator_init.F90 $(SRCROOT)/am5_phys/atmos_param/microphysics/micro_mg2_utils.F90 $(SRCROOT)/am5_phys/atmos_param/gw_beres/gw_common_utils.F90 $(SRCROOT)/am5_phys/atmos_param/cosp/quickbeam/radar_simulator_types.F90 $(SRCROOT)/am5_phys/atmos_shared/tracer_driver/tropchem/AM4_chem/moz.mat.F90 $(SRCROOT)/am5_phys/atmos_param/cg_drag/cg_drag.F90 $(SRCROOT)/am5_phys/atmos_shared/tracer_driver/atmos_age_tracer.F90 $(SRCROOT)/am5_phys/atmos_param/cosp/quickbeam/mrgrnk.F90 $(SRCROOT)/am5_phys/atmos_param/cosp/cosp_rttov_simulator.F90 $(SRCROOT)/am5_phys/atmos_shared/tracer_driver/atmos_regional_tracer_driver.F90 $(SRCROOT)/am5_phys/atmos_shared/tracer_driver/tropchem/isoropiaIIcode.F $(SRCROOT)/am5_phys/atmos_shared/atmos_cmip_diag/atmos_cmip_diag.F90 $(SRCROOT)/am5_phys/atmos_shared/tracer_driver/atmos_convection_tracer.F90 $(SRCROOT)/am5_phys/atmos_param/vert_diff_driver/vert_diff_driver.F90 $(SRCROOT)/am5_phys/atmos_shared/tracer_driver/atmos_sulfur_hex.F90 $(SRCROOT)/am5_phys/atmos_param/cosp/icarus-scops-3.7/icarus.F90 $(SRCROOT)/am5_phys/atmos_param/radiation/aerosols/aerosolrad_driver.F90 $(SRCROOT)/am5_phys/atmos_param/radiation/clouds/cloudrad_diagnostics.F90 $(SRCROOT)/am5_phys/atmos_param/radiation/driver/rad_output_file.F90 $(SRCROOT)/am5_phys/atmos_param/lscloud_driver/lscloud_netcdf.F90 $(SRCROOT)/am5_phys/atmos_shared/tracer_driver/tropchem/mo_imp_slv.F90 $(SRCROOT)/am5_phys/atmos_shared/tracer_driver/tropchem/strat_chem_utilities.F90 $(SRCROOT)/am5_phys/atmos_param/gw_beres/gw_beres.F90 $(SRCROOT)/am5_phys/atmos_param/cosp/llnl/cosp_radar.F90 $(SRCROOT)/am5_phys/atmos_param/radiation/driver/radiation_driver_diag.F90 $(SRCROOT)/am5_phys/atmos_param/cosp/icarus-scops-3.7/scops.F $(SRCROOT)/am5_phys/atmos_param/radiation/shortwave/esfsw_parameters.F90 $(SRCROOT)/am5_phys/atmos_param/radiation/clouds/microphys_rad.F90 $(SRCROOT)/am5_phys/atmos_param/moist_processes/moist_processes.F90 $(SRCROOT)/am5_phys/atmos_shared/tracer_driver/atmos_radon.F90 $(SRCROOT)/am5_phys/atmos_param/cosp/quickbeam/dsd.F90 $(SRCROOT)/am5_phys/atmos_param/radiation/longwave/longwave_utilities.F90 $(SRCROOT)/am5_phys/atmos_shared/tracer_driver/cloud_chem.F90 $(SRCROOT)/am5_phys/atmos_shared/tracer_driver/tropchem/mo_chemdr.F90 $(SRCROOT)/am5_phys/atmos_shared/tracer_driver/mo_setsox.F90 $(SRCROOT)/am5_phys/atmos_param/cloud_generator/betaDistribution.F90 $(SRCROOT)/am5_phys/atmos_shared/tracer_driver/atmos_ch3i.F90 $(SRCROOT)/am5_phys/atmos_param/radiation/driver/radiation_types.F90 $(SRCROOT)/am5_phys/atmos_param/moist_processes/moist_processes_utils.F90 $(SRCROOT)/am5_phys/atmos_shared/tracer_driver/aer_ccn_act/aer_ccn_act_k.F90 $(SRCROOT)/am5_phys/atmos_param/cosp/cosp_types.F90 $(SRCROOT)/am5_phys/atmos_param/shallow_cu/conv_plumes.F90 $(SRCROOT)/am5_phys/atmos_param/convection_driver/convection_driver.F90 $(SRCROOT)/am5_phys/atmos_shared/tracer_driver/atmos_ch4.F90 $(SRCROOT)/am5_phys/atmos_param/mg_drag/mg_drag.F90 $(SRCROOT)/am5_phys/atmos_param/cosp/MODIS_simulator/modis_simulator.F90 $(SRCROOT)/am5_phys/atmos_param/cosp/actsim/lidar_simulator.F90 $(SRCROOT)/am5_phys/atmos_param/cosp/actsim/lmd_ipsl_stats.F90 $(SRCROOT)/am5_phys/atmos_param/radiation/driver/shortwave_driver.F90 $(SRCROOT)/am5_phys/atmos_param/cosp/quickbeam/format_input.F90 $(SRCROOT)/am5_phys/atmos_shared/tracer_driver/tropchem/mo_hook.F90 $(SRCROOT)/am5_phys/atmos_param/cosp/quickbeam/atmos_lib.F90 $(SRCROOT)/am5_phys/atmos_shared/tracer_driver/atmos_tracer_driver.F90 $(SRCROOT)/am5_phys/atmos_param/shallow_cu/conv_utilities.F90 $(SRCROOT)/am5_phys/atmos_shared/tracer_driver/aer_ccn_act/aerosol_params.F90 $(SRCROOT)/am5_phys/atmos_shared/tracer_driver/tropchem/mo_read_sim_chm.F90 $(SRCROOT)/am5_phys/atmos_param/cosp/cosp_diagnostics.F90 $(SRCROOT)/am5_phys/atmos_shared/tracer_driver/tropchem/mo_usrrxt.F90 $(SRCROOT)/am5_phys/atmos_param/han_edmf/hanedmf.F90 $(SRCROOT)/am5_phys/atmos_shared/tracer_driver/aer_ccn_act/ice_nucl.F90 $(SRCROOT)/am5_phys/atmos_param/shallow_cu/conv_closures.F90 $(SRCROOT)/am5_phys/atmos_param/radiation/shortwave/esfsw_utilities.F90 $(SRCROOT)/am5_phys/atmos_shared/tracer_driver/tropchem/mo_jpl.F90 $(SRCROOT)/am5_phys/atmos_param/radiation/shortwave/esfsw_driver.F90 $(SRCROOT)/am5_phys/atmos_param/cosp/cosp_modis_simulator.F90 $(SRCROOT)/am5_phys/atmos_shared/tracer_driver/atmos_soa.F90 $(SRCROOT)/am5_phys/atmos_param/lscloud_driver/lscloud_driver.F90 $(SRCROOT)/am5_phys/atmos_shared/tracer_driver/tropchem/mo_exp_slv.F90 $(SRCROOT)/am5_phys/atmos_param/radiation/aerosols/aerosolrad_types.F90 $(SRCROOT)/am5_phys/atmos_param/physics_driver/physics_driver.F90 $(SRCROOT)/am5_phys/atmos_param/lscloud_driver/polysvp.F90 $(SRCROOT)/am5_phys/atmos_param/cosp/cosp.F90 $(SRCROOT)/am5_phys/atmos_shared/tracer_driver/tropchem/mo_setinv.F90 $(SRCROOT)/am5_phys/atmos_param/vert_diff/vert_diff.F90 $(SRCROOT)/am5_phys/atmos_shared/tracer_driver/aer_ccn_act/aer_ccn_act.F90 $(SRCROOT)/am5_phys/atmos_param/microphysics/ls_cloud_microphysics.F90 $(SRCROOT)/am5_phys/atmos_param/radiation/aerosols/aerosolrad_package.F90 $(SRCROOT)/am5_phys/atmos_shared/tracer_driver/tropchem/AM4_chem/moz.mods.F90 $(SRCROOT)/am5_phys/atmos_shared/tracer_driver/tropchem/AM4_chem/moz.subs.F90 $(SRCROOT)/am5_phys/atmos_shared/aerosol/aerosol_types.F90 $(SRCROOT)/am5_phys/atmos_param/radiation/longwave/longwave_fluxes.F90 $(SRCROOT)/am5_phys/atmos_shared/vert_advection/vert_advection.F90 $(SRCROOT)/am5_phys/atmos_shared/tracer_driver/tropchem/AM4_chem/mo_fastjx.F90 $(SRCROOT)/am5_phys/atmos_shared/tracer_driver/stratchem/strat_chem_model.F90 $(SRCROOT)/am5_phys/atmos_param/cosp/llnl/llnl_stats.F90 $(SRCROOT)/am5_phys/atmos_param/cosp/cosp_utils.F90 $(SRCROOT)/am5_phys/atmos_param/shallow_cu/conv_plumes_k.F90 $(SRCROOT)/am5_phys/atmos_param/cosp/quickbeam/array_lib.F90 $(SRCROOT)/am5_phys/atmos_param/vert_turb_driver/vert_turb_driver.F90 $(SRCROOT)/am5_phys/atmos_param/cosp/cosp_misr_simulator.F90 $(SRCROOT)/am5_phys/atmos_param/radiation/longwave/gas_tf.F90 $(SRCROOT)/am5_phys/atmos_param/cosp/llnl/pf_to_mr.F $(SRCROOT)/am5_phys/atmos_param/radiation/clouds/cloudrad_types.F90 $(SRCROOT)/am5_phys/atmos_shared/tracer_driver/xactive_bvoc.F90 $(SRCROOT)/am5_phys/atmos_param/radiation/driver/longwave_driver.F90 $(SRCROOT)/am5_phys/atmos_param/han_edmf/satmedmfvdifq.F90 $(SRCROOT)/am5_phys/atmos_param/shallow_cu/conv_utilities_k.F90 $(SRCROOT)/am5_phys/atmos_shared/tracer_driver/atmos_tropopause.F90 $(SRCROOT)/am5_phys/atmos_param/lscloud_driver/check_nan.F90 $(SRCROOT)/am5_phys/atmos_param/convection_driver/detr_ice_num.F90 $(SRCROOT)/am5_phys/atmos_param/han_edmf/tridi.F90 $(SRCROOT)/am5_phys/atmos_shared/tracer_driver/tropchem/mo_chemini.F90 $(SRCROOT)/am5_phys/atmos_shared/tracer_driver/stratchem/strat_chem_driver.F90 $(SRCROOT)/am5_phys/atmos_param/radiation/shortwave/shortwave_types.F90 $(SRCROOT)/am5_phys/atmos_shared/aerosol/aerosol.F90 $(SRCROOT)/am5_phys/atmos_param/radiation/radiative_gases/ozone.F90 $(SRCROOT)/am5_phys/atmos_param/radiation/driver/radiation_driver.F90 $(SRCROOT)/am5_phys/atmos_param/han_edmf/mfpbltq.F90 $(SRCROOT)/am5_phys/atmos_param/shallow_cu/uw_conv.F90 $(SRCROOT)/am5_phys/atmos_param/radiation/clouds/isccp_clouds.F90 $(SRCROOT)/am5_phys/atmos_param/macrophysics/ls_cloud_macrophysics.F90 $(SRCROOT)/am5_phys/atmos_param/lscloud_driver/lscloud_constants.F90 $(SRCROOT)/am5_phys/atmos_param/han_edmf/mfscuq.F90 $(SRCROOT)/am5_phys/atmos_shared/tracer_driver/atmos_co2.F90 $(SRCROOT)/am5_phys/atmos_shared/atmos_cmip_diag/atmos_global_diag.F90 $(SRCROOT)/am5_phys/atmos_param/radiation/longwave/longwave_params.F90 $(SRCROOT)/am5_phys/atmos_param/cosp/icarus-scops-3.7/isccp_cloud_types.F $(SRCROOT)/am5_phys/atmos_param/physics_radiation_exch/physics_radiation_exch.F90 $(SRCROOT)/am5_phys/atmos_param/convection_driver/convection_utilities.F90 $(SRCROOT)/am5_phys/atmos_shared/tracer_driver/atmos_carbon_aerosol.F90 $(SRCROOT)/am5_phys/atmos_shared/tracer_driver/tropchem/mo_chem_utls.F90 $(SRCROOT)/am5_phys/atmos_param/radiation/clouds/cloud_spec.F90 $(SRCROOT)/am5_phys/atmos_shared/tracer_driver/tropchem/mo_photo.F90 $(SRCROOT)/am5_phys/atmos_shared/atmos_nudge/atmos_nudge.F90 $(SRCROOT)/am5_phys/atmos_param/radiation/longwave/lw_gases_stdtf.F90 $(SRCROOT)/am5_phys/atmos_param/cosp/quickbeam/scale_LUTs_io.F90 $(SRCROOT)/am5_phys/atmos_param/lscloud_driver/lscloud_debug.F90 $(SRCROOT)/am5_phys/atmos_param/cosp/quickbeam/gases.F90 $(SRCROOT)/am5_phys/atmos_param/cosp/llnl/prec_scops.F $(SRCROOT)/am5_phys/atmos_param/microphysics/rotstayn_klein_mp.F90 $(SRCROOT)/am5_phys/atmos_param/radiation/longwave/sealw99.F90 $(SRCROOT)/am5_phys/atmos_shared/tracer_driver/atmos_sea_salt.F90 $(SRCROOT)/am5_phys/atmos_param/microphysics/micro_mg2.F90 $(SRCROOT)/am5_phys/atmos_param/cosp/quickbeam/zeff.F90 $(SRCROOT)/am5_phys/atmos_param/radiation/clouds/cloudrad_package.F90 $(SRCROOT)/am5_phys/atmos_param/macrophysics/tiedtke_macro.F90 $(SRCROOT)/am5_phys/atmos_shared/tracer_driver/tropchem_driver.F90 $(SRCROOT)/am5_phys/atmos_param/cosp/cosp_stats.F90 $(SRCROOT)/am5_phys/atmos_param/cosp/MISR_simulator/MISR_simulator.F $(SRCROOT)/am5_phys/atmos_param/cosp/cosp_simulator.F90 $(SRCROOT)/am5_phys/atmos_param/aerosol_cloud/aerosol_cloud.F90 $(SRCROOT)/am5_phys/atmos_param/stable_bl_turb/stable_bl_turb.F90 $(SRCROOT)/am5_phys/atmos_param/lscale_cond/lscale_cond.F90 $(SRCROOT)/am5_phys/atmos_param/cosp/quickbeam/calc_Re.F90 $(SRCROOT)/am5_phys/atmos_shared/tracer_driver/atmos_dust.F90 $(SRCROOT)/am5_phys/atmos_param/radiation/shortwave/esfsw_bands.F90 $(SRCROOT)/am5_phys/atmos_param/cosp/quickbeam/radar_simulator.F90 $(SRCROOT)/am5_phys/atmos_param/radiation/radiative_gases/radiative_gases.F90 $(SRCROOT)/am5_phys/atmos_param/cosp/cosp_isccp_simulator.F90 $(SRCROOT)/am5_phys/atmos_param/cloud_generator/cloud_generator.F90 $(SRCROOT)/am5_phys/atmos_shared/tracer_driver/tropchem/mo_rodas_slv.F90 $(SRCROOT)/am5_phys/atmos_param/cosp/cosp_driver.F90 $(SRCROOT)/am5_phys/atmos_shared/tracer_driver/tropchem/AM4_chem/mo_fphoto.F90 $(SRCROOT)/am5_phys/atmos_param/radiation/util/radiation_diag.F90 $(SRCROOT)/am5_phys/atmos_shared/tracer_driver/atmos_ozone_tracer.F90 $(SRCROOT)/am5_phys/atmos_param/microphysics/gamma_mg.F90 $(SRCROOT)/am5_phys/atmos_param/radiation/longwave/longwave_types.F90 $(SRCROOT)/am5_phys/atmos_param/shallow_cu/deep_conv.F90 $(SRCROOT)/am5_phys/atmos_param/damping_driver/damping_driver.F90 $(SRCROOT)/am5_phys/atmos_param/radiation/clouds/cloudrad_driver.F90 $(SRCROOT)/am5_phys/atmos_shared/tracer_driver/tropchem/tropchem_types_mod.F90 $(SRCROOT)/am5_phys/atmos_param/cosp/cosp_io.F90 $(SRCROOT)/am5_phys/atmos_param/cloud_rad/cloud_rad.F90 $(SRCROOT)/am5_phys/atmos_shared/tracer_driver/aer_ccn_act/aer_in_act.F90 $(SRCROOT)/am5_phys/atmos_param/radiation/clouds/get_random_number_stream.F90 $(SRCROOT)/am5_phys/atmos_shared/tracer_driver/atmos_sulfate.F90 $(SRCROOT)/am5_phys/atmos_param/radiation/longwave/longwave_clouds.F90 $(SRCROOT)/am5_phys/atmos_param/lscloud_driver/lscloud_types.F90 $(SRCROOT)/am5_phys/atmos_param/cosp/cosp_lidar.F90 $(SRCROOT)/am5_phys/atmos_param/cu_mo_trans/cu_mo_trans.F90 $(SRCROOT)/am5_phys/atmos_param/radiation/clouds/strat_clouds_W.F90 $(SRCROOT)/am5_phys/atmos_param/cosp/quickbeam/optics_lib.F90 $(SRCROOT)/am5_phys/atmos_shared/tracer_driver/tropchem/aerosol_thermodynamics.F90 $(SRCROOT)/am5_phys/atmos_param/cosp/cosp_constants.F90 $(SRCROOT)/am5_phys/atmos_param/cosp/quickbeam/math_lib.F90 $(SRCROOT)/am5_phys/atmos_param/radiation/driver/radiation_driver_types.F90 $(SRCROOT)/am5_phys/atmos_param/shallow_cu/wet_deposition_0D.F90 $(SRCROOT)/am5_phys/atmos_shared/tracer_driver/atmos_nh3_tag.F90 $(SRCROOT)/am5_phys/atmos_param/topo_drag/topo_drag.F90 $(SRCROOT)/am5_phys/atmos_param/radiation/clouds/uw_clouds_W.F90 $(SRCROOT)/am5_phys/atmos_param/physics_driver/physics_types.F90 $(SRCROOT)/am5_phys/atmos_param/entrain/entrain.F90 $(SRCROOT)/am5_phys/atmos_param/cosp/icarus-scops-3.7/congvec.H $(SRCROOT)/am5_phys/atmos_shared/tracer_driver/tropchem/isrpia.inc $(SRCROOT)FMS/include/fms_platform.h $(SRCROOT)/am5_phys/atmos_param/cosp/cosp_defs.H $(SRCROOT)/am5_phys/atmos_shared/atmos_cmip_diag/atmos_cmip_interp.inc
OBJ = radiative_gases_types.o solar_data_driver.o atmos_tracer_utilities.o optical_path.o m_tracname.o longwave_tables.o radar_simulator_init.o micro_mg2_utils.o gw_common_utils.o radar_simulator_types.o moz.mat.o cg_drag.o atmos_age_tracer.o mrgrnk.o cosp_rttov_simulator.o atmos_regional_tracer_driver.o isoropiaIIcode.o atmos_cmip_diag.o atmos_convection_tracer.o vert_diff_driver.o atmos_sulfur_hex.o icarus.o aerosolrad_driver.o cloudrad_diagnostics.o rad_output_file.o lscloud_netcdf.o mo_imp_slv.o strat_chem_utilities.o gw_beres.o cosp_radar.o radiation_driver_diag.o scops.o esfsw_parameters.o microphys_rad.o moist_processes.o atmos_radon.o dsd.o longwave_utilities.o cloud_chem.o mo_chemdr.o mo_setsox.o betaDistribution.o atmos_ch3i.o radiation_types.o moist_processes_utils.o aer_ccn_act_k.o cosp_types.o conv_plumes.o convection_driver.o atmos_ch4.o mg_drag.o modis_simulator.o lidar_simulator.o lmd_ipsl_stats.o shortwave_driver.o format_input.o mo_hook.o atmos_lib.o atmos_tracer_driver.o conv_utilities.o aerosol_params.o mo_read_sim_chm.o cosp_diagnostics.o mo_usrrxt.o hanedmf.o ice_nucl.o conv_closures.o esfsw_utilities.o mo_jpl.o esfsw_driver.o cosp_modis_simulator.o atmos_soa.o lscloud_driver.o mo_exp_slv.o aerosolrad_types.o physics_driver.o polysvp.o cosp.o mo_setinv.o vert_diff.o aer_ccn_act.o ls_cloud_microphysics.o aerosolrad_package.o moz.mods.o moz.subs.o aerosol_types.o longwave_fluxes.o vert_advection.o mo_fastjx.o strat_chem_model.o llnl_stats.o cosp_utils.o conv_plumes_k.o array_lib.o vert_turb_driver.o cosp_misr_simulator.o gas_tf.o pf_to_mr.o cloudrad_types.o xactive_bvoc.o longwave_driver.o satmedmfvdifq.o conv_utilities_k.o atmos_tropopause.o check_nan.o detr_ice_num.o tridi.o mo_chemini.o strat_chem_driver.o shortwave_types.o aerosol.o ozone.o radiation_driver.o mfpbltq.o uw_conv.o isccp_clouds.o ls_cloud_macrophysics.o lscloud_constants.o mfscuq.o atmos_co2.o atmos_global_diag.o longwave_params.o isccp_cloud_types.o physics_radiation_exch.o convection_utilities.o atmos_carbon_aerosol.o mo_chem_utls.o cloud_spec.o mo_photo.o atmos_nudge.o lw_gases_stdtf.o scale_LUTs_io.o lscloud_debug.o gases.o prec_scops.o rotstayn_klein_mp.o sealw99.o atmos_sea_salt.o micro_mg2.o zeff.o cloudrad_package.o tiedtke_macro.o tropchem_driver.o cosp_stats.o MISR_simulator.o cosp_simulator.o aerosol_cloud.o stable_bl_turb.o lscale_cond.o calc_Re.o atmos_dust.o esfsw_bands.o radar_simulator.o radiative_gases.o cosp_isccp_simulator.o cloud_generator.o mo_rodas_slv.o cosp_driver.o mo_fphoto.o radiation_diag.o atmos_ozone_tracer.o gamma_mg.o longwave_types.o deep_conv.o damping_driver.o cloudrad_driver.o tropchem_types_mod.o cosp_io.o cloud_rad.o aer_in_act.o get_random_number_stream.o atmos_sulfate.o longwave_clouds.o lscloud_types.o cosp_lidar.o cu_mo_trans.o strat_clouds_W.o optics_lib.o aerosol_thermodynamics.o cosp_constants.o math_lib.o radiation_driver_types.o wet_deposition_0D.o atmos_nh3_tag.o topo_drag.o uw_clouds_W.o physics_types.o entrain.o
OFF = $(SRCROOT)/am5_phys/atmos_shared/tracer_driver/atmos_dust.F90 $(SRCROOT)/am5_phys/atmos_param/cosp/cosp_stats.F90 $(SRCROOT)/am5_phys/atmos_param/radiation/clouds/get_random_number_stream.F90 $(SRCROOT)/am5_phys/atmos_param/cosp/cosp.F90 $(SRCROOT)/am5_phys/atmos_param/lscale_cond/lscale_cond.F90 $(SRCROOT)/am5_phys/atmos_param/cosp/llnl/prec_scops.F $(SRCROOT)/am5_phys/atmos_param/radiation/driver/radiation_types.F90 $(SRCROOT)/am5_phys/atmos_param/han_edmf/satmedmfvdifq.F90 $(SRCROOT)/am5_phys/atmos_param/han_edmf/hanedmf.F90 $(SRCROOT)/am5_phys/atmos_param/shallow_cu/wet_deposition_0D.F90 $(SRCROOT)/am5_phys/atmos_param/shallow_cu/conv_utilities_k.F90 $(SRCROOT)/am5_phys/atmos_param/radiation/longwave/optical_path.F90 $(SRCROOT)/am5_phys/atmos_param/cosp/cosp_io.F90 $(SRCROOT)/am5_phys/atmos_param/gw_beres/gw_common_utils.F90 $(SRCROOT)/am5_phys/atmos_param/radiation/radiative_gases/ozone.F90 $(SRCROOT)/am5_phys/atmos_param/radiation/driver/radiation_driver.F90 $(SRCROOT)/am5_phys/atmos_param/radiation/clouds/cloudrad_driver.F90 $(SRCROOT)/am5_phys/atmos_param/radiation/shortwave/esfsw_bands.F90 $(SRCROOT)/am5_phys/atmos_shared/tracer_driver/atmos_ozone_tracer.F90 $(SRCROOT)/am5_phys/atmos_param/radiation/driver/shortwave_driver.F90 $(SRCROOT)/am5_phys/atmos_param/microphysics/rotstayn_klein_mp.F90 $(SRCROOT)/am5_phys/atmos_param/lscloud_driver/lscloud_netcdf.F90 $(SRCROOT)/am5_phys/atmos_param/vert_diff/vert_diff.F90 $(SRCROOT)/am5_phys/atmos_shared/tracer_driver/atmos_convection_tracer.F90 $(SRCROOT)/am5_phys/atmos_shared/tracer_driver/atmos_regional_tracer_driver.F90 $(SRCROOT)/am5_phys/atmos_shared/tracer_driver/aer_ccn_act/aer_ccn_act_k.F90 $(SRCROOT)/am5_phys/atmos_param/lscloud_driver/lscloud_types.F90 $(SRCROOT)/am5_phys/atmos_param/lscloud_driver/check_nan.F90 $(SRCROOT)/am5_phys/atmos_shared/tracer_driver/tropchem/isoropiaIIcode.F $(SRCROOT)/am5_phys/atmos_shared/tracer_driver/tropchem/mo_chemdr.F90 $(SRCROOT)/am5_phys/atmos_shared/tracer_driver/tropchem/AM4_chem/moz.mat.F90 $(SRCROOT)/am5_phys/atmos_param/radiation/driver/solar_data_driver.F90 $(SRCROOT)/am5_phys/atmos_shared/tracer_driver/atmos_sea_salt.F90 $(SRCROOT)/am5_phys/atmos_shared/tracer_driver/atmos_carbon_aerosol.F90 $(SRCROOT)/am5_phys/atmos_shared/tracer_driver/atmos_sulfur_hex.F90 $(SRCROOT)/am5_phys/atmos_param/cosp/quickbeam/mrgrnk.F90 $(SRCROOT)/am5_phys/atmos_param/cosp/quickbeam/scale_LUTs_io.F90 $(SRCROOT)/am5_phys/atmos_param/stable_bl_turb/stable_bl_turb.F90 $(SRCROOT)/am5_phys/atmos_param/cloud_generator/betaDistribution.F90 $(SRCROOT)/am5_phys/atmos_param/cosp/cosp_driver.F90 $(SRCROOT)/am5_phys/atmos_shared/tracer_driver/atmos_tracer_utilities.F90 $(SRCROOT)/am5_phys/atmos_shared/tracer_driver/tropchem/mo_imp_slv.F90 $(SRCROOT)/am5_phys/atmos_param/macrophysics/tiedtke_macro.F90 $(SRCROOT)/am5_phys/atmos_shared/atmos_cmip_diag/atmos_cmip_interp.inc $(SRCROOT)/am5_phys/atmos_param/radiation/util/radiation_diag.F90 $(SRCROOT)/am5_phys/atmos_param/shallow_cu/conv_utilities.F90 $(SRCROOT)/am5_phys/atmos_param/radiation/shortwave/esfsw_utilities.F90 $(SRCROOT)/am5_phys/atmos_param/cosp/quickbeam/gases.F90 $(SRCROOT)/am5_phys/atmos_param/convection_driver/convection_driver.F90 $(SRCROOT)/am5_phys/atmos_param/radiation/clouds/isccp_clouds.F90 $(SRCROOT)/am5_phys/atmos_param/cosp/quickbeam/array_lib.F90 $(SRCROOT)/am5_phys/atmos_shared/tracer_driver/atmos_sulfate.F90 $(SRCROOT)/am5_phys/atmos_shared/tracer_driver/tropchem/mo_hook.F90 $(SRCROOT)/am5_phys/atmos_param/shallow_cu/conv_closures.F90 $(SRCROOT)/am5_phys/atmos_shared/tracer_driver/tropchem/mo_jpl.F90 $(SRCROOT)/am5_phys/atmos_param/radiation/driver/radiation_driver_types.F90 $(SRCROOT)/am5_phys/atmos_param/cosp/llnl/llnl_stats.F90 $(SRCROOT)/am5_phys/atmos_param/cosp/icarus-scops-3.7/isccp_cloud_types.F $(SRCROOT)/am5_phys/atmos_param/cosp/quickbeam/optics_lib.F90 $(SRCROOT)/am5_phys/atmos_param/radiation/shortwave/esfsw_parameters.F90 $(SRCROOT)/am5_phys/atmos_param/topo_drag/topo_drag.F90 $(SRCROOT)/am5_phys/atmos_param/microphysics/micro_mg2.F90 $(SRCROOT)/am5_phys/atmos_shared/tracer_driver/tropchem/m_tracname.F90 $(SRCROOT)/am5_phys/atmos_param/lscloud_driver/lscloud_driver.F90 $(SRCROOT)/am5_phys/atmos_param/cloud_rad/cloud_rad.F90 $(SRCROOT)/am5_phys/atmos_param/radiation/aerosols/aerosolrad_types.F90 $(SRCROOT)/am5_phys/atmos_param/cosp/cosp_simulator.F90 $(SRCROOT)/am5_phys/atmos_param/cosp/llnl/cosp_radar.F90 $(SRCROOT)/am5_phys/atmos_param/cosp/cosp_isccp_simulator.F90 $(SRCROOT)/am5_phys/atmos_shared/tracer_driver/tropchem/mo_setinv.F90 $(SRCROOT)/am5_phys/atmos_param/shallow_cu/deep_conv.F90 $(SRCROOT)/am5_phys/atmos_shared/vert_advection/vert_advection.F90 $(SRCROOT)/am5_phys/atmos_param/cosp/cosp_modis_simulator.F90 $(SRCROOT)/am5_phys/atmos_param/moist_processes/moist_processes.F90 $(SRCROOT)/am5_phys/atmos_param/radiation/longwave/longwave_types.F90 $(SRCROOT)/am5_phys/atmos_param/han_edmf/tridi.F90 $(SRCROOT)/am5_phys/atmos_param/cosp/quickbeam/zeff.F90 $(SRCROOT)/am5_phys/atmos_param/cosp/cosp_rttov_simulator.F90 $(SRCROOT)/am5_phys/atmos_shared/tracer_driver/mo_setsox.F90 $(SRCROOT)/am5_phys/atmos_param/radiation/radiative_gases/radiative_gases_types.F90 $(SRCROOT)/am5_phys/atmos_param/damping_driver/damping_driver.F90 $(SRCROOT)/am5_phys/atmos_param/shallow_cu/uw_conv.F90 $(SRCROOT)/am5_phys/atmos_shared/aerosol/aerosol.F90 $(SRCROOT)/am5_phys/atmos_shared/tracer_driver/atmos_nh3_tag.F90 $(SRCROOT)/am5_phys/atmos_param/cosp/quickbeam/radar_simulator.F90 $(SRCROOT)/am5_phys/atmos_shared/tracer_driver/cloud_chem.F90 $(SRCROOT)/am5_phys/atmos_param/radiation/longwave/longwave_utilities.F90 $(SRCROOT)/am5_phys/atmos_shared/tracer_driver/tropchem/mo_usrrxt.F90 $(SRCROOT)/am5_phys/atmos_shared/tracer_driver/stratchem/strat_chem_driver.F90 $(SRCROOT)/am5_phys/atmos_param/vert_diff_driver/vert_diff_driver.F90 $(SRCROOT)/am5_phys/atmos_param/cosp/cosp_utils.F90 $(SRCROOT)/am5_phys/atmos_shared/tracer_driver/tropchem/AM4_chem/moz.mods.F90 $(SRCROOT)/am5_phys/atmos_param/radiation/clouds/cloud_spec.F90 $(SRCROOT)/am5_phys/atmos_param/aerosol_cloud/aerosol_cloud.F90 $(SRCROOT)/am5_phys/atmos_param/lscloud_driver/lscloud_constants.F90 $(SRCROOT)/am5_phys/atmos_shared/tracer_driver/tropchem/mo_chem_utls.F90 $(SRCROOT)/am5_phys/atmos_param/radiation/clouds/cloudrad_types.F90 $(SRCROOT)/am5_phys/atmos_shared/tracer_driver/tropchem/mo_read_sim_chm.F90 $(SRCROOT)/am5_phys/atmos_shared/atmos_cmip_diag/atmos_global_diag.F90 $(SRCROOT)/am5_phys/atmos_shared/tracer_driver/tropchem/mo_exp_slv.F90 $(SRCROOT)/am5_phys/atmos_shared/tracer_driver/tropchem/isrpia.inc $(SRCROOT)/am5_phys/atmos_param/radiation/longwave/longwave_params.F90 $(SRCROOT)/am5_phys/atmos_shared/aerosol/aerosol_types.F90 $(SRCROOT)/am5_phys/atmos_param/radiation/driver/longwave_driver.F90 $(SRCROOT)/am5_phys/atmos_shared/tracer_driver/xactive_bvoc.F90 $(SRCROOT)/am5_phys/atmos_shared/tracer_driver/stratchem/strat_chem_model.F90 $(SRCROOT)/am5_phys/atmos_param/radiation/aerosols/aerosolrad_package.F90 $(SRCROOT)/am5_phys/atmos_param/lscloud_driver/polysvp.F90 $(SRCROOT)/am5_phys/atmos_shared/tracer_driver/atmos_age_tracer.F90 $(SRCROOT)/am5_phys/atmos_param/cosp/quickbeam/radar_simulator_types.F90 $(SRCROOT)/am5_phys/atmos_param/mg_drag/mg_drag.F90 $(SRCROOT)/am5_phys/atmos_shared/tracer_driver/tropchem/AM4_chem/moz.subs.F90 $(SRCROOT)/am5_phys/atmos_shared/tracer_driver/atmos_ch4.F90 $(SRCROOT)/am5_phys/atmos_param/radiation/longwave/longwave_clouds.F90 $(SRCROOT)/am5_phys/atmos_param/convection_driver/detr_ice_num.F90 $(SRCROOT)/am5_phys/atmos_param/cosp/icarus-scops-3.7/icarus.F90 $(SRCROOT)/am5_phys/atmos_param/radiation/radiative_gases/radiative_gases.F90 $(SRCROOT)/am5_phys/atmos_param/cosp/llnl/pf_to_mr.F $(SRCROOT)/am5_phys/atmos_shared/tracer_driver/aer_ccn_act/aer_ccn_act.F90 $(SRCROOT)/am5_phys/atmos_param/radiation/longwave/sealw99.F90 $(SRCROOT)/am5_phys/atmos_param/han_edmf/mfscuq.F90 $(SRCROOT)/am5_phys/atmos_shared/tracer_driver/atmos_tracer_driver.F90 $(SRCROOT)/am5_phys/atmos_param/cosp/MISR_simulator/MISR_simulator.F $(SRCROOT)/am5_phys/atmos_param/radiation/aerosols/aerosolrad_driver.F90 $(SRCROOT)/am5_phys/atmos_param/radiation/clouds/cloudrad_package.F90 $(SRCROOT)/am5_phys/atmos_param/cosp/quickbeam/format_input.F90 $(SRCROOT)/am5_phys/atmos_shared/atmos_nudge/atmos_nudge.F90 $(SRCROOT)/am5_phys/atmos_param/radiation/driver/rad_output_file.F90 $(SRCROOT)/am5_phys/atmos_param/radiation/longwave/lw_gases_stdtf.F90 $(SRCROOT)/am5_phys/atmos_param/cosp/cosp_constants.F90 $(SRCROOT)/am5_phys/atmos_shared/tracer_driver/aer_ccn_act/ice_nucl.F90 $(SRCROOT)/am5_phys/atmos_shared/tracer_driver/aer_ccn_act/aer_in_act.F90 $(SRCROOT)/am5_phys/atmos_param/radiation/clouds/microphys_rad.F90 $(SRCROOT)/am5_phys/atmos_param/cloud_generator/cloud_generator.F90 $(SRCROOT)/am5_phys/atmos_shared/tracer_driver/tropchem/AM4_chem/mo_fastjx.F90 $(SRCROOT)/am5_phys/atmos_shared/tracer_driver/tropchem/mo_chemini.F90 $(SRCROOT)/am5_phys/atmos_param/cosp/quickbeam/math_lib.F90 $(SRCROOT)/am5_phys/atmos_param/radiation/shortwave/esfsw_driver.F90 $(SRCROOT)/am5_phys/atmos_param/han_edmf/mfpbltq.F90 $(SRCROOT)/am5_phys/atmos_param/cosp/cosp_lidar.F90 $(SRCROOT)/am5_phys/atmos_shared/tracer_driver/atmos_co2.F90 $(SRCROOT)FMS/include/fms_platform.h $(SRCROOT)/am5_phys/atmos_param/vert_turb_driver/vert_turb_driver.F90 $(SRCROOT)/am5_phys/atmos_param/cosp/cosp_defs.H $(SRCROOT)/am5_phys/atmos_shared/tracer_driver/tropchem/strat_chem_utilities.F90 $(SRCROOT)/am5_phys/atmos_shared/tracer_driver/atmos_ch3i.F90 $(SRCROOT)/am5_phys/atmos_param/cosp/actsim/lidar_simulator.F90 $(SRCROOT)/am5_phys/atmos_param/cosp/actsim/lmd_ipsl_stats.F90 $(SRCROOT)/am5_phys/atmos_param/shallow_cu/conv_plumes.F90 $(SRCROOT)/am5_phys/atmos_param/shallow_cu/conv_plumes_k.F90 $(SRCROOT)/am5_phys/atmos_param/cosp/icarus-scops-3.7/scops.F $(SRCROOT)/am5_phys/atmos_param/radiation/shortwave/shortwave_types.F90 $(SRCROOT)/am5_phys/atmos_param/cosp/MODIS_simulator/modis_simulator.F90 $(SRCROOT)/am5_phys/atmos_param/cu_mo_trans/cu_mo_trans.F90 $(SRCROOT)/am5_phys/atmos_param/lscloud_driver/lscloud_debug.F90 $(SRCROOT)/am5_phys/atmos_param/radiation/clouds/strat_clouds_W.F90 $(SRCROOT)/am5_phys/atmos_param/radiation/clouds/uw_clouds_W.F90 $(SRCROOT)/am5_phys/atmos_shared/tracer_driver/tropchem/tropchem_types_mod.F90 $(SRCROOT)/am5_phys/atmos_param/cg_drag/cg_drag.F90 $(SRCROOT)/am5_phys/atmos_param/cosp/quickbeam/atmos_lib.F90 $(SRCROOT)/am5_phys/atmos_shared/tracer_driver/atmos_radon.F90 $(SRCROOT)/am5_phys/atmos_param/microphysics/micro_mg2_utils.F90 $(SRCROOT)/am5_phys/atmos_param/radiation/longwave/longwave_fluxes.F90 $(SRCROOT)/am5_phys/atmos_param/gw_beres/gw_beres.F90 $(SRCROOT)/am5_phys/atmos_param/moist_processes/moist_processes_utils.F90 $(SRCROOT)/am5_phys/atmos_param/macrophysics/ls_cloud_macrophysics.F90 $(SRCROOT)/am5_phys/atmos_param/microphysics/ls_cloud_microphysics.F90 $(SRCROOT)/am5_phys/atmos_param/entrain/entrain.F90 $(SRCROOT)/am5_phys/atmos_param/cosp/quickbeam/dsd.F90 $(SRCROOT)/am5_phys/atmos_param/cosp/cosp_misr_simulator.F90 $(SRCROOT)/am5_phys/atmos_param/physics_driver/physics_types.F90 $(SRCROOT)/am5_phys/atmos_shared/tracer_driver/tropchem/AM4_chem/mo_fphoto.F90 $(SRCROOT)/am5_phys/atmos_param/radiation/longwave/gas_tf.F90 $(SRCROOT)/am5_phys/atmos_shared/tracer_driver/atmos_soa.F90 $(SRCROOT)/am5_phys/atmos_shared/tracer_driver/tropchem/aerosol_thermodynamics.F90 $(SRCROOT)/am5_phys/atmos_param/cosp/icarus-scops-3.7/congvec.H $(SRCROOT)/am5_phys/atmos_param/physics_radiation_exch/physics_radiation_exch.F90 $(SRCROOT)/am5_phys/atmos_shared/tracer_driver/aer_ccn_act/aerosol_params.F90 $(SRCROOT)/am5_phys/atmos_param/cosp/quickbeam/calc_Re.F90 $(SRCROOT)/am5_phys/atmos_param/cosp/quickbeam/radar_simulator_init.F90 $(SRCROOT)/am5_phys/atmos_param/physics_driver/physics_driver.F90 $(SRCROOT)/am5_phys/atmos_param/cosp/cosp_diagnostics.F90 $(SRCROOT)/am5_phys/atmos_param/cosp/cosp_types.F90 $(SRCROOT)/am5_phys/atmos_shared/atmos_cmip_diag/atmos_cmip_diag.F90 $(SRCROOT)/am5_phys/atmos_param/radiation/driver/radiation_driver_diag.F90 $(SRCROOT)/am5_phys/atmos_shared/tracer_driver/tropchem_driver.F90 $(SRCROOT)/am5_phys/atmos_param/radiation/clouds/cloudrad_diagnostics.F90 $(SRCROOT)/am5_phys/atmos_shared/tracer_driver/tropchem/mo_photo.F90 $(SRCROOT)/am5_phys/atmos_shared/tracer_driver/tropchem/mo_rodas_slv.F90 $(SRCROOT)/am5_phys/atmos_param/microphysics/gamma_mg.F90 $(SRCROOT)/am5_phys/atmos_param/radiation/longwave/longwave_tables.F90 $(SRCROOT)/am5_phys/atmos_shared/tracer_driver/atmos_tropopause.F90 $(SRCROOT)/am5_phys/atmos_param/convection_driver/convection_utilities.F90
clean: neat
	-rm -f .libam5_phys.a.cppdefs $(OBJ) libam5_phys.a

neat:
	-rm -f $(TMPFILES)

localize: $(OFF)
	cp $(OFF) .

TAGS: $(SRC)
	etags $(SRC)

tags: $(SRC)
	ctags $(SRC)

libam5_phys.a: $(OBJ)
	$(AR) $(ARFLAGS) libam5_phys.a $(OBJ)

