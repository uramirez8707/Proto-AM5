# Makefile created by mkmf 2022.01.00

CPPDEFS = -DINTERNAL_FILE_NML -DMAX_FIELDS_=100 -DNOT_SET_AFFINITY -D_USE_MOM6_DIAG -D_USE_GENERIC_TRACER  -DUSE_PRECISION=2
OTHERFLAGS = -I$(BUILDROOT)/FMS 
include $(MK_TEMPLATE)


.DEFAULT:
	-echo $@ does not exist.
all: libMOM6.a
BFB_initialization.o: $(SRCROOT)/MOM6/src/user/BFB_initialization.F90 $(SRCROOT)/MOM6/config_src/memory/dynamic_nonsymmetric/MOM_memory.h $(SRCROOT)/MOM6/src/framework/MOM_memory_macros.h $(SRCROOT)/MOM6/src/framework/version_variable.h $(SRCROOT)MOM6/src/framework/version_variable.h MOM_error_handler.o MOM_file_parser.o MOM_get_input.o MOM_grid.o MOM_sponge.o MOM_tracer_registry.o MOM_unit_scaling.o MOM_variables.o MOM_verticalGrid.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)/MOM6/config_src/memory/dynamic_nonsymmetric -I$(SRCROOT)/MOM6/src/framework -I$(SRCROOT)MOM6/src/framework	$(SRCROOT)/MOM6/src/user/BFB_initialization.F90
BFB_surface_forcing.o: $(SRCROOT)/MOM6/src/user/BFB_surface_forcing.F90 $(SRCROOT)/MOM6/src/framework/version_variable.h MOM_diag_mediator.o MOM_domains.o MOM_error_handler.o MOM_file_parser.o MOM_forcing_type.o MOM_grid.o MOM_safe_alloc.o MOM_time_manager.o MOM_tracer_flow_control.o MOM_unit_scaling.o MOM_variables.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)/MOM6/src/framework	$(SRCROOT)/MOM6/src/user/BFB_surface_forcing.F90
DOME2d_initialization.o: $(SRCROOT)/MOM6/src/user/DOME2d_initialization.F90 $(SRCROOT)/MOM6/config_src/memory/dynamic_nonsymmetric/MOM_memory.h $(SRCROOT)/MOM6/src/framework/MOM_memory_macros.h $(SRCROOT)/MOM6/src/framework/version_variable.h MOM_ALE_sponge.o MOM_dyn_horgrid.o MOM_error_handler.o MOM_file_parser.o MOM_get_input.o MOM_grid.o MOM_interface_heights.o MOM_sponge.o MOM_unit_scaling.o MOM_variables.o MOM_verticalGrid.o regrid_consts.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)/MOM6/config_src/memory/dynamic_nonsymmetric -I$(SRCROOT)/MOM6/src/framework	$(SRCROOT)/MOM6/src/user/DOME2d_initialization.F90
DOME_initialization.o: $(SRCROOT)/MOM6/src/user/DOME_initialization.F90 $(SRCROOT)/MOM6/config_src/memory/dynamic_nonsymmetric/MOM_memory.h $(SRCROOT)/MOM6/src/framework/MOM_memory_macros.h $(SRCROOT)/MOM6/src/framework/version_variable.h MOM_sponge.o MOM_dyn_horgrid.o MOM_error_handler.o MOM_file_parser.o MOM_get_input.o MOM_grid.o MOM_open_boundary.o MOM_tracer_registry.o MOM_unit_scaling.o MOM_variables.o MOM_verticalGrid.o MOM_EOS.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)/MOM6/config_src/memory/dynamic_nonsymmetric -I$(SRCROOT)/MOM6/src/framework	$(SRCROOT)/MOM6/src/user/DOME_initialization.F90
DOME_tracer.o: $(SRCROOT)/MOM6/src/tracer/DOME_tracer.F90 $(SRCROOT)/MOM6/config_src/memory/dynamic_nonsymmetric/MOM_memory.h $(SRCROOT)/MOM6/src/framework/MOM_memory_macros.h $(SRCROOT)/MOM6/src/framework/version_variable.h MOM_coupler_types.o MOM_diag_mediator.o MOM_error_handler.o MOM_file_parser.o MOM_forcing_type.o MOM_hor_index.o MOM_grid.o MOM_io.o MOM_open_boundary.o MOM_restart.o MOM_sponge.o MOM_time_manager.o MOM_tracer_registry.o MOM_tracer_diabatic.o MOM_unit_scaling.o MOM_variables.o MOM_verticalGrid.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)/MOM6/config_src/memory/dynamic_nonsymmetric -I$(SRCROOT)/MOM6/src/framework	$(SRCROOT)/MOM6/src/tracer/DOME_tracer.F90
FMS_coupler_util.o: $(SRCROOT)/MOM6/config_src/external/GFDL_ocean_BGC/FMS_coupler_util.F90
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)/MOM6/config_src/external/GFDL_ocean_BGC/FMS_coupler_util.F90
ISOMIP_initialization.o: $(SRCROOT)/MOM6/src/user/ISOMIP_initialization.F90 $(SRCROOT)/MOM6/config_src/memory/dynamic_nonsymmetric/MOM_memory.h $(SRCROOT)/MOM6/src/framework/MOM_memory_macros.h $(SRCROOT)/MOM6/src/framework/version_variable.h MOM_ALE_sponge.o MOM_sponge.o MOM_dyn_horgrid.o MOM_error_handler.o MOM_file_parser.o MOM_get_input.o MOM_grid.o MOM_interface_heights.o MOM_io.o MOM_unit_scaling.o MOM_variables.o MOM_verticalGrid.o MOM_EOS.o regrid_consts.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)/MOM6/config_src/memory/dynamic_nonsymmetric -I$(SRCROOT)/MOM6/src/framework	$(SRCROOT)/MOM6/src/user/ISOMIP_initialization.F90
ISOMIP_tracer.o: $(SRCROOT)/MOM6/src/tracer/ISOMIP_tracer.F90 $(SRCROOT)/MOM6/config_src/memory/dynamic_nonsymmetric/MOM_memory.h $(SRCROOT)/MOM6/src/framework/MOM_memory_macros.h $(SRCROOT)/MOM6/src/framework/version_variable.h MOM_coms.o MOM_coupler_types.o MOM_diag_mediator.o MOM_error_handler.o MOM_file_parser.o MOM_forcing_type.o MOM_hor_index.o MOM_grid.o MOM_io.o MOM_open_boundary.o MOM_restart.o MOM_ALE_sponge.o MOM_time_manager.o MOM_tracer_registry.o MOM_tracer_diabatic.o MOM_unit_scaling.o MOM_variables.o MOM_verticalGrid.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)/MOM6/config_src/memory/dynamic_nonsymmetric -I$(SRCROOT)/MOM6/src/framework	$(SRCROOT)/MOM6/src/tracer/ISOMIP_tracer.F90
Idealized_Hurricane.o: $(SRCROOT)/MOM6/src/user/Idealized_Hurricane.F90 $(SRCROOT)/MOM6/config_src/memory/dynamic_nonsymmetric/MOM_memory.h $(SRCROOT)/MOM6/src/framework/MOM_memory_macros.h $(SRCROOT)/MOM6/src/framework/version_variable.h $(SRCROOT)MOM6/src/framework/version_variable.h MOM_error_handler.o MOM_file_parser.o MOM_forcing_type.o MOM_grid.o MOM_safe_alloc.o MOM_time_manager.o MOM_unit_scaling.o MOM_variables.o MOM_verticalGrid.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)/MOM6/config_src/memory/dynamic_nonsymmetric -I$(SRCROOT)/MOM6/src/framework -I$(SRCROOT)MOM6/src/framework	$(SRCROOT)/MOM6/src/user/Idealized_Hurricane.F90
Kelvin_initialization.o: $(SRCROOT)/MOM6/src/user/Kelvin_initialization.F90 $(SRCROOT)/MOM6/config_src/memory/dynamic_nonsymmetric/MOM_memory.h $(SRCROOT)/MOM6/src/framework/MOM_memory_macros.h $(SRCROOT)/MOM6/src/framework/version_variable.h MOM_dyn_horgrid.o MOM_error_handler.o MOM_file_parser.o MOM_grid.o MOM_open_boundary.o MOM_unit_scaling.o MOM_verticalGrid.o MOM_time_manager.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)/MOM6/config_src/memory/dynamic_nonsymmetric -I$(SRCROOT)/MOM6/src/framework	$(SRCROOT)/MOM6/src/user/Kelvin_initialization.F90
MOM.o: $(SRCROOT)/MOM6/src/core/MOM.F90 $(SRCROOT)/MOM6/config_src/memory/dynamic_nonsymmetric/MOM_memory.h $(SRCROOT)/MOM6/src/framework/MOM_memory_macros.h $(SRCROOT)/MOM6/src/framework/version_variable.h MOM_array_transform.o MOM_debugging.o MOM_checksum_packages.o MOM_coms.o MOM_cpu_clock.o MOM_diag_mediator.o MOM_domains.o MOM_error_handler.o MOM_file_parser.o MOM_forcing_type.o MOM_get_input.o MOM_io.o MOM_obsolete_params.o MOM_restart.o MOM_spatial_means.o MOM_time_manager.o MOM_unit_tests.o MOM_ALE.o MOM_ALE_sponge.o MOM_barotropic.o MOM_boundary_update.o MOM_check_scaling.o MOM_coord_initialization.o MOM_diabatic_driver.o MOM_stochastics.o MOM_diagnostics.o MOM_dynamics_unsplit.o MOM_dynamics_split_RK2.o MOM_dynamics_unsplit_RK2.o MOM_dyn_horgrid.o MOM_EOS.o MOM_fixed_initialization.o MOM_grid.o MOM_hor_index.o MOM_interface_heights.o MOM_interface_filter.o MOM_lateral_mixing_coeffs.o MOM_MEKE.o MOM_MEKE_types.o MOM_mixed_layer_restrat.o MOM_obsolete_diagnostics.o MOM_open_boundary.o MOM_porous_barriers.o MOM_set_viscosity.o MOM_shared_initialization.o MOM_sponge.o MOM_state_initialization.o MOM_stoch_eos.o MOM_sum_output.o MOM_thickness_diffuse.o MOM_tracer_advect.o MOM_tracer_hor_diff.o MOM_tracer_registry.o MOM_tracer_flow_control.o MOM_transcribe_grid.o MOM_unit_scaling.o MOM_variables.o MOM_verticalGrid.o MOM_wave_interface.o MOM_database_comms.o MOM_oda_driver.o MOM_oda_incupd.o MOM_offline_main.o MOM_ice_shelf.o MOM_particles.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)/MOM6/config_src/memory/dynamic_nonsymmetric -I$(SRCROOT)/MOM6/src/framework	$(SRCROOT)/MOM6/src/core/MOM.F90
MOM_ALE.o: $(SRCROOT)/MOM6/src/ALE/MOM_ALE.F90 $(SRCROOT)/MOM6/config_src/memory/dynamic_nonsymmetric/MOM_memory.h $(SRCROOT)/MOM6/src/framework/MOM_memory_macros.h MOM_debugging.o MOM_diag_mediator.o MOM_domains.o MOM_error_handler.o MOM_hybgen_unmix.o MOM_hybgen_regrid.o MOM_file_parser.o MOM_interface_heights.o MOM_open_boundary.o MOM_regridding.o MOM_remapping.o MOM_string_functions.o MOM_tracer_registry.o MOM_unit_scaling.o MOM_variables.o MOM_verticalGrid.o regrid_consts.o regrid_edge_values.o PLM_functions.o PPM_functions.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)/MOM6/config_src/memory/dynamic_nonsymmetric -I$(SRCROOT)/MOM6/src/framework	$(SRCROOT)/MOM6/src/ALE/MOM_ALE.F90
MOM_ALE_sponge.o: $(SRCROOT)/MOM6/src/parameterizations/vertical/MOM_ALE_sponge.F90 $(SRCROOT)/MOM6/config_src/memory/dynamic_nonsymmetric/MOM_memory.h $(SRCROOT)/MOM6/src/framework/MOM_memory_macros.h $(SRCROOT)/MOM6/src/framework/version_variable.h $(SRCROOT)MOM6/src/framework/version_variable.h MOM_array_transform.o MOM_coms.o MOM_diag_mediator.o MOM_domains.o MOM_error_handler.o MOM_file_parser.o MOM_grid.o MOM_horizontal_regridding.o MOM_interpolate.o MOM_remapping.o MOM_spatial_means.o MOM_time_manager.o MOM_unit_scaling.o MOM_verticalGrid.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)/MOM6/config_src/memory/dynamic_nonsymmetric -I$(SRCROOT)/MOM6/src/framework -I$(SRCROOT)MOM6/src/framework	$(SRCROOT)/MOM6/src/parameterizations/vertical/MOM_ALE_sponge.F90
MOM_CFC_cap.o: $(SRCROOT)/MOM6/src/tracer/MOM_CFC_cap.F90 $(SRCROOT)/MOM6/config_src/memory/dynamic_nonsymmetric/MOM_memory.h $(SRCROOT)/MOM6/src/framework/MOM_memory_macros.h $(SRCROOT)/MOM6/src/framework/version_variable.h MOM_coms.o MOM_diag_mediator.o MOM_error_handler.o MOM_file_parser.o MOM_forcing_type.o MOM_hor_index.o MOM_grid.o MOM_CVMix_KPP.o MOM_io.o MOM_tracer_registry.o MOM_open_boundary.o MOM_restart.o MOM_spatial_means.o MOM_time_manager.o MOM_interpolate.o MOM_tracer_types.o MOM_tracer_diabatic.o MOM_tracer_Z_init.o MOM_unit_scaling.o MOM_variables.o MOM_verticalGrid.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)/MOM6/config_src/memory/dynamic_nonsymmetric -I$(SRCROOT)/MOM6/src/framework	$(SRCROOT)/MOM6/src/tracer/MOM_CFC_cap.F90
MOM_CVMix_KPP.o: $(SRCROOT)/MOM6/src/parameterizations/vertical/MOM_CVMix_KPP.F90 $(SRCROOT)/MOM6/config_src/memory/dynamic_nonsymmetric/MOM_memory.h $(SRCROOT)/MOM6/src/framework/MOM_memory_macros.h $(SRCROOT)/MOM6/src/framework/version_variable.h MOM_coms.o MOM_debugging.o MOM_diag_mediator.o MOM_error_handler.o MOM_EOS.o MOM_file_parser.o MOM_grid.o MOM_unit_scaling.o MOM_variables.o MOM_verticalGrid.o MOM_wave_interface.o MOM_domains.o MOM_cpu_clock.o MOM_tracer_types.o cvmix_kpp.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)/MOM6/config_src/memory/dynamic_nonsymmetric -I$(SRCROOT)/MOM6/src/framework	$(SRCROOT)/MOM6/src/parameterizations/vertical/MOM_CVMix_KPP.F90
MOM_CVMix_conv.o: $(SRCROOT)/MOM6/src/parameterizations/vertical/MOM_CVMix_conv.F90 $(SRCROOT)/MOM6/config_src/memory/dynamic_nonsymmetric/MOM_memory.h $(SRCROOT)/MOM6/src/framework/MOM_memory_macros.h $(SRCROOT)/MOM6/src/framework/version_variable.h MOM_debugging.o MOM_diag_mediator.o MOM_EOS.o MOM_error_handler.o MOM_file_parser.o MOM_grid.o MOM_unit_scaling.o MOM_variables.o MOM_verticalGrid.o cvmix_convection.o cvmix_kpp.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)/MOM6/config_src/memory/dynamic_nonsymmetric -I$(SRCROOT)/MOM6/src/framework	$(SRCROOT)/MOM6/src/parameterizations/vertical/MOM_CVMix_conv.F90
MOM_CVMix_ddiff.o: $(SRCROOT)/MOM6/src/parameterizations/vertical/MOM_CVMix_ddiff.F90 $(SRCROOT)/MOM6/config_src/memory/dynamic_nonsymmetric/MOM_memory.h $(SRCROOT)/MOM6/src/framework/MOM_memory_macros.h $(SRCROOT)/MOM6/src/framework/version_variable.h MOM_diag_mediator.o MOM_EOS.o MOM_error_handler.o MOM_file_parser.o MOM_debugging.o MOM_grid.o MOM_unit_scaling.o MOM_variables.o MOM_verticalGrid.o cvmix_ddiff.o cvmix_kpp.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)/MOM6/config_src/memory/dynamic_nonsymmetric -I$(SRCROOT)/MOM6/src/framework	$(SRCROOT)/MOM6/src/parameterizations/vertical/MOM_CVMix_ddiff.F90
MOM_CVMix_shear.o: $(SRCROOT)/MOM6/src/parameterizations/vertical/MOM_CVMix_shear.F90 $(SRCROOT)/MOM6/config_src/memory/dynamic_nonsymmetric/MOM_memory.h $(SRCROOT)/MOM6/src/framework/MOM_memory_macros.h $(SRCROOT)/MOM6/src/framework/version_variable.h MOM_diag_mediator.o MOM_error_handler.o MOM_file_parser.o MOM_grid.o MOM_unit_scaling.o MOM_variables.o MOM_verticalGrid.o MOM_EOS.o cvmix_shear.o MOM_kappa_shear.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)/MOM6/config_src/memory/dynamic_nonsymmetric -I$(SRCROOT)/MOM6/src/framework	$(SRCROOT)/MOM6/src/parameterizations/vertical/MOM_CVMix_shear.F90
MOM_CoriolisAdv.o: $(SRCROOT)/MOM6/src/core/MOM_CoriolisAdv.F90 $(SRCROOT)/MOM6/config_src/memory/dynamic_nonsymmetric/MOM_memory.h $(SRCROOT)/MOM6/src/framework/MOM_memory_macros.h $(SRCROOT)/MOM6/src/framework/version_variable.h MOM_diag_mediator.o MOM_error_handler.o MOM_file_parser.o MOM_grid.o MOM_open_boundary.o MOM_string_functions.o MOM_unit_scaling.o MOM_variables.o MOM_verticalGrid.o MOM_wave_interface.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)/MOM6/config_src/memory/dynamic_nonsymmetric -I$(SRCROOT)/MOM6/src/framework	$(SRCROOT)/MOM6/src/core/MOM_CoriolisAdv.F90
MOM_EOS.o: $(SRCROOT)/MOM6/src/equation_of_state/MOM_EOS.F90 $(SRCROOT)/MOM6/src/framework/version_variable.h MOM_EOS_linear.o MOM_EOS_Wright.o MOM_EOS_Wright_full.o MOM_EOS_Wright_red.o MOM_EOS_Jackett06.o MOM_EOS_UNESCO.o MOM_EOS_Roquet_rho.o MOM_EOS_Roquet_SpV.o MOM_EOS_TEOS10.o MOM_temperature_convert.o MOM_TFreeze.o MOM_error_handler.o MOM_file_parser.o MOM_hor_index.o MOM_io.o MOM_string_functions.o MOM_unit_scaling.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)/MOM6/src/framework	$(SRCROOT)/MOM6/src/equation_of_state/MOM_EOS.F90
MOM_EOS_Jackett06.o: $(SRCROOT)/MOM6/src/equation_of_state/MOM_EOS_Jackett06.F90 MOM_hor_index.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)/MOM6/src/equation_of_state/MOM_EOS_Jackett06.F90
MOM_EOS_Roquet_SpV.o: $(SRCROOT)/MOM6/src/equation_of_state/MOM_EOS_Roquet_SpV.F90
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)/MOM6/src/equation_of_state/MOM_EOS_Roquet_SpV.F90
MOM_EOS_Roquet_rho.o: $(SRCROOT)/MOM6/src/equation_of_state/MOM_EOS_Roquet_rho.F90
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)/MOM6/src/equation_of_state/MOM_EOS_Roquet_rho.F90
MOM_EOS_TEOS10.o: $(SRCROOT)/MOM6/src/equation_of_state/MOM_EOS_TEOS10.F90 gsw_mod_toolbox.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)/MOM6/src/equation_of_state/MOM_EOS_TEOS10.F90
MOM_EOS_UNESCO.o: $(SRCROOT)/MOM6/src/equation_of_state/MOM_EOS_UNESCO.F90
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)/MOM6/src/equation_of_state/MOM_EOS_UNESCO.F90
MOM_EOS_Wright.o: $(SRCROOT)/MOM6/src/equation_of_state/MOM_EOS_Wright.F90 MOM_hor_index.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)/MOM6/src/equation_of_state/MOM_EOS_Wright.F90
MOM_EOS_Wright_full.o: $(SRCROOT)/MOM6/src/equation_of_state/MOM_EOS_Wright_full.F90 MOM_hor_index.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)/MOM6/src/equation_of_state/MOM_EOS_Wright_full.F90
MOM_EOS_Wright_red.o: $(SRCROOT)/MOM6/src/equation_of_state/MOM_EOS_Wright_red.F90 MOM_hor_index.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)/MOM6/src/equation_of_state/MOM_EOS_Wright_red.F90
MOM_EOS_linear.o: $(SRCROOT)/MOM6/src/equation_of_state/MOM_EOS_linear.F90 MOM_hor_index.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)/MOM6/src/equation_of_state/MOM_EOS_linear.F90
MOM_MEKE.o: $(SRCROOT)/MOM6/src/parameterizations/lateral/MOM_MEKE.F90 $(SRCROOT)/MOM6/config_src/memory/dynamic_nonsymmetric/MOM_memory.h $(SRCROOT)/MOM6/src/framework/MOM_memory_macros.h $(SRCROOT)/MOM6/src/framework/version_variable.h MOM_coms.o MOM_database_comms.o MOM_debugging.o MOM_cpu_clock.o MOM_diag_mediator.o MOM_domains.o MOM_error_handler.o MOM_file_parser.o MOM_grid.o MOM_hor_index.o MOM_interface_heights.o MOM_interpolate.o MOM_io.o MOM_isopycnal_slopes.o MOM_restart.o MOM_string_functions.o MOM_time_manager.o MOM_unit_scaling.o MOM_variables.o MOM_verticalGrid.o MOM_MEKE_types.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)/MOM6/config_src/memory/dynamic_nonsymmetric -I$(SRCROOT)/MOM6/src/framework	$(SRCROOT)/MOM6/src/parameterizations/lateral/MOM_MEKE.F90
MOM_MEKE_types.o: $(SRCROOT)/MOM6/src/parameterizations/lateral/MOM_MEKE_types.F90
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)/MOM6/src/parameterizations/lateral/MOM_MEKE_types.F90
MOM_OCMIP2_CFC.o: $(SRCROOT)/MOM6/src/tracer/MOM_OCMIP2_CFC.F90 $(SRCROOT)/MOM6/config_src/memory/dynamic_nonsymmetric/MOM_memory.h $(SRCROOT)/MOM6/src/framework/MOM_memory_macros.h $(SRCROOT)/MOM6/src/framework/version_variable.h MOM_coms.o MOM_coupler_types.o MOM_diag_mediator.o MOM_error_handler.o MOM_file_parser.o MOM_forcing_type.o MOM_hor_index.o MOM_grid.o MOM_io.o MOM_open_boundary.o MOM_restart.o MOM_spatial_means.o MOM_sponge.o MOM_time_manager.o MOM_tracer_registry.o MOM_tracer_diabatic.o MOM_tracer_Z_init.o MOM_unit_scaling.o MOM_variables.o MOM_verticalGrid.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)/MOM6/config_src/memory/dynamic_nonsymmetric -I$(SRCROOT)/MOM6/src/framework	$(SRCROOT)/MOM6/src/tracer/MOM_OCMIP2_CFC.F90
MOM_PointAccel.o: $(SRCROOT)/MOM6/src/diagnostics/MOM_PointAccel.F90 $(SRCROOT)/MOM6/config_src/memory/dynamic_nonsymmetric/MOM_memory.h $(SRCROOT)/MOM6/src/framework/MOM_memory_macros.h $(SRCROOT)/MOM6/src/framework/version_variable.h MOM_diag_mediator.o MOM_domains.o MOM_error_handler.o MOM_file_parser.o MOM_get_input.o MOM_grid.o MOM_io.o MOM_time_manager.o MOM_unit_scaling.o MOM_variables.o MOM_verticalGrid.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)/MOM6/config_src/memory/dynamic_nonsymmetric -I$(SRCROOT)/MOM6/src/framework	$(SRCROOT)/MOM6/src/diagnostics/MOM_PointAccel.F90
MOM_PressureForce.o: $(SRCROOT)/MOM6/src/core/MOM_PressureForce.F90 $(SRCROOT)/MOM6/config_src/memory/dynamic_nonsymmetric/MOM_memory.h $(SRCROOT)/MOM6/src/framework/MOM_memory_macros.h $(SRCROOT)/MOM6/src/framework/version_variable.h MOM_diag_mediator.o MOM_error_handler.o MOM_file_parser.o MOM_grid.o MOM_PressureForce_FV.o MOM_PressureForce_Montgomery.o MOM_tidal_forcing.o MOM_unit_scaling.o MOM_variables.o MOM_verticalGrid.o MOM_ALE.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)/MOM6/config_src/memory/dynamic_nonsymmetric -I$(SRCROOT)/MOM6/src/framework	$(SRCROOT)/MOM6/src/core/MOM_PressureForce.F90
MOM_PressureForce_FV.o: $(SRCROOT)/MOM6/src/core/MOM_PressureForce_FV.F90 $(SRCROOT)/MOM6/config_src/memory/dynamic_nonsymmetric/MOM_memory.h $(SRCROOT)/MOM6/src/framework/MOM_memory_macros.h $(SRCROOT)/MOM6/src/framework/version_variable.h MOM_diag_mediator.o MOM_error_handler.o MOM_file_parser.o MOM_grid.o MOM_PressureForce_Montgomery.o MOM_tidal_forcing.o MOM_unit_scaling.o MOM_variables.o MOM_verticalGrid.o MOM_EOS.o MOM_density_integrals.o MOM_ALE.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)/MOM6/config_src/memory/dynamic_nonsymmetric -I$(SRCROOT)/MOM6/src/framework	$(SRCROOT)/MOM6/src/core/MOM_PressureForce_FV.F90
MOM_PressureForce_Montgomery.o: $(SRCROOT)/MOM6/src/core/MOM_PressureForce_Montgomery.F90 $(SRCROOT)/MOM6/config_src/memory/dynamic_nonsymmetric/MOM_memory.h $(SRCROOT)/MOM6/src/framework/MOM_memory_macros.h $(SRCROOT)/MOM6/src/framework/version_variable.h MOM_density_integrals.o MOM_diag_mediator.o MOM_error_handler.o MOM_file_parser.o MOM_grid.o MOM_tidal_forcing.o MOM_unit_scaling.o MOM_variables.o MOM_verticalGrid.o MOM_EOS.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)/MOM6/config_src/memory/dynamic_nonsymmetric -I$(SRCROOT)/MOM6/src/framework	$(SRCROOT)/MOM6/src/core/MOM_PressureForce_Montgomery.F90
MOM_TFreeze.o: $(SRCROOT)/MOM6/src/equation_of_state/MOM_TFreeze.F90 gsw_mod_toolbox.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)/MOM6/src/equation_of_state/MOM_TFreeze.F90
MOM_Zanna_Bolton.o: $(SRCROOT)/MOM6/src/parameterizations/lateral/MOM_Zanna_Bolton.F90 $(SRCROOT)/MOM6/config_src/memory/dynamic_nonsymmetric/MOM_memory.h $(SRCROOT)/MOM6/src/framework/MOM_memory_macros.h $(SRCROOT)/MOM6/src/framework/version_variable.h MOM_grid.o MOM_verticalGrid.o MOM_diag_mediator.o MOM_file_parser.o MOM_unit_scaling.o MOM_domains.o MOM_coms.o MOM_error_handler.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)/MOM6/config_src/memory/dynamic_nonsymmetric -I$(SRCROOT)/MOM6/src/framework	$(SRCROOT)/MOM6/src/parameterizations/lateral/MOM_Zanna_Bolton.F90
MOM_array_transform.o: $(SRCROOT)/MOM6/src/framework/MOM_array_transform.F90
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)/MOM6/src/framework/MOM_array_transform.F90
MOM_barotropic.o: $(SRCROOT)/MOM6/src/core/MOM_barotropic.F90 $(SRCROOT)/MOM6/config_src/memory/dynamic_nonsymmetric/MOM_memory.h $(SRCROOT)/MOM6/src/framework/MOM_memory_macros.h $(SRCROOT)/MOM6/src/framework/version_variable.h MOM_debugging.o MOM_cpu_clock.o MOM_diag_mediator.o MOM_domains.o MOM_error_handler.o MOM_file_parser.o MOM_forcing_type.o MOM_grid.o MOM_hor_index.o MOM_io.o MOM_open_boundary.o MOM_restart.o MOM_tidal_forcing.o MOM_time_manager.o MOM_unit_scaling.o MOM_variables.o MOM_verticalGrid.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)/MOM6/config_src/memory/dynamic_nonsymmetric -I$(SRCROOT)/MOM6/src/framework	$(SRCROOT)/MOM6/src/core/MOM_barotropic.F90
MOM_bkgnd_mixing.o: $(SRCROOT)/MOM6/src/parameterizations/vertical/MOM_bkgnd_mixing.F90 $(SRCROOT)/MOM6/config_src/memory/dynamic_nonsymmetric/MOM_memory.h $(SRCROOT)/MOM6/src/framework/MOM_memory_macros.h $(SRCROOT)/MOM6/src/framework/version_variable.h MOM_debugging.o MOM_diag_mediator.o MOM_error_handler.o MOM_file_parser.o MOM_forcing_type.o MOM_grid.o MOM_unit_scaling.o MOM_verticalGrid.o MOM_variables.o MOM_intrinsic_functions.o cvmix_background.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)/MOM6/config_src/memory/dynamic_nonsymmetric -I$(SRCROOT)/MOM6/src/framework	$(SRCROOT)/MOM6/src/parameterizations/vertical/MOM_bkgnd_mixing.F90
MOM_boundary_update.o: $(SRCROOT)/MOM6/src/core/MOM_boundary_update.F90 $(SRCROOT)/MOM6/config_src/memory/dynamic_nonsymmetric/MOM_memory.h $(SRCROOT)/MOM6/src/framework/MOM_memory_macros.h $(SRCROOT)/MOM6/src/framework/version_variable.h MOM_cpu_clock.o MOM_diag_mediator.o MOM_error_handler.o MOM_file_parser.o MOM_grid.o MOM_dyn_horgrid.o MOM_open_boundary.o MOM_unit_scaling.o MOM_tracer_registry.o MOM_variables.o MOM_verticalGrid.o DOME_initialization.o tidal_bay_initialization.o Kelvin_initialization.o shelfwave_initialization.o dyed_channel_initialization.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)/MOM6/config_src/memory/dynamic_nonsymmetric -I$(SRCROOT)/MOM6/src/framework	$(SRCROOT)/MOM6/src/core/MOM_boundary_update.F90
MOM_bulk_mixed_layer.o: $(SRCROOT)/MOM6/src/parameterizations/vertical/MOM_bulk_mixed_layer.F90 $(SRCROOT)/MOM6/config_src/memory/dynamic_nonsymmetric/MOM_memory.h $(SRCROOT)/MOM6/src/framework/MOM_memory_macros.h $(SRCROOT)/MOM6/src/framework/version_variable.h MOM_cpu_clock.o MOM_diag_mediator.o MOM_domains.o MOM_error_handler.o MOM_file_parser.o MOM_forcing_type.o MOM_grid.o MOM_opacity.o MOM_unit_scaling.o MOM_variables.o MOM_verticalGrid.o MOM_EOS.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)/MOM6/config_src/memory/dynamic_nonsymmetric -I$(SRCROOT)/MOM6/src/framework	$(SRCROOT)/MOM6/src/parameterizations/vertical/MOM_bulk_mixed_layer.F90
MOM_check_scaling.o: $(SRCROOT)/MOM6/src/core/MOM_check_scaling.F90 MOM_error_handler.o MOM_unique_scales.o MOM_unit_scaling.o MOM_verticalGrid.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)/MOM6/src/core/MOM_check_scaling.F90
MOM_checksum_packages.o: $(SRCROOT)/MOM6/src/core/MOM_checksum_packages.F90 $(SRCROOT)/MOM6/config_src/memory/dynamic_nonsymmetric/MOM_memory.h $(SRCROOT)/MOM6/src/framework/MOM_memory_macros.h MOM_coms.o MOM_debugging.o MOM_error_handler.o MOM_grid.o MOM_unit_scaling.o MOM_variables.o MOM_verticalGrid.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)/MOM6/config_src/memory/dynamic_nonsymmetric -I$(SRCROOT)/MOM6/src/framework	$(SRCROOT)/MOM6/src/core/MOM_checksum_packages.F90
MOM_checksums.o: $(SRCROOT)/MOM6/src/framework/MOM_checksums.F90 $(SRCROOT)/MOM6/src/framework/version_variable.h MOM_array_transform.o MOM_coms.o MOM_error_handler.o MOM_file_parser.o MOM_hor_index.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)/MOM6/src/framework	$(SRCROOT)/MOM6/src/framework/MOM_checksums.F90
MOM_coms.o: $(SRCROOT)/MOM6/src/framework/MOM_coms.F90 MOM_coms_infra.o MOM_error_handler.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)/MOM6/src/framework/MOM_coms.F90
MOM_coms_infra.o: $(SRCROOT)/MOM6/config_src/infra/FMS2/MOM_coms_infra.F90
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)/MOM6/config_src/infra/FMS2/MOM_coms_infra.F90
MOM_constants.o: $(SRCROOT)/MOM6/config_src/infra/FMS2/MOM_constants.F90
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)/MOM6/config_src/infra/FMS2/MOM_constants.F90
MOM_continuity.o: $(SRCROOT)/MOM6/src/core/MOM_continuity.F90 $(SRCROOT)/MOM6/config_src/memory/dynamic_nonsymmetric/MOM_memory.h $(SRCROOT)/MOM6/src/framework/MOM_memory_macros.h $(SRCROOT)/MOM6/src/framework/version_variable.h MOM_continuity_PPM.o MOM_diag_mediator.o MOM_error_handler.o MOM_file_parser.o MOM_string_functions.o MOM_grid.o MOM_open_boundary.o MOM_unit_scaling.o MOM_variables.o MOM_verticalGrid.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)/MOM6/config_src/memory/dynamic_nonsymmetric -I$(SRCROOT)/MOM6/src/framework	$(SRCROOT)/MOM6/src/core/MOM_continuity.F90
MOM_continuity_PPM.o: $(SRCROOT)/MOM6/src/core/MOM_continuity_PPM.F90 $(SRCROOT)/MOM6/config_src/memory/dynamic_nonsymmetric/MOM_memory.h $(SRCROOT)/MOM6/src/framework/MOM_memory_macros.h $(SRCROOT)/MOM6/src/framework/version_variable.h MOM_cpu_clock.o MOM_diag_mediator.o MOM_error_handler.o MOM_file_parser.o MOM_grid.o MOM_open_boundary.o MOM_unit_scaling.o MOM_variables.o MOM_verticalGrid.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)/MOM6/config_src/memory/dynamic_nonsymmetric -I$(SRCROOT)/MOM6/src/framework	$(SRCROOT)/MOM6/src/core/MOM_continuity_PPM.F90
MOM_controlled_forcing.o: $(SRCROOT)/MOM6/src/user/MOM_controlled_forcing.F90 $(SRCROOT)/MOM6/config_src/memory/dynamic_nonsymmetric/MOM_memory.h $(SRCROOT)/MOM6/src/framework/MOM_memory_macros.h $(SRCROOT)/MOM6/src/framework/version_variable.h MOM_diag_mediator.o MOM_domains.o MOM_error_handler.o MOM_file_parser.o MOM_forcing_type.o MOM_grid.o MOM_restart.o MOM_time_manager.o MOM_unit_scaling.o MOM_variables.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)/MOM6/config_src/memory/dynamic_nonsymmetric -I$(SRCROOT)/MOM6/src/framework	$(SRCROOT)/MOM6/src/user/MOM_controlled_forcing.F90
MOM_coord_initialization.o: $(SRCROOT)/MOM6/src/initialization/MOM_coord_initialization.F90 $(SRCROOT)/MOM6/src/framework/version_variable.h MOM_debugging.o MOM_EOS.o MOM_error_handler.o MOM_file_parser.o MOM_io.o MOM_string_functions.o MOM_unit_scaling.o MOM_variables.o MOM_verticalGrid.o user_initialization.o BFB_initialization.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)/MOM6/src/framework	$(SRCROOT)/MOM6/src/initialization/MOM_coord_initialization.F90
MOM_coupler_types.o: $(SRCROOT)/MOM6/src/framework/MOM_coupler_types.F90 MOM_couplertype_infra.o MOM_domain_infra.o MOM_time_manager.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)/MOM6/src/framework/MOM_coupler_types.F90
MOM_couplertype_infra.o: $(SRCROOT)/MOM6/config_src/infra/FMS2/MOM_couplertype_infra.F90 MOM_domain_infra.o MOM_time_manager.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)/MOM6/config_src/infra/FMS2/MOM_couplertype_infra.F90
MOM_cpu_clock.o: $(SRCROOT)/MOM6/src/framework/MOM_cpu_clock.F90 MOM_cpu_clock_infra.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)/MOM6/src/framework/MOM_cpu_clock.F90
MOM_cpu_clock_infra.o: $(SRCROOT)/MOM6/config_src/infra/FMS2/MOM_cpu_clock_infra.F90
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)/MOM6/config_src/infra/FMS2/MOM_cpu_clock_infra.F90
MOM_data_override.o: $(SRCROOT)/MOM6/src/framework/MOM_data_override.F90 MOM_data_override_infra.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)/MOM6/src/framework/MOM_data_override.F90
MOM_data_override_infra.o: $(SRCROOT)/MOM6/config_src/infra/FMS2/MOM_data_override_infra.F90 MOM_domain_infra.o MOM_time_manager.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)/MOM6/config_src/infra/FMS2/MOM_data_override_infra.F90
MOM_database_comms.o: $(SRCROOT)/MOM6/config_src/external/database_comms/MOM_database_comms.F90 MOM_file_parser.o MOM_error_handler.o database_client_interface.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)/MOM6/config_src/external/database_comms/MOM_database_comms.F90
MOM_debugging.o: $(SRCROOT)/MOM6/src/diagnostics/MOM_debugging.F90 $(SRCROOT)/MOM6/src/framework/version_variable.h MOM_checksums.o MOM_coms.o MOM_domains.o MOM_error_handler.o MOM_file_parser.o MOM_grid.o MOM_hor_index.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)/MOM6/src/framework	$(SRCROOT)/MOM6/src/diagnostics/MOM_debugging.F90
MOM_density_integrals.o: $(SRCROOT)/MOM6/src/core/MOM_density_integrals.F90 $(SRCROOT)/MOM6/config_src/memory/dynamic_nonsymmetric/MOM_memory.h $(SRCROOT)/MOM6/src/framework/MOM_memory_macros.h MOM_EOS.o MOM_error_handler.o MOM_hor_index.o MOM_string_functions.o MOM_variables.o MOM_unit_scaling.o MOM_verticalGrid.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)/MOM6/config_src/memory/dynamic_nonsymmetric -I$(SRCROOT)/MOM6/src/framework	$(SRCROOT)/MOM6/src/core/MOM_density_integrals.F90
MOM_diabatic_aux.o: $(SRCROOT)/MOM6/src/parameterizations/vertical/MOM_diabatic_aux.F90 $(SRCROOT)/MOM6/config_src/memory/dynamic_nonsymmetric/MOM_memory.h $(SRCROOT)/MOM6/src/framework/MOM_memory_macros.h $(SRCROOT)/MOM6/src/framework/version_variable.h MOM_cpu_clock.o MOM_diag_mediator.o MOM_EOS.o MOM_error_handler.o MOM_file_parser.o MOM_forcing_type.o MOM_grid.o MOM_interpolate.o MOM_io.o MOM_opacity.o MOM_tracer_flow_control.o MOM_unit_scaling.o MOM_variables.o MOM_verticalGrid.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)/MOM6/config_src/memory/dynamic_nonsymmetric -I$(SRCROOT)/MOM6/src/framework	$(SRCROOT)/MOM6/src/parameterizations/vertical/MOM_diabatic_aux.F90
MOM_diabatic_driver.o: $(SRCROOT)/MOM6/src/parameterizations/vertical/MOM_diabatic_driver.F90 $(SRCROOT)/MOM6/config_src/memory/dynamic_nonsymmetric/MOM_memory.h $(SRCROOT)/MOM6/src/framework/MOM_memory_macros.h $(SRCROOT)/MOM6/src/framework/version_variable.h $(SRCROOT)MOM6/src/framework/version_variable.h MOM_bulk_mixed_layer.o MOM_debugging.o MOM_checksum_packages.o MOM_cpu_clock.o MOM_CVMix_shear.o MOM_CVMix_ddiff.o MOM_diabatic_aux.o MOM_diag_mediator.o MOM_diapyc_energy_req.o MOM_CVMix_conv.o MOM_domains.o MOM_energetic_PBL.o MOM_entrain_diffusive.o MOM_EOS.o MOM_error_handler.o MOM_file_parser.o MOM_forcing_type.o MOM_geothermal.o MOM_grid.o MOM_internal_tide_input.o MOM_interface_heights.o MOM_internal_tides.o MOM_kappa_shear.o MOM_CVMix_KPP.o MOM_oda_incupd.o MOM_opacity.o MOM_open_boundary.o MOM_regularize_layers.o MOM_set_diffusivity.o MOM_sponge.o MOM_ALE_sponge.o MOM_time_manager.o MOM_tracer_flow_control.o MOM_tracer_diabatic.o MOM_unit_scaling.o MOM_variables.o MOM_verticalGrid.o MOM_wave_interface.o MOM_stochastics.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)/MOM6/config_src/memory/dynamic_nonsymmetric -I$(SRCROOT)/MOM6/src/framework -I$(SRCROOT)MOM6/src/framework	$(SRCROOT)/MOM6/src/parameterizations/vertical/MOM_diabatic_driver.F90
MOM_diag_manager_infra.o: $(SRCROOT)/MOM6/config_src/infra/FMS2/MOM_diag_manager_infra.F90 MOM_time_manager.o MOM_domain_infra.o MOM_error_infra.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)/MOM6/config_src/infra/FMS2/MOM_diag_manager_infra.F90
MOM_diag_mediator.o: $(SRCROOT)/MOM6/src/framework/MOM_diag_mediator.F90 $(SRCROOT)/MOM6/src/framework/version_variable.h MOM_checksums.o MOM_coms.o MOM_cpu_clock.o MOM_diag_manager_infra.o MOM_diag_remap.o MOM_EOS.o MOM_error_handler.o MOM_file_parser.o MOM_grid.o MOM_io.o MOM_safe_alloc.o MOM_string_functions.o MOM_time_manager.o MOM_unit_scaling.o MOM_verticalGrid.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)/MOM6/src/framework	$(SRCROOT)/MOM6/src/framework/MOM_diag_mediator.F90
MOM_diag_remap.o: $(SRCROOT)/MOM6/src/framework/MOM_diag_remap.F90 MOM_coms.o MOM_error_handler.o MOM_debugging.o MOM_diag_manager_infra.o MOM_file_parser.o MOM_string_functions.o MOM_grid.o MOM_unit_scaling.o MOM_verticalGrid.o MOM_EOS.o MOM_remapping.o MOM_regridding.o regrid_consts.o coord_zlike.o coord_sigma.o coord_rho.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)/MOM6/src/framework/MOM_diag_remap.F90
MOM_diagnostics.o: $(SRCROOT)/MOM6/src/diagnostics/MOM_diagnostics.F90 $(SRCROOT)/MOM6/config_src/memory/dynamic_nonsymmetric/MOM_memory.h $(SRCROOT)/MOM6/src/framework/MOM_memory_macros.h $(SRCROOT)/MOM6/src/framework/version_variable.h MOM_coms.o MOM_coupler_types.o MOM_density_integrals.o MOM_diag_mediator.o MOM_domains.o MOM_EOS.o MOM_error_handler.o MOM_file_parser.o MOM_grid.o MOM_interface_heights.o MOM_spatial_means.o MOM_tracer_registry.o MOM_unit_scaling.o MOM_variables.o MOM_verticalGrid.o MOM_wave_speed.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)/MOM6/config_src/memory/dynamic_nonsymmetric -I$(SRCROOT)/MOM6/src/framework	$(SRCROOT)/MOM6/src/diagnostics/MOM_diagnostics.F90
MOM_diapyc_energy_req.o: $(SRCROOT)/MOM6/src/parameterizations/vertical/MOM_diapyc_energy_req.F90 $(SRCROOT)/MOM6/src/framework/version_variable.h MOM_diag_mediator.o MOM_error_handler.o MOM_file_parser.o MOM_grid.o MOM_unit_scaling.o MOM_variables.o MOM_verticalGrid.o MOM_EOS.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)/MOM6/src/framework	$(SRCROOT)/MOM6/src/parameterizations/vertical/MOM_diapyc_energy_req.F90
MOM_document.o: $(SRCROOT)/MOM6/src/framework/MOM_document.F90 MOM_time_manager.o MOM_error_handler.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)/MOM6/src/framework/MOM_document.F90
MOM_domain_infra.o: $(SRCROOT)/MOM6/config_src/infra/FMS2/MOM_domain_infra.F90 MOM_coms_infra.o MOM_cpu_clock_infra.o MOM_error_infra.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)/MOM6/config_src/infra/FMS2/MOM_domain_infra.F90
MOM_domains.o: $(SRCROOT)/MOM6/src/framework/MOM_domains.F90 $(SRCROOT)/MOM6/src/framework/version_variable.h MOM_coms_infra.o MOM_domain_infra.o MOM_error_handler.o MOM_file_parser.o MOM_io_infra.o MOM_string_functions.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)/MOM6/src/framework	$(SRCROOT)/MOM6/src/framework/MOM_domains.F90
MOM_dyn_horgrid.o: $(SRCROOT)/MOM6/src/framework/MOM_dyn_horgrid.F90 MOM_array_transform.o MOM_domains.o MOM_error_handler.o MOM_hor_index.o MOM_unit_scaling.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)/MOM6/src/framework/MOM_dyn_horgrid.F90
MOM_dynamics_split_RK2.o: $(SRCROOT)/MOM6/src/core/MOM_dynamics_split_RK2.F90 $(SRCROOT)/MOM6/config_src/memory/dynamic_nonsymmetric/MOM_memory.h $(SRCROOT)/MOM6/src/framework/MOM_memory_macros.h $(SRCROOT)/MOM6/src/framework/version_variable.h MOM_variables.o MOM_forcing_type.o MOM_checksum_packages.o MOM_cpu_clock.o MOM_diag_mediator.o MOM_domains.o MOM_debugging.o MOM_error_handler.o MOM_file_parser.o MOM_get_input.o MOM_io.o MOM_restart.o MOM_time_manager.o MOM_ALE.o MOM_barotropic.o MOM_boundary_update.o MOM_continuity.o MOM_CoriolisAdv.o MOM_grid.o MOM_hor_index.o MOM_hor_visc.o MOM_interface_heights.o MOM_lateral_mixing_coeffs.o MOM_MEKE_types.o MOM_open_boundary.o MOM_PressureForce.o MOM_set_viscosity.o MOM_thickness_diffuse.o MOM_tidal_forcing.o MOM_unit_scaling.o MOM_vert_friction.o MOM_verticalGrid.o MOM_wave_interface.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)/MOM6/config_src/memory/dynamic_nonsymmetric -I$(SRCROOT)/MOM6/src/framework	$(SRCROOT)/MOM6/src/core/MOM_dynamics_split_RK2.F90
MOM_dynamics_unsplit.o: $(SRCROOT)/MOM6/src/core/MOM_dynamics_unsplit.F90 $(SRCROOT)/MOM6/config_src/memory/dynamic_nonsymmetric/MOM_memory.h $(SRCROOT)/MOM6/src/framework/MOM_memory_macros.h $(SRCROOT)/MOM6/src/framework/version_variable.h MOM_variables.o MOM_forcing_type.o MOM_checksum_packages.o MOM_cpu_clock.o MOM_diag_mediator.o MOM_domains.o MOM_error_handler.o MOM_file_parser.o MOM_get_input.o MOM_time_manager.o MOM_ALE.o MOM_barotropic.o MOM_boundary_update.o MOM_continuity.o MOM_CoriolisAdv.o MOM_debugging.o MOM_grid.o MOM_hor_index.o MOM_hor_visc.o MOM_interface_heights.o MOM_lateral_mixing_coeffs.o MOM_MEKE_types.o MOM_open_boundary.o MOM_PressureForce.o MOM_set_viscosity.o MOM_tidal_forcing.o MOM_unit_scaling.o MOM_vert_friction.o MOM_verticalGrid.o MOM_wave_interface.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)/MOM6/config_src/memory/dynamic_nonsymmetric -I$(SRCROOT)/MOM6/src/framework	$(SRCROOT)/MOM6/src/core/MOM_dynamics_unsplit.F90
MOM_dynamics_unsplit_RK2.o: $(SRCROOT)/MOM6/src/core/MOM_dynamics_unsplit_RK2.F90 $(SRCROOT)/MOM6/config_src/memory/dynamic_nonsymmetric/MOM_memory.h $(SRCROOT)/MOM6/src/framework/MOM_memory_macros.h $(SRCROOT)/MOM6/src/framework/version_variable.h MOM_variables.o MOM_forcing_type.o MOM_checksum_packages.o MOM_cpu_clock.o MOM_diag_mediator.o MOM_domains.o MOM_error_handler.o MOM_file_parser.o MOM_get_input.o MOM_time_manager.o MOM_ALE.o MOM_boundary_update.o MOM_barotropic.o MOM_continuity.o MOM_CoriolisAdv.o MOM_debugging.o MOM_grid.o MOM_hor_index.o MOM_hor_visc.o MOM_interface_heights.o MOM_lateral_mixing_coeffs.o MOM_MEKE_types.o MOM_open_boundary.o MOM_PressureForce.o MOM_set_viscosity.o MOM_tidal_forcing.o MOM_unit_scaling.o MOM_vert_friction.o MOM_verticalGrid.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)/MOM6/config_src/memory/dynamic_nonsymmetric -I$(SRCROOT)/MOM6/src/framework	$(SRCROOT)/MOM6/src/core/MOM_dynamics_unsplit_RK2.F90
MOM_energetic_PBL.o: $(SRCROOT)/MOM6/src/parameterizations/vertical/MOM_energetic_PBL.F90 $(SRCROOT)/MOM6/config_src/memory/dynamic_nonsymmetric/MOM_memory.h $(SRCROOT)/MOM6/src/framework/MOM_memory_macros.h $(SRCROOT)/MOM6/src/framework/version_variable.h MOM_cpu_clock.o MOM_coms.o MOM_diag_mediator.o MOM_domains.o MOM_error_handler.o MOM_file_parser.o MOM_forcing_type.o MOM_grid.o MOM_string_functions.o MOM_unit_scaling.o MOM_variables.o MOM_verticalGrid.o MOM_wave_interface.o MOM_stochastics.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)/MOM6/config_src/memory/dynamic_nonsymmetric -I$(SRCROOT)/MOM6/src/framework	$(SRCROOT)/MOM6/src/parameterizations/vertical/MOM_energetic_PBL.F90
MOM_ensemble_manager.o: $(SRCROOT)/MOM6/src/framework/MOM_ensemble_manager.F90 MOM_ensemble_manager_infra.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)/MOM6/src/framework/MOM_ensemble_manager.F90
MOM_ensemble_manager_infra.o: $(SRCROOT)/MOM6/config_src/infra/FMS2/MOM_ensemble_manager_infra.F90
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)/MOM6/config_src/infra/FMS2/MOM_ensemble_manager_infra.F90
MOM_entrain_diffusive.o: $(SRCROOT)/MOM6/src/parameterizations/vertical/MOM_entrain_diffusive.F90 $(SRCROOT)/MOM6/config_src/memory/dynamic_nonsymmetric/MOM_memory.h $(SRCROOT)/MOM6/src/framework/MOM_memory_macros.h $(SRCROOT)/MOM6/src/framework/version_variable.h MOM_diag_mediator.o MOM_error_handler.o MOM_file_parser.o MOM_forcing_type.o MOM_grid.o MOM_unit_scaling.o MOM_variables.o MOM_verticalGrid.o MOM_EOS.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)/MOM6/config_src/memory/dynamic_nonsymmetric -I$(SRCROOT)/MOM6/src/framework	$(SRCROOT)/MOM6/src/parameterizations/vertical/MOM_entrain_diffusive.F90
MOM_error_handler.o: $(SRCROOT)/MOM6/src/framework/MOM_error_handler.F90 MOM_coms_infra.o MOM_error_infra.o posix.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)/MOM6/src/framework/MOM_error_handler.F90
MOM_error_infra.o: $(SRCROOT)/MOM6/config_src/infra/FMS2/MOM_error_infra.F90
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)/MOM6/config_src/infra/FMS2/MOM_error_infra.F90
MOM_file_parser.o: $(SRCROOT)/MOM6/src/framework/MOM_file_parser.F90 $(SRCROOT)/MOM6/src/framework/version_variable.h MOM_coms.o MOM_error_handler.o MOM_time_manager.o MOM_document.o MOM_string_functions.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)/MOM6/src/framework	$(SRCROOT)/MOM6/src/framework/MOM_file_parser.F90
MOM_file_parser_tests.o: $(SRCROOT)/MOM6/src/framework/testing/MOM_file_parser_tests.F90 posix.o MOM_file_parser.o MOM_time_manager.o MOM_error_handler.o MOM_unit_testing.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)/MOM6/src/framework/testing/MOM_file_parser_tests.F90
MOM_fixed_initialization.o: $(SRCROOT)/MOM6/src/initialization/MOM_fixed_initialization.F90 $(SRCROOT)/MOM6/src/framework/version_variable.h MOM_debugging.o MOM_domains.o MOM_dyn_horgrid.o MOM_error_handler.o MOM_file_parser.o MOM_io.o MOM_grid_initialize.o MOM_open_boundary.o MOM_shared_initialization.o MOM_unit_scaling.o user_initialization.o DOME_initialization.o ISOMIP_initialization.o basin_builder.o benchmark_initialization.o Neverworld_initialization.o DOME2d_initialization.o Kelvin_initialization.o sloshing_initialization.o seamount_initialization.o dumbbell_initialization.o shelfwave_initialization.o Phillips_initialization.o dense_water_initialization.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)/MOM6/src/framework	$(SRCROOT)/MOM6/src/initialization/MOM_fixed_initialization.F90
MOM_forcing_type.o: $(SRCROOT)/MOM6/src/core/MOM_forcing_type.F90 $(SRCROOT)/MOM6/config_src/memory/dynamic_nonsymmetric/MOM_memory.h $(SRCROOT)/MOM6/src/framework/MOM_memory_macros.h MOM_array_transform.o MOM_coupler_types.o MOM_cpu_clock.o MOM_debugging.o MOM_diag_mediator.o MOM_EOS.o MOM_error_handler.o MOM_file_parser.o MOM_grid.o MOM_interface_heights.o MOM_opacity.o MOM_spatial_means.o MOM_unit_scaling.o MOM_variables.o MOM_verticalGrid.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)/MOM6/config_src/memory/dynamic_nonsymmetric -I$(SRCROOT)/MOM6/src/framework	$(SRCROOT)/MOM6/src/core/MOM_forcing_type.F90
MOM_full_convection.o: $(SRCROOT)/MOM6/src/parameterizations/vertical/MOM_full_convection.F90 $(SRCROOT)/MOM6/config_src/memory/dynamic_nonsymmetric/MOM_memory.h $(SRCROOT)/MOM6/src/framework/MOM_memory_macros.h MOM_grid.o MOM_unit_scaling.o MOM_variables.o MOM_verticalGrid.o MOM_EOS.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)/MOM6/config_src/memory/dynamic_nonsymmetric -I$(SRCROOT)/MOM6/src/framework	$(SRCROOT)/MOM6/src/parameterizations/vertical/MOM_full_convection.F90
MOM_generic_tracer.o: $(SRCROOT)/MOM6/src/tracer/MOM_generic_tracer.F90 $(SRCROOT)/MOM6/config_src/memory/dynamic_nonsymmetric/MOM_memory.h $(SRCROOT)/MOM6/src/framework/MOM_memory_macros.h $(SRCROOT)/MOM6/src/framework/version_variable.h $(SRCROOT)MOM6/src/framework/version_variable.h generic_tracer.o generic_tracer_utils.o MOM_ALE_sponge.o MOM_coms.o MOM_diag_mediator.o MOM_error_handler.o MOM_file_parser.o MOM_forcing_type.o MOM_grid.o MOM_hor_index.o MOM_io.o MOM_open_boundary.o MOM_restart.o MOM_spatial_means.o MOM_sponge.o MOM_time_manager.o MOM_tracer_diabatic.o MOM_tracer_registry.o MOM_tracer_Z_init.o MOM_tracer_initialization_from_Z.o MOM_unit_scaling.o MOM_variables.o MOM_verticalGrid.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)/MOM6/config_src/memory/dynamic_nonsymmetric -I$(SRCROOT)/MOM6/src/framework -I$(SRCROOT)MOM6/src/framework	$(SRCROOT)/MOM6/src/tracer/MOM_generic_tracer.F90
MOM_geothermal.o: $(SRCROOT)/MOM6/src/parameterizations/vertical/MOM_geothermal.F90 $(SRCROOT)/MOM6/config_src/memory/dynamic_nonsymmetric/MOM_memory.h $(SRCROOT)/MOM6/src/framework/MOM_memory_macros.h $(SRCROOT)/MOM6/src/framework/version_variable.h MOM_diag_mediator.o MOM_domains.o MOM_error_handler.o MOM_file_parser.o MOM_io.o MOM_grid.o MOM_unit_scaling.o MOM_variables.o MOM_verticalGrid.o MOM_EOS.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)/MOM6/config_src/memory/dynamic_nonsymmetric -I$(SRCROOT)/MOM6/src/framework	$(SRCROOT)/MOM6/src/parameterizations/vertical/MOM_geothermal.F90
MOM_get_input.o: $(SRCROOT)/MOM6/src/framework/MOM_get_input.F90 MOM_error_handler.o MOM_file_parser.o MOM_io.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)/MOM6/src/framework/MOM_get_input.F90
MOM_grid.o: $(SRCROOT)/MOM6/src/core/MOM_grid.F90 $(SRCROOT)/MOM6/config_src/memory/dynamic_nonsymmetric/MOM_memory.h $(SRCROOT)/MOM6/src/framework/MOM_memory_macros.h $(SRCROOT)/MOM6/src/framework/version_variable.h MOM_hor_index.o MOM_domains.o MOM_error_handler.o MOM_file_parser.o MOM_unit_scaling.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)/MOM6/config_src/memory/dynamic_nonsymmetric -I$(SRCROOT)/MOM6/src/framework	$(SRCROOT)/MOM6/src/core/MOM_grid.F90
MOM_grid_initialize.o: $(SRCROOT)/MOM6/src/initialization/MOM_grid_initialize.F90 $(SRCROOT)/MOM6/src/framework/version_variable.h MOM_checksums.o MOM_domains.o MOM_dyn_horgrid.o MOM_error_handler.o MOM_file_parser.o MOM_io.o MOM_unit_scaling.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)/MOM6/src/framework	$(SRCROOT)/MOM6/src/initialization/MOM_grid_initialize.F90
MOM_hor_index.o: $(SRCROOT)/MOM6/src/framework/MOM_hor_index.F90 $(SRCROOT)/MOM6/src/framework/version_variable.h MOM_domains.o MOM_error_handler.o MOM_file_parser.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)/MOM6/src/framework	$(SRCROOT)/MOM6/src/framework/MOM_hor_index.F90
MOM_hor_visc.o: $(SRCROOT)/MOM6/src/parameterizations/lateral/MOM_hor_visc.F90 $(SRCROOT)/MOM6/config_src/memory/dynamic_nonsymmetric/MOM_memory.h $(SRCROOT)/MOM6/src/framework/MOM_memory_macros.h $(SRCROOT)/MOM6/src/framework/version_variable.h MOM_checksums.o MOM_coms.o MOM_diag_mediator.o MOM_domains.o MOM_error_handler.o MOM_file_parser.o MOM_grid.o MOM_lateral_mixing_coeffs.o MOM_barotropic.o MOM_thickness_diffuse.o MOM_io.o MOM_MEKE_types.o MOM_open_boundary.o MOM_unit_scaling.o MOM_verticalGrid.o MOM_variables.o MOM_Zanna_Bolton.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)/MOM6/config_src/memory/dynamic_nonsymmetric -I$(SRCROOT)/MOM6/src/framework	$(SRCROOT)/MOM6/src/parameterizations/lateral/MOM_hor_visc.F90
MOM_horizontal_regridding.o: $(SRCROOT)/MOM6/src/framework/MOM_horizontal_regridding.F90 $(SRCROOT)/MOM6/config_src/memory/dynamic_nonsymmetric/MOM_memory.h $(SRCROOT)/MOM6/src/framework/MOM_memory_macros.h MOM_debugging.o MOM_coms.o MOM_cpu_clock.o MOM_domains.o MOM_error_handler.o MOM_file_parser.o MOM_grid.o MOM_interpolate.o MOM_interp_infra.o MOM_time_manager.o MOM_io.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)/MOM6/config_src/memory/dynamic_nonsymmetric -I$(SRCROOT)/MOM6/src/framework	$(SRCROOT)/MOM6/src/framework/MOM_horizontal_regridding.F90
MOM_hybgen_regrid.o: $(SRCROOT)/MOM6/src/ALE/MOM_hybgen_regrid.F90 $(SRCROOT)/MOM6/config_src/memory/dynamic_nonsymmetric/MOM_memory.h $(SRCROOT)/MOM6/src/framework/MOM_memory_macros.h MOM_EOS.o MOM_error_handler.o MOM_file_parser.o MOM_io.o MOM_string_functions.o MOM_unit_scaling.o MOM_variables.o MOM_verticalGrid.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)/MOM6/config_src/memory/dynamic_nonsymmetric -I$(SRCROOT)/MOM6/src/framework	$(SRCROOT)/MOM6/src/ALE/MOM_hybgen_regrid.F90
MOM_hybgen_remap.o: $(SRCROOT)/MOM6/src/ALE/MOM_hybgen_remap.F90
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)/MOM6/src/ALE/MOM_hybgen_remap.F90
MOM_hybgen_unmix.o: $(SRCROOT)/MOM6/src/ALE/MOM_hybgen_unmix.F90 $(SRCROOT)/MOM6/config_src/memory/dynamic_nonsymmetric/MOM_memory.h $(SRCROOT)/MOM6/src/framework/MOM_memory_macros.h MOM_EOS.o MOM_error_handler.o MOM_file_parser.o MOM_hybgen_regrid.o MOM_tracer_registry.o MOM_unit_scaling.o MOM_variables.o MOM_verticalGrid.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)/MOM6/config_src/memory/dynamic_nonsymmetric -I$(SRCROOT)/MOM6/src/framework	$(SRCROOT)/MOM6/src/ALE/MOM_hybgen_unmix.F90
MOM_ice_shelf.o: $(SRCROOT)/MOM6/src/ice_shelf/MOM_ice_shelf.F90 $(SRCROOT)/MOM6/config_src/memory/dynamic_nonsymmetric/MOM_memory.h $(SRCROOT)/MOM6/src/framework/MOM_memory_macros.h $(SRCROOT)/MOM6/src/framework/version_variable.h MOM_array_transform.o MOM_constants.o MOM_cpu_clock.o MOM_coms.o MOM_data_override.o MOM_diag_mediator.o MOM_ice_shelf_diag_mediator.o MOM_domains.o MOM_dyn_horgrid.o MOM_error_handler.o MOM_file_parser.o MOM_grid.o MOM_grid_initialize.o MOM_hor_index.o MOM_fixed_initialization.o user_initialization.o MOM_io.o MOM_restart.o MOM_time_manager.o MOM_transcribe_grid.o MOM_unit_scaling.o MOM_variables.o MOM_forcing_type.o MOM_get_input.o MOM_EOS.o MOM_ice_shelf_dynamics.o MOM_ice_shelf_initialize.o MOM_ice_shelf_state.o user_shelf_init.o MOM_spatial_means.o MOM_checksums.o MOM_interpolate.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)/MOM6/config_src/memory/dynamic_nonsymmetric -I$(SRCROOT)/MOM6/src/framework	$(SRCROOT)/MOM6/src/ice_shelf/MOM_ice_shelf.F90
MOM_ice_shelf_diag_mediator.o: $(SRCROOT)/MOM6/src/ice_shelf/MOM_ice_shelf_diag_mediator.F90 $(SRCROOT)/MOM6/src/framework/version_variable.h MOM_coms.o MOM_diag_manager_infra.o MOM_error_handler.o MOM_file_parser.o MOM_grid.o MOM_safe_alloc.o MOM_string_functions.o MOM_time_manager.o MOM_unit_scaling.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)/MOM6/src/framework	$(SRCROOT)/MOM6/src/ice_shelf/MOM_ice_shelf_diag_mediator.F90
MOM_ice_shelf_dynamics.o: $(SRCROOT)/MOM6/src/ice_shelf/MOM_ice_shelf_dynamics.F90 $(SRCROOT)/MOM6/config_src/memory/dynamic_nonsymmetric/MOM_memory.h $(SRCROOT)/MOM6/src/framework/MOM_memory_macros.h $(SRCROOT)/MOM6/src/framework/version_variable.h MOM_cpu_clock.o MOM_ice_shelf_diag_mediator.o MOM_domains.o MOM_error_handler.o MOM_file_parser.o MOM_grid.o MOM_io.o MOM_restart.o MOM_time_manager.o MOM_unit_scaling.o MOM_ice_shelf_state.o MOM_coms.o MOM_checksums.o MOM_ice_shelf_initialize.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)/MOM6/config_src/memory/dynamic_nonsymmetric -I$(SRCROOT)/MOM6/src/framework	$(SRCROOT)/MOM6/src/ice_shelf/MOM_ice_shelf_dynamics.F90
MOM_ice_shelf_initialize.o: $(SRCROOT)/MOM6/src/ice_shelf/MOM_ice_shelf_initialize.F90 $(SRCROOT)/MOM6/config_src/memory/dynamic_nonsymmetric/MOM_memory.h $(SRCROOT)/MOM6/src/framework/MOM_memory_macros.h MOM_grid.o MOM_array_transform.o MOM_hor_index.o MOM_file_parser.o MOM_io.o MOM_error_handler.o MOM_unit_scaling.o user_shelf_init.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)/MOM6/config_src/memory/dynamic_nonsymmetric -I$(SRCROOT)/MOM6/src/framework	$(SRCROOT)/MOM6/src/ice_shelf/MOM_ice_shelf_initialize.F90
MOM_ice_shelf_state.o: $(SRCROOT)/MOM6/src/ice_shelf/MOM_ice_shelf_state.F90 MOM_cpu_clock.o MOM_dyn_horgrid.o MOM_error_handler.o MOM_file_parser.o MOM_grid.o MOM_get_input.o MOM_coms.o MOM_checksums.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)/MOM6/src/ice_shelf/MOM_ice_shelf_state.F90
MOM_interface_filter.o: $(SRCROOT)/MOM6/src/parameterizations/lateral/MOM_interface_filter.F90 $(SRCROOT)/MOM6/config_src/memory/dynamic_nonsymmetric/MOM_memory.h $(SRCROOT)/MOM6/src/framework/MOM_memory_macros.h $(SRCROOT)/MOM6/src/framework/version_variable.h MOM_debugging.o MOM_diag_mediator.o MOM_domains.o MOM_error_handler.o MOM_file_parser.o MOM_grid.o MOM_interface_heights.o MOM_unit_scaling.o MOM_variables.o MOM_verticalGrid.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)/MOM6/config_src/memory/dynamic_nonsymmetric -I$(SRCROOT)/MOM6/src/framework	$(SRCROOT)/MOM6/src/parameterizations/lateral/MOM_interface_filter.F90
MOM_interface_heights.o: $(SRCROOT)/MOM6/src/core/MOM_interface_heights.F90 $(SRCROOT)/MOM6/config_src/memory/dynamic_nonsymmetric/MOM_memory.h $(SRCROOT)/MOM6/src/framework/MOM_memory_macros.h MOM_density_integrals.o MOM_debugging.o MOM_error_handler.o MOM_EOS.o MOM_file_parser.o MOM_grid.o MOM_unit_scaling.o MOM_variables.o MOM_verticalGrid.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)/MOM6/config_src/memory/dynamic_nonsymmetric -I$(SRCROOT)/MOM6/src/framework	$(SRCROOT)/MOM6/src/core/MOM_interface_heights.F90
MOM_internal_tide_input.o: $(SRCROOT)/MOM6/src/parameterizations/vertical/MOM_internal_tide_input.F90 $(SRCROOT)/MOM6/config_src/memory/dynamic_nonsymmetric/MOM_memory.h $(SRCROOT)/MOM6/src/framework/MOM_memory_macros.h $(SRCROOT)/MOM6/src/framework/version_variable.h MOM_cpu_clock.o MOM_diag_mediator.o MOM_debugging.o MOM_error_handler.o MOM_file_parser.o MOM_forcing_type.o MOM_grid.o MOM_io.o MOM_isopycnal_slopes.o MOM_time_manager.o MOM_unit_scaling.o MOM_variables.o MOM_verticalGrid.o MOM_EOS.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)/MOM6/config_src/memory/dynamic_nonsymmetric -I$(SRCROOT)/MOM6/src/framework	$(SRCROOT)/MOM6/src/parameterizations/vertical/MOM_internal_tide_input.F90
MOM_internal_tides.o: $(SRCROOT)/MOM6/src/parameterizations/lateral/MOM_internal_tides.F90 $(SRCROOT)/MOM6/config_src/memory/dynamic_nonsymmetric/MOM_memory.h $(SRCROOT)/MOM6/src/framework/MOM_memory_macros.h $(SRCROOT)/MOM6/src/framework/version_variable.h MOM_debugging.o MOM_diag_mediator.o MOM_domains.o MOM_error_handler.o MOM_file_parser.o MOM_grid.o MOM_io.o MOM_restart.o MOM_spatial_means.o MOM_time_manager.o MOM_unit_scaling.o MOM_variables.o MOM_verticalGrid.o MOM_wave_speed.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)/MOM6/config_src/memory/dynamic_nonsymmetric -I$(SRCROOT)/MOM6/src/framework	$(SRCROOT)/MOM6/src/parameterizations/lateral/MOM_internal_tides.F90
MOM_interp_infra.o: $(SRCROOT)/MOM6/config_src/infra/FMS2/MOM_interp_infra.F90 MOM_domain_infra.o MOM_io.o MOM_time_manager.o MOM_error_handler.o MOM_string_functions.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)/MOM6/config_src/infra/FMS2/MOM_interp_infra.F90
MOM_interpolate.o: $(SRCROOT)/MOM6/src/framework/MOM_interpolate.F90 MOM_array_transform.o MOM_error_handler.o MOM_interp_infra.o MOM_time_manager.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)/MOM6/src/framework/MOM_interpolate.F90
MOM_intrinsic_functions.o: $(SRCROOT)/MOM6/src/framework/MOM_intrinsic_functions.F90
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)/MOM6/src/framework/MOM_intrinsic_functions.F90
MOM_io.o: $(SRCROOT)/MOM6/src/framework/MOM_io.F90 $(SRCROOT)/MOM6/src/framework/version_variable.h MOM_array_transform.o MOM_domains.o MOM_dyn_horgrid.o MOM_ensemble_manager.o MOM_error_handler.o MOM_file_parser.o MOM_grid.o MOM_io_infra.o MOM_io_file.o MOM_string_functions.o MOM_verticalGrid.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)/MOM6/src/framework	$(SRCROOT)/MOM6/src/framework/MOM_io.F90
MOM_io_file.o: $(SRCROOT)/MOM6/src/framework/MOM_io_file.F90 MOM_domains.o MOM_io_infra.o MOM_hor_index.o MOM_netcdf.o MOM_error_handler.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)/MOM6/src/framework/MOM_io_file.F90
MOM_io_infra.o: $(SRCROOT)/MOM6/config_src/infra/FMS2/MOM_io_infra.F90 MOM_domain_infra.o MOM_error_infra.o MOM_string_functions.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)/MOM6/config_src/infra/FMS2/MOM_io_infra.F90
MOM_isopycnal_slopes.o: $(SRCROOT)/MOM6/src/core/MOM_isopycnal_slopes.F90 $(SRCROOT)/MOM6/config_src/memory/dynamic_nonsymmetric/MOM_memory.h $(SRCROOT)/MOM6/src/framework/MOM_memory_macros.h MOM_debugging.o MOM_grid.o MOM_unit_scaling.o MOM_variables.o MOM_verticalGrid.o MOM_EOS.o MOM_open_boundary.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)/MOM6/config_src/memory/dynamic_nonsymmetric -I$(SRCROOT)/MOM6/src/framework	$(SRCROOT)/MOM6/src/core/MOM_isopycnal_slopes.F90
MOM_kappa_shear.o: $(SRCROOT)/MOM6/src/parameterizations/vertical/MOM_kappa_shear.F90 $(SRCROOT)/MOM6/config_src/memory/dynamic_nonsymmetric/MOM_memory.h $(SRCROOT)/MOM6/src/framework/MOM_memory_macros.h $(SRCROOT)/MOM6/src/framework/version_variable.h MOM_cpu_clock.o MOM_diag_mediator.o MOM_debugging.o MOM_error_handler.o MOM_file_parser.o MOM_grid.o MOM_unit_scaling.o MOM_variables.o MOM_verticalGrid.o MOM_EOS.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)/MOM6/config_src/memory/dynamic_nonsymmetric -I$(SRCROOT)/MOM6/src/framework	$(SRCROOT)/MOM6/src/parameterizations/vertical/MOM_kappa_shear.F90
MOM_lateral_boundary_diffusion.o: $(SRCROOT)/MOM6/src/tracer/MOM_lateral_boundary_diffusion.F90 $(SRCROOT)/MOM6/config_src/memory/dynamic_nonsymmetric/MOM_memory.h $(SRCROOT)/MOM6/src/framework/MOM_memory_macros.h $(SRCROOT)/MOM6/src/framework/version_variable.h MOM_cpu_clock.o MOM_checksums.o MOM_domains.o MOM_diag_mediator.o MOM_error_handler.o MOM_file_parser.o MOM_grid.o MOM_remapping.o MOM_spatial_means.o MOM_tracer_registry.o MOM_unit_scaling.o MOM_verticalGrid.o MOM_CVMix_KPP.o MOM_energetic_PBL.o MOM_diabatic_driver.o MOM_io.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)/MOM6/config_src/memory/dynamic_nonsymmetric -I$(SRCROOT)/MOM6/src/framework	$(SRCROOT)/MOM6/src/tracer/MOM_lateral_boundary_diffusion.F90
MOM_lateral_mixing_coeffs.o: $(SRCROOT)/MOM6/src/parameterizations/lateral/MOM_lateral_mixing_coeffs.F90 $(SRCROOT)/MOM6/config_src/memory/dynamic_nonsymmetric/MOM_memory.h $(SRCROOT)/MOM6/src/framework/MOM_memory_macros.h $(SRCROOT)/MOM6/src/framework/version_variable.h MOM_debugging.o MOM_error_handler.o MOM_diag_mediator.o MOM_domains.o MOM_file_parser.o MOM_interface_heights.o MOM_isopycnal_slopes.o MOM_grid.o MOM_unit_scaling.o MOM_variables.o MOM_verticalGrid.o MOM_wave_speed.o MOM_open_boundary.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)/MOM6/config_src/memory/dynamic_nonsymmetric -I$(SRCROOT)/MOM6/src/framework	$(SRCROOT)/MOM6/src/parameterizations/lateral/MOM_lateral_mixing_coeffs.F90
MOM_load_love_numbers.o: $(SRCROOT)/MOM6/src/parameterizations/lateral/MOM_load_love_numbers.F90
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)/MOM6/src/parameterizations/lateral/MOM_load_love_numbers.F90
MOM_marine_ice.o: $(SRCROOT)/MOM6/src/ice_shelf/MOM_marine_ice.F90 $(SRCROOT)/MOM6/config_src/memory/dynamic_nonsymmetric/MOM_memory.h $(SRCROOT)/MOM6/src/framework/MOM_memory_macros.h $(SRCROOT)/MOM6/src/framework/version_variable.h MOM_constants.o MOM_diag_mediator.o MOM_domains.o MOM_error_handler.o MOM_file_parser.o MOM_forcing_type.o MOM_grid.o MOM_time_manager.o MOM_unit_scaling.o MOM_variables.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)/MOM6/config_src/memory/dynamic_nonsymmetric -I$(SRCROOT)/MOM6/src/framework	$(SRCROOT)/MOM6/src/ice_shelf/MOM_marine_ice.F90
MOM_mixed_layer_restrat.o: $(SRCROOT)/MOM6/src/parameterizations/lateral/MOM_mixed_layer_restrat.F90 $(SRCROOT)/MOM6/config_src/memory/dynamic_nonsymmetric/MOM_memory.h $(SRCROOT)/MOM6/src/framework/MOM_memory_macros.h $(SRCROOT)/MOM6/src/framework/version_variable.h MOM_debugging.o MOM_diag_mediator.o MOM_domains.o MOM_error_handler.o MOM_file_parser.o MOM_forcing_type.o MOM_grid.o MOM_hor_index.o MOM_lateral_mixing_coeffs.o MOM_restart.o MOM_unit_scaling.o MOM_variables.o MOM_verticalGrid.o MOM_EOS.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)/MOM6/config_src/memory/dynamic_nonsymmetric -I$(SRCROOT)/MOM6/src/framework	$(SRCROOT)/MOM6/src/parameterizations/lateral/MOM_mixed_layer_restrat.F90
MOM_netcdf.o: $(SRCROOT)/MOM6/src/framework/MOM_netcdf.F90 MOM_error_handler.o MOM_io_infra.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)/MOM6/src/framework/MOM_netcdf.F90
MOM_neutral_diffusion.o: $(SRCROOT)/MOM6/src/tracer/MOM_neutral_diffusion.F90 $(SRCROOT)/MOM6/config_src/memory/dynamic_nonsymmetric/MOM_memory.h $(SRCROOT)/MOM6/src/framework/MOM_memory_macros.h $(SRCROOT)/MOM6/src/framework/version_variable.h MOM_cpu_clock.o MOM_domains.o MOM_diag_mediator.o MOM_EOS.o MOM_error_handler.o MOM_file_parser.o MOM_grid.o MOM_remapping.o MOM_tracer_registry.o MOM_unit_scaling.o MOM_verticalGrid.o polynomial_functions.o PPM_functions.o regrid_edge_values.o MOM_CVMix_KPP.o MOM_energetic_PBL.o MOM_diabatic_driver.o MOM_lateral_boundary_diffusion.o MOM_io.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)/MOM6/config_src/memory/dynamic_nonsymmetric -I$(SRCROOT)/MOM6/src/framework	$(SRCROOT)/MOM6/src/tracer/MOM_neutral_diffusion.F90
MOM_obsolete_diagnostics.o: $(SRCROOT)/MOM6/src/diagnostics/MOM_obsolete_diagnostics.F90 $(SRCROOT)/MOM6/config_src/memory/dynamic_nonsymmetric/MOM_memory.h $(SRCROOT)/MOM6/src/framework/MOM_memory_macros.h $(SRCROOT)/MOM6/src/framework/version_variable.h MOM_diag_mediator.o MOM_error_handler.o MOM_file_parser.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)/MOM6/config_src/memory/dynamic_nonsymmetric -I$(SRCROOT)/MOM6/src/framework	$(SRCROOT)/MOM6/src/diagnostics/MOM_obsolete_diagnostics.F90
MOM_obsolete_params.o: $(SRCROOT)/MOM6/src/diagnostics/MOM_obsolete_params.F90 $(SRCROOT)/MOM6/config_src/memory/dynamic_nonsymmetric/MOM_memory.h $(SRCROOT)/MOM6/src/framework/MOM_memory_macros.h $(SRCROOT)/MOM6/src/framework/version_variable.h MOM_error_handler.o MOM_file_parser.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)/MOM6/config_src/memory/dynamic_nonsymmetric -I$(SRCROOT)/MOM6/src/framework	$(SRCROOT)/MOM6/src/diagnostics/MOM_obsolete_params.F90
MOM_oda_driver.o: $(SRCROOT)/MOM6/src/ocean_data_assim/MOM_oda_driver.F90 $(SRCROOT)/MOM6/config_src/memory/dynamic_nonsymmetric/MOM_memory.h $(SRCROOT)/MOM6/src/framework/MOM_memory_macros.h MOM_coms.o MOM_domains.o MOM_diag_mediator.o MOM_ensemble_manager.o MOM_error_handler.o MOM_io.o MOM_interp_infra.o MOM_interpolate.o MOM_remapping.o MOM_time_manager.o MOM_cpu_clock.o MOM_horizontal_regridding.o ocean_da_types.o ocean_da_core.o kdtree.o MOM_get_input.o MOM_grid.o MOM_grid_initialize.o MOM_hor_index.o MOM_dyn_horgrid.o MOM_transcribe_grid.o MOM_fixed_initialization.o MOM_coord_initialization.o MOM_file_parser.o MOM_string_functions.o MOM_ALE.o MOM_regridding.o MOM_unit_scaling.o MOM_variables.o MOM_verticalGrid.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)/MOM6/config_src/memory/dynamic_nonsymmetric -I$(SRCROOT)/MOM6/src/framework	$(SRCROOT)/MOM6/src/ocean_data_assim/MOM_oda_driver.F90
MOM_oda_incupd.o: $(SRCROOT)/MOM6/src/ocean_data_assim/MOM_oda_incupd.F90 $(SRCROOT)/MOM6/config_src/memory/dynamic_nonsymmetric/MOM_memory.h $(SRCROOT)/MOM6/src/framework/MOM_memory_macros.h $(SRCROOT)/MOM6/src/framework/version_variable.h MOM_array_transform.o MOM_coms.o MOM_diag_mediator.o MOM_domains.o MOM_error_handler.o MOM_file_parser.o MOM_get_input.o MOM_grid.o MOM_io.o MOM_remapping.o MOM_restart.o MOM_spatial_means.o MOM_time_manager.o MOM_unit_scaling.o MOM_variables.o MOM_verticalGrid.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)/MOM6/config_src/memory/dynamic_nonsymmetric -I$(SRCROOT)/MOM6/src/framework	$(SRCROOT)/MOM6/src/ocean_data_assim/MOM_oda_incupd.F90
MOM_offline_aux.o: $(SRCROOT)/MOM6/src/tracer/MOM_offline_aux.F90 $(SRCROOT)/MOM6/config_src/memory/dynamic_nonsymmetric/MOM_memory.h $(SRCROOT)/MOM6/src/framework/MOM_memory_macros.h MOM_debugging.o MOM_domains.o MOM_diag_mediator.o MOM_error_handler.o MOM_file_parser.o MOM_forcing_type.o MOM_grid.o MOM_io.o MOM_opacity.o MOM_time_manager.o MOM_unit_scaling.o MOM_variables.o MOM_verticalGrid.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)/MOM6/config_src/memory/dynamic_nonsymmetric -I$(SRCROOT)/MOM6/src/framework	$(SRCROOT)/MOM6/src/tracer/MOM_offline_aux.F90
MOM_offline_main.o: $(SRCROOT)/MOM6/src/tracer/MOM_offline_main.F90 $(SRCROOT)/MOM6/config_src/memory/dynamic_nonsymmetric/MOM_memory.h $(SRCROOT)/MOM6/src/framework/MOM_memory_macros.h $(SRCROOT)/MOM6/src/framework/version_variable.h MOM_ALE.o MOM_checksums.o MOM_coms.o MOM_cpu_clock.o MOM_diabatic_aux.o MOM_diabatic_driver.o MOM_diag_mediator.o MOM_domains.o MOM_error_handler.o MOM_file_parser.o MOM_forcing_type.o MOM_grid.o MOM_interface_heights.o MOM_io.o MOM_offline_aux.o MOM_opacity.o MOM_open_boundary.o MOM_time_manager.o MOM_tracer_advect.o MOM_tracer_diabatic.o MOM_tracer_flow_control.o MOM_tracer_registry.o MOM_unit_scaling.o MOM_variables.o MOM_verticalGrid.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)/MOM6/config_src/memory/dynamic_nonsymmetric -I$(SRCROOT)/MOM6/src/framework	$(SRCROOT)/MOM6/src/tracer/MOM_offline_main.F90
MOM_opacity.o: $(SRCROOT)/MOM6/src/parameterizations/vertical/MOM_opacity.F90 $(SRCROOT)/MOM6/config_src/memory/dynamic_nonsymmetric/MOM_memory.h $(SRCROOT)/MOM6/src/framework/MOM_memory_macros.h $(SRCROOT)/MOM6/src/framework/version_variable.h MOM_diag_mediator.o MOM_error_handler.o MOM_file_parser.o MOM_string_functions.o MOM_grid.o MOM_unit_scaling.o MOM_variables.o MOM_verticalGrid.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)/MOM6/config_src/memory/dynamic_nonsymmetric -I$(SRCROOT)/MOM6/src/framework	$(SRCROOT)/MOM6/src/parameterizations/vertical/MOM_opacity.F90
MOM_open_boundary.o: $(SRCROOT)/MOM6/src/core/MOM_open_boundary.F90 $(SRCROOT)/MOM6/config_src/memory/dynamic_nonsymmetric/MOM_memory.h $(SRCROOT)/MOM6/src/framework/MOM_memory_macros.h $(SRCROOT)/MOM6/src/framework/version_variable.h $(SRCROOT)MOM6/src/framework/version_variable.h MOM_array_transform.o MOM_coms.o MOM_cpu_clock.o MOM_diag_mediator.o MOM_domains.o MOM_error_handler.o MOM_file_parser.o MOM_grid.o MOM_dyn_horgrid.o MOM_io.o MOM_restart.o MOM_obsolete_params.o MOM_string_functions.o MOM_tidal_forcing.o MOM_time_manager.o MOM_tracer_registry.o MOM_interpolate.o MOM_remapping.o MOM_regridding.o MOM_unit_scaling.o MOM_variables.o MOM_verticalGrid.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)/MOM6/config_src/memory/dynamic_nonsymmetric -I$(SRCROOT)/MOM6/src/framework -I$(SRCROOT)MOM6/src/framework	$(SRCROOT)/MOM6/src/core/MOM_open_boundary.F90
MOM_particles.o: $(SRCROOT)/MOM6/config_src/external/drifters/MOM_particles.F90 MOM_grid.o MOM_time_manager.o MOM_variables.o MOM_particles_types.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)/MOM6/config_src/external/drifters/MOM_particles.F90
MOM_particles_types.o: $(SRCROOT)/MOM6/config_src/external/drifters/MOM_particles_types.F90 MOM_grid.o MOM_domains.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)/MOM6/config_src/external/drifters/MOM_particles_types.F90
MOM_porous_barriers.o: $(SRCROOT)/MOM6/src/core/MOM_porous_barriers.F90 $(SRCROOT)/MOM6/config_src/memory/dynamic_nonsymmetric/MOM_memory.h $(SRCROOT)/MOM6/src/framework/MOM_memory_macros.h $(SRCROOT)/MOM6/src/framework/version_variable.h MOM_cpu_clock.o MOM_error_handler.o MOM_grid.o MOM_unit_scaling.o MOM_variables.o MOM_verticalGrid.o MOM_interface_heights.o MOM_time_manager.o MOM_diag_mediator.o MOM_file_parser.o MOM_debugging.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)/MOM6/config_src/memory/dynamic_nonsymmetric -I$(SRCROOT)/MOM6/src/framework	$(SRCROOT)/MOM6/src/core/MOM_porous_barriers.F90
MOM_random.o: $(SRCROOT)/MOM6/src/framework/MOM_random.F90 MOM_hor_index.o MOM_time_manager.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)/MOM6/src/framework/MOM_random.F90
MOM_regridding.o: $(SRCROOT)/MOM6/src/ALE/MOM_regridding.F90 $(SRCROOT)/MOM6/config_src/memory/dynamic_nonsymmetric/MOM_memory.h $(SRCROOT)/MOM6/src/framework/MOM_memory_macros.h MOM_error_handler.o MOM_file_parser.o MOM_io.o MOM_unit_scaling.o MOM_variables.o MOM_verticalGrid.o MOM_EOS.o MOM_string_functions.o MOM_remapping.o regrid_consts.o regrid_interp.o coord_zlike.o coord_sigma.o coord_rho.o coord_hycom.o coord_adapt.o MOM_hybgen_regrid.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)/MOM6/config_src/memory/dynamic_nonsymmetric -I$(SRCROOT)/MOM6/src/framework	$(SRCROOT)/MOM6/src/ALE/MOM_regridding.F90
MOM_regularize_layers.o: $(SRCROOT)/MOM6/src/parameterizations/vertical/MOM_regularize_layers.F90 $(SRCROOT)/MOM6/config_src/memory/dynamic_nonsymmetric/MOM_memory.h $(SRCROOT)/MOM6/src/framework/MOM_memory_macros.h $(SRCROOT)/MOM6/src/framework/version_variable.h MOM_cpu_clock.o MOM_diag_mediator.o MOM_domains.o MOM_error_handler.o MOM_file_parser.o MOM_grid.o MOM_unit_scaling.o MOM_variables.o MOM_verticalGrid.o MOM_EOS.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)/MOM6/config_src/memory/dynamic_nonsymmetric -I$(SRCROOT)/MOM6/src/framework	$(SRCROOT)/MOM6/src/parameterizations/vertical/MOM_regularize_layers.F90
MOM_remapping.o: $(SRCROOT)/MOM6/src/ALE/MOM_remapping.F90 MOM_error_handler.o MOM_io.o MOM_string_functions.o regrid_edge_values.o remapping_attic.o PCM_functions.o PLM_functions.o PPM_functions.o PQM_functions.o MOM_hybgen_remap.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)/MOM6/src/ALE/MOM_remapping.F90
MOM_restart.o: $(SRCROOT)/MOM6/src/framework/MOM_restart.F90 $(SRCROOT)/MOM6/src/framework/version_variable.h MOM_checksums.o MOM_domains.o MOM_error_handler.o MOM_file_parser.o MOM_grid.o MOM_io.o MOM_string_functions.o MOM_time_manager.o MOM_verticalGrid.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)/MOM6/src/framework	$(SRCROOT)/MOM6/src/framework/MOM_restart.F90
MOM_safe_alloc.o: $(SRCROOT)/MOM6/src/framework/MOM_safe_alloc.F90
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)/MOM6/src/framework/MOM_safe_alloc.F90
MOM_set_diffusivity.o: $(SRCROOT)/MOM6/src/parameterizations/vertical/MOM_set_diffusivity.F90 $(SRCROOT)/MOM6/config_src/memory/dynamic_nonsymmetric/MOM_memory.h $(SRCROOT)/MOM6/src/framework/MOM_memory_macros.h $(SRCROOT)/MOM6/src/framework/version_variable.h MOM_bkgnd_mixing.o MOM_cpu_clock.o MOM_CVMix_ddiff.o MOM_CVMix_shear.o MOM_diag_mediator.o MOM_debugging.o MOM_EOS.o MOM_error_handler.o MOM_file_parser.o MOM_forcing_type.o MOM_full_convection.o MOM_grid.o MOM_internal_tides.o MOM_intrinsic_functions.o MOM_io.o MOM_isopycnal_slopes.o MOM_kappa_shear.o MOM_open_boundary.o MOM_string_functions.o MOM_tidal_mixing.o MOM_unit_scaling.o MOM_variables.o MOM_verticalGrid.o user_change_diffusivity.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)/MOM6/config_src/memory/dynamic_nonsymmetric -I$(SRCROOT)/MOM6/src/framework	$(SRCROOT)/MOM6/src/parameterizations/vertical/MOM_set_diffusivity.F90
MOM_set_viscosity.o: $(SRCROOT)/MOM6/src/parameterizations/vertical/MOM_set_viscosity.F90 $(SRCROOT)/MOM6/config_src/memory/dynamic_nonsymmetric/MOM_memory.h $(SRCROOT)/MOM6/src/framework/MOM_memory_macros.h $(SRCROOT)/MOM6/src/framework/version_variable.h MOM_ALE.o MOM_cpu_clock.o MOM_debugging.o MOM_diag_mediator.o MOM_domains.o MOM_error_handler.o MOM_file_parser.o MOM_forcing_type.o MOM_grid.o MOM_hor_index.o MOM_io.o MOM_kappa_shear.o MOM_CVMix_shear.o MOM_CVMix_conv.o MOM_CVMix_ddiff.o MOM_restart.o MOM_safe_alloc.o MOM_unit_scaling.o MOM_variables.o MOM_verticalGrid.o MOM_EOS.o MOM_open_boundary.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)/MOM6/config_src/memory/dynamic_nonsymmetric -I$(SRCROOT)/MOM6/src/framework	$(SRCROOT)/MOM6/src/parameterizations/vertical/MOM_set_viscosity.F90
MOM_shared_initialization.o: $(SRCROOT)/MOM6/src/initialization/MOM_shared_initialization.F90 $(SRCROOT)/MOM6/src/framework/version_variable.h MOM_coms.o MOM_domains.o MOM_dyn_horgrid.o MOM_error_handler.o MOM_file_parser.o MOM_io.o MOM_string_functions.o MOM_unit_scaling.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)/MOM6/src/framework	$(SRCROOT)/MOM6/src/initialization/MOM_shared_initialization.F90
MOM_spatial_means.o: $(SRCROOT)/MOM6/src/diagnostics/MOM_spatial_means.F90 $(SRCROOT)/MOM6/config_src/memory/dynamic_nonsymmetric/MOM_memory.h $(SRCROOT)/MOM6/src/framework/MOM_memory_macros.h MOM_coms.o MOM_error_handler.o MOM_file_parser.o MOM_grid.o MOM_verticalGrid.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)/MOM6/config_src/memory/dynamic_nonsymmetric -I$(SRCROOT)/MOM6/src/framework	$(SRCROOT)/MOM6/src/diagnostics/MOM_spatial_means.F90
MOM_spherical_harmonics.o: $(SRCROOT)/MOM6/src/parameterizations/lateral/MOM_spherical_harmonics.F90 $(SRCROOT)/MOM6/config_src/memory/dynamic_nonsymmetric/MOM_memory.h $(SRCROOT)/MOM6/src/framework/MOM_memory_macros.h $(SRCROOT)/MOM6/src/framework/version_variable.h MOM_cpu_clock.o MOM_error_handler.o MOM_file_parser.o MOM_grid.o MOM_coms_infra.o MOM_coms.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)/MOM6/config_src/memory/dynamic_nonsymmetric -I$(SRCROOT)/MOM6/src/framework	$(SRCROOT)/MOM6/src/parameterizations/lateral/MOM_spherical_harmonics.F90
MOM_sponge.o: $(SRCROOT)/MOM6/src/parameterizations/vertical/MOM_sponge.F90 $(SRCROOT)/MOM6/config_src/memory/dynamic_nonsymmetric/MOM_memory.h $(SRCROOT)/MOM6/src/framework/MOM_memory_macros.h $(SRCROOT)/MOM6/src/framework/version_variable.h MOM_coms.o MOM_diag_mediator.o MOM_error_handler.o MOM_file_parser.o MOM_grid.o MOM_spatial_means.o MOM_time_manager.o MOM_unit_scaling.o MOM_variables.o MOM_verticalGrid.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)/MOM6/config_src/memory/dynamic_nonsymmetric -I$(SRCROOT)/MOM6/src/framework	$(SRCROOT)/MOM6/src/parameterizations/vertical/MOM_sponge.F90
MOM_state_initialization.o: $(SRCROOT)/MOM6/src/initialization/MOM_state_initialization.F90 $(SRCROOT)/MOM6/config_src/memory/dynamic_nonsymmetric/MOM_memory.h $(SRCROOT)/MOM6/src/framework/MOM_memory_macros.h $(SRCROOT)/MOM6/src/framework/version_variable.h $(SRCROOT)MOM6/src/framework/version_variable.h MOM_debugging.o MOM_density_integrals.o MOM_coms.o MOM_cpu_clock.o MOM_domains.o MOM_error_handler.o MOM_file_parser.o MOM_get_input.o MOM_grid.o MOM_interface_heights.o MOM_io.o MOM_open_boundary.o MOM_grid_initialize.o MOM_restart.o MOM_sponge.o MOM_ALE_sponge.o MOM_string_functions.o MOM_time_manager.o MOM_tracer_registry.o MOM_unit_scaling.o MOM_variables.o MOM_verticalGrid.o MOM_EOS.o user_initialization.o DOME_initialization.o ISOMIP_initialization.o RGC_initialization.o baroclinic_zone_initialization.o benchmark_initialization.o Neverworld_initialization.o circle_obcs_initialization.o lock_exchange_initialization.o external_gwave_initialization.o DOME2d_initialization.o adjustment_initialization.o sloshing_initialization.o seamount_initialization.o dumbbell_initialization.o Phillips_initialization.o Rossby_front_2d_initialization.o SCM_CVMix_tests.o dyed_channel_initialization.o dyed_obcs_initialization.o supercritical_initialization.o soliton_initialization.o BFB_initialization.o dense_water_initialization.o MOM_tracer_Z_init.o MOM_ALE.o MOM_regridding.o MOM_remapping.o MOM_horizontal_regridding.o MOM_oda_incupd.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)/MOM6/config_src/memory/dynamic_nonsymmetric -I$(SRCROOT)/MOM6/src/framework -I$(SRCROOT)MOM6/src/framework	$(SRCROOT)/MOM6/src/initialization/MOM_state_initialization.F90
MOM_stoch_eos.o: $(SRCROOT)/MOM6/src/core/MOM_stoch_eos.F90 $(SRCROOT)/MOM6/config_src/memory/dynamic_nonsymmetric/MOM_memory.h $(SRCROOT)/MOM6/src/framework/MOM_memory_macros.h MOM_diag_mediator.o MOM_error_handler.o MOM_file_parser.o MOM_grid.o MOM_hor_index.o MOM_isopycnal_slopes.o MOM_random.o MOM_restart.o MOM_time_manager.o MOM_unit_scaling.o MOM_variables.o MOM_verticalGrid.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)/MOM6/config_src/memory/dynamic_nonsymmetric -I$(SRCROOT)/MOM6/src/framework	$(SRCROOT)/MOM6/src/core/MOM_stoch_eos.F90
MOM_stochastics.o: $(SRCROOT)/MOM6/src/parameterizations/stochastic/MOM_stochastics.F90 $(SRCROOT)/MOM6/config_src/memory/dynamic_nonsymmetric/MOM_memory.h $(SRCROOT)/MOM6/src/framework/MOM_memory_macros.h $(SRCROOT)/MOM6/src/framework/version_variable.h MOM_diag_mediator.o MOM_grid.o MOM_verticalGrid.o MOM_error_handler.o MOM_file_parser.o MOM_domains.o MOM_coms.o stochastic_physics.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)/MOM6/config_src/memory/dynamic_nonsymmetric -I$(SRCROOT)/MOM6/src/framework	$(SRCROOT)/MOM6/src/parameterizations/stochastic/MOM_stochastics.F90
MOM_string_functions.o: $(SRCROOT)/MOM6/src/framework/MOM_string_functions.F90
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)/MOM6/src/framework/MOM_string_functions.F90
MOM_sum_output.o: $(SRCROOT)/MOM6/src/diagnostics/MOM_sum_output.F90 $(SRCROOT)/MOM6/config_src/memory/dynamic_nonsymmetric/MOM_memory.h $(SRCROOT)/MOM6/src/framework/MOM_memory_macros.h $(SRCROOT)/MOM6/src/framework/version_variable.h MOM_checksums.o MOM_coms.o MOM_error_handler.o MOM_file_parser.o MOM_forcing_type.o MOM_grid.o MOM_interface_heights.o MOM_io.o MOM_time_manager.o MOM_tracer_flow_control.o MOM_unit_scaling.o MOM_variables.o MOM_verticalGrid.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)/MOM6/config_src/memory/dynamic_nonsymmetric -I$(SRCROOT)/MOM6/src/framework	$(SRCROOT)/MOM6/src/diagnostics/MOM_sum_output.F90
MOM_surface_forcing_gfdl.o: $(SRCROOT)/MOM6/config_src/drivers/FMS_cap/MOM_surface_forcing_gfdl.F90 $(SRCROOT)/MOM6/config_src/memory/dynamic_nonsymmetric/MOM_memory.h $(SRCROOT)/MOM6/src/framework/MOM_memory_macros.h $(SRCROOT)/MOM6/src/framework/version_variable.h MOM_coms.o MOM_constants.o MOM_coupler_types.o MOM_cpu_clock.o MOM_data_override.o MOM_diag_mediator.o MOM_domains.o MOM_error_handler.o MOM_file_parser.o MOM_forcing_type.o MOM_get_input.o MOM_grid.o MOM_interpolate.o MOM_io.o MOM_restart.o MOM_string_functions.o MOM_spatial_means.o MOM_unit_scaling.o MOM_variables.o user_revise_forcing.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)/MOM6/config_src/memory/dynamic_nonsymmetric -I$(SRCROOT)/MOM6/src/framework	$(SRCROOT)/MOM6/config_src/drivers/FMS_cap/MOM_surface_forcing_gfdl.F90
MOM_temperature_convert.o: $(SRCROOT)/MOM6/src/equation_of_state/MOM_temperature_convert.F90
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)/MOM6/src/equation_of_state/MOM_temperature_convert.F90
MOM_thickness_diffuse.o: $(SRCROOT)/MOM6/src/parameterizations/lateral/MOM_thickness_diffuse.F90 $(SRCROOT)/MOM6/config_src/memory/dynamic_nonsymmetric/MOM_memory.h $(SRCROOT)/MOM6/src/framework/MOM_memory_macros.h $(SRCROOT)/MOM6/src/framework/version_variable.h MOM_debugging.o MOM_diag_mediator.o MOM_domains.o MOM_error_handler.o MOM_EOS.o MOM_file_parser.o MOM_grid.o MOM_io.o MOM_interface_heights.o MOM_isopycnal_slopes.o MOM_lateral_mixing_coeffs.o MOM_MEKE_types.o MOM_unit_scaling.o MOM_variables.o MOM_verticalGrid.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)/MOM6/config_src/memory/dynamic_nonsymmetric -I$(SRCROOT)/MOM6/src/framework	$(SRCROOT)/MOM6/src/parameterizations/lateral/MOM_thickness_diffuse.F90
MOM_tidal_forcing.o: $(SRCROOT)/MOM6/src/parameterizations/lateral/MOM_tidal_forcing.F90 $(SRCROOT)/MOM6/config_src/memory/dynamic_nonsymmetric/MOM_memory.h $(SRCROOT)/MOM6/src/framework/MOM_memory_macros.h $(SRCROOT)/MOM6/src/framework/version_variable.h MOM_cpu_clock.o MOM_domains.o MOM_error_handler.o MOM_file_parser.o MOM_grid.o MOM_io.o MOM_time_manager.o MOM_unit_scaling.o MOM_spherical_harmonics.o MOM_load_love_numbers.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)/MOM6/config_src/memory/dynamic_nonsymmetric -I$(SRCROOT)/MOM6/src/framework	$(SRCROOT)/MOM6/src/parameterizations/lateral/MOM_tidal_forcing.F90
MOM_tidal_mixing.o: $(SRCROOT)/MOM6/src/parameterizations/vertical/MOM_tidal_mixing.F90 $(SRCROOT)/MOM6/config_src/memory/dynamic_nonsymmetric/MOM_memory.h $(SRCROOT)/MOM6/src/framework/MOM_memory_macros.h $(SRCROOT)/MOM6/src/framework/version_variable.h MOM_diag_mediator.o MOM_debugging.o MOM_error_handler.o MOM_file_parser.o MOM_grid.o MOM_io.o MOM_internal_tides.o MOM_remapping.o MOM_string_functions.o MOM_unit_scaling.o MOM_variables.o MOM_verticalGrid.o cvmix_tidal.o cvmix_kinds_and_types.o cvmix_put_get.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)/MOM6/config_src/memory/dynamic_nonsymmetric -I$(SRCROOT)/MOM6/src/framework	$(SRCROOT)/MOM6/src/parameterizations/vertical/MOM_tidal_mixing.F90
MOM_time_manager.o: $(SRCROOT)/MOM6/config_src/infra/FMS2/MOM_time_manager.F90
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)/MOM6/config_src/infra/FMS2/MOM_time_manager.F90
MOM_tracer_Z_init.o: $(SRCROOT)/MOM6/src/tracer/MOM_tracer_Z_init.F90 $(SRCROOT)/MOM6/config_src/memory/dynamic_nonsymmetric/MOM_memory.h $(SRCROOT)/MOM6/src/framework/MOM_memory_macros.h $(SRCROOT)/MOM6/src/framework/version_variable.h $(SRCROOT)MOM6/src/framework/version_variable.h MOM_error_handler.o MOM_file_parser.o MOM_grid.o MOM_io.o MOM_EOS.o MOM_unit_scaling.o MOM_verticalGrid.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)/MOM6/config_src/memory/dynamic_nonsymmetric -I$(SRCROOT)/MOM6/src/framework -I$(SRCROOT)MOM6/src/framework	$(SRCROOT)/MOM6/src/tracer/MOM_tracer_Z_init.F90
MOM_tracer_advect.o: $(SRCROOT)/MOM6/src/tracer/MOM_tracer_advect.F90 $(SRCROOT)/MOM6/config_src/memory/dynamic_nonsymmetric/MOM_memory.h $(SRCROOT)/MOM6/src/framework/MOM_memory_macros.h $(SRCROOT)/MOM6/src/framework/version_variable.h MOM_cpu_clock.o MOM_diag_mediator.o MOM_domains.o MOM_error_handler.o MOM_file_parser.o MOM_grid.o MOM_open_boundary.o MOM_tracer_registry.o MOM_unit_scaling.o MOM_verticalGrid.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)/MOM6/config_src/memory/dynamic_nonsymmetric -I$(SRCROOT)/MOM6/src/framework	$(SRCROOT)/MOM6/src/tracer/MOM_tracer_advect.F90
MOM_tracer_diabatic.o: $(SRCROOT)/MOM6/src/tracer/MOM_tracer_diabatic.F90 $(SRCROOT)/MOM6/config_src/memory/dynamic_nonsymmetric/MOM_memory.h $(SRCROOT)/MOM6/src/framework/MOM_memory_macros.h MOM_grid.o MOM_verticalGrid.o MOM_forcing_type.o MOM_error_handler.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)/MOM6/config_src/memory/dynamic_nonsymmetric -I$(SRCROOT)/MOM6/src/framework	$(SRCROOT)/MOM6/src/tracer/MOM_tracer_diabatic.F90
MOM_tracer_flow_control.o: $(SRCROOT)/MOM6/src/tracer/MOM_tracer_flow_control.F90 $(SRCROOT)/MOM6/config_src/memory/dynamic_nonsymmetric/MOM_memory.h $(SRCROOT)/MOM6/src/framework/MOM_memory_macros.h $(SRCROOT)/MOM6/src/framework/version_variable.h MOM_coms.o MOM_diag_mediator.o MOM_error_handler.o MOM_file_parser.o MOM_forcing_type.o MOM_get_input.o MOM_grid.o MOM_hor_index.o MOM_CVMix_KPP.o MOM_open_boundary.o MOM_restart.o MOM_sponge.o MOM_ALE_sponge.o MOM_tracer_registry.o MOM_unit_scaling.o MOM_variables.o MOM_verticalGrid.o tracer_example.o DOME_tracer.o ISOMIP_tracer.o RGC_tracer.o ideal_age_example.o dye_example.o MOM_OCMIP2_CFC.o MOM_CFC_cap.o oil_tracer.o advection_test_tracer.o dyed_obc_tracer.o MOM_generic_tracer.o pseudo_salt_tracer.o boundary_impulse_tracer.o nw2_tracers.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)/MOM6/config_src/memory/dynamic_nonsymmetric -I$(SRCROOT)/MOM6/src/framework	$(SRCROOT)/MOM6/src/tracer/MOM_tracer_flow_control.F90
MOM_tracer_hor_diff.o: $(SRCROOT)/MOM6/src/tracer/MOM_tracer_hor_diff.F90 $(SRCROOT)/MOM6/config_src/memory/dynamic_nonsymmetric/MOM_memory.h $(SRCROOT)/MOM6/src/framework/MOM_memory_macros.h $(SRCROOT)/MOM6/src/framework/version_variable.h MOM_cpu_clock.o MOM_diag_mediator.o MOM_domains.o MOM_debugging.o MOM_diabatic_driver.o MOM_EOS.o MOM_error_handler.o MOM_file_parser.o MOM_grid.o MOM_lateral_mixing_coeffs.o MOM_MEKE_types.o MOM_neutral_diffusion.o MOM_lateral_boundary_diffusion.o MOM_tracer_registry.o MOM_unit_scaling.o MOM_variables.o MOM_verticalGrid.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)/MOM6/config_src/memory/dynamic_nonsymmetric -I$(SRCROOT)/MOM6/src/framework	$(SRCROOT)/MOM6/src/tracer/MOM_tracer_hor_diff.F90
MOM_tracer_initialization_from_Z.o: $(SRCROOT)/MOM6/src/initialization/MOM_tracer_initialization_from_Z.F90 $(SRCROOT)/MOM6/config_src/memory/dynamic_nonsymmetric/MOM_memory.h $(SRCROOT)/MOM6/src/framework/MOM_memory_macros.h $(SRCROOT)/MOM6/src/framework/version_variable.h MOM_debugging.o MOM_cpu_clock.o MOM_domains.o MOM_error_handler.o MOM_file_parser.o MOM_grid.o MOM_horizontal_regridding.o MOM_interface_heights.o MOM_remapping.o MOM_unit_scaling.o MOM_verticalGrid.o MOM_ALE.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)/MOM6/config_src/memory/dynamic_nonsymmetric -I$(SRCROOT)/MOM6/src/framework	$(SRCROOT)/MOM6/src/initialization/MOM_tracer_initialization_from_Z.F90
MOM_tracer_registry.o: $(SRCROOT)/MOM6/src/tracer/MOM_tracer_registry.F90 $(SRCROOT)/MOM6/config_src/memory/dynamic_nonsymmetric/MOM_memory.h $(SRCROOT)/MOM6/src/framework/MOM_memory_macros.h $(SRCROOT)/MOM6/src/framework/version_variable.h MOM_coms.o MOM_debugging.o MOM_diag_mediator.o MOM_error_handler.o MOM_file_parser.o MOM_hor_index.o MOM_grid.o MOM_io.o MOM_restart.o MOM_string_functions.o MOM_time_manager.o MOM_unit_scaling.o MOM_verticalGrid.o MOM_tracer_types.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)/MOM6/config_src/memory/dynamic_nonsymmetric -I$(SRCROOT)/MOM6/src/framework	$(SRCROOT)/MOM6/src/tracer/MOM_tracer_registry.F90
MOM_tracer_types.o: $(SRCROOT)/MOM6/src/tracer/MOM_tracer_types.F90 $(SRCROOT)/MOM6/config_src/memory/dynamic_nonsymmetric/MOM_memory.h $(SRCROOT)/MOM6/src/framework/MOM_memory_macros.h
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)/MOM6/config_src/memory/dynamic_nonsymmetric -I$(SRCROOT)/MOM6/src/framework	$(SRCROOT)/MOM6/src/tracer/MOM_tracer_types.F90
MOM_transcribe_grid.o: $(SRCROOT)/MOM6/src/core/MOM_transcribe_grid.F90 MOM_array_transform.o MOM_domains.o MOM_dyn_horgrid.o MOM_error_handler.o MOM_grid.o MOM_unit_scaling.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)/MOM6/src/core/MOM_transcribe_grid.F90
MOM_unique_scales.o: $(SRCROOT)/MOM6/src/framework/MOM_unique_scales.F90 MOM_error_handler.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)/MOM6/src/framework/MOM_unique_scales.F90
MOM_unit_scaling.o: $(SRCROOT)/MOM6/src/framework/MOM_unit_scaling.F90 $(SRCROOT)/MOM6/src/framework/version_variable.h MOM_error_handler.o MOM_file_parser.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)/MOM6/src/framework	$(SRCROOT)/MOM6/src/framework/MOM_unit_scaling.F90
MOM_unit_testing.o: $(SRCROOT)/MOM6/src/framework/MOM_unit_testing.F90 posix.o MOM_coms.o MOM_error_handler.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)/MOM6/src/framework/MOM_unit_testing.F90
MOM_unit_tests.o: $(SRCROOT)/MOM6/src/core/MOM_unit_tests.F90 MOM_error_handler.o MOM_string_functions.o MOM_remapping.o MOM_neutral_diffusion.o MOM_random.o MOM_lateral_boundary_diffusion.o MOM_CFC_cap.o MOM_EOS.o MOM_mixed_layer_restrat.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)/MOM6/src/core/MOM_unit_tests.F90
MOM_variables.o: $(SRCROOT)/MOM6/src/core/MOM_variables.F90 $(SRCROOT)/MOM6/config_src/memory/dynamic_nonsymmetric/MOM_memory.h $(SRCROOT)/MOM6/src/framework/MOM_memory_macros.h MOM_array_transform.o MOM_coupler_types.o MOM_debugging.o MOM_domains.o MOM_EOS.o MOM_error_handler.o MOM_grid.o MOM_unit_scaling.o MOM_verticalGrid.o MOM_tracer_types.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)/MOM6/config_src/memory/dynamic_nonsymmetric -I$(SRCROOT)/MOM6/src/framework	$(SRCROOT)/MOM6/src/core/MOM_variables.F90
MOM_vert_friction.o: $(SRCROOT)/MOM6/src/parameterizations/vertical/MOM_vert_friction.F90 $(SRCROOT)/MOM6/config_src/memory/dynamic_nonsymmetric/MOM_memory.h $(SRCROOT)/MOM6/src/framework/MOM_memory_macros.h $(SRCROOT)/MOM6/src/framework/version_variable.h MOM_domains.o MOM_diag_mediator.o MOM_debugging.o MOM_error_handler.o MOM_file_parser.o MOM_forcing_type.o MOM_get_input.o MOM_grid.o MOM_io.o MOM_open_boundary.o MOM_PointAccel.o MOM_time_manager.o MOM_unit_scaling.o MOM_variables.o MOM_verticalGrid.o MOM_wave_interface.o MOM_lateral_mixing_coeffs.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)/MOM6/config_src/memory/dynamic_nonsymmetric -I$(SRCROOT)/MOM6/src/framework	$(SRCROOT)/MOM6/src/parameterizations/vertical/MOM_vert_friction.F90
MOM_verticalGrid.o: $(SRCROOT)/MOM6/src/core/MOM_verticalGrid.F90 $(SRCROOT)/MOM6/config_src/memory/dynamic_nonsymmetric/MOM_memory.h $(SRCROOT)/MOM6/src/framework/MOM_memory_macros.h $(SRCROOT)/MOM6/src/framework/version_variable.h MOM_error_handler.o MOM_file_parser.o MOM_unit_scaling.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)/MOM6/config_src/memory/dynamic_nonsymmetric -I$(SRCROOT)/MOM6/src/framework	$(SRCROOT)/MOM6/src/core/MOM_verticalGrid.F90
MOM_wave_interface.o: $(SRCROOT)/MOM6/src/user/MOM_wave_interface.F90 $(SRCROOT)/MOM6/config_src/memory/dynamic_nonsymmetric/MOM_memory.h $(SRCROOT)/MOM6/src/framework/MOM_memory_macros.h $(SRCROOT)/MOM6/src/framework/version_variable.h MOM_data_override.o MOM_diag_mediator.o MOM_domains.o MOM_error_handler.o MOM_file_parser.o MOM_forcing_type.o MOM_grid.o MOM_hor_index.o MOM_io.o MOM_safe_alloc.o MOM_time_manager.o MOM_unit_scaling.o MOM_variables.o MOM_verticalGrid.o MOM_restart.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)/MOM6/config_src/memory/dynamic_nonsymmetric -I$(SRCROOT)/MOM6/src/framework	$(SRCROOT)/MOM6/src/user/MOM_wave_interface.F90
MOM_wave_speed.o: $(SRCROOT)/MOM6/src/diagnostics/MOM_wave_speed.F90 $(SRCROOT)/MOM6/config_src/memory/dynamic_nonsymmetric/MOM_memory.h $(SRCROOT)/MOM6/src/framework/MOM_memory_macros.h $(SRCROOT)/MOM6/src/framework/version_variable.h MOM_diag_mediator.o MOM_error_handler.o MOM_file_parser.o MOM_grid.o MOM_remapping.o MOM_unit_scaling.o MOM_variables.o MOM_verticalGrid.o MOM_EOS.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)/MOM6/config_src/memory/dynamic_nonsymmetric -I$(SRCROOT)/MOM6/src/framework	$(SRCROOT)/MOM6/src/diagnostics/MOM_wave_speed.F90
MOM_write_cputime.o: $(SRCROOT)/MOM6/src/framework/MOM_write_cputime.F90 $(SRCROOT)/MOM6/src/framework/version_variable.h MOM_coms.o MOM_error_handler.o MOM_io.o MOM_file_parser.o MOM_time_manager.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)/MOM6/src/framework	$(SRCROOT)/MOM6/src/framework/MOM_write_cputime.F90
Neverworld_initialization.o: $(SRCROOT)/MOM6/src/user/Neverworld_initialization.F90 $(SRCROOT)/MOM6/config_src/memory/dynamic_nonsymmetric/MOM_memory.h $(SRCROOT)/MOM6/src/framework/MOM_memory_macros.h $(SRCROOT)/MOM6/src/framework/version_variable.h MOM_sponge.o MOM_dyn_horgrid.o MOM_error_handler.o MOM_file_parser.o MOM_get_input.o MOM_grid.o MOM_tracer_registry.o MOM_unit_scaling.o MOM_variables.o MOM_verticalGrid.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)/MOM6/config_src/memory/dynamic_nonsymmetric -I$(SRCROOT)/MOM6/src/framework	$(SRCROOT)/MOM6/src/user/Neverworld_initialization.F90
P1M_functions.o: $(SRCROOT)/MOM6/src/ALE/P1M_functions.F90 regrid_edge_values.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)/MOM6/src/ALE/P1M_functions.F90
P3M_functions.o: $(SRCROOT)/MOM6/src/ALE/P3M_functions.F90 regrid_edge_values.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)/MOM6/src/ALE/P3M_functions.F90
PCM_functions.o: $(SRCROOT)/MOM6/src/ALE/PCM_functions.F90
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)/MOM6/src/ALE/PCM_functions.F90
PLM_functions.o: $(SRCROOT)/MOM6/src/ALE/PLM_functions.F90
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)/MOM6/src/ALE/PLM_functions.F90
PPM_functions.o: $(SRCROOT)/MOM6/src/ALE/PPM_functions.F90 regrid_edge_values.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)/MOM6/src/ALE/PPM_functions.F90
PQM_functions.o: $(SRCROOT)/MOM6/src/ALE/PQM_functions.F90 regrid_edge_values.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)/MOM6/src/ALE/PQM_functions.F90
Phillips_initialization.o: $(SRCROOT)/MOM6/src/user/Phillips_initialization.F90 $(SRCROOT)/MOM6/config_src/memory/dynamic_nonsymmetric/MOM_memory.h $(SRCROOT)/MOM6/src/framework/MOM_memory_macros.h $(SRCROOT)/MOM6/src/framework/version_variable.h MOM_error_handler.o MOM_dyn_horgrid.o MOM_file_parser.o MOM_get_input.o MOM_grid.o MOM_sponge.o MOM_tracer_registry.o MOM_unit_scaling.o MOM_variables.o MOM_verticalGrid.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)/MOM6/config_src/memory/dynamic_nonsymmetric -I$(SRCROOT)/MOM6/src/framework	$(SRCROOT)/MOM6/src/user/Phillips_initialization.F90
RGC_initialization.o: $(SRCROOT)/MOM6/src/user/RGC_initialization.F90 $(SRCROOT)/MOM6/config_src/memory/dynamic_nonsymmetric/MOM_memory.h $(SRCROOT)/MOM6/src/framework/MOM_memory_macros.h MOM_ALE_sponge.o MOM_domains.o MOM_dyn_horgrid.o MOM_error_handler.o MOM_file_parser.o MOM_get_input.o MOM_grid.o MOM_io.o MOM_sponge.o MOM_unit_scaling.o MOM_variables.o MOM_verticalGrid.o MOM_EOS.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)/MOM6/config_src/memory/dynamic_nonsymmetric -I$(SRCROOT)/MOM6/src/framework	$(SRCROOT)/MOM6/src/user/RGC_initialization.F90
RGC_tracer.o: $(SRCROOT)/MOM6/src/tracer/RGC_tracer.F90 $(SRCROOT)/MOM6/config_src/memory/dynamic_nonsymmetric/MOM_memory.h $(SRCROOT)/MOM6/src/framework/MOM_memory_macros.h $(SRCROOT)/MOM6/src/framework/version_variable.h MOM_diag_mediator.o MOM_error_handler.o MOM_file_parser.o MOM_forcing_type.o MOM_hor_index.o MOM_grid.o MOM_io.o MOM_restart.o MOM_ALE_sponge.o MOM_sponge.o MOM_time_manager.o MOM_tracer_registry.o MOM_tracer_diabatic.o MOM_unit_scaling.o MOM_variables.o MOM_open_boundary.o MOM_verticalGrid.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)/MOM6/config_src/memory/dynamic_nonsymmetric -I$(SRCROOT)/MOM6/src/framework	$(SRCROOT)/MOM6/src/tracer/RGC_tracer.F90
Rossby_front_2d_initialization.o: $(SRCROOT)/MOM6/src/user/Rossby_front_2d_initialization.F90 $(SRCROOT)/MOM6/config_src/memory/dynamic_nonsymmetric/MOM_memory.h $(SRCROOT)/MOM6/src/framework/MOM_memory_macros.h $(SRCROOT)/MOM6/src/framework/version_variable.h MOM_error_handler.o MOM_file_parser.o MOM_get_input.o MOM_grid.o MOM_unit_scaling.o MOM_variables.o MOM_verticalGrid.o regrid_consts.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)/MOM6/config_src/memory/dynamic_nonsymmetric -I$(SRCROOT)/MOM6/src/framework	$(SRCROOT)/MOM6/src/user/Rossby_front_2d_initialization.F90
SCM_CVMix_tests.o: $(SRCROOT)/MOM6/src/user/SCM_CVMix_tests.F90 $(SRCROOT)/MOM6/config_src/memory/dynamic_nonsymmetric/MOM_memory.h $(SRCROOT)/MOM6/src/framework/MOM_memory_macros.h $(SRCROOT)/MOM6/src/framework/version_variable.h $(SRCROOT)MOM6/src/framework/version_variable.h MOM_domains.o MOM_error_handler.o MOM_file_parser.o MOM_forcing_type.o MOM_grid.o MOM_verticalGrid.o MOM_safe_alloc.o MOM_unit_scaling.o MOM_time_manager.o MOM_variables.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)/MOM6/config_src/memory/dynamic_nonsymmetric -I$(SRCROOT)/MOM6/src/framework -I$(SRCROOT)MOM6/src/framework	$(SRCROOT)/MOM6/src/user/SCM_CVMix_tests.F90
adjustment_initialization.o: $(SRCROOT)/MOM6/src/user/adjustment_initialization.F90 $(SRCROOT)/MOM6/config_src/memory/dynamic_nonsymmetric/MOM_memory.h $(SRCROOT)/MOM6/src/framework/MOM_memory_macros.h $(SRCROOT)/MOM6/src/framework/version_variable.h MOM_error_handler.o MOM_file_parser.o MOM_get_input.o MOM_grid.o MOM_unit_scaling.o MOM_variables.o MOM_verticalGrid.o regrid_consts.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)/MOM6/config_src/memory/dynamic_nonsymmetric -I$(SRCROOT)/MOM6/src/framework	$(SRCROOT)/MOM6/src/user/adjustment_initialization.F90
advection_test_tracer.o: $(SRCROOT)/MOM6/src/tracer/advection_test_tracer.F90 $(SRCROOT)/MOM6/config_src/memory/dynamic_nonsymmetric/MOM_memory.h $(SRCROOT)/MOM6/src/framework/MOM_memory_macros.h $(SRCROOT)/MOM6/src/framework/version_variable.h MOM_coms.o MOM_coupler_types.o MOM_diag_mediator.o MOM_error_handler.o MOM_file_parser.o MOM_forcing_type.o MOM_grid.o MOM_hor_index.o MOM_io.o MOM_open_boundary.o MOM_restart.o MOM_spatial_means.o MOM_sponge.o MOM_time_manager.o MOM_tracer_registry.o MOM_tracer_diabatic.o MOM_unit_scaling.o MOM_variables.o MOM_verticalGrid.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)/MOM6/config_src/memory/dynamic_nonsymmetric -I$(SRCROOT)/MOM6/src/framework	$(SRCROOT)/MOM6/src/tracer/advection_test_tracer.F90
baroclinic_zone_initialization.o: $(SRCROOT)/MOM6/src/user/baroclinic_zone_initialization.F90 $(SRCROOT)/MOM6/config_src/memory/dynamic_nonsymmetric/MOM_memory.h $(SRCROOT)/MOM6/src/framework/MOM_memory_macros.h $(SRCROOT)/MOM6/src/framework/version_variable.h MOM_file_parser.o MOM_grid.o MOM_unit_scaling.o MOM_verticalGrid.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)/MOM6/config_src/memory/dynamic_nonsymmetric -I$(SRCROOT)/MOM6/src/framework	$(SRCROOT)/MOM6/src/user/baroclinic_zone_initialization.F90
basin_builder.o: $(SRCROOT)/MOM6/src/user/basin_builder.F90 $(SRCROOT)/MOM6/config_src/memory/dynamic_nonsymmetric/MOM_memory.h $(SRCROOT)/MOM6/src/framework/MOM_memory_macros.h $(SRCROOT)/MOM6/src/framework/version_variable.h MOM_dyn_horgrid.o MOM_error_handler.o MOM_file_parser.o MOM_get_input.o MOM_grid.o MOM_string_functions.o MOM_unit_scaling.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)/MOM6/config_src/memory/dynamic_nonsymmetric -I$(SRCROOT)/MOM6/src/framework	$(SRCROOT)/MOM6/src/user/basin_builder.F90
benchmark_initialization.o: $(SRCROOT)/MOM6/src/user/benchmark_initialization.F90 $(SRCROOT)/MOM6/config_src/memory/dynamic_nonsymmetric/MOM_memory.h $(SRCROOT)/MOM6/src/framework/MOM_memory_macros.h $(SRCROOT)/MOM6/src/framework/version_variable.h $(SRCROOT)MOM6/src/framework/version_variable.h MOM_sponge.o MOM_dyn_horgrid.o MOM_error_handler.o MOM_file_parser.o MOM_get_input.o MOM_grid.o MOM_tracer_registry.o MOM_unit_scaling.o MOM_variables.o MOM_verticalGrid.o MOM_EOS.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)/MOM6/config_src/memory/dynamic_nonsymmetric -I$(SRCROOT)/MOM6/src/framework -I$(SRCROOT)MOM6/src/framework	$(SRCROOT)/MOM6/src/user/benchmark_initialization.F90
boundary_impulse_tracer.o: $(SRCROOT)/MOM6/src/tracer/boundary_impulse_tracer.F90 $(SRCROOT)/MOM6/config_src/memory/dynamic_nonsymmetric/MOM_memory.h $(SRCROOT)/MOM6/src/framework/MOM_memory_macros.h $(SRCROOT)/MOM6/src/framework/version_variable.h MOM_coms.o MOM_coupler_types.o MOM_diag_mediator.o MOM_error_handler.o MOM_file_parser.o MOM_forcing_type.o MOM_grid.o MOM_hor_index.o MOM_io.o MOM_open_boundary.o MOM_restart.o MOM_spatial_means.o MOM_sponge.o MOM_time_manager.o MOM_tracer_registry.o MOM_tracer_diabatic.o MOM_tracer_Z_init.o MOM_unit_scaling.o MOM_variables.o MOM_verticalGrid.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)/MOM6/config_src/memory/dynamic_nonsymmetric -I$(SRCROOT)/MOM6/src/framework	$(SRCROOT)/MOM6/src/tracer/boundary_impulse_tracer.F90
circle_obcs_initialization.o: $(SRCROOT)/MOM6/src/user/circle_obcs_initialization.F90 $(SRCROOT)/MOM6/config_src/memory/dynamic_nonsymmetric/MOM_memory.h $(SRCROOT)/MOM6/src/framework/MOM_memory_macros.h $(SRCROOT)/MOM6/src/framework/version_variable.h MOM_sponge.o MOM_error_handler.o MOM_file_parser.o MOM_get_input.o MOM_grid.o MOM_tracer_registry.o MOM_unit_scaling.o MOM_variables.o MOM_verticalGrid.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)/MOM6/config_src/memory/dynamic_nonsymmetric -I$(SRCROOT)/MOM6/src/framework	$(SRCROOT)/MOM6/src/user/circle_obcs_initialization.F90
coord_adapt.o: $(SRCROOT)/MOM6/src/ALE/coord_adapt.F90 $(SRCROOT)/MOM6/config_src/memory/dynamic_nonsymmetric/MOM_memory.h $(SRCROOT)/MOM6/src/framework/MOM_memory_macros.h MOM_EOS.o MOM_error_handler.o MOM_unit_scaling.o MOM_variables.o MOM_verticalGrid.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)/MOM6/config_src/memory/dynamic_nonsymmetric -I$(SRCROOT)/MOM6/src/framework	$(SRCROOT)/MOM6/src/ALE/coord_adapt.F90
coord_hycom.o: $(SRCROOT)/MOM6/src/ALE/coord_hycom.F90 MOM_error_handler.o MOM_remapping.o MOM_EOS.o regrid_interp.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)/MOM6/src/ALE/coord_hycom.F90
coord_rho.o: $(SRCROOT)/MOM6/src/ALE/coord_rho.F90 MOM_error_handler.o MOM_remapping.o MOM_EOS.o regrid_interp.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)/MOM6/src/ALE/coord_rho.F90
coord_sigma.o: $(SRCROOT)/MOM6/src/ALE/coord_sigma.F90 MOM_error_handler.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)/MOM6/src/ALE/coord_sigma.F90
coord_zlike.o: $(SRCROOT)/MOM6/src/ALE/coord_zlike.F90 MOM_error_handler.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)/MOM6/src/ALE/coord_zlike.F90
cvmix_background.o: $(SRCROOT)/MOM6/src/parameterizations/CVmix/cvmix_background.F90 cvmix_kinds_and_types.o cvmix_put_get.o cvmix_utils.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)/MOM6/src/parameterizations/CVmix/cvmix_background.F90
cvmix_convection.o: $(SRCROOT)/MOM6/src/parameterizations/CVmix/cvmix_convection.F90 cvmix_kinds_and_types.o cvmix_utils.o cvmix_put_get.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)/MOM6/src/parameterizations/CVmix/cvmix_convection.F90
cvmix_ddiff.o: $(SRCROOT)/MOM6/src/parameterizations/CVmix/cvmix_ddiff.F90 cvmix_kinds_and_types.o cvmix_put_get.o cvmix_utils.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)/MOM6/src/parameterizations/CVmix/cvmix_ddiff.F90
cvmix_kinds_and_types.o: $(SRCROOT)/MOM6/src/parameterizations/CVmix/cvmix_kinds_and_types.F90
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)/MOM6/src/parameterizations/CVmix/cvmix_kinds_and_types.F90
cvmix_kpp.o: $(SRCROOT)/MOM6/src/parameterizations/CVmix/cvmix_kpp.F90 cvmix_kinds_and_types.o cvmix_math.o cvmix_put_get.o cvmix_utils.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)/MOM6/src/parameterizations/CVmix/cvmix_kpp.F90
cvmix_math.o: $(SRCROOT)/MOM6/src/parameterizations/CVmix/cvmix_math.F90 cvmix_kinds_and_types.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)/MOM6/src/parameterizations/CVmix/cvmix_math.F90
cvmix_put_get.o: $(SRCROOT)/MOM6/src/parameterizations/CVmix/cvmix_put_get.F90 cvmix_kinds_and_types.o cvmix_utils.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)/MOM6/src/parameterizations/CVmix/cvmix_put_get.F90
cvmix_shear.o: $(SRCROOT)/MOM6/src/parameterizations/CVmix/cvmix_shear.F90 cvmix_kinds_and_types.o cvmix_put_get.o cvmix_utils.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)/MOM6/src/parameterizations/CVmix/cvmix_shear.F90
cvmix_tidal.o: $(SRCROOT)/MOM6/src/parameterizations/CVmix/cvmix_tidal.F90 cvmix_kinds_and_types.o cvmix_utils.o cvmix_put_get.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)/MOM6/src/parameterizations/CVmix/cvmix_tidal.F90
cvmix_utils.o: $(SRCROOT)/MOM6/src/parameterizations/CVmix/cvmix_utils.F90 cvmix_kinds_and_types.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)/MOM6/src/parameterizations/CVmix/cvmix_utils.F90
database_client_interface.o: $(SRCROOT)/MOM6/config_src/external/database_comms/database_client_interface.F90
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)/MOM6/config_src/external/database_comms/database_client_interface.F90
dense_water_initialization.o: $(SRCROOT)/MOM6/src/user/dense_water_initialization.F90 $(SRCROOT)/MOM6/config_src/memory/dynamic_nonsymmetric/MOM_memory.h $(SRCROOT)/MOM6/src/framework/MOM_memory_macros.h MOM_ALE_sponge.o MOM_dyn_horgrid.o MOM_EOS.o MOM_error_handler.o MOM_file_parser.o MOM_interface_heights.o MOM_grid.o MOM_sponge.o MOM_unit_scaling.o MOM_variables.o MOM_verticalGrid.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)/MOM6/config_src/memory/dynamic_nonsymmetric -I$(SRCROOT)/MOM6/src/framework	$(SRCROOT)/MOM6/src/user/dense_water_initialization.F90
dumbbell_initialization.o: $(SRCROOT)/MOM6/src/user/dumbbell_initialization.F90 $(SRCROOT)/MOM6/config_src/memory/dynamic_nonsymmetric/MOM_memory.h $(SRCROOT)/MOM6/src/framework/MOM_memory_macros.h $(SRCROOT)/MOM6/src/framework/version_variable.h MOM_domains.o MOM_dyn_horgrid.o MOM_error_handler.o MOM_file_parser.o MOM_get_input.o MOM_grid.o MOM_interface_heights.o MOM_sponge.o MOM_tracer_registry.o MOM_unit_scaling.o MOM_variables.o MOM_verticalGrid.o regrid_consts.o MOM_ALE_sponge.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)/MOM6/config_src/memory/dynamic_nonsymmetric -I$(SRCROOT)/MOM6/src/framework	$(SRCROOT)/MOM6/src/user/dumbbell_initialization.F90
dumbbell_surface_forcing.o: $(SRCROOT)/MOM6/src/user/dumbbell_surface_forcing.F90 $(SRCROOT)/MOM6/src/framework/version_variable.h MOM_diag_mediator.o MOM_domains.o MOM_error_handler.o MOM_file_parser.o MOM_forcing_type.o MOM_grid.o MOM_safe_alloc.o MOM_time_manager.o MOM_tracer_flow_control.o MOM_unit_scaling.o MOM_variables.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)/MOM6/src/framework	$(SRCROOT)/MOM6/src/user/dumbbell_surface_forcing.F90
dye_example.o: $(SRCROOT)/MOM6/src/tracer/dye_example.F90 $(SRCROOT)/MOM6/config_src/memory/dynamic_nonsymmetric/MOM_memory.h $(SRCROOT)/MOM6/src/framework/MOM_memory_macros.h $(SRCROOT)/MOM6/src/framework/version_variable.h MOM_coms.o MOM_coupler_types.o MOM_diag_mediator.o MOM_error_handler.o MOM_file_parser.o MOM_forcing_type.o MOM_grid.o MOM_hor_index.o MOM_io.o MOM_open_boundary.o MOM_restart.o MOM_spatial_means.o MOM_sponge.o MOM_time_manager.o MOM_tracer_registry.o MOM_tracer_diabatic.o MOM_tracer_Z_init.o MOM_unit_scaling.o MOM_variables.o MOM_verticalGrid.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)/MOM6/config_src/memory/dynamic_nonsymmetric -I$(SRCROOT)/MOM6/src/framework	$(SRCROOT)/MOM6/src/tracer/dye_example.F90
dyed_channel_initialization.o: $(SRCROOT)/MOM6/src/user/dyed_channel_initialization.F90 $(SRCROOT)/MOM6/config_src/memory/dynamic_nonsymmetric/MOM_memory.h $(SRCROOT)/MOM6/src/framework/MOM_memory_macros.h MOM_dyn_horgrid.o MOM_error_handler.o MOM_file_parser.o MOM_get_input.o MOM_grid.o MOM_open_boundary.o MOM_time_manager.o MOM_tracer_registry.o MOM_unit_scaling.o MOM_variables.o MOM_verticalGrid.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)/MOM6/config_src/memory/dynamic_nonsymmetric -I$(SRCROOT)/MOM6/src/framework	$(SRCROOT)/MOM6/src/user/dyed_channel_initialization.F90
dyed_obc_tracer.o: $(SRCROOT)/MOM6/src/tracer/dyed_obc_tracer.F90 $(SRCROOT)/MOM6/config_src/memory/dynamic_nonsymmetric/MOM_memory.h $(SRCROOT)/MOM6/src/framework/MOM_memory_macros.h $(SRCROOT)/MOM6/src/framework/version_variable.h MOM_coupler_types.o MOM_diag_mediator.o MOM_error_handler.o MOM_file_parser.o MOM_forcing_type.o MOM_hor_index.o MOM_grid.o MOM_io.o MOM_open_boundary.o MOM_restart.o MOM_time_manager.o MOM_tracer_registry.o MOM_tracer_diabatic.o MOM_unit_scaling.o MOM_variables.o MOM_verticalGrid.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)/MOM6/config_src/memory/dynamic_nonsymmetric -I$(SRCROOT)/MOM6/src/framework	$(SRCROOT)/MOM6/src/tracer/dyed_obc_tracer.F90
dyed_obcs_initialization.o: $(SRCROOT)/MOM6/src/user/dyed_obcs_initialization.F90 $(SRCROOT)/MOM6/config_src/memory/dynamic_nonsymmetric/MOM_memory.h $(SRCROOT)/MOM6/src/framework/MOM_memory_macros.h MOM_dyn_horgrid.o MOM_error_handler.o MOM_file_parser.o MOM_get_input.o MOM_grid.o MOM_open_boundary.o MOM_tracer_registry.o MOM_variables.o MOM_verticalGrid.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)/MOM6/config_src/memory/dynamic_nonsymmetric -I$(SRCROOT)/MOM6/src/framework	$(SRCROOT)/MOM6/src/user/dyed_obcs_initialization.F90
external_gwave_initialization.o: $(SRCROOT)/MOM6/src/user/external_gwave_initialization.F90 $(SRCROOT)/MOM6/config_src/memory/dynamic_nonsymmetric/MOM_memory.h $(SRCROOT)/MOM6/src/framework/MOM_memory_macros.h $(SRCROOT)/MOM6/src/framework/version_variable.h MOM_error_handler.o MOM_file_parser.o MOM_get_input.o MOM_grid.o MOM_tracer_registry.o MOM_unit_scaling.o MOM_variables.o MOM_verticalGrid.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)/MOM6/config_src/memory/dynamic_nonsymmetric -I$(SRCROOT)/MOM6/src/framework	$(SRCROOT)/MOM6/src/user/external_gwave_initialization.F90
generic_tracer.o: $(SRCROOT)/MOM6/config_src/external/GFDL_ocean_BGC/generic_tracer.F90 generic_tracer_utils.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)/MOM6/config_src/external/GFDL_ocean_BGC/generic_tracer.F90
generic_tracer_utils.o: $(SRCROOT)/MOM6/config_src/external/GFDL_ocean_BGC/generic_tracer_utils.F90 MOM_diag_mediator.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)/MOM6/config_src/external/GFDL_ocean_BGC/generic_tracer_utils.F90
get_stochy_pattern.o: $(SRCROOT)/MOM6/config_src/external/stochastic_physics/get_stochy_pattern.F90
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)/MOM6/config_src/external/stochastic_physics/get_stochy_pattern.F90
gsw_chem_potential_water_t_exact.o: $(SRCROOT)/MOM6/src/equation_of_state/TEOS10/gsw_chem_potential_water_t_exact.f90 gsw_mod_teos10_constants.o gsw_mod_kinds.o
	$(FC) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)/MOM6/src/equation_of_state/TEOS10/gsw_chem_potential_water_t_exact.f90
gsw_ct_freezing_exact.o: $(SRCROOT)/MOM6/src/equation_of_state/TEOS10/gsw_ct_freezing_exact.f90 gsw_mod_toolbox.o gsw_mod_kinds.o
	$(FC) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)/MOM6/src/equation_of_state/TEOS10/gsw_ct_freezing_exact.f90
gsw_ct_freezing_poly.o: $(SRCROOT)/MOM6/src/equation_of_state/TEOS10/gsw_ct_freezing_poly.f90 gsw_mod_teos10_constants.o gsw_mod_freezing_poly_coefficients.o gsw_mod_kinds.o
	$(FC) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)/MOM6/src/equation_of_state/TEOS10/gsw_ct_freezing_poly.f90
gsw_ct_from_pt.o: $(SRCROOT)/MOM6/src/equation_of_state/TEOS10/gsw_ct_from_pt.f90 gsw_mod_teos10_constants.o gsw_mod_kinds.o
	$(FC) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)/MOM6/src/equation_of_state/TEOS10/gsw_ct_from_pt.f90
gsw_ct_from_t.o: $(SRCROOT)/MOM6/src/equation_of_state/TEOS10/gsw_ct_from_t.f90 gsw_mod_toolbox.o gsw_mod_kinds.o
	$(FC) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)/MOM6/src/equation_of_state/TEOS10/gsw_ct_from_t.f90
gsw_entropy_part.o: $(SRCROOT)/MOM6/src/equation_of_state/TEOS10/gsw_entropy_part.f90 gsw_mod_teos10_constants.o gsw_mod_kinds.o
	$(FC) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)/MOM6/src/equation_of_state/TEOS10/gsw_entropy_part.f90
gsw_entropy_part_zerop.o: $(SRCROOT)/MOM6/src/equation_of_state/TEOS10/gsw_entropy_part_zerop.f90 gsw_mod_teos10_constants.o gsw_mod_kinds.o
	$(FC) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)/MOM6/src/equation_of_state/TEOS10/gsw_entropy_part_zerop.f90
gsw_gibbs.o: $(SRCROOT)/MOM6/src/equation_of_state/TEOS10/gsw_gibbs.f90 gsw_mod_teos10_constants.o gsw_mod_kinds.o
	$(FC) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)/MOM6/src/equation_of_state/TEOS10/gsw_gibbs.f90
gsw_gibbs_ice.o: $(SRCROOT)/MOM6/src/equation_of_state/TEOS10/gsw_gibbs_ice.f90 gsw_mod_teos10_constants.o gsw_mod_gibbs_ice_coefficients.o gsw_mod_kinds.o
	$(FC) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)/MOM6/src/equation_of_state/TEOS10/gsw_gibbs_ice.f90
gsw_gibbs_pt0_pt0.o: $(SRCROOT)/MOM6/src/equation_of_state/TEOS10/gsw_gibbs_pt0_pt0.f90 gsw_mod_teos10_constants.o gsw_mod_kinds.o
	$(FC) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)/MOM6/src/equation_of_state/TEOS10/gsw_gibbs_pt0_pt0.f90
gsw_mod_error_functions.o: $(SRCROOT)/MOM6/src/equation_of_state/TEOS10/gsw_mod_error_functions.f90 gsw_mod_kinds.o
	$(FC) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)/MOM6/src/equation_of_state/TEOS10/gsw_mod_error_functions.f90
gsw_mod_freezing_poly_coefficients.o: $(SRCROOT)/MOM6/src/equation_of_state/TEOS10/gsw_mod_freezing_poly_coefficients.f90 gsw_mod_kinds.o
	$(FC) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)/MOM6/src/equation_of_state/TEOS10/gsw_mod_freezing_poly_coefficients.f90
gsw_mod_gibbs_ice_coefficients.o: $(SRCROOT)/MOM6/src/equation_of_state/TEOS10/gsw_mod_gibbs_ice_coefficients.f90 gsw_mod_kinds.o
	$(FC) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)/MOM6/src/equation_of_state/TEOS10/gsw_mod_gibbs_ice_coefficients.f90
gsw_mod_kinds.o: $(SRCROOT)/MOM6/src/equation_of_state/TEOS10/gsw_mod_kinds.f90
	$(FC) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)/MOM6/src/equation_of_state/TEOS10/gsw_mod_kinds.f90
gsw_mod_specvol_coefficients.o: $(SRCROOT)/MOM6/src/equation_of_state/TEOS10/gsw_mod_specvol_coefficients.f90 gsw_mod_kinds.o
	$(FC) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)/MOM6/src/equation_of_state/TEOS10/gsw_mod_specvol_coefficients.f90
gsw_mod_teos10_constants.o: $(SRCROOT)/MOM6/src/equation_of_state/TEOS10/gsw_mod_teos10_constants.f90 gsw_mod_kinds.o
	$(FC) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)/MOM6/src/equation_of_state/TEOS10/gsw_mod_teos10_constants.f90
gsw_mod_toolbox.o: $(SRCROOT)/MOM6/src/equation_of_state/TEOS10/gsw_mod_toolbox.f90 gsw_mod_kinds.o
	$(FC) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)/MOM6/src/equation_of_state/TEOS10/gsw_mod_toolbox.f90
gsw_pt0_from_t.o: $(SRCROOT)/MOM6/src/equation_of_state/TEOS10/gsw_pt0_from_t.f90 gsw_mod_toolbox.o gsw_mod_teos10_constants.o gsw_mod_kinds.o
	$(FC) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)/MOM6/src/equation_of_state/TEOS10/gsw_pt0_from_t.f90
gsw_pt_from_ct.o: $(SRCROOT)/MOM6/src/equation_of_state/TEOS10/gsw_pt_from_ct.f90 gsw_mod_toolbox.o gsw_mod_teos10_constants.o gsw_mod_kinds.o
	$(FC) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)/MOM6/src/equation_of_state/TEOS10/gsw_pt_from_ct.f90
gsw_pt_from_t.o: $(SRCROOT)/MOM6/src/equation_of_state/TEOS10/gsw_pt_from_t.f90 gsw_mod_toolbox.o gsw_mod_teos10_constants.o gsw_mod_kinds.o
	$(FC) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)/MOM6/src/equation_of_state/TEOS10/gsw_pt_from_t.f90
gsw_rho.o: $(SRCROOT)/MOM6/src/equation_of_state/TEOS10/gsw_rho.f90 gsw_mod_toolbox.o gsw_mod_kinds.o
	$(FC) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)/MOM6/src/equation_of_state/TEOS10/gsw_rho.f90
gsw_rho_first_derivatives.o: $(SRCROOT)/MOM6/src/equation_of_state/TEOS10/gsw_rho_first_derivatives.f90 gsw_mod_teos10_constants.o gsw_mod_specvol_coefficients.o gsw_mod_kinds.o
	$(FC) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)/MOM6/src/equation_of_state/TEOS10/gsw_rho_first_derivatives.f90
gsw_rho_second_derivatives.o: $(SRCROOT)/MOM6/src/equation_of_state/TEOS10/gsw_rho_second_derivatives.f90 gsw_mod_toolbox.o gsw_mod_kinds.o
	$(FC) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)/MOM6/src/equation_of_state/TEOS10/gsw_rho_second_derivatives.f90
gsw_sp_from_sr.o: $(SRCROOT)/MOM6/src/equation_of_state/TEOS10/gsw_sp_from_sr.f90 gsw_mod_teos10_constants.o gsw_mod_kinds.o
	$(FC) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)/MOM6/src/equation_of_state/TEOS10/gsw_sp_from_sr.f90
gsw_specvol.o: $(SRCROOT)/MOM6/src/equation_of_state/TEOS10/gsw_specvol.f90 gsw_mod_teos10_constants.o gsw_mod_specvol_coefficients.o gsw_mod_kinds.o
	$(FC) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)/MOM6/src/equation_of_state/TEOS10/gsw_specvol.f90
gsw_specvol_first_derivatives.o: $(SRCROOT)/MOM6/src/equation_of_state/TEOS10/gsw_specvol_first_derivatives.f90 gsw_mod_teos10_constants.o gsw_mod_specvol_coefficients.o gsw_mod_kinds.o
	$(FC) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)/MOM6/src/equation_of_state/TEOS10/gsw_specvol_first_derivatives.f90
gsw_specvol_second_derivatives.o: $(SRCROOT)/MOM6/src/equation_of_state/TEOS10/gsw_specvol_second_derivatives.f90 gsw_mod_teos10_constants.o gsw_mod_specvol_coefficients.o gsw_mod_kinds.o
	$(FC) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)/MOM6/src/equation_of_state/TEOS10/gsw_specvol_second_derivatives.f90
gsw_sr_from_sp.o: $(SRCROOT)/MOM6/src/equation_of_state/TEOS10/gsw_sr_from_sp.f90 gsw_mod_teos10_constants.o gsw_mod_kinds.o
	$(FC) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)/MOM6/src/equation_of_state/TEOS10/gsw_sr_from_sp.f90
gsw_t_deriv_chem_potential_water_t_exact.o: $(SRCROOT)/MOM6/src/equation_of_state/TEOS10/gsw_t_deriv_chem_potential_water_t_exact.f90 gsw_mod_teos10_constants.o gsw_mod_kinds.o
	$(FC) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)/MOM6/src/equation_of_state/TEOS10/gsw_t_deriv_chem_potential_water_t_exact.f90
gsw_t_freezing_exact.o: $(SRCROOT)/MOM6/src/equation_of_state/TEOS10/gsw_t_freezing_exact.f90 gsw_mod_teos10_constants.o gsw_mod_toolbox.o gsw_mod_kinds.o
	$(FC) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)/MOM6/src/equation_of_state/TEOS10/gsw_t_freezing_exact.f90
gsw_t_freezing_poly.o: $(SRCROOT)/MOM6/src/equation_of_state/TEOS10/gsw_t_freezing_poly.f90 gsw_mod_teos10_constants.o gsw_mod_freezing_poly_coefficients.o gsw_mod_toolbox.o gsw_mod_kinds.o
	$(FC) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)/MOM6/src/equation_of_state/TEOS10/gsw_t_freezing_poly.f90
gsw_t_from_ct.o: $(SRCROOT)/MOM6/src/equation_of_state/TEOS10/gsw_t_from_ct.f90 gsw_mod_toolbox.o gsw_mod_kinds.o
	$(FC) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)/MOM6/src/equation_of_state/TEOS10/gsw_t_from_ct.f90
ideal_age_example.o: $(SRCROOT)/MOM6/src/tracer/ideal_age_example.F90 $(SRCROOT)/MOM6/config_src/memory/dynamic_nonsymmetric/MOM_memory.h $(SRCROOT)/MOM6/src/framework/MOM_memory_macros.h $(SRCROOT)/MOM6/src/framework/version_variable.h MOM_coms.o MOM_coupler_types.o MOM_diag_mediator.o MOM_error_handler.o MOM_file_parser.o MOM_forcing_type.o MOM_grid.o MOM_hor_index.o MOM_io.o MOM_open_boundary.o MOM_restart.o MOM_spatial_means.o MOM_sponge.o MOM_time_manager.o MOM_tracer_registry.o MOM_tracer_diabatic.o MOM_tracer_Z_init.o MOM_unit_scaling.o MOM_variables.o MOM_verticalGrid.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)/MOM6/config_src/memory/dynamic_nonsymmetric -I$(SRCROOT)/MOM6/src/framework	$(SRCROOT)/MOM6/src/tracer/ideal_age_example.F90
kdtree.o: $(SRCROOT)/MOM6/config_src/external/ODA_hooks/kdtree.f90
	$(FC) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)/MOM6/config_src/external/ODA_hooks/kdtree.f90
lock_exchange_initialization.o: $(SRCROOT)/MOM6/src/user/lock_exchange_initialization.F90 $(SRCROOT)/MOM6/config_src/memory/dynamic_nonsymmetric/MOM_memory.h $(SRCROOT)/MOM6/src/framework/MOM_memory_macros.h $(SRCROOT)/MOM6/src/framework/version_variable.h MOM_error_handler.o MOM_file_parser.o MOM_get_input.o MOM_grid.o MOM_tracer_registry.o MOM_unit_scaling.o MOM_variables.o MOM_verticalGrid.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)/MOM6/config_src/memory/dynamic_nonsymmetric -I$(SRCROOT)/MOM6/src/framework	$(SRCROOT)/MOM6/src/user/lock_exchange_initialization.F90
nw2_tracers.o: $(SRCROOT)/MOM6/src/tracer/nw2_tracers.F90 $(SRCROOT)/MOM6/config_src/memory/dynamic_nonsymmetric/MOM_memory.h $(SRCROOT)/MOM6/src/framework/MOM_memory_macros.h $(SRCROOT)/MOM6/src/framework/version_variable.h MOM_diag_mediator.o MOM_error_handler.o MOM_file_parser.o MOM_forcing_type.o MOM_grid.o MOM_hor_index.o MOM_io.o MOM_restart.o MOM_time_manager.o MOM_tracer_registry.o MOM_tracer_diabatic.o MOM_unit_scaling.o MOM_variables.o MOM_verticalGrid.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)/MOM6/config_src/memory/dynamic_nonsymmetric -I$(SRCROOT)/MOM6/src/framework	$(SRCROOT)/MOM6/src/tracer/nw2_tracers.F90
ocean_da_core.o: $(SRCROOT)/MOM6/config_src/external/ODA_hooks/ocean_da_core.F90 MOM_domains.o MOM_time_manager.o ocean_da_types.o kdtree.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)/MOM6/config_src/external/ODA_hooks/ocean_da_core.F90
ocean_da_types.o: $(SRCROOT)/MOM6/config_src/external/ODA_hooks/ocean_da_types.F90 MOM_time_manager.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)/MOM6/config_src/external/ODA_hooks/ocean_da_types.F90
ocean_model_MOM.o: $(SRCROOT)/MOM6/config_src/drivers/FMS_cap/ocean_model_MOM.F90 $(SRCROOT)/MOM6/config_src/memory/dynamic_nonsymmetric/MOM_memory.h $(SRCROOT)/MOM6/src/framework/MOM_memory_macros.h $(SRCROOT)/MOM6/src/framework/version_variable.h MOM.o MOM_coms.o MOM_constants.o MOM_coupler_types.o MOM_diag_mediator.o MOM_domains.o MOM_error_handler.o MOM_EOS.o MOM_file_parser.o MOM_forcing_type.o MOM_get_input.o MOM_grid.o MOM_io.o MOM_marine_ice.o MOM_restart.o MOM_string_functions.o MOM_surface_forcing_gfdl.o MOM_time_manager.o MOM_tracer_flow_control.o MOM_unit_scaling.o MOM_variables.o MOM_verticalGrid.o MOM_ice_shelf.o MOM_wave_interface.o MOM_generic_tracer.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)/MOM6/config_src/memory/dynamic_nonsymmetric -I$(SRCROOT)/MOM6/src/framework	$(SRCROOT)/MOM6/config_src/drivers/FMS_cap/ocean_model_MOM.F90
oil_tracer.o: $(SRCROOT)/MOM6/src/tracer/oil_tracer.F90 $(SRCROOT)/MOM6/config_src/memory/dynamic_nonsymmetric/MOM_memory.h $(SRCROOT)/MOM6/src/framework/MOM_memory_macros.h $(SRCROOT)/MOM6/src/framework/version_variable.h MOM_coms.o MOM_coupler_types.o MOM_diag_mediator.o MOM_error_handler.o MOM_file_parser.o MOM_forcing_type.o MOM_grid.o MOM_hor_index.o MOM_io.o MOM_open_boundary.o MOM_restart.o MOM_spatial_means.o MOM_sponge.o MOM_time_manager.o MOM_tracer_registry.o MOM_tracer_diabatic.o MOM_tracer_Z_init.o MOM_unit_scaling.o MOM_variables.o MOM_verticalGrid.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)/MOM6/config_src/memory/dynamic_nonsymmetric -I$(SRCROOT)/MOM6/src/framework	$(SRCROOT)/MOM6/src/tracer/oil_tracer.F90
polynomial_functions.o: $(SRCROOT)/MOM6/src/ALE/polynomial_functions.F90
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)/MOM6/src/ALE/polynomial_functions.F90
posix.o: $(SRCROOT)/MOM6/src/framework/posix.F90 $(SRCROOT)/MOM6/src/framework/posix.h
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)/MOM6/src/framework	$(SRCROOT)/MOM6/src/framework/posix.F90
pseudo_salt_tracer.o: $(SRCROOT)/MOM6/src/tracer/pseudo_salt_tracer.F90 $(SRCROOT)/MOM6/config_src/memory/dynamic_nonsymmetric/MOM_memory.h $(SRCROOT)/MOM6/src/framework/MOM_memory_macros.h $(SRCROOT)/MOM6/src/framework/version_variable.h MOM_coms.o MOM_debugging.o MOM_diag_mediator.o MOM_error_handler.o MOM_file_parser.o MOM_forcing_type.o MOM_grid.o MOM_CVMix_KPP.o MOM_hor_index.o MOM_io.o MOM_open_boundary.o MOM_restart.o MOM_spatial_means.o MOM_sponge.o MOM_time_manager.o MOM_tracer_registry.o MOM_tracer_diabatic.o MOM_tracer_Z_init.o MOM_unit_scaling.o MOM_variables.o MOM_verticalGrid.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)/MOM6/config_src/memory/dynamic_nonsymmetric -I$(SRCROOT)/MOM6/src/framework	$(SRCROOT)/MOM6/src/tracer/pseudo_salt_tracer.F90
regrid_consts.o: $(SRCROOT)/MOM6/src/ALE/regrid_consts.F90 MOM_error_handler.o MOM_string_functions.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)/MOM6/src/ALE/regrid_consts.F90
regrid_edge_values.o: $(SRCROOT)/MOM6/src/ALE/regrid_edge_values.F90 MOM_error_handler.o regrid_solvers.o polynomial_functions.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)/MOM6/src/ALE/regrid_edge_values.F90
regrid_interp.o: $(SRCROOT)/MOM6/src/ALE/regrid_interp.F90 MOM_error_handler.o MOM_string_functions.o regrid_edge_values.o PLM_functions.o PPM_functions.o PQM_functions.o P1M_functions.o P3M_functions.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)/MOM6/src/ALE/regrid_interp.F90
regrid_solvers.o: $(SRCROOT)/MOM6/src/ALE/regrid_solvers.F90 MOM_error_handler.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)/MOM6/src/ALE/regrid_solvers.F90
remapping_attic.o: $(SRCROOT)/MOM6/src/ALE/remapping_attic.F90 MOM_error_handler.o MOM_io.o PPM_functions.o regrid_edge_values.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)/MOM6/src/ALE/remapping_attic.F90
seamount_initialization.o: $(SRCROOT)/MOM6/src/user/seamount_initialization.F90 $(SRCROOT)/MOM6/config_src/memory/dynamic_nonsymmetric/MOM_memory.h $(SRCROOT)/MOM6/src/framework/MOM_memory_macros.h MOM_domains.o MOM_dyn_horgrid.o MOM_error_handler.o MOM_file_parser.o MOM_get_input.o MOM_grid.o MOM_sponge.o MOM_tracer_registry.o MOM_unit_scaling.o MOM_variables.o MOM_verticalGrid.o regrid_consts.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)/MOM6/config_src/memory/dynamic_nonsymmetric -I$(SRCROOT)/MOM6/src/framework	$(SRCROOT)/MOM6/src/user/seamount_initialization.F90
shelfwave_initialization.o: $(SRCROOT)/MOM6/src/user/shelfwave_initialization.F90 $(SRCROOT)/MOM6/config_src/memory/dynamic_nonsymmetric/MOM_memory.h $(SRCROOT)/MOM6/src/framework/MOM_memory_macros.h MOM_domains.o MOM_dyn_horgrid.o MOM_error_handler.o MOM_file_parser.o MOM_grid.o MOM_open_boundary.o MOM_time_manager.o MOM_unit_scaling.o MOM_verticalGrid.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)/MOM6/config_src/memory/dynamic_nonsymmetric -I$(SRCROOT)/MOM6/src/framework	$(SRCROOT)/MOM6/src/user/shelfwave_initialization.F90
sloshing_initialization.o: $(SRCROOT)/MOM6/src/user/sloshing_initialization.F90 $(SRCROOT)/MOM6/config_src/memory/dynamic_nonsymmetric/MOM_memory.h $(SRCROOT)/MOM6/src/framework/MOM_memory_macros.h $(SRCROOT)/MOM6/src/framework/version_variable.h MOM_domains.o MOM_dyn_horgrid.o MOM_error_handler.o MOM_file_parser.o MOM_get_input.o MOM_grid.o MOM_sponge.o MOM_tracer_registry.o MOM_unit_scaling.o MOM_variables.o MOM_verticalGrid.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)/MOM6/config_src/memory/dynamic_nonsymmetric -I$(SRCROOT)/MOM6/src/framework	$(SRCROOT)/MOM6/src/user/sloshing_initialization.F90
soliton_initialization.o: $(SRCROOT)/MOM6/src/user/soliton_initialization.F90 $(SRCROOT)/MOM6/config_src/memory/dynamic_nonsymmetric/MOM_memory.h $(SRCROOT)/MOM6/src/framework/MOM_memory_macros.h MOM_error_handler.o MOM_file_parser.o MOM_get_input.o MOM_grid.o MOM_unit_scaling.o MOM_variables.o MOM_verticalGrid.o regrid_consts.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)/MOM6/config_src/memory/dynamic_nonsymmetric -I$(SRCROOT)/MOM6/src/framework	$(SRCROOT)/MOM6/src/user/soliton_initialization.F90
stochastic_physics.o: $(SRCROOT)/MOM6/config_src/external/stochastic_physics/stochastic_physics.F90 MOM_error_handler.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)/MOM6/config_src/external/stochastic_physics/stochastic_physics.F90
supercritical_initialization.o: $(SRCROOT)/MOM6/src/user/supercritical_initialization.F90 $(SRCROOT)/MOM6/config_src/memory/dynamic_nonsymmetric/MOM_memory.h $(SRCROOT)/MOM6/src/framework/MOM_memory_macros.h MOM_dyn_horgrid.o MOM_error_handler.o MOM_file_parser.o MOM_grid.o MOM_open_boundary.o MOM_time_manager.o MOM_unit_scaling.o MOM_verticalGrid.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)/MOM6/config_src/memory/dynamic_nonsymmetric -I$(SRCROOT)/MOM6/src/framework	$(SRCROOT)/MOM6/src/user/supercritical_initialization.F90
tidal_bay_initialization.o: $(SRCROOT)/MOM6/src/user/tidal_bay_initialization.F90 $(SRCROOT)/MOM6/config_src/memory/dynamic_nonsymmetric/MOM_memory.h $(SRCROOT)/MOM6/src/framework/MOM_memory_macros.h MOM_coms.o MOM_dyn_horgrid.o MOM_error_handler.o MOM_file_parser.o MOM_grid.o MOM_open_boundary.o MOM_unit_scaling.o MOM_verticalGrid.o MOM_time_manager.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)/MOM6/config_src/memory/dynamic_nonsymmetric -I$(SRCROOT)/MOM6/src/framework	$(SRCROOT)/MOM6/src/user/tidal_bay_initialization.F90
tracer_example.o: $(SRCROOT)/MOM6/src/tracer/tracer_example.F90 $(SRCROOT)/MOM6/config_src/memory/dynamic_nonsymmetric/MOM_memory.h $(SRCROOT)/MOM6/src/framework/MOM_memory_macros.h $(SRCROOT)/MOM6/src/framework/version_variable.h MOM_coms.o MOM_coupler_types.o MOM_diag_mediator.o MOM_error_handler.o MOM_file_parser.o MOM_forcing_type.o MOM_grid.o MOM_hor_index.o MOM_io.o MOM_open_boundary.o MOM_restart.o MOM_spatial_means.o MOM_sponge.o MOM_time_manager.o MOM_tracer_registry.o MOM_unit_scaling.o MOM_variables.o MOM_verticalGrid.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)/MOM6/config_src/memory/dynamic_nonsymmetric -I$(SRCROOT)/MOM6/src/framework	$(SRCROOT)/MOM6/src/tracer/tracer_example.F90
user_change_diffusivity.o: $(SRCROOT)/MOM6/src/user/user_change_diffusivity.F90 $(SRCROOT)/MOM6/config_src/memory/dynamic_nonsymmetric/MOM_memory.h $(SRCROOT)/MOM6/src/framework/MOM_memory_macros.h $(SRCROOT)/MOM6/src/framework/version_variable.h MOM_diag_mediator.o MOM_error_handler.o MOM_file_parser.o MOM_grid.o MOM_unit_scaling.o MOM_variables.o MOM_verticalGrid.o MOM_EOS.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)/MOM6/config_src/memory/dynamic_nonsymmetric -I$(SRCROOT)/MOM6/src/framework	$(SRCROOT)/MOM6/src/user/user_change_diffusivity.F90
user_initialization.o: $(SRCROOT)/MOM6/src/user/user_initialization.F90 $(SRCROOT)/MOM6/config_src/memory/dynamic_nonsymmetric/MOM_memory.h $(SRCROOT)/MOM6/src/framework/MOM_memory_macros.h $(SRCROOT)/MOM6/src/framework/version_variable.h MOM_error_handler.o MOM_dyn_horgrid.o MOM_file_parser.o MOM_get_input.o MOM_grid.o MOM_open_boundary.o MOM_sponge.o MOM_tracer_registry.o MOM_unit_scaling.o MOM_variables.o MOM_verticalGrid.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)/MOM6/config_src/memory/dynamic_nonsymmetric -I$(SRCROOT)/MOM6/src/framework	$(SRCROOT)/MOM6/src/user/user_initialization.F90
user_revise_forcing.o: $(SRCROOT)/MOM6/src/user/user_revise_forcing.F90 $(SRCROOT)/MOM6/src/framework/version_variable.h MOM_domains.o MOM_error_handler.o MOM_file_parser.o MOM_forcing_type.o MOM_grid.o MOM_io.o MOM_restart.o MOM_time_manager.o MOM_tracer_flow_control.o MOM_variables.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)/MOM6/src/framework	$(SRCROOT)/MOM6/src/user/user_revise_forcing.F90
user_shelf_init.o: $(SRCROOT)/MOM6/src/ice_shelf/user_shelf_init.F90 $(SRCROOT)/MOM6/config_src/memory/dynamic_nonsymmetric/MOM_memory.h $(SRCROOT)/MOM6/src/framework/MOM_memory_macros.h MOM_error_handler.o MOM_file_parser.o MOM_grid.o MOM_time_manager.o MOM_unit_scaling.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)/MOM6/config_src/memory/dynamic_nonsymmetric -I$(SRCROOT)/MOM6/src/framework	$(SRCROOT)/MOM6/src/ice_shelf/user_shelf_init.F90
write_ocean_obs.o: $(SRCROOT)/MOM6/config_src/external/ODA_hooks/write_ocean_obs.F90 ocean_da_types.o MOM_time_manager.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)/MOM6/config_src/external/ODA_hooks/write_ocean_obs.F90
./MOM_CVMix_ddiff.F90: $(SRCROOT)/MOM6/src/parameterizations/vertical/MOM_CVMix_ddiff.F90
	cp $(SRCROOT)/MOM6/src/parameterizations/vertical/MOM_CVMix_ddiff.F90 .
./BFB_initialization.F90: $(SRCROOT)/MOM6/src/user/BFB_initialization.F90
	cp $(SRCROOT)/MOM6/src/user/BFB_initialization.F90 .
./gsw_pt0_from_t.f90: $(SRCROOT)/MOM6/src/equation_of_state/TEOS10/gsw_pt0_from_t.f90
	cp $(SRCROOT)/MOM6/src/equation_of_state/TEOS10/gsw_pt0_from_t.f90 .
./PLM_functions.F90: $(SRCROOT)/MOM6/src/ALE/PLM_functions.F90
	cp $(SRCROOT)/MOM6/src/ALE/PLM_functions.F90 .
./MOM_sum_output.F90: $(SRCROOT)/MOM6/src/diagnostics/MOM_sum_output.F90
	cp $(SRCROOT)/MOM6/src/diagnostics/MOM_sum_output.F90 .
./advection_test_tracer.F90: $(SRCROOT)/MOM6/src/tracer/advection_test_tracer.F90
	cp $(SRCROOT)/MOM6/src/tracer/advection_test_tracer.F90 .
./gsw_t_from_ct.f90: $(SRCROOT)/MOM6/src/equation_of_state/TEOS10/gsw_t_from_ct.f90
	cp $(SRCROOT)/MOM6/src/equation_of_state/TEOS10/gsw_t_from_ct.f90 .
./MOM_regridding.F90: $(SRCROOT)/MOM6/src/ALE/MOM_regridding.F90
	cp $(SRCROOT)/MOM6/src/ALE/MOM_regridding.F90 .
./MOM_diapyc_energy_req.F90: $(SRCROOT)/MOM6/src/parameterizations/vertical/MOM_diapyc_energy_req.F90
	cp $(SRCROOT)/MOM6/src/parameterizations/vertical/MOM_diapyc_energy_req.F90 .
./gsw_gibbs.f90: $(SRCROOT)/MOM6/src/equation_of_state/TEOS10/gsw_gibbs.f90
	cp $(SRCROOT)/MOM6/src/equation_of_state/TEOS10/gsw_gibbs.f90 .
./MOM_density_integrals.F90: $(SRCROOT)/MOM6/src/core/MOM_density_integrals.F90
	cp $(SRCROOT)/MOM6/src/core/MOM_density_integrals.F90 .
./MOM_EOS_Wright_full.F90: $(SRCROOT)/MOM6/src/equation_of_state/MOM_EOS_Wright_full.F90
	cp $(SRCROOT)/MOM6/src/equation_of_state/MOM_EOS_Wright_full.F90 .
./DOME2d_initialization.F90: $(SRCROOT)/MOM6/src/user/DOME2d_initialization.F90
	cp $(SRCROOT)/MOM6/src/user/DOME2d_initialization.F90 .
./MOM_ice_shelf.F90: $(SRCROOT)/MOM6/src/ice_shelf/MOM_ice_shelf.F90
	cp $(SRCROOT)/MOM6/src/ice_shelf/MOM_ice_shelf.F90 .
./MOM_stochastics.F90: $(SRCROOT)/MOM6/src/parameterizations/stochastic/MOM_stochastics.F90
	cp $(SRCROOT)/MOM6/src/parameterizations/stochastic/MOM_stochastics.F90 .
./gsw_mod_teos10_constants.f90: $(SRCROOT)/MOM6/src/equation_of_state/TEOS10/gsw_mod_teos10_constants.f90
	cp $(SRCROOT)/MOM6/src/equation_of_state/TEOS10/gsw_mod_teos10_constants.f90 .
./MOM_unit_testing.F90: $(SRCROOT)/MOM6/src/framework/MOM_unit_testing.F90
	cp $(SRCROOT)/MOM6/src/framework/MOM_unit_testing.F90 .
./MOM_database_comms.F90: $(SRCROOT)/MOM6/config_src/external/database_comms/MOM_database_comms.F90
	cp $(SRCROOT)/MOM6/config_src/external/database_comms/MOM_database_comms.F90 .
./external_gwave_initialization.F90: $(SRCROOT)/MOM6/src/user/external_gwave_initialization.F90
	cp $(SRCROOT)/MOM6/src/user/external_gwave_initialization.F90 .
./version_variable.h: $(SRCROOT)MOM6/src/framework/version_variable.h
	cp $(SRCROOT)MOM6/src/framework/version_variable.h .
./MOM_io.F90: $(SRCROOT)/MOM6/src/framework/MOM_io.F90
	cp $(SRCROOT)/MOM6/src/framework/MOM_io.F90 .
./RGC_tracer.F90: $(SRCROOT)/MOM6/src/tracer/RGC_tracer.F90
	cp $(SRCROOT)/MOM6/src/tracer/RGC_tracer.F90 .
./MOM_random.F90: $(SRCROOT)/MOM6/src/framework/MOM_random.F90
	cp $(SRCROOT)/MOM6/src/framework/MOM_random.F90 .
./MOM_PressureForce.F90: $(SRCROOT)/MOM6/src/core/MOM_PressureForce.F90
	cp $(SRCROOT)/MOM6/src/core/MOM_PressureForce.F90 .
./MOM_tracer_flow_control.F90: $(SRCROOT)/MOM6/src/tracer/MOM_tracer_flow_control.F90
	cp $(SRCROOT)/MOM6/src/tracer/MOM_tracer_flow_control.F90 .
./tidal_bay_initialization.F90: $(SRCROOT)/MOM6/src/user/tidal_bay_initialization.F90
	cp $(SRCROOT)/MOM6/src/user/tidal_bay_initialization.F90 .
./ISOMIP_initialization.F90: $(SRCROOT)/MOM6/src/user/ISOMIP_initialization.F90
	cp $(SRCROOT)/MOM6/src/user/ISOMIP_initialization.F90 .
./MOM_load_love_numbers.F90: $(SRCROOT)/MOM6/src/parameterizations/lateral/MOM_load_love_numbers.F90
	cp $(SRCROOT)/MOM6/src/parameterizations/lateral/MOM_load_love_numbers.F90 .
./lock_exchange_initialization.F90: $(SRCROOT)/MOM6/src/user/lock_exchange_initialization.F90
	cp $(SRCROOT)/MOM6/src/user/lock_exchange_initialization.F90 .
./MOM_verticalGrid.F90: $(SRCROOT)/MOM6/src/core/MOM_verticalGrid.F90
	cp $(SRCROOT)/MOM6/src/core/MOM_verticalGrid.F90 .
./MOM_error_infra.F90: $(SRCROOT)/MOM6/config_src/infra/FMS2/MOM_error_infra.F90
	cp $(SRCROOT)/MOM6/config_src/infra/FMS2/MOM_error_infra.F90 .
./MOM_EOS_linear.F90: $(SRCROOT)/MOM6/src/equation_of_state/MOM_EOS_linear.F90
	cp $(SRCROOT)/MOM6/src/equation_of_state/MOM_EOS_linear.F90 .
./regrid_edge_values.F90: $(SRCROOT)/MOM6/src/ALE/regrid_edge_values.F90
	cp $(SRCROOT)/MOM6/src/ALE/regrid_edge_values.F90 .
./ISOMIP_tracer.F90: $(SRCROOT)/MOM6/src/tracer/ISOMIP_tracer.F90
	cp $(SRCROOT)/MOM6/src/tracer/ISOMIP_tracer.F90 .
./gsw_mod_specvol_coefficients.f90: $(SRCROOT)/MOM6/src/equation_of_state/TEOS10/gsw_mod_specvol_coefficients.f90
	cp $(SRCROOT)/MOM6/src/equation_of_state/TEOS10/gsw_mod_specvol_coefficients.f90 .
./MOM_dyn_horgrid.F90: $(SRCROOT)/MOM6/src/framework/MOM_dyn_horgrid.F90
	cp $(SRCROOT)/MOM6/src/framework/MOM_dyn_horgrid.F90 .
./dyed_obcs_initialization.F90: $(SRCROOT)/MOM6/src/user/dyed_obcs_initialization.F90
	cp $(SRCROOT)/MOM6/src/user/dyed_obcs_initialization.F90 .
./MOM_grid.F90: $(SRCROOT)/MOM6/src/core/MOM_grid.F90
	cp $(SRCROOT)/MOM6/src/core/MOM_grid.F90 .
./user_change_diffusivity.F90: $(SRCROOT)/MOM6/src/user/user_change_diffusivity.F90
	cp $(SRCROOT)/MOM6/src/user/user_change_diffusivity.F90 .
./MOM_memory.h: $(SRCROOT)/MOM6/config_src/memory/dynamic_nonsymmetric/MOM_memory.h
	cp $(SRCROOT)/MOM6/config_src/memory/dynamic_nonsymmetric/MOM_memory.h .
./MOM_grid_initialize.F90: $(SRCROOT)/MOM6/src/initialization/MOM_grid_initialize.F90
	cp $(SRCROOT)/MOM6/src/initialization/MOM_grid_initialize.F90 .
./get_stochy_pattern.F90: $(SRCROOT)/MOM6/config_src/external/stochastic_physics/get_stochy_pattern.F90
	cp $(SRCROOT)/MOM6/config_src/external/stochastic_physics/get_stochy_pattern.F90 .
./MOM_transcribe_grid.F90: $(SRCROOT)/MOM6/src/core/MOM_transcribe_grid.F90
	cp $(SRCROOT)/MOM6/src/core/MOM_transcribe_grid.F90 .
./MOM_array_transform.F90: $(SRCROOT)/MOM6/src/framework/MOM_array_transform.F90
	cp $(SRCROOT)/MOM6/src/framework/MOM_array_transform.F90 .
./MOM_interp_infra.F90: $(SRCROOT)/MOM6/config_src/infra/FMS2/MOM_interp_infra.F90
	cp $(SRCROOT)/MOM6/config_src/infra/FMS2/MOM_interp_infra.F90 .
./gsw_specvol_first_derivatives.f90: $(SRCROOT)/MOM6/src/equation_of_state/TEOS10/gsw_specvol_first_derivatives.f90
	cp $(SRCROOT)/MOM6/src/equation_of_state/TEOS10/gsw_specvol_first_derivatives.f90 .
./user_shelf_init.F90: $(SRCROOT)/MOM6/src/ice_shelf/user_shelf_init.F90
	cp $(SRCROOT)/MOM6/src/ice_shelf/user_shelf_init.F90 .
./cvmix_math.F90: $(SRCROOT)/MOM6/src/parameterizations/CVmix/cvmix_math.F90
	cp $(SRCROOT)/MOM6/src/parameterizations/CVmix/cvmix_math.F90 .
./benchmark_initialization.F90: $(SRCROOT)/MOM6/src/user/benchmark_initialization.F90
	cp $(SRCROOT)/MOM6/src/user/benchmark_initialization.F90 .
./MOM_ice_shelf_state.F90: $(SRCROOT)/MOM6/src/ice_shelf/MOM_ice_shelf_state.F90
	cp $(SRCROOT)/MOM6/src/ice_shelf/MOM_ice_shelf_state.F90 .
./MOM_bkgnd_mixing.F90: $(SRCROOT)/MOM6/src/parameterizations/vertical/MOM_bkgnd_mixing.F90
	cp $(SRCROOT)/MOM6/src/parameterizations/vertical/MOM_bkgnd_mixing.F90 .
./version_variable.h: $(SRCROOT)/MOM6/src/framework/version_variable.h
	cp $(SRCROOT)/MOM6/src/framework/version_variable.h .
./MOM_unit_scaling.F90: $(SRCROOT)/MOM6/src/framework/MOM_unit_scaling.F90
	cp $(SRCROOT)/MOM6/src/framework/MOM_unit_scaling.F90 .
./MOM_coms.F90: $(SRCROOT)/MOM6/src/framework/MOM_coms.F90
	cp $(SRCROOT)/MOM6/src/framework/MOM_coms.F90 .
./MOM_tracer_diabatic.F90: $(SRCROOT)/MOM6/src/tracer/MOM_tracer_diabatic.F90
	cp $(SRCROOT)/MOM6/src/tracer/MOM_tracer_diabatic.F90 .
./MOM_safe_alloc.F90: $(SRCROOT)/MOM6/src/framework/MOM_safe_alloc.F90
	cp $(SRCROOT)/MOM6/src/framework/MOM_safe_alloc.F90 .
./MOM_marine_ice.F90: $(SRCROOT)/MOM6/src/ice_shelf/MOM_marine_ice.F90
	cp $(SRCROOT)/MOM6/src/ice_shelf/MOM_marine_ice.F90 .
./MOM_netcdf.F90: $(SRCROOT)/MOM6/src/framework/MOM_netcdf.F90
	cp $(SRCROOT)/MOM6/src/framework/MOM_netcdf.F90 .
./MOM_hybgen_regrid.F90: $(SRCROOT)/MOM6/src/ALE/MOM_hybgen_regrid.F90
	cp $(SRCROOT)/MOM6/src/ALE/MOM_hybgen_regrid.F90 .
./MOM_oda_incupd.F90: $(SRCROOT)/MOM6/src/ocean_data_assim/MOM_oda_incupd.F90
	cp $(SRCROOT)/MOM6/src/ocean_data_assim/MOM_oda_incupd.F90 .
./basin_builder.F90: $(SRCROOT)/MOM6/src/user/basin_builder.F90
	cp $(SRCROOT)/MOM6/src/user/basin_builder.F90 .
./MOM_dynamics_split_RK2.F90: $(SRCROOT)/MOM6/src/core/MOM_dynamics_split_RK2.F90
	cp $(SRCROOT)/MOM6/src/core/MOM_dynamics_split_RK2.F90 .
./MOM_wave_interface.F90: $(SRCROOT)/MOM6/src/user/MOM_wave_interface.F90
	cp $(SRCROOT)/MOM6/src/user/MOM_wave_interface.F90 .
./MOM_obsolete_diagnostics.F90: $(SRCROOT)/MOM6/src/diagnostics/MOM_obsolete_diagnostics.F90
	cp $(SRCROOT)/MOM6/src/diagnostics/MOM_obsolete_diagnostics.F90 .
./gsw_pt_from_t.f90: $(SRCROOT)/MOM6/src/equation_of_state/TEOS10/gsw_pt_from_t.f90
	cp $(SRCROOT)/MOM6/src/equation_of_state/TEOS10/gsw_pt_from_t.f90 .
./MOM_ice_shelf_diag_mediator.F90: $(SRCROOT)/MOM6/src/ice_shelf/MOM_ice_shelf_diag_mediator.F90
	cp $(SRCROOT)/MOM6/src/ice_shelf/MOM_ice_shelf_diag_mediator.F90 .
./gsw_rho_second_derivatives.f90: $(SRCROOT)/MOM6/src/equation_of_state/TEOS10/gsw_rho_second_derivatives.f90
	cp $(SRCROOT)/MOM6/src/equation_of_state/TEOS10/gsw_rho_second_derivatives.f90 .
./MOM_diagnostics.F90: $(SRCROOT)/MOM6/src/diagnostics/MOM_diagnostics.F90
	cp $(SRCROOT)/MOM6/src/diagnostics/MOM_diagnostics.F90 .
./MOM_set_diffusivity.F90: $(SRCROOT)/MOM6/src/parameterizations/vertical/MOM_set_diffusivity.F90
	cp $(SRCROOT)/MOM6/src/parameterizations/vertical/MOM_set_diffusivity.F90 .
./MOM_open_boundary.F90: $(SRCROOT)/MOM6/src/core/MOM_open_boundary.F90
	cp $(SRCROOT)/MOM6/src/core/MOM_open_boundary.F90 .
./DOME_initialization.F90: $(SRCROOT)/MOM6/src/user/DOME_initialization.F90
	cp $(SRCROOT)/MOM6/src/user/DOME_initialization.F90 .
./dumbbell_initialization.F90: $(SRCROOT)/MOM6/src/user/dumbbell_initialization.F90
	cp $(SRCROOT)/MOM6/src/user/dumbbell_initialization.F90 .
./MOM_EOS_Wright.F90: $(SRCROOT)/MOM6/src/equation_of_state/MOM_EOS_Wright.F90
	cp $(SRCROOT)/MOM6/src/equation_of_state/MOM_EOS_Wright.F90 .
./boundary_impulse_tracer.F90: $(SRCROOT)/MOM6/src/tracer/boundary_impulse_tracer.F90
	cp $(SRCROOT)/MOM6/src/tracer/boundary_impulse_tracer.F90 .
./gsw_ct_freezing_exact.f90: $(SRCROOT)/MOM6/src/equation_of_state/TEOS10/gsw_ct_freezing_exact.f90
	cp $(SRCROOT)/MOM6/src/equation_of_state/TEOS10/gsw_ct_freezing_exact.f90 .
./MOM_data_override_infra.F90: $(SRCROOT)/MOM6/config_src/infra/FMS2/MOM_data_override_infra.F90
	cp $(SRCROOT)/MOM6/config_src/infra/FMS2/MOM_data_override_infra.F90 .
./gsw_gibbs_pt0_pt0.f90: $(SRCROOT)/MOM6/src/equation_of_state/TEOS10/gsw_gibbs_pt0_pt0.f90
	cp $(SRCROOT)/MOM6/src/equation_of_state/TEOS10/gsw_gibbs_pt0_pt0.f90 .
./MOM_tidal_mixing.F90: $(SRCROOT)/MOM6/src/parameterizations/vertical/MOM_tidal_mixing.F90
	cp $(SRCROOT)/MOM6/src/parameterizations/vertical/MOM_tidal_mixing.F90 .
./MOM_regularize_layers.F90: $(SRCROOT)/MOM6/src/parameterizations/vertical/MOM_regularize_layers.F90
	cp $(SRCROOT)/MOM6/src/parameterizations/vertical/MOM_regularize_layers.F90 .
./MOM_sponge.F90: $(SRCROOT)/MOM6/src/parameterizations/vertical/MOM_sponge.F90
	cp $(SRCROOT)/MOM6/src/parameterizations/vertical/MOM_sponge.F90 .
./cvmix_tidal.F90: $(SRCROOT)/MOM6/src/parameterizations/CVmix/cvmix_tidal.F90
	cp $(SRCROOT)/MOM6/src/parameterizations/CVmix/cvmix_tidal.F90 .
./MOM_CVMix_conv.F90: $(SRCROOT)/MOM6/src/parameterizations/vertical/MOM_CVMix_conv.F90
	cp $(SRCROOT)/MOM6/src/parameterizations/vertical/MOM_CVMix_conv.F90 .
./MOM_entrain_diffusive.F90: $(SRCROOT)/MOM6/src/parameterizations/vertical/MOM_entrain_diffusive.F90
	cp $(SRCROOT)/MOM6/src/parameterizations/vertical/MOM_entrain_diffusive.F90 .
./pseudo_salt_tracer.F90: $(SRCROOT)/MOM6/src/tracer/pseudo_salt_tracer.F90
	cp $(SRCROOT)/MOM6/src/tracer/pseudo_salt_tracer.F90 .
./gsw_mod_freezing_poly_coefficients.f90: $(SRCROOT)/MOM6/src/equation_of_state/TEOS10/gsw_mod_freezing_poly_coefficients.f90
	cp $(SRCROOT)/MOM6/src/equation_of_state/TEOS10/gsw_mod_freezing_poly_coefficients.f90 .
./MOM_hor_index.F90: $(SRCROOT)/MOM6/src/framework/MOM_hor_index.F90
	cp $(SRCROOT)/MOM6/src/framework/MOM_hor_index.F90 .
./MOM_EOS_Roquet_rho.F90: $(SRCROOT)/MOM6/src/equation_of_state/MOM_EOS_Roquet_rho.F90
	cp $(SRCROOT)/MOM6/src/equation_of_state/MOM_EOS_Roquet_rho.F90 .
./MOM_porous_barriers.F90: $(SRCROOT)/MOM6/src/core/MOM_porous_barriers.F90
	cp $(SRCROOT)/MOM6/src/core/MOM_porous_barriers.F90 .
./MOM_intrinsic_functions.F90: $(SRCROOT)/MOM6/src/framework/MOM_intrinsic_functions.F90
	cp $(SRCROOT)/MOM6/src/framework/MOM_intrinsic_functions.F90 .
./MOM_EOS.F90: $(SRCROOT)/MOM6/src/equation_of_state/MOM_EOS.F90
	cp $(SRCROOT)/MOM6/src/equation_of_state/MOM_EOS.F90 .
./gsw_gibbs_ice.f90: $(SRCROOT)/MOM6/src/equation_of_state/TEOS10/gsw_gibbs_ice.f90
	cp $(SRCROOT)/MOM6/src/equation_of_state/TEOS10/gsw_gibbs_ice.f90 .
./MOM_ice_shelf_dynamics.F90: $(SRCROOT)/MOM6/src/ice_shelf/MOM_ice_shelf_dynamics.F90
	cp $(SRCROOT)/MOM6/src/ice_shelf/MOM_ice_shelf_dynamics.F90 .
./gsw_entropy_part.f90: $(SRCROOT)/MOM6/src/equation_of_state/TEOS10/gsw_entropy_part.f90
	cp $(SRCROOT)/MOM6/src/equation_of_state/TEOS10/gsw_entropy_part.f90 .
./MOM_EOS_Jackett06.F90: $(SRCROOT)/MOM6/src/equation_of_state/MOM_EOS_Jackett06.F90
	cp $(SRCROOT)/MOM6/src/equation_of_state/MOM_EOS_Jackett06.F90 .
./MOM_TFreeze.F90: $(SRCROOT)/MOM6/src/equation_of_state/MOM_TFreeze.F90
	cp $(SRCROOT)/MOM6/src/equation_of_state/MOM_TFreeze.F90 .
./gsw_ct_from_t.f90: $(SRCROOT)/MOM6/src/equation_of_state/TEOS10/gsw_ct_from_t.f90
	cp $(SRCROOT)/MOM6/src/equation_of_state/TEOS10/gsw_ct_from_t.f90 .
./P1M_functions.F90: $(SRCROOT)/MOM6/src/ALE/P1M_functions.F90
	cp $(SRCROOT)/MOM6/src/ALE/P1M_functions.F90 .
./gsw_rho.f90: $(SRCROOT)/MOM6/src/equation_of_state/TEOS10/gsw_rho.f90
	cp $(SRCROOT)/MOM6/src/equation_of_state/TEOS10/gsw_rho.f90 .
./MOM_io_file.F90: $(SRCROOT)/MOM6/src/framework/MOM_io_file.F90
	cp $(SRCROOT)/MOM6/src/framework/MOM_io_file.F90 .
./PPM_functions.F90: $(SRCROOT)/MOM6/src/ALE/PPM_functions.F90
	cp $(SRCROOT)/MOM6/src/ALE/PPM_functions.F90 .
./MOM_coord_initialization.F90: $(SRCROOT)/MOM6/src/initialization/MOM_coord_initialization.F90
	cp $(SRCROOT)/MOM6/src/initialization/MOM_coord_initialization.F90 .
./user_initialization.F90: $(SRCROOT)/MOM6/src/user/user_initialization.F90
	cp $(SRCROOT)/MOM6/src/user/user_initialization.F90 .
./generic_tracer_utils.F90: $(SRCROOT)/MOM6/config_src/external/GFDL_ocean_BGC/generic_tracer_utils.F90
	cp $(SRCROOT)/MOM6/config_src/external/GFDL_ocean_BGC/generic_tracer_utils.F90 .
./gsw_chem_potential_water_t_exact.f90: $(SRCROOT)/MOM6/src/equation_of_state/TEOS10/gsw_chem_potential_water_t_exact.f90
	cp $(SRCROOT)/MOM6/src/equation_of_state/TEOS10/gsw_chem_potential_water_t_exact.f90 .
./MOM_fixed_initialization.F90: $(SRCROOT)/MOM6/src/initialization/MOM_fixed_initialization.F90
	cp $(SRCROOT)/MOM6/src/initialization/MOM_fixed_initialization.F90 .
./MOM_error_handler.F90: $(SRCROOT)/MOM6/src/framework/MOM_error_handler.F90
	cp $(SRCROOT)/MOM6/src/framework/MOM_error_handler.F90 .
./Rossby_front_2d_initialization.F90: $(SRCROOT)/MOM6/src/user/Rossby_front_2d_initialization.F90
	cp $(SRCROOT)/MOM6/src/user/Rossby_front_2d_initialization.F90 .
./MOM_cpu_clock.F90: $(SRCROOT)/MOM6/src/framework/MOM_cpu_clock.F90
	cp $(SRCROOT)/MOM6/src/framework/MOM_cpu_clock.F90 .
./MOM_Zanna_Bolton.F90: $(SRCROOT)/MOM6/src/parameterizations/lateral/MOM_Zanna_Bolton.F90
	cp $(SRCROOT)/MOM6/src/parameterizations/lateral/MOM_Zanna_Bolton.F90 .
./circle_obcs_initialization.F90: $(SRCROOT)/MOM6/src/user/circle_obcs_initialization.F90
	cp $(SRCROOT)/MOM6/src/user/circle_obcs_initialization.F90 .
./gsw_ct_freezing_poly.f90: $(SRCROOT)/MOM6/src/equation_of_state/TEOS10/gsw_ct_freezing_poly.f90
	cp $(SRCROOT)/MOM6/src/equation_of_state/TEOS10/gsw_ct_freezing_poly.f90 .
./MOM_PressureForce_FV.F90: $(SRCROOT)/MOM6/src/core/MOM_PressureForce_FV.F90
	cp $(SRCROOT)/MOM6/src/core/MOM_PressureForce_FV.F90 .
./MOM_variables.F90: $(SRCROOT)/MOM6/src/core/MOM_variables.F90
	cp $(SRCROOT)/MOM6/src/core/MOM_variables.F90 .
./MOM_CFC_cap.F90: $(SRCROOT)/MOM6/src/tracer/MOM_CFC_cap.F90
	cp $(SRCROOT)/MOM6/src/tracer/MOM_CFC_cap.F90 .
./MOM.F90: $(SRCROOT)/MOM6/src/core/MOM.F90
	cp $(SRCROOT)/MOM6/src/core/MOM.F90 .
./gsw_t_freezing_poly.f90: $(SRCROOT)/MOM6/src/equation_of_state/TEOS10/gsw_t_freezing_poly.f90
	cp $(SRCROOT)/MOM6/src/equation_of_state/TEOS10/gsw_t_freezing_poly.f90 .
./dyed_channel_initialization.F90: $(SRCROOT)/MOM6/src/user/dyed_channel_initialization.F90
	cp $(SRCROOT)/MOM6/src/user/dyed_channel_initialization.F90 .
./MOM_continuity_PPM.F90: $(SRCROOT)/MOM6/src/core/MOM_continuity_PPM.F90
	cp $(SRCROOT)/MOM6/src/core/MOM_continuity_PPM.F90 .
./MOM_EOS_UNESCO.F90: $(SRCROOT)/MOM6/src/equation_of_state/MOM_EOS_UNESCO.F90
	cp $(SRCROOT)/MOM6/src/equation_of_state/MOM_EOS_UNESCO.F90 .
./cvmix_kinds_and_types.F90: $(SRCROOT)/MOM6/src/parameterizations/CVmix/cvmix_kinds_and_types.F90
	cp $(SRCROOT)/MOM6/src/parameterizations/CVmix/cvmix_kinds_and_types.F90 .
./MOM_energetic_PBL.F90: $(SRCROOT)/MOM6/src/parameterizations/vertical/MOM_energetic_PBL.F90
	cp $(SRCROOT)/MOM6/src/parameterizations/vertical/MOM_energetic_PBL.F90 .
./database_client_interface.F90: $(SRCROOT)/MOM6/config_src/external/database_comms/database_client_interface.F90
	cp $(SRCROOT)/MOM6/config_src/external/database_comms/database_client_interface.F90 .
./tracer_example.F90: $(SRCROOT)/MOM6/src/tracer/tracer_example.F90
	cp $(SRCROOT)/MOM6/src/tracer/tracer_example.F90 .
./gsw_mod_kinds.f90: $(SRCROOT)/MOM6/src/equation_of_state/TEOS10/gsw_mod_kinds.f90
	cp $(SRCROOT)/MOM6/src/equation_of_state/TEOS10/gsw_mod_kinds.f90 .
./MOM_thickness_diffuse.F90: $(SRCROOT)/MOM6/src/parameterizations/lateral/MOM_thickness_diffuse.F90
	cp $(SRCROOT)/MOM6/src/parameterizations/lateral/MOM_thickness_diffuse.F90 .
./MOM_checksum_packages.F90: $(SRCROOT)/MOM6/src/core/MOM_checksum_packages.F90
	cp $(SRCROOT)/MOM6/src/core/MOM_checksum_packages.F90 .
./MOM_remapping.F90: $(SRCROOT)/MOM6/src/ALE/MOM_remapping.F90
	cp $(SRCROOT)/MOM6/src/ALE/MOM_remapping.F90 .
./MOM_diag_remap.F90: $(SRCROOT)/MOM6/src/framework/MOM_diag_remap.F90
	cp $(SRCROOT)/MOM6/src/framework/MOM_diag_remap.F90 .
./MOM_hybgen_remap.F90: $(SRCROOT)/MOM6/src/ALE/MOM_hybgen_remap.F90
	cp $(SRCROOT)/MOM6/src/ALE/MOM_hybgen_remap.F90 .
./coord_sigma.F90: $(SRCROOT)/MOM6/src/ALE/coord_sigma.F90
	cp $(SRCROOT)/MOM6/src/ALE/coord_sigma.F90 .
./MOM_oda_driver.F90: $(SRCROOT)/MOM6/src/ocean_data_assim/MOM_oda_driver.F90
	cp $(SRCROOT)/MOM6/src/ocean_data_assim/MOM_oda_driver.F90 .
./MOM_unit_tests.F90: $(SRCROOT)/MOM6/src/core/MOM_unit_tests.F90
	cp $(SRCROOT)/MOM6/src/core/MOM_unit_tests.F90 .
./stochastic_physics.F90: $(SRCROOT)/MOM6/config_src/external/stochastic_physics/stochastic_physics.F90
	cp $(SRCROOT)/MOM6/config_src/external/stochastic_physics/stochastic_physics.F90 .
./MOM_ensemble_manager.F90: $(SRCROOT)/MOM6/src/framework/MOM_ensemble_manager.F90
	cp $(SRCROOT)/MOM6/src/framework/MOM_ensemble_manager.F90 .
./cvmix_ddiff.F90: $(SRCROOT)/MOM6/src/parameterizations/CVmix/cvmix_ddiff.F90
	cp $(SRCROOT)/MOM6/src/parameterizations/CVmix/cvmix_ddiff.F90 .
./sloshing_initialization.F90: $(SRCROOT)/MOM6/src/user/sloshing_initialization.F90
	cp $(SRCROOT)/MOM6/src/user/sloshing_initialization.F90 .
./MOM_diag_manager_infra.F90: $(SRCROOT)/MOM6/config_src/infra/FMS2/MOM_diag_manager_infra.F90
	cp $(SRCROOT)/MOM6/config_src/infra/FMS2/MOM_diag_manager_infra.F90 .
./gsw_specvol.f90: $(SRCROOT)/MOM6/src/equation_of_state/TEOS10/gsw_specvol.f90
	cp $(SRCROOT)/MOM6/src/equation_of_state/TEOS10/gsw_specvol.f90 .
./Idealized_Hurricane.F90: $(SRCROOT)/MOM6/src/user/Idealized_Hurricane.F90
	cp $(SRCROOT)/MOM6/src/user/Idealized_Hurricane.F90 .
./MOM_coms_infra.F90: $(SRCROOT)/MOM6/config_src/infra/FMS2/MOM_coms_infra.F90
	cp $(SRCROOT)/MOM6/config_src/infra/FMS2/MOM_coms_infra.F90 .
./gsw_mod_gibbs_ice_coefficients.f90: $(SRCROOT)/MOM6/src/equation_of_state/TEOS10/gsw_mod_gibbs_ice_coefficients.f90
	cp $(SRCROOT)/MOM6/src/equation_of_state/TEOS10/gsw_mod_gibbs_ice_coefficients.f90 .
./MOM_dynamics_unsplit.F90: $(SRCROOT)/MOM6/src/core/MOM_dynamics_unsplit.F90
	cp $(SRCROOT)/MOM6/src/core/MOM_dynamics_unsplit.F90 .
./gsw_mod_error_functions.f90: $(SRCROOT)/MOM6/src/equation_of_state/TEOS10/gsw_mod_error_functions.f90
	cp $(SRCROOT)/MOM6/src/equation_of_state/TEOS10/gsw_mod_error_functions.f90 .
./cvmix_kpp.F90: $(SRCROOT)/MOM6/src/parameterizations/CVmix/cvmix_kpp.F90
	cp $(SRCROOT)/MOM6/src/parameterizations/CVmix/cvmix_kpp.F90 .
./ideal_age_example.F90: $(SRCROOT)/MOM6/src/tracer/ideal_age_example.F90
	cp $(SRCROOT)/MOM6/src/tracer/ideal_age_example.F90 .
./cvmix_put_get.F90: $(SRCROOT)/MOM6/src/parameterizations/CVmix/cvmix_put_get.F90
	cp $(SRCROOT)/MOM6/src/parameterizations/CVmix/cvmix_put_get.F90 .
./coord_zlike.F90: $(SRCROOT)/MOM6/src/ALE/coord_zlike.F90
	cp $(SRCROOT)/MOM6/src/ALE/coord_zlike.F90 .
./gsw_sp_from_sr.f90: $(SRCROOT)/MOM6/src/equation_of_state/TEOS10/gsw_sp_from_sr.f90
	cp $(SRCROOT)/MOM6/src/equation_of_state/TEOS10/gsw_sp_from_sr.f90 .
./gsw_mod_toolbox.f90: $(SRCROOT)/MOM6/src/equation_of_state/TEOS10/gsw_mod_toolbox.f90
	cp $(SRCROOT)/MOM6/src/equation_of_state/TEOS10/gsw_mod_toolbox.f90 .
./MOM_domain_infra.F90: $(SRCROOT)/MOM6/config_src/infra/FMS2/MOM_domain_infra.F90
	cp $(SRCROOT)/MOM6/config_src/infra/FMS2/MOM_domain_infra.F90 .
./MOM_get_input.F90: $(SRCROOT)/MOM6/src/framework/MOM_get_input.F90
	cp $(SRCROOT)/MOM6/src/framework/MOM_get_input.F90 .
./cvmix_utils.F90: $(SRCROOT)/MOM6/src/parameterizations/CVmix/cvmix_utils.F90
	cp $(SRCROOT)/MOM6/src/parameterizations/CVmix/cvmix_utils.F90 .
./user_revise_forcing.F90: $(SRCROOT)/MOM6/src/user/user_revise_forcing.F90
	cp $(SRCROOT)/MOM6/src/user/user_revise_forcing.F90 .
./MOM_particles_types.F90: $(SRCROOT)/MOM6/config_src/external/drifters/MOM_particles_types.F90
	cp $(SRCROOT)/MOM6/config_src/external/drifters/MOM_particles_types.F90 .
./MOM_internal_tide_input.F90: $(SRCROOT)/MOM6/src/parameterizations/vertical/MOM_internal_tide_input.F90
	cp $(SRCROOT)/MOM6/src/parameterizations/vertical/MOM_internal_tide_input.F90 .
./MOM_interface_heights.F90: $(SRCROOT)/MOM6/src/core/MOM_interface_heights.F90
	cp $(SRCROOT)/MOM6/src/core/MOM_interface_heights.F90 .
./P3M_functions.F90: $(SRCROOT)/MOM6/src/ALE/P3M_functions.F90
	cp $(SRCROOT)/MOM6/src/ALE/P3M_functions.F90 .
./BFB_surface_forcing.F90: $(SRCROOT)/MOM6/src/user/BFB_surface_forcing.F90
	cp $(SRCROOT)/MOM6/src/user/BFB_surface_forcing.F90 .
./MOM_hybgen_unmix.F90: $(SRCROOT)/MOM6/src/ALE/MOM_hybgen_unmix.F90
	cp $(SRCROOT)/MOM6/src/ALE/MOM_hybgen_unmix.F90 .
./MOM_string_functions.F90: $(SRCROOT)/MOM6/src/framework/MOM_string_functions.F90
	cp $(SRCROOT)/MOM6/src/framework/MOM_string_functions.F90 .
./MOM_diag_mediator.F90: $(SRCROOT)/MOM6/src/framework/MOM_diag_mediator.F90
	cp $(SRCROOT)/MOM6/src/framework/MOM_diag_mediator.F90 .
./remapping_attic.F90: $(SRCROOT)/MOM6/src/ALE/remapping_attic.F90
	cp $(SRCROOT)/MOM6/src/ALE/remapping_attic.F90 .
./cvmix_convection.F90: $(SRCROOT)/MOM6/src/parameterizations/CVmix/cvmix_convection.F90
	cp $(SRCROOT)/MOM6/src/parameterizations/CVmix/cvmix_convection.F90 .
./PQM_functions.F90: $(SRCROOT)/MOM6/src/ALE/PQM_functions.F90
	cp $(SRCROOT)/MOM6/src/ALE/PQM_functions.F90 .
./MOM_obsolete_params.F90: $(SRCROOT)/MOM6/src/diagnostics/MOM_obsolete_params.F90
	cp $(SRCROOT)/MOM6/src/diagnostics/MOM_obsolete_params.F90 .
./soliton_initialization.F90: $(SRCROOT)/MOM6/src/user/soliton_initialization.F90
	cp $(SRCROOT)/MOM6/src/user/soliton_initialization.F90 .
./MOM_CoriolisAdv.F90: $(SRCROOT)/MOM6/src/core/MOM_CoriolisAdv.F90
	cp $(SRCROOT)/MOM6/src/core/MOM_CoriolisAdv.F90 .
./MOM_constants.F90: $(SRCROOT)/MOM6/config_src/infra/FMS2/MOM_constants.F90
	cp $(SRCROOT)/MOM6/config_src/infra/FMS2/MOM_constants.F90 .
./MOM_barotropic.F90: $(SRCROOT)/MOM6/src/core/MOM_barotropic.F90
	cp $(SRCROOT)/MOM6/src/core/MOM_barotropic.F90 .
./MOM_internal_tides.F90: $(SRCROOT)/MOM6/src/parameterizations/lateral/MOM_internal_tides.F90
	cp $(SRCROOT)/MOM6/src/parameterizations/lateral/MOM_internal_tides.F90 .
./posix.F90: $(SRCROOT)/MOM6/src/framework/posix.F90
	cp $(SRCROOT)/MOM6/src/framework/posix.F90 .
./MOM_io_infra.F90: $(SRCROOT)/MOM6/config_src/infra/FMS2/MOM_io_infra.F90
	cp $(SRCROOT)/MOM6/config_src/infra/FMS2/MOM_io_infra.F90 .
./MOM_state_initialization.F90: $(SRCROOT)/MOM6/src/initialization/MOM_state_initialization.F90
	cp $(SRCROOT)/MOM6/src/initialization/MOM_state_initialization.F90 .
./regrid_solvers.F90: $(SRCROOT)/MOM6/src/ALE/regrid_solvers.F90
	cp $(SRCROOT)/MOM6/src/ALE/regrid_solvers.F90 .
./MOM_diabatic_aux.F90: $(SRCROOT)/MOM6/src/parameterizations/vertical/MOM_diabatic_aux.F90
	cp $(SRCROOT)/MOM6/src/parameterizations/vertical/MOM_diabatic_aux.F90 .
./MOM_offline_main.F90: $(SRCROOT)/MOM6/src/tracer/MOM_offline_main.F90
	cp $(SRCROOT)/MOM6/src/tracer/MOM_offline_main.F90 .
./MOM_ALE.F90: $(SRCROOT)/MOM6/src/ALE/MOM_ALE.F90
	cp $(SRCROOT)/MOM6/src/ALE/MOM_ALE.F90 .
./Kelvin_initialization.F90: $(SRCROOT)/MOM6/src/user/Kelvin_initialization.F90
	cp $(SRCROOT)/MOM6/src/user/Kelvin_initialization.F90 .
./MOM_PressureForce_Montgomery.F90: $(SRCROOT)/MOM6/src/core/MOM_PressureForce_Montgomery.F90
	cp $(SRCROOT)/MOM6/src/core/MOM_PressureForce_Montgomery.F90 .
./MOM_tracer_registry.F90: $(SRCROOT)/MOM6/src/tracer/MOM_tracer_registry.F90
	cp $(SRCROOT)/MOM6/src/tracer/MOM_tracer_registry.F90 .
./seamount_initialization.F90: $(SRCROOT)/MOM6/src/user/seamount_initialization.F90
	cp $(SRCROOT)/MOM6/src/user/seamount_initialization.F90 .
./MOM_full_convection.F90: $(SRCROOT)/MOM6/src/parameterizations/vertical/MOM_full_convection.F90
	cp $(SRCROOT)/MOM6/src/parameterizations/vertical/MOM_full_convection.F90 .
./MOM_controlled_forcing.F90: $(SRCROOT)/MOM6/src/user/MOM_controlled_forcing.F90
	cp $(SRCROOT)/MOM6/src/user/MOM_controlled_forcing.F90 .
./gsw_sr_from_sp.f90: $(SRCROOT)/MOM6/src/equation_of_state/TEOS10/gsw_sr_from_sp.f90
	cp $(SRCROOT)/MOM6/src/equation_of_state/TEOS10/gsw_sr_from_sp.f90 .
./MOM_shared_initialization.F90: $(SRCROOT)/MOM6/src/initialization/MOM_shared_initialization.F90
	cp $(SRCROOT)/MOM6/src/initialization/MOM_shared_initialization.F90 .
./MOM_ice_shelf_initialize.F90: $(SRCROOT)/MOM6/src/ice_shelf/MOM_ice_shelf_initialize.F90
	cp $(SRCROOT)/MOM6/src/ice_shelf/MOM_ice_shelf_initialize.F90 .
./MOM_particles.F90: $(SRCROOT)/MOM6/config_src/external/drifters/MOM_particles.F90
	cp $(SRCROOT)/MOM6/config_src/external/drifters/MOM_particles.F90 .
./cvmix_shear.F90: $(SRCROOT)/MOM6/src/parameterizations/CVmix/cvmix_shear.F90
	cp $(SRCROOT)/MOM6/src/parameterizations/CVmix/cvmix_shear.F90 .
./MOM_ensemble_manager_infra.F90: $(SRCROOT)/MOM6/config_src/infra/FMS2/MOM_ensemble_manager_infra.F90
	cp $(SRCROOT)/MOM6/config_src/infra/FMS2/MOM_ensemble_manager_infra.F90 .
./MOM_geothermal.F90: $(SRCROOT)/MOM6/src/parameterizations/vertical/MOM_geothermal.F90
	cp $(SRCROOT)/MOM6/src/parameterizations/vertical/MOM_geothermal.F90 .
./MOM_temperature_convert.F90: $(SRCROOT)/MOM6/src/equation_of_state/MOM_temperature_convert.F90
	cp $(SRCROOT)/MOM6/src/equation_of_state/MOM_temperature_convert.F90 .
./regrid_consts.F90: $(SRCROOT)/MOM6/src/ALE/regrid_consts.F90
	cp $(SRCROOT)/MOM6/src/ALE/regrid_consts.F90 .
./MOM_CVMix_shear.F90: $(SRCROOT)/MOM6/src/parameterizations/vertical/MOM_CVMix_shear.F90
	cp $(SRCROOT)/MOM6/src/parameterizations/vertical/MOM_CVMix_shear.F90 .
./MOM_tracer_initialization_from_Z.F90: $(SRCROOT)/MOM6/src/initialization/MOM_tracer_initialization_from_Z.F90
	cp $(SRCROOT)/MOM6/src/initialization/MOM_tracer_initialization_from_Z.F90 .
./MOM_continuity.F90: $(SRCROOT)/MOM6/src/core/MOM_continuity.F90
	cp $(SRCROOT)/MOM6/src/core/MOM_continuity.F90 .
./MOM_memory_macros.h: $(SRCROOT)/MOM6/src/framework/MOM_memory_macros.h
	cp $(SRCROOT)/MOM6/src/framework/MOM_memory_macros.h .
./MOM_coupler_types.F90: $(SRCROOT)/MOM6/src/framework/MOM_coupler_types.F90
	cp $(SRCROOT)/MOM6/src/framework/MOM_coupler_types.F90 .
./baroclinic_zone_initialization.F90: $(SRCROOT)/MOM6/src/user/baroclinic_zone_initialization.F90
	cp $(SRCROOT)/MOM6/src/user/baroclinic_zone_initialization.F90 .
./MOM_surface_forcing_gfdl.F90: $(SRCROOT)/MOM6/config_src/drivers/FMS_cap/MOM_surface_forcing_gfdl.F90
	cp $(SRCROOT)/MOM6/config_src/drivers/FMS_cap/MOM_surface_forcing_gfdl.F90 .
./supercritical_initialization.F90: $(SRCROOT)/MOM6/src/user/supercritical_initialization.F90
	cp $(SRCROOT)/MOM6/src/user/supercritical_initialization.F90 .
./MOM_tracer_types.F90: $(SRCROOT)/MOM6/src/tracer/MOM_tracer_types.F90
	cp $(SRCROOT)/MOM6/src/tracer/MOM_tracer_types.F90 .
./shelfwave_initialization.F90: $(SRCROOT)/MOM6/src/user/shelfwave_initialization.F90
	cp $(SRCROOT)/MOM6/src/user/shelfwave_initialization.F90 .
./MOM_lateral_mixing_coeffs.F90: $(SRCROOT)/MOM6/src/parameterizations/lateral/MOM_lateral_mixing_coeffs.F90
	cp $(SRCROOT)/MOM6/src/parameterizations/lateral/MOM_lateral_mixing_coeffs.F90 .
./coord_rho.F90: $(SRCROOT)/MOM6/src/ALE/coord_rho.F90
	cp $(SRCROOT)/MOM6/src/ALE/coord_rho.F90 .
./MOM_file_parser.F90: $(SRCROOT)/MOM6/src/framework/MOM_file_parser.F90
	cp $(SRCROOT)/MOM6/src/framework/MOM_file_parser.F90 .
./MOM_couplertype_infra.F90: $(SRCROOT)/MOM6/config_src/infra/FMS2/MOM_couplertype_infra.F90
	cp $(SRCROOT)/MOM6/config_src/infra/FMS2/MOM_couplertype_infra.F90 .
./MOM_dynamics_unsplit_RK2.F90: $(SRCROOT)/MOM6/src/core/MOM_dynamics_unsplit_RK2.F90
	cp $(SRCROOT)/MOM6/src/core/MOM_dynamics_unsplit_RK2.F90 .
./MOM_MEKE.F90: $(SRCROOT)/MOM6/src/parameterizations/lateral/MOM_MEKE.F90
	cp $(SRCROOT)/MOM6/src/parameterizations/lateral/MOM_MEKE.F90 .
./MOM_document.F90: $(SRCROOT)/MOM6/src/framework/MOM_document.F90
	cp $(SRCROOT)/MOM6/src/framework/MOM_document.F90 .
./MOM_spatial_means.F90: $(SRCROOT)/MOM6/src/diagnostics/MOM_spatial_means.F90
	cp $(SRCROOT)/MOM6/src/diagnostics/MOM_spatial_means.F90 .
./MOM_hor_visc.F90: $(SRCROOT)/MOM6/src/parameterizations/lateral/MOM_hor_visc.F90
	cp $(SRCROOT)/MOM6/src/parameterizations/lateral/MOM_hor_visc.F90 .
./MOM_debugging.F90: $(SRCROOT)/MOM6/src/diagnostics/MOM_debugging.F90
	cp $(SRCROOT)/MOM6/src/diagnostics/MOM_debugging.F90 .
./MOM_forcing_type.F90: $(SRCROOT)/MOM6/src/core/MOM_forcing_type.F90
	cp $(SRCROOT)/MOM6/src/core/MOM_forcing_type.F90 .
./MOM_spherical_harmonics.F90: $(SRCROOT)/MOM6/src/parameterizations/lateral/MOM_spherical_harmonics.F90
	cp $(SRCROOT)/MOM6/src/parameterizations/lateral/MOM_spherical_harmonics.F90 .
./dyed_obc_tracer.F90: $(SRCROOT)/MOM6/src/tracer/dyed_obc_tracer.F90
	cp $(SRCROOT)/MOM6/src/tracer/dyed_obc_tracer.F90 .
./gsw_t_deriv_chem_potential_water_t_exact.f90: $(SRCROOT)/MOM6/src/equation_of_state/TEOS10/gsw_t_deriv_chem_potential_water_t_exact.f90
	cp $(SRCROOT)/MOM6/src/equation_of_state/TEOS10/gsw_t_deriv_chem_potential_water_t_exact.f90 .
./ocean_da_types.F90: $(SRCROOT)/MOM6/config_src/external/ODA_hooks/ocean_da_types.F90
	cp $(SRCROOT)/MOM6/config_src/external/ODA_hooks/ocean_da_types.F90 .
./gsw_pt_from_ct.f90: $(SRCROOT)/MOM6/src/equation_of_state/TEOS10/gsw_pt_from_ct.f90
	cp $(SRCROOT)/MOM6/src/equation_of_state/TEOS10/gsw_pt_from_ct.f90 .
./gsw_ct_from_pt.f90: $(SRCROOT)/MOM6/src/equation_of_state/TEOS10/gsw_ct_from_pt.f90
	cp $(SRCROOT)/MOM6/src/equation_of_state/TEOS10/gsw_ct_from_pt.f90 .
./MOM_cpu_clock_infra.F90: $(SRCROOT)/MOM6/config_src/infra/FMS2/MOM_cpu_clock_infra.F90
	cp $(SRCROOT)/MOM6/config_src/infra/FMS2/MOM_cpu_clock_infra.F90 .
./MOM_checksums.F90: $(SRCROOT)/MOM6/src/framework/MOM_checksums.F90
	cp $(SRCROOT)/MOM6/src/framework/MOM_checksums.F90 .
./MOM_kappa_shear.F90: $(SRCROOT)/MOM6/src/parameterizations/vertical/MOM_kappa_shear.F90
	cp $(SRCROOT)/MOM6/src/parameterizations/vertical/MOM_kappa_shear.F90 .
./FMS_coupler_util.F90: $(SRCROOT)/MOM6/config_src/external/GFDL_ocean_BGC/FMS_coupler_util.F90
	cp $(SRCROOT)/MOM6/config_src/external/GFDL_ocean_BGC/FMS_coupler_util.F90 .
./MOM_EOS_Roquet_SpV.F90: $(SRCROOT)/MOM6/src/equation_of_state/MOM_EOS_Roquet_SpV.F90
	cp $(SRCROOT)/MOM6/src/equation_of_state/MOM_EOS_Roquet_SpV.F90 .
./MOM_EOS_TEOS10.F90: $(SRCROOT)/MOM6/src/equation_of_state/MOM_EOS_TEOS10.F90
	cp $(SRCROOT)/MOM6/src/equation_of_state/MOM_EOS_TEOS10.F90 .
./MOM_wave_speed.F90: $(SRCROOT)/MOM6/src/diagnostics/MOM_wave_speed.F90
	cp $(SRCROOT)/MOM6/src/diagnostics/MOM_wave_speed.F90 .
./kdtree.f90: $(SRCROOT)/MOM6/config_src/external/ODA_hooks/kdtree.f90
	cp $(SRCROOT)/MOM6/config_src/external/ODA_hooks/kdtree.f90 .
./Phillips_initialization.F90: $(SRCROOT)/MOM6/src/user/Phillips_initialization.F90
	cp $(SRCROOT)/MOM6/src/user/Phillips_initialization.F90 .
./MOM_tracer_Z_init.F90: $(SRCROOT)/MOM6/src/tracer/MOM_tracer_Z_init.F90
	cp $(SRCROOT)/MOM6/src/tracer/MOM_tracer_Z_init.F90 .
./MOM_EOS_Wright_red.F90: $(SRCROOT)/MOM6/src/equation_of_state/MOM_EOS_Wright_red.F90
	cp $(SRCROOT)/MOM6/src/equation_of_state/MOM_EOS_Wright_red.F90 .
./MOM_offline_aux.F90: $(SRCROOT)/MOM6/src/tracer/MOM_offline_aux.F90
	cp $(SRCROOT)/MOM6/src/tracer/MOM_offline_aux.F90 .
./MOM_interface_filter.F90: $(SRCROOT)/MOM6/src/parameterizations/lateral/MOM_interface_filter.F90
	cp $(SRCROOT)/MOM6/src/parameterizations/lateral/MOM_interface_filter.F90 .
./MOM_interpolate.F90: $(SRCROOT)/MOM6/src/framework/MOM_interpolate.F90
	cp $(SRCROOT)/MOM6/src/framework/MOM_interpolate.F90 .
./MOM_CVMix_KPP.F90: $(SRCROOT)/MOM6/src/parameterizations/vertical/MOM_CVMix_KPP.F90
	cp $(SRCROOT)/MOM6/src/parameterizations/vertical/MOM_CVMix_KPP.F90 .
./MOM_data_override.F90: $(SRCROOT)/MOM6/src/framework/MOM_data_override.F90
	cp $(SRCROOT)/MOM6/src/framework/MOM_data_override.F90 .
./gsw_rho_first_derivatives.f90: $(SRCROOT)/MOM6/src/equation_of_state/TEOS10/gsw_rho_first_derivatives.f90
	cp $(SRCROOT)/MOM6/src/equation_of_state/TEOS10/gsw_rho_first_derivatives.f90 .
./MOM_tracer_advect.F90: $(SRCROOT)/MOM6/src/tracer/MOM_tracer_advect.F90
	cp $(SRCROOT)/MOM6/src/tracer/MOM_tracer_advect.F90 .
./MOM_time_manager.F90: $(SRCROOT)/MOM6/config_src/infra/FMS2/MOM_time_manager.F90
	cp $(SRCROOT)/MOM6/config_src/infra/FMS2/MOM_time_manager.F90 .
./MOM_MEKE_types.F90: $(SRCROOT)/MOM6/src/parameterizations/lateral/MOM_MEKE_types.F90
	cp $(SRCROOT)/MOM6/src/parameterizations/lateral/MOM_MEKE_types.F90 .
./MOM_isopycnal_slopes.F90: $(SRCROOT)/MOM6/src/core/MOM_isopycnal_slopes.F90
	cp $(SRCROOT)/MOM6/src/core/MOM_isopycnal_slopes.F90 .
./oil_tracer.F90: $(SRCROOT)/MOM6/src/tracer/oil_tracer.F90
	cp $(SRCROOT)/MOM6/src/tracer/oil_tracer.F90 .
./gsw_entropy_part_zerop.f90: $(SRCROOT)/MOM6/src/equation_of_state/TEOS10/gsw_entropy_part_zerop.f90
	cp $(SRCROOT)/MOM6/src/equation_of_state/TEOS10/gsw_entropy_part_zerop.f90 .
./PCM_functions.F90: $(SRCROOT)/MOM6/src/ALE/PCM_functions.F90
	cp $(SRCROOT)/MOM6/src/ALE/PCM_functions.F90 .
./MOM_write_cputime.F90: $(SRCROOT)/MOM6/src/framework/MOM_write_cputime.F90
	cp $(SRCROOT)/MOM6/src/framework/MOM_write_cputime.F90 .
./gsw_t_freezing_exact.f90: $(SRCROOT)/MOM6/src/equation_of_state/TEOS10/gsw_t_freezing_exact.f90
	cp $(SRCROOT)/MOM6/src/equation_of_state/TEOS10/gsw_t_freezing_exact.f90 .
./DOME_tracer.F90: $(SRCROOT)/MOM6/src/tracer/DOME_tracer.F90
	cp $(SRCROOT)/MOM6/src/tracer/DOME_tracer.F90 .
./dumbbell_surface_forcing.F90: $(SRCROOT)/MOM6/src/user/dumbbell_surface_forcing.F90
	cp $(SRCROOT)/MOM6/src/user/dumbbell_surface_forcing.F90 .
./MOM_tidal_forcing.F90: $(SRCROOT)/MOM6/src/parameterizations/lateral/MOM_tidal_forcing.F90
	cp $(SRCROOT)/MOM6/src/parameterizations/lateral/MOM_tidal_forcing.F90 .
./coord_hycom.F90: $(SRCROOT)/MOM6/src/ALE/coord_hycom.F90
	cp $(SRCROOT)/MOM6/src/ALE/coord_hycom.F90 .
./MOM_generic_tracer.F90: $(SRCROOT)/MOM6/src/tracer/MOM_generic_tracer.F90
	cp $(SRCROOT)/MOM6/src/tracer/MOM_generic_tracer.F90 .
./MOM_file_parser_tests.F90: $(SRCROOT)/MOM6/src/framework/testing/MOM_file_parser_tests.F90
	cp $(SRCROOT)/MOM6/src/framework/testing/MOM_file_parser_tests.F90 .
./polynomial_functions.F90: $(SRCROOT)/MOM6/src/ALE/polynomial_functions.F90
	cp $(SRCROOT)/MOM6/src/ALE/polynomial_functions.F90 .
./MOM_unique_scales.F90: $(SRCROOT)/MOM6/src/framework/MOM_unique_scales.F90
	cp $(SRCROOT)/MOM6/src/framework/MOM_unique_scales.F90 .
./MOM_boundary_update.F90: $(SRCROOT)/MOM6/src/core/MOM_boundary_update.F90
	cp $(SRCROOT)/MOM6/src/core/MOM_boundary_update.F90 .
./MOM_lateral_boundary_diffusion.F90: $(SRCROOT)/MOM6/src/tracer/MOM_lateral_boundary_diffusion.F90
	cp $(SRCROOT)/MOM6/src/tracer/MOM_lateral_boundary_diffusion.F90 .
./MOM_check_scaling.F90: $(SRCROOT)/MOM6/src/core/MOM_check_scaling.F90
	cp $(SRCROOT)/MOM6/src/core/MOM_check_scaling.F90 .
./MOM_neutral_diffusion.F90: $(SRCROOT)/MOM6/src/tracer/MOM_neutral_diffusion.F90
	cp $(SRCROOT)/MOM6/src/tracer/MOM_neutral_diffusion.F90 .
./ocean_model_MOM.F90: $(SRCROOT)/MOM6/config_src/drivers/FMS_cap/ocean_model_MOM.F90
	cp $(SRCROOT)/MOM6/config_src/drivers/FMS_cap/ocean_model_MOM.F90 .
./nw2_tracers.F90: $(SRCROOT)/MOM6/src/tracer/nw2_tracers.F90
	cp $(SRCROOT)/MOM6/src/tracer/nw2_tracers.F90 .
./MOM_opacity.F90: $(SRCROOT)/MOM6/src/parameterizations/vertical/MOM_opacity.F90
	cp $(SRCROOT)/MOM6/src/parameterizations/vertical/MOM_opacity.F90 .
./MOM_restart.F90: $(SRCROOT)/MOM6/src/framework/MOM_restart.F90
	cp $(SRCROOT)/MOM6/src/framework/MOM_restart.F90 .
./MOM_domains.F90: $(SRCROOT)/MOM6/src/framework/MOM_domains.F90
	cp $(SRCROOT)/MOM6/src/framework/MOM_domains.F90 .
./MOM_stoch_eos.F90: $(SRCROOT)/MOM6/src/core/MOM_stoch_eos.F90
	cp $(SRCROOT)/MOM6/src/core/MOM_stoch_eos.F90 .
./posix.h: $(SRCROOT)/MOM6/src/framework/posix.h
	cp $(SRCROOT)/MOM6/src/framework/posix.h .
./MOM_horizontal_regridding.F90: $(SRCROOT)/MOM6/src/framework/MOM_horizontal_regridding.F90
	cp $(SRCROOT)/MOM6/src/framework/MOM_horizontal_regridding.F90 .
./MOM_ALE_sponge.F90: $(SRCROOT)/MOM6/src/parameterizations/vertical/MOM_ALE_sponge.F90
	cp $(SRCROOT)/MOM6/src/parameterizations/vertical/MOM_ALE_sponge.F90 .
./ocean_da_core.F90: $(SRCROOT)/MOM6/config_src/external/ODA_hooks/ocean_da_core.F90
	cp $(SRCROOT)/MOM6/config_src/external/ODA_hooks/ocean_da_core.F90 .
./MOM_bulk_mixed_layer.F90: $(SRCROOT)/MOM6/src/parameterizations/vertical/MOM_bulk_mixed_layer.F90
	cp $(SRCROOT)/MOM6/src/parameterizations/vertical/MOM_bulk_mixed_layer.F90 .
./regrid_interp.F90: $(SRCROOT)/MOM6/src/ALE/regrid_interp.F90
	cp $(SRCROOT)/MOM6/src/ALE/regrid_interp.F90 .
./MOM_set_viscosity.F90: $(SRCROOT)/MOM6/src/parameterizations/vertical/MOM_set_viscosity.F90
	cp $(SRCROOT)/MOM6/src/parameterizations/vertical/MOM_set_viscosity.F90 .
./write_ocean_obs.F90: $(SRCROOT)/MOM6/config_src/external/ODA_hooks/write_ocean_obs.F90
	cp $(SRCROOT)/MOM6/config_src/external/ODA_hooks/write_ocean_obs.F90 .
./dye_example.F90: $(SRCROOT)/MOM6/src/tracer/dye_example.F90
	cp $(SRCROOT)/MOM6/src/tracer/dye_example.F90 .
./coord_adapt.F90: $(SRCROOT)/MOM6/src/ALE/coord_adapt.F90
	cp $(SRCROOT)/MOM6/src/ALE/coord_adapt.F90 .
./MOM_mixed_layer_restrat.F90: $(SRCROOT)/MOM6/src/parameterizations/lateral/MOM_mixed_layer_restrat.F90
	cp $(SRCROOT)/MOM6/src/parameterizations/lateral/MOM_mixed_layer_restrat.F90 .
./SCM_CVMix_tests.F90: $(SRCROOT)/MOM6/src/user/SCM_CVMix_tests.F90
	cp $(SRCROOT)/MOM6/src/user/SCM_CVMix_tests.F90 .
./gsw_specvol_second_derivatives.f90: $(SRCROOT)/MOM6/src/equation_of_state/TEOS10/gsw_specvol_second_derivatives.f90
	cp $(SRCROOT)/MOM6/src/equation_of_state/TEOS10/gsw_specvol_second_derivatives.f90 .
./MOM_vert_friction.F90: $(SRCROOT)/MOM6/src/parameterizations/vertical/MOM_vert_friction.F90
	cp $(SRCROOT)/MOM6/src/parameterizations/vertical/MOM_vert_friction.F90 .
./generic_tracer.F90: $(SRCROOT)/MOM6/config_src/external/GFDL_ocean_BGC/generic_tracer.F90
	cp $(SRCROOT)/MOM6/config_src/external/GFDL_ocean_BGC/generic_tracer.F90 .
./Neverworld_initialization.F90: $(SRCROOT)/MOM6/src/user/Neverworld_initialization.F90
	cp $(SRCROOT)/MOM6/src/user/Neverworld_initialization.F90 .
./RGC_initialization.F90: $(SRCROOT)/MOM6/src/user/RGC_initialization.F90
	cp $(SRCROOT)/MOM6/src/user/RGC_initialization.F90 .
./MOM_PointAccel.F90: $(SRCROOT)/MOM6/src/diagnostics/MOM_PointAccel.F90
	cp $(SRCROOT)/MOM6/src/diagnostics/MOM_PointAccel.F90 .
./MOM_OCMIP2_CFC.F90: $(SRCROOT)/MOM6/src/tracer/MOM_OCMIP2_CFC.F90
	cp $(SRCROOT)/MOM6/src/tracer/MOM_OCMIP2_CFC.F90 .
./dense_water_initialization.F90: $(SRCROOT)/MOM6/src/user/dense_water_initialization.F90
	cp $(SRCROOT)/MOM6/src/user/dense_water_initialization.F90 .
./MOM_diabatic_driver.F90: $(SRCROOT)/MOM6/src/parameterizations/vertical/MOM_diabatic_driver.F90
	cp $(SRCROOT)/MOM6/src/parameterizations/vertical/MOM_diabatic_driver.F90 .
./adjustment_initialization.F90: $(SRCROOT)/MOM6/src/user/adjustment_initialization.F90
	cp $(SRCROOT)/MOM6/src/user/adjustment_initialization.F90 .
./cvmix_background.F90: $(SRCROOT)/MOM6/src/parameterizations/CVmix/cvmix_background.F90
	cp $(SRCROOT)/MOM6/src/parameterizations/CVmix/cvmix_background.F90 .
./MOM_tracer_hor_diff.F90: $(SRCROOT)/MOM6/src/tracer/MOM_tracer_hor_diff.F90
	cp $(SRCROOT)/MOM6/src/tracer/MOM_tracer_hor_diff.F90 .
SRC = $(SRCROOT)/MOM6/config_src/external/database_comms/database_client_interface.F90 $(SRCROOT)/MOM6/src/ALE/PCM_functions.F90 $(SRCROOT)/MOM6/src/equation_of_state/TEOS10/gsw_gibbs_pt0_pt0.f90 $(SRCROOT)/MOM6/src/framework/MOM_safe_alloc.F90 $(SRCROOT)/MOM6/src/equation_of_state/TEOS10/gsw_gibbs.f90 $(SRCROOT)/MOM6/src/initialization/MOM_state_initialization.F90 $(SRCROOT)/MOM6/src/user/dumbbell_surface_forcing.F90 $(SRCROOT)/MOM6/src/user/basin_builder.F90 $(SRCROOT)/MOM6/src/parameterizations/vertical/MOM_diabatic_aux.F90 $(SRCROOT)/MOM6/src/equation_of_state/TEOS10/gsw_pt_from_t.f90 $(SRCROOT)/MOM6/src/parameterizations/lateral/MOM_tidal_forcing.F90 $(SRCROOT)/MOM6/src/equation_of_state/MOM_EOS_UNESCO.F90 $(SRCROOT)/MOM6/src/user/supercritical_initialization.F90 $(SRCROOT)/MOM6/src/user/benchmark_initialization.F90 $(SRCROOT)/MOM6/config_src/infra/FMS2/MOM_io_infra.F90 $(SRCROOT)/MOM6/src/core/MOM_open_boundary.F90 $(SRCROOT)/MOM6/src/ALE/remapping_attic.F90 $(SRCROOT)/MOM6/src/core/MOM_boundary_update.F90 $(SRCROOT)/MOM6/src/parameterizations/vertical/MOM_bkgnd_mixing.F90 $(SRCROOT)/MOM6/src/equation_of_state/TEOS10/gsw_pt_from_ct.f90 $(SRCROOT)/MOM6/src/parameterizations/vertical/MOM_energetic_PBL.F90 $(SRCROOT)/MOM6/src/framework/MOM_diag_remap.F90 $(SRCROOT)/MOM6/src/core/MOM_barotropic.F90 $(SRCROOT)/MOM6/src/initialization/MOM_tracer_initialization_from_Z.F90 $(SRCROOT)/MOM6/src/framework/MOM_string_functions.F90 $(SRCROOT)/MOM6/src/tracer/RGC_tracer.F90 $(SRCROOT)/MOM6/src/tracer/tracer_example.F90 $(SRCROOT)/MOM6/src/ice_shelf/MOM_ice_shelf_diag_mediator.F90 $(SRCROOT)/MOM6/src/core/MOM_transcribe_grid.F90 $(SRCROOT)/MOM6/src/core/MOM_dynamics_split_RK2.F90 $(SRCROOT)/MOM6/src/framework/MOM_coupler_types.F90 $(SRCROOT)/MOM6/src/user/BFB_initialization.F90 $(SRCROOT)/MOM6/src/core/MOM.F90 $(SRCROOT)/MOM6/src/framework/MOM_hor_index.F90 $(SRCROOT)/MOM6/src/ALE/coord_zlike.F90 $(SRCROOT)/MOM6/src/user/adjustment_initialization.F90 $(SRCROOT)/MOM6/config_src/drivers/FMS_cap/MOM_surface_forcing_gfdl.F90 $(SRCROOT)/MOM6/src/parameterizations/lateral/MOM_MEKE.F90 $(SRCROOT)/MOM6/config_src/drivers/FMS_cap/ocean_model_MOM.F90 $(SRCROOT)/MOM6/src/parameterizations/vertical/MOM_CVMix_shear.F90 $(SRCROOT)/MOM6/src/user/BFB_surface_forcing.F90 $(SRCROOT)/MOM6/src/user/circle_obcs_initialization.F90 $(SRCROOT)/MOM6/src/tracer/MOM_tracer_diabatic.F90 $(SRCROOT)/MOM6/src/user/SCM_CVMix_tests.F90 $(SRCROOT)/MOM6/config_src/infra/FMS2/MOM_cpu_clock_infra.F90 $(SRCROOT)/MOM6/src/parameterizations/vertical/MOM_full_convection.F90 $(SRCROOT)/MOM6/src/framework/testing/MOM_file_parser_tests.F90 $(SRCROOT)/MOM6/config_src/infra/FMS2/MOM_interp_infra.F90 $(SRCROOT)/MOM6/src/user/soliton_initialization.F90 $(SRCROOT)/MOM6/src/core/MOM_stoch_eos.F90 $(SRCROOT)/MOM6/src/ice_shelf/MOM_ice_shelf_state.F90 $(SRCROOT)/MOM6/src/tracer/nw2_tracers.F90 $(SRCROOT)/MOM6/src/parameterizations/lateral/MOM_lateral_mixing_coeffs.F90 $(SRCROOT)/MOM6/src/tracer/dye_example.F90 $(SRCROOT)/MOM6/src/parameterizations/vertical/MOM_CVMix_conv.F90 $(SRCROOT)/MOM6/src/parameterizations/vertical/MOM_CVMix_KPP.F90 $(SRCROOT)/MOM6/src/tracer/MOM_tracer_advect.F90 $(SRCROOT)/MOM6/src/user/dyed_obcs_initialization.F90 $(SRCROOT)/MOM6/src/user/Rossby_front_2d_initialization.F90 $(SRCROOT)/MOM6/src/parameterizations/CVmix/cvmix_convection.F90 $(SRCROOT)/MOM6/src/parameterizations/CVmix/cvmix_background.F90 $(SRCROOT)/MOM6/src/parameterizations/CVmix/cvmix_kinds_and_types.F90 $(SRCROOT)/MOM6/src/equation_of_state/MOM_TFreeze.F90 $(SRCROOT)/MOM6/src/equation_of_state/TEOS10/gsw_sr_from_sp.f90 $(SRCROOT)/MOM6/config_src/external/stochastic_physics/get_stochy_pattern.F90 $(SRCROOT)/MOM6/src/framework/MOM_data_override.F90 $(SRCROOT)/MOM6/src/user/seamount_initialization.F90 $(SRCROOT)/MOM6/src/parameterizations/vertical/MOM_regularize_layers.F90 $(SRCROOT)/MOM6/src/core/MOM_variables.F90 $(SRCROOT)/MOM6/src/framework/MOM_interpolate.F90 $(SRCROOT)/MOM6/src/framework/posix.F90 $(SRCROOT)/MOM6/src/framework/MOM_write_cputime.F90 $(SRCROOT)/MOM6/src/diagnostics/MOM_wave_speed.F90 $(SRCROOT)/MOM6/src/core/MOM_density_integrals.F90 $(SRCROOT)/MOM6/src/equation_of_state/MOM_EOS_Roquet_rho.F90 $(SRCROOT)/MOM6/src/tracer/pseudo_salt_tracer.F90 $(SRCROOT)/MOM6/src/user/external_gwave_initialization.F90 $(SRCROOT)/MOM6/src/diagnostics/MOM_obsolete_diagnostics.F90 $(SRCROOT)/MOM6/src/tracer/MOM_CFC_cap.F90 $(SRCROOT)/MOM6/src/parameterizations/lateral/MOM_MEKE_types.F90 $(SRCROOT)/MOM6/src/equation_of_state/TEOS10/gsw_pt0_from_t.f90 $(SRCROOT)/MOM6/src/framework/MOM_random.F90 $(SRCROOT)/MOM6/src/user/ISOMIP_initialization.F90 $(SRCROOT)/MOM6/src/framework/MOM_io_file.F90 $(SRCROOT)/MOM6/src/user/MOM_wave_interface.F90 $(SRCROOT)/MOM6/src/equation_of_state/MOM_EOS_Roquet_SpV.F90 $(SRCROOT)/MOM6/src/parameterizations/vertical/MOM_entrain_diffusive.F90 $(SRCROOT)/MOM6/src/parameterizations/vertical/MOM_sponge.F90 $(SRCROOT)/MOM6/src/tracer/oil_tracer.F90 $(SRCROOT)/MOM6/config_src/infra/FMS2/MOM_domain_infra.F90 $(SRCROOT)/MOM6/src/framework/MOM_file_parser.F90 $(SRCROOT)/MOM6/src/equation_of_state/TEOS10/gsw_t_freezing_poly.f90 $(SRCROOT)/MOM6/src/equation_of_state/TEOS10/gsw_rho_second_derivatives.f90 $(SRCROOT)/MOM6/src/ALE/regrid_edge_values.F90 $(SRCROOT)/MOM6/src/parameterizations/vertical/MOM_set_viscosity.F90 $(SRCROOT)/MOM6/src/user/baroclinic_zone_initialization.F90 $(SRCROOT)/MOM6/src/equation_of_state/TEOS10/gsw_entropy_part.f90 $(SRCROOT)/MOM6/src/ALE/P3M_functions.F90 $(SRCROOT)/MOM6/src/equation_of_state/TEOS10/gsw_ct_from_t.f90 $(SRCROOT)/MOM6/src/parameterizations/CVmix/cvmix_kpp.F90 $(SRCROOT)/MOM6/config_src/infra/FMS2/MOM_diag_manager_infra.F90 $(SRCROOT)/MOM6/src/user/user_change_diffusivity.F90 $(SRCROOT)/MOM6/src/tracer/MOM_tracer_flow_control.F90 $(SRCROOT)/MOM6/config_src/external/drifters/MOM_particles_types.F90 $(SRCROOT)/MOM6/src/tracer/MOM_neutral_diffusion.F90 $(SRCROOT)/MOM6/src/equation_of_state/TEOS10/gsw_specvol_second_derivatives.f90 $(SRCROOT)/MOM6/src/parameterizations/vertical/MOM_opacity.F90 $(SRCROOT)/MOM6/src/parameterizations/vertical/MOM_geothermal.F90 $(SRCROOT)/MOM6/src/framework/MOM_document.F90 $(SRCROOT)/MOM6/src/ice_shelf/user_shelf_init.F90 $(SRCROOT)/MOM6/src/tracer/MOM_OCMIP2_CFC.F90 $(SRCROOT)/MOM6/src/core/MOM_isopycnal_slopes.F90 $(SRCROOT)/MOM6/src/ALE/PQM_functions.F90 $(SRCROOT)/MOM6/src/framework/MOM_array_transform.F90 $(SRCROOT)/MOM6/src/core/MOM_forcing_type.F90 $(SRCROOT)/MOM6/src/ice_shelf/MOM_marine_ice.F90 $(SRCROOT)/MOM6/src/initialization/MOM_fixed_initialization.F90 $(SRCROOT)/MOM6/src/user/Kelvin_initialization.F90 $(SRCROOT)/MOM6/src/diagnostics/MOM_obsolete_params.F90 $(SRCROOT)/MOM6/config_src/infra/FMS2/MOM_couplertype_infra.F90 $(SRCROOT)/MOM6/src/equation_of_state/TEOS10/gsw_rho.f90 $(SRCROOT)/MOM6/src/user/DOME_initialization.F90 $(SRCROOT)/MOM6/src/equation_of_state/MOM_EOS_Jackett06.F90 $(SRCROOT)/MOM6/src/parameterizations/vertical/MOM_ALE_sponge.F90 $(SRCROOT)/MOM6/src/tracer/advection_test_tracer.F90 $(SRCROOT)/MOM6/src/core/MOM_interface_heights.F90 $(SRCROOT)/MOM6/src/equation_of_state/TEOS10/gsw_mod_gibbs_ice_coefficients.f90 $(SRCROOT)/MOM6/src/equation_of_state/TEOS10/gsw_entropy_part_zerop.f90 $(SRCROOT)/MOM6/src/ice_shelf/MOM_ice_shelf_dynamics.F90 $(SRCROOT)/MOM6/src/diagnostics/MOM_spatial_means.F90 $(SRCROOT)/MOM6/src/framework/MOM_dyn_horgrid.F90 $(SRCROOT)/MOM6/src/framework/MOM_unit_testing.F90 $(SRCROOT)/MOM6/src/core/MOM_verticalGrid.F90 $(SRCROOT)/MOM6/config_src/infra/FMS2/MOM_error_infra.F90 $(SRCROOT)/MOM6/src/parameterizations/CVmix/cvmix_math.F90 $(SRCROOT)/MOM6/src/framework/MOM_unique_scales.F90 $(SRCROOT)/MOM6/src/tracer/ideal_age_example.F90 $(SRCROOT)/MOM6/src/initialization/MOM_shared_initialization.F90 $(SRCROOT)/MOM6/config_src/infra/FMS2/MOM_coms_infra.F90 $(SRCROOT)/MOM6/src/core/MOM_unit_tests.F90 $(SRCROOT)/MOM6/src/equation_of_state/TEOS10/gsw_t_deriv_chem_potential_water_t_exact.f90 $(SRCROOT)/MOM6/src/framework/MOM_diag_mediator.F90 $(SRCROOT)/MOM6/config_src/external/GFDL_ocean_BGC/FMS_coupler_util.F90 $(SRCROOT)/MOM6/src/ocean_data_assim/MOM_oda_incupd.F90 $(SRCROOT)/MOM6/src/ALE/MOM_hybgen_unmix.F90 $(SRCROOT)/MOM6/src/parameterizations/CVmix/cvmix_tidal.F90 $(SRCROOT)/MOM6/src/tracer/MOM_tracer_hor_diff.F90 $(SRCROOT)/MOM6/src/parameterizations/lateral/MOM_thickness_diffuse.F90 $(SRCROOT)/MOM6/config_src/external/stochastic_physics/stochastic_physics.F90 $(SRCROOT)/MOM6/src/core/MOM_check_scaling.F90 $(SRCROOT)/MOM6/src/equation_of_state/MOM_EOS_TEOS10.F90 $(SRCROOT)/MOM6/src/equation_of_state/TEOS10/gsw_mod_specvol_coefficients.f90 $(SRCROOT)/MOM6/src/parameterizations/vertical/MOM_internal_tide_input.F90 $(SRCROOT)/MOM6/src/core/MOM_PressureForce_FV.F90 $(SRCROOT)/MOM6/src/ALE/coord_rho.F90 $(SRCROOT)/MOM6/src/parameterizations/CVmix/cvmix_shear.F90 $(SRCROOT)/MOM6/src/equation_of_state/MOM_EOS.F90 $(SRCROOT)/MOM6/src/parameterizations/vertical/MOM_set_diffusivity.F90 $(SRCROOT)/MOM6/src/core/MOM_checksum_packages.F90 $(SRCROOT)/MOM6/src/parameterizations/lateral/MOM_load_love_numbers.F90 $(SRCROOT)/MOM6/src/parameterizations/CVmix/cvmix_utils.F90 $(SRCROOT)/MOM6/src/equation_of_state/TEOS10/gsw_specvol_first_derivatives.f90 $(SRCROOT)/MOM6/src/equation_of_state/TEOS10/gsw_mod_freezing_poly_coefficients.f90 $(SRCROOT)/MOM6/config_src/infra/FMS2/MOM_ensemble_manager_infra.F90 $(SRCROOT)/MOM6/src/equation_of_state/MOM_EOS_Wright_red.F90 $(SRCROOT)/MOM6/src/tracer/ISOMIP_tracer.F90 $(SRCROOT)/MOM6/src/equation_of_state/TEOS10/gsw_t_freezing_exact.f90 $(SRCROOT)/MOM6/src/equation_of_state/TEOS10/gsw_chem_potential_water_t_exact.f90 $(SRCROOT)/MOM6/src/parameterizations/vertical/MOM_CVMix_ddiff.F90 $(SRCROOT)/MOM6/src/framework/MOM_io.F90 $(SRCROOT)/MOM6/src/parameterizations/stochastic/MOM_stochastics.F90 $(SRCROOT)/MOM6/src/user/dumbbell_initialization.F90 $(SRCROOT)/MOM6/src/tracer/MOM_lateral_boundary_diffusion.F90 $(SRCROOT)/MOM6/src/framework/MOM_get_input.F90 $(SRCROOT)/MOM6/src/tracer/MOM_generic_tracer.F90 $(SRCROOT)/MOM6/src/parameterizations/lateral/MOM_interface_filter.F90 $(SRCROOT)/MOM6/src/framework/MOM_horizontal_regridding.F90 $(SRCROOT)/MOM6/src/user/user_revise_forcing.F90 $(SRCROOT)/MOM6/src/core/MOM_dynamics_unsplit_RK2.F90 $(SRCROOT)/MOM6/config_src/infra/FMS2/MOM_time_manager.F90 $(SRCROOT)/MOM6/src/diagnostics/MOM_sum_output.F90 $(SRCROOT)/MOM6/src/tracer/MOM_offline_aux.F90 $(SRCROOT)/MOM6/config_src/external/ODA_hooks/write_ocean_obs.F90 $(SRCROOT)/MOM6/src/equation_of_state/MOM_temperature_convert.F90 $(SRCROOT)/MOM6/src/core/MOM_continuity.F90 $(SRCROOT)/MOM6/config_src/external/GFDL_ocean_BGC/generic_tracer.F90 $(SRCROOT)/MOM6/src/user/Idealized_Hurricane.F90 $(SRCROOT)/MOM6/src/framework/MOM_checksums.F90 $(SRCROOT)/MOM6/src/equation_of_state/TEOS10/gsw_specvol.f90 $(SRCROOT)/MOM6/src/user/shelfwave_initialization.F90 $(SRCROOT)/MOM6/src/parameterizations/vertical/MOM_diabatic_driver.F90 $(SRCROOT)/MOM6/src/equation_of_state/TEOS10/gsw_mod_teos10_constants.f90 $(SRCROOT)/MOM6/src/initialization/MOM_grid_initialize.F90 $(SRCROOT)/MOM6/src/equation_of_state/TEOS10/gsw_mod_error_functions.f90 $(SRCROOT)/MOM6/src/parameterizations/vertical/MOM_tidal_mixing.F90 $(SRCROOT)/MOM6/config_src/external/database_comms/MOM_database_comms.F90 $(SRCROOT)/MOM6/src/user/dense_water_initialization.F90 $(SRCROOT)/MOM6/src/user/MOM_controlled_forcing.F90 $(SRCROOT)/MOM6/src/parameterizations/CVmix/cvmix_ddiff.F90 $(SRCROOT)/MOM6/src/tracer/MOM_tracer_types.F90 $(SRCROOT)/MOM6/config_src/external/ODA_hooks/ocean_da_types.F90 $(SRCROOT)/MOM6/src/ALE/coord_sigma.F90 $(SRCROOT)/MOM6/src/framework/MOM_cpu_clock.F90 $(SRCROOT)/MOM6/src/diagnostics/MOM_PointAccel.F90 $(SRCROOT)/MOM6/src/ALE/PLM_functions.F90 $(SRCROOT)/MOM6/src/diagnostics/MOM_debugging.F90 $(SRCROOT)/MOM6/src/ALE/coord_adapt.F90 $(SRCROOT)/MOM6/src/parameterizations/vertical/MOM_bulk_mixed_layer.F90 $(SRCROOT)/MOM6/src/user/tidal_bay_initialization.F90 $(SRCROOT)/MOM6/src/equation_of_state/MOM_EOS_Wright_full.F90 $(SRCROOT)/MOM6/src/framework/MOM_intrinsic_functions.F90 $(SRCROOT)/MOM6/src/equation_of_state/TEOS10/gsw_sp_from_sr.f90 $(SRCROOT)/MOM6/src/equation_of_state/TEOS10/gsw_rho_first_derivatives.f90 $(SRCROOT)/MOM6/src/ALE/MOM_regridding.F90 $(SRCROOT)/MOM6/src/tracer/dyed_obc_tracer.F90 $(SRCROOT)/MOM6/src/equation_of_state/TEOS10/gsw_gibbs_ice.f90 $(SRCROOT)/MOM6/src/user/Neverworld_initialization.F90 $(SRCROOT)/MOM6/src/core/MOM_continuity_PPM.F90 $(SRCROOT)/MOM6/src/equation_of_state/MOM_EOS_Wright.F90 $(SRCROOT)/MOM6/src/equation_of_state/MOM_EOS_linear.F90 $(SRCROOT)/MOM6/src/framework/MOM_error_handler.F90 $(SRCROOT)/MOM6/src/tracer/MOM_offline_main.F90 $(SRCROOT)/MOM6/src/equation_of_state/TEOS10/gsw_t_from_ct.f90 $(SRCROOT)/MOM6/src/equation_of_state/TEOS10/gsw_ct_freezing_exact.f90 $(SRCROOT)/MOM6/src/ALE/regrid_consts.F90 $(SRCROOT)/MOM6/config_src/external/GFDL_ocean_BGC/generic_tracer_utils.F90 $(SRCROOT)/MOM6/config_src/external/ODA_hooks/ocean_da_core.F90 $(SRCROOT)/MOM6/src/core/MOM_CoriolisAdv.F90 $(SRCROOT)/MOM6/src/parameterizations/lateral/MOM_Zanna_Bolton.F90 $(SRCROOT)/MOM6/src/framework/MOM_restart.F90 $(SRCROOT)/MOM6/src/framework/MOM_unit_scaling.F90 $(SRCROOT)/MOM6/src/core/MOM_dynamics_unsplit.F90 $(SRCROOT)/MOM6/src/user/sloshing_initialization.F90 $(SRCROOT)/MOM6/src/parameterizations/lateral/MOM_spherical_harmonics.F90 $(SRCROOT)/MOM6/src/diagnostics/MOM_diagnostics.F90 $(SRCROOT)/MOM6/src/ALE/MOM_remapping.F90 $(SRCROOT)/MOM6/src/initialization/MOM_coord_initialization.F90 $(SRCROOT)/MOM6/src/parameterizations/lateral/MOM_mixed_layer_restrat.F90 $(SRCROOT)/MOM6/src/ALE/regrid_solvers.F90 $(SRCROOT)/MOM6/src/equation_of_state/TEOS10/gsw_ct_freezing_poly.f90 $(SRCROOT)/MOM6/src/equation_of_state/TEOS10/gsw_mod_toolbox.f90 $(SRCROOT)/MOM6/src/parameterizations/lateral/MOM_hor_visc.F90 $(SRCROOT)/MOM6/src/user/lock_exchange_initialization.F90 $(SRCROOT)/MOM6/src/core/MOM_PressureForce.F90 $(SRCROOT)/MOM6/config_src/external/drifters/MOM_particles.F90 $(SRCROOT)/MOM6/src/tracer/DOME_tracer.F90 $(SRCROOT)/MOM6/src/parameterizations/vertical/MOM_kappa_shear.F90 $(SRCROOT)/MOM6/src/tracer/MOM_tracer_Z_init.F90 $(SRCROOT)/MOM6/src/framework/MOM_coms.F90 $(SRCROOT)/MOM6/src/core/MOM_PressureForce_Montgomery.F90 $(SRCROOT)/MOM6/src/parameterizations/vertical/MOM_vert_friction.F90 $(SRCROOT)/MOM6/src/core/MOM_porous_barriers.F90 $(SRCROOT)/MOM6/src/ALE/MOM_hybgen_remap.F90 $(SRCROOT)/MOM6/src/ice_shelf/MOM_ice_shelf_initialize.F90 $(SRCROOT)/MOM6/src/ALE/P1M_functions.F90 $(SRCROOT)/MOM6/src/ALE/polynomial_functions.F90 $(SRCROOT)/MOM6/src/ice_shelf/MOM_ice_shelf.F90 $(SRCROOT)/MOM6/src/framework/MOM_ensemble_manager.F90 $(SRCROOT)/MOM6/src/user/Phillips_initialization.F90 $(SRCROOT)/MOM6/src/parameterizations/vertical/MOM_diapyc_energy_req.F90 $(SRCROOT)/MOM6/src/core/MOM_grid.F90 $(SRCROOT)/MOM6/src/equation_of_state/TEOS10/gsw_ct_from_pt.f90 $(SRCROOT)/MOM6/src/parameterizations/lateral/MOM_internal_tides.F90 $(SRCROOT)/MOM6/src/user/RGC_initialization.F90 $(SRCROOT)/MOM6/src/ALE/coord_hycom.F90 $(SRCROOT)/MOM6/src/parameterizations/CVmix/cvmix_put_get.F90 $(SRCROOT)/MOM6/config_src/infra/FMS2/MOM_constants.F90 $(SRCROOT)/MOM6/src/framework/MOM_domains.F90 $(SRCROOT)/MOM6/src/ALE/MOM_hybgen_regrid.F90 $(SRCROOT)/MOM6/src/equation_of_state/TEOS10/gsw_mod_kinds.f90 $(SRCROOT)/MOM6/src/ocean_data_assim/MOM_oda_driver.F90 $(SRCROOT)/MOM6/src/tracer/boundary_impulse_tracer.F90 $(SRCROOT)/MOM6/config_src/infra/FMS2/MOM_data_override_infra.F90 $(SRCROOT)/MOM6/src/user/DOME2d_initialization.F90 $(SRCROOT)/MOM6/src/ALE/regrid_interp.F90 $(SRCROOT)/MOM6/config_src/external/ODA_hooks/kdtree.f90 $(SRCROOT)/MOM6/src/ALE/MOM_ALE.F90 $(SRCROOT)/MOM6/src/ALE/PPM_functions.F90 $(SRCROOT)/MOM6/src/tracer/MOM_tracer_registry.F90 $(SRCROOT)/MOM6/src/framework/MOM_netcdf.F90 $(SRCROOT)/MOM6/src/user/dyed_channel_initialization.F90 $(SRCROOT)/MOM6/src/user/user_initialization.F90 $(SRCROOT)/MOM6/config_src/memory/dynamic_nonsymmetric/MOM_memory.h $(SRCROOT)/MOM6/src/framework/posix.h $(SRCROOT)/MOM6/src/framework/version_variable.h $(SRCROOT)MOM6/src/framework/version_variable.h $(SRCROOT)/MOM6/src/framework/MOM_memory_macros.h
OBJ = database_client_interface.o PCM_functions.o gsw_gibbs_pt0_pt0.o MOM_safe_alloc.o gsw_gibbs.o MOM_state_initialization.o dumbbell_surface_forcing.o basin_builder.o MOM_diabatic_aux.o gsw_pt_from_t.o MOM_tidal_forcing.o MOM_EOS_UNESCO.o supercritical_initialization.o benchmark_initialization.o MOM_io_infra.o MOM_open_boundary.o remapping_attic.o MOM_boundary_update.o MOM_bkgnd_mixing.o gsw_pt_from_ct.o MOM_energetic_PBL.o MOM_diag_remap.o MOM_barotropic.o MOM_tracer_initialization_from_Z.o MOM_string_functions.o RGC_tracer.o tracer_example.o MOM_ice_shelf_diag_mediator.o MOM_transcribe_grid.o MOM_dynamics_split_RK2.o MOM_coupler_types.o BFB_initialization.o MOM.o MOM_hor_index.o coord_zlike.o adjustment_initialization.o MOM_surface_forcing_gfdl.o MOM_MEKE.o ocean_model_MOM.o MOM_CVMix_shear.o BFB_surface_forcing.o circle_obcs_initialization.o MOM_tracer_diabatic.o SCM_CVMix_tests.o MOM_cpu_clock_infra.o MOM_full_convection.o MOM_file_parser_tests.o MOM_interp_infra.o soliton_initialization.o MOM_stoch_eos.o MOM_ice_shelf_state.o nw2_tracers.o MOM_lateral_mixing_coeffs.o dye_example.o MOM_CVMix_conv.o MOM_CVMix_KPP.o MOM_tracer_advect.o dyed_obcs_initialization.o Rossby_front_2d_initialization.o cvmix_convection.o cvmix_background.o cvmix_kinds_and_types.o MOM_TFreeze.o gsw_sr_from_sp.o get_stochy_pattern.o MOM_data_override.o seamount_initialization.o MOM_regularize_layers.o MOM_variables.o MOM_interpolate.o posix.o MOM_write_cputime.o MOM_wave_speed.o MOM_density_integrals.o MOM_EOS_Roquet_rho.o pseudo_salt_tracer.o external_gwave_initialization.o MOM_obsolete_diagnostics.o MOM_CFC_cap.o MOM_MEKE_types.o gsw_pt0_from_t.o MOM_random.o ISOMIP_initialization.o MOM_io_file.o MOM_wave_interface.o MOM_EOS_Roquet_SpV.o MOM_entrain_diffusive.o MOM_sponge.o oil_tracer.o MOM_domain_infra.o MOM_file_parser.o gsw_t_freezing_poly.o gsw_rho_second_derivatives.o regrid_edge_values.o MOM_set_viscosity.o baroclinic_zone_initialization.o gsw_entropy_part.o P3M_functions.o gsw_ct_from_t.o cvmix_kpp.o MOM_diag_manager_infra.o user_change_diffusivity.o MOM_tracer_flow_control.o MOM_particles_types.o MOM_neutral_diffusion.o gsw_specvol_second_derivatives.o MOM_opacity.o MOM_geothermal.o MOM_document.o user_shelf_init.o MOM_OCMIP2_CFC.o MOM_isopycnal_slopes.o PQM_functions.o MOM_array_transform.o MOM_forcing_type.o MOM_marine_ice.o MOM_fixed_initialization.o Kelvin_initialization.o MOM_obsolete_params.o MOM_couplertype_infra.o gsw_rho.o DOME_initialization.o MOM_EOS_Jackett06.o MOM_ALE_sponge.o advection_test_tracer.o MOM_interface_heights.o gsw_mod_gibbs_ice_coefficients.o gsw_entropy_part_zerop.o MOM_ice_shelf_dynamics.o MOM_spatial_means.o MOM_dyn_horgrid.o MOM_unit_testing.o MOM_verticalGrid.o MOM_error_infra.o cvmix_math.o MOM_unique_scales.o ideal_age_example.o MOM_shared_initialization.o MOM_coms_infra.o MOM_unit_tests.o gsw_t_deriv_chem_potential_water_t_exact.o MOM_diag_mediator.o FMS_coupler_util.o MOM_oda_incupd.o MOM_hybgen_unmix.o cvmix_tidal.o MOM_tracer_hor_diff.o MOM_thickness_diffuse.o stochastic_physics.o MOM_check_scaling.o MOM_EOS_TEOS10.o gsw_mod_specvol_coefficients.o MOM_internal_tide_input.o MOM_PressureForce_FV.o coord_rho.o cvmix_shear.o MOM_EOS.o MOM_set_diffusivity.o MOM_checksum_packages.o MOM_load_love_numbers.o cvmix_utils.o gsw_specvol_first_derivatives.o gsw_mod_freezing_poly_coefficients.o MOM_ensemble_manager_infra.o MOM_EOS_Wright_red.o ISOMIP_tracer.o gsw_t_freezing_exact.o gsw_chem_potential_water_t_exact.o MOM_CVMix_ddiff.o MOM_io.o MOM_stochastics.o dumbbell_initialization.o MOM_lateral_boundary_diffusion.o MOM_get_input.o MOM_generic_tracer.o MOM_interface_filter.o MOM_horizontal_regridding.o user_revise_forcing.o MOM_dynamics_unsplit_RK2.o MOM_time_manager.o MOM_sum_output.o MOM_offline_aux.o write_ocean_obs.o MOM_temperature_convert.o MOM_continuity.o generic_tracer.o Idealized_Hurricane.o MOM_checksums.o gsw_specvol.o shelfwave_initialization.o MOM_diabatic_driver.o gsw_mod_teos10_constants.o MOM_grid_initialize.o gsw_mod_error_functions.o MOM_tidal_mixing.o MOM_database_comms.o dense_water_initialization.o MOM_controlled_forcing.o cvmix_ddiff.o MOM_tracer_types.o ocean_da_types.o coord_sigma.o MOM_cpu_clock.o MOM_PointAccel.o PLM_functions.o MOM_debugging.o coord_adapt.o MOM_bulk_mixed_layer.o tidal_bay_initialization.o MOM_EOS_Wright_full.o MOM_intrinsic_functions.o gsw_sp_from_sr.o gsw_rho_first_derivatives.o MOM_regridding.o dyed_obc_tracer.o gsw_gibbs_ice.o Neverworld_initialization.o MOM_continuity_PPM.o MOM_EOS_Wright.o MOM_EOS_linear.o MOM_error_handler.o MOM_offline_main.o gsw_t_from_ct.o gsw_ct_freezing_exact.o regrid_consts.o generic_tracer_utils.o ocean_da_core.o MOM_CoriolisAdv.o MOM_Zanna_Bolton.o MOM_restart.o MOM_unit_scaling.o MOM_dynamics_unsplit.o sloshing_initialization.o MOM_spherical_harmonics.o MOM_diagnostics.o MOM_remapping.o MOM_coord_initialization.o MOM_mixed_layer_restrat.o regrid_solvers.o gsw_ct_freezing_poly.o gsw_mod_toolbox.o MOM_hor_visc.o lock_exchange_initialization.o MOM_PressureForce.o MOM_particles.o DOME_tracer.o MOM_kappa_shear.o MOM_tracer_Z_init.o MOM_coms.o MOM_PressureForce_Montgomery.o MOM_vert_friction.o MOM_porous_barriers.o MOM_hybgen_remap.o MOM_ice_shelf_initialize.o P1M_functions.o polynomial_functions.o MOM_ice_shelf.o MOM_ensemble_manager.o Phillips_initialization.o MOM_diapyc_energy_req.o MOM_grid.o gsw_ct_from_pt.o MOM_internal_tides.o RGC_initialization.o coord_hycom.o cvmix_put_get.o MOM_constants.o MOM_domains.o MOM_hybgen_regrid.o gsw_mod_kinds.o MOM_oda_driver.o boundary_impulse_tracer.o MOM_data_override_infra.o DOME2d_initialization.o regrid_interp.o kdtree.o MOM_ALE.o PPM_functions.o MOM_tracer_registry.o MOM_netcdf.o dyed_channel_initialization.o user_initialization.o
OFF = $(SRCROOT)/MOM6/src/parameterizations/vertical/MOM_CVMix_ddiff.F90 $(SRCROOT)/MOM6/src/user/BFB_initialization.F90 $(SRCROOT)/MOM6/src/equation_of_state/TEOS10/gsw_pt0_from_t.f90 $(SRCROOT)/MOM6/src/ALE/PLM_functions.F90 $(SRCROOT)/MOM6/src/diagnostics/MOM_sum_output.F90 $(SRCROOT)/MOM6/src/tracer/advection_test_tracer.F90 $(SRCROOT)/MOM6/src/equation_of_state/TEOS10/gsw_t_from_ct.f90 $(SRCROOT)/MOM6/src/ALE/MOM_regridding.F90 $(SRCROOT)/MOM6/src/parameterizations/vertical/MOM_diapyc_energy_req.F90 $(SRCROOT)/MOM6/src/equation_of_state/TEOS10/gsw_gibbs.f90 $(SRCROOT)/MOM6/src/core/MOM_density_integrals.F90 $(SRCROOT)/MOM6/src/equation_of_state/MOM_EOS_Wright_full.F90 $(SRCROOT)/MOM6/src/user/DOME2d_initialization.F90 $(SRCROOT)/MOM6/src/ice_shelf/MOM_ice_shelf.F90 $(SRCROOT)/MOM6/src/parameterizations/stochastic/MOM_stochastics.F90 $(SRCROOT)/MOM6/src/equation_of_state/TEOS10/gsw_mod_teos10_constants.f90 $(SRCROOT)/MOM6/src/framework/MOM_unit_testing.F90 $(SRCROOT)/MOM6/config_src/external/database_comms/MOM_database_comms.F90 $(SRCROOT)/MOM6/src/user/external_gwave_initialization.F90 $(SRCROOT)MOM6/src/framework/version_variable.h $(SRCROOT)/MOM6/src/framework/MOM_io.F90 $(SRCROOT)/MOM6/src/tracer/RGC_tracer.F90 $(SRCROOT)/MOM6/src/framework/MOM_random.F90 $(SRCROOT)/MOM6/src/core/MOM_PressureForce.F90 $(SRCROOT)/MOM6/src/tracer/MOM_tracer_flow_control.F90 $(SRCROOT)/MOM6/src/user/tidal_bay_initialization.F90 $(SRCROOT)/MOM6/src/user/ISOMIP_initialization.F90 $(SRCROOT)/MOM6/src/parameterizations/lateral/MOM_load_love_numbers.F90 $(SRCROOT)/MOM6/src/user/lock_exchange_initialization.F90 $(SRCROOT)/MOM6/src/core/MOM_verticalGrid.F90 $(SRCROOT)/MOM6/config_src/infra/FMS2/MOM_error_infra.F90 $(SRCROOT)/MOM6/src/equation_of_state/MOM_EOS_linear.F90 $(SRCROOT)/MOM6/src/ALE/regrid_edge_values.F90 $(SRCROOT)/MOM6/src/tracer/ISOMIP_tracer.F90 $(SRCROOT)/MOM6/src/equation_of_state/TEOS10/gsw_mod_specvol_coefficients.f90 $(SRCROOT)/MOM6/src/framework/MOM_dyn_horgrid.F90 $(SRCROOT)/MOM6/src/user/dyed_obcs_initialization.F90 $(SRCROOT)/MOM6/src/core/MOM_grid.F90 $(SRCROOT)/MOM6/src/user/user_change_diffusivity.F90 $(SRCROOT)/MOM6/config_src/memory/dynamic_nonsymmetric/MOM_memory.h $(SRCROOT)/MOM6/src/initialization/MOM_grid_initialize.F90 $(SRCROOT)/MOM6/config_src/external/stochastic_physics/get_stochy_pattern.F90 $(SRCROOT)/MOM6/src/core/MOM_transcribe_grid.F90 $(SRCROOT)/MOM6/src/framework/MOM_array_transform.F90 $(SRCROOT)/MOM6/config_src/infra/FMS2/MOM_interp_infra.F90 $(SRCROOT)/MOM6/src/equation_of_state/TEOS10/gsw_specvol_first_derivatives.f90 $(SRCROOT)/MOM6/src/ice_shelf/user_shelf_init.F90 $(SRCROOT)/MOM6/src/parameterizations/CVmix/cvmix_math.F90 $(SRCROOT)/MOM6/src/user/benchmark_initialization.F90 $(SRCROOT)/MOM6/src/ice_shelf/MOM_ice_shelf_state.F90 $(SRCROOT)/MOM6/src/parameterizations/vertical/MOM_bkgnd_mixing.F90 $(SRCROOT)/MOM6/src/framework/version_variable.h $(SRCROOT)/MOM6/src/framework/MOM_unit_scaling.F90 $(SRCROOT)/MOM6/src/framework/MOM_coms.F90 $(SRCROOT)/MOM6/src/tracer/MOM_tracer_diabatic.F90 $(SRCROOT)/MOM6/src/framework/MOM_safe_alloc.F90 $(SRCROOT)/MOM6/src/ice_shelf/MOM_marine_ice.F90 $(SRCROOT)/MOM6/src/framework/MOM_netcdf.F90 $(SRCROOT)/MOM6/src/ALE/MOM_hybgen_regrid.F90 $(SRCROOT)/MOM6/src/ocean_data_assim/MOM_oda_incupd.F90 $(SRCROOT)/MOM6/src/user/basin_builder.F90 $(SRCROOT)/MOM6/src/core/MOM_dynamics_split_RK2.F90 $(SRCROOT)/MOM6/src/user/MOM_wave_interface.F90 $(SRCROOT)/MOM6/src/diagnostics/MOM_obsolete_diagnostics.F90 $(SRCROOT)/MOM6/src/equation_of_state/TEOS10/gsw_pt_from_t.f90 $(SRCROOT)/MOM6/src/ice_shelf/MOM_ice_shelf_diag_mediator.F90 $(SRCROOT)/MOM6/src/equation_of_state/TEOS10/gsw_rho_second_derivatives.f90 $(SRCROOT)/MOM6/src/diagnostics/MOM_diagnostics.F90 $(SRCROOT)/MOM6/src/parameterizations/vertical/MOM_set_diffusivity.F90 $(SRCROOT)/MOM6/src/core/MOM_open_boundary.F90 $(SRCROOT)/MOM6/src/user/DOME_initialization.F90 $(SRCROOT)/MOM6/src/user/dumbbell_initialization.F90 $(SRCROOT)/MOM6/src/equation_of_state/MOM_EOS_Wright.F90 $(SRCROOT)/MOM6/src/tracer/boundary_impulse_tracer.F90 $(SRCROOT)/MOM6/src/equation_of_state/TEOS10/gsw_ct_freezing_exact.f90 $(SRCROOT)/MOM6/config_src/infra/FMS2/MOM_data_override_infra.F90 $(SRCROOT)/MOM6/src/equation_of_state/TEOS10/gsw_gibbs_pt0_pt0.f90 $(SRCROOT)/MOM6/src/parameterizations/vertical/MOM_tidal_mixing.F90 $(SRCROOT)/MOM6/src/parameterizations/vertical/MOM_regularize_layers.F90 $(SRCROOT)/MOM6/src/parameterizations/vertical/MOM_sponge.F90 $(SRCROOT)/MOM6/src/parameterizations/CVmix/cvmix_tidal.F90 $(SRCROOT)/MOM6/src/parameterizations/vertical/MOM_CVMix_conv.F90 $(SRCROOT)/MOM6/src/parameterizations/vertical/MOM_entrain_diffusive.F90 $(SRCROOT)/MOM6/src/tracer/pseudo_salt_tracer.F90 $(SRCROOT)/MOM6/src/equation_of_state/TEOS10/gsw_mod_freezing_poly_coefficients.f90 $(SRCROOT)/MOM6/src/framework/MOM_hor_index.F90 $(SRCROOT)/MOM6/src/equation_of_state/MOM_EOS_Roquet_rho.F90 $(SRCROOT)/MOM6/src/core/MOM_porous_barriers.F90 $(SRCROOT)/MOM6/src/framework/MOM_intrinsic_functions.F90 $(SRCROOT)/MOM6/src/equation_of_state/MOM_EOS.F90 $(SRCROOT)/MOM6/src/equation_of_state/TEOS10/gsw_gibbs_ice.f90 $(SRCROOT)/MOM6/src/ice_shelf/MOM_ice_shelf_dynamics.F90 $(SRCROOT)/MOM6/src/equation_of_state/TEOS10/gsw_entropy_part.f90 $(SRCROOT)/MOM6/src/equation_of_state/MOM_EOS_Jackett06.F90 $(SRCROOT)/MOM6/src/equation_of_state/MOM_TFreeze.F90 $(SRCROOT)/MOM6/src/equation_of_state/TEOS10/gsw_ct_from_t.f90 $(SRCROOT)/MOM6/src/ALE/P1M_functions.F90 $(SRCROOT)/MOM6/src/equation_of_state/TEOS10/gsw_rho.f90 $(SRCROOT)/MOM6/src/framework/MOM_io_file.F90 $(SRCROOT)/MOM6/src/ALE/PPM_functions.F90 $(SRCROOT)/MOM6/src/initialization/MOM_coord_initialization.F90 $(SRCROOT)/MOM6/src/user/user_initialization.F90 $(SRCROOT)/MOM6/config_src/external/GFDL_ocean_BGC/generic_tracer_utils.F90 $(SRCROOT)/MOM6/src/equation_of_state/TEOS10/gsw_chem_potential_water_t_exact.f90 $(SRCROOT)/MOM6/src/initialization/MOM_fixed_initialization.F90 $(SRCROOT)/MOM6/src/framework/MOM_error_handler.F90 $(SRCROOT)/MOM6/src/user/Rossby_front_2d_initialization.F90 $(SRCROOT)/MOM6/src/framework/MOM_cpu_clock.F90 $(SRCROOT)/MOM6/src/parameterizations/lateral/MOM_Zanna_Bolton.F90 $(SRCROOT)/MOM6/src/user/circle_obcs_initialization.F90 $(SRCROOT)/MOM6/src/equation_of_state/TEOS10/gsw_ct_freezing_poly.f90 $(SRCROOT)/MOM6/src/core/MOM_PressureForce_FV.F90 $(SRCROOT)/MOM6/src/core/MOM_variables.F90 $(SRCROOT)/MOM6/src/tracer/MOM_CFC_cap.F90 $(SRCROOT)/MOM6/src/core/MOM.F90 $(SRCROOT)/MOM6/src/equation_of_state/TEOS10/gsw_t_freezing_poly.f90 $(SRCROOT)/MOM6/src/user/dyed_channel_initialization.F90 $(SRCROOT)/MOM6/src/core/MOM_continuity_PPM.F90 $(SRCROOT)/MOM6/src/equation_of_state/MOM_EOS_UNESCO.F90 $(SRCROOT)/MOM6/src/parameterizations/CVmix/cvmix_kinds_and_types.F90 $(SRCROOT)/MOM6/src/parameterizations/vertical/MOM_energetic_PBL.F90 $(SRCROOT)/MOM6/config_src/external/database_comms/database_client_interface.F90 $(SRCROOT)/MOM6/src/tracer/tracer_example.F90 $(SRCROOT)/MOM6/src/equation_of_state/TEOS10/gsw_mod_kinds.f90 $(SRCROOT)/MOM6/src/parameterizations/lateral/MOM_thickness_diffuse.F90 $(SRCROOT)/MOM6/src/core/MOM_checksum_packages.F90 $(SRCROOT)/MOM6/src/ALE/MOM_remapping.F90 $(SRCROOT)/MOM6/src/framework/MOM_diag_remap.F90 $(SRCROOT)/MOM6/src/ALE/MOM_hybgen_remap.F90 $(SRCROOT)/MOM6/src/ALE/coord_sigma.F90 $(SRCROOT)/MOM6/src/ocean_data_assim/MOM_oda_driver.F90 $(SRCROOT)/MOM6/src/core/MOM_unit_tests.F90 $(SRCROOT)/MOM6/config_src/external/stochastic_physics/stochastic_physics.F90 $(SRCROOT)/MOM6/src/framework/MOM_ensemble_manager.F90 $(SRCROOT)/MOM6/src/parameterizations/CVmix/cvmix_ddiff.F90 $(SRCROOT)/MOM6/src/user/sloshing_initialization.F90 $(SRCROOT)/MOM6/config_src/infra/FMS2/MOM_diag_manager_infra.F90 $(SRCROOT)/MOM6/src/equation_of_state/TEOS10/gsw_specvol.f90 $(SRCROOT)/MOM6/src/user/Idealized_Hurricane.F90 $(SRCROOT)/MOM6/config_src/infra/FMS2/MOM_coms_infra.F90 $(SRCROOT)/MOM6/src/equation_of_state/TEOS10/gsw_mod_gibbs_ice_coefficients.f90 $(SRCROOT)/MOM6/src/core/MOM_dynamics_unsplit.F90 $(SRCROOT)/MOM6/src/equation_of_state/TEOS10/gsw_mod_error_functions.f90 $(SRCROOT)/MOM6/src/parameterizations/CVmix/cvmix_kpp.F90 $(SRCROOT)/MOM6/src/tracer/ideal_age_example.F90 $(SRCROOT)/MOM6/src/parameterizations/CVmix/cvmix_put_get.F90 $(SRCROOT)/MOM6/src/ALE/coord_zlike.F90 $(SRCROOT)/MOM6/src/equation_of_state/TEOS10/gsw_sp_from_sr.f90 $(SRCROOT)/MOM6/src/equation_of_state/TEOS10/gsw_mod_toolbox.f90 $(SRCROOT)/MOM6/config_src/infra/FMS2/MOM_domain_infra.F90 $(SRCROOT)/MOM6/src/framework/MOM_get_input.F90 $(SRCROOT)/MOM6/src/parameterizations/CVmix/cvmix_utils.F90 $(SRCROOT)/MOM6/src/user/user_revise_forcing.F90 $(SRCROOT)/MOM6/config_src/external/drifters/MOM_particles_types.F90 $(SRCROOT)/MOM6/src/parameterizations/vertical/MOM_internal_tide_input.F90 $(SRCROOT)/MOM6/src/core/MOM_interface_heights.F90 $(SRCROOT)/MOM6/src/ALE/P3M_functions.F90 $(SRCROOT)/MOM6/src/user/BFB_surface_forcing.F90 $(SRCROOT)/MOM6/src/ALE/MOM_hybgen_unmix.F90 $(SRCROOT)/MOM6/src/framework/MOM_string_functions.F90 $(SRCROOT)/MOM6/src/framework/MOM_diag_mediator.F90 $(SRCROOT)/MOM6/src/ALE/remapping_attic.F90 $(SRCROOT)/MOM6/src/parameterizations/CVmix/cvmix_convection.F90 $(SRCROOT)/MOM6/src/ALE/PQM_functions.F90 $(SRCROOT)/MOM6/src/diagnostics/MOM_obsolete_params.F90 $(SRCROOT)/MOM6/src/user/soliton_initialization.F90 $(SRCROOT)/MOM6/src/core/MOM_CoriolisAdv.F90 $(SRCROOT)/MOM6/config_src/infra/FMS2/MOM_constants.F90 $(SRCROOT)/MOM6/src/core/MOM_barotropic.F90 $(SRCROOT)/MOM6/src/parameterizations/lateral/MOM_internal_tides.F90 $(SRCROOT)/MOM6/src/framework/posix.F90 $(SRCROOT)/MOM6/config_src/infra/FMS2/MOM_io_infra.F90 $(SRCROOT)/MOM6/src/initialization/MOM_state_initialization.F90 $(SRCROOT)/MOM6/src/ALE/regrid_solvers.F90 $(SRCROOT)/MOM6/src/parameterizations/vertical/MOM_diabatic_aux.F90 $(SRCROOT)/MOM6/src/tracer/MOM_offline_main.F90 $(SRCROOT)/MOM6/src/ALE/MOM_ALE.F90 $(SRCROOT)/MOM6/src/user/Kelvin_initialization.F90 $(SRCROOT)/MOM6/src/core/MOM_PressureForce_Montgomery.F90 $(SRCROOT)/MOM6/src/tracer/MOM_tracer_registry.F90 $(SRCROOT)/MOM6/src/user/seamount_initialization.F90 $(SRCROOT)/MOM6/src/parameterizations/vertical/MOM_full_convection.F90 $(SRCROOT)/MOM6/src/user/MOM_controlled_forcing.F90 $(SRCROOT)/MOM6/src/equation_of_state/TEOS10/gsw_sr_from_sp.f90 $(SRCROOT)/MOM6/src/initialization/MOM_shared_initialization.F90 $(SRCROOT)/MOM6/src/ice_shelf/MOM_ice_shelf_initialize.F90 $(SRCROOT)/MOM6/config_src/external/drifters/MOM_particles.F90 $(SRCROOT)/MOM6/src/parameterizations/CVmix/cvmix_shear.F90 $(SRCROOT)/MOM6/config_src/infra/FMS2/MOM_ensemble_manager_infra.F90 $(SRCROOT)/MOM6/src/parameterizations/vertical/MOM_geothermal.F90 $(SRCROOT)/MOM6/src/equation_of_state/MOM_temperature_convert.F90 $(SRCROOT)/MOM6/src/ALE/regrid_consts.F90 $(SRCROOT)/MOM6/src/parameterizations/vertical/MOM_CVMix_shear.F90 $(SRCROOT)/MOM6/src/initialization/MOM_tracer_initialization_from_Z.F90 $(SRCROOT)/MOM6/src/core/MOM_continuity.F90 $(SRCROOT)/MOM6/src/framework/MOM_memory_macros.h $(SRCROOT)/MOM6/src/framework/MOM_coupler_types.F90 $(SRCROOT)/MOM6/src/user/baroclinic_zone_initialization.F90 $(SRCROOT)/MOM6/config_src/drivers/FMS_cap/MOM_surface_forcing_gfdl.F90 $(SRCROOT)/MOM6/src/user/supercritical_initialization.F90 $(SRCROOT)/MOM6/src/tracer/MOM_tracer_types.F90 $(SRCROOT)/MOM6/src/user/shelfwave_initialization.F90 $(SRCROOT)/MOM6/src/parameterizations/lateral/MOM_lateral_mixing_coeffs.F90 $(SRCROOT)/MOM6/src/ALE/coord_rho.F90 $(SRCROOT)/MOM6/src/framework/MOM_file_parser.F90 $(SRCROOT)/MOM6/config_src/infra/FMS2/MOM_couplertype_infra.F90 $(SRCROOT)/MOM6/src/core/MOM_dynamics_unsplit_RK2.F90 $(SRCROOT)/MOM6/src/parameterizations/lateral/MOM_MEKE.F90 $(SRCROOT)/MOM6/src/framework/MOM_document.F90 $(SRCROOT)/MOM6/src/diagnostics/MOM_spatial_means.F90 $(SRCROOT)/MOM6/src/parameterizations/lateral/MOM_hor_visc.F90 $(SRCROOT)/MOM6/src/diagnostics/MOM_debugging.F90 $(SRCROOT)/MOM6/src/core/MOM_forcing_type.F90 $(SRCROOT)/MOM6/src/parameterizations/lateral/MOM_spherical_harmonics.F90 $(SRCROOT)/MOM6/src/tracer/dyed_obc_tracer.F90 $(SRCROOT)/MOM6/src/equation_of_state/TEOS10/gsw_t_deriv_chem_potential_water_t_exact.f90 $(SRCROOT)/MOM6/config_src/external/ODA_hooks/ocean_da_types.F90 $(SRCROOT)/MOM6/src/equation_of_state/TEOS10/gsw_pt_from_ct.f90 $(SRCROOT)/MOM6/src/equation_of_state/TEOS10/gsw_ct_from_pt.f90 $(SRCROOT)/MOM6/config_src/infra/FMS2/MOM_cpu_clock_infra.F90 $(SRCROOT)/MOM6/src/framework/MOM_checksums.F90 $(SRCROOT)/MOM6/src/parameterizations/vertical/MOM_kappa_shear.F90 $(SRCROOT)/MOM6/config_src/external/GFDL_ocean_BGC/FMS_coupler_util.F90 $(SRCROOT)/MOM6/src/equation_of_state/MOM_EOS_Roquet_SpV.F90 $(SRCROOT)/MOM6/src/equation_of_state/MOM_EOS_TEOS10.F90 $(SRCROOT)/MOM6/src/diagnostics/MOM_wave_speed.F90 $(SRCROOT)/MOM6/config_src/external/ODA_hooks/kdtree.f90 $(SRCROOT)/MOM6/src/user/Phillips_initialization.F90 $(SRCROOT)/MOM6/src/tracer/MOM_tracer_Z_init.F90 $(SRCROOT)/MOM6/src/equation_of_state/MOM_EOS_Wright_red.F90 $(SRCROOT)/MOM6/src/tracer/MOM_offline_aux.F90 $(SRCROOT)/MOM6/src/parameterizations/lateral/MOM_interface_filter.F90 $(SRCROOT)/MOM6/src/framework/MOM_interpolate.F90 $(SRCROOT)/MOM6/src/parameterizations/vertical/MOM_CVMix_KPP.F90 $(SRCROOT)/MOM6/src/framework/MOM_data_override.F90 $(SRCROOT)/MOM6/src/equation_of_state/TEOS10/gsw_rho_first_derivatives.f90 $(SRCROOT)/MOM6/src/tracer/MOM_tracer_advect.F90 $(SRCROOT)/MOM6/config_src/infra/FMS2/MOM_time_manager.F90 $(SRCROOT)/MOM6/src/parameterizations/lateral/MOM_MEKE_types.F90 $(SRCROOT)/MOM6/src/core/MOM_isopycnal_slopes.F90 $(SRCROOT)/MOM6/src/tracer/oil_tracer.F90 $(SRCROOT)/MOM6/src/equation_of_state/TEOS10/gsw_entropy_part_zerop.f90 $(SRCROOT)/MOM6/src/ALE/PCM_functions.F90 $(SRCROOT)/MOM6/src/framework/MOM_write_cputime.F90 $(SRCROOT)/MOM6/src/equation_of_state/TEOS10/gsw_t_freezing_exact.f90 $(SRCROOT)/MOM6/src/tracer/DOME_tracer.F90 $(SRCROOT)/MOM6/src/user/dumbbell_surface_forcing.F90 $(SRCROOT)/MOM6/src/parameterizations/lateral/MOM_tidal_forcing.F90 $(SRCROOT)/MOM6/src/ALE/coord_hycom.F90 $(SRCROOT)/MOM6/src/tracer/MOM_generic_tracer.F90 $(SRCROOT)/MOM6/src/framework/testing/MOM_file_parser_tests.F90 $(SRCROOT)/MOM6/src/ALE/polynomial_functions.F90 $(SRCROOT)/MOM6/src/framework/MOM_unique_scales.F90 $(SRCROOT)/MOM6/src/core/MOM_boundary_update.F90 $(SRCROOT)/MOM6/src/tracer/MOM_lateral_boundary_diffusion.F90 $(SRCROOT)/MOM6/src/core/MOM_check_scaling.F90 $(SRCROOT)/MOM6/src/tracer/MOM_neutral_diffusion.F90 $(SRCROOT)/MOM6/config_src/drivers/FMS_cap/ocean_model_MOM.F90 $(SRCROOT)/MOM6/src/tracer/nw2_tracers.F90 $(SRCROOT)/MOM6/src/parameterizations/vertical/MOM_opacity.F90 $(SRCROOT)/MOM6/src/framework/MOM_restart.F90 $(SRCROOT)/MOM6/src/framework/MOM_domains.F90 $(SRCROOT)/MOM6/src/core/MOM_stoch_eos.F90 $(SRCROOT)/MOM6/src/framework/posix.h $(SRCROOT)/MOM6/src/framework/MOM_horizontal_regridding.F90 $(SRCROOT)/MOM6/src/parameterizations/vertical/MOM_ALE_sponge.F90 $(SRCROOT)/MOM6/config_src/external/ODA_hooks/ocean_da_core.F90 $(SRCROOT)/MOM6/src/parameterizations/vertical/MOM_bulk_mixed_layer.F90 $(SRCROOT)/MOM6/src/ALE/regrid_interp.F90 $(SRCROOT)/MOM6/src/parameterizations/vertical/MOM_set_viscosity.F90 $(SRCROOT)/MOM6/config_src/external/ODA_hooks/write_ocean_obs.F90 $(SRCROOT)/MOM6/src/tracer/dye_example.F90 $(SRCROOT)/MOM6/src/ALE/coord_adapt.F90 $(SRCROOT)/MOM6/src/parameterizations/lateral/MOM_mixed_layer_restrat.F90 $(SRCROOT)/MOM6/src/user/SCM_CVMix_tests.F90 $(SRCROOT)/MOM6/src/equation_of_state/TEOS10/gsw_specvol_second_derivatives.f90 $(SRCROOT)/MOM6/src/parameterizations/vertical/MOM_vert_friction.F90 $(SRCROOT)/MOM6/config_src/external/GFDL_ocean_BGC/generic_tracer.F90 $(SRCROOT)/MOM6/src/user/Neverworld_initialization.F90 $(SRCROOT)/MOM6/src/user/RGC_initialization.F90 $(SRCROOT)/MOM6/src/diagnostics/MOM_PointAccel.F90 $(SRCROOT)/MOM6/src/tracer/MOM_OCMIP2_CFC.F90 $(SRCROOT)/MOM6/src/user/dense_water_initialization.F90 $(SRCROOT)/MOM6/src/parameterizations/vertical/MOM_diabatic_driver.F90 $(SRCROOT)/MOM6/src/user/adjustment_initialization.F90 $(SRCROOT)/MOM6/src/parameterizations/CVmix/cvmix_background.F90 $(SRCROOT)/MOM6/src/tracer/MOM_tracer_hor_diff.F90
clean: neat
	-rm -f .libMOM6.a.cppdefs $(OBJ) libMOM6.a

neat:
	-rm -f $(TMPFILES)

localize: $(OFF)
	cp $(OFF) .

TAGS: $(SRC)
	etags $(SRC)

tags: $(SRC)
	ctags $(SRC)

libMOM6.a: $(OBJ)
	$(AR) $(ARFLAGS) libMOM6.a $(OBJ)

