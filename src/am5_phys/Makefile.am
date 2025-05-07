ACLOCAL_AMFLAGS = -I m4
AM_CPPFLAGS = -I${top_srcdir}/atmos_param/cosp -DINTERNAL_FILE_NML

lib_LTLIBRARIES = libatmos_phys.la
libatmos_phys_la_SOURCES = ${top_srcdir}/atmos_shared/atmos_nudge/atmos_nudge.F90 \
		${top_srcdir}/atmos_shared/atmos_cmip_diag/atmos_global_diag.F90 \
		${top_srcdir}/atmos_shared/atmos_cmip_diag/atmos_cmip_diag.F90 \
		${top_srcdir}/atmos_shared/vert_advection/vert_advection.F90 \
		${top_srcdir}/atmos_shared/aerosol/aerosol.F90 \
		${top_srcdir}/atmos_shared/aerosol/aerosol_types.F90 \
		${top_srcdir}/atmos_shared/tracer_driver/atmos_convection_tracer.F90 \
		${top_srcdir}/atmos_shared/tracer_driver/atmos_sulfur_hex.F90 \
		${top_srcdir}/atmos_shared/tracer_driver/atmos_tracer_driver.F90 \
		${top_srcdir}/atmos_shared/tracer_driver/atmos_sea_salt.F90 \
		${top_srcdir}/atmos_shared/tracer_driver/atmos_tracer_utilities.F90 \
		${top_srcdir}/atmos_shared/tracer_driver/atmos_co2.F90 \
		${top_srcdir}/atmos_shared/tracer_driver/tropchem_driver.F90 \
		${top_srcdir}/atmos_shared/tracer_driver/atmos_dust.F90 \
		${top_srcdir}/atmos_shared/tracer_driver/atmos_tropopause.F90 \
		${top_srcdir}/atmos_shared/tracer_driver/atmos_sulfate.F90 \
		${top_srcdir}/atmos_shared/tracer_driver/atmos_radon.F90 \
		${top_srcdir}/atmos_shared/tracer_driver/atmos_nh3_tag.F90 \
		${top_srcdir}/atmos_shared/tracer_driver/atmos_soa.F90 \
		${top_srcdir}/atmos_shared/tracer_driver/atmos_ch4.F90 \
		${top_srcdir}/atmos_shared/tracer_driver/atmos_regional_tracer_driver.F90 \
		${top_srcdir}/atmos_shared/tracer_driver/cloud_chem.F90 \
		${top_srcdir}/atmos_shared/tracer_driver/atmos_age_tracer.F90 \
		${top_srcdir}/atmos_shared/tracer_driver/mo_setsox.F90 \
		${top_srcdir}/atmos_shared/tracer_driver/atmos_ch3i.F90 \
		${top_srcdir}/atmos_shared/tracer_driver/xactive_bvoc.F90 \
		${top_srcdir}/atmos_shared/tracer_driver/atmos_carbon_aerosol.F90 \
		${top_srcdir}/atmos_shared/tracer_driver/stratchem/strat_chem_model.F90 \
		${top_srcdir}/atmos_shared/tracer_driver/stratchem/strat_chem_driver.F90 \
		${top_srcdir}/atmos_shared/tracer_driver/tropchem/mo_chemini.F90 \
		${top_srcdir}/atmos_shared/tracer_driver/tropchem/mo_jpl.F90 \
		${top_srcdir}/atmos_shared/tracer_driver/tropchem/mo_exp_slv.F90 \
		${top_srcdir}/atmos_shared/tracer_driver/tropchem/mo_chem_utls.F90 \
		${top_srcdir}/atmos_shared/tracer_driver/tropchem/tropchem_types_mod.F90 \
		${top_srcdir}/atmos_shared/tracer_driver/tropchem/mo_usrrxt.F90 \
		${top_srcdir}/atmos_shared/tracer_driver/tropchem/mo_photo.F90 \
		${top_srcdir}/atmos_shared/tracer_driver/tropchem/aerosol_thermodynamics.F90 \
		${top_srcdir}/atmos_shared/tracer_driver/tropchem/mo_chemdr.F90 \
		${top_srcdir}/atmos_shared/tracer_driver/tropchem/strat_chem_utilities.F90 \
		${top_srcdir}/atmos_shared/tracer_driver/tropchem/mo_imp_slv.F90 \
		${top_srcdir}/atmos_shared/tracer_driver/tropchem/isoropiaIIcode.F \
		${top_srcdir}/atmos_shared/tracer_driver/tropchem/mo_hook.F90 \
		${top_srcdir}/atmos_shared/tracer_driver/tropchem/mo_setinv.F90 \
		${top_srcdir}/atmos_shared/tracer_driver/tropchem/mo_rodas_slv.F90 \
		${top_srcdir}/atmos_shared/tracer_driver/tropchem/mo_read_sim_chm.F90 \
		${top_srcdir}/atmos_shared/tracer_driver/tropchem/m_tracname.F90 \
		${top_srcdir}/atmos_shared/tracer_driver/tropchem/AM4_chem/moz.mods.F90 \
		${top_srcdir}/atmos_shared/tracer_driver/tropchem/AM4_chem/mo_fphoto.F90 \
		${top_srcdir}/atmos_shared/tracer_driver/tropchem/AM4_chem/moz.mat.F90 \
		${top_srcdir}/atmos_shared/tracer_driver/tropchem/AM4_chem/moz.subs.F90 \
		${top_srcdir}/atmos_shared/tracer_driver/tropchem/AM4_chem/mo_fastjx.F90 \
		${top_srcdir}/atmos_shared/tracer_driver/aer_ccn_act/aer_ccn_act.F90 \
		${top_srcdir}/atmos_shared/tracer_driver/aer_ccn_act/aer_in_act.F90 \
		${top_srcdir}/atmos_shared/tracer_driver/aer_ccn_act/aer_ccn_act_k.F90 \
		${top_srcdir}/atmos_shared/tracer_driver/aer_ccn_act/ice_nucl.F90 \
		${top_srcdir}/atmos_shared/tracer_driver/aer_ccn_act/aerosol_params.F90 \
		${top_srcdir}/atmos_shared/tracer_driver/atmos_ozone_tracer.F90 \
		${top_srcdir}/atmos_param/physics_driver/physics_types.F90 \
		${top_srcdir}/atmos_param/physics_driver/physics_driver.F90 \
		${top_srcdir}/atmos_param/topo_drag/topo_drag.F90 \
		${top_srcdir}/atmos_param/macrophysics/ls_cloud_macrophysics.F90 \
		${top_srcdir}/atmos_param/macrophysics/tiedtke_macro.F90 \
		${top_srcdir}/atmos_param/radiation/clouds/cloudrad_types.F90 \
		${top_srcdir}/atmos_param/radiation/clouds/microphys_rad.F90 \
		${top_srcdir}/atmos_param/radiation/clouds/strat_clouds_W.F90 \
		${top_srcdir}/atmos_param/radiation/clouds/isccp_clouds.F90 \
		${top_srcdir}/atmos_param/radiation/clouds/get_random_number_stream.F90 \
		${top_srcdir}/atmos_param/radiation/clouds/uw_clouds_W.F90 \
		${top_srcdir}/atmos_param/radiation/clouds/cloudrad_driver.F90 \
		${top_srcdir}/atmos_param/radiation/clouds/cloud_spec.F90 \
		${top_srcdir}/atmos_param/radiation/clouds/cloudrad_diagnostics.F90 \
		${top_srcdir}/atmos_param/radiation/clouds/cloudrad_package.F90 \
		${top_srcdir}/atmos_param/radiation/longwave/longwave_tables.F90 \
		${top_srcdir}/atmos_param/radiation/longwave/longwave_clouds.F90 \
		${top_srcdir}/atmos_param/radiation/longwave/gas_tf.F90 \
		${top_srcdir}/atmos_param/radiation/longwave/longwave_params.F90 \
		${top_srcdir}/atmos_param/radiation/longwave/sealw99.F90 \
		${top_srcdir}/atmos_param/radiation/longwave/optical_path.F90 \
		${top_srcdir}/atmos_param/radiation/longwave/longwave_types.F90 \
		${top_srcdir}/atmos_param/radiation/longwave/longwave_utilities.F90 \
		${top_srcdir}/atmos_param/radiation/longwave/longwave_fluxes.F90 \
		${top_srcdir}/atmos_param/radiation/longwave/lw_gases_stdtf.F90 \
		${top_srcdir}/atmos_param/radiation/aerosols/aerosolrad_types.F90 \
		${top_srcdir}/atmos_param/radiation/aerosols/aerosolrad_package.F90 \
		${top_srcdir}/atmos_param/radiation/aerosols/aerosolrad_driver.F90 \
		${top_srcdir}/atmos_param/radiation/shortwave/esfsw_driver.F90 \
		${top_srcdir}/atmos_param/radiation/shortwave/esfsw_utilities.F90 \
		${top_srcdir}/atmos_param/radiation/shortwave/shortwave_types.F90 \
		${top_srcdir}/atmos_param/radiation/shortwave/esfsw_bands.F90 \
		${top_srcdir}/atmos_param/radiation/shortwave/esfsw_parameters.F90 \
		${top_srcdir}/atmos_param/radiation/radiative_gases/radiative_gases.F90 \
		${top_srcdir}/atmos_param/radiation/radiative_gases/radiative_gases_types.F90 \
		${top_srcdir}/atmos_param/radiation/radiative_gases/ozone.F90 \
		${top_srcdir}/atmos_param/radiation/util/radiation_diag.F90 \
		${top_srcdir}/atmos_param/radiation/driver/longwave_driver.F90 \
		${top_srcdir}/atmos_param/radiation/driver/rad_output_file.F90 \
		${top_srcdir}/atmos_param/radiation/driver/radiation_driver_types.F90 \
		${top_srcdir}/atmos_param/radiation/driver/radiation_types.F90 \
		${top_srcdir}/atmos_param/radiation/driver/solar_data_driver.F90 \
		${top_srcdir}/atmos_param/radiation/driver/shortwave_driver.F90 \
		${top_srcdir}/atmos_param/radiation/driver/radiation_driver_diag.F90 \
		${top_srcdir}/atmos_param/radiation/driver/radiation_driver.F90 \
		${top_srcdir}/atmos_param/convection_driver/convection_utilities.F90 \
		${top_srcdir}/atmos_param/convection_driver/detr_ice_num.F90 \
		${top_srcdir}/atmos_param/convection_driver/convection_driver.F90 \
		${top_srcdir}/atmos_param/physics_radiation_exch/physics_radiation_exch.F90 \
		${top_srcdir}/atmos_param/cloud_generator/betaDistribution.F90 \
		${top_srcdir}/atmos_param/cloud_generator/cloud_generator.F90 \
		${top_srcdir}/atmos_param/lscloud_driver/lscloud_debug.F90 \
		${top_srcdir}/atmos_param/lscloud_driver/lscloud_constants.F90 \
		${top_srcdir}/atmos_param/lscloud_driver/lscloud_netcdf.F90 \
		${top_srcdir}/atmos_param/lscloud_driver/lscloud_types.F90 \
		${top_srcdir}/atmos_param/lscloud_driver/lscloud_driver.F90 \
		${top_srcdir}/atmos_param/lscloud_driver/polysvp.F90 \
		${top_srcdir}/atmos_param/lscloud_driver/check_nan.F90 \
		${top_srcdir}/atmos_param/stable_bl_turb/stable_bl_turb.F90 \
		${top_srcdir}/atmos_param/microphysics/gamma_mg.F90 \
		${top_srcdir}/atmos_param/microphysics/micro_mg2.F90 \
		${top_srcdir}/atmos_param/microphysics/micro_mg2_utils.F90 \
		${top_srcdir}/atmos_param/microphysics/rotstayn_klein_mp.F90 \
		${top_srcdir}/atmos_param/microphysics/ls_cloud_microphysics.F90 \
		${top_srcdir}/atmos_param/damping_driver/damping_driver.F90 \
		${top_srcdir}/atmos_param/vert_diff_driver/vert_diff_driver.F90 \
		${top_srcdir}/atmos_param/moist_processes/moist_processes.F90 \
		${top_srcdir}/atmos_param/moist_processes/moist_processes_utils.F90 \
		${top_srcdir}/atmos_param/cosp/cosp_modis_simulator.F90 \
		${top_srcdir}/atmos_param/cosp/cosp_driver.F90 \
		${top_srcdir}/atmos_param/cosp/cosp_misr_simulator.F90 \
		${top_srcdir}/atmos_param/cosp/cosp_stats.F90 \
		${top_srcdir}/atmos_param/cosp/cosp_utils.F90 \
		${top_srcdir}/atmos_param/cosp/cosp_types.F90 \
		${top_srcdir}/atmos_param/cosp/cosp_rttov_simulator.F90 \
		${top_srcdir}/atmos_param/cosp/cosp_isccp_simulator.F90 \
		${top_srcdir}/atmos_param/cosp/cosp.F90 \
		${top_srcdir}/atmos_param/cosp/cosp_constants.F90 \
		${top_srcdir}/atmos_param/cosp/cosp_diagnostics.F90 \
		${top_srcdir}/atmos_param/cosp/cosp_simulator.F90 \
		${top_srcdir}/atmos_param/cosp/cosp_lidar.F90 \
		${top_srcdir}/atmos_param/cosp/cosp_io.F90 \
		${top_srcdir}/atmos_param/cosp/actsim/lidar_simulator.F90 \
		${top_srcdir}/atmos_param/cosp/actsim/lmd_ipsl_stats.F90 \
		${top_srcdir}/atmos_param/cosp/MISR_simulator/MISR_simulator.F \
		${top_srcdir}/atmos_param/cosp/quickbeam/radar_simulator_types.F90 \
		${top_srcdir}/atmos_param/cosp/quickbeam/array_lib.F90 \
		${top_srcdir}/atmos_param/cosp/quickbeam/dsd.F90 \
		${top_srcdir}/atmos_param/cosp/quickbeam/optics_lib.F90 \
		${top_srcdir}/atmos_param/cosp/quickbeam/gases.F90 \
		${top_srcdir}/atmos_param/cosp/quickbeam/format_input.F90 \
		${top_srcdir}/atmos_param/cosp/quickbeam/mrgrnk.F90 \
		${top_srcdir}/atmos_param/cosp/quickbeam/radar_simulator.F90 \
		${top_srcdir}/atmos_param/cosp/quickbeam/scale_LUTs_io.F90 \
		${top_srcdir}/atmos_param/cosp/quickbeam/zeff.F90 \
		${top_srcdir}/atmos_param/cosp/quickbeam/atmos_lib.F90 \
		${top_srcdir}/atmos_param/cosp/quickbeam/calc_Re.F90 \
		${top_srcdir}/atmos_param/cosp/quickbeam/radar_simulator_init.F90 \
		${top_srcdir}/atmos_param/cosp/quickbeam/math_lib.F90 \
		${top_srcdir}/atmos_param/cosp/llnl/prec_scops.F \
		${top_srcdir}/atmos_param/cosp/llnl/cosp_radar.F90 \
		${top_srcdir}/atmos_param/cosp/llnl/llnl_stats.F90 \
		${top_srcdir}/atmos_param/cosp/llnl/pf_to_mr.F \
		${top_srcdir}/atmos_param/cosp/icarus-scops-3.7/icarus.F90 \
		${top_srcdir}/atmos_param/cosp/icarus-scops-3.7/scops.F \
		${top_srcdir}/atmos_param/cosp/icarus-scops-3.7/isccp_cloud_types.F \
		${top_srcdir}/atmos_param/cosp/MODIS_simulator/modis_simulator.F90 \
		${top_srcdir}/atmos_param/lscale_cond/lscale_cond.F90 \
		${top_srcdir}/atmos_param/entrain/entrain.F90 \
		${top_srcdir}/atmos_param/cg_drag/cg_drag.F90 \
		${top_srcdir}/atmos_param/cu_mo_trans/cu_mo_trans.F90 \
		${top_srcdir}/atmos_param/cloud_rad/cloud_rad.F90 \
    ${top_srcdir}/atmos_param/gw_beres/gw_beres.F90 \
    ${top_srcdir}/atmos_param/gw_beres/gw_common_utils.F90 \
		${top_srcdir}/atmos_param/vert_diff/vert_diff.F90 \
		${top_srcdir}/atmos_param/aerosol_cloud/aerosol_cloud.F90 \
		${top_srcdir}/atmos_param/mg_drag/mg_drag.F90 \
		${top_srcdir}/atmos_param/vert_turb_driver/vert_turb_driver.F90 \
    ${top_srcdir}/atmos_param/shallow_cu/wet_deposition_0D.F90 \
		${top_srcdir}/atmos_param/shallow_cu/conv_plumes_k.F90 \
		${top_srcdir}/atmos_param/shallow_cu/deep_conv.F90 \
		${top_srcdir}/atmos_param/shallow_cu/conv_closures.F90 \
		${top_srcdir}/atmos_param/shallow_cu/conv_plumes.F90 \
		${top_srcdir}/atmos_param/shallow_cu/conv_utilities_k.F90 \
		${top_srcdir}/atmos_param/shallow_cu/conv_utilities.F90 \
		${top_srcdir}/atmos_param/shallow_cu/uw_conv.F90 \
    ${top_srcdir}/atmos_param/han_edmf/hanedmf.F90 \
    ${top_srcdir}/atmos_param/han_edmf/mfpbltq.F90 \
    ${top_srcdir}/atmos_param/han_edmf/mfscuq.F90 \
    ${top_srcdir}/atmos_param/han_edmf/satmedmfvdifq.F90 \
    ${top_srcdir}/atmos_param/han_edmf/tridi.F90

${top_builddir}/atmos_shared/atmos_nudge/atmos_nudge.lo: ${top_srcdir}/atmos_shared/atmos_nudge/atmos_nudge.F90
${top_builddir}/atmos_shared/atmos_cmip_diag/atmos_global_diag.lo: ${top_srcdir}/atmos_shared/atmos_cmip_diag/atmos_global_diag.F90
${top_builddir}/atmos_shared/atmos_cmip_diag/atmos_cmip_diag.lo: ${top_srcdir}/atmos_shared/atmos_cmip_diag/atmos_cmip_diag.F90
${top_builddir}/atmos_shared/vert_advection/vert_advection.lo: ${top_srcdir}/atmos_shared/vert_advection/vert_advection.F90
${top_builddir}/atmos_shared/aerosol/aerosol.lo: ${top_srcdir}/atmos_shared/aerosol/aerosol.F90 \
		 ${top_srcdir}/atmos_shared/aerosol/aerosol_types.lo
${top_builddir}/atmos_shared/aerosol/aerosol_types.lo: ${top_srcdir}/atmos_shared/aerosol/aerosol_types.F90
${top_builddir}/atmos_shared/tracer_driver/atmos_convection_tracer.lo: ${top_srcdir}/atmos_shared/tracer_driver/atmos_convection_tracer.F90 \
		 ${top_srcdir}/atmos_shared/tracer_driver/atmos_tracer_utilities.lo
${top_builddir}/atmos_shared/tracer_driver/atmos_sulfur_hex.lo: ${top_srcdir}/atmos_shared/tracer_driver/atmos_sulfur_hex.F90 \
		 ${top_srcdir}/atmos_shared/tracer_driver/atmos_tracer_utilities.lo
${top_builddir}/atmos_shared/tracer_driver/atmos_tracer_driver.lo: ${top_srcdir}/atmos_shared/tracer_driver/atmos_tracer_driver.F90 \
		 ${top_srcdir}/atmos_shared/atmos_cmip_diag/atmos_cmip_diag.lo \
		 ${top_srcdir}/atmos_shared/tracer_driver/atmos_tracer_utilities.lo \
		 ${top_srcdir}/atmos_shared/tracer_driver/atmos_radon.lo \
		 ${top_srcdir}/atmos_shared/tracer_driver/atmos_carbon_aerosol.lo \
		 ${top_srcdir}/atmos_shared/tracer_driver/atmos_convection_tracer.lo \
		 ${top_srcdir}/atmos_shared/tracer_driver/atmos_sulfur_hex.lo \
		 ${top_srcdir}/atmos_shared/tracer_driver/atmos_ch3i.lo \
		 ${top_srcdir}/atmos_shared/tracer_driver/atmos_sea_salt.lo \
		 ${top_srcdir}/atmos_shared/tracer_driver/atmos_dust.lo \
		 ${top_srcdir}/atmos_shared/tracer_driver/atmos_sulfate.lo \
		 ${top_srcdir}/atmos_shared/tracer_driver/atmos_soa.lo \
		 ${top_srcdir}/atmos_shared/tracer_driver/atmos_nh3_tag.lo \
		 ${top_srcdir}/atmos_shared/tracer_driver/tropchem_driver.lo \
		 ${top_srcdir}/atmos_shared/tracer_driver/atmos_regional_tracer_driver.lo \
		 ${top_srcdir}/atmos_shared/tracer_driver/stratchem/strat_chem_driver.lo \
		 ${top_srcdir}/atmos_shared/tracer_driver/atmos_age_tracer.lo \
		 ${top_srcdir}/atmos_shared/tracer_driver/atmos_co2.lo \
		 ${top_srcdir}/atmos_shared/tracer_driver/atmos_ch4.lo \
		 ${top_srcdir}/atmos_shared/tracer_driver/atmos_tropopause.lo \
		 ${top_srcdir}/atmos_shared/tracer_driver/xactive_bvoc.lo \
		 ${top_srcdir}/atmos_shared/tracer_driver/atmos_ozone_tracer.lo
${top_builddir}/atmos_shared/tracer_driver/atmos_sea_salt.lo: ${top_srcdir}/atmos_shared/tracer_driver/atmos_sea_salt.F90 \
		 ${top_srcdir}/atmos_shared/atmos_cmip_diag/atmos_cmip_diag.lo \
		 ${top_srcdir}/atmos_shared/tracer_driver/atmos_tracer_utilities.lo
${top_builddir}/atmos_shared/tracer_driver/atmos_tracer_utilities.lo: ${top_srcdir}/atmos_shared/tracer_driver/atmos_tracer_utilities.F90 \
		 ${top_srcdir}/atmos_shared/atmos_cmip_diag/atmos_cmip_diag.lo
${top_builddir}/atmos_shared/tracer_driver/atmos_co2.lo: ${top_srcdir}/atmos_shared/tracer_driver/atmos_co2.F90
${top_builddir}/atmos_shared/tracer_driver/tropchem_driver.lo: ${top_srcdir}/atmos_shared/tracer_driver/tropchem_driver.F90 \
		 ${top_srcdir}/atmos_shared/tracer_driver/tropchem/tropchem_types_mod.lo \
		 ${top_srcdir}/atmos_shared/atmos_cmip_diag/atmos_cmip_diag.lo \
		 ${top_srcdir}/atmos_shared/tracer_driver/atmos_tracer_utilities.lo \
		 ${top_srcdir}/atmos_shared/tracer_driver/tropchem/mo_chemdr.lo \
		 ${top_srcdir}/atmos_shared/tracer_driver/mo_setsox.lo \
		 ${top_srcdir}/atmos_shared/tracer_driver/tropchem/mo_chemini.lo \
		 ${top_srcdir}/atmos_shared/tracer_driver/tropchem/m_tracname.lo \
		 ${top_srcdir}/atmos_shared/tracer_driver/tropchem/AM4_chem/moz.mods.lo \
		 ${top_srcdir}/atmos_shared/tracer_driver/tropchem/mo_hook.lo \
		 ${top_srcdir}/atmos_shared/tracer_driver/tropchem/strat_chem_utilities.lo \
		 ${top_srcdir}/atmos_shared/tracer_driver/tropchem/mo_chem_utls.lo \
		 ${top_srcdir}/atmos_shared/tracer_driver/atmos_sulfate.lo \
		 ${top_srcdir}/atmos_param/radiation/driver/shortwave_driver.lo \
		 ${top_srcdir}/atmos_shared/tracer_driver/cloud_chem.lo \
		 ${top_srcdir}/atmos_shared/tracer_driver/tropchem/aerosol_thermodynamics.lo \
		 ${top_srcdir}/atmos_shared/tracer_driver/tropchem/mo_usrrxt.lo
${top_builddir}/atmos_shared/tracer_driver/atmos_dust.lo: ${top_srcdir}/atmos_shared/tracer_driver/atmos_dust.F90 \
		 ${top_srcdir}/atmos_shared/atmos_cmip_diag/atmos_cmip_diag.lo \
		 ${top_srcdir}/atmos_shared/tracer_driver/atmos_tracer_utilities.lo
${top_builddir}/atmos_shared/tracer_driver/atmos_tropopause.lo: ${top_srcdir}/atmos_shared/tracer_driver/atmos_tropopause.F90 \
		 ${top_srcdir}/atmos_shared/atmos_cmip_diag/atmos_cmip_diag.lo
${top_builddir}/atmos_shared/tracer_driver/atmos_sulfate.lo: ${top_srcdir}/atmos_shared/tracer_driver/atmos_sulfate.F90 \
		 ${top_srcdir}/atmos_shared/atmos_cmip_diag/atmos_cmip_diag.lo \
		 ${top_srcdir}/atmos_shared/tracer_driver/cloud_chem.lo
${top_builddir}/atmos_shared/tracer_driver/atmos_radon.lo: ${top_srcdir}/atmos_shared/tracer_driver/atmos_radon.F90 \
		 ${top_srcdir}/atmos_shared/tracer_driver/atmos_tracer_utilities.lo
${top_builddir}/atmos_shared/tracer_driver/atmos_nh3_tag.lo: ${top_srcdir}/atmos_shared/tracer_driver/atmos_nh3_tag.F90 \
		 ${top_srcdir}/atmos_shared/atmos_cmip_diag/atmos_cmip_diag.lo \
		 ${top_srcdir}/atmos_shared/tracer_driver/atmos_tracer_utilities.lo
${top_builddir}/atmos_shared/tracer_driver/atmos_soa.lo: ${top_srcdir}/atmos_shared/tracer_driver/atmos_soa.F90 \
		 ${top_srcdir}/atmos_shared/atmos_cmip_diag/atmos_cmip_diag.lo \
		 ${top_srcdir}/atmos_shared/tracer_driver/xactive_bvoc.lo
${top_builddir}/atmos_shared/tracer_driver/atmos_ch4.lo: ${top_srcdir}/atmos_shared/tracer_driver/atmos_ch4.F90
${top_builddir}/atmos_shared/tracer_driver/atmos_regional_tracer_driver.lo: ${top_srcdir}/atmos_shared/tracer_driver/atmos_regional_tracer_driver.F90
${top_builddir}/atmos_shared/tracer_driver/cloud_chem.lo: ${top_srcdir}/atmos_shared/tracer_driver/cloud_chem.F90
${top_builddir}/atmos_shared/tracer_driver/atmos_age_tracer.lo: ${top_srcdir}/atmos_shared/tracer_driver/atmos_age_tracer.F90
${top_builddir}/atmos_shared/tracer_driver/mo_setsox.lo: ${top_srcdir}/atmos_shared/tracer_driver/mo_setsox.F90 \
		 ${top_srcdir}/atmos_shared/tracer_driver/tropchem/mo_chem_utls.lo \
		 ${top_srcdir}/atmos_shared/tracer_driver/tropchem/tropchem_types_mod.lo \
		 ${top_srcdir}/atmos_shared/tracer_driver/cloud_chem.lo \
		 ${top_srcdir}/atmos_shared/tracer_driver/tropchem/aerosol_thermodynamics.lo
${top_builddir}/atmos_shared/tracer_driver/atmos_ch3i.lo: ${top_srcdir}/atmos_shared/tracer_driver/atmos_ch3i.F90
${top_builddir}/atmos_shared/tracer_driver/xactive_bvoc.lo: ${top_srcdir}/atmos_shared/tracer_driver/xactive_bvoc.F90 \
		 ${top_srcdir}/atmos_shared/tracer_driver/tropchem/AM4_chem/moz.mods.lo \
		 ${top_srcdir}/atmos_shared/tracer_driver/tropchem/m_tracname.lo
${top_builddir}/atmos_shared/tracer_driver/atmos_carbon_aerosol.lo: ${top_srcdir}/atmos_shared/tracer_driver/atmos_carbon_aerosol.F90 \
		 ${top_srcdir}/atmos_shared/atmos_cmip_diag/atmos_cmip_diag.lo \
		 ${top_srcdir}/atmos_shared/tracer_driver/atmos_tracer_utilities.lo
${top_builddir}/atmos_shared/tracer_driver/stratchem/strat_chem_model.lo: ${top_srcdir}/atmos_shared/tracer_driver/stratchem/strat_chem_model.F90
${top_builddir}/atmos_shared/tracer_driver/stratchem/strat_chem_driver.lo: ${top_srcdir}/atmos_shared/tracer_driver/stratchem/strat_chem_driver.F90 \
		 ${top_srcdir}/atmos_shared/tracer_driver/stratchem/strat_chem_model.lo
${top_builddir}/atmos_shared/tracer_driver/tropchem/mo_chemini.lo: ${top_srcdir}/atmos_shared/tracer_driver/tropchem/mo_chemini.F90 \
		 ${top_srcdir}/atmos_shared/tracer_driver/tropchem/mo_photo.lo \
		 ${top_srcdir}/atmos_shared/tracer_driver/tropchem/mo_chem_utls.lo \
		 ${top_srcdir}/atmos_shared/tracer_driver/tropchem/mo_usrrxt.lo \
		 ${top_srcdir}/atmos_shared/tracer_driver/tropchem/AM4_chem/moz.mods.lo \
		 ${top_srcdir}/atmos_shared/tracer_driver/tropchem/mo_exp_slv.lo \
		 ${top_srcdir}/atmos_shared/tracer_driver/tropchem/mo_imp_slv.lo \
		 ${top_srcdir}/atmos_shared/tracer_driver/tropchem/mo_rodas_slv.lo \
		 ${top_srcdir}/atmos_shared/tracer_driver/tropchem/mo_read_sim_chm.lo \
		 ${top_srcdir}/atmos_shared/tracer_driver/tropchem/AM4_chem/mo_fphoto.lo \
		 ${top_srcdir}/atmos_shared/tracer_driver/tropchem/tropchem_types_mod.lo
${top_builddir}/atmos_shared/tracer_driver/tropchem/mo_jpl.lo: ${top_srcdir}/atmos_shared/tracer_driver/tropchem/mo_jpl.F90
${top_builddir}/atmos_shared/tracer_driver/tropchem/mo_exp_slv.lo: ${top_srcdir}/atmos_shared/tracer_driver/tropchem/mo_exp_slv.F90 \
		 ${top_srcdir}/atmos_shared/tracer_driver/tropchem/AM4_chem/moz.mods.lo \
		 ${top_srcdir}/atmos_shared/tracer_driver/tropchem/mo_chem_utls.lo \
		 ${top_srcdir}/atmos_shared/tracer_driver/tropchem/AM4_chem/moz.mat.lo
${top_builddir}/atmos_shared/tracer_driver/tropchem/mo_chem_utls.lo: ${top_srcdir}/atmos_shared/tracer_driver/tropchem/mo_chem_utls.F90 \
		 ${top_srcdir}/atmos_shared/tracer_driver/tropchem/AM4_chem/moz.mods.lo \
		 ${top_srcdir}/atmos_shared/tracer_driver/tropchem/m_tracname.lo
${top_builddir}/atmos_shared/tracer_driver/tropchem/tropchem_types_mod.lo: ${top_srcdir}/atmos_shared/tracer_driver/tropchem/tropchem_types_mod.F90
${top_builddir}/atmos_shared/tracer_driver/tropchem/mo_usrrxt.lo: ${top_srcdir}/atmos_shared/tracer_driver/tropchem/mo_usrrxt.F90 \
		 ${top_srcdir}/atmos_shared/tracer_driver/tropchem/aerosol_thermodynamics.lo \
		 ${top_srcdir}/atmos_shared/tracer_driver/tropchem/strat_chem_utilities.lo \
		 ${top_srcdir}/atmos_shared/tracer_driver/tropchem/tropchem_types_mod.lo \
		 ${top_srcdir}/atmos_shared/tracer_driver/tropchem/mo_chem_utls.lo \
		 ${top_srcdir}/atmos_shared/tracer_driver/tropchem/AM4_chem/moz.mods.lo
${top_builddir}/atmos_shared/tracer_driver/tropchem/mo_photo.lo: ${top_srcdir}/atmos_shared/tracer_driver/tropchem/mo_photo.F90 \
		 ${top_srcdir}/atmos_shared/tracer_driver/tropchem/mo_chem_utls.lo \
		 ${top_srcdir}/atmos_shared/tracer_driver/tropchem/AM4_chem/moz.mods.lo
${top_builddir}/atmos_shared/tracer_driver/tropchem/aerosol_thermodynamics.lo: ${top_srcdir}/atmos_shared/tracer_driver/tropchem/aerosol_thermodynamics.F90 \
		 ${top_srcdir}/atmos_shared/tracer_driver/tropchem/tropchem_types_mod.lo
${top_builddir}/atmos_shared/tracer_driver/tropchem/mo_chemdr.lo: ${top_srcdir}/atmos_shared/tracer_driver/tropchem/mo_chemdr.F90 \
		 ${top_srcdir}/atmos_shared/tracer_driver/tropchem/tropchem_types_mod.lo \
		 ${top_srcdir}/atmos_shared/tracer_driver/tropchem/mo_chem_utls.lo \
		 ${top_srcdir}/atmos_shared/tracer_driver/tropchem/AM4_chem/moz.mods.lo \
		 ${top_srcdir}/atmos_shared/tracer_driver/tropchem/mo_photo.lo \
		 ${top_srcdir}/atmos_shared/tracer_driver/tropchem/mo_exp_slv.lo \
		 ${top_srcdir}/atmos_shared/tracer_driver/tropchem/mo_imp_slv.lo \
		 ${top_srcdir}/atmos_shared/tracer_driver/tropchem/mo_rodas_slv.lo \
		 ${top_srcdir}/atmos_shared/tracer_driver/tropchem/mo_usrrxt.lo \
		 ${top_srcdir}/atmos_shared/tracer_driver/tropchem/mo_setinv.lo \
		 ${top_srcdir}/atmos_shared/tracer_driver/tropchem/AM4_chem/moz.subs.lo \
		 ${top_srcdir}/atmos_shared/tracer_driver/mo_setsox.lo \
		 ${top_srcdir}/atmos_shared/tracer_driver/tropchem/AM4_chem/mo_fphoto.lo \
		 ${top_srcdir}/atmos_shared/tracer_driver/tropchem/strat_chem_utilities.lo
${top_builddir}/atmos_shared/tracer_driver/tropchem/strat_chem_utilities.lo: ${top_srcdir}/atmos_shared/tracer_driver/tropchem/strat_chem_utilities.F90
${top_builddir}/atmos_shared/tracer_driver/tropchem/mo_imp_slv.lo: ${top_srcdir}/atmos_shared/tracer_driver/tropchem/mo_imp_slv.F90 \
		 ${top_srcdir}/atmos_shared/tracer_driver/tropchem/AM4_chem/moz.mods.lo \
		 ${top_srcdir}/atmos_shared/tracer_driver/tropchem/mo_chem_utls.lo \
		 ${top_srcdir}/atmos_shared/tracer_driver/tropchem/m_tracname.lo \
		 ${top_srcdir}/atmos_shared/tracer_driver/tropchem/AM4_chem/moz.mat.lo
${top_builddir}/atmos_shared/tracer_driver/tropchem/isoropiaIIcode.lo: ${top_srcdir}/atmos_shared/tracer_driver/tropchem/isoropiaIIcode.F
${top_builddir}/atmos_shared/tracer_driver/tropchem/mo_hook.lo: ${top_srcdir}/atmos_shared/tracer_driver/tropchem/mo_hook.F90 \
		 ${top_srcdir}/atmos_shared/atmos_cmip_diag/atmos_cmip_diag.lo
${top_builddir}/atmos_shared/tracer_driver/tropchem/mo_setinv.lo: ${top_srcdir}/atmos_shared/tracer_driver/tropchem/mo_setinv.F90
${top_builddir}/atmos_shared/tracer_driver/tropchem/mo_rodas_slv.lo: ${top_srcdir}/atmos_shared/tracer_driver/tropchem/mo_rodas_slv.F90 \
		 ${top_srcdir}/atmos_shared/tracer_driver/tropchem/AM4_chem/moz.mods.lo \
		 ${top_srcdir}/atmos_shared/tracer_driver/tropchem/mo_chem_utls.lo \
		 ${top_srcdir}/atmos_shared/tracer_driver/tropchem/AM4_chem/moz.mat.lo
${top_builddir}/atmos_shared/tracer_driver/tropchem/mo_read_sim_chm.lo: ${top_srcdir}/atmos_shared/tracer_driver/tropchem/mo_read_sim_chm.F90 \
		 ${top_srcdir}/atmos_shared/tracer_driver/tropchem/AM4_chem/moz.mods.lo \
		 ${top_srcdir}/atmos_shared/tracer_driver/tropchem/m_tracname.lo
${top_builddir}/atmos_shared/tracer_driver/tropchem/m_tracname.lo: ${top_srcdir}/atmos_shared/tracer_driver/tropchem/m_tracname.F90 \
		 ${top_srcdir}/atmos_shared/tracer_driver/tropchem/AM4_chem/moz.mods.lo
${top_builddir}/atmos_shared/tracer_driver/tropchem/AM4_chem/moz.mods.lo: ${top_srcdir}/atmos_shared/tracer_driver/tropchem/AM4_chem/moz.mods.F90
${top_builddir}/atmos_shared/tracer_driver/tropchem/AM4_chem/mo_fphoto.lo: ${top_srcdir}/atmos_shared/tracer_driver/tropchem/AM4_chem/mo_fphoto.F90 \
		 ${top_srcdir}/atmos_shared/tracer_driver/tropchem/AM4_chem/mo_fastjx.lo \
		 ${top_srcdir}/atmos_shared/tracer_driver/tropchem/mo_chem_utls.lo \
		 ${top_srcdir}/atmos_shared/tracer_driver/tropchem/AM4_chem/moz.mods.lo
${top_builddir}/atmos_shared/tracer_driver/tropchem/AM4_chem/moz.mat.lo: ${top_srcdir}/atmos_shared/tracer_driver/tropchem/AM4_chem/moz.mat.F90 \
		 ${top_srcdir}/atmos_shared/tracer_driver/tropchem/AM4_chem/moz.mods.lo
${top_builddir}/atmos_shared/tracer_driver/tropchem/AM4_chem/moz.subs.lo: ${top_srcdir}/atmos_shared/tracer_driver/tropchem/AM4_chem/moz.subs.F90 \
		 ${top_srcdir}/atmos_shared/tracer_driver/tropchem/AM4_chem/moz.mods.lo \
		 ${top_srcdir}/atmos_shared/tracer_driver/tropchem/mo_jpl.lo
${top_builddir}/atmos_shared/tracer_driver/tropchem/AM4_chem/mo_fastjx.lo: ${top_srcdir}/atmos_shared/tracer_driver/tropchem/AM4_chem/mo_fastjx.F90 \
		 ${top_srcdir}/atmos_shared/tracer_driver/tropchem/tropchem_types_mod.lo
${top_builddir}/atmos_shared/tracer_driver/aer_ccn_act/aer_ccn_act.lo: ${top_srcdir}/atmos_shared/tracer_driver/aer_ccn_act/aer_ccn_act.F90 \
		 ${top_srcdir}/atmos_shared/tracer_driver/aer_ccn_act/aer_ccn_act_k.lo
${top_builddir}/atmos_shared/tracer_driver/aer_ccn_act/aer_in_act.lo: ${top_srcdir}/atmos_shared/tracer_driver/aer_ccn_act/aer_in_act.F90
${top_builddir}/atmos_shared/tracer_driver/aer_ccn_act/aer_ccn_act_k.lo: ${top_srcdir}/atmos_shared/tracer_driver/aer_ccn_act/aer_ccn_act_k.F90
${top_builddir}/atmos_shared/tracer_driver/aer_ccn_act/ice_nucl.lo: ${top_srcdir}/atmos_shared/tracer_driver/aer_ccn_act/ice_nucl.F90 \
		 ${top_srcdir}/atmos_shared/tracer_driver/aer_ccn_act/aer_ccn_act_k.lo \
		 ${top_srcdir}/atmos_shared/tracer_driver/aer_ccn_act/aerosol_params.lo
${top_builddir}/atmos_shared/tracer_driver/aer_ccn_act/aerosol_params.lo: ${top_srcdir}/atmos_shared/tracer_driver/aer_ccn_act/aerosol_params.F90
${top_builddir}/atmos_shared/tracer_driver/atmos_ozone_tracer.lo: ${top_srcdir}/atmos_shared/tracer_driver/atmos_ozone_tracer.F90
${top_builddir}/atmos_param/physics_driver/physics_types.lo: ${top_srcdir}/atmos_param/physics_driver/physics_types.F90
${top_builddir}/atmos_param/physics_driver/physics_driver.lo: ${top_srcdir}/atmos_param/physics_driver/physics_driver.F90 \
		 ${top_srcdir}/atmos_shared/tracer_driver/atmos_tracer_driver.lo \
		 ${top_srcdir}/atmos_shared/atmos_cmip_diag/atmos_cmip_diag.lo \
		 ${top_srcdir}/atmos_shared/aerosol/aerosol_types.lo \
		 ${top_srcdir}/atmos_param/physics_radiation_exch/physics_radiation_exch.lo \
		 ${top_srcdir}/atmos_param/physics_driver/physics_types.lo \
		 ${top_srcdir}/atmos_param/moist_processes/moist_processes_utils.lo \
		 ${top_srcdir}/atmos_shared/aerosol/aerosol.lo \
		 ${top_srcdir}/atmos_param/cosp/cosp_driver.lo \
		 ${top_srcdir}/atmos_param/moist_processes/moist_processes.lo \
		 ${top_srcdir}/atmos_param/vert_turb_driver/vert_turb_driver.lo \
		 ${top_srcdir}/atmos_param/vert_diff_driver/vert_diff_driver.lo \
		 ${top_srcdir}/atmos_param/damping_driver/damping_driver.lo
${top_builddir}/atmos_param/topo_drag/topo_drag.lo: ${top_srcdir}/atmos_param/topo_drag/topo_drag.F90
${top_builddir}/atmos_param/macrophysics/ls_cloud_macrophysics.lo: ${top_srcdir}/atmos_param/macrophysics/ls_cloud_macrophysics.F90 \
		 ${top_srcdir}/atmos_param/physics_driver/physics_types.lo \
		 ${top_srcdir}/atmos_param/macrophysics/tiedtke_macro.lo \
		 ${top_srcdir}/atmos_param/aerosol_cloud/aerosol_cloud.lo \
		 ${top_srcdir}/atmos_param/lscloud_driver/lscloud_types.lo \
		 ${top_srcdir}/atmos_shared/aerosol/aerosol_types.lo \
		 ${top_srcdir}/atmos_param/moist_processes/moist_processes_utils.lo \
		 ${top_srcdir}/atmos_param/physics_radiation_exch/physics_radiation_exch.lo
${top_builddir}/atmos_param/macrophysics/tiedtke_macro.lo: ${top_srcdir}/atmos_param/macrophysics/tiedtke_macro.F90 \
		 ${top_srcdir}/atmos_param/cloud_generator/betaDistribution.lo \
		 ${top_srcdir}/atmos_param/lscloud_driver/lscloud_types.lo \
		 ${top_srcdir}/atmos_param/lscloud_driver/lscloud_debug.lo \
		 ${top_srcdir}/atmos_param/lscloud_driver/polysvp.lo \
		 ${top_srcdir}/atmos_param/moist_processes/moist_processes_utils.lo \
		 ${top_srcdir}/atmos_param/physics_radiation_exch/physics_radiation_exch.lo
${top_builddir}/atmos_param/radiation/clouds/cloudrad_types.lo: ${top_srcdir}/atmos_param/radiation/clouds/cloudrad_types.F90
${top_builddir}/atmos_param/radiation/clouds/microphys_rad.lo: ${top_srcdir}/atmos_param/radiation/clouds/microphys_rad.F90 \
		 ${top_srcdir}/atmos_param/radiation/clouds/cloudrad_types.lo \
		 ${top_srcdir}/atmos_param/radiation/longwave/sealw99.lo \
		 ${top_srcdir}/atmos_param/radiation/shortwave/esfsw_driver.lo
${top_builddir}/atmos_param/radiation/clouds/strat_clouds_W.lo: ${top_srcdir}/atmos_param/radiation/clouds/strat_clouds_W.F90 \
		 ${top_srcdir}/atmos_shared/aerosol/aerosol_types.lo \
		 ${top_srcdir}/atmos_param/physics_radiation_exch/physics_radiation_exch.lo \
		 ${top_srcdir}/atmos_param/cloud_rad/cloud_rad.lo \
		 ${top_srcdir}/atmos_param/radiation/clouds/cloudrad_types.lo \
		 ${top_srcdir}/atmos_param/radiation/clouds/get_random_number_stream.lo \
		 ${top_srcdir}/atmos_param/cloud_generator/cloud_generator.lo
${top_builddir}/atmos_param/radiation/clouds/isccp_clouds.lo: ${top_srcdir}/atmos_param/radiation/clouds/isccp_clouds.F90
${top_builddir}/atmos_param/radiation/clouds/get_random_number_stream.lo: ${top_srcdir}/atmos_param/radiation/clouds/get_random_number_stream.F90 \
		 ${top_srcdir}/atmos_param/radiation/clouds/cloudrad_types.lo
${top_builddir}/atmos_param/radiation/clouds/uw_clouds_W.lo: ${top_srcdir}/atmos_param/radiation/clouds/uw_clouds_W.F90 \
		 ${top_srcdir}/atmos_param/physics_radiation_exch/physics_radiation_exch.lo \
		 ${top_srcdir}/atmos_param/radiation/clouds/cloudrad_types.lo
${top_builddir}/atmos_param/radiation/clouds/cloudrad_driver.lo: ${top_srcdir}/atmos_param/radiation/clouds/cloudrad_driver.F90 \
		 ${top_srcdir}/atmos_shared/aerosol/aerosol_types.lo \
		 ${top_srcdir}/atmos_param/physics_radiation_exch/physics_radiation_exch.lo \
		 ${top_srcdir}/atmos_param/radiation/clouds/cloudrad_types.lo \
		 ${top_srcdir}/atmos_param/radiation/clouds/cloud_spec.lo \
		 ${top_srcdir}/atmos_param/radiation/clouds/cloudrad_package.lo
${top_builddir}/atmos_param/radiation/clouds/cloud_spec.lo: ${top_srcdir}/atmos_param/radiation/clouds/cloud_spec.F90 \
		 ${top_srcdir}/atmos_param/radiation/clouds/get_random_number_stream.lo \
		 ${top_srcdir}/atmos_shared/aerosol/aerosol_types.lo \
		 ${top_srcdir}/atmos_param/physics_radiation_exch/physics_radiation_exch.lo \
		 ${top_srcdir}/atmos_param/radiation/clouds/cloudrad_types.lo \
		 ${top_srcdir}/atmos_param/radiation/clouds/strat_clouds_W.lo \
		 ${top_srcdir}/atmos_param/radiation/clouds/uw_clouds_W.lo
${top_builddir}/atmos_param/radiation/clouds/cloudrad_diagnostics.lo: ${top_srcdir}/atmos_param/radiation/clouds/cloudrad_diagnostics.F90 \
		 ${top_srcdir}/atmos_shared/atmos_cmip_diag/atmos_cmip_diag.lo \
		 ${top_srcdir}/atmos_param/radiation/clouds/cloudrad_types.lo \
		 ${top_srcdir}/atmos_param/radiation/clouds/microphys_rad.lo \
		 ${top_srcdir}/atmos_param/radiation/clouds/isccp_clouds.lo
${top_builddir}/atmos_param/radiation/clouds/cloudrad_package.lo: ${top_srcdir}/atmos_param/radiation/clouds/cloudrad_package.F90 \
		 ${top_srcdir}/atmos_param/radiation/clouds/cloudrad_types.lo \
		 ${top_srcdir}/atmos_param/radiation/clouds/cloudrad_diagnostics.lo \
		 ${top_srcdir}/atmos_param/radiation/clouds/microphys_rad.lo
${top_builddir}/atmos_param/radiation/longwave/longwave_tables.lo: ${top_srcdir}/atmos_param/radiation/longwave/longwave_tables.F90 \
		 ${top_srcdir}/atmos_param/radiation/longwave/longwave_utilities.lo \
		 ${top_srcdir}/atmos_param/radiation/longwave/longwave_params.lo
${top_builddir}/atmos_param/radiation/longwave/longwave_clouds.lo: ${top_srcdir}/atmos_param/radiation/longwave/longwave_clouds.F90 \
		 ${top_srcdir}/atmos_param/radiation/longwave/longwave_utilities.lo \
		 ${top_srcdir}/atmos_param/radiation/longwave/longwave_types.lo
${top_builddir}/atmos_param/radiation/longwave/gas_tf.lo: ${top_srcdir}/atmos_param/radiation/longwave/gas_tf.F90 \
		 ${top_srcdir}/atmos_param/radiation/longwave/longwave_utilities.lo \
		 ${top_srcdir}/atmos_param/radiation/longwave/longwave_params.lo
${top_builddir}/atmos_param/radiation/longwave/longwave_params.lo: ${top_srcdir}/atmos_param/radiation/longwave/longwave_params.F90
${top_builddir}/atmos_param/radiation/longwave/sealw99.lo: ${top_srcdir}/atmos_param/radiation/longwave/sealw99.F90 \
		 ${top_srcdir}/atmos_param/radiation/longwave/longwave_utilities.lo \
		 ${top_srcdir}/atmos_param/radiation/longwave/longwave_types.lo \
		 ${top_srcdir}/atmos_param/radiation/longwave/longwave_params.lo \
		 ${top_srcdir}/atmos_param/radiation/longwave/longwave_clouds.lo \
		 ${top_srcdir}/atmos_param/radiation/longwave/longwave_fluxes.lo \
		 ${top_srcdir}/atmos_param/radiation/longwave/longwave_tables.lo \
		 ${top_srcdir}/atmos_param/radiation/longwave/optical_path.lo \
		 ${top_srcdir}/atmos_param/radiation/longwave/gas_tf.lo \
		 ${top_srcdir}/atmos_param/radiation/longwave/lw_gases_stdtf.lo
${top_builddir}/atmos_param/radiation/longwave/optical_path.lo: ${top_srcdir}/atmos_param/radiation/longwave/optical_path.F90 \
		 ${top_srcdir}/atmos_param/radiation/longwave/longwave_utilities.lo \
		 ${top_srcdir}/atmos_param/radiation/longwave/longwave_params.lo \
		 ${top_srcdir}/atmos_param/radiation/longwave/lw_gases_stdtf.lo
${top_builddir}/atmos_param/radiation/longwave/longwave_types.lo: ${top_srcdir}/atmos_param/radiation/longwave/longwave_types.F90
${top_builddir}/atmos_param/radiation/longwave/longwave_utilities.lo: ${top_srcdir}/atmos_param/radiation/longwave/longwave_utilities.F90
${top_builddir}/atmos_param/radiation/longwave/longwave_fluxes.lo: ${top_srcdir}/atmos_param/radiation/longwave/longwave_fluxes.F90 \
		 ${top_srcdir}/atmos_param/radiation/longwave/longwave_types.lo
${top_builddir}/atmos_param/radiation/longwave/lw_gases_stdtf.lo: ${top_srcdir}/atmos_param/radiation/longwave/lw_gases_stdtf.F90 \
		 ${top_srcdir}/atmos_param/radiation/longwave/longwave_utilities.lo \
		 ${top_srcdir}/atmos_param/radiation/longwave/gas_tf.lo
${top_builddir}/atmos_param/radiation/aerosols/aerosolrad_types.lo: ${top_srcdir}/atmos_param/radiation/aerosols/aerosolrad_types.F90
${top_builddir}/atmos_param/radiation/aerosols/aerosolrad_package.lo: ${top_srcdir}/atmos_param/radiation/aerosols/aerosolrad_package.F90 \
		 ${top_srcdir}/atmos_shared/aerosol/aerosol_types.lo \
		 ${top_srcdir}/atmos_param/radiation/aerosols/aerosolrad_types.lo \
		 ${top_srcdir}/atmos_param/radiation/longwave/sealw99.lo \
		 ${top_srcdir}/atmos_param/radiation/shortwave/esfsw_driver.lo
${top_builddir}/atmos_param/radiation/aerosols/aerosolrad_driver.lo: ${top_srcdir}/atmos_param/radiation/aerosols/aerosolrad_driver.F90 \
		 ${top_srcdir}/atmos_shared/aerosol/aerosol_types.lo \
		 ${top_srcdir}/atmos_shared/aerosol/aerosol.lo \
		 ${top_srcdir}/atmos_param/radiation/aerosols/aerosolrad_types.lo \
		 ${top_srcdir}/atmos_param/radiation/aerosols/aerosolrad_package.lo
${top_builddir}/atmos_param/radiation/shortwave/esfsw_driver.lo: ${top_srcdir}/atmos_param/radiation/shortwave/esfsw_driver.F90 \
		 ${top_srcdir}/atmos_param/radiation/shortwave/esfsw_parameters.lo \
		 ${top_srcdir}/atmos_param/radiation/shortwave/esfsw_bands.lo \
		 ${top_srcdir}/atmos_param/radiation/shortwave/esfsw_utilities.lo \
		 ${top_srcdir}/atmos_param/radiation/shortwave/shortwave_types.lo
${top_builddir}/atmos_param/radiation/shortwave/esfsw_utilities.lo: ${top_srcdir}/atmos_param/radiation/shortwave/esfsw_utilities.F90
${top_builddir}/atmos_param/radiation/shortwave/shortwave_types.lo: ${top_srcdir}/atmos_param/radiation/shortwave/shortwave_types.F90
${top_builddir}/atmos_param/radiation/shortwave/esfsw_bands.lo: ${top_srcdir}/atmos_param/radiation/shortwave/esfsw_bands.F90 \
		 ${top_srcdir}/atmos_param/radiation/shortwave/esfsw_utilities.lo
${top_builddir}/atmos_param/radiation/shortwave/esfsw_parameters.lo: ${top_srcdir}/atmos_param/radiation/shortwave/esfsw_parameters.F90
${top_builddir}/atmos_param/radiation/radiative_gases/radiative_gases.lo: ${top_srcdir}/atmos_param/radiation/radiative_gases/radiative_gases.F90 \
		 ${top_srcdir}/atmos_param/radiation/longwave/sealw99.lo \
		 ${top_srcdir}/atmos_param/radiation/radiative_gases/radiative_gases_types.lo \
		 ${top_srcdir}/atmos_param/radiation/radiative_gases/ozone.lo
${top_builddir}/atmos_param/radiation/radiative_gases/radiative_gases_types.lo: ${top_srcdir}/atmos_param/radiation/radiative_gases/radiative_gases_types.F90
${top_builddir}/atmos_param/radiation/radiative_gases/ozone.lo: ${top_srcdir}/atmos_param/radiation/radiative_gases/ozone.F90 \
		 ${top_srcdir}/atmos_param/radiation/radiative_gases/radiative_gases_types.lo
${top_builddir}/atmos_param/radiation/util/radiation_diag.lo: ${top_srcdir}/atmos_param/radiation/util/radiation_diag.F90 \
		 ${top_srcdir}/atmos_param/radiation/driver/radiation_driver.lo
${top_builddir}/atmos_param/radiation/driver/longwave_driver.lo: ${top_srcdir}/atmos_param/radiation/driver/longwave_driver.F90 \
		 ${top_srcdir}/atmos_param/radiation/driver/radiation_driver_types.lo \
		 ${top_srcdir}/atmos_param/radiation/aerosols/aerosolrad_types.lo \
		 ${top_srcdir}/atmos_param/radiation/longwave/sealw99.lo \
		 ${top_srcdir}/atmos_param/radiation/longwave/longwave_types.lo \
		 ${top_srcdir}/atmos_param/radiation/radiative_gases/radiative_gases.lo \
		 ${top_srcdir}/atmos_param/radiation/radiative_gases/radiative_gases_types.lo
${top_builddir}/atmos_param/radiation/driver/rad_output_file.lo: ${top_srcdir}/atmos_param/radiation/driver/rad_output_file.F90 \
		 ${top_srcdir}/atmos_shared/aerosol/aerosol_types.lo \
		 ${top_srcdir}/atmos_shared/atmos_cmip_diag/atmos_cmip_diag.lo \
		 ${top_srcdir}/atmos_param/radiation/driver/radiation_driver_types.lo \
		 ${top_srcdir}/atmos_param/radiation/aerosols/aerosolrad_types.lo \
		 ${top_srcdir}/atmos_param/radiation/clouds/cloudrad_types.lo \
		 ${top_srcdir}/atmos_param/radiation/shortwave/esfsw_driver.lo \
		 ${top_srcdir}/atmos_param/radiation/shortwave/shortwave_types.lo \
		 ${top_srcdir}/atmos_param/radiation/longwave/longwave_types.lo
${top_builddir}/atmos_param/radiation/driver/radiation_driver_types.lo: ${top_srcdir}/atmos_param/radiation/driver/radiation_driver_types.F90
${top_builddir}/atmos_param/radiation/driver/radiation_types.lo: ${top_srcdir}/atmos_param/radiation/driver/radiation_types.F90 \
		 ${top_srcdir}/atmos_shared/tracer_driver/atmos_co2.lo \
		 ${top_srcdir}/atmos_shared/tracer_driver/atmos_ch4.lo
${top_builddir}/atmos_param/radiation/driver/solar_data_driver.lo: ${top_srcdir}/atmos_param/radiation/driver/solar_data_driver.F90
${top_builddir}/atmos_param/radiation/driver/shortwave_driver.lo: ${top_srcdir}/atmos_param/radiation/driver/shortwave_driver.F90 \
		 ${top_srcdir}/atmos_param/radiation/driver/radiation_driver_types.lo \
		 ${top_srcdir}/atmos_param/radiation/shortwave/esfsw_driver.lo \
		 ${top_srcdir}/atmos_param/radiation/shortwave/shortwave_types.lo \
		 ${top_srcdir}/atmos_param/radiation/radiative_gases/radiative_gases_types.lo \
		 ${top_srcdir}/atmos_param/radiation/driver/solar_data_driver.lo
${top_builddir}/atmos_param/radiation/driver/radiation_driver_diag.lo: ${top_srcdir}/atmos_param/radiation/driver/radiation_driver_diag.F90 \
		 ${top_srcdir}/atmos_shared/atmos_cmip_diag/atmos_cmip_diag.lo \
		 ${top_srcdir}/atmos_shared/atmos_cmip_diag/atmos_global_diag.lo \
		 ${top_srcdir}/atmos_param/radiation/driver/radiation_driver_types.lo \
		 ${top_srcdir}/atmos_param/radiation/aerosols/aerosolrad_types.lo \
		 ${top_srcdir}/atmos_param/radiation/shortwave/shortwave_types.lo \
		 ${top_srcdir}/atmos_param/radiation/longwave/longwave_types.lo \
		 ${top_srcdir}/atmos_param/radiation/radiative_gases/radiative_gases_types.lo \
		 ${top_srcdir}/atmos_param/radiation/driver/shortwave_driver.lo
${top_builddir}/atmos_param/radiation/driver/radiation_driver.lo: ${top_srcdir}/atmos_param/radiation/driver/radiation_driver.F90 \
		 ${top_srcdir}/atmos_param/radiation/driver/radiation_types.lo \
		 ${top_srcdir}/atmos_param/physics_radiation_exch/physics_radiation_exch.lo \
		 ${top_srcdir}/atmos_param/radiation/driver/radiation_driver_types.lo \
		 ${top_srcdir}/atmos_shared/aerosol/aerosol_types.lo \
		 ${top_srcdir}/atmos_param/radiation/aerosols/aerosolrad_types.lo \
		 ${top_srcdir}/atmos_param/radiation/clouds/cloudrad_types.lo \
		 ${top_srcdir}/atmos_param/radiation/driver/longwave_driver.lo \
		 ${top_srcdir}/atmos_param/radiation/longwave/longwave_types.lo \
		 ${top_srcdir}/atmos_param/radiation/driver/shortwave_driver.lo \
		 ${top_srcdir}/atmos_param/radiation/shortwave/shortwave_types.lo \
		 ${top_srcdir}/atmos_param/radiation/driver/rad_output_file.lo \
		 ${top_srcdir}/atmos_param/radiation/aerosols/aerosolrad_driver.lo \
		 ${top_srcdir}/atmos_param/radiation/clouds/cloudrad_driver.lo \
		 ${top_srcdir}/atmos_param/radiation/clouds/cloudrad_diagnostics.lo \
		 ${top_srcdir}/atmos_param/radiation/radiative_gases/radiative_gases.lo \
		 ${top_srcdir}/atmos_param/radiation/radiative_gases/radiative_gases_types.lo \
		 ${top_srcdir}/atmos_param/radiation/driver/radiation_driver_diag.lo
${top_builddir}/atmos_param/convection_driver/convection_utilities.lo: ${top_srcdir}/atmos_param/convection_driver/convection_utilities.F90
${top_builddir}/atmos_param/convection_driver/detr_ice_num.lo: ${top_srcdir}/atmos_param/convection_driver/detr_ice_num.F90 \
		 ${top_srcdir}/atmos_param/lscloud_driver/lscloud_constants.lo
${top_builddir}/atmos_param/convection_driver/convection_driver.lo: ${top_srcdir}/atmos_param/convection_driver/convection_driver.F90 \
		 ${top_srcdir}/atmos_shared/atmos_cmip_diag/atmos_global_diag.lo \
		 ${top_srcdir}/atmos_shared/atmos_cmip_diag/atmos_cmip_diag.lo \
		 ${top_srcdir}/atmos_param/physics_driver/physics_types.lo \
		 ${top_srcdir}/atmos_param/vert_diff_driver/vert_diff_driver.lo \
		 ${top_srcdir}/atmos_param/physics_radiation_exch/physics_radiation_exch.lo \
		 ${top_srcdir}/atmos_param/shallow_cu/uw_conv.lo \
		 ${top_srcdir}/atmos_param/convection_driver/detr_ice_num.lo \
		 ${top_srcdir}/atmos_param/cu_mo_trans/cu_mo_trans.lo \
		 ${top_srcdir}/atmos_shared/tracer_driver/tropchem/mo_hook.lo \
		 ${top_srcdir}/atmos_shared/aerosol/aerosol_types.lo \
		 ${top_srcdir}/atmos_param/moist_processes/moist_processes_utils.lo \
		 ${top_srcdir}/atmos_param/convection_driver/convection_utilities.lo \
		 ${top_srcdir}/atmos_shared/tracer_driver/atmos_tracer_utilities.lo
${top_builddir}/atmos_param/physics_radiation_exch/physics_radiation_exch.lo: ${top_srcdir}/atmos_param/physics_radiation_exch/physics_radiation_exch.F90
${top_builddir}/atmos_param/cloud_generator/betaDistribution.lo: ${top_srcdir}/atmos_param/cloud_generator/betaDistribution.F90
${top_builddir}/atmos_param/cloud_generator/cloud_generator.lo: ${top_srcdir}/atmos_param/cloud_generator/cloud_generator.F90 \
		 ${top_srcdir}/atmos_param/cloud_generator/betaDistribution.lo
${top_builddir}/atmos_param/lscloud_driver/lscloud_debug.lo: ${top_srcdir}/atmos_param/lscloud_driver/lscloud_debug.F90 \
		 ${top_srcdir}/atmos_param/lscloud_driver/lscloud_types.lo \
		 ${top_srcdir}/atmos_param/moist_processes/moist_processes_utils.lo \
		 ${top_srcdir}/atmos_param/lscloud_driver/check_nan.lo
${top_builddir}/atmos_param/lscloud_driver/lscloud_constants.lo: ${top_srcdir}/atmos_param/lscloud_driver/lscloud_constants.F90
${top_builddir}/atmos_param/lscloud_driver/lscloud_netcdf.lo: ${top_srcdir}/atmos_param/lscloud_driver/lscloud_netcdf.F90 \
		 ${top_srcdir}/atmos_param/lscloud_driver/lscloud_types.lo \
		 ${top_srcdir}/atmos_shared/atmos_cmip_diag/atmos_cmip_diag.lo
${top_builddir}/atmos_param/lscloud_driver/lscloud_types.lo: ${top_srcdir}/atmos_param/lscloud_driver/lscloud_types.F90 \
		 ${top_srcdir}/atmos_shared/atmos_cmip_diag/atmos_cmip_diag.lo
${top_builddir}/atmos_param/lscloud_driver/lscloud_driver.lo: ${top_srcdir}/atmos_param/lscloud_driver/lscloud_driver.F90 \
		 ${top_srcdir}/atmos_param/physics_driver/physics_types.lo \
		 ${top_srcdir}/atmos_param/lscale_cond/lscale_cond.lo \
		 ${top_srcdir}/atmos_param/lscloud_driver/lscloud_debug.lo \
		 ${top_srcdir}/atmos_param/lscloud_driver/lscloud_types.lo \
		 ${top_srcdir}/atmos_param/lscloud_driver/polysvp.lo \
		 ${top_srcdir}/atmos_param/lscloud_driver/lscloud_netcdf.lo \
		 ${top_srcdir}/atmos_param/physics_radiation_exch/physics_radiation_exch.lo \
		 ${top_srcdir}/atmos_param/macrophysics/ls_cloud_macrophysics.lo \
		 ${top_srcdir}/atmos_param/microphysics/ls_cloud_microphysics.lo \
		 ${top_srcdir}/atmos_param/aerosol_cloud/aerosol_cloud.lo \
		 ${top_srcdir}/atmos_param/moist_processes/moist_processes_utils.lo \
		 ${top_srcdir}/atmos_shared/aerosol/aerosol_types.lo \
		 ${top_srcdir}/atmos_shared/tracer_driver/atmos_tracer_utilities.lo \
		 ${top_srcdir}/atmos_shared/atmos_cmip_diag/atmos_global_diag.lo \
		 ${top_srcdir}/atmos_shared/atmos_cmip_diag/atmos_cmip_diag.lo
${top_builddir}/atmos_param/lscloud_driver/polysvp.lo: ${top_srcdir}/atmos_param/lscloud_driver/polysvp.F90 \
		 ${top_srcdir}/atmos_param/lscloud_driver/lscloud_constants.lo \
		 ${top_srcdir}/atmos_param/lscloud_driver/lscloud_types.lo
${top_builddir}/atmos_param/lscloud_driver/check_nan.lo: ${top_srcdir}/atmos_param/lscloud_driver/check_nan.F90
${top_builddir}/atmos_param/stable_bl_turb/stable_bl_turb.lo: ${top_srcdir}/atmos_param/stable_bl_turb/stable_bl_turb.F90
${top_builddir}/atmos_param/microphysics/gamma_mg.lo: ${top_srcdir}/atmos_param/microphysics/gamma_mg.F90 \
		 ${top_srcdir}/atmos_param/lscloud_driver/lscloud_constants.lo
${top_builddir}/atmos_param/microphysics/micro_mg2.lo: ${top_srcdir}/atmos_param/microphysics/micro_mg2.F90 \
		 ${top_srcdir}/atmos_param/microphysics/micro_mg2_utils.lo \
		 ${top_srcdir}/atmos_param/microphysics/gamma_mg.lo \
		 ${top_srcdir}/atmos_param/lscloud_driver/lscloud_types.lo \
		 ${top_srcdir}/atmos_param/physics_radiation_exch/physics_radiation_exch.lo
${top_builddir}/atmos_param/microphysics/micro_mg2_utils.lo: ${top_srcdir}/atmos_param/microphysics/micro_mg2_utils.F90 \
		 ${top_srcdir}/atmos_param/microphysics/gamma_mg.lo
${top_builddir}/atmos_param/microphysics/rotstayn_klein_mp.lo: ${top_srcdir}/atmos_param/microphysics/rotstayn_klein_mp.F90 \
		 ${top_srcdir}/atmos_param/cloud_generator/cloud_generator.lo \
		 ${top_srcdir}/atmos_shared/tracer_driver/aer_ccn_act/aer_in_act.lo \
		 ${top_srcdir}/atmos_param/lscloud_driver/polysvp.lo \
		 ${top_srcdir}/atmos_param/lscloud_driver/lscloud_types.lo \
		 ${top_srcdir}/atmos_param/lscloud_driver/lscloud_debug.lo \
		 ${top_srcdir}/atmos_param/physics_radiation_exch/physics_radiation_exch.lo
${top_builddir}/atmos_param/microphysics/ls_cloud_microphysics.lo: ${top_srcdir}/atmos_param/microphysics/ls_cloud_microphysics.F90 \
		 ${top_srcdir}/atmos_param/physics_driver/physics_types.lo \
		 ${top_srcdir}/atmos_param/lscloud_driver/lscloud_types.lo \
		 ${top_srcdir}/atmos_shared/aerosol/aerosol_types.lo \
		 ${top_srcdir}/atmos_param/physics_radiation_exch/physics_radiation_exch.lo \
		 ${top_srcdir}/atmos_param/moist_processes/moist_processes_utils.lo \
		 ${top_srcdir}/atmos_param/lscloud_driver/lscloud_debug.lo \
		 ${top_srcdir}/atmos_param/microphysics/rotstayn_klein_mp.lo \
		 ${top_srcdir}/atmos_param/microphysics/micro_mg2.lo
${top_builddir}/atmos_param/damping_driver/damping_driver.lo: ${top_srcdir}/atmos_param/damping_driver/damping_driver.F90 \
		 ${top_srcdir}/atmos_param/mg_drag/mg_drag.lo \
		 ${top_srcdir}/atmos_param/cg_drag/cg_drag.lo \
		 ${top_srcdir}/atmos_param/topo_drag/topo_drag.lo \
		 ${top_srcdir}/atmos_shared/atmos_cmip_diag/atmos_cmip_diag.lo
${top_builddir}/atmos_param/vert_diff_driver/vert_diff_driver.lo: ${top_srcdir}/atmos_param/vert_diff_driver/vert_diff_driver.F90 \
		 ${top_srcdir}/atmos_param/vert_diff/vert_diff.lo \
		 ${top_srcdir}/atmos_shared/atmos_cmip_diag/atmos_cmip_diag.lo
${top_builddir}/atmos_param/moist_processes/moist_processes.lo: ${top_srcdir}/atmos_param/moist_processes/moist_processes.F90 \
		 ${top_srcdir}/atmos_param/physics_driver/physics_types.lo \
		 ${top_srcdir}/atmos_param/physics_radiation_exch/physics_radiation_exch.lo \
		 ${top_srcdir}/atmos_param/lscloud_driver/lscloud_driver.lo \
		 ${top_srcdir}/atmos_param/convection_driver/convection_driver.lo \
		 ${top_srcdir}/atmos_param/convection_driver/convection_utilities.lo \
		 ${top_srcdir}/atmos_shared/atmos_cmip_diag/atmos_global_diag.lo \
		 ${top_srcdir}/atmos_param/vert_diff_driver/vert_diff_driver.lo \
		 ${top_srcdir}/atmos_shared/aerosol/aerosol_types.lo \
		 ${top_srcdir}/atmos_shared/tracer_driver/atmos_tracer_utilities.lo \
		 ${top_srcdir}/atmos_param/moist_processes/moist_processes_utils.lo \
		 ${top_srcdir}/atmos_shared/tracer_driver/atmos_dust.lo \
		 ${top_srcdir}/atmos_shared/tracer_driver/atmos_tracer_driver.lo \
		 ${top_srcdir}/atmos_shared/tracer_driver/atmos_sea_salt.lo \
		 ${top_srcdir}/atmos_shared/atmos_cmip_diag/atmos_cmip_diag.lo
${top_builddir}/atmos_param/moist_processes/moist_processes_utils.lo: ${top_srcdir}/atmos_param/moist_processes/moist_processes_utils.F90 \
		 ${top_srcdir}/atmos_param/lscloud_driver/lscloud_types.lo
${top_builddir}/atmos_param/cosp/cosp_modis_simulator.lo: ${top_srcdir}/atmos_param/cosp/cosp_modis_simulator.F90 \
		 ${top_srcdir}/atmos_param/cosp/cosp_constants.lo \
		 ${top_srcdir}/atmos_param/cosp/cosp_types.lo \
		 ${top_srcdir}/atmos_param/cosp/MODIS_simulator/modis_simulator.lo
${top_builddir}/atmos_param/cosp/cosp_driver.lo: ${top_srcdir}/atmos_param/cosp/cosp_driver.F90 \
		 ${top_srcdir}/atmos_param/physics_radiation_exch/physics_radiation_exch.lo \
		 ${top_srcdir}/atmos_param/physics_driver/physics_types.lo \
		 ${top_srcdir}/atmos_param/cosp/cosp_types.lo \
		 ${top_srcdir}/atmos_param/cosp/cosp.lo \
		 ${top_srcdir}/atmos_param/cosp/cosp_io.lo \
		 ${top_srcdir}/atmos_param/cosp/cosp_constants.lo \
		 ${top_srcdir}/atmos_param/cosp/cosp_modis_simulator.lo \
		 ${top_srcdir}/atmos_param/cosp/cosp_diagnostics.lo \
		 ${top_srcdir}/atmos_param/cosp/cosp_utils.lo
${top_builddir}/atmos_param/cosp/cosp_misr_simulator.lo: ${top_srcdir}/atmos_param/cosp/cosp_misr_simulator.F90 \
		 ${top_srcdir}/atmos_param/cosp/cosp_constants.lo \
		 ${top_srcdir}/atmos_param/cosp/cosp_types.lo
${top_builddir}/atmos_param/cosp/cosp_stats.lo: ${top_srcdir}/atmos_param/cosp/cosp_stats.F90 \
		 ${top_srcdir}/atmos_param/cosp/cosp_constants.lo \
		 ${top_srcdir}/atmos_param/cosp/cosp_types.lo \
		 ${top_srcdir}/atmos_param/cosp/llnl/llnl_stats.lo \
		 ${top_srcdir}/atmos_param/cosp/actsim/lmd_ipsl_stats.lo
${top_builddir}/atmos_param/cosp/cosp_utils.lo: ${top_srcdir}/atmos_param/cosp/cosp_utils.F90 \
		 ${top_srcdir}/atmos_param/cosp/cosp_constants.lo
${top_builddir}/atmos_param/cosp/cosp_types.lo: ${top_srcdir}/atmos_param/cosp/cosp_types.F90 \
		 ${top_srcdir}/atmos_param/cosp/cosp_constants.lo \
		 ${top_srcdir}/atmos_param/cosp/cosp_utils.lo \
		 ${top_srcdir}/atmos_param/cosp/quickbeam/radar_simulator_types.lo \
		 ${top_srcdir}/atmos_param/cosp/quickbeam/scale_LUTs_io.lo
${top_builddir}/atmos_param/cosp/cosp_rttov_simulator.lo: ${top_srcdir}/atmos_param/cosp/cosp_rttov_simulator.F90 \
		 ${top_srcdir}/atmos_param/cosp/cosp_constants.lo \
		 ${top_srcdir}/atmos_param/cosp/cosp_types.lo
${top_builddir}/atmos_param/cosp/cosp_isccp_simulator.lo: ${top_srcdir}/atmos_param/cosp/cosp_isccp_simulator.F90 \
		 ${top_srcdir}/atmos_param/cosp/cosp_constants.lo \
		 ${top_srcdir}/atmos_param/cosp/cosp_types.lo
${top_builddir}/atmos_param/cosp/cosp.lo: ${top_srcdir}/atmos_param/cosp/cosp.F90 \
		 ${top_srcdir}/atmos_param/cosp/cosp_types.lo \
		 ${top_srcdir}/atmos_param/cosp/cosp_simulator.lo \
		 ${top_srcdir}/atmos_param/cosp/cosp_modis_simulator.lo
${top_builddir}/atmos_param/cosp/cosp_constants.lo: ${top_srcdir}/atmos_param/cosp/cosp_constants.F90
${top_builddir}/atmos_param/cosp/cosp_diagnostics.lo: ${top_srcdir}/atmos_param/cosp/cosp_diagnostics.F90 \
		 ${top_srcdir}/atmos_param/cosp/cosp_types.lo \
		 ${top_srcdir}/atmos_param/cosp/cosp_io.lo \
		 ${top_srcdir}/atmos_param/cosp/cosp_constants.lo \
		 ${top_srcdir}/atmos_param/cosp/actsim/lmd_ipsl_stats.lo \
		 ${top_srcdir}/atmos_param/cosp/cosp_modis_simulator.lo \
		 ${top_srcdir}/atmos_param/cosp/MODIS_simulator/modis_simulator.lo \
		 ${top_srcdir}/atmos_param/cosp/cosp_utils.lo
${top_builddir}/atmos_param/cosp/cosp_simulator.lo: ${top_srcdir}/atmos_param/cosp/cosp_simulator.F90 \
		 ${top_srcdir}/atmos_param/cosp/cosp_constants.lo \
		 ${top_srcdir}/atmos_param/cosp/cosp_types.lo \
		 ${top_srcdir}/atmos_param/cosp/llnl/cosp_radar.lo \
		 ${top_srcdir}/atmos_param/cosp/cosp_lidar.lo \
		 ${top_srcdir}/atmos_param/cosp/cosp_isccp_simulator.lo \
		 ${top_srcdir}/atmos_param/cosp/cosp_modis_simulator.lo \
		 ${top_srcdir}/atmos_param/cosp/cosp_misr_simulator.lo \
		 ${top_srcdir}/atmos_param/cosp/cosp_rttov_simulator.lo \
		 ${top_srcdir}/atmos_param/cosp/cosp_stats.lo
${top_builddir}/atmos_param/cosp/cosp_lidar.lo: ${top_srcdir}/atmos_param/cosp/cosp_lidar.F90 \
		 ${top_srcdir}/atmos_param/cosp/cosp_constants.lo \
		 ${top_srcdir}/atmos_param/cosp/cosp_types.lo
${top_builddir}/atmos_param/cosp/cosp_io.lo: ${top_srcdir}/atmos_param/cosp/cosp_io.F90 \
		 ${top_srcdir}/atmos_param/cosp/cosp_constants.lo \
		 ${top_srcdir}/atmos_param/cosp/cosp_types.lo \
		 ${top_srcdir}/atmos_param/cosp/cosp_modis_simulator.lo
${top_builddir}/atmos_param/cosp/actsim/lidar_simulator.lo: ${top_srcdir}/atmos_param/cosp/actsim/lidar_simulator.F90
${top_builddir}/atmos_param/cosp/actsim/lmd_ipsl_stats.lo: ${top_srcdir}/atmos_param/cosp/actsim/lmd_ipsl_stats.F90 \
		 ${top_srcdir}/atmos_param/cosp/llnl/llnl_stats.lo \
		 ${top_srcdir}/atmos_param/cosp/cosp_constants.lo
${top_builddir}/atmos_param/cosp/MISR_simulator/MISR_simulator.lo: ${top_srcdir}/atmos_param/cosp/MISR_simulator/MISR_simulator.F
${top_builddir}/atmos_param/cosp/quickbeam/radar_simulator_types.lo: ${top_srcdir}/atmos_param/cosp/quickbeam/radar_simulator_types.F90
${top_builddir}/atmos_param/cosp/quickbeam/array_lib.lo: ${top_srcdir}/atmos_param/cosp/quickbeam/array_lib.F90 \
		 ${top_srcdir}/atmos_param/cosp/quickbeam/mrgrnk.lo
${top_builddir}/atmos_param/cosp/quickbeam/dsd.lo: ${top_srcdir}/atmos_param/cosp/quickbeam/dsd.F90 \
		 ${top_srcdir}/atmos_param/cosp/quickbeam/array_lib.lo \
		 ${top_srcdir}/atmos_param/cosp/quickbeam/math_lib.lo
${top_builddir}/atmos_param/cosp/quickbeam/optics_lib.lo: ${top_srcdir}/atmos_param/cosp/quickbeam/optics_lib.F90
${top_builddir}/atmos_param/cosp/quickbeam/gases.lo: ${top_srcdir}/atmos_param/cosp/quickbeam/gases.F90
${top_builddir}/atmos_param/cosp/quickbeam/format_input.lo: ${top_srcdir}/atmos_param/cosp/quickbeam/format_input.F90 \
		 ${top_srcdir}/atmos_param/cosp/quickbeam/array_lib.lo
${top_builddir}/atmos_param/cosp/quickbeam/mrgrnk.lo: ${top_srcdir}/atmos_param/cosp/quickbeam/mrgrnk.F90
${top_builddir}/atmos_param/cosp/quickbeam/radar_simulator.lo: ${top_srcdir}/atmos_param/cosp/quickbeam/radar_simulator.F90 \
		 ${top_srcdir}/atmos_param/cosp/quickbeam/mrgrnk.lo \
		 ${top_srcdir}/atmos_param/cosp/quickbeam/array_lib.lo \
		 ${top_srcdir}/atmos_param/cosp/quickbeam/math_lib.lo \
		 ${top_srcdir}/atmos_param/cosp/quickbeam/optics_lib.lo \
		 ${top_srcdir}/atmos_param/cosp/quickbeam/radar_simulator_types.lo \
		 ${top_srcdir}/atmos_param/cosp/quickbeam/scale_LUTs_io.lo
${top_builddir}/atmos_param/cosp/quickbeam/scale_LUTs_io.lo: ${top_srcdir}/atmos_param/cosp/quickbeam/scale_LUTs_io.F90 \
		 ${top_srcdir}/atmos_param/cosp/quickbeam/radar_simulator_types.lo
${top_builddir}/atmos_param/cosp/quickbeam/zeff.lo: ${top_srcdir}/atmos_param/cosp/quickbeam/zeff.F90 \
		 ${top_srcdir}/atmos_param/cosp/quickbeam/math_lib.lo \
		 ${top_srcdir}/atmos_param/cosp/quickbeam/optics_lib.lo
${top_builddir}/atmos_param/cosp/quickbeam/atmos_lib.lo: ${top_srcdir}/atmos_param/cosp/quickbeam/atmos_lib.F90
${top_builddir}/atmos_param/cosp/quickbeam/calc_Re.lo: ${top_srcdir}/atmos_param/cosp/quickbeam/calc_Re.F90 \
		 ${top_srcdir}/atmos_param/cosp/quickbeam/math_lib.lo
${top_builddir}/atmos_param/cosp/quickbeam/radar_simulator_init.lo: ${top_srcdir}/atmos_param/cosp/quickbeam/radar_simulator_init.F90 \
		 ${top_srcdir}/atmos_param/cosp/quickbeam/radar_simulator_types.lo
${top_builddir}/atmos_param/cosp/quickbeam/math_lib.lo: ${top_srcdir}/atmos_param/cosp/quickbeam/math_lib.F90 \
		 ${top_srcdir}/atmos_param/cosp/quickbeam/mrgrnk.lo \
		 ${top_srcdir}/atmos_param/cosp/quickbeam/array_lib.lo
${top_builddir}/atmos_param/cosp/llnl/prec_scops.lo: ${top_srcdir}/atmos_param/cosp/llnl/prec_scops.F
${top_builddir}/atmos_param/cosp/llnl/cosp_radar.lo: ${top_srcdir}/atmos_param/cosp/llnl/cosp_radar.F90 \
		 ${top_srcdir}/atmos_param/cosp/cosp_constants.lo \
		 ${top_srcdir}/atmos_param/cosp/cosp_types.lo \
		 ${top_srcdir}/atmos_param/cosp/cosp_utils.lo \
		 ${top_srcdir}/atmos_param/cosp/quickbeam/radar_simulator_types.lo \
		 ${top_srcdir}/atmos_param/cosp/quickbeam/array_lib.lo \
		 ${top_srcdir}/atmos_param/cosp/quickbeam/atmos_lib.lo \
		 ${top_srcdir}/atmos_param/cosp/quickbeam/format_input.lo \
		 ${top_srcdir}/atmos_param/cosp/quickbeam/mrgrnk.lo \
		 ${top_srcdir}/atmos_param/cosp/quickbeam/math_lib.lo \
		 ${top_srcdir}/atmos_param/cosp/quickbeam/optics_lib.lo
${top_builddir}/atmos_param/cosp/llnl/llnl_stats.lo: ${top_srcdir}/atmos_param/cosp/llnl/llnl_stats.F90 \
		 ${top_srcdir}/atmos_param/cosp/cosp_constants.lo
${top_builddir}/atmos_param/cosp/llnl/pf_to_mr.lo: ${top_srcdir}/atmos_param/cosp/llnl/pf_to_mr.F
${top_builddir}/atmos_param/cosp/icarus-scops-3.7/icarus.lo: ${top_srcdir}/atmos_param/cosp/icarus-scops-3.7/icarus.F90
${top_builddir}/atmos_param/cosp/icarus-scops-3.7/scops.lo: ${top_srcdir}/atmos_param/cosp/icarus-scops-3.7/scops.F
${top_builddir}/atmos_param/cosp/icarus-scops-3.7/isccp_cloud_types.lo: ${top_srcdir}/atmos_param/cosp/icarus-scops-3.7/isccp_cloud_types.F
${top_builddir}/atmos_param/cosp/MODIS_simulator/modis_simulator.lo: ${top_srcdir}/atmos_param/cosp/MODIS_simulator/modis_simulator.F90 \
		 ${top_srcdir}/atmos_param/cosp/cosp_types.lo
${top_builddir}/atmos_param/lscale_cond/lscale_cond.lo: ${top_srcdir}/atmos_param/lscale_cond/lscale_cond.F90
${top_builddir}/atmos_param/entrain/entrain.lo: ${top_srcdir}/atmos_param/entrain/entrain.F90
${top_builddir}/atmos_param/cg_drag/cg_drag.lo: ${top_srcdir}/atmos_param/cg_drag/cg_drag.F90
${top_builddir}/atmos_param/cu_mo_trans/cu_mo_trans.lo: ${top_srcdir}/atmos_param/cu_mo_trans/cu_mo_trans.F90 \
		 ${top_srcdir}/atmos_param/convection_driver/convection_utilities.lo \
		 ${top_srcdir}/atmos_param/moist_processes/moist_processes_utils.lo
${top_builddir}/atmos_param/cloud_rad/cloud_rad.lo: ${top_srcdir}/atmos_param/cloud_rad/cloud_rad.F90 \
		 ${top_srcdir}/atmos_param/microphysics/gamma_mg.lo \
		 ${top_srcdir}/atmos_param/lscloud_driver/lscloud_constants.lo \
		 ${top_srcdir}/atmos_shared/aerosol/aerosol_types.lo \
		 ${top_srcdir}/atmos_param/physics_radiation_exch/physics_radiation_exch.lo
${top_builddir}/atmos_param/NCARgwd/gw_utils.lo: ${top_srcdir}/atmos_param/NCARgwd/gw_utils.F90
${top_builddir}/atmos_param/NCARgwd/gw_common.lo: ${top_srcdir}/atmos_param/NCARgwd/gw_common.F90
${top_builddir}/atmos_param/NCARgwd/gw_drag.lo: ${top_srcdir}/atmos_param/NCARgwd/gw_drag.F90
${top_builddir}/atmos_param/NCARgwd/gw_convect_utils.lo: ${top_srcdir}/atmos_param/NCARgwd/gw_convect_utils.F90 \
		 ${top_builddir}/atmos_param/NCARgwd/gw_common.lo
${top_builddir}/atmos_param/NCARgwd/gw_convect.lo: ${top_srcdir}/atmos_param/NCARgwd/gw_convect.F90 \
		 ${top_builddir}/atmos_param/NCARgwd/gw_common.lo \
     ${top_builddir}/atmos_param/NCARgwd/gw_convect_utils.lo \
     ${top_builddir}/atmos_param/NCARgwd/gw_utils.lo
${top_builddir}/atmos_param/vert_diff/vert_diff.lo: ${top_srcdir}/atmos_param/vert_diff/vert_diff.F90
${top_builddir}/atmos_param/aerosol_cloud/aerosol_cloud.lo: ${top_srcdir}/atmos_param/aerosol_cloud/aerosol_cloud.F90 \
		 ${top_srcdir}/atmos_shared/aerosol/aerosol_types.lo \
		 ${top_srcdir}/atmos_param/lscloud_driver/polysvp.lo \
		 ${top_srcdir}/atmos_shared/tracer_driver/aer_ccn_act/aerosol_params.lo \
		 ${top_srcdir}/atmos_shared/tracer_driver/aer_ccn_act/aer_ccn_act.lo \
		 ${top_srcdir}/atmos_shared/tracer_driver/aer_ccn_act/ice_nucl.lo \
		 ${top_srcdir}/atmos_param/lscloud_driver/lscloud_types.lo \
		 ${top_srcdir}/atmos_param/lscloud_driver/lscloud_debug.lo \
		 ${top_srcdir}/atmos_param/moist_processes/moist_processes_utils.lo \
		 ${top_srcdir}/atmos_param/physics_radiation_exch/physics_radiation_exch.lo
${top_builddir}/atmos_param/mg_drag/mg_drag.lo: ${top_srcdir}/atmos_param/mg_drag/mg_drag.F90
${top_builddir}/atmos_param/vert_turb_driver/vert_turb_driver.lo: ${top_srcdir}/atmos_param/vert_turb_driver/vert_turb_driver.F90 \
		 ${top_srcdir}/atmos_param/physics_radiation_exch/physics_radiation_exch.lo \
		 ${top_srcdir}/atmos_param/physics_driver/physics_types.lo \
		 ${top_srcdir}/atmos_param/stable_bl_turb/stable_bl_turb.lo \
		 ${top_srcdir}/atmos_param/entrain/entrain.lo \
		 ${top_srcdir}/atmos_param/moist_processes/moist_processes_utils.lo \
		 ${top_srcdir}/atmos_shared/atmos_cmip_diag/atmos_cmip_diag.lo \
     ${top_srcdir}/atmos_param/han_edmf/hanedmf.lo \
     ${top_srcdir}/atmos_shared/aerosol/aerosol_types.lo
${top_builddir}/atmos_param/donner_deep/wet_deposition_0D.lo: ${top_srcdir}/atmos_param/donner_deep/wet_deposition_0D.F90
${top_builddir}/atmos_param/shallow_cu/conv_plumes_k.lo: ${top_srcdir}/atmos_param/shallow_cu/conv_plumes_k.F90 \
		 ${top_srcdir}/atmos_shared/tracer_driver/aer_ccn_act/aer_ccn_act_k.lo \
		 ${top_srcdir}/atmos_param/shallow_cu/conv_utilities_k.lo
${top_builddir}/atmos_param/shallow_cu/deep_conv.lo: ${top_srcdir}/atmos_param/shallow_cu/deep_conv.F90 \
		 ${top_srcdir}/atmos_param/shallow_cu/conv_utilities.lo \
		 ${top_srcdir}/atmos_param/shallow_cu/conv_utilities_k.lo \
		 ${top_srcdir}/atmos_param/shallow_cu/conv_plumes_k.lo \
		 ${top_srcdir}/atmos_param/shallow_cu/conv_closures.lo
${top_builddir}/atmos_param/shallow_cu/conv_closures.lo: ${top_srcdir}/atmos_param/shallow_cu/conv_closures.F90 \
		 ${top_srcdir}/atmos_param/shallow_cu/conv_utilities_k.lo \
		 ${top_srcdir}/atmos_param/shallow_cu/conv_plumes_k.lo
${top_builddir}/atmos_param/shallow_cu/conv_plumes.lo: ${top_srcdir}/atmos_param/shallow_cu/conv_plumes.F90
${top_builddir}/atmos_param/shallow_cu/conv_utilities_k.lo: ${top_srcdir}/atmos_param/shallow_cu/conv_utilities_k.F90
${top_builddir}/atmos_param/shallow_cu/conv_utilities.lo: ${top_srcdir}/atmos_param/shallow_cu/conv_utilities.F90 \
		 ${top_srcdir}/atmos_param/shallow_cu/conv_utilities_k.lo
${top_builddir}/atmos_param/shallow_cu/uw_conv.lo: ${top_srcdir}/atmos_param/shallow_cu/uw_conv.F90 \
		 ${top_srcdir}/atmos_shared/atmos_cmip_diag/atmos_cmip_diag.lo \
		 ${top_srcdir}/atmos_shared/tracer_driver/atmos_tracer_utilities.lo \
		 ${top_srcdir}/atmos_param/moist_processes/moist_processes_utils.lo \
		 ${top_srcdir}/atmos_shared/aerosol/aerosol_types.lo \
		 ${top_srcdir}/atmos_shared/tracer_driver/aer_ccn_act/aer_ccn_act.lo \
		 ${top_srcdir}/atmos_param/shallow_cu/conv_utilities.lo \
		 ${top_srcdir}/atmos_param/shallow_cu/conv_utilities_k.lo \
		 ${top_srcdir}/atmos_param/shallow_cu/conv_plumes_k.lo \
		 ${top_srcdir}/atmos_param/shallow_cu/conv_closures.lo \
		 ${top_srcdir}/atmos_param/shallow_cu/deep_conv.lo
${top_builddir}/atmos_param/han_edmf/hanedmf.lo: ${top_srcdir}/atmos_param/han_edmf/hanedmf.F90 \
		 ${top_srcdir}/atmos_shared/aerosol/aerosol_types.lo \
		 ${top_srcdir}/atmos_shared/tracer_driver/aer_ccn_act/aer_ccn_act.lo \
     ${top_srcdir}/atmos_param/han_edmf/satmedmfvdifq.lo \
     ${top_srcdir}/atmos_param/han_edmf/mfpbltq.lo \
     ${top_srcdir}/atmos_param/han_edmf/mfscuq.lo
${top_builddir}/atmos_param/han_edmf/tridi.lo: ${top_srcdir}/atmos_param/han_edmf/tridi.F90
${top_builddir}/atmos_param/han_edmf/mfpbltq.lo: ${top_srcdir}/atmos_param/han_edmf/mfpbltq.F90 \
     ${top_srcdir}/atmos_shared/tracer_driver/aer_ccn_act/aer_ccn_act_k.lo
${top_builddir}/atmos_param/han_edmf/mfscuq.lo: ${top_srcdir}/atmos_param/han_edmf/mfscuq.F90
${top_builddir}/atmos_param/han_edmf/satmedmfvdifq.lo: ${top_srcdir}/atmos_param/han_edmf/satmedmfvdifq.F90 \
     ${top_srcdir}/atmos_param/han_edmf/tridi.lo
atmos_nudge_mod.mod: ${top_srcdir}/atmos_shared/atmos_nudge/atmos_nudge.lo
atmos_global_diag_mod.mod: ${top_srcdir}/atmos_shared/atmos_cmip_diag/atmos_global_diag.lo
atmos_cmip_diag_mod.mod: ${top_srcdir}/atmos_shared/atmos_cmip_diag/atmos_cmip_diag.lo
vert_advection_mod.mod: ${top_srcdir}/atmos_shared/vert_advection/vert_advection.lo
aerosol_mod.mod: ${top_srcdir}/atmos_shared/aerosol/aerosol.lo
aerosol_types_mod.mod: ${top_srcdir}/atmos_shared/aerosol/aerosol_types.lo
atmos_convection_tracer_mod.mod: ${top_srcdir}/atmos_shared/tracer_driver/atmos_convection_tracer.lo
atmos_sulfur_hex_mod.mod: ${top_srcdir}/atmos_shared/tracer_driver/atmos_sulfur_hex.lo
atmos_tracer_driver_mod.mod: ${top_srcdir}/atmos_shared/tracer_driver/atmos_tracer_driver.lo
atmos_sea_salt_mod.mod: ${top_srcdir}/atmos_shared/tracer_driver/atmos_sea_salt.lo
atmos_tracer_utilities_mod.mod: ${top_srcdir}/atmos_shared/tracer_driver/atmos_tracer_utilities.lo
atmos_co2_mod.mod: ${top_srcdir}/atmos_shared/tracer_driver/atmos_co2.lo
tropchem_driver_mod.mod: ${top_srcdir}/atmos_shared/tracer_driver/tropchem_driver.lo
atmos_dust_mod.mod: ${top_srcdir}/atmos_shared/tracer_driver/atmos_dust.lo
atmos_tropopause_mod.mod: ${top_srcdir}/atmos_shared/tracer_driver/atmos_tropopause.lo
atmos_sulfate_mod.mod: ${top_srcdir}/atmos_shared/tracer_driver/atmos_sulfate.lo
atmos_radon_mod.mod: ${top_srcdir}/atmos_shared/tracer_driver/atmos_radon.lo
atmos_nh3_tag_mod.mod: ${top_srcdir}/atmos_shared/tracer_driver/atmos_nh3_tag.lo
atmos_soa_mod.mod: ${top_srcdir}/atmos_shared/tracer_driver/atmos_soa.lo
atmos_ch4_mod.mod: ${top_srcdir}/atmos_shared/tracer_driver/atmos_ch4.lo
atmos_regional_tracer_driver_mod.mod: ${top_srcdir}/atmos_shared/tracer_driver/atmos_regional_tracer_driver.lo
cloud_chem.mod: ${top_srcdir}/atmos_shared/tracer_driver/cloud_chem.lo
atmos_age_tracer_mod.mod: ${top_srcdir}/atmos_shared/tracer_driver/atmos_age_tracer.lo
mo_setsox_mod.mod: ${top_srcdir}/atmos_shared/tracer_driver/mo_setsox.lo
atmos_ch3i_mod.mod: ${top_srcdir}/atmos_shared/tracer_driver/atmos_ch3i.lo
xactive_bvoc_mod.mod: ${top_srcdir}/atmos_shared/tracer_driver/xactive_bvoc.lo
atmos_carbon_aerosol_mod.mod: ${top_srcdir}/atmos_shared/tracer_driver/atmos_carbon_aerosol.lo
strat_chem_mod.mod: ${top_srcdir}/atmos_shared/tracer_driver/stratchem/strat_chem_model.lo
strat_chem_driver_mod.mod: ${top_srcdir}/atmos_shared/tracer_driver/stratchem/strat_chem_driver.lo
mo_chemini_mod.mod: ${top_srcdir}/atmos_shared/tracer_driver/tropchem/mo_chemini.lo
mo_jpl_mod.mod: ${top_srcdir}/atmos_shared/tracer_driver/tropchem/mo_jpl.lo
mo_exp_sol_mod.mod: ${top_srcdir}/atmos_shared/tracer_driver/tropchem/mo_exp_slv.lo
mo_chem_utls_mod.mod: ${top_srcdir}/atmos_shared/tracer_driver/tropchem/mo_chem_utls.lo
tropchem_types_mod.mod: ${top_srcdir}/atmos_shared/tracer_driver/tropchem/tropchem_types_mod.lo
mo_usrrxt_mod.mod: ${top_srcdir}/atmos_shared/tracer_driver/tropchem/mo_usrrxt.lo
mo_photo_mod.mod: ${top_srcdir}/atmos_shared/tracer_driver/tropchem/mo_photo.lo
aerosol_thermodynamics.mod: ${top_srcdir}/atmos_shared/tracer_driver/tropchem/aerosol_thermodynamics.lo
mo_chemdr_mod.mod: ${top_srcdir}/atmos_shared/tracer_driver/tropchem/mo_chemdr.lo
strat_chem_utilities_mod.mod: ${top_srcdir}/atmos_shared/tracer_driver/tropchem/strat_chem_utilities.lo
mo_imp_sol_mod.mod: ${top_srcdir}/atmos_shared/tracer_driver/tropchem/mo_imp_slv.lo
moz_hook_mod.mod: ${top_srcdir}/atmos_shared/tracer_driver/tropchem/mo_hook.lo
mo_setinv_mod.mod: ${top_srcdir}/atmos_shared/tracer_driver/tropchem/mo_setinv.lo
mo_rodas_sol_mod.mod: ${top_srcdir}/atmos_shared/tracer_driver/tropchem/mo_rodas_slv.lo
mo_read_sim_chm_mod.mod: ${top_srcdir}/atmos_shared/tracer_driver/tropchem/mo_read_sim_chm.lo
m_tracname_mod.mod: ${top_srcdir}/atmos_shared/tracer_driver/tropchem/m_tracname.lo
mo_grid_mod.mod: ${top_srcdir}/atmos_shared/tracer_driver/tropchem/AM4_chem/moz.mods.lo
chem_mods_mod.mod: ${top_srcdir}/atmos_shared/tracer_driver/tropchem/AM4_chem/moz.mods.lo
m_spc_id_mod.mod: ${top_srcdir}/atmos_shared/tracer_driver/tropchem/AM4_chem/moz.mods.lo
m_rxt_id_mod.mod: ${top_srcdir}/atmos_shared/tracer_driver/tropchem/AM4_chem/moz.mods.lo
m_het_id_mod.mod: ${top_srcdir}/atmos_shared/tracer_driver/tropchem/AM4_chem/moz.mods.lo
mo_fphoto_mod.mod: ${top_srcdir}/atmos_shared/tracer_driver/tropchem/AM4_chem/mo_fphoto.lo
mo_exp_prod_loss_mod.mod: ${top_srcdir}/atmos_shared/tracer_driver/tropchem/AM4_chem/moz.mat.lo
mo_imp_prod_loss_mod.mod: ${top_srcdir}/atmos_shared/tracer_driver/tropchem/AM4_chem/moz.mat.lo
mo_rodas_prod_loss_mod.mod: ${top_srcdir}/atmos_shared/tracer_driver/tropchem/AM4_chem/moz.mat.lo
mo_indprd_mod.mod: ${top_srcdir}/atmos_shared/tracer_driver/tropchem/AM4_chem/moz.mat.lo
mo_imp_lin_matrix_mod.mod: ${top_srcdir}/atmos_shared/tracer_driver/tropchem/AM4_chem/moz.mat.lo
mo_rod_lin_matrix_mod.mod: ${top_srcdir}/atmos_shared/tracer_driver/tropchem/AM4_chem/moz.mat.lo
mo_imp_nln_matrix_mod.mod: ${top_srcdir}/atmos_shared/tracer_driver/tropchem/AM4_chem/moz.mat.lo
mo_rod_nln_matrix_mod.mod: ${top_srcdir}/atmos_shared/tracer_driver/tropchem/AM4_chem/moz.mat.lo
mo_imp_factor_mod.mod: ${top_srcdir}/atmos_shared/tracer_driver/tropchem/AM4_chem/moz.mat.lo
mo_rod_factor_mod.mod: ${top_srcdir}/atmos_shared/tracer_driver/tropchem/AM4_chem/moz.mat.lo
mo_imp_solve_mod.mod: ${top_srcdir}/atmos_shared/tracer_driver/tropchem/AM4_chem/moz.mat.lo
mo_rod_solve_mod.mod: ${top_srcdir}/atmos_shared/tracer_driver/tropchem/AM4_chem/moz.mat.lo
mo_setrxt_mod.mod: ${top_srcdir}/atmos_shared/tracer_driver/tropchem/AM4_chem/moz.subs.lo
mo_adjrxt_mod.mod: ${top_srcdir}/atmos_shared/tracer_driver/tropchem/AM4_chem/moz.subs.lo
mo_phtadj_mod.mod: ${top_srcdir}/atmos_shared/tracer_driver/tropchem/AM4_chem/moz.subs.lo
mo_rxt_mod.mod: ${top_srcdir}/atmos_shared/tracer_driver/tropchem/AM4_chem/moz.subs.lo
mo_make_grp_vmr_mod.mod: ${top_srcdir}/atmos_shared/tracer_driver/tropchem/AM4_chem/moz.subs.lo
mo_fastjx_mod.mod: ${top_srcdir}/atmos_shared/tracer_driver/tropchem/AM4_chem/mo_fastjx.lo
aer_ccn_act_mod.mod: ${top_srcdir}/atmos_shared/tracer_driver/aer_ccn_act/aer_ccn_act.lo
aer_in_act_mod.mod: ${top_srcdir}/atmos_shared/tracer_driver/aer_ccn_act/aer_in_act.lo
aer_ccn_act_k_mod.mod: ${top_srcdir}/atmos_shared/tracer_driver/aer_ccn_act/aer_ccn_act_k.lo
ice_nucl_mod.mod: ${top_srcdir}/atmos_shared/tracer_driver/aer_ccn_act/ice_nucl.lo
aerosol_params_mod.mod: ${top_srcdir}/atmos_shared/tracer_driver/aer_ccn_act/aerosol_params.lo
atmos_ozone_tracer_mod.mod: ${top_srcdir}/atmos_shared/tracer_driver/atmos_ozone_tracer.lo
physics_types_mod.mod: ${top_srcdir}/atmos_param/physics_driver/physics_types.lo
physics_driver_mod.mod: ${top_srcdir}/atmos_param/physics_driver/physics_driver.lo
topo_drag_mod.mod: ${top_srcdir}/atmos_param/topo_drag/topo_drag.lo
ls_cloud_macrophysics_mod.mod: ${top_srcdir}/atmos_param/macrophysics/ls_cloud_macrophysics.lo
tiedtke_macro_mod.mod: ${top_srcdir}/atmos_param/macrophysics/tiedtke_macro.lo
cloudrad_types_mod.mod: ${top_srcdir}/atmos_param/radiation/clouds/cloudrad_types.lo
microphys_rad_mod.mod: ${top_srcdir}/atmos_param/radiation/clouds/microphys_rad.lo
strat_clouds_w_mod.mod: ${top_srcdir}/atmos_param/radiation/clouds/strat_clouds_W.lo
isccp_clouds_mod.mod: ${top_srcdir}/atmos_param/radiation/clouds/isccp_clouds.lo
random_number_streams_mod.mod: ${top_srcdir}/atmos_param/radiation/clouds/get_random_number_stream.lo
uw_clouds_w_mod.mod: ${top_srcdir}/atmos_param/radiation/clouds/uw_clouds_W.lo
cloudrad_driver_mod.mod: ${top_srcdir}/atmos_param/radiation/clouds/cloudrad_driver.lo
cloud_spec_mod.mod: ${top_srcdir}/atmos_param/radiation/clouds/cloud_spec.lo
cloudrad_diagnostics_mod.mod: ${top_srcdir}/atmos_param/radiation/clouds/cloudrad_diagnostics.lo
cloudrad_package_mod.mod: ${top_srcdir}/atmos_param/radiation/clouds/cloudrad_package.lo
longwave_tables_mod.mod: ${top_srcdir}/atmos_param/radiation/longwave/longwave_tables.lo
longwave_clouds_mod.mod: ${top_srcdir}/atmos_param/radiation/longwave/longwave_clouds.lo
gas_tf_mod.mod: ${top_srcdir}/atmos_param/radiation/longwave/gas_tf.lo
longwave_params_mod.mod: ${top_srcdir}/atmos_param/radiation/longwave/longwave_params.lo
sealw99_mod.mod: ${top_srcdir}/atmos_param/radiation/longwave/sealw99.lo
optical_path_mod.mod: ${top_srcdir}/atmos_param/radiation/longwave/optical_path.lo
longwave_types_mod.mod: ${top_srcdir}/atmos_param/radiation/longwave/longwave_types.lo
longwave_utilities_mod.mod: ${top_srcdir}/atmos_param/radiation/longwave/longwave_utilities.lo
longwave_fluxes_mod.mod: ${top_srcdir}/atmos_param/radiation/longwave/longwave_fluxes.lo
lw_gases_stdtf_mod.mod: ${top_srcdir}/atmos_param/radiation/longwave/lw_gases_stdtf.lo
aerosolrad_types_mod.mod: ${top_srcdir}/atmos_param/radiation/aerosols/aerosolrad_types.lo
aerosolrad_package_mod.mod: ${top_srcdir}/atmos_param/radiation/aerosols/aerosolrad_package.lo
aerosolrad_driver_mod.mod: ${top_srcdir}/atmos_param/radiation/aerosols/aerosolrad_driver.lo
esfsw_driver_mod.mod: ${top_srcdir}/atmos_param/radiation/shortwave/esfsw_driver.lo
esfsw_utilities_mod.mod: ${top_srcdir}/atmos_param/radiation/shortwave/esfsw_utilities.lo
shortwave_types_mod.mod: ${top_srcdir}/atmos_param/radiation/shortwave/shortwave_types.lo
esfsw_bands_mod.mod: ${top_srcdir}/atmos_param/radiation/shortwave/esfsw_bands.lo
esfsw_parameters_mod.mod: ${top_srcdir}/atmos_param/radiation/shortwave/esfsw_parameters.lo
radiative_gases_mod.mod: ${top_srcdir}/atmos_param/radiation/radiative_gases/radiative_gases.lo
radiative_gases_types_mod.mod: ${top_srcdir}/atmos_param/radiation/radiative_gases/radiative_gases_types.lo
ozone_mod.mod: ${top_srcdir}/atmos_param/radiation/radiative_gases/ozone.lo
radiation_diag_mod.mod: ${top_srcdir}/atmos_param/radiation/util/radiation_diag.lo
longwave_driver_mod.mod: ${top_srcdir}/atmos_param/radiation/driver/longwave_driver.lo
rad_output_file_mod.mod: ${top_srcdir}/atmos_param/radiation/driver/rad_output_file.lo
radiation_driver_types_mod.mod: ${top_srcdir}/atmos_param/radiation/driver/radiation_driver_types.lo
radiation_types_mod.mod: ${top_srcdir}/atmos_param/radiation/driver/radiation_types.lo
solar_data_driver_mod.mod: ${top_srcdir}/atmos_param/radiation/driver/solar_data_driver.lo
shortwave_driver_mod.mod: ${top_srcdir}/atmos_param/radiation/driver/shortwave_driver.lo
radiation_driver_diag_mod.mod: ${top_srcdir}/atmos_param/radiation/driver/radiation_driver_diag.lo
radiation_driver_mod.mod: ${top_srcdir}/atmos_param/radiation/driver/radiation_driver.lo
convection_utilities_mod.mod: ${top_srcdir}/atmos_param/convection_driver/convection_utilities.lo
detr_ice_num_mod.mod: ${top_srcdir}/atmos_param/convection_driver/detr_ice_num.lo
convection_driver_mod.mod: ${top_srcdir}/atmos_param/convection_driver/convection_driver.lo
physics_radiation_exch_mod.mod: ${top_srcdir}/atmos_param/physics_radiation_exch/physics_radiation_exch.lo
beta_dist_mod.mod: ${top_srcdir}/atmos_param/cloud_generator/betaDistribution.lo
cloud_generator_mod.mod: ${top_srcdir}/atmos_param/cloud_generator/cloud_generator.lo
lscloud_debug_mod.mod: ${top_srcdir}/atmos_param/lscloud_driver/lscloud_debug.lo
lscloud_constants_mod.mod: ${top_srcdir}/atmos_param/lscloud_driver/lscloud_constants.lo
lscloud_netcdf_mod.mod: ${top_srcdir}/atmos_param/lscloud_driver/lscloud_netcdf.lo
lscloud_types_mod.mod: ${top_srcdir}/atmos_param/lscloud_driver/lscloud_types.lo
lscloud_driver_mod.mod: ${top_srcdir}/atmos_param/lscloud_driver/lscloud_driver.lo
polysvp_mod.mod: ${top_srcdir}/atmos_param/lscloud_driver/polysvp.lo
check_nan_mod.mod: ${top_srcdir}/atmos_param/lscloud_driver/check_nan.lo
stable_bl_turb_mod.mod: ${top_srcdir}/atmos_param/stable_bl_turb/stable_bl_turb.lo
gamma_mg_mod.mod: ${top_srcdir}/atmos_param/microphysics/gamma_mg.lo
micro_mg2_mod.mod: ${top_srcdir}/atmos_param/microphysics/micro_mg2.lo
micro_mg2_utils.mod: ${top_srcdir}/atmos_param/microphysics/micro_mg2_utils.lo
rotstayn_klein_mp_mod.mod: ${top_srcdir}/atmos_param/microphysics/rotstayn_klein_mp.lo
ls_cloud_microphysics_mod.mod: ${top_srcdir}/atmos_param/microphysics/ls_cloud_microphysics.lo
damping_driver_mod.mod: ${top_srcdir}/atmos_param/damping_driver/damping_driver.lo
vert_diff_driver_mod.mod: ${top_srcdir}/atmos_param/vert_diff_driver/vert_diff_driver.lo
moist_processes_mod.mod: ${top_srcdir}/atmos_param/moist_processes/moist_processes.lo
moist_proc_utils_mod.mod: ${top_srcdir}/atmos_param/moist_processes/moist_processes_utils.lo
mod_cosp_modis_simulator.mod: ${top_srcdir}/atmos_param/cosp/cosp_modis_simulator.lo
cosp_driver_mod.mod: ${top_srcdir}/atmos_param/cosp/cosp_driver.lo
mod_cosp_misr_simulator.mod: ${top_srcdir}/atmos_param/cosp/cosp_misr_simulator.lo
mod_cosp_stats.mod: ${top_srcdir}/atmos_param/cosp/cosp_stats.lo
mod_cosp_utils.mod: ${top_srcdir}/atmos_param/cosp/cosp_utils.lo
mod_cosp_types.mod: ${top_srcdir}/atmos_param/cosp/cosp_types.lo
mod_cosp_rttov_simulator.mod: ${top_srcdir}/atmos_param/cosp/cosp_rttov_simulator.lo
mod_cosp_isccp_simulator.mod: ${top_srcdir}/atmos_param/cosp/cosp_isccp_simulator.lo
mod_cosp.mod: ${top_srcdir}/atmos_param/cosp/cosp.lo
mod_cosp_constants.mod: ${top_srcdir}/atmos_param/cosp/cosp_constants.lo
cosp_diagnostics_mod.mod: ${top_srcdir}/atmos_param/cosp/cosp_diagnostics.lo
mod_cosp_simulator.mod: ${top_srcdir}/atmos_param/cosp/cosp_simulator.lo
mod_cosp_lidar.mod: ${top_srcdir}/atmos_param/cosp/cosp_lidar.lo
mod_cosp_io.mod: ${top_srcdir}/atmos_param/cosp/cosp_io.lo
mod_lmd_ipsl_stats.mod: ${top_srcdir}/atmos_param/cosp/actsim/lmd_ipsl_stats.lo
radar_simulator_types.mod: ${top_srcdir}/atmos_param/cosp/quickbeam/radar_simulator_types.lo
array_lib.mod: ${top_srcdir}/atmos_param/cosp/quickbeam/array_lib.lo
optics_lib.mod: ${top_srcdir}/atmos_param/cosp/quickbeam/optics_lib.lo
format_input.mod: ${top_srcdir}/atmos_param/cosp/quickbeam/format_input.lo
m_mrgrnk.mod: ${top_srcdir}/atmos_param/cosp/quickbeam/mrgrnk.lo
scale_luts_io.mod: ${top_srcdir}/atmos_param/cosp/quickbeam/scale_LUTs_io.lo
atmos_lib.mod: ${top_srcdir}/atmos_param/cosp/quickbeam/atmos_lib.lo
math_lib.mod: ${top_srcdir}/atmos_param/cosp/quickbeam/math_lib.lo
mod_cosp_radar.mod: ${top_srcdir}/atmos_param/cosp/llnl/cosp_radar.lo
mod_llnl_stats.mod: ${top_srcdir}/atmos_param/cosp/llnl/llnl_stats.lo
mod_modis_sim.mod: ${top_srcdir}/atmos_param/cosp/MODIS_simulator/modis_simulator.lo
lscale_cond_mod.mod: ${top_srcdir}/atmos_param/lscale_cond/lscale_cond.lo
entrain_mod.mod: ${top_srcdir}/atmos_param/entrain/entrain.lo
cg_drag_mod.mod: ${top_srcdir}/atmos_param/cg_drag/cg_drag.lo
cu_mo_trans_mod.mod: ${top_srcdir}/atmos_param/cu_mo_trans/cu_mo_trans.lo
cloud_rad_mod.mod: ${top_srcdir}/atmos_param/cloud_rad/cloud_rad.lo
gw_common.mod : ${top_srcdir}/atmos_param/NCARgwd/gw_common.lo
gw_drag.mod : ${top_srcdir}/atmos_param/NCARgwd/gw_drag.lo
gw_convect.mod : ${top_srcdir}/atmos_param/NCARgwd/gw_convect.lo
gw_convect_utils.mod : ${top_srcdir}/atmos_param/NCARgwd/gw_convect_utils.lo
gw_utils.mod : ${top_srcdir}/atmos_param/NCARgwd/gw_utils.lo
vert_diff_mod.mod: ${top_srcdir}/atmos_param/vert_diff/vert_diff.lo
aerosol_cloud_mod.mod: ${top_srcdir}/atmos_param/aerosol_cloud/aerosol_cloud.lo
mg_drag_mod.mod: ${top_srcdir}/atmos_param/mg_drag/mg_drag.lo
vert_turb_driver_mod.mod: ${top_srcdir}/atmos_param/vert_turb_driver/vert_turb_driver.lo
conv_plumes_k_mod.mod: ${top_srcdir}/atmos_param/shallow_cu/conv_plumes_k.lo
deep_conv_mod.mod: ${top_srcdir}/atmos_param/shallow_cu/deep_conv.lo
conv_closures_mod.mod: ${top_srcdir}/atmos_param/shallow_cu/conv_closures.lo
conv_plumes_mod.mod: ${top_srcdir}/atmos_param/shallow_cu/conv_plumes.lo
conv_utilities_k_mod.mod: ${top_srcdir}/atmos_param/shallow_cu/conv_utilities_k.lo
conv_utilities_mod.mod: ${top_srcdir}/atmos_param/shallow_cu/conv_utilities.lo
uw_conv_mod.mod: ${top_srcdir}/atmos_param/shallow_cu/uw_conv.lo
hanedmf_mod.mod: ${top_srcdir}/atmos_param/han_edmf/hanedmf.lo
satmedmfvdifq_mod.mod: ${top_srcdir}/atmos_param/han_edmf/satmedmfvdifq.lo
tridi_mod.mod: ${top_srcdir}/atmos_param/han_edmf/tridi.lo
include_HEADERS = atmos_nudge_mod.mod \
		atmos_global_diag_mod.mod \
		atmos_cmip_diag_mod.mod \
		vert_advection_mod.mod \
		aerosol_mod.mod \
		aerosol_types_mod.mod \
		atmos_convection_tracer_mod.mod \
		atmos_sulfur_hex_mod.mod \
		atmos_tracer_driver_mod.mod \
		atmos_sea_salt_mod.mod \
		atmos_tracer_utilities_mod.mod \
		atmos_co2_mod.mod \
		tropchem_driver_mod.mod \
		atmos_dust_mod.mod \
		atmos_tropopause_mod.mod \
		atmos_sulfate_mod.mod \
		atmos_radon_mod.mod \
		atmos_nh3_tag_mod.mod \
		atmos_soa_mod.mod \
		atmos_ch4_mod.mod \
		atmos_regional_tracer_driver_mod.mod \
		cloud_chem.mod \
		atmos_age_tracer_mod.mod \
		mo_setsox_mod.mod \
		atmos_ch3i_mod.mod \
		xactive_bvoc_mod.mod \
		atmos_carbon_aerosol_mod.mod \
		strat_chem_mod.mod \
		strat_chem_driver_mod.mod \
		mo_chemini_mod.mod \
		mo_jpl_mod.mod \
		mo_exp_sol_mod.mod \
		mo_chem_utls_mod.mod \
		tropchem_types_mod.mod \
		mo_usrrxt_mod.mod \
		mo_photo_mod.mod \
		aerosol_thermodynamics.mod \
		mo_chemdr_mod.mod \
		strat_chem_utilities_mod.mod \
		mo_imp_sol_mod.mod \
		moz_hook_mod.mod \
		mo_setinv_mod.mod \
		mo_rodas_sol_mod.mod \
		mo_read_sim_chm_mod.mod \
		m_tracname_mod.mod \
		mo_grid_mod.mod \
		chem_mods_mod.mod \
		m_spc_id_mod.mod \
		m_rxt_id_mod.mod \
		m_het_id_mod.mod \
		mo_fphoto_mod.mod \
		mo_exp_prod_loss_mod.mod \
		mo_imp_prod_loss_mod.mod \
		mo_rodas_prod_loss_mod.mod \
		mo_indprd_mod.mod \
		mo_imp_lin_matrix_mod.mod \
		mo_rod_lin_matrix_mod.mod \
		mo_imp_nln_matrix_mod.mod \
		mo_rod_nln_matrix_mod.mod \
		mo_imp_factor_mod.mod \
		mo_rod_factor_mod.mod \
		mo_imp_solve_mod.mod \
		mo_rod_solve_mod.mod \
		mo_setrxt_mod.mod \
		mo_adjrxt_mod.mod \
		mo_phtadj_mod.mod \
		mo_rxt_mod.mod \
		mo_make_grp_vmr_mod.mod \
		mo_fastjx_mod.mod \
		aer_ccn_act_mod.mod \
		aer_in_act_mod.mod \
		aer_ccn_act_k_mod.mod \
		ice_nucl_mod.mod \
		aerosol_params_mod.mod \
		physics_types_mod.mod \
		physics_driver_mod.mod \
		topo_drag_mod.mod \
		ls_cloud_macrophysics_mod.mod \
		tiedtke_macro_mod.mod \
		cloudrad_types_mod.mod \
		microphys_rad_mod.mod \
		strat_clouds_w_mod.mod \
		isccp_clouds_mod.mod \
		random_number_streams_mod.mod \
		uw_clouds_w_mod.mod \
		cloudrad_driver_mod.mod \
		cloud_spec_mod.mod \
		cloudrad_diagnostics_mod.mod \
		cloudrad_package_mod.mod \
		longwave_tables_mod.mod \
		longwave_clouds_mod.mod \
		gas_tf_mod.mod \
		longwave_params_mod.mod \
		sealw99_mod.mod \
		optical_path_mod.mod \
		longwave_types_mod.mod \
		longwave_utilities_mod.mod \
		longwave_fluxes_mod.mod \
		lw_gases_stdtf_mod.mod \
		aerosolrad_types_mod.mod \
		aerosolrad_package_mod.mod \
		aerosolrad_driver_mod.mod \
		esfsw_driver_mod.mod \
		esfsw_utilities_mod.mod \
		shortwave_types_mod.mod \
		esfsw_bands_mod.mod \
		esfsw_parameters_mod.mod \
		radiative_gases_mod.mod \
		radiative_gases_types_mod.mod \
		ozone_mod.mod \
		radiation_diag_mod.mod \
		longwave_driver_mod.mod \
		rad_output_file_mod.mod \
		radiation_driver_types_mod.mod \
		radiation_types_mod.mod \
		solar_data_driver_mod.mod \
		shortwave_driver_mod.mod \
		radiation_driver_diag_mod.mod \
		radiation_driver_mod.mod \
		convection_utilities_mod.mod \
		detr_ice_num_mod.mod \
		convection_driver_mod.mod \
		physics_radiation_exch_mod.mod \
		beta_dist_mod.mod \
		cloud_generator_mod.mod \
		lscloud_debug_mod.mod \
		lscloud_constants_mod.mod \
		lscloud_netcdf_mod.mod \
		lscloud_types_mod.mod \
		lscloud_driver_mod.mod \
		polysvp_mod.mod \
		check_nan_mod.mod \
		stable_bl_turb_mod.mod \
		gamma_mg_mod.mod \
		micro_mg2_mod.mod \
		micro_mg2_utils.mod \
		rotstayn_klein_mp_mod.mod \
		ls_cloud_microphysics_mod.mod \
		damping_driver_mod.mod \
		vert_diff_driver_mod.mod \
		moist_processes_mod.mod \
		moist_proc_utils_mod.mod \
		mod_cosp_modis_simulator.mod \
		cosp_driver_mod.mod \
		mod_cosp_misr_simulator.mod \
		mod_cosp_stats.mod \
		mod_cosp_utils.mod \
		mod_cosp_types.mod \
		mod_cosp_rttov_simulator.mod \
		mod_cosp_isccp_simulator.mod \
		mod_cosp.mod \
		mod_cosp_constants.mod \
		cosp_diagnostics_mod.mod \
		mod_cosp_simulator.mod \
		mod_cosp_lidar.mod \
		mod_cosp_io.mod \
		mod_lmd_ipsl_stats.mod \
		radar_simulator_types.mod \
		array_lib.mod \
		optics_lib.mod \
		format_input.mod \
		m_mrgrnk.mod \
		scale_luts_io.mod \
		atmos_lib.mod \
		math_lib.mod \
		mod_cosp_radar.mod \
		mod_llnl_stats.mod \
		mod_modis_sim.mod \
		lscale_cond_mod.mod \
		entrain_mod.mod \
		cg_drag_mod.mod \
		cu_mo_trans_mod.mod \
		cloud_rad_mod.mod \
		vert_diff_mod.mod \
		aerosol_cloud_mod.mod \
		mg_drag_mod.mod \
		vert_turb_driver_mod.mod \
		conv_plumes_k_mod.mod \
		deep_conv_mod.mod \
		conv_closures_mod.mod \
		conv_plumes_mod.mod \
		conv_utilities_k_mod.mod \
		conv_utilities_mod.mod \
		uw_conv_mod.mod \
		atmos_ozone_tracer_mod.mod \
    gw_beres_mod.mod \
    gw_common_utils.mod \
    hanedmf_mod.mod \
    tridi_mod.mod \
    satmedmfvdifq_mod.mod
CLEANFILES = *.mod
