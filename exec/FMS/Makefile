# Makefile created by mkmf 2022.01.00

CPPDEFS =  -DHAVE_GETTID -DINTERNAL_FILE_NML -Duse_libMPI -Duse_netCDF -DMAXFIELDMETHODS_=500
include $(MK_TEMPLATE)


.DEFAULT:
	-echo $@ does not exist.
all: libFMS.a
affinity.o: $(SRCROOT)/FMS/affinity/affinity.c
	$(CC) $(CPPDEFS) $(CPPFLAGS) $(CFLAGS) $(OTHERFLAGS) $(OTHER_CFLAGS) -c	$(SRCROOT)/FMS/affinity/affinity.c
amip_interp.o: $(SRCROOT)/FMS/amip_interp/amip_interp.F90 $(SRCROOT)/FMS/include/file_version.h $(SRCROOT)/FMS/amip_interp/include/amip_interp_r4.fh $(SRCROOT)/FMS/amip_interp/include/amip_interp.inc $(SRCROOT)/FMS/amip_interp/include/amip_interp_r8.fh time_interp.o time_manager.o get_cal_time.o horiz_interp.o fms.o constants.o platform.o mpp.o fms2_io.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)/FMS/include -I$(SRCROOT)/FMS/amip_interp/include	$(SRCROOT)/FMS/amip_interp/amip_interp.F90
astronomy.o: $(SRCROOT)/FMS/astronomy/astronomy.F90 $(SRCROOT)/FMS/include/file_version.h $(SRCROOT)/FMS/astronomy/include/astronomy_r4.fh $(SRCROOT)/FMS/astronomy/include/astronomy.inc $(SRCROOT)/FMS/astronomy/include/astronomy_r8.fh fms.o time_manager.o fms_io_utils.o constants.o mpp.o platform.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)/FMS/include -I$(SRCROOT)/FMS/astronomy/include	$(SRCROOT)/FMS/astronomy/astronomy.F90
atmos_ocean_fluxes.o: $(SRCROOT)/FMS/coupler/atmos_ocean_fluxes.F90 $(SRCROOT)/FMS/include/file_version.h mpp.o fms.o coupler_types.o field_manager.o fm_util.o fms_io_utils.o platform.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)/FMS/include	$(SRCROOT)/FMS/coupler/atmos_ocean_fluxes.F90
axis_utils.o: $(SRCROOT)/FMS/axis_utils/axis_utils.F90 $(SRCROOT)/FMS/include/file_version.h mpp_io.o mpp.o fms.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)/FMS/include	$(SRCROOT)/FMS/axis_utils/axis_utils.F90
axis_utils2.o: $(SRCROOT)/FMS/axis_utils/axis_utils2.F90 $(SRCROOT)/FMS/include/file_version.h $(SRCROOT)/FMS/axis_utils/include/axis_utils2_r4.fh $(SRCROOT)/FMS/axis_utils/include/axis_utils2.inc $(SRCROOT)/FMS/axis_utils/include/axis_utils2_r8.fh mpp.o fms.o fms2_io.o platform.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)/FMS/include -I$(SRCROOT)/FMS/axis_utils/include	$(SRCROOT)/FMS/axis_utils/axis_utils2.F90
blackboxio.o: $(SRCROOT)/FMS/fms2_io/blackboxio.F90 mpp_domains.o fms_io_utils.o netcdf_io.o fms_netcdf_domain_io.o fms_netcdf_unstructured_domain_io.o mpp.o platform.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)/FMS/fms2_io/blackboxio.F90
block_control.o: $(SRCROOT)/FMS/block_control/block_control.F90 mpp.o mpp_domains.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)/FMS/block_control/block_control.F90
cloud_interpolator.o: $(SRCROOT)/FMS/drifters/cloud_interpolator.F90 $(SRCROOT)/FMS/include/file_version.h
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)/FMS/include	$(SRCROOT)/FMS/drifters/cloud_interpolator.F90
column_diagnostics.o: $(SRCROOT)/FMS/column_diagnostics/column_diagnostics.F90 $(SRCROOT)/FMS/include/file_version.h $(SRCROOT)/FMS/column_diagnostics/include/column_diagnostics_r4.fh $(SRCROOT)/FMS/column_diagnostics/include/column_diagnostics.inc $(SRCROOT)/FMS/column_diagnostics/include/column_diagnostics_r8.fh fms.o time_manager.o constants.o mpp.o platform.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)/FMS/include -I$(SRCROOT)/FMS/column_diagnostics/include	$(SRCROOT)/FMS/column_diagnostics/column_diagnostics.F90
constants.o: $(SRCROOT)/FMS/constants/constants.F90 fmsconstants.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)/FMS/constants/constants.F90
constantsr4.o: $(SRCROOT)/FMS/constants4/constantsr4.F90 fmsconstantsr4.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)/FMS/constants4/constantsr4.F90
coupler_types.o: $(SRCROOT)/FMS/coupler/coupler_types.F90 $(SRCROOT)/FMS/include/file_version.h $(SRCROOT)/FMS/coupler/include/coupler_types_r4.fh $(SRCROOT)/FMS/coupler/include/coupler_types.inc $(SRCROOT)/FMS/coupler/include/coupler_types_r8.fh fms.o fms2_io.o fms_io.o time_manager.o diag_manager.o data_override.o mpp_domains.o mpp.o fms_string_utils.o platform.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)/FMS/include -I$(SRCROOT)/FMS/coupler/include	$(SRCROOT)/FMS/coupler/coupler_types.F90
create_xgrid.o: $(SRCROOT)/FMS/mosaic/create_xgrid.c $(SRCROOT)/FMS/mosaic/mosaic_util.h $(SRCROOT)/FMS/mosaic/create_xgrid.h $(SRCROOT)/FMS/mosaic/constant.h
	$(CC) $(CPPDEFS) $(CPPFLAGS) $(CFLAGS) $(OTHERFLAGS) $(OTHER_CFLAGS) -c -I$(SRCROOT)/FMS/mosaic	$(SRCROOT)/FMS/mosaic/create_xgrid.c
data_override.o: $(SRCROOT)/FMS/data_override/data_override.F90 $(SRCROOT)/FMS/data_override/include/data_override_r4.fh $(SRCROOT)/FMS/data_override/include/data_override.inc $(SRCROOT)/FMS/include/file_version.h $(SRCROOT)/FMS/data_override/include/data_override_r8.fh platform.o mpp.o mpp_domains.o time_manager.o yaml_parser.o constants.o horiz_interp.o time_interp_external2.o fms.o axis_utils2.o fms2_io.o get_grid_version.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)/FMS/data_override/include -I$(SRCROOT)/FMS/include	$(SRCROOT)/FMS/data_override/data_override.F90
diag_axis.o: $(SRCROOT)/FMS/diag_manager/diag_axis.F90 $(SRCROOT)/FMS/include/file_version.h platform.o mpp_domains.o fms.o diag_data.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)/FMS/include	$(SRCROOT)/FMS/diag_manager/diag_axis.F90
diag_data.o: $(SRCROOT)/FMS/diag_manager/diag_data.F90 $(SRCROOT)/FMS/include/file_version.h platform.o time_manager.o mpp_domains.o fms.o fms_diag_bbox.o fms2_io.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)/FMS/include	$(SRCROOT)/FMS/diag_manager/diag_data.F90
diag_grid.o: $(SRCROOT)/FMS/diag_manager/diag_grid.F90 $(SRCROOT)/FMS/include/file_version.h platform.o constants.o fms.o mpp.o mpp_domains.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)/FMS/include	$(SRCROOT)/FMS/diag_manager/diag_grid.F90
diag_integral.o: $(SRCROOT)/FMS/diag_integral/diag_integral.F90 $(SRCROOT)/FMS/include/file_version.h $(SRCROOT)/FMS/diag_integral/include/diag_integral_r4.fh $(SRCROOT)/FMS/diag_integral/include/diag_integral.inc $(SRCROOT)/FMS/diag_integral/include/diag_integral_r8.fh platform.o time_manager.o mpp.o fms.o fms2_io.o constants.o ensemble_manager.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)/FMS/include -I$(SRCROOT)/FMS/diag_integral/include	$(SRCROOT)/FMS/diag_integral/diag_integral.F90
diag_manager.o: $(SRCROOT)/FMS/diag_manager/diag_manager.F90 $(SRCROOT)/FMS/include/file_version.h platform.o time_manager.o mpp.o fms.o diag_axis.o diag_util.o diag_data.o diag_table.o diag_output.o diag_grid.o constants.o fms_diag_outfield.o fms_diag_fieldbuff_update.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)/FMS/include	$(SRCROOT)/FMS/diag_manager/diag_manager.F90
diag_output.o: $(SRCROOT)/FMS/diag_manager/diag_output.F90 platform.o mpp_domains.o mpp.o diag_axis.o diag_data.o time_manager.o fms.o fms2_io.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)/FMS/diag_manager/diag_output.F90
diag_table.o: $(SRCROOT)/FMS/diag_manager/diag_table.F90 fms2_io.o fms.o time_manager.o constants.o diag_data.o diag_util.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)/FMS/diag_manager/diag_table.F90
diag_util.o: $(SRCROOT)/FMS/diag_manager/diag_util.F90 $(SRCROOT)/FMS/include/file_version.h platform.o diag_data.o diag_axis.o diag_output.o diag_grid.o fms.o mpp_domains.o time_manager.o mpp.o constants.o fms2_io.o fms_diag_bbox.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)/FMS/include	$(SRCROOT)/FMS/diag_manager/diag_util.F90
drifters.o: $(SRCROOT)/FMS/drifters/drifters.F90 $(SRCROOT)/FMS/drifters/fms_switches.h $(SRCROOT)/FMS/include/file_version.h $(SRCROOT)/FMS/drifters/drifters_push.fh $(SRCROOT)/FMS/drifters/drifters_set_field.fh $(SRCROOT)/FMS/drifters/drifters_compute_k.fh mpp.o mpp_domains.o drifters_core.o drifters_input.o drifters_io.o drifters_comm.o cloud_interpolator.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)/FMS/drifters -I$(SRCROOT)/FMS/include	$(SRCROOT)/FMS/drifters/drifters.F90
drifters_comm.o: $(SRCROOT)/FMS/drifters/drifters_comm.F90 $(SRCROOT)/FMS/drifters/fms_switches.h mpp.o mpp_domains.o drifters_core.o drifters_input.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)/FMS/drifters	$(SRCROOT)/FMS/drifters/drifters_comm.F90
drifters_core.o: $(SRCROOT)/FMS/drifters/drifters_core.F90 $(SRCROOT)/FMS/include/file_version.h platform.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)/FMS/include	$(SRCROOT)/FMS/drifters/drifters_core.F90
drifters_input.o: $(SRCROOT)/FMS/drifters/drifters_input.F90 $(SRCROOT)/FMS/include/file_version.h
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)/FMS/include	$(SRCROOT)/FMS/drifters/drifters_input.F90
drifters_io.o: $(SRCROOT)/FMS/drifters/drifters_io.F90 $(SRCROOT)/FMS/include/file_version.h
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)/FMS/include	$(SRCROOT)/FMS/drifters/drifters_io.F90
ensemble_manager.o: $(SRCROOT)/FMS/coupler/ensemble_manager.F90 fms.o mpp.o fms2_io.o fms_io.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)/FMS/coupler/ensemble_manager.F90
field_manager.o: $(SRCROOT)/FMS/field_manager/field_manager.F90 $(SRCROOT)/FMS/include/file_version.h $(SRCROOT)/FMS/field_manager/parse.inc $(SRCROOT)/FMS/field_manager/include/field_manager_r4.fh $(SRCROOT)/FMS/field_manager/include/field_manager.inc $(SRCROOT)/FMS/field_manager/include/field_manager_r8.fh mpp.o fms.o fms2_io.o platform.o fm_yaml.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)/FMS/include -I$(SRCROOT)/FMS/field_manager -I$(SRCROOT)/FMS/field_manager/include	$(SRCROOT)/FMS/field_manager/field_manager.F90
fm_util.o: $(SRCROOT)/FMS/field_manager/fm_util.F90 $(SRCROOT)/FMS/include/file_version.h $(SRCROOT)/FMS/field_manager/include/fm_util_r4.fh $(SRCROOT)/FMS/field_manager/include/fm_util.inc $(SRCROOT)/FMS/field_manager/include/fm_util_r8.fh field_manager.o fms.o mpp.o platform.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)/FMS/include -I$(SRCROOT)/FMS/field_manager/include	$(SRCROOT)/FMS/field_manager/fm_util.F90
fm_yaml.o: $(SRCROOT)/FMS/field_manager/fm_yaml.F90 yaml_parser.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)/FMS/field_manager/fm_yaml.F90
fms.o: $(SRCROOT)/FMS/fms/fms.F90 $(SRCROOT)/FMS/include/file_version.h $(SRCROOT)/FMS/fms/include/fms_r4.fh $(SRCROOT)/FMS/fms/include/fms.inc $(SRCROOT)/FMS/fms/include/fms_r8.fh mpp.o mpp_domains.o mpp_io.o fms_io.o fms2_io.o memutils.o grid2.o fms_string_utils.o platform.o constants.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)/FMS/include -I$(SRCROOT)/FMS/fms/include	$(SRCROOT)/FMS/fms/fms.F90
fms2_io.o: $(SRCROOT)/FMS/fms2_io/fms2_io.F90 fms_io_utils.o netcdf_io.o fms_netcdf_domain_io.o fms_netcdf_unstructured_domain_io.o blackboxio.o mpp.o mpp_domains.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)/FMS/fms2_io/fms2_io.F90
fms_affinity.o: $(SRCROOT)/FMS/affinity/fms_affinity.F90 $(SRCROOT)/FMS/include/file_version.h mpp.o fms.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)/FMS/include	$(SRCROOT)/FMS/affinity/fms_affinity.F90
fms_diag_bbox.o: $(SRCROOT)/FMS/diag_manager/fms_diag_bbox.F90 fms.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)/FMS/diag_manager/fms_diag_bbox.F90
fms_diag_elem_weight_procs.o: $(SRCROOT)/FMS/diag_manager/fms_diag_elem_weight_procs.F90 platform.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)/FMS/diag_manager/fms_diag_elem_weight_procs.F90
fms_diag_fieldbuff_update.o: $(SRCROOT)/FMS/diag_manager/fms_diag_fieldbuff_update.F90 $(SRCROOT)/FMS/diag_manager/include/fms_diag_fieldbuff_update.inc $(SRCROOT)/FMS/diag_manager/include/fms_diag_fieldbuff_update.fh platform.o mpp.o time_manager.o fms.o diag_data.o fms_diag_outfield.o diag_util.o fms_diag_time_reduction.o fms_diag_elem_weight_procs.o fms_diag_bbox.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)/FMS/diag_manager/include	$(SRCROOT)/FMS/diag_manager/fms_diag_fieldbuff_update.F90
fms_diag_outfield.o: $(SRCROOT)/FMS/diag_manager/fms_diag_outfield.F90 platform.o mpp.o fms.o diag_data.o fms_diag_time_reduction.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)/FMS/diag_manager/fms_diag_outfield.F90
fms_diag_time_reduction.o: $(SRCROOT)/FMS/diag_manager/fms_diag_time_reduction.F90 diag_data.o fms.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)/FMS/diag_manager/fms_diag_time_reduction.F90
fms_io.o: $(SRCROOT)/FMS/fms/fms_io.F90 $(SRCROOT)/FMS/include/fms_platform.h $(SRCROOT)/FMS/include/file_version.h $(SRCROOT)/FMS/fms/read_data_2d.inc $(SRCROOT)/FMS/fms/read_data_3d.inc $(SRCROOT)/FMS/fms/read_data_4d.inc $(SRCROOT)/FMS/fms/write_data.inc $(SRCROOT)/FMS/fms/fms_io_unstructured_register_restart_axis.inc $(SRCROOT)/FMS/fms/fms_io_unstructured_setup_one_field.inc $(SRCROOT)/FMS/fms/fms_io_unstructured_register_restart_field.inc $(SRCROOT)/FMS/fms/fms_io_unstructured_save_restart.inc $(SRCROOT)/FMS/fms/fms_io_unstructured_read.inc $(SRCROOT)/FMS/fms/fms_io_unstructured_get_file_name.inc $(SRCROOT)/FMS/fms/fms_io_unstructured_get_file_unit.inc $(SRCROOT)/FMS/fms/fms_io_unstructured_file_unit.inc $(SRCROOT)/FMS/fms/fms_io_unstructured_get_field_size.inc $(SRCROOT)/FMS/fms/fms_io_unstructured_field_exist.inc mpp_io.o mpp_domains.o mpp.o platform.o mpp_parameter.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)/FMS/include -I$(SRCROOT)/FMS/fms	$(SRCROOT)/FMS/fms/fms_io.F90
fms_io_utils.o: $(SRCROOT)/FMS/fms2_io/fms_io_utils.F90 $(SRCROOT)/FMS/fms2_io/include/array_utils.inc $(SRCROOT)/FMS/fms2_io/include/array_utils_char.inc $(SRCROOT)/FMS/fms2_io/include/get_data_type_string.inc mpp.o mpp_domains.o platform.o fms_string_utils.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)/FMS/fms2_io/include	$(SRCROOT)/FMS/fms2_io/fms_io_utils.F90
fms_netcdf_domain_io.o: $(SRCROOT)/FMS/fms2_io/fms_netcdf_domain_io.F90 $(SRCROOT)/FMS/fms2_io/include/register_domain_restart_variable.inc $(SRCROOT)/FMS/fms2_io/include/domain_read.inc $(SRCROOT)/FMS/fms2_io/include/domain_write.inc $(SRCROOT)/FMS/fms2_io/include/compute_global_checksum.inc mpp.o mpp_domains.o fms_io_utils.o netcdf_io.o platform.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)/FMS/fms2_io/include	$(SRCROOT)/FMS/fms2_io/fms_netcdf_domain_io.F90
fms_netcdf_unstructured_domain_io.o: $(SRCROOT)/FMS/fms2_io/fms_netcdf_unstructured_domain_io.F90 $(SRCROOT)/FMS/fms2_io/include/register_unstructured_domain_restart_variable.inc $(SRCROOT)/FMS/fms2_io/include/unstructured_domain_read.inc $(SRCROOT)/FMS/fms2_io/include/unstructured_domain_write.inc mpp_domains.o fms_io_utils.o netcdf_io.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)/FMS/fms2_io/include	$(SRCROOT)/FMS/fms2_io/fms_netcdf_unstructured_domain_io.F90
fms_stacksize.o: $(SRCROOT)/FMS/fms/fms_stacksize.c
	$(CC) $(CPPDEFS) $(CPPFLAGS) $(CFLAGS) $(OTHERFLAGS) $(OTHER_CFLAGS) -c	$(SRCROOT)/FMS/fms/fms_stacksize.c
fms_string_utils.o: $(SRCROOT)/FMS/string_utils/fms_string_utils.F90 $(SRCROOT)/FMS/string_utils/include/fms_string_utils_r4.fh $(SRCROOT)/FMS/string_utils/include/fms_string_utils.inc $(SRCROOT)/FMS/string_utils/include/fms_string_utils_r8.fh platform.o mpp.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)/FMS/string_utils/include	$(SRCROOT)/FMS/string_utils/fms_string_utils.F90
fms_string_utils_binding.o: $(SRCROOT)/FMS/string_utils/fms_string_utils_binding.c
	$(CC) $(CPPDEFS) $(CPPFLAGS) $(CFLAGS) $(OTHERFLAGS) $(OTHER_CFLAGS) -c	$(SRCROOT)/FMS/string_utils/fms_string_utils_binding.c
fms_yaml_output.o: $(SRCROOT)/FMS/parser/fms_yaml_output.F90 fms_string_utils.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)/FMS/parser/fms_yaml_output.F90
fmsconstants.o: $(SRCROOT)/FMS/constants/fmsconstants.F90 $(SRCROOT)/FMS/constants/gfdl_constants.fh $(SRCROOT)/FMS/constants/gfs_constants.fh $(SRCROOT)/FMS/constants/geos_constants.fh platform.o mpp.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)/FMS/constants	$(SRCROOT)/FMS/constants/fmsconstants.F90
fmsconstantsr4.o: $(SRCROOT)/FMS/constants4/fmsconstantsr4.F90 $(SRCROOT)/FMS/constants4/gfdl_constantsR4.fh $(SRCROOT)/FMS/constants4/gfs_constantsR4.fh $(SRCROOT)/FMS/constants4/geos_constantsR4.fh platform.o mpp.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)/FMS/constants4	$(SRCROOT)/FMS/constants4/fmsconstantsr4.F90
gaussian_topog.o: $(SRCROOT)/FMS/topography/gaussian_topog.F90 $(SRCROOT)/FMS/include/file_version.h $(SRCROOT)/FMS/topography/include/gaussian_topog_r4.fh $(SRCROOT)/FMS/topography/include/gaussian_topog.inc $(SRCROOT)/FMS/topography/include/gaussian_topog_r8.fh fms.o constants.o mpp.o platform.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)/FMS/include -I$(SRCROOT)/FMS/topography/include	$(SRCROOT)/FMS/topography/gaussian_topog.F90
get_cal_time.o: $(SRCROOT)/FMS/time_manager/get_cal_time.F90 $(SRCROOT)/FMS/include/file_version.h fms.o time_manager.o mpp.o platform.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)/FMS/include	$(SRCROOT)/FMS/time_manager/get_cal_time.F90
get_grid_version.o: $(SRCROOT)/FMS/data_override/get_grid_version.F90 $(SRCROOT)/FMS/data_override/include/get_grid_version_r4.fh $(SRCROOT)/FMS/data_override/include/get_grid_version.inc $(SRCROOT)/FMS/data_override/include/get_grid_version_r8.fh constants.o platform.o mpp.o mpp_domains.o fms2_io.o mosaic2.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)/FMS/data_override/include	$(SRCROOT)/FMS/data_override/get_grid_version.F90
gradient.o: $(SRCROOT)/FMS/mosaic/gradient.F90 $(SRCROOT)/FMS/include/file_version.h mpp.o constants.o platform.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)/FMS/include	$(SRCROOT)/FMS/mosaic/gradient.F90
gradient_c2l.o: $(SRCROOT)/FMS/mosaic/gradient_c2l.c $(SRCROOT)/FMS/mosaic/constant.h $(SRCROOT)/FMS/mosaic/mosaic_util.h $(SRCROOT)/FMS/mosaic/gradient_c2l.h
	$(CC) $(CPPDEFS) $(CPPFLAGS) $(CFLAGS) $(OTHERFLAGS) $(OTHER_CFLAGS) -c -I$(SRCROOT)/FMS/mosaic	$(SRCROOT)/FMS/mosaic/gradient_c2l.c
grid.o: $(SRCROOT)/FMS/mosaic/grid.F90 $(SRCROOT)/FMS/include/file_version.h mpp.o constants.o fms_io.o mosaic.o mpp_domains.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)/FMS/include	$(SRCROOT)/FMS/mosaic/grid.F90
grid2.o: $(SRCROOT)/FMS/mosaic2/grid2.F90 $(SRCROOT)/FMS/include/file_version.h $(SRCROOT)/FMS/mosaic2/include/grid2_r4.fh $(SRCROOT)/FMS/mosaic2/include/grid2.inc $(SRCROOT)/FMS/mosaic2/include/grid2_r8.fh mpp.o constants.o fms2_io.o fms_string_utils.o mosaic2.o mpp_domains.o platform.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)/FMS/include -I$(SRCROOT)/FMS/mosaic2/include	$(SRCROOT)/FMS/mosaic2/grid2.F90
horiz_interp.o: $(SRCROOT)/FMS/horiz_interp/horiz_interp.F90 $(SRCROOT)/FMS/include/file_version.h $(SRCROOT)/FMS/horiz_interp/include/horiz_interp_r4.fh $(SRCROOT)/FMS/horiz_interp/include/horiz_interp.inc $(SRCROOT)/FMS/horiz_interp/include/horiz_interp_r8.fh fms.o mpp.o constants.o horiz_interp_type.o horiz_interp_conserve.o horiz_interp_bilinear.o horiz_interp_bicubic.o horiz_interp_spherical.o platform.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)/FMS/include -I$(SRCROOT)/FMS/horiz_interp/include	$(SRCROOT)/FMS/horiz_interp/horiz_interp.F90
horiz_interp_bicubic.o: $(SRCROOT)/FMS/horiz_interp/horiz_interp_bicubic.F90 $(SRCROOT)/FMS/include/file_version.h $(SRCROOT)/FMS/horiz_interp/include/horiz_interp_bicubic_r4.fh $(SRCROOT)/FMS/horiz_interp/include/horiz_interp_bicubic.inc $(SRCROOT)/FMS/horiz_interp/include/horiz_interp_bicubic_r8.fh mpp.o fms.o horiz_interp_type.o constants.o platform.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)/FMS/include -I$(SRCROOT)/FMS/horiz_interp/include	$(SRCROOT)/FMS/horiz_interp/horiz_interp_bicubic.F90
horiz_interp_bilinear.o: $(SRCROOT)/FMS/horiz_interp/horiz_interp_bilinear.F90 $(SRCROOT)/FMS/include/file_version.h $(SRCROOT)/FMS/horiz_interp/include/horiz_interp_bilinear_r4.fh $(SRCROOT)/FMS/horiz_interp/include/horiz_interp_bilinear.inc $(SRCROOT)/FMS/horiz_interp/include/horiz_interp_bilinear_r8.fh mpp.o fms.o constants.o horiz_interp_type.o platform.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)/FMS/include -I$(SRCROOT)/FMS/horiz_interp/include	$(SRCROOT)/FMS/horiz_interp/horiz_interp_bilinear.F90
horiz_interp_conserve.o: $(SRCROOT)/FMS/horiz_interp/horiz_interp_conserve.F90 $(SRCROOT)/FMS/include/file_version.h $(SRCROOT)/FMS/horiz_interp/include/horiz_interp_conserve_r4.fh $(SRCROOT)/FMS/horiz_interp/include/horiz_interp_conserve.inc $(SRCROOT)/FMS/horiz_interp/include/horiz_interp_conserve_r8.fh platform.o mpp.o fms.o grid2.o constants.o horiz_interp_type.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)/FMS/include -I$(SRCROOT)/FMS/horiz_interp/include	$(SRCROOT)/FMS/horiz_interp/horiz_interp_conserve.F90
horiz_interp_spherical.o: $(SRCROOT)/FMS/horiz_interp/horiz_interp_spherical.F90 $(SRCROOT)/FMS/include/file_version.h $(SRCROOT)/FMS/horiz_interp/include/horiz_interp_spherical_r4.fh $(SRCROOT)/FMS/horiz_interp/include/horiz_interp_spherical.inc $(SRCROOT)/FMS/horiz_interp/include/horiz_interp_spherical_r8.fh platform.o mpp.o fms.o constants.o horiz_interp_type.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)/FMS/include -I$(SRCROOT)/FMS/horiz_interp/include	$(SRCROOT)/FMS/horiz_interp/horiz_interp_spherical.F90
horiz_interp_type.o: $(SRCROOT)/FMS/horiz_interp/horiz_interp_type.F90 $(SRCROOT)/FMS/horiz_interp/include/horiz_interp_type_r4.fh $(SRCROOT)/FMS/horiz_interp/include/horiz_interp_type.inc $(SRCROOT)/FMS/horiz_interp/include/horiz_interp_type_r8.fh mpp.o platform.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)/FMS/horiz_interp/include	$(SRCROOT)/FMS/horiz_interp/horiz_interp_type.F90
interp.o: $(SRCROOT)/FMS/mosaic/interp.c $(SRCROOT)/FMS/mosaic/mosaic_util.h $(SRCROOT)/FMS/mosaic/interp.h $(SRCROOT)/FMS/mosaic/create_xgrid.h
	$(CC) $(CPPDEFS) $(CPPFLAGS) $(CFLAGS) $(OTHERFLAGS) $(OTHER_CFLAGS) -c -I$(SRCROOT)/FMS/mosaic	$(SRCROOT)/FMS/mosaic/interp.c
interpolator.o: $(SRCROOT)/FMS/interpolator/interpolator.F90 $(SRCROOT)/FMS/include/file_version.h $(SRCROOT)/FMS/interpolator/include/interpolator_r4.fh $(SRCROOT)/FMS/interpolator/include/interpolator.inc $(SRCROOT)/FMS/interpolator/include/interpolator_r8.fh mpp.o mpp_domains.o diag_manager.o fms.o fms2_io.o horiz_interp.o time_manager.o time_interp.o constants.o platform.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)/FMS/include -I$(SRCROOT)/FMS/interpolator/include	$(SRCROOT)/FMS/interpolator/interpolator.F90
libFMS.o: $(SRCROOT)/FMS/libFMS.F90 $(SRCROOT)/FMS/include/file_version.h fms_affinity.o amip_interp.o astronomy.o axis_utils2.o block_control.o column_diagnostics.o coupler_types.o ensemble_manager.o atmos_ocean_fluxes.o data_override.o diag_integral.o diag_manager.o xgrid.o stock_constants.o field_manager.o fm_util.o fms2_io.o fms.o horiz_interp.o horiz_interp_type.o interpolator.o memutils.o monin_obukhov.o monin_obukhov_inter.o mosaic2.o grid2.o gradient.o mpp.o mpp_parameter.o mpp_utilities.o mpp_memutils.o mpp_efp.o mpp_domains.o yaml_parser.o platform.o random_numbers.o sat_vapor_pres.o fms_string_utils.o time_interp.o time_interp_external2.o time_manager.o get_cal_time.o gaussian_topog.o topography.o tracer_manager.o tridiagonal.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)/FMS/include	$(SRCROOT)/FMS/libFMS.F90
memutils.o: $(SRCROOT)/FMS/memutils/memutils.F90 mpp.o mpp_memutils.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)/FMS/memutils/memutils.F90
mersennetwister.o: $(SRCROOT)/FMS/random_numbers/mersennetwister.F90 platform.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)/FMS/random_numbers/mersennetwister.F90
monin_obukhov.o: $(SRCROOT)/FMS/monin_obukhov/monin_obukhov.F90 $(SRCROOT)/FMS/include/file_version.h $(SRCROOT)/FMS/monin_obukhov/include/monin_obukhov_r4.fh $(SRCROOT)/FMS/monin_obukhov/include/monin_obukhov.inc $(SRCROOT)/FMS/monin_obukhov/include/monin_obukhov_r8.fh constants.o mpp.o fms.o monin_obukhov_inter.o platform.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)/FMS/include -I$(SRCROOT)/FMS/monin_obukhov/include	$(SRCROOT)/FMS/monin_obukhov/monin_obukhov.F90
monin_obukhov_inter.o: $(SRCROOT)/FMS/monin_obukhov/monin_obukhov_inter.F90 $(SRCROOT)/FMS/monin_obukhov/include/monin_obukhov_inter_r4.fh $(SRCROOT)/FMS/monin_obukhov/include/monin_obukhov_inter.inc $(SRCROOT)/FMS/monin_obukhov/include/monin_obukhov_inter_r8.fh platform.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)/FMS/monin_obukhov/include	$(SRCROOT)/FMS/monin_obukhov/monin_obukhov_inter.F90
mosaic.o: $(SRCROOT)/FMS/mosaic/mosaic.F90 $(SRCROOT)/FMS/include/file_version.h mpp.o mpp_io.o fms_io.o constants.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)/FMS/include	$(SRCROOT)/FMS/mosaic/mosaic.F90
mosaic2.o: $(SRCROOT)/FMS/mosaic2/mosaic2.F90 $(SRCROOT)/FMS/include/file_version.h $(SRCROOT)/FMS/mosaic2/include/mosaic2_r4.fh $(SRCROOT)/FMS/mosaic2/include/mosaic2.inc $(SRCROOT)/FMS/mosaic2/include/mosaic2_r8.fh mpp.o mpp_domains.o constants.o fms2_io.o platform.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)/FMS/include -I$(SRCROOT)/FMS/mosaic2/include	$(SRCROOT)/FMS/mosaic2/mosaic2.F90
mosaic_util.o: $(SRCROOT)/FMS/mosaic/mosaic_util.c $(SRCROOT)/FMS/mosaic/mosaic_util.h $(SRCROOT)/FMS/mosaic/constant.h
	$(CC) $(CPPDEFS) $(CPPFLAGS) $(CFLAGS) $(OTHERFLAGS) $(OTHER_CFLAGS) -c -I$(SRCROOT)/FMS/mosaic	$(SRCROOT)/FMS/mosaic/mosaic_util.c
mpp.o: $(SRCROOT)/FMS/mpp/mpp.F90 $(SRCROOT)/FMS/include/file_version.h $(SRCROOT)/FMS/mpp/include/system_clock.fh $(SRCROOT)/FMS/mpp/include/mpp_util.inc $(SRCROOT)/FMS/mpp/include/mpp_util_mpi.inc $(SRCROOT)/FMS/mpp/include/mpp_util_nocomm.inc $(SRCROOT)/FMS/mpp/include/mpp_error_a_a.fh $(SRCROOT)FMS/mpp/include/mpp_error_a_a.fh $(SRCROOT)/FMS/mpp/include/mpp_error_a_s.fh $(SRCROOT)FMS/mpp/include/mpp_error_a_s.fh $(SRCROOT)/FMS/mpp/include/mpp_error_s_a.fh $(SRCROOT)FMS/mpp/include/mpp_error_s_a.fh $(SRCROOT)/FMS/mpp/include/mpp_error_s_s.fh $(SRCROOT)FMS/mpp/include/mpp_error_s_s.fh $(SRCROOT)FMS/include/file_version.h $(SRCROOT)/FMS/mpp/include/mpp_comm.inc $(SRCROOT)/FMS/mpp/include/mpp_comm_mpi.inc $(SRCROOT)/FMS/mpp/include/mpp_transmit_mpi.fh $(SRCROOT)/FMS/mpp/include/mpp_transmit.inc $(SRCROOT)FMS/mpp/include/mpp_transmit_mpi.fh $(SRCROOT)FMS/mpp/include/mpp_transmit.inc $(SRCROOT)/FMS/mpp/include/mpp_reduce_mpi.fh $(SRCROOT)FMS/mpp/include/mpp_reduce_mpi.fh $(SRCROOT)/FMS/mpp/include/mpp_sum_mpi.fh $(SRCROOT)/FMS/mpp/include/mpp_sum.inc $(SRCROOT)FMS/mpp/include/mpp_sum_mpi.fh $(SRCROOT)FMS/mpp/include/mpp_sum.inc $(SRCROOT)/FMS/mpp/include/mpp_sum_mpi_ad.fh $(SRCROOT)/FMS/mpp/include/mpp_sum_ad.inc $(SRCROOT)FMS/mpp/include/mpp_sum_mpi_ad.fh $(SRCROOT)FMS/mpp/include/mpp_sum_ad.inc $(SRCROOT)/FMS/mpp/include/mpp_alltoall_mpi.fh $(SRCROOT)FMS/mpp/include/mpp_alltoall_mpi.fh $(SRCROOT)/FMS/mpp/include/mpp_type_mpi.fh $(SRCROOT)FMS/mpp/include/mpp_type_mpi.fh $(SRCROOT)/FMS/mpp/include/mpp_comm_nocomm.inc $(SRCROOT)/FMS/mpp/include/mpp_transmit_nocomm.fh $(SRCROOT)FMS/mpp/include/mpp_transmit_nocomm.fh $(SRCROOT)/FMS/mpp/include/mpp_reduce_nocomm.fh $(SRCROOT)FMS/mpp/include/mpp_reduce_nocomm.fh $(SRCROOT)/FMS/mpp/include/mpp_sum_nocomm.fh $(SRCROOT)FMS/mpp/include/mpp_sum_nocomm.fh $(SRCROOT)/FMS/mpp/include/mpp_sum_nocomm_ad.fh $(SRCROOT)FMS/mpp/include/mpp_sum_nocomm_ad.fh $(SRCROOT)/FMS/mpp/include/mpp_alltoall_nocomm.fh $(SRCROOT)FMS/mpp/include/mpp_alltoall_nocomm.fh $(SRCROOT)/FMS/mpp/include/mpp_type_nocomm.fh $(SRCROOT)FMS/mpp/include/mpp_type_nocomm.fh $(SRCROOT)/FMS/mpp/include/mpp_chksum_int.fh $(SRCROOT)FMS/mpp/include/mpp_chksum_int.fh $(SRCROOT)/FMS/mpp/include/mpp_chksum_scalar.fh $(SRCROOT)/FMS/mpp/include/mpp_chksum.fh $(SRCROOT)FMS/mpp/include/mpp_chksum.fh $(SRCROOT)FMS/mpp/include/mpp_chksum_scalar.fh $(SRCROOT)/FMS/mpp/include/mpp_gather.fh $(SRCROOT)FMS/mpp/include/mpp_gather.fh $(SRCROOT)/FMS/mpp/include/mpp_scatter.fh $(SRCROOT)FMS/mpp/include/mpp_scatter.fh mpp_parameter.o mpp_data.o platform.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)/FMS/include -I$(SRCROOT)/FMS/mpp/include -I$(SRCROOT)FMS/mpp/include -I$(SRCROOT)FMS/include	$(SRCROOT)/FMS/mpp/mpp.F90
mpp_data.o: $(SRCROOT)/FMS/mpp/mpp_data.F90 $(SRCROOT)/FMS/include/file_version.h $(SRCROOT)/FMS/mpp/include/mpp_data_mpi.inc $(SRCROOT)/FMS/mpp/include/mpp_data_nocomm.inc mpp_parameter.o platform.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)/FMS/include -I$(SRCROOT)/FMS/mpp/include	$(SRCROOT)/FMS/mpp/mpp_data.F90
mpp_domains.o: $(SRCROOT)/FMS/mpp/mpp_domains.F90 $(SRCROOT)/FMS/include/file_version.h $(SRCROOT)/FMS/mpp/include/mpp_define_nest_domains.inc $(SRCROOT)/FMS/mpp/include/mpp_domains_util.inc $(SRCROOT)/FMS/mpp/include/mpp_domains_comm.inc $(SRCROOT)/FMS/mpp/include/mpp_domains_define.inc $(SRCROOT)/FMS/mpp/include/mpp_domains_misc.inc $(SRCROOT)/FMS/mpp/include/mpp_update_domains2D.fh $(SRCROOT)FMS/mpp/include/mpp_update_domains2D.fh $(SRCROOT)/FMS/mpp/include/mpp_update_domains2D_nonblock.fh $(SRCROOT)FMS/mpp/include/mpp_update_domains2D_nonblock.fh $(SRCROOT)/FMS/mpp/include/mpp_do_update_nonblock.fh $(SRCROOT)/FMS/mpp/include/mpp_do_updateV_nonblock.fh $(SRCROOT)FMS/mpp/include/mpp_do_update_nonblock.fh $(SRCROOT)FMS/mpp/include/mpp_do_updateV_nonblock.fh $(SRCROOT)/FMS/mpp/include/mpp_do_update.fh $(SRCROOT)/FMS/mpp/include/mpp_do_updateV.fh $(SRCROOT)FMS/mpp/include/mpp_do_update.fh $(SRCROOT)FMS/mpp/include/mpp_do_updateV.fh $(SRCROOT)/FMS/mpp/include/mpp_do_check.fh $(SRCROOT)/FMS/mpp/include/mpp_do_checkV.fh $(SRCROOT)FMS/mpp/include/mpp_do_check.fh $(SRCROOT)FMS/mpp/include/mpp_do_checkV.fh $(SRCROOT)/FMS/mpp/include/mpp_update_nest_domains.fh $(SRCROOT)FMS/mpp/include/mpp_update_nest_domains.fh $(SRCROOT)/FMS/mpp/include/mpp_do_update_nest.fh $(SRCROOT)FMS/mpp/include/mpp_do_update_nest.fh $(SRCROOT)/FMS/mpp/include/mpp_update_domains2D_ad.fh $(SRCROOT)FMS/mpp/include/mpp_update_domains2D_ad.fh $(SRCROOT)/FMS/mpp/include/mpp_do_update_ad.fh $(SRCROOT)/FMS/mpp/include/mpp_do_updateV_ad.fh $(SRCROOT)FMS/mpp/include/mpp_do_update_ad.fh $(SRCROOT)FMS/mpp/include/mpp_do_updateV_ad.fh $(SRCROOT)/FMS/mpp/include/mpp_do_redistribute.fh $(SRCROOT)FMS/mpp/include/mpp_do_redistribute.fh $(SRCROOT)/FMS/mpp/include/mpp_get_boundary.fh $(SRCROOT)/FMS/mpp/include/mpp_get_boundary_ad.fh $(SRCROOT)FMS/mpp/include/mpp_get_boundary.fh $(SRCROOT)FMS/mpp/include/mpp_get_boundary_ad.fh $(SRCROOT)/FMS/mpp/include/mpp_do_get_boundary.fh $(SRCROOT)/FMS/mpp/include/mpp_do_get_boundary_ad.fh $(SRCROOT)FMS/mpp/include/mpp_do_get_boundary.fh $(SRCROOT)FMS/mpp/include/mpp_do_get_boundary_ad.fh $(SRCROOT)/FMS/mpp/include/mpp_group_update.fh $(SRCROOT)/FMS/mpp/include/group_update_pack.inc $(SRCROOT)/FMS/mpp/include/group_update_unpack.inc $(SRCROOT)FMS/mpp/include/group_update_pack.inc $(SRCROOT)FMS/mpp/include/group_update_unpack.inc $(SRCROOT)FMS/mpp/include/mpp_group_update.fh $(SRCROOT)/FMS/mpp/include/mpp_domains_reduce.inc $(SRCROOT)/FMS/mpp/include/mpp_global_reduce.fh $(SRCROOT)FMS/mpp/include/mpp_global_reduce.fh $(SRCROOT)/FMS/mpp/include/mpp_global_sum.fh $(SRCROOT)FMS/mpp/include/mpp_global_sum.fh $(SRCROOT)/FMS/mpp/include/mpp_global_sum_tl.fh $(SRCROOT)FMS/mpp/include/mpp_global_sum_tl.fh $(SRCROOT)/FMS/mpp/include/mpp_global_sum_ad.fh $(SRCROOT)FMS/mpp/include/mpp_global_sum_ad.fh $(SRCROOT)/FMS/mpp/include/mpp_global_field.fh $(SRCROOT)FMS/mpp/include/mpp_global_field.fh $(SRCROOT)/FMS/mpp/include/mpp_global_field_ad.fh $(SRCROOT)FMS/mpp/include/mpp_global_field_ad.fh $(SRCROOT)/FMS/mpp/include/mpp_do_global_field.fh $(SRCROOT)FMS/mpp/include/mpp_do_global_field.fh $(SRCROOT)/FMS/mpp/include/mpp_do_global_field_ad.fh $(SRCROOT)FMS/mpp/include/mpp_do_global_field_ad.fh $(SRCROOT)/FMS/mpp/include/mpp_unstruct_domain.inc $(SRCROOT)/FMS/mpp/include/mpp_unstruct_pass_data.fh $(SRCROOT)FMS/mpp/include/mpp_unstruct_pass_data.fh $(SRCROOT)/FMS/mpp/include/mpp_global_field_ug.fh $(SRCROOT)FMS/mpp/include/mpp_global_field_ug.fh mpp_parameter.o mpp_data.o mpp.o mpp_memutils.o mpp_efp.o platform.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)/FMS/include -I$(SRCROOT)/FMS/mpp/include -I$(SRCROOT)FMS/mpp/include	$(SRCROOT)/FMS/mpp/mpp_domains.F90
mpp_efp.o: $(SRCROOT)/FMS/mpp/mpp_efp.F90 mpp.o platform.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)/FMS/mpp/mpp_efp.F90
mpp_io.o: $(SRCROOT)/FMS/mpp/mpp_io.F90 $(SRCROOT)/FMS/include/file_version.h $(SRCROOT)/FMS/mpp/include/mpp_io_util.inc $(SRCROOT)/FMS/mpp/include/mpp_io_misc.inc $(SRCROOT)/FMS/mpp/include/mpp_io_connect.inc $(SRCROOT)/FMS/mpp/include/mpp_io_read.inc $(SRCROOT)/FMS/mpp/include/mpp_read_2Ddecomp.fh $(SRCROOT)FMS/mpp/include/mpp_read_2Ddecomp.fh $(SRCROOT)/FMS/mpp/include/mpp_read_compressed.fh $(SRCROOT)FMS/mpp/include/mpp_read_compressed.fh $(SRCROOT)/FMS/mpp/include/mpp_read_distributed_ascii.inc $(SRCROOT)/FMS/mpp/include/mpp_read_distributed_ascii.fh $(SRCROOT)FMS/mpp/include/mpp_read_distributed_ascii.fh $(SRCROOT)/FMS/mpp/include/mpp_io_write.inc $(SRCROOT)/FMS/mpp/include/mpp_write_2Ddecomp.fh $(SRCROOT)FMS/mpp/include/mpp_write_2Ddecomp.fh $(SRCROOT)/FMS/mpp/include/mpp_write_compressed.fh $(SRCROOT)FMS/mpp/include/mpp_write_compressed.fh $(SRCROOT)/FMS/mpp/include/mpp_write_unlimited_axis.fh $(SRCROOT)/FMS/mpp/include/mpp_write.fh $(SRCROOT)FMS/mpp/include/mpp_write.fh $(SRCROOT)/FMS/mpp/include/mpp_io_unstructured_write.inc $(SRCROOT)/FMS/mpp/include/mpp_io_unstructured_read.inc mpp_parameter.o mpp.o mpp_domains.o platform.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)/FMS/include -I$(SRCROOT)/FMS/mpp/include -I$(SRCROOT)FMS/mpp/include	$(SRCROOT)/FMS/mpp/mpp_io.F90
mpp_memuse.o: $(SRCROOT)/FMS/mpp/mpp_memuse.c
	$(CC) $(CPPDEFS) $(CPPFLAGS) $(CFLAGS) $(OTHERFLAGS) $(OTHER_CFLAGS) -c	$(SRCROOT)/FMS/mpp/mpp_memuse.c
mpp_memutils.o: $(SRCROOT)/FMS/mpp/mpp_memutils.F90 mpp.o platform.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)/FMS/mpp/mpp_memutils.F90
mpp_parameter.o: $(SRCROOT)/FMS/mpp/mpp_parameter.F90 $(SRCROOT)/FMS/include/file_version.h platform.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)/FMS/include	$(SRCROOT)/FMS/mpp/mpp_parameter.F90
mpp_pset.o: $(SRCROOT)/FMS/mpp/mpp_pset.F90 $(SRCROOT)/FMS/include/fms_platform.h mpp.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)/FMS/include	$(SRCROOT)/FMS/mpp/mpp_pset.F90
mpp_utilities.o: $(SRCROOT)/FMS/mpp/mpp_utilities.F90 $(SRCROOT)/FMS/include/file_version.h mpp.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)/FMS/include	$(SRCROOT)/FMS/mpp/mpp_utilities.F90
netcdf_io.o: $(SRCROOT)/FMS/fms2_io/netcdf_io.F90 $(SRCROOT)/FMS/fms2_io/include/netcdf_add_restart_variable.inc $(SRCROOT)/FMS/fms2_io/include/netcdf_read_data.inc $(SRCROOT)/FMS/fms2_io/include/netcdf_write_data.inc $(SRCROOT)/FMS/fms2_io/include/register_global_attribute.inc $(SRCROOT)/FMS/fms2_io/include/register_variable_attribute.inc $(SRCROOT)/FMS/fms2_io/include/get_global_attribute.inc $(SRCROOT)/FMS/fms2_io/include/get_variable_attribute.inc $(SRCROOT)/FMS/fms2_io/include/compressed_write.inc $(SRCROOT)/FMS/fms2_io/include/compressed_read.inc $(SRCROOT)/FMS/fms2_io/include/scatter_data_bc.inc $(SRCROOT)/FMS/fms2_io/include/gather_data_bc.inc $(SRCROOT)/FMS/fms2_io/include/unpack_data.inc mpp.o fms_io_utils.o platform.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)/FMS/fms2_io/include	$(SRCROOT)/FMS/fms2_io/netcdf_io.F90
platform.o: $(SRCROOT)/FMS/platform/platform.F90 $(SRCROOT)/FMS/include/fms_platform.h
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)/FMS/include	$(SRCROOT)/FMS/platform/platform.F90
quicksort.o: $(SRCROOT)/FMS/drifters/quicksort.F90
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)/FMS/drifters/quicksort.F90
random_numbers.o: $(SRCROOT)/FMS/random_numbers/random_numbers.F90 $(SRCROOT)/FMS/random_numbers/include/random_numbers_r4.fh $(SRCROOT)/FMS/random_numbers/include/random_numbers.inc $(SRCROOT)/FMS/random_numbers/include/random_numbers_r8.fh mersennetwister.o time_manager.o platform.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)/FMS/random_numbers/include	$(SRCROOT)/FMS/random_numbers/random_numbers.F90
read_mosaic.o: $(SRCROOT)/FMS/mosaic/read_mosaic.c $(SRCROOT)/FMS/mosaic/read_mosaic.h $(SRCROOT)/FMS/mosaic/constant.h $(SRCROOT)/FMS/mosaic/mosaic_util.h
	$(CC) $(CPPDEFS) $(CPPFLAGS) $(CFLAGS) $(OTHERFLAGS) $(OTHER_CFLAGS) -c -I$(SRCROOT)/FMS/mosaic	$(SRCROOT)/FMS/mosaic/read_mosaic.c
sat_vapor_pres.o: $(SRCROOT)/FMS/sat_vapor_pres/sat_vapor_pres.F90 $(SRCROOT)/FMS/include/file_version.h $(SRCROOT)/FMS/sat_vapor_pres/include/sat_vapor_pres_r4.fh $(SRCROOT)/FMS/sat_vapor_pres/include/sat_vapor_pres.inc $(SRCROOT)/FMS/sat_vapor_pres/include/sat_vapor_pres_r8.fh constants.o fms.o mpp.o sat_vapor_pres_k.o platform.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)/FMS/include -I$(SRCROOT)/FMS/sat_vapor_pres/include	$(SRCROOT)/FMS/sat_vapor_pres/sat_vapor_pres.F90
sat_vapor_pres_k.o: $(SRCROOT)/FMS/sat_vapor_pres/sat_vapor_pres_k.F90 $(SRCROOT)/FMS/include/file_version.h $(SRCROOT)/FMS/sat_vapor_pres/include/sat_vapor_pres_k_r4.fh $(SRCROOT)/FMS/sat_vapor_pres/include/sat_vapor_pres_k.inc $(SRCROOT)/FMS/sat_vapor_pres/include/sat_vapor_pres_k_r8.fh fms.o platform.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)/FMS/include -I$(SRCROOT)/FMS/sat_vapor_pres/include	$(SRCROOT)/FMS/sat_vapor_pres/sat_vapor_pres_k.F90
stock_constants.o: $(SRCROOT)/FMS/exchange/stock_constants.F90 $(SRCROOT)/FMS/include/file_version.h mpp.o fms.o time_manager.o diag_manager.o platform.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)/FMS/include	$(SRCROOT)/FMS/exchange/stock_constants.F90
time_interp.o: $(SRCROOT)/FMS/time_interp/time_interp.F90 $(SRCROOT)/FMS/include/file_version.h $(SRCROOT)/FMS/time_interp/include/time_interp_r4.fh $(SRCROOT)/FMS/time_interp/include/time_interp.inc $(SRCROOT)/FMS/time_interp/include/time_interp_r8.fh time_manager.o fms.o mpp.o platform.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)/FMS/include -I$(SRCROOT)/FMS/time_interp/include	$(SRCROOT)/FMS/time_interp/time_interp.F90
time_interp_external.o: $(SRCROOT)/FMS/time_interp/time_interp_external.F90 $(SRCROOT)/FMS/include/fms_platform.h $(SRCROOT)/FMS/include/file_version.h fms.o mpp.o mpp_io.o time_manager.o get_cal_time.o mpp_domains.o time_interp.o axis_utils.o platform.o horiz_interp.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)/FMS/include	$(SRCROOT)/FMS/time_interp/time_interp_external.F90
time_interp_external2.o: $(SRCROOT)/FMS/time_interp/time_interp_external2.F90 $(SRCROOT)/FMS/include/file_version.h $(SRCROOT)/FMS/time_interp/include/time_interp_external2_r4.fh $(SRCROOT)/FMS/time_interp/include/time_interp_external2.inc $(SRCROOT)/FMS/time_interp/include/time_interp_external2_r8.fh platform.o fms.o mpp.o time_manager.o get_cal_time.o mpp_domains.o time_interp.o axis_utils2.o horiz_interp.o fms2_io.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)/FMS/include -I$(SRCROOT)/FMS/time_interp/include	$(SRCROOT)/FMS/time_interp/time_interp_external2.F90
time_manager.o: $(SRCROOT)/FMS/time_manager/time_manager.F90 $(SRCROOT)/FMS/include/file_version.h platform.o constants.o fms.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)/FMS/include	$(SRCROOT)/FMS/time_manager/time_manager.F90
topography.o: $(SRCROOT)/FMS/topography/topography.F90 $(SRCROOT)/FMS/include/file_version.h $(SRCROOT)/FMS/topography/include/topography_r4.fh $(SRCROOT)/FMS/topography/include/topography.inc $(SRCROOT)/FMS/topography/include/topography_r8.fh gaussian_topog.o horiz_interp.o fms.o fms2_io.o constants.o mpp.o platform.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)/FMS/include -I$(SRCROOT)/FMS/topography/include	$(SRCROOT)/FMS/topography/topography.F90
tracer_manager.o: $(SRCROOT)/FMS/tracer_manager/tracer_manager.F90 $(SRCROOT)/FMS/include/file_version.h $(SRCROOT)/FMS/tracer_manager/include/tracer_manager_r4.fh $(SRCROOT)/FMS/tracer_manager/include/tracer_manager.inc $(SRCROOT)/FMS/tracer_manager/include/tracer_manager_r8.fh mpp.o fms.o field_manager.o platform.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)/FMS/include -I$(SRCROOT)/FMS/tracer_manager/include	$(SRCROOT)/FMS/tracer_manager/tracer_manager.F90
tridiagonal.o: $(SRCROOT)/FMS/tridiagonal/tridiagonal.F90 $(SRCROOT)/FMS/tridiagonal/include/tridiagonal_r4.fh $(SRCROOT)/FMS/tridiagonal/include/tridiagonal.inc $(SRCROOT)/FMS/tridiagonal/include/tridiagonal_r8.fh platform.o mpp.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)/FMS/tridiagonal/include	$(SRCROOT)/FMS/tridiagonal/tridiagonal.F90
xgrid.o: $(SRCROOT)/FMS/exchange/xgrid.F90 $(SRCROOT)/FMS/include/file_version.h fms.o mpp.o mpp_domains.o constants.o mosaic2.o stock_constants.o gradient.o fms2_io.o platform.o time_manager.o diag_manager.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)/FMS/include	$(SRCROOT)/FMS/exchange/xgrid.F90
yaml_output_functions.o: $(SRCROOT)/FMS/parser/yaml_output_functions.c
	$(CC) $(CPPDEFS) $(CPPFLAGS) $(CFLAGS) $(OTHERFLAGS) $(OTHER_CFLAGS) -c	$(SRCROOT)/FMS/parser/yaml_output_functions.c
yaml_parser.o: $(SRCROOT)/FMS/parser/yaml_parser.F90 fms.o fms_string_utils.o platform.o mpp.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)/FMS/parser/yaml_parser.F90
yaml_parser_binding.o: $(SRCROOT)/FMS/parser/yaml_parser_binding.c
	$(CC) $(CPPDEFS) $(CPPFLAGS) $(CFLAGS) $(OTHERFLAGS) $(OTHER_CFLAGS) -c	$(SRCROOT)/FMS/parser/yaml_parser_binding.c
./mpp_read_compressed.fh: $(SRCROOT)FMS/mpp/include/mpp_read_compressed.fh
	cp $(SRCROOT)FMS/mpp/include/mpp_read_compressed.fh .
./mpp_domains_define.inc: $(SRCROOT)/FMS/mpp/include/mpp_domains_define.inc
	cp $(SRCROOT)/FMS/mpp/include/mpp_domains_define.inc .
./drifters_compute_k.fh: $(SRCROOT)/FMS/drifters/drifters_compute_k.fh
	cp $(SRCROOT)/FMS/drifters/drifters_compute_k.fh .
./mpp_do_global_field.fh: $(SRCROOT)/FMS/mpp/include/mpp_do_global_field.fh
	cp $(SRCROOT)/FMS/mpp/include/mpp_do_global_field.fh .
./grid.F90: $(SRCROOT)/FMS/mosaic/grid.F90
	cp $(SRCROOT)/FMS/mosaic/grid.F90 .
./fms_diag_bbox.F90: $(SRCROOT)/FMS/diag_manager/fms_diag_bbox.F90
	cp $(SRCROOT)/FMS/diag_manager/fms_diag_bbox.F90 .
./mpp_read_distributed_ascii.fh: $(SRCROOT)/FMS/mpp/include/mpp_read_distributed_ascii.fh
	cp $(SRCROOT)/FMS/mpp/include/mpp_read_distributed_ascii.fh .
./mpp_do_update_nest.fh: $(SRCROOT)/FMS/mpp/include/mpp_do_update_nest.fh
	cp $(SRCROOT)/FMS/mpp/include/mpp_do_update_nest.fh .
./amip_interp.F90: $(SRCROOT)/FMS/amip_interp/amip_interp.F90
	cp $(SRCROOT)/FMS/amip_interp/amip_interp.F90 .
./diag_output.F90: $(SRCROOT)/FMS/diag_manager/diag_output.F90
	cp $(SRCROOT)/FMS/diag_manager/diag_output.F90 .
./mpp_sum_mpi_ad.fh: $(SRCROOT)/FMS/mpp/include/mpp_sum_mpi_ad.fh
	cp $(SRCROOT)/FMS/mpp/include/mpp_sum_mpi_ad.fh .
./sat_vapor_pres_r4.fh: $(SRCROOT)/FMS/sat_vapor_pres/include/sat_vapor_pres_r4.fh
	cp $(SRCROOT)/FMS/sat_vapor_pres/include/sat_vapor_pres_r4.fh .
./system_clock.fh: $(SRCROOT)/FMS/mpp/include/system_clock.fh
	cp $(SRCROOT)/FMS/mpp/include/system_clock.fh .
./fms_diag_fieldbuff_update.F90: $(SRCROOT)/FMS/diag_manager/fms_diag_fieldbuff_update.F90
	cp $(SRCROOT)/FMS/diag_manager/fms_diag_fieldbuff_update.F90 .
./time_interp.inc: $(SRCROOT)/FMS/time_interp/include/time_interp.inc
	cp $(SRCROOT)/FMS/time_interp/include/time_interp.inc .
./mpp_domains_reduce.inc: $(SRCROOT)/FMS/mpp/include/mpp_domains_reduce.inc
	cp $(SRCROOT)/FMS/mpp/include/mpp_domains_reduce.inc .
./diag_data.F90: $(SRCROOT)/FMS/diag_manager/diag_data.F90
	cp $(SRCROOT)/FMS/diag_manager/diag_data.F90 .
./diag_integral_r8.fh: $(SRCROOT)/FMS/diag_integral/include/diag_integral_r8.fh
	cp $(SRCROOT)/FMS/diag_integral/include/diag_integral_r8.fh .
./get_grid_version_r8.fh: $(SRCROOT)/FMS/data_override/include/get_grid_version_r8.fh
	cp $(SRCROOT)/FMS/data_override/include/get_grid_version_r8.fh .
./parse.inc: $(SRCROOT)/FMS/field_manager/parse.inc
	cp $(SRCROOT)/FMS/field_manager/parse.inc .
./horiz_interp_bicubic_r4.fh: $(SRCROOT)/FMS/horiz_interp/include/horiz_interp_bicubic_r4.fh
	cp $(SRCROOT)/FMS/horiz_interp/include/horiz_interp_bicubic_r4.fh .
./mpp_chksum_scalar.fh: $(SRCROOT)FMS/mpp/include/mpp_chksum_scalar.fh
	cp $(SRCROOT)FMS/mpp/include/mpp_chksum_scalar.fh .
./monin_obukhov_inter_r8.fh: $(SRCROOT)/FMS/monin_obukhov/include/monin_obukhov_inter_r8.fh
	cp $(SRCROOT)/FMS/monin_obukhov/include/monin_obukhov_inter_r8.fh .
./amip_interp.inc: $(SRCROOT)/FMS/amip_interp/include/amip_interp.inc
	cp $(SRCROOT)/FMS/amip_interp/include/amip_interp.inc .
./mpp_global_sum.fh: $(SRCROOT)/FMS/mpp/include/mpp_global_sum.fh
	cp $(SRCROOT)/FMS/mpp/include/mpp_global_sum.fh .
./horiz_interp_spherical_r4.fh: $(SRCROOT)/FMS/horiz_interp/include/horiz_interp_spherical_r4.fh
	cp $(SRCROOT)/FMS/horiz_interp/include/horiz_interp_spherical_r4.fh .
./mpp_get_boundary_ad.fh: $(SRCROOT)/FMS/mpp/include/mpp_get_boundary_ad.fh
	cp $(SRCROOT)/FMS/mpp/include/mpp_get_boundary_ad.fh .
./column_diagnostics_r8.fh: $(SRCROOT)/FMS/column_diagnostics/include/column_diagnostics_r8.fh
	cp $(SRCROOT)/FMS/column_diagnostics/include/column_diagnostics_r8.fh .
./fms_yaml_output.F90: $(SRCROOT)/FMS/parser/fms_yaml_output.F90
	cp $(SRCROOT)/FMS/parser/fms_yaml_output.F90 .
./fms_diag_time_reduction.F90: $(SRCROOT)/FMS/diag_manager/fms_diag_time_reduction.F90
	cp $(SRCROOT)/FMS/diag_manager/fms_diag_time_reduction.F90 .
./mpp_scatter.fh: $(SRCROOT)FMS/mpp/include/mpp_scatter.fh
	cp $(SRCROOT)FMS/mpp/include/mpp_scatter.fh .
./memutils.F90: $(SRCROOT)/FMS/memutils/memutils.F90
	cp $(SRCROOT)/FMS/memutils/memutils.F90 .
./mpp_comm.inc: $(SRCROOT)/FMS/mpp/include/mpp_comm.inc
	cp $(SRCROOT)/FMS/mpp/include/mpp_comm.inc .
./diag_util.F90: $(SRCROOT)/FMS/diag_manager/diag_util.F90
	cp $(SRCROOT)/FMS/diag_manager/diag_util.F90 .
./data_override_r4.fh: $(SRCROOT)/FMS/data_override/include/data_override_r4.fh
	cp $(SRCROOT)/FMS/data_override/include/data_override_r4.fh .
./mpp_do_update_nonblock.fh: $(SRCROOT)/FMS/mpp/include/mpp_do_update_nonblock.fh
	cp $(SRCROOT)/FMS/mpp/include/mpp_do_update_nonblock.fh .
./get_global_attribute.inc: $(SRCROOT)/FMS/fms2_io/include/get_global_attribute.inc
	cp $(SRCROOT)/FMS/fms2_io/include/get_global_attribute.inc .
./time_interp_external2.inc: $(SRCROOT)/FMS/time_interp/include/time_interp_external2.inc
	cp $(SRCROOT)/FMS/time_interp/include/time_interp_external2.inc .
./mpp_update_domains2D_nonblock.fh: $(SRCROOT)FMS/mpp/include/mpp_update_domains2D_nonblock.fh
	cp $(SRCROOT)FMS/mpp/include/mpp_update_domains2D_nonblock.fh .
./data_override.inc: $(SRCROOT)/FMS/data_override/include/data_override.inc
	cp $(SRCROOT)/FMS/data_override/include/data_override.inc .
./fms_io_unstructured_get_file_unit.inc: $(SRCROOT)/FMS/fms/fms_io_unstructured_get_file_unit.inc
	cp $(SRCROOT)/FMS/fms/fms_io_unstructured_get_file_unit.inc .
./mpp_do_global_field_ad.fh: $(SRCROOT)FMS/mpp/include/mpp_do_global_field_ad.fh
	cp $(SRCROOT)FMS/mpp/include/mpp_do_global_field_ad.fh .
./fms_io_unstructured_get_file_name.inc: $(SRCROOT)/FMS/fms/fms_io_unstructured_get_file_name.inc
	cp $(SRCROOT)/FMS/fms/fms_io_unstructured_get_file_name.inc .
./mpp_transmit_nocomm.fh: $(SRCROOT)FMS/mpp/include/mpp_transmit_nocomm.fh
	cp $(SRCROOT)FMS/mpp/include/mpp_transmit_nocomm.fh .
./mpp_update_domains2D.fh: $(SRCROOT)FMS/mpp/include/mpp_update_domains2D.fh
	cp $(SRCROOT)FMS/mpp/include/mpp_update_domains2D.fh .
./mpp_pset.F90: $(SRCROOT)/FMS/mpp/mpp_pset.F90
	cp $(SRCROOT)/FMS/mpp/mpp_pset.F90 .
./sat_vapor_pres_k.inc: $(SRCROOT)/FMS/sat_vapor_pres/include/sat_vapor_pres_k.inc
	cp $(SRCROOT)/FMS/sat_vapor_pres/include/sat_vapor_pres_k.inc .
./mosaic.F90: $(SRCROOT)/FMS/mosaic/mosaic.F90
	cp $(SRCROOT)/FMS/mosaic/mosaic.F90 .
./constantsr4.F90: $(SRCROOT)/FMS/constants4/constantsr4.F90
	cp $(SRCROOT)/FMS/constants4/constantsr4.F90 .
./gradient_c2l.h: $(SRCROOT)/FMS/mosaic/gradient_c2l.h
	cp $(SRCROOT)/FMS/mosaic/gradient_c2l.h .
./mpp_unstruct_pass_data.fh: $(SRCROOT)/FMS/mpp/include/mpp_unstruct_pass_data.fh
	cp $(SRCROOT)/FMS/mpp/include/mpp_unstruct_pass_data.fh .
./fm_util.inc: $(SRCROOT)/FMS/field_manager/include/fm_util.inc
	cp $(SRCROOT)/FMS/field_manager/include/fm_util.inc .
./mpp_write_compressed.fh: $(SRCROOT)/FMS/mpp/include/mpp_write_compressed.fh
	cp $(SRCROOT)/FMS/mpp/include/mpp_write_compressed.fh .
./mpp_reduce_mpi.fh: $(SRCROOT)FMS/mpp/include/mpp_reduce_mpi.fh
	cp $(SRCROOT)FMS/mpp/include/mpp_reduce_mpi.fh .
./mpp_reduce_nocomm.fh: $(SRCROOT)FMS/mpp/include/mpp_reduce_nocomm.fh
	cp $(SRCROOT)FMS/mpp/include/mpp_reduce_nocomm.fh .
./mpp_io_unstructured_write.inc: $(SRCROOT)/FMS/mpp/include/mpp_io_unstructured_write.inc
	cp $(SRCROOT)/FMS/mpp/include/mpp_io_unstructured_write.inc .
./mpp_chksum.fh: $(SRCROOT)/FMS/mpp/include/mpp_chksum.fh
	cp $(SRCROOT)/FMS/mpp/include/mpp_chksum.fh .
./constant.h: $(SRCROOT)/FMS/mosaic/constant.h
	cp $(SRCROOT)/FMS/mosaic/constant.h .
./scatter_data_bc.inc: $(SRCROOT)/FMS/fms2_io/include/scatter_data_bc.inc
	cp $(SRCROOT)/FMS/fms2_io/include/scatter_data_bc.inc .
./drifters_comm.F90: $(SRCROOT)/FMS/drifters/drifters_comm.F90
	cp $(SRCROOT)/FMS/drifters/drifters_comm.F90 .
./mpp_global_sum_ad.fh: $(SRCROOT)/FMS/mpp/include/mpp_global_sum_ad.fh
	cp $(SRCROOT)/FMS/mpp/include/mpp_global_sum_ad.fh .
./stock_constants.F90: $(SRCROOT)/FMS/exchange/stock_constants.F90
	cp $(SRCROOT)/FMS/exchange/stock_constants.F90 .
./astronomy.inc: $(SRCROOT)/FMS/astronomy/include/astronomy.inc
	cp $(SRCROOT)/FMS/astronomy/include/astronomy.inc .
./mpp_error_s_s.fh: $(SRCROOT)/FMS/mpp/include/mpp_error_s_s.fh
	cp $(SRCROOT)/FMS/mpp/include/mpp_error_s_s.fh .
./mpp_global_reduce.fh: $(SRCROOT)/FMS/mpp/include/mpp_global_reduce.fh
	cp $(SRCROOT)/FMS/mpp/include/mpp_global_reduce.fh .
./mersennetwister.F90: $(SRCROOT)/FMS/random_numbers/mersennetwister.F90
	cp $(SRCROOT)/FMS/random_numbers/mersennetwister.F90 .
./register_domain_restart_variable.inc: $(SRCROOT)/FMS/fms2_io/include/register_domain_restart_variable.inc
	cp $(SRCROOT)/FMS/fms2_io/include/register_domain_restart_variable.inc .
./column_diagnostics_r4.fh: $(SRCROOT)/FMS/column_diagnostics/include/column_diagnostics_r4.fh
	cp $(SRCROOT)/FMS/column_diagnostics/include/column_diagnostics_r4.fh .
./read_mosaic.h: $(SRCROOT)/FMS/mosaic/read_mosaic.h
	cp $(SRCROOT)/FMS/mosaic/read_mosaic.h .
./grid2_r4.fh: $(SRCROOT)/FMS/mosaic2/include/grid2_r4.fh
	cp $(SRCROOT)/FMS/mosaic2/include/grid2_r4.fh .
./mpp_write_2Ddecomp.fh: $(SRCROOT)FMS/mpp/include/mpp_write_2Ddecomp.fh
	cp $(SRCROOT)FMS/mpp/include/mpp_write_2Ddecomp.fh .
./fms_io_unstructured_file_unit.inc: $(SRCROOT)/FMS/fms/fms_io_unstructured_file_unit.inc
	cp $(SRCROOT)/FMS/fms/fms_io_unstructured_file_unit.inc .
./register_variable_attribute.inc: $(SRCROOT)/FMS/fms2_io/include/register_variable_attribute.inc
	cp $(SRCROOT)/FMS/fms2_io/include/register_variable_attribute.inc .
./mpp_global_sum_tl.fh: $(SRCROOT)/FMS/mpp/include/mpp_global_sum_tl.fh
	cp $(SRCROOT)/FMS/mpp/include/mpp_global_sum_tl.fh .
./array_utils_char.inc: $(SRCROOT)/FMS/fms2_io/include/array_utils_char.inc
	cp $(SRCROOT)/FMS/fms2_io/include/array_utils_char.inc .
./tracer_manager_r4.fh: $(SRCROOT)/FMS/tracer_manager/include/tracer_manager_r4.fh
	cp $(SRCROOT)/FMS/tracer_manager/include/tracer_manager_r4.fh .
./mpp_global_field.fh: $(SRCROOT)/FMS/mpp/include/mpp_global_field.fh
	cp $(SRCROOT)/FMS/mpp/include/mpp_global_field.fh .
./mpp_unstruct_pass_data.fh: $(SRCROOT)FMS/mpp/include/mpp_unstruct_pass_data.fh
	cp $(SRCROOT)FMS/mpp/include/mpp_unstruct_pass_data.fh .
./time_interp_external2_r8.fh: $(SRCROOT)/FMS/time_interp/include/time_interp_external2_r8.fh
	cp $(SRCROOT)/FMS/time_interp/include/time_interp_external2_r8.fh .
./unstructured_domain_read.inc: $(SRCROOT)/FMS/fms2_io/include/unstructured_domain_read.inc
	cp $(SRCROOT)/FMS/fms2_io/include/unstructured_domain_read.inc .
./fmsconstantsr4.F90: $(SRCROOT)/FMS/constants4/fmsconstantsr4.F90
	cp $(SRCROOT)/FMS/constants4/fmsconstantsr4.F90 .
./mpp_memuse.c: $(SRCROOT)/FMS/mpp/mpp_memuse.c
	cp $(SRCROOT)/FMS/mpp/mpp_memuse.c .
./mpp_read_distributed_ascii.inc: $(SRCROOT)/FMS/mpp/include/mpp_read_distributed_ascii.inc
	cp $(SRCROOT)/FMS/mpp/include/mpp_read_distributed_ascii.inc .
./mpp_sum_ad.inc: $(SRCROOT)/FMS/mpp/include/mpp_sum_ad.inc
	cp $(SRCROOT)/FMS/mpp/include/mpp_sum_ad.inc .
./drifters_input.F90: $(SRCROOT)/FMS/drifters/drifters_input.F90
	cp $(SRCROOT)/FMS/drifters/drifters_input.F90 .
./mpp_global_sum_ad.fh: $(SRCROOT)FMS/mpp/include/mpp_global_sum_ad.fh
	cp $(SRCROOT)FMS/mpp/include/mpp_global_sum_ad.fh .
./mpp_do_updateV.fh: $(SRCROOT)/FMS/mpp/include/mpp_do_updateV.fh
	cp $(SRCROOT)/FMS/mpp/include/mpp_do_updateV.fh .
./mpp_do_updateV_nonblock.fh: $(SRCROOT)FMS/mpp/include/mpp_do_updateV_nonblock.fh
	cp $(SRCROOT)FMS/mpp/include/mpp_do_updateV_nonblock.fh .
./astronomy_r4.fh: $(SRCROOT)/FMS/astronomy/include/astronomy_r4.fh
	cp $(SRCROOT)/FMS/astronomy/include/astronomy_r4.fh .
./field_manager.inc: $(SRCROOT)/FMS/field_manager/include/field_manager.inc
	cp $(SRCROOT)/FMS/field_manager/include/field_manager.inc .
./read_mosaic.c: $(SRCROOT)/FMS/mosaic/read_mosaic.c
	cp $(SRCROOT)/FMS/mosaic/read_mosaic.c .
./mpp_type_mpi.fh: $(SRCROOT)/FMS/mpp/include/mpp_type_mpi.fh
	cp $(SRCROOT)/FMS/mpp/include/mpp_type_mpi.fh .
./field_manager_r8.fh: $(SRCROOT)/FMS/field_manager/include/field_manager_r8.fh
	cp $(SRCROOT)/FMS/field_manager/include/field_manager_r8.fh .
./mpp_group_update.fh: $(SRCROOT)FMS/mpp/include/mpp_group_update.fh
	cp $(SRCROOT)FMS/mpp/include/mpp_group_update.fh .
./mpp_do_redistribute.fh: $(SRCROOT)/FMS/mpp/include/mpp_do_redistribute.fh
	cp $(SRCROOT)/FMS/mpp/include/mpp_do_redistribute.fh .
./mpp_io_unstructured_read.inc: $(SRCROOT)/FMS/mpp/include/mpp_io_unstructured_read.inc
	cp $(SRCROOT)/FMS/mpp/include/mpp_io_unstructured_read.inc .
./amip_interp_r8.fh: $(SRCROOT)/FMS/amip_interp/include/amip_interp_r8.fh
	cp $(SRCROOT)/FMS/amip_interp/include/amip_interp_r8.fh .
./horiz_interp_conserve.inc: $(SRCROOT)/FMS/horiz_interp/include/horiz_interp_conserve.inc
	cp $(SRCROOT)/FMS/horiz_interp/include/horiz_interp_conserve.inc .
./topography_r8.fh: $(SRCROOT)/FMS/topography/include/topography_r8.fh
	cp $(SRCROOT)/FMS/topography/include/topography_r8.fh .
./mpp_do_update_ad.fh: $(SRCROOT)/FMS/mpp/include/mpp_do_update_ad.fh
	cp $(SRCROOT)/FMS/mpp/include/mpp_do_update_ad.fh .
./constants.F90: $(SRCROOT)/FMS/constants/constants.F90
	cp $(SRCROOT)/FMS/constants/constants.F90 .
./horiz_interp_conserve.F90: $(SRCROOT)/FMS/horiz_interp/horiz_interp_conserve.F90
	cp $(SRCROOT)/FMS/horiz_interp/horiz_interp_conserve.F90 .
./fms_diag_outfield.F90: $(SRCROOT)/FMS/diag_manager/fms_diag_outfield.F90
	cp $(SRCROOT)/FMS/diag_manager/fms_diag_outfield.F90 .
./horiz_interp_type_r4.fh: $(SRCROOT)/FMS/horiz_interp/include/horiz_interp_type_r4.fh
	cp $(SRCROOT)/FMS/horiz_interp/include/horiz_interp_type_r4.fh .
./fms_string_utils.inc: $(SRCROOT)/FMS/string_utils/include/fms_string_utils.inc
	cp $(SRCROOT)/FMS/string_utils/include/fms_string_utils.inc .
./mpp_do_updateV_ad.fh: $(SRCROOT)FMS/mpp/include/mpp_do_updateV_ad.fh
	cp $(SRCROOT)FMS/mpp/include/mpp_do_updateV_ad.fh .
./horiz_interp_type.F90: $(SRCROOT)/FMS/horiz_interp/horiz_interp_type.F90
	cp $(SRCROOT)/FMS/horiz_interp/horiz_interp_type.F90 .
./mpp_do_get_boundary.fh: $(SRCROOT)/FMS/mpp/include/mpp_do_get_boundary.fh
	cp $(SRCROOT)/FMS/mpp/include/mpp_do_get_boundary.fh .
./sat_vapor_pres.F90: $(SRCROOT)/FMS/sat_vapor_pres/sat_vapor_pres.F90
	cp $(SRCROOT)/FMS/sat_vapor_pres/sat_vapor_pres.F90 .
./mpp_util.inc: $(SRCROOT)/FMS/mpp/include/mpp_util.inc
	cp $(SRCROOT)/FMS/mpp/include/mpp_util.inc .
./fms_io_unstructured_read.inc: $(SRCROOT)/FMS/fms/fms_io_unstructured_read.inc
	cp $(SRCROOT)/FMS/fms/fms_io_unstructured_read.inc .
./mpp_global_field.fh: $(SRCROOT)FMS/mpp/include/mpp_global_field.fh
	cp $(SRCROOT)FMS/mpp/include/mpp_global_field.fh .
./astronomy.F90: $(SRCROOT)/FMS/astronomy/astronomy.F90
	cp $(SRCROOT)/FMS/astronomy/astronomy.F90 .
./mpp_global_field_ad.fh: $(SRCROOT)/FMS/mpp/include/mpp_global_field_ad.fh
	cp $(SRCROOT)/FMS/mpp/include/mpp_global_field_ad.fh .
./time_interp.F90: $(SRCROOT)/FMS/time_interp/time_interp.F90
	cp $(SRCROOT)/FMS/time_interp/time_interp.F90 .
./drifters_push.fh: $(SRCROOT)/FMS/drifters/drifters_push.fh
	cp $(SRCROOT)/FMS/drifters/drifters_push.fh .
./topography_r4.fh: $(SRCROOT)/FMS/topography/include/topography_r4.fh
	cp $(SRCROOT)/FMS/topography/include/topography_r4.fh .
./interp.c: $(SRCROOT)/FMS/mosaic/interp.c
	cp $(SRCROOT)/FMS/mosaic/interp.c .
./random_numbers.inc: $(SRCROOT)/FMS/random_numbers/include/random_numbers.inc
	cp $(SRCROOT)/FMS/random_numbers/include/random_numbers.inc .
./domain_write.inc: $(SRCROOT)/FMS/fms2_io/include/domain_write.inc
	cp $(SRCROOT)/FMS/fms2_io/include/domain_write.inc .
./compute_global_checksum.inc: $(SRCROOT)/FMS/fms2_io/include/compute_global_checksum.inc
	cp $(SRCROOT)/FMS/fms2_io/include/compute_global_checksum.inc .
./mpp_read_distributed_ascii.fh: $(SRCROOT)FMS/mpp/include/mpp_read_distributed_ascii.fh
	cp $(SRCROOT)FMS/mpp/include/mpp_read_distributed_ascii.fh .
./fms_diag_fieldbuff_update.fh: $(SRCROOT)/FMS/diag_manager/include/fms_diag_fieldbuff_update.fh
	cp $(SRCROOT)/FMS/diag_manager/include/fms_diag_fieldbuff_update.fh .
./quicksort.F90: $(SRCROOT)/FMS/drifters/quicksort.F90
	cp $(SRCROOT)/FMS/drifters/quicksort.F90 .
./data_override_r8.fh: $(SRCROOT)/FMS/data_override/include/data_override_r8.fh
	cp $(SRCROOT)/FMS/data_override/include/data_override_r8.fh .
./mpp_sum_mpi.fh: $(SRCROOT)FMS/mpp/include/mpp_sum_mpi.fh
	cp $(SRCROOT)FMS/mpp/include/mpp_sum_mpi.fh .
./random_numbers_r4.fh: $(SRCROOT)/FMS/random_numbers/include/random_numbers_r4.fh
	cp $(SRCROOT)/FMS/random_numbers/include/random_numbers_r4.fh .
./mosaic_util.c: $(SRCROOT)/FMS/mosaic/mosaic_util.c
	cp $(SRCROOT)/FMS/mosaic/mosaic_util.c .
./mosaic2.F90: $(SRCROOT)/FMS/mosaic2/mosaic2.F90
	cp $(SRCROOT)/FMS/mosaic2/mosaic2.F90 .
./mpp_global_sum.fh: $(SRCROOT)FMS/mpp/include/mpp_global_sum.fh
	cp $(SRCROOT)FMS/mpp/include/mpp_global_sum.fh .
./mpp_do_updateV_nonblock.fh: $(SRCROOT)/FMS/mpp/include/mpp_do_updateV_nonblock.fh
	cp $(SRCROOT)/FMS/mpp/include/mpp_do_updateV_nonblock.fh .
./mpp_scatter.fh: $(SRCROOT)/FMS/mpp/include/mpp_scatter.fh
	cp $(SRCROOT)/FMS/mpp/include/mpp_scatter.fh .
./fms_io_unstructured_get_field_size.inc: $(SRCROOT)/FMS/fms/fms_io_unstructured_get_field_size.inc
	cp $(SRCROOT)/FMS/fms/fms_io_unstructured_get_field_size.inc .
./gfdl_constants.fh: $(SRCROOT)/FMS/constants/gfdl_constants.fh
	cp $(SRCROOT)/FMS/constants/gfdl_constants.fh .
./get_cal_time.F90: $(SRCROOT)/FMS/time_manager/get_cal_time.F90
	cp $(SRCROOT)/FMS/time_manager/get_cal_time.F90 .
./mpp_do_redistribute.fh: $(SRCROOT)FMS/mpp/include/mpp_do_redistribute.fh
	cp $(SRCROOT)FMS/mpp/include/mpp_do_redistribute.fh .
./fms_diag_elem_weight_procs.F90: $(SRCROOT)/FMS/diag_manager/fms_diag_elem_weight_procs.F90
	cp $(SRCROOT)/FMS/diag_manager/fms_diag_elem_weight_procs.F90 .
./coupler_types.F90: $(SRCROOT)/FMS/coupler/coupler_types.F90
	cp $(SRCROOT)/FMS/coupler/coupler_types.F90 .
./read_data_2d.inc: $(SRCROOT)/FMS/fms/read_data_2d.inc
	cp $(SRCROOT)/FMS/fms/read_data_2d.inc .
./fms_affinity.F90: $(SRCROOT)/FMS/affinity/fms_affinity.F90
	cp $(SRCROOT)/FMS/affinity/fms_affinity.F90 .
./fms_r4.fh: $(SRCROOT)/FMS/fms/include/fms_r4.fh
	cp $(SRCROOT)/FMS/fms/include/fms_r4.fh .
./write_data.inc: $(SRCROOT)/FMS/fms/write_data.inc
	cp $(SRCROOT)/FMS/fms/write_data.inc .
./mpp_alltoall_mpi.fh: $(SRCROOT)/FMS/mpp/include/mpp_alltoall_mpi.fh
	cp $(SRCROOT)/FMS/mpp/include/mpp_alltoall_mpi.fh .
./fms_netcdf_domain_io.F90: $(SRCROOT)/FMS/fms2_io/fms_netcdf_domain_io.F90
	cp $(SRCROOT)/FMS/fms2_io/fms_netcdf_domain_io.F90 .
./yaml_output_functions.c: $(SRCROOT)/FMS/parser/yaml_output_functions.c
	cp $(SRCROOT)/FMS/parser/yaml_output_functions.c .
./fms_platform.h: $(SRCROOT)/FMS/include/fms_platform.h
	cp $(SRCROOT)/FMS/include/fms_platform.h .
./mpp_define_nest_domains.inc: $(SRCROOT)/FMS/mpp/include/mpp_define_nest_domains.inc
	cp $(SRCROOT)/FMS/mpp/include/mpp_define_nest_domains.inc .
./topography.inc: $(SRCROOT)/FMS/topography/include/topography.inc
	cp $(SRCROOT)/FMS/topography/include/topography.inc .
./file_version.h: $(SRCROOT)FMS/include/file_version.h
	cp $(SRCROOT)FMS/include/file_version.h .
./mpp_sum_nocomm.fh: $(SRCROOT)FMS/mpp/include/mpp_sum_nocomm.fh
	cp $(SRCROOT)FMS/mpp/include/mpp_sum_nocomm.fh .
./monin_obukhov_r4.fh: $(SRCROOT)/FMS/monin_obukhov/include/monin_obukhov_r4.fh
	cp $(SRCROOT)/FMS/monin_obukhov/include/monin_obukhov_r4.fh .
./grid2.inc: $(SRCROOT)/FMS/mosaic2/include/grid2.inc
	cp $(SRCROOT)/FMS/mosaic2/include/grid2.inc .
./mpp_write.fh: $(SRCROOT)/FMS/mpp/include/mpp_write.fh
	cp $(SRCROOT)/FMS/mpp/include/mpp_write.fh .
./group_update_pack.inc: $(SRCROOT)FMS/mpp/include/group_update_pack.inc
	cp $(SRCROOT)FMS/mpp/include/group_update_pack.inc .
./diag_manager.F90: $(SRCROOT)/FMS/diag_manager/diag_manager.F90
	cp $(SRCROOT)/FMS/diag_manager/diag_manager.F90 .
./affinity.c: $(SRCROOT)/FMS/affinity/affinity.c
	cp $(SRCROOT)/FMS/affinity/affinity.c .
./mpp_update_domains2D_ad.fh: $(SRCROOT)/FMS/mpp/include/mpp_update_domains2D_ad.fh
	cp $(SRCROOT)/FMS/mpp/include/mpp_update_domains2D_ad.fh .
./mosaic2_r4.fh: $(SRCROOT)/FMS/mosaic2/include/mosaic2_r4.fh
	cp $(SRCROOT)/FMS/mosaic2/include/mosaic2_r4.fh .
./mpp_global_field_ug.fh: $(SRCROOT)FMS/mpp/include/mpp_global_field_ug.fh
	cp $(SRCROOT)FMS/mpp/include/mpp_global_field_ug.fh .
./tridiagonal.inc: $(SRCROOT)/FMS/tridiagonal/include/tridiagonal.inc
	cp $(SRCROOT)/FMS/tridiagonal/include/tridiagonal.inc .
./fms_string_utils.F90: $(SRCROOT)/FMS/string_utils/fms_string_utils.F90
	cp $(SRCROOT)/FMS/string_utils/fms_string_utils.F90 .
./time_interp_external2.F90: $(SRCROOT)/FMS/time_interp/time_interp_external2.F90
	cp $(SRCROOT)/FMS/time_interp/time_interp_external2.F90 .
./mpp_read_compressed.fh: $(SRCROOT)/FMS/mpp/include/mpp_read_compressed.fh
	cp $(SRCROOT)/FMS/mpp/include/mpp_read_compressed.fh .
./mpp_error_a_a.fh: $(SRCROOT)/FMS/mpp/include/mpp_error_a_a.fh
	cp $(SRCROOT)/FMS/mpp/include/mpp_error_a_a.fh .
./mpp_get_boundary.fh: $(SRCROOT)FMS/mpp/include/mpp_get_boundary.fh
	cp $(SRCROOT)FMS/mpp/include/mpp_get_boundary.fh .
./group_update_pack.inc: $(SRCROOT)/FMS/mpp/include/group_update_pack.inc
	cp $(SRCROOT)/FMS/mpp/include/group_update_pack.inc .
./mpp_transmit_mpi.fh: $(SRCROOT)/FMS/mpp/include/mpp_transmit_mpi.fh
	cp $(SRCROOT)/FMS/mpp/include/mpp_transmit_mpi.fh .
./drifters.F90: $(SRCROOT)/FMS/drifters/drifters.F90
	cp $(SRCROOT)/FMS/drifters/drifters.F90 .
./fms_string_utils_r4.fh: $(SRCROOT)/FMS/string_utils/include/fms_string_utils_r4.fh
	cp $(SRCROOT)/FMS/string_utils/include/fms_string_utils_r4.fh .
./mpp_data_nocomm.inc: $(SRCROOT)/FMS/mpp/include/mpp_data_nocomm.inc
	cp $(SRCROOT)/FMS/mpp/include/mpp_data_nocomm.inc .
./horiz_interp_spherical_r8.fh: $(SRCROOT)/FMS/horiz_interp/include/horiz_interp_spherical_r8.fh
	cp $(SRCROOT)/FMS/horiz_interp/include/horiz_interp_spherical_r8.fh .
./mpp_do_update.fh: $(SRCROOT)FMS/mpp/include/mpp_do_update.fh
	cp $(SRCROOT)FMS/mpp/include/mpp_do_update.fh .
./field_manager_r4.fh: $(SRCROOT)/FMS/field_manager/include/field_manager_r4.fh
	cp $(SRCROOT)/FMS/field_manager/include/field_manager_r4.fh .
./axis_utils2.inc: $(SRCROOT)/FMS/axis_utils/include/axis_utils2.inc
	cp $(SRCROOT)/FMS/axis_utils/include/axis_utils2.inc .
./mpp.F90: $(SRCROOT)/FMS/mpp/mpp.F90
	cp $(SRCROOT)/FMS/mpp/mpp.F90 .
./horiz_interp_r8.fh: $(SRCROOT)/FMS/horiz_interp/include/horiz_interp_r8.fh
	cp $(SRCROOT)/FMS/horiz_interp/include/horiz_interp_r8.fh .
./time_interp_r8.fh: $(SRCROOT)/FMS/time_interp/include/time_interp_r8.fh
	cp $(SRCROOT)/FMS/time_interp/include/time_interp_r8.fh .
./time_interp_external2_r4.fh: $(SRCROOT)/FMS/time_interp/include/time_interp_external2_r4.fh
	cp $(SRCROOT)/FMS/time_interp/include/time_interp_external2_r4.fh .
./tridiagonal_r4.fh: $(SRCROOT)/FMS/tridiagonal/include/tridiagonal_r4.fh
	cp $(SRCROOT)/FMS/tridiagonal/include/tridiagonal_r4.fh .
./mpp_do_get_boundary_ad.fh: $(SRCROOT)FMS/mpp/include/mpp_do_get_boundary_ad.fh
	cp $(SRCROOT)FMS/mpp/include/mpp_do_get_boundary_ad.fh .
./mpp_sum.inc: $(SRCROOT)FMS/mpp/include/mpp_sum.inc
	cp $(SRCROOT)FMS/mpp/include/mpp_sum.inc .
./fms_io_unstructured_register_restart_field.inc: $(SRCROOT)/FMS/fms/fms_io_unstructured_register_restart_field.inc
	cp $(SRCROOT)/FMS/fms/fms_io_unstructured_register_restart_field.inc .
./mpp_global_sum_tl.fh: $(SRCROOT)FMS/mpp/include/mpp_global_sum_tl.fh
	cp $(SRCROOT)FMS/mpp/include/mpp_global_sum_tl.fh .
./mpp_gather.fh: $(SRCROOT)FMS/mpp/include/mpp_gather.fh
	cp $(SRCROOT)FMS/mpp/include/mpp_gather.fh .
./mpp_sum_nocomm_ad.fh: $(SRCROOT)/FMS/mpp/include/mpp_sum_nocomm_ad.fh
	cp $(SRCROOT)/FMS/mpp/include/mpp_sum_nocomm_ad.fh .
./mpp_transmit.inc: $(SRCROOT)FMS/mpp/include/mpp_transmit.inc
	cp $(SRCROOT)FMS/mpp/include/mpp_transmit.inc .
./mpp_do_get_boundary_ad.fh: $(SRCROOT)/FMS/mpp/include/mpp_do_get_boundary_ad.fh
	cp $(SRCROOT)/FMS/mpp/include/mpp_do_get_boundary_ad.fh .
./netcdf_read_data.inc: $(SRCROOT)/FMS/fms2_io/include/netcdf_read_data.inc
	cp $(SRCROOT)/FMS/fms2_io/include/netcdf_read_data.inc .
./horiz_interp_type_r8.fh: $(SRCROOT)/FMS/horiz_interp/include/horiz_interp_type_r8.fh
	cp $(SRCROOT)/FMS/horiz_interp/include/horiz_interp_type_r8.fh .
./mpp_sum_nocomm_ad.fh: $(SRCROOT)FMS/mpp/include/mpp_sum_nocomm_ad.fh
	cp $(SRCROOT)FMS/mpp/include/mpp_sum_nocomm_ad.fh .
./gaussian_topog_r8.fh: $(SRCROOT)/FMS/topography/include/gaussian_topog_r8.fh
	cp $(SRCROOT)/FMS/topography/include/gaussian_topog_r8.fh .
./domain_read.inc: $(SRCROOT)/FMS/fms2_io/include/domain_read.inc
	cp $(SRCROOT)/FMS/fms2_io/include/domain_read.inc .
./mpp_io_write.inc: $(SRCROOT)/FMS/mpp/include/mpp_io_write.inc
	cp $(SRCROOT)/FMS/mpp/include/mpp_io_write.inc .
./mpp_update_domains2D.fh: $(SRCROOT)/FMS/mpp/include/mpp_update_domains2D.fh
	cp $(SRCROOT)/FMS/mpp/include/mpp_update_domains2D.fh .
./mpp_global_reduce.fh: $(SRCROOT)FMS/mpp/include/mpp_global_reduce.fh
	cp $(SRCROOT)FMS/mpp/include/mpp_global_reduce.fh .
./mpp_do_checkV.fh: $(SRCROOT)/FMS/mpp/include/mpp_do_checkV.fh
	cp $(SRCROOT)/FMS/mpp/include/mpp_do_checkV.fh .
./mpp_do_global_field_ad.fh: $(SRCROOT)/FMS/mpp/include/mpp_do_global_field_ad.fh
	cp $(SRCROOT)/FMS/mpp/include/mpp_do_global_field_ad.fh .
./horiz_interp_conserve_r8.fh: $(SRCROOT)/FMS/horiz_interp/include/horiz_interp_conserve_r8.fh
	cp $(SRCROOT)/FMS/horiz_interp/include/horiz_interp_conserve_r8.fh .
./mpp_type_mpi.fh: $(SRCROOT)FMS/mpp/include/mpp_type_mpi.fh
	cp $(SRCROOT)FMS/mpp/include/mpp_type_mpi.fh .
./mpp_util_mpi.inc: $(SRCROOT)/FMS/mpp/include/mpp_util_mpi.inc
	cp $(SRCROOT)/FMS/mpp/include/mpp_util_mpi.inc .
./netcdf_write_data.inc: $(SRCROOT)/FMS/fms2_io/include/netcdf_write_data.inc
	cp $(SRCROOT)/FMS/fms2_io/include/netcdf_write_data.inc .
./random_numbers_r8.fh: $(SRCROOT)/FMS/random_numbers/include/random_numbers_r8.fh
	cp $(SRCROOT)/FMS/random_numbers/include/random_numbers_r8.fh .
./mpp_util_nocomm.inc: $(SRCROOT)/FMS/mpp/include/mpp_util_nocomm.inc
	cp $(SRCROOT)/FMS/mpp/include/mpp_util_nocomm.inc .
./grid2_r8.fh: $(SRCROOT)/FMS/mosaic2/include/grid2_r8.fh
	cp $(SRCROOT)/FMS/mosaic2/include/grid2_r8.fh .
./time_manager.F90: $(SRCROOT)/FMS/time_manager/time_manager.F90
	cp $(SRCROOT)/FMS/time_manager/time_manager.F90 .
./diag_table.F90: $(SRCROOT)/FMS/diag_manager/diag_table.F90
	cp $(SRCROOT)/FMS/diag_manager/diag_table.F90 .
./mpp_sum.inc: $(SRCROOT)/FMS/mpp/include/mpp_sum.inc
	cp $(SRCROOT)/FMS/mpp/include/mpp_sum.inc .
./tracer_manager.inc: $(SRCROOT)/FMS/tracer_manager/include/tracer_manager.inc
	cp $(SRCROOT)/FMS/tracer_manager/include/tracer_manager.inc .
./axis_utils2_r8.fh: $(SRCROOT)/FMS/axis_utils/include/axis_utils2_r8.fh
	cp $(SRCROOT)/FMS/axis_utils/include/axis_utils2_r8.fh .
./mpp_error_s_a.fh: $(SRCROOT)FMS/mpp/include/mpp_error_s_a.fh
	cp $(SRCROOT)FMS/mpp/include/mpp_error_s_a.fh .
./monin_obukhov.inc: $(SRCROOT)/FMS/monin_obukhov/include/monin_obukhov.inc
	cp $(SRCROOT)/FMS/monin_obukhov/include/monin_obukhov.inc .
./gaussian_topog.inc: $(SRCROOT)/FMS/topography/include/gaussian_topog.inc
	cp $(SRCROOT)/FMS/topography/include/gaussian_topog.inc .
./field_manager.F90: $(SRCROOT)/FMS/field_manager/field_manager.F90
	cp $(SRCROOT)/FMS/field_manager/field_manager.F90 .
./mpp_io_misc.inc: $(SRCROOT)/FMS/mpp/include/mpp_io_misc.inc
	cp $(SRCROOT)/FMS/mpp/include/mpp_io_misc.inc .
./topography.F90: $(SRCROOT)/FMS/topography/topography.F90
	cp $(SRCROOT)/FMS/topography/topography.F90 .
./block_control.F90: $(SRCROOT)/FMS/block_control/block_control.F90
	cp $(SRCROOT)/FMS/block_control/block_control.F90 .
./fms_diag_fieldbuff_update.inc: $(SRCROOT)/FMS/diag_manager/include/fms_diag_fieldbuff_update.inc
	cp $(SRCROOT)/FMS/diag_manager/include/fms_diag_fieldbuff_update.inc .
./fm_util_r8.fh: $(SRCROOT)/FMS/field_manager/include/fm_util_r8.fh
	cp $(SRCROOT)/FMS/field_manager/include/fm_util_r8.fh .
./coupler_types_r4.fh: $(SRCROOT)/FMS/coupler/include/coupler_types_r4.fh
	cp $(SRCROOT)/FMS/coupler/include/coupler_types_r4.fh .
./mpp_do_updateV.fh: $(SRCROOT)FMS/mpp/include/mpp_do_updateV.fh
	cp $(SRCROOT)FMS/mpp/include/mpp_do_updateV.fh .
./tracer_manager.F90: $(SRCROOT)/FMS/tracer_manager/tracer_manager.F90
	cp $(SRCROOT)/FMS/tracer_manager/tracer_manager.F90 .
./get_variable_attribute.inc: $(SRCROOT)/FMS/fms2_io/include/get_variable_attribute.inc
	cp $(SRCROOT)/FMS/fms2_io/include/get_variable_attribute.inc .
./tridiagonal.F90: $(SRCROOT)/FMS/tridiagonal/tridiagonal.F90
	cp $(SRCROOT)/FMS/tridiagonal/tridiagonal.F90 .
./mpp_sum_mpi_ad.fh: $(SRCROOT)FMS/mpp/include/mpp_sum_mpi_ad.fh
	cp $(SRCROOT)FMS/mpp/include/mpp_sum_mpi_ad.fh .
./register_global_attribute.inc: $(SRCROOT)/FMS/fms2_io/include/register_global_attribute.inc
	cp $(SRCROOT)/FMS/fms2_io/include/register_global_attribute.inc .
./interpolator.F90: $(SRCROOT)/FMS/interpolator/interpolator.F90
	cp $(SRCROOT)/FMS/interpolator/interpolator.F90 .
./geos_constantsR4.fh: $(SRCROOT)/FMS/constants4/geos_constantsR4.fh
	cp $(SRCROOT)/FMS/constants4/geos_constantsR4.fh .
./sat_vapor_pres_k_r8.fh: $(SRCROOT)/FMS/sat_vapor_pres/include/sat_vapor_pres_k_r8.fh
	cp $(SRCROOT)/FMS/sat_vapor_pres/include/sat_vapor_pres_k_r8.fh .
./file_version.h: $(SRCROOT)/FMS/include/file_version.h
	cp $(SRCROOT)/FMS/include/file_version.h .
./mpp_io_read.inc: $(SRCROOT)/FMS/mpp/include/mpp_io_read.inc
	cp $(SRCROOT)/FMS/mpp/include/mpp_io_read.inc .
./get_grid_version.F90: $(SRCROOT)/FMS/data_override/get_grid_version.F90
	cp $(SRCROOT)/FMS/data_override/get_grid_version.F90 .
./fms_switches.h: $(SRCROOT)/FMS/drifters/fms_switches.h
	cp $(SRCROOT)/FMS/drifters/fms_switches.h .
./interpolator.inc: $(SRCROOT)/FMS/interpolator/include/interpolator.inc
	cp $(SRCROOT)/FMS/interpolator/include/interpolator.inc .
./mpp_sum_mpi.fh: $(SRCROOT)/FMS/mpp/include/mpp_sum_mpi.fh
	cp $(SRCROOT)/FMS/mpp/include/mpp_sum_mpi.fh .
./drifters_io.F90: $(SRCROOT)/FMS/drifters/drifters_io.F90
	cp $(SRCROOT)/FMS/drifters/drifters_io.F90 .
./mosaic2_r8.fh: $(SRCROOT)/FMS/mosaic2/include/mosaic2_r8.fh
	cp $(SRCROOT)/FMS/mosaic2/include/mosaic2_r8.fh .
./mpp_unstruct_domain.inc: $(SRCROOT)/FMS/mpp/include/mpp_unstruct_domain.inc
	cp $(SRCROOT)/FMS/mpp/include/mpp_unstruct_domain.inc .
./register_unstructured_domain_restart_variable.inc: $(SRCROOT)/FMS/fms2_io/include/register_unstructured_domain_restart_variable.inc
	cp $(SRCROOT)/FMS/fms2_io/include/register_unstructured_domain_restart_variable.inc .
./create_xgrid.c: $(SRCROOT)/FMS/mosaic/create_xgrid.c
	cp $(SRCROOT)/FMS/mosaic/create_xgrid.c .
./horiz_interp_spherical.inc: $(SRCROOT)/FMS/horiz_interp/include/horiz_interp_spherical.inc
	cp $(SRCROOT)/FMS/horiz_interp/include/horiz_interp_spherical.inc .
./cloud_interpolator.F90: $(SRCROOT)/FMS/drifters/cloud_interpolator.F90
	cp $(SRCROOT)/FMS/drifters/cloud_interpolator.F90 .
./gradient.F90: $(SRCROOT)/FMS/mosaic/gradient.F90
	cp $(SRCROOT)/FMS/mosaic/gradient.F90 .
./monin_obukhov_r8.fh: $(SRCROOT)/FMS/monin_obukhov/include/monin_obukhov_r8.fh
	cp $(SRCROOT)/FMS/monin_obukhov/include/monin_obukhov_r8.fh .
./fmsconstants.F90: $(SRCROOT)/FMS/constants/fmsconstants.F90
	cp $(SRCROOT)/FMS/constants/fmsconstants.F90 .
./mpp_error_s_s.fh: $(SRCROOT)FMS/mpp/include/mpp_error_s_s.fh
	cp $(SRCROOT)FMS/mpp/include/mpp_error_s_s.fh .
./unstructured_domain_write.inc: $(SRCROOT)/FMS/fms2_io/include/unstructured_domain_write.inc
	cp $(SRCROOT)/FMS/fms2_io/include/unstructured_domain_write.inc .
./mpp_alltoall_mpi.fh: $(SRCROOT)FMS/mpp/include/mpp_alltoall_mpi.fh
	cp $(SRCROOT)FMS/mpp/include/mpp_alltoall_mpi.fh .
./mpp_reduce_nocomm.fh: $(SRCROOT)/FMS/mpp/include/mpp_reduce_nocomm.fh
	cp $(SRCROOT)/FMS/mpp/include/mpp_reduce_nocomm.fh .
./gaussian_topog_r4.fh: $(SRCROOT)/FMS/topography/include/gaussian_topog_r4.fh
	cp $(SRCROOT)/FMS/topography/include/gaussian_topog_r4.fh .
./yaml_parser_binding.c: $(SRCROOT)/FMS/parser/yaml_parser_binding.c
	cp $(SRCROOT)/FMS/parser/yaml_parser_binding.c .
./fms_string_utils_r8.fh: $(SRCROOT)/FMS/string_utils/include/fms_string_utils_r8.fh
	cp $(SRCROOT)/FMS/string_utils/include/fms_string_utils_r8.fh .
./get_data_type_string.inc: $(SRCROOT)/FMS/fms2_io/include/get_data_type_string.inc
	cp $(SRCROOT)/FMS/fms2_io/include/get_data_type_string.inc .
./interp.h: $(SRCROOT)/FMS/mosaic/interp.h
	cp $(SRCROOT)/FMS/mosaic/interp.h .
./netcdf_io.F90: $(SRCROOT)/FMS/fms2_io/netcdf_io.F90
	cp $(SRCROOT)/FMS/fms2_io/netcdf_io.F90 .
./mpp_do_update.fh: $(SRCROOT)/FMS/mpp/include/mpp_do_update.fh
	cp $(SRCROOT)/FMS/mpp/include/mpp_do_update.fh .
./read_data_4d.inc: $(SRCROOT)/FMS/fms/read_data_4d.inc
	cp $(SRCROOT)/FMS/fms/read_data_4d.inc .
./mpp_error_s_a.fh: $(SRCROOT)/FMS/mpp/include/mpp_error_s_a.fh
	cp $(SRCROOT)/FMS/mpp/include/mpp_error_s_a.fh .
./mpp_chksum_int.fh: $(SRCROOT)/FMS/mpp/include/mpp_chksum_int.fh
	cp $(SRCROOT)/FMS/mpp/include/mpp_chksum_int.fh .
./gfdl_constantsR4.fh: $(SRCROOT)/FMS/constants4/gfdl_constantsR4.fh
	cp $(SRCROOT)/FMS/constants4/gfdl_constantsR4.fh .
./monin_obukhov.F90: $(SRCROOT)/FMS/monin_obukhov/monin_obukhov.F90
	cp $(SRCROOT)/FMS/monin_obukhov/monin_obukhov.F90 .
./mpp_memutils.F90: $(SRCROOT)/FMS/mpp/mpp_memutils.F90
	cp $(SRCROOT)/FMS/mpp/mpp_memutils.F90 .
./platform.F90: $(SRCROOT)/FMS/platform/platform.F90
	cp $(SRCROOT)/FMS/platform/platform.F90 .
./horiz_interp_bilinear_r4.fh: $(SRCROOT)/FMS/horiz_interp/include/horiz_interp_bilinear_r4.fh
	cp $(SRCROOT)/FMS/horiz_interp/include/horiz_interp_bilinear_r4.fh .
./horiz_interp_bicubic_r8.fh: $(SRCROOT)/FMS/horiz_interp/include/horiz_interp_bicubic_r8.fh
	cp $(SRCROOT)/FMS/horiz_interp/include/horiz_interp_bicubic_r8.fh .
./mpp_group_update.fh: $(SRCROOT)/FMS/mpp/include/mpp_group_update.fh
	cp $(SRCROOT)/FMS/mpp/include/mpp_group_update.fh .
./horiz_interp_type.inc: $(SRCROOT)/FMS/horiz_interp/include/horiz_interp_type.inc
	cp $(SRCROOT)/FMS/horiz_interp/include/horiz_interp_type.inc .
./mpp_error_a_s.fh: $(SRCROOT)/FMS/mpp/include/mpp_error_a_s.fh
	cp $(SRCROOT)/FMS/mpp/include/mpp_error_a_s.fh .
./horiz_interp.F90: $(SRCROOT)/FMS/horiz_interp/horiz_interp.F90
	cp $(SRCROOT)/FMS/horiz_interp/horiz_interp.F90 .
./group_update_unpack.inc: $(SRCROOT)/FMS/mpp/include/group_update_unpack.inc
	cp $(SRCROOT)/FMS/mpp/include/group_update_unpack.inc .
./atmos_ocean_fluxes.F90: $(SRCROOT)/FMS/coupler/atmos_ocean_fluxes.F90
	cp $(SRCROOT)/FMS/coupler/atmos_ocean_fluxes.F90 .
./mpp_gather.fh: $(SRCROOT)/FMS/mpp/include/mpp_gather.fh
	cp $(SRCROOT)/FMS/mpp/include/mpp_gather.fh .
./get_grid_version.inc: $(SRCROOT)/FMS/data_override/include/get_grid_version.inc
	cp $(SRCROOT)/FMS/data_override/include/get_grid_version.inc .
./fms_io_unstructured_save_restart.inc: $(SRCROOT)/FMS/fms/fms_io_unstructured_save_restart.inc
	cp $(SRCROOT)/FMS/fms/fms_io_unstructured_save_restart.inc .
./column_diagnostics.F90: $(SRCROOT)/FMS/column_diagnostics/column_diagnostics.F90
	cp $(SRCROOT)/FMS/column_diagnostics/column_diagnostics.F90 .
./mpp_error_a_s.fh: $(SRCROOT)FMS/mpp/include/mpp_error_a_s.fh
	cp $(SRCROOT)FMS/mpp/include/mpp_error_a_s.fh .
./create_xgrid.h: $(SRCROOT)/FMS/mosaic/create_xgrid.h
	cp $(SRCROOT)/FMS/mosaic/create_xgrid.h .
./diag_grid.F90: $(SRCROOT)/FMS/diag_manager/diag_grid.F90
	cp $(SRCROOT)/FMS/diag_manager/diag_grid.F90 .
./fms_io.F90: $(SRCROOT)/FMS/fms/fms_io.F90
	cp $(SRCROOT)/FMS/fms/fms_io.F90 .
./tridiagonal_r8.fh: $(SRCROOT)/FMS/tridiagonal/include/tridiagonal_r8.fh
	cp $(SRCROOT)/FMS/tridiagonal/include/tridiagonal_r8.fh .
./horiz_interp_r4.fh: $(SRCROOT)/FMS/horiz_interp/include/horiz_interp_r4.fh
	cp $(SRCROOT)/FMS/horiz_interp/include/horiz_interp_r4.fh .
./mpp_write_unlimited_axis.fh: $(SRCROOT)/FMS/mpp/include/mpp_write_unlimited_axis.fh
	cp $(SRCROOT)/FMS/mpp/include/mpp_write_unlimited_axis.fh .
./mpp_efp.F90: $(SRCROOT)/FMS/mpp/mpp_efp.F90
	cp $(SRCROOT)/FMS/mpp/mpp_efp.F90 .
./mpp_read_2Ddecomp.fh: $(SRCROOT)FMS/mpp/include/mpp_read_2Ddecomp.fh
	cp $(SRCROOT)FMS/mpp/include/mpp_read_2Ddecomp.fh .
./mpp_do_check.fh: $(SRCROOT)FMS/mpp/include/mpp_do_check.fh
	cp $(SRCROOT)FMS/mpp/include/mpp_do_check.fh .
./yaml_parser.F90: $(SRCROOT)/FMS/parser/yaml_parser.F90
	cp $(SRCROOT)/FMS/parser/yaml_parser.F90 .
./mpp_update_nest_domains.fh: $(SRCROOT)FMS/mpp/include/mpp_update_nest_domains.fh
	cp $(SRCROOT)FMS/mpp/include/mpp_update_nest_domains.fh .
./data_override.F90: $(SRCROOT)/FMS/data_override/data_override.F90
	cp $(SRCROOT)/FMS/data_override/data_override.F90 .
./mpp_alltoall_nocomm.fh: $(SRCROOT)FMS/mpp/include/mpp_alltoall_nocomm.fh
	cp $(SRCROOT)FMS/mpp/include/mpp_alltoall_nocomm.fh .
./sat_vapor_pres_k.F90: $(SRCROOT)/FMS/sat_vapor_pres/sat_vapor_pres_k.F90
	cp $(SRCROOT)/FMS/sat_vapor_pres/sat_vapor_pres_k.F90 .
./axis_utils2_r4.fh: $(SRCROOT)/FMS/axis_utils/include/axis_utils2_r4.fh
	cp $(SRCROOT)/FMS/axis_utils/include/axis_utils2_r4.fh .
./grid2.F90: $(SRCROOT)/FMS/mosaic2/grid2.F90
	cp $(SRCROOT)/FMS/mosaic2/grid2.F90 .
./mpp_utilities.F90: $(SRCROOT)/FMS/mpp/mpp_utilities.F90
	cp $(SRCROOT)/FMS/mpp/mpp_utilities.F90 .
./mpp_data.F90: $(SRCROOT)/FMS/mpp/mpp_data.F90
	cp $(SRCROOT)/FMS/mpp/mpp_data.F90 .
./mpp_data_mpi.inc: $(SRCROOT)/FMS/mpp/include/mpp_data_mpi.inc
	cp $(SRCROOT)/FMS/mpp/include/mpp_data_mpi.inc .
./mpp_get_boundary.fh: $(SRCROOT)/FMS/mpp/include/mpp_get_boundary.fh
	cp $(SRCROOT)/FMS/mpp/include/mpp_get_boundary.fh .
./mpp_get_boundary_ad.fh: $(SRCROOT)FMS/mpp/include/mpp_get_boundary_ad.fh
	cp $(SRCROOT)FMS/mpp/include/mpp_get_boundary_ad.fh .
./fms_stacksize.c: $(SRCROOT)/FMS/fms/fms_stacksize.c
	cp $(SRCROOT)/FMS/fms/fms_stacksize.c .
./fm_yaml.F90: $(SRCROOT)/FMS/field_manager/fm_yaml.F90
	cp $(SRCROOT)/FMS/field_manager/fm_yaml.F90 .
./gaussian_topog.F90: $(SRCROOT)/FMS/topography/gaussian_topog.F90
	cp $(SRCROOT)/FMS/topography/gaussian_topog.F90 .
./sat_vapor_pres_k_r4.fh: $(SRCROOT)/FMS/sat_vapor_pres/include/sat_vapor_pres_k_r4.fh
	cp $(SRCROOT)/FMS/sat_vapor_pres/include/sat_vapor_pres_k_r4.fh .
./mpp_do_update_nest.fh: $(SRCROOT)FMS/mpp/include/mpp_do_update_nest.fh
	cp $(SRCROOT)FMS/mpp/include/mpp_do_update_nest.fh .
./tracer_manager_r8.fh: $(SRCROOT)/FMS/tracer_manager/include/tracer_manager_r8.fh
	cp $(SRCROOT)/FMS/tracer_manager/include/tracer_manager_r8.fh .
./mpp_global_field_ug.fh: $(SRCROOT)/FMS/mpp/include/mpp_global_field_ug.fh
	cp $(SRCROOT)/FMS/mpp/include/mpp_global_field_ug.fh .
./coupler_types.inc: $(SRCROOT)/FMS/coupler/include/coupler_types.inc
	cp $(SRCROOT)/FMS/coupler/include/coupler_types.inc .
./mpp_comm_mpi.inc: $(SRCROOT)/FMS/mpp/include/mpp_comm_mpi.inc
	cp $(SRCROOT)/FMS/mpp/include/mpp_comm_mpi.inc .
./interpolator_r4.fh: $(SRCROOT)/FMS/interpolator/include/interpolator_r4.fh
	cp $(SRCROOT)/FMS/interpolator/include/interpolator_r4.fh .
./read_data_3d.inc: $(SRCROOT)/FMS/fms/read_data_3d.inc
	cp $(SRCROOT)/FMS/fms/read_data_3d.inc .
./compressed_read.inc: $(SRCROOT)/FMS/fms2_io/include/compressed_read.inc
	cp $(SRCROOT)/FMS/fms2_io/include/compressed_read.inc .
./mpp_do_update_nonblock.fh: $(SRCROOT)FMS/mpp/include/mpp_do_update_nonblock.fh
	cp $(SRCROOT)FMS/mpp/include/mpp_do_update_nonblock.fh .
./mpp_sum_nocomm.fh: $(SRCROOT)/FMS/mpp/include/mpp_sum_nocomm.fh
	cp $(SRCROOT)/FMS/mpp/include/mpp_sum_nocomm.fh .
./monin_obukhov_inter.inc: $(SRCROOT)/FMS/monin_obukhov/include/monin_obukhov_inter.inc
	cp $(SRCROOT)/FMS/monin_obukhov/include/monin_obukhov_inter.inc .
./horiz_interp_spherical.F90: $(SRCROOT)/FMS/horiz_interp/horiz_interp_spherical.F90
	cp $(SRCROOT)/FMS/horiz_interp/horiz_interp_spherical.F90 .
./mpp_transmit_nocomm.fh: $(SRCROOT)/FMS/mpp/include/mpp_transmit_nocomm.fh
	cp $(SRCROOT)/FMS/mpp/include/mpp_transmit_nocomm.fh .
./mpp_update_domains2D_nonblock.fh: $(SRCROOT)/FMS/mpp/include/mpp_update_domains2D_nonblock.fh
	cp $(SRCROOT)/FMS/mpp/include/mpp_update_domains2D_nonblock.fh .
./drifters_core.F90: $(SRCROOT)/FMS/drifters/drifters_core.F90
	cp $(SRCROOT)/FMS/drifters/drifters_core.F90 .
./gradient_c2l.c: $(SRCROOT)/FMS/mosaic/gradient_c2l.c
	cp $(SRCROOT)/FMS/mosaic/gradient_c2l.c .
./mpp_type_nocomm.fh: $(SRCROOT)FMS/mpp/include/mpp_type_nocomm.fh
	cp $(SRCROOT)FMS/mpp/include/mpp_type_nocomm.fh .
./mpp_parameter.F90: $(SRCROOT)/FMS/mpp/mpp_parameter.F90
	cp $(SRCROOT)/FMS/mpp/mpp_parameter.F90 .
./fms_r8.fh: $(SRCROOT)/FMS/fms/include/fms_r8.fh
	cp $(SRCROOT)/FMS/fms/include/fms_r8.fh .
./mpp_domains.F90: $(SRCROOT)/FMS/mpp/mpp_domains.F90
	cp $(SRCROOT)/FMS/mpp/mpp_domains.F90 .
./fms_io_unstructured_setup_one_field.inc: $(SRCROOT)/FMS/fms/fms_io_unstructured_setup_one_field.inc
	cp $(SRCROOT)/FMS/fms/fms_io_unstructured_setup_one_field.inc .
./diag_integral_r4.fh: $(SRCROOT)/FMS/diag_integral/include/diag_integral_r4.fh
	cp $(SRCROOT)/FMS/diag_integral/include/diag_integral_r4.fh .
./mpp_do_get_boundary.fh: $(SRCROOT)FMS/mpp/include/mpp_do_get_boundary.fh
	cp $(SRCROOT)FMS/mpp/include/mpp_do_get_boundary.fh .
./sat_vapor_pres.inc: $(SRCROOT)/FMS/sat_vapor_pres/include/sat_vapor_pres.inc
	cp $(SRCROOT)/FMS/sat_vapor_pres/include/sat_vapor_pres.inc .
./mpp_reduce_mpi.fh: $(SRCROOT)/FMS/mpp/include/mpp_reduce_mpi.fh
	cp $(SRCROOT)/FMS/mpp/include/mpp_reduce_mpi.fh .
./random_numbers.F90: $(SRCROOT)/FMS/random_numbers/random_numbers.F90
	cp $(SRCROOT)/FMS/random_numbers/random_numbers.F90 .
./time_interp_r4.fh: $(SRCROOT)/FMS/time_interp/include/time_interp_r4.fh
	cp $(SRCROOT)/FMS/time_interp/include/time_interp_r4.fh .
./horiz_interp_bilinear_r8.fh: $(SRCROOT)/FMS/horiz_interp/include/horiz_interp_bilinear_r8.fh
	cp $(SRCROOT)/FMS/horiz_interp/include/horiz_interp_bilinear_r8.fh .
./time_interp_external.F90: $(SRCROOT)/FMS/time_interp/time_interp_external.F90
	cp $(SRCROOT)/FMS/time_interp/time_interp_external.F90 .
./column_diagnostics.inc: $(SRCROOT)/FMS/column_diagnostics/include/column_diagnostics.inc
	cp $(SRCROOT)/FMS/column_diagnostics/include/column_diagnostics.inc .
./drifters_set_field.fh: $(SRCROOT)/FMS/drifters/drifters_set_field.fh
	cp $(SRCROOT)/FMS/drifters/drifters_set_field.fh .
./mpp_write_2Ddecomp.fh: $(SRCROOT)/FMS/mpp/include/mpp_write_2Ddecomp.fh
	cp $(SRCROOT)/FMS/mpp/include/mpp_write_2Ddecomp.fh .
./interpolator_r8.fh: $(SRCROOT)/FMS/interpolator/include/interpolator_r8.fh
	cp $(SRCROOT)/FMS/interpolator/include/interpolator_r8.fh .
./unpack_data.inc: $(SRCROOT)/FMS/fms2_io/include/unpack_data.inc
	cp $(SRCROOT)/FMS/fms2_io/include/unpack_data.inc .
./fms_string_utils_binding.c: $(SRCROOT)/FMS/string_utils/fms_string_utils_binding.c
	cp $(SRCROOT)/FMS/string_utils/fms_string_utils_binding.c .
./fms_io_unstructured_field_exist.inc: $(SRCROOT)/FMS/fms/fms_io_unstructured_field_exist.inc
	cp $(SRCROOT)/FMS/fms/fms_io_unstructured_field_exist.inc .
./diag_axis.F90: $(SRCROOT)/FMS/diag_manager/diag_axis.F90
	cp $(SRCROOT)/FMS/diag_manager/diag_axis.F90 .
./xgrid.F90: $(SRCROOT)/FMS/exchange/xgrid.F90
	cp $(SRCROOT)/FMS/exchange/xgrid.F90 .
./mpp_chksum_int.fh: $(SRCROOT)FMS/mpp/include/mpp_chksum_int.fh
	cp $(SRCROOT)FMS/mpp/include/mpp_chksum_int.fh .
./horiz_interp_bicubic.inc: $(SRCROOT)/FMS/horiz_interp/include/horiz_interp_bicubic.inc
	cp $(SRCROOT)/FMS/horiz_interp/include/horiz_interp_bicubic.inc .
./mpp_io_connect.inc: $(SRCROOT)/FMS/mpp/include/mpp_io_connect.inc
	cp $(SRCROOT)/FMS/mpp/include/mpp_io_connect.inc .
./mpp_type_nocomm.fh: $(SRCROOT)/FMS/mpp/include/mpp_type_nocomm.fh
	cp $(SRCROOT)/FMS/mpp/include/mpp_type_nocomm.fh .
./sat_vapor_pres_r8.fh: $(SRCROOT)/FMS/sat_vapor_pres/include/sat_vapor_pres_r8.fh
	cp $(SRCROOT)/FMS/sat_vapor_pres/include/sat_vapor_pres_r8.fh .
./mpp_io.F90: $(SRCROOT)/FMS/mpp/mpp_io.F90
	cp $(SRCROOT)/FMS/mpp/mpp_io.F90 .
./mpp_domains_util.inc: $(SRCROOT)/FMS/mpp/include/mpp_domains_util.inc
	cp $(SRCROOT)/FMS/mpp/include/mpp_domains_util.inc .
./monin_obukhov_inter.F90: $(SRCROOT)/FMS/monin_obukhov/monin_obukhov_inter.F90
	cp $(SRCROOT)/FMS/monin_obukhov/monin_obukhov_inter.F90 .
./mpp_update_nest_domains.fh: $(SRCROOT)/FMS/mpp/include/mpp_update_nest_domains.fh
	cp $(SRCROOT)/FMS/mpp/include/mpp_update_nest_domains.fh .
./mpp_chksum.fh: $(SRCROOT)FMS/mpp/include/mpp_chksum.fh
	cp $(SRCROOT)FMS/mpp/include/mpp_chksum.fh .
./amip_interp_r4.fh: $(SRCROOT)/FMS/amip_interp/include/amip_interp_r4.fh
	cp $(SRCROOT)/FMS/amip_interp/include/amip_interp_r4.fh .
./netcdf_add_restart_variable.inc: $(SRCROOT)/FMS/fms2_io/include/netcdf_add_restart_variable.inc
	cp $(SRCROOT)/FMS/fms2_io/include/netcdf_add_restart_variable.inc .
./astronomy_r8.fh: $(SRCROOT)/FMS/astronomy/include/astronomy_r8.fh
	cp $(SRCROOT)/FMS/astronomy/include/astronomy_r8.fh .
./mpp_sum_ad.inc: $(SRCROOT)FMS/mpp/include/mpp_sum_ad.inc
	cp $(SRCROOT)FMS/mpp/include/mpp_sum_ad.inc .
./mpp_do_checkV.fh: $(SRCROOT)FMS/mpp/include/mpp_do_checkV.fh
	cp $(SRCROOT)FMS/mpp/include/mpp_do_checkV.fh .
./mosaic_util.h: $(SRCROOT)/FMS/mosaic/mosaic_util.h
	cp $(SRCROOT)/FMS/mosaic/mosaic_util.h .
./horiz_interp_bicubic.F90: $(SRCROOT)/FMS/horiz_interp/horiz_interp_bicubic.F90
	cp $(SRCROOT)/FMS/horiz_interp/horiz_interp_bicubic.F90 .
./mpp_do_update_ad.fh: $(SRCROOT)FMS/mpp/include/mpp_do_update_ad.fh
	cp $(SRCROOT)FMS/mpp/include/mpp_do_update_ad.fh .
./fm_util.F90: $(SRCROOT)/FMS/field_manager/fm_util.F90
	cp $(SRCROOT)/FMS/field_manager/fm_util.F90 .
./fms_io_utils.F90: $(SRCROOT)/FMS/fms2_io/fms_io_utils.F90
	cp $(SRCROOT)/FMS/fms2_io/fms_io_utils.F90 .
./diag_integral.F90: $(SRCROOT)/FMS/diag_integral/diag_integral.F90
	cp $(SRCROOT)/FMS/diag_integral/diag_integral.F90 .
./mpp_domains_misc.inc: $(SRCROOT)/FMS/mpp/include/mpp_domains_misc.inc
	cp $(SRCROOT)/FMS/mpp/include/mpp_domains_misc.inc .
./gfs_constantsR4.fh: $(SRCROOT)/FMS/constants4/gfs_constantsR4.fh
	cp $(SRCROOT)/FMS/constants4/gfs_constantsR4.fh .
./mpp_comm_nocomm.inc: $(SRCROOT)/FMS/mpp/include/mpp_comm_nocomm.inc
	cp $(SRCROOT)/FMS/mpp/include/mpp_comm_nocomm.inc .
./mpp_global_field_ad.fh: $(SRCROOT)FMS/mpp/include/mpp_global_field_ad.fh
	cp $(SRCROOT)FMS/mpp/include/mpp_global_field_ad.fh .
./group_update_unpack.inc: $(SRCROOT)FMS/mpp/include/group_update_unpack.inc
	cp $(SRCROOT)FMS/mpp/include/group_update_unpack.inc .
./diag_integral.inc: $(SRCROOT)/FMS/diag_integral/include/diag_integral.inc
	cp $(SRCROOT)/FMS/diag_integral/include/diag_integral.inc .
./compressed_write.inc: $(SRCROOT)/FMS/fms2_io/include/compressed_write.inc
	cp $(SRCROOT)/FMS/fms2_io/include/compressed_write.inc .
./mpp_alltoall_nocomm.fh: $(SRCROOT)/FMS/mpp/include/mpp_alltoall_nocomm.fh
	cp $(SRCROOT)/FMS/mpp/include/mpp_alltoall_nocomm.fh .
./libFMS.F90: $(SRCROOT)/FMS/libFMS.F90
	cp $(SRCROOT)/FMS/libFMS.F90 .
./mpp_chksum_scalar.fh: $(SRCROOT)/FMS/mpp/include/mpp_chksum_scalar.fh
	cp $(SRCROOT)/FMS/mpp/include/mpp_chksum_scalar.fh .
./mpp_do_global_field.fh: $(SRCROOT)FMS/mpp/include/mpp_do_global_field.fh
	cp $(SRCROOT)FMS/mpp/include/mpp_do_global_field.fh .
./coupler_types_r8.fh: $(SRCROOT)/FMS/coupler/include/coupler_types_r8.fh
	cp $(SRCROOT)/FMS/coupler/include/coupler_types_r8.fh .
./mpp_read_2Ddecomp.fh: $(SRCROOT)/FMS/mpp/include/mpp_read_2Ddecomp.fh
	cp $(SRCROOT)/FMS/mpp/include/mpp_read_2Ddecomp.fh .
./mpp_write_compressed.fh: $(SRCROOT)FMS/mpp/include/mpp_write_compressed.fh
	cp $(SRCROOT)FMS/mpp/include/mpp_write_compressed.fh .
./monin_obukhov_inter_r4.fh: $(SRCROOT)/FMS/monin_obukhov/include/monin_obukhov_inter_r4.fh
	cp $(SRCROOT)/FMS/monin_obukhov/include/monin_obukhov_inter_r4.fh .
./fms.F90: $(SRCROOT)/FMS/fms/fms.F90
	cp $(SRCROOT)/FMS/fms/fms.F90 .
./mpp_transmit.inc: $(SRCROOT)/FMS/mpp/include/mpp_transmit.inc
	cp $(SRCROOT)/FMS/mpp/include/mpp_transmit.inc .
./fms_io_unstructured_register_restart_axis.inc: $(SRCROOT)/FMS/fms/fms_io_unstructured_register_restart_axis.inc
	cp $(SRCROOT)/FMS/fms/fms_io_unstructured_register_restart_axis.inc .
./horiz_interp_bilinear.inc: $(SRCROOT)/FMS/horiz_interp/include/horiz_interp_bilinear.inc
	cp $(SRCROOT)/FMS/horiz_interp/include/horiz_interp_bilinear.inc .
./mpp_transmit_mpi.fh: $(SRCROOT)FMS/mpp/include/mpp_transmit_mpi.fh
	cp $(SRCROOT)FMS/mpp/include/mpp_transmit_mpi.fh .
./fms_netcdf_unstructured_domain_io.F90: $(SRCROOT)/FMS/fms2_io/fms_netcdf_unstructured_domain_io.F90
	cp $(SRCROOT)/FMS/fms2_io/fms_netcdf_unstructured_domain_io.F90 .
./mpp_io_util.inc: $(SRCROOT)/FMS/mpp/include/mpp_io_util.inc
	cp $(SRCROOT)/FMS/mpp/include/mpp_io_util.inc .
./horiz_interp.inc: $(SRCROOT)/FMS/horiz_interp/include/horiz_interp.inc
	cp $(SRCROOT)/FMS/horiz_interp/include/horiz_interp.inc .
./mpp_error_a_a.fh: $(SRCROOT)FMS/mpp/include/mpp_error_a_a.fh
	cp $(SRCROOT)FMS/mpp/include/mpp_error_a_a.fh .
./mpp_domains_comm.inc: $(SRCROOT)/FMS/mpp/include/mpp_domains_comm.inc
	cp $(SRCROOT)/FMS/mpp/include/mpp_domains_comm.inc .
./gather_data_bc.inc: $(SRCROOT)/FMS/fms2_io/include/gather_data_bc.inc
	cp $(SRCROOT)/FMS/fms2_io/include/gather_data_bc.inc .
./array_utils.inc: $(SRCROOT)/FMS/fms2_io/include/array_utils.inc
	cp $(SRCROOT)/FMS/fms2_io/include/array_utils.inc .
./mosaic2.inc: $(SRCROOT)/FMS/mosaic2/include/mosaic2.inc
	cp $(SRCROOT)/FMS/mosaic2/include/mosaic2.inc .
./mpp_do_check.fh: $(SRCROOT)/FMS/mpp/include/mpp_do_check.fh
	cp $(SRCROOT)/FMS/mpp/include/mpp_do_check.fh .
./gfs_constants.fh: $(SRCROOT)/FMS/constants/gfs_constants.fh
	cp $(SRCROOT)/FMS/constants/gfs_constants.fh .
./mpp_update_domains2D_ad.fh: $(SRCROOT)FMS/mpp/include/mpp_update_domains2D_ad.fh
	cp $(SRCROOT)FMS/mpp/include/mpp_update_domains2D_ad.fh .
./fm_util_r4.fh: $(SRCROOT)/FMS/field_manager/include/fm_util_r4.fh
	cp $(SRCROOT)/FMS/field_manager/include/fm_util_r4.fh .
./ensemble_manager.F90: $(SRCROOT)/FMS/coupler/ensemble_manager.F90
	cp $(SRCROOT)/FMS/coupler/ensemble_manager.F90 .
./axis_utils.F90: $(SRCROOT)/FMS/axis_utils/axis_utils.F90
	cp $(SRCROOT)/FMS/axis_utils/axis_utils.F90 .
./geos_constants.fh: $(SRCROOT)/FMS/constants/geos_constants.fh
	cp $(SRCROOT)/FMS/constants/geos_constants.fh .
./get_grid_version_r4.fh: $(SRCROOT)/FMS/data_override/include/get_grid_version_r4.fh
	cp $(SRCROOT)/FMS/data_override/include/get_grid_version_r4.fh .
./axis_utils2.F90: $(SRCROOT)/FMS/axis_utils/axis_utils2.F90
	cp $(SRCROOT)/FMS/axis_utils/axis_utils2.F90 .
./horiz_interp_conserve_r4.fh: $(SRCROOT)/FMS/horiz_interp/include/horiz_interp_conserve_r4.fh
	cp $(SRCROOT)/FMS/horiz_interp/include/horiz_interp_conserve_r4.fh .
./mpp_write.fh: $(SRCROOT)FMS/mpp/include/mpp_write.fh
	cp $(SRCROOT)FMS/mpp/include/mpp_write.fh .
./mpp_do_updateV_ad.fh: $(SRCROOT)/FMS/mpp/include/mpp_do_updateV_ad.fh
	cp $(SRCROOT)/FMS/mpp/include/mpp_do_updateV_ad.fh .
./fms.inc: $(SRCROOT)/FMS/fms/include/fms.inc
	cp $(SRCROOT)/FMS/fms/include/fms.inc .
./fms2_io.F90: $(SRCROOT)/FMS/fms2_io/fms2_io.F90
	cp $(SRCROOT)/FMS/fms2_io/fms2_io.F90 .
./horiz_interp_bilinear.F90: $(SRCROOT)/FMS/horiz_interp/horiz_interp_bilinear.F90
	cp $(SRCROOT)/FMS/horiz_interp/horiz_interp_bilinear.F90 .
./blackboxio.F90: $(SRCROOT)/FMS/fms2_io/blackboxio.F90
	cp $(SRCROOT)/FMS/fms2_io/blackboxio.F90 .
SRC = $(SRCROOT)/FMS/mosaic/mosaic.F90 $(SRCROOT)/FMS/diag_manager/diag_axis.F90 $(SRCROOT)/FMS/coupler/atmos_ocean_fluxes.F90 $(SRCROOT)/FMS/field_manager/fm_yaml.F90 $(SRCROOT)/FMS/horiz_interp/horiz_interp.F90 $(SRCROOT)/FMS/mosaic/read_mosaic.c $(SRCROOT)/FMS/fms2_io/fms2_io.F90 $(SRCROOT)/FMS/time_interp/time_interp_external.F90 $(SRCROOT)/FMS/mpp/mpp_pset.F90 $(SRCROOT)/FMS/diag_manager/fms_diag_time_reduction.F90 $(SRCROOT)/FMS/monin_obukhov/monin_obukhov.F90 $(SRCROOT)/FMS/time_interp/time_interp.F90 $(SRCROOT)/FMS/diag_manager/diag_manager.F90 $(SRCROOT)/FMS/tracer_manager/tracer_manager.F90 $(SRCROOT)/FMS/horiz_interp/horiz_interp_bicubic.F90 $(SRCROOT)/FMS/horiz_interp/horiz_interp_conserve.F90 $(SRCROOT)/FMS/astronomy/astronomy.F90 $(SRCROOT)/FMS/time_interp/time_interp_external2.F90 $(SRCROOT)/FMS/coupler/coupler_types.F90 $(SRCROOT)/FMS/fms/fms_io.F90 $(SRCROOT)/FMS/column_diagnostics/column_diagnostics.F90 $(SRCROOT)/FMS/mpp/mpp_parameter.F90 $(SRCROOT)/FMS/data_override/data_override.F90 $(SRCROOT)/FMS/affinity/affinity.c $(SRCROOT)/FMS/libFMS.F90 $(SRCROOT)/FMS/topography/topography.F90 $(SRCROOT)/FMS/drifters/drifters_core.F90 $(SRCROOT)/FMS/horiz_interp/horiz_interp_type.F90 $(SRCROOT)/FMS/sat_vapor_pres/sat_vapor_pres.F90 $(SRCROOT)/FMS/mosaic/gradient_c2l.c $(SRCROOT)/FMS/affinity/fms_affinity.F90 $(SRCROOT)/FMS/mosaic/grid.F90 $(SRCROOT)/FMS/diag_manager/diag_table.F90 $(SRCROOT)/FMS/field_manager/fm_util.F90 $(SRCROOT)/FMS/field_manager/field_manager.F90 $(SRCROOT)/FMS/constants/fmsconstants.F90 $(SRCROOT)/FMS/parser/yaml_parser_binding.c $(SRCROOT)/FMS/mosaic2/mosaic2.F90 $(SRCROOT)/FMS/axis_utils/axis_utils2.F90 $(SRCROOT)/FMS/interpolator/interpolator.F90 $(SRCROOT)/FMS/mpp/mpp_utilities.F90 $(SRCROOT)/FMS/parser/fms_yaml_output.F90 $(SRCROOT)/FMS/mpp/mpp.F90 $(SRCROOT)/FMS/mosaic/gradient.F90 $(SRCROOT)/FMS/drifters/drifters_input.F90 $(SRCROOT)/FMS/fms/fms_stacksize.c $(SRCROOT)/FMS/horiz_interp/horiz_interp_spherical.F90 $(SRCROOT)/FMS/exchange/xgrid.F90 $(SRCROOT)/FMS/fms2_io/fms_netcdf_unstructured_domain_io.F90 $(SRCROOT)/FMS/diag_manager/fms_diag_fieldbuff_update.F90 $(SRCROOT)/FMS/fms2_io/netcdf_io.F90 $(SRCROOT)/FMS/drifters/drifters_comm.F90 $(SRCROOT)/FMS/fms/fms.F90 $(SRCROOT)/FMS/constants4/fmsconstantsr4.F90 $(SRCROOT)/FMS/drifters/quicksort.F90 $(SRCROOT)/FMS/diag_manager/fms_diag_elem_weight_procs.F90 $(SRCROOT)/FMS/drifters/drifters_io.F90 $(SRCROOT)/FMS/diag_manager/diag_data.F90 $(SRCROOT)/FMS/random_numbers/mersennetwister.F90 $(SRCROOT)/FMS/parser/yaml_parser.F90 $(SRCROOT)/FMS/diag_manager/diag_output.F90 $(SRCROOT)/FMS/fms2_io/blackboxio.F90 $(SRCROOT)/FMS/fms2_io/fms_netcdf_domain_io.F90 $(SRCROOT)/FMS/platform/platform.F90 $(SRCROOT)/FMS/fms2_io/fms_io_utils.F90 $(SRCROOT)/FMS/diag_integral/diag_integral.F90 $(SRCROOT)/FMS/tridiagonal/tridiagonal.F90 $(SRCROOT)/FMS/parser/yaml_output_functions.c $(SRCROOT)/FMS/diag_manager/diag_grid.F90 $(SRCROOT)/FMS/axis_utils/axis_utils.F90 $(SRCROOT)/FMS/block_control/block_control.F90 $(SRCROOT)/FMS/sat_vapor_pres/sat_vapor_pres_k.F90 $(SRCROOT)/FMS/time_manager/time_manager.F90 $(SRCROOT)/FMS/data_override/get_grid_version.F90 $(SRCROOT)/FMS/constants4/constantsr4.F90 $(SRCROOT)/FMS/mpp/mpp_data.F90 $(SRCROOT)/FMS/diag_manager/fms_diag_bbox.F90 $(SRCROOT)/FMS/mpp/mpp_efp.F90 $(SRCROOT)/FMS/mpp/mpp_io.F90 $(SRCROOT)/FMS/diag_manager/fms_diag_outfield.F90 $(SRCROOT)/FMS/memutils/memutils.F90 $(SRCROOT)/FMS/drifters/drifters.F90 $(SRCROOT)/FMS/mpp/mpp_memutils.F90 $(SRCROOT)/FMS/mosaic/interp.c $(SRCROOT)/FMS/string_utils/fms_string_utils.F90 $(SRCROOT)/FMS/time_manager/get_cal_time.F90 $(SRCROOT)/FMS/monin_obukhov/monin_obukhov_inter.F90 $(SRCROOT)/FMS/mosaic/create_xgrid.c $(SRCROOT)/FMS/exchange/stock_constants.F90 $(SRCROOT)/FMS/horiz_interp/horiz_interp_bilinear.F90 $(SRCROOT)/FMS/mpp/mpp_memuse.c $(SRCROOT)/FMS/mosaic/mosaic_util.c $(SRCROOT)/FMS/random_numbers/random_numbers.F90 $(SRCROOT)/FMS/topography/gaussian_topog.F90 $(SRCROOT)/FMS/string_utils/fms_string_utils_binding.c $(SRCROOT)/FMS/amip_interp/amip_interp.F90 $(SRCROOT)/FMS/drifters/cloud_interpolator.F90 $(SRCROOT)/FMS/mpp/mpp_domains.F90 $(SRCROOT)/FMS/mosaic2/grid2.F90 $(SRCROOT)/FMS/coupler/ensemble_manager.F90 $(SRCROOT)/FMS/constants/constants.F90 $(SRCROOT)/FMS/diag_manager/diag_util.F90 $(SRCROOT)/FMS/time_interp/include/time_interp_r8.fh $(SRCROOT)/FMS/fms2_io/include/compressed_write.inc $(SRCROOT)/FMS/horiz_interp/include/horiz_interp_r8.fh $(SRCROOT)/FMS/axis_utils/include/axis_utils2.inc $(SRCROOT)/FMS/field_manager/include/field_manager_r4.fh $(SRCROOT)/FMS/diag_integral/include/diag_integral.inc $(SRCROOT)/FMS/mpp/include/mpp_chksum_scalar.fh $(SRCROOT)FMS/mpp/include/mpp_do_get_boundary_ad.fh $(SRCROOT)FMS/mpp/include/mpp_sum.inc $(SRCROOT)/FMS/tridiagonal/include/tridiagonal_r4.fh $(SRCROOT)/FMS/mpp/include/mpp_alltoall_nocomm.fh $(SRCROOT)/FMS/time_interp/include/time_interp_external2_r4.fh $(SRCROOT)/FMS/mpp/include/mpp_read_2Ddecomp.fh $(SRCROOT)/FMS/coupler/include/coupler_types_r8.fh $(SRCROOT)/FMS/fms/fms_io_unstructured_register_restart_field.inc $(SRCROOT)FMS/mpp/include/mpp_do_global_field.fh $(SRCROOT)/FMS/fms2_io/include/netcdf_read_data.inc $(SRCROOT)FMS/mpp/include/mpp_transmit.inc $(SRCROOT)/FMS/mpp/include/mpp_do_get_boundary_ad.fh $(SRCROOT)FMS/mpp/include/mpp_gather.fh $(SRCROOT)/FMS/mpp/include/mpp_sum_nocomm_ad.fh $(SRCROOT)/FMS/monin_obukhov/include/monin_obukhov_inter_r4.fh $(SRCROOT)FMS/mpp/include/mpp_global_sum_tl.fh $(SRCROOT)FMS/mpp/include/mpp_write_compressed.fh $(SRCROOT)FMS/mpp/include/mpp_global_field_ug.fh $(SRCROOT)/FMS/mosaic2/include/mosaic2_r4.fh $(SRCROOT)FMS/mpp/include/mpp_do_update_ad.fh $(SRCROOT)/FMS/mpp/include/mpp_error_a_a.fh $(SRCROOT)/FMS/mpp/include/mpp_read_compressed.fh $(SRCROOT)/FMS/tridiagonal/include/tridiagonal.inc $(SRCROOT)/FMS/string_utils/include/fms_string_utils_r4.fh $(SRCROOT)/FMS/mpp/include/mpp_domains_misc.inc $(SRCROOT)/FMS/mpp/include/mpp_transmit_mpi.fh $(SRCROOT)FMS/mpp/include/mpp_get_boundary.fh $(SRCROOT)/FMS/mpp/include/group_update_pack.inc $(SRCROOT)/FMS/horiz_interp/include/horiz_interp_spherical_r8.fh $(SRCROOT)FMS/mpp/include/group_update_unpack.inc $(SRCROOT)FMS/mpp/include/mpp_do_update.fh $(SRCROOT)FMS/mpp/include/mpp_global_field_ad.fh $(SRCROOT)/FMS/constants4/gfs_constantsR4.fh $(SRCROOT)/FMS/mpp/include/mpp_data_nocomm.inc $(SRCROOT)/FMS/mpp/include/mpp_comm_nocomm.inc $(SRCROOT)/FMS/data_override/include/get_grid_version_r4.fh $(SRCROOT)/FMS/tracer_manager/include/tracer_manager.inc $(SRCROOT)/FMS/constants/geos_constants.fh $(SRCROOT)/FMS/mpp/include/mpp_sum.inc $(SRCROOT)/FMS/field_manager/include/fm_util_r4.fh $(SRCROOT)/FMS/mosaic2/include/grid2_r8.fh $(SRCROOT)/FMS/monin_obukhov/include/monin_obukhov.inc $(SRCROOT)FMS/mpp/include/mpp_error_s_a.fh $(SRCROOT)FMS/mpp/include/mpp_write.fh $(SRCROOT)/FMS/axis_utils/include/axis_utils2_r8.fh $(SRCROOT)/FMS/horiz_interp/include/horiz_interp_conserve_r4.fh $(SRCROOT)/FMS/mpp/include/mpp_do_updateV_ad.fh $(SRCROOT)/FMS/topography/include/gaussian_topog.inc $(SRCROOT)/FMS/mpp/include/mpp_io_misc.inc $(SRCROOT)/FMS/fms/include/fms.inc $(SRCROOT)/FMS/topography/include/gaussian_topog_r8.fh $(SRCROOT)/FMS/horiz_interp/include/horiz_interp_bilinear.inc $(SRCROOT)FMS/mpp/include/mpp_sum_nocomm_ad.fh $(SRCROOT)/FMS/fms/fms_io_unstructured_register_restart_axis.inc $(SRCROOT)/FMS/horiz_interp/include/horiz_interp_type_r8.fh $(SRCROOT)/FMS/mpp/include/mpp_transmit.inc $(SRCROOT)/FMS/mpp/include/mpp_do_global_field_ad.fh $(SRCROOT)/FMS/mpp/include/mpp_io_util.inc $(SRCROOT)/FMS/mpp/include/mpp_do_checkV.fh $(SRCROOT)FMS/mpp/include/mpp_global_reduce.fh $(SRCROOT)/FMS/mpp/include/mpp_update_domains2D.fh $(SRCROOT)/FMS/fms2_io/include/domain_read.inc $(SRCROOT)/FMS/mpp/include/mpp_io_write.inc $(SRCROOT)FMS/mpp/include/mpp_transmit_mpi.fh $(SRCROOT)/FMS/fms2_io/include/gather_data_bc.inc $(SRCROOT)/FMS/mpp/include/mpp_domains_comm.inc $(SRCROOT)/FMS/fms2_io/include/netcdf_write_data.inc $(SRCROOT)/FMS/mpp/include/mpp_util_mpi.inc $(SRCROOT)FMS/mpp/include/mpp_error_a_a.fh $(SRCROOT)FMS/mpp/include/mpp_type_mpi.fh $(SRCROOT)/FMS/horiz_interp/include/horiz_interp.inc $(SRCROOT)/FMS/horiz_interp/include/horiz_interp_conserve_r8.fh $(SRCROOT)FMS/mpp/include/mpp_update_domains2D_ad.fh $(SRCROOT)/FMS/mpp/include/mpp_util_nocomm.inc $(SRCROOT)/FMS/mpp/include/mpp_do_check.fh $(SRCROOT)/FMS/constants/gfs_constants.fh $(SRCROOT)/FMS/mosaic2/include/mosaic2.inc $(SRCROOT)/FMS/fms2_io/include/array_utils.inc $(SRCROOT)/FMS/random_numbers/include/random_numbers_r8.fh $(SRCROOT)/FMS/topography/include/topography_r4.fh $(SRCROOT)/FMS/horiz_interp/include/horiz_interp_bilinear_r8.fh $(SRCROOT)/FMS/drifters/drifters_push.fh $(SRCROOT)/FMS/column_diagnostics/include/column_diagnostics.inc $(SRCROOT)/FMS/time_interp/include/time_interp_r4.fh $(SRCROOT)/FMS/mpp/include/mpp_reduce_mpi.fh $(SRCROOT)/FMS/sat_vapor_pres/include/sat_vapor_pres.inc $(SRCROOT)/FMS/fms2_io/include/unpack_data.inc $(SRCROOT)/FMS/interpolator/include/interpolator_r8.fh $(SRCROOT)/FMS/mpp/include/mpp_write_2Ddecomp.fh $(SRCROOT)/FMS/drifters/drifters_set_field.fh $(SRCROOT)/FMS/fms2_io/include/domain_write.inc $(SRCROOT)/FMS/random_numbers/include/random_numbers.inc $(SRCROOT)FMS/mpp/include/mpp_sum_mpi.fh $(SRCROOT)/FMS/data_override/include/data_override_r8.fh $(SRCROOT)/FMS/diag_manager/include/fms_diag_fieldbuff_update.fh $(SRCROOT)FMS/mpp/include/mpp_read_distributed_ascii.fh $(SRCROOT)/FMS/fms2_io/include/compute_global_checksum.inc $(SRCROOT)/FMS/mpp/include/mpp_do_update_ad.fh $(SRCROOT)/FMS/topography/include/topography_r8.fh $(SRCROOT)/FMS/monin_obukhov/include/monin_obukhov_inter.inc $(SRCROOT)/FMS/horiz_interp/include/horiz_interp_conserve.inc $(SRCROOT)/FMS/mpp/include/mpp_update_domains2D_nonblock.fh $(SRCROOT)/FMS/mpp/include/mpp_transmit_nocomm.fh $(SRCROOT)/FMS/fms/fms_io_unstructured_setup_one_field.inc $(SRCROOT)/FMS/mpp/include/mpp_do_get_boundary.fh $(SRCROOT)FMS/mpp/include/mpp_do_updateV_ad.fh $(SRCROOT)/FMS/fms/include/fms_r8.fh $(SRCROOT)FMS/mpp/include/mpp_type_nocomm.fh $(SRCROOT)/FMS/horiz_interp/include/horiz_interp_type_r4.fh $(SRCROOT)/FMS/string_utils/include/fms_string_utils.inc $(SRCROOT)/FMS/mpp/include/mpp_global_field_ad.fh $(SRCROOT)FMS/mpp/include/mpp_do_get_boundary.fh $(SRCROOT)/FMS/diag_integral/include/diag_integral_r4.fh $(SRCROOT)FMS/mpp/include/mpp_global_field.fh $(SRCROOT)/FMS/mpp/include/mpp_util.inc $(SRCROOT)/FMS/fms/fms_io_unstructured_read.inc $(SRCROOT)/FMS/amip_interp/include/amip_interp_r4.fh $(SRCROOT)FMS/mpp/include/mpp_chksum.fh $(SRCROOT)/FMS/mpp/include/mpp_alltoall_mpi.fh $(SRCROOT)/FMS/fms/write_data.inc $(SRCROOT)/FMS/topography/include/topography.inc $(SRCROOT)FMS/include/file_version.h $(SRCROOT)FMS/mpp/include/mpp_sum_nocomm.fh $(SRCROOT)/FMS/fms2_io/include/netcdf_add_restart_variable.inc $(SRCROOT)/FMS/mpp/include/mpp_define_nest_domains.inc $(SRCROOT)/FMS/include/fms_platform.h $(SRCROOT)FMS/mpp/include/group_update_pack.inc $(SRCROOT)/FMS/mpp/include/mpp_write.fh $(SRCROOT)/FMS/mosaic2/include/grid2.inc $(SRCROOT)FMS/mpp/include/mpp_do_checkV.fh $(SRCROOT)FMS/mpp/include/mpp_sum_ad.inc $(SRCROOT)/FMS/astronomy/include/astronomy_r8.fh $(SRCROOT)/FMS/monin_obukhov/include/monin_obukhov_r4.fh $(SRCROOT)/FMS/mpp/include/mpp_update_domains2D_ad.fh $(SRCROOT)/FMS/mosaic/mosaic_util.h $(SRCROOT)FMS/mpp/include/mpp_chksum_int.fh $(SRCROOT)FMS/mpp/include/mpp_global_sum.fh $(SRCROOT)/FMS/fms/fms_io_unstructured_field_exist.inc $(SRCROOT)/FMS/random_numbers/include/random_numbers_r4.fh $(SRCROOT)/FMS/fms/fms_io_unstructured_get_field_size.inc $(SRCROOT)/FMS/mpp/include/mpp_scatter.fh $(SRCROOT)/FMS/horiz_interp/include/horiz_interp_bicubic.inc $(SRCROOT)/FMS/mpp/include/mpp_do_updateV_nonblock.fh $(SRCROOT)FMS/mpp/include/mpp_do_redistribute.fh $(SRCROOT)/FMS/sat_vapor_pres/include/sat_vapor_pres_r8.fh $(SRCROOT)/FMS/mpp/include/mpp_type_nocomm.fh $(SRCROOT)/FMS/mpp/include/mpp_io_connect.inc $(SRCROOT)/FMS/constants/gfdl_constants.fh $(SRCROOT)/FMS/mpp/include/mpp_update_nest_domains.fh $(SRCROOT)/FMS/fms/include/fms_r4.fh $(SRCROOT)/FMS/fms/read_data_2d.inc $(SRCROOT)/FMS/mpp/include/mpp_domains_util.inc $(SRCROOT)/FMS/mpp/include/mpp_global_reduce.fh $(SRCROOT)/FMS/fms/fms_io_unstructured_save_restart.inc $(SRCROOT)/FMS/mpp/include/mpp_error_s_s.fh $(SRCROOT)/FMS/data_override/include/get_grid_version.inc $(SRCROOT)/FMS/astronomy/include/astronomy.inc $(SRCROOT)/FMS/mpp/include/mpp_gather.fh $(SRCROOT)/FMS/mosaic/read_mosaic.h $(SRCROOT)/FMS/mosaic/create_xgrid.h $(SRCROOT)FMS/mpp/include/mpp_error_a_s.fh $(SRCROOT)/FMS/column_diagnostics/include/column_diagnostics_r4.fh $(SRCROOT)/FMS/fms2_io/include/register_domain_restart_variable.inc $(SRCROOT)/FMS/tridiagonal/include/tridiagonal_r8.fh $(SRCROOT)FMS/mpp/include/mpp_read_2Ddecomp.fh $(SRCROOT)/FMS/mpp/include/mpp_global_sum_tl.fh $(SRCROOT)/FMS/fms2_io/include/register_variable_attribute.inc $(SRCROOT)/FMS/mpp/include/mpp_write_unlimited_axis.fh $(SRCROOT)/FMS/fms/fms_io_unstructured_file_unit.inc $(SRCROOT)/FMS/horiz_interp/include/horiz_interp_r4.fh $(SRCROOT)/FMS/mosaic2/include/grid2_r4.fh $(SRCROOT)FMS/mpp/include/mpp_write_2Ddecomp.fh $(SRCROOT)/FMS/mpp/include/mpp_write_compressed.fh $(SRCROOT)/FMS/field_manager/include/fm_util.inc $(SRCROOT)/FMS/constants4/gfdl_constantsR4.fh $(SRCROOT)/FMS/mpp/include/mpp_unstruct_pass_data.fh $(SRCROOT)/FMS/mosaic/gradient_c2l.h $(SRCROOT)/FMS/mpp/include/mpp_chksum.fh $(SRCROOT)/FMS/mpp/include/mpp_io_unstructured_write.inc $(SRCROOT)/FMS/horiz_interp/include/horiz_interp_bilinear_r4.fh $(SRCROOT)FMS/mpp/include/mpp_reduce_nocomm.fh $(SRCROOT)FMS/mpp/include/mpp_reduce_mpi.fh $(SRCROOT)/FMS/fms2_io/include/scatter_data_bc.inc $(SRCROOT)/FMS/mpp/include/mpp_error_a_s.fh $(SRCROOT)/FMS/mpp/include/mpp_group_update.fh $(SRCROOT)/FMS/horiz_interp/include/horiz_interp_type.inc $(SRCROOT)/FMS/horiz_interp/include/horiz_interp_bicubic_r8.fh $(SRCROOT)/FMS/mosaic/constant.h $(SRCROOT)/FMS/mpp/include/mpp_global_sum_ad.fh $(SRCROOT)/FMS/mpp/include/group_update_unpack.inc $(SRCROOT)/FMS/field_manager/include/field_manager.inc $(SRCROOT)/FMS/astronomy/include/astronomy_r4.fh $(SRCROOT)FMS/mpp/include/mpp_do_updateV_nonblock.fh $(SRCROOT)/FMS/tracer_manager/include/tracer_manager_r8.fh $(SRCROOT)/FMS/mpp/include/mpp_global_field_ug.fh $(SRCROOT)FMS/mpp/include/mpp_do_update_nest.fh $(SRCROOT)/FMS/mpp/include/mpp_do_updateV.fh $(SRCROOT)/FMS/sat_vapor_pres/include/sat_vapor_pres_k_r4.fh $(SRCROOT)/FMS/mpp/include/mpp_comm_mpi.inc $(SRCROOT)/FMS/coupler/include/coupler_types.inc $(SRCROOT)/FMS/fms2_io/include/compressed_read.inc $(SRCROOT)/FMS/fms/read_data_3d.inc $(SRCROOT)/FMS/field_manager/include/field_manager_r8.fh $(SRCROOT)/FMS/interpolator/include/interpolator_r4.fh $(SRCROOT)/FMS/mpp/include/mpp_type_mpi.fh $(SRCROOT)/FMS/mpp/include/mpp_sum_nocomm.fh $(SRCROOT)/FMS/amip_interp/include/amip_interp_r8.fh $(SRCROOT)FMS/mpp/include/mpp_do_update_nonblock.fh $(SRCROOT)/FMS/mpp/include/mpp_do_redistribute.fh $(SRCROOT)/FMS/mpp/include/mpp_io_unstructured_read.inc $(SRCROOT)FMS/mpp/include/mpp_group_update.fh $(SRCROOT)/FMS/fms2_io/include/array_utils_char.inc $(SRCROOT)FMS/mpp/include/mpp_do_check.fh $(SRCROOT)/FMS/axis_utils/include/axis_utils2_r4.fh $(SRCROOT)FMS/mpp/include/mpp_unstruct_pass_data.fh $(SRCROOT)FMS/mpp/include/mpp_alltoall_nocomm.fh $(SRCROOT)/FMS/mpp/include/mpp_global_field.fh $(SRCROOT)/FMS/tracer_manager/include/tracer_manager_r4.fh $(SRCROOT)FMS/mpp/include/mpp_update_nest_domains.fh $(SRCROOT)/FMS/fms2_io/include/unstructured_domain_read.inc $(SRCROOT)/FMS/mpp/include/mpp_data_mpi.inc $(SRCROOT)/FMS/time_interp/include/time_interp_external2_r8.fh $(SRCROOT)FMS/mpp/include/mpp_global_sum_ad.fh $(SRCROOT)FMS/mpp/include/mpp_get_boundary_ad.fh $(SRCROOT)/FMS/mpp/include/mpp_get_boundary.fh $(SRCROOT)/FMS/mpp/include/mpp_sum_ad.inc $(SRCROOT)/FMS/mpp/include/mpp_read_distributed_ascii.inc $(SRCROOT)/FMS/sat_vapor_pres/include/sat_vapor_pres_r4.fh $(SRCROOT)/FMS/mpp/include/mpp_sum_mpi_ad.fh $(SRCROOT)/FMS/fms2_io/include/register_global_attribute.inc $(SRCROOT)/FMS/diag_integral/include/diag_integral_r8.fh $(SRCROOT)/FMS/mpp/include/mpp_domains_reduce.inc $(SRCROOT)/FMS/drifters/fms_switches.h $(SRCROOT)/FMS/include/file_version.h $(SRCROOT)/FMS/time_interp/include/time_interp.inc $(SRCROOT)/FMS/mpp/include/mpp_io_read.inc $(SRCROOT)/FMS/sat_vapor_pres/include/sat_vapor_pres_k_r8.fh $(SRCROOT)/FMS/constants4/geos_constantsR4.fh $(SRCROOT)/FMS/mpp/include/system_clock.fh $(SRCROOT)/FMS/mpp/include/mpp_unstruct_domain.inc $(SRCROOT)/FMS/mosaic2/include/mosaic2_r8.fh $(SRCROOT)/FMS/interpolator/include/interpolator.inc $(SRCROOT)/FMS/horiz_interp/include/horiz_interp_bicubic_r4.fh $(SRCROOT)/FMS/mpp/include/mpp_sum_mpi.fh $(SRCROOT)FMS/mpp/include/mpp_chksum_scalar.fh $(SRCROOT)/FMS/field_manager/parse.inc $(SRCROOT)/FMS/data_override/include/get_grid_version_r8.fh $(SRCROOT)/FMS/horiz_interp/include/horiz_interp_spherical_r4.fh $(SRCROOT)/FMS/mpp/include/mpp_global_sum.fh $(SRCROOT)/FMS/monin_obukhov/include/monin_obukhov_inter_r8.fh $(SRCROOT)/FMS/amip_interp/include/amip_interp.inc $(SRCROOT)/FMS/drifters/drifters_compute_k.fh $(SRCROOT)/FMS/mpp/include/mpp_domains_define.inc $(SRCROOT)FMS/mpp/include/mpp_read_compressed.fh $(SRCROOT)/FMS/field_manager/include/fm_util_r8.fh $(SRCROOT)/FMS/coupler/include/coupler_types_r4.fh $(SRCROOT)/FMS/mpp/include/mpp_do_global_field.fh $(SRCROOT)/FMS/diag_manager/include/fms_diag_fieldbuff_update.inc $(SRCROOT)FMS/mpp/include/mpp_do_updateV.fh $(SRCROOT)FMS/mpp/include/mpp_sum_mpi_ad.fh $(SRCROOT)/FMS/fms2_io/include/get_variable_attribute.inc $(SRCROOT)/FMS/mpp/include/mpp_read_distributed_ascii.fh $(SRCROOT)/FMS/mpp/include/mpp_do_update_nest.fh $(SRCROOT)/FMS/time_interp/include/time_interp_external2.inc $(SRCROOT)/FMS/fms2_io/include/get_global_attribute.inc $(SRCROOT)/FMS/mpp/include/mpp_reduce_nocomm.fh $(SRCROOT)/FMS/mpp/include/mpp_do_update_nonblock.fh $(SRCROOT)/FMS/fms/fms_io_unstructured_get_file_unit.inc $(SRCROOT)/FMS/data_override/include/data_override.inc $(SRCROOT)/FMS/topography/include/gaussian_topog_r4.fh $(SRCROOT)FMS/mpp/include/mpp_update_domains2D_nonblock.fh $(SRCROOT)/FMS/mosaic/interp.h $(SRCROOT)FMS/mpp/include/mpp_transmit_nocomm.fh $(SRCROOT)/FMS/fms2_io/include/get_data_type_string.inc $(SRCROOT)/FMS/string_utils/include/fms_string_utils_r8.fh $(SRCROOT)/FMS/fms/fms_io_unstructured_get_file_name.inc $(SRCROOT)FMS/mpp/include/mpp_do_global_field_ad.fh $(SRCROOT)/FMS/mpp/include/mpp_chksum_int.fh $(SRCROOT)/FMS/mpp/include/mpp_error_s_a.fh $(SRCROOT)/FMS/sat_vapor_pres/include/sat_vapor_pres_k.inc $(SRCROOT)/FMS/fms/read_data_4d.inc $(SRCROOT)/FMS/mpp/include/mpp_do_update.fh $(SRCROOT)FMS/mpp/include/mpp_update_domains2D.fh $(SRCROOT)/FMS/horiz_interp/include/horiz_interp_spherical.inc $(SRCROOT)/FMS/column_diagnostics/include/column_diagnostics_r8.fh $(SRCROOT)/FMS/mpp/include/mpp_get_boundary_ad.fh $(SRCROOT)/FMS/fms2_io/include/register_unstructured_domain_restart_variable.inc $(SRCROOT)/FMS/monin_obukhov/include/monin_obukhov_r8.fh $(SRCROOT)/FMS/mpp/include/mpp_comm.inc $(SRCROOT)FMS/mpp/include/mpp_error_s_s.fh $(SRCROOT)FMS/mpp/include/mpp_scatter.fh $(SRCROOT)/FMS/data_override/include/data_override_r4.fh $(SRCROOT)FMS/mpp/include/mpp_alltoall_mpi.fh $(SRCROOT)/FMS/fms2_io/include/unstructured_domain_write.inc
OBJ = mosaic.o diag_axis.o atmos_ocean_fluxes.o fm_yaml.o horiz_interp.o read_mosaic.o fms2_io.o time_interp_external.o mpp_pset.o fms_diag_time_reduction.o monin_obukhov.o time_interp.o diag_manager.o tracer_manager.o horiz_interp_bicubic.o horiz_interp_conserve.o astronomy.o time_interp_external2.o coupler_types.o fms_io.o column_diagnostics.o mpp_parameter.o data_override.o affinity.o libFMS.o topography.o drifters_core.o horiz_interp_type.o sat_vapor_pres.o gradient_c2l.o fms_affinity.o grid.o diag_table.o fm_util.o field_manager.o fmsconstants.o yaml_parser_binding.o mosaic2.o axis_utils2.o interpolator.o mpp_utilities.o fms_yaml_output.o mpp.o gradient.o drifters_input.o fms_stacksize.o horiz_interp_spherical.o xgrid.o fms_netcdf_unstructured_domain_io.o fms_diag_fieldbuff_update.o netcdf_io.o drifters_comm.o fms.o fmsconstantsr4.o quicksort.o fms_diag_elem_weight_procs.o drifters_io.o diag_data.o mersennetwister.o yaml_parser.o diag_output.o blackboxio.o fms_netcdf_domain_io.o platform.o fms_io_utils.o diag_integral.o tridiagonal.o yaml_output_functions.o diag_grid.o axis_utils.o block_control.o sat_vapor_pres_k.o time_manager.o get_grid_version.o constantsr4.o mpp_data.o fms_diag_bbox.o mpp_efp.o mpp_io.o fms_diag_outfield.o memutils.o drifters.o mpp_memutils.o interp.o fms_string_utils.o get_cal_time.o monin_obukhov_inter.o create_xgrid.o stock_constants.o horiz_interp_bilinear.o mpp_memuse.o mosaic_util.o random_numbers.o gaussian_topog.o fms_string_utils_binding.o amip_interp.o cloud_interpolator.o mpp_domains.o grid2.o ensemble_manager.o constants.o diag_util.o
OFF = $(SRCROOT)FMS/mpp/include/mpp_read_compressed.fh $(SRCROOT)/FMS/mpp/include/mpp_domains_define.inc $(SRCROOT)/FMS/drifters/drifters_compute_k.fh $(SRCROOT)/FMS/mpp/include/mpp_do_global_field.fh $(SRCROOT)/FMS/mosaic/grid.F90 $(SRCROOT)/FMS/diag_manager/fms_diag_bbox.F90 $(SRCROOT)/FMS/mpp/include/mpp_read_distributed_ascii.fh $(SRCROOT)/FMS/mpp/include/mpp_do_update_nest.fh $(SRCROOT)/FMS/amip_interp/amip_interp.F90 $(SRCROOT)/FMS/diag_manager/diag_output.F90 $(SRCROOT)/FMS/mpp/include/mpp_sum_mpi_ad.fh $(SRCROOT)/FMS/sat_vapor_pres/include/sat_vapor_pres_r4.fh $(SRCROOT)/FMS/mpp/include/system_clock.fh $(SRCROOT)/FMS/diag_manager/fms_diag_fieldbuff_update.F90 $(SRCROOT)/FMS/time_interp/include/time_interp.inc $(SRCROOT)/FMS/mpp/include/mpp_domains_reduce.inc $(SRCROOT)/FMS/diag_manager/diag_data.F90 $(SRCROOT)/FMS/diag_integral/include/diag_integral_r8.fh $(SRCROOT)/FMS/data_override/include/get_grid_version_r8.fh $(SRCROOT)/FMS/field_manager/parse.inc $(SRCROOT)/FMS/horiz_interp/include/horiz_interp_bicubic_r4.fh $(SRCROOT)FMS/mpp/include/mpp_chksum_scalar.fh $(SRCROOT)/FMS/monin_obukhov/include/monin_obukhov_inter_r8.fh $(SRCROOT)/FMS/amip_interp/include/amip_interp.inc $(SRCROOT)/FMS/mpp/include/mpp_global_sum.fh $(SRCROOT)/FMS/horiz_interp/include/horiz_interp_spherical_r4.fh $(SRCROOT)/FMS/mpp/include/mpp_get_boundary_ad.fh $(SRCROOT)/FMS/column_diagnostics/include/column_diagnostics_r8.fh $(SRCROOT)/FMS/parser/fms_yaml_output.F90 $(SRCROOT)/FMS/diag_manager/fms_diag_time_reduction.F90 $(SRCROOT)FMS/mpp/include/mpp_scatter.fh $(SRCROOT)/FMS/memutils/memutils.F90 $(SRCROOT)/FMS/mpp/include/mpp_comm.inc $(SRCROOT)/FMS/diag_manager/diag_util.F90 $(SRCROOT)/FMS/data_override/include/data_override_r4.fh $(SRCROOT)/FMS/mpp/include/mpp_do_update_nonblock.fh $(SRCROOT)/FMS/fms2_io/include/get_global_attribute.inc $(SRCROOT)/FMS/time_interp/include/time_interp_external2.inc $(SRCROOT)FMS/mpp/include/mpp_update_domains2D_nonblock.fh $(SRCROOT)/FMS/data_override/include/data_override.inc $(SRCROOT)/FMS/fms/fms_io_unstructured_get_file_unit.inc $(SRCROOT)FMS/mpp/include/mpp_do_global_field_ad.fh $(SRCROOT)/FMS/fms/fms_io_unstructured_get_file_name.inc $(SRCROOT)FMS/mpp/include/mpp_transmit_nocomm.fh $(SRCROOT)FMS/mpp/include/mpp_update_domains2D.fh $(SRCROOT)/FMS/mpp/mpp_pset.F90 $(SRCROOT)/FMS/sat_vapor_pres/include/sat_vapor_pres_k.inc $(SRCROOT)/FMS/mosaic/mosaic.F90 $(SRCROOT)/FMS/constants4/constantsr4.F90 $(SRCROOT)/FMS/mosaic/gradient_c2l.h $(SRCROOT)/FMS/mpp/include/mpp_unstruct_pass_data.fh $(SRCROOT)/FMS/field_manager/include/fm_util.inc $(SRCROOT)/FMS/mpp/include/mpp_write_compressed.fh $(SRCROOT)FMS/mpp/include/mpp_reduce_mpi.fh $(SRCROOT)FMS/mpp/include/mpp_reduce_nocomm.fh $(SRCROOT)/FMS/mpp/include/mpp_io_unstructured_write.inc $(SRCROOT)/FMS/mpp/include/mpp_chksum.fh $(SRCROOT)/FMS/mosaic/constant.h $(SRCROOT)/FMS/fms2_io/include/scatter_data_bc.inc $(SRCROOT)/FMS/drifters/drifters_comm.F90 $(SRCROOT)/FMS/mpp/include/mpp_global_sum_ad.fh $(SRCROOT)/FMS/exchange/stock_constants.F90 $(SRCROOT)/FMS/astronomy/include/astronomy.inc $(SRCROOT)/FMS/mpp/include/mpp_error_s_s.fh $(SRCROOT)/FMS/mpp/include/mpp_global_reduce.fh $(SRCROOT)/FMS/random_numbers/mersennetwister.F90 $(SRCROOT)/FMS/fms2_io/include/register_domain_restart_variable.inc $(SRCROOT)/FMS/column_diagnostics/include/column_diagnostics_r4.fh $(SRCROOT)/FMS/mosaic/read_mosaic.h $(SRCROOT)/FMS/mosaic2/include/grid2_r4.fh $(SRCROOT)FMS/mpp/include/mpp_write_2Ddecomp.fh $(SRCROOT)/FMS/fms/fms_io_unstructured_file_unit.inc $(SRCROOT)/FMS/fms2_io/include/register_variable_attribute.inc $(SRCROOT)/FMS/mpp/include/mpp_global_sum_tl.fh $(SRCROOT)/FMS/fms2_io/include/array_utils_char.inc $(SRCROOT)/FMS/tracer_manager/include/tracer_manager_r4.fh $(SRCROOT)/FMS/mpp/include/mpp_global_field.fh $(SRCROOT)FMS/mpp/include/mpp_unstruct_pass_data.fh $(SRCROOT)/FMS/time_interp/include/time_interp_external2_r8.fh $(SRCROOT)/FMS/fms2_io/include/unstructured_domain_read.inc $(SRCROOT)/FMS/constants4/fmsconstantsr4.F90 $(SRCROOT)/FMS/mpp/mpp_memuse.c $(SRCROOT)/FMS/mpp/include/mpp_read_distributed_ascii.inc $(SRCROOT)/FMS/mpp/include/mpp_sum_ad.inc $(SRCROOT)/FMS/drifters/drifters_input.F90 $(SRCROOT)FMS/mpp/include/mpp_global_sum_ad.fh $(SRCROOT)/FMS/mpp/include/mpp_do_updateV.fh $(SRCROOT)FMS/mpp/include/mpp_do_updateV_nonblock.fh $(SRCROOT)/FMS/astronomy/include/astronomy_r4.fh $(SRCROOT)/FMS/field_manager/include/field_manager.inc $(SRCROOT)/FMS/mosaic/read_mosaic.c $(SRCROOT)/FMS/mpp/include/mpp_type_mpi.fh $(SRCROOT)/FMS/field_manager/include/field_manager_r8.fh $(SRCROOT)FMS/mpp/include/mpp_group_update.fh $(SRCROOT)/FMS/mpp/include/mpp_do_redistribute.fh $(SRCROOT)/FMS/mpp/include/mpp_io_unstructured_read.inc $(SRCROOT)/FMS/amip_interp/include/amip_interp_r8.fh $(SRCROOT)/FMS/horiz_interp/include/horiz_interp_conserve.inc $(SRCROOT)/FMS/topography/include/topography_r8.fh $(SRCROOT)/FMS/mpp/include/mpp_do_update_ad.fh $(SRCROOT)/FMS/constants/constants.F90 $(SRCROOT)/FMS/horiz_interp/horiz_interp_conserve.F90 $(SRCROOT)/FMS/diag_manager/fms_diag_outfield.F90 $(SRCROOT)/FMS/horiz_interp/include/horiz_interp_type_r4.fh $(SRCROOT)/FMS/string_utils/include/fms_string_utils.inc $(SRCROOT)FMS/mpp/include/mpp_do_updateV_ad.fh $(SRCROOT)/FMS/horiz_interp/horiz_interp_type.F90 $(SRCROOT)/FMS/mpp/include/mpp_do_get_boundary.fh $(SRCROOT)/FMS/sat_vapor_pres/sat_vapor_pres.F90 $(SRCROOT)/FMS/mpp/include/mpp_util.inc $(SRCROOT)/FMS/fms/fms_io_unstructured_read.inc $(SRCROOT)FMS/mpp/include/mpp_global_field.fh $(SRCROOT)/FMS/astronomy/astronomy.F90 $(SRCROOT)/FMS/mpp/include/mpp_global_field_ad.fh $(SRCROOT)/FMS/time_interp/time_interp.F90 $(SRCROOT)/FMS/drifters/drifters_push.fh $(SRCROOT)/FMS/topography/include/topography_r4.fh $(SRCROOT)/FMS/mosaic/interp.c $(SRCROOT)/FMS/random_numbers/include/random_numbers.inc $(SRCROOT)/FMS/fms2_io/include/domain_write.inc $(SRCROOT)/FMS/fms2_io/include/compute_global_checksum.inc $(SRCROOT)FMS/mpp/include/mpp_read_distributed_ascii.fh $(SRCROOT)/FMS/diag_manager/include/fms_diag_fieldbuff_update.fh $(SRCROOT)/FMS/drifters/quicksort.F90 $(SRCROOT)/FMS/data_override/include/data_override_r8.fh $(SRCROOT)FMS/mpp/include/mpp_sum_mpi.fh $(SRCROOT)/FMS/random_numbers/include/random_numbers_r4.fh $(SRCROOT)/FMS/mosaic/mosaic_util.c $(SRCROOT)/FMS/mosaic2/mosaic2.F90 $(SRCROOT)FMS/mpp/include/mpp_global_sum.fh $(SRCROOT)/FMS/mpp/include/mpp_do_updateV_nonblock.fh $(SRCROOT)/FMS/mpp/include/mpp_scatter.fh $(SRCROOT)/FMS/fms/fms_io_unstructured_get_field_size.inc $(SRCROOT)/FMS/constants/gfdl_constants.fh $(SRCROOT)/FMS/time_manager/get_cal_time.F90 $(SRCROOT)FMS/mpp/include/mpp_do_redistribute.fh $(SRCROOT)/FMS/diag_manager/fms_diag_elem_weight_procs.F90 $(SRCROOT)/FMS/coupler/coupler_types.F90 $(SRCROOT)/FMS/fms/read_data_2d.inc $(SRCROOT)/FMS/affinity/fms_affinity.F90 $(SRCROOT)/FMS/fms/include/fms_r4.fh $(SRCROOT)/FMS/fms/write_data.inc $(SRCROOT)/FMS/mpp/include/mpp_alltoall_mpi.fh $(SRCROOT)/FMS/fms2_io/fms_netcdf_domain_io.F90 $(SRCROOT)/FMS/parser/yaml_output_functions.c $(SRCROOT)/FMS/include/fms_platform.h $(SRCROOT)/FMS/mpp/include/mpp_define_nest_domains.inc $(SRCROOT)/FMS/topography/include/topography.inc $(SRCROOT)FMS/include/file_version.h $(SRCROOT)FMS/mpp/include/mpp_sum_nocomm.fh $(SRCROOT)/FMS/monin_obukhov/include/monin_obukhov_r4.fh $(SRCROOT)/FMS/mosaic2/include/grid2.inc $(SRCROOT)/FMS/mpp/include/mpp_write.fh $(SRCROOT)FMS/mpp/include/group_update_pack.inc $(SRCROOT)/FMS/diag_manager/diag_manager.F90 $(SRCROOT)/FMS/affinity/affinity.c $(SRCROOT)/FMS/mpp/include/mpp_update_domains2D_ad.fh $(SRCROOT)/FMS/mosaic2/include/mosaic2_r4.fh $(SRCROOT)FMS/mpp/include/mpp_global_field_ug.fh $(SRCROOT)/FMS/tridiagonal/include/tridiagonal.inc $(SRCROOT)/FMS/string_utils/fms_string_utils.F90 $(SRCROOT)/FMS/time_interp/time_interp_external2.F90 $(SRCROOT)/FMS/mpp/include/mpp_read_compressed.fh $(SRCROOT)/FMS/mpp/include/mpp_error_a_a.fh $(SRCROOT)FMS/mpp/include/mpp_get_boundary.fh $(SRCROOT)/FMS/mpp/include/group_update_pack.inc $(SRCROOT)/FMS/mpp/include/mpp_transmit_mpi.fh $(SRCROOT)/FMS/drifters/drifters.F90 $(SRCROOT)/FMS/string_utils/include/fms_string_utils_r4.fh $(SRCROOT)/FMS/mpp/include/mpp_data_nocomm.inc $(SRCROOT)/FMS/horiz_interp/include/horiz_interp_spherical_r8.fh $(SRCROOT)FMS/mpp/include/mpp_do_update.fh $(SRCROOT)/FMS/field_manager/include/field_manager_r4.fh $(SRCROOT)/FMS/axis_utils/include/axis_utils2.inc $(SRCROOT)/FMS/mpp/mpp.F90 $(SRCROOT)/FMS/horiz_interp/include/horiz_interp_r8.fh $(SRCROOT)/FMS/time_interp/include/time_interp_r8.fh $(SRCROOT)/FMS/time_interp/include/time_interp_external2_r4.fh $(SRCROOT)/FMS/tridiagonal/include/tridiagonal_r4.fh $(SRCROOT)FMS/mpp/include/mpp_do_get_boundary_ad.fh $(SRCROOT)FMS/mpp/include/mpp_sum.inc $(SRCROOT)/FMS/fms/fms_io_unstructured_register_restart_field.inc $(SRCROOT)FMS/mpp/include/mpp_global_sum_tl.fh $(SRCROOT)FMS/mpp/include/mpp_gather.fh $(SRCROOT)/FMS/mpp/include/mpp_sum_nocomm_ad.fh $(SRCROOT)FMS/mpp/include/mpp_transmit.inc $(SRCROOT)/FMS/mpp/include/mpp_do_get_boundary_ad.fh $(SRCROOT)/FMS/fms2_io/include/netcdf_read_data.inc $(SRCROOT)/FMS/horiz_interp/include/horiz_interp_type_r8.fh $(SRCROOT)FMS/mpp/include/mpp_sum_nocomm_ad.fh $(SRCROOT)/FMS/topography/include/gaussian_topog_r8.fh $(SRCROOT)/FMS/fms2_io/include/domain_read.inc $(SRCROOT)/FMS/mpp/include/mpp_io_write.inc $(SRCROOT)/FMS/mpp/include/mpp_update_domains2D.fh $(SRCROOT)FMS/mpp/include/mpp_global_reduce.fh $(SRCROOT)/FMS/mpp/include/mpp_do_checkV.fh $(SRCROOT)/FMS/mpp/include/mpp_do_global_field_ad.fh $(SRCROOT)/FMS/horiz_interp/include/horiz_interp_conserve_r8.fh $(SRCROOT)FMS/mpp/include/mpp_type_mpi.fh $(SRCROOT)/FMS/mpp/include/mpp_util_mpi.inc $(SRCROOT)/FMS/fms2_io/include/netcdf_write_data.inc $(SRCROOT)/FMS/random_numbers/include/random_numbers_r8.fh $(SRCROOT)/FMS/mpp/include/mpp_util_nocomm.inc $(SRCROOT)/FMS/mosaic2/include/grid2_r8.fh $(SRCROOT)/FMS/time_manager/time_manager.F90 $(SRCROOT)/FMS/diag_manager/diag_table.F90 $(SRCROOT)/FMS/mpp/include/mpp_sum.inc $(SRCROOT)/FMS/tracer_manager/include/tracer_manager.inc $(SRCROOT)/FMS/axis_utils/include/axis_utils2_r8.fh $(SRCROOT)FMS/mpp/include/mpp_error_s_a.fh $(SRCROOT)/FMS/monin_obukhov/include/monin_obukhov.inc $(SRCROOT)/FMS/topography/include/gaussian_topog.inc $(SRCROOT)/FMS/field_manager/field_manager.F90 $(SRCROOT)/FMS/mpp/include/mpp_io_misc.inc $(SRCROOT)/FMS/topography/topography.F90 $(SRCROOT)/FMS/block_control/block_control.F90 $(SRCROOT)/FMS/diag_manager/include/fms_diag_fieldbuff_update.inc $(SRCROOT)/FMS/field_manager/include/fm_util_r8.fh $(SRCROOT)/FMS/coupler/include/coupler_types_r4.fh $(SRCROOT)FMS/mpp/include/mpp_do_updateV.fh $(SRCROOT)/FMS/tracer_manager/tracer_manager.F90 $(SRCROOT)/FMS/fms2_io/include/get_variable_attribute.inc $(SRCROOT)/FMS/tridiagonal/tridiagonal.F90 $(SRCROOT)FMS/mpp/include/mpp_sum_mpi_ad.fh $(SRCROOT)/FMS/fms2_io/include/register_global_attribute.inc $(SRCROOT)/FMS/interpolator/interpolator.F90 $(SRCROOT)/FMS/constants4/geos_constantsR4.fh $(SRCROOT)/FMS/sat_vapor_pres/include/sat_vapor_pres_k_r8.fh $(SRCROOT)/FMS/include/file_version.h $(SRCROOT)/FMS/mpp/include/mpp_io_read.inc $(SRCROOT)/FMS/data_override/get_grid_version.F90 $(SRCROOT)/FMS/drifters/fms_switches.h $(SRCROOT)/FMS/interpolator/include/interpolator.inc $(SRCROOT)/FMS/mpp/include/mpp_sum_mpi.fh $(SRCROOT)/FMS/drifters/drifters_io.F90 $(SRCROOT)/FMS/mosaic2/include/mosaic2_r8.fh $(SRCROOT)/FMS/mpp/include/mpp_unstruct_domain.inc $(SRCROOT)/FMS/fms2_io/include/register_unstructured_domain_restart_variable.inc $(SRCROOT)/FMS/mosaic/create_xgrid.c $(SRCROOT)/FMS/horiz_interp/include/horiz_interp_spherical.inc $(SRCROOT)/FMS/drifters/cloud_interpolator.F90 $(SRCROOT)/FMS/mosaic/gradient.F90 $(SRCROOT)/FMS/monin_obukhov/include/monin_obukhov_r8.fh $(SRCROOT)/FMS/constants/fmsconstants.F90 $(SRCROOT)FMS/mpp/include/mpp_error_s_s.fh $(SRCROOT)/FMS/fms2_io/include/unstructured_domain_write.inc $(SRCROOT)FMS/mpp/include/mpp_alltoall_mpi.fh $(SRCROOT)/FMS/mpp/include/mpp_reduce_nocomm.fh $(SRCROOT)/FMS/topography/include/gaussian_topog_r4.fh $(SRCROOT)/FMS/parser/yaml_parser_binding.c $(SRCROOT)/FMS/string_utils/include/fms_string_utils_r8.fh $(SRCROOT)/FMS/fms2_io/include/get_data_type_string.inc $(SRCROOT)/FMS/mosaic/interp.h $(SRCROOT)/FMS/fms2_io/netcdf_io.F90 $(SRCROOT)/FMS/mpp/include/mpp_do_update.fh $(SRCROOT)/FMS/fms/read_data_4d.inc $(SRCROOT)/FMS/mpp/include/mpp_error_s_a.fh $(SRCROOT)/FMS/mpp/include/mpp_chksum_int.fh $(SRCROOT)/FMS/constants4/gfdl_constantsR4.fh $(SRCROOT)/FMS/monin_obukhov/monin_obukhov.F90 $(SRCROOT)/FMS/mpp/mpp_memutils.F90 $(SRCROOT)/FMS/platform/platform.F90 $(SRCROOT)/FMS/horiz_interp/include/horiz_interp_bilinear_r4.fh $(SRCROOT)/FMS/horiz_interp/include/horiz_interp_bicubic_r8.fh $(SRCROOT)/FMS/mpp/include/mpp_group_update.fh $(SRCROOT)/FMS/horiz_interp/include/horiz_interp_type.inc $(SRCROOT)/FMS/mpp/include/mpp_error_a_s.fh $(SRCROOT)/FMS/horiz_interp/horiz_interp.F90 $(SRCROOT)/FMS/mpp/include/group_update_unpack.inc $(SRCROOT)/FMS/coupler/atmos_ocean_fluxes.F90 $(SRCROOT)/FMS/mpp/include/mpp_gather.fh $(SRCROOT)/FMS/data_override/include/get_grid_version.inc $(SRCROOT)/FMS/fms/fms_io_unstructured_save_restart.inc $(SRCROOT)/FMS/column_diagnostics/column_diagnostics.F90 $(SRCROOT)FMS/mpp/include/mpp_error_a_s.fh $(SRCROOT)/FMS/mosaic/create_xgrid.h $(SRCROOT)/FMS/diag_manager/diag_grid.F90 $(SRCROOT)/FMS/fms/fms_io.F90 $(SRCROOT)/FMS/tridiagonal/include/tridiagonal_r8.fh $(SRCROOT)/FMS/horiz_interp/include/horiz_interp_r4.fh $(SRCROOT)/FMS/mpp/include/mpp_write_unlimited_axis.fh $(SRCROOT)/FMS/mpp/mpp_efp.F90 $(SRCROOT)FMS/mpp/include/mpp_read_2Ddecomp.fh $(SRCROOT)FMS/mpp/include/mpp_do_check.fh $(SRCROOT)/FMS/parser/yaml_parser.F90 $(SRCROOT)FMS/mpp/include/mpp_update_nest_domains.fh $(SRCROOT)/FMS/data_override/data_override.F90 $(SRCROOT)FMS/mpp/include/mpp_alltoall_nocomm.fh $(SRCROOT)/FMS/sat_vapor_pres/sat_vapor_pres_k.F90 $(SRCROOT)/FMS/axis_utils/include/axis_utils2_r4.fh $(SRCROOT)/FMS/mosaic2/grid2.F90 $(SRCROOT)/FMS/mpp/mpp_utilities.F90 $(SRCROOT)/FMS/mpp/mpp_data.F90 $(SRCROOT)/FMS/mpp/include/mpp_data_mpi.inc $(SRCROOT)/FMS/mpp/include/mpp_get_boundary.fh $(SRCROOT)FMS/mpp/include/mpp_get_boundary_ad.fh $(SRCROOT)/FMS/fms/fms_stacksize.c $(SRCROOT)/FMS/field_manager/fm_yaml.F90 $(SRCROOT)/FMS/topography/gaussian_topog.F90 $(SRCROOT)/FMS/sat_vapor_pres/include/sat_vapor_pres_k_r4.fh $(SRCROOT)FMS/mpp/include/mpp_do_update_nest.fh $(SRCROOT)/FMS/tracer_manager/include/tracer_manager_r8.fh $(SRCROOT)/FMS/mpp/include/mpp_global_field_ug.fh $(SRCROOT)/FMS/coupler/include/coupler_types.inc $(SRCROOT)/FMS/mpp/include/mpp_comm_mpi.inc $(SRCROOT)/FMS/interpolator/include/interpolator_r4.fh $(SRCROOT)/FMS/fms/read_data_3d.inc $(SRCROOT)/FMS/fms2_io/include/compressed_read.inc $(SRCROOT)FMS/mpp/include/mpp_do_update_nonblock.fh $(SRCROOT)/FMS/mpp/include/mpp_sum_nocomm.fh $(SRCROOT)/FMS/monin_obukhov/include/monin_obukhov_inter.inc $(SRCROOT)/FMS/horiz_interp/horiz_interp_spherical.F90 $(SRCROOT)/FMS/mpp/include/mpp_transmit_nocomm.fh $(SRCROOT)/FMS/mpp/include/mpp_update_domains2D_nonblock.fh $(SRCROOT)/FMS/drifters/drifters_core.F90 $(SRCROOT)/FMS/mosaic/gradient_c2l.c $(SRCROOT)FMS/mpp/include/mpp_type_nocomm.fh $(SRCROOT)/FMS/mpp/mpp_parameter.F90 $(SRCROOT)/FMS/fms/include/fms_r8.fh $(SRCROOT)/FMS/mpp/mpp_domains.F90 $(SRCROOT)/FMS/fms/fms_io_unstructured_setup_one_field.inc $(SRCROOT)/FMS/diag_integral/include/diag_integral_r4.fh $(SRCROOT)FMS/mpp/include/mpp_do_get_boundary.fh $(SRCROOT)/FMS/sat_vapor_pres/include/sat_vapor_pres.inc $(SRCROOT)/FMS/mpp/include/mpp_reduce_mpi.fh $(SRCROOT)/FMS/random_numbers/random_numbers.F90 $(SRCROOT)/FMS/time_interp/include/time_interp_r4.fh $(SRCROOT)/FMS/horiz_interp/include/horiz_interp_bilinear_r8.fh $(SRCROOT)/FMS/time_interp/time_interp_external.F90 $(SRCROOT)/FMS/column_diagnostics/include/column_diagnostics.inc $(SRCROOT)/FMS/drifters/drifters_set_field.fh $(SRCROOT)/FMS/mpp/include/mpp_write_2Ddecomp.fh $(SRCROOT)/FMS/interpolator/include/interpolator_r8.fh $(SRCROOT)/FMS/fms2_io/include/unpack_data.inc $(SRCROOT)/FMS/string_utils/fms_string_utils_binding.c $(SRCROOT)/FMS/fms/fms_io_unstructured_field_exist.inc $(SRCROOT)/FMS/diag_manager/diag_axis.F90 $(SRCROOT)/FMS/exchange/xgrid.F90 $(SRCROOT)FMS/mpp/include/mpp_chksum_int.fh $(SRCROOT)/FMS/horiz_interp/include/horiz_interp_bicubic.inc $(SRCROOT)/FMS/mpp/include/mpp_io_connect.inc $(SRCROOT)/FMS/mpp/include/mpp_type_nocomm.fh $(SRCROOT)/FMS/sat_vapor_pres/include/sat_vapor_pres_r8.fh $(SRCROOT)/FMS/mpp/mpp_io.F90 $(SRCROOT)/FMS/mpp/include/mpp_domains_util.inc $(SRCROOT)/FMS/monin_obukhov/monin_obukhov_inter.F90 $(SRCROOT)/FMS/mpp/include/mpp_update_nest_domains.fh $(SRCROOT)FMS/mpp/include/mpp_chksum.fh $(SRCROOT)/FMS/amip_interp/include/amip_interp_r4.fh $(SRCROOT)/FMS/fms2_io/include/netcdf_add_restart_variable.inc $(SRCROOT)/FMS/astronomy/include/astronomy_r8.fh $(SRCROOT)FMS/mpp/include/mpp_sum_ad.inc $(SRCROOT)FMS/mpp/include/mpp_do_checkV.fh $(SRCROOT)/FMS/mosaic/mosaic_util.h $(SRCROOT)/FMS/horiz_interp/horiz_interp_bicubic.F90 $(SRCROOT)FMS/mpp/include/mpp_do_update_ad.fh $(SRCROOT)/FMS/field_manager/fm_util.F90 $(SRCROOT)/FMS/fms2_io/fms_io_utils.F90 $(SRCROOT)/FMS/diag_integral/diag_integral.F90 $(SRCROOT)/FMS/mpp/include/mpp_domains_misc.inc $(SRCROOT)/FMS/constants4/gfs_constantsR4.fh $(SRCROOT)/FMS/mpp/include/mpp_comm_nocomm.inc $(SRCROOT)FMS/mpp/include/mpp_global_field_ad.fh $(SRCROOT)FMS/mpp/include/group_update_unpack.inc $(SRCROOT)/FMS/diag_integral/include/diag_integral.inc $(SRCROOT)/FMS/fms2_io/include/compressed_write.inc $(SRCROOT)/FMS/mpp/include/mpp_alltoall_nocomm.fh $(SRCROOT)/FMS/libFMS.F90 $(SRCROOT)/FMS/mpp/include/mpp_chksum_scalar.fh $(SRCROOT)FMS/mpp/include/mpp_do_global_field.fh $(SRCROOT)/FMS/coupler/include/coupler_types_r8.fh $(SRCROOT)/FMS/mpp/include/mpp_read_2Ddecomp.fh $(SRCROOT)FMS/mpp/include/mpp_write_compressed.fh $(SRCROOT)/FMS/monin_obukhov/include/monin_obukhov_inter_r4.fh $(SRCROOT)/FMS/fms/fms.F90 $(SRCROOT)/FMS/mpp/include/mpp_transmit.inc $(SRCROOT)/FMS/fms/fms_io_unstructured_register_restart_axis.inc $(SRCROOT)/FMS/horiz_interp/include/horiz_interp_bilinear.inc $(SRCROOT)FMS/mpp/include/mpp_transmit_mpi.fh $(SRCROOT)/FMS/fms2_io/fms_netcdf_unstructured_domain_io.F90 $(SRCROOT)/FMS/mpp/include/mpp_io_util.inc $(SRCROOT)/FMS/horiz_interp/include/horiz_interp.inc $(SRCROOT)FMS/mpp/include/mpp_error_a_a.fh $(SRCROOT)/FMS/mpp/include/mpp_domains_comm.inc $(SRCROOT)/FMS/fms2_io/include/gather_data_bc.inc $(SRCROOT)/FMS/fms2_io/include/array_utils.inc $(SRCROOT)/FMS/mosaic2/include/mosaic2.inc $(SRCROOT)/FMS/mpp/include/mpp_do_check.fh $(SRCROOT)/FMS/constants/gfs_constants.fh $(SRCROOT)FMS/mpp/include/mpp_update_domains2D_ad.fh $(SRCROOT)/FMS/field_manager/include/fm_util_r4.fh $(SRCROOT)/FMS/coupler/ensemble_manager.F90 $(SRCROOT)/FMS/axis_utils/axis_utils.F90 $(SRCROOT)/FMS/constants/geos_constants.fh $(SRCROOT)/FMS/data_override/include/get_grid_version_r4.fh $(SRCROOT)/FMS/axis_utils/axis_utils2.F90 $(SRCROOT)/FMS/horiz_interp/include/horiz_interp_conserve_r4.fh $(SRCROOT)FMS/mpp/include/mpp_write.fh $(SRCROOT)/FMS/mpp/include/mpp_do_updateV_ad.fh $(SRCROOT)/FMS/fms/include/fms.inc $(SRCROOT)/FMS/fms2_io/fms2_io.F90 $(SRCROOT)/FMS/horiz_interp/horiz_interp_bilinear.F90 $(SRCROOT)/FMS/fms2_io/blackboxio.F90
clean: neat
	-rm -f .libFMS.a.cppdefs $(OBJ) libFMS.a

neat:
	-rm -f $(TMPFILES)

localize: $(OFF)
	cp $(OFF) .

TAGS: $(SRC)
	etags $(SRC)

tags: $(SRC)
	ctags $(SRC)

libFMS.a: $(OBJ)
	$(AR) $(ARFLAGS) libFMS.a $(OBJ)

