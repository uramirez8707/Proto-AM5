# RELEASE 2024.03
This 2024.03 is **not expected** to reproduce previous releases because of the inclusion of aerosol wet scavenging by ice flux.
1. micro_mg2.F90 has been updated to include aerosol wet scavenging by ice flux.


# RELEASE 2024.02
This 2024.02 is **not expected** to reproduce previous releases because of the changes to sulfate, DMS and SOA.
1. DMS chemistry is updated so that DMS is conserved and pH is replaced by H_cloud in atmos_sulfate.F90,
2. Diurnal variation for radicals is fixed in atmos_soa.F90

# RELEASE 2024.01
This release contains a new PBL scheme (NCEP TKE-based EDMF) along with other changes. This 2024.01 version is **expected** to reprod 2023.07 if the default namelist parameters are used.
1. The NCEP TKE-based EDMF package is incorporated under the path atmos_param/han_edmf.
2. Various changes are made to the other modules including physics_driver, vert_turb_driver, shallow_cu (UW convection scheme) for consistency with NCEP EDMF.
3. Additional changes related to the Tiedtke cloud fraction scheme are made for future work on PBL-convection-cloud coupling.
4. All updates above are disabled with the default namelist parameters, and thus the default runs reproduce 2023.07.

# RELEASE 2023.07
This release contains namelist pruning for the sea salt module. This 2023.07 version is **NOT** expected to reprod 2023.06 due to the changes in namelist.
1. The default value for the namelist `ssalt_nml::use_zieger_growth` was changed from `.FALSE.` to `.TRUE.`
2. The namelist `ssalt_nml::min_scale_marthenson` has been removed. The code will use a `min_scale_marthenson` of 0
3. The namelist `ssalt_nml::use_tsurf_for_scaling` have been removed. The code will follow the `use_tsurf_for_scaling=.true.` path

# RELEASE 2023.06
This release is to fix the negative vapor in cloud macro- and micro-physics. This 2023.06 version is expected to reprod 2023.04 if the default values of the new namelist parameters are used.
1. A maximum limit for the large-scale condensation/depostion was added in ls_cloud_microphysics.F90.


# RELEASE 2023.05
This release contains updated to the MG2 cloud microphysics. This 2023.05 version is expected to reprod 2023.04 if the default values of the new namelist parameters are used.
1. A relative humidity threshold for rain evaporation and snow sublimation was added for MG2 cloud microphysics
2. A nml variable "debug_cld_microphysics" was added to turn on/off debugging printout.

# RELEASE 2023.04
This release brings in updates from atmos_phys. This 2023.04 version is expected to reproduce 2023.03. 
1. The sim_data_flsp and beta distribution files were updated to be read with Fortran's open, read, and close
2. These updates are going to be needed to run with the upcoming version of FMS (2023.02)

# RELEASE 2023.03
This release contains chemistry namelist and code pruning. This 2023.03 version is not expected to reproduce 2023.02 answers due to the removal of bug flags.

1. The following ununsed macros have been removed:
    - _ALLOC
    - USE_MPI
    - test_aerosol
2. The following namelists were updated:
    - **atmos_tracer_driver_nml**:  
        - *prevent_flux_through_ice* default value has been changed from `false` to `.true`
        - *do_cmip6_bug_diag* has been removed following the `.false.` path (The default value had been set to .true.)
    - **tropchem_driver_nml**
        - *solar_flux_bugfix* has been removed following the `.true.` path (The default value had been set to .false.)
        - *het_chem_bug1* has been removed following the `.false.` path (The default value had been set to .true.)
    - **ice_nucl_nml**
        - *retain_ice_nucl_bug* has been removed following the `.true.` path (The default value had been set to .true.)
    - **atmos_nh3_tag_nml**
        - *debug* default value has been changed from `true` to `false`

# RELEASE 2023.02
This release contains namelist pruning, and the addition of a new linear ozone scheme. This 2023.02 version does not change AM5 answers, but it will require removing the pruned namelist from the xml.

1. Namelist pruning:
    - The default values of nml variables have been changed to match those specified in the AM5 xml. Details can be found [here](https://docs.google.com/spreadsheets/d/1K4A6FXj-vv056AFFox_qllfzHfLtfEh8Og4QkrVLyLA/edit#gid=0)
    - The following namelist variables have been pruned:
        <details><summary> rotstayn_klein_nml</summary>
        do_old_snowmelt (code now follows the .false. path) <br/>
        retain_cm3_bug (code now follows the .false. path) <br/>
        use_inconsistent_lh (code now follows the .true. path) <br/>
        </details>
        <details><summary> moist_processes_nml</summary>
        do_simple  (code now follows the .false. path) <br/>
        </details>
        <details><summary> uw_conv_nml</summary>
        reproduce_old_version  (code now follows the .false. path) <br/>
        </details>
        <details><summary> cloudrad_diagnostics_nml</summary>
        do_outdated_isccp (code now follows the .false. path) <br/>
        </details>
        <details><summary> cloud_spec_nml</summary>
        reproduce_ulm (code now follows the .true. path) <br/>
        </details>
        <details><summary> get_random_number_stream_nml</summary>
        do_legacy_seed_generation (code now follows the .false. path)  <br/>
        </details>
        <details><summary> microphys_rad_nml</summary>
        remain_hu_bug (code now follows the .false. path)  <br/>
        </details>
        <details><summary> uw_clouds_W_nml</summary>
        preserve_inconsistency (code now follows the .false. path)  <br/>
        </details>
        <details><summary> lw_gases_stdtf_nml</summary>
        do_co2_bug (code now follows the .false. path)  <br/>
        </details>
        <details><summary> esfsw_driver_nml</summary>
        reproduce_ulm (code now follows the .false. path) <br/>
        remain_rayleigh_bug (code now follows the .false. path) <br/>
        </details>
        <details><summary> lscloud_driver_nml</summary>
        use_cf_metadata (code now follows the .true. path) <br/>
        </details>
2. New linear ozone scheme is added:
    - This is a option to represent stratospheric ozone that is parallel to prescribing ozone concentration or full stratospheric chemistry. In this scheme, ozone is treated as a tracer with linearized chemical tendency. By default this option is turned off. See https://gitlab.gfdl.noaa.gov/fms/am5_phys/-/merge_requests/48 for details on how to turn this option on.

# RELEASE 2023.01
This release contains minor chemistry pruning updates, updates brought in from atmos_phys, and updates MG2 microphysics. This 2023.01 version does not change AM5 answers if **MG2 microphysics is not turned on**.

1. Chemistry pruning:
    -  All `AM3_CHEM macro` statements have been removed.
    -  The tracer_driver/tropchem/AM3_chem directory has been removed.
    -  The namelist, `tropchem_driver_nml :: retain_cm3_bugs` and the buggy `retain_cm3_bugs=.true.` code have been removed
    -  The legacy `f1p_bug` and `f1p_bug2` options have been removed from the cloud_chem_solver
    -  The `CLOUD_CHEM_LEGACY`, `CLOUD_CHEM_F1P_BUG`, AND `CLOUD_CHEM_F1P_BUG2` options have been removed from the `cloud_chem_type` derived type
    -  The optional argument, `do_am3_bug`, for subroutine cloud_so2_chem has been removed
    -  The default nml values for `cloud_chem_type` and `cloud_chem_solver` have been changed to `f1p`
    -  The html markups in `atmos_param` have been changed to regular comment sections.
    -  All `_ALLOCATABLE` and `_ALLOCATED` usage have been changed to allocatable and allocated.
    -  The `_NULL` macro has been removed.
    -  `qsrout_3d` and `qrout_3dvariables` associated with the pruned MG scheme have been removed
    -  The namelist, `ssalt_nml :: ulm_ssalt_deposition` and `ssalt_nml :: ssalt_debug` have been removed
2.  Atmos_phys updates:
    - https://gitlab.gfdl.noaa.gov/fms/atmos_phys/-/commit/709fdd120bb7e260474a65146a8c76ad7ee9f3ec Mixed mode updates (needed for FMS 2022.03 + versions)
    - https://gitlab.gfdl.noaa.gov/fms/atmos_phys/-/commit/46705feafd85940faafefd4cf715ee92f1ece0bc Bug fix for CH4 diagnostics (changes AM4.1/2 answers)
    - https://gitlab.gfdl.noaa.gov/fms/atmos_phys/-/commit/1f959c4953fde06cecc0f30b47c99058d32d6d80 Update to receive alkalinity deposition flux from Ocean_BGC
    - https://gitlab.gfdl.noaa.gov/fms/atmos_phys/-/commit/bdf29bdcc7684650e2b58bef95d03bcd6e128fae Bug fix for biogenic terpene emissions (changes AM4.1/2 answers unless bug flag is specified)
    - https://gitlab.gfdl.noaa.gov/fms/atmos_phys/-/commit/bbbb4f985439755ccebe67cf792576bcd44a1d5d Bug in MEGAN BVOC routines
    - https://gitlab.gfdl.noaa.gov/fms/atmos_phys/-/commit/bbbb4f985439755ccebe67cf792576bcd44a1d5d Fixes a bug causing memory leak/crash related to use of horiz_interp_new
3.  MG2 microphysics updates:
    - MG2 microphysics was updated to strictly conserve total mass and enthalpy.
    - Changes answers if MG2 microphysics is turned on

# RELEASE 2022.01
In this release, portions of `atmos_param` that are unused in the `AM4 model` have been pruned.  This 2022.01 version reproduces all the CHECKSUMS and restart numerical data for AM4, CM4, ESM4, and SPEAR.  Users will however see that 
1.  the variable **restart_version** no longer exists in cg_drag.res.nc 
2.  the variable **doing_edt** no longer exists in physics_driver.res.nc
3.  the variable **tke** no longer exists in physics_driver.res.tile[1-6].nc

#### DIRECTORIES
In details, the following directories, referred to here as they are referred in the code, have been pruned:
- betts_miller 
- cloud_obs
- clouds
- cloud_zonal
- CLUBB
- diag_cloud
- diag_cloud_radiation
- diffusivity
- donner_deep
- dry_adj
- edt
- grey_radiation
- lin_cloud_microphysics
- my25_turb
- qe_moist_convection
- ras
- rh_clouds
- shallow_conv
- shallow_physics
- strat_cloud
- tke_turb
- two_stream_gray_rad 

### FILES
In addition, the following files (modules) have been pruned
- bulkphys_rad_mod
- cldwat2m_micro_mod
- donner_deep_clouds_W_mod
- micro_mg_mod
- morrison_gettelman_microp_mod
- simple_pdf_mod 

### NML VARIABLES
Removal of the abovementioned directories and files required the following nml options to be pruned:
<details><summary> cg_drag_nml</summary>
        calculate_ked <br/>
        num_diag_pts_ij <br/>
        num_diag_pts_latlon <br/>
        i_coords_gl <br/>
        j_coords_gl <br/>
        lat_coords_gl <br/>
        lon_coords_gl <br/>
        Bt_eq<br/>
        Bt_eq_width<br/>
</details>

<details><summary> cloud_rad_nml</summary>
        clubb_error <br />
        prog_ccn
</details>

<details><summary> convection_driver_nml</summary>
       do_limit_donner <br/> 
       do_unified_convective_cloure <br/>
       do_donner_before_uw <br/>
       use_updated_profiles_for_uw <br/>
       use_updated_profiles_for_donner <br/>
       only_one_conv_scheme_per_column <br/>
       force_donner_moist_conserv <br/>
       do_donner_conservation_checks <br/>
       do_donner_mca <br/> 
       conv_frac_max <br/> 
       cmt_mass_flux_source = 'donner', 'donner_and_ras', 'donner_and_uw', 'ras_and_uw', 'donner_and_ras_anduw'
       remain_detrain_bug <br/> 
       keep_icenum_detrain_bug <br/>
</details>

<details><summary> ls_cloud_driver_nml </summary>
        do_legacy_strat_cloud <br/>
        microphys_scheme = 'lin','morrison_gettelman', 'mg_ncar', 'ncar'
</details>

<details><summary> ls_cloud_macrophysics_nml </summary>
        use_updated_profiles_for_clubb
</details>

<details><summary> ls_coud_microphysics_nml </summary>
        lin_microphys_top_press <br/>
        override_liq_num <br/>
        override_ice_num <br/> 
        use_Meyers <br/> 
        use_Cooper <br/> 
        micro_begin_sec <br/> 
        min_precip_needing_adjustment
</details>

<details><summary> rotstayn_klein_mp_nml</summary>
       use_inconsistent_lh
</details>

<details><summary> moist_processes_nml</summary>
        do_mca <br/> 
        do_ras <br/>
        do_donner_deep <br/> 
        do_dryadj <br/> 
        do_bm <br/> 
        do_bmmass <br/> 
        do_bmomp <br/> 
        do_simple <br/> 
        do_rh_clouds <br/> 
        include_donmca_in_cosp
</details>

<details><summary> physics_driver_nml </summary>
        do_clubb <br/>
        donner_meso_is_largescale <br/> 
        do_grey_radiation <br/> 
        R1, R2, R3, R4 <br/> 
        l_host_applies_sfc_fluxes
</details>

<details><summary> cloud_spec_nml </summary>
        ignore_donner_cells <br/>
        cloud_type_form = 'rh', 'deep', 'stratdeep', or 'stratdeepuw', 'deepuw'
</details>

<details><summary> microphys_rad_nml </summary>
        do_orig_donner_stoch <br/> 
        ignore_donner_cells
</details>
    
<details><summary> uw_conv_nml </summary>
        use_turb_tke
</details>
  
<details><summary> vert_turb_driver_nml </summary>
        do_shallow_conv <br/> 
        do_mellor_yamada <br/> 
        do_tke_turb <br/> 
        do_diffusivity <br/> 
        do_molecular_diffusion <br/>
        do_edt <br/> 
        do_simple
</details>
