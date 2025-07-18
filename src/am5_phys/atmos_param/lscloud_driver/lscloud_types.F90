                    module lscloud_types_mod

use  fms_mod,        only : write_version_number
use  atmos_cmip_diag_mod, only : cmip_diag_id_type

implicit none
private

!-------------------------------------------------------------------------
!--interfaces-------------------------------------------------------------

public lscloud_types_init

public diag_id_type, diag_pt_type, lscloud_debug_type, &
       lscloud_nml_type, atmos_state_type, particles_type,   &
       cloud_state_type, precip_state_type, cloud_processes_type,    &
       lsc_constants_type

!----------------------------------------------------------------------
!----version number----------------------------------------------------

Character(len=128) :: Version = '$Id:  $'
Character(len=128) :: Tagname = '$Name: $'

logical  :: module_is_initialized = .false.


!########################################################################

TYPE diag_id_type

!  cloud area variables

  integer :: aall, aliq, aice, cf_liq_init, cf_ice_init, aauto
  integer :: SA3d, qadt_lsform, qadt_lsdiss, qadt_rhred, qadt_eros,  &
             qadt_fill, qadt_super, qadt_destr, qadt_limits, qadt_ahuco, &
             SA_imb
  integer :: SA2d, qa_lsform_col, qa_lsdiss_col, qa_rhred_col,  &
             qa_eros_col, qa_fill_col, qa_super_col, qa_destr_col,     &
             qa_limits_col, qa_ahuco_col, SA_imb_col

!  cloud liquid variables

  integer :: SL3d, qldt_cond, qldt_evap, qldt_eros, qldt_berg, qldt_freez,&
             liq_adj, qldt_rime, qldt_accr, qldt_auto, qldt_fill, qldt_tiny, &
             qldt_destr, qldt_freez2, qldt_sedi, qldt_accrs, qldt_bergs, &
             qldt_HM_splinter, SL_imb
  integer :: SL2d, ql_cond_col, ql_evap_col, ql_eros_col, ql_berg_col,   &
             ql_freez_col, liq_adj_col, ql_rime_col, ql_accr_col,  &
             ql_auto_col, ql_fill_col, ql_tiny_col, ql_destr_col, ql_freez2_col,  &
             ql_sedi_col, ql_accrs_col, ql_bergs_col, ql_HM_splinter_col, &
             SL_imb_col

!  cloud ice variables

  integer :: SI3d, qidt_dep, qidt_subl, qidt_fall, qidt_eros, qidt_melt, &
             qidt_melt2, qidt_fill, qidt_tiny, qidt_destr, qidt_qvdep, qidt_auto,  &
             qidt_accr, qidt_accrs, ice_adj,  qidt_rain2ice, SI_imb
  integer :: SI2d, qi_dep_col, qi_subl_col, qi_fall_col, qi_eros_col,  &
             qi_melt_col, qi_melt2_col, qi_fill_col, qi_tiny_col, qi_destr_col,  &
             qi_qvdep_col, qi_auto_col, qi_accr_col, qi_accrs_col,  &
           ice_adj_col, qi_rain2ice_col,  SI_imb_col

!  rain variables
  integer :: qrdt_fill, qr_fill_col,  qrdt_destr, qr_destr_col, qrdt_tiny, qr_tiny_col

! snow variables
  integer :: qsdt_fill, qs_fill_col, qsdt_destr, qs_destr_col, qsdt_tiny, qs_tiny_col

!  rain number variables
  integer :: qnrdt_fill, qnr_fill_col, qnrdt_destr, qnr_destr_col, qnrdt_tiny, qnr_tiny_col

! snow number variables
  integer :: qnsdt_fill, qns_fill_col, qnsdt_destr, qns_destr_col, qnsdt_tiny, qns_tiny_col

  integer :: SR3d, SR2d, SNR3d, SNR2d, SS3d, SS2d, SNS3d, SNS2d

!  cloud droplet variables

  integer :: droplets_col250, gb_droplets_col, potential_droplets, &
             droplets, droplets_wtd, ql_wt, droplets_col, rvolume
  integer :: SN3d, qndt_cond , qndt_evap, qndt_fill, qndt_tiny, qndt_berg,  qndt_rime, &  !h1g, 2014-07-24
             qndt_destr, qndt_super, qndt_freez, qndt_sacws, qndt_sacws_o, &
             qndt_eros, qndt_pra, qndt_auto, qndt_nucclim, qndt_sedi, &
             qndt_melt, qndt_ihom, qndt_size_adj, qndt_fill2,   &
             qndt_contact_frz, qndt_cleanup, qndt_cleanup2, SN_imb
  integer :: SN2d, qn_cond_col, qn_evap_col, qn_fill_col, qn_tiny_col, qn_berg_col, qn_rime_col, & !h1g, 2014-07-24
             qn_destr_col, qn_super_col, qn_freez_col, qn_sacws_col,  &
             qn_sacws_o_col, qn_eros_col, qn_pra_col, qn_auto_col,    &
             qn_nucclim_col, qn_sedi_col, qn_melt_col, qn_ihom_col,  &
             qn_size_adj_col, qn_fill2_col, qn_contact_frz_col,  &
             qn_cleanup_col, qn_cleanup2_col, SN_imb_col

!  cloud ice particle variables

  integer :: nice, nice_col, gb_nice_col, potential_crystals
  integer :: SNi3d, qnidt_fill, qnidt_tiny, qnidt_nnuccd, qnidt_nsubi,  &
             qnidt_nerosi, qnidt_auto, qnidt_accr, qnidt_nprci, qnidt_nprai, &
             qnidt_nucclim1, qnidt_nucclim2, qnidt_sedi,             &
             qnidt_melt, qnidt_size_adj, qnidt_fill2,                &
             qnidt_super, qnidt_ihom, qnidt_destr, qnidt_rain2ice,   &
             qnidt_cleanup, qnidt_cleanup2, qnidt_nsacwi, SNi_imb
  integer :: SNi2d, qni_fill_col, qni_tiny_col, qni_nnuccd_col, qni_nsubi_col, &
             qni_nerosi_col, qni_auto_col, qni_accr_col, qni_nprci_col, qni_nprai_col, &
             qni_nucclim1_col, qni_nucclim2_col, qni_sedi_col, &
             qni_melt_col, qni_size_adj_col, qni_fill2_col, &
             qni_super_col, qni_ihom_col, qni_destr_col, qni_rain2ice_col, &
             qni_cleanup_col, qni_cleanup2_col,                    &
             qni_nsacwi_col, SNi_imb_col

!  aerosol diagnostics

  integer :: delta_cf, sulfate, seasalt_sub, seasalt_sup, om,   &
             rhcrit, rhcrit_min, rhiin, rhlin, cfin, imass7,     &
             ni_dust, ni_sulf, ni_bc, ndust1, ndust2, ndust3,  &
             ndust4, ndust5, dust_berg_flag, subgrid_w_variance

!  rain diagnostics

  integer :: rain3d, qrout, rain_clr, rain_cld, a_rain_clr, a_rain_cld, &
             rain_evap, rain_freeze, srfrain_accrs, srfrain_freez,  &
             srfrain_evap, rain_inst, rain_sedi, &
             rain_evap_col, rain_freeze_col,  &
             srfrain_accrs_col, srfrain_freez_col, srfrain_evap_col, &
             rain_inst_col, rain_sedi_col, &
             rain_mass_conv, rain_imb, rain_imb_col, cld_liq_imb,  &
             cld_liq_imb_col, neg_rain, qrout_col

!  rain number diagnostics
  integer :: rain_num_inst,     rain_num_sedi,     rain_num_adj,      rain_num2snow, &
             rain_num_evap,     rain_num_freez,    rain_num_selfcoll,                &
             rain_num_inst_col, rain_num_sedi_col, rain_num_adj_col,  rain_num2snow_col, &
             rain_num_evap_col, rain_num_freez_col,rain_num_selfcoll_col

!  snow diagnostics

  integer :: snow3d, qsout, snow_clr, snow_cld, a_snow_clr, a_snow_cld, &
             snow_melt, snow_inst, snow_sedi, &
             snow_melt_col, snow_mass_conv, sedi_ice, snow_imb, &
             snow_inst_col, snow_sedi_col, &
             snow_imb_col, cld_ice_imb, cld_ice_imb_col, neg_snow, qsout_col

!  snow number diagnostics
  integer :: snow_num_inst,     snow_num_sedi,     snow_num_melt,     snow_num_adj,    &
             snow_num_inst_col, snow_num_sedi_col, snow_num_melt_col, snow_num_adj_col


!  total precip diagnostics

  integer :: a_precip_cld, a_precip_clr, sedi_sfc

!  temperature diagnostics

  integer ::  ST3d, ST_imb
  integer ::  ST2d, ST_imb_col

!  vapor diagnostics

  integer :: SQ3d, qdt_liquid_init, qdt_ice_init, qdt_tiny, qdt_rain_evap,   &
             qdt_cond, qdt_deposition, qdt_eros_l, qdt_eros_i,        &
             qdt_qv_on_qi, qdt_sedi_ice2vapor, qdt_sedi_liquid2vapor,  &
             qdt_super_sat_rm, qdt_destr, qdt_cleanup_liquid,  &
             qdt_cleanup_ice, qdt_snow_sublim, qdt_snow2vapor, SQ_imb
  integer :: SQ2d, q_liquid_init_col, q_ice_init_col, q_tiny_col,  q_rain_evap_col, &
             q_cond_col, q_deposition_col, q_eros_l_col, q_eros_i_col, &
             q_qv_on_qi_col, q_sedi_ice2vapor_col, q_sedi_liquid2vapor_col,&
             q_super_sat_rm_col, q_destr_col, q_cleanup_liquid_col, &
             q_cleanup_ice_col, q_snow_sublim_col, q_snow2vapor_col,  &
             SQ_imb_col

!   miscellaneous diagnostics

  integer :: f_snow_berg, f_snow_berg_col, &
             lsf_strat, lcf_strat, mfls_strat, &
             dcond, vfall

!   CMIP diagnostics (3d)
  type(cmip_diag_id_type) :: cdnc

END TYPE diag_id_type


!#########################################################################

TYPE diag_pt_type

!  cloud area variables

  integer :: aall, aliq, aice, cf_liq_init, cf_ice_init, aauto
  integer :: SA3d, qadt_lsform, qadt_lsdiss, qadt_rhred, qadt_eros,  &
             qadt_fill, qadt_super, qadt_destr, qadt_limits, qadt_ahuco, &
             SA_imb

!  cloud liquid variables

  integer :: SL3d, qldt_cond, qldt_evap, qldt_eros, qldt_berg, qldt_freez,&
             liq_adj, qldt_rime, qldt_accr, qldt_auto, qldt_fill, qldt_tiny,  &
             qldt_destr, qldt_freez2, qldt_sedi, qldt_accrs, qldt_bergs, &
             qldt_HM_splinter, SL_imb

!  cloud ice variables

  integer :: SI3d, qidt_dep, qidt_subl, qidt_fall, qidt_eros, qidt_melt, &
             qidt_melt2, qidt_fill, qidt_tiny, qidt_destr, qidt_qvdep, qidt_auto,  &
             qidt_accr, qidt_accrs, ice_adj, qidt_rain2ice,  SI_imb

!  rain variables
  integer :: qrdt_fill, qrdt_destr, qrdt_tiny

! snow variables
  integer :: qsdt_fill, qsdt_destr, qsdt_tiny

!  rain number variables
  integer :: qnrdt_fill, qnrdt_destr, qnrdt_tiny

! snow number variables
  integer :: qnsdt_fill, qnsdt_destr, qnsdt_tiny

  integer :: SR3d, SR2d, SNR3d, SNR2d, SS3d, SS2d, SNS3d, SNS2d


!  cloud droplet variables

  integer :: droplets_col250, gb_droplets_col, potential_droplets, &
             droplets, droplets_wtd, ql_wt, droplets_col, rvolume
  integer :: SN3d, qndt_cond , qndt_evap, qndt_fill, qndt_tiny, qndt_berg, qndt_rime, &  !h1g, 2014-07-24
             qndt_destr, qndt_super, qndt_freez, qndt_sacws, qndt_sacws_o, &
             qndt_eros, qndt_pra, qndt_auto, qndt_nucclim, qndt_sedi, &
             qndt_melt, qndt_ihom, qndt_size_adj, qndt_fill2,   &
             qndt_contact_frz, qndt_cleanup, qndt_cleanup2, SN_imb

!  cloud ice particle variables

  integer :: nice, nice_col, gb_nice_col, potential_crystals
  integer :: SNi3d, qnidt_fill,  qnidt_tiny, qnidt_nnuccd, qnidt_nsubi,  &
             qnidt_nerosi, qnidt_auto, qnidt_accr, qnidt_nprci, qnidt_nprai,  &
             qnidt_nucclim1, qnidt_nucclim2, qnidt_sedi,             &
             qnidt_melt, qnidt_size_adj, qnidt_fill2,                &
             qnidt_super, qnidt_ihom, qnidt_destr,  qnidt_rain2ice,  &
             qnidt_cleanup, qnidt_cleanup2, qnidt_nsacwi, SNi_imb

!  aerosol diagnostics

  integer :: delta_cf, sulfate, seasalt_sub, seasalt_sup, om,   &
             rhcrit, rhcrit_min, rhiin, rhlin, cfin, imass7,     &
             ni_dust, ni_sulf, ni_bc, ndust1, ndust2, ndust3,  &
             ndust4, ndust5, dust_berg_flag, subgrid_w_variance

!  rain diagnostics
  integer :: rain3d, qrout, rain_clr, rain_cld, a_rain_clr, a_rain_cld, &
             rain_evap, rain_freeze, srfrain_accrs, srfrain_freez,  &
             rain_inst, rain_sedi, &
             srfrain_evap, rain_mass_conv, rain_imb, cld_liq_imb, neg_rain

!  rain number diagnostics
  integer :: rain_num_inst, rain_num_sedi,  rain_num_adj, rain_num2snow,  &
             rain_num_evap, rain_num_freez, rain_num_selfcoll

!  snow diagnostics
  integer :: snow3d, qsout, snow_clr, snow_cld, a_snow_clr, a_snow_cld, &
             snow_inst, snow_sedi,&
             snow_melt, snow_mass_conv, sedi_ice, snow_imb, cld_ice_imb, &
             neg_snow

!  snow number diagnostics
  integer :: snow_num_inst,  snow_num_sedi, snow_num_melt, snow_num_adj



!  total precip diagnostics

  integer :: a_precip_cld, a_precip_clr, sedi_sfc

!  temperature diagnostics

  integer ::  ST3d, ST_imb

!  vapor diagnostics

  integer :: SQ3d, qdt_liquid_init, qdt_ice_init, qdt_tiny, qdt_rain_evap,   &
             qdt_cond, qdt_deposition, qdt_eros_l, qdt_eros_i,        &
             qdt_qv_on_qi, qdt_sedi_ice2vapor, qdt_sedi_liquid2vapor,  &
             qdt_super_sat_rm, qdt_destr, qdt_cleanup_liquid,  &
             qdt_cleanup_ice, qdt_snow_sublim, qdt_snow2vapor, SQ_imb

!   miscellaneous diagnostics

  integer :: f_snow_berg, &
             lsf_strat, lcf_strat, mfls_strat, &
             dcond, vfall

END TYPE diag_pt_type


!#########################################################################

type lscloud_debug_type

!-----------------------------------------------------------------------
!    variables related to debugging options.
!-----------------------------------------------------------------------
! otun               ! file where debug output is written
! debugo
! debugo0 = .false.  ! small output
! debugo4 = .false.  ! when true, nrefuse will be output
! ncall   = 1        ! timestep counter of calls to strat_cloud
! nrefuse
! isamp
! jsamp
! ksamp
!-----------------------------------------------------------------------


 integer :: otun
 integer :: ncall
 integer :: isamp
 integer :: jsamp
 integer :: ksamp
 integer :: nrefuse
 logical :: debugo
 logical :: debugo0
 logical :: debugo4
 integer :: num_strat_pts
 integer, dimension(:,:), pointer :: strat_pts=>NULL()

end type lscloud_debug_type

!#########################################################################

type lscloud_nml_type

  real    :: Dmin
  real    :: cfact
  integer :: super_ice_opt
  logical :: do_ice_nucl_wpdf
  logical :: do_dust_berg
  logical :: do_pdf_clouds
  logical :: pdf_org
  integer :: betaP
  real :: qthalfwidth
  integer :: nsublevels
  integer :: kmap
  integer :: kord

end type lscloud_nml_type


!#########################################################################

type atmos_state_type

!       airdens        air density                     kg air/(m*m*m)
!       qs             saturation specific humidity    kg vapor/kg air
!       dqsdT          T derivative of qs              kg vapor/kg air/K
!       gamma          (L/cp)*dqsdT                    dimensionless
!       deltpg         pressure thickness of grid box  kg air/(m*m)
!                      divided by gravity
!       U_ca           grid box relative humidity      fraction

  real, dimension(:,:,:), pointer ::  &
                                        airdens        =>NULL(), &
                                        tn             =>NULL(), &
                                        qvn            =>NULL(), &
                                        qs             =>NULL(), &
                                        dqsdT          =>NULL(), &
                                        qsi            =>NULL(), &
                                        qsl            =>NULL(), &
                                        rh_crit        =>NULL(), &
                                        rh_crit_min    =>NULL(), &
                                        gamma          =>NULL(), &
                                        esat0          =>NULL(), &
                                        U_ca           =>NULL(), &
                                        delp           =>NULL(), &
                                        U01            =>NULL(), &
                                        pthickness     =>NULL()

end type atmos_state_type


!##########################################################################

type  particles_type

! drop1           number conc                     [1/cm^3]
! drop2           mass concentration              [1/kg]

  real, dimension(:,:,:), pointer ::  &
                                        concen_dust_sub=>NULL(), &
                                        drop1          =>NULL(), &
                                        drop2          =>NULL(), &
                                        crystal1       =>NULL(), &
                                        N3D            =>NULL(), &
                                        N3Di           =>NULL(), &
                                        rbar_dust      =>NULL(), &
                                        ndust          =>NULL(), &
                                        hom            =>NULL()
  real, dimension(:,:,:,:), pointer ::  &
                                        imass1         => NULL(), &
                                        totalmass1     => NULL()

end type particles_type


!########################################################################

type cloud_state_type

!       ql_upd         updated value of ql             kg condensate/
!                                                      kg air
!
!       qi_upd         updated value of qi             kg condensate/
!                                                      kg air
!
!       qa_upd         updated value of qa             fraction
!
!       qa_mean        qa + SA; semi-implicit          fraction
!                      saturated volume fraction
!       ql_mean        ql + positive increment         kg condensate/
!                      of ql; i.e. a sort of           kg air
!                      intermediate ql
!
!       qi_mean        ql + positive increment         kg condensate/
!                      of qi; i.e. a sort of           kg air
!                      intermediate qi
  real, dimension(:,:,:), pointer ::   &
                                        ql_upd         =>NULL(), &
                                        qr_upd         =>NULL(), &
                                        qi_upd         =>NULL(), &
                                        qs_upd         =>NULL(), &
                                        qg_upd         =>NULL(), &
                                        qa_upd         =>NULL(), &
                                        qn_upd         =>NULL(), &
                                        qni_upd        =>NULL(), &
                                        qnr_upd        =>NULL(), &
                                        qns_upd        =>NULL(), &
                                        ql_mean        =>NULL(), &
                                        qr_mean        =>NULL(), &
                                        qi_mean        =>NULL(), &
                                        qs_mean        =>NULL(), &
                                        qg_mean        =>NULL(), &
                                        qa_mean        =>NULL(), &
                                        qn_mean        =>NULL(), &
                                        qni_mean       =>NULL(), &
                                        qnr_mean       =>NULL(), &
                                        qns_mean       =>NULL(), &
                                        ql_in          =>NULL(), &
                                        qr_in          =>NULL(), &
                                        qi_in          =>NULL(), &
                                        qs_in          =>NULL(), &
                                        qg_in          =>NULL(), &
                                        qa_in          =>NULL(), &
                                        qn_in          =>NULL(), &
                                        qni_in         =>NULL(), &
                                        qnr_in         =>NULL(), &
                                        qns_in         =>NULL(), &
                                        SL_out         =>NULL(), &
                                        SI_out         =>NULL(), &
                                        SA_out         =>NULL(), &
                                        SN_out         =>NULL(), &
                                        SNi_out        =>NULL(), &
                                        SR_out         =>NULL(), &
                                        SS_out         =>NULL(), &
                                        SNR_out        =>NULL(), &
                                        SNS_out        =>NULL(), &
                                        SA_0           =>NULL(), &
                                        qa_upd_0       =>NULL(), &
                                        relvarn        =>NULL()

end type cloud_state_type


!#########################################################################

type precip_state_type

  real, dimension(:,:,:), pointer ::   &
                                        lsc_snow       =>NULL(), &
                                        lsc_rain       =>NULL(), &
                                        lsc_snow_size  =>NULL(), &
                                        lsc_rain_size  =>NULL()

  real, dimension(:,:), pointer   ::   &

                                        precip         =>NULL(), &
                                        surfrain       =>NULL(), &
                                        surfsnow       =>NULL()

end type precip_state_type


!#########################################################################

type cloud_processes_type

!       da_ls          change in saturated volume      fraction
!                      fraction due to large-scale
!                      processes
!       D_eros         Sink in ql, qi and qa equation  dimensionless
!                      due to turbulent erosion of
!                      cloud sides
!       dcond_ls       change in condensate due to     kg condensate/
!                      non-convective condensation.    kg air
!                      After phase determination,
!                      this variable refers only to
!                      liquid condensation.
!
!       dcond_ls_ice   change in ice due to            kg condensate/
!                      non-convective condensation.    kg air
!       qvg            equilibrium value of water      kg vapor /
!                      vapor in the clear portion      kg air
!                      of the grid box that PDF
!                      clouds wants
!

  real, dimension(:,:,:), pointer ::   &
                                        da_ls          =>NULL(), &
                                        D_eros         =>NULL(), &
                                        qvg            =>NULL(), &
                                        dcond_ls       =>NULL(), &
                                        dcond_ls_ice   =>NULL(), &
                                        dcond_ls_liquid   =>NULL(), &
                                        dcond_ls_tot   =>NULL(), &
                                        delta_cf       =>NULL(), &
                                        f_snow_berg    =>NULL()

end type cloud_processes_type



!##########################################################################

type lsc_constants_type

  logical                         ::                             &
                                     do_rk_microphys,            &
                                     do_ncar_MG2,                &
                                     tiedtke_macrophysics,       &
                                     dqa_activation,             &
                                     total_activation

end type lsc_constants_type



                             CONTAINS




!########################################################################

subroutine lscloud_types_init

!------------------------------------------------------------------------

      if (module_is_initialized) return

!------------------------------------------------------------------------
!    write version number to output file.
!------------------------------------------------------------------------
      call write_version_number (version, tagname)

!------------------------------------------------------------------------
!    mark this module as initialized.
!------------------------------------------------------------------------
      module_is_initialized = .true.

!------------------------------------------------------------------------

end subroutine lscloud_types_init



!########################################################################



                    end module lscloud_types_mod
