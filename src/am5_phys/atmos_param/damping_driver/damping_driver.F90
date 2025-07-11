module damping_driver_mod

!-----------------------------------------------------------------------
!
!       This module controls five functions:
!
!   (1) rayleigh friction applied to momentum fields at levels
!       1 to kbot (i.e., momentum is damped toward zero).
!
!   (2) mountain gravity wave drag module may be called
!
!   (3) Garner topo_drag module may be called
!
!   (4) Alexander-Dunkerton gravity wave drag may be called
!
!   (5) Beres non-orographic gravity wave scheme may be called, which 
!   can be combined with AD scheme or used by itself.

!-----------------------------------------------------------------------

 use      mg_drag_mod, only:  mg_drag, mg_drag_init, mg_drag_end, &
                              mg_drag_restart
 use      cg_drag_mod, only:  cg_drag_init, cg_drag_calc, cg_drag_end, &
                              cg_drag_time_vary, cg_drag_endts, &
                              cg_drag_restart
 use    topo_drag_mod, only:  topo_drag_init, topo_drag, topo_drag_end, &
                              topo_drag_restart
 use    gw_beres_mod, only: gw_beres_init, gw_beres_ifc, gw_beres_end
 use    tracer_manager_mod, only: get_tracer_index
 use    field_manager_mod, only: MODEL_ATMOS
 use          mpp_mod, only:  input_nml_file
 use  mpp_domains_mod, only:  domain2D
 use          fms_mod, only:  mpp_pe, mpp_root_pe, stdlog, &
                              write_version_number, &
                              error_mesg, &
                              check_nml_error,                   &
                              FATAL
 use diag_manager_mod, only:  register_diag_field,  &
                              register_static_field, send_data
 use atmos_cmip_diag_mod, only: register_cmip_diag_field_3d, &
                                send_cmip_data_3d, cmip_diag_id_type, &
                                query_cmip_diag_id
 use time_manager_mod, only:  time_type
 use    constants_mod, only:  cp_air, grav

 implicit none
 private

 public   damping_driver, damping_driver_init, damping_driver_end
 public   damping_driver_time_vary, damping_driver_endts
 public   damping_driver_restart

!-----------------------------------------------------------------------
!---------------------- namelist ---------------------------------------

   real     :: trayfric = 0.
   integer  :: nlev_rayfric = 1
   logical  :: do_mg_drag = .false.
   logical  :: do_cg_drag = .true.
   logical  :: do_beres_gw = .true. !p1l
   logical  :: beres_diagonly = .true. ! p1l
   logical  :: do_geos_bck = .false. ! p1l
   logical  :: do_topo_drag = .false., use_topo_drag = .true.
   logical  :: do_conserve_energy = .true.

   integer  :: kstart = 0.0         ! rjw 

   namelist /damping_driver_nml/  trayfric, nlev_rayfric,  &
                                  do_cg_drag, do_topo_drag,  &
                                  do_mg_drag, do_conserve_energy, &
                                  use_topo_drag, kstart, &               !stg
                                  do_beres_gw, beres_diagonly, do_geos_bck  !p1l

!
!   trayfric = damping time in seconds for rayleigh damping momentum
!              in the top nlev_rayfric layers (if trayfric < 0 then time
!              in days)
!                 [real, default: trayfric=0.]
!
!   nlev_rayfric = number of levels at the top of the model where
!                  rayleigh friction of momentum is performed, if
!                  trayfric=0. then nlev_rayfric has no effect
!                    [integer, default: nlev_rayfric=1]
!
!-----------------------------------------------------------------------
!----- id numbers for diagnostic fields -----

integer :: id_udt_rdamp,  id_vdt_rdamp,   &
           id_udt_gwd,    id_vdt_gwd,     &
                          id_sgsmtn,      &
           id_udt_cgwd,   id_vdt_cgwd,    &
           id_taubx,      id_tauby,       &
           id_taus

integer :: id_tdt_diss_rdamp,  id_diss_heat_rdamp, &
           id_tdt_diss_gwd,    id_diss_heat_gwd,   &
           id_tdt_diss_topo,   id_diss_heat_topo,  &
           id_mom_flux,  id_diss_heat_cgwd                 !stg

integer :: id_udt_topo,   id_vdt_topo,    &
           id_udtnp_topo, id_vdtnp_topo,  &
           id_taubx_topo, id_tauby_topo,  &
           id_taus_topo

type(cmip_diag_id_type) :: ID_utendogw, ID_utendnogw, &
                           ID_vtendogw, ID_vtendnogw, &
                           ID_tntogw,   ID_tntnogw

!----- missing value for all fields ------

real :: missing_value = -999.

character(len=7) :: mod_name = 'damping'

!-----------------------------------------------------------------------

 logical :: do_rayleigh

 real, parameter ::  daypsec=1./86400.
 logical :: module_is_initialized =.false.

 real :: rfactr

!   note:  
!     rfactr = coeff. for damping momentum at the top level

 character(len=128) :: version = '$Id$'
 character(len=128) :: tagname = '$Name$'

!-----------------------------------------------------------------------

contains

!#######################################################################

 subroutine damping_driver (is, js, lat, Time, delt, area, pfull, phalf, zfull, zhalf, &
                            u, v, t, q, r, u_ref, v_ref, z_pbl, udt, vdt, tdt, &
                            dqcdt_gw, tten_gw, cqa_gw, & !p1l
                             mask, kbot)
 
!-----------------------------------------------------------------------
 integer,         intent(in)                :: is, js
 real, dimension(:,:), intent(in)           :: lat
 type(time_type), intent(in)                :: Time
 real,            intent(in)                :: delt
 real,    intent(in),    dimension(:,:,:)   :: pfull, phalf, &
                                               zfull, zhalf, &
                                               u, v, t, q
 real,    intent(in),    dimension(:,:)     :: u_ref, v_ref !bqx
 real,    intent(in),    dimension(:,:,:,:) :: r
 real,    intent(inout), dimension(:,:,:)   :: udt,vdt,tdt
 real,    intent(in),    dimension(:,:)     :: z_pbl, area 
 class(*),    intent(in),    dimension(:,:,:), optional :: mask
 integer, intent(in),    dimension(:,:),   optional :: kbot
 real,    intent(in),   dimension(:,:,:)  :: dqcdt_gw ! p1l: total condenstation [kg/kg/s]
 real,    intent(in),   dimension(:,:,:)  :: tten_gw ! p1l: grid-mean convective heating rates to be used in gw_beres [K/s]
 real,    intent(in),   dimension(:,:,:)  :: cqa_gw ! p1l: convection updraft area to be used in gw_beres [1]

!-----------------------------------------------------------------------
 real, dimension(size(udt,1),size(udt,2))             :: diag2
 real, dimension(size(udt,1),size(udt,2))             :: taubx, tauby
 real, dimension(size(udt,1),size(udt,2),size(udt,3)) :: taus
 real, dimension(size(udt,1),size(udt,2),size(udt,3)) :: utnd, vtnd, &
                                                         utnd_np, vtnd_np, & !bqx
                                                         ttnd, pmass, &
                                                         p2, uxv            !stg
 real, dimension(size(udt,1),size(udt,2),size(udt,3)) :: utnd_ad, vtnd_ad, & !p1l
                             ttnd_ad, utnd_beres, vtnd_beres, ttnd_beres !p1l
 real, dimension(size(udt,1),size(udt,2),size(udt,3)+1) :: lphalf
 integer :: k
 logical :: used

!-----------------------------------------------------------------------

   if (.not.module_is_initialized) call error_mesg ('damping_driver',  &
                     'damping_driver_init must be called first', FATAL)

!-----------------------------------------------------------------------
! log(phalf) may be needed for diagnostics (pressure level interp)
! pre-compute here

   lphalf = log(phalf)

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!----------------- r a y l e i g h   d a m p i n g ---------------------
!-----------------------------------------------------------------------
   if (do_rayleigh) then

       p2 = pfull * pfull
       call rayleigh (delt, p2, u, v, utnd, vtnd, ttnd)
       udt = udt + utnd
       vdt = vdt + vtnd
       tdt = tdt + ttnd

!----- diagnostics -----

       if ( id_udt_rdamp > 0 ) then
            used = send_data ( id_udt_rdamp, utnd, Time, is, js, 1, &
                               rmask=mask )
       endif

       if ( id_vdt_rdamp > 0 ) then
            used = send_data ( id_vdt_rdamp, vtnd, Time, is, js, 1, &
                               rmask=mask )
       endif

       if ( id_tdt_diss_rdamp > 0 ) then
            used = send_data ( id_tdt_diss_rdamp, ttnd, Time, is, js, 1, &
                               rmask=mask )
       endif

       if ( id_diss_heat_rdamp > 0 ) then
            do k = 1,size(u,3)
              pmass(:,:,k) = phalf(:,:,k+1)-phalf(:,:,k)
            enddo
            diag2 = cp_air/grav * sum(ttnd*pmass,3)
            used = send_data ( id_diss_heat_rdamp, diag2, Time, is, js )
       endif

   endif
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!--------- m t n   g r a v i t y   w a v e   d r a g -------------------
!-----------------------------------------------------------------------
   if (do_mg_drag) then

       call mg_drag (is, js, delt, area, u, v, t, pfull, phalf, zfull, zhalf,  &
                     utnd, vtnd, ttnd, taubx,tauby,taus,        kbot)
       udt = udt + utnd
       vdt = vdt + vtnd
       tdt = tdt + ttnd

!----- diagnostics -----

       if ( id_udt_gwd > 0 ) then
            used = send_data ( id_udt_gwd, utnd, Time, is, js, 1, &
                               rmask=mask )
       endif

       if ( id_vdt_gwd > 0 ) then
            used = send_data ( id_vdt_gwd, vtnd, Time, is, js, 1, &
                               rmask=mask )
       endif

       if ( id_taubx > 0 ) then
            used = send_data ( id_taubx, taubx, Time, is, js )
       endif

       if ( id_tauby > 0 ) then
            used = send_data ( id_tauby, tauby, Time, is, js )
       endif

       if ( id_taus > 0 ) then
           used = send_data ( id_taus, taus, Time, is, js, 1, &
                              rmask=mask )
       endif

       if ( id_tdt_diss_gwd > 0 ) then
            used = send_data ( id_tdt_diss_gwd, ttnd, Time, is, js, 1, &
                               rmask=mask )
       endif

       if ( id_diss_heat_gwd > 0 ) then
            do k = 1,size(u,3)
              pmass(:,:,k) = phalf(:,:,k+1)-phalf(:,:,k)
            enddo
            diag2 = cp_air/grav * sum(ttnd*pmass,3)
            used = send_data ( id_diss_heat_gwd, diag2, Time, is, js )
       endif

       !--- cmip fields (could pre-compute log(phalf) ---
       if (query_cmip_diag_id(ID_utendogw)) then
          used = send_cmip_data_3d (ID_utendogw, utnd, Time, is, js, 1, &
                                    phalf=lphalf) !, rmask=mask )
       endif
       if (query_cmip_diag_id(ID_vtendogw)) then
          used = send_cmip_data_3d (ID_vtendogw, vtnd, Time, is, js, 1, &
                                    phalf=lphalf) !, rmask=mask )
       endif
       if (query_cmip_diag_id(ID_tntogw)) then
          used = send_cmip_data_3d (ID_tntogw, ttnd, Time, is, js, 1, &
                                    phalf=lphalf) !, rmask=mask )
       endif

   endif

!-----------------------------------------------------------------------
!---------topographic   w a v e   d r a g -------------------
!-----------------------------------------------------------------------
   if (do_topo_drag) then

     call topo_drag ( is, js, delt, u, v, t, pfull, phalf, zfull, zhalf, & 
                      lat, u_ref, v_ref, z_pbl,                          & !bqx+
             utnd, vtnd, utnd_np, vtnd_np, ttnd, taubx, tauby, taus, kbot )

     if (use_topo_drag) then  
         if ( kstart > 0 ) then
           do k = kstart, size(u,3)
                utnd(:,:,k)= 0.0*utnd(:,:,k)
                vtnd(:,:,k)= 0.0*vtnd(:,:,k)
                utnd_np(:,:,k)= 0.0*utnd_np(:,:,k)
                vtnd_np(:,:,k)= 0.0*vtnd_np(:,:,k)
           enddo
         endif 

       udt = udt + utnd
       vdt = vdt + vtnd
       tdt = tdt + ttnd  
     endif

!----- diagnostics -----

     if ( id_udt_topo > 0 ) then
        used = send_data ( id_udt_topo, utnd, Time, is, js, 1, &
                           rmask=mask )
     endif

     if ( id_vdt_topo > 0 ) then
        used = send_data ( id_vdt_topo, vtnd, Time, is, js, 1, &
                           rmask=mask )
     endif
!bqx+
     if ( id_udtnp_topo > 0 ) then
        used = send_data ( id_udtnp_topo, utnd_np, Time, is, js, 1, &
                           rmask=mask )
     endif

     if ( id_vdtnp_topo > 0 ) then
        used = send_data ( id_vdtnp_topo, vtnd_np, Time, is, js, 1, &
                           rmask=mask )
     endif
!bqx

     if ( id_taubx_topo > 0 ) then
       used = send_data ( id_taubx_topo, taubx, Time, is, js )
     endif

     if ( id_tauby_topo > 0 ) then
        used = send_data ( id_tauby_topo, tauby, Time, is, js )
     endif

     if ( id_taus_topo > 0 ) then
        used = send_data ( id_taus_topo, taus, Time, is, js, 1, &
                           rmask=mask )
     endif

     if ( id_tdt_diss_topo > 0 ) then
        used = send_data ( id_tdt_diss_topo, ttnd, Time, is, js, 1, &
                               rmask=mask )
     endif

     if ( id_diss_heat_topo > 0 ) then
          do k = 1,size(u,3)
             pmass(:,:,k) = phalf(:,:,k+1)-phalf(:,:,k)
          enddo
          diag2 = cp_air/grav * sum(ttnd*pmass,3)
          used = send_data ( id_diss_heat_topo, diag2, Time, is, js )
     endif

     !--- cmip fields (could pre-compute log(phalf) ---
     if (query_cmip_diag_id(ID_utendogw)) then
        used = send_cmip_data_3d (ID_utendogw, utnd, Time, is, js, 1, &
                                  phalf=lphalf) !, rmask=mask )
     endif
     if (query_cmip_diag_id(ID_vtendogw)) then
        used = send_cmip_data_3d (ID_vtendogw, vtnd, Time, is, js, 1, &
                                  phalf=lphalf) !, rmask=mask )
     endif
     if (query_cmip_diag_id(ID_tntogw)) then
        used = send_cmip_data_3d (ID_tntogw, ttnd, Time, is, js, 1, &
                                  phalf=lphalf) !, rmask=mask )
     endif

 endif


!rjw         Save vertically-integrated momemtum flux 
     uxv = u*v   !stg
!!!     vxt = v*t   !rjw 

     if ( id_mom_flux > 0 ) then      !stg
          do k = 1,size(u,3)
            pmass(:,:,k) = phalf(:,:,k+1)-phalf(:,:,k)
          enddo
          diag2 = sum(uxv*pmass,3)/grav
          used = send_data ( id_mom_flux, diag2, Time, is, js )
     endif

!-----------------------------------------------------------------------
!---------non-orographic  wave drag: AD99--------------------------
!-----------------------------------------------------------------------
!   Alexander-Dunkerton gravity wave drag

   if (do_cg_drag) then

     call cg_drag_calc (is, js, lat, pfull, zfull, t, u, v, Time,    &
                        delt, utnd_ad, vtnd_ad, ttnd_ad )
   else
   utnd_ad=0.0
   vtnd_ad=0.0
   ttnd_ad=0.0
   endif

!-----------------------------------------------------------------------
!---------non-orographic  wave drag: Beres--------------------------
!-----------------------------------------------------------------------
   if (do_beres_gw) then
      if (do_geos_bck) then
          call gw_beres_ifc(is, js, lat, Time, delt, pfull, phalf, zfull, u, v, t, tten_gw, &
        cqa_gw, utnd_beres, vtnd_beres, ttnd_beres, dqcdt_in = dqcdt_gw)
      else
        call gw_beres_ifc(is, js, lat, Time, delt, pfull, phalf, zfull, u, v, t, tten_gw, &
        cqa_gw, utnd_beres, vtnd_beres, ttnd_beres)
      endif
   else
       utnd_beres=0.0
       vtnd_beres=0.0
       ttnd_beres=0.0
   endif

! Combining the two CG scheme
 if (beres_diagonly) then
     utnd=utnd_ad
     vtnd=vtnd_ad
     ttnd=ttnd_ad
 else
     utnd=utnd_ad+utnd_beres
     vtnd=vtnd_ad +vtnd_beres
     ttnd=ttnd_ad + ttnd_beres
 endif

     udt =  udt + utnd
     vdt =  vdt + vtnd
     tdt =  tdt + ttnd

!----- diagnostics -----

       if ( id_udt_cgwd > 0 ) then
          used = send_data ( id_udt_cgwd, utnd, Time, is, js, 1, &
                            rmask=mask )
       endif
       if ( id_vdt_cgwd > 0 ) then
          used = send_data ( id_vdt_cgwd, vtnd, Time, is, js, 1, &
                            rmask=mask )
       endif

       if ( id_diss_heat_cgwd > 0 ) then
           do k = 1,size(u,3)
             pmass(:,:,k) = phalf(:,:,k+1)-phalf(:,:,k)
           enddo
           diag2 = cp_air/grav * sum(ttnd*pmass,3)
           used = send_data ( id_diss_heat_cgwd, diag2, Time, is, js )
       endif

       !--- cmip fields (could pre-compute log(phalf) ---
       if (query_cmip_diag_id(ID_utendnogw)) then
          used = send_cmip_data_3d (ID_utendnogw, utnd, Time, is, js, 1, &
                                    phalf=lphalf) !, rmask=mask )
       endif
       if (query_cmip_diag_id(ID_vtendnogw)) then
          used = send_cmip_data_3d (ID_vtendnogw, vtnd, Time, is, js, 1, &
                                    phalf=lphalf) !, rmask=mask )
       endif
       if (query_cmip_diag_id(ID_tntnogw)) then
          used = send_cmip_data_3d (ID_tntnogw, ttnd, Time, is, js, 1, &
                                    phalf=lphalf) !, rmask=mask )
       endif
!-----------------------------------------------------------------------

 end subroutine damping_driver

!#######################################################################

 subroutine damping_driver_init ( domain, lonb, latb, pref, axes, Time, sgsmtn)

 type(domain2D), target, intent(in)  :: domain !< Atmosphere domain
 real,            intent(in) :: lonb(:,:), latb(:,:), pref(:)
 integer,         intent(in) :: axes(4)
 type(time_type), intent(in) :: Time
 real, dimension(:,:), intent(out) :: sgsmtn
!-----------------------------------------------------------------------
!     lonb  = longitude in radians of the grid box corners
!     latb  = latitude  in radians of the grid box corners
!     axes  = axis indices, (/x,y,pf,ph/)
!               (returned from diag axis manager)
!     Time  = current time (time_type)
!     sgsmtn = subgrid scale topography variance
!-----------------------------------------------------------------------
 integer :: unit, ierr, io, logunit
 logical :: used

!-----------------------------------------------------------------------
!----------------- namelist (read & write) -----------------------------

   read (input_nml_file, nml=damping_driver_nml, iostat=io)
   ierr = check_nml_error(io,"damping_driver_nml")

   call write_version_number(version, tagname)
   logunit = stdlog()
   if(mpp_pe() == mpp_root_pe() ) then
        write (logunit,nml=damping_driver_nml)
   endif

!-----------------------------------------------------------------------
!--------- both mg_drag and topo_drag can not be active ------

   if (do_mg_drag .and. do_topo_drag) call error_mesg ('damping_driver',  &
                 'do_mg_drag and do_topo_drag can not both be true', FATAL)

!-----------------------------------------------------------------------
!--------- rayleigh friction ----------

   do_rayleigh=.false.

   if (abs(trayfric) > 0.0001 .and. nlev_rayfric > 0) then
      if (trayfric > 0.0) then
         rfactr=(1./trayfric)
      else
         rfactr=(1./abs(trayfric))*daypsec
      endif
         do_rayleigh=.true.
   else
         rfactr=0.0
   endif

!-----------------------------------------------------------------------
!----- mountain gravity wave drag -----

   if (do_mg_drag) call mg_drag_init (domain, lonb, latb, sgsmtn)

!--------------------------------------------------------------------
!----- Alexander-Dunkerton gravity wave drag -----
 
   if (do_cg_drag)  then
     call cg_drag_init (domain, lonb, latb, pref, Time=Time, axes=axes)
   endif

!--------------------------------------------------------------------
!----- Beres gravity wave drag -----

   if (do_beres_gw)  then
     call gw_beres_init(lonb, latb, pref, Time, axes)
   endif

!-----------------------------------------------------------------------
!----- initialize diagnostic fields -----

if (do_rayleigh) then

   id_udt_rdamp = &
   register_diag_field ( mod_name, 'udt_rdamp', axes(1:3), Time,       &
                       'u wind tendency for Rayleigh damping', 'm/s2', &
                        missing_value=missing_value               )

   id_vdt_rdamp = &
   register_diag_field ( mod_name, 'vdt_rdamp', axes(1:3), Time,       &
                       'v wind tendency for Rayleigh damping', 'm/s2', &
                        missing_value=missing_value               )

   id_tdt_diss_rdamp = &
   register_diag_field ( mod_name, 'tdt_diss_rdamp', axes(1:3), Time,  &
                      'Dissipative heating from Rayleigh damping',&
                             'deg_k/s', missing_value=missing_value   )
       
   id_diss_heat_rdamp = &
   register_diag_field ( mod_name, 'diss_heat_rdamp', axes(1:2), Time,   &
                'Integrated dissipative heating from Rayleigh damping',&
                  'W/m2' )
endif

if (do_mg_drag) then

 ! register and send static field
   id_sgsmtn = &
   register_static_field ( mod_name, 'sgsmtn', axes(1:2), &
               'sub-grid scale topography for gravity wave drag', 'm')
   if (id_sgsmtn > 0) used = send_data (id_sgsmtn, sgsmtn, Time)

 ! register non-static field
   id_udt_gwd = &
   register_diag_field ( mod_name, 'udt_gwd', axes(1:3), Time,        &
                     'u wind tendency for gravity wave drag', 'm/s2', &
                        missing_value=missing_value               )

   id_vdt_gwd = &
   register_diag_field ( mod_name, 'vdt_gwd', axes(1:3), Time,        &
                     'v wind tendency for gravity wave drag', 'm/s2', &
                        missing_value=missing_value               )

   id_taubx = &
   register_diag_field ( mod_name, 'taubx', axes(1:2), Time,        &
                         'x base flux for grav wave drag', 'N/m2',  &
                         missing_value=missing_value )

   id_tauby = &
   register_diag_field ( mod_name, 'tauby', axes(1:2), Time,        &
                         'y base flux for grav wave drag', 'N/m2',  &
                         missing_value=missing_value  )

   id_taus = &
   register_diag_field ( mod_name, 'taus', axes(1:3), Time,             &
                       'saturation flux for gravity wave drag', 'N/m2', &
                       missing_value=missing_value    )

   id_tdt_diss_gwd = &
   register_diag_field ( mod_name, 'tdt_diss_gwd', axes(1:3), Time,     &
                          'Dissipative heating from gravity wave drag', &
                              'deg_K/s', missing_value=missing_value   )
       
   id_diss_heat_gwd = &
   register_diag_field ( mod_name, 'diss_heat_gwd', axes(1:2), Time,      &
                'Integrated dissipative heating from gravity wave drag',  &
                                 'W/m2' )
endif

   if (do_cg_drag) then

    id_udt_cgwd = &
    register_diag_field ( mod_name, 'udt_cgwd', axes(1:3), Time,     &
                 'u wind tendency for cg gravity wave drag', 'm/s2', &
                      missing_value=missing_value               )


    id_vdt_cgwd = &
    register_diag_field ( mod_name, 'vdt_cgwd', axes(1:3), Time,     &
                 'v wind tendency for cg gravity wave drag', 'm/s2', &
                      missing_value=missing_value               )

   id_diss_heat_cgwd = &
   register_diag_field ( mod_name, 'diss_heat_cgwd', axes(1:2), Time,      &
                'Integrated dissipative heating from convective gravity wave drag',  &
                                 'W/m2' )

     !----- cmip diagnostics for non-orographic drag -----
     ID_utendnogw = register_cmip_diag_field_3d ( mod_name, 'utendnogw', Time, &
                           'U-tendency Nonorographic Gravity Wave Drag', 'm s-2', &
              standard_name='tendency_of_eastward_wind_due_to_nonorographic_gravity_wave_drag')

     ID_vtendnogw = register_cmip_diag_field_3d ( mod_name, 'vtendnogw', Time, &
                           'V-tendency Nonorographic Gravity Wave Drag', 'm s-2', &
              standard_name='tendency_of_northward_wind_due_to_nonorographic_gravity_wave_drag')

     ID_tntnogw = register_cmip_diag_field_3d ( mod_name, 'tntnogw', Time, &
                           'Temperature Tendency due to Nonorographic Gravity Wave Dissipation', 'K s-1', &
              standard_name='temperature_tendency_due_to_dissipation_nonorographic_gravity_wave_drag')
   endif

!-----------------------------------------------------------------------
!----- topo wave drag -----



  if (do_topo_drag) then
          call topo_drag_init (domain, lonb, latb)
          sgsmtn(:,:) = -99999.
  endif

  if (do_topo_drag) then

   id_udt_topo = &
   register_diag_field ( mod_name, 'udt_topo', axes(1:3), Time,        &
                         'u wind tendency for topo wave drag', 'm/s2', &
                         missing_value=missing_value )

   id_vdt_topo = &
   register_diag_field ( mod_name, 'vdt_topo', axes(1:3), Time,        &
                         'v wind tendency for topo wave drag', 'm/s2', &
                         missing_value=missing_value )

   id_udtnp_topo = &
   register_diag_field ( mod_name, 'udtnp_topo', axes(1:3), Time,        &
                         'u wind tendency for non-propagating topo wave drag', 'm/s2', &
                         missing_value=missing_value )

   id_vdtnp_topo = &
   register_diag_field ( mod_name, 'vdtnp_topo', axes(1:3), Time,        &
                         'v wind tendency for non-propagating topo wave drag', 'm/s2', &
                         missing_value=missing_value )

   id_taubx_topo = &
   register_diag_field ( mod_name, 'taubx_topo', axes(1:2), Time,      &
                         'x base flux for topo wave drag', 'kg/m/s2',  &
                         missing_value=missing_value )

   id_tauby_topo = &
   register_diag_field ( mod_name, 'tauby_topo', axes(1:2), Time,      &
                         'y base flux for topo wave drag', 'kg/m/s2',  &
                         missing_value=missing_value )

   id_taus_topo = &
   register_diag_field ( mod_name, 'taus_topo', axes(1:3), Time,          &
                         'saturation flux for topo wave drag', 'kg/m/s2', &
                         missing_value=missing_value )

   id_tdt_diss_topo = &
   register_diag_field ( mod_name, 'tdt_diss_topo', axes(1:3), Time,   &
                         'Dissipative heating from topo wave drag',    &
                         'deg_k/s', missing_value=missing_value )
       
   id_diss_heat_topo = &
   register_diag_field ( mod_name, 'diss_heat_topo', axes(1:2), Time,          &
                         'Integrated dissipative heating from topo wave drag', &
                         'W/m2' )


 endif

! rjw     Save vertically-integrated momentum flux 

   id_mom_flux = &                                                            !stg
   register_diag_field ( mod_name, 'mom_flux', axes(1:2), Time,                         &
                'Integrated meridional flux of zonal momentum from topo wave drag',     &
                  'J/m2' )

!-----------------------------------------------------------------------
!----- cmip diagnostics for orographic drag -----

  if (do_mg_drag .or. do_topo_drag) then

     ID_utendogw = register_cmip_diag_field_3d ( mod_name, 'utendogw', Time, &
                           'U-tendency Orographic Gravity Wave Drag', 'm s-2', &
              standard_name='tendency_of_eastward_wind_due_to_orographic_gravity_wave_drag')

     ID_vtendogw = register_cmip_diag_field_3d ( mod_name, 'vtendogw', Time, &
                           'V-tendency Orographic Gravity Wave Drag', 'm s-2', &
              standard_name='tendency_of_northward_wind_due_to_orographic_gravity_wave_drag')

     ID_tntogw = register_cmip_diag_field_3d ( mod_name, 'tntogw', Time, &
                           'Temperature Tendency due to Orographic Gravity Wave Dissipation', 'K s-1', &
              standard_name='temperature_tendency_due_to_dissipation_orographic_gravity_wave_drag')
   endif

!-----------------------------------------------------------------------

   module_is_initialized =.true.

!******************** end of initialization ****************************
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

 end subroutine damping_driver_init


!#####################################################################

subroutine damping_driver_time_vary (delt)

real, intent(in) :: delt


       call cg_drag_time_vary (delt)

end subroutine damping_driver_time_vary



!#####################################################################

subroutine damping_driver_endts


     call cg_drag_endts

end subroutine damping_driver_endts



!######################################################################
!#######################################################################

 subroutine damping_driver_end

     if (do_mg_drag)   call mg_drag_end
     if (do_cg_drag)   call cg_drag_end
     if (do_topo_drag) call topo_drag_end
     if (do_beres_gw) call gw_beres_end

     module_is_initialized =.false.


 end subroutine damping_driver_end

!#######################################################################
!
subroutine damping_driver_restart(timestamp)
  character(len=*), intent(in), optional :: timestamp
!-----------------------------------------------------------------------
! write out restart file.
! Arguments: 
!   timestamp (optional, intent(in)) : A character string that represents the model time, 
!                                      used for writing restart. timestamp will append to
!                                      the any restart file name as a prefix. 
!-----------------------------------------------------------------------
     if (do_mg_drag)   call mg_drag_restart(timestamp)
     if (do_cg_drag)   call cg_drag_restart(timestamp)
     if (do_topo_drag) call topo_drag_restart(timestamp)

end subroutine damping_driver_restart

!#######################################################################

 subroutine rayleigh (dt, p2, u, v, udt, vdt, tdt)

  real,    intent(in)                      :: dt
  real,    intent(in),  dimension(:,:,:)   :: p2, u, v
  real,    intent(out), dimension(:,:,:)   :: udt, vdt, tdt

  real, dimension(size(u,1),size(u,2)) :: fact
  integer :: k
!-----------------------------------------------------------------------
!--------------rayleigh damping of momentum (to zero)-------------------

   do k = 1, nlev_rayfric
     fact(:,:) = rfactr*(1.+(p2(:,:,1)-p2(:,:,k))/(p2(:,:,1)+p2(:,:,k)))
     udt(:,:,k) = -u(:,:,k)*fact(:,:)
     vdt(:,:,k) = -v(:,:,k)*fact(:,:)
   enddo

   do k = nlev_rayfric+1, size(u,3)
     udt(:,:,k) = 0.0
     vdt(:,:,k) = 0.0
   enddo

!  total energy conservation
!  compute temperature change loss due to ke dissipation

   if (do_conserve_energy) then
       do k = 1, nlev_rayfric
          tdt(:,:,k) = -((u(:,:,k)+.5*dt*udt(:,:,k))*udt(:,:,k) +  &
                         (v(:,:,k)+.5*dt*vdt(:,:,k))*vdt(:,:,k)) / cp_air
       enddo
       do k = nlev_rayfric+1, size(u,3)
          tdt(:,:,k) = 0.0
       enddo
   else
       tdt = 0.0
   endif

!-----------------------------------------------------------------------

 end subroutine rayleigh

!#######################################################################

end module damping_driver_mod
