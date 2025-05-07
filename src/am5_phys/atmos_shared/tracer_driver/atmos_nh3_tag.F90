module atmos_nh3_tag_mod

  use  astronomy_mod,          only : astronomy_init, diurnal_solar, universal_time
  use  atmos_cmip_diag_mod,    only : register_cmip_diag_field_2d
  use  atmos_ocean_fluxes_mod, only : aof_set_coupler_flux
  use atmos_tracer_utilities_mod, only :                      &
       get_cmip_param, get_chem_param
  use  constants_mod,          only : grav, rdgas, WTMAIR, WTMH2O, AVOGNO, &
       PI, DEG_TO_RAD, SECONDS_PER_DAY, epsln, WTMN
  use  coupler_types_mod, only : coupler_2d_bc_type, ind_pcair
  use  data_override_mod, only : data_override
  use  diag_manager_mod,  only : send_data,            &
       register_diag_field,  &
       register_static_field, &
       get_base_time
  use  field_manager_mod, only : MODEL_ATMOS, MODEL_LAND, parse
  use  fms_mod, only : write_version_number, &
       mpp_pe,  &
       mpp_root_pe, &
       lowercase,   &
       uppercase, &
       stdlog, &
       check_nml_error, &
       error_mesg, &
       FATAL, &
       WARNING, &
       NOTE
  use  fms2_io_mod, only: file_exists
  use  interpolator_mod, only : interpolate_type,     &
       interpolate_type_eq,  &
       interpolator_init,    &
       obtain_interpolator_time_slices, &
       unset_interpolator_time_flag, &
       interpolator_end,     &
       interpolator,         &
       query_interpolator,   &
       init_clim_diag,       &
       CONSTANT,             &
       INTERP_WEIGHTED_P
  use  mpp_mod, only : input_nml_file
  use  time_manager_mod, only : time_type, &
       get_date, &
       set_date, &
       set_time, &
       days_in_year, &
       real_to_time_type, &
       time_type_to_real, &
       operator(+), operator(-)
  use  tracer_manager_mod, only : get_tracer_index, &
       query_method

  implicit none
  private

  public :: atmos_nh3_tag_driver,          &
       atmos_nh3_tag_init,                 &
       atmos_nh3_tag_time_vary,            &
       atmos_nh3_tag_end,                  &
       atmos_nh3_tag_endts,                &
       atmos_nh3_tag_adjust,               &
       atmos_nh3_tag_gather_data,          &
       atmos_nh3_tag_flux_init,            &
       is_nh3_tag_tracer

  character(len=7), parameter :: mod_name = 'nh3_tag'
  character(len=64), parameter:: sub_name = 'atmos_nh3_tag_init'      

  character(len=80)     :: file_emis    = ' '
  character(len=80)     :: file_emis3d  = ' '
  logical               :: debug        = .false.
  integer               :: nb_tag_max   = 10 
  logical               :: read_nml = .true.

!---- version number -----
  character(len=128) :: version = '$Id$'
  character(len=128) :: tagname = '$Name$'


  namelist /atmos_nh3_tag_nml/ &
       file_emis,file_emis3d,  &
       debug,nb_tag_max


  character(len=7), parameter :: module_name = 'tracers'

  type        :: field_init_type
     character(len=64), pointer :: field_names(:)
     character(len=64)          :: keep_field
     logical                    :: is_ag
     character(len=64)          :: exclude_field(3)
     character(len=64)          :: mask_name
     logical                    :: do_mask_emis        
  end type field_init_type

  type(field_init_type),  allocatable :: emis_field_names(:), &
       emis3d_field_names(:)

  type(interpolate_type), allocatable :: inter_emis(:), inter_emis3d(:)

  integer              :: nnh3,nnh4,noh, nb_tag
  integer, allocatable :: id_loss_mol(:)
  integer, allocatable :: id_emis(:)
  integer, allocatable :: id_emis3d(:)  
  integer, allocatable :: id_emis_tag(:)  

  real, allocatable    :: mask_emis(:,:,:),mask_emis3d(:,:,:)

  logical, allocatable ::  is_nh3_tag(:)
  logical, allocatable :: has_emis(:), has_emis3d(:), &
       land_does_emission(:), &
       diurnal_emis(:), diurnal_emis3d(:)

  real, parameter :: g_to_kg    = 1.e-3,    & !conversion factor (kg/g)
       m2_to_cm2  = 1.e4,     & !conversion factor (cm2/m2)
       twopi      = 2.*PI
  real, parameter :: emis_cons = WTMAIR * g_to_kg * m2_to_cm2 / AVOGNO

  integer, allocatable     :: nnh3_tag(:), nnh4_tag(:)
  integer, allocatable     :: ind_nh3_flux(:)

  logical :: module_is_initialized=.false., first = .true.

contains

  subroutine atmos_nh3_tag_driver( lon, lat, land, ocn_flx_fraction, pwt, r, rm, rdt, &
       Time, phalf, pfull, t, is, ie, js, je, dt,         &
       z_half, z_full, coszen,  &
       half_day,                              &
       Time_next, kbot, emis_ag_masage )

    real, intent(in),    dimension(:,:)            :: lon, lat
    real, intent(in),    dimension(:,:)            :: land    ! land fraction
    real, intent(in),    dimension(:,:)            :: ocn_flx_fraction ! grid box fraction over which DMS flux from ocean occurs
    real, intent(in),    dimension(:,:,:)          :: pwt
    real, intent(in),    dimension(:,:,:,:)        :: r, rm
    real, intent(inout), dimension(:,:,:,:)        :: rdt
    type(time_type), intent(in)                    :: Time, Time_next
    integer, intent(in)                            :: is, ie, js, je
    real, intent(in),    dimension(:,:,:)          :: phalf,pfull,t
    real, intent(in)                               :: dt      ! timestep (s)
    real, intent(in),    dimension(:,:,:)          :: z_half  ! height in meters at half levels
    real, intent(in),    dimension(:,:,:)          :: z_full  ! height in meters at full levels
    real, intent(in),    dimension(:,:)            :: coszen  ! cosine of solar zenith angle
    real, intent(in), dimension(:,:)               :: half_day! half-day length (0 to pi)
    integer, intent(in),  dimension(:,:), optional :: kbot
    real, intent(in),    dimension(:,:), optional  :: emis_ag_masage    ! 

    !local variables
    real, dimension(size(r,1),size(r,2),size(r,3),nb_tag)  :: emis_source, oh_chem_dt
    real, dimension(size(r,1),size(r,2),nb_tag)            :: emisz
    real, dimension(size(r,1),size(r,2),size(r,3))         :: cnh4,cnhx,f_aero
    real, dimension(size(r,1),size(r,2))                   :: emis!, mask_emis
    real, dimension(size(r,1),size(r,2),size(r,3))         :: emis3d

    integer :: i,j,k,n,kb,id,jd,kd
    real :: air_dens,k_nh3_oh

    logical :: used

    if (.not. module_is_initialized)  &
         call error_mesg ('atmos_nh3_tag_driver','atmos_nh3_tag_init must be called first.', FATAL)

    id=size(r,1); jd=size(r,2); kd=size(r,3)

    emis_source(:,:,:,:) = 0.0
    emisz(:,:,:)         = 0.0

    do n = 1, nb_tag

       emis3d(:,:,:)        = 0.0
       emis(:,:)            = 0.0

       !read emissions
       if (has_emis(n)) then
          call read_2D_emis_data( inter_emis(n), emis, Time, Time_next, &
               emis_field_names(n)%field_names, &
               diurnal_emis(n), coszen, half_day, lon, &
               is, js, emis_field_names(n)%keep_field, emis_field_names(n)%exclude_field)

          if (emis_field_names(n)%is_ag .and. present(emis_ag_masage)) then
             emis = emis+emis_ag_masage
          end if

          if (emis_field_names(n)%do_mask_emis) then
             !             call data_override('ATM', trim(emis_field_names(n)%mask_name), mask_emis, time, override=used, is_in=is, ie_in=ie, js_in=js, je_in=je)
             do j=1,jd
                do i=1,id                
                   emis(i,j) = emis(i,j)*mask_emis(i+is-1,j+js-1,n)
                end do
             end do
          end if


          if (id_emis(n) > 0) then
             used = send_data(id_emis(n),emis,Time_next,is_in=is,js_in=js)
          end if
       end if

       !read emissions 3d
       if (has_emis3d(n)) then
          call read_3D_emis_data( inter_emis3d(n), emis3d, Time, Time_next,phalf, &
               emis3d_field_names(n)%field_names, &
               diurnal_emis3d(n), coszen, half_day, lon, &
               is, js )

          if (emis3d_field_names(n)%do_mask_emis) then
             !             call data_override('ATM', trim(emis3d_field_names(n)%mask_name), mask_emis, time, override=used, is_in=is, ie_in=ie, js_in=js, je_in=je)
             do k=1,size(emis3d,3)
                do j=1,jd
                   do i=1,id                
                      emis3d(i,j,k) = emis3d(i,j,k)*mask_emis3d(i+is-1,j+js-1,n)
                   end do
                end do
             end do
          end if
          if (id_emis3d(n) > 0) then
             used = send_data(id_emis3d(n),emis3d,Time_next,is_in=is,js_in=js)
          end if
       end if

       if (present(kbot)) then
          do j=1,jd
             do i=1,id
                kb=kbot(i,j)
                emis_source(i,j,kb,n) = emis(i,j)/pwt(i,j,kb) * emis_cons
             end do
          end do
       else
          emis_source(:,:,kd,n) = emis(:,:)/pwt(:,:,kd) * emis_cons
       end if

       emis_source(:,:,:,n) = emis_source(:,:,:,n) &
            + emis3d(:,:,:)/pwt(:,:,:) * emis_cons       

       emisz(:,:,n) = emis(:,:)
       do k=1, size(emis3d,3)
          emisz(:,:,n) = emisz(:,:,n) + emis3d(:,:,k)
       end do

    end do

    !chemical loss by OH
    do i=1,id
       do j=1,jd
          do k=1,kd
             air_dens = pfull(i,j,k)/(rdgas*wtmair*1e-3*t(i,j,k)) * AVOGNO * 1e-6 ! molec/cm3
             k_nh3_oh = 1.7e-12*exp(-710./t(i,j,k))    
             do n=1,nb_tag
                oh_chem_dt(i,j,k,n) = -k_nh3_oh*r(i,j,k,noh)*r(i,j,k,nnh3_tag(n))*air_dens
             end do
          end do
       end do
    end do

    do n=1,nb_tag
       if(id_loss_mol(n)>0) then
          used = send_data(id_loss_mol(n),-oh_chem_dt(:,:,:,n)*pwt(:,:,:)*1.e3/WTMAIR, &
               Time_next,is_in=is,js_in=js)
       end if
       if(id_emis_tag(n)>0) then
          used = send_data(id_emis_tag(n), emisz(:,:,n)*1.0e04*0.017/AVOGNO, Time_next, is_in=is,js_in=js)
       endif
    end do

    !update nh3 tendency
    do n = 1, nb_tag
       rdt(:,:,:,nnh3_tag(n)) = emis_source(:,:,:,n) + oh_chem_dt(:,:,:,n) + rdt(:,:,:,nnh3_tag(n))
    end do

    !adjustment due to partioning
    cNH4 = rm(:,:,:,nnh4)+rdt(:,:,:,nnh4)*dt
    cNHx = rm(:,:,:,nnh3)+rdt(:,:,:,nnh3)*dt+cNH4
    f_aero =  min(max(cNH4/(max(cNHx,epsln)),0.),1.)

    do i=1,nb_tag
       cNH4 = rm(:,:,:,nnh4_tag(i))+rdt(:,:,:,nnh4_tag(i))*dt
       cNHx = rm(:,:,:,nnh3_tag(i))+rdt(:,:,:,nnh3_tag(i))*dt+cNH4

       rdt(:,:,:,nnh4_tag(i)) = (cNHx*f_aero-rm(:,:,:,nnh4_tag(i)))/dt
       rdt(:,:,:,nnh3_tag(i)) = (cNHx*(1.-f_aero)-rm(:,:,:,nnh3_tag(i)))/dt
    end do


  end subroutine atmos_nh3_tag_driver

  function atmos_nh3_tag_init(nt,axes,Time,lonb_mod,latb_mod,do_nh3_atm_ocean_exchange,do_masage) result(do_nh3_tag)

    integer, intent(in) :: nt
    integer        , intent(in) :: axes(4)
    type(time_type), intent(in) :: Time
    real, intent(in), dimension(:,:) :: lonb_mod
    real, intent(in), dimension(:,:) :: latb_mod
    logical, intent(in)     :: do_nh3_atm_ocean_exchange
    logical, intent(in), optional :: do_masage
    logical                       :: do_nh3_tag

    character*9 :: fld !nh3_tag01
    integer :: ierr, io, logunit
    integer :: i

    character(len=256) :: cmip_name

    if (module_is_initialized) return

    if (debug .and. mpp_pe().eq.mpp_root_pe()) write(*,*) 'starting atmos_nh3_tag'

    !-----------------------------------------------------------------------
    !     ... read namelist
    !-----------------------------------------------------------------------
    call read_nml_file()

    nNH4      = get_tracer_index(MODEL_ATMOS,'nh4')
    nNH3      = get_tracer_index(MODEL_ATMOS,'nh3')
    nOH       = get_tracer_index(MODEL_ATMOS,'oh')

    allocate(is_nh3_tag(nt))
    allocate(nnh3_tag(nb_tag_max))
    allocate(nnh4_tag(nb_tag_max))

    is_nh3_tag(:) = .false.

    nb_tag = 0

    do_nh3_tag = .false.

    if (nnh4.gt.0.and.nnh3.gt.0) then
       ! find out if there are tag nh3 tracers dust tracers
       do i=1,nb_tag_max
          write(fld,'(A7,I2.2)') 'nh3_tag',i
          nNH3_tag(i) = get_tracer_index(MODEL_ATMOS,fld)
          if (debug .and. mpp_root_pe().eq.mpp_pe() .and. nNH3_tag(i).gt.0) write(*,*) fld,nNH3_tag(i)
          write(fld,'(A7,I2.2)') 'nh4_tag',i
          nNH4_tag(i) = get_tracer_index(MODEL_ATMOS,fld)
          if (debug .and. mpp_root_pe().eq.mpp_pe() .and. nNH4_tag(i).gt.0) write(*,*) fld,nNH4_tag(i)
          if (nNH3_tag(i).gt.0.and.nNH4_tag(i).gt.0) then
             do_nh3_tag=.true.
             nb_tag = nb_tag+1
             is_nh3_tag(nnh3_tag(i)) = .true.
          elseif (nNH3_tag(i).gt.0.and.nNH4_tag(i).le.0) then
             call error_mesg ('Tracer_driver','NH3 tag without NH4 tag', FATAL)             
          elseif (nNH4_tag(i).gt.0.and.nNH3_tag(i).le.0) then
             call error_mesg ('Tracer_driver','NH4 tag without NH3 tag', FATAL)             
          end if
       end do

       if (mpp_root_pe().eq.mpp_pe()) write(*,*) 'nb_nh3_tag=',nb_tag

       if (do_nh3_tag) then
          allocate(inter_emis(nb_tag))
          allocate(inter_emis3d(nb_tag))
          allocate(has_emis(nb_tag))
          allocate(has_emis3d(nb_tag))
          allocate(land_does_emission(nb_tag))
          allocate(diurnal_emis(nb_tag))
          allocate(diurnal_emis3d(nb_tag))
          allocate(emis_field_names(nb_tag))
          allocate(emis3d_field_names(nb_tag))
          allocate(id_loss_mol(nb_tag))
          allocate(id_emis(nb_tag))
          allocate(id_emis3d(nb_tag))          
          allocate(id_emis_tag(nb_tag))
          allocate(mask_emis(size(lonb_mod,1)-1,size(latb_mod,2)-1,nb_tag))
          allocate(mask_emis3d(size(lonb_mod,1)-1,size(latb_mod,2)-1,nb_tag))

          mask_emis    = 0.
          mask_emis3d  = 0.

          do i = 1,nb_tag

             write(fld,'(A7,I2.2)') 'nh3_tag',i
             if (mpp_root_pe().eq.mpp_pe()) write(*,*) fld

             call init_emis_data( inter_emis(i), MODEL_ATMOS, 'emissions', nNH3_tag(i),  &
                  lonb_mod, latb_mod, emis_field_names(i), &
                  has_emis(i), diurnal_emis(i), axes, Time, land_does_emis=land_does_emission(i), do_nh3_atm_ocean_exchange=do_nh3_atm_ocean_exchange, do_masage=do_masage )
             call init_emis_data( inter_emis3d(i), MODEL_ATMOS, 'emissions3d', nNH3_tag(i),  &
                  lonb_mod, latb_mod, emis3d_field_names(i), &
                  has_emis3d(i), diurnal_emis3d(i), axes, Time )

             write(fld,'(A7,I2.2)') 'nh3_tag',i
             id_loss_mol(i) = register_diag_field( module_name, fld//'_lossm', axes(1:3), &
                  Time, fld//'_lossm','mole/m2/s')

             if( has_emis(i) ) then
                id_emis(i) = register_diag_field( module_name, fld//'_emis', axes(1:2), &
                     Time, fld//'_emis', 'molec/cm2/s')
             else
                id_emis(i) = 0
             end if
             if( has_emis3d(i) ) then
                id_emis3d(i) = register_diag_field( module_name, fld//'_emis3d', axes(1:3), &
                     Time, fld//'_emis3d', 'molec/cm2/s')
             else
                id_emis3d(i) = 0
             end if

             !register regardless 
             call  get_cmip_param (nNH3_tag(i), cmip_name=cmip_name)

             id_emis_tag(i) = register_diag_field( module_name, 'emi'//fld, axes(1:2), &
                  Time, 'emi'//fld, 'kg/m2/s', &
                  standard_name='emission '//trim(cmip_name))

          end do
       end if
    end if

    module_is_initialized = .true.

  end function atmos_nh3_tag_init

  subroutine atmos_nh3_tag_end

    if (allocated(inter_emis))         deallocate(inter_emis)
    if (allocated(inter_emis3d))       deallocate(inter_emis3d)
    if (allocated(has_emis))           deallocate(has_emis)
    if (allocated(has_emis3d))         deallocate(has_emis3d)
    if (allocated(land_does_emission)) deallocate(land_does_emission)
    if (allocated(diurnal_emis))       deallocate(diurnal_emis)
    if (allocated(diurnal_emis3d))     deallocate(diurnal_emis3d)
    if (allocated(emis_field_names))   deallocate(emis_field_names)
    if (allocated(emis3d_field_names)) deallocate(emis3d_field_names)
    if (allocated(id_emis))            deallocate(id_emis)
    if (allocated(id_emis3d))          deallocate(id_emis3d)
    if (allocated(id_loss_mol))        deallocate(id_loss_mol)
    if (allocated(nnh3_tag))           deallocate(nnh3_tag)
    if (allocated(nnh4_tag))           deallocate(nnh4_tag)
    if (allocated(mask_emis))          deallocate(mask_emis)
    if (allocated(mask_emis3d))        deallocate(mask_emis3d) 
    if (allocated(ind_nh3_flux))       deallocate(ind_nh3_flux)

    module_is_initialized = .false.

  end subroutine atmos_nh3_tag_end

  subroutine atmos_nh3_tag_endts

    integer :: n

    do n=1,size(inter_emis,1)
       if (has_emis(n)) then
          call unset_interpolator_time_flag(inter_emis(n))
       end if
    end do
    do n=1,size(inter_emis3d,1)
       if (has_emis(n)) then
          call unset_interpolator_time_flag(inter_emis3d(n))
       end if
    end do
  end subroutine atmos_nh3_tag_endts

  subroutine atmos_nh3_tag_time_vary(Time)

    type(time_type), intent(in) :: Time
    integer :: n
    logical :: used

    do n=1, size(inter_emis,1)
       if (has_emis(n)) then
          call obtain_interpolator_time_slices (inter_emis(n), Time)
       endif
    end do

    do n=1, size(inter_emis3d,1)
       if (has_emis3d(n)) then
          call obtain_interpolator_time_slices (inter_emis3d(n), Time)
       endif
    end do

    if (first) then
       do n = 1,nb_tag
          if (emis_field_names(n)%do_mask_emis) then
             call data_override('ATM', trim(emis_field_names(n)%mask_name), mask_emis(:,:,n), time, override=used)
          end if
          if (emis3d_field_names(n)%do_mask_emis) then
             call data_override('ATM', trim(emis3d_field_names(n)%mask_name), mask_emis3d(:,:,n), time, override=used)
          end if
       end do
       first=.false.
    end if

  end subroutine atmos_nh3_tag_time_vary

  !copied from tropchem_driver
  subroutine read_2D_emis_data( emis_type, emis, Time, Time_next, &
       field_names, &
       Ldiurnal, coszen, half_day, lon, &
       is, js, keep_field,  skip_field )

    type(interpolate_type),intent(inout) :: emis_type
    real, dimension(:,:),intent(out) :: emis
    type(time_type),intent(in) :: Time, Time_next
    character(len=*),dimension(:), intent(in) :: field_names
    character(len=*), intent(in) :: keep_field
    character(len=*), intent(in) :: skip_field(:)
    logical, intent(in) :: Ldiurnal
    real, dimension(:,:), intent(in) :: coszen, half_day, lon
    integer, intent(in) :: is, js
    integer :: i, j, k, k2
    logical :: used
    real, dimension(size(emis,1),size(emis,2)) :: temp_data
    real :: diurnal_scale_factor, gmt, iso_on, iso_off, dayfrac
    real :: local_angle, factor_tmp
    logical :: found

    emis(:,:) = 0.
    temp_data(:,:) = 0.

    do k = 1,size(field_names)
       temp_data(:,:) = 0.                      
       found=.false.

       do k2 = 1,size(skip_field)
          if (trim(field_names(k)).eq.trim(skip_field(k2))) found=.true.
       end do

       if (.not. found) then
          if (trim(keep_field).eq.'N/A' .or. trim(field_names(k)).eq.trim(keep_field)) then
             call interpolator(emis_type,Time,temp_data,field_names(k),is,js)
          end if
       end if
       emis(:,:) = emis(:,:) + temp_data(:,:)
    end do

    if (Ldiurnal) then
       do j=1,size(emis,2)
          do i=1,size(emis,1)
             if( coszen(i,j) < 0. ) then
                diurnal_scale_factor = 0.
             else
                iso_off = .8 * half_day(i,j)
                iso_on  = -iso_off
                dayfrac = iso_off/PI
                gmt = universal_time(Time)
                local_angle = gmt + lon(i,j) - PI
                if (local_angle >= PI) local_angle = local_angle - twopi
                if (local_angle < -PI) local_angle = local_angle + twopi
                if( local_angle >= iso_off .or. local_angle <= iso_on ) then
                   diurnal_scale_factor = 0.
                else
                   factor_tmp = local_angle - iso_on
                   factor_tmp = factor_tmp / MAX(2.*iso_off,1.e-6)
                   diurnal_scale_factor = 2. / dayfrac * (sin(PI*factor_tmp))**2
                end if
             end if
             emis(i,j) = emis(i,j) * diurnal_scale_factor
          end do
       end do
    end if
  end subroutine read_2D_emis_data


  subroutine read_3D_emis_data( emis_type, emis, Time, Time_next, phalf, &
       field_names, &
       Ldiurnal, coszen, half_day, lon, &
       is, js)

    type(interpolate_type),intent(inout) :: emis_type
    real, dimension(:,:,:),intent(in) :: phalf
    real, dimension(:,:,:),intent(out) :: emis
    type(time_type),intent(in) :: Time, Time_next
    character(len=*),dimension(:), intent(in) :: field_names
    logical, intent(in) :: Ldiurnal
    real, dimension(:,:), intent(in) :: coszen, half_day, lon
    integer, intent(in) :: is, js

    integer :: i, j, k
    logical :: used
    real, dimension(size(emis,1),size(emis,2),size(emis,3)) :: temp_data
    real :: diurnal_scale_factor, gmt, iso_on, iso_off, dayfrac
    real :: local_angle, factor_tmp

    emis(:,:,:) = 0.
    temp_data(:,:,:) = 0.
    do k = 1,size(field_names)
       call interpolator(emis_type,Time,phalf,temp_data,field_names(k),is,js)
       emis(:,:,:) = emis(:,:,:) + temp_data(:,:,:)
    end do
    if (Ldiurnal) then
       do j=1,size(emis,2)
          do i=1,size(emis,1)
             if( coszen(i,j) < 0. ) then
                diurnal_scale_factor = 0.
             else
                iso_off = .8 * half_day(i,j)
                iso_on  = -iso_off
                dayfrac = iso_off/PI
                gmt = universal_time(Time)
                local_angle = gmt + lon(i,j) + PI
                if (local_angle >= PI) local_angle = local_angle - twopi
                if (local_angle < -PI) local_angle = local_angle + twopi
                if( local_angle >= iso_off .or. local_angle <= iso_on ) then
                   diurnal_scale_factor = 0.
                else
                   factor_tmp = local_angle - iso_on
                   factor_tmp = factor_tmp / MAX(2.*iso_off,1.e-6)
                   diurnal_scale_factor = 2. / dayfrac * (sin(PI*factor_tmp))**2
                end if
             end if
             emis(i,j,:) = emis(i,j,:) * diurnal_scale_factor
          end do
       end do
    end if

  end subroutine read_3D_emis_data


  subroutine init_emis_data( emis_type, model, method_type, pos, &
       lonb_mod, latb_mod, field_type, flag, diurnal, &
       axes, Time, land_does_emis, do_nh3_atm_ocean_exchange, do_masage )

    type(interpolate_type),intent(inout) :: emis_type
    integer, intent(in) :: model,pos
    character(len=*),intent(in) :: method_type
    real,intent(in),dimension(:,:) :: lonb_mod,latb_mod
    type(field_init_type),intent(out) :: field_type
    logical, intent(out) :: flag, diurnal
    logical, intent(out), optional :: land_does_emis
    logical, intent(in), optional :: do_nh3_atm_ocean_exchange, do_masage
    integer        , intent(in)  :: axes(4)
    type(time_type), intent(in)  :: Time
    character(len=64) ::file_name

    character(len=64) :: name, control, mask
    integer :: nfields
    integer :: flag_file, flag_diurnal, flag_name, flag_field, flag_mask
    character(len=64) ::  emis_file, control_diurnal, exclude_field, keep_field

    logical   :: used

    flag = .false.
    diurnal = .false.
    control = ''

    if( query_method(trim(method_type),model,pos,name,control) ) then
       if( trim(name(1:4)) == 'file' ) then
          flag = .true.
          flag_file = parse(control, 'file', emis_file)
          flag_diurnal = parse(control, 'diurnal', control_diurnal)
          if(flag_file > 0) then
             file_name = trim(emis_file)
          else
             select case (trim(method_type))
             case ('emissions3d')
                file_name  = trim(file_emis3d)
             case default
                file_name  = trim(file_emis)
             end select
          end if
          flag_field = parse(control, 'keep_field', keep_field)
          if (flag_field>0) then
             field_type%keep_field = keep_field
             if (trim(keep_field) .eq. 'agriculture') then
                field_type%is_ag = .true.
             end if
          else
             field_type%keep_field = 'N/A'
          end if
          flag_field = parse(control, 'exclude_field',exclude_field)
          if (flag_field>0) then
             field_type%exclude_field(1) = exclude_field
          else
             field_type%exclude_field(1) = 'N/A'
          end if
          if (present(do_nh3_atm_ocean_exchange) .and. do_nh3_atm_ocean_exchange) then
             field_type%exclude_field(2) = 'ocean'
          else
             field_type%exclude_field(2) = 'N/A'
          end if
          if (present(do_masage) .and. do_masage) then
             field_type%exclude_field(3) = 'agriculture'
          else
             field_type%exclude_field(3) = 'N/A'
          end if
          diurnal = (flag_diurnal > 0)          
          call interpolator_init( emis_type, trim(file_name), &
               lonb_mod, latb_mod,  &
               data_out_of_bounds=(/CONSTANT/), &
               vert_interp=(/INTERP_WEIGHTED_P/) )
          call query_interpolator(emis_type,nfields=nfields)
          allocate(field_type%field_names(nfields))
          call query_interpolator(emis_type,field_names=field_type%field_names)
       end if
       if ( present(land_does_emis) )  land_does_emis  = (index(lowercase(name),'land:lm3')>0)

       if (mpp_root_pe().eq.mpp_pe()) then
          write(*,*) 'file=              ',file_name
          write(*,*) 'field_names        ',field_type%field_names
          write(*,*) 'exclude field_names',field_type%exclude_field
          write(*,*) 'keep field_names   ',field_type%keep_field
       end if
    else
       field_type%exclude_field='N/A'
       field_type%keep_field='N/A'
    end if

    flag_mask = parse(control, 'mask', mask)
    if (flag_mask.le.0) then
       field_type%do_mask_emis = .false.
       field_type%mask_name    =  'N/A'
    else
       field_type%do_mask_emis = .true.
       field_type%mask_name    =  trim(mask)
       if (mpp_root_pe().eq.mpp_pe()) write(*,*) 'mask_name=',trim(mask)
    end if

  end subroutine init_emis_data

  function is_nh3_tag_tracer(tr) result(ret)
    logical :: ret
    integer, intent(in) :: tr
    ret = (is_nh3_tag(tr))
  end function is_nh3_tag_tracer

  subroutine atmos_nh3_tag_adjust()

  end subroutine atmos_nh3_tag_adjust

  subroutine atmos_nh3_tag_gather_data(gas_fields,tr_bot)
    type(coupler_2d_bc_type), intent(inout) :: gas_fields
    real, dimension(:,:,:), intent(in)      :: tr_bot

    integer :: n

    do n=1,nb_tag
       if (ind_nh3_flux(n) .gt. 0) then
          gas_fields%bc(ind_nh3_flux(n))%field(ind_pcair)%values(:,:) = tr_bot(:,:,nnh3_tag(n))
       end if
    end do

  end subroutine atmos_nh3_tag_gather_data

  subroutine atmos_nh3_tag_flux_init

    integer             :: n
    integer             :: nnh3_tag_for_flux
    character*14        :: fld !nh3_tag01

    if (mpp_root_pe().eq.mpp_pe()) write(*,*) 'setting up nh3_tag_flux (atmos)'

    call read_nml_file()
    allocate(ind_nh3_flux(nb_tag_max))
    ind_nh3_flux = 0

    do n=1,nb_tag_max
       write(fld,'(A7,I2.2)') 'nh3_tag',n
       nnh3_tag_for_flux = get_tracer_index(MODEL_ATMOS,fld)
       if (nnh3_tag_for_flux.gt.0) then
          write(fld,'(A7,I2.2,A5)') 'nh3_tag',n,'_flux'
          ind_nh3_flux(n) = aof_set_coupler_flux(fld,                                    &
               flux_type = 'air_sea_gas_flux_generic', implementation = 'johnson',       &
               atm_tr_index = nnh3_tag_for_flux,                                               &
               mol_wt = WTMN, param = (/ 17.,25. /),                                     &
               caller = trim(mod_name) // '(' // trim(sub_name) // ')')
       end if
    end do
       
  end subroutine atmos_nh3_tag_flux_init

  subroutine read_nml_file()
    integer :: io
    integer :: ierr
    integer :: logunit
    if (read_nml) then
       if (file_exists('input.nml')) then
          read (input_nml_file, nml=atmos_nh3_tag_nml, iostat=io)
          ierr = check_nml_error(io,'atmos_nh3_tag_nml')
       endif
       !--------- write version and namelist to standard log ------------
       call write_version_number(version,tagname)
       logunit = stdlog()
       if (mpp_pe() .eq. mpp_root_pe()) then
          write(logunit,nml=atmos_nh3_tag_nml)
       endif
       read_nml = .false.
    endif
  end subroutine read_nml_file

end module atmos_nh3_tag_mod

