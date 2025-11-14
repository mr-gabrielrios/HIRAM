
module ice_model_mod

use   ice_albedo_mod, only:  ice_albedo_init, ice_albedo
use ocean_albedo_mod, only:  compute_ocean_albedo_new
use  ocean_rough_mod, only:  compute_ocean_roughness, fixed_ocean_roughness

use  amip_interp_mod, only: amip_interp_type, amip_interp_new, &
                            get_amip_ice, get_amip_sst
use time_manager_mod, only: time_type, get_time, operator(+)
use diag_manager_mod, only: diag_axis_init, register_diag_field, send_data
use    constants_mod, only: HLV, HLF, TFREEZE, pi

use          fms_mod, only: file_exist, open_namelist_file, open_restart_file, &
                            close_file, mpp_pe, mpp_root_pe, mpp_npes,         &
                            write_version_number, stdlog, error_mesg, FATAL,   &
                            check_nml_error, read_data, write_data, NOTE,      &
                            set_domain, nullify_domain
use       fms_io_mod, only: get_restart_io_mode

use       mpp_io_mod, only: mpp_open, mpp_close, mpp_get_info, mpp_read,   &
                            MPP_RDONLY, MPP_NETCDF, MPP_MULTI, MPP_SINGLE, &
                            fieldtype, mpp_get_atts, mpp_get_fields

use  mpp_domains_mod, only: domain2d,  mpp_get_layout,  &
                            mpp_get_global_domain, mpp_get_compute_domain
use mpp_mod, only: mpp_min, mpp_max

!miz
use atmos_model_mod,  only: atmos_data_type
use diag_manager_mod, only: register_diag_field, send_data
use data_override_mod,only: data_override !for overriding flux_oh
!miz

implicit none
private
real :: cmin,cmax

public :: ice_data_type, atmos_ice_boundary_type,               &
          ice_model_init, ice_model_end, update_ice_model_fast, &
          update_ice_model_slow
!miz	  
character(len=9) :: mod_name = 'ice_model'	  
integer :: id_flux_oh
logical :: used
!miz
!----------------------------------------------------------------
!--------- NAMELIST INTERFACE ---------
!
!  use_climo_ice           = use monthly climatological amip ice mask
!  use_annual_ice          = use annual climatology amip ice mask
!  temp_ice_freeze         = temperature at which sea ice melts
!  heat_capacity_ocean     = heat capacity for ocean in simple mixed layer model
!
real    :: diff                     = 2.092
real    :: thickness_min            = 0.10
real    :: specified_ice_thickness  = 2.0
real    :: heat_capacity_ocean      = 1.e07
real    :: temp_ice_freeze          = -1.66
real    :: roughness_ice            = 1.e-4
logical :: mixed_layer_ocean        = .false.
logical :: use_climo_ice            = .false.
logical :: use_annual_ice           = .false.
logical :: use_climo_sst            = .false.
logical :: use_annual_sst           = .false.
character(len=64) :: ice_method = 'prognostic' ! none, uniform, or prognostic
character(len=64) :: sst_method = 'specified'  ! specified, uniform, or mixed_layer
                                               ! Additional sst specifications: 'aqua_planet_#'
                                               !   aqua_planet_1 = Control profile
                                               !   aqua_planet_2 = Peaked
                                               !   aqua_planet_3 = Flat
                                               !   aqua_planet_4 = Qobs
                                               !   aqua_planet_5 = Control-5N
                                               !   aqua_planet_6 = 1KEQ
                                               !   aqua_planet_7 = 3KEQ
                                               !   aqua_planet_8 = 3KW1
real              :: temp_ice = 270.      ! used when ice_method = 'uniform'
real              :: temp_sst = 280.      ! used when sst_method = 'uniform'
!miz/tmm
real    :: ohflux_hemi_asym_width=50.
real    :: ohflux_hemi_asym_amp=0.
real    :: trop_ohflux_amp = 0.0
real    :: trop_ohflux_width = 16.0  ! width of tropqflux region in degrees
real    :: walker_ohflux_amp = 0.0
real    :: rose_amp_oh = 0.0
real    :: itcz_ohflux_amp = 0.0 ! can be negative
real    :: itcz_ohflux_width = 6.0 ! half width in degrees
real    :: itcz_ohflux_phi0 = 0.0 ! lat in degrees of center Qflux
logical :: use_wide_itcz_ohflux  = .false. ! spread compensating perturbation Q-flux over wider region
logical :: use_simple_ice_albedo = .false. ! assume different albedo below temp_ice
real    :: simple_ice_albedo = 0.5
!miz/tmm
character(len=64) :: interp_method  = "bilinear" ! conservative or bilinear
logical :: do_netcdf_restart = .true.

real    :: e_folding_width = 5.
real    :: sine_exponent = 4.
real    :: latitude_of_maximum_SST = 10.
real    :: max_SST = 29 ! maximum temperature in Celsius
real    :: min_SST = 0 ! minimum temperature in Celsius
real    :: sine_wavenumber = 1.25

! GR: slab-ocean parameters (Q-flux = ocean energy flux divergence)
logical :: use_TC_ohflux  = .false. ! 
real    :: TC_ohflux_a = 1
real    :: TC_ohflux_b = 1
real    :: TC_ohflux_c = 1
real    :: q_flux_width = 5 ! e-folding width of positive Q-flux in degrees (keep between 1 and 10)
real    :: q_flux_center_latitude = 0 ! meridional offset of Q-flux profile from Equator
real    :: q_flux_0 = 40 ! magnitude of Q-flux

namelist /ice_model_nml/ diff, thickness_min, specified_ice_thickness,        &
                         heat_capacity_ocean, temp_ice_freeze, roughness_ice, &
                         ice_method, use_climo_ice, use_annual_ice, temp_ice, &
                         sst_method, use_climo_sst, use_annual_sst, temp_sst, &
                         interp_method, do_netcdf_restart, ohflux_hemi_asym_width,&
                         ohflux_hemi_asym_amp, trop_ohflux_amp, trop_ohflux_width, &
                         walker_ohflux_amp, rose_amp_oh, itcz_ohflux_amp, &
                         itcz_ohflux_width, itcz_ohflux_phi0, use_wide_itcz_ohflux, &
                         use_simple_ice_albedo, simple_ice_albedo,  & !miz/tmm
                         e_folding_width, sine_exponent, latitude_of_maximum_SST, &
                         max_SST, min_SST, sine_wavenumber, q_flux_width, q_flux_center_latitude, & 
                         q_flux_0, use_TC_ohflux, TC_ohflux_a, TC_ohflux_b, TC_ohflux_c

!----------------------------------------------------------------

type ice_data_type
  type(domain2d),pointer                :: Domain

   real,    pointer, dimension(:,:)     :: glon_bnd =>NULL(), &
                                           glat_bnd =>NULL(), &
                                           lon_bnd =>NULL() , &
                                           lat_bnd =>NULL()

   real,    pointer, dimension(:,:)     :: glon =>NULL(), &
                                           glat =>NULL(), &
                                           lon =>NULL(), &
                                           lat =>NULL()

   logical, pointer, dimension(:,:)     :: gmask =>NULL(), &
                                           mask =>NULL(), &
                                           ice_mask =>NULL()

   real,    pointer, dimension(:,:)   :: t_surf =>NULL(), &
                                         albedo =>NULL(), &
                                         albedo_vis_dir =>NULL(), &
                                         albedo_nir_dir =>NULL(), &
                                         albedo_vis_dif =>NULL(), &
                                         albedo_nir_dif =>NULL(), &
                                         rough_mom =>NULL(),&
                                         rough_heat =>NULL(), &
                                         rough_moist =>NULL(), &
                                         thickness =>NULL()

   type (time_type)                   :: Time_Init, Time,  &
                                         Time_step_fast,   &
                                         Time_step_slow
end type ice_data_type

!----------------------------------------------------------------

type :: atmos_ice_boundary_type
  real, dimension(:,:), pointer :: u_star =>NULL(), &
                                   t_flux =>NULL(), &
                                   q_flux =>NULL(), &
                                   lw_flux =>NULL(), &
                                   sw_flux =>NULL(), &
                                   lprec =>NULL(), &
                                   fprec =>NULL()
  real, dimension(:,:), pointer :: dhdt =>NULL(), &
                                   dedt =>NULL(), &
                                   drdt =>NULL(), &
                                   coszen =>NULL(), &
                                   data =>NULL()
  integer :: xtype
end type atmos_ice_boundary_type

!----------------------------------------------------------------

integer :: is, ie, js, je
type(amip_interp_type), save :: Amip_ice, Amip_sst
logical :: module_is_initialized = .false.
character(len=64) :: fname = 'INPUT/ice_model.res.nc'

character(len=128) :: version = '$Id: ice_model.F90,v 1.1.1.1 2011/07/11 21:09:17 miz Exp $'
character(len=128) :: tagname = '$Name:  $'

real, parameter :: LATENT = HLV + HLF

contains

!######################################################################

 subroutine update_ice_model_fast( Atmos_boundary, Ice )
 type(atmos_ice_boundary_type), intent(in)    :: Atmos_boundary
 type (ice_data_type),          intent(inout) :: Ice

 real, dimension(is:ie,js:je) :: ts_new, gamma, flux_i, t_dt_surf, &
                 flux_t_new, flux_q_new, flux_lw_new, flux_sw_new, &
                 flux_u_new, flux_v_new, lprec_new,   fprec_new,   &
                 deriv, flux_oh, flux_oh_itcz, &
                 q_flux_rios, q_flux_pos, q_flux_neg

 logical, dimension(is:ie,js:je,2) :: mask_ice
 real,    dimension(is:ie,js:je,2) :: thickness_ice, t_surf_ice, albedo
 real,    dimension(is:ie,js:je)   :: albedo_vis_dir, albedo_nir_dir, &
                                      albedo_vis_dif, albedo_nir_dif
 integer :: dt
!miz
 real    :: lat, lon, pi, freq, rad_trop_width
 real    :: q_flux_width_radians, q_flux_center_latitude_radians
 real    :: meridional_offset, root_min, root_max, period
 integer :: i, j
 pi = 4. * atan(1.0)
!miz

!-----------------------------------------------------------------------
!
!   updates ice model on the atmospheric (fast) time step 
!   averages input quantities to be seen by the ocean
!
!    flux_u  = zonal wind stress
!    flux_v  = meridional wind stress
!    flux_sw = net shortwave radiation (down-up) 
!    flux_sw_vis = net visible shortwave radiation (down-up)
!    flux_sw_dir = net direct shortwave radiation (down-up)
!    flux_sw_dif = net diffuse shortwave radiation (down-up)
!    flux_sw_vis_dir = net visible direct shortwave radiation (down-up)
!    flux_sw_vis_dif = net visible diffuse shortwave radiation (down-up)
!    flux_lw = net longwave radiation (down-up) 
!    flux_t  = sensible heat flux
!    flux_q  = specific humidity flux
!    lprec   = liquid precipitiation rate (kg/m2/s)
!    fprec   = frozen precipitiation rate (kg/m2/s)
!    coszen  = cosine of the zenith angle
!
!-----------------------------------------------------------------------
!----- implicit update of ice surface temperature -----

if (trim(ice_method) == 'prognostic') then

    where (Ice%ice_mask)
      gamma = diff / max(Ice%thickness,thickness_min)
      flux_i = gamma * (TFREEZE + temp_ice_freeze - Ice%t_surf)

      t_dt_surf = (flux_i + Atmos_boundary%lw_flux  + Atmos_boundary%sw_flux - &
                            Atmos_boundary%t_flux - Atmos_boundary%q_flux*LATENT)  &
                            / (Atmos_boundary%dhdt + Atmos_boundary%dedt*LATENT +  &
                               Atmos_boundary%drdt + gamma)
      ts_new = Ice%t_surf + t_dt_surf
    elsewhere
      ts_new = 0.
    endwhere

  ! update sea ice temperature
  ! new temperature cannot be greater than freezing temperature
    where (Ice%ice_mask .and. ts_new > TFREEZE )
      Ice%t_surf  = TFREEZE
    endwhere

    where (Ice%ice_mask .and. ts_new <= TFREEZE)
      Ice%t_surf  = ts_new
    endwhere

endif


if (trim(ice_method) == 'uniform') then
    where (Ice%ice_mask)
      !t_dt_surf = temp_ice - Ice%t_surf
       Ice%t_surf = temp_ice
    endwhere
endif

!-----------------------------------------------------------------------
!----- simple mixed layer ocean ------

   if (trim(sst_method) == 'mixed_layer') then

       call get_time ( Ice%Time_step_fast, dt )

       ! hemispherically asymmetric Qflux as in Kang et al 2008
       do i = is, ie
          do j = js, je
             lat  = Ice%lat(i,j)*180/pi; 
       
             if ( -90. < lat .and. lat < -90. + ohflux_hemi_asym_width ) then
                flux_oh(i,j) = ohflux_hemi_asym_amp * sin( (lat + 90.-ohflux_hemi_asym_width)*pi / ohflux_hemi_asym_width )
             elseif ( 90.-ohflux_hemi_asym_width < lat .and. lat < 90. ) then
                flux_oh(i,j) = ohflux_hemi_asym_amp * sin( (lat - 90.+ohflux_hemi_asym_width)*pi / ohflux_hemi_asym_width )
             else
                flux_oh(i,j) = 0.0
             endif
          enddo
       enddo
         
       where (Ice%mask .and. .not. Ice%ice_mask)
          flux_i = ( Atmos_boundary%lw_flux + Atmos_boundary%sw_flux + flux_oh - &
                     Atmos_boundary%t_flux - Atmos_boundary%q_flux*HLV - &
                     Atmos_boundary%fprec*HLF ) * real(dt)/heat_capacity_ocean
          deriv = -( Atmos_boundary%dhdt + Atmos_boundary%dedt*HLV + &
                     Atmos_boundary%drdt) * real(dt)/heat_capacity_ocean 
          t_dt_surf = flux_i/(1.0 -deriv) 
          ts_new = Ice%t_surf + t_dt_surf
       endwhere

     ! update sea surface temperature
     ! note: temperatures allowed below freezing for conservation
       where (Ice%mask .and. .not.Ice%ice_mask)
         Ice%t_surf  = ts_new
       endwhere

   endif

   ! Rose and Ferreira 2012 Q-flux formulation
   if (trim(sst_method) == 'mixed_layer_rose') then

       call get_time ( Ice%Time_step_fast, dt )

       ! lowest order Qflux as in Rose and Ferreira 2012
       ! equivalent to second legendre polynomial
       flux_oh = - rose_amp_oh * ( cos(Ice%lat)**2 - 2*sin(Ice%lat)**2 )

       where (Ice%mask .and. .not. Ice%ice_mask)
          flux_i = ( Atmos_boundary%lw_flux + Atmos_boundary%sw_flux + flux_oh - &
                     Atmos_boundary%t_flux - Atmos_boundary%q_flux*HLV - &
                     Atmos_boundary%fprec*HLF ) * real(dt)/heat_capacity_ocean
          deriv = -( Atmos_boundary%dhdt + Atmos_boundary%dedt*HLV + &
                     Atmos_boundary%drdt) * real(dt)/heat_capacity_ocean 
          t_dt_surf = flux_i/(1.0 -deriv) 
          ts_new = Ice%t_surf + t_dt_surf
       endwhere

     ! update sea surface temperature
     ! note: temperatures allowed below freezing for conservation
       where (Ice%mask .and. .not.Ice%ice_mask)
         Ice%t_surf  = ts_new
       endwhere

   endif
   
   ! Rios Q-flux formulation
   if (trim(sst_method) == 'mixed_layer_rios') then

       call get_time ( Ice%Time_step_fast, dt )

       ! Convert latitudes into radians
       q_flux_width_radians = q_flux_width * pi / 180 ! e-folding width of the Q-flux width
       q_flux_center_latitude_radians = q_flux_center_latitude * pi / 180 ! meridional offset of the Q-flux maximum from the Equator

       ! Q-flux is a sum of a strong positive component and weaker negative component
       q_flux_pos = exp(-(Ice%lat - q_flux_center_latitude_radians)**2 / (2 * q_flux_width_radians **2)) 
       q_flux_neg = -0.5 * exp(-(Ice%lat - 2 * q_flux_center_latitude_radians)**2 / (2 * q_flux_width_radians **2))
       q_flux_rios = 2 * q_flux_0 * (q_flux_pos + q_flux_neg)

       where (Ice%mask .and. .not. Ice%ice_mask)
          flux_i = ( Atmos_boundary%lw_flux + Atmos_boundary%sw_flux + q_flux_rios - &
                     Atmos_boundary%t_flux - Atmos_boundary%q_flux*HLV - &
                     Atmos_boundary%fprec*HLF ) * real(dt)/heat_capacity_ocean
          deriv = -( Atmos_boundary%dhdt + Atmos_boundary%dedt*HLV + &
                     Atmos_boundary%drdt) * real(dt)/heat_capacity_ocean 
          t_dt_surf = flux_i/(1.0 -deriv) 
          ts_new = Ice%t_surf + t_dt_surf
       endwhere
     
       ! Update sea surface temperature
       ! Note: temperatures allowed below freezing for conservation
       where (Ice%mask .and. .not.Ice%ice_mask)
         Ice%t_surf  = ts_new
       endwhere

   endif
   
   if (trim(sst_method) == 'mixed_layer_rios-sine') then

       call get_time ( Ice%Time_step_fast, dt )

       ! Constants for SST profile shift
       meridional_offset = q_flux_center_latitude * pi / 180. ! convert from degrees to radians
       root_min = (-(pi / 2) / sine_wavenumber + meridional_offset) ! calculate minimum root (where SST = 0 degC)
       root_max = ((pi / 2) / sine_wavenumber + meridional_offset) ! calculate maximum root for filtering (where SST = 0 degC)

       q_flux_rios = merge(q_flux_0 * (1 - sin(sine_wavenumber * Ice%lat - sine_wavenumber * meridional_offset)**sine_exponent) + min_SST, &
                           min_SST, &
                           ((Ice%lat .ge. root_min) .and. (Ice%lat .le. root_max))) 

       where (Ice%mask .and. .not. Ice%ice_mask)
          flux_i = ( Atmos_boundary%lw_flux + Atmos_boundary%sw_flux + q_flux_rios - &
                     Atmos_boundary%t_flux - Atmos_boundary%q_flux*HLV - &
                     Atmos_boundary%fprec*HLF ) * real(dt)/heat_capacity_ocean
          deriv = -( Atmos_boundary%dhdt + Atmos_boundary%dedt*HLV + &
                     Atmos_boundary%drdt) * real(dt)/heat_capacity_ocean 
          t_dt_surf = flux_i/(1.0 -deriv) 
          ts_new = Ice%t_surf + t_dt_surf
       endwhere
     
       ! Update sea surface temperature
       ! Note: temperatures allowed below freezing for conservation
       where (Ice%mask .and. .not.Ice%ice_mask)
         Ice%t_surf  = ts_new
       endwhere

   endif
   
   if (trim(sst_method) == 'mixed_layer_rios-implied_ocean_qflux') then

       call get_time ( Ice%Time_step_fast, dt )

       ! Constants for SST profile shift
       q_flux_rios = -0.1343 * exp(-(Ice%lat - 0.6048)**2 / (2 * (0.3924**2)))
        
       where (Ice%mask .and. .not. Ice%ice_mask)
          flux_i = ( Atmos_boundary%lw_flux + Atmos_boundary%sw_flux + q_flux_rios - &
                     Atmos_boundary%t_flux - Atmos_boundary%q_flux*HLV - &
                     Atmos_boundary%fprec*HLF ) * real(dt)/heat_capacity_ocean
          deriv = -( Atmos_boundary%dhdt + Atmos_boundary%dedt*HLV + &
                     Atmos_boundary%drdt) * real(dt)/heat_capacity_ocean 
          t_dt_surf = flux_i/(1.0 -deriv) 
          ts_new = Ice%t_surf + t_dt_surf
       endwhere
     
       ! Update sea surface temperature
       ! Note: temperatures allowed below freezing for conservation

       ! Constants for SST profile shift
       min_SST = 0. ! minimum SST is chosen to be 7K lower than the peak, per Hsieh et al. (2020) 
       max_SST = 29. ! taken from HadISST zonal mean SST maximum, see hadisst_sst.clim.1986-2005.nc
       period = 1.25
       meridional_offset = 15 * pi / 180. ! convert from degrees to radians
       root_min = (-(pi / 2) / period + meridional_offset) ! calculate minimum root (where SST = 0 degC)
       root_max = ((pi / 2) / period + meridional_offset) ! calculate maximum root for filtering (where SST = 0 degC)
       
       do j = js, je
           do i = is, ie
               Ice%t_surf(i, j) = merge((max_SST - min_SST) * (1 - sin(period * Ice%lat(i, j) - meridional_offset)**2) + min_SST + TFREEZE, & 
                                         min_SST + TFREEZE, &
                                        ((Ice%lat(i, j) .ge. root_min) .and. (Ice%lat(i, j) .le. root_max))) 
           enddo
       enddo
       where (Ice%mask .and. .not.Ice%ice_mask)
         Ice%t_surf  = ts_new
       endwhere

   endif

   ! Merlis and Schneider 2011 Q-flux formulation
   if (trim(sst_method) == 'mixed_layer_merlis') then

       call get_time ( Ice%Time_step_fast, dt )

       ! tropical Qflux as in Bordoni 2007
       ! with Merlis et al 2012 correction
       rad_trop_width = trop_ohflux_width*pi/180.
       flux_oh =  -trop_ohflux_amp*(1-2.*Ice%lat**2/rad_trop_width**2) * &
                  exp(- ((Ice%lat)**2/(rad_trop_width)**2)) / cos(Ice%lat)

       ! parameters for zonally asymmetric perturbation hard coded
       flux_oh = flux_oh - &
          walker_ohflux_amp*exp( - ( (Ice%lon - 270.0*pi/180.)**2/((30.0*pi/180.)**2) + &
           (Ice%lat)**2/((7.0*pi/180.)**2) ) ) + &
          walker_ohflux_amp*exp( - ( (Ice%lon - 90.0*pi/180.)**2/((30.0*pi/180.)**2) + &
           (Ice%lat)**2/((7.0*pi/180)**2) ) )

       where (Ice%mask .and. .not. Ice%ice_mask)
          flux_i = ( Atmos_boundary%lw_flux + Atmos_boundary%sw_flux + flux_oh - &
                     Atmos_boundary%t_flux - Atmos_boundary%q_flux*HLV - &
                     Atmos_boundary%fprec*HLF ) * real(dt)/heat_capacity_ocean
          deriv = -( Atmos_boundary%dhdt + Atmos_boundary%dedt*HLV + &
                     Atmos_boundary%drdt) * real(dt)/heat_capacity_ocean 
          t_dt_surf = flux_i/(1.0 -deriv) 
          ts_new = Ice%t_surf + t_dt_surf
       endwhere

     ! update sea surface temperature
     ! note: temperatures allowed below freezing for conservation
       where (Ice%mask .and. .not.Ice%ice_mask)
         Ice%t_surf  = ts_new
       endwhere

   endif

   !miz
   if (trim(sst_method) == 'mixed_layer_data') then

       call get_time ( Ice%Time_step_fast, dt )
        
       call data_override('ICE', 'flux_oh', flux_oh, Ice%Time)   

       where (Ice%mask .and. .not. Ice%ice_mask)
          flux_i = ( Atmos_boundary%lw_flux + Atmos_boundary%sw_flux + flux_oh - &
                     Atmos_boundary%t_flux - Atmos_boundary%q_flux*HLV - &
                     Atmos_boundary%fprec*HLF ) * real(dt)/heat_capacity_ocean
          deriv = -( Atmos_boundary%dhdt + Atmos_boundary%dedt*HLV + &
                     Atmos_boundary%drdt) * real(dt)/heat_capacity_ocean 
          t_dt_surf = flux_i/(1.0 -deriv) 
          ts_new = Ice%t_surf + t_dt_surf
       endwhere

     ! update sea surface temperature
     ! note: temperatures allowed below freezing for conservation
       where (Ice%mask .and. .not.Ice%ice_mask)
         Ice%t_surf  = ts_new
       endwhere

   endif

   if (trim(sst_method) == 'uniform') then
       where (Ice%mask .and. .not.Ice%ice_mask)
          Ice%t_surf = temp_sst
       endwhere
   endif

!-----------------------------------------------------------------------
!------ update ocean/ice surface parameters --------

   !---- over all ocean points (regardless of ice) ----

    call compute_ocean_roughness ( Ice%mask,        &
                        Atmos_boundary%u_star(:,:), &
                              Ice%rough_mom  (:,:), &
                              Ice%rough_heat (:,:), &
                              Ice%rough_moist(:,:) )

    call compute_ocean_albedo_new ( Ice%mask, Atmos_boundary%coszen(:,:),  &
               albedo_vis_dir, albedo_vis_dif, albedo_nir_dir, albedo_nir_dif, Ice%lat )

   !---- over only ice points -----

    where (Ice%ice_mask)
      Ice%rough_mom   = roughness_ice
      Ice%rough_heat  = roughness_ice
      Ice%rough_moist = roughness_ice
    endwhere

    mask_ice     (:,:,2) = Ice%ice_mask
    thickness_ice(:,:,2) = Ice%thickness
    t_surf_ice   (:,:,2) = Ice%t_surf
    albedo       (:,:,2) = Ice%albedo
    albedo       (:,:,1) = Ice%albedo
    call ice_albedo (mask_ice, thickness_ice, t_surf_ice, albedo)

    where (Ice%ice_mask)
       Ice%albedo = albedo(:,:,2)
       Ice%albedo_vis_dir = albedo(:,:,2)
       Ice%albedo_nir_dir = albedo(:,:,2)
       Ice%albedo_vis_dif = albedo(:,:,2)
       Ice%albedo_nir_dif = albedo(:,:,2)
    elsewhere
       Ice%albedo = albedo(:,:,1)
       Ice%albedo_vis_dir = albedo_vis_dir
       Ice%albedo_nir_dir = albedo_vis_dir
       Ice%albedo_vis_dif = albedo_vis_dir
       Ice%albedo_nir_dif = albedo_vis_dir
    endwhere

    ! Tim M. mod Aug. 8, 2013       
    if ( use_simple_ice_albedo ) then
       where (Ice%t_surf < temp_ice .and. .not. Ice%ice_mask)
    	     Ice%albedo = simple_ice_albedo
	     Ice%albedo_vis_dir = simple_ice_albedo
	     Ice%albedo_nir_dir = simple_ice_albedo
	     Ice%albedo_vis_dif = simple_ice_albedo
	     Ice%albedo_nir_dif = simple_ice_albedo
       endwhere	
    endif	

!miz
   if ( id_flux_oh > 0 ) used = send_data ( id_flux_oh, flux_oh, Ice%Time )
!miz
!-----------------------------------------------------------------------
!--------- advance time -----------------

 Ice%Time = Ice%Time + Ice%Time_step_fast

!-----------------------------------------------------------------------

 end subroutine update_ice_model_fast

!######################################################################

 subroutine update_ice_model_slow( Atmos_boundary, Ice )
 type(atmos_ice_boundary_type), intent(in)    :: Atmos_boundary
 type(ice_data_type),           intent(inout) :: Ice

 !---- get the specified sea-ice fraction -----

      if (trim(ice_method) == 'prognostic' .or. trim(ice_method) == 'uniform') then

          call prognostic_ice ( Ice )

      endif

 !---- get the specified ocean temperature -----

      if (trim(sst_method) == 'specified') then

          call prognostic_sst ( Ice )

      endif


 end subroutine update_ice_model_slow

!######################################################################

 subroutine prognostic_ice ( Ice )
 type(ice_data_type),           intent(inout) :: Ice

 real, dimension(is:ie, js:je) :: ice_frac

 !---- get the specified sea-ice fraction -----

    call get_amip_ice (Ice%Time, Amip_ice, ice_frac )

  ! determine which grid boxes have ice coverage
    where ( Ice%mask(:,:) .and. ice_frac > 0.5 )
        Ice%thickness(:,:) = specified_ice_thickness
        Ice%ice_mask (:,:) = .true.
        Ice%t_surf   (:,:) = MIN( Ice%t_surf(:,:), TFREEZE )
    elsewhere
        Ice%thickness(:,:) = 0.0
        Ice%ice_mask (:,:) = .false.
        Ice%t_surf   (:,:) = MAX( Ice%t_surf(:,:), TFREEZE )
    endwhere

 end subroutine prognostic_ice

!######################################################################

 subroutine prognostic_sst ( Ice )
 type(ice_data_type),           intent(inout) :: Ice

 real, dimension(is:ie, js:je) :: sea_temp

 !---- get the specified ocean temperature -----

    call get_amip_sst (Ice%Time, Amip_sst, sea_temp )

  ! determine which grid boxes have open ocean coverage ----
    where ( Ice%mask(:,:) .and. .not.Ice%ice_mask(:,:))
        Ice%t_surf   (:,:) = sea_temp
    endwhere

 end subroutine prognostic_sst

!######################################################################

! subroutine update_ice_model_slow_up ( Ocean_boundary, Ice )
! type(ocean_ice_boundary_type), intent(in)    :: Ocean_boundary
! type(ice_data_type),           intent(inout) :: Ice

! set temp at open ocean points

! end subroutine update_ice_model_slow_up

!######################################################################

 subroutine ice_model_init ( Ice, Atm, Time_Init, Time, &   !miz
                             Time_step_fast, Time_step_slow, &
                             glon_bnd, glat_bnd, Atmos_domain )
 type(ice_data_type), intent(inout) :: Ice
 type(time_type)    , intent(in)    :: Time_Init, Time, &
                                       Time_step_fast, Time_step_slow
!miz
 type(atmos_data_type), intent(in)  :: Atm
!miz
 real               , intent(in)    :: glon_bnd(:,:), glat_bnd(:,:)
 type(domain2d), intent(in), target :: Atmos_domain

 real :: lon0, lond, latd, amp, t_control, dellon
 ! GR: custom SST profile constants
 real :: period, meridional_offset, root_min, root_max

 integer :: isg, ieg, jsg, jeg
 integer :: unit, ierr, io, i, j
 integer :: ndim, nvar, natt, ntime, nlon, nlat, mlon, mlat, layout(2)
 type(fieldtype), allocatable :: Fields(:)
 logical :: need_ic

 if (module_is_initialized) then
     return
 endif

 if ( file_exist( 'input.nml' ) ) then
    unit = open_namelist_file ( )
    ierr = 1
    do while ( ierr /= 0 )       
       read ( unit,  nml = ice_model_nml, iostat = io, end = 10 )
       ierr = check_nml_error ( io, 'ice_model_nml' )
    enddo
 10 continue
    call close_file (unit)       
 endif

 call get_restart_io_mode(do_netcdf_restart)

 call write_version_number (version, tagname)
 if ( mpp_pe() == mpp_root_pe() ) then
    write (stdlog(), nml=ice_model_nml)
 endif

!---- error checks ----

  if ( trim(ice_method) /= 'none'    .and. &
       trim(ice_method) /= 'uniform' .and. &
       trim(ice_method) /= 'prognostic' ) call error_mesg &
     ('ice_model_init', 'namelist variable ice_method has invalid value', FATAL)

  if ( trim(sst_method) /= 'specified' .and. &
       trim(sst_method) /= 'uniform'   .and. &
       trim(sst_method) /= 'aqua_planet_1'  .and. &
       trim(sst_method) /= 'aqua_planet_2'  .and. &
       trim(sst_method) /= 'aqua_planet_3'  .and. &
       trim(sst_method) /= 'aqua_planet_4'  .and. &
       trim(sst_method) /= 'aqua_planet_5'  .and. &
       trim(sst_method) /= 'aqua_planet_6'  .and. &
       trim(sst_method) /= 'aqua_planet_7'  .and. &
       trim(sst_method) /= 'aqua_planet_8'  .and. &
       trim(sst_method) /= 'aqua_planet_9'  .and. &
       trim(sst_method) /= 'aqua_planet_9b'  .and. &
       trim(sst_method) /= 'aqua_planet_11'  .and. &
       trim(sst_method) /= 'aqua_planet_12'  .and. &
       trim(sst_method) /= 'aqua_planet_10N' .and. &
       trim(sst_method) /= 'aqua_planet_15N' .and. &
       trim(sst_method) /= 'aqua_planet_20N' .and. &
       trim(sst_method) /= 'aqua_planet_sine'  .and. &
       trim(sst_method) /= 'aqua_planet_gauss' .and. &
       trim(sst_method) /= 'mixed_layer'    .and. &               !miz
       trim(sst_method) /= 'mixed_layer_data' .and. &             !miz
       trim(sst_method) /= 'mixed_layer_merlis' .and. &             !tmm
       trim(sst_method) /= 'mixed_layer_rios' .and. &             
       trim(sst_method) /= 'mixed_layer_rios-sine' .and. &             
       trim(sst_method) /= 'mixed_layer_rios-implied_ocean_qflux' .and. &             
       trim(sst_method) /= 'mixed_layer_rose' ) call error_mesg & !miz
     ('ice_model_init', 'namelist variable sst_method has invalid value', FATAL)

!----------------------------------------------------------
!--- open the grid_spec file ---

! call mpp_open ( unit, 'INPUT/grid_spec.nc', MPP_RDONLY, MPP_NETCDF, &
!                 threading=MPP_MULTI, fileset = MPP_SINGLE )
! call mpp_get_info (unit, ndim, nvar, natt, ntime)
! allocate (Fields(nvar))
! call mpp_get_fields (unit, Fields)

! call get_grid_size ( Fields, nlon, nlat )
  nlon = size(glon_bnd,1)-1
  nlat = size(glon_bnd,2)-1

!----------------------------------------------------------
!--- set up domain type ---

! if (present(Atmos_domain)) then
!     call mpp_get_layout (Atmos_domain, layout)
! else
!     call mpp_define_layout  ( (/1,nlon,1,nlat/), mpp_npes(), layout )
! endif

! call mpp_define_domains ( (/1,nlon,1,nlat/), layout, Ice%Domain, &
!                           name='ice grid')

  Ice%Domain => Atmos_domain

!----------------------------------------------------------
! get global domain indices 
! this assumes that domain2d type has been assigned

  call mpp_get_global_domain ( Ice%Domain, isg, ieg, jsg, jeg )

  allocate ( Ice%glon_bnd (isg:ieg+1,jsg:jeg+1), &
             Ice%glat_bnd (isg:ieg+1,jsg:jeg+1), &
             Ice%glon     (isg:ieg,jsg:jeg), &
             Ice%glat     (isg:ieg,jsg:jeg), &
             Ice%gmask    (isg:ieg,jsg:jeg)  )

!----------------------------------------------------------
!--- read grid info ----

   Ice%glon_bnd = glon_bnd
   Ice%glat_bnd = glat_bnd
  !do j = js, je+1
  !do i = is, ie+1
  !   Ice%glon_bnd(i,j) = glon_bnd(i-isg+1,j-jsg+1)
  !   Ice%glat_bnd(i,j) = glat_bnd(i-isg+1,j-jsg+1)
  !enddo
  !enddo

  ! read the land mask from a file (land=1)
   if (file_exist('INPUT/land_mask.nc')) then
      call read_data ('INPUT/land_mask.nc', 'land_mask', Ice%glon, no_domain=.true.) !, Ice%Domain)
      where (Ice%glon > 0.50)
         Ice%gmask = .false.
      elsewhere
         Ice%gmask = .true.
      endwhere
   else
      Ice%gmask = .true.  ! aqua-planet
   endif

   ! mid-point of grid box
  !if (is_latlon(glon_bnd,latb_out)) then
  !   do j = js, je
  !   do i = is, ie
  !      Ice%glon(i,j) = (Ice%glon_bnd(i,j  )+Ice%glon_bnd(i+1,j  )+ &
  !                       Ice%glon_bnd(i,j+1)+Ice%glon_bnd(i+1,j+1))*0.25
  !      Ice%glat(i,j) = (Ice%glat_bnd(i,j  )+Ice%glat_bnd(i+1,j  )+ &
  !                       Ice%glat_bnd(i,j+1)+Ice%glat_bnd(i+1,j+1))*0.25
  !   enddo
  !   enddo
  !else
      call get_cell_center (Ice%glon_bnd, Ice%glat_bnd, Ice%glon, Ice%glat)
  !endif

!  call read_grid_data ( unit, Fields, Ice%glon_bnd, Ice%glat_bnd, &
!                        Ice%glon, Ice%glat, Ice%gmask )
!  call mpp_close(unit)
!  deallocate (Fields)

!----------------------------------------------------------
! get compute domain indices 

  call mpp_get_compute_domain ( Ice%Domain, is, ie, js, je )

  allocate ( Ice%lon_bnd        (is:ie+1,js:je+1), &
             Ice%lat_bnd        (is:ie+1,js:je+1), &
             Ice%lon            (is:ie, js:je)   , &
             Ice%lat            (is:ie, js:je)   , &
             Ice%ice_mask       (is:ie, js:je)   , &
             Ice%t_surf         (is:ie, js:je)   , &
             Ice%albedo         (is:ie, js:je)   , &
             Ice%albedo_vis_dir (is:ie, js:je)   , &
             Ice%albedo_nir_dir (is:ie, js:je)   , &
             Ice%albedo_vis_dif (is:ie, js:je)   , &
             Ice%albedo_nir_dif (is:ie, js:je)   , &
             Ice%rough_mom      (is:ie, js:je)   , &
             Ice%rough_heat     (is:ie, js:je)   , &
             Ice%rough_moist    (is:ie, js:je)   , &
             Ice%thickness      (is:ie, js:je)   , &
             Ice%mask           (is:ie, js:je)     )

    Ice%lon_bnd    = Ice%glon_bnd(is:ie+1,js:je+1)
    Ice%lat_bnd    = Ice%glat_bnd(is:ie+1,js:je+1)
    Ice%lon        = Ice%glon(is:ie, js:je)
    Ice%lat        = Ice%glat(is:ie, js:je)
    Ice%mask       = Ice%gmask(is:ie, js:je)
    Ice%Time           = Time
    Ice%Time_init      = Time_init
    Ice%Time_step_fast = Time_step_fast
    Ice%Time_step_slow = Time_step_slow

!----------------------------------------------------------
!----------- read restart -------------

need_ic = .false.

if (file_exist('INPUT/ice_model.res.nc', domain=Ice%Domain )) then
   if (mpp_pe() == mpp_root_pe()) call error_mesg ('ice_model_mod', &
            'Reading NetCDF formatted restart file: INPUT/ice_model.res.nc', NOTE)
   call read_data(fname, 'mlon', mlon, Ice%Domain)
   call read_data(fname, 'mlat', mlat, Ice%Domain)
   if (mlon /= nlon .or. mlat /= nlat )  &
        call error_mesg ('ice_model_init',           &       
                        'incorrect resolution on restart', FATAL)
   call read_data ( fname, 't_surf',         Ice%t_surf,         Ice%Domain )
   call read_data ( fname, 'thickness',      Ice%thickness,      Ice%Domain )
   call read_data ( fname, 'albedo',         Ice%albedo,         Ice%Domain )
   call read_data ( fname, 'albedo_vis_dir', Ice%albedo_vis_dir, Ice%Domain )
   call read_data ( fname, 'albedo_nir_dir', Ice%albedo_nir_dir, Ice%Domain )
   call read_data ( fname, 'albedo_vis_dif', Ice%albedo_vis_dif, Ice%Domain )
   call read_data ( fname, 'albedo_nir_dif', Ice%albedo_nir_dif, Ice%Domain )
   call read_data ( fname, 'rough_mom',      Ice%rough_mom,      Ice%Domain )
   call read_data ( fname, 'rough_heat',     Ice%rough_heat,     Ice%Domain )
   call read_data ( fname, 'rough_moist',    Ice%rough_moist,    Ice%Domain)
else
   if (file_exist('INPUT/ice_model.res')) then
      if (mpp_pe() == mpp_root_pe()) call error_mesg ('ice_model_mod', &
            'Reading native formatted restart file.', NOTE)
      call set_domain (Ice%Domain)
      unit = open_restart_file ('INPUT/ice_model.res', 'read')
      read  (unit) mlon, mlat

    ! restart resolution must be consistent with grid spec
      if (mlon /= nlon .or. mlat /= nlat) then
           call error_mesg ('ice_model_init',           &
            'incorrect resolution on restart', FATAL)
      endif

      call read_data ( unit, Ice%t_surf        )
      call read_data ( unit, Ice%thickness     )
      call read_data ( unit, Ice%albedo        )
      call read_data ( unit, Ice%albedo_vis_dir)
      call read_data ( unit, Ice%albedo_nir_dir)
      call read_data ( unit, Ice%albedo_vis_dif)
      call read_data ( unit, Ice%albedo_nir_dif)
      call read_data ( unit, Ice%rough_mom     )
      call read_data ( unit, Ice%rough_heat    )
      call read_data ( unit, Ice%rough_moist   )
      call close_file (unit)
      call nullify_domain ()

  !--- if no restart then no ice ---
   else

      need_ic = .true.
      Ice%t_surf      = TFREEZE   ! + temp_ice_freeze
      Ice%thickness   = 0.0       ! no ice initially
      Ice%albedo      = 0.14
      Ice%albedo_vis_dir = 0.14
      Ice%albedo_nir_dir = 0.14
      Ice%albedo_vis_dif = 0.14
      Ice%albedo_nir_dif = 0.14
      Ice%rough_mom   = 0.0004
      Ice%rough_heat  = 0.0004
      Ice%rough_moist = 0.0004
    ! Ice%ice_mask    = .false. ! no ice initially

    ! fixed roughness
      call fixed_ocean_roughness ( Ice%mask, Ice%rough_mom, &
                                   Ice%rough_heat, Ice%rough_moist )
   endif
endif

  ! initialize mask where ice exists
    Ice%ice_mask = Ice%mask .and. Ice%thickness .ge. thickness_min

    call ice_albedo_init (TFREEZE)
    
  ! analytic distribution with no ice
  ! melt all ice
    if (sst_method == "aqua_planet_1") then
        ice_method = 'none'
        Ice%ice_mask = .false.
        where (Ice%mask)
            Ice%t_surf = 27.*(1.-sin(max(min(1.5*Ice%lat,pi*0.5),-pi*0.5))**2) + TFREEZE
        endwhere
    else if (sst_method == "aqua_planet_2") then
        ice_method = 'none'
        Ice%ice_mask = .false.
        where (Ice%mask)
            Ice%t_surf = 27.*(1.-min(3.*abs(Ice%lat)/pi,1.)) + TFREEZE
        endwhere
    else if (sst_method == "aqua_planet_3") then
        ice_method = 'none'
        Ice%ice_mask = .false.
        where (Ice%mask)
            Ice%t_surf = 27.*(1.-sin(max(min(1.5*Ice%lat,pi*0.5),-pi*0.5))**4) + TFREEZE
        endwhere
    else if (sst_method == "aqua_planet_4") then
        ice_method = 'none'
        Ice%ice_mask = .false.
        where (Ice%mask)
            Ice%t_surf = max(min(1.5*Ice%lat,pi*0.5),-pi*0.5) ! use t_surf as work array
            Ice%t_surf = 27.*(1.-0.5*(sin(Ice%t_surf)**2+sin(Ice%t_surf)**4)) + TFREEZE
        endwhere
    else if (sst_method == "aqua_planet_5") then
        ice_method = 'none'
        Ice%ice_mask = .false.
        do j = js, je
        do i = is, ie
           if (Ice%mask(i,j)) then
              if (Ice%lat(i,j) > pi/36.) then
                  Ice%t_surf(i,j) = 27.*(1.-sin(min(90./55.*(Ice%lat(i,j)-pi/36.),pi*0.5))**2) + TFREEZE
              else
                  Ice%t_surf(i,j) = 27.*(1.-sin(max(90./65.*(Ice%lat(i,j)-pi/36.),-pi*0.5))**2) + TFREEZE
              endif
           endif
        enddo
        enddo
    else if (sst_method == "aqua_planet_6" .or. sst_method == "aqua_planet_7") then
        ice_method = 'none'
        Ice%ice_mask = .false.
        ! constants
        lon0 = 0.
        lond = 30.*pi/180.
        latd = 15.*pi/180.
        if (sst_method == "aqua_planet_6") amp = 1.
        if (sst_method == "aqua_planet_7") amp = 3.
        do j = js, je
        do i = is, ie
           if (Ice%mask(i,j)) then
               t_control = 27.*(1.-sin(max(min(1.5*Ice%lat(i,j),pi*0.5),-pi*0.5))**2) + TFREEZE
               dellon = Ice%lon(i,j)-lon0
               if (dellon >  pi) dellon = dellon - 2.*pi
               if (dellon < -pi) dellon = dellon + 2.*pi
               Ice%t_surf(i,j) = t_control + amp * cos(0.5*pi*min(max(dellon/lond,-1.),1.))**2 * &
                                                   cos(0.5*pi*min(max(Ice%lat(i,j)/latd,-1.),1.))**2
            endif
        enddo
        enddo
    else if (sst_method == "aqua_planet_8") then
        ice_method = 'none'
        Ice%ice_mask = .false.
        ! constants
        lon0 = 0.
        latd = 30.*pi/180.
        amp = 3.
        where (Ice%mask)
            Ice%t_surf = 27.*(1.-sin(max(min(1.5*Ice%lat,pi*0.5),-pi*0.5))**2) + TFREEZE + &
                         amp * cos(Ice%lon-lon0) * cos(0.5*pi*min(max(Ice%lat/latd,-1.),1.))**2
        endwhere
    
    ! -------------------------------------------------------------------------------------------
    ! GR: begin definition of custom SST profiles
    ! -------------------------------------------------------------------------------------------

    ! SST maximum centered at equator
    else if (sst_method == "aqua_planet_9") then
        ice_method = 'none'
        Ice%ice_mask = .false.

        ! Constants for SST profile shift
        min_SST = 0. ! minimum SST is chosen to be 7K lower than the peak, per Hsieh et al. (2020) 
        max_SST = 29. ! taken from HadISST zonal mean SST maximum, see hadisst_sst.clim.1986-2005.nc
        meridional_offset = 0 * pi / 180. ! convert from degrees to radians
        root_min = (-(pi / 2) + meridional_offset) / sine_wavenumber ! calculate minimum root (where SST = 0 degC)
        root_max = ((pi / 2) + meridional_offset) / sine_wavenumber ! calculate maximum root for filtering (where SST = 0 degC)
       
        do j = js, je
            do i = is, ie
                Ice%t_surf(i, j) = merge((max_SST - min_SST) * (1 - sin(sine_wavenumber * Ice%lat(i, j) - meridional_offset)**4) + min_SST + TFREEZE, & 
                                          min_SST + TFREEZE, &
                                         ((Ice%lat(i, j) .ge. root_min) .and. (Ice%lat(i, j) .le. root_max))) 
            enddo
        enddo
    
    ! SST maximum centered at equator, narrower SST maximum
    else if (sst_method == "aqua_planet_9b") then
        ice_method = 'none'
        Ice%ice_mask = .false.

        ! Constants for SST profile shift
        min_SST = 0. ! minimum SST is chosen to be 7K lower than the peak, per Hsieh et al. (2020) 
        max_SST = 29. ! taken from HadISST zonal mean SST maximum, see hadisst_sst.clim.1986-2005.nc
        
        meridional_offset = 0 * pi / 180. ! convert from degrees to radians
        root_min = (-(pi / 2) / sine_wavenumber + meridional_offset) ! calculate minimum root (where SST = 0 degC)
        root_max = ((pi / 2) / sine_wavenumber + meridional_offset) ! calculate maximum root for filtering (where SST = 0 degC)
       
        do j = js, je
            do i = is, ie
                Ice%t_surf(i, j) = merge((max_SST - min_SST) * (1 - sin(sine_wavenumber * Ice%lat(i, j) - meridional_offset)**2) + min_SST + TFREEZE, & 
                                          min_SST + TFREEZE, &
                                         ((Ice%lat(i, j) .ge. root_min) .and. (Ice%lat(i, j) .le. root_max))) 
            enddo
        enddo
    
    ! SST maximum centered at 10 N based on a sinusoid
    else if (sst_method == "aqua_planet_sine") then
        ice_method = 'none'
        Ice%ice_mask = .false.

        ! Constants for SST profile shift
        meridional_offset = latitude_of_maximum_SST * pi / 180. ! convert from degrees to radians
        root_min = (-(pi / 2) / sine_wavenumber + meridional_offset) ! calculate minimum root (where SST = 0 degC)
        root_max = ((pi / 2) / sine_wavenumber + meridional_offset) ! calculate maximum root for filtering (where SST = 0 degC)
        do j = js, je
            do i = is, ie
                Ice%t_surf(i, j) = merge((max_SST - min_SST) * (1 - sin(sine_wavenumber * (Ice%lat(i, j) - meridional_offset))**sine_exponent) + min_SST + TFREEZE, &
                                          min_SST + TFREEZE, &
                                          ((Ice%lat(i, j) .ge. root_min) .and. (Ice%lat(i, j) .le. root_max)))
            enddo
        enddo

    ! SST maximum centered at 10 N
    else if (sst_method == "aqua_planet_10N") then
        ice_method = 'none'
        Ice%ice_mask = .false.

        ! Constants for SST profile shift
        min_SST = 0. ! minimum SST is chosen to be 7K lower than the peak, per Hsieh et al. (2020) 
        max_SST = 29. ! taken from HadISST zonal mean SST maximum, see hadisst_sst.clim.1986-2005.nc
        period = 1.25
        meridional_offset = 10 * pi / 180. ! convert from degrees to radians
        root_min = (-(pi / 2) / period + meridional_offset) ! calculate minimum root (where SST = 0 degC)
        root_max = ((pi / 2) / period + meridional_offset) ! calculate maximum root for filtering (where SST = 0 degC)
       
        do j = js, je
            do i = is, ie
                Ice%t_surf(i, j) = merge((max_SST - min_SST) * (1 - sin(period * Ice%lat(i, j) - meridional_offset)**2) + min_SST + TFREEZE, & 
                                          min_SST + TFREEZE, &
                                         ((Ice%lat(i, j) .ge. root_min) .and. (Ice%lat(i, j) .le. root_max))) 
            enddo
        enddo
        
    ! SST maximum centered at 15 N
    else if (sst_method == "aqua_planet_15N") then
        ice_method = 'none'
        Ice%ice_mask = .false.

        ! Constants for SST profile shift
        min_SST = 0. ! minimum SST is chosen to be 7K lower than the peak, per Hsieh et al. (2020) 
        max_SST = 29. ! taken from HadISST zonal mean SST maximum, see hadisst_sst.clim.1986-2005.nc
        period = 1.25
        meridional_offset = 15 * pi / 180. ! convert from degrees to radians
        root_min = (-(pi / 2) / period + meridional_offset) ! calculate minimum root (where SST = 0 degC)
        root_max = ((pi / 2) / period + meridional_offset) ! calculate maximum root for filtering (where SST = 0 degC)
       
        do j = js, je
            do i = is, ie
                Ice%t_surf(i, j) = merge((max_SST - min_SST) * (1 - sin(period * Ice%lat(i, j) - meridional_offset)**2) + min_SST + TFREEZE, & 
                                          min_SST + TFREEZE, &
                                         ((Ice%lat(i, j) .ge. root_min) .and. (Ice%lat(i, j) .le. root_max))) 
            enddo
        enddo

    ! SST maximum centered at 20 N
    else if (sst_method == "aqua_planet_20N") then
        ice_method = 'none'
        Ice%ice_mask = .false.

        ! Constants for SST profile shift
        min_SST = 0. ! minimum SST is chosen to be 7K lower than the peak, per Hsieh et al. (2020) 
        max_SST = 29. ! taken from HadISST zonal mean SST maximum, see hadisst_sst.clim.1986-2005.nc
        period = 1.25
        meridional_offset = 20 * pi / 180. ! convert from degrees to radians
        root_min = (-(pi / 2) + meridional_offset) / period ! calculate minimum root (where SST = 0 degC)
        root_max = ((pi / 2) + meridional_offset) / period ! calculate maximum root for filtering (where SST = 0 degC)
       
        do j = js, je
            do i = is, ie
                Ice%t_surf(i, j) = merge((max_SST - min_SST) * (1 - sin(period * Ice%lat(i, j) - meridional_offset)**4) + min_SST + TFREEZE, & 
                                          min_SST + TFREEZE, &
                                         ((Ice%lat(i, j) .ge. root_min) .and. (Ice%lat(i, j) .le. root_max))) 
            enddo
        enddo
    
    ! Pseudo-constant SST
    else if (sst_method == "aqua_planet_11") then
        ice_method = 'none'
        Ice%ice_mask = .false.

        ! Constants for SST profile shift
        min_SST = 22. ! minimum SST is chosen to be 7K lower than the peak, per Hsieh et al. (2020) 
        max_SST = 29. ! taken from HadISST zonal mean SST maximum, see hadisst_sst.clim.1986-2005.nc
        period = 1.1
        meridional_offset = 0 * pi / 180. ! convert from degrees to radians
        root_min = (-(pi / 2) + meridional_offset) / period ! calculate minimum root (where SST = 0 degC)
        root_max = ((pi / 2) + meridional_offset) / period ! calculate maximum root for filtering (where SST = 0 degC)
       
        do j = js, je
            do i = is, ie
                Ice%t_surf(i, j) = merge((max_SST - min_SST) * (1 - sin(period * Ice%lat(i, j) - meridional_offset)**100) + min_SST + TFREEZE, & 
                                          min_SST + TFREEZE, &
                                         ((Ice%lat(i, j) .ge. root_min) .and. (Ice%lat(i, j) .le. root_max))) 
            enddo
        enddo
    
    ! SST flat profile with constricted warm pool
    else if (sst_method == "aqua_planet_12") then
        ice_method = 'none'
        Ice%ice_mask = .false.

        ! Constants for SST profile shift
        min_SST = -2. ! minimum SST is chosen to be 7K lower than the peak, per Hsieh et al. (2020) 
        max_SST = 29. ! taken from HadISST zonal mean SST maximum, see hadisst_sst.clim.1986-2005.nc
        period = 1.5
        meridional_offset = 15 * pi / 180. ! convert from degrees to radians
        root_min = (-(pi / 2) + meridional_offset) / period ! calculate minimum root (where SST = 0 degC)
        root_max = ((pi / 2) + meridional_offset) / period ! calculate maximum root for filtering (where SST = 0 degC)
       
        do j = js, je
            do i = is, ie
                Ice%t_surf(i, j) = merge((max_SST - min_SST) * (1 - sin(period * Ice%lat(i, j) - meridional_offset)**4) + min_SST + TFREEZE, & 
                                          min_SST + TFREEZE, &
                                         ((Ice%lat(i, j) .ge. root_min) .and. (Ice%lat(i, j) .le. root_max))) 
            enddo
        enddo

    ! Gaussian, 5-degree e-folding
    else if (sst_method == "aqua_planet_gauss") then
        ice_method = 'none'
        Ice%ice_mask = .false.

        ! Constants for SST profile shift
        meridional_offset = latitude_of_maximum_SST * pi / 180

        do j = js, je
            do i = is, ie
                Ice%t_surf(i, j) = max_SST * exp((-(Ice%lat(i, j) - meridional_offset)**2) / (2 * (e_folding_width * pi / 180) ** 2)) + TFREEZE
            enddo
        enddo
    endif


!----------------------------------------------------------

  if (trim(ice_method) == 'prognostic' .or. &
      trim(ice_method) == 'uniform') then
      if (trim(interp_method) == "conservative") then
          Amip_ice = amip_interp_new ( Ice%lon_bnd(:,1),     Ice%lat_bnd(1,:),  &
                         Ice%mask(:,:), interp_method = interp_method, &
                    use_climo=use_climo_ice, use_annual=use_annual_ice )
      else if(trim(interp_method) == "bilinear") then
          Amip_ice = amip_interp_new ( Ice%lon,     Ice%lat,          &
                         Ice%mask(:,:), interp_method = interp_method, &
                    use_climo=use_climo_ice, use_annual=use_annual_ice )
      else
          call error_mesg ('ice_model_init', 'interp_method should be '// &
                           'conservative or bilinear', FATAL)
      endif
      ! initialize ice (if needed)
      if (need_ic) then
         call prognostic_ice ( Ice )
      endif
  endif

!----------------------------------------------------------

  if (trim(sst_method) == 'specified') then
      if (trim(interp_method) == "conservative") then
          Amip_sst = amip_interp_new ( Ice%lon_bnd(:,1),     Ice%lat_bnd(1,:),  &
                         Ice%mask(:,:), interp_method = interp_method, &
                    use_climo=use_climo_sst, use_annual=use_annual_sst )
      else if(trim(interp_method) == "bilinear") then
          Amip_sst = amip_interp_new ( Ice%lon,     Ice%lat,          &
                         Ice%mask(:,:), interp_method = interp_method, &
                    use_climo=use_climo_sst, use_annual=use_annual_sst )
      else
          call error_mesg ('ice_model_init', 'interp_method should be '// &
                           'conservative or bilinear', FATAL)
      endif
      ! initialize sst (if needed)
      if (need_ic) then
         call prognostic_sst ( Ice )
      endif
  endif

print *, 'pe,count(ice,all,ocean)=',mpp_pe(),count(Ice%ice_mask),count(Ice%mask),count(Ice%mask .and. .not.Ice%ice_mask)

!miz
   id_flux_oh = register_diag_field ( mod_name, 'flux_oh', Atm%axes(1:2), &
                                      Ice%Time, 'flux_oh', 'W/m2' )
!miz
!----------------------------------------------------------

  module_is_initialized = .true.

!----------------------------------------------------------

 end subroutine ice_model_init

!######################################################################
!######## netcdf interface routines #########
!######################################################################

 subroutine get_grid_size ( Fields, nlon, nlat )
 type(fieldtype), intent(in)  :: Fields(:)
 integer,         intent(out) :: nlon, nlat
 integer :: i, j, dimsiz(4)
 character(len=128) :: name

  nlon = 0; nlat = 0
  do i = 1, size(Fields(:))
     do j=1,128; name(j:j)=' '; enddo
     call mpp_get_atts (Fields(i), name=name, siz=dimsiz)
       select case (trim(name))
          case ('geolon_t')
             nlon= dimsiz(1); nlat = dimsiz(2)
       end select
  enddo

 end subroutine get_grid_size

!----------------------------------------------------------------------

 subroutine read_grid_data ( unit, Fields, glonb, glatb, glon, glat, gmask )
 integer, intent(in) :: unit
 type(fieldtype), intent(in)  :: Fields(:)
 real,    intent(out) :: glonb(:), glatb(:)
 real,    intent(out) :: glon(:,:), glat(:,:)
 logical, intent(out) :: gmask(:,:)
 
 integer :: i, m, n
 character(len=128) :: name
 real, dimension(size(glon,1)+1,size(glon,2)+1) :: data2d

      m = size(glon,1);  n= size(glon,2)

      do i = 1, size(Fields(:))
         call mpp_get_atts(Fields(i), name=name)
         select case (trim(name))
            case ('geolon_t')
               call mpp_read(unit,Fields(i),glon)
               glon = glon*pi/180.
            case ('geolat_t')
               call mpp_read(unit,Fields(i),glat)
               glat = glat*pi/180.
            case ('geolon_vert_t')
               call mpp_read(unit,Fields(i),data2d)
               glonb = data2d(:,1)*pi/180.
            case('geolat_vert_t')
               call mpp_read(unit,Fields(i),data2d)
               glatb = data2d(1,:)*pi/180.
            case('wet')
               call mpp_read(unit,Fields(i),data2d(1:m,1:n))
               gmask = data2d(1:m,1:n) .gt. 0.50
         end select
      enddo

 end subroutine read_grid_data

!######################################################################

 subroutine ice_model_end ( Ice )
 type(ice_data_type), intent(inout) :: Ice
 integer :: unit
 character(len=64) :: fname='RESTART/ice_model.res.nc'

 if (.not.module_is_initialized) return
 if( do_netcdf_restart) then

    if(mpp_pe() == mpp_root_pe() ) then
       call error_mesg ('ice_model_mod', 'Writing NetCDF formatted restart file: RESTART/ice_model.res.nc', NOTE)
    endif   
    call write_data(fname, 'mlon', size(Ice%gmask,1), Ice%Domain)
    call write_data(fname, 'mlat', size(Ice%gmask,2), Ice%Domain)
    
    call write_data ( fname, 't_surf',         Ice%t_surf,         Ice%Domain )
    call write_data ( fname, 'thickness',      Ice%thickness,      Ice%Domain )
    call write_data ( fname, 'albedo',         Ice%albedo,         Ice%Domain )
    call write_data ( fname, 'albedo_vis_dir', Ice%albedo_vis_dir, Ice%Domain )
    call write_data ( fname, 'albedo_nir_dir', Ice%albedo_nir_dir, Ice%Domain )
    call write_data ( fname, 'albedo_vis_dif', Ice%albedo_vis_dif, Ice%Domain )
    call write_data ( fname, 'albedo_nir_dif', Ice%albedo_nir_dif, Ice%Domain )
    call write_data ( fname, 'rough_mom',      Ice%rough_mom,      Ice%Domain )
    call write_data ( fname, 'rough_heat',     Ice%rough_heat,     Ice%Domain )
    call write_data ( fname, 'rough_moist',    Ice%rough_moist,    Ice%Domain )
 else
    if (mpp_pe() == mpp_root_pe()) then
       call error_mesg ('ice_model_mod', 'Writing native formatted restart file.', NOTE)
    endif
    unit = open_restart_file ('RESTART/ice_model.res', 'write')
    if ( mpp_pe() == mpp_root_pe() ) then
       write (unit) size(Ice%gmask,1), size(Ice%gmask,2)
    endif

    call set_domain (Ice%Domain)
    call write_data ( unit, Ice%t_surf        )
    call write_data ( unit, Ice%thickness     )
    call write_data ( unit, Ice%albedo        )
    call write_data ( unit, Ice%albedo_vis_dir)
    call write_data ( unit, Ice%albedo_nir_dir)
    call write_data ( unit, Ice%albedo_vis_dif)
    call write_data ( unit, Ice%albedo_nir_dif)
    call write_data ( unit, Ice%rough_mom     )
    call write_data ( unit, Ice%rough_heat    )
    call write_data ( unit, Ice%rough_moist   )
    call close_file ( unit )
    call nullify_domain ()
 endif


  deallocate ( Ice%glon_bnd, Ice%glat_bnd, Ice%glon, Ice%glat, Ice%gmask )
  deallocate ( Ice%lon_bnd, Ice%lat_bnd, Ice%lon, Ice%lat, Ice%ice_mask,       &
               Ice%t_surf, Ice%albedo, Ice%albedo_vis_dir, Ice%albedo_nir_dir, &
               Ice%albedo_vis_dif, Ice%albedo_nir_dif, Ice%rough_mom,          &
               Ice%rough_heat, Ice%rough_moist, Ice%thickness, Ice%mask        )

  module_is_initialized = .false.

 end subroutine ice_model_end

!######################################################################
!               Routines added for computing then
!            mid-point of grid boxes with cubed sphere
!######################################################################

function is_latlon ( lon, lat )
real, intent(in) :: lon(:,:), lat(:,:)

! Determines if the latitude/longitude values form a
! regular latitude/longitude grid (or something else).
!

logical :: is_latlon
integer :: i, j 

  is_latlon = .true.

  do j = 2, size(lon,2)
  do i = 1, size(lon,1)
     if (lon(i,j) .ne. lon(i,1)) then
        is_latlon = .false.
        return
     endif
  enddo
  enddo

  do j = 1, size(lat,2)
  do i = 2, size(lat,1)
     if (lat(i,j) .ne. lat(1,j)) then
        is_latlon = .false.
        return
     endif
  enddo
  enddo

end function is_latlon

!########################################################################

  subroutine get_cell_center (lonb, latb, lon, lat)
    !------------------------------------------------------------------!
    ! calculate cell center (lon,lat)                                  !
    ! by averaging cell corner locations (lonb,latb) in Cartesian coor !
    !------------------------------------------------------------------!
    real, dimension(:,:), intent(in)  :: lonb, latb
    real, dimension(:,:), intent(out) :: lon, lat
    !------------------------------------------------------------------!
    ! local variables                                                  !
    !------------------------------------------------------------------!
    real :: sph_coor(2), xyz_corner(3,size(lonb,1),size(latb,2)), xyz_center(3), abs_center
    integer :: i,j, nx,ny

    nx = size(lonb,1)
    ny = size(latb,2)

    do j=1,ny
       do i=1,nx
          sph_coor(1)=lonb(i,j)
          sph_coor(2)=latb(i,j)
          call latlon2xyz(sph_coor, xyz_corner(1,i,j))
       enddo
    enddo
    do j=1,ny-1
       do i=1,nx-1
          xyz_center(:)=0.25*(xyz_corner(:,i  ,j  )+xyz_corner(:,i+1,j  )  &
                             +xyz_corner(:,i  ,j+1)+xyz_corner(:,i+1,j+1))
          abs_center=xyz_center(1)*xyz_center(1)                           &
                    +xyz_center(2)*xyz_center(2)                           &
                    +xyz_center(3)*xyz_center(3)
          xyz_center(:)=xyz_center(:)/abs_center
          call xyz2latlon(xyz_center, sph_coor)
          lon(i,j)=sph_coor(1)
          lat(i,j)=sph_coor(2)
       enddo
    enddo

  end subroutine get_cell_center

!########################################################################

  subroutine latlon2xyz(sph_coor, xyz_coor)
    !------------------------------------------------------------------!
    ! calculate cartesian coordinates from spherical coordinates       !
    !                                                                  !
    ! input:                                                           !
    ! sph_coor [rad]   latlon coordinate                               !
    !                                                                  !
    ! output:                                                          !
    ! xyz_coor [1]     normalized cartesian vector                     !
    !------------------------------------------------------------------!
    real, dimension(2), intent(in)    :: sph_coor
    real, dimension(3), intent(inout) :: xyz_coor

    xyz_coor(1) = cos(sph_coor(2)) * cos(sph_coor(1))
    xyz_coor(2) = cos(sph_coor(2)) * sin(sph_coor(1))
    xyz_coor(3) = sin(sph_coor(2))

  end subroutine latlon2xyz
  !====================================================================!
  subroutine xyz2latlon(xyz_coor, sph_coor)
    !------------------------------------------------------------------!
    ! calculate spherical coordinates from cartesian coordinates       !
    !                                                                  !
    ! input:                                                           !
    ! xyz_coor [1]     normalized cartesian vector                     !
    !                                                                  !
    ! output:                                                          !
    ! sph_coor [rad]   latlon coordinate                               !
    !------------------------------------------------------------------!

    real, dimension(3), intent(in)    :: xyz_coor
    real, dimension(2), intent(inout) :: sph_coor
    !------------------------------------------------------------------!
    ! local variables                                                  !
    !------------------------------------------------------------------!
    real :: coslat, radius
    real, parameter :: epsilon=1.e-7
    integer :: i,j


    radius=sqrt(xyz_coor(1)*xyz_coor(1) &
               +xyz_coor(2)*xyz_coor(2) &
               +xyz_coor(3)*xyz_coor(3))

    sph_coor(1)=atan2(xyz_coor(2),xyz_coor(1))
    if (sph_coor(1)<0.) sph_coor(1)=sph_coor(1)+2*PI
    sph_coor(2)=asin(xyz_coor(3)/radius)

  end subroutine xyz2latlon

!######################################################################
 
end module ice_model_mod

