#include <fms_platform.h>

! ============================================================================
! This null land model was created from land_lad, not land_lad2.
! It contains all inferfaces expected by the simple coupler.
! Upon initialization, it allocates fields referenced by coupler
! and fills them with non-NAN values.
! All routines called after initialization are do-nothing routines
! that immediately return.
! This version was not tested with the full coupler,
! nor was it tested as a replacement for land_lad2.
! ============================================================================

module land_model_mod

  use time_manager_mod, only: time_type

  use mpp_mod,          only: mpp_chksum

  use mpp_domains_mod,  only: domain2d, mpp_get_layout, mpp_define_layout, &
                              mpp_define_domains, mpp_get_compute_domain,  &
                              CYCLIC_GLOBAL_DOMAIN, mpp_define_mosaic

  use fms_mod,          only: write_version_number, error_mesg, FATAL, &
                              mpp_npes, stdout

  use constants_mod,   only:  radius, pi 

  use tracer_manager_mod, only: NO_TRACER

implicit none
private

! ==== public interfaces =====================================================
public land_model_init          ! initialize the land model
public land_model_end           ! finish land model calculations
public land_model_restart       ! dummy routines
public update_land_model_fast   ! fast time-scale update of the land
public update_land_model_slow   ! slow time-scale update of the land

public Lnd_stock_pe             ! calculate and return total amount of requested quantitiy per PE
public land_data_type_chksum, atm_lnd_bnd_type_chksum
! ==== end of public interfaces ==============================================

! ==== public data type =====================================================
public land_type, land_data_type, atmos_land_boundary_type

! ==== generic interfaces ====================================================
! <INTERFACE NAME="land_model_init">
!   <OVERVIEW>
!     Initializes the state of the land model.
!   </OVERVIEW>
!   <DESCRIPTION>
!    Allocates fields referenced by coupler and fills them with non-NAN values.
!   </DESCRIPTION>
interface land_model_init
   module procedure init_land_with_area
   module procedure init_land_with_xgrid
end interface
! </INTERFACE>


! <TYPE NAME="land_type">
!   <DESCRIPTION>
!     The type describing the state of the land model. It is private to this
!     module and is basically a remnant of the previous design. There is only one
!     variable of this type, called "theLand", which is used in the model to
!     hold the data.
!     Strictly speaking, this type is not necessary, but there is no harm in
!     having it and it is possible to image a situation where a single entity
!     describing the land-surface model state would be useful.
!   </DESCRIPTION>
type land_type

   private  ! it's an opaque type: all data members are private
   
!   <DATA NAME="domain" TYPE="domain2d" DIM="2">
!     Domain of calculations
!   </DATA>
!   <DATA NAME="is" TYPE="integer">
!     Domain bound
!   </DATA>
!   <DATA NAME="ie" TYPE="integer">
!     Domain bound
!   </DATA>
!   <DATA NAME="js" TYPE="integer">
!     Domain bound
!   </DATA>
!   <DATA NAME="je" TYPE="integer">
!     Domain bound
!   </DATA>
!   <DATA NAME="n_tiles" TYPE="integer">
!     Domain bound
!   </DATA>
   type(domain2d)         :: domain     ! domain of calculations
   integer :: is,ie,js,je,n_tiles       ! domain bounds, for convenience

   ! the values below is just a quick fix to start running Land with original
   ! exchange grid, since it assumes something like that to be present. It also 
   ! assumes the number of tiles to be fixed.
   !   Of course, in general case there is no such thing as "latitude/longitude"
   ! boundaries of the cell, since cell boundaries do not have to be parallel
   ! to the lat/lon coordinate axes

!   <DATA NAME="blon" TYPE="real, pointer" DIM="2">
!     Longitude domain corners
!   </DATA>
!   <DATA NAME="blat" TYPE="real, pointer" DIM="2">
!     Latitude domain corners
!   </DATA>
!   <DATA NAME="mask" TYPE="logical, pointer" DIM="2">
!     Land mask, true where there is land
!   </DATA>
   real, _ALLOCATABLE          :: blon(:,:) _NULL  ! longitude corners of our domain
   real, _ALLOCATABLE          :: blat(:,:) _NULL  ! latitude corners of our domain 
   logical, _ALLOCATABLE       :: mask(:,:) _NULL  ! land mask (true where ther is some land)

end type land_type
! </TYPE>


type :: atmos_land_boundary_type
   real, dimension(:,:,:), pointer :: &
        t_flux =>NULL(),  &
        lw_flux =>NULL(), &
        lwdn_flux => NULL(), &
        sw_flux =>NULL(), &
        sw_flux_down_vis_dir =>NULL(), &
        sw_flux_down_total_dir =>NULL(), &
        sw_flux_down_vis_dif =>NULL(), &
        sw_flux_down_total_dif =>NULL(), &
        lprec =>NULL(),   &
        fprec =>NULL(),   &
        tprec =>NULL()

   real, dimension(:,:,:), pointer :: &
        dhdt =>NULL(),    &
        drdt =>NULL()

   real, dimension(:,:,:), pointer :: &
        cd_m => NULL(),      &
        cd_t => NULL(),      &
        ustar => NULL(),     &
        bstar => NULL(),     &
        wind => NULL(),      &
        z_bot => NULL(),     &
        drag_q =>NULL(),     &
        p_surf =>NULL()

   real, dimension(:,:,:,:), pointer :: & ! (lon, lat, tile,tracer)
        tr_flux => NULL(), &
        dfdtr   => NULL()

   integer :: xtype

end type atmos_land_boundary_type

type :: land_data_type

   type(domain2d) :: domain

   real, pointer, dimension(:,:,:)   :: &  ! (lon, lat, tile)
        tile_size =>NULL(),       &
        t_surf =>NULL(),          &
        t_ca =>NULL(),            &
        albedo =>NULL(),          &
        albedo_vis_dir =>NULL(),  &
        albedo_nir_dir =>NULL(),  &
        albedo_vis_dif =>NULL(),  &
        albedo_nir_dif =>NULL(),  &
        rough_mom =>NULL(),       &
        rough_heat =>NULL(),      &
        rough_scale =>NULL()

   real, pointer, dimension(:,:,:,:) :: tr => NULL() ! (lon, lat, tile, tracer)

   real, pointer, dimension(:,:) :: &  ! (lon, lat)
        discharge           => NULL(),  &
        discharge_heat      => NULL(),  &
        discharge_snow      => NULL(),  &
        discharge_snow_heat => NULL()

   logical, pointer, dimension(:,:,:) :: mask =>NULL()  ! true if land

   logical, pointer, dimension(:,:) :: maskmap =>NULL()

   integer :: axes(2)      ! axes IDs for diagnostics  
   logical :: pe

end type land_data_type

! ==== some names, for information only ======================================
logical                       :: module_is_initialized = .FALSE.
character(len=*),   parameter :: module_name = 'land_mod'
character(len=128), parameter :: version     = '$Id: land_model.F90,v 19.0.6.1 2012/07/18 14:55:25 pjp Exp $'
character(len=128), parameter :: tagname     = '$Name: land_lad_null_for_simple_coupler_pjp $'


! ==== local module variables ================================================
integer                 :: n_tiles = 1  ! number of tiles
type(land_type),save    :: theLand  ! internal state of the model
integer :: isphum ! number of specific humidity in the tracer array, or NO_TRACER

contains ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-


! <SUBROUTINE NAME="init_land_with_xgrid" INTERFACE="land_model_init">
!   <OVERVIEW>
!     Initializes the land model
!   </OVERVIEW>
!   <DESCRIPTION>
!     The state of the land model is initialized using the grid description
!     file as an input. This routine reads land grid corners and area of
!     land from a grid description file.
!   </DESCRIPTION>
!   <PUBLICROUTINE INTERFACE="land_model_init">
subroutine init_land_with_xgrid &
     (atmos2land, bnd, time_init, time, dt_fast, dt_slow, &
      glonb, glatb, atmos_domain)

  type(atmos_land_boundary_type), intent(inout) :: atmos2land ! land boundary data
  type(land_data_type), intent(inout) :: bnd ! land boundary data
  type(time_type), intent(in)    :: time_init ! initial time of simulation (?)
  type(time_type), intent(in)    :: time      ! current time
  type(time_type), intent(in)    :: dt_fast   ! fast time step
  type(time_type), intent(in)    :: dt_slow   ! slow time step
  real,            intent(in), optional    :: glonb(:,:), glatb(:,:)
  type(domain2d),  intent(in), target, optional :: atmos_domain ! domain of computations
!   </PUBLICROUTINE>

  ! ---- local vars ----------------------------------------------------------
  real, allocatable :: glon(:,:)    ! lon centers of global grid
  real, allocatable :: glat(:,:)    ! lat centers of global grid
  real, allocatable :: gfrac    (:,:) ! fraction of land in cells
  real, allocatable :: gcellarea(:,:) ! area of land cells
  real, allocatable :: tmpx(:,:), tmpy(:,:)
  real, allocatable :: tmp_latb(:), tmp_lat(:)
  real, allocatable :: geo_lonv(:,:)
  real, allocatable :: geo_latv(:,:)
  real, allocatable :: xgrid_area(:)
  integer, allocatable :: i1(:), j1(:)
  integer, allocatable :: i2(:), j2(:)

  integer :: nlonb, nlatb  ! number of the cell bounds in lon and lat directions
  integer :: nlon,  nlat   ! number of cells in lon and lat directions
  integer :: siz(4)        ! used to store field size
  integer,save :: ntiles=6 ! number of tiles in land mosaic, should be 1 or 6.
  integer :: nfile_axl     ! number of atmosXland exchange grid file. Can be more than one.
  integer :: nxgrid        ! number of exchange grid cell.
  integer :: i, j, n, m, tile, ind, digit

  character(len=256) :: grid_file   = "INPUT/grid_spec.nc"
  character(len=256) :: tile_file
  character(len=256) :: land_mosaic   ! land mosaic file
  character(len=256) :: axl_file      ! atmosXland exchange grid file

  
! initialize grid info from input arguments
  nlonb = size(glonb,1)
  nlatb = size(glonb,2)
  nlon = nlonb-1
  nlat = nlatb-1
  allocate( glon(nlon,nlat), glat(nlon,nlat), gcellarea(nlon,nlat), gfrac(nlon,nlat) )

! approximate coordinates of grid cell centers
  glon = (glonb(1:nlon,1:nlat)+glonb(2:nlon+1,1:nlat)+glonb(1:nlon,2:nlat+1)+glonb(2:nlon+1,2:nlat+1))/4.
  glat = (glatb(1:nlon,1:nlat)+glatb(2:nlon+1,1:nlat)+glatb(1:nlon,2:nlat+1)+glatb(2:nlon+1,2:nlat+1))/4.

! initialize with no land
  gfrac = 0.0 ! aqua-planet

! uniform area
  gcellarea = 1./real(6*(nlonb-1)*(nlatb-1))
  ! convert relative area to absolute value, m2
  gcellarea = gcellarea*4*pi*radius**2

! initialize the land model using the global grid we just obtained
  call init_land_with_area &
       ( atmos2land, bnd, glonb, glatb, gcellarea, gfrac, time, dt_fast, dt_slow, &
       atmos_domain, glon=glon, glat=glat, numtiles=ntiles )

! deallocate memory that we are not going to use any more
  deallocate ( glon, glat, gfrac, gcellarea )


end subroutine init_land_with_xgrid
! </SUBROUTINE>


! <SUBROUTINE NAME="init_land_with_area" INTERFACE="land_model_init">
!   <OVERVIEW>
!     Initializes the land model
!   </OVERVIEW>
!   <DESCRIPTION>
!     Initializes the land state using land area for each of the grid points.
!   </DESCRIPTION>
!   <PUBLICROUTINE INTERFACE="land_model_init">
subroutine init_land_with_area &
     ( atmos2land, bnd, gblon, gblat, gcellarea, gfrac, time, dt_fast, dt_slow, domain, &
     glon, glat, numtiles )

  type(atmos_land_boundary_type), intent(inout) :: atmos2land ! land boundary data
  type(land_data_type), intent(inout) :: bnd        ! state to update
  real,              intent(in) :: gblon(:,:)! lon corners of the grid cells
  real,              intent(in) :: gblat(:,:)! lat corners of the grid cells
  real,              intent(in) :: gcellarea(:,:) ! full area of the grid cells
  real,              intent(in) :: gfrac(:,:)     ! fraction of land in the grid cell
  type(time_type),   intent(in) :: time    ! current time
  type(time_type),   intent(in) :: dt_fast ! fast time step
  type(time_type),   intent(in) :: dt_slow ! slow time step
  type(domain2d),    intent(in), optional :: domain  ! domain of computations
  real,              intent(in), optional :: glon(:,:), glat(:,:) ! centers
  integer,           intent(in), optional :: numtiles
                          ! of the grid cells
!   </PUBLICROUTINE>

  ! ---- local vars ----------------------------------------------------------
  integer :: layout(2) ! layout of our domains
  integer :: nlon, nlat
  integer :: is,ie,js,je,k
  integer :: ntprog   ! numbers of tracers
  integer :: ntiles ! number of tiles for land mosaic

  module_is_initialized = .TRUE.

  ! write the version and tagname to the logfile
  call write_version_number(version, tagname)

! null version of land has no tracers
  isphum = NO_TRACER

  ! get the size of the global grid
  nlon = size(gfrac, 1)
  nlat = size(gfrac, 2)
  ! number of land mosaic tiles
  ntiles = 1
  if (present(numtiles)) ntiles = numtiles
  if (mod(mpp_npes(),ntiles) .NE. 0) call error_mesg('land_model_init', &
      'Number of processors must be a multiple of number of mosaic tiles', FATAL)
  ! get the processor layout information, either from upper-layer module
  ! decomposition, or define according to the grid size
  if(present(domain)) then
     call mpp_get_layout (domain, layout)
  else
     call mpp_define_layout &
          ((/1,nlon, 1,nlat/),mpp_npes()/ntiles,layout)
  endif

  ! create our computation domain according to obtained processor layout 
  if (ntiles .EQ. 1) then
     call mpp_define_domains ( &
          (/1,nlon, 1, nlat/),           &  ! global grid size
          layout,                        &  ! layout for our domains
          theLand%domain,                &  ! domain to define
          xflags = CYCLIC_GLOBAL_DOMAIN, name = 'LAND MODEL' )
  else if (ntiles .EQ. 6) then
     call define_cube_mosaic ( theLand%domain, nlon, nlat, layout )
  else
     call error_mesg('land_model_init',  &
        'number of mosaic tiles should be 1 or 6', FATAL)
  endif

  ! get the size of our computation domain
  call mpp_get_compute_domain ( theLand%domain, is,ie,js,je )

  theLand%is=is; theLand%ie=ie; theLand%js=js; theLand%je=je
  theLand%n_tiles=n_tiles

  ! allocate storage for the boundary data 

  allocate ( &
       bnd % tile_size      (is:ie,js:je,n_tiles), &
       bnd % t_surf         (is:ie,js:je,n_tiles), &
       bnd % t_ca           (is:ie,js:je,n_tiles), &
       bnd % albedo         (is:ie,js:je,n_tiles), &
       bnd % albedo_vis_dir (is:ie,js:je,n_tiles), &
       bnd % albedo_nir_dir (is:ie,js:je,n_tiles), &
       bnd % albedo_vis_dif (is:ie,js:je,n_tiles), &
       bnd % albedo_nir_dif (is:ie,js:je,n_tiles), &
       bnd % rough_mom      (is:ie,js:je,n_tiles), &
       bnd % rough_heat     (is:ie,js:je,n_tiles), &
       bnd % rough_scale    (is:ie,js:je,n_tiles), &
       bnd % mask           (is:ie,js:je,n_tiles), &
       bnd % discharge             (is:ie,js:je),  &
       bnd % discharge_heat        (is:ie,js:je),  &
       bnd % discharge_snow        (is:ie,js:je),  &
       bnd % discharge_snow_heat   (is:ie,js:je),  &
       bnd % tr             (is:ie,js:je,n_tiles,1))

  bnd % tile_size      = 0.0
  bnd % t_surf         = 280.
  bnd % t_ca           = 280.
  bnd % albedo         = 0.0
  bnd % albedo_vis_dir = 0.0
  bnd % albedo_nir_dir = 0.0
  bnd % albedo_vis_dif = 0.0
  bnd % albedo_nir_dif = 0.0
  bnd % rough_mom      = 0.0
  bnd % rough_heat     = 0.0
  bnd % rough_scale    = 0.0
  bnd % mask           = .FALSE.
  bnd % discharge      = 0.0
  bnd % discharge_heat = 0.0
  bnd % discharge_snow = 0.0
  bnd % discharge_snow_heat = 0.0
  bnd % tr             = 0.0

  allocate( &
       atmos2land % t_flux  (is:ie,js:je,n_tiles), &
       atmos2land % lw_flux (is:ie,js:je,n_tiles), &
       atmos2land % sw_flux (is:ie,js:je,n_tiles), &
       atmos2land % lprec   (is:ie,js:je,n_tiles), &
       atmos2land % fprec   (is:ie,js:je,n_tiles), &
       atmos2land % tprec   (is:ie,js:je,n_tiles), &
       atmos2land % dhdt    (is:ie,js:je,n_tiles), &
       atmos2land % drdt    (is:ie,js:je,n_tiles), &
       atmos2land % drag_q  (is:ie,js:je,n_tiles), &
       atmos2land % p_surf  (is:ie,js:je,n_tiles), &
       atmos2land % sw_flux_down_vis_dir   (is:ie,js:je,n_tiles), &
       atmos2land % sw_flux_down_total_dir (is:ie,js:je,n_tiles), &
       atmos2land % sw_flux_down_vis_dif   (is:ie,js:je,n_tiles), &
       atmos2land % sw_flux_down_total_dif (is:ie,js:je,n_tiles), &
       atmos2land % tr_flux(is:ie,js:je,n_tiles,1), &
       atmos2land % dfdtr(is:ie,js:je,n_tiles,1), &
       atmos2land % lwdn_flux(is:ie,js:je,n_tiles), &
       atmos2land % cd_m(is:ie,js:je,n_tiles), &
       atmos2land % cd_t(is:ie,js:je,n_tiles), &
       atmos2land % bstar(is:ie,js:je,n_tiles), &
       atmos2land % ustar(is:ie,js:je,n_tiles), &
       atmos2land % wind(is:ie,js:je,n_tiles), &
       atmos2land % z_bot(is:ie,js:je,n_tiles) )

       atmos2land % t_flux  = 0.0
       atmos2land % lw_flux = 0.0
       atmos2land % sw_flux = 0.0
       atmos2land % lprec   = 0.0
       atmos2land % fprec   = 0.0
       atmos2land % tprec   = 0.0
       atmos2land % dhdt    = 0.0
       atmos2land % drdt    = 0.0
       atmos2land % drag_q  = 0.0
       atmos2land % p_surf  = 1.0e5  
       atmos2land % sw_flux_down_vis_dir   = 0.0
       atmos2land % sw_flux_down_total_dir = 0.0
       atmos2land % sw_flux_down_vis_dif   = 0.0
       atmos2land % sw_flux_down_total_dif = 0.0
       atmos2land % tr_flux = 0.0
       atmos2land % dfdtr   = 0.0
       atmos2land % lwdn_flux = 0.0
       atmos2land % cd_m  = 0.0
       atmos2land % cd_t  = 0.0
       atmos2land % bstar = 0.0
       atmos2land % ustar = 0.0
       atmos2land % wind  = 0.0
       atmos2land % z_bot = 0.0

  ! set up boundary values that we know already
  bnd % domain = theLand % domain
  do k = 1, size(bnd%tile_size,3)
     where (gfrac(is:ie,js:je)>0)
        bnd % tile_size (:,:,k) = 1.0/n_tiles
        bnd % mask      (:,:,k) = .true.
     elsewhere
        bnd % tile_size (:,:,k) = 0.0
        bnd % mask      (:,:,k) = .false.
     endwhere
  enddo

  ! init grid variables: this is a quick fix since it is not yet
  ! clear how to initialize land properties with arbitrary grid
  allocate(theLand%blon(is:ie+1,js:je+1), theLand%blat(is:ie+1,js:je+1), theLand%mask(is:ie,js:je))
  theLand%blon(:,:)   = gblon(is:ie+1,js:je+1)
  theLand%blat(:,:)   = gblat(is:ie+1,js:je+1)
  theLand%mask(:,:)   = any( bnd % mask, DIM=3 )

end subroutine init_land_with_area
! </SUBROUTINE>


subroutine land_model_end ( atmos2land, bnd )
  type(atmos_land_boundary_type), intent(inout) :: atmos2land
  type(land_data_type), intent(inout) :: bnd

  module_is_initialized = .FALSE.
  deallocate ( theLand%blon, theLand%blat, theLand%mask )

  call deallocate_boundary_data ( atmos2land, bnd )
  
end subroutine land_model_end

subroutine land_model_restart(timestamp)
  character(*), intent(in), optional :: timestamp ! timestamp to add to the file name

  return ! do nothing routine

end subroutine land_model_restart

subroutine update_land_model_fast ( atmos2land, bnd )

  type(atmos_land_boundary_type), intent(inout) :: atmos2land
  type(land_data_type),  intent(inout) :: bnd

  return ! do nothing routine

end subroutine update_land_model_fast


subroutine update_land_model_slow ( atmos2land, bnd )
  type(atmos_land_boundary_type), intent(inout) :: atmos2land
  type(land_data_type), intent(inout) :: bnd

  return ! do nothing routine

end subroutine update_land_model_slow
! </SUBROUTINE>

! ============================================================================
subroutine Lnd_stock_pe(bnd,index,value)
  type(land_data_type), intent(in)  :: bnd
  integer             , intent(in)  :: index ! ID of the stock to calculate
  real                , intent(out) :: value ! calculated value of the stock

  real :: soil_value, vegn_value

  value = 0.0

end subroutine Lnd_stock_pe
! ============================================================================

! <SUBROUTINE NAME="define_cube_mosaic">
subroutine define_cube_mosaic ( Domain, nx, ny, layout )
type(domain2d), intent(inout) :: Domain
integer, intent(in) :: nx, ny, layout(2)
integer, parameter :: ntiles = 6, num_contact = 12, ng = 0
integer :: global_indices(4,ntiles), layout_2d(2,ntiles)
integer :: pe_start(ntiles), pe_end(ntiles)
integer :: east(4), west(4), north(4), south(4)
integer :: tile1(num_contact), tile2(num_contact),  &
           index1(num_contact,4), index2(num_contact,4)
integer :: n

    do n = 1, ntiles
       global_indices(:,n) = (/ 1, nx, 1, ny /)
       layout_2d     (:,n) = layout
       pe_start        (n) = (n-1)*layout(1)*layout(2)
       pe_end          (n) =     n*layout(1)*layout(2) - 1
    enddo

  ! istart, iend, jstart, jend for each face
    west =  (/  1,  1,  1, ny /)
    east =  (/ nx, nx,  1, ny /)
    south = (/  1, nx,  1,  1 /)
    north = (/  1, nx, ny, ny /)

  ! define the 12 contact lines bewteen the faces
    tile1( 1) = 1; index1( 1,:) = east;   tile2( 1) = 2; index2( 1,:) = west
    tile1( 2) = 1; index1( 2,:) = north;  tile2( 2) = 3; index2( 2,:) = west
    tile1( 3) = 1; index1( 3,:) = west;   tile2( 3) = 5; index2( 3,:) = north
    tile1( 4) = 1; index1( 4,:) = south;  tile2( 4) = 6; index2( 4,:) = north
    tile1( 5) = 2; index1( 5,:) = north;  tile2( 5) = 3; index2( 5,:) = south
    tile1( 6) = 2; index1( 6,:) = east;   tile2( 6) = 4; index2( 6,:) = south
    tile1( 7) = 2; index1( 7,:) = south;  tile2( 7) = 6; index2( 7,:) = east
    tile1( 8) = 3; index1( 8,:) = east;   tile2( 8) = 4; index2( 8,:) = west
    tile1( 9) = 3; index1( 9,:) = north;  tile2( 9) = 5; index2( 9,:) = west
    tile1(10) = 4; index1(10,:) = north;  tile2(10) = 5; index2(10,:) = south
    tile1(11) = 4; index1(11,:) = east;   tile2(11) = 6; index2(11,:) = south
    tile1(12) = 5; index1(12,:) = east;   tile2(12) = 6; index2(12,:) = west

  ! create the domain2d variable
    call mpp_define_mosaic ( global_indices, layout_2d, Domain,                  &
                             ntiles, num_contact, tile1, tile2,                  &
                             index1(:,1), index1(:,2), index1(:,3), index1(:,4), &
                             index2(:,1), index2(:,2), index2(:,3), index2(:,4), &
                             pe_start=pe_start, pe_end=pe_end, symmetry=.true.,  &
                             shalo = ng, nhalo = ng, whalo = ng, ehalo = ng,     &
                             name = "Land Model Cubic-Grid" )

end subroutine define_cube_mosaic
! </SUBROUTINE>
! ============================================================================
subroutine land_data_type_chksum(id, timestep, land)

    character(len=*), intent(in) :: id
    integer         , intent(in) :: timestep
    type(land_data_type), intent(in) :: land
    integer :: outunit
    
    outunit = stdout()
        
    write(outunit,*) "BEGIN CHECKSUM(land_type):: ", id, timestep
    write(outunit,100) 'land%tile_size         ',mpp_chksum(land%tile_size)
    write(outunit,100) 'land%t_surf            ',mpp_chksum(land%t_surf)
    write(outunit,100) 'land%t_ca              ',mpp_chksum(land%t_ca)
    write(outunit,100) 'land%albedo            ',mpp_chksum(land%albedo)
    write(outunit,100) 'land%albedo_vis_dir    ',mpp_chksum(land%albedo_vis_dir)
    write(outunit,100) 'land%albedo_nir_dir    ',mpp_chksum(land%albedo_nir_dir)
    write(outunit,100) 'land%albedo_vis_dif    ',mpp_chksum(land%albedo_vis_dif)
    write(outunit,100) 'land%albedo_nir_dif    ',mpp_chksum(land%albedo_nir_dif)
    write(outunit,100) 'land%rough_mom         ',mpp_chksum(land%rough_mom)
    write(outunit,100) 'land%rough_heat        ',mpp_chksum(land%rough_heat)
    write(outunit,100) 'land%rough_scale       ',mpp_chksum(land%rough_scale)
    write(outunit,100) 'land%discharge         ',mpp_chksum(land%discharge)
    write(outunit,100) 'land%discharge_snow    ',mpp_chksum(land%discharge_snow)

100 FORMAT("   CHECKSUM::",A32," = ",Z20)

end subroutine land_data_type_chksum
! ============================================================================
subroutine atm_lnd_bnd_type_chksum(id, timestep, albt)

    character(len=*), intent(in) :: id
    integer         , intent(in) :: timestep
    type(atmos_land_boundary_type), intent(in) :: albt
    integer :: outunit

    outunit = stdout()

    write(outunit,*) "BEGIN CHECKSUM(atmos_land_boundary_type):: ", id, timestep
    write(outunit,100) 'albt%t_flux                ', mpp_chksum( albt%t_flux)
    write(outunit,100) 'albt%lw_flux               ', mpp_chksum( albt%lw_flux)
    write(outunit,100) 'albt%lwdn_flux             ', mpp_chksum( albt%lwdn_flux)
    write(outunit,100) 'albt%sw_flux               ', mpp_chksum( albt%sw_flux)
    write(outunit,100) 'albt%sw_flux_down_vis_dir  ', mpp_chksum( albt%sw_flux_down_vis_dir)
    write(outunit,100) 'albt%sw_flux_down_total_dir', mpp_chksum( albt%sw_flux_down_total_dir)
    write(outunit,100) 'albt%sw_flux_down_vis_dif  ', mpp_chksum( albt%sw_flux_down_vis_dif)
    write(outunit,100) 'albt%sw_flux_down_total_dif', mpp_chksum( albt%sw_flux_down_total_dif)
    write(outunit,100) 'albt%lprec                 ', mpp_chksum( albt%lprec)
    write(outunit,100) 'albt%fprec                 ', mpp_chksum( albt%fprec)
    write(outunit,100) 'albt%tprec                 ', mpp_chksum( albt%tprec)
    write(outunit,100) 'albt%dhdt                  ', mpp_chksum( albt%dhdt)
    write(outunit,100) 'albt%drdt                  ', mpp_chksum( albt%drdt)
    write(outunit,100) 'albt%cd_m                  ', mpp_chksum( albt%cd_m)
    write(outunit,100) 'albt%cd_t                  ', mpp_chksum( albt%cd_t)
    write(outunit,100) 'albt%ustar                 ', mpp_chksum( albt%ustar)
    write(outunit,100) 'albt%bstar                 ', mpp_chksum( albt%bstar)
    write(outunit,100) 'albt%wind                  ', mpp_chksum( albt%wind)
    write(outunit,100) 'albt%z_bot                 ', mpp_chksum( albt%z_bot)
    write(outunit,100) 'albt%drag_q                ', mpp_chksum( albt%drag_q)
    write(outunit,100) 'albt%p_surf                ', mpp_chksum( albt%p_surf)

100 FORMAT("   CHECKSUM::",A32," = ",Z20)

end subroutine atm_lnd_bnd_type_chksum
! ============================================================================

subroutine deallocate_boundary_data ( a2l, bnd )
  type(atmos_land_boundary_type), intent(inout) :: a2l
  type(land_data_type), intent(inout) :: bnd

  deallocate ( &
       bnd % tile_size  , & 
       bnd % t_surf     , &
       bnd % t_ca       , &
       bnd % albedo     , & 
       bnd % albedo_vis_dir , &
       bnd % albedo_nir_dir , &
       bnd % albedo_vis_dif , &
       bnd % albedo_nir_dif , &
       bnd % rough_mom  , & 
       bnd % rough_heat , & 
       bnd % rough_scale, &
       bnd % discharge      ,      & 
       bnd % discharge_heat ,      & 
       bnd % discharge_snow ,      &
       bnd % discharge_snow_heat , &
       bnd % mask        )
  if(associated( bnd%tr )) deallocate( bnd%tr )

  deallocate( a2l % t_flux   )
  deallocate( a2l % lw_flux  )
  deallocate( a2l % sw_flux  )
  deallocate( a2l % lprec    )
  deallocate( a2l % fprec    )
  deallocate( a2l % tprec    )
  deallocate( a2l % dhdt     )
  deallocate( a2l % drdt     )
  deallocate( a2l % drag_q   )
  deallocate( a2l % p_surf   )
  deallocate( a2l % sw_flux_down_vis_dir    )
  deallocate( a2l % sw_flux_down_total_dir  )
  deallocate( a2l % sw_flux_down_vis_dif    )
  deallocate( a2l % sw_flux_down_total_dif  )
  deallocate( a2l % tr_flux )
  deallocate( a2l % dfdtr )

end subroutine deallocate_boundary_data
! ============================================================================

end module land_model_mod
