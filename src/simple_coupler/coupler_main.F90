
program coupler_main

!-----------------------------------------------------------------------
!
!   program that couples component models for the atmosphere,
!   ocean (amip), land, and sea-ice using the exchange module
!
!-----------------------------------------------------------------------

use time_manager_mod, only: time_type, set_calendar_type, set_time,  &
                            set_date, days_in_month, month_name,     &
                            operator(+), operator (<), operator (>), &
                            operator (/=), operator (/), get_date,   &
                            operator (*), THIRTY_DAY_MONTHS, JULIAN, &
                            NOLEAP, NO_CALENDAR

use  atmos_model_mod, only: atmos_model_init, atmos_model_end, &
                            update_atmos_model_down,           &
                            update_atmos_model_up,             &
                            atmos_data_type,                   &
                            land_ice_atmos_boundary_type

use   land_model_mod, only: land_model_init, land_model_end, &
                            update_land_model_fast,          &
                            update_land_model_slow,          &
                            land_data_type,                  &
                            atmos_land_boundary_type

use    ice_model_mod, only: ice_model_init, ice_model_end,  &
                            update_ice_model_fast,          &
                            update_ice_model_slow,          &
                            ice_data_type,                  &
                            atmos_ice_boundary_type
                           !land_ice_boundary_type

use constants_mod, only:    constants_init
use       fms_mod, only: open_namelist_file, file_exist, check_nml_error,  &
                         error_mesg, fms_init, fms_end, close_file,        &
                         write_version_number, uppercase
use    fms_io_mod, only: fms_io_exit

use mpp_mod, only: mpp_init, mpp_pe, mpp_root_pe, &
                   stdlog, mpp_error, NOTE, FATAL, WARNING
use mpp_mod, only: mpp_clock_id, mpp_clock_begin, mpp_clock_end

use mpp_io_mod, only: mpp_open, mpp_close, &
                      MPP_NATIVE, MPP_RDONLY, MPP_DELETE

use mpp_domains_mod, only: mpp_get_global_domain, mpp_global_field, CORNER

use flux_exchange_mod, only: flux_exchange_init,   &
                             sfc_boundary_layer,   &
                             flux_down_from_atmos, &
                             flux_up_to_atmos,     &
                             flux_exchange_end       ! may not be used?

use  diag_manager_mod, only: diag_manager_init, diag_manager_end, &
                             get_base_date, diag_manager_set_time_end

use data_override_mod, only: data_override_init


implicit none

!-----------------------------------------------------------------------

character(len=128) :: version = '$Id: coupler_main.F90,v 19.0 2012/01/06 22:06:48 fms Exp $'
character(len=128) :: tag = '$Name: siena $'

!-----------------------------------------------------------------------
!---- model defined-types ----

 type (atmos_data_type) :: Atm
 type  (land_data_type) :: Land
 type   (ice_data_type) :: Ice

 type(atmos_land_boundary_type)     :: Atmos_land_boundary
 type(atmos_ice_boundary_type)      :: Atmos_ice_boundary
 type(land_ice_atmos_boundary_type) :: Land_ice_atmos_boundary

!-----------------------------------------------------------------------
!---- storage for fluxes ----

 real, allocatable, dimension(:,:)   ::                         &
    t_surf_atm, albedo_atm, land_frac_atm, dt_t_atm, dt_q_atm,  &
    flux_u_atm, flux_v_atm, dtaudv_atm, u_star_atm, b_star_atm, &
    rough_mom_atm

!-----------------------------------------------------------------------
! ----- coupled model time -----

   type (time_type) :: Time_atmos, Time_init, Time_end,  &
                       Time_step_atmos, Time_step_ocean
   integer :: num_cpld_calls, num_atmos_calls, nc, na

! ----- coupled model initial date -----

   integer :: date_init(6)
   integer :: calendar_type = -99

! ----- timing flags -----

   integer :: initClock, mainClock, termClock
   integer, parameter :: timing_level = 1

!-----------------------------------------------------------------------

      integer, dimension(6) :: current_date = (/ 0, 0, 0, 0, 0, 0 /)
      character(len=17) :: calendar = '                 '
      logical :: force_date_from_namelist = .false.  ! override restart values for date
      integer :: months=0, days=0, hours=0, minutes=0, seconds=0
      integer :: dt_atmos = 0
      integer :: dt_ocean = 0

      namelist /coupler_nml/ current_date, calendar, force_date_from_namelist, &
                             months, days, hours, minutes, seconds,  &
                             dt_atmos, dt_ocean

!#######################################################################

 call fms_init()
 call mpp_init()
 initClock = mpp_clock_id( 'Initialization' )
 mainClock = mpp_clock_id( 'Main loop' )
 termClock = mpp_clock_id( 'Termination' )
 call mpp_clock_begin (initClock)
  
 call fms_init
 call constants_init

 call coupler_init

 call mpp_clock_end (initClock) !end initialization
 call mpp_clock_begin(mainClock) !begin main loop


!------ ocean/slow-ice integration loop ------

 do nc = 1, num_cpld_calls

!------ atmos/fast-land/fast-ice integration loop -------

 do na = 1, num_atmos_calls

       Time_atmos = Time_atmos + Time_step_atmos

       call sfc_boundary_layer (real(dt_atmos), Time_atmos, Atm, Land, Ice, &
                                Land_ice_atmos_boundary                     )

       call update_atmos_model_down( Land_ice_atmos_boundary, Atm )

       call flux_down_from_atmos( Time_atmos, Atm, Land, Ice, &
                                  Land_ice_atmos_boundary,    &
                                  Atmos_land_boundary,        &
                                  Atmos_ice_boundary          )

     !--- land and ice models ---

      call update_land_model_fast ( Atmos_land_boundary, Land )
      call update_ice_model_fast  ( Atmos_ice_boundary,  Ice  )

     !--- atmosphere up ---

       call flux_up_to_atmos( Time_atmos, Land, Ice, Land_ice_atmos_boundary )

       call update_atmos_model_up( Land_ice_atmos_boundary, Atm )

 enddo

    !--- call land slow for diagnostics ---
      call update_land_model_slow ( Atmos_land_boundary, Land )

    ! need flux call to put runoff and p_surf on ice grid

    ! call flux_land_to_ice ( Time_atmos, Land, Ice, Land_ice_boundary )

    !----- slow-ice/ocean model ------

       call update_ice_model_slow ( Atmos_ice_boundary, Ice )

  enddo

!-----------------------------------------------------------------------

 call mpp_clock_end(mainClock)
 call mpp_clock_begin(termClock)

 call coupler_end
 call mpp_clock_end(termClock)

 call fms_end

!-----------------------------------------------------------------------

 stop

contains

!#######################################################################

   subroutine coupler_init

!-----------------------------------------------------------------------
!   initialize all defined exchange grids and all boundary maps
!-----------------------------------------------------------------------
    integer :: total_days, total_seconds, unit, ierr, io
    integer :: n, gnlon, gnlat
    integer :: date(6), flags
    type (time_type) :: Run_length
    character(len=9) :: month
    logical :: use_namelist
    
    logical, allocatable, dimension(:,:) :: mask
    real,    allocatable, dimension(:,:) :: glon_bnd, glat_bnd
!-----------------------------------------------------------------------
!----- initialization timing identifiers ----

!----- read namelist -------
!----- for backwards compatibilty read from file coupler.nml -----

   if (file_exist('input.nml')) then
      unit = open_namelist_file ()
   else
      call error_mesg ('program coupler',  &
                       'namelist file input.nml does not exist', FATAL)
   endif
   
      ierr=1
      do while (ierr /= 0)
          read  (unit, nml=coupler_nml, iostat=io, end=10)
          ierr = check_nml_error (io, 'coupler_nml')
      enddo
10    call close_file (unit)

!----- write namelist to logfile -----

   call write_version_number (version, tag)
   if (mpp_pe() == mpp_root_pe()) write(stdlog(),nml=coupler_nml)

!----- read restart file -----

   if (file_exist('INPUT/coupler.res')) then
       call mpp_open( unit, 'INPUT/coupler.res', action=MPP_RDONLY )
       read (unit,*,err=999) calendar_type
       read (unit,*) date_init
       read (unit,*) date
       goto 998 !back to fortran-4
     ! read old-style coupler.res
   999 call mpp_close (unit)
       call mpp_open (unit, 'INPUT/coupler.res', action=MPP_RDONLY, form=MPP_NATIVE)
       read (unit) calendar_type
       read (unit) date
   998 call mpp_close(unit)
   else
       force_date_from_namelist = .true.
   endif       

!----- use namelist value (either no restart or override flag on) ---

 if ( force_date_from_namelist ) then

    if ( sum(current_date) <= 0 ) then
         call error_mesg ('program coupler',  &
              'no namelist value for current_date', FATAL)
    else
         date      = current_date
    endif

!----- override calendar type with namelist value -----

        select case( uppercase(trim(calendar)) )
        case( 'JULIAN' )
            calendar_type = JULIAN
        case( 'NOLEAP' )
            calendar_type = NOLEAP
        case( 'THIRTY_DAY' )
            calendar_type = THIRTY_DAY_MONTHS
        case( 'NO_CALENDAR' )
            calendar_type = NO_CALENDAR
        case default
            call mpp_error ( FATAL, 'COUPLER_MAIN: coupler_nml entry calendar must '// &
                                    'be one of JULIAN|NOLEAP|THIRTY_DAY|NO_CALENDAR.' )
        end select

 endif

    call set_calendar_type (calendar_type)

!----- write current/initial date actually used to logfile file -----

    if ( mpp_pe() == mpp_root_pe() ) then
      write (stdlog(),16) date(1),trim(month_name(date(2))),date(3:6)
    endif

 16 format ('  current date used = ',i4,1x,a,2i3,2(':',i2.2),' gmt') 

!-----------------------------------------------------------------------
!------ initialize diagnostics manager ------

      call diag_manager_init

!----- always override initial/base date with diag_manager value -----

 call get_base_date ( date_init(1), date_init(2), date_init(3), &
                      date_init(4), date_init(5), date_init(6)  )

!----- use current date if no base date ------

    if ( date_init(1) == 0 ) date_init = date

!----- set initial and current time types ------

    Time_init  = set_date (date_init(1), date_init(2), date_init(3), &
                           date_init(4), date_init(5), date_init(6))

    Time_atmos = set_date (date(1), date(2), date(3),  &
                           date(4), date(5), date(6))

!-----------------------------------------------------------------------
!----- compute the ending time (compute days in each month first) -----
!
!   (NOTE: if run length in months then starting day must be <= 28)

    if ( months > 0 .and. date(3) > 28 )     &
        call error_mesg ('program coupler',  &
       'if run length in months then starting day must be <= 28', FATAL)

    Time_end = Time_atmos
    total_days = 0
    do n = 1, months
       total_days = total_days + days_in_month(Time_end)
       Time_end = Time_atmos + set_time (0,total_days)
    enddo

    total_days    = total_days + days
    total_seconds = hours*3600 + minutes*60 + seconds
    Run_length    = set_time (total_seconds,total_days)
    Time_end      = Time_atmos + Run_length

    !Need to pass Time_end into diag_manager for multiple thread case.
    call diag_manager_set_time_end(Time_end)


!-----------------------------------------------------------------------
!----- write time stamps (for start time and end time) ------

      call mpp_open( unit, 'time_stamp.out', nohdrs=.TRUE. )

      month = month_name(date(2))
      if ( mpp_pe() == mpp_root_pe() ) write (unit,20) date, month(1:3)

      call get_date (Time_end, date(1), date(2), date(3),  &
                               date(4), date(5), date(6))
      month = month_name(date(2))
      if ( mpp_pe() == mpp_root_pe() ) write (unit,20) date, month(1:3)

      call mpp_close (unit)

  20  format (6i4,2x,a3)

!-----------------------------------------------------------------------
!----- compute the time steps ------

Time_step_atmos = set_time (dt_atmos,0)
Time_step_ocean = set_time (dt_ocean,0)
num_cpld_calls  = Run_length / Time_step_ocean
num_atmos_calls = Time_step_ocean / Time_step_atmos

!-----------------------------------------------------------------------
!------------------- some error checks ---------------------------------

!----- initial time cannot be greater than current time -------

    if ( Time_init > Time_atmos ) call error_mesg ('program coupler',  &
                    'initial time is greater than current time', FATAL)

!----- make sure run length is a multiple of ocean time step ------

    if ( num_cpld_calls * Time_step_ocean /= Run_length )  &
         call error_mesg ('program coupler',  &
         'run length must be multiple of ocean time step', FATAL)

! ---- make sure cpld time step is a multiple of atmos time step ----

    if ( num_atmos_calls * Time_step_atmos /= Time_step_ocean )  &
         call error_mesg ('program coupler',   &
         'atmos time step is not a multiple of the ocean time step', FATAL)


!------ initialize component models ------

      call  atmos_model_init (Atm,  Time_init, Time_atmos, Time_step_atmos)

      call mpp_get_global_domain(Atm%Domain, xsize=gnlon, ysize=gnlat)
      allocate ( glon_bnd(gnlon+1,gnlat+1), glat_bnd(gnlon+1,gnlat+1) )
      call mpp_global_field(Atm%Domain, Atm%lon_bnd, glon_bnd, position=CORNER)
      call mpp_global_field(Atm%Domain, Atm%lat_bnd, glat_bnd, position=CORNER)

      call   land_model_init (Atmos_land_boundary, Land, &
                              Time_init, Time_atmos, Time_step_atmos, Time_step_ocean, &
                              glon_bnd, glat_bnd, atmos_domain=Atm%Domain)
      call    ice_model_init (Ice,  Time_init, Time_atmos, Time_step_atmos, Time_step_ocean, &
                              glon_bnd, glat_bnd, atmos_domain=Atm%Domain)

      call data_override_init ( ) ! Atm_domain_in  = Atm%domain, &
                                  ! Ice_domain_in  = Ice%domain, &
                                  ! Land_domain_in = Land%domain )

!------------------------------------------------------------------------
!---- setup allocatable storage for fluxes exchanged between models ----
!---- use local grids -----

   call flux_exchange_init (Time_atmos, Atm, Land, Ice, &
                     !!!!!  atmos_land_boundary,        &
                            atmos_ice_boundary,         &
                            land_ice_atmos_boundary     )




!-----------------------------------------------------------------------
!---- open and close dummy file in restart dir to check if dir exists --

    call mpp_open( unit, 'RESTART/file' )
    call mpp_close(unit, MPP_DELETE)

!-----------------------------------------------------------------------

   end subroutine coupler_init

!#######################################################################

   subroutine coupler_end

   integer :: unit, date(6)
!-----------------------------------------------------------------------

      call atmos_model_end (Atm)
      call  land_model_end (Atmos_land_boundary, Land)
      call   ice_model_end (Ice)

      call  fms_io_exit

! call flux_exchange_end (Atm)

!----- compute current date ------

      call get_date (Time_atmos, date(1), date(2), date(3),  &
                                 date(4), date(5), date(6))

!----- check time versus expected ending time ----

      if (Time_atmos /= Time_end) call error_mesg ('program coupler',  &
              'final time does not match expected ending time', WARNING)

!----- write restart file ------

    call mpp_open( unit, 'RESTART/coupler.res', nohdrs=.TRUE. )
    if (mpp_pe() == mpp_root_pe())then
        write( unit, '(i6,8x,a)' )calendar_type, &
             '(Calendar: no_calendar=0, thirty_day_months=1, julian=2, gregorian=3, noleap=4)'

        write( unit, '(6i6,8x,a)' )date_init, &
             'Model start time:   year, month, day, hour, minute, second'
        write( unit, '(6i6,8x,a)' )date, &
             'Current model time: year, month, day, hour, minute, second'
    endif
    call mpp_close(unit)


!----- final output of diagnostic fields ----

   call diag_manager_end (Time_atmos)

!-----------------------------------------------------------------------

   end subroutine coupler_end

!#######################################################################

end program coupler_main

