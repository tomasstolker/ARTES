program artes

  use omp_lib

  implicit none

  ! Constants
  integer,  parameter       :: dp    = kind(0.d0)                     ! Double precision
  real(dp), parameter       :: pi    = 4._dp*datan(1._dp)             ! Surface area of a disk with radius unity
  real(dp), parameter       :: k_b   = 1.3806488e-23_dp               ! Boltzmann constant [m2 kg s-2 K-1]
  real(dp), parameter       :: sb    = 5.670373e-8_dp                 ! Stefan-Boltzmann constant [J s-1 m-2 K-4]
  real(dp), parameter       :: hh    = 6.62606957e-34_dp              ! Planck constant [m2 kg s-1]
  real(dp), parameter       :: cc    = 2.99792458e8_dp                ! Speed of light [m s-1]
  real(dp), parameter       :: r_sun = 6.95500e8_dp                   ! Solar radius [m]
  real(dp), parameter       :: pc    = 3.08572e16_dp                  ! Parsec [m]
  real(dp), parameter       :: au    = 1.49598e11_dp                  ! Astronomical Unit [m]

  ! Input files
  logical                   :: log_file                               ! Log file on or output to the screen
  integer                   :: photon_source                          ! Photon emission point (1=star, 2=planet)
  integer                   :: det_type                               ! 1. image 2. spectrum 3. phase curve
  real(dp)                  :: wavelength                             ! Photon wavelength
  integer(16)               :: packages                               ! Number of photon packages
  real(dp)                  :: orbit                                  ! Star-planet distance [m]
  real(dp)                  :: distance_planet                        ! Distance from the observer to the planet [m]
  real(dp)                  :: fstop                                  ! Photon stop parameter [0:1]
  real(dp)                  :: photon_minimum                         ! Minimum photon energy possible (0<=photon_minimum<=1), otherwise removed
  logical                   :: thermal_weight                         ! Weight cell luminosity
  logical                   :: photon_scattering                      ! Scattering on or off
  integer                   :: photon_emission                        ! Isotropic (1) or biased (2) emission
  real(dp)                  :: photon_bias                            ! Bias parameter (0 <= bias < 1)
  real(dp)                  :: photon_walk                            ! Modified random walk parameter (0 <= walk)
  integer                   :: photon_number                          ! Photon number for debug output
  real(dp)                  :: t_star                                 ! Stellar effective temperature [K]
  real(dp)                  :: r_star                                 ! Stellar radius [m]
  real(dp)                  :: star_theta                             ! Stellar location, theta direction
  real(dp)                  :: star_phi                               ! Stellar location, phi direction
  real(dp)                  :: sin_star_theta, cos_star_theta         ! Sine and cosine of star theta
  real(dp)                  :: sin_star_phi, cos_star_phi             ! Sine and cosine of star phi
  real(dp)                  :: surface_albedo                         ! Surface albedo [0,1] (0=black, 1=Lambertian)
  real(dp)                  :: oblateness                             ! Planet oblateness (0-1)
  real(dp)                  :: oblate_x, oblate_y, oblate_z           ! Planet aspect_ratio (=R_equator/R_poles)
  integer                   :: nr                                     ! Number of grid points in r-direction
  integer                   :: ntheta                                 ! Number of grid points in theta-direction
  integer                   :: nphi                                   ! Number of grid points in phi-direction
  real(dp)                  :: det_theta, det_phi                     ! Detector direction [rad]
  real(dp)                  :: det_angle                              ! Detector position angle [rad]
  real(dp)                  :: sin_det_angle, cos_det_angle           ! Sine and cosine of detector position angle
  real(dp)                  :: theta_phase, phi_phase                 ! Theta and phi coordinate of phase curve detector
  integer                   :: npix                                   ! Number of pixels in x and y direction
  real(dp)                  :: x_max                                  ! Image plane size ranges from [-x_max:x_max]
  real(dp)                  :: y_max                                  ! Image plane size ranges from [-y_max:y_max]
  character(100)            :: email                                  ! Send email when job is finished
  logical                   :: global                                 ! Global energy transport output (on/off)
  logical                   :: latitudinal                            ! Latitudinal energy transport output (on/off)
  logical                   :: ring                                   ! Ring system (on/off)
  real(dp)                  :: grid_min                               ! Minimum fractional distance to cell boundary
  logical                   :: debug                                  ! Output debug errors (on/off)
  real(dp)                  :: tau_max                                ! Maximum vertical optical depth
  logical                   :: ranseed                                ! Random seed for random number generator (on/off)

  ! Arrays
  real(dp), allocatable     :: rfront(:)                              ! Radial front point of the grid cell
  real(dp), allocatable     :: thetafront(:)                          ! Theta front point of the grid cell
  integer,  allocatable     :: thetaplane(:)                          ! Type of theta surface (1=normal/cone, 2=flat/xy-plane)
  real(dp), allocatable     :: phifront(:)                            ! Phi front point of the grid cell
  real(dp), allocatable     :: cell_density(:,:,:)                    ! Cell density [kg/m3]
  real(dp), allocatable     :: cell_temperature(:,:,:)                ! Cell temperature [K]
  real(dp), allocatable     :: cell_scattering_opacity(:,:,:,:)       ! Cell scattering opacity [m-1]
  real(dp), allocatable     :: cell_absorption_opacity(:,:,:,:)       ! Cell absorption opacity [m-1]
  real(dp), allocatable     :: cell_opacity(:,:,:,:)                  ! Cell total opacity [m-1]
  real(dp), allocatable     :: cell_asymmetry(:,:,:,:)                ! Cell asymmetry parameter
  real(dp), allocatable     :: cell_luminosity(:,:,:)                 ! Cell thermal luminosity [W]
  real(dp), allocatable     :: cell_weight(:,:,:)                     ! Cell luminosity weight
  real(dp), allocatable     :: cell_scatter_matrix(:,:,:,:,:,:)       ! Cell scatter matrix number
  real(dp), allocatable     :: cell_albedo(:,:,:,:)                   ! Cell albedo
  real(dp), allocatable     :: cell_volume(:,:,:)                     ! Cell volume [m3]
  real(dp), allocatable     :: cell_p11_int(:,:,:,:)                  ! Integral of the P11 scattering matrix element
  real(dp), allocatable     :: cell_p12_int(:,:,:,:)                  ! Integral of the P12 scattering matrix element
  real(dp), allocatable     :: cell_p13_int(:,:,:,:)                  ! Integral of the P13 scattering matrix element
  real(dp), allocatable     :: cell_p14_int(:,:,:,:)                  ! Integral of the P14 scattering matrix element
  real(dp), allocatable     :: theta_grid_cos(:)                      ! Cosine of theta boundaries
  real(dp), allocatable     :: theta_grid_sin(:)                      ! Sine of theta boundaries
  real(dp), allocatable     :: theta_grid_tan(:)                      ! Tangent of theta boundaries
  real(dp), allocatable     :: phi_grid_sin(:)                        ! Sine of phi boundaries
  real(dp), allocatable     :: phi_grid_cos(:)                        ! Cosine of phi boundaries
  real(dp), allocatable     :: wavelengths(:)                         ! Wavelength points [micron]
  real(dp), allocatable     :: cell_flow_global(:,:,:,:,:)            ! Global energy transport through each cell
  real(dp), allocatable     :: cell_flow(:,:,:,:,:)                   ! Energy transport through cell boundary (1=upward, 2=downward, 3=south, 4=north)
  real(dp), allocatable     :: emissivity_cumulative(:,:,:)           ! Cumulative emissivity in each grid cell CDF
  real(dp), allocatable     :: detector(:,:,:,:)                      ! Array(flux/flux^2/photons, Stokes I/Q/U/V, npix, npix)
  real(dp), allocatable     :: detector_thread(:,:,:,:,:)             ! Array(thread_id, flux/flux^2/photons, Stokes I/Q/U/V, npix, npix)
  real(dp), allocatable     :: flux_emitted(:)                        ! Thermal flux emitted in the atmosphere
  real(dp), allocatable     :: flux_exit(:)                           ! Flux exiting from the atmosphere
  real(dp)                  :: sinbeta(360)                           ! sin(beta)
  real(dp)                  :: sin2beta(360), cos2beta(360)           ! sin(2*beta), cos(2*beta)
  real(dp)                  :: photometry(11)                         ! Detector integrated values
  real(dp)                  :: det_dir(5)                             ! Detector direction

  ! Global variables
  real(dp)                  :: phase_observer                         ! Phase angle of observation [deg]
  real(dp)                  :: sin_det_theta, cos_det_theta
  real(dp)                  :: sin_det_phi, cos_det_phi
  integer                   :: out_unit                               ! Output unit, 10=logfile, 6=screen
  integer                   :: cores                                  ! Number of computer cores
  integer                   :: threads                                ! Number of threads in use
  integer                   :: cells                                  ! Number of grid cells
  integer                   :: cell_depth                             ! Maximum cell depth to emit thermal photons
  real(dp)                  :: t_start, t_end
  character(100)            :: output_name
  character(100)            :: input_file                             ! Input parameter file
  real(dp)                  :: package_energy                         ! Photon package energy [J]
  real(dp)                  :: x_fov, y_fov                           ! Field of view [mas]
  real(dp)                  :: pixel_scale                            ! Pixel size [mas pixel-1]
  integer                   :: n_wavelength                           ! Number of wavelengths
  integer                   :: wl_count                               ! Wavelength point
  character(300)            :: error_log                              ! Filename of error log
  integer, parameter        :: ns = 4                                 ! Random number generator
  integer, parameter        :: default_seed(ns) = (/ 521288629, 362436069, 16163801, 1131199299 /)
  integer, allocatable      :: state(:,:)                             ! Random number generator
  real(dp), allocatable     :: y_mrw(:), phi_mrw(:)                   ! Modified random walk
  integer                   :: ny_mrw                                 ! Modified random walk

  call run

contains

  subroutine run

    integer  :: i
    
    call initialize
    call grid_initialize(1)
    call output(1)

    wl_count = 1

    if (det_type.eq.2) then

       do

          wavelength = wavelengths(wl_count)

          call array_start
          call grid_initialize(2)

          if (wl_count.eq.1) call output(2)

          write (6,fmt="(a1,a,t13,f7.3,a)", advance="no") achar(13), &
               "Wavelength: ", wavelength*1.e6_dp, " micron"

          call radiative_transfer
          call write_output

          call grid_finished(1)
          call grid_finished(2)

          if (wl_count.eq.size(wavelengths)) then

             write (6,fmt="(a1,a,t13,f6.3,a)", advance="yes") achar(13), ""
             if (.not.log_file) write (6,'(a)') ""
             call output(3)
             exit

          end if

          wl_count = wl_count + 1

       end do

       call grid_finished(3)
       call output(4)

    else if (det_type.eq.1.or.det_type.eq.3) then

       wavelength = wavelengths(wl_count)
       
       call grid_initialize(2)
       call output(2)

       if (det_type.eq.3) then

          do i = 1,74

             if (i.eq.1) then
                det_phi = 1.e-5_dp*pi/180._dp
             else if (i.eq.2) then
                det_phi = 2.5_dp*pi/180._dp
             else if (i.eq.73) then
                det_phi = (180._dp-1e-5_dp)*pi/180._dp
             else if (i.eq.74) then
                write (6,fmt="(a1,a,t14,f6.1,a)", advance="yes") achar(13), ""
                if (.not.log_file) write (6,'(a)') ""
                exit
             else
                det_phi = det_phi + 2.5_dp*pi/180._dp
             end if

             det_dir(5) = det_phi
             cos_det_phi = cos(det_dir(5))
             sin_det_phi = sqrt(1._dp-cos_det_phi*cos_det_phi)

             call spherical_cartesian(1._dp, det_theta, det_phi, det_dir(1), det_dir(2), det_dir(3))

             write (6,fmt="(a1,a,t14,f6.1,a)", advance="no") achar(13), &
                  "Phase angle: ", det_phi*180._dp/pi, " deg"

             call array_start
             call radiative_transfer
             call write_output
             call grid_finished(2)

          end do

          call output(3)
          call grid_finished(1)
          call grid_finished(3)
          call output(4)

       else if (det_type.eq.1) then

          call array_start
          call radiative_transfer
          call write_output
          call output(3)
          call grid_finished(1)
          call grid_finished(2)
          call grid_finished(3)
          call output(4)

       end if

    end if

  end subroutine run

  subroutine initialize

    integer                     :: i, j, nlines, nmax_mrw
    logical                     :: file_exists
    character(1)                :: first_char
    character(100)              :: key_word, key_value
    character(500), allocatable :: data(:)
    real(dp)                    :: xi

    ! ARTES input file
    
    log_file = .false.
    email = ""
    ranseed = .true.

    photon_source = 1
    packages = 100000
    fstop = 1.e-5_dp
    photon_minimum = 1.e-20_dp
    thermal_weight = .false.
    photon_scattering = .true.
    photon_emission = 1
    photon_bias = 0.8_dp
    photon_walk = -1._dp
    photon_number = 0

    t_star = 5800._dp
    r_star = r_sun
    star_theta = pi/2._dp
    star_phi = 0._dp

    surface_albedo = 0._dp
    oblateness = 0._dp
    orbit = 5._dp*au
    ring = .false.
    grid_min = 1.e-6_dp
    tau_max = -1._dp

    det_type = 1
    det_theta = pi/2._dp
    det_phi = pi/2._dp
    npix = 25
    det_angle = -1._dp
    distance_planet = 10._dp*pc

    debug = .false.
    global = .false.
    latitudinal = .false.

    ! Other parameters
    
    det_dir = 0._dp
    phase_observer = 0._dp
    input_file = ""
    package_energy = 0._dp
    i = 0
    j = 0
    nlines = 0
    first_char = ""
    file_exists = .false.
    wavelength = 0._dp
    x_max = 0._dp
    y_max = 0._dp
    oblate_x = 1._dp
    oblate_y = 1._dp
    oblate_z = 1._dp
    cell_depth = 0

    inquire(file="/proc/cpuinfo", exist=file_exists)

    if (file_exists) then

       ! Check number of cores on a Linux computer

       call system("grep -c ^processor /proc/cpuinfo > cores.txt")
       open (100, file='cores.txt')
       read (100,*) cores
       close (100, status='delete')

    else if (.not.file_exists) then

       ! Check number of cores on a Apple computer

       call system("sysctl hw.ncpu | awk '{print $2}' > cores.txt")
       open (100, file='cores.txt')
       read (100,*) cores
       close (100, status='delete')

    end if
    
    ! Number of threads in use
    !$omp parallel
    threads = omp_get_num_threads()
    !$omp end parallel

    ! CPU start time
    t_start = omp_get_wtime()

    ! Get input directory and number of photons
    call argument_input(input_file)

    ! Check if input file exists
    inquire (file=input_file, exist=file_exists)
    if (.not.file_exists) then
       write (6,'(a)') "Input file does not exist!"
       call exit(0)
    end if

    ! Read the input file
    call readfile(input_file, data, nlines)

    ! Get the keywords and values, ignore comment lines
    do i=1,nlines

       ! Get first character of string
       first_char = data(i)(1:1)

       ! Check if the line is not a comment line
       if (first_char.ne."*".and.first_char.ne."-".and.first_char.ne."=".and.len_trim(data(i)).ne.0) then

          call get_key_value(data(i), key_word, key_value)
          call input_parameters(key_word, key_value)

       end if

    enddo

    ! Get command line keywords
    call argument_keywords
    
    ! Sine and cosine arrays

    sinbeta  = 0._dp
    cos2beta = 0._dp
    sin2beta = 0._dp

    do i=1,180

       sinbeta(i)      = ( sin(dble(i)*pi/180._dp) + sin(dble(i-1)*pi/180._dp) ) / 2._dp
       sinbeta(i+180)  = -sinbeta(i)
       cos2beta(i)     = ( cos(2._dp*dble(i)*pi/180._dp) + cos(2._dp*dble(i-1)*pi/180._dp) ) / 2._dp
       cos2beta(i+180) = cos2beta(i)
       sin2beta(i)     = ( sin(2._dp*dble(i)*pi/180._dp) + sin(2._dp*dble(i-1)*pi/180._dp) ) / 2._dp
       sin2beta(i+180) = sin2beta(i)

    end do

    ! Get atmospheric structure
    call get_atmosphere

    ! Open error log
    if (debug) then
       error_log = trim(output_name)//'/error.log'
       open (11, file=trim(error_log), status="new")
       close (11)
    end if

    ! Open output log
    if (log_file) open (10, file=trim(output_name)//'/output.log', status='new')

    ! Random number generator

    call init_random_seed

    allocate (state(threads,ns))
    state = 0
	
    do j=1,ns
       state(:,j) = default_seed(j)
    end do

    if (ranseed) then
       do i=1,threads
          call random_number(xi)
          state(i,1) = int(xi*1.e6_dp)
       end do
    end if

    ! Detector

    if (det_type.eq.2) then
       npix = 1
    else if (det_type.eq.3) then
       npix = 1
       det_theta = pi/2._dp
       det_phi = 1.e-5_dp
    end if
    
    ! Oblateness

    oblate_x = 1._dp / (1._dp-oblateness)
    oblate_y = oblate_x
    oblate_z = 1._dp

    ! Detector field of view

    x_max = 1.3_dp*rfront(nr)
    y_max = 1.3_dp*rfront(nr)

    x_max = ( oblateness + 1._dp ) * x_max
    y_max = ( oblateness + 1._dp ) * y_max

    ! Field of view [mas]
 
    x_fov = 2._dp*atan(x_max/distance_planet)*3600._dp*180._dp/pi*1000._dp
    y_fov = 2._dp*atan(y_max/distance_planet)*3600._dp*180._dp/pi*1000._dp
    
    ! Pixel size [mas pixel-1]
    
    pixel_scale = x_fov/npix

    ! Detector

    if (abs(det_phi).lt.1.e-3_dp.or.det_phi.gt.2._dp*pi-1.e-3_dp) det_phi = 1.e-3_dp
    if (det_phi.gt.pi-1.e-3_dp.and.det_phi.lt.pi+1.e-3_dp) det_phi = pi-1.e-3_dp

    call spherical_cartesian(1._dp, det_theta, det_phi, det_dir(1), det_dir(2), det_dir(3))
    det_dir(4) = det_theta
    det_dir(5) = det_phi

    cos_det_theta = cos(det_dir(4))
    sin_det_theta = sqrt(1._dp-cos_det_theta*cos_det_theta)

    cos_det_phi = cos(det_dir(5))
    sin_det_phi = sqrt(1._dp-cos_det_phi*cos_det_phi)
    if (det_dir(5).gt.pi.and.det_dir(5).lt.2._dp*pi) sin_det_phi = -sin_det_phi
    
    if (det_angle.gt.0._dp) then
       cos_det_angle = cos(det_angle)
       sin_det_angle = sqrt(1._dp-cos_det_angle*cos_det_angle)
       if (det_angle.gt.pi.and.det_angle.lt.2._dp*pi) sin_det_angle = -sin_det_angle
    else
       sin_det_angle = 0._dp
       cos_det_angle = 0._dp
    end if

    ! Star

    cos_star_theta = cos(star_theta)
    sin_star_theta = sqrt(1._dp-cos_star_theta*cos_star_theta)

    cos_star_phi = cos(star_phi)
    sin_star_phi = sqrt(1._dp-cos_star_phi*cos_star_phi)
    if (star_phi.gt.pi) sin_star_phi = -sin_star_phi

    ! Observer phase angle

    if (det_type.eq.1.or.det_type.eq.2) then

       phase_observer = sin_star_theta*cos_star_phi*sin_det_theta*cos_det_phi + &
            sin_star_theta*sin_star_phi*sin_det_theta*sin_det_phi + &
            cos_star_theta*cos_det_theta

       phase_observer = acos(phase_observer) * 180._dp/pi

    end if

    ! Modified random walk
    
    nmax_mrw = 1000
    ny_mrw = 1000

    allocate (phi_mrw(ny_mrw))
    phi_mrw = 0._dp

    allocate (y_mrw(ny_mrw))
    y_mrw = 0._dp
	
    y_mrw(1) = 1._dp
    phi_mrw(1) = 1._dp

    do i=2,ny_mrw

       y_mrw(i) = dble(ny_mrw-i+1) / dble(ny_mrw-1) * 0.75_dp
       
       phi_mrw(i) = 0._dp
       do j=1,nmax_mrw
          phi_mrw(i) = phi_mrw(i) + (-1._dp)**(j+1)*y_mrw(i)**(j**2)
       enddo
       phi_mrw(i) = 2._dp*phi_mrw(i)

    enddo
    
  end subroutine initialize

  subroutine radiative_transfer

    integer(16) :: i
    integer     :: j, k, l, thread_id, cell(3), cell_out(3), face(2), next_face(2)
    integer     :: face_check(2), cell_check(3), buffer(13), status, iy
    logical     :: grid_exit, surface, check1, check2, check3, check4, cell_error, walk
    real(dp)    :: face_distance, tau, tau_cell, tau_run, xi, x_out, y_out, z_out, tau_first
    real(dp)    :: stokes(4), stokes_new(4), alpha, beta, scatter(16), gamma, c2b
    real(dp)    :: direction(3), direction_new(3), x, y, z, s, x_check, y_check, z_check
    real(dp)    :: r_disk, phi_disk, dpi, dummy, bias_weight, d_mrw, r_min, acos_alpha
    
    check1 = .false.
    check2 = .false.
    check3 = .false.
    check4 = .false.

    !$omp parallel do default(none) private(surface) &

    !$omp& private(x, y, z, s, x_out, y_out, z_out, stokes, stokes_new, direction, direction_new, d_mrw, iy, r_min) &
    !$omp& private(face, next_face, face_distance, alpha, beta, acos_alpha, cell, cell_out, tau, tau_run, tau_first) &
    !$omp& private(scatter, xi, theta_phase, phi_phase, grid_exit, tau_cell, cell_error, thread_id, bias_weight) &
    !$omp& private(x_check, y_check, z_check, face_check, cell_check, buffer, status, r_disk, phi_disk, walk, c2b) &

    !$omp& shared(cell_opacity, cell_absorption_opacity, cell_albedo, nr, gamma, photon_source, log_file, photon_minimum) &
    !$omp& shared(surface_albedo, packages, threads, cell_depth, fstop, det_type, y_mrw, phi_mrw, photon_walk, ny_mrw) &
    !$omp& shared(check1, check2, check3, check4, det_phi, det_theta, output_name, wl_count, rfront, cell_scattering_opacity) &
    !$omp& shared(global, latitudinal, package_energy, cell_weight, flux_emitted, flux_exit, error_log, photon_scattering) &
    !$omp& shared(cell_asymmetry, debug)
    
    do i=1,packages

       if (debug) then
          
          ! Check if error file is not becoming larger than 100 MB
          call stat(trim(output_name)//'/error.log', buffer, status)
          if (buffer(8).gt.1.e8_dp) then
             open (11, file=trim(error_log), position="append")
             write (11,*) "error 001"
             close (11)
             call exit(0)
          end if

       end if
          
       surface = .false.
       grid_exit = .false.
       cell_error = .false.
       walk = .false.
       cell = 0
       cell_out = 0
       face = 0
       next_face = 0
       x = 0._dp
       y = 0._dp
       z = 0._dp
       bias_weight = 0._dp

       thread_id = OMP_GET_THREAD_NUM()

       if (thread_id.eq.0.and.det_type.eq.1) then

          if (dble(i)/(dble(packages)/dble(threads)).gt.0.2_dp.and..not.check1) then
             write (6,fmt="(a3)", advance="no") "20%"
             check1 = .true.
          else if (dble(i)/(dble(packages)/dble(threads)).gt.0.4_dp.and..not.check2) then
             write (6,fmt="(a9)", advance="no") "  --  40%"
             check2 = .true.
          else if (dble(i)/(dble(packages)/dble(threads)).gt.0.6_dp.and..not.check3) then
             write (6,fmt="(a9)", advance="no") "  --  60%"
             check3 = .true.
          else if (dble(i)/(dble(packages)/dble(threads)).gt.0.8_dp.and..not.check4) then
             write (6,fmt="(a9)", advance="no") "  --  80%"
             check4 = .true.
          else if (i.eq.packages/threads) then
             write (6,fmt='(a10)', advance="yes") "  --  100%"
             if (.not.log_file) write (6,'(a)') ""
          end if

       end if

       call emit_photon(thread_id, x, y, z, direction, face, cell, r_disk, phi_disk, bias_weight)

       stokes(1) = 1._dp
       stokes(2) = 0._dp
       stokes(3) = 0._dp
       stokes(4) = 0._dp

       if (photon_source.eq.2) then

          ! Correction of the photon energy for the cell emission probability
          ! stokes(1) = stokes(1) * exp( cell_absorption_opacity(cell(1),cell(2),cell(3),wl_count) * &
          !     ( rfront(cell(1)+1) - rfront(cell(1)) ) )

          stokes(1) = stokes(1) * bias_weight / cell_weight(cell(1),cell(2),cell(3))

          flux_emitted(thread_id+1) = flux_emitted(thread_id+1) + stokes(1)

          call peel_thermal(thread_id, i, x, y, z, stokes, cell, face, cell_error)

          if (cell_error) then
             if (debug) then
                open (11, file=trim(error_log), position="append")
                write (11,*) "error 047"
                close (11)
             end if
             cycle
          end if
          
       end if

       ! Get optical depth to grid boundary or surface

       tau_first = 0._dp

       x_check = x
       y_check = y
       z_check = z
       face_check = face
       cell_check = cell

       do

          call cell_face(i, x_check, y_check, z_check, direction, face_check, &
               next_face, face_distance, grid_exit, cell_check, cell_out, cell_error)

          tau_cell = face_distance*cell_opacity(cell_check(1), cell_check(2), cell_check(3), wl_count)
          tau_first = tau_first + tau_cell

          x_check = x_check + face_distance*direction(1)
          y_check = y_check + face_distance*direction(2)
          z_check = z_check + face_distance*direction(3)

          if (cell_error.or.grid_exit.or.(next_face(1).eq.1.and.next_face(2).eq.cell_depth)) exit

          face_check = next_face
          cell_check = cell_out

       end do

       if (cell_error) then
          if (debug) then
             open (11, file=trim(error_log), position="append")
             write (11,*) "error 071"
             close (11)
          end if
          cycle
       end if

       ! First optical depth

       if (next_face(1).eq.1.and.next_face(2).eq.cell_depth.and.surface_albedo.gt.0._dp) then

          ! Zero optical depth, photon hits planet surface
          call random(thread_id,xi)
          tau = -log(1._dp-xi)

       else

          if (tau_first.ge.1.e-6_dp) then

             ! Sample optical depth, force first interaction, weight photon
             call random(thread_id,xi)
             if (tau_first.lt.50._dp) then
                tau = -log(1._dp-xi*(1._dp-exp(-tau_first)))
                stokes = stokes * (1._dp-exp(-tau_first))
             else
                tau = -log(1._dp-xi)
             end if
          
          else if (tau_first.lt.1.e-6_dp) then

             ! New photon in case optical depth is zero and not crossing the surface
             cycle

          end if


       end if

       ! Do loop to start crossing cells till next interaction point

       tau_run = 0._dp

       do

          call cell_face(i, x, y, z, direction, face, next_face, face_distance, grid_exit, cell, cell_out, cell_error)

          if (cell_error) exit
          
          tau_cell = face_distance*cell_opacity(cell(1), cell(2), cell(3), wl_count)

          ! Check if next interaction happens in the upcoming cell

          if (tau_run+tau_cell.gt.tau) then
             
             grid_exit = .false.

             s = (tau-tau_run) / cell_opacity(cell(1), cell(2), cell(3), wl_count)

             x = x + s*direction(1)
             y = y + s*direction(2)
             z = z + s*direction(3)

             if (global) call add_flow_global(thread_id, x, y, z, direction, stokes(1), s, cell)
             
             face = 0
             next_face = 0

             exit

          else

             x = x + face_distance*direction(1)
             y = y + face_distance*direction(2)
             z = z + face_distance*direction(3)

             if (global) call add_flow_global(thread_id, x, y, z, direction, stokes(1), face_distance, cell)

             if (latitudinal) then
                if (next_face(1).eq.1) then
                   if (cell_out(1).gt.cell(1)) then
                      call add_flow(thread_id, 1, stokes(1), cell)
                   else if (cell_out(1).lt.cell(1)) then
                      call add_flow(thread_id, 2, stokes(1), cell)
                   end if
                else if (next_face(1).eq.2) then
                   if (cell_out(2).gt.cell(2)) then
                      call add_flow(thread_id, 3, stokes(1), cell)
                   else if (cell_out(2).lt.cell(2)) then
                      call add_flow(thread_id, 4, stokes(1), cell)
                   end if
                end if
             end if
             
             face = next_face
             cell = cell_out

             if (grid_exit) exit

          end if

          ! Photon is crossing the planet surface

          if (next_face(1).eq.1.and.next_face(2).eq.cell_depth) then

             call random(thread_id,xi)

             if (xi.gt.surface_albedo) then

                surface = .true.
                exit

             else if (xi.le.surface_albedo) then

                call lambertian(thread_id, x, y, z, stokes, direction)
                call peel_surface(thread_id, i, x, y, z, stokes, cell, face)

                cell(1) = cell(1)+1

                if (cell_error) exit

             end if

          end if

          tau_run = tau_run + tau_cell

       end do

       if (grid_exit) flux_exit(thread_id+1) = flux_exit(thread_id+1) + stokes(1)
       
       ! Emit new photon when current photon leaves the atmosphere or is absorbed by the planet surface

       if (surface.or.grid_exit.or.cell_error) cycle

       ! Let the photon scatter in the atmosphere until it either exits or is absorbed
       
       do

          ! No scattering
          if (.not.photon_scattering) exit
          
          call random(thread_id,xi)

          if (xi.lt.fstop) then

             ! Absorption
             exit

          else

             if (cell_albedo(cell(1),cell(2),cell(3),wl_count).lt.1._dp.and. &
                  cell_albedo(cell(1),cell(2),cell(3),wl_count).gt.0._dp) then

                ! Weighted scattering
                gamma = cell_albedo(cell(1),cell(2),cell(3),wl_count)/(1._dp-fstop)
                stokes = gamma*stokes

             end if

             ! Remove photon when Stokes I becomes too small
             if (stokes(1).le.photon_minimum) then
                cell_error = .true.
                exit
             end if

             call peel_photon(thread_id, i, x, y, z, stokes, direction, cell, face, cell_error)

             if (cell_error) exit

             call scattering_sampling(thread_id, stokes, alpha, beta, acos_alpha, c2b, cell)
             call get_scatter(cell, scatter, acos_alpha)
             call direction_cosine(alpha, beta, c2b, direction, direction_new)
             call stokes_rotation(alpha, beta, c2b, stokes, scatter, direction, direction_new, stokes_new, .false.)

             stokes = stokes_new
             direction = direction_new

          end if

          if (photon_walk.gt.0._dp.and.face(1).eq.0) then

             ! Minimum distance to cell face
             call min_distance(x, y, z, cell, r_min)

             if (r_min.gt.1e100_dp.or.r_min.le.0._dp) then
                if (debug) then
                   open (11, file=trim(error_log), position="append")
                   write (11,*) "error 057"
                   write (11,*) r_min
                   close (11)
                end if
                r_min = 0._dp
             end if

             if (r_min.gt.photon_walk/cell_opacity(cell(1), cell(2), cell(3), wl_count)) then

                ! Modified random walk

                walk = .true.
                
                iy = 1

                call random(thread_id, xi)
                call hunt(phi_mrw, xi, iy)

                d_mrw = ( 1._dp - cell_asymmetry(cell(1),cell(2),cell(3),wl_count) ) * &
                     cell_scattering_opacity(cell(1), cell(2), cell(3), wl_count) + &
                     cell_absorption_opacity(cell(1), cell(2), cell(3), wl_count)

                d_mrw = -3._dp * d_mrw * log(y_mrw(iy)) * r_min**2 / (pi**2)
                
                call random_direction(thread_id, direction)
                x = x + r_min*direction(1)
                y = y + r_min*direction(2)
                z = z + r_min*direction(3)

                stokes = stokes * exp(-cell_absorption_opacity(cell(1), cell(2), cell(3), wl_count)*d_mrw)

                if (global) call add_flow_global(thread_id, x, y, z, direction, stokes(1), r_min, cell)

                if (stokes(1).le.photon_minimum) exit

             else

                walk = .false.

             end if

          end if

          if (.not.walk) then

             ! Sample optical depth
             call random(thread_id,xi)
             tau = -log(1._dp-xi)

             tau_run = 0._dp

             do

                call cell_face(i, x, y, z, direction, face, next_face, face_distance, grid_exit, cell, cell_out, cell_error)

                if (cell_error) exit
                
                tau_cell = face_distance*cell_opacity(cell(1), cell(2), cell(3), wl_count)

                if (tau_run+tau_cell.gt.tau) then

                   ! Next interaction point is in the same grid cell
                   
                   grid_exit = .false.

                   s = (tau-tau_run) / cell_opacity(cell(1), cell(2), cell(3), wl_count)

                   x = x + s*direction(1)
                   y = y + s*direction(2)
                   z = z + s*direction(3)

                   if (global) call add_flow_global(thread_id, x, y, z, direction, stokes(1), s, cell)
                
                   face = 0
                   next_face = 0

                   exit

                else

                   ! Photon crosses a cell face before next interaction point

                   x = x + face_distance*direction(1)
                   y = y + face_distance*direction(2)
                   z = z + face_distance*direction(3)

                   if (global) call add_flow_global(thread_id, x, y, z, direction, stokes(1), face_distance, cell)

                   if (latitudinal) then
                      if (next_face(1).eq.1) then
                         if (cell_out(1).gt.cell(1)) then
                            call add_flow(thread_id, 1, stokes(1), cell)
                         else if (cell_out(1).lt.cell(1)) then
                            call add_flow(thread_id, 2, stokes(1), cell)
                         end if
                      else if (next_face(1).eq.2) then
                         if (cell_out(2).gt.cell(2)) then
                            call add_flow(thread_id, 3, stokes(1), cell)
                         else if (cell_out(2).lt.cell(2)) then
                            call add_flow(thread_id, 4, stokes(1), cell)
                         end if
                      end if
                   end if
                
                   face = next_face
                   cell = cell_out

                   if (grid_exit) exit

                end if

                ! Photon crosses the surface

                if (next_face(1).eq.1.and.next_face(2).eq.cell_depth) then

                   call random(thread_id,xi)

                   if (xi.gt.surface_albedo) then
                      
                      surface = .true.
                      exit

                   else if (xi.le.surface_albedo) then

                      call lambertian(thread_id, x, y, z, stokes, direction)
                      call peel_surface(thread_id, i, x, y, z, stokes, cell, face)

                      cell(1) = cell(1)+1

                      if (cell_error) exit

                   end if

                end if

                tau_run = tau_run + tau_cell
                face = next_face

             end do

             if (grid_exit.or.surface.or.cell_error) exit

          end if

       end do

       if (grid_exit) flux_exit(thread_id+1) = flux_exit(thread_id+1) + stokes(1)
       
    end do

    call photon_package
    
    do l = 1,3
       do k = 1,4
          do j = 1,npix
             do i = 1,npix

                if (l.eq.1) then
                   detector(i,j,k,l) = sum(detector_thread(i,j,k,l,:)) * package_energy
                else if (l.eq.2) then
                   detector(i,j,k,l) = sum(detector_thread(i,j,k,l,:)) * package_energy * package_energy
                else if (l.eq.3) then
                   detector(i,j,k,l) = sum(detector_thread(i,j,k,l,:))
                end if

             end do
          end do
       end do
    end do

    photometry = 0._dp
    
    photometry(1) = sum(detector(:,:,1,1)) ! Stokes I
    photometry(3) = sum(detector(:,:,2,1)) ! Stokes Q
    photometry(5) = sum(detector(:,:,3,1)) ! Stokes U
    photometry(7) = sum(detector(:,:,4,1)) ! Stokes V
    photometry(9) = sqrt( sum(detector(:,:,2,1))**2 + sum(detector(:,:,3,1))**2 ) ! PI
    if (photometry(1).gt.0._dp) photometry(10) = photometry(9) / photometry(1) ! PI/I

    do i=1,4

       if (sum(detector(:,:,i,3)).gt.0._dp) then
       
          dummy = ( sum(detector(:,:,i,2)) / sum(detector(:,:,i,3)) ) - ( sum(detector(:,:,i,1)) / sum(detector(:,:,i,3)) )**2
          if (dummy.gt.0._dp) photometry(i*2) = sqrt(dummy) * sqrt(sum(detector(:,:,i,3)))

       end if
       
    end do

    ! Degree of polarization error

    if (photometry(3)**2+photometry(5)**2.gt.0._dp) then
    
       dpi = sqrt( ( (photometry(3)*photometry(4))**2 + (photometry(5)*photometry(6))**2 ) / &
            ( 2._dp*(photometry(3)**2+photometry(5)**2 ) ) )

       photometry(11) = photometry(10) * sqrt( (dpi/photometry(9))**2 + (photometry(2)/photometry(1))**2  )

    end if

  end subroutine radiative_transfer

  subroutine emit_photon(thread_id, x_emission, y_emission, z_emission, emission_direction, &
       face_emission, cell_emission, r_disk, phi_disk, bias_weight)

    ! Photon emission

    integer,  intent(in)  :: thread_id
    integer,  intent(out) :: cell_emission(3), face_emission(2)
    integer               :: i, j, k
    real(dp), intent(out) :: x_emission, y_emission, z_emission, emission_direction(3), r_disk, phi_disk, bias_weight
    real(dp)              :: alpha, beta, rot_matrix(3,3), disk(2), emissivity_cumulative_sampled, cos_phi_sampled
    real(dp)              :: r_sampled, phi_sampled, phase_obs, sin_theta_sampled, cos_theta_sampled, emis_prev, xi
    real(dp)              :: theta_direction, phi_direction, x_temp, y_temp, z_temp, sin_phi_sampled
    real(dp)              :: radial_unit(3), y_bias, cos_phi_disk, sin_phi_disk
    logical               :: solution, grid_exit, cell_ok

    solution = .false.
    grid_exit = .false.
    bias_weight = 1._dp

    if (photon_source.eq.1) then

       ! Stellar photons

       face_emission(1) = 1
       face_emission(2) = nr

       r_disk = 0._dp
       phi_disk = 0._dp

       ! http://mathworld.wolfram.com/DiskPointPicking.html
          
       if (det_type.eq.3.and.det_phi*180._dp/pi.ge.170._dp) then

          do
             call random(thread_id,xi)
             r_disk = sqrt(xi)
             if (r_disk.gt.0.9_dp) exit
          end do
             
       else

          call random(thread_id,xi)
          r_disk = sqrt(xi)

       end if

       call random(thread_id,xi)
       phi_disk = 2._dp * pi * xi

       cos_phi_disk = cos(phi_disk)
       sin_phi_disk = sqrt(1._dp-cos_phi_disk*cos_phi_disk)

       if (phi_disk.gt.pi) sin_phi_disk = -sin_phi_disk

       disk(1) = rfront(nr) * r_disk * sin_phi_disk
       disk(2) = rfront(nr) * r_disk * cos_phi_disk

       emission_direction(1) = -1._dp
       emission_direction(2) = 0._dp
       emission_direction(3) = 0._dp

       x_emission = sqrt( rfront(nr)*rfront(nr) - disk(1)*disk(1) - disk(2)*disk(2) )
       y_emission = disk(1)
       z_emission = disk(2)

       if (abs(star_theta-pi/2._dp).gt.1.e-10_dp.or.abs(star_phi).gt.1.e-10_dp) then

          ! Rotate coordinates

          ! Rotate coordinates around y-axis
          call rotation_matrix(2, -(pi/2._dp-star_theta), rot_matrix)

          x_temp = x_emission*rot_matrix(1,1) + y_emission*rot_matrix(1,2) + z_emission*rot_matrix(1,3)
          y_temp = x_emission*rot_matrix(2,1) + y_emission*rot_matrix(2,2) + z_emission*rot_matrix(2,3)
          z_temp = x_emission*rot_matrix(3,1) + y_emission*rot_matrix(3,2) + z_emission*rot_matrix(3,3)

          ! Rotate coordinates around z-axis
          call rotation_matrix(3, star_phi, rot_matrix)

          x_emission = x_temp*rot_matrix(1,1) + y_temp*rot_matrix(1,2) + z_temp*rot_matrix(1,3)
          y_emission = x_temp*rot_matrix(2,1) + y_temp*rot_matrix(2,2) + z_temp*rot_matrix(2,3)
          z_emission = x_temp*rot_matrix(3,1) + y_temp*rot_matrix(3,2) + z_temp*rot_matrix(3,3)

          ! Rotate photon direction

          theta_direction = pi - star_theta
          phi_direction = pi + star_phi

          if (theta_direction.lt.0._dp) theta_direction = theta_direction + 2._dp*pi
          if (theta_direction.gt.2._dp*pi) theta_direction = theta_direction - 2._dp*pi
          if (phi_direction.lt.0._dp) phi_direction = phi_direction + 2._dp*pi
          if (phi_direction.gt.2._dp*pi) phi_direction = phi_direction - 2._dp*pi

          call spherical_cartesian(1._dp, theta_direction, phi_direction, emission_direction(1), &
               emission_direction(2), emission_direction(3))

       end if

       call initial_cell(x_emission, y_emission, z_emission, cell_emission)

    else if (photon_source.eq.2) then

       ! Thermal photon emission from the planet

       face_emission(1) = 0
       face_emission(2) = 0

       ! Sample grid cell

       call random(thread_id,xi)
       emissivity_cumulative_sampled = xi*emissivity_cumulative(nr-1,ntheta-1,nphi-1)

       emis_prev = 0._dp
       cell_ok = .false.
       
       do i=cell_depth,nr-1
          do j=0,ntheta-1
             do k=0,nphi-1

                if ( emissivity_cumulative_sampled.ge.emis_prev.and. &
                     emissivity_cumulative_sampled.le.emissivity_cumulative(i,j,k) ) then

                   cell_emission(1) = i
                   cell_emission(2) = j
                   cell_emission(3) = k

                   cell_ok = .true.
                   
                   exit

                end if

                emis_prev = emissivity_cumulative(i,j,k)

             end do
             if (cell_ok) exit
          end do
          if (cell_ok) exit
       end do

       ! Sample location in grid cell

       call random(thread_id,xi)
       r_sampled = xi * ( rfront(cell_emission(1)+1) - rfront(cell_emission(1)) )
       r_sampled = rfront(cell_emission(1)) + r_sampled

       call random(thread_id,xi)
       cos_theta_sampled = xi * ( theta_grid_cos(cell_emission(2)+1) - theta_grid_cos(cell_emission(2)) )
       cos_theta_sampled = theta_grid_cos(cell_emission(2)) + cos_theta_sampled
       sin_theta_sampled = sqrt(1._dp-cos_theta_sampled*cos_theta_sampled)

       if (nphi.eq.1) then

          call random(thread_id,xi)
          phi_sampled = 2._dp*pi*xi
          
       else if (nphi.gt.1) then

          if (cell_emission(3).lt.nphi-1) then

             call random(thread_id,xi)
             phi_sampled = xi * ( phifront(cell_emission(3)+1) - phifront(cell_emission(3)) )
             phi_sampled = phifront(cell_emission(3)) + phi_sampled

          else if (cell_emission(3).eq.nphi-1) then

             call random(thread_id,xi)
             phi_sampled = xi * ( 2._dp*pi - phifront(cell_emission(3)) )
             phi_sampled = phifront(cell_emission(3)) + phi_sampled

          end if

       end if

       cos_phi_sampled = cos(phi_sampled)
       sin_phi_sampled = sqrt(1._dp-cos_phi_sampled*cos_phi_sampled)
       if (phi_sampled.gt.pi) sin_phi_sampled = -sin_phi_sampled
       
       ! Emission coordinates
       
       x_emission = r_sampled*sin_theta_sampled*cos_phi_sampled
       y_emission = r_sampled*sin_theta_sampled*sin_phi_sampled
       z_emission = r_sampled*cos_theta_sampled

       ! Scale spherical coordinates to oblate spheroid coordinates

       x_emission = oblate_x * x_emission
       y_emission = oblate_y * y_emission
       z_emission = oblate_z * z_emission

       phase_obs = sin_theta_sampled*cos_phi_sampled*sin_det_theta*cos_det_phi + &
            sin_theta_sampled*sin_phi_sampled*sin_det_theta*sin_det_phi + cos_theta_sampled*cos_det_theta

       ! Sample emission direction

       if (photon_emission.eq.1) then

          ! Isotropic
          call random_direction(thread_id, emission_direction)
          
       else if (photon_emission.eq.2) then

          ! Biased upward (Gordon 1987)
          ! http://www.oceanopticsbook.info/view/monte_carlo_simulation/importance_sampling

          call random(thread_id,xi)
          y_bias = (1._dp+photon_bias) * tan(pi*xi/2._dp) / sqrt(1._dp-photon_bias*photon_bias)
          alpha = (1._dp-y_bias*y_bias) / (1._dp+y_bias*y_bias)

          call random(thread_id,xi)
          beta = 2._dp*pi*xi

          ! Local radial vector, normal vector to sphere
          radial_unit(1) = x_emission / ( oblate_x * oblate_x )
          radial_unit(2) = y_emission / ( oblate_y * oblate_y )
          radial_unit(3) = z_emission / ( oblate_z * oblate_z )
          radial_unit = radial_unit / sqrt(radial_unit(1)**2+radial_unit(2)**2+radial_unit(3)**2)

          ! Photon emission direction
          call direction_cosine(-alpha, beta, cos(2._dp*beta), radial_unit, emission_direction)

          ! Weight factor
          bias_weight = ( pi * sqrt(1._dp-alpha*alpha) * (1._dp+photon_bias*alpha) ) / &
               ( 2._dp * sqrt(1._dp-photon_bias*photon_bias) )

       end if

    end if

  end subroutine emit_photon
  
  subroutine rotation_matrix(rotation_axis, rotation_angle, matrix)

    ! Create rotation matrix for rotation aroung x-axis (1), y-axis (2) or z-axis (3)

    integer,  intent(in)  :: rotation_axis
    real(dp), intent(in)  :: rotation_angle
    real(dp), intent(out) :: matrix(3,3)
    real(dp)              :: cos_rot_angle, sin_rot_angle

    cos_rot_angle = cos(rotation_angle)
    sin_rot_angle = sin(rotation_angle)
    
    if (rotation_axis.eq.1) then

       matrix(1,1) = 1._dp
       matrix(1,2) = 0._dp
       matrix(1,3) = 0._dp

       matrix(2,1) = 0._dp
       matrix(2,2) = cos_rot_angle
       matrix(2,3) = -sin_rot_angle

       matrix(3,1) = 0._dp
       matrix(3,2) = sin_rot_angle
       matrix(3,3) = cos_rot_angle

    else if (rotation_axis.eq.2) then

       matrix(1,1) = cos_rot_angle
       matrix(1,2) = 0._dp
       matrix(1,3) = sin_rot_angle

       matrix(2,1) = 0._dp
       matrix(2,2) = 1._dp
       matrix(2,3) = 0._dp

       matrix(3,1) = -sin_rot_angle
       matrix(3,2) = 0._dp
       matrix(3,3) = cos_rot_angle

    else if (rotation_axis.eq.3) then

       matrix(1,1) = cos_rot_angle
       matrix(1,2) = -sin_rot_angle
       matrix(1,3) = 0._dp

       matrix(2,1) = sin_rot_angle
       matrix(2,2) = cos_rot_angle
       matrix(2,3) = 0._dp

       matrix(3,1) = 0._dp
       matrix(3,2) = 0._dp
       matrix(3,3) = 1._dp

    end if

  end subroutine rotation_matrix
  
  subroutine planck_function(temperature, flux)

    real(dp), intent(in)  :: temperature
    real(dp), intent(out) :: flux

    if (photon_source.eq.1) then

       ! Planck function [W m-2 m-1]
       flux = (2._dp*pi*hh*cc*cc/(wavelength**5._dp)) / ( exp(hh*cc/(wavelength*k_b*temperature)) - 1._dp )
       
    else if (photon_source.eq.2) then

       ! Planck function [W m-2 m-1 sr-1]
       flux = (2._dp*hh*cc*cc/(wavelength**5._dp)) / ( exp(hh*cc/(wavelength*k_b*temperature)) - 1._dp )

    end if

  end subroutine planck_function

  subroutine lambertian(thread_id, x, y, z, stokes, direction)

    ! Isotropic Lambertian reflection

    integer,  intent(in)    :: thread_id
    real(dp), intent(in)    :: x, y, z
    real(dp), intent(inout) :: stokes(4)
    real(dp), intent(out)   :: direction(3)
    real(dp)                :: xi, alpha, beta, surface_normal(3), norm

    ! Surface normal on spherical surface at point of reflection and scale for oblate surface
    surface_normal(1) = x / ( oblate_x * oblate_x )
    surface_normal(2) = y / ( oblate_y * oblate_y )
    surface_normal(3) = z / ( oblate_z * oblate_z )

    ! Make the vector a unit vector
    norm = sqrt(surface_normal(1)*surface_normal(1)+surface_normal(2)*surface_normal(2)+surface_normal(3)*surface_normal(3))
    surface_normal = surface_normal / norm

    ! Sample emission direction
    call random(thread_id,xi)
    alpha = sqrt(xi)
    call random(thread_id,xi)
    beta = 2._dp*pi*xi

    ! Calculate photon emission direction on a spherical planet
    call direction_cosine(alpha, beta, cos(2._dp*beta), surface_normal, direction)

    ! Lambertian surface depolarizes the light 100%
    stokes(2) = 0._dp
    stokes(3) = 0._dp
    stokes(4) = 0._dp

  end subroutine lambertian

  subroutine cartesian_spherical(x, y, z, r, theta, phi)

    ! Transform from Cartesian to spherical coordinates
    ! Note that atan2 returns a value in the range [-pi:pi]

    real(dp), intent(in)  :: x, y, z
    real(dp), intent(out) :: r, theta, phi

    r = sqrt(x*x+y*y+z*z)

    theta = acos(z/r)
    phi = atan2(y,x)

    if (phi.lt.0._dp) phi = phi+2._dp*pi

  end subroutine cartesian_spherical

  subroutine spherical_cartesian(r,theta,phi,x,y,z)

    ! Transform from spherical to Cartesian coordinates

    real(dp), intent(in)  :: r, theta, phi
    real(dp), intent(out) :: x, y, z

    x = r*sin(theta)*cos(phi)
    y = r*sin(theta)*sin(phi)
    z = r*cos(theta)

  end subroutine spherical_cartesian

  subroutine get_scatter(cell, scatter, acos_alpha)

    integer,  intent(in)  :: cell(3)
    integer               :: i, angle_upper, angle_lower
    real(dp), intent(in)  :: acos_alpha
    real(dp), intent(out) :: scatter(16)
    real(dp)              :: x0(16), x1(16), y0, y1, y_inter

    ! Scattering angle
    
    if (mod(acos_alpha*180._dp/pi,1._dp).gt.0.5_dp) then
       angle_upper = int(acos_alpha*180._dp/pi) + 2
       angle_lower = int(acos_alpha*180._dp/pi) + 1
    else
       angle_upper = int(acos_alpha*180._dp/pi) + 1
       angle_lower = int(acos_alpha*180._dp/pi)
    end if

    ! Scattering matrix
    
    if (angle_upper.eq.1) then

       scatter = cell_scatter_matrix(cell(1),cell(2),cell(3),wl_count,:,1)

    else if (angle_lower.eq.180) then

       scatter = cell_scatter_matrix(cell(1),cell(2),cell(3),wl_count,:,180)
       
    else

       do i=1,16

          x0(i) = cell_scatter_matrix(cell(1),cell(2),cell(3),wl_count,i,angle_lower)
          x1(i) = cell_scatter_matrix(cell(1),cell(2),cell(3),wl_count,i,angle_upper)
          y0 = dble(angle_lower) - 0.5_dp
          y1 = dble(angle_upper) - 0.5_dp
          y_inter = acos_alpha*180._dp/pi
          scatter(i) = (x1(i)-x0(i)) * (y_inter-y0) / (y1-y0) + x0(i)

       end do

    end if

  end subroutine get_scatter

  subroutine scattering_sampling(thread_id, stokes, alpha, beta, acos_alpha, c2b, cell)

    integer               :: i
    integer,  intent(in)  :: cell(3), thread_id
    real(dp), intent(in)  :: stokes(4)
    real(dp), intent(out) :: alpha, beta, acos_alpha, c2b
    real(dp)              :: intensity(180), intensity_cum(0:180), intensity_cum_sampled, xi, s2b

    ! Azimuthal angle sampling

    intensity_cum(0) = 0._dp

    do i=1,180

       ! Calculate the intensity for different azimuthal angles
       
       intensity(i) = cell_p11_int(cell(1),cell(2),cell(3),wl_count)*stokes(1) + &
            cell_p12_int(cell(1),cell(2),cell(3),wl_count)*stokes(2)*cos2beta(i) + &
            cell_p12_int(cell(1),cell(2),cell(3),wl_count)*stokes(3)*sin2beta(i) - &
            cell_p13_int(cell(1),cell(2),cell(3),wl_count)*stokes(2)*sin2beta(i) + &
            cell_p13_int(cell(1),cell(2),cell(3),wl_count)*stokes(3)*cos2beta(i) + &
            cell_p14_int(cell(1),cell(2),cell(3),wl_count)*stokes(4)

       intensity_cum(i) = intensity_cum(i-1) + intensity(i)

    end do

    call random(thread_id,xi)
    intensity_cum_sampled = xi*intensity_cum(180)

    do i=1,180

       if (intensity_cum_sampled.ge.intensity_cum(i-1).and.intensity_cum_sampled.le.intensity_cum(i)) then

          beta = (dble(i)-dble(i-1)) * (intensity_cum_sampled-intensity_cum(i-1)) / &
               (intensity_cum(i)-intensity_cum(i-1)) + dble(i-1)
          beta = beta*pi/180._dp
          exit

       end if

    end do

    call random(thread_id,xi)
    if (xi.gt.0.5_dp) beta = beta + pi
    
    if (beta.ge.2._dp*pi) beta = 2._dp*pi - 1.e-10_dp
    if (beta.le.0._dp) beta = 1.e-10_dp
    
    ! Scattering angle sampling

    c2b = cos(2._dp*beta)
    s2b = sqrt(1._dp-c2b*c2b)

    if (beta.gt.pi/2._dp.and.beta.lt.pi) then
       s2b = -s2b
    else if (beta.gt.3._dp*pi/2._dp.and.beta.lt.2._dp*pi) then
       s2b = -s2b
    else if (beta.lt.0._dp.or.beta.gt.2._dp*pi) then
       if (debug) then
          open (11, file=trim(error_log), position="append")
          write (11,*) "error 059"
          close (11)
       end if
    end if
    
    do i=1,180

       intensity(i) = cell_scatter_matrix(cell(1),cell(2),cell(3),wl_count,1,i)*stokes(1) + &
            cell_scatter_matrix(cell(1),cell(2),cell(3),wl_count,2,i)*c2b*stokes(2) + &
            cell_scatter_matrix(cell(1),cell(2),cell(3),wl_count,2,i)*s2b*stokes(3) - &
            cell_scatter_matrix(cell(1),cell(2),cell(3),wl_count,3,i)*s2b*stokes(2) + &
            cell_scatter_matrix(cell(1),cell(2),cell(3),wl_count,3,i)*c2b*stokes(3) + &
            cell_scatter_matrix(cell(1),cell(2),cell(3),wl_count,4,i)*stokes(4)

       intensity(i) = intensity(i) * sinbeta(i) * pi/180._dp

       intensity_cum(i) = intensity_cum(i-1) + intensity(i)

    end do

    call random(thread_id,xi)
    intensity_cum_sampled = xi*intensity_cum(180)

    do i=1,180

       if (intensity_cum_sampled.ge.intensity_cum(i-1).and.intensity_cum_sampled.le.intensity_cum(i)) then

          acos_alpha = (dble(i)-dble(i-1)) * (intensity_cum_sampled-intensity_cum(i-1)) / &
               (intensity_cum(i)-intensity_cum(i-1)) + dble(i-1)
          acos_alpha = acos_alpha*pi/180._dp
          alpha = cos(acos_alpha)
          exit

       end if

    end do

    if (alpha.ge.1._dp) alpha = 1._dp - 1.e-10_dp
    if (alpha.le.-1._dp) alpha = -1._dp + 1.e-10_dp

  end subroutine scattering_sampling

  subroutine stokes_rotation(alpha, beta, c2b, stokes_in, scatter, direction, direction_new, stokes_out, peeling)

    logical,  intent(in)  :: peeling
    real(dp), intent(in)  :: alpha, beta, c2b, stokes_in(4), scatter(16), direction_new(3), direction(3)
    real(dp), intent(out) :: stokes_out(4)
    real(dp)              :: stokes_rot(4), stokes_scattered(4), norm, cos_beta_2, s2b, s2b2, c2b2, mueller(4)
    real(dp)              :: num_check_1, beta_check

    if (debug) then
    
       ! Check if sampled beta and beta_check are the same
       
       if (abs(alpha).lt.1._dp.and.abs(direction(3)).lt.1._dp) then

          num_check_1 = (direction_new(3) - direction(3)*alpha) / &
               ( sqrt(1._dp-alpha*alpha)*sqrt(1._dp-direction(3)*direction(3)) )

          if (num_check_1.ge.1._dp)  then
             num_check_1 = 1._dp - 1.e-10_dp
          else if (num_check_1.le.-1._dp) then
             num_check_1 = -1._dp + 1.e-10_dp
          end if

          beta_check = acos(num_check_1)

          if (beta.ge.pi.and.beta.le.2._dp*pi) beta_check = 2._dp*pi - beta_check

          if (abs(beta-beta_check)*180._dp/pi.gt.1.e-2_dp) then

             if (debug) then
                open (11, file=trim(error_log), position="append")
                write (11,*) "error 009"
                write (11,*) beta, beta_check
                close (11)
             end if

          end if

       else

          if (debug) then
             open (11, file=trim(error_log), position="append")
             write (11,*) "error 010"
             close (11)
          end if

       end if

    end if

    ! Rotation from scattering plane to meridiane plane with the spherical cosine rule

    if (abs(direction_new(3)).lt.1._dp) then

       cos_beta_2 = (direction(3) - direction_new(3)*alpha) / &
            ( sqrt(1._dp-alpha*alpha) * sqrt(1._dp-direction_new(3)*direction_new(3)) )

       if (cos_beta_2.ge.1._dp) then
          cos_beta_2 = 1._dp - 1.e-10_dp
       else if (cos_beta_2.le.-1._dp) then
          cos_beta_2 = -1._dp + 1.e-10_dp
       end if

       s2b = sqrt(1._dp-c2b*c2b)

       if (beta.gt.pi/2._dp.and.beta.lt.pi) then
          s2b = -s2b
       else if (beta.gt.3._dp*pi/2._dp.and.beta.lt.2._dp*pi) then
          s2b = -s2b
       else if (beta.lt.0._dp.or.beta.gt.2._dp*pi) then
          if (debug) then
             open (11, file=trim(error_log), position="append")
             write (11,*) "error 060"
             close (11)
          end if
       end if

       mueller(1) = c2b
       mueller(2) = s2b
       mueller(3) = -s2b
       mueller(4) = c2b

       ! Rotate the Stokes vector from the meridian plane to the scattering plane

       stokes_rot(1) = stokes_in(1)
       stokes_rot(2) = mueller(1)*stokes_in(2) + mueller(2)*stokes_in(3)
       stokes_rot(3) = mueller(3)*stokes_in(2) + mueller(4)*stokes_in(3)
       stokes_rot(4) = stokes_in(4)

       ! Re-normalized the polarized intensity to keep it constant before and after the Stokes rotation

       if (sqrt(stokes_rot(2)**2+stokes_rot(3)**2+stokes_rot(4)**2).gt.0._dp) then
          norm = sqrt(stokes_in(2)**2+stokes_in(3)**2+stokes_in(4)**2) / sqrt(stokes_rot(2)**2+stokes_rot(3)**2+stokes_rot(4)**2)
       else
          norm = 1._dp
       end if

       if (norm.lt.1._dp.or.norm.gt.1._dp) then

          stokes_rot(2) = stokes_rot(2) * norm
          stokes_rot(3) = stokes_rot(3) * norm
          stokes_rot(4) = stokes_rot(4) * norm
          
       end if

       ! Apply scattering matrix

       stokes_scattered(1) = scatter(1)*stokes_rot(1) + scatter(2)*stokes_rot(2) + &
            scatter(3)*stokes_rot(3) + scatter(4)*stokes_rot(4)

       stokes_scattered(2) = scatter(5)*stokes_rot(1) + scatter(6)*stokes_rot(2) + &
            scatter(7)*stokes_rot(3) + scatter(8)*stokes_rot(4)

       stokes_scattered(3) = scatter(9)*stokes_rot(1) + scatter(10)*stokes_rot(2) + &
            scatter(11)*stokes_rot(3) + scatter(12)*stokes_rot(4)

       stokes_scattered(4) = scatter(13)*stokes_rot(1) + scatter(14)*stokes_rot(2) + &
            scatter(15)*stokes_rot(3) + scatter(16)*stokes_rot(4)

       ! Normalize the stokes vector such that intensity is conserved

       if (.not.peeling.and.stokes_scattered(1).gt.0._dp) then

          norm = stokes_rot(1) / stokes_scattered(1)
          stokes_scattered = norm * stokes_scattered

       end if

       ! Obtain Mueller matrix

       c2b2 = 2._dp*cos_beta_2*cos_beta_2-1._dp
       s2b2 = sqrt(1._dp-c2b2*c2b2)

       if (cos_beta_2.lt.0._dp) s2b2 = -s2b2

       if (beta.ge.0._dp.and.beta.lt.pi) then

          mueller(1) = c2b2
          mueller(2) = s2b2
          mueller(3) = -s2b2
          mueller(4) = c2b2

       else if (beta.ge.pi.and.beta.lt.2._dp*pi) then

          mueller(1) = c2b2
          mueller(2) = -s2b2
          mueller(3) = s2b2
          mueller(4) = c2b2
          
       end if
       
       ! Rotate the new Stokes vector into the meridian plane

       stokes_out(1) = stokes_scattered(1)
       stokes_out(2) = mueller(1)*stokes_scattered(2) + mueller(2)*stokes_scattered(3)
       stokes_out(3) = mueller(3)*stokes_scattered(2) + mueller(4)*stokes_scattered(3)
       stokes_out(4) = stokes_scattered(4)

       ! Re-normalized the polarized intensity to keep it constant before and after the Stokes rotation

       if (sqrt(stokes_out(2)**2+stokes_out(3)**2+stokes_out(4)**2).gt.0._dp) then
          norm = sqrt(stokes_scattered(2)**2+stokes_scattered(3)**2+stokes_scattered(4)**2) / &
               sqrt(stokes_out(2)**2+stokes_out(3)**2+stokes_out(4)**2)
       else
          norm = 1._dp
       end if

       stokes_out(2) = stokes_out(2) * norm
       stokes_out(3) = stokes_out(3) * norm
       stokes_out(4) = stokes_out(4) * norm

    else

       if (debug) then
          open (11, file=trim(error_log), position="append")
          write (11,*) "error 061"
          close (11)
       end if

    end if

  end subroutine stokes_rotation

  subroutine direction_cosine(alpha, beta, c2b, direction, direction_new)

    real(dp), intent(in)  :: alpha, beta, direction(3), c2b
    real(dp), intent(out) :: direction_new(3)
    real(dp)              :: cos_theta_old, sin_theta_old, cos_theta_new, sin_theta_new
    real(dp)              :: cos_phi_new, sin_phi_new, phi_old, phi_new, num_check_2, cos_beta

    ! Initital photon direction
    
    cos_theta_old = direction(3)/sqrt(direction(1)*direction(1)+direction(2)*direction(2)+direction(3)*direction(3))

    if (cos_theta_old.ge.1._dp) then
       cos_theta_old = 1._dp - 1.e-10_dp
    else if (cos_theta_old.le.-1._dp) then
       cos_theta_old = -1._dp + 1.e-10_dp
    end if
    
    sin_theta_old = sqrt(1._dp-cos_theta_old*cos_theta_old)

    phi_old = atan2(direction(2),direction(1))
    if (phi_old.lt.0._dp) phi_old = phi_old+2._dp*pi

    ! New theta direction from spherical cosine rule

    cos_beta = sqrt((c2b+1._dp)/2._dp)
    if (beta.gt.pi/2._dp.and.beta.lt.3._dp*pi/2._dp) cos_beta = -cos_beta

    cos_theta_new = cos_theta_old*alpha + sin_theta_old*sqrt(1._dp-alpha*alpha)*cos_beta
    sin_theta_new = sqrt(1._dp-cos_theta_new*cos_theta_new)

    ! New phi direction from spherical cosine rule

    num_check_2 = (alpha-cos_theta_new*cos_theta_old) / (sin_theta_new*sin_theta_old)

    if (num_check_2.ge.1._dp) then
       num_check_2 = 1._dp - 1.e-10_dp
    else if (num_check_2.le.-1._dp) then
       num_check_2 = -1._dp + 1.e-10_dp
    end if
    
    if (beta.ge.pi.and.beta.lt.2._dp*pi) then
       phi_new = phi_old - acos( num_check_2  )
    else if (beta.ge.0._dp.and.beta.lt.pi) then
       phi_new = phi_old + acos( num_check_2 )
    else
       if (debug) then
          open (11, file=trim(error_log), position="append")
          write (11,*) "error 019"
          close (11)
       end if
    end if

    if (phi_new.lt.0._dp) phi_new = phi_new + 2._dp*pi
    if (phi_new.gt.2._dp*pi) phi_new = phi_new - 2._dp*pi

    ! New photon direction
    
    cos_phi_new = cos(phi_new)
    
    if (phi_new.ge.0._dp.and.phi_new.lt.pi) then
       sin_phi_new = sqrt(1._dp-cos_phi_new*cos_phi_new)
    else if (phi_new.ge.pi.and.phi_new.le.2._dp*pi) then
       sin_phi_new = -sqrt(1._dp-cos_phi_new*cos_phi_new)
    else
       if (debug) then
          open (11, file=trim(error_log), position="append")
          write (11,*) "error 021"
          close (11)
       end if
    end if

    direction_new(1) = sin_theta_new*cos_phi_new
    direction_new(2) = sin_theta_new*sin_phi_new
    direction_new(3) = cos_theta_new

    if (direction_new(1).ge.1._dp) then
       direction_new(1) = 1._dp - 1.e-10_dp
       if (debug) then
          open (11, file=trim(error_log), position="append")
          write (11,*) "error 064"
          close (11)
       end if
    end if
    
    if (direction_new(2).ge.1._dp) then
       direction_new(2) = 1._dp - 1.e-10_dp
       if (debug) then
          open (11, file=trim(error_log), position="append")
          write (11,*) "error 065"
          close (11)
       end if
    end if
    
    if (direction_new(3).ge.1._dp) then
       direction_new(3) = 1._dp - 1.e-10_dp
       if (debug) then
          open (11, file=trim(error_log), position="append")
          write (11,*) "error 066"
          close (11)
       end if
    end if
    
    direction_new = direction_new / sqrt(direction_new(1)**2+direction_new(2)**2+direction_new(3)**2)
    
  end subroutine direction_cosine

  subroutine get_atmosphere

    integer               :: unit, status, readwrite, hdutype, nfound, blocksize, i, j, k, m, n
    integer, allocatable  :: naxes(:)
    real(dp), allocatable :: temp(:)
    logical               :: anynul

    status = 0
    readwrite = 0

    call ftgiou(unit, status)
    call ftopen(unit, "atmosphere.fits", readwrite, blocksize, status)
    call ftmahd(unit, 1, hdutype, status)

    ! Radial grid [m]

    call ftgknj(unit,'NAXIS',1,1,nr,nfound,status)
    allocate (temp(nr))
    temp = 0._dp
    call ftgpvd(unit,1,1,nr,-999._dp,temp,anynul,status)
    nr = nr-1
    allocate (rfront(0:nr))
    rfront = 0._dp
    do i = 0,nr
       rfront(i) = temp(i+1)
    end do
    deallocate (temp)

    ! Theta grid [rad]

    call ftmrhd(unit,1,hdutype,status)    
    call ftgknj(unit,'NAXIS',1,1,ntheta,nfound,status)
    allocate (temp(ntheta))
    temp = 0._dp
    call ftgpvd(unit,1,1,ntheta,-999._dp,temp,anynul,status)
    ntheta = ntheta-1
    allocate (thetafront(0:ntheta))
    thetafront = 0._dp
    allocate (thetaplane(0:ntheta))
    thetaplane = 0
    do i = 0,ntheta
       thetafront(i) = temp(i+1)
       if (thetafront(i).lt.90._dp-1.e-6_dp.or.thetafront(i).gt.90._dp+1.e-6_dp) then
          thetaplane(i) = 1
       else
          thetaplane(i) = 2
       end if
    end do
    deallocate (temp)
    thetafront = thetafront*pi/180.

    ! Phi grid [rad]

    call ftmrhd(unit,1,hdutype,status)
    call ftgknj(unit,'NAXIS',1,1,nphi,nfound,status)
    allocate (temp(nphi))
    temp = 0._dp
    call ftgpvd(unit,1,1,nphi,-999._dp,temp,anynul,status)
    nphi = nphi
    allocate (phifront(0:nphi-1))
    phifront = 0._dp
    do i = 0,nphi-1
       phifront(i) = temp(i+1)
    end do
    deallocate (temp)
    phifront = phifront*pi/180.

    ! Wavelengths [micron]

    call ftmrhd(unit,1,hdutype,status)
    call ftgknj(unit,'NAXIS',1,1,n_wavelength,nfound,status)
    allocate (wavelengths(n_wavelength))
    wavelengths = 0._dp
    call ftgpvd(unit,1,1,n_wavelength,-999._dp,wavelengths,anynul,status)
    wavelengths = wavelengths * 1.e-6_dp ! [micron] to [m]

    ! Density [kg m-3]

    allocate (cell_density(0:nr-1,0:ntheta-1,0:nphi-1))
    cell_density = 0._dp
    allocate (naxes(3))
    naxes = 0
    call ftmrhd(unit, 1, hdutype, status)
    call ftgknj(unit,'NAXIS',1,3,naxes,nfound,status)
    call ftgpvd(unit, 1, 1, naxes(1)*naxes(2)*naxes(3), -999._dp, cell_density, anynul, status)
    deallocate (cell_density)
    
    ! Temperature [K]

    allocate (cell_temperature(0:nr-1,0:ntheta-1,0:nphi-1))
    cell_temperature = 0._dp
    call ftmrhd(unit, 1, hdutype, status)
    call ftgknj(unit,'NAXIS',1,3,naxes,nfound,status)
    call ftgpvd(unit, 1, 1, naxes(1)*naxes(2)*naxes(3), -999._dp, cell_temperature, anynul, status)
    deallocate (naxes)
    
    ! Scattering opacity [m-1]

    allocate (naxes(4))
    naxes = 0
    allocate (cell_scattering_opacity(0:nr-1,0:ntheta-1,0:nphi-1,n_wavelength))
    cell_scattering_opacity = 0._dp
    call ftmrhd(unit, 1, hdutype, status)
    call ftgknj(unit,'NAXIS',1,4,naxes,nfound,status)
    call ftgpvd(unit, 1, 1, naxes(1)*naxes(2)*naxes(3)*naxes(4), -999._dp, cell_scattering_opacity, anynul, status)

    ! Absorption opacity [m-1]

    allocate (cell_absorption_opacity(0:nr-1,0:ntheta-1,0:nphi-1,n_wavelength))
    cell_absorption_opacity = 0._dp
    call ftmrhd(unit, 1, hdutype, status)
    call ftgknj(unit,'NAXIS',1,4,naxes,nfound,status)
    call ftgpvd(unit, 1, 1, naxes(1)*naxes(2)*naxes(3)*naxes(4), -999._dp, cell_absorption_opacity, anynul, status)
    deallocate (naxes)

    ! Opacity [m-1]

    allocate (cell_opacity(0:nr-1,0:ntheta-1,0:nphi-1,n_wavelength))
    allocate (cell_albedo(0:nr-1,0:ntheta-1,0:nphi-1,n_wavelength))
    cell_opacity = 0._dp
    cell_albedo = 0._dp
    do i=0,nr-1
       do j=0,ntheta-1
          do k=0,nphi-1
             do m=1,n_wavelength
                cell_opacity(i,j,k,m) = cell_scattering_opacity(i,j,k,m) + cell_absorption_opacity(i,j,k,m)
                if (cell_opacity(i,j,k,m).gt.0._dp) cell_albedo(i,j,k,m) = cell_scattering_opacity(i,j,k,m) / cell_opacity(i,j,k,m)
                if (cell_albedo(i,j,k,m).lt.1.e-20_dp) cell_albedo(i,j,k,m) = 1.e-20_dp
             end do
          end do
       end do
    end do
    
    ! Scattering matrix
    
    allocate (naxes(6))
    naxes = 0
    call ftmrhd(unit, 1, hdutype, status)
    call ftgknj(unit,'NAXIS',1,6,naxes,nfound,status)
    allocate (cell_scatter_matrix(0:nr-1,0:ntheta-1,0:nphi-1,n_wavelength,16,180))
    cell_scatter_matrix = 0._dp
    call ftgpvd(unit, 1, 1, naxes(1)*naxes(2)*naxes(3)*naxes(4)*naxes(5)*naxes(6), -999._dp, cell_scatter_matrix, anynul, status)
    deallocate (naxes)
    
    ! Asymmetry parameter

    allocate (naxes(4))
    naxes = 0
    allocate (cell_asymmetry(0:nr-1,0:ntheta-1,0:nphi-1,n_wavelength))
    cell_asymmetry = 0._dp
    call ftmrhd(unit, 1, hdutype, status)
    call ftgknj(unit,'NAXIS',1,4,naxes,nfound,status)
    call ftgpvd(unit, 1, 1, naxes(1)*naxes(2)*naxes(3)*naxes(4), -999._dp, cell_asymmetry, anynul, status)
    deallocate (naxes)

    call ftclos(unit, status)
    call ftfiou(unit, status)
    
    ! Integral of first row elements
    
    allocate (cell_p11_int(0:nr-1,0:ntheta-1,0:nphi-1,n_wavelength))
    allocate (cell_p12_int(0:nr-1,0:ntheta-1,0:nphi-1,n_wavelength))
    allocate (cell_p13_int(0:nr-1,0:ntheta-1,0:nphi-1,n_wavelength))
    allocate (cell_p14_int(0:nr-1,0:ntheta-1,0:nphi-1,n_wavelength))

    cell_p11_int = 0._dp
    cell_p12_int = 0._dp
    cell_p13_int = 0._dp
    cell_p14_int = 0._dp

    do i=0,nr-1
       do j=0,ntheta-1
          do k=0,nphi-1
             do m=1,n_wavelength
                do n=1,180

                   cell_p11_int(i,j,k,m) = cell_p11_int(i,j,k,m) + cell_scatter_matrix(i,j,k,m,1,n)*sinbeta(n)*pi/180._dp
                   cell_p12_int(i,j,k,m) = cell_p12_int(i,j,k,m) + cell_scatter_matrix(i,j,k,m,2,n)*sinbeta(n)*pi/180._dp
                   cell_p13_int(i,j,k,m) = cell_p13_int(i,j,k,m) + cell_scatter_matrix(i,j,k,m,3,n)*sinbeta(n)*pi/180._dp
                   cell_p14_int(i,j,k,m) = cell_p14_int(i,j,k,m) + cell_scatter_matrix(i,j,k,m,4,n)*sinbeta(n)*pi/180._dp

                end do
             end do
          end do
       end do
    end do

    ! Total number of cells
    cells = nr*ntheta*nphi

  end subroutine get_atmosphere

  subroutine grid_initialize(mode)

    ! Atmospheric grid initialization

    integer, intent(in) :: mode
    integer             :: i, j, k, cell_max, grid_out
    real(dp)            :: planck_flux, optical_depth_total1, optical_depth_total2, optical_depth_total
    real(dp)            :: weight_norm, emissivity_total
    logical             :: exist

    if (tau_max.lt.0._dp) then
       if (photon_source.eq.1) then
          tau_max = 30._dp
       else if (photon_source.eq.2) then
          tau_max = 5._dp
       end if
    end if

    if (mode.eq.1) then

       ! Wavelength independent arrays

       allocate (theta_grid_cos(0:ntheta))
       theta_grid_cos = 0._dp
       allocate (theta_grid_sin(0:ntheta))
       theta_grid_sin = 0._dp
       allocate (theta_grid_tan(0:ntheta))
       theta_grid_tan = 0._dp
       allocate (phi_grid_cos(0:nphi-1))
       phi_grid_cos = 0._dp
       allocate (phi_grid_sin(0:nphi-1))
       phi_grid_sin = 0._dp

       ! Set sine and cosine values of theta faces
       do i=0,ntheta
          theta_grid_cos(i) = cos(thetafront(i))
          theta_grid_sin(i) = sin(thetafront(i))
          theta_grid_tan(i) = tan(thetafront(i))
       end do

       ! Set sine and cosine values of phi faces
       do i=0,nphi-1
          phi_grid_cos(i) = cos(phifront(i))
          phi_grid_sin(i) = sin(phifront(i))
       end do

       ! Cell volume

       allocate (cell_volume(0:nr-1,0:ntheta-1,0:nphi-1))
       cell_volume = 0._dp

       do k=0,nphi-1
          do j=0,ntheta-1
             do i=0,nr-1

                ! Cell volume
                ! Volume = a * b * c * (1/3) * (r_out^3-r_in^3) * (cos(theta_in)-cos(theta_out)) * (phi_out-phi_in)

                if (nphi.eq.1) then
                
                   cell_volume(i,j,k) = oblate_x * oblate_y * oblate_z * (1._dp/3._dp) * &
                        (rfront(i+1)**3-rfront(i)**3) * (theta_grid_cos(j)-theta_grid_cos(j+1)) * 2._dp*pi

                else if (nphi.gt.1) then

                   if (k.lt.nphi-1) then

                      cell_volume(i,j,k) = oblate_x * oblate_y * oblate_z * (1._dp/3._dp) * &
                           (rfront(i+1)**3-rfront(i)**3) * (theta_grid_cos(j)-theta_grid_cos(j+1)) * (phifront(k+1)-phifront(k))

                   else if (k.eq.nphi-1) then

                      cell_volume(i,j,k) = oblate_x * oblate_y * oblate_z * (1._dp/3._dp) * &
                           (rfront(i+1)**3-rfront(i)**3) * (theta_grid_cos(j)-theta_grid_cos(j+1)) * (2._dp*pi-phifront(k))

                   end if
                      
                end if

             end do
          end do
       end do

       ! Cell energy transport
       
       if (global) then
       
          allocate (cell_flow_global(threads,3,0:nr-1,0:ntheta-1,0:nphi-1))
          cell_flow_global = 0._dp

       end if

       if (latitudinal) then

          allocate (cell_flow(threads,4,0:nr-1,0:ntheta-1,0:nphi-1))
          cell_flow = 0._dp

       end if

    else if (mode.eq.2) then

       ! Wavelength dependent arrays

       ! Check at which radial depth the optical thickness is tau_max or larger
       ! Determine at which theta/phi location this is deepest in the atmosphere
       
       if (photon_source.eq.1) then

          if (photon_walk.lt.0._dp) then
          
             cell_max = 1000000

             do j=0,ntheta-1
                do k=0,nphi-1

                   optical_depth_total = 0._dp

                   do i=0,nr-1

                      optical_depth_total = optical_depth_total + cell_opacity(nr-i-1,j,k,wl_count) * &
                           ( rfront(nr-i)-rfront(nr-i-1) )

                      cell_depth = nr-i-1

                      if (optical_depth_total.gt.tau_max) exit

                   end do

                   if (cell_depth.lt.cell_max) cell_max = cell_depth

                end do
             end do

             cell_depth = cell_max

          end if

       else if (photon_source.eq.2) then

          cell_max = 1000000

          if (ring) then
             grid_out = 2
          else if (.not.ring) then
             grid_out = 0
          end if
          
          do j=0,ntheta-1
             do k=0,nphi-1

                optical_depth_total = 0._dp

                do i=grid_out,nr-1

                   optical_depth_total = optical_depth_total + cell_absorption_opacity(nr-i-1,j,k,wl_count) * &
                        ( rfront(nr-i)-rfront(nr-i-1) )

                   cell_depth = nr-i-1

                   if (optical_depth_total.gt.tau_max) exit

                end do

                if (cell_depth.lt.cell_max) cell_max = cell_depth

             end do
          end do

          cell_depth = cell_max

          ! Thermal cell luminosity
          
          weight_norm = 0._dp

          do i=cell_depth,nr-1
             do j=0,ntheta-1
                do k=0,nphi-1

                   if (cell_temperature(i,j,k).gt.0._dp) then
                      
                      call planck_function(cell_temperature(i,j,k), planck_flux)
                      weight_norm = weight_norm + cell_absorption_opacity(i,j,k,wl_count)*planck_flux*cell_volume(i,j,k)

                   end if

                end do
             end do
          end do

          allocate (cell_weight(0:nr-1,0:ntheta-1,0:nphi-1))
          cell_weight = 0._dp

          allocate (cell_luminosity(0:nr-1,0:ntheta-1,0:nphi-1))
          cell_luminosity = 0._dp

          allocate (emissivity_cumulative(0:nr-1,0:ntheta-1,0:nphi-1))
          emissivity_cumulative = 0._dp

          emissivity_total = 0._dp
          
          do i=cell_depth,nr-1
             do j=0,ntheta-1
                do k=0,nphi-1

                   if (cell_temperature(i,j,k).gt.0._dp.and.cell_absorption_opacity(i,j,k,wl_count).gt.0._dp) then
                      
                      call planck_function(cell_temperature(i,j,k), planck_flux)

                      if (thermal_weight) then
                         cell_weight(i,j,k) = weight_norm / (cell_volume(i,j,k)*cell_absorption_opacity(i,j,k,wl_count)*planck_flux)
                      else if (.not.thermal_weight) then
                         cell_weight(i,j,k) = 1._dp
                      end if
                      
                      cell_luminosity(i,j,k) = 4._dp*pi*cell_volume(i,j,k)*cell_absorption_opacity(i,j,k,wl_count)*planck_flux ! [W m-1]

                      emissivity_cumulative(i,j,k) = emissivity_total + cell_luminosity(i,j,k) * cell_weight(i,j,k)

                      emissivity_total = emissivity_cumulative(i,j,k)

                   else

                      emissivity_cumulative(i,j,k) = emissivity_total

                   end if

                end do
             end do
          end do

       end if

       ! Optical depth output

       if (det_type.eq.2) then

          inquire(file=trim(output_name)//'/output/optical_depth.dat', exist=exist)

          if (exist) then

             open (100, file=trim(output_name)//'/output/optical_depth.dat', status='old', position='append', action='write')

          else if (.not.exist) then

             open (100, file=trim(output_name)//'/output/optical_depth.dat', status='new', action='write')

             write (100,*) "# Wavelength [micron] - Total optical depth - Absorption optical depth - Scattering optical depth"
             write (100,*)

          end if

          optical_depth_total = 0._dp
          optical_depth_total1 = 0._dp
          optical_depth_total2 = 0._dp

          do i=0,nr-1
             optical_depth_total = optical_depth_total + ( (rfront(i+1)-rfront(i)) * cell_opacity(i,0,0,wl_count) )
             optical_depth_total1 = optical_depth_total1 + ( (rfront(i+1)-rfront(i)) * cell_scattering_opacity(i,0,0,wl_count) )
             optical_depth_total2 = optical_depth_total2 + ( (rfront(i+1)-rfront(i)) * cell_absorption_opacity(i,0,0,wl_count) )
          end do

          ! Wavelength [micron] - Total radial optical depth total - absorption - scattering
          write (100,*) wavelength*1.e6_dp, optical_depth_total, optical_depth_total2, optical_depth_total1

          close (100)

       end if

       if (photon_source.eq.2) then
       
          allocate (flux_emitted(threads))
          flux_emitted = 0._dp

       end if
          
       allocate (flux_exit(threads))
       flux_exit = 0._dp
          
    end if

  end subroutine grid_initialize

  subroutine photon_package

    ! Photon package energy

    real(dp) :: planck_flux

    if (photon_source.eq.1) then
    
       ! L_star = 4 * pi * R_star^2 * F_lambda = 4 * pi * D^2 * F_lambda,planet
       ! F_lambda,planet = F_lambda * R_star^2/D^2
       ! F_lambda,planet * pi * R_planet^2 = d^2 * F_obs
       ! F_obs = pi * R_planet^2 * R_star^2 * F_lambda / ( D^2 * d^2 )
    
       call planck_function(t_star, planck_flux)

       package_energy = pi * planck_flux * rfront(nr) * rfront(nr) * r_star * r_star / &
            ( orbit * orbit * distance_planet * distance_planet * dble(packages) )

       if (det_type.eq.3.and.det_phi*180._dp/pi.ge.170._dp) then
          
          package_energy = package_energy * (pi*r_star*r_star-0.9_dp*0.9_dp*pi*r_star*r_star)/(pi*r_star*r_star)

       end if

    else if (photon_source.eq.2) then
       
       package_energy = emissivity_cumulative(nr-1,ntheta-1,nphi-1) / ( distance_planet * distance_planet * dble(packages) )

    end if

  end subroutine photon_package

  subroutine array_start

    allocate (detector(npix,npix,4,3))
    detector = 0._dp
    allocate (detector_thread(npix,npix,4,3,threads))
    detector_thread = 0._dp

  end subroutine array_start

  subroutine grid_finished(i)

    integer, intent(in) :: i

    if (i.eq.1) then

       ! Wavelength dependent arrays

       if (photon_source.eq.2) then
       
          deallocate (cell_luminosity)
          deallocate (cell_weight)
          deallocate (emissivity_cumulative)
          deallocate (flux_emitted)

       end if

       deallocate (flux_exit)
       
    else if (i.eq.2) then

       ! Detector arrays

       deallocate (detector)
       deallocate (detector_thread)

    else if (i.eq.3) then

       ! Global grid arrays

       deallocate (cell_absorption_opacity)
       deallocate (cell_scattering_opacity)
       deallocate (cell_asymmetry)
       deallocate (cell_opacity)
       deallocate (cell_albedo)
       deallocate (cell_scatter_matrix)
       deallocate (cell_p11_int)
       deallocate (cell_p12_int)
       deallocate (cell_p13_int)
       deallocate (cell_p14_int)
       deallocate (rfront)
       deallocate (thetafront)
       deallocate (thetaplane)
       deallocate (phifront)
       deallocate (theta_grid_cos)
       deallocate (theta_grid_sin)
       deallocate (theta_grid_tan)
       deallocate (phi_grid_cos)
       deallocate (phi_grid_sin)
       deallocate (cell_volume)
       deallocate (wavelengths)
       if (global) deallocate (cell_flow_global)
       if (latitudinal) deallocate (cell_flow)
       deallocate (phi_mrw)
       deallocate (y_mrw)
       
    end if

  end subroutine grid_finished

  subroutine initial_cell(x, y, z, cell)

    ! Initial grid cell that is crossed with reflected light

    integer               :: j
    real(dp), intent(in)  :: x, y, z
    real(dp)              :: theta, phi
    integer,  intent(out) :: cell(3)

    ! Radial grid cell
    
    cell(1) = nr-1
    cell(2) = 0
    cell(3) = 0

    ! Polar grid cell
    
    if (ntheta.gt.1.) then
       
       theta = acos(z/sqrt(x*x+y*y+z*z))
       
       do j=0,ntheta-1
          
          if (theta.gt.thetafront(j).and.theta.lt.thetafront(j+1)) then
             
             cell(2) = j
             exit
             
          end if
          
       end do
       
    end if

    ! Azimuthal grid cell

    if (nphi.gt.1) then
       
       phi = atan2(y,x)

       if (phi.lt.0._dp) phi = phi+2._dp*pi
       
       do j=0,nphi-1
          
          if (j.lt.nphi-1) then
             
             if (phi.gt.phifront(j).and.phi.lt.phifront(j+1)) then
                
                cell(3) = j
                exit
                
             end if
             
          else if (j.eq.nphi-1) then
             
             if (phi.gt.phifront(j).and.phi.lt.2._dp*pi) then
                
                cell(3) = j
                exit
                
             end if
             
          end if
          
       end do
       
    end if
       
  end subroutine initial_cell

  subroutine next_cell(face, next_face, cell_in, cell_out)

    integer, intent(in)  :: face(2), next_face(2), cell_in(3)
    integer, intent(out) :: cell_out(3)

    cell_out = cell_in

    ! Next face is radial
    if (next_face(1).eq.1) then

       if (face(1).eq.1.and.next_face(2).eq.face(2)) then

          ! Next cell is radially outward, photon crosses the same radial surface
          cell_out(1) = cell_in(1)+1

       else if (next_face(2).eq.cell_in(1)) then

          ! Next cell is radially inward
          cell_out(1) = cell_in(1)-1

       else if (next_face(2).eq.cell_in(1)+1) then

          ! Next cell is radially outward
          cell_out(1) = cell_in(1)+1

       else

          if (debug) then
             open (11, file=trim(error_log), position="append")
             write (11,*) "error 022"
             close (11)
          end if

       end if

    end if

    ! Next face is polar
    if (next_face(1).eq.2) then

       if (face(1).eq.2.and.next_face(2).eq.face(2).and.thetafront(next_face(2)).lt.pi/2._dp) then

          ! Next cell is theta outward, photon crosses the same theta surface
          cell_out(2) = cell_in(2)+1

       else if (face(1).eq.2.and.next_face(2).eq.face(2).and.thetafront(next_face(2)).gt.pi/2._dp) then

          ! Next cell is theta inward, photon crosses the same theta surface
          cell_out(2) = cell_in(2)-1

       else if (next_face(2).eq.cell_in(2)) then

          ! Next cell is theta inward
          cell_out(2) = cell_in(2)-1

       else if (next_face(2).eq.cell_in(2)+1) then

          ! Next cell is theta outward
          cell_out(2) = cell_in(2)+1

       else

          if (debug) then
             open (11, file=trim(error_log), position="append")
             write (11,*) "error 023"
             close (11)
          end if

       end if

    end if

    ! Change in phi grid cell
    if (next_face(1).eq.3) then

       if (cell_in(3).eq.nphi-1.and.next_face(2).eq.0) then

          ! Next cell is azimuthally outward
          cell_out(3) = 0

       else if (cell_in(3).eq.0.and.next_face(2).eq.0) then

          ! Next cell is azimuthally inward
          cell_out(3) = nphi-1

       else if (next_face(2).eq.cell_in(3)+1) then

          ! Next cell is azimuthally outward
          cell_out(3) = cell_in(3)+1

       else if (next_face(2).eq.cell_in(3)) then

          ! Next cell is azimuthally inward
          cell_out(3) = cell_in(3)-1

       else

          if (debug) then
             open (11, file=trim(error_log), position="append")
             write (11,*) "error 024"
             close (11)
          end if

       end if

    end if

  end subroutine next_cell

  subroutine cell_face(photon, x, y, z, n, face, next_face, face_distance, grid_exit, cell, cell_out, cell_error)

    integer,  intent(in)    :: face(2), cell(3)
    integer(16), intent(in) :: photon
    integer,  intent(out)   :: next_face(2), cell_out(3)
    integer                 :: r_in, r_out, r_same, t_in, t_out, t_same, p_in, p_out
    real(dp), intent(in)    :: x, y, z, n(3)
    real(dp), intent(out)   :: face_distance
    real(dp)                :: sol_r_in(2), sol_r_out(2), sol_r_same(2), sol_t_in(2), sol_t_out(2), sol_t_same(2)
    real(dp)                :: sol_p_in, sol_p_out, a, b, c, qa, qb, qc, z_test, tan_theta, aa, bb, cc
    logical,  intent(out)   :: grid_exit, cell_error

    grid_exit = .false.
    cell_error = .false.
    next_face = 0

    face_distance = 1.e20_dp
    
    sol_r_in = 0._dp
    sol_r_out = 0._dp
    sol_r_same = 0._dp

    sol_t_in = 0._dp
    sol_t_out = 0._dp
    sol_t_same = 0._dp

    sol_p_in = 0._dp
    sol_p_out = 0._dp

    a = 1._dp/oblate_x
    b = 1._dp/oblate_y
    c = 1._dp/oblate_z

    r_in = cell(1)
    r_out = cell(1) + 1

    t_in = cell(2)
    t_out = cell(2) + 1

    p_in = cell(3)
    p_out = cell(3) + 1
    if (p_out.eq.nphi) p_out = 0

    if (face(1).eq.1) then

       r_in = face(2)-1
       r_out = face(2)+1
       r_same = face(2)

    else if (face(1).eq.2) then

       t_in = face(2)-1
       t_out = face(2)+1
       t_same = face(2)

    else if (face(1).eq.3) then

       if (face(2).eq.0) then
          p_in = nphi-1
       else
          p_in = face(2)-1
       end if
       
       if (face(2).eq.nphi-1) then
          p_out = 0
       else
          p_out = face(2)+1
       end if

    end if

    ! Radial
    
    qa = a*a*n(1)*n(1) + b*b*n(2)*n(2) + c*c*n(3)*n(3)
    qb = 2._dp*(a*a*x*n(1) + b*b*y*n(2) + c*c*z*n(3))
    qc = a*a*x*x + b*b*y*y + c*c*z*z
    
    if (face(1).eq.1) then

       ! Current cell face is a radial face
       
       if (cell(1).eq.face(2)-1) then

          ! Radial inward face, photon comes from radial outward cell
          call quadratic_equation(qa, qb, qc-rfront(r_in)*rfront(r_in), sol_r_in)

          ! Radial same face, only possible for radially inward moving photon
          call quadratic_equation(qa, qb, qc-rfront(r_same)*rfront(r_same), sol_r_same)

       else if (cell(1).eq.face(2)) then

          ! Radial outward face, photon comes from radial inward cell
          call quadratic_equation(qa, qb, qc-rfront(r_out)*rfront(r_out), sol_r_out)

       end if

    else if (face(1).ne.1) then

       ! Cell face is not a radial cell face
       
       ! Radial inward face
       call quadratic_equation(qa, qb, qc-rfront(r_in)*rfront(r_in), sol_r_in)

       ! Radial outward face
       call quadratic_equation(qa, qb, qc-rfront(r_out)*rfront(r_out), sol_r_out)

    end if

    if (sol_r_in(1).gt.grid_min.and.sol_r_in(1).lt.1.e20_dp.and.sol_r_in(1).lt.face_distance) then
       face_distance = sol_r_in(1)
       next_face(1) = 1
       next_face(2) = r_in
    end if

    if (sol_r_in(2).gt.grid_min.and.sol_r_in(2).lt.1.e20_dp.and.sol_r_in(2).lt.face_distance) then
       face_distance = sol_r_in(2)
       next_face(1) = 1
       next_face(2) = r_in
    end if

    if (sol_r_out(1).gt.grid_min.and.sol_r_out(1).lt.1.e20_dp.and.sol_r_out(1).lt.face_distance) then
       face_distance = sol_r_out(1)
       next_face(1) = 1
       next_face(2) = r_out
    end if

    if (sol_r_out(2).gt.grid_min.and.sol_r_out(2).lt.1.e20_dp.and.sol_r_out(2).lt.face_distance) then
       face_distance = sol_r_out(2)
       next_face(1) = 1
       next_face(2) = r_out
    end if

    if (sol_r_same(1).gt.1.e-3_dp.and.sol_r_same(1).lt.1.e20_dp.and.sol_r_same(1).lt.face_distance) then
       face_distance = sol_r_same(1)
       next_face(1) = 1
       next_face(2) = r_same
    end if

    if (sol_r_same(2).gt.1.e-3_dp.and.sol_r_same(2).lt.1.e20_dp.and.sol_r_same(2).lt.face_distance) then
       face_distance = sol_r_same(2)
       next_face(1) = 1
       next_face(2) = r_same
    end if

    ! Polar

    if (ntheta.gt.1) then

       qa = a*a*n(1)*n(1) + b*b*n(2)*n(2)
       qb = 2._dp*(a*a*x*n(1) + b*b*y*n(2))
       qc = a*a*x*x + b*b*y*y

       aa = c*c*n(3)*n(3)
       bb = 2._dp*c*c*z*n(3)
       cc = c*c*z*z
       
       if (face(1).eq.2) then

          ! Current cell face is a theta face
          
          if (cell(2).eq.face(2)-1.and.t_in.ne.0) then

             ! Theta inward face, photon coming from theta outward face

             if (thetaplane(t_in).eq.1) then

                tan_theta = theta_grid_tan(t_in)*theta_grid_tan(t_in)
                call quadratic_equation(qa-aa*tan_theta, qb-bb*tan_theta, qc-cc*tan_theta, sol_t_in)
                
                if (sol_t_in(1).gt.grid_min) then

                   z_test = z+sol_t_in(1)*n(3)

                   ! if (abs(sign(1._dp,z_test)-sign(1._dp,pi/2._dp-thetafront(t_in))).gt.1.e-10_dp)

                   if ( (z_test.gt.0._dp.and.thetafront(t_in).gt.pi/2._dp) .or. &
                        (z_test.lt.0._dp.and.thetafront(t_in).lt.pi/2._dp) ) sol_t_in(1) = 0._dp
                      
                end if

                if (sol_t_in(2).gt.grid_min) then

                   z_test = z+sol_t_in(2)*n(3)

                   if ( (z_test.gt.0._dp.and.thetafront(t_in).gt.pi/2._dp) .or. &
                        (z_test.lt.0._dp.and.thetafront(t_in).lt.pi/2._dp) ) sol_t_in(2) = 0._dp

                end if

             else if (thetaplane(t_in).eq.2) then

                if (abs(n(3)).gt.0._dp) sol_t_in(1) = -z/n(3)

             end if

          else if (cell(2).eq.face(2).and.t_out.ne.ntheta) then

             ! Theta outward face, photon coming from inward face

             if (thetaplane(t_out).eq.1) then

                tan_theta = theta_grid_tan(t_out)*theta_grid_tan(t_out)
                call quadratic_equation(qa-aa*tan_theta, qb-bb*tan_theta, qc-cc*tan_theta, sol_t_out)
                
                if (sol_t_out(1).gt.grid_min) then

                   z_test = z+sol_t_out(1)*n(3)

                   if ( (z_test.gt.0._dp.and.thetafront(t_out).gt.pi/2._dp) .or. &
                        (z_test.lt.0._dp.and.thetafront(t_out).lt.pi/2._dp) ) sol_t_out(1) = 0._dp

                end if

                if (sol_t_out(2).gt.grid_min) then

                   z_test = z+sol_t_out(2)*n(3)

                   if ( (z_test.gt.0._dp.and.thetafront(t_out).gt.pi/2._dp) .or. &
                        (z_test.lt.0._dp.and.thetafront(t_out).lt.pi/2._dp) ) sol_t_out(2) = 0._dp

                end if

             else if (thetaplane(t_out).eq.2) then

                if (abs(n(3)).gt.0._dp) sol_t_out(1) = -z/n(3)

             end if

          end if

          if ( (thetafront(t_same).lt.pi/2._dp.and.cell(2).eq.face(2)-1) .or. &
               (thetafront(t_same).gt.pi/2._dp.and.cell(2).eq.face(2)) ) then

             ! Theta same face

             if (thetaplane(t_same).eq.1) then

                tan_theta = theta_grid_tan(t_same)*theta_grid_tan(t_same)
                call quadratic_equation(qa-aa*tan_theta, qb-bb*tan_theta, qc-cc*tan_theta, sol_t_same)

                if (sol_t_same(1).gt.grid_min) then

                   z_test = z+sol_t_same(1)*n(3)

                   ! if (sign(1._dp,z_test)-sign(1._dp,pi/2._dp-thetafront(t_same)).gt.1.e-10_dp) then

                   if ( (z_test.gt.0._dp.and.thetafront(t_same).gt.pi/2._dp) .or. &
                        (z_test.lt.0._dp.and.thetafront(t_same).lt.pi/2._dp) ) sol_t_same(1) = 0._dp

                end if

                if (sol_t_same(2).gt.grid_min) then

                   z_test = z+sol_t_same(2)*n(3)

                   if ( (z_test.gt.0._dp.and.thetafront(t_same).gt.pi/2._dp) .or. &
                        (z_test.lt.0._dp.and.thetafront(t_same).lt.pi/2._dp) ) sol_t_same(2) = 0._dp

                end if

             end if

          end if

       else if (face(1).ne.2) then

          ! Current face is not a theta face

          ! Check if the inward theta face exists

          if (t_in.ne.0) then

             ! Theta inward face

             if (thetaplane(t_in).eq.1) then

                tan_theta = theta_grid_tan(t_in)*theta_grid_tan(t_in)
                call quadratic_equation(qa-aa*tan_theta, qb-bb*tan_theta, qc-cc*tan_theta, sol_t_in)

                if (sol_t_in(1).gt.grid_min) then

                   z_test = z+sol_t_in(1)*n(3)

                   if ( (z_test.gt.0._dp.and.thetafront(t_in).gt.pi/2._dp) .or. &
                        (z_test.lt.0._dp.and.thetafront(t_in).lt.pi/2._dp) ) sol_t_in(1) = 0._dp

                end if

                if (sol_t_in(2).gt.grid_min) then

                   z_test = z+sol_t_in(2)*n(3)

                   if ( (z_test.gt.0._dp.and.thetafront(t_in).gt.pi/2._dp) .or. &
                        (z_test.lt.0._dp.and.thetafront(t_in).lt.pi/2._dp) ) sol_t_in(2) = 0._dp

                end if

             else if (thetaplane(t_in).eq.2) then

                if (abs(n(3)).gt.0._dp) sol_t_in(1) = -z/n(3)

             end if

          end if

          if (t_out.ne.ntheta) then

             ! Theta outward face

             if (thetaplane(t_out).eq.1) then

                tan_theta = theta_grid_tan(t_out)*theta_grid_tan(t_out)
                call quadratic_equation(qa-aa*tan_theta, qb-bb*tan_theta, qc-cc*tan_theta, sol_t_out)

                if (sol_t_out(1).gt.grid_min) then

                   z_test = z+sol_t_out(1)*n(3)

                   if ( (z_test.gt.0._dp.and.thetafront(t_out).gt.pi/2._dp) .or. &
                        (z_test.lt.0._dp.and.thetafront(t_out).lt.pi/2._dp) ) sol_t_out(1) = 0._dp

                end if

                if (sol_t_out(2).gt.grid_min) then

                   z_test = z+sol_t_out(2)*n(3)

                   if ( (z_test.gt.0._dp.and.thetafront(t_out).gt.pi/2._dp) .or. &
                        (z_test.lt.0._dp.and.thetafront(t_out).lt.pi/2._dp) ) sol_t_out(2) = 0._dp

                end if

             else if (thetaplane(t_out).eq.2) then

                if (abs(n(3)).gt.0._dp) sol_t_out(1) = -z/n(3)

             end if

          end if

       end if

       if (sol_t_in(1).gt.grid_min.and.sol_t_in(1).lt.1.e20_dp.and.sol_t_in(1).lt.face_distance) then
          face_distance = sol_t_in(1)
          next_face(1) = 2
          next_face(2) = t_in
       end if

       if (sol_t_in(2).gt.grid_min.and.sol_t_in(2).lt.1.e20_dp.and.sol_t_in(2).lt.face_distance) then
          face_distance = sol_t_in(2)
          next_face(1) = 2
          next_face(2) = t_in
       end if

       if (sol_t_out(1).gt.grid_min.and.sol_t_out(1).lt.1.e20_dp.and.sol_t_out(1).lt.face_distance) then
          face_distance = sol_t_out(1)
          next_face(1) = 2
          next_face(2) = t_out
       end if

       if (sol_t_out(2).gt.grid_min.and.sol_t_out(2).lt.1.e20_dp.and.sol_t_out(2).lt.face_distance) then
          face_distance = sol_t_out(2)
          next_face(1) = 2
          next_face(2) = t_out
       end if

       if (sol_t_same(1).gt.1.e-3_dp.and.sol_t_same(1).lt.1.e20_dp.and.sol_t_same(1).lt.face_distance) then
          face_distance = sol_t_same(1)
          next_face(1) = 2
          next_face(2) = t_same
       end if

       if (sol_t_same(2).gt.1.e-3_dp.and.sol_t_same(2).lt.1.e20_dp.and.sol_t_same(2).lt.face_distance) then
          face_distance = sol_t_same(2)
          next_face(1) = 2
          next_face(2) = t_same
       end if
       
    end if

    ! Azimuthal
    
    if (nphi.gt.1) then
       
       if (face(1).eq.3) then

          ! Current face is a phi face

          if (cell(3).eq.face(2)-1.or.(cell(3).eq.nphi-1.and.face(2).eq.0)) then

             ! Photon is travelling inward in phi direction

             if (abs(b*n(2)*phi_grid_cos(p_in)-a*n(1)*phi_grid_sin(p_in)).gt.0._dp) then

                sol_p_in = (a*x*phi_grid_sin(p_in)-b*y*phi_grid_cos(p_in)) / &
                     (b*n(2)*phi_grid_cos(p_in)-a*n(1)*phi_grid_sin(p_in))

             end if

          else if (cell(3).eq.face(2)) then

             ! Photon is travelling outward in phi direction

             if (abs((b*n(2)*phi_grid_cos(p_out)-a*n(1)*phi_grid_sin(p_out))).gt.0._dp) then

                sol_p_out = (a*x*phi_grid_sin(p_out)-b*y*phi_grid_cos(p_out)) / &
                     (b*n(2)*phi_grid_cos(p_out)-a*n(1)*phi_grid_sin(p_out))

             end if

          end if

       else if (face(1).ne.3.and.nphi.gt.1) then

          ! Current face in not a phi face
          ! Inward phi direction

          if (abs(b*n(2)*phi_grid_cos(p_in)-a*n(1)*phi_grid_sin(p_in)).gt.0._dp) then

             sol_p_in = (a*x*phi_grid_sin(p_in)-b*y*phi_grid_cos(p_in)) / &
                  (b*n(2)*phi_grid_cos(p_in)-a*n(1)*phi_grid_sin(p_in))

          end if

          ! Outward phi direction

          if (abs(n(2)*phi_grid_cos(p_out)-n(1)*phi_grid_sin(p_out)).gt.0._dp) then

             sol_p_out = (a*x*phi_grid_sin(p_out)-b*y*phi_grid_cos(p_out)) / &
                  (b*n(2)*phi_grid_cos(p_out)-a*n(1)*phi_grid_sin(p_out))

          end if

       end if

       if (sol_p_in.gt.grid_min.and.sol_p_in.lt.1.e20_dp.and.sol_p_in.lt.face_distance) then
          face_distance = sol_p_in
          next_face(1) = 3
          next_face(2) = p_in
       end if

       if (sol_p_in.gt.grid_min.and.sol_p_in.lt.1.e20_dp.and.sol_p_in.lt.face_distance) then
          face_distance = sol_p_in
          next_face(1) = 3
          next_face(2) = p_in
       end if

       if (sol_p_out.gt.grid_min.and.sol_p_out.lt.1.e20_dp.and.sol_p_out.lt.face_distance) then
          face_distance = sol_p_out
          next_face(1) = 3
          next_face(2) = p_out
       end if

       if (sol_p_out.gt.grid_min.and.sol_p_out.lt.1.e20_dp.and.sol_p_out.lt.face_distance) then
          face_distance = sol_p_out
          next_face(1) = 3
          next_face(2) = p_out
       end if

    end if

    if (next_face(1).eq.0) then

       if (debug) then
          open (11, file=trim(error_log), position="append")
          write (11,*) "error 031"
          write (11,*) photon
          write (11,*) sol_r_in
          write (11,*) sol_r_out
          write (11,*) sol_r_same
          write (11,*) sol_t_in
          write (11,*) sol_t_out
          write (11,*) sol_t_same
          write (11,*) sol_p_in
          write (11,*) sol_p_out
          write (11,*) face
          write (11,*) cell
          write (11,*) next_face
          write (11,*) n
          write (11,*) x,y,z, sqrt(x*x+y*y+z*z)
          close (11)
          call exit(0)
       end if
       cell_error = .true.

    else

       call next_cell(face, next_face, cell, cell_out)
       
    end if

    if (next_face(1).eq.1.and.next_face(2).eq.nr) grid_exit = .true.

    if (face(1).eq.1.and.face(2).eq.cell_depth.and.&
         next_face(1).eq.1.and.next_face(2).eq.cell_depth) then

       if (debug) then
          open (11, file=trim(error_log), position="append")
          write (11,*) "error 034"
          close (11)
       end if
       cell_error = .true.

    else if (cell_out(1).eq.nr.and..not.grid_exit) then

       if (debug) then
          open (11, file=trim(error_log), position="append")
          write (11,*) "error 035"
          close (11)
       end if
       cell_error = .true.

    else if (cell_out(2).eq.ntheta) then

       if (debug) then
          open (11, file=trim(error_log), position="append")
          write (11,*) "error 036"
          close (11)
       end if
       cell_error = .true.

    else if (cell(1).eq.cell_out(1).and.cell(2).eq.cell_out(2).and.&
         cell(3).eq.cell_out(3).and..not.grid_exit) then

       if (debug) then
          open (11, file=trim(error_log), position="append")
          write (11,*) "error 037"
          close (11)
       end if
       cell_error = .true.

    end if

    if (cell_error) then

       if (debug) then
          open (11, file=trim(error_log), position="append")
          write (11,*) "error 058"
          close (11)
       end if

    end if

  end subroutine cell_face

  subroutine write_output

    integer               :: i, j, k, m, n
    real(dp)              :: planck_flux, error(npix,npix,5), norm, e_pack, dummy, pol, dpol
    real(dp), allocatable :: transport(:,:,:,:)
    logical               :: exist
    
    error = 0._dp

    ! Stokes error
    
    do k=1,4
       do j=1,npix
          do i=1,npix

             ! Error: stdev = sqrt( sum_x2/n - mean^2 )
             ! mean = sum_x / n

             if (detector(i,j,k,3).gt.0._dp) then

                dummy = (detector(i,j,k,2)/detector(i,j,k,3)) - (detector(i,j,k,1)/detector(i,j,k,3))**2
                if (dummy.gt.0._dp) error(i,j,k) = sqrt(dummy) * sqrt(detector(i,j,k,3))

             end if

          end do
       end do
    end do

    ! Degree of polarization error
    
    do j=1,npix
       do i=1,npix

          if (detector(i,j,2,1)**2+detector(i,j,3,1)**2.gt.0._dp) then
          
             pol = sqrt(detector(i,j,2,1)**2+detector(i,j,3,1)**2)

             dpol = sqrt( ( (detector(i,j,2,1)*error(i,j,2))**2 + (detector(i,j,3,1)*error(i,j,3))**2 ) / &
                  ( 2._dp*(detector(i,j,2,1)**2+detector(i,j,3,1)**2) ) )

          end if

          if (detector(i,j,1,1).gt.0._dp) error(i,j,5) = &
               (pol/detector(i,j,1,1)) * sqrt( (dpol/pol)**2 + (error(i,j,1)/detector(i,j,1,1))**2  )
          
       end do
    end do

    if (det_type.eq.3) then

       ! Phase curve

       inquire(file=trim(output_name)//'/output/phase.dat', exist=exist)

       if (exist) then

          open (100, file=trim(output_name)//'/output/phase.dat', status='old', position='append', action='write')

       else if (.not.exist) then

          open (100, file=trim(output_name)//'/output/phase.dat', status='new', action='write')
          write (100,*) "# Phase angle [deg] - Stokes I, Q, U, V [W m-2 micron-1]"
          write (100,*)

       end if

       ! Phase [deg] - Stokes I, I error, Q, Q error, U, U error, V, V error [W m-2 micron-1]
       
       if (det_phi*180._dp/pi.lt.1._dp) then

          write (100,*) 0._dp, detector(1,1,1,1)*1.e-6_dp, error(1,1,1)*1.e-6_dp, &
               detector(1,1,2,1)*1.e-6_dp, error(1,1,2)*1.e-6_dp, detector(1,1,3,1)*1.e-6_dp, &
               error(1,1,3)*1.e-6_dp, detector(1,1,4,1)*1.e-6_dp, error(1,1,4)*1.e-6_dp

       else if (det_phi*180._dp/pi.gt.179_dp) then

          write (100,*) 180._dp, detector(1,1,1,1)*1.e-6_dp, error(1,1,1)*1.e-6_dp, &
               detector(1,1,2,1)*1.e-6_dp, error(1,1,2)*1.e-6_dp, detector(1,1,3,1)*1.e-6_dp, &
               error(1,1,3)*1.e-6_dp, detector(1,1,4,1)*1.e-6_dp, error(1,1,4)*1.e-6_dp
          
       else
          
          write (100,*) det_phi*180._dp/pi, detector(1,1,1,1)*1.e-6_dp, error(1,1,1)*1.e-6_dp, &
               detector(1,1,2,1)*1.e-6_dp, error(1,1,2)*1.e-6_dp, detector(1,1,3,1)*1.e-6_dp, &
               error(1,1,3)*1.e-6_dp, detector(1,1,4,1)*1.e-6_dp, error(1,1,4)*1.e-6_dp

       end if

       close (100)

    else if (det_type.eq.1) then
          
       ! Write FITS images

       call write_fits_3D('stokes.fits', detector(:,:,:,1)*1.e-6_dp/(pixel_scale*pixel_scale), npix, npix, 4)
       call write_fits_3D('error.fits', error, npix, npix, 5)
       
       ! Write photometry

       if (det_type.eq.1) then

          open (100, file=trim(output_name)//'/output/photometry.dat', status='new', action='write')

          write (100,*) "# Wavelength [micron] - Stokes I, Q, U, V [W m-2 micron-1]"
          write (100,*)

          ! Wavelength [micron] - Stokes I, I error, Q, Q error, U, U error, V, V error [W m-2 micron-1]
          write (100,*) wavelength*1.e6_dp, 1.e-6_dp*photometry(1), 1.e-6_dp*photometry(2), 1.e-6_dp*photometry(3), &
               1.e-6_dp*photometry(4), 1.e-6_dp*photometry(5), 1.e-6_dp*photometry(6), &
               1.e-6_dp*photometry(7), 1.e-6_dp*photometry(8)

          close (100)

       end if
       
    else if (det_type.eq.2) then

       ! Spectrum
       
       inquire(file=trim(output_name)//'/output/spectrum.dat', exist=exist)

       if (exist) then

          open (100, file=trim(output_name)//'/output/spectrum.dat', status='old', position='append', action='write')

       else if (.not.exist) then

          open (100, file=trim(output_name)//'/output/spectrum.dat', status='new', action='write')

          write (100,*) "# Wavelength [micron] - Stokes I, Q, U, V [W m-2 micron-1]"
          write (100,*)

       end if

       ! Wavelength [micron] - Stokes I, Q, U, V [W m-2 micron-1]
       write (100,*) wavelength*1.e6_dp, 1.e-6_dp*detector(1,1,1,1), 1.e-6_dp*detector(1,1,2,1), &
            1.e-6_dp*detector(1,1,3,1), 1.e-6_dp*detector(1,1,4,1)

       close (100)

    end if

    if (photon_source.eq.1) then

       ! Flux normalization constants

       if ((det_type.eq.3.and.det_phi.lt.pi/180._dp).or.det_type.eq.1.or.det_type.eq.2) then

          inquire(file=trim(output_name)//'/output/normalization.dat', exist=exist)

          if (exist) then

             open (100, file=trim(output_name)//'/output/normalization.dat', status='old', position='append', action='write')

          else if (.not.exist) then

             open (100, file=trim(output_name)//'/output/normalization.dat', status='new', action='write')
             write (100,*) "# Wavelength [micron] - Norm constant 1 [W m-2 micron-1] - Norm constant 2 [W m-2 micron-1]"
             write (100,*)

          end if

          call planck_function(t_star, planck_flux)

          ! Wavelength [micron] - Norm1 [W m-2 micron-1] - Norm2 [W m-2 micron-1]
          write (100,*) wavelength*1.e6_dp, 1.e-6_dp*planck_flux*r_star*r_star / &
               (distance_planet*distance_planet), 1.e-6_dp*planck_flux*rfront(nr)*rfront(nr)*r_star*r_star / &
               (orbit*orbit*distance_planet*distance_planet)

          close (100)

       end if

    else if (photon_source.eq.2) then

       ! Write cell luminosity
       
       if (det_type.eq.1) call write_fits_3D("cell_luminosity.fits", cell_luminosity, nr, ntheta, nphi)

       ! Write emitted and emergent luminosity
       
       e_pack = emissivity_cumulative(nr-1,ntheta-1,nphi-1)/dble(packages)

       inquire(file=trim(output_name)//'/output/luminosity.dat', exist=exist)

       if (exist) then
          
          open (100, file=trim(output_name)//'/output/luminosity.dat', status='old', position='append', action='write')
          
       else if (.not.exist) then

          open (100, file=trim(output_name)//'/output/luminosity.dat', status='new', action='write')

          write (100,*) "# Wavelength [deg] - Emitted luminosity [W micron-1] - &
               & Emergent luminosity [W micron-1] - Emergent luminosity [a.u.]"
          write (100,*)
          
       end if

       write (100,*) wavelength, sum(flux_emitted)*e_pack*1.e-6_dp, sum(flux_exit)*e_pack*1.e-6_dp, sum(flux_exit)

       close (100)
       
    end if

    ! Write cell depth

    if (det_type.eq.1.or.det_type.eq.2) then
    
       inquire(file=trim(output_name)//'/output/cell_depth.dat', exist=exist)

       if (exist) then
          
          open (100, file=trim(output_name)//'/output/cell_depth.dat', status='old', position='append', action='write')
          
       else if (.not.exist) then

          open (100, file=trim(output_name)//'/output/cell_depth.dat', status='new', action='write')

          write (100,*) "# Wavelength [micron] - Cell depth"
          write (100,*)
          
       end if

       write (100,*) wavelength*1.e6_dp, cell_depth

       close (100)

    end if
    
    ! Write flow output

    if (global) then
    
       allocate (transport(3,0:nr-1,0:ntheta-1,0:nphi-1))
       transport = 0._dp

       do k=0,nphi-1
          do j=0,ntheta-1
             do i=cell_depth,nr-1
                do m=1,3
                   do n=1,threads
                      transport(m,i,j,k) = transport(m,i,j,k) + cell_flow_global(n,m,i,j,k)
                   end do
                end do
                norm = sqrt(transport(1,i,j,k)**2+transport(2,i,j,k)**2+transport(3,i,j,k)**2)
                if (norm.gt.0._dp) then
                   transport(1,i,j,k) = transport(1,i,j,k)/norm
                   transport(2,i,j,k) = transport(2,i,j,k)/norm
                   transport(3,i,j,k) = transport(3,i,j,k)/norm
                end if
             end do
          end do
       end do

       call write_fits_4D("global.fits", transport, 3, nr, ntheta, nphi)

       deallocate (transport)

    end if

    if (latitudinal) then

       ! Write latitudinal energy transport, normalized to the total emergent flux

       allocate (transport(4,0:nr-1,0:ntheta-1,0:nphi-1))
       transport = 0._dp

       do k=0,nphi-1
          do j=0,ntheta-1
             do i=cell_depth,nr-1

                transport(1,i,j,k) = sum(cell_flow(:,1,i,j,k))
                transport(2,i,j,k) = sum(cell_flow(:,2,i,j,k))
                transport(3,i,j,k) = sum(cell_flow(:,3,i,j,k))
                transport(4,i,j,k) = sum(cell_flow(:,4,i,j,k))
                
             end do
          end do
       end do
       
       if (sum(flux_exit).gt.0._dp) then
          transport = transport / sum(flux_exit)
          call write_fits_4D("latitudinal.fits", transport, 4, nr, ntheta, nphi)
       end if

       deallocate (transport)

    end if

  end subroutine write_output

  subroutine write_fits_3D(filename, array, n1, n2, n3)

    ! Write 3D array to FITS file

    integer,      intent(in) :: n1, n2, n3
    integer                  :: status, unit, blocksize, bitpix, naxis, naxes(3), group, fpixel
    real(dp),     intent(in) :: array(n1,n2,n3)
    character(*), intent(in) :: filename
    logical                  :: simple, extend
    character(300)           :: fits_path

    status = 0
    blocksize = 1
    simple = .true.
    naxis = 3
    naxes(1) = n1
    naxes(2) = n2
    naxes(3) = n3
    extend = .true.
    group = 1
    fpixel = 1
    bitpix = -64

    write (fits_path, "((a),'/output/',(a))") trim(output_name), trim(filename)

    call ftgiou(unit, status)
    call ftinit(unit, trim(fits_path), blocksize, status)
    call ftphpr(unit, simple, bitpix, naxis, naxes, 0, 1, extend, status)
    call ftpprd(unit, group, fpixel, n1*n2*n3, array, status)
    call ftclos(unit, status)
    call ftfiou(unit, status)

  end subroutine write_fits_3D

  subroutine write_fits_4D(filename, array, n1, n2, n3, n4)

    ! Write 4D array to FITS file

    integer,      intent(in) :: n1, n2, n3, n4
    integer                  :: status, unit, blocksize, bitpix, naxis, naxes(4), group, fpixel
    real(dp),     intent(in) :: array(n1,n2,n3,n4)
    character(*), intent(in) :: filename
    logical                  :: simple, extend
    character(300)           :: fits_path

    status = 0
    blocksize = 1
    simple = .true.
    naxis = 4
    naxes(1) = n1
    naxes(2) = n2
    naxes(3) = n3
    naxes(4) = n4
    extend = .true.
    group = 1
    fpixel = 1
    bitpix = -64

    write (fits_path, "((a),'/output/',(a))") trim(output_name), trim(filename)

    call ftgiou(unit, status)
    call ftinit(unit, trim(fits_path), blocksize, status)
    call ftphpr(unit, simple, bitpix, naxis, naxes, 0, 1, extend, status)
    call ftpprd(unit, group, fpixel, n1*n2*n3*n4, array, status)
    call ftclos(unit, status)
    call ftfiou(unit, status)

  end subroutine write_fits_4D
  
  subroutine output(log_number)

    integer, intent(in) :: log_number
    integer             :: i, j, k, hours, minutes, seconds, log_size
    real(dp)            :: angular_size, total_optical_depth, cpu_total, e_pack, norm, norm2, planck_flux
    logical             :: file_exists
    character(200)      :: hostname

    if (log_number.eq.1) then

       angular_size = 2._dp*atan(rfront(nr)/distance_planet)*3600._dp*180._dp/pi*1000._dp

       write (out_unit,'(a)') "########################################################"
       write (out_unit,'(a)') "              _     ___   _____   ___   ___             "
       write (out_unit,'(a)') "             /_\   | _ \ |_   _| | __| / __|            "
       write (out_unit,'(a)') "            / _ \  |   /   | |   | _|  \__ \            "
       write (out_unit,'(a)') "           /_/ \_\ |_|_\   |_|   |___| |___/            "
       write (out_unit,'(a)') "                                                        "
       write (out_unit,'(a)') "  Atmospheric Radiative Transfer for Exoplanet Science  "
       write (out_unit,'(a)') "                                                        "
       write (out_unit,'(a)') "    Please cite Stolker et al., 2017, A&A, 607, A42     "
       write (out_unit,'(a)') " whenever results from ARTES are used in a publication. "
       write (out_unit,'(a)') ""
       write (out_unit,'(a)') "--------------------------------------------------------"
       write (out_unit,'(a)') ""
       write (out_unit,'(a)') "--> Initialization"
       write (out_unit,'(a)') ""
       write (out_unit,'(a,a)') "Input file: ", input_file
       write (out_unit,'(a,i0)') "Computer cores: ", cores
       write (out_unit,'(a,i0)') "Threads in use: ", threads
       write (out_unit,'(a)') ""
       write (out_unit,'(a)') "--------------------------------------------------------"
       write (out_unit,'(a)') ""
       write (out_unit,'(a)') "--> Build planet atmosphere"
       write (out_unit,'(a)') ""
       write (out_unit,'(a,es8.2)') "Planet radius [km]: ", rfront(nr)/1000._dp
       write (out_unit,'(a,es8.2)') "Atmosphere height [km]: ", (rfront(nr)-rfront(0))/1000._dp
       write (out_unit,'(a)') ""
       write (out_unit,'(a,es8.2)') "Oblateness: ", oblateness
       if (oblateness.gt.0._dp) write (out_unit,'(a,es8.2)') "Eccentricity: ", &
            sqrt( 1._dp - (oblate_z*oblate_z) / (oblate_x*oblate_x) )
       write (out_unit,'(a)') ""
       write (out_unit,'(a,es8.2)') "Surface albedo: ", surface_albedo
       write (out_unit,'(a)') ""
       write (out_unit,'(a,i0)') "Radial grid cells: ", nr
       write (out_unit,'(a,i0)') "Latitudinal grid cells: ", ntheta
       write (out_unit,'(a,i0)') "Longitudial grid cells: ", nphi
       write (out_unit,'(a)') ""
       write (out_unit,'(a,es8.2,a,es8.2)') "Field of view [mas x mas]: ", x_fov, " x ", y_fov
       write (out_unit,'(a,es8.2,a,es8.2)') "Planet angular diameter [mas]: ", angular_size
       write (out_unit,'(a,es8.2,a,es8.2)') "Pixel scale [mas pixel-1]: ", pixel_scale
       write (out_unit,'(a)') ""

    else if (log_number.eq.2) then

       write (out_unit,'(a)') "--------------------------------------------------------"
       write (out_unit,'(a)') ""
       write (out_unit,'(a)') "--> Photon transfer"
       write (out_unit,'(a)') ""
       if (photon_source.eq.1) then
          write (out_unit,'(a,a)') "Photon source: star"
       else if (photon_source.eq.2) then
          write (out_unit,'(a,a)') "Photon source: planet"
       end if
       write (out_unit,'(a,es8.2)') "Emitted photons: ", dble(packages)
       write (out_unit,'(a)') ""
       if (photon_source.eq.1.and.(det_type.eq.1.or.det_type.eq.2)) then
          write (out_unit,'(a,es8.2)') "Phase angle [deg]: ", phase_observer
          write (out_unit,'(a)') ""
       end if
       if (det_type.eq.1.or.det_type.eq.3) then
          write (out_unit,'(a,es8.2)') "Wavelength [micron]: ", wavelength*1.e6_dp
       else if (det_type.eq.2) then
          write (out_unit,'(a,es8.2,a,es8.2)') "Spectrum range [micron]: ", wavelengths(1)*1.e6, &
               " - ", wavelengths(size(wavelengths))*1.e6
       end if
       write (out_unit,'(a)') ""
       write (out_unit,'(a,es8.2)') "Stellar luminosity [W]: ", 4._dp*pi*r_star*r_star*sb*t_star**4
       write (out_unit,'(a)') ""
       if (photon_walk.lt.0._dp) then
          write (out_unit,'(a)') "Modified random walk: off"
       else if (photon_walk.gt.0._dp) then
          write (out_unit,'(a)') "Modified random walk: on"
       end if
       write (out_unit,'(a)') ""
       if (det_type.eq.1.or.det_type.eq.3) then
          write (out_unit,'(a)') "Total optical depth:"
          do i=0,ntheta-1
             do j=0,nphi-1
                total_optical_depth = 0._dp
                do k=cell_depth,nr-1
                   total_optical_depth = total_optical_depth + ( (rfront(k+1)-rfront(k)) * &
                        cell_opacity(k,i,j,wl_count) )
                end do
                write (out_unit,'(a,i0,a,i0,a,es10.4)') "[Theta, phi] = [", i, ", ", j, "] --> ", total_optical_depth
             end do
          end do
          write (out_unit,'(a)') ""
          write (out_unit,'(a)') "Scattering optical depth:"
          do i=0,ntheta-1
             do j=0,nphi-1
                total_optical_depth = 0._dp
                do k=cell_depth,nr-1
                   total_optical_depth = total_optical_depth + ( (rfront(k+1)-rfront(k)) * &
                        cell_scattering_opacity(k,i,j,wl_count) )
                end do
                write (out_unit,'(a,i0,a,i0,a,es10.4)') "[Theta, phi] = [", i, ", ", j, "] --> ", total_optical_depth
             end do
          end do
          write (out_unit,'(a)') ""
          write (out_unit,'(a)') "Absorption optical depth:"
          do i=0,ntheta-1
             do j=0,nphi-1
                total_optical_depth = 0._dp
                do k=cell_depth,nr-1
                   total_optical_depth = total_optical_depth + ( (rfront(k+1)-rfront(k)) * &
                        cell_absorption_opacity(k,i,j,wl_count) )
                end do
                write (out_unit,'(a,i0,a,i0,a,es10.4)') "[Theta, phi] = [", i, ", ", j, "] --> ", total_optical_depth
             end do
          end do
          write (out_unit,'(a)') ""
       end if

    else if (log_number.eq.3) then

       write (out_unit,'(a)') "--------------------------------------------------------"
       write (out_unit,'(a)') ""

       if (det_type.eq.1) then

          if (photon_source.eq.1) then

             call planck_function(t_star, planck_flux)

             norm = planck_flux*rfront(nr)*rfront(nr)*r_star*r_star / (orbit*orbit*distance_planet*distance_planet)
             norm2 = planck_flux*r_star*r_star/(distance_planet*distance_planet)

             if (photometry(1).gt.0._dp) then

                write (out_unit,'(a)') "Planet integrated flux"
                write (out_unit,'(a)') ""
                write (out_unit,'(a,es9.2)') "Stokes I [W m-2 micron-1]: ", photometry(1)*1.e-6_dp
                write (out_unit,'(a,es9.2)') "Stokes Q [W m-2 micron-1]: ", photometry(3)*1.e-6_dp
                write (out_unit,'(a,es9.2)') "Stokes U [W m-2 micron-1]: ", photometry(5)*1.e-6_dp
                write (out_unit,'(a,es9.2)') "Stokes V [W m-2 micron-1]: ", photometry(7)*1.e-6_dp
                write (out_unit,'(a)') ""
                write (out_unit,'(a,es9.2)') "Normalized Stokes I: ", photometry(1)/norm
                write (out_unit,'(a,es9.2)') "Normalized Stokes Q: ", photometry(3)/norm
                write (out_unit,'(a,es9.2)') "Normalized Stokes U: ", photometry(5)/norm
                write (out_unit,'(a,es9.2)') "Normalized Stokes V: ", photometry(7)/norm
                write (out_unit,'(a)') ""
                write (out_unit,'(a,es9.2)') "Stellar normalized Stokes I: ", photometry(1)/norm2
                write (out_unit,'(a,es9.2)') "Stellar normalized Stokes Q: ", photometry(3)/norm2
                write (out_unit,'(a,es9.2)') "Stellar normalized Stokes U: ", photometry(5)/norm2
                write (out_unit,'(a,es9.2)') "Stellar normalized Stokes V: ", photometry(7)/norm2
                write (out_unit,'(a)') ""
                write (out_unit,'(a,es9.2)') "-Q/I: ", -photometry(3)/photometry(1)
                write (out_unit,'(a,es9.2)') " U/I: ", photometry(5)/photometry(1)
                write (out_unit,'(a,es9.2)') " V/I: ", photometry(7)/photometry(1)
                write (out_unit,'(a)') ""
                write (out_unit,'(a,es8.2,a,es8.2)') "Degree of polarization [%]: ", &
                     100._dp * photometry(10), " +/- ", 100._dp * photometry(11)
                write (out_unit,'(a,es9.2)') "Direction of polarization [deg]: ", &
                     0.5*atan2(photometry(5),photometry(3)) * 180._dp/pi
                write (out_unit,'(a)') ""

             else

                write (out_unit,'(a)') "Error in Stokes I"
                write (out_unit,'(a)') ""

             end if

             write (out_unit,'(a)') "--------------------------------------------------------"
             write (out_unit,'(a)') ""

          else if (photon_source.eq.2) then

             if (photometry(1).gt.0._dp) then

                write (out_unit,'(a)') "Planet integrated flux"
                write (out_unit,'(a)') ""
                write (out_unit,'(a,es9.2)') "Stokes I [W m-2 micron-1]: ", photometry(1)*1.e-6_dp
                write (out_unit,'(a,es9.2)') "Stokes Q [W m-2 micron-1]: ", photometry(3)*1.e-6_dp
                write (out_unit,'(a,es9.2)') "Stokes U [W m-2 micron-1]: ", photometry(5)*1.e-6_dp
                write (out_unit,'(a,es9.2)') "Stokes V [W m-2 micron-1]: ", photometry(7)*1.e-6_dp
                write (out_unit,'(a)') ""
                write (out_unit,'(a,es9.2)') "-Q/I: ", -photometry(3)/photometry(1)
                write (out_unit,'(a,es9.2)') " U/I: ", photometry(5)/photometry(1)
                write (out_unit,'(a,es9.2)') " V/I: ", photometry(7)/photometry(1)
                write (out_unit,'(a)') ""
                write (out_unit,'(a,es8.2,a,es8.2)') "Degree of polarization [%]: ", &
                     100._dp * photometry(10), " +/- ", 100._dp * photometry(11)
                write (out_unit,'(a,es9.2)') "Direction of polarization [deg]: ", &
                     0.5*atan2(photometry(5),photometry(3)) * 180._dp/pi
                write (out_unit,'(a)') ""

             else

                write (out_unit,'(a)') "Error: Stokes I is zero"
                write (out_unit,'(a)') ""

             end if

             write (out_unit,'(a)') "--------------------------------------------------------"
             write (out_unit,'(a)') ""

             if (det_type.eq.1.and.photon_source.eq.2) then

                e_pack = emissivity_cumulative(nr-1,ntheta-1,nphi-1)/dble(packages)
                
                write (out_unit,'(a,es9.2)') "Total emitted flux [W micron-1] =", sum(flux_emitted)*e_pack*1.e-6_dp
                write (out_unit,'(a,es9.2)') "Total emergent flux [W micron-1] =", sum(flux_exit)*e_pack*1.e-6_dp
                write (out_unit,'(a)') ""
                write (out_unit,'(a)') "--------------------------------------------------------"
                write (out_unit,'(a)') ""

             end if

          end if

       end if

    else if (log_number.eq.4) then

       t_end = omp_get_wtime()
       cpu_total = t_end - t_start

       hours   = int( cpu_total / 3600._dp )
       minutes = int( cpu_total / 60._dp - dble(hours)*60._dp )
       seconds = int( cpu_total - dble(hours)*3600._dp - dble(minutes)*60._dp )

       write (out_unit,'(a,i0.2,a,i0.2,a,i0.2)') "CPU time [hour:min:sec]: ", hours, ":", minutes, ":", seconds
       write (out_unit,'(a)') ""

       if (debug) then
          inquire (file=trim(error_log),size=log_size)
          if (log_size.ne.0) then
             write (out_unit,'(a)') "WARNING: check error log!"
             write (out_unit,'(a)') ""
          end if
       end if
       
       write (out_unit,'(a)') "########################################################"

       if (log_file) close (10)

       if (len(trim(email)).ne.0) then

          call system('hostname > hostname.txt')
          open (100, file='hostname.txt')
          read (100,*) hostname
          close (100, status='delete')

          inquire(file="/usr/bin/mail", exist=file_exists)

          if (file_exists) then

             open (100, file='mail.txt')
             write (100,'(a,a,a,a,a)') "Job with input ", adjustl(trim(input_file)), " on ", &
                  adjustl(trim(hostname)), " is finished."
             write (100,'(a)')
             write (100,'(a)') "Have a nice day!"

             call system('mail -s "ARTES is finished" ' // adjustl(trim(email)) // ' < mail.txt')

             close (100, status='delete')

          else if (.not.file_exists) then

             inquire(file="/usr/sbin/ssmtp", exist=file_exists)

             if (file_exists) then

                open (100, file='mail.txt')

                write (100,'(a,a)') "To:" // adjustl(trim(email))
                write (100,'(a)') "From:ARTES"
                write (100,'(a)') "Subject: ARTES is finished"
                write (100,'(a)') ""
                write (100,'(a,a,a,a,a)') "Job with input ", adjustl(trim(input_file)), " on ", &
                     adjustl(trim(hostname)), " is finished."
                write (100,'(a)')
                write (100,'(a)') "Have a nice day!"

                call system('ssmtp ' // adjustl(trim(email)) // ' < mail.txt')

                close (100, status='delete')

             else

                if (debug) then
                   open (11, file=trim(error_log), position="append")
                   write (11,*) "error 038"
                   close (11)
                end if

             end if

          end if

       end if

       close (11)

    end if

  end subroutine output

  subroutine quadratic_equation(a, b, c, solutions)

    ! Solve a quadratic equation

    real(dp), intent(in)  :: a, b, c
    real(dp), intent(out) :: solutions(2)
    real(dp)              :: discriminant, q

    solutions = 0._dp
    discriminant = b*b - 4._dp*a*c

    if (discriminant.ge.0._dp) then

       q = -0.5_dp * ( b + sign(1._dp,b) * sqrt(discriminant)  )
       if (abs(a).gt.1.e-100_dp) solutions(1) = q/a
       if (abs(q).gt.1.e-100_dp) solutions(2) = c/q

    end if

  end subroutine quadratic_equation

  subroutine init_random_seed

    ! Initialize a random seed for the random number generator

    integer               :: clock, i, n
    integer, allocatable  :: seed(:)

    call random_seed(size = n)

    allocate (seed(n))
    seed = 0

    call system_clock(count=clock)

    seed = clock + 37 * (/ (i - 1, i = 1, n) /)

    call random_seed(put=seed)

    deallocate (seed)

  end subroutine init_random_seed

  subroutine random(thread_id, xi)

    ! Marsaglia & Zaman (1994)

    integer,  intent(in)  :: thread_id
    integer               :: imz
    real(dp), intent(out) :: xi

    imz = state(thread_id+1,1) - state(thread_id+1,3)

    if (imz.lt.0) imz = imz + 2147483579

    state(thread_id+1,1) = state(thread_id+1,2)
    state(thread_id+1,2) = state(thread_id+1,3)
    state(thread_id+1,3) = imz
    state(thread_id+1,4) = 69069 * state(thread_id+1,4) + 1013904243

    imz = imz + state(thread_id+1,4)

    xi = 0.5_dp + 0.23283064e-9_dp * imz

    if (xi.le.0._dp.or.xi.ge.1._dp) then
       if (debug) then
          open (11, file=trim(error_log), position="append")
          write (11,*) "error 055"
          close (11)
       end if
    end if
       
  end subroutine random

  subroutine argument_input(input_file)

    ! Check the command line arguments

    character(100), intent(out) :: input_file
    character(100)              :: dummy_char
    real(dp)                    :: dummy_dble

    input_file = ""

    if (iargc().le.3) then
       
       write (6,'(a)') "How to run ARTES:"
       write (6,'(a)') "./path/to/bin/artes [input_file] [photons] -o [output_folder] -k [keyword]=[value]"
       call exit(0)
    
    else

       call getarg(1, input_file)
       call getarg(2, dummy_char)

       read (dummy_char,*) dummy_dble
       packages = int(dummy_dble,kind=16)

    end if
    
  end subroutine argument_input

  subroutine argument_keywords

    integer        :: i
    character(100) :: arg, keyword_dummy, key_word, key_value
    logical        :: exists

    do i = 1, command_argument_count()

       call get_command_argument(i, arg)

       select case (arg)
       case ('-o')
          
          call get_command_argument(i+1, output_name)

          ! Make output directories
          call system('rm -rf ' // trim(output_name))
          call system('mkdir -p ' // trim(output_name) )
          call system('mkdir -p ' // trim(output_name) // '/input' )
          call system('mkdir -p ' // trim(output_name) // '/output' )
          call system('mkdir -p ' // trim(output_name) // '/plot' )

          ! Copy input files to output directory
          call system('cp '//adjustl(trim(input_file))//' '//adjustl(trim(output_name))//'/input/')
          call system('cp '//trim(input_file(1:len_trim(input_file)-8)//'atmosphere.fits')// &
               ' '//adjustl(trim(output_name))//'/input/')
          inquire(file=trim(input_file(1:len_trim(input_file)-8))//'atmosphere.dat', exist=exists)
          if (exists) call system('cp '//trim(input_file(1:len_trim(input_file)-8)//'atmosphere.dat')// &
               ' '//adjustl(trim(output_name))//'/input/')
          inquire(file=trim(input_file(1:len_trim(input_file)-8))//'pressure_temperature.dat', exist=exists)
          if (exists) call system('cp '//trim(input_file(1:len_trim(input_file)-8)//'pressure_temperature.dat')// &
               ' '//adjustl(trim(output_name))//'/input/')
          
       case ('-k')
          
          call get_command_argument(i+1, keyword_dummy)
          call get_key_value(keyword_dummy, key_word, key_value)
          call input_parameters(key_word, key_value)

          open (100, file=trim(output_name)//'/input/'//adjustl(trim(input_file)), status='old', position='append', action='write')
          write (100,*) new_line('a')//trim(keyword_dummy)
          close (100)
          
       end select
       
    end do

  end subroutine argument_keywords
  
  subroutine readfile(input_file, data, nlines)

    ! Check number of lines in a data file and read into array

    integer,                     intent(out) :: nlines
    integer                                  :: j, line_counter, ios
    integer, parameter                       :: max_lines = 10000
    character(100),              intent(in)  :: input_file 
    character(500), allocatable, intent(out) :: data(:)
    character(100)                           :: junk

    line_counter = 0

    open (100,file=input_file)

    do j=1,max_lines

       read (100,'(a500)',iostat=ios) junk

       if (ios.ne.0) exit

       if (j.eq.max_lines) then

          if (debug) then
             open (11, file=trim(error_log), position="append")
             write (11,*) "error 039"
             close (11)
          end if
          stop

       end if

       line_counter = line_counter + 1

    end do

    rewind(100)

    allocate (data(line_counter))

    do j=1,line_counter

       read (100,'(a500)') data(j)

    enddo

    nlines = line_counter

    close (100)

  end subroutine readfile

  subroutine input_parameters(key_word, key_value)

    integer                     :: key_value_int
    real(dp)                    :: key_value_double
    character(100), intent(in)  :: key_word, key_value
    character(100)              :: key_value_char

    key_value_int = 0
    key_value_double = 0._dp
    key_value_char = ""

    select case(key_word)

    case("general:log")
       if (key_value.eq."on") then
          log_file = .true.
          out_unit = 10
       else if (key_value.eq."off") then
          log_file = .false.
          out_unit = 6
       end if
    case("general:email")
       read (key_value, '(a100)') key_value_char
       email = key_value_char
    case("general:random")
       if (key_value.eq."on") then
          ranseed = .true.
       else if (key_value.eq."off") then
          ranseed = .false.
       end if
    case("photon:source")
       if (key_value.eq."star") then
          photon_source = 1
       else if (key_value.eq."planet") then
          photon_source = 2
       end if
    case("photon:fstop")
       read (key_value, '(e50.0)') key_value_double
       fstop = key_value_double
    case("photon:minimum")
       read (key_value, '(e50.0)') key_value_double
       photon_minimum = key_value_double
    case("photon:weight")
       if (key_value.eq."on") then
          thermal_weight = .true.
       else if (key_value.eq."off") then
          thermal_weight = .false.
       end if
    case("photon:scattering")
       if (key_value.eq."on") then
          photon_scattering = .true.
       else if (key_value.eq."off") then
          photon_scattering = .false.
       end if
    case("photon:emission")
       if (key_value.eq."isotropic") then
          photon_emission = 1
       else if (key_value.eq."biased") then
          photon_emission = 2
       end if
    case("photon:bias")
       read (key_value, '(e50.0)') key_value_double
       photon_bias = key_value_double
    case("photon:walk")
       read (key_value, '(e50.0)') key_value_double
       photon_walk = key_value_double
    case("photon:number")
       read (key_value, '(i50)') key_value_int
       photon_number = key_value_int
    case("star:temperature")
       read (key_value, '(e50.0)') key_value_double
       t_star = key_value_double
    case("star:radius")
       read (key_value, '(e50.0)') key_value_double
       r_star = key_value_double*r_sun
    case("star:theta")
       read (key_value, '(e50.0)') key_value_double
       star_theta = key_value_double*pi/180._dp
       if (star_theta.lt.1.e-6) star_theta = 1.e-6_dp
       if (star_theta.gt.pi-1.e-6) star_theta = pi-1.e-6
    case("star:phi")
       read (key_value, '(e50.0)') key_value_double
       star_phi = key_value_double*pi/180._dp
       if (star_phi.lt.0._dp) star_phi = star_phi + 2._dp*pi
       if (star_phi.ge.2._dp*pi) star_phi = star_phi - 2._dp*pi
    case("planet:albedo")
       read (key_value, '(e50.0)') key_value_double
       surface_albedo = key_value_double
    case("planet:oblateness")
       read (key_value, '(e50.0)') key_value_double
       oblateness = key_value_double
    case("planet:orbit")
       read (key_value, '(e50.0)') key_value_double
       orbit = key_value_double*au
    case("planet:ring")
       if (key_value.eq."on") then
          ring = .true.
       else if (key_value.eq."off") then
          ring = .false.
       end if
    case("planet:grid")
       read (key_value, '(e50.0)') key_value_double
       grid_min = key_value_double
    case("planet:tau")
       read (key_value, '(e50.0)') key_value_double
       tau_max = key_value_double
    case("detector:type")
       if (key_value.eq."imaging") then
          det_type = 1
       else if (key_value.eq."spectrum") then
          det_type = 2
       else if (key_value.eq."phase") then
          det_type = 3
       end if
    case("detector:theta")
       read (key_value, '(e50.0)') key_value_double
       det_theta = key_value_double*pi/180._dp
       if (det_theta.lt.1.e-3) det_theta = 1.e-3_dp
       if (det_theta.gt.pi-1.e-3) det_theta = pi-1.e-3
    case("detector:phi")
       read (key_value, '(e50.0)') key_value_double
       det_phi = key_value_double*pi/180._dp
       if (det_phi.lt.0._dp) det_phi = det_phi + 2._dp*pi
       if (det_phi.ge.2._dp*pi) det_phi = det_phi - 2._dp*pi
    case("detector:pixel")
       read (key_value, '(i50)') key_value_int
       npix = key_value_int
    case("detector:distance")
       read (key_value, '(e50.0)') key_value_double
       distance_planet = key_value_double*pc
    case("detector:angle")
       read (key_value, '(e50.0)') key_value_double
       det_angle = key_value_double*pi/180._dp
       if (det_angle.lt.0._dp) det_angle = det_angle + 2._dp*pi
       if (det_angle.ge.2._dp*pi) det_angle = det_angle - 2._dp*pi
    case("output:debug")
       if (key_value.eq."on") then
          debug = .true.
       else if (key_value.eq."off") then
          debug = .false.
       end if
    case("output:global")
       if (key_value.eq."on") then
          global = .true.
       else if (key_value.eq."off") then
          global = .false.
       end if
    case("output:latitudinal")
       if (key_value.eq."on") then
          latitudinal = .true.
       else if (key_value.eq."off") then
          latitudinal = .false.
       end if
    case default
       write (6,*) "Wrong keyword found in input file: ", key_word
       call exit(0)

    end select

  end subroutine input_parameters

  subroutine get_key_value(line,key_word,key_value)

    ! Seperate a key_word and key_value component of a string given key_word=key_value syntax

    character(100), intent(in)  :: line
    character(100), intent(out) :: key_word, key_value

    key_word = ""
    key_value = ""

    key_word=line(1:index(line,'=')-1)
    key_value=line(index(line,'=')+1:len_trim(line))

    if(key_value(1:1).eq.'"'.or.key_value(1:1).eq."'") key_value=key_value(2:len_trim(key_value)-1)

  end subroutine get_key_value

  subroutine peel_thermal(thread_id, photon, x_photon, y_photon, z_photon, stokes_peel_in, cell_in, face_in, cell_error)

    integer,  intent(in)    :: cell_in(3), face_in(2), thread_id
    integer(16), intent(in) :: photon
    integer                 :: cell(3), cell_out(3), face(2), next_face(2), ix, iy
    real(dp), intent(in)    :: x_photon, y_photon, z_photon, stokes_peel_in(4)
    real(dp)                :: x, y, z, face_distance, photon_weight, tau_cell, tau_total, x_im, y_im, det(3), x_new, y_new
    logical                 :: grid_exit
    logical, intent(out)    :: cell_error

    x = x_photon
    y = y_photon
    z = z_photon

    det(1) = det_dir(1)
    det(2) = det_dir(2)
    det(3) = det_dir(3)

    face = face_in
    cell = cell_in

    tau_total = 0._dp
    photon_weight = 0._dp

    do

       call cell_face(photon, x, y, z, det, face, next_face, face_distance, grid_exit, cell, cell_out, cell_error)
       
       tau_cell = face_distance*cell_opacity(cell(1), cell(2), cell(3), wl_count)
       tau_total = tau_total + tau_cell

       x = x + face_distance*det(1)
       y = y + face_distance*det(2)
       z = z + face_distance*det(3)

       if (cell_error.or.grid_exit.or.(next_face(1).eq.1.and.next_face(2).eq.cell_depth)) exit

       face = next_face
       cell = cell_out

    end do

    if (grid_exit.and.tau_total.lt.50._dp.and..not.cell_error) then

       photon_weight = exp(-tau_total) / (4._dp*pi)

       x_im = y_photon*cos_det_phi - x_photon*sin_det_phi
       y_im = z_photon*sin_det_theta - y_photon*cos_det_theta*sin_det_phi - x_photon*cos_det_theta*cos_det_phi
       
       if (det_angle.gt.0._dp) then

          x_new = x_im*cos_det_angle - y_im*sin_det_angle
          y_new = x_im*sin_det_angle + y_im*cos_det_angle

          x_im = x_new
          y_im = y_new
          
       end if

       ix = int(npix*(x_im+x_max)/(2._dp*x_max))+1
       iy = int(npix*(y_im+y_max)/(2._dp*y_max))+1

       if (photon_weight*stokes_peel_in(1).gt.0._dp.and.photon_weight*stokes_peel_in(1).lt.1.e100_dp) then
       
          detector_thread(ix,iy,1,1,thread_id+1) = detector_thread(ix,iy,1,1,thread_id+1) + photon_weight*stokes_peel_in(1)
          detector_thread(ix,iy,1,2,thread_id+1) = detector_thread(ix,iy,1,2,thread_id+1) + (photon_weight*stokes_peel_in(1))**2
          detector_thread(ix,iy,1,3,thread_id+1) = detector_thread(ix,iy,1,3,thread_id+1) + 1._dp

       else

          if (debug) then
             open (11, file=trim(error_log), position="append")
             write (11,*) "error 051"
             close (11)
          end if

       end if
       
    end if

  end subroutine peel_thermal

  subroutine peel_surface(thread_id, photon, x_photon, y_photon, z_photon, stokes_peel_in, cell_in, face_in)

    integer,  intent(in)    :: cell_in(3), face_in(2), thread_id
    integer(16), intent(in) :: photon
    integer                 :: cell(3), cell_out(3), face(2), next_face(2), ix, iy
    real(dp), intent(in)    :: x_photon, y_photon, z_photon, stokes_peel_in(4)
    real(dp)                :: x, y, z, face_distance, photon_weight, surface_normal(3), det_sphere(3), normal_sphere(3)
    real(dp)                :: tau_cell, tau_total, x_im, y_im, det(3), cos_angle, norm, x_new, y_new
    logical                 :: grid_exit, cell_error

    ! Surface normal on spherical surface at point of reflection and scale for oblate surface
    surface_normal(1) = x_photon / ( oblate_x * oblate_x )
    surface_normal(2) = y_photon / ( oblate_y * oblate_y )
    surface_normal(3) = z_photon / ( oblate_z * oblate_z )

    ! Make the vector a unit vector
    norm = sqrt(surface_normal(1)*surface_normal(1)+surface_normal(2)*surface_normal(2)+surface_normal(3)*surface_normal(3))
    surface_normal = surface_normal / norm

    x = x_photon
    y = y_photon
    z = z_photon

    det(1) = det_dir(1)
    det(2) = det_dir(2)
    det(3) = det_dir(3)

    ! Determine angle between surface normal and detector direction

    call cartesian_spherical(det(1), det(2), det(3), det_sphere(1), det_sphere(2), det_sphere(3))
    call cartesian_spherical(surface_normal(1), surface_normal(2), surface_normal(3), &
         normal_sphere(1), normal_sphere(2), normal_sphere(3))

    cos_angle = sin(det_sphere(2))*cos(det_sphere(3))*sin(normal_sphere(2))*cos(normal_sphere(3)) + &
         sin(det_sphere(2))*sin(det_sphere(3))*sin(normal_sphere(2))*sin(normal_sphere(3)) + &
         cos(det_sphere(2))*cos(normal_sphere(2))

    ! Angle should be smaller than 90 deg otherwise photon is directed into the surface

    if (cos_angle.gt.0._dp) then

       face = face_in

       cell(1) = cell_in(1) + 1
       cell(2) = cell_in(2)
       cell(3) = cell_in(3)

       ! Calculate optical depth

       tau_total = 0._dp
       photon_weight = 0._dp

       do

          call cell_face(photon, x, y, z, det, face, next_face, face_distance, grid_exit, cell, cell_out, cell_error)

          tau_cell = face_distance*cell_opacity(cell(1), cell(2), cell(3), wl_count)
          tau_total = tau_total + tau_cell

          x = x + face_distance*det(1)
          y = y + face_distance*det(2)
          z = z + face_distance*det(3)

          if (cell_error.or.grid_exit.or.(next_face(1).eq.1.and.next_face(2).eq.cell_depth)) exit

          face = next_face
          cell = cell_out

       end do

       if (grid_exit.and.tau_total.lt.50._dp.and..not.cell_error) then

          ! Lambertian scattering probability: P(mu,phi) dmu dphi = 2mu dmu dphi/(2pi) -> P(mu,phi) = mu/pi

          ! Weight photon for optical depth and scattering angle

          photon_weight = exp(-tau_total) * cos_angle / pi

          x_im = y_photon*cos_det_phi - x_photon*sin_det_phi
          y_im = z_photon*sin_det_theta - y_photon*cos_det_theta*sin_det_phi - x_photon*cos_det_theta*cos_det_phi

          if (det_angle.gt.0._dp) then

             x_new = x_im*cos_det_angle - y_im*sin_det_angle
             y_new = x_im*sin_det_angle + y_im*cos_det_angle

             x_im = x_new
             y_im = y_new
          
          end if
          
          ix = int(npix*(x_im+x_max)/(2._dp*x_max))+1
          iy = int(npix*(y_im+y_max)/(2._dp*y_max))+1

          if (photon_weight*stokes_peel_in(1).gt.0._dp.and.photon_weight*stokes_peel_in(1).lt.1.e100_dp) then

             detector_thread(ix,iy,1,1,thread_id+1) = detector_thread(ix,iy,1,1,thread_id+1) + photon_weight*stokes_peel_in(1)
             detector_thread(ix,iy,1,2,thread_id+1) = detector_thread(ix,iy,1,2,thread_id+1) + (photon_weight*stokes_peel_in(1))**2
             detector_thread(ix,iy,1,3,thread_id+1) = detector_thread(ix,iy,1,3,thread_id+1) + 1._dp
             
          else

             if (debug) then
                open (11, file=trim(error_log), position="append")
                write (11,*) "error 052"
                close (11)
             end if

          end if

       end if

    end if

  end subroutine peel_surface

  subroutine peel_photon(thread_id, photon, x_in, y_in, z_in, stokes_peel_in, direction_in, cell_in, face_in, cell_error)

    integer,  intent(in)    :: cell_in(3), face_in(2), thread_id
    integer(16), intent(in) :: photon
    integer                 :: cell(3), cell_out(3), face(2), next_face(2)
    integer                 :: i, angle_upper, angle_lower, ix, iy
    real(dp), intent(in)    :: x_in, y_in, z_in, stokes_peel_in(4), direction_in(3)
    real(dp)                :: x, y, z, face_distance, stokes_out(4), x_im, y_im, phi_old, x_new, y_new
    real(dp)                :: mu_scatter, phi_scatter, scatter_dummy(16), scatter(16), cos_phi, phi_new, c2p
    real(dp)                :: x0(16), x1(16), y0, y1, y_inter, tau_cell, tau_total, photon_weight, det(3), acos_mu_scatter
    logical, intent(out)    :: cell_error
    logical                 :: grid_exit

    x = x_in
    y = y_in
    z = z_in

    det(1) = det_dir(1)
    det(2) = det_dir(2)
    det(3) = det_dir(3)

    face = face_in
    cell = cell_in

    scatter = 0._dp

    ! Optical depth

    tau_total = 0._dp
    photon_weight = 0._dp
    cell_error = .false.

    do

       call cell_face(photon, x, y, z, det, face, next_face, face_distance, grid_exit, cell, cell_out, cell_error)

       tau_cell = face_distance*cell_opacity(cell(1), cell(2), cell(3), wl_count)
       tau_total = tau_total + tau_cell

       x = x + face_distance*det(1)
       y = y + face_distance*det(2)
       z = z + face_distance*det(3)

       if (grid_exit.or.cell_error.or.(next_face(1).eq.1.and.next_face(2).eq.cell_depth)) exit

       face = next_face
       cell = cell_out

    end do

    if (.not.cell_error) then

       if (grid_exit.and.tau_total.lt.50._dp) then

          photon_weight = exp(-tau_total)

          mu_scatter = direction_in(1)*det(1) + direction_in(2)*det(2) + direction_in(3)*det(3)

          if (mu_scatter.ge.1._dp) then
             mu_scatter = 1._dp - 1.e-10_dp
          else if (mu_scatter.le.-1._dp) then
             mu_scatter = -1._dp + 1.e-10_dp
          end if

          ! Scattering matrix

          acos_mu_scatter = acos(mu_scatter)

          if (mod(acos_mu_scatter*180._dp/pi,1._dp).gt.0.5_dp) then

             angle_upper = int(acos_mu_scatter*180._dp/pi) + 2
             angle_lower = int(acos_mu_scatter*180._dp/pi) + 1

          else

             angle_upper = int(acos_mu_scatter*180._dp/pi) + 1
             angle_lower = int(acos_mu_scatter*180._dp/pi)

          end if
          
          if (angle_upper.eq.1) then

             scatter = cell_scatter_matrix(cell_in(1),cell_in(2),cell_in(3),wl_count,:,1)
             
          else if (angle_lower.eq.180) then

             scatter = cell_scatter_matrix(cell_in(1),cell_in(2),cell_in(3),wl_count,:,180)
             
          else

             do i=1,16

                x0(i) = cell_scatter_matrix(cell_in(1),cell_in(2),cell_in(3),wl_count,i,angle_lower)
                x1(i) = cell_scatter_matrix(cell_in(1),cell_in(2),cell_in(3),wl_count,i,angle_upper)
                y0 = dble(angle_lower) - 0.5_dp
                y1 = dble(angle_upper) - 0.5_dp
                y_inter = acos_mu_scatter*180._dp/pi
                scatter_dummy(i) = (x1(i)-x0(i)) * (y_inter-y0) / (y1-y0) + x0(i)
             
             end do

             scatter = scatter_dummy

          end if

          ! Check phi difference between old and new direction

          phi_old = atan2(direction_in(2),direction_in(1))

          if (phi_old.lt.0._dp) phi_old = phi_old + 2._dp*pi
          if (phi_old.gt.2._dp*pi) phi_old = phi_old - 2._dp*pi

          phi_new = det_dir(5)
          
          if (phi_new.lt.0._dp) phi_new = phi_new + 2._dp*pi
          if (phi_new.gt.2._dp*pi) phi_new = phi_new - 2._dp*pi

          ! Azimuthal scattering angle

          if (abs(direction_in(3)).lt.1._dp) then

             cos_phi = (det(3) - direction_in(3)*mu_scatter) / &
                  ( sqrt(1._dp-mu_scatter*mu_scatter)*sqrt(1._dp-direction_in(3)*direction_in(3)) )

             if (cos_phi.ge.1._dp) cos_phi = 1._dp - 1.e-6_dp
             if (cos_phi.le.-1._dp) cos_phi = -1._dp + 1.e-6_dp
             
             phi_scatter = acos(cos_phi)

             if (phi_old-phi_new.ge.0._dp.and.phi_old-phi_new.lt.pi) then

                phi_scatter = 2._dp*pi - phi_scatter

             end if

             if (2._dp*pi+phi_old-phi_new.ge.0._dp.and.2._dp*pi+phi_old-phi_new.lt.pi) then

                phi_scatter = 2._dp*pi - phi_scatter

             end if

             if (phi_scatter.lt.0._dp) phi_scatter = phi_scatter + 2._dp*pi

             c2p = 2._dp*cos_phi*cos_phi-1._dp

             call stokes_rotation(mu_scatter, phi_scatter, c2p, stokes_peel_in, scatter, direction_in, &
                     det, stokes_out, .true.)

          else

             if (debug) then
                open (11, file=trim(error_log), position="append")
                write (11,*) "error 045"
                close (11)
             end if

          end if

          if (.not.cell_error) then

             x_im = y_in*cos_det_phi - x_in*sin_det_phi
             y_im = z_in*sin_det_theta - y_in*cos_det_theta*sin_det_phi - x_in*cos_det_theta*cos_det_phi

             if (det_angle.gt.0._dp) then

                x_new = x_im*cos_det_angle - y_im*sin_det_angle
                y_new = x_im*sin_det_angle + y_im*cos_det_angle
                
                x_im = x_new
                y_im = y_new
          
             end if
             
             ix = int(npix*(x_im+x_max)/(2._dp*x_max))+1
             iy = int(npix*(y_im+y_max)/(2._dp*y_max))+1

             if (photon_weight*stokes_out(1).gt.0._dp.and.photon_weight*stokes_out(1).lt.1.e100_dp) then

                detector_thread(ix,iy,1,1,thread_id+1) = detector_thread(ix,iy,1,1,thread_id+1) + photon_weight*stokes_out(1)
                detector_thread(ix,iy,2,1,thread_id+1) = detector_thread(ix,iy,2,1,thread_id+1) - photon_weight*stokes_out(2)
                detector_thread(ix,iy,3,1,thread_id+1) = detector_thread(ix,iy,3,1,thread_id+1) + photon_weight*stokes_out(3)
                detector_thread(ix,iy,4,1,thread_id+1) = detector_thread(ix,iy,4,1,thread_id+1) + photon_weight*stokes_out(4)

                detector_thread(ix,iy,1,2,thread_id+1) = &
                     detector_thread(ix,iy,1,2,thread_id+1) + (photon_weight*stokes_out(1))**2
                detector_thread(ix,iy,2,2,thread_id+1) = &
                     detector_thread(ix,iy,2,2,thread_id+1) + (photon_weight*stokes_out(2))**2
                detector_thread(ix,iy,3,2,thread_id+1) = &
                     detector_thread(ix,iy,3,2,thread_id+1) + (photon_weight*stokes_out(3))**2
                detector_thread(ix,iy,4,2,thread_id+1) = &
                     detector_thread(ix,iy,4,2,thread_id+1) + (photon_weight*stokes_out(4))**2

                detector_thread(ix,iy,1,3,thread_id+1) = detector_thread(ix,iy,1,3,thread_id+1) + 1._dp
                detector_thread(ix,iy,2,3,thread_id+1) = detector_thread(ix,iy,2,3,thread_id+1) + 1._dp
                detector_thread(ix,iy,3,3,thread_id+1) = detector_thread(ix,iy,3,3,thread_id+1) + 1._dp
                detector_thread(ix,iy,4,3,thread_id+1) = detector_thread(ix,iy,4,3,thread_id+1) + 1._dp
             
             else

                if (debug) then
                   open (11, file=trim(error_log), position="append")
                   write (11,*) "error 053"
                   write (11,*) photon_weight
                   write (11,*) stokes_out(1)
                   write (11,*) mu_scatter
                   write (11,*) phi_scatter
                   write (11,*) c2p
                   write (11,*) stokes_peel_in
                   write (11,*) scatter
                   write (11,*) direction_in
                   write (11,*) det
                   write (11,*) stokes_out
                   close (11)
                end if

             end if

          end if

       end if

    end if

  end subroutine peel_photon

  subroutine add_flow_global(thread_id, x, y, z, direction, energy, face_distance, cell)

    integer,  intent(in) :: thread_id, cell(3)
    real(dp), intent(in) :: x, y, z, face_distance, direction(3), energy
    real(dp)             :: rho, r_dir, theta_dir, phi_dir, cos_theta, sin_theta, cos_phi, sin_phi

    rho = sqrt(x*x+y*y)
    
    cos_theta = z/sqrt(x*x+y*y+z*z)
    sin_theta = sqrt(1._dp-cos_theta*cos_theta)

    cos_phi = x/rho
    sin_phi = y/rho

    r_dir = sin_theta*cos_phi*direction(1) + sin_theta*sin_phi*direction(2) + cos_theta*direction(3)
    theta_dir = cos_theta*cos_phi*direction(1) + cos_theta*sin_phi*direction(2) - sin_theta*direction(3)
    phi_dir = -sin_phi*direction(1) + cos_phi*direction(2)

    cell_flow_global(thread_id+1,1,cell(1),cell(2),cell(3)) = &
         cell_flow_global(thread_id+1,1,cell(1),cell(2),cell(3)) + r_dir*face_distance*energy

    cell_flow_global(thread_id+1,2,cell(1),cell(2),cell(3)) = &
         cell_flow_global(thread_id+1,2,cell(1),cell(2),cell(3)) + theta_dir*face_distance*energy

    cell_flow_global(thread_id+1,3,cell(1),cell(2),cell(3)) = &
         cell_flow_global(thread_id+1,3,cell(1),cell(2),cell(3)) + phi_dir*face_distance*energy

  end subroutine add_flow_global

  subroutine add_flow(thread_id, direction, energy, cell)

    integer,  intent(in) :: thread_id, direction, cell(3)
    real(dp), intent(in) :: energy

    if (direction.eq.1) then

       ! Energy flowing upward
       cell_flow(thread_id+1,1,cell(1),cell(2),cell(3)) = cell_flow(thread_id+1,1,cell(1),cell(2),cell(3)) + energy
       
    else if (direction.eq.2) then

       ! Energy flowing downward
       cell_flow(thread_id+1,2,cell(1),cell(2),cell(3)) = cell_flow(thread_id+1,2,cell(1),cell(2),cell(3)) + energy

    else if (direction.eq.3) then

       ! Energy flowing towards the south
       cell_flow(thread_id+1,3,cell(1),cell(2),cell(3)) = cell_flow(thread_id+1,3,cell(1),cell(2),cell(3)) + energy

    else if (direction.eq.4) then

       ! Energy flowing towards the north
       cell_flow(thread_id+1,4,cell(1),cell(2),cell(3)) = cell_flow(thread_id+1,4,cell(1),cell(2),cell(3)) + energy

    end if

  end subroutine add_flow

  subroutine hunt(xx, x, jlo)

    !> Table lookup subroutine `hunt` from Numerical Recipes.
    !!
    !! Given an array `xx(1:N)`, and given a value `x`, returns a value `jlo` such 
    !! that `x` is between `xx(jlo)` and `xx(jlo+1)`. `xx` must be monotonic, 
    !! either increasing or decreasing. `jlo=0` or `jlo=N` is returned to indicate 
    !! that `x` is out of range. `jlo` on input is taken as the initial guess for 
    !! `jlo` on output.
    !!
    !! @author W.H. Press, S.A. Teukolsky, W.T. Vetterling, B.P. Flannery: 
    !! Numerical Recipes in Fortran 90
    
    integer, intent(inout) :: jlo
    real(dp), intent(in)   :: x
    real(dp), intent(in)   :: xx(:)
    integer                :: n,inc,jhi,jm
    logical                :: ascnd

    n=size(xx)
    ascnd = (xx(n) >= xx(1))

    if (jlo <= 0 .or. jlo > n) then
       jlo=0
       jhi=n+1
    else
       inc=1
       if (x >= xx(jlo) .eqv. ascnd) then
          do
             jhi=jlo+inc
             if (jhi > n) then
                jhi=n+1
                exit
             else
                if (x < xx(jhi) .eqv. ascnd) exit
                jlo=jhi
                inc=inc+inc
             end if
          end do
       else
          jhi=jlo
          do
             jlo=jhi-inc
             if (jlo < 1) then
                jlo=0
                exit
             else
                if (x >= xx(jlo) .eqv. ascnd) exit
                jhi=jlo
                inc=inc+inc
             end if
          end do
       end if
    end if

    do
       if (jhi-jlo <= 1) then
          ! if (x == xx(n)) jlo=n-1
          ! if (x == xx(1)) jlo=1
          if (abs(x-xx(n)).lt.1.e-6_dp) jlo=n-1
          if (abs(x-xx(1)).lt.1.e-6_dp) jlo=1
          exit
       else
          jm=(jhi+jlo)/2
          if (x >= xx(jm) .eqv. ascnd) then
             jlo=jm
          else
             jhi=jm
          end if
       end if
    end do

  end subroutine hunt

  subroutine random_direction(thread_id, direction)

    integer, intent(in)   :: thread_id
    real(dp), intent(out) :: direction(3)
    real(dp)              :: r, xi, x, y, z

    call random(thread_id, xi)
    x = 2._dp*xi-1._dp
    call random(thread_id, xi)
    y = 2._dp*xi-1._dp
    call random(thread_id, xi)
    z = 2._dp*xi-1._dp

    r = sqrt(x**2+y**2+z**2)
    
    direction(1) = x/r
    direction(2) = y/r
    direction(3) = z/r

  end subroutine random_direction

  subroutine min_distance(x, y, z, cell, r_min)

    integer, intent(in)   :: cell(3)
    real(dp), intent(in)  :: x, y, z
    real(dp), intent(out) :: r_min
    real(dp)              :: r, rho, phi, r1, r2, t1, t2, p1, p2, cos_theta, sin_theta, cos_phi, sin_phi, sin_diff

    r_min = 1e100_dp

    rho = sqrt(x*x+y*y)
    r = sqrt(x*x+y*y+z*z)
    cos_theta = z/r
    
    r1 = abs(r-rfront(cell(1)))
    r2 = abs(rfront(cell(1)+1)-r)
    if (r1.lt.r_min) r_min = r1
    if (r2.lt.r_min) r_min = r2

    if (ntheta.gt.1) then

       ! Trigonometry product and sum formulas

       sin_theta = sqrt(1._dp-cos_theta*cos_theta)
       
       if (thetaplane(cell(2)).eq.1) then
          ! sin_diff = sin(theta-thetafront(cell(2)))
          sin_diff = sin_theta*theta_grid_cos(cell(2)) - cos_theta*theta_grid_sin(cell(2))
          t1 = abs(r*sin_diff)
       else if (thetaplane(cell(2)).eq.2) then
          t1 = abs(z)
       end if

       if (thetaplane(cell(2)+1).eq.1) then
          ! sin_diff = sin(thetafront(cell(2)+1)-theta)
          sin_diff = cos_theta*theta_grid_sin(cell(2)+1) - sin_theta*theta_grid_cos(cell(2)+1)
          t2 = abs(r*sin_diff)
       else if (thetaplane(cell(2)+1).eq.2) then
          t2 = abs(z)
       end if

       if (t1.lt.r_min) r_min = t1
       if (t2.lt.r_min) r_min = t2

    end if

    if (nphi.gt.1.) then

       sin_phi = y/rho
       cos_phi = x/rho

       sin_diff = sin_phi*phi_grid_cos(cell(3)) - cos_phi*phi_grid_sin(cell(3))
       p1 = abs(rho*sin_diff)

       if (cell(3)+1.lt.nphi) then
          sin_diff = cos_phi*phi_grid_sin(cell(3)+1) - sin_phi*phi_grid_cos(cell(3)+1)
       else if (cell(3)+1.eq.nphi) then
          sin_diff = cos_phi*phi_grid_sin(0) - sin_phi*phi_grid_cos(0)
       else
          if (debug) then
             open (11, file=trim(error_log), position="append")
             write (11,*) "error 067"
             close (11)
          end if
       end if

       p2 = abs(rho*sin_diff)

       if (p1.lt.r_min) then
          
          phi = atan2(y,x)

          if (phi.lt.0._dp) phi = phi+2._dp*pi
          
          if (phi-phifront(cell(3)).gt.pi/2._dp.or.phi-phifront(cell(3)).lt.0._dp) then
             if (debug) then
                open (11, file=trim(error_log), position="append")
                write (11,*) "error 068"
                write (11,*) phifront(cell(3)), phi
                close (11)
             end if
          else
             r_min = p1
          end if
          
       end if
       
       if (p2.lt.r_min) then

          phi = atan2(y,x)
          if (phi.lt.0._dp) phi = phi+2._dp*pi

          if (cell(3)+1.lt.nphi) then
             if (phifront(cell(3)+1)-phi.gt.pi/2._dp.or.phifront(cell(3)+1)-phi.lt.0._dp) then
                if (debug) then
                   open (11, file=trim(error_log), position="append")
                   write (11,*) "error 069"
                   write (11,*) phifront(cell(3)+1), phi
                   close (11)
                end if
             else
                r_min = p2
             end if
          elseif (cell(3)+1.eq.nphi) then
             if (2._dp*pi-phi.gt.pi/2._dp.or.2._dp*pi-phi.lt.0._dp) then
                if (debug) then
                   open (11, file=trim(error_log), position="append")
                   write (11,*) "error 070"
                   write (11,*) 2._dp*pi, phi
                   close (11)
                end if
             else
                r_min = p2
             end if
          end if

       end if
       
    end if

  end subroutine min_distance

end program artes
