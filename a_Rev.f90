module mod_timeseries
  !
  ! tools for handling time series data
  !
  implicit none

  private
  public :: timeseries, timeseries_read, timeseries_interp

  type :: timeseries
    real, allocatable :: time(:)
    real, allocatable :: vars(:,:)
    integer :: nt, nv, i_curr
  end type timeseries

contains

  subroutine count_elements(input, nel)

    character(len=*), intent(in) :: input
    integer, intent(out) :: nel

    character(len=512) :: temp
    integer :: i, len_input

    len_input = len_trim(input)
    temp = adjustl(input)

    nel = 0
    do i = 1, len_input
      if (temp(i:i) == ' ') then
        if (i > 1 .and. temp(i-1:i-1) /= ' ') nel = nel + 1
      endif
    enddo

    if (len_trim(temp) > 0 .and. temp(len_input:len_input) /= ' ') &
      nel = nel + 1

  end subroutine count_elements


  subroutine timeseries_read(fname, data)

    character(len=*), intent(in) :: fname
    type(timeseries), intent(out) :: data

    integer :: i, nt, nv, ios
    real, allocatable :: tmp(:)
    character(len=512) :: line

    ! count numbers of variables
    open(unit=10, file=trim(fname), status='old', action='read')
    read(10, '(a)', iostat=ios) ! skip header
    read(10, '(a)', iostat=ios) line
    if (ios /= 0) then
      print *, 'Error: cannot read ', fname
      stop
    endif

    call count_elements(line, nv)
    nv = nv - 1  ! because of the left column is time
    rewind(10)

    ! count numbers of steps
    nt = 0
    do
      read(10, *, iostat=ios)
      if (ios /= 0) exit
      nt = nt + 1
    enddo
    nt = nt - 1 ! because of one line header
    rewind(10)

    ! allocate arrays
    allocate (data%time(nt), data%vars(nt,nv))

    ! read data and store in timeseries
    data%i_curr = 1
    data%nt = nt
    data%nv = nv

    allocate (tmp(nv+1))
    read(10,*) ! skip header
    do i=1, nt
      read(10,*) tmp
      data%time(i) = tmp(1)
      data%vars(i,:) = tmp(2:)
    enddo
    deallocate (tmp)

    close(10)

  end subroutine timeseries_read


  subroutine timeseries_interp(data, time_interp, data_interp)

    implicit none
    type(timeseries), intent(inout) :: data
    real, intent(in) :: time_interp
    real, intent(out) :: data_interp(:)

    real :: w_
    integer :: i, j

    if (time_interp < data%time(1) .or. time_interp > data%time(data%nt)) then
      print *, 'error: current time is out of range of time series data'
      stop
    endif

    do i=data%i_curr, data%nt - 1
      if (data%time(i) <= time_interp .and. &
        time_interp <= data%time(i+1)) exit
    enddo

    data%i_curr = i

    w_ = (data%time(i+1) - time_interp)/(data%time(i+1) - data%time(i))

    do j=1, data%nv
      data_interp(j) = w_*data%vars(i,j) + (1.0-w_)*data%vars(i+1,j)
    enddo

  end subroutine timeseries_interp

end module mod_timeseries


module mod_parameter
  !
  ! constant parameters
  !
  implicit none

  ! math constant
  real, parameter :: pi=4*atan(1.0)

  ! physical constant
  real, parameter ::  &
    gravity=9.81,         & ! gravity, m/s2
    k_boltzmann=1.171e-6    ! Stefan Boltzmann's constant, kcal/m2/day/K4

  ! water
  real, parameter ::  &
    rho_water=998.20, &  ! reference water density(20 deg-C), kg/m3
    mu_water=1.00e-3, &  ! reference water viscosity(20 deg-C), Pa s
    c_water=1.0,      &  ! specific heat of water, kcal/kg/K
    l_vapor=595.         ! latent heat of water evaporation, kcal/kg

  ! suspended solid particle
  integer :: nps         ! number of particle size
  real :: rho_ss=1650.   ! supended solid particle density, kg/m3
  real, allocatable, dimension(:) :: &
    d_ss,             &  ! particle sizes, m/s
    w_ss                 ! settling velocity for each particle sizes, m/s

contains

  real function water_density(t)
    ! water density, kg/m3
    real, intent(in) :: t ! temperature, deg-C
    water_density = 0.0000400*t**3 - 0.00784*t**2 &
      + 0.0581700*t + 999.85497
  end function water_density

  real function water_viscosity(t)
    ! water viscosity, Pa s
    real, intent(in) :: t ! temperature, deg-C
    water_viscosity = 4.486e-7*t**2 - 4.597e-5*t + 0.00175574
  end function water_viscosity

  real function water_thermal_conductivity(t)
    ! water viscosity, Pa s
    real, intent(in) :: t ! temperature, deg-C
    water_thermal_conductivity = 4.1e-8*t**4 &
      - 5.2e-6*t**3 + 0.00021*t**2 - 0.00177*t + 0.58176
  end function water_thermal_conductivity


  real function turbid_density(t, c)
    ! turbid(mixture) water density, kg/m3
    real, intent(in) :: t ! temperature, deg-C
    real, intent(in) :: c ! suspended solids conc., g/m3
    real :: rho
    rho = water_density(t)
    turbid_density = rho + c/1000*(1. - rho/rho_ss)    
  end function turbid_density


  real function psat(t)
    ! saturation vapour pressure, mmHg
    real, intent(in) :: t    ! temperature, deg-C

    ! from KCC code
    !psat = 0.00045*t**3 + 0.00363*t**2 + 0.39636*t + 4.47110

    ! Murray (1966)
    ! F. W. Murray, On the computation of saturation vapor pressure, 
    ! Journal of Applied Meteorology, vol.6, pp.203-204. (1966)
    psat = 6.1078*exp(17.2693882*t/(t + 237.3)) ! mb
    psat = psat * 0.75006 ! mmHg
 
  end function psat

end module mod_parameter


module mod_domain
  !
  ! module for domain
  !
  use mod_timeseries, only : timeseries

  implicit none

  type :: domain

    !-- domin

    character(len=128) :: name      ! domain name (ex. 'd01', 'd02', ..)
    integer :: id=0                 ! domain id (= 1, 2, 3, ..)
    character(len=128) :: fname_nml ! namelist file of domain (ex. namelist.d01)

    !-- geometric variables

    integer :: nx, nz ! number of cells

    real, pointer, dimension(:) :: &
      x, z,      & ! axis of cell faces, m
      xc, zc,    & ! axis of cell centers, m
      dx, dz,    & ! cell size, m
      dxs, dzs,  & ! dx at u-position, dz at w-position, m
      z_bed        ! bed elevation, m

    real, pointer, dimension(:,:) :: &
      b,         & ! resevoir width, m
      au,        & ! cell face area at u-position, m
      aw,        & ! cell face area at w-position, m2
      aus,       & ! au at cell center, m2
      aws,       & ! aw at cell vertex, m2
      vol,       & ! cell volume, m3
      vols         ! cell volume at u-position, m3

    ! original variables dependent on water surface elevation
    real, pointer, dimension(:) :: &
      dz0,       & ! cell size, m
      dzs0         ! dx at u-position, dz at w-position, m
    !
    real, pointer, dimension(:,:) :: &
      au0,       & ! cell face area at u-position, m
      aw0,       & ! cell face area at w-position, m2
      aus0,      & ! au at cell center, m2
      aws0,      & ! aw at cell vertex, m2
      vol0,      & ! cell volume, m3
      vols0        ! cell volume at u-position, m3

    real, pointer, dimension(:) :: &
      vol_hgt,   & ! total volume with respect to height, m3
      area_hgt,  & ! surface area with respect to height, m2
      len_hgt,   & ! reservoir length with respect to height, m
      vol_lay      ! volume of each layers, m3

    integer, pointer, dimension(:) :: &
      k_bot,     & ! bottom k at cell center
      kc_bot       ! bottom k at cell face

    ! variables related to water surface
    integer :: k_srf   ! index k of the surface layer
    real :: z_srf      ! water surface height, m
    real :: w_srf      ! vertical velocity at water surface, m/s
    real :: dz_srf     ! depth of the surface layer, m
    real :: total_vol  ! total water volume, m3
    real :: q_total_vol! flowrate due to water surface level fluctuations, m3
    integer :: log_i_inlet

    !-- primary variables

    real, pointer, dimension(:,:) :: &
      u,      &  ! horizontal velocity, m/s
      w,      &  ! vertical velocity, m/s
      p,      &  ! pressure, Pa
      t          ! temperature, deg-C

    real, pointer, dimension(:,:,:) :: &
      c          ! concentraion of suspended solids, g/m3 (=mg/l)

    !-- model variavles

    real, pointer, dimension(:,:) :: &
      rho,    &  ! water density, kg/m3
      GRAV,   &  ! gravity from river gradient, 
      PGX,    &  ! pressure gradient or filterd pressure gradient in time filter
      PGXraw, &  ! new pressure gradient in time filter, 
      PGXold     ! old pressure gradient in time filter,

    real, pointer, dimension(:) :: &
      rho_avg    ! horizontary averaged water density at w-pos, kg/m3


    real, pointer, dimension(:,:) :: &
      c_sed      ! sedimentation rate of suspended solids, g/s

    real :: &
      dmx0,  &  ! coeff of dmx, 1/day
      dmz0,  &  ! coeff of dmz, 1/day
      dhx0,  &  ! coeff of dhx, 1/day
      dhz0,  &  ! coeff of dhz, 1/day
      dcx0,  &  ! coeff of dcx, 1/day
      dcz0,  &  ! coeff of dcz, 1/day
      ll,    &  ! strength param of the expornent
      mm,    &  ! strength param of the expornent
      nn,    &  ! strength param of the expornent
      dmix      ! vertical diffusivity of mixing zone, m2/s

    real, pointer, dimension(:,:) :: &
      dmx,  &  ! horizontal eddy viscosity, m2/s
      dmz,  &  ! vertical eddy viscosity, m2/s
      dhx,  &  ! horizontal eddy diffusivity of heat, m2/s
      dhz,  &  ! vertical eddy diffusivity of heat, m2/s
      dcx,  &  ! horizontal eddy diffusivity of conc., m2/s
      dcz      ! vertical eddy diffusivity of conc., m2/s

    ! k-epsilon model variables
    real, pointer, dimension(:,:) :: &
      tke,     & ! turbulence kinetic energy
      td_eps,  & ! turbulence dissipation rate
      nut        ! eddy viscosity from k-epsilon, m2/s
    logical :: use_kepsilon = .false.
    logical :: semi_implicit = .false.
    logical :: freeslip = .false.

    !-- boundary conditions

    ! column flowrate
    real, pointer :: q_col(:)

    ! notes: 
    !   internal Froude number is 0.25 for 2d-flow, 
    !   0.134 for axisymmetric flow, 0.324 for surface and bottom layer flow

    ! inlet flow
    integer :: id_up      ! domain id of upstream (=0: read from fname_in)
    real ::          &
      q_in,          &    ! inlet flowrate, m3/s
      t_in,          &    ! inlet temperature, deg-C
      rho_in,        &    ! inlet density, kg/m3
      z_in,          &    ! median level of inlet flow, m
      delta_in,      &    ! depth of inlet flow, m
      fr_in,         &    ! internal Froude number of inflow
      b_in,          &    ! inflow width, m
      z_in_low            ! lower limit of inflow elevation, m
    real, pointer :: &
      c_in(:)             ! inlet SS concentration, g/m3
    type(timeseries), pointer :: &
      ts_in               ! time series of inlet
    integer, pointer :: &
      i_inlet(:)          ! index i of inlet face
    real :: qtot_in       ! total flowrate at inlet

    ! outlet flow
    integer :: n_out      ! number of outlets
    logical, pointer :: &
      surf_outs(:)        ! outlet at water surface(=T) or in-water(=F)
    real, pointer :: &
      q_outs(:),     &    ! outlet flowrates
      t_outs(:),     &    ! outlet temperature, deg-C
      c_outs(:,:),   &    ! outlet SS concentration, kg/m3
      z_outs(:),     &    ! median level of outlet flow, m
      fr_outs(:),    &    ! internal Froude number of outflow
      phi_outs(:),   &    ! outflow width, m
      out_heights(:), &    ! USE ONLY "DDD" and "EDO", outflow gate height, m
      zKTSWs(:),       &    ! USE ONLY "DDD", KTOP limitation level, m
      zKBSWs(:)             ! USE ONLY "DDD", KBOT limitation level, m
    type(timeseries), pointer :: &
      ts_out              ! time series of outlets
    real, pointer :: &
      t_out, c_out(:)     ! averaged values
    real :: qtot_out      ! total flowrate at outlet

    ! tributary flow
    integer :: n_trb      ! number of tributaries
    real :: &
      z_trb_low=-999.     ! lower limit of tributary elevation, m
    real, pointer :: &
      q_trbs(:),     &    ! tributary flowrate, m3/s
      t_trbs(:),     &    ! tributary temperature, deg-C
      c_trbs(:,:),   &    ! tributary SS concentration, g/m3
      fr_trbs(:),    &    ! internal Froude number of tributary flow
      b_trbs(:),     &    ! triburary width, m
      theta_trbs(:), &    ! tributary injection angle, deg
      z_trbs(:)           ! median level of tributary flow, m
    type(timeseries), pointer :: &
      ts_trbs(:)          ! time series of tributaries
    integer, pointer :: &
      i_trbs(:)           ! index i of tributaries
    real, pointer :: &
      q_trb(:,:),    &    ! q in x-z space, m3/s
      u_trb(:,:),    &    ! u in x-z space, m/s
      t_trb(:,:),    &    ! t in x-z space, deg-C
      c_trb(:,:,:)        ! c in x-z space, g/m3
    real :: qtot_trb      ! total flowrate of tributaries

    ! cnfluence
    integer :: n_cnf      ! number of confluences
    integer, pointer :: &
      id_cnfs(:)          ! domain id of confluence connections
    real, pointer :: &
      q_cnfs(:),     &    ! confluence flowrate, m3/s
      t_cnfs(:),     &    ! confluence temperature, deg-C
      c_cnfs(:,:),   &    ! confluence SS concentration, g/m3
      theta_cnfs(:)       ! confluence injection angle, deg
    integer, pointer :: &
      i_cnfs(:)           ! index i of confluences
    real, pointer :: &
      q_cnf(:,:),    &    ! q in x-z space, m3/s
      u_cnf(:,:),    &    ! u in x-z space, m/s
      t_cnf(:,:),    &    ! t in x-z space, deg-C
      c_cnf(:,:,:)        ! c in x-z space, g/m3
    real :: qtot_cnf      ! total flowrate of confluences

    ! water pipe flow
    integer :: n_wtp      ! number of water pipes
    integer, pointer :: &
      id_wtps(:)          ! domain id of water pipe connections
    real, pointer :: &
      q_wtps(:),     &    ! water pipe flowrate
      t_wtps(:),     &    ! water pipe temperature
      c_wtps(:,:),   &    ! water pipe SS concentration, g/m3
      fr_wtps(:),    &    ! internal Froude number of water pipe flow
      phi_wtps(:),   &    ! aperture angle of water pipe, deg
      theta_wtps(:), &    ! water pipe injection angle
      z_wtps(:)           ! median level of pipe flow, m
    type(timeseries), pointer :: &
      ts_wtps(:)          ! time series of water pipes
    integer, pointer :: &
      i_wtps(:)           ! index i of water pipes
    real, pointer :: &
      q_wtp(:,:),    &    ! q in x-z space
      u_wtp(:,:),    &    ! u in x-z space
      t_wtp(:,:),    &    ! t in x-z space
      c_wtp(:,:,:)        ! c in x-z space
    real :: qtot_wtp      ! total flowrate of water pipes

!---
    ! point inflow
    integer :: n_pin      ! number of point inflow
    real, pointer :: &
      q_pins(:),     &    ! point inflow flowrate
      t_pins(:),     &    ! point inflow temperature
      c_pins(:,:),   &    ! point inflow SS concentration, g/m3
      fr_pins(:),    &    ! internal Froude number of point inflow
      phi_pins(:),   &    ! aperture angle of point inflow, deg
      theta_pins(:), &    ! point inflow injection angle
      z_pins(:)           ! median level of point inflow, m
    type(timeseries), pointer :: &
      ts_pins(:)          ! time series of point inflow
    integer, pointer :: &
      i_pins(:),     &    ! index i of point inflow
      k_pins(:)           ! index k of point inflow
    real, pointer :: &
      q_pin(:,:),    &    ! q in x-z space
      u_pin(:,:),    &    ! u in x-z space
      t_pin(:,:),    &    ! t in x-z space
      c_pin(:,:,:)        ! c in x-z space
    real :: qtot_pin      ! total flowrate of point inflows

    ! point outflow
    integer :: n_pout     ! number of point outflow
    real, pointer :: &
      q_pouts(:),    &    ! point outflow flowrate
      t_pouts(:),    &    ! point outflow temperature
      c_pouts(:,:),  &    ! point outflow SS concentration, g/m3
      fr_pouts(:),   &    ! internal Froude number of point outflow
      phi_pouts(:),  &    ! aperture angle of point outflow, deg
      z_pouts(:)          ! median level of point outflow, m
    type(timeseries), pointer :: &
      ts_pout             ! time series of point outflow
    integer, pointer :: &
      i_pouts(:),    &    ! index i of point outflow
      k_pouts(:)          ! index k of point outflow
    real, pointer :: &
      q_pout(:,:)         ! q in x-z space
    real :: qtot_pout     ! total flowrate of point outflows
!---

    ! fence
    integer :: n_fnc      ! number of fences
    logical, pointer :: &
      type_fncs(:)        ! fence type (=T: floating, =F: fixed at bed)
    real, pointer ::    &
      width_fncs(:)       ! fence width
    integer, pointer :: &
      i_fncs(:),        & ! i-indices of fences
      k_fncs(:,:)         ! k-indices of fences

    !-- probe

    integer :: n_prb      ! number of probing points
    integer, pointer :: &
      i_prb(:),         & ! index i of probes, m
      k_prb(:)            ! index k of probes, m
    real, pointer ::    &
      x_prb(:),         & ! x of probes, m
      z_prb(:)            ! z of probes, m, (negative value is water surface)

    logical, pointer :: &
      exceed_logged(:,:) ! exceed logging flag

#ifdef SCALAR
    real, pointer ::     &
      s_in,              &  ! inlet scalar, -
      s_out,             &  ! outlet scalar, -
      s(:,:),            &  ! scalar field, -
      s_outs(:),         &  ! outlet scalar, -
      s_trbs(:),         &  ! tributary scalar, -
      s_trb(:,:),        &  ! s in x-z space, -
      s_cnfs(:),         &  ! confluence scalar, -
      s_cnf(:,:),        &  ! s in x-z space, -
      s_wtps(:),         &  ! water pipe scalar, -
      s_wtp(:,:),        &  ! s in x-z space, -
      s_pins(:),         &  ! point inflow scalar, -
      s_pin(:,:),        &  ! s in x-z space, -
      s_pouts(:),        &  ! point inflow scalar, -
      s_pout(:,:)           ! s in x-z space, -
#endif

  end type domain

  ! all domains
  type(domain), pointer :: doms(:)


contains

  subroutine set_geometry(dom)
    !
    ! input  : fname_geo
    ! output : x, z, z_bed, b, 
    !          xc, zc, dx, dz, dxs, dzs, au, aw, aus, aws, vol, vols,
    !          vol_hgt, area_hgt, len_hgt, k_bot, kc_bot,
    !          dz0, dzs0, au0, aw0, aus0, aws0, vol0, vols0
    !
    implicit none

    type(domain), intent(inout) :: dom

    integer :: nx, nz
    real, pointer, dimension(:) :: x, z, xc, zc, dx, dz, dxs, dzs, z_bed
    real, pointer, dimension(:,:) :: b, au, aw, aus, aws, vol, vols
    real, pointer, dimension(:) :: dz0, dzs0
    real, pointer, dimension(:,:) :: au0, aw0, aus0, aws0, vol0, vols0
    real, pointer, dimension(:) :: vol_hgt, area_hgt, len_hgt, vol_lay
    integer, pointer, dimension(:) :: k_bot, kc_bot

    integer :: i, k
    real :: tmp

    character(len=128) :: fname_geo=''
    namelist /geometry/ fname_geo

    ! loading from namelist file
    open(10, file=trim(dom%fname_nml), status='old')
    read(10, nml=geometry)
    close(10)

    write(6, '(a)') trim(dom%name)
    write(6, nml=geometry)

    !-- loading geometric data --
    open(10, file=trim(fname_geo), status='old')

    read(10,*)
    read(10,*) nx, nz

    allocate(x(0:nx), z(0:nz), z_bed(0:nx), b(0:nx,0:nz))

    read(10,*)
    read(10,*)

    do i=0, nx
      read(10,*) x(i), z_bed(i)
    enddo

    read(10,*)
    read(10,*)

    do k=0, nz
      read(10,*) z(k)
    enddo

    read(10,*)
    read(10,*)

    do k=0, nz
      read(10,*) (b(i,k), i=0, nx)
    enddo

    close(10)

    do i=0, nx
      tmp = 0.0
      do k=nz, 0, -1
        if (z(k) >= z_bed(i)) tmp = b(i,k)
        if (z(k) < z_bed(i)) b(i,k) = tmp 
      enddo
    enddo
    !---------------------

    allocate (xc(nx), zc(nz), &
      dx(nx), dz(nz), dxs(nx-1), dzs(nz-1), &
      au(0:nx,nz), aw(nx,0:nz), aus(nx,nz), aws(nx-1,0:nz), &
      vol(nx,nz), vols(nx-1,nz), vol_hgt(0:nz), & 
      area_hgt(0:nz), len_hgt(0:nz), vol_lay(nz), &
      k_bot(0:nx), kc_bot(nx))

    allocate (dz0(nz), dzs0(nz-1), &
      au0(0:nx,nz), aw0(nx,0:nz), aus0(nx,nz), aws0(nx-1,0:nz), &
      vol0(nx,nz), vols0(nx-1,nz))

    ! cell center position, cell size 
    do i=1, nx
      xc(i) = (x(i) + x(i-1))/2
      dx(i) = x(i) - x(i-1)
    enddo

    do i=1, nx-1
      dxs(i) = (dx(i) + dx(i+1))/2
    enddo

    do k=1, nz
      zc(k) = (z(k) + z(k-1))/2
      dz(k) = z(k) - z(k-1)
    enddo

    do k=1, nz-1
      dzs(k) = (dz(k) + dz(k+1))/2
    enddo

    ! cell face area at u-formation
    do i=0, nx
      do k=1, nz
        au(i,k) = (b(i,k-1) + b(i,k))/2 * dz(k)
      enddo
    enddo

    ! cell face area at w-formation
    do i=1, nx
      do k=0, nz
        aw(i,k) = (b(i-1,k) + b(i,k))/2 * dx(i)
      enddo
    enddo

    ! cell face area at u-stag
    do i=1, nx
      do k=1, nz
        aus(i,k) = (au(i-1,k) + au(i,k))/2
      enddo
    enddo

    ! cell face area at w-stag
    do i=1, nx-1
      do k=0, nz
        aws(i,k) = (aw(i,k) + aw(i+1,k))/2
      enddo
    enddo

    ! cell volume
    do i=1, nx
      do k=1, nz
        vol(i,k) = (aw(i,k-1) + aw(i,k))/2 * dz(k)
      enddo
    enddo

    ! cell volume at u-formation
    do i=1, nx-1
      do k=1, nz
        vols(i,k) = (vol(i,k) + vol(i+1,k))/2
      enddo
    enddo

    ! bottom k at cell face
    do i=0, nx
      k_bot(i) = 1
      do k=2, nz
        if (zc(k-1) <= z_bed(i) .and. z_bed(i) <= zc(k)) &
          k_bot(i) = k
      enddo
    enddo

    ! bottom k at cell center
    do i=1, nx
      kc_bot(i) = 1
      do k=1, nz
        kc_bot(i) = min(k_bot(i-1), k_bot(i))
      enddo
      if (i==nx) k_bot(nx) = kc_bot(nx) !!!
    enddo

    ! layer volume
    vol_lay(:) = 0.
    do i=1, nx
      do k=kc_bot(i), nz
        vol_lay(k) = vol_lay(k) + vol(i,k)
      enddo
    enddo

    ! total water volume with respect to height
    vol_hgt(0) = 0.
    do k=1, nz
      vol_hgt(k) = vol_hgt(k-1) + vol_lay(k)
    enddo

    ! surface area with respect to height
    area_hgt = 0.
    do i=1, nx
      do k=kc_bot(i)-1, nz
        area_hgt(k) = area_hgt(k) + aw(i,k)
      enddo
    enddo

    ! total length with respect to height
    len_hgt(:) = x(nx)

    do k=0, minval(k_bot)
      len_hgt(k) = 0.
    enddo

    do k=minval(k_bot), nz
      do i=0, nx-1
        if (z_bed(i) > z(k) .and. z(k) >= z_bed(i+1)) then
          tmp = (z_bed(i+1) - z(k))/(z_bed(i+1) - z_bed(i))
          len_hgt(k) = x(nx) - (tmp*x(i) + (1.-tmp)*x(i+1))
          exit
        endif
      enddo
    enddo

    ! save the original geometry
    dz0 = dz;  dzs0 = dzs
    au0 = au;  aw0 = aw;  aus0 = aus;  aws0 = aws
    vol0 = vol;  vols0 = vols

    !--- pointer ----------------------------------------
    dom%nx=nx; dom%nz=nz
    dom%x=>x; dom%z=>z; dom%z_bed=>z_bed; dom%b=>b
    dom%xc=>xc; dom%zc=>zc; dom%dx=>dx; dom%dz=>dz; dom%dxs=>dxs; dom%dzs=>dzs
    dom%au=>au; dom%aw=>aw; dom%aus=>aus; dom%aws=>aws
    dom%vol=>vol; dom%vols=>vols; dom%vol_hgt=>vol_hgt
    dom%area_hgt=>area_hgt; dom%len_hgt=>len_hgt; dom%vol_lay=>vol_lay
    dom%k_bot=>k_bot; dom%kc_bot=>kc_bot

    dom%dz0=>dz0; dom%dzs0=>dzs0
    dom%au0=>au0; dom%aw0=>aw0; dom%aus0=>aus0; dom%aws0=>aws0
    dom%vol0=>vol0; dom%vols0=>vols0
    !----------------------------------------------------

  end subroutine set_geometry


  subroutine allocate_variables(dom)

    use mod_parameter, only : nps

    implicit none

    type(domain), intent(inout) :: dom

    integer :: nx, nz
  
    ! allocate variables
    nx=dom%nx; nz=dom%nz

    ! primary variables
    allocate ( &
      dom%u(0:nx,nz), dom%w(nx,0:nz), dom%t(0:nx+1,0:nz+1), &
      dom%c(0:nx+1,0:nz+1,nps))

    ! model variables
    allocate ( &
      dom%p(0:nx+1,0:nz+1), dom%rho(0:nx+1,0:nz+1), dom%rho_avg(0:nz) )
 
    allocate ( &
      dom%dmx(nx,nz), dom%dmz(nx,0:nz), &
      dom%dhx(nx,nz), dom%dhz(nx,0:nz), &
      dom%dcx(nx,nz), dom%dcz(nx,0:nz))

    allocate (dom%tke(0:nx+1,0:nz+1), dom%td_eps(0:nx+1,0:nz+1))
    allocate (dom%nut(nx,nz))
    allocate (dom%PGX(nx,nz), dom%PGXold(nx,nz), dom%PGXraw(nx,nz))

    allocate (dom%exceed_logged(nx,0:nz))
    dom%exceed_logged = .false.

    allocate (dom%c_sed(0:nx+1,nps))
    dom%c_sed=0.0

#ifdef SCALAR
    allocate (dom%s(0:nx+1,0:nz+1))  ! scalar field
#endif

  end subroutine allocate_variables


  subroutine deallocate_variables(dom)

    implicit none

    type(domain), intent(inout) :: dom

    ! geometry variables
    nullify ( &
      dom%x, dom%z, dom%xc, dom%zc, dom%dx, dom%dz, dom%dxs, dom%dzs, dom%z_bed, &
      dom%b, dom%au, dom%aw, dom%aus, dom%aws, dom%vol, dom%vols, &
      dom%dz0, dom%dzs0, &
      dom%au0, dom%aw0, dom%aus0, dom%aws0, dom%vol0, dom%vols0, &
      dom%vol_hgt, dom%area_hgt, dom%len_hgt, &
      dom%k_bot, dom%kc_bot)
    
    ! primary variables
    nullify (dom%u, dom%w, dom%t, dom%c)

    ! model variables
    nullify (dom%p, dom%rho, dom%rho_avg, dom%PGX, dom%PGXraw, dom%PGXold)
    nullify (dom%dmx, dom%dmz, dom%dhx, dom%dhz, dom%dcx, dom%dcz)
    nullify (dom%tke, dom%td_eps, dom%nut)
    nullify (dom%c_sed)
    nullify (dom%exceed_logged, dom%c_sed)

#ifdef SCALAR
    nullify (dom%s)
#endif

  end subroutine deallocate_variables


  subroutine update_surface_layer(dom)
    !
    ! input: z_srf
    ! output: k_srf, dz_srf, dz, dzs, au, aus, aw, aws, vol, vols
    !
    implicit none

    type(domain), intent(inout) :: dom

    integer :: nx, nz, k_srf
    real :: z_srf, dz_srf
    real, pointer :: z(:), zc(:), dx(:), dz(:), b(:,:)

    integer :: i, k
    real, allocatable, dimension(:) :: b_srf, au_srf, aus_srf, &
      aw_srf, aws_srf, vol_srf, vols_srf

    !--- pointer ----------------------------------------
    nx=dom%nx; nz=dom%nz; z_srf=dom%z_srf
    z=>dom%z; zc=>dom%zc; dx=>dom%dx; dz=>dom%dz; b=>dom%b
    !----------------------------------------------------

    allocate (b_srf(0:nx), au_srf(0:nx), aus_srf(nx), &
      aw_srf(nx), aws_srf(nx-1), vol_srf(nx), vols_srf(nx-1))

    ! k_srf, dz_srf
    if (z(nz) < z_srf) then
      print *, 'Error: z_srf exceeds the upper limit of the domain.'
      stop
    endif

    if (z_srf <= zc(1)) then
      print *, 'Error: z_srf exceeds the lower limit of the domain.'
      stop
    endif

    k_srf = 0
    do k=1, nz-1
      if (zc(k) <= z_srf .and. z_srf < zc(k+1)) k_srf = k
    enddo
    if (k_srf==0) k_srf = nz

    dz_srf = z_srf - z(k_srf-1)

    ! b_srf, au_srf
    do i=0, nx
      if (z_srf <= z(k_srf)) then
        b_srf(i) = b(i,k_srf-1) + &
          dz_srf/dz(k_srf) * (b(i,k_srf) - b(i,k_srf-1))
        au_srf(i) = (b(i,k_srf-1) + b_srf(i))/2 * dz_srf
      else
        b_srf(i) = b(i,k_srf) + &
          (z_srf - z(k_srf))/dz(k_srf+1) * (b(i,k_srf+1) - b(i,k_srf))
        au_srf(i) = (b(i,k_srf-1) + b(i,k_srf))/2 * dz(k_srf) + &
          (b(i,k_srf) + b_srf(i))/2 * (z_srf - z(k_srf))
      endif
    enddo

    do i=1, nx
      aus_srf(i) = (au_srf(i-1) + au_srf(i))/2
      aw_srf(i) = (b_srf(i-1) + b_srf(i))/2 * dx(i)
      vol_srf(i) = aus_srf(i) * dx(i)
    enddo

    do i=1, nx-1
      aws_srf(i) = (aw_srf(i) + aw_srf(i+1))/2
      vols_srf(i) = (vol_srf(i) + vol_srf(i+1))/2
    enddo

    ! set surface layer to geometry
    dom%dz = dom%dz0;  dom%dzs = dom%dzs0
    dom%au = dom%au0;  dom%aus = dom%aus0
    dom%aw = dom%aw0;  dom%aws = dom%aws0
    dom%vol = dom%vol0;  dom%vols = dom%vols0

    k = k_srf
    dom%k_srf = k_srf
    dom%dz_srf = dz_srf

    dom%dz(k) = dz_srf;  dom%dzs(k) = dz_srf
    dom%au(:,k) = au_srf(:);  dom%aus(:,k) = aus_srf(:)
    dom%aw(:,k) = aw_srf(:);  dom%aws(:,k) = aws_srf(:)
    dom%vol(:,k) = vol_srf(:);  dom%vols(:,k) = vols_srf(:)

    deallocate (b_srf, au_srf, aus_srf, aw_srf, vol_srf, vols_srf)

  end subroutine update_surface_layer


  subroutine get_total_volume(dom)
    !
    ! inputs:  z_srf, vol_hgt
    ! outputs: total_vol
    !
    implicit none

    type(domain), intent(inout) :: dom

    integer :: nz
    real :: z_srf, total_vol
    real, pointer :: z(:), vol_hgt(:)

    integer :: k

    !--- pointer ---
    nz=dom%nz; z=>dom%z
    z_srf=dom%z_srf; vol_hgt=>dom%vol_hgt
    !---------------

    total_vol = 0.0

    do k=1, nz
      if (z(k-1) <= z_srf .and. z_srf <= z(k)) then
        total_vol = vol_hgt(k-1) + &
          (vol_hgt(k) - vol_hgt(k-1))*(z_srf - z(k-1))/(z(k) - z(k-1))
      endif
    enddo

    if (total_vol > vol_hgt(nz)) then
      print *, 'Error: total_vol > vol_hgt(nz)'
      stop
    endif

    if (total_vol < 1.e-3) then
      print *, 'Error: total_vol=0'
      stop
    endif

    dom%total_vol = total_vol

  end subroutine get_total_volume


  subroutine get_surface_height(nz, z, vol_hgt, total_vol, z_srf)

    implicit none

    integer, intent(in) :: nz
    real, intent(in) :: z(0:nz), vol_hgt(0:nz), total_vol
    real, intent(out) :: z_srf

    integer :: k

    if (vol_hgt(nz) < total_vol) then
      print *, 'water volume has exceeded reservoir capacity' 
      stop
    endif

    z_srf = -999.

    do k=1, nz
      if (vol_hgt(k-1) < total_vol .and. total_vol <= vol_hgt(k)) then
        z_srf = z(k-1) + (z(k) - z(k-1))/(vol_hgt(k) - vol_hgt(k-1)) * &
          (total_vol - vol_hgt(k-1))
      endif
    enddo

    if (z_srf > z(nz)) then
      print *, 'Error: z_srf > z(nz)'
      stop
    endif

    if (z_srf < -900.) then
      print *, 'Error: z_srf is out of bounds'
      stop
    endif

  end subroutine get_surface_height


  subroutine set_initial(dom)
    !
    ! input  : fname_nml
    ! output : z_srf, u, w, t, c
    !
    use mod_parameter, only : nps, turbid_density

    implicit none

    type(domain), intent(inout) :: dom

    integer :: i, j, k, io_stat, k_, n_
    real :: w_
    real, allocatable :: z_(:), u_(:), t_(:), c_(:,:)
#ifdef SCALAR
    real, allocatable :: s_(:)
#endif

    real :: z_srf=40.      ! initial water surface level, m
    logical :: uniform=.false. ! uniform(T) or read from an file(F)
    real :: u_init=0.0     ! uniform velocity, m/s (used when uniform=T)
    real :: t_init=20.     ! uniform temp., m/s (used when uniform=T)
    real :: c_init(20)=0.0 ! uniform SS conc, g/m3 (used when uniform=T)
#ifdef SCALAR
    real :: s_init=0.      ! uniform scalar, - (used when uniform=T)
#endif
    character(len=128) :: fname_init=''

#ifdef SCALAR
    namelist /initial/ z_srf, uniform, u_init, t_init, c_init, s_init, fname_init
#else
    namelist /initial/ z_srf, uniform, u_init, t_init, c_init, fname_init
#endif

    ! loading from namelist file
    open(10, file=trim(dom%fname_nml), status='old')
    read(10, nml=initial)
    close(10)

    write(6, '(a)') trim(dom%name)
    write(6, nml=initial)

    ! set initial surface height, m
    dom%z_srf = z_srf

    ! set uniform distribution
    if (uniform) then
      dom%u = u_init
      dom%w = 0.0
      dom%t = t_init
      do j=1, nps
        dom%c(:,:,j) = c_init(j)
      enddo
#ifdef SCALAR
      dom%s = s_init
#endif
      return
    endif

    ! read from a file

    ! count number of data
    open(10, file=trim(fname_init), status='old')
    i = 0
    do
      read(10, '(a)', iostat=io_stat)
      if (io_stat /= 0) exit
      i = i + 1
    enddo
    n_ = i - 1  ! skip header line
    close(10)

    ! allocate arrays
#ifdef SCALAR
    allocate (z_(n_), u_(n_), t_(n_), c_(n_,nps), s_(n_))
#else
    allocate (z_(n_), u_(n_), t_(n_), c_(n_,nps))
#endif

    ! read data
    open(10, file=trim(fname_init), status='old')
    read(10,*)  ! skip header line
    do i=1, n_
#ifdef SCALAR
      read(10,*) z_(i), u_(i), t_(i), (c_(i,j), j=1, nps), s_(i)
#else
      read(10,*) z_(i), u_(i), t_(i), (c_(i,j), j=1, nps)
#endif
    enddo
    close(10)

    ! horizontal velocity, u
    do k=1, dom%nz

      do k_=1, n_-1
        if (z_(k_) <= dom%zc(k) .and. dom%zc(k) <= z_(k_+1)) then
          w_ = (z_(k_+1) - dom%zc(k))/(z_(k_+1) - z_(k_))
          dom%u(:,k) = w_*u_(k_) + (1.0-w_)*u_(k_+1)
          dom%t(:,k) = w_*t_(k_) + (1.0-w_)*t_(k_+1)
          do j=1, nps
            dom%c(:,k,j) = w_*c_(k_,j) + (1.0-w_)*c_(k_+1,j)
          enddo
#ifdef SCALAR
          dom%s(:,k) = w_*s_(k_) + (1.0-w_)*s_(k_+1)
#endif
        endif
      enddo

      dom%t(:,0) = dom%t(:,1)
      dom%t(:,dom%nz+1) = dom%t(:,dom%nz)
      do j=1, nps
        dom%c(:,0,j) = dom%c(:,1,j)
        dom%c(:,dom%nz+1,j) = dom%c(:,dom%nz,j)
      enddo
#ifdef SCALAR
      dom%s(:,0) = dom%s(:,1)
      dom%s(:,dom%nz+1) = dom%s(:,dom%nz)
#endif

    enddo

    ! vertical velocity, w
    dom%w = 0.0

    ! turbulence quantities (k-epsilon)
    dom%tke = 1.25e-7
    dom%td_eps = 1.0e-9
    dom%nut = 0.09 * dom%tke**2 / dom%td_eps

    deallocate (z_, u_, t_, c_)
#ifdef SCALAR
    deallocate (s_)
#endif

  end subroutine set_initial


  subroutine write_geo(dom)
    !
    ! input  : nx, nz, x, z, k_bot, kc_bot, i_inlet, z_bed, b, 
    !          vol_hgt, area_hgt, len_hgt
    ! output : 'out/d01/geo.dat'
    !
    implicit none

    type(domain), intent(in) :: dom

    open(10, file='out/' // trim(dom%name) // '/geo.dat' , status='unknown')

    write(10,'(2i5)') dom%nx, dom%nz
    write(10,'(1000e15.7)') dom%x
    write(10,'(1000e15.7)') dom%z
    write(10,'(1000i5)') dom%k_bot
    write(10,'(1000i5)') dom%kc_bot
    write(10,'(1000i5)') dom%i_inlet
    write(10,'(1000e15.7)') dom%z_bed
    write(10,'(1000000e15.7)') dom%b
    write(10,'(1000e15.7)') dom%vol_hgt
    write(10,'(1000e15.7)') dom%area_hgt
    write(10,'(1000e15.7)') dom%len_hgt

    close(10)

  end subroutine write_geo


  subroutine update_density(dom)
    !
    ! input  : t, c
    ! output : p, rho, rho_avg
    !
    use mod_parameter, only : gravity, turbid_density

    implicit none

    type(domain), intent(inout) :: dom

    integer :: nx, nz
    integer :: i, k, k_srf
    real, allocatable :: volay(:), rhoav(:)

    !--- pointer ---
    nx=dom%nx; nz=dom%nz
    !----------------

    if (dom%k_srf==0) then
      k_srf = dom%nz
    else
      k_srf = dom%k_srf
    endif

    do i=1, nx

      ! density, kg/m3
      do k=dom%kc_bot(i), k_srf
        dom%rho(i,k) = turbid_density(dom%t(i,k), sum(dom%c(i,k,:)))
      enddo

      ! ghost cell
      dom%rho(i,k_srf+1) = dom%rho(i,k_srf)

      ! hydrostatic pressure, Pa
      dom%p(i,k_srf) = dom%rho(i,k_srf) * gravity * (dom%z_srf - dom%zc(k_srf))
      do k=k_srf-1, dom%kc_bot(i), -1
        dom%p(i,k) = dom%p(i,k+1) &
          + (dom%rho(i,k) + dom%rho(i,k+1))/2 * gravity * dom%dzs(k)
      enddo

    enddo

    ! horizontaly averaged density

    ! rho level
    allocate (volay(nz), rhoav(nz))

    volay = 1e-6
    rhoav = 0.0

    do i=1, nx
      do k=dom%kc_bot(i), k_srf
        volay(k) = volay(k) + dom%vol(i,k)
        rhoav(k) = rhoav(k) + dom%vol(i,k)*dom%rho(i,k)
      enddo
    enddo

    do k=1, k_srf
      rhoav(k) = rhoav(k)/volay(k)
    enddo

    deallocate (volay)

    ! w level
    dom%rho_avg(0) = rhoav(1)
    dom%rho_avg(k_srf) = rhoav(k_srf)
    dom%rho_avg(min(k_srf+1,nz)) = rhoav(k_srf)

    do k=1, k_srf-1
      dom%rho_avg(k) = (rhoav(k) + rhoav(k+1))/2
    enddo

  end subroutine update_density


  subroutine set_turbulence(dom)

    implicit none

    type(domain), intent(inout) :: dom

    ! namelist
    real ::        &
      dmx0 = 5.0,  &  ! coeff of dmx, 1/day
      dmz0 = 1e-3, &  ! coeff of dmz, 1/day
      dhx0 = 5.0,  &  ! coeff of dhx, 1/day
      dhz0 = 1e-3, &  ! coeff of dhz, 1/day
      dcx0 = 5.0,  &  ! coeff of dcx, 1/day
      dcz0 = 1e-3, &  ! coeff of dcz, 1/day
      ll = 0.5,    &  ! strength param of the expornent
      mm = 0.5,    &  ! strength param of the expornent
      nn = 0.5,    &  ! strength param of the expornent
      dmix = 1e-3     ! vertical diffusivity of mixing zone, m2/s
    logical :: use_kepsilon = .false.   ! impose(T) or expose(F) k-epsilon 
    logical :: semi_implicit= .false.   ! impose(T) or expose(F) semi_implicit method (TRIDIAG) 
    logical :: freeslip = .false.       ! impose(T) or expose(F) freeslip in k-epsilon 

    namelist /turbulence/ dmx0, dmz0, dhx0, dhz0, dcx0, dcz0, &
      ll, mm, nn, dmix, use_kepsilon, semi_implicit, freeslip

    ! loading from namelist file
    open(10, file=trim(dom%fname_nml), status='old')
    read(10, nml=turbulence)
    close(10)

    write(6, '(a)') trim(dom%name)
    write(6, nml=turbulence)

    dom%dmx0=dmx0; dom%dmz0=dmz0
    dom%dhx0=dhx0; dom%dhz0=dhz0
    dom%dcx0=dcx0; dom%dcz0=dcz0
    dom%ll=ll; dom%mm=mm; dom%nn=nn 
    dom%dmix=dmix
    dom%use_kepsilon=use_kepsilon
    dom%semi_implicit=semi_implicit
    dom%freeslip=freeslip

  end subroutine set_turbulence


  subroutine update_turbulence(dom, dt_sec)
    !
    ! input  : u, t, rho
    ! output : dmx, dmz, dhx, dhz
    !
    use mod_parameter, only : gravity, rho_water, c_water, &
      water_viscosity, water_density, water_thermal_conductivity

    implicit none

    type(domain), intent(inout) :: dom
    real, intent(in) :: dt_sec

    integer :: i, k, nx, nz, k_srf
    integer, pointer, dimension(:) :: k_bot, kc_bot
    real, pointer, dimension(:) :: dx, dz, dxs, dzs, x, z_bed
    real, pointer, dimension(:,:) :: au, aw, aus, aws, vols
    real, pointer, dimension(:,:) :: rho, p, dmx, dmz, tke, td_eps
    real :: ri_, u1_, u2_, uz_, rhoz_, rho_, &
      t_, nu_, nut_, alpha_, alphat_, beta_, betat_
    ! k-e variables
    real :: Prdk, Prhe, Prhk, Unst, Unse
    real :: buoy_k, buoy_e, flux_up_tke, flux_down_tke, diff_tke, flux_up_e, flux_down_e, diff_e, dz_up, dz_down
    real :: u_starS, u_starBS, u_starB, kappa=0.41, z0b=1.0e-4, WIND10, CZ, GC2, FRIC, SIG1, SIG2 

    real, allocatable :: u(:,:), w(:,:), AT(:,:), CT(:,:), VT(:,:), DT(:,:), BTA1(:,:), GMA1(:,:)
    real, allocatable :: depth(:)

    !--- pointer, allocation ----------------------------
    nx=dom%nx; nz=dom%nz
    k_bot=>dom%k_bot; kc_bot=>dom%kc_bot; k_srf=dom%k_srf
    dx=>dom%dx; dxs=>dom%dxs; dz=>dom%dz; dzs=>dom%dzs; x=>dom%x; z_bed=>dom%z_bed
    rho=>dom%rho; p=>dom%p; dmx=>dom%dmx; dmz=>dom%dmz
    au=>dom%au; aw=>dom%aw; aus=>dom%aus; aws=>dom%aws; vols=>dom%vols
    tke=>dom%tke; td_eps=>dom%td_eps
    allocate (u(0:nx,nz), w(nx,0:nz), AT(nx,nz), CT(nx,nz), VT(nx,nz), DT(nx,nz), BTA1(nx,nz), GMA1(nx,nz))
    allocate(depth(nx))
    !----------------------------------------------------

    ! horizontal diffusivity
    do i=1, nx
      do k=1, nz

        ! temperature, deg-C
        t_ = dom%t(i,k)

        ! kinematic viscosity, m2/s
        nu_ = water_viscosity(t_)/water_density(t_)

        ! thermal diffusivity, m2/s
        alpha_ = water_thermal_conductivity(t_)/(rho_water*c_water*4186)

        ! concentration diffusivity, m2/s
        beta_ = nu_

        ! eddy viscosity, m2/s
        nut_ = (dom%dmx0/(3600*24))*dom%dx(i)**2

        ! eddy thermal diffusivity, m2/s
        alphat_ = (dom%dhx0/(3600*24))*dom%dx(i)**2

        ! eddy concentraion diffusivity, m2/s
        betat_ = nut_

        ! effective viscosity, m2/s
        dom%dmx(i,k) = nu_ + nut_

        ! effective thermal diffusivity, m2/s
        dom%dhx(i,k) = alpha_ + alphat_

        ! effective concentration diffusivity, m2/s
        dom%dcx(i,k) = beta_ + betat_

      enddo
    enddo

    if (dom%use_kepsilon) then
      ! --- from k-epsilon model ---
      do i = 1, nx
        ! WIND10 using Power Law
        ! WIND10 = wind*(10/2)**0.40
        WIND10 = 2.0D0

        ! surface shear velocity for wind10
        IF(WIND10 >= 15.0)THEN
          CZ = 0.0026D0
        ELSEIF(WIND10 >= 4.0)THEN
          CZ = 0.0005D0*SQRT(WIND10)
        ELSEIF(WIND10 >= 0.5)THEN
          CZ= 0.0044D0*WIND10**(-1.15D0)
        ELSE
          CZ= 0.01D0
        ENDIF

        u_starS = SQRT(1.25 * CZ * wind10**2 / rho(i,dom%k_srf))

        ! bottom shear velocity using simple model
        FRIC = 100.0
        GC2 = gravity / (FRIC**2)
        u_starBS = SQRT(GC2) * ABS(0.5* (u(i, dom%k_srf) + u(i-1, dom%k_srf)))

        ! set surface tke and td_eps
        depth(i) = dom%z_srf - dom%z_bed(i)
        tke(i, dom%k_srf) = ( 3.33*(u_starS**2 + u_starBS**2) )
        td_eps(i, dom%k_srf) = (u_starS**3 + u_starBS**3)*5.0 / depth(i)

        ! intermediate tke and td_eps
        do k = dom%kc_bot(i)+1, dom%k_srf-1
          ! shear production and buoyoncy
          buoy_k = MAX( dom%dmz(i,k)*gravity*(dom%rho(i,k-1) - dom%rho(i,k))/(rho_water*dom%dz(k)), 0.0 )
          ! select Prdk calc 
          Prdk = dom%dmz(i,k)*( 0.5*(dom%u(i,k) - dom%u(i,k-1) + dom%u(i-1,k)  - dom%u(i-1,k-1)) &
          / (dom%dz(k)*0.5 + dom%dz(k-1)*0.5) )**2
          !Prdk = dom%dmz(i,k)*( (dom%u(i,k) - dom%u(i,k-1))/ (dom%dz(k)*0.5 + dom%dz(k-1)*0.5) )**2
          Prhe = 10.0* GC2**1.25 *ABS( 0.5*(dom%u(i,k) + dom%u(i-1,k)))**4.0 &
          / (0.5*dom%au(i,k)/dom%dz(k))**2
          Prhk = GC2 / (0.5*dom%au(i,k)/dom%dz(k)) * ABS(0.5 * (dom%u(i,k) + dom%u(i-1,k)))**3

          ! update (Production - Dissipasion)
          dom%tke(i,k)    = max(dom%tke(i,k),    1.25d-7)
          dom%td_eps(i,k) = max(dom%td_eps(i,k), 1.0d-9)
          Unst = Prdk - dom%td_eps(i,k)
          Unse = 1.44 * dom%td_eps(i,k)/dom%tke(i,k)*Prdk - 1.92 * (dom%td_eps(i,k)**2/dom%tke(i,k))

          dom%tke(i,k) = dom%tke(i,k) + dt_sec*(Unst + Prhk - buoy_k)
          dom%td_eps(i,k) = dom%td_eps(i,k) + dt_sec*(Unse + Prhe)
          !print '(2i5, 8f14.10)',  i, k, buoy_k, Prdk, Prhe, Prhk, Unst, Unse, dom%tke(i, k), dom%td_eps(i,k)
        end do

        ! bottom tke and td_eps
        u_starB = SQRT(GC2) * ABS( 0.5 * (u(i,dom%kc_bot(i)) + u(i-1,dom%kc_bot(i))) )
        dom%tke(i,dom%kc_bot(i)) = 0.5 *( 3.33*u_starB**2 + tke(i,dom%kc_bot(i)))
        dom%td_eps(i,dom%kc_bot(i)) = 0.5 *( u_starB**3*5.0/dom%dzs(dom%kc_bot(i)) + td_eps(i,dom%kc_bot(i)) )

        dom%tke(i,dom%kc_bot(i)-1) = dom%tke(i,dom%kc_bot(i))
        dom%td_eps(i,dom%kc_bot(i)-1) = dom%td_eps(i,dom%kc_bot(i))

      enddo

      ! semi implicit using tridiagonal method estimate TKE and epsilon
      if (dom%semi_implicit) then
        SIG1 = 1.0      ! turbulence model coefficient for tke
        SIG2 = 1.3      ! turbulence model coefficient for epsilon

        ! update tke using tridiagonal method
        do i = 1, nx
          k = k_srf      
          AT(i,k) = 0.0
          CT(i,k) = 0.0
          VT(i,k) = 1.0
          DT(i,k) = dom%tke(i,k)

          do k = kc_bot(i), k_srf-1
            AT(i,k) = - dt_sec/dom%au(i,k)*dom%dx(i)/SIG1*dom%dmz(i,k-1)/((dom%dz(k) + dom%dz(k-1))*0.5)
            CT(i,k) = - dt_sec/dom%au(i,k)*dom%dx(i)/SIG1*dom%dmz(i,k)/((dom%dz(k) + dom%dz(k+1))*0.5)
            VT(i,k) = 1.0 - AT(i,k) - CT(i,k)
            DT(i,k) = dom%tke(i,k)
          end do 

          k = kc_bot(i)
          AT(i,k) = 0.0
          CT(i,k) = 0.0
          VT(i,k) = 1.0
          DT(i,k) = dom%tke(i,k) 

          ! TRIDIAG for tke
          k= k_srf
          BTA1(i, k) = VT(i,k)
          GMA1(i, k) = DT(i,k)

          do k= k_srf-1, kc_bot(i), -1
            BTA1(i, k) = VT(i,k) - AT(i,k)/BTA1(i,k+1)*CT(i,k+1)
            GMA1(i, k) = DT(i,k) - AT(i,k)/BTA1(i,k+1)*GMA1(i,k+1)
          end do


          do k= kc_bot(i), k_srf-1
            BTA1(i,k) = max(BTA1(i,k), 1.0e-10)
          end do

          dom%tke(i,kc_bot(i)) = GMA1(i, kc_bot(i))/BTA1(i,kc_bot(i))

          do k= kc_bot(i)+1, k_srf
            dom%tke(i,k) = ( GMA1(i,k) - CT(i,k)*dom%tke(i,k-1))/BTA1(i,k)

            !print '(2i5, 4f20.15)',  i, k, dom%nut(i,k), dom%td_eps(i, k), dom%tke(i, k), dom%dmz(i, k) 

          end do
        end do

        ! update td_eps using tridiagonal method
        do i = 1, nx
          k = k_srf      
          AT(i,k) = 0.0
          CT(i,k) = 0.0
          VT(i,k) = 1.0
          DT(i,k) = dom%td_eps(i,k)

          do k = kc_bot(i), k_srf-1
            AT(i,k) = - dt_sec/dom%au(i,k)*dom%dx(i)/SIG2*dom%dmz(i,k-1)/((dom%dz(k) + dom%dz(k-1))*0.5)
            CT(i,k) = - dt_sec/dom%au(i,k)*dom%dx(i)/SIG2*dom%dmz(i,k)/((dom%dz(k) + dom%dz(k+1))*0.5)
            VT(i,k) = 1.0 - AT(i,k) - CT(i,k)
            DT(i,k) = dom%td_eps(i,k)
          end do 

          k = kc_bot(i)
          AT(i,k) = 0.0
          CT(i,k) = 0.0
          VT(i,k) = 1.0
          DT(i,k) = dom%td_eps(i,k) 

          ! TRIDIAG for td_eps 
          k= k_srf
          BTA1(i, k) = VT(i,k)
          GMA1(i, k) = DT(i,k)

          do k= k_srf-1, kc_bot(i), -1
            BTA1(i, k) = VT(i,k) - AT(i,k)/BTA1(i,k+1)*CT(i,k+1)
            GMA1(i, k) = DT(i,k) - AT(i,k)/BTA1(i,k+1)*GMA1(i,k+1)
          end do

          do k= kc_bot(i), k_srf
            BTA1(i,k) = max(BTA1(i,k), 1.0e-10)
          end do

          dom%td_eps(i,kc_bot(i)) = GMA1(i, kc_bot(i))/BTA1(i,kc_bot(i))

          do k= kc_bot(i)+1, k_srf
            dom%td_eps(i,k) = ( GMA1(i,k) - CT(i,k)*dom%td_eps(i,k-1))/BTA1(i,k)

            !print '(2i5, 5f20.15)',  i, k, dom%au(i,k), dom%dx(i), dom%dz(k), BTA1(i,k), GMA1(i, k) 
            !print '(2i5, 4f20.15)',  i, k, dom%nut(i,k), dom%td_eps(i, k), dom%tke(i, k), dom%dmz(i, k) 

          end do
        end do

        deallocate(AT, CT, VT, DT, BTA1, GMA1)

      else    ! upwind method calc TKE and epsilon
        ! update vertical diffusion
        do i = 1, nx
          do k = dom%kc_bot(i), dom%k_srf-1
            flux_up_tke = dom%dmz(i,k) * (dom%tke(i,k+1) - dom%tke(i,k))/( 0.5*(dom%dz(k+1)+dom%dz(k)) )
            flux_down_tke = dom%dmz(i,k-1) * (dom%tke(i,k) - dom%tke(i,k-1))/( 0.5*(dom%dz(k)+dom%dz(k-1)) )
            flux_up_e = dom%dmz(i,k) * (dom%td_eps(i,k+1) - dom%td_eps(i,k))/( 0.5*(dom%dz(k+1)+dom%dz(k)) )
            flux_down_e = dom%dmz(i,k-1) * (dom%td_eps(i,k) - dom%td_eps(i,k-1))/( 0.5*(dom%dz(k)+dom%dz(k-1)) )
            diff_tke = (flux_up_tke - flux_down_tke) / dom%aw(i,k)
            diff_e = (flux_up_e - flux_down_e) / dom%aw(i,k)

            dom%tke(i,k) = dom%tke(i,k) + dt_sec*diff_tke
            dom%td_eps(i,k) = dom%td_eps(i,k) + dt_sec*diff_e

            ! avoid NaNs and non-positive values
            if (dom%tke(i,k) /= dom%tke(i,k) .or. dom%tke(i,k) <= 0.0) then
              dom%tke(i,k) = 1.25e-7
            end if
            if (dom%td_eps(i,k) /= dom%td_eps(i,k) .or. dom%td_eps(i,k) <= 0.0) then
              dom%td_eps(i,k) = 1.0e-9
            end if

          !print '(2i5, 3f30.10)',  i, k, dom%nut(i,k), dom%td_eps(i, k), dom%tke(i, k)

          enddo
        enddo

      end if

      ! update eddy viscousity coefficient
      do i = 1, nx
        do k = dom%kc_bot(i), dom%k_srf-1
          dom%tke(i,k) = MAX(dom%tke(i,k), 1.25D-7)
          dom%td_eps(i,k) = MAX(dom%td_eps(i,k), 1.0e-9)
          ! WARNING forced control maximum tke due to prevent a sudden increase
          dom%tke(i,k) = MIN(dom%tke(i,k), 10.0)
          ! WARNING forced control maximum turbulent dissipation eps due to prevent a sudden increase
          !dom%td_eps(i,k) = MIN(dom%td_eps(i,k), 1.0)

          dom%nut(i,k) = min(0.09*dom%tke(i,k)**2/max(dom%td_eps(i,k),1.0d-9), 0.2)
        !print '(2i5, 3f30.10)',  i, k, dom%nut(i,k), dom%td_eps(i, k), dom%tke(i, k)
        enddo
      enddo

      ! --- new dmz, dhz, dcz ---
      do i = 1, nx
        do k = dom%kc_bot(i), dom%k_srf-1
          dom%dmz(i,k) = 0.5*(dom%nut(i,k) + dom%nut(i,k-1))
          dom%dmz(i,k) = MAX(1.4e-6, dom%dmz(i,k))
          dom%dmz(i,k) = MIN(2.0e-1, dom%dmz(i,k))
          dom%dhz(i,k) = MAX(1.4e-7, 0.14*dom%dmz(i,k))
          dom%dcz(i,k) = MAX(1.4e-7, 0.14*dom%dmz(i,k))

          ! mixing for unstable zone
          if (dom%rho(i,k+1) > dom%rho(i,k)) then
          !  dom%dmz(i,k) = dom%dmix
            dom%dhz(i,k) = dom%dmix
            dom%dcz(i,k) = dom%dmix
          endif

        ! print '(2i5, 3f30.10)',  i, k, dom%dmz(i, k)
        enddo
      enddo

      do i=1, nx
        dom%dmz(i,dom%kc_bot(i)-1) = dom%dmz(i,dom%kc_bot(i))
        dom%dhz(i,dom%kc_bot(i)-1) = dom%dhz(i,dom%kc_bot(i))
        dom%dcz(i,dom%kc_bot(i)-1) = dom%dcz(i,dom%kc_bot(i))

        dom%dmz(i,dom%k_srf) = dom%dmz(i,dom%k_srf-1)
        dom%dhz(i,dom%k_srf) = dom%dhz(i,dom%k_srf-1)
        dom%dcz(i,dom%k_srf) = dom%dcz(i,dom%k_srf-1)
      enddo

      if (dom%freeslip) then
        do i=1, nx
          dom%dmz(i,dom%kc_bot(i)) = 0.0
          dom%dhz(i,dom%kc_bot(i)) = 0.0
          dom%dcz(i,dom%kc_bot(i)) = 0.0

          dom%dmz(i,dom%k_srf) = 0.0
          dom%dhz(i,dom%k_srf) = 0.0
          dom%dcz(i,dom%k_srf) = 0.0
        enddo        
      endif 

    else
      ! --- from Richardson number function ---
      ! vertical diffusivity
      dom%dmz(:,:) = 0.0
      dom%dhz(:,:) = 0.0
      dom%dcz(:,:) = 0.0

      do i=1, nx
        do k=dom%kc_bot(i), dom%k_srf-1

          ! temperature, deg-C
          t_ = (dom%t(i,k) + dom%t(i,k+1))/2

          ! kinematic viscosity, m2/s
          nu_ = water_viscosity(t_)/water_density(t_)

          ! thermal diffusivity, m2/s
          alpha_ = water_thermal_conductivity(t_)/(rho_water*c_water*4186)

          ! concentration diffusivity, m2/s
          beta_ = nu_

          ! local Richardson number
          u1_ = (dom%u(i-1,k) + dom%u(i,k))/2
          u2_ = (dom%u(i-1,k+1) + dom%u(i,k+1))/2
          uz_ = abs(u2_ - u1_)/dom%dzs(k) + 1.e-6
          rhoz_ = (dom%rho(i,k+1) - dom%rho(i,k))/dom%dzs(k)
          rho_ = (dom%rho(i,k) + dom%rho(i,k+1))/2
          ri_ = min(max(0.0, -gravity*rhoz_/(rho_*uz_**2)), 15.)

          ! eddy viscosity, m2/s
          nut_ = dom%dmz0*exp(-dom%ll*ri_)

          ! eddy thermal diffusivity, m2/s
          alphat_ = dom%dhz0*exp(-dom%mm*ri_)

          ! eddy concentration diffusivity, m2/s
          betat_ = dom%dcz0*exp(-dom%nn*ri_)

          ! effective diffusivity
          dom%dmz(i,k) = nu_ + nut_
          dom%dhz(i,k) = alpha_ + alphat_
          dom%dcz(i,k) = beta_ + betat_

          ! mixing for unstable zone
          if (dom%rho(i,k+1) > dom%rho(i,k)) then
          !  dom%dmz(i,k) = dom%dmix
            dom%dhz(i,k) = dom%dmix
            dom%dcz(i,k) = dom%dmix
          endif

        enddo

        dom%dmz(i,dom%kc_bot(i)-1) = dom%dmz(i,dom%kc_bot(i))
        dom%dhz(i,dom%kc_bot(i)-1) = dom%dhz(i,dom%kc_bot(i))
        dom%dcz(i,dom%kc_bot(i)-1) = dom%dcz(i,dom%kc_bot(i))

        dom%dmz(i,dom%k_srf) = dom%dmz(i,dom%k_srf-1)
        dom%dhz(i,dom%k_srf) = dom%dhz(i,dom%k_srf-1)
        dom%dcz(i,dom%k_srf) = dom%dhz(i,dom%k_srf-1)


      enddo
    endif

    deallocate (depth)
    deallocate (u, w)

  end subroutine update_turbulence


  subroutine write_exceed_point(dom, it_step, time_day)
    !
    ! write log when exceeds threshold
    !
    implicit none

    type(domain), intent(in) :: dom
    integer, intent(in) :: it_step
    real, intent(in) :: time_day

    real, parameter :: nut_threshold = 10.0
    real, parameter :: dmz_threshold = 0.01
    real, parameter :: PGX_threshold = 1.0
    integer, parameter :: unit_no = 9100
    logical, save :: first_call = .true.
    character(len=128), save :: fname

    integer :: i, k, nx, k_srf

    if (first_call) then
      fname = 'out/' // trim(dom%name) // '/exceed_points.csv'
      open(unit_no, file=fname, status='unknown')
      write(unit_no,'(a)') &
        'STEP, DAYS, I, K, U, W, P, NUT, TKE, EPS, DMX, DMZ, PGX'
      first_call = .false.
    endif

    nx = dom%nx; k_srf = dom%k_srf

    do i=1, nx
      do k=dom%kc_bot(i), k_srf-1
        if (dom%nut(i,k) <= nut_threshold) cycle
        if (.not. dom%exceed_logged(i,k)) then
          write(unit_no,'(i8,a,f10.4,a,i6,a,i6,a,e15.7,a,e15.7,a,e15.7,a,e15.7,a,e15.7,a,e15.7,a,e15.7,a,e15.7,a,e15.7)') &
            it_step, ',', time_day, ',', i, ',', k, ',', dom%u(i,k), ',', dom%w(i,k), ',', dom%p(i,k), ',', &
            dom%nut(i,k), ',', dom%tke(i,k), ',', dom%td_eps(i,k), ',', dom%dmx(i,k), ',', dom%dmz(i,k), ',', dom%PGX(i,k)
          dom%exceed_logged(i,k) = .true.
        end if
      enddo
    enddo

  end subroutine write_exceed_point


end module mod_domain


module mod_boundary
  !
  ! module for boundary conditions
  !
  use mod_timeseries, only : timeseries, timeseries_read, timeseries_interp
  use mod_domain, only : domain, doms

  implicit none

contains

  subroutine set_boundary(dom)

    use mod_domain, only : get_total_volume, update_surface_layer
    use mod_parameter, only : nps

    implicit none

    type(domain), intent(inout) :: dom

    integer :: nx, nz, k_srf
    real, pointer, dimension(:) :: x, z, xc, zc, z_bed
    integer, pointer :: k_bot(:), kc_bot(:)

    integer :: i, j, k
    real :: tmp
    logical, allocatable :: east(:,:)  ! cell faces on the east side

    ! namelist
    integer :: id_up
    real :: fr_in, b_in, z_in_low
    character(len=128) :: fname_in

    integer :: n_out
    logical, dimension(10) :: surf_out
    real, dimension(10) :: z_out, fr_out, phi_out, out_height, zKTSW, zKBSW
    character(len=128) :: fname_out

    integer :: n_trb
    real, dimension(10) :: x_trb, fr_trb, b_trb, theta_trb
    real :: z_trb_low
    character(len=128) :: fname_trb(10)

    integer :: n_cnf
    integer :: id_cnf(10)
    real, dimension(10) :: x_cnf, theta_cnf

    integer :: n_wtp
    integer :: id_wtp(10)
    real, dimension(10) :: x_wtp, z_wtp, fr_wtp, phi_wtp, theta_wtp
    character(len=128) :: fname_wtp(10)

    integer :: n_pin
    real, dimension(10) :: x_pin, z_pin, fr_pin, phi_pin, theta_pin
    character(len=128) :: fname_pin(10)

    integer :: n_pout
    real, dimension(10) :: x_pout, z_pout, fr_pout, phi_pout
    character(len=128) :: fname_pout

    integer :: n_fnc
    logical :: type_fnc(10)
    real :: x_fnc(10), width_fnc(10)
    type(domain), pointer :: trb

    namelist /boundary/ &
      id_up, fr_in, b_in, z_in_low, fname_in, &
      n_out, surf_out, z_out, fr_out, phi_out, fname_out, out_height, zKTSW, zKBSW,&
      n_trb, x_trb, fr_trb, b_trb, theta_trb, z_trb_low, fname_trb, &
      n_cnf, id_cnf, x_cnf, theta_cnf, &
      n_wtp, id_wtp, x_wtp, z_wtp, fr_wtp, phi_wtp, theta_wtp, fname_wtp, &
      n_pin, x_pin, z_pin, fr_pin, phi_pin, theta_pin, fname_pin, &
      n_pout, x_pout, z_pout, fr_pout, phi_pout, fname_pout, &
      n_fnc, type_fnc, x_fnc, width_fnc

    ! default parameters
    id_up=0; fr_in=0.25; b_in=100.; z_in_low=-999.; fname_in=''
    n_out=1; surf_out=.false.; z_out=-999.; fr_out=0.134; phi_out=90.; fname_out=''
    out_height= 4.0; zKTSW= 999.0; zKBSW=-999.0
    n_trb=0; x_trb=-999.; fr_trb=0.25; b_trb=100.; theta_trb=90.; z_trb_low=-999.
    fname_trb=''
    n_cnf=0; id_cnf=0; x_cnf=-999.; theta_cnf=90.
    n_wtp=0; id_wtp=0; x_wtp=-999.; z_wtp=-999.; fr_wtp=0.134; phi_wtp=90.
    theta_wtp=90.; fname_wtp=''
    n_pin=0; x_pin=-999.; z_pin=-999.; fr_pin=0.134; phi_pin=90.
    theta_pin=90.; fname_pin=''
    n_pout=0; x_pout=-999.; z_pout=-999.; fr_pout=0.134; phi_pout=90.; fname_pout=''
    n_fnc=0; type_fnc=.false.; x_fnc=-999.; width_fnc=10.

    ! read namelist parameters
    open(10, file=trim(dom%fname_nml), status='old')
    read(10, nml=boundary)
    close(10)

    ! write namelsit parameters
    write(6, '(a)') trim(dom%name)
    write(6, nml=boundary)

    !--- pointer ----------------------------------------
    nx=dom%nx;  nz=dom%nz
    x=>dom%x; z=>dom%z; xc=>dom%xc; zc=>dom%zc; z_bed=>dom%z_bed
    k_bot=>dom%k_bot; kc_bot=>dom%kc_bot; k_srf=dom%k_srf
    !----------------------------------------------------

    !-- inlet
    !
    allocate (dom%c_in(nps), dom%ts_in, dom%i_inlet(nz))
#ifdef SCALAR
    allocate (dom%s_in)
#endif

    ! find inlet cell faces
    allocate (east(0:nx,nz))  ! east cell face (u-position)

    east(:,:) = .false.
    do i=1, nx-1
      do k=1, nz
        if (kc_bot(i) > k .and. k >= kc_bot(i+1)) east(i,k) = .true.
      enddo
    enddo

    do k=kc_bot(1), nz
      east(0,k) = .true.
    enddo

    dom%i_inlet(:) = nx

    do i=0, nx-1
      do k=1, nz
        if (east(i,k)) dom%i_inlet(k) = i
      enddo
    enddo

    deallocate (east)

    ! setup time series data
    if (id_up == 0) call timeseries_read(fname_in, dom%ts_in)

    ! parameter setting
    dom%id_up=id_up; dom%fr_in=fr_in; dom%b_in=b_in; dom%z_in_low=z_in_low

    !-- outlet
    !
    allocate ( &
      dom%surf_outs(n_out), dom%z_outs(n_out), dom%fr_outs(n_out), dom%phi_outs(n_out), &
      dom%out_heights(n_out), dom%zKTSWs(n_out), dom%zKBSWs(n_out), &
      dom%ts_out, dom%q_outs(n_out), dom%t_outs(n_out), dom%c_outs(n_out,nps), &
      dom%t_out, dom%c_out(nps))
#ifdef SCALAR
    allocate (dom%s_outs(n_out), dom%s_out)
#endif

    ! setup time series data
    if (n_out > 0) call timeseries_read(fname_out, dom%ts_out)

    ! parameter setting
    dom%n_out=n_out; dom%surf_outs=surf_out(1:n_out)
    dom%z_outs=z_out(1:n_out); dom%fr_outs=fr_out(1:n_out)
    dom%phi_outs=phi_out(1:n_out)
    dom%out_heights=out_height(1:n_out); dom%zKTSWs=zKTSW(1:n_out); dom%zKBSWs=zKBSW(1:n_out)

    !-- tributary
    !
    allocate ( &
      dom%q_trbs(n_trb), dom%t_trbs(n_trb), dom%c_trbs(n_trb,nps), &
      dom%fr_trbs(n_trb), dom%b_trbs(n_trb), dom%theta_trbs(n_trb), &
      dom%z_trbs(n_trb), dom%i_trbs(n_trb), dom%ts_trbs(n_trb), &
      dom%q_trb(nx,nz), dom%u_trb(nx,nz), &
      dom%t_trb(nx,nz), dom%c_trb(nx,nz,nps))
#ifdef SCALAR
    allocate (dom%s_trbs(n_trb), dom%s_trb(nx,nz))
#endif

    ! find tributary index, i
    do j=1, n_trb
      do i=1, nx
        if (x(i-1) <= x_trb(j) .and. x_trb(j) < x(i)) then
          dom%i_trbs(j) = i
          exit
        endif
      enddo
    enddo

    ! setup time series data
    do i=1, n_trb
      call timeseries_read(fname_trb(i), dom%ts_trbs(i))
    enddo

    ! parameter setting
    dom%n_trb=n_trb
    dom%fr_trbs=fr_trb(1:n_trb); dom%b_trbs=b_trb(1:n_trb)
    dom%z_trb_low=z_trb_low; dom%theta_trbs=theta_trb(1:n_trb)

    !-- confluence
    !
    allocate ( &
      dom%id_cnfs(n_cnf), &
      dom%q_cnfs(n_cnf), dom%t_cnfs(n_cnf), dom%c_cnfs(n_cnf,nps), &
      dom%theta_cnfs(n_trb), dom%i_cnfs(n_cnf), &
      dom%q_cnf(nx,nz), dom%u_cnf(nx,nz), &
      dom%t_cnf(nx,nz), dom%c_cnf(nx,nz,nps))
#ifdef SCALAR
    allocate (dom%s_cnfs(n_trb), dom%s_cnf(nx,nz))
#endif

    ! check outflow conditions in tributary domains
    do j=1, n_cnf
      trb => doms(id_cnf(j))
      if (trb%n_out > 0) then
        trb%n_out = 0
        print '(a,a,a)', 'warning: overwrite n_out in domain ', &
          trim(trb%name), ' to 0.'
      endif
    enddo

    ! find confluence index, i
    do j=1, n_cnf
      do i=1, nx
        if (x(i-1) <= x_cnf(j) .and. x_cnf(j) < x(i)) then
          dom%i_cnfs(j) = i
          exit
        endif
      enddo
    enddo

    ! parameter setting
    dom%n_cnf=n_cnf
    dom%id_cnfs=id_cnf(1:n_cnf); dom%theta_cnfs=theta_cnf(1:n_cnf)

    !-- water pipe
    !
    allocate ( &
      dom%id_wtps(n_wtp), &
      dom%q_wtps(n_wtp), dom%t_wtps(n_wtp), dom%c_wtps(n_wtp,nps), &
      dom%fr_wtps(n_wtp), dom%phi_wtps(n_wtp), dom%theta_wtps(n_wtp), &
      dom%z_wtps(n_wtp), dom%i_wtps(n_wtp), dom%ts_wtps(n_wtp), &
      dom%q_wtp(nx,nz), dom%u_wtp(nx,nz), dom%t_wtp(nx,nz), dom%c_wtp(nx,nz,nps))
#ifdef SCALAR
    allocate (dom%s_wtps(n_wtp), dom%s_wtp(nx,nz))
#endif

    ! find water pipe index, i
    do j=1, n_wtp
      do i=1, nx
        if (x(i-1) <= x_wtp(j) .and. x_wtp(j) < x(i)) then
          dom%i_wtps(j) = i
          exit
        endif
      enddo
    enddo

    ! setup time series data
    if (n_wtp > 0) then
      do i=1, n_wtp
        if (trim(fname_wtp(i)) /= '') &
          call timeseries_read(fname_wtp(i), dom%ts_wtps(i))
      enddo
    endif

    ! parameter setting
    dom%n_wtp=n_wtp; dom%id_wtps=id_wtp(1:n_wtp)
    dom%fr_wtps=fr_wtp(1:n_wtp); dom%phi_wtps=phi_wtp(1:n_wtp)
    dom%theta_wtps=theta_wtp(1:n_wtp); dom%z_wtps=z_wtp(1:n_wtp)

    !-- point inflow
    !
    allocate ( &
      dom%q_pins(n_pin), dom%t_pins(n_pin), dom%c_pins(n_pin,nps), &
      dom%fr_pins(n_pin), dom%phi_pins(n_pin), dom%theta_pins(n_pin), &
      dom%z_pins(n_pin), dom%i_pins(n_pin), dom%k_pins(n_pin), dom%ts_pins(n_pin), &
      dom%q_pin(nx,nz), dom%u_pin(nx,nz), dom%t_pin(nx,nz), dom%c_pin(nx,nz,nps))
#ifdef SCALAR
    allocate (dom%s_pins(n_pin), dom%s_pin(nx,nz))
#endif

    ! find point inflow index, i
    do j=1, n_pin
      do i=1, nx
        if (x(i-1) <= x_pin(j) .and. x_pin(j) < x(i)) then
          dom%i_pins(j) = i
          exit
        endif
      enddo
      do k=1, nz
        if (z(k-1) <= z_pin(j) .and. z_pin(j) < z(k)) then
          dom%k_pins(j) = k
          exit
        endif
      enddo
    enddo

    ! setup time series data
    if (n_pin > 0) then
      do i=1, n_pin
        if (trim(fname_pin(i)) /= '') &
          call timeseries_read(fname_pin(i), dom%ts_pins(i))
      enddo
    endif

    ! parameter setting
    dom%n_pin=n_pin
    dom%fr_pins=fr_pin(1:n_pin); dom%phi_pins=phi_pin(1:n_pin)
    dom%theta_pins=theta_pin(1:n_pin); dom%z_pins=z_pin(1:n_pin)

    !-- point outflow
    !
    allocate ( &
      dom%q_pouts(n_pout), dom%t_pouts(n_pout), dom%c_pouts(n_pout,nps), &
      dom%fr_pouts(n_pout), dom%phi_pouts(n_pout), dom%z_pouts(n_pout), &
      dom%i_pouts(n_pout), dom%k_pouts(n_pout), dom%ts_pout, dom%q_pout(nx,nz))
#ifdef SCALAR
    allocate (dom%s_pouts(n_pout), dom%s_pout(nx,nz))
#endif

    ! find point inflow index, i
    do j=1, n_pout
      do i=1, nx
        if (x(i-1) <= x_pout(j) .and. x_pout(j) < x(i)) then
          dom%i_pouts(j) = i
          exit
        endif
      enddo
      do k=1, nz
        if (z(k-1) <= z_pout(j) .and. z_pout(j) < z(k)) then
          dom%k_pouts(j) = k
          exit
        endif
      enddo
    enddo

    ! setup time series data
    if (n_pout > 0) call timeseries_read(fname_pout, dom%ts_pout)

    ! parameter setting
    dom%n_pout=n_pout
    dom%fr_pouts=fr_pout(1:n_pout); dom%phi_pouts=phi_pout(1:n_pout)
    dom%z_pouts=z_pout(1:n_pout)

    !-- fence
    !
    allocate ( &
      dom%type_fncs(n_fnc), dom%width_fncs(n_fnc), &
      dom%i_fncs(n_fnc), dom%k_fncs(2,n_fnc))

    ! position indices of fences
    do j=1, n_fnc

      do i=1, nx-1
        if (xc(i) <= x_fnc(j) .and. x_fnc(j) < xc(i+1)) then
          dom%i_fncs(j) = i
          exit
        endif
      enddo

      if (type_fnc(j)) then ! floating type
      else                  ! fixed type
        dom%k_fncs(1,j) = k_bot(dom%i_fncs(j))
        tmp = z_bed(dom%i_fncs(j)) + width_fnc(j)
        do k=1, nz-1
          if (z(k) < tmp .and. tmp <= z(k+1)) then
            dom%k_fncs(2,j) = k
            exit
          endif
        enddo
      endif

    enddo

    ! parameter setting
    dom%n_fnc=n_fnc
    dom%type_fncs=type_fnc(1:n_fnc); dom%width_fncs=width_fnc(1:n_fnc)

    !-- column flowrate, m3/s
    !
    allocate (dom%q_col(nx))
    dom%q_col=0.0

    !-- surface layer
    !
    ! update total water volume by new z_srf
    call get_total_volume(dom)

    ! update surface layer by new z_srf
    call update_surface_layer(dom)

  end subroutine set_boundary


  subroutine update_inflow(dom, time_day)
    !
    ! inputs:  ts_in, time_day, rho_avg
    ! outputs: z_in, q_in, t_in, c_in, rho_in, u, t, c
    !
    use mod_parameter, only : gravity, nps, turbid_density, rho_water

    implicit none

    type(domain), intent(inout) :: dom
    real, intent(in) :: time_day

    integer :: nz, k_srf
    real :: z_srf, rho_in, z_in
    integer, pointer :: kc_bot(:), k_bot(:)
    real, pointer :: z(:), zc(:), dz(:), rho_(:), z_bed(:)
    type(domain), pointer :: dom_up

    integer :: i, k, kc, kc_min
    real :: tmp, eps, delta, zeta
    real, allocatable :: array(:), f(:)

    character(3) :: inflowtype
    ! namelist
    namelist /inflow/ inflowtype

    ! loading from namelist file
    open(10, file='namelist.in', status='old')
    read(10, nml=inflow)
    close(10)

    !--- pointer and allocate arrays ---------------
    nz=dom%nz; k_srf=dom%k_srf; z_srf=dom%z_srf; kc_bot=>dom%kc_bot; k_bot=>dom%k_bot
    z=>dom%z; zc=>dom%zc; dz=>dom%dz; rho_=>dom%rho_avg; z_bed=>dom%z_bed
    allocate (f(dom%nz))
    !-----------------------------------------------

    if (dom%id_up > 0) then

      ! one-way coupling up- and down-stream domains.
      ! inflow of downstream domain to outflow of upstream domain
      dom_up => doms(dom%id_up)
      dom%q_in = sum(dom_up%q_outs)
      dom%t_in = dom_up%t_out
      dom%c_in = dom_up%c_out
#ifdef SCALAR
      dom%s_in = dom_up%s_out
#endif

    else

      ! individual inflows to up- and down-stream domains.
      ! set inlet values at current time
      allocate (array(dom%ts_in%nv))
      call timeseries_interp(dom%ts_in, time_day, array)
      dom%q_in = array(1)  ! inlet flow, m3/s
      dom%t_in = array(2)  ! inlet temp., deg-C
      do i=1, nps
        dom%c_in(i) = array(2+i)  ! inlet SS, g/m3
      enddo
#ifdef SCALAR
      dom%s_in = array(2+nps+1)  ! inlet scalar, -
#endif
      deallocate (array)

    endif

    ! inlet density
    rho_in = turbid_density(dom%t_in, sum(dom%c_in))

    ! inlet height as equivalent density
    kc_min = minval(kc_bot)

    kc = -1
    z_in = -999.

  SELECT CASE (inflowtype)

  CASE ('EDI')        ! Equidensity Inflow, normal distributed
    ! search equidensity zone from layer averaga density
    if (rho_in <= rho_(k_srf-1)) then
      z_in = z_srf
      kc = k_srf
    else if (rho_(kc_min) <= rho_in) then
      z_in = zc(kc_min)
      kc = kc_min
    else
      do k=k_srf-1, kc_min+1, -1
        if (rho_(k-1) >= rho_in .and. rho_in >= rho_(k)) then
          tmp = (rho_in - rho_(k))/(rho_(k-1) - rho_(k) + 1.e-8)
          z_in = tmp*z(k-1) + (1.0 - tmp)*z(k)
          kc = k
          exit
        endif
      enddo
    endif

    if (z_in < dom%z_in_low) then
      z_in = dom%z_in_low
      dom%z_in = z_in
      do k=kc_min, k_srf
        if (zc(k-1) <= z_in .and. z_in <= zc(k)) then
          kc = k
          exit
        endif
      enddo
    endif

    if (kc < 0) then
      print *, 'error: cant find kc, in subroutine update_inflow'
      stop
    endif

    ! normalized density gradient, 1/m
    eps = -(rho_(kc) - rho_(kc-1))/(z(kc) - z(kc-1))/rho_water
    eps = max(1.e-6, eps)

    ! flow depth of 2d jet, m
    delta = abs(dom%q_in) / (dom%fr_in*dom%b_in*sqrt(eps*gravity))
    delta = min(max(dz(kc), sqrt(delta)), z(nz)-z(0))

    ! velocity profile function
    f(:) = 0.0
    do k=1, k_srf
      zeta = (zc(k) - z_in)/delta
      if (-0.5 <= zeta .and. zeta <= 0.5) &
        f(k) = exp(-0.5*(zeta*3.92)**2)
    enddo

    ! inlet velocity
    tmp = 1.e-8
    do k=1, k_srf
      i = dom%i_inlet(k)
      tmp = tmp + dom%au(i,k)*f(k)
    enddo
    do k=1, k_srf
      i = dom%i_inlet(k)
      dom%u(i,k) = dom%q_in*f(k)/tmp
    enddo

    ! inlet temerature and SS
    do k=1, k_srf
      i = dom%i_inlet(k)
      dom%t(i,k) = dom%t_in
      dom%c(i,k,:) = dom%c_in(:)
#ifdef SCALAR
      dom%s(i,k) = dom%s_in
#endif
    enddo

    ! total flowrate
    dom%qtot_in = dom%q_in

  CASE ('RID')      ! Distributed River Inflow
    ! Set river inflow level at average depth
    i = dom%i_inlet(k_srf)
    z_in = z_srf - (z_srf - z_bed(i))/2
    kc = (k_srf + kc_bot(i))/2

    ! normalized density gradient, 1/m
    eps = -(rho_(kc) - rho_(kc-1))/(z(kc) - z(kc-1))/rho_water
    eps = max(1.e-6, eps)

    ! flow depth of 2d jet, m
    delta = abs(dom%q_in) / (dom%fr_in*dom%b_in*sqrt(eps*gravity))
    delta = min(max(dz(kc), sqrt(delta)), z(nz)-z(0))

    ! velocity profile function
    f(:) = 0.0
    do k=1, k_srf
      zeta = (zc(k) - z_in)/delta
      if (-0.5 <= zeta .and. zeta <= 0.5) &
        f(k) = exp(-0.5*(zeta*3.92)**2)
    enddo

    ! inlet velocity
    tmp = 1.e-8
    do k=1, k_srf
      i = dom%i_inlet(k)
      tmp = tmp + dom%au(i,k)*f(k)
    enddo
    do k=1, k_srf
      i = dom%i_inlet(k)
      dom%u(i,k) = dom%q_in*f(k)/tmp
    enddo

    ! inlet temerature and SS
    do k=1, k_srf
      i = dom%i_inlet(k)
      dom%t(i,k) = dom%t_in
      dom%c(i,k,:) = dom%c_in(:)
#ifdef SCALAR
      dom%s(i,k) = dom%s_in
#endif
    enddo

    ! total flowrate
    dom%qtot_in = dom%q_in

  CASE ('RI1')   ! 1 cell River Inflow due to density selection
    ! determine inflow layer from inlet colomn density
    kc = k_srf
    i = dom%i_inlet(k_srf)
    do while (rho_in > rho_(kc) .and. kc > k_bot(i))
      kc = kc - 1
      i = dom%i_inlet(k_srf)
    enddo
    z_in = z(kc)
    
    ! clear and update horizontal velocity into selected layer
    i = dom%i_inlet(k_srf)
    do k=kc_bot(i), k_srf
      dom%u(i,k) = 0.0
    enddo

    ! ONLY surface inflow, cell height check
    if(kc == k_srf) then
      if(dom%dz(k_srf)/dom%dz0(k_srf) < 0.7) then
        dom%u(i, k_srf) = dom%q_in*0.5/dom%au(i, k_srf)
        dom%u(i, k_srf-1) = dom%q_in*0.5/dom%au(i, k_srf-1)
      endif
    else
      dom%u(i,kc) = dom%q_in/dom%au(i,kc)
    endif

    !print '(4i4, 4f15.5)', i, dom%i_inlet(kc), k_srf, kc, dom%au(i,kc), dom%q_in, dom%u(i,kc), dom%u(i,kc-1)

    ! inlet temerature and SS
    do k=1, k_srf
      i = dom%i_inlet(k)
      dom%t(i,k) = dom%t_in
      dom%c(i,k,:) = dom%c_in(:)
#ifdef SCALAR
      dom%s(i,k) = dom%s_in
#endif
    enddo

    ! total flowrate
    dom%qtot_in = dom%q_in

  END SELECT

    !--- pointer and de-allocate arrays ---
    dom%z_in=z_in; dom%delta_in=delta/1.00; dom%rho_in=rho_in
    deallocate (f)
    !--------------------------------------

  end subroutine update_inflow

  subroutine update_outflow(dom, time_day)
    !
    ! inputs:  ts_out, time_day, rho_avg
    ! outputs: q_outs, u, t_out, c_out
    !
    use mod_parameter, only : pi, gravity, nps, rho_water, turbid_density

    implicit none

    type(domain), intent(inout) :: dom
    real, intent(in) :: time_day

    integer :: nx, nz, k_srf
    real :: z_srf
    integer, pointer :: kc_bot(:)
    real, pointer :: z(:), zc(:), dz(:), rho_(:)
    real, pointer, dimension(:,:) :: rho
    integer :: i, j, k, l, kc, kb, kc_up, kc_dw, kc_thick, KTOP, KBOT, SETT, SETB, fncnx, rhonx
    real :: tmp, eps, delta, zeta, q_, z_out, lev_outtop, lev_outbot, zKTSW, zKBSW
    real :: HT, HB, OUTCOEF, RHOFT, HSWT, RHOFB, HSWB, delta_top, delta_bot, DLRHOT, DLRHOB, DLRHOMAX, VSUM
    real, allocatable :: array(:), f(:), c_(:), u(:), q(:), rho_outavg(:), volay(:)

    character(3) :: outflowtype
    logical :: ktswkbsw = .false.   ! impose(T) or expose(F) KTSW and KBSW to calc KTOP and KBOT 
    real ::         &
      out_height= 3.0        ! outlet gate height

    ! namelist
    namelist /outflow/ outflowtype, ktswkbsw, out_height

    ! loading from namelist file
    open(10, file='namelist.in', status='old')
    read(10, nml=outflow)
    close(10)

    !--- pointer and allocate arrays ---------------
    nx=dom%nx; nz=dom%nz; k_srf=dom%k_srf; z_srf=dom%z_srf; kc_bot=>dom%kc_bot
    z=>dom%z; zc=>dom%zc; dz=>dom%dz; rho_=>dom%rho_avg; rho=>dom%rho
    allocate (f(dom%nz), c_(nps), u(dom%nz), q(dom%nz), rho_outavg(dom%nz), volay(dom%nz))
    !-----------------------------------------------

    ! bottom k at outlet
    kb = dom%k_bot(nx)

    if (dom%n_out > 0) then  ! dam with some outlets

      ! get current outlet flowrate
      allocate (array(dom%ts_out%nv))
      call timeseries_interp(dom%ts_out, time_day, array)
      do i=1, dom%ts_out%nv
        dom%q_outs(i) = array(i)  ! outlet flow, m3/s
      enddo
      deallocate (array)

      ! outlet variables
      dom%u(nx,:)=0.0


    SELECT CASE (outflowtype)

    CASE ('NDO')        ! normal distributed outflow
      do j=1, dom%n_out

        kc = -1

        if (dom%surf_outs(j)) then ! surface outflow

          z_out = z_srf
          kc = k_srf

        else                     ! in-water outflow

          z_out = dom%z_outs(j)

          ! find kc
          if (dom%z_outs(j) >= z(k_srf-1)) then
            z_out = z_srf
            kc = k_srf
          else if (z(kb) >= dom%z_outs(j)) then
            z_out = z(kb)
            kc = kb
          else
            do k=k_srf-1, kb+1, -1
              if (z(k-1) <= dom%z_outs(j) .and. dom%z_outs(j) <= z(k)) then
                z_out = dom%z_outs(j)
                kc = k
                exit
              endif
            enddo
          endif

        endif

        if (kc < 0) then
          print *, 'error: cant find kc, in subroutine update_outflow'
          stop
       endif

        ! normalized density gradient, 1/m
        eps = -(rho_(kc) - rho_(kc-1))/(z(kc) - z(kc-1))/rho_water
        eps = max(1.e-8, eps)

        ! flow depth of axi-symmetric jet, m
        delta = abs(dom%q_outs(j)) / &
          (dom%fr_outs(j)*(dom%phi_outs(j)*pi/180)*sqrt(eps*gravity))
        delta = min(max(2*dz(kc), delta**0.333), z(nz)-z(0))

        ! velocity profile function
        f = 0.0
        do k=kb, k_srf
          zeta = (zc(k) - z_out)/delta
          if (-0.5 <= zeta .and. zeta <= 0.5) &
            f(k) = exp(-0.5*(zeta*3.92)**2)
        enddo

        ! outlet velocity, temperature and suspended solids
        tmp = sum(dom%au(nx,kb:k_srf)*f(kb:k_srf))
        f(kb:k_srf) = f(kb:k_srf)/tmp
        u(kb:k_srf) = dom%q_outs(j)*f(kb:k_srf)
        q(kb:k_srf) = u(kb:k_srf)*dom%au(nx,kb:k_srf)
        dom%u(nx,kb:k_srf) = dom%u(nx,kb:k_srf) + u(kb:k_srf)
        q_ = dom%q_outs(j) + 1.e-8
        dom%t_outs(j) = sum(q(kb:k_srf)*dom%t(nx,kb:k_srf))/q_
        do i=1, nps
          dom%c_outs(j,i) = sum(q(kb:k_srf)*dom%c(nx,kb:k_srf,i))/q_
        enddo
#ifdef SCALAR
        dom%s_outs(j) = sum(q(kb:k_srf)*dom%s(nx,kb:k_srf))/q_
#endif

      enddo             ! End outflow number j normal distributed in outflow height

    CASE ('EDO')        ! equal distributed in outflow height
      do j=1, dom%n_out
        kc = -1

        if (dom%surf_outs(j)) then ! surface outflow

          z_out = z_srf
          kc = k_srf

        else                     ! in-water outflow

          z_out = dom%z_outs(j)

          ! find kc
          if (dom%z_outs(j) >= z(k_srf-1)) then
            z_out = z_srf
            kc = k_srf
          else if (z(kb) >= dom%z_outs(j)) then
            z_out = z(kb)
            kc = kb
          else
            do k=k_srf-1, kb+1, -1
              if (z(k-1) <= dom%z_outs(j) .and. dom%z_outs(j) <= z(k)) then
                z_out = dom%z_outs(j)
                kc = k
                exit
              endif
            enddo
          endif

        endif

        if (kc < 0) then
          print *, 'error: cant find kc, in subroutine update_outflow'
          stop
        endif

        ! check distribution thickness
        lev_outtop = z_out + out_height*0.5
        lev_outbot = z_out - out_height*0.5

        do k=k_srf-1, kb+1, -1
          if (z(k-1) <= lev_outtop .and. lev_outtop <= z(k)) then
            kc_up = k
            exit
          endif
        enddo

        do k=k_srf-1, kb+1, -1
          if (z(k-1) <= lev_outbot .and. lev_outbot <= z(k)) then
            kc_dw = k
            exit
          endif
        enddo

        kc_thick = MAX((kc_up - kc_dw), 1)

        do l= 1, int(out_height)
          if (kc_up > k_srf) then
            kc_up = kc_up - 1
            kc_thick = kc_thick - 1
          endif
        enddo

        do l= 1, int(out_height)
          if (kc_dw < kc_bot(nx)) then
            kc_dw = kc_dw - 1
            kc_thick = kc_thick - 1
          endif
        enddo

        kc_thick = MAX(kc_thick, 1)  

        ! outlet velocity, temperature and suspended solids
        ! equal distribution
        u(:) = 0.0
        if(kc_thick > 1) then
          do k= kc_dw, kc_up
            u(k) = dom%q_outs(j)/dom%au(nx,k)/kc_thick
          enddo
        else
          u(kc) = dom%q_outs(j)/dom%au(nx,kc)
        endif

        !print '(3i4, 4f15.5)',  kc_dw, kc_up, kc_thick, dom%q_outs(j), dom%au(nx,k), u(kc_dw), u(kc)

        q(kb:k_srf) = u(kb:k_srf)*dom%au(nx,kb:k_srf)
        dom%u(nx,kb:k_srf) = dom%u(nx,kb:k_srf) + u(kb:k_srf)
        q_ = dom%q_outs(j) + 1.e-8
        dom%t_outs(j) = sum(q(kb:k_srf)*dom%t(nx,kb:k_srf))/q_
        do i=1, nps
          dom%c_outs(j,i) = sum(q(kb:k_srf)*dom%c(nx,kb:k_srf,i))/q_
        enddo
#ifdef SCALAR
        dom%s_outs(j) = sum(q(kb:k_srf)*dom%s(nx,kb:k_srf))/q_
#endif

      enddo             ! End outflow number j loop equal distributed in outflow height

    CASE ('DDD')        ! Density denpended distributed outflow
      ! average density include SS in front of dam
      volay(:) = 1.0e-10
      rho_outavg(:) = 0.0

      ! fence position check
      fncnx = dom%i_fncs(1)
      do j = 2, dom%n_fnc
        if(dom%i_fncs(1) < dom%i_fncs(j)) then
          fncnx = dom%i_fncs(j)
        end if
      end do

      ! average density near dam calc position check
      rhonx = nx-5 
      if(rhonx <= fncnx) then
        rhonx = fncnx + 1
      endif

      print '(8i4, 4f15.5)', dom%i_fncs(1), dom%i_fncs(2), dom%n_fnc, fncnx, rhonx, nx-5, nx

      ! average density near dam calc each layer
      do k = 1, k_srf
        do i=rhonx, nx
          if(k >= dom%kc_bot(i)) then
            volay(k) = volay(k) + dom%vol(i,k)
            rho_outavg(k) = rho_outavg(k) + dom%vol(i,k)*turbid_density(dom%t(i,k), sum(dom%c(i,k,:)))
          end if

          !print '(4i4, 7f15.5)', i, k, k_srf, kc_bot(nx), volay(k), rho_outavg(k)
        end do
      end do

      do k=1, k_srf
        rho_outavg(k) = rho_outavg(k)/volay(k)
      enddo

      ! find kc
      do j=1, dom%n_out

        kc = -1

        if (dom%surf_outs(j)) then ! surface outflow

          z_out = z_srf
          kc = k_srf

        else                     ! in-water outflow

          z_out = dom%z_outs(j)

          ! find kc
          if (dom%z_outs(j) >= z(k_srf-1)) then
            z_out = z_srf
            kc = k_srf
          else if (z(kb) >= dom%z_outs(j)) then
            z_out = z(kb)
            kc = kb
          else
            do k=k_srf-1, kb+1, -1
              if (z(k-1) <= dom%z_outs(j) .and. dom%z_outs(j) <= z(k)) then
                z_out = dom%z_outs(j)
                kc = k
                exit
              endif
            enddo
          endif

        endif

        if (kc < 0) then
          print *, 'error: cant find kc, in subroutine update_outflow'
          stop
        endif

        ! density based velocity profile
        f(:) = 0.0
        HSWT = 0.0
        HSWB = 0.0

        ! user diffinition KTSW and KBSW
        if (ktswkbsw) then
          KTOP = k_srf
          KBOT = kc_bot(nx)

          zKTSW = dom%zKTSWs(j)
          zKTSW = MIN(zKTSW, z_srf)
          do k= kc, k_srf
            if (zc(k) >= zKTSW) then
              KTOP = k-1; exit
            end if
          end do

          zKBSW = dom%zKBSWs(j)
          zKBSW = max(zKBSW, z(kc_bot(nx)))
          do k= kc_bot(nx), kc
            if (zKBSW <= zc(k)) then
              KBOT = k-1; exit
            end if
          end do

          !print '(4i4, 7f15.5)', i, k, KBOT, KTOP, zKBSW, ZKTSW
          
        else   ! if no input KTSW and KBSW, KTOP and KBOT calc. from flow depth of axi-symmetric jet
          if(abs(rho_outavg(k_srf) - rho_outavg(kc_bot(nx))) > 0.2) then
            eps = -(rho_outavg(kc+1) - rho_outavg(kc-1))/(z(kc+1) - z(kc-1))/rho_water
          else 
            eps = -(rho_outavg(k_srf) - rho_outavg(kc_bot(nx)))/(z(k_srf) - z(kc-1))/rho_water
          endif
            eps = max(1.e-8, eps)

          ! flow depth of axi-symmetric jet, m
          delta = (abs(dom%q_outs(j)) / (dom%fr_outs(j)*(dom%phi_outs(j)*pi/180)*sqrt(eps*gravity)))**0.33333
          delta = min(max(2*dz(kc), delta), z(nz)-z(0))

          KTOP = kc+1
          do k= kc, k_srf
            if (zc(k) >= z_out + delta/2) then
              KTOP = k; exit
            end if
          end do

          KBOT = kc-1
          do k= kc_bot(nx), kc
            if (z_out - delta/2 <= zc(k)) then
              KBOT = k; exit
            end if
            !print '(4i4, 7f15.5)', i, k, kc, KBOT, zc(k), z_out - delta/2
          end do
        endif

        ! coefficient for outflow thickness near the water surface or bottom
        OUTCOEF = 1.0 
        if ((z_out/z_srf) > 0.9) then
          OUTCOEF = 2.0
        endif

        ! outflow zone above structure
        do k= kc-1, KTOP
          ! Density frequency
          HT = z(k) - z_out
          RHOFT = MAX( SQRT((ABS(rho_outavg(k) - rho_outavg(kc)))/(HT*rho_outavg(kc) + 1.0E-10)*gravity), 1.0E-10 )

          ! Thickness of point sink flow
          HSWT = (OUTCOEF*dom%q_outs(j)/RHOFT)**0.333333
          ! Thickness of distributed sink flow
          !HSWT = SQRT(2.0*OUTCOEF*dom%q_outs(j)/(out_height*RHOFT))

          if (HT >= HSWT) then
            KTOP = k; exit
          end if
        
          !print '(4i4, 7f15.5)', k, k_srf, kc_bot(nx), KTOP, z(k), z_out, rho_outavg(k), rho_outavg(kc), HT, RHOFT, HSWT

        enddo

        ! Reference density of above structure
        if ((z_out + HSWT) < z_srf) THEN
          DLRHOT = ABS(rho_outavg(kc) - rho_outavg(KTOP))
          SETT=1
        else if (z_srf == z_out) THEN
          DLRHOT = 1.0e-10
          SETT=2
        else
          DLRHOT = ABS(rho_outavg(kc) - rho_outavg(k_srf))*HSWT/(z_srf-z_out)
          SETT=3
        endif
        DLRHOT = MAX(DLRHOT, 1.0e-10)

        ! outflow zone below structure
        do k= kc-1, kc_bot(nx), -1
          ! Density frequency
          HB = z_out - z(k)
          RHOFB = MAX( SQRT((ABS(rho_outavg(k) - rho_outavg(kc)))/(HB*rho_outavg(kc) + 1.0E-10)*gravity), 1.0E-10 )

          ! Thickness of point sink flow
          HSWB = (OUTCOEF*dom%q_outs(j)/RHOFB)**0.333333
          ! Thickness of distributed sink flow
          !HSWB = SQRT(2.0*OUTCOEF*dom%q_outs(j)/(out_height*RHOFB))

          if (HB >= HSWB) then
            KBOT = k; exit
          end if
          
          !print '(4i4, 7f15.5)', k, k_srf, kc_bot(nx), KBOT, z(k), z_out, rho(nx,k), rho(nx,kc), HB, RHOFB, HSWB
        enddo

        ! Reference density of below structure
        if ((z_out - HSWB) > z(KBOT-1)) then
          DLRHOB = ABS(rho_outavg(kc) - rho_outavg(KBOT))
          SETB=1
        else if (z(KBOT-1) == z_out) then
          DLRHOB = 1.0e-10
          SETB=2
        else
          DLRHOB = ABS(rho_outavg(kc) - rho_outavg(KBOT))*HSWB/(z_out - z(KBOT-1))
          SETB=3
        endif
        DLRHOB = MAX(DLRHOB, 1.0e-10)

        !print '(3i4, 7f10.5)', k, kc_bot(nx), KBOT, OUTCOEF, HB, RHOFB, rho(nx,k), rho(nx,kc), HSWB, DLRHOB

        ! Velocity profile 
        VSUM = 0.0
        do k=KBOT, KTOP
          if(k > kc) then
            DLRHOMAX = MAX(DLRHOT, 1.0e-10)   
          else
            DLRHOMAX = MAX(DLRHOB, 1.0e-10)
          endif

          f(k) = 1.0 - ((rho_outavg(kc) - rho_outavg(k))/DLRHOMAX)**2
          IF(f(k).GT.1.0) f(k)=1.0  
          IF(f(k).LT.0.0) f(k)=0.0
          !print '(5i4, 7f11.5)', k, KBOT, KTOP, SETT, SETB, DLRHOB, DLRHOT, DLRHOMAX, rho(nx,kc), rho(nx,k), f(k), ((rho(nx,kc) - rho(nx,k))/DLRHOMAX)**2

          f(k) = f(k)*dom%au(nx,k)
          VSUM = VSUM + f(k)
        end do

        ! outlet velocity, temperature and suspended solids
        !u(kb:k_srf) = dom%q_outs(j)*f(kb:k_srf)/dom%au(nx,kb:k_srf)
        u(kb:k_srf) = dom%q_outs(j)*f(kb:k_srf)/VSUM/dom%au(nx,kb:k_srf)
        q(kb:k_srf) = dom%q_outs(j)*f(kb:k_srf)
        dom%u(nx,kb:k_srf) = dom%u(nx,kb:k_srf) + u(kb:k_srf)
        q_ = dom%q_outs(j) + 1.e-8
        dom%t_outs(j) = sum(q(kb:k_srf)*dom%t(nx,kb:k_srf))/q_
        do i=1, nps
          dom%c_outs(j,i) = sum(q(kb:k_srf)*dom%c(nx,kb:k_srf,i))/q_
        enddo
#ifdef SCALAR
        dom%s_outs(j) = sum(q(kb:k_srf)*dom%s(nx,kb:k_srf))/q_
#endif

      enddo            ! End DO loop Density denpended distributed outflow

    END SELECT

    else if (dom%n_out == 0) then  ! open boundary

      dom%u(nx,kb:k_srf) = dom%u(nx-1,kb:k_srf)
      dom%t(nx+1,kb:k_srf) = dom%t(nx,kb:k_srf)
      do l=1, nps
        dom%c(nx+1,kb:k_srf,l) = dom%c(nx,kb:k_srf,l)
      enddo
#ifdef SCALAR
      dom%s(nx+1,kb:k_srf) = dom%s(nx,kb:k_srf)
#endif

    endif

    ! outlet temperature, suspended solids
    q = dom%au(nx,:)*dom%u(nx,:)
    tmp = sum(q(kb:k_srf)) + 1e-8
    dom%t_out = sum(q(kb:k_srf)*dom%t(nx,kb:k_srf))/tmp
    do i=1, nps
      dom%c_out(i) = sum(q(kb:k_srf)*dom%c(nx,kb:k_srf,i))/tmp
    enddo
#ifdef SCALAR
    dom%s_out = sum(q(kb:k_srf)*dom%s(nx,kb:k_srf))/tmp
#endif

    ! total flowrate
    if (dom%n_out == 0) then
      kb = dom%kc_bot(nx)
      dom%qtot_out = sum(dom%au(nx-1,kb:k_srf)*dom%u(nx-1,kb:k_srf))
    else
      dom%qtot_out = sum(dom%q_outs)
    endif


    !--- de-allocate arrays ---
    deallocate (f, c_, u, q, rho_outavg, volay)
    !--------------------------

  end subroutine update_outflow


  subroutine update_tributary(dom, time_day)
    !
    ! inputs:  ts_trbs, time_day, rho_avg
    ! outputs: q_trb, u_trb, t_trb, c_trb
    !
    use mod_parameter, only : pi, gravity, nps, turbid_density, rho_water

    implicit none

    type(domain), intent(inout) :: dom
    real, intent(in) :: time_day

    integer :: nz, k_srf
    real :: z_srf, rho_trb, z_trb
    real, pointer :: dx(:), z(:), zc(:), dz(:), rho_(:)

    integer :: i, j, k, l, kc, kb
    real :: tmp, eps, zeta, delta, cos_
    real, allocatable :: array(:), f(:), q(:)

    !--- pointer and allocate arrays ---------------
    nz=dom%nz; k_srf=dom%k_srf; z_srf=dom%z_srf
    dx=>dom%dx; z=>dom%z; zc=>dom%zc; dz=>dom%dz
    rho_=>dom%rho_avg
    allocate (f(dom%nz), q(dom%nz))
    !-----------------------------------------------

    ! get current tributary variables
    do i=1, dom%n_trb
      allocate (array(dom%ts_trbs(i)%nv))
      call timeseries_interp(dom%ts_trbs(i), time_day, array)
      dom%q_trbs(i) = array(1)  ! tributary flow, m3/s
      dom%t_trbs(i) = array(2)  ! tributary temp., deg-C
      do j=1, nps
        dom%c_trbs(i,j) = array(2+j)  ! tributary SS, g/m3
      enddo
#ifdef SCALAR
      dom%s_trbs(i) = array(2+nps+1)  ! tributary scalar, -
#endif
      deallocate (array)
    enddo

    ! tributary source term for momentum-eq, heat-eq and SS-eq
    dom%q_trb=0.0; dom%u_trb=0.0; dom%t_trb=0.0; dom%c_trb=0.0
#ifdef SCALAR
    ! tributary source term for scalar-eq
    dom%s_trb=0.0
#endif

    do j=1, dom%n_trb

      ! tributary position index
      i = dom%i_trbs(j)
      kb = dom%kc_bot(i)

      ! tributary flow density
      rho_trb = turbid_density(dom%t_trbs(j), sum(dom%c_trbs(j,:)))

      ! tributary height as equivalent density
      kc = -1

      if (rho_trb <= rho_(k_srf-1)) then
        z_trb = z_srf
        kc = k_srf
      else if (rho_(kb) <= rho_trb) then
        z_trb = zc(kb)
        kc = kb
      else
        do k=k_srf-1, kb+1, -1
          if (rho_(k-1) >= rho_trb .and. rho_trb >= rho_(k)) then
            tmp = (rho_trb - rho_(k))/(rho_(k-1) - rho_(k) + 1.e-8)
            z_trb = tmp*z(k-1) + (1.0 - tmp)*z(k)
            kc = k
            exit
          endif
        enddo
      endif

      if (kc < 0) then
        print *, 'error: cant find kc, in subroutine update_tributary'
        stop
      endif

      dom%z_trbs(j) = z_trb

      ! normalized density gradient, 1/m
      eps = -(rho_(kc) - rho_(kc-1))/(z(kc) - z(kc-1))/rho_water
      eps = max(1.e-6, eps)

      ! flow depth of 2d jet, m
      delta = abs(dom%q_trbs(j))/(dom%fr_trbs(j)*dom%b_trbs(j)*sqrt(eps*gravity))
      delta = min(max(2*dz(kc), sqrt(delta)), z(nz)-z(0))

      ! velocity profile function
      f = 0.0
      do k=1, k_srf
        zeta = (zc(k) - dom%z_trbs(j))/delta
        if (-0.5 <= zeta .and. zeta <= 0.5) &
          f(k) = exp(-0.5*(zeta*3.92)**2)
      enddo
      f(kb:k_srf) = f(kb:k_srf)/sum(f(kb:k_srf))

      ! source term
      cos_ = cos(dom%theta_trbs(j)*pi/180)
      q(kb:k_srf) = dom%q_trbs(j)*f(kb:k_srf)

      dom%q_trb(i,kb:k_srf) = q(kb:k_srf)
      dom%u_trb(i,kb:k_srf) = cos_*(q(kb:k_srf)/(dz(kb:k_srf)*dx(i)))
      dom%t_trb(i,kb:k_srf) = dom%t_trbs(j)
      do l=1, nps
        dom%c_trb(i,kb:k_srf,l) = dom%c_trbs(j,l)
      enddo
#ifdef SCALAR
      dom%s_trb(i,kb:k_srf) = dom%s_trbs(j)
#endif
    enddo

    ! total flowrate
    dom%qtot_trb = sum(dom%q_trbs)

    !--- de-allocate arrays ---
    deallocate (f, q)
    !--------------------------

  end subroutine update_tributary


  subroutine update_confluence(dom)
    !
    ! inputs:  dom%id_cnfs, time_day
    ! outputs: q_cnf, u_cnf, t_cnf, c_cnf, q_cnfs, t_cnfs, c_cnfs
    !
    use mod_parameter, only : pi, nps

    implicit none

    type(domain), intent(inout) :: dom

    integer :: i, j, l, kb, k_srf, nx
    real :: cos_, qtot_
    real, allocatable :: q(:)
    type(domain), pointer :: trb

    !--- pointer and allocate arrays ---------------
    allocate (q(dom%nz))
    !-----------------------------------------------

    ! confluence source term for momentum-eq, heat-eq and SS-eq
    dom%q_cnf=0.0; dom%u_cnf=0.0; dom%t_cnf=0.0; dom%c_cnf=0.0
#ifdef SCALAR
    dom%s_cnf=0.0
#endif

    do j=1, dom%n_cnf

      ! upstream tributary domain
      trb => doms(dom%id_cnfs(j))

      ! confluence position index
      i = dom%i_cnfs(j)

      ! check mesh consistency
      if (dom%nz /= trb%nz .or. dom%kc_bot(i) /= trb%kc_bot(trb%nx)) then
        print *, 'error: the confluent domains mesh is not consistent.'
        print *, 'dom%nz, trb%nz = ', dom%nz, trb%nz
        print *, 'i, dom%kc_bot(i), trb%kc_bot(trb%nx) = ', &
          i, dom%kc_bot(i), trb%kc_bot(trb%nx)
        stop
      endif

      ! confluecne angle
      cos_ = cos(dom%theta_cnfs(j)*pi/180)

      nx = trb%nx
      kb = dom%kc_bot(i)
      k_srf = dom%k_srf

      q(kb:k_srf) = trb%u(nx-1,kb:k_srf)*trb%au(nx-1,kb:k_srf)
      qtot_ = sum(q(kb:k_srf)) + 1e-8

      dom%q_cnf(i,kb:k_srf) = q(kb:k_srf) !*(trb%q_col(dom%nx)/qtot_)
      dom%u_cnf(i,kb:k_srf) = cos_*(q(kb:k_srf)/(dom%dz(kb:k_srf)*dom%dx(i)))
      dom%t_cnf(i,kb:k_srf) = trb%t(nx,kb:k_srf)
      do l=1, nps
        dom%c_cnf(i,kb:k_srf,l) = trb%c(nx,kb:k_srf,l)
      enddo
#ifdef SCALAR
      dom%s_cnf(i,kb:k_srf) = trb%s(nx,kb:k_srf)
#endif

      dom%q_cnfs(j) = qtot_
      dom%t_cnfs(j) = sum(q(kb:k_srf)*dom%t_cnf(i,kb:k_srf))/qtot_
      do l=1, nps
        dom%c_cnfs(j,l) = sum(q(kb:k_srf)*dom%c_cnf(i,kb:k_srf,l))/qtot_
      enddo
#ifdef SCALAR
      dom%t_cnfs(j) = sum(q(kb:k_srf)*dom%t_cnf(i,kb:k_srf))/qtot_
#endif
    enddo

    ! outflow boundary conditions of tributary domain
    do j=1, dom%n_cnf
      trb => doms(dom%id_cnfs(j))
      kb = trb%k_bot(trb%nx)
      k_srf = trb%k_srf

      trb%t(trb%nx+1,kb:k_srf) = dom%t(dom%nx,kb:k_srf)
      trb%c(trb%nx+1,kb:k_srf,:) = dom%c(dom%nx,kb:k_srf,:)
#ifdef SCALAR
      trb%s(trb%nx+1,kb:k_srf) = dom%s(dom%nx,kb:k_srf)
#endif
    enddo

    ! total flowrate
    dom%qtot_cnf = sum(dom%q_cnfs)

    !--- de-allocate arrays ---
    deallocate (q)
    !--------------------------

  end subroutine update_confluence


  subroutine update_waterpipe(dom, time_day)
    !
    ! inputs:  ts_wtps, time_day, rho_avg
    ! outputs: q_wtp, u_wtp, t_wtp, c_wtp
    !
    use mod_parameter, only : pi, gravity, nps, rho_water

    implicit none

    type(domain), intent(inout) :: dom
    real, intent(in) :: time_day

    integer :: nz, k_srf
    real :: z_srf
    real, pointer :: z(:), zc(:), dz(:), rho_(:)
    type(domain), pointer :: wtp

    integer :: i, j, k, l, kc, kb
    real :: tmp, eps, zeta, delta, cos_
    real, allocatable :: array(:), f(:,:), q(:)

    !--- pointer and allocate arrays ---------------
    nz=dom%nz; k_srf=dom%k_srf; z_srf=dom%z_srf
    z=>dom%z; zc=>dom%zc; dz=>dom%dz; rho_=>dom%rho_avg
    allocate (f(dom%n_wtp,dom%nz), q(dom%nz))
    !-----------------------------------------------

    do j=1, dom%n_wtp

      if (dom%id < dom%id_wtps(j)) then

        ! get current flowrate
        allocate (array(dom%ts_wtps(j)%nv))
        call timeseries_interp(dom%ts_wtps(j), time_day, array)
        dom%q_wtps(j) = array(1)  ! water pipe flow, m3/s
        deallocate (array)

      else

        ! set flowrate as negative of upstreams value
        wtp => doms(dom%id_wtps(j))  ! water pipe connection domain
        dom%q_wtps(j) = -wtp%q_wtps(j)

      endif

    enddo

    ! water pipe source term for momentum-eq, heat-eq and SS-eq
    dom%q_wtp=0.0; dom%u_wtp=0.0; dom%t_wtp=0.0; dom%c_wtp=0.0
#ifdef SCALAR
    dom%s_wtp=0.0
#endif

    do j=1, dom%n_wtp

      i = dom%i_wtps(j)    ! position index of water pipe
      kb = dom%kc_bot(i)   ! bottom k of water pipe

      ! find kc
      kc = -1

      if (dom%z_wtps(j) >= z(k_srf-1)) then
        kc = k_srf
      else if (z(kb) >= dom%z_wtps(j)) then
        kc = 1
      else
        do k=k_srf-1, kb+1, -1
          if (z(k-1) <= dom%z_wtps(j) .and. dom%z_wtps(j) <= z(k)) then
            kc = k
            exit
          endif
        enddo
      endif

      if (kc < 0) then
        print *, 'error: cant find kc, in subroutine update_waterpipe'
        stop
      endif

      ! normalized density gradient, 1/m
      eps = -(rho_(kc) - rho_(kc-1))/(z(kc) - z(kc-1))/rho_water
      eps = max(1.e-6, eps)

      ! flow depth of axi-symmetric jet, m
      delta = abs(dom%q_wtps(j)) / &
        (dom%fr_wtps(j)*(dom%phi_wtps(j)*pi/180)*sqrt(eps*gravity))
      delta = min(max(2*dz(kc), delta**0.333), z(nz)-z(0))

      ! velocity profile function
      f(j,:) = 0.0
      do k=kb, k_srf
        zeta = (zc(k) - dom%z_wtps(j))/delta
        if (-0.5 <= zeta .and. zeta <= 0.5) &
          f(j,k) = exp(-0.5*(zeta*3.92)**2)
      enddo

      wtp => doms(dom%id_wtps(j))           ! water pipe connection domain
      cos_ = cos(dom%theta_wtps(j)*pi/180)  ! water pipe angle

      f(j,kb:k_srf) = f(j,kb:k_srf)/sum(f(j,kb:k_srf))
      q(kb:k_srf) = dom%q_wtps(j)*f(j,kb:k_srf)

      ! source term
      dom%q_wtp(i,kb:k_srf) = q(kb:k_srf)
      dom%u_wtp(i,kb:k_srf) = cos_*(dom%q_wtps(j)/dom%aus(i,kb:k_srf))
      dom%t_wtp(i,kb:k_srf) = wtp%t_wtps(j)
      do l=1, nps
        dom%c_wtp(i,kb:k_srf,l) = wtp%c_wtps(j,l)
      enddo
#ifdef SCALAR
      dom%s_wtp(i,kb:k_srf) = wtp%s_wtps(j)
#endif
    enddo

    ! water pipe temperature and suspended solids
    do j=1, dom%n_wtp
      i = dom%i_wtps(j)
      kb = dom%kc_bot(i)

      f(j,kb:k_srf) = f(j,kb:k_srf)*dom%vol(i,kb:k_srf)
      tmp = sum(f(j,kb:k_srf))

      dom%t_wtps(j) = sum(f(j,kb:k_srf)*dom%t(i,kb:k_srf))/tmp
      do l=1, nps
        dom%c_wtps(j,l) = sum(f(j,kb:k_srf)*dom%c(i,kb:k_srf,l))/tmp
      enddo
#ifdef SCALAR
      dom%s_wtps(j) = sum(f(j,kb:k_srf)*dom%s(i,kb:k_srf))/tmp
#endif
    enddo

    ! total flowrate
    dom%qtot_wtp = sum(dom%q_wtps)

    !--- de-allocate arrays ---
    deallocate (f, q)
    !--------------------------

  end subroutine update_waterpipe


  subroutine update_pointin(dom, time_day)
    !
    ! inputs:  ts_pins, time_day, rho_avg
    ! outputs: q_pin, u_pin, t_pin, c_pin
    !
    use mod_parameter, only : pi, gravity, nps, rho_water

    implicit none

    type(domain), intent(inout) :: dom
    real, intent(in) :: time_day

    integer :: nz, k_srf
    real :: z_srf
    real, pointer :: z(:), zc(:), dz(:), rho_(:)

    integer :: i, j, k, l, kc, kb
    real :: eps, zeta, delta, cos_
    real, allocatable :: array(:), f(:,:), q(:)

    !--- pointer and allocate arrays ---------------
    nz=dom%nz; k_srf=dom%k_srf; z_srf=dom%z_srf
    z=>dom%z; zc=>dom%zc; dz=>dom%dz; rho_=>dom%rho_avg
    allocate (f(dom%n_pin,dom%nz), q(dom%nz))
    !-----------------------------------------------

    do j=1, dom%n_pin

      ! get current flowrate
      allocate (array(dom%ts_pins(j)%nv))
      call timeseries_interp(dom%ts_pins(j), time_day, array)
      dom%q_pins(j) = array(1)  ! flow rate, m3/s
      dom%t_pins(j) = array(2)  ! temperature, deg-C
      do i=1, nps
        dom%c_pins(j,i) = array(2+i)  ! SS, g/m3
      enddo
#ifdef SCALAR
      dom%s_pins(j) = array(2+nps+1)  ! scalar, -
#endif
      deallocate (array)

    enddo

    ! water pipe source term for momentum-eq, heat-eq and SS-eq
    dom%q_pin=0.0; dom%u_pin=0.0; dom%t_pin=0.0; dom%c_pin=0.0
#ifdef SCALAR
    dom%s_pin=0.0
#endif

    do j=1, dom%n_pin

      i = dom%i_pins(j)    ! position index of water pipe
      kb = dom%kc_bot(i)   ! bottom k of water pipe

      ! find kc
      kc = dom%k_pins(j)

      ! normalized density gradient, 1/m
      eps = -(rho_(kc) - rho_(kc-1))/(z(kc) - z(kc-1))/rho_water
      eps = max(1.e-6, eps)

      ! flow depth of axi-symmetric jet, m
      delta = abs(dom%q_pins(j)) / &
        (dom%fr_pins(j)*(dom%phi_pins(j)*pi/180)*sqrt(eps*gravity))
      delta = min(max(2*dz(kc), delta**0.333), z(nz)-z(0))

      ! velocity profile function
      f(j,:) = 0.0
      do k=kb, k_srf
        zeta = (zc(k) - dom%z_pins(j))/delta
        if (-0.5 <= zeta .and. zeta <= 0.5) &
          f(j,k) = exp(-0.5*(zeta*3.92)**2)
      enddo

      cos_ = cos(dom%theta_pins(j)*pi/180)  ! point inflow angle

      f(j,kb:k_srf) = f(j,kb:k_srf)/sum(f(j,kb:k_srf))
      q(kb:k_srf) = dom%q_pins(j)*f(j,kb:k_srf)

      ! source term
      dom%q_pin(i,kb:k_srf) = dom%q_pin(i,kb:k_srf) + q(kb:k_srf)
      dom%u_pin(i,kb:k_srf) = dom%u_pin(i,kb:k_srf) + cos_*(q(kb:k_srf)/dom%aus(i,kb:k_srf))
      dom%t_pin(i,kb:k_srf) = dom%t_pins(j)
      do l=1, nps
        dom%c_pin(i,kb:k_srf,l) = dom%c_pins(j,l)
      enddo
#ifdef SCALAR
      dom%s_pin(i,kb:k_srf) = dom%s_pins(j)
#endif
    enddo

    ! point inflow temperature and suspended solids
    do j=1, dom%n_pin
      i = dom%i_pins(j)
      kb = dom%kc_bot(i)

      f(j,kb:k_srf) = f(j,kb:k_srf)*dom%vol(i,kb:k_srf)
      f(j,kb:k_srf) = f(j,kb:k_srf)/sum(f(j,kb:k_srf))

      dom%t_pins(j) = sum(f(j,kb:k_srf)*dom%t(i,kb:k_srf))
      do l=1, nps
        dom%c_pins(j,l) = sum(f(j,kb:k_srf)*dom%c(i,kb:k_srf,l))
      enddo
#ifdef SCALAR
      dom%s_pins(j) = sum(f(j,kb:k_srf)*dom%s(i,kb:k_srf))
#endif
    enddo

    ! total flowrate
    dom%qtot_pin = sum(dom%q_pins)

    !--- de-allocate arrays ---
    deallocate (f, q)
    !--------------------------

  end subroutine update_pointin


  subroutine update_pointout(dom, time_day)
    !
    ! inputs:  ts_pout, time_day, rho_avg
    ! outputs: q_pout
    !
    use mod_parameter, only : pi, gravity, nps, rho_water

    implicit none

    type(domain), intent(inout) :: dom
    real, intent(in) :: time_day

    integer :: nz, k_srf
    real :: z_srf
    real, pointer :: z(:), zc(:), dz(:), rho_(:)

    integer :: i, j, k, l, kc, kb
    real :: eps, zeta, delta
    real, allocatable :: array(:), f(:,:), q(:)

    !--- pointer and allocate arrays ---------------
    nz=dom%nz; k_srf=dom%k_srf; z_srf=dom%z_srf
    z=>dom%z; zc=>dom%zc; dz=>dom%dz; rho_=>dom%rho_avg
    allocate (f(dom%n_pout,dom%nz), q(dom%nz))
    !-----------------------------------------------

    ! get current flowrate
    allocate (array(dom%ts_pout%nv))
    call timeseries_interp(dom%ts_pout, time_day, array)
    do i=1, dom%ts_pout%nv
      dom%q_pouts(i) = (-1.)*array(i)  ! flow rate, m3/s
    enddo
    deallocate (array)

    ! water pipe source term for momentum-eq, heat-eq and SS-eq
    dom%q_pout=0.0

    do j=1, dom%n_pout

      i = dom%i_pouts(j)    ! position index of water pipe
      kb = dom%kc_bot(i)    ! bottom k of water pipe

      ! find kc
      kc = dom%k_pouts(j)

      ! normalized density gradient, 1/m
      eps = -(rho_(kc) - rho_(kc-1))/(z(kc) - z(kc-1))/rho_water
      eps = max(1.e-6, eps)

      ! flow depth of axi-symmetric jet, m
      delta = abs(dom%q_pouts(j)) / &
        (dom%fr_pouts(j)*(dom%phi_pouts(j)*pi/180)*sqrt(eps*gravity))
      delta = min(max(2*dz(kc), delta**0.333), z(nz)-z(0))

      ! velocity profile function
      f(j,:) = 0.0
      do k=kb, k_srf
        zeta = (zc(k) - dom%z_pouts(j))/delta
        if (-0.5 <= zeta .and. zeta <= 0.5) &
          f(j,k) = exp(-0.5*(zeta*3.92)**2)
      enddo

      f(j,kb:k_srf) = f(j,kb:k_srf)/sum(f(j,kb:k_srf))
      q(kb:k_srf) = dom%q_pouts(j)*f(j,kb:k_srf)

      ! source term
      dom%q_pout(i,kb:k_srf) = dom%q_pout(i,kb:k_srf) + q(kb:k_srf)
    enddo

    ! point outflow temperature and suspended solids
    do j=1, dom%n_pout
      i = dom%i_pouts(j)
      kb = dom%kc_bot(i)

      f(j,kb:k_srf) = f(j,kb:k_srf)*dom%vol(i,kb:k_srf)
      f(j,kb:k_srf) = f(j,kb:k_srf)/sum(f(j,kb:k_srf))

      dom%t_pouts(j) = sum(f(j,kb:k_srf)*dom%t(i,kb:k_srf))
      do l=1, nps
        dom%c_pouts(j,l) = sum(f(j,kb:k_srf)*dom%c(i,kb:k_srf,l))
      enddo
#ifdef SCALAR
      dom%s_pouts(j) = sum(f(j,kb:k_srf)*dom%s(i,kb:k_srf))
#endif
    enddo

    ! total flowrate
    dom%qtot_pout = sum(dom%q_pouts)

    !--- de-allocate arrays ---
    deallocate (f, q)
    !--------------------------

  end subroutine update_pointout


  subroutine update_column_flowrate(dom)
    !
    ! inputs:  u, w_srf, q_trbs, q_cnfs, q_wtps
    ! outputs: q_col
    !
    implicit none

    type(domain), intent(inout) :: dom

    integer :: nx, k_srf
    real, pointer :: q_col(:)

    integer :: i, j, k, ii
    real :: tmp, tmp1

    !--- pointer  ---------------
    nx=dom%nx; k_srf=dom%k_srf; q_col=>dom%q_col
    !----------------------------

    ! zero clear
    q_col(:) = 0.0

    ! inlet
    do i=1, nx
      do k=1, k_srf
        ii = dom%i_inlet(k)
        if (ii < i) q_col(i) = q_col(i) + dom%au(ii,k)*dom%u(ii,k)
      enddo
    enddo

    ! surface
    tmp1 = 1.e-20
    do i=1, nx
      tmp1 = tmp1 + dom%w_srf*dom%aw(i,k_srf)
    enddo
    tmp1 = dom%q_total_vol/tmp1

    tmp = 0.0
    do i=1, nx
      tmp = tmp + dom%w_srf*dom%aw(i,k_srf)*tmp1
      q_col(i) = q_col(i) - tmp
    enddo

    ! tributary
    do j=1, dom%n_trb
      do i=dom%i_trbs(j), nx
        q_col(i) = q_col(i) + dom%q_trbs(j)
      enddo
    enddo

    ! confluence
    do j=1, dom%n_cnf
      do i=dom%i_cnfs(j), nx
        q_col(i) = q_col(i) + dom%q_cnfs(j)
      enddo
    enddo

    ! water pipe
    do j=1, dom%n_wtp
      do i=dom%i_wtps(j), nx
        q_col(i) = q_col(i) + dom%q_wtps(j)
      enddo
    enddo

    ! point inflow
    do j=1, dom%n_pin
      do i=dom%i_pins(j), nx
        q_col(i) = q_col(i) + dom%q_pins(j)
      enddo
    enddo

    ! point outflow
    do j=1, dom%n_pout
      do i=dom%i_pouts(j), nx
        q_col(i) = q_col(i) + dom%q_pouts(j)
      enddo
    enddo

    ! because of numerical errors in the flowrate due to surface variation
    if (dom%n_out > 0) q_col(nx) = sum(dom%q_outs)

  end subroutine update_column_flowrate


  subroutine update_floating_fence(dom)
    !
    ! inputs:  z_srf, width_fncs
    ! outputs: k_fnc
    !
    implicit none

    type(domain), intent(inout) :: dom

    integer :: k_srf
    integer, pointer :: k_bot(:)
    real, pointer :: zc(:)

    integer :: i, j, k
    real :: z_fbot

    !--- pointer  ---------------------------------
    k_srf=dom%k_srf; k_bot=>dom%k_bot; zc=>dom%zc
    !----------------------------------------------

    do j=1, dom%n_fnc

      if (dom%type_fncs(j)) then

        ! bottom height of floating fence
        z_fbot = dom%z_srf - dom%width_fncs(j)

        ! lower bound of k_fnc
        i = dom%i_fncs(j)
        do k=k_bot(i), k_srf
          if (zc(k) < z_fbot .and. z_fbot <= zc(k+1)) then
            dom%k_fncs(1,j) = k
            exit
          endif
        enddo

        ! upper bound of k_fnc
        dom%k_fncs(2,j) = k_srf

      endif

    enddo

  end subroutine update_floating_fence


  subroutine update_boundary(dom, time_day)

    use mod_parameter, only : turbid_density

    implicit none

    type(domain), intent(inout) :: dom
    real, intent(in) :: time_day

    !-- update inflow
    call update_inflow(dom, time_day)

    !-- update outflow
    call update_outflow(dom, time_day)

    !-- update tributary
    dom%q_trb=0.0; dom%u_trb=0.0; dom%t_trb=0.0; dom%c_trb=0.0
#ifdef SCALAR
    dom%s_trb=0.0
#endif
    if (dom%n_trb > 0) call update_tributary(dom, time_day)

    !-- update confluence
    dom%q_cnf=0.0; dom%u_cnf=0.0; dom%t_cnf=0.0; dom%c_cnf=0.0
#ifdef SCALAR
    dom%s_cnf=0.0
#endif
    if (dom%n_cnf > 0) call update_confluence(dom)

    !-- update water pipe
    dom%q_wtp=0.0; dom%u_wtp=0.0; dom%t_wtp=0.0; dom%c_wtp=0.0
#ifdef SCALAR
    dom%s_wtp=0.0
#endif
    if (dom%n_wtp > 0) call update_waterpipe(dom, time_day)

    !-- update point inflow
    dom%q_pin=0.0; dom%u_pin=0.0; dom%t_pin=0.0; dom%c_pin=0.0
#ifdef SCALAR
    dom%s_pin=0.0
#endif
    if (dom%n_pin > 0) call update_pointin(dom, time_day)

    !-- update point outflow
    dom%q_pout=0.0
    if (dom%n_pout > 0) call update_pointout(dom, time_day)

  end subroutine update_boundary


  subroutine update_surface(ndom, dt_sec)

    use mod_domain, only : update_surface_layer, get_surface_height
    implicit none

    integer, intent(in) :: ndom
    real, intent(in) :: dt_sec

    integer :: i, j
    real :: z_srf, z_srf_old, total_vol
    real, allocatable :: vol_hgt(:)
    type(domain), pointer :: dom

    do i=1, ndom
      dom => doms(i)

      !-- update total volume
      dom%q_total_vol = &
        dom%qtot_in + dom%qtot_trb + dom%qtot_wtp + dom%qtot_cnf &
        + dom%qtot_pin + dom%qtot_pout - dom%qtot_out
      dom%total_vol = dom%total_vol + dt_sec*dom%q_total_vol

      !-- update column flow rate (q_col)
      call update_column_flowrate(dom)

    enddo

    do i=1, ndom
      dom => doms(i)

      z_srf_old = dom%z_srf

      ! estimate water surface elevation, z_srf
      if (dom%n_out==0) then  ! no z_srf estimation
      else if (dom%n_cnf > 0) then  ! z_srf shared by all confluent domains

        allocate (vol_hgt(0:dom%nz))

        total_vol = dom%total_vol
        vol_hgt = dom%vol_hgt

        do j=1, dom%n_cnf
          total_vol = total_vol + doms(dom%id_cnfs(j))%total_vol
          vol_hgt = vol_hgt + doms(dom%id_cnfs(j))%vol_hgt
        enddo

        call get_surface_height(dom%nz, dom%z, vol_hgt, total_vol, z_srf)

        dom%z_srf = z_srf
        dom%w_srf = (z_srf - z_srf_old)/dt_sec

        do j=1, dom%n_cnf
          doms(dom%id_cnfs(j))%z_srf = z_srf
          doms(dom%id_cnfs(j))%w_srf = (z_srf - z_srf_old)/dt_sec 
        enddo

        deallocate (vol_hgt)

      else  ! z_srf of single domain

        call get_surface_height(dom%nz, dom%z, dom%vol_hgt, &
          dom%total_vol, dom%z_srf)

        dom%w_srf = (dom%z_srf - z_srf_old)/dt_sec

      endif

    enddo

    do i=1, ndom
      dom => doms(i)

      ! update surface geometry by new z_srf
      call update_surface_layer(dom)

      !-- update floating fence
      if (dom%n_fnc > 0) call update_floating_fence(dom)

    enddo

  end subroutine update_surface


  subroutine deallocate_boundary(dom)

    use mod_domain, only : domain

    implicit none

    type(domain), intent(inout) :: dom

    nullify (dom%c_in, dom%ts_in, dom%i_inlet)

    nullify ( &
      dom%surf_outs, dom%z_outs, dom%fr_outs, dom%phi_outs, &
      dom%out_heights, dom%zKTSWs, dom%zKBSWs, &
      dom%ts_out, dom%q_outs, dom%t_outs, dom%c_outs, &
      dom%t_out, dom%c_out)

    nullify ( &
      dom%q_trbs, dom%t_trbs, dom%c_trbs, &
      dom%fr_trbs, dom%b_trbs, dom%theta_trbs, &
      dom%i_trbs, dom%ts_trbs, &
      dom%q_trb, dom%u_trb, dom%t_trb, dom%c_trb)

    nullify ( &
      dom%q_cnfs, & !dom%t_cnfs, dom%c_cnfs, &
      dom%theta_cnfs, dom%i_cnfs, &
      dom%q_cnf, dom%u_cnf, dom%t_cnf, dom%c_cnf)

    nullify ( &
      dom%id_wtps, dom%q_wtps, dom%t_wtps, dom%c_wtps, &
      dom%fr_wtps, dom%phi_wtps, dom%theta_wtps, &
      dom%i_wtps, dom%ts_wtps, &
      dom%q_wtp, dom%u_wtp, dom%t_wtp, dom%c_wtp)

    nullify ( &
      dom%q_pins, dom%t_pins, dom%c_pins, &
      dom%fr_pins, dom%phi_pins, dom%theta_pins, &
      dom%i_pins, dom%k_pins, dom%ts_pins, &
      dom%q_pin, dom%u_pin, dom%t_pin, dom%c_pin)

    nullify ( &
      dom%q_pouts, dom%t_pouts, dom%c_pouts, &
      dom%fr_pouts, dom%phi_pouts, &
      dom%i_pouts, dom%k_pouts, dom%ts_pout, dom%q_pout)

#ifdef SCALAR
    nullify ( &
      dom%s_in, dom%s_outs, dom%s_out, dom%s_trbs, dom%s_trb, &
      dom%s_cnfs, dom%s_cnf, dom%s_wtps, dom%s_wtp)
#endif

  end subroutine deallocate_boundary

end module mod_boundary


module mod_utility
  !
  ! module for utility tools
  !
  use mod_domain, only : domain

  implicit none

contains

  subroutine check_dt(dom, dt_sec, dt_cfl)

    use mod_parameter, only : nps, w_ss

    implicit none

    type(domain), intent(in) :: dom
    real, intent(in) :: dt_sec
    real, intent(out) :: dt_cfl

    integer :: i, k
    real :: uc_, wc_, dc_

    real, parameter :: eps_u=1.e-9, eps_d=1.e-9
    real :: dt_cx=1e9, dt_dx=1e9, dt_cz=1e9, dt_dz=1e9

    do i=1, dom%nx
      do k=dom%kc_bot(i), dom%k_srf
        uc_ = (dom%u(i-1,k) + dom%u(i,k))/2 + eps_u
        dc_ = dom%dmx(i,k) + eps_d
        dt_cx = min(dt_cx, dom%dx(i)/abs(uc_))
        dt_dx = min(dt_dx, dom%dx(i)**2/(2*dc_))

        wc_ = (dom%w(i,k-1) + dom%w(i,k) - w_ss(nps))/2 + eps_u
        dc_ = (dom%dmz(i,k-1) + dom%dmz(i,k))/2 + eps_d
        dt_cz = min(dt_cz, dom%dz(k)/abs(wc_))
        dt_dz = min(dt_dz, dom%dz(k)**2/(2*dc_))
      enddo
    enddo

    dt_cfl = min(dt_cx, dt_dx, dt_cz, dt_dz)

    if (dt_sec > dt_cfl) then
      print '(a,a)', 'Error: dt does not satisfy the CFL condition in domain ', &
        trim(dom%name)
      print '(a,4e10.3)', 'allowable dt limits: cx, dx, cz, dz =', &
        dt_cx, dt_dx, dt_cz, dt_dz
      stop
    endif

  end subroutine check_dt


  subroutine set_probe(dom)
    !
    ! inputs:  namelist.d01, ..
    ! outputs: n_prb, x_prb, z_prb, i_prb, k_prb
    !
    implicit none

    type(domain), intent(inout) :: dom

    integer :: i, j, k

    integer :: nx, nz
    real, pointer :: x(:), z(:)

    integer :: n_prb
    real :: x_prb(50)=-999., z_prb(50)=-999.
    namelist /probe/ n_prb, x_prb, z_prb

    !--- pointer ---
    nx=dom%nx; nz=dom%nz
    x=>dom%x; z=>dom%z
    !---------------

    ! loading from namelist file
    open(10, file=trim(dom%fname_nml), status='old')
    read(10, nml=probe)
    close(10)

    write(6, '(a)') trim(dom%name)
    write(6, nml=probe)

    allocate ( &
      dom%x_prb(n_prb), dom%z_prb(n_prb), &
      dom%i_prb(n_prb), dom%k_prb(n_prb))

    dom%n_prb = n_prb
    dom%x_prb = x_prb(1:n_prb)
    dom%z_prb = z_prb(1:n_prb)

    do j=1, n_prb

      dom%i_prb(j) = 0
      dom%k_prb(j) = 0

      do i=1, nx

        if (x(i-1) <= x_prb(j) .and. x_prb(j) <= x(i)) then
          dom%i_prb(j) = i
          exit
        endif

        if (dom%z_prb(j) < -900.) then
          dom%k_prb(j) = dom%k_srf
        else
          do k=dom%kc_bot(i), dom%k_srf
            if (z(k-1) <= z_prb(j) .and. z_prb(j) <= z(k)) then
              dom%k_prb(j) = k
              exit
            endif
          enddo
        endif

      enddo

      if (dom%i_prb(j)==0 .or. dom%k_prb(j)==0) then
        write (*,'(a)', advance='no') 'error: probing point is out of bounds, '
        write (*,'(a,a,i1)') trim(dom%name), ' point-', j
        print *, j, dom%x_prb(j), dom%z_prb(j)
        stop
      endif

    enddo

  end subroutine set_probe


  subroutine write_timeseries(dom, idom, it_snap, time_day)

    use mod_parameter, only : nps

    implicit none

    type(domain), intent(in) :: dom
    integer, intent(in) :: idom, it_snap
    real, intent(in) :: time_day

    logical, save :: first_call(10)=.true.
    character(len=128) :: fname
    character(len=200) :: header 
    character(len=50), save :: &
      fmt, fmt_in, fmt_out, fmt_sed, fmt_trb, fmt_cnf, fmt_wtp, &
      fmt_pin, fmt_pout, fmt_prb
    integer :: i, j, k, l
    real :: u_, w_

    if (first_call(idom)) then

      ! water body
      header = 'DAYS Z_SRF(m) W_SRF(m/s) TOTAL_VOL(m3)'
      write(fname, '(a,i3.3)') 'ts_waterbody.', it_snap
      fname = 'out/' // trim(dom%name) // '/' // trim(fname)
      open(100+idom, file=fname, status='unknown')
      write(100+idom, '(a)') trim(header)
      write(fmt, '(a)') '(4e15.7)'

      ! inflow
      header = 'DAYS Q(m3/s) T(deg-C)'
      do l=1, nps
        write(header(len_trim(header)+1:), '(a, "C", i1, "(g/m3)")') ' ', l
      enddo
#ifdef SCALAR
      header = trim(header) // ' S(days)'
#endif
      write(fname, '(a,i3.3)') 'ts_inflow.', it_snap
      fname = 'out/' // trim(dom%name) // '/' // trim(fname)
      open(110+idom, file=fname, status='unknown')
      write(110+idom, '(a)') trim(header)
      write(fmt_in, '(a,i0,a)') '(3e15.7, ', nps, '(e15.7))'

      ! outflow
      header = 'DAYS Q(m3/s) T(deg-C)'
      do l=1, nps
         write(header(len_trim(header)+1:), '(a, "C", i1, "(g/m3)")') ' ', l
      enddo
#ifdef SCALAR
      header = trim(header) // ' S(days)'
#endif
      write(fname, '(a,i3.3)') 'ts_outflow.', it_snap
      fname = 'out/' // trim(dom%name) // '/' // trim(fname)
      open(120+idom, file=fname, status='unknown')
      write(120+idom, '(a)') trim(header)
      write(fmt_out, '(a,i0,a)') '(3e15.7, ', nps, '(e15.7))'

      ! sedimentation
      header = 'DAYS '
      do l=1, nps
        write(header(len_trim(header)+1:), '(a, "C", i1, "(g/s)")') ' ', l
      enddo
      write(fname, '(a,i3.3)') 'ts_sediment.', it_snap
      fname = 'out/' // trim(dom%name) // '/' // trim(fname)
      open(130+idom, file=fname, status='unknown')
      write(130+idom, '(a)') trim(header)
      write(fmt_sed, '(a,i0,a)') '(3e15.7, ', nps, '(e15.7))'

      ! tributary
      if (dom%n_trb > 0) then
        header = 'DAYS Q(m3/s) T(deg-C)'
        do l=1, nps
          write(header(len_trim(header)+1:), '(a, "C", i1, "(g/m3)")') ' ', l
        enddo
#ifdef SCALAR
        header = trim(header) // ' S(days)'
#endif
        do j=1, dom%n_trb
          write(fname, '(a,i0,a,i3.3)') 'ts_tributary.', j, '.', it_snap
          fname = 'out/' // trim(dom%name) // '/' // trim(fname)
          open(200+10*idom+j, file=fname, status='unknown')
          write(200+10*idom+j, '(a)') trim(header)
        enddo
        write(fmt_trb, '(a,i0,a)') '(3e15.7, ', nps, '(e15.7))'
      endif

      ! confluence
      if (dom%n_cnf > 0) then
        header = 'DAYS Q(m3/s) T(deg-C)'
        do l=1, nps
          write(header(len_trim(header)+1:), '(a, "C", i1, "(g/m3)")') ' ', l
        enddo
#ifdef SCALAR
        header = trim(header) // ' S(days)'
#endif
        do j=1, dom%n_cnf
          write(fname, '(a,i0,a,i3.3)') 'ts_confluence.', j, '.', it_snap
          fname = 'out/' // trim(dom%name) // '/' // trim(fname)
          open(300+10*idom+j, file=fname, status='unknown')
          write(300+10*idom+j, '(a)') trim(header)
        enddo
        write(fmt_cnf, '(a,i0,a)') '(3e15.7, ', nps, '(e15.7))'
      endif

      ! water pipe
      if (dom%n_wtp > 0) then
        header = 'DAYS Q(m3/s) T(deg-C)'
        do l=1, nps
          write(header(len_trim(header)+1:), '(a, "C", i1, "(g/m3)")') ' ', l
        enddo
#ifdef SCALAR
        header = trim(header) // ' S(days)'
#endif
        do j=1, dom%n_wtp
          write(fname, '(a,i0,a,i3.3)') 'ts_waterpipe.', j, '.', it_snap
          fname = 'out/' // trim(dom%name) // '/' // trim(fname)
          open(400+10*idom+j, file=fname, status='unknown')
          write(400+10*idom+j, '(a)') trim(header)
        enddo
         write(fmt_wtp, '(a,i0,a)') '(3e15.7, ', nps, '(e15.7))'
      endif

      ! point inflow
      if (dom%n_pin > 0) then
        header = 'DAYS Q(m3/s) T(deg-C)'
        do l=1, nps
          write(header(len_trim(header)+1:), '(a, "C", i1, "(g/m3)")') ' ', l
        enddo
#ifdef SCALAR
        header = trim(header) // ' S(days)'
#endif
        do j=1, dom%n_pin
          write(fname, '(a,i0,a,i3.3)') 'ts_pointin.', j, '.', it_snap
          fname = 'out/' // trim(dom%name) // '/' // trim(fname)
          open(500+10*idom+j, file=fname, status='unknown')
          write(500+10*idom+j, '(a)') trim(header)
        enddo
         write(fmt_pin, '(a,i0,a)') '(3e15.7, ', nps, '(e15.7))'
      endif

      ! point outflow
      if (dom%n_pout > 0) then
        header = 'DAYS Q(m3/s) T(deg-C)'
        do l=1, nps
          write(header(len_trim(header)+1:), '(a, "C", i1, "(g/m3)")') ' ', l
        enddo
#ifdef SCALAR
        header = trim(header) // ' S(days)'
#endif
        do j=1, dom%n_pout
          write(fname, '(a,i0,a,i3.3)') 'ts_pointout.', j, '.', it_snap
          fname = 'out/' // trim(dom%name) // '/' // trim(fname)
          open(600+10*idom+j, file=fname, status='unknown')
          write(600+10*idom+j, '(a)') trim(header)
        enddo
         write(fmt_pout, '(a,i0,a)') '(3e15.7, ', nps, '(e15.7))'
      endif

      ! probe
      header = 'DAYS U(m/s) W(m/s) T(deg-C)'
      do l=1, nps
        write(header(len_trim(header)+1:), '(a, "C", i1, "(g/m3)")') ' ', l
      enddo
#ifdef SCALAR
      header = trim(header) // ' S(days)'
#endif
      do j=1, dom%n_prb
        write(fname, '(a,i0,a,i3.3)') 'ts_probe.', j, '.', it_snap
        fname = 'out/' // trim(dom%name) // '/' // trim(fname)
        open(1000+100*idom+j, file=fname, status='unknown')
        write(1000+100*idom+j, '(a)') trim(header)
      enddo
      write(fmt_prb, '(a,i0,a)') '(4e15.7, ', nps, '(e15.7))'

      first_call(idom) = .false.
    endif

    ! water body
    write(100+idom, fmt) time_day, dom%z_srf, dom%w_srf, dom%total_vol

    ! inflow
    write(110+idom, fmt_in, advance='no') time_day, &
      dom%q_in, dom%t_in, (dom%c_in(l), l=1, nps)
#ifdef SCALAR
    write(110+idom, '(e15.7)', advance='no') dom%s_in
#endif
    write(110+idom, *) ''

    ! outflow
    write(120+idom, fmt_out, advance='no') time_day, &
      dom%qtot_out, dom%t_out, (dom%c_out(l), l=1, nps)
#ifdef SCALAR
    write(120+idom, '(e15.7)', advance='no') dom%s_out
#endif
    write(120+idom, *) ''

    ! sedimentation
    write(130+idom, fmt_sed) time_day, (sum(dom%c_sed(:,l)), l=1, nps)

    ! tributary
    do j=1, dom%n_trb
      write(200+10*idom+j, fmt_trb, advance='no') time_day, &
        dom%q_trbs(j), dom%t_trbs(j), (dom%c_trbs(j,l), l=1, nps)
#ifdef SCALAR
      write(200+10*idom+j, '(e15.7)', advance='no') dom%s_trbs(j)
#endif
      write(200+10*idom+j, *) ''
    enddo

    ! confluence
    do j=1, dom%n_cnf
      write(300+10*idom+j, fmt_cnf, advance='no') time_day, &
        dom%q_cnfs(j), dom%t_cnfs(j), (dom%c_cnfs(j,l), l=1, nps)
#ifdef SCALAR
      write(300+10*idom+j, '(e15.7)', advance='no') dom%s_cnfs(j)
#endif
      write(300+10*idom+j, *) ''
    enddo

    ! water pipe
    do j=1, dom%n_wtp
      write(400+10*idom+j, fmt_wtp, advance='no') time_day, &
        dom%q_wtps(j), dom%t_wtps(j), (dom%c_wtps(j,l), l=1, nps)
#ifdef SCALAR
      write(400+10*idom+j, '(e15.7)', advance='no') dom%s_wtps(j)
#endif
      write(400+10*idom+j, *) ''
    enddo

    ! point inflow
    do j=1, dom%n_pin
      write(500+10*idom+j, fmt_pin, advance='no') time_day, &
        dom%q_pins(j), dom%t_pins(j), (dom%c_pins(j,l), l=1, nps)
#ifdef SCALAR
      write(500+10*idom+j, '(e15.7)', advance='no') dom%s_pins(j)
#endif
      write(500+10*idom+j, *) ''
    enddo

    ! point outflow
    do j=1, dom%n_pout
      write(600+10*idom+j, fmt_pout, advance='no') time_day, &
        dom%q_pouts(j), dom%t_pouts(j), (dom%c_pouts(j,l), l=1, nps)
#ifdef SCALAR
      write(600+10*idom+j, '(e15.7)', advance='no') dom%s_pouts(j)
#endif
      write(600+10*idom+j, *) ''
    enddo

    ! probe
    do j=1, dom%n_prb

      i = dom%i_prb(j)

      if (dom%z_prb(j) < -900.) then ! surface value
        k = dom%k_srf
      else
        k = dom%k_prb(j)
      endif

      if (k > dom%k_srf) then ! dry up
        write(1000+100*idom+j, fmt_prb, advance='no') time_day, &
          -999., -999., -999., (-999., l=1, nps)
#ifdef SCALAR
        write(1000+100*idom+j, '(e15.7)', advance='no') -999.
#endif
        write(1000+100*idom+j, *) ''
      else
        u_ = (dom%u(i-1,k) + dom%u(i,k))/2
        w_ = (dom%w(i,k-1) + dom%w(i,k))/2
        write(1000+100*idom+j, fmt_prb, advance='no') time_day, &
          u_, w_, dom%t(i,k), (dom%c(i,k,l), l=1, nps)
#ifdef SCALAR
        write(1000+100*idom+j, '(e15.7)', advance='no') dom%s(i,k)
#endif
        write(1000+100*idom+j, *) ''
      endif

    enddo

  end subroutine write_timeseries


  subroutine write_snap(dom, it_snap, time_day)
    !
    ! input  : it_snap, time_day, total_vol, z_srf, &
    !          k_srf, n_fnc, i_fncs, k_fncs, u, w, t, c
    ! output : 'out/d01/snap.<it_snap>', ..
    !
    use mod_parameter, only : nps, d_ss

    implicit none

    type(domain), intent(in) :: dom
    integer, intent(in) :: it_snap
    real, intent(in) :: time_day

    integer :: i, j, k, nx, nz
    character(len=128) :: fname, fmt

    write(fname, '(a,i4.4)') 'snap.', it_snap
    open(10, file='out/'//trim(dom%name)//'/'//trim(fname), status='unknown')

    write(10,'(3e15.7,30i5)') time_day, dom%total_vol, dom%z_srf, dom%k_srf, &
      dom%n_fnc, dom%i_fncs, ((dom%k_fncs(i,j),i=1,2),j=1,dom%n_fnc)

    write(10,'(i5,20e15.7)') nps, d_ss(1:nps)

    nx=dom%nx; nz=dom%nz
    write(fmt, '(a,i0,a)') '(', (nx+2)*(nz+2)*nps, '(e15.7))'

do i=0, nx
  do k=1, nz
    if (abs(dom%u(i,k)) < 1e-20) dom%u(i,k) = 0.0
  enddo
enddo

do i=1, nx
  do k=0, nz
    if (abs(dom%u(i,k)) < 1e-20) dom%w(i,k) = 0.0
  enddo
enddo

do i=0, nx+1
  do k=0, nz+1
    do j=1,nps
      if (abs(dom%c(i,k,j)) < 1e-20) dom%c(i,k,j) = 0.0
    enddo
  enddo
enddo

    write(10,fmt) ((dom%u(i,k),i=0,nx),k=1,nz)
    write(10,fmt) ((dom%w(i,k),i=1,nx),k=0,nz)
    write(10,fmt) ((dom%rho(i,k),i=0,nx+1),k=0,nz+1)
    write(10,fmt) ((dom%p(i,k),i=0,nx+1),k=0,nz+1)
    write(10,fmt) ((dom%t(i,k),i=0,nx+1),k=0,nz+1)
    write(10,fmt) (((dom%c(i,k,j),i=0,nx+1),k=0,nz+1),j=1,nps)
    write(10,fmt) ((dom%c_sed(i,j),i=0,nx+1),j=1,nps)
    write(10,fmt) ((dom%dhz(i,k),i=1,nx),k=0,nz)
    write(10,fmt) ((dom%dmx(i,k),i=1,nx),k=1,nz)
    write(10,fmt) ((dom%dmz(i,k),i=1,nx),k=0,nz)
#ifdef SCALAR
    write(10,fmt) ((dom%s(i,k),i=0,nx+1),k=0,nz+1)
#endif

    close(10)

  end subroutine write_snap


  subroutine read_snap(dom, it_snap, time_day)
    !
    ! input  : it_snap, "out/snap.<it_snap>", ..
    ! outout : time_day, total_vol, z_srf, k_srf, n_fnc, i_fncs, k_fncs, 
    !          u, w, t, c
    !
    use mod_parameter, only : nps, d_ss

    implicit none

    type(domain), intent(inout) :: dom
    integer, intent(in) :: it_snap
    real, intent(out) :: time_day

    character(len=128) :: fname

    write(fname, '(a,i4.4)') 'snap.', it_snap
    open(10, file='out/'//trim(dom%name)//'/'//trim(fname), status='unknown')

    read(10,*) time_day, dom%total_vol, dom%z_srf, dom%k_srf!, &
!      dom%n_fnc, dom%i_fncs, dom%k_fncs

    read(10,*) nps, d_ss(1:nps)

    read(10,*) dom%u
    read(10,*) dom%w
    read(10,*) dom%rho
    read(10,*) dom%p
    read(10,*) dom%t
    read(10,*) dom%c
    read(10,*) dom%c_sed
    read(10,*) dom%dhz
    read(10,*) dom%dmx
    read(10,*) dom%dmz
#ifdef SCALAR
    read(10,*) dom%s    
#endif

    close(10)

  end subroutine read_snap


  subroutine write_monitor(dom, time_day)
    !
    ! output monitoring data of turbulence variables
    ! keep logs only for the latest ten timesteps
    !
    use mod_parameter, only : rho_water

    implicit none

    type(domain), intent(in) :: dom
    real, intent(in) :: time_day

    integer, parameter :: nkeep = 3
    integer :: i, k, nx, nz, k_srf
    real, pointer :: dxs(:)
    real :: pg
    logical, save :: first_call = .true.
    integer, save :: step_count = 0
    character(len=128), save :: base_fname
    character(len=128) :: fname, old_fname
    logical :: ex

    if (first_call) then
      base_fname = 'out/'//trim(dom%name)//'/monitor'
      first_call = .false.
    endif

    step_count = step_count + 1
    write(fname,'(a,i6.6,a)') trim(base_fname), step_count, '.log'
    open(9000, file=fname, status='replace')
    write(9000,'(a)') 'TIME(days) I K PRESS_GRAD(Pa/m) U W NUT TKE TD_EPS'

    nx = dom%nx; nz = dom%nz; k_srf = dom%k_srf
    dxs => dom%dxs

    do i=1, nx-1
      do k=dom%kc_bot(i), k_srf
        pg = -(dom%p(i+1,k) - dom%p(i,k))/(dxs(i)*rho_water)
        write(9000,'(f10.4,2i6,6e15.7)') time_day, i, k, pg, &
          dom%u(i,k), dom%w(i,k), &
          dom%nut(i,k), dom%tke(i,k), dom%td_eps(i,k)
      enddo
    enddo

    close(9000)

    if (step_count > nkeep) then
      write(old_fname,'(a,i6.6,a)') trim(base_fname), step_count-nkeep, '.log'
      inquire(file=trim(old_fname), exist=ex)
      if (ex) then
        open(9001, file=trim(old_fname), status='old')
        close(9001, status='delete')
      endif
    endif

  end subroutine write_monitor


end module mod_utility


module mod_momentum

  implicit none

  logical :: pressure_grad=.false.  ! impose(T) or expose(F) pressure gradient
  logical :: smooth_PGX=.false.     ! impose(T) or expose(F) smoothing pressure gradient
  logical :: gravity_slope=.false.  ! impose(T) or expose(F) gravity from river slope

  real ::         &
    threshold_PGX= 10.0,   &   ! correction factor for suspended solid diffusion dcx and dcz
    beta_PGX= 0.50             ! time filter coefficient 0-1 (0: no filter, 1: almost old)

contains

  subroutine set_momentum()

    implicit none

    ! namelist
    namelist /momentum/ pressure_grad, smooth_PGX, gravity_slope, threshold_PGX, beta_PGX

    ! loading from namelist file
    open(10, file='namelist.in', status='old')
    read(10, nml=momentum)
    close(10)

    write(6, nml=momentum)

  end subroutine set_momentum

  subroutine update_momentum(dom, dt_sec)
    !
    ! inputs:  u, w, dmx, dmz, dt_sec
    ! outputs: u, w
    !
    use mod_parameter, only : rho_water, gravity
    use mod_domain, only : domain

    implicit none

    type(domain), intent(inout) :: dom
    real, intent(in) :: dt_sec

    integer :: nx, nz, k_srf
    integer, pointer, dimension(:) :: k_bot, kc_bot
    real, pointer, dimension(:) :: dx, dzs, x, z_bed
    real, pointer, dimension(:,:) :: au, aw, aus, aws, vols
    real, pointer, dimension(:,:) :: rho, p, dmx, dmz, PGX, PGXraw, PGXold

    integer :: i, j, k
    real :: rhs, fp, fm, fxm, fxp, fzm, fzp, q_col_, ini_slope, dmzmax, PG
    real, allocatable :: u(:,:), w(:,:)
    real, allocatable :: slope(:)

    !--- pointer, allocation ----------------------------
    nx=dom%nx; nz=dom%nz
    k_bot=>dom%k_bot; kc_bot=>dom%kc_bot; k_srf=dom%k_srf
    dx=>dom%dx; dzs=>dom%dzs; x=>dom%x; z_bed=>dom%z_bed
    rho=>dom%rho; p=>dom%p; dmx=>dom%dmx; dmz=>dom%dmz
    au=>dom%au; aw=>dom%aw; aus=>dom%aus; aws=>dom%aws; vols=>dom%vols
    allocate (u(0:nx,nz), w(nx,0:nz), slope(nx-1))
    do i=1, nx-1
      slope(i) = (z_bed(i) - z_bed(i+1)) / (x(i+1) - x(i))
    enddo
    !----------------------------------------------------

    u=dom%u;  w=dom%w

    do i=1, nx-1
      do k=k_bot(i), k_srf

        ! advection-x
        fm = (u(i-1,k) + u(i,k))/2*dom%aus(i,k)
        if (fm >= 0.0) then
          fm = fm*u(i-1,k)
        else
          fm = fm*u(i,k)
        endif

        fp = (u(i,k) + u(i+1,k))/2*dom%aus(i+1,k)
        if (fp >= 0.0) then
          fp = fp*u(i,k)
        else
          fp = fp*u(i+1,k)
        endif

        fxm = -fm
        fxp = -fp

        ! advection-z
        fm = (w(i,k-1) + w(i+1,k-1))/2*aws(i,k-1)
        if (k == k_bot(i)) then
          fm = fm*0.0
        else
          if (fm >= 0.0) then
            fm = fm*u(i,k-1)
          else
            fm = fm*u(i,k)
          endif
        endif

        fp = (w(i,k) + w(i+1,k))/2*aws(i,k)
        if (k == k_srf) then
          fp = fp*u(i,k)
        else
          if (fp >= 0.0) then
            fp = fp*u(i,k)
          else
            fp = fp*u(i,k+1)
          endif
        endif

        fzm = -fm
        fzp = -fp

        ! diffusion-x
        fm = dmx(i,k)*(u(i,k) - u(i-1,k))/dx(i)*aus(i,k)
        fp = dmx(i+1,k)*(u(i+1,k) - u(i,k))/dx(i)*aus(i+1,k)

        fxm = fxm + fm
        fxp = fxp + fp

        ! diffusion-z 
        if (k == k_bot(i)) then
          fm = 0.0
        else
          fm = dmz(i,k-1)*(u(i,k) - u(i,k-1))/dzs(k-1)*aws(i,k-1)
        endif

        if (k == k_srf) then
          fp = 0.0
        else
          fp = dmz(i,k)*(u(i,k+1) - u(i,k))/dzs(k)*aws(i,k)
        endif

        fzm = fzm + fm
        fzp = fzp + fp

        ! right hand side term
        rhs = fxp - fxm + fzp - fzm

        ! pressure gradient
        if (pressure_grad) then
          if (k == k_srf) then
            dom%PGX(i,k) = (p(i+1,k) - p(i,k))/rho_water/dom%dxs(i)*vols(i,k)
            if (smooth_PGX) then
              if (abs(dom%PGX(i,k)) > threshold_PGX) then
                if (i > 1 .and. k > k_bot(i-1)) then
                  dom%PGXraw(i,k) = ((p(i+1,k) - p(i-1,k))/2)/rho_water/dom%dxs(i)*vols(i,k)
                  dom%PGX(i,k) = beta_PGX*dom%PGXold(i,k) + (1.0 - beta_PGX)*dom%PGXraw(i,k)
                endif
              endif
            endif
            rhs = rhs - dom%PGX(i,k)
          else
            dom%PGX(i,k) = ( (p(i+1,k) - p(i,k)) + (p(i+1,k+1) - p(i,k+1)) )/(rho_water*2*dom%dxs(i))*vols(i,k)
            if (smooth_PGX) then
              if (abs(dom%PGX(i,k)) > threshold_PGX) then
                if (i > 1 .and. k > k_bot(i-1) .and. k+1 > k_bot(i-1)) then
                  dom%PGXraw(i,k) = ( (p(i+1,k) - p(i-1,k))/2 + (p(i+1,k+1) - p(i-1,k+1))/2 )/(rho_water*2*dom%dxs(i))*vols(i,k)
                  dom%PGX(i,k) = beta_PGX*dom%PGXold(i,k) + (1.0 - beta_PGX)*dom%PGXraw(i,k)
                endif
              endif
            endif
            rhs = rhs - dom%PGX(i,k)
          endif
          dom%PGXold(i,k) = dom%PGX(i,k)
          !print '(2i4, 3f15.5)',  i, k, p(i+1,k), p(i,k), dom%PGX(i,k)
        endif

        ! Gravity from river slope
        if (gravity_slope) then
          ini_slope = 0.001
          if (slope(i) > ini_slope) then
            rhs = rhs + gravity*sin(atan(slope(i)))*vols(i,k)
          else
            rhs = rhs + gravity*sin(atan(ini_slope))*vols(i,k)
          endif
        endif
        
        ! tributary
        if (dom%n_trb > 0) then
          if (dom%u_trb(i,k) >= 0.) &
            rhs = rhs + (dom%u_trb(i,k) - u(i,k))*dom%q_trb(i,k)
          if (dom%u_trb(i+1,k) < 0.) &
            rhs = rhs + (dom%u_trb(i+1,k) - u(i,k))*dom%q_trb(i+1,k)
        endif

        ! confluence
        if (dom%n_cnf > 0) then
          if (dom%u_cnf(i,k) >= 0.) &
            rhs = rhs + (dom%u_cnf(i,k) - u(i,k))*dom%q_cnf(i,k)
          if (dom%u_cnf(i+1,k) < 0.) &
            rhs = rhs + (dom%u_cnf(i+1,k) - u(i,k))*dom%q_cnf(i+1,k)
        endif

        ! water pipe
        if (dom%n_wtp > 0) then
          if (dom%u_wtp(i,k) >= 0.) &
            rhs = rhs + (dom%u_wtp(i,k) - u(i,k))*dom%q_wtp(i,k)
          if (dom%u_wtp(i+1,k) < 0.) &
            rhs = rhs + (dom%u_wtp(i+1,k) - u(i,k))*dom%q_wtp(i+1,k)
        endif

        ! point inflow
        if (dom%n_pin > 0) then
          if (dom%u_pin(i,k) >= 0.) &
            rhs = rhs + (dom%u_pin(i,k) - u(i,k))*dom%q_pin(i,k)
          if (dom%u_pin(i+1,k) < 0.) &
            rhs = rhs + (dom%u_pin(i+1,k) - u(i,k))*dom%q_pin(i+1,k)
        endif

        ! new u
        dom%u(i,k) = u(i,k) + dt_sec*rhs/vols(i,k)

      enddo
    enddo

    ! fence
    do j=1, dom%n_fnc
      dom%u(dom%i_fncs(j),dom%k_fncs(1,j):dom%k_fncs(2,j)) = 0.0
    enddo

    ! impose continuity on each column
    do i=1, nx-1
      q_col_ = 1e-10
      do k=k_bot(i), k_srf
        q_col_ = q_col_ + au(i,k)*dom%u(i,k)
      enddo
      do k=k_bot(i), k_srf
        dom%u(i,k) = dom%u(i,k)*abs(dom%q_col(i)/q_col_)
      enddo
    enddo

    ! new w by continuity eq.
    do i=1, nx
      dom%w(i,kc_bot(i)-1) = 0.0
      do k=kc_bot(i), k_srf
        dom%w(i,k) = ( &
          dom%q_trb(i,k) + dom%q_cnf(i,k) + dom%q_wtp(i,k) + &
          dom%q_pin(i,k) + dom%q_pout(i,k) + &
          dom%w(i,k-1)*aw(i,k-1) - &
          dom%u(i,k)*au(i,k) + dom%u(i-1,k)*au(i-1,k)) / aw(i,k)
      enddo
    enddo


    ! ghost cell
    if (k_srf < nz) then
      do i=1, nx
        dom%u(i,k_srf+1) = dom%u(i,k_srf)
        dom%w(i,k_srf+1) = dom%w(i,k_srf)
      enddo
    endif

    deallocate (u, w, slope)

  end subroutine update_momentum

end module mod_momentum


module mod_heat

  use mod_timeseries, only : timeseries, timeseries_read, timeseries_interp
  use mod_parameter, only : psat, rho_water, c_water, l_vapor, k_boltzmann

  implicit none

  logical :: &
    radiation=.true.  ! include(T) or exclude(F) solar radiation and SHX

  real ::         &
    ar=0.06,      &   ! refraction rate of solar radiation on water surface
    beta=0.5,     &   ! absorption parameter [0.4:0.6]
    eta=0.5,      &   ! decay parameter [0.3:1.5]
    z_wind=1.0,   &   ! measured height of wind speed, m
    alpha_heat=1.0    ! correction factor for dhx and dhz

  real, parameter :: &
    z0=1.e-4          ! roughness length at water surface, m

  character(len=128) :: &
    fname_meteo=''    ! meteorological data file

  type(timeseries) :: &
    ts_meteo          ! time series of met data

  real ::  &
    time,  &          ! days
    solar, &          ! solar radiation, kcal/m2/day
    cloud, &          ! cloud coverage, [0:1]
    t_air, &          ! air temperature, deg-C
    rh,    &          ! relative humidity, [0:1]
    wind,  &          ! wind speed at 0.15m, m/s
    rain              ! precipitation, mm/hour

contains

  subroutine set_heat()

    implicit none

    ! namelist (heat)
    namelist /heat/ radiation, ar, beta, eta, z_wind, alpha_heat, fname_meteo

    ! loading from namelist file
    open(10, file='namelist.in', status='old')
    read(10, nml=heat)
    close(10)

    write(6, nml=heat)

    ! setup meteorological data
    call timeseries_read(fname_meteo, ts_meteo)

  end subroutine set_heat


  subroutine update_meteorology(time_day)

    implicit none
    real, intent(in) :: time_day
    real, allocatable :: array(:)

    allocate (array(ts_meteo%nv))

    call timeseries_interp(ts_meteo, time_day, array)

    solar = array(1)    ! solar radiation, kcal/m2/day
    cloud = array(2)    ! cloud coverage, [0:1]
    t_air = array(3)    ! air temperature, deg-C
    rh    = array(4)    ! relative humidity, [0:1]
    wind  = array(5)    ! wind speed at 0.15m, m/s
    rain  = array(6)    ! precipitation, mm/hour

    deallocate (array)

  end subroutine update_meteorology


  real function surface_heat_flux(t_srf)
    ! water surface heat echange, Km/s
    implicit none

    real, intent(in) :: t_srf ! water surface temperatrure, deg-C
    real :: phi_0, phi_s, phi_ec, phi_ra, tmp, ts, ta
    real :: ce0, cc0, ce, cc, s, w2

    ! evaporation and conduction (Rohwer's eq.)
    tmp = psat(t_srf) - rh*psat(t_air)

    if (.true.) then

      ! Rohwer's eq.
      wind = wind * log(0.15/z0)/log(z_wind/z0)  ! wind speed at 0.15m

      phi_ec = (3.08e-4 + 1.85e-4*wind*0.635)*rho_water &
        * (tmp*(l_vapor + (c_water-0.54)*t_srf) + 269.1*(t_srf - t_air))

    else

      ! Kondo's eq.
      wind = wind * log(10./z0)/log(z_wind/z0)  ! wind speed at 10m
      w2 = wind**2

      ce0 = 1.49
      cc0 = 1.43
      if (0.3 <= wind .and. wind < 2.3) then
        ce0 = 1.230*wind**(-0.160)
        cc0 = 1.185*wind**(-0.157)
      else if (2.3 <= wind .and. wind < 5.0) then
        ce0 = 0.969 + 0.0521*wind
        cc0 = 0.927 + 0.0546*wind
      else if (5.0 <= wind .and. wind < 8.0) then
        ce0 = 1.18 + 0.01*wind
        cc0 = 1.15 + 0.01*wind
      else if (8.0 <= wind) then
        ce0 = 1.170 + 0.0144*wind - 0.004*w2
        cc0 = 1.141 + 0.0147*wind - 0.00045*w2
      endif

      if (abs(t_srf - t_air) < 1.e-3) then
        ce = ce0
        cc = cc0
      else
        ! stability parameter, s
        if (wind > 0.0) then
          s = (t_srf - t_air)/w2/(1.0 + 0.01+w2/abs(t_srf - t_air))
        else
          s = 100*(t_srf - t_air)
        endif
        if (s < 0.) then ! unstable
          if (-3.3 < s) then
            ce = ce0*(0.1 + 0.03*s + 0.9*exp(4.8*s))
            cc = cc0*(0.1 + 0.03*s + 0.9*exp(4.8*s))
          else
            ce = 0.0
            cc = 0.0
          endif
        else
          ce = ce0*(1.0 + 0.63+sqrt(s))
          cc = cc0*(1.0 + 0.63+sqrt(s))
        endif
      endif
      
      phi_ec = 50.33*ce*tmp*wind + 24.88*cc*(t_srf - t_air)*wind

    endif

    if (tmp < 0.0) phi_ec = 0.0

    ! atmospheric radiation (Swinbank's eq.)
    ts = t_srf + 273.15
    ta = t_air + 273.15
    phi_ra = 0.97*k_boltzmann*(ts**4 - 0.937e-5*ta**6*(1.+0.17*cloud**2))

    ! solar radiation, kcal/m2/day
    phi_0 = (1.0 - ar)*solar

    ! water surface heat exchange, kcal/m2/day
    phi_s = phi_0 - (phi_ec + phi_ra) ! kcc code

    ! kcal/m2/day --> Km/s
    surface_heat_flux = phi_s/(rho_water*c_water)/(3600*24)

  end function surface_heat_flux


  real function radiation_heat_flux(depth)
    ! radiation heat flux in water body, Km/s
    implicit none

    real, intent(in) :: depth ! water depth, m
    real :: phi_0, phi

    ! penetration of solar radiation, kcal/m2/day
    phi_0 = (1.0 - ar)*solar

    ! heat flux in water
    phi = (1.0 - beta)*phi_0*exp(-eta*depth)

    ! kcal/m2/day --> Km/s
    radiation_heat_flux = phi/(rho_water*c_water)/(3600*24)

  end function radiation_heat_flux


  subroutine update_heat(dom, dt_sec)
    !
    ! inputs:  t, u, w, dhx, dhz, dt_sec
    ! outputs: t
    !
    use mod_domain, only : domain
    
    implicit none

    type(domain), intent(inout) :: dom
    real, intent(in) :: dt_sec

    integer :: nx, nz, k_srf
    integer, pointer, dimension(:) :: k_bot, kc_bot
    real, pointer, dimension(:) :: dxs, dzs
    real, pointer, dimension(:,:) :: u, w, dhx, dhz
    real, pointer, dimension(:,:) :: au, aw, aus, aws, vol

    integer :: i, j, k
    real :: fm, fp, fxm, fxp, fzm, fzp, rhs
    real, allocatable :: t(:,:)

    !--- pointer, allocation ----------------------------
    nx=dom%nx; nz=dom%nz
    k_bot=>dom%k_bot; kc_bot=>dom%kc_bot; k_srf=dom%k_srf
    dxs=>dom%dxs; dzs=>dom%dzs
    u=>dom%u; w=>dom%w; dhx=>dom%dhx; dhz=>dom%dhz
    au=>dom%au; aw=>dom%aw; aus=>dom%aus; aws=>dom%aws; vol=>dom%vol
    allocate (t(0:nx+1,0:nz+1))
    !----------------------------------------------------

    t=dom%t

    do i=1, nx
      do k=kc_bot(i), k_srf

        ! advection-x
        fm = u(i-1,k)*au(i-1,k)
        if (fm > 0.0) then
          fm = fm*t(i-1,k)
        else
          fm = fm*t(i,k)
        endif

        fp = u(i,k)*au(i,k)
        if (fp > 0.0) then
          fp = fp*t(i,k)
        else
          fp = fp*t(i+1,k)
        endif

        fxm = -fm
        fxp = -fp

        ! advection-z
        fm = w(i,k-1)*aw(i,k-1)
        if (k == kc_bot(i)) then
          fm = fm*0.0
        else
          if (fm > 0.0) then
            fm = fm*t(i,k-1)
          else
            fm = fm*t(i,k)
          endif
        endif

        fp = w(i,k)*aw(i,k)
        if (k == k_srf) then
          fp = fp*t(i,k)
        else
          if (fp > 0.0) then
            fp = fp*t(i,k)
          else
            fp = fp*t(i,k+1)
          endif
        endif

        fzm = -fm
        fzp = -fp

        ! diffusion-x
        if (i == 1 .or. k <= k_bot(i-1)) then
          fm = 0.0
        else
          fm = alpha_heat*dhx(i-1,k)*(t(i,k) - t(i-1,k))/dxs(i-1)*au(i-1,k)
        endif

        if (i == nx .or. k <= k_bot(i)) then
          fp = 0.0
        else
          fp = alpha_heat*dhx(i,k)*(t(i+1,k) - t(i,k))/dxs(i)*au(i,k)
        endif

        fxm = fxm + fm
        fxp = fxp + fp

        ! diffusion-z
        if (k == kc_bot(i)) then
          fm = fm + 0.0
        else
          fm = fm + alpha_heat*dhz(i,k-1)*(t(i,k) - t(i,k-1))/dzs(k-1)*aw(i,k-1)
        endif

        if (k == k_srf) then
          fp = fp + 0.0
        else
          fp = fp + alpha_heat*dhz(i,k)*(t(i,k+1) - t(i,k))/dzs(k)*aw(i,k)
        endif

        ! radiation heat flux
        if (radiation) then
          fm = fm + radiation_heat_flux(dom%z_srf - dom%z(k-1))*aw(i,k-1)
          if (k == k_srf) then
            fp = fp + surface_heat_flux(t(i,k_srf))*aw(i,k)
          else
            fp = fp + radiation_heat_flux(dom%z_srf - dom%z(k))*aw(i,k)
          endif
        endif

        fzm = fzm + fm
        fzp = fzp + fp

        ! fence
        do j=1, dom%n_fnc
          if (dom%k_fncs(1,j) <= k .and. k <= dom%k_fncs(2,j)) then
            if (i-1 == dom%i_fncs(j)) fxm = 0.0
            if (i == dom%i_fncs(j)) fxp = 0.0
          endif
        enddo

        ! right hand side term
        rhs = fxp - fxm + fzp - fzm

        ! tributary
        if (dom%n_trb > 0) rhs = rhs + (dom%t_trb(i,k) - t(i,k))*dom%q_trb(i,k)

        ! confluence
        if (dom%n_cnf > 0) rhs = rhs + (dom%t_cnf(i,k) - t(i,k))*dom%q_cnf(i,k)

        ! water pipe
        if (dom%n_wtp > 0) rhs = rhs + (dom%t_wtp(i,k) - t(i,k))*dom%q_wtp(i,k)

        ! point inflow
        if (dom%n_pin > 0) rhs = rhs + (dom%t_pin(i,k) - t(i,k))*dom%q_pin(i,k)

        ! divergence correction
        rhs = rhs + (u(i,k)*au(i,k) - u(i-1,k)*au(i-1,k) &
          + w(i,k)*aw(i,k) - w(i,k-1)*aw(i,k-1))*t(i,k)

        ! new t
        dom%t(i,k) = t(i,k) + dt_sec*rhs/vol(i,k)

      enddo
    enddo

    ! ghost cell
    if (k_srf < nz) then
      do i=1, nx
        dom%t(i,k_srf+1) = dom%t(i,k_srf)
      enddo
    endif

    deallocate (t)

  end subroutine update_heat

end module mod_heat


module mod_suspended_solids

  implicit none

  real ::         &
    alpha_ss=1.0    ! correction factor for suspended solid diffusion dcx and dcz

contains

  subroutine set_suspended_solids()

    use mod_parameter, only : gravity, rho_water, mu_water, &
      nps, rho_ss, d_ss, w_ss

    implicit none

    ! namelist (particle)
    real :: particle_size(20)=0.0   ! particle size for each bins, m
    namelist /particle/ nps, particle_size, rho_ss, alpha_ss

    integer :: l
    real :: rep

    ! loading from namelist file
    open(10, file='namelist.in', status='old')
    read(10, nml=particle)
    close(10)

    write(6, nml=particle)

    ! particle size, settling velocity(Stokes' law) 
    ! and particle Reynolds number
    allocate (d_ss(nps), w_ss(nps))

    do l=1, nps

      d_ss(l) = particle_size(l)
      w_ss(l) = gravity*(rho_ss - rho_water)*d_ss(l)**2/(18*mu_water)

      rep = w_ss(l)*d_ss(l)*rho_water/mu_water ! particle Re
      if (rep > 1.5) then
        print *, 'Error: particle Reynolds number is too large, Rep=', rep
        stop
      endif

    enddo

  end subroutine set_suspended_solids


  subroutine update_suspended_solids(dom, dt_sec)
    !
    ! inputs:  c, u, w, dcx, dcz, dt_sec
    ! outputs: c, c_sed
    !
    use mod_parameter, only : nps, w_ss
    use mod_domain, only : domain
    
    implicit none

    type(domain), intent(inout) :: dom
    real, intent(in) :: dt_sec

    integer :: nx, nz, k_srf
    integer, pointer, dimension(:) :: k_bot, kc_bot
    real, pointer, dimension(:) :: dxs, dzs
    real, pointer, dimension(:,:) :: u, w, dcx, dcz
    real, pointer, dimension(:,:) :: au, aw, aus, aws, vol

    integer :: i, j, k, l
    real :: fm, fp, fxm, fxp, fzm, fzp, rhs
    real, allocatable :: c(:,:,:)

    !--- pointer, allocation -----------------------------
    nx=dom%nx; nz=dom%nz
    k_bot=>dom%k_bot; kc_bot=>dom%kc_bot; k_srf=dom%k_srf
    dxs=>dom%dxs; dzs=>dom%dzs; dcx=>dom%dcx; dcz=>dom%dcz
    u=>dom%u; w=>dom%w
    au=>dom%au; aw=>dom%aw; aus=>dom%aus; aws=>dom%aws; vol=>dom%vol
    allocate (c(0:nx+1,0:nz+1,nps))
    !----------------------------------------------------

    c = dom%c
    dom%c_sed = 0.0  ! deposition of suspended solids, g/s

    do l=1, nps  ! particle-size loop

    do i=1, nx
      do k=kc_bot(i), k_srf

        ! advection-x
        fm = u(i-1,k)*au(i-1,k)
        if (fm > 0.0) then
          fm = fm*c(i-1,k,l)
        else
          fm = fm*c(i,k,l)
        endif

        fp = u(i,k)*au(i,k)
        if (fp > 0.0) then
          fp = fp*c(i,k,l)
        else
          fp = fp*c(i+1,k,l)
        endif

        fxm = -fm
        fxp = -fp

        ! advection-z
        fm = (w(i,k-1) - w_ss(l))*aw(i,k-1)
        if (k == kc_bot(i)) then
          fm = fm*c(i,k,l)
          dom%c_sed(i,l) = dom%c_sed(i,l) - fm ! sedimentation rate, g/s
        else
          if (fm > 0.0) then
            fm = fm*c(i,k-1,l)
          else
            fm = fm*c(i,k,l)
          endif
        endif

        if (k == k_srf) then
          fp = w(i,k)*aw(i,k)*c(i,k,l)
        else
          fp = (w(i,k) - w_ss(l))*aw(i,k)
          if (fp > 0.0) then
            fp = fp*c(i,k,l)
          else
            fp = fp*c(i,k+1,l)
          endif
        endif

        fzm = -fm
        fzp = -fp

        ! diffusion-x
        if (i == 1 .or. k <= k_bot(i-1)) then
          fm = 0.0
        else
          fm = alpha_ss*dcx(i-1,k)*(c(i,k,l) - c(i-1,k,l))/dxs(i-1)*au(i-1,k)
        endif

        if (i == nx .or. k <= k_bot(i)) then
          fp = 0.0
        else
          fp = alpha_ss*dcx(i,k)*(c(i+1,k,l) - c(i,k,l))/dxs(i)*au(i,k)
        endif

        fxm = fxm + fm
        fxp = fxp + fp

        ! diffusion-z
        if (k == kc_bot(i)) then
          fm = fm + 0.0
        else
          fm = fm + alpha_ss*dcz(i,k-1)*(c(i,k,l) - c(i,k-1,l))/dzs(k-1)*aw(i,k-1)
        endif

        if (k == k_srf) then
          fp = fp + 0.0
        else
          fp = fp + alpha_ss*dcz(i,k)*(c(i,k+1,l) - c(i,k,l))/dzs(k)*aw(i,k)
        endif

        fzm = fzm + fm
        fzp = fzp + fp

        ! fence
        do j=1, dom%n_fnc
          if (dom%k_fncs(1,j) <= k .and. k <= dom%k_fncs(2,j)) then
            if (i-1 == dom%i_fncs(j)) fxm = 0.0
            if (i == dom%i_fncs(j)) fxp = 0.0
          endif
        enddo

        ! right hand side term
        rhs = fxp - fxm + fzp - fzm

        ! tributary
        if (dom%q_trb(i,k) > 0.) &
          rhs = rhs + (dom%c_trb(i,k,l) - c(i,k,l))*dom%q_trb(i,k)

        ! confluence
        if (dom%q_cnf(i,k) > 0.) &
          rhs = rhs + (dom%c_cnf(i,k,l) - c(i,k,l))*dom%q_cnf(i,k)

        ! water pipe
        if (dom%q_wtp(i,k) > 0.) &
          rhs = rhs + (dom%c_wtp(i,k,l) - c(i,k,l))*dom%q_wtp(i,k)

        ! point inflow
        if (dom%q_pin(i,k) > 0.) &
          rhs = rhs + (dom%c_pin(i,k,l) - c(i,k,l))*dom%q_pin(i,k)

        ! divergence correction
        rhs = rhs + (u(i,k)*au(i,k) - u(i-1,k)*au(i-1,k) &
          + w(i,k)*aw(i,k) - w(i,k-1)*aw(i,k-1))*c(i,k,l)

        ! new c
        dom%c(i,k,l) = c(i,k,l) + dt_sec*rhs/vol(i,k)
      enddo
    enddo

    ! ghost cell
    if (k_srf < nz) then
      do i=1, nx
        dom%c(i,k_srf+1,l) = dom%c(i,k_srf,l)
      enddo
    endif

    enddo  ! particle-size loop

    deallocate (c)

  end subroutine update_suspended_solids

end module mod_suspended_solids


#ifdef SCALAR

module mod_scalar

  implicit none

contains

  subroutine update_scalar(dom, dt_sec)
    !
    ! inputs:  s, u, w, dcx, dcz, dt_sec
    ! outputs: s
    !
    use mod_domain, only : domain
    
    implicit none

    type(domain), intent(inout) :: dom
    real, intent(in) :: dt_sec

    integer :: nx, nz, k_srf
    integer, pointer, dimension(:) :: k_bot, kc_bot
    real, pointer, dimension(:) :: dxs, dzs
    real, pointer, dimension(:,:) :: u, w, dcx, dcz
    real, pointer, dimension(:,:) :: au, aw, aus, aws, vol

    integer :: i, j, k
    real :: fm, fp, fxm, fxp, fzm, fzp, rhs
    real, allocatable :: s(:,:)

    !--- pointer, allocation -----------------------------
    nx=dom%nx; nz=dom%nz
    k_bot=>dom%k_bot; kc_bot=>dom%kc_bot; k_srf=dom%k_srf
    dxs=>dom%dxs; dzs=>dom%dzs
    u=>dom%u; w=>dom%w; dcx=>dom%dcx; dcz=>dom%dcz
    au=>dom%au; aw=>dom%aw; aus=>dom%aus; aws=>dom%aws; vol=>dom%vol
    allocate (s(0:nx+1,0:nz+1))
    !----------------------------------------------------

    s = dom%s

    do i=1, nx
      do k=kc_bot(i), k_srf

        ! advection-x
        fm = u(i-1,k)*au(i-1,k)
        if (fm > 0.0) then
          fm = fm*s(i-1,k)
        else
          fm = fm*s(i,k)
        endif

        fp = u(i,k)*au(i,k)
        if (fp > 0.0) then
          fp = fp*s(i,k)
        else
          fp = fp*s(i+1,k)
        endif

        fxm = -fm
        fxp = -fp

        ! advection-z
        fm = w(i,k-1)*aw(i,k-1)
        if (k == kc_bot(i)) then
          fm = fm*0.0
        else
          if (fm > 0.0) then
            fm = fm*s(i,k-1)
          else
            fm = fm*s(i,k)
          endif
        endif

        fp = w(i,k)*aw(i,k)
        if (k == k_srf) then
          fp = fp*s(i,k)
        else
          if (fp > 0.0) then
            fp = fp*s(i,k)
          else
            fp = fp*s(i,k+1)
          endif
        endif

        fzm = -fm
        fzp = -fp

        ! diffusion-x
        if (i == 1 .or. k <= k_bot(i-1)) then
          fm = 0.0
        else
          fm = dcx(i-1,k)*(s(i,k) - s(i-1,k))/dxs(i-1)*au(i-1,k)
        endif

        if (i == nx .or. k <= k_bot(i)) then
          fp = 0.0
        else
          fp = dcx(i,k)*(s(i+1,k) - s(i,k))/dxs(i)*au(i,k)
        endif

        fxm = fxm + fm
        fxp = fxp + fp

        ! diffusion-z
        if (k == kc_bot(i)) then
          fm = fm + 0.0
        else
          fm = fm + dcz(i,k-1)*(s(i,k) - s(i,k-1))/dzs(k-1)*aw(i,k-1)
        endif

        if (k == k_srf) then
          fp = fp + 0.0
        else
          fp = fp + dcz(i,k)*(s(i,k+1) - s(i,k))/dzs(k)*aw(i,k)
        endif

        fzm = fzm + fm
        fzp = fzp + fp

        ! fence
        do j=1, dom%n_fnc
          if (dom%k_fncs(1,j) <= k .and. k <= dom%k_fncs(2,j)) then
            if (i-1 == dom%i_fncs(j)) fxm = 0.0
            if (i == dom%i_fncs(j)) fxp = 0.0
          endif
        enddo

        ! right hand side term
        rhs = fxp - fxm + fzp - fzm

        ! tributary
        if (dom%n_trb > 0) rhs = rhs + (dom%s_trb(i,k) - s(i,k))*dom%q_trb(i,k)

        ! confluence
        if (dom%n_cnf > 0) rhs = rhs + (dom%s_cnf(i,k) - s(i,k))*dom%q_cnf(i,k)

        ! water pipe
        if (dom%n_wtp > 0) rhs = rhs + (dom%s_wtp(i,k) - s(i,k))*dom%q_wtp(i,k)

        ! point inflow
        if (dom%n_pin > 0) rhs = rhs + (dom%s_pin(i,k) - s(i,k))*dom%q_pin(i,k)

        ! divergence correction
        rhs = rhs + (u(i,k)*au(i,k) - u(i-1,k)*au(i-1,k) &
          + w(i,k)*aw(i,k) - w(i,k-1)*aw(i,k-1))*s(i,k)

        ! time generator for age estimation, m3/days
        if (.true.) rhs = rhs + vol(i,k)/(3600*24)

        ! new s
        dom%s(i,k) = s(i,k) + dt_sec*rhs/vol(i,k)

      enddo
    enddo

    ! ghost cell
    if (k_srf < nz) then
      do i=1, nx
        dom%s(i,k_srf+1) = dom%s(i,k_srf)
      enddo
    endif

    deallocate (s)

  end subroutine update_scalar

end module mod_scalar

#endif


program main

  use mod_domain, only : domain, &
    doms, &
    set_geometry, allocate_variables, deallocate_variables, &
    set_initial, update_surface_layer, update_density, &
    set_turbulence, update_turbulence, write_geo, write_exceed_point

  use mod_boundary, only : &
    set_boundary, update_boundary, update_surface, deallocate_boundary

  use mod_heat, only : &
    set_heat, update_meteorology, radiation, update_heat

  use mod_momentum, only : &
    set_momentum, update_momentum

  use mod_suspended_solids, only : &
    set_suspended_solids, update_suspended_solids

  use mod_utility, only : &
    check_dt, set_probe, write_timeseries, read_snap, write_snap !, write_monitor

#ifdef SCALAR
  use mod_scalar, only : update_scalar
#endif

  implicit none

  type(domain), pointer :: dom

  ! namelist (control)
  integer :: ndom             ! number of domains
  character(len=128) :: names(10)='' ! name for each domains
  logical :: restart=.false.  ! restart(T) or cold start(F)
  integer :: it_rst=10        ! restart file number
  logical :: fixed_dt=.true.  ! fixed(T) or Variable(F) time step
  real    :: dt_sec=10.       ! initial dt, s
  integer :: it_max=1000      ! max time step
  integer :: it_out_fld=100   ! output interval of 2d field 
  integer :: it_out_ts=10     ! output interval of timeseries

  namelist /control/ ndom, names, restart, it_rst, &
    fixed_dt, dt_sec, it_max, it_out_fld, it_out_ts

  ! misc
  integer :: i, it, it_snap
  real :: time_day, dt_cfl

  ! loading from namelist file
  open(10, file='namelist.in', status='old')
  read(10, nml=control)
  close(10)

  write(6, nml=control)

  ! allocate domains
  allocate (doms(ndom))

  ! read namelist params
  do i=1, ndom
    doms(i)%id = i
    doms(i)%name = trim(names(i))
    doms(i)%fname_nml = trim('namelist' // '.' // doms(i)%name)
  enddo

  ! setupd geometry data (grid, bed, width, volume, ..)
  do i=1, ndom
    call set_geometry(doms(i))
  enddo

  ! setup momentum, heat and suspended solids eq.
  call set_momentum()
  call set_heat()
  call set_suspended_solids()

  ! allocate field variables
  do i=1, ndom
    call allocate_variables(doms(i))
  enddo

  ! initialization
  if (restart) then
    it_snap = it_rst
    do i=1, ndom
      call read_snap(doms(i), it_snap, time_day)
    enddo
  else
    time_day = 0.0
    it_snap = 0
    do i=1, ndom
      call set_initial(doms(i))
    enddo
  endif

  do i=1, ndom
    call set_boundary(doms(i))
  enddo

  do i=1, ndom
    dom => doms(i)
    call update_density(dom)
    call update_boundary(dom, time_day)
    call set_turbulence(dom)
    call update_turbulence(dom, dt_sec)
    call write_exceed_point(dom, 0, time_day)
    call set_probe(doms(i))
  enddo

  call update_surface(ndom, dt_sec)

  do i=1, ndom
    dom => doms(i)
    call write_geo(dom)
    if (.not. restart) call write_snap(dom, it_snap, time_day)
    call write_timeseries(dom, i, it_snap, time_day)
  enddo


  ! std log
  write(6,'(a)')
  write(6,'(a)',advance='no') 'days, dom1:z_srf z_in del_in maxT maxC rho_in i_inlet k_srf'
  if (ndom > 1) then
    write(6,'(a)') ', dom2:z_srf z_in del_in maxT maxC rho_in i_inlet k_srf'
  else
    write(6,*) ''
  endif

  write(6,'(f6.1,a,5f6.1,4f7.1,2i4)', advance='no') time_day, ', ', &
    doms(1)%z_srf, doms(1)%z_in, doms(1)%delta_in*2, &
    maxval(doms(1)%t), maxval(sum(doms(1)%c,dim=3)), &
    doms(1)%rho_in, doms(1)%log_i_inlet, doms(1)%k_srf
  if (ndom > 1) then
    write(6,'(a,5f5.1)') ', ', &
      doms(2)%z_srf, doms(2)%z_in, doms(2)%delta_in*2, &
      maxval(doms(2)%t), maxval(sum(doms(2)%c,dim=3))
  else
    write(6,*) ''
  endif

  ! time loop
  do it=1, it_max

    do i=1, ndom
      call check_dt(doms(i), dt_sec, dt_cfl)
    enddo

    time_day = time_day + dt_sec/(3600*24)

    ! update meteorological data
    if (radiation) call update_meteorology(time_day)

    ! boundary conditions
    do i=1, ndom
      call update_boundary(doms(i), time_day)
    enddo

    ! update water surface
    call update_surface(ndom, dt_sec)

    ! primary variables
    do i=1, ndom
      dom => doms(i)
      call update_momentum(dom, dt_sec)
      call update_heat(dom, dt_sec)
      call update_suspended_solids(dom, dt_sec)
#ifdef SCALAR
      call update_scalar(dom, dt_sec)
#endif
    enddo

    ! model variables
    do i=1, ndom
      call update_density(doms(i))
      call update_turbulence(doms(i), dt_sec)
      call write_exceed_point(doms(i), it, time_day)
    !  call write_monitor(doms(i), time_day)
    enddo

    ! output 2d field
    if (mod(it, it_out_fld) == 0) then
      it_snap = it_snap + 1
      do i=1, ndom
        call write_snap(doms(i), it_snap, time_day)
      enddo

      ! std log
      write(6,'(f6.1,a,5f6.1,4f7.1,2i4)', advance='no') time_day, ', ', &
        doms(1)%z_srf, doms(1)%z_in, doms(1)%delta_in*2, &
        maxval(doms(1)%t), maxval(sum(doms(1)%c,dim=3)), &
        doms(1)%rho_in, doms(1)%log_i_inlet, doms(1)%k_srf
      if (ndom > 1) then
        write(6,'(a,5f5.1)') ', ', &
          doms(2)%z_srf, doms(2)%z_in, doms(2)%delta_in*2, &
          maxval(doms(2)%t), maxval(sum(doms(2)%c,dim=3))
      else
        write(6,*) ''
      endif
    endif

    ! output timeseries
    if (mod(it, it_out_ts) == 0) then
      do i=1, ndom
        call write_timeseries(doms(i), i, it_snap, time_day)
      enddo
    endif

  enddo

  ! deallocate arrays
  do i=1, ndom
    call deallocate_boundary(doms(i))
    call deallocate_variables(doms(i))
  enddo

  deallocate (doms)

  print *, 'finish normaly'

end program main
