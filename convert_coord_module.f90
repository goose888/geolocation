module convert_coord_module
  implicit none

! ==========================================================================
!  ezlhconv.for - FORTRAN routines for conversion of azimuthal 
! 		equal area and equal area cylindrical grid coordinates
! 
! 	30-Jan.-1992 H.Maybee
! 	20-Mar-1992 Ken Knowles  303-492-0644  knowles@kryos.colorado.edu
!        16-Dec-1993 MJ Brodzik   303-492-8263  brodzik@jokull.colorado.edu
!                    Copied from nsmconv.f, changed resolutions from 
!                    40-20-10 km to 25-12.5 km
!        21-Dec-1993 MJ Brodzik   303-492-8263  brodzik@jokull.colorado.edu
!                    Fixed sign of Southern latitudes in ease_inverse.
! 	12-Sep-1994 David Hoogstrate 303-492-4116 hoogstra@jokull.colorado.edu
! 		    Changed grid cell size. Changed "c","f" to "l","h"
! 	25-Oct-1994 David Hoogstrate 303-492-4116 hoogstra@jokull.colorado.edu
! 		    Changed row size from 587 to 586 for Mercator projection
! 		    Changed function names to "ezlh-.."
!       Shijie: F90 version
! $Log: ezlhconv.f,v $
! Revision 1.3  1994/11/01 23:40:43  brodzik
! Replaced all references to 'ease' with 'ezlh'
! 
! ==========================================================================

  public :: ezlh_convert
  public :: ezlh_inverse
  private :: deg
  private :: rad

  ! radius of the earth (km), authalic sphere based on International datum 
  real, parameter :: RE_km = 6371.228
  ! nominal cell size in kilometers
  real, parameter :: CELL_km = 25.067525
  ! scale factor for standard paralles at +/-30.00 degrees
  real, parameter :: COS_PHI1 = .866025403
  real, parameter :: PI = 3.141592653589793

contains

!--------------------------------------------------------------------------
  subroutine ezlh_convert (grid, lat, lon, r, s)
    character(len=2), intent(in) :: grid
    real, intent(in) :: lat, lon
    real, intent(out) :: r, s
!
!	convert geographic coordinates (spherical earth) to 
!	azimuthal equal area or equal area cylindrical grid coordinates
!
!	status = ezlh_convert (grid, lat, lon, r, s)
!
!	input : grid - projection name '[NSM][lh]'
!               where l = "low"  = 25km resolution
!                     h = "high" = 12.5km resolution
!		lat, lon - geo. coords. (decimal degrees)
!
!	output: r, s - column, row coordinates
!
!	result: status = 0 indicates normal successful completion
!			-1 indicates error status (point not on grid)
!
!--------------------------------------------------------------------------
    ! Local Vars
    integer :: cols, rows, scale
    real :: Rg, phi, lam, rho
    real :: r0, s0, c


    if ((grid(1:1).eq.'N').or.(grid(1:1).eq.'S')) then
      cols = 721
      rows = 721
    else if (grid(1:1).eq.'M') then
      cols = 1383
      rows = 586
    else
      print *, 'ezlh_convert: unknown projection: ', grid
      return
    endif
    
    if (grid(2:2).eq.'l') then
      scale = 1
    else if (grid(2:2).eq.'h') then
      scale = 2
    else
      print *, 'ezlh_convert: unknown projection: ', grid
      return
    endif
    
    Rg = scale * RE_km/CELL_km

    !
    ! r0,s0 are defined such that cells at all scales
    ! have coincident center points
    !
    r0 = (cols-1)/2. * scale
    s0 = (rows-1)/2. * scale

    phi = rad(lat)
    lam = rad(lon)
    
    if (grid(1:1).eq.'N') then
      rho = 2 * Rg * sin(PI/4. - phi/2.)
      r = r0 + rho * sin(lam)
      s = s0 + rho * cos(lam)
    
    else if (grid(1:1).eq.'S') then
      rho = 2 * Rg * cos(PI/4. - phi/2.)
      r = r0 + rho * sin(lam)
      s = s0 - rho * cos(lam)
    
    else if (grid(1:1).eq.'M') then
      r = r0 + Rg * lam * COS_PHI1
      s = s0 - Rg * sin(phi) / COS_PHI1
    
    endif
  
  end subroutine ezlh_convert

!--------------------------------------------------------------------------
  subroutine ezlh_inverse (grid, r, s, lat, lon)
    character(len=2) :: grid
    real, intent(in) :: r, s
    real, intent(out) :: lat, lon
!
!	convert azimuthal equal area or equal area cylindrical 
!	grid coordinates to geographic coordinates (spherical earth)
!
!	status = ezlh_inverse (grid, r, s, lat, lon)
!
!	input : grid - projection name '[NSM][lh]'
!               where l = "low"  = 25km resolution
!                     h = "high" = 12.5km resolution
!		r, s - grid column and row coordinates
!
!	output: lat, lon - geo. coords. (decimal degrees)
!
!	result: status = 0 indicates normal successful completion
!			-1 indicates error status (point not on grid)
!
!--------------------------------------------------------------------------
    ! Local
    integer :: cols, rows, scale
    real :: Rg, phi, lam, rho
    real :: gamma, beta, epsilon, x, y
    real :: sinphi1, cosphi1
    real :: r0, s0, c

    if ((grid(1:1).eq.'N').or.(grid(1:1).eq.'S')) then
      cols = 721
      rows = 721
    else if (grid(1:1).eq.'M') then
      cols = 1383
      rows = 586
    else
      print *, 'ezlh_inverse: unknown projection: ', grid
      return
    endif
    
    if (grid(2:2).eq.'l') then
      scale = 1
    else if (grid(2:2).eq.'h') then
      scale = 2
    else
      print *, 'ezlh_inverse: unknown projection: ', grid
      return
    endif
    
    Rg = scale * RE_km/CELL_km
    
    r0 = (cols-1)/2. * scale
    s0 = (rows-1)/2. * scale
    
    x = r - r0
    y = -(s - s0)
    
    if ((grid(1:1).eq.'N').or.(grid(1:1).eq.'S')) then 
      rho = sqrt(x*x + y*y)
      if (rho.eq.0.0) then
        if (grid(1:1).eq.'N') lat = 90.0 
        if (grid(1:1).eq.'S') lat = -90.0 
        lon = 0.0
      else
        if (grid(1:1).eq.'N') then
          sinphi1 = sin(PI/2.)
          cosphi1 = cos(PI/2.)
          if (y.eq.0.) then
            if (r.le.r0) lam = -PI/2.
            if (r.gt.r0) lam = PI/2.
          else
            lam = atan2(x,-y)
          endif
        else if (grid(1:1).eq.'S') then
          sinphi1 = sin(-PI/2.)
          cosphi1 = cos(-PI/2.)
          if (y.eq.0.) then
            if (r.le.r0) lam = -PI/2.
            if (r.gt.r0) lam = PI/2.
          else
            lam = atan2(x,y)
          endif
        endif
        gamma = rho/(2 * Rg)
        if (abs(gamma).gt.1.) return
        c = 2 * asin(gamma)
        beta = cos(c) * sinphi1 + y * sin(c) * (cosphi1/rho)
        if (abs(beta).gt.1.) return
        phi = asin(beta)
        lat = deg(phi)
        lon = deg(lam)
      endif
    
    else if (grid(1:1).eq.'M') then
    !
    !  allow .5 cell tolerance in arcsin function
    !  so that grid coordinates which are less than .5 cells
    !  above 90.00N or below 90.00S are given a lat of 90.00
    !
      epsilon = 1 + 0.5/Rg
      beta = y*COS_PHI1/Rg
      if (abs(beta).gt.epsilon) return
      if (beta.le.-1.) then
        phi = -PI/2.
      else if (beta.ge.1.) then
        phi = PI/2.
      else
        phi = asin(beta)
      endif
      lam = x/COS_PHI1/Rg
      lat = deg(phi)
      lon = deg(lam)
    endif
    
  end subroutine

  real function rad(t)
    real :: t
    rad = t*PI/180.
    return
  end function rad


  real function deg(t)
    real :: t

    deg = t*180./PI
    return

  end function deg

end module convert_coord_module
