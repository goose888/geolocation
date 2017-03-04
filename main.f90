program easelatlon
  use convert_coord_module
  implicit none

  ! ==========================================================================
  !  MAIN FORTRAN routines for conversion of azimuthal 
  ! 		equal area and equal area cylindrical grid coordinates
  !  Shijie: F90 version
  !
  ! ==========================================================================

  !--------------------------------------------------------------------------
  !
  ! convert geographic coordinates (spherical earth) to 
  ! azimuthal equal area or equal area cylindrical grid coordinates
  !
  ! status = ezlh_convert (grid, lat, lon, r, s)
  !
  ! input : grid - projection name '[NSM][lh]'
  !               where l = "low"  = 25km resolution
  !                     h = "high" = 12.5km resolution
  !         lat, lon - geo. coords. (decimal degrees)
  !
  !         output: r, s - column, row coordinates
  !
  !         result: status = 0 indicates normal successful completion
  !                    -1 indicates error status (point not on grid)
  !
  !--------------------------------------------------------------------------

  ! Local Vars
  integer :: ierr, i
  character(1) :: region          ! 'N', 'S', 'M'
  character(1) :: res             ! 'h', 'l'
  character(256) :: latfile       ! Name of ascii file storing latitudes
  character(256) :: lonfile       ! Name of ascii file storing longitudes
  character(256) :: outlat        ! Name of ascii file storing latitudes after transformation
  character(256) :: outlon        ! Name of ascii file storing longitudes after transformation
  integer :: numrows, numcols     ! 'number of rows and cols'
  logical :: tolatlon             ! .true.  -> transfer from row/col to lon/lat
                                  ! .false. -> transfer from lon/lat to col/row
  real :: temp

  real, dimension(:), allocatable :: lons
  real, dimension(:), allocatable :: lats
  real, dimension(:), allocatable :: reslons
  real, dimension(:), allocatable :: reslats

  !----------------------------------
  ! Namelist variables for ISAM_STANDALONE
  !----------------------------------
  namelist /easecfg/ region, res, numrows, numcols, tolatlon, &
                  & latfile, lonfile, outlat, outlon

  !----------------------------------
  ! Set Defaults for ISAM_STANDALONE
  !----------------------------------
  region = 'M'
  res = 'l'
  numrows = 586
  numcols = 1383
  tolatlon = .true.
  latfile = 'lat.dat'
  lonfile = 'lon.dat'
  outlat  = 'outlat.dat'
  outlon  = 'outlon.dat'

  !----------------------------------
  ! Read Namelist for ISAM_STANDALONE
  !----------------------------------
  read (*, easecfg, iostat=ierr)
  if (ierr /= 0) then
     print *, 'Main Program: Error reading namelist.  Aborting'
     stop 2
  end if

  !--------------------------------------
  ! Allocate the size of the coordinates
  !--------------------------------------
  allocate(lats(numrows), stat=ierr)
  if (ierr /= 0) then
     print *, 'Main Program: Error allocating lats.  Aborting'
     stop 2
  end if

  allocate(lons(numcols), stat=ierr)
  if (ierr /= 0) then
     print *, 'Main Program: Error allocating lons.  Aborting'
     stop 2
  end if

  allocate(reslats(numrows), stat=ierr)
  if (ierr /= 0) then
     print *, 'Main Program: Error allocating reslats.  Aborting'
     stop 2
  end if
  reslats = -99.

  allocate(reslons(numcols), stat=ierr)
  if (ierr /= 0) then
     print *, 'Main Program: Error allocating reslons.  Aborting'
     stop 2
  end if
  reslons = -99.

  !--------------------------------------
  ! Read in lats and lons
  !--------------------------------------
  open(unit=08, file=trim(latfile), status='OLD', iostat=ierr)
  if (ierr /= 0) then
     print *, 'Main Program: Error reading latfile.  Aborting'
     stop 2
  end if
  do i = 1, numrows
     read(08, *) lats(i)
  end do
  close(08)

  open(unit=09, file=trim(lonfile), status='OLD', iostat=ierr)
  if (ierr /= 0) then
     print *, 'Main Program: Error reading lonfile.  Aborting'
     stop 2
  end if
  do i = 1, numcols
     read(09, *) lons(i)
  end do
  close(09)

  !--------------------------------------
  ! Transform the coordination system
  !--------------------------------------
  do i = 1, numrows
     if(tolatlon) then
        call ezlh_inverse(region//res, lons(1), lats(i), reslats(i), temp)!, temp)
     else
        call ezlh_convert(region//res, lats(i), lons(1), temp, reslats(i))
     end if
  end do
  
  do i = 1, numcols
     if(tolatlon) then
        call ezlh_inverse(region//res, lons(i), lats(1), temp, reslons(i))
     else
        call ezlh_convert(region//res, lats(1), lons(i), reslons(i), temp)
     end if
  end do

  !--------------------------------------
  ! Write the output
  !--------------------------------------
  open(unit=08, file=trim(outlat), status='unknown', iostat=ierr)
  if (ierr /= 0) then
     print *, 'Main Program: Error writing outlat file.  Aborting'
     stop 2
  end if
  do i = 1, numrows
     write(08, *) reslats(i)
  end do
  close(08)

  open(unit=09, file=trim(outlon), status='unknown', iostat=ierr)
  if (ierr /= 0) then
     print *, 'Main Program: Error writing outlon file.  Aborting'
     stop 2
  end if
  do i = 1, numcols
     write(09, *) reslons(i)
  end do
  close(09)

  print *, 'Finish transformation!'

end program easelatlon
