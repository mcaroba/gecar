program get_bands

  use printouts

  implicit none

! norbs should be changed to accept any number, depending on the format of the PROCAR file !!!!!!!!!!!!!!!!!!!!!!!
  integer :: nk, nbands, kstart, nions, norbs, orb, k, band, j, iostatus, ijunk, ion, i
  integer :: band2, band_match, band3, kend, ktotal
  real*8, allocatable :: proj(:,:,:,:), energy(:,:), mx(:,:,:,:), my(:,:,:,:), mz(:,:,:,:)
  real*8, allocatable :: kpoint(:,:), dist_matrix(:,:)
  character*64 :: cjunk, filename, output_filename
  integer, allocatable :: correct_band(:,:), ordered_band(:,:)
  real*8 :: ecut, dist, dist_prev, dE, dE_prev, n1, n2, dE_weight, ene, dk(1:3), dE2, dk2(1:3), dk3, dE3, dEdk
  logical, allocatable :: band_assigned(:), band2_assigned(:)
  logical :: spin_orbit = .false.
  real*8 :: proj_w = 1.d0
  real*8, allocatable :: proj_total(:,:,:)
  integer :: iostatus2, nkpts, n_lines, line
  integer, allocatable :: linestart(:), n_anchor_points(:), anchor_point_position(:,:)
  character*64 :: keyword
  character*1 :: keyword_first, crap
  character*1, allocatable :: labels(:, :)
  character*16 :: routine
  character*1024 :: gnuplot_tags(1:5)
  logical :: gnuplot, gnuplot_tags_logical(1:5) = .false., use_derivative
  integer :: max_anchor_points, band_pick, band2_pick


!  read(*,*) kstart, kend, filename, ecut, dE_weight
!  nk = kend - kstart + 1

! Read input
  filename = 'PROCAR'
  open(unit=10, file='checkpoint', status='old')
  read(10, *)
  read(10, *) n_lines, nkpts
  allocate( linestart(1:n_lines) )
  allocate( n_anchor_points(1:n_lines) )
  read(10,*) linestart(1:n_lines)
  read(10,*) ecut, dE_weight
  read(10,*) max_anchor_points, n_anchor_points(1:n_lines)
  allocate( anchor_point_position( 1:max_anchor_points, 1:n_lines) )
  allocate( labels( 1:max_anchor_points, 1:n_lines) )
  do i = 1, n_lines
    read(10, *) anchor_point_position(1:max_anchor_points, i)
  end do
  do i = 1, n_lines
    do j = 1, n_anchor_points(i)
      read(10,'(A1)',advance='no',iostat=iostatus) labels(j, i)
    end do
  end do
  close(10)



! This subroutine opens the input file as unit 10
  routine = 'get_bands'
  call check_input_file(routine)
  iostatus = 0
  do while(iostatus==0)
    read(10, *, iostat=iostatus) keyword
    keyword = trim(keyword)
    if(iostatus/=0)then
      exit
    end if
    keyword_first = keyword
    if( keyword_first=='#' .or. keyword_first=='!' )then
      continue
    else if( keyword == 'gnuplot' )then
      backspace(10)
      read(10,*) crap, crap, gnuplot
    else if( keyword == 'gnuplot_title' )then
      backspace(10)
      read(10,*) crap, crap, gnuplot_tags(1)
      gnuplot_tags_logical(1) = .true.
    else if( keyword == 'gnuplot_xlabel' )then
      backspace(10)
      read(10,*) crap, crap, gnuplot_tags(2)
      gnuplot_tags_logical(2) = .true.
    else if( keyword == 'gnuplot_ylabel' )then
      backspace(10)
      read(10,*) crap, crap, gnuplot_tags(3)
      gnuplot_tags_logical(3) = .true.
    else if( keyword == 'gnuplot_style' )then
      backspace(10)
      read(10,*) crap, crap, gnuplot_tags(4)
      gnuplot_tags_logical(4) = .true.
    else if( keyword == 'gnuplot_size' )then
      backspace(10)
      read(10,*) crap, crap, gnuplot_tags(5)
      gnuplot_tags_logical(5) = .true.
    end if
  end do
  close(10)



  do line = 1, n_lines

  kstart = linestart(line)
  if( n_lines == 1 .or. line == n_lines )then
    kend = nkpts
  else if( line /= n_lines )then
    kend = linestart(line+1) - 1
  end if
  nk = kend - kstart + 1

  open(unit=10, file=filename, status='old')

  read(10,*) cjunk, cjunk
  if( cjunk == 'new' )then
    norbs = 3
  else
    norbs = 9
  end if
  read(10,*) cjunk, cjunk, cjunk, ijunk, cjunk, cjunk, cjunk, nbands, cjunk, cjunk, cjunk, nions

  allocate( proj(1:norbs, 1:nions, 1:nk, 1:nbands) )
  allocate( proj_total(1:norbs, 1:nk, 1:nbands) )
  allocate( mx(1:norbs, 1:nions, 1:nk, 1:nbands) )
  allocate( my(1:norbs, 1:nions, 1:nk, 1:nbands) )
  allocate( mz(1:norbs, 1:nions, 1:nk, 1:nbands) )
  proj = 0.d0
  proj_total = 0.d0
  proj_w = 1.d0
  mx = 0.d0
  my = 0.d0
  mz = 0.d0

  allocate( kpoint(1:nk, 3) )

  allocate( correct_band(1:nk, 1:nbands) )
  allocate( ordered_band(1:nk, 1:nbands) )
  correct_band = 0
  ordered_band = 0
  allocate( energy(1:nk, 1:nbands) )
  energy = 0.d0
  allocate( band_assigned(1:nbands) )
  allocate( band2_assigned(1:nbands) )
  allocate( dist_matrix(1:nbands, 1:nbands) )
  band_assigned = .false.
  band2_assigned = .false.

  iostatus = 0
  k = 0
  ktotal = 0
  do while( iostatus /= -1 )
    read(10, *, iostat=iostatus) cjunk, ijunk
    if( cjunk == "k-point" )then
      ktotal = ktotal + 1
      if( ktotal > kend )then
        exit
      end if
      k = ktotal - kstart + 1
      if( k > 0)then
        backspace(10)
        read(10, *) cjunk, ijunk, cjunk, kpoint(k, 1:3)
      end if
    end if
    if( cjunk == "band" .and. k > 0 )then
      backspace(10)
      read(10, *) cjunk, band, cjunk, cjunk, energy(k, band)
      read(10, *)
      read(10, *)
      do ion = 1, nions
        read(10, *) ijunk, proj(1:norbs, ion, k, band)
      end do
      read(10,*) cjunk, proj_total(1:norbs, k, band)
      ijunk = 0
      read(10, *, iostat=iostatus2) ijunk
      backspace(10)
      if( ijunk == 1 )then
        spin_orbit = .true.
        do ion = 1, nions
          read(10, *) ijunk, mx(1:norbs, ion, k, band)
        end do
        read(10, *)
        do ion = 1, nions
          read(10, *) ijunk, my(1:norbs, ion, k, band)
        end do
        read(10, *)
        do ion = 1, nions
          read(10, *) ijunk, mz(1:norbs, ion, k, band)
        end do
        read(10, *)
      end if
    end if
  end do
  close(10)


! If the magnetic projections are available then forget about the total projections and
! use those instead
  if( spin_orbit )then
!    proj = 0.d0
    proj_w = 0.d0
  end if


! Keep order consistent with band order at first k point
  do band = 1, nbands
    correct_band(1, band) = band
    ordered_band(1, band) = band
  end do

! Now we check for the following k points which band resembles the most the bands at the
! previous k point
  do k = 2, nk
!   Predict energy based on derivative dE/dk. Here we obtain dk just to check whether there
!   has been a change in k-space direction
    if( k == 2 )then
      use_derivative = .false.
    else
      dk = (kpoint(k, 1:3) - kpoint(k-1, 1:3)) / dsqrt(dot_product(kpoint(k, 1:3)-kpoint(k-1, 1:3), &
                                                                   kpoint(k, 1:3)-kpoint(k-1, 1:3)))
      dk2 = (kpoint(k-1, 1:3) - kpoint(k-2, 1:3)) / dsqrt(dot_product(kpoint(k-1, 1:3)-kpoint(k-2, 1:3), &
                                                                   kpoint(k-1, 1:3)-kpoint(k-2, 1:3)))
      if( dot_product(dk, dk2) < 0.999999d0 )then
        use_derivative = .false.
      else
      use_derivative = .true.
      end if
    end if
!   Create a distance matrix
    band_assigned = .false.
    band2_assigned = .false.
    dist_matrix = 1.d99
    do band = 1, nbands
      do band2 = 1, nbands
!       Predict the expected energy
        if( use_derivative )then
          band3 = ordered_band(k-1, band2)
!         Since dk is fixed, E + dE/dk * dk has this expression:
          ene = 2.d0*energy(k-1, band2) - energy(k-2, correct_band(k-2, band3))
        else
          ene = energy(k-1, band2)
        end if
        dE = dsqrt( (ene - energy(k, band))**2 )
        if( dE < ecut)then
          dist = 0.d0
          do ion = 1, nions
            do orb = 1, norbs
              dist = dist + proj_w*(proj(orb, ion, k-1, band2) - proj(orb, ion, k, band))**2 + &
                            (mx(orb, ion, k-1, band2) - mx(orb, ion, k, band))**2 + &
                            (my(orb, ion, k-1, band2) - my(orb, ion, k, band))**2 + &
                            (mz(orb, ion, k-1, band2) - mz(orb, ion, k, band))**2
            end do
          end do
          dist = dist + dE_weight*dE**2
           dist_matrix(band, band2) = dist
        end if
      end do
    end do
!   Now we find the minima in the distance matrix to assign the bands
!    do while( .not. all(band_assigned) )
    do i = 1, nbands
      dist_prev = 1.d100
      do band = 1, nbands
        if( .not. band_assigned(band) )then
          do band2 = 1, nbands
            if( .not. band2_assigned(band2) )then
              dist = dist_matrix(band, band2)
              if( dist < dist_prev )then
                band_pick = band
                band2_pick = band2
                dist_prev = dist
              end if
            end if
          end do
        end if
      end do
      band_match = ordered_band(k-1, band2_pick)
      correct_band(k, band_match) = band_pick
      ordered_band(k, band_pick) = band_match
      band_assigned(band_pick) = .true.
      band2_assigned(band2_pick) = .true.
    end do
  end do



! Check that all bands are assigned
  do k = 1, nk
    do band = 1, nbands
      if( correct_band(k, band) == 0 .or. correct_band(k, band) > nbands )then
        write(*,*)'                                       |'
        write(*,*)'ERROR: One or more bands could not be  |  <-- ERROR'
        write(*,*)'assigned. This error is usually        |'
        write(*,*)'related to small search windows (i.e., |'
        write(*,*)'try increasing the value of ecut).     |'
        write(*,*)'_______________________________________/'
        stop
      end if
    end do
  end do


! Output
  write(output_filename, '(A,I0,A)') 'bands_', line, '.dat'
  open(unit=10, file=output_filename, status='unknown')
  do band = 1, nbands
    write(10,*) "# Band no.", band
    do k = 1, nk
      write(10,*) k-1 + kstart, kpoint(k,1:3), energy(k, correct_band(k, band)), &
                  proj_total(1:norbs, k, correct_band(k, band))
    end do
    write(10,*)
  end do
  close(10)




  deallocate( proj, proj_total, mx, my, mz, kpoint, correct_band, ordered_band, energy, band_assigned, &
              band2_assigned, dist_matrix )

  end do



! Generate gnuplot script
 if( gnuplot )then
   call create_gnuplot_script(n_lines, nbands, n_anchor_points, anchor_point_position, nkpts, &
                              labels, linestart, gnuplot_tags, gnuplot_tags_logical)
 end if



end program
