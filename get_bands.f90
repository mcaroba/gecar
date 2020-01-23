program get_bands

  implicit none

! norbs should be changed to accept any number, depending on the format of the PROCAR file !!!!!!!!!!!!!!!!!!!!!!!
  integer :: nk, nbands, kstart, nions, norbs, orb, k, band, j, iostatus, ijunk, ion
  integer :: band2, band_match, band3, kend, ktotal
  real*8, allocatable :: proj(:,:,:,:), energy(:,:), mx(:,:,:,:), my(:,:,:,:), mz(:,:,:,:)
  real*8, allocatable :: kpoint(:,:)
  character*64 :: cjunk, filename, output_filename
  integer, allocatable :: correct_band(:,:), ordered_band(:,:)
  real*8 :: ecut, dist, dist_prev, dE, dE_prev, n1, n2, dE_weight, ene, dk, dE2, dk2, dk3, dE3
  logical, allocatable :: band_assigned(:)
  logical :: spin_orbit = .false.
  real*8 :: proj_w = 1.d0
  real*8, allocatable :: proj_total(:,:,:)
  integer :: iostatus2, nkpts, n_lines, line
  integer, allocatable :: linestart(:)

!  read(*,*) kstart, kend, filename, ecut, dE_weight
!  nk = kend - kstart + 1

! Read input
  filename = 'PROCAR'
  open(unit=10, file='checkpoint', status='old')
  read(10, *)
  read(10, *) n_lines, nkpts
  allocate( linestart(1:n_lines) )
  read(10,*) linestart(1:n_lines)
  read(10,*) ecut, dE_weight
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
  band_assigned = .false.

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
      read(10, fmt="(I1)", iostat=iostatus2) ijunk
      if( ijunk == 1 )then
        spin_orbit = .true.
        backspace(10)
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
    band_assigned = .false.
    do band = 1, nbands
!     For each band at the present k point, we compute the similarity with every band in the previous
!     k point provided that the two energy values are less than ecut apart
      dist_prev = 1.d100
      do band2 = 1, nbands
        ene = energy(k-1, band2)
        dE = dsqrt( (ene - energy(k, band))**2 )
        if( dE < ecut)then
          dist = 0.d0
!          n1 = 0.d0
!          n2 = 0.d0
          do orb = 1, norbs
            do ion = 1, nions
!             This overlap scheme gives problems
!              n1 = n1 + proj_w*proj(orb, ion, k-1, band2)**2 + mx(orb, ion, k-1, band2)**2 + &
!                        my(orb, ion, k-1, band2)**2 + mz(orb, ion, k-1, band2)**2
!              n2 = n2 + proj_w*proj(orb, ion, k, band)**2 + mx(orb, ion, k, band)**2 + &
!                        my(orb, ion, k, band)**2 + mz(orb, ion, k, band)**2
!              dist = dist + proj_w*proj(orb, ion, k-1, band2) * proj(orb, ion, k, band) + &
!                            mx(orb, ion, k-1, band2) * mx(orb, ion, k, band) + &
!                            my(orb, ion, k-1, band2) * my(orb, ion, k, band) + &
!                            mz(orb, ion, k-1, band2) * mz(orb, ion, k, band)
              dist = dist + proj_w*(proj(orb, ion, k-1, band2) - proj(orb, ion, k, band))**2 + &
                            (mx(orb, ion, k-1, band2) - mx(orb, ion, k, band))**2 + &
                            (my(orb, ion, k-1, band2) * my(orb, ion, k, band))**2 + &
                            (mz(orb, ion, k-1, band2) * mz(orb, ion, k, band))**2
            end do
          end do
          dist = dist + dE_weight*dE
!          if( n1 > 0.d0 .and. n2 > 0.d0 )then
!            dist = 1.d0 - dist/dsqrt(n1)/dsqrt(n2) + dE_weight*dE
!          else
!            dist = 1.d0 + dE_weight*dE
!          end if
          band3 = ordered_band(k-1, band2)
          if( dist < dist_prev .and. (.not. band_assigned(band3)) )then
            dist_prev = dist
            dE_prev = dE
            band_match = band3
          end if
        end if
      end do
      correct_band(k, band_match) = band
      ordered_band(k, band) = band_match
      band_assigned(band_match) = .true.
    end do
  end do


! Check that all bands are assigned
  do k = 1, nk
    do band = 1, nbands
      if( correct_band(k, band) == 0 )then
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




  deallocate( proj, proj_total, mx, my, mz, kpoint, correct_band, ordered_band, energy, band_assigned )

  end do

end program
