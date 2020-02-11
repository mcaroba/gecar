program gen_kpts

  use printouts
  use crystal_special

  implicit none


! Variable definitions
  character*1 :: keyword_first, crap
  character*1, allocatable :: labels(:, :)
  character*64 :: keyword, crystal
  character*16 :: routine
  character*1024 :: string1, string2, input_file, output_file
  integer :: iostatus, n_lines, max_anchor_points, n_anchor, iostatus2, i, ibeg, iend, j, nkpts, extra_kpts, k
  integer, allocatable :: n_anchor_points(:), linestart(:), anchor_point_position(:,:)
  real*8, allocatable :: anchor_point(:, :, :), kpts(:, :)
  logical :: do_special_lines = .false.
  real*8 :: kstep = 0.1d0, kdist, kdist_vec(1:3), ecut = 1.d0, dE_weight = 0.d0

  call print_welcome()

! This subroutine opens the input file as unit 10
  routine = 'gen_kpts'
  call check_input_file(routine)


! Read in options from input_kpts. On the first round we only look for the number of lines
  n_lines = 0
  max_anchor_points = 0
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
    else if( keyword=='special_lines' .and. .not. do_special_lines )then
      do_special_lines = .true.
      backspace(10)
      read(10, *) crap, crap, crystal
      call is_crystal_implemented(crystal)
    else if( keyword=='line' .and. .not. do_special_lines )then
      n_lines = n_lines + 1
      backspace(10)
      n_anchor = 1
      iostatus2 = 0
      do while(iostatus2 /= -2 ) 
        read(10, '(A1)', advance='no', iostat=iostatus2) keyword_first
        if( keyword_first == '-' )then
          n_anchor = n_anchor + 1
        end if
      end do
      if( n_anchor > max_anchor_points )then
        max_anchor_points = n_anchor
      end if
    end if
  end do
  if( n_lines == 0 .and. .not. do_special_lines )then
    write(*,*)'                                       |'
    write(*,*)'ERROR: you must define at least one    |  <-- ERROR'
    write(*,*)'line along k space                     |'
    write(*,*)'_______________________________________/'
    stop
  end if
  if( max_anchor_points < 2 .and. .not. do_special_lines )then
    write(*,*)'                                       |'
    write(*,*)'ERROR: you must define at least two    |  <-- ERROR'
    write(*,*)'anchor points per line along k space   |'
    write(*,*)'_______________________________________/'
    stop
  end if



! Choose between special lines (if implemented) and explicit line definitions
  if( do_special_lines )then
    call special_lines(crystal, n_lines, max_anchor_points, n_anchor_points, anchor_point, anchor_point_position, labels)
  else
    allocate( n_anchor_points(1:n_lines) )
    allocate( anchor_point(1:3, 1:max_anchor_points, 1:n_lines) )
    allocate( anchor_point_position( 1:max_anchor_points, 1:n_lines) )
    allocate( labels( 1:max_anchor_points, 1:n_lines) )
  end if


! Now we read in all the options
! Read in options from input_kpts
  rewind(10)
  iostatus = 0
  if( .not. do_special_lines )then
    n_lines = 0
  end if
  do while(iostatus==0)
    read(10, *, iostat=iostatus) keyword
    keyword = trim(keyword)
    if(iostatus/=0)then
      exit
    end if
    keyword_first = keyword
    if(keyword_first=='#' .or. keyword_first=='!')then
      continue
    else if( keyword == 'kstep' )then
      backspace(10)
      read(10,*) crap, crap, kstep
    else if( keyword == 'ecut' )then
      backspace(10)
      read(10,*) crap, crap, ecut
    else if( keyword == 'dE_weight' )then
      backspace(10)
      read(10,*) crap, crap, dE_weight
    else if( keyword == 'input_file' )then
      backspace(10)
      read(10,*) crap, crap, input_file
    else if( keyword == 'output_file' )then
      backspace(10)
      read(10,*) crap, crap, output_file
    else if(keyword=='line' .and. .not. do_special_lines )then
      n_lines = n_lines + 1
      backspace(10)
      i = 0
      do
        i = i + 1
        read(10, '(A1)', advance='no') keyword_first
        if( keyword_first == '=' )then
          ibeg = i + 1
          backspace(10)
          exit
        end if
      end do
      read(10,'(A)') string1
      i = 0
      n_anchor = 1
      read(string1(ibeg:), *) anchor_point(1:3, n_anchor, n_lines)
      do while( i < 1024 )
        i = i + 1
        read(string1(i:i), '(A1)') keyword_first
        if( keyword_first == '-' )then
          ibeg = i + 1
          n_anchor = n_anchor + 1
          read(string1(ibeg:), *) anchor_point(1:3, n_anchor, n_lines)
        end if
      end do
      n_anchor_points(n_lines) = n_anchor
      read(10, *, iostat=iostatus2) keyword
      backspace(10)
      if( keyword == 'labels' )then
        read(10, *) crap, crap, labels(1:n_anchor, n_lines)
      else
        continue
      end if
    end if
  end do
  close(10)




! For debugging
  if( .false. )then
    do i = 1, n_lines
      write(*,*) labels(1:n_anchor_points(i), i)
      do j = 1, n_anchor_points(i)
        write(*,*) anchor_point(1:3, j, i)
      end do
    end do
  end if






  write(*,*)'                                       |'
  write(*,*)'Generating k points along specified    |'
  write(*,*)'lines...                               |'
  write(*,*)'                                       |'
  write(*,*)'.......................................|'

  allocate( linestart(1:n_lines) )
  linestart(1) = 1
! This is just to count the number of kpts
  nkpts = 0
  do i = 1, n_lines
    do j = 1, n_anchor_points(i)-1
      nkpts = nkpts+1
      anchor_point_position(j, i) = nkpts
      kdist = sqrt( (anchor_point(1, j+1, i) - anchor_point(1, j, i))**2 + &
                    (anchor_point(2, j+1, i) - anchor_point(2, j, i))**2 + &
                    (anchor_point(3, j+1, i) - anchor_point(3, j, i))**2 )
      extra_kpts = nint(kdist/kstep) - 1
      if( extra_kpts > 0 )then
        nkpts = nkpts + extra_kpts
      end if
    end do
    nkpts = nkpts + 1
    anchor_point_position(n_anchor_points(i), i) = nkpts
    if( i < n_lines )then
      linestart(i+1) = nkpts + 1
    end if
  end do
! This is to compute the kpoint coordinates
  allocate( kpts(1:3, 1:nkpts) )
  nkpts = 0
  do i = 1, n_lines
    do j = 1, n_anchor_points(i)-1
      nkpts = nkpts+1
      kpts(1:3, nkpts) = anchor_point(1:3, j, i)
      kdist = sqrt( (anchor_point(1, j+1, i) - anchor_point(1, j, i))**2 + &
                    (anchor_point(2, j+1, i) - anchor_point(2, j, i))**2 + &
                    (anchor_point(3, j+1, i) - anchor_point(3, j, i))**2 )
      extra_kpts = nint(kdist/kstep) - 1
      if( extra_kpts > 0 )then
        kdist_vec(1:3) = (/ anchor_point(1, j+1, i) - anchor_point(1, j, i), &
                            anchor_point(2, j+1, i) - anchor_point(2, j, i), &
                            anchor_point(3, j+1, i) - anchor_point(3, j, i) /) &
                         / dfloat(extra_kpts + 1)

        do k = 1, extra_kpts
          nkpts = nkpts + 1
          kpts(1:3, nkpts) = anchor_point(1:3, j, i) + dfloat(k) * kdist_vec(1:3)
        end do
      end if
    end do
    nkpts = nkpts+1
    kpts(1:3, nkpts) = anchor_point(1:3, n_anchor_points(i), i)
  end do




! For debugging
  if( .false. )then
    write(*,*) linestart
    do i=1,nkpts
      write(*,*) kpts(1:3,i)
    end do
  end if






! Create output
  open(unit=10, file=input_file, status='old')
  open(unit=20, file=output_file, status='unknown')
  read(10,*)
  write(20,'(A)') 'This KPOINTS file was generated with GECAR''s gen_kpts code'
  read(10,*) k
  write(20,'(I8)') k+nkpts
  do i = 1, k + 1
    read(10,'(A)') string1
    write(20,'(A)') trim(string1)
  end do
  do i = 1, nkpts
    write(20,'(F20.14,F20.14,F20.14,I14)') kpts(1:3, i), 0
  end do
  close(10)
  close(20)





! Check point for the get_bands program to take over
  open(unit=10, file='checkpoint', status='unknown')
  write(10,*) 'Do not delete this file! It''s needed for the second GECAR pass (get_bands)!'
  write(10, *) n_lines, nkpts+k
  write(10, *) linestart(1:n_lines)+k
  write(10, *) ecut, dE_weight
  write(10, *) max_anchor_points, n_anchor_points(1:n_lines)
  do i = 1, n_lines
    write(10, *) anchor_point_position(1:max_anchor_points, i) + k
  end do
  do i = 1, n_lines
    do j = 1, n_anchor_points(i)
      write(10,'(A1)',advance='no') labels(j, i)
    end do
  end do
  close(10)











  write(*,*)'                                       |'
  write(*,*)'Code gen_kpts successfully executed.   |'
  write(*,*)'                                       |'
  write(*,*)'If you want to do the band disentan-   |'
  write(*,*)'glement based on band projection, run  |'
  write(*,*)'a VASP calculation with the output     |'
  write(*,*)'file as KPOINTS and then run get_bands |'
  write(*,*)'with the resulting PROCAR (LORBIT = 10 |'
  write(*,*)'or LORBIT = 11 in your INCAR file).    |'
  write(*,*)'_______________________________________/'



end program gen_kpts
