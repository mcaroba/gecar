module crystal_special

  contains





  subroutine is_crystal_implemented(crystal)

    implicit none

    character*64, intent(in) :: crystal
    character*64 :: crystal_library(1:1)
    integer :: i

!   List of implemented crystals
    crystal_library(1) = 'zincblende'

    do i = 1, len(crystal_library)
      if( crystal == crystal_library(i) )then
        return
      end if
    end do

!   If the crystal is not implemented we print an error and exit the program
    write(*,*)'                                       |'
    write(*,*)'ERROR: invalid special_lines argument  |  <-- ERROR'
    write(*,*)'                                       |'
    write(*,*)'.......................................|'
    stop

  end subroutine






  subroutine special_lines(crystal, n_lines, max_anchor_points, n_anchor_points, anchor_point, labels)

    implicit none

    character*64, intent(in) :: crystal
    integer, intent(inout) :: n_lines, max_anchor_points
    integer, allocatable, intent(inout) :: n_anchor_points(:)
    real*8, allocatable, intent(inout) :: anchor_point(:, :, :)
    character*1, allocatable, intent(inout) :: labels(:, :)

    if( crystal == 'zincblende' )then
!     For zincblende, the usual band plot is:
!       First line scan L - G - X - U
!         [[1./2., 1./2., 1./2.],  # L
!          [   0.,    0.,    0.],  # G
!          [   0., 1./2., 1./2.],  # X
!          [1./4., 5./8., 5./8.]], # U
!       Second line scan K - G
!         [[3./8., 3./4., 3./8.],  # K
!          [   0.,    0.,    0.]]  # G
      n_lines = 2
      max_anchor_points = 4
      allocate( n_anchor_points(1:n_lines) )
      allocate( anchor_point(1:3, 1:max_anchor_points, 1:n_lines) )
      allocate( labels(1:max_anchor_points, 1:n_lines) )
      n_anchor_points(1) = 4
      labels(1:4, 1) = (/ 'L', 'G', 'X', 'U' /)
      anchor_point( 1:3, 1, 1 ) = (/ 0.5d0, 0.5d0, 0.5d0 /)
      anchor_point( 1:3, 2, 1 ) = (/ 0.d0, 0.d0, 0.d0 /)
      anchor_point( 1:3, 3, 1 ) = (/ 0.d0, 0.5d0, 0.5d0 /)
      anchor_point( 1:3, 4, 1 ) = (/ 0.25d0, 0.625d0, 0.625d0 /)
      n_anchor_points(2) = 2
      labels(1:2, 2) = (/ 'K', 'G' /)
      anchor_point( 1:3, 1, 2 ) = (/ 0.375d0, 0.75d0, 0.375d0 /)
      anchor_point( 1:3, 2, 2 ) = (/ 0.d0, 0.d0, 0.d0 /)
    end if

  end subroutine

end module
