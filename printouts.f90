module printouts

  contains

  subroutine print_welcome()

    implicit none

    write(*,*)'_________________________________________________________________ '
    write(*,*)'                                                                 \'
    write(*,*)'                             GECAR                               |'
    write(*,*)'           ( GEt Crossings and Anticrossings Right )             |'
    write(*,*)'                                                                 |'
    write(*,*)'    A light-weight tool to plot electronic band diagrams from    |'
    write(*,*)'    VASP output, preserving the band characted across adjacent   |'
    write(*,*)'            k points, according to band projections              |'
    write(*,*)'                                                                 |'
    write(*,*)'                     Welcome to GECAR v0.1                       |'
    write(*,*)'                    Written and maintained by                    |'
    write(*,*)'                                                                 |'
    write(*,*)'                         Miguel A. Caro                          |'
    write(*,*)'                       mcaroba@gmail.com                         |'
    write(*,*)'                      miguel.caro@aalto.fi                       |'
    write(*,*)'                                                                 |'
    write(*,*)'       Department of Electrical Engineering and Automation       |'
    write(*,*)'                     Aalto University, Finland                   |'
    write(*,*)'                                                                 |'
    write(*,*)'                     Last updated: Jan. 2020                     |'
    write(*,*)'                                        _________________________/'
    write(*,*)'.......................................|'

  end subroutine






  subroutine check_input_file(routine)

    implicit none

    character*16 :: routine
    character*128 :: input_file
    integer :: iostatus

    if( routine == 'gen_kpts' )then
      input_file = 'input_kpts'
    else if( routine == 'get_bands' )then
!      input_file = 'input_bands'
      input_file = 'input_kpts'
    else
      write(*,*) 'Error in check_input_file() subroutine'
    end if

!   Print message about mode of operation
    if( routine == 'gen_kpts' )then
      write(*,*)'                                       |'
      write(*,*)'You''re running the gen_kpts routine:   |'
    else if( routine == 'get_bands' )then
      write(*,*)'                                       |'
      write(*,*)'You''re running the get_bands routine:  |'
    end if

!   Read input file
    open(unit=10, file=input_file, status='old', iostat=iostatus)
!   Check for existence of input file
    write(*,*)'                                       |'
    write(*,*)'Checking input file...                 |'
    if(iostatus/=0)then
      close(10)
      write(*,*)'                                       |'
      write(*,*)'ERROR: input file could not be found   |  <-- ERROR'
      write(*,*)'_______________________________________/'
      stop
  end if

  end subroutine





  subroutine create_gnuplot_script(n_lines, nbands, n_anchor_points, anchor_point_position, nkpts, labels, &
                                   linestart, gnuplot_tags, gnuplot_tags_logical)

    implicit none

    integer, intent(in) :: n_lines, nkpts, linestart(:), n_anchor_points(:), nbands, anchor_point_position(:,:)
    character*1, intent(in) :: labels(:,:)
    character*1024, intent(in) :: gnuplot_tags(:)
    logical, intent(in) :: gnuplot_tags_logical(:)
    character*1024 :: title, ylabel, xlabel, style, term_size
    real*8 :: left_margin, right_margin, bottom_margin, top_margin, space, line_spacer, space_per_k_point
    integer :: i, j, border

    if( gnuplot_tags_logical(1) )then
      title = adjustl(trim(gnuplot_tags(1)))
    end if
    if( gnuplot_tags_logical(2) )then
      xlabel = adjustl(trim(gnuplot_tags(2)))
    end if
    if( gnuplot_tags_logical(3) )then
      ylabel = adjustl(trim(gnuplot_tags(3)))
    end if
    if( gnuplot_tags_logical(4) )then
      style = adjustl(trim(gnuplot_tags(4)))
    else
      style = 'plain'
    end if
    if( gnuplot_tags_logical(5) )then
      term_size = adjustl(trim(gnuplot_tags(5)))
    else
      term_size = '640,480'
    end if

    left_margin = 0.12d0
    right_margin = 0.98d0
    top_margin = 0.98d0
    bottom_margin = 0.12d0
!    line_spacer = 0.02d0
    line_spacer = 0.d0
    space = right_margin - left_margin
    space = space - line_spacer * dfloat(n_lines-1)
    space_per_k_point = space/dfloat(nkpts-n_lines)

    write(*,*)'                                       |'
    write(*,*)'Generating plot_bands.gnuplot script...|'
    write(*,*)'                                       |'
    write(*,*)'.......................................|'

    open(unit=10, file='plot_bands.gnuplot', status='unknown')
    write(10, '(A)') 'set term pngcairo'
    write(10, '(A)') 'set output "bands.png"'
    write(10, '(A,1X,A)') 'set size', adjustl(trim(term_size))
    write(10, '(A,1X,A,A,A)') 'set title', '"', adjustl(trim(title)), '"'
    write(10, '(A,1X,A,A,A,F4.2,A,F4.2,A)') 'set label 1', '"', adjustl(trim(xlabel)), '" at screen ', &
          (left_margin+right_margin)/2.d0, ',', bottom_margin/3.d0, ' center'
    write(10, '(A,1X,A,A,A)') 'set ylabel', '"', adjustl(trim(ylabel)), '"'
    write(10, '(A)') 'set multiplot'
    do i = 1, n_lines
      do j = 1, n_anchor_points(i)
        if( j == 1 )then
          if( i > 1 .and. line_spacer == 0.d0 )then
            write(10, '(A,A1,A,A1,A,I0,A)') 'set xtics ("', labels(j, i), ',', labels(n_anchor_points(i-1), i-1), &
                                       '" ', anchor_point_position(j,i), ')'
          else
            write(10, '(A,A1,A,I0,A)') 'set xtics ("', labels(j, i), '" ', anchor_point_position(j,i), ')'
          end if
        else if( j == n_anchor_points(i) .and. i < n_lines .and. line_spacer == 0.d0 )then
          continue
        else
          write(10, '(A,A1,A,I0,A)') 'set xtics add ("', labels(j, i), '" ', anchor_point_position(j,i), ')'
        end if
      end do
      if( i == 2 )then
        write(10, '(A)') 'unset ylabel'
        write(10, '(A)') 'unset xlabel'
        write(10, '(A)') 'unset label 1'
        write(10, '(A)') 'set format y ""'
      end if
!     b + l + t + r
      border = 1+2+4+8
      if( line_spacer == 0.d0 )then
        if( i < n_lines )then
          border = border - 8
        end if
        if( i > 1 )then  
          border = border - 2
        end if
        write(10, '(A,I0)') 'set border ', border
      end if
      if( i == n_lines )then
        write(10, '(A,I0,A,I0,A)') 'set xrange [', linestart(i), ':', nkpts, ']'
        right_margin = left_margin + dfloat(nkpts - linestart(i)) * space_per_k_point
      else
        write(10, '(A,I0,A,I0,A)') 'set xrange [', linestart(i), ':', linestart(i+1)-1, ']'
        right_margin = left_margin + dfloat(linestart(i+1) - linestart(i) - 1) * space_per_k_point
      end if
      write(10, '(A,1X,F4.2)') 'set lmargin at screen', left_margin
      write(10, '(A,1X,F4.2)') 'set rmargin at screen', right_margin
      write(10, '(A,1X,F4.2)') 'set tmargin at screen', top_margin
      write(10, '(A,1X,F4.2)') 'set bmargin at screen', bottom_margin
      left_margin = right_margin + line_spacer
      write(10, '(A,I0,A,I0,A)') 'plot for [i=0:', nbands-1, '] "bands_', i, '.dat" every :::i::i u 1:5 w l lc i not'
    end do

    close(10)

  end subroutine



end module
