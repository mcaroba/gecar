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
      input_file = 'input_bands'
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





  subroutine create_gnuplot_script(n_lines, n_anchor_points, nkpts, labels, linestart, gnuplot_tags, gnuplot_tags_logical)

    implicit none

    integer, intent(in) :: n_lines, nkpts, linestart(:), n_anchor_points(:)
    character*1, intent(in) :: labels(:,:)
    character*1024, intent(in) :: gnuplot_tags(:)
    logical, intent(in) :: gnuplot_tags_logical(:)
    character*1024 :: title, ylabel, xlabel, style, term_size


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

    write(*,*)'                                       |'
    write(*,*)'Generating plot_bands.gnuplot script...|'
    write(*,*)'                                       |'
    write(*,*)'.......................................|'

    open(unit=10, file='plot_bands.gnuplot', status='unknown')
    write(10, '(A)') 'set term pngcairo'
    write(10, '(A)') 'set output bands.png'
    write(10, '(A,1X,A)') 'set size', adjustl(trim(term_size))
    write(10, '(A,1X,A,A,A)') 'set title', '"', adjustl(trim(title)), '"'
    write(10, '(A,1X,A,A,A)') 'set xlabel', '"', adjustl(trim(xlabel)), '"'
    write(10, '(A,1X,A,A,A)') 'set ylabel', '"', adjustl(trim(ylabel)), '"'

    close(10)

  end subroutine



end module
