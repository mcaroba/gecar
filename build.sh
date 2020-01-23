#debug="-g -fcheck=all"
gfortran -c printouts.f90 crystal_special.f90
gfortran -o gen_kpts gen_kpts.f90 *.o
gfortran $debug -o get_bands get_bands.f90
