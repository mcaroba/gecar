debug="-g -fcheck=all"
mkdir -p bin/
gfortran -c printouts.f90 crystal_special.f90 functions.f90
gfortran $debug -o bin/gen_kpts gen_kpts.f90 *.o
gfortran $debug -o bin/get_bands get_bands.f90 *.o
