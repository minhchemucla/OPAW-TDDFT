## compiler options
FC   	= mpif90 -O3 -g -fcheck=all -fbacktrace #-Wall
#FC   	= mpif90 -g -ffpe-trap=zero,invalid,overflow,underflow -fcheck=all -fbacktrace
#FC	= mpif90 -Ofast
CC		= gcc 
#FC      = /home/wenfei/codes/openmpi-4.0.1/build-lahey/bin/mpif90 --chkglobal
FFTFLG 	= -lfftw3 -lfftw3f
MPIFLG  = -DMPI  
libpaw = libpaw/0_m_libpaw_defs.o
lib_ = lib/0_library_mpi_module.o

#FC = /home/samuelh/app-src/openmpi-1.6/build-lahey/bin/mpif90 -DMPI -lfftw3
#LIBXC = -lxcf90 -lxc
XCID  = XCI/*f90   # XCI includes a few interface files to LIBXC

LIBXC = -I/home/minh/codes/LIBXC/include -L/home/minh/codes/LIBXC/lib -lxcf90 -lxc
#XC=XCID/*.f90 
#LIBXC = ' '

# compile
all :  $(libpaw) $(lib_) main clean

$(libpaw) :
	cd libpaw  && $(CC) -c libpaw_libxc.c && $(FC) $(MPIFLG) -c *.F90 *.f90 && cd ..

$(lib_) :
	cd lib && $(FC) $(MPIFLG) -c *.f90 *.f   && cd ..

main :
	$(FC) -c -I ./libpaw 1_libpaw_mod.f90
	$(FC) -I ./lib ./lib/*.o -I ./libpaw ./libpaw/*.o *.f90 vloc/*.f90 -o opaw_tddft.x \
        $(XCID) $(FFTFLG) -lblas -llapack $(LIBXC) 
       

clean :
	touch empty.o empty.mod 
	rm *.o *mod 
