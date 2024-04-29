## compiler options
FC   	= mpif90 
CC		= gcc 
FFTFLG 	= -lfftw3 -lfftw3f
BLASFLG   = -lblas
LAPACKFLG   = -llapack
LIBXC = -I/home/minh/codes/LIBXC/include -L/home/minh/codes/LIBXC/lib -lxcf90 -lxc
libs = $(FFTFLG) $(BLASFLG) $(LAPACKFLG) $(LIBXC)
MPIFLG  = -DMPI -O3 -g -fcheck=all -fbacktrace #-Wall  
XCI  = XCI/*f90  

libpaw_ = libpaw/0_m_libpaw_defs.o
lib_ = lib/0_library_mpi_module.o

# compile
all :  $(libpaw_) $(lib_) main clean

$(libpaw_) :
	cd libpaw  && $(CC) -c libpaw_libxc.c && $(FC) $(MPIFLG) -c *.F90 *.f90 && cd ..

$(lib_) :
	cd lib && $(FC) $(MPIFLG) -c *.f90 *.f   && cd ..

main :
	$(FC) -c -I ./libpaw 1_libpaw_mod.f90
	$(FC) -I ./lib ./lib/*.o -I ./libpaw ./libpaw/*.o *.f90 vloc/*.f90 -o opaw_tddft.x \
        $(XCI) $(libs) 

cleanlib :
	touch lib/a.o lib/a.mod
	rm lib/*.o lib/*.mod
          
cleanlibpaw :
	touch libpaw/a.o libpaw/a.mod
	rm $(libpaw)/*.o $(libpaw)/*.mod
       
clean :
	touch empty.o empty.mod 
	rm *.o *mod 
