# OPAW_TDDFT

The projector augmented wave (PAW) method linearly maps the highly oscillatory "all-electron" wavefunctions to smoother pseudowavefunctions that are cheaper to represent on a real-space grid. However, PAW has non-orthogonal pseudowavefunctions complicating its implementation in theory and code. We have developed orthogonal PAW (OPAW) to remedy this issue and successfully implemented in the real-space plane-wave DFT framework in Wenfei Li and Daniel Neuhausers' [paper](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.102.195118).

Further work was taken by Minh Nguyen, Tim Duong, and Daniel Neuhauser to develop [OPAW-TDDFT](https://pubs.aip.org/aip/jcp/article-abstract/160/14/144101/3281117/Time-dependent-density-functional-theory-with-the?redirectedFrom=fulltext) using 4th order Runge-Kutta for the time propagation. This code was used to generate the data in the OPAW-TDDFT paper and is now available here on GitHub mainly for other Neuhauser members to generate OPAW-DFT wavefunctions. 

For the OPAW-TDDFT, the code will produce a dipole-dipole correlation function: 

$$\Delta n_i(r,t)=\frac{\hat{r}}{\gamma}\Big(n^{\gamma}(r,t)-n^{\gamma=0}(r,t)\Big).$$

where $\hat{r}$ is the direction of perturbation. This function is then used to calculate the absorption spectra of molecules. Refer to the OPAW-DFT and OPAW-TDDFT papers for theory.

The code has the following features

 - gamma-point (# k points = 1).
 - close-shelled (`nspin=1`).
 - LDA and PBE exchange-correlation functionals.
 - Periodic and nonperiodic ([Martyna-Tuckermann](https://pubs.aip.org/aip/jcp/article-abstract/110/6/2810/474725/A-reciprocal-space-based-method-for-treating-long?redirectedFrom=fulltext) method)

The parallelization scheme in this code is over the number of grid points, so the number of MPI processes has to be a factor of it.

# Chapters

  1. [Prerequisites](#prequisites) 
  2. [OPAW Input Files](#input)
  3.  [Compiling](#compiling)  
  4. [How to run](#how_to)
  5. [Acknowledgements](#acknowledgements)

##      <a id="prequisites"></a>     1.   Prerequisites
This software is written in FORTRAN with MPI parallelization. Compile this code using the `mpif90` compiler wrapper or any other  MPI wrapper. There is a `libpaw_libxc.c` C file in the `opawlib/libpaw` directory so a C compiler like `gcc` will be need too.

The following are scientific and mathematical libraries that are necessary. 

  1. **LIBXC**: Library for exchange-correlation functions. More information on how to install LIBXC can be found at this [link](https://libxc.gitlab.io/).

  2. **FFTW3** : Fastest Fourier Transform in the West version 3.x.x.

  3. **BLAS**: Basic Linear Algebra Subprograms.

  4. **LAPACK**: Linear Algebra Package.

The last three are typically available on supercomputers and Linux package tools such as `apt-get`.

##    <a id="input"></a>        2.     OPAW Input Files

The following input files are needed for this code to work. 

- **dim_paw.inp**: Input file explained in the next [chapter](#how_to).
- **counter.inp**: Related to random number generation. Changing the integer will change the random numbers.
- **random.inp**: Random number seed. You can change the numbers in here to also change the seed.
- **pawfiles:** This file should contain a list of elements and the corresponding PAW potential.   For example:

	    H     H.LDA_PW-JTH.xml
	    C     C.LDA_PW-JTH.xml

- **PAW potentials:** The PAW potentials in the pawfiles should be in the same directory as the pawfiles itself.
- **cnt.ini:** Coordinate file of the system in Bohr. An example of Naphthalene is provided in the `example` directory.
- 
##   <a id="compiling"></a> 3. Compiling

The `src` directory has the following subdirectories:
	
 - `vloc`
 - `libpaw`
 - `lib`
 - `XCI`

The `lib` directory contains various subroutines used across Neuhauser group software. The `vloc` directory contains routines to calculate the local part of the PAW potential and `libpaw` contains many ABINIT subroutines. The `XCI` directory contains subroutines to interface with LIBXC. In the compilation of the example program, there are 3 steps:

 1. Compile the ABINIT subroutines in `opawlib/libpaw`.
 2. Compile the Neuhauser library in `opawlib/lib`.
 3. Compile the OPAW-TDDFT code.

In the `makefile` the following variables are defined
	 
	libpaw_ = libpaw/0_m_libpaw_defs.o
	lib_ = lib/0_library_mpi_module.o
	
With the following rules.

	$(libpaw_) :
	  cd libpaw  && $(CC) -c libpaw_libxc.c && $(FC) $(MPIFLG) -c *.F90 *.f90 && cd ..

	$(lib_) :
	  cd lib && $(FC) $(MPIFLG) -c *.f90 *.f   && cd ..

	cleanlib :
	  touch lib/a.o lib/a.mod
	  rm lib/*.o lib/*.mod

	cleanlibpaw :
	  touch libpaw/a.o libpaw/a.mod
	  rm $(libpaw)/*.o $(libpaw)/*.mod

The `libpaw_` and `lib_` variables are there to tell GNU Make to not compile `libpaw` and `lib` if `0_m_libpaw_defs.o` and `0_library_mpi_module.o` are present i.e. if the object files in `libpaw` and `lib` have been compile already. If you make modifications to any files in `libpaw` or `lib`, you will have to run `make cleanlibpaw` and `make cleanlib` and then recompile.

Below here are other make variables defining the MPI Fortran and C compilers, LIBXC, FFTW2, BLAS, and LAPACK  libraries, and some MPI flags.

	FC    = mpif90
	CC    = gcc
	FFTFLG  = -lfftw3 -lfftw3f
	BLASFLG   = -lblas
	LAPACKFLG   = -llapack
	LIBXC = -I/home/minh/codes/LIBXC/include -L/home/minh/codes/LIBXC/lib -lxcf90 -lxc
	libs = $(FFTFLG) $(BLASFLG) $(LAPACKFLG) $(LIBXC)
	MPIFLG  = -DMPI -O3 -g -fcheck=all -fbacktrace #-Wall
	XCI  = XCI/*f90

The `-g`, ` -fcheck=all` and `-fbacktrace` are debugging flags. 

To compile the library with the main program, base your `main` rule off the following example:

	main :
	  $(FC) -c -I ./libpaw 1_libpaw_mod.f90
	  $(FC) -I ./lib ./lib/*.o -I ./libpaw ./libpaw/*.o *.f90 vloc/*.f90 -o opaw_tddft.x \
        $(XCI) $(libs)

Then define the following all rule.

	all :  $(libpaw_) $(lib_) main clean

Simply compile with `make` while in the `src` directory.


##  <a id="how_to"></a>   4. How to run

### dim_paw.inp
This section is to explain the `dim_paw.inp` input file that controls the parameters of the OPAW-DFT calculation and OPAW-TDDFT simulation.

	nx  56         !number of grid points in x direction
	ny  52         !number of grid points in y direction
	nz  32         !number of grid points in z direction
	box_x 28       !Size of box in x direction. dx=box_x/nx
	box_y 26       !Size of box in y direction. dy=box_y/ny
	box_z 16       !Size of box in z direction. dx=box_z/nz
	periodic F     !Use periodic or non-periodic (Martyna-Tuckermann approach)

	dmua 0.1       !Positive shift to mu in energy filtering

	mix_diis 0.15  !DIIS mixing parameter
	mix_diis1 -1   !if < 0 then no diis for dij
	nscf 30        ! number of SCF cycles
	funct 0        !0:pwlda, 1:pbe

	iscf_hminmax 5 !Will update H_max and H_min up to iscf_hminmax. if <0 will not update.

	hmax 70        !Maximum of energy range of Hamiltonian
	hmin -10			 !Minimum of energy range of Hamiltonian

	finegrid 0.15  !nfovnr=grid spacing/finegrid
	ekcut -1d0     !Kinetic energy cutoff. Negative values will have PAW subroutines generate
	h_type 1       !0=s^-1h psuedo wfs output => wf.txt/bin, 1=s^-1/2hs^-1/2 orthogonal wfs => wf_bar.txt/bin
	flg_bin F      ! wavefunction output to have .txt or more compact .bin format.

	tddft 0 			 !-1 just opaw-dft, 0 opaw-dft followed by tddft, 1 read in wfs from wf(_bar).txt/bin then do tddft
	nt  10000			 ! number of time-steps
	dt 0.05        ! length of time-step
	exct_pol 1     ! Direction of perturbation 1,2,3 -> x,y,z
	strength 1e-3  ! gamma
	theory 3       !1=static Ham, 2=RPA, 3=TDDFT

Included in this repository is an Naphthalene example. You can run the code with `mpirun`. For example

	mpirun -n 4 opaw_tddft.x 
 
The code will output

- A `wf.txt`, `wf_bar.txt`, `wf.bin`, or `wf_bar.bin` file depending on `flg_bin` and `h_type` that contains either. The `wf.txt/bin` are the PAW pseudowavefunctions and `wf_bar.txt` are the OPAW pseudowavefunctions. For developing OPAW-sGW and OPAW-sBSE, you should generally use `h_type=1` for get the orthogonal wavefunctions.
- `eig` (`h_type=0`) or `eig_bar` (`h_type=1`) that contains the eigenvalues.
- `dip.dat` - The perturbed and unperturbed portions of the dipole-dipole correlation function in the introduction of this README.
- `d_dip.dat` - The dipole-dipole correlation function.
- `fort.xxx` - Various debugging outputs from the code and ABINIT subroutines. You can safely delete them after the calculation is finished.


##   <a id="acknowledgements"></a> 5.  Acknowledgments

This code is supported by the U.S. Department of Energy, Office of Science, Office of Advanced Scientific Computing Research and Office of Basic Energy Sciences, Scientific Discovery through Advanced Computing (SciDAC) program under Award Number DE-SC0022198.

Much of the PAW capabilities were adapted from the open source [ABINIT](https://www.abinit.org/) software version 8.0.8. The PAW algorithm in ABINIT was based on the work of [Torrent, Marc, et al. "_Implementation of the projector augmented-wave method in the ABINIT code: Application to the study of iron under pressure_" Computational Materials Science **42**, 337-351 (2008)](https://doi.org/10.1016/j.commatsci.2007.07.020).

This code authored by Minh Nguyen based off the OPAW-DFT work done by Wenfei Li and Daniel Neuhauser.
