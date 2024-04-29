module paw_mod
    implicit none
    save
    !for chebyshev    
    integer :: nfovnr ! change to 4 or 6, DN    !for spline !nfine over nrough
    ! tells you how dense fine grid is (2 means twice as dense, etc)

    type paw_element !datatype for storing paw information for each element
!general info    
        real*8  :: Zat, core, val !Z, core electron, valence electron
        character(len=2 ) :: symbol !atomic symbol
        character(len=20) :: filename !name of xml file
    
! 1 dimensional
        real*8  :: rcut !radius of augmentatin sphere
        integer :: nrcut !radial grid corresponding to rcut
        integer :: nr_int !carry out integration up to this point nr_int > nrcut

        integer :: nr, nstates !# of radial grid pts, number of states
        real*8,  allocatable :: rr(:), dr(:) 
        !phi is the 1d radial all electron partial waves. \phi_i^(a) in the paper
        !phitilde are the 1d radial psuedo partial waves. \tilde\phi_i^(a)  in the paper
        !ptilde are the 1d radial projector functions 
        ! They are all on the rough grid
        real*8,  allocatable :: phi(:,:),phitilde(:,:),ptilde(:,:) 
        !phi1, phitilde1, ptilde1 are transformed orthogonal version (rough grid)
        real*8,  allocatable :: phi1(:,:),phitilde1(:,:),ptilde1(:,:) 
        real*8,  allocatable :: y2a_p(:,:) !used for spline of radial part
        real*8,  allocatable :: y2a_p1(:,:) !used for spline of radial part
        real*8,  allocatable :: vloc(:), ncoretilde(:) 
        real*8,  allocatable :: gl(:,:), qijlm(:,:,:) !M,mi,mj
        !gl = shape function; qijlm Eq (14)
        real*8,  allocatable :: sij(:,:)
        real*8,  allocatable :: ssqinvij(:)
        real*8,  allocatable :: occ_l(:) 
        !read from xml, radial grid, derivative, projector, ae wf, pseudo wf,
        !local ionic potential, pseudo core density
        !shape function, qijkl: equation (14)
!        real*8,  allocatable :: ca(:) !coefficients,for testing 1d only

        integer :: nl !how many different l
        
        integer, allocatable :: lstate(:) !l of each state
        integer, allocatable :: lstate_diff(:) !what are the different l's
       
        integer :: mstates !The total number of states. For each l state, there will be 2l+1 m states
        integer, allocatable :: mstate(:) !m of each state
        integer, allocatable :: ms_ls(:)!each 3d state correspond to which 1d state 
                                        ! for example 2 s states and 2 p states would be
                                        ! ms_ls=(\1,2,3,3,3,4,4,4\)

        real*8,  allocatable :: mat_t(:,:)! for checking only. Initially proportional sij matrix Eq. (4)
        real*8,  allocatable :: mat_sp(:,:)! for checking only. Initially the overlap matrix L_{ij}^(a) under Eq. (A.1)

        integer :: flg_comp,nrcomp !compensation charge type, integer
        real*8  :: rcomp !compensation charge radius
        
        integer :: nrough(3),nfine(3) !number of grid points for rough/fine grid
        real*8  :: dxf(3) !fine grid spacing 
        real*8,  allocatable :: bx(:,:),by(:,:),bz(:,:) ! Ono-Hirose method- linear projection matrices
        !the b matrices also include the dxf/dx, dyf/dy, and dzf/dz elements
        !the indices are reverse of the paper i.e b(i,i_f) instead of b(i_f,i)


        ! TDDFT Additions - Minh
        real*8,  allocatable :: ek_fine(:,:,:,:) !fine grid kinetic energy
        ! tptilde1 = T ptilde1
        real*8,  allocatable :: tptilde1(:,:)
    end type paw_element

! allocate after counting number of different elements in molecule    
    
! other stuff
    logical :: file_exists
    integer :: stat
    integer :: whereis

end module      

