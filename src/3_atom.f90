module atom_mod
    use paw_mod
    implicit none
    save
    character(len=50) :: atomfile
    integer :: natom, ntype

    integer, allocatable :: atom_z(:), diff_z(:) !z of atoms, different zs, from small to large
    integer, allocatable :: atom_map(:) !type of each atom
    integer, allocatable :: natom_type(:) !how many atoms of each type
    integer, allocatable :: atindx1(:),atindx(:) !used for libpaw
    integer, allocatable :: ngrid(:) !count how many grid points in local grid

!    type atom_old
!        real*8, allocatable :: rhoij
!    end type atom_old
    type atom
        ! wenfei's stuff for dft ------------------------------------------------------
        real*8  :: coord(3) !location
        
        real*8,  allocatable :: rhoij(:,:) !ca*ca
        real*8,  allocatable :: dij(:,:) 

        real*8,  allocatable :: local_g3d(:,:) !shapefunction*spherical harmonics
        integer, allocatable :: local_grid(:,:)
        !ngrid*3, stores ix,iy,iz of each point in local grid

        integer :: ir_start(3) !starting point of rough grid
        logical :: edge=.false. !for periodic

        real*8,  allocatable :: local_p3d(:,:,:,:) !projector functions
        complex*16,  allocatable :: local_p3d_c(:,:,:,:,:)

        real*8,  allocatable :: local_aep3d(:,:,:,:) !all-electron partial wfs
        complex*16,  allocatable :: local_aep3d_c(:,:,:,:,:)

        real*8,  allocatable :: local_pp3d(:,:,:,:) !pseudo partial wfs
        complex*16,  allocatable :: local_pp3d_c(:,:,:,:,:)

        real*8,  allocatable :: local_p3d1(:,:,:,:) !orthogonal projector functions
        complex*16,  allocatable :: local_p3d1_c(:,:,:,:,:)

        real*8,  allocatable :: s(:,:),sinv(:,:),ssqinv(:,:) !overlap matrices

        real*8  :: zval
            
        real*8,  allocatable :: dijfockhat(:,:,:,:) !dij(mu,nu) = -sum_LM int([v_x(r)]_(mu,nu) Q_ij^LM,a(r))
                                    ![v_x(r)]_(mu,nu) = int((n^~_mu,nu(r')+nhat_mu,nu(r')/|r-r'|)
        real*8,  allocatable :: dijfock(:,:) !dij(mu,nu) = -sum_LM int([v_x(r)]_(mu,nu) Q_ij^LM,a(r))
                                    ![v_x(r)]_(mu,nu) = int((n^~_mu,nu(r')+nhat_mu,nu(r')/|r-r'|)
        real*8,  allocatable :: dijfock_vv(:,:), dijfock_cv(:,:) 
!tddft stuff added by minh ------------------------------------------------------------

        real*8,  allocatable :: rhoij_pert(:,:) !perturbation stuff
        real*8,  allocatable :: dij_pert(:,:) !perturbation stuff
        real*8,  allocatable :: dij_pert_old(:,:) !perturbation stuff
        real*8,  allocatable :: dij_old(:,:) 
        real*8,  allocatable :: rho0ij(:,:) 

   !tbar variables and arrays-----------------------------------------------------------
        !!tbar = s^{-1/2}ts^{-1/2} = t + sum_i |b_i>w_i<b_i|
        !real*8, allocatable  :: local_b(:,:,:,:) 
        !complex*16, allocatable  :: local_b_c(:,:,:,:)
        !complex*16, allocatable  :: w(:) 

        !tbar_2 = s^-1 t = t + sum_i |local_p3d_i> y_i <local_p3d_i| t 
        !       = t + sum_i |local_p3d_i> y_i <c_i| 
        !  |c_i> = t|local_p3d_i>  
        !  <c_i| = (t|local_p3d_i>)^dagger  !we can do this since t is hermitian
        real*8, allocatable :: local_c(:,:,:,:) !|c_i>
        complex*16, allocatable :: local_c_c(:,:,:,:) ! complex |c_i>
        complex*16, allocatable :: local_c_c_test(:,:,:,:) ! complex |c_i>
        !real*8, allocatable :: mat_transform_c(:,:),vec_transform_c(:,:)

        ! c1 = a^t c is the biorthogonal basis set to p3d1_c
        ! <p3d1_i|c1_j> = delta_ij
        ! thus sum_i |local_p3d_i> y_i <c_i| = sum_i |local_p3d_i> y_iw^dagger_ij<c1_j|
        ! w(i,j) = <p3d_i|c_j>, a = inv(w)
        complex*16, allocatable :: wmat(:,:), amat(:,:)
        complex*16, allocatable :: local_c1_c(:,:,:,:)
        !taylor expansion for c_ij= sinv_iw^dagger_ij 
        complex*16, allocatable :: c_taylor(:,:)
        complex*16, allocatable :: c_taylor_terms(:,:,:)
        complex*16, allocatable :: cmat(:,:)
        complex*16, allocatable :: cbarmat(:,:)

   !vbar variables and arrays-----------------------------------------------------------
   !Note to Minh from Minh: Please change the variable names later to be more clear
        real*8, allocatable :: local_cv(:,:,:,:) !|c_i>
        complex*16, allocatable :: local_cv_c(:,:,:,:) ! complex |c_i>
        complex*16, allocatable :: local_cv_c_test(:,:,:,:) ! complex |c_i>
        complex*16, allocatable :: wvmat(:,:), avmat(:,:)
        complex*16, allocatable :: local_cv1_c(:,:,:,:)
        complex*16, allocatable :: cv_taylor(:,:)
        complex*16, allocatable :: cv_taylor_terms(:,:,:)
        complex*16, allocatable :: cvmat(:,:)
        complex*16, allocatable :: cvbarmat(:,:)

    end type atom

    type(atom), allocatable :: atominfo(:)!,atominfo_old(:)
    type(paw_element), allocatable :: pawinfo(:)
end module atom_mod      
