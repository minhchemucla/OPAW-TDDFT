module main_mod
    implicit none
    save

    integer :: nx,ny,nz,nn!values are for testing
    real*8  :: box_x,box_y,box_z
    real*8  :: dx,dy,dz,dv
    real*8  :: kkx,kky,kkz
    real*8  :: rnel
    real*8  :: ekcut
    integer :: nocc,nb !occupied states,total #.states
    integer :: scale_vh
    integer :: nscf,iscf
    logical :: periodic
    logical :: read_hminmax=.false.
    logical :: diis_dij=.true.
    logical :: ekread=.false.

    real*8  :: p_fg!parameter for fine grid
    real*8  :: mix_diis
    real*8 :: mix_diis1

    !number of grids; number of bands
    real*8  :: xmax,ymax,zmax

    real*8  :: dmua,mu0
    real*8  :: hmin_0,hmax_0
    real*8, allocatable :: mu(:),hmin(:),hmax(:)

    real*8, allocatable :: kpt(:,:) !3*#. k points
    real*8, allocatable :: wk(:) !#. kpoints
    integer :: nkpt,nk_loc

    integer :: funct=0 !0:lda,1:pbe
    integer :: h_type=0 !0: S^-1Htilde, 1: S^-1/2(Htilde)S^-1/2
    integer :: funct_x,funct_c

    real*8, allocatable :: vk(:)
    real*8, allocatable :: ek(:,:,:,:) !nx,ny,nz,number of k points
    real*8, allocatable :: vloc_tot(:,:,:) !local ionic potential
    real*8, allocatable :: nhat(:,:,:) !compensation charge
    real*8, allocatable :: ncoret(:,:,:) !pseudo core density
    real*8, allocatable :: vxc(:,:),vks(:,:,:),vh(:,:,:),vks_old(:,:,:)
    real*8, allocatable :: dijall(:)
    integer :: siz_dijall
!    real*8, allocatable :: phi(:,:,:,:) !stores occupied orbitals
    complex*16, allocatable :: expk(:,:,:,:) !nx,ny,nz,number of k points
    complex*16, allocatable :: phit(:,:,:,:,:) !phibar=s^1/2 phit
    complex*16, allocatable :: phit_tot(:,:,:,:,:) !phibar=s^1/2 phit
    real*8, allocatable :: dens(:)!,dens_p(:) !the electron density of this and previous iteration
end module main_mod
