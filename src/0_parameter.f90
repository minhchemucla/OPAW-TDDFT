module param_mod
      implicit none
      save

      logical :: k_flg = .false.
      integer :: nch=10
      integer :: niter=25
      integer :: ndi=4
      integer :: ndi1=4
      integer :: iscf_hminmax
      integer :: iscf_start

      real*8  :: fctr=1.2d0
      real*8  :: rpad=1.0d0 !for spline
      real*8  :: rpad_r=2d0
      real*8  :: ek_factor=1d0
      real*8  :: param_dij_max=500d0
      real*8  :: tollsij=0.01

      real*8  :: rnel_1=1.3
      integer :: rnel_2=5,rnel_3=20
      !nb=max(ceiling(nocc*r_1)+r_2,r_3)
        
      logical :: flg_bin=.false.
end module param_mod      
