!interface with libpaw subroutines
!ixc, ecut, ecutpaw will be subject to change in actual implementation
!other flags are fixed
module libpaw_mod 
    use m_pawpsp 
    use m_pawtab
    use m_pawrad
    use m_pawxmlps
    use m_pawang
    use m_paw_init
    use m_pawrhoij
    use m_paw_ij
!    use m_libpaw_libxc_funcs 
    use m_paw_an
    use m_pawfgrtab
    use m_pawdij
    use m_paw_denpot
    use m_paw_finegrid
    use m_paw_nhat
    implicit none

    save
    type(paw_setup_t) :: pawsetup
    type(pawpsp_header_type) :: pawpsp_header
    type(pawang_type) :: pawang
    
    type(pawrad_type),allocatable :: pawrad(:)
    type(pawtab_type),allocatable :: pawtab(:)
    type(pawrhoij_type),allocatable :: pawrhoij(:)
    type(paw_ij_type),allocatable :: paw_ij(:)
    type(paw_an_type),allocatable :: paw_an(:)
    type(pawfgrtab_type),allocatable :: pawfgrtab(:)

    real*8 ,allocatable :: znucl(:)
    integer,allocatable :: typat(:)
    integer,allocatable :: nattyp(:)
    integer,allocatable :: lexexch(:),lpawu(:)
    integer,allocatable :: l_size_atm(:)
    real*8 ,allocatable :: xred(:,:)
    

    character(264)     :: filename

    integer           :: lloc,lmax,pspcod,pspxc
    real*8            :: r2well,zion
    real*8            :: xcccrc

    integer           :: mqgrid_ff=3001,lnmax=6,ipsp=1,ixc,mqgrid_vl=3001
    integer           :: xclevel,usexcnhat=0,xcdev=1,usewvl=0,icoulomb=0
    integer           :: gnt_option=1,lcutdens=10,lmix=10,mpsang,nphi=13,nsym=1,ntheta=12
    integer           :: pawspnorb=0,usepotzero=0,nspden=1,ntypat
    integer           :: cplex=1,nspinor=1,nsppol=1
    integer           :: enunit=0,ipert=0,nfft,nfftot
    integer           :: pawprtvol=0
    integer           :: ngfft(18)
    integer           :: ngfft1(18)
    
    real*8            :: gsqcut_eff,hyb_range_fock=0.0d0
    real*8            :: ecutpaw,denpos=1d-14
    real*8            :: gsqcut,gsqcutdg,dq,qmax
    real*8            :: spnorbscl=1d0,ucvol,charge=0.0d0
    real*8            :: gprimd(3,3)=0d0,qphon(3)=0d0

    !real*8,allocatable:: vxc_r(:,:),vtrial_r(:) !for checking only
    real*8,allocatable:: ffspl(:,:,:),vlspl(:,:,:)
    real*8,allocatable:: qgrid_ff(:),qgrid_vl(:)

    real*8            :: epsatm, scccrc

    real*8            :: gmet(3,3)
    real*8            :: rprimd(3,3)
    real*8            :: spinat(3,1)

    real*8            :: tol10=1d-10

    real*8            :: compch_sph,epaw,epawdc 
    real*8            :: nucdipmom(3,1)
    
    integer           :: ii

end module libpaw_mod      
