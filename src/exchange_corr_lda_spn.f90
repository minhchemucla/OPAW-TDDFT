subroutine get_vx_lda(nn, ng, vx)
  implicit none

  integer ng
  real*8 nn(ng), vx(ng)

  vx = -0.738558766382022046 * 4.d0/3.d0* (max(nn,1e-20))**(1.d0/3.d0)
end subroutine get_vx_lda

subroutine get_vcn_spn_lda(nn, ng, nspn, vc)
  implicit none
  integer ng, nspn, ig, is
  real*8 nn(ng, nspn), vc(ng, nspn), n1, n2, dn, nnn, ec0a, de1, de2
  real*8, external :: ec_pw

  select case(nspn)
  case(1)
     do ig=1,ng
        n1 = max(1e-20,nn(ig,1)/2d0 )
        n2 = n1
        nnn = n1+n2
        ec0a  = ec_pw(n1,n2)
        dn = 1e-5 * nnn
        !de1 = (ec_pw(n1+dn,n2   )-ec0a)/dn
        de1 = (ec_pw(n1+dn/2d0,n2+dn/2d0)-ec0a)/dn
        vc(ig,1) = ec0a + nnn*de1
     enddo
  case(2)
      stop 'not using this yet'
!     do ig=1,ng
!        !
!        ! vc(ig,1) = d( e(n1,n2)*(nnn) )/ dn1 = e(n1,n2) + nnn*de(n1,n2)/dn1
!        !
!        n1 = max(1e-20,nn(ig,1))
!        n2 = max(1e-20,nn(ig,2))
!        nnn = n1+n2
!        ec0a  = ec_pw(n1,n2)
!        dn = 1e-5 * nnn
!        de1 = (ec_pw(n1+dn,n2   )-ec0a)/dn
!        de2 = (ec_pw(n1,   n2+dn)-ec0a)/dn
!        vc(ig,1) = ec0a + nnn*de1
!        vc(ig,2) = ec0a + nnn*de2
!     end do
  case default
     stop ' nspn in vcn_spn '
  end select
end subroutine get_vcn_spn_lda

