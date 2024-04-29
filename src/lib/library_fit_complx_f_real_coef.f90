!
!  Given fexp(j) [= (eps-1)*w] = complex (exp. found) function,
!
!  Let's fit to N terms, p=1,...,N 
!
!  f(j) = sum_p cp g(ap,bp,w)= sum_p cp gp(w).  The form of g will be given later, but we want to be general.
!
!  For our case, cp is real and gp is complex.
!  Also, the n'th term is what we call +i*sigma  (times w/w =1). It can be a part of these terms, with g(p=N)= +i
!
!  Let's first find the best c(:) set.
!
!  Ob=sum_j h(j) |f(j)-fexp(j)|^2 + eta* sum_p cp^2      
!
!  {A side note: the weight h(j) could be simply dw, 
!  and the way to do it is : for j=1 take 0.5(w(2)-w(1)), for j=Nj take 0.5*(w(Nj)-w(nj-1)), for the rest take
!  h(j)=(w(j+1)-w(j-1))/2.d0 1<j<Nj .}
!
!  Ob = sum(h |f|^2)+ sum(h fexp f* ) + sum(h fexp* f) + sum(|fexp|^2) + eta sum(cp^2)  -- cp real, remember!
!
!  Ob = sum_wpq (h(w) c(p) c(q) gp*(w) gq(w)) + sum_wp (h(w) fexp*(w) cp gp(w) ) + c.c.of 2nd term 
!    + sum(|fexp|^2) 
!    + eta sum_p cp^2
!
!  Ob = sum_pq M(p,q) c(p) c(q) - sum_p L*(p) c(p) + c.c.of 2nd term +   
!    + sum(|fexp|^2) 
!    + eta sum_p cp^2
!
!  where  
!    sum_pq M(p,q) = sum_pq Mprimitive(p,q) c(p) c(q),  and 
!  
!  Mprimitive(p,q)= sum (h(w) gp*(w) gq(w))
!
!
!  and since Mprimitve is hermitian, and the c's are real, we can replace cT Mprimitive c by cT M c, where
!  M is the real part of Mprimitive, which is a symmetric matrix:
!
!           M = real(Mprimiitive)
!
!  Also
!      Lp* = sum_w (h(w) fexp*(w) gp(w))  
!  so
!      Lp  = sum_w h(w) fexp(w) gp*(w)
!
!
!
!
!  So let's diff.:
!
!  0 = dOb_dc = 2Mc - L* - L + 2eta c = 0
!   Mc + Re L + eta c = 0
! 
!   c = (M+eta*Unity)^-1  ReL
! 
!


  subroutine Linear_fit_complex_func_with_real_coef(nw,np,ww,h,fexp,gsub,anlp,nlp,cc,eta,obj)
    use mat_module
    implicit none
    integer        nw, np, nlp   ! nlp: # of nonlinear parameters for each "ip".  
    real*8         ww(nw),  h(nw)
    real*8         anlp(nlp,np), cc(np) ! cc is output
    complex*16     fexp(nw)
    external gsub  
    complex*16     g(np)
    real*8         eta

    real*8   Lr(np)
    real*8   M(    np, np)
    real*8   Minv( np, np)
    real*8   Obj

    integer  ip, iq, iw

    M = 0.d0
    Lr = 0.d0

    do iw=1, nw
       call gsub(ww(iw), anlp,nlp, np, g )
       
       do ip=1,np
          do iq=1,np
             M(ip,iq) = M(ip,iq) + dble(h(iw)* conjg(g(ip)) * g(iq))
          enddo
       enddo
       
       Lr(:) = Lr(:) + dble(h(iw) * fexp(iw) * conjg(g(:)))
    enddo
    
    ! begin check

    !!!! call check_M

    ! end check

    do ip=1,np
       M(ip,ip) = M(ip,ip)+eta
    enddo

    call mat_inv(M, Minv)

    cc =  matmul(Minv, Lr)

    Obj = Obj_all()
    
  contains
    real*8 function Obj_all()
      implicit none
      real*8 term
      complex*16 ff

      term  = 0.d0
      do iw=1, nw
         call gsub(ww(iw), anlp, nlp, np, g )
         ff = sum(cc*g)
         term = term+ h(iw)* abs((ff-fexp(iw))**2)
      enddo

      Obj_all = term

    end function Obj_all
  end subroutine Linear_fit_complex_func_with_real_coef


    
    
    
