SUBROUTINE LPN(N,X,PN,PD)
  !
  !       ===============================================
  !       Purpose: Compute Legendre polynomials Pn(x)
  !                and their derivatives Pn'(x)
  !       Input :  x --- Argument of Pn(x)
  !                n --- Degree of Pn(x) ( n = 0,1,...)
  !       Output:  PN(n) --- Pn(x)
  !                PD(n) --- Pn'(x)
  !       ===============================================
  !
  IMPLICIT DOUBLE PRECISION (P,X)
  DIMENSION PN(0:N),PD(0:N)
  PN(0)=1.0D0
  PN(1)=X
  PD(0)=0.0D0
  PD(1)=1.0D0
  P0=1.0D0
  P1=X
  DO K=2,N
     PF=(2.0D0*K-1.0D0)/K*X*P1-(K-1.0D0)/K*P0
     PN(K)=PF
     IF (DABS(X).EQ.1.0D0) THEN
        PD(K)=0.5D0*X**(K+1)*K*(K+1.0D0)
     ELSE
        PD(K)=K*(P1-X*PF)/(1.0D0-X*X)
     ENDIF
     P0=P1
     P1=PF           
  END DO
end SUBROUTINE LPN
