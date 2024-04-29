C
C PAGE 326-328: NUMERICAL MATHEMATICS AND COMPUTING, CHENEY/KINCAID, 1985
C
C FILE: RK45AD.FOR
C
C ADAPTIVE SCHEME BASED ON RUNGE-KUTTA-FEHLBERG METHOD (RK45AD,RK45,F)
C
      subroutine crk45ad_driver
      implicit none

      EXTERNAL Scfsimp
      integer, parameter ::nn=2
      real*8 ta, h, tb, emin,emax,hmin,hmax,tdiff
      complex*16 xa(nn)
      integer itmax, iflag
      complex*16, parameter :: ci = (0.d0,1.d0)
      
      TA = 1.0d0
      XA = 2.0d0
      H = 7.8125d-3 
      TB = 1.5625d0 

      tdiff = tb-ta


      ITMAX = 100 
      EMIN = 1.0d-12
      EMAX = 1.0d-10
      HMIN = 1.0d-6 
      HMAX = 1.0d0
      CALL CRK45AD(Scfsimp,TA,XA,H,TB,ITMAX,EMIN,EMAX,HMIN,HMAX,
     $     IFLAG,nn)
      PRINT 3,H,TA,XA(1:min(nn,5))
 3     FORMAT(5X,E10.3,2(5X,E22.14))
      PRINT 5,IFLAG
 5     FORMAT('0 IFLAG =',I2)
       
      print *,' xa_start* exp ',2.d0*exp(-ci*(tdiff))
      STOP
      END 
  
      subroutine Scfsimp(F,T,X, nn) 
      implicit none
      integer nn
      complex*16 x(nn), f(nn)
      complex*16, parameter :: ci=(0.d0,1.d0)
      real*8 t
!      F = 2.d0+(X-T-1.d0)**2  
      f = -ci*X
         
      RETURN
      END 
  
      SUBROUTINE CRK45AD(SF,T,X,H,TB,ITMAX,EMIN,EMAX,HMIN,HMAX,IFLAG, 
     $     nn) 
      implicit none
      integer nn, itmax,iflag, nstep
      real*8 emin,emax, hmin,hmax, t, tb,h,  tsave,
     $     epsi, est, dt
      complex*16 x(nn), xsave(nn)


      data epsi/1d-6/
      external SF
      IFLAG = 2
      NSTEP = 0   
!      PRINT 4     
!      PRINT 3,H,T,X(1:min(nn,5)) 
 2      IF(ABS(H) .LT. HMIN)  H = HMIN
      IF(ABS(H) .GT. HMAX)  H = HMAX
      DT = ABS(TB - T)      
      IF(DT .LE. ABS(H))  THEN
        IFLAG = 0 
        IF(DT .LE. EPSI*max(ABS(TB),ABS(T)))  RETURN    
        H = SIGN(DT,H)      
      END IF      
      XSAVE = X   
      TSAVE = T   
       CALL CRK45(SF,T,X,H,EST,nn)
      NSTEP = NSTEP + 1     
!      PRINT 3,H,T,X(1:min(nn,5)),EST     
      IF(IFLAG .EQ. 0)  RETURN
      IF(NSTEP .GT. ITMAX)  THEN
        IFLAG = 1 
        RETURN
      END IF      
      IF((EMIN .LE. EST) .AND. (EST .LE. EMAX))  GO TO 2  
      IF(EST .LT. EMIN)  THEN 
        H = 2.d0*H 
      ELSE
        H = 0.5d0*H 
      END IF      
      X = XSAVE   
      T = TSAVE   
      GO TO 2 
 3      FORMAT(5X,E10.3,2(5X,E22.14),5X,E10.3)    
 4        FORMAT(9X,'H',22X,'T',26X,'X',18X,'EST')  
      END 
      SUBROUTINE CRK45(SF,T,X,H,EST,nn)    
      implicit none
      real*8 C21,C31,C32, C41,C42,C43, C51,C52,C53,C54,    
     A      C61,C62,C63,C64,C65, A1,A3,A4,A5, B1,B3,B4,B5, B6, C40  

      integer nn
      external SF

      real*8 est,h,t
      complex*16 x(nn),f1(nn), f2(nn), f3(nn), f4(nn),f5(nn),
     $     x5(nn),f6(nn)

      DATA  C21,C31,C32, C41,C42,C43, C51,C52,C53,C54
     B  /0.25d0,0.09375d0,0.28125d0,
     C  0.87938097405553d0,-3.2771961766045d0,3.3208921256258d0,
     D  2.0324074074074d0,-8.0d0,7.1734892787524d0,-0.20589668615984d0/
      
      DATA   C61,C62,C63,C64,C65, A1,A3,A4,A5, B1,B3,B4,B5, B6, C40 
     E / -0.2962962962963d0,2.0d0,-1.3816764132554d0, 
     $     0.45297270955166d0,-0.275d0,
     F  0.11574074074074d0,0.54892787524366d0,0.5353313840156d0,-0.2d0,     
     G  0.11851851851852d0,0.51898635477583d0,0.50613149034201d0,
     $     -0.18d0,   
     H  0.036363636363636d0, 0.92307692307692d0/    
      call SF(F1,T,         X,nn)
      call SF(F2,T+ 0.25d0* H,X + C21*H*F1, nn)  
      call SF(F3,T+ 0.375d0*H,X + C31*H*F1 + C32*H*F2, nn)   
      call SF(F4,T+C40*H,   X + C41*H*F1 + C42*H*F2 + 
     $                          C43*H*F3,nn)  
      call SF(F5,T+H    ,   X + C51*H*F1 + C52*H*F2 + 
     $                          C53*H*F3 + C54*H*F4,nn)  
      call SF(F6,T+0.5d0*H,   X + C61*H*F1 + C62*H*F2 + 
     $                          C63*H*F3 + C64*H*F4 
     $                                   + C65*H*F5, nn)
      X5 = X + B1*H*F1 + B3*H*F3 + B4*H*F4 + B5*H*F5 + B6*H*F6      
      X  = X + A1*H*F1 + A3*H*F3 + A4*H*F4 + A5*H*F5    
      T = T + H   
      EST = sum(ABS(X - X5))
      RETURN      
      END 
