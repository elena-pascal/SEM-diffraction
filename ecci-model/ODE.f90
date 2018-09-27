MODULE odes
   !USE betavalue
   USE doublep
   IMPLICIT NONE
   complex(kind=dp), parameter :: im = CMPLX(0.0_dp,1.0_dp)
   real(kind=dp), parameter :: PI = 3.14159265358979d0
   CONTAINS

     SUBROUTINE FEX (NEQ, T, Y, YDOT, RPAR, IPAR)
     !-------------------------------------------------------------------
     !The ODE function  YDOT=dY/dT=F(T,Y)
     !                  INPUT  : T, Y
     !                  OUTPUT : YDOT
     !RPAR is an array for passing real parameters (9 elements):
     !                  RPAR(1) = w 
     !                  RPAR(2) = Xsig/Xsig'
     !                  RPAR(3) = Xsi0'/Xsig'
     !                  RPAR(4) = gb 
     !                  RPAR(5) = nu
     !                  RPAR(6) = thetaB
     !                  RPAR(7) = X   !current x value
     !                  RPAR(8) = Y   !current Y value
     !                  RPAR(9) = thickness (=-TEND)
     !IPAR is an array for passing integer parameters
     !-------------------------------------------------------------------
         USE doublep
         IMPLICIT NONE
         integer,                          intent(in) :: NEQ
         real(kind=dp),                    intent(in) :: T
         real(kind=dp),                    intent(in) :: RPAR(9) !(:) doesnt do what I expect it to do !!!!!
         integer,                          intent(in) :: IPAR
         complex(kind=dp), dimension(NEQ), intent(in) :: Y
         complex(kind=dp), dimension(NEQ), intent(out):: YDOT
         real(kind=dp)   :: beta
         !-------no absorbtion----------------------------------
         ! YDOT(1) = PI*im*Y(2)
         ! YDOT(2) = PI*im*Y(1) + (2.0_dp*PI*im*RPAR*Y(2))
         
         !------- with absorbtion and imperfections-------------

         !call beta_ECCI(RPAR(7), RPAR(8), T, RPAR(5), RPAR(6), RPAR(4), beta)
         
         !beta_ECCI(x, y, z, betaH)
		 call beta_ECCI(RPAR(7), RPAR(8), T*30., beta)
         
         YDOT(1) = -PI*Y(1)*RPAR(3)/RPAR(2) + PI*Y(2)*(im-RPAR(3))
         YDOT(2) = PI*Y(1)*(im-RPAR(3)) + PI*Y(2)*(2.0_dp*im*RPAR(1)-(RPAR(3)/RPAR(2))+2.0_dp*im*beta)
     END SUBROUTINE FEX

  
  
  
     SUBROUTINE JEX (NEQ, T, Y, ML, MU, PD, NRPD, RPAR, IPAR)
     !-------------------------------------------------------------------
     !The Jacobian of the system df/dy that I can't get rid of
     !               INPUT  : T, Y
     !               OUTPUT : PD
     !ML, MU lower and upper band width for Jacobian
     !RPAR is an array for passing real parameters (9 elements):
     !                  RPAR(1) = w 
     !                  RPAR(2) = Xsig/Xsig'
     !                  RPAR(3) = Xsi0'/Xsig'
     !                  RPAR(4) = gb 
     !                  RPAR(5) = nu
     !                  RPAR(6) = thetaB
     !                  RPAR(7) = X   !current x value
     !                  RPAR(8) = Y   !current Y value 
     !                  RPAR(9) = thickness (=TEND)
     !IPAR is an array for passing integer parameters
     !-------------------------------------------------------------------
         USE doublep
         IMPLICIT NONE
         integer,                                intent(in) :: NEQ, NRPD
         integer,                                intent(in) :: ML, MU
         real(kind=dp),                          intent(in) :: T
         real(kind=dp),                          intent(in) :: RPAR(9)
         integer,                                intent(in) :: IPAR
         complex(kind=dp), dimension(NEQ),       intent(in) :: Y
         complex(kind=dp), dimension(NRPD, NEQ), intent(out):: PD
         real(kind=dp)   :: beta
         !------no absorbtion--------------------------------
         ! PD(1,2) = PI*im
         ! PD(2,1) = PI*im
         ! PD(2,2) = 2.0_dp*PI*im*RPAR
     
         !------- with absorbtion and imperfections-------------
         

         !call beta_ECCI(RPAR(7), RPAR(8), T, RPAR(5), RPAR(6), RPAR(4), beta)
         
         !beta_ECCI(x, y, z, betaH)
		 call beta_ECCI(RPAR(7), RPAR(8), T*30., beta)

                 
         PD(1,1) = -PI*RPAR(3)/RPAR(2) 
         PD(1,2) = PI*(im-RPAR(3))
         PD(2,1) = PI*(im-RPAR(3))
         PD(2,2) = 2*PI*im*RPAR(1) - PI*RPAR(3)/RPAR(2) + 2.0_dp*im*beta
     END SUBROUTINE JEX
  
END MODULE odes
