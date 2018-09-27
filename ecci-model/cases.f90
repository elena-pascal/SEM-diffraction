MODULE cases
     !USE betavalue
     USE odes
     IMPLICIT NONE
   
     integer, parameter :: NEQ=2
!----zvode parameters-------------------------------------
     integer, parameter       :: LZW = 24, LRW = 22, LIW = 32 !array legths
     integer, parameter       :: ITOL = 1, MF = 21, ITASK = 1
     integer, parameter       :: IOPT = 1                     !optional input? yes
     real(kind=dp), parameter :: RTOL = 1.E-9_dp, ATOL = 1.E-8_dp 
     !
     integer          :: IPAR
     integer          :: IWORK(LIW) 
     real(kind=dp)    :: RWORK(LRW)
     real(kind=dp)    :: RPAR(3)              !keeps w, Xsig/Xsig' and Xsi0'/Xsig' 
     complex(kind=dp) :: ZWORK(LZW)
!--------------------------------------------------------
    
   CONTAINS
   
!--------interface for subroutine ZVODE called from zvode.f---
!   interface
!      subroutine ZVODE(F,NEQ,Y,T,TOUT,ITOL,RTOL,ATOL,ITASK,ISTATE,IOPT, &
!          & ZWORK,LZW,RWORK,LRW,IWORK,LIW,JAC,MF,RPAR,IPAR)
!          use double
!          external F, JAC
!          integer,          intent(in)   :: NEQ !number of first order ODEs of the system 
!          integer,          intent(in)   :: ITASK !index specifying the task to be performed 
!          integer,          intent(in)   :: IOPT !flag to specfy if optional input is being used
!          integer,          intent(in)   :: MF !flag for Jacobian strategy
!          integer,          intent(in)   :: ITOL !type of error control
!          integer,          intent(in)   :: LZW, LRW, LIW !length of arrays: ZWORK, RWORK, ZWORK
!          integer,          intent(inout):: ISTATE !index to specify the state of the calculation 
!          integer, optional              :: IWORK(LIW) !integer working array for opt params
!          integer, optional              :: IPAR !array of user specified ineger parameters
!          real(kind=dp),    intent(inout):: T !the independent variable
!          real(kind=dp),    intent(in)   :: TOUT !the next value of t
!          real(kind=dp),    intent(in)   :: RTOL, ATOL !relative and abs error tolerance parameters
!          real(kind=dp), optional        :: RWORK(LRW) !real woring array for opt params
!          complex(kind=dp), intent(inout):: Y(NEQ) !the vector of dependent variables
!          complex(kind=dp), optional     :: ZWORK(LZW) !complex working array
!          complex(kind=dp), optional     :: RPAR(:) !array of user specified real/complex parameters
!      end subroutine ZVODE
!   end interface
!-------------------------------------------------------------

    
     SUBROUTINE fringes(Y, T, TEND, param, NOUT, ISTATE, OUTFILE)
         IMPLICIT NONE
         complex(kind=dp), intent(inout) :: Y(NEQ)
         real(kind=dp),    intent(in)    :: T, TEND
         real(kind=dp),    intent(in)    :: param(3)
         integer,          intent(in)    :: NOUT
         integer,          intent(inout) :: ISTATE
         character(30),    intent(in)    :: OUTFILE(:)
         
         integer       :: IOUT
         real(kind=dp) :: TPRINT(NOUT)
         real(kind=dp) :: DTOUT, TOUT           
         real(kind=dp) :: Int0(NOUT), Intg(NOUT)
         real(kind=dp) :: Realg(NOUT), Img(NOUT)
         
         IWORK(6) = 1000 !maximum number of steps allowed during one call
         DTOUT = (TEND-T)/(NOUT-1)   !T step size    
         TOUT = DTOUT           
         RPAR = param        
         
         MainLoopT: DO IOUT = 1,NOUT 
         !start or continue integration to new TEND, T->TEND
             CALL ZVODE(FEX,NEQ,Y,T,TOUT,ITOL,RTOL,ATOL,ITASK,ISTATE,IOPT, &
                   & ZWORK,LZW,RWORK,LRW,IWORK,LIW,JEX,MF,RPAR,IPAR)
     
             TPRINT(IOUT) = T
           
             Realg(IOUT) = REAL(Y(2))              !calculate amplitude and
             Img(IOUT) = IMAG(Y(2))                !phase for phig

             Int0(IOUT) = ABS(Y(1))**2             !
             Intg(IOUT) = ABS(Y(2))**2             !calculate intensities

             IF (ISTATE .LT. 0) THEN               !if integration fails exit
                 EXIT MainLoopT
             ENDIF
    
             TOUT = TOUT + DTOUT
         END DO MainLoopT
         
         !write output file(s)
         open(80, file = OUTFILE(1))
         open(120, file = OUTFILE(2))
         write(80,*) " T ", "                 Int0 ", "               Intg "
         write(120,*) " T ", "     Amplitude(realG) ", "      Phase(imagG) "
         
         DO IOUT=1, NOUT
             write(80,*) TPRINT(IOUT), Int0(IOUT), Intg(IOUT)
             write(120,*) TPRINT(IOUT), Realg(IOUT), Img(IOUT)
         END DO
     END SUBROUTINE


     SUBROUTINE bands(Y, T0, TEND, param, WEND, NOUT, ISTATE, OUTFILE)
         IMPLICIT NONE
         complex(kind=dp), intent(inout) :: Y(NEQ)
         real(kind=dp),    intent(in)    :: T0, TEND
         real(kind=dp),    intent(in)    :: param(3)
         real(kind=dp),    intent(in)    :: WEND
         integer,          intent(in)    :: NOUT
         integer,          intent(inout) :: ISTATE
         character(30),    intent(in)    :: OUTFILE(:)
         
         integer       :: IOUT
         real(kind=dp) :: WPRINT(NOUT), RPAR(3)
         real(kind=dp) :: W, DW, WOUT, TLOOP    !working values      
         real(kind=dp) :: Int0(NOUT), Intg(NOUT)


         IWORK(6) = 3000 !maximum number of steps allowed during one call
         W = param(1) 
         DW = (WEND-W)/(NOUT-1)   !W step size                   
         WOUT = W
         
         MainLoopW: DO IOUT = 1,NOUT 
             TLOOP = T0
             RPAR = (/WOUT, param(2), param(3)/)
             ISTATE = 1                            !start from the begining
             Y(1) = DCMPLX(1.0_dp, 0.0_dp)         !reininitialise Y
             Y(2) = DCMPLX(0.0_dp, 0.0_dp)         !every loop

             CALL ZVODE(FEX,NEQ,Y,TLOOP,TEND,ITOL,RTOL,ATOL,ITASK,ISTATE,IOPT, &
              & ZWORK,LZW,RWORK,LRW,IWORK,LIW,JEX,MF,RPAR,IPAR)
         
             WPrint(IOUT) = WOUT

             Int0(IOUT) = ABS(Y(1))**2
             Intg(IOUT) = ABS(Y(2))**2  
     
             IF (ISTATE .LT. 0) THEN
                 !EXIT MainLoopW
                 RETURN
             ENDIF
             
             WOUT = WOUT + dW
         END DO MainLoopW
         
         !write output file(s)
         open(90, file = OUTFILE(1))
         !write(90,*) " W ", "                      Int0 ", "                      Intg "
         
         DO IOUT=1, NOUT
             write(90,*) WPRINT(IOUT), Int0(IOUT), Intg(IOUT)
         END DO

     END SUBROUTINE
     
!============================== main routine ===========================================     
     SUBROUTINE contrast(phi, T0, TEND, param, Nmesh, SizeMesh, ISTATE, OUTFILE)
         IMPLICIT NONE
         integer,          intent(in)    :: Nmesh(2)
         integer,          intent(inout) :: ISTATE
         real(kind=dp),    intent(in)    :: T0, TEND
         real(kind=dp),    intent(in)    :: param(6)
         real(kind=dp),    intent(in)    :: SizeMesh(2)
         complex(kind=dp), intent(inout) :: phi(NEQ)
         character(30),    intent(in)    :: OUTFILE(:)
         !local variables
         integer       :: i, j, counter, counterw, t
         integer       :: failureGridX(3), failureGridY(3)
         real(kind=dp) :: Xmesh, Ymesh, TRUN
         real(kind=dp) :: XOUT(Nmesh(1)), YOUT(Nmesh(2))
         real(kind=dp) :: dX, dY, dt
         real(kind=dp) :: beta, Average
         real(kind=dp) :: RPAR(9)
         real(kind=dp) :: Int0(Nmesh(1),Nmesh(2)), Intg(Nmesh(1),Nmesh(2))
         real(kind=dp) :: Int0Ref, IntgRef
         real(kind=dp) :: Int0Norm(Nmesh(1),Nmesh(2)), IntgNorm(Nmesh(1),Nmesh(2))  
         real(kind=dp) :: TotalInt0(Nmesh(1),Nmesh(2)), TotalIntg(Nmesh(1),Nmesh(2))
         real(kind=dp) :: BSEInt(Nmesh(1),Nmesh(2)), BSENorm(Nmesh(1),Nmesh(2))
         real(kind=dp) :: BSEConv(Nmesh(1),Nmesh(2)), BSEConvN(Nmesh(1),Nmesh(2))
         real(kind=dp) :: BSE_BF(Nmesh(1),Nmesh(2)), BSE_DF(Nmesh(1),Nmesh(2))
         real(kind=dp) :: BSE_BFNorm(Nmesh(1),Nmesh(2)), BSE_DFNorm(Nmesh(1),Nmesh(2))
         real(kind=dp) :: refBSE, ECCI_BFbgd, ECCI_DFbgd
         
         real(kind=dp) :: betaD, sigmaBF, sigmaDF, TENDrun, T0run, Tstart, w, TENDT, tilt
         integer       :: IntSteps
         
                  
         !bs parameters
         real(kind=dp) :: Z_Ga, Z_N, DWF_Ga, DWF_N, betaBF, det_sa, Bragg, k, Bohr_r, A_BF, A_DF
         real(kind=dp), parameter :: PI = 3.14159265358979d0
         
         !------------------tilt------------------------
         tilt = -0.87   !-49.6 rad
                
         IWORK(6) = 150000 !maximum number of steps allowed during one call
         RPAR(1:6) = param(1:6) 
         IntSteps = 100 ! intensity integration steps
         dt = TEND*cos(tilt)/IntSteps
         !dt = TEND/0.94_dp/IntSteps !inclined sample TEM
         !dt = TEND/0.34_dp/IntSteps
         
         BSEConv(:,:) = 0.0_dp
         w = RPAR(1)
         DO counterw = 1, 1  !vary w
           RPAR(1) = w - 0.002_dp + (0.004_dp*counterw/2)
           
         TotalInt0(:,:) = 0.0_dp
         TotalIntg(:,:) = 0.0_dp
         BSEInt(:,:) = 0.0_dp

         !-----set mesh---------------------
         Xmesh = -SizeMesh(1)                !start mapping contrast at 
         Ymesh = -SizeMesh(2)                !bottom left corner
             
         dX = 2.0_dp*(-Xmesh)/(Nmesh(1)-1)   !define mesh step
         dY = 2.0_dp*(-Ymesh)/(Nmesh(2)-1)   !
         
         counter = 0
         DO i = Nmesh(1)/2-1, Nmesh(1)/2+1
			counter = counter + 1
			failureGridX(counter) = i
         END DO
         
         counter = 0
         DO  j = Nmesh(2)/2-1, Nmesh(2)/2+1
			counter = counter + 1
			failureGridY(counter) = j
         END DO
         
        !For  the backscattered signal the relative contributions of dark and bright field 
        !are weighted by the cross section for impact ionisation sigma
        !see Rossouw Phil Mag A, 1994
    
        !Dynamical scattering parameters
        Z_Ga = 31.0_dp               !atomic numbers of species in sample
        Z_N = 7.0_dp
    
        DWF_Ga = 0.9969_dp          !exp of unitless DW factors at room temperature
        DWF_N = 0.9973_dp
    
        !betaBF = 2.443_dp           !detector angle=140
        betaBF = 0.5236_dp
        det_sa = 0.087_dp           !detector "solid angle"
        Bragg = 0.05_dp           !betaDF = betaBF+2Bragg from Twigg&Picard 2009
    
        k = 88.73                    !k vector for 30 keV electrons in Angstroms^-1
        Bohr_r = 0.53                !Bohr radius in Angstroms
    
        !number of unit cells N
        !integrate over atoms?
    

       A_BF = 2*PI/(Bohr_r**2 * k**3) * (cos(betaBF)-cos(betaBF+det_sa))&
       &/((1-cos(betaBF))*(1-cos(betaBF+det_sa))) 
       sigmaBF =  A_BF * ((Z_Ga**2 * DWF_Ga) + (Z_N**2 * DWF_N))!!
     
       A_DF = 2*PI/(Bohr_r**2 * k**3) * (cos(betaBF+(2*Bragg))-cos(betaBF+det_sa*(2*Bragg)))&
       &/((1-cos(betaBF+(2*Bragg)))*(1-cos(betaBF+det_sa+2*(Bragg)))) 
       sigmaDF =  A_DF * ((Z_Ga**2 * DWF_Ga) + (Z_N**2 * DWF_N))!!  
                
                
                
         XLoop: DO i = 1, Nmesh(1)
			 print *, Xmesh
             Ymesh = -SizeMesh(2)                     !reset Y
             YLoop: DO j = 1, Nmesh(2)
                 ISTATE = 1                           !start from the begining
                
                 !different surface position
                 TRUN = T0 -Ymesh*tan(tilt)             !70 degrees tilt 
                 
                 phi(1) = DCMPLX(1.0_dp, 0.0_dp)      !reininitialise phi
                 phi(2) = DCMPLX(0.0_dp, 0.0_dp)      !every loop
                 RPAR(7) = Xmesh
                 RPAR(8) = Ymesh
                 RPAR(9) = TEND                       !pass thickness info to ODEs 
                 
                 !avoid convergence failure close to dislocation
                 IF ((ANY(failureGridX.EQ.i)).AND.(ANY(failureGridY.EQ.j))) THEN
                     Int0(i,j) = -1000.0_dp
                     Intg(i,j) = -1000.0_dp
                     print *, i, j
                 ELSE  
					  Tstart = TRUN
        IntDepthIntgr:DO t = 1, IntSteps  
                        TENDrun = Tstart + t*dt
                        
                        CALL ZVODE(FEX,NEQ,phi,TRUN,TENDrun,ITOL,RTOL,ATOL,ITASK,ISTATE,IOPT, &
                       & ZWORK,LZW,RWORK,LRW,IWORK,LIW,JEX,MF,RPAR,IPAR)

                        TotalInt0(i,j) = TotalInt0(i,j) + ABS(phi(1))**2
                        TotalIntg(i,j) = TotalIntg(i,j) + ABS(phi(2))**2 
                   
                           
                      END DO IntDepthIntgr  
                      
                   BSE_BF(i,j) = A_BF*TotalInt0(i,j)
                   BSE_DF(i,j) = A_DF*TotalIntg(i,j)      
                   BSEInt(i,j) = BSE_BF(i,j) + BSE_DF(i,j)
                   BSEConv(i,j) = BSEConv(i,j) + BSEInt(i,j)
                 ENDIF
                   
                 IF (ISTATE .LT. 0) THEN
                     write(*,*)'Grid point where failure',i,j
                     EXIT XLoop
                 ENDIF
             
             
                 YOUT(j) = Ymesh /cos(tilt)
                 Ymesh = Ymesh + dY
             END DO YLoop
             XOUT(i) = Xmesh
             Xmesh = Xmesh + dX
         END DO XLoop
         
             
         !write output file(s)
         !open(101, file = OUTFILE(1)) ! BS BF
         !open(102, file = OUTFILE(2)) ! BS DF
         !open(103, file = OUTFILE(3)) ! Backscattered signal
         open(104, file = OUTFILE(4)) ! convergent beam 

       
       
       print *, "W=", RPAR(1)
     END DO ! Vary w
    
      DO i=1, Nmesh(1)
             DO j=1, Nmesh(2)
                   BSENorm(i,j) = (BSEConv(i,j)-BSEConv(1,1)) *100 / BSEConv(1,1)          
                    write(104,*) XOUT(i), YOUT(j), BSENorm(i,j)   
             END DO
        write(104,*)
          
      END DO
         
     END SUBROUTINE contrast
END MODULE cases


















