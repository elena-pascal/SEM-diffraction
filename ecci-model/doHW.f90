PROGRAM doHW
     USE cases
     IMPLICIT NONE
 
     integer                          :: tORw, counter
     integer                          :: NOUT, Nmesh(2)
     integer                          :: ISTATE
     real(kind=dp)                    :: T, TEND, W, WEND
     real(kind=dp)                    :: strParam(6)
     real(kind=dp)                    :: SizeMesh(2)
     real(kind=dp)                    :: gb
     complex(kind=dp), dimension(NEQ) :: Y
     character(30)                    :: OUTFILE(4), string
	 character(*), parameter :: fileplace = "Edge/"
   
     ! initialise everything   
     Y = (/(1.0_dp, 0.0_dp), (0.0_dp, 0.0_dp)/)    !Y initial conditions
     T = 0.0_dp                                    !T start
     ISTATE = 1                                    !start from begining
     NOUT = 200                                    !how many points of integration
     strParam(2) = 0.165_dp                        !Xsig/Xsig' ratio 
     strParam(3) = 0.163_dp                        !Xsi0'/Xsig' ratio
     strParam(5) = 0.26_dp                         !nu for wurtzite GaN
     !
     tORw = 3                                      !choose case

     SELECT CASE(tORw)
         CASE (1)                                  ! VARY T 
             strParam(1) = 0.0_dp                  !define w
             TEND = 20.0_dp                         !max thickness
             OUTFILE(1) = 'intensitiesT.dat'       !name outfiles
             OUTFILE(2) = 'phase.dat'              !
          
             CALL fringes(Y, T, TEND, strParam, NOUT, ISTATE, OUTFILE)
         
             IF(ISTATE.LT.0) THEN
                 print *, '***** Error halt.  ISTATE =',ISTATE 
             ELSE
                 print *, OUTFILE(1), ' written successfully'
                 print *, OUTFILE(2), ' written successfully'
             ENDIF
             
         CASE(2)                                   ! VARY W
             strParam(1) = 0.0_dp                 !W start
             WEND = 3.0_dp                         !W end
             TEND = 4.0_dp                         !max thickness
             OUTFILE(1) = 'intensitiesW.dat'       !name outfiles
          
             CALL bands(Y, T, TEND, strParam, WEND, NOUT, ISTATE, OUTFILE)
             
             IF(ISTATE.LT.0) THEN
                 print *, '***** Error halt.  ISTATE =',ISTATE 
             ELSE
                 print *, OUTFILE(1), ' written succesfully'    
             END IF
     
!========================= HW on a mesh =====================================
         CASE(3)                                   ! VARY X,Y
         
 ! GaN parameters
             ! g//b//x
              strParam(1) = 0_dp                !define w
              !strParam(4) = 0.39_dp             
              strParam(4) = 1_dp                  !define g.b -edge
              !strParam(4) = 3.25_dp							!screw
              !TEND = 34_dp                        !define thickness t
              TEND = 4_dp                         !in the direction of rinc in units of extintion distance
              strParam(6) = 0.05_dp              !thetaB
!-----mesh set-up------------------------------------------

 !           SizeMesh = (/0.02_dp,0.02_dp/)    !define HW sampling mesh

             SizeMesh = (/300_dp,300_dp/)    ! xlonger
             Nmesh = (/250, 250/)                  !number of sampling points
                                                   !even number please
				write(string, '(F4.1)')  
				
				OUTFILE(4) = fileplace//'edge_g1_b0'
				
				!contrast(phi, T0, TEND, param, Nmesh, SizeMesh, ISTATE, OUTFILE)
				CALL contrast(Y, T, TEND, strParam, Nmesh, SizeMesh, ISTATE, OUTFILE)
             
				IF(ISTATE.LT.0) THEN
					print *, '***** Error halt.  ISTATE =',ISTATE 
				ELSE
					print *, OUTFILE(:), ' written succesfully'   
				ENDIF
     END SELECT
         
 
END PROGRAM doHW

