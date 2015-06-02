 ! 
      SUBROUTINE ROS2(N,T,Tnext,Hmin,Hmax,Hstart,
     +                   y,Fix,Rconst,AbsTol,RelTol,
     +                   Info)
     
      INCLUDE 'saprcnov_stem_s.h'                                                                                                  
      
C  INPUT ARGUMENTS:
C     y = Vector of (NVAR) concentrations, contains the
C         initial values on input
C     [T, Tnext] = the integration interval
C     Hmin, Hmax = lower and upper bounds for the selected step-size.
C          Note that for Step = Hmin the current computed
C          solution is unconditionally accepted by the error
C          control mechanism.
C     AbsTol, RelTol = (NVAR) dimensional vectors of
C          componentwise absolute and relative tolerances.
C     FUN = name of routine of derivatives. KPP syntax.
C          See the header below.
C     JAC_SP = name of routine that computes the Jacobian, in
C          sparse format. KPP syntax. See the header below.
C     Info(1) = 1  for  autonomous   system
C             = 0  for nonautonomous system
C     Info(2) = 1  for third order embedded formula
C             = 0  for first order embedded formula 
C
C   Note: Stage 3 used to build strongly A-stable order 3 formula for error control
C            Embed3 = (Info(2).EQ.1)
C         if Embed3 = .true. then the third order embedded formula is used
C                     .false. then a first order embedded  formula is used
C
C
C  OUTPUT ARGUMENTS:
C     y = the values of concentrations at Tend.
C     T = equals Tend on output.
C     Info(2) = # of FUN calls.
C     Info(3) = # of JAC_SP calls.
C     Info(4) = # of accepted steps.
C     Info(5) = # of rejected steps.
 
      REAL*8     K1(NVAR), K2(NVAR), K3(NVAR)
      REAL*8     F1(NVAR), JAC(LU_NONZERO_V)
      REAL*8     DFDT(NVAR)
      REAL*8     Fix(NFIX), Rconst(NREACT), Rad(1)
      REAL*8     Hmin,Hmax,Hstart,ghinv,uround      
      REAL*8     y(NVAR), ynew(NVAR)
      REAL*8     AbsTol(NVAR), RelTol(NVAR)
      REAL*8     T, Tnext, H, Hold, Tplus
      REAL*8     ERR, factor, facmax
      INTEGER    n,nfcn,njac,Naccept,Nreject,i,j
      INTEGER    Info(5)
      LOGICAL    IsReject, Autonom, Embed3
      EXTERNAL    FUN, JAC_SP                                                                                                

      REAL*8 gamma, m1, m2, alpha, beta1, beta2, delta, w, e
 
c     Initialization of counters, etc.
      Autonom = Info(1) .EQ. 1
      Embed3  = Info(2) .EQ. 1
      uround  = 1.d-15
      dround  = dsqrt(uround)
      H = DMAX1(Hstart,DMAX1(1.d-8, Hmin))
      Tplus = T
      IsReject = .false.
      Naccept  = 0
      Nreject  = 0
      Nfcn     = 0
      Njac     = 0
      
C     Method Parameters      
      gamma = 1.d0 + 1.d0/sqrt(2.d0)
      a21   = - 1.d0/gamma
      m1 = -3.d0/(2.d0*gamma)
      m2 = -1.d0/(2.d0*gamma)
      c31 = -1.0D0/gamma**2*(1.0D0-7.0D0*gamma+9.0D0*gamma**2)
     &         /(-1.0D0+2.0D0*gamma)
      c32 = -1.0D0/gamma**2*(1.0D0-6.0D0*gamma+6.0D0*gamma**2)
     &          /(-1.0D0+2.0D0*gamma)/2
      gamma3 = 0.5D0 - 2*gamma
      d1 = ((-9.0D0*gamma+8.0D0*gamma**2+2.0D0)/gamma**2/
     &          (-1.0D0+2*gamma))/6.0D0
      d2 = ((-1.0D0+3.0D0*gam)/gamma**2/
     &          (-1.0D0+2.0D0*gamma))/6.0D0
      d3 = -1.0D0/(3.0D0*gamma)
 
C === Starting the time loop ===
 10    CONTINUE
       Tplus = T + H
       if ( Tplus .gt. Tnext ) then
          H = Tnext - T
          Tplus = Tnext
       end if
 
       call JacVar_SP( Y, Rad, Fix, Rconst, JAC )
 
       Njac = Njac+1
       ghinv = -1.0d0/(gamma*H)
       DO 20 j=1,NVAR
         JAC(LU_DIAG_V(j)) = JAC(LU_DIAG_V(j)) + ghinv
 20    CONTINUE
       call KppDecomp (NVAR, JAC, ier)
 
       if (ier.ne.0) then
         if ( H.gt.Hmin) then
            H = 5.0d-1*H
            go to 10
         else
            print *,'IER <> 0, H=',H
            stop
         end if
       end if
 
       call FunVar( Y, Rad, Fix, Rconst, F1 )
 
C ====== NONAUTONOMOUS CASE ===============
       IF (.not. Autonom) THEN
         tau = dsign(dround*dmax1( 1.0d-6, dabs(T) ), T)
         call FunVar( Y, Rad, Fix, Rconst, K2 )
         nfcn=nfcn+1
         DO 30 j = 1,NVAR
           DFDT(j) = ( K2(j)-F1(j) )/tau
 30      CONTINUE
       END IF ! .NOT.Autonom
 
C ----- STAGE 1 -----
       delta = gamma*H
       DO 40 j = 1,NVAR
          K1(j) =  F1(j) 
 40    CONTINUE
       IF (.NOT.Autonom) THEN
         DO 45 j = 1,NVAR
           K1(j) = K1(j) + delta*DFDT(j)
 45      CONTINUE       
       END IF ! .NOT.Autonom
       call KppSolve (JAC, K1)
 
C ----- STAGE 2  -----
       DO 50 j = 1,NVAR
         ynew(j) = y(j) + a21*K1(j)
 50    CONTINUE
       call FunVar( Ynew, Rad, Fix, Rconst, F1 )
       nfcn=nfcn+1
       beta = 2.d0/(gamma*H)
       delta = -gamma*H
       DO 55 j = 1,NVAR
         K2(j) = F1(j) + beta*K1(j) 
 55    CONTINUE
       IF (.NOT.Autonom) THEN
         DO 56 j = 1,NVAR
           K2(j) = K2(j) + delta*DFDT(j)
 56      CONTINUE       
       END IF ! .NOT.Autonom
       call KppSolve (JAC, K2)
 
C ----- STAGE 3  -----
       IF (Embed3) THEN
       beta1 = -c31/H
       beta2 = -c32/H
       delta = gamma3*H
       DO 57 j = 1,NVAR
         K3(j) = F1(j) + beta1*K1(j) + beta2*K2(j) 
 57    CONTINUE
       IF (.NOT.Autonom) THEN
         DO 58 j = 1,NVAR
           K3(j) = K3(j) + delta*DFDT(j)
 58      CONTINUE       
       END IF ! .NOT.Autonom
       CALL KppSolve (JAC, K3)
       END IF ! Embed3
 


C ---- The Solution ---
       DO 120 j = 1,NVAR
         ynew(j) = y(j) + m1*K1(j) + m2*K2(j)
 120   CONTINUE
 
 
C ====== Error estimation ========
 
        ERR=0.d0
        DO 130 i=1,NVAR
           w = AbsTol(i) + RelTol(i)*DMAX1(DABS(y(i)),DABS(ynew(i)))
	   IF ( Embed3 ) THEN
	     e = d1*K1(i) + d2*K2(i) + d3*K3(i)
	   ELSE
             e = 1.d0/(2.d0*gamma)*(K1(i)+K2(i))
	   END IF  ! Embed3 
           ERR = ERR + ( e/w )**2
 130    CONTINUE
        ERR = DMAX1( uround, DSQRT( ERR/NVAR ) )
	
C ======= Choose the stepsize ===============================
 
        IF ( Embed3 ) THEN
            elo = 3.0D0 ! estimator local order
	ELSE
	    elo = 2.0D0
	END IF    
        factor = DMAX1(2.0D-1,DMIN1(6.0D0,ERR**(1.0D0/elo)/.9D0))
        Hnew   = DMIN1(Hmax,DMAX1(Hmin, H/factor))
        Hold = H
         
C ======= Rejected/Accepted Step ============================
 
        IF ( (ERR.gt.1).and.(H.gt.Hmin) ) THEN
          IsReject = .true.
	  H = DMIN1(H/10,Hnew)
          Nreject  = Nreject+1
        ELSE
          DO 140 i=1,NVAR
             y(i)  = ynew(i)
 140      CONTINUE
          T = Tplus
	  IF (.NOT.IsReject) THEN
	      H = Hnew   ! Do not increase stepsize if previous step was rejected
	  END IF    
          IsReject = .false.
          Naccept = Naccept+1
!	  CALL TRAJISTORE(y,hold)    
        end if        
C ======= End of the time loop ===============================
      if ( T .lt. Tnext ) go to 10 
 
C ======= Output Information =================================
      Info(2) = Nfcn
      Info(3) = Njac
      Info(4) = Naccept
      Info(5) = Nreject      
     
      RETURN
      END
 

! 
      SUBROUTINE ROS2_cadj(N,T,Tnext,Hmin,Hmax,Hstart,
     +                   y,Lambda,Fix,Rconst,AbsTol,RelTol,
     +                   Info)
      INCLUDE 'saprcnov_stem_s.h'                                                                                                  
      
C  INPUT ARGUMENTS:
C     y = Vector of (NVAR) concentrations, contains the
C         initial values on input
C     [T, Tnext] = the integration interval
C     Hmin, Hmax = lower and upper bounds for the selected step-size.
C          Note that for Step = Hmin the current computed
C          solution is unconditionally accepted by the error
C          control mechanism.
C     AbsTol, RelTol = (NVAR) dimensional vectors of
C          componentwise absolute and relative tolerances.
C     FUN = name of routine of derivatives. KPP syntax.
C          See the header below.
C     JAC_SP = name of routine that computes the Jacobian, in
C          sparse format. KPP syntax. See the header below.
C     Info(1) = 1  for  autonomous   system
C             = 0  for nonautonomous system
C     Info(2) = 1  for third order embedded formula
C             = 0  for first order embedded formula 
C
C   Note: Stage 3 used to build strongly A-stable order 3 formula for error control
C            Embed3 = (Info(2).EQ.1)
C         if Embed3 = .true. then the third order embedded formula is used
C                     .false. then a first order embedded  formula is used
C
C
C  OUTPUT ARGUMENTS:
C     y = the values of concentrations at Tend.
C     T = equals Tend on output.
C     Info(2) = # of FUN calls.
C     Info(3) = # of JAC_SP calls.
C     Info(4) = # of accepted steps.
C     Info(5) = # of rejected steps.

      integer, parameter :: max_no_steps = 200      
      real*8 :: Trajectory(NVAR,max_no_steps)
      real*8 :: StepSize(max_no_steps)
 
      REAL*8     K1(NVAR), K2(NVAR), K3(NVAR)
      REAL*8     F1(NVAR), JAC(LU_NONZERO_V)
      REAL*8     DFDT(NVAR), Rad(NRAD)
      REAL*8     Fix(NFIX), Rconst(NREACT)
      REAL*8     Hmin,Hmax,Hstart,ghinv,uround      
      REAL*8     y(NVAR), ynew(NVAR)
      REAL*8     AbsTol(NVAR), RelTol(NVAR)
      REAL*8     T, Tnext, H, Hold, Tplus
      REAL*8     ERR, factor, facmax
      REAL*8     Lambda(NVAR), K11(NVAR), JAC1(LU_NONZERO_V)
      INTEGER    n,nfcn,njac,Naccept,Nreject,i,j
      INTEGER    Info(5)
      LOGICAL    IsReject, Autonom, Embed3
      EXTERNAL    FUN, JAC_SP                                                                                                

      REAL*8 gamma, m1, m2, alpha, beta1, beta2, delta, w, e
      REAL*8 ginv
c     Initialization of counters, etc.
      Autonom = Info(1) .EQ. 1
      Embed3  = Info(2) .EQ. 1
      uround  = 1.d-15
      dround  = dsqrt(uround)
      H = DMAX1(Hstart,DMAX1(1.d-8, Hmin))
      Tplus = T
      IsReject = .false.
      Naccept  = 0
      Nreject  = 0
      Nfcn     = 0
      Njac     = 0
      
C     Method Parameters      
      gamma = 1.d0 + 1.d0/sqrt(2.d0)
      a21   = - 1.d0/gamma
      m1 = -3.d0/(2.d0*gamma)
      m2 = -1.d0/(2.d0*gamma)
      c31 = -1.0D0/gamma**2*(1.0D0-7.0D0*gamma+9.0D0*gamma**2)
     &         /(-1.0D0+2.0D0*gamma)
      c32 = -1.0D0/gamma**2*(1.0D0-6.0D0*gamma+6.0D0*gamma**2)
     &          /(-1.0D0+2.0D0*gamma)/2
      gamma3 = 0.5D0 - 2*gamma
      d1 = ((-9.0D0*gamma+8.0D0*gamma**2+2.0D0)/gamma**2/
     &          (-1.0D0+2*gamma))/6.0D0
      d2 = ((-1.0D0+3.0D0*gam)/gamma**2/
     &          (-1.0D0+2.0D0*gamma))/6.0D0
      d3 = -1.0D0/(3.0D0*gamma)

      Trajectory(1:NVAR,1) = ynew(1)
 
C === Starting the time loop ===
 10    CONTINUE
       Tplus = T + H
       if ( Tplus .gt. Tnext ) then
          H = Tnext - T
          Tplus = Tnext
       end if
 
       call JacVar_SP( Y, Rad, Fix, Rconst, JAC )
 
       Njac = Njac+1
       ghinv = -1.0d0/(gamma*H)
       DO 20 j=1,NVAR
         JAC(LU_DIAG_V(j)) = JAC(LU_DIAG_V(j)) + ghinv
 20    CONTINUE
       call KppDecomp (NVAR, JAC, ier)
 
       if (ier.ne.0) then
         if ( H.gt.Hmin) then
            H = 5.0d-1*H
            go to 10
         else
            print *,'IER <> 0, H=',H
            stop
         end if
       end if
 
       call FunVar( Y, Rad, Fix, Rconst, F1 )
 
C ====== NONAUTONOMOUS CASE ===============
       IF (.not. Autonom) THEN
         tau = dsign(dround*dmax1( 1.0d-6, dabs(T) ), T)
         call FunVar( Y, Rad, Fix, Rconst, K2 )
         nfcn=nfcn+1
         DO 30 j = 1,NVAR
           DFDT(j) = ( K2(j)-F1(j) )/tau
 30      CONTINUE
       END IF ! .NOT.Autonom
 
C ----- STAGE 1 -----
       delta = gamma*H
       DO 40 j = 1,NVAR
          K1(j) =  F1(j) 
 40    CONTINUE
       IF (.NOT.Autonom) THEN
         DO 45 j = 1,NVAR
           K1(j) = K1(j) + delta*DFDT(j)
 45      CONTINUE       
       END IF ! .NOT.Autonom
       call KppSolve (JAC, K1)
 
C ----- STAGE 2  -----
       DO 50 j = 1,NVAR
         ynew(j) = y(j) + a21*K1(j)
 50    CONTINUE
       call FunVar( Ynew, Rad, Fix, Rconst, F1 )
       nfcn=nfcn+1
       beta = 2.d0/(gamma*H)
       delta = -gamma*H
       DO 55 j = 1,NVAR
         K2(j) = F1(j) + beta*K1(j) 
 55    CONTINUE
       IF (.NOT.Autonom) THEN
         DO 56 j = 1,NVAR
           K2(j) = K2(j) + delta*DFDT(j)
 56      CONTINUE       
       END IF ! .NOT.Autonom
       call KppSolve (JAC, K2)
 
C ----- STAGE 3  -----
       IF (Embed3) THEN
       beta1 = -c31/H
       beta2 = -c32/H
       delta = gamma3*H
       DO 57 j = 1,NVAR
         K3(j) = F1(j) + beta1*K1(j) + beta2*K2(j) 
 57    CONTINUE
       IF (.NOT.Autonom) THEN
         DO 58 j = 1,NVAR
           K3(j) = K3(j) + delta*DFDT(j)
 58      CONTINUE       
       END IF ! .NOT.Autonom
       CALL KppSolve (JAC, K3)
       END IF ! Embed3
 


C ---- The Solution ---
       DO 120 j = 1,NVAR
         ynew(j) = y(j) + m1*K1(j) + m2*K2(j)
 120   CONTINUE
 
 
C ====== Error estimation ========
 
        ERR=0.d0
        DO 130 i=1,NVAR
           w = AbsTol(i) + RelTol(i)*DMAX1(DABS(y(i)),DABS(ynew(i)))
	   IF ( Embed3 ) THEN
	     e = d1*K1(i) + d2*K2(i) + d3*K3(i)
	   ELSE
             e = 1.d0/(2.d0*gamma)*(K1(i)+K2(i))
	   END IF  ! Embed3 
           ERR = ERR + ( e/w )**2
 130    CONTINUE
        ERR = DMAX1( uround, DSQRT( ERR/NVAR ) )
	
C ======= Choose the stepsize ===============================
 
        IF ( Embed3 ) THEN
            elo = 3.0D0 ! estimator local order
	ELSE
	    elo = 2.0D0
	END IF    
        factor = DMAX1(2.0D-1,DMIN1(6.0D0,ERR**(1.0D0/elo)/.9D0))
        Hnew   = DMIN1(Hmax,DMAX1(Hmin, H/factor))
        Hold = H
         
C ======= Rejected/Accepted Step ============================
 
        IF ( (ERR.gt.1).and.(H.gt.Hmin) ) THEN
          IsReject = .true.
	  H = DMIN1(H/10,Hnew)
          Nreject  = Nreject+1
        ELSE
          DO 140 i=1,NVAR
             y(i)  = ynew(i)
 140      CONTINUE
          T = Tplus
	  IF (.NOT.IsReject) THEN
	      H = Hnew   ! Do not increase stepsize if previous step was rejected
	  END IF    
          IsReject = .false.
          Naccept = Naccept+1
	  if (Naccept+1>max_no_steps) then
	    print*,'Error in Ros2: more steps than allowed'
	    stop
	  end if  
	  Trajectory(1:NVAR,Naccept+1) = ynew(1:NVAR)
	  StepSize(Naccept) = Hold
!	  CALL TRAJISTORE(y,hold)    
        end if        
C ======= End of the time loop ===============================
      if ( T .lt. Tnext ) go to 10 
 
C ======= Output Information =================================
      Info(2) = Nfcn
      Info(3) = Njac
      Info(4) = Naccept
      Info(5) = Nreject      

      ginv = 1.d0/gamma
C -- The backwards loop for the CONTINUOUS ADJOINT    
      do istep = Naccept,1,-1
        
	h = StepSize(istep)
	y(1:NVAR) = Trajectory(1:NVAR,istep+1)
	gHinv = -ginv/H  

	call JacVar_SP(Y, Rad, Fix, Rconst, JAC)
        JAC1(1:LU_NONZERO_V)=JAC(1:LU_NONZERO_V)
        do j=1,NVAR
          JAC(lu_diag_v(j)) = JAC(lu_diag_v(j)) + gHinv
        end do
        call KppDecomp (NVAR,JAC,ier) 
ccc equivalent to function evaluation in forward integration  
ccc is J^T*Lambda in backward integration    
        call JacVarTR_SP_Vec ( JAC1, Lambda, F1)      
      
C ----- STAGE 1 (AUTONOMOUS) -----
        K11(1:NVAR) =  F1(1:NVAR) 
        call KppSolveTR (JAC,K11,K1) 
C ----- STAGE 2 (AUTONOMOUS) -----       
	y(1:NVAR) = Trajectory(1:NVAR,istep)
        call JacVar_SP(Y, Rad, Fix, Rconst, JAC1)
        Ynew(1:NVAR) = Lambda(1:NVAR) - ginv*K1(1:NVAR)
        call JacVarTR_SP_Vec ( JAC1, Ynew, F1) 
        beta = -2.d0*ghinv
        K11(1:NVAR) = F1(1:NVAR) + beta*K1(1:NVAR)
        call KppSolveTR (JAC,K11,K2)
c ---- The solution
        Lambda(1:NVAR) = Lambda(1:NVAR)+m1*K1(1:NVAR)+m2*K2(1:NVAR)                

      end do ! istep
     
     
      RETURN
      END
 
