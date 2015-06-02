ccc discrete adjoint modeling with ros2
ccc implements the algorithm presented by Daescu et al. (JCP 165 (2), 2000)
  
      subroutine dadj_ros2(y,u,h,t)
      
      include 'saprcnov_stem_s.h' 
          
      real*8 y(nvar),u(nvar),v(nvar)
      real*8 omega(nvar),omega1(nvar),theta(nvar)
             
      REAL*8    K1(NVAR),K2(NVAR)
      REAL*8    gamma,F1(NVAR),thetaj(nvar)
      REAL*8    JAC(LU_NONZERO_V),JAC1(LU_NONZERO_V)
      REAL*8    JAC2(LU_NONZERO_V)
      REAL*8    ghinv,m1,m2,beta,alpha
      REAL*8    ynew(NVAR)
      REAL*8    t,H
      EXTERNAL   FUN, JAC_SP  
      
      integer   i,j,ier     
C HESS - Hessian (derivative of Jacvar)                            
      REAL*8 HESS(NHESS+1)
      REAL*8 HTU1(NVAR),HTU2(NVAR)
      
            
c     Initialization of counters, etc.
      gamma = (2.d0+DSQRT(2.d0))/2.d0        
      m1 = -3.d0/(2.d0*gamma)
      m2 = -1.d0/(2.d0*gamma)
      alpha = - 1.d0/gamma
      beta = 2.d0/(gamma*H)
      gHinv = -1.0d0/(gamma*H)
      
C ===== Start the time loop ======
           
       call JAC_SP(NVAR, T, y, JAC)
c save the jacobian value in jac1      
       do i=1,LU_NONZERO_V
        JAC1(i)=JAC(i)
       enddo 
       
       do 20 j=1,NVAR
         JAC(lu_diag_v(j)) = JAC(lu_diag_v(j)) + gHinv
 20    continue

        
       call KppDecomp (NVAR,JAC,ier)

       call FUN(NVAR, T, y, F1)

C ----- STAGE 1 (AUTONOMOUS) -----
         do 70 j = 1,NVAR
           K1(j) =  F1(j)
 70      continue

       call KppSolve (JAC,K1)
       
C ----- STAGE 2 (AUTONOMOUS) -----
         do 80 j = 1,NVAR
           ynew(j) = y(j) + alpha*K1(j) 
 80      continue
         call FUN(NVAR, T+H, ynew, F1)
         do 90 j = 1,NVAR
           K2(j) = F1(j) + beta*K1(j)
 90      continue

         call KppSolve (JAC,K2)
c set adjoint var zero if the fwd values were set zero	
c uncomment the lines below if fwd integration used clipping 
c        do i=1,nvar
c         if(y(i)+ m1*K1(j) + m2*K2(j).le.0.0d0)then
c           u(i)=0.0d0
c         endif
c        enddo  


       call KppSolveTR(jac,u,v)
       call JAC_SP(NVAR, T+h, ynew, JAC2)
       call JacVarTR_SP_Vec ( JAC2, v, omega )     
       
       do i=1,nvar
         omega1(i)=m2*(alpha*omega(i)+beta*v(i))
       enddo
       
       call KppSolveTR (jac,omega1,theta)

       do i=1,nvar
        theta(i) =theta(i) + m1*v(i) 
       enddo 
       
cccc       call JACVARTR_SP_VEC ( JAC1, theta, thetaj ) 
cccc   replaced using properties of theta by 
       do i=1,nvar 
         thetaj(i) = omega1(i)+ m1*u(i) - ghinv*theta(i)
       enddo	       
cccc       
       call HessVar ( y, RAD, FIX, RCONST, HESS )
       call HessVarTR_Vec ( HESS, theta, k1, HTU1 )
       call HessVarTR_Vec ( HESS, v, k2, HTU2 )       
       do i=1,nvar
         u(i)=u(i)+thetaj(i)+m2*omega(i)-htu1(i)-m2*htu2(i)        
       enddo
	        
      return
      end


ccc continuous adjoint modeling with ros2

      subroutine cadj_ros2(y,advar,h,t)      
      INCLUDE 'saprcnov_stem_s.h'     
      real*8  y(nvar),advar(nvar)
      REAL*8  K1(NVAR)
      REAL*8  K2(NVAR),K11(nvar)
      REAL*8  F1(NVAR)
      REAL*8  JAC(LU_NONZERO_V),JAC1(LU_NONZERO_V)
      REAL*8  gamma,ginv,ghinv,m1,m2,beta
      REAL*8  ynew(NVAR)
      REAL*8  t,H
      
      integer   i,j,ier     
      
c     Initialization of counters, etc.
       
      gamma = 0.5d0*(2.0d0 + dsqrt(2.d0)) 
      ginv  = 2.0d0 - dsqrt(2.d0) 
      gHinv = -ginv/H  
      m1 = -1.5d0*ginv
      m2 = -0.5d0*ginv     
            
      call JAC_SP(NVAR, T, y, JAC)
      do i=1,LU_NONZERO_V
        JAC1(i)=JAC(i)
      enddo 
            
      do 20 j=1,NVAR
         JAC(lu_diag_v(j)) = JAC(lu_diag_v(j)) + gHinv
 20   continue 
      call KppDecomp (NVAR,JAC,ier) 
ccc equivalent to function evaluation in forward integration  
ccc is J^T*advar in backward integration    
      call JacVarTR_SP_Vec ( JAC1, advar, F1)      
      
C ----- STAGE 1 (AUTONOMOUS) -----

      do 70 j = 1,NVAR
           K11(j) =  F1(j)
 70   continue
 
      call KppSolveTR (jac,K11,K1) 
C ----- STAGE 2 (AUTONOMOUS) -----       
      do 80 j = 1,NVAR
           ynew(j) = advar(j) - ginv*K1(j)
 80   continue 
      call JacVarTR_SP_Vec ( JAC1, ynew, F1) 
      beta = -2.d0*ghinv
      do 90 j = 1,NVAR
           K11(j) = F1(j) + beta*K1(j)
 90   continue 
      call KppSolveTR (jac,K11,K2)
      do 110 j = 1,NVAR
         advar(j) =   advar(j) + m1*K1(j) + m2*K2(j)                
 110  continue 
 
       return
       end
                    
                      
