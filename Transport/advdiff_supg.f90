! Transport Routines
! Note: The concentration is expected to be molec/cm^3 or similar
!       Do not use with relative fractions (ppb etc).


!*************************************************************************
      subroutine tranx_supg(ix,iy,iz,nstart,nend,s1,u,kh,sx,dt,dx)
      use ParallelDataMap
!**************************************************************************
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7--
      include 'aqmax.param'
      !include 'aqindx.cmm'
      integer :: ix, iy, iz, nstart, nend
      real :: s1(ix,iy,iz,*), u(ix,iy,*), kh(ix,iy,*)
      real :: sx(iy,iz,2,*), dx(*), dt
      !real :: bndf(2,mxspg),bndv(2,mxspg)
      !real :: gen(mxgr*mxspg)
      real*8 :: wind(1:ix), dif(1:ix), air(1:ix), conc(1:ix,mxspg)
      real*8 :: sb(1:2,mxspg), Deltax, Deltat
      real*8 :: SurfaceEm(mxspg), VolumeEm(ix,mxspg)
      real*8 :: Depvel(mxspg),X(ix)
      integer :: num, i, j, k, Nsteps
!
!                           transport in x-direction

      num = nend-nstart+1         ! number of transported species
      
      Deltax = dx(1)       ! grid size: this works for uniform grids only
      ! max step = 300.0
      Nsteps = int(dt/300.0)+1
      Deltat = dt/dble(Nsteps)          ! time interval
      
      do i=1,ix
        X(i) = (i-1)*DeltaX
      end do	
      SurfaceEm(1:num) = 0.d0
      VolumeEm(1:ix,1:num) = 0.d0
      DepVel(1:num)        = 0.d0      
      
      do j=1,no_of_xslices(MyId)  ! 1,iy
      do k=1,iz
      
      
      Conc(1:ix,1:num) = s1(1:ix,j,k,nstart:nend)
      Wind(1:ix) = u(1:ix,j,k)
      Dif(1:ix) = kh(1:ix,j,k)
      ! Air(1:ix) = s1(1:ix,j,k,iair)
      Sb(1:2,1:num) = sx(j,k,1:2,1:num)
      ! CALL ADVDIFF_U(Deltat,ix,num,Deltax,Wind,Dif,Conc,Sb)
      CALL ADVDIFF_SUPG(Deltat, Nsteps, ix, num, X, Wind, Dif, Conc, &
                     Sb, SurfaceEm, VolumeEm, DepVel)
      s1(1:ix,j,k,nstart:nend) = max(Conc(1:ix,1:num),0.d0)
      end do ! iz
      end do ! iy
      return
      end

!*********************************************************************
      subroutine trany_supg(ix,iy,iz,nstart,nend,s1,v,kh,sy,dt,dy)
      use ParallelDataMap
!*********************************************************************
      include 'aqmax.param'
      !include 'aqindx.cmm'
      integer :: ix, iy, iz, nstart, nend
      real ::  s1(ix,iy,iz,*), v(ix,iy,*), kh(ix,iy,*)
      real ::  sy(ix,iz,2,*), dy(*), dt
      !real ::  bndf(2,mxspg),bndv(2,mxspg)
      !real ::  gen(mxgr*mxspg)
      real*8 :: Wind(1:iy), Dif(1:iy), Air(1:iy), Conc(1:iy,mxspg)
      real*8 :: Sb(1:2,mxspg), Deltay, Deltat
      integer :: num, i, j, k, Nsteps
      real*8 :: SurfaceEm(mxspg),VolumeEm(iy,mxspg)
      real*8 :: DepVel(mxspg),Y(iy)
!      
      num=nend-nstart+1    ! number of transported species
      Deltay = dy(1)       ! grid size: this works for uniform grids only
      ! max step = 300.0
      Nsteps = int(dt/300.0)+1
      Deltat = dt/dble(Nsteps)          ! time interval

      do i=1,iy
        Y(i) = (i-1)*DeltaY
      end do	
      SurfaceEm(1:num) = 0.d0
      VolumeEm(1:iy,1:num) = 0.d0
      DepVel(1:num)        = 0.d0      
      
      do i=1,no_of_yslices(MyId)   ! 1,ix
      do k=1,iz
      
      Conc(1:iy,1:num) = s1(i,1:iy,k,nstart:nend)
      Wind(1:iy) = v(i,1:iy,k)
      Dif(1:iy) = kh(i,1:iy,k)
      ! Air(1:iy) = s1(i,1:iy,k,iair)
      Sb(1:2,1:num) = sy(i,k,1:2,1:num)
      ! CALL ADVDIFF_U(Deltat,iy,num,Deltay,Wind,Dif,Conc,Sb)
      CALL ADVDIFF_SUPG(Deltat, Nsteps, iy, num, Y, Wind, Dif, &
                   Conc, Sb, SurfaceEm, VolumeEm, DepVel)
      s1(i,1:iy,k,nstart:nend) = max(Conc(1:iy,1:num),0.d0)
      
      end do ! iz
      end do ! ix
      return
      end
      
      
!*********************************************************************
      subroutine tranz_supg(ix,iy,iz,nstart,nend,s1,w,kv,q,em,vg,sz,dt,dz)
      use ParallelDataMap
!*********************************************************************
      include 'aqmax.param'
      !include 'aqindx.cmm'
!      Arguments:      
      integer :: ix, iy, iz, nstart, nend
      real :: s1(ix,iy,iz,*), w(ix,iy,iz), kv(ix,iy,iz)
      real :: q(ix,iy,*), em(ix,iy,iz,*), vg(ix,iy,*), sz(ix,iy,*)
      real :: dt, dz(ix,iy,*)
! Local variables:      
      real*8 :: Deltat, Z(iz), Wind(iz), Dif(iz)
      real*8 :: C(iz,mxspg), Bdry(2,mxspg) 
      real*8 :: SurfaceEm(mxspg), VolumeEm(iz,mxspg), DepVel(mxspg) 
      integer :: i, j, k, Nsteps    

      num = nend-nstart+1 ! number of species
      ! max step = 300.0
      Nsteps = int(dt/300.0)+1
      Deltat = dt/dble(Nsteps)          ! time interval
      
      do i=1,no_of_yslices(MyId)  ! i=1,ix
      do j=1,iy
      
      ! The vertical grid
      Z(1) = 0.d0
      do k = 1, iz-1
        Z(k+1) = Z(k) + dz(i,j,k)
      end do	
      
      Wind(1:iz)  = w(i,j,1:iz)
      Dif(1:iz)   = kv(i,j,1:iz)
      C(1:iz,1:num) = s1(i,j,1:iz,nstart:nend)
      Bdry(1,1:num) = s1(i,j,1,nstart:nend)
      Bdry(2,1:num) = sz(i,j,nstart:nend)
      SurfaceEm(1:num) = q(i,j,nstart:nend)
      VolumeEm(1:iz,1:num) = em(i,j,1:iz,nstart:nend)
      DepVel(1:num)        = vg(i,j,nstart:nend)
        
      !call ADVDIFF_FVZ(Deltat, Nsteps, iz, num, Z, Wind, Dif, &
      !              C, BDRY, SurfaceEm, VolumeEm, DepVel)   
      call ADVDIFF_SUPG(Deltat, Nsteps, iz, num, Z, Wind, Dif, &
                    Conc, Bdry, SurfaceEm, VolumeEm, Depvel)

      s1(i,j,1:iz,nstart:nend) = max(C(1:iz,1:num),0.d0)
		       
      end do ! j
      end do ! i
!
      return
      end
      

!*************************************************************************
      subroutine tranx_mf(ix,iy,iz,nstart,nend,s1,u,kh,sx,dt,dx)
      use ParallelDataMap
!**************************************************************************
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7--
      include 'aqmax.param'
      !include 'aqindx.cmm'
      integer :: ix, iy, iz, nstart, nend
      real :: s1(ix,iy,iz,*), u(ix,iy,*), kh(ix,iy,*)
      real :: sx(iy,iz,2,*), dx(*), dt
      !real :: bndf(2,mxspg),bndv(2,mxspg)
      !real :: gen(mxgr*mxspg)
      real*8 :: wind(1:ix), dif(1:ix), air(1:ix), conc(1:ix,mxspg)
      real*8 :: sb(1:2,mxspg), Deltax, Deltat
      real*8 :: SurfaceEm(mxspg), VolumeEm(ix,mxspg),Vd(mxspg),X(ix)
      integer :: num, i, j, k
!
!                           transport in x-direction

      num = nend-nstart+1         ! number of transported species
      
      Deltax = dx(1)       ! grid size: this works for uniform grids only
      Deltat = dt          ! time interval
      
      do i=1,ix
        X(i) = (i-1)*DeltaX
      end do	
      SurfaceEm(1:num) = 0.d0
      VolumeEm(1:ix,1:num) = 0.d0
      Vd(1:num)        = 0.d0      
      
      do j=1,no_of_xslices(MyId)  ! 1,iy
      do k=1,iz
      
      
      Conc(1:ix,1:num) = s1(1:ix,j,k,nstart:nend)
      Wind(1:ix) = u(1:ix,j,k)
      Dif(1:ix) = kh(1:ix,j,k)
      Air(1:ix) = s1(1:ix,j,k,iair)
      Sb(1:2,1:num) = sx(j,k,1:2,1:num)
      ! CALL ADVDIFF_U(Deltat,ix,num,Deltax,Wind,Dif,Conc,Sb)
      CALL ADVDIFF_SUPG_MF(Deltat, 1, ix, num, X, Wind, Dif, &
                        Air, Conc, Sb, SurfaceEm, VolumeEm, Vd)
      end do ! iz
      end do ! iy
      return
      end

!*********************************************************************
      subroutine trany_mf(ix,iy,iz,nstart,nend,s1,v,kh,sy,dt,dy)
      use ParallelDataMap
!*********************************************************************
      include 'aqmax.param'
      !include 'aqindx.cmm'
      integer :: ix, iy, iz, nstart, nend
      real ::  s1(ix,iy,iz,*), v(ix,iy,*), kh(ix,iy,*)
      real ::  sy(ix,iz,2,*), dy(*), dt
      !real ::  bndf(2,mxspg),bndv(2,mxspg)
      !real ::  gen(mxgr*mxspg)
      real*8 :: Wind(1:iy), Dif(1:iy), Air(1:iy), Conc(1:iy,mxspg)
      real*8 :: Sb(1:2,mxspg), Deltay, Deltat
      integer :: num, i, j, k
      real*8 :: SurfaceEm(mxspg),VolumeEm(iy,mxspg),Vd(mxspg),Y(iy)
!      
      num=nend-nstart+1    ! number of transported species
      Deltay = dy(1)       ! grid size: this works for uniform grids only
      Deltat = dt          ! time interval

      do i=1,iy
        Y(i) = (i-1)*DeltaY
      end do	
      SurfaceEm(1:num) = 0.d0
      VolumeEm(1:iy,1:num) = 0.d0
      Vd(1:num)        = 0.d0      
      
      do i=1,no_of_yslices(MyId)   ! 1,ix
      do k=1,iz
      
      Conc(1:iy,1:num) = s1(i,1:iy,k,nstart:nend)
      Wind(1:iy) = v(i,1:iy,k)
      Dif(1:iy) = kh(i,1:iy,k)
      Air(1:iy) = s1(i,1:iy,k,iair)
      Sb(1:2,1:num) = sy(i,k,1:2,1:num)
      ! CALL ADVDIFF_U(Deltat,iy,num,Deltay,Wind,Dif,Conc,Sb)
      CALL ADVDIFF_SUPG_MF(Deltat, 1, iy, num, Y, Wind, Dif, &
                    Air, Conc, Sb, SurfaceEm, VolumeEm, Vd)
      s1(i,1:iy,k,nstart:nend) = max(Conc(1:iy,1:num),0.d0)
      
      end do ! iz
      end do ! ix
      return
      end
      
      
!*********************************************************************
      subroutine tranz_mf(ix,iy,iz,nstart,nend,s1,w,kv,q,em,vg,sz,dt,dz)
      use ParallelDataMap
!*********************************************************************
      include 'aqmax.param'
      !include 'aqindx.cmm'
!      Arguments:      
      integer :: ix, iy, iz, nstart, nend
      real :: s1(ix,iy,iz,*), w(ix,iy,iz), kv(ix,iy,iz)
      real :: q(ix,iy,*), em(ix,iy,iz,*), vg(ix,iy,*), sz(ix,iy,*)
      real :: dt, dz(ix,iy,*)
! Local variables:      
      real*8 :: Deltat, Z(iz), Wind(iz), Dif(iz), Air(iz)
      real*8 :: C(iz,mxspg), Bdry(2,mxspg) 
      real*8 :: SurfaceEm(mxspg), VolumeEm(iz,mxspg), DepVel(mxspg)     

      num = nend-nstart+1 ! number of species
      Deltat = dt         ! The time step
      Nstep = 1           ! No. of time steps
      
      do i=1,no_of_yslices(MyId)  ! i=1,ix
      do j=1,iy
      
      ! The vertical grid
      Z(1) = 0.d0
      do k = 1, iz-1
        Z(k+1) = Z(k) + dz(i,j,k)
      end do	
      
      Wind(1:iz)  = w(i,j,1:iz)
      Dif(1:iz)   = kv(i,j,1:iz)
      C(1:iz,1:num) = s1(i,j,1:iz,nstart:nend)
      Air(1:iz) = s1(i,j,1:iz,iair)
      Bdry(1,1:num) = s1(i,j,1,nstart:nend)
      Bdry(2,1:num) = sz(i,j,nstart:nend)
      SurfaceEm(1:num) = q(i,j,nstart:nend)
      VolumeEm(1:iz,1:num) = em(i,j,1:iz,nstart:nend)
      DepVel(1:num)        = vg(i,j,nstart:nend)
        
      call ADVDIFF_SUPG_MF(Deltat, 1, iz, num, Z, Wind, Dif, &
                    Air, Conc, Bdry, SurfaceEm, VolumeEm, Depvel)

      s1(i,j,1:iz,nstart:nend) = max(C(1:iz,1:num),0.d0)
		       
      end do ! j
      end do ! i
!
      return
      end

     
!*********************************************************************
      subroutine transigma(ix,iy,iz,nstart,nend,s1,w,kv,q, &
                        em,vg,sz,dt,sigma)
      use ParallelDataMap
!*********************************************************************
      include 'aqmax.param'
      !include 'aqindx.cmm'
!      Arguments:      
      integer :: ix, iy, iz, nstart, nend
      real :: s1(ix,iy,iz,*), w(ix,iy,iz), kv(ix,iy,iz)
      real :: q(ix,iy,*), em(ix,iy,iz,*), vg(ix,iy,*), sz(ix,iy,*)
      real :: dt, sigma(*)
! Local variables:      
      real*8 :: Deltat, Z(iz), Wind(iz), Dif(iz)
      real*8 :: C(iz,mxspg), Bdry(2,mxspg) 
      real*8 :: SurfaceEm(mxspg), VolumeEm(iz,mxspg), DepVel(mxspg)     

      num = nend-nstart+1 ! number of species
      Deltat = dt         ! The time step
      Nstep = 1           ! No. of time steps
      
      do i=1,no_of_yslices(MyId)  ! i=1,ix
      do j=1,iy
      
      ! The vertical grid
      Z(1:iz) = sigma(1:iz)	
      
      Wind(1:iz)  = w(i,j,1:iz)
      Dif(1:iz)   = kv(i,j,1:iz)
      C(1:iz,1:num) = s1(i,j,1:iz,nstart:nend)
      Bdry(1,1:num) = s1(i,j,1,nstart:nend)
      Bdry(2,1:num) = sz(i,j,nstart:nend)
      SurfaceEm(1:num) = q(i,j,nstart:nend)
      VolumeEm(1:iz,1:num) = em(i,j,1:iz,nstart:nend)
      DepVel(1:num)        = vg(i,j,nstart:nend)
        
      call ADVDIFF_FVZ(Deltat, Nstep, iz, num, Z, Wind, Dif, &
                    C, BDRY, SurfaceEm, VolumeEm, DepVel)   

      s1(i,j,1:iz,nstart:nend) = max(C(1:iz,1:num),0.d0)
		       
      end do ! j
      end do ! i
!
      return
      end

!-----------------------------------------------------------------------
! Timesteps the advection-diffusion equation on [0,T]
! It takes the minimal number of time steps needed to satisfy CFL<CFLMAX
! Note: this routine works ONLY for equidistant grids
!       Please check the grid is uniform before 
!-----------------------------------------------------------------------
  SUBROUTINE ADVDIFF_U(T,N,Nspec,DX,W,K,C,BDRY)
  IMPLICIT NONE
! -- Input arguments:  
      INTEGER, INTENT(IN) :: N, Nspec
      REAL*8, INTENT(IN)  :: T, W(N), K(N), BDRY(2,Nspec)
! -- Input/Output arguments:  
      REAL*8, INTENT(INOUT)  :: C(N,Nspec)
! -- Local variables:
      REAL*8, PARAMETER :: CFLMAX = 0.5d0
      REAL*8  :: AFLUX(N), S1(N), S2(N)
      INTEGER :: Istep, Nsteps, Ispec
      REAL*8  :: DA(N), DB(N), DC(N), Gamma
      REAL*8  :: DT, DX, CFL, C1(N)
!-----------------------------------------------------------------------
!       T       = Time interval for integration is [0,T]
!       N       = No. of Grid Points
!       Nspec   = No. of species 
!       X(I)    = CENTER OF GRID I (not explicitly used)
!       DX      = Length of cell I (equal for all cells)
!       C(i,j)  = CONCENTRATION of species j AT X(i) 
!                 - at time=0 on input
!                 - at time=T on output   
!       W(I)    = WIND SPEED AT X(I)
!       K(I)    = DIFFUSION COEFFICIENT AT X(I)
!       SRC(I)  = SOURCE STRENGTH AT X(I)
!
!       AFLUX(I) = TOTAL FLUX DIFFERENCE AT X(I)
!-----------------------------------------------------------------------
! --- Check that we have uniform grid


! Estimate the Courant number	
      CFL = MAXVAL(T/DX*W(1:N))
      IF ( CFL > CFLMAX ) THEN
	 NSTEPS = INT(CFL/CFLMAX)+1
	 DT = T/DBLE(NSTEPS)  
      ELSE
	 NSTEPS = 1 
	 DT = T  
      END IF 

! Build the diffusion flux Jacobian J
      CALL DIFF_FLUX_U_JAC(N,DX,W,K,DA,DB,DC)
! Build I-h*gamma*J
      Gamma = 1.0D0 + SQRT(2.0D0)/2.0D0
      DA(1:N) = 1.0D0 - (DT*Gamma)*DA(1:N)
      DB(1:N) =       - (DT*Gamma)*DB(1:N)
      DC(1:N) =       - (DT*Gamma)*DC(1:N)
! Factorize this matrix
      CALL TRI_FA(N, DA, DB, DC)
      

! Start the time loop
time:  DO Istep = 1, NSTEPS  
! Start the species loop
spec:   DO Ispec = 1, Nspec  
! ---  RK Stage 1 ---
          CALL AD_FLUX_U(N,DX,W,K,C(1,Ispec),BDRY(1,Ispec),AFLUX)
	  S1(1:N) = AFLUX(1:N)
	  S2(1:N) = S1(1:N)
	  CALL TRI_SOL(N, DA, DB, DC, S1)	  
	  C1(1:N) = C(1:N,Ispec) + DT*S1(1:N)
! ---  RK Stage 2 ---
          CALL AD_FLUX_U(N,DX,W,K,C1,BDRY(1,Ispec),AFLUX)
	  S2(1:N) = AFLUX(1:N)+2*(S2(1:N)-S1(1:N))
	  CALL TRI_SOL(N, DA, DB, DC, S2)	  
	  C1(1:N) = C1(1:N) + DT*S2(1:N)
! ---  RK Solution ---
	  C(1:N,Ispec) = ( C(1:N,Ispec) + C1(1:N) )/2
        END DO spec
       END DO time
! End of the time loop

       END SUBROUTINE ADVDIFF_U
!-----------------------------------------------------------------------
       
       
       
       
!-----------------------------------------------------------------------
!     Computes the limited advective flux plus the 
!     diffusive flux on a uniform grid
!-----------------------------------------------------------------------
  SUBROUTINE AD_FLUX_U(N,DX,W,K,C,BDRY,AFLUX)
!
  IMPLICIT NONE
!
! Arguments:
  INTEGER, INTENT(IN) :: N
  REAL*8, INTENT(IN) :: DX, W(N), K(N), C(N), BDRY(2) 
  REAL*8, INTENT(OUT) :: AFLUX(N)
! Local Variables:
  REAL*8 :: R, PHI_U, Wind, F(N), FLUX(N+1)
  REAL*8 :: Eps = 1.d-14
  INTEGER :: I

! ----------------------------      
! Input Arguments:
!       N       = No. OF Grid Points
!       X(I)    = CENTER OF GRID I, not explicitly given
!       C(I)    = CONCENTRATION  AT X(I)    
!       W(I)    = WIND SPEED AT X(I)
!       K(I)    = DIFFUSION COEF. AT X(I)
! ----------------------------
! Output Arguments:
!       AFLUX(I) = TOTAL FLUX DIFFERENCE OVER DX, CELL I
! ----------------------------
! Local Variables:
!       FLUX(I) = TOTAL FLUX AT X(I-1/2) interface
!       F(I)    = ADVECTIVE FLUX AT X(I)
! ----------------------------  

! --- The advective flux values at the cell centers
  F(1:N) = C(1:N)*W(1:N)
	
  DO I=1,N+1
!
!  The wind strength through the I-1/2 interface is:
     IF (I==1) THEN
        Wind = W(1)
     ELSEIF (I==N+1) THEN
        Wind = W(N)	
     ELSE
        Wind = ( W(I)+W(I-1) )/2.0d0
     END IF	
!   
!  Compute the total (advective + diffusive) flux through I-1/2 interface.
!
     IF ( Wind >= 0.d0 ) THEN !========================================>
!	 
          IF (I == 1) THEN       ! I=1 is inflow boundary
	    FLUX(1) = W(1)*BDRY(1)  &
	              - K(1)*( C(1)-BDRY(1) )/DX
	  ELSEIF (I == 2) THEN
	    R = ( C(2)-C(1)+Eps )/( C(2)-BDRY(1)+Eps )
	    FLUX(2) = (F(1) +  PHI_U(R)*(F(1)-BDRY(1)*W(1))) &
                      - (K(1)+K(2))/2.d0*(C(2)-C(1))/DX 
	  ELSEIF (I == N+1) THEN  ! I=N is outflow boundary
	    FLUX(N+1) = (1.5D0*F(N)-0.5D0*F(N-1)) - 0.0D0 
	  ELSE
	    R = ( C(I)-C(I-1)+Eps )/( C(I-1)-C(I-2)+Eps )
	    FLUX(I) = F(I-1) + PHI_U(R)*(F(I-1)-F(I-2)) &
	              - (K(I)+K(I-1))/2.d0*(C(I)-C(I-1))/DX
          END IF
!
     ELSE  ! Wind < 0          !========================================>
!	  
          IF (I == 1) THEN        ! I=1 is outflow boundary
	    FLUX(1) = (1.5D0*F(1)-0.5D0*F(2)) &
	              - 0.0D0
	  ELSEIF (I == N) THEN
	    R = ( C(I-1)-C(I)+Eps )/( C(I)-BDRY(2)+Eps )
	    FLUX(I) = (F(N) + PHI_U(R)*(F(N)-BDRY(2)*W(I))) &
	              - (K(N)+K(N-1))/2.d0*(C(N)-C(N-1))/DX 
	  ELSEIF (I == N+1) THEN  ! I=N is inflow boundary
	    FLUX(N+1) = W(N)*BDRY(2) &
	                - K(N)*(BDRY(2)-C(N))/DX
	  ELSE
	    R = ( C(I-1)-C(I)+Eps )/( C(I)-C(I+1)+Eps )
	    FLUX(I) = F(I) + PHI_U(R)*(F(I)-F(I+1)) &
	              - (K(I)+K(I-1))/2.d0*(C(I)-C(I-1))/DX 
          END IF 
!	  
     END IF                    !========================================> 
!	    
  END DO 
!
!     
! Compute Flux differences for each cell over DX.
! This gives the time derivative of the mean concentration in each cell.	
  DO I = 1, N
      AFLUX(I) = - ( FLUX(I+1) - FLUX(I) )/DX 
  END DO  


  END SUBROUTINE AD_FLUX_U
!-----------------------------------------------------------------------



      
!-----------------------------------------------------------------------
! The Limiter Function
      DOUBLE PRECISION FUNCTION PHI_U ( R )
      DOUBLE PRECISION, INTENT(IN) ::  R
      PHI_U = 5.0D-1*DMAX1(0.D0,DMIN1(2.D0*R, &
             DMIN1(2.D0,(1.D0+2.D0*R)/3.D0)))
      END FUNCTION PHI_U
!-----------------------------------------------------------------------



      
      
!-----------------------------------------------------------------------
!     Computes the Jacobian of the diffusive flux on a uniform grid
!-----------------------------------------------------------------------
      SUBROUTINE DIFF_FLUX_U_JAC(N,DX,W,K,DA,DB,DC)
!
      IMPLICIT NONE
!
! Input Arguments:
      INTEGER, INTENT(IN) :: N
      REAL*8, INTENT(IN) :: DX, W(N), K(N)
! Output Arguments:     
      REAL*8, INTENT(OUT) :: DA(N), DB(N), DC(N)
! Local Variables:      
      INTEGER :: I
      REAL*8 :: DF1(N+1), DF2(N+1)

!
!       N       = No. OF Grid Points
!       X(I)    = CENTER OF GRID I, does not appear explicitly
!       DX      = uniform grid spacing
!       K(I)    = Diffusion at AT X(I)
!       FLUX(I) = FLUX AT X(I-1/2)
! DA=DIAGONAL, DB=SUPRA, DC=SUB, INDEX=LINE NO. entries in the Jacobian
! DF1(i) = D_flux(i)/d_c(i-1),  DF2(i) = D_flux(i)/d_c(i)
      
!   Diffusive fluxes through the cell boundaries I-1/2
        DF1(1:N+1) = 0.0D0
	DF2(1:N+1) = 0.0D0
        DO I = 2, N
	    ! FLUX(I) = ( K(I)*C(I) - K(I-1)*C(I-1) )/DX
	    DF1(I)  = - (K(I-1)+K(I))/2.d0/DX 
	    DF2(I)  =   (K(I-1)+K(I))/2.d0/DX
	END DO  
	IF ( W(1) > 0.0D0 ) THEN ! I=1 is inflow boundary
	    ! FLUX(1) = ( K(1)*C(1) - K(1)*BDRY(1) )/DX
	    DF1(1) = 0.d0; DF2(1) = K(1)/DX
	ELSE                     ! I=1 is outflow boundary
	    ! FLUX(1) = 0.0D0
	    DF1(1) = 0.d0; DF2(1) = 0.d0
	END IF
	IF ( W(N) < 0.0D0 ) THEN ! I=N is inflow boundary
	    ! FLUX(N+1) = ( K(N)*BDRY(2) - K(N)*C(N) )/DX
	    DF1(N+1) = -K(N)/DX; DF2(N+1) = 0.d0
	ELSE
	    ! FLUX(N+1) = 0.0D0
	    DF1(N+1) = 0.0D0; DF2(N+1) = 0.0D0
	END IF
        
!  Diffusive flux differences     
        DO I = 1, N
	  ! DFLUX(I) = (FLUX(I+1) - FLUX(I))/DX
	  DA(I) = ( DF1(I+1) - DF2(I) )/DX  ! d_dflux(i)/d_c(i)
	  DB(I) = DF2(I+1)/DX               ! d_dflux(i)/d_c(i+1)
	  DC(I) = - DF1(I)/DX	            ! d_dflux(i)/d_c(i-1)  
	END DO  

      END SUBROUTINE DIFF_FLUX_U_JAC
!-----------------------------------------------------------------------



!-----------------------------------------------------------------------
subroutine ADVDIFF_SUPG(DT, Nstep, NGP, Nspec, X, U, K, C, &
                    BDRY, SurfaceEm, VolumeEm, Vd)
!-----------------------------------------------------------------------
!  Performs Nstep timesteps of length DT
!      to solve the adv_diff equation
! 
!  NGP   = no. of grid points
!  Nspec = no. of chemical species
!  Nstep = no of time steps
!  X(1:NGP) = grid point coordinates
!  U(1:NGP) = wind speeds
!  K(1:NGP) = diffusion coefficients
!  SurfaceEm  = Surface Emission intensity
!  VolumeEm   = Elevated Emission intensity
!  Vd    = deposition velocity
!  C     = concentration of each species
!-----------------------------------------------------------------------

  implicit none
  integer, intent(in) :: NGP, Nstep, Nspec
  double precision, intent(in)  :: DT, X(NGP), U(NGP), K(NGP), &
                          Vd(Nspec), SurfaceEm(Nspec), &
			  Bdry(2,Nspec), VolumeEm(NGP, Nspec)
  double precision, intent(inout) :: C(NGP,Nspec)

!  Local Variables
  ! Mass Matrix 
  double precision :: MA(NGP), MB(NGP), MC(NGP)
  ! Convection Matrix 
  double precision :: CA(NGP), CB(NGP), CC(NGP)
  ! Diffusion Matrix 
  double precision :: DA(NGP), DB(NGP), DC(NGP)
  ! Buffer Matrices
  double precision :: FA(NGP), FB(NGP), FC(NGP), &
                GA(NGP), GB(NGP), GC(NGP)
  double precision :: rhs(NGP), B(NGP), x1, x2, w1, w2
  double precision :: DU2(NGP)
  integer :: istep, ispec, info, ipiv(NGP)

! Build the matrices  
  call mass_matrix      (NGP,X,U,K,MA,MB,MC)
  call convection_matrix(NGP,X,U,K,CA,CB,CC)
  call diffusion_matrix (NGP,X,U,K,DA,DB,DC)
  
125  format(1000(E12.5,2x))
 
! Build the LHS matrix and factorize it
  GA(1:NGP) = MA(1:NGP) + DT/2.d0*(CA(1:NGP)+DA(1:NGP))
  GB(1:NGP) = MB(1:NGP) + DT/2.d0*(CB(1:NGP)+DB(1:NGP))
  GC(1:NGP) = MC(1:NGP) + DT/2.d0*(CC(1:NGP)+DC(1:NGP))
  call TRI_FA(NGP, GA, GB, GC)
  ! call DGTTRF( NGP, GC(2), GA, GB, DU2, IPIV, INFO )
 
! Build M - dt/2*(Conv+Diff)  
  FA(1:NGP) = MA(1:NGP) - DT/2.d0*(CA(1:NGP)+DA(1:NGP))
  FB(1:NGP) = MB(1:NGP) - DT/2.d0*(CB(1:NGP)+DB(1:NGP))
  FC(1:NGP) = MC(1:NGP) - DT/2.d0*(CC(1:NGP)+DC(1:NGP))
   
   
time:  do istep = 1, Nstep
! Build the rhs  
spc:   do ispec = 1,Nspec
          call MxV(NGP,FA,FB,FC,C(1,ispec),rhs)
          call boundary( NGP,X,U,K,C(1,ispec),Bdry(1,ispec), &
                 Vd(ispec),SurfaceEm(ispec),B )
          rhs = rhs + DT*B 
          ! Solve the system and advance the solution  
          call TRI_SOL( NGP, GA, GB, GC, rhs )
          ! call DGTTRS( 'N', NGP, 1, GC(2), GA, GB, DU2, &
	  !            IPIV, RHS, NGP, INFO )
          ! Add the volume emissions  
          C(1:NGP,ispec) = rhs(1:NGP) + DT*VolumeEm(1:NGP,ispec)
     end do spc
  end do time

end subroutine ADVDIFF_SUPG



subroutine ADVDIFF_SUPG_MF(DT, Nstep, NGP, Nspec, X, U, K, Air, C, &
                    BDRY, SurfaceEm, Vd)
!  Performs Nstep timesteps of length DT
!      to solve the adv_diff equation
! 
!  NGP   = no. of grid points
!  Nspec = no. of chemical species
!  Nstep = no of time steps
!  X(1:NGP) = grid point coordinates
!  U(1:NGP) = wind speeds
!  K(1:NGP) = diffusion coefficients
!  SurfaceEm  = Surface Emission intensity
!  Vd    = deposition velocity
!  C     = concentration of each species IN MOLE FRACTION
!  Air = Air Concentration

  implicit none
  integer, intent(in) :: NGP, Nstep, Nspec
  double precision, intent(in)  :: DT, X(NGP), U(NGP), K(NGP), &
                          Vd(Nspec), SurfaceEm(Nspec), &
			  Bdry(2,Nspec), Air(NGP)
  double precision, intent(inout) :: C(NGP,Nspec)

!  Local Variables
  ! Mass Matrix 
  double precision :: MA(NGP), MB(NGP), MC(NGP)
  ! Convection Matrix 
  double precision :: CA(NGP), CB(NGP), CC(NGP)
  ! Diffusion Matrix 
  double precision :: DA(NGP), DB(NGP), DC(NGP)
  ! Buffer Matrices
  double precision :: FA(NGP), FB(NGP), FC(NGP), &
                GA(NGP), GB(NGP), GC(NGP)
  double precision :: rhs(NGP), B(NGP), x1, x2, w1, w2
  integer :: istep, ispec

! Build the matrices  
  call mass_matrix         (NGP,X,U,K,MA,MB,MC)
  call convection_matrix_mf(NGP,X,U,K,CA,CB,CC)
  call diffusion_matrix_mf (NGP,X,U,K,Air,DA,DB,DC)
  
125  format(1000(E12.5,2x))
 
! Build the LHS matrix and factorize it
  GA(1:NGP) = MA(1:NGP) + DT/2.d0*(CA(1:NGP)+DA(1:NGP))
  GB(1:NGP) = MB(1:NGP) + DT/2.d0*(CB(1:NGP)+DB(1:NGP))
  GC(1:NGP) = MC(1:NGP) + DT/2.d0*(CC(1:NGP)+DC(1:NGP))
  call TRI_FA(NGP, GA, GB, GC)
 
! Build M - dt/2*(Conv+Diff)  
  FA(1:NGP) = MA(1:NGP) - DT/2.d0*(CA(1:NGP)+DA(1:NGP))
  FB(1:NGP) = MB(1:NGP) - DT/2.d0*(CB(1:NGP)+DB(1:NGP))
  FC(1:NGP) = MC(1:NGP) - DT/2.d0*(CC(1:NGP)+DC(1:NGP))
   
   
time:  do istep = 1, Nstep
! Build the rhs  
spc:   do ispec = 1,Nspec
          call MxV(NGP,FA,FB,FC,C(1,ispec),rhs)
          call boundary(NGP,X,U,K,C(1,ispec),Bdry(1,ispec), &
                 Vd(ispec),SurfaceEm(ispec),B)
          rhs = rhs + DT*B
          ! Solve the system and advance the solution  
          call TRI_SOL(NGP, GA, GB, GC, rhs)
          C(1:NGP,ispec) = rhs(1:NGP)
     end do spc
  end do time

end subroutine ADVDIFF_SUPG_MF


subroutine MxV(N,A,B,C,X,Y)
!  Matrix A times vector X equals vector Y
! The matrix is tridiagonal, given by 
! A=DIAGONAL, B=SUPRA, C=SUB, INDEX=LINE NO.
  implicit none
  integer, intent(in) :: N
  double precision, intent(in)  :: A(N), B(N), C(N), X(N) 
  double precision, intent(out) :: Y(N)
  integer :: i

  Y(1) = A(1)*X(1)+B(1)*X(2)
  do i=2,N-1
     Y(i) = C(i)*X(i-1) + A(i)*X(i) + B(i)*X(i+1)
  end do   
  Y(N) = C(N)*X(N-1)+A(N)*X(N)

end subroutine MxV


subroutine boundary(N,X,U,K,C,Bdry,Vd,E,B)
!
!  U  =  wind speed
! Vd  =  deposition velocity
! C   =  concentration
! Cabove = concentration above the domain upper boundary
! B_j =  the value of (K dc/dz)(z_j) csi_j
! Bdry = boundary values (for Dirichlet condition)

  implicit none
  
  integer, intent(in) :: N
  double precision, intent(in)  :: X(N), U(N), K(N), &
                           C(N), Vd, E, Bdry(2)
  double precision, intent(out) :: B(N)
  double precision, external :: xi

   B(1:N) = 0.d0
! Bottom    
   ! B(1) = - K dc/dz at z(1)
   !     = (-U(1)*(Bdry(1)-C(1)) - Vd*C(1) + E)*xi(N,X,U,K,X(1),1)
   if ( U(1)>0.d0 ) then
     B(1) = ( U(1)*(Bdry(1)-C(1)) - Vd*C(1) + E )*xi(N,X,U,K,X(1),1)
   else
     B(1) = ( - Vd*C(1) + E )*xi(N,X,U,K,X(1),1)
   end if    
! Top of the domain
   if ( U(N)<0 ) then ! inflow
      B(N) = -U(N)*(Bdry(2) - C(N))*xi(N,X,U,K,X(N),N)
   else ! outflow
      B(N) = 0.d0
   end if      
  
end subroutine boundary


! The convection matrix
subroutine diffusion_matrix(N,X,U,K,DA,DB,DC)
! Builds Convection matrix D_{ij} = integral K d(phi_j)/dz d(xi_i)/dz 
! A=DIAGONAL, B=SUPRA, C=SUB, INDEX=LINE NO.
!
!  X(N) = grid points
! K(N) = diffusion coefficient

  implicit none
  
  integer, intent(in) :: N
  double precision, intent(in)  :: X(N), U(N), K(N)
  double precision, intent(out) :: DA(N), DB(N), DC(N)
  
  integer :: i
  double precision :: xa, xb, xc, xd, wa, wb, wc, wd
  double precision, external :: phi, phi_prime, xi, xi_prime

! Gauss Points 

  call gauss( X(1),X(2),xc,xd,wc,wd )
  ! integral phi_1 xi_1 from x_{1} to x_{2}
  DA(1) =   wc*phi_prime(N,X,xc,1)*xi_prime(N,X,U,K,xc,1)  & 
              *(K(1)*phi(N,X,xc,1)+K(2)*phi(N,X,xc,2)) & 
	  + wd*phi_prime(N,X,xd,1)*xi_prime(N,X,U,K,xd,1)  &
              *(K(1)*phi(N,X,xd,1)+K(2)*phi(N,X,xd,2))  
  ! xi_{1} phi_{2}	    
  DB(1) =   wc*phi_prime(N,X,xc,2)*xi_prime(N,X,U,K,xc,1)  & 
              *(K(1)*phi(N,X,xc,1)+K(2)*phi(N,X,xc,2)) & 
	  + wd*phi_prime(N,X,xd,2)*xi_prime(N,X,U,K,xd,1)  &
              *(K(1)*phi(N,X,xd,1)+K(2)*phi(N,X,xd,2))  

  do i=2,N-1
    call gauss( X(i-1),X(i),xa,xb,wa,wb )
    call gauss( X(i),X(i+1),xc,xd,wc,wd )
    ! integral phi_j xi_i from x_{i-1} to x_{i+1}
    DA(i) =   wa*phi_prime(N,X,xa,i)*xi_prime(N,X,U,K,xa,i)      & 
                *(K(i)*phi(N,X,xa,i-1)+K(i)*phi(N,X,xa,i))   & 
            + wb*phi_prime(N,X,xb,i)*xi_prime(N,X,U,K,xb,i)      & 
                *(K(i)*phi(N,X,xb,i-1)+K(i)*phi(N,X,xb,i))   & 
	    + wc*phi_prime(N,X,xc,i)*xi_prime(N,X,U,K,xc,i)      & 
                *(K(i)*phi(N,X,xc,i)+K(i+1)*phi(N,X,xc,i+1)) & 
	    + wd*phi_prime(N,X,xd,i)*xi_prime(N,X,U,K,xd,i)      &
                *(K(i)*phi(N,X,xd,i)+K(i+1)*phi(N,X,xd,i+1)) 
    ! xi_{i} phi_{i+1}	    
    DB(i) =   wc*phi_prime(N,X,xc,i+1)*xi_prime(N,X,U,K,xc,i)    & 
                *(K(i)*phi(N,X,xc,i)+K(i+1)*phi(N,X,xc,i+1)) & 
	    + wd*phi_prime(N,X,xd,i+1)*xi_prime(N,X,U,K,xd,i)    &
                *(K(i)*phi(N,X,xd,i)+K(i+1)*phi(N,X,xd,i+1))  
    ! xi_{i} phi_{i-1}	    
    DC(i) =   wa*phi_prime(N,X,xa,i-1)*xi_prime(N,X,U,K,xa,i)    & 
                *(K(i)*phi(N,X,xa,i-1)+K(i)*phi(N,X,xa,i))   & 
            + wb*phi_prime(N,X,xb,i-1)*xi_prime(N,X,U,K,xb,i)    &
                *(K(i-1)*phi(N,X,xb,i-1)+K(i)*phi(N,X,xb,i))  
  end do

  call gauss( X(N-1),X(N),xa,xb,wa,wb )
  ! integral phi_j xi_i from x_{i-1} to x_{i+1}
  DA(N) =   wa*phi_prime(N,X,xa,N)*xi_prime(N,X,U,K,xa,N)        & 
              *(K(N-1)*phi(N,X,xa,N-1)+K(N)*phi(N,X,xa,N))   & 
          + wb*phi_prime(N,X,xb,N)*xi_prime(N,X,U,K,xb,N)        & 
              *(K(N-1)*phi(N,X,xb,N-1)+K(N)*phi(N,X,xb,N))  
  ! xi_{i} phi_{i-1}	    
  DC(N) =   wa*phi_prime(N,X,xa,N-1)*xi_prime(N,X,U,K,xa,N)      & 
              *(K(N-1)*phi(N,X,xa,N-1)+K(N)*phi(N,X,xa,N))   & 
          + wb*phi_prime(N,X,xb,N-1)*xi_prime(N,X,U,K,xb,N)      &
              *(K(N-1)*phi(N,X,xb,N-1)+K(N)*phi(N,X,xb,N))    

end subroutine diffusion_matrix



! The convection matrix
subroutine convection_matrix(N,X,U,K,CA,CB,CC)
! Builds Convection matrix C_{ij} = integral d(phi_j w)/dz xi_i 
! A=DIAGONAL, B=SUPRA, C=SUB, INDEX=LINE NO.
  implicit none
  
  integer, intent(in) :: N
  double precision, intent(in)  :: X(N), U(N), K(N)
  double precision, intent(out) :: CA(N), CB(N), CC(N)
  
  integer :: i
  double precision :: xa, xb, xc, xd, wa, wb, wc, wd
  double precision, external :: phi, phi_prime, xi, xi_prime

! Gauss Points 

  call gauss( X(1),X(2),xc,xd,wc,wd )
  ! integral phi_1 xi_1 from x_{1} to x_{2}
  CA(1) =   wc*xi(N,X,U,K,xc,1)                                 & 
              *( U(2)*( phi_prime(N,X,xc,2)*phi(N,X,xc,1)   &
		+phi_prime(N,X,xc,1)*phi(N,X,xc,2) )        &
	      +2.d0*U(1)*phi(N,X,xc,1)*phi_prime(N,X,xc,1) )&
	  + wd*xi(N,X,U,K,xd,1)                                 &
              *( U(2)*( phi_prime(N,X,xd,2)*phi(N,X,xd,1)   &
		+phi_prime(N,X,xd,1)*phi(N,X,xd,2) )        &
	      +2.d0*U(1)*phi(N,X,xd,1)*phi_prime(N,X,xd,1) )
  ! xi_{1} phi_{2}	    
  CB(1) =   wc*xi(N,X,U,K,xc,1)                                 & 
              *(U(i)*( phi_prime(N,X,xc,1)*phi(N,X,xc,2)    &
	      + phi(N,X,xc,1)*phi_prime(N,X,xc,2) )         & 
	      + 2.d0*U(2)*phi_prime(N,X,xc,2)*phi(N,X,xc,2))& 
	  + wd*xi(N,X,U,K,xd,1)                                 &
              *(U(1)*( phi_prime(N,X,xd,1)*phi(N,X,xd,2)    &
	      + phi(N,X,xd,1)*phi_prime(N,X,xd,2) )         & 
	      + 2.d0*U(2)*phi_prime(N,X,xd,2)*phi(N,X,xd,2))  

  do i=2,N-1
    call gauss( X(i-1),X(i),xa,xb,wa,wb )
    call gauss( X(i),X(i+1),xc,xd,wc,wd )
    ! integral phi_j xi_i from x_{i-1} to x_{i+1}
    CA(i) =   wa*xi(N,X,U,K,xa,i)                                   & 
                *( U(i-1)*( phi_prime(N,X,xa,i-1)*phi(N,X,xa,i) &
		  +phi_prime(N,X,xa,i)*phi(N,X,xa,i-1) )        &
		  +2.d0*U(i)*phi(N,X,xa,i)*phi_prime(N,X,xa,i) )&
            + wb*xi(N,X,U,K,xb,i)                                   & 
                *( U(i-1)*( phi_prime(N,X,xb,i-1)*phi(N,X,xb,i) &
		  +phi_prime(N,X,xb,i)*phi(N,X,xb,i-1) )        &
		  +2.d0*U(i)*phi(N,X,xb,i)*phi_prime(N,X,xb,i) )&
	    + wc*xi(N,X,U,K,xc,i)                                   & 
                *( U(i+1)*( phi_prime(N,X,xc,i+1)*phi(N,X,xc,i) &
		  +phi_prime(N,X,xc,i)*phi(N,X,xc,i+1) )        &
		  +2.d0*U(i)*phi(N,X,xc,i)*phi_prime(N,X,xc,i) )&
	    + wd*xi(N,X,U,K,xd,i)                                  &
                *( U(i+1)*( phi_prime(N,X,xd,i+1)*phi(N,X,xd,i) &
		  +phi_prime(N,X,xd,i)*phi(N,X,xd,i+1) )        &
		  +2.d0*U(i)*phi(N,X,xd,i)*phi_prime(N,X,xd,i) )
    ! xi_{i} phi_{i+1}	    
    CB(i) =   wc*xi(N,X,U,K,xc,i)                                      & 
                *(U(i)*( phi_prime(N,X,xc,i)*phi(N,X,xc,i+1)         &
		+ phi(N,X,xc,i)*phi_prime(N,X,xc,i+1) )              & 
		+ 2.d0*U(i+1)*phi_prime(N,X,xc,i+1)*phi(N,X,xc,i+1)) & 
	    + wd*xi(N,X,U,K,xd,i)                      &
                *(U(i)*( phi_prime(N,X,xd,i)*phi(N,X,xd,i+1)         &
		+ phi(N,X,xd,i)*phi_prime(N,X,xd,i+1) )              & 
		+ 2.d0*U(i+1)*phi_prime(N,X,xd,i+1)*phi(N,X,xd,i+1))  
    ! xi_{i} phi_{i-1}	    
    CC(i) =   wa*xi(N,X,U,K,xa,i)                                      & 
                *(U(i)*( phi_prime(N,X,xa,i)*phi(N,X,xa,i-1)         &
		+ phi(N,X,xa,i)*phi_prime(N,X,xa,i-1) )              & 
		+ 2.d0*U(i-1)*phi_prime(N,X,xa,i-1)*phi(N,X,xa,i-1)) & 
            + wb*xi(N,X,U,K,xb,i)                                      & 
                *(U(i)*( phi_prime(N,X,xb,i)*phi(N,X,xb,i-1)         &
		+ phi(N,X,xb,i)*phi_prime(N,X,xb,i-1) )              & 
		+ 2.d0*U(i-1)*phi_prime(N,X,xb,i-1)*phi(N,X,xb,i-1))  
  end do

  call gauss( X(N-1),X(N),xa,xb,wa,wb )
  ! integral phi_j xi_i from x_{i-1} to x_{i+1}
  CA(N) =   wa*xi(N,X,U,K,xa,N)                                   & 
                *( U(N-1)*( phi_prime(N,X,xa,N-1)*phi(N,X,xa,N) &
		  +phi_prime(N,X,xa,N)*phi(N,X,xa,N-1) )        &
		  +2.d0*U(N)*phi(N,X,xa,N)*phi_prime(N,X,xa,N) )&
          + wb*xi(N,X,U,K,xb,N)                                   & 
                *( U(N-1)*( phi_prime(N,X,xb,N-1)*phi(N,X,xb,N) &
		  +phi_prime(N,X,xb,N)*phi(N,X,xb,N-1) )        &
		  +2.d0*U(N)*phi(N,X,xb,N)*phi_prime(N,X,xb,N) )
  ! xi_{i} phi_{i-1}	    
  CC(N) =   wa*xi(N,X,U,K,xa,N)                                        & 
                *(U(N)*( phi_prime(N,X,xa,N)*phi(N,X,xa,N-1)         &
		+ phi(N,X,xa,N)*phi_prime(N,X,xa,N-1) )              & 
		+ 2.d0*U(N-1)*phi_prime(N,X,xa,N-1)*phi(N,X,xa,N-1)) & 
          + wb*xi(N,X,U,K,xb,N)                                        &
                *(U(N)*( phi_prime(N,X,xb,N)*phi(N,X,xb,N-1)         &
		+ phi(N,X,xb,N)*phi_prime(N,X,xb,N-1) )              & 
		+ 2.d0*U(N-1)*phi_prime(N,X,xb,N-1)*phi(N,X,xb,N-1))  

end subroutine convection_matrix



! The convection matrix in mole fraction formulation
subroutine convection_matrix_mf(N,X,U,K,CA,CB,CC)
! Builds Convection matrix C_{ij} = integral d(phi_j w)/dz xi_i 
! A=DIAGONAL, B=SUPRA, C=SUB, INDEX=LINE NO.
  implicit none
  
  integer, intent(in) :: N
  double precision, intent(in)  :: X(N), U(N), K(N)
  double precision, intent(out) :: CA(N), CB(N), CC(N)
  
  integer :: i
  double precision :: xa, xb, xc, xd, wa, wb, wc, wd
  double precision, external :: phi, phi_prime, xi, xi_prime


  do i=1,N
    ! Gauss Points 
    if ( i>=2 ) then   
      call gauss( X(i-1),X(i),xa,xb,wa,wb )
    end if
    if ( i<= N-1 ) then  
      call gauss( X(i),X(i+1),xc,xd,wc,wd )
    end if  
    ! integral phi_j xi_i from x_{i-1} to x_{i+1}
    CA(i) = 0.d0
    if ( i>=2 ) then 
        CA(i) = CA(i) &
            + wa*xi(N,X,U,K,xa,i)*phi_prime(N,X,xa,i)           & 
                *( U(i-1)*phi(N,X,xa,i-1)+U(i)*phi(N,X,xa,i) )  &
            + wb*xi(N,X,U,K,xb,i)*phi_prime(N,X,xb,i)           &                          
                *( U(i-1)*phi(N,X,xb,i-1)+U(i)*phi(N,X,xb,i) )  
    end if
    if ( i<= N-1 ) then	     
        CA(i) = CA(i) &
	    + wc*xi(N,X,U,K,xc,i)*phi_prime(N,X,xc,i)           &       
                *( U(i+1)*phi(N,X,xc,i+1)+U(i)*phi(N,X,xc,i) )  &
	    + wd*xi(N,X,U,K,xd,i)*phi_prime(N,X,xd,i)           &                       
                *( U(i+1)*phi(N,X,xd,i+1)+U(i)*phi(N,X,xd,i) ) 
    end if
    ! xi_{i} phi_{i+1}	    
    if ( i<=N-1 ) then
      CB(i) = wc*xi(N,X,U,K,xc,i)*phi_prime(N,X,xc,i+1)         & 
                *( U(i)*phi(N,X,xc,i)+U(i+1)*phi(N,X,xc,i+1) )  & 
	    + wd*xi(N,X,U,K,xd,i)*phi_prime(N,X,xd,i+1)         &
                *( U(i)*phi(N,X,xd,i)+U(i+1)*phi(N,X,xd,i+1) )
    else
      CB(i) = 0.0d0
    end if   		
    ! xi_{i} phi_{i-1}	
    if ( i>=2 ) then   
      CC(i) = wa*xi(N,X,U,K,xa,i)*phi_prime(N,X,xa,i-1)      & 
                *(U(i)*phi(N,X,xa,i)+U(i-1)*phi(N,X,xa,i-1)) & 
            + wb*xi(N,X,U,K,xb,i)*phi_prime(N,X,xb,i-1)      & 
                *(U(i)*phi(N,X,xb,i)+U(i-1)*phi(N,X,xb,i-1))  
    else
      CC(i) = 0.0d0
    end if   		
  end do


end subroutine convection_matrix_mf


! The diffusion matrix in mole-fraction formulation
subroutine diffusion_matrix_mf(N,X,U,K,Air,DA,DB,DC)
! Builds Convection matrix D_{ij} = integral K d(phi_j)/dz d(xi_i)/dz 
! A=DIAGONAL, B=SUPRA, C=SUB, INDEX=LINE NO.
!
!  X(N) = grid points
!  K(N) = diffusion coefficient
!  Air = air density

  implicit none
  
  integer, intent(in) :: N
  double precision, intent(in)  :: X(N), U(N), K(N), Air(N)
  double precision, intent(out) :: DA(N), DB(N), DC(N)
  
  integer :: i
  double precision :: xa, xb, xc, xd, wa, wb, wc, wd
  double precision, external :: phi, phi_prime, xi, xi_prime

  do i=1,N
    ! Gauss Points 
    if ( i>=2 ) then   
      call gauss( X(i-1),X(i),xa,xb,wa,wb )
    end if
    if ( i<= N-1 ) then  
      call gauss( X(i),X(i+1),xc,xd,wc,wd )
    end if  
    ! integral phi_j xi_i from x_{i-1} to x_{i+1}
    DA(i) = 0.d0
    if ( i>=2 ) then 
        DA(i) = DA(i) &
	   + wa*phi_prime(N,X,xa,i)                                & 
                *( Air(i)*K(i)*phi(N,X,xa,i)+Air(i-1)*K(i-1)*phi(N,X,xa,i-1) )    &
		*( xi_prime(N,X,U,K,xa,i)*( 1.d0/Air(i)*phi(N,X,xa,i) &
		      +1.d0/Air(i-1)*phi(N,X,xa,i-1) ) +              &
		   xi(N,X,U,K,xa,i)*( 1.d0/Air(i)*phi_prime(N,X,xa,i)& 
		      +1.d0/Air(i-1)*phi_prime(N,X,xa,i-1) ) )        & 
	    + wb*phi_prime(N,X,xb,i)                                & 
                *( Air(i)*K(i)*phi(N,X,xb,i)+Air(i-1)*K(i-1)*phi(N,X,xb,i-1) )    &
		*( xi_prime(N,X,U,K,xb,i)*( 1.d0/Air(i)*phi(N,X,xb,i) &
		      +1.d0/Air(i-1)*phi(N,X,xb,i-1) ) +              &
		   xi(N,X,U,K,xb,i)*( 1.d0/Air(i)*phi_prime(N,X,xb,i)& 
		      +1.d0/Air(i-1)*phi_prime(N,X,xb,i-1) ) ) 
    end if
    if ( i<= N-1 ) then	     
        DA(i) = DA(i) &
	    + wc*phi_prime(N,X,xc,i)                                & 
                *( Air(i)*K(i)*phi(N,X,xc,i)+Air(i+1)*K(i+1)*phi(N,X,xc,i+1) )    &
		*( xi_prime(N,X,U,K,xc,i)*( 1.d0/Air(i)*phi(N,X,xc,i) &
		      +1.d0/Air(i+1)*phi(N,X,xc,i+1) ) +              &
		   xi(N,X,U,K,xc,i)*( 1.d0/Air(i)*phi_prime(N,X,xc,i)& 
		      +1.d0/Air(i+1)*phi_prime(N,X,xc,i+1) ) )        & 
	    + wd*phi_prime(N,X,xd,i)                                & 
                *( Air(i)*K(i)*phi(N,X,xd,i)+Air(i+1)*K(i+1)*phi(N,X,xd,i+1) )    &
		*( xi_prime(N,X,U,K,xd,i)*( 1.d0/Air(i)*phi(N,X,xd,i) &
		      +1.d0/Air(i+1)*phi(N,X,xd,i+1) ) +              &
		   xi(N,X,U,K,xd,i)*( 1.d0/Air(i)*phi_prime(N,X,xd,i)& 
		      +1.d0/Air(i+1)*phi_prime(N,X,xd,i+1) ) )
    ! xi_{i} phi_{i+1}	    
    end if
    if ( i<=N-1 ) then
      DB(i) = wc*phi_prime(N,X,xc,i+1)                                & 
                *( Air(i)*K(i)*phi(N,X,xc,i)+Air(i+1)*K(i+1)*phi(N,X,xc,i+1) )    &
		*( xi_prime(N,X,U,K,xc,i)*( 1.d0/Air(i)*phi(N,X,xc,i) &
		      +1.d0/Air(i+1)*phi(N,X,xc,i+1) ) +              &
		   xi(N,X,U,K,xc,i)*( 1.d0/Air(i)*phi_prime(N,X,xc,i)& 
		      +1.d0/Air(i+1)*phi_prime(N,X,xc,i+1) ) )        & 
	    + wd*phi_prime(N,X,xd,i+1)                                & 
                *( Air(i)*K(i)*phi(N,X,xd,i)+Air(i+1)*K(i+1)*phi(N,X,xd,i+1) )    &
		*( xi_prime(N,X,U,K,xd,i)*( 1.d0/Air(i)*phi(N,X,xd,i) &
		      +1.d0/Air(i+1)*phi(N,X,xd,i+1) ) +              &
		   xi(N,X,U,K,xd,i)*( 1.d0/Air(i)*phi_prime(N,X,xd,i)& 
		      +1.d0/Air(i+1)*phi_prime(N,X,xd,i+1) ) )
    else
      DB(i) = 0.0d0
    end if   		
    ! xi_{i} phi_{i-1}	
    if ( i>=2 ) then   
      DC(i) = wa*phi_prime(N,X,xa,i-1)                                & 
                *( Air(i)*K(i)*phi(N,X,xa,i)+Air(i-1)*K(i-1)*phi(N,X,xa,i-1) )    &
		*( xi_prime(N,X,U,K,xa,i)*( 1.d0/Air(i)*phi(N,X,xa,i) &
		      +1.d0/Air(i-1)*phi(N,X,xa,i-1) ) +              &
		   xi(N,X,U,K,xa,i)*( 1.d0/Air(i)*phi_prime(N,X,xa,i)& 
		      +1.d0/Air(i-1)*phi_prime(N,X,xa,i-1) ) )        & 
	    + wb*phi_prime(N,X,xb,i-1)                                & 
                *( Air(i)*K(i)*phi(N,X,xb,i)+Air(i-1)*K(i-1)*phi(N,X,xb,i-1) )    &
		*( xi_prime(N,X,U,K,xb,i)*( 1.d0/Air(i)*phi(N,X,xb,i) &
		      +1.d0/Air(i-1)*phi(N,X,xb,i-1) ) +              &
		   xi(N,X,U,K,xb,i)*( 1.d0/Air(i)*phi_prime(N,X,xb,i)& 
		      +1.d0/Air(i-1)*phi_prime(N,X,xb,i-1) ) ) 
    else
      DC(i) = 0.0d0
    end if   		
  end do

end subroutine diffusion_matrix_mf



! The mass matrix
subroutine mass_matrix(N,X,U,K,MA,MB,MC)
! Builds mass matrix A_{ij} = integral phi_j xi_i 
! A=DIAGONAL, B=SUPRA, C=SUB, INDEX=LINE NO.
  implicit none
  
  integer, intent(in) :: N
  double precision, intent(in)  :: X(N), U(N), K(N)
  double precision, intent(out) :: MA(N), MB(N), MC(N)
  
  integer :: i
  double precision :: xa, xb, xc, xd, wa, wb, wc, wd
  double precision, external :: phi, xi

! Gauss Points 

  call gauss( X(1),X(2),xc,xd,wc,wd )
  ! integral phi_1 xi_1 from x_{1} to x_{2}
  MA(1) =   wc*phi(N,X,xc,1)*xi(N,X,U,K,xc,1) & 
	  + wd*phi(N,X,xd,1)*xi(N,X,U,K,xd,1)
  ! xi_{1} phi_{2}	    
  MB(1) =   wc*phi(N,X,xc,2)*xi(N,X,U,K,xc,1) & 
	  + wd*phi(N,X,xd,2)*xi(N,X,U,K,xd,1)
  MC(1) = 0.d0	  

  do i=2,N-1
    call gauss( X(i-1),X(i),xa,xb,wa,wb )
    call gauss( X(i),X(i+1),xc,xd,wc,wd )
    ! integral phi_j xi_i from x_{i-1} to x_{i+1}
    MA(i) =   wa*phi(N,X,xa,i)*xi(N,X,U,K,xa,i) & 
            + wb*phi(N,X,xb,i)*xi(N,X,U,K,xb,i) & 
	    + wc*phi(N,X,xc,i)*xi(N,X,U,K,xc,i) & 
	    + wd*phi(N,X,xd,i)*xi(N,X,U,K,xd,i)
    ! xi_{i} phi_{i+1}	    
    MB(i) =   wc*phi(N,X,xc,i+1)*xi(N,X,U,K,xc,i) & 
	    + wd*phi(N,X,xd,i+1)*xi(N,X,U,K,xd,i)
    ! xi_{i} phi_{i-1}	    
    MC(i) =   wa*phi(N,X,xa,i-1)*xi(N,X,U,K,xa,i) & 
            + wb*phi(N,X,xb,i-1)*xi(N,X,U,K,xb,i)
  end do

  call gauss( X(N-1),X(N),xa,xb,wa,wb )
  ! integral phi_j xi_i from x_{i-1} to x_{i+1}
  MA(N) =   wa*phi(N,X,xa,N)*xi(N,X,U,K,xa,N) & 
          + wb*phi(N,X,xb,N)*xi(N,X,U,K,xb,N)
  MB(N) = 0.d0	   
  ! xi_{i} phi_{i-1}	    
  MC(N) =   wa*phi(N,X,xa,N-1)*xi(N,X,U,K,xa,N) & 
          + wb*phi(N,X,xb,N-1)*xi(N,X,U,K,xb,N)

end subroutine mass_matrix

! Gauss quadrature points


! The gauss points x1,x2 and weights w1,w2 on [a,b]
subroutine gauss(a,b,x1,x2,w1,w2)
implicit none
double precision, intent(in) :: a,b
double precision, intent(out) :: x1,x2,w1,w2
double precision, parameter :: t1 = 0.78867513459481d0
double precision, parameter :: t2 = 0.21132486540519d0
double precision, parameter :: sqrt3 = 1.73205080756888d0
  w1 = (b-a)/2.d0
  w2 = w1
  x1 = t1*(2*b+a-sqrt3*b)
  x2 = t2*(2*b+a+sqrt3*b)
end subroutine gauss




! The basis function
double precision function phi(N,X,y,i)
  implicit none
  integer :: N, i
  double precision :: X(N), y
  
  if (i==1) then
     if  ( (y>=X(1)).and.(y<=X(2)) ) then
        phi = (X(2)-y)/(X(2)-X(1))
     else
        phi = 0.d0
     end if
  else if (i==N) then
     if  ( (y>X(N-1)).and.(y<=X(N)) ) then
        phi = (y-X(N-1))/(X(N)-X(N-1))
     else
        phi = 0.d0
     end if
  else
     if     ( (y>X(i-1)).and.(y<=X(i)) ) then
        phi = (y-X(i-1))/(X(i)-X(i-1))
     else if ( (y>=X(i))  .and.(y<X(i+1)) ) then
        phi = (X(i+1)-y)/(X(i+1)-X(i))
     else
        phi = 0.d0
     end if
  end if

end function phi


! The basis function
double precision function phi_prime(N,X,y,i)
  implicit none
  integer :: N, i
  double precision :: X(N), y
  
  if (i==1) then
     if  ( (y>=X(1)).and.(y<=X(2)) ) then
        phi_prime = -1.d0/(X(2)-X(1))
     else
        phi_prime = 0.d0
     end if
  else if (i==N) then
     if  ( (y>X(N-1)).and.(y<=X(N)) ) then
        phi_prime = 1.d0/(X(N)-X(N-1))
     else
        phi_prime = 0.d0
     end if
  else
     if     ( (y>=X(i-1)).and.(y<=X(i)) ) then
        phi_prime = 1.d0/(X(i)-X(i-1))
     else if ( (y>X(i))  .and.(y<X(i+1)) ) then
        phi_prime = -1.d0/(X(i+1)-X(i))
     else
        phi_prime = 0.d0
     end if
  end if

end function phi_prime



! The basis function
double precision function xi(N,X,U,K,y,i)
  implicit none
  integer, intent(IN) :: N, i
  double precision, intent(IN) :: X(N), U(N), K(N), y
  double precision, external :: phi, phi_prime, uc, uc_prime
  double precision :: dx, alpha, zeta
  
  if ( i==1 ) then
     dx = X(2)-X(1)
  else if ( i==N ) then
     dx = X(N)-X(N-1)
  else
     dx = ( X(i+1)-X(i-1) )/2.d0
  end if
  
!  if ( K(i) > 0.d0 ) then
  if ( .false. ) then
    alpha = U(i)*dx/(2.d0*K(i))  
    ! Critical approximation  
    if ( alpha <= -1.d0 ) then
      zeta = -1.d0 - 1.d0/alpha
    else if ( alpha >= 1.d0 ) then
      zeta = 1.d0 - 1.d0/alpha
    else
      zeta = 0.d0
    end if  
  else
    if ( U(i) <= 0.d0 ) then
      zeta = -1.d0
    else 
      zeta = 1.d0
    end if  
  end if               

! Classical SUPG/Discontinuous
    xi = phi( N,X,y,i ) + (dx*zeta/2)*phi_prime ( N,X,y,i )
! Christie SUPG/Continuous
!    xi = phi( N,X,y,i ) + (dx*zeta/2)*uc( N,X,y,i )
  
end function xi


! The basis function
double precision function xi_prime(N,X,U,K,y,i)
  implicit none
  integer, intent(IN) :: N, i
  double precision, intent(IN) :: X(N), U(N), K(N), y
  double precision, external :: phi, phi_prime, uc, uc_prime
  double precision :: dx, alpha, zeta
     
  if ( i==1 ) then
     dx = X(2)-X(1)
  else if ( i==N ) then
     dx = X(N)-X(N-1)
  else
     dx = ( X(i+1)-X(i-1) )/2.d0
  end if
  
  if ( K(i) > 0.d0 ) then
    alpha = U(i)*dx/(2.d0*K(i))  
    ! Critical approximation  
    if ( alpha <= -1.d0 ) then
      zeta = -1.d0 - 1.d0/alpha
    else if ( alpha >= 1.d0 ) then
      zeta = 1.d0 - 1.d0/alpha
    else
      zeta = 0.d0
    end if  
  else
    if ( U(i) <= 0.d0 ) then
      zeta = -1.d0
    else 
      zeta = 1.d0
    end if  
  end if               

! Classical SUPG/Discontinuous
  xi_prime = phi_prime ( N,X,y,i )
! Christie SUPG/Continuous
!  xi_prime = phi_prime( N,X,y,i ) + (dx*zeta/2)*uc_prime( N,X,y,i )

end function xi_prime


! The upwind correction function
double precision function uc(N,X,y,i)
  implicit none
  integer, intent(IN) :: N, i
  double precision, intent(IN) :: X(N), y
   
  if ( (y>=X(i-1)).and.(y<X(i)) ) then 
      uc = -(y-X(i-1))*(y-X(i))
  else if ( (y>=X(i)).and.(y<X(i+1)) ) then
      uc = +(y-X(i))*(y-X(i+1))  
  else
      uc = 0.d0    
  end if
  
end function uc

! The upwind correction function
double precision function uc_prime(N,X,y,i)
  implicit none
  integer, intent(IN) :: N, i
  double precision, intent(IN) :: X(N), y
   
  if ( (y>=X(i-1)).and.(y<X(i)) ) then 
      uc_prime = -(y-X(i-1))-(y-X(i))
  else if ( (y>=X(i)).and.(y<X(i+1)) ) then
      uc_prime = +(y-X(i))+(y-X(i+1))  
  else
      uc_prime = 0.d0    
  end if
  

end function uc_prime



SUBROUTINE ADVDIFF_FVZ(DT, Nstep, N, Nspec, Z, W, K, C, &
                    BDRY, SurfaceEm, VolumeEm, Vd)
!  Performs Nstep timesteps of length DT
!      to solve the adv_diff equation in vertical direction
!      using finite volume method
! 
!  N     = no. of grid points
!  Nspec = no. of chemical species
!  Nstep = no of time steps
!  Z(1:N) = grid point coordinates
!  W(1:N) = wind speeds
!  K(1:N) = diffusion coefficients
!  SurfaceEm  = Surface Emission intensity
!  VolumeEm   = Elevated Emission intensity
!  Vd    = deposition velocity
!  C     = concentration of each species
!
! Note: it uses Ros-2 with positive implementation
!
  implicit none
  integer, intent(in) :: N, Nstep, Nspec
  double precision, intent(in)  :: DT, Z(N), W(N), K(N), &
                          Vd(Nspec), SurfaceEm(Nspec), &
			  Bdry(2,Nspec), VolumeEm(N, Nspec)
  double precision, intent(inout) :: C(N,Nspec)

!  Local Variables
  ! Jacobian Matrix 
  double precision :: DCL(N), DCD(N), DCU(N)
  double precision :: S1(N), S2(N)
  double precision :: C1(N), DU2(N), gam
  integer :: istep, ispec, info, ipiv(N)

  gam = 1 + sqrt(2.d0)/2.d0
!  The Jacobian
  call ADVDIFF_JAC_FVZ(N,Z,W,K,Vd,C,DCL,DCD,DCU)  
!  I - DT*J
  DCL(1:N) =      - DT*gam*DCL(1:N)
  DCD(1:N) = 1.d0 - DT*gam*DCD(1:N)
  DCU(1:N) =      - DT*gam*DCU(1:N)
  call TRI_FA(N, DCD, DCU, DCL)
  ! call DGTTRF( N, DCL(2), DCD, DCU, DU2, IPIV, INFO )
  
  
time:  DO istep = 1, Nstep
! Build the rhs  
spc:   DO ispec = 1,Nspec
          ! Stage 1
          call advdiff_flux_fvz(N,Z,W,K,Bdry(1,ispec), &
	         SurfaceEm(ispec),Vd(ispec),C(1,ispec),S1)
	  S1(1:N) = S1(1:N) + VolumeEm(1:N,ispec)	 
          call TRI_SOL(N, DCD, DCU, DCL, S1)          
	  !call DGTTRS( 'N', N, 1, DCL(2), DCD, DCU, DU2, &
	  !           IPIV, S1, N, INFO )
          ! Stage 2
	  C1(1:N) = C(1:N,ispec) + DT*S1(1:N)	     
          call advdiff_flux_fvz(N,Z,W,K,Bdry(1,ispec), &
	         SurfaceEm(ispec),Vd(ispec),C1,S2)
	  S2 = S2 + VolumeEm(1:N,ispec) - 2*S1	 
          call TRI_SOL(N, DCD, DCU, DCL, S2)          
	  !call DGTTRS( 'N', N, 1, DCL(2), DCD, DCU, DU2, &
	  !           IPIV, S2, N, INFO )
          ! The Solution 
          C(1:N,ispec) = C(1:N,ispec) + DT*1.5*S1(1:N) &
	                              + DT*0.5*S2(1:N)
     END DO spc
  END DO time

END SUBROUTINE ADVDIFF_FVZ


! Advection-diffusion derivative by finite volumes
! Advection is discretized by simple upwind
subroutine ADVDIFF_FLUX_FVZ(N,Z,W,K,Bdry,SurfEm,Vd,C,DC)
!
implicit none
!
integer, intent(in) :: N
double precision, intent(in)  :: Z(N), W(N), K(N), &
                          Vd, SurfEm, Bdry(2),C(N)
! Time derivative of the concentration
real*8, intent(out) :: DC(N)

! difflux/advflux = diffusive/advective fluxes through i-1/2
integer :: i
real*8 :: difflux, advflux
! flux = total flux through i-1/2
real*8 :: flux(N+1)

! Leftmost boundary
   ! B(1) = - K dc/dz at z(1)
   !     = (-U(1)*(Bdry(1)-C(1)) - Vd*C(1) + E)*xi(N,X,U,K,X(1),1)
if ( W(1)>0.d0 ) then
     advflux = 0 ! W(1)*Bdry(1)
     difflux = Vd*C(1) - SurfEm 
else
     advflux = 0 ! W(1)*C(1)
     difflux = Vd*C(1) - SurfEm 
end if  
flux(1) =  difflux - advflux 

! Intermediate Boundaries
do i=2,N
  if (W(i)>=0) then
    advflux = W(i-1)*C(i-1)
  else
    advflux = W(i)*C(i)
  end if  
  difflux = (K(i)+K(i-1))/2.d0*(C(i)-C(i-1))/(Z(i)-Z(i-1))
  flux(i) = difflux - advflux
end do

! Top of the domain
if ( W(N)<0 ) then ! inflow
    advflux = W(N)*Bdry(2)
    difflux = 0.d0
else ! outflow
    advflux = W(N)*C(N)
    difflux = 0.d0
end if      
flux(N+1) =  difflux - advflux 

! Time derivatives
DC(1) = (flux(2)-flux(1))/(z(2)-z(1))
do i=2,N
   DC(i) = (flux(i+1)-flux(i))/(z(i+1)-z(i-1))*2.d0
end do
DC(N) = (flux(N+1)-flux(N))/(z(N)-z(N-1))

end subroutine ADVDIFF_FLUX_FVZ


! Advection-diffusion derivative by finite volumes
! Advection is discretized by simple upwind
subroutine ADVDIFF_JAC_FVZ(N,Z,W,K,Vd,C,DCL,DCD,DCU)
!
implicit none
!
integer, intent(in) :: N
double precision, intent(in)  :: Z(N), W(N), K(N), Vd, C(N)
! Time derivative of the concentration
real*8, intent(out) :: DCL(N), DCD(N), DCU(N)

! difflux/advflux = diffusive/advective fluxes through i-1/2
integer :: i
real*8 :: difflux, advflux
! flux = total flux through i-1/2
real*8 :: dflux(N+1), lflux(N+1)

! Leftmost boundary
if ( W(1)>0.d0 ) then
    ! flux(1) = Vd*C(1) - SurfEm - W(1)*Bdry(1)
    dflux(1) = Vd
    lflux(1) = 0.d0
else
    ! flux(1) = Vd*C(1) - SurfEm - W(1)*C(1)
    dflux(1) = Vd-W(1)
    lflux(1) = 0.d0
end if  

! Intermediate Boundaries
do i=2,N
  if (W(i)>=0) then
    ! flux(i) = (K(i)+K(i-1))/2.d0*(C(i)-C(i-1))/(Z(i)-Z(i-1)) - W(i-1)*C(i-1)
    dflux(i) =(K(i)+K(i-1))/2.d0/(Z(i)-Z(i-1)) 
    lflux(i) = -(K(i)+K(i-1))/2.d0/(Z(i)-Z(i-1)) - W(i-1)
  else
    ! flux(i) = (K(i)+K(i-1))/2.d0*(C(i)-C(i-1))/(Z(i)-Z(i-1)) - W(i)*C(i)
    dflux(i) =(K(i)+K(i-1))/2.d0/(Z(i)-Z(i-1))  - W(i)
    lflux(i) = -(K(i)+K(i-1))/2.d0/(Z(i)-Z(i-1)) 
  end if  
end do

! Top of the domain
if ( W(N)<0 ) then ! inflow
    ! flux(N+1) = - W(N)*Bdry(2)
    dflux(N+1) = 0.d0
    lflux(N+1) = 0.d0
else ! outflow
    ! flux(N+1) = - W(N)*C(N)
    dflux(N+1) = 0.d0
    lflux(N+1) = -W(N)
end if      

! Time derivatives
! DC(1) = (flux(2)-flux(1))/(z(2)-z(1))
DCL(1) = 0.d0
DCD(1) = (lflux(2)-dflux(1))/(z(2)-z(1))
DCU(1) = (dflux(2))/(z(2)-z(1))
do i=2,N
   ! DC(i) = (flux(i+1)-flux(i))/(z(i+1)-z(i-1))*2.d0
   DCL(i) = (-lflux(i))/(z(i+1)-z(i-1))*2.d0
   DCD(i) = (lflux(i+1)-dflux(i))/(z(i+1)-z(i-1))*2.d0
   DCU(i) = (dflux(i+1))/(z(i+1)-z(i-1))*2.d0
end do
! DC(N) = (flux(N+1)-flux(N))/(z(N)-z(N-1))
DCL(N) = (-lflux(N))/(z(N)-z(N-1))
DCD(N) = (lflux(N+1)-dflux(N))/(z(N)-z(N-1))
DCU(N) = 0.d0

end subroutine ADVDIFF_JAC_FVZ



      SUBROUTINE TRI_FA(N, A, B, C)
! FACTORIZES TRIDIAGONAL SYSTEM 
! A=DIAGONAL, B=SUPRA, C=SUB, INDEX=LINE NO.
      INTEGER N 
      REAL*8 A(N), B(N), C(N)
!      
      DO I=1, N-1
	C(I)   = - C(I+1)/A(I)
        A(I+1) = A(I+1) + C(I)*B(I)
      END DO      
!
      RETURN
      END

      SUBROUTINE TRI_SOL(N, A, B, C, Y)
! SOLVES DIFFUSION (TRIDIAGONAL SYSTEM) 
! A=DIAGONAL, B=SUPRA, C=SUB, INDEX=LINE NO. 
      REAL*8 Y(N), A(N), B(N), C(N)
!       
      DO I=1, N-1
        Y(I+1) = Y(I+1) + C(I)*Y(I)
      ENDDO      
!
      Y(N) = Y(N)/A(N)
!
      DO I=N-1,1,-1
        Y(I) = ( Y(I) - B(I)*Y(I+1) )/A(I)
      END DO
!
      RETURN
      END      

