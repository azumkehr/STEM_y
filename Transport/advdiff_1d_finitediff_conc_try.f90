!------------------------------------------------------------------------------
! Contains the solution for advection-diffusion eqns by finite differences
! ***** CONCENTRATIONS ARE ABSOLUTE, E.G. MOLEC/CM3 *****
!------------------------------------------------------------------------------


!*************************************************************************
      subroutine tranx_fd(ix,iy,iz,bounds,s1,u,kh,sx,dt,dx)
!**************************************************************************
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7--
      include 'aqmax.param'
      include 'aqindx.cmm'
      integer :: ix, iy, iz
      integer :: bounds(8)
      real :: s1(ix,iy,iz,*), u(ix,iy,*), kh(ix,iy,*)
      real :: sx(iy,iz,2,*), dx(*), dt
      !real :: bndf(2,mxspg),bndv(2,mxspg)
      !real :: gen(mxgr*mxspg)
      real*8 :: wind(1:ix), dif(1:ix), air(1:ix), conc(1:ix,mxspg)
      real*8 :: sb(1:2,mxspg), Deltax, Deltat
      integer :: num, i, j, k, Nsteps
      integer :: xstart,xend,ystart,yend
      integer :: zstart,zend,nstart,nend
!
!    The computational bounds 
      xstart = bounds(1); xend = bounds(2)
      ystart = bounds(3); yend = bounds(4)
      zstart = bounds(5); zend = bounds(6)
      nstart = bounds(7); nend = bounds(8)

!
!                           transport in x-direction

      num = nend-nstart+1         ! number of transported species
            
      Deltax = dx(1)       ! grid size: this works for uniform grids only
      ! Max time step is 900 sec
      Nsteps = int(dt/900.0)+1
      Deltat = dt/Dble(Nsteps)        
      
      do j=ystart,yend
      do k=zstart,zend
      
      
      Conc(1:ix,1:num) = s1(1:ix,j,k,nstart:nend)
      Wind(1:ix) = u(1:ix,j,k)
      Dif(1:ix) = kh(1:ix,j,k)
      ! Air(1:ix) = s1(1:ix,j,k,iair)
      Sb(1:2,1:num) = sx(j,k,1:2,1:num)
! For air
      Conc(1:ix,num+1) = s1(1:ix,j,k,iair)
      Sb(1:2,num+1)    = sx(j,k,1:2,iair)      
      
      call ADVDIFF_FDH_U(Deltat, Nsteps, ix, num+1, Deltax, &
                         Wind, Dif, Conc, Sb)
			 
      s1(1:ix,j,k,nstart:nend) = max(Conc(1:ix,1:num),0.d0)
      s1(1:ix,j,k,iair) =  max(Conc(1:ix,num+1),0.d0)
      
      end do ! iz
      end do ! iy      
      
      return
      end

!*********************************************************************
      subroutine trany_fd(ix,iy,iz,bounds,s1,v,kh,sy,dt,dy)
!*********************************************************************
      include 'aqmax.param'
      include 'aqindx.cmm'
      integer :: ix, iy, iz
      integer :: bounds(8)
      real ::  s1(ix,iy,iz,*), v(ix,iy,*), kh(ix,iy,*)
      real ::  sy(ix,iz,2,*), dy(*), dt
      !real ::  bndf(2,mxspg),bndv(2,mxspg)
      !real ::  gen(mxgr*mxspg)
      real*8 :: Wind(1:iy), Dif(1:iy), Air(1:iy), Conc(1:iy,mxspg)
      real*8 :: Sb(1:2,mxspg), Deltay, Deltat
      integer :: num, i, j, k, Nsteps
      integer :: xstart,xend,ystart,yend
      integer :: zstart,zend,nstart,nend
!
!    The computational bounds 
      xstart = bounds(1); xend = bounds(2)
      ystart = bounds(3); yend = bounds(4)
      zstart  = bounds(5); zend = bounds(6)
      nstart = bounds(7); nend = bounds(8)
 
      num=nend-nstart+1    ! number of transported species
      Deltay = dy(1)       ! grid size: this works for uniform grids only
      ! Max time step is 300 sec
      Nsteps = int(dt/300.0)+1
      Deltat = dt/Dble(Nsteps)        

      do i=xstart,xend
      do k=zstart,zend
      
      Conc(1:iy,1:num) = s1(i,1:iy,k,nstart:nend)
      Wind(1:iy) = v(i,1:iy,k)
      Dif(1:iy) = kh(i,1:iy,k)
      ! Air(1:iy) = s1(i,1:iy,k,iair)
      Sb(1:2,1:num) = sy(i,k,1:2,1:num)
! For air
      Conc(1:iy,num+1) = s1(i,1:iy,k,iair)
      Sb(1:2,num+1) = sy(i,k,1:2,iair)
      
      call ADVDIFF_FDH_U(Deltat, Nsteps, iy, num+1, Deltay, &
                         Wind, Dif, Conc, Sb)
			 
      s1(i,1:iy,k,nstart:nend) = max(Conc(1:iy,1:num),0.d0)
      s1(i,1:iy,k,iair) = max(Conc(1:iy,num+1),0.d0)
      
      end do ! iz
      end do ! ix

      return
      end
      
      
!*********************************************************************
      subroutine tranz_fd(ix,iy,iz,bounds,s1,w,kv,q,em,vg,sz,dt,dz)
!*********************************************************************
      include 'aqmax.param'
      include 'aqindx.cmm'
!      Arguments:      
      integer :: ix, iy, iz
      integer :: bounds(8)
      real :: s1(ix,iy,iz,*), w(ix,iy,iz), kv(ix,iy,iz)
      real :: q(ix,iy,*), em(ix,iy,iz,*), vg(ix,iy,*), sz(ix,iy,*)
      real :: dt, dz(ix,iy,*)
! Local variables:      
      real*8 :: Deltat, Z(iz), Wind(iz), Dif(iz)
      real*8 :: Conc(iz,mxspg), Sb(2,mxspg), Air(iz) 
      real*8 :: SurfaceEm(mxspg), VolumeEm(iz,mxspg), DepVel(mxspg)     
      real*8 :: Conc1(iz,mxspg)
      integer :: i, j, Nsteps
      integer :: xstart,xend,ystart,yend
      integer :: zstart,zend,nstart,nend
!
!    The computational bounds 
      xstart = bounds(1); xend = bounds(2)
      ystart = bounds(3); yend = bounds(4)
      zstart = bounds(5); zend = bounds(6)
      nstart = bounds(7); nend = bounds(8)

      num = nend-nstart+1 ! number of species
      ! Max time step is 300 sec
      Nsteps = int(dt/300.0)+1
      Deltat = dt/Dble(Nsteps)        
      
      do i=xstart,xend
      do j=ystart,yend
      
      ! The vertical grid
      Z(1) = 0.d0
      do k = 1, iz-1
        Z(k+1) = Z(k) + dz(i,j,k)
      end do	
      
      
      Wind(1:iz)  = w(i,j,1:iz)
      Dif(1:iz)   = kv(i,j,1:iz)
      Conc(1:iz,1:num) = s1(i,j,1:iz,nstart:nend)
      ! Air(1:iz) = s1(i,j,1:iz,iair)
      Sb(1,1:num) = s1(i,j,1,nstart:nend)
      Sb(2,1:num) = sz(i,j,nstart:nend)
      SurfaceEm(1:num) = q(i,j,nstart:nend)
      VolumeEm(1:iz,1:num) = em(i,j,1:iz,nstart:nend)
      DepVel(1:num)        = vg(i,j,nstart:nend)
! For Air:
      Conc(1:iz,num+1) = s1(i,j,1:iz,iair)
      Sb(1,num+1) = s1(i,j,1,iair)
      Sb(2,num+1) = sz(i,j,iair)
      SurfaceEm(num+1) = 0.d0
      VolumeEm(1:iz,num+1) = 0.d0
      DepVel(num+1)        = 0.d0
              
      call ADVDIFF_FDZ(Deltat, Nsteps, iz, num+1, Z, Wind, Dif, &
                    Conc,Sb, SurfaceEm, VolumeEm, DepVel)   
			 
      s1(i,j,1:iz,nstart:nend) = max(Conc(1:iz,1:num),0.d0)
      s1(i,j,1:iz,iair) = max(Conc(1:iz,num+1),0.d0)
		       
      end do ! j
      end do ! i
!
      return
      end
  
      


!*************************************************************************
! X-Transport in the sigma vertical coordinate
      subroutine tranx_sigma_fd(ix,iy,iz,bounds,s1,u,kh,&
                            sx,dt,dx,deltah)
!**************************************************************************
!---+----1----+----2----+----3----+----4----+----5----+----6----+----7--
      include 'aqmax.param'
      include 'aqindx.cmm'
      integer :: ix, iy, iz
      integer :: bounds(8)
      real :: s1(ix,iy,iz,*), u(ix,iy,*), kh(ix,iy,*)
      real :: sx(iy,iz,2,*), dx(*), dt, deltah(ix,iy)
      !real :: bndf(2,mxspg),bndv(2,mxspg)
      !real :: gen(mxgr*mxspg)
      real*8 :: wind(1:ix), dif(1:ix), air(1:ix), conc(1:ix,mxspg)
      real*8 :: sb(1:2,mxspg), Deltax, Deltat, dh(1:ix)
      integer :: num, i, j, k, Nsteps
      integer :: xstart,xend,ystart,yend
      integer :: zstart,zend,nstart,nend
!
!    The computational bounds 
      xstart = bounds(1); xend = bounds(2)
      ystart = bounds(3); yend = bounds(4)
      zstart = bounds(5); zend = bounds(6)
      nstart = bounds(7); nend = bounds(8)

      num = nend-nstart+1         ! number of transported species
            
      Deltax = dx(1)       ! grid size: this works for uniform grids only
      ! Max time step is 300 sec
      Nsteps = int(dt/300.0)+1
      Deltat = dt/Dble(Nsteps)        
      
      do j=ystart,yend
      do k=zstart,zend
      
      ! A scaled version of terrain height
      dh(1:ix) = deltah(1:ix,j) !/maxval( deltah(1:ix,j) )
      
      do m=1,num
        Conc(1:ix,m) = s1(1:ix,j,k,nstart+m-1)*dh(1:ix)
      end do
      Wind(1:ix) = u(1:ix,j,k)
      Dif(1:ix)  = kh(1:ix,j,k)
      ! Air(1:ix) = s1(1:ix,j,k,iair)
      Sb(1,1:num) = sx(j,k,1,1:num)*dh(1)
      Sb(2,1:num) = sx(j,k,2,1:num)*dh(ix)
! For air
      Conc(1:ix,num+1) = s1(1:ix,j,k,iair)*dh(1:ix)
      Sb(1,num+1)      = sx(j,k,1,iair)*dh(1)      
      Sb(2,num+1)      = sx(j,k,2,iair)*dh(ix)      
     
      call ADVDIFF_FDH_U(Deltat, Nsteps, ix, num+1, Deltax, &
                         Wind, Dif, Conc, Sb)
			
      do m=1,num 
        s1(1:ix,j,k,nstart+m-1) = max(Conc(1:ix,m)/dh(1:ix),0.d0)
      end do
      s1(1:ix,j,k,iair) =  max(Conc(1:ix,num+1)/dh(1:ix),0.d0)
      
      end do ! iz
      end do ! iy      
      
      return
      end

!*********************************************************************
! Y-Transport in the sigma vertical coordinate
      subroutine trany_sigma_fd(ix,iy,iz,bounds,s1,v,kh, &
                             sy,dt,dy,deltah)
!*********************************************************************
      include 'aqmax.param'
      include 'aqindx.cmm'
      integer :: ix, iy, iz
      integer :: bounds(8)
      real ::  s1(ix,iy,iz,*), v(ix,iy,*), kh(ix,iy,*)
      real ::  sy(ix,iz,2,*), dy(*), dt, deltah(ix,iy)
      real*8 :: Wind(1:iy), Dif(1:iy), Air(1:iy), Conc(1:iy,mxspg)
      real*8 :: Sb(1:2,mxspg), Deltay, Deltat, dh(1:iy)
      integer :: num, i, j, k, Nsteps
      integer :: xstart,xend,ystart,yend
      integer :: zstart,zend,nstart,nend
!
!    The computational bounds 
      xstart = bounds(1); xend = bounds(2)
      ystart = bounds(3); yend = bounds(4)
      zstart = bounds(5); zend = bounds(6)
      nstart = bounds(7); nend = bounds(8)

!      
      num=nend-nstart+1    ! number of transported species
      Deltay = dy(1)       ! grid size: this works for uniform grids only
      ! Max time step is 300 sec
      Nsteps = int(dt/300.0)+1
      Deltat = dt/Dble(Nsteps)        

      do i=xstart,xend
      do k=zstart,zend
      
      ! A scaled version of terrain height
      dh(1:iy) = deltah(i,1:iy) ! /maxval( deltah(i,1:iy) )
     
      do m=1,num 
        Conc(1:iy,m) = s1(i,1:iy,k,nstart+m-1)*dh(1:iy)
      end do
      Wind(1:iy) = v(i,1:iy,k)
      Dif(1:iy) = kh(i,1:iy,k)
      ! Air(1:iy) = s1(i,1:iy,k,iair)
      Sb(1,1:num) = sy(i,k,1,1:num)*dh(1)
      Sb(2,1:num) = sy(i,k,2,1:num)*dh(iy)
! For air
      Conc(1:iy,num+1) = s1(i,1:iy,k,iair)*dh(1:iy)
      Sb(1,num+1)      = sy(i,k,1,iair)*dh(1)
      Sb(2,num+1)      = sy(i,k,2,iair)*dh(iy)
      
      call ADVDIFF_FDH_U(Deltat, Nsteps, iy, num+1, Deltay, &
                         Wind, Dif, Conc, Sb)
	
      do m=1,num		 
          s1(i,1:iy,k,nstart+m-1) = max(Conc(1:iy,m)/dh(1:iy),0.d0)
      end do
      s1(i,1:iy,k,iair) = max(Conc(1:iy,num+1)/dh(1:iy),0.d0)
      
      end do ! iz
      end do ! ix

      return
      end

      
!*********************************************************************
! Z-Transport in the sigma vertical coordinate
      subroutine tranz_sigma_fd(ix,iy,iz,bounds,s1,w,kv,&
                             q,em,vg,sz,dt,sigma)
!*********************************************************************
      include 'aqmax.param'
      include 'aqindx.cmm'
!      Arguments:      
      integer :: ix, iy, iz
      integer :: bounds(8)
      real :: s1(ix,iy,iz,*), w(ix,iy,iz), kv(ix,iy,iz)
      real :: q(ix,iy,*), em(ix,iy,iz,*), vg(ix,iy,*), sz(ix,iy,*)
      real :: dt, sigma(iz)
! Local variables:      
      real*8 :: Deltat, Z(iz), Wind(iz), Dif(iz)
      real*8 :: Conc(iz,mxspg), Sb(2,mxspg), Air(iz) 
      real*8 :: SurfaceEm(mxspg), VolumeEm(iz,mxspg), DepVel(mxspg)     
      real*8 :: Conc1(iz,mxspg)
      integer :: i, j, Nsteps
      integer :: xstart,xend,ystart,yend
      integer :: zstart,zend,nstart,nend
!
!    The computational bounds 
      xstart = bounds(1); xend = bounds(2)
      ystart = bounds(3); yend = bounds(4)
      zstart = bounds(5); zend = bounds(6)
      nstart = bounds(7); nend = bounds(8)

      num = nend-nstart+1 ! number of species
      ! Max time step is 300 sec
      Nsteps = int(dt/300.0)+1
      Deltat = dt/Dble(Nsteps)        
      
      do i=xstart,xend
      do j=ystart,yend
      
      ! The vertical grid
      Z(1:iz) = sigma(1:iz)
      
      Wind(1:iz)  = w(i,j,1:iz)
      Dif(1:iz)   = kv(i,j,1:iz)
      Conc(1:iz,1:num) = s1(i,j,1:iz,nstart:nend)
      ! Air(1:iz) = s1(i,j,1:iz,iair)
      Sb(1,1:num) = s1(i,j,1,nstart:nend)
      Sb(2,1:num) = sz(i,j,nstart:nend)
      SurfaceEm(1:num) = q(i,j,nstart:nend)
      VolumeEm(1:iz,1:num) = em(i,j,1:iz,nstart:nend)
      DepVel(1:num)        = vg(i,j,nstart:nend)
! For Air:
      Conc(1:iz,num+1) = s1(i,j,1:iz,iair)
      Sb(1,num+1) = s1(i,j,1,iair)
      Sb(2,num+1) = sz(i,j,iair)
      SurfaceEm(num+1) = 0.d0
      VolumeEm(1:iz,num+1) = 0.d0
      DepVel(num+1)        = 0.d0
             
      call ADVDIFF_FDZ(Deltat, Nsteps, iz, num+1, Z, Wind, Dif, &
                    Conc,Sb, SurfaceEm, VolumeEm, DepVel)   

      s1(i,j,1:iz,nstart:nend) = max(Conc(1:iz,1:num),0.d0)
      s1(i,j,1:iz,iair) = max(Conc(1:iz,num+1),0.d0)
		       
      end do ! j
      end do ! i
!
      return
      end
  
      

subroutine ADVDIFF_FDH_U(DT, Nstep, N, Nspec, DX, U, K, C, BDRY)
! -----------------------------------------------------------------------------
!  Performs Nstep timesteps of length DT
!      to solve the adv_diff equation using linearly implicit midpoint
! 
!  N   = no. of grid points
!  Nspec = no. of chemical species
!  Nstep = no of time steps
!  X(1:N) = grid point coordinates
!  U(1:N) = wind speeds
!  K(1:N) = diffusion coefficients
!  SurfaceEm  = Surface Emission intensity
!  VolumeEm   = Elevated Emission intensity
!  Vd    = deposition velocity
!  C     = concentration of each species
! -----------------------------------------------------------------------------

  implicit none
  integer, intent(in) :: N, Nstep, Nspec
  double precision, intent(in)  :: DT, DX, U(N), K(N), &
			  Bdry(2,Nspec)
  double precision, intent(inout) :: C(N,Nspec)

!  Local Variables
  integer, parameter :: kl=2, ku=2
  integer, parameter :: ldjac=kl+ku+1, lda = 2*kl+ku+1
  double precision, parameter :: alpha = 1.d0, beta = 1.d0
  double precision :: Jac(ldjac,N),  A(lda,N)
  double precision :: C1(N), B(N), D(2), S1(N) 
  integer :: istep, ispec, info, ipiv(N)

!  The Jacobian
  call ADVDIFF_JAC_FDH(N,DX,U,K,Jac)  
!  A = I - DT*gam*Jac
  A(kl+1:lda,1:N) = -DT/2.d0*Jac(1:ldjac,1:N)
  A(kl+ku+1,1:N)  = 1.d0 + A(kl+ku+1,1:N)  ! add 1 to diagonal terms
  call DGBTRF( N, N, kl, ku, A, lda, IPIV, INFO )
  if (INFO.ne.0) then
     print*,'In ADVDIFF_FDH_U INFO = ',INFO
  end if
            
time:do istep = 1, Nstep
spc:   do ispec = 1, Nspec
          C1(1:N) = C(1:N,ispec)
	  D(1:2)  = Bdry(1:2,ispec)
          ! The free term
	  call ADVDIFF_FREE_FDH(N,DX,U,K,D,B)
	  C1(1:N) = C1(1:N) + DT*B(1:N)
	  call DGBMV('N', N, N, kl, ku, DT/2.d0, Jac, ldjac,&
                    C(1,ispec), 1, beta, C1, 1)
          !OR: call Banded_times_vector(N,kl,ku,Jac,C1,S1); S1 = S1 + B
          call DGBTRS( 'N', N, kl, ku, 1, A, lda, IPIV, &
	              C1, N, INFO )
          C(1:N,ispec) = C1(1:N)
       end do spc
     end do   time
   
end subroutine ADVDIFF_FDH_U

     

     

subroutine ADJ_ADVDIFF_FDH(DT, Nstep, N, Nspec, DX, U, K, Lam)
! -----------------------------------------------------------------------------
!  The adjoint of the above function
! 
!  N   = no. of grid points
!  Nspec = no. of chemical species
!  Nstep = no of time steps
!  X(1:N) = grid point coordinates
!  U(1:N) = wind speeds
!  K(1:N) = diffusion coefficients
!  SurfaceEm  = Surface Emission intensity
!  VolumeEm   = Elevated Emission intensity
!  Vd      = deposition velocity
!  Lam     = adjoint of concentration of each species
! -----------------------------------------------------------------------------

  implicit none
  integer, intent(in) :: N, Nstep, Nspec
  double precision, intent(in)  :: DT, DX, U(N), K(N)
  double precision, intent(inout) :: Lam(N,Nspec)

!  Local Variables
  integer, parameter :: kl=2, ku=2
  integer, parameter :: ldjac=kl+ku+1, lda = 2*kl+ku+1
  double precision, parameter :: alpha = 1.d0, beta = 1.d0
  double precision :: Jac(ldjac,N),  A(lda,N)
  double precision :: C1(N), B(N), D(2), S1(N) 
  integer :: istep, ispec, info, ipiv(N)

!  The Jacobian
  call ADVDIFF_JAC_FDH(N,DX,U,K,Jac)  
!  A = I - DT*gam*Jac
  A(kl+1:lda,1:N) = -DT/2.d0*Jac(1:ldjac,1:N)
  A(kl+ku+1,1:N)  = 1.d0 + A(kl+ku+1,1:N)  ! add 1 to diagonal terms
  call DGBTRF( N, N, kl, ku, A, lda, IPIV, INFO )
  if (INFO.ne.0) then
     print*,'In ADJ_ADVDIFF_FDH_MF INFO = ',INFO
  end if
            
time:do istep = 1, Nstep
spc:   do ispec = 1, Nspec
          call DGBTRS( 'T', N, kl, ku, 1, A, lda, IPIV, &
	              Lam(1,ispec), N, INFO )
	  call DGBMV('T', N, N, kl, ku, DT/2.d0, Jac, ldjac,&
                    Lam(1,ispec), 1, beta, Lam(1,ispec), 1)
       end do spc
     end do   time
   
end subroutine ADJ_ADVDIFF_FDH




subroutine ADVDIFF_FDH_ROS2(DT, Nstep, N, Nspec, DX, U, K, C, BDRY)
! -----------------------------------------------------------------------------
!  Performs Nstep timesteps of length DT
!      to solve the adv_diff equation using Ros2
! 
!  N   = no. of grid points
!  Nspec = no. of chemical species
!  Nstep = no of time steps
!  X(1:N) = grid point coordinates
!  U(1:N) = wind speeds
!  K(1:N) = diffusion coefficients
!  SurfaceEm  = Surface Emission intensity
!  VolumeEm   = Elevated Emission intensity
!  Vd    = deposition velocity
!  C     = concentration of each species
! -----------------------------------------------------------------------------

  implicit none
  integer, intent(in) :: N, Nstep, Nspec
  double precision, intent(in)  :: DT, DX, U(N), K(N), &
			  Bdry(2,Nspec)
  double precision, intent(inout) :: C(N,Nspec)

!  Local Variables
  integer, parameter :: kl=2, ku=2
  integer, parameter :: ldjac=kl+ku+1, lda = 2*kl+ku+1
  double precision, parameter :: alpha = 1.d0, beta = 1.d0
  double precision :: Jac(ldjac,N),  A(lda,N)
  double precision :: gam, C1(N), B(N), S1(N), S2(N)
  integer :: istep, ispec, info, ipiv(N)

  gam = 1 + sqrt(2.d0)/2.d0
!  The Jacobian
  call ADVDIFF_JAC_FDH(N,DX,U,K,Jac)  
!  A = I - DT*gam*Jac
  A(kl+1:lda,1:N) = -DT*gam*Jac(1:ldjac,1:N)
  A(kl+ku+1,1:N)  = 1.d0 + A(kl+ku+1,1:N)  ! add 1 to diagonal terms
  call DGBTRF( N, N, kl, ku, A, lda, IPIV, INFO )
  if (INFO.ne.0) then
     print*,'In ADVDIFF_FDH_ROS2 INFO = ',INFO
  end if
            
time:do istep = 1, Nstep
spc:   do ispec = 1, Nspec
	  call ADVDIFF_FUN_FDH(N,DX,U,K,Bdry(1,ispec),C(1,ispec),S1)
	  ! Stage 1: 
          call DGBTRS( 'N', N, kl, ku, 1, A, lda, IPIV, &
	              S1, N, INFO )
          ! Stage 2: S2 = A\(Jac*C1 + B - 2*S1)
	  C1(1:N) = C(1:N,ispec) + DT*S1(1:N)	     
	  call ADVDIFF_FUN_FDH(N,DX,U,K,Bdry(1,ispec),C1,S2)	
	  S2 = S2 - 2*S1
	  call DGBTRS( 'N', N, kl, ku, 1, A, lda, IPIV, &
	              S2, N, INFO )
          ! Next-step solution 
          C(1:N,ispec) = C(1:N,ispec) + 1.5d0*DT*S1(1:N) &
	                              + 0.5d0*DT*S2(1:N) 
       end do spc
     end do   time
   
end subroutine ADVDIFF_FDH_ROS2

! -----------------------------------------------------------------------------
! Advection-diffusion derivative by finite differences
! Advection is discretized by third order, unlimited upwind
! Dirichlet b.c. and uniform grid
! Horizontal transport in STEM
! -----------------------------------------------------------------------------
subroutine ADVDIFF_FUN_FDH(N,DX,U,K,Bdry,C,DC)
!
implicit none
!
integer, intent(in) :: N
double precision, intent(in)  :: DX, U(N), K(N), Bdry(2), C(N)
! Time derivative of the concentration
real*8, intent(out) :: DC(N)

! difflux/advflux = diffusive/advective fluxes through i-1/2
integer :: i
real*8 :: difflux, advflux
! Concentration and boundaries in a single vector
real*8 :: F(0:N+1)

!
real*8, parameter :: ap = -1.0d0/6.0d0, bp = 1.0d0, &
               cp = -1.0d0/2.0d0, dp = -1.0d0/3.0d0, &
	       an = 1.d0/3.d0, bn = 1.d0/2.d0, &
	       cn = -1.d0, dn = 1.d0/6.d0

F(1:N) = C(1:N)*U(1:N)
F(0)   = Bdry(1)*U(1)
F(N+1) = Bdry(2)*U(N)


! The advection discretization
if ( U(1) >= 0.d0 ) then   !  inflow
   DC(1) = ( F(0)-F(1) )/DX
else  !  outflow
   DC(1) = ( F(1)-F(2) )/DX
end if 
!
do i=2,N-1
  if ( U(i) >= 0.d0 ) then  !  inflow
    DC(i) = ( ap*F(i-2)+bp*F(i-1)+cp*F(i)+dp*F(i+1) )/DX
  else
    DC(i) = ( an*F(i-1)+bn*F(i)+cn*F(i+1)+dn*F(i+2) )/DX
  end if 
end do
!
if ( U(N) >= 0.d0 ) then  ! outflow
  DC(N) = ( F(N-1)-F(N) )/DX
else   ! inflow
  DC(N) = ( F(N)-F(N+1) )/DX
end if 

! The diffusion part
if ( U(1) >= 0.d0 ) then   !  inflow
  DC(1) = DC(1) + ( (K(1)+K(2))*(C(2)-C(1))  &
                -   2*K(1)*(C(1)-Bdry(1)) &
	          )/(2*DX**2)
else  !  outflow
  DC(1) = DC(1) + ( (K(1)+K(2))*(C(2)-C(1))  &
	         )/(2*DX**2)
end if 
!
do i=2,N-1
  DC(i) = DC(i) + ( (K(i+1)+K(i))*(C(i+1)-C(i))  &
                -   (K(i)+K(i-1))*(C(i)-C(i-1))  &
	           )/(2*DX**2)
end do
!
if ( U(N) >= 0.d0 ) then  ! outflow
  DC(N) = DC(N) + (                      &
                -  (K(N)+ K(N-1))*(C(N)-C(N-1)) &
	           )/(2*DX**2)
else   ! inflow
  DC(N) = DC(N) + ( 2*K(N)*(Bdry(2)-C(N))  &
                -   (K(N)+K(N-1))*(C(N)-C(N-1)) &
	           )/(2*DX**2)
end if 
!
end subroutine ADVDIFF_FUN_FDH


! -----------------------------------------------------------------------------
! Jacobian of Advection-diffusion derivative by finite differences
! Advection is discretized by third order, unlimited upwind
! Dirichlet b.c. and uniform grid
! The Jacobian is pentadiagonal in Blas banded representation:
!      Jac( KU + 1 - J + I, J ) = matrix( I, J ), J = 1, N,  I = MAX( 1, J - KU ), MIN( M, J + KL )
! -----------------------------------------------------------------------------

subroutine ADVDIFF_JAC_FDH(N,DX,U,K,Jac)
!
implicit none
!
integer, intent(in) :: N
double precision, intent(in)  :: DX, U(N), K(N)
! Time derivative of the concentration
integer, parameter :: ku=2, kl=2
real*8, intent(out) :: Jac(ku+kl+1,N)

integer :: i

!
real*8, parameter :: ap = -1.0d0/6.0d0, bp = 1.0d0, &
               cp = -1.0d0/2.0d0, dp = -1.0d0/3.0d0, &
	       an = 1.d0/3.d0, bn = 1.d0/2.d0, &
	       cn = -1.d0, dn = 1.d0/6.d0

!  Initialize Jacobian to zzero
Jac(1:ku+kl+1,1:N) = 0.d0

! The advection discretization
if ( U(1) >= 0.d0 ) then
   ! DC(1) = 1.d0/DX*( F(0)-F(1) )
   Jac(jrow(1,1),1) = -U(1)/DX
else
   !   DC(1) = 1.d0/DX*( F(1)-F(2) ) 
    Jac(jrow(1,1),1) =  U(1)/DX
    Jac(jrow(1,2),2) = -U(2)/DX
end if 
!
do i=2,N-1
  if ( U(i) >= 0.d0 ) then
    ! DC(i) = 1.d0/DX*( ap*F(i-2)+bp*F(i-1)+cp*F(i)+dp*F(i+1) )
    if (i>2) Jac(jrow(i,i-2),i-2) = ap*U(i-2)/DX
    Jac(jrow(i,i-1),i-1) = bp*U(i-1)/DX
    Jac(jrow(i,i  ),i  ) = cp*U(i  )/DX
    Jac(jrow(i,i+1),i+1) = dp*U(i+1)/DX
  else
    ! DC(i) = 1.d0/DX*( an*F(i-1)+bn*F(i)+cn*F(i+1)+dn*F(i+2) )
    Jac(jrow(i,i-1),i-1) = an*U(i-1)/DX
    Jac(jrow(i,i  ),i  ) = bn*U(i  )/DX
    Jac(jrow(i,i+1),i+1) = cn*U(i+1)/DX
    if (i<N-1) Jac(jrow(i,i+2),i+2) = dn*U(i+2)/DX
  end if 
end do
!
if ( U(N) >= 0.d0 ) then
  ! DC(N) = 1.d0/DX*( F(N-1)-F(N) )
  Jac(jrow(N,N-1),N-1) =  U(N-1)/DX
  Jac(jrow(N,N)  ,N  ) = -U(N)  /DX
else
  ! DC(N) =1.d0/DX*( F(N)-F(N+1) )
  Jac(jrow(N,N)  ,N  ) =  U(N)/DX
end if 

! The diffusion part
if ( U(1) >= 0.d0 ) then   !  inflow
  !DC(1) = DC(1) + ( (K(1)+K(2))*(C(2)-C(1))  &
  !              -   2*K(1)*(C(1)-Bdry(1))  )/(2*DX**2)
  Jac(jrow(1,1),1) = Jac(jrow(1,1),1) &
                     - (K(2)+3*K(1))/(2*DX**2)
  Jac(jrow(1,2),2) = Jac(jrow(1,2),2) &
                     +  (K(1)+K(2))/(2*DX**2)		  
else  !  outflow
  !DC(1) = DC(1) + ( (K(1)+K(2))*(C(2)-C(1)) )/(2*DX**2)
  Jac(jrow(1,1),1) = Jac(jrow(1,1),1) - (K(1)+K(2))/(2*DX**2)
  Jac(jrow(1,2),2) = Jac(jrow(1,2),2) + (K(1)+K(2))/(2*DX**2)		  
end if 
!
do i=2,N-1
  !DC(i) = DC(i) + ( K(i+1)*(C(i+1)-C(i))  &
  !              -   K(i-1)*(C(i)-C(i-1))  &
  !	           )/(2*DX**2)
  Jac(jrow(i,i-1),i-1) = Jac(jrow(i,i-1),i-1)   &
             + (K(i)+K(i-1))/(2*DX**2)
  Jac(jrow(i,i)  ,i)   = Jac(jrow(i,i) ,i)      &
             - (K(i+1)+2*K(i)+K(i-1))/(2*DX**2)
  Jac(jrow(i,i+1),i+1) = Jac(jrow(i,i+1),i+1)   &
             +  (K(i+1)+K(i))/(2*DX**2)
end do
!
if ( U(N) >= 0.d0 ) then  ! outflow
  !DC(N) = DC(N) + (                      &
  !              -  (K(N)+ K(N-1))*(C(N)-C(N-1)) &
!	           )/(2*DX**2)
  Jac(jrow(N,N-1),N-1) = Jac(jrow(N,N-1),N-1) &
                          + (K(N)+K(N-1))/(2*DX**2)
  Jac(jrow(N,N)  ,N)   = Jac(jrow(N,N),N)     &
                          - (K(N)+K(N-1))/(2*DX**2)
else   ! inflow
  !DC(N) = DC(N) + ( 2*K(N)*(Bdry(2)-C(N))  &
  !              -   (K(N)+K(N-1))*(C(N)-C(N-1)) &
!	           )/(2*DX**2)
  Jac(jrow(N,N-1),N-1) = Jac(jrow(N,N-1),N-1)  &
                          + (K(N)+K(N-1))/(2*DX**2)
  Jac(jrow(N,N)  ,N)   = Jac(jrow(N,N),N)      &
                          - (3*K(N)+K(N-1))/(2*DX**2)
end if 

contains
  
  integer function jrow(i,j)
  ! gives the row of the Blas banded format for pentadiagonal Jacobian
  integer :: i, j
  integer, parameter :: kl=2, ku=2
  if ( (i<=0) .or. (j<=0) ) then
     print*,'Error in ADVDIFF_JAC_FDH. i,j=',i,j
     stop
  end if
  jrow = ku + 1 + i - j
  end function jrow

end subroutine ADVDIFF_JAC_FDH



! -----------------------------------------------------------------------------
! Advection-diffusion derivative by finite differences
! Advection is discretized by third order, unlimited upwind
! Dirichlet b.c. and uniform grid
! Free term B such that: c' = Fun_fdh = Jac_fdh*c + B
! -----------------------------------------------------------------------------
subroutine ADVDIFF_FREE_FDH(N,DX,U,K,Bdry,B)
!
implicit none
!
integer, intent(in) :: N
double precision, intent(in)  :: DX, U(N), K(N), Bdry(2)
! Time derivative of the concentration
real*8, intent(out) :: B(N)

! difflux/advflux = diffusive/advective fluxes through i-1/2
integer :: i
real*8 :: difflux, advflux
! Concentration and boundaries in a single vector
real*8 :: F(0:N+1)

!
real*8, parameter :: ap = -1.0d0/6.0d0, bp = 1.0d0, &
               cp = -1.0d0/2.0d0, dp = -1.0d0/3.0d0, &
	       an = 1.d0/3.d0, bn = 1.d0/2.d0, &
	       cn = -1.d0, dn = 1.d0/6.d0

B(1:N) = 0.d0

! The advection discretization
if ( U(1) >= 0.d0 ) then   
   B(1) = Bdry(1)*U(1)/DX + 2*K(1)*Bdry(1)/(2*DX**2)
end if 
!
if ( U(2) >= 0.d0 ) then  
   B(2) = ap*Bdry(1)*U(1)/DX
end if 
!
if ( U(N-1) < 0.d0 ) then  
   B(N-1) = dn*Bdry(2)*U(N)/DX
end if 
!
if ( U(N) < 0.d0 ) then  
  B(N) = -Bdry(2)*U(N)/DX + 2*K(N)*Bdry(2)/(2*DX**2)
end if 
!
end subroutine ADVDIFF_FREE_FDH



subroutine ADVDIFF_FDZ(DT, Nstep, N, Nspec, X, U, K, C, &
                    BDRY, SurfaceEm, VolumeEm, Vd)
! -----------------------------------------------------------------------------
!  Performs Nstep timesteps of length DT
!      to solve the adv_diff equation in vertical direction
!      using finite volume method and Ros2
! 
!  N     = no. of grid points
!  Nspec = no. of chemical species
!  Nstep = no of time steps
!  X(1:NGP) = grid point coordinates
!  U(1:NGP) = wind speeds
!  K(1:NGP) = diffusion coefficients
!  SurfaceEm  = Surface Emission intensity
!  VolumeEm   = Elevated Emission intensity
!  Vd    = deposition velocity
!  C     = concentration of each species
!
! Note: it uses Midpoint rule
! -----------------------------------------------------------------------------
!
  implicit none
  integer, intent(in) :: N, Nstep, Nspec
  double precision, intent(in)  :: DT, X(N), U(N), K(N), &
                          Vd(Nspec), SurfaceEm(Nspec), &
			  Bdry(2,Nspec), VolumeEm(N, Nspec)
  double precision, intent(inout) :: C(N,Nspec)

!  Local Variables
  integer, parameter :: kl=1, ku=1
  integer, parameter :: ldjac=kl+ku+1, lda = 2*kl+ku+1
  double precision, parameter :: alpha = 1.d0, beta = 1.d0
  double precision :: Jac(ldjac,N),  A(lda,N)
  double precision :: C1(N), B(N), D(2)
  integer :: istep, ispec, info, ipiv(N)

!  The Jacobian
  call ADVDIFF_JAC_FDZ(N,X,U,K,Jac)  
!  A = I - DT*gam*Jac
  A(kl+1:lda,1:N) = -DT/2.d0*Jac(1:ldjac,1:N)
  A(kl+ku+1,1:N)  = 1.d0 + A(kl+ku+1,1:N)  ! add 1 to diagonal terms
  call DGBTRF( N, N, kl, ku, A, lda, IPIV, INFO )
  if (INFO.ne.0) then
     print*,'In ADVDIFF_FDZ INFO = ',INFO
  end if
            
time:do istep = 1, Nstep
spc:   do ispec = 1, Nspec
          C1(1:N) = C(1:N,ispec)
	  D(1:2)  = Bdry(1:2,ispec)
          ! The free term
	  call ADVDIFF_FREE_FDZ(N,X,U,K,D,SurfaceEm(ispec), &
	           Vd(ispec),C1,B)
	  ! Stage 1: S1 = A\(Jac*C1 + B)
	  C1(1:N) = C1(1:N) + DT*B(1:N) + DT*VolumeEm(1:N,ispec)
	  call DGBMV('N', N, N, kl, ku, DT/2.d0, Jac, ldjac,&
                    C(1:N,ispec), 1, beta, C1, 1)
          !OR: call Banded_times_vector(N,kl,ku,Jac,C1,S1); S1 = S1 + B
          call DGBTRS( 'N', N, kl, ku, 1, A, lda, IPIV, &
	              C1, N, INFO )
          C(1:N,ispec) = C1(1:N)
       end do spc
     end do   time
   
end subroutine ADVDIFF_FDZ




subroutine ADJ_ADVDIFF_FDZ(DT, Nstep, N, Nspec, X, U, K, Lam, &
                       BDRY, SurfaceEm, VolumeEm, Vd)
! -----------------------------------------------------------------------------
!  The adjoint of the above
! 
!  N     = no. of grid points
!  Nspec = no. of chemical species
!  Nstep = no of time steps
!  X(1:NGP) = grid point coordinates
!  U(1:NGP) = wind speeds
!  K(1:NGP) = diffusion coefficients
!  SurfaceEm  = Surface Emission intensity
!  VolumeEm   = Elevated Emission intensity
!  Vd    = deposition velocity
!  Lam     = adjoint of concentration of each species
!
! Note: it uses Midpoint rule
! -----------------------------------------------------------------------------
!
  implicit none
  integer, intent(in) :: N, Nstep, Nspec
  double precision, intent(in)  :: DT, X(N), U(N), K(N), &
                          Vd(Nspec), SurfaceEm(Nspec), &
			  Bdry(2,Nspec), VolumeEm(N, Nspec)
  double precision, intent(inout) :: Lam(N,Nspec)

!  Local Variables
  integer, parameter :: kl=1, ku=1
  integer, parameter :: ldjac=kl+ku+1, lda = 2*kl+ku+1
  double precision, parameter :: alpha = 1.d0, beta = 1.d0
  double precision :: Jac(ldjac,N),  A(lda,N)
  double precision :: C1(N), B(N), D(2)
  integer :: istep, ispec, info, ipiv(N)

!  The Jacobian
  call ADVDIFF_JAC_FDZ(N,X,U,K,Jac)  
!  A = I - DT*gam*Jac
  A(kl+1:lda,1:N) = -DT/2.d0*Jac(1:ldjac,1:N)
  A(kl+ku+1,1:N)  = 1.d0 + A(kl+ku+1,1:N)  ! add 1 to diagonal terms
  call DGBTRF( N, N, kl, ku, A, lda, IPIV, INFO )
  if (INFO.ne.0) then
     print*,'In ADJ_ADVDIFF_FDZ INFO = ',INFO
  end if
            
time:do istep = 1, Nstep
spc:   do ispec = 1, Nspec
          call DGBTRS( 'T', N, kl, ku, 1, A, lda, IPIV, &
	              Lam(1,ispec), N, INFO )
          C1(1:N) = Lam(1:N,ispec) 
	  Lam(1,ispec) = (1.d0 - DT*Vd(ispec)/(X(2)-X(1)))*C1(1)
	  call DGBMV('T', N, N, kl, ku, DT/2.d0, Jac, ldjac,&
                    C1, 1, beta, Lam(1,ispec), 1)
       end do spc
     end do   time
   
end subroutine ADJ_ADVDIFF_FDZ




! Advection-diffusion derivative by finite volumes
! Advection is discretized by simple upwind
subroutine ADVDIFF_FUN_FDZ(N,Z,W,K,Bdry,SurfEm,Vd,C,DC)
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
     advflux = W(1)*Bdry(1)
     difflux = Vd*C(1) - SurfEm 
else
     advflux = W(1)*C(1)
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
    advflux =  W(N)*Bdry(2)
    difflux = -W(N)*(Bdry(2)-C(N))
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

end subroutine ADVDIFF_FUN_FDZ




! -----------------------------------------------------------------------------
! Advection-diffusion derivative by finite volumes
! Advection is discretized by simple upwind
! Jac is in Blas-banded-format:  Jac(1:3,j) = A(j-1:j+1,j)
! -----------------------------------------------------------------------------
subroutine ADVDIFF_JAC_FDZ(N,Z,W,K,Jac)
!
implicit none
!
integer, intent(in) :: N
double precision, intent(in)  :: Z(N), W(N), K(N)
! Jacobian for time derivative of the concentration
integer, parameter :: kl=1, ku=1
real*8, intent(out) :: Jac(kl+ku+1,N)

! difflux/advflux = diffusive/advective fluxes through i-1/2
integer :: i
real*8 :: difflux, advflux
! flux = total flux through i-1/2
real*8 :: dflux(N+1), lflux(N+1)

! Leftmost boundary
! Note: The -Vd*C(1) term is NOT incorporated into the Jacobian
!       consequently the Jacobian is species-independent
if ( W(1)>0.d0 ) then
    ! flux(1) = Vd*C(1) - SurfEm - W(1)*Bdry(1)
    dflux(1) = 0.d0 ! Vd
    lflux(1) = 0.d0
else
    ! flux(1) = Vd*C(1) - SurfEm - W(1)*C(1)
    dflux(1) = -W(1) ! Vd-W(1)
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
    dflux(N+1) =  0.d0
    lflux(N+1) =  W(N)
else ! outflow
    ! flux(N+1) = - W(N)*C(N)
    dflux(N+1) =  0.d0
    lflux(N+1) = -W(N)
end if      

Jac = 0.d0
! DC(1) = (flux(2)-flux(1))/(z(2)-z(1))
Jac(jrow(1,1),1) = (lflux(2)-dflux(1))/(z(2)-z(1))
Jac(jrow(1,2),2) = (dflux(2))/(z(2)-z(1))
do i=2,N-1
   ! DC(i) = (flux(i+1)-flux(i))/(z(i+1)-z(i-1))*2.d0
   Jac(jrow(i,i-1),i-1) = (-lflux(i))/(z(i+1)-z(i-1))*2.d0
   Jac(jrow(i,i),i)     = (lflux(i+1)-dflux(i))/(z(i+1)-z(i-1))*2.d0
   Jac(jrow(i,i+1),i+1) = (dflux(i+1))/(z(i+1)-z(i-1))*2.d0
end do
! DC(N) = (flux(N+1)-flux(N))/(z(N)-z(N-1))
Jac(jrow(N,N-1),N-1) = (-lflux(N))/(z(N)-z(N-1))
Jac(jrow(N,N),N)     = (lflux(N+1)-dflux(N))/(z(N)-z(N-1))

contains
  
  integer function jrow(i,j)
  ! gives the row of the Blas banded format for pentadiagonal Jacobian
  integer :: i, j
  integer, parameter :: kl=1, ku=1
  if ( (i<=0) .or. (j<=0) ) then
     print*,'Error in ADVDIFF_JAC_FDZ. i,j=',i,j
     stop
  end if
  jrow = ku + 1 + i - j
  end function jrow

end subroutine ADVDIFF_JAC_FDZ



!-------------------------------------------------------------------------------
! Advection-diffusion derivative by finite volumes
! Advection is discretized by simple upwind
! The free vertical term B such that c' = fun_fdz = jac_fdz * c + B
!-------------------------------------------------------------------------------
subroutine ADVDIFF_FREE_FDZ(N,Z,W,K,Bdry,SurfEm,Vd,C,B)
!
implicit none
!
integer, intent(in) :: N
double precision, intent(in)  :: Z(N), W(N), K(N), &
                          Vd, SurfEm, Bdry(2),C(N)
! The free term
real*8, intent(out) :: B(N)

! difflux/advflux = diffusive/advective fluxes through i-1/2
integer :: i
real*8 :: difflux, advflux

B(1:N) = 0.d0
! Leftmost boundary
   ! B(1) = - K dc/dz at z(1)
   !     = (-U(1)*(Bdry(1)-C(1)) - Vd*C(1) + E)*xi(N,X,U,K,X(1),1)
if ( W(1)>0.d0 ) then
     advflux = W(1)*Bdry(1)
     difflux = Vd*C(1) - SurfEm 
else
     advflux = 0.d0
     difflux = Vd*C(1) - SurfEm 
end if  
B(1) =  - (difflux - advflux)/(z(2)-z(1)) 

! Top of the domain
if ( W(N)<0 ) then ! inflow
    advflux = W(N)*Bdry(2)
    difflux = K(N)*Bdry(2)/(Z(N)-Z(N-1)) !-W(N)*Bdry(2)
else ! outflow
    advflux = 0.d0
    difflux = 0.d0
end if      
B(N) = (difflux - advflux)/(z(N)-z(N-1))


end subroutine ADVDIFF_FREE_FDZ

! -----------------------------------------------------------------------------
! Jacobian of Advection-diffusion derivative by finite differences
! Advection is discretized by third order, unlimited upwind
! The grid is general grid
! The Jacobian is pentadiagonal in Blas banded representation:
!      Jac( KU + 1 - J + I, J ) = matrix( I, J ), J = 1, N,  I = MAX( 1, J - KU ), MIN( M, J + KL )
! -----------------------------------------------------------------------------

subroutine ADVDIFF_JAC_FD3(N,X,U,K,Jac)
!
implicit none
!
integer, intent(in) :: N
double precision, intent(in)  :: X(N), U(N), K(N)
! Time derivative of the concentration
integer, parameter :: ku=2, kl=2
real*8, intent(out) :: Jac(ku+kl+1,N)
!
integer :: i
!
real*8 :: a, b, c, d

!  Initialize Jacobian to zzero
Jac(1:ku+kl+1,1:N) = 0.d0

! The advection discretization
if ( U(1) >= 0.d0 ) then
   ! DC(1) = 1.d0/DX*( F(0)-F(1) )
   Jac(jrow(1,1),1) = -U(1)/(X(2)-X(1))
else
   !   DC(1) = 1.d0/DX*( F(1)-F(2) ) 
    Jac(jrow(1,1),1) =  U(1)/(X(2)-X(1))
    Jac(jrow(1,2),2) = -U(2)/(X(2)-X(1))
end if 
!
i=2
if ( U(2) >= 0.d0 ) then
    call coeff_fd3(+1,X(1),X(i-1),X(i),X(i+1),a,b,c,d)
    Jac(jrow(i,i-1),i-1) = b*U(i-1)
    Jac(jrow(i,i  ),i  ) = c*U(i  )
    Jac(jrow(i,i+1),i+1) = d*U(i+1)
else
    call coeff_fd3(-1,X(i-1),X(i),X(i+1),X(i+2),a,b,c,d)
    Jac(jrow(i,i-1),i-1) = a*U(i-1)
    Jac(jrow(i,i  ),i  ) = b*U(i  )
    Jac(jrow(i,i+1),i+1) = c*U(i+1)
    if (i<N-1) Jac(jrow(i,i+2),i+2) = d*U(i+2)
end if 
!
do i=3,N-2
  if ( U(i) >= 0.d0 ) then
    call coeff_fd3(+1,X(i-2),X(i-1),X(i),X(i+1),a,b,c,d)
    Jac(jrow(i,i-2),i-2) = a*U(i-2)
    Jac(jrow(i,i-1),i-1) = b*U(i-1)
    Jac(jrow(i,i  ),i  ) = c*U(i  )
    Jac(jrow(i,i+1),i+1) = d*U(i+1)
  else
    call coeff_fd3(-1,X(i-1),X(i),X(i+1),X(i+2),a,b,c,d)
    Jac(jrow(i,i-1),i-1) = a*U(i-1)
    Jac(jrow(i,i  ),i  ) = b*U(i  )
    Jac(jrow(i,i+1),i+1) = c*U(i+1)
    Jac(jrow(i,i+2),i+2) = d*U(i+2)
  end if 
end do
!
i=N-1
if ( U(i) >= 0.d0 ) then
    call coeff_fd3(+1,X(i-2),X(i-1),X(i),X(i+1),a,b,c,d)
    Jac(jrow(i,i-2),i-2) = a*U(i-2)
    Jac(jrow(i,i-1),i-1) = b*U(i-1)
    Jac(jrow(i,i  ),i  ) = c*U(i  )
    Jac(jrow(i,i+1),i+1) = d*U(i+1)
else
    call coeff_fd3(-1,X(N-2),X(N-1),X(N),X(N),a,b,c,d)
    Jac(jrow(i,i-1),i-1) = a*U(i-1)
    Jac(jrow(i,i  ),i  ) = b*U(i  )
    Jac(jrow(i,i+1),i+1) = c*U(i+1)
end if 
!
if ( U(N) >= 0.d0 ) then
  ! DC(N) = 1.d0/DX*( F(N-1)-F(N) )
  Jac(jrow(N,N-1),N-1) =  U(N-1)/(X(N)-X(N-1))
  Jac(jrow(N,N)  ,N  ) = -U(N)  /(X(N)-X(N-1))
else
  ! DC(N) =1.d0/DX*( F(N)-F(N+1) )
  Jac(jrow(N,N)  ,N  ) =  U(N)/(X(N)-X(N-1))
end if 


contains
  
  integer function jrow(i,j)
  ! gives the row of the Blas banded format for pentadiagonal Jacobian
  integer :: i, j
  integer, parameter :: kl=2, ku=2
  if ( (i<=0) .or. (j<=0) ) then
     print*,'Error in ADVDIFF_JAC_FDH. i,j=',i,j
     stop
  end if
  jrow = ku + 1 + i - j
  end function jrow

end subroutine ADVDIFF_JAC_FD3

!--------------------------------------------------------------------
! Finite difference coefficients
! isign = +1: positive wind 
!             a f(x1) + b f(x2) + c f(x3) + d f(x4 ) ~= - df/dx (x3)
! isign = -1: negative wind 
!             a f(x1) + b f(x2) + c f(x3) + d f(x4 ) ~= - df/dx (x2)
!--------------------------------------------------------------------
subroutine coeff_fd3(isign,x1,x2,x3,x4,a,b,c,d)
implicit none
integer, intent(in) :: isign
real*8, intent(in) :: x1,x2,x3,x4
real*8, intent(out) :: a,b,c,d
real*8 :: h1, h2, h3, z

h1 = x2-x1; h2 = x3-x2; h3 = x4-x3
z = h1*(h1+h2)*(h1+h2+h3)

if ( isign == 1 ) then
    a = -h2*h3/z
    b = (h1+h2)*h3/(h1*h2*(h2+h3))
    d = -h2*(h1+h2)/z
    c = -a-b-d
else if (isign == -1) then
    a = h2*(h2+h3)/z
    c = -h1*(h2+h3)/(h2*h3*(h1+h2))
    d = h1*h2/z
    b = -a-c-d
end if
 
end subroutine coeff_fd3

!--------------------------------------------------------------------
! Blas-format banded matrix times vector:  y = A*x
! Note: A(KU+1+i-j) = A_normal(i,j)
!--------------------------------------------------------------------
subroutine Banded_times_vector(N,KL,KU,A,X,Y)
implicit none
integer, intent(in) :: N, KL, KU
real*8, intent(in) :: A(KL+KU+1,N), X(N)
real*8, intent(out) :: Y(N)
!
integer :: i, j

Y(1:N) = 0.d0
do j=1,N
  do i = max(j-KU,1), min(j+KL,N)
   Y(i) = Y(i) + X(j)*A(KU+1+i-j,j)
  end do 
end do

end subroutine Banded_times_vector


subroutine ADVDIFF_JAC_OLD(N,DX,U,K,C,Jac)
!
implicit none
!
integer, intent(in) :: N
double precision, intent(in)  :: DX, U(N), K(N), C(N)
! Time derivative of the concentration
real*8, intent(out) :: Jac(7,N)

integer :: i

!
real*8, parameter :: ap = -1.0d0/6.0d0, bp = 1.0d0, &
               cp = -1.0d0/2.0d0, dp = -1.0d0/3.0d0, &
	       an = 1.d0/3.d0, bn = 1.d0/2.d0, &
	       cn = -1.d0, dn = 1.d0/6.d0

!  Initialize Jacobian to zzero
Jac(1:5,1:N) = 0.d0

! The advection discretization
if ( U(1) >= 0.d0 ) then
   ! DC(1) = 1.d0/DX*( F(0)-F(1) )
   Jac(jrow(1,1),1) = -U(1)/DX
else
   !   DC(1) = 1.d0/DX*( F(1)-F(2) ) 
    Jac(jrow(1,1),1) =  U(1)/DX
    Jac(jrow(1,2),2) = -U(2)/DX
end if 
!
do i=2,N-1
  if ( U(i) >= 0.d0 ) then
    ! DC(i) = 1.d0/DX*( ap*F(i-2)+bp*F(i-1)+cp*F(i)+dp*F(i+1) )
    if (i>2) Jac(jrow(i,i-2),i-2) = ap*U(i-2)/DX
    Jac(jrow(i,i-1),i-1) = bp*U(i-1)/DX
    Jac(jrow(i,i  ),i  ) = cp*U(i  )/DX
    Jac(jrow(i,i+1),i+1) = dp*U(i+1)/DX
  else
    ! DC(i) = 1.d0/DX*( an*F(i-1)+bn*F(i)+cn*F(i+1)+dn*F(i+2) )
    Jac(jrow(i,i-1),i-1) = an*U(i-1)/DX
    Jac(jrow(i,i  ),i  ) = bn*U(i  )/DX
    Jac(jrow(i,i+1),i+1) = cn*U(i+1)/DX
    if (i<N-1) Jac(jrow(i,i+2),i+2) = dn*U(i+2)/DX
  end if 
end do
!
if ( U(N) >= 0.d0 ) then
  ! DC(N) = 1.d0/DX*( F(N-1)-F(N) )
  Jac(jrow(N,N-1),N-1) =  U(N-1)/DX
  Jac(jrow(N,N)  ,N  ) = -U(N)  /DX
else
  ! DC(N) =1.d0/DX*( F(N)-F(N+1) )
  Jac(jrow(N,N)  ,N  ) =  U(N)/DX
end if 

! The diffusion part
if ( U(1) >= 0.d0 ) then   !  inflow
  !DC(1) = DC(1) + ( K(2)*(C(2)-C(1))  &
  !              -   K(1)*(C(1)-Bdry(1)) &
  !	          )/(2*DX**2)
  Jac(jrow(1,1),1) = Jac(jrow(1,1),1) - (K(2)+K(1))/(2*DX**2)
  Jac(jrow(1,2),2) = Jac(jrow(1,2),2) +  K(2)/(2*DX**2)		  
else  !  outflow
  ! DC(1) = DC(1) + ( K(2)*(C(2)-C(1))  &
  ! 	         )/(2*DX**2)
  Jac(jrow(1,1),1) = Jac(jrow(1,1),1) - K(2)/(2*DX**2)
  Jac(jrow(1,2),2) = Jac(jrow(1,2),2) + K(2)/(2*DX**2)		  
end if 
!
do i=2,N-1
  !DC(i) = DC(i) + ( K(i+1)*(C(i+1)-C(i))  &
  !              -   K(i-1)*(C(i)-C(i-1))  &
  !	           )/(2*DX**2)
  Jac(jrow(i,i-1),i-1) = Jac(jrow(i,i-1),i-1) + K(i-1)/(2*DX**2)
  Jac(jrow(i,i)  ,i)   = Jac(jrow(i,i) ,i) - (K(i+1)+K(i-1))/(2*DX**2)
  Jac(jrow(i,i+1),i+1) = Jac(jrow(i,i+1),i+1) +  K(i+1)/(2*DX**2)
end do
if ( U(N) >= 0.d0 ) then  ! outflow
  !DC(N) = DC(N) + (                      &
  !              -   K(N-1)*(C(N)-C(N-1)) &
  !	           )/(2*DX**2)
  Jac(jrow(N,N-1),N-1) = Jac(jrow(N,N-1),N-1) + K(N)/(2*DX**2)
  Jac(jrow(N,N)  ,N)   = Jac(jrow(N,N),N) - (K(N)+K(N-1))/(2*DX**2)
else   ! inflow
  !DC(N) = DC(N) + ( K(N)*(Bdry(2)-C(N))  &
  !              -   K(N-1)*(C(N)-C(N-1)) &
  !	           )/(2*DX**2)
  Jac(jrow(N,N-1),N-1) = Jac(jrow(N,N-1),N-1) + K(N-1)/(2*DX**2)
  Jac(jrow(N,N)  ,N)   = Jac(jrow(N,N),N) - (K(N)+K(N-1))/(2*DX**2)
end if 

contains
  
  integer function jrow(i,j)
  ! gives the row of the compressed format pentadiagonal Jacobian
  integer :: i, j
  integer, parameter :: kl=2, ku=2
  if ( (i<=0) .or. (j<=0) ) then
     print*,'Error in ADVDIFF_JAC_FDH. i,j=',i,j
     stop
  end if
  jrow = ku + 1 + i - j
  end function jrow

end subroutine ADVDIFF_JAC_OLD


