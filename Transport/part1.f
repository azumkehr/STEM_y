subroutine ADJ_ADVDIFF_FDZ_MF(DT, Nstep, N, Nspec, X, U, K, Air, &
                       Lam, BDRY, SurfaceEm, VolumeEm, Vd)
! -----------------------------------------------------------------------------
!  Performs Nstep timesteps of length DT
!      to solve the adv_diff equation in vertical direction
!      using finite volume method and Ros2
! 
!  N     = no. of grid points
!  Nspec = no. of chemical species
!  Nstep = no of time steps
!  X(1:N) = grid point coordinates
!  U(1:N) = wind speeds
!  K(1:N) = diffusion coefficients
! Air(1:N) = Air density
!  SurfaceEm  = Surface Emission intensity
!  VolumeEm(1:N)   = Elevated Emission intensity
!  Vd    = deposition velocity
!  Lam(1:N)     = concentration of each species
!
! Note: it uses Midpoint rule
! -----------------------------------------------------------------------------
!
  implicit none
  integer, intent(in) :: N, Nstep, Nspec
  double precision, intent(in)  :: DT, X(N), U(N), K(N), Air(N), &
                          Vd(Nspec), SurfaceEm(Nspec),   &
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
  call ADVDIFF_JAC_FDZ_MF(N,X,U,K,Air,Jac)  
!  A = I - DT*gam*Jac
  A(kl+1:lda,1:N) = -DT/2.d0*Jac(1:ldjac,1:N)
  A(kl+ku+1,1:N)  = 1.d0 + A(kl+ku+1,1:N)  ! add 1 to diagonal terms
  call DGBTRF( N, N, kl, ku, A, lda, IPIV, INFO )
  if (INFO.ne.0) then
     print*,'In ADVDIFF_FDZ_MF INFO = ',INFO
  end if
            
time:do istep = 1, Nstep
spc:   do ispec = 1, Nspec
          call DGBTRS('T', N, kl, ku, 1, A, lda, IPIV, &
          	       Lam(1,ispec), N, INFO )
          C1(1:N) = Lam(1:N,ispec) 
	  Lam(1,ispec) = (1.d0 - DT*Vd(ispec)/(X(2)-X(1)))*C1(1)
	  call DGBMV('T', N, N, kl, ku, DT/2.d0, Jac, ldjac,&
                       C1, 1, beta, Lam(1,ispec), 1)
       end do spc
     end do   time
   
end subroutine ADJ_ADVDIFF_FDZ_MF
