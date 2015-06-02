!
! Changes U, V, W, Kh and Kv from Cartesian to sigma coordinates
!
! Note: the modified equation will use (DeltaH*C) as the variable.
!
!
!
!
subroutine cartesian2sigma(Nx,Ny,Nz,dx,dy,sigma,Ground,DeltaH, &
                           U,V,W,Kh,Kv)
!
! The domain is Nx * Ny * Nz
! Ground(Nx,Ny) = Ground level altitude
! Top(Nx,Ny)   = Top of the domain altitude
!

!
integer, intent(in) :: Nx, Ny, Nz
real, intent(in)   :: dx, dy, sigma(Nz)
real, intent(in) :: Ground(Nx, Ny),DeltaH(Nx, Ny)
real, intent(in) :: Kh(Nx,Ny,Nz)
!
real, intent(inout) :: U(Nx,Ny,Nz), V(Nx,Ny,Nz), W(Nx,Ny,Nz)
real, intent(inout) :: Kv(Nx,Ny,Nz)
!
real :: Ground_x(Nx,Ny), Ground_y(Nx,Ny)
real :: DeltaH_x(Nx,Ny), DeltaH_y(Nx,Ny)

! x and y Derivatives
!
do j = 1, Ny
  Ground_x(1,j)   = (Ground(2,j)-Ground(1,j))/(dx)
  Ground_x(Nx,j)  = (Ground(Nx,j)-Ground(Nx-1,j))/(dx)
  DeltaH_x(1,j)   = (DeltaH(2,j)-DeltaH(1,j))/(dx)
  DeltaH_x(Nx,j)  = (DeltaH(Nx,j)-DeltaH(Nx-1,j))/(dx)
  do i = 2, Nx-1
    Ground_x(i,j) = (Ground(i+1,j)-Ground(i-1,j))/(2.d0*dx)
    DeltaH_x(i,j) = (DeltaH(i+1,j)-DeltaH(i-1,j))/(2.d0*dx)
  end do
end do
!
do i = 1, Nx
  Ground_y(i,1)   = (Ground(i,2)-Ground(i,1))/(dy)
  Ground_y(i,Ny)  = (Ground(i,Ny)-Ground(i,Ny-1))/(dy)
  DeltaH_y(i,1)   = (DeltaH(i,2)-DeltaH(i,1))/(dy)
  DeltaH_y(i,Ny)  = (DeltaH(i,Ny)-DeltaH(i,Ny-1))/(dy)
  do j = 2, Ny-1
    Ground_y(i,j) = (Ground(i,j+1)-Ground(i,j-1))/(2.d0*dy)
    DeltaH_y(i,j) = (DeltaH(i,j+1)-DeltaH(i,j-1))/(2.d0*dy)
  end do
end do

! The new Quantities
do i = 1, Nx
  do j=1, Ny
     do k=1, Nz
	 W(i,j,k) = ( W(i,j,k) &
	            -(Ground_x(i,j)+sigma(k)*DeltaH_x(i,j))*U(i,j,k) &
		    -(Ground_y(i,j)+sigma(k)*DeltaH_y(i,j))*V(i,j,k) &
		     )/DeltaH(i,j)
     end do
  end do   
end do

! The new Quantities
do i = 1, Nx
  do j=1, Ny
     do k=1, Nz
         U(i,j,k) =  U(i,j,k) + Kh(i,j,k)*DeltaH_x(i,j)/DeltaH(i,j)
     end do
  end do
end do   

! The new Quantities
do i = 1, Nx
  do j=1, Ny
     do k=1, Nz
         V(i,j,k) =  V(i,j,k) + Kh(i,j,k)*DeltaH_y(i,j)/DeltaH(i,j)
     end do
  end do   
end do

! The new Quantities
do i = 1, Nx
  do j=1, Ny
     do k=1, Nz
	 Kv(i,j,k) = Kv(i,j,k)/DeltaH(i,j)**2	     
     end do
  end do   
end do

end subroutine cartesian2sigma

! Coordinate transformation for mole fraction
subroutine cartesian2sigma_mf(Nx,Ny,Nz,dx,dy,sigma,Ground,DeltaH, &
                           U,V,W,Kh,Kv)
!
! The domain is Nx * Ny * Nz
! Ground(Nx,Ny) = Ground level altitude
! Top(Nx,Ny)   = Top of the domain altitude
!

!
integer, intent(in) :: Nx, Ny, Nz
real, intent(in)   :: dx, dy, sigma(Nz)
real, intent(in) :: Ground(Nx, Ny),DeltaH(Nx, Ny)
real, intent(in) :: Kh(Nx,Ny,Nz)
!
real, intent(inout) :: U(Nx,Ny,Nz), V(Nx,Ny,Nz), W(Nx,Ny,Nz)
real, intent(inout) :: Kv(Nx,Ny,Nz)
!
real :: Ground_x(Nx,Ny), Ground_y(Nx,Ny)
real :: DeltaH_x(Nx,Ny), DeltaH_y(Nx,Ny)

! x and y Derivatives
!
do j = 1, Ny
  Ground_x(1,j)   = (Ground(2,j)-Ground(1,j))/(dx)
  Ground_x(Nx,j)  = (Ground(Nx,j)-Ground(Nx-1,j))/(dx)
  DeltaH_x(1,j)   = (DeltaH(2,j)-DeltaH(1,j))/(dx)
  DeltaH_x(Nx,j)  = (DeltaH(Nx,j)-DeltaH(Nx-1,j))/(dx)
  do i = 2, Nx-1
    Ground_x(i,j) = (Ground(i+1,j)-Ground(i-1,j))/(2.d0*dx)
    DeltaH_x(i,j) = (DeltaH(i+1,j)-DeltaH(i-1,j))/(2.d0*dx)
  end do
end do
!
do i = 1, Nx
  Ground_y(i,1)   = (Ground(i,2)-Ground(i,1))/(dy)
  Ground_y(i,Ny)  = (Ground(i,Ny)-Ground(i,Ny-1))/(dy)
  DeltaH_y(i,1)   = (DeltaH(i,2)-DeltaH(i,1))/(dy)
  DeltaH_y(i,Ny)  = (DeltaH(i,Ny)-DeltaH(i,Ny-1))/(dy)
  do j = 2, Ny-1
    Ground_y(i,j) = (Ground(i,j+1)-Ground(i,j-1))/(2.d0*dy)
    DeltaH_y(i,j) = (DeltaH(i,j+1)-DeltaH(i,j-1))/(2.d0*dy)
  end do
end do

! The new Quantities
do i = 1, Nx
  do j=1, Ny
     do k=1, Nz
	 W(i,j,k) = ( W(i,j,k) &
	            -(Ground_x(i,j)+sigma(k)*DeltaH_x(i,j))*U(i,j,k) &
		    -(Ground_y(i,j)+sigma(k)*DeltaH_y(i,j))*V(i,j,k) &
		     )/DeltaH(i,j)
     end do
  end do   
end do

! The new Quantities
do i = 1, Nx
  do j=1, Ny
     do k=1, Nz
	 Kv(i,j,k) = Kv(i,j,k)/DeltaH(i,j)**2	     
     end do
  end do   
end do


end subroutine cartesian2sigma_mf
