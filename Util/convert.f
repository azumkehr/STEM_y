! Converts the first L species

      subroutine convert_4D(Nx,Ny,Nz,Ns,iair,num,s,how)
      implicit none
      integer :: Nx,Ny,Nz,Ns,iair,num,L
      real :: s(Nx,Ny,Nz,Ns)
      character(len=9) :: how
      
      if( how == 'conc2frac' ) then
        do L=1,num
	  s(1:Nx,1:Ny,1:Nz,L) = s(1:Nx,1:Ny,1:Nz,L)
     &                           /s(1:Nx,1:Ny,1:Nz,iair)	  
        end do
      else if (how == 'frac2conc' ) then
        do L=1,num
	  s(1:Nx,1:Ny,1:Nz,L) = s(1:Nx,1:Ny,1:Nz,L)
     &                           *s(1:Nx,1:Ny,1:Nz,iair)	  
        end do      
      else
        print*,'Error: Unknown conversion type ',how
	stop
      end if      
      
      end subroutine convert_4D
