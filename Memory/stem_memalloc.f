      module StemMemAlloc
      
      contains

c********************************************************************
      subroutine MemAlloc( N_x, N_y, N_z, 
     &                   N_gas, N_liquid, N_particle,
     &                   sg1, sl1, sp1,
     &                   u, v, w,
     &                   kh, kv, t, dz,
     &                   wc, wr, sprc, rvel,
     &                   sx, sy, sz, q, em,
     &                   vg, fz, hdz,
     &                   h, deltah, tlon, tlat, cldod, 
     &                   kctop, ccover, dobson)
c-----------------------------------------------------------------------c
c
c                  A. Sandu, Dec. 2000  
c
c-----------------------------------------------------------------------c
c Allocates memory for the following arrays:
c
c Nx, Ny, Nz = dimension of the grid   
c N_gas      = no. of gas species
c N_liquid   = no. of liquid species
c N_particle = no of particulate species  
c sg1        = gas concentrations field
c sl1        = liquid concentrations field
c sp1        = particulate matter concentration field
c u, v, w    = wind field components
c kh, kv     = diffusivity coefficients
c t          = global time
c dz         = vertical resolution
c wc, wr     = photolysis rate data
c sprc       = unknown
c rvel       = removal velocity
c sx, sy, sz = boundary concentrations
c q          = surface emission rates 
c em         = elevated emission rates 
c vg         = deposition velocity 
c fz         = unknown
c hdz        = height of each layer 
c h          = ground level height
c deltah     = computational domain height
c tlon, tlat = latitude and longitude 
c
c********************************************************************
      implicit none
      real, dimension(:,:,:,:), pointer :: sg1, sl1, sp1
      real, dimension(:,:,:),   pointer :: u, v, w
      real, dimension(:,:,:),   pointer :: kh, kv, t
      real, dimension(:,:,:),   pointer :: dz
      real, dimension(:,:,:),   pointer :: wc, wr
      real, dimension(:,:),     pointer :: sprc
      real, dimension(:,:,:),   pointer :: rvel
      real, dimension(:,:,:,:), pointer :: sx, sy
      real, dimension(:,:,:),   pointer :: sz
      real, dimension(:,:,:),   pointer :: q
      real, dimension(:,:,:,:), pointer :: em
      real, dimension(:,:,:),   pointer :: vg, fz
      real, dimension(:,:,:),   pointer :: hdz, cldod, ccover
      real, dimension(:,:),     pointer :: h, deltah, tlon, tlat, dobson 

      real,dimension(:,:),   pointer :: kctop
      integer :: N_x, N_y, N_z, N_gas, N_liquid, N_particle
      integer :: ierr
c
c-----dimension sg1(ixm*iym*izm*iLm)
      allocate( sg1(N_x, N_y, N_z, N_gas), STAT=ierr)
      call AllocErrorCheck(ierr,"sg1")
c
c-----dimension sl1(ixm*iym*izm*iLm)
c      allocate( sl1(N_x, N_y, N_z, N_liquid), STAT=ierr)
      allocate( sl1(1,1,1,1), STAT=ierr)
      call AllocErrorCheck(ierr,"sl1")
c
c-----dimension sp1(ixm*iym*izm*ilptm)
c      allocate( sp1(N_x, N_y, N_z, N_particle), STAT=ierr)
      allocate( sp1(1,1,1,1), STAT=ierr)
      call AllocErrorCheck(ierr,"sp1")
c
c-----dimension u(ixm*iym*izm)
      allocate( u(N_x, N_y, N_z), STAT=ierr)
      call AllocErrorCheck(ierr,"u")
c
c-----dimension v(ixm*iym*izm)
      allocate( v(N_x, N_y, N_z), STAT=ierr)
      call AllocErrorCheck(ierr,"v")
c
c-----dimension w(ixm*iym*izm)
      allocate( w(N_x, N_y, N_z), STAT=ierr)
      call AllocErrorCheck(ierr,"w")
c
c-----dimension kh(ixm*iym*izm)
      allocate( kh(N_x, N_y, N_z), STAT=ierr)
      call AllocErrorCheck(ierr,"kh")
c
c-----dimension kv(ixm*iym*izm)
      allocate( kv(N_x, N_y, N_z), STAT=ierr)
      call AllocErrorCheck(ierr,"kv")
c
c-----dimension t(ixm*iym*izm)
      allocate( t(N_x, N_y, N_z), STAT=ierr)
      call AllocErrorCheck(ierr,"t")
c
c-----dimension dz(ixm*iym*izm)
      allocate( dz(N_x, N_y, N_z), STAT=ierr)
      call AllocErrorCheck(ierr,"dz")
c
c-----dimension wc(ixm*iym*izm)
      allocate( wc(N_x, N_y, N_z), STAT=ierr)
      call AllocErrorCheck(ierr,"wc")
c
c-----dimension wr(ixm*iym*izm)
      allocate( wr(N_x, N_y, N_z), STAT=ierr)
      call AllocErrorCheck(ierr,"wr")
c
c-----dimension sprc(ixm*iym)
      allocate( sprc(N_x, N_y), STAT=ierr)
      call AllocErrorCheck(ierr,"sprc")
c
c-----dimension rvel(ixm*iym*izm)
      allocate( rvel(N_x, N_y, N_z), STAT=ierr)
      call AllocErrorCheck(ierr,"rvel")
c
c-----dimension sx(iym*izm*2*iLm)          
      allocate( sx(N_y, N_z, 2, N_gas), STAT=ierr)
      call AllocErrorCheck(ierr,"sx")
c
c-----dimension sy(ixm*izm*2*iLm)          
      allocate( sy(N_x, N_z, 2, N_gas), STAT=ierr)
      call AllocErrorCheck(ierr,"sy")
c
c-----dimension sz(ixm*iym*iLm)            
      allocate( sz(N_x, N_y, N_gas), STAT=ierr)
      call AllocErrorCheck(ierr,"sz")
c
c-----dimension q(ixm*iym*iLm)
      allocate( q(N_x, N_y, N_gas), STAT=ierr)
      call AllocErrorCheck(ierr,"q")
c
c-----dimension em(ixm*iym*izm*iLm)
      allocate( em(N_x, N_y, N_z, N_gas), STAT=ierr)
      call AllocErrorCheck(ierr,"em")
c
c-----dimension vg(ixm*iym*iLm)
      allocate( vg(N_x, N_y, N_gas), STAT=ierr)
      call AllocErrorCheck(ierr,"vg")
c
c-----dimension fz(ixm*iym*iLm)
      allocate( fz(N_x, N_y, N_gas), STAT=ierr)
      call AllocErrorCheck(ierr,"fz")
c
c-----dimension hdz(ixm*iym*izm)
      allocate( hdz(N_x, N_y, N_z), STAT=ierr)
      call AllocErrorCheck(ierr,"hdz")

c-----dimension cldod(ixm*iym*izm)
      allocate( cldod(N_x, N_y, N_z), STAT=ierr)
      call AllocErrorCheck(ierr,"cldod")

c-----dimension ccover(ixm*iym*izm)
      allocate( ccover(N_x, N_y, N_z), STAT=ierr)
      call AllocErrorCheck(ierr,"cldod")

c
c-----dimension h(ixm*iym)
      allocate( h(N_x, N_y), STAT=ierr)
      call AllocErrorCheck(ierr,"h")
c
c-----dimension deltah(ixm*iym)
      allocate( deltah(N_x, N_y), STAT=ierr)
      call AllocErrorCheck(ierr,"deltah")
c
c-----dimension tlon(ixm*iym)
      allocate( Tlon(N_x, N_y), STAT=ierr)
      call AllocErrorCheck(ierr,"Tlon")
c
c-----dimension tlat(ixm*iym)
      allocate( Tlat(N_x, N_y), STAT=ierr)
      call AllocErrorCheck(ierr,"Tlat")

c-----dimension kctop(ixm*iym)
      allocate( Kctop(N_x, N_y), STAT=ierr)
      call AllocErrorCheck(ierr,"Kctop")

c-----dimension dobson(ixm*iym)
      allocate( dobson(N_x, N_y), STAT=ierr)
      call AllocErrorCheck(ierr,"Kctop")
c
      return
c
      contains
c
        subroutine AllocErrorCheck(ierr,s)
	  integer :: ierr
	  character(LEN=*) :: s
          if (ierr .NE. 0) then
            print*,"Error in MemAlloc: allocation for ",s," failed"
	    stop
          end if 
	end subroutine AllocErrorCheck
c
      end subroutine MemAlloc


      end module StemMemAlloc


