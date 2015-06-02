      module XYStemCommunication
      
      use XYParallelCommunication
      
      contains
      
      subroutine distrib_conc_xy( ix,iy,iz,ixloc,iyloc,
     &                   N_gas, N_liquid, N_particle,
     &                   sg1, sg1_x, sg1_y,
     &                   sl1, sl1_x, sl1_y,
     &                   sp1, sp1_x, sp1_y)
c      
      implicit none
c      
      integer :: ix, iy, iz, ixloc, iyloc  
      integer :: N_gas, N_liquid, N_particle 
c         
      real, dimension(ix,iy,iz,N_gas)         :: sg1
      real, dimension(ix,iyloc,iz,N_gas)      :: sg1_x
      real, dimension(ixloc,iy,iz,N_gas)      :: sg1_y
      real, dimension(ix,iy,iz,N_liquid)      :: sl1
      real, dimension(ix,iyloc,iz,N_liquid)   :: sl1_x
      real, dimension(ixloc,iy,iz,N_liquid)   :: sl1_y
      real, dimension(ix,iy,iz,N_particle)    :: sp1
      real, dimension(ix,iyloc,iz,N_particle) :: sp1_x
      real, dimension(ixloc,iy,iz,N_particle) :: sp1_y
     
c 
      call distrib_x_4D(ix,iy,iz,iyloc,N_gas,sg1,sg1_x)
c  (Not Needed )    call distrib_x_4D(ix,iy,iz,iyloc,N_liquid,sl1,sl1_x)
c  (Not Needed)     call distrib_x_4D(ix,iy,iz,iyloc,N_particle,sp1,sp1_x)
     
      end subroutine distrib_conc_xy


      subroutine distrib_topo_xy( ix,iy,iz,ixloc,iyloc,
     &                   hdz, hdz_x, hdz_y, 
     &                   h,   h_x,  h_y, 
     &                   deltah, deltah_x, deltah_y,
     &                   tlon, tlon_x, tlon_y,
     &                   tlat, tlat_x, tlat_y, 
     &                   dz, dz_x, dz_y )
     
      implicit none
c      
      integer :: ix, iy, iz, ixloc, iyloc  
     
      real, dimension(ix,iy,iz)    :: hdz, dz
      real, dimension(ix,iyloc,iz) :: hdz_x, dz_x
      real, dimension(ixloc,iy,iz) :: hdz_y, dz_y
      real, dimension(ix,iy)    :: h, deltah, tlon, tlat
      real, dimension(ix,iyloc) :: h_x, deltah_x, tlon_x,tlat_x
      real, dimension(ixloc,iy) :: h_y, deltah_y, tlon_y,tlat_y
c     
      call distrib_x_3D(ix,iy,iz,iyloc,hdz,hdz_x)
      call distrib_y_3D(ix,iy,iz,ixloc,hdz,hdz_y)
     
      call distrib_x_2D(ix,iy,iyloc,h,h_x)
      call distrib_y_2D(ix,iy,ixloc,h,h_y)
      call distrib_x_2D(ix,iy,iyloc,deltah,deltah_x)
      call distrib_y_2D(ix,iy,ixloc,deltah,deltah_y)
      call distrib_x_2D(ix,iy,iyloc,tlon,tlon_x)
      call distrib_y_2D(ix,iy,ixloc,tlon,tlon_y)
      call distrib_x_2D(ix,iy,iyloc,tlat,tlat_x)
      call distrib_y_2D(ix,iy,ixloc,tlat,tlat_y)

      call distrib_x_3D(ix,iy,iz,iyloc,dz,dz_x)
      call distrib_y_3D(ix,iy,iz,ixloc,dz,dz_y)
      
      end subroutine distrib_topo_xy     




      subroutine distrib_xy( ix,iy,iz,ixloc,iyloc,
     &                   N_gas, N_liquid, N_particle,
     &                   u, v, w, u_x, v_x, w_x, u_y, v_y, w_y,
     &                   kh, kv, t, kh_x, kv_x, t_x, kh_y, kv_y, t_y,
     &                   wc, wc_x, wc_y, wr, wr_x, wr_y, 
     &                   rvel, rvel_x, rvel_y,
     &                   q, q_x, q_y, sprc, sprc_x, sprc_y,    
     &                   em, em_x, em_y,     
     &                   vg, vg_x, vg_y, fz, fz_x, fz_y,     
     &                   sx, sx_x, 
     &                   sy, sy_y, 
     &                   sz, sz_x, sz_y,
     &                   cldod, cldod_x, cldod_y,
     &                   kctop, kctop_x, kctop_y,
     &                   ccover, ccover_x, ccover_y,
     &                   dobson, dobson_x, dobson_y )
     
c  Master distributes data in x-slice format to all workers
c  A. Sandu, Jan 2001
c-----------------------------------------------------------------------------------------------------------------------------
      implicit none
c      
      integer :: ix, iy, iz, ixloc, iyloc  
      integer :: N_gas, N_liquid, N_particle          
c
      real, dimension(ix,iy,iz,N_gas)    :: em
      real, dimension(ix,iyloc,iz,N_gas) :: em_x
      real, dimension(ixloc,iy,iz,N_gas) :: em_y
c
      real, dimension(ix,iy,iz)    :: u, v, w, kh, kv, t, wc, wr,
     &                             rvel, cldod,  ccover
      real, dimension(ix,iyloc,iz) :: u_x, v_x, w_x, kh_x, kv_x,t_x,
     &       wc_x, wr_x, rvel_x, cldod_x, ccover_x
      real, dimension(ixloc,iy,iz) :: u_y, v_y, w_y, kh_y, kv_y, t_y, 
     &       wc_y, wr_y, rvel_y, cldod_y, ccover_y

c      
      real, dimension(ix,iy,N_gas)    :: vg,   fz,   q
      real, dimension(ix,iyloc,N_gas) :: vg_x, fz_x, q_x
      real, dimension(ixloc,iy,N_gas) :: vg_y, fz_y, q_y
c      
      real, dimension(ix,iy)    :: sprc,   dobson,   kctop
      real, dimension(ix,iyloc) :: sprc_x, dobson_x, kctop_x
      real, dimension(ixloc,iy) :: sprc_y, dobson_y, kctop_y
c      
      real, dimension(iy,iz,2,N_gas)    :: sx
      real, dimension(iyloc,iz,2,N_gas) :: sx_x
      real, dimension(ix,iz,2,N_gas)    :: sy
      real, dimension(ixloc,iz,2,N_gas) :: sy_y
      real, dimension(ix,iy,N_gas)    :: sz
      real, dimension(ix,iyloc,N_gas) :: sz_x
      real, dimension(ixloc,iy,N_gas) :: sz_y

c      external  distrib_y_3d

c 
c      call distrib_x_4D(ix,iy,iz,iyloc,N_gas,sg1,sg1_x)
c  (Not Needed )    call distrib_x_4D(ix,iy,iz,iyloc,N_liquid,sl1,sl1_x)
c  (Not Needed)     call distrib_x_4D(ix,iy,iz,iyloc,N_particle,sp1,sp1_x)
c       
      call distrib_x_3D(ix,iy,iz,iyloc,u,u_x)
c  (Not Needed)      call distrib_x_3D(ix,iy,iz,iyloc,v,v_x)
      call distrib_x_3D(ix,iy,iz,iyloc,w,w_x)
c  (Not Needed)      call distrib_y_3D(ix,iy,iz,ixloc,u,u_y)
      call distrib_y_3D(ix,iy,iz,ixloc,v,v_y)
      call distrib_y_3D(ix,iy,iz,ixloc,w,w_y)
c
      call distrib_x_3D(ix,iy,iz,iyloc,kh,kh_x)
      call distrib_x_3D(ix,iy,iz,iyloc,kv,kv_x)
      call distrib_x_3D(ix,iy,iz,iyloc,t,t_x)
      call distrib_y_3D(ix,iy,iz,ixloc,kh,kh_y)
      call distrib_y_3D(ix,iy,iz,ixloc,kv,kv_y)
      call distrib_y_3D(ix,iy,iz,ixloc,t,t_y)
c
c  (Not Needed)      call distrib_x_3D(ix,iy,iz,iyloc,wc,wc_x)
c  (Not Needed)      call distrib_x_3D(ix,iy,iz,iyloc,wr,wr_x)
      call distrib_y_3D(ix,iy,iz,ixloc,wc,wc_y)
      call distrib_y_3D(ix,iy,iz,ixloc,wr,wr_y)
c
c  (Not Needed)     call distrib_x_3D(ix,iy,iz,iyloc,rvel,rvel_x)
      call distrib_y_3D(ix,iy,iz,ixloc,rvel,rvel_y)
c 
c  (Not Needed)   call distrib_x_4D(ix,iy,iz,iyloc,N_gas,em,em_x)
      call distrib_y_4D(ix,iy,iz,ixloc,N_gas,em,em_y)
c
c   (Not Needed)    call distrib_x_BDz(ix,iy,iyloc,N_gas,vg,vg_x)
      call distrib_y_BDz(ix,iy,ixloc,N_gas,vg,vg_y)
c   (Not Needed)      call distrib_x_BDz(ix,iy,iyloc,N_gas,fz,fz_x)
      call distrib_y_BDz(ix,iy,ixloc,N_gas,fz,fz_y)
      call distrib_x_BDz(ix,iy,iyloc,N_gas,q,q_x)
      call distrib_y_BDz(ix,iy,ixloc,N_gas,q,q_y)
c
C  (not needed)      call distrib_x_3D(ix,iy,iz,iyloc,cldod,cldod_x)
      call distrib_y_3D(ix,iy,iz,ixloc,cldod,cldod_y)
C  (not needed)      call distrib_x_3D(ix,iy,iz,iyloc,ccover,ccover_x)
      call distrib_y_3D(ix,iy,iz,ixloc,ccover,ccover_y)
c      
      call distrib_x_2D(ix,iy,iyloc,sprc,sprc_x)
      call distrib_y_2D(ix,iy,ixloc,sprc,sprc_y)
c
C (not needed)      call distrib_x_2D_INT(ix,iy,iyloc,kctop,kctop_x)
      call distrib_y_2D(ix,iy,ixloc,kctop,kctop_y)
c
      call distrib_x_2D(ix,iy,iyloc,dobson,dobson_x)
      call distrib_y_2D(ix,iy,ixloc,dobson,dobson_y)
c    
      call distrib_x_BD(iy,iz,iyloc,N_gas,sx,sx_x)
      call distrib_y_BD(ix,iz,ixloc,N_gas,sy,sy_y)
      call distrib_x_BDz(ix,iy,iyloc,N_gas,sz,sz_x)
      call distrib_y_BDz(ix,iy,ixloc,N_gas,sz,sz_y)

c 
      end subroutine distrib_xy
 
c-----------------------------------------------------------------------------------------------------------------------------
      subroutine int_distrib(numl,nbin,ixtrn,iytrn,iztrn,
     &                     irxng,irxnl,num,mdt,idate,iend) 
c----------------------
c Master distributes several integer variables to all workers
c  A. Sandu, Jan 2001
c-----------------------------------------------------------------------------------------------------------------------------
         implicit none
      include 'mpif.h'
c      
      integer ::  numl(3,4)
      integer :: nbin,ixtrn,iytrn,iztrn,irxng,irxnl,num,mdt
      integer :: idate(3), iend
      integer, parameter ::  Nbuf=24   
      integer ::  i, j, buf(Nbuf), Ierr, rc    
c
      print*, 'P[',MyId,'] in int_distrib 1'
      if (Master) then
        do j=1,4
	  do i=1,3
           buf((j-1)*3+i) = numl(i,j)
	  enddo
	enddo  
	buf(13:15) = idate(1:3) 
	buf(16) = nbin
	buf(17) = ixtrn
	buf(18) = iytrn
	buf(19) = iztrn
	buf(20) = irxng
	buf(21) = irxnl  
	buf(22) = num  
	buf(23) = mdt  
	buf(24) = iend
      endif     
c
      print*, 'P[',MyId,'] in int_distrib 2'
      call MPI_BCAST(buf,Nbuf,MPI_INTEGER,0,MPI_COMM_WORLD,Ierr)
      if (Ierr.ne.MPI_SUCCESS) then
        print*,'BCAST failed inside int_distrib'
	call MPI_FINALIZE(rc)
      endif
      print*, 'P[',MyId,'] in int_distrib 3'
c      
      if (XWorker.or.YWorker) then
        do j=1,4
	  do i=1,3
           numl(i,j) = buf((j-1)*3+i)
	  enddo
	enddo  
      print*, 'P[',MyId,'] in int_distrib 4'
	idate(1:3) = buf(13:15)
	nbin  = buf(16)
	ixtrn = buf(17)
	iytrn = buf(18)
	iztrn = buf(19)
	irxng = buf(20)
	irxnl = buf(21)  
	num   = buf(22)  
	mdt   = buf(23)  
	iend  = buf(24)
      endif     
c       
      print*, 'P[',MyId,'] in int_distrib 5'
      end  subroutine int_distrib
 
c-----------------------------------------------------------------------------------------------------------------------------


      subroutine int_distrib1(iend) 
c-----------------------
c  Master distributes another bunch of integers to Workers
c  A. Sandu, Jan 2001
c-----------------------------------------------------------------------------------------------------------------------------
        implicit none
      include 'mpif.h'
c      
      include 'aqmax.param'
      include 'aqindx.cmm'
c      
      integer :: iend
      integer, parameter ::  Nbuf=35
      integer ::  i, j, k, buf(Nbuf), Ierr, status(MPI_STATUS_SIZE)    
c
      if (Master) then
! arguments
	buf(1) = iend
!  /aqspid/ in aqindx.cmm stuff
	buf(2) = iair
	buf(3) = ih2o
	buf(4) = io2
	buf(5) = ico 
	buf(6) = ino2
	buf(7) = iho2
	buf(8) = iso2
	buf(9) = io3
	buf(10)= ich4
	buf(11)= ico2
	buf(12)= ih2
	buf(13)= in2
	buf(14)= itrace
	k=15
	buf(k:k+9) = ispg_idx(1:10); k=k+10
	buf(k:k+9) = ispl_idx(1:10); k=k+10
      endif   !   (Master)
c
      call MPI_BCAST(buf,Nbuf,MPI_INTEGER,0,MPI_COMM_WORLD,Ierr)
c      
      if (XWorker.or.YWorker) then
	iend  = buf(1)
! /aqspid/ in aqindx.cmm stuff
	iair  = buf(2)
	ih2o  = buf(3)
	io2   = buf(4)
	ico   = buf(5)
	ino2  = buf(6)
	iho2  = buf(7) 
	iso2  = buf(8)
	io3   = buf(9)
	ich4  = buf(10)
	ico2  = buf(11)
	ih2   = buf(12)
	in2   = buf(13)
	itrace= buf(14)
	k=15
	ispg_idx(1:10) = buf(k:k+9); k=k+10
	ispl_idx(1:10) = buf(k:k+9); k=k+10
      endif  !    (Worker)
c       
      end  subroutine int_distrib1
      



c-----------------------------------------------------------------------------------------------------------------------------
      subroutine int_distrib2(numl,nbin,ixtrn,iytrn,iztrn,
     &                      irxng,irxnl,num,mdt,idate,iend) 
c -----------------------------
c  Master distributes yet another bunch of integers to Workers
c  A. Sandu, Jan 2001
c-----------------------------------------------------------------------------------------------------------------------------
c      
          implicit none
      include 'mpif.h'
c      
      integer ::  numl(3,4)
      integer :: nbin,ixtrn,iytrn,iztrn,irxng,irxnl,num,mdt
      integer :: idate(3), iend
      integer, parameter ::  Nbuf=24   
      integer ::  i, j, buf(Nbuf), Ierr, rc, status(MPI_STATUS_SIZE)    
c
      if (Master) then
        do j=1,4
	  do i=1,3
            buf((j-1)*3+i) = numl(i,j)
	  enddo
	enddo  
	do i=1,3
            buf(12+i) = idate(i)
	enddo
	buf(16) = nbin
	buf(17) = ixtrn
	buf(18) = iytrn
	buf(19) = iztrn
	buf(20) = irxng
	buf(21) = irxnl  
	buf(22) = num  
	buf(23) = mdt  
	buf(24) = iend
      endif   !   (Master)
c
      call MPI_BCAST(buf,Nbuf,MPI_INTEGER,0,MPI_COMM_WORLD,Ierr)
c      
      if (XWorker.or.YWorker) then
        do j=1,4
	  do i=1,3
           numl(i,j) = buf((j-1)*3+i)
	  enddo
	enddo  
	do i=1,3
           idate(i) = buf(12+i)
	enddo
	nbin  = buf(16)
	ixtrn = buf(17)
	iytrn = buf(18)
	iztrn = buf(19)
	irxng = buf(20)
	irxnl = buf(21)  
	num   = buf(22)  
	mdt   = buf(23)  
	iend  = buf(24)
      endif  !    (Worker)
c       
      end  subroutine int_distrib2
      
c-----------------------------------------------------------------------------------------------------------------------------
      subroutine real_distrib(ix,iy,iz,is,dx,dy,
     &                 sigmaz,dht,baseh,rmw,dt,ut) 
c -----------------------------
c  Master distributes a bunch of reals to Workers
c  A. Sandu, Jan 2001
c-----------------------------------------------------------------------------------------------------------------------------
      implicit none
      include 'mpif.h'
c      
      integer :: ix,iy,iz,is,k
      real   :: dx(ix), dy(iy), sigmaz(iz)
      real   :: dht, baseh, rmw(is), dt, ut
      integer :: Nbuf, Ierr
      real   :: buf(ix+iy+iz+is+4)
c
      Nbuf = ix+iy+iz+is+4
      if (Master) then 
        k=0
	buf(k+1:k+ix)=dx(1:ix) ; k=k+ix
	buf(k+1:k+iy)=dy(1:iy) ; k=k+iy
	buf(k+1:k+iz)=sigmaz(1:iz) ; k=k+iz
	buf(k+1)=dht ; k=k+1
	buf(k+1)=baseh ; k=k+1	
	buf(k+1:k+is)=rmw(1:is) ; k=k+is
	buf(k+1)=dt ; k=k+1
	buf(k+1)=ut ; k=k+1
	if (k.ne.Nbuf) then
	  print*, 'Error in real_distrib. Nbuf=',Nbuf,
     &             '   needed=',k
          stop   	  
	endif	
      endif     
c
      call MPI_BCAST(buf(1),Nbuf,MPI_REAL,0,MPI_COMM_WORLD,Ierr)
c      
      if (XWorker.or.YWorker) then
        k=0
	dx(1:ix)=buf(1:ix) ; k=ix
	dy(1:iy)=buf(k+1:k+iy) ; k=k+iy
	sigmaz(1:iz)=buf(k+1:k+iz) ; k=k+iz
	dht=buf(k+1) ; k=k+1
	baseh=buf(k+1) ; k=k+1	
	rmw(1:is)=buf(k+1:k+is); k=k+is
	dt=buf(k+1); k=k+1
	ut=buf(k+1); k=k+1
      endif     
c       
      end  subroutine real_distrib


c-----------------------------------------------------------------------------------------------------------------------------
      subroutine real_distrib2(ix,iy,iz,is,dx,dy,
     &                        sigmaz,dht,baseh,rmw,dt,ut) 
c -----------------------------
c  Master distributes a bunch of reals to Workers
c   (non-broadcast version of real_distrib)
c  A. Sandu, Jan 2001
c-----------------------------------------------------------------------------------------------------------------------------
      implicit none
      include 'mpif.h'
c      
      integer :: ix,iy,iz,is,k
      real   :: dx(ix), dy(iy), sigmaz(iz)
      real   :: dht, baseh, rmw(is), dt, ut
      integer :: Nbuf, Ierr, status(MPI_STATUS_SIZE) 
      real   :: buf(ix+iy+iz+is+4)
      integer :: i
c
      Nbuf = ix+iy+iz+is+4
      if (Master) then 
        k=0
	buf(1:ix)=dx(1:ix) ; k=k+ix
	buf(k+1:k+iy)=dy(1:iy) ; k=k+iy
	buf(k+1:k+iz)=sigmaz(1:iz) ; k=k+iz
	buf(k+1)=dht ; k=k+1
	buf(k+1)=baseh ; k=k+1	
	buf(k+1:k+is)=rmw(1:is) ; k=k+is
	buf(k+1)=dt ; k=k+1
	buf(k+1)=ut ; k=k+1
	if (k.ne.Nbuf) then
	  print*, 'Error in real_distrib. Nbuf=',Nbuf,
     &             '   needed=',k
          stop   	  
	endif	
      endif     
c
      call MPI_BCAST(buf,Nbuf,MPI_REAL,0,MPI_COMM_WORLD,Ierr)
c      
      if (XWorker.or.YWorker) then
        k=0
	dx(1:ix)=buf(k+1:k+ix) ; k=k+ix
	dy(1:iy)=buf(k+1:k+iy) ; k=k+iy
	sigmaz(1:iz)=buf(k+1:k+iz) ; k=k+iz
	dht=buf(k+1) ; k=k+1
	baseh=buf(k+1) ; k=k+1	
	rmw(1:is)=buf(k+1:k+is); k=k+is
	dt=buf(k+1); k=k+1
	ut=buf(k+1); k=k+1
      endif     
c       
      end  subroutine real_distrib2

c-----------------------------------------------------------------------------------------------------------------------------
      subroutine cmm_distrib() 
c -----------------------------
c  Master distributes all common block data to Workers
c               (non-broadcast version)
c  A. Sandu, Jan 2001
c-----------------------------------------------------------------------------------------------------------------------------
      implicit none
      include 'mpif.h'
      include 'aqmax.param'
      !1.include 'aqcon1.cmm'
      real aqcon1,aqcon2 
      common /aqcon1/aqcon1(22*mxspg)      
      common /aqcon2/aqcon2(2000)
      integer iaqcon2
      common /iaqcon2/iaqcon2(22*mxspg+2135)
      !2.include 'aqcon2.cmm'
      double precision aqcon3
      common /aqcon3/aqcon3(3660)
      integer iaqcon3
      common /iaqcon3/iaqcon3(1862)
      !3.include 'aqcon5.cmm'
      double precision aqcon5
      common /aqcon5/aqcon5(13*mxspg+4*mxrxn+1)
      !4.include 'aqcont.cmm'
      integer aqcont, aq_unit
      common /aqcont/aqcont(24)
      common /aq_unit/aq_unit(6)      
      !6.include 'aqindx.cmm'
      integer iaqspid, iaqmatr
      common /aqspid/iaqspid(33)
      common /aqmatr/iaqmatr(mxspg*20+1)
      !include 'aqrxn.cmm'
      double precision aqrxng
      common /aqrxng/  aqrxng(66*mxrxn)
      real raqrxng
      common /raqrxng/  raqrxng(2*mxrxn)
      integer iaqrxng
      common /iaqrxng/ iaqrxng(64*mxrxn+1)
      double precision aqrxnl
      common /aqrxnl/  aqrxnl(63*mxrxn)
      integer iaqrxnl
      common /iaqrxnl/ iaqrxnl(64*mxrxn+1)
      real  aqphoto
      common /aqphoto/  aqphoto(3*30*500)
      integer iaqphoto
      common /iaqphoto/ iaqphoto(61)
      !include 'aqspec.cmm'
      character specinfo
      common /specinfo/specinfo(16*4*mxspg)
      integer ispecinfo
      integer, parameter :: mxal=5000
      common /ispecinfo/ispecinfo(4) 
      real  rspecinfo
      common /rspecinfo/rspecinfo(4*mxal)
      !include 'aqsymb.cmm'
      character aqsymb
      common /aqsymb/aqsymb( 16*(3*mxspg+1000) )
                 
c      
      integer :: Ierr, status(MPI_STATUS_SIZE) 
      integer :: Ndbuf, Nrbuf, Nibuf, Ncbuf 
      integer :: kd, kr, ki, kc, i
      double precision  :: dbuf(3660+13*mxspg+4*mxrxn+1+66*mxrxn
     &                      +63*mxrxn)
      real         :: rbuf(22*mxspg+2*mxrxn+2000
     &                      +4*mxal+3*30*500)
      integer       :: ibuf(22*mxspg+2135+1862+24+6
     &                      +33+mxspg*20+1+64*mxrxn
     &                      +1+64*mxrxn+1+61+4)
      character      :: cbuf(16*4*mxspg+16*(3*mxspg+1000))
c
      Ndbuf=+3660+13*mxspg+4*mxrxn+1+66*mxrxn+63*mxrxn 
      Nrbuf=+22*mxspg+2*mxrxn+2000+4*mxal+3*30*500
      Nibuf=+22*mxspg+2135+1862+24+6+33+mxspg*20+1
     &        +64*mxrxn+1+64*mxrxn+1+61+4 
      Ncbuf=+16*4*mxspg+16*(3*mxspg+1000) 
c
      if (Master) then 
          kd=0; kr=0; ki=0; kc=0;
c  aqcon1.cmm  
          rbuf(kr+1:kr+22*mxspg)=aqcon1(1:22*mxspg); 
	  kr=kr+22*mxspg     
          rbuf(kr+1:kr+2000)=aqcon2(1:2000)    ;     
	  kr=kr+2000
          ibuf(ki+1:ki+22*mxspg+2135)=
     &              iaqcon2(1:22*mxspg+2135); 
	         ki=ki+22*mxspg+2135
c  aqcon2.cmm  
          dbuf(kd+1:kd+3660)=aqcon3(1:3660);  kd=kd+3660
          ibuf(ki+1:ki+1862)=iaqcon3(1:1862); ki=ki+1862	  
c  aqcon5.cmm  
          dbuf(kd+1:kd+13*mxspg+4*mxrxn+1)=aqcon5(1:13*mxspg+4*mxrxn+1); 
	  kd=kd+13*mxspg+4*mxrxn+1
c  aqcont.cmm  
          ibuf(ki+1:ki+24)=aqcont(1:24); ki=ki+24
          ibuf(ki+1:ki+6)=aq_unit(1:6);  ki=ki+6
c  aqindx.cmm  
          ibuf(ki+1:ki+33)=iaqspid(1:33); ki=ki+33
          ibuf(ki+1:ki+mxspg*20+1)=iaqmatr(1:mxspg*20+1); 
	  ki=ki+mxspg*20+1
c  aqrxn.cmm
          dbuf(kd+1:kd+66*mxrxn)=aqrxng(1:66*mxrxn); 
	  kd=kd+66*mxrxn
          rbuf(kr+1:kr+2*mxrxn)=raqrxng(1:2*mxrxn); 
	  kr=kr+2*mxrxn
          ibuf(ki+1:ki+64*mxrxn+1)=iaqrxng(1:64*mxrxn+1); 
	  ki=ki+64*mxrxn+1
          dbuf(kd+1:kd+63*mxrxn)=aqrxnl(1:63*mxrxn); 
	  kd=kd+63*mxrxn
          ibuf(ki+1:ki+64*mxrxn+1)=iaqrxnl(1:64*mxrxn+1); 
	  ki=ki+64*mxrxn+1
          rbuf(kr+1:kr+3*30*500)=aqphoto(1:3*30*500); 
	  kr=kr+3*30*500
          ibuf(ki+1:ki+61)=iaqphoto(1:61); ki=ki+61
c   aqspec.cmm 
          cbuf(kc+1:kc+16*4*mxspg)=specinfo(1:16*4*mxspg); 
	  kc=kc+16*4*mxspg
          ibuf(ki+1:ki+4)=ispecinfo(1:4); ki=ki+4 
          rbuf(kr+1:kr+4*mxal)=rspecinfo(1:4*mxal); kr=kr+4*mxal
c   qsymb.cmm
          cbuf(kc+1:kc+16*(3*mxspg+1000))=
     &            aqsymb(1:16*(3*mxspg+1000));
	        kc=kc+16*(3*mxspg+1000)	  
c
	if (kd.ne.Ndbuf) then
	  print*, 'Error in double cmm_distrib. Nbuf=',Ndbuf,
     &             '   needed=',kd
          stop   	  
	else if (kr.ne.Nrbuf) then
	  print*, 'Error in real cmm_distrib. Nbuf=',Nrbuf,
     &             '   needed=',kr
          stop   	  
	else if (ki.ne.Nibuf) then
	  print*, 'Error in integer cmm_distrib. Nbuf=',Nibuf,
     &             '   needed=',ki
          stop   	  
	else if (kc.ne.Ncbuf) then
	  print*, 'Error in character cmm_distrib. Nbuf=',Ncbuf,
     &             '   needed=',kc
          stop   	  
	endif
		
      endif     
c
      call MPI_BCAST(dbuf,Ndbuf,MPI_DOUBLE_PRECISION,0,
     &                  MPI_COMM_WORLD,Ierr)
      call MPI_BCAST(rbuf,Nrbuf,MPI_REAL,0,MPI_COMM_WORLD,Ierr)
      call MPI_BCAST(ibuf,Nibuf,MPI_INTEGER,0,MPI_COMM_WORLD,Ierr)
      call MPI_BCAST(cbuf,Ncbuf,MPI_CHARACTER,0,MPI_COMM_WORLD,Ierr)

c      
      if (XWorker.or.YWorker) then
c
          kd=0; kr=0; ki=0; kc=0;
c  aqcon1.cmm  
          aqcon1(1:22*mxspg)=rbuf(kr+1:kr+22*mxspg); kr=kr+22*mxspg     
          aqcon2(1:2000)=rbuf(kr+1:kr+2000)    ;     kr=kr+2000
          iaqcon2(1:22*mxspg+2135)=ibuf(ki+1:ki+22*mxspg+2135); 
	                                  ki=ki+22*mxspg+2135
c  aqcon2.cmm  
          aqcon3(1:3660)=dbuf(kd+1:kd+3660);  kd=kd+3660
          iaqcon3(1:1862)=ibuf(ki+1:ki+1862); ki=ki+1862	  
c  aqcon5.cmm  
          aqcon5(1:13*mxspg+4*mxrxn+1)=dbuf(kd+1:kd+13*mxspg+4*mxrxn+1); 
	  kd=kd+13*mxspg+4*mxrxn+1
c  aqcont.cmm  
          aqcont(1:24)=ibuf(ki+1:ki+24); ki=ki+24
          aq_unit(1:6)=ibuf(ki+1:ki+6);  ki=ki+6
c  aqindx.cmm  
          iaqspid(1:33)=ibuf(ki+1:ki+33); ki=ki+33
          iaqmatr(1:mxspg*20+1)=ibuf(ki+1:ki+mxspg*20+1); 
	  ki=ki+mxspg*20+1
c  aqrxn.cmm
          aqrxng(1:66*mxrxn)=dbuf(kd+1:kd+66*mxrxn); 
	  kd=kd+66*mxrxn
          raqrxng(1:2*mxrxn)=rbuf(kr+1:kr+2*mxrxn); 
	  kr=kr+2*mxrxn
          iaqrxng(1:64*mxrxn+1)=ibuf(ki+1:ki+64*mxrxn+1); 
	  ki=ki+64*mxrxn+1
          aqrxnl(1:63*mxrxn)=dbuf(kd+1:kd+63*mxrxn); 
	  kd=kd+63*mxrxn
          iaqrxnl(1:64*mxrxn+1)=ibuf(ki+1:ki+64*mxrxn+1); 
	  ki=ki+64*mxrxn+1
          aqphoto(1:3*30*500)=rbuf(kr+1:kr+3*30*500); 
	  kr=kr+3*30*500
          iaqphoto(1:61)=ibuf(ki+1:ki+61); ki=ki+61
c   aqspec.cmm 
          specinfo(1:16*4*mxspg)=cbuf(kc+1:kc+16*4*mxspg); 
	  kc=kc+16*4*mxspg
          ispecinfo(1:4)=ibuf(ki+1:ki+4); ki=ki+4 
          rspecinfo(1:4*mxal)=rbuf(kr+1:kr+4*mxal); 
	  kr=kr+4*mxal
c   qsymb.cmm
          aqsymb(1:16*(3*mxspg+1000))=
     &           cbuf(kc+1:kc+16*(3*mxspg+1000));
	       kc=kc+16*(3*mxspg+1000)	  
c
      endif !  MyId.ne.0   
c           
      end  subroutine cmm_distrib


      subroutine cmm_distrib2(Nworkers)
c -----------------------------
c  Master distributes all common block data to Workers
c               (non-broadcast version)
c  A. Sandu, Jan 2001
c-----------------------------------------------------------------------------------------------------------------------------
      implicit none
      include 'mpif.h'
      include 'aqmax.param'
      !1.include 'aqcon1.cmm'
      integer Nworkers
      real aqcon1,aqcon2
      common /aqcon1/aqcon1(22*mxspg)
      common /aqcon2/aqcon2(2000)
      integer iaqcon2
      common /iaqcon2/iaqcon2(22*mxspg+2135)
      !2.include 'aqcon2.cmm'
      double precision aqcon3
      common /aqcon3/aqcon3(3660)
      integer iaqcon3
      common /iaqcon3/iaqcon3(1862)
      !3.include 'aqcon5.cmm'
      double precision aqcon5
      common /aqcon5/aqcon5(13*mxspg+4*mxrxn+1)
      !4.include 'aqcont.cmm'
      integer aqcont, aq_unit
      common /aqcont/aqcont(24)
      common /aq_unit/aq_unit(6)
      !6.include 'aqindx.cmm'
      integer iaqspid, iaqmatr
      common /aqspid/iaqspid(33)
      common /aqmatr/iaqmatr(mxspg*20+1)
      !include 'aqrxn.cmm'
      double precision aqrxng
      common /aqrxng/  aqrxng(66*mxrxn)
      real raqrxng
      common /raqrxng/  raqrxng(2*mxrxn)
      integer iaqrxng
      common /iaqrxng/ iaqrxng(64*mxrxn+1)
      double precision aqrxnl
      common /aqrxnl/  aqrxnl(63*mxrxn)
      integer iaqrxnl
      common /iaqrxnl/ iaqrxnl(64*mxrxn+1)
      real  aqphoto
      common /aqphoto/  aqphoto(3*30*500)
      integer iaqphoto
      common /iaqphoto/ iaqphoto(61)
      !include 'aqspec.cmm'
      character specinfo
      common /specinfo/specinfo(16*4*mxspg)
      integer ispecinfo
      integer, parameter :: mxal=5000
      common /ispecinfo/ispecinfo(4)
      real  rspecinfo
      common /rspecinfo/rspecinfo(4*mxal)
      !include 'aqsymb.cmm'
      character aqsymb
      common /aqsymb/aqsymb( 16*(3*mxspg+1000) )

c      
      integer :: Ierr, status(MPI_STATUS_SIZE)
      integer :: Ndbuf, Nrbuf, Nibuf, Ncbuf
      integer :: kd, kr, ki, kc, i
      double precision  :: dbuf(3660+13*mxspg+4*mxrxn+1+66*mxrxn
     &                      +63*mxrxn)
      real         :: rbuf(22*mxspg+2*mxrxn+2000
     &                      +4*mxal+3*30*500)
      integer       :: ibuf(22*mxspg+2135+1862+24+6
     &                      +33+mxspg*20+1+64*mxrxn
     &                      +1+64*mxrxn+1+61+4)
      character      :: cbuf(16*4*mxspg+16*(3*mxspg+1000))
c
      Ndbuf=+3660+13*mxspg+4*mxrxn+1+66*mxrxn+63*mxrxn
      Nrbuf=+22*mxspg+2*mxrxn+2000+4*mxal+3*30*500
      Nibuf=+22*mxspg+2135+1862+24+6+33+mxspg*20+1
     &        +64*mxrxn+1+64*mxrxn+1+61+4
      Ncbuf=+16*4*mxspg+16*(3*mxspg+1000)
c
      if (Master) then
          kd=0; kr=0; ki=0; kc=0;
c  aqcon1.cmm  
          rbuf(kr+1:kr+22*mxspg)=aqcon1(1:22*mxspg);
          kr=kr+22*mxspg
          rbuf(kr+1:kr+2000)=aqcon2(1:2000)    ;
          kr=kr+2000
          ibuf(ki+1:ki+22*mxspg+2135)=
     &              iaqcon2(1:22*mxspg+2135);
                 ki=ki+22*mxspg+2135
c  aqcon2.cmm  
          dbuf(kd+1:kd+3660)=aqcon3(1:3660);  kd=kd+3660
          ibuf(ki+1:ki+1862)=iaqcon3(1:1862); ki=ki+1862        
c  aqcon5.cmm  
          dbuf(kd+1:kd+13*mxspg+4*mxrxn+1)=aqcon5(1:13*mxspg+4*mxrxn+1); 
          kd=kd+13*mxspg+4*mxrxn+1
c  aqcont.cmm  
          ibuf(ki+1:ki+24)=aqcont(1:24); ki=ki+24
          ibuf(ki+1:ki+6)=aq_unit(1:6);  ki=ki+6
c  aqindx.cmm  
          ibuf(ki+1:ki+33)=iaqspid(1:33); ki=ki+33
          ibuf(ki+1:ki+mxspg*20+1)=iaqmatr(1:mxspg*20+1);
          ki=ki+mxspg*20+1
c  aqrxn.cmm
          dbuf(kd+1:kd+66*mxrxn)=aqrxng(1:66*mxrxn);
          kd=kd+66*mxrxn
          rbuf(kr+1:kr+2*mxrxn)=raqrxng(1:2*mxrxn);
          kr=kr+2*mxrxn
          ibuf(ki+1:ki+64*mxrxn+1)=iaqrxng(1:64*mxrxn+1);
          ki=ki+64*mxrxn+1
          dbuf(kd+1:kd+63*mxrxn)=aqrxnl(1:63*mxrxn);
          kd=kd+63*mxrxn
          ibuf(ki+1:ki+64*mxrxn+1)=iaqrxnl(1:64*mxrxn+1);
          ki=ki+64*mxrxn+1
          rbuf(kr+1:kr+3*30*500)=aqphoto(1:3*30*500);
          kr=kr+3*30*500
          ibuf(ki+1:ki+61)=iaqphoto(1:61); ki=ki+61
c   aqspec.cmm 
          cbuf(kc+1:kc+16*4*mxspg)=specinfo(1:16*4*mxspg);
          kc=kc+16*4*mxspg
          ibuf(ki+1:ki+4)=ispecinfo(1:4); ki=ki+4
          rbuf(kr+1:kr+4*mxal)=rspecinfo(1:4*mxal); kr=kr+4*mxal
c   qsymb.cmm
          cbuf(kc+1:kc+16*(3*mxspg+1000))=
     &            aqsymb(1:16*(3*mxspg+1000));
                kc=kc+16*(3*mxspg+1000) 
c
        if (kd.ne.Ndbuf) then
          print*, 'Error in double cmm_distrib. Nbuf=',Ndbuf,
     &             '   needed=',kd
          stop          
        else if (kr.ne.Nrbuf) then
          print*, 'Error in real cmm_distrib. Nbuf=',Nrbuf,
     &             '   needed=',kr
          stop          
        else if (ki.ne.Nibuf) then
          print*, 'Error in integer cmm_distrib. Nbuf=',Nibuf,
     &             '   needed=',ki
          stop          
        else if (kc.ne.Ncbuf) then
          print*, 'Error in character cmm_distrib. Nbuf=',Ncbuf,
     &             '   needed=',kc
          stop          
        endif
                
        do i=1,Nworkers
          call MPI_SEND(dbuf(1), Ndbuf, MPI_DOUBLE_PRECISION, i,
     &          i,  MPI_COMM_WORLD, Ierr)
          call MPI_SEND(rbuf(1), Nrbuf, MPI_REAL, i,
     &          i,  MPI_COMM_WORLD, Ierr)
          call MPI_SEND(ibuf(1), Nibuf, MPI_INTEGER, i,
     &          i,  MPI_COMM_WORLD, Ierr)
          call MPI_SEND(cbuf(1), Ncbuf, MPI_CHARACTER, i,
     &          i,  MPI_COMM_WORLD, Ierr)
        enddo
      endif
c
c      call MPI_BCAST(buf,Nbuf,MPI_REAL,0,MPI_COMM_WORLD,Ierr)
c      
      if (XWorker.or.YWorker) then
        call MPI_RECV(dbuf, Ndbuf, MPI_DOUBLE_PRECISION, 0, MyId,
     &                 MPI_COMM_WORLD, status, ierr)
        call MPI_RECV(rbuf, Nrbuf, MPI_REAL, 0, MyId,
     &                 MPI_COMM_WORLD, status, ierr)
        call MPI_RECV(ibuf, Nibuf, MPI_INTEGER, 0, MyId,
     &                 MPI_COMM_WORLD, status, ierr)
        call MPI_RECV(cbuf, Ncbuf, MPI_CHARACTER, 0, MyId,
     &                 MPI_COMM_WORLD, status, ierr)
                   kd=0; kr=0; ki=0; kc=0;
c  aqcon1.cmm  
          aqcon1(1:22*mxspg)=rbuf(kr+1:kr+22*mxspg); kr=kr+22*mxspg
          aqcon2(1:2000)=rbuf(kr+1:kr+2000)    ;     kr=kr+2000
          iaqcon2(1:22*mxspg+2135)=ibuf(ki+1:ki+22*mxspg+2135);
                                          ki=ki+22*mxspg+2135
c  aqcon2.cmm  
          aqcon3(1:3660)=dbuf(kd+1:kd+3660);  kd=kd+3660
          iaqcon3(1:1862)=ibuf(ki+1:ki+1862); ki=ki+1862        
c  aqcon5.cmm  
          aqcon5(1:13*mxspg+4*mxrxn+1)=dbuf(kd+1:kd+13*mxspg+4*mxrxn+1); 
          kd=kd+13*mxspg+4*mxrxn+1
c  aqcont.cmm  
          aqcont(1:24)=ibuf(ki+1:ki+24); ki=ki+24
          aq_unit(1:6)=ibuf(ki+1:ki+6);  ki=ki+6
c  aqindx.cmm  
          iaqspid(1:33)=ibuf(ki+1:ki+33); ki=ki+33
          iaqmatr(1:mxspg*20+1)=ibuf(ki+1:ki+mxspg*20+1);
          ki=ki+mxspg*20+1
c  aqrxn.cmm
          aqrxng(1:66*mxrxn)=dbuf(kd+1:kd+66*mxrxn);
          kd=kd+66*mxrxn
          raqrxng(1:2*mxrxn)=rbuf(kr+1:kr+2*mxrxn);
          kr=kr+2*mxrxn
          iaqrxng(1:64*mxrxn+1)=ibuf(ki+1:ki+64*mxrxn+1);
          ki=ki+64*mxrxn+1
          aqrxnl(1:63*mxrxn)=dbuf(kd+1:kd+63*mxrxn);
          kd=kd+63*mxrxn
          iaqrxnl(1:64*mxrxn+1)=ibuf(ki+1:ki+64*mxrxn+1);
          ki=ki+64*mxrxn+1
          aqphoto(1:3*30*500)=rbuf(kr+1:kr+3*30*500);
          kr=kr+3*30*500
          iaqphoto(1:61)=ibuf(ki+1:ki+61); ki=ki+61
c   aqspec.cmm 
          specinfo(1:16*4*mxspg)=cbuf(kc+1:kc+16*4*mxspg);
          kc=kc+16*4*mxspg
          ispecinfo(1:4)=ibuf(ki+1:ki+4); ki=ki+4
          rspecinfo(1:4*mxal)=rbuf(kr+1:kr+4*mxal);
          kr=kr+4*mxal
c   qsymb.cmm
          aqsymb(1:16*(3*mxspg+1000))=
     &           cbuf(kc+1:kc+16*(3*mxspg+1000));
               kc=kc+16*(3*mxspg+1000)  
      endif !  MyId.ne.0   
c           
      end  subroutine cmm_distrib2



      end module XYStemCommunication




      module HVStemCommunication
      
      use HVParallelCommunication
           

      contains
 
      subroutine distrib_conc_hv( ix,iy,iz,izloc,icloc,
     &                   N_gas, N_liquid, N_particle,
     &                   sg1, sg1_h, sg1_v,
     &                   sl1, sl1_h, sl1_v,
     &                   sp1, sp1_h, sp1_v)
c      
      implicit none
c      
c      
      integer :: ix, iy, iz, izloc, icloc  
      integer :: N_gas, N_liquid, N_particle 
c         
      real, dimension(ix,iy,iz,N_gas)         :: sg1
      real, dimension(ix,iy,izloc,N_gas)      :: sg1_h
      real, dimension(1,icloc,iz,N_gas)       :: sg1_v
      real, dimension(ix,iy,iz,N_liquid)      :: sl1
      real, dimension(ix,iy,izloc,N_liquid)   :: sl1_h
      real, dimension(1,icloc,iz,N_liquid)    :: sl1_v
      real, dimension(ix,iy,iz,N_particle)    :: sp1
      real, dimension(ix,iy,izloc,N_particle) :: sp1_h
      real, dimension(1,icloc,iz,N_particle)  :: sp1_v
     
c 
      call distrib_h_4D(ix,iy,iz,izloc,N_gas,sg1,sg1_h)
c  (Not Needed )    call distrib_h_4D(ix,iy,iz,izloc,N_liquid,sl1,sl1_h)
c  (Not Needed)     call distrib_h_4D(ix,iy,iz,izloc,N_particle,sp1,sp1_h)
     
      end subroutine distrib_conc_hv

      
c-------------------------------------------------------------------------------------      

      subroutine distrib_topo_hv( ix,iy,iz,izloc,icloc,
     &                   hdz, hdz_h, hdz_v, 
     &                   h,   h_v, 
     &                   deltah, deltah_v,
     &                   tlon, tlon_v,
     &                   tlat, tlat_v, 
     &                   dz, dz_h, dz_v )
     
      implicit none
c      
      integer :: ix, iy, iz, izloc, icloc  
     

      real, dimension(ix,iy,iz)    :: hdz,   dz
      real, dimension(ix,iy,izloc) :: hdz_h, dz_h
      real, dimension(1,icloc,iz)  :: hdz_v, dz_v

      real, dimension(ix,iy)    :: h, deltah, tlon, tlat
      real, dimension(1,icloc) :: h_v, deltah_v, tlon_v, tlat_v
c     

      call distrib_v_3D(ix,iy,iz,icloc,hdz,hdz_v)
      call distrib_v_3D(ix,iy,iz,icloc,dz,dz_v)
      call distrib_h_3D(ix,iy,iz,izloc,hdz,hdz_h)
      call distrib_h_3D(ix,iy,iz,izloc,dz,dz_h)

      call distrib_v_2D(ix,iy,icloc,h,h_v)
      !call distrib_h_2D(ix,iy,h)
      
      call distrib_v_2D(ix,iy,icloc,deltah,deltah_v)
      !call distrib_h_2D(ix,iy,deltah)
      
      call distrib_v_2D(ix,iy,icloc,tlon,tlon_v)
      !call distrib_h_2D(ix,iy,tlon)
     
      call distrib_v_2D(ix,iy,icloc,tlat,tlat_v)
      !call distrib_h_2D(ix,iy,tlat)
     
      end subroutine distrib_topo_hv     


c-----------------------------------------------------------------------------------------------------------------------------
      subroutine distrib_hv( ix,iy,iz,izloc,icloc,
     &                   N_gas, N_liquid, N_particle,
     &                   u, v, w, u_h, v_h, w_v, 
     &                   kh, kv, t, kh_h, kv_h, t_h, kh_v, kv_v, t_v,
     &                   wc, wc_h, wc_v, wr, wr_h, wr_v,
     &                   rvel, rvel_h, rvel_v,
     &                   q, q_v, sprc, sprc_v,    
     &                   em, em_h, em_v,
     &                   vg, vg_h, vg_v, fz, fz_h, fz_v,     
     &                   sx, sx_h, sy, sy_h, sz, sz_v,
     &                   cldod, cldod_h, cldod_v,
     &                   kctop, kctop_h, kctop_v,
     &                   ccover, ccover_h, ccover_v,
     &                   dobson, dobson_h, dobson_v )
c--------
c  Master distributes data in H-slice format to all workers
c-----------------------------------------------------------------------------------------------------------------------------
c      
      implicit none
c      
      integer :: ix, iy, iz, izloc, icloc  
      integer :: N_gas, N_liquid, N_particle 
c         
c
      real, dimension(ix,iy,iz,N_gas)    :: em
      real, dimension(1,icloc,iz,N_gas) :: em_v
      real, dimension(ix,iy,izloc, N_gas) :: em_h
c
      real, dimension(ix,iy,iz)    :: u, v, w, kh, kv, t, wc, wr,
     &                             rvel, hdz,   dz
      real, dimension(ix,iy,izloc) :: u_h, v_h, kh_h, kv_h, t_h,
     &                             wc_h, wr_h, rvel_h
      real, dimension(1,icloc,iz)  :: w_v, kv_v, t_v, kh_v,
     &                             wc_v, wr_v, rvel_v, hdz_v, dz_v

c      
      real, dimension(ix,iy,N_gas)    :: vg,   fz,   sz,   q
      real, dimension(1,icloc,N_gas)  :: vg_v, fz_v, sz_v, q_v
      real, dimension(ix,iy,N_gas)    :: vg_h,fz_h,q_h
c      
      real, dimension(ix,iy)    :: sprc
      real, dimension(1,icloc) :: sprc_v
c      
      real, dimension(iy,iz,2,N_gas)    :: sx
      real, dimension(iy,izloc,2,N_gas) :: sx_h
      real, dimension(ix,iz,2,N_gas)    :: sy
      real, dimension(ix,izloc,2,N_gas) :: sy_h
      
      real, dimension(ix,iy)    :: dobson,   kctop
      real, dimension(ix,iy   ) :: dobson_h, kctop_h
      real, dimension(1,icloc)  :: dobson_v, kctop_v
      real, dimension(ix,iy,iz)    :: cldod,  ccover
      real, dimension(ix,iy,izloc) :: cldod_h, ccover_h
      real, dimension(1,icloc,iz)  :: cldod_v, ccover_v
c      external  distrib_y_3d

c 
c       
      call distrib_h_3D(ix,iy,iz,izloc,u,u_h)
      call distrib_h_3D(ix,iy,iz,izloc,v,v_h)
      call distrib_v_3D(ix,iy,iz,icloc,w,w_v)
c
      call distrib_h_3D(ix,iy,iz,izloc,kh,kh_h)
      call distrib_v_3D(ix,iy,iz,icloc,kv,kv_v)
      call distrib_v_3D(ix,iy,iz,icloc,t,t_v)

      call distrib_v_3D(ix,iy,iz,icloc,wc,wc_v)
      call distrib_v_3D(ix,iy,iz,icloc,wr,wr_v)
c
      call distrib_v_3D(ix,iy,iz,icloc,rvel,rvel_v)
c 
      call distrib_v_4D(ix,iy,iz,icloc,N_gas,em,em_v)
c
      call distrib_v_BDz(ix,iy,icloc,N_gas,vg,vg_v)
      call distrib_v_BDz(ix,iy,icloc,N_gas,fz,fz_v)
      call distrib_v_BDz(ix,iy,icloc,N_gas,sz,sz_v)
      call distrib_v_BDz(ix,iy,icloc,N_gas,q,q_v)
c
      call distrib_v_2D(ix,iy,icloc,sprc,sprc_v)
c
c    
      call distrib_yh_BD(iy,iz,izloc,N_gas,sx,sx_h)
      call distrib_xh_BD(ix,iz,izloc,N_gas,sy,sy_h)
c 
      call distrib_v_3D(ix,iy,iz,icloc,cldod,cldod_v)
      call distrib_v_3D(ix,iy,iz,icloc,ccover,ccover_v)
c
C (not needed)      call distrib_x_2D_INT(ix,iy,iyloc,kctop,kctop_x)
      call distrib_v_2D(ix,iy,icloc,kctop,kctop_v)
c
      call distrib_v_2D(ix,iy,icloc,dobson,dobson_v)
      
      end subroutine distrib_hv
      

 
c-----------------------------------------------------------------------------------------------------------------------------
      subroutine int_distrib(numl,nbin,ixtrn,iytrn,iztrn,
     &                     irxng,irxnl,num,mdt,idate,iend) 
c----------------------
c Master distributes several integer variables to all workers
c  A. Sandu, Jan 2001
c-----------------------------------------------------------------------------------------------------------------------------
         implicit none
      include 'mpif.h'
c      
      integer ::  numl(3,4)
      integer :: nbin,ixtrn,iytrn,iztrn,irxng,irxnl,num,mdt
      integer :: idate(3), iend
      integer, parameter ::  Nbuf=24   
      integer ::  i, j, buf(Nbuf), Ierr, rc    
c
      print*, 'P[',MyId,'] in int_distrib 1'
      if (Master) then
        do j=1,4
	  do i=1,3
           buf((j-1)*3+i) = numl(i,j)
	  enddo
	enddo  
	buf(13:15) = idate(1:3) 
	buf(16) = nbin
	buf(17) = ixtrn
	buf(18) = iytrn
	buf(19) = iztrn
	buf(20) = irxng
	buf(21) = irxnl  
	buf(22) = num  
	buf(23) = mdt  
	buf(24) = iend
      endif     
c
      print*, 'P[',MyId,'] in int_distrib 2'
      call MPI_BCAST(buf,Nbuf,MPI_INTEGER,0,MPI_COMM_WORLD,Ierr)
      if (Ierr.ne.MPI_SUCCESS) then
        print*,'BCAST failed inside int_distrib'
	call MPI_FINALIZE(rc)
      endif
      print*, 'P[',MyId,'] in int_distrib 3'
c      
      if (HWorker.or.VWorker) then
        do j=1,4
	  do i=1,3
           numl(i,j) = buf((j-1)*3+i)
	  enddo
	enddo  
      print*, 'P[',MyId,'] in int_distrib 4'
	idate(1:3) = buf(13:15)
	nbin  = buf(16)
	ixtrn = buf(17)
	iytrn = buf(18)
	iztrn = buf(19)
	irxng = buf(20)
	irxnl = buf(21)  
	num   = buf(22)  
	mdt   = buf(23)  
	iend  = buf(24)
      endif     
c       
      print*, 'P[',MyId,'] in int_distrib 5'
      end  subroutine int_distrib
 
c-----------------------------------------------------------------------------------------------------------------------------


      subroutine int_distrib1(iend) 
c-----------------------
c  Master distributes another bunch of integers to Workers
c  A. Sandu, Jan 2001
c-----------------------------------------------------------------------------------------------------------------------------
        implicit none
      include 'mpif.h'
c      
      include 'aqmax.param'
      include 'aqindx.cmm'
c      
      integer :: iend
      integer, parameter ::  Nbuf=35
      integer ::  i, j, k, buf(Nbuf), Ierr, status(MPI_STATUS_SIZE)    
c
      if (Master) then
! arguments
	buf(1) = iend
!  /aqspid/ in aqindx.cmm stuff
	buf(2) = iair
	buf(3) = ih2o
	buf(4) = io2
	buf(5) = ico 
	buf(6) = ino2
	buf(7) = iho2
	buf(8) = iso2
	buf(9) = io3
	buf(10)= ich4
	buf(11)= ico2
	buf(12)= ih2
	buf(13)= in2
	buf(14)= itrace
	k=15
	buf(k:k+9) = ispg_idx(1:10); k=k+10
	buf(k:k+9) = ispl_idx(1:10); k=k+10
      endif   !   (Master)
c
       call MPI_BCAST(buf,Nbuf,MPI_INTEGER,0,MPI_COMM_WORLD,Ierr)     

c
        if (HWorker.or.VWorker) then
	iend  = buf(1)
! /aqspid/ in aqindx.cmm stuff
	iair  = buf(2)
	ih2o  = buf(3)
	io2   = buf(4)
	ico   = buf(5)
	ino2  = buf(6)
	iho2  = buf(7) 
	iso2  = buf(8)
	io3   = buf(9)
	ich4  = buf(10)
	ico2  = buf(11)
	ih2   = buf(12)
	in2   = buf(13)
	itrace= buf(14)
	k=15
	ispg_idx(1:10) = buf(k:k+9); k=k+10
	ispl_idx(1:10) = buf(k:k+9); k=k+10
      endif  !    (Worker)
c       
      end  subroutine int_distrib1
      



c-----------------------------------------------------------------------------------------------------------------------------
      subroutine int_distrib2(numl,nbin,ixtrn,iytrn,iztrn,
     &                      irxng,irxnl,num,mdt,idate,iend) 
c -----------------------------
c  Master distributes yet another bunch of integers to Workers
c  A. Sandu, Jan 2001
c-----------------------------------------------------------------------------------------------------------------------------
c      
          implicit none
      include 'mpif.h'
c      
      integer ::  numl(3,4)
      integer :: nbin,ixtrn,iytrn,iztrn,irxng,irxnl,num,mdt
      integer :: idate(3), iend
      integer, parameter ::  Nbuf=24   
      integer ::  i, j, buf(Nbuf), Ierr, rc, status(MPI_STATUS_SIZE)    
c
      if (Master) then
        do j=1,4
	  do i=1,3
            buf((j-1)*3+i) = numl(i,j)
	  enddo
	enddo  
	do i=1,3
            buf(12+i) = idate(i)
	enddo
	buf(16) = nbin
	buf(17) = ixtrn
	buf(18) = iytrn
	buf(19) = iztrn
	buf(20) = irxng
	buf(21) = irxnl  
	buf(22) = num  
	buf(23) = mdt  
	buf(24) = iend
      endif   !   (Master)
c
      call MPI_BCAST(buf,Nbuf,MPI_INTEGER,0,MPI_COMM_WORLD,Ierr)
c      
      if (HWorker.or.VWorker) then
        do j=1,4
	  do i=1,3
           numl(i,j) = buf((j-1)*3+i)
	  enddo
	enddo  
	do i=1,3
           idate(i) = buf(12+i)
	enddo
	nbin  = buf(16)
	ixtrn = buf(17)
	iytrn = buf(18)
	iztrn = buf(19)
	irxng = buf(20)
	irxnl = buf(21)  
	num   = buf(22)  
	mdt   = buf(23)  
	iend  = buf(24)
      endif  !    (Worker)
c       
      end  subroutine int_distrib2
      
c-----------------------------------------------------------------------------------------------------------------------------
      subroutine real_distrib(ix,iy,iz,is,dx,dy,
     &                 sigmaz,dht,baseh,rmw,dt,ut) 
c -----------------------------
c  Master distributes a bunch of reals to Workers
c  A. Sandu, Jan 2001
c-----------------------------------------------------------------------------------------------------------------------------
      implicit none
      include 'mpif.h'
c      
      integer :: ix,iy,iz,is,k
      real   :: dx(ix), dy(iy), sigmaz(iz)
      real   :: dht, baseh, rmw(is), dt, ut
      integer :: Nbuf, Ierr
      real   :: buf(ix+iy+iz+is+4)
c
      Nbuf = ix+iy+iz+is+4
      if (Master) then 
        k=0
	buf(k+1:k+ix)=dx(1:ix) ; k=k+ix
	buf(k+1:k+iy)=dy(1:iy) ; k=k+iy
	buf(k+1:k+iz)=sigmaz(1:iz) ; k=k+iz
	buf(k+1)=dht ; k=k+1
	buf(k+1)=baseh ; k=k+1	
	buf(k+1:k+is)=rmw(1:is) ; k=k+is
	buf(k+1)=dt ; k=k+1
	buf(k+1)=ut ; k=k+1
	if (k.ne.Nbuf) then
	  print*, 'Error in real_distrib. Nbuf=',Nbuf,
     &             '   needed=',k
          stop   	  
	endif	
      endif     
c
      call MPI_BCAST(buf(1),Nbuf,MPI_REAL,0,MPI_COMM_WORLD,Ierr)
c      
      if (HWorker.or.VWorker) then
        k=0
	dx(1:ix)=buf(1:ix) ; k=ix
	dy(1:iy)=buf(k+1:k+iy) ; k=k+iy
	sigmaz(1:iz)=buf(k+1:k+iz) ; k=k+iz
	dht=buf(k+1) ; k=k+1
	baseh=buf(k+1) ; k=k+1	
	rmw(1:is)=buf(k+1:k+is); k=k+is
	dt=buf(k+1); k=k+1
	ut=buf(k+1); k=k+1
      endif     
c       
      end  subroutine real_distrib


c-----------------------------------------------------------------------------------------------------------------------------
      subroutine real_distrib2(ix,iy,iz,is,dx,dy,
     &                        sigmaz,dht,baseh,rmw,dt,ut) 
c -----------------------------
c  Master distributes a bunch of reals to Workers
c   (non-broadcast version of real_distrib)
c  A. Sandu, Jan 2001
c-----------------------------------------------------------------------------------------------------------------------------
      implicit none
      include 'mpif.h'
c      
      integer :: ix,iy,iz,is,k
      real   :: dx(ix), dy(iy), sigmaz(iz)
      real   :: dht, baseh, rmw(is), dt, ut
      integer :: Nbuf, Ierr, status(MPI_STATUS_SIZE) 
      real   :: buf(ix+iy+iz+is+4)
      integer :: i
c
      Nbuf = ix+iy+iz+is+4
      if (Master) then 
        k=0
	buf(1:ix)=dx(1:ix) ; k=k+ix
	buf(k+1:k+iy)=dy(1:iy) ; k=k+iy
	buf(k+1:k+iz)=sigmaz(1:iz) ; k=k+iz
	buf(k+1)=dht ; k=k+1
	buf(k+1)=baseh ; k=k+1	
	buf(k+1:k+is)=rmw(1:is) ; k=k+is
	buf(k+1)=dt ; k=k+1
	buf(k+1)=ut ; k=k+1
	if (k.ne.Nbuf) then
	  print*, 'Error in real_distrib. Nbuf=',Nbuf,
     &             '   needed=',k
          stop   	  
	endif	
      endif     
c
      call MPI_BCAST(buf,Nbuf,MPI_REAL,0,MPI_COMM_WORLD,Ierr)
c      
      if (HWorker.or.VWorker) then
        k=0
	dx(1:ix)=buf(k+1:k+ix) ; k=k+ix
	dy(1:iy)=buf(k+1:k+iy) ; k=k+iy
	sigmaz(1:iz)=buf(k+1:k+iz) ; k=k+iz
	dht=buf(k+1) ; k=k+1
	baseh=buf(k+1) ; k=k+1	
	rmw(1:is)=buf(k+1:k+is); k=k+is
	dt=buf(k+1); k=k+1
	ut=buf(k+1); k=k+1
      endif     
c       
      end  subroutine real_distrib2

c-----------------------------------------------------------------------------------------------------------------------------
      subroutine cmm_distrib() 
c -----------------------------
c  Master distributes all common block data to Workers
c               (non-broadcast version)
c  A. Sandu, Jan 2001
c-----------------------------------------------------------------------------------------------------------------------------
      implicit none
      include 'mpif.h'
      include 'aqmax.param'
      !1.include 'aqcon1.cmm'
      real aqcon1,aqcon2 
      common /aqcon1/aqcon1(22*mxspg)      
      common /aqcon2/aqcon2(2000)
      integer iaqcon2
      common /iaqcon2/iaqcon2(22*mxspg+2135)
      !2.include 'aqcon2.cmm'
      double precision aqcon3
      common /aqcon3/aqcon3(3660)
      integer iaqcon3
      common /iaqcon3/iaqcon3(1862)
      !3.include 'aqcon5.cmm'
      double precision aqcon5
      common /aqcon5/aqcon5(13*mxspg+4*mxrxn+1)
      !4.include 'aqcont.cmm'
      integer aqcont, aq_unit
      common /aqcont/aqcont(24)
      common /aq_unit/aq_unit(6)      
      !6.include 'aqindx.cmm'
      integer iaqspid, iaqmatr
      common /aqspid/iaqspid(33)
      common /aqmatr/iaqmatr(mxspg*20+1)
      !include 'aqrxn.cmm'
      double precision aqrxng
      common /aqrxng/  aqrxng(66*mxrxn)
      real raqrxng
      common /raqrxng/  raqrxng(2*mxrxn)
      integer iaqrxng
      common /iaqrxng/ iaqrxng(64*mxrxn+1)
      double precision aqrxnl
      common /aqrxnl/  aqrxnl(63*mxrxn)
      integer iaqrxnl
      common /iaqrxnl/ iaqrxnl(64*mxrxn+1)
      real  aqphoto
      common /aqphoto/  aqphoto(3*30*500)
      integer iaqphoto
      common /iaqphoto/ iaqphoto(61)
      !include 'aqspec.cmm'
      character specinfo
      common /specinfo/specinfo(16*4*mxspg)
      integer ispecinfo
      integer, parameter :: mxal=5000
      common /ispecinfo/ispecinfo(4) 
      real  rspecinfo
      common /rspecinfo/rspecinfo(4*mxal)
      !include 'aqsymb.cmm'
      character aqsymb
      common /aqsymb/aqsymb( 16*(3*mxspg+1000) )
                 
c      
      integer :: Ierr, status(MPI_STATUS_SIZE) 
      integer :: Ndbuf, Nrbuf, Nibuf, Ncbuf 
      integer :: kd, kr, ki, kc, i
      double precision  :: dbuf(3660+13*mxspg+4*mxrxn+1+66*mxrxn
     &                      +63*mxrxn)
      real         :: rbuf(22*mxspg+2*mxrxn+2000
     &                      +4*mxal+3*30*500)
      integer       :: ibuf(22*mxspg+2135+1862+24+6
     &                      +33+mxspg*20+1+64*mxrxn
     &                      +1+64*mxrxn+1+61+4)
      character      :: cbuf(16*4*mxspg+16*(3*mxspg+1000))
c
      Ndbuf=+3660+13*mxspg+4*mxrxn+1+66*mxrxn+63*mxrxn 
      Nrbuf=+22*mxspg+2*mxrxn+2000+4*mxal+3*30*500
      Nibuf=+22*mxspg+2135+1862+24+6+33+mxspg*20+1
     &        +64*mxrxn+1+64*mxrxn+1+61+4 
      Ncbuf=+16*4*mxspg+16*(3*mxspg+1000) 
c
      if (Master) then 
          kd=0; kr=0; ki=0; kc=0;
c  aqcon1.cmm  
          rbuf(kr+1:kr+22*mxspg)=aqcon1(1:22*mxspg); 
	  kr=kr+22*mxspg     
          rbuf(kr+1:kr+2000)=aqcon2(1:2000)    ;     
	  kr=kr+2000
          ibuf(ki+1:ki+22*mxspg+2135)=
     &              iaqcon2(1:22*mxspg+2135); 
	         ki=ki+22*mxspg+2135
c  aqcon2.cmm  
          dbuf(kd+1:kd+3660)=aqcon3(1:3660);  kd=kd+3660
          ibuf(ki+1:ki+1862)=iaqcon3(1:1862); ki=ki+1862	  
c  aqcon5.cmm  
          dbuf(kd+1:kd+13*mxspg+4*mxrxn+1)=aqcon5(1:13*mxspg+4*mxrxn+1); 
	  kd=kd+13*mxspg+4*mxrxn+1
c  aqcont.cmm  
          ibuf(ki+1:ki+24)=aqcont(1:24); ki=ki+24
          ibuf(ki+1:ki+6)=aq_unit(1:6);  ki=ki+6
c  aqindx.cmm  
          ibuf(ki+1:ki+33)=iaqspid(1:33); ki=ki+33
          ibuf(ki+1:ki+mxspg*20+1)=iaqmatr(1:mxspg*20+1); 
	  ki=ki+mxspg*20+1
c  aqrxn.cmm
          dbuf(kd+1:kd+66*mxrxn)=aqrxng(1:66*mxrxn); 
	  kd=kd+66*mxrxn
          rbuf(kr+1:kr+2*mxrxn)=raqrxng(1:2*mxrxn); 
	  kr=kr+2*mxrxn
          ibuf(ki+1:ki+64*mxrxn+1)=iaqrxng(1:64*mxrxn+1); 
	  ki=ki+64*mxrxn+1
          dbuf(kd+1:kd+63*mxrxn)=aqrxnl(1:63*mxrxn); 
	  kd=kd+63*mxrxn
          ibuf(ki+1:ki+64*mxrxn+1)=iaqrxnl(1:64*mxrxn+1); 
	  ki=ki+64*mxrxn+1
          rbuf(kr+1:kr+3*30*500)=aqphoto(1:3*30*500); 
	  kr=kr+3*30*500
          ibuf(ki+1:ki+61)=iaqphoto(1:61); ki=ki+61
c   aqspec.cmm 
          cbuf(kc+1:kc+16*4*mxspg)=specinfo(1:16*4*mxspg); 
	  kc=kc+16*4*mxspg
          ibuf(ki+1:ki+4)=ispecinfo(1:4); ki=ki+4 
          rbuf(kr+1:kr+4*mxal)=rspecinfo(1:4*mxal); kr=kr+4*mxal
c   qsymb.cmm
          cbuf(kc+1:kc+16*(3*mxspg+1000))=
     &            aqsymb(1:16*(3*mxspg+1000));
	        kc=kc+16*(3*mxspg+1000)	  
c
	if (kd.ne.Ndbuf) then
	  print*, 'Error in double cmm_distrib. Nbuf=',Ndbuf,
     &             '   needed=',kd
          stop   	  
	else if (kr.ne.Nrbuf) then
	  print*, 'Error in real cmm_distrib. Nbuf=',Nrbuf,
     &             '   needed=',kr
          stop   	  
	else if (ki.ne.Nibuf) then
	  print*, 'Error in integer cmm_distrib. Nbuf=',Nibuf,
     &             '   needed=',ki
          stop   	  
	else if (kc.ne.Ncbuf) then
	  print*, 'Error in character cmm_distrib. Nbuf=',Ncbuf,
     &             '   needed=',kc
          stop   	  
	endif
		
      endif     
c

      call MPI_BCAST(dbuf,Ndbuf,MPI_DOUBLE_PRECISION,0,
     &                  MPI_COMM_WORLD,Ierr)
      call MPI_BCAST(rbuf,Nrbuf,MPI_REAL,0,MPI_COMM_WORLD,Ierr)
      call MPI_BCAST(ibuf,Nibuf,MPI_INTEGER,0,MPI_COMM_WORLD,Ierr)
      call MPI_BCAST(cbuf,Ncbuf,MPI_CHARACTER,0,MPI_COMM_WORLD,Ierr)

c      
      if (HWorker.or.VWorker) then
c
          kd=0; kr=0; ki=0; kc=0;
c  aqcon1.cmm  
          aqcon1(1:22*mxspg)=rbuf(kr+1:kr+22*mxspg); kr=kr+22*mxspg     
          aqcon2(1:2000)=rbuf(kr+1:kr+2000)    ;     kr=kr+2000
          iaqcon2(1:22*mxspg+2135)=ibuf(ki+1:ki+22*mxspg+2135); 
	                                  ki=ki+22*mxspg+2135
c  aqcon2.cmm  
          aqcon3(1:3660)=dbuf(kd+1:kd+3660);  kd=kd+3660
          iaqcon3(1:1862)=ibuf(ki+1:ki+1862); ki=ki+1862	  
c  aqcon5.cmm  
          aqcon5(1:13*mxspg+4*mxrxn+1)=dbuf(kd+1:kd+13*mxspg+4*mxrxn+1); 
	  kd=kd+13*mxspg+4*mxrxn+1
c  aqcont.cmm  
          aqcont(1:24)=ibuf(ki+1:ki+24); ki=ki+24
          aq_unit(1:6)=ibuf(ki+1:ki+6);  ki=ki+6
c  aqindx.cmm  
          iaqspid(1:33)=ibuf(ki+1:ki+33); ki=ki+33
          iaqmatr(1:mxspg*20+1)=ibuf(ki+1:ki+mxspg*20+1); 
	  ki=ki+mxspg*20+1
c  aqrxn.cmm
          aqrxng(1:66*mxrxn)=dbuf(kd+1:kd+66*mxrxn); 
	  kd=kd+66*mxrxn
          raqrxng(1:2*mxrxn)=rbuf(kr+1:kr+2*mxrxn); 
	  kr=kr+2*mxrxn
          iaqrxng(1:64*mxrxn+1)=ibuf(ki+1:ki+64*mxrxn+1); 
	  ki=ki+64*mxrxn+1
          aqrxnl(1:63*mxrxn)=dbuf(kd+1:kd+63*mxrxn); 
	  kd=kd+63*mxrxn
          iaqrxnl(1:64*mxrxn+1)=ibuf(ki+1:ki+64*mxrxn+1); 
	  ki=ki+64*mxrxn+1
          aqphoto(1:3*30*500)=rbuf(kr+1:kr+3*30*500); 
	  kr=kr+3*30*500
          iaqphoto(1:61)=ibuf(ki+1:ki+61); ki=ki+61
c   aqspec.cmm 
          specinfo(1:16*4*mxspg)=cbuf(kc+1:kc+16*4*mxspg); 
	  kc=kc+16*4*mxspg
          ispecinfo(1:4)=ibuf(ki+1:ki+4); ki=ki+4 
          rspecinfo(1:4*mxal)=rbuf(kr+1:kr+4*mxal); 
	  kr=kr+4*mxal
c   qsymb.cmm
          aqsymb(1:16*(3*mxspg+1000))=
     &           cbuf(kc+1:kc+16*(3*mxspg+1000));
	       kc=kc+16*(3*mxspg+1000)	  
c
      endif !  MyId.ne.0   
c           
      end  subroutine cmm_distrib

      end module HVStemCommunication
