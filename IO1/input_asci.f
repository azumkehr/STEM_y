c**********************************************************************
      subroutine input1(istday2,ut,iperiod2,ix,iy,iz,sg1,sl1,sp1,
     1    tlon,tlat,h,hdz,dx,dy,dz,dht,sigmaz,baseh,dt2,iunit,iout)
c****************************************************************************
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-----
      parameter (mxspg=250,n_var2=3)
      include 'aqcon1.cmm'
      include 'aqsymb.cmm'
      include 'aqcont.cmm'
cc
      dimension sg1(ix,iy,iz,1),sl1(ix,iy,iz,1)
      dimension sp1(ix,iy,iz,1)
      dimension tlon(ix,iy),tlat(ix,iy),h(ix,iy)
      dimension hdz(ix,iy,1),dz(ix,iy,1),sigmaz(1),istday(3)
      dimension iunit(25),iout(25),idum1(3),dx(1),dy(1)
      dimension istday2(3)
      dimension var_2d(n_var2)
      character*16 chdum,var_2d
      data var_2d/'LONGITUDE','LATITUDE','TOPOGRAPHY'/
cc
      namelist /aqtime/istday,isthr,iperiod,iprnt,dt
c------------------------------------------------------------------------
c     read namelist input file
      open(11,file='aqms.mif',status='old')
      read(11,aqtime)
      do i=1,3
        istday2(i)=istday(i)
      enddo
      iperiod2=iperiod
      dt2=dt 
      close(11)
      ut=isthr
c     call datec
c     call headline
c-------------------------------------------------------------------------
c     open all the data file
      call aq_open(iunit,iout)
c-------------------------------------------------------------------------
c        read geological data: longitude, latitude, terrain base ht
      if(iunit(8).ne.0) then
	 do 8 iter=1,n_var2
         call aq_readhd(iunit(8),1)
         read(iunit(8),'(a)',end=8) chdum
         call aq_chop(16,chdum)
         call aq_find(n_var2,chdum,var_2d,lpsec,iflag)
         if(iflag.ne.0) then
            write(6,*) '** critical error in input2 **',chdum
            stop
         else
            if(lpsec.eq.1) call aq_read(iunit(8),ix*iy,tlon)
            if(lpsec.eq.2) call aq_read(iunit(8),ix*iy,tlat)
            if(lpsec.eq.3) call aq_read(iunit(8),ix*iy,h)
         endif
8        continue
	 if(iz.gt.2) then
	    do 10 i=1,ix
	    do 10 j=1,iy
	    do k=1,iz
	    hdz(i,j,k)=sigmaz(k)*(dht-h(i,j))+baseh+h(i,j)
	    enddo
	    do 10 k=2,iz
	    dz(i,j,k-1)=hdz(i,j,k)-hdz(i,j,k-1)
10          continue
         endif
      else 
	 write(6,*) '** critical error: no lat-lon file **'
	 stop
      endif
      call input_locate(ix,iy,iz,iunit,sg1)
c-------------------------------------------------------------------------
c                                 initial condition
      call aq_readhd(iunit(3),1)
      call aq_read_conc(iunit(3),ix*iy*iz,numl(1,3),numl(1,3),
     1                 sname(1,1),sg1,
     2 '** error in reading gas initial concentrations **')
c---------------------------------------------------------------------
      if(iunit(15).ne.0) then
      do 20 i=1,ix
      do 20 j=1,iy
      do 20 k=1,iz
      do 20 ibin=1,nbin
c     do l=1,numl(3,1)
      sp1(i,j,k,ibin)=10.**(-7+ibin)
      if(sp1(i,j,k,ibin).gt.10.0) sp1(i,j,k,ibin)=10.0
20    continue
c     call aq_readhd(iunit(15),1)
c     call aq_read_conc(iunit(15),ix*iy*iz,
c    1           numl(3,3),numl(3,3),sname(3,1),sp1,
c    2 '** error in reading particle initial concentrations **')
      endif
      return
      end
c**********************************************************************
      subroutine input2(ix,iy,iz,numsp,it,idate,sg1,u,v,w,kh,kv,t,
     1  wc,wr,rvel,q,em,vg,fz,sprc,sx,sy,sz,dx,dy,hdz,iunit)
c**********************************************************************
      parameter (mxspg=250)
      include 'aqcont.cmm'
      include 'aqindx.cmm'
      include 'aqsymb.cmm'
      include 'aqspec.cmm'
      parameter (n_var3=11)
      dimension u(ix,iy,1),v(ix,iy,1),w(ix,iy,1),kh(ix,iy,1)
      dimension kv(ix,iy,1),t(ix,iy,1),wc(ix,iy,1),wr(ix,iy,1)
      dimension sg1(ix,iy,iz,1),q(ix,iy,1),em(ix,iy,iz,1)
      dimension rvel(ix,iy,1),idate(1)
      dimension vg(ix,iy,1),fz(ix,iy,1),sprc(1),dx(1),dy(1),hdz(1)
      dimension sx(iy,iz,2,1),sy(ix,iz,2,1),sz(ix,iy,1)
      dimension iunit(25)
      character*16 chdum,chvar3(n_var3)
      real kh, kv
      data chvar3/'U','V','W','KH','KV','T','P','WC','WR','RVEL','H2O'/
      kount=it-1
c---------------------------------------------------------------------
c      read  2-d field
      if(mod(kount,ifreq2(1)).eq.0.and.iunit(9).ne.0) then
         read(iunit(9),'(a)') chdum
         call aq_read(iunit(9),ix*iy,sprc)
      endif
c------------------------------------------------------------------
c      read 3-d field
      numt=ix*iy*iz
      if(mod(kount,ifreq3(1)).eq.0.and.iunit(10).ne.0) then
         do iter=1,nvar3(1)
         read(iunit(10),'(a)') chdum
         call aq_chop(16,chdum)
         call aq_find(n_var3,chdum,chvar3,lpsec,iflag)
         if(iflag.ne.0) then
            write(6,*) '** critical error in input2 **',chdum
            stop
         else
            if(lpsec.eq.1) call aq_read(iunit(10),numt,u)
            if(lpsec.eq.2) call aq_read(iunit(10),numt,v)
            if(lpsec.eq.3) call aq_read(iunit(10),numt,w)
            if(lpsec.eq.4) call aq_read(iunit(10),numt,kh)
            if(lpsec.eq.5) call aq_read(iunit(10),numt,kv)
            if(lpsec.eq.6) call aq_read(iunit(10),numt,t)
            if(lpsec.eq.7) 
     1          call aq_read(iunit(10),numt,sg1(1,1,1,iair)) 
            if(lpsec.eq.8) call aq_read(iunit(10),numt,wc)
            if(lpsec.eq.9) call aq_read(iunit(10),numt,wr)
            if(lpsec.eq.10)call aq_read(iunit(10),numt,rvel)
            if(lpsec.eq.11) 
     1          call aq_read(iunit(10),numt,sg1(1,1,1,ih2o))
         endif
         enddo
         do 10 i=1,ix
         do 10 j=1,iy
         do 10 k=1,iz
         sg1(i,j,k,iair)=2.687e+19*sg1(i,j,k,iair)
     1                  /1013.*273/t(i,j,k)   
		 ! 2.687e+19=avo/22400=6.02e+23/22400 molecules/cm3
         sg1(i,j,k,ih2o)=sg1(i,j,k,ih2o)*22.84/18.0
     1                  *sg1(i,j,k,iair)
10       continue
      endif
c---------------------------------------------------------------------
c     Read surface emission
      if(mod(kount,iqm).eq.0.and.iunit(5).ne.0.and.iz.gt.1) then
	 write(6,*) '** reading emission data**'
         call aq_readhd(iunit(5),1)
         call aq_read_conc(iunit(5),ix*iy,emnum,nqm,emname,q,
     2                 '** error in reading surface emissions **')
         call aq_speci(ix*iy,numsp,emnum,emal,q)
         call aq_mult_val(ix*iy*numsp,q,0.01)
      endif
c---------------------------------------------------------------------
c     Read Dry Deposition Velocities
      if(mod(kount,ivg).eq.0.and.iunit(6).ne.0.and.iz.gt.1) then
         call aq_readhd(iunit(6),1)
         call aq_read_conc(iunit(6),ix*iy,vgnum,nvg,vgname,vg,
     2                 '** error in reading Dry Dep. velocity **')
         call aq_speci(ix*iy,numsp,vgnum,vgal,vg)
      endif
c---------------------------------------------------------------------
c     Read Boundar Conditions
      if(mod(kount,ibd).eq.0.and.iunit(7).ne.0) then
         call aq_readhd(iunit(7),1)
         call aq_read_bc(iunit(7),ix,iy,iz,numsp,nbd,sname(1,1),
     1     sx,sy,sz,'** error in reading boundary conditions **')
      endif
      do i=1,ix
      do j=1,iy
      sz(i,j,7)=30.e-09*sg1(i,j,iz,iair)
      enddo
      enddo
      return
      end
c**************************************************************************
      subroutine input_locate(ix,iy,iz,iunit,sg1)
c**************************************************************************
      dimension iunit(25),sg1(1),istday(3)
      include 'aqcont.cmm'
      namelist /aqtime/istday,isthr,iperiod,iprnt,dt
      open(11,file='aqms.mif',status='old')
      read(11,aqtime)
      close(11)
c-----------------------------------------------------------------------
c     Initialize time dependent input files
      call aq_locat_time(ix*iy*iz,iunit(4),istday,isthr,iem,nem,sg1,
     1                   'in elevated emission file' )
      call aq_locat_time(ix*iy*iz,iunit(5),istday,isthr,iqm,nqm,sg1,
     1                   'in surface emission file'  )
      call aq_locat_time(ix*iy*iz,iunit(6),istday,isthr,ivg,nvg,sg1,
     1                   'in dry deposition file'  )
      call aq_locat_time(ix*iy*iz,iunit(7),istday,isthr,ibd,nbd,sg1,
     1                   'in boundary file'  )
      call aq_locat_time(ix*iy*iz,iunit(9),istday,isthr,
     1     ifreq2(1),nvar2(1),sg1,'in two dimension file ' )
      call aq_locat_time(ix*iy*iz,iunit(10),istday,isthr,
     1     ifreq3(1),nvar3(1),sg1,'in three dimension file' )
c----------------------------------------------------------------------
      return
      end
c***********************************************************************
      subroutine aq_open(iunit,iout)
c***********************************************************************
      implicit none
      dimension iunit(25),iout(25) 
      include 'aqfile.cmm'
      character*80 old,emarea
      integer iunit,iout,i
c-----------------------------------------------------------------------
c     initialize file name
      call aq_blank( 80,specif)
      call aq_blank( 80,gmechf)
      call aq_blank( 80,blank)
      call aq_blank( 80,initf)
      call aq_blank( 80,emif)
      call aq_blank( 80,emarea)
      call aq_blank( 80,depvel)
      call aq_blank( 80,bdf)
      call aq_blank( 80,domain)
      call aq_blank( 80,meteo2d)
      call aq_blank( 80,meteo3d)
      call aq_blank( 80,init_lqf)
      call aq_blank( 80,init_ptf)
      call aq_blank( 80,emarea_pt)
      call aq_blank( 80,lmechf)
c------------------------------------------------------------------------
c     read namelist input file
      open(11,file='aqms.mif',status='old')
      read(11,aqinpf)
      close(11)
cc
      do i=1,22
      iunit(i)=i+19
      enddo
      old='old'
      call aq_openfi(iunit(1),specif,old)
      call aq_openfi(iunit(2),gmechf,old)
      call aq_openfi(iunit(3),initf,old)
      call aq_openfi(iunit(4),emif,old)
      call aq_openfi(iunit(5),emarea,old)
      call aq_openfi(iunit(6),depvel,old)
      call aq_openfi(iunit(7),bdf,old)
      call aq_openfi(iunit(8),domain,old)
      call aq_openfi(iunit(9),meteo2d,old)
      call aq_openfi(iunit(10),meteo3d,old)
      call aq_openfi(iunit(14),init_lqf,old)
      call aq_openfi(iunit(15),init_ptf,old)
      call aq_openfi(iunit(16),emarea_pt,old)
      call aq_openfi(iunit(21),lmechf,old)
      return
      end
