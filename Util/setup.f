c*****************************************************************************
      subroutine aq_setp(ix,iy,iz,numl,nbin,ixtrn,iytrn,iztrn,irxng,
     1  irxnl,ixm,iym,izm,ilm,dx,dy,dht,sigmaz,baseh,sg1,sl1,sp1,
     2  u,v,w,kh,kv,t,wc,wr,sprc,q,em,vg,fz,sx,sy,sz)
c****************************************************************************
      dimension numl(3,4),dx(1),dy(1),sigmaz(1),sg1(1),sl1(1),u(1),v(1)
      dimension w(1),kh(1),kv(1),t(1),wc(1),wr(1),sprc(1),em(1),vg(1)
      dimension fz(1),sx(1),sy(1),sz(1),sp1(1)
      real kv,kh
c-------------------------------------------------------------------------------------
      call aq_chem
      call aq_aero_int
      call aq_jobc
     1 (ix,iy,iz,numl,ixtrn,iytrn,iztrn,irxng,irxnl,ixm,iym,izm,ilm)
      call aq_grid(ix,iy,iz,dx,dy,dht,sigmaz,baseh)
      call aq_init(ix,iy,iz,numl,nbin,sg1,sl1,sp1,u,v,w,kh,kv,t,
     1     wc,wr,sprc,q,em,vg,fz,sx,sy,sz)
c
      call aq_info(numl(1,1))
c
      return
      end
c*****************************************************************************
      subroutine aq_chem
c*****************************************************************************
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-----
      parameter (mxrxn=500,mxspg=250)
      include 'aqcon1.cmm'
      include 'aqrxn.cmm'
      include 'aqsymb.cmm'
cc
      character*50 dname(20),old
      dname(1)='gas phase species'
      dname(2)='gas phase peroperty'
      dname(3)='gas phase reactions'
      dname(4)='photolysis(g)'
      dname(5)='liquid phase species'
      dname(6)='Equilibrium Constants'
      dname(7)='Liquid phase reactions'
      dname(8)='particle phase species'
c-------------------------------------------------------------------------
      old='old' 
      call aq_openfi(14,
     &  'aq_property.dat',
     &   old)
      call aq_locate(14,dname(1),iflag)
C    Longlive gas species, longlive + shortlive gas species,
c       longlive + shortlive + thirdbody gas species  
      read(14,*) numl(1,1),numl(1,2),numl(1,3) 
      do l=1,numl(1,3)
      read(14,'(a)') sname(l,1)  
      enddo
      rewind(14)
cc
      call aq_locate(14,dname(2),iflag)
      do l=1,numl(1,3)
      read(14,*) dgas(l),rmw(l),(hencof(l,j),j=1,2)
      enddo
      rewind(14)
cc
      call aq_locate(14,dname(3),iflag)
      read(14,*) nrxng
      do l=1,nrxng
      read(14,*) (aa(j,l),ee(j,l),bb(j,l),j=1,2),ff(l),
     1  fnn(l),(coerg(j,l),coesg(j,l),nus(j,l),nur(j,l),j=1,30),
     2  irs(l),irr(l),irt(l),irc(l)
      enddo
      rewind(14)
cc 
      call aq_locate(14,dname(4),iflag)
      read(14,*) npht
      do l=1,npht
      read(14,*) (ipht(l,j),j=1,2)
      read(14,*) (wls(l,j),absorb(l,j),qy(l,j),j=1,ipht(l,2))
      enddo
      rewind(14)
      call aq_specid_g
      nhion=0
c-------------------------------------------------------------------
      call aq_locate(14,dname(5),iflag)
      if(iflag.eq.0) then    ! liquid species
c--- longlive, longlive + shortlive, longlive + shortlive + thirdbody     
c---- longlive + shortlive + thirdbody + longlive counter
      read(14,*) numl(2,1),numl(2,2),numl(2,3),numl(2,4)
      do l=1,numl(2,3)
      read(14,'(a)') sname(l,2)  
      enddo
      do l=1,numl(2,4)
      read(14,'(a)') lname1(l)
      enddo
      read(14,*) (mgas(l),l=1,numl(2,3))
      read(14,*) (mliq(l),l=1,numl(1,3))
      read(14,*) ((l1to2(l,j),j=1,2),l=1,numl(2,4))
      read(14,*) ((l2to1(l,j),j=1,3),l=1,numl(2,3))
      read(14,*) (liqnum(l),l=1,numl(2,3))
      read(14,*) ((icharg(l,j),j=1,3),l=1,numl(2,3))
      rewind(14)
cc
      call aq_locate(14,dname(6),iflag)
      read(14,*) 
     1         (((eqkcof(l,j,k),l=1,numl(2,3)),j=1,3),k=1,3)
      rewind(14)
cc
      call aq_locate(14,dname(7),iflag)
      read(14,*) nrxnl
      do l=1,nrxnl
      read(14,*) aal(l),eel(l),bbl(l),
     1      (coerl(j,l),coesl(j,l),lus(j,l),lur(j,l),j=1,30),
     2      lrs(l),lrr(l),lrt(l),lrc(l)
      enddo
      call aq_specid_l
      endif
      rewind(14)
c------------------------------------------------------------------
      call aq_locate(14,dname(8),iflag)
      if(iflag.eq.0) then          ! solid particles
c----  longlive, longlive + shortlive, longlive + shortlive + thirdbody    
         read(14,*) numl(3,1),numl(3,2),numl(3,3)
         do l=1,numl(3,3)
         read(14,'(a)') sname(l,3)  
         enddo
      endif
c------------------------------------------------------------------
      close(14)
      return
      end
c***********************************************************************
      subroutine aq_aero_int
c**********************************************************************
      include 'aqmax.param'
      include 'aqsymb.cmm'
      include 'aqcon1.cmm'
      namelist /aqaero/nbin,daero,daero_min,daero_max
      nbin=0.
      daero_min=-10.
      daero_max=-10.
      open(14,file=
     & 'aqms.mif',
     &  status='old')
      read(14,NML=aqaero)
      close(14)
      if(nbin.eq.0) return
      if(numl(3,3).eq.0) then
	 write(6,*) '** critical error in species data **'
	 write(6,*) '** no particle species with nbin .ne. 0 **'
	 stop
      endif
      if(daero_min.gt.0.and.daero_max.gt.daero_min) then
	 dxpt=(alog10(daero_max)-alog10(daero_min))/float(nbin-1)
	 do iter=1,nbin
	 daero(iter)=alog10(daero_min)+float(iter-1)*dxpt
	 enddo
      endif
      return
      end
c**********************************************************************
      subroutine aq_specid_g
c**********************************************************************
      parameter(mxspg=250)
      include 'aqsymb.cmm'
      include 'aqcon1.cmm'
      include 'aqindx.cmm'
      dimension iadum(12)
      character*16 adumch(12)
      data adumch/'AIR','H2O','O2','CO','NO2','HO2','SO2','O3','CH4',
     1 'CO2','H2','N2'/
cc
      do i=1,12
      iadum(i)=-1
      enddo 
      do i=1,12 
      call aq_find(numl(1,3),adumch(i),sname(1,1),iadum(i),iflag)
      enddo
cc
      iair=iadum(1)
      ih2o=iadum(2)
      io2=iadum(3)
      ico=iadum(4)
      ino2=iadum(5)
      iho2=iadum(6)      
      iso2=iadum(7)
      io3=iadum(8)
      ich4=iadum(9)
      ico2=iadum(10)
      ih2=iadum(11)
      in2=iadum(12)
      itrace=numl(1,3)+2    ! gas tracer in sg1

      if(iair.le.0) then
	 write(6,*) '** critical error: no air in speci.dat **'
	 write(6,*) iadum(1),adumch(1)
	 stop
      endif
      if(ih2o.le.0) write(6,*) '**Warning: no water species '
c-------------------------------------------------------------------
c                     set parameters for Lurman Mechanism
      do i=1,2
      call aq_find(numl(1,3),adumch(i+4),sname(1,1),ispg_idx(i),iflag)
      if(iflag.ne.0) then
	 ispg_idx(i)=0
	 write(6,*) '** can not find **',adumch(i+4),'for lurman mech'
	 if(i.eq.1)  write(6,*) '** it is also required in TUV **'
      endif
      enddo
cc---------------------------------------------------------------------
c                    set parameters for TUV
      do i=3,4
      call aq_find(numl(1,3),adumch(i+4),sname(1,1),ispg_idx(i),iflag)
      if(iflag.ne.0) then
	 ispg_idx(i)=0
	 write(6,*) '** can not find **',adumch(i+4),'for TUV'
      endif
      enddo
c
      return
      end
c**********************************************************************
      subroutine aq_specid_l  
c**********************************************************************
      parameter(mxspg=250)
      include 'aqsymb.cmm'
      include 'aqcon1.cmm'
      include 'aqindx.cmm'
      dimension iadum(10)
      character*16 adumch(2)
      data adumch/'H+','OH-'/
cc
      do i=1,1
      iadum(i)=-1
      enddo 
      do i=1,1 
      call aq_find(numl(2,3),adumch(i),sname(1,2),iadum(i),iflag)
      enddo
cc
      nhion=iadum(1)
      nohion=iadum(2)
      if(nhion.le.0.) then
	 write(6,*) '** critical error: No H+ liquid phase species **'
	 write(6,*) 'nhion',nhion
	 stop
      endif
      call aq_find(numl(1,3),adumch(1),sname(1,1),iadum(1),iflag)
      if(iflag.ne.0) then
	 write(6,*) '** critical error: No H+ gas phase species **'
	 stop
      endif
      return
      end           
c*****************************************************************
      subroutine aq_jobc
C    1  (ix,iy,iz,numl,ixtrn,iytrn,iztrn,irxng,irxnl,ixm,iym,izm,ilm)
     1  (ix2,iy2,iz2,numl,ixtrn2,iytrn2,iztrn2,
     2   irxng2,irxnl2,ixm,iym,izm,ilm)
c*******************************************************************
      dimension numl(3,4)
      include 'aqcont.cmm'
      namelist /aqjob/ix,iy,iz,ixtrn,iytrn,iztrn,irxng,irxnl
      namelist /aqunit/unit_ini,unit_em,unit_qm,unit_bd,unit_out,
     1 unit_aero 
c-------------------------------------------------------------------------------------
      open(11,file=
     &   'aqms.mif',
     &    status='old')
      read(11,aqjob)
      ix2=ix
      iy2=iy
      iz2=iz
      ixtrn2=ixtrn
      iytrn2=iytrn
      iztrn2=iztrn
      irxng2=irxng
      irxnl2=irxnl
      close(11)
      unit_ini=0
      unit_em=0
      unit_qm=0
      unit_bd=0
      unit_out=0
      open(11,file=
     &   'aqms.mif',
     &    status='old')
      read(11,aqunit)
      close(11)
      if(ix.gt.ixm.or.iy.gt.iym.or.iz.gt.izm.or.numl(1,3).gt.iLm)then
	  write(6,*) '** critical error in number of grids or **'
	  write(6,*) '** number of chemical species **'
	  write(6,*) '** ix,iy,iz,nsp **',ix,iy,iz,numl(1,3)
	  write(6,*) '** ixm,iym,izm,iLm **',ixm,iym,izm,iLm
	  stop
      endif
      ierror=1 
      if(ixtrn.eq.0.and.ix.le.1) then 
         write(6,*) '** error in ixtrn or ix',ixtrn,ix 
      else if(iytrn.eq.0.and.iy.le.1) then
         write(6,*) '** error in iytrn or iy',iytrn,iy 
      else if(iztrn.eq.0.and.iz.le.1) then
         write(6,*) '** error in iztrn or iz',iytrn,iy 
      else
	 ierror=0
      endif
      if(ierror.ne.0) stop
      return
      end
c*******************************************************************
      subroutine aq_grid(ix,iy,iz,dx,dy,dht2,sigmaz,baseh2)
c****************************************************************
      include 'PARMS3.EXT'      ! i/o API
      include 'FDESC3.EXT'      ! i/o API
      include 'IODECL3.EXT'     ! i/o API 
      character*50 filex,filey,filez
      namelist /aqgrid/flagx,flagy,flagz,filex,filey,filez,
     1                 dxc,dyc,dht,baseh
CSandu      namelist /aqgrid/flagx,flagy,flagz,filex,filey,filez,
CSandu     1                 dxc,dyc,dht,baseh,dx,dy,sigmaz
      dimension dx2(1),dy2(1),sigmaz2(1)
      dimension dx(ix), dy(iy), sigmaz(iz)
c-----------------------------------------------------------------------------------------
      open(11,file=
     &   'aqms.mif',
     &    status='old')
      read(11,aqgrid)
      baseh2=baseh
      dht2=dht 
      close(11)
      if(ix.gt.2) then
	 if(flagx.ge.0.) then
	    do i=1,(ix-1)
	    dx(i)=dxc
	    enddo
	 call aq_test_val_r4((ix-1),dx,0.0,'**error:dx(i)=0',iflag)
         endif
      endif
cc
      if(iy.gt.2) then
	 if(flagy.ge.0.) then
	    do j=1,(iy-1)
	    dy(j)=dyc
	    enddo
	 call aq_test_val_r4((iy-1),dy,0.0,'**error:dy(i)=0',iflag)
         endif
      endif
cc
      if(iz.gt.2) then 
ctyh         if(flagz.gt.0.1) then
ctyh            dum1=float(iz-1)**flagz
ctyh            do k=1,iz
ctyh	       sigmaz(k)=float(k-1)**flagz/dum1
ctyh	    enddo
ctyh         endif
       if(.not.open3('METEO3D',FSREAD3,'aq_open')) then
        print*, 'failed to open METEO3D in aq_grid'
        stop
       endif

       if (.not. DESC3('METEO3D') ) then   ! get grid information from meteorological 3d file to fill the description of 3d chemical output
        print*, 'Error getting grid info from METEO3D in aq_grid' 
        stop
       endif
       dht = vgtop3d               ! meteorological model top
       do k=1,iz
        sigmaz(k)=vglvs3d(k)/vgtop3d
       enddo
       io_log=close3('METEO3D')
       call aq_test_val_r4(1,dht,0.0,'**error:dht=0',iflag)
       call aq_test_val_r4
     1          ((iz-1),sigmaz(2),0.0,'**error:z(i)=0',iflag)
      endif
      return
      end
c************************************************************************
      subroutine aq_hori_dt
     1 (ix,iy,iz,ixtrn,iytrn,u,v,dx,dy,dtmax,dt,mdt)
c***********************************************************************
      dimension u(1),v(1),dx(1),dy(1)
      if(ixtrn*iytrn.ne.0) then
         mdt=3600/int(dtmax+0.5)
      else
         cocr=0.4
         if(ixtrn.eq.0) cox=courant_num(ix,iy,iz,dx,u,dtmax,1)
         if(iytrn.eq.0) coy=courant_num(ix,iy,iz,dy,v,dtmax,2)
         comax=amax1(cox,coy)
         if(comax.le.cocr) then
            mdt=3600/int(dtmax+0.5)
         else
	    dt=dtmax*cocr/comax
	    mdt=int(3600./dt+0.99)
         endif
      endif
      dt=float(3600/mdt)
      return
      end
c************************************************************************
      subroutine aq_hori_dt_Strang
     1 (ix,iy,iz,ixtrn,iytrn,u,v,dx,dy,dtmax,dt,mdt)
c***********************************************************************
      dimension u(1),v(1),dx(1),dy(1)

      if(ixtrn*iytrn.ne.0) then
         mdt=3600/int(dtmax+0.5)
      else
         cocr=0.4
         if(ixtrn.eq.0) cox=courant_num(ix,iy,iz,dx,u,dtmax,1)
         if(iytrn.eq.0) coy=courant_num(ix,iy,iz,dy,v,dtmax,2)
         comax=amax1(cox,coy)
         if(comax.le.cocr) then
            mdt=3600/int(dtmax+0.5)
         else
	    dt=dtmax*cocr/comax
	    mdt=int(3600./dt+0.99)
         endif
      endif
      if(mod(mdt,2).eq.0) then
        mdt=mdt/2
      else
        mdt=(mdt+1)/2
      endif		
      dt=float(3600)/(2*mdt)
 
      return
      end
c*****************************************************************
      subroutine input_check(ix,iy,iz,numsp,sg1,
     1                       ixtrn,iytrn,iztrn,irxng )
c****************************************************************
      include 'aqmax.param'
      include 'aqcon1.cmm'
      include 'aqrxn.cmm'
      include 'aqindx.cmm'
      include 'aqcont.cmm'
      dimension sg1(ix*iy*iz,1)
      call aq_test_val_r4
     1          (ix*iy*iz,sg1(1,iair),0.0,'*error:air=0 **',iflag)
      if(iflag.eq.1) then
	 sair=6.02d+23/(22.4*1000.)   ! molecules/cm3
	 call aq_const_r4(ix*iy*iz,sg1(1,iair),sair)
	 write(6,*) '** air concentration is assumed to ',sg1(1,iair)
      endif
c--------------------------------------------------------------------------
      if(irxng.ne.0.and.nrxng.eq.0) then
	 write(6,*) '** critical error: no gas phase rxns **'
         stop
      endif
c--------------------------------------------------------------------------
      if(unit_ini.ne.0) 
     1  call conv_conc(ix*iy*iz,numsp,sg1,sg1(1,iair),rmw,iair,unit_ini)
      return
      end
c**************************************************************************
      subroutine aq_init(ix,iy,iz,numl,nbin,sg1,sl1,sp1,u,v,w,
     1     kh,kv,t,wc,wr,sprc,q,em,vg,fz,sx,sy,sz)
c**************************************************************************
      dimension sg1(1),sl1(1),u(1),v(1),w(1),kh(1),kv(1),t(1)
      dimension wc(1),wr(1),q(1),em(1),vg(1),fz(1),sprc(1)
      dimension sx(1),sy(1),sz(1),numl(3,4)
      common /aq_tuv2/nz,ini_tuv,nw0
      real kh,kv
      numt=ix*iy*iz
      call aq_zero_r4(numt*numl(1,3),sg1)
      if(nbin.ne.0) call aq_zero_r4(numt*numl(3,1)*nbin,sp1)
      call aq_zero_r4(numt,u)
      call aq_zero_r4(numt,v)
      call aq_zero_r4(numt,w)
      call aq_zero_r4(numt,kv)
      call aq_zero_r4(numt,kh)
      call aq_zero_r4(numt,t)
      call aq_zero_r4(numt,wc)
      call aq_zero_r4(numt,wr)
      call aq_zero_r4(ix*iy,sprc)
      call aq_zero_r4(ix*iy*numl(1,3),q)
      call aq_zero_r4(ix*iy*iz*numl(1,3),em)
      call aq_zero_r4(ix*iy*numl(1,3),vg)
      call aq_zero_r4(ix*iy*numl(1,3),fz)
      call aq_zero_r4(ix*iy*numl(1,3),sz)
      call aq_zero_r4(iy*iz*2*numl(1,3),sx)
      call aq_zero_r4(ix*iz*2*numl(1,3),sy)
      return
      end
c*************************************************************
      subroutine aq_speci_info
     1        (filen,numsp,sname,name1,numv,alpha,mxal)
c*************************************************************
      parameter (mxspg=250)
      dimension sname(1),char(mxspg),name1(1)
      dimension val(mxspg),alpha(numsp,1)
      character*16 name1,sname,char,blank,filen
cc
      call aq_blank(16,blank)
      numv=0
      if(filen.eq.blank) then
	 do l=1,numsp
	 name1(l)=sname(l)
	 enddo
	 numv=numsp
         return
      endif
      write(6,*) '***reading   ',filen
      open(14,file=filen,status='old')
      call aq_zero_r4(mxal,alpha)
      do 10 iter=1,mxspg
         call aq_readhd(14,1)
         read(14,*,end=98,err=99) name1(iter),num
         call aq_readhd(14,1)
         read(14,*,end=98,err=99) (char(i),val(i),i=1,num)
         do i=1,num
         call aq_find(numsp,char(i),sname,lpsec,iflag)
         if(iflag.eq.0) then
            alpha(lpsec,iter)=val(i)
         else
	    write(6,*) '** critical error in aq_speciate_info**'
	    write(6,*) 'name=',name1(iter),'  var=',char(i)
	    stop
         endif
         enddo
	 if(iter*numsp.gt.mxal) then
	    write(6,*) '** critical error: too many species **'
	    write(6,*) '** in **',filen
	    stop
         endif
10    continue
      write(6,*) '** too many data in aq_speciate_info **'
98    continue
      numv=iter-1
      close(14)
      return
99    write(6,*) '** error in the vicinity of ',name1(iter)
      stop
      end
c***********************************************************************
      subroutine aq_speci(npoints,numsp,numv,alpha,val)
c***********************************************************************
      dimension val(npoints,1),alpha(numsp,1)
      dimension speci(1000) ,lump(1000)
      real lump
      if(numv.eq.0) return
      do iter=1,npoints
         call aq_zero_r4(numsp,speci)
	 do j=1,numv             ! original species number
	 lump(j)=val(iter,j)
	 enddo
         do i=1,numsp
         do j=1,numv
         speci(i)=speci(i)+alpha(i,j)*lump(j)
         enddo
         enddo
	 do i=1,numsp         ! transport species
	 val(iter,i)=speci(i)
	 enddo
      enddo
      return
      end
c**********************************************************************
      subroutine aq_info(numsp)
c**********************************************************************
      parameter (mxspg=250)
      include 'aqsymb.cmm'
      include 'aqfile.cmm'
      include 'aqspec.cmm'
      call aq_blank( 16,spec_emf)
      call aq_blank( 16,spec_icf)
      call aq_blank( 16,spec_bcf)
      call aq_blank( 16,spec_vgf)
c------------------------------------------------------------------------
c     read namelist input file
      open(11,file=
     &   'aqms.mif',
     &     status='old')
      read(11,aqinpf)
      close(11)
cc
      call aq_speci_info
     1    (spec_emf,numsp,sname(1,1),emname,emnum,emal,mxal)
      call aq_speci_info
     1    (spec_icf,numsp,sname(1,1),icname,icnum,ical,mxal)
      call aq_speci_info
     1    (spec_bcf,numsp,sname(1,1),bcname,bcnum,bcal,mxal)
      call aq_speci_info
     1    (spec_vgf,numsp,sname(1,1),vgname,vgnum,vgal,mxal)
      return
      end
