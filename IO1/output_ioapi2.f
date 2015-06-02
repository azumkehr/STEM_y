C********************************************************************
      subroutine prtini
c********************************************************************
      parameter (mxspg=250)
      include 'aqcon1.cmm'
      include 'aqsymb.cmm'
      do i=1,numl(1,3)
      write(*,'(2x,i3,2x,a20)') i,sname(i,1)
      enddo
      return
      end
c*********************************************************************
      subroutine prtemp(idate,ut,ix,iy,iz,sg1,sl1,iout)
c*********************************************************************
      parameter (mxspg=250,mxgr=200)
      include 'aqcon1.cmm'      
      include 'aqcont.cmm'      
      include 'aqsymb.cmm'
      include 'aqindx.cmm'
      dimension sg1(ix,iy,iz,1),sl1(ix,iy,iz,1),iout(1),idate(3)
      real dum1(mxgr*mxgr*7)
      call open_ioapw("AQRST")
      ihr=nint(ut)

      !Process H2SO4=H2SO4+SO41+SO42+SO41+SO43
      !Add BC1+BC2=> Ca1, OC1+OC2=>Ca2
      do i=1,ix
      do j=1,iy
      do k=1,iz
         sg1(i,j,k,11)=sg1(i,j,k,11)+sg1(i,j,k,104)+
     &   sg1(i,j,k,105)+sg1(i,j,k,106) !H2SO4=H2SO4+SO41+SO42+SO41+SO43
         sg1(i,j,k,72)=sg1(i,j,k,58)+sg1(i,j,k,59) ! BC1+BC2=> Ca1
         sg1(i,j,k,73)=sg1(i,j,k,60)+sg1(i,j,k,61) ! OC1+OC2=> Ca2
!         sg1(i,j,k,62)=sg1(i,j,k,62)+sg1(i,j,k,63)
!     &    +sg1(i,j,k,64) ! SS1-3 => SSF
      enddo
      enddo
      enddo

      do l=1,numl(1,2)
        call aq_copy_mat(ix*iy*iz,sg1(1,1,1,L),dum1)
        print *, '***************l,  sg1(10,10,10,l)  **********'
	print *,  l, sg1(10,10,10,l)
!        if(unit_out.ne.0) 	
!     1   call conv_conc(ix*iy*iz,1,dum1,sg1(1,1,1,iair),
!     1                  rmw,iair,unit_out)
        print *,  l, sg1(10,10,10,l)
      call write_out_ioapi("AQRST",idate(1),idate(2),idate(3),ihr,
     1      sname(l,1),dum1)
         print *,  l, sg1(10,10,10,l) 
       print *, '**************____________________***********'
      enddo
      lstat= close3("AQRST")
      return
      end
c**********************************************************************
      subroutine prtout(idate,ut,ix,iy,iz,sg1,sl1,iout)
c*********************************************************************
      parameter (mxspg=250,mxgr=200)
      include 'aqcon1.cmm'
      include 'aqsymb.cmm'
      include 'aqcont.cmm'
      include 'aqindx.cmm'
      integer aqspout,joutindex,aeroindex
      common /aqodx/aqspout(mxspg),noutsp,izout,joutindex(mxspg),
     1 joutsp,jzout,aeroindex(mxspg),nout_aero,izaero

      dimension sg1(ix,iy,iz,1),sl1(ix,iy,iz,1),iout(1),idate(3)
      dimension dum1(mxgr*mxgr*7)
      call open_ioapw("AQOUT")
      ihr=nint(ut)
      if(mxgr*mxgr*5.lt.ix*iy*izout) then
	 write(6,*) '** critical error in prtout **'
	 stop
      endif
     
      print *, 'idate,ut,ix,iy,iz,iout'
      print *, idate,ut,ix,iy,iz,iout 
      print *, 'mxgr,mxgr*mxgr*7, dum1(1001)'
      print *, mxgr,mxgr*mxgr*7, dum1(55009)
      print *, 'Before writing, Before writing'
       
      !Process H2SO4=H2SO4+SO41+SO42+SO41+SO43
      !Add BC1+BC2=> BC1, OC1+OC2=>OC1
      do i=1,ix
      do j=1,iy
      do k=1,iz
!         sg1(i,j,k,11)=sg1(i,j,k,11)+sg1(i,j,k,104)+
!     &   sg1(i,j,k,105)+sg1(i,j,k,106) !H2SO4=H2SO4+SO41+SO42+SO41+SO43
!         sg1(i,j,k,72)=sg1(i,j,k,58)+sg1(i,j,k,59) ! BC1+BC2=> BC1          
!         sg1(i,j,k,73)=sg1(i,j,k,60)+sg1(i,j,k,61) ! OC1+OC2=>OC1
      enddo
      enddo
      enddo      
 
      do L=1,55            ! gas-phase
      print *, l, sname(L,1)
      call aq_copy_mat(ix*iy*izout,sg1(1,1,1,L),dum1)

      call write_out_ioapi("AQOUT",idate(1),idate(2),idate(3),ihr,
     1    sname(L,1),dum1)
      enddo
     
      !BC
      call aq_copy_mat(ix*iy*izout,sg1(1,1,1,58),dum1)
      call write_out_ioapi("AQOUT",idate(1),idate(2),idate(3),ihr,
     1    sname(72,1),dum1)
      !OC
      call aq_copy_mat(ix*iy*izout,sg1(1,1,1,60),dum1)
      call write_out_ioapi("AQOUT",idate(1),idate(2),idate(3),ihr,
     1    sname(73,1),dum1)   
      !SSF
!      call aq_copy_mat(ix*iy*izout,sg1(1,1,1,62),dum1)
!      call write_out_ioapi("AQOUT",idate(1),idate(2),idate(3),ihr,
!     1    'SSF',dum1)
      !SSC
!      call aq_copy_mat(ix*iy*izout,sg1(1,1,1,65),dum1)
!      call write_out_ioapi("AQOUT",idate(1),idate(2),idate(3),ihr,
!     1    'SSC',dum1)     

      do L=66,70            ! OPM10,OPM25,DUST1-3
      lsp=L
      print *, lsp, sname(lsp,1)
      call aq_copy_mat(ix*iy*izout,sg1(1,1,1,lsp),dum1)

      call write_out_ioapi("AQOUT",idate(1),idate(2),idate(3),ihr,
     1    sname(lsp,1),dum1)
      enddo

 
      do L=112,135           ! gas-phase
      lsp=L
      print *, lsp, sname(lsp,1)
      call aq_copy_mat(ix*iy*izout,sg1(1,1,1,lsp),dum1)

      call write_out_ioapi("AQOUT",idate(1),idate(2),idate(3),ihr,
     1    sname(lsp,1),dum1)
      enddo



      lstat= close3("AQOUT")

      return
      end


c**********************************************************************
      subroutine prtchkp1(idate,ut,ix,iy,iz,nspec,sg1)
c
c     Print Checkpoint Level 1
c
c*********************************************************************
      parameter (mxspg=250,mxgr=200)
      include 'aqcon1.cmm'
      include 'aqsymb.cmm'
      include 'aqcont.cmm'
      include 'aqindx.cmm'
      integer aqspout,joutindex,aeroindex
      common /aqodx/aqspout(mxspg),noutsp,izout,joutindex(mxspg),
     1 joutsp,jzout,aeroindex(mxspg),nout_aero,izaero

      real sg1(ix,iy,iz,1),sl1(ix,iy,iz,1),iout(1),idate(3)
      real buffer(ix,iy,iz)
      call open_ioapw("AQCHKP1")
      ihr=nint(ut)

      do L=1,nspec             ! gas-phase
        print *, 'CHKP1:',L, sname(L,1)
        buffer(1:ix,1:iy,1:iz) = sg1(1:ix,1:iy,1:iz,L)
        call write_out_chkp("AQCHKP1",idate(1),idate(2),idate(3),ut,
     1    sname(L,1),buffer)
      end do
      lstat= close3("AQCHKP1")

      return
      end

            subroutine WRITE_OUT_CHKP
     +  (loutname, year, month, day, realhour, varname, field)
c------------------------------------------------------------------------------
c A. Sandu, 2002, after G. Calori, 1.11.99
c
c LAST REV:
c
c PURPOSE: Write a field to an AQMS output file
c          (EDSS/Models-3 I/O API format)
c
c PRECONDITIONS: 
c 'create_out_ioapi' must has been previously called
c
c INPUT:
c <loutname>            C  logical name of output file
c year, month, day      I  current year (yyyy), month, day
c hour                  I  current hour
c varname            C*16  name of variable to be written
c field(nwords)         R  field to be written
c 
c OUTPUT:
c on <loutname> file 
c
c CALLS: I/O API library
c------------------------------------------------------------------------------
      implicit none

c      include 'netcdf.inc'
      include 'PARMS3.EXT'      ! i/o API
      include 'FDESC3.EXT'      ! i/o API
      include 'IODECL3.EXT'     ! i/o API

      character*(*) loutname
      integer day, month, year, hour, minute, second
      real realhour
      integer timesec, restsec
      character*16 varname
      real field(1)

      integer julian, trimlen

      integer io_msg, ios
      integer date3d, time3d, jul3d
      character*256 outfile
      character*17 subname
      character*80 msg3d

      data io_msg /6/
      data subname /'write_out_chkp'/

c       convert to i/o api time

      timesec = nint(realhour*3600.0)
      hour    = timesec/3600 
      restsec = mod(timesec,3600)
      minute  = restsec/60
      second  = mod(restsec,60)

      jul3d = JULIAN (year, month, day)
      date3d = 1000 * year + jul3d    ! file current date YYYYDDD
      ! file current time HHMMSS
      time3d = 10000 * hour  + 100*minute + second

      if(time3d.eq.240000) then
        time3d = 0
        date3d = date3d + 1
      end if

c       write field

      if (.not.
     +    WRITE3(loutname(:TRIMLEN(loutname)),
     +           varname, date3d, time3d, field)
     +   ) then
        call ENVSTR (loutname(:TRIMLEN(loutname)),
     +               'AQMS output file', 'dummy_out.nc',
     +               outfile, ios)
        msg3d = 'Error writing var ' // varname
c     +     //  ' at date&time ', date3d, time3d,
c     +    ' on AQMS output file ' // outfile
        call M3EXIT (subname, date3d, time3d, msg3d, 1)
        stop
      endif

      end subroutine WRITE_OUT_CHKP


c**********************************************************************
      subroutine praero(idate,ut,ix,iy,iz,sg1,sl1,iout)
c*********************************************************************
      parameter (mxspg=250,mxgr=200)
      include 'aqcon1.cmm'
      include 'aqsymb.cmm'
      include 'aqcont.cmm'
      include 'aqindx.cmm'
      integer aqspout,joutindex,aeroindex
      character outaeronam(mxspg)*16
      common /aqodx/aqspout(mxspg),noutsp,izout,joutindex(mxspg),
     1 joutsp,jzout,aeroindex(mxspg),nout_aero,izaero,outaeronam

      dimension sg1(ix,iy,iz,1),sl1(ix,iy,iz,1),iout(1),idate(3)
      dimension dum1(mxgr*mxgr*7)

      call open_ioapw("AEROOUT")  ! open aerosol output file
      ihr=nint(ut)
      if(mxgr*mxgr*5.lt.ix*iy*izout) then
	 write(6,*) '** critical error in prtout **'
	 stop
      endif
      
      do L=1,nout_aero              !aerosol  
      lsp=aeroindex(L)
      print *, lsp, outaeronam(L)
      call aq_copy_mat(ix*iy*izaero,sg1(1,1,1,lsp),dum1)
      if(unit_aero.ne.0.and.lsp.le.numl(1,3))        ! AOE storage index is beyond numl(1,3)
     1   call conv_conc(ix*iy*izaero,1,dum1,sg1(1,1,1,iair),
     1                  rmw(lsp),iair,unit_aero)
      call write_out_ioapi("AEROOUT",idate(1),idate(2),idate(3),ihr,
     1    outaeronam(L),dum1)
      enddo
      lstat= close3("AEROOUT")
      return
      end


c**********************************************************************
      subroutine prjout(idate,ut,ix,iy,iz,sg1,sl1,iout)
c*********************************************************************
      include 'aqmax.param'
      include 'aqcon1.cmm'
      include 'aqsymb.cmm'
      include 'aqcont.cmm'
      include 'aqindx.cmm'
      include 'aqrxn.cmm'
      integer aqspout,joutindex,aeroindex
      common /aqodx/aqspout(mxspg),noutsp,izout,joutindex(mxspg),
     1 joutsp,jzout,aeroindex(mxspg),nout_aero,izaero

      dimension sg1(ix,iy,iz,1),sl1(ix,iy,iz,1),iout(1),idate(3)
      dimension dum1(mxgr*mxgr*7)
      call open_ioapw("JOUT")
      ihr=nint(ut)
      if(mxgr*mxgr*5.lt.ix*iy*izout) then
	 write(6,*) '** critical error in prtout **'
	 stop
      endif
      
      print*, joutsp,joutindex(1:joutsp)
      do l=1,joutsp
      lsp=joutindex(l)
      print *, lsp, jvname(lsp)
      call aq_copy_mat(ix*iy*jzout,sg1(1,1,1,lsp+jno2-1),dum1)
      call write_out_ioapi("JOUT",idate(1),idate(2),idate(3),ihr,
     1    jvname(lsp),dum1)
      enddo
      lstat= close3("JOUT")
      return
      end
c******************************************************************
      subroutine close_ioap
c******************************************************************
      lstat= close3("AQOUT")
      lstat= close3("AQRST")
      lstat= close3("JOUT")
      lstat= close3("AEROOUT")
      return
      end
