c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-----
c*********************************************************************
      subroutine aq_readhd(iunit,iflag)
c***********************************************************************
      character*1 char(80)
      if(iflag.eq.1) then      ! formatted
	 do iter=1,1000
	 read(iunit,100,end=8)  (char(i),i=1,80)
	 if(char(1).eq.'$') then
	    write(6,200) iunit,(char(i),i=1,80)
         else if(char(1).eq.'#') then
	 else
	    backspace iunit
	    return
         endif
	 enddo
      else
	 do iter=1,1000
	 read(iunit,end=8)  (char(i),i=1,80)
	 if(char(1).eq.'$') then
	    write(6,200) iunit,(char(i),i=1,80)
         else if(char(1).eq.'#') then
	 else
	    backspace iunit
	    return
         endif
	 enddo
      endif
8     continue
100   format(80a1)
200   format(2x,'iunit=',i3,2x,80a1)
      return
      end
c**********************************************************************
      subroutine aq_openfi(iunit,fname,type)
c**********************************************************************
      character*(*) fname,type
      character*50 blank
      call aq_blank(50,blank)
      if(fname.eq.blank) then
	 iunit=0
      else
         open(iunit,file=fname,status=type)
      endif
      return
      end
c**********************************************************************
      subroutine aq_locate(iunit,char,iflag)
c***********************************************************************
      character*(*) char
      character*80 dum1
      nchar=len(char)
      iflag=0
      do iter=1,10000
      read(iunit,'(a)',end=98) dum1(1:nchar)
      if(dum1(1:nchar).eq.char) return
      enddo
98    iflag=1
      return
      end
c**********************************************************************
      subroutine aq_locat_time
     1    (n,iunit,ist1,ist2,ifreq,nvar,rwk,cherr)
c**********************************************************************
      dimension rwk(1),ist1(3),irt1(3)
      character*(*) cherr
      if(iunit.eq.0) then 
	 ifreq=1
      else
         read(iunit,*) ifreq,irt1,irt2,nvar
         if(ist1(1).eq.irt1(1).and.ist1(2).eq.irt1(2).and.
     1      ist1(3).eq.irt1(3).and.ist2.eq.irt2) return
         write(6,*) '** inconsistent data file **'
	 write(6,*) ist1,ist2,irt1,irt2
	 write(6,*) cherr
	 stop
      endif
      return
      end
c**********************************************************************
      subroutine aq_find(num,cdum1,sname,lpsec,iflag)
c***********************************************************************
      dimension sname(1)
      character*(*) sname,cdum1
      iflag=0
      do 15 l=1,num
      if(cdum1.ne.sname(l)) go to 15
      lpsec=l
      return
15    continue
      iflag=1
      return
      end
c***********************************************************************
      subroutine aq_zero_i4(ntot,iy)
c***********************************************************************
      dimension iy(ntot)
      do i=1,ntot
      iy(i)=0.
      enddo
      return
      end
c***********************************************************************
      subroutine aq_zero_r4(ntot,y)
c***********************************************************************
      dimension y(ntot)
      do i=1,ntot
      y(i)=0.
      enddo
      return
      end
c********************************************************************
      subroutine aq_const_r8(ntot,y,val)
c*******************************************************************
      double precision y(1),val
      do i=1,ntot
      y(i)=val
      enddo
      return
      end
c********************************************************************
      subroutine aq_const_r4(n,val1,value)
c*******************************************************************
      dimension val1(1)
      do i=1,n
      val1(i)=value
      enddo
      return
      end
c********************************************************************
      subroutine aq_const_i4(ntot,y,val)
c*******************************************************************
      integer y(1),val
      do i=1,ntot
      y(i)=val
      enddo
      return
      end
c***********************************************************************
      subroutine aq_zero_r8(ntot,y)
c***********************************************************************
      double precision y
      dimension y(1)
      do i=1,ntot
      y(i)=0.
      enddo
      return
      end
c***********************************************************************
      subroutine aq_blank(ntot,y)
c***********************************************************************
      character*1 y(ntot)
      do i=1,ntot
      y(i)=' '
      enddo
      return
      end
c***********************************************************************
      subroutine aq_pos(ntot,y)
c***********************************************************************
      dimension y(ntot)
      do i=1,ntot
      if(y(i).lt.0.0) y(i)=0.0
      enddo
      return
      end
c***********************************************************************
      subroutine aq_read(nunit,ntot,y)
c***********************************************************************
      dimension y(ntot)
      call aq_readhd(nunit,1)
      if(nunit.gt.0) then
         read(nunit,*) y
      else
	 read(nunit) y
      endif
      return
      end
c***********************************************************************
      subroutine aq_write(nunit,ntot,y)
c***********************************************************************
      dimension y(ntot)
      if(nunit.gt.0) then
         write(nunit,'(1x,6e13.5)') y
      else
	 write(nunit) y
      endif
      return
      end
c***********************************************************************
      subroutine aq_phyval
c***********************************************************************
      common /aqsphy/avm,pi
      avm=6.02e+23
      pi=3.141592
      return
      end
C*****************************************************************
      SUBROUTINE aq_TOTAL2(ix,iy,il,CONC,TOT,DX,DY,I0,I1,J0,J1)
C*****************************************************************
      parameter (mxgr=400)
      DIMENSION CONC(IX,IY,IL),conc3(mxgr,mxgr),conc1(100)
      DIMENSION TOT(1),dx(1),dy(1)
cc
      DO 10 L=1,IL
cc
      do 12 i=1,ix
      do 12 j=1,iy
      conc3(i,j)=conc(i,j,L)
12    continue
cc
      do 15 j=1,iy
         do i=i0,i1
         conc1(i)=conc(i,j,l)
         enddo
         conc3(1,j)=sumint(i0,i1,dx,conc1)
15    continue
         do j=j0,j1
         conc1(j)=conc3(1,j)
         enddo
	 tot(l)=sumint(j0,j1,dy,conc1)
17    continue
10    CONTINUE
      RETURN
      END
C*****************************************************************
      SUBROUTINE aq_TOTAL3
     1         (ix,iy,iz,il,CONC,TOT,DX,DY,DZ,I0,I1,J0,J1,k0,k1)
C*****************************************************************
      parameter (mxgr=400)
      DIMENSION CONC(IX,IY,IZ,IL),conc3(mxgr,mxgr),conc1(mxgr)
      DIMENSION TOT(1),dx(1),dy(1),dz(1)
cc
      DO 10 L=1,IL
cc
      do 15 j=1,iy
      do 15 k=1,iz
         do i=i0,i1
         conc1(i)=conc(i,j,k,l)
         enddo
         conc3(j,k)=sumint(i0,i1,dx,conc1)
15    continue
      do 17 k=1,iz
         do j=j0,j1
         conc1(j)=conc3(j,k)
         enddo
	 conc3(1,k)=sumint(j0,j1,dy,conc1)
17    continue
      do k=k0,k1
      conc1(k)=conc3(1,k)
      enddo
      tot(l)=sumint(k0,k1,dz,conc1)
10    CONTINUE
      RETURN
      END
C**********************************************************************
      function sumint(i0,i1,dx,conc)
c************************************************************************
      dimension dx(1),conc(1)
      sumint=0.
      do i=i0,(i1-1)
      sumint=sumint+0.5*(conc(i)+conc(i+1))*dx(i)
      enddo
      return
      end
c ***************************************************************
c  subroutine trid is the solver for tri-diagonal matrix.
c
      subroutine trid(n,a,b,c,f,x)
c ***************************************************************
      dimension x(1),a(1),b(1),c(1),f(1),a1(200),c1(200),y(200)
      n1=n-1
      n2=n-2
      a1(1)=a(1)
      c1(1)=c(1)/a1(1)
      do 20 i=1,n2
      i1p=i+1
      a1(i1p)=a(i1p)-b(i1p)*c1(i)
   20 c1(i1p)=c(i1p)/a1(i1p)
      a1(n)=a(n)-b(n)*c1(n1)
cc
      y(1)=f(1)/a1(1)
      do 30 i=1,n1
      i1p=i+1
   30 y(i1p)=(f(i1p)-b(i1p)*y(i))/a1(i1p)
      x(n)=y(n)
      do 40 j=1,n1
   40 x(n-j)=y(n-j)-c1(n-j)*x(n+1-j)
      return
      end
c**************************************************************
c  subroutine decomp do a lu decompositin
c
      subroutine decomp(n,a,b,c)
c**************************************************************
      dimension a(1),b(1),c(1),a1(300),c1(300)
      n1=n-1
      n2=n-2
      a1(1)=a(1)
      c1(1)=c(1)/a1(1)
      do 20 i=1,n2
      i1p=i+1
      a1(i1p)=a(i1p)-b(i1p)*c1(i)
   20 c1(i1p)=c(i1p)/a1(i1p)
      a1(n)=a(n)-b(n)*c1(n1)
      do 10 i=1,n
      a(i)=a1(i)
      c(i)=c1(i)
10    continue
      return
      end
c ***************************************************************
c            subroutine trid1 is the solver
c     for tri-diagonal matrix which is already lu decomposed.
c
      subroutine trid1(n,a,b,c,f,x)
c ***************************************************************
      dimension x(1),a(1),b(1),c(1),f(1),y(300)
      n1=n-1
      n2=n-2
      y(1)=f(1)/a(1)
      do 30 i=1,n1
      i1p=i+1
   30 y(i1p)=(f(i1p)-b(i1p)*y(i))/a(i1p)
      x(n)=y(n)
      do 40 j=1,n1
   40 x(n-j)=y(n-j)-c(n-j)*x(n+1-j)
      return
      end
c*******************************************************************
      subroutine conv_conc2(nt,il,conc,sg4,wm,iair,iflag)
c********************************************************************
c     iflag=7:  convert molefraction to molecules/cm3
c     iflag=8:  convert molecules/cm3 to molefraction
c------------------------------------------------------------------
      data avm/6.02e+23/
      
      dimension conc(*),wm(*),sg4(*)
      print*,'conv_conc',nt,il,iair,sg4(1)
      if(iflag.eq.7) then
	 do 40 l=1,il
	   if(l.eq.iair) go to 40
           do i=1,nt
	     n=i+(l-1)*nt
             conc(n)=conc(n)*sg4(i)
	   enddo
40       continue
      else if(iflag.eq.8) then
	 do 45 l=1,il
	   if(l.eq.iair) go to 45
           do i=1,nt
	     n=i+(l-1)*nt
             conc(n)=conc(n)/sg4(i)
	   enddo
45       continue
      else
	 write(6,*) '** unknown flag in conc_convt: flag=',iflag
	 stop
      endif
      return
      end
c*******************************************************************
      subroutine conv_conc(nt,il,conc,sg4,wm,iair,iflag)
c********************************************************************
c     iflag=1:  convert moelcules/cm3 to ton/m3
c     iflag=2:  convert moelcules/cm3 to microgram/m3
cc
c     iflag=3:  convert ppb to molecules/cm3
c     iflag=4:  convert molecules/cm3 to ppb
cc
c     iflag=5:  convert m*molecules/cm3 to ton/m2
cc
c     iflag=7:  convert molefraction to molecules/cm3
c     iflag=8:  convert molecules/cm3 to molefraction
c------------------------------------------------------------------
       data avm/6.02e+23/
      dimension conc(1),wm(1),sg4(1)
      if(iflag.eq.1) then
	 do 10 l=1,il
	 if(l.eq.iair) go to 10
         do i=1,nt
	 n=i+(l-1)*nt
	 conc(n)=conc(n)/avm*wm(l)  !  unit=*gm/cm3=ton/m3  
	 enddo
10       continue
      else if(iflag.eq.2) then
	 do 15 l=1,il
	 if(l.eq.iair) go to 15
         do i=1,nt
	 n=i+(l-1)*nt
	 conc(n)=conc(n)/avm*wm(l)*1.E+12  !  unit=gm/cm3=ton/m3  
	 enddo
15       continue                              !  1.E+12 ton -> micrgram
      else if(iflag.eq.3) then
	 do 20 l=1,il
	 if(l.eq.iair) go to 20
         do i=1,nt
	 n=i+(l-1)*nt
         conc(n)=conc(n)*1.e-09*sg4(i)
	 enddo
20       continue
      else if(iflag.eq.4) then
	 do 25 l=1,il
	 if(l.eq.iair) go to 25
         do i=1,nt
	 n=i+(l-1)*nt
         conc(n)=conc(n)/sg4(i)*1.e+09
	 enddo
25       continue
      else if(iflag.eq.5) then
	 do 30 l=1,il
	 if(l.eq.iair) go to 30
         do i=1,nt
	 n=i+(l-1)*nt
	 conc(n)=conc(n)/avm*wm(l)     !  unit=m-gm/cm3= ton/m2  
	 enddo
30       continue
      else if(iflag.eq.7) then
	 do 40 l=1,il
	 if(l.eq.iair) go to 40
         do i=1,nt
	 n=i+(l-1)*nt
         conc(n)=conc(n)*sg4(i)
	 enddo
40       continue
      else if(iflag.eq.8) then
	 do 45 l=1,il
	 if(l.eq.iair) go to 45
         do i=1,nt
	 n=i+(l-1)*nt
         conc(n)=conc(n)/sg4(i)
	 enddo
45       continue
      else
	 write(6,*) '** unknown flag in conc_convt: flag=',iflag
	 stop
      endif
      return
      end
c********************************************************************
      subroutine copy_val(n,val1,val2)
c*******************************************************************
      dimension val1(1),val2(1)
      do i=1,n
      val2(i)=val1(i)
      enddo
      return
      end
c********************************************************************
      subroutine aq_copy_mat(n,val1,val2)
c*******************************************************************
      dimension val1(1),val2(1)
      do i=1,n
        val2(i)=val1(i)
      enddo
      return
      end
c********************************************************************
      subroutine aq_copy_mat_i4(n,val1,val2)
c*******************************************************************
      dimension val1(1),val2(1)
      integer val1,val2
      do i=1,n
      val2(i)=val1(i)
      enddo
      return
      end
c********************************************************************
      subroutine aq_copy_mat_r8(n,val1,val2)
c*******************************************************************
      dimension val1(1),val2(1)
      double precision val1,val2
      do i=1,n
      val2(i)=val1(i)
      enddo
      return
      end
c********************************************************************
      subroutine mult_val(n,val1,factor)
c*******************************************************************
      dimension val1(1)
      do i=1,n
      val1(i)=factor*val1(i)
      enddo
      return
      end
c********************************************************************
      subroutine aq_mult_val(n,val1,factor)
c*******************************************************************
      dimension val1(1)
      do i=1,n
      val1(i)=factor*val1(i)
      enddo
      return
      end
c*********************************************************************
      subroutine aq_add_val(n,val1,val2)
c*********************************************************************
      dimension val1(1)
      do i=1,n
      val1(i)=val1(i)+val2
      enddo
      return
      end
c********************************************************************
      subroutine aq_add_mat(n,val1,val2)
c*******************************************************************
      dimension val1(1),val2(1)
      do i=1,n
      val1(i)=val1(i)+val2(i)
      enddo
      return
      end
c********************************************************************
      subroutine aq_subs_mat(n,val1,val2)
c*******************************************************************
      dimension val1(1),val2(1)
      do i=1,n
      val1(i)=val1(i)-val2(i)
      enddo
      return
      end
c**************************************************************************
      subroutine d3tod1(ix,iy,iz,il,i1,i2,three,one,idir)
c**************************************************************************
      dimension three(ix,iy,iz,1),one(1)
      if(idir.eq.1) then
	 do 10 l=1,il
	 do 10 i=1,ix
	 num=i+(l-1)*ix
	 one(num)=three(i,i1,i2,l)
10       continue
         return
      else if(idir.eq.2) then
	 do 20 l=1,il
	 do 20 j=1,iy
	 num=j+(l-1)*iy
	 one(num)=three(i1,j,i2,l)
20       continue
         return
      else if(idir.eq.3) then
	 do 30 l=1,il
	 do 30 k=1,iz
	 num=k+(l-1)*iz
	 one(num)=three(i1,i2,k,l)
30       continue
         return
      else
	 write(6,*) '**** fatal error in d3tod1 idir='
	 stop
      endif
      return
      end
c*************************************************************************
      subroutine d3tod2(ix,iy,iz,il,i1,three,two,idir)
c*************************************************************************
      dimension three(ix,iy,iz,1),two(1)
      if(idir.eq.1) then
	 do 10 l=1,il
	 do 10 j=1,iy
	 do 10 k=1,iz
	 num=j+(k-1)*iy+(l-1)*iy*iz
	 two(num)=three(i1,j,k,l)
10       continue
         return
      else if(idir.eq.2) then
	 do 20 l=1,il
	 do 20 i=1,ix
	 do 20 k=1,iz
	 num=i+(k-1)*ix+(l-1)*ix*iz
	 two(num)=three(i,i1,k,l)
20       continue
         return
      else if(idir.eq.3) then
	 do 30 l=1,il
	 do 30 i=1,ix
	 do 30 j=1,iy
	 num=i+(j-1)*ix+(l-1)*ix*iy
	 two(num)=three(i,j,i1,l)
30       continue
         return
      else
	 write(6,*) '**** fatal error in d3tod2 idir='
	 stop
      endif
      return
      end
c**************************************************************************
      subroutine d1tod3(ix,iy,iz,il,i1,i2,three,one,idir)
c**************************************************************************
      dimension three(ix,iy,iz,1),one(1)
      if(idir.eq.1) then
	 do 10 l=1,il
	 do 10 i=1,ix
	 num=i+(l-1)*ix
	 three(i,i1,i2,l)=one(num)
10       continue
         return
      else if(idir.eq.2) then
	 do 20 l=1,il
	 do 20 j=1,iy
	 num=j+(l-1)*iy
	 three(i1,j,i2,l)=one(num)
20       continue
         return
      else if(idir.eq.3) then
	 do 30 l=1,il
	 do 30 k=1,iz
	 num=k+(l-1)*iz
         three(i1,i2,k,l)=one(num)
30       continue
         return
      else
	 write(6,*) '**** fatal error in d1tod3 idir='
	 stop
      endif
      return
      end
c**************************************************************************
      subroutine d2tod0(ix,iy,il,i1,i2,two,zero)
c**************************************************************************
      dimension two(ix,iy,1),zero(1)
      do l=1,il
      zero(l)=two(i1,i2,l)
      enddo
      return
      end
c**************************************************************************
      subroutine d3tod0(ix,iy,iz,il,i1,i2,i3,three,zero)
c**************************************************************************
      dimension three(ix,iy,iz,1),zero(1)
      do l=1,il
      zero(l)=three(i1,i2,i3,l)
      enddo
      return
      end
c**************************************************************************
      subroutine d0tod3(ix,iy,iz,il,i1,i2,i3,three,zero)
c**************************************************************************
      dimension three(ix,iy,iz,1),zero(1)
      do l=1,il
      three(i1,i2,i3,l)=zero(l)
      enddo
      return
      end
c**************************************************************************
      function courant_num(ix,iy,iz,dx,wind,dt,idir)
c**************************************************************************
      dimension wind(ix,iy,iz),dx(1) 
      courant_num=0.
      if(idir.eq.1) then
	 do 10 i=1,(ix-1)
	 do 10 j=1,iy
	 do 10 k=1,iz
	 dum1=abs((wind(i,j,k)+wind(i+1,j,k))*dt/(2.*dx(i)))
	 courant_num=amax1(dum1,courant_num)
10       continue
      else if(idir.eq.2) then
	 do 20 i=1,ix
	 do 20 j=1,(iy-1)
	 do 20 k=1,iz
	 dum1=abs((wind(i,j,k)+wind(i,j+1,k))*dt/(2.*dx(j)))
	 courant_num=amax1(dum1,courant_num)
20       continue
      else 
	 do 30 i=1,ix
	 do 30 j=1,iy
	 do 30 k=1,(iz-1)
	 dum1=abs((wind(i,j,k)+wind(i,j,k+1))*dt/(2.*dx(k)))
	 courant_num=amax1(dum1,courant_num)
30       continue
      endif
      return
      end
C*****************************************************************
      SUBROUTINE DATEC(IDATE,IYR,MO,IDY,IHR)
C*****************************************************************
      IYR=IDATE/1000000
      MO=(IDATE-IYR*1000000)/10000
      IDY=(IDATE-IYR*1000000-MO*10000)/100
      IHR=(IDATE-IYR*1000000-MO*10000-IDY*100)
      WRITE(6,*) 'IYR,MO,IDY,IHR',IYR,MO,IDY,IHR
      RETURN
      END
C*****************************************************************
      SUBROUTINE DATEC1(IDATE,IYR,MO,IDY)
C*****************************************************************
      IYR=IDATE/10000
      MO=(IDATE-IYR*10000)/100
      IDY=IDATE-IYR*10000-MO*100
      RETURN
      END
c*****************************************************************
      subroutine dater1(idate,iyr,mo,idy)
c*****************************************************************
      idate=iyr*10000+mo*100+idy
      return
      end
c******************************************************************
      function julday(iyr,mo,idy)
c******************************************************************
      dimension month(12) 
      data month/0,31,59,90,120,151,180,211,242,272,303,334/
      julday=month(mo)+idy
      if(mod(iyr,4).eq.0.and.mo.ge.3) julday=julday+1
      return
      end
c*******************************************************************
      subroutine aq_clock_j(jday,time)
c*******************************************************************
      do iter=1,100
      if(time.ge.24.0) then
	 time=time-24.0
         jday=jday+1
      else
	 return
      endif
      enddo
      return
      end
c*******************************************************************
      subroutine aq_clock_n(idate,time)
c*******************************************************************
      dimension month(12) ,idate(3)
      data month/31,28,31,30,31,30,31,31,30,31,30,31/
      if(mod(idate(1),4).eq.0) month(2)=29

      if(abs(nint(time)-time)<1e-4) time=real(nint(time)) ! remove time cut error
      do iter=1,100
         if(time.ge.24.0) then
	    time=time-24.0
            idate(3)=idate(3)+1
            if(idate(3).gt.month(idate(2))) then
	       idate(3)=idate(3)-month(idate(2))
	       idate(2)=idate(2)+1
	       if(idate(2).gt.12) then
                  idate(2)=idate(2)-12
	          idate(1)=idate(1)+1
	       endif
            endif
         else if (time<0.0) then
            time=time+24.0
            idate(3)=idate(3)-1
            if(idate(3)<=0.0) then
               idate(2)=idate(2)-1
               if(idate(2)<1) then
                  idate(2)=idate(2)+12
                  idate(1)=idate(1)-1
               endif
               idate(3) = month(idate(2))
            endif
         else
	    return
         endif
      enddo
      return
      end
c****************************************************************************
      subroutine aq_int1d_r8(ix,il,conc,num,xp,concp)
c****************************************************************************
      implicit double precision(a-h,o-z)
      dimension conc(ix,1),xp(1),concp(il,1)
      do iter=1,num
      ip=int(xp(iter))
      if(ip.lt.1.or.ip.gt.ix) then
	 write(6,*) '** critical error in inpt2d **'
	 write(6,*) '** out of bound: ip=',ip
	 return
      endif
      frx=xp(iter)-ip
      do l=1,il
      concp(l,iter)=conc(ip,l)+frx*(conc(ip+1,l)-conc(ip,l))
      enddo
      enddo
      return
      end
c****************************************************************************
      subroutine aq_int1v(numsp,num1,conc,xval,num2,xp,concp)
c****************************************************************************
      dimension conc(num1,1),xval(1),xp(1),concp(num2,1)
      do 10 iter=1,num2
      if(xp(iter).le.xval(1)) then
         do l=1,numsp
         concp(iter,l)=conc(1,l)
         enddo
	 go to 10
      else if(xp(iter).ge.xval(num1)) then
         do l=1,numsp
         concp(iter,l)=conc(num1,l)
         enddo
	 go to 10
      else
         do it1=1,(num1-1)
         if(xp(iter).ge.xval(it1).and.xp(iter).le.xval(it1+1))then
	    delta=xval(it1+1)-xval(it1)
	    fract=(xp(iter)-xval(it1))/delta
            do l=1,numsp
            concp(iter,l)=conc(it1,l)+fract*(conc(it1+1,l)-conc(it1,l))
            enddo
	 endif
         enddo
	 go to 10
      endif
      write(6,*) '** critical error in int1v',xp(iter),xval(num1)
10    continue
      return
      end
c****************************************************************************
      subroutine aq_int1v2(numsp,num1,conc,xval,num2,xp,concp)
c****************************************************************************
      dimension conc(numsp,1),xval(1),xp(1),concp(numsp,1)
      do 10 iter=1,num2
      if(xp(iter).le.xval(1)) then
         do l=1,numsp
         concp(l,iter)=0.0
         enddo
	 go to 10
      else if(xp(iter).ge.xval(num1)) then
         do l=1,numsp
         concp(l,iter)=0.0
         enddo
	 go to 10
      else
         do it1=1,(num1-1)
         if(xp(iter).ge.xval(it1).and.xp(iter).le.xval(it1+1))then
	    delta=xval(it1+1)-xval(it1)
	    fract=(xp(iter)-xval(it1))/delta
            do l=1,numsp
            concp(l,iter)=conc(l,it1)+fract*(conc(l,it1+1)-conc(l,it1))
            enddo
	 endif
         enddo
	 go to 10
      endif
      write(6,*) '** critical error in int1v',xp(iter),xval(num1)
10    continue
      return
      end
c****************************************************************************
      subroutine aq_int1m(numsp,num,conc,xval,xp,concp)
c****************************************************************************
      dimension conc(numsp,1),xval(1),concp(1)
      if(num.le.1) then
	 call aq_copy_mat(numsp,conc,concp)
	 return
      endif
      if(xp.le.xval(1)) then
	 call aq_copy_mat(numsp,conc(1,1),concp)
         return
      else if(xp.ge.xval(num)) then
	 call aq_copy_mat(numsp,conc(1,num),concp)
         return
      endif
      do iter=1,(num-1)
      if(xp.ge.xval(iter).and.xp.le.xval(iter+1))then
	 delta=xval(iter+1)-xval(iter)
	 fract=(xp-xval(iter))/delta
         do l=1,numsp
         concp(l)=conc(l,iter)+fract*(conc(l,iter+1)-conc(l,iter))
         enddo
	 return
      endif
      enddo
      write(6,*) '** critical error in aq_intm **'
      return
      end
c****************************************************************************
      subroutine aq_int2d(ix,iy,il,conc,num,xp,yp,concp)
c****************************************************************************
c     purpose:  find the value at (x,y) from the 2 dimensional field
c----------------------------------------------------------------------------
      dimension conc(ix,iy,1),xp(1),yp(1),concp(il,1)
      do iter=1,num
cc
         ip=int(xp(iter))
         jp=int(yp(iter))
	 if(ip.lt.1.or.ip.gt.ix.or.jp.lt.1.or.jp.gt.iy) then
	    write(6,*) '** critical error in intp2d **' 
	    write(6,*) '** out of the bound: ip=',ip,' jp=',jp
	    ip=1
	    jp=1
         endif
         frx=xp(iter)-ip
         fry=yp(iter)-jp
cc
          do l=1,il
          conc1=conc(ip,jp,l)
     1     +frx*(conc(ip+1,jp,l)  -conc(ip,jp,l))
          conc2=conc(ip,jp+1,l)
     1     +frx*(conc(ip+1,jp+1,l)-conc(ip,jp+1,l))
          concp(l,iter)=conc1+fry*(conc2-conc1)
          enddo
      enddo
      return
      end
c******************************************************************
      subroutine aq_r4toi2(n,rval,ival,iround)
c******************************************************************
      integer*2 ival
      dimension rval(1),ival(1)
      if(iround.eq.0) then
	 round=0.5               !  rounding
      else 
	 round=0.0               !  chopping
      endif
      do i=1,n
      ival(i)=int(rval(i)+0.5)
      enddo
      return
      end
c******************************************************************
      subroutine aq_r4tor8(n,rval,dval)
c******************************************************************
      dimension rval(1),dval(1)
ctyh      real*8 dval
      double precision dval
      do i=1,n
      dval(i)=rval(i)
      enddo
      return
      end
c******************************************************************
      subroutine aq_r8tor4(n,dval,rval)
c******************************************************************
      dimension rval(1),dval(1)
c      real*8 dval
      double precision dval
      do i=1,n
      rval(i)=dval(i)
      enddo
      return
      end
c***************************************************************
      subroutine aq_cor(n,x,y,corval)
c***************************************************************
      dimension x(1),y(1)
      call aq_ave(n,x,xave)
      call aq_ave(n,y,yave)
      deno1=0.
      deno2=0.
      rnum=0.
      do i=1,n
         deno1=deno1+(x(i)-xave)*(x(i)-xave)
         deno2=deno2+(y(i)-yave)*(y(i)-yave)
         rnum=rnum+(x(i)-xave)*(y(i)-yave)
      enddo
      deno=sqrt(deno1*deno2)
      corval=rnum/deno
c  
      return
      end
c*************************************************************
      subroutine aq_ave(n,x,xave)
c************************************************************
      dimension x(n)
      xave=0.
      do i=1,n
         xave=xave+x(i)
      enddo
      xave=xave/float(n)
      return
      end
c*************************************************************
      subroutine aq_max_i4(n,x,xmax)
c************************************************************
      dimension x(n)
      integer xmax,x
      xmax=x(1)
      do i=2,n
      if(x(i).gt.xmax) xmax=x(i)
      enddo
      return
      end
c*************************************************************
      subroutine aq_max(n,x,xmax)
c************************************************************
      dimension x(n)
      xmax=x(1)
      do i=2,n
      if(x(i).gt.xmax) xmax=x(i)
      enddo
      return
      end
c*************************************************************
      subroutine aq_min(n,x,xmin)
c************************************************************
      dimension x(n)
      xmax=x(1)
      do i=2,n
      if(x(i).lt.xmin) xmin=x(i)
      enddo
      return
      end
c****************************************************
      subroutine aq_i2toi4(n,i2conc,i4conc)
c****************************************************
      dimension i2conc(1),i4conc(1)
      integer*2 i2conc
      do i=1,n
      i4conc(i)=i2conc(i)
      enddo
      return
      end
c****************************************************
      subroutine aq_i2tor4(n,i2conc,r4conc)
c****************************************************
      dimension i2conc(1),r4conc(1)
      integer*2 i2conc
      do i=1,n
      r4conc(i)=i2conc(i)
      enddo
      return
      end
c*******************************************************
      subroutine aq_mat_swap(n,m,a,b)
c*******************************************************
      dimension a(n,m),b(m,n)
      do i=1,m
      do j=1,n
c     a(i,j)=b(j,i)
      b(i,j)=a(j,i)
      enddo
      enddo
      return
      end
c*******************************************************
      subroutine aq_fudge(rconc,ist,iend,conc,iflag,minpt)
c*******************************************************
c     minpt: minimum number of valid data points 
c---------------------------------------------------------------
      dimension rconc(1),conc(1),idx(100),temp(100)
      real rconc,isum
      if(iflag.ne.0) return
      isum=0.
      num=0
      do i=ist,iend
      if(rconc(i).ne.9999) then
	 isum=isum+rconc(i)
         num=num+1
      end if
      enddo
      if(isum.eq.0.and.num.lt.minpt) then
         iflag=2
         return
      endif
      ave=isum/float(num)
c
      if(rconc( ist).eq.9999) rconc( ist)=rconc( ist+1)
      if(rconc(iend).eq.9999) rconc(iend)=rconc(iend-1)
      if(rconc( ist).eq.9999) rconc( ist)=ave
      if(rconc(iend).eq.9999) rconc(iend)=ave
      conc(1)=rconc(ist)
      conc(iend-ist+1)=rconc(iend)
c 
c     Interpolation Begins Here
      num=0
      do i=ist,iend
      if(rconc(i).ne.9999) then
	 num=num+1
         temp(num)=rconc(i)
         idx(num)=i-ist+1
      end if
      enddo
      do i=2,(iend-ist)
      do j=1,(num-1)
      if(i.ge.idx(j).and.i.lt.idx(J+1)) then
	 factor=float(i-idx(j))/float(idx(j+1)-idx(j))
	 conc(i)=temp(j)+(temp(j+1)-temp(j))*factor
      end if
      enddo
      enddo
      return
      end
c***********************************************************************
      subroutine aq_chop(n,ch_string)
c*************************************************************
      character*1 cdum(80),ch_string(1)
      do i=1,n
      cdum(i)=ch_string(i)
      enddo
      call aq_blank(n,ch_string)
      k=0
      do i=1,n
      if(cdum(i).ne.' '.and.cdum(i).ne."'".
     1                  and.cdum(i).ne.'"') then
	 k=k+1
	 ch_string(k)=cdum(i)
      endif
      enddo
      return
      end
c*****************************************************************************
      subroutine aq_read_conc(iunit,num,numsp,numv,sname,conc,cherr)
c****************************************************************************
      character*(*) cherr,sname
      character*200 cdum1
      dimension conc(num,1),sname(1) 
c-------------------------------------------------------------------------
      nchar=len(sname(1))
      do l=1,numv
         read(iunit,*,end=35) cdum1(1:nchar)
         call aq_find(numsp,cdum1(1:nchar),sname,lpsec,iflag)
cc
	 if(iflag.eq.1) then
	    write(6,*) cherr
	    write(6,*) 'The chem. species does not have  ',
     1                                           cdum1(1:nchar)
	    stop
         else
            call aq_read(iunit,num,conc(1,lpsec))
	 endif    
      enddo
35    continue
c---------------------------------------------------------------------
      return
      end
c**********************************************************************
      subroutine aq_read_bc
     1        (iunit,ix,iy,iz,numsp,numv,sname,sx,sy,sz,cherr)
c****************************************************************************
      character*(*) cherr,sname
      character*100 cdum1
      dimension sx(iy,iz,2,1),sy(ix,iz,2,1),sz(ix,iy,1),sname(1) 
c-------------------------------------------------------------------------
      nchar=len(sname(1))
      do l=1,numv
         read(iunit,*,end=35) cdum1(1:nchar)
         call aq_find(numsp,cdum1(1:nchar),sname,lpsec,iflag)
cc
	 if(iflag.eq.1) then
	    write(6,*) cherr
	    write(6,*) 'The chem species file does not have  ',cdum1
	    stop
         else
            if(ix.gt.1) 
     1              call aq_read(iunit,iy*iz*2,sx(1,1,1,lpsec))
            if(iy.gt.1) 
     1              call aq_read(iunit,ix*iz*2,sy(1,1,1,lpsec))
            if(iz.gt.1) 
     1              call aq_read(iunit,ix*iy,sz(1,1,lpsec))
	 endif    
      enddo
35    continue
c---------------------------------------------------------------------
      return
      end
c**********************************************************************
      subroutine aq_test_val_r4(n,val,crit,cherr,iflag)
c********************************************************************
      dimension val(1)
      character*(*) cherr
      iflag=0
      do i=1,n
      if(abs(val(i)-crit).le.1.e-30) then
	 write(6,*) cherr,i
	 iflag=1
         return 
      endif
      enddo
      return
      end
c***************************************************************
      subroutine prtspec(iunit,n,spname,val)
c**************************************************************
      character*(*) spname(1)
      dimension val(1)
      icol=6
      last=mod(n,icol)
      num=n/icol
      do l=1,num
      write(iunit,100) (spname((l-1)*icol+i),i=1,icol)
      write(iunit,200) (val((l-1)*icol+i),i=1,icol)
      enddo
      write(iunit,100) (spname(num*icol+i),i=1,last)
      write(iunit,200) (val(num*icol+i),i=1,last)
100   format(3x,6a14)
200   format(6(e11.3,3x))
      return
      end
c***************************************************************
      subroutine prtspec_r8(iunit,n,spname,val)
c**************************************************************
      character*(*) spname(1)
      dimension val(1)
      double precision val
      icol=6
      last=mod(n,icol)
      num=n/icol
      do l=1,num
      write(iunit,100) (spname((l-1)*icol+i),i=1,icol)
      write(iunit,200) (val((l-1)*icol+i),i=1,icol)
      enddo
      write(iunit,100) (spname(num*icol+i),i=1,last)
      write(iunit,200) (val(num*icol+i),i=1,last)
100   format(3x,6a14)
200   format(6(e11.3,3x))
      return
      end
c*****************************************************************
      subroutine chem_err(i,j,k,nw,valj,zen)
c****************************************************************
      include 'aqmax.param'
      include 'aqcon5.cmm'
      include 'aqrxn.cmm'
      include 'tuv.params'
      common /chemer/chem_method,istate
      common /aq_tuv/wl(kw),wc(kw),wu(kw),f(kw),albedo(kw),xso3(kw),
     1      s226(kw),s263(kw),s298(kw),xsso2(kw),xsno2(kw),tlev(kz),
     2      tlay(kz),airlev(kz),cz_air(kz),cz_ozo(kz),cz_no2(kz),
     3      cz_so2(kz),cz_aero(kz),cz_cloud(kz),z(kz),dtaer(kz,kw),
     4      omaer(kz,kw),gaer(kz,kw),dtno2(kz,kw),dtso2(kz,kw),
     5      dto3(kz,kw),dtcld(kz,kw),omcld(kz,kw),gcld(kz,kw),
     6      dtrl(kz,kw)
      common /aq_tuv2/nz,ini_tuv,nw0
      character*80 chem_method
      dimension valj(1)
c     write(6,100) (valj(i1),i1=1,nw)
      
      if(istate.gt.0) then
c      write(6,'(1x,6e10.3)') (kg(i1),i1=1,nrxng)
      return
      endif
      
      write(6,*) '***  critical error in chemistry solver **'
      write(6,*) istate
      write(6,*) 'i,j,k',i,j,k
      write(6,*) 'rate constant'
      write(6,'(1x,6e10.3)') (kg(i1),i1=1,nrxng)
      write(6,*) 'photo-related data'
      write(6,*) 'nw=',nw,'zen=',zen
      write(6,100) (valj(i1),i1=1,nw)
      do i1=1,nw
      if(valj(i1).gt.1.e+25.or.valj(i1).le.1.e-25)  
     1               write(6,*) 'valj inf', i1
      enddo
      if(ini_tuv.ne.59819) then
       print*,'chem_err stop'
       stop
      endif 
      sumaer  =0.0
      sumomaer=0.0
      sumgaer =0.0
      sumdtno2=0.0
      sumdtso2=0.0
      sumdto3 =0.0
      sumdtcld=0.0
      sumomcld=0.0
      sumgcld =0.0
      sumdtrl =0.0
      do i1=1,nz
      do i2=1,(nw-1)
      sumaer  =sumaer  +dtaer(i1,i2)
      sumomaer=sumomaer+omaer(i1,i2)
      sumgaer =sumgaer +gaer(i1,i2)
      sumdtno2=sumdtno2+dtno2(i1,i2)
      sumdtso2=sumdtso2+dtso2(i1,i2)
      sumdto3 =sumdto3 +dto3(i1,i2)
      sumdtcld=sumdtcld+dtcld(i1,i2)
      sumomcld=sumomcld+omcld(i1,i2)
      sumgcld =sumgcld +gcld(i1,i2)
      sumdtrl =sumdtrl +dtrl(i1,i2)
      enddo
      enddo
      write(6,100) sumaer,sumomaer,sumdtno2,sumdtso2
      write(6,100) sumdto3,sumdtcld,sumomcld,sumgcld
      write(6,100) sumdtrl,sumgaer
100   format(1x,6e10.3)
      stop
      end
c***********************************************************
      subroutine set_nan_r4(n,val)
c************************************************************
      dimension val(1)
      iopt=1
      do 10  i=1,n 
      if(abs(val(i)).le.1.e-31) go to 10
      dum1=abs(val(i))*0.1
      dum2=abs(val(i))*0.1
      if(dum1.ne.dum2) then
	 if(iopt.eq.0.or.i.eq.1.or.i.eq.n) then
       	   val(i)=0.0
         else
	    val(i)=(val(i-1)+val(i+1))/2.0
         endif
      endif
10    continue
      return
      end
c***********************************************************
      subroutine find_nan_1d(n,nan,val)
c************************************************************
      dimension val(1)
      nan=0
      do 10  i=1,n 
      if(abs(val(i)).le.1.e-31) go to 10
      dum1=abs(val(i))*0.1
      dum2=abs(val(i))*0.1
      if(dum1.ne.dum2) nan=i
10    continue
      return
      end
