c***********************************************************************
      subroutine rxn(ix,iy,iz,bounds,idate,ut,dt,tlon,tlat,h,hdz,
     1    sg1,sl1,sp1,t,wc,wrc,sprc,rvel,cldod,kctop,ccover,dobson)
c**********************************************************************
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-----
      include 'aqmax.param'
      include 'aqcon1.cmm'
      include 'aqcon5.cmm'
      include 'aqindx.cmm'
      include 'aqrxn.cmm'
      include 'tuv.params'
cc
      integer :: bounds(8)
      dimension sg1(ix,iy,iz,*),t(ix,iy,*),wc(ix,iy,*),wrc(ix,iy,*)
      dimension sl1(ix,iy,iz,*),tlon(ix,*),tlat(ix,*),sprc(ix,*)
      dimension sp1(ix,iy,iz,*),adum1(mxspg),adum2(mxspg),toms(1)
      dimension rvel(ix,iy,*),idate(3),oneh(mxgr),onet(mxgr)
      real  cldod(ix,iy,iz),ccover(ix,iy,iz),dobson(ix,iy),fixwet(4)
      dimension h(ix,iy),hdz(ix,iy,*),one(mxgr*mxspg),yt(mxspg)
      real  valj(kj,kz),zstem(iz),o3stem(iz),so2stem(iz),aerod(iz,6),
     1 no2stem(iz),airstem(iz),tempstem(iz),cld1d(iz),cc1d(iz),
     2 mssaero(iz,5),rhu(iz)
      
      real kctop(ix,iy)
      integer jtdate 
c      real co_kg(ix,iy,iz), dtw
      dimension wval(kj),wint(kj)   
      real*8 fixspec(5)

      !real c303,c302
      !parameter(C303=19.83,C302=5417.4)

      integer :: xstart,xend,ystart,yend
      integer :: zstart,zend,nstart,nend
!
 
      !ESAT(TEMK)=.611*EXP(C303-C302/TEMK)       ! for calculating saturated water vapor pressure  
      !QSAT(ESAT1,PCB)=ESAT1*.622/(PCB-ESAT1)    ! TEMK is ambient temperature in K, PCB is the pressue in KPa
                                                ! QSAT is the saturated humidity in kg/kg

c
!    The computational bounds
      xstart = bounds(1); xend = bounds(2)
      ystart = bounds(3); yend = bounds(4)
      zstart = bounds(5); zend = bounds(6)
      nstart = bounds(7); nend = bounds(8)




!      call aq_pos(ix*iy*iz*numl(1,1),sg1)    ! remove negative concentrations


      jtdate=mod(idate(1),100)*10000+idate(2)*100+idate(3)   ! YYMMDD
      
      
      do i=xstart,xend
      do j=ystart,yend
cc    
      
        call update_photolysis(ix,iy,iz,i,j,idate,ut,dt,tlon,tlat,h,hdz,
     1    sg1,sl1,sp1,t,wc,wrc,sprc,rvel,cldod,kctop,ccover,dobson,valj)

      do lj=1,npht
      ! sg1(i,j,k,jno2+lj-1)=kg(ipht(lj,1))  ! store J-value
      enddo

      nrain=nint(kctop(i,j))
      

      do k=zstart,zend

      call rkg1(ix,iy,iz,i,j,k,sg1,wc,wrc,t(i,j,k),valj(1,k))

      
      yt(1:numl(1,3))=max(0.0,sg1(i,j,k,1:numl(1,3)))


      fixspec(1)=sg1(i,j,k,iair)
      fixspec(2)=sg1(i,j,k,io2)
      fixspec(3)=sg1(i,j,k,ih2o)
      fixspec(4)=sg1(i,j,k,ih2)
      fixspec(5)=sg1(i,j,k,ich4)


      call chem_box_gas(numl(1,2),dt,yt,fixspec,kg)

      sg1(i,j,k,1:numl(1,2))=max(0.0,yt(1:numl(1,2)))

      end do ! k
      end do ! j
      end do ! i

      return
      end

c***********************************************************************
      subroutine update_photolysis(ix,iy,iz,i,j,idate,ut,
     &     dt,tlon,tlat,h,hdz,
     1    sg1,sl1,sp1,t,wc,wrc,sprc,rvel,cldod,kctop,ccover,dobson,valj)
c**********************************************************************
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-----
      include 'aqmax.param'
      include 'aqcon1.cmm'
      include 'aqcon5.cmm'
      include 'aqindx.cmm'
      include 'aqrxn.cmm'
      include 'tuv.params'
cc
      integer :: bounds(8)
      dimension sg1(ix,iy,iz,*),t(ix,iy,*),wc(ix,iy,*),wrc(ix,iy,*)
      dimension sl1(ix,iy,iz,*),tlon(ix,*),tlat(ix,*),sprc(ix,*)
      dimension sp1(ix,iy,iz,*),adum1(mxspg),adum2(mxspg),toms(1)
      dimension rvel(ix,iy,*),idate(3),oneh(mxgr),onet(mxgr)
      real  cldod(ix,iy,iz),ccover(ix,iy,iz),dobson(ix,iy),fixwet(4)
      dimension h(ix,iy),hdz(ix,iy,*),one(mxgr*mxspg),yt(mxspg)
      real  valj(kj,kz),zstem(iz),o3stem(iz),so2stem(iz),aerod(iz,6),
     1 no2stem(iz),airstem(iz),tempstem(iz),cld1d(iz),cc1d(iz),
     2 mssaero(iz,5),rhu(iz)
      
      integer kctop(ix,iy),jtdate
      
c      real co_kg(ix,iy,iz), dtw
      dimension wval(kj),wint(kj)   

      real c303,c302
      parameter(C303=19.83,C302=5417.4)

      integer :: xstart,xend,ystart,yend
      integer :: zstart,zend,nstart,nend
!
 
      ESAT(TEMK)=.611*EXP(C303-C302/TEMK)       ! for calculating saturated water vapor pressure  
      QSAT(ESAT1,PCB)=ESAT1*.622/(PCB-ESAT1)    ! TEMK is ambient temperature in K, PCB is the pressue in KPa
                                                ! QSAT is the saturated humidity in kg/kg


!      call aq_pos(ix*iy*iz*numl(1,1),sg1)    ! remove negative concentrations


      jtdate=mod(idate(1),100)*10000+idate(2)*100+idate(3)   ! YYMMDD
      
cc    
      o3stem(1:iz)=sg1(i,j,1:iz,io3)         ! load STEM variable to one dimension
      so2stem(1:iz)=sg1(i,j,1:iz,iso2)
      no2stem(1:iz)=sg1(i,j,1:iz,ino2)
      airstem(1:iz)=sg1(i,j,1:iz,iair)
      zstem(1:iz)=hdz(i,j,1:iz)/1000.  ! convert to km
      tempstem(1:iz)=t(i,j,1:iz)
      cld1d(1:iz)=cldod(i,j,1:iz)
      cc1d(1:iz)=ccover(i,j,1:iz)

      do k=1,iz                  ! calculate relative humidity
       wvapor=sg1(i,j,k,ih2o)/sg1(i,j,k,iair)*18./28.9          ! water vapor in Kg/Kg
       airpress=sg1(i,j,k,iair)/2.687e+19*101300./273.*t(i,j,k) ! pressure in Pa
       rhu(k)=100.*wvapor/QSAT(ESAT(         ! computing relative humudity
     1		t(i,j,k)),airpress/1000)         !convert to KPa
       rhu(k)=amin1(rhu(k),100.)
       
       so2=sg1(i,j,k,iso2)/sg1(i,j,k,iair)*1e9   ! so2 in ppbv
       mssaero(k,1)=sg1(i,j,k,idust1)+ sg1(i,j,k,idust2)/125.
                                                     ! load aerosol for TUV, currently consider fine dust only 						
       mssaero(k,2)=sg1(i,j,k,iso4)/amin1(exp(so2/30.),2.)      ! 1=dust,2=watersoluble (50%sulfate),
       mssaero(k,3)=sg1(i,j,k,ibc)              !  3=black carbon, 4=sea salt,   5=organic carbon
c     1     + amax1((sg1(i,j,k,ico)-sg1(i,j,k,iair)*1e-7)/50.,0.)  ! use difference between current CO and background
                                                                 ! as an indicator of backgroup	BC 
       mssaero(k,4)=sg1(i,j,k,issf)+sg1(i,j,k,issc)/125.  
       mssaero(k,5)=sg1(i,j,k,ioc)+sg1(i,j,k,iopm25)+      ! opm25 and opm10 merged inte OC
     1   sg1(i,j,k,iopm10)/125.                           ! for computing optical depth 
      enddo
  
      call tuvsub(iz,zstem,o3stem,so2stem,no2stem,airstem,tempstem,
     1 cld1d,cc1d,mssaero,rhu,aerod,ut,jtdate,tlat(i,j),tlon(i,j),
     2  dobson(i,j),valj)                  ! valj include 30 J-values for all layers               

      do k=1,iz              ! check J-value
       do lj=1,npht
       	if(abs(valj(lj,k)).gt.1) then
	  print*, 'wrong J values LJ,K, valj(LJ,k), jtdate, ut=',
     1	  LJ,K,valj(LJ,k), jtdate, ut
          print*,'mssaero=',mssaero
	  print*,'aerod=',aerod
         stop
	endif  
       enddo
      enddo 	

      end subroutine update_photolysis



c***********************************************************************
      subroutine chem_box_gas(numeq,dt,yt,fixconc,rkg)
c**********************************************************************
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-----
      include 'saprcnov_stem.h'
      real :: yt(numeq), dt
      double precision :: rkg(NREACT)
      double precision :: abstol(NVAR),reltol(NVAR)
      double precision :: conc(NVAR), fixconc(NFIX)
      double precision :: t, tout, hmin, hmax, hstart 
      integer :: info(5)
!      integer, dimension(NVAR), parameter :: permute = (/
!     * 76, 28, 81, 85, 82, 31, 32, 56, 37, 20,
!     *  1, 53, 70, 69, 73, 55, 74, 13, 41, 30,
!     * 14, 15, 59, 54, 23, 46, 49, 66, 64, 68,
!     * 75, 48, 45, 44, 57, 61, 22, 29, 47, 58,
!     * 33, 42, 35, 43, 36, 65, 67, 63, 71, 51,
!     * 52, 24, 25, 26, 27,  2,  3,  4,  5,  6,
!     *  7,  9, 10, 11,  8, 12, 16, 88, 17, 18,
!     * 19, 72, 21, 79, 78, 83, 38, 50, 87, 62,
!     * 86, 39, 84, 77, 40, 60, 80, 34 /)
      integer, dimension(NVAR), parameter :: permute = (/
     & I_O3, I_H2O2, I_NO, I_NO2, I_NO3, I_N2O5, I_HONO,              
     & I_HNO3, I_HNO4, I_SO2, I_H2SO4, I_CO , I_HCHO ,             
     & I_CCHO, I_RCHO, I_ACET, I_MEK, I_HCOOH, I_MEOH,              
     & I_ETOH, I_CCO_OH, I_RCO_OH, I_GLY, I_MGLY, I_BACL,              
     & I_CRES, I_BALD, I_ISOPROD, I_METHACRO, I_MVK, I_PROD2,            
     & I_DCB1, I_DCB2, I_DCB3, I_ETHENE, I_ISOPRENE, I_C2H6,             
     & I_C3H8, I_C2H2, I_C3H6, I_ALK3 , I_ALK4, I_ALK5, I_ARO1,             
     & I_ARO2, I_OLE1, I_OLE2, I_TERP, I_RNO3, I_NPHE, I_PHEN,             
     & I_PAN, I_PAN2, I_PBZN, I_MA_PAN, I_BC, I_OC, I_SSF, I_SSC,              
     & I_PM10, I_PM25, I_DST1, I_DST2, I_DST3, I_DMS, I_CO2,              
     & I_CCO_OOH, I_RCO_O2, I_RCO_OOH, I_XN, I_XC, I_O3P, I_O1D,              
     & I_OH, I_HO2, I_C_O2, I_COOH, I_ROOH, I_RO2_R, I_R2O2, I_RO2_N,           
     & I_HOCOO, I_CCO_O2, I_BZCO_O2, I_BZNO2_O, I_BZ_O, I_MA_RCO3,         
     & I_TBU_O    /)       

c---------------------------------------------------------------------
c    input parameter for lsode
      t =  0.d0
      tout = dt
      reltol(1:NVAR) = 1.d-02
      abstol(1:NVAR) = 1.d+04
      hmin=1.d-06
      hmax=tout-t
      hstart=1.d0
      info(1)=1 ! 1=autonomous
      info(2)=0 ! 1 = 3rd order embedded, 0 = 1st order embedded

      
      do i=1,NVAR
        conc(permute(i)) = yt(i)
      end do
c
      call ros2(numeq,t,tout,hmin,hmax,hstart,
     &           conc,fixconc,rkg,
     &           abstol,reltol,info)

      do i=1,NVAR
        yt(i) = conc(permute(i))
      end do

      return
      end


c***********************************************************************
      subroutine diff_g(num,t,y,dkg,pkg)
c***********************************************************************
      include 'aqmax.param'
      include 'aqcon1.cmm'
      include 'aqcon5.cmm'
      include 'aqrxn.cmm'
      dimension y(1),wrg(mxrxn),dkg(1),pkg(1),yg(mxspg)
      double precision wrg,y,dkg,pkg,t,yg
cc
      call aq_copy_mat_r8(numl(1,2),y,yg)
      call aq_copy_mat_r8(numl(1,3)-numl(1,2),
     1          y(numl(1,2)+numl(3,1)*nbin+1),yg(numl(1,2)+1))
      call wdotig(num,nrxng,kg,wrg,yg,irs,nus)
      call rmolg(num,nrxng,wrg,dkg,pkg,irs,irr,irc,
     1           nus,nur,coesg,coerg)
cc
      if(nbin.ne.0) call rmolp(dkg,pkg,yg,y)
c      
      return
      end
c***********************************************************************
      subroutine wdotig(num,nrxng,kg,wr,yg,irs,nus)
c**********************************************************************
      implicit double precision(a-h,o-z)
      dimension kg(1),wr(1),yg(1),irs(1),nus(30,1)
      double precision kg
      do l=1,nrxng
         wr(l)=kg(l)
         i1=irs(l)
         do i=1,i1
         wr(l)=wr(l)*yg(nus(i,l))
         enddo
      enddo
      return
      end
c**********************************************************************
      subroutine rmolg(num,nrxng,wr,dkg,pkg,irs,irr,irc,
     1     nus,nur,coes,coer)
c**********************************************************************
      implicit double precision(a-h,o-z)
      parameter (mxspg=250,mxrxn=500)
      dimension wr(1),dkg(1),pkg(1)
      dimension irs(1),irr(1),nus(30,1),nur(30,1),irc(1)
      dimension coes(30,1),coer(30,1)
c
      call aq_zero_r8(num,dkg)
      call aq_zero_r8(num,pkg)
c
      do 20 l=1,nrxng
c
      if(irc(l).eq.0) then
         i1=irs(l)
         do i=1,i1
	 l1=nus(i,l)
         dkg(l1)=dkg(l1)+wr(l)
	 enddo

         i1=irr(l)
         do i=1,i1
	 l1=nur(i,l)
         pkg(l1)=pkg(l1)+wr(l)
	 enddo
      else
         i1=irs(l)
         do i=1,i1
	 l1=nus(i,l)
         dkg(l1)=dkg(l1)+coes(i,l)*wr(l)
	 enddo

         i1=irr(l)
         do i=1,i1
	 l1=nur(i,l)
         pkg(l1)=pkg(l1)+coer(i,l)*wr(l)
	 enddo
      endif
cc
20    continue
      return
      end
c******************************************************************
      subroutine rmolp(dkg,pkg,yg,y)
c******************************************************************
      include 'aqmax.param'
      include 'aqcon1.cmm'
      include 'aqcon5.cmm'
      double precision dkg,pkg,yg,y,rc,hetero
      dimension dkg(1),pkg(1),yg(1),y(1),rc(1000)
      call get_aero_rc(numl(3,1),y(numl(1,2)+1),rc)
      do ibin=1,nbin
         call dfcons_aero(numl(1,2),rc(ibin))
	 dkg(numl(1,2)+ibin)=0.0
	 pkg(numl(1,2)+ibin)=0.0
         do l=1,numl(1,2)
	     hetero=dif(l)*aero_num(ibin)
c            if(l.eq.1) write(6,'(a,i3,3e10.3)')  
c    1           'ibin',ibin,hetero,dif(l),aero_num(ibin)
             dkg(l)=dkg(l)+hetero*yg(l)
             pkg(numl(1,2)+ibin)=pkg(numl(1,2)+ibin)+
     1               hetero*yg(l)/6.02d+23*rmw(l)*1.d+12   ! unit=microgam/m3
         enddo
      enddo
      return
      end
c**********************************************************************
      subroutine senjac_g(num,t,y,pd)           ! for lsode
c***********************************************************************
      parameter (mxspg=250,mxrxn=500)
      include 'aqrxn.cmm'
      include 'aqcon5.cmm'
      dimension pd(num,1),y(1)
      dimension dwdc(mxrxn,mxspg)
      double precision t,y,pd,dwdc,prod
      call aq_zero_r8(num*num,pd)
      do i=1,nrxng
        do j=1,num
         dwdc(i,j)=0.d0
	enddo
      enddo	  
 
      do l=1,nrxng
       i1=irs(l)
       do i=1,i1
        j=nus(i,l)
        if(j.le.num)  then
          prod=kg(l)
          do j1=1,i1
           if(i.ne.j1) prod=prod*y(nus(j1,l))
	  enddo
          dwdc(l,j)=dwdc(l,j)+prod
        endif
       enddo
      enddo

c
      do 60 l=1,nrxng
      i1=irs(l)
      do 65 i=1,i1
      l1=nus(i,l)
      if(l1.le.num) then
         do l2=1,i1
	 l3=nus(l2,l)
	 if(l3.le.num) pd(l1,l3)=pd(l1,l3)-coesg(i,l)*dwdc(l,l3)
	 enddo
      endif
 65   continue
c
      i2=irr(l)
      do 67 i=1,i2
      l1=nur(i,l)
      if(l1.le.num) then
         do l2=1,i1
         l3=nus(l2,l)
         if(l3.le.num) pd(l1,l3)=pd(l1,l3)+coerg(i,l)*dwdc(l,l3)
	 enddo
      endif
 67   continue
 60   continue
      return
      end  
c***********************************************************************
      subroutine solpam(nw,npht,wval,val,wls,absorb,qy,ipht)
c***********************************************************************
      dimension wval(1),val(30,1),wls(30,1),absorb(30,1),qy(30,1)
      dimension ipht(30,2)
      do l=1,npht
      num=ipht(l,2)-1
      iter0=1
      do int=1,nw
      wlm=wval(int)
      if(wlm.lt.wls(l,1).or.wlm.gt.wls(l,num)) then 
	 val(l,int)=0.0
      else
         do 10 iter=iter0,(num-1)
         if(wlm.ge.wls(l,iter).and.wlm.le.wls(l,iter+1)) then
            delta=wls(l,iter+1)-wls(l,iter)
            dev=wlm-wls(l,iter)
            qval=qy(l,iter)+(qy(l,iter)-qy(l,iter+1))*dev/delta
            qabs=absorb(l,iter)+(absorb(l,iter)-absorb(l,iter+1))
     1  	    *dev/delta
            val(l,int)=qval*qabs
	    iter0=iter
	    go to 20
         end if
10       continue
         write(6,*) '** critical error in solpam **'
      endif
20    continue
      enddo
      enddo
      return
      end
c***********************************************************************
      subroutine rkg1(ix,iy,iz,i,j,k,sg1,wc,wrc,tp1,valj)
c***********************************************************************
      include 'aqmax.param'
      include 'aqcon1.cmm'
      include 'aqcon5.cmm'
      include 'aqindx.cmm'
      include 'aqrxn.cmm'
      real sg1(ix,iy,iz,*), wc(ix,iy,iz), wrc(ix,iy,iz)
      dimension yt(1),wval(1),wint(1),valj(1),val(30)
      double precision rk0,rk2,rk3,rkinf
cc
      ruc=0.0019872
	 do l=1,npht
           kg(ipht(l,1))=valj(l)
	 enddo

      yt(1:numl(1,3))=sg1(i,j,k,1:numl(1,3))
cc    
      do 10 l=1,nrxng
      if(irt(l).eq.1) then           ! K(T)=A*(T/300)**B*exp(-Ea/R/T) 
	 kg(l)=aa(1,l)*(tp1/300.)**bb(1,l)*dexp(-ee(1,l)/tp1)
      else if(irt(l).eq.3) then      ! K(T)=K0+K2*[M]
	 rk0  =aa(1,l)*(tp1/300.)**bb(1,l)*dexp(-ee(1,l)/(tp1))
	 rkinf=aa(2,l)*(tp1/300.)**bb(2,l)*dexp(-ee(2,l)/(tp1))
	 kg(l)=rk0 + rkinf*yt(iair)
      else if(irt(l).eq.4) then      ! falloff reaction
	 rk0  =aa(1,l)*(tp1/300.)**bb(1,l)*dexp(-ee(1,l)/(tp1))
	 rkinf=aa(2,l)*(tp1/300.)**bb(2,l)*dexp(-ee(2,l)/(tp1))
	 dum1=dlog10(rk0*yt(iair)/rkinf)/fnn(l)
	 gf=1./(1.+dum1*dum1)
	 kg(l)=(rk0*yt(iair)/(1.+(rk0*yt(iair)/rkinf)))*ff(l)**gf
      end if
10    continue

      rk0=7.2d-15*dexp(785.d0/tp1)
      rk2=4.1d-16*dexp(1440.d0/tp1)
      rk3=1.9d-33*dexp(725.d0/tp1)*yt(iair)
      kg(27)=rk0+rk3/(1+rk3/rk2)

      kg(219) = 2.7d-6*(1+                                      ! additional reaction from SO2 to SO4
     1  (wc(i,j,k)+wrc(i,j,k))*28.966/6.022e23*1e6*
     &      sg1(i,j,k,iair)   ! convert Cloud and Rain water from Kg/Kg to g/m3
     2   *1000.*0.082*tp1                               ! based on equation 6.75 in Seinfeld book 
     3   /4.4e5 )                                                 ! "atmospheric chemistry and physics"
          ! 4.4e5= P(pa)/temp(k)*1.2027e3 (ppbv->moles/L)
                                                        
      kg(220) = sg1(i,j,k,nwetso2)/100   ! HO2 wet removing
      kg(221) = sg1(i,j,k,nwetso2)       ! SO2 wet removing      
      kg(222) = sg1(i,j,k,nwetso4)/2.    ! SO4 wet removing      
      kg(223) = sg1(i,j,k,nwethno3)      ! HNO3 wet removing
      kg(224) = sg1(i,j,k,nweth2o2)      ! H2O2 wet removing
      kg(225) = 0.                       ! BC wet removing, using sulfate data
      kg(226) = 0.                       ! OC wet removing
      kg(227) = sg1(i,j,k,nwetso4)       ! SSF wet removing
      kg(228) = sg1(i,j,k,nwetso4)       ! SSC wet removing
      kg(229) = 0.                       ! PM10 wet removing
      kg(230) = 0.                       ! PM25 wet removing
      kg(231) = sg1(i,j,k,nwetso4)/2.    ! Dust1 wet removing
      kg(232) = sg1(i,j,k,nwetso4)/2.    ! Dust2 wet removing
      kg(233) = sg1(i,j,k,nwetso4)/2.    ! Dust3 wet removing
      kg(234) = 0.                       ! DMS wet removing                              

      kg(235) = 0.d0                       ! empty reaction of CO2

 
      return
      end
c********************************************************************
      subroutine aq_rain_top(ix,iy,iz,i,j,sprc,wr,nrain)
c**************************************************************************
      dimension wr(ix,iy,1)
      nrain=0
      if(sprc.gt.1.e-04) then
         do k=iz,1,-1
         if(wr(i,j,k).gt.1.e-10) then
            nrain=k
            if(k.eq.iz) then
	    write(6,*) '** critical error: wr(iz) .ne.0 **'
            stop
            endif
	    go to 10
         endif
         enddo
         write(6,*) '** critical error: wr=0.0 **'
      endif
10    continue
      return
      end
c**************************************************************************
      subroutine aq_wash(ix,iy,iz,i,j,dt,sprc,rvel,hdz,nrain,nwash)
c**************************************************************************
      dimension sprc(ix,1),rvel(ix,iy,1)
      dimension hdz(ix,iy,1)
      co_nm=0.2
      nwash=1
      if(sprc(i,j).gt.1.e-04) then
	 uave=0.
	 do k=1,nrain
	 uave=uave+rvel(i,j,k)
	 enddo
         uave=uave/nrain
	 if(uave.gt.0.01) then
	 delt=co_nm*hdz(i,j,nrain)/uave
	 nwash=int(dt/delt+0.95)
	 endif
      endif
      if(nwash.lt.1) nwash=1
10    continue
      return
      end
c***************************************************************
      subroutine get_aero_num(ndim,sp)
c************************************************************
      include 'aqmax.param'
      include 'aqcon1.cmm'
      integer ndim
      real sp(ndim,1)
      aero_den=1.7
      do ibin=1,nbin
      rc3=(10.**daero(ibin))**3/8.0
      airmas=0.0
      do l=1,numl(3,1)
      airmas=airmas+sp(l,ibin)
      enddo
      airmas=airmas*1.e-12    ! convert unit to gm/cm3
      aero_num(ibin)=airmas/(aero_den*1.33333*3.141592*rc3)
      enddo
      return
      end
c*****************************************************************
      subroutine get_aero_rc(ndim,sp,rc)
c**************************************************************
      include 'aqmax.param'
      include 'aqcon1.cmm'
      integer ndim
      dimension sp(ndim,1),rc(1)
      double precision sp,rc,airmas,aero_den
      aero_den=1.7
      do ibin=1,nbin
         airmas=0.0
         do l=1,numl(3,1)
         airmas=airmas+sp(l,ibin)
         enddo
         airmas=airmas*1.d-12   ! conver unit to gm/cm3
         if(airmas.le.1.d-30.or.aero_num(ibin).le.1.d-10) then
	    rc(ibin)=10.**daero(ibin)/2.0
         else
            rc(ibin)=airmas/
     1          (1.33333*3.141592*aero_den*aero_num(ibin))
	    rc(ibin)=rc(ibin)**0.33333
         endif
      enddo
      return
      end
c**************************************************************
      subroutine get_aero_rc_r4(ndim,sp,rc)
c**************************************************************
      include 'aqmax.param'
      include 'aqcon1.cmm'
      integer ndim
      dimension sp(ndim,1),rc(1)
      aero_den=1.7
      do ibin=1,nbin
         airmas=0.0
         do l=1,numl(3,1)
         airmas=airmas+sp(l,ibin)
         enddo
         airmas=airmas*1.e-12   ! conver unit to gm/cm3
         if(airmas.le.1.e-30) then
	    rc(ibin)=10.**daero(ibin)/2.0
         else
            rc(ibin)=airmas
     1        /(1.33333*3.141592*aero_den*aero_num(ibin))
            rc(ibin)=rc(ibin)**0.33333
         endif
      enddo
      return
      end
c************************************************************
      subroutine distr_aero(ndim,sp) 
c*************************************************************
*   Eqns:  n1v1+n2v2=nv,  n1+n2=n  --> n2v2/nv=(v-v1)/(v2-v1)(v2/v)
      include 'aqmax.param'
      include 'aqcon1.cmm'
      integer ndim
      dimension sp(ndim,1),adum1(500*100)
      dimension rc(1000),xval(1000),xp(1000)
      if(nbin.eq.0) return
      call get_aero_rc_r4(numl(3,1),sp,rc)
      call aq_zero_r4(nbin*numl(3,1),adum1)
      do i=1,nbin
         xp(i  )=3.141592*(10.**daero(i))**3/6.0
	 xval(i)=3.141592*rc(i)**3.0*1.333333
      enddo
      do 10 i=1,nbin
      if(aero_num(i).le.1.e-20) go to 10
      if(xval(i).le.xp(1)) then
	 do l=1,numl(3,1)
	 adum1(l)=sp(l,i)
	 enddo
      else 
	 do iter=1,nbin
         if(xval(i).gt.xp(iter).and.xval(i).le.xp(iter+1)) then
	    fract=(xval(i)-xp(iter))/(xp(iter+1)-xp(iter))*
     1            xp(iter+1)/xval(i)
	    do l=1,numl(3,1)
	    adum1(numl(3,1)*iter+l)=adum1(numl(3,1)*iter+l)
     1                             +sp(l,i)*fract
	    adum1(numl(3,1)*(iter-1)+l)
     1       =adum1(numl(3,1)*(iter-1)+1)+sp(l,i)*(1.-fract)
	    enddo
	 go to 10
	 endif
	 enddo
	 do l=1,numl(3,1)
	 adum1(numl(3,1)*(nbin-1)+l)=sp(l,i)
	 enddo
      endif
10    continue       
cc
c     write(6,*) 'distribut'
c     do i=1,nbin
c     write(6,'(i4,4e11.3)') i,xval(i),sp(1,i),xp(i),adum1(i)
c     enddo
      call aq_copy_mat(numl(3,1)*nbin,adum1,sp)
      return
      end
c*********************************************************************
      subroutine dfcons_aero(num,rc) ! unit of rc=cm
c*********************************************************************
      include 'aqmax.param'
      include 'aqcon1.cmm'
      include 'aqcon5.cmm'
      double precision rc
      path=0.651D-02*1.e-04                ! mean free path in centi-meter
      rkn=path/rc
      do l=1,num
      accom(l)=1.e-05
      ff=(1.333+0.71/rkn)/(1.+1./rkn)+4.*(1.-accom(l))/(3.*accom(l))
      dif(l)=dgas(l)/(1.+ff*rkn)*4.*3.141592*rc
      enddo
      return
      end
