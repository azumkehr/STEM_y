c***********************************************************************
      subroutine rxn_adjoint(ix,iy,iz,bounds,idate,ut,dt,
     &            tlon,tlat,h,hdz,
     &            sg1,Adjoint,sl1,sp1,t,wc,wrc,sprc,rvel,
     &            cldod,kctop,ccover,dobson)
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
      real  Adjoint(ix,iy,iz,*)
      dimension sl1(ix,iy,iz,*),tlon(ix,*),tlat(ix,*),sprc(ix,*)
      dimension sp1(ix,iy,iz,*),adum1(mxspg),adum2(mxspg),toms(1)
      dimension rvel(ix,iy,*),idate(3),oneh(mxgr),onet(mxgr)
      real  cldod(ix,iy,iz),ccover(ix,iy,iz),dobson(ix,iy),fixwet(4)
      dimension h(ix,iy),hdz(ix,iy,*),one(mxgr*mxspg),yt(mxspg)
      real Lambda(mxspg)
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



      fixspec(1)=sg1(i,j,k,iair)
      fixspec(2)=sg1(i,j,k,io2)
      fixspec(3)=sg1(i,j,k,ih2o)
      fixspec(4)=sg1(i,j,k,ih2)
      fixspec(5)=sg1(i,j,k,ich4)

      
      yt(1:numl(1,2))=sg1(i,j,k,1:numl(1,2))
      Lambda(1:numl(1,2))=Adjoint(i,j,k,1:numl(1,2))

      call chem_box_gas_adjoint(numl(1,2),dt,yt,Lambda,fixspec,kg)

      Adjoint(i,j,k,1:numl(1,2))=Lambda(1:numl(1,2))

      end do ! k
      end do ! j
      end do ! i

      return
      end




c***********************************************************************
      subroutine chem_box_gas_adjoint(numeq,dt,yt,adj,fixconc,rkg)
c**********************************************************************
c---+----1----+----2----+----3----+----4----+----5----+----6----+----7-----
      include 'saprcnov_stem.h'
      real :: yt(NVAR), adj(NVAR), dt
      double precision :: rkg(NREACT)
      double precision :: abstol(NVAR),reltol(NVAR)
      double precision :: conc(NVAR), fixconc(NFIX)
      double precision :: Lambda(NVAR)
      double precision :: t, tout, hmin, hmax, hstart 
      integer :: info(5)
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
        if ( yt(i)>0.0 ) then
          conc(permute(i))   = yt(i)
          Lambda(permute(i)) = adj(i)
	else
          conc(permute(i))   = 0.d0
          Lambda(permute(i)) = adj(i) !0.d0	
	end if  
      end do
c
      call ros2_cadj(numeq,t,tout,hmin,hmax,hstart,
     &           conc,Lambda,fixconc,rkg,
     &           abstol,reltol,info)

      do i=1,NVAR
        if ( conc(permute(i))>0.d0 ) then
          yt(i)  = conc(permute(i))
          adj(i) = Lambda(permute(i))
	else
          yt(i)  = 0.0
          adj(i) = Lambda(permute(i)) !0.0	
	end if  
      end do

      return
      end
