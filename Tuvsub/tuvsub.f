      subroutine tuvsub(nzstem,zstem,o3stem,so2stem,no2stem,airstem,
     1 tempstem,cldod,ccover,mssaero,rhu,aerod,ut,idate,alat,along,
     2 dobnew,valj)
*_______________________________________________________________________
*     Tropospheric Ultraviolet-Visible (TUV) radiation model
*     Version 4.1
*     July 2000 by Madronich et al.
*_______________________________________________________________________
*-----------------------------------------------------------------------------*
*= This program is free software;  you can redistribute it and/or modify     =*
*= it under the terms of the GNU General Public License as published by the  =*
*= Free Software Foundation;  either version 2 of the license, or (at your   =*
*= option) any later version.                                                =*
*= The TUV package is distributed in the hope that it will be useful, but    =*
*= WITHOUT ANY WARRANTY;  without even the implied warranty of MERCHANTIBI-  =*
*= LITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public     =*
*= License for more details.                                                 =*
*= To obtain a copy of the GNU General Public License, write to:             =*
*= Free Software Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.   =*
*-----------------------------------------------------------------------------*
*= To contact the authors, please mail to:                                   =*
*= Sasha Madronich, NCAR/ACD, P.O.Box 3000, Boulder, CO, 80307-3000, USA  or =*
*= send email to:  sasha@ucar.edu                                            =*
*-----------------------------------------------------------------------------*
*= Copyright (C) 1994,95,96,97,98,99  University Corporation for Atmospheric =*
*= Research                                                                  =*
*-----------------------------------------------------------------------------*

      IMPLICIT NONE

* include parameter file

      INCLUDE 'tuv.params'

* ___ SECTION 1: VARIABLES AND PARAMETERS ______________________________

      real cclean, cblock, cfull
      data cclean,cblock,cfull/0.15,0.5,0.85/ 
           !cloud coverage for defining clean sky, broken cloud with blocked cloud
	   ! and near-full cloud coverage

* input dimension from STEM
      integer nzstem              ! stem vertical grid number
      real zstem(nzstem),o3stem(nzstem),so2stem(nzstem),no2stem(nzstem),
     1 airstem(nzstem),tempstem(nzstem),cldod(nzstem),ccover(nzstem),
     2 mssaero(nzstem,5), rhu(nzstem) 
        ! stem Z-grid (km), o3(mole/cm3),so2, no2, air desity(mole/cm3), 
	! temperature (K), cloud optic depth, cloud coverage
	! aerosol concentration in molecular/cm3 (ug/m3),  1=dust,2=watersoluble(50%sulfate),
        !          3=black carbon, 4=sea salt,   5=organic carbon
        ! relative humidity 
         
* additional output dimension	
      real aerod(nzstem,6)  ! aerosol extinction coefficient
                            ! 6 represent total aerosol depth counted from top layer
			    ! 7 is the Single Scatter Albedo
       
* wavelength grid:

      INTEGER nw, iw
      REAL wl(kw), wc(kw), wu(kw)

* altitude grid

      INTEGER nz, iz
      REAL z(kz)

* solar zenith angle and azimuth
* slant pathlengths in spherical geometry

      REAL zen, azim
      INTEGER nid(0:kz)
      REAL dsdh(0:kz,kz)

* extra terrestrial solar flux and earth-Sun distance ^-2

      REAL f(kw), etf(kw)
      REAL esfact

* ozone absorption cross section and ozone optical depth:

      REAL xso3(kw), s226(kw), s263(kw), s298(kw)

* O2 absorption cross section

      REAL xso2(kz,kw)

* SO2 absorption cross section
     
      REAL xsso2(kw)

* NO2 absorption cross section
     
      REAL xsno2(kw)

* atmospheric optical parameters:

      REAL tlev(kz), tlay(kz)
      REAL airlev(kz), colinc(kz), vcol(kz), scol(kz)
      REAL dtrl(kz,kw)
      REAL dto3(kz,kw), dto2(kz,kw), dtso2(kz,kw), dtno2(kz,kw)
      REAL dtcld(kz,kw), omcld(kz,kw), gcld(kz,kw)
      REAL dtaer(kz,kw), omaer(kz,kw), gaer(kz,kw)
      REAL albedo(kw)

* spectral irradiance and actinic flux (scalar irradiance):

      REAL edir(kz), edn(kz), eup(kz), edir2(kz),edn2(kz), eup2(kz)
      REAL fdir(kz), fdn(kz), fup(kz), fdir2(kz),fdn2(kz), fup2(kz)

* spectral weighting functions and weighted radiation:

      INTEGER ns, is
      REAL sw(ks,kw), rate(ks), dose(ks)
      REAL sirrad, sprate 
      CHARACTER*40 label(ks)

*! j-values:

      INTEGER nj, ij
      REAL sj(kj,kz,kw), valj(kj,kz)
      REAL saflux, deltaj
      CHARACTER*40 jlabel(kj)

* new sea level pressure, surface dobson, etc.

      REAL pmbnew, dobnew, so2new, no2new

* Location and time

      REAL alat, along 
      INTEGER idate
      REAL dtime, ut, ut0

* Commonly used looping indices

      integer i, j
      integer  idat, idob, itime, izen

* Other user-defined variables here:

      integer nzen
      real  sav(100,kz), sza(100)
      
      logical first
      save first
      data first/.true./

***


* ___ SECTION 2: SET GRIDS _________________________________________________

* wavelengths

      if(first) CALL gridw(nw,wl,wc,wu)
      
* for version B only (no sr bands):  limit  to wavelengths longer than 250 nm.

C      IF (wl(iw) .LE. 250.) STOP

* altitudes (set the surface elevation, z(1), in km asl).

       z(1:nzstem)=zstem(1:nzstem)
       nz=nzstem
       do while(z(nz).le.78)  ! set top to 82km
        nz=nz+1
	z(nz)=z(nz-1)+2         ! set vertical interval 2km while above stem layers
       enddo 	
c       write(kout,"(' total z grid number:',i2)")nz
c       write(kout,"('Z grid ', 20f7.3)")z(1:nz)  
c       z(1) = 0.
c       CALL gridz(nz,z)

* ___ SECTION 3: SPECTRAL DATA ____________________________

* read (and grid) extra terrestrial flux data:
      
      if(first) then
      CALL rdetfl(nw,wl,f)

* read cross section data for 
*    ozone (temperature-dependent)
*    SO2 
*    NO2

      CALL rdo3xs(nw,wl,xso3,s226,s263,s298)
      CALL rdso2xs(nw,wl,xsso2)
      CALL rdno2xs(nw,wl,xsno2)
      endif

* ___ SECTION 4: SET MODEL ATMOSPHERE __________________________________

* temperature profile

      CALL settmp(nz,z,
     $     tlev,tlay,tempstem,nzstem)

*  air profile and Rayleigh optical depths

      pmbnew = -999.
      CALL setair(pmbnew,
     $     nz,z,nw,wl,
     $     airlev,dtrl,colinc,airstem,nzstem)

* Photo-chemical and photo-biological weigting functions. 
* For pchem, need to know temperature and pressure profiles.
* Output:
* from pbiol:  s(ks,kw) - for each weigting function label(ks)
* from pchem:  sj(kj,kz,kw) - for each reaction jlabel(kj)

      is = 0
c      CALL pbiol1(nw,wl,wc,is,sw,label)
      ns = is

      CALL pchem(nw,wl,nz,tlev,airlev,
     $     nj,sj,jlabel)

* ozone optical depths (must give temperature)

c      dobnew = 270.
      CALL setozo(dobnew,
     $     nz,z,nw,wl,
     $     xso3,s226,s263,s298,tlay,
     $     dto3,o3stem,nzstem)

* SO2 optical depth (also has temperature correction)
* so2new = new column SO2, in Dobson Units

      so2new =  0.
      CALL setso2(so2new,
     $     nz,z,nw,wl,
     $     xsso2, tlay,
     $     dtso2,so2stem,nzstem)

* NO2 optical depth (also has temperature correction)
* no2new = new column NO2, in Dobson Units

      no2new =  0.
      CALL setno2(no2new,
     $     nz,z,nw,wl,
     $     xsno2, tlay,
     $     dtno2,no2stem,nzstem)

*  cloud and aerosol optical depths:

      CALL setcld(nz,z,nw,wl,      ! clean sky
     $     dtcld,omcld,gcld)

      CALL setaer(nz,z,nw,wl,
     $     dtaer,omaer,gaer,mssaero,rhu,nzstem,aerod)

* surface albedo:

      CALL setalb(nw,wl,albedo)

* ___ SECTION 5: TIME AND LOCATION _____________________________________

* specify date and compute earth-sun distance correction

c      idate = 010307
      CALL sundis(idate,esfact)
c      WRITE(kout,*) 'idate = ', idate,' esfact = ', esfact
      do iw = 1, nw-1
         etf(iw) = f(iw) * esfact
      enddo

* specify latitude and longitude

c      alat = 22.3
c      along = 113.9
c      WRITE(kout,*)'lat = ',alat,' long = ',along

* below, can  chose between specific time (solar zenith angle is calculated) 
* or can set zenith angle to arbitrary value(s).  Loop DO 20 allows calculation 
* at multiple solar zenith angles  (or multiple times).

* Set starting time (ut = Universal Time, hrs.), and 
* time increment (dtime, in seconds)

c      ut0 = 0.
c      dtime = 3600./4.

* initalize time-integrated quantities

cC      call zero1(dose,ks)

* Loop over time (alternatively, can loop over solar zenith angle)

c      DO 20, itime = 1, 96
C         ut = ut0 + (dtime/3600.) * FLOAT(itime-1)

c         ut = 2.28889

* solar zenith angle calculation:

       CALL zenith(alat,along,idate,ut,azim,zen)
c       WRITE(kout,*) 'ut = ', ut, 'azimuth = ', azim, ' zen = ', zen

* ____ SECTION 6: CALCULATE ZENITH-ANGLE DEPENDENT QUANTITIES __________

* slant path lengths for spherical geometry

       CALL sphers(nz, z, zen, dsdh, nid)
       CALL airmas(nz, z, zen, dsdh, nid, colinc,
     $      vcol, scol)

* effective O2 optical depth (SR bands, must know zenith angle!)
* assign O2 cross section to sj(1,*,*)

       CALL seto2(nz,z,nw,wl,colinc,vcol,scol,dto2,xso2)
c       CALL sjo2(nz,nw,xso2,1,sj)  ! do not output O2 J-value

* ____ SECTION 7: WAVELENGTH LOOP ______________________________________

* initialize for wavelength integration

       CALL zero1(rate,ks)
       CALL zero2(valj,kj,kz)

** Main wavelength loop:

       DO 10, iw = 1, nw-1

** monochromatic radiative transfer:
*  outputs are  edir(iz), eup, fdir, fdn, fup

         CALL rtlink(nz,z,
     $        iw, albedo(iw), zen,
     $        dsdh,nid,
     $        dtrl,
     $        dto3,
     $        dto2,
     $        dtso2,
     $        dtno2,
     $        dtcld, omcld, gcld,
     $        dtaer,omaer,gaer,
     $        edir, edn, eup, fdir, fdn, fup)  ! edir, edn: direct beam, down-welling
                                               ! diffuse light 
        if(ccover(1).gt.cclean) then	    ! cloud takes effect
	  dtcld(1:nzstem,iw)=cldod(1:nzstem)   ! load cloud optical depth from STEM
	  CALL rtlink(nz,z,
     $        iw, albedo(iw), zen,
     $        dsdh,nid,
     $        dtrl,
     $        dto3,
     $        dto2,
     $        dtso2,
     $        dtno2,
     $        dtcld, omcld, gcld,
     $        dtaer,omaer,gaer,
     $        edir2, edn2, eup2, fdir2, fdn2, fup2)
          do iz=1,nzstem
c	   if(ccover(iz).gt.cblock) then  ! if blocked by cloud
c	    edir(iz)=edir2(iz)
c	    fdir(iz)=fdir2(iz)
c	   endif
c	   edn(iz)=(1-ccover(iz))*edn(iz)+ccover(iz)*edn2(iz)
c	   fdn(iz)=(1-ccover(iz))*fdn(iz)+ccover(iz)*fdn2(iz)
c	   fup(iz)=fup2(iz)
           edir(iz)=edir2(iz)
	   fdir(iz)=fdir2(iz)
	   edn(iz)=edn2(iz)
	   fdn(iz)=fdn2(iz)
	   fup(iz)=fup2(iz)
	  enddo
	endif   
	   
** surface irradiance and weighted radiation

         iz = 1 
         sirrad = etf(iw) * (edir(iz) + edn(iz))

         DO 15, is = 1, ns
            sprate = sirrad * sw(is,iw) 
            rate(is) = rate(is) + sprate * (wu(iw) - wl(iw))
 15      CONTINUE

** spherical irradiance (actinic flux)
* as a function of altitude
* convert to quanta s-1 nm-1 cm-2
* ( 1.e-4 * (wc*1e-9) / (hc = 6.62E-34 * 2.998E8) )

         DO 17 iz = 1, nz
            saflux = etf(iw)* 5.039e11 * wc(iw) *
     $           (fdir(iz) + fdn(iz) + fup(iz))

            DO 16, ij = 1, nj
               deltaj = saflux * sj(ij,iz,iw)
               valj(ij,iz) = valj(ij,iz) + deltaj * (wu(iw) - wl(iw))

	     if(abs(valj(ij,iz)).gt.1) then
	       print*, 'wrong J values iJ, iz, valj(iJ,iz), idate, ut=',
     1	          iJ,iz,valj(iJ,iz), idate, ut
               print*, 'saflux,sj(ij,iz,iw)=',saflux,sj(ij,iz,iw),
     2		etf(iw),wc(iw),fdir(iz),fdn(iz),fup(iz)
               print*,'dtaer = ', dtaer(1:nz,iw)
	       print*,'nz,z= ',nz,z(1:nz)
	       stop
	     endif 
 16         CONTINUE
            
 17      CONTINUE

 10   CONTINUE

*^^^^^^^^^^^^^^^^ end wavelength loop

 20   CONTINUE

** some examples of output:

* dose rates weighted by specific action spectra:

c      DO 35, is = 1, ns
c         WRITE(kout,99) is, label(is), rate(is)
c 35   CONTINUE

* photolysis rate coefficients (j-values) at surface

c      iz = 1
      do iz=1,15
c      write(kout,*)'Z=',iz
c      DO 36, ij = 1, nj
c         WRITE(kout,99) ij, jlabel(ij), valj(ij,iz)
c 36   CONTINUE
      enddo
 99   FORMAT(I4,1X,A40,1X,1PE10.3)

*^^^^^^^^^^^^^^^^ end zenith loop
      first=.false. 

*_______________________________________________________________________

      END


      SUBROUTINE addpnt ( x, y, ld, n, xnew, ynew )

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Add a point <xnew,ynew> to a set of data pairs <x,y>.  x must be in      =*
*=  ascending order                                                          =*
*-----------------------------------------------------------------------------*
*=  PARAMETERS:                                                              =*
*=  X    - REAL vector of length LD, x-coordinates                       (IO)=*
*=  Y    - REAL vector of length LD, y-values                            (IO)=*
*=  LD   - INTEGER, dimension of X, Y exactly as declared in the calling  (I)=*
*=         program                                                           =*
*=  N    - INTEGER, number of elements in X, Y.  On entry, it must be:   (IO)=*
*=         N < LD.  On exit, N is incremented by 1.                          =*
*=  XNEW - REAL, x-coordinate at which point is to be added               (I)=*
*=  YNEW - REAL, y-value of point to be added                             (I)=*
*-----------------------------------------------------------------------------*
*=  EDIT HISTORY:                                                            =*
*=  02/95  Original                                                          =*
*-----------------------------------------------------------------------------*
*= This program is free software;  you can redistribute it and/or modify     =*
*= it under the terms of the GNU General Public License as published by the  =*
*= Free Software Foundation;  either version 2 of the license, or (at your   =*
*= option) any later version.                                                =*
*= The TUV package is distributed in the hope that it will be useful, but    =*
*= WITHOUT ANY WARRANTY;  without even the implied warranty of MERCHANTIBI-  =*
*= LITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public     =*
*= License for more details.                                                 =*
*= To obtain a copy of the GNU General Public License, write to:             =*
*= Free Software Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.   =*
*-----------------------------------------------------------------------------*
*= To contact the authors, please mail to:                                   =*
*= Sasha Madronich, NCAR/ACD, P.O.Box 3000, Boulder, CO, 80307-3000, USA  or =*
*= send email to:  sasha@ucar.edu                                            =*
*-----------------------------------------------------------------------------*
*= Copyright (C) 1994,95,96  University Corporation for Atmospheric Research =*
*-----------------------------------------------------------------------------*

      IMPLICIT NONE

C calling parameters

      INTEGER ld, n
      REAL x(ld), y(ld)
      REAL xnew, ynew
      INTEGER ierr

C local variables

      INTEGER insert
      INTEGER i

C-----------------------------------------------------------------------

* initialize error flag

      ierr = 0

* check n<ld to make sure x will hold another point

      IF (n .GE. ld) THEN
         WRITE(0,*) '>>> ERROR (ADDPNT) <<<  Cannot expand array '
         WRITE(0,*) '                        All elements used.'
         STOP
      ENDIF

      insert = 1
      i = 2

* check, whether x is already sorted.
* also, use this loop to find the point at which xnew needs to be inserted
* into vector x, if x is sorted.

 10   CONTINUE
      IF (i .LT. n) THEN
        IF (x(i) .LT. x(i-1)) THEN
           WRITE(*,*) '>>> ERROR (ADDPNT) <<<  x-data must be '//
     >                'in ascending order!',i,n
           write(*,*)x
           STOP
        ELSE
           IF (xnew .GT. x(i)) insert = i + 1
        ENDIF
        i = i+1
        GOTO 10
      ENDIF

* if <xnew,ynew> needs to be appended at the end, just do so,
* otherwise, insert <xnew,ynew> at position INSERT

      IF ( xnew .GT. x(n) ) THEN
 
         x(n+1) = xnew
         y(n+1) = ynew
  
      ELSE

* shift all existing points one index up

         DO i = n, insert, -1
           x(i+1) = x(i)
           y(i+1) = y(i)
         ENDDO

* insert new point

         x(insert) = xnew
         y(insert) = ynew
  
      ENDIF

* increase total number of elements in x, y

      n = n+1

      END
      SUBROUTINE airmas(nz, z, zen, dsdh, nid, cz,
     $      vcol, scol)

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Calculate vertical and slant air columns, in spherical geometry, as a    =*
*=  function of altitude.                                                    =*
*-----------------------------------------------------------------------------*
*=  PARAMETERS:                                                              =*
*=  NZ      - INTEGER, number of specified altitude levels in the working (I)=*
*=            grid                                                           =*
*=  Z       - REAL, specified altitude working grid (km)                  (I)=*
*=  ZEN     - REAL, solar zenith angle (degrees)                          (I)=*
*=  DSDH    - REAL, slant path of direct beam through each layer crossed  (O)=*
*=            when travelling from the top of the atmosphere to layer i;     =*
*=            DSDH(i,j), i = 0..NZ-1, j = 1..NZ-1                            =*
*=  NID     - INTEGER, number of layers crossed by the direct beam when   (O)=*
*=            travelling from the top of the atmosphere to layer i;          =*
*=            NID(i), i = 0..NZ-1                                            =*
*=  VCOL    - REAL, output, vertical air column, molec cm-2, above level iz  =*
*=  SCOL    - REAL, output, slant air column in direction of sun, above iz   =*
*=            also in molec cm-2                                             =*
*-----------------------------------------------------------------------------*
*=  EDIT HISTORY:                                                            =*
*=  Original                                                                 =*
*-----------------------------------------------------------------------------*
*= This program is free software;  you can redistribute it and/or modify     =*
*= it under the terms of the GNU General Public License as published by the  =*
*= Free Software Foundation;  either version 2 of the license, or (at your   =*
*= option) any later version.                                                =*
*= The TUV package is distributed in the hope that it will be useful, but    =*
*= WITHOUT ANY WARRANTY;  without even the implied warranty of MERCHANTIBI-  =*
*= LITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public     =*
*= License for more details.                                                 =*
*= To obtain a copy of the GNU General Public License, write to:             =*
*= Free Software Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.   =*
*-----------------------------------------------------------------------------*
*= To contact the authors, please mail to:                                   =*
*= Sasha Madronich, NCAR/ACD, P.O.Box 3000, Boulder, CO, 80307-3000, USA  or =*
*= send email to:  sasha@ucar.edu                                            =*
*-----------------------------------------------------------------------------*
*= Copyright (C) 1994-2000   University Corporation for Atmospheric Research =*
*-----------------------------------------------------------------------------*

      IMPLICIT NONE
      INCLUDE 'tuv.params'

* Input:

      INTEGER nz
      REAL z(kz)
      REAL zen
      INTEGER nid(0:kz)
      REAL dsdh(0:kz,kz)
      REAL cz(kz)

* output: 

      REAL vcol(kz), scol(kz)

* internal:

      INTEGER id, j
      REAL sum, ssum, vsum

* calculate vertical and slant column from each level:
* work downward

      vsum = 0.
      ssum = 0.
      DO id = 0, nz - 1
         vsum = vsum + cz(nz-id)
         vcol(nz-id) = vsum
         sum = 0.
         IF(nid(id) .LT. 0) THEN
            sum = largest
         ELSE

* single pass layers:

            DO j = 1, MIN(nid(id), id)
               sum = sum + cz(nz-j)*dsdh(id,j)
c	       if(cz(nz-j).lt.0.or.dsdh(id,j).lt.0) then
c	        print*,'wrong cz or dsdh =',id,nz-j,cz(nz-j),dsdh(id,j)
c     1		  ,sum
c		stop
c	       endif	
            ENDDO

* double pass layers:

            DO j = MIN(nid(id),id)+1, nid(id)
               sum = sum + 2.*cz(nz-j)*dsdh(id,j)
            ENDDO

         ENDIF
         scol(nz - id) = sum
	 	 
c         if(id.gt.0.and.scol(nz-id+1).gt.sum) then
c	  print*,'scol in wrong sequence ',nz-id, sum, scol(nz-id+1),
c     1	    id, dsdh(id,j-1), dsdh(id,j),cz(nz-j),cz(nz-j+1)
c          print*,'nid=',nid
c	  stop
c	 endif   
      ENDDO
ctyh      
       DO id = nz-2, nz-1
        if(scol(id).le.scol(id+1)) 
     1   scol(id+1)=scol(id)-0.5*(scol(id-1)-scol(id))  ! correct the unknown error
       enddo	
      RETURN
      END




      FUNCTION fery(w)

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Calculate the action spectrum value for erythema at a given wavelength   =*
*=  according to: McKinlay, A.F and B.L.Diffey, A reference action spectrum  =*
*=  for ultraviolet induced erythema in human skin, CIE Journal, vol 6,      =*
*=  pp 17-22, 1987.                                                          =*
*=  Value at 300 nm = 0.6486                                                 =*
*-----------------------------------------------------------------------------*
*=  PARAMETERS:                                                              =*
*=  W - REAL, wavelength (nm)                                             (I)=*
*-----------------------------------------------------------------------------*
*=  EDIT HISTORY:                                                            =*
*=  Original                                                                 =*
*-----------------------------------------------------------------------------*
*= This program is free software;  you can redistribute it and/or modify     =*
*= it under the terms of the GNU General Public License as published by the  =*
*= Free Software Foundation;  either version 2 of the license, or (at your   =*
*= option) any later version.                                                =*
*= The TUV package is distributed in the hope that it will be useful, but    =*
*= WITHOUT ANY WARRANTY;  without even the implied warranty of MERCHANTIBI-  =*
*= LITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public     =*
*= License for more details.                                                 =*
*= To obtain a copy of the GNU General Public License, write to:             =*
*= Free Software Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.   =*
*-----------------------------------------------------------------------------*
*= To contact the authors, please mail to:                                   =*
*= Sasha Madronich, NCAR/ACD, P.O.Box 3000, Boulder, CO, 80307-3000, USA  or =*
*= send email to:  sasha@ucar.edu                                            =*
*-----------------------------------------------------------------------------*
*= Copyright (C) 1994,95,96  University Corporation for Atmospheric Research =*
*-----------------------------------------------------------------------------*

      IMPLICIT NONE

* input:
      REAL w 

* function value:
      REAL fery
*_______________________________________________________________________

      IF (w .LT. 250.) THEN
          fery = 1.
C outside the ery spectrum range
      ELSEIF ((w .GE. 250.) .AND. (w .LT. 298)) THEN
          fery = 1.
      ELSEIF ((w .GE. 298.) .AND. (w .LT. 328.)) THEN
          fery = 10.**( 0.094*(298.-w) )
      ELSEIF ((w .GE. 328.) .AND. (w .LT. 400.)) THEN
          fery = 10.**( 0.015*(139.-w) )
      ELSE
         fery = 1.E-36
C outside the ery spectrum range
      ENDIF

*_______________________________________________________________________

      RETURN
      END
      FUNCTION fsum(n,x)

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Compute the sum of the first N elements of a floating point vector.      =*
*-----------------------------------------------------------------------------*
*=  PARAMETERS:                                                              =*
*=  N  - INTEGER, number of elements to sum                               (I)=*
*=  X  - REAL, vector whose components are to be summed                   (I)=*
*-----------------------------------------------------------------------------*
*=  EDIT HISTORY:                                                            =*
*=  Original                                                                 =*
*-----------------------------------------------------------------------------*
*= This program is free software;  you can redistribute it and/or modify     =*
*= it under the terms of the GNU General Public License as published by the  =*
*= Free Software Foundation;  either version 2 of the license, or (at your   =*
*= option) any later version.                                                =*
*= The TUV package is distributed in the hope that it will be useful, but    =*
*= WITHOUT ANY WARRANTY;  without even the implied warranty of MERCHANTIBI-  =*
*= LITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public     =*
*= License for more details.                                                 =*
*= To obtain a copy of the GNU General Public License, write to:             =*
*= Free Software Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.   =*
*-----------------------------------------------------------------------------*
*= To contact the authors, please mail to:                                   =*
*= Sasha Madronich, NCAR/ACD, P.O.Box 3000, Boulder, CO, 80307-3000, USA  or =*
*= send email to:  sasha@ucar.edu                                            =*
*-----------------------------------------------------------------------------*
*= Copyright (C) 1994,95,96  University Corporation for Atmospheric Research =*
*-----------------------------------------------------------------------------*

      IMPLICIT NONE

* input:
      INTEGER n
      REAL x(n)

* function value:
      REAL fsum

* local:
      INTEGER i
*_______________________________________________________________________

      fsum = 0.
      DO 10, i = 1, n
         fsum=fsum+x(i)
   10 CONTINUE
*_______________________________________________________________________

      RETURN
      END
      FUNCTION futr(w)

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Calculate the action spectrum value for skin cancer of albino hairless   =*
*=  mice at a given wavelength according to:  deGRuijl, F.R., H.J.C.M.Steren-=*
*=  borg, P.D.Forbes, R.E.Davies, C.Colse, G.Kelfkens, H.vanWeelden,         =*
*=  and J.C.van der Leun, Wavelength dependence of skin cancer induction by  =*
*=  ultraviolet irradiation of albino hairless mice, Cancer Research, vol 53,=*
*=  pp. 53-60, 1993                                                          =*
*=  (Action spectrum for carcinomas)                                         =*
*-----------------------------------------------------------------------------*
*=  PARAMETERS:                                                              =*
*=  W  - REAL, wavelength (nm)                                            (I)=*
*-----------------------------------------------------------------------------*
*=  EDIT HISTORY:                                                            =*
*=  Original                                                                 =*
*-----------------------------------------------------------------------------*
*= This program is free software;  you can redistribute it and/or modify     =*
*= it under the terms of the GNU General Public License as published by the  =*
*= Free Software Foundation;  either version 2 of the license, or (at your   =*
*= option) any later version.                                                =*
*= The TUV package is distributed in the hope that it will be useful, but    =*
*= WITHOUT ANY WARRANTY;  without even the implied warranty of MERCHANTIBI-  =*
*= LITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public     =*
*= License for more details.                                                 =*
*= To obtain a copy of the GNU General Public License, write to:             =*
*= Free Software Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.   =*
*-----------------------------------------------------------------------------*
*= To contact the authors, please mail to:                                   =*
*= Sasha Madronich, NCAR/ACD, P.O.Box 3000, Boulder, CO, 80307-3000, USA  or =*
*= send email to:  sasha@ucar.edu                                            =*
*-----------------------------------------------------------------------------*
*= Copyright (C) 1994,95,96  University Corporation for Atmospheric Research =*
*-----------------------------------------------------------------------------*

      IMPLICIT NONE

* input:
      REAL w

* function value:
      REAL futr

* local:
      REAL a1, a2, a3, a4, a5,
     >     x1, x2, x3, x4, x5,
     >     t1, t2, t3, t4, t5,
     >     b1, b2, b3, b4, b5,
     >     p
*_______________________________________________________________________

      a1 = -10.91
      a2 = - 0.86
      a3 = - 8.60
      a4 = - 9.36
      a5 = -13.15

      x1 = 270.
      x2 = 302.
      x3 = 334.
      x4 = 367.
      x5 = 400.

      t1 = (w-x2)*(w-x3)*(w-x4)*(w-x5)
      t2 = (w-x1)*(w-x3)*(w-x4)*(w-x5)
      t3 = (w-x1)*(w-x2)*(w-x4)*(w-x5)
      t4 = (w-x1)*(w-x2)*(w-x3)*(w-x5)
      t5 = (w-x1)*(w-x2)*(w-x3)*(w-x4)

      b1 = (x1-x2)*(x1-x3)*(x1-x4)*(x1-x5)
      b2 = (x2-x1)*(x2-x3)*(x2-x4)*(x2-x5)
      b3 = (x3-x1)*(x3-x2)*(x3-x4)*(x3-x5)
      b4 = (x4-x1)*(x4-x2)*(x4-x3)*(x4-x5)
      b5 = (x5-x1)*(x5-x2)*(x5-x3)*(x5-x4)

      p = a1*t1/b1 + a2*t2/b2 + a3*t3/b3 + a4*t4/b4 + a5*t5/b5

      futr  = EXP(p)
*_______________________________________________________________________

      RETURN
      END
      SUBROUTINE gridck(k,n,x,ok)

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Check a grid X for various improperties.  The values in X have to comply =*
*=  with the following rules:                                                =*
*=  1) Number of actual points cannot exceed declared length of X            =*
*=  2) Number of actual points has to be greater than or equal to 2          =*
*=  3) X-values must be non-negative                                         =*
*=  4) X-values must be unique                                               =*
*=  5) X-values must be in ascending order                                   =*
*-----------------------------------------------------------------------------*
*=  PARAMETERS:                                                              =*
*=  K  - INTEGER, length of X as declared in the calling program          (I)=*
*=  N  - INTEGER, number of actual points in X                            (I)=*
*=  X  - REAL, vector (grid) to be checked                                (I)=*
*=  OK - LOGICAL, .TRUE. -> X agrees with rules 1)-5)                     (O)=*
*=                .FALSE.-> X violates at least one of 1)-5)                 =*
*-----------------------------------------------------------------------------*
*=  EDIT HISTORY:                                                            =*
*=  Original                                                                 =*
*-----------------------------------------------------------------------------*
*= This program is free software;  you can redistribute it and/or modify     =*
*= it under the terms of the GNU General Public License as published by the  =*
*= Free Software Foundation;  either version 2 of the license, or (at your   =*
*= option) any later version.                                                =*
*= The TUV package is distributed in the hope that it will be useful, but    =*
*= WITHOUT ANY WARRANTY;  without even the implied warranty of MERCHANTIBI-  =*
*= LITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public     =*
*= License for more details.                                                 =*
*= To obtain a copy of the GNU General Public License, write to:             =*
*= Free Software Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.   =*
*-----------------------------------------------------------------------------*
*= To contact the authors, please mail to:                                   =*
*= Sasha Madronich, NCAR/ACD, P.O.Box 3000, Boulder, CO, 80307-3000, USA  or =*
*= send email to:  sasha@ucar.edu                                            =*
*-----------------------------------------------------------------------------*
*= Copyright (C) 1994,95,96  University Corporation for Atmospheric Research =*
*-----------------------------------------------------------------------------*

      IMPLICIT NONE

      INCLUDE 'tuv.params'

* input:
      INTEGER k, n
      REAL x(k)

* output:
      LOGICAL ok

* local:
      INTEGER i
*_______________________________________________________________________

      ok = .TRUE.

* check if dimension meaningful and within bounds

      IF (n .GT. k) THEN
         ok = .false.
c         WRITE(kout,100)
         RETURN
      ENDIF         
  100 FORMAT('Number of data exceeds dimension')

      IF (n .LT. 2) THEN
         ok = .FALSE.
c         WRITE(kout,101)
         RETURN
      ENDIF
  101 FORMAT('Too few data, number of data points must be >= 2')

* disallow negative grid values

      IF(x(1) .LT. 0.) THEN
         ok = .FALSE.
c         WRITE(kout,105)
         RETURN
      ENDIF
  105 FORMAT('Grid cannot start below zero')

* check sorting

      DO 10, i = 2, n
         IF( x(i) .LE. x(i-1)) THEN
            ok = .FALSE.
c            WRITE(kout,110)
            RETURN
         ENDIF
   10 CONTINUE
  110 FORMAT('Grid is not sorted or contains multiple values')
*_______________________________________________________________________

      RETURN
      END
      SUBROUTINE gridw(nw,wl,wc,wu)

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Create the altitude grid for all interpolations and radiative transfer   =*
*=  calculations.  Grid may be irregularly spaced.  Wavelengths are in nm.   =*
*=  No gaps are allowed within the wavelength grid.                          =*
*-----------------------------------------------------------------------------*
*=  PARAMETERS:                                                              =*
*=  NW  - INTEGER, number of wavelength grid _points_                     (O)=*
*=  WL  - REAL, vector carrying the lower limit of each wavel. interval   (O)=*
*=  WC  - REAL, vector carrying the center wavel of each wavel. interval  (O)=*
*=              (wc(i) = 0.5*(wl(i)+wu(i), i = 1..NW-1)                      =*
*=  WU  - REAL, vector carrying the upper limit of each wavel. interval   (O)=*
*-----------------------------------------------------------------------------*
*=  EDIT HISTORY:                                                            =*
*=  Original                                                                 =*
*-----------------------------------------------------------------------------*
*= This program is free software;  you can redistribute it and/or modify     =*
*= it under the terms of the GNU General Public License as published by the  =*
*= Free Software Foundation;  either version 2 of the license, or (at your   =*
*= option) any later version.                                                =*
*= The TUV package is distributed in the hope that it will be useful, but    =*
*= WITHOUT ANY WARRANTY;  without even the implied warranty of MERCHANTIBI-  =*
*= LITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public     =*
*= License for more details.                                                 =*
*= To obtain a copy of the GNU General Public License, write to:             =*
*= Free Software Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.   =*
*-----------------------------------------------------------------------------*
*= To contact the authors, please mail to:                                   =*
*= Sasha Madronich, NCAR/ACD, P.O.Box 3000, Boulder, CO, 80307-3000, USA  or =*
*= send email to:  sasha@ucar.edu                                            =*
*-----------------------------------------------------------------------------*
*= Copyright (C) 1994,95,96  University Corporation for Atmospheric Research =*
*-----------------------------------------------------------------------------*

      IMPLICIT NONE

      INCLUDE 'tuv.params'

* output:

      REAL wl(kw), wc(kw), wu(kw)
      INTEGER nw

* local:
      integer nkock,nisak
      parameter(nkock=16,nisak=120)
      real wl_kock(nkock),wu_kock(nkock),wc_isak(nisak),wl_isak(nisak),
     1  wu_isak(nisak)

      data wl_kock(1:nkock)/175.439,176.991,178.571,180.180,181.818,
     1 183.486,185.185, 186.916,188.679,190.476,192.308,194.175,
     2 196.078,198.020,200.000,202.020/
       data wu_kock(1:nkock)/176.991,178.571,180.180,181.818,183.486,
     1  185.185,186.916,188.679,190.476,192.308,194.175,196.078,198.020,
     2  200.000, 202.020,204.082/
       data wc_isak(1:nisak)/205.128,207.254,209.424,211.640,213.904,
     1 216.216,218.579,220.994,223.464,225.989,228.571,231.214,233.918,
     2 236.686,239.521,242.424,245.399,248.447,251.572,254.777,258.065,
     3 261.438,264.901,268.456,272.109,275.862,279.720,283.688,
     4 287.770,291.971,296.296,300.500,303.000,304.000,305.000,
     5 306.000,307.000,308.000,309.000,310.000,311.000,312.000,
     6 313.000,314.000,316.000,320.000,325.000,330.000,335.000,
     7 340.000,345.000,350.000,355.000,360.000,365.000,370.000,
     8 375.000,380.000,385.000,390.000,395.000,400.000,405.000,
     9 410.000,415.000,420.000,425.000,430.000,435.000,440.000,
     a 445.000,450.000,455.000,460.000,465.000,470.000,475.000,
     b 480.000,485.000,490.000,495.000,500.000,505.000,510.000,
     c 515.000,520.000,525.000,530.000,535.000,540.000,545.000,
     d 550.000,555.000,560.000,565.000,570.000,575.000,580.000,
     e 585.000,590.000,595.000,600.000,605.000,610.000,615.000,
     f 620.000,625.000,630.000,635.000,640.000,644.800,651.050,
     g 660.000,670.000,680.000,690.000,700.000,710.000,720.000,
     h 730.000/
       data wl_isak(1:nisak)/204.082,206.186,208.333,210.526,212.766,
     1 215.054,217.391,219.780,222.222,224.719,227.273,229.885,232.558,
     2 235.294,238.095,240.964,243.902,246.914,250.000,253.165,256.410,
     3 259.740,263.158,266.667,270.270,273.973,277.778,281.690,
     4 285.714,289.855,294.118,298.500,302.500,303.500,304.500,
     5 305.500,306.500,307.500,308.500,309.500,310.500,311.500,
     6 312.500,313.500,314.500,317.500,322.500,327.500,332.500,
     7 337.500,342.500,347.500,352.500,357.500,362.500,367.500,
     8 372.500,377.500,382.500,387.500,392.500,397.500,402.500,
     9 407.500,412.500,417.500,422.500,427.500,432.500,437.500,
     a 442.500,447.500,452.500,457.500,462.500,467.500,472.500,
     b 477.500,482.500,487.500,492.500,497.500,502.500,507.500,
     c 512.500,517.500,522.500,527.500,532.500,537.500,542.500,
     d 547.500,552.500,557.500,562.500,567.500,572.500,577.500,
     e 582.500,587.500,592.500,597.500,602.500,607.500,612.500,
     f 617.500,622.500,627.500,632.500,637.500,642.500,647.100,
     g 655.000,665.000,675.000,685.000,695.000,705.000,715.000,
     h 725.000/
       data wu_isak(1:nisak)/206.186,208.333,210.526,212.766,215.054,
     1 217.391,219.780,222.222,224.719,227.273,229.885,232.558,235.294,
     2 238.095,240.964,243.902,246.914,250.000,253.165,256.410,259.740,
     3 263.158,266.667,270.270,273.973,277.778,281.690,285.714,
     4 289.855,294.118,298.500,302.500,303.500,304.500,305.500,
     5 306.500,307.500,308.500,309.500,310.500,311.500,312.500,
     6 313.500,314.500,317.500,322.500,327.500,332.500,337.500,
     7 342.500,347.500,352.500,357.500,362.500,367.500,372.500,
     8 377.500,382.500,387.500,392.500,397.500,402.500,407.500,
     9 412.500,417.500,422.500,427.500,432.500,437.500,442.500,
     a 447.500,452.500,457.500,462.500,467.500,472.500,477.500,
     b 482.500,487.500,492.500,497.500,502.500,507.500,512.500,
     c 517.500,522.500,527.500,532.500,537.500,542.500,547.500,
     d 552.500,557.500,562.500,567.500,572.500,577.500,582.500,
     e 587.500,592.500,597.500,602.500,607.500,612.500,617.500,
     f 622.500,627.500,632.500,637.500,642.500,647.100,655.000,
     g 665.000,675.000,685.000,695.000,705.000,715.000,725.000,
     h 735.000/

      REAL wincr
      INTEGER iw
      LOGICAL ok
      INTEGER idum
      REAL dum
      INTEGER mopt

*_______________________________________________________________________

**** chose wavelengths

* some pre-set options
*     mopt = 1    equal spacing
*     mopt = 2    Isaksen's grid
*     mopt = 3    combined Kockarts/Isaksen grid + Lyman-Alpha
*     mopt = 4    user-defined

      mopt = 3
      IF (mopt .EQ. 1) GO TO 1
      IF (mopt .EQ. 2) GO TO 2
      IF (mopt .EQ. 3) GO TO 3
      IF (mopt .EQ. 4) GO TO 4

 1    CONTINUE
      nw = 140 + 1
      wincr = 1.0
      DO 10, iw = 1, nw-1
         wl(iw) = 280. + wincr*FLOAT(iw-1)
         wu(iw) = wl(iw) + wincr
         wc(iw) = ( wl(iw) + wu(iw) )/2.
   10 CONTINUE
      wl(nw) = wu(nw-1)
      GO TO 9

 2    CONTINUE
      nw = 0
      OPEN(unit=kin,file='DATAE1/GRIDS/isaksen.grid',status='old')
      DO iw = 1, 2 
         READ(kin,*)
      ENDDO
      DO iw = 1, 130
         nw = nw + 1
         READ(kin,*) idum, dum, wc(nw), wl(nw), wu(nw)
      ENDDO
      CLOSE(kin)
      nw = nw + 1
      wl(nw) = wu(nw-1)
      GO TO 9

*** grid for strat photolysis calculations, extended at short wavelengths

 3    CONTINUE
      nw = 1
* include Lyman-Alpha wavelengths ([120.,121.4],[121.4,121.9],[123.,-])
      wl(nw) = 120.0
      wu(nw) = 121.4
      wc(nw) = (wl(nw)+wu(nw))*0.5
      nw = nw+1
      wl(nw) = wu(nw-1)
      wu(nw) = 121.9
      wc(nw) = (wl(nw)+wu(nw))*0.5
      nw = nw+1
      wl(nw) = wu(nw-1)
      wu(nw) = 123.0
      wc(nw) = (wl(nw)+wu(nw))*0.5
      nw = nw+1
      wl(nw) = wu(nw-1)

      DO iw = 1, nkock
         nw = nw + 1
         wl(nw)=wl_kock(iw)
	 wu(nw)=wu_kock(iw)
         wc(nw) = ( wl(nw) + wu(nw) ) / 2.
         IF (iw	.eq. 1) THEN
             wu(nw-1) = wl(nw)
             wc(nw-1) = (wl(nw-1) + wu(nw-1))*0.5
         ENDIF
      ENDDO

      DO iw = 1, nisak
         nw = nw + 1
         wl(nw)=wl_isak(iw)
	 wu(nw)=wu_isak(iw)
         wc(nw) = wc_isak(iw)
      ENDDO
      nw = nw + 1
      wl(nw) = wu(nw-1)
      GO TO 9

 4    CONTINUE
* define wavelength intervals of width 1 nm from 150 - 198 nm:
      nw = 1
      wl(1) = 150.
      DO iw = 151, 198
        wu(nw) = Float(iw)
        wc(nw) = (wl(nw) + wu(nw))/2.
        nw = nw+1
        wl(nw) = Float(iw)
      ENDDO
c      write(kout,*), 'nw,wu,wc,wl=',nw-1,wu(nw-1),wc(nw-1),wl(nw-1)
* define wavelength intervals of width 0.5 nm from 198 - 494 nm:

      do while(wu(nw-1).le.494)
        wu(nw) = wu(nw-1)+0.5
        wc(nw) = (wl(nw) + wu(nw))/2.
        nw = nw+1
        wl(nw) = wl(nw-1)+0.5
      ENDDO
c      write(kout,*), 'nw,wu,wc,wl=',nw-1,wu(nw-1),wc(nw-1),wl(nw-1)      
* define wavelength intervals of width 1 nm from 494 - 900 nm:

      do while(wu(nw-1).le.900)
        wu(nw) = wu(nw-1)+ 1.
        wc(nw) = (wl(nw) + wu(nw))/2.
        nw = nw+1
        wl(nw) = wl(nw-1)+1.
      ENDDO
c      write(kout,*), 'nw,wu,wc,wl=',nw-1,wu(nw-1),wc(nw-1),wl(nw-1)

 9    CONTINUE

***
* write to record

c      WRITE(kout,*)'w-grid:',nw,wl(1),wl(nw)

* check grid for assorted improprieties:

      CALL gridck(kw,nw,wl,ok)

      IF (.NOT. ok) THEN
c         WRITE(kout,*)'STOP in GRIDW:  The w-grid does not make sense'
         STOP
      ENDIF

*_______________________________________________________________________

      RETURN
      END
      SUBROUTINE gridz(nz,z)

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Create the altitude grid for all interpolations and radiative transfer   =*
*=  calculations.  Grid may be irregularly spaced.  All altitudes are in     =*
*=  kilometers (km).  The altitude at index 1 specifies the surface elevation=*
*-----------------------------------------------------------------------------*
*=  PARAMETERS:                                                              =*
*=  NZ  - INTEGER, number of altitude points (levels)                     (O)=*
*=  Z   - REAL, vector of altitude levels (in km)                         (O)=*
*-----------------------------------------------------------------------------*
*=  EDIT HISTORY:                                                            =*
*=  Original                                                                 =*
*-----------------------------------------------------------------------------*
*= This program is free software;  you can redistribute it and/or modify     =*
*= it under the terms of the GNU General Public License as published by the  =*
*= Free Software Foundation;  either version 2 of the license, or (at your   =*
*= option) any later version.                                                =*
*= The TUV package is distributed in the hope that it will be useful, but    =*
*= WITHOUT ANY WARRANTY;  without even the implied warranty of MERCHANTIBI-  =*
*= LITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public     =*
*= License for more details.                                                 =*
*= To obtain a copy of the GNU General Public License, write to:             =*
*= Free Software Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.   =*
*-----------------------------------------------------------------------------*
*= To contact the authors, please mail to:                                   =*
*= Sasha Madronich, NCAR/ACD, P.O.Box 3000, Boulder, CO, 80307-3000, USA  or =*
*= send email to:  sasha@ucar.edu                                            =*
*-----------------------------------------------------------------------------*
*= Copyright (C) 1994,95,96  University Corporation for Atmospheric Research =*
*-----------------------------------------------------------------------------*

      IMPLICIT NONE

      INCLUDE 'tuv.params'

* output: altitude working grid:

      REAL z(kz)
      INTEGER nz

* local:

      REAL zincr
      INTEGER i
      LOGICAL ok
*_______________________________________________________________________

* set vertical grid of the atmosphere.  All values should be in km.
* User specifies upright grid (surface at lowest km value, increasing
* upwards:
*     -  NZ = total number of user levels
*     -  Z(I) = altitude in km for each level.
* Note "levels" are vertical points
*      "layers" are vertical distances between levels

* set atmospheric level altitudes (in real km), including 
* top-most level.
* non-uniform spacing is possible 
* z(1) is the elevation of the surface (km asl), and can be specified either
* here or in the main progarm.

      nz = 80 + 1
      zincr = 1.
      DO 10, i = 2, nz
         z(i) =  z(1) + FLOAT(i-1) * zincr
   10 CONTINUE

* write to record:

c      WRITE(kout,*)'z-grid:',nz,z(1),z(nz)

* check grid for assorted improprieties:

      CALL gridck(kz,nz,z,ok)

      IF (.NOT. ok) THEN
c         WRITE(kout,*)'STOP in GRIDZ:  The z-grid does not make sense'
         STOP
      ENDIF
*_______________________________________________________________________

      RETURN
      END
      SUBROUTINE inter1(ng,xg,yg, n,x,y)

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Map input data given on single, discrete points, onto a discrete target  =*
*=  grid.                                                                    =*
*=  The original input data are given on single, discrete points of an       =*
*=  arbitrary grid and are being linearly interpolated onto a specified      =*
*=  discrete target grid.  A typical example would be the re-gridding of a   =*
*=  given data set for the vertical temperature profile to match the speci-  =*
*=  fied altitude grid.                                                      =*
*=  Some caution should be used near the end points of the grids.  If the    =*
*=  input data set does not span the range of the target grid, the remaining =*
*=  points will be set to zero, as extrapolation is not permitted.           =*
*=  If the input data does not encompass the target grid, use ADDPNT to      =*
*=  expand the input array.                                                  =*
*-----------------------------------------------------------------------------*
*=  PARAMETERS:                                                              =*
*=  NG  - INTEGER, number of points in the target grid                    (I)=*
*=  XG  - REAL, target grid (e.g. altitude grid)                          (I)=*
*=  YG  - REAL, y-data re-gridded onto XG                                 (O)=*
*=  N   - INTEGER, number of points in the input data set                 (I)=*
*=  X   - REAL, grid on which input data are defined                      (I)=*
*=  Y   - REAL, input y-data                                              (I)=*
*-----------------------------------------------------------------------------*
*=  EDIT HISTORY:                                                            =*
*=  01/95  Loop 10 restructured                                              =*
*-----------------------------------------------------------------------------*
*= This program is free software;  you can redistribute it and/or modify     =*
*= it under the terms of the GNU General Public License as published by the  =*
*= Free Software Foundation;  either version 2 of the license, or (at your   =*
*= option) any later version.                                                =*
*= The TUV package is distributed in the hope that it will be useful, but    =*
*= WITHOUT ANY WARRANTY;  without even the implied warranty of MERCHANTIBI-  =*
*= LITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public     =*
*= License for more details.                                                 =*
*= To obtain a copy of the GNU General Public License, write to:             =*
*= Free Software Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.   =*
*-----------------------------------------------------------------------------*
*= To contact the authors, please mail to:                                   =*
*= Sasha Madronich, NCAR/ACD, P.O.Box 3000, Boulder, CO, 80307-3000, USA  or =*
*= send email to:  sasha@ucar.edu                                            =*
*-----------------------------------------------------------------------------*
*= Copyright (C) 1994,95,96  University Corporation for Atmospheric Research =*
*-----------------------------------------------------------------------------*

      IMPLICIT NONE

* input:
      INTEGER n, ng
      REAL xg(ng)
      REAL x(n), y(n)

* output:
      REAL yg(ng)

* local:
      REAL slope
      INTEGER jsave, i, j
*_______________________________________________________________________

      jsave = 1
      DO 20, i = 1, ng
         yg(i) = 0.
         j = jsave
   10    CONTINUE
            IF ((x(j) .GT. xg(i)) .OR. (xg(i) .GE. x(j+1))) THEN
               j = j+1
               IF (j .LE. n-1) GOTO 10
*        ---- end of loop 10 ----
            ELSE
               slope = (y(j+1)-y(j)) / (x(j+1)-x(j))
               yg(i) = y(j) + slope * (xg(i) - x(j))
               jsave = j
             ENDIF
   20 CONTINUE
*_______________________________________________________________________

      RETURN
      END
      SUBROUTINE inter2(ng,xg,yg,n,x,y,ierr)

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Map input data given on single, discrete points onto a set of target     =*
*=  bins.                                                                    =*
*=  The original input data are given on single, discrete points of an       =*
*=  arbitrary grid and are being linearly interpolated onto a specified set  =*
*=  of target bins.  In general, this is the case for most of the weighting  =*
*=  functions (action spectra, molecular cross section, and quantum yield    =*
*=  data), which have to be matched onto the specified wavelength intervals. =*
*=  The average value in each target bin is found by averaging the trapezoi- =*
*=  dal area underneath the input data curve (constructed by linearly connec-=*
*=  ting the discrete input values).                                         =*
*=  Some caution should be used near the endpoints of the grids.  If the     =*
*=  input data set does not span the range of the target grid, an error      =*
*=  message is printed and the execution is stopped, as extrapolation of the =*
*=  data is not permitted.                                                   =*
*=  If the input data does not encompass the target grid, use ADDPNT to      =*
*=  expand the input array.                                                  =*
*-----------------------------------------------------------------------------*
*=  PARAMETERS:                                                              =*
*=  NG  - INTEGER, number of bins + 1 in the target grid                  (I)=*
*=  XG  - REAL, target grid (e.g., wavelength grid);  bin i is defined    (I)=*
*=        as [XG(i),XG(i+1)] (i = 1..NG-1)                                   =*
*=  YG  - REAL, y-data re-gridded onto XG, YG(i) specifies the value for  (O)=*
*=        bin i (i = 1..NG-1)                                                =*
*=  N   - INTEGER, number of points in input grid                         (I)=*
*=  X   - REAL, grid on which input data are defined                      (I)=*
*=  Y   - REAL, input y-data                                              (I)=*
*-----------------------------------------------------------------------------*
*=  EDIT HISTORY:                                                            =*
*=  01/95  Major restucturing of the entire subroutine                       =*
*-----------------------------------------------------------------------------*
*= This program is free software;  you can redistribute it and/or modify     =*
*= it under the terms of the GNU General Public License as published by the  =*
*= Free Software Foundation;  either version 2 of the license, or (at your   =*
*= option) any later version.                                                =*
*= The TUV package is distributed in the hope that it will be useful, but    =*
*= WITHOUT ANY WARRANTY;  without even the implied warranty of MERCHANTIBI-  =*
*= LITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public     =*
*= License for more details.                                                 =*
*= To obtain a copy of the GNU General Public License, write to:             =*
*= Free Software Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.   =*
*-----------------------------------------------------------------------------*
*= To contact the authors, please mail to:                                   =*
*= Sasha Madronich, NCAR/ACD, P.O.Box 3000, Boulder, CO, 80307-3000, USA  or =*
*= send email to:  sasha@ucar.edu                                            =*
*-----------------------------------------------------------------------------*
*= Copyright (C) 1994,95,96  University Corporation for Atmospheric Research =*
*-----------------------------------------------------------------------------*

      IMPLICIT NONE

* input:
      INTEGER ng, n
      REAL x(n), y(n), xg(ng)

* output:
      REAL yg(ng)

* local:
      REAL area, xgl, xgu
      REAL darea, slope
      REAL a1, a2, b1, b2
      INTEGER ngintv
      INTEGER i, k, jstart
      INTEGER ierr
*_______________________________________________________________________

      ierr = 0

*  test for correct ordering of data, by increasing value of x

      DO 10, i = 2, n
         IF (x(i) .LE. x(i-1)) THEN
            ierr = 1
            WRITE(*,*)'data not sorted',i
	    write(*,*)x
            RETURN
         ENDIF
   10 CONTINUE     

      DO i = 2, ng
        IF (xg(i) .LE. xg(i-1)) THEN
           ierr = 2
          WRITE(0,*) '>>> ERROR (inter2) <<<  xg-grid not sorted!'
          RETURN
        ENDIF
      ENDDO

* check for xg-values outside the x-range

      IF ( (x(1) .GT. xg(1)) .OR. (x(n) .LT. xg(ng)) ) THEN
          WRITE(0,*) '>>> ERROR (inter2) <<<  Data do not span '//
     >               'grid.  '
          WRITE(0,*) 'x(1), xg(1), x(n), xg(ng)= ',x(1), xg(1), x(n),
     1	   xg(ng)

          WRITE(0,*) '                        Use ADDPNT to '//
     >               'expand data and re-run.'
          STOP
      ENDIF

*  find the integral of each grid interval and use this to 
*  calculate the average y value for the interval      
*  xgl and xgu are the lower and upper limits of the grid interval

      jstart = 1
      ngintv = ng - 1
      DO 50, i = 1,ngintv

* initalize:

            area = 0.0
            xgl = xg(i)
            xgu = xg(i+1)

*  discard data before the first grid interval and after the 
*  last grid interval
*  for internal grid intervals, start calculating area by interpolating
*  between the last point which lies in the previous interval and the
*  first point inside the current interval

            k = jstart
            IF (k .LE. n-1) THEN

*  if both points are before the first grid, go to the next point
   30         CONTINUE
                IF (x(k+1) .LE. xgl) THEN
                   jstart = k - 1
                   k = k+1
                   IF (k .LE. n-1) GO TO 30
                ENDIF


*  if the last point is beyond the end of the grid, complete and go to the next
*  grid
   40         CONTINUE
                 IF ((k .LE. n-1) .AND. (x(k) .LT. xgu)) THEN          

                    jstart = k-1

* compute x-coordinates of increment

                    a1 = MAX(x(k),xgl)
                    a2 = MIN(x(k+1),xgu)

*  if points coincide, contribution is zero

                    IF (x(k+1).EQ.x(k)) THEN
                       darea = 0.e0
                    ELSE
                       slope = (y(k+1) - y(k))/(x(k+1) - x(k))
                       b1 = y(k) + slope*(a1 - x(k))
                       b2 = y(k) + slope*(a2 - x(k))
                       darea = (a2 - a1)*(b2 + b1)/2.
                    ENDIF


*  find the area under the trapezoid from a1 to a2

                    area = area + darea

* go to next point
              
                    k = k+1
                    GO TO 40

                ENDIF

            ENDIF

*  calculate the average y after summing the areas in the interval
            yg(i) = area/(xgu - xgl)

   50 CONTINUE
*_______________________________________________________________________

      RETURN
      END

      SUBROUTINE inter3(ng,xg,yg, n,x,y, FoldIn)

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Map input data given on a set of bins onto a different set of target     =*
*=  bins.                                                                    =*
*=  The input data are given on a set of bins (representing the integral     =*
*=  of the input quantity over the range of each bin) and are being matched  =*
*=  onto another set of bins (target grid).  A typical example would be an   =*
*=  input data set spcifying the extra-terrestrial flux on wavelength inter- =*
*=  vals, that has to be matched onto the working wavelength grid.           =*
*=  The resulting area in a given bin of the target grid is calculated by    =*
*=  simply adding all fractional areas of the input data that cover that     =*
*=  particular target bin.                                                   =*
*=  Some caution should be used near the endpoints of the grids.  If the     =*
*=  input data do not span the full range of the target grid, the area in    =*
*=  the "missing" bins will be assumed to be zero.  If the input data extend =*
*=  beyond the upper limit of the target grid, the user has the option to    =*
*=  integrate the "overhang" data and fold the remaining area back into the  =*
*=  last target bin.  Using this option is recommended when re-gridding      =*
*=  vertical profiles that directly affect the total optical depth of the    =*
*=  model atmosphere.                                                        =*
*-----------------------------------------------------------------------------*
*=  PARAMETERS:                                                              =*
*=  NG     - INTEGER, number of bins + 1 in the target grid               (I)=*
*=  XG     - REAL, target grid (e.g. working wavelength grid);  bin i     (I)=*
*=           is defined as [XG(i),XG(i+1)] (i = 1..NG-1)                     =*
*=  YG     - REAL, y-data re-gridded onto XG;  YG(i) specifies the        (O)=*
*=           y-value for bin i (i = 1..NG-1)                                 =*
*=  N      - INTEGER, number of bins + 1 in the input grid                (I)=*
*=  X      - REAL, input grid (e.g. data wavelength grid);  bin i is      (I)=*
*=           defined as [X(i),X(i+1)] (i = 1..N-1)                           =*
*=  Y      - REAL, input y-data on grid X;  Y(i) specifies the            (I)=*
*=           y-value for bin i (i = 1..N-1)                                  =*
*=  FoldIn - Switch for folding option of "overhang" data                 (I)=*
*=           FoldIn = 0 -> No folding of "overhang" data                     =*
*=           FoldIn = 1 -> Integerate "overhang" data and fold back into     =*
*=                         last target bin                                   =*
*-----------------------------------------------------------------------------*
*=  EDIT HISTORY:                                                            =*
*=  06/96  Added FoldIn switch                                               =*
*-----------------------------------------------------------------------------*
*= This program is free software;  you can redistribute it and/or modify     =*
*= it under the terms of the GNU General Public License as published by the  =*
*= Free Software Foundation;  either version 2 of the license, or (at your   =*
*= option) any later version.                                                =*
*= The TUV package is distributed in the hope that it will be useful, but    =*
*= WITHOUT ANY WARRANTY;  without even the implied warranty of MERCHANTIBI-  =*
*= LITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public     =*
*= License for more details.                                                 =*
*= To obtain a copy of the GNU General Public License, write to:             =*
*= Free Software Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.   =*
*-----------------------------------------------------------------------------*
*= To contact the authors, please mail to:                                   =*
*= Sasha Madronich, NCAR/ACD, P.O.Box 3000, Boulder, CO, 80307-3000, USA  or =*
*= send email to:  sasha@ucar.edu                                            =*
*-----------------------------------------------------------------------------*
*= Copyright (C) 1994,95,96  University Corporation for Atmospheric Research =*
*-----------------------------------------------------------------------------*

      IMPLICIT NONE
      
* input:
      INTEGER n, ng
      REAL xg(ng)
      REAL x(n), y(n)

      INTEGER FoldIn

* output:
      REAL yg(ng)

* local:
      REAL a1, a2, sum
      REAL tail
      INTEGER jstart, i, j, k
*_______________________________________________________________________

* check whether flag given is legal
      IF ((FoldIn .NE. 0) .AND. (FoldIn .NE. 1)) THEN
         WRITE(0,*) '>>> ERROR (inter3) <<<  Value for FOLDIN invalid. '
         WRITE(0,*) '                        Must be 0 or 1'
         STOP
      ENDIF

* do interpolation

      jstart = 1

      DO 30, i = 1, ng - 1

         yg(i) = 0.
         sum = 0.
         j = jstart

         IF (j .LE. n-1) THEN

   20      CONTINUE

             IF (x(j+1) .LT. xg(i)) THEN
                jstart = j
                j = j+1
                IF (j .LE. n-1) GO TO 20
             ENDIF               

   25      CONTINUE

             IF ((x(j) .LE. xg(i+1)) .AND. (j .LE. n-1)) THEN

                a1 = AMAX1(x(j),xg(i))
                a2 = AMIN1(x(j+1),xg(i+1))

                sum = sum + y(j) * (a2-a1)/(x(j+1)-x(j))
                j = j+1
                GO TO 25

             ENDIF

           yg(i) = sum 

         ENDIF

   30 CONTINUE


* if wanted, integrate data "overhang" and fold back into last bin

      IF (FoldIn .EQ. 1) THEN

         j = j-1
         a1 = xg(ng)     ! upper limit of last interpolated bin
         a2 = x(j+1)     ! upper limit of last input bin considered

*        do folding only if grids don't match up and there is more input 
         IF ((a2 .GT. a1) .OR. (j+1 .LT. n)) THEN
           tail = y(j) * (a2-a1)/(x(j+1)-x(j))
           DO k = j+1, n-1
              tail = tail + y(k) * (x(k+1)-x(k))
           ENDDO
           yg(ng-1) = yg(ng-1) + tail
         ENDIF

      ENDIF
*_______________________________________________________________________

      RETURN
      END
      SUBROUTINE inter4(ng,xg,yg, n,x,y, FoldIn)

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Map input data given on a set of bins onto a different set of target     =*
*=  bins.                                                                    =*
*=  The input data are given on a set of bins (representing the integral     =*
*=  of the input quantity over the range of each bin) and are being matched  =*
*=  onto another set of bins (target grid).  A typical example would be an   =*
*=  input data set spcifying the extra-terrestrial flux on wavelength inter- =*
*=  vals, that has to be matched onto the working wavelength grid.           =*
*=  The resulting area in a given bin of the target grid is calculated by    =*
*=  simply adding all fractional areas of the input data that cover that     =*
*=  particular target bin.                                                   =*
*=  Some caution should be used near the endpoints of the grids.  If the     =*
*=  input data do not span the full range of the target grid, the area in    =*
*=  the "missing" bins will be assumed to be zero.  If the input data extend =*
*=  beyond the upper limit of the target grid, the user has the option to    =*
*=  integrate the "overhang" data and fold the remaining area back into the  =*
*=  last target bin.  Using this option is recommended when re-gridding      =*
*=  vertical profiles that directly affect the total optical depth of the    =*
*=  model atmosphere.                                                        =*
*-----------------------------------------------------------------------------*
*=  PARAMETERS:                                                              =*
*=  NG     - INTEGER, number of bins + 1 in the target grid               (I)=*
*=  XG     - REAL, target grid (e.g. working wavelength grid);  bin i     (I)=*
*=           is defined as [XG(i),XG(i+1)] (i = 1..NG-1)                     =*
*=  YG     - REAL, y-data re-gridded onto XG;  YG(i) specifies the        (O)=*
*=           y-value for bin i (i = 1..NG-1)                                 =*
*=  N      - INTEGER, number of bins + 1 in the input grid                (I)=*
*=  X      - REAL, input grid (e.g. data wavelength grid);  bin i is      (I)=*
*=           defined as [X(i),X(i+1)] (i = 1..N-1)                           =*
*=  Y      - REAL, input y-data on grid X;  Y(i) specifies the            (I)=*
*=           y-value for bin i (i = 1..N-1)                                  =*
*=  FoldIn - Switch for folding option of "overhang" data                 (I)=*
*=           FoldIn = 0 -> No folding of "overhang" data                     =*
*=           FoldIn = 1 -> Integerate "overhang" data and fold back into     =*
*=                         last target bin                                   =*
*-----------------------------------------------------------------------------*
*=  EDIT HISTORY:                                                            =*
*=  06/96  Added FoldIn switch                                               =*
*-----------------------------------------------------------------------------*
*= This program is free software;  you can redistribute it and/or modify     =*
*= it under the terms of the GNU General Public License as published by the  =*
*= Free Software Foundation;  either version 2 of the license, or (at your   =*
*= option) any later version.                                                =*
*= The TUV package is distributed in the hope that it will be useful, but    =*
*= WITHOUT ANY WARRANTY;  without even the implied warranty of MERCHANTIBI-  =*
*= LITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public     =*
*= License for more details.                                                 =*
*= To obtain a copy of the GNU General Public License, write to:             =*
*= Free Software Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.   =*
*-----------------------------------------------------------------------------*
*= To contact the authors, please mail to:                                   =*
*= Sasha Madronich, NCAR/ACD, P.O.Box 3000, Boulder, CO, 80307-3000, USA  or =*
*= send email to:  sasha@ucar.edu                                            =*
*-----------------------------------------------------------------------------*
*= Copyright (C) 1994,95,96  University Corporation for Atmospheric Research =*
*-----------------------------------------------------------------------------*

      IMPLICIT NONE
      
* input:
      INTEGER n, ng
      REAL xg(ng)
      REAL x(n), y(n)

      INTEGER FoldIn

* output:
      REAL yg(ng)

* local:
      REAL a1, a2, sum
      REAL tail
      INTEGER jstart, i, j, k
*_______________________________________________________________________

* check whether flag given is legal
      IF ((FoldIn .NE. 0) .AND. (FoldIn .NE. 1)) THEN
         WRITE(0,*) '>>> ERROR (inter3) <<<  Value for FOLDIN invalid. '
         WRITE(0,*) '                        Must be 0 or 1'
         STOP
      ENDIF

* do interpolation

      jstart = 1

      DO 30, i = 1, ng - 1

         yg(i) = 0.
         sum = 0.
         j = jstart

         IF (j .LE. n-1) THEN

   20      CONTINUE

             IF (x(j+1) .LT. xg(i)) THEN
                jstart = j
                j = j+1
                IF (j .LE. n-1) GO TO 20
             ENDIF               

   25      CONTINUE

           IF ((x(j) .LE. xg(i+1)) .AND. (j .LE. n-1)) THEN

              a1 = AMAX1(x(j),xg(i))
              a2 = AMIN1(x(j+1),xg(i+1))

              sum = sum + y(j) * (a2-a1)

              j = j+1
              GO TO 25

           ENDIF

           yg(i) = sum /(xg(i+1)-xg(i))

        ENDIF

 30   CONTINUE


* if wanted, integrate data "overhang" and fold back into last bin

      IF (FoldIn .EQ. 1) THEN

         j = j-1
         a1 = xg(ng)     ! upper limit of last interpolated bin
         a2 = x(j+1)     ! upper limit of last input bin considered

*        do folding only if grids don't match up and there is more input 
         IF ((a2 .GT. a1) .OR. (j+1 .LT. n)) THEN
           tail = y(j) * (a2-a1)/(x(j+1)-x(j))
           DO k = j+1, n-1
              tail = tail + y(k) * (x(k+1)-x(k))
           ENDDO
           yg(ng-1) = yg(ng-1) + tail
         ENDIF

      ENDIF
*_______________________________________________________________________

      RETURN
      END
      SUBROUTINE lymana(nz,o2col,secchi,dto2la,xso2la)

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Calculate the effective absorption cross section of O2 in the Lyman-Alpha=*
*=  bands and an effective O2 optical depth at all altitudes.  Parameterized =*
*=  after:  Chabrillat, S., and G. Kockarts, Simple parameterization of the  =*
*=  absorption of the solar Lyman-Alpha line, Geophysical Research Letters,  =*
*=  Vol.24, No.21, pp 2659-2662, 1997.                                       =*
*-----------------------------------------------------------------------------*
*=  PARAMETERS:                                                              =*
*=  NZ      - INTEGER, number of specified altitude levels in the working (I)=*
*=            grid                                                           =*
*=  O2COL   - REAL, slant overhead O2 column (molec/cc) at each specified (I)=*
*=            altitude                                                       =*
*=  DTO2LA  - REAL, optical depth due to O2 absorption at each specified  (O)=*
*=            vertical layer                                                 =*
*=  XSO2LA  - REAL, molecular absorption cross section in LA bands        (O)=*
*-----------------------------------------------------------------------------*
*=  EDIT HISTORY:                                                            =*
*=  01/98  Original                                                          =*
*-----------------------------------------------------------------------------*
*= This program is free software;  you can redistribute it and/or modify     =*
*= it under the terms of the GNU General Public License as published by the  =*
*= Free Software Foundation;  either version 2 of the license, or (at your   =*
*= option) any later version.                                                =*
*= The TUV package is distributed in the hope that it will be useful, but    =*
*= WITHOUT ANY WARRANTY;  without even the implied warranty of MERCHANTIBI-  =*
*= LITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public     =*
*= License for more details.                                                 =*
*= To obtain a copy of the GNU General Public License, write to:             =*
*= Free Software Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.   =*
*-----------------------------------------------------------------------------*
*= To contact the authors, please mail to:                                   =*
*= Sasha Madronich, NCAR/ACD, P.O.Box 3000, Boulder, CO, 80307-3000, USA  or =*
*= send email to:  sasha@ucar.edu                                            =*
*-----------------------------------------------------------------------------*
*= Copyright (C) 1994 - 1998 University Corporation for Atmospheric Research =*
*-----------------------------------------------------------------------------*

      IMPLICIT NONE
      INCLUDE 'tuv.params'

      INTEGER nz
      REAL o2col(kz)
      REAL secchi(kz)

      REAL dto2la(kz,*), xso2la(kz,*)

* local variables
      DOUBLE PRECISION RM(kz), RO2(kz)
      DOUBLE PRECISION b(3), c(3), d(3), e(3)
      DATA b/ 6.8431D-01, 2.29841D-01,  8.65412D-02/,
     >     c/8.22114D-21, 1.77556D-20,  8.22112D-21/,
     >     d/ 6.0073D-21, 4.28569D-21,  1.28059D-20/,
     >     e/8.21666D-21, 1.63296D-20,  4.85121D-17/

      INTEGER iz, i

*------------------------------------------------------------------------------*

C calculate reduction factors at every altitude

      DO iz = 1, nz
        RM(iz) = 0.D+00
        RO2(iz) = 0.D+00
        DO i = 1, 3
          RM(iz) = RM(iz) + b(i) * DEXP(-c(i) * DBLE(o2col(iz)))
          RO2(iz) = RO2(iz) + d(i) * DEXP(-e(i) * DBLE(o2col(iz)))
        ENDDO
      ENDDO

C calculate effective O2 optical depths and effective O2 cross sections

      DO iz = 1, nz-1

         IF (rm(iz) .GT. 1.0D-100) THEN
            IF (ro2(iz) .GT. 1.D-100) THEN
               xso2la(iz,1) = ro2(iz)/rm(iz)
            ELSE
               xso2la(iz,1) = 0.               
            ENDIF

            IF (rm(iz+1) .GT. 0.) THEN

               dto2la(iz,1) = LOG(rm(iz+1)) / secchi(iz+1) 
     $                      - LOG(rm(iz))   / secchi(iz)

            ELSE
               dto2la(iz,1) = 1000.
            ENDIF
         ELSE
            dto2la(iz,1) = 1000.
            xso2la(iz,1) = 0.
         ENDIF

      ENDDO

C do top layer separate

      dto2la(nz,1) = 0.
      IF(rm(nz) .GT. 1.D-100) THEN
         xso2la(nz,1) = RO2(nz)/RM(nz)
      ELSE
         xso2la(nz,1) = 0.
      ENDIF

*------------------------------------------------------------------------------*

      END
      SUBROUTINE pbiol1(nw,wl,wc,j,s,label)

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Create or read various weighting functions, e.g. biological action       =*
*=  spectra, instrument responses etc.                                       =*
*-----------------------------------------------------------------------------*
*=  PARAMETERS:                                                              =*
*=  NW     - INTEGER, number of specified intervals + 1 in working        (I)=*
*=           wavelength grid                                                 =*
*=  WL     - REAL, vector of lower limits of wavelength intervals in      (I)=*
*=           working wavelength grid                                         =*
*=  WC     - REAL, vector of central wavelength of wavelength intervals    I)=*
*=           in working wavelength grid                                      =*
*=  J      - INTEGER, counter for number of weighting functions defined  (IO)=*
*=  S      - REAL, value of each defined weighting function at each       (O)=*
*=           defined wavelength                                              =*
*=  LABEL  - CHARACTER*40, string identifier for each weighting function  (O)=*
*=           defined                                                         =*
*-----------------------------------------------------------------------------*
*=  EDIT HISTORY:                                                            =*
*=  02/97  Changed offset for grid-end interpolation to relative number      =*
*=         (x * (1 +- deltax))                                               =*
*=  05/96  Renamed from LOADW1 to WSPEC1                                     =*
*=  03/95  Added vis+ and UV index                                           =*
*-----------------------------------------------------------------------------*
*= This program is free software;  you can redistribute it and/or modify     =*
*= it under the terms of the GNU General Public License as published by the  =*
*= Free Software Foundation;  either version 2 of the license, or (at your   =*
*= option) any later version.                                                =*
*= The TUV package is distributed in the hope that it will be useful, but    =*
*= WITHOUT ANY WARRANTY;  without even the implied warranty of MERCHANTIBI-  =*
*= LITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public     =*
*= License for more details.                                                 =*
*= To obtain a copy of the GNU General Public License, write to:             =*
*= Free Software Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.   =*
*-----------------------------------------------------------------------------*
*= To contact the authors, please mail to:                                   =*
*= Sasha Madronich, NCAR/ACD, P.O.Box 3000, Boulder, CO, 80307-3000, USA  or =*
*= send email to:  sasha@ucar.edu                                            =*
*-----------------------------------------------------------------------------*
*= Copyright (C) 1994,95,96  University Corporation for Atmospheric Research =*
*-----------------------------------------------------------------------------*

      IMPLICIT NONE
      INCLUDE 'tuv.params'

      INTEGER kdata
      PARAMETER(kdata=1000)

* input:
      REAL wl(kw), wc(kw)
      INTEGER nw

* input/output:
      INTEGER j

* output: (weighting functions and labels)
      REAL s(ks,kw)
      CHARACTER*40 label(ks)

* internal:
      REAL x1(kdata)
      REAL y1(kdata)
      REAL yg(kw)

      REAL fery, futr
      EXTERNAL fery, futr
      INTEGER i, iw, n

      INTEGER ierr

      INTEGER idum
      REAL dum1, dum2
      REAL em, a, b, c
      REAL sum

*_______________________________________________________________________

********* UV-B (280-315 nm)
 
      j = j + 1
      label(j) = 'UV-B, 280-315 nm'
      DO iw = 1, nw-1
         IF (wc(iw) .GT. 280. .AND. wc(iw) .LT. 315.) THEN
            s(j,iw) = 1.
         ELSE
            s(j,iw) = 0.
         ENDIF
      ENDDO

********* UV-B* (280-320 nm)
 
      j = j + 1
      label(j) = 'UV-B*, 280-320 nm'
      DO iw = 1, nw-1
         IF (wc(iw) .GT. 280. .AND. wc(iw) .LT. 320.) THEN
            s(j,iw) = 1.
         ELSE
            s(j,iw) = 0.
         ENDIF
      ENDDO

********* UV-A (315-400 nm)
 
      j = j + 1
      label(j) = 'UV-A, 315-400 nm'
      DO iw = 1, nw-1
         IF (wc(iw) .GT. 315. .AND. wc(iw) .LT. 400.) THEN
            s(j,iw) = 1.
         ELSE
            s(j,iw) = 0.
         ENDIF
      ENDDO

********* visible+ (> 400 nm)
 
      j = j + 1
      label(j) = 'vis+, > 400 nm'
      DO iw = 1, nw-1
         IF (wc(iw) .GT. 400.) THEN
            s(j,iw) = 1.
         ELSE
            s(j,iw) = 0.
         ENDIF
      ENDDO

********** unity raf constant slope:  

      j = j + 1
      label(j) = 'decay, 14 nm/10'
      DO iw = 1, nw-1
         s(j,iw) = 10.**(-(wc(iw) -300.)/14.)
      ENDDO

************ DNA damage action spectrum
* from: Setlow, R. B., The wavelengths in sunlight effective in 
*       producing skin cancer: a theoretical analysis, Proceedings 
*       of the National Academy of Science, 71, 3363 -3366, 1974.
* normalize to unity at 300 nm
* Data read from original hand-drawn plot by Setlow
* received from R. Setlow in May 1995
* data is per quantum (confirmed with R. Setlow in May 1995).  
* Therefore must put on energy basis if irradiance is is energy
* (rather than quanta) units.

      j = j + 1
      label(j) = 'Setlow dna.new'
      OPEN(UNIT=kin,FILE='DATAS1/dna.setlow.new',STATUS='old')
      do i = 1, 11
         read(kin,*)
      enddo
      n = 55
      DO i = 1, n
         READ(kin,*) x1(i), y1(i)
         y1(i) = y1(i) / 2.4E-02  *  x1(i)/300.
      ENDDO
      CLOSE (kin)

      CALL addpnt(x1,y1,kdata,n,x1(1)*(1.-deltax),y1(1))
      CALL addpnt(x1,y1,kdata,n,          0.,y1(1))
      CALL addpnt(x1,y1,kdata,n,x1(n)*(1.+deltax),   0.)
      CALL addpnt(x1,y1,kdata,n,      1.e+38,   0.)
      CALL inter2(nw,wl,yg,n,x1,y1,ierr)
      IF (ierr .NE. 0) THEN
         WRITE(*,*) ierr, label(j)
         STOP
      ENDIF

      DO iw = 1, nw-1
         s(j,iw) = yg(iw)
      ENDDO

********* skin cancer in mice,  Utrecht/Phildelphia study
*from de Gruijl, F. R., H. J. C. M. Sterenborg, P. D. Forbes, 
*     R. E. Davies, C. Cole, G. Kelfkens, H. van Weelden, H. Slaper,
*     and J. C. van der Leun, Wavelength dependence of skin cancer 
*     induction by ultraviolet irradiation of albino hairless mice, 
*     Cancer Res., 53, 53-60, 1993.
* normalize at 300 nm.

      j = j + 1
      label(j) = 'SCUP-m'
      DO iw = 1, nw-1
         s(j,iw) =  futr(wc(iw)) / futr(300.)
      ENDDO
         
*********** Utrecht mice spectrum corrected for humans skin.

      j = j + 1
      label(j) = 'SCUP-h'
      OPEN(UNIT=kin,FILE='DATAS1/SCUP-h',STATUS='old')
      n = 28
      DO i = 1, n
         READ(kin,*) x1(i), y1(i)
      ENDDO

      CALL addpnt(x1,y1,kdata,n,x1(1)*(1.-deltax),y1(1))
      CALL addpnt(x1,y1,kdata,n,          0.,y1(1))
      CALL addpnt(x1,y1,kdata,n,x1(n)*(1.+deltax),   0.)
      CALL addpnt(x1,y1,kdata,n,      1.e+38,   0.)
      CALL inter2(nw,wl,yg,n,x1,y1,ierr)
      IF (ierr .NE. 0) THEN
         WRITE(*,*) ierr, label(j)
         STOP
      ENDIF
            
      DO iw = 1, nw-1
         s(j,iw) = yg(iw)
      ENDDO
      CLOSE (kin)
      
***************** compute ery action spectrum
*from
* McKinlay, A. F., and B. L. Diffey, A reference action spectrum for 
* ultraviolet induced erythema in human skin, in Human Exposure to 
* Ultraviolet Radiation: Risks and Regulations, W. R. Passchler 
* and B. F. M. Bosnajokovic, (eds.), Elsevier, Amsterdam, 1987.

      j = j + 1
      label(j) = 'CIE hum erythema'
      DO iw = 1, nw-1
         s(j,iw) = fery(wc(iw))
      ENDDO

***************** UV index (Canadian - WMO/WHO)
* based on erythema

      j = j + 1
      label(j) = 'UV index'
      DO iw = 1, nw-1
         s(j,iw) = 40. * fery(wc(iw))
      ENDDO

************* erythema - Anders et al.
* for skin types II and III, from Anders et al., Photochem. and
* Photobiol., 61, 200-203, 1995. Units are J m-2.

      j = j + 1
      label(j) = 'ery.anders'
      OPEN(UNIT=kin,FILE='DATAS1/ery.anders',STATUS='old')
      do i = 1, 3
         read(kin,*)
      enddo
      n = 28
      DO i = 1, n
         READ(kin,*) x1(i), y1(i)
         y1(i) = 1./y1(i)
      ENDDO

      CALL addpnt(x1,y1,kdata,n,x1(1)*(1.-deltax),y1(1))
      CALL addpnt(x1,y1,kdata,n,          0.,y1(1))
      CALL addpnt(x1,y1,kdata,n,x1(n)*(1.+deltax),   0.)
      CALL addpnt(x1,y1,kdata,n,      1.e+38,   0.)
      CALL inter2(nw,wl,yg,n,x1,y1,ierr)
      IF (ierr .NE. 0) THEN
         WRITE(*,*) ierr, label(j)
         STOP
      ENDIF
            
      DO iw = 1, nw-1
         s(j,iw) = yg(iw)
      ENDDO
      CLOSE (kin)

********* read 1991-92 acgih threshold limit values
* from
* ACGIH, 1991-1992 Threshold Limit Values, American Conference 
*  of Governmental and Industrial Hygienists, 1992.

      j = j + 1
      label(j) = 'ACGIH 1992 TLVs'
      OPEN(UNIT=kin,FILE='DATAS1/acgih.1992',STATUS='old')
      n = 56
      DO i = 1, n
         READ(kin,*) x1(i), y1(i)
         y1(i) = y1(i)
      ENDDO

      CALL addpnt(x1,y1,kdata,n,x1(1)*(1.-deltax),y1(1))
      CALL addpnt(x1,y1,kdata,n,          0.,y1(1))
      CALL addpnt(x1,y1,kdata,n,x1(n)*(1.+deltax),   0.)
      CALL addpnt(x1,y1,kdata,n,      1.e+38,   0.)
      CALL inter2(nw,wl,yg,n,x1,y1,ierr)
      IF (ierr .NE. 0) THEN
         WRITE(*,*) ierr, label(j)
         STOP
      ENDIF
            
      DO iw = 1, nw-1
         s(j,iw) = yg(iw)
      ENDDO
      CLOSE (kin)

********** RB Meter, model 501
*  private communication, M. Morys (Solar Light Co.), 1994.
* From: morys@omni.voicenet.com (Marian Morys)
* Received: from acd.ucar.edu by sasha.acd.ucar.edu (AIX 3.2/UCB 5.64/4.03)
*          id AA17274; Wed, 21 Sep 1994 11:35:44 -0600

      j = j + 1
      label(j) = 'RB Meter, model 501'
      OPEN(UNIT=kin,FILE='DATAS1/rbm.501',STATUS='old')
      n = 57
      DO i = 1, n
         READ(kin,*) x1(i), y1(i)
      ENDDO
      
      CALL addpnt(x1,y1,kdata,n,x1(1)*(1.-deltax),y1(1))
      CALL addpnt(x1,y1,kdata,n,          0.,y1(1))
      CALL addpnt(x1,y1,kdata,n,x1(n)*(1.+deltax),   0.)
      CALL addpnt(x1,y1,kdata,n,      1.e+38,   0.)
      CALL inter2(nw,wl,yg,n,x1,y1,ierr)
      IF (ierr .NE. 0) THEN
         WRITE(*,*) ierr, label(j)
         STOP
      ENDIF
            
      DO iw = 1, nw-1
         s(j,iw) = yg(iw)
      ENDDO
      CLOSE (kin)

********* phytoplankton, Boucher et al. (1994) 
* from Boucher, N., Prezelin, B.B., Evens, T., Jovine, R., Kroon, B., Moline, M.A.,
* and Schofield, O., Icecolors '93: Biological weighting function for the ultraviolet
*  inhibition  of carbon fixation in a natural antarctic phytoplankton community, 
* Antarctic Journal, Review 1994, pp. 272-275, 1994.
* In original paper, value of b and m (em below are given as positive.  Correct values
* are negative. Also, limit to positive values.

      j = j + 1
      label(j) = 'phytoplankton, Boucher et al. (1994)'
      a = 112.5
      b = -6.223E-01
      c = 7.670E-04
      em = -3.17E-06
      DO iw = 1, nw-1
         IF (wc(iw) .GT. 290. .AND. wc(iw) .LT. 400.) THEN
            s(j,iw) = em + EXP(a+b*wc(iw)+c*wc(iw)*wc(iw))
         ELSE
            s(j,iw) = 0.
         ENDIF
         s(j,iw) = max(s(j,iw),0.)
      ENDDO

********* phytoplankton, Cullen et al.
* Cullen, J.J., Neale, P.J., and Lesser, M.P., Biological weighting function for the  
*  inhibition of phytoplankton photosynthesis by ultraviolet radiation, Science, 25,
*  646-649, 1992.
* phaeo

      j = j + 1
      label(j) = 'Cullen, phaeo'
      OPEN(UNIT=kin,FILE='DATAS1/phaeo.bio',STATUS='old')
      n = 106
      DO i = 1, n
         READ(kin,*) idum, dum1, dum2, y1(i)
         x1(i) = (dum1+dum2)/2.
      ENDDO

      CALL addpnt(x1,y1,kdata,n,x1(1)*(1.-deltax),y1(1))
      CALL addpnt(x1,y1,kdata,n,          0.,y1(1))
      CALL addpnt(x1,y1,kdata,n,x1(n)*(1.+deltax),   0.)
      CALL addpnt(x1,y1,kdata,n,      1.e+38,   0.)
      CALL inter2(nw,wl,yg,n,x1,y1,ierr)
      IF (ierr .NE. 0) THEN
         WRITE(*,*) ierr, label(j)
         STOP
      ENDIF

      DO iw = 1, nw-1
         s(j,iw) = yg(iw)
      ENDDO
      CLOSE(kin)

* proro

      j = j + 1
      label(j) = 'Cullen, proro'
      OPEN(UNIT=kin,FILE='DATAS1/proro.bio',STATUS='old')
      n = 100
      DO i = 1, n
         READ(kin,*) idum, dum1, dum2, y1(i)
         x1(i) = (dum1+dum2)/2.
      ENDDO

      CALL addpnt(x1,y1,kdata,n,x1(1)*(1.-deltax),y1(1))
      CALL addpnt(x1,y1,kdata,n,          0.,y1(1))
      CALL addpnt(x1,y1,kdata,n,x1(n)*(1.+deltax),   0.)
      CALL addpnt(x1,y1,kdata,n,      1.e+38,   0.)
      CALL inter2(nw,wl,yg,n,x1,y1,ierr)
      IF (ierr .NE. 0) THEN
         WRITE(*,*) ierr, label(j)
         STOP
      ENDIF

      DO iw = 1, nw-1
         s(j,iw) = yg(iw)
      ENDDO
      CLOSE (kin)

*  Gaussians

      j = j + 1
      label(j) = 'Gaussian, 305 nm, 10 nm FWHM'
      sum = 0.
      DO iw = 1, nw-1
         s(j,iw) = exp(- ( log(2.) * ((wc(iw)-305.)/(5.))**2) )
         sum = sum + s(j,iw)
      ENDDO
      DO iw = 1, nw-1
         s(j,iw) = s(j,iw)/sum
      ENDDO

      j = j + 1
      label(j) = 'Gaussian, 320 nm, 10 nm FWHM'
      sum = 0.
      DO iw = 1, nw-1
         s(j,iw) = exp(- ( log(2.) * ((wc(iw)-320.)/(5.))**2) )
         sum = sum + s(j,iw)
      ENDDO
      DO iw = 1, nw-1
         s(j,iw) = s(j,iw)/sum
      ENDDO

      j = j + 1
      label(j) = 'Gaussian, 340 nm, 10 nm FWHM'
      sum = 0.
      DO iw = 1, nw-1
         s(j,iw) = exp(- ( log(2.) * ((wc(iw)-340.)/(5.))**2) )
         sum = sum + s(j,iw)
      ENDDO
      DO iw = 1, nw-1
         s(j,iw) = s(j,iw)/sum
      ENDDO

      j = j + 1
      label(j) = 'Gaussian, 380 nm, 10 nm FWHM'
      sum = 0.
      DO iw = 1, nw-1
         s(j,iw) = exp(- ( log(2.) * ((wc(iw)-380.)/(5.))**2) )
         sum = sum + s(j,iw)
      ENDDO
      DO iw = 1, nw-1
         s(j,iw) = s(j,iw)/sum
      ENDDO

****************************************************************
****************************************************************

*_______________________________________________________________________

      IF (j .GT. ks) STOP '1001'
*_______________________________________________________________________

      RETURN
      END
      SUBROUTINE pchem(nw,wl,nz,tlev,airlev,
     $     j,sq,jlabel)

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Load various "weighting functions" (products of cross section and        =*
*=  quantum yield at each altitude and for wavelength).  The altitude        =*
*=  dependence is necessary to ensure the consideration of pressure and      =*
*=  temperature dependence of the cross sections or quantum yields.          =*
*=  The actual reading, evaluation and interpolation is done is separate     =*
*=  subroutines for ease of management and manipulation.  Please refer to    =*
*=  the inline documentation of the specific subroutines for detail          =*
*=  information.                                                             =*
*-----------------------------------------------------------------------------*
*=  PARAMETERS:                                                              =*
*=  NW     - INTEGER, number of specified intervals + 1 in working        (I)=*
*=           wavelength grid                                                 =*
*=  WL     - REAL, vector of lower limits of wavelength intervals in      (I)=*
*=           working wavelength grid                                         =*
*=  NZ     - INTEGER, number of altitude levels in working altitude grid  (I)=*
*=  TLEV   - REAL, temperature (K) at each specified altitude level       (I)=*
*=  AIRLEV - REAL, air density (molec/cc) at each altitude level          (I)=*
*=  J      - INTEGER, counter for number of weighting functions defined  (IO)=*
*=  SQ     - REAL, cross section x quantum yield (cm^2) for each          (O)=*
*=           photolysis reaction defined, at each defined wavelength and     =*
*=           at each defined altitude level                                  =*
*=  JLABEL - CHARACTER*40, string identifier for each photolysis reaction (O)=*
*=           defined                                                         =*
*-----------------------------------------------------------------------------*
*=  EDIT HISTORY:                                                            =*
*=  05/98  Original; adapted from the "old" JSPEC1 routine                   =*
*-----------------------------------------------------------------------------*
*= This program is free software;  you can redistribute it and/or modify     =*
*= it under the terms of the GNU General Public License as published by the  =*
*= Free Software Foundation;  either version 2 of the license, or (at your   =*
*= option) any later version.                                                =*
*= The TUV package is distributed in the hope that it will be useful, but    =*
*= WITHOUT ANY WARRANTY;  without even the implied warranty of MERCHANTIBI-  =*
*= LITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public     =*
*= License for more details.                                                 =*
*= To obtain a copy of the GNU General Public License, write to:             =*
*= Free Software Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.   =*
*-----------------------------------------------------------------------------*
*= To contact the authors, please mail to:                                   =*
*= Sasha Madronich, NCAR/ACD, P.O.Box 3000, Boulder, CO, 80307-3000, USA  or =*
*= send email to:  sasha@ucar.edu                                            =*
*-----------------------------------------------------------------------------*
*= Copyright (C) 1994,95,96  University Corporation for Atmospheric Research =*
*-----------------------------------------------------------------------------*

      IMPLICIT NONE
      INCLUDE 'tuv.params'

* input

      INTEGER nw
      REAL wl(kw)
      
      INTEGER nz

      REAL tlev(kz)
      REAL airlev(kz)

* weighting functions

      CHARACTER*40 jlabel(kj),jlabeltmp
      REAL sq(kj,kz,kw),sqtmp(kz,kw)

* input/output:
      INTEGER j

* local:
      REAL wc(kw), wu(kw)
      INTEGER iw
*_______________________________________________________________________

* complete wavelength grid

      DO 5, iw = 1, nw - 1
         wc(iw) = (wl(iw) + wl(iw+1))/2.
         wu(iw) =  wl(iw+1)
 5    CONTINUE

*____________________________________________________________________________

C O2 + hv -> O + O
* reserve first position.  Cross section parameterization in Schumman-Runge and 
* Lyman-alpha regions are zenith-angle dependent, will be written in 
* subroutine sto2xs(nz,nw,xso2,nj,sj).
 
      j = 0

C NO2 + hv -> NO + O(3P)
      CALL rn_NO2(nw,wl,wc,nz,tlev,airlev,j,sq,jlabel)

C NO3 + hv ->  (both channels)
      CALL rn_N03(nw,wl,wc,nz,tlev,airlev,j,sq,jlabel)

C O3 + hv ->  (both channels)
      CALL rn_O3(nw,wl,wc,nz,tlev,airlev,j,sq,jlabel)

C HNO2 + hv -> OH + NO
C HNO2 + hv -> HO2 + NO2 
      CALL rn_HNO2(nw,wl,wc,nz,tlev,airlev,j,sq,jlabel)

C HNO3 + hv -> OH + NO2
      CALL rn_HNO3(nw,wl,wc,nz,tlev,airlev,j,sq,jlabel)

C HNO4 + hv -> HO2 + NO2
      CALL rn_HNO4(nw,wl,wc,nz,tlev,airlev,j,sq,jlabel)

C H2O2 + hv -> 2 OH
      CALL rn_H2O2(nw,wl,wc,nz,tlev,airlev,j,sq,jlabel)

C CH2O + hv -> (both channels)
      CALL rn_HCHO(nw,wl,wc,nz,tlev,airlev,j,sq,jlabel)
     
C CH3CHO + hv -> CO + HO2 + C-O2
      CALL rn_CCHO(nw,wl,wc,nz,tlev,airlev,j,sq,jlabel)

C C2CHO + hv -> CCHO + RO2_R + CO + HO2
      CALL rn_RCHO(nw,wl,wc,nz,tlev,airlev,j,sq,jlabel)

C CH3COCH3 + hv -> Products
      CALL rn_acet(nw,wl,wc,nz,tlev,airlev,j,sq,jlabel)

C MEK + hv -> Products and PROD2
      CALL rn_mek(nw,wl,wc,nz,tlev,airlev,j,sq,jlabel)
      jlabeltmp=jlabel(j)      ! move PROD2 to the 27th
      do iw=1,nw-1
        sqtmp(1:nz,iw)=sq(j,1:nz,iw)
      enddo 	
      j=j-1

C CH3OOH + hv -> CH3O + OH and ROOH
      CALL rn_COOH(nw,wl,wc,nz,tlev,airlev,j,sq,jlabel)

C GLY + hv -> Products (both channels)
      CALL rn_GLY(nw,wl,wc,nz,tlev,airlev,j,sq,jlabel)

C MGLY + hv -> Products 
      CALL rn_MGLY(nw,wl,wc,nz,tlev,airlev,j,sq,jlabel)

C BACL + hv -> Products (two channels)
      CALL rn_BACL(nw,wl,wc,nz,tlev,airlev,j,sq,jlabel)

C ACROLEIN (four channels)
      CALL rn_ACROLEIN(nw,wl,wc,nz,tlev,airlev,j,sq,jlabel)
      jlabel(30)=jlabel(j)              ! move DCB3 to the last 
      do iw=1,nw-1       
       sq(30,1:nz,iw)=sq(j,1:nz,iw)
       sq(27,1:nz,iw)=sqtmp(1:nz,iw)
      enddo 
      jlabel(27)=jlabeltmp              ! move PROD2 to the 27th
      j=27
      
C RNO3 -> products
      CALL rn_RNO3(nw,wl,wc,nz,tlev,airlev,j,sq,jlabel)

C DCB2 -> products
      CALL rn_DCB2(nw,wl,wc,nz,tlev,airlev,j,sq,jlabel)
      
      j=30
       
****************************************************************

      IF (j .GT. kj) STOP '1002'
      RETURN
      END
      SUBROUTINE ps2str(nlevel,zen,rsfc,tauu,omu,gu,
     $     dsdh, nid, delta,
     $     fdr, fup, fdn, edr, eup, edn)
*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Solve two-stream equations for multiple layers.  The subroutine is based =*
*=  on equations from:  Toon et al., J.Geophys.Res., v94 (D13), Nov 20, 1989.=*
*=  It contains 9 two-stream methods to choose from.  A pseudo-spherical     =*
*=  correction has also been added.                                          =*
*-----------------------------------------------------------------------------*
*=  PARAMETERS:                                                              =*
*=  NLEVEL  - INTEGER, number of specified altitude levels in the working (I)=*
*=            grid                                                           =*
*=  ZEN     - REAL, solar zenith angle (degrees)                          (I)=*
*=  RSFC    - REAL, surface albedo at current wavelength                  (I)=*
*=  TAUU    - REAL, unscaled optical depth of each layer                  (I)=*
*=  OMU     - REAL, unscaled single scattering albedo of each layer       (I)=*
*=  GU      - REAL, unscaled asymmetry parameter of each layer            (I)=*
*=  DSDH    - REAL, slant path of direct beam through each layer crossed  (I)=*
*=            when travelling from the top of the atmosphere to layer i;     =*
*=            DSDH(i,j), i = 0..NZ-1, j = 1..NZ-1                            =*
*=  NID     - INTEGER, number of layers crossed by the direct beam when   (I)=*
*=            travelling from the top of the atmosphere to layer i;          =*
*=            NID(i), i = 0..NZ-1                                            =*
*=  DELTA   - LOGICAL, switch to use delta-scaling                        (I)=*
*=            .TRUE. -> apply delta-scaling                                  =*
*=            .FALSE.-> do not apply delta-scaling                           =*
*=  FDR     - REAL, contribution of the direct component to the total     (O)=*
*=            actinic flux at each altitude level                            =*
*=  FUP     - REAL, contribution of the diffuse upwelling component to    (O)=*
*=            the total actinic flux at each altitude level                  =*
*=  FDN     - REAL, contribution of the diffuse downwelling component to  (O)=*
*=            the total actinic flux at each altitude level                  =*
*=  EDR     - REAL, contribution of the direct component to the total     (O)=*
*=            spectral irradiance at each altitude level                     =*
*=  EUP     - REAL, contribution of the diffuse upwelling component to    (O)=*
*=            the total spectral irradiance at each altitude level           =*
*=  EDN     - REAL, contribution of the diffuse downwelling component to  (O)=*
*=            the total spectral irradiance at each altitude level           =*
*-----------------------------------------------------------------------------*
*=  EDIT HISTORY:                                                            =*
*=  02/96  Added pseudo-spherical correction                                 =*
*=  05/94  Added various two-stream methods                                  =*
*-----------------------------------------------------------------------------*
*= This program is free software;  you can redistribute it and/or modify     =*
*= it under the terms of the GNU General Public License as published by the  =*
*= Free Software Foundation;  either version 2 of the license, or (at your   =*
*= option) any later version.                                                =*
*= The TUV package is distributed in the hope that it will be useful, but    =*
*= WITHOUT ANY WARRANTY;  without even the implied warranty of MERCHANTIBI-  =*
*= LITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public     =*
*= License for more details.                                                 =*
*= To obtain a copy of the GNU General Public License, write to:             =*
*= Free Software Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.   =*
*-----------------------------------------------------------------------------*
*= To contact the authors, please mail to:                                   =*
*= Sasha Madronich, NCAR/ACD, P.O.Box 3000, Boulder, CO, 80307-3000, USA  or =*
*= send email to:  sasha@ucar.edu                                            =*
*-----------------------------------------------------------------------------*
*= Copyright (C) 1994,95,96  University Corporation for Atmospheric Research =*
*-----------------------------------------------------------------------------*

      IMPLICIT NONE
      INCLUDE 'tuv.params'

      INTEGER nrows
      PARAMETER(nrows=2*kz)

*******
* input:
*******
      INTEGER nlevel
      REAL zen, rsfc
      REAL tauu(kz), omu(kz), gu(kz)
      REAL dsdh(0:kz,kz)
      INTEGER nid(0:kz)
      LOGICAL delta

*******
* output:
*******
      REAL fup(kz),fdn(kz),fdr(kz)
      REAL eup(kz),edn(kz),edr(kz)

*******
* local:
*******
      REAL tausla(0:kz), tauc(0:kz)
      REAL mu2(0:kz), mu, sum

* internal coefficients and matrix
      REAL lam(kz),taun(kz),bgam(kz)
      REAL e1(kz),e2(kz),e3(kz),e4(kz)
      REAL cup(kz),cdn(kz),cuptn(kz),cdntn(kz)
      REAL mu1(kz)
      INTEGER row
      REAL a(nrows),b(nrows),d(nrows),e(nrows),y(nrows)

*******
* other:
*******
      REAL pifs, fdn0
      REAL gi(kz), omi(kz), tempg
      REAL f, g, om
      REAL gam1, gam2, gam3, gam4

* For calculations of Associated Legendre Polynomials for GAMA1,2,3,4
* in delta-function, modified quadrature, hemispheric constant,
* Hybrid modified Eddington-delta function metods, p633,Table1.
* W.E.Meador and W.R.Weaver, GAS,1980,v37,p.630
* W.J.Wiscombe and G.W. Grams, GAS,1976,v33,p2440, 
* uncomment the following two lines and the appropriate statements further
* down.
C     REAL YLM0, YLM2, YLM4, YLM6, YLM8, YLM10, YLM12, YLMS, BETA0,
C    >     BETA1, BETAn, amu1, subd

      REAL expon, expon0, expon1, divisr, temp, up, dn
      REAL ssfc
      INTEGER nlayer, mrows, lev

      INTEGER i, j

* Some additional program constants:
      REAL eps, precis
      PARAMETER (eps = 1.E-3, precis = 1.E-7) 
*_______________________________________________________________________

* MU = cosine of solar zenith angle
* RSFC = surface albedo
* TAUU =  unscaled optical depth of each layer
* OMU  =  unscaled single scattering albedo
* GU   =  unscaled asymmetry factor
* KLEV = max dimension of number of layers in atmosphere
* NLAYER = number of layers in the atmosphere
* NLEVEL = nlayer + 1 = number of levels

* initial conditions:  pi*solar flux = 1;  diffuse incidence = 0

      pifs = 1.      
      fdn0 = 0.

      nlayer = nlevel - 1

       mu = COS(zen*pi/180.)

************** compute coefficients for each layer:
* GAM1 - GAM4 = 2-stream coefficients, different for different approximations
* EXPON0 = calculation of e when TAU is zero
* EXPON1 = calculation of e when TAU is TAUN
* CUP and CDN = calculation when TAU is zero
* CUPTN and CDNTN = calc. when TAU is TAUN
* DIVISR = prevents division by zero

        do j = 0, kz
           tauc(j) = 0.
           tausla(j) = 0.
           mu2(j) = 1./SQRT(largest)

        end do

       IF( .NOT. delta ) THEN
         DO i = 1, nlayer
           gi(i) = gu(i)
           omi(i) = omu(i)
           taun(i) = tauu(i)
         ENDDO
       ELSE 

* delta-scaling. Have to be done for delta-Eddington approximation, 
* delta discrete ordinate, Practical Improved Flux Method, delta function,
* and Hybrid modified Eddington-delta function methods approximations

         DO i = 1, nlayer
           f = gu(i)*gu(i)
           gi(i) = (gu(i) - f)/(1 - f)
           omi(i) = (1 - f)*omu(i)/(1 - omu(i)*f)       
           taun(i) = (1 - omu(i)*f)*tauu(i)
         ENDDO
        END IF

*
* calculate slant optical depth at the top of the atmosphere when zen>90.
* in this case, higher altitude of the top layer is recommended which can 
* be easily changed in gridz.f.
*
         IF(zen .GT. 90.0) THEN
           IF(nid(0) .LT. 0) THEN
             tausla(0) = largest
           ELSE
             sum = 0.0
             DO j = 1, nid(0)
              sum = sum + 2.*taun(j)*dsdh(0,j)
             END DO
             tausla(0) = sum 
           END IF
         END IF
  
*
        DO 11, i = 1, nlayer

         g = gi(i)
         om = omi(i)
         tauc(i) = tauc(i-1) + taun(i)

* stay away from 1 by precision.  For g, also stay away from -1

         tempg = AMIN1(abs(g),1. - precis)
         g = SIGN(tempg,g)
         om = AMIN1(om,1.-precis)


* calculate slant optical depth
*              
          IF(nid(i) .LT. 0) THEN
            tausla(i) = largest
          ELSE
            sum = 0.0
            DO j = 1, MIN(nid(i),i)
               sum = sum + taun(j)*dsdh(i,j)
            ENDDO
            DO j = MIN(nid(i),i)+1,nid(i)
               sum = sum + 2.*taun(j)*dsdh(i,j)
            ENDDO
            tausla(i) = sum 
            IF(tausla(i) .EQ. tausla(i-1)) THEN
              mu2(i) = SQRT(largest)
            ELSE
              mu2(i) = (tauc(i)-tauc(i-1))/(tausla(i)-tausla(i-1))
              mu2(i) = SIGN( AMAX1(ABS(mu2(i)),1./SQRT(largest)),
     $                     mu2(i) )
            END IF
          END IF
*
*** the following gamma equations are from pg 16,289, Table 1
*** save mu1 for each approx. for use in converting irradiance to actinic flux

* Eddington approximation(Joseph et al., 1976, JAS, 33, 2452):

        gam1 =  (7. - om*(4. + 3.*g))/4.
        gam2 = -(1. - om*(4. - 3.*g))/4.
        gam3 = (2. - 3.*g*mu)/4.
        gam4 = 1. - gam3
        mu1(i) = 0.5

* quadrature (Liou, 1973, JAS, 30, 1303-1326; 1974, JAS, 31, 1473-1475):

c          gam1 = 1.7320508*(2. - om*(1. + g))/2.
c          gam2 = 1.7320508*om*(1. - g)/2.
c          gam3 = (1. - 1.7320508*g*mu)/2.
c          gam4 = 1. - gam3
c          mu1(i) = 1./sqrt(3.)
         
* hemispheric mean (Toon et al., 1089, JGR, 94, 16287):

c          gam1 = 2. - om*(1. + g)
c          gam2 = om*(1. - g)
c          gam3 = (2. - g*mu)/4.
c          gam4 = 1. - gam3
c          mu1(i) = 0.5

* PIFM  (Zdunkovski et al.,1980, Conrib.Atmos.Phys., 53, 147-166):
c         GAM1 = 0.25*(8. - OM*(5. + 3.*G))
c         GAM2 = 0.75*OM*(1.-G)
c         GAM3 = 0.25*(2.-3.*G*MU)
c         GAM4 = 1. - GAM3
c         mu1(i) = 0.5

* delta discrete ordinates  (Schaller, 1979, Contrib.Atmos.Phys, 52, 17-26):
c         GAM1 = 0.5*1.7320508*(2. - OM*(1. + G))
c         GAM2 = 0.5*1.7320508*OM*(1.-G)
c         GAM3 = 0.5*(1.-1.7320508*G*MU)
c         GAM4 = 1. - GAM3
c         mu1(i) = 1./sqrt(3.)

* Calculations of Associated Legendre Polynomials for GAMA1,2,3,4
* in delta-function, modified quadrature, hemispheric constant,
* Hybrid modified Eddington-delta function metods, p633,Table1.
* W.E.Meador and W.R.Weaver, GAS,1980,v37,p.630
* W.J.Wiscombe and G.W. Grams, GAS,1976,v33,p2440
c      YLM0 = 2.
c      YLM2 = -3.*G*MU
c      YLM4 = 0.875*G**3*MU*(5.*MU**2-3.)
c      YLM6=-0.171875*G**5*MU*(15.-70.*MU**2+63.*MU**4)
c     YLM8=+0.073242*G**7*MU*(-35.+315.*MU**2-693.*MU**4
c    *+429.*MU**6)
c     YLM10=-0.008118*G**9*MU*(315.-4620.*MU**2+18018.*MU**4
c    *-25740.*MU**6+12155.*MU**8)
c     YLM12=0.003685*G**11*MU*(-693.+15015.*MU**2-90090.*MU**4
c    *+218790.*MU**6-230945.*MU**8+88179.*MU**10)
c      YLMS=YLM0+YLM2+YLM4+YLM6+YLM8+YLM10+YLM12
c      YLMS=0.25*YLMS
c      BETA0 = YLMS
c
c         amu1=1./1.7320508
c      YLM0 = 2.
c      YLM2 = -3.*G*amu1
c      YLM4 = 0.875*G**3*amu1*(5.*amu1**2-3.)
c      YLM6=-0.171875*G**5*amu1*(15.-70.*amu1**2+63.*amu1**4)
c     YLM8=+0.073242*G**7*amu1*(-35.+315.*amu1**2-693.*amu1**4
c    *+429.*amu1**6)
c     YLM10=-0.008118*G**9*amu1*(315.-4620.*amu1**2+18018.*amu1**4
c    *-25740.*amu1**6+12155.*amu1**8)
c     YLM12=0.003685*G**11*amu1*(-693.+15015.*amu1**2-90090.*amu1**4
c    *+218790.*amu1**6-230945.*amu1**8+88179.*amu1**10)
c      YLMS=YLM0+YLM2+YLM4+YLM6+YLM8+YLM10+YLM12
c      YLMS=0.25*YLMS
c      BETA1 = YLMS
c
c         BETAn = 0.25*(2. - 1.5*G-0.21875*G**3-0.085938*G**5
c    *-0.045776*G**7)


* Hybrid modified Eddington-delta function(Meador and Weaver,1980,JAS,37,630):
c         subd=4.*(1.-G*G*(1.-MU))
c         GAM1 = (7.-3.*G*G-OM*(4.+3.*G)+OM*G*G*(4.*BETA0+3.*G))/subd
c         GAM2 =-(1.-G*G-OM*(4.-3.*G)-OM*G*G*(4.*BETA0+3.*G-4.))/subd
c         GAM3 = BETA0
c         GAM4 = 1. - GAM3
c         mu1(i) = (1. - g*g*(1.- mu) )/(2. - g*g)

*****
* delta function  (Meador, and Weaver, 1980, JAS, 37, 630):
c         GAM1 = (1. - OM*(1. - beta0))/MU
c         GAM2 = OM*BETA0/MU
c         GAM3 = BETA0
c         GAM4 = 1. - GAM3
c         mu1(i) = mu
*****
* modified quadrature (Meador, and Weaver, 1980, JAS, 37, 630):
c         GAM1 = 1.7320508*(1. - OM*(1. - beta1))
c         GAM2 = 1.7320508*OM*beta1
c         GAM3 = BETA0
c         GAM4 = 1. - GAM3
c         mu1(i) = 1./sqrt(3.)

* hemispheric constant (Toon et al., 1989, JGR, 94, 16287):
c         GAM1 = 2.*(1. - OM*(1. - betan))
c         GAM2 = 2.*OM*BETAn
c         GAM3 = BETA0
c         GAM4 = 1. - GAM3
c         mu1(i) = 0.5

*****

* lambda = pg 16,290 equation 21
* big gamma = pg 16,290 equation 22
 
         lam(i) = sqrt(gam1*gam1 - gam2*gam2)
	 if( gam2.ne.0) then
          bgam(i) = (gam1 - lam(i))/gam2
	 else
	  bgam(i) = 0.
	 endif  

         expon = EXP(-lam(i)*taun(i))

* e1 - e4 = pg 16,292 equation 44
         
         e1(i) = 1. + bgam(i)*expon
         e2(i) = 1. - bgam(i)*expon
         e3(i) = bgam(i) + expon
         e4(i) = bgam(i) - expon
	 if(abs(e1(i)).gt.1e34.or.abs(e2(i)).gt.1e34.or.
     1     abs(e3(i)).gt.1e34.or.abs(e4(i).gt.1e34)) then
          write(*,*)' E data overflowed in ps2str, i,bgam,expon=',i,
     1	    bgam(i),expon
          write(*,*)' e =', e1(i),e2(i),e3(i),e4(i)
	  write(*,*)' lam, taun, gam1, gam2, gam3, gam4=',lam(i),
     1	  taun(i),gam1,gam2,gam3,gam4
          write(*,*)'g,om,mu, mu2=',g,om,mu,mu2(i)
	  stop 2001
	 endif 

* the following sets up for the C equations 23, and 24
* found on page 16,290
* prevent division by zero (if LAMBDA=1/MU, shift 1/MU^2 by EPS = 1.E-3
* which is approx equiv to shifting MU by 0.5*EPS* (MU)**3

         expon0 = EXP(-tausla(i-1))
         expon1 = EXP(-tausla(i))
          
         divisr = lam(i)*lam(i) - 1./(mu2(i)*mu2(i))
         temp = AMAX1(eps,abs(divisr))
         divisr = SIGN(temp,divisr)

         up = om*pifs*((gam1 - 1./mu2(i))*gam3 + gam4*gam2)/divisr
         dn = om*pifs*((gam1 + 1./mu2(i))*gam4 + gam2*gam3)/divisr
         
* cup and cdn are when tau is equal to zero
* cuptn and cdntn are when tau is equal to taun

         cup(i) = up*expon0
         cdn(i) = dn*expon0
         cuptn(i) = up*expon1
         cdntn(i) = dn*expon1
c	if(abs(cup(i)).gt.largest.or.abs(cdn(i)).gt.largest.or.
c     1	 abs(cuptn(i)).gt.largest.or.abs(cdntn(i)).gt.largest) then
c         print*,'something overflow in ps2str, i, nlayer=',i, nlayer
c	 print*,'cup(i),cdn(i),cuptn(i),cdntn(i)=',cup(i),cdn(i),
c     1	   cuptn(i),cdntn(i)
c         print*,'taun=',taun 
c	 print*,'dsdh=',dsdh(i,1:kz)
c	 print*,'omu=',omu
c	 print*,'gu=',gu
c         print*,'expon0, expon1, nid(i), tausla(i), tausla=',expon0, 
c     1	 expon1,nid(i), tausla(i),tausla
c	 print*,'up,dn,divisr,om,pifs,gam1,gam2,gam3,gam4,mu2(i)=',
c     1	  up,dn,divisr,om,pifs,gam1,gam2,gam3,gam4,mu2(i)
c         stop
c	endif  
   11 CONTINUE

***************** set up matrix ******
* ssfc = pg 16,292 equation 37  where pi Fs is one (unity).

      ssfc = rsfc*mu*EXP(-tausla(nlayer))*pifs

* MROWS = the number of rows in the matrix

      mrows = 2*nlayer     
      
* the following are from pg 16,292  equations 39 - 43.
* set up first row of matrix:

      i = 1
      a(1) = 0.
      b(1) = e1(i)
      d(1) = -e2(i)
      e(1) = fdn0 - cdn(i)

      row=1

* set up odd rows 3 thru (MROWS - 1):

      i = 0
      DO 20, row = 3, mrows - 1, 2
         i = i + 1
         a(row) = e2(i)*e3(i) - e4(i)*e1(i)
         b(row) = e1(i)*e1(i + 1) - e3(i)*e3(i + 1)
         d(row) = e3(i)*e4(i + 1) - e1(i)*e2(i + 1)
         e(row) = e3(i)*(cup(i + 1) - cuptn(i)) + 
     $        e1(i)*(cdntn(i) - cdn(i + 1))
   20 CONTINUE

* set up even rows 2 thru (MROWS - 2): 

      i = 0
      DO 30, row = 2, mrows - 2, 2
         i = i + 1
         a(row) = e2(i + 1)*e1(i) - e3(i)*e4(i + 1)
         b(row) = e2(i)*e2(i + 1) - e4(i)*e4(i + 1)
         d(row) = e1(i + 1)*e4(i + 1) - e2(i + 1)*e3(i + 1)
         e(row) = (cup(i + 1) - cuptn(i))*e2(i + 1) - 
     $        (cdn(i + 1) - cdntn(i))*e4(i + 1)
   30 CONTINUE

* set up last row of matrix at MROWS:

      row = mrows
      i = nlayer
      
      a(row) = e1(i) - rsfc*e3(i)
      b(row) = e2(i) - rsfc*e4(i)
      d(row) = 0.
      e(row) = ssfc - cuptn(i) + rsfc*cdntn(i)

* solve tri-diagonal matrix:

      CALL tridag(a, b, d, e, y, mrows)
c      do row=1,mrows     ! check y value
c       if(abs(y(row)).gt.100) then
c         print*,'Y overflowed after tridag in TUV',row
c	 print*,'a=',a
c	 print*,'b=',b
c	 print*,'d=',d
c	 print*,'e=',e
c	 print*,'y=',y
c	 stop
c	endif
c       enddo	 
**** unfold solution of matrix, compute output fluxes:

      row = 1 
      lev = 1
      j = 1
      
* the following equations are from pg 16,291  equations 31 & 32

      fdr(lev) = EXP( -tausla(0) )
      edr(lev) = mu * fdr(lev)
      edn(lev) = fdn0
      eup(lev) =  y(row)*e3(j) - y(row + 1)*e4(j) + cup(j)
      fdn(lev) = edn(lev)/mu1(lev)
      fup(lev) = eup(lev)/mu1(lev)

      DO 60, lev = 2, nlayer + 1
         fdr(lev) = EXP(-tausla(lev-1))
         edr(lev) =  mu *fdr(lev)
         edn(lev) =  y(row)*e3(j) + y(row + 1)*e4(j) + cdntn(j)
         eup(lev) =  y(row)*e1(j) + y(row + 1)*e2(j) + cuptn(j)
         fdn(lev) = edn(lev)/mu1(j)
         fup(lev) = eup(lev)/mu1(j)
c         if(fdr(lev).gt.100.or.fdn(lev).gt.100.or.fup(lev).gt.100) then
c	  print*,'fdr=',fdr
c	  print*,'fdn=',fdn
c	  print*,'fup=',fup
c	  print*,'mrows,row,lev,y(row),y(row+1)=',mrows,row,lev,y(row),
c     1	   y(row+1)
c	  print*,'j,e1(j),e2(j),e3(j),e4(j),cuptn(j),cdntn(j),mu1(j)=',
c     1	   j,e1(j),e2(j),e3(j),e4(j),cuptn(j),cdntn(j),mu1(j)
c          stop
c	 endif 
         row = row + 2
         j = j + 1
   60 CONTINUE
*_______________________________________________________________________

      RETURN
      END

************************************************************************

      SUBROUTINE tridag(a,b,c,r,u,n)
*_______________________________________________________________________
* solves tridiagonal system.  From Numerical Recipies, p. 40
*_______________________________________________________________________

      IMPLICIT NONE

* input:
      INTEGER n
      REAL a, b, c, r
      DIMENSION a(n),b(n),c(n),r(n)

* output:
      REAL u
      DIMENSION u(n)

* local:
      INTEGER j

      INCLUDE 'tuv.params'
      REAL bet, gam
      DIMENSION gam(2*kz)
*_______________________________________________________________________

      IF (b(1) .EQ. 0.) STOP 1001
      bet   = b(1)
      u(1) = r(1)/bet
      DO 11, j = 2, n   
         gam(j) = c(j - 1)/bet
         bet = b(j) - a(j)*gam(j)
         IF (bet .EQ. 0.) STOP 2002 
         u(j) = (r(j) - a(j)*u(j - 1))/bet
   11 CONTINUE
      DO 12, j = n - 1, 1, -1  
         u(j) = u(j) - gam(j + 1)*u(j + 1)
   12 CONTINUE
*_______________________________________________________________________

      RETURN
      END


c  Note: CDIR$ and CFPP$ comment lines are relevant only when running
c        on Cray computers.  They cause better optimization of loops
c        immediately following.


C      SUBROUTINE DISORT( dsdh, nid,
      SUBROUTINE PSNDO( dsdh, nid,
     &                   NLYR, DTAUC, SSALB, PMOM, TEMPER, WVNMLO,
     &                   WVNMHI, USRTAU, NTAU, UTAU, NSTR, USRANG, NUMU,
     &                   UMU, CWT, NPHI, PHI, IBCND, FBEAM, UMU0, PHI0,
     &                   FISOT, LAMBER, ALBEDO, HL, BTEMP, TTEMP, TEMIS,
     &                   DELTAM, PLANK, ONLYFL, ACCUR, PRNT, HEADER,
     &                   MAXCLY, MAXULV, MAXUMU, MAXCMU, MAXPHI, RFLDIR,
     &                   RFLDN, FLUP, DFDT, UAVG, UU, U0U, ALBMED,
     &                   TRNMED, uavgso, uavgup, uavgdn, 
     &                   sindir, sinup, sindn )

c Improved handling of numerical instabilities. Bernhard Mayer on 5/3/99.
c  disort seems to produce unstable results for certain combinations
c  of single scattering albedo and phase function. A temporary fix has been 
c  introduced to avoid this problem: The original instability check in 
c  UPBEAM fails on certain compiler/machine combinations (e.g., gcc/LINUX, 
c  or xlf/IBM RS6000). This check has therefore been replaced by a new one.
c  Whenever UPBEAM reports an instability, the single scattering albedo 
c  of the respective layer is changed by a small amount, and the 
c  calculation is repeated until numerically stable conditions are reached 
c  (all the necessary changes are confined to the new subroutine SOLVEC 
c  and the slighly changed subroutine UPBEAM). To check for potential 
c  instabilities, the variable 'RCOND' returned by SGECO is compared to
c  a machine-dependent constant, 'MINRCOND'. The value of this constant 
c  determines (a) if really all instabilities are caught; and (b) the 
c  amount by which the single scattering albedo has to be changed. The 
c  value of 'MINRCOND' is therefore a compromise between numerical 
c  stability on the one hand and uncertainties introduced by changing 
c  the atmospheric conditions and increased computational time on the 
c  other hand (an increase of MINRCOND will lead to the detection of
c  more potential numerical instabilities, and thus to an increase in 
c  computational time; by changing the atmospheric conditions, that is,
c  the single scattering albedo, the result might however be changed 
c  unfavourably, if the change is too large). From a limited number 
c  of experiments we found that 'MINRCOND = 5000. * R1MACH(4)' seems 
c  to be a good choice if high accuracy is required (more tests are 
c  definitely neccessary!). If an instability is encountered, a message 
c  is printed telling about neccessary changes to the single scattering 
c  albedo. This message may be switched off by setting 'DEBUG = .FALSE.' 
c  in subroutine SOLVEC. 
c
c
c modified to calculate sine-weighted intensities. Bernhard Mayer on 2/12/99.
c modified to handle some numerical instabilities. Chris Fischer on 1/22/99.
c modified by adding pseudo-spherical correction. Jun Zeng on 3/11/97.
c dsdh: slant path of direct beam through each layer crossed  
c       when travelling from the top of the atmosphere to layer i;    
c       dsdh(i,j), i = 0..nlyr, j = 1..nlyr;
c nid:  number of layers crossed by the direct beam when   
c       travelling from the top of the atmosphere to layer i; 
c       NID(i), i = 0..nlyr.
c uavgso, uvagup, and uvagdn are direct, downward diffuse, and upward
c diffuse actinic flux (mean intensity).
c u0u is the azimuthally averaged intensity, check DISORT.doc for details.
c *******************************************************************
c       Plane-parallel discrete ordinates radiative transfer program
c                      V E R S I O N    1.1
c             ( see DISORT.DOC for complete documentation )
c *******************************************************************


c +------------------------------------------------------------------+
c  Calling Tree (omitting calls to ERRMSG):
c  (routines in parentheses are not in this file)

c  DISORT-+-(R1MACH)
c         +-ZEROIT
c         +-CHEKIN-+-(WRTBAD)
c         |        +-(WRTDIM)
c         |        +-DREF
c         +-ZEROAL
c         +-SETDIS-+-QGAUSN (1)-+-(D1MACH)
c         +-PRTINP
c         +-LEPOLY see 2
c         +-SURFAC-+-QGAUSN see 1
c         |        +-LEPOLY see 2
c         |        +-ZEROIT
c         +-SOLEIG see 3
c         +-UPBEAM-+-(SGECO)
c         |        +-(SGESL)
c         +-TERPEV
c         +-TERPSO
c         +-SETMTX see 4
c         +-SOLVE0-+-ZEROIT
c         |        +-(SGBCO)
c         |        +-(SGBSL)
c         +-FLUXES--ZEROIT
c         +-PRAVIN
c         +-RATIO--(R1MACH)
c         +-PRTINT

c *** Intrinsic Functions used in DISORT package which take
c     non-negligible amount of time:

c    EXP :  Called by- ALBTRN, ALTRIN, CMPINT, FLUXES, SETDIS,
c                      SETMTX, SPALTR, USRINT, PLKAVG

c    SQRT : Called by- ASYMTX, LEPOLY, SOLEIG

c +-------------------------------------------------------------------+

c  Index conventions (for all DO-loops and all variable descriptions):

c     IU     :  for user polar angles

c  IQ,JQ,KQ  :  for computational polar angles ('quadrature angles')

c   IQ/2     :  for half the computational polar angles (just the ones
c               in either 0-90 degrees, or 90-180 degrees)

c     J      :  for user azimuthal angles

c     K,L    :  for Legendre expansion coefficients or, alternatively,
c               subscripts of associated Legendre polynomials

c     LU     :  for user levels

c     LC     :  for computational layers (each having a different
c               single-scatter albedo and/or phase function)

c    LEV     :  for computational levels

c    MAZIM   :  for azimuthal components in Fourier cosine expansion
c               of intensity and phase function

c +------------------------------------------------------------------+

c               I N T E R N A L    V A R I A B L E S

c   AMB(IQ/2,IQ/2)    First matrix factor in reduced eigenvalue problem
c                     of Eqs. SS(12), STWJ(8E)  (used only in SOLEIG)

c   APB(IQ/2,IQ/2)    Second matrix factor in reduced eigenvalue problem
c                     of Eqs. SS(12), STWJ(8E)  (used only in SOLEIG)

c   ARRAY(IQ,IQ)      Scratch matrix for SOLEIG, UPBEAM and UPISOT
c                     (see each subroutine for definition)

c   B()               Right-hand side vector of Eq. SC(5) going into
c                     SOLVE0,1;  returns as solution vector
c                     vector  L, the constants of integration

c   BDR(IQ/2,0:IQ/2)  Bottom-boundary bidirectional reflectivity for a
c                     given azimuthal component.  First index always
c                     refers to a computational angle.  Second index:
c                     if zero, refers to incident beam angle UMU0;
c                     if non-zero, refers to a computational angle.

c   BEM(IQ/2)         Bottom-boundary directional emissivity at compu-
c                     tational angles.

c   BPLANK            Intensity emitted from bottom boundary

c   CBAND()           Matrix of left-hand side of the linear system
c                     Eq. SC(5), scaled by Eq. SC(12);  in banded
c                     form required by LINPACK solution routines

c   CC(IQ,IQ)         C-sub-IJ in Eq. SS(5)

c   CMU(IQ)           Computational polar angles (Gaussian)

c   CWT(IQ)           Quadrature weights corresponding to CMU

c   DELM0             Kronecker delta, delta-sub-M0, where M = MAZIM
c                     is the number of the Fourier component in the
c                     azimuth cosine expansion

c   DITHER            Small quantity subtracted from single-scattering
c                     albedos of unity, in order to avoid using special
c                     case formulas;  prevents an eigenvalue of exactly
c                     zero from occurring, which would cause an
c                     immediate overflow

c   DTAUCP(LC)        Computational-layer optical depths (delta-M-scaled
c                     if DELTAM = TRUE, otherwise equal to DTAUC)

c   EMU(IU)           Bottom-boundary directional emissivity at user
c                     angles.

c   EVAL(IQ)          Temporary storage for eigenvalues of Eq. SS(12)

c   EVECC(IQ,IQ)      Complete eigenvectors of SS(7) on return from
c                     SOLEIG; stored permanently in  GC

c   EXPBEA(LC)        Transmission of direct beam in delta-M optical
c                     depth coordinates

c   FLYR(LC)          Truncated fraction in delta-M method

c   GL(K,LC)          Phase function Legendre polynomial expansion
c                     coefficients, calculated from PMOM by
c                     including single-scattering albedo, factor
c                     2K+1, and (if DELTAM=TRUE) the delta-M
c                     scaling

c   GC(IQ,IQ,LC)      Eigenvectors at polar quadrature angles,
c                     g  in Eq. SC(1)

c   GU(IU,IQ,LC)      Eigenvectors interpolated to user polar angles
c                     ( g  in Eqs. SC(3) and S1(8-9), i.e.
c                       G without the L factor )

c   HLPR()            Legendre coefficients of bottom bidirectional
c                     reflectivity (after inclusion of 2K+1 factor)

c   IPVT(LC*IQ)       Integer vector of pivot indices for LINPACK
c                     routines

c   KK(IQ,LC)         Eigenvalues of coeff. matrix in Eq. SS(7)

c   KCONV             Counter in azimuth convergence test

c   LAYRU(LU)         Computational layer in which user output level
c                     UTAU(LU) is located

c   LL(IQ,LC)         Constants of integration L in Eq. SC(1),
c                     obtained by solving scaled version of Eq. SC(5)

c   LYRCUT            TRUE, radiation is assumed zero below layer
c                     NCUT because of almost complete absorption

c   NAZ               Number of azimuthal components considered

c   NCUT              Computational layer number in which absorption
c                     optical depth first exceeds ABSCUT

c   OPRIM(LC)         Single scattering albedo after delta-M scaling

c   PASS1             TRUE on first entry, FALSE thereafter

c   PKAG(0:LC)        Integrated Planck function for internal emission

c   PSI(IQ)           Sum just after square bracket in  Eq. SD(9)

c   RMU(IU,0:IQ)      Bottom-boundary bidirectional reflectivity for a
c                     given azimuthal component.  First index always
c                     refers to a user angle.  Second index:
c                     if zero, refers to incident beam angle UMU0;
c                     if non-zero, refers to a computational angle.

c   TAUC(0:LC)        Cumulative optical depth (un-delta-M-scaled)

c   TAUCPR(0:LC)      Cumulative optical depth (delta-M-scaled if
c                     DELTAM = TRUE, otherwise equal to TAUC)

c   TPLANK            Intensity emitted from top boundary

c   UUM(IU,LU)        Expansion coefficients when the intensity
c                     (u-super-M) is expanded in Fourier cosine series
c                     in azimuth angle

c   U0C(IQ,LU)        Azimuthally-averaged intensity

c   UTAUPR(LU)        Optical depths of user output levels in delta-M
c                     coordinates;  equal to  UTAU(LU) if no delta-M

c   WK()              scratch array

c   XR0(LC)           X-sub-zero in expansion of thermal source func-
c                     tion preceding Eq. SS(14) (has no mu-dependence)

c   XR1(LC)           X-sub-one in expansion of thermal source func-
c                     tion;  see  Eqs. SS(14-16)

c   YLM0(L)           Normalized associated Legendre polynomial
c                     of subscript L at the beam angle (not saved
c                     as function of superscipt M)

c   YLMC(L,IQ)        Normalized associated Legendre polynomial
c                     of subscript L at the computational angles
c                     (not saved as function of superscipt M)

c   YLMU(L,IU)        Normalized associated Legendre polynomial
c                     of subscript L at the user angles
c                     (not saved as function of superscipt M)

c   Z()               scratch array used in  SOLVE0,1  to solve a
c                     linear system for the constants of integration

c   Z0(IQ)            Solution vectors Z-sub-zero of Eq. SS(16)

c   Z0U(IU,LC)        Z-sub-zero in Eq. SS(16) interpolated to user
c                     angles from an equation derived from SS(16)

c   Z1(IQ)            Solution vectors Z-sub-one  of Eq. SS(16)

c   Z1U(IU,LC)        Z-sub-one in Eq. SS(16) interpolated to user
c                     angles from an equation derived from SS(16)

c   ZBEAM(IU,LC)      Particular solution for beam source

c   ZJ(IQ)            Right-hand side vector  X-sub-zero in
c                     Eq. SS(19), also the solution vector
c                     Z-sub-zero after solving that system

c   ZZ(IQ,LC)         Permanent storage for the beam source vectors ZJ

c   ZPLK0(IQ,LC)      Permanent storage for the thermal source
c                     vectors  Z0  obtained by solving  Eq. SS(16)

c   ZPLK1(IQ,LC)      Permanent storage for the thermal source
c                     vectors  Z1  obtained by solving  Eq. SS(16)

c +-------------------------------------------------------------------+

c  LOCAL SYMBOLIC DIMENSIONS (have big effect on storage requirements):

c       MXCLY  = Max no. of computational layers
c       MXULV  = Max no. of output levels
c       MXCMU  = Max no. of computation polar angles
c       MXUMU  = Max no. of output polar angles
c       MXPHI  = Max no. of output azimuthal angles

c +-------------------------------------------------------------------+

      INCLUDE 'tuv.params'

c     .. Parameters ..

      INTEGER   MXCLY, MXULV, MXCMU, MXUMU, MXPHI, MI, MI9M2, NNLYRI
      PARAMETER ( MXCLY = 101, MXULV = 101, MXCMU = 32, MXUMU = 32,
     &          MXPHI = 3, MI = MXCMU / 2, MI9M2 = 9*MI - 2,
     &          NNLYRI = MXCMU*MXCLY )
c     ..
c     .. Scalar Arguments ..

      CHARACTER HEADER*127
      LOGICAL   DELTAM, LAMBER, ONLYFL, PLANK, USRANG, USRTAU
      INTEGER   IBCND, MAXCLY, MAXCMU, MAXPHI, MAXULV, MAXUMU, NLYR,
     &          NPHI, NSTR, NTAU, NUMU
      REAL      ACCUR, ALBEDO, BTEMP, FBEAM, FISOT, PHI0, TEMIS, TTEMP,
     &          UMU0, WVNMHI, WVNMLO

c     sherical geometry
      REAL dsdh(0:kz,kz)
      INTEGER nid(0:kz)
      REAL tausla(0:kz), tauslau(0:kz), mu2(0:kz)
c     ..
c     .. Array Arguments ..

      LOGICAL   PRNT( 7 )
      REAL      ALBMED( MAXUMU ), DFDT( MAXULV ), DTAUC( MAXCLY ),
     &          FLUP( MAXULV ), HL( 0:MAXCMU ), PHI( MAXPHI ),
     &          PMOM( 0:MAXCMU, MAXCLY ), RFLDIR( MAXULV ),
     &          RFLDN( MAXULV ), SSALB( MAXCLY ), TEMPER( 0:MAXCLY ),
     &          TRNMED( MAXUMU ), U0U( MAXUMU, MAXULV ), UAVG( MAXULV ),
     &          UMU( MAXUMU ), CWT( MAXCMU ), UTAU( MAXULV ),
     &          UU( MAXUMU, MAXULV, MAXPHI ), 
     &          uavgso( maxulv ), uavgup( maxulv ), uavgdn( maxulv ),
     &          sindir( maxulv ), sinup( maxulv ),  sindn ( maxulv )
c     ..
c     .. Local Scalars ..

      LOGICAL   COMPAR, LYRCUT, PASS1
      INTEGER   IQ, IU, J, KCONV, L, LC, LEV, LU, MAZIM, NAZ, NCOL,
     &          NCOS, NCUT, NN
      REAL      ANGCOS, AZERR, AZTERM, BPLANK, COSPHI, DELM0, DITHER,
     &          DUM, RPD, SGN, TPLANK
c     ..
c     .. Local Arrays ..

      INTEGER   IPVT( NNLYRI ), LAYRU( MXULV )

      REAL      AMB( MI, MI ), APB( MI, MI ), ARRAY( MXCMU, MXCMU ),
     &          B( NNLYRI ), BDR( MI, 0:MI ), BEM( MI ),
     &          CBAND( MI9M2, NNLYRI ), CC( MXCMU, MXCMU ),
     &          CMU( MXCMU ), DTAUCP( MXCLY ),
     &          EMU( MXUMU ), EVAL( MI ), EVECC( MXCMU, MXCMU ),
     &          EXPBEA( 0:MXCLY ), FLDIR( MXULV ), FLDN( MXULV ),
     &          FLYR( MXCLY ), GC( MXCMU, MXCMU, MXCLY ),
     &          GL( 0:MXCMU, MXCLY ), GU( MXUMU, MXCMU, MXCLY ),
     &          HLPR( 0:MXCMU ), KK( MXCMU, MXCLY ), LL( MXCMU, MXCLY ),
     &          OPRIM( MXCLY ), PHIRAD( MXPHI ), PKAG( 0:MXCLY ),
     &          PSI( MXCMU ), RMU( MXUMU, 0:MI ), TAUC( 0:MXCLY ),
     &          TAUCPR( 0:MXCLY ), U0C( MXCMU, MXULV ), UTAUPR( MXULV ),
     &          UUM( MXUMU, MXULV ), WK( MXCMU ), XR0( MXCLY ),
     &          XR1( MXCLY ), YLM0( 0:MXCMU ), YLMC( 0:MXCMU, MXCMU ),
     &          YLMU( 0:MXCMU, MXUMU ), Z( NNLYRI ), Z0( MXCMU ),
     &          Z0U( MXUMU, MXCLY ), Z1( MXCMU ), Z1U( MXUMU, MXCLY ),
     &          ZBEAM( MXUMU, MXCLY ), ZJ( MXCMU ),
     &          ZPLK0( MXCMU, MXCLY ), ZPLK1( MXCMU, MXCLY ),
     &          ZZ( MXCMU, MXCLY )

cgy added glsave and dgl to allow adjustable dimensioning in SOLVEC
      REAL GLSAVE( 0:MXCMU ), DGL( 0:MXCMU )

      DOUBLE PRECISION AAD( MI, MI ), EVALD( MI ), EVECCD( MI, MI ),
     &                 WKD( MXCMU )
c     ..
c     .. External Functions ..

      REAL      PLKAVG, R1MACH, RATIO
      EXTERNAL  PLKAVG, R1MACH, RATIO
c     ..
c     .. External Subroutines ..

      EXTERNAL  CHEKIN, FLUXES, LEPOLY, PRAVIN, PRTINP,
     &          PRTINT, SETDIS, SETMTX, SOLEIG, SOLVE0, SURFAC,
     &          UPBEAM, ZEROAL, ZEROIT
c     ..
c     .. Intrinsic Functions ..

      INTRINSIC ABS, ASIN, COS, LEN, MAX
c     ..
      SAVE      PASS1, DITHER, RPD
      DATA      PASS1 / .TRUE. /



      IF( PASS1 ) THEN

         DITHER = 10.*R1MACH( 4 )

c                            ** Must dither more on Cray (14-digit prec)

         IF( DITHER.LT.1.E-10 ) DITHER = 10.*DITHER

         RPD  = PI / 180.0
         PASS1 = .FALSE.
      END IF
 
   10 CONTINUE

c                                  ** Calculate cumulative optical depth
c                                     and dither single-scatter albedo
c                                     to improve numerical behavior of
c                                     eigenvalue/vector computation
      CALL ZEROIT( TAUC, MXCLY + 1 )

      DO 20 LC = 1, NLYR

         IF( SSALB( LC ).EQ.1.0 ) SSALB( LC ) = 1.0 - DITHER
         TAUC( LC ) = TAUC( LC - 1 ) + DTAUC( LC )

   20 CONTINUE
c                                ** Check input dimensions and variables

      CALL CHEKIN( NLYR, DTAUC, SSALB, PMOM, TEMPER, WVNMLO, WVNMHI,
     &             USRTAU, NTAU, UTAU, NSTR, USRANG, NUMU, UMU, NPHI,
     &             PHI, IBCND, FBEAM, UMU0, PHI0, FISOT, LAMBER, ALBEDO,
     &             HL, BTEMP, TTEMP, TEMIS, PLANK, ONLYFL, ACCUR, TAUC,
     &             MAXCLY, MAXULV, MAXUMU, MAXCMU, MAXPHI, MXCLY, MXULV,
     &             MXUMU, MXCMU, MXPHI )

c                                 ** Zero internal and output arrays

      CALL  ZEROAL( MXCLY, EXPBEA(1), FLYR, OPRIM, TAUCPR(1), XR0, XR1,
     $              MXCMU, CMU, CWT, PSI, WK, Z0, Z1, ZJ,
     $              MXCMU+1, HLPR, YLM0,
     $              MXCMU**2, ARRAY, CC, EVECC,
     $              (MXCMU+1)*MXCLY, GL,
     $              (MXCMU+1)*MXCMU, YLMC,
     $              (MXCMU+1)*MXUMU, YLMU,
     $              MXCMU*MXCLY, KK, LL, ZZ, ZPLK0, ZPLK1,
     $              MXCMU**2*MXCLY, GC,
     $              MXULV, LAYRU, UTAUPR,
     $              MXUMU*MXCMU*MXCLY, GU,
     $              MXUMU*MXCLY, Z0U, Z1U, ZBEAM,
     $              MI, EVAL,
     $              MI**2, AMB, APB,
     $              NNLYRI, IPVT, Z,
     $              MAXULV, RFLDIR, RFLDN, FLUP, UAVG, DFDT,
     $              MAXUMU, ALBMED, TRNMED,
     $              MAXUMU*MAXULV, U0U,
     $              MAXUMU*MAXULV*MAXPHI, UU )

c                                 ** Perform various setup operations

      CALL SETDIS( dsdh, nid, tausla, tauslau, mu2,
     &             CMU, CWT, DELTAM, DTAUC, DTAUCP, EXPBEA, FBEAM, FLYR,
     &             GL, HL, HLPR, IBCND, LAMBER, LAYRU, LYRCUT, MAXUMU,
     &             MAXCMU, MXCMU, NCUT, NLYR, NTAU, NN, NSTR, PLANK,
     &             NUMU, ONLYFL, OPRIM, PMOM, SSALB, TAUC, TAUCPR, UTAU,
     &             UTAUPR, UMU, UMU0, USRTAU, USRANG )

c                                 ** Print input information
      IF ( PRNT(1) )
     $     CALL PRTINP( NLYR, DTAUC, DTAUCP, SSALB, PMOM, TEMPER,
     $                  WVNMLO, WVNMHI, NTAU, UTAU, NSTR, NUMU, UMU,
     $                  NPHI, PHI, IBCND, FBEAM, UMU0, PHI0, FISOT,
     $                  LAMBER, ALBEDO, HL, BTEMP, TTEMP, TEMIS,
     $                  DELTAM, PLANK, ONLYFL, ACCUR, FLYR, LYRCUT,
     $                  OPRIM, TAUC, TAUCPR, MAXCMU, PRNT(7) )

c                              ** Handle special case for getting albedo
c                                 and transmissivity of medium for many
c                                 beam angles at once
c                                   ** Calculate Planck functions

         BPLANK = 0.0
         TPLANK = 0.0
         CALL ZEROIT( PKAG, MXCLY + 1 )

c ========  BEGIN LOOP TO SUM AZIMUTHAL COMPONENTS OF INTENSITY  =======
c           (EQ STWJ 5)

      KCONV  = 0
      NAZ  = NSTR - 1
c                                    ** Azimuth-independent case

      IF( FBEAM.EQ.0.0 .OR. ( 1.- UMU0 ).LT.1.E-5 .OR. ONLYFL .OR.
     &      ( NUMU.EQ.1 .AND. ( 1.- UMU(1) ).LT.1.E-5 ) )
     &   NAZ = 0

      DO 160 MAZIM = 0, NAZ

         IF( MAZIM.EQ.0 ) DELM0  = 1.0
         IF( MAZIM.GT.0 ) DELM0  = 0.0

c                             ** Get normalized associated Legendre
c                                polynomials for
c                                (a) incident beam angle cosine
c                                (b) computational and user polar angle
c                                    cosines
         IF( FBEAM.GT.0.0 ) THEN

            NCOS   = 1
            ANGCOS = -UMU0

            CALL LEPOLY( NCOS, MAZIM, MXCMU, NSTR - 1, ANGCOS, YLM0 )

         END IF


         IF( .NOT.ONLYFL .AND. USRANG )
     &       CALL LEPOLY( NUMU, MAZIM, MXCMU, NSTR-1, UMU, YLMU )

         CALL LEPOLY( NN, MAZIM, MXCMU, NSTR-1, CMU, YLMC )

c                       ** Get normalized associated Legendre polys.
c                          with negative arguments from those with
c                          positive arguments; Dave/Armstrong Eq. (15)
         SGN  = - 1.0

         DO 50 L = MAZIM, NSTR - 1

            SGN  = - SGN

            DO 40 IQ = NN + 1, NSTR
               YLMC( L, IQ ) = SGN*YLMC( L, IQ - NN )
   40       CONTINUE

   50    CONTINUE
c                                 ** Specify users bottom reflectivity
c                                    and emissivity properties
      IF ( .NOT.LYRCUT )
     $   CALL  SURFAC( ALBEDO, DELM0, FBEAM, HLPR, LAMBER,
     $                 MI, MAZIM, MXCMU, MXUMU, NN, NUMU, NSTR, ONLYFL,
     $                 UMU, USRANG, YLM0, YLMC, YLMU, BDR, EMU, BEM,
     $                 RMU )


c ===================  BEGIN LOOP ON COMPUTATIONAL LAYERS  =============

         DO 60 LC = 1, NCUT

            CALL SOLVEC( AMB, APB, ARRAY, CMU, CWT, GL( 0,LC ), MI,
     &           MAZIM, MXCMU, NN, NSTR, YLM0, YLMC, CC, 
     &           EVECC, EVAL, KK( 1,LC ), GC( 1,1,LC ), AAD, EVECCD, 
     &           EVALD, WK, WKD, DELM0, FBEAM, IPVT, PI, UMU0,
     &           ZJ, ZZ(1,LC), OPRIM(LC), LC, DITHER, mu2(lc),
     &           glsave, dgl)
cgy added glsave and dgl to call to allow adjustable dimensioning

 60      CONTINUE


c ===================  END LOOP ON COMPUTATIONAL LAYERS  ===============


c                      ** Set coefficient matrix of equations combining
c                         boundary and layer interface conditions

         CALL SETMTX( BDR, CBAND, CMU, CWT, DELM0, DTAUCP, GC, KK,
     &                LAMBER, LYRCUT, MI, MI9M2, MXCMU, NCOL, NCUT,
     &                NNLYRI, NN, NSTR, TAUCPR, WK )

c                      ** Solve for constants of integration in homo-
c                         geneous solution (general boundary conditions)

         CALL SOLVE0( B, BDR, BEM, BPLANK, CBAND, CMU, CWT, EXPBEA,
     &                FBEAM, FISOT, IPVT, LAMBER, LL, LYRCUT, MAZIM, MI,
     &                MI9M2, MXCMU, NCOL, NCUT, NN, NSTR, NNLYRI, PI,
     &                TPLANK, TAUCPR, UMU0, Z, ZZ, ZPLK0, ZPLK1 )

c                                  ** Compute upward and downward fluxes

      IF ( MAZIM.EQ.0 )
     $     CALL FLUXES( tausla, tauslau,
     $                  CMU, CWT, FBEAM, GC, KK, LAYRU, LL, LYRCUT,
     $                  MAXULV, MXCMU, MXULV, NCUT, NN, NSTR, NTAU,
     $                  PI, PRNT, SSALB, TAUCPR, UMU0, UTAU, UTAUPR,
     $                  XR0, XR1, ZZ, ZPLK0, ZPLK1, DFDT, FLUP,
     $                  FLDN, FLDIR, RFLDIR, RFLDN, UAVG, U0C,
     $                  uavgso, uavgup, uavgdn,
     $                  sindir, sinup, sindn)

         IF( ONLYFL ) THEN

            IF( MAXUMU.GE.NSTR ) THEN
c                                     ** Save azimuthal-avg intensities
c                                        at quadrature angles
               DO 80 LU = 1, NTAU

                  DO 70 IQ = 1, NSTR
                     U0U( IQ, LU ) = U0C( IQ, LU )
   70             CONTINUE

   80          CONTINUE

            END IF

            GO TO  170

         END IF


         CALL ZEROIT( UUM, MXUMU*MXULV )

         IF( MAZIM.EQ.0 ) THEN
c                               ** Save azimuthally averaged intensities

            DO 110 LU = 1, NTAU

               DO 100 IU = 1, NUMU
                  U0U( IU, LU ) = UUM( IU, LU )

                  DO 90 J = 1, NPHI
                     UU( IU, LU, J ) = UUM( IU, LU )
   90             CONTINUE

  100          CONTINUE

  110       CONTINUE
c                              ** Print azimuthally averaged intensities
c                                 at user angles

            IF( PRNT( 4 ) ) CALL PRAVIN( UMU, NUMU, MAXUMU, UTAU, NTAU,
     &                                   U0U )
            IF( NAZ.GT.0 ) THEN

               CALL ZEROIT( PHIRAD, MXPHI )
               DO 120 J = 1, NPHI
                  PHIRAD( J ) = RPD*( PHI( J ) - PHI0 )
  120          CONTINUE

            END IF


         ELSE
c                                ** Increment intensity by current
c                                   azimuthal component (Fourier
c                                   cosine series);  Eq SD(2)
            AZERR  = 0.0

            DO 150 J = 1, NPHI

               COSPHI = COS( MAZIM*PHIRAD( J ) )

               DO 140 LU = 1, NTAU

                  DO 130 IU = 1, NUMU
                     AZTERM = UUM( IU, LU )*COSPHI
                     UU( IU, LU, J ) = UU( IU, LU, J ) + AZTERM
                     AZERR = MAX( AZERR,
     &                       RATIO( ABS(AZTERM), ABS(UU(IU,LU,J)) ) )
  130             CONTINUE

  140          CONTINUE

  150       CONTINUE

            IF( AZERR.LE.ACCUR ) KCONV  = KCONV + 1

            IF( KCONV.GE.2 ) GO TO  170

         END IF

  160 CONTINUE

c ===================  END LOOP ON AZIMUTHAL COMPONENTS  ===============


c                                          ** Print intensities
  170 CONTINUE
      IF( PRNT( 5 ) .AND. .NOT.ONLYFL ) CALL PRTINT( UU, UTAU, NTAU,
     &    UMU, NUMU, PHI, NPHI, MAXULV, MAXUMU )

      END

      SUBROUTINE ASYMTX( AA, EVEC, EVAL, M, IA, IEVEC, IER, WKD, AAD,
     &                   EVECD, EVALD )

c    =======  D O U B L E    P R E C I S I O N    V E R S I O N  ======

c       Solves eigenfunction problem for real asymmetric matrix
c       for which it is known a priori that the eigenvalues are real.

c       This is an adaptation of a subroutine EIGRF in the IMSL
c       library to use real instead of complex arithmetic, accounting
c       for the known fact that the eigenvalues and eigenvectors in
c       the discrete ordinate solution are real.  Other changes include
c       putting all the called subroutines in-line, deleting the
c       performance index calculation, updating many DO-loops
c       to Fortran77, and in calculating the machine precision
c       TOL instead of specifying it in a data statement.

c       EIGRF is based primarily on EISPACK routines.  The matrix is
c       first balanced using the Parlett-Reinsch algorithm.  Then
c       the Martin-Wilkinson algorithm is applied.

c       References:
c          Dongarra, J. and C. Moler, EISPACK -- A Package for Solving
c             Matrix Eigenvalue Problems, in Cowell, ed., 1984:
c             Sources and Development of Mathematical Software,
c             Prentice-Hall, Englewood Cliffs, NJ
c         Parlett and Reinsch, 1969: Balancing a Matrix for Calculation
c             of Eigenvalues and Eigenvectors, Num. Math. 13, 293-304
c         Wilkinson, J., 1965: The Algebraic Eigenvalue Problem,
c             Clarendon Press, Oxford

c   I N P U T    V A R I A B L E S:

c       AA    :  input asymmetric matrix, destroyed after solved
c        M    :  order of  AA
c       IA    :  first dimension of  AA
c    IEVEC    :  first dimension of  EVEC

c   O U T P U T    V A R I A B L E S:

c       EVEC  :  (unnormalized) eigenvectors of  AA
c                   ( column J corresponds to EVAL(J) )

c       EVAL  :  (unordered) eigenvalues of AA ( dimension at least M )

c       IER   :  if .NE. 0, signals that EVAL(IER) failed to converge;
c                   in that case eigenvalues IER+1,IER+2,...,M  are
c                   correct but eigenvalues 1,...,IER are set to zero.

c   S C R A T C H   V A R I A B L E S:

c       WKD   :  work area ( dimension at least 2*M )
c       AAD   :  double precision stand-in for AA
c       EVECD :  double precision stand-in for EVEC
c       EVALD :  double precision stand-in for EVAL

c   Called by- SOLEIG
c   Calls- D1MACH, ERRMSG
c +-------------------------------------------------------------------+

c     .. Scalar Arguments ..

      INTEGER   IA, IER, IEVEC, M
c     ..
c     .. Array Arguments ..

      REAL      AA( IA, M ), EVAL( M ), EVEC( IEVEC, M )
      DOUBLE PRECISION AAD( IA, M ), EVALD( M ), EVECD( IA, M ),
     &                 WKD( * )
c     ..
c     .. Local Scalars ..

      LOGICAL   NOCONV, NOTLAS
      INTEGER   I, II, IN, J, K, KA, KKK, L, LB, LLL, N, N1, N2
      DOUBLE PRECISION C1, C2, C3, C4, C5, C6, COL, DISCRI, F, G, H,
     &                 ONE, P, Q, R, REPL, RNORM, ROW, S, SCALE, SGN, T,
     &                 TOL, UU, VV, W, X, Y, Z, ZERO
c     ..
c     .. External Functions ..

      DOUBLE PRECISION D1MACH
      EXTERNAL  D1MACH
c     ..
c     .. External Subroutines ..

      EXTERNAL  ERRMSG
c     ..
c     .. Intrinsic Functions ..

      INTRINSIC ABS, DBLE, MIN, SIGN, SQRT
c     ..
      DATA      C1 / 0.4375D0 / , C2 / 0.5D0 / , C3 / 0.75D0 / ,
     &          C4 / 0.95D0 / , C5 / 16.D0 / , C6 / 256.D0 / ,
     &          ZERO / 0.D0 / , ONE / 1.D0 /


      IER  = 0
      TOL  = D1MACH( 4 )

      IF( M.LT.1 .OR. IA.LT.M .OR. IEVEC.LT.M )
     &    CALL ERRMSG( 'ASYMTX--bad input variable(s)', .TRUE. )

c                           ** Handle 1x1 and 2x2 special cases

      IF( M.EQ.1 ) THEN

         EVAL( 1 )    = AA( 1, 1 )
         EVEC( 1, 1 ) = 1.0
         RETURN

      ELSE IF( M.EQ.2 ) THEN

         DISCRI = ( AA( 1,1 ) - AA( 2,2 ) )**2 +
     &              4.*AA( 1, 2 )*AA( 2, 1 )

         IF( DISCRI.LT.0.0 )
     &       CALL ERRMSG( 'ASYMTX--complex evals in 2x2 case',.TRUE. )

         SGN  = 1.0

         IF( AA( 1,1 ).LT.AA( 2,2 ) ) SGN  = - 1.0

         EVAL( 1 ) = 0.5*( AA( 1,1 ) + AA( 2,2 ) + SGN*SQRT( DISCRI ) )
         EVAL( 2 ) = 0.5*( AA( 1,1 ) + AA( 2,2 ) - SGN*SQRT( DISCRI ) )
         EVEC( 1, 1 ) = 1.0
         EVEC( 2, 2 ) = 1.0

         IF( AA( 1,1 ).EQ.AA( 2,2 ) .AND.
     &       ( AA( 2,1 ).EQ.0.0 .OR. AA( 1,2 ).EQ.0.0 ) ) THEN

            RNORM  = ABS( AA( 1,1 ) ) + ABS( AA( 1,2 ) ) +
     &               ABS( AA( 2,1 ) ) + ABS( AA( 2,2 ) )
            W  = TOL*RNORM
            EVEC( 2, 1 ) =   AA( 2, 1 ) / W
            EVEC( 1, 2 ) = - AA( 1, 2 ) / W

         ELSE

            EVEC( 2, 1 ) = AA( 2, 1 ) / ( EVAL( 1 ) - AA( 2,2 ) )
            EVEC( 1, 2 ) = AA( 1, 2 ) / ( EVAL( 2 ) - AA( 1,1 ) )

         END IF

         RETURN

      END IF
c                               ** Put s.p. matrix into d.p. matrix
      DO 20 J = 1, M

         DO 10 K = 1, M
            AAD( J, K ) = DBLE( AA( J,K ) )
   10    CONTINUE

   20 CONTINUE

c                                ** Initialize output variables
      IER  = 0

      DO 40 I = 1, M
         EVALD( I ) = ZERO

         DO 30 J = 1, M
            EVECD( I, J ) = ZERO
   30    CONTINUE

         EVECD( I, I ) = ONE
   40 CONTINUE

c                  ** Balance the input matrix and reduce its norm by
c                     diagonal similarity transformation stored in WK;
c                     then search for rows isolating an eigenvalue
c                     and push them down
      RNORM  = ZERO
      L  = 1
      K  = M

   50 CONTINUE
      KKK  = K

      DO 90 J = KKK, 1, -1

         ROW  = ZERO

         DO 60 I = 1, K

            IF( I.NE.J ) ROW  = ROW + ABS( AAD( J,I ) )

   60    CONTINUE

         IF( ROW.EQ.ZERO ) THEN

            WKD( K ) = J

            IF( J.NE.K ) THEN

               DO 70 I = 1, K
                  REPL        = AAD( I, J )
                  AAD( I, J ) = AAD( I, K )
                  AAD( I, K ) = REPL
   70          CONTINUE

               DO 80 I = L, M
                  REPL        = AAD( J, I )
                  AAD( J, I ) = AAD( K, I )
                  AAD( K, I ) = REPL
   80          CONTINUE

            END IF

            K  = K - 1
            GO TO  50

         END IF

   90 CONTINUE
c                                ** Search for columns isolating an
c                                   eigenvalue and push them left
  100 CONTINUE
      LLL  = L

      DO 140 J = LLL, K

         COL  = ZERO

         DO 110 I = L, K

            IF( I.NE.J ) COL  = COL + ABS( AAD( I,J ) )

  110    CONTINUE

         IF( COL.EQ.ZERO ) THEN

            WKD( L ) = J

            IF( J.NE.L ) THEN

               DO 120 I = 1, K
                  REPL        = AAD( I, J )
                  AAD( I, J ) = AAD( I, L )
                  AAD( I, L ) = REPL
  120          CONTINUE

               DO 130 I = L, M
                  REPL        = AAD( J, I )
                  AAD( J, I ) = AAD( L, I )
                  AAD( L, I ) = REPL
  130          CONTINUE

            END IF

            L  = L + 1
            GO TO  100

         END IF

  140 CONTINUE

c                           ** Balance the submatrix in rows L through K
      DO 150 I = L, K
         WKD( I ) = ONE
  150 CONTINUE

  160 CONTINUE
      NOCONV = .FALSE.

      DO 220 I = L, K

         COL  = ZERO
         ROW  = ZERO

         DO 170 J = L, K

            IF( J.NE.I ) THEN

               COL  = COL + ABS( AAD( J,I ) )
               ROW  = ROW + ABS( AAD( I,J ) )

            END IF

  170    CONTINUE

         F  = ONE
         G  = ROW / C5
         H  = COL + ROW

  180    CONTINUE
         IF( COL.LT.G ) THEN

            F    = F*C5
            COL  = COL*C6
            GO TO  180

         END IF

         G  = ROW*C5

  190    CONTINUE
         IF( COL.GE.G ) THEN

            F    = F / C5
            COL  = COL / C6
            GO TO  190

         END IF
c                                                ** Now balance
         IF( ( COL + ROW ) / F.LT.C4*H ) THEN

            WKD( I ) = WKD( I )*F
            NOCONV = .TRUE.

            DO 200 J = L, M
               AAD( I, J ) = AAD( I, J ) / F
  200       CONTINUE

            DO 210 J = 1, K
               AAD( J, I ) = AAD( J, I )*F
  210       CONTINUE

         END IF

  220 CONTINUE


      IF( NOCONV ) GO TO  160
c                                   ** Is A already in Hessenberg form?
      IF( K-1 .LT. L+1 ) GO TO  370

c                                   ** Transfer A to a Hessenberg form
      DO 310 N = L + 1, K - 1

         H  = ZERO
         WKD( N + M ) = ZERO
         SCALE  = ZERO
c                                                 ** Scale column
         DO 230 I = N, K
            SCALE  = SCALE + ABS( AAD( I,N - 1 ) )
  230    CONTINUE

         IF( SCALE.NE.ZERO ) THEN

            DO 240 I = K, N, -1
               WKD( I + M ) = AAD( I, N - 1 ) / SCALE
               H  = H + WKD( I + M )**2
  240       CONTINUE

            G    = - SIGN( SQRT( H ), WKD( N + M ) )
            H    = H - WKD( N + M )*G
            WKD( N + M ) = WKD( N + M ) - G
c                                            ** Form (I-(U*UT)/H)*A
            DO 270 J = N, M

               F  = ZERO

               DO 250 I = K, N, -1
                  F  = F + WKD( I + M )*AAD( I, J )
  250          CONTINUE

               DO 260 I = N, K
                  AAD( I, J ) = AAD( I, J ) - WKD( I + M )*F / H
  260          CONTINUE

  270       CONTINUE
c                                    ** Form (I-(U*UT)/H)*A*(I-(U*UT)/H)
            DO 300 I = 1, K

               F  = ZERO

               DO 280 J = K, N, -1
                  F  = F + WKD( J + M )*AAD( I, J )
  280          CONTINUE

               DO 290 J = N, K
                  AAD( I, J ) = AAD( I, J ) - WKD( J + M )*F / H
  290          CONTINUE

  300       CONTINUE

            WKD( N + M ) = SCALE*WKD( N + M )
            AAD( N, N - 1 ) = SCALE*G

         END IF

  310 CONTINUE


      DO 360 N = K - 2, L, -1

         N1   = N + 1
         N2   = N + 2
         F  = AAD( N + 1, N )

         IF( F.NE.ZERO ) THEN

            F  = F*WKD( N + 1 + M )

            DO 320 I = N + 2, K
               WKD( I + M ) = AAD( I, N )
  320       CONTINUE

            IF( N + 1.LE.K ) THEN

               DO 350 J = 1, M

                  G  = ZERO

                  DO 330 I = N + 1, K
                     G  = G + WKD( I + M )*EVECD( I, J )
  330             CONTINUE

                  G  = G / F

                  DO 340 I = N + 1, K
                     EVECD( I, J ) = EVECD( I, J ) + G*WKD( I + M )
  340             CONTINUE

  350          CONTINUE

            END IF

         END IF

  360 CONTINUE


  370 CONTINUE

      N  = 1

      DO 390 I = 1, M

         DO 380 J = N, M
            RNORM  = RNORM + ABS( AAD( I,J ) )
  380    CONTINUE

         N  = I

         IF( I.LT.L .OR. I.GT.K ) EVALD( I ) = AAD( I, I )

  390 CONTINUE

      N  = K
      T  = ZERO

c                                      ** Search for next eigenvalues
  400 CONTINUE
      IF( N.LT.L ) GO TO  550

      IN  = 0
      N1  = N - 1
      N2  = N - 2
c                          ** Look for single small sub-diagonal element
  410 CONTINUE

      DO 420 I = L, N
         LB  = N + L - I

         IF( LB.EQ.L ) GO TO  430

         S  = ABS( AAD( LB - 1,LB - 1 ) ) + ABS( AAD( LB,LB ) )

         IF( S.EQ.ZERO ) S  = RNORM

         IF( ABS( AAD( LB, LB-1 ) ).LE. TOL*S ) GO TO  430

  420 CONTINUE


  430 CONTINUE
      X  = AAD( N, N )

      IF( LB.EQ.N ) THEN
c                                        ** One eigenvalue found
         AAD( N, N ) = X + T
         EVALD( N ) = AAD( N, N )
         N  = N1
         GO TO  400

      END IF

C next line has been included to avoid run time error caused by xlf

      IF ( ( N1.LE.0 ).OR.( N.LE.0 ) ) THEN
        WRITE(0,*) 'Subscript out of bounds in ASYMTX'
        STOP 9999
      ENDIF

      Y  = AAD( N1, N1 )
      W  = AAD( N, N1 )*AAD( N1, N )

      IF( LB.EQ.N1 ) THEN
c                                        ** Two eigenvalues found
         P  = ( Y - X )*C2
         Q  = P**2 + W
         Z  = SQRT( ABS( Q ) )
         AAD( N, N ) = X + T
         X  = AAD( N, N )
         AAD( N1, N1 ) = Y + T
c                                        ** Real pair
         Z  = P + SIGN( Z, P )
         EVALD( N1 ) = X + Z
         EVALD( N ) = EVALD( N1 )

         IF( Z.NE.ZERO ) EVALD( N ) = X - W / Z

         X  = AAD( N, N1 )
c                                  ** Employ scale factor in case
c                                     X and Z are very small
         R  = SQRT( X*X + Z*Z )
         P  = X / R
         Q  = Z / R
c                                             ** Row modification
         DO 440 J = N1, M
            Z  = AAD( N1, J )
            AAD( N1, J ) = Q*Z + P*AAD( N, J )
            AAD( N, J ) = Q*AAD( N, J ) - P*Z
  440    CONTINUE
c                                             ** Column modification
         DO 450 I = 1, N
            Z  = AAD( I, N1 )
            AAD( I, N1 ) = Q*Z + P*AAD( I, N )
            AAD( I, N ) = Q*AAD( I, N ) - P*Z
  450    CONTINUE
c                                          ** Accumulate transformations
         DO 460 I = L, K
            Z  = EVECD( I, N1 )
            EVECD( I, N1 ) = Q*Z + P*EVECD( I, N )
            EVECD( I, N ) = Q*EVECD( I, N ) - P*Z
  460    CONTINUE

         N  = N2
         GO TO  400

      END IF


      IF( IN.EQ.30 ) THEN

c                    ** No convergence after 30 iterations; set error
c                       indicator to the index of the current eigenvalue
         IER  = N
         GO TO  700

      END IF
c                                                  ** Form shift
      IF( IN.EQ.10 .OR. IN.EQ.20 ) THEN

         T  = T + X

         DO 470 I = L, N
            AAD( I, I ) = AAD( I, I ) - X
  470    CONTINUE

         S  = ABS( AAD( N,N1 ) ) + ABS( AAD( N1,N2 ) )
         X  = C3*S
         Y  = X
         W  = -C1*S**2

      END IF


      IN  = IN + 1

c                ** Look for two consecutive small sub-diagonal elements

C inhibit vectorization by CF77, as this will cause a run time error

CDIR$ NEXTSCALAR
      DO 480 J = LB, N2
         I  = N2 + LB - J
         Z  = AAD( I, I )
         R  = X - Z
         S  = Y - Z
         P  = ( R*S - W ) / AAD( I + 1, I ) + AAD( I, I + 1 )
         Q  = AAD( I + 1, I + 1 ) - Z - R - S
         R  = AAD( I + 2, I + 1 )
         S  = ABS( P ) + ABS( Q ) + ABS( R )
         P  = P / S
         Q  = Q / S
         R  = R / S

         IF( I.EQ.LB ) GO TO  490

         UU   = ABS( AAD( I, I-1 ) )*( ABS( Q ) + ABS( R ) )
         VV   = ABS( P ) * ( ABS( AAD( I-1, I-1 ) ) + ABS( Z ) +
     &                       ABS( AAD( I+1, I+1 ) ) )

         IF( UU .LE. TOL*VV ) GO TO  490

  480 CONTINUE

  490 CONTINUE
      AAD( I+2, I ) = ZERO

c                      ** fpp vectorization of this loop triggers
c                         array bounds errors, so inhibit
CFPP$ NOVECTOR L
      DO 500 J = I + 3, N
         AAD( J, J - 2 ) = ZERO
         AAD( J, J - 3 ) = ZERO
  500 CONTINUE

c             ** Double QR step involving rows K to N and columns M to N

      DO 540 KA = I, N1

         NOTLAS = KA.NE.N1

         IF( KA.EQ.I ) THEN

            S  = SIGN( SQRT( P*P + Q*Q + R*R ), P )

            IF( LB.NE.I ) AAD( KA, KA - 1 ) = -AAD( KA, KA - 1 )

         ELSE

            P  = AAD( KA, KA - 1 )
            Q  = AAD( KA + 1, KA - 1 )
            R  = ZERO

            IF( NOTLAS ) R  = AAD( KA + 2, KA - 1 )

            X  = ABS( P ) + ABS( Q ) + ABS( R )

            IF( X.EQ.ZERO ) GO TO  540

            P  = P / X
            Q  = Q / X
            R  = R / X
            S  = SIGN( SQRT( P*P + Q*Q + R*R ), P )
            AAD( KA, KA - 1 ) = -S*X

         END IF

         P  = P + S
         X  = P / S
         Y  = Q / S
         Z  = R / S
         Q  = Q / P
         R  = R / P
c                                              ** Row modification
         DO 510 J = KA, M

            P  = AAD( KA, J ) + Q*AAD( KA + 1, J )

            IF( NOTLAS ) THEN

               P  = P + R*AAD( KA + 2, J )
               AAD( KA + 2, J ) = AAD( KA + 2, J ) - P*Z

            END IF

            AAD( KA + 1, J ) = AAD( KA + 1, J ) - P*Y
            AAD( KA, J ) = AAD( KA, J ) - P*X
  510    CONTINUE
c                                                 ** Column modification
         DO 520 II = 1, MIN( N, KA + 3 )

            P  = X*AAD( II, KA ) + Y*AAD( II, KA + 1 )

            IF( NOTLAS ) THEN

               P  = P + Z*AAD( II, KA + 2 )
               AAD( II, KA + 2 ) = AAD( II, KA + 2 ) - P*R

            END IF

            AAD( II, KA + 1 ) = AAD( II, KA + 1 ) - P*Q
            AAD( II, KA ) = AAD( II, KA ) - P
  520    CONTINUE
c                                          ** Accumulate transformations
         DO 530 II = L, K

            P  = X*EVECD( II, KA ) + Y*EVECD( II, KA + 1 )

            IF( NOTLAS ) THEN

               P  = P + Z*EVECD( II, KA + 2 )
               EVECD( II, KA + 2 ) = EVECD( II, KA + 2 ) - P*R

            END IF

            EVECD( II, KA + 1 ) = EVECD( II, KA + 1 ) - P*Q
            EVECD( II, KA ) = EVECD( II, KA ) - P
  530    CONTINUE

  540 CONTINUE

      GO TO  410
c                     ** All evals found, now backsubstitute real vector
  550 CONTINUE

      IF( RNORM.NE.ZERO ) THEN

         DO 580 N = M, 1, -1
            N2   = N
            AAD( N, N ) = ONE

            DO 570 I = N - 1, 1, -1
               W  = AAD( I, I ) - EVALD( N )

               IF( W.EQ.ZERO ) W  = TOL*RNORM

               R  = AAD( I, N )

               DO 560 J = N2, N - 1
                  R  = R + AAD( I, J )*AAD( J, N )
  560          CONTINUE

               AAD( I, N ) = -R / W
               N2   = I
  570       CONTINUE

  580    CONTINUE
c                      ** End backsubstitution vectors of isolated evals
         DO 600 I = 1, M

            IF( I.LT.L .OR. I.GT.K ) THEN

               DO 590 J = I, M
                  EVECD( I, J ) = AAD( I, J )
  590          CONTINUE

            END IF

  600    CONTINUE
c                                   ** Multiply by transformation matrix
         IF( K.NE.0 ) THEN

            DO 630 J = M, L, -1

               DO 620 I = L, K
                  Z  = ZERO

                  DO 610 N = L, MIN( J, K )
                     Z  = Z + EVECD( I, N )*AAD( N, J )
  610             CONTINUE

                  EVECD( I, J ) = Z
  620          CONTINUE

  630       CONTINUE

         END IF

      END IF


      DO 650 I = L, K

         DO 640 J = 1, M
            EVECD( I, J ) = EVECD( I, J )*WKD( I )
  640    CONTINUE
  650 CONTINUE

c                           ** Interchange rows if permutations occurred
      DO 670 I = L-1, 1, -1

         J  = WKD( I )

         IF( I.NE.J ) THEN

            DO 660 N = 1, M
               REPL   = EVECD( I, N )
               EVECD( I, N ) = EVECD( J, N )
               EVECD( J, N ) = REPL
  660       CONTINUE

         END IF

  670 CONTINUE


      DO 690 I = K + 1, M

         J  = WKD( I )

         IF( I.NE.J ) THEN

            DO 680 N = 1, M
               REPL   = EVECD( I, N )
               EVECD( I, N ) = EVECD( J, N )
               EVECD( J, N ) = REPL
  680       CONTINUE

         END IF

  690 CONTINUE

c                         ** Put results into output arrays
  700 CONTINUE

      DO 720 J = 1, M

         EVAL( J ) = EVALD( J )

         DO 710 K = 1, M
            EVEC( J, K ) = EVECD( J, K )
  710    CONTINUE

  720 CONTINUE


      END

      SUBROUTINE CHEKIN( NLYR, DTAUC, SSALB, PMOM, TEMPER, WVNMLO,
     &                   WVNMHI, USRTAU, NTAU, UTAU, NSTR, USRANG, NUMU,
     &                   UMU, NPHI, PHI, IBCND, FBEAM, UMU0, PHI0,
     &                   FISOT, LAMBER, ALBEDO, HL, BTEMP, TTEMP, TEMIS,
     &                   PLANK, ONLYFL, ACCUR, TAUC, MAXCLY, MAXULV,
     &                   MAXUMU, MAXCMU, MAXPHI, MXCLY, MXULV, MXUMU,
     &                   MXCMU, MXPHI )

c           Checks the input dimensions and variables

c   Calls- WRTBAD, WRTDIM, DREF, ERRMSG
c   Called by- DISORT
c --------------------------------------------------------------------

c     .. Scalar Arguments ..

      LOGICAL   LAMBER, ONLYFL, PLANK, USRANG, USRTAU
      INTEGER   IBCND, MAXCLY, MAXCMU, MAXPHI, MAXULV, MAXUMU, MXCLY,
     &          MXCMU, MXPHI, MXULV, MXUMU, NLYR, NPHI, NSTR, NTAU, NUMU
      REAL      ACCUR, ALBEDO, BTEMP, FBEAM, FISOT, PHI0, TEMIS, TTEMP,
     &          UMU0, WVNMHI, WVNMLO
c     ..
c     .. Array Arguments ..

      REAL      DTAUC( MAXCLY ), HL( 0:MAXCMU ), PHI( MAXPHI ),
     &          PMOM( 0:MAXCMU, MAXCLY ), SSALB( MAXCLY ),
     &          TAUC( 0:MXCLY ), TEMPER( 0:MAXCLY ), UMU( MAXUMU ),
     &          UTAU( MAXULV )
c     ..
c     .. Local Scalars ..

      LOGICAL   INPERR
      INTEGER   IRMU, IU, J, K, LC, LU
      REAL      FLXALB, RMU
c     ..
c     .. External Functions ..

      LOGICAL   WRTBAD, WRTDIM
      REAL      DREF
      EXTERNAL  WRTBAD, WRTDIM, DREF
c     ..
c     .. External Subroutines ..

      EXTERNAL  ERRMSG
c     ..
c     .. Intrinsic Functions ..

      INTRINSIC ABS, MOD
c     ..


      INPERR = .FALSE.

      IF( NLYR.LT.1 ) INPERR = WRTBAD( 'NLYR' )

      IF( NLYR.GT.MAXCLY ) INPERR = WRTBAD( 'MAXCLY' )

      DO 20 LC = 1, NLYR

         IF( DTAUC( LC ).LT.0.0 ) INPERR = WRTBAD( 'DTAUC' )

         IF( SSALB( LC ).LT.0.0 .OR. SSALB( LC ).GT.1.0 )
     &       INPERR = WRTBAD( 'SSALB' )

         IF( PLANK .AND. IBCND.NE.1 ) THEN

            IF( LC.EQ.1 .AND. TEMPER( 0 ).LT.0.0 )
     &          INPERR = WRTBAD( 'TEMPER' )

            IF( TEMPER( LC ).LT.0.0 ) INPERR = WRTBAD( 'TEMPER' )

         END IF

         DO 10 K = 0, NSTR

            IF( PMOM( K,LC ).LT.-1.0 .OR. PMOM( K,LC ).GT.1.0 )
     &          INPERR = WRTBAD( 'PMOM' )

   10    CONTINUE

   20 CONTINUE


      IF( IBCND.EQ.1 ) THEN

         IF( MAXULV.LT.2 ) INPERR = WRTBAD( 'MAXULV' )

      ELSE IF( USRTAU ) THEN

         IF( NTAU.LT.1 ) INPERR = WRTBAD( 'NTAU' )

         IF( MAXULV.LT.NTAU ) INPERR = WRTBAD( 'MAXULV' )

         DO 30 LU = 1, NTAU

            IF( ABS( UTAU( LU ) - TAUC( NLYR ) ).LE. 1.E-4 )
     &          UTAU( LU ) = TAUC( NLYR )

            IF( UTAU( LU ).LT.0.0 .OR. UTAU( LU ).GT. TAUC( NLYR ) )
     &          INPERR = WRTBAD( 'UTAU' )

   30    CONTINUE

      ELSE

         IF( MAXULV.LT.NLYR + 1 ) INPERR = WRTBAD( 'MAXULV' )

      END IF


      IF( NSTR.LT.2 .OR. MOD( NSTR,2 ).NE.0 ) INPERR = WRTBAD( 'NSTR' )

c     IF( NSTR.EQ.2 )
c    &    CALL ERRMSG( 'CHEKIN--2 streams not recommended;'//
c    &                 ' use specialized 2-stream code instead',.False.)

      IF( NSTR.GT.MAXCMU ) INPERR = WRTBAD( 'MAXCMU' )

      IF( USRANG ) THEN

         IF( NUMU.LT.0 ) INPERR = WRTBAD( 'NUMU' )

         IF( .NOT.ONLYFL .AND. NUMU.EQ.0 ) INPERR = WRTBAD( 'NUMU' )

         IF( NUMU.GT.MAXUMU ) INPERR = WRTBAD( 'MAXUMU' )

         IF( IBCND.EQ.1 .AND. 2*NUMU.GT.MAXUMU )
     &       INPERR = WRTBAD( 'MAXUMU' )

         DO 40 IU = 1, NUMU

            IF( UMU( IU ).LT.-1.0 .OR. UMU( IU ).GT.1.0 .OR.
     &          UMU( IU ).EQ.0.0 ) INPERR = WRTBAD( 'UMU' )

            IF( IBCND.EQ.1 .AND. UMU( IU ).LT.0.0 )
     &          INPERR = WRTBAD( 'UMU' )

            IF( IU.GT.1 ) THEN

               IF( UMU( IU ).LT.UMU( IU-1 ) ) INPERR = WRTBAD( 'UMU' )

            END IF

   40    CONTINUE

      ELSE

         IF( MAXUMU.LT.NSTR ) INPERR = WRTBAD( 'MAXUMU' )

      END IF


      IF( .NOT.ONLYFL .AND. IBCND.NE.1 ) THEN

         IF( NPHI.LE.0 ) INPERR = WRTBAD( 'NPHI' )

         IF( NPHI.GT.MAXPHI ) INPERR = WRTBAD( 'MAXPHI' )

         DO 50 J = 1, NPHI

            IF( PHI( J ).LT.0.0 .OR. PHI( J ).GT.360.0 )
     &          INPERR = WRTBAD( 'PHI' )

   50    CONTINUE

      END IF


      IF( IBCND.LT.0 .OR. IBCND.GT.1 ) INPERR = WRTBAD( 'IBCND' )

      IF( IBCND.EQ.0 ) THEN

         IF( FBEAM.LT.0.0 ) INPERR = WRTBAD( 'FBEAM' )

         IF( FBEAM.GT.0.0 .AND. abs(UMU0).GT.1.0 )
     &       INPERR = WRTBAD( 'UMU0' )

         IF( FBEAM.GT.0.0 .AND. ( PHI0.LT.0.0 .OR.PHI0.GT.360.0 ) )
     &       INPERR = WRTBAD( 'PHI0' )

         IF( FISOT.LT.0.0 ) INPERR = WRTBAD( 'FISOT' )

         IF( LAMBER ) THEN

            IF( ALBEDO.LT.0.0 .OR. ALBEDO.GT.1.0 )
     &          INPERR = WRTBAD( 'ALBEDO' )

         ELSE
c                    ** Make sure flux albedo at dense mesh of incident
c                       angles does not assume unphysical values

            DO 60 IRMU = 0, 100
               RMU  = IRMU*0.01
               FLXALB = DREF( RMU, HL, NSTR )

               IF( FLXALB.LT.0.0 .OR. FLXALB.GT.1.0 )
     &             INPERR = WRTBAD( 'HL' )

   60       CONTINUE

         END IF


      ELSE IF( IBCND.EQ.1 ) THEN

         IF( ALBEDO.LT.0.0 .OR. ALBEDO.GT.1.0 )
     &       INPERR = WRTBAD( 'ALBEDO' )

      END IF


      IF( PLANK .AND. IBCND.NE.1 ) THEN

         IF( WVNMLO.LT.0.0 .OR. WVNMHI.LE.WVNMLO )
     &       INPERR = WRTBAD( 'WVNMLO,HI' )

         IF( TEMIS.LT.0.0 .OR. TEMIS.GT.1.0 ) INPERR = WRTBAD( 'TEMIS' )

         IF( BTEMP.LT.0.0 ) INPERR = WRTBAD( 'BTEMP' )

         IF( TTEMP.LT.0.0 ) INPERR = WRTBAD( 'TTEMP' )

      END IF


      IF( ACCUR.LT.0.0 .OR. ACCUR.GT.1.E-2 ) INPERR = WRTBAD( 'ACCUR' )

      IF( MXCLY.LT.NLYR ) INPERR = WRTDIM( 'MXCLY', NLYR )

      IF( IBCND.NE.1 ) THEN

         IF( USRTAU .AND. MXULV.LT.NTAU )
     &       INPERR = WRTDIM( 'MXULV',NTAU )

         IF( .NOT.USRTAU .AND. MXULV .LT. NLYR + 1 )
     &       INPERR = WRTDIM( 'MXULV', NLYR + 1 )

      ELSE

         IF( MXULV.LT.2 ) INPERR = WRTDIM( 'MXULV', 2 )

      END IF

      IF( MXCMU.LT.NSTR ) INPERR = WRTDIM( 'MXCMU', NSTR )

      IF( USRANG .AND. MXUMU.LT.NUMU ) INPERR = WRTDIM( 'MXUMU', NUMU )

      IF( USRANG .AND. IBCND.EQ.1 .AND.MXUMU.LT.2*NUMU )
     &    INPERR = WRTDIM( 'MXUMU', NUMU )

      IF( .NOT.USRANG .AND. MXUMU.LT.NSTR )
     &    INPERR = WRTDIM( 'MXUMU', NSTR )

      IF( .NOT.ONLYFL .AND. IBCND.NE.1 .AND. MXPHI.LT.NPHI )
     &    INPERR = WRTDIM( 'MXPHI', NPHI )

      IF( INPERR )
     &    CALL ERRMSG( 'DISORT--input and/or dimension errors',.True.)

      IF( PLANK ) THEN

         DO 70 LC = 1, NLYR

            IF( ABS( TEMPER( LC ) - TEMPER( LC-1 ) ).GT. 20.0 )
     &          CALL ERRMSG('CHEKIN--vertical temperature step may'
     &                      // ' be too large for good accuracy',
     &                      .False.)
   70    CONTINUE

      END IF

      END

      SUBROUTINE FLUXES( tausla, tauslau,
     &                   CMU, CWT, FBEAM, GC, KK, LAYRU, LL, LYRCUT,
     &                   MAXULV, MXCMU, MXULV, NCUT, NN, NSTR, NTAU, PI,
     &                   PRNT, SSALB, TAUCPR, UMU0, UTAU, UTAUPR, XR0,
     &                   XR1, ZZ, ZPLK0, ZPLK1, DFDT, FLUP, FLDN, FLDIR,
     &                   RFLDIR, RFLDN, UAVG, U0C,
     &                   uavgso, uavgup, uavgdn,
     $                   sindir, sinup, sindn)

c       Calculates the radiative fluxes, mean intensity, and flux
c       derivative with respect to optical depth from the m=0 intensity
c       components (the azimuthally-averaged intensity)

c    I N P U T     V A R I A B L E S:

c       CMU      :  Abscissae for Gauss quadrature over angle cosine
c       CWT      :  Weights for Gauss quadrature over angle cosine
c       GC       :  Eigenvectors at polar quadrature angles, SC(1)
c       KK       :  Eigenvalues of coeff. matrix in Eq. SS(7)
c       LAYRU    :  Layer number of user level UTAU
c       LL       :  Constants of integration in Eq. SC(1), obtained
c                     by solving scaled version of Eq. SC(5);
c                     exponential term of Eq. SC(12) not included
c       LYRCUT   :  Logical flag for truncation of comput. layer
c       NN       :  Order of double-Gauss quadrature (NSTR/2)
c       NCUT     :  Number of computational layer where absorption
c                     optical depth exceeds ABSCUT
c       TAUCPR   :  Cumulative optical depth (delta-M-scaled)
c       UTAUPR   :  Optical depths of user output levels in delta-M
c                     coordinates;  equal to UTAU if no delta-M
c       XR0      :  Expansion of thermal source function in Eq. SS(14)
c       XR1      :  Expansion of thermal source function Eqs. SS(16)
c       ZZ       :  Beam source vectors in Eq. SS(19)
c       ZPLK0    :  Thermal source vectors Z0, by solving Eq. SS(16)
c       ZPLK1    :  Thermal source vectors Z1, by solving Eq. SS(16)
c       (remainder are DISORT input variables)


c                   O U T P U T     V A R I A B L E S:

c       U0C      :  Azimuthally averaged intensities
c                   ( at polar quadrature angles )
c       (RFLDIR, RFLDN, FLUP, DFDT, UAVG are DISORT output variables)


c                   I N T E R N A L       V A R I A B L E S:

c       DIRINT   :  Direct intensity attenuated
c       FDNTOT   :  Total downward flux (direct + diffuse)
c       FLDIR    :  Direct-beam flux (delta-M scaled)
c       FLDN     :  Diffuse down-flux (delta-M scaled)
c       FNET     :  Net flux (total-down - diffuse-up)
c       FACT     :  EXP( - UTAUPR / UMU0 )
c       PLSORC   :  Planck source function (thermal)
c       ZINT     :  Intensity of m = 0 case, in Eq. SC(1)

c   Called by- DISORT
c   Calls- ZEROIT
c +-------------------------------------------------------------------+

c     .. Scalar Arguments ..

      LOGICAL   LYRCUT
      INTEGER   MAXULV, MXCMU, MXULV, NCUT, NN, NSTR, NTAU
      REAL      FBEAM, PI, UMU0
c     ..
c     .. Array Arguments ..

      LOGICAL   PRNT( * )
      INTEGER   LAYRU( MXULV )
      REAL      CMU( MXCMU ), CWT( MXCMU ), DFDT( MAXULV ),
     &          FLDIR( MXULV ), FLDN( MXULV ), FLUP( MAXULV ),
     &          GC( MXCMU, MXCMU, * ), KK( MXCMU, * ), LL( MXCMU, * ),
     &          RFLDIR( MAXULV ), RFLDN( MAXULV ), SSALB( * ),
     &          TAUCPR( 0:* ), U0C( MXCMU, MXULV ), UAVG( MAXULV ),
     &          UTAU( MAXULV ), UTAUPR( MXULV ), XR0( * ), XR1( * ),
     &          ZPLK0( MXCMU, * ), ZPLK1( MXCMU, * ), ZZ( MXCMU, * ),
     &          uavgso(*),uavgup(*), uavgdn(*),
     &          sindir(*),sinup(*), sindn(*)
      REAL tausla(0:*), tauslau(0:*)
c     ..
c     .. Local Scalars ..

      INTEGER   IQ, JQ, LU, LYU
      REAL      ANG1, ANG2, DIRINT, FACT, FDNTOT, FNET, PLSORC, ZINT
c     ..
c     .. External Subroutines ..

      EXTERNAL  ZEROIT
c     ..
c     .. Intrinsic Functions ..

      INTRINSIC ACOS, EXP
c     ..


      IF( PRNT( 2 ) ) WRITE( *, 9000 )
c                                          ** Zero DISORT output arrays
      CALL ZEROIT( U0C, MXULV*MXCMU )
      CALL ZEROIT( FLDIR, MXULV )
      CALL ZEROIT( FLDN, MXULV )
      call  zeroit( uavgso,   maxulv )
      call  zeroit( uavgup,   maxulv )
      call  zeroit( uavgdn,   maxulv )
      call  zeroit( sindir,   maxulv )
      call  zeroit( sinup,    maxulv )
      call  zeroit( sindn,    maxulv )

c                                        ** Loop over user levels
      DO 80 LU = 1, NTAU

         LYU  = LAYRU( LU )

         IF( LYRCUT .AND. LYU.GT.NCUT ) THEN
c                                                ** No radiation reaches
c                                                ** this level
            FDNTOT = 0.0
            FNET   = 0.0
            PLSORC = 0.0
            GO TO  70

         END IF

         IF( FBEAM.GT.0.0 ) THEN
 
            FACT  = EXP( - tausla(LU-1) )
            DIRINT       = FBEAM*FACT
            FLDIR( LU )  = UMU0*( FBEAM*FACT )
            RFLDIR( LU ) = UMU0*FBEAM * EXP( -tauslau(lu-1) )
            sindir( LU ) = SQRT(1.-UMU0*UMU0)*FBEAM * 
     $                     EXP( -tauslau(lu-1) )

         ELSE

            DIRINT       = 0.0
            FLDIR( LU )  = 0.0
            RFLDIR( LU ) = 0.0
            sindir( LU ) = 0.0

         END IF


         DO 30 IQ = 1, NN

            ZINT   = 0.0

            DO 10 JQ = 1, NN
               ZINT   = ZINT + GC( IQ, JQ, LYU )*LL( JQ, LYU )*
     &                  EXP( -KK( JQ,LYU )*( UTAUPR( LU ) -
     &                  TAUCPR( LYU ) ) )
   10       CONTINUE

            DO 20 JQ = NN + 1, NSTR
               ZINT   = ZINT + GC( IQ, JQ, LYU )*LL( JQ, LYU )*
     &                  EXP( -KK( JQ,LYU )*( UTAUPR( LU ) -
     &                  TAUCPR( LYU - 1 ) ) )
   20       CONTINUE

            U0C( IQ, LU ) = ZINT

            IF( FBEAM.GT.0.0 ) U0C( IQ, LU ) = ZINT + ZZ( IQ, LYU )*FACT

            U0C( IQ, LU ) = U0C( IQ, LU ) + ZPLK0( IQ, LYU ) +
     &                      ZPLK1( IQ, LYU )*UTAUPR( LU )
            UAVG( LU ) = UAVG( LU ) + CWT( NN + 1 - IQ )*U0C( IQ, LU )
            uavgdn(lu) = uavgdn(lu) + cwt(nn+1-iq) * u0c( iq,lu )
            sindn(lu)  = sindn(lu)  + cwt(nn+1-iq) * 
     &                   SQRT(1.-CMU(NN+1-IQ)*CMU(NN+1-IQ))*
     &                   U0C( IQ, LU )
            FLDN( LU ) = FLDN( LU ) + CWT( NN + 1 - IQ )*
     &                   CMU( NN + 1 - IQ )*U0C( IQ, LU )
   30    CONTINUE


         DO 60 IQ = NN + 1, NSTR

            ZINT   = 0.0

            DO 40 JQ = 1, NN
               ZINT   = ZINT + GC( IQ, JQ, LYU )*LL( JQ, LYU )*
     &                  EXP( -KK( JQ,LYU )*( UTAUPR( LU ) -
     &                  TAUCPR( LYU ) ) )
   40       CONTINUE

            DO 50 JQ = NN + 1, NSTR
               ZINT   = ZINT + GC( IQ, JQ, LYU )*LL( JQ, LYU )*
     &                  EXP( -KK( JQ,LYU )*( UTAUPR( LU ) -
     &                  TAUCPR( LYU - 1 ) ) )
   50       CONTINUE

            U0C( IQ, LU ) = ZINT

            IF( FBEAM.GT.0.0 ) U0C( IQ, LU ) = ZINT + ZZ( IQ, LYU )*FACT

            U0C( IQ, LU ) = U0C( IQ, LU ) + ZPLK0( IQ, LYU ) +
     &                      ZPLK1( IQ, LYU )*UTAUPR( LU )
            UAVG( LU ) = UAVG( LU ) + CWT( IQ - NN )*U0C( IQ, LU )
            uavgup(lu) = uavgup(lu) + cwt(iq-nn) * u0c( iq,lu )
            sinup (lu) = sinup(lu)  + cwt(iq-nn) * 
     &                   SQRT(1.-CMU(IQ-NN)*CMU(IQ-NN))*
     &                   U0C( IQ, LU )
            FLUP( LU ) = FLUP( LU ) + CWT( IQ - NN )*CMU( IQ - NN )*
     &                   U0C( IQ, LU )
   60    CONTINUE


         FLUP( LU )  = 2.*PI*FLUP( LU )
         FLDN( LU )  = 2.*PI*FLDN( LU )
         FDNTOT      = FLDN( LU ) + FLDIR( LU )
         FNET        = FDNTOT - FLUP( LU )
         RFLDN( LU ) = FDNTOT - RFLDIR( LU )
         UAVG( LU )  = ( 2.*PI*UAVG( LU ) + DIRINT ) / ( 4.*PI )
         uavgso( lu ) = dirint / (4.*pi)
         uavgup( lu ) = (2.0 * pi * uavgup(lu) )/ (4.*pi)
         uavgdn( lu)  = (2.0 * pi * uavgdn(lu) )/ (4.*pi)
         sindn ( lu ) = 2.*PI*sindn ( LU )
         sinup ( lu ) = 2.*PI*sinup ( LU )

         PLSORC      = XR0( LYU ) + XR1( LYU )*UTAUPR( LU )
         DFDT( LU )  = ( 1.- SSALB( LYU ) ) * 4.*PI *
     &                 ( UAVG( LU ) - PLSORC )

   70    CONTINUE
         IF( PRNT( 2 ) ) WRITE( *, FMT = 9010 ) UTAU( LU ), LYU,
     &       RFLDIR( LU ), RFLDN( LU ), FDNTOT, FLUP( LU ), FNET,
     &       UAVG( LU ), PLSORC, DFDT( LU )

   80 CONTINUE


      IF( PRNT( 3 ) ) THEN

         WRITE( *, FMT = 9020 )

         DO 100 LU = 1, NTAU

            WRITE( *, FMT = 9030 ) UTAU( LU )

            DO 90 IQ = 1, NN
               ANG1   = 180./ PI* ACOS( CMU( 2*NN - IQ + 1 ) )
               ANG2   = 180./ PI* ACOS( CMU( IQ ) )
               WRITE( *, 9040 ) ANG1, CMU(2*NN-IQ+1), U0C(IQ,LU),
     $                          ANG2, CMU(IQ),        U0C(IQ+NN,LU)
   90       CONTINUE

  100    CONTINUE

      END IF


 9000 FORMAT( //, 21X,
     $ '<----------------------- FLUXES ----------------------->', /,
     $ '   Optical  Compu    Downward    Downward    Downward     ',
     $ ' Upward                    Mean      Planck   d(Net Flux)', /,
     $ '     Depth  Layer      Direct     Diffuse       Total     ',
     $ 'Diffuse         Net   Intensity      Source   / d(Op Dep)', / )
 9010 FORMAT( F10.4, I7, 1P, 7E12.3, E14.3 )
 9020 FORMAT( / , / , ' ******** AZIMUTHALLY AVERAGED INTENSITIES',
     &      ' ( at polar quadrature angles ) *******' )
 9030 FORMAT( /, ' Optical depth =', F10.4, //,
     $  '     Angle (deg)   cos(Angle)     Intensity',
     $  '     Angle (deg)   cos(Angle)     Intensity' )
 9040 FORMAT( 2( 0P,F16.4,F13.5,1P,E14.3 ) )

      END

      SUBROUTINE LEPOLY( NMU, M, MAXMU, TWONM1, MU, YLM )

c       Computes the normalized associated Legendre polynomial,
c       defined in terms of the associated Legendre polynomial
c       Plm = P-sub-l-super-m as

c             Ylm(MU) = sqrt( (l-m)!/(l+m)! ) * Plm(MU)

c       for fixed order m and all degrees from l = m to TWONM1.
c       When m.GT.0, assumes that Y-sub(m-1)-super(m-1) is available
c       from a prior call to the routine.

c       REFERENCE: Dave, J.V. and B.H. Armstrong, Computations of
c                  High-Order Associated Legendre Polynomials,
c                  J. Quant. Spectrosc. Radiat. Transfer 10,
c                  557-562, 1970.  (hereafter D/A)

c       METHOD: Varying degree recurrence relationship.

c       NOTE 1: The D/A formulas are transformed by
c               setting  M = n-1; L = k-1.
c       NOTE 2: Assumes that routine is called first with  M = 0,
c               then with  M = 1, etc. up to  M = TWONM1.
c       NOTE 3: Loops are written in such a way as to vectorize.

c  I N P U T     V A R I A B L E S:

c       NMU    :  Number of arguments of YLM
c       M      :  Order of YLM
c       MAXMU  :  First dimension of YLM
c       TWONM1 :  Max degree of YLM
c       MU(i)  :  Arguments of YLM (i = 1 to NMU)

c       If M.GT.0, YLM(M-1,i) for i = 1 to NMU is assumed to exist
c       from a prior call.

c  O U T P U T     V A R I A B L E:

c       YLM(l,i) :  l = M to TWONM1, normalized associated Legendre
c                   polynomials evaluated at argument MU(i)

c   Called by- DISORT, ALBTRN, SURFAC
c   Calls- ERRMSG
c +-------------------------------------------------------------------+

c     .. Parameters ..

      INTEGER   MAXSQT
      PARAMETER ( MAXSQT = 1000 )
c     ..
c     .. Scalar Arguments ..

      INTEGER   M, MAXMU, NMU, TWONM1
c     ..
c     .. Array Arguments ..

      REAL      MU( * ), YLM( 0:MAXMU, * )
c     ..
c     .. Local Scalars ..

      LOGICAL   PASS1
      INTEGER   I, L, NS
      REAL      TMP1, TMP2
c     ..
c     .. Local Arrays ..

      REAL      SQT( MAXSQT )
c     ..
c     .. External Subroutines ..

      EXTERNAL  ERRMSG
c     ..
c     .. Intrinsic Functions ..

      INTRINSIC FLOAT, SQRT
c     ..
      SAVE      SQT, PASS1
      DATA      PASS1 / .TRUE. /


      IF( PASS1 ) THEN

         PASS1  = .FALSE.

         DO 10 NS = 1, MAXSQT
            SQT( NS ) = SQRT( FLOAT( NS ) )
   10    CONTINUE

      END IF

      IF( 2*TWONM1.GT.MAXSQT )
     &    CALL ERRMSG('LEPOLY--need to increase param MAXSQT',.True.)


      IF( M.EQ.0 ) THEN
c                             ** Upward recurrence for ordinary
c                                Legendre polynomials
         DO 20 I = 1, NMU
            YLM( 0, I ) = 1.0
            YLM( 1, I ) = MU( I )
   20    CONTINUE


         DO 40 L = 2, TWONM1

            DO 30 I = 1, NMU
               YLM( L, I ) = ( ( 2*L - 1 )*MU( I )*YLM( L - 1, I ) -
     &                         ( L - 1 )*YLM( L - 2, I ) ) / L
   30       CONTINUE

   40    CONTINUE


      ELSE

         DO 50 I = 1, NMU
c                               ** Y-sub-m-super-m; derived from
c                               ** D/A Eqs. (11,12)

            YLM( M, I ) = - SQT( 2*M - 1 ) / SQT( 2*M )*
     &                      SQRT( 1.- MU(I)**2 )*YLM( M - 1, I )

c                              ** Y-sub-(m+1)-super-m; derived from
c                              ** D/A Eqs.(13,14) using Eqs.(11,12)

            YLM( M + 1, I ) = SQT( 2*M + 1 )*MU( I )*YLM( M, I )

   50    CONTINUE

c                                   ** Upward recurrence; D/A EQ.(10)
         DO 70 L = M + 2, TWONM1

            TMP1  = SQT( L - M )*SQT( L + M )
            TMP2  = SQT( L - M - 1 )*SQT( L + M - 1 )

            DO 60 I = 1, NMU
               YLM( L, I ) = ( ( 2*L - 1 )*MU( I )*YLM( L-1, I ) -
     &                         TMP2*YLM( L-2, I ) ) / TMP1
   60       CONTINUE

   70    CONTINUE

      END IF


      END

      SUBROUTINE PRAVIN( UMU, NUMU, MAXUMU, UTAU, NTAU, U0U )

c        Print azimuthally averaged intensities at user angles

c   Called by- DISORT

c     LENFMT   Max number of polar angle cosines UMU that can be
c                printed on one line, as set in FORMAT statement
c --------------------------------------------------------------------

c     .. Scalar Arguments ..

      INTEGER   MAXUMU, NTAU, NUMU
c     ..
c     .. Array Arguments ..

      REAL      U0U( MAXUMU, NTAU ), UMU( NUMU ), UTAU( NTAU )
c     ..
c     .. Local Scalars ..

      INTEGER   IU, IUMAX, IUMIN, LENFMT, LU, NP, NPASS
c     ..
c     .. Intrinsic Functions ..

      INTRINSIC MIN
c     ..


      IF( NUMU.LT.1 )  RETURN

      WRITE( *, '(//,A)' )
     &   ' *******  AZIMUTHALLY AVERAGED INTENSITIES ' //
     &   '(at user polar angles)  ********'

      LENFMT = 8
      NPASS  = 1 + (NUMU-1) / LENFMT

      WRITE( *,'(/,A,/,A)') '   Optical   Polar Angle Cosines',
     &                      '     Depth'

      DO 20 NP = 1, NPASS

         IUMIN  = 1 + LENFMT * ( NP - 1 )
         IUMAX  = MIN( LENFMT*NP, NUMU )
         WRITE( *,'(/,10X,8F14.5)') ( UMU(IU), IU = IUMIN, IUMAX )

         DO 10 LU = 1, NTAU
            WRITE( *, '(0P,F10.4,1P,8E14.4)' ) UTAU( LU ),
     &           ( U0U( IU,LU ), IU = IUMIN, IUMAX )
   10    CONTINUE

   20 CONTINUE


      END

      SUBROUTINE PRTINP( NLYR, DTAUC, DTAUCP, SSALB, PMOM, TEMPER,
     &                   WVNMLO, WVNMHI, NTAU, UTAU, NSTR, NUMU, UMU,
     &                   NPHI, PHI, IBCND, FBEAM, UMU0, PHI0, FISOT,
     &                   LAMBER, ALBEDO, HL, BTEMP, TTEMP, TEMIS,
     &                   DELTAM, PLANK, ONLYFL, ACCUR, FLYR, LYRCUT,
     &                   OPRIM, TAUC, TAUCPR, MAXCMU, PRTMOM )

c        Print values of input variables

c   Called by- DISORT
c --------------------------------------------------------------------

c     .. Scalar Arguments ..

      LOGICAL   DELTAM, LAMBER, LYRCUT, ONLYFL, PLANK, PRTMOM
      INTEGER   IBCND, MAXCMU, NLYR, NPHI, NSTR, NTAU, NUMU
      REAL      ACCUR, ALBEDO, BTEMP, FBEAM, FISOT, PHI0, TEMIS, TTEMP,
     &          UMU0, WVNMHI, WVNMLO
c     ..
c     .. Array Arguments ..

      REAL      DTAUC( * ), DTAUCP( * ), FLYR( * ), HL( 0:MAXCMU ),
     &          OPRIM( * ), PHI( * ), PMOM( 0:MAXCMU, * ), SSALB( * ),
     &          TAUC( 0:* ), TAUCPR( 0:* ), TEMPER( 0:* ), UMU( * ),
     &          UTAU( * )
c     ..
c     .. Local Scalars ..

      INTEGER   IU, J, K, LC, LU
      REAL      YESSCT
c     ..


      WRITE( *, '(/,A,I4,A,I4)' ) ' No. streams =', NSTR,
     &       '     No. computational layers =', NLYR

      IF( IBCND.NE.1 ) WRITE( *, '(I4,A,10F10.4,/,(26X,10F10.4))' )
     &    NTAU,' User optical depths :', ( UTAU(LU), LU = 1, NTAU )

      IF( .NOT.ONLYFL ) WRITE( *, '(I4,A,10F9.5,/,(31X,10F9.5))' )
     &    NUMU,' User polar angle cosines :',( UMU(IU), IU = 1, NUMU )

      IF( .NOT.ONLYFL .AND. IBCND.NE.1 )
     &    WRITE( *, '(I4,A,10F9.2,/,(28X,10F9.2))' )
     &           NPHI,' User azimuthal angles :',( PHI(J), J = 1, NPHI )

      IF( .NOT.PLANK .OR. IBCND.EQ.1 )
     &    WRITE( *, '(A)' ) ' No thermal emission'


      WRITE( *, '(A,I2)' ) ' Boundary condition flag: IBCND =', IBCND

      IF( IBCND.EQ.0 ) THEN

         WRITE( *, '(A,1P,E11.3,A,0P,F8.5,A,F7.2,/,A,1P,E11.3)' )
     &          '    Incident beam with intensity =', FBEAM,
     &          ' and polar angle cosine = ', UMU0,
     &          '  and azimuth angle =', PHI0,
     &          '    plus isotropic incident intensity =', FISOT

         IF( LAMBER ) WRITE( *, '(A,0P,F8.4)' )
     &                '    Bottom albedo (Lambertian) =', ALBEDO

         IF( .NOT.LAMBER ) WRITE( *, '(A,/,(10X,10F9.5))' )
     &     '    Legendre coeffs of bottom bidirectional reflectivity :',
     &         ( HL( K ), K = 0, NSTR )

         IF( PLANK ) WRITE( *, '(A,2F14.4,/,A,F10.2,A,F10.2,A,F8.4)' )
     &       '    Thermal emission in wavenumber interval :', WVNMLO,
     &       WVNMHI,
     &       '    Bottom temperature =', BTEMP,
     &       '    Top temperature =', TTEMP,
     &       '    Top emissivity =',TEMIS

      ELSE IF( IBCND.EQ.1 ) THEN

         WRITE(*,'(A)') '    Isotropic illumination from top and bottom'
         WRITE( *, '(A,0P,F8.4)' )
     &          '    Bottom albedo (Lambertian) =', ALBEDO
      END IF


      IF( DELTAM ) WRITE( *, '(A)' ) ' Uses delta-M method'
      IF( .NOT.DELTAM ) WRITE( *, '(A)' ) ' Does not use delta-M method'


      IF( IBCND.EQ.1 ) THEN

         WRITE( *, '(A)' ) ' Calculate albedo and transmissivity of'//
     &                     ' medium vs. incident beam angle'

      ELSE IF( ONLYFL ) THEN

         WRITE( *, '(A)' )
     &          ' Calculate fluxes and azim-averaged intensities only'

      ELSE

         WRITE( *, '(A)' ) ' Calculate fluxes and intensities'

      END IF


      WRITE( *, '(A,1P,E11.2)' )
     &       ' Relative convergence criterion for azimuth series =',
     &       ACCUR

      IF( LYRCUT ) WRITE( *, '(A)' )
     &    ' Sets radiation = 0 below absorption optical depth 10'


c                                        ** Print layer variables
      IF( PLANK ) WRITE( *, FMT = 9180 )
      IF( .NOT.PLANK ) WRITE( *, FMT = 9190 )

      YESSCT = 0.0

      DO 10 LC = 1, NLYR

         YESSCT = YESSCT + SSALB( LC )

         IF( PLANK )
     &       WRITE(*,'(I4,2F10.4,F10.5,F12.5,2F10.4,F10.5,F9.4,F14.3)')
     &             LC, DTAUC( LC ), TAUC( LC ), SSALB( LC ), FLYR( LC ),
     &             DTAUCP( LC ), TAUCPR( LC ), OPRIM( LC ), PMOM(1,LC),
     &             TEMPER( LC-1 )

         IF( .NOT.PLANK )
     &       WRITE(*,'(I4,2F10.4,F10.5,F12.5,2F10.4,F10.5,F9.4)')
     &             LC, DTAUC( LC ), TAUC( LC ), SSALB( LC ), FLYR( LC ),
     &             DTAUCP( LC ), TAUCPR( LC ), OPRIM( LC ), PMOM( 1,LC )
   10 CONTINUE

      IF( PLANK ) WRITE( *, '(85X,F14.3)' ) TEMPER( NLYR )


      IF( PRTMOM .AND. YESSCT.GT.0.0 ) THEN

         WRITE( *, '(/,A)' ) ' Layer   Phase Function Moments'

         DO 20 LC = 1, NLYR

            IF( SSALB( LC ).GT.0.0 )
     &          WRITE( *, '(I6,10F11.6,/,(6X,10F11.6))' )
     &                 LC, ( PMOM( K, LC ), K = 0, NSTR )
   20    CONTINUE

      END IF

c                ** (Read every other line in these formats)

 9180 FORMAT( /, 37X, '<------------- Delta-M --------------->', /,
     &'                   Total    Single                           ',
     &               'Total    Single', /,
     &'       Optical   Optical   Scatter   Truncated   ',
     &   'Optical   Optical   Scatter    Asymm', /,
     &'         Depth     Depth    Albedo    Fraction     ',
     &     'Depth     Depth    Albedo   Factor   Temperature' )
 9190 FORMAT( /, 37X, '<------------- Delta-M --------------->', /,
     &'                   Total    Single                           ',
     &               'Total    Single', /,
     &'       Optical   Optical   Scatter   Truncated   ',
     &   'Optical   Optical   Scatter    Asymm', /,
     &'         Depth     Depth    Albedo    Fraction     ',
     &     'Depth     Depth    Albedo   Factor' )

      END

      SUBROUTINE PRTINT( UU, UTAU, NTAU, UMU, NUMU, PHI, NPHI, MAXULV,
     &                   MAXUMU )

c         Prints the intensity at user polar and azimuthal angles

c     All arguments are DISORT input or output variables

c   Called by- DISORT

c     LENFMT   Max number of azimuth angles PHI that can be printed
c                on one line, as set in FORMAT statement
c +-------------------------------------------------------------------+


c     .. Scalar Arguments ..

      INTEGER   MAXULV, MAXUMU, NPHI, NTAU, NUMU
c     ..
c     .. Array Arguments ..

      REAL      PHI( * ), UMU( * ), UTAU( * ), UU( MAXUMU, MAXULV, * )
c     ..
c     .. Local Scalars ..

      INTEGER   IU, J, JMAX, JMIN, LENFMT, LU, NP, NPASS
c     ..
c     .. Intrinsic Functions ..

      INTRINSIC MIN
c     ..


      IF( NPHI.LT.1 )  RETURN

      WRITE( *, '(//,A)' )
     &   ' *********  I N T E N S I T I E S  *********'

      LENFMT = 10
      NPASS  = 1 + (NPHI-1) / LENFMT

      WRITE( *, '(/,A,/,A,/,A)' )
     &   '             Polar   Azimuth angles (degrees)',
     &   '   Optical   Angle',
     &   '    Depth   Cosine'

      DO 30 LU = 1, NTAU

         DO 20 NP = 1, NPASS

            JMIN   = 1 + LENFMT * ( NP - 1 )
            JMAX   = MIN( LENFMT*NP, NPHI )

            WRITE( *, '(/,18X,10F11.2)' ) ( PHI(J), J = JMIN, JMAX )

            IF( NP.EQ.1 ) WRITE( *, '(F10.4,F8.4,1P,10E11.3)' )
     &             UTAU(LU), UMU(1), (UU(1, LU, J), J = JMIN, JMAX)
            IF( NP.GT.1 ) WRITE( *, '(10X,F8.4,1P,10E11.3)' )
     &                       UMU(1), (UU(1, LU, J), J = JMIN, JMAX)

            DO 10 IU = 2, NUMU
               WRITE( *, '(10X,F8.4,1P,10E11.3)' ) 
     &                 UMU( IU ), ( UU( IU, LU, J ), J = JMIN, JMAX )
   10       CONTINUE

   20    CONTINUE

   30 CONTINUE


      END

      SUBROUTINE QGAUSN( M, GMU, GWT )

c       Compute weights and abscissae for ordinary Gaussian quadrature
c       on the interval (0,1);  that is, such that

c           sum(i=1 to M) ( GWT(i) f(GMU(i)) )

c       is a good approximation to

c           integral(0 to 1) ( f(x) dx )

c   INPUT :    M       order of quadrature rule

c   OUTPUT :  GMU(I)   array of abscissae (I = 1 TO M)
c             GWT(I)   array of weights (I = 1 TO M)

c   REFERENCE:  Davis, P.J. and P. Rabinowitz, Methods of Numerical
c                   Integration, Academic Press, New York, pp. 87, 1975

c   METHOD:  Compute the abscissae as roots of the Legendre
c            polynomial P-sub-M using a cubically convergent
c            refinement of Newton's method.  Compute the
c            weights from EQ. 2.7.3.8 of Davis/Rabinowitz.  Note
c            that Newton's method can very easily diverge; only a
c            very good initial guess can guarantee convergence.
c            The initial guess used here has never led to divergence
c            even for M up to 1000.

c   ACCURACY:  relative error no better than TOL or computer
c              precision (machine epsilon), whichever is larger

c   INTERNAL VARIABLES:

c    ITER      : number of Newton Method iterations
c    MAXIT     : maximum allowed iterations of Newton Method
c    PM2,PM1,P : 3 successive Legendre polynomials
c    PPR       : derivative of Legendre polynomial
c    P2PRI     : 2nd derivative of Legendre polynomial
c    TOL       : convergence criterion for Legendre poly root iteration
c    X,XI      : successive iterates in cubically-convergent version
c                of Newtons Method (seeking roots of Legendre poly.)

c   Called by- SETDIS, SURFAC
c   Calls- D1MACH, ERRMSG
c +-------------------------------------------------------------------+

c     .. Scalar Arguments ..

      INTEGER   M
c     ..
c     .. Array Arguments ..

      REAL      GMU( M ), GWT( M )
c     ..
c     .. Local Scalars ..

      INTEGER   ITER, K, LIM, MAXIT, NN, NP1
      REAL      CONA, PI, T
      DOUBLE PRECISION EN, NNP1, ONE, P, P2PRI, PM1, PM2, PPR, PROD,
     &                 TMP, TOL, TWO, X, XI
c     ..
c     .. External Functions ..

      DOUBLE PRECISION D1MACH
      EXTERNAL  D1MACH
c     ..
c     .. External Subroutines ..

      EXTERNAL  ERRMSG
c     ..
c     .. Intrinsic Functions ..

      INTRINSIC ABS, ASIN, COS, FLOAT, MOD, TAN
c     ..
      SAVE      PI, TOL

      DATA      PI / 0.0 / , MAXIT / 1000 / , ONE / 1.D0 / ,
     &          TWO / 2.D0 /


      IF( PI.EQ.0.0 ) THEN

         PI   = 2.*ASIN( 1.0 )
         TOL  = 10.*D1MACH( 4 )

      END IF


      IF( M.LT.1 ) CALL ERRMSG( 'QGAUSN--Bad value of M',.True.)

      IF( M.EQ.1 ) THEN

         GMU( 1 ) = 0.5
         GWT( 1 ) = 1.0
         RETURN

      END IF

      EN   = M
      NP1  = M + 1
      NNP1 = M*NP1
      CONA = FLOAT( M - 1 ) / ( 8*M**3 )

      LIM  = M / 2

      DO 30 K = 1, LIM
c                                        ** Initial guess for k-th root
c                                           of Legendre polynomial, from
c                                           Davis/Rabinowitz (2.7.3.3a)
         T  = ( 4*K - 1 )*PI / ( 4*M + 2 )
         X  = COS( T + CONA / TAN( T ) )
         ITER = 0
c                                        ** Upward recurrence for
c                                           Legendre polynomials
   10    CONTINUE
         ITER   = ITER + 1
         PM2    = ONE
         PM1    = X

         DO 20 NN = 2, M
            P    = ( ( 2*NN - 1 )*X*PM1 - ( NN - 1 )*PM2 ) / NN
            PM2  = PM1
            PM1  = P
   20    CONTINUE
c                                              ** Newton Method
         TMP    = ONE / ( ONE - X**2 )
         PPR    = EN*( PM2 - X*P )*TMP
         P2PRI  = ( TWO*X*PPR - NNP1*P )*TMP
         XI     = X - ( P / PPR )*( ONE +
     &            ( P / PPR )*P2PRI / ( TWO*PPR ) )

c                                              ** Check for convergence
         IF( ABS( XI - X ).GT.TOL ) THEN

            IF( ITER.GT.MAXIT )
     &          CALL ERRMSG( 'QGAUSN--max iteration count',.True.)

            X  = XI
            GO TO  10

         END IF
c                             ** Iteration finished--calculate weights,
c                                abscissae for (-1,1)
         GMU( K ) = -X
         GWT( K ) = TWO / ( TMP*( EN*PM2 )**2 )
         GMU( NP1 - K ) = -GMU( K )
         GWT( NP1 - K ) = GWT( K )
   30 CONTINUE
c                                    ** Set middle abscissa and weight
c                                       for rules of odd order
      IF( MOD( M,2 ).NE.0 ) THEN

         GMU( LIM + 1 ) = 0.0
         PROD   = ONE

         DO 40 K = 3, M, 2
            PROD   = PROD * K / ( K - 1 )
   40    CONTINUE

         GWT( LIM + 1 ) = TWO / PROD**2
      END IF

c                                        ** Convert from (-1,1) to (0,1)
      DO 50 K = 1, M
         GMU( K ) = 0.5*GMU( K ) + 0.5
         GWT( K ) = 0.5*GWT( K )
   50 CONTINUE


      END

      SUBROUTINE SETDIS( dsdh, nid, tausla, tauslau, mu2,
     &                   CMU, CWT, DELTAM, DTAUC, DTAUCP, EXPBEA, FBEAM,
     &                   FLYR, GL, HL, HLPR, IBCND, LAMBER, LAYRU,
     &                   LYRCUT, MAXUMU, MAXCMU, MXCMU, NCUT, NLYR,
     &                   NTAU, NN, NSTR, PLANK, NUMU, ONLYFL, OPRIM,
     &                   PMOM, SSALB, TAUC, TAUCPR, UTAU, UTAUPR, UMU,
     &                   UMU0, USRTAU, USRANG )

c          Perform miscellaneous setting-up operations

c       INPUT :  all are DISORT input variables (see DOC file)

c       OUTPUT:  NTAU,UTAU   if USRTAU = FALSE
c                NUMU,UMU    if USRANG = FALSE
c                CMU,CWT     computational polar angles and
c                               corresponding quadrature weights
c                EXPBEA      transmission of direct beam
c                FLYR        truncated fraction in delta-M method
c                GL          phase function Legendre coefficients multi-
c                              plied by (2L+1) and single-scatter albedo
c                HLPR        Legendre moments of surface bidirectional
c                              reflectivity, times 2K+1
c                LAYRU       Computational layer in which UTAU falls
c                LYRCUT      flag as to whether radiation will be zeroed
c                              below layer NCUT
c                NCUT        computational layer where absorption
c                              optical depth first exceeds  ABSCUT
c                NN          NSTR / 2
c                OPRIM       delta-M-scaled single-scatter albedo
c                TAUCPR      delta-M-scaled optical depth
c                UTAUPR      delta-M-scaled version of  UTAU

c   Called by- DISORT
c   Calls- QGAUSN, ERRMSG
c ----------------------------------------------------------------------

      INCLUDE 'tuv.params'

c     .. Scalar Arguments ..

      LOGICAL   DELTAM, LAMBER, LYRCUT, ONLYFL, PLANK, USRANG, USRTAU
      INTEGER   IBCND, MAXCMU, MAXUMU, MXCMU, NCUT, NLYR, NN, NSTR,
     &          NTAU, NUMU
      REAL      FBEAM, UMU0

c geometry
      REAL dsdh(0:kz,kz)
      INTEGER nid(0:kz)
      REAL tausla(0:kz), tauslau(0:kz), mu2(0:kz)
      REAL sum, sumu
c     ..
c     .. Array Arguments ..

      INTEGER   LAYRU( * )
      REAL      CMU( MXCMU ), CWT( MXCMU ), DTAUC( * ), DTAUCP( * ),
     &          EXPBEA( 0:* ), FLYR( * ), GL( 0:MXCMU, * ),
     &          HL( 0:MAXCMU ), HLPR( 0:MXCMU ), OPRIM( * ),
     &          PMOM( 0:MAXCMU, * ), SSALB( * ), TAUC( 0:* ),
     &          TAUCPR( 0:* ), UMU( MAXUMU ), UTAU( * ), UTAUPR( * )
c     ..
c     .. Local Scalars ..

      INTEGER   IQ, IU, K, LC, LU
      REAL      ABSCUT, ABSTAU, F

      REAL      R1MACH
      EXTERNAL  R1MACH
c     ..
c     .. External Subroutines ..

      EXTERNAL  ERRMSG, QGAUSN
c     ..
c     .. Intrinsic Functions ..

      INTRINSIC ABS, EXP
c     ..
      DATA      ABSCUT / 30. /


      IF( .NOT.USRTAU ) THEN
c                              ** Set output levels at computational
c                                 layer boundaries
         NTAU  = NLYR + 1

         DO 10 LC = 0, NTAU - 1
            UTAU( LC + 1 ) = TAUC( LC )
   10    CONTINUE

      END IF
c                        ** Apply delta-M scaling and move description
c                           of computational layers to local variables
      EXPBEA( 0 ) = 1.0
      TAUCPR( 0 ) = 0.0
      ABSTAU      = 0.0
      do i = 0, kz
       tausla( i ) = 0.0
       tauslau( i ) = 0.0
       mu2(i) = 1./largest
      end do

      DO 40 LC = 1, NLYR

         PMOM( 0, LC ) = 1.0

         IF( ABSTAU.LT.ABSCUT ) NCUT  = LC

         ABSTAU = ABSTAU + ( 1.- SSALB( LC ) )*DTAUC( LC )

         IF( .NOT.DELTAM ) THEN

            OPRIM( LC )  = SSALB( LC )
            DTAUCP( LC ) = DTAUC( LC )
            TAUCPR( LC ) = TAUC( LC )

            DO 20 K = 0, NSTR - 1
               GL( K, LC ) = ( 2*K + 1 )*OPRIM( LC )*PMOM( K, LC )
   20       CONTINUE

            F  = 0.0


         ELSE
c                                    ** Do delta-M transformation

            F  = PMOM( NSTR, LC )
            OPRIM(LC) = SSALB(LC) * ( 1.- F ) / ( 1.- F * SSALB(LC) )
            DTAUCP( LC ) = ( 1.- F*SSALB( LC ) )*DTAUC( LC )
            TAUCPR( LC ) = TAUCPR( LC-1 ) + DTAUCP( LC )

            DO 30 K = 0, NSTR - 1
               GL( K, LC ) = ( 2*K + 1 ) * OPRIM( LC ) *
     &                       ( PMOM( K,LC ) - F ) / ( 1.- F )
   30       CONTINUE

         END IF

         FLYR( LC )   = F
         EXPBEA( LC ) = 0.0

   40 CONTINUE
c 
* calculate slant optical depth
*              
         IF(umu0 .LT. 0.0) THEN
           IF(nid(0) .LT. 0) THEN
             tausla(0) = largest
             tauslau(0) = largest
           ELSE
             sum = 0.0
             sumu = 0.0
             DO lc = 1, nid(0)
               sum = sum + 2.*dtaucp(lc)*dsdh(0,lc)
               sumu = sumu + 2.*dtauc(lc)*dsdh(0,lc)
             END DO
             tausla(0) = sum 
             tauslau(0) = sumu 
           END IF
         END IF

         expbea( 0 ) = EXP( -tausla( 0 ) )

*
         DO 41, lc = 1, nlyr
          IF(nid(lc) .LT. 0) THEN
            tausla(lc) = largest
            tauslau(lc) = largest
          ELSE
            sum = 0.0
            sumu = 0.0
            DO lu = 1, MIN(nid(lc),lc)
               sum = sum + dtaucp(lu)*dsdh(lc,lu)
               sumu = sumu + dtauc(lu)*dsdh(lc,lu)
            ENDDO
            DO lu = MIN(nid(lc),lc)+1,nid(lc)
               sum = sum + 2.*dtaucp(lu)*dsdh(lc,lu)
               sumu = sumu + 2.*dtauc(lu)*dsdh(lc,lu)
            ENDDO
            tausla(lc) = sum 
            tauslau(lc) = sumu 
            IF(tausla(lc) .EQ. tausla(lc-1)) THEN
              mu2(lc) = largest
            ELSE
              mu2(lc) = (taucpr(lc)-taucpr(lc-1))
     $                         /(tausla(lc)-tausla(lc-1))
              mu2(lc) = SIGN( AMAX1(ABS(mu2(lc)),1./largest),
     $                     mu2(lc) )
            END IF
          END IF
          expbea(lc) = EXP( -tausla( lc ) )
 41      CONTINUE

c                      ** If no thermal emission, cut off medium below
c                         absorption optical depth = ABSCUT ( note that
c                         delta-M transformation leaves absorption
c                         optical depth invariant ).  Not worth the
c                         trouble for one-layer problems, though.
      LYRCUT = .FALSE.

      IF( ABSTAU.GE.ABSCUT .AND. .NOT.PLANK .AND. IBCND.NE.1 .AND.
     &    NLYR.GT.1 ) LYRCUT = .TRUE.

      IF( .NOT.LYRCUT ) NCUT   = NLYR

c                             ** Set arrays defining location of user
c                             ** output levels within delta-M-scaled
c                             ** computational mesh
      DO 70 LU = 1, NTAU

         DO 50 LC = 1, NLYR

            IF( UTAU( LU ).GE.TAUC( LC - 1 ) .AND.
     &          UTAU( LU ).LE.TAUC( LC ) ) GO TO  60

   50    CONTINUE
         LC   = NLYR

   60    CONTINUE
         UTAUPR( LU ) = UTAU( LU )
         IF( DELTAM ) UTAUPR( LU ) = TAUCPR( LC - 1 ) +
     &                               ( 1.- SSALB( LC )*FLYR( LC ) )*
     &                               ( UTAU( LU ) - TAUC( LC-1 ) )
         LAYRU( LU ) = LC

   70 CONTINUE
c                      ** Calculate computational polar angle cosines
c                         and associated quadrature weights for Gaussian
c                         quadrature on the interval (0,1) (upward)
      NN   = NSTR / 2

      CALL QGAUSN( NN, CMU, CWT )
c                                  ** Downward (neg) angles and weights
      DO 80 IQ = 1, NN
         CMU( IQ + NN ) = - CMU( IQ )
         CWT( IQ + NN ) = CWT( IQ )
   80 CONTINUE


c     IF( FBEAM.GT.0.0 ) THEN
c                               ** Compare beam angle to comput. angles
         DO 90 IQ = 1, NN

C                      ** Dither mu2 if it is close to one of the 
C                         quadrature angles.

         DO  lc = 1, nlyr
          IF (  ABS(mu2(lc)) .lt. 1.E5 ) THEN
            IF( ABS( 1. - ABS(mu2(lc))/CMU( IQ ) ) .LT. 0.05 ) 
     &           mu2(lc) = mu2(lc)*0.999
          ENDIF
         END DO

   90    CONTINUE

c     END IF

      IF( .NOT.USRANG .OR. ( ONLYFL .AND. MAXUMU.GE.NSTR ) ) THEN

c                                   ** Set output polar angles to
c                                      computational polar angles
         NUMU   = NSTR

         DO 100 IU = 1, NN
            UMU( IU ) = - CMU( NN + 1 - IU )
  100    CONTINUE

         DO 110 IU = NN + 1, NSTR
            UMU( IU ) = CMU( IU - NN )
  110    CONTINUE

      END IF


      IF( USRANG .AND. IBCND.EQ.1 ) THEN

c                               ** Shift positive user angle cosines to
c                                  upper locations and put negatives
c                                  in lower locations
         DO 120 IU = 1, NUMU
            UMU( IU + NUMU ) = UMU( IU )
  120    CONTINUE

         DO 130 IU = 1, NUMU
            UMU( IU ) = -UMU( 2*NUMU + 1 - IU )
  130    CONTINUE

         NUMU   = 2*NUMU

      END IF


      IF( .NOT.LYRCUT .AND. .NOT.LAMBER ) THEN

         DO 140 K = 0, NSTR
            HLPR( K ) = ( 2*K + 1 )*HL( K )
  140    CONTINUE

      END IF


      END

      SUBROUTINE SETMTX( BDR, CBAND, CMU, CWT, DELM0, DTAUCP, GC, KK,
     &                   LAMBER, LYRCUT, MI, MI9M2, MXCMU, NCOL, NCUT,
     &                   NNLYRI, NN, NSTR, TAUCPR, WK )

c        Calculate coefficient matrix for the set of equations
c        obtained from the boundary conditions and the continuity-
c        of-intensity-at-layer-interface equations;  store in the
c        special banded-matrix format required by LINPACK routines

c     I N P U T      V A R I A B L E S:

c       BDR      :  Surface bidirectional reflectivity
c       CMU      :  Abscissae for Gauss quadrature over angle cosine
c       CWT      :  Weights for Gauss quadrature over angle cosine
c       DELM0    :  Kronecker delta, delta-sub-m0
c       GC       :  Eigenvectors at polar quadrature angles, SC(1)
c       KK       :  Eigenvalues of coeff. matrix in Eq. SS(7)
c       LYRCUT   :  Logical flag for truncation of comput. layer
c       NN       :  Number of streams in a hemisphere (NSTR/2)
c       NCUT     :  Total number of computational layers considered
c       TAUCPR   :  Cumulative optical depth (delta-M-scaled)
c       (remainder are DISORT input variables)

c   O U T P U T     V A R I A B L E S:

c       CBAND    :  Left-hand side matrix of linear system Eq. SC(5),
c                      scaled by Eq. SC(12); in banded form required
c                      by LINPACK solution routines
c       NCOL     :  Counts of columns in CBAND

c   I N T E R N A L    V A R I A B L E S:

c       IROW     :  Points to row in CBAND
c       JCOL     :  Points to position in layer block
c       LDA      :  Row dimension of CBAND
c       NCD      :  Number of diagonals below or above main diagonal
c       NSHIFT   :  For positioning number of rows in band storage
c       WK       :  Temporary storage for EXP evaluations

c   Called by- DISORT, ALBTRN
c   Calls- ZEROIT
c +--------------------------------------------------------------------+


c     .. Scalar Arguments ..

      LOGICAL   LAMBER, LYRCUT
      INTEGER   MI, MI9M2, MXCMU, NCOL, NCUT, NN, NNLYRI, NSTR
      REAL      DELM0
c     ..
c     .. Array Arguments ..

      REAL      BDR( MI, 0:MI ), CBAND( MI9M2, NNLYRI ), CMU( MXCMU ),
     &          CWT( MXCMU ), DTAUCP( * ), GC( MXCMU, MXCMU, * ),
     &          KK( MXCMU, * ), TAUCPR( 0:* ), WK( MXCMU )
c     ..
c     .. Local Scalars ..

      INTEGER   IQ, IROW, JCOL, JQ, K, LC, LDA, NCD, NNCOL, NSHIFT
      REAL      EXPA, SUM
c     ..
c     .. External Subroutines ..

      EXTERNAL  ZEROIT
c     ..
c     .. Intrinsic Functions ..

      INTRINSIC EXP
c     ..


      CALL ZEROIT( CBAND, MI9M2*NNLYRI )

      NCD    = 3*NN - 1
      LDA    = 3*NCD + 1
      NSHIFT = LDA - 2*NSTR + 1
      NCOL   = 0
c                         ** Use continuity conditions of Eq. STWJ(17)
c                            to form coefficient matrix in STWJ(20);
c                            employ scaling transformation STWJ(22)
      DO 60 LC = 1, NCUT

         DO 10 IQ = 1, NN
            WK( IQ ) = EXP( KK( IQ,LC )*DTAUCP( LC ) )
   10    CONTINUE

         JCOL  = 0

         DO 30 IQ = 1, NN

            NCOL  = NCOL + 1
            IROW  = NSHIFT - JCOL

            DO 20 JQ = 1, NSTR
               CBAND( IROW + NSTR, NCOL ) =   GC( JQ, IQ, LC )
               CBAND( IROW, NCOL )        = - GC( JQ, IQ, LC )*WK( IQ )
               IROW  = IROW + 1
   20       CONTINUE

            JCOL  = JCOL + 1

   30    CONTINUE


         DO 50 IQ = NN + 1, NSTR

            NCOL  = NCOL + 1
            IROW  = NSHIFT - JCOL

            DO 40 JQ = 1, NSTR
               CBAND( IROW + NSTR, NCOL ) =   GC( JQ, IQ, LC )*
     &                                          WK( NSTR + 1 - IQ )
               CBAND( IROW, NCOL )        = - GC( JQ, IQ, LC )
               IROW  = IROW + 1
   40       CONTINUE

            JCOL  = JCOL + 1

   50    CONTINUE

   60 CONTINUE
c                  ** Use top boundary condition of STWJ(20a) for
c                     first layer

      JCOL  = 0

      DO 80 IQ = 1, NN

         EXPA  = EXP( KK( IQ,1 )*TAUCPR( 1 ) )
         IROW  = NSHIFT - JCOL + NN

         DO 70 JQ = NN, 1, -1
            CBAND( IROW, JCOL + 1 ) = GC( JQ, IQ, 1 )*EXPA
            IROW  = IROW + 1
   70    CONTINUE

         JCOL  = JCOL + 1

   80 CONTINUE


      DO 100 IQ = NN + 1, NSTR

         IROW  = NSHIFT - JCOL + NN

         DO 90 JQ = NN, 1, -1
            CBAND( IROW, JCOL + 1 ) = GC( JQ, IQ, 1 )
            IROW  = IROW + 1
   90    CONTINUE

         JCOL  = JCOL + 1

  100 CONTINUE
c                           ** Use bottom boundary condition of
c                              STWJ(20c) for last layer

      NNCOL = NCOL - NSTR
      JCOL  = 0

      DO 130 IQ = 1, NN

         NNCOL  = NNCOL + 1
         IROW   = NSHIFT - JCOL + NSTR

         DO 120 JQ = NN + 1, NSTR

            IF( LYRCUT .OR. ( LAMBER .AND. DELM0.EQ.0 ) ) THEN

c                          ** No azimuthal-dependent intensity if Lam-
c                             bert surface; no intensity component if
c                             truncated bottom layer

               CBAND( IROW, NNCOL ) = GC( JQ, IQ, NCUT )

            ELSE

               SUM  = 0.0

               DO 110 K = 1, NN
                  SUM  = SUM + CWT( K )*CMU( K )*BDR( JQ - NN, K )*
     &                     GC( NN + 1 - K, IQ, NCUT )
  110          CONTINUE

               CBAND( IROW, NNCOL ) = GC( JQ, IQ, NCUT ) -
     &                                ( 1.+ DELM0 )*SUM
            END IF

            IROW  = IROW + 1

  120    CONTINUE

         JCOL  = JCOL + 1

  130 CONTINUE


      DO 160 IQ = NN + 1, NSTR

         NNCOL  = NNCOL + 1
         IROW   = NSHIFT - JCOL + NSTR
         EXPA   = WK( NSTR + 1 - IQ )

         DO 150 JQ = NN + 1, NSTR

            IF( LYRCUT .OR. ( LAMBER .AND. DELM0.EQ.0 ) ) THEN

               CBAND( IROW, NNCOL ) = GC( JQ, IQ, NCUT )*EXPA

            ELSE

               SUM  = 0.0

               DO 140 K = 1, NN
                  SUM  = SUM + CWT( K )*CMU( K )*BDR( JQ - NN, K )*
     &                         GC( NN + 1 - K, IQ, NCUT )
  140          CONTINUE

               CBAND( IROW, NNCOL ) = ( GC( JQ,IQ,NCUT ) -
     &                                ( 1.+ DELM0 )*SUM )*EXPA
            END IF

            IROW  = IROW + 1

  150    CONTINUE

         JCOL  = JCOL + 1

  160 CONTINUE

      END


      SUBROUTINE SOLEIG( AMB, APB, ARRAY, CMU, CWT, GL, MI, MAZIM,
     &                   MXCMU, NN, NSTR, YLMC, CC, EVECC, EVAL, KK, GC,
     &                   AAD, EVECCD, EVALD, WKD )

c         Solves eigenvalue/vector problem necessary to construct
c         homogeneous part of discrete ordinate solution; STWJ(8b)
c         ** NOTE ** Eigenvalue problem is degenerate when single
c                    scattering albedo = 1;  present way of doing it
c                    seems numerically more stable than alternative
c                    methods that we tried

c   I N P U T     V A R I A B L E S:

c       GL     :  Delta-M scaled Legendre coefficients of phase function
c                    (including factors 2l+1 and single-scatter albedo)
c       CMU    :  Computational polar angle cosines
c       CWT    :  Weights for quadrature over polar angle cosine
c       MAZIM  :  Order of azimuthal component
c       NN     :  Half the total number of streams
c       YLMC   :  Normalized associated Legendre polynomial
c                    at the quadrature angles CMU
c       (remainder are DISORT input variables)

c   O U T P U T    V A R I A B L E S:

c       CC     :  C-sub-ij in Eq. SS(5); needed in SS(15&18)
c       EVAL   :  NN eigenvalues of Eq. SS(12) on return from ASYMTX
c                    but then square roots taken
c       EVECC  :  NN eigenvectors  (G+) - (G-)  on return
c                    from ASYMTX ( column j corresponds to EVAL(j) )
c                    but then  (G+) + (G-)  is calculated from SS(10),
c                    G+  and  G-  are separated, and  G+  is stacked on
c                    top of  G-  to form NSTR eigenvectors of SS(7)
c       GC     :  Permanent storage for all NSTR eigenvectors, but
c                    in an order corresponding to KK
c       KK     :  Permanent storage for all NSTR eigenvalues of SS(7),
c                    but re-ordered with negative values first ( square
c                    roots of EVAL taken and negatives added )

c   I N T E R N A L   V A R I A B L E S:

c       AMB,APB :  Matrices (alpha-beta), (alpha+beta) in reduced
c                    eigenvalue problem
c       ARRAY   :  Complete coefficient matrix of reduced eigenvalue
c                    problem: (alfa+beta)*(alfa-beta)
c       GPPLGM  :  (G+) + (G-) (cf. Eqs. SS(10-11))
c       GPMIGM  :  (G+) - (G-) (cf. Eqs. SS(10-11))
c       WKD     :  Scratch array required by ASYMTX

c   Called by- DISORT, ALBTRN
c   Calls- ASYMTX, ERRMSG
c +-------------------------------------------------------------------+


c     .. Scalar Arguments ..

      INTEGER   MAZIM, MI, MXCMU, NN, NSTR
c     ..
c     .. Array Arguments ..

      REAL      AMB( MI, MI ), APB( MI, MI ), ARRAY( MI, * ),
     &          CC( MXCMU, MXCMU ), CMU( MXCMU ), CWT( MXCMU ),
     &          EVAL( MI ), EVECC( MXCMU, MXCMU ), GC( MXCMU, MXCMU ),
     &          GL( 0:MXCMU ), KK( MXCMU ), YLMC( 0:MXCMU, MXCMU )
      DOUBLE PRECISION AAD( MI, MI ), EVALD( MI ), EVECCD( MI, MI ),
     &                 WKD( MXCMU )
c     ..
c     .. Local Scalars ..

      INTEGER   IER, IQ, JQ, KQ, L
      REAL      ALPHA, BETA, GPMIGM, GPPLGM, SUM
c     ..
c     .. External Subroutines ..

      EXTERNAL  ASYMTX, ERRMSG
c     ..
c     .. Intrinsic Functions ..

      INTRINSIC ABS, SQRT
c     ..

c                             ** Calculate quantities in Eqs. SS(5-6)
      DO 40 IQ = 1, NN

         DO 20 JQ = 1, NSTR

            SUM  = 0.0
            DO 10 L = MAZIM, NSTR - 1
               SUM  = SUM + GL( L )*YLMC( L, IQ )*YLMC( L, JQ )
   10       CONTINUE

            CC( IQ, JQ ) = 0.5*SUM*CWT( JQ )

   20    CONTINUE

         DO 30 JQ = 1, NN
c                             ** Fill remainder of array using symmetry
c                                relations  C(-mui,muj) = C(mui,-muj)
c                                and        C(-mui,-muj) = C(mui,muj)

            CC( IQ + NN, JQ ) = CC( IQ, JQ + NN )
            CC( IQ + NN, JQ + NN ) = CC( IQ, JQ )

c                                       ** Get factors of coeff. matrix
c                                          of reduced eigenvalue problem

            ALPHA  = CC( IQ, JQ ) / CMU( IQ )
            BETA   = CC( IQ, JQ + NN ) / CMU( IQ )
            AMB( IQ, JQ ) = ALPHA - BETA
            APB( IQ, JQ ) = ALPHA + BETA

   30    CONTINUE

         AMB( IQ, IQ ) = AMB( IQ, IQ ) - 1.0 / CMU( IQ )
         APB( IQ, IQ ) = APB( IQ, IQ ) - 1.0 / CMU( IQ )

   40 CONTINUE
c                      ** Finish calculation of coefficient matrix of
c                         reduced eigenvalue problem:  get matrix
c                         product (alfa+beta)*(alfa-beta); SS(12)
      DO 70 IQ = 1, NN

         DO 60 JQ = 1, NN

            SUM  = 0.
            DO 50 KQ = 1, NN
               SUM  = SUM + APB( IQ, KQ )*AMB( KQ, JQ )
   50       CONTINUE

            ARRAY( IQ, JQ ) = SUM

   60    CONTINUE

   70 CONTINUE
c                      ** Find (real) eigenvalues and eigenvectors

      CALL ASYMTX( ARRAY, EVECC, EVAL, NN, MI, MXCMU, IER, WKD, AAD,
     &             EVECCD, EVALD )

      IF( IER.GT.0 ) THEN

         WRITE( *, FMT = '(//,A,I4,A)' ) ' ASYMTX--eigenvalue no. ',
     &      IER, '  didnt converge.  Lower-numbered eigenvalues wrong.'

         CALL ERRMSG( 'ASYMTX--convergence problems',.True.)

      END IF

CDIR$ IVDEP
      DO 80 IQ = 1, NN
         EVAL( IQ )    = SQRT( ABS( EVAL( IQ ) ) )
         KK( IQ + NN ) = EVAL( IQ )
c                                      ** Add negative eigenvalue
         KK( NN + 1 - IQ ) = -EVAL( IQ )
   80 CONTINUE

c                          ** Find eigenvectors (G+) + (G-) from SS(10)
c                             and store temporarily in APB array
      DO 110 JQ = 1, NN

         DO 100 IQ = 1, NN

            SUM  = 0.
            DO 90 KQ = 1, NN
               SUM  = SUM + AMB( IQ, KQ )*EVECC( KQ, JQ )
   90       CONTINUE

            APB( IQ, JQ ) = SUM / EVAL( JQ )

  100    CONTINUE

  110 CONTINUE


      DO 130 JQ = 1, NN
CDIR$ IVDEP
         DO 120 IQ = 1, NN

            GPPLGM = APB( IQ, JQ )
            GPMIGM = EVECC( IQ, JQ )
c                                ** Recover eigenvectors G+,G- from
c                                   their sum and difference; stack them
c                                   to get eigenvectors of full system
c                                   SS(7) (JQ = eigenvector number)

            EVECC( IQ,      JQ ) = 0.5*( GPPLGM + GPMIGM )
            EVECC( IQ + NN, JQ ) = 0.5*( GPPLGM - GPMIGM )

c                                ** Eigenvectors corresponding to
c                                   negative eigenvalues (corresp. to
c                                   reversing sign of 'k' in SS(10) )
            GPPLGM = - GPPLGM
            EVECC(IQ,   JQ+NN) = 0.5 * ( GPPLGM + GPMIGM )
            EVECC(IQ+NN,JQ+NN) = 0.5 * ( GPPLGM - GPMIGM )
            GC( IQ+NN,   JQ+NN )   = EVECC( IQ,    JQ )
            GC( NN+1-IQ, JQ+NN )   = EVECC( IQ+NN, JQ )
            GC( IQ+NN,   NN+1-JQ ) = EVECC( IQ,    JQ+NN )
            GC( NN+1-IQ, NN+1-JQ ) = EVECC( IQ+NN, JQ+NN )

  120    CONTINUE

  130 CONTINUE


      END

      SUBROUTINE SOLVE0( B, BDR, BEM, BPLANK, CBAND, CMU, CWT, EXPBEA,
     &                   FBEAM, FISOT, IPVT, LAMBER, LL, LYRCUT, MAZIM,
     &                   MI, MI9M2, MXCMU, NCOL, NCUT, NN, NSTR, NNLYRI,
     &                   PI, TPLANK, TAUCPR, UMU0, Z, ZZ, ZPLK0, ZPLK1 )

c        Construct right-hand side vector B for general boundary
c        conditions STWJ(17) and solve system of equations obtained
c        from the boundary conditions and the continuity-of-
c        intensity-at-layer-interface equations.
c        Thermal emission contributes only in azimuthal independence.

c     I N P U T      V A R I A B L E S:

c       BDR      :  Surface bidirectional reflectivity
c       BEM      :  Surface bidirectional emissivity
c       BPLANK   :  Bottom boundary thermal emission
c       CBAND    :  Left-hand side matrix of linear system Eq. SC(5),
c                   scaled by Eq. SC(12); in banded form required
c                   by LINPACK solution routines
c       CMU      :  Abscissae for Gauss quadrature over angle cosine
c       CWT      :  Weights for Gauss quadrature over angle cosine
c       EXPBEA   :  Transmission of incident beam, EXP(-TAUCPR/UMU0)
c       LYRCUT   :  Logical flag for truncation of comput. layer
c       MAZIM    :  Order of azimuthal component
c       ncol     :  Counts of columns in CBAND
c       NN       :  Order of double-Gauss quadrature (NSTR/2)
c       NCUT     :  Total number of computational layers considered
c       TPLANK   :  Top boundary thermal emission
c       TAUCPR   :  Cumulative optical depth (delta-M-scaled)
c       ZZ       :  Beam source vectors in Eq. SS(19)
c       ZPLK0    :  Thermal source vectors Z0, by solving Eq. SS(16)
c       ZPLK1    :  Thermal source vectors Z1, by solving Eq. SS(16)
c       (remainder are DISORT input variables)

c   O U T P U T     V A R I A B L E S:

c       B        :  Right-hand side vector of Eq. SC(5) going into
c                   SGBSL; returns as solution vector of Eq. SC(12),
c                   constants of integration without exponential term
c
c      LL        :  Permanent storage for B, but re-ordered

c   I N T E R N A L    V A R I A B L E S:

c       IPVT     :  Integer vector of pivot indices
c       IT       :  Pointer for position in  B
c       NCD      :  Number of diagonals below or above main diagonal
c       RCOND    :  Indicator of singularity for CBAND
c       Z        :  Scratch array required by SGBCO

c   Called by- DISORT
c   Calls- ZEROIT, SGBCO, ERRMSG, SGBSL
c +-------------------------------------------------------------------+


c     .. Scalar Arguments ..

      LOGICAL   LAMBER, LYRCUT
      INTEGER   MAZIM, MI, MI9M2, MXCMU, NCOL, NCUT, NN, NNLYRI, NSTR
      REAL      BPLANK, FBEAM, FISOT, PI, TPLANK, UMU0
c     ..
c     .. Array Arguments ..

      INTEGER   IPVT( * )
      REAL      B( NNLYRI ), BDR( MI, 0:MI ), BEM( MI ),
     &          CBAND( MI9M2, NNLYRI ), CMU( MXCMU ), CWT( MXCMU ),
     &          EXPBEA( 0:* ), LL( MXCMU, * ), TAUCPR( 0:* ),
     &          Z( NNLYRI ), ZPLK0( MXCMU, * ), ZPLK1( MXCMU, * ),
     &          ZZ( MXCMU, * )
c     ..
c     .. Local Scalars ..

      INTEGER   IPNT, IQ, IT, JQ, LC, NCD
      REAL      RCOND, SUM
c     ..
c     .. External Subroutines ..

      EXTERNAL  ERRMSG, SGBCO, SGBSL, ZEROIT
c     ..


      CALL ZEROIT( B, NNLYRI )
c                              ** Construct B,  STWJ(20a,c) for
c                                 parallel beam + bottom reflection +
c                                 thermal emission at top and/or bottom

      IF( MAZIM.GT.0 .AND. FBEAM.GT.0.0 ) THEN

c                                         ** Azimuth-dependent case
c                                            (never called if FBEAM = 0)
         IF( LYRCUT .OR. LAMBER ) THEN

c               ** No azimuthal-dependent intensity for Lambert surface;
c                  no intensity component for truncated bottom layer

            DO 10 IQ = 1, NN
c                                                  ** Top boundary
               B( IQ ) = - ZZ( NN + 1 - IQ, 1 )*EXPBEA( 0 )
c                                                  ** Bottom boundary

               B( NCOL - NN + IQ ) = -ZZ( IQ + NN, NCUT )*EXPBEA( NCUT )

   10       CONTINUE


         ELSE

            DO 30 IQ = 1, NN

               B( IQ ) = - ZZ( NN + 1 - IQ, 1 )*EXPBEA( 0 )

               SUM  = 0.
               DO 20 JQ = 1, NN
                  SUM  = SUM + CWT( JQ )*CMU( JQ )*BDR( IQ, JQ )*
     &                         ZZ( NN + 1 - JQ, NCUT )*EXPBEA( NCUT )
   20          CONTINUE

               B( NCOL - NN + IQ ) = SUM
               IF( FBEAM.GT.0.0 ) B( NCOL - NN + IQ ) = SUM +
     &             ( BDR( IQ,0 )*UMU0*FBEAM / PI - ZZ( IQ + NN,NCUT ) )*
     &             EXPBEA( NCUT )

   30       CONTINUE

         END IF
c                             ** Continuity condition for layer
c                                interfaces of Eq. STWJ(20b)
         IT   = NN

         DO 50 LC = 1, NCUT - 1

            DO 40 IQ = 1, NSTR
               IT   = IT + 1
               B( IT ) = ( ZZ( IQ, LC+1 ) - ZZ( IQ, LC ) )*EXPBEA( LC )
   40       CONTINUE

   50    CONTINUE


      ELSE
c                                   ** Azimuth-independent case

         IF( FBEAM.EQ.0.0 ) THEN

            DO 60 IQ = 1, NN
c                                      ** Top boundary

               B( IQ ) = -ZPLK0( NN + 1 - IQ, 1 ) + FISOT + TPLANK

   60       CONTINUE


            IF( LYRCUT ) THEN
c                               ** No intensity component for truncated
c                                  bottom layer
               DO 70 IQ = 1, NN
c                                      ** Bottom boundary

                  B( NCOL - NN + IQ ) = - ZPLK0( IQ + NN, NCUT ) -
     &                                    ZPLK1( IQ + NN, NCUT )*
     &                                    TAUCPR( NCUT )
   70          CONTINUE


            ELSE

               DO 90 IQ = 1, NN

                  SUM  = 0.
                  DO 80 JQ = 1, NN
                     SUM  = SUM + CWT( JQ )*CMU( JQ )*BDR( IQ, JQ )*
     &                            ( ZPLK0( NN + 1 - JQ,NCUT ) +
     &                        ZPLK1( NN + 1 - JQ,NCUT )*TAUCPR( NCUT ) )
   80             CONTINUE

                  B( NCOL - NN + IQ ) = 2.*SUM + BEM( IQ )*BPLANK -
     &                                  ZPLK0( IQ + NN, NCUT ) -
     &                                  ZPLK1( IQ + NN, NCUT )*
     &                                  TAUCPR( NCUT )
   90          CONTINUE

            END IF
c                             ** Continuity condition for layer
c                                interfaces, STWJ(20b)
            IT   = NN
            DO 110 LC = 1, NCUT - 1

               DO 100 IQ = 1, NSTR
                  IT   = IT + 1
                  B( IT ) =   ZPLK0( IQ, LC + 1 ) - ZPLK0( IQ, LC ) +
     &                      ( ZPLK1( IQ, LC + 1 ) - ZPLK1( IQ, LC ) )*
     &                      TAUCPR( LC )
  100          CONTINUE

  110       CONTINUE


         ELSE

            DO 120 IQ = 1, NN
               B( IQ ) = - ZZ( NN + 1 - IQ, 1 )*EXPBEA( 0 ) -
     &                   ZPLK0( NN + 1 - IQ, 1 ) + FISOT + TPLANK
  120       CONTINUE

            IF( LYRCUT ) THEN

               DO 130 IQ = 1, NN
                  B(NCOL-NN+IQ) = - ZZ(IQ+NN, NCUT) * EXPBEA(NCUT)
     &                            - ZPLK0(IQ+NN, NCUT)
     &                            - ZPLK1(IQ+NN, NCUT) * TAUCPR(NCUT)
  130          CONTINUE


            ELSE

               DO 150 IQ = 1, NN

                  SUM  = 0.
                  DO 140 JQ = 1, NN
                     SUM = SUM + CWT(JQ) * CMU(JQ) * BDR(IQ,JQ)
     &                          * ( ZZ(NN+1-JQ, NCUT) * EXPBEA(NCUT)
     &                            + ZPLK0(NN+1-JQ, NCUT)
     &                            + ZPLK1(NN+1-JQ, NCUT) * TAUCPR(NCUT))
  140             CONTINUE

                  B(NCOL-NN+IQ) = 2.*SUM + ( BDR(IQ,0) * UMU0*FBEAM/PI
     &                                - ZZ(IQ+NN, NCUT) ) * EXPBEA(NCUT)
     &                            + BEM(IQ) * BPLANK
     &                            - ZPLK0(IQ+NN, NCUT)
     &                            - ZPLK1(IQ+NN, NCUT) * TAUCPR(NCUT)
  150          CONTINUE

            END IF


            IT   = NN

            DO 170 LC = 1, NCUT - 1

               DO 160 IQ = 1, NSTR

                  IT   = IT + 1
                  B(IT) = ( ZZ(IQ,LC+1) - ZZ(IQ,LC) ) * EXPBEA(LC)
     &                    + ZPLK0(IQ,LC+1) - ZPLK0(IQ,LC) +
     &                    ( ZPLK1(IQ,LC+1) - ZPLK1(IQ,LC) ) * TAUCPR(LC)
  160          CONTINUE

  170       CONTINUE

         END IF

      END IF
c                     ** Find L-U (lower/upper triangular) decomposition
c                        of band matrix CBAND and test if it is nearly
c                        singular (note: CBAND is destroyed)
c                        (CBAND is in LINPACK packed format)
      RCOND  = 0.0
      NCD    = 3*NN - 1

      CALL SGBCO( CBAND, MI9M2, NCOL, NCD, NCD, IPVT, RCOND, Z )

      IF( 1.0 + RCOND.EQ.1.0 )
     &    CALL ERRMSG('SOLVE0--SGBCO says matrix near singular',.FALSE.)

c                   ** Solve linear system with coeff matrix CBAND
c                      and R.H. side(s) B after CBAND has been L-U
c                      decomposed.  Solution is returned in B.

      CALL SGBSL( CBAND, MI9M2, NCOL, NCD, NCD, IPVT, B, 0 )

c                   ** Zero CBAND (it may contain 'foreign'
c                      elements upon returning from LINPACK);
c                      necessary to prevent errors

      CALL ZEROIT( CBAND, MI9M2*NNLYRI )

      DO 190 LC = 1, NCUT

         IPNT  = LC*NSTR - NN

         DO 180 IQ = 1, NN
            LL( NN + 1 - IQ, LC ) = B( IPNT + 1 - IQ )
            LL( IQ + NN,     LC ) = B( IQ + IPNT )
  180    CONTINUE

  190 CONTINUE

      RETURN
      END

      SUBROUTINE SURFAC( ALBEDO, DELM0, FBEAM, HLPR, LAMBER, MI, MAZIM,
     &                   MXCMU, MXUMU, NN, NUMU, NSTR, ONLYFL, UMU,
     &                   USRANG, YLM0, YLMC, YLMU, BDR, EMU, BEM, RMU )

c       Specifies user's surface bidirectional properties, STWJ(21)

c   I N P U T     V A R I A B L E S:

c       DELM0  :  Kronecker delta, delta-sub-m0
c       HLPR   :  Legendre moments of surface bidirectional reflectivity
c                    (with 2K+1 factor included)
c       MAZIM  :  Order of azimuthal component
c       NN     :  Order of double-Gauss quadrature (NSTR/2)
c       YLM0   :  Normalized associated Legendre polynomial
c                 at the beam angle
c       YLMC   :  Normalized associated Legendre polynomials
c                 at the quadrature angles
c       YLMU   :  Normalized associated Legendre polynomials
c                 at the user angles
c       (remainder are DISORT input variables)

c    O U T P U T     V A R I A B L E S:

c       BDR :  Surface bidirectional reflectivity (computational angles)
c       RMU :  Surface bidirectional reflectivity (user angles)
c       BEM :  Surface directional emissivity (computational angles)
c       EMU :  Surface directional emissivity (user angles)

c    I N T E R N A L     V A R I A B L E S:

c       DREF      Directional reflectivity
c       NMUG   :  Number of angle cosine quadrature points on (0,1) for
c                   integrating bidirectional reflectivity to get
c                   directional emissivity (it is necessary to use a
c                   quadrature set distinct from the computational
c                   angles, because the computational angles may not be
c                   dense enough--NSTR may be too small--to give an
c                   accurate approximation for the integration).
c       GMU    :  The NMUG angle cosine quadrature points on (0,1)
c       GWT    :  The NMUG angle cosine quadrature weights on (0,1)
c       YLMG   :  Normalized associated Legendre polynomials
c                   at the NMUG quadrature angles

c   Called by- DISORT
c   Calls- QGAUSN, LEPOLY, ZEROIT, ERRMSG
c +-------------------------------------------------------------------+

c     .. Parameters ..

      INTEGER   NMUG, MAXSTR
      PARAMETER ( NMUG = 10, MAXSTR = 100 )
c     ..
c     .. Scalar Arguments ..

      LOGICAL   LAMBER, ONLYFL, USRANG
      INTEGER   MAZIM, MI, MXCMU, MXUMU, NN, NSTR, NUMU
      REAL      ALBEDO, DELM0, FBEAM
c     ..
c     .. Array Arguments ..

      REAL      BDR( MI, 0:MI ), BEM( MI ), EMU( MXUMU ),
     &          HLPR( 0:MXCMU ), RMU( MXUMU, 0:MI ), UMU( * ),
     &          YLM0( 0:MXCMU ), YLMC( 0:MXCMU, MXCMU ),
     &          YLMU( 0:MXCMU, MXUMU )
c     ..
c     .. Local Scalars ..

      LOGICAL   PASS1
      INTEGER   IQ, IU, JG, JQ, K
      REAL      DREF, SGN, SUM
c     ..
c     .. Local Arrays ..

      REAL      GMU( NMUG ), GWT( NMUG ), YLMG( 0:MAXSTR, NMUG )
c     ..
c     .. External Subroutines ..

      EXTERNAL  ERRMSG, LEPOLY, QGAUSN, ZEROIT
c     ..
      SAVE      PASS1, GMU, GWT, YLMG
      DATA      PASS1 / .TRUE. /


      IF( PASS1 ) THEN

         PASS1  = .FALSE.

         CALL QGAUSN( NMUG, GMU, GWT )

         CALL LEPOLY( NMUG, 0, MAXSTR, MAXSTR, GMU, YLMG )

c                       ** Convert Legendre polys. to negative GMU
         SGN  = - 1.0

         DO 20 K = 0, MAXSTR

            SGN  = - SGN

            DO 10 JG = 1, NMUG
               YLMG( K, JG ) = SGN*YLMG( K, JG )
   10       CONTINUE

   20    CONTINUE

      END IF


      CALL ZEROIT( BDR, MI*( MI + 1 ) )
      CALL ZEROIT( BEM, MI )

      IF( LAMBER .AND. MAZIM.EQ.0 ) THEN

         DO 40 IQ = 1, NN

            BEM( IQ ) = 1.- ALBEDO

            DO 30 JQ = 0, NN
               BDR( IQ, JQ ) = ALBEDO
   30       CONTINUE

   40    CONTINUE


      ELSE IF( .NOT.LAMBER ) THEN
c                                  ** Compute surface bidirectional
c                                     properties at computational angles
         DO 80 IQ = 1, NN

            DO 60 JQ = 1, NN

               SUM  = 0.0
               DO 50 K = MAZIM, NSTR - 1
                  SUM  = SUM + HLPR( K )*YLMC( K, IQ )*
     &                         YLMC( K, JQ + NN )
   50          CONTINUE

               BDR( IQ, JQ ) = ( 2.- DELM0 )*SUM

   60       CONTINUE


            IF( FBEAM.GT.0.0 ) THEN

               SUM  = 0.0
               DO 70 K = MAZIM, NSTR - 1
                  SUM  = SUM + HLPR( K )*YLMC( K, IQ )*YLM0( K )
   70          CONTINUE

               BDR( IQ, 0 ) = ( 2.- DELM0 )*SUM

            END IF

   80    CONTINUE


         IF( MAZIM.EQ.0 ) THEN

            IF( NSTR.GT.MAXSTR )
     &          CALL ERRMSG('SURFAC--parameter MAXSTR too small',.True.)

c                              ** Integrate bidirectional reflectivity
c                                 at reflection polar angles CMU and
c                                 incident angles GMU to get
c                                 directional emissivity at
c                                 computational angles CMU.
            DO 110 IQ = 1, NN

               DREF  = 0.0

               DO 100 JG = 1, NMUG

                  SUM  = 0.0
                  DO 90 K = 0, NSTR - 1
                     SUM  = SUM + HLPR( K )*YLMC( K, IQ )*
     &                            YLMG( K, JG )
   90             CONTINUE

                  DREF  = DREF + 2.*GWT( JG )*GMU( JG )*SUM

  100          CONTINUE

               BEM( IQ ) = 1.- DREF

  110       CONTINUE

         END IF

      END IF
c                                       ** Compute surface bidirectional
c                                          properties at user angles

      IF( .NOT.ONLYFL .AND. USRANG ) THEN

         CALL ZEROIT( EMU, MXUMU )
         CALL ZEROIT( RMU, MXUMU*( MI + 1 ) )

         DO 180 IU = 1, NUMU

            IF( UMU( IU ).GT.0.0 ) THEN

               IF( LAMBER .AND. MAZIM.EQ.0 ) THEN

                  DO 120 IQ = 0, NN
                     RMU( IU, IQ ) = ALBEDO
  120             CONTINUE

                  EMU( IU ) = 1.- ALBEDO


               ELSE IF( .NOT.LAMBER ) THEN

                  DO 140 IQ = 1, NN

                     SUM  = 0.0
                     DO 130 K = MAZIM, NSTR - 1
                        SUM  = SUM + HLPR( K )*YLMU( K, IU )*
     &                               YLMC( K, IQ + NN )
  130                CONTINUE

                     RMU( IU, IQ ) = ( 2.- DELM0 )*SUM

  140             CONTINUE


                  IF( FBEAM.GT.0.0 ) THEN

                     SUM  = 0.0
                     DO 150 K = MAZIM, NSTR - 1
                        SUM  = SUM + HLPR( K )*YLMU( K, IU )*YLM0( K )
  150                CONTINUE

                     RMU( IU, 0 ) = ( 2.- DELM0 )*SUM

                  END IF


                  IF( MAZIM.EQ.0 ) THEN

c                               ** Integrate bidirectional reflectivity
c                                  at reflection angles UMU and
c                                  incident angles GMU to get
c                                  directional emissivity at
c                                  user angles UMU.
                     DREF  = 0.0

                     DO 170 JG = 1, NMUG

                        SUM  = 0.0
                        DO 160 K = 0, NSTR - 1
                           SUM  = SUM + HLPR( K )*YLMU( K, IU )*
     &                                  YLMG( K, JG )
  160                   CONTINUE

                        DREF  = DREF + 2.*GWT( JG )*GMU( JG )*SUM

  170                CONTINUE

                     EMU( IU ) = 1.- DREF

                  END IF

               END IF

            END IF

  180    CONTINUE

      END IF

      END


*bm  SOLVEC calls SOLEIG and UPBEAM; if UPBEAM reports a potenially 
*bm  unstable solution, the calculation is repeated with a slightly 
*bm  changed single scattering albedo; this process is iterates 
*bm  until a stable solution is found; as stable solutions may be 
*bm  reached either by increasing or by decreasing the single 
*bm  scattering albedo, both directions are explored ('upward' and
*bm  'downward' iteration); the solution which required the smaller 
*bm  change in the single scattering albedo is finally returned 
*bm  by SOLVEC.

      SUBROUTINE SOLVEC( AMB, APB, ARRAY, CMU, CWT, GL, MI,
     &     MAZIM, MXCMU, NN, NSTR, YLM0, YLMC, CC, 
     &     EVECC, EVAL, KK, GC, AAD, EVECCD, EVALD,
     &     WK, WKD, DELM0, FBEAM, IPVT, PI, UMU0, ZJ, ZZ,
     &     OPRIM, LC, DITHER, mu2, glsave, dgl)

cgy added glsave and dgl to call to allow adjustable dimensioning


c     .. Scalar Arguments ..

      INTEGER   MAZIM, MI, MXCMU, NN, NSTR, LC
      REAL      DELM0, FBEAM, PI, UMU0, OPRIM, DITHER
      REAL      mu2

c     ..
c     .. Array Arguments ..

      INTEGER   IPVT( * )
      
      REAL      AMB( MI, MI ), APB( MI, MI ), ARRAY( MI, * ),
     &     CC( MXCMU, MXCMU ), CMU( MXCMU ), CWT( MXCMU ),
     &     EVAL( MI ), EVECC( MXCMU, MXCMU ), GC( MXCMU, MXCMU ),
     &     GL( 0:MXCMU ), KK( MXCMU ), 
     &     YLM0( 0:MXCMU ), YLMC( 0:MXCMU, MXCMU ),
     &     WK( MXCMU ), ZJ( MXCMU ), ZZ( MXCMU )

      DOUBLE PRECISION AAD( MI, MI ), EVALD( MI ), EVECCD( MI, MI ),
     &                 WKD( MXCMU )

*bm   Variables for instability fix
      
      INTEGER UAGAIN, DAGAIN
      REAL MINRCOND, ADD, UADD, DADD, SSA, DSSA, FACTOR
      REAL GLSAVE( 0:MXCMU ), DGL( 0:MXCMU )
      
      LOGICAL  DONE, NOUP, NODN, DEBUG, INSTAB
      
*bm   reset parameters

      DONE = .FALSE.
      NOUP = .FALSE.
      NODN = .FALSE.


*bm   flag for printing debugging output      
*      DEBUG  = .TRUE.
      DEBUG  = .FALSE.

*bm   instability parameter; the solution is considered 
*bm   unstable, if the RCOND reported by SGECO is smaller 
*bm   than MINRCOND
      MINRCOND = 5000. * R1MACH(4)

*bm   if an instability is detected, the single scattering albedo
*bm   is iterated downwards in steps of DADD and upwards in steps 
*bm   of UADD; in practice, MINRCOND and -MINRCOND should 
*bm   be reasonable choices for these parameters
      DADD    = -MINRCOND
      UADD    = MINRCOND

      UAGAIN = 0
      DAGAIN = 0
      ADD   = DADD
      

*bm   save array GL( ) because it will be 
*bm   changed if an iteration should be neccessary
      DO K = MAZIM, NSTR - 1
         GLSAVE( K ) =  GL( K )
      ENDDO
      
      SSA = OPRIM


*bm   in case of an instability reported by UPBEAM (INSTAB)
*bm   the single scattering albedo will be changed by a small 
*bm   amount (ADD); this is indicated by DAGAIN or UAGAIN 
*bm   being larger than 0; a change in the single scattering 
*bm   albedo is equivalent to scaling the array GL( )

 666  IF ( DAGAIN .GT. 0 .OR. UAGAIN .GT. 0)  THEN
         FACTOR = (SSA + ADD) / SSA
         DO K = MAZIM, NSTR - 1
            GL( K ) =  GL( K ) * FACTOR
         ENDDO

         SSA = SSA + ADD
         
*bm   if the single scattering albedo is now smaller than 0
*bm   the downward iteration is stopped and upward iteration 
*bm   is forced instead

         IF( SSA .LT. DITHER) THEN
            NODN = .TRUE.
            DAGAIN = -1
            goto 778
         ENDIF

*bm   if the single scattering albedo is now larger than its maximum 
*bm   allowed value (1.0 - DITHER), the upward iteration is 
*bm   stopped and downward iteration is forced instead

         IF( SSA .GT. 1.0 - DITHER) THEN
            NOUP = .TRUE.
            UAGAIN = -1
            goto 888
         ENDIF
      ENDIF


c     ** Solve eigenfunction problem in Eq. STWJ(8B);
c        return eigenvalues and eigenvectors

 777     CALL SOLEIG( AMB, APB, ARRAY, CMU, CWT, GL, MI,
     &     MAZIM, MXCMU, NN, NSTR, YLMC, CC, EVECC, EVAL,
     &     KK, GC, AAD, EVECCD, EVALD,
     &     WKD )

c     ** Calculate particular solutions of
c        q.SS(18) for incident beam source

      IF ( FBEAM.GT.0.0 ) THEN
         CALL  UPBEAM( mu2,
     $        ARRAY, CC, CMU, DELM0, FBEAM, GL,
     $        IPVT, MAZIM, MXCMU, NN, NSTR, PI, UMU0, WK,
     $        YLM0, YLMC, ZJ, ZZ, MINRCOND, INSTAB)
      ENDIF
      
c     ** Calculate particular solutions of
c        Eq. SS(15) for thermal emission source
c        (not available in psndo.f)
      
*bm   finished if the result is stable on the first try
      IF ( (.NOT. INSTAB) .AND. 
     $     (UAGAIN .EQ. 0) .AND. (DAGAIN .EQ. 0)) THEN
         goto 999
      ENDIF

*bm   downward iteration
      IF( INSTAB .AND. UAGAIN .EQ. 0 )  THEN
         DAGAIN = DAGAIN + 1
         GOTO 666
      ENDIF
      
*bm   upward iteration
      IF( INSTAB .AND. UAGAIN .GT. 0 )  THEN
         UAGAIN = UAGAIN + 1
         GOTO 666
      ENDIF


*bm   ( DAGAIN .NE. 0 ) at this place means that the downward
*bm   iteration is finished 

 778  IF (DAGAIN .NE. 0 .AND. UAGAIN .EQ. 0) THEN
         
*bm   save downward iteration data for later use and 
*bm   restore original input data
         DO K = MAZIM, NSTR - 1
            DGL( K ) =  GL( K )
            GL( K ) =  GLSAVE( K )
         ENDDO

         DSSA = SSA
         SSA = OPRIM

*bm   start upward iteration
         ADD = UADD
         UAGAIN = UAGAIN + 1
         GOTO 666
      ENDIF

*bm   both iterations finished
 888  IF (DONE) THEN
         goto 998
      ENDIF


*bm  if neither upward nor downward iteration converged, the 
*bm  original conditions are restored and SOLEIG/UPBEAM 
*bm  is called for the last time 
         
      IF (NOUP .AND. NODN) THEN
         
         DO K = MAZIM, NSTR - 1
            GL( K ) =  GLSAVE( K )
         ENDDO
         
         SSA = OPRIM
         
         IF (DEBUG) THEN
            write (*,*) '! *** Neither upward nor downward iteration'
            write (*,*) '! *** converged; using original result.'
         ENDIF

         DONE = .TRUE.
         GOTO 777
      ENDIF

*bm  if upward iteration did not converge, the stable downward conditions
*bm  are restored and SOLEIG/UPBEAM is called for the last time
      IF (NOUP) THEN
         DO K = MAZIM, NSTR - 1
            GL( K ) =  DGL( K )
         ENDDO
         
         SSA = DSSA
         
         IF (DEBUG) THEN
            write (*,*) '! *** The upward iteration did not converge.'
            write (*,*) '! *** Had to iterate ', DAGAIN,
     $           ' times in layer LC =', LC,';'
            write (*,*) '! *** changed SSA from ',
     $           OPRIM, ' to ', SSA,','
            write (*,*) '! *** by a factor of ', SSA/OPRIM
         ENDIF

         DONE = .TRUE.
         GOTO 777
      ENDIF

*bm  if downward iteration did not converge, we are done 
*bm  (the result of the upward iteration will be used)
      IF (NODN) THEN
         IF (DEBUG) THEN
            write (*,*) '! *** The downward iteration did not converge.'
            write (*,*) '! *** Had to iterate ', UAGAIN,
     $           ' times in layer LC =', LC,';'
            write (*,*) '! *** changed SSA from ',
     $           OPRIM, ' to ', SSA,','
            write (*,*) '! *** by a factor of ', SSA/OPRIM
         ENDIF
         
         DONE = .TRUE.
         GOTO 998
      ENDIF

      
*bm   if both iterations converged, and if the upward iteration 
*bm   required more steps than the downward iteration, the stable 
*bm   downward conditions are restored and SOLEIG/UPBEAM is 
*bm   called for the last time 
         
      IF (UAGAIN .GT. DAGAIN) THEN
         DO K = MAZIM, NSTR - 1
            GL( K ) =  DGL( K )
         ENDDO
         
         SSA = DSSA
         
         IF (DEBUG) THEN
            write (*,*) '! *** Both iterations converged;',
     $           ' using downward.'
            write (*,*) '! *** Had to iterate ', DAGAIN,
     $        ' times in layer LC =', LC,';'
            write (*,*) '! *** changed SSA from ',
     $           OPRIM, ' to ', SSA,','
            write (*,*) '! *** by a factor of ', SSA/OPRIM
         ENDIF

         DONE = .TRUE.
         GOTO 777
      ELSE
         
         IF (DEBUG) THEN
            write (*,*) '! *** Both iterations converged;',
     $           ' using upward.'
            write (*,*) '! *** Had to iterate ', UAGAIN,
     $        ' times in layer LC =', LC,';'
            write (*,*) '! *** changed SSA from ',
     $           OPRIM, ' to ', SSA,','
            write (*,*) '! *** by a factor of ', SSA/OPRIM
         ENDIF

         DONE = .TRUE.
         goto 998
      ENDIF
      
*bm   finally restore original input data
 998  DO K = MAZIM, NSTR - 1
         GL( K ) =  GLSAVE( K )
      ENDDO
      
 999  CONTINUE
      END



      SUBROUTINE UPBEAM( mu2,
     &                   ARRAY, CC, CMU, DELM0, FBEAM, GL, IPVT, MAZIM,
     &                   MXCMU, NN, NSTR, PI, UMU0, WK, YLM0, YLMC, ZJ,
     &                   ZZ, MINRCOND, INSTAB )

c         Finds the incident-beam particular solution of SS(18)

c   I N P U T    V A R I A B L E S:

c       CC     :  C-sub-ij in Eq. SS(5)
c       CMU    :  Abscissae for Gauss quadrature over angle cosine
c       DELM0  :  Kronecker delta, delta-sub-m0
c       GL     :  Delta-M scaled Legendre coefficients of phase function
c                    (including factors 2L+1 and single-scatter albedo)
c       MAZIM  :  Order of azimuthal component
c       YLM0   :  Normalized associated Legendre polynomial
c                    at the beam angle
c       YLMC   :  Normalized associated Legendre polynomial
c                    at the quadrature angles
c       (remainder are DISORT input variables)

c   O U T P U T    V A R I A B L E S:

c       ZJ     :  Right-hand side vector X-sub-zero in SS(19); also the
c                 solution vector Z-sub-zero after solving that system

c       ZZ     :  Permanent storage for ZJ, but re-ordered

c   I N T E R N A L    V A R I A B L E S:

c       ARRAY  :  Coefficient matrix in left-hand side of Eq. SS(19)
c       IPVT   :  Integer vector of pivot indices required by LINPACK
c       WK     :  Scratch array required by LINPACK

c   Called by- DISORT
c   Calls- SGECO, ERRMSG, SGESL
c +-------------------------------------------------------------------+


c     .. Scalar Arguments ..

      INTEGER   MAZIM, MXCMU, NN, NSTR
      LOGICAL   INSTAB
      REAL      MINRCOND
      REAL      DELM0, FBEAM, PI, UMU0
      REAL mu2
c     ..
c     .. Array Arguments ..

      INTEGER   IPVT( * )
      REAL      ARRAY( MXCMU, MXCMU ), CC( MXCMU, MXCMU ), CMU( MXCMU ),
     &          GL( 0:MXCMU ), WK( MXCMU ), YLM0( 0:MXCMU ),
     &          YLMC( 0:MXCMU, * ), ZJ( MXCMU ), ZZ( MXCMU )
c     ..
c     .. Local Scalars ..

      INTEGER   IQ, JOB, JQ, K
      REAL      RCOND, SUM
c     ..
c     .. External Subroutines ..

      EXTERNAL  ERRMSG, SGECO, SGESL
c     ..


      DO 30 IQ = 1, NSTR

         DO 10 JQ = 1, NSTR
            ARRAY( IQ, JQ ) = -CC( IQ, JQ )
   10    CONTINUE

         ARRAY( IQ, IQ ) = 1.+ CMU( IQ ) / mu2 + ARRAY( IQ, IQ )

         SUM  = 0.
         DO 20 K = MAZIM, NSTR - 1
            SUM  = SUM + GL( K )*YLMC( K, IQ )*YLM0( K )
   20    CONTINUE

         ZJ( IQ ) = ( 2.- DELM0 )*FBEAM*SUM / ( 4.*PI )
   30 CONTINUE

c                  ** Find L-U (lower/upper triangular) decomposition
c                     of ARRAY and see if it is nearly singular
c                     (NOTE:  ARRAY is destroyed)
      RCOND  = 0.0

      CALL SGECO( ARRAY, MXCMU, NSTR, IPVT, RCOND, WK )

*bm      IF( 1.0 + RCOND.EQ.1.0 )
*bm     &    CALL ERRMSG('UPBEAM--SGECO says matrix near singular',.FALSE.)
*bm
*bm   replaced original check of RCOND by the following:

      INSTAB = .FALSE.
      IF( ABS(RCOND) .LT. MINRCOND )  THEN
         INSTAB = .TRUE.
         RETURN
      ENDIF

c                ** Solve linear system with coeff matrix ARRAY
c                   (assumed already L-U decomposed) and R.H. side(s)
c                   ZJ;  return solution(s) in ZJ
      JOB  = 0

      CALL SGESL( ARRAY, MXCMU, NSTR, IPVT, ZJ, JOB )

CDIR$ IVDEP
      DO 40 IQ = 1, NN
         ZZ( IQ + NN )     = ZJ( IQ )
         ZZ( NN + 1 - IQ ) = ZJ( IQ + NN )
   40 CONTINUE

      END


      SUBROUTINE ZEROAL( ND1, EXPBEA, FLYR, OPRIM, TAUCPR, XR0, XR1,
     &                    ND2, CMU, CWT, PSI, WK, Z0, Z1, ZJ,
     &                    ND3, HLPR, YLM0,
     &                    ND4, ARRAY, CC, EVECC,
     &                    ND5, GL,
     &                    ND6, YLMC,
     &                    ND7, YLMU,
     &                    ND8, KK, LL, ZZ, ZPLK0, ZPLK1,
     &                    ND9, GC,
     &                    ND10, LAYRU, UTAUPR,
     &                    ND11, GU,
     &                    ND12, Z0U, Z1U, ZBEAM,
     &                    ND13, EVAL,
     &                    ND14, AMB, APB,
     &                    ND15, IPVT, Z,
     &                    ND16, RFLDIR, RFLDN, FLUP, UAVG, DFDT,
     &                    ND17, ALBMED, TRNMED,
     &                    ND18, U0U,
     &                    ND19, UU )

c         ZERO ARRAYS; NDn is dimension of all arrays following
c         it in the argument list

c   Called by- DISORT
c --------------------------------------------------------------------

c     .. Scalar Arguments ..

      INTEGER   ND1, ND10, ND11, ND12, ND13, ND14, ND15, ND16, ND17,
     &          ND18, ND19, ND2, ND3, ND4, ND5, ND6, ND7, ND8, ND9
c     ..
c     .. Array Arguments ..

      INTEGER   IPVT( * ), LAYRU( * )
      REAL      ALBMED( * ), AMB( * ), APB( * ), ARRAY( * ), CC( * ),
     &          CMU( * ), CWT( * ), DFDT( * ), EVAL( * ), EVECC( * ),
     &          EXPBEA( * ), FLUP( * ), FLYR( * ), GC( * ), GL( * ),
     &          GU( * ), HLPR( * ), KK( * ), LL( * ), OPRIM( * ),
     &          PSI( * ), RFLDIR( * ), RFLDN( * ), TAUCPR( * ),
     &          TRNMED( * ), U0U( * ), UAVG( * ), UTAUPR( * ), UU( * ),
     &          WK( * ), XR0( * ), XR1( * ), YLM0( * ), YLMC( * ),
     &          YLMU( * ), Z( * ), Z0( * ), Z0U( * ), Z1( * ), Z1U( * ),
     &          ZBEAM( * ), ZJ( * ), ZPLK0( * ), ZPLK1( * ), ZZ( * )
c     ..
c     .. Local Scalars ..

      INTEGER   N
c     ..


      DO 10 N = 1, ND1
         EXPBEA( N ) = 0.0
         FLYR( N )   = 0.0
         OPRIM( N )  = 0.0
         TAUCPR( N ) = 0.0
         XR0( N )    = 0.0
         XR1( N )    = 0.0
   10 CONTINUE

      DO 20 N = 1, ND2
         CMU( N ) = 0.0
         CWT( N ) = 0.0
         PSI( N ) = 0.0
         WK( N )  = 0.0
         Z0( N )  = 0.0
         Z1( N )  = 0.0
         ZJ( N )  = 0.0
   20 CONTINUE

      DO 30 N = 1, ND3
         HLPR( N ) = 0.0
         YLM0( N ) = 0.0
   30 CONTINUE

      DO 40 N = 1, ND4
         ARRAY( N ) = 0.0
         CC( N )    = 0.0
         EVECC( N ) = 0.0
   40 CONTINUE

      DO 50 N = 1, ND5
         GL( N ) = 0.0
   50 CONTINUE

      DO 60 N = 1, ND6
         YLMC( N ) = 0.0
   60 CONTINUE

      DO 70 N = 1, ND7
         YLMU( N ) = 0.0
   70 CONTINUE

      DO 80 N = 1, ND8
         KK( N )    = 0.0
         LL( N )    = 0.0
         ZZ( N )    = 0.0
         ZPLK0( N ) = 0.0
         ZPLK1( N ) = 0.0
   80 CONTINUE

      DO 90 N = 1, ND9
         GC( N ) = 0.0
   90 CONTINUE

      DO 100 N = 1, ND10
         LAYRU( N )  = 0
         UTAUPR( N ) = 0.0
  100 CONTINUE

      DO 110 N = 1, ND11
         GU( N ) = 0.0
  110 CONTINUE

      DO 120 N = 1, ND12
         Z0U( N )   = 0.0
         Z1U( N )   = 0.0
         ZBEAM( N ) = 0.0
  120 CONTINUE

      DO 130 N = 1, ND13
         EVAL( N ) = 0.0
  130 CONTINUE

      DO 140 N = 1, ND14
         AMB( N ) = 0.0
         APB( N ) = 0.0
  140 CONTINUE

      DO 150 N = 1, ND15
         IPVT( N ) = 0
         Z( N )    = 0.0
  150 CONTINUE

      DO 160 N = 1, ND16
         RFLDIR( N ) = 0.
         RFLDN( N )  = 0.
         FLUP( N )   = 0.
         UAVG( N )   = 0.
         DFDT( N )   = 0.
  160 CONTINUE

      DO 170 N = 1, ND17
         ALBMED( N ) = 0.
         TRNMED( N ) = 0.
  170 CONTINUE

      DO 180 N = 1, ND18
         U0U( N ) = 0.
  180 CONTINUE

      DO 190 N = 1, ND19
         UU( N ) = 0.
  190 CONTINUE


      END

      SUBROUTINE ZEROIT( A, LENGTH )

c         Zeros a real array A having LENGTH elements
c --------------------------------------------------------------------

c     .. Scalar Arguments ..

      INTEGER   LENGTH
c     ..
c     .. Array Arguments ..

      REAL      A( LENGTH )
c     ..
c     .. Local Scalars ..

      INTEGER   L
c     ..

      DO 10 L = 1, LENGTH
         A( L ) = 0.0
   10 CONTINUE

      END

      REAL FUNCTION DREF( MU, HL, NSTR )

c        Exact flux albedo for given angle of incidence, given
c        a bidirectional reflectivity characterized by its
c        Legendre coefficients ( NOTE** these will only agree
c        with bottom-boundary albedos calculated by DISORT in
c        the limit as number of streams go to infinity, because
c        DISORT evaluates the integral 'CL' only approximately,
c        by quadrature, while this routine calculates it exactly.)

c  INPUT :   MU     Cosine of incidence angle
c            HL     Legendre coefficients of bidirectional reflectivity
c          NSTR     Number of elements of HL to consider

c  INTERNAL VARIABLES (P-sub-L is the L-th Legendre polynomial) :

c       CL      Integral from 0 to 1 of  MU * P-sub-L(MU)
c                   (vanishes for  L = 3, 5, 7, ... )
c       PL      P-sub-L
c       PLM1    P-sub-(L-1)
c       PLM2    P-sub-(L-2)

c   Called by- CHEKIN
c   Calls- ERRMSG
c +-------------------------------------------------------------------+

c     .. Parameters ..

      INTEGER   MAXTRM
      PARAMETER ( MAXTRM = 100 )
c     ..
c     .. Scalar Arguments ..

      INTEGER   NSTR
      REAL      MU
c     ..
c     .. Array Arguments ..

      REAL      HL( 0:NSTR )
c     ..
c     .. Local Scalars ..

      LOGICAL   PASS1
      INTEGER   L
      REAL      CL, PL, PLM1, PLM2
c     ..
c     .. Local Arrays ..

      REAL      C( MAXTRM )
c     ..
c     .. External Subroutines ..

      EXTERNAL  ERRMSG
c     ..
c     .. Intrinsic Functions ..

      INTRINSIC MOD
c     ..
      SAVE      PASS1, C
      DATA      PASS1 / .TRUE. /
c     ..


      IF( PASS1 ) THEN

         PASS1  = .FALSE.
         CL     = 0.125
         C( 2 ) = 10.*CL

         DO 10 L = 4, MAXTRM, 2
            CL     = - CL*( L - 3 ) / ( L + 2 )
            C( L ) = 2.*( 2*L + 1 )*CL
   10    CONTINUE

      END IF


      IF( NSTR.LT.2 .OR. ABS(MU).GT.1.0 )
     &    CALL ERRMSG( 'DREF--input argument error(s)',.True. )

      IF( NSTR.GT.MAXTRM )
     &    CALL ERRMSG( 'DREF--parameter MAXTRM too small',.True. )


      DREF  = HL( 0 ) - 2.*HL( 1 )*MU
      PLM2  = 1.0
      PLM1  = - MU

      DO 20 L = 2, NSTR - 1
c                                ** Legendre polynomial recurrence

         PL = ( ( 2*L - 1 )*( -MU )*PLM1 - ( L-1 )*PLM2 ) / L

         IF( MOD( L,2 ).EQ.0 ) DREF   = DREF + C( L )*HL( L )*PL

         PLM2  = PLM1
         PLM1  = PL

   20 CONTINUE

      IF( DREF.LT.0.0 .OR. DREF.GT.1.0 )
     &    CALL ERRMSG( 'DREF--albedo value not in (0,1)',.False. )

      END

      REAL FUNCTION RATIO( A, B )

c        Calculate ratio  A/B  with over- and under-flow protection
c        (thanks to Prof. Jeff Dozier for some suggestions here).
c        Since this routine takes two logs, it is no speed demon,
c        but it is invaluable for comparing results from two runs
c        of a program under development.

c        NOTE:  In Fortran90, built-in functions TINY and HUGE
c               can replace the R1MACH calls.
c ---------------------------------------------------------------

c     .. Scalar Arguments ..

      REAL      A, B
c     ..
c     .. Local Scalars ..

      LOGICAL   PASS1
      REAL      ABSA, ABSB, HUGE, POWA, POWB, POWMAX, POWMIN, TINY
c     ..
c     .. External Functions ..

      REAL      R1MACH
      EXTERNAL  R1MACH
c     ..
c     .. Intrinsic Functions ..

      INTRINSIC ABS, LOG10, SIGN
c     ..
      SAVE      PASS1, TINY, HUGE, POWMAX, POWMIN
      DATA      PASS1 / .TRUE. /
c     ..


      IF( PASS1 ) THEN

         TINY   = R1MACH( 1 )
         HUGE   = R1MACH( 2 )
         POWMAX = LOG10( HUGE )
         POWMIN = LOG10( TINY )
         PASS1  = .FALSE.

      END IF


      IF( A.EQ.0.0 ) THEN

         IF( B.EQ.0.0 ) THEN

            RATIO  = 1.0

         ELSE

            RATIO  = 0.0

         END IF


      ELSE IF( B.EQ.0.0 ) THEN

         RATIO  = SIGN( HUGE, A )

      ELSE

         ABSA   = ABS( A )
         ABSB   = ABS( B )
         POWA   = LOG10( ABSA )
         POWB   = LOG10( ABSB )

         IF( ABSA.LT.TINY .AND. ABSB.LT.TINY ) THEN

            RATIO  = 1.0

         ELSE IF( POWA - POWB.GE.POWMAX ) THEN

            RATIO  = HUGE

         ELSE IF( POWA - POWB.LE.POWMIN ) THEN

            RATIO  = TINY

         ELSE

            RATIO  = ABSA / ABSB

         END IF
c                      ** DONT use old trick of determining sign
c                      ** from A*B because A*B may (over/under)flow

         IF( ( A.GT.0.0 .AND. B.LT.0.0 ) .OR.
     &       ( A.LT.0.0 .AND. B.GT.0.0 ) ) RATIO = -RATIO

      END IF

      END
      SUBROUTINE  ErrMsg( MESSAG, FATAL )

c        Print out a warning or error message;  abort if error
c        after making symbolic dump (machine-specific)

      LOGICAL       FATAL, MsgLim, Cray
      CHARACTER*(*) MESSAG
      INTEGER       MaxMsg, NumMsg
      SAVE          MaxMsg, NumMsg, MsgLim
      DATA NumMsg / 0 /,  MaxMsg / 100 /,  MsgLim / .FALSE. /


      IF ( FATAL )  THEN
         WRITE ( *, '(//,2A,//)' )  ' ******* ERROR >>>>>>  ', MESSAG
         STOP
      END IF

      NumMsg = NumMsg + 1
      IF( MsgLim )  RETURN

      IF ( NumMsg.LE.MaxMsg )  THEN
         WRITE ( *, '(/,2A,/)' )  ' ******* WARNING >>>>>>  ', MESSAG
      ELSE
         WRITE ( *,99 )
         MsgLim = .True.
      ENDIF

      RETURN

   99 FORMAT( //,' >>>>>>  TOO MANY WARNING MESSAGES --  ',
     $   'They will no longer be printed  <<<<<<<', // )
      END

      LOGICAL FUNCTION  WrtBad ( VarNam )

c          Write names of erroneous variables and return 'TRUE'

c      INPUT :   VarNam = Name of erroneous variable to be written
c                         ( CHARACTER, any length )

      CHARACTER*(*)  VarNam
      INTEGER        MaxMsg, NumMsg
      SAVE  NumMsg, MaxMsg
      DATA  NumMsg / 0 /,  MaxMsg / 50 /


      WrtBad = .TRUE.
      NumMsg = NumMsg + 1
      WRITE ( *, '(3A)' )  ' ****  Input variable  ', VarNam,
     $                     '  in error  ****'
      IF ( NumMsg.EQ.MaxMsg )
     $   CALL  ErrMsg ( 'Too many input errors.  Aborting...', .TRUE. )
      RETURN
      END

      LOGICAL FUNCTION  WrtDim ( DimNam, MinVal )

c          Write name of too-small symbolic dimension and
c          the value it should be increased to;  return 'TRUE'

c      INPUT :  DimNam = Name of symbolic dimension which is too small
c                        ( CHARACTER, any length )
c               Minval = Value to which that dimension should be
c                        increased (at least)

      CHARACTER*(*)  DimNam
      INTEGER        MinVal


      WRITE ( *, '(3A,I7)' )  ' ****  Symbolic dimension  ', DimNam,
     $                     '  should be increased to at least ', MinVal
      WrtDim = .TRUE.
      RETURN
      END

      LOGICAL FUNCTION  TstBad( VarNam, RelErr )

c       Write name (VarNam) of variable failing self-test and its
c       percent error from the correct value;  return  'FALSE'.

      CHARACTER*(*)  VarNam
      REAL           RelErr


      TstBad = .FALSE.
      WRITE( *, '(/,3A,1P,E11.2,A)' )
     $       ' Output variable ', VarNam,' differed by ', 100.*RelErr,
     $       ' per cent from correct value.  Self-test failed.'
      RETURN
      END
	SUBROUTINE  SGBCO( ABD, LDA, N, ML, MU, IPVT, RCOND, Z )

C         FACTORS A REAL BAND MATRIX BY GAUSSIAN ELIMINATION 
C         AND ESTIMATES THE CONDITION OF THE MATRIX.

C         REVISION DATE:  8/1/82
C         AUTHOR:  MOLER, C. B., (U. OF NEW MEXICO)

C     IF  RCOND  IS NOT NEEDED, SGBFA IS SLIGHTLY FASTER.
C     TO SOLVE  A*X = B , FOLLOW SBGCO BY SGBSL.

C     INPUT:

C        ABD     REAL(LDA, N)
C                CONTAINS THE MATRIX IN BAND STORAGE.  THE COLUMNS
C                OF THE MATRIX ARE STORED IN THE COLUMNS OF  ABD  AND
C                THE DIAGONALS OF THE MATRIX ARE STORED IN ROWS
C                ML+1 THROUGH 2*ML+MU+1 OF  ABD .
C                SEE THE COMMENTS BELOW FOR DETAILS.

C        LDA     INTEGER
C                THE LEADING DIMENSION OF THE ARRAY  ABD .
C                LDA MUST BE .GE. 2*ML + MU + 1 .

C        N       INTEGER
C                THE ORDER OF THE ORIGINAL MATRIX.

C        ML      INTEGER
C                NUMBER OF DIAGONALS BELOW THE MAIN DIAGONAL.
C                0 .LE. ML .LT. N .

C        MU      INTEGER
C                NUMBER OF DIAGONALS ABOVE THE MAIN DIAGONAL.
C                0 .LE. MU .LT. N .
C                MORE EFFICIENT IF  ML .LE. MU .

C     ON RETURN

C        ABD     AN UPPER TRIANGULAR MATRIX IN BAND STORAGE AND
C                THE MULTIPLIERS WHICH WERE USED TO OBTAIN IT.
C                THE FACTORIZATION CAN BE WRITTEN  A = L*U  WHERE
C                L  IS A PRODUCT OF PERMUTATION AND UNIT LOWER
C                TRIANGULAR MATRICES AND  U  IS UPPER TRIANGULAR.

C        IPVT    INTEGER(N)
C                AN INTEGER VECTOR OF PIVOT INDICES.

C        RCOND   REAL
C                AN ESTIMATE OF THE RECIPROCAL CONDITION OF  A .
C                FOR THE SYSTEM  A*X = B , RELATIVE PERTURBATIONS
C                IN  A  AND  B  OF SIZE  EPSILON  MAY CAUSE
C                RELATIVE PERTURBATIONS IN  X  OF SIZE  EPSILON/RCOND .
C                IF  RCOND  IS SO SMALL THAT THE LOGICAL EXPRESSION
C                           1.0 + RCOND .EQ. 1.0
C                IS TRUE, THEN  A  MAY BE SINGULAR TO WORKING
C                PRECISION.  IN PARTICULAR,  RCOND  IS ZERO  IF
C                EXACT SINGULARITY IS DETECTED OR THE ESTIMATE
C                UNDERFLOWS.

C        Z       REAL(N)
C                A WORK VECTOR WHOSE CONTENTS ARE USUALLY UNIMPORTANT.
C                IF  A  IS CLOSE TO A SINGULAR MATRIX, THEN  Z  IS
C                AN APPROXIMATE NULL VECTOR IN THE SENSE THAT
C                NORM(A*Z) = RCOND*NORM(A)*NORM(Z) .

C     BAND STORAGE

C           IF  A  IS A BAND MATRIX, THE FOLLOWING PROGRAM SEGMENT
C           WILL SET UP THE INPUT.

C                   ML = (BAND WIDTH BELOW THE DIAGONAL)
C                   MU = (BAND WIDTH ABOVE THE DIAGONAL)
C                   M = ML + MU + 1
C                   DO 20 J = 1, N
C                      I1 = MAX0(1, J-MU)
C                      I2 = MIN0(N, J+ML)
C                      DO 10 I = I1, I2
C                         K = I - J + M
C                         ABD(K,J) = A(I,J)
C                10    CONTINUE
C                20 CONTINUE

C           THIS USES ROWS  ML+1  THROUGH  2*ML+MU+1  OF  ABD .
C           IN ADDITION, THE FIRST  ML  ROWS IN  ABD  ARE USED FOR
C           ELEMENTS GENERATED DURING THE TRIANGULARIZATION.
C           THE TOTAL NUMBER OF ROWS NEEDED IN  ABD  IS  2*ML+MU+1 .
C           THE  ML+MU BY ML+MU  UPPER LEFT TRIANGLE AND THE
C           ML BY ML  LOWER RIGHT TRIANGLE ARE NOT REFERENCED.

C     EXAMPLE:  IF THE ORIGINAL MATRIX IS

C           11 12 13  0  0  0
C           21 22 23 24  0  0
C            0 32 33 34 35  0
C            0  0 43 44 45 46
C            0  0  0 54 55 56
C            0  0  0  0 65 66

C      THEN  N = 6, ML = 1, MU = 2, LDA .GE. 5  AND ABD SHOULD CONTAIN

C            *  *  *  +  +  +  , * = NOT USED
C            *  * 13 24 35 46  , + = USED FOR PIVOTING
C            * 12 23 34 45 56
C           11 22 33 44 55 66
C           21 32 43 54 65  *


C     ROUTINES CALLED:  FROM LINPACK: SGBFA
C                       FROM BLAS:    SAXPY, SDOT, SSCAL, SASUM
C                       FROM FORTRAN: ABS, AMAX1, MAX0, MIN0, SIGN

	INTEGER  LDA, N, ML, MU, IPVT(*)
	REAL     ABD(LDA,*), Z(*)
	REAL     RCOND

	REAL     SDOT, EK, T, WK, WKM
	REAL     ANORM, S, SASUM, SM, YNORM
	INTEGER  IS, INFO, J, JU, K, KB, KP1, L, LA, LM, LZ, M, MM


C                       ** COMPUTE 1-NORM OF A
	ANORM = 0.0E0
	L = ML + 1
	IS = L + MU
	DO 10 J = 1, N
	   ANORM = AMAX1(ANORM, SASUM(L,ABD(IS,J), 1))
	   IF (IS .GT. ML + 1) IS = IS - 1
	   IF (J .LE. MU) L = L + 1
	   IF (J .GE. N - ML) L = L - 1
   10 CONTINUE
C                                               ** FACTOR
	CALL SGBFA(ABD, LDA, N, ML, MU, IPVT, INFO)

C     RCOND = 1/(NORM(A)*(ESTIMATE OF NORM(INVERSE(A)))) .
C     ESTIMATE = NORM(Z)/NORM(Y) WHERE  A*Z = Y  AND  TRANS(A)*Y = E .
C     TRANS(A)  IS THE TRANSPOSE OF A .  THE COMPONENTS OF  E  ARE
C     CHOSEN TO CAUSE MAXIMUM LOCAL GROWTH IN THE ELEMENTS OF W  WHERE
C     TRANS(U)*W = E .  THE VECTORS ARE FREQUENTLY RESCALED TO AVOID
C     OVERFLOW.

C                     ** SOLVE TRANS(U)*W = E
	EK = 1.0E0
	DO 20 J = 1, N
	   Z(J) = 0.0E0
   20 CONTINUE

	M = ML + MU + 1
	JU = 0
	DO 100 K = 1, N
	   IF (Z(K) .NE. 0.0E0) EK = SIGN(EK, -Z(K))
	   IF (ABS(EK-Z(K)) .GT. ABS(ABD(M,K))) THEN
	      S = ABS(ABD(M,K))/ABS(EK-Z(K))
	      CALL SSCAL(N, S, Z, 1)
	      EK = S*EK
	   ENDIF
	   WK = EK - Z(K)
	   WKM = -EK - Z(K)
	   S = ABS(WK)
	   SM = ABS(WKM)
	   IF (ABD(M,K) .NE. 0.0E0) THEN
	      WK  = WK /ABD(M,K)
	      WKM = WKM/ABD(M,K)
	   ELSE
	      WK  = 1.0E0
	      WKM = 1.0E0
	   ENDIF
	   KP1 = K + 1
	   JU = MIN0(MAX0(JU, MU+IPVT(K)), N)
	   MM = M
	   IF (KP1 .LE. JU) THEN
	      DO 60 J = KP1, JU
	         MM = MM - 1
	         SM = SM + ABS(Z(J)+WKM*ABD(MM,J))
	         Z(J) = Z(J) + WK*ABD(MM,J)
	         S = S + ABS(Z(J))
   60       CONTINUE
	      IF (S .LT. SM) THEN
	         T = WKM - WK
	         WK = WKM
	         MM = M
	         DO 70 J = KP1, JU
	            MM = MM - 1
	            Z(J) = Z(J) + T*ABD(MM,J)
   70          CONTINUE
	      ENDIF
	   ENDIF
	   Z(K) = WK
  100 CONTINUE

	S = 1.0E0 / SASUM(N, Z, 1)
	CALL SSCAL(N, S, Z, 1)

C                         ** SOLVE TRANS(L)*Y = W
	DO 120 KB = 1, N
	   K = N + 1 - KB
	   LM = MIN0(ML, N-K)
	   IF (K .LT. N) Z(K) = Z(K) + SDOT(LM, ABD(M+1,K), 1, Z(K+1), 1)
	   IF (ABS(Z(K)) .GT. 1.0E0) THEN
	      S = 1.0E0 / ABS(Z(K))
	      CALL SSCAL(N, S, Z, 1)
	   ENDIF
	   L = IPVT(K)
	   T = Z(L)
	   Z(L) = Z(K)
	   Z(K) = T
  120 CONTINUE

	S = 1.0E0 / SASUM(N, Z, 1)
	CALL SSCAL(N, S, Z, 1)

	YNORM = 1.0E0
C                         ** SOLVE L*V = Y
	DO 140 K = 1, N
	   L = IPVT(K)
	   T = Z(L)
	   Z(L) = Z(K)
	   Z(K) = T
	   LM = MIN0(ML, N-K)
	   IF (K .LT. N) CALL SAXPY(LM, T, ABD(M+1,K), 1, Z(K+1), 1)
	   IF (ABS(Z(K)) .GT. 1.0E0) THEN
	      S = 1.0E0 / ABS(Z(K))
	      CALL SSCAL(N, S, Z, 1)
	      YNORM = S*YNORM
	   ENDIF
  140 CONTINUE

	S = 1.0E0/SASUM(N, Z, 1)
	CALL SSCAL(N, S, Z, 1)
	YNORM = S*YNORM
C                           ** SOLVE  U*Z = W
	DO 160 KB = 1, N
	   K = N + 1 - KB
	   IF (ABS(Z(K)) .GT. ABS(ABD(M,K))) THEN
	      S = ABS(ABD(M,K)) / ABS(Z(K))
	      CALL SSCAL(N, S, Z, 1)
	      YNORM = S*YNORM
	   ENDIF
	   IF (ABD(M,K) .NE. 0.0E0) Z(K) = Z(K)/ABD(M,K)
	   IF (ABD(M,K) .EQ. 0.0E0) Z(K) = 1.0E0
	   LM = MIN0(K, M) - 1
	   LA = M - LM
	   LZ = K - LM
	   T = -Z(K)
	   CALL SAXPY(LM, T, ABD(LA,K), 1, Z(LZ), 1)
  160 CONTINUE
C                              ** MAKE ZNORM = 1.0
	S = 1.0E0 / SASUM(N, Z, 1)
	CALL SSCAL(N, S, Z, 1)
	YNORM = S*YNORM

	IF (ANORM .NE. 0.0E0) RCOND = YNORM/ANORM
	IF (ANORM .EQ. 0.0E0) RCOND = 0.0E0
	RETURN
	END
	SUBROUTINE  SGBFA( ABD, LDA, N, ML, MU, IPVT, INFO )

C         FACTORS A REAL BAND MATRIX BY ELIMINATION.

C         REVISION DATE:  8/1/82
C         AUTHOR:  MOLER, C. B., (U. OF NEW MEXICO)

C     SGBFA IS USUALLY CALLED BY SBGCO, BUT IT CAN BE CALLED
C     DIRECTLY WITH A SAVING IN TIME IF  RCOND  IS NOT NEEDED.

C     INPUT:  SAME AS 'SGBCO'

C     ON RETURN:

C        ABD,IPVT    SAME AS 'SGBCO'

C        INFO    INTEGER
C                = 0  NORMAL VALUE.
C                = K  IF  U(K,K) .EQ. 0.0 .  THIS IS NOT AN ERROR
C                     CONDITION FOR THIS SUBROUTINE, BUT IT DOES
C                     INDICATE THAT SGBSL WILL DIVIDE BY ZERO IF
C                     CALLED.  USE  RCOND  IN SBGCO FOR A RELIABLE
C                     INDICATION OF SINGULARITY.

C     (SEE 'SGBCO' FOR DESCRIPTION OF BAND STORAGE MODE)

C     ROUTINES CALLED:  FROM BLAS:    SAXPY, SSCAL, ISAMAX
C                       FROM FORTRAN: MAX0, MIN0

	INTEGER  LDA, N, ML, MU, IPVT(*), INFO
	REAL     ABD(LDA,*)

	REAL     T
	INTEGER  I,ISAMAX,I0,J,JU,JZ,J0,J1,K,KP1,L,LM,M,MM,NM1


	M = ML + MU + 1
	INFO = 0
C                        ** ZERO INITIAL FILL-IN COLUMNS
	J0 = MU + 2
	J1 = MIN0(N, M) - 1
	DO 20 JZ = J0, J1
	   I0 = M + 1 - JZ
	   DO 10 I = I0, ML
	      ABD(I,JZ) = 0.0E0
   10    CONTINUE
   20 CONTINUE
	JZ = J1
	JU = 0

C                       ** GAUSSIAN ELIMINATION WITH PARTIAL PIVOTING
	NM1 = N - 1
	DO 120 K = 1, NM1
	   KP1 = K + 1
C                                  ** ZERO NEXT FILL-IN COLUMN
	   JZ = JZ + 1
	   IF (JZ .LE. N) THEN
	      DO 40 I = 1, ML
	         ABD(I,JZ) = 0.0E0
   40       CONTINUE
	   ENDIF
C                                  ** FIND L = PIVOT INDEX
	   LM = MIN0(ML, N-K)
	   L = ISAMAX(LM+1, ABD(M,K), 1) + M - 1
	   IPVT(K) = L + K - M

	   IF (ABD(L,K) .EQ. 0.0E0) THEN
C                                      ** ZERO PIVOT IMPLIES THIS COLUMN 
C                                      ** ALREADY TRIANGULARIZED
	      INFO = K
	   ELSE
C                                ** INTERCHANGE IF NECESSARY
	      IF (L .NE. M) THEN
	         T = ABD(L,K)
	         ABD(L,K) = ABD(M,K)
	         ABD(M,K) = T
	      ENDIF
C                                   ** COMPUTE MULTIPLIERS
	      T = -1.0E0 / ABD(M,K)
	      CALL SSCAL(LM, T, ABD(M+1,K), 1)

C                               ** ROW ELIMINATION WITH COLUMN INDEXING

	      JU = MIN0(MAX0(JU, MU+IPVT(K)), N)
	      MM = M
	      DO 80 J = KP1, JU
	         L = L - 1
	         MM = MM - 1
	         T = ABD(L,J)
	         IF (L .NE. MM) THEN
	            ABD(L,J) = ABD(MM,J)
	            ABD(MM,J) = T
	         ENDIF
	         CALL SAXPY(LM, T, ABD(M+1,K), 1, ABD(MM+1,J), 1)
   80       CONTINUE

	   ENDIF

  120 CONTINUE

	IPVT(N) = N
	IF (ABD(M,N) .EQ. 0.0E0) INFO = N
	RETURN
	END
	SUBROUTINE  SGBSL( ABD, LDA, N, ML, MU, IPVT, B, JOB )

C         SOLVES THE REAL BAND SYSTEM
C            A * X = B  OR  TRANSPOSE(A) * X = B
C         USING THE FACTORS COMPUTED BY SBGCO OR SGBFA.

C         REVISION DATE:  8/1/82
C         AUTHOR:  MOLER, C. B., (U. OF NEW MEXICO)

C     INPUT:

C        ABD     REAL(LDA, N)
C                THE OUTPUT FROM SBGCO OR SGBFA.

C        LDA     INTEGER
C                THE LEADING DIMENSION OF THE ARRAY  ABD .

C        N       INTEGER
C                THE ORDER OF THE ORIGINAL MATRIX.

C        ML      INTEGER
C                NUMBER OF DIAGONALS BELOW THE MAIN DIAGONAL.

C        MU      INTEGER
C                NUMBER OF DIAGONALS ABOVE THE MAIN DIAGONAL.

C        IPVT    INTEGER(N)
C                THE PIVOT VECTOR FROM SBGCO OR SGBFA.

C        B       REAL(N)
C                THE RIGHT HAND SIDE VECTOR.

C        JOB     INTEGER
C                = 0         TO SOLVE  A*X = B ,
C                = NONZERO   TO SOLVE  TRANS(A)*X = B , WHERE
C                            TRANS(A)  IS THE TRANSPOSE.

C     ON RETURN

C        B       THE SOLUTION VECTOR  X .

C     ERROR CONDITION

C        A DIVISION BY ZERO WILL OCCUR IF THE INPUT FACTOR CONTAINS A
C        ZERO ON THE DIAGONAL.  TECHNICALLY, THIS INDICATES SINGULARITY,
C        BUT IT IS OFTEN CAUSED BY IMPROPER ARGUMENTS OR IMPROPER
C        SETTING OF LDA .  IT WILL NOT OCCUR IF THE SUBROUTINES ARE
C        CALLED CORRECTLY AND IF SBGCO HAS SET RCOND .GT. 0.0
C        OR SGBFA HAS SET INFO .EQ. 0 .

C     TO COMPUTE  INVERSE(A) * C  WHERE  C  IS A MATRIX
C     WITH  P  COLUMNS
C           CALL SGBCO(ABD,LDA,N,ML,MU,IPVT,RCOND,Z)
C           IF (RCOND IS TOO SMALL) GO TO ...
C           DO 10 J = 1, P
C              CALL SGBSL(ABD,LDA,N,ML,MU,IPVT,C(1,J),0)
C        10 CONTINUE

C     ROUTINES CALLED:  FROM BLAS:    SAXPY, SDOT
C                       FROM FORTRAN: MIN0

	INTEGER  LDA, N, ML, MU, IPVT(*), JOB
	REAL     ABD(LDA,*), B(*)

	REAL     SDOT,T
	INTEGER  K,KB,L,LA,LB,LM,M,NM1


	M = MU + ML + 1
	NM1 = N - 1
	IF (JOB .EQ. 0) THEN
C                               ** JOB = 0 , SOLVE  A * X = B
C                               ** FIRST SOLVE L*Y = B
	   IF (ML .NE. 0) THEN
	      DO 20 K = 1, NM1
	         LM = MIN0(ML, N-K)
	         L = IPVT(K)
	         T = B(L)
	         IF (L .NE. K) THEN
	            B(L) = B(K)
	            B(K) = T
	         ENDIF
	         CALL SAXPY( LM, T, ABD(M+1,K), 1, B(K+1), 1 )
   20       CONTINUE
	   ENDIF
C                           ** NOW SOLVE  U*X = Y
	   DO 40 KB = 1, N
	      K = N + 1 - KB
	      B(K) = B(K) / ABD(M,K)
	      LM = MIN0(K, M) - 1
	      LA = M - LM
	      LB = K - LM
	      T = -B(K)
	      CALL SAXPY(LM, T, ABD(LA,K), 1, B(LB), 1)
   40    CONTINUE

	ELSE
C                          ** JOB = NONZERO, SOLVE  TRANS(A) * X = B
C                                  ** FIRST SOLVE  TRANS(U)*Y = B
	   DO 60 K = 1, N
	      LM = MIN0(K, M) - 1
	      LA = M - LM
	      LB = K - LM
	      T = SDOT(LM, ABD(LA,K), 1, B(LB), 1)
	      B(K) = (B(K) - T)/ABD(M,K)
   60    CONTINUE
C                                  ** NOW SOLVE TRANS(L)*X = Y
	   IF (ML .NE. 0) THEN
	      DO 80 KB = 1, NM1
	         K = N - KB
	         LM = MIN0(ML, N-K)
	         B(K) = B(K) + SDOT(LM, ABD(M+1,K), 1, B(K+1), 1)
	         L = IPVT(K)
	         IF (L .NE. K) THEN
	            T = B(L)
	            B(L) = B(K)
	            B(K) = T
	         ENDIF
   80       CONTINUE
	   ENDIF

	ENDIF

	RETURN
	END
	SUBROUTINE  SGECO( A, LDA, N,IPVT, RCOND, Z )

C         FACTORS A REAL MATRIX BY GAUSSIAN ELIMINATION
C         AND ESTIMATES THE CONDITION OF THE MATRIX.

C         REVISION DATE:  8/1/82
C         AUTHOR:  MOLER, C. B., (U. OF NEW MEXICO)

C         IF  RCOND  IS NOT NEEDED, SGEFA IS SLIGHTLY FASTER.
C         TO SOLVE  A*X = B , FOLLOW SGECO BY SGESL.

C     ON ENTRY

C        A       REAL(LDA, N)
C                THE MATRIX TO BE FACTORED.

C        LDA     INTEGER
C                THE LEADING DIMENSION OF THE ARRAY  A .

C        N       INTEGER
C                THE ORDER OF THE MATRIX  A .

C     ON RETURN

C        A       AN UPPER TRIANGULAR MATRIX AND THE MULTIPLIERS
C                WHICH WERE USED TO OBTAIN IT.
C                THE FACTORIZATION CAN BE WRITTEN  A = L*U , WHERE
C                L  IS A PRODUCT OF PERMUTATION AND UNIT LOWER
C                TRIANGULAR MATRICES AND  U  IS UPPER TRIANGULAR.

C        IPVT    INTEGER(N)
C                AN INTEGER VECTOR OF PIVOT INDICES.

C        RCOND   REAL
C                AN ESTIMATE OF THE RECIPROCAL CONDITION OF  A .
C                FOR THE SYSTEM  A*X = B , RELATIVE PERTURBATIONS
C                IN  A  AND  B  OF SIZE  EPSILON  MAY CAUSE
C                RELATIVE PERTURBATIONS IN  X  OF SIZE  EPSILON/RCOND .
C                IF  RCOND  IS SO SMALL THAT THE LOGICAL EXPRESSION
C                           1.0 + RCOND .EQ. 1.0
C                IS TRUE, THEN  A  MAY BE SINGULAR TO WORKING
C                PRECISION.  IN PARTICULAR,  RCOND  IS ZERO  IF
C                EXACT SINGULARITY IS DETECTED OR THE ESTIMATE
C                UNDERFLOWS.

C        Z       REAL(N)
C                A WORK VECTOR WHOSE CONTENTS ARE USUALLY UNIMPORTANT.
C                IF  A  IS CLOSE TO A SINGULAR MATRIX, THEN  Z  IS
C                AN APPROXIMATE NULL VECTOR IN THE SENSE THAT
C                NORM(A*Z) = RCOND*NORM(A)*NORM(Z) .

C     ROUTINES CALLED:  FROM LINPACK: SGEFA
C                       FROM BLAS:    SAXPY, SDOT, SSCAL, SASUM
C                       FROM FORTRAN: ABS, AMAX1, SIGN

	INTEGER  LDA, N, IPVT(*)
	REAL     A(LDA,*), Z(*)
	REAL     RCOND

	REAL     SDOT,EK,T,WK,WKM
	REAL     ANORM,S,SASUM,SM,YNORM
	INTEGER  INFO,J,K,KB,KP1,L


C                        ** COMPUTE 1-NORM OF A
	ANORM = 0.0E0
	DO 10 J = 1, N
	   ANORM = AMAX1( ANORM, SASUM(N,A(1,J),1) )
   10 CONTINUE
C                                      ** FACTOR
	CALL SGEFA(A,LDA,N,IPVT,INFO)

C     RCOND = 1/(NORM(A)*(ESTIMATE OF NORM(INVERSE(A)))) .
C     ESTIMATE = NORM(Z)/NORM(Y) WHERE  A*Z = Y  AND  TRANS(A)*Y = E .
C     TRANS(A)  IS THE TRANSPOSE OF A .  THE COMPONENTS OF  E  ARE
C     CHOSEN TO CAUSE MAXIMUM LOCAL GROWTH IN THE ELEMENTS OF W  WHERE
C     TRANS(U)*W = E .  THE VECTORS ARE FREQUENTLY RESCALED TO AVOID
C     OVERFLOW.

C                        ** SOLVE TRANS(U)*W = E
	EK = 1.0E0
	DO 20 J = 1, N
	   Z(J) = 0.0E0
   20 CONTINUE

	DO 100 K = 1, N
	   IF (Z(K) .NE. 0.0E0) EK = SIGN(EK, -Z(K))
	   IF (ABS(EK-Z(K)) .GT. ABS(A(K,K))) THEN
	      S = ABS(A(K,K)) / ABS(EK-Z(K))
	      CALL SSCAL(N, S, Z, 1)
	      EK = S*EK
	   ENDIF
	   WK = EK - Z(K)
	   WKM = -EK - Z(K)
	   S = ABS(WK)
	   SM = ABS(WKM)
	   IF (A(K,K) .NE. 0.0E0) THEN
	      WK  = WK  / A(K,K)
	      WKM = WKM / A(K,K)
	   ELSE
	      WK  = 1.0E0
	      WKM = 1.0E0
	   ENDIF
	   KP1 = K + 1
	   IF (KP1 .LE. N) THEN
	      DO 60 J = KP1, N
	         SM = SM + ABS(Z(J)+WKM*A(K,J))
	         Z(J) = Z(J) + WK*A(K,J)
	         S = S + ABS(Z(J))
   60       CONTINUE
	      IF (S .LT. SM) THEN
	         T = WKM - WK
	         WK = WKM
	         DO 70 J = KP1, N
	            Z(J) = Z(J) + T*A(K,J)
   70          CONTINUE
	      ENDIF
	   ENDIF
	   Z(K) = WK
  100 CONTINUE

	S = 1.0E0 / SASUM(N, Z, 1)
	CALL SSCAL(N, S, Z, 1)
C                                ** SOLVE TRANS(L)*Y = W
	DO 120 KB = 1, N
	   K = N + 1 - KB
	   IF (K .LT. N) Z(K) = Z(K) + SDOT(N-K, A(K+1,K), 1, Z(K+1), 1)
	   IF (ABS(Z(K)) .GT. 1.0E0) THEN
	      S = 1.0E0/ABS(Z(K))
	      CALL SSCAL(N, S, Z, 1)
	   ENDIF
	   L = IPVT(K)
	   T = Z(L)
	   Z(L) = Z(K)
	   Z(K) = T
  120 CONTINUE

	S = 1.0E0 / SASUM(N, Z, 1)
	CALL SSCAL(N, S, Z, 1)
C                                 ** SOLVE L*V = Y
	YNORM = 1.0E0
	DO 140 K = 1, N
	   L = IPVT(K)
	   T = Z(L)
	   Z(L) = Z(K)
	   Z(K) = T
	   IF (K .LT. N) CALL SAXPY(N-K, T, A(K+1,K), 1, Z(K+1), 1)
	   IF (ABS(Z(K)) .GT. 1.0E0) THEN
	      S = 1.0E0/ABS(Z(K))
	      CALL SSCAL(N, S, Z, 1)
	      YNORM = S*YNORM
	   ENDIF
  140 CONTINUE

	S = 1.0E0 / SASUM(N, Z, 1)
	CALL SSCAL(N, S, Z, 1)
C                                  ** SOLVE  U*Z = V
	YNORM = S*YNORM
	DO 160 KB = 1, N
	   K = N + 1 - KB
	   IF (ABS(Z(K)) .GT. ABS(A(K,K))) THEN
	      S = ABS(A(K,K))/ABS(Z(K))
	      CALL SSCAL(N, S, Z, 1)
	      YNORM = S*YNORM
	   ENDIF
	   IF (A(K,K) .NE. 0.0E0) Z(K) = Z(K)/A(K,K)
	   IF (A(K,K) .EQ. 0.0E0) Z(K) = 1.0E0
	   T = -Z(K)
	   CALL SAXPY(K-1, T, A(1,K), 1, Z(1), 1)
  160 CONTINUE
C                                   ** MAKE ZNORM = 1.0
	S = 1.0E0 / SASUM(N, Z, 1)
	CALL SSCAL(N, S, Z, 1)
	YNORM = S*YNORM

	IF (ANORM .NE. 0.0E0) RCOND = YNORM/ANORM
	IF (ANORM .EQ. 0.0E0) RCOND = 0.0E0
	RETURN
	END
	SUBROUTINE  SGEFA( A, LDA, N, IPVT, INFO )

C         FACTORS A REAL MATRIX BY GAUSSIAN ELIMINATION.

C         REVISION DATE:  8/1/82
C         AUTHOR:  MOLER, C. B., (U. OF NEW MEXICO)

C     SGEFA IS USUALLY CALLED BY SGECO, BUT IT CAN BE CALLED
C     DIRECTLY WITH A SAVING IN TIME IF  RCOND  IS NOT NEEDED.
C     (TIME FOR SGECO) = (1 + 9/N)*(TIME FOR SGEFA) .

C     INPUT:  SAME AS 'SGECO'

C     ON RETURN:

C        A,IPVT  SAME AS 'SGECO'

C        INFO    INTEGER
C                = 0  NORMAL VALUE.
C                = K  IF  U(K,K) .EQ. 0.0 .  THIS IS NOT AN ERROR
C                     CONDITION FOR THIS SUBROUTINE, BUT IT DOES
C                     INDICATE THAT SGESL OR SGEDI WILL DIVIDE BY ZERO
C                     IF CALLED.  USE  RCOND  IN SGECO FOR A RELIABLE
C                     INDICATION OF SINGULARITY.

C     ROUTINES CALLED:  FROM BLAS:    SAXPY, SSCAL, ISAMAX

	INTEGER  LDA, N, IPVT(*), INFO
	REAL     A(LDA,*)

	REAL     T
	INTEGER  ISAMAX,J,K,KP1,L,NM1


C                      ** GAUSSIAN ELIMINATION WITH PARTIAL PIVOTING
	INFO = 0
	NM1 = N - 1
	DO 60 K = 1, NM1
	   KP1 = K + 1
C                                            ** FIND L = PIVOT INDEX
	   L = ISAMAX( N-K+1, A(K,K), 1) + K-1
	   IPVT(K) = L

	   IF (A(L,K) .EQ. 0.0E0) THEN
C                                     ** ZERO PIVOT IMPLIES THIS COLUMN 
C                                     ** ALREADY TRIANGULARIZED
	      INFO = K
	   ELSE
C                                     ** INTERCHANGE IF NECESSARY
	      IF (L .NE. K) THEN
	         T = A(L,K)
	         A(L,K) = A(K,K)
	         A(K,K) = T
	      ENDIF
C                                     ** COMPUTE MULTIPLIERS
	      T = -1.0E0 / A(K,K)
	      CALL SSCAL( N-K, T, A(K+1,K), 1 )

C                              ** ROW ELIMINATION WITH COLUMN INDEXING
	      DO 30 J = KP1, N
	         T = A(L,J)
	         IF (L .NE. K) THEN
	            A(L,J) = A(K,J)
	            A(K,J) = T
	         ENDIF
	         CALL SAXPY( N-K, T, A(K+1,K), 1, A(K+1,J), 1 )
   30       CONTINUE

	   ENDIF

   60 CONTINUE

	IPVT(N) = N
	IF (A(N,N) .EQ. 0.0E0) INFO = N
	RETURN
	END
	SUBROUTINE  SGESL( A, LDA, N,IPVT, B, JOB )

C         SOLVES THE REAL SYSTEM
C            A * X = B  OR  TRANS(A) * X = B
C         USING THE FACTORS COMPUTED BY SGECO OR SGEFA.

C         REVISION DATE:  8/1/82
C         AUTHOR:  MOLER, C. B., (U. OF NEW MEXICO)

C     ON ENTRY

C        A       REAL(LDA, N)
C                THE OUTPUT FROM SGECO OR SGEFA.

C        LDA     INTEGER
C                THE LEADING DIMENSION OF THE ARRAY  A .

C        N       INTEGER
C                THE ORDER OF THE MATRIX  A .

C        IPVT    INTEGER(N)
C                THE PIVOT VECTOR FROM SGECO OR SGEFA.

C        B       REAL(N)
C                THE RIGHT HAND SIDE VECTOR.

C        JOB     INTEGER
C                = 0         TO SOLVE  A*X = B ,
C                = NONZERO   TO SOLVE  TRANS(A)*X = B  WHERE
C                            TRANS(A)  IS THE TRANSPOSE.

C     ON RETURN

C        B       THE SOLUTION VECTOR  X .

C     ERROR CONDITION

C        A DIVISION BY ZERO WILL OCCUR IF THE INPUT FACTOR CONTAINS A
C        ZERO ON THE DIAGONAL.  TECHNICALLY, THIS INDICATES SINGULARITY,
C        BUT IT IS OFTEN CAUSED BY IMPROPER ARGUMENTS OR IMPROPER
C        SETTING OF LDA .  IT WILL NOT OCCUR IF THE SUBROUTINES ARE
C        CALLED CORRECTLY AND IF SGECO HAS SET RCOND .GT. 0.0
C        OR SGEFA HAS SET INFO .EQ. 0 .

C     TO COMPUTE  INVERSE(A) * C  WHERE  C  IS A MATRIX
C     WITH  P  COLUMNS
C           CALL SGECO(A,LDA,N,IPVT,RCOND,Z)
C           IF (RCOND IS TOO SMALL) GO TO ...
C           DO 10 J = 1, P
C              CALL SGESL(A,LDA,N,IPVT,C(1,J),0)
C        10 CONTINUE


C     ROUTINES CALLED:  FROM BLAS:    SAXPY, SDOT

	INTEGER  LDA, N, IPVT(*), JOB
	REAL     A(LDA,*), B(*)

	REAL     SDOT,T
	INTEGER  K,KB,L,NM1


	NM1 = N - 1
	IF (JOB .EQ. 0) THEN
C                                 ** JOB = 0 , SOLVE  A * X = B
C                                     ** FIRST SOLVE  L*Y = B
	   DO 20 K = 1, NM1
	      L = IPVT(K)
	      T = B(L)
	      IF (L .NE. K) THEN
	         B(L) = B(K)
	         B(K) = T
	      ENDIF
	      CALL SAXPY( N-K, T, A(K+1,K), 1, B(K+1), 1 )
   20    CONTINUE
C                                    ** NOW SOLVE  U*X = Y
	   DO 40 KB = 1, N
	      K = N + 1 - KB
	      B(K) = B(K) / A(K,K)
	      T = -B(K)
	      CALL SAXPY( K-1, T, A(1,K), 1, B(1), 1 )
   40    CONTINUE

	ELSE
C                         ** JOB = NONZERO, SOLVE  TRANS(A) * X = B
C                                    ** FIRST SOLVE  TRANS(U)*Y = B
	   DO 60 K = 1, N
	      T = SDOT( K-1, A(1,K), 1, B(1), 1 )
	      B(K) = (B(K) - T) / A(K,K)
   60    CONTINUE
C                                    ** NOW SOLVE  TRANS(L)*X = Y
	   DO 80 KB = 1, NM1
	      K = N - KB
	      B(K) = B(K) + SDOT( N-K, A(K+1,K), 1, B(K+1), 1 )
	      L = IPVT(K)
	      IF (L .NE. K) THEN
	         T = B(L)
	         B(L) = B(K)
	         B(K) = T
	      ENDIF
   80    CONTINUE

	ENDIF

	RETURN
	END
	REAL FUNCTION  SASUM( N, SX, INCX )

C  --INPUT--  N  NUMBER OF ELEMENTS IN VECTOR TO BE SUMMED
C            SX  SING-PREC ARRAY, LENGTH 1+(N-1)*INCX, CONTAINING VECTOR
C          INCX  SPACING OF VECTOR ELEMENTS IN 'SX'

C --OUTPUT-- SASUM   SUM FROM 0 TO N-1 OF  ABS(SX(1+I*INCX))

	REAL SX(*)


	SASUM = 0.0
	IF( N.LE.0 )  RETURN
	IF( INCX.NE.1 ) THEN
C                                          ** NON-UNIT INCREMENTS
	    DO 10 I = 1, 1+(N-1)*INCX, INCX
	       SASUM = SASUM + ABS(SX(I))
   10     CONTINUE
	ELSE
C                                          ** UNIT INCREMENTS
	   M = MOD(N,6)
	   IF( M.NE.0 ) THEN
C                             ** CLEAN-UP LOOP SO REMAINING VECTOR 
C                             ** LENGTH IS A MULTIPLE OF 6.
	      DO 30  I = 1, M
	        SASUM = SASUM + ABS(SX(I))
   30       CONTINUE
	   ENDIF
C                              ** UNROLL LOOP FOR SPEED
	   DO 50  I = M+1, N, 6
	     SASUM = SASUM + ABS(SX(I))   + ABS(SX(I+1)) + ABS(SX(I+2))
     $                   + ABS(SX(I+3)) + ABS(SX(I+4)) + ABS(SX(I+5))
   50    CONTINUE
	ENDIF

	RETURN
	END
	SUBROUTINE     SAXPY( N, SA, SX, INCX, SY, INCY )

C          Y = A*X + Y  (X, Y = VECTORS, A = SCALAR)

C  --INPUT--
C        N  NUMBER OF ELEMENTS IN INPUT VECTORS 'X' AND 'Y'
C       SA  SINGLE PRECISION SCALAR MULTIPLIER 'A'
C       SX  SING-PREC ARRAY CONTAINING VECTOR 'X'
C     INCX  SPACING OF ELEMENTS OF VECTOR 'X' IN 'SX'
C       SY  SING-PREC ARRAY CONTAINING VECTOR 'Y'
C     INCY  SPACING OF ELEMENTS OF VECTOR 'Y' IN 'SY'

C --OUTPUT--
C       SY   FOR I = 0 TO N-1, OVERWRITE  SY(LY+I*INCY) WITH 
C                 SA*SX(LX+I*INCX) + SY(LY+I*INCY), 
C            WHERE LX = 1          IF INCX .GE. 0,
C                     = (-INCX)*N  IF INCX .LT. 0
C            AND LY IS DEFINED IN A SIMILAR WAY USING INCY.

	REAL SX(*), SY(*), SA


	IF( N.LE.0 .OR. SA.EQ.0.0 ) RETURN

	IF ( INCX.EQ.INCY .AND. INCX.GT.1 )  THEN

	    DO 10  I = 1, 1+(N-1)*INCX, INCX
	       SY(I) = SY(I) + SA * SX(I)
   10     CONTINUE

	ELSE IF ( INCX.EQ.INCY .AND. INCX.EQ.1 )  THEN

C                                        ** EQUAL, UNIT INCREMENTS
	   M = MOD(N,4)
	   IF( M .NE. 0 ) THEN
C                            ** CLEAN-UP LOOP SO REMAINING VECTOR LENGTH
C                            ** IS A MULTIPLE OF 4.
	      DO 20  I = 1, M
	        SY(I) = SY(I) + SA * SX(I)
   20       CONTINUE
	   ENDIF
C                              ** UNROLL LOOP FOR SPEED
	   DO 30  I = M+1, N, 4
	      SY(I)   = SY(I)   + SA * SX(I)
	      SY(I+1) = SY(I+1) + SA * SX(I+1)
	      SY(I+2) = SY(I+2) + SA * SX(I+2)
	      SY(I+3) = SY(I+3) + SA * SX(I+3)
   30    CONTINUE

	ELSE
C               ** NONEQUAL OR NONPOSITIVE INCREMENTS.
	   IX = 1
	   IY = 1
	   IF( INCX.LT.0 )  IX = 1 + (N-1)*(-INCX)
	   IF( INCY.LT.0 )  IY = 1 + (N-1)*(-INCY)
	   DO 40  I = 1, N
	      SY(IY) = SY(IY) + SA*SX(IX)
	      IX = IX + INCX
	      IY = IY + INCY
   40    CONTINUE

	ENDIF

	RETURN
	END
	REAL FUNCTION  SDOT( N, SX, INCX, SY, INCY )

C          S.P. DOT PRODUCT OF VECTORS  'X'  AND  'Y'

C  --INPUT--
C        N  NUMBER OF ELEMENTS IN INPUT VECTORS 'X' AND 'Y'
C       SX  SING-PREC ARRAY CONTAINING VECTOR 'X'
C     INCX  SPACING OF ELEMENTS OF VECTOR 'X' IN 'SX'
C       SY  SING-PREC ARRAY CONTAINING VECTOR 'Y'
C     INCY  SPACING OF ELEMENTS OF VECTOR 'Y' IN 'SY'

C --OUTPUT--
C     SDOT   SUM FOR I = 0 TO N-1 OF  SX(LX+I*INCX) * SY(LY+I*INCY),
C            WHERE  LX = 1          IF INCX .GE. 0, 
C                      = (-INCX)*N  IF INCX .LT. 0,
C            AND LY IS DEFINED IN A SIMILAR WAY USING INCY.

	REAL SX(*), SY(*)


	SDOT = 0.0
	IF( N.LE.0 )  RETURN

	IF ( INCX.EQ.INCY .AND. INCX.GT.1 )  THEN

	    DO 10  I = 1, 1+(N-1)*INCX, INCX
	       SDOT = SDOT + SX(I) * SY(I)
   10     CONTINUE

	ELSE IF ( INCX.EQ.INCY .AND. INCX.EQ.1 )  THEN

C                                        ** EQUAL, UNIT INCREMENTS
	   M = MOD(N,5)
	   IF( M .NE. 0 ) THEN
C                            ** CLEAN-UP LOOP SO REMAINING VECTOR LENGTH
C                            ** IS A MULTIPLE OF 4.
	      DO 20  I = 1, M
	         SDOT = SDOT + SX(I) * SY(I)
   20       CONTINUE
	   ENDIF
C                              ** UNROLL LOOP FOR SPEED
	   DO 30  I = M+1, N, 5
	      SDOT = SDOT + SX(I)*SY(I)     + SX(I+1)*SY(I+1)
     $                  + SX(I+2)*SY(I+2) + SX(I+3)*SY(I+3)
     $                  + SX(I+4)*SY(I+4)
   30    CONTINUE

	ELSE
C               ** NONEQUAL OR NONPOSITIVE INCREMENTS.
	   IX = 1
	   IY = 1
	   IF( INCX.LT.0 )  IX = 1 + (N-1)*(-INCX)
	   IF( INCY.LT.0 )  IY = 1 + (N-1)*(-INCY)
	   DO 40  I = 1, N
	      SDOT = SDOT + SX(IX) * SY(IY)
	      IX = IX + INCX
	      IY = IY + INCY
   40    CONTINUE

	ENDIF

	RETURN
	END
	SUBROUTINE     SSCAL( N, SA, SX, INCX )

C         CALCULATE  X = A*X  (X = VECTOR, A = SCALAR)

C  --INPUT--  N  NUMBER OF ELEMENTS IN VECTOR
C            SA  SINGLE PRECISION SCALE FACTOR
C            SX  SING-PREC ARRAY, LENGTH 1+(N-1)*INCX, CONTAINING VECTOR
C          INCX  SPACING OF VECTOR ELEMENTS IN 'SX'

C --OUTPUT-- SX  REPLACE  SX(1+I*INCX)  WITH  SA * SX(1+I*INCX) 
C                FOR I = 0 TO N-1

	REAL SA, SX(*)


	IF( N.LE.0 ) RETURN

	IF( INCX.NE.1 ) THEN

	    DO 10  I = 1, 1+(N-1)*INCX, INCX
	       SX(I) = SA * SX(I)
   10     CONTINUE

	ELSE

	   M = MOD(N,5)
	   IF( M.NE.0 ) THEN
C                           ** CLEAN-UP LOOP SO REMAINING VECTOR LENGTH
C                           ** IS A MULTIPLE OF 5.
	      DO 30  I = 1, M
	         SX(I) = SA * SX(I)
   30       CONTINUE
	   ENDIF
C                             ** UNROLL LOOP FOR SPEED
	   DO 50  I = M+1, N, 5
	      SX(I)   = SA * SX(I)
	      SX(I+1) = SA * SX(I+1)
	      SX(I+2) = SA * SX(I+2)
	      SX(I+3) = SA * SX(I+3)
	      SX(I+4) = SA * SX(I+4)
   50    CONTINUE

	ENDIF

	RETURN
	END
	SUBROUTINE     SSWAP( N, SX, INCX, SY, INCY )

C          INTERCHANGE S.P VECTORS  X  AND  Y

C  --INPUT--
C        N  NUMBER OF ELEMENTS IN INPUT VECTORS 'X' AND 'Y'
C       SX  SING-PREC ARRAY CONTAINING VECTOR 'X'
C     INCX  SPACING OF ELEMENTS OF VECTOR 'X' IN 'SX'
C       SY  SING-PREC ARRAY CONTAINING VECTOR 'Y'
C     INCY  SPACING OF ELEMENTS OF VECTOR 'Y' IN 'SY'

C --OUTPUT--
C       SX  INPUT VECTOR SY (UNCHANGED IF N .LE. 0)
C       SY  INPUT VECTOR SX (UNCHANGED IF N .LE. 0)

C     FOR I = 0 TO N-1, INTERCHANGE  SX(LX+I*INCX) AND SY(LY+I*INCY),
C     WHERE LX = 1          IF INCX .GE. 0, 
C              = (-INCX)*N  IF INCX .LT. 0
C     AND LY IS DEFINED IN A SIMILAR WAY USING INCY.

	REAL SX(*), SY(*), STEMP1, STEMP2, STEMP3


	IF( N.LE.0 ) RETURN

	IF ( INCX.EQ.INCY .AND. INCX.GT.1 )  THEN

	    DO 10  I = 1, 1+(N-1)*INCX, INCX
	       STEMP1 = SX(I)
	       SX(I) = SY(I)
	       SY(I) = STEMP1
   10     CONTINUE

	ELSE IF ( INCX.EQ.INCY .AND. INCX.EQ.1 )  THEN

C                                        ** EQUAL, UNIT INCREMENTS
	   M = MOD(N,3)
	   IF( M .NE. 0 ) THEN
C                            ** CLEAN-UP LOOP SO REMAINING VECTOR LENGTH
C                            ** IS A MULTIPLE OF 3.
	      DO 20  I = 1, M
	         STEMP1 = SX(I)
	         SX(I) = SY(I)
	         SY(I) = STEMP1
   20       CONTINUE
	   ENDIF
C                              ** UNROLL LOOP FOR SPEED
	   DO 30  I = M+1, N, 3
	      STEMP1  = SX(I)
	      STEMP2  = SX(I+1)
	      STEMP3  = SX(I+2)
	      SX(I)   = SY(I)
	      SX(I+1) = SY(I+1)
	      SX(I+2) = SY(I+2)
	      SY(I)   = STEMP1
	      SY(I+1) = STEMP2
	      SY(I+2) = STEMP3
   30    CONTINUE

	ELSE
C               ** NONEQUAL OR NONPOSITIVE INCREMENTS.
	   IX = 1
	   IY = 1
	   IF( INCX.LT.0 )  IX = 1 + (N-1)*(-INCX)
	   IF( INCY.LT.0 )  IY = 1 + (N-1)*(-INCY)
	   DO 40  I = 1, N
	      STEMP1 = SX(IX)
	      SX(IX) = SY(IY)
	      SY(IY) = STEMP1
	      IX = IX + INCX
	      IY = IY + INCY
   40    CONTINUE

	ENDIF

	RETURN
	END
	INTEGER FUNCTION  ISAMAX( N, SX, INCX )

C  --INPUT--  N  NUMBER OF ELEMENTS IN VECTOR OF INTEREST
C            SX  SING-PREC ARRAY, LENGTH 1+(N-1)*INCX, CONTAINING VECTOR
C          INCX  SPACING OF VECTOR ELEMENTS IN 'SX'

C --OUTPUT-- ISAMAX   FIRST I, I = 1 TO N, TO MAXIMIZE
C                         ABS(SX(1+(I-1)*INCX))

	REAL SX(*), SMAX, XMAG


	IF( N.LE.0 ) THEN
	   ISAMAX = 0
	ELSE IF( N.EQ.1 ) THEN
	   ISAMAX = 1
	ELSE
	   SMAX = 0.0
	   II = 1
	   DO 20  I = 1, 1+(N-1)*INCX, INCX
	      XMAG = ABS(SX(I))
	      IF( SMAX.LT.XMAG ) THEN
	         SMAX = XMAG
	         ISAMAX = II
	      ENDIF
	      II = II + 1
   20    CONTINUE
	ENDIF

	RETURN
	END
      FUNCTION D1MACH(i)

*-----------------------------------------------------------------------------*
*= PURPOSE:                                                                  =*
*= D1MACH calculates various machine constants in single precision.          =*
*-----------------------------------------------------------------------------*
*= PARAMETERS:                                                               =*
*=   I       -  INTEGER, identifies the machine constant (0<I<5)         (I) =*
*=   D1MACH  -  REAL, machine constant in single precision               (O) =*
*=      I=1     - the smallest non-vanishing normalized floating-point       =*
*=                power of the radix, i.e., D1MACH=FLOAT(IBETA)**MINEXP      =*
*=      I=2     - the largest finite floating-point number.  In              =*
*=                particular D1MACH=(1.0-EPSNEG)*FLOAT(IBETA)**MAXEXP        =*
*=                Note - on some machines D1MACH will be only the            =*
*=                second, or perhaps third, largest number, being            =*
*=                too small by 1 or 2 units in the last digit of             =*
*=                the significand.                                           =*
*=      I=3     - A small positive floating-point number such that           =*
*=                1.0-D1MACH .NE. 1.0. In particular, if IBETA = 2           =*
*=                or  IRND = 0, D1MACH = FLOAT(IBETA)**NEGEPS.               =*
*=                Otherwise,  D1MACH = (IBETA**NEGEPS)/2.  Because           =*
*=                NEGEPS is bounded below by -(IT+3), D1MACH may not         =*
*=                be the smallest number that can alter 1.0 by               =*
*=                subtraction.                                               =*
*=      I=4     - the smallest positive floating-point number such           =*
*=                that  1.0+D1MACH .NE. 1.0. In particular, if either        =*
*=                IBETA = 2  or  IRND = 0, D1MACH=FLOAT(IBETA)**MACHEP.      =*
*=                Otherwise, D1MACH=(FLOAT(IBETA)**MACHEP)/2                 =*
*=  (see routine T665D for more information on different constants)          =*
*-----------------------------------------------------------------------------*

      DOUBLE PRECISION d1mach
      INTEGER i
   
      LOGICAL doinit
      DATA doinit/.TRUE./
      SAVE doinit

      DOUBLE PRECISION dmach(4) 
      SAVE dmach

      IF (( i .GE. 1 ) .AND. ( i .LE. 4 )) THEN
* compute constants at first call only
        IF (doinit) THEN
           CALL t665d(dmach)
           doinit = .FALSE.
        ENDIF
        d1mach = dmach(i)
      ELSE
        WRITE(0,*) '>>> ERROR (D1MACH) <<<  invalid argument'
        STOP
      ENDIF

      END


C      ALGORITHM 665, COLLECTED ALGORITHMS FROM ACM.
C      THIS WORK PUBLISHED IN TRANSACTIONS ON MATHEMATICAL SOFTWARE,
C      VOL. 14, NO. 4, PP. 303-311.
      SUBROUTINE T665D(DMACH)
C-----------------------------------------------------------------------
C This subroutine is a double precision version of subroutine T665R.
C See code of T665R for detailed comments and explanation
C-----------------------------------------------------------------------
      DOUBLE PRECISION DMACH(4)
      INTEGER I,IBETA,IEXP,IRND,IT,ITEMP,IZ,J,K,MACHEP,MAXEXP,
     1        MINEXP,MX,NEGEP,NGRD,NXRES
CS    REAL A,B,BETA,BETAIN,BETAH,CONV,EPS,EPSNEG,ONE,T,TEMP,TEMPA,
CS   1     TEMP1,TWO,XMAX,XMIN,Y,Z,ZERO
      DOUBLE PRECISION A,B,BETA,BETAIN,BETAH,CONV,EPS,EPSNEG,ONE,
     1                 T,TEMP,TEMPA,TEMP1,TWO,XMAX,XMIN,Y,Z,ZERO
C-----------------------------------------------------------------------
CS    CONV(I) = REAL(I)
      CONV(I) = DBLE(I)
      ONE = CONV(1)
      TWO = ONE + ONE
      ZERO = ONE - ONE
C-----------------------------------------------------------------------
C  Determine IBETA, BETA ala Malcolm.
C-----------------------------------------------------------------------
      A = ONE
   10 A = A + A
         TEMP = A+ONE
         TEMP1 = TEMP-A
         IF (TEMP1-ONE .EQ. ZERO) GO TO 10
      B = ONE
   20 B = B + B
         TEMP = A+B
         ITEMP = INT(TEMP-A)
         IF (ITEMP .EQ. 0) GO TO 20
      IBETA = ITEMP
      BETA = CONV(IBETA)
C-----------------------------------------------------------------------
C  Determine IT, IRND.
C-----------------------------------------------------------------------
      IT = 0
      B = ONE
  100 IT = IT + 1
         B = B * BETA
         TEMP = B+ONE
         TEMP1 = TEMP-B
         IF (TEMP1-ONE .EQ. ZERO) GO TO 100
      IRND = 0
      BETAH = BETA / TWO
      TEMP = A+BETAH
      IF (TEMP-A .NE. ZERO) IRND = 1
      TEMPA = A + BETA
      TEMP = TEMPA+BETAH
      IF ((IRND .EQ. 0) .AND. (TEMP-TEMPA .NE. ZERO)) IRND = 2
C-----------------------------------------------------------------------
C  Determine NEGEP, EPSNEG.
C-----------------------------------------------------------------------
      NEGEP = IT + 3
      BETAIN = ONE / BETA
      A = ONE
      DO 200 I = 1, NEGEP
         A = A * BETAIN
  200 CONTINUE
      B = A
  210 TEMP = ONE-A
         IF (TEMP-ONE .NE. ZERO) GO TO 220
         A = A * BETA
         NEGEP = NEGEP - 1
      GO TO 210
  220 NEGEP = -NEGEP
      EPSNEG = A
      IF ((IBETA .EQ. 2) .OR. (IRND .EQ. 0)) GO TO 300
      A = (A*(ONE+A)) / TWO
      TEMP = ONE-A
      IF (TEMP-ONE .NE. ZERO) EPSNEG = A
C-----------------------------------------------------------------------
C  Determine MACHEP, EPS.
C-----------------------------------------------------------------------
  300 MACHEP = -IT - 3
      A = B
  310 TEMP = ONE+A
         IF (TEMP-ONE .NE. ZERO) GO TO 320
         A = A * BETA
         MACHEP = MACHEP + 1
      GO TO 310
  320 EPS = A
      TEMP = TEMPA+BETA*(ONE+EPS)
      IF ((IBETA .EQ. 2) .OR. (IRND .EQ. 0)) GO TO 350
      A = (A*(ONE+A)) / TWO
      TEMP = ONE+A
      IF (TEMP-ONE .NE. ZERO) EPS = A
C-----------------------------------------------------------------------
C  Determine NGRD.
C-----------------------------------------------------------------------
  350 NGRD = 0
      TEMP = ONE+EPS
      IF ((IRND .EQ. 0) .AND. (TEMP*ONE-ONE .NE. ZERO)) NGRD = 1
C-----------------------------------------------------------------------
C  Determine IEXP, MINEXP, XMIN.
C
C  Loop to determine largest I and K = 2**I such that
C         (1/BETA) ** (2**(I))
C  does not underflow.
C  Exit from loop is signaled by an underflow.
C-----------------------------------------------------------------------
      I = 0
      K = 1
      Z = BETAIN
      T = ONE + EPS
      NXRES = 0
  400 Y = Z
         Z = Y * Y
C-----------------------------------------------------------------------
C  Check for underflow here.
C-----------------------------------------------------------------------
         A = Z * ONE
         TEMP = Z * T
         IF ((A+A .EQ. ZERO) .OR. (ABS(Z) .GE. Y)) GO TO 410
         TEMP1 = TEMP * BETAIN
         IF (TEMP1*BETA .EQ. Z) GO TO 410
         I = I + 1
         K = K + K
      GO TO 400
  410 IF (IBETA .EQ. 10) GO TO 420
      IEXP = I + 1
      MX = K + K
      GO TO 450
C-----------------------------------------------------------------------
C  This segment is for decimal machines only.
C-----------------------------------------------------------------------
  420 IEXP = 2
      IZ = IBETA
  430 IF (K .LT. IZ) GO TO 440
         IZ = IZ * IBETA
         IEXP = IEXP + 1
      GO TO 430
  440 MX = IZ + IZ - 1
C-----------------------------------------------------------------------
C  Loop to determine MINEXP, XMIN.
C  Exit from loop is signaled by an underflow.
C-----------------------------------------------------------------------
  450 XMIN = Y
         Y = Y * BETAIN
C-----------------------------------------------------------------------
C  Check for underflow here.
C-----------------------------------------------------------------------
         A = Y * ONE
         TEMP = Y * T
         IF (((A+A) .EQ. ZERO) .OR. (ABS(Y) .GE. XMIN)) GO TO 460
         K = K + 1
         TEMP1 = TEMP * BETAIN
         IF (TEMP1*BETA .NE. Y) GO TO 450
      NXRES = 3
      XMIN = Y
  460 MINEXP = -K
C-----------------------------------------------------------------------
C  Determine MAXEXP, XMAX.
C-----------------------------------------------------------------------
      IF ((MX .GT. K+K-3) .OR. (IBETA .EQ. 10)) GO TO 500
      MX = MX + MX
      IEXP = IEXP + 1
  500 MAXEXP = MX + MINEXP
C-----------------------------------------------------------------
C  Adjust IRND to reflect partial underflow.
C-----------------------------------------------------------------
      IRND = IRND + NXRES
C-----------------------------------------------------------------
C  Adjust for IEEE-style machines.
C-----------------------------------------------------------------
      IF ((IRND .EQ. 2) .OR. (IRND .EQ. 5)) MAXEXP = MAXEXP - 2
C-----------------------------------------------------------------
C  Adjust for non-IEEE machines with partial underflow.
C-----------------------------------------------------------------
      IF ((IRND .EQ. 3) .OR. (IRND .EQ. 4)) MAXEXP = MAXEXP - IT
C-----------------------------------------------------------------
C  Adjust for machines with implicit leading bit in binary
C  significand, and machines with radix point at extreme
C  right of significand.
C-----------------------------------------------------------------
      I = MAXEXP + MINEXP
      IF ((IBETA .EQ. 2) .AND. (I .EQ. 0)) MAXEXP = MAXEXP - 1
      IF (I .GT. 20) MAXEXP = MAXEXP - 1
      IF (A .NE. Y) MAXEXP = MAXEXP - 2
      XMAX = ONE - EPSNEG
      IF (XMAX*ONE .NE. XMAX) XMAX = ONE - BETA * EPSNEG
      XMAX = XMAX / (BETA * BETA * BETA * XMIN)
      I = MAXEXP + MINEXP + 3
      IF (I .LE. 0) GO TO 520
      DO 510 J = 1, I
          IF (IBETA .EQ. 2) XMAX = XMAX + XMAX
          IF (IBETA .NE. 2) XMAX = XMAX * BETA
  510 CONTINUE
      DMACH(1) = XMIN
      DMACH(2) = XMAX
      DMACH(3) = EPSNEG
      DMACH(4) = EPS
  520 RETURN
C---------- LAST CARD OF T665D ----------
      END


      FUNCTION R1MACH(i)

*-----------------------------------------------------------------------------*
*= PURPOSE:                                                                  =*
*= R1MACH calculates various machine constants in single precision.          =*
*-----------------------------------------------------------------------------*
*= PARAMETERS:                                                               =*
*=   I       -  INTEGER, identifies the machine constant (0<I<5)         (I) =*
*=   R1MACH  -  REAL, machine constant in single precision               (O) =*
*=      I=1     - the smallest non-vanishing normalized floating-point       =*
*=                power of the radix, i.e., R1MACH=FLOAT(IBETA)**MINEXP      =*
*=      I=2     - the largest finite floating-point number.  In              =*
*=                particular R1MACH=(1.0-EPSNEG)*FLOAT(IBETA)**MAXEXP        =*
*=                Note - on some machines R1MACH will be only the            =*
*=                second, or perhaps third, largest number, being            =*
*=                too small by 1 or 2 units in the last digit of             =*
*=                the significand.                                           =*
*=      I=3     - A small positive floating-point number such that           =*
*=                1.0-R1MACH .NE. 1.0. In particular, if IBETA = 2           =*
*=                or  IRND = 0, R1MACH = FLOAT(IBETA)**NEGEPS.               =*
*=                Otherwise,  R1MACH = (IBETA**NEGEPS)/2.  Because           =*
*=                NEGEPS is bounded below by -(IT+3), R1MACH may not         =*
*=                be the smallest number that can alter 1.0 by               =*
*=                subtraction.                                               =*
*=      I=4     - the smallest positive floating-point number such           =*
*=                that  1.0+R1MACH .NE. 1.0. In particular, if either        =*
*=                IBETA = 2  or  IRND = 0, R1MACH=FLOAT(IBETA)**MACHEP.      =*
*=                Otherwise, R1MACH=(FLOAT(IBETA)**MACHEP)/2                 =*
*=  (see routine T665R for more information on different constants)          =*
*-----------------------------------------------------------------------------*

      REAL r1mach
      INTEGER i
   
      LOGICAL doinit
      DATA doinit/.TRUE./
      SAVE doinit

      REAL rmach(4) 
      SAVE rmach

      IF (( i .GE. 1 ) .AND. ( i .LE. 4 )) THEN
* compute constants at first call only
        IF (doinit) THEN
           CALL t665r(rmach)
           doinit = .FALSE.
        ENDIF
        r1mach = rmach(i)
      ELSE
        WRITE(0,*) '>>> ERROR (R1MACH) <<<  invalid argument'
        STOP
      ENDIF

      END


C      ALGORITHM 665, COLLECTED ALGORITHMS FROM ACM.
C      THIS WORK PUBLISHED IN TRANSACTIONS ON MATHEMATICAL SOFTWARE,
C      VOL. 14, NO. 4, PP. 303-311.
      SUBROUTINE T665R(RMACH)
C-----------------------------------------------------------------------
C  This Fortran 77 subroutine is intended to determine the parameters
C   of the floating-point arithmetic system specified below.  The
C   determination of the first three uses an extension of an algorithm
C   due to M. Malcolm, CACM 15 (1972), pp. 949-951, incorporating some,
C   but not all, of the improvements suggested by M. Gentleman and S.
C   Marovich, CACM 17 (1974), pp. 276-277.  An earlier version of this
C   program was published in the book Software Manual for the
C   Elementary Functions by W. J. Cody and W. Waite, Prentice-Hall,
C   Englewood Cliffs, NJ, 1980.
C
C  The program as given here must be modified before compiling.  If
C   a single (double) precision version is desired, change all
C   occurrences of CS (CD) in columns 1 and 2 to blanks.
C
C  Parameter values reported are as follows:
C
C       IBETA   - the radix for the floating-point representation
C       IT      - the number of base IBETA digits in the floating-point
C                 significand
C       IRND    - 0 if floating-point addition chops
C                 1 if floating-point addition rounds, but not in the
C                   IEEE style
C                 2 if floating-point addition rounds in the IEEE style
C                 3 if floating-point addition chops, and there is
C                   partial underflow
C                 4 if floating-point addition rounds, but not in the
C                   IEEE style, and there is partial underflow
C                 5 if floating-point addition rounds in the IEEE style,
C                   and there is partial underflow
C       NGRD    - the number of guard digits for multiplication with
C                 truncating arithmetic.  It is
C                 0 if floating-point arithmetic rounds, or if it
C                   truncates and only  IT  base  IBETA digits
C                   participate in the post-normalization shift of the
C                   floating-point significand in multiplication;
C                 1 if floating-point arithmetic truncates and more
C                   than  IT  base  IBETA  digits participate in the
C                   post-normalization shift of the floating-point
C                   significand in multiplication.
C       MACHEP  - the largest negative integer such that
C                 1.0+FLOAT(IBETA)**MACHEP .NE. 1.0, except that
C                 MACHEP is bounded below by  -(IT+3)
C       NEGEPS  - the largest negative integer such that
C                 1.0-FLOAT(IBETA)**NEGEPS .NE. 1.0, except that
C                 NEGEPS is bounded below by  -(IT+3)
C       IEXP    - the number of bits (decimal places if IBETA = 10)
C                 reserved for the representation of the exponent
C                 (including the bias or sign) of a floating-point
C                 number
C       MINEXP  - the largest in magnitude negative integer such that
C                 FLOAT(IBETA)**MINEXP is positive and normalized
C       MAXEXP  - the smallest positive power of  BETA  that overflows
C       EPS     - the smallest positive floating-point number such
C                 that  1.0+EPS .NE. 1.0. In particular, if either
C                 IBETA = 2  or  IRND = 0, EPS = FLOAT(IBETA)**MACHEP.
C                 Otherwise,  EPS = (FLOAT(IBETA)**MACHEP)/2
C       EPSNEG  - A small positive floating-point number such that
C                 1.0-EPSNEG .NE. 1.0. In particular, if IBETA = 2
C                 or  IRND = 0, EPSNEG = FLOAT(IBETA)**NEGEPS.
C                 Otherwise,  EPSNEG = (IBETA**NEGEPS)/2.  Because
C                 NEGEPS is bounded below by -(IT+3), EPSNEG may not
C                 be the smallest number that can alter 1.0 by
C                 subtraction.
C       XMIN    - the smallest non-vanishing normalized floating-point
C                 power of the radix, i.e.,  XMIN = FLOAT(IBETA)**MINEXP
C       XMAX    - the largest finite floating-point number.  In
C                 particular  XMAX = (1.0-EPSNEG)*FLOAT(IBETA)**MAXEXP
C                 Note - on some machines  XMAX  will be only the
C                 second, or perhaps third, largest number, being
C                 too small by 1 or 2 units in the last digit of
C                 the significand.
C
C     Latest revision - April 20, 1987
C
C     Author - W. J. Cody
C              Argonne National Laboratory
C
C-----------------------------------------------------------------------
      REAL rmach(4)
      INTEGER I,IBETA,IEXP,IRND,IT,ITEMP,IZ,J,K,MACHEP,MAXEXP,
     1        MINEXP,MX,NEGEP,NGRD,NXRES
      REAL A,B,BETA,BETAIN,BETAH,CONV,EPS,EPSNEG,ONE,T,TEMP,TEMPA,
     1     TEMP1,TWO,XMAX,XMIN,Y,Z,ZERO
CD    DOUBLE PRECISION A,B,BETA,BETAIN,BETAH,CONV,EPS,EPSNEG,ONE,
CD   1                 T,TEMP,TEMPA,TEMP1,TWO,XMAX,XMIN,Y,Z,ZERO
C-----------------------------------------------------------------------
      CONV(I) = REAL(I)
CD    CONV(I) = DBLE(I)
      ONE = CONV(1)
      TWO = ONE + ONE
      ZERO = ONE - ONE
C-----------------------------------------------------------------------
C  Determine IBETA, BETA ala Malcolm.
C-----------------------------------------------------------------------
      A = ONE
   10 A = A + A
         TEMP = A+ONE
         TEMP1 = TEMP-A
         IF (TEMP1-ONE .EQ. ZERO) GO TO 10
      B = ONE
   20 B = B + B
         TEMP = A+B
         ITEMP = INT(TEMP-A)
         IF (ITEMP .EQ. 0) GO TO 20
      IBETA = ITEMP
      BETA = CONV(IBETA)
C-----------------------------------------------------------------------
C  Determine IT, IRND.
C-----------------------------------------------------------------------
      IT = 0
      B = ONE
  100 IT = IT + 1
         B = B * BETA
         TEMP = B+ONE
         TEMP1 = TEMP-B
         IF (TEMP1-ONE .EQ. ZERO) GO TO 100
      IRND = 0
      BETAH = BETA / TWO
      TEMP = A+BETAH
      IF (TEMP-A .NE. ZERO) IRND = 1
      TEMPA = A + BETA
      TEMP = TEMPA+BETAH
      IF ((IRND .EQ. 0) .AND. (TEMP-TEMPA .NE. ZERO)) IRND = 2
C-----------------------------------------------------------------------
C  Determine NEGEP, EPSNEG.
C-----------------------------------------------------------------------
      NEGEP = IT + 3
      BETAIN = ONE / BETA
      A = ONE
      DO 200 I = 1, NEGEP
         A = A * BETAIN
  200 CONTINUE
      B = A
  210 TEMP = ONE-A
         IF (TEMP-ONE .NE. ZERO) GO TO 220
         A = A * BETA
         NEGEP = NEGEP - 1
      GO TO 210
  220 NEGEP = -NEGEP
      EPSNEG = A
      IF ((IBETA .EQ. 2) .OR. (IRND .EQ. 0)) GO TO 300
      A = (A*(ONE+A)) / TWO
      TEMP = ONE-A
      IF (TEMP-ONE .NE. ZERO) EPSNEG = A
C-----------------------------------------------------------------------
C  Determine MACHEP, EPS.
C-----------------------------------------------------------------------
  300 MACHEP = -IT - 3
      A = B
  310 TEMP = ONE+A
         IF (TEMP-ONE .NE. ZERO) GO TO 320
         A = A * BETA
         MACHEP = MACHEP + 1
      GO TO 310
  320 EPS = A
      TEMP = TEMPA+BETA*(ONE+EPS)
      IF ((IBETA .EQ. 2) .OR. (IRND .EQ. 0)) GO TO 350
      A = (A*(ONE+A)) / TWO
      TEMP = ONE+A
      IF (TEMP-ONE .NE. ZERO) EPS = A
C-----------------------------------------------------------------------
C  Determine NGRD.
C-----------------------------------------------------------------------
  350 NGRD = 0
      TEMP = ONE+EPS
      IF ((IRND .EQ. 0) .AND. (TEMP*ONE-ONE .NE. ZERO)) NGRD = 1
C-----------------------------------------------------------------------
C  Determine IEXP, MINEXP, XMIN.
C
C  Loop to determine largest I and K = 2**I such that
C         (1/BETA) ** (2**(I))
C  does not underflow.
C  Exit from loop is signaled by an underflow.
C-----------------------------------------------------------------------
      I = 0
      K = 1
      Z = BETAIN
      T = ONE + EPS
      NXRES = 0
  400 Y = Z
         Z = Y * Y
C-----------------------------------------------------------------------
C  Check for underflow here.
C-----------------------------------------------------------------------
         A = Z * ONE
         TEMP = Z * T
         IF ((A+A .EQ. ZERO) .OR. (ABS(Z) .GE. Y)) GO TO 410
         TEMP1 = TEMP * BETAIN
         IF (TEMP1*BETA .EQ. Z) GO TO 410
         I = I + 1
         K = K + K
      GO TO 400
  410 IF (IBETA .EQ. 10) GO TO 420
      IEXP = I + 1
      MX = K + K
      GO TO 450
C-----------------------------------------------------------------------
C  This segment is for decimal machines only.
C-----------------------------------------------------------------------
  420 IEXP = 2
      IZ = IBETA
  430 IF (K .LT. IZ) GO TO 440
         IZ = IZ * IBETA
         IEXP = IEXP + 1
      GO TO 430
  440 MX = IZ + IZ - 1
C-----------------------------------------------------------------------
C  Loop to determine MINEXP, XMIN.
C  Exit from loop is signaled by an underflow.
C-----------------------------------------------------------------------
  450 XMIN = Y
         Y = Y * BETAIN
C-----------------------------------------------------------------------
C  Check for underflow here.
C-----------------------------------------------------------------------
         A = Y * ONE
         TEMP = Y * T
         IF (((A+A) .EQ. ZERO) .OR. (ABS(Y) .GE. XMIN)) GO TO 460
         K = K + 1
         TEMP1 = TEMP * BETAIN
         IF (TEMP1*BETA .NE. Y) GO TO 450
      NXRES = 3
      XMIN = Y
  460 MINEXP = -K
C-----------------------------------------------------------------------
C  Determine MAXEXP, XMAX.
C-----------------------------------------------------------------------
      IF ((MX .GT. K+K-3) .OR. (IBETA .EQ. 10)) GO TO 500
      MX = MX + MX
      IEXP = IEXP + 1
  500 MAXEXP = MX + MINEXP
C-----------------------------------------------------------------
C  Adjust IRND to reflect partial underflow.
C-----------------------------------------------------------------
      IRND = IRND + NXRES
C-----------------------------------------------------------------
C  Adjust for IEEE-style machines.
C-----------------------------------------------------------------
      IF ((IRND .EQ. 2) .OR. (IRND .EQ. 5)) MAXEXP = MAXEXP - 2
C-----------------------------------------------------------------
C  Adjust for non-IEEE machines with partial underflow.
C-----------------------------------------------------------------
      IF ((IRND .EQ. 3) .OR. (IRND .EQ. 4)) MAXEXP = MAXEXP - IT
C-----------------------------------------------------------------
C  Adjust for machines with implicit leading bit in binary
C  significand, and machines with radix point at extreme
C  right of significand.
C-----------------------------------------------------------------
      I = MAXEXP + MINEXP
      IF ((IBETA .EQ. 2) .AND. (I .EQ. 0)) MAXEXP = MAXEXP - 1
      IF (I .GT. 20) MAXEXP = MAXEXP - 1
      IF (A .NE. Y) MAXEXP = MAXEXP - 2
      XMAX = ONE - EPSNEG
      IF (XMAX*ONE .NE. XMAX) XMAX = ONE - BETA * EPSNEG
      XMAX = XMAX / (BETA * BETA * BETA * XMIN)
      I = MAXEXP + MINEXP + 3
      IF (I .LE. 0) GO TO 520
      DO 510 J = 1, I
          IF (IBETA .EQ. 2) XMAX = XMAX + XMAX
          IF (IBETA .NE. 2) XMAX = XMAX * BETA
  510 CONTINUE
      RMACH(1) = XMIN
      RMACH(2) = XMAX
      RMACH(3) = EPSNEG
      RMACH(4) = EPS
  520 RETURN
C---------- LAST CARD OF T665R ----------
      END

      SUBROUTINE rdetfl(nw,wl,f)

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Read and re-grid extra-terrestrial flux data.                            =*
*-----------------------------------------------------------------------------*
*=  PARAMETERS:                                                              =*
*=  NW     - INTEGER, number of specified intervals + 1 in working        (I)=*
*=           wavelength grid                                                 =*
*=  WL     - REAL, vector of lower limits of wavelength intervals in      (I)=*
*=           working wavelength grid                                         =*
*=  F      - REAL, spectral irradiance at the top of the atmosphere at    (O)=*
*=           each specified wavelength                                       =*
*-----------------------------------------------------------------------------*
*=  EDIT HISTORY:                                                            =*
*=  02/97  Changed offset for grid-end interpolation to relative number      =*
*=         (x * (1 +- deltax))                                               =*
*=  05/96  Put in different preset options                                   =*
*-----------------------------------------------------------------------------*
*= This program is free software;  you can redistribute it and/or modify     =*
*= it under the terms of the GNU General Public License as published by the  =*
*= Free Software Foundation;  either version 2 of the license, or (at your   =*
*= option) any later version.                                                =*
*= The TUV package is distributed in the hope that it will be useful, but    =*
*= WITHOUT ANY WARRANTY;  without even the implied warranty of MERCHANTIBI-  =*
*= LITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public     =*
*= License for more details.                                                 =*
*= To obtain a copy of the GNU General Public License, write to:             =*
*= Free Software Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.   =*
*-----------------------------------------------------------------------------*
*= To contact the authors, please mail to:                                   =*
*= Sasha Madronich, NCAR/ACD, P.O.Box 3000, Boulder, CO, 80307-3000, USA  or =*
*= send email to:  sasha@ucar.edu                                            =*
*-----------------------------------------------------------------------------*
*= Copyright (C) 1994,95,96  University Corporation for Atmospheric Research =*
*-----------------------------------------------------------------------------*

      IMPLICIT NONE
      INCLUDE 'tuv.params'

      integer kdata
      parameter(kdata=20000)

* input: (wavelength grid)
      INTEGER nw
      REAL wl(kw)
      INTEGER iw

* output: (extra terrestrial solar flux)
      REAL f(kw)

* INTERNAL:

* work arrays for input data files:

      CHARACTER*40 fil
      REAL x1(kdata), x2(kdata)
      REAL y1(kdata), y2(kdata)
      INTEGER nhead, n, i, ierr
      REAL dum

* data gridded onto wl(kw) grid:

      REAL yg1(kw)
      REAL yg2(kw)
      REAL yg3(kw)

      REAL hc
      PARAMETER(hc = 6.62E-34 * 2.998E8)

      INTEGER msun

      REAL refrac
      EXTERNAL refrac


*_______________________________________________________________________
* select desired extra-terrestrial solar irradiance, using msun:

*  1 =   extsol.flx:  De Luisi, JGR 80, 345-354, 1975
*                     280-400 nm, 1 nm steps.
*  2 =   lowsun3.flx:  Lowtran (John Bahr, priv. comm.)
*                      173.974-500000 nm, ca. 0.1 nm steps in UV-B
*  3 =   modtran1.flx:  Modtran (Gail Anderson, priv. comm.)
*                       200.55-949.40, 0.05 nm steps
*  4 =   nicolarv.flx:  wvl<300 nm from Nicolet, Plan. Sp. Sci., 29,  951-974, 1981.
*                       wvl>300 nm supplied by Thekaekera, Arvesen Applied Optics 8, 
*                       11, 2215-2232, 1969 (also see Thekaekera, Applied Optics, 13,
*                       3, 518, 1974) but with corrections recommended by:
*                       Nicolet, Plan. Sp. Sci., 37, 1249-1289, 1989.
*                       270.0-299.0 nm in 0.5 nm steps
*                       299.6-340.0 nm in ca. 0.4 nm steps
*                       340.0-380.0 nm in ca. 0.2 nm steps
*                       380.0-470.0 nm in ca. 0.1 nm steps   
*  5 =  solstice.flx:  From:   MX%"ROTTMAN@virgo.hao.ucar.edu" 12-OCT-1994 13:03:01.62
*                      Original data gave Wavelength in vacuum
*                      (Converted to wavelength in air using Pendorf, 1967, J. Opt. Soc. Am.)
*                      279.5 to 420 nm, 0.24 nm spectral resolution, approx 0.07 nm steps
*  6 =  suntoms.flx: (from TOMS CD-ROM).  280-340 nm, 0.05 nm steps.
*  7 =  neckel.flx:  H.Neckel and D.Labs, "The Solar Radiation Between 3300 and 12500 A",
*                    Solar Physics v.90, pp.205-258 (1984).
*                    1 nm between 330.5 and 529.5 nm
*                    2 nm between 631.0 and 709.0 nm
*                    5 nm between 872.5 and 1247.4 nm
*                    Units: must convert to W m-2 nm-1 from photons cm-2 s-1 nm-1
*  8 =  atlas3.flx:  ATLAS3-SUSIM 13 Nov 94 high resolution (0.15 nm FWHM)
*                    available by ftp from susim.nrl.navy.mil
*                    atlas3_1994_317_a.dat, downloaded 30 Sept 98.
*                    150-407.95 nm, in 0.05 nm steps
*                    (old version from Dianne Prinz through Jim Slusser)
*                    orig wavelengths in vac, correct here to air.
*  9 =  solstice.flx:  solstice 1991-1996, average
*                    119.5-420.5 nm in 1 nm steps

* 10 =  susim_hi.flx:  SUSIM SL2 high resolution
*                      120.5-400.0 in 0.05 nm intervals (0.15 nm resolution)
* 11 =  wmo85.flx: from WMO 1995 Ozone Atmospehric Ozone (report no. 16)
*                  on variable-size bins.  Original values are per bin, not
*                  per nm.
* 12 = combine susim_hi.flx for .lt. 350 nm, neckel.flx for .gt. 350 nm.
*
* 13 = combine 
*     for wl(iw) .lt. 150.01                                susim_hi.flx
*     for wl(iw) .ge. 150.01 and wl(iw) .le. 400            atlas3.flx 
*     for wl(iw) .gt. 400                                   Neckel & Labs 

      msun = 13

* simple files are read and interpolated here in-line. Reading of 
* more complex files may be done with longer code in a read#.f subroutine.

      IF (msun .EQ. 1) THEN
         fil = 'DATAE1/SUN/extsol.flx'
c         write(kout,*) fil
         OPEN(UNIT=kin,FILE=fil,STATUS='old')
         nhead = 3
         n =121
         DO i = 1, nhead
            READ(kin,*)
         ENDDO
         DO i = 1, n
            READ(kin,*) x1(i), y1(i)
         ENDDO
         CLOSE (kin)
         CALL addpnt(x1,y1,kdata,n,x1(1)*(1.-deltax),0.)
         CALL addpnt(x1,y1,kdata,n,          0.,0.)
         CALL addpnt(x1,y1,kdata,n,x1(n)*(1.+deltax),0.)
         CALL addpnt(x1,y1,kdata,n,      1.e+38,0.)
         CALL inter2(nw,wl,yg1,n,x1,y1,ierr)
         IF (ierr .NE. 0) THEN
            WRITE(*,*) ierr, fil
            STOP
         ENDIF         
         DO iw = 1, nw-1
            f(iw) = yg1(iw)
         ENDDO

      ELSEIF (msun .EQ. 2) THEN
         fil = 'DATAE1/SUN/lowsun3.flx'
c         write(kout,*) fil
         OPEN(UNIT=kin,FILE=fil,STATUS='old')
         nhead = 3
         n = 4327
         DO i = 1, nhead
            READ(kin,*)
         ENDDO
         DO i = 1, n
            READ(kin,*) x1(i), y1(i)
         ENDDO
         CLOSE (kin)
         CALL addpnt(x1,y1,kdata,n,x1(1)*(1.-deltax),0.)
         CALL addpnt(x1,y1,kdata,n,          0.,0.)
         CALL addpnt(x1,y1,kdata,n,x1(n)*(1.+deltax),0.)
         CALL addpnt(x1,y1,kdata,n,      1.e+38,0.)
         CALL inter2(nw,wl,yg1,n,x1,y1,ierr)
         IF (ierr .NE. 0) THEN
            WRITE(*,*) ierr, fil
            STOP
         ENDIF         
         DO iw = 1, nw-1
            f(iw) = yg1(iw)
         ENDDO

      ELSEIF (msun .EQ. 3) THEN
         fil = 'DATAE1/SUN/modtran1.flx'
c         write(kout,*) fil
         OPEN(UNIT=kin,FILE=fil,STATUS='old')
         nhead = 6
         n = 14980
         DO i = 1, nhead
            READ(kin,*)
         ENDDO
         DO i = 1, n
            READ(kin,*) x1(i), y1(i)
         ENDDO
         CLOSE (kin)
         CALL addpnt(x1,y1,kdata,n,x1(1)*(1.-deltax),0.)
         CALL addpnt(x1,y1,kdata,n,          0.,0.)
         CALL addpnt(x1,y1,kdata,n,x1(n)*(1.+deltax),0.)
         CALL addpnt(x1,y1,kdata,n,      1.e+38,0.)
         CALL inter2(nw,wl,yg1,n,x1,y1,ierr)
         IF (ierr .NE. 0) THEN
            WRITE(*,*) ierr, fil
            STOP
         ENDIF         
         DO iw = 1, nw-1
            f(iw) = yg1(iw)
         ENDDO

      ELSEIF (msun .EQ. 4) THEN
         fil = 'DATAE1/SUN/nicolarv.flx'
c         write(kout,*) fil
         OPEN(UNIT=kin,FILE=fil,STATUS='old')
         nhead = 8
         n = 1260
         DO i = 1, nhead
            READ(kin,*)
         ENDDO
         DO i = 1, n
            READ(kin,*) x1(i), y1(i)
         ENDDO
         CLOSE (kin)
         CALL addpnt(x1,y1,kdata,n,x1(1)*(1.-deltax),0.)
         CALL addpnt(x1,y1,kdata,n,          0.,0.)
         CALL addpnt(x1,y1,kdata,n,x1(n)*(1.+deltax),0.)
         CALL addpnt(x1,y1,kdata,n,      1.e+38,0.)
         CALL inter2(nw,wl,yg1,n,x1,y1,ierr)
         IF (ierr .NE. 0) THEN
            WRITE(*,*) ierr, fil
            STOP
         ENDIF         
         DO iw = 1, nw-1
            f(iw) = yg1(iw)
         ENDDO

      ELSEIF (msun .EQ. 5) THEN
* unofficial - do not use
         fil = 'DATAE2/SUN/solstice.flx'
c         write(kout,*) fil
         OPEN(UNIT=kin,FILE=fil,STATUS='old')
         nhead = 11
         n = 2047
         DO i = 1, nhead
            READ(kin,*)
         ENDDO
         DO i = 1, n
            READ(kin,*) x1(i), y1(i)
         ENDDO
         CLOSE (kin)
         CALL addpnt(x1,y1,kdata,n,x1(1)*(1.-deltax),0.)
         CALL addpnt(x1,y1,kdata,n,          0.,0.)
         CALL addpnt(x1,y1,kdata,n,x1(n)*(1.+deltax),0.)
         CALL addpnt(x1,y1,kdata,n,      1.e+38,0.)
         CALL inter2(nw,wl,yg1,n,x1,y1,ierr)
         IF (ierr .NE. 0) THEN
            WRITE(*,*) ierr, fil
            STOP
         ENDIF         
         DO iw = 1, nw-1
            f(iw) = yg1(iw)
         ENDDO

      ELSEIF (msun .EQ. 6) THEN
* unofficial - do not use
         fil = 'DATAE2/SUN/suntoms.flx'
c         write(kout,*) fil
         OPEN(UNIT=kin,FILE=fil,STATUS='old')
         nhead = 3
         n = 1200
         DO i = 1, nhead
            READ(kin,*)
         ENDDO
         DO i = 1, n
            READ(kin,*) x1(i), y1(i)
            y1(i) = y1(i)* 1.e-3
         ENDDO
         CLOSE (kin)
         CALL addpnt(x1,y1,kdata,n,x1(1)*(1.-deltax),0.)
         CALL addpnt(x1,y1,kdata,n,          0.,0.)
         CALL addpnt(x1,y1,kdata,n,x1(n)*(1.+deltax),0.)
         CALL addpnt(x1,y1,kdata,n,      1.e+38,0.)
         CALL inter2(nw,wl,yg1,n,x1,y1,ierr)
         IF (ierr .NE. 0) THEN
            WRITE(*,*) ierr, fil
            STOP
         ENDIF         
         DO iw = 1, nw-1
            f(iw) = yg1(iw)
         ENDDO

      ELSEIF (msun .EQ. 7) THEN
         fil = 'DATAE1/SUN/neckel.flx'
c         write(kout,*) fil
         OPEN(UNIT=kin,FILE=fil,STATUS='old')
         nhead = 11
         n = 496
         DO i = 1, nhead
            READ(kin,*)
         ENDDO
         DO i = 1, n
            READ(kin,*) dum, y1(i)
            if (dum .lt. 630.0) x1(i) = dum - 0.5
            if (dum .gt. 630.0 .and. dum .lt. 870.0) x1(i) = dum - 1.0
            if (dum .gt. 870.0) x1(i) = dum - 2.5
            y1(i) = y1(i) * 1.E4 * hc / (dum * 1.E-9)
         ENDDO
         CLOSE (kin)
         x1(n+1) = x1(n) + 2.5
         do i = 1, n
            y1(i) = y1(i) * (x1(i+1)-x1(i))
         enddo
         call inter3(nw,wl,yg2,n+1,x1,y1,0)
         do iw = 1, nw-1
            yg1(iw) = yg1(iw) / (wl(iw+1)-wl(iw))
         enddo
         DO iw = 1, nw-1
            f(iw) = yg1(iw)
         ENDDO

      ELSEIF (msun .EQ. 8) THEN
         nhead = 5
         fil = 'DATAE1/SUN/atlas3_1994_317_a.dat'
c         write(kout,*) fil
         OPEN(UNIT=kin,FILE=fil,STATUS='old')
         nhead = 13
         n = 5160
         DO i = 1, nhead
            READ(kin,*)
         ENDDO
         DO i = 1, n
            READ(kin,*) x1(i), y1(i)
            y1(i) = y1(i) * 1.E-3
            x1(i) = x1(i)/refrac(x1(i))
         ENDDO
         CLOSE (kin)
         CALL addpnt(x1,y1,kdata,n,x1(1)*(1.-deltax),0.)
         CALL addpnt(x1,y1,kdata,n,          0.,0.)
         CALL addpnt(x1,y1,kdata,n,x1(n)*(1.+deltax),0.)
         CALL addpnt(x1,y1,kdata,n,      1.e+38,0.)
         CALL inter2(nw,wl,yg1,n,x1,y1,ierr)
         IF (ierr .NE. 0) THEN
            WRITE(*,*) ierr, fil
            STOP
         ENDIF         
         DO iw = 1, nw-1
            f(iw) = yg1(iw)
         ENDDO

      ELSEIF (msun .EQ. 9) THEN
         fil = 'DATAE1/SUN/solstice.flx'
c         write(kout,*) fil
         OPEN(UNIT=kin,FILE=fil,STATUS='old')
         nhead = 2
         n = 302
         DO i = 1, nhead
            READ(kin,*)
         ENDDO
         DO i = 1, n
            READ(kin,*) x1(i), y1(i)
         ENDDO
         CLOSE (kin)
         CALL addpnt(x1,y1,kdata,n,x1(1)*(1.-deltax),0.)
         CALL addpnt(x1,y1,kdata,n,          0.,0.)
         CALL addpnt(x1,y1,kdata,n,x1(n)*(1.+deltax),0.)
         CALL addpnt(x1,y1,kdata,n,      1.e+38,0.)
         CALL inter2(nw,wl,yg1,n,x1,y1,ierr)
         IF (ierr .NE. 0) THEN
            WRITE(*,*) ierr, fil
            STOP
         ENDIF         
         DO iw = 1, nw-1
            f(iw) = yg1(iw)
         ENDDO

      ELSEIF (msun .EQ. 10) THEN
c         WRITE(kout,*) 'DATAE1/SUN/susim_hi.flx'
         CALL read1(nw,wl,yg1)
         DO iw = 1, nw-1
            f(iw) = yg1(iw)
         ENDDO


      ELSEIF (msun .EQ. 11) THEN
c         WRITE(kout,*) 'DATAE1/SUN/wmo85.flx'
         CALL read2(nw,wl,yg1)
         DO iw = 1, nw-1
            f(iw) = yg1(iw)
         ENDDO

      ELSEIF (msun .EQ. 12) THEN
c         WRITE(kout,*) 'DATAE1/SUN/susim_hi.flx'
         CALL read1(nw,wl,yg1)
         fil = 'DATAE1/SUN/neckel.flx'
c         write(kout,*) fil
         OPEN(UNIT=kin,FILE=fil,STATUS='old')
         nhead = 11
         n = 496
         DO i = 1, nhead
            READ(kin,*)
         ENDDO
         DO i = 1, n
            READ(kin,*) dum, y1(i)
            if (dum .lt. 630.0) x1(i) = dum - 0.5
            if (dum .gt. 630.0 .and. dum .lt. 870.0) x1(i) = dum - 1.0
            if (dum .gt. 870.0) x1(i) = dum - 2.5
            y1(i) = y1(i) * 1.E4 * hc / (dum * 1.E-9)
         ENDDO
         CLOSE (kin)
         x1(n+1) = x1(n) + 2.5
         do i = 1, n
            y1(i) = y1(i) * (x1(i+1)-x1(i))
         enddo
         call inter3(nw,wl,yg2,n+1,x1,y1,0)
         do iw = 1, nw-1
            yg2(iw) = yg2(iw) / (wl(iw+1)-wl(iw))
         enddo

         DO iw = 1, nw-1
            IF (wl(iw) .GT. 350.) THEN
               f(iw) = yg2(iw)
            ELSE
               f(iw) = yg1(iw)
            ENDIF
         ENDDO

      ELSEIF (msun .EQ. 13) THEN

c         WRITE(kout,*) 'DATAE1/SUN/susim_hi.flx'
         CALL read1(nw,wl,yg1)

         fil = 'DATAE1/SUN/atlas3_1994_317_a.dat'
	 call read_atlas(x1,y1,kdata,n)         ! read atlas extra-terrestrial data
         DO i = 1, n
            y1(i) = y1(i) * 1.E-3
            x1(i) = x1(i)/refrac(x1(i))
         ENDDO

         CALL addpnt(x1,y1,kdata,n,x1(1)*(1.-deltax),0.)
         CALL addpnt(x1,y1,kdata,n,          0.,0.)
         CALL addpnt(x1,y1,kdata,n,x1(n)*(1.+deltax),0.)
         CALL addpnt(x1,y1,kdata,n,      1.e+38,0.)
         CALL inter2(nw,wl,yg2,n,x1,y1,ierr)
         IF (ierr .NE. 0) THEN
            WRITE(*,*) ierr, fil
            STOP
         ENDIF         

         fil = 'DATAE1/SUN/neckel.flx'
	 call read_neck(x1,y1,kdata,n)         ! read neckel extra-terrestrial data

         x1(n+1) = x1(n) + 2.5
         call inter4(nw,wl,yg3,n+1,x1,y1,0)

         DO iw = 1, nw-1

            IF (wl(iw) .LT. 150.01) THEN
               f(iw) = yg1(iw)
            ELSE IF ((wl(iw) .GE. 150.01) .AND. wl(iw) .LE. 400.) THEN
               f(iw) = yg2(iw)
            ELSE IF (wl(iw) .GT. 400.) THEN
               f(iw) = yg3(iw)
            ENDIF

         ENDDO


      ENDIF


      
*_______________________________________________________________________

      RETURN
      END

      SUBROUTINE rdno2xs(nw,wl,xsno2)

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Read NO2 molecular absorption cross section.  Re-grid data to match      =*
*=  specified wavelength working grid.                                       =*
*-----------------------------------------------------------------------------*
*=  PARAMETERS:                                                              =*
*=  NW     - INTEGER, number of specified intervals + 1 in working        (I)=*
*=           wavelength grid                                                 =*
*=  WL     - REAL, vector of lower limits of wavelength intervals in      (I)=*
*=           working wavelength grid                                         =*
*=  XSNO2  - REAL, molecular absoprtion cross section (cm^2) of NO2 at    (O)=*
*=           each specified wavelength                                       =*
*-----------------------------------------------------------------------------*
*=  EDIT HISTORY:                                                            =*
*=  02/97  Changed offset for grid-end interpolation to relative number      =*
*=         (x * (1 +- deltax)                                                =*
*-----------------------------------------------------------------------------*
*= This program is free software;  you can redistribute it and/or modify     =*
*= it under the terms of the GNU General Public License as published by the  =*
*= Free Software Foundation;  either version 2 of the license, or (at your   =*
*= option) any later version.                                                =*
*= The TUV package is distributed in the hope that it will be useful, but    =*
*= WITHOUT ANY WARRANTY;  without even the implied warranty of MERCHANTIBI-  =*
*= LITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public     =*
*= License for more details.                                                 =*
*= To obtain a copy of the GNU General Public License, write to:             =*
*= Free Software Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.   =*
*-----------------------------------------------------------------------------*
*= To contact the authors, please mail to:                                   =*
*= Sasha Madronich, NCAR/ACD, P.O.Box 3000, Boulder, CO, 80307-3000, USA  or =*
*= send email to:  sasha@ucar.edu                                            =*
*-----------------------------------------------------------------------------*
*= Copyright (C) 1994,95,96  University Corporation for Atmospheric Research =*
*-----------------------------------------------------------------------------*

      IMPLICIT NONE
      INCLUDE 'tuv.params'

      INTEGER kdata
      PARAMETER(kdata=1000)

* input: (altitude working grid)
      INTEGER nw
      REAL wl(kw)

* output:

      REAL xsno2(kw)

* local:
      integer nno2
      parameter(nno2=750)
      REAL x1(kdata)
      REAL y1(kdata),y2(kdata),y3(kdata)       
      data x1(1:nno2)/263.80,264.30,264.80,265.30,265.80,
     A 266.40,266.90,267.40,267.90,268.40,268.90,269.40,270.00,270.50,
     A 271.00,271.50,272.00,272.50,273.00,273.60,274.10,274.60,275.10,
     A 275.60,276.10,276.60,277.20,277.70,278.20,278.70,279.20,279.70,
     A 280.20,280.80,281.30,281.80,282.30,282.80,283.30,283.80,284.40,
     A 284.90,285.40,285.90,286.40,286.90,287.40,287.90,288.50,289.00,
     A 289.50,290.00,290.50,291.00,291.50,292.10,292.60,293.10,293.60,
     A 294.10,294.60,295.10,295.70,296.20,296.70,297.20,297.70,298.20,
     A 298.70,299.30,299.80,300.30,300.80,301.30,301.80,302.30,302.90,
     A 303.40,303.90,304.40,304.90,305.40,305.90,306.50,307.00,307.50,
     A 308.00,308.50,309.00,309.50,310.10,310.60,311.10,311.60,312.10,
     A 312.60,313.10,313.60,314.20,314.70,315.20,315.70,316.20,316.70,
     A 317.20,317.80,318.30,318.80,319.30,319.80,320.30,320.80,321.40,
     A 321.90,322.40,322.90,323.40,323.90,324.40,325.00,325.50,326.00,
     A 326.50,327.00,327.50,328.00,328.60,329.10,329.60,330.10,330.60,
     A 331.10,331.60,332.20,332.70,333.20,333.70,334.20,334.70,335.20,
     A 335.70,336.30,336.80,337.30,337.80,338.30,338.80,339.30,339.90,
     A 340.40,340.90,341.40,341.90,342.40,342.90,343.50,344.00,344.50,
     A 345.00,345.50,346.00,346.50,347.10,347.60,348.10,348.60,349.10,
     A 349.60,350.10,350.70,351.20,351.70,352.20,352.70,353.20,353.70,
     A 354.30,354.80,355.30,355.80,356.30,356.80,357.30,357.80,358.40,
     A 358.90,359.40,359.90,360.40,360.90,361.40,362.00,362.50,363.00,
     A 363.50,364.00,364.50,365.00,365.60,366.10,366.60,367.10,367.60,
     A 368.10,368.60,369.20,369.70,370.20,370.70,371.20,371.70,372.20,
     A 372.80,373.30,373.80,374.30,374.80,375.30,375.80,376.40,376.90,
     A 377.40,377.90,378.40,378.90,379.40,379.90,380.50,381.00,381.50,
     A 382.00,382.50,383.00,383.50,384.10,384.60,385.10,385.60,386.10,
     A 386.60,387.10,387.70,388.20,388.70,389.20,389.70,390.20,390.70,
     A 391.30,391.80,392.30,392.80,393.30,393.80,394.30,394.90,395.40,
     A 395.90,396.40,396.90,397.40,397.90,398.50,399.00,399.50,400.00,
     A 400.50,401.00,401.50,402.10,402.60,403.10,403.60,404.10,404.60,
     A 405.10,405.60,406.20,406.70,407.20,407.70,408.20,408.70,409.20,
     A 409.80,410.30,410.80,411.30,411.80,412.30,412.80,413.40,413.90,
     A 414.40,414.90,415.40,415.90,416.40,417.00,417.50,418.00,418.50,
     A 419.00,419.50,420.00,420.60,421.10,421.60,422.10,422.60,423.10,
     A 423.60,424.20,424.70,425.20,425.70,426.20,426.70,427.20,427.70,
     A 428.30,428.80,429.30,429.80,430.30,430.80,431.30,431.90,432.40,
     A 432.90,433.40,433.90,434.40,434.90,435.50,436.00,436.50,437.00,
     A 437.50,438.00,438.50,439.10,439.60,440.10,440.60,441.10,441.60,
     A 442.10,442.70,443.20,443.70,444.20,444.70,445.20,445.70,446.30,
     A 446.80,447.30,447.80,448.30,448.80,449.30,449.80,450.40,450.90,
     A 451.40,451.90,452.40,452.90,453.40,454.00,454.50,455.00,455.50,
     A 456.00,456.50,457.00,457.60,458.10,458.60,459.10,459.60,460.10,
     A 460.60,461.20,461.70,462.20,462.70,463.20,463.70,464.20,464.80,
     A 465.30,465.80,466.30,466.80,467.30,467.80,468.40,468.90,469.40,
     A 469.90,470.40,470.90,471.40,472.00,472.50,473.00,473.50,474.00,
     A 474.50,475.00,475.50,476.10,476.60,477.10,477.60,478.10,478.60,
     A 479.10,479.70,480.20,480.70,481.20,481.70,482.20,482.70,483.30,
     A 483.80,484.30,484.80,485.30,485.80,486.30,486.90,487.40,487.90,
     A 488.40,488.90,489.40,489.90,490.50,491.00,491.50,492.00,492.50,
     A 493.00,493.50,494.10,494.60,495.10,495.60,496.10,496.60,497.10,
     A 497.60,498.20,498.70,499.20,499.70,500.20,500.70,501.20,501.80,
     A 502.30,502.80,503.30,503.80,504.30,504.80,505.40,505.90,506.40,
     A 506.90,507.40,507.90,508.40,509.00,509.50,510.00,510.50,511.00,
     A 511.50,512.00,512.60,513.10,513.60,514.10,514.60,515.10,515.60,
     A 516.20,516.70,517.20,517.70,518.20,518.70,519.20,519.70,520.30,
     A 520.80,521.30,521.80,522.30,522.80,523.30,523.90,524.40,524.90,
     A 525.40,525.90,526.40,526.90,527.50,528.00,528.50,529.00,529.50,
     A 530.00,530.50,531.10,531.60,532.10,532.60,533.10,533.60,534.10,
     A 534.70,535.20,535.70,536.20,536.70,537.20,537.70,538.30,538.80,
     A 539.30,539.80,540.30,540.80,541.30,541.80,542.40,542.90,543.40,
     A 543.90,544.40,544.90,545.40,546.00,546.50,547.00,547.50,548.00,
     A 548.50,549.00,549.60,550.10,550.60,551.10,551.60,552.10,552.60,
     A 553.20,553.70,554.20,554.70,555.20,555.70,556.20,556.80,557.30,
     A 557.80,558.30,558.80,559.30,559.80,560.40,560.90,561.40,561.90,
     A 562.40,562.90,563.40,564.00,564.50,565.00,565.50,566.00,566.50,
     A 567.00,567.50,568.10,568.60,569.10,569.60,570.10,570.60,571.10,
     A 571.70,572.20,572.70,573.20,573.70,574.20,574.70,575.30,575.80,
     A 576.30,576.80,577.30,577.80,578.30,578.90,579.40,579.90,580.40,
     A 580.90,581.40,581.90,582.50,583.00,583.50,584.00,584.50,585.00,
     A 585.50,586.10,586.60,587.10,587.60,588.10,588.60,589.10,589.60,
     A 590.20,590.70,591.20,591.70,592.20,592.70,593.20,593.80,594.30,
     A 594.80,595.30,595.80,596.30,596.80,597.40,597.90,598.40,598.90,
     A 599.40,599.90,600.40,601.00,601.50,602.00,602.50,603.00,603.50,
     A 604.00,604.60,605.10,605.60,606.10,606.60,607.10,607.60,608.20,
     A 608.70,609.20,609.70,610.20,610.70,611.20,611.70,612.30,612.80,
     A 613.30,613.80,614.30,614.80,615.30,615.90,616.40,616.90,617.40,
     A 617.90,618.40,618.90,619.50,620.00,620.50,621.00,621.50,622.00,
     A 622.50,623.10,623.60,624.10,624.60,625.10,625.60,626.10,626.70,
     A 627.20,627.70,628.20,628.70,629.20,629.70,630.30,630.80,631.30,
     A 631.80,632.30,632.80,633.30,633.80,634.40,634.90,635.40,635.90,
     A 636.40,636.90,637.40,638.00,638.50,639.00,639.50,640.00,640.50,
     A 641.00,641.60,642.10,642.60,643.10,643.60,644.10,644.60,645.20,
     A 645.70,646.20,646.70,647.20,647.70,648.20,648.80/
       data y1(1:nno2)/.2988E-19,.3611E-19,.3128E-19,.3475E-19,
     A .3233E-19,.3636E-19,.3357E-19,.3625E-19,.3531E-19,.3962E-19,
     A .3703E-19,.3937E-19,.3787E-19,.4042E-19,.3785E-19,.4137E-19,
     A .4089E-19,.4396E-19,.4232E-19,.4590E-19,.4489E-19,.4728E-19,
     A .4567E-19,.4781E-19,.4781E-19,.5054E-19,.4884E-19,.5155E-19,
     A .5202E-19,.5575E-19,.5494E-19,.5747E-19,.5749E-19,.6014E-19,
     A .5869E-19,.6071E-19,.6055E-19,.6432E-19,.6528E-19,.6801E-19,
     A .6905E-19,.7337E-19,.7329E-19,.7505E-19,.7423E-19,.7668E-19,
     A .7646E-19,.7897E-19,.8037E-19,.8438E-19,.8523E-19,.8826E-19,
     A .9027E-19,.9445E-19,.9527E-19,.9710E-19,.9648E-19,.9868E-19,
     A .9894E-19,.1018E-18,.1028E-18,.1066E-18,.1093E-18,.1157E-18,
     A .1194E-18,.1229E-18,.1221E-18,.1227E-18,.1231E-18,.1257E-18,
     A .1251E-18,.1274E-18,.1298E-18,.1363E-18,.1399E-18,.1450E-18,
     A .1475E-18,.1508E-18,.1515E-18,.1541E-18,.1548E-18,.1579E-18,
     A .1592E-18,.1631E-18,.1644E-18,.1678E-18,.1707E-18,.1767E-18,
     A .1799E-18,.1831E-18,.1840E-18,.1887E-18,.1921E-18,.1975E-18,
     A .1979E-18,.1990E-18,.1997E-18,.2051E-18,.2081E-18,.2127E-18,
     A .2136E-18,.2169E-18,.2207E-18,.2278E-18,.2322E-18,.2356E-18,
     A .2349E-18,.2390E-18,.2424E-18,.2477E-18,.2478E-18,.2496E-18,
     A .2535E-18,.2602E-18,.2614E-18,.2647E-18,.2654E-18,.2716E-18,
     A .2759E-18,.2816E-18,.2837E-18,.2877E-18,.2882E-18,.2934E-18,
     A .2932E-18,.2988E-18,.3015E-18,.3069E-18,.3067E-18,.3099E-18,
     A .3103E-18,.3177E-18,.3259E-18,.3377E-18,.3368E-18,.3340E-18,
     A .3307E-18,.3372E-18,.3386E-18,.3440E-18,.3449E-18,.3484E-18,
     A .3457E-18,.3542E-18,.3656E-18,.3774E-18,.3757E-18,.3787E-18,
     A .3816E-18,.3912E-18,.3886E-18,.3832E-18,.3742E-18,.3772E-18,
     A .3794E-18,.3871E-18,.3915E-18,.4001E-18,.4045E-18,.4154E-18,
     A .4187E-18,.4253E-18,.4296E-18,.4435E-18,.4394E-18,.4296E-18,
     A .4188E-18,.4234E-18,.4275E-18,.4354E-18,.4314E-18,.4309E-18,
     A .4267E-18,.4382E-18,.4549E-18,.4735E-18,.4733E-18,.4683E-18,
     A .4627E-18,.4750E-18,.4860E-18,.4938E-18,.4852E-18,.4771E-18,
     A .4658E-18,.4679E-18,.4700E-18,.4828E-18,.4883E-18,.4952E-18,
     A .4913E-18,.4955E-18,.4925E-18,.4982E-18,.5016E-18,.5178E-18,
     A .5279E-18,.5310E-18,.5167E-18,.5150E-18,.5155E-18,.5239E-18,
     A .5200E-18,.5209E-18,.5193E-18,.5224E-18,.5167E-18,.5224E-18,
     A .5314E-18,.5478E-18,.5506E-18,.5487E-18,.5374E-18,.5361E-18,
     A .5314E-18,.5412E-18,.5536E-18,.5727E-18,.5697E-18,.5591E-18,
     A .5413E-18,.5390E-18,.5384E-18,.5543E-18,.5683E-18,.5815E-18,
     A .5757E-18,.5742E-18,.5652E-18,.5620E-18,.5549E-18,.5601E-18,
     A .5662E-18,.5787E-18,.5793E-18,.5809E-18,.5731E-18,.5717E-18,
     A .5696E-18,.5790E-18,.5812E-18,.5858E-18,.5864E-18,.5929E-18,
     A .5897E-18,.5931E-18,.5909E-18,.6003E-18,.6008E-18,.5956E-18,
     A .5794E-18,.5725E-18,.5665E-18,.5738E-18,.5830E-18,.5961E-18,
     A .5928E-18,.5857E-18,.5786E-18,.5922E-18,.6069E-18,.6094E-18,
     A .5995E-18,.6099E-18,.6199E-18,.6242E-18,.6135E-18,.6005E-18,
     A .5800E-18,.5706E-18,.5680E-18,.5803E-18,.5889E-18,.6001E-18,
     A .5951E-18,.5812E-18,.5609E-18,.5524E-18,.5581E-18,.5881E-18,
     A .6091E-18,.6159E-18,.6089E-18,.6110E-18,.6083E-18,.5958E-18,
     A .5749E-18,.5724E-18,.5863E-18,.6136E-18,.6172E-18,.5997E-18,
     A .5786E-18,.5632E-18,.5422E-18,.5333E-18,.5334E-18,.5452E-18,
     A .5535E-18,.5627E-18,.5644E-18,.5748E-18,.5818E-18,.5872E-18,
     A .5853E-18,.5857E-18,.5858E-18,.5878E-18,.5829E-18,.5872E-18,
     A .5898E-18,.5891E-18,.5698E-18,.5429E-18,.5145E-18,.5070E-18,
     A .5201E-18,.5460E-18,.5533E-18,.5427E-18,.5167E-18,.5022E-18,
     A .5102E-18,.5356E-18,.5462E-18,.5366E-18,.5123E-18,.5030E-18,
     A .5127E-18,.5340E-18,.5475E-18,.5669E-18,.5756E-18,.5647E-18,
     A .5367E-18,.5069E-18,.4751E-18,.4627E-18,.4775E-18,.5197E-18,
     A .5446E-18,.5355E-18,.4979E-18,.4633E-18,.4362E-18,.4249E-18,
     A .4261E-18,.4381E-18,.4513E-18,.4760E-18,.5040E-18,.5190E-18,
     A .5045E-18,.4818E-18,.4608E-18,.4677E-18,.4991E-18,.5282E-18,
     A .5153E-18,.4824E-18,.4520E-18,.4426E-18,.4413E-18,.4406E-18,
     A .4388E-18,.4448E-18,.4437E-18,.4374E-18,.4228E-18,.4173E-18,
     A .4200E-18,.4206E-18,.4132E-18,.4041E-18,.3991E-18,.4159E-18,
     A .4351E-18,.4400E-18,.4269E-18,.4201E-18,.4234E-18,.4322E-18,
     A .4261E-18,.4125E-18,.4013E-18,.4131E-18,.4364E-18,.4468E-18,
     A .4307E-18,.4092E-18,.3795E-18,.3540E-18,.3447E-18,.3620E-18,
     A .3790E-18,.3838E-18,.3726E-18,.3570E-18,.3403E-18,.3324E-18,
     A .3290E-18,.3303E-18,.3277E-18,.3327E-18,.3359E-18,.3442E-18,
     A .3614E-18,.3898E-18,.4038E-18,.4048E-18,.3893E-18,.3692E-18,
     A .3484E-18,.3353E-18,.3257E-18,.3221E-18,.3197E-18,.3271E-18,
     A .3393E-18,.3457E-18,.3380E-18,.3310E-18,.3221E-18,.3075E-18,
     A .2873E-18,.2706E-18,.2588E-18,.2586E-18,.2626E-18,.2688E-18,
     A .2705E-18,.2707E-18,.2662E-18,.2665E-18,.2760E-18,.3004E-18,
     A .3267E-18,.3359E-18,.3227E-18,.3016E-18,.2768E-18,.2590E-18,
     A .2514E-18,.2516E-18,.2538E-18,.2635E-18,.2746E-18,.2859E-18,
     A .2895E-18,.2883E-18,.2848E-18,.2798E-18,.2648E-18,.2448E-18,
     A .2243E-18,.2088E-18,.1954E-18,.1861E-18,.1809E-18,.1821E-18,
     A .1867E-18,.1946E-18,.2063E-18,.2215E-18,.2279E-18,.2268E-18,
     A .2318E-18,.2455E-18,.2486E-18,.2384E-18,.2255E-18,.2160E-18,
     A .2064E-18,.2018E-18,.2034E-18,.2117E-18,.2205E-18,.2277E-18,
     A .2275E-18,.2242E-18,.2190E-18,.2154E-18,.2118E-18,.2051E-18,
     A .1922E-18,.1811E-18,.1719E-18,.1643E-18,.1575E-18,.1536E-18,
     A .1490E-18,.1450E-18,.1439E-18,.1504E-18,.1600E-18,.1660E-18,
     A .1673E-18,.1670E-18,.1626E-18,.1560E-18,.1493E-18,.1465E-18,
     A .1456E-18,.1468E-18,.1469E-18,.1517E-18,.1631E-18,.1763E-18,
     A .1822E-18,.1826E-18,.1770E-18,.1689E-18,.1629E-18,.1588E-18,
     A .1522E-18,.1467E-18,.1426E-18,.1428E-18,.1418E-18,.1411E-18,
     A .1403E-18,.1402E-18,.1372E-18,.1314E-18,.1237E-18,.1168E-18,
     A .1081E-18,.1002E-18,.9535E-19,.9400E-19,.9250E-19,.9364E-19,
     A .9891E-19,.1062E-18,.1090E-18,.1072E-18,.1007E-18,.9574E-19,
     A .9610E-19,.1035E-18,.1125E-18,.1197E-18,.1230E-18,.1227E-18,
     A .1195E-18,.1196E-18,.1215E-18,.1194E-18,.1129E-18,.1089E-18,
     A .1097E-18,.1133E-18,.1140E-18,.1139E-18,.1123E-18,.1081E-18,
     A .1014E-18,.9592E-19,.9093E-19,.8928E-19,.9012E-19,.9309E-19,
     A .9346E-19,.9064E-19,.8508E-19,.7872E-19,.7264E-19,.6851E-19,
     A .6581E-19,.6528E-19,.6448E-19,.6436E-19,.6304E-19,.6064E-19,
     A .5848E-19,.5868E-19,.6007E-19,.6240E-19,.6438E-19,.6707E-19,
     A .7045E-19,.7660E-19,.8343E-19,.8884E-19,.9050E-19,.8780E-19,
     A .8041E-19,.7178E-19,.6550E-19,.6453E-19,.6664E-19,.6981E-19,
     A .7006E-19,.6959E-19,.6998E-19,.7328E-19,.7694E-19,.7922E-19,
     A .7851E-19,.7591E-19,.7018E-19,.6453E-19,.5861E-19,.5399E-19,
     A .4973E-19,.4637E-19,.4234E-19,.3960E-19,.3764E-19,.3778E-19,
     A .3936E-19,.4266E-19,.4567E-19,.4760E-19,.4722E-19,.4568E-19,
     A .4252E-19,.3999E-19,.3776E-19,.3691E-19,.3693E-19,.3799E-19,
     A .3991E-19,.4316E-19,.4613E-19,.4810E-19,.4804E-19,.4728E-19,
     A .4547E-19,.4365E-19,.4129E-19,.4084E-19,.4186E-19,.4474E-19,
     A .4708E-19,.4937E-19,.5042E-19,.5343E-19,.5558E-19,.5637E-19,
     A .5360E-19,.4890E-19,.4300E-19,.3848E-19,.3421E-19,.3249E-19,
     A .3136E-19,.3107E-19,.3089E-19,.3137E-19,.3211E-19,.3354E-19,
     A .3472E-19,.3630E-19,.3648E-19,.3653E-19,.3544E-19,.3409E-19,
     A .3110E-19,.2777E-19,.2416E-19,.2212E-19,.2032E-19,.2016E-19,
     A .1989E-19,.2009E-19,.1990E-19,.2029E-19,.2050E-19,.2091E-19,
     A .2102E-19,.2139E-19,.2113E-19,.2185E-19,.2327E-19,.2687E-19,
     A .3025E-19,.3433E-19,.3509E-19,.3329E-19,.3093E-19,.3006E-19,
     A .2844E-19,.2882E-19,.2981E-19,.3168E-19,.3156E-19,.3037E-19,
     A .2704E-19,.2535E-19,.2430E-19,.2540E-19,.2559E-19,.2545E-19,
     A .2337E-19,.2319E-19,.2153E-19,.2240E-19,.2065E-19,.2123E-19,
     A .1851E-19,.2083E-19,.1809E-19,.1948E-19,.1759E-19,.1951E-19,
     A .1765E-19,.1853E-19,.1513E-19,.1739E-19,.1593E-19,.1785E-19,
     A .1512E-19,.1705E-19,.1485E-19,.1655E-19,.1372E-19,.1616E-19,
     A .1367E-19,.1604E-19,.1288E-19,.1596E-19,.1288E-19,.1559E-19,
     A .1189E-19,.1527E-19,.1276E-19,.1512E-19,.1188E-19,.1527E-19,
     A .1218E-19,.1502E-19,.1151E-19,.1483E-19,.1127E-19,.1474E-19,
     A .1106E-19,.1463E-19,.1180E-19,.1466E-19,.1176E-19,.1454E-19,
     A .1080E-19,.1419E-19,.1087E-19,.1470E-19,.1151E-19,.1405E-19,
     A .1120E-19,.1446E-19/
      REAL yg(kw)
      REAL a1, a2, dum
      INTEGER ierr
      INTEGER i, l, n, idum
      CHARACTER*40 fil
*_______________________________________________________________________

************* absorption cross sections:
*     measurements of Davidson et al. (198x) at 273K
*     from 263.8 to 648.8 nm in approximately 0.5 nm intervals

      fil = 'DATAE1/NO2/NO2_ncar_00.abs'
      n = nno2

      CALL addpnt(x1,y1,kdata,n,x1(1)*(1.-deltax),0.)
      CALL addpnt(x1,y1,kdata,n,          0.,0.)
      CALL addpnt(x1,y1,kdata,n,x1(n)*(1.+deltax),0.)
      CALL addpnt(x1,y1,kdata,n,      1.e+38,0.)
      CALL inter2(nw,wl,yg,n,x1,y1,ierr)
      IF (ierr .NE. 0) THEN
         WRITE(*,*) ierr, fil
         STOP
      ENDIF

      DO 13, l = 1, nw-1
         xsno2(l) = yg(l)
   13 CONTINUE

*_______________________________________________________________________

      RETURN
      END
      SUBROUTINE rdo3xs(nw,wl,xso3,s226,s263,s298)

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Read ozone molecular absorption cross section.  Re-grid data to match    =*
*=  specified wavelength working grid.                                       =*
*-----------------------------------------------------------------------------*
*=  PARAMETERS:                                                              =*
*=  NW     - INTEGER, number of specified intervals + 1 in working        (I)=*
*=           wavelength grid                                                 =*
*=  WL     - REAL, vector of lower limits of wavelength intervals in      (I)=*
*=           working wavelength grid                                         =*
*=  XSO3   - REAL, molecular absoprtion cross section (cm^2) of O3 at     (O)=*
*=           each specified wavelength (WMO value at 273)                    =*
*=  S226   - REAL, molecular absoprtion cross section (cm^2) of O3 at     (O)=*
*=           each specified wavelength (value from Molina and Molina at 226K)=*
*=  S263   - REAL, molecular absoprtion cross section (cm^2) of O3 at     (O)=*
*=           each specified wavelength (value from Molina and Molina at 263K)=*
*=  S298   - REAL, molecular absoprtion cross section (cm^2) of O3 at     (O)=*
*=           each specified wavelength (value from Molina and Molina at 298K)=*
*-----------------------------------------------------------------------------*
*=  EDIT HISTORY:                                                            =*
*=  02/97  Changed offset for grid-end interpolation to relative number      =*
*=         (x * (1 +- deltax))                                               =*
*-----------------------------------------------------------------------------*
*= This program is free software;  you can redistribute it and/or modify     =*
*= it under the terms of the GNU General Public License as published by the  =*
*= Free Software Foundation;  either version 2 of the license, or (at your   =*
*= option) any later version.                                                =*
*= The TUV package is distributed in the hope that it will be useful, but    =*
*= WITHOUT ANY WARRANTY;  without even the implied warranty of MERCHANTIBI-  =*
*= LITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public     =*
*= License for more details.                                                 =*
*= To obtain a copy of the GNU General Public License, write to:             =*
*= Free Software Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.   =*
*-----------------------------------------------------------------------------*
*= To contact the authors, please mail to:                                   =*
*= Sasha Madronich, NCAR/ACD, P.O.Box 3000, Boulder, CO, 80307-3000, USA  or =*
*= send email to:  sasha@ucar.edu                                            =*
*-----------------------------------------------------------------------------*
*= Copyright (C) 1994,95,96  University Corporation for Atmospheric Research =*
*-----------------------------------------------------------------------------*

      IMPLICIT NONE
      INCLUDE 'tuv.params'

      INTEGER kdata
      PARAMETER(kdata=1000)

* input: (altitude working grid)
      INTEGER nw
      REAL wl(kw)

* output:
* ozone absorption cross section at three different 
* temperatures: 226, 263, 298 Kelvin.  Can interpolate
* to different temperatures. Units are cm2 molecule-1

      REAL xso3(kw), s226(kw),s263(kw),s298(kw)

* local:
      integer nwmo,ns
      parameter(nwmo=158,ns=220) 
      REAL x1(kdata),x2(kdata),x3(kdata),xs1wmo(kdata),ys1wmo(kdata)
      REAL xs1(kdata),xs2(kdata),xs3(kdata),ys1(kdata),
     1 ys2(kdata),ys3(kdata)       
      
      data xs1wmo(1:nwmo)/176.2,177.8,179.4,181.0,182.7,184.3,185.7,
     1 187.4,189.6,191.4,193.2,195.1,197.0,199.0,
     2 201.0,203.1,205.1,207.3,209.4,211.6,213.9,
     3 216.2,218.6,221.0,223.5,226.0,228.6,231.2,
     4 233.9,236.7,239.5,242.4,245.4,248.5,251.6,
     5 254.8,258.1,261.4,264.9,268.5,272.1,275.9,
     6 279.7,283.7,287.8,292.0,296.3,300.8,305.4,
     7 310.1,315.0,320.0,325.0,330.0,335.0,340.0,
     8 345.0,350.0,355.0,360.0,365.0,370.0,375.0,
     9 380.0,385.0,390.0,395.0,400.0,405.0,410.0,
     a 415.0,420.0,425.0,430.0,435.0,440.0,445.0,
     b 450.0,455.0,460.0,465.0,470.0,475.0,480.0,
     c 485.0,490.0,495.0,500.0,505.0,510.0,515.0,
     d 520.0,525.0,530.0,535.0,540.0,545.0,550.0,
     e 555.0,560.0,565.0,570.0,575.0,580.0,585.0,
     f 590.0,595.0,600.0,605.0,610.0,615.0,620.0,
     g 625.0,630.0,635.0,640.0,645.0,650.0,655.0,
     h 660.0,665.0,670.0,675.0,680.0,685.0,690.0,
     i 695.0,700.0,705.0,710.0,715.0,720.0,725.0,
     j 730.0,735.0,740.0,745.0,750.0,755.0,760.0,
     k 765.0,770.0,775.0,780.0,785.0,790.0,795.0,
     l 800.0,805.0,810.0,815.0,820.0,825.0,830.0,
     m 835.0,840.0,845.0,850.0/
       data ys1wmo(1:nwmo)/.811E-18,.799E-18,.786E-18,.763E-18,.729E-18,
     1  .688E-18,.640E-18,.588E-18,.531E-18,.480E-18,.438E-18,.411E-18,
     2 .369E-18,.330E-18,.326E-18,.326E-18,.351E-18,.411E-18,.484E-18,
     3 .626E-18,.857E-18,.117E-17,.152E-17,.197E-17,.255E-17,.324E-17,
     4 .400E-17,.483E-17,.579E-17,.686E-17,.797E-17,.900E-17,.100E-16,
     5 .108E-16,.113E-16,.115E-16,.112E-16,.106E-16,.965E-17,.834E-17,
     6 .692E-17,.542E-17,.402E-17,.277E-17,.179E-17,.109E-17,.624E-18,
     7 .343E-18,.185E-18,.980E-19,.501E-19,.249E-19,.120E-19,.617E-20,
     8 .274E-20,.117E-20,.588E-21,.266E-21,.109E-21,.549E-22,.000E+00,
     9 .000E+00,.000E+00,.000E+00,.000E+00,.000E+00,.000E+00,.000E+00,
     a .000E+00,.291E-22,.314E-22,.399E-22,.654E-22,.683E-22,.866E-22,
     b .125E-21,.149E-21,.171E-21,.212E-21,.357E-21,.368E-21,.406E-21,
     c .489E-21,.711E-21,.843E-21,.828E-21,.909E-21,.122E-20,.162E-20,
     d .158E-20,.160E-20,.178E-20,.207E-20,.255E-20,.274E-20,.288E-20,
     e .307E-20,.317E-20,.336E-20,.388E-20,.431E-20,.467E-20,.475E-20,
     f .455E-20,.435E-20,.442E-20,.461E-20,.489E-20,.484E-20,.454E-20,
     g .424E-20,.390E-20,.360E-20,.343E-20,.317E-20,.274E-20,.261E-20,
     h .242E-20,.220E-20,.202E-20,.185E-20,.167E-20,.154E-20,.142E-20,
     i .125E-20,.112E-20,.102E-20,.920E-21,.840E-21,.770E-21,.690E-21,
     j .630E-21,.570E-21,.525E-21,.475E-21,.447E-21,.420E-21,.375E-21,
     k .325E-21,.292E-21,.276E-21,.270E-21,.280E-21,.285E-21,.252E-21,
     l .220E-21,.182E-21,.163E-21,.175E-21,.190E-21,.185E-21,.170E-21,
     m .152E-21,.142E-21,.140E-21,.140E-21,.142E-21,.145E-21/
       data xs1(1:ns)/240.5,241.0,241.5,242.0,242.5,243.0,243.5,
     1  244.0,244.5,245.0,245.5,246.0,246.5,247.0,
     2  247.5,248.0,248.5,249.0,249.5,250.0,250.5,
     3  251.0,251.5,252.0,252.5,253.0,253.5,254.0,
     4  254.5,255.0,255.5,256.0,256.5,257.0,257.5,
     5  258.0,258.5,259.0,259.5,260.0,260.5,261.0,
     6  261.5,262.0,262.5,263.0,263.5,264.0,264.5,
     7  265.0,265.5,266.0,266.5,267.0,267.5,268.0,
     8  268.5,269.0,269.5,270.0,270.5,271.0,271.5,
     9  272.0,272.5,273.0,273.5,274.0,274.5,275.0,
     a  275.5,276.0,276.5,277.0,277.5,278.0,278.5,
     b  279.0,279.5,280.0,280.5,281.0,281.5,282.0,
     c  282.5,283.0,283.5,284.0,284.5,285.0,285.5,
     d  286.0,286.5,287.0,287.5,288.0,288.5,289.0,
     e  289.5,290.0,290.5,291.0,291.5,292.0,292.5,
     f  293.0,293.5,294.0,294.5,295.0,295.5,296.0,
     g  296.5,297.0,297.5,298.0,298.5,299.0,299.5,
     h  300.0,300.5,301.0,301.5,302.0,302.5,303.0,
     i  303.5,304.0,304.5,305.0,305.5,306.0,306.5,
     j  307.0,307.5,308.0,308.5,309.0,309.5,310.0,
     k  310.5,311.0,311.5,312.0,312.5,313.0,313.5,
     l  314.0,314.5,315.0,315.5,316.0,316.5,317.0,
     m  317.5,318.0,318.5,319.0,319.5,320.0,320.5,
     n  321.0,321.5,322.0,322.5,323.0,323.5,324.0,
     o  324.5,325.0,325.5,326.0,326.5,327.0,327.5,
     p  328.0,328.5,329.0,329.5,330.0,330.5,331.0,
     q  331.5,332.0,332.5,333.0,333.5,334.0,334.5,
     r  335.0,335.5,336.0,336.5,337.0,337.5,338.0,
     s  338.5,339.0,339.5,340.0,340.5,341.0,341.5,
     t  342.0,342.5,343.0,343.5,344.0,344.5,345.0,
     u  345.5,346.0,346.5,347.0,347.5,348.0,348.5,
     v  349.0,349.5,350.0/
       data ys1(1:ns)/ 842.4, 864.0, 891.8, 901.6, 916.3, 939.6, 962.8,
     1 975.2, 977.6,1007.0,1040.0,1042.0,1039.0,1058.0,
     2 1080.0,1079.0,1094.0,1124.0,1121.0,1134.0,1118.0,
     3 1123.0,1164.0,1165.0,1143.0,1149.0,1167.0,1169.0,
     4 1145.0,1174.0,1184.0,1158.0,1132.0,1147.0,1134.0,
     5 1130.0,1172.0,1151.0,1084.0,1086.0,1128.0,1095.0,
     6 1056.0,1064.0,1072.0,1016.0, 994.0,1013.0,1013.0,
     7 961.2, 968.3, 946.7, 920.6, 877.9, 900.5, 872.2,
     8 829.3, 803.8, 820.8, 796.1, 754.4, 736.7, 724.8,
     9 710.7, 684.2, 666.0, 640.1, 601.3, 592.7, 587.1,
     a 562.9, 537.9, 519.6, 504.0, 492.3, 461.7, 445.7,
     b 424.5, 412.9, 398.3, 371.8, 364.2, 333.7, 322.0,
     c 310.0, 299.0, 293.5, 264.4, 252.4, 239.8, 225.4,
     d 219.5, 205.3, 197.8, 184.2, 168.7, 162.0, 153.5,
     e 142.9, 136.5, 129.9, 122.5, 111.8, 104.3,  99.6900,
     f  93.7500, 86.3700, 81.4800, 78.2400, 72.7, 65.9800, 62.0200,
     g  58.6500, 54.6900, 50.0200, 47.1400, 44.4900, 42.3100, 38.0,
     h  36.1600, 33.9600, 31.2800,  28.9500, 27.9600, 26.0, 23.35,
     i  21.9800, 21.5600, 19.8400,  17.7100, 16.3100, 16.0300, 15.49,
     j  13.6700,  12.0,  12.0200,  11.1800,  10.6400,9.4250,8.6370,
     k   8.3140,7.9250,7.8310,6.6970,6.0740,5.6910,6.1690,
     l   5.3340,4.2940,4.1860,4.8370,3.8960,3.3510,3.3510,
     m   3.4200,3.0630,2.3400,1.9990,2.7120,2.8590,1.9630,
     n   1.3680,1.3610,2.1170,1.7930,1.5290,0.8902,0.7852,
     o   0.9488,1.4860,1.2410,0.7276,0.4945,0.6158,0.6491,
     p   1.1580,0.6374,0.3460,0.2666,0.2854,0.5695,0.7130,
     q   0.4584,0.2415,0.1567,0.2187,0.3693,0.3788,0.2117,
     r   0.1274,0.0876,0.0871,0.1837,0.2464,0.2590,0.1311,
     s   0.0695,0.0546,0.0713,0.1549,0.1133,0.0594,0.0352,
     t   0.0285,0.0236,0.0372,0.0441,0.1271,0.0730,0.0390,
     u   0.0269,0.0195,0.0178,0.0283,0.0295,0.0181,0.0096,
     v   0.0083,0.0100,0.0076/
       data ys2(1:ns)/ 837.0, 858.6, 886.0, 895.1, 911.1, 932.0, 957.3,
     1  968.4, 970.7, 994.8,1028.0,1032.0,1035.0,1051.0,
     2  1069.0,1072.0,1082.0,1110.0,1110.0,1121.0,1110.0,
     3  1115.0,1153.0,1153.0,1137.0,1142.0,1155.0,1159.0,
     4  1136.0,1167.0,1176.0,1151.0,1125.0,1139.0,1126.0,
     5  1122.0,1164.0,1143.0,1083.0,1085.0,1120.0,1090.0,
     6  1054.0,1061.0,1067.0,1015.0, 994.7,1010.0,1012.0,
     7  960.5, 965.6, 943.4, 920.0, 878.0, 895.1, 871.1,
     8  830.1, 803.7, 822.0, 798.2, 755.5, 736.5, 726.6,
     9  711.4, 687.6, 675.4, 642.1, 609.9, 597.5, 586.1,
     a  561.0, 539.5, 523.3, 508.9, 490.3, 465.0, 452.1,
     b  429.8, 413.8, 395.2, 374.7, 367.7, 337.6, 323.9,
     c  311.1, 301.1, 283.3, 267.4, 257.7, 244.5, 227.6,
     d  220.5, 209.4, 201.0, 185.8, 172.6, 164.1, 155.2,
     e  147.0, 139.6, 131.9, 124.9, 114.5, 107.4, 102.8,
     f  96.59,  88.31, 84.47, 79.43, 74.92, 68.53, 64.66,
     g  60.1, 56.96,  52.3100,  48.7500,  45.9700,  43.4200,  39.59,
     h  37.4, 35.4,  32.2700,  30.1800,  28.9100,  26.7800,  24.52,
     i  23.0, 22.26,  20.8100,  18.6200,  17.3200,  16.7800,  16.12,
     j  14.4900,  13.0600,  12.6700,  11.7300,  11.4500,9.7150,9.3030,
     k   8.8030,8.4800,8.2430,7.2,6.5750,6.1990,6.4070,
     l   5.6700,4.7760,4.5850,5.0280,4.2260,3.7500,3.7400,
     m   3.6620,3.3220,2.6510,2.3500,2.8800,3.0220,2.2230,
     n   1.6700,1.6150,2.2420,1.9750,1.7600,1.1410,0.9893,
     o   1.1040,1.5930,1.3590,0.9125,0.6579,0.7214,0.7691,
     p   1.2210,0.7731,0.4785,0.3750,0.3770,0.6381,0.7511,
     q   0.5535,0.3283,0.2295,0.2731,0.4033,0.4474,0.2788,
     r   0.1730,0.1347,0.1312,0.2216,0.2686,0.2978,0.1641,
     s   0.1002,0.0836,0.1127,0.1670,0.1434,0.0863,0.0537,
     t   0.0488,0.0414,0.0640,0.0633,0.1336,0.0865,0.0590,
     u   0.0439,0.0384,0.0271,0.0323,0.0429,0.0249,0.0217,
     v   0.0163,0.0251,0.0193/
       data ys3(1:ns)/ 839.8, 860.3, 886.1, 897.1, 912.6, 933.3, 956.0,
     1  971.7, 975.1, 993.2,1025.0,1033.0,1034.0,1047.0,
     2  1068.0,1071.0,1082.0,1112.0,1112.0,1124.0,1113.0,
     3  1114.0,1146.0,1155.0,1139.0,1140.0,1155.0,1159.0,
     4  1140.0,1161.0,1173.0,1154.0,1130.0,1139.0,1129.0,
     5  1124.0,1157.0,1145.0,1090.0,1080.0,1114.0,1094.0,
     6  1057.0,1057.0,1066.0,1022.0, 995.9,1006.0,1007.0,
     7  965.7, 965.0, 948.5, 922.1, 884.1, 893.8, 875.4,
     8  836.3, 810.4, 815.1, 798.0, 763.6, 741.5, 726.6,
     9  714.7, 689.9, 669.8, 646.0, 614.0, 596.8, 591.3,
     a  560.6, 545.0, 522.1, 509.6, 490.6, 466.8, 451.6,
     b  432.6, 415.8, 400.1, 372.7, 367.3, 340.4, 325.0,
     c  312.2, 302.5, 296.3, 271.2, 259.8, 246.5, 232.2,
     d  223.8, 210.3, 203.4, 191.2, 175.0, 167.3, 158.5,
     e  151.2, 141.8, 135.8, 128.5, 119.4, 111.1, 105.6,
     f  100.2,  92.1500,  87.1100,  81.9300,  77.53,  71.55,  67.27,
     g  62.7400,  59.5500,  55.2,  51.2400,  48.4,  45.51,  41.87,
     h  39.6400,  37.2300,  34.6300,  32.2300,  30.73,  28.69,  26.5,
     i  24.5800,  24.0100,  21.8800,  20.15,  18.79,  18.08,  17.11,
     j  15.6500,  14.2100,  13.6400,  12.9,  12.43,  10.97,  10.2,
     k   9.7050,9.2600,8.9900,7.9470,7.2090,6.8830,7.0320,
     l   6.2940,5.4,5.1990,5.5290,4.7920,4.2510,4.1460,
     m   4.0980,3.7570,3.1270,2.7650,3.2130,3.2430,2.5610,
     n   2.0410,1.9220,2.4350,2.2,1.9830,1.4210,1.2500,
     o   1.3210,1.7270,1.5130,1.1050,0.8550,0.8875,0.9160,
     p   1.3,0.9452,0.6504,0.5180,0.4923,0.7302,0.8328,
     q   0.6572,0.4347,0.3194,0.3528,0.4665,0.5343,0.2550,
     r   0.2434,0.1955,0.1875,0.2753,0.3236,0.3299,0.2086,
     s   0.1413,0.1343,0.1662,0.2082,0.1784,0.1134,0.0808,
     t   0.0776,0.0714,0.0977,0.0936,0.1419,0.0916,0.0644,
     u   0.0698,0.0641,0.0539,0.0467,0.0503,0.0386,0.0285,
     v   0.0271,0.0405,0.0294/
      REAL yg(kw)
      REAL a1, a2, dum
      INTEGER ierr
      INTEGER i, iw, n, idum, n1, n2, n3

      character*40 fil

*_______________________________________________________________________


************ from WMO 1985 Ozone Assessment
* from 175.439 to 847.500 nm
* use value at 273 K

      n = nwmo

      CALL addpnt(xs1wmo,ys1wmo,kdata,n,xs1wmo(1)*(1.-deltax),0.)
      CALL addpnt(xs1wmo,ys1wmo,kdata,n,          0.,0.)
      CALL addpnt(xs1wmo,ys1wmo,kdata,n,xs1wmo(n)*(1.+deltax),0.)
      CALL addpnt(xs1wmo,ys1wmo,kdata,n,      1.e+38,0.)
      CALL inter2(nw,wl,yg,n,xs1wmo,ys1wmo,ierr)
      IF (ierr .NE. 0) THEN
         WRITE(*,*) ierr, fil
         STOP
      ENDIF
            
      DO 13, iw = 1, nw-1
         xso3(iw) = yg(iw)
   13 CONTINUE

************* ozone absorption cross sections:
* For Hartley aand Huggins bands, use temperature-dependent values from
* Molina, L. T., and M. J. Molina, Absolute absorption cross sections
* of ozone in the 185- to 350-nm wavelength range,
* J. Geophys. Res., vol. 91, 14501-14508, 1986.

      n1 = ns
      n2 = ns
      n3 = ns
      DO  i = 1, n1
         xs2(i) = xs1(i)
         xs3(i) = xs1(i)
      ENDDO

      CALL addpnt(xs1,ys1,kdata,n1,xs1(1)*(1.-deltax),0.)
      CALL addpnt(xs1,ys1,kdata,n1,               0.,0.)
      CALL addpnt(xs1,ys1,kdata,n1,xs1(n1)*(1.+deltax),0.)
      CALL addpnt(xs1,ys1,kdata,n1,            1.e+38,0.)
      CALL inter2(nw,wl,yg,n1,xs1,ys1,ierr)
      IF (ierr .NE. 0) THEN
         WRITE(*,*) ierr, 'rdo3xs'
         STOP
      ENDIF
      DO iw = 1, nw-1
         s226(iw) = yg(iw)*1.E-20
      ENDDO

      CALL addpnt(xs2,ys2,kdata,n2,xs2(1)*(1.-deltax),0.)
      CALL addpnt(xs2,ys2,kdata,n2,               0.,0.)
      CALL addpnt(xs2,ys2,kdata,n2,xs2(n2)*(1.+deltax),0.)
      CALL addpnt(xs2,ys2,kdata,n2,            1.e+38,0.)
      CALL inter2(nw,wl,yg,n2,xs2,ys2,ierr)
      IF (ierr .NE. 0) THEN
         WRITE(*,*) ierr, 'rdo3xs'
         STOP
      ENDIF
      DO iw = 1, nw-1
         s263(iw) = yg(iw)*1.E-20
      ENDDO

      CALL addpnt(xs3,ys3,kdata,n3,xs3(1)*(1.-deltax),0.)
      CALL addpnt(xs3,ys3,kdata,n3,               0.,0.)
      CALL addpnt(xs3,ys3,kdata,n3,xs3(n3)*(1.+deltax),0.)
      CALL addpnt(xs3,ys3,kdata,n3,            1.e+38,0.)
      CALL inter2(nw,wl,yg,n3,xs3,ys3,ierr)
      IF (ierr .NE. 0) THEN
         WRITE(*,*) ierr, 'rdo3xs'
         STOP
      ENDIF
      DO iw = 1, nw-1
         s298(iw) = yg(iw)*1.E-20
      ENDDO

*_______________________________________________________________________

      RETURN
      END
      SUBROUTINE rdso2xs(nw,wl,xsso2)

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Read SO2 molecular absorption cross section.  Re-grid data to match      =*
*=  specified wavelength working grid.                                       =*
*-----------------------------------------------------------------------------*
*=  PARAMETERS:                                                              =*
*=  NW     - INTEGER, number of specified intervals + 1 in working        (I)=*
*=           wavelength grid                                                 =*
*=  WL     - REAL, vector of lower limits of wavelength intervals in      (I)=*
*=           working wavelength grid                                         =*
*=  XSSO2  - REAL, molecular absoprtion cross section (cm^2) of SO2 at    (O)=*
*=           each specified wavelength                                       =*
*-----------------------------------------------------------------------------*
*=  EDIT HISTORY:                                                            =*
*=  02/97  Changed offset for grid-end interpolation to relative number      =*
*=         (x * (1 +- deltax)                                                =*
*-----------------------------------------------------------------------------*
*= This program is free software;  you can redistribute it and/or modify     =*
*= it under the terms of the GNU General Public License as published by the  =*
*= Free Software Foundation;  either version 2 of the license, or (at your   =*
*= option) any later version.                                                =*
*= The TUV package is distributed in the hope that it will be useful, but    =*
*= WITHOUT ANY WARRANTY;  without even the implied warranty of MERCHANTIBI-  =*
*= LITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public     =*
*= License for more details.                                                 =*
*= To obtain a copy of the GNU General Public License, write to:             =*
*= Free Software Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.   =*
*-----------------------------------------------------------------------------*
*= To contact the authors, please mail to:                                   =*
*= Sasha Madronich, NCAR/ACD, P.O.Box 3000, Boulder, CO, 80307-3000, USA  or =*
*= send email to:  sasha@ucar.edu                                            =*
*-----------------------------------------------------------------------------*
*= Copyright (C) 1994,95,96  University Corporation for Atmospheric Research =*
*-----------------------------------------------------------------------------*

      IMPLICIT NONE
      INCLUDE 'tuv.params'

      INTEGER kdata
      PARAMETER(kdata=1000)

* input: (altitude working grid)
      INTEGER nw
      REAL wl(kw)

* output:

      REAL xsso2(kw)

* local:
      integer nso2
      parameter(nso2=704)
      REAL x1(kdata)
      REAL y1(kdata),y2(kdata),y3(kdata)
      data x1(1:nso2)/280.00,281.00,281.90,282.50,283.10,
     A 283.90,284.90,286.10,287.00,288.30,289.00,290.20,290.50,290.80,
     A 291.50,292.20,293.20,294.00,295.40,296.30,297.80,298.10,299.20,
     A 300.00,300.03,300.06,300.09,300.12,300.15,300.18,300.21,300.24,
     A 300.27,300.30,300.33,300.36,300.39,300.42,300.45,300.48,300.51,
     A 300.54,300.57,300.60,300.63,300.66,300.69,300.72,300.75,300.78,
     A 300.81,300.84,300.87,300.90,300.93,300.96,300.99,301.02,301.05,
     A 301.08,301.11,301.14,301.17,301.20,301.23,301.26,301.29,301.32,
     A 301.35,301.38,301.41,301.44,301.47,301.50,301.53,301.56,301.59,
     A 301.62,301.65,301.68,301.71,301.74,301.77,301.80,301.83,301.86,
     A 301.89,301.92,301.95,301.98,302.01,302.04,302.07,302.10,302.13,
     A 302.16,302.19,302.22,302.25,302.28,302.31,302.34,302.37,302.40,
     A 302.43,302.46,302.49,302.52,302.55,302.58,302.61,302.64,302.67,
     A 302.70,302.73,302.76,302.79,302.82,302.85,302.88,302.91,302.94,
     A 302.97,303.00,303.03,303.06,303.09,303.12,303.15,303.18,303.21,
     A 303.24,303.27,303.30,303.33,303.36,303.39,303.42,303.45,303.48,
     A 303.51,303.54,303.57,303.60,303.63,303.66,303.69,303.72,303.75,
     A 303.78,303.81,303.84,303.87,303.90,303.93,303.96,303.99,304.02,
     A 304.05,304.08,304.11,304.14,304.17,304.20,304.23,304.26,304.29,
     A 304.32,304.35,304.38,304.41,304.44,304.47,304.50,304.53,304.56,
     A 304.59,304.62,304.65,304.68,304.71,304.74,304.77,304.80,304.83,
     A 304.86,304.89,304.92,304.95,304.98,305.01,305.04,305.07,305.10,
     A 305.13,305.16,305.19,305.22,305.25,305.28,305.31,305.34,305.37,
     A 305.40,305.43,305.46,305.49,305.52,305.55,305.58,305.61,305.64,
     A 305.67,305.70,305.73,305.76,305.79,305.82,305.85,305.88,305.91,
     A 305.94,305.97,306.00,306.03,306.06,306.09,306.12,306.15,306.18,
     A 306.21,306.24,306.27,306.30,306.33,306.36,306.39,306.42,306.45,
     A 306.48,306.51,306.54,306.57,306.60,306.63,306.66,306.69,306.72,
     A 306.75,306.78,306.81,306.84,306.87,306.90,306.93,306.96,306.99,
     A 307.02,307.05,307.08,307.11,307.14,307.17,307.20,307.23,307.26,
     A 307.29,307.32,307.35,307.38,307.41,307.44,307.47,307.50,307.53,
     A 307.56,307.59,307.62,307.65,307.68,307.71,307.74,307.77,307.80,
     A 307.83,307.86,307.89,307.92,307.95,307.98,308.01,308.04,308.07,
     A 308.10,308.13,308.16,308.19,308.22,308.25,308.28,308.31,308.34,
     A 308.37,308.40,308.43,308.46,308.49,308.52,308.55,308.58,308.61,
     A 308.64,308.67,308.70,308.73,308.76,308.79,308.82,308.85,308.88,
     A 308.91,308.94,308.97,309.00,309.03,309.06,309.09,309.12,309.15,
     A 309.18,309.21,309.24,309.27,309.30,309.33,309.36,309.39,309.42,
     A 309.45,309.48,309.51,309.54,309.57,309.60,309.63,309.66,309.69,
     A 309.72,309.75,309.78,309.81,309.84,309.87,309.90,309.93,309.96,
     A 309.99,310.02,310.05,310.08,310.11,310.14,310.17,310.20,310.23,
     A 310.26,310.29,310.32,310.35,310.38,310.41,310.44,310.47,310.50,
     A 310.53,310.56,310.59,310.62,310.65,310.68,310.71,310.74,310.77,
     A 310.80,310.83,310.86,310.89,310.92,310.95,310.98,311.01,311.04,
     A 311.07,311.10,311.13,311.16,311.19,311.22,311.25,311.28,311.31,
     A 311.34,311.37,311.40,311.43,311.46,311.49,311.52,311.55,311.58,
     A 311.61,311.64,311.67,311.70,311.73,311.76,311.79,311.82,311.85,
     A 311.88,311.91,311.94,311.97,312.00,312.03,312.06,312.09,312.12,
     A 312.15,312.18,312.21,312.24,312.27,312.30,312.33,312.36,312.39,
     A 312.42,312.45,312.48,312.51,312.54,312.57,312.60,312.63,312.66,
     A 312.69,312.72,312.75,312.78,312.81,312.84,312.87,312.90,312.93,
     A 312.96,312.99,313.02,313.05,313.08,313.11,313.14,313.17,313.20,
     A 313.23,313.26,313.29,313.32,313.35,313.38,313.41,313.44,313.47,
     A 313.50,313.53,313.56,313.59,313.62,313.65,313.68,313.71,313.74,
     A 313.77,313.80,313.83,313.86,313.89,313.92,313.95,313.98,314.01,
     A 314.04,314.07,314.10,314.13,314.16,314.19,314.22,314.25,314.28,
     A 314.31,314.34,314.37,314.40,314.43,314.46,314.49,314.52,314.55,
     A 314.58,314.61,314.64,314.67,314.70,314.73,314.76,314.79,314.82,
     A 314.85,314.88,314.91,314.94,314.97,315.00,315.03,315.06,315.09,
     A 315.12,315.15,315.18,315.21,315.24,315.27,315.30,315.33,315.36,
     A 315.39,315.42,315.45,315.48,315.51,315.54,315.57,315.60,315.63,
     A 315.66,315.69,315.72,315.75,315.78,315.81,315.84,315.87,315.90,
     A 315.93,315.96,315.99,316.02,316.05,316.08,316.11,316.14,316.17,
     A 316.20,316.23,316.26,316.29,316.32,316.35,316.38,316.41,316.44,
     A 316.47,316.50,316.53,316.56,316.59,316.62,316.65,316.68,316.71,
     A 316.74,316.77,316.80,316.83,316.86,316.89,316.92,316.95,316.98,
     A 317.01,317.04,317.07,317.10,317.13,317.16,317.19,317.22,317.25,
     A 317.28,317.31,317.34,317.37,317.40,317.43,317.46,317.49,317.52,
     A 317.55,317.58,317.61,317.64,317.67,317.70,317.73,317.76,317.79,
     A 317.82,317.85,317.88,317.91,317.94,317.97,318.00,318.03,318.06,
     A 318.09,318.12,318.15,318.18,318.21,318.24,318.27,318.30,318.33,
     A 318.36,318.39,318.42,318.45,318.48,318.51,318.54,318.57,318.60,
     A 318.63,318.66,318.69,318.72,318.75,318.78,318.81,318.84,318.87,
     A 318.90,318.93,318.96,318.99,319.02,319.05,319.08,319.11,319.14,
     A 319.17,319.20,319.23,319.26,319.29,319.32,319.35,319.38,319.41,
     A 319.44,319.47,319.50,319.53,319.56,319.59,319.62,319.65,319.68,
     A 319.71,319.74,319.77,319.80,319.83,319.86,319.89,319.92,319.95,
     A 319.98,320.01,320.04,320.07,320.10,320.13,320.16,320.19,320.22,
     A 320.25,320.28,320.31,320.34,320.37,320.40/
       data y1(1:nso2)/.9744E-18,.7714E-18,.1137E-17,.6496E-18,
     A .1015E-17,.7308E-18,.1137E-17,.6902E-18,.1096E-17,.6496E-18,
     A .1137E-17,.9338E-18,.8932E-18,.1096E-17,.6090E-18,.1056E-17,
     A .5684E-18,.1056E-17,.4872E-18,.1096E-17,.3938E-18,.1116E-17,
     A .3248E-18,.1458E-17,.1449E-17,.1352E-17,.1293E-17,.1274E-17,
     A .1095E-17,.1061E-17,.1042E-17,.9794E-18,.9533E-18,.9039E-18,
     A .8546E-18,.8444E-18,.7709E-18,.7653E-18,.7301E-18,.6703E-18,
     A .6352E-18,.6045E-18,.5789E-18,.5550E-18,.5114E-18,.4705E-18,
     A .4562E-18,.4603E-18,.4473E-18,.3968E-18,.3548E-18,.3526E-18,
     A .3446E-18,.3489E-18,.3374E-18,.3338E-18,.2977E-18,.2881E-18,
     A .2868E-18,.3043E-18,.2594E-18,.2520E-18,.2287E-18,.2235E-18,
     A .2106E-18,.2362E-18,.2385E-18,.2075E-18,.1778E-18,.1708E-18,
     A .1978E-18,.1899E-18,.1931E-18,.1842E-18,.1898E-18,.1690E-18,
     A .1885E-18,.2037E-18,.2711E-18,.2854E-18,.3126E-18,.2943E-18,
     A .3641E-18,.4414E-18,.5524E-18,.6181E-18,.7100E-18,.7583E-18,
     A .8579E-18,.9400E-18,.1033E-17,.1014E-17,.1038E-17,.1090E-17,
     A .1123E-17,.1078E-17,.1032E-17,.1019E-17,.9648E-18,.9248E-18,
     A .9110E-18,.9316E-18,.8568E-18,.8183E-18,.7679E-18,.7300E-18,
     A .6639E-18,.6722E-18,.5949E-18,.5323E-18,.5066E-18,.4520E-18,
     A .4728E-18,.4386E-18,.4305E-18,.3980E-18,.3732E-18,.3609E-18,
     A .3485E-18,.3288E-18,.3047E-18,.2854E-18,.2789E-18,.2632E-18,
     A .2373E-18,.2374E-18,.2334E-18,.2199E-18,.2104E-18,.1813E-18,
     A .1793E-18,.1798E-18,.1681E-18,.1778E-18,.1578E-18,.1606E-18,
     A .1594E-18,.1701E-18,.1668E-18,.1909E-18,.1811E-18,.1787E-18,
     A .1789E-18,.1787E-18,.2002E-18,.2004E-18,.1984E-18,.2056E-18,
     A .2196E-18,.2215E-18,.2461E-18,.2352E-18,.2640E-18,.3089E-18,
     A .3675E-18,.4289E-18,.5174E-18,.5418E-18,.6201E-18,.7430E-18,
     A .9542E-18,.1120E-17,.1096E-17,.1007E-17,.9071E-18,.8700E-18,
     A .8569E-18,.8047E-18,.7922E-18,.7688E-18,.6876E-18,.6362E-18,
     A .6718E-18,.6355E-18,.5775E-18,.5716E-18,.5420E-18,.4814E-18,
     A .4234E-18,.4290E-18,.3784E-18,.3329E-18,.3244E-18,.3423E-18,
     A .3283E-18,.2904E-18,.2772E-18,.2757E-18,.2826E-18,.2675E-18,
     A .2384E-18,.2102E-18,.2058E-18,.2083E-18,.2029E-18,.1932E-18,
     A .1713E-18,.1633E-18,.1585E-18,.1525E-18,.1555E-18,.1492E-18,
     A .1330E-18,.1192E-18,.1299E-18,.1292E-18,.1311E-18,.1481E-18,
     A .1453E-18,.1398E-18,.1470E-18,.1567E-18,.1515E-18,.1545E-18,
     A .1617E-18,.1576E-18,.1713E-18,.1908E-18,.2079E-18,.2318E-18,
     A .2526E-18,.2783E-18,.3572E-18,.3672E-18,.3880E-18,.3938E-18,
     A .4158E-18,.4691E-18,.5290E-18,.5420E-18,.5497E-18,.5327E-18,
     A .6065E-18,.6702E-18,.6761E-18,.6910E-18,.7491E-18,.7652E-18,
     A .7794E-18,.7809E-18,.7682E-18,.7191E-18,.6557E-18,.6069E-18,
     A .5960E-18,.5482E-18,.5198E-18,.4610E-18,.4286E-18,.4287E-18,
     A .4067E-18,.3578E-18,.3414E-18,.3291E-18,.3126E-18,.3028E-18,
     A .3022E-18,.2865E-18,.2688E-18,.2658E-18,.2743E-18,.2555E-18,
     A .2591E-18,.2346E-18,.2308E-18,.2028E-18,.1932E-18,.1810E-18,
     A .1770E-18,.1842E-18,.1664E-18,.1720E-18,.1652E-18,.1655E-18,
     A .1686E-18,.1696E-18,.1941E-18,.1729E-18,.1885E-18,.1912E-18,
     A .1754E-18,.1822E-18,.1793E-18,.2015E-18,.1846E-18,.1818E-18,
     A .1647E-18,.1654E-18,.1746E-18,.1852E-18,.1913E-18,.1911E-18,
     A .1985E-18,.2187E-18,.2432E-18,.2793E-18,.2843E-18,.2669E-18,
     A .2661E-18,.2456E-18,.2744E-18,.2750E-18,.3477E-18,.4052E-18,
     A .4230E-18,.4685E-18,.5388E-18,.5477E-18,.6264E-18,.6552E-18,
     A .6611E-18,.6292E-18,.5853E-18,.5260E-18,.5168E-18,.4613E-18,
     A .4562E-18,.4397E-18,.4352E-18,.3723E-18,.3511E-18,.3293E-18,
     A .3018E-18,.2908E-18,.2772E-18,.2518E-18,.2369E-18,.2291E-18,
     A .2308E-18,.2107E-18,.2041E-18,.1905E-18,.1878E-18,.1852E-18,
     A .1779E-18,.1685E-18,.1612E-18,.1359E-18,.1326E-18,.1404E-18,
     A .1379E-18,.1487E-18,.1486E-18,.1423E-18,.1537E-18,.1560E-18,
     A .1581E-18,.1594E-18,.1620E-18,.1502E-18,.1517E-18,.1504E-18,
     A .1549E-18,.1375E-18,.1444E-18,.1538E-18,.1492E-18,.1476E-18,
     A .1758E-18,.1754E-18,.1876E-18,.1745E-18,.1755E-18,.1796E-18,
     A .1792E-18,.1849E-18,.1857E-18,.1846E-18,.1840E-18,.1906E-18,
     A .1964E-18,.2434E-18,.2803E-18,.3220E-18,.3313E-18,.3651E-18,
     A .3880E-18,.4281E-18,.4248E-18,.4158E-18,.3897E-18,.3472E-18,
     A .3378E-18,.3261E-18,.2987E-18,.2757E-18,.2644E-18,.2553E-18,
     A .2713E-18,.2820E-18,.2575E-18,.2511E-18,.2584E-18,.2367E-18,
     A .2161E-18,.2150E-18,.2150E-18,.2155E-18,.2063E-18,.1884E-18,
     A .1746E-18,.1626E-18,.1611E-18,.1561E-18,.1317E-18,.1244E-18,
     A .1222E-18,.1220E-18,.1102E-18,.9686E-19,.9890E-19,.1048E-18,
     A .1025E-18,.9563E-19,.9875E-19,.9610E-19,.9905E-19,.1066E-18,
     A .1106E-18,.1206E-18,.1275E-18,.1279E-18,.1393E-18,.1331E-18,
     A .1337E-18,.1420E-18,.1398E-18,.1250E-18,.1214E-18,.1216E-18,
     A .1211E-18,.1152E-18,.1193E-18,.1238E-18,.1220E-18,.1233E-18,
     A .1211E-18,.1317E-18,.1376E-18,.1431E-18,.1435E-18,.1470E-18,
     A .1526E-18,.1498E-18,.1432E-18,.1390E-18,.1332E-18,.1363E-18,
     A .1508E-18,.1647E-18,.1868E-18,.1982E-18,.2044E-18,.2296E-18,
     A .2532E-18,.2811E-18,.2942E-18,.3056E-18,.3000E-18,.2854E-18,
     A .2703E-18,.2559E-18,.2502E-18,.2338E-18,.2266E-18,.2187E-18,
     A .2086E-18,.1958E-18,.1807E-18,.1652E-18,.1548E-18,.1426E-18,
     A .1339E-18,.1270E-18,.1201E-18,.1165E-18,.1103E-18,.1064E-18,
     A .1056E-18,.1011E-18,.9842E-19,.1222E-18,.1386E-18,.1241E-18,
     A .1191E-18,.1107E-18,.1057E-18,.9674E-19,.9306E-19,.8353E-19,
     A .8314E-19,.7906E-19,.7599E-19,.8161E-19,.7207E-19,.6425E-19,
     A .6574E-19,.7085E-19,.6545E-19,.6826E-19,.7997E-19,.8297E-19,
     A .7719E-19,.7352E-19,.8277E-19,.1023E-18,.1070E-18,.1057E-18,
     A .1093E-18,.1080E-18,.1054E-18,.1058E-18,.1044E-18,.1140E-18,
     A .1202E-18,.1162E-18,.1194E-18,.1191E-18,.1018E-18,.9313E-19,
     A .9804E-19,.9138E-19,.9043E-19,.9508E-19,.9708E-19,.9832E-19,
     A .1045E-18,.1158E-18,.1410E-18,.1558E-18,.1530E-18,.1536E-18,
     A .1518E-18,.1426E-18,.1413E-18,.1316E-18,.1283E-18,.1263E-18,
     A .1114E-18,.1060E-18,.1135E-18,.1321E-18,.1208E-18,.1195E-18,
     A .1211E-18,.1457E-18,.1516E-18,.1533E-18,.1518E-18,.1399E-18,
     A .1304E-18,.1126E-18,.1028E-18,.1083E-18,.9960E-19,.9418E-19,
     A .8286E-19,.8049E-19,.8219E-19,.7231E-19,.6265E-19,.6607E-19,
     A .6718E-19,.6381E-19,.6012E-19,.6105E-19,.5813E-19,.5854E-19,
     A .6199E-19,.5408E-19,.5615E-19,.6364E-19,.5939E-19,.6273E-19,
     A .6861E-19,.7715E-19,.8488E-19,.9438E-19,.9284E-19,.8905E-19,
     A .8386E-19,.8536E-19,.9176E-19,.9322E-19,.8907E-19,.8161E-19,
     A .8337E-19,.9234E-19,.8779E-19,.8792E-19,.8532E-19,.8274E-19,
     A .8855E-19,.9435E-19,.9637E-19,.1000E-18,.1087E-18,.1066E-18,
     A .1040E-18,.9610E-19,.9075E-19,.8583E-19,.7762E-19,.7748E-19,
     A .7170E-19,.7259E-19,.6860E-19,.6743E-19,.6486E-19,.6145E-19,
     A .6533E-19,.7024E-19,.7039E-19,.7014E-19,.6979E-19,.6606E-19,
     A .6728E-19,.6761E-19,.7415E-19,.7468E-19,.7456E-19,.7851E-19,
     A .8124E-19,.9392E-19,.1115E-18,.1171E-18,.1118E-18,.1048E-18,
     A .9506E-19,.9720E-19,.8582E-19,.7631E-19,.8253E-19,.7262E-19,
     A .7373E-19,.7645E-19,.7342E-19,.6153E-19,.6094E-19,.6180E-19,
     A .5933E-19,.5025E-19,.4567E-19,.4381E-19,.4209E-19,.3942E-19,
     A .3871E-19,.3328E-19,.3158E-19,.2936E-19,.3242E-19,.2941E-19,
     A .2903E-19,.2743E-19,.2840E-19,.3129E-19,.3039E-19,.3098E-19,
     A .2791E-19,.3234E-19,.3246E-19,.3368E-19,.3365E-19,.3243E-19,
     A .3230E-19,.3195E-19,.3437E-19,.3723E-19,.4169E-19,.4617E-19,
     A .4482E-19,.5742E-19,.5677E-19,.5761E-19,.6062E-19,.6868E-19,
     A .6984E-19,.6796E-19,.6442E-19,.5720E-19,.5497E-19,.5883E-19,
     A .6610E-19,.6322E-19,.7501E-19,.6537E-19,.5746E-19,.5044E-19,
     A .4889E-19,.4719E-19,.4620E-19,.4637E-19,.4503E-19,.3860E-19,
     A .4036E-19,.3601E-19,.3343E-19,.3329E-19,.2985E-19,.2973E-19,
     A .2787E-19,.2718E-19,.2767E-19,.2498E-19/

      REAL yg(kw)
      REAL a1, a2, dum
      INTEGER ierr
      INTEGER i, l, n, idum
      CHARACTER*40 fil
*_______________________________________________________________________

************* absorption cross sections:
* SO2 absorption cross sections from J. Quant. Spectrosc. Radiat. Transfer
* 37, 165-182, 1987, T. J. McGee and J. Burris Jr.
* Angstrom vs. cm2/molecule, value at 221 K

      fil = 'DATA/McGee87'

      n = nso2

      CALL addpnt(x1,y1,kdata,n,x1(1)*(1.-deltax),0.)
      CALL addpnt(x1,y1,kdata,n,          0.,0.)
      CALL addpnt(x1,y1,kdata,n,x1(n)*(1.+deltax),0.)
      CALL addpnt(x1,y1,kdata,n,      1.e+38,0.)
      CALL inter2(nw,wl,yg,n,x1,y1,ierr)
      IF (ierr .NE. 0) THEN
         WRITE(*,*) ierr, fil
         STOP
      ENDIF
      
      DO 13, l = 1, nw-1
         xsso2(l) = yg(l)
   13 CONTINUE

*_______________________________________________________________________

      RETURN
      END

      SUBROUTINE read1(nw,wl,f)

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Read extra-terrestrial flux data.  Re-grid data to match specified       =*
*=  working wavelength grid.                                                 =*
*-----------------------------------------------------------------------------*
*=  PARAMETERS:                                                              =*
*=  NW     - INTEGER, number of specified intervals + 1 in working        (I)=*
*=           wavelength grid                                                 =*
*=  WL     - REAL, vector of lower limits of wavelength intervals in      (I)=*
*=           working wavelength grid                                         =*
*=  F      - REAL, spectral irradiance at the top of the atmosphere at    (O)=*
*=           each specified wavelength                                       =*
*-----------------------------------------------------------------------------*
*=  EDIT HISTORY:                                                            =*
*=  02/97  Changed offset for grid-end interpolation to relative number      =*
*=         (x * (1 +- deltax))                                               =*
*-----------------------------------------------------------------------------*
*= This program is free software;  you can redistribute it and/or modify     =*
*= it under the terms of the GNU General Public License as published by the  =*
*= Free Software Foundation;  either version 2 of the license, or (at your   =*
*= option) any later version.                                                =*
*= The TUV package is distributed in the hope that it will be useful, but    =*
*= WITHOUT ANY WARRANTY;  without even the implied warranty of MERCHANTIBI-  =*
*= LITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public     =*
*= License for more details.                                                 =*
*= To obtain a copy of the GNU General Public License, write to:             =*
*= Free Software Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.   =*
*-----------------------------------------------------------------------------*
*= To contact the authors, please mail to:                                   =*
*= Sasha Madronich, NCAR/ACD, P.O.Box 3000, Boulder, CO, 80307-3000, USA  or =*
*= send email to:  sasha@ucar.edu                                            =*
*-----------------------------------------------------------------------------*
*= Copyright (C) 1994,95,96  University Corporation for Atmospheric Research =*
*-----------------------------------------------------------------------------*

      IMPLICIT NONE
      INCLUDE 'tuv.params'

* input: (wavelength grid)
      INTEGER nw
      REAL wl(kw)
      INTEGER iw

* output: (extra terrestrial solar flux)
      REAL f(kw)

* local:
      integer nsusim
      parameter(nsusim=5590)
      REAL lambda_hi(10000),irrad_hi(10000)
       data irrad_hi(1:nsusim)/.6525E-01,.2662E+00,.5396E+00,.5506E+00,
     A .3325E+00,.5456E-01,.5860E-01,.6201E-01,.5934E-01,.6588E-01,
     A .6787E-01,.8138E-01,.8514E-01,.1009E+00,.1272E+00,.1409E+00,
     A .1875E+00,.3042E+00,.5199E+00,.1469E+01,.6272E+01,.1400E+02,
     A .2645E+02,.2873E+02,.1863E+02,.6863E+01,.1381E+01,.4207E+00,
     A .3229E+00,.2668E+00,.1765E+00,.1098E+00,.8889E-01,.8038E-01,
     A .6954E-01,.5756E-01,.4989E-01,.4755E-01,.4633E-01,.3926E-01,
     A .3735E-01,.3612E-01,.3001E-01,.3038E-01,.2856E-01,.2802E-01,
     A .2865E-01,.2814E-01,.2651E-01,.2626E-01,.2577E-01,.2378E-01,
     A .2152E-01,.2070E-01,.2084E-01,.2112E-01,.2122E-01,.1988E-01,
     A .2006E-01,.1868E-01,.1921E-01,.1972E-01,.2055E-01,.2020E-01,
     A .1961E-01,.1664E-01,.3151E-01,.5076E-01,.6939E-01,.4579E-01,
     A .2100E-01,.1886E-01,.1938E-01,.2038E-01,.3138E-01,.4421E-01,
     A .4507E-01,.2842E-01,.1492E-01,.1453E-01,.1477E-01,.1654E-01,
     A .1496E-01,.1512E-01,.1931E-01,.1759E-01,.1478E-01,.1514E-01,
     A .1533E-01,.1526E-01,.1865E-01,.2031E-01,.1695E-01,.1254E-01,
     A .1234E-01,.1606E-01,.2126E-01,.2602E-01,.2098E-01,.1556E-01,
     A .1403E-01,.1755E-01,.2363E-01,.2201E-01,.1753E-01,.1818E-01,
     A .2366E-01,.3283E-01,.3711E-01,.4311E-01,.5027E-01,.4084E-01,
     A .2634E-01,.1812E-01,.1471E-01,.9788E-02,.1084E-01,.3914E-01,
     A .8849E-01,.1015E+00,.4607E-01,.2981E-01,.2373E-01,.2186E-01,
     A .1755E-01,.1231E-01,.1225E-01,.1165E-01,.1072E-01,.1255E-01,
     A .1221E-01,.1257E-01,.1345E-01,.1050E-01,.1163E-01,.1095E-01,
     A .1139E-01,.1562E-01,.1928E-01,.2403E-01,.2547E-01,.2053E-01,
     A .2124E-01,.2473E-01,.2479E-01,.1914E-01,.1561E-01,.1617E-01,
     A .1914E-01,.2501E-01,.2105E-01,.1789E-01,.1548E-01,.1397E-01,
     A .1888E-01,.1915E-01,.1883E-01,.2021E-01,.1900E-01,.1753E-01,
     A .1462E-01,.1410E-01,.1486E-01,.1229E-01,.1376E-01,.2190E-01,
     A .3000E-01,.2177E-01,.2008E-01,.1981E-01,.1648E-01,.1494E-01,
     A .1392E-01,.1306E-01,.1302E-01,.1233E-01,.1151E-01,.1445E-01,
     A .1617E-01,.1631E-01,.2081E-01,.2509E-01,.2367E-01,.1748E-01,
     A .1444E-01,.1478E-01,.1952E-01,.2902E-01,.2574E-01,.2263E-01,
     A .2445E-01,.3403E-01,.9062E-01,.1774E+00,.2312E+00,.1893E+00,
     A .8462E-01,.3443E-01,.1105E+00,.2424E+00,.2968E+00,.3120E+00,
     A .2960E+00,.1898E+00,.6908E-01,.1702E-01,.1476E-01,.3080E-01,
     A .5714E-01,.7034E-01,.4775E-01,.2319E-01,.2641E-01,.3111E-01,
     A .2686E-01,.2085E-01,.2240E-01,.2512E-01,.1888E-01,.1369E-01,
     A .1459E-01,.1775E-01,.2639E-01,.3201E-01,.3118E-01,.2551E-01,
     A .2162E-01,.2280E-01,.2998E-01,.2733E-01,.1349E-01,.1325E-01,
     A .1163E-01,.1107E-01,.1132E-01,.1346E-01,.1299E-01,.1736E-01,
     A .2097E-01,.1822E-01,.1716E-01,.1118E-01,.1105E-01,.1455E-01,
     A .1795E-01,.1827E-01,.1440E-01,.1587E-01,.2382E-01,.3352E-01,
     A .2879E-01,.1945E-01,.1503E-01,.1534E-01,.1475E-01,.1599E-01,
     A .1237E-01,.1136E-01,.6356E-01,.2897E+00,.5284E+00,.5910E+00,
     A .6725E+00,.4339E+00,.9271E-01,.3915E-01,.1607E-01,.1720E-01,
     A .1534E-01,.1237E-01,.1478E-01,.1589E-01,.1583E-01,.1665E-01,
     A .1980E-01,.2100E-01,.1761E-01,.1828E-01,.2007E-01,.1679E-01,
     A .1470E-01,.1749E-01,.1789E-01,.1846E-01,.2181E-01,.1990E-01,
     A .1958E-01,.1739E-01,.1937E-01,.2132E-01,.1939E-01,.1855E-01,
     A .3103E-01,.5343E-01,.5397E-01,.3402E-01,.2322E-01,.2148E-01,
     A .2398E-01,.3002E-01,.5305E-01,.8904E-01,.8065E-01,.4746E-01,
     A .3437E-01,.3091E-01,.4061E-01,.5062E-01,.4852E-01,.4130E-01,
     A .2762E-01,.2379E-01,.2958E-01,.3144E-01,.2788E-01,.2173E-01,
     A .2212E-01,.2968E-01,.3500E-01,.3352E-01,.2461E-01,.2137E-01,
     A .2201E-01,.2200E-01,.2469E-01,.2523E-01,.2710E-01,.2611E-01,
     A .2707E-01,.2811E-01,.2769E-01,.3155E-01,.3084E-01,.3007E-01,
     A .3166E-01,.2961E-01,.2698E-01,.2663E-01,.2842E-01,.2880E-01,
     A .2577E-01,.2518E-01,.2551E-01,.2764E-01,.2751E-01,.2540E-01,
     A .2707E-01,.2631E-01,.2791E-01,.3022E-01,.2816E-01,.2760E-01,
     A .3082E-01,.3316E-01,.2903E-01,.2685E-01,.2794E-01,.2532E-01,
     A .2608E-01,.3049E-01,.2658E-01,.2732E-01,.2561E-01,.2894E-01,
     A .3353E-01,.2849E-01,.2801E-01,.2895E-01,.2640E-01,.2777E-01,
     A .2919E-01,.2700E-01,.2631E-01,.2870E-01,.3879E-01,.7091E-01,
     A .1449E+00,.2934E+00,.2967E+00,.9113E-01,.4227E-01,.3434E-01,
     A .3272E-01,.3456E-01,.3210E-01,.3025E-01,.3266E-01,.3620E-01,
     A .3877E-01,.3953E-01,.4595E-01,.6574E-01,.6243E-01,.7330E-01,
     A .1252E+00,.1783E+00,.1213E+00,.5062E-01,.4826E-01,.4723E-01,
     A .4512E-01,.4325E-01,.4120E-01,.3491E-01,.3862E-01,.4219E-01,
     A .3752E-01,.3801E-01,.3576E-01,.3334E-01,.3763E-01,.3993E-01,
     A .4072E-01,.4044E-01,.4178E-01,.4262E-01,.4033E-01,.3888E-01,
     A .3604E-01,.3475E-01,.3812E-01,.4050E-01,.4254E-01,.4256E-01,
     A .3853E-01,.3535E-01,.3447E-01,.3752E-01,.4052E-01,.3776E-01,
     A .3507E-01,.3588E-01,.3890E-01,.3917E-01,.3823E-01,.3813E-01,
     A .4022E-01,.4424E-01,.6461E-01,.7802E-01,.6645E-01,.5313E-01,
     A .4047E-01,.4078E-01,.4182E-01,.4167E-01,.4066E-01,.3998E-01,
     A .4041E-01,.4000E-01,.4118E-01,.4309E-01,.4869E-01,.6081E-01,
     A .6687E-01,.6283E-01,.6106E-01,.5276E-01,.4347E-01,.4997E-01,
     A .5148E-01,.4691E-01,.5192E-01,.5378E-01,.5141E-01,.4960E-01,
     A .4481E-01,.4254E-01,.4232E-01,.4276E-01,.4821E-01,.5080E-01,
     A .4835E-01,.4845E-01,.5251E-01,.5307E-01,.4949E-01,.5244E-01,
     A .5429E-01,.5302E-01,.5303E-01,.5110E-01,.5599E-01,.5393E-01,
     A .5057E-01,.5013E-01,.5094E-01,.4768E-01,.4971E-01,.5072E-01,
     A .4695E-01,.4598E-01,.4685E-01,.4622E-01,.4953E-01,.5204E-01,
     A .5209E-01,.5407E-01,.5689E-01,.5587E-01,.5751E-01,.5748E-01,
     A .5295E-01,.4759E-01,.5185E-01,.5410E-01,.5584E-01,.5889E-01,
     A .5693E-01,.5355E-01,.5541E-01,.5950E-01,.5565E-01,.5174E-01,
     A .5672E-01,.6284E-01,.6859E-01,.6612E-01,.5806E-01,.5729E-01,
     A .5822E-01,.5940E-01,.6289E-01,.7520E-01,.8892E-01,.9071E-01,
     A .8056E-01,.6741E-01,.6243E-01,.6809E-01,.7086E-01,.6596E-01,
     A .6634E-01,.7142E-01,.9957E-01,.1382E+00,.1422E+00,.1414E+00,
     A .1141E+00,.9059E-01,.6669E-01,.6287E-01,.6457E-01,.6491E-01,
     A .6415E-01,.6301E-01,.6312E-01,.5980E-01,.5974E-01,.5944E-01,
     A .6259E-01,.8395E-01,.1109E+00,.1235E+00,.1030E+00,.9986E-01,
     A .9631E-01,.7233E-01,.7046E-01,.7452E-01,.7807E-01,.9309E-01,
     A .9472E-01,.8904E-01,.8616E-01,.7531E-01,.6988E-01,.7000E-01,
     A .7012E-01,.6890E-01,.7054E-01,.7273E-01,.7591E-01,.7103E-01,
     A .8591E-01,.8651E-01,.7515E-01,.7627E-01,.7531E-01,.7568E-01,
     A .7408E-01,.7428E-01,.7596E-01,.7428E-01,.7395E-01,.7905E-01,
     A .8245E-01,.7662E-01,.7932E-01,.8283E-01,.8404E-01,.8406E-01,
     A .8030E-01,.8875E-01,.8763E-01,.8456E-01,.8382E-01,.8541E-01,
     A .8794E-01,.8510E-01,.8572E-01,.8681E-01,.8622E-01,.8505E-01,
     A .8472E-01,.8301E-01,.8501E-01,.8714E-01,.8786E-01,.9348E-01,
     A .8528E-01,.1100E+00,.1098E+00,.1054E+00,.9757E-01,.8746E-01,
     A .9173E-01,.9614E-01,.8344E-01,.8795E-01,.9245E-01,.9174E-01,
     A .9108E-01,.8811E-01,.8758E-01,.9153E-01,.9323E-01,.9232E-01,
     A .9658E-01,.9831E-01,.8579E-01,.9206E-01,.9464E-01,.9631E-01,
     A .9530E-01,.9868E-01,.1033E+00,.1054E+00,.1016E+00,.9994E-01,
     A .9918E-01,.1032E+00,.1458E+00,.2007E+00,.2054E+00,.1649E+00,
     A .1056E+00,.1110E+00,.1103E+00,.1049E+00,.9826E-01,.1031E+00,
     A .1088E+00,.1126E+00,.1091E+00,.1443E+00,.1950E+00,.2183E+00,
     A .1894E+00,.1350E+00,.1146E+00,.1062E+00,.1056E+00,.1160E+00,
     A .1461E+00,.1415E+00,.1238E+00,.1041E+00,.9938E-01,.1006E+00,
     A .1111E+00,.1300E+00,.1462E+00,.1408E+00,.1348E+00,.1115E+00,
     A .9608E-01,.9929E-01,.1040E+00,.1153E+00,.1247E+00,.1193E+00,
     A .1184E+00,.9577E-01,.1766E+00,.4678E+00,.6565E+00,.6417E+00,
     A .3351E+00,.1256E+00,.2067E+00,.3868E+00,.4074E+00,.2281E+00,
     A .1591E+00,.1361E+00,.1059E+00,.1035E+00,.9743E-01,.1057E+00,
     A .1231E+00,.1194E+00,.1130E+00,.1200E+00,.1161E+00,.1232E+00,
     A .1641E+00,.2015E+00,.2270E+00,.2461E+00,.2910E+00,.3284E+00,
     A .3568E+00,.2963E+00,.1710E+00,.1209E+00,.1400E+00,.1671E+00,
     A .1582E+00,.1281E+00,.1218E+00,.1165E+00,.1299E+00,.1456E+00,
     A .1573E+00,.1550E+00,.1513E+00,.1374E+00,.1439E+00,.1838E+00,
     A .1975E+00,.1715E+00,.1609E+00,.1494E+00,.1366E+00,.1444E+00,
     A .1583E+00,.1868E+00,.2309E+00,.2407E+00,.2184E+00,.1762E+00,
     A .1437E+00,.1495E+00,.1518E+00,.1431E+00,.1335E+00,.1310E+00,
     A .1202E+00,.1264E+00,.1635E+00,.1872E+00,.1713E+00,.1404E+00,
     A .1241E+00,.1351E+00,.1460E+00,.1499E+00,.1677E+00,.1835E+00,
     A .1769E+00,.1623E+00,.1551E+00,.1502E+00,.1546E+00,.1842E+00,
     A .2099E+00,.1886E+00,.1619E+00,.1642E+00,.1681E+00,.1650E+00,
     A .1615E+00,.1665E+00,.1702E+00,.1675E+00,.1657E+00,.1714E+00,
     A .1843E+00,.1956E+00,.1908E+00,.1738E+00,.1529E+00,.1384E+00,
     A .1452E+00,.1572E+00,.1576E+00,.1496E+00,.1447E+00,.1596E+00,
     A .1770E+00,.1743E+00,.1714E+00,.1911E+00,.2208E+00,.2199E+00,
     A .1971E+00,.1709E+00,.1628E+00,.1664E+00,.1734E+00,.1819E+00,
     A .1792E+00,.1743E+00,.1762E+00,.2029E+00,.2240E+00,.2100E+00,
     A .1932E+00,.1867E+00,.2125E+00,.2562E+00,.2474E+00,.2218E+00,
     A .2417E+00,.2763E+00,.2908E+00,.2700E+00,.2375E+00,.2121E+00,
     A .2032E+00,.1950E+00,.1949E+00,.2041E+00,.2101E+00,.2142E+00,
     A .2293E+00,.2267E+00,.2048E+00,.2119E+00,.2137E+00,.2166E+00,
     A .2516E+00,.2698E+00,.2634E+00,.2699E+00,.2670E+00,.2462E+00,
     A .2296E+00,.2404E+00,.2999E+00,.3263E+00,.2816E+00,.2469E+00,
     A .2544E+00,.2479E+00,.2420E+00,.2707E+00,.2861E+00,.2738E+00,
     A .2541E+00,.2283E+00,.2272E+00,.2611E+00,.2812E+00,.2961E+00,
     A .3076E+00,.3116E+00,.3064E+00,.2644E+00,.2318E+00,.2375E+00,
     A .2627E+00,.2954E+00,.3174E+00,.3064E+00,.2614E+00,.2553E+00,
     A .3234E+00,.4168E+00,.4697E+00,.4097E+00,.3063E+00,.2711E+00,
     A .2740E+00,.3024E+00,.3398E+00,.3477E+00,.3245E+00,.2990E+00,
     A .2941E+00,.2788E+00,.2788E+00,.2788E+00,.2742E+00,.2845E+00,
     A .2840E+00,.2896E+00,.3080E+00,.3035E+00,.2790E+00,.2710E+00,
     A .2676E+00,.2694E+00,.2785E+00,.2727E+00,.2798E+00,.3125E+00,
     A .3523E+00,.3390E+00,.3658E+00,.4814E+00,.7587E+00,.9372E+00,
     A .1042E+01,.1032E+01,.8773E+00,.6644E+00,.5012E+00,.4159E+00,
     A .3792E+00,.3311E+00,.3402E+00,.3410E+00,.3413E+00,.3826E+00,
     A .3832E+00,.3730E+00,.3380E+00,.3232E+00,.3406E+00,.3709E+00,
     A .3971E+00,.3717E+00,.3475E+00,.3385E+00,.3335E+00,.3288E+00,
     A .3570E+00,.3813E+00,.4438E+00,.5991E+00,.5330E+00,.3716E+00,
     A .3159E+00,.3220E+00,.3457E+00,.4062E+00,.4433E+00,.4368E+00,
     A .3886E+00,.3545E+00,.3842E+00,.4072E+00,.4134E+00,.4070E+00,
     A .3807E+00,.3800E+00,.3977E+00,.4109E+00,.4049E+00,.4146E+00,
     A .4427E+00,.4398E+00,.4093E+00,.4033E+00,.4183E+00,.4193E+00,
     A .4187E+00,.4211E+00,.4326E+00,.4695E+00,.5398E+00,.5902E+00,
     A .5437E+00,.4858E+00,.4952E+00,.5187E+00,.5073E+00,.5029E+00,
     A .5147E+00,.5271E+00,.5897E+00,.6094E+00,.5730E+00,.5633E+00,
     A .5669E+00,.5723E+00,.5957E+00,.6360E+00,.6464E+00,.6178E+00,
     A .5964E+00,.6616E+00,.7383E+00,.7260E+00,.6429E+00,.6291E+00,
     A .6683E+00,.6818E+00,.6574E+00,.6479E+00,.6859E+00,.7953E+00,
     A .8769E+00,.8611E+00,.7538E+00,.6823E+00,.7078E+00,.7135E+00,
     A .6779E+00,.6572E+00,.6760E+00,.6815E+00,.6712E+00,.6663E+00,
     A .6960E+00,.7535E+00,.7777E+00,.7774E+00,.7856E+00,.7606E+00,
     A .7388E+00,.7333E+00,.6500E+00,.6655E+00,.7213E+00,.7198E+00,
     A .6573E+00,.6142E+00,.6315E+00,.6659E+00,.7177E+00,.7702E+00,
     A .8016E+00,.8240E+00,.8204E+00,.7694E+00,.7518E+00,.7381E+00,
     A .7534E+00,.8217E+00,.8665E+00,.8323E+00,.7579E+00,.7254E+00,
     A .7267E+00,.7209E+00,.7035E+00,.7810E+00,.8567E+00,.8750E+00,
     A .8148E+00,.8074E+00,.8255E+00,.8252E+00,.8229E+00,.8083E+00,
     A .7613E+00,.6811E+00,.6619E+00,.6772E+00,.6900E+00,.7333E+00,
     A .7512E+00,.7401E+00,.7637E+00,.7732E+00,.7770E+00,.8001E+00,
     A .7910E+00,.7880E+00,.8165E+00,.8233E+00,.8115E+00,.8138E+00,
     A .8226E+00,.8344E+00,.8623E+00,.8655E+00,.8682E+00,.8687E+00,
     A .8779E+00,.9340E+00,.9895E+00,.9841E+00,.9885E+00,.9906E+00,
     A .9905E+00,.9572E+00,.1031E+01,.1023E+01,.1015E+01,.1043E+01,
     A .9634E+00,.8697E+00,.8767E+00,.9202E+00,.9493E+00,.9523E+00,
     A .9575E+00,.1035E+01,.1076E+01,.1109E+01,.1121E+01,.1102E+01,
     A .1080E+01,.1113E+01,.1153E+01,.1195E+01,.1242E+01,.1229E+01,
     A .1202E+01,.1183E+01,.1176E+01,.1218E+01,.1218E+01,.1213E+01,
     A .1282E+01,.1253E+01,.1229E+01,.1240E+01,.1253E+01,.1318E+01,
     A .1342E+01,.1280E+01,.1141E+01,.1126E+01,.1153E+01,.1211E+01,
     A .1157E+01,.1085E+01,.1043E+01,.1158E+01,.1409E+01,.1544E+01,
     A .1492E+01,.1372E+01,.1302E+01,.1326E+01,.1351E+01,.1449E+01,
     A .1533E+01,.1626E+01,.1613E+01,.1696E+01,.1822E+01,.1788E+01,
     A .1716E+01,.1557E+01,.1341E+01,.1217E+01,.1228E+01,.1306E+01,
     A .1365E+01,.1371E+01,.1399E+01,.1470E+01,.1620E+01,.1677E+01,
     A .1637E+01,.1707E+01,.1731E+01,.1651E+01,.1510E+01,.1479E+01,
     A .1555E+01,.1631E+01,.1646E+01,.1646E+01,.1630E+01,.1645E+01,
     A .1712E+01,.1687E+01,.1640E+01,.1645E+01,.1653E+01,.1669E+01,
     A .1743E+01,.1848E+01,.1810E+01,.1822E+01,.1868E+01,.1857E+01,
     A .1796E+01,.1613E+01,.1532E+01,.1460E+01,.1435E+01,.1473E+01,
     A .1523E+01,.1603E+01,.1653E+01,.1602E+01,.1599E+01,.1651E+01,
     A .1686E+01,.1696E+01,.1708E+01,.1720E+01,.1738E+01,.1771E+01,
     A .1824E+01,.1863E+01,.1873E+01,.1887E+01,.1850E+01,.1880E+01,
     A .1854E+01,.1856E+01,.1920E+01,.1974E+01,.1976E+01,.1922E+01,
     A .1774E+01,.1880E+01,.2577E+01,.2967E+01,.2749E+01,.2253E+01,
     A .1867E+01,.1902E+01,.2020E+01,.1982E+01,.1926E+01,.1898E+01,
     A .1931E+01,.1893E+01,.1818E+01,.1731E+01,.1785E+01,.1877E+01,
     A .1822E+01,.2386E+01,.3597E+01,.4543E+01,.4393E+01,.3024E+01,
     A .2392E+01,.2227E+01,.2245E+01,.2074E+01,.1995E+01,.2130E+01,
     A .2269E+01,.2268E+01,.2301E+01,.2471E+01,.2582E+01,.2676E+01,
     A .2602E+01,.2498E+01,.2372E+01,.2254E+01,.2237E+01,.2248E+01,
     A .2179E+01,.2141E+01,.2204E+01,.2299E+01,.2326E+01,.2250E+01,
     A .2215E+01,.2248E+01,.2336E+01,.2410E+01,.2459E+01,.2400E+01,
     A .2397E+01,.2505E+01,.2629E+01,.2651E+01,.2593E+01,.2446E+01,
     A .2358E+01,.2456E+01,.2521E+01,.2509E+01,.2566E+01,.2741E+01,
     A .2666E+01,.2380E+01,.2262E+01,.1996E+01,.1945E+01,.2135E+01,
     A .2303E+01,.2360E+01,.2308E+01,.2297E+01,.2347E+01,.2198E+01,
     A .2069E+01,.2079E+01,.2092E+01,.1946E+01,.1730E+01,.1751E+01,
     A .1920E+01,.2190E+01,.2474E+01,.2335E+01,.1879E+01,.1880E+01,
     A .1957E+01,.1997E+01,.2026E+01,.2234E+01,.2395E+01,.2761E+01,
     A .2898E+01,.2906E+01,.2818E+01,.2882E+01,.2782E+01,.2561E+01,
     A .2488E+01,.2500E+01,.2508E+01,.2446E+01,.2424E+01,.2518E+01,
     A .2540E+01,.2529E+01,.2592E+01,.2635E+01,.2754E+01,.2807E+01,
     A .2691E+01,.2664E+01,.2556E+01,.2554E+01,.2502E+01,.2659E+01,
     A .2942E+01,.3032E+01,.3240E+01,.3277E+01,.3238E+01,.3257E+01,
     A .3233E+01,.3270E+01,.3233E+01,.3416E+01,.3407E+01,.3286E+01,
     A .3134E+01,.3084E+01,.3197E+01,.3288E+01,.3231E+01,.3128E+01,
     A .3185E+01,.3166E+01,.3151E+01,.3140E+01,.3155E+01,.3114E+01,
     A .3058E+01,.3015E+01,.2974E+01,.2994E+01,.2961E+01,.2948E+01,
     A .2930E+01,.3045E+01,.3425E+01,.3520E+01,.3453E+01,.3544E+01,
     A .3372E+01,.3378E+01,.3440E+01,.3397E+01,.3443E+01,.3498E+01,
     A .3413E+01,.3342E+01,.3431E+01,.3630E+01,.3641E+01,.3585E+01,
     A .3595E+01,.3773E+01,.4015E+01,.4191E+01,.4086E+01,.3403E+01,
     A .3633E+01,.3721E+01,.3948E+01,.3776E+01,.3709E+01,.3717E+01,
     A .3758E+01,.3613E+01,.3348E+01,.3208E+01,.3119E+01,.3209E+01,
     A .3459E+01,.3616E+01,.3562E+01,.3330E+01,.3331E+01,.3572E+01,
     A .3649E+01,.3787E+01,.3993E+01,.3958E+01,.3828E+01,.3850E+01,
     A .3688E+01,.3726E+01,.3841E+01,.3925E+01,.3936E+01,.4136E+01,
     A .4336E+01,.4405E+01,.4231E+01,.3660E+01,.3757E+01,.3749E+01,
     A .3878E+01,.4035E+01,.3971E+01,.3967E+01,.4044E+01,.4086E+01,
     A .4259E+01,.4457E+01,.4687E+01,.4400E+01,.4385E+01,.4301E+01,
     A .4131E+01,.4058E+01,.4049E+01,.4166E+01,.4296E+01,.4453E+01,
     A .4598E+01,.4802E+01,.4796E+01,.4522E+01,.4398E+01,.4557E+01,
     A .4513E+01,.4324E+01,.4304E+01,.4226E+01,.4424E+01,.4449E+01,
     A .4258E+01,.4286E+01,.4351E+01,.4474E+01,.4499E+01,.4326E+01,
     A .4092E+01,.3703E+01,.3319E+01,.2958E+01,.3020E+01,.2891E+01,
     A .3051E+01,.3362E+01,.3288E+01,.3341E+01,.3363E+01,.2995E+01,
     A .2935E+01,.2899E+01,.2730E+01,.3279E+01,.3709E+01,.4033E+01,
     A .4485E+01,.4735E+01,.4813E+01,.4838E+01,.4901E+01,.5228E+01,
     A .5828E+01,.5869E+01,.5936E+01,.6014E+01,.5722E+01,.5535E+01,
     A .5403E+01,.5569E+01,.5347E+01,.5209E+01,.5170E+01,.5432E+01,
     A .6086E+01,.6212E+01,.6138E+01,.6118E+01,.5714E+01,.5395E+01,
     A .5466E+01,.5260E+01,.5080E+01,.4777E+01,.4946E+01,.5203E+01,
     A .5832E+01,.5888E+01,.5656E+01,.5380E+01,.5368E+01,.5557E+01,
     A .5843E+01,.5783E+01,.5397E+01,.5555E+01,.5896E+01,.6068E+01,
     A .6044E+01,.5838E+01,.5640E+01,.5273E+01,.5289E+01,.5392E+01,
     A .5454E+01,.5495E+01,.5544E+01,.5933E+01,.6400E+01,.6526E+01,
     A .6666E+01,.6699E+01,.6771E+01,.6787E+01,.6806E+01,.6848E+01,
     A .6958E+01,.6912E+01,.6813E+01,.6508E+01,.6184E+01,.6435E+01,
     A .6313E+01,.6274E+01,.6379E+01,.6418E+01,.6171E+01,.6124E+01,
     A .6126E+01,.6192E+01,.6278E+01,.6252E+01,.6079E+01,.5925E+01,
     A .6297E+01,.6388E+01,.6135E+01,.6141E+01,.6150E+01,.6095E+01,
     A .6313E+01,.6654E+01,.6940E+01,.6916E+01,.6142E+01,.6216E+01,
     A .6252E+01,.6331E+01,.6558E+01,.6505E+01,.5841E+01,.5678E+01,
     A .5564E+01,.5749E+01,.6015E+01,.5743E+01,.5609E+01,.5670E+01,
     A .5965E+01,.6203E+01,.6154E+01,.6016E+01,.6228E+01,.6926E+01,
     A .7091E+01,.7309E+01,.7658E+01,.7591E+01,.7121E+01,.6972E+01,
     A .7072E+01,.7024E+01,.7120E+01,.6959E+01,.6623E+01,.6618E+01,
     A .7038E+01,.7169E+01,.6869E+01,.6966E+01,.7029E+01,.6899E+01,
     A .6976E+01,.7318E+01,.7459E+01,.7704E+01,.7459E+01,.7521E+01,
     A .7214E+01,.7186E+01,.7488E+01,.7976E+01,.8012E+01,.7856E+01,
     A .7503E+01,.7498E+01,.7422E+01,.7495E+01,.7606E+01,.7716E+01,
     A .8136E+01,.7572E+01,.7019E+01,.7275E+01,.7633E+01,.7990E+01,
     A .8493E+01,.8364E+01,.8065E+01,.7902E+01,.7826E+01,.7777E+01,
     A .7805E+01,.8002E+01,.8249E+01,.8537E+01,.8971E+01,.9138E+01,
     A .8742E+01,.8551E+01,.8535E+01,.8449E+01,.8050E+01,.7704E+01,
     A .7714E+01,.7756E+01,.8013E+01,.8682E+01,.8850E+01,.8128E+01,
     A .8037E+01,.6694E+01,.5917E+01,.6216E+01,.7735E+01,.8646E+01,
     A .8882E+01,.9033E+01,.9075E+01,.8754E+01,.8756E+01,.9332E+01,
     A .9735E+01,.9981E+01,.9977E+01,.9596E+01,.8920E+01,.8206E+01,
     A .7821E+01,.7909E+01,.8462E+01,.8916E+01,.9725E+01,.1012E+02,
     A .1001E+02,.1009E+02,.9855E+01,.9632E+01,.8873E+01,.8399E+01,
     A .8536E+01,.9288E+01,.1010E+02,.1029E+02,.1107E+02,.1122E+02,
     A .1082E+02,.1027E+02,.9827E+01,.9780E+01,.1016E+02,.1092E+02,
     A .1147E+02,.1112E+02,.1059E+02,.1011E+02,.9945E+01,.1043E+02,
     A .1061E+02,.1039E+02,.1062E+02,.1069E+02,.1076E+02,.1092E+02,
     A .1099E+02,.1118E+02,.1153E+02,.1132E+02,.1127E+02,.1083E+02,
     A .1008E+02,.9602E+01,.1010E+02,.1090E+02,.1132E+02,.1091E+02,
     A .1017E+02,.9872E+01,.1095E+02,.1144E+02,.1159E+02,.1172E+02,
     A .1077E+02,.9720E+01,.9707E+01,.9918E+01,.1049E+02,.1115E+02,
     A .1129E+02,.1169E+02,.1150E+02,.1116E+02,.1092E+02,.1130E+02,
     A .1196E+02,.1215E+02,.1177E+02,.1155E+02,.1159E+02,.1147E+02,
     A .1169E+02,.1206E+02,.1247E+02,.1227E+02,.1212E+02,.1248E+02,
     A .1258E+02,.1304E+02,.1315E+02,.1319E+02,.1343E+02,.1373E+02,
     A .1403E+02,.1428E+02,.1439E+02,.1409E+02,.1324E+02,.1321E+02,
     A .1347E+02,.1438E+02,.1407E+02,.1366E+02,.1324E+02,.1311E+02,
     A .1284E+02,.1297E+02,.1365E+02,.1359E+02,.1322E+02,.1289E+02,
     A .1392E+02,.1598E+02,.1687E+02,.1672E+02,.1623E+02,.1696E+02,
     A .1871E+02,.1919E+02,.1738E+02,.1601E+02,.1563E+02,.1637E+02,
     A .1786E+02,.1865E+02,.1821E+02,.2169E+02,.2534E+02,.2251E+02,
     A .1843E+02,.1807E+02,.1996E+02,.2141E+02,.2258E+02,.2372E+02,
     A .2509E+02,.2695E+02,.2675E+02,.2448E+02,.2137E+02,.1902E+02,
     A .2034E+02,.2364E+02,.2890E+02,.3427E+02,.3444E+02,.3026E+02,
     A .2252E+02,.1749E+02,.2141E+02,.2832E+02,.3551E+02,.3866E+02,
     A .3502E+02,.3369E+02,.3365E+02,.3302E+02,.2819E+02,.2140E+02,
     A .1948E+02,.1938E+02,.2020E+02,.2437E+02,.3386E+02,.3965E+02,
     A .4051E+02,.3484E+02,.3067E+02,.2919E+02,.2619E+02,.2367E+02,
     A .2425E+02,.3158E+02,.4100E+02,.4505E+02,.4665E+02,.4032E+02,
     A .3528E+02,.3616E+02,.4101E+02,.4888E+02,.5173E+02,.4501E+02,
     A .3860E+02,.3514E+02,.3100E+02,.2687E+02,.1992E+02,.1644E+02,
     A .1414E+02,.1343E+02,.1962E+02,.2663E+02,.3207E+02,.3640E+02,
     A .3624E+02,.3639E+02,.3693E+02,.3734E+02,.3387E+02,.2820E+02,
     A .2469E+02,.2613E+02,.2963E+02,.3011E+02,.2947E+02,.2866E+02,
     A .3089E+02,.3129E+02,.2784E+02,.2471E+02,.2589E+02,.2957E+02,
     A .3288E+02,.4055E+02,.4277E+02,.3857E+02,.3189E+02,.2519E+02,
     A .2204E+02,.2652E+02,.3357E+02,.4303E+02,.5199E+02,.5068E+02,
     A .4700E+02,.4630E+02,.4846E+02,.4929E+02,.4666E+02,.4359E+02,
     A .3841E+02,.3443E+02,.3309E+02,.3247E+02,.3318E+02,.3593E+02,
     A .3996E+02,.4196E+02,.4632E+02,.4349E+02,.3583E+02,.3064E+02,
     A .3057E+02,.3282E+02,.3540E+02,.4124E+02,.4622E+02,.4793E+02,
     A .4559E+02,.4120E+02,.3860E+02,.3623E+02,.3337E+02,.2864E+02,
     A .2403E+02,.2167E+02,.2071E+02,.2235E+02,.2366E+02,.2866E+02,
     A .3877E+02,.4378E+02,.4233E+02,.4117E+02,.3848E+02,.4004E+02,
     A .3925E+02,.3706E+02,.3672E+02,.3667E+02,.3162E+02,.2548E+02,
     A .1994E+02,.1759E+02,.2095E+02,.2954E+02,.3132E+02,.2687E+02,
     A .2626E+02,.2965E+02,.3378E+02,.3582E+02,.3595E+02,.3789E+02,
     A .3922E+02,.3757E+02,.3473E+02,.2954E+02,.2628E+02,.2477E+02,
     A .2651E+02,.3143E+02,.3424E+02,.2933E+02,.2449E+02,.2453E+02,
     A .2919E+02,.4161E+02,.5075E+02,.5359E+02,.5413E+02,.5500E+02,
     A .5717E+02,.5920E+02,.5793E+02,.5385E+02,.4900E+02,.4475E+02,
     A .4305E+02,.4786E+02,.5352E+02,.4556E+02,.3061E+02,.2712E+02,
     A .3221E+02,.3937E+02,.4239E+02,.4006E+02,.3797E+02,.4199E+02,
     A .4654E+02,.4460E+02,.3991E+02,.3677E+02,.5015E+02,.6416E+02,
     A .6217E+02,.5898E+02,.5229E+02,.4485E+02,.3850E+02,.3293E+02,
     A .3542E+02,.4118E+02,.4784E+02,.5604E+02,.5988E+02,.5937E+02,
     A .5716E+02,.4821E+02,.3456E+02,.3185E+02,.3530E+02,.4430E+02,
     A .6032E+02,.6485E+02,.6933E+02,.6571E+02,.5878E+02,.5785E+02,
     A .5930E+02,.5601E+02,.4768E+02,.4447E+02,.4123E+02,.3685E+02,
     A .3321E+02,.3888E+02,.3729E+02,.3060E+02,.2186E+02,.1832E+02,
     A .2074E+02,.2358E+02,.4319E+02,.5046E+02,.4638E+02,.4425E+02,
     A .4610E+02,.4701E+02,.4012E+02,.3112E+02,.2467E+02,.2242E+02,
     A .2517E+02,.2963E+02,.3784E+02,.4886E+02,.5275E+02,.5207E+02,
     A .5346E+02,.5692E+02,.5754E+02,.5074E+02,.4878E+02,.4724E+02,
     A .4909E+02,.5227E+02,.5066E+02,.4890E+02,.4939E+02,.4794E+02,
     A .4381E+02,.4175E+02,.4303E+02,.4712E+02,.5351E+02,.5411E+02,
     A .5109E+02,.4744E+02,.4240E+02,.4175E+02,.4987E+02,.6633E+02,
     A .7564E+02,.8278E+02,.8027E+02,.7441E+02,.7416E+02,.7661E+02,
     A .7836E+02,.7665E+02,.7122E+02,.6706E+02,.6250E+02,.5236E+02,
     A .5202E+02,.6074E+02,.6796E+02,.6798E+02,.6239E+02,.5498E+02,
     A .5203E+02,.5148E+02,.5096E+02,.4906E+02,.4649E+02,.4467E+02,
     A .5205E+02,.5901E+02,.6298E+02,.6406E+02,.6760E+02,.7116E+02,
     A .7471E+02,.6832E+02,.5798E+02,.5256E+02,.5250E+02,.5314E+02,
     A .5099E+02,.4848E+02,.5203E+02,.5852E+02,.5854E+02,.5176E+02,
     A .4520E+02,.4613E+02,.5281E+02,.5179E+02,.4557E+02,.4381E+02,
     A .4893E+02,.5158E+02,.5112E+02,.5172E+02,.5466E+02,.5336E+02,
     A .4616E+02,.3774E+02,.3770E+02,.4279E+02,.4756E+02,.5019E+02,
     A .4897E+02,.4336E+02,.3863E+02,.3678E+02,.3358E+02,.3135E+02,
     A .3247E+02,.3332E+02,.3371E+02,.3257E+02,.3169E+02,.3132E+02,
     A .3031E+02,.2967E+02,.2939E+02,.3071E+02,.3196E+02,.3278E+02,
     A .3458E+02,.3601E+02,.3379E+02,.3297E+02,.3632E+02,.3810E+02,
     A .3584E+02,.3571E+02,.3915E+02,.4219E+02,.4431E+02,.4502E+02,
     A .4638E+02,.4656E+02,.4496E+02,.4377E+02,.4762E+02,.5507E+02,
     A .6000E+02,.6321E+02,.6355E+02,.5702E+02,.4982E+02,.5124E+02,
     A .5819E+02,.6478E+02,.6521E+02,.5445E+02,.4148E+02,.3561E+02,
     A .3825E+02,.4686E+02,.5251E+02,.4900E+02,.4530E+02,.4265E+02,
     A .4067E+02,.4205E+02,.4469E+02,.4582E+02,.4776E+02,.4998E+02,
     A .5138E+02,.5088E+02,.4978E+02,.5339E+02,.5975E+02,.6183E+02,
     A .5489E+02,.4067E+02,.3146E+02,.2964E+02,.3145E+02,.3496E+02,
     A .3831E+02,.3842E+02,.3760E+02,.4087E+02,.4579E+02,.4974E+02,
     A .5079E+02,.4371E+02,.3919E+02,.4605E+02,.5268E+02,.6102E+02,
     A .7170E+02,.7651E+02,.7724E+02,.7191E+02,.6556E+02,.6169E+02,
     A .5727E+02,.5152E+02,.4601E+02,.4277E+02,.3830E+02,.3398E+02,
     A .3204E+02,.3218E+02,.3264E+02,.3166E+02,.3598E+02,.4717E+02,
     A .5749E+02,.6353E+02,.6208E+02,.5782E+02,.5439E+02,.5206E+02,
     A .5279E+02,.5679E+02,.5704E+02,.5076E+02,.4108E+02,.3341E+02,
     A .3428E+02,.3703E+02,.4463E+02,.5686E+02,.6880E+02,.8081E+02,
     A .8040E+02,.7224E+02,.6028E+02,.4556E+02,.3779E+02,.3559E+02,
     A .3658E+02,.3864E+02,.4320E+02,.5169E+02,.5805E+02,.5874E+02,
     A .5672E+02,.5211E+02,.4462E+02,.3906E+02,.3659E+02,.3342E+02,
     A .3245E+02,.3661E+02,.4129E+02,.4495E+02,.5102E+02,.5767E+02,
     A .6080E+02,.6015E+02,.5224E+02,.3978E+02,.2969E+02,.2666E+02,
     A .3006E+02,.3410E+02,.3781E+02,.4075E+02,.4436E+02,.4899E+02,
     A .5024E+02,.4800E+02,.4097E+02,.3038E+02,.2437E+02,.2355E+02,
     A .2577E+02,.2826E+02,.2945E+02,.3171E+02,.3470E+02,.3390E+02,
     A .2915E+02,.2720E+02,.3283E+02,.4343E+02,.5085E+02,.5245E+02,
     A .5244E+02,.5410E+02,.5868E+02,.6292E+02,.6385E+02,.6234E+02,
     A .5612E+02,.4702E+02,.4446E+02,.4973E+02,.5886E+02,.6500E+02,
     A .6527E+02,.6527E+02,.6252E+02,.5225E+02,.3999E+02,.3314E+02,
     A .2977E+02,.3016E+02,.3640E+02,.4211E+02,.4566E+02,.5217E+02,
     A .5829E+02,.6006E+02,.5310E+02,.4339E+02,.3903E+02,.3995E+02,
     A .4225E+02,.4160E+02,.3993E+02,.4054E+02,.4223E+02,.4436E+02,
     A .4955E+02,.5463E+02,.5567E+02,.5838E+02,.6398E+02,.6850E+02,
     A .6907E+02,.5959E+02,.4447E+02,.3612E+02,.3557E+02,.3857E+02,
     A .4195E+02,.4519E+02,.5135E+02,.6169E+02,.6812E+02,.6524E+02,
     A .5908E+02,.4935E+02,.4094E+02,.3946E+02,.3972E+02,.3691E+02,
     A .3071E+02,.2488E+02,.2201E+02,.2188E+02,.2264E+02,.2603E+02,
     A .3064E+02,.3614E+02,.4422E+02,.5103E+02,.5464E+02,.5649E+02,
     A .5474E+02,.4962E+02,.4161E+02,.3201E+02,.3115E+02,.3796E+02,
     A .4816E+02,.6073E+02,.6393E+02,.5979E+02,.5937E+02,.5918E+02,
     A .5841E+02,.5687E+02,.5025E+02,.3607E+02,.2423E+02,.2218E+02,
     A .2806E+02,.3596E+02,.4117E+02,.4539E+02,.4349E+02,.3709E+02,
     A .3225E+02,.3697E+02,.4714E+02,.5169E+02,.5171E+02,.4937E+02,
     A .4821E+02,.4834E+02,.4750E+02,.3983E+02,.2988E+02,.2663E+02,
     A .2845E+02,.3121E+02,.3210E+02,.3182E+02,.3309E+02,.3960E+02,
     A .4773E+02,.5014E+02,.4939E+02,.4216E+02,.3243E+02,.2798E+02,
     A .2870E+02,.3498E+02,.4021E+02,.3894E+02,.3657E+02,.3881E+02,
     A .4332E+02,.4924E+02,.5398E+02,.5569E+02,.5874E+02,.6107E+02,
     A .6060E+02,.6162E+02,.6579E+02,.6845E+02,.6860E+02,.6998E+02,
     A .7344E+02,.7635E+02,.7934E+02,.8448E+02,.8334E+02,.7361E+02,
     A .6613E+02,.6352E+02,.6172E+02,.6320E+02,.6766E+02,.7159E+02,
     A .7672E+02,.8503E+02,.8976E+02,.8421E+02,.7356E+02,.6665E+02,
     A .6723E+02,.7282E+02,.7835E+02,.8036E+02,.7799E+02,.7257E+02,
     A .6709E+02,.6275E+02,.6335E+02,.6904E+02,.7032E+02,.6429E+02,
     A .5687E+02,.5525E+02,.5809E+02,.6120E+02,.6676E+02,.7443E+02,
     A .7295E+02,.6531E+02,.6086E+02,.5693E+02,.6208E+02,.7033E+02,
     A .7561E+02,.8117E+02,.8217E+02,.7985E+02,.7655E+02,.7022E+02,
     A .6507E+02,.6494E+02,.6600E+02,.6596E+02,.6510E+02,.6212E+02,
     A .5864E+02,.5575E+02,.5110E+02,.5150E+02,.5572E+02,.5586E+02,
     A .5398E+02,.5231E+02,.4972E+02,.4810E+02,.4740E+02,.4757E+02,
     A .4737E+02,.4643E+02,.4892E+02,.5118E+02,.5260E+02,.5313E+02,
     A .5132E+02,.4938E+02,.4988E+02,.5573E+02,.5906E+02,.5514E+02,
     A .5262E+02,.5231E+02,.5150E+02,.4948E+02,.4716E+02,.4511E+02,
     A .4343E+02,.4113E+02,.4041E+02,.4434E+02,.4908E+02,.5207E+02,
     A .5404E+02,.5398E+02,.5315E+02,.5303E+02,.5363E+02,.5641E+02,
     A .5882E+02,.5971E+02,.5984E+02,.5862E+02,.5767E+02,.5774E+02,
     A .5999E+02,.6488E+02,.6193E+02,.5195E+02,.4505E+02,.4671E+02,
     A .5766E+02,.6599E+02,.6748E+02,.6656E+02,.6539E+02,.6352E+02,
     A .6246E+02,.6158E+02,.5753E+02,.5393E+02,.5007E+02,.4454E+02,
     A .4635E+02,.5340E+02,.5773E+02,.6006E+02,.5696E+02,.4569E+02,
     A .3787E+02,.3370E+02,.2970E+02,.3556E+02,.4419E+02,.4797E+02,
     A .5384E+02,.5520E+02,.5048E+02,.4688E+02,.3982E+02,.3500E+02,
     A .3627E+02,.3782E+02,.3613E+02,.3271E+02,.3065E+02,.3129E+02,
     A .3444E+02,.4504E+02,.5436E+02,.5659E+02,.6184E+02,.7029E+02,
     A .7636E+02,.7678E+02,.7212E+02,.6880E+02,.7120E+02,.7327E+02,
     A .7309E+02,.7228E+02,.7420E+02,.7935E+02,.8143E+02,.7546E+02,
     A .6405E+02,.5367E+02,.5090E+02,.5738E+02,.6508E+02,.6884E+02,
     A .7218E+02,.7499E+02,.7413E+02,.6846E+02,.5938E+02,.5072E+02,
     A .4230E+02,.4184E+02,.5109E+02,.6022E+02,.6793E+02,.7309E+02,
     A .7383E+02,.6866E+02,.5937E+02,.5580E+02,.5639E+02,.5654E+02,
     A .5693E+02,.5249E+02,.4594E+02,.4203E+02,.4034E+02,.4001E+02,
     A .3940E+02,.3650E+02,.3635E+02,.3978E+02,.4084E+02,.4130E+02,
     A .4133E+02,.4040E+02,.4177E+02,.4708E+02,.5303E+02,.5460E+02,
     A .5229E+02,.4783E+02,.4019E+02,.3589E+02,.3557E+02,.3635E+02,
     A .3869E+02,.4314E+02,.4891E+02,.5143E+02,.4938E+02,.4481E+02,
     A .4062E+02,.3652E+02,.3321E+02,.3362E+02,.3473E+02,.3957E+02,
     A .5134E+02,.6117E+02,.6401E+02,.6208E+02,.5981E+02,.5771E+02,
     A .5771E+02,.5856E+02,.5836E+02,.5797E+02,.5477E+02,.5217E+02,
     A .5206E+02,.5422E+02,.5774E+02,.5868E+02,.5674E+02,.5567E+02,
     A .5914E+02,.6413E+02,.6217E+02,.5660E+02,.5513E+02,.5677E+02,
     A .6149E+02,.6724E+02,.7036E+02,.7068E+02,.6854E+02,.6657E+02,
     A .6270E+02,.5832E+02,.5696E+02,.5796E+02,.6065E+02,.6216E+02,
     A .5915E+02,.5143E+02,.4967E+02,.5482E+02,.6454E+02,.8056E+02,
     A .9225E+02,.9754E+02,.9172E+02,.8305E+02,.8766E+02,.9682E+02,
     A .9822E+02,.9203E+02,.8617E+02,.8607E+02,.8684E+02,.8364E+02,
     A .8022E+02,.8483E+02,.9440E+02,.1052E+03,.1082E+03,.1017E+03,
     A .9925E+02,.9792E+02,.9008E+02,.7938E+02,.7072E+02,.6644E+02,
     A .6972E+02,.7795E+02,.8791E+02,.9812E+02,.1133E+03,.1284E+03,
     A .1300E+03,.1240E+03,.1225E+03,.1235E+03,.1230E+03,.1277E+03,
     A .1392E+03,.1472E+03,.1424E+03,.1348E+03,.1374E+03,.1390E+03,
     A .1409E+03,.1512E+03,.1608E+03,.1566E+03,.1370E+03,.1138E+03,
     A .9358E+02,.7825E+02,.7040E+02,.7876E+02,.9311E+02,.1084E+03,
     A .1250E+03,.1332E+03,.1484E+03,.1633E+03,.1712E+03,.1872E+03,
     A .2046E+03,.1876E+03,.1446E+03,.1237E+03,.1272E+03,.1326E+03,
     A .1277E+03,.1092E+03,.8116E+02,.6270E+02,.6331E+02,.7469E+02,
     A .9700E+02,.1110E+03,.1212E+03,.1392E+03,.1655E+03,.1911E+03,
     A .1890E+03,.1556E+03,.1234E+03,.1154E+03,.1229E+03,.1153E+03,
     A .9265E+02,.8096E+02,.8601E+02,.1019E+03,.1169E+03,.1159E+03,
     A .1131E+03,.1293E+03,.1293E+03,.1041E+03,.7759E+02,.6047E+02,
     A .4956E+02,.4779E+02,.5559E+02,.8173E+02,.1150E+03,.1349E+03,
     A .1452E+03,.1490E+03,.1444E+03,.1369E+03,.1328E+03,.1232E+03,
     A .1032E+03,.8856E+02,.7862E+02,.6880E+02,.6403E+02,.7288E+02,
     A .9040E+02,.1058E+03,.1176E+03,.1213E+03,.1184E+03,.1033E+03,
     A .8193E+02,.6890E+02,.6368E+02,.6655E+02,.7402E+02,.7074E+02,
     A .8113E+02,.1019E+03,.1189E+03,.1336E+03,.1371E+03,.1292E+03,
     A .1053E+03,.8846E+02,.8997E+02,.9772E+02,.1168E+03,.1402E+03,
     A .1420E+03,.1301E+03,.1196E+03,.1203E+03,.1356E+03,.1446E+03,
     A .1449E+03,.1498E+03,.1542E+03,.1486E+03,.1323E+03,.1065E+03,
     A .9578E+02,.1028E+03,.1103E+03,.1087E+03,.9454E+02,.9051E+02,
     A .9969E+02,.1031E+03,.1004E+03,.7754E+02,.6147E+02,.7087E+02,
     A .8395E+02,.1002E+03,.1201E+03,.1417E+03,.1686E+03,.1868E+03,
     A .1871E+03,.1865E+03,.2019E+03,.2253E+03,.2403E+03,.2706E+03,
     A .3145E+03,.3228E+03,.2980E+03,.2773E+03,.2622E+03,.2516E+03,
     A .2483E+03,.2678E+03,.3063E+03,.3481E+03,.3783E+03,.3706E+03,
     A .3283E+03,.2950E+03,.2765E+03,.2616E+03,.2520E+03,.2497E+03,
     A .2391E+03,.2160E+03,.1945E+03,.1897E+03,.2065E+03,.2210E+03,
     A .2258E+03,.2259E+03,.2275E+03,.2395E+03,.2541E+03,.2631E+03,
     A .2649E+03,.2803E+03,.3078E+03,.3381E+03,.3476E+03,.3232E+03,
     A .3029E+03,.3085E+03,.3150E+03,.2979E+03,.2770E+03,.2820E+03,
     A .3016E+03,.3109E+03,.3081E+03,.2951E+03,.2720E+03,.2515E+03,
     A .2559E+03,.2618E+03,.2459E+03,.2224E+03,.2327E+03,.2708E+03,
     A .2928E+03,.2883E+03,.2502E+03,.2363E+03,.2750E+03,.2968E+03,
     A .2861E+03,.2673E+03,.2435E+03,.2279E+03,.2509E+03,.3004E+03,
     A .3190E+03,.2744E+03,.2223E+03,.1915E+03,.2240E+03,.2916E+03,
     A .3291E+03,.3434E+03,.3412E+03,.3344E+03,.3177E+03,.2845E+03,
     A .2653E+03,.2519E+03,.2317E+03,.2198E+03,.2203E+03,.2298E+03,
     A .2457E+03,.2680E+03,.2840E+03,.2892E+03,.2873E+03,.2790E+03,
     A .2690E+03,.2585E+03,.2509E+03,.2536E+03,.2665E+03,.2843E+03,
     A .2979E+03,.3004E+03,.2910E+03,.2761E+03,.2677E+03,.2513E+03,
     A .2303E+03,.2157E+03,.2202E+03,.2363E+03,.2609E+03,.2760E+03,
     A .2706E+03,.2714E+03,.2873E+03,.2948E+03,.2857E+03,.2568E+03,
     A .2093E+03,.1929E+03,.2116E+03,.2484E+03,.2619E+03,.2448E+03,
     A .2358E+03,.2507E+03,.2976E+03,.3629E+03,.3780E+03,.3287E+03,
     A .2768E+03,.2674E+03,.2853E+03,.3107E+03,.3097E+03,.3251E+03,
     A .3576E+03,.3589E+03,.3346E+03,.2889E+03,.2596E+03,.2616E+03,
     A .2723E+03,.2782E+03,.2664E+03,.2610E+03,.2768E+03,.2854E+03,
     A .2806E+03,.2698E+03,.2612E+03,.2826E+03,.3302E+03,.3452E+03,
     A .3030E+03,.2491E+03,.2073E+03,.1990E+03,.2165E+03,.2384E+03,
     A .2715E+03,.2814E+03,.2523E+03,.2160E+03,.1589E+03,.1181E+03,
     A .1192E+03,.1270E+03,.1274E+03,.1336E+03,.1644E+03,.1948E+03,
     A .2206E+03,.2122E+03,.1939E+03,.1912E+03,.1953E+03,.2058E+03,
     A .2182E+03,.2255E+03,.2251E+03,.2158E+03,.1965E+03,.2001E+03,
     A .2309E+03,.3166E+03,.3702E+03,.3342E+03,.2900E+03,.2762E+03,
     A .2781E+03,.2815E+03,.2764E+03,.2463E+03,.2128E+03,.1983E+03,
     A .2036E+03,.2124E+03,.2233E+03,.2127E+03,.1695E+03,.1425E+03,
     A .1442E+03,.1598E+03,.1784E+03,.1670E+03,.1343E+03,.1359E+03,
     A .1542E+03,.1648E+03,.1628E+03,.1480E+03,.1269E+03,.1067E+03,
     A .1013E+03,.1110E+03,.1316E+03,.1659E+03,.1842E+03,.1587E+03,
     A .1251E+03,.1095E+03,.1184E+03,.1456E+03,.1555E+03,.1326E+03,
     A .9572E+02,.9431E+02,.1157E+03,.1599E+03,.1989E+03,.2351E+03,
     A .2655E+03,.2560E+03,.2288E+03,.2058E+03,.1855E+03,.1671E+03,
     A .1310E+03,.1089E+03,.1272E+03,.1590E+03,.2003E+03,.2596E+03,
     A .3242E+03,.3423E+03,.3119E+03,.2898E+03,.2883E+03,.2873E+03,
     A .2619E+03,.2116E+03,.2109E+03,.2270E+03,.2530E+03,.2899E+03,
     A .3246E+03,.3405E+03,.3251E+03,.2990E+03,.2744E+03,.2658E+03,
     A .2664E+03,.2609E+03,.2390E+03,.2307E+03,.2431E+03,.2716E+03,
     A .3036E+03,.3129E+03,.2879E+03,.2613E+03,.2670E+03,.2881E+03,
     A .3076E+03,.3176E+03,.3098E+03,.2886E+03,.2654E+03,.2474E+03,
     A .2372E+03,.2357E+03,.2222E+03,.1967E+03,.1841E+03,.1808E+03,
     A .1745E+03,.1711E+03,.1732E+03,.1751E+03,.1820E+03,.1870E+03,
     A .1817E+03,.1730E+03,.1711E+03,.1795E+03,.1873E+03,.1891E+03,
     A .1910E+03,.1910E+03,.1885E+03,.1802E+03,.1589E+03,.1426E+03,
     A .1408E+03,.1362E+03,.1255E+03,.1181E+03,.1109E+03,.1044E+03,
     A .1002E+03,.9661E+02,.9089E+02,.8352E+02,.7073E+02,.5916E+02,
     A .7047E+02,.9171E+02,.9386E+02,.8182E+02,.6926E+02,.6167E+02,
     A .6164E+02,.6576E+02,.7122E+02,.7687E+02,.8358E+02,.8562E+02,
     A .8306E+02,.6994E+02,.6656E+02,.7956E+02,.8586E+02,.8225E+02,
     A .7571E+02,.7207E+02,.7545E+02,.8706E+02,.1006E+03,.1111E+03,
     A .1175E+03,.1279E+03,.1382E+03,.1488E+03,.1607E+03,.1723E+03,
     A .1823E+03,.1846E+03,.1831E+03,.1834E+03,.1827E+03,.1825E+03,
     A .1858E+03,.1800E+03,.2196E+03,.2548E+03,.2773E+03,.2866E+03,
     A .2819E+03,.2750E+03,.2669E+03,.2587E+03,.2649E+03,.2839E+03,
     A .2993E+03,.3103E+03,.3213E+03,.3294E+03,.3295E+03,.3240E+03,
     A .3200E+03,.3177E+03,.3119E+03,.3348E+03,.3685E+03,.3828E+03,
     A .3645E+03,.3176E+03,.2930E+03,.3074E+03,.3291E+03,.3391E+03,
     A .3409E+03,.3408E+03,.3539E+03,.3770E+03,.3939E+03,.3951E+03,
     A .3935E+03,.3930E+03,.3757E+03,.3651E+03,.3745E+03,.3859E+03,
     A .3897E+03,.3855E+03,.3577E+03,.3072E+03,.2737E+03,.2676E+03,
     A .3155E+03,.3775E+03,.4021E+03,.3995E+03,.3931E+03,.3863E+03,
     A .3575E+03,.3336E+03,.3235E+03,.3353E+03,.3630E+03,.3704E+03,
     A .3504E+03,.2953E+03,.2649E+03,.2943E+03,.3242E+03,.3163E+03,
     A .2855E+03,.2595E+03,.2453E+03,.2329E+03,.2112E+03,.1897E+03,
     A .1696E+03,.1453E+03,.1250E+03,.1073E+03,.9557E+02,.8250E+02,
     A .6286E+02,.5116E+02,.5302E+02,.7098E+02,.9571E+02,.1193E+03,
     A .1354E+03,.1494E+03,.1682E+03,.1870E+03,.2051E+03,.2231E+03,
     A .2352E+03,.2467E+03,.2661E+03,.3030E+03,.3617E+03,.3735E+03,
     A .3717E+03,.3793E+03,.3810E+03,.3781E+03,.3742E+03,.3776E+03,
     A .4033E+03,.4246E+03,.4108E+03,.4053E+03,.4194E+03,.4161E+03,
     A .3851E+03,.3612E+03,.3547E+03,.3478E+03,.3350E+03,.3235E+03,
     A .3550E+03,.4212E+03,.4821E+03,.5118E+03,.5101E+03,.4948E+03,
     A .4727E+03,.4451E+03,.4159E+03,.4053E+03,.4190E+03,.4222E+03,
     A .3979E+03,.3727E+03,.3714E+03,.4111E+03,.4126E+03,.3689E+03,
     A .3294E+03,.3017E+03,.2876E+03,.2594E+03,.2061E+03,.1477E+03,
     A .1299E+03,.1635E+03,.2204E+03,.2889E+03,.3385E+03,.4040E+03,
     A .4645E+03,.4864E+03,.4814E+03,.4654E+03,.4558E+03,.4499E+03,
     A .4416E+03,.4432E+03,.4544E+03,.4467E+03,.4457E+03,.4503E+03,
     A .4440E+03,.4516E+03,.4631E+03,.4656E+03,.4561E+03,.4541E+03,
     A .4608E+03,.4585E+03,.4901E+03,.5528E+03,.5924E+03,.6015E+03,
     A .6000E+03,.5916E+03,.5836E+03,.5858E+03,.5967E+03,.6232E+03,
     A .6872E+03,.7497E+03,.7735E+03,.6963E+03,.6229E+03,.6305E+03,
     A .6754E+03,.7332E+03,.7521E+03,.7541E+03,.7556E+03,.7386E+03,
     A .7010E+03,.6685E+03,.6650E+03,.6547E+03,.6338E+03,.6228E+03,
     A .6220E+03,.6255E+03,.6258E+03,.6189E+03,.6190E+03,.6247E+03,
     A .6379E+03,.7039E+03,.7719E+03,.7539E+03,.6959E+03,.6517E+03,
     A .5990E+03,.5727E+03,.5685E+03,.6011E+03,.6247E+03,.6135E+03,
     A .5725E+03,.5719E+03,.6088E+03,.6286E+03,.6226E+03,.6144E+03,
     A .6264E+03,.6328E+03,.6384E+03,.6360E+03,.5826E+03,.5284E+03,
     A .5132E+03,.5264E+03,.5372E+03,.5230E+03,.4884E+03,.4585E+03,
     A .4931E+03,.5422E+03,.5559E+03,.4955E+03,.4040E+03,.4824E+03,
     A .6277E+03,.6846E+03,.6937E+03,.7139E+03,.7336E+03,.7081E+03,
     A .6690E+03,.6356E+03,.6337E+03,.6593E+03,.7323E+03,.7341E+03,
     A .6070E+03,.4500E+03,.3654E+03,.3827E+03,.4296E+03,.4666E+03,
     A .4894E+03,.5111E+03,.5735E+03,.5737E+03,.5257E+03,.5026E+03,
     A .5457E+03,.6131E+03,.6603E+03,.6436E+03,.5990E+03,.5561E+03,
     A .5433E+03,.5852E+03,.6995E+03,.7153E+03,.5625E+03,.4262E+03,
     A .3585E+03,.3265E+03,.3315E+03,.3821E+03,.4704E+03,.5371E+03,
     A .6074E+03,.6758E+03,.7216E+03,.7510E+03,.6842E+03,.5373E+03,
     A .4673E+03,.4976E+03,.5826E+03,.6730E+03,.6637E+03,.5782E+03,
     A .5224E+03,.5502E+03,.6214E+03,.7022E+03,.7251E+03,.6708E+03,
     A .6142E+03,.5926E+03,.5933E+03,.6279E+03,.6956E+03,.7497E+03,
     A .7792E+03,.7400E+03,.6747E+03,.5669E+03,.4528E+03,.4062E+03,
     A .3770E+03,.3170E+03,.2849E+03,.3256E+03,.4630E+03,.5246E+03,
     A .4898E+03,.4208E+03,.3549E+03,.3966E+03,.4756E+03,.5367E+03,
     A .5470E+03,.4656E+03,.3822E+03,.3836E+03,.5090E+03,.6174E+03,
     A .6692E+03,.6678E+03,.6437E+03,.6474E+03,.7495E+03,.8532E+03,
     A .8393E+03,.7648E+03,.6714E+03,.6171E+03,.6006E+03,.5756E+03,
     A .5227E+03,.4177E+03,.4022E+03,.4246E+03,.3966E+03,.3327E+03,
     A .3068E+03,.3067E+03,.3126E+03,.3207E+03,.3557E+03,.4128E+03,
     A .4281E+03,.4311E+03,.4642E+03,.5170E+03,.5732E+03,.6277E+03,
     A .6790E+03,.7162E+03,.7304E+03,.6858E+03,.6026E+03,.5851E+03,
     A .5794E+03,.5406E+03,.4595E+03,.3722E+03,.3411E+03,.3872E+03,
     A .5542E+03,.6190E+03,.6707E+03,.6961E+03,.6849E+03,.6545E+03,
     A .5530E+03,.4204E+03,.3275E+03,.2749E+03,.2824E+03,.3119E+03,
     A .3135E+03,.3000E+03,.2913E+03,.3321E+03,.4368E+03,.5262E+03,
     A .5727E+03,.5937E+03,.5603E+03,.5021E+03,.4482E+03,.3793E+03,
     A .3481E+03,.3612E+03,.3766E+03,.4064E+03,.4555E+03,.5534E+03,
     A .6326E+03,.5669E+03,.5126E+03,.5450E+03,.5766E+03,.5936E+03,
     A .5935E+03,.5807E+03,.5822E+03,.5846E+03,.5872E+03,.5616E+03,
     A .5135E+03,.4702E+03,.4330E+03,.4009E+03,.3762E+03,.3504E+03,
     A .2771E+03,.2039E+03,.1893E+03,.2506E+03,.4046E+03,.5678E+03,
     A .6379E+03,.6556E+03,.6123E+03,.5842E+03,.5211E+03,.4488E+03,
     A .4180E+03,.4283E+03,.5383E+03,.7271E+03,.8315E+03,.8111E+03,
     A .7513E+03,.7047E+03,.6775E+03,.6631E+03,.6043E+03,.5850E+03,
     A .7033E+03,.8131E+03,.8318E+03,.8240E+03,.7936E+03,.7725E+03,
     A .7894E+03,.8333E+03,.8060E+03,.6765E+03,.5387E+03,.4080E+03,
     A .3984E+03,.5187E+03,.6402E+03,.6969E+03,.6996E+03,.6609E+03,
     A .6021E+03,.5369E+03,.4912E+03,.5516E+03,.6423E+03,.7014E+03,
     A .7414E+03,.7609E+03,.7675E+03,.7874E+03,.8067E+03,.7091E+03,
     A .5473E+03,.4406E+03,.4334E+03,.5362E+03,.6434E+03,.7570E+03,
     A .7709E+03,.7099E+03,.6905E+03,.7672E+03,.9233E+03,.9538E+03,
     A .8485E+03,.7441E+03,.6730E+03,.6070E+03,.6469E+03,.7221E+03,
     A .7272E+03,.6678E+03,.5636E+03,.4541E+03,.3773E+03,.3533E+03,
     A .3457E+03,.3942E+03,.4645E+03,.5769E+03,.6954E+03,.7935E+03,
     A .8060E+03,.7956E+03,.7577E+03,.6998E+03,.6165E+03,.5748E+03,
     A .6072E+03,.6508E+03,.6498E+03,.5855E+03,.4962E+03,.4629E+03,
     A .5252E+03,.6396E+03,.7167E+03,.7756E+03,.8016E+03,.8020E+03,
     A .7916E+03,.7606E+03,.7167E+03,.6829E+03,.6951E+03,.7394E+03,
     A .7751E+03,.7315E+03,.6272E+03,.5763E+03,.6327E+03,.7445E+03,
     A .7845E+03,.7548E+03,.7015E+03,.6616E+03,.6987E+03,.7798E+03,
     A .8141E+03,.7815E+03,.6927E+03,.6164E+03,.5726E+03,.5773E+03,
     A .5903E+03,.6337E+03,.7153E+03,.8030E+03,.8682E+03,.9011E+03,
     A .8877E+03,.7664E+03,.7037E+03,.6370E+03,.6181E+03,.6627E+03,
     A .6951E+03,.6768E+03,.6521E+03,.6004E+03,.5338E+03,.5089E+03,
     A .4874E+03,.4366E+03,.3731E+03,.4023E+03,.4606E+03,.5417E+03,
     A .6954E+03,.7336E+03,.6850E+03,.5909E+03,.4785E+03,.4843E+03,
     A .5764E+03,.6433E+03,.5839E+03,.4738E+03,.3914E+03,.3439E+03,
     A .3335E+03,.3537E+03,.4017E+03,.4859E+03,.6303E+03,.7653E+03,
     A .8363E+03,.8541E+03,.8358E+03,.7996E+03,.8264E+03,.8532E+03,
     A .8756E+03,.9034E+03,.9492E+03,.9970E+03,.9868E+03,.9312E+03,
     A .8579E+03,.8234E+03,.8498E+03,.8795E+03,.9259E+03,.9776E+03,
     A .9689E+03,.8908E+03,.8329E+03,.8445E+03,.8903E+03,.9301E+03,
     A .8919E+03,.8256E+03,.7547E+03,.7152E+03,.7046E+03,.6996E+03,
     A .6914E+03,.6571E+03,.6511E+03,.6631E+03,.6625E+03,.6527E+03,
     A .6515E+03,.7076E+03,.7662E+03,.8387E+03,.8815E+03,.8146E+03,
     A .7287E+03,.6780E+03,.6783E+03,.7299E+03,.7763E+03,.7984E+03,
     A .7855E+03,.7676E+03,.7904E+03,.8086E+03,.7830E+03,.7540E+03,
     A .7596E+03,.7564E+03,.7484E+03,.7479E+03,.7511E+03,.7333E+03,
     A .6948E+03,.6823E+03,.7379E+03,.7834E+03,.7746E+03,.7714E+03,
     A .8147E+03,.8874E+03,.9323E+03,.9641E+03,.9493E+03,.9078E+03,
     A .8766E+03,.8816E+03,.9114E+03,.9307E+03,.9234E+03,.8787E+03,
     A .7428E+03,.5448E+03,.4210E+03,.3742E+03,.4127E+03,.5326E+03,
     A .6258E+03,.6819E+03,.6885E+03,.6889E+03,.7094E+03,.7785E+03,
     A .8748E+03,.8797E+03,.8261E+03,.8080E+03,.8104E+03,.8121E+03,
     A .8037E+03,.8060E+03,.8234E+03,.8384E+03,.8166E+03,.7640E+03,
     A .7568E+03,.7771E+03,.7648E+03,.7284E+03,.6115E+03,.5317E+03,
     A .4587E+03,.4409E+03,.4501E+03,.4896E+03,.5145E+03,.5106E+03,
     A .5085E+03,.5346E+03,.5978E+03,.6647E+03,.7431E+03,.7684E+03,
     A .7341E+03,.7050E+03,.7146E+03,.7116E+03,.7082E+03,.7437E+03,
     A .8141E+03,.8201E+03,.7600E+03,.7197E+03,.7185E+03,.7547E+03,
     A .8505E+03,.9173E+03,.9928E+03,.1022E+04,.1020E+04,.1027E+04,
     A .1012E+04,.9744E+03,.9388E+03,.8927E+03,.8767E+03,.9222E+03,
     A .9943E+03,.1045E+04,.1035E+04,.9307E+03,.8190E+03,.6919E+03,
     A .6118E+03,.5615E+03,.5636E+03,.5802E+03,.5687E+03,.5601E+03,
     A .5822E+03,.6111E+03,.6302E+03,.6618E+03,.6777E+03,.6892E+03,
     A .7299E+03,.7659E+03,.7761E+03,.7660E+03,.7772E+03,.7883E+03,
     A .8022E+03,.8795E+03,.9949E+03,.1021E+04,.9417E+03,.8277E+03,
     A .7197E+03,.6464E+03,.5527E+03,.5001E+03,.5484E+03,.6125E+03,
     A .6740E+03,.7253E+03,.7492E+03,.7518E+03,.7417E+03,.7222E+03,
     A .7069E+03,.7413E+03,.8504E+03,.1006E+04,.9532E+03,.8549E+03,
     A .7946E+03,.7727E+03,.8435E+03,.9094E+03,.9116E+03,.8770E+03,
     A .8682E+03,.9100E+03,.9887E+03,.1022E+04,.9863E+03,.9605E+03,
     A .9688E+03,.9894E+03,.9988E+03,.9911E+03,.9317E+03,.8242E+03,
     A .7832E+03,.7657E+03,.7352E+03,.7008E+03,.6493E+03,.6093E+03,
     A .6268E+03,.6514E+03,.6378E+03,.6100E+03,.5905E+03,.6020E+03,
     A .6653E+03,.7848E+03,.7866E+03,.7564E+03,.7297E+03,.7569E+03,
     A .8617E+03,.9601E+03,.9717E+03,.9339E+03,.9360E+03,.9902E+03,
     A .9955E+03,.8989E+03,.7768E+03,.7148E+03,.7135E+03,.7919E+03,
     A .8970E+03,.9020E+03,.8474E+03,.7424E+03,.6921E+03,.6927E+03,
     A .7054E+03,.6768E+03,.6427E+03,.6174E+03,.5922E+03,.5729E+03,
     A .6001E+03,.6710E+03,.7225E+03,.7871E+03,.8554E+03,.8541E+03,
     A .8195E+03,.7831E+03,.6778E+03,.6693E+03,.7144E+03,.7347E+03,
     A .7102E+03,.6685E+03,.6880E+03,.7648E+03,.8396E+03,.8262E+03,
     A .7613E+03,.7261E+03,.8088E+03,.9148E+03,.9372E+03,.8909E+03,
     A .8671E+03,.8773E+03,.8708E+03,.8372E+03,.8252E+03,.8612E+03,
     A .9233E+03,.9651E+03,.9942E+03,.9526E+03,.8797E+03,.8070E+03,
     A .7601E+03,.7416E+03,.7811E+03,.8155E+03,.8133E+03,.8052E+03,
     A .8102E+03,.8344E+03,.8647E+03,.8726E+03,.8848E+03,.8821E+03,
     A .8403E+03,.8852E+03,.1025E+04,.1101E+04,.1120E+04,.1106E+04,
     A .1069E+04,.1036E+04,.1045E+04,.1086E+04,.1147E+04,.1157E+04,
     A .1144E+04,.1144E+04,.1099E+04,.1065E+04,.1098E+04,.1109E+04,
     A .1037E+04,.9687E+03,.9103E+03,.8830E+03,.9174E+03,.9925E+03,
     A .1099E+04,.1135E+04,.1127E+04,.1132E+04,.1126E+04,.1103E+04,
     A .1088E+04,.1066E+04,.1004E+04,.9318E+03,.9172E+03,.9605E+03,
     A .1021E+04,.1080E+04,.1049E+04,.1038E+04,.1078E+04,.1113E+04,
     A .1092E+04,.1016E+04,.9788E+03,.9677E+03,.1009E+04,.1011E+04,
     A .9856E+03,.9415E+03,.9279E+03,.9628E+03,.9993E+03,.1021E+04,
     A .1001E+04,.9411E+03,.9174E+03,.9399E+03,.9995E+03,.1034E+04,
     A .1013E+04,.9906E+03,.9713E+03,.9159E+03,.8082E+03,.7532E+03,
     A .7633E+03,.8201E+03,.8779E+03,.9591E+03,.1080E+04,.1157E+04,
     A .1160E+04,.1142E+04,.1118E+04,.1086E+04,.1101E+04,.1143E+04,
     A .1184E+04,.1248E+04,.1236E+04,.1166E+04,.1140E+04,.1155E+04,
     A .1186E+04,.1222E+04,.1219E+04,.1161E+04,.1142E+04,.1134E+04,
     A .1234E+04,.1313E+04,.1373E+04,.1357E+04,.1319E+04,.1143E+04,
     A .1032E+04,.9970E+03,.1046E+04,.1092E+04,.1085E+04,.1020E+04,
     A .8963E+03,.7890E+03,.7833E+03,.8250E+03,.8920E+03,.9770E+03,
     A .1045E+04,.1077E+04,.1085E+04,.1068E+04,.1051E+04,.1081E+04,
     A .1086E+04,.1053E+04,.1005E+04,.1005E+04,.9691E+03,.9290E+03,
     A .8801E+03,.8889E+03,.9533E+03,.1011E+04,.1048E+04,.1068E+04,
     A .1091E+04,.1106E+04,.1106E+04,.1098E+04,.1060E+04,.1047E+04,
     A .1037E+04,.1052E+04,.1063E+04,.9932E+03,.8901E+03,.8298E+03,
     A .8172E+03,.8653E+03,.9427E+03,.1053E+04,.1131E+04,.1173E+04,
     A .1119E+04,.1034E+04,.1019E+04,.1018E+04,.9991E+03,.8909E+03,
     A .8653E+03,.8979E+03,.1008E+04,.1069E+04,.1042E+04,.1031E+04,
     A .1061E+04,.1092E+04,.1056E+04,.1028E+04,.9854E+03,.9114E+03,
     A .8247E+03,.7559E+03,.7423E+03,.7618E+03,.8294E+03,.9280E+03,
     A .9957E+03,.1029E+04,.1017E+04,.9955E+03,.9890E+03,.1026E+04,
     A .9719E+03,.9259E+03,.9338E+03,.9936E+03,.1064E+04,.1105E+04,
     A .1125E+04,.1116E+04,.1132E+04,.1130E+04,.1109E+04,.1102E+04,
     A .1090E+04,.1048E+04,.9006E+03,.7469E+03,.7069E+03,.7740E+03,
     A .9309E+03,.1058E+04,.1102E+04,.1118E+04,.1135E+04,.1131E+04,
     A .1090E+04,.1052E+04,.1051E+04,.1099E+04,.1167E+04,.1190E+04,
     A .1168E+04,.1159E+04,.1148E+04,.1062E+04,.9271E+03,.7499E+03,
     A .6573E+03,.5781E+03,.5270E+03,.4542E+03,.4842E+03,.6016E+03,
     A .7689E+03,.9158E+03,.9853E+03,.1024E+04,.1094E+04,.1130E+04,
     A .1116E+04,.1060E+04,.9897E+03,.9055E+03,.8300E+03,.8093E+03,
     A .7777E+03,.7122E+03,.6774E+03,.6996E+03,.7345E+03,.7602E+03,
     A .7384E+03,.7072E+03,.6697E+03,.7399E+03,.8167E+03,.8490E+03,
     A .8735E+03,.9491E+03,.1057E+04,.1119E+04,.1155E+04,.1168E+04,
     A .1165E+04,.1150E+04,.1052E+04,.9005E+03,.8347E+03,.7346E+03,
     A .6491E+03,.6813E+03,.8528E+03,.9684E+03,.1008E+04,.9998E+03,
     A .9364E+03,.9200E+03,.9668E+03,.1034E+04,.1039E+04,.1096E+04,
     A .1167E+04,.1142E+04,.9953E+03,.9278E+03,.9171E+03,.9548E+03,
     A .1061E+04,.1168E+04,.1183E+04,.1060E+04,.9662E+03,.9561E+03,
     A .8532E+03,.7725E+03,.7327E+03,.7798E+03,.8634E+03,.9591E+03,
     A .1012E+04,.1029E+04,.1049E+04,.1083E+04,.1102E+04,.1138E+04,
     A .1168E+04,.1141E+04,.1134E+04,.1203E+04,.1312E+04,.1338E+04,
     A .1267E+04,.1164E+04,.1110E+04,.1098E+04,.1116E+04,.1103E+04,
     A .1043E+04,.1002E+04,.1045E+04,.1083E+04,.1075E+04,.1003E+04,
     A .9351E+03,.9303E+03,.9755E+03,.9893E+03,.9586E+03,.9581E+03,
     A .9925E+03,.1052E+04,.1097E+04,.1154E+04,.1113E+04,.9980E+03,
     A .8584E+03,.7977E+03,.7189E+03,.7008E+03,.8253E+03,.9513E+03,
     A .1000E+04,.1028E+04,.1034E+04,.1021E+04,.1034E+04,.1116E+04,
     A .1195E+04,.1247E+04,.1244E+04,.1197E+04,.1180E+04,.1160E+04,
     A .1078E+04,.9836E+03,.9300E+03,.8984E+03,.8900E+03,.9501E+03,
     A .1064E+04,.1162E+04,.1069E+04,.9925E+03,.9774E+03,.9924E+03,
     A .1022E+04,.1108E+04,.1230E+04,.1287E+04,.1308E+04,.1318E+04,
     A .1273E+04,.1178E+04,.1098E+04,.1039E+04,.9480E+03,.8249E+03,
     A .9924E+03,.1173E+04,.1219E+04,.1210E+04,.1150E+04,.1050E+04,
     A .9768E+03,.9545E+03,.9695E+03,.1004E+04,.9750E+03,.8814E+03,
     A .7129E+03,.5449E+03,.5008E+03,.5372E+03,.5945E+03,.6527E+03,
     A .6720E+03,.6364E+03,.6594E+03,.7645E+03,.8510E+03,.8639E+03,
     A .8217E+03,.8025E+03,.8334E+03,.9330E+03,.1008E+04,.1017E+04,
     A .1013E+04,.1047E+04,.1114E+04,.1186E+04,.1230E+04,.1183E+04,
     A .1042E+04,.8819E+03,.8776E+03,.9541E+03,.1071E+04,.1194E+04,
     A .1207E+04,.1184E+04,.1138E+04,.1095E+04,.1094E+04,.1043E+04,
     A .9420E+03,.9052E+03,.9003E+03,.8856E+03,.8506E+03,.8381E+03,
     A .7898E+03,.7798E+03,.8186E+03,.9291E+03,.1038E+04,.1128E+04,
     A .1194E+04,.1207E+04,.1132E+04,.9761E+03,.8859E+03,.9411E+03,
     A .1045E+04,.1140E+04,.1176E+04,.1150E+04,.1107E+04,.1119E+04,
     A .1187E+04,.1258E+04,.1187E+04,.1052E+04,.9768E+03,.9674E+03,
     A .9769E+03,.9648E+03,.9330E+03,.8765E+03,.7811E+03,.6947E+03,
     A .6887E+03,.6963E+03,.7940E+03,.9251E+03,.1009E+04,.1053E+04,
     A .1068E+04,.1065E+04,.1050E+04,.1042E+04,.1033E+04,.1049E+04,
     A .1015E+04,.9457E+03,.8814E+03,.8770E+03,.8988E+03,.8993E+03,
     A .8990E+03,.9182E+03,.9992E+03,.1122E+04,.1215E+04,.1196E+04,
     A .1068E+04,.9415E+03,.8856E+03,.8354E+03,.7622E+03,.7254E+03,
     A .8053E+03,.9193E+03,.8713E+03,.8091E+03,.8466E+03,.9209E+03,
     A .1008E+04,.1073E+04,.9992E+03,.9288E+03,.8571E+03,.8089E+03,
     A .7568E+03,.7555E+03,.7987E+03,.1022E+04,.1190E+04,.1253E+04,
     A .1238E+04,.1204E+04,.1206E+04,.1188E+04,.1151E+04,.1180E+04,
     A .1279E+04,.1357E+04,.1327E+04,.1208E+04,.1119E+04,.1078E+04,
     A .1047E+04,.1043E+04,.1095E+04,.1151E+04,.1177E+04,.1209E+04,
     A .1176E+04,.1055E+04,.9336E+03,.9402E+03,.1020E+04,.1047E+04,
     A .1038E+04,.1002E+04,.9280E+03,.8527E+03,.8297E+03,.8027E+03,
     A .8371E+03,.9577E+03,.1085E+04,.1147E+04,.1258E+04,.1328E+04,
     A .1302E+04,.1287E+04,.1253E+04,.1181E+04,.1114E+04,.1040E+04,
     A .9523E+03,.1008E+04,.1042E+04,.1008E+04,.9261E+03,.8276E+03,
     A .7536E+03,.7527E+03,.7993E+03,.7549E+03,.6793E+03,.7109E+03,
     A .8622E+03,.1076E+04,.1168E+04,.1190E+04,.1175E+04,.1159E+04,
     A .1176E+04,.1205E+04,.1145E+04,.1052E+04,.9855E+03,.9669E+03,
     A .9992E+03,.1137E+04,.1254E+04,.1279E+04,.1245E+04,.1212E+04,
     A .1178E+04,.1166E+04,.1143E+04,.1111E+04,.1147E+04,.1259E+04,
     A .1412E+04,.1369E+04,.1238E+04,.1155E+04,.1110E+04,.1079E+04,
     A .1086E+04,.1156E+04,.1241E+04,.1311E+04,.1309E+04,.1263E+04,
     A .1227E+04,.1233E+04,.1281E+04,.1303E+04,.1229E+04,.1141E+04,
     A .1174E+04,.1244E+04,.1291E+04,.1319E+04,.1301E+04,.1269E+04,
     A .1246E+04,.1228E+04,.1168E+04,.1122E+04,.1076E+04,.9952E+03,
     A .9487E+03,.9355E+03,.1040E+04,.1182E+04,.1170E+04,.1204E+04,
     A .1207E+04,.1130E+04,.1056E+04,.1035E+04,.1076E+04,.1129E+04,
     A .1153E+04,.1165E+04,.1191E+04,.1218E+04,.1281E+04,.1390E+04,
     A .1440E+04,.1398E+04,.1200E+04,.1050E+04,.7454E+03,.6720E+03,
     A .6744E+03,.8102E+03,.9951E+03,.1118E+04,.1145E+04,.1041E+04,
     A .6716E+03,.6455E+03,.6765E+03,.8132E+03,.8618E+03,.8894E+03,
     A .1047E+04,.1099E+04,.1100E+04,.1089E+04,.1094E+04,.1088E+04,
     A .1097E+04,.1146E+04,.1221E+04,.1260E+04,.1203E+04,.1108E+04,
     A .1008E+04,.9660E+03,.9809E+03,.9428E+03,.7967E+03,.5580E+03,
     A .5084E+03,.5801E+03,.6627E+03,.7615E+03,.8821E+03,.9240E+03,
     A .8525E+03,.7750E+03,.7106E+03,.6571E+03,.6236E+03,.6295E+03,
     A .6817E+03,.7762E+03,.8384E+03,.8663E+03,.8596E+03,.9790E+03,
     A .1068E+04,.1212E+04,.1314E+04,.1333E+04,.1325E+04,.1208E+04,
     A .1079E+04,.1123E+04,.1152E+04,.1196E+04,.1249E+04,.1323E+04,
     A .1390E+04,.1406E+04,.1389E+04,.1327E+04,.1324E+04,.1290E+04,
     A .1316E+04,.1382E+04,.1469E+04,.1494E+04,.1426E+04,.1191E+04,
     A .1075E+04,.1059E+04,.1044E+04,.1053E+04,.1058E+04,.1010E+04,
     A .9457E+03,.1022E+04,.1079E+04,.1139E+04,.1161E+04,.1078E+04,
     A .8932E+03,.7463E+03,.7295E+03,.7597E+03,.8151E+03,.9494E+03,
     A .1102E+04,.1188E+04,.1167E+04,.1083E+04,.1017E+04,.9977E+03,
     A .1052E+04,.1159E+04,.1253E+04,.1286E+04,.1281E+04,.1204E+04,
     A .1054E+04,.8409E+03,.6454E+03,.5571E+03,.5852E+03,.7230E+03,
     A .8864E+03,.9705E+03,.9936E+03,.1125E+04,.1260E+04,.1277E+04,
     A .1202E+04,.1136E+04,.1125E+04,.1172E+04,.1294E+04,.1440E+04,
     A .1521E+04,.1544E+04,.1512E+04,.1481E+04,.1467E+04,.1487E+04,
     A .1437E+04,.1290E+04,.1006E+04,.8299E+03,.7005E+03,.6814E+03,
     A .8116E+03,.9694E+03,.1085E+04,.1062E+04,.9785E+03,.1040E+04,
     A .1114E+04,.1128E+04,.1134E+04,.1145E+04,.1151E+04,.1163E+04,
     A .1239E+04,.1303E+04,.1346E+04,.1353E+04,.1362E+04,.1319E+04,
     A .1308E+04,.1272E+04,.1233E+04,.1182E+04,.1157E+04,.1131E+04,
     A .1060E+04,.9930E+03,.1005E+04,.1086E+04,.1123E+04,.1044E+04,
     A .8932E+03,.8424E+03,.9081E+03,.9603E+03,.1023E+04,.1075E+04,
     A .1137E+04,.1139E+04,.1175E+04,.1294E+04,.1425E+04,.1438E+04,
     A .1425E+04,.1419E+04,.1410E+04,.1385E+04,.1347E+04,.1385E+04,
     A .1434E+04,.1484E+04,.1444E+04,.1497E+04,.1561E+04,.1592E+04,
     A .1571E+04,.1529E+04,.1513E+04,.1535E+04,.1568E+04,.1534E+04,
     A .1499E+04,.1446E+04,.1401E+04,.1324E+04,.1311E+04,.1404E+04,
     A .1491E+04,.1463E+04,.1390E+04,.1381E+04,.1426E+04,.1442E+04,
     A .1387E+04,.1289E+04,.1209E+04,.1186E+04,.1237E+04,.1322E+04,
     A .1445E+04,.1572E+04,.1570E+04,.1517E+04,.1424E+04,.1396E+04,
     A .1427E+04,.1473E+04,.1525E+04,.1521E+04,.1460E+04,.1316E+04,
     A .1200E+04,.1196E+04,.1275E+04,.1281E+04,.1165E+04,.1104E+04,
     A .1153E+04,.1206E+04,.1270E+04,.1287E+04,.1298E+04,.1346E+04,
     A .1406E+04,.1434E+04,.1398E+04,.1305E+04,.1239E+04,.1226E+04,
     A .1188E+04,.1091E+04,.1028E+04,.1052E+04,.1118E+04,.1181E+04,
     A .1272E+04,.1406E+04,.1519E+04,.1587E+04,.1674E+04,.1681E+04,
     A .1612E+04,.1454E+04,.1345E+04,.1322E+04,.1348E+04,.1404E+04,
     A .1419E+04,.1429E+04,.1429E+04,.1415E+04,.1421E+04,.1465E+04,
     A .1535E+04,.1605E+04,.1682E+04,.1703E+04,.1662E+04,.1603E+04,
     A .1580E+04,.1553E+04,.1555E+04,.1448E+04,.1341E+04,.1230E+04,
     A .1169E+04,.1034E+04,.9538E+03,.1031E+04,.1105E+04,.1110E+04,
     A .1130E+04,.1155E+04,.1124E+04,.1058E+04,.1078E+04,.1214E+04,
     A .1362E+04,.1446E+04,.1504E+04,.1565E+04,.1605E+04,.1670E+04,
     A .1806E+04,.1873E+04,.1843E+04,.1748E+04,.1689E+04,.1697E+04,
     A .1720E+04,.1734E+04,.1691E+04,.1532E+04,.1361E+04,.1107E+04,
     A .8051E+03,.7811E+03,.7817E+03,.8148E+03,.8157E+03,.8107E+03,
     A .8819E+03,.1046E+04,.1222E+04,.1325E+04,.1406E+04,.1524E+04,
     A .1641E+04,.1689E+04,.1561E+04,.1282E+04,.1099E+04,.1147E+04,
     A .1336E+04,.1503E+04,.1555E+04,.1491E+04,.1380E+04,.1347E+04,
     A .1331E+04,.1287E+04,.1199E+04,.1064E+04,.9209E+03,.7933E+03,
     A .6282E+03,.6147E+03,.6570E+03,.6910E+03,.7247E+03,.7646E+03,
     A .8230E+03,.9333E+03,.1105E+04,.1213E+04,.1314E+04,.1435E+04,
     A .1519E+04,.1564E+04,.1544E+04,.1449E+04,.1347E+04,.1159E+04,
     A .1080E+04,.1129E+04,.1098E+04,.1017E+04,.8800E+03,.8049E+03,
     A .8752E+03,.9807E+03,.9965E+03,.8883E+03,.7312E+03,.6077E+03,
     A .5726E+03,.6666E+03,.8569E+03,.1066E+04,.1237E+04,.1378E+04,
     A .1398E+04,.1399E+04,.1431E+04,.1511E+04,.1594E+04,.1699E+04,
     A .1749E+04,.1650E+04,.1476E+04,.1277E+04,.1078E+04,.9043E+03,
     A .8110E+03,.8613E+03,.1003E+04,.1216E+04,.1297E+04,.1288E+04,
     A .1296E+04,.1316E+04,.1306E+04,.1236E+04,.1089E+04,.1044E+04,
     A .1180E+04,.1314E+04,.1360E+04,.1329E+04,.1198E+04,.1045E+04,
     A .1011E+04,.1113E+04,.1260E+04,.1355E+04,.1345E+04,.1247E+04,
     A .1169E+04,.1177E+04,.1273E+04,.1394E+04,.1500E+04,.1563E+04,
     A .1562E+04,.1562E+04,.1550E+04,.1528E+04,.1516E+04,.1535E+04,
     A .1582E+04,.1609E+04,.1620E+04,.1579E+04,.1550E+04,.1566E+04,
     A .1638E+04,.1772E+04,.1817E+04,.1796E+04,.1780E+04,.1758E+04,
     A .1758E+04,.1761E+04,.1809E+04,.1915E+04,.1966E+04,.1917E+04,
     A .1701E+04,.1463E+04,.1362E+04,.1340E+04,.1348E+04,.1361E+04,
     A .1358E+04,.1312E+04,.1237E+04,.1211E+04,.1266E+04,.1384E+04,
     A .1500E+04,.1516E+04,.1458E+04,.1365E+04,.1319E+04,.1266E+04,
     A .1177E+04,.1124E+04,.1219E+04,.1273E+04,.1274E+04,.1184E+04,
     A .1059E+04,.9639E+03,.9365E+03,.9434E+03,.1024E+04,.1195E+04,
     A .1389E+04,.1427E+04,.1396E+04,.1439E+04,.1524E+04,.1611E+04,
     A .1673E+04,.1700E+04,.1670E+04,.1640E+04,.1626E+04,.1548E+04,
     A .1426E+04,.1329E+04,.1325E+04,.1428E+04,.1537E+04,.1613E+04,
     A .1684E+04,.1748E+04,.1760E+04,.1743E+04,.1734E+04,.1692E+04,
     A .1584E+04,.1400E+04,.1300E+04,.1286E+04,.1239E+04,.1165E+04,
     A .1004E+04,.9337E+03,.1063E+04,.1255E+04,.1390E+04,.1485E+04,
     A .1500E+04,.1384E+04,.1159E+04,.9092E+03,.7824E+03,.8276E+03,
     A .9830E+03,.1179E+04,.1353E+04,.1407E+04,.1409E+04,.1275E+04,
     A .1055E+04,.9449E+03,.8593E+03,.7391E+03,.7232E+03,.8005E+03,
     A .8223E+03,.8026E+03,.8433E+03,.8475E+03,.8288E+03,.8856E+03,
     A .1046E+04,.1030E+04,.9680E+03,.8231E+03,.7172E+03,.6958E+03,
     A .7088E+03,.6896E+03,.6861E+03,.7303E+03,.8134E+03,.8977E+03,
     A .9574E+03,.9766E+03,.9505E+03,.8389E+03,.6674E+03,.6716E+03,
     A .7375E+03,.7487E+03,.7074E+03,.7244E+03,.9065E+03,.1076E+04,
     A .1263E+04,.1317E+04,.1319E+04,.1369E+04,.1387E+04,.1375E+04,
     A .1337E+04,.1329E+04,.1320E+04,.1365E+04,.1462E+04,.1529E+04,
     A .1532E+04,.1458E+04,.1290E+04,.1132E+04,.1099E+04,.1169E+04,
     A .1325E+04,.1490E+04,.1474E+04,.1440E+04,.1391E+04,.1292E+04,
     A .1255E+04,.1261E+04,.1146E+04,.1066E+04,.1109E+04,.1197E+04,
     A .1223E+04,.1166E+04,.1045E+04,.8941E+03,.7599E+03,.7428E+03,
     A .8013E+03,.8925E+03,.1008E+04,.1168E+04,.1370E+04,.1469E+04,
     A .1464E+04,.1450E+04,.1389E+04,.1301E+04,.1240E+04,.1257E+04,
     A .1286E+04,.1279E+04,.1248E+04,.1259E+04,.1243E+04,.1205E+04,
     A .1140E+04,.1056E+04,.9471E+03,.8610E+03,.8292E+03,.8148E+03,
     A .8160E+03,.8821E+03,.9795E+03,.1103E+04,.1314E+04,.1410E+04,
     A .1428E+04,.1422E+04,.1385E+04,.1323E+04,.1172E+04,.1010E+04,
     A .8595E+03,.9283E+03,.1050E+04,.1121E+04,.1177E+04,.1168E+04,
     A .1067E+04,.9670E+03,.8663E+03,.9637E+03,.1221E+04,.1289E+04,
     A .1338E+04,.1257E+04,.1088E+04,.9100E+03,.8546E+03,.9448E+03,
     A .1042E+04,.1040E+04,.9700E+03,.9269E+03,.9954E+03,.1114E+04,
     A .1247E+04,.1368E+04,.1469E+04,.1518E+04,.1497E+04,.1457E+04,
     A .1457E+04,.1503E+04,.1499E+04,.1408E+04,.1302E+04,.1327E+04,
     A .1491E+04,.1642E+04,.1617E+04,.1570E+04,.1556E+04,.1503E+04,
     A .1395E+04,.1357E+04,.1433E+04,.1546E+04,.1660E+04,.1621E+04,
     A .1449E+04,.1315E+04,.1301E+04,.1375E+04,.1410E+04,.1309E+04,
     A .1161E+04,.1090E+04,.1160E+04,.1312E+04,.1461E+04,.1562E+04,
     A .1633E+04,.1668E+04,.1657E+04,.1592E+04,.1548E+04,.1594E+04,
     A .1646E+04,.1714E+04,.1780E+04,.1750E+04,.1666E+04,.1602E+04,
     A .1605E+04,.1663E+04,.1648E+04,.1584E+04,.1522E+04,.1536E+04,
     A .1585E+04,.1566E+04,.1507E+04,.1490E+04,.1450E+04,.1358E+04,
     A .1230E+04,.1144E+04,.1175E+04,.1236E+04,.1244E+04,.1243E+04,
     A .1323E+04,.1408E+04,.1444E+04,.1434E+04,.1346E+04,.1267E+04,
     A .1256E+04,.1221E+04,.1102E+04,.9681E+03,.8773E+03,.8266E+03,
     A .7843E+03,.7218E+03,.6724E+03,.6265E+03,.5740E+03,.5256E+03,
     A .4636E+03,.4048E+03,.3421E+03,.3047E+03,.3290E+03,.3833E+03,
     A .4618E+03,.5208E+03,.6259E+03,.6878E+03,.7329E+03,.7867E+03,
     A .8414E+03,.9267E+03,.1041E+04,.1121E+04,.1165E+04,.1192E+04,
     A .1212E+04,.1224E+04,.1229E+04,.1216E+04,.1150E+04,.1046E+04,
     A .9997E+03,.1084E+04,.1307E+04,.1510E+04,.1534E+04,.1537E+04,
     A .1488E+04,.1418E+04,.1366E+04,.1372E+04,.1443E+04,.1539E+04,
     A .1626E+04,.1689E+04,.1720E+04,.1665E+04,.1428E+04,.1363E+04,
     A .1421E+04,.1614E+04,.1690E+04,.1700E+04,.1646E+04,.1476E+04,
     A .1423E+04,.1437E+04,.1466E+04,.1500E+04,.1545E+04,.1611E+04,
     A .1675E+04,.1684E+04,.1529E+04,.1228E+04,.1073E+04,.9796E+03,
     A .9243E+03,.9428E+03,.9888E+03,.1004E+04,.1003E+04,.9634E+03,
     A .8960E+03,.7893E+03,.6855E+03,.5908E+03,.4950E+03,.4087E+03,
     A .3496E+03,.3154E+03,.3322E+03,.4038E+03,.4827E+03,.5917E+03,
     A .6922E+03,.7639E+03,.8546E+03,.9448E+03,.9788E+03,.1006E+04,
     A .1057E+04,.1164E+04,.1281E+04,.1353E+04,.1382E+04,.1406E+04,
     A .1520E+04,.1578E+04,.1640E+04,.1723E+04,.1791E+04,.1866E+04,
     A .1924E+04,.1942E+04,.1894E+04,.1824E+04,.1784E+04,.1770E+04,
     A .1747E+04,.1706E+04,.1702E+04,.1717E+04,.1780E+04,.1814E+04,
     A .1725E+04,.1648E+04,.1708E+04,.1820E+04,.1912E+04,.1978E+04,
     A .1925E+04,.1852E+04,.1806E+04,.1829E+04,.1893E+04,.1923E+04,
     A .1948E+04,.1982E+04,.1991E+04,.1980E+04,.1941E+04,.1896E+04,
     A .1887E+04,.1891E+04,.1887E+04,.1826E+04,.1747E+04,.1715E+04,
     A .1712E+04,.1811E+04,.1972E+04,.2095E+04,.2160E+04,.2114E+04/      
      REAL lambda
      INTEGER ierr
      INTEGER i, j, n
      CHARACTER*40 FIL

*_______________________________________________________________________

******* SUSIM irradiance 
*_______________________________________________________________________
* VanHoosier, M. E., J.-D. F. Bartoe, G. E. Brueckner, and
* D. K. Prinz, Absolute solar spectral irradiance 120 nm -
* 400 nm (Results from the Solar Ultraviolet Spectral Irradiance
* Monitor - SUSIM- Experiment on board Spacelab 2), 
* Astro. Lett. and Communications, 1988, vol. 27, pp. 163-168.
*     SUSIM SL2 high resolution (0.15nm) Solar Irridance data.
*     Irradiance values are given in milliwatts/m^2/nanomenters
*     and are listed at 0.05nm intervals.  The wavelength given is
*     the center wavelength of the 0.15nm triangular bandpass.
*     Normalized to 1 astronomical unit.
*  DATA for wavelengths > 350 nm are unreliable
* (Van Hoosier, personal communication, 1994).
*_______________________________________________________________________




* compute wavelengths, convert from mW to W

      n = nsusim
      DO 13, i = 1, n
         lambda_hi(i)=120.5 + FLOAT(i-1)*.05
         irrad_hi(i) = irrad_hi(i)  /  1000.
   13 CONTINUE
*_______________________________________________________________________

      CALL addpnt(lambda_hi,irrad_hi,10000,n,
     >            lambda_hi(1)*(1.-deltax),0.)
      CALL addpnt(lambda_hi,irrad_hi,10000,n,                 0.,0.)
      CALL addpnt(lambda_hi,irrad_hi,10000,n,
     >            lambda_hi(n)*(1.+deltax),0.)
      CALL addpnt(lambda_hi,irrad_hi,10000,n,              1.e38,0.)
      CALL inter2(nw,wl,f,n,lambda_hi,irrad_hi,ierr)
      IF (ierr .NE. 0) THEN
         WRITE(*,*) ierr, fil
         STOP
      ENDIF

*_______________________________________________________________________

      RETURN
      END
      SUBROUTINE read2(nw,wl,f)

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Read extra-terrestrial flux data.  Re-grid data to match specified       =*
*=  working wavelength grid.                                                 =*
*-----------------------------------------------------------------------------*
*=  PARAMETERS:                                                              =*
*=  NW     - INTEGER, number of specified intervals + 1 in working        (I)=*
*=           wavelength grid                                                 =*
*=  WL     - REAL, vector of lower limits of wavelength intervals in      (I)=*
*=           working wavelength grid                                         =*
*=  F      - REAL, spectral irradiance at the top of the atmosphere at    (O)=*
*=           each specified wavelength                                       =*
*-----------------------------------------------------------------------------*
*=  EDIT HISTORY:                                                            =*
*=  Original                                                                 =*
*-----------------------------------------------------------------------------*
*= This program is free software;  you can redistribute it and/or modify     =*
*= it under the terms of the GNU General Public License as published by the  =*
*= Free Software Foundation;  either version 2 of the license, or (at your   =*
*= option) any later version.                                                =*
*= The TUV package is distributed in the hope that it will be useful, but    =*
*= WITHOUT ANY WARRANTY;  without even the implied warranty of MERCHANTIBI-  =*
*= LITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public     =*
*= License for more details.                                                 =*
*= To obtain a copy of the GNU General Public License, write to:             =*
*= Free Software Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.   =*
*-----------------------------------------------------------------------------*
*= To contact the authors, please mail to:                                   =*
*= Sasha Madronich, NCAR/ACD, P.O.Box 3000, Boulder, CO, 80307-3000, USA  or =*
*= send email to:  sasha@ucar.edu                                            =*
*-----------------------------------------------------------------------------*
*= Copyright (C) 1994,95,96  University Corporation for Atmospheric Research =*
*-----------------------------------------------------------------------------*

      IMPLICIT NONE
      INCLUDE 'tuv.params'

* input: (wavelength grid)
      INTEGER nw
      REAL wl(kw)
      REAL yg(kw)

*
      INTEGER iw

* output: (extra terrestrial solar flux)
      REAL f(kw)

* local:

      REAL x1(1000), y1(1000) 
      REAL x2(1000), y2(1000)
      REAL x3(1000), y3(1000)
      INTEGER i, n
      REAL DUM
      INTEGER IDUM

*_______________________________________________________________________

*********WMO 85 irradiance

      OPEN(UNIT=kin,FILE='DATAE1/SUN/wmo85.flx',STATUS='old')
      DO 11, i = 1, 3
         READ(kin,*)
   11 CONTINUE
      n = 158
      DO 12, i = 1, n
         READ(kin,*) idum, x1(i),x2(i),y1(i), dum, dum, dum
         x3(i) = 0.5 * (x1(i) + x2(i))

C average value needs to be calculated only if inter2 is
C used to interpolate onto wavelength grid (see below)
C        y1(i) =  y1(i) / (x2(i) - x1(i)) 

   12 CONTINUE
      CLOSE (kin)

      x1(n+1) = x2(n)

C inter2: INPUT : average value in each bin 
C         OUTPUT: average value in each bin
C inter3: INPUT : total area in each bin
C         OUTPUT: total area in each bin

      CALL inter3(nw,wl,yg, n+1,x1,y1,0)
C      CALL inter2(nw,wl,yg,n,x3,y1,ierr)

      DO 10,  iw = 1, nw-1
* from quanta s-1 cm-2 bin-1 to  watts m-2 nm-1
* 1.e4 * ([hc =] 6.62E-34 * 2.998E8)/(wc*1e-9) 
         
C the scaling by bin width needs to be done only if
C inter3 is used for interpolation

         yg(iw) = yg(iw) / (wl(iw+1)-wl(iw))
         f(iw) = yg(iw) * 1.e4 * (6.62E-34 * 2.998E8) / 
     $        ( 0.5 * (wl(iw+1)+wl(iw)) * 1.e-9)

   10 CONTINUE
      
*_______________________________________________________________________

      RETURN
      END
      FUNCTION refrac(w)

      IMPLICIT NONE

* input vacuum wavelength, nm

      REAL w

* output refractive index for standard air
* (dry air at 15 deg. C, 101.325 kPa, 0.03% CO2)

      REAL refrac

* internal

      REAL sig,  dum

* from CRC Handbook, originally from Edlen, B., Metrologia, 2, 71, 1966.
* valid from 200 nm to 2000 nm
* beyond this range, use constant value

      sig = 1.E3/w

      IF (w .LT. 200.) sig = 1.E3/200.
      IF (w .GT. 2000.) sig = 1.E3/2000.

      dum = 8342.13 + 2406030./(130. - sig*sig) + 
     $     15997./(38.9 - sig*sig)
      refrac = 1. + 1.E-8 * dum


      RETURN
      END

      Subroutine read_atlas(x1,y1,kdata,n) ! data atlas extra-terrestrial flux data
      parameter (natlas=5160)
      real x1(kdata),y1(kdata),x1_atlas(natlas),y1_atlas(natlas)
      
       data x1_atlas/150.01,150.06,150.11,150.16,150.21,
     A 150.26,150.31,150.36,150.41,150.46,150.51,150.56,150.61,150.66,
     A 150.71,150.76,150.81,150.86,150.91,150.96,151.01,151.06,151.11,
     A 151.16,151.21,151.26,151.31,151.36,151.41,151.46,151.51,151.56,
     A 151.61,151.66,151.71,151.76,151.81,151.86,151.91,151.96,152.01,
     A 152.06,152.11,152.16,152.21,152.26,152.31,152.36,152.41,152.46,
     A 152.51,152.56,152.61,152.66,152.71,152.76,152.81,152.86,152.91,
     A 152.96,153.01,153.06,153.11,153.16,153.21,153.26,153.31,153.36,
     A 153.41,153.46,153.51,153.56,153.61,153.66,153.71,153.76,153.81,
     A 153.86,153.91,153.96,154.01,154.06,154.11,154.16,154.21,154.26,
     A 154.31,154.36,154.41,154.46,154.51,154.56,154.61,154.66,154.71,
     A 154.76,154.81,154.86,154.91,154.96,155.01,155.06,155.11,155.16,
     A 155.21,155.26,155.31,155.36,155.41,155.46,155.51,155.56,155.61,
     A 155.66,155.71,155.76,155.81,155.86,155.91,155.96,156.01,156.06,
     A 156.11,156.16,156.21,156.26,156.31,156.36,156.41,156.46,156.51,
     A 156.56,156.61,156.66,156.71,156.76,156.81,156.86,156.91,156.96,
     A 157.01,157.06,157.11,157.16,157.21,157.26,157.31,157.36,157.41,
     A 157.46,157.51,157.56,157.61,157.66,157.71,157.76,157.81,157.86,
     A 157.91,157.96,158.01,158.06,158.11,158.16,158.21,158.26,158.31,
     A 158.36,158.41,158.46,158.51,158.56,158.61,158.66,158.71,158.76,
     A 158.81,158.86,158.91,158.96,159.01,159.06,159.11,159.16,159.21,
     A 159.26,159.31,159.36,159.41,159.46,159.51,159.56,159.61,159.66,
     A 159.71,159.76,159.81,159.86,159.91,159.96,160.01,160.06,160.11,
     A 160.16,160.21,160.26,160.31,160.36,160.41,160.46,160.51,160.56,
     A 160.61,160.66,160.71,160.76,160.81,160.86,160.91,160.96,161.01,
     A 161.06,161.11,161.16,161.21,161.26,161.31,161.36,161.41,161.46,
     A 161.51,161.56,161.61,161.66,161.71,161.76,161.81,161.86,161.91,
     A 161.96,162.01,162.06,162.11,162.16,162.21,162.26,162.31,162.36,
     A 162.41,162.46,162.51,162.56,162.61,162.66,162.71,162.76,162.81,
     A 162.86,162.91,162.96,163.01,163.06,163.11,163.16,163.21,163.26,
     A 163.31,163.36,163.41,163.46,163.51,163.56,163.61,163.66,163.71,
     A 163.76,163.81,163.86,163.91,163.96,164.01,164.06,164.11,164.16,
     A 164.21,164.26,164.31,164.36,164.41,164.46,164.51,164.56,164.61,
     A 164.66,164.71,164.76,164.81,164.86,164.91,164.96,165.01,165.06,
     A 165.11,165.16,165.21,165.26,165.31,165.36,165.41,165.46,165.51,
     A 165.56,165.61,165.66,165.71,165.76,165.81,165.86,165.91,165.96,
     A 166.01,166.06,166.11,166.16,166.21,166.26,166.31,166.36,166.41,
     A 166.46,166.51,166.56,166.61,166.66,166.71,166.76,166.81,166.86,
     A 166.91,166.96,167.01,167.06,167.11,167.16,167.21,167.26,167.31,
     A 167.36,167.41,167.46,167.51,167.56,167.61,167.66,167.71,167.76,
     A 167.81,167.86,167.91,167.96,168.01,168.06,168.11,168.16,168.21,
     A 168.26,168.31,168.36,168.41,168.46,168.51,168.56,168.61,168.66,
     A 168.71,168.76,168.81,168.86,168.91,168.96,169.01,169.06,169.11,
     A 169.16,169.21,169.26,169.31,169.36,169.41,169.46,169.51,169.56,
     A 169.61,169.66,169.71,169.76,169.81,169.86,169.91,169.96,170.01,
     A 170.06,170.11,170.16,170.21,170.26,170.31,170.36,170.41,170.46,
     A 170.51,170.56,170.61,170.66,170.71,170.76,170.81,170.86,170.91,
     A 170.96,171.01,171.06,171.11,171.16,171.21,171.26,171.31,171.36,
     A 171.41,171.46,171.51,171.56,171.61,171.66,171.71,171.76,171.81,
     A 171.86,171.91,171.96,172.01,172.06,172.11,172.16,172.21,172.26,
     A 172.31,172.36,172.41,172.46,172.51,172.56,172.61,172.66,172.71,
     A 172.76,172.81,172.86,172.91,172.96,173.01,173.06,173.11,173.16,
     A 173.21,173.26,173.31,173.36,173.41,173.46,173.51,173.56,173.61,
     A 173.66,173.71,173.76,173.81,173.86,173.91,173.96,174.01,174.06,
     A 174.11,174.16,174.21,174.26,174.31,174.36,174.41,174.46,174.51,
     A 174.56,174.61,174.66,174.71,174.76,174.81,174.86,174.91,174.96,
     A 175.01,175.06,175.11,175.16,175.21,175.26,175.31,175.36,175.41,
     A 175.46,175.51,175.56,175.61,175.66,175.71,175.76,175.81,175.86,
     A 175.91,175.96,176.01,176.06,176.11,176.16,176.21,176.26,176.31,
     A 176.36,176.41,176.46,176.51,176.56,176.61,176.66,176.71,176.76,
     A 176.81,176.86,176.91,176.96,177.01,177.06,177.11,177.16,177.21,
     A 177.26,177.31,177.36,177.41,177.46,177.51,177.56,177.61,177.66,
     A 177.71,177.76,177.81,177.86,177.91,177.96,178.01,178.06,178.11,
     A 178.16,178.21,178.26,178.31,178.36,178.41,178.46,178.51,178.56,
     A 178.61,178.66,178.71,178.76,178.81,178.86,178.91,178.96,179.01,
     A 179.06,179.11,179.16,179.21,179.26,179.31,179.36,179.41,179.46,
     A 179.51,179.56,179.61,179.66,179.71,179.76,179.81,179.86,179.91,
     A 179.96,180.01,180.06,180.11,180.16,180.21,180.26,180.31,180.36,
     A 180.41,180.46,180.51,180.56,180.61,180.66,180.71,180.76,180.81,
     A 180.86,180.91,180.96,181.01,181.06,181.11,181.16,181.21,181.26,
     A 181.31,181.36,181.41,181.46,181.51,181.56,181.61,181.66,181.71,
     A 181.76,181.81,181.86,181.91,181.96,182.01,182.06,182.11,182.16,
     A 182.21,182.26,182.31,182.36,182.41,182.46,182.51,182.56,182.61,
     A 182.66,182.71,182.76,182.81,182.86,182.91,182.96,183.01,183.06,
     A 183.11,183.16,183.21,183.26,183.31,183.36,183.41,183.46,183.51,
     A 183.56,183.61,183.66,183.71,183.76,183.81,183.86,183.91,183.96,
     A 184.01,184.06,184.11,184.16,184.21,184.26,184.31,184.36,184.41,
     A 184.46,184.51,184.56,184.61,184.66,184.71,184.76,184.81,184.86,
     A 184.91,184.96,185.01,185.06,185.11,185.16,185.21,185.26,185.31,
     A 185.36,185.41,185.46,185.51,185.56,185.61,185.66,185.71,185.76,
     A 185.81,185.86,185.91,185.96,186.01,186.06,186.11,186.16,186.21,
     A 186.26,186.31,186.36,186.41,186.46,186.51,186.56,186.61,186.66,
     A 186.71,186.76,186.81,186.86,186.91,186.96,187.01,187.06,187.11,
     A 187.16,187.21,187.26,187.31,187.36,187.41,187.46,187.51,187.56,
     A 187.61,187.66,187.71,187.76,187.81,187.86,187.91,187.96,188.01,
     A 188.06,188.11,188.16,188.21,188.26,188.31,188.36,188.41,188.46,
     A 188.51,188.56,188.61,188.66,188.71,188.76,188.81,188.86,188.91,
     A 188.96,189.01,189.06,189.11,189.16,189.21,189.26,189.31,189.36,
     A 189.41,189.46,189.51,189.56,189.61,189.66,189.71,189.76,189.81,
     A 189.86,189.91,189.96,190.01,190.06,190.11,190.16,190.21,190.26,
     A 190.31,190.36,190.41,190.46,190.51,190.56,190.61,190.66,190.71,
     A 190.76,190.81,190.86,190.91,190.96,191.01,191.06,191.11,191.16,
     A 191.21,191.26,191.31,191.36,191.41,191.46,191.51,191.56,191.61,
     A 191.66,191.71,191.76,191.81,191.86,191.91,191.96,192.01,192.06,
     A 192.11,192.16,192.21,192.26,192.31,192.36,192.41,192.46,192.51,
     A 192.56,192.61,192.66,192.71,192.76,192.81,192.86,192.91,192.96,
     A 193.01,193.06,193.11,193.16,193.21,193.26,193.31,193.36,193.41,
     A 193.46,193.51,193.56,193.61,193.66,193.71,193.76,193.81,193.86,
     A 193.91,193.96,194.01,194.06,194.11,194.16,194.21,194.26,194.31,
     A 194.36,194.41,194.46,194.51,194.56,194.61,194.66,194.71,194.76,
     A 194.81,194.86,194.91,194.96,195.01,195.06,195.11,195.16,195.21,
     A 195.26,195.31,195.36,195.41,195.46,195.51,195.56,195.61,195.66,
     A 195.71,195.76,195.81,195.86,195.91,195.96,196.01,196.06,196.11,
     A 196.16,196.21,196.26,196.31,196.36,196.41,196.46,196.51,196.56,
     A 196.61,196.66,196.71,196.76,196.81,196.86,196.91,196.96,197.01,
     A 197.06,197.11,197.16,197.21,197.26,197.31,197.36,197.41,197.46,
     A 197.51,197.56,197.61,197.66,197.71,197.76,197.81,197.86,197.91,
     A 197.96,198.01,198.06,198.11,198.16,198.21,198.26,198.31,198.36,
     A 198.41,198.46,198.51,198.56,198.61,198.66,198.71,198.76,198.81,
     A 198.86,198.91,198.96,199.01,199.06,199.11,199.16,199.21,199.26,
     A 199.31,199.36,199.41,199.46,199.51,199.56,199.61,199.66,199.71,
     A 199.76,199.81,199.86,199.91,199.96,200.01,200.06,200.11,200.16,
     A 200.21,200.26,200.31,200.36,200.41,200.46,200.51,200.56,200.61,
     A 200.66,200.71,200.76,200.81,200.86,200.91,200.96,201.01,201.06,
     A 201.11,201.16,201.21,201.26,201.31,201.36,201.41,201.46,201.51,
     A 201.56,201.61,201.66,201.71,201.76,201.81,201.86,201.91,201.96,
     A 202.01,202.06,202.11,202.16,202.21,202.26,202.31,202.36,202.41,
     A 202.46,202.51,202.56,202.61,202.66,202.71,202.76,202.81,202.86,
     A 202.91,202.96,203.01,203.06,203.11,203.16,203.21,203.26,203.31,
     A 203.36,203.41,203.46,203.51,203.56,203.61,203.66,203.71,203.76,
     A 203.81,203.86,203.91,203.96,204.01,204.06,204.11,204.16,204.21,
     A 204.26,204.31,204.36,204.41,204.46,204.51,204.56,204.61,204.66,
     A 204.71,204.76,204.81,204.86,204.91,204.96,205.01,205.06,205.11,
     A 205.16,205.21,205.26,205.31,205.36,205.41,205.46,205.51,205.56,
     A 205.61,205.66,205.71,205.76,205.81,205.86,205.91,205.96,206.01,
     A 206.06,206.11,206.16,206.21,206.26,206.31,206.36,206.41,206.46,
     A 206.51,206.56,206.61,206.66,206.71,206.76,206.81,206.86,206.91,
     A 206.96,207.01,207.06,207.11,207.16,207.21,207.26,207.31,207.36,
     A 207.41,207.46,207.51,207.56,207.61,207.66,207.71,207.76,207.81,
     A 207.86,207.91,207.96,208.01,208.06,208.11,208.16,208.21,208.26,
     A 208.31,208.36,208.41,208.46,208.51,208.56,208.61,208.66,208.71,
     A 208.76,208.81,208.86,208.91,208.96,209.01,209.06,209.11,209.16,
     A 209.21,209.26,209.31,209.36,209.41,209.46,209.51,209.56,209.61,
     A 209.66,209.71,209.76,209.81,209.86,209.91,209.96,210.01,210.06,
     A 210.11,210.16,210.21,210.26,210.31,210.36,210.41,210.46,210.51,
     A 210.56,210.61,210.66,210.71,210.76,210.81,210.86,210.91,210.96,
     A 211.01,211.06,211.11,211.16,211.21,211.26,211.31,211.36,211.41,
     A 211.46,211.51,211.56,211.61,211.66,211.71,211.76,211.81,211.86,
     A 211.91,211.96,212.01,212.06,212.11,212.16,212.21,212.26,212.31,
     A 212.36,212.41,212.46,212.51,212.56,212.61,212.66,212.71,212.76,
     A 212.81,212.86,212.91,212.96,213.01,213.06,213.11,213.16,213.21,
     A 213.26,213.31,213.36,213.41,213.46,213.51,213.56,213.61,213.66,
     A 213.71,213.76,213.81,213.86,213.91,213.96,214.01,214.06,214.11,
     A 214.16,214.21,214.26,214.31,214.36,214.41,214.46,214.51,214.56,
     A 214.61,214.66,214.71,214.76,214.81,214.86,214.91,214.96,215.01,
     A 215.06,215.11,215.16,215.21,215.26,215.31,215.36,215.41,215.46,
     A 215.51,215.56,215.61,215.66,215.71,215.76,215.81,215.86,215.91,
     A 215.96,216.01,216.06,216.11,216.16,216.21,216.26,216.31,216.36,
     A 216.41,216.46,216.51,216.56,216.61,216.66,216.71,216.76,216.81,
     A 216.86,216.91,216.96,217.01,217.06,217.11,217.16,217.21,217.26,
     A 217.31,217.36,217.41,217.46,217.51,217.56,217.61,217.66,217.71,
     A 217.76,217.81,217.86,217.91,217.96,218.01,218.06,218.11,218.16,
     A 218.21,218.26,218.31,218.36,218.41,218.46,218.51,218.56,218.61,
     A 218.66,218.71,218.76,218.81,218.86,218.91,218.96,219.01,219.06,
     A 219.11,219.16,219.21,219.26,219.31,219.36,219.41,219.46,219.51,
     A 219.56,219.61,219.66,219.71,219.76,219.81,219.86,219.91,219.96,
     A 220.01,220.06,220.11,220.16,220.21,220.26,220.31,220.36,220.41,
     A 220.46,220.51,220.56,220.61,220.66,220.71,220.76,220.81,220.86,
     A 220.91,220.96,221.01,221.06,221.11,221.16,221.21,221.26,221.31,
     A 221.36,221.41,221.46,221.51,221.56,221.61,221.66,221.71,221.76,
     A 221.81,221.86,221.91,221.96,222.01,222.06,222.11,222.16,222.21,
     A 222.26,222.31,222.36,222.41,222.46,222.51,222.56,222.61,222.66,
     A 222.71,222.76,222.81,222.86,222.91,222.96,223.01,223.06,223.11,
     A 223.16,223.21,223.26,223.31,223.36,223.41,223.46,223.51,223.56,
     A 223.61,223.66,223.71,223.76,223.81,223.86,223.91,223.96,224.01,
     A 224.06,224.11,224.16,224.21,224.26,224.31,224.36,224.41,224.46,
     A 224.51,224.56,224.61,224.66,224.71,224.76,224.81,224.86,224.91,
     A 224.96,225.01,225.06,225.11,225.16,225.21,225.26,225.31,225.36,
     A 225.41,225.46,225.51,225.56,225.61,225.66,225.71,225.76,225.81,
     A 225.86,225.91,225.96,226.01,226.06,226.11,226.16,226.21,226.26,
     A 226.31,226.36,226.41,226.46,226.51,226.56,226.61,226.66,226.71,
     A 226.76,226.81,226.86,226.91,226.96,227.01,227.06,227.11,227.16,
     A 227.21,227.26,227.31,227.36,227.41,227.46,227.51,227.56,227.61,
     A 227.66,227.71,227.76,227.81,227.86,227.91,227.96,228.01,228.06,
     A 228.11,228.16,228.21,228.26,228.31,228.36,228.41,228.46,228.51,
     A 228.56,228.61,228.66,228.71,228.76,228.81,228.86,228.91,228.96,
     A 229.01,229.06,229.11,229.16,229.21,229.26,229.31,229.36,229.41,
     A 229.46,229.51,229.56,229.61,229.66,229.71,229.76,229.81,229.86,
     A 229.91,229.96,230.01,230.06,230.11,230.16,230.21,230.26,230.31,
     A 230.36,230.41,230.46,230.51,230.56,230.61,230.66,230.71,230.76,
     A 230.81,230.86,230.91,230.96,231.01,231.06,231.11,231.16,231.21,
     A 231.26,231.31,231.36,231.41,231.46,231.51,231.56,231.61,231.66,
     A 231.71,231.76,231.81,231.86,231.91,231.96,232.01,232.06,232.11,
     A 232.16,232.21,232.26,232.31,232.36,232.41,232.46,232.51,232.56,
     A 232.61,232.66,232.71,232.76,232.81,232.86,232.91,232.96,233.01,
     A 233.06,233.11,233.16,233.21,233.26,233.31,233.36,233.41,233.46,
     A 233.51,233.56,233.61,233.66,233.71,233.76,233.81,233.86,233.91,
     A 233.96,234.01,234.06,234.11,234.16,234.21,234.26,234.31,234.36,
     A 234.41,234.46,234.51,234.56,234.61,234.66,234.71,234.76,234.81,
     A 234.86,234.91,234.96,235.01,235.06,235.11,235.16,235.21,235.26,
     A 235.31,235.36,235.41,235.46,235.51,235.56,235.61,235.66,235.71,
     A 235.76,235.81,235.86,235.91,235.96,236.01,236.06,236.11,236.16,
     A 236.21,236.26,236.31,236.36,236.41,236.46,236.51,236.56,236.61,
     A 236.66,236.71,236.76,236.81,236.86,236.91,236.96,237.01,237.06,
     A 237.11,237.16,237.21,237.26,237.31,237.36,237.41,237.46,237.51,
     A 237.56,237.61,237.66,237.71,237.76,237.81,237.86,237.91,237.96,
     A 238.01,238.06,238.11,238.16,238.21,238.26,238.31,238.36,238.41,
     A 238.46,238.51,238.56,238.61,238.66,238.71,238.76,238.81,238.86,
     A 238.91,238.96,239.01,239.06,239.11,239.16,239.21,239.26,239.31,
     A 239.36,239.41,239.46,239.51,239.56,239.61,239.66,239.71,239.76,
     A 239.81,239.86,239.91,239.96,240.01,240.06,240.11,240.16,240.21,
     A 240.26,240.31,240.36,240.41,240.46,240.51,240.56,240.61,240.66,
     A 240.71,240.76,240.81,240.86,240.91,240.96,241.01,241.06,241.11,
     A 241.16,241.21,241.26,241.31,241.36,241.41,241.46,241.51,241.56,
     A 241.61,241.66,241.71,241.76,241.81,241.86,241.91,241.96,242.01,
     A 242.06,242.11,242.16,242.21,242.26,242.31,242.36,242.41,242.46,
     A 242.51,242.56,242.61,242.66,242.71,242.76,242.81,242.86,242.91,
     A 242.96,243.01,243.06,243.11,243.16,243.21,243.26,243.31,243.36,
     A 243.41,243.46,243.51,243.56,243.61,243.66,243.71,243.76,243.81,
     A 243.86,243.91,243.96,244.01,244.06,244.11,244.16,244.21,244.26,
     A 244.31,244.36,244.41,244.46,244.51,244.56,244.61,244.66,244.71,
     A 244.76,244.81,244.86,244.91,244.96,245.01,245.06,245.11,245.16,
     A 245.21,245.26,245.31,245.36,245.41,245.46,245.51,245.56,245.61,
     A 245.66,245.71,245.76,245.81,245.86,245.91,245.96,246.01,246.06,
     A 246.11,246.16,246.21,246.26,246.31,246.36,246.41,246.46,246.51,
     A 246.56,246.61,246.66,246.71,246.76,246.81,246.86,246.91,246.96,
     A 247.01,247.06,247.11,247.16,247.21,247.26,247.31,247.36,247.41,
     A 247.46,247.51,247.56,247.61,247.66,247.71,247.76,247.81,247.86,
     A 247.91,247.96,248.01,248.06,248.11,248.16,248.21,248.26,248.31,
     A 248.36,248.41,248.46,248.51,248.56,248.61,248.66,248.71,248.76,
     A 248.81,248.86,248.91,248.96,249.01,249.06,249.11,249.16,249.21,
     A 249.26,249.31,249.36,249.41,249.46,249.51,249.56,249.61,249.66,
     A 249.71,249.76,249.81,249.86,249.91,249.96,250.01,250.06,250.11,
     A 250.16,250.21,250.26,250.31,250.36,250.41,250.46,250.51,250.56,
     A 250.61,250.66,250.71,250.76,250.81,250.86,250.91,250.96,251.01,
     A 251.06,251.11,251.16,251.21,251.26,251.31,251.36,251.41,251.46,
     A 251.51,251.56,251.61,251.66,251.71,251.76,251.81,251.86,251.91,
     A 251.96,252.01,252.06,252.11,252.16,252.21,252.26,252.31,252.36,
     A 252.41,252.46,252.51,252.56,252.61,252.66,252.71,252.76,252.81,
     A 252.86,252.91,252.96,253.01,253.06,253.11,253.16,253.21,253.26,
     A 253.31,253.36,253.41,253.46,253.51,253.56,253.61,253.66,253.71,
     A 253.76,253.81,253.86,253.91,253.96,254.01,254.06,254.11,254.16,
     A 254.21,254.26,254.31,254.36,254.41,254.46,254.51,254.56,254.61,
     A 254.66,254.71,254.76,254.81,254.86,254.91,254.96,255.01,255.06,
     A 255.11,255.16,255.21,255.26,255.31,255.36,255.41,255.46,255.51,
     A 255.56,255.61,255.66,255.71,255.76,255.81,255.86,255.91,255.96,
     A 256.01,256.06,256.11,256.16,256.21,256.26,256.31,256.36,256.41,
     A 256.46,256.51,256.56,256.61,256.66,256.71,256.76,256.81,256.86,
     A 256.91,256.96,257.01,257.06,257.11,257.16,257.21,257.26,257.31,
     A 257.36,257.41,257.46,257.51,257.56,257.61,257.66,257.71,257.76,
     A 257.81,257.86,257.91,257.96,258.01,258.06,258.11,258.16,258.21,
     A 258.26,258.31,258.36,258.41,258.46,258.51,258.56,258.61,258.66,
     A 258.71,258.76,258.81,258.86,258.91,258.96,259.01,259.06,259.11,
     A 259.16,259.21,259.26,259.31,259.36,259.41,259.46,259.51,259.56,
     A 259.61,259.66,259.71,259.76,259.81,259.86,259.91,259.96,260.01,
     A 260.06,260.11,260.16,260.21,260.26,260.31,260.36,260.41,260.46,
     A 260.51,260.56,260.61,260.66,260.71,260.76,260.81,260.86,260.91,
     A 260.96,261.01,261.06,261.11,261.16,261.21,261.26,261.31,261.36,
     A 261.41,261.46,261.51,261.56,261.61,261.66,261.71,261.76,261.81,
     A 261.86,261.91,261.96,262.01,262.06,262.11,262.16,262.21,262.26,
     A 262.31,262.36,262.41,262.46,262.51,262.56,262.61,262.66,262.71,
     A 262.76,262.81,262.86,262.91,262.96,263.01,263.06,263.11,263.16,
     A 263.21,263.26,263.31,263.36,263.41,263.46,263.51,263.56,263.61,
     A 263.66,263.71,263.76,263.81,263.86,263.91,263.96,264.01,264.06,
     A 264.11,264.16,264.21,264.26,264.31,264.36,264.41,264.46,264.51,
     A 264.56,264.61,264.66,264.71,264.76,264.81,264.86,264.91,264.96,
     A 265.01,265.06,265.11,265.16,265.21,265.26,265.31,265.36,265.41,
     A 265.46,265.51,265.56,265.61,265.66,265.71,265.76,265.81,265.86,
     A 265.91,265.96,266.01,266.06,266.11,266.16,266.21,266.26,266.31,
     A 266.36,266.41,266.46,266.51,266.56,266.61,266.66,266.71,266.76,
     A 266.81,266.86,266.91,266.96,267.01,267.06,267.11,267.16,267.21,
     A 267.26,267.31,267.36,267.41,267.46,267.51,267.56,267.61,267.66,
     A 267.71,267.76,267.81,267.86,267.91,267.96,268.01,268.06,268.11,
     A 268.16,268.21,268.26,268.31,268.36,268.41,268.46,268.51,268.56,
     A 268.61,268.66,268.71,268.76,268.81,268.86,268.91,268.96,269.01,
     A 269.06,269.11,269.16,269.21,269.26,269.31,269.36,269.41,269.46,
     A 269.51,269.56,269.61,269.66,269.71,269.76,269.81,269.86,269.91,
     A 269.96,270.01,270.06,270.11,270.16,270.21,270.26,270.31,270.36,
     A 270.41,270.46,270.51,270.56,270.61,270.66,270.71,270.76,270.81,
     A 270.86,270.91,270.96,271.01,271.06,271.11,271.16,271.21,271.26,
     A 271.31,271.36,271.41,271.46,271.51,271.56,271.61,271.66,271.71,
     A 271.76,271.81,271.86,271.91,271.96,272.01,272.06,272.11,272.16,
     A 272.21,272.26,272.31,272.36,272.41,272.46,272.51,272.56,272.61,
     A 272.66,272.71,272.76,272.81,272.86,272.91,272.96,273.01,273.06,
     A 273.11,273.16,273.21,273.26,273.31,273.36,273.41,273.46,273.51,
     A 273.56,273.61,273.66,273.71,273.76,273.81,273.86,273.91,273.96,
     A 274.01,274.06,274.11,274.16,274.21,274.26,274.31,274.36,274.41,
     A 274.46,274.51,274.56,274.61,274.66,274.71,274.76,274.81,274.86,
     A 274.91,274.96,275.01,275.06,275.11,275.16,275.21,275.26,275.31,
     A 275.36,275.41,275.46,275.51,275.56,275.61,275.66,275.71,275.76,
     A 275.81,275.86,275.91,275.96,276.01,276.06,276.11,276.16,276.21,
     A 276.26,276.31,276.36,276.41,276.46,276.51,276.56,276.61,276.66,
     A 276.71,276.76,276.81,276.86,276.91,276.96,277.01,277.06,277.11,
     A 277.16,277.21,277.26,277.31,277.36,277.41,277.46,277.51,277.56,
     A 277.61,277.66,277.71,277.76,277.81,277.86,277.91,277.96,278.01,
     A 278.06,278.11,278.16,278.21,278.26,278.31,278.36,278.41,278.46,
     A 278.51,278.56,278.61,278.66,278.71,278.76,278.81,278.86,278.91,
     A 278.96,279.01,279.06,279.11,279.16,279.21,279.26,279.31,279.36,
     A 279.41,279.46,279.51,279.56,279.61,279.66,279.71,279.76,279.81,
     A 279.86,279.91,279.96,280.01,280.06,280.11,280.16,280.21,280.26,
     A 280.31,280.36,280.41,280.46,280.51,280.56,280.61,280.66,280.71,
     A 280.76,280.81,280.86,280.91,280.96,281.01,281.06,281.11,281.16,
     A 281.21,281.26,281.31,281.36,281.41,281.46,281.51,281.56,281.61,
     A 281.66,281.71,281.76,281.81,281.86,281.91,281.96,282.01,282.06,
     A 282.11,282.16,282.21,282.26,282.31,282.36,282.41,282.46,282.51,
     A 282.56,282.61,282.66,282.71,282.76,282.81,282.86,282.91,282.96,
     A 283.01,283.06,283.11,283.16,283.21,283.26,283.31,283.36,283.41,
     A 283.46,283.51,283.56,283.61,283.66,283.71,283.76,283.81,283.86,
     A 283.91,283.96,284.01,284.06,284.11,284.16,284.21,284.26,284.31,
     A 284.36,284.41,284.46,284.51,284.56,284.61,284.66,284.71,284.76,
     A 284.81,284.86,284.91,284.96,285.01,285.06,285.11,285.16,285.21,
     A 285.26,285.31,285.36,285.41,285.46,285.51,285.56,285.61,285.66,
     A 285.71,285.76,285.81,285.86,285.91,285.96,286.01,286.06,286.11,
     A 286.16,286.21,286.26,286.31,286.36,286.41,286.46,286.51,286.56,
     A 286.61,286.66,286.71,286.76,286.81,286.86,286.91,286.96,287.01,
     A 287.06,287.11,287.16,287.21,287.26,287.31,287.36,287.41,287.46,
     A 287.51,287.56,287.61,287.66,287.71,287.76,287.81,287.86,287.91,
     A 287.96,288.01,288.06,288.11,288.16,288.21,288.26,288.31,288.36,
     A 288.41,288.46,288.51,288.56,288.61,288.66,288.71,288.76,288.81,
     A 288.86,288.91,288.96,289.01,289.06,289.11,289.16,289.21,289.26,
     A 289.31,289.36,289.41,289.46,289.51,289.56,289.61,289.66,289.71,
     A 289.76,289.81,289.86,289.91,289.96,290.01,290.06,290.11,290.16,
     A 290.21,290.26,290.31,290.36,290.41,290.46,290.51,290.56,290.61,
     A 290.66,290.71,290.76,290.81,290.86,290.91,290.96,291.01,291.06,
     A 291.11,291.16,291.21,291.26,291.31,291.36,291.41,291.46,291.51,
     A 291.56,291.61,291.66,291.71,291.76,291.81,291.86,291.91,291.96,
     A 292.01,292.06,292.11,292.16,292.21,292.26,292.31,292.36,292.41,
     A 292.46,292.51,292.56,292.61,292.66,292.71,292.76,292.81,292.86,
     A 292.91,292.96,293.01,293.06,293.11,293.16,293.21,293.26,293.31,
     A 293.36,293.41,293.46,293.51,293.56,293.61,293.66,293.71,293.76,
     A 293.81,293.86,293.91,293.96,294.01,294.06,294.11,294.16,294.21,
     A 294.26,294.31,294.36,294.41,294.46,294.51,294.56,294.61,294.66,
     A 294.71,294.76,294.81,294.86,294.91,294.96,295.01,295.06,295.11,
     A 295.16,295.21,295.26,295.31,295.36,295.41,295.46,295.51,295.56,
     A 295.61,295.66,295.71,295.76,295.81,295.86,295.91,295.96,296.01,
     A 296.06,296.11,296.16,296.21,296.26,296.31,296.36,296.41,296.46,
     A 296.51,296.56,296.61,296.66,296.71,296.76,296.81,296.86,296.91,
     A 296.96,297.01,297.06,297.11,297.16,297.21,297.26,297.31,297.36,
     A 297.41,297.46,297.51,297.56,297.61,297.66,297.71,297.76,297.81,
     A 297.86,297.91,297.96,298.01,298.06,298.11,298.16,298.21,298.26,
     A 298.31,298.36,298.41,298.46,298.51,298.56,298.61,298.66,298.71,
     A 298.76,298.81,298.86,298.91,298.96,299.01,299.06,299.11,299.16,
     A 299.21,299.26,299.31,299.36,299.41,299.46,299.51,299.56,299.61,
     A 299.66,299.71,299.76,299.81,299.86,299.91,299.96,300.01,300.06,
     A 300.11,300.16,300.21,300.26,300.31,300.36,300.41,300.46,300.51,
     A 300.56,300.61,300.66,300.71,300.76,300.81,300.86,300.91,300.96,
     A 301.01,301.06,301.11,301.16,301.21,301.26,301.31,301.36,301.41,
     A 301.46,301.51,301.56,301.61,301.66,301.71,301.76,301.81,301.86,
     A 301.91,301.96,302.01,302.06,302.11,302.16,302.21,302.26,302.31,
     A 302.36,302.41,302.46,302.51,302.56,302.61,302.66,302.71,302.76,
     A 302.81,302.86,302.91,302.96,303.01,303.06,303.11,303.16,303.21,
     A 303.26,303.31,303.36,303.41,303.46,303.51,303.56,303.61,303.66,
     A 303.71,303.76,303.81,303.86,303.91,303.96,304.01,304.06,304.11,
     A 304.16,304.21,304.26,304.31,304.36,304.41,304.46,304.51,304.56,
     A 304.61,304.66,304.71,304.76,304.81,304.86,304.91,304.96,305.01,
     A 305.06,305.11,305.16,305.21,305.26,305.31,305.36,305.41,305.46,
     A 305.51,305.56,305.61,305.66,305.71,305.76,305.81,305.86,305.91,
     A 305.96,306.01,306.06,306.11,306.16,306.21,306.26,306.31,306.36,
     A 306.41,306.46,306.51,306.56,306.61,306.66,306.71,306.76,306.81,
     A 306.86,306.91,306.96,307.01,307.06,307.11,307.16,307.21,307.26,
     A 307.31,307.36,307.41,307.46,307.51,307.56,307.61,307.66,307.71,
     A 307.76,307.81,307.86,307.91,307.96,308.01,308.06,308.11,308.16,
     A 308.21,308.26,308.31,308.36,308.41,308.46,308.51,308.56,308.61,
     A 308.66,308.71,308.76,308.81,308.86,308.91,308.96,309.01,309.06,
     A 309.11,309.16,309.21,309.26,309.31,309.36,309.41,309.46,309.51,
     A 309.56,309.61,309.66,309.71,309.76,309.81,309.86,309.91,309.96,
     A 310.01,310.06,310.11,310.16,310.21,310.26,310.31,310.36,310.41,
     A 310.46,310.51,310.56,310.61,310.66,310.71,310.76,310.81,310.86,
     A 310.91,310.96,311.01,311.06,311.11,311.16,311.21,311.26,311.31,
     A 311.36,311.41,311.46,311.51,311.56,311.61,311.66,311.71,311.76,
     A 311.81,311.86,311.91,311.96,312.01,312.06,312.11,312.16,312.21,
     A 312.26,312.31,312.36,312.41,312.46,312.51,312.56,312.61,312.66,
     A 312.71,312.76,312.81,312.86,312.91,312.96,313.01,313.06,313.11,
     A 313.16,313.21,313.26,313.31,313.36,313.41,313.46,313.51,313.56,
     A 313.61,313.66,313.71,313.76,313.81,313.86,313.91,313.96,314.01,
     A 314.06,314.11,314.16,314.21,314.26,314.31,314.36,314.41,314.46,
     A 314.51,314.56,314.61,314.66,314.71,314.76,314.81,314.86,314.91,
     A 314.96,315.01,315.06,315.11,315.16,315.21,315.26,315.31,315.36,
     A 315.41,315.46,315.51,315.56,315.61,315.66,315.71,315.76,315.81,
     A 315.86,315.91,315.96,316.01,316.06,316.11,316.16,316.21,316.26,
     A 316.31,316.36,316.41,316.46,316.51,316.56,316.61,316.66,316.71,
     A 316.76,316.81,316.86,316.91,316.96,317.01,317.06,317.11,317.16,
     A 317.21,317.26,317.31,317.36,317.41,317.46,317.51,317.56,317.61,
     A 317.66,317.71,317.76,317.81,317.86,317.91,317.96,318.01,318.06,
     A 318.11,318.16,318.21,318.26,318.31,318.36,318.41,318.46,318.51,
     A 318.56,318.61,318.66,318.71,318.76,318.81,318.86,318.91,318.96,
     A 319.01,319.06,319.11,319.16,319.21,319.26,319.31,319.36,319.41,
     A 319.46,319.51,319.56,319.61,319.66,319.71,319.76,319.81,319.86,
     A 319.91,319.96,320.01,320.06,320.11,320.16,320.21,320.26,320.31,
     A 320.36,320.41,320.46,320.51,320.56,320.61,320.66,320.71,320.76,
     A 320.81,320.86,320.91,320.96,321.01,321.06,321.11,321.16,321.21,
     A 321.26,321.31,321.36,321.41,321.46,321.51,321.56,321.61,321.66,
     A 321.71,321.76,321.81,321.86,321.91,321.96,322.01,322.06,322.11,
     A 322.16,322.21,322.26,322.31,322.36,322.41,322.46,322.51,322.56,
     A 322.61,322.66,322.71,322.76,322.81,322.86,322.91,322.96,323.01,
     A 323.06,323.11,323.16,323.21,323.26,323.31,323.36,323.41,323.46,
     A 323.51,323.56,323.61,323.66,323.71,323.76,323.81,323.86,323.91,
     A 323.96,324.01,324.06,324.11,324.16,324.21,324.26,324.31,324.36,
     A 324.41,324.46,324.51,324.56,324.61,324.66,324.71,324.76,324.81,
     A 324.86,324.91,324.96,325.01,325.06,325.11,325.16,325.21,325.26,
     A 325.31,325.36,325.41,325.46,325.51,325.56,325.61,325.66,325.71,
     A 325.76,325.81,325.86,325.91,325.96,326.01,326.06,326.11,326.16,
     A 326.21,326.26,326.31,326.36,326.41,326.46,326.51,326.56,326.61,
     A 326.66,326.71,326.76,326.81,326.86,326.91,326.96,327.01,327.06,
     A 327.11,327.16,327.21,327.26,327.31,327.36,327.41,327.46,327.51,
     A 327.56,327.61,327.66,327.71,327.76,327.81,327.86,327.91,327.96,
     A 328.01,328.06,328.11,328.16,328.21,328.26,328.31,328.36,328.41,
     A 328.46,328.51,328.56,328.61,328.66,328.71,328.76,328.81,328.86,
     A 328.91,328.96,329.01,329.06,329.11,329.16,329.21,329.26,329.31,
     A 329.36,329.41,329.46,329.51,329.56,329.61,329.66,329.71,329.76,
     A 329.81,329.86,329.91,329.96,330.01,330.06,330.11,330.16,330.21,
     A 330.26,330.31,330.36,330.41,330.46,330.51,330.56,330.61,330.66,
     A 330.71,330.76,330.81,330.86,330.91,330.96,331.01,331.06,331.11,
     A 331.16,331.21,331.26,331.31,331.36,331.41,331.46,331.51,331.56,
     A 331.61,331.66,331.71,331.76,331.81,331.86,331.91,331.96,332.01,
     A 332.06,332.11,332.16,332.21,332.26,332.31,332.36,332.41,332.46,
     A 332.51,332.56,332.61,332.66,332.71,332.76,332.81,332.86,332.91,
     A 332.96,333.01,333.06,333.11,333.16,333.21,333.26,333.31,333.36,
     A 333.41,333.46,333.51,333.56,333.61,333.66,333.71,333.76,333.81,
     A 333.86,333.91,333.96,334.01,334.06,334.11,334.16,334.21,334.26,
     A 334.31,334.36,334.41,334.46,334.51,334.56,334.61,334.66,334.71,
     A 334.76,334.81,334.86,334.91,334.96,335.01,335.06,335.11,335.16,
     A 335.21,335.26,335.31,335.36,335.41,335.46,335.51,335.56,335.61,
     A 335.66,335.71,335.76,335.81,335.86,335.91,335.96,336.01,336.06,
     A 336.11,336.16,336.21,336.26,336.31,336.36,336.41,336.46,336.51,
     A 336.56,336.61,336.66,336.71,336.76,336.81,336.86,336.91,336.96,
     A 337.01,337.06,337.11,337.16,337.21,337.26,337.31,337.36,337.41,
     A 337.46,337.51,337.56,337.61,337.66,337.71,337.76,337.81,337.86,
     A 337.91,337.96,338.01,338.06,338.11,338.16,338.21,338.26,338.31,
     A 338.36,338.41,338.46,338.51,338.56,338.61,338.66,338.71,338.76,
     A 338.81,338.86,338.91,338.96,339.01,339.06,339.11,339.16,339.21,
     A 339.26,339.31,339.36,339.41,339.46,339.51,339.56,339.61,339.66,
     A 339.71,339.76,339.81,339.86,339.91,339.96,340.01,340.06,340.11,
     A 340.16,340.21,340.26,340.31,340.36,340.41,340.46,340.51,340.56,
     A 340.61,340.66,340.71,340.76,340.81,340.86,340.91,340.96,341.01,
     A 341.06,341.11,341.16,341.21,341.26,341.31,341.36,341.41,341.46,
     A 341.51,341.56,341.61,341.66,341.71,341.76,341.81,341.86,341.91,
     A 341.96,342.01,342.06,342.11,342.16,342.21,342.26,342.31,342.36,
     A 342.41,342.46,342.51,342.56,342.61,342.66,342.71,342.76,342.81,
     A 342.86,342.91,342.96,343.01,343.06,343.11,343.16,343.21,343.26,
     A 343.31,343.36,343.41,343.46,343.51,343.56,343.61,343.66,343.71,
     A 343.76,343.81,343.86,343.91,343.96,344.01,344.06,344.11,344.16,
     A 344.21,344.26,344.31,344.36,344.41,344.46,344.51,344.56,344.61,
     A 344.66,344.71,344.76,344.81,344.86,344.91,344.96,345.01,345.06,
     A 345.11,345.16,345.21,345.26,345.31,345.36,345.41,345.46,345.51,
     A 345.56,345.61,345.66,345.71,345.76,345.81,345.86,345.91,345.96,
     A 346.01,346.06,346.11,346.16,346.21,346.26,346.31,346.36,346.41,
     A 346.46,346.51,346.56,346.61,346.66,346.71,346.76,346.81,346.86,
     A 346.91,346.96,347.01,347.06,347.11,347.16,347.21,347.26,347.31,
     A 347.36,347.41,347.46,347.51,347.56,347.61,347.66,347.71,347.76,
     A 347.81,347.86,347.91,347.96,348.01,348.06,348.11,348.16,348.21,
     A 348.26,348.31,348.36,348.41,348.46,348.51,348.56,348.61,348.66,
     A 348.71,348.76,348.81,348.86,348.91,348.96,349.01,349.06,349.11,
     A 349.16,349.21,349.26,349.31,349.36,349.41,349.46,349.51,349.56,
     A 349.61,349.66,349.71,349.76,349.81,349.86,349.91,349.96,350.01,
     A 350.06,350.11,350.16,350.21,350.26,350.31,350.36,350.41,350.46,
     A 350.51,350.56,350.61,350.66,350.71,350.76,350.81,350.86,350.91,
     A 350.96,351.01,351.06,351.11,351.16,351.21,351.26,351.31,351.36,
     A 351.41,351.46,351.51,351.56,351.61,351.66,351.71,351.76,351.81,
     A 351.86,351.91,351.96,352.01,352.06,352.11,352.16,352.21,352.26,
     A 352.31,352.36,352.41,352.46,352.51,352.56,352.61,352.66,352.71,
     A 352.76,352.81,352.86,352.91,352.96,353.01,353.06,353.11,353.16,
     A 353.21,353.26,353.31,353.36,353.41,353.46,353.51,353.56,353.61,
     A 353.66,353.71,353.76,353.81,353.86,353.91,353.96,354.01,354.06,
     A 354.11,354.16,354.21,354.26,354.31,354.36,354.41,354.46,354.51,
     A 354.56,354.61,354.66,354.71,354.76,354.81,354.86,354.91,354.96,
     A 355.01,355.06,355.11,355.16,355.21,355.26,355.31,355.36,355.41,
     A 355.46,355.51,355.56,355.61,355.66,355.71,355.76,355.81,355.86,
     A 355.91,355.96,356.01,356.06,356.11,356.16,356.21,356.26,356.31,
     A 356.36,356.41,356.46,356.51,356.56,356.61,356.66,356.71,356.76,
     A 356.81,356.86,356.91,356.96,357.01,357.06,357.11,357.16,357.21,
     A 357.26,357.31,357.36,357.41,357.46,357.51,357.56,357.61,357.66,
     A 357.71,357.76,357.81,357.86,357.91,357.96,358.01,358.06,358.11,
     A 358.16,358.21,358.26,358.31,358.36,358.41,358.46,358.51,358.56,
     A 358.61,358.66,358.71,358.76,358.81,358.86,358.91,358.96,359.01,
     A 359.06,359.11,359.16,359.21,359.26,359.31,359.36,359.41,359.46,
     A 359.51,359.56,359.61,359.66,359.71,359.76,359.81,359.86,359.91,
     A 359.96,360.01,360.06,360.11,360.16,360.21,360.26,360.31,360.36,
     A 360.41,360.46,360.51,360.56,360.61,360.66,360.71,360.76,360.81,
     A 360.86,360.91,360.96,361.01,361.06,361.11,361.16,361.21,361.26,
     A 361.31,361.36,361.41,361.46,361.51,361.56,361.61,361.66,361.71,
     A 361.76,361.81,361.86,361.91,361.96,362.01,362.06,362.11,362.16,
     A 362.21,362.26,362.31,362.36,362.41,362.46,362.51,362.56,362.61,
     A 362.66,362.71,362.76,362.81,362.86,362.91,362.96,363.01,363.06,
     A 363.11,363.16,363.21,363.26,363.31,363.36,363.41,363.46,363.51,
     A 363.56,363.61,363.66,363.71,363.76,363.81,363.86,363.91,363.96,
     A 364.01,364.06,364.11,364.16,364.21,364.26,364.31,364.36,364.41,
     A 364.46,364.51,364.56,364.61,364.66,364.71,364.76,364.81,364.86,
     A 364.91,364.96,365.01,365.06,365.11,365.16,365.21,365.26,365.31,
     A 365.36,365.41,365.46,365.51,365.56,365.61,365.66,365.71,365.76,
     A 365.81,365.86,365.91,365.96,366.01,366.06,366.11,366.16,366.21,
     A 366.26,366.31,366.36,366.41,366.46,366.51,366.56,366.61,366.66,
     A 366.71,366.76,366.81,366.86,366.91,366.96,367.01,367.06,367.11,
     A 367.16,367.21,367.26,367.31,367.36,367.41,367.46,367.51,367.56,
     A 367.61,367.66,367.71,367.76,367.81,367.86,367.91,367.96,368.01,
     A 368.06,368.11,368.16,368.21,368.26,368.31,368.36,368.41,368.46,
     A 368.51,368.56,368.61,368.66,368.71,368.76,368.81,368.86,368.91,
     A 368.96,369.01,369.06,369.11,369.16,369.21,369.26,369.31,369.36,
     A 369.41,369.46,369.51,369.56,369.61,369.66,369.71,369.76,369.81,
     A 369.86,369.91,369.96,370.01,370.06,370.11,370.16,370.21,370.26,
     A 370.31,370.36,370.41,370.46,370.51,370.56,370.61,370.66,370.71,
     A 370.76,370.81,370.86,370.91,370.96,371.01,371.06,371.11,371.16,
     A 371.21,371.26,371.31,371.36,371.41,371.46,371.51,371.56,371.61,
     A 371.66,371.71,371.76,371.81,371.86,371.91,371.96,372.01,372.06,
     A 372.11,372.16,372.21,372.26,372.31,372.36,372.41,372.46,372.51,
     A 372.56,372.61,372.66,372.71,372.76,372.81,372.86,372.91,372.96,
     A 373.01,373.06,373.11,373.16,373.21,373.26,373.31,373.36,373.41,
     A 373.46,373.51,373.56,373.61,373.66,373.71,373.76,373.81,373.86,
     A 373.91,373.96,374.01,374.06,374.11,374.16,374.21,374.26,374.31,
     A 374.36,374.41,374.46,374.51,374.56,374.61,374.66,374.71,374.76,
     A 374.81,374.86,374.91,374.96,375.01,375.06,375.11,375.16,375.21,
     A 375.26,375.31,375.36,375.41,375.46,375.51,375.56,375.61,375.66,
     A 375.71,375.76,375.81,375.86,375.91,375.96,376.01,376.06,376.11,
     A 376.16,376.21,376.26,376.31,376.36,376.41,376.46,376.51,376.56,
     A 376.61,376.66,376.71,376.76,376.81,376.86,376.91,376.96,377.01,
     A 377.06,377.11,377.16,377.21,377.26,377.31,377.36,377.41,377.46,
     A 377.51,377.56,377.61,377.66,377.71,377.76,377.81,377.86,377.91,
     A 377.96,378.01,378.06,378.11,378.16,378.21,378.26,378.31,378.36,
     A 378.41,378.46,378.51,378.56,378.61,378.66,378.71,378.76,378.81,
     A 378.86,378.91,378.96,379.01,379.06,379.11,379.16,379.21,379.26,
     A 379.31,379.36,379.41,379.46,379.51,379.56,379.61,379.66,379.71,
     A 379.76,379.81,379.86,379.91,379.96,380.01,380.06,380.11,380.16,
     A 380.21,380.26,380.31,380.36,380.41,380.46,380.51,380.56,380.61,
     A 380.66,380.71,380.76,380.81,380.86,380.91,380.96,381.01,381.06,
     A 381.11,381.16,381.21,381.26,381.31,381.36,381.41,381.46,381.51,
     A 381.56,381.61,381.66,381.71,381.76,381.81,381.86,381.91,381.96,
     A 382.01,382.06,382.11,382.16,382.21,382.26,382.31,382.36,382.41,
     A 382.46,382.51,382.56,382.61,382.66,382.71,382.76,382.81,382.86,
     A 382.91,382.96,383.01,383.06,383.11,383.16,383.21,383.26,383.31,
     A 383.36,383.41,383.46,383.51,383.56,383.61,383.66,383.71,383.76,
     A 383.81,383.86,383.91,383.96,384.01,384.06,384.11,384.16,384.21,
     A 384.26,384.31,384.36,384.41,384.46,384.51,384.56,384.61,384.66,
     A 384.71,384.76,384.81,384.86,384.91,384.96,385.01,385.06,385.11,
     A 385.16,385.21,385.26,385.31,385.36,385.41,385.46,385.51,385.56,
     A 385.61,385.66,385.71,385.76,385.81,385.86,385.91,385.96,386.01,
     A 386.06,386.11,386.16,386.21,386.26,386.31,386.36,386.41,386.46,
     A 386.51,386.56,386.61,386.66,386.71,386.76,386.81,386.86,386.91,
     A 386.96,387.01,387.06,387.11,387.16,387.21,387.26,387.31,387.36,
     A 387.41,387.46,387.51,387.56,387.61,387.66,387.71,387.76,387.81,
     A 387.86,387.91,387.96,388.01,388.06,388.11,388.16,388.21,388.26,
     A 388.31,388.36,388.41,388.46,388.51,388.56,388.61,388.66,388.71,
     A 388.76,388.81,388.86,388.91,388.96,389.01,389.06,389.11,389.16,
     A 389.21,389.26,389.31,389.36,389.41,389.46,389.51,389.56,389.61,
     A 389.66,389.71,389.76,389.81,389.86,389.91,389.96,390.01,390.06,
     A 390.11,390.16,390.21,390.26,390.31,390.36,390.41,390.46,390.51,
     A 390.56,390.61,390.66,390.71,390.76,390.81,390.86,390.91,390.96,
     A 391.01,391.06,391.11,391.16,391.21,391.26,391.31,391.36,391.41,
     A 391.46,391.51,391.56,391.61,391.66,391.71,391.76,391.81,391.86,
     A 391.91,391.96,392.01,392.06,392.11,392.16,392.21,392.26,392.31,
     A 392.36,392.41,392.46,392.51,392.56,392.61,392.66,392.71,392.76,
     A 392.81,392.86,392.91,392.96,393.01,393.06,393.11,393.16,393.21,
     A 393.26,393.31,393.36,393.41,393.46,393.51,393.56,393.61,393.66,
     A 393.71,393.76,393.81,393.86,393.91,393.96,394.01,394.06,394.11,
     A 394.16,394.21,394.26,394.31,394.36,394.41,394.46,394.51,394.56,
     A 394.61,394.66,394.71,394.76,394.81,394.86,394.91,394.96,395.01,
     A 395.06,395.11,395.16,395.21,395.26,395.31,395.36,395.41,395.46,
     A 395.51,395.56,395.61,395.66,395.71,395.76,395.81,395.86,395.91,
     A 395.96,396.01,396.06,396.11,396.16,396.21,396.26,396.31,396.36,
     A 396.41,396.46,396.51,396.56,396.61,396.66,396.71,396.76,396.81,
     A 396.86,396.91,396.96,397.01,397.06,397.11,397.16,397.21,397.26,
     A 397.31,397.36,397.41,397.46,397.51,397.56,397.61,397.66,397.71,
     A 397.76,397.81,397.86,397.91,397.96,398.01,398.06,398.11,398.16,
     A 398.21,398.26,398.31,398.36,398.41,398.46,398.51,398.56,398.61,
     A 398.66,398.71,398.76,398.81,398.86,398.91,398.96,399.01,399.06,
     A 399.11,399.16,399.21,399.26,399.31,399.36,399.41,399.46,399.51,
     A 399.56,399.61,399.66,399.71,399.76,399.81,399.86,399.91,399.96,
     A 400.01,400.06,400.11,400.16,400.21,400.26,400.31,400.36,400.41,
     A 400.46,400.51,400.56,400.61,400.66,400.71,400.76,400.81,400.86,
     A 400.91,400.96,401.01,401.06,401.11,401.16,401.21,401.26,401.31,
     A 401.36,401.41,401.46,401.51,401.56,401.61,401.66,401.71,401.76,
     A 401.81,401.86,401.91,401.96,402.01,402.06,402.11,402.16,402.21,
     A 402.26,402.31,402.36,402.41,402.46,402.51,402.56,402.61,402.66,
     A 402.71,402.76,402.81,402.86,402.91,402.96,403.01,403.06,403.11,
     A 403.16,403.21,403.26,403.31,403.36,403.41,403.46,403.51,403.56,
     A 403.61,403.66,403.71,403.76,403.81,403.86,403.91,403.96,404.01,
     A 404.06,404.11,404.16,404.21,404.26,404.31,404.36,404.41,404.46,
     A 404.51,404.56,404.61,404.66,404.71,404.76,404.81,404.86,404.91,
     A 404.96,405.01,405.06,405.11,405.16,405.21,405.26,405.31,405.36,
     A 405.41,405.46,405.51,405.56,405.61,405.66,405.71,405.76,405.81,
     A 405.86,405.91,405.96,406.01,406.06,406.11,406.16,406.21,406.26,
     A 406.31,406.36,406.41,406.46,406.51,406.56,406.61,406.66,406.71,
     A 406.76,406.81,406.86,406.91,406.96,407.01,407.06,407.11,407.16,
     A 407.21,407.26,407.31,407.36,407.41,407.46,407.51,407.56,407.61,
     A 407.66,407.71,407.76,407.81,407.86,407.91,407.96/
     
       data y1_atlas/.1016E+00,.9724E-01,.9506E-01,.9816E-01,
     A .9452E-01,.9815E-01,.1016E+00,.9245E-01,.9415E-01,.9548E-01,
     A .9849E-01,.9356E-01,.9330E-01,.9741E-01,.9714E-01,.9723E-01,
     A .9314E-01,.9572E-01,.9805E-01,.1020E+00,.1158E+00,.1210E+00,
     A .1245E+00,.1112E+00,.1020E+00,.1028E+00,.9836E-01,.9947E-01,
     A .9813E-01,.1006E+00,.9890E-01,.9947E-01,.1057E+00,.1001E+00,
     A .1006E+00,.1061E+00,.1010E+00,.1045E+00,.1008E+00,.1052E+00,
     A .1058E+00,.1051E+00,.1081E+00,.1086E+00,.1104E+00,.1124E+00,
     A .1157E+00,.1159E+00,.1134E+00,.1155E+00,.1199E+00,.1204E+00,
     A .1564E+00,.2108E+00,.1850E+00,.1407E+00,.1172E+00,.1126E+00,
     A .1212E+00,.1169E+00,.1158E+00,.1189E+00,.1255E+00,.1204E+00,
     A .1180E+00,.1426E+00,.2028E+00,.2214E+00,.1804E+00,.1506E+00,
     A .1265E+00,.1214E+00,.1261E+00,.1323E+00,.1469E+00,.1495E+00,
     A .1310E+00,.1167E+00,.1141E+00,.1217E+00,.1463E+00,.1560E+00,
     A .1615E+00,.1585E+00,.1409E+00,.1288E+00,.1202E+00,.1200E+00,
     A .1315E+00,.1383E+00,.1496E+00,.1445E+00,.1243E+00,.1320E+00,
     A .1507E+00,.3318E+00,.7653E+00,.6647E+00,.3737E+00,.1668E+00,
     A .2461E+00,.4521E+00,.5043E+00,.3600E+00,.1846E+00,.1462E+00,
     A .1291E+00,.1265E+00,.1209E+00,.1275E+00,.1363E+00,.1416E+00,
     A .1447E+00,.1423E+00,.1325E+00,.1364E+00,.1652E+00,.2228E+00,
     A .2578E+00,.2466E+00,.2693E+00,.3636E+00,.3754E+00,.3435E+00,
     A .2486E+00,.1598E+00,.1545E+00,.1878E+00,.1846E+00,.1628E+00,
     A .1408E+00,.1390E+00,.1457E+00,.1624E+00,.1738E+00,.1798E+00,
     A .1724E+00,.1619E+00,.1718E+00,.2029E+00,.2217E+00,.2084E+00,
     A .1831E+00,.1698E+00,.1561E+00,.1593E+00,.1678E+00,.2086E+00,
     A .2436E+00,.2589E+00,.2470E+00,.1890E+00,.1644E+00,.1615E+00,
     A .1740E+00,.1628E+00,.1453E+00,.1477E+00,.1421E+00,.1522E+00,
     A .1794E+00,.2062E+00,.1996E+00,.1602E+00,.1481E+00,.1410E+00,
     A .1548E+00,.1617E+00,.1844E+00,.1965E+00,.1967E+00,.1792E+00,
     A .1725E+00,.1690E+00,.1741E+00,.1989E+00,.2114E+00,.2041E+00,
     A .1794E+00,.1819E+00,.1811E+00,.1778E+00,.1838E+00,.1785E+00,
     A .1836E+00,.1790E+00,.1751E+00,.1904E+00,.1936E+00,.2113E+00,
     A .2030E+00,.1866E+00,.1798E+00,.1665E+00,.1677E+00,.1720E+00,
     A .1786E+00,.1686E+00,.1692E+00,.1787E+00,.1920E+00,.1909E+00,
     A .2046E+00,.2215E+00,.2271E+00,.2173E+00,.2034E+00,.1927E+00,
     A .1840E+00,.1917E+00,.1888E+00,.1942E+00,.1900E+00,.1918E+00,
     A .1965E+00,.2146E+00,.2332E+00,.2386E+00,.2010E+00,.1928E+00,
     A .2151E+00,.2490E+00,.2733E+00,.2421E+00,.2398E+00,.2639E+00,
     A .2975E+00,.2870E+00,.2444E+00,.2448E+00,.2145E+00,.2094E+00,
     A .2147E+00,.2261E+00,.2298E+00,.2226E+00,.2409E+00,.2595E+00,
     A .2469E+00,.2300E+00,.2328E+00,.2353E+00,.2542E+00,.2733E+00,
     A .2763E+00,.2861E+00,.2949E+00,.2865E+00,.2515E+00,.2605E+00,
     A .2840E+00,.3372E+00,.3211E+00,.2706E+00,.2700E+00,.2785E+00,
     A .2656E+00,.2743E+00,.2998E+00,.3088E+00,.2848E+00,.2608E+00,
     A .2524E+00,.2580E+00,.2740E+00,.3140E+00,.3138E+00,.3211E+00,
     A .3346E+00,.3074E+00,.2697E+00,.2666E+00,.2659E+00,.2878E+00,
     A .3302E+00,.3439E+00,.3290E+00,.2851E+00,.2898E+00,.3626E+00,
     A .4689E+00,.5131E+00,.3929E+00,.3056E+00,.2932E+00,.2999E+00,
     A .3380E+00,.3758E+00,.3646E+00,.3294E+00,.3146E+00,.3053E+00,
     A .2960E+00,.2956E+00,.3028E+00,.3061E+00,.2987E+00,.2990E+00,
     A .3174E+00,.3295E+00,.3071E+00,.2906E+00,.2882E+00,.2965E+00,
     A .2922E+00,.2956E+00,.2910E+00,.3110E+00,.3676E+00,.3785E+00,
     A .3582E+00,.3886E+00,.6873E+00,.8989E+00,.9837E+00,.1034E+01,
     A .9174E+00,.6692E+00,.5701E+00,.5033E+00,.4063E+00,.3548E+00,
     A .3696E+00,.3622E+00,.3597E+00,.3769E+00,.4057E+00,.4121E+00,
     A .3792E+00,.3431E+00,.3438E+00,.3648E+00,.3917E+00,.4038E+00,
     A .3706E+00,.3579E+00,.3490E+00,.3440E+00,.3571E+00,.3749E+00,
     A .4162E+00,.5798E+00,.5779E+00,.4519E+00,.3319E+00,.3389E+00,
     A .3765E+00,.4198E+00,.4381E+00,.4563E+00,.4125E+00,.3820E+00,
     A .3790E+00,.4150E+00,.4280E+00,.4043E+00,.3929E+00,.3921E+00,
     A .4190E+00,.4334E+00,.4213E+00,.4244E+00,.4586E+00,.4461E+00,
     A .4017E+00,.4052E+00,.4177E+00,.4331E+00,.4329E+00,.4303E+00,
     A .4375E+00,.5000E+00,.5797E+00,.6068E+00,.5653E+00,.4928E+00,
     A .5083E+00,.5403E+00,.5281E+00,.5059E+00,.5432E+00,.5681E+00,
     A .5967E+00,.5972E+00,.5681E+00,.5609E+00,.5562E+00,.5802E+00,
     A .6201E+00,.6364E+00,.6333E+00,.6115E+00,.6060E+00,.7001E+00,
     A .7351E+00,.7093E+00,.6146E+00,.6432E+00,.6808E+00,.6873E+00,
     A .6568E+00,.6472E+00,.7005E+00,.7783E+00,.8337E+00,.7682E+00,
     A .6870E+00,.6870E+00,.6753E+00,.6840E+00,.6715E+00,.6532E+00,
     A .6662E+00,.6700E+00,.6667E+00,.6594E+00,.6959E+00,.7610E+00,
     A .7808E+00,.7715E+00,.7563E+00,.7345E+00,.7154E+00,.6969E+00,
     A .6460E+00,.6498E+00,.7316E+00,.6971E+00,.6363E+00,.6097E+00,
     A .6426E+00,.6797E+00,.7151E+00,.7590E+00,.8054E+00,.7977E+00,
     A .7620E+00,.7558E+00,.7498E+00,.7381E+00,.7831E+00,.8223E+00,
     A .8448E+00,.7603E+00,.7240E+00,.7028E+00,.7028E+00,.6890E+00,
     A .7340E+00,.8046E+00,.8422E+00,.8225E+00,.8033E+00,.8090E+00,
     A .7976E+00,.7924E+00,.8038E+00,.7899E+00,.7439E+00,.6926E+00,
     A .6528E+00,.6727E+00,.6943E+00,.7305E+00,.7444E+00,.7406E+00,
     A .7603E+00,.7747E+00,.7646E+00,.7748E+00,.7795E+00,.7948E+00,
     A .8122E+00,.8041E+00,.8023E+00,.7901E+00,.7995E+00,.8277E+00,
     A .8502E+00,.8461E+00,.8464E+00,.8569E+00,.9070E+00,.9231E+00,
     A .9465E+00,.9472E+00,.9539E+00,.9574E+00,.9611E+00,.1000E+01,
     A .9809E+00,.9520E+00,.9672E+00,.9701E+00,.8925E+00,.8475E+00,
     A .8738E+00,.8996E+00,.9403E+00,.9413E+00,.9545E+00,.1005E+01,
     A .1068E+01,.1103E+01,.1083E+01,.1074E+01,.1075E+01,.1095E+01,
     A .1118E+01,.1124E+01,.1189E+01,.1169E+01,.1115E+01,.1138E+01,
     A .1158E+01,.1187E+01,.1171E+01,.1170E+01,.1210E+01,.1202E+01,
     A .1187E+01,.1191E+01,.1225E+01,.1279E+01,.1278E+01,.1243E+01,
     A .1141E+01,.1112E+01,.1164E+01,.1166E+01,.1170E+01,.1052E+01,
     A .1033E+01,.1180E+01,.1312E+01,.1461E+01,.1382E+01,.1298E+01,
     A .1240E+01,.1258E+01,.1355E+01,.1415E+01,.1467E+01,.1528E+01,
     A .1577E+01,.1663E+01,.1739E+01,.1686E+01,.1585E+01,.1374E+01,
     A .1248E+01,.1175E+01,.1236E+01,.1308E+01,.1310E+01,.1332E+01,
     A .1363E+01,.1479E+01,.1579E+01,.1600E+01,.1598E+01,.1666E+01,
     A .1657E+01,.1567E+01,.1496E+01,.1461E+01,.1529E+01,.1577E+01,
     A .1593E+01,.1577E+01,.1583E+01,.1603E+01,.1649E+01,.1636E+01,
     A .1590E+01,.1594E+01,.1598E+01,.1607E+01,.1735E+01,.1785E+01,
     A .1748E+01,.1742E+01,.1778E+01,.1786E+01,.1707E+01,.1615E+01,
     A .1468E+01,.1410E+01,.1398E+01,.1439E+01,.1494E+01,.1581E+01,
     A .1616E+01,.1598E+01,.1560E+01,.1602E+01,.1664E+01,.1670E+01,
     A .1650E+01,.1704E+01,.1755E+01,.1772E+01,.1804E+01,.1811E+01,
     A .1804E+01,.1807E+01,.1818E+01,.1805E+01,.1813E+01,.1870E+01,
     A .1929E+01,.1902E+01,.1862E+01,.1788E+01,.1792E+01,.2242E+01,
     A .2891E+01,.2453E+01,.2085E+01,.1836E+01,.1904E+01,.1985E+01,
     A .1915E+01,.1887E+01,.1854E+01,.1899E+01,.1858E+01,.1809E+01,
     A .1741E+01,.1778E+01,.1848E+01,.1897E+01,.2315E+01,.3849E+01,
     A .4406E+01,.3851E+01,.3053E+01,.2220E+01,.2220E+01,.2170E+01,
     A .2090E+01,.1962E+01,.2066E+01,.2224E+01,.2277E+01,.2279E+01,
     A .2407E+01,.2493E+01,.2545E+01,.2539E+01,.2415E+01,.2296E+01,
     A .2210E+01,.2173E+01,.2235E+01,.2128E+01,.2073E+01,.2084E+01,
     A .2277E+01,.2290E+01,.2237E+01,.2181E+01,.2208E+01,.2350E+01,
     A .2415E+01,.2388E+01,.2357E+01,.2373E+01,.2455E+01,.2509E+01,
     A .2579E+01,.2489E+01,.2354E+01,.2318E+01,.2374E+01,.2513E+01,
     A .2448E+01,.2476E+01,.2640E+01,.2550E+01,.2273E+01,.2227E+01,
     A .2085E+01,.1975E+01,.2082E+01,.2216E+01,.2296E+01,.2295E+01,
     A .2251E+01,.2300E+01,.2248E+01,.2072E+01,.2029E+01,.2101E+01,
     A .1984E+01,.1829E+01,.1738E+01,.1999E+01,.2254E+01,.2391E+01,
     A .2250E+01,.1929E+01,.1912E+01,.1980E+01,.2034E+01,.2015E+01,
     A .2115E+01,.2348E+01,.2628E+01,.2851E+01,.2913E+01,.2788E+01,
     A .2791E+01,.2808E+01,.2672E+01,.2472E+01,.2422E+01,.2513E+01,
     A .2438E+01,.2389E+01,.2440E+01,.2545E+01,.2580E+01,.2534E+01,
     A .2584E+01,.2706E+01,.2850E+01,.2675E+01,.2700E+01,.2541E+01,
     A .2550E+01,.2543E+01,.2655E+01,.2840E+01,.3075E+01,.3212E+01,
     A .3213E+01,.3222E+01,.3245E+01,.3268E+01,.3297E+01,.3231E+01,
     A .3341E+01,.3338E+01,.3190E+01,.3144E+01,.3122E+01,.3254E+01,
     A .3356E+01,.3316E+01,.3184E+01,.3148E+01,.3150E+01,.3183E+01,
     A .3109E+01,.3124E+01,.3093E+01,.3066E+01,.3056E+01,.3030E+01,
     A .2994E+01,.3036E+01,.3022E+01,.3016E+01,.3198E+01,.3440E+01,
     A .3485E+01,.3553E+01,.3532E+01,.3510E+01,.3539E+01,.3484E+01,
     A .3475E+01,.3535E+01,.3511E+01,.3501E+01,.3410E+01,.3514E+01,
     A .3610E+01,.3713E+01,.3770E+01,.3719E+01,.3835E+01,.3958E+01,
     A .4212E+01,.4035E+01,.3650E+01,.3588E+01,.3811E+01,.4020E+01,
     A .3936E+01,.3793E+01,.3728E+01,.3833E+01,.3731E+01,.3581E+01,
     A .3417E+01,.3279E+01,.3381E+01,.3609E+01,.3704E+01,.3717E+01,
     A .3469E+01,.3506E+01,.3618E+01,.3841E+01,.3981E+01,.4102E+01,
     A .4160E+01,.4122E+01,.4007E+01,.3906E+01,.3861E+01,.3875E+01,
     A .4057E+01,.4148E+01,.4206E+01,.4348E+01,.4558E+01,.4349E+01,
     A .4008E+01,.3864E+01,.3889E+01,.4079E+01,.4210E+01,.4197E+01,
     A .4183E+01,.4075E+01,.4220E+01,.4358E+01,.4493E+01,.4704E+01,
     A .4627E+01,.4481E+01,.4402E+01,.4335E+01,.4198E+01,.4208E+01,
     A .4280E+01,.4458E+01,.4626E+01,.4763E+01,.4787E+01,.4906E+01,
     A .4840E+01,.4665E+01,.4665E+01,.4567E+01,.4567E+01,.4451E+01,
     A .4430E+01,.4512E+01,.4599E+01,.4418E+01,.4367E+01,.4452E+01,
     A .4651E+01,.4650E+01,.4525E+01,.4336E+01,.3989E+01,.3514E+01,
     A .3239E+01,.3183E+01,.3147E+01,.3265E+01,.3441E+01,.3542E+01,
     A .3590E+01,.3521E+01,.3374E+01,.3254E+01,.3075E+01,.2988E+01,
     A .3226E+01,.3533E+01,.3826E+01,.4449E+01,.4839E+01,.4975E+01,
     A .5011E+01,.4975E+01,.5177E+01,.5594E+01,.5885E+01,.6024E+01,
     A .6088E+01,.6029E+01,.5909E+01,.5681E+01,.5531E+01,.5647E+01,
     A .5482E+01,.5395E+01,.5404E+01,.6044E+01,.6372E+01,.6329E+01,
     A .6277E+01,.5985E+01,.5663E+01,.5555E+01,.5489E+01,.5237E+01,
     A .4952E+01,.4923E+01,.5026E+01,.5734E+01,.6115E+01,.6018E+01,
     A .5783E+01,.5443E+01,.5661E+01,.5952E+01,.5954E+01,.5704E+01,
     A .5447E+01,.5649E+01,.5916E+01,.6161E+01,.6005E+01,.5731E+01,
     A .5574E+01,.5467E+01,.5472E+01,.5545E+01,.5623E+01,.5685E+01,
     A .5914E+01,.6206E+01,.6482E+01,.6523E+01,.6742E+01,.6828E+01,
     A .6896E+01,.6894E+01,.6908E+01,.6975E+01,.7019E+01,.6980E+01,
     A .6867E+01,.6455E+01,.6344E+01,.6412E+01,.6344E+01,.6429E+01,
     A .6541E+01,.6628E+01,.6403E+01,.6265E+01,.6327E+01,.6368E+01,
     A .6438E+01,.6337E+01,.6060E+01,.6050E+01,.6286E+01,.6359E+01,
     A .6151E+01,.6283E+01,.6273E+01,.6298E+01,.6536E+01,.7039E+01,
     A .7166E+01,.7007E+01,.6584E+01,.6399E+01,.6475E+01,.6613E+01,
     A .6873E+01,.6607E+01,.6265E+01,.5902E+01,.5887E+01,.6191E+01,
     A .6295E+01,.6026E+01,.5964E+01,.6231E+01,.6452E+01,.6503E+01,
     A .6445E+01,.6460E+01,.6690E+01,.6928E+01,.7454E+01,.7723E+01,
     A .7185E+01,.7028E+01,.6893E+01,.7015E+01,.7033E+01,.6970E+01,
     A .6666E+01,.6697E+01,.7196E+01,.7468E+01,.7233E+01,.7069E+01,
     A .7167E+01,.7147E+01,.7012E+01,.7168E+01,.7659E+01,.7637E+01,
     A .7627E+01,.7720E+01,.7593E+01,.7256E+01,.7402E+01,.7635E+01,
     A .8142E+01,.7902E+01,.7757E+01,.7465E+01,.7520E+01,.7679E+01,
     A .7690E+01,.7833E+01,.8122E+01,.8209E+01,.7675E+01,.7246E+01,
     A .7409E+01,.8272E+01,.8713E+01,.8810E+01,.8517E+01,.8184E+01,
     A .8039E+01,.8036E+01,.7962E+01,.8101E+01,.8346E+01,.8790E+01,
     A .9012E+01,.9256E+01,.8948E+01,.8645E+01,.8723E+01,.8774E+01,
     A .8559E+01,.8236E+01,.7856E+01,.7755E+01,.8223E+01,.8920E+01,
     A .9122E+01,.8974E+01,.8261E+01,.7110E+01,.6572E+01,.6153E+01,
     A .7002E+01,.8460E+01,.9075E+01,.9287E+01,.9265E+01,.8997E+01,
     A .8949E+01,.9225E+01,.9759E+01,.1030E+02,.1028E+02,.9968E+01,
     A .9471E+01,.8660E+01,.8113E+01,.8142E+01,.8387E+01,.9127E+01,
     A .9909E+01,.1022E+02,.1026E+02,.1019E+02,.1022E+02,.1011E+02,
     A .9539E+01,.8931E+01,.8792E+01,.9008E+01,.9691E+01,.1026E+02,
     A .1081E+02,.1124E+02,.1131E+02,.1109E+02,.1031E+02,.9922E+01,
     A .1005E+02,.1036E+02,.1116E+02,.1153E+02,.1116E+02,.1056E+02,
     A .9972E+01,.1034E+02,.1067E+02,.1073E+02,.1060E+02,.1096E+02,
     A .1113E+02,.1104E+02,.1105E+02,.1119E+02,.1133E+02,.1153E+02,
     A .1152E+02,.1143E+02,.1102E+02,.1071E+02,.9844E+01,.1014E+02,
     A .1120E+02,.1143E+02,.1117E+02,.1001E+02,.9846E+01,.1059E+02,
     A .1139E+02,.1161E+02,.1162E+02,.1142E+02,.1072E+02,.9948E+01,
     A .1027E+02,.1121E+02,.1151E+02,.1143E+02,.1171E+02,.1182E+02,
     A .1139E+02,.1114E+02,.1163E+02,.1228E+02,.1235E+02,.1231E+02,
     A .1198E+02,.1180E+02,.1181E+02,.1195E+02,.1241E+02,.1286E+02,
     A .1264E+02,.1246E+02,.1239E+02,.1289E+02,.1342E+02,.1342E+02,
     A .1341E+02,.1359E+02,.1397E+02,.1414E+02,.1445E+02,.1435E+02,
     A .1409E+02,.1363E+02,.1345E+02,.1362E+02,.1427E+02,.1462E+02,
     A .1413E+02,.1371E+02,.1343E+02,.1309E+02,.1314E+02,.1365E+02,
     A .1424E+02,.1363E+02,.1321E+02,.1343E+02,.1563E+02,.1722E+02,
     A .1735E+02,.1724E+02,.1703E+02,.1958E+02,.2027E+02,.1867E+02,
     A .1617E+02,.1631E+02,.1695E+02,.1793E+02,.1877E+02,.1931E+02,
     A .2340E+02,.2493E+02,.2421E+02,.1861E+02,.1809E+02,.2086E+02,
     A .2320E+02,.2520E+02,.2353E+02,.2323E+02,.2459E+02,.2771E+02,
     A .2511E+02,.2311E+02,.2106E+02,.1958E+02,.2163E+02,.2775E+02,
     A .3252E+02,.3564E+02,.2944E+02,.2265E+02,.1773E+02,.2043E+02,
     A .3096E+02,.3732E+02,.4047E+02,.3824E+02,.3494E+02,.3365E+02,
     A .3320E+02,.2857E+02,.2347E+02,.1976E+02,.2011E+02,.2106E+02,
     A .2363E+02,.3299E+02,.3607E+02,.4067E+02,.4017E+02,.3151E+02,
     A .2890E+02,.2874E+02,.2708E+02,.2474E+02,.3116E+02,.3853E+02,
     A .4318E+02,.4679E+02,.4353E+02,.3881E+02,.3488E+02,.3790E+02,
     A .4765E+02,.5191E+02,.5268E+02,.4398E+02,.3387E+02,.3118E+02,
     A .2888E+02,.2559E+02,.1729E+02,.1381E+02,.1336E+02,.1613E+02,
     A .2603E+02,.3257E+02,.3605E+02,.3768E+02,.3723E+02,.3743E+02,
     A .3767E+02,.3844E+02,.3413E+02,.2724E+02,.2481E+02,.2633E+02,
     A .3072E+02,.3030E+02,.2944E+02,.3010E+02,.3334E+02,.3113E+02,
     A .2678E+02,.2518E+02,.2706E+02,.3222E+02,.3657E+02,.4307E+02,
     A .4383E+02,.3415E+02,.2568E+02,.2239E+02,.2425E+02,.3319E+02,
     A .4729E+02,.5324E+02,.5320E+02,.4998E+02,.4741E+02,.4926E+02,
     A .5037E+02,.4932E+02,.4522E+02,.4063E+02,.3584E+02,.3468E+02,
     A .3393E+02,.3443E+02,.3572E+02,.3743E+02,.4190E+02,.4536E+02,
     A .4721E+02,.4856E+02,.3432E+02,.3011E+02,.3281E+02,.3574E+02,
     A .3873E+02,.4605E+02,.5042E+02,.5043E+02,.4410E+02,.4005E+02,
     A .3801E+02,.3524E+02,.3363E+02,.2562E+02,.2291E+02,.2201E+02,
     A .2164E+02,.2453E+02,.2951E+02,.3591E+02,.4250E+02,.4619E+02,
     A .4174E+02,.3935E+02,.3955E+02,.4161E+02,.3859E+02,.3710E+02,
     A .3756E+02,.3693E+02,.2917E+02,.2230E+02,.1895E+02,.1940E+02,
     A .2693E+02,.3345E+02,.3244E+02,.3037E+02,.2751E+02,.3308E+02,
     A .3747E+02,.3758E+02,.3706E+02,.3983E+02,.4140E+02,.4172E+02,
     A .3470E+02,.2888E+02,.2540E+02,.2575E+02,.2843E+02,.3541E+02,
     A .3671E+02,.3336E+02,.2944E+02,.2537E+02,.3208E+02,.3938E+02,
     A .4938E+02,.5531E+02,.5595E+02,.5669E+02,.5753E+02,.5920E+02,
     A .6140E+02,.5830E+02,.5047E+02,.4155E+02,.4151E+02,.4990E+02,
     A .5399E+02,.5305E+02,.3465E+02,.2456E+02,.2775E+02,.3901E+02,
     A .4484E+02,.3862E+02,.3719E+02,.4297E+02,.4836E+02,.4424E+02,
     A .3994E+02,.3754E+02,.4875E+02,.6473E+02,.6742E+02,.6465E+02,
     A .5485E+02,.4375E+02,.3655E+02,.3354E+02,.3591E+02,.4483E+02,
     A .5429E+02,.5875E+02,.6214E+02,.6047E+02,.5496E+02,.4795E+02,
     A .3638E+02,.3063E+02,.4021E+02,.4959E+02,.5851E+02,.6654E+02,
     A .7084E+02,.6634E+02,.6185E+02,.5814E+02,.6132E+02,.5558E+02,
     A .5099E+02,.4522E+02,.4192E+02,.3541E+02,.3360E+02,.3702E+02,
     A .4058E+02,.3081E+02,.2503E+02,.2139E+02,.2124E+02,.2951E+02,
     A .3991E+02,.4828E+02,.5213E+02,.4494E+02,.4516E+02,.4853E+02,
     A .4748E+02,.3712E+02,.2604E+02,.2372E+02,.2468E+02,.3016E+02,
     A .3754E+02,.4480E+02,.5306E+02,.5436E+02,.5392E+02,.5708E+02,
     A .5984E+02,.5785E+02,.5080E+02,.4516E+02,.4708E+02,.5451E+02,
     A .5271E+02,.5030E+02,.5014E+02,.5029E+02,.4566E+02,.4316E+02,
     A .4325E+02,.4676E+02,.5635E+02,.5758E+02,.5534E+02,.5137E+02,
     A .4633E+02,.4348E+02,.4471E+02,.5293E+02,.7094E+02,.8425E+02,
     A .8534E+02,.8348E+02,.7560E+02,.7749E+02,.8078E+02,.8101E+02,
     A .7453E+02,.6808E+02,.6331E+02,.5976E+02,.5278E+02,.6308E+02,
     A .7235E+02,.7254E+02,.6903E+02,.5591E+02,.5339E+02,.5338E+02,
     A .5392E+02,.5267E+02,.4974E+02,.4713E+02,.4661E+02,.5589E+02,
     A .6481E+02,.6614E+02,.6632E+02,.6943E+02,.7600E+02,.7545E+02,
     A .7015E+02,.5350E+02,.5296E+02,.5608E+02,.5566E+02,.5194E+02,
     A .5157E+02,.6011E+02,.6452E+02,.5897E+02,.4505E+02,.4683E+02,
     A .5527E+02,.5780E+02,.5153E+02,.4562E+02,.4689E+02,.5083E+02,
     A .5478E+02,.5342E+02,.5475E+02,.5873E+02,.5246E+02,.3948E+02,
     A .3758E+02,.4209E+02,.4795E+02,.5289E+02,.5080E+02,.4838E+02,
     A .4280E+02,.4018E+02,.3560E+02,.3447E+02,.3281E+02,.3426E+02,
     A .3572E+02,.3500E+02,.3364E+02,.3323E+02,.3290E+02,.3170E+02,
     A .3104E+02,.3162E+02,.3147E+02,.3298E+02,.3403E+02,.3519E+02,
     A .3733E+02,.3846E+02,.3748E+02,.3492E+02,.3738E+02,.4067E+02,
     A .3797E+02,.3668E+02,.3967E+02,.4467E+02,.4643E+02,.4678E+02,
     A .4765E+02,.4880E+02,.4791E+02,.4603E+02,.4946E+02,.5925E+02,
     A .6368E+02,.6616E+02,.6667E+02,.5919E+02,.5033E+02,.5166E+02,
     A .5482E+02,.6662E+02,.6907E+02,.6237E+02,.5087E+02,.3983E+02,
     A .3830E+02,.4658E+02,.5374E+02,.5573E+02,.4894E+02,.4377E+02,
     A .4255E+02,.4408E+02,.4795E+02,.4963E+02,.5068E+02,.5247E+02,
     A .5451E+02,.5323E+02,.5268E+02,.5521E+02,.6284E+02,.6423E+02,
     A .5730E+02,.4531E+02,.3313E+02,.3176E+02,.3321E+02,.3743E+02,
     A .4125E+02,.4120E+02,.4004E+02,.4094E+02,.4860E+02,.5302E+02,
     A .5165E+02,.4631E+02,.4284E+02,.4975E+02,.5755E+02,.6434E+02,
     A .7432E+02,.8220E+02,.7850E+02,.7698E+02,.6929E+02,.6490E+02,
     A .5521E+02,.5595E+02,.5058E+02,.4424E+02,.4066E+02,.3579E+02,
     A .3448E+02,.3508E+02,.3577E+02,.3521E+02,.3513E+02,.4367E+02,
     A .5348E+02,.6757E+02,.6996E+02,.6480E+02,.5648E+02,.5538E+02,
     A .5537E+02,.5874E+02,.6225E+02,.5468E+02,.4673E+02,.3766E+02,
     A .3574E+02,.4032E+02,.4137E+02,.4864E+02,.6411E+02,.8043E+02,
     A .8713E+02,.8609E+02,.7458E+02,.5667E+02,.4339E+02,.3768E+02,
     A .3731E+02,.4008E+02,.4179E+02,.4769E+02,.5530E+02,.6155E+02,
     A .6202E+02,.5813E+02,.5416E+02,.4698E+02,.3927E+02,.3750E+02,
     A .3586E+02,.3369E+02,.3955E+02,.4632E+02,.4940E+02,.5429E+02,
     A .6164E+02,.6440E+02,.5854E+02,.4549E+02,.3306E+02,.2793E+02,
     A .2776E+02,.3046E+02,.3663E+02,.4120E+02,.4477E+02,.4725E+02,
     A .5222E+02,.5324E+02,.4725E+02,.3981E+02,.2743E+02,.2379E+02,
     A .2589E+02,.2979E+02,.3105E+02,.3249E+02,.3622E+02,.3782E+02,
     A .3559E+02,.2891E+02,.3002E+02,.3308E+02,.4599E+02,.5355E+02,
     A .5477E+02,.5430E+02,.5589E+02,.6123E+02,.6625E+02,.6717E+02,
     A .6595E+02,.5709E+02,.4902E+02,.4539E+02,.4879E+02,.6297E+02,
     A .6814E+02,.6781E+02,.6669E+02,.6518E+02,.5525E+02,.4140E+02,
     A .3508E+02,.3039E+02,.3306E+02,.4200E+02,.4561E+02,.4797E+02,
     A .5603E+02,.6438E+02,.6244E+02,.5214E+02,.4195E+02,.3867E+02,
     A .4161E+02,.4549E+02,.4300E+02,.4015E+02,.4248E+02,.4433E+02,
     A .4639E+02,.5067E+02,.5604E+02,.5819E+02,.5765E+02,.6292E+02,
     A .6953E+02,.7197E+02,.6601E+02,.4988E+02,.3598E+02,.3511E+02,
     A .3887E+02,.4188E+02,.4655E+02,.5357E+02,.5923E+02,.6884E+02,
     A .6639E+02,.6062E+02,.5043E+02,.4219E+02,.3940E+02,.4034E+02,
     A .3860E+02,.3363E+02,.2629E+02,.2370E+02,.2254E+02,.2409E+02,
     A .2764E+02,.3359E+02,.3857E+02,.4353E+02,.5165E+02,.5641E+02,
     A .5655E+02,.5542E+02,.4922E+02,.3734E+02,.3035E+02,.2953E+02,
     A .3739E+02,.5329E+02,.6187E+02,.6522E+02,.6034E+02,.5708E+02,
     A .5840E+02,.5878E+02,.5558E+02,.4813E+02,.4009E+02,.3277E+02,
     A .2517E+02,.2509E+02,.3256E+02,.3900E+02,.4520E+02,.4840E+02,
     A .4277E+02,.3623E+02,.3306E+02,.3733E+02,.4906E+02,.5318E+02,
     A .5342E+02,.5017E+02,.4932E+02,.4956E+02,.4942E+02,.4270E+02,
     A .2958E+02,.2618E+02,.3064E+02,.3467E+02,.3328E+02,.3320E+02,
     A .3776E+02,.4525E+02,.5103E+02,.4929E+02,.4906E+02,.4243E+02,
     A .3275E+02,.2694E+02,.3302E+02,.3736E+02,.4220E+02,.4052E+02,
     A .3828E+02,.4176E+02,.4893E+02,.5442E+02,.5660E+02,.5757E+02,
     A .6159E+02,.6287E+02,.6165E+02,.6327E+02,.6919E+02,.7013E+02,
     A .6900E+02,.7203E+02,.7462E+02,.7695E+02,.8121E+02,.8453E+02,
     A .8394E+02,.7299E+02,.6376E+02,.6243E+02,.6164E+02,.6343E+02,
     A .6803E+02,.7547E+02,.8241E+02,.8973E+02,.8435E+02,.7442E+02,
     A .7017E+02,.6517E+02,.6824E+02,.7542E+02,.7758E+02,.7907E+02,
     A .7606E+02,.6735E+02,.6221E+02,.6078E+02,.6516E+02,.6919E+02,
     A .6534E+02,.5371E+02,.5308E+02,.5814E+02,.6076E+02,.6406E+02,
     A .7192E+02,.7439E+02,.6681E+02,.6334E+02,.5602E+02,.5864E+02,
     A .6658E+02,.7258E+02,.7678E+02,.7883E+02,.7827E+02,.7369E+02,
     A .6655E+02,.6194E+02,.6367E+02,.6500E+02,.6552E+02,.6271E+02,
     A .6004E+02,.5518E+02,.5323E+02,.4977E+02,.5042E+02,.5478E+02,
     A .5593E+02,.5263E+02,.5011E+02,.4819E+02,.4746E+02,.4551E+02,
     A .4649E+02,.4742E+02,.4656E+02,.4662E+02,.4907E+02,.5316E+02,
     A .5282E+02,.5027E+02,.4809E+02,.5038E+02,.5338E+02,.5700E+02,
     A .5552E+02,.5094E+02,.5055E+02,.5095E+02,.4990E+02,.4645E+02,
     A .4501E+02,.4205E+02,.3961E+02,.4030E+02,.4564E+02,.4883E+02,
     A .5208E+02,.5292E+02,.5167E+02,.5156E+02,.5074E+02,.5211E+02,
     A .5595E+02,.5737E+02,.5753E+02,.5750E+02,.5635E+02,.5521E+02,
     A .5540E+02,.5878E+02,.6227E+02,.5458E+02,.5044E+02,.4242E+02,
     A .4952E+02,.6341E+02,.6597E+02,.6484E+02,.6443E+02,.6304E+02,
     A .6110E+02,.5983E+02,.5969E+02,.5601E+02,.5225E+02,.4961E+02,
     A .4345E+02,.4555E+02,.5302E+02,.6125E+02,.5995E+02,.5257E+02,
     A .4027E+02,.3543E+02,.3001E+02,.3244E+02,.4503E+02,.4755E+02,
     A .5538E+02,.5621E+02,.5192E+02,.4877E+02,.4323E+02,.3602E+02,
     A .3515E+02,.3682E+02,.3679E+02,.3442E+02,.3033E+02,.3017E+02,
     A .3339E+02,.4278E+02,.4992E+02,.5560E+02,.5636E+02,.6028E+02,
     A .6771E+02,.7413E+02,.7496E+02,.6884E+02,.6441E+02,.6748E+02,
     A .7059E+02,.7001E+02,.6774E+02,.7046E+02,.7444E+02,.7735E+02,
     A .7051E+02,.5140E+02,.4618E+02,.5086E+02,.6004E+02,.6514E+02,
     A .6768E+02,.6980E+02,.7199E+02,.6471E+02,.5804E+02,.4960E+02,
     A .4310E+02,.3771E+02,.4229E+02,.5056E+02,.6145E+02,.6706E+02,
     A .6950E+02,.6877E+02,.6482E+02,.5381E+02,.5229E+02,.5433E+02,
     A .5501E+02,.5338E+02,.4846E+02,.4438E+02,.3924E+02,.3891E+02,
     A .3912E+02,.3588E+02,.3428E+02,.3697E+02,.4037E+02,.4044E+02,
     A .3989E+02,.3934E+02,.4013E+02,.4673E+02,.5147E+02,.5319E+02,
     A .5001E+02,.4500E+02,.4026E+02,.3685E+02,.3536E+02,.3572E+02,
     A .3671E+02,.4036E+02,.4587E+02,.4910E+02,.5003E+02,.4837E+02,
     A .4227E+02,.3648E+02,.3433E+02,.3269E+02,.3417E+02,.4412E+02,
     A .5043E+02,.5908E+02,.6269E+02,.6077E+02,.5835E+02,.5733E+02,
     A .5702E+02,.5623E+02,.5729E+02,.5607E+02,.5317E+02,.5058E+02,
     A .5161E+02,.5324E+02,.5584E+02,.5867E+02,.5473E+02,.5377E+02,
     A .5748E+02,.6353E+02,.5885E+02,.5407E+02,.5380E+02,.5782E+02,
     A .6113E+02,.6742E+02,.6993E+02,.6832E+02,.6529E+02,.6284E+02,
     A .6124E+02,.5712E+02,.5451E+02,.5886E+02,.6178E+02,.6190E+02,
     A .5583E+02,.5012E+02,.4821E+02,.5230E+02,.6253E+02,.7651E+02,
     A .8560E+02,.9392E+02,.9243E+02,.7992E+02,.8516E+02,.9309E+02,
     A .9752E+02,.8947E+02,.8220E+02,.8345E+02,.8601E+02,.8143E+02,
     A .7935E+02,.8424E+02,.9403E+02,.1052E+03,.1058E+03,.9903E+02,
     A .9685E+02,.9574E+02,.8828E+02,.7864E+02,.7098E+02,.6508E+02,
     A .6887E+02,.7736E+02,.8624E+02,.9606E+02,.1132E+03,.1284E+03,
     A .1283E+03,.1224E+03,.1185E+03,.1201E+03,.1195E+03,.1226E+03,
     A .1352E+03,.1426E+03,.1413E+03,.1331E+03,.1319E+03,.1367E+03,
     A .1364E+03,.1428E+03,.1550E+03,.1503E+03,.1336E+03,.1089E+03,
     A .8521E+02,.6881E+02,.7376E+02,.8922E+02,.1028E+03,.1187E+03,
     A .1345E+03,.1390E+03,.1488E+03,.1595E+03,.1689E+03,.1772E+03,
     A .1953E+03,.1820E+03,.1290E+03,.1137E+03,.1244E+03,.1320E+03,
     A .1081E+03,.8765E+02,.7623E+02,.6129E+02,.7103E+02,.9396E+02,
     A .1110E+03,.1180E+03,.1340E+03,.1676E+03,.1864E+03,.1947E+03,
     A .1613E+03,.1362E+03,.1154E+03,.1160E+03,.1255E+03,.1104E+03,
     A .9631E+02,.7680E+02,.9063E+02,.1160E+03,.1255E+03,.1218E+03,
     A .1142E+03,.1031E+03,.8989E+02,.7208E+02,.5276E+02,.4122E+02,
     A .4110E+02,.4725E+02,.6363E+02,.9057E+02,.1119E+03,.1218E+03,
     A .1272E+03,.1268E+03,.1223E+03,.1155E+03,.1121E+03,.1002E+03,
     A .8104E+02,.7191E+02,.6451E+02,.5719E+02,.6288E+02,.7703E+02,
     A .8859E+02,.9974E+02,.1057E+03,.1081E+03,.1067E+03,.8739E+02,
     A .6538E+02,.5545E+02,.5773E+02,.6674E+02,.6894E+02,.6506E+02,
     A .6918E+02,.9121E+02,.1164E+03,.1268E+03,.1243E+03,.1155E+03,
     A .9024E+02,.7794E+02,.8261E+02,.1005E+03,.1284E+03,.1340E+03,
     A .1263E+03,.1156E+03,.1102E+03,.1122E+03,.1306E+03,.1383E+03,
     A .1355E+03,.1388E+03,.1412E+03,.1292E+03,.1035E+03,.8556E+02,
     A .9319E+02,.1057E+03,.1086E+03,.9213E+02,.8300E+02,.8730E+02,
     A .9729E+02,.9734E+02,.8167E+02,.6837E+02,.6072E+02,.6770E+02,
     A .8646E+02,.1059E+03,.1277E+03,.1500E+03,.1760E+03,.1785E+03,
     A .1759E+03,.1803E+03,.2094E+03,.2221E+03,.2381E+03,.2723E+03,
     A .2986E+03,.2988E+03,.2852E+03,.2571E+03,.2408E+03,.2352E+03,
     A .2329E+03,.2506E+03,.3148E+03,.3534E+03,.3528E+03,.3357E+03,
     A .2920E+03,.2688E+03,.2634E+03,.2453E+03,.2373E+03,.2367E+03,
     A .2311E+03,.2104E+03,.1877E+03,.1808E+03,.2010E+03,.2154E+03,
     A .2209E+03,.2191E+03,.2216E+03,.2287E+03,.2396E+03,.2565E+03,
     A .2570E+03,.2573E+03,.2922E+03,.3301E+03,.3230E+03,.2914E+03,
     A .2757E+03,.2944E+03,.3014E+03,.2868E+03,.2670E+03,.2610E+03,
     A .2897E+03,.3009E+03,.2960E+03,.2820E+03,.2447E+03,.2307E+03,
     A .2440E+03,.2557E+03,.2231E+03,.2066E+03,.2201E+03,.2690E+03,
     A .2880E+03,.2566E+03,.2120E+03,.2117E+03,.2651E+03,.2898E+03,
     A .2636E+03,.2510E+03,.2229E+03,.2351E+03,.2778E+03,.3029E+03,
     A .2928E+03,.2026E+03,.1805E+03,.2013E+03,.2139E+03,.3016E+03,
     A .3203E+03,.3282E+03,.3214E+03,.3088E+03,.2970E+03,.2753E+03,
     A .2466E+03,.2364E+03,.2119E+03,.2069E+03,.2129E+03,.2331E+03,
     A .2505E+03,.2632E+03,.2717E+03,.2726E+03,.2665E+03,.2596E+03,
     A .2513E+03,.2461E+03,.2379E+03,.2468E+03,.2601E+03,.2760E+03,
     A .2807E+03,.2791E+03,.2735E+03,.2539E+03,.2399E+03,.2272E+03,
     A .2131E+03,.2041E+03,.2141E+03,.2305E+03,.2468E+03,.2581E+03,
     A .2603E+03,.2611E+03,.2686E+03,.2814E+03,.2629E+03,.2089E+03,
     A .1796E+03,.1851E+03,.2112E+03,.2460E+03,.2523E+03,.2347E+03,
     A .2260E+03,.2323E+03,.2990E+03,.3365E+03,.3458E+03,.3047E+03,
     A .2585E+03,.2496E+03,.2761E+03,.2997E+03,.3008E+03,.2953E+03,
     A .3116E+03,.3484E+03,.3202E+03,.2681E+03,.2490E+03,.2461E+03,
     A .2653E+03,.2651E+03,.2546E+03,.2468E+03,.2646E+03,.2738E+03,
     A .2693E+03,.2623E+03,.2478E+03,.2731E+03,.2901E+03,.3255E+03,
     A .3047E+03,.2247E+03,.1904E+03,.1860E+03,.2101E+03,.2524E+03,
     A .2579E+03,.2650E+03,.2469E+03,.2045E+03,.1554E+03,.1262E+03,
     A .1129E+03,.1255E+03,.1279E+03,.1391E+03,.1547E+03,.2003E+03,
     A .2133E+03,.2149E+03,.1968E+03,.1785E+03,.1853E+03,.1947E+03,
     A .2083E+03,.2180E+03,.2115E+03,.2023E+03,.1866E+03,.1985E+03,
     A .2574E+03,.2786E+03,.3316E+03,.3308E+03,.2793E+03,.2491E+03,
     A .2562E+03,.2663E+03,.2592E+03,.2297E+03,.1883E+03,.1811E+03,
     A .1957E+03,.2094E+03,.2107E+03,.2075E+03,.1830E+03,.1439E+03,
     A .1243E+03,.1481E+03,.1778E+03,.1591E+03,.1322E+03,.1180E+03,
     A .1406E+03,.1647E+03,.1621E+03,.1479E+03,.1171E+03,.9899E+02,
     A .9858E+02,.1044E+03,.1293E+03,.1640E+03,.1847E+03,.1704E+03,
     A .1180E+03,.1055E+03,.1144E+03,.1421E+03,.1529E+03,.1247E+03,
     A .9368E+02,.8662E+02,.1060E+03,.1468E+03,.1947E+03,.2245E+03,
     A .2408E+03,.2370E+03,.2119E+03,.1921E+03,.1745E+03,.1543E+03,
     A .1216E+03,.1040E+03,.1084E+03,.1468E+03,.1973E+03,.2502E+03,
     A .2912E+03,.3203E+03,.2862E+03,.2654E+03,.2652E+03,.2686E+03,
     A .2364E+03,.2010E+03,.1795E+03,.2009E+03,.2451E+03,.2685E+03,
     A .2832E+03,.3107E+03,.3159E+03,.2935E+03,.2627E+03,.2498E+03,
     A .2486E+03,.2430E+03,.2303E+03,.2202E+03,.2219E+03,.2522E+03,
     A .2793E+03,.2912E+03,.2890E+03,.2622E+03,.2407E+03,.2528E+03,
     A .2887E+03,.3035E+03,.2929E+03,.2784E+03,.2629E+03,.2466E+03,
     A .2262E+03,.2201E+03,.2261E+03,.2098E+03,.1788E+03,.1741E+03,
     A .1723E+03,.1702E+03,.1689E+03,.1708E+03,.1726E+03,.1802E+03,
     A .1840E+03,.1735E+03,.1676E+03,.1679E+03,.1797E+03,.1848E+03,
     A .1859E+03,.1853E+03,.1827E+03,.1811E+03,.1708E+03,.1490E+03,
     A .1388E+03,.1351E+03,.1339E+03,.1241E+03,.1161E+03,.1073E+03,
     A .1016E+03,.1009E+03,.9008E+02,.8593E+02,.7438E+02,.7001E+02,
     A .6397E+02,.7339E+02,.9721E+02,.1008E+03,.8343E+02,.6165E+02,
     A .6161E+02,.6489E+02,.7001E+02,.7997E+02,.8649E+02,.8916E+02,
     A .8558E+02,.7638E+02,.6955E+02,.7610E+02,.8231E+02,.8930E+02,
     A .7931E+02,.7377E+02,.7451E+02,.8829E+02,.1034E+03,.1114E+03,
     A .1151E+03,.1249E+03,.1355E+03,.1488E+03,.1561E+03,.1668E+03,
     A .1777E+03,.1825E+03,.1809E+03,.1781E+03,.1782E+03,.1805E+03,
     A .1812E+03,.1817E+03,.2038E+03,.2429E+03,.2501E+03,.2740E+03,
     A .2669E+03,.2599E+03,.2515E+03,.2475E+03,.2512E+03,.2649E+03,
     A .2752E+03,.2824E+03,.2981E+03,.3102E+03,.3101E+03,.3074E+03,
     A .2977E+03,.2945E+03,.2908E+03,.3043E+03,.3279E+03,.3524E+03,
     A .3593E+03,.3177E+03,.2822E+03,.2856E+03,.3011E+03,.3146E+03,
     A .3180E+03,.3187E+03,.3266E+03,.3354E+03,.3549E+03,.3701E+03,
     A .3611E+03,.3608E+03,.3566E+03,.3459E+03,.3288E+03,.3338E+03,
     A .3561E+03,.3453E+03,.3177E+03,.2807E+03,.2436E+03,.2471E+03,
     A .3227E+03,.3567E+03,.3717E+03,.3684E+03,.3579E+03,.3481E+03,
     A .3465E+03,.3234E+03,.2884E+03,.2845E+03,.3110E+03,.3401E+03,
     A .3229E+03,.2623E+03,.2215E+03,.2229E+03,.2786E+03,.3099E+03,
     A .2802E+03,.2551E+03,.2280E+03,.2166E+03,.1977E+03,.1871E+03,
     A .1581E+03,.1408E+03,.1264E+03,.1184E+03,.1010E+03,.8693E+02,
     A .6969E+02,.5835E+02,.4953E+02,.5788E+02,.7748E+02,.9271E+02,
     A .1105E+03,.1311E+03,.1439E+03,.1551E+03,.1696E+03,.1972E+03,
     A .2220E+03,.2245E+03,.2268E+03,.2583E+03,.3102E+03,.3312E+03,
     A .3396E+03,.3302E+03,.3368E+03,.3412E+03,.3367E+03,.3324E+03,
     A .3574E+03,.3717E+03,.3853E+03,.3719E+03,.3629E+03,.3785E+03,
     A .3807E+03,.3544E+03,.3165E+03,.3145E+03,.3123E+03,.3012E+03,
     A .3027E+03,.3476E+03,.3959E+03,.4310E+03,.4630E+03,.4622E+03,
     A .4377E+03,.4234E+03,.3810E+03,.3570E+03,.3664E+03,.3860E+03,
     A .3760E+03,.3404E+03,.3288E+03,.3442E+03,.3752E+03,.3683E+03,
     A .3390E+03,.3011E+03,.2630E+03,.2412E+03,.2094E+03,.1604E+03,
     A .1367E+03,.1292E+03,.1784E+03,.2215E+03,.2630E+03,.3218E+03,
     A .4017E+03,.4283E+03,.4560E+03,.4417E+03,.4267E+03,.4171E+03,
     A .4118E+03,.4063E+03,.4180E+03,.4244E+03,.4129E+03,.4055E+03,
     A .4173E+03,.4150E+03,.4141E+03,.4236E+03,.4266E+03,.4265E+03,
     A .4247E+03,.4262E+03,.4420E+03,.4698E+03,.5109E+03,.5446E+03,
     A .5536E+03,.5450E+03,.5460E+03,.5374E+03,.5371E+03,.5589E+03,
     A .6109E+03,.6543E+03,.6931E+03,.6761E+03,.6239E+03,.5735E+03,
     A .5718E+03,.6440E+03,.6836E+03,.6945E+03,.6932E+03,.6820E+03,
     A .6635E+03,.6362E+03,.6091E+03,.6119E+03,.5980E+03,.5902E+03,
     A .5739E+03,.5765E+03,.5831E+03,.5918E+03,.5841E+03,.5800E+03,
     A .5935E+03,.6028E+03,.6370E+03,.7088E+03,.7283E+03,.6645E+03,
     A .6002E+03,.5658E+03,.5320E+03,.5516E+03,.5667E+03,.5878E+03,
     A .5576E+03,.5266E+03,.5228E+03,.5693E+03,.5919E+03,.5948E+03,
     A .5821E+03,.5829E+03,.6027E+03,.6011E+03,.5931E+03,.5665E+03,
     A .4992E+03,.4823E+03,.4962E+03,.5076E+03,.4796E+03,.4408E+03,
     A .4419E+03,.4591E+03,.5053E+03,.4928E+03,.4255E+03,.3932E+03,
     A .4395E+03,.6074E+03,.6589E+03,.6638E+03,.6531E+03,.6821E+03,
     A .6699E+03,.6407E+03,.5983E+03,.6037E+03,.6557E+03,.6883E+03,
     A .6982E+03,.5994E+03,.4302E+03,.3431E+03,.3516E+03,.4227E+03,
     A .4759E+03,.4804E+03,.4784E+03,.5302E+03,.5575E+03,.5203E+03,
     A .4724E+03,.4823E+03,.6050E+03,.6422E+03,.6193E+03,.5562E+03,
     A .5131E+03,.5186E+03,.5709E+03,.6615E+03,.6829E+03,.5779E+03,
     A .4601E+03,.3460E+03,.3070E+03,.3481E+03,.3948E+03,.4546E+03,
     A .5415E+03,.6070E+03,.6480E+03,.6789E+03,.6944E+03,.6719E+03,
     A .5511E+03,.4434E+03,.4566E+03,.5951E+03,.6437E+03,.6513E+03,
     A .5979E+03,.5002E+03,.4904E+03,.5248E+03,.6584E+03,.6877E+03,
     A .6532E+03,.5960E+03,.5563E+03,.5645E+03,.5870E+03,.6399E+03,
     A .6858E+03,.7378E+03,.7198E+03,.6762E+03,.5990E+03,.4404E+03,
     A .3660E+03,.3327E+03,.3107E+03,.2817E+03,.3633E+03,.4463E+03,
     A .5014E+03,.5012E+03,.3955E+03,.3346E+03,.3179E+03,.4131E+03,
     A .5371E+03,.5631E+03,.5162E+03,.4162E+03,.3368E+03,.4378E+03,
     A .5345E+03,.6456E+03,.6557E+03,.6141E+03,.6060E+03,.6855E+03,
     A .8106E+03,.8338E+03,.7915E+03,.6751E+03,.5880E+03,.5659E+03,
     A .5569E+03,.5113E+03,.4115E+03,.3679E+03,.4052E+03,.4293E+03,
     A .3605E+03,.2919E+03,.2888E+03,.3069E+03,.3168E+03,.3545E+03,
     A .4143E+03,.4364E+03,.4282E+03,.4402E+03,.4881E+03,.5396E+03,
     A .5977E+03,.6790E+03,.7260E+03,.7113E+03,.6773E+03,.6046E+03,
     A .5762E+03,.5757E+03,.5665E+03,.4365E+03,.3412E+03,.3105E+03,
     A .3744E+03,.5318E+03,.6173E+03,.6478E+03,.6810E+03,.7001E+03,
     A .6205E+03,.5251E+03,.4172E+03,.3286E+03,.2879E+03,.2836E+03,
     A .2931E+03,.3257E+03,.3107E+03,.2909E+03,.2918E+03,.4177E+03,
     A .5294E+03,.5869E+03,.5987E+03,.5936E+03,.4916E+03,.4515E+03,
     A .3964E+03,.3660E+03,.3545E+03,.3858E+03,.4276E+03,.4724E+03,
     A .5786E+03,.6383E+03,.5885E+03,.5199E+03,.5041E+03,.5760E+03,
     A .5958E+03,.5880E+03,.5798E+03,.5738E+03,.5769E+03,.5688E+03,
     A .5602E+03,.5031E+03,.4645E+03,.4098E+03,.3871E+03,.3596E+03,
     A .3395E+03,.2891E+03,.2240E+03,.1755E+03,.2960E+03,.3169E+03,
     A .5179E+03,.6025E+03,.6408E+03,.6164E+03,.5608E+03,.5126E+03,
     A .4503E+03,.3964E+03,.4582E+03,.5913E+03,.7519E+03,.8129E+03,
     A .8330E+03,.7430E+03,.6791E+03,.6367E+03,.6280E+03,.6052E+03,
     A .5913E+03,.6721E+03,.7162E+03,.7974E+03,.8081E+03,.7616E+03,
     A .7480E+03,.7709E+03,.7931E+03,.7986E+03,.7621E+03,.5723E+03,
     A .4194E+03,.4217E+03,.4380E+03,.5779E+03,.6705E+03,.6461E+03,
     A .6250E+03,.6052E+03,.5210E+03,.4966E+03,.5534E+03,.5330E+03,
     A .6664E+03,.7387E+03,.7518E+03,.7425E+03,.7586E+03,.7465E+03,
     A .7544E+03,.5987E+03,.4935E+03,.4320E+03,.4723E+03,.5859E+03,
     A .7036E+03,.7289E+03,.6667E+03,.6448E+03,.7237E+03,.8742E+03,
     A .8306E+03,.8642E+03,.7996E+03,.5934E+03,.5643E+03,.6223E+03,
     A .6363E+03,.6681E+03,.6655E+03,.5398E+03,.3972E+03,.3480E+03,
     A .3338E+03,.3492E+03,.3822E+03,.4433E+03,.6129E+03,.6784E+03,
     A .7147E+03,.7505E+03,.7413E+03,.7151E+03,.6960E+03,.5862E+03,
     A .5345E+03,.5748E+03,.6135E+03,.6146E+03,.4967E+03,.4543E+03,
     A .4170E+03,.4616E+03,.5525E+03,.6434E+03,.7115E+03,.7460E+03,
     A .7447E+03,.7263E+03,.7271E+03,.6417E+03,.6404E+03,.6708E+03,
     A .6909E+03,.7060E+03,.7113E+03,.6082E+03,.5593E+03,.5268E+03,
     A .6264E+03,.7289E+03,.7337E+03,.6894E+03,.6315E+03,.6599E+03,
     A .7494E+03,.7740E+03,.7442E+03,.7021E+03,.5904E+03,.5391E+03,
     A .5278E+03,.5444E+03,.5904E+03,.5840E+03,.7286E+03,.8250E+03,
     A .8448E+03,.8154E+03,.8040E+03,.6977E+03,.6107E+03,.5828E+03,
     A .5830E+03,.6268E+03,.6576E+03,.6275E+03,.5976E+03,.5535E+03,
     A .4920E+03,.4653E+03,.4531E+03,.3950E+03,.3828E+03,.3920E+03,
     A .5111E+03,.6473E+03,.6747E+03,.5894E+03,.5554E+03,.4616E+03,
     A .4274E+03,.5223E+03,.5915E+03,.6119E+03,.5349E+03,.3989E+03,
     A .3112E+03,.2902E+03,.3153E+03,.3540E+03,.4488E+03,.5404E+03,
     A .6738E+03,.7725E+03,.7880E+03,.7711E+03,.7502E+03,.7316E+03,
     A .7665E+03,.7769E+03,.8233E+03,.8616E+03,.9103E+03,.9126E+03,
     A .9008E+03,.8157E+03,.7629E+03,.7268E+03,.7633E+03,.8125E+03,
     A .8557E+03,.8829E+03,.8895E+03,.8128E+03,.7448E+03,.7520E+03,
     A .8064E+03,.8276E+03,.7839E+03,.6732E+03,.6465E+03,.6523E+03,
     A .6398E+03,.6393E+03,.6285E+03,.6094E+03,.6032E+03,.6072E+03,
     A .6091E+03,.5980E+03,.6196E+03,.6964E+03,.7318E+03,.7963E+03,
     A .8055E+03,.7412E+03,.6640E+03,.6159E+03,.6359E+03,.6979E+03,
     A .7349E+03,.7512E+03,.7347E+03,.7268E+03,.7339E+03,.7565E+03,
     A .7378E+03,.7069E+03,.7003E+03,.7070E+03,.7022E+03,.7016E+03,
     A .6936E+03,.6711E+03,.6370E+03,.6778E+03,.7311E+03,.7522E+03,
     A .7298E+03,.7455E+03,.7820E+03,.8210E+03,.8846E+03,.9064E+03,
     A .8719E+03,.8385E+03,.8145E+03,.8315E+03,.8745E+03,.8675E+03,
     A .8637E+03,.6584E+03,.5241E+03,.4359E+03,.3730E+03,.3695E+03,
     A .5342E+03,.6023E+03,.6355E+03,.6615E+03,.6460E+03,.6676E+03,
     A .6842E+03,.8125E+03,.8391E+03,.8086E+03,.7858E+03,.7482E+03,
     A .7698E+03,.7606E+03,.7591E+03,.7768E+03,.7840E+03,.7875E+03,
     A .7475E+03,.7220E+03,.7256E+03,.7274E+03,.6846E+03,.6652E+03,
     A .5111E+03,.4472E+03,.4358E+03,.4282E+03,.4395E+03,.4757E+03,
     A .4912E+03,.4901E+03,.4865E+03,.5175E+03,.5718E+03,.6448E+03,
     A .6933E+03,.7143E+03,.6955E+03,.6528E+03,.5316E+03,.6799E+03,
     A .6972E+03,.6989E+03,.7470E+03,.7715E+03,.7141E+03,.6800E+03,
     A .6951E+03,.7639E+03,.8585E+03,.9064E+03,.9505E+03,.9666E+03,
     A .9467E+03,.9562E+03,.9384E+03,.9085E+03,.8604E+03,.8401E+03,
     A .8347E+03,.8812E+03,.9446E+03,.9393E+03,.9483E+03,.8606E+03,
     A .7224E+03,.6490E+03,.5800E+03,.5229E+03,.5281E+03,.5539E+03,
     A .5470E+03,.5367E+03,.5559E+03,.5954E+03,.6038E+03,.6280E+03,
     A .6522E+03,.6662E+03,.6912E+03,.7158E+03,.7510E+03,.7378E+03,
     A .7291E+03,.7514E+03,.7787E+03,.8587E+03,.9643E+03,.9757E+03,
     A .8688E+03,.7628E+03,.6855E+03,.6427E+03,.5367E+03,.4773E+03,
     A .5401E+03,.5766E+03,.6606E+03,.7134E+03,.7185E+03,.7217E+03,
     A .7199E+03,.6906E+03,.6958E+03,.7924E+03,.8242E+03,.9081E+03,
     A .9575E+03,.8809E+03,.7491E+03,.7630E+03,.8399E+03,.8974E+03,
     A .8753E+03,.8458E+03,.8754E+03,.8994E+03,.9703E+03,.9845E+03,
     A .9255E+03,.9060E+03,.9358E+03,.9572E+03,.9616E+03,.9166E+03,
     A .8395E+03,.7650E+03,.7419E+03,.6974E+03,.7072E+03,.6443E+03,
     A .6039E+03,.5974E+03,.5925E+03,.6025E+03,.6138E+03,.5632E+03,
     A .5601E+03,.6134E+03,.6927E+03,.7540E+03,.7609E+03,.7311E+03,
     A .6956E+03,.7505E+03,.8491E+03,.9217E+03,.8997E+03,.8708E+03,
     A .8603E+03,.9111E+03,.8653E+03,.7610E+03,.6794E+03,.6575E+03,
     A .6909E+03,.8178E+03,.8486E+03,.8311E+03,.7478E+03,.6404E+03,
     A .6167E+03,.6373E+03,.6508E+03,.6122E+03,.5734E+03,.5556E+03,
     A .5442E+03,.5457E+03,.6071E+03,.6183E+03,.7068E+03,.7885E+03,
     A .8002E+03,.8018E+03,.7580E+03,.6847E+03,.5937E+03,.6012E+03,
     A .6348E+03,.6942E+03,.6544E+03,.6185E+03,.6570E+03,.7563E+03,
     A .7954E+03,.8078E+03,.7138E+03,.6984E+03,.8442E+03,.8943E+03,
     A .9037E+03,.8434E+03,.8180E+03,.8088E+03,.7876E+03,.7903E+03,
     A .7766E+03,.8367E+03,.8710E+03,.9011E+03,.9282E+03,.8693E+03,
     A .8438E+03,.7417E+03,.7017E+03,.7183E+03,.7443E+03,.7698E+03,
     A .7732E+03,.7647E+03,.7737E+03,.8087E+03,.8347E+03,.8425E+03,
     A .8357E+03,.8265E+03,.8055E+03,.9406E+03,.1005E+04,.1082E+04,
     A .1055E+04,.1029E+04,.9936E+03,.9878E+03,.9994E+03,.1080E+04,
     A .1110E+04,.1111E+04,.1104E+04,.1052E+04,.1050E+04,.1007E+04,
     A .1011E+04,.1057E+04,.1070E+04,.1017E+04,.9238E+03,.8702E+03,
     A .8732E+03,.9255E+03,.9877E+03,.1085E+04,.1093E+04,.1073E+04,
     A .1080E+04,.1073E+04,.1065E+04,.1040E+04,.1017E+04,.9685E+03,
     A .8825E+03,.8864E+03,.9564E+03,.9962E+03,.1055E+04,.1013E+04,
     A .1014E+04,.1035E+04,.1076E+04,.1055E+04,.1012E+04,.9393E+03,
     A .9481E+03,.9591E+03,.9959E+03,.9750E+03,.9384E+03,.9181E+03,
     A .9285E+03,.9690E+03,.1004E+04,.9723E+03,.9537E+03,.9061E+03,
     A .9495E+03,.9849E+03,.1010E+04,.1021E+04,.9360E+03,.9111E+03,
     A .8382E+03,.7967E+03,.7368E+03,.7951E+03,.8579E+03,.9216E+03,
     A .1010E+04,.1105E+04,.1159E+04,.1140E+04,.1106E+04,.1080E+04,
     A .1066E+04,.1086E+04,.1114E+04,.1186E+04,.1196E+04,.1190E+04,
     A .1130E+04,.1106E+04,.1129E+04,.1144E+04,.1164E+04,.1153E+04,
     A .1126E+04,.1101E+04,.1113E+04,.1192E+04,.1303E+04,.1343E+04,
     A .1280E+04,.1184E+04,.1061E+04,.9664E+03,.9691E+03,.1039E+04,
     A .1087E+04,.1007E+04,.1005E+04,.8930E+03,.7938E+03,.7762E+03,
     A .8475E+03,.9245E+03,.1023E+04,.1060E+04,.1073E+04,.1084E+04,
     A .1061E+04,.1049E+04,.1064E+04,.1064E+04,.1037E+04,.9907E+03,
     A .9871E+03,.9564E+03,.9282E+03,.8761E+03,.8760E+03,.9163E+03,
     A .1002E+04,.1058E+04,.1065E+04,.1073E+04,.1093E+04,.1089E+04,
     A .1093E+04,.1066E+04,.1032E+04,.1025E+04,.1005E+04,.1035E+04,
     A .9213E+03,.8364E+03,.7921E+03,.7941E+03,.8952E+03,.9676E+03,
     A .1072E+04,.1102E+04,.1130E+04,.1054E+04,.9987E+03,.9781E+03,
     A .9965E+03,.9521E+03,.9051E+03,.8625E+03,.8725E+03,.9932E+03,
     A .1055E+04,.1082E+04,.1027E+04,.1024E+04,.1059E+04,.1090E+04,
     A .1018E+04,.9326E+03,.8855E+03,.8027E+03,.7631E+03,.7400E+03,
     A .7991E+03,.8355E+03,.9292E+03,.9910E+03,.1004E+04,.9533E+03,
     A .9775E+03,.9746E+03,.9982E+03,.9633E+03,.9481E+03,.9044E+03,
     A .1018E+04,.1070E+04,.1103E+04,.1098E+04,.1094E+04,.1102E+04,
     A .1116E+04,.1091E+04,.1076E+04,.1059E+04,.1043E+04,.8910E+03,
     A .7495E+03,.6864E+03,.7627E+03,.8717E+03,.1069E+04,.1102E+04,
     A .1111E+04,.1105E+04,.1089E+04,.1035E+04,.1018E+04,.1042E+04,
     A .1109E+04,.1175E+04,.1153E+04,.1148E+04,.1120E+04,.1081E+04,
     A .1016E+04,.8497E+03,.6928E+03,.5711E+03,.5098E+03,.4860E+03,
     A .4643E+03,.5135E+03,.6580E+03,.7003E+03,.8984E+03,.9819E+03,
     A .1023E+04,.1043E+04,.1104E+04,.1069E+04,.9964E+03,.9353E+03,
     A .8974E+03,.8080E+03,.7741E+03,.7656E+03,.7056E+03,.6721E+03,
     A .7047E+03,.7429E+03,.7588E+03,.7260E+03,.6942E+03,.6899E+03,
     A .7160E+03,.8145E+03,.8865E+03,.8933E+03,.9954E+03,.1094E+04,
     A .1152E+04,.1154E+04,.1145E+04,.1135E+04,.1045E+04,.9241E+03,
     A .8769E+03,.7815E+03,.6515E+03,.6408E+03,.6814E+03,.8950E+03,
     A .1011E+04,.1010E+04,.9785E+03,.9158E+03,.9287E+03,.9662E+03,
     A .1035E+04,.1054E+04,.1097E+04,.1154E+04,.1112E+04,.1070E+04,
     A .8917E+03,.9006E+03,.1040E+04,.1118E+04,.1170E+04,.1094E+04,
     A .9832E+03,.9500E+03,.8107E+03,.7553E+03,.7099E+03,.7316E+03,
     A .8256E+03,.9218E+03,.9631E+03,.1002E+04,.1021E+04,.1050E+04,
     A .1088E+04,.1111E+04,.1121E+04,.1138E+04,.1143E+04,.1126E+04,
     A .1223E+04,.1291E+04,.1248E+04,.1077E+04,.1083E+04,.1060E+04,
     A .1073E+04,.1080E+04,.1029E+04,.1007E+04,.1007E+04,.1007E+04,
     A .1047E+04,.1024E+04,.9110E+03,.8848E+03,.9219E+03,.9529E+03,
     A .9492E+03,.9338E+03,.9503E+03,.1011E+04,.1090E+04,.1117E+04,
     A .1052E+04,.1057E+04,.8313E+03,.8326E+03,.7518E+03,.6958E+03,
     A .7294E+03,.8314E+03,.8769E+03,.1007E+04,.1006E+04,.9759E+03,
     A .9818E+03,.1022E+04,.1126E+04,.1183E+04,.1205E+04,.1192E+04,
     A .1133E+04,.1102E+04,.1097E+04,.9892E+03,.9002E+03,.8646E+03,
     A .8604E+03,.9113E+03,.1048E+04,.1093E+04,.1088E+04,.1038E+04,
     A .9117E+03,.9311E+03,.9963E+03,.1003E+04,.1121E+04,.1226E+04,
     A .1273E+04,.1244E+04,.1213E+04,.1107E+04,.1043E+04,.9732E+03,
     A .9228E+03,.8196E+03,.8918E+03,.1011E+04,.1168E+04,.1187E+04,
     A .1132E+04,.1067E+04,.9934E+03,.9234E+03,.9294E+03,.9616E+03,
     A .9753E+03,.8413E+03,.7828E+03,.6104E+03,.4720E+03,.5256E+03,
     A .5850E+03,.6583E+03,.6704E+03,.6415E+03,.6686E+03,.7419E+03,
     A .7797E+03,.8460E+03,.8015E+03,.7803E+03,.8138E+03,.9126E+03,
     A .9742E+03,.9839E+03,.9962E+03,.1010E+04,.1091E+04,.1169E+04,
     A .1206E+04,.1096E+04,.1006E+04,.8781E+03,.8375E+03,.8681E+03,
     A .1036E+04,.1168E+04,.1197E+04,.1145E+04,.1050E+04,.1036E+04,
     A .1010E+04,.1003E+04,.9028E+03,.8341E+03,.8316E+03,.8443E+03,
     A .8153E+03,.7868E+03,.7700E+03,.7461E+03,.7861E+03,.8809E+03,
     A .1013E+04,.1041E+04,.1120E+04,.1140E+04,.1075E+04,.1005E+04,
     A .8238E+03,.8591E+03,.9326E+03,.1055E+04,.1113E+04,.1122E+04,
     A .1042E+04,.1039E+04,.1131E+04,.1196E+04,.1203E+04,.1161E+04,
     A .9649E+03,.9205E+03,.9317E+03,.9283E+03,.8868E+03,.8197E+03,
     A .7285E+03,.6884E+03,.6650E+03,.6830E+03,.7510E+03,.8687E+03,
     A .9534E+03,.1038E+04,.1066E+04,.1050E+04,.1046E+04,.1028E+04,
     A .1023E+04,.1038E+04,.1008E+04,.9482E+03,.8673E+03,.8508E+03,
     A .8652E+03,.9006E+03,.8945E+03,.9018E+03,.9460E+03,.1099E+04,
     A .1195E+04,.1213E+04,.1103E+04,.1007E+04,.8244E+03,.8003E+03,
     A .7551E+03,.7137E+03,.8081E+03,.9349E+03,.9077E+03,.8166E+03,
     A .7817E+03,.9307E+03,.1018E+04,.1047E+04,.9362E+03,.8742E+03,
     A .8240E+03,.7578E+03,.7125E+03,.7520E+03,.8566E+03,.9949E+03,
     A .1196E+04,.1216E+04,.1194E+04,.1166E+04,.1146E+04,.1147E+04,
     A .1132E+04,.1169E+04,.1278E+04,.1323E+04,.1294E+04,.1118E+04,
     A .1066E+04,.1028E+04,.1026E+04,.1026E+04,.1082E+04,.1149E+04,
     A .1190E+04,.1183E+04,.1133E+04,.1012E+04,.8714E+03,.9166E+03,
     A .9888E+03,.1028E+04,.1007E+04,.9780E+03,.8809E+03,.8249E+03,
     A .7947E+03,.7765E+03,.8229E+03,.9806E+03,.1082E+04,.1146E+04,
     A .1200E+04,.1277E+04,.1254E+04,.1223E+04,.1184E+04,.1102E+04,
     A .1037E+04,.9724E+03,.9266E+03,.9205E+03,.9747E+03,.9708E+03,
     A .9339E+03,.7991E+03,.7295E+03,.7199E+03,.7212E+03,.7020E+03,
     A .6495E+03,.6804E+03,.8295E+03,.9953E+03,.1091E+04,.1169E+04,
     A .1148E+04,.1111E+04,.1091E+04,.1131E+04,.1115E+04,.1075E+04,
     A .9344E+03,.8783E+03,.9406E+03,.1037E+04,.1184E+04,.1221E+04,
     A .1185E+04,.1148E+04,.1102E+04,.1080E+04,.1064E+04,.1054E+04,
     A .1092E+04,.1234E+04,.1307E+04,.1319E+04,.1205E+04,.1072E+04,
     A .1005E+04,.1005E+04,.1043E+04,.1100E+04,.1155E+04,.1208E+04,
     A .1227E+04,.1219E+04,.1178E+04,.1163E+04,.1188E+04,.1200E+04,
     A .1140E+04,.1128E+04,.1075E+04,.1142E+04,.1245E+04,.1241E+04,
     A .1235E+04,.1195E+04,.1161E+04,.1150E+04,.1124E+04,.1070E+04,
     A .1012E+04,.9741E+03,.9220E+03,.9011E+03,.9871E+03,.1146E+04,
     A .1136E+04,.1114E+04,.1103E+04,.1082E+04,.1012E+04,.9636E+03,
     A .1011E+04,.1038E+04,.1060E+04,.1069E+04,.1062E+04,.1104E+04,
     A .1172E+04,.1250E+04,.1292E+04,.1278E+04,.1167E+04,.9438E+03,
     A .6831E+03,.5447E+03,.5923E+03,.7113E+03,.8794E+03,.1008E+04,
     A .9964E+03,.8966E+03,.6964E+03,.5406E+03,.5312E+03,.6215E+03,
     A .7429E+03,.8141E+03,.8502E+03,.8959E+03,.9533E+03,.9358E+03,
     A .9399E+03,.9278E+03,.9520E+03,.9788E+03,.1009E+04,.1048E+04,
     A .1066E+04,.9974E+03,.8694E+03,.8249E+03,.8119E+03,.8154E+03,
     A .7711E+03,.6182E+03,.4364E+03,.4451E+03,.5049E+03,.6126E+03,
     A .7555E+03,.8065E+03,.7970E+03,.7189E+03,.5693E+03,.4925E+03,
     A .5422E+03,.5608E+03,.5946E+03,.6666E+03,.7374E+03,.7900E+03,
     A .7822E+03,.7897E+03,.8388E+03,.9248E+03,.1097E+04,.1193E+04,
     A .1169E+04,.1052E+04,.9886E+03,.9445E+03,.1006E+04,.1049E+04,
     A .1098E+04,.1178E+04,.1225E+04,.1265E+04,.1259E+04,.1198E+04,
     A .1179E+04,.1149E+04,.1167E+04,.1207E+04,.1315E+04,.1346E+04,
     A .1351E+04,.1184E+04,.1029E+04,.9685E+03,.9340E+03,.9442E+03,
     A .9493E+03,.9450E+03,.8985E+03,.8619E+03,.9223E+03,.9927E+03,
     A .1041E+04,.1023E+04,.8902E+03,.7030E+03,.6237E+03,.6570E+03,
     A .7418E+03,.8072E+03,.9245E+03,.1038E+04,.1078E+04,.1020E+04,
     A .9741E+03,.8957E+03,.9290E+03,.1020E+04,.1081E+04,.1160E+04,
     A .1171E+04,.1123E+04,.1069E+04,.8842E+03,.6725E+03,.5089E+03,
     A .4872E+03,.6021E+03,.7666E+03,.8928E+03,.8931E+03,.9406E+03,
     A .1073E+04,.1157E+04,.1116E+04,.1085E+04,.9982E+03,.1029E+04,
     A .1174E+04,.1196E+04,.1342E+04,.1366E+04,.1389E+04,.1341E+04,
     A .1308E+04,.1307E+04,.1314E+04,.1264E+04,.1116E+04,.7879E+03,
     A .6907E+03,.5603E+03,.6743E+03,.8468E+03,.9622E+03,.9752E+03,
     A .9486E+03,.8859E+03,.9639E+03,.1016E+04,.1059E+04,.1059E+04,
     A .1050E+04,.1076E+04,.1110E+04,.1189E+04,.1251E+04,.1270E+04,
     A .1282E+04,.1267E+04,.1265E+04,.1209E+04,.1187E+04,.1107E+04,
     A .1097E+04,.1067E+04,.1033E+04,.9507E+03,.9430E+03,.9473E+03,
     A .1033E+04,.1053E+04,.8858E+03,.8575E+03,.7471E+03,.8767E+03,
     A .9570E+03,.9023E+03,.9218E+03,.1044E+04,.1102E+04,.1106E+04,
     A .1157E+04,.1299E+04,.1355E+04,.1285E+04,.1289E+04,.1278E+04,
     A .1263E+04,.1236E+04,.1269E+04,.1340E+04,.1352E+04,.1325E+04,
     A .1332E+04,.1405E+04,.1437E+04,.1416E+04,.1373E+04,.1365E+04,
     A .1395E+04,.1396E+04,.1385E+04,.1312E+04,.1252E+04,.1196E+04,
     A .1179E+04,.1219E+04,.1318E+04,.1350E+04,.1290E+04,.1269E+04,
     A .1254E+04,.1287E+04,.1313E+04,.1240E+04,.1123E+04,.1072E+04,
     A .1070E+04,.1153E+04,.1242E+04,.1391E+04,.1437E+04,.1434E+04,
     A .1357E+04,.1275E+04,.1274E+04,.1286E+04,.1350E+04,.1384E+04,
     A .1378E+04,.1226E+04,.1179E+04,.1045E+04,.1075E+04,.1137E+04,
     A .1126E+04,.1060E+04,.9821E+03,.1009E+04,.1064E+04,.1173E+04,
     A .1187E+04,.1183E+04,.1207E+04,.1287E+04,.1291E+04,.1252E+04,
     A .1171E+04,.1088E+04,.1079E+04,.1058E+04,.9797E+03,.9062E+03,
     A .9680E+03,.1029E+04,.1079E+04,.1184E+04,.1263E+04,.1359E+04,
     A .1448E+04,.1498E+04,.1497E+04,.1422E+04,.1297E+04,.1202E+04,
     A .1175E+04,.1200E+04,.1275E+04,.1290E+04,.1266E+04,.1264E+04,
     A .1263E+04,.1291E+04,.1321E+04,.1365E+04,.1439E+04,.1504E+04,
     A .1509E+04,.1455E+04,.1415E+04,.1375E+04,.1382E+04,.1364E+04,
     A .1339E+04,.1210E+04,.1060E+04,.9972E+03,.9610E+03,.8558E+03,
     A .8618E+03,.9548E+03,.9995E+03,.9913E+03,.9767E+03,.9705E+03,
     A .9276E+03,.9136E+03,.1082E+04,.1208E+04,.1279E+04,.1299E+04,
     A .1338E+04,.1388E+04,.1483E+04,.1512E+04,.1592E+04,.1596E+04,
     A .1519E+04,.1441E+04,.1359E+04,.1394E+04,.1429E+04,.1427E+04,
     A .1320E+04,.1152E+04,.9414E+03,.8051E+03,.6224E+03,.6505E+03,
     A .7033E+03,.7125E+03,.7098E+03,.7966E+03,.8969E+03,.1104E+04,
     A .1153E+04,.1268E+04,.1302E+04,.1416E+04,.1418E+04,.1221E+04,
     A .1045E+04,.9016E+03,.9799E+03,.1184E+04,.1338E+04,.1399E+04,
     A .1315E+04,.1203E+04,.1135E+04,.1142E+04,.1118E+04,.9892E+03,
     A .8524E+03,.7399E+03,.6875E+03,.5583E+03,.5063E+03,.5843E+03,
     A .6212E+03,.6261E+03,.6352E+03,.6809E+03,.7743E+03,.1001E+04,
     A .1128E+04,.1194E+04,.1224E+04,.1360E+04,.1348E+04,.1343E+04,
     A .1301E+04,.1112E+04,.9566E+03,.9186E+03,.9318E+03,.9440E+03,
     A .7892E+03,.7076E+03,.6966E+03,.7813E+03,.9161E+03,.8236E+03,
     A .7522E+03,.6064E+03,.5095E+03,.5481E+03,.6151E+03,.8432E+03,
     A .1053E+04,.1184E+04,.1249E+04,.1254E+04,.1252E+04,.1289E+04,
     A .1374E+04,.1445E+04,.1539E+04,.1559E+04,.1522E+04,.1446E+04,
     A .1212E+04,.9425E+03,.7658E+03,.7354E+03,.8211E+03,.1008E+04,
     A .1135E+04,.1178E+04,.1192E+04,.1181E+04,.1221E+04,.1197E+04,
     A .1064E+04,.9627E+03,.1003E+04,.1121E+04,.1224E+04,.1254E+04,
     A .1158E+04,.1004E+04,.9442E+03,.9525E+03,.1159E+04,.1161E+04,
     A .1249E+04,.1242E+04,.1080E+04,.1060E+04,.1084E+04,.1181E+04,
     A .1325E+04,.1382E+04,.1422E+04,.1431E+04,.1404E+04,.1384E+04,
     A .1369E+04,.1350E+04,.1358E+04,.1402E+04,.1430E+04,.1415E+04,
     A .1368E+04,.1370E+04,.1390E+04,.1478E+04,.1565E+04,.1588E+04,
     A .1556E+04,.1529E+04,.1501E+04,.1504E+04,.1527E+04,.1556E+04,
     A .1633E+04,.1672E+04,.1569E+04,.1355E+04,.1152E+04,.1120E+04,
     A .1151E+04,.1168E+04,.1161E+04,.1150E+04,.1093E+04,.1068E+04,
     A .1066E+04,.1145E+04,.1217E+04,.1316E+04,.1296E+04,.1237E+04,
     A .1182E+04,.1116E+04,.1051E+04,.9967E+03,.9802E+03,.1006E+04,
     A .1099E+04,.1106E+04,.1042E+04,.8935E+03,.8225E+03,.8103E+03,
     A .8386E+03,.9439E+03,.1097E+04,.1239E+04,.1264E+04,.1238E+04,
     A .1271E+04,.1340E+04,.1404E+04,.1459E+04,.1454E+04,.1422E+04,
     A .1400E+04,.1386E+04,.1290E+04,.1173E+04,.1098E+04,.1091E+04,
     A .1225E+04,.1334E+04,.1379E+04,.1410E+04,.1463E+04,.1473E+04,
     A .1453E+04,.1427E+04,.1369E+04,.1236E+04,.1100E+04,.1057E+04,
     A .1034E+04,.1007E+04,.9270E+03,.8178E+03,.7721E+03,.9013E+03,
     A .1046E+04,.1166E+04,.1226E+04,.1203E+04,.1104E+04,.9165E+03,
     A .7222E+03,.6028E+03,.6667E+03,.7821E+03,.9722E+03,.1156E+04,
     A .1119E+04,.1087E+04,.9148E+03,.7431E+03,.6894E+03,.6283E+03,
     A .6077E+03,.6586E+03,.7114E+03,.7116E+03,.7018E+03,.7242E+03,
     A .7378E+03,.7429E+03,.8136E+03,.9002E+03,.8662E+03,.7804E+03,
     A .6987E+03,.5854E+03,.6082E+03,.6162E+03,.6120E+03,.6473E+03,
     A .7198E+03,.8205E+03,.8445E+03,.8589E+03,.8301E+03,.6959E+03,
     A .6322E+03,.5805E+03,.6415E+03,.6773E+03,.6865E+03,.6290E+03,
     A .6744E+03,.8807E+03,.1000E+04,.1154E+04,.1151E+04,.1139E+04,
     A .1173E+04,.1203E+04,.1192E+04,.1165E+04,.1142E+04,.1134E+04,
     A .1198E+04,.1268E+04,.1316E+04,.1302E+04,.1084E+04,.1013E+04,
     A .9354E+03,.9204E+03,.1045E+04,.1234E+04,.1306E+04,.1280E+04,
     A .1205E+04,.1127E+04,.1091E+04,.1072E+04,.1072E+04,.1014E+04,
     A .9254E+03,.9235E+03,.1033E+04,.1066E+04,.9938E+03,.8652E+03,
     A .7097E+03,.6158E+03,.5987E+03,.6297E+03,.7898E+03,.9760E+03,
     A .1119E+04,.1300E+04,.1349E+04,.1327E+04,.1284E+04,.1273E+04,
     A .1197E+04,.1121E+04,.1157E+04,.1174E+04,.1183E+04,.1150E+04,
     A .1136E+04,.1135E+04,.1115E+04,.1044E+04,.9334E+03,.8875E+03,
     A .7909E+03,.7574E+03,.7667E+03,.7839E+03,.8638E+03,.1021E+04,
     A .1153E+04,.1271E+04,.1356E+04,.1369E+04,.1319E+04,.1317E+04,
     A .1203E+04,.1028E+04,.8211E+03,.8321E+03,.8714E+03,.1057E+04,
     A .1128E+04,.1128E+04,.1075E+04,.9737E+03,.8498E+03,.8241E+03,
     A .8891E+03,.1099E+04,.1291E+04,.1294E+04,.1194E+04,.1011E+04,
     A .8193E+03,.7686E+03,.8519E+03,.9850E+03,.9389E+03,.8978E+03,
     A .8636E+03,.9215E+03,.1049E+04,.1170E+04,.1231E+04,.1339E+04,
     A .1373E+04,.1349E+04,.1325E+04,.1307E+04,.1323E+04,.1282E+04,
     A .1249E+04,.1174E+04,.1229E+04,.1417E+04,.1482E+04,.1457E+04,
     A .1361E+04,.1340E+04,.1295E+04,.1228E+04,.1219E+04,.1309E+04,
     A .1403E+04,.1449E+04,.1345E+04,.1166E+04,.1126E+04,.1152E+04,
     A .1228E+04,.1200E+04,.1148E+04,.1006E+04,.9786E+03,.1164E+04,
     A .1207E+04,.1352E+04,.1416E+04,.1452E+04,.1470E+04,.1437E+04,
     A .1373E+04,.1373E+04,.1401E+04,.1475E+04,.1563E+04,.1535E+04,
     A .1536E+04,.1450E+04,.1408E+04,.1457E+04,.1474E+04,.1456E+04,
     A .1377E+04,.1354E+04,.1378E+04,.1404E+04,.1392E+04,.1303E+04,
     A .1289E+04,.1258E+04,.1165E+04,.1022E+04,.1004E+04,.1050E+04,
     A .1121E+04,.1120E+04,.1120E+04,.1163E+04,.1276E+04,.1271E+04,
     A .1217E+04,.1155E+04,.1062E+04,.1037E+04,.9766E+03,.9336E+03,
     A .7666E+03,.6929E+03,.6690E+03,.6394E+03,.5889E+03,.5242E+03,
     A .4950E+03,.4691E+03,.4256E+03,.3553E+03,.3031E+03,.2631E+03,
     A .2516E+03,.2938E+03,.3419E+03,.3898E+03,.4409E+03,.5326E+03,
     A .6067E+03,.6318E+03,.6838E+03,.7927E+03,.8823E+03,.9404E+03,
     A .1002E+04,.1014E+04,.1014E+04,.1012E+04,.1054E+04,.1078E+04,
     A .1025E+04,.9437E+03,.8893E+03,.8983E+03,.1059E+04,.1178E+04,
     A .1337E+04,.1412E+04,.1356E+04,.1283E+04,.1236E+04,.1225E+04,
     A .1294E+04,.1396E+04,.1450E+04,.1510E+04,.1532E+04,.1482E+04,
     A .1416E+04,.1238E+04,.1238E+04,.1333E+04,.1467E+04,.1552E+04,
     A .1510E+04,.1427E+04,.1327E+04,.1242E+04,.1293E+04,.1347E+04,
     A .1394E+04,.1449E+04,.1517E+04,.1509E+04,.1458E+04,.1245E+04,
     A .9816E+03,.8695E+03,.8172E+03,.8490E+03,.8989E+03,.9140E+03,
     A .9110E+03,.8610E+03,.8048E+03,.7170E+03,.6353E+03,.5451E+03,
     A .4482E+03,.3753E+03,.3407E+03,.2967E+03,.2898E+03,.3393E+03,
     A .4312E+03,.4787E+03,.6118E+03,.6815E+03,.7431E+03,.8003E+03,
     A .8591E+03,.9084E+03,.9458E+03,.9973E+03,.1134E+04,.1199E+04,
     A .1207E+04,.1229E+04,.1270E+04,.1354E+04,.1397E+04,.1447E+04,
     A .1517E+04,.1582E+04,.1616E+04,.1640E+04,.1624E+04,.1542E+04,
     A .1521E+04,.1512E+04,.1509E+04,.1471E+04,.1462E+04,.1508E+04,
     A .1568E+04,.1583E+04,.1547E+04,.1483E+04,.1452E+04,.1658E+04,
     A .1739E+04,.1763E+04,.1699E+04,.1630E+04,.1608E+04,.1626E+04,
     A .1693E+04,.1732E+04,.1740E+04,.1760E+04,.1782E+04,.1737E+04,
     A .1699E+04,.1674E+04,.1659E+04,.1659E+04,.1648E+04,.1613E+04,
     A .1552E+04,.1476E+04,.1507E+04,.1697E+04,.1718E+04,.1852E+04,
     A .1859E+04,.1787E+04,.1740E+04,.1761E+04,.1791E+04,.1876E+04,
     A .1922E+04,.1898E+04,.1812E+04,.1639E+04,.1371E+04,.1324E+04,
     A .1403E+04,.1483E+04,.1613E+04,.1713E+04,.1835E+04,.1824E+04,
     A .1787E+04,.1758E+04,.1771E+04,.1786E+04,.1816E+04,.1818E+04,
     A .1815E+04,.1840E+04,.1850E+04,.1808E+04,.1766E+04,.1786E+04,
     A .1840E+04,.1885E+04,.1841E+04,.1787E+04,.1680E+04,.1622E+04,
     A .1706E+04,.1816E+04,.1904E+04,.1913E+04,.1825E+04,.1823E+04,
     A .1793E+04,.1803E+04,.1809E+04,.1851E+04,.1859E+04,.1841E+04,
     A .1773E+04,.1707E+04,.1708E+04,.1769E+04,.1844E+04,.1893E+04,
     A .1905E+04,.1934E+04,.1951E+04,.1959E+04,.1947E+04,.1776E+04,
     A .1584E+04,.1443E+04,.1458E+04,.1548E+04,.1553E+04,.1486E+04,
     A .1484E+04,.1514E+04,.1575E+04,.1600E+04,.1648E+04,.1706E+04,
     A .1817E+04,.1926E+04,.1989E+04,.2023E+04,.2031E+04,.2013E+04,
     A .1982E+04,.1949E+04,.1836E+04,.1761E+04,.1771E+04,.1811E+04,
     A .1976E+04,.2058E+04,.2000E+04,.1863E+04,.1633E+04,.1384E+04,
     A .1164E+04,.1056E+04,.1043E+04,.1413E+04,.1679E+04,.1816E+04,
     A .1909E+04,.1865E+04,.1851E+04,.1881E+04,.1964E+04,.2016E+04,
     A .1990E+04,.1928E+04,.1804E+04,.1750E+04,.1779E+04,.1793E+04,
     A .1750E+04,.1676E+04,.1660E+04,.1691E+04,.1759E+04,.1739E+04,
     A .1616E+04,.1541E+04,.1519E+04,.1563E+04,.1656E+04,.1727E+04,
     A .1856E+04,.1922E+04,.1905E+04,.1862E+04,.1726E+04,.1492E+04,
     A .1253E+04,.1166E+04,.1243E+04,.1529E+04,.1701E+04,.1756E+04,
     A .1695E+04,.1594E+04,.1554E+04,.1585E+04,.1745E+04,.1884E+04,
     A .1961E+04,.1969E+04,.1885E+04,.1651E+04,.1365E+04,.1286E+04,
     A .1314E+04,.1629E+04,.1828E+04,.1875E+04,.1852E+04,.1789E+04,
     A .1740E+04,.1710E+04,.1591E+04,.1444E+04,.1370E+04,.1375E+04,
     A .1452E+04,.1568E+04/
      
        n=natlas
	x1(1:n)=x1_atlas(1:n)
	y1(1:n)=y1_atlas(1:n)
      end 	

      Subroutine read_neck(x1,y1,kdata,n) ! data neckel extra-terrestrial flux data
      parameter (neckel=496)
      real x1(kdata),y1(kdata),x1_neck(neckel),y1_neck(neckel)
       data x1_neck/ 330.0, 331.0, 332.0, 333.0, 334.0,
     A  335.0, 336.0, 337.0, 338.0, 339.0, 340.0, 341.0, 342.0, 343.0,
     A  344.0, 345.0, 346.0, 347.0, 348.0, 349.0, 350.0, 351.0, 352.0,
     A  353.0, 354.0, 355.0, 356.0, 357.0, 358.0, 359.0, 360.0, 361.0,
     A  362.0, 363.0, 364.0, 365.0, 366.0, 367.0, 368.0, 369.0, 370.0,
     A  371.0, 372.0, 373.0, 374.0, 375.0, 376.0, 377.0, 378.0, 379.0,
     A  380.0, 381.0, 382.0, 383.0, 384.0, 385.0, 386.0, 387.0, 388.0,
     A  389.0, 390.0, 391.0, 392.0, 393.0, 394.0, 395.0, 396.0, 397.0,
     A  398.0, 399.0, 400.0, 401.0, 402.0, 403.0, 404.0, 405.0, 406.0,
     A  407.0, 408.0, 409.0, 410.0, 411.0, 412.0, 413.0, 414.0, 415.0,
     A  416.0, 417.0, 418.0, 419.0, 420.0, 421.0, 422.0, 423.0, 424.0,
     A  425.0, 426.0, 427.0, 428.0, 429.0, 430.0, 431.0, 432.0, 433.0,
     A  434.0, 435.0, 436.0, 437.0, 438.0, 439.0, 440.0, 441.0, 442.0,
     A  443.0, 444.0, 445.0, 446.0, 447.0, 448.0, 449.0, 450.0, 451.0,
     A  452.0, 453.0, 454.0, 455.0, 456.0, 457.0, 458.0, 459.0, 460.0,
     A  461.0, 462.0, 463.0, 464.0, 465.0, 466.0, 467.0, 468.0, 469.0,
     A  470.0, 471.0, 472.0, 473.0, 474.0, 475.0, 476.0, 477.0, 478.0,
     A  479.0, 480.0, 481.0, 482.0, 483.0, 484.0, 485.0, 486.0, 487.0,
     A  488.0, 489.0, 490.0, 491.0, 492.0, 493.0, 494.0, 495.0, 496.0,
     A  497.0, 498.0, 499.0, 500.0, 501.0, 502.0, 503.0, 504.0, 505.0,
     A  506.0, 507.0, 508.0, 509.0, 510.0, 511.0, 512.0, 513.0, 514.0,
     A  515.0, 516.0, 517.0, 518.0, 519.0, 520.0, 521.0, 522.0, 523.0,
     A  524.0, 525.0, 526.0, 527.0, 528.0, 529.0, 530.0, 531.0, 532.0,
     A  533.0, 534.0, 535.0, 536.0, 537.0, 538.0, 539.0, 540.0, 541.0,
     A  542.0, 543.0, 544.0, 545.0, 546.0, 547.0, 548.0, 549.0, 550.0,
     A  551.0, 552.0, 553.0, 554.0, 555.0, 556.0, 557.0, 558.0, 559.0,
     A  560.0, 561.0, 562.0, 563.0, 564.0, 565.0, 566.0, 567.0, 568.0,
     A  569.0, 570.0, 571.0, 572.0, 573.0, 574.0, 575.0, 576.0, 577.0,
     A  578.0, 579.0, 580.0, 581.0, 582.0, 583.0, 584.0, 585.0, 586.0,
     A  587.0, 588.0, 589.0, 590.0, 591.0, 592.0, 593.0, 594.0, 595.0,
     A  596.0, 597.0, 598.0, 599.0, 600.0, 601.0, 602.0, 603.0, 604.0,
     A  605.0, 606.0, 607.0, 608.0, 609.0, 610.0, 611.0, 612.0, 613.0,
     A  614.0, 615.0, 616.0, 617.0, 618.0, 619.0, 620.0, 621.0, 622.0,
     A  623.0, 624.0, 625.0, 626.0, 627.0, 628.0, 629.0, 630.0, 632.0,
     A  634.0, 636.0, 638.0, 640.0, 642.0, 644.0, 646.0, 648.0, 650.0,
     A  652.0, 654.0, 656.0, 658.0, 660.0, 662.0, 664.0, 666.0, 668.0,
     A  670.0, 672.0, 674.0, 676.0, 678.0, 680.0, 682.0, 684.0, 686.0,
     A  688.0, 690.0, 692.0, 694.0, 696.0, 698.0, 700.0, 702.0, 704.0,
     A  706.0, 708.0, 710.0, 712.0, 714.0, 716.0, 718.0, 720.0, 722.0,
     A  724.0, 726.0, 728.0, 730.0, 732.0, 734.0, 736.0, 738.0, 740.0,
     A  742.0, 744.0, 746.0, 748.0, 750.0, 752.0, 754.0, 756.0, 758.0,
     A  760.0, 762.0, 764.0, 766.0, 768.0, 770.0, 772.0, 774.0, 776.0,
     A  778.0, 780.0, 782.0, 784.0, 786.0, 788.0, 790.0, 792.0, 794.0,
     A  796.0, 798.0, 800.0, 802.0, 804.0, 806.0, 808.0, 810.0, 812.0,
     A  814.0, 816.0, 818.0, 820.0, 822.0, 824.0, 826.0, 828.0, 830.0,
     A  832.0, 834.0, 836.0, 838.0, 840.0, 842.0, 844.0, 846.0, 848.0,
     A  850.0, 852.0, 854.0, 856.0, 858.0, 860.0, 862.0, 864.0, 866.0,
     A  868.0, 870.0, 875.0, 880.0, 885.0, 890.0, 895.0, 900.0, 905.0,
     A  910.0, 915.0, 920.0, 925.0, 930.0, 935.0, 940.0, 945.0, 950.0,
     A  955.0, 960.0, 965.0, 970.0, 975.0, 980.0, 985.0, 990.0, 995.0,
     A 1000.0,1005.0,1010.0,1015.0,1020.0,1025.0,1030.0,1035.0,1040.0,
     A 1045.0,1050.0,1055.0,1060.0,1065.0,1070.0,1075.0,1080.0,1085.0,
     A 1090.0,1095.0,1100.0,1105.0,1110.0,1115.0,1120.0,1125.0,1130.0,
     A 1135.0,1140.0,1145.0,1150.0,1155.0,1160.0,1165.0,1170.0,1175.0,
     A 1180.0,1185.0,1190.0,1195.0,1200.0,1205.0,1210.0,1215.0,1220.0,
     A 1225.0,1230.0,1235.0,1240.0,1245.0/
     
       data y1_neck/.1006E+01,.9682E+00,.9213E+00,.9053E+00,
     A .9402E+00,.9822E+00,.7654E+00,.8663E+00,.9163E+00,.9372E+00,
     A .9922E+00,.9362E+00,.9952E+00,.9852E+00,.7194E+00,.9672E+00,
     A .9193E+00,.9023E+00,.9482E+00,.8653E+00,.1120E+01,.9932E+00,
     A .8713E+00,.1116E+01,.1134E+01,.1058E+01,.9382E+00,.8913E+00,
     A .6275E+00,.1137E+01,.9792E+00,.8943E+00,.1176E+01,.9582E+00,
     A .1015E+01,.1264E+01,.1250E+01,.1215E+01,.1089E+01,.1332E+01,
     A .1076E+01,.1308E+01,.1065E+01,.8383E+00,.8783E+00,.1142E+01,
     A .1102E+01,.1292E+01,.1342E+01,.1000E+01,.1290E+01,.1097E+01,
     A .7334E+00,.6844E+00,.1027E+01,.9542E+00,.1071E+01,.9662E+00,
     A .9123E+00,.1228E+01,.1224E+01,.1399E+01,.9552E+00,.4896E+00,
     A .1102E+01,.1379E+01,.6505E+00,.1040E+01,.1539E+01,.1656E+01,
     A .1650E+01,.1798E+01,.1805E+01,.1659E+01,.1603E+01,.1673E+01,
     A .1625E+01,.1546E+01,.1826E+01,.1707E+01,.1503E+01,.1821E+01,
     A .1793E+01,.1759E+01,.1740E+01,.1737E+01,.1846E+01,.1668E+01,
     A .1687E+01,.1704E+01,.1761E+01,.1801E+01,.1585E+01,.1714E+01,
     A .1771E+01,.1698E+01,.1701E+01,.1572E+01,.1590E+01,.1478E+01,
     A .1137E+01,.1689E+01,.1649E+01,.1734E+01,.1673E+01,.1726E+01,
     A .1932E+01,.1810E+01,.1570E+01,.1829E+01,.1716E+01,.1934E+01,
     A .1983E+01,.1912E+01,.1976E+01,.1825E+01,.1894E+01,.2080E+01,
     A .1976E+01,.2030E+01,.2147E+01,.2112E+01,.1944E+01,.1973E+01,
     A .1982E+01,.2037E+01,.2080E+01,.2103E+01,.1974E+01,.2012E+01,
     A .2043E+01,.2058E+01,.2107E+01,.2043E+01,.1979E+01,.2045E+01,
     A .1924E+01,.2018E+01,.1997E+01,.1993E+01,.1880E+01,.2021E+01,
     A .2044E+01,.1994E+01,.2054E+01,.2019E+01,.1959E+01,.2078E+01,
     A .2012E+01,.2079E+01,.2038E+01,.2093E+01,.2026E+01,.2022E+01,
     A .1972E+01,.1834E+01,.1628E+01,.1834E+01,.1917E+01,.1963E+01,
     A .2010E+01,.1899E+01,.1899E+01,.1891E+01,.2061E+01,.1929E+01,
     A .2020E+01,.2021E+01,.1869E+01,.1973E+01,.1860E+01,.1816E+01,
     A .1897E+01,.1937E+01,.1872E+01,.1996E+01,.1964E+01,.1909E+01,
     A .1922E+01,.1919E+01,.1950E+01,.2000E+01,.1870E+01,.1864E+01,
     A .1877E+01,.1903E+01,.1672E+01,.1729E+01,.1657E+01,.1832E+01,
     A .1835E+01,.1909E+01,.1827E+01,.1897E+01,.1961E+01,.1933E+01,
     A .1677E+01,.1832E+01,.1900E+01,.1921E+01,.1955E+01,.1966E+01,
     A .1774E+01,.1926E+01,.1861E+01,.1993E+01,.1874E+01,.1885E+01,
     A .1907E+01,.1836E+01,.1773E+01,.1884E+01,.1829E+01,.1882E+01,
     A .1882E+01,.1904E+01,.1882E+01,.1837E+01,.1866E+01,.1898E+01,
     A .1865E+01,.1874E+01,.1849E+01,.1885E+01,.1901E+01,.1900E+01,
     A .1825E+01,.1849E+01,.1791E+01,.1812E+01,.1846E+01,.1828E+01,
     A .1853E+01,.1864E+01,.1857E+01,.1802E+01,.1833E+01,.1890E+01,
     A .1814E+01,.1863E+01,.1773E+01,.1827E+01,.1895E+01,.1879E+01,
     A .1870E+01,.1834E+01,.1849E+01,.1860E+01,.1787E+01,.1832E+01,
     A .1842E+01,.1856E+01,.1876E+01,.1860E+01,.1863E+01,.1787E+01,
     A .1834E+01,.1851E+01,.1753E+01,.1615E+01,.1817E+01,.1791E+01,
     A .1812E+01,.1800E+01,.1777E+01,.1786E+01,.1809E+01,.1784E+01,
     A .1761E+01,.1778E+01,.1749E+01,.1754E+01,.1722E+01,.1791E+01,
     A .1780E+01,.1767E+01,.1763E+01,.1761E+01,.1746E+01,.1747E+01,
     A .1706E+01,.1749E+01,.1708E+01,.1686E+01,.1716E+01,.1716E+01,
     A .1612E+01,.1710E+01,.1727E+01,.1710E+01,.1737E+01,.1693E+01,
     A .1716E+01,.1669E+01,.1659E+01,.1635E+01,.1700E+01,.1700E+01,
     A .1700E+01,.1680E+01,.1642E+01,.1654E+01,.1659E+01,.1657E+01,
     A .1654E+01,.1617E+01,.1624E+01,.1630E+01,.1606E+01,.1561E+01,
     A .1609E+01,.1602E+01,.1535E+01,.1387E+01,.1552E+01,.1574E+01,
     A .1558E+01,.1563E+01,.1538E+01,.1549E+01,.1519E+01,.1524E+01,
     A .1513E+01,.1511E+01,.1501E+01,.1495E+01,.1482E+01,.1458E+01,
     A .1470E+01,.1464E+01,.1451E+01,.1451E+01,.1439E+01,.1419E+01,
     A .1428E+01,.1389E+01,.1391E+01,.1418E+01,.1403E+01,.1387E+01,
     A .1388E+01,.1376E+01,.1369E+01,.1356E+01,.1330E+01,.1333E+01,
     A .1350E+01,.1352E+01,.1348E+01,.1321E+01,.1328E+01,.1320E+01,
     A .1311E+01,.1309E+01,.1280E+01,.1260E+01,.1288E+01,.1281E+01,
     A .1285E+01,.1272E+01,.1264E+01,.1261E+01,.1257E+01,.1250E+01,
     A .1242E+01,.1239E+01,.1243E+01,.1223E+01,.1187E+01,.1205E+01,
     A .1206E+01,.1210E+01,.1190E+01,.1198E+01,.1189E+01,.1189E+01,
     A .1178E+01,.1182E+01,.1179E+01,.1176E+01,.1160E+01,.1145E+01,
     A .1136E+01,.1154E+01,.1137E+01,.1144E+01,.1131E+01,.1117E+01,
     A .1122E+01,.1097E+01,.1116E+01,.1117E+01,.1109E+01,.1106E+01,
     A .1065E+01,.1082E+01,.1075E+01,.1077E+01,.1077E+01,.1074E+01,
     A .1069E+01,.1034E+01,.1053E+01,.1052E+01,.1042E+01,.1045E+01,
     A .1028E+01,.1033E+01,.1025E+01,.9712E+00,.1003E+01,.9732E+00,
     A .8773E+00,.1011E+01,.9972E+00,.9972E+00,.9992E+00,.9702E+00,
     A .8803E+00,.9672E+00,.9972E+00,.9882E+00,.9782E+00,.9682E+00,
     A .9572E+00,.9472E+00,.9372E+00,.9272E+00,.9173E+00,.9083E+00,
     A .8983E+00,.8883E+00,.8793E+00,.8693E+00,.8603E+00,.8513E+00,
     A .8423E+00,.8333E+00,.8243E+00,.8153E+00,.8063E+00,.7974E+00,
     A .7894E+00,.7804E+00,.7724E+00,.7644E+00,.7554E+00,.7464E+00,
     A .7374E+00,.7294E+00,.7204E+00,.7124E+00,.7044E+00,.6954E+00,
     A .6874E+00,.6794E+00,.6715E+00,.6625E+00,.6545E+00,.6485E+00,
     A .6435E+00,.6375E+00,.6315E+00,.6265E+00,.6205E+00,.6145E+00,
     A .6095E+00,.6045E+00,.5985E+00,.5935E+00,.5885E+00,.5825E+00,
     A .5775E+00,.5725E+00,.5675E+00,.5625E+00,.5575E+00,.5526E+00,
     A .5476E+00,.5426E+00,.5376E+00,.5326E+00,.5286E+00,.5246E+00,
     A .5196E+00,.5156E+00,.5116E+00,.5076E+00,.5036E+00,.4996E+00,
     A .4956E+00,.4916E+00,.4876E+00,.4826E+00,.4776E+00,.4726E+00/
        n=neckel
	x1(1:n)=x1_neck(1:n)
	y1(1:n)=y1_neck(1:n)
      end 	


      SUBROUTINE rn_O3(nw,wl,wc,nz,tlev,airlev,j,sq,jlabel)

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Provide the product of (cross section) x (quantum yield) for the two     =*
*=  O3 photolysis reactions:                                                 =*
*=             (a) O3 + hv -> O2 + O(1D)                                     =*
*=             (b) O3 + hv -> O2 + O(3P)                                     =*
*=  Cross section:  Combined data from WMO 85 Ozone Assessment (use 273K     =*
*=                  value from 175.439-847.5 nm) and data from Molina and    =*
*=                  Molina (use in Hartley and Huggins bans (240.5-350 nm)   =*
*=  Quantum yield:  Choice between                                           =*
*=                   (1) data from Michelsen et al, 1994                     =*
*=                   (2) JPL 87 recommendation                               =*
*=                   (3) JPL 90/92 recommendation (no "tail")                =*
*=                   (4) data from Shetter et al., 1996                      =*
*=                   (5) JPL 97 recommendation                               =*
*-----------------------------------------------------------------------------*
*=  PARAMETERS:                                                              =*
*=  NW     - INTEGER, number of specified intervals + 1 in working        (I)=*
*=           wavelength grid                                                 =*
*=  WL     - REAL, vector of lower limits of wavelength intervals in      (I)=*
*=           working wavelength grid                                         =*
*=  WC     - REAL, vector of center points of wavelength intervals in     (I)=*
*=           working wavelength grid                                         =*
*=  NZ     - INTEGER, number of altitude levels in working altitude grid  (I)=*
*=  TLEV   - REAL, temperature (K) at each specified altitude level       (I)=*
*=  AIRLEV - REAL, air density (molec/cc) at each altitude level          (I)=*
*=  J      - INTEGER, counter for number of weighting functions defined  (IO)=*
*=  SQ     - REAL, cross section x quantum yield (cm^2) for each          (O)=*
*=           photolysis reaction defined, at each defined wavelength and     =*
*=           at each defined altitude level                                  =*
*=  JLABEL - CHARACTER*40, string identifier for each photolysis reaction (O)=*
*=           defined                                                         =*
*-----------------------------------------------------------------------------*
*=  EDIT HISTORY:                                                            =*
*=  05/98  Original, adapted from former JSPEC1 subroutine                   =*
*-----------------------------------------------------------------------------*
*= This program is free software;  you can redistribute it and/or modify     =*
*= it under the terms of the GNU General Public License as published by the  =*
*= Free Software Foundation;  either version 2 of the license, or (at your   =*
*= option) any later version.                                                =*
*= The TUV package is distributed in the hope that it will be useful, but    =*
*= WITHOUT ANY WARRANTY;  without even the implied warranty of MERCHANTIBI-  =*
*= LITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public     =*
*= License for more details.                                                 =*
*= To obtain a copy of the GNU General Public License, write to:             =*
*= Free Software Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.   =*
*-----------------------------------------------------------------------------*

      IMPLICIT NONE
      INCLUDE 'tuv.params'

* input

      INTEGER nw
      REAL wl(kw), wc(kw)

      INTEGER nz

      REAL tlev(kz)
      REAL airlev(kz)

* weighting functions

      CHARACTER*40 jlabel(kj)
      REAL sq(kj,kz,kw)

* input/output:
      INTEGER j

* data arrays

      INTEGER n1, n2, n3, n4, n5
      INTEGER kdata,ny,ns,nwmo
      PARAMETER (kdata=250,ny=21,nwmo=158,ns=220)
      REAL xs1(kdata), xs2(kdata), xs3(kdata), xs1wmo(kdata)
      REAL ys1(kdata), ys2(kdata), ys3(kdata), ys1wmo(kdata)
      real xy1(kdata), xy2(kdata), yy1(kdata), yy2(kdata), 
     1 xwork(kdata), ywork(kdata)

       data xs1wmo(1:nwmo)/176.2,177.8,179.4,181.0,182.7,184.3,185.7,
     1 187.4,189.6,191.4,193.2,195.1,197.0,199.0,
     2 201.0,203.1,205.1,207.3,209.4,211.6,213.9,
     3 216.2,218.6,221.0,223.5,226.0,228.6,231.2,
     4 233.9,236.7,239.5,242.4,245.4,248.5,251.6,
     5 254.8,258.1,261.4,264.9,268.5,272.1,275.9,
     6 279.7,283.7,287.8,292.0,296.3,300.8,305.4,
     7 310.1,315.0,320.0,325.0,330.0,335.0,340.0,
     8 345.0,350.0,355.0,360.0,365.0,370.0,375.0,
     9 380.0,385.0,390.0,395.0,400.0,405.0,410.0,
     a 415.0,420.0,425.0,430.0,435.0,440.0,445.0,
     b 450.0,455.0,460.0,465.0,470.0,475.0,480.0,
     c 485.0,490.0,495.0,500.0,505.0,510.0,515.0,
     d 520.0,525.0,530.0,535.0,540.0,545.0,550.0,
     e 555.0,560.0,565.0,570.0,575.0,580.0,585.0,
     f 590.0,595.0,600.0,605.0,610.0,615.0,620.0,
     g 625.0,630.0,635.0,640.0,645.0,650.0,655.0,
     h 660.0,665.0,670.0,675.0,680.0,685.0,690.0,
     i 695.0,700.0,705.0,710.0,715.0,720.0,725.0,
     j 730.0,735.0,740.0,745.0,750.0,755.0,760.0,
     k 765.0,770.0,775.0,780.0,785.0,790.0,795.0,
     l 800.0,805.0,810.0,815.0,820.0,825.0,830.0,
     m 835.0,840.0,845.0,850.0/
       data ys1wmo(1:nwmo)/.811E-18,.799E-18,.786E-18,.763E-18,.729E-18,
     1  .688E-18,.640E-18,.588E-18,.531E-18,.480E-18,.438E-18,.411E-18,
     2 .369E-18,.330E-18,.326E-18,.326E-18,.351E-18,.411E-18,.484E-18,
     3 .626E-18,.857E-18,.117E-17,.152E-17,.197E-17,.255E-17,.324E-17,
     4 .400E-17,.483E-17,.579E-17,.686E-17,.797E-17,.900E-17,.100E-16,
     5 .108E-16,.113E-16,.115E-16,.112E-16,.106E-16,.965E-17,.834E-17,
     6 .692E-17,.542E-17,.402E-17,.277E-17,.179E-17,.109E-17,.624E-18,
     7 .343E-18,.185E-18,.980E-19,.501E-19,.249E-19,.120E-19,.617E-20,
     8 .274E-20,.117E-20,.588E-21,.266E-21,.109E-21,.549E-22,.000E+00,
     9 .000E+00,.000E+00,.000E+00,.000E+00,.000E+00,.000E+00,.000E+00,
     a .000E+00,.291E-22,.314E-22,.399E-22,.654E-22,.683E-22,.866E-22,
     b .125E-21,.149E-21,.171E-21,.212E-21,.357E-21,.368E-21,.406E-21,
     c .489E-21,.711E-21,.843E-21,.828E-21,.909E-21,.122E-20,.162E-20,
     d .158E-20,.160E-20,.178E-20,.207E-20,.255E-20,.274E-20,.288E-20,
     e .307E-20,.317E-20,.336E-20,.388E-20,.431E-20,.467E-20,.475E-20,
     f .455E-20,.435E-20,.442E-20,.461E-20,.489E-20,.484E-20,.454E-20,
     g .424E-20,.390E-20,.360E-20,.343E-20,.317E-20,.274E-20,.261E-20,
     h .242E-20,.220E-20,.202E-20,.185E-20,.167E-20,.154E-20,.142E-20,
     i .125E-20,.112E-20,.102E-20,.920E-21,.840E-21,.770E-21,.690E-21,
     j .630E-21,.570E-21,.525E-21,.475E-21,.447E-21,.420E-21,.375E-21,
     k .325E-21,.292E-21,.276E-21,.270E-21,.280E-21,.285E-21,.252E-21,
     l .220E-21,.182E-21,.163E-21,.175E-21,.190E-21,.185E-21,.170E-21,
     m .152E-21,.142E-21,.140E-21,.140E-21,.142E-21,.145E-21/
       data xs1(1:ns)/240.5,241.0,241.5,242.0,242.5,243.0,243.5,
     1  244.0,244.5,245.0,245.5,246.0,246.5,247.0,
     2  247.5,248.0,248.5,249.0,249.5,250.0,250.5,
     3  251.0,251.5,252.0,252.5,253.0,253.5,254.0,
     4  254.5,255.0,255.5,256.0,256.5,257.0,257.5,
     5  258.0,258.5,259.0,259.5,260.0,260.5,261.0,
     6  261.5,262.0,262.5,263.0,263.5,264.0,264.5,
     7  265.0,265.5,266.0,266.5,267.0,267.5,268.0,
     8  268.5,269.0,269.5,270.0,270.5,271.0,271.5,
     9  272.0,272.5,273.0,273.5,274.0,274.5,275.0,
     a  275.5,276.0,276.5,277.0,277.5,278.0,278.5,
     b  279.0,279.5,280.0,280.5,281.0,281.5,282.0,
     c  282.5,283.0,283.5,284.0,284.5,285.0,285.5,
     d  286.0,286.5,287.0,287.5,288.0,288.5,289.0,
     e  289.5,290.0,290.5,291.0,291.5,292.0,292.5,
     f  293.0,293.5,294.0,294.5,295.0,295.5,296.0,
     g  296.5,297.0,297.5,298.0,298.5,299.0,299.5,
     h  300.0,300.5,301.0,301.5,302.0,302.5,303.0,
     i  303.5,304.0,304.5,305.0,305.5,306.0,306.5,
     j  307.0,307.5,308.0,308.5,309.0,309.5,310.0,
     k  310.5,311.0,311.5,312.0,312.5,313.0,313.5,
     l  314.0,314.5,315.0,315.5,316.0,316.5,317.0,
     m  317.5,318.0,318.5,319.0,319.5,320.0,320.5,
     n  321.0,321.5,322.0,322.5,323.0,323.5,324.0,
     o  324.5,325.0,325.5,326.0,326.5,327.0,327.5,
     p  328.0,328.5,329.0,329.5,330.0,330.5,331.0,
     q  331.5,332.0,332.5,333.0,333.5,334.0,334.5,
     r  335.0,335.5,336.0,336.5,337.0,337.5,338.0,
     s  338.5,339.0,339.5,340.0,340.5,341.0,341.5,
     t  342.0,342.5,343.0,343.5,344.0,344.5,345.0,
     u  345.5,346.0,346.5,347.0,347.5,348.0,348.5,
     v  349.0,349.5,350.0/
       data ys1(1:ns)/ 842.4, 864.0, 891.8, 901.6, 916.3, 939.6, 962.8,
     1 975.2, 977.6,1007.0,1040.0,1042.0,1039.0,1058.0,
     2 1080.0,1079.0,1094.0,1124.0,1121.0,1134.0,1118.0,
     3 1123.0,1164.0,1165.0,1143.0,1149.0,1167.0,1169.0,
     4 1145.0,1174.0,1184.0,1158.0,1132.0,1147.0,1134.0,
     5 1130.0,1172.0,1151.0,1084.0,1086.0,1128.0,1095.0,
     6 1056.0,1064.0,1072.0,1016.0, 994.0,1013.0,1013.0,
     7 961.2, 968.3, 946.7, 920.6, 877.9, 900.5, 872.2,
     8 829.3, 803.8, 820.8, 796.1, 754.4, 736.7, 724.8,
     9 710.7, 684.2, 666.0, 640.1, 601.3, 592.7, 587.1,
     a 562.9, 537.9, 519.6, 504.0, 492.3, 461.7, 445.7,
     b 424.5, 412.9, 398.3, 371.8, 364.2, 333.7, 322.0,
     c 310.0, 299.0, 293.5, 264.4, 252.4, 239.8, 225.4,
     d 219.5, 205.3, 197.8, 184.2, 168.7, 162.0, 153.5,
     e 142.9, 136.5, 129.9, 122.5, 111.8, 104.3,  99.6900,
     f  93.7500, 86.3700, 81.4800, 78.2400, 72.7, 65.9800, 62.0200,
     g  58.6500, 54.6900, 50.0200, 47.1400, 44.4900, 42.3100, 38.0,
     h  36.1600, 33.9600, 31.2800,  28.9500, 27.9600, 26.0, 23.35,
     i  21.9800, 21.5600, 19.8400,  17.7100, 16.3100, 16.0300, 15.49,
     j  13.6700,  12.0,  12.0200,  11.1800,  10.6400,9.4250,8.6370,
     k   8.3140,7.9250,7.8310,6.6970,6.0740,5.6910,6.1690,
     l   5.3340,4.2940,4.1860,4.8370,3.8960,3.3510,3.3510,
     m   3.4200,3.0630,2.3400,1.9990,2.7120,2.8590,1.9630,
     n   1.3680,1.3610,2.1170,1.7930,1.5290,0.8902,0.7852,
     o   0.9488,1.4860,1.2410,0.7276,0.4945,0.6158,0.6491,
     p   1.1580,0.6374,0.3460,0.2666,0.2854,0.5695,0.7130,
     q   0.4584,0.2415,0.1567,0.2187,0.3693,0.3788,0.2117,
     r   0.1274,0.0876,0.0871,0.1837,0.2464,0.2590,0.1311,
     s   0.0695,0.0546,0.0713,0.1549,0.1133,0.0594,0.0352,
     t   0.0285,0.0236,0.0372,0.0441,0.1271,0.0730,0.0390,
     u   0.0269,0.0195,0.0178,0.0283,0.0295,0.0181,0.0096,
     v   0.0083,0.0100,0.0076/
       data ys2(1:ns)/ 837.0, 858.6, 886.0, 895.1, 911.1, 932.0, 957.3,
     1  968.4, 970.7, 994.8,1028.0,1032.0,1035.0,1051.0,
     2  1069.0,1072.0,1082.0,1110.0,1110.0,1121.0,1110.0,
     3  1115.0,1153.0,1153.0,1137.0,1142.0,1155.0,1159.0,
     4  1136.0,1167.0,1176.0,1151.0,1125.0,1139.0,1126.0,
     5  1122.0,1164.0,1143.0,1083.0,1085.0,1120.0,1090.0,
     6  1054.0,1061.0,1067.0,1015.0, 994.7,1010.0,1012.0,
     7  960.5, 965.6, 943.4, 920.0, 878.0, 895.1, 871.1,
     8  830.1, 803.7, 822.0, 798.2, 755.5, 736.5, 726.6,
     9  711.4, 687.6, 675.4, 642.1, 609.9, 597.5, 586.1,
     a  561.0, 539.5, 523.3, 508.9, 490.3, 465.0, 452.1,
     b  429.8, 413.8, 395.2, 374.7, 367.7, 337.6, 323.9,
     c  311.1, 301.1, 283.3, 267.4, 257.7, 244.5, 227.6,
     d  220.5, 209.4, 201.0, 185.8, 172.6, 164.1, 155.2,
     e  147.0, 139.6, 131.9, 124.9, 114.5, 107.4, 102.8,
     f  96.59,  88.31, 84.47, 79.43, 74.92, 68.53, 64.66,
     g  60.1, 56.96,  52.3100,  48.7500,  45.9700,  43.4200,  39.59,
     h  37.4, 35.4,  32.2700,  30.1800,  28.9100,  26.7800,  24.52,
     i  23.0, 22.26,  20.8100,  18.6200,  17.3200,  16.7800,  16.12,
     j  14.4900,  13.0600,  12.6700,  11.7300,  11.4500,9.7150,9.3030,
     k   8.8030,8.4800,8.2430,7.2,6.5750,6.1990,6.4070,
     l   5.6700,4.7760,4.5850,5.0280,4.2260,3.7500,3.7400,
     m   3.6620,3.3220,2.6510,2.3500,2.8800,3.0220,2.2230,
     n   1.6700,1.6150,2.2420,1.9750,1.7600,1.1410,0.9893,
     o   1.1040,1.5930,1.3590,0.9125,0.6579,0.7214,0.7691,
     p   1.2210,0.7731,0.4785,0.3750,0.3770,0.6381,0.7511,
     q   0.5535,0.3283,0.2295,0.2731,0.4033,0.4474,0.2788,
     r   0.1730,0.1347,0.1312,0.2216,0.2686,0.2978,0.1641,
     s   0.1002,0.0836,0.1127,0.1670,0.1434,0.0863,0.0537,
     t   0.0488,0.0414,0.0640,0.0633,0.1336,0.0865,0.0590,
     u   0.0439,0.0384,0.0271,0.0323,0.0429,0.0249,0.0217,
     v   0.0163,0.0251,0.0193/
       data ys3(1:ns)/ 839.8, 860.3, 886.1, 897.1, 912.6, 933.3, 956.0,
     1  971.7, 975.1, 993.2,1025.0,1033.0,1034.0,1047.0,
     2  1068.0,1071.0,1082.0,1112.0,1112.0,1124.0,1113.0,
     3  1114.0,1146.0,1155.0,1139.0,1140.0,1155.0,1159.0,
     4  1140.0,1161.0,1173.0,1154.0,1130.0,1139.0,1129.0,
     5  1124.0,1157.0,1145.0,1090.0,1080.0,1114.0,1094.0,
     6  1057.0,1057.0,1066.0,1022.0, 995.9,1006.0,1007.0,
     7  965.7, 965.0, 948.5, 922.1, 884.1, 893.8, 875.4,
     8  836.3, 810.4, 815.1, 798.0, 763.6, 741.5, 726.6,
     9  714.7, 689.9, 669.8, 646.0, 614.0, 596.8, 591.3,
     a  560.6, 545.0, 522.1, 509.6, 490.6, 466.8, 451.6,
     b  432.6, 415.8, 400.1, 372.7, 367.3, 340.4, 325.0,
     c  312.2, 302.5, 296.3, 271.2, 259.8, 246.5, 232.2,
     d  223.8, 210.3, 203.4, 191.2, 175.0, 167.3, 158.5,
     e  151.2, 141.8, 135.8, 128.5, 119.4, 111.1, 105.6,
     f  100.2,  92.1500,  87.1100,  81.9300,  77.53,  71.55,  67.27,
     g  62.7400,  59.5500,  55.2,  51.2400,  48.4,  45.51,  41.87,
     h  39.6400,  37.2300,  34.6300,  32.2300,  30.73,  28.69,  26.5,
     i  24.5800,  24.0100,  21.8800,  20.15,  18.79,  18.08,  17.11,
     j  15.6500,  14.2100,  13.6400,  12.9,  12.43,  10.97,  10.2,
     k   9.7050,9.2600,8.9900,7.9470,7.2090,6.8830,7.0320,
     l   6.2940,5.4,5.1990,5.5290,4.7920,4.2510,4.1460,
     m   4.0980,3.7570,3.1270,2.7650,3.2130,3.2430,2.5610,
     n   2.0410,1.9220,2.4350,2.2,1.9830,1.4210,1.2500,
     o   1.3210,1.7270,1.5130,1.1050,0.8550,0.8875,0.9160,
     p   1.3,0.9452,0.6504,0.5180,0.4923,0.7302,0.8328,
     q   0.6572,0.4347,0.3194,0.3528,0.4665,0.5343,0.2550,
     r   0.2434,0.1955,0.1875,0.2753,0.3236,0.3299,0.2086,
     s   0.1413,0.1343,0.1662,0.2082,0.1784,0.1134,0.0808,
     t   0.0776,0.0714,0.0977,0.0936,0.1419,0.0916,0.0644,
     u   0.0698,0.0641,0.0539,0.0467,0.0503,0.0386,0.0285,
     v   0.0271,0.0405,0.0294/

       data xy1(1:ny)/305.,306.,307.,308.,309.,310.,311.,
     1  312.,313.,314.,315.,316.,317.,318.,
     2 319.,320.,321.,322.,323.,324.,325./
       data yy1(1:ny)/ 0.96, 0.96, 1.00, 1.09, 1.32, 1.80, 2.78,
     1  4.63, 7.80,12.60,16.70,19.40,17.10,20.70,
     2  17.20,16.30, 7.59,10.90,13.60,10.20,11.20/
       data yy2(1:ny)/5.659,16.560,47.61, 114.200,230.1,392.1,586.9,
     1 793.3 , 981.7 ,1139., 1225. , 1300. , 1295.,1365.,
     2 1282. ,1534. ,1395. ,1728. ,1701. ,1657. ,2065./

* local

      REAL yg(kw), yg1(kw), yg2(kw), yg3(kw), yg4(kw)
      REAL xso3(kw),s226(kw),s263(kw),s298(kw)
      REAL qy1d, qy3p
      REAL tau, tau2, tau3
      REAL a, b, c
      REAL a0, a1, a2, a3, a4, a5, a6, a7
      REAL xl, xl0
      REAL so3
      REAL dum
      INTEGER myld
      INTEGER kmich, kjpl87, kjpl92, kshet, kjpl97
      INTEGER i, iw, n, idum
      INTEGER ierr

****************************************************************

*************       jlabel(j) = 'O3 -> O2 + O(1D)'
*************       jlabel(j) = 'O3 -> O2 + O(3P)'

      j = j + 1
      jlabel(j) = 'O3 -> O2 + O(3P)'
      
      j = j + 1
      jlabel(j) = 'O3 -> O2 + O(1D)'
      
* cross sections:
* from WMO 1985 Ozone Assessment
* from 175.439 to 847.500 nm
* use value at 273 K

      n = nwmo
      xwork(1:n)=xs1wmo(1:n)
      ywork(1:n)=ys1wmo(1:n)
      
      CALL addpnt(xwork,ywork,kdata,n,xwork(1)*(1.-deltax),0.)
      CALL addpnt(xwork,ywork,kdata,n,               0.,0.)
      CALL addpnt(xwork,ywork,kdata,n,xwork(n)*(1.+deltax),0.)
      CALL addpnt(xwork,ywork,kdata,n,           1.e+38,0.)
      CALL inter2(nw,wl,yg,n,xwork,ywork,ierr)
      IF (ierr .NE. 0) THEN
         WRITE(*,*) ierr, jlabel(j)
         STOP
      ENDIF

      DO iw = 1, nw-1
         xso3(iw) = yg(iw)
      ENDDO

* For Hartley and Huggins bands, use temperature-dependent values from
* Molina, L. T., and M. J. Molina, Absolute absorption cross sections
* of ozone in the 185- to 350-nm wavelength range,
* J. Geophys. Res., vol. 91, 14501-14508, 1986.

      n1 = ns
      n2 = ns
      n3 = ns
      DO i = 1, n1
         xs2(i) = xs1(i)
         xs3(i) = xs1(i)
      ENDDO
      
      xwork(1:n1)=xs1(1:n1)
      ywork(1:n1)=ys1(1:n1)      
      CALL addpnt(xwork,ywork,kdata,n1,xwork(1)*(1.-deltax),0.)
      CALL addpnt(xwork,ywork,kdata,n1,               0.,0.)
      CALL addpnt(xwork,ywork,kdata,n1,xwork(n1)*(1.+deltax),0.)
      CALL addpnt(xwork,ywork,kdata,n1,            1.e+38,0.)
      CALL inter2(nw,wl,yg,n1,xwork,ywork,ierr)
      IF (ierr .NE. 0) THEN
         WRITE(*,*) ierr, jlabel(j)
         STOP 100
      ENDIF
      DO iw = 1, nw-1
         s226(iw) = yg(iw)*1.E-20
      ENDDO

      xwork(1:n2)=xs2(1:n2)
      ywork(1:n2)=ys2(1:n2)
      CALL addpnt(xwork,ywork,kdata,n2,xwork(1)*(1.-deltax),0.)
      CALL addpnt(xwork,ywork,kdata,n2,               0.,0.)
      CALL addpnt(xwork,ywork,kdata,n2,xwork(n2)*(1.+deltax),0.)
      CALL addpnt(xwork,ywork,kdata,n2,            1.e+38,0.)
      CALL inter2(nw,wl,yg,n2,xwork,ywork,ierr)
      IF (ierr .NE. 0) THEN
         WRITE(*,*) ierr, jlabel(j)
         STOP 101
      ENDIF
      DO iw = 1, nw-1
         s263(iw) = yg(iw)*1.E-20
      ENDDO

      xwork(1:n3)=xs3(1:n3)
      ywork(1:n3)=ys3(1:n3)
      CALL addpnt(xwork,ywork,kdata,n3,xwork(1)*(1.-deltax),0.)
      CALL addpnt(xwork,ywork,kdata,n3,               0.,0.)
      CALL addpnt(xwork,ywork,kdata,n3,xwork(n3)*(1.+deltax),0.)
      CALL addpnt(xwork,ywork,kdata,n3,            1.e+38,0.)
      CALL inter2(nw,wl,yg,n3,xwork,ywork,ierr)
      IF (ierr .NE. 0) THEN
         WRITE(*,*) ierr, jlabel(j)
         STOP 102
      ENDIF
      DO iw = 1, nw-1
         s298(iw) = yg(iw)*1.E-20
      ENDDO

* quantum yield:
*    kjpl97:  JPL recommendation 1997, includes tail, similar to Shetter et al.
* read parameters from JPL'97

        n1 = ny
        n2 = n1
        DO i = 1, n1
           xy2(i) = xy1(i)
        ENDDO

        xwork(1:n1)=xy1(1:n1)
        ywork(1:n1)=yy1(1:n1)        
        CALL addpnt(xwork,ywork,kdata,n1,xwork(1)*(1.-deltax),ywork(1))
        CALL addpnt(xwork,ywork,kdata,n1,               0.,ywork(1))
        CALL addpnt(xwork,ywork,kdata,n1,xwork(n1)*(1.+deltax),
     1     ywork(n1))
        CALL addpnt(xwork,ywork,kdata,n1,            1.e+38,ywork(n1))
        CALL inter2(nw,wl,yg1,n1,xwork,ywork,ierr)
        IF (ierr .NE. 0) THEN
           WRITE(*,*) ierr, jlabel(j)
           STOP
        ENDIF

        xwork(1:n2)=xy2(1:n2)
        ywork(1:n2)=yy2(1:n2)        
        CALL addpnt(xwork,ywork,kdata,n2,xwork(1)*(1.-deltax),ywork(1))
        CALL addpnt(xwork,ywork,kdata,n2,               0.,ywork(1))
        CALL addpnt(xwork,ywork,kdata,n2,xwork(n2)*(1.+deltax),
     1	 ywork(n2))
        CALL addpnt(xwork,ywork,kdata,n2,            1.e+38,ywork(n2))
        CALL inter2(nw,wl,yg2,n2,xwork,ywork,ierr)
        IF (ierr .NE. 0) THEN
           WRITE(*,*) ierr, jlabel(j)
           STOP
        ENDIF

* compute cross sections and yields at different wavelengths, altitudes:

      DO iw = 1, nw-1

         so3 = xso3(iw)

         DO i = 1, nz

            IF ( wl(iw) .GT. 240.5  .AND. wl(iw+1) .LT. 350. ) THEN
               IF (tlev(i) .LT. 263.) THEN
                  so3 = s226(iw) + (s263(iw)-s226(iw)) / (263.-226.) *
     $                 (tlev(i)-226.)
               ELSE
                  so3 = s263(iw) + (s298(iw)-s263(iw)) / (298.-263.) *
     $              (tlev(i)-263.)
               ENDIF
            ENDIF

             IF (wc(iw) .LT. 271.) THEN
                qy1d = 0.87
             ELSE IF (wc(iw) .GE. 271. .AND. wc(iw) .LT. 290.) THEN
                qy1d = 0.87 + (wc(iw)-271.)*(.95-.87)/(290.-271.)
             ELSE IF (wc(iw) .GE. 290. .AND. wc(iw) .LT. 305.) THEN
                qy1d = 0.95
             ELSE IF (wc(iw) .GE. 305. .AND. wc(iw) .LE. 325.) THEN
                qy1d = yg1(iw) * EXP ( -yg2(iw) /tlev(i) )
             ELSE
                qy1d = 0.
             ENDIF

           sq(j,i,iw) = qy1d*so3
           qy3p = 1.0 - qy1d
           sq(j-1,i,iw) = qy3p*so3

         ENDDO
      ENDDO

      END

      SUBROUTINE rn_NO2(nw,wl,wc,nz,tlev,airlev,j,sq,jlabel)

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Provide the product (cross section) x (quantum yield) for NO2            =*
*=  photolysis:                                                              =*
*=         NO2 + hv -> NO + O(3P)                                            =*
*=  Cross section from JPL94 (can also have Davidson et al.)                 =*
*=  Quantum yield from Gardiner, Sperry, and Calvert                         =*
*-----------------------------------------------------------------------------*
*=  PARAMETERS:                                                              =*
*=  NW     - INTEGER, number of specified intervals + 1 in working        (I)=*
*=           wavelength grid                                                 =*
*=  WL     - REAL, vector of lower limits of wavelength intervals in      (I)=*
*=           working wavelength grid                                         =*
*=  WC     - REAL, vector of center points of wavelength intervals in     (I)=*
*=           working wavelength grid                                         =*
*=  NZ     - INTEGER, number of altitude levels in working altitude grid  (I)=*
*=  TLEV   - REAL, temperature (K) at each specified altitude level       (I)=*
*=  AIRLEV - REAL, air density (molec/cc) at each altitude level          (I)=*
*=  J      - INTEGER, counter for number of weighting functions defined  (IO)=*
*=  SQ     - REAL, cross section x quantum yield (cm^2) for each          (O)=*
*=           photolysis reaction defined, at each defined wavelength and     =*
*=           at each defined altitude level                                  =*
*=  JLABEL - CHARACTER*40, string identifier for each photolysis reaction (O)=*
*=           defined                                                         =*
*-----------------------------------------------------------------------------*
*=  EDIT HISTORY:                                                            =*
*=  05/98  Original, adapted from former JSPEC1 subroutine                   =*
*-----------------------------------------------------------------------------*
*= This program is free software;  you can redistribute it and/or modify     =*
*= it under the terms of the GNU General Public License as published by the  =*
*= Free Software Foundation;  either version 2 of the license, or (at your   =*
*= option) any later version.                                                =*
*= The TUV package is distributed in the hope that it will be useful, but    =*
*= WITHOUT ANY WARRANTY;  without even the implied warranty of MERCHANTIBI-  =*
*= LITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public     =*
*= License for more details.                                                 =*
*= To obtain a copy of the GNU General Public License, write to:             =*
*= Free Software Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.   =*
*-----------------------------------------------------------------------------*


      IMPLICIT NONE
      INCLUDE 'tuv.params'

* input

      INTEGER nw
      REAL wl(kw), wc(kw)
      
      INTEGER nz

      REAL tlev(kz)
      REAL airlev(kz)

* weighting functions

      CHARACTER*40 jlabel(kj)
      REAL sq(kj,kz,kw)

* input/output:
      INTEGER j

* data arrays

      INTEGER kdata,ns,ny   ! number for obsorbtion, quantum yield coefficient
      PARAMETER(kdata=200,ns=57,ny=66)

      INTEGER n1
      REAL xs1(kdata), xs2(kdata), xs3(kdata), xy1(kdata)
      REAL ys1(kdata), ys2(kdata),yy1(kdata),xwork(kdata),ywork(kdata),
     1 ywork2(kdata) 
      data xs1(1:ns)/202.02,204.08,206.19,208.33,210.53,212.77,215.09,
     1   217.39,219.78,222.22,224.72,227.27,229.89,232.56,
     2   235.29,238.09,240.96,243.90,246.91,250.00,253.17,
     3   256.41,259.74,263.16,266.67,270.27,273.97,277.78,
     4   281.69,285.71,289.85,294.12,298.51,303.03,307.69,
     5   312.50,317.50,322.50,327.50,332.50,337.50,342.50,
     6   347.50,352.50,357.50,362.50,367.50,372.50,377.50,
     7   382.50,387.50,392.50,397.50,402.50,407.50,412.50,
     8   417.50/
       data xs3(1:ns)/204.08,206.19,208.33,210.53,212.77,215.09,217.39,
     1 219.78,222.22,224.72,227.27,229.89,232.56,235.29,
     2 238.09,240.96,243.90,246.91,250.00,253.17,256.41,
     3 259.74,263.16,266.67,270.27,273.97,277.78,281.69,
     4 285.71,289.85,294.12,298.51,303.03,307.69,312.50,
     5 317.50,322.50,327.50,332.50,337.50,342.50,347.50,
     6 352.50,357.50,362.50,367.50,372.50,377.50,382.50,
     7 387.50,392.50,397.50,402.50,407.50,412.50,417.50,
     8 422.50/
       data ys1(1:ns)/41.450,44.780,44.540,46.410,48.660,48.180,50.220,
     1 44.410,47.130,37.720,39.290,27.400,27.780,16.890,
     2 16.180, 8.812, 7.472, 3.909, 2.753, 2.007, 1.973,
     3 2.111, 2.357, 2.698, 3.247, 3.785, 5.030, 5.880,
     4 7.000, 8.150, 9.720,11.540,13.440,15.890,18.670,
     5 21.530,24.770,28.070,31.330,34.250,37.980,40.650,
     6 43.130,47.170,48.330,51.660,53.150,55.080,56.440,
     7 57.570,59.270,58.450,60.210,57.810,59.990,56.510,
     8 58.120/
       data ys2(1:ns)/0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     1 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     2 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
     3 0.000, 0.000, 0.000, 0.000, 0.000, 0.075, 0.082,
     4 -0.053,-0.043,-0.031,-0.162,-0.284,-0.357,-0.536,
     5 -0.686,-0.786,-1.105,-1.355,-1.277,-1.612,-1.890,
     6 -1.219,-1.921,-1.095,-1.322,-1.102,-0.806,-0.867,
     7 -0.945,-0.923,-0.738,-0.599,-0.545,-1.129, 0.001,
     8 -1.208/
       data xy1(1:ny)/1.,285.,290.,295.,300.,305.,310.,
     1 315.,320.,325.,330.,335.,340.,345.,
     2 350.,355.,360.,365.,370.,375.,380.,
     3 381.,382.,383.,384.,385.,386.,387.,
     4 388.,389.,390.,391.,392.,393.,394.,
     5 395.,396.,397.,398.,399.,400.,401.,
     6 402.,403.,404.,405.,406.,407.,408.,
     7 409.,410.,411.,412.,413.,414.,415.,
     8 416.,417.,418.,419.,420.,421.,422.,
     9 423.,424.,450./
       data yy1(1:ny)/1.000,1.000,0.999,0.998,0.997,0.996,0.995,
     1 0.994,0.993,0.992,0.991,0.990,0.989,0.988,
     2 0.987,0.986,0.984,0.983,0.981,0.979,0.975,
     3 0.974,0.973,0.972,0.971,0.969,0.967,0.966,
     4 0.964,0.962,0.960,0.959,0.957,0.953,0.950,
     5 0.942,0.922,0.870,0.820,0.760,0.695,0.635,
     6 0.560,0.485,0.425,0.350,0.290,0.225,0.185,
     7 0.153,0.130,0.110,0.094,0.083,0.070,0.059,
     8 0.048,0.039,0.030,0.023,0.018,0.012,0.008,
     9 0.004,0.000,0.000/

* local

      REAL yg1(kw), yg2(kw)
      REAL xsno2(kz,kw)
      REAL dum
      INTEGER i, iw, n, idum, ierr
      integer mabs


**************** NO2 photodissociation

      j = j + 1
      jlabel(j) = 'NO2 -> NO + O(3P)'

* cross section
*      NO2_jpl94.abs  (same as JPL97)

* read in wavelength bins, cross section at T0 and temperature correction
* coefficient a;  see input file for details.
* data need to be scaled to total area per bin so that they can be used with
* inter3
* cross section data from JPL 94 recommendation
* JPL 97 recommendation is identical

         n=ns                ! obsorption coefficient 
         DO i = 1, n
            ywork(i) = (xs3(i)-xs1(i)) * ys1(i)*1.E-20
            ywork2(i) = (xs3(i)-xs1(i)) * ys2(i)*1.E-22
            xs2(i) = xs1(i) 
         ENDDO
         CLOSE(kin)

         xs1(n+1) = xs3(n)
         xs2(n+1) = xs3(n)
         n = n+1
         n1 = n

	 xwork(1:n)=xs1(1:n)	  
         CALL inter3(nw,wl,yg1,n,xwork,ywork,0)
	 xwork(1:n1)=xs2(1:n1)
         CALL inter3(nw,wl,yg2,n1,xwork,ywork2,0)

* yg1, yg2 are per nm, so rescale by bin widths

         DO iw = 1, nw-1
            yg1(iw) = yg1(iw)/(wl(iw+1)-wl(iw))
            yg2(iw) = yg2(iw)/(wl(iw+1)-wl(iw))
         ENDDO

         DO iw = 1, nw-1
            DO i = 1, nz
               xsno2(i,iw) = yg1(iw) + yg2(iw)*(tlev(i)-273.15)
            ENDDO
         ENDDO 


* quantum yield
* from Gardiner, Sperry, and Calvert

      n = ny
      xwork(1:n)=xy1(1:n)
      ywork(1:n)=yy1(1:n)
      CALL addpnt(xwork,ywork,kdata,n,xwork(1)*(1.-deltax),ywork(1))
      CALL addpnt(xwork,ywork,kdata,n,               0.,ywork(1))
      CALL addpnt(xwork,ywork,kdata,n,xwork(n)*(1.+deltax),   0.)
      CALL addpnt(xwork,ywork,kdata,n,           1.e+38,   0.)
      CALL inter2(nw,wl,yg1,n,xwork,ywork,ierr)
      IF (ierr .NE. 0) THEN
         WRITE(*,*) ierr, jlabel(j)
         STOP
      ENDIF

* combine

      DO iw = 1, nw - 1
         DO i = 1, nz
            sq(j,i,iw) = xsno2(i,iw)*yg1(iw)
         ENDDO
      ENDDO

      END

      SUBROUTINE rn_N03(nw,wl,wc,nz,tlev,airlev,j,sq,jlabel)

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Provide the product (absorptioon cross section) x (quantum yield) for    =*
*=  both channels of NO3 photolysis:                                         =*
*=          (a) NO3 + hv -> NO + O2                                          =*
*=          (b) NO3 + hv -> NO2 + O(3P)                                      =*
*=  Cross section combined from Graham and Johnston (<600 nm) and JPL 94     =*
*=  Quantum yield from Madronich (1988)                                      =*
*-----------------------------------------------------------------------------*
*=  PARAMETERS:                                                              =*
*=  NW     - INTEGER, number of specified intervals + 1 in working        (I)=*
*=           wavelength grid                                                 =*
*=  WL     - REAL, vector of lower limits of wavelength intervals in      (I)=*
*=           working wavelength grid                                         =*
*=  WC     - REAL, vector of center points of wavelength intervals in     (I)=*
*=           working wavelength grid                                         =*
*=  NZ     - INTEGER, number of altitude levels in working altitude grid  (I)=*
*=  TLEV   - REAL, temperature (K) at each specified altitude level       (I)=*
*=  AIRLEV - REAL, air density (molec/cc) at each altitude level          (I)=*
*=  J      - INTEGER, counter for number of weighting functions defined  (IO)=*
*=  SQ     - REAL, cross section x quantum yield (cm^2) for each          (O)=*
*=           photolysis reaction defined, at each defined wavelength and     =*
*=           at each defined altitude level                                  =*
*=  JLABEL - CHARACTER*40, string identifier for each photolysis reaction (O)=*
*=           defined                                                         =*
*-----------------------------------------------------------------------------*
*=  EDIT HISTORY:                                                            =*
*=  05/98  Original, adapted from former JSPEC1 subroutine                   =*
*-----------------------------------------------------------------------------*
*= This program is free software;  you can redistribute it and/or modify     =*
*= it under the terms of the GNU General Public License as published by the  =*
*= Free Software Foundation;  either version 2 of the license, or (at your   =*
*= option) any later version.                                                =*
*= The TUV package is distributed in the hope that it will be useful, but    =*
*= WITHOUT ANY WARRANTY;  without even the implied warranty of MERCHANTIBI-  =*
*= LITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public     =*
*= License for more details.                                                 =*
*= To obtain a copy of the GNU General Public License, write to:             =*
*= Free Software Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.   =*
*-----------------------------------------------------------------------------*

      IMPLICIT NONE
      INCLUDE 'tuv.params'

* input

      INTEGER nw
      REAL wl(kw), wc(kw)
      
      INTEGER nz

      REAL tlev(kz)
      REAL airlev(kz)

* weighting functions

      CHARACTER*40 jlabel(kj)
      REAL sq(kj,kz,kw)

* input/output:
      INTEGER j

* data arrays

      INTEGER kdata,ngj,njpl
      PARAMETER(kdata=350,ngj=305,njpl=71)

      REAL x1(kdata),x1jpl(kdata),xwork(kdata),ywork(kdata)
      REAL y1gj(kdata),y1jpl(kdata)
        data y1gj(1:ngj)/0.0,  0.1,  0.1,  0.3,  0.2,  0.5,  0.3,
     1   0.1,  0.3,  0.5,  0.6,  0.5,  0.3,  0.7,
     2   0.7,  0.6,  0.3,  0.4,  0.6,  0.9,  0.9,
     3   0.9,  0.8,  1.0,  1.2,  1.3,  0.9,  0.8,
     4   1.2,  1.2,  1.2,  1.5,  1.4,  1.5,  1.7,
     5   2.1,  2.1,  1.8,  1.8,  2.1,  1.9,  1.9,
     6   2.0,  1.9,  2.1,  2.3,  2.3,  2.5,  2.8,
     7   2.8,  2.7,  2.8,  3.1,  3.2,  3.4,  3.5,
     8   3.2,  3.4,  3.7,  3.9,  3.9,  3.6,  3.5,
     9   3.8,  4.1,  4.5,  4.5,  4.8,  5.0,  5.2,
     a   4.9,  5.0,  5.4,  5.5,  5.6,  5.9,  6.4,
     b   6.8,  6.6,  6.4,  6.4,  6.5,  6.3,  6.1,
     c   6.2,  6.6,  7.4,  8.0,  8.0,  8.6,  9.3,
     d   9.2,  8.9,  8.9,  8.8,  9.1, 10.4, 11.2,
     e  10.8, 10.3,  9.8,  9.4,  9.1,  9.5, 10.5,
     f  11.6, 11.9, 11.4, 10.6, 11.2, 13.0, 15.1,
     g  16.1, 15.1, 14.1, 14.0, 14.0, 13.0, 12.1,
     h  12.8, 14.4, 15.8, 17.2, 16.6, 15.0, 13.8,
     i  13.7, 15.1, 17.9, 21.0, 20.9, 19.1, 18.1,
     j  17.3, 17.7, 20.2, 23.2, 23.8, 21.1, 18.8,
     k  18.1, 16.8, 16.8, 14.3, 13.9, 16.2, 20.4,
     l  25.6, 27.5, 24.9, 22.4, 21.4, 21.6, 22.2,
     m  24.5, 27.8, 29.5, 30.0, 31.7, 34.3, 32.3,
     n   28.5, 26.8, 25.9, 24.8, 24.7, 25.8, 25.5,
     o   25.7, 26.3, 25.3, 25.1, 24.8, 24.7, 25.5,
     p   27.0, 29.2, 30.5, 30.3, 29.4, 29.9, 32.0,
     q   31.0, 26.8, 24.7, 24.6, 27.5, 34.8, 44.8,
     r   55.2, 56.7, 51.9, 48.3, 43.2, 39.2, 39.1,
     s   41.6, 40.9, 35.4, 28.9, 24.5, 24.5, 28.4,
     t   33.9, 40.0, 41.8, 33.8, 23.2, 15.9, 13.3,
     u   13.5, 14.3, 16.9, 21.7, 22.4, 19.9, 17.4,
     v   16.7, 18.3, 20.2, 24.7, 39.8, 76.1,120.4,
     x  116.6, 86.5, 70.0, 69.0, 68.9, 67.0, 64.1,
     y  50.2, 32.7, 19.9, 13.2, 10.6, 12.3, 16.4,
     z 17.6, 13.4,  9.8,  7.8,  6.8,  6.9,  7.1,
     A  6.7,  5.6,  4.9,  4.8,  3.7,  3.2,  3.3,
     B  3.9,  4.7,  5.7,  6.9,  8.9, 11.8, 16.8,
     C 27.6, 51.2,101.5,170.8,170.4,115.4, 73.5,
     D 48.6, 29.7, 17.5, 10.7,  7.5,  6.0,  5.7,
     E  4.7,  3.6,  3.0,  3.1,  4.0,  5.5,  5.9,
     F  4.9,  3.5,  2.5,  1.6,  0.9,  0.5,  0.3,
     G  0.2,  0.4,  0.2,  0.1,  0.0,  0.0,  0.1,
     H  0.1,  0.2,  0.4,  0.4,  0.4,  0.4,  0.3,
     I  0.2,  0.2,  0.1,  0.0/
       data x1jpl(1:njpl)/600.,601.,602.,603.,604.,605.,606.,
     1 607.,608.,609.,610.,611.,612.,613.,
     2 614.,615.,616.,617.,618.,619.,620.,
     3 621.,622.,623.,624.,625.,626.,627.,
     4 628.,629.,630.,631.,632.,633.,634.,
     5 635.,636.,637.,638.,639.,640.,641.,
     6 642.,643.,644.,645.,646.,647.,648.,
     7 649.,650.,651.,652.,653.,654.,655.,
     8 656.,657.,658.,659.,660.,661.,662.,
     9 663.,664.,665.,666.,667.,668.,669.,
     a 670./
        data y1jpl(1:njpl)/ 258., 263., 302., 351., 413., 415., 322.,
     1  225., 170., 153., 192., 171., 202., 241.,
     2  242., 210., 190., 189., 208., 229., 292.,
     3  450., 941.,1407.,1139., 796., 703., 715.,
     4  702., 672., 638., 470., 344., 194., 142.,
     5  128., 159., 191., 193., 162., 121.,  99.,
     6   91.,  93.,  92.,  85.,  72.,  69.,  60.,
     7   51.,  49.,  52.,  55.,  61.,  76.,  93.,
     8  131., 172., 222., 356., 658.,1308.,2000.,
     9 1742.,1110., 752., 463., 254., 163., 113.,
     a   85./     
     
* local

      REAL yg(kw), yg1(kw)
      REAL qy
      INTEGER irow, icol
      INTEGER i, iw, n, idum
      INTEGER ierr

****************      jlabel(j) = 'NO3 -> NO + O2'
****************      jlabel(j) = 'NO3 -> NO2 + O(3P)'


* cross section
*     measurements of Graham and Johnston 1978
      
      n=ngj
       
      DO i = 1, ngj
         ywork(i) =  y1gj(i) * 1.E-19
         x1(i) = 400. + 1.*FLOAT(i-1)
      ENDDO

      CALL addpnt(x1,ywork,kdata,n,x1(1)*(1.-deltax),0.)
      CALL addpnt(x1,ywork,kdata,n,               0.,0.)
      CALL addpnt(x1,ywork,kdata,n,x1(n)*(1.+deltax),0.)
      CALL addpnt(x1,ywork,kdata,n,           1.e+38,0.)
      CALL inter2(nw,wl,yg,n,x1,ywork,ierr)
      IF (ierr .NE. 0) THEN
         WRITE(*,*) ierr, jlabel(j)
         STOP
      ENDIF

*     cross section from JPL94:
      
      n=njpl
      DO i = 1, njpl
        ywork(i) = y1jpl(i)*1E-20	
	xwork(i) = x1jpl(i)
      ENDDO 

      CALL addpnt(xwork,ywork,kdata,n,xwork(1)*(1.-deltax),0.)
      CALL addpnt(xwork,ywork,kdata,n,               0.,0.)
      CALL addpnt(xwork,ywork,kdata,n,xwork(n)*(1.+deltax),0.)
      CALL addpnt(xwork,ywork,kdata,n,           1.e+38,0.)
      CALL inter2(nw,wl,yg1,n,xwork,ywork,ierr)
      IF (ierr .NE. 0) THEN
         WRITE(*,*) ierr, jlabel(j)
         STOP
      ENDIF

* use JPL94 for wavelengths longer than 600 nm

      DO iw = 1, nw-1
         IF(wl(iw) .GT. 600.) yg(iw) = yg1(iw)
      ENDDO

* quantum yield:
* from Madronich (1988) see CEC NO3 book.

* for   NO3 ->NO+O2

      j = j + 1
      jlabel(j) = 'NO3 -> NO + O2'
      DO iw = 1, nw - 1
         IF (wc(iw).LT.584.) THEN 
            qy = 0.
         ELSEIF (wc(iw).GE.640.) THEN
            qy = 0.
         ELSEIF (wc(iw).GE.595.) THEN 
            qy = 0.35*(1.-(wc(iw)-595.)/45.)
         ELSE
            qy = 0.35*(wc(iw)-584.)/11.
         ENDIF
         DO i = 1, nz
            sq(j,i,iw) = yg(iw)*qy
         ENDDO
      ENDDO

* for  NO3 ->NO2+O
      j = j + 1
      jlabel(j) = 'NO3 -> NO2 + O(3P)'
      DO iw = 1, nw - 1
         IF (wc(iw).LT.584.) THEN
            qy = 1.
         ELSEIF (wc(iw).GT.640.) THEN
            qy = 0.
         ELSEIF (wc(iw).GT.595.) THEN
            qy = 0.65*(1-(wc(iw)-595.)/45.)
         ELSE
            qy = 1.-0.35*(wc(iw)-584.)/11.
         ENDIF
         DO i = 1, nz
            sq(j,i,iw) = yg(iw)*qy
         ENDDO
      ENDDO

      END


      SUBROUTINE rn_HNO2(nw,wl,wc,nz,tlev,airlev,j,sq,jlabel)

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Provide the product (cross section) x (quantum yield) for HNO2 photolysis=*
*=     HNO2 + hv -> NO + OH                                                  =*
*=     HNO2 + hv -> NO2 + HO2                                                =*
*=  Cross section:  from JPL97                                               =*
*=  Quantum yield:  assumed to be unity                                      =*
*-----------------------------------------------------------------------------*
*=  PARAMETERS:                                                              =*
*=  NW     - INTEGER, number of specified intervals + 1 in working        (I)=*
*=           wavelength grid                                                 =*
*=  WL     - REAL, vector of lower limits of wavelength intervals in      (I)=*
*=           working wavelength grid                                         =*
*=  WC     - REAL, vector of center points of wavelength intervals in     (I)=*
*=           working wavelength grid                                         =*
*=  NZ     - INTEGER, number of altitude levels in working altitude grid  (I)=*
*=  TLEV   - REAL, temperature (K) at each specified altitude level       (I)=*
*=  AIRLEV - REAL, air density (molec/cc) at each altitude level          (I)=*
*=  J      - INTEGER, counter for number of weighting functions defined  (IO)=*
*=  SQ     - REAL, cross section x quantum yield (cm^2) for each          (O)=*
*=           photolysis reaction defined, at each defined wavelength and     =*
*=           at each defined altitude level                                  =*
*=  JLABEL - CHARACTER*40, string identifier for each photolysis reaction (O)=*
*=           defined                                                         =*
*-----------------------------------------------------------------------------*
*=  EDIT HISTORY:                                                            =*
*=  05/98  Original, adapted from former JSPEC1 subroutine                   =*
*=  07/01  tyh add HNO2 + hv -> NO2 + HO2                                    =* 
*-----------------------------------------------------------------------------*
*= This program is free software;  you can redistribute it and/or modify     =*
*= it under the terms of the GNU General Public License as published by the  =*
*= Free Software Foundation;  either version 2 of the license, or (at your   =*
*= option) any later version.                                                =*
*= The TUV package is distributed in the hope that it will be useful, but    =*
*= WITHOUT ANY WARRANTY;  without even the implied warranty of MERCHANTIBI-  =*
*= LITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public     =*
*= License for more details.                                                 =*
*= To obtain a copy of the GNU General Public License, write to:             =*
*= Free Software Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.   =*
*-----------------------------------------------------------------------------*

      IMPLICIT NONE
      INCLUDE 'tuv.params'

* input

      INTEGER nw
      REAL wl(kw), wc(kw)
      
      INTEGER nz

      REAL tlev(kz)
      REAL airlev(kz)

* weighting functions

      CHARACTER*40 jlabel(kj)
      REAL sq(kj,kz,kw)

* input/output:
      INTEGER j

* data arrays

      INTEGER kdata,nno,nno2
      PARAMETER(kdata=100,nno=89,nno2=58)

      REAL x1(kdata),x1no(kdata),x1no2(kdata),s1no(kdata),s1no2(kdata)
      REAL y1no(kdata),y1no2(kdata),xwork(kdata),ywork(kdata)

       data x1no(1:nno)/309.0,310.0,311.0,312.0,313.0,314.0,315.0,
     1 316.0,317.0,318.0,319.0,320.0,321.0,322.0,
     2 323.0,324.0,325.0,326.0,327.0,328.0,329.0,
     3 330.0,331.0,332.0,333.0,334.0,335.0,336.0,
     4 337.0,338.0,339.0,340.0,341.0,342.0,343.0,
     5 344.0,345.0,346.0,347.0,348.0,349.0,350.0,
     6 351.0,352.0,353.0,354.0,355.0,356.0,357.0,
     7 358.0,359.0,360.0,361.0,362.0,363.0,364.0,
     8 365.0,366.0,367.0,368.0,369.0,370.0,371.0,
     9 372.0,373.0,374.0,375.0,376.0,377.0,378.0,
     a 379.0,380.0,381.0,382.0,383.0,384.0,385.0,
     b 386.0,387.0,388.0,389.0,390.0,391.0,392.0,
     c 393.0,394.0,395.0,396.0,400.0/
       data s1no(1:nno)/.000E+00,.130E-19,.190E-19,.280E-19,.220E-19,
     1 .360E-19,.300E-19,.140E-19,.310E-19,.560E-19,.360E-19,.490E-19,
     2 .780E-19,.490E-19,.510E-19,.710E-19,.500E-19,.290E-19,.660E-19,
     3 .117E-18,.610E-19,.111E-18,.179E-18,.870E-19,.760E-19,.960E-19,
     4 .960E-19,.720E-19,.530E-19,.100E-18,.188E-18,.100E-18,.170E-18,
     5 .386E-18,.149E-18,.970E-19,.109E-18,.123E-18,.104E-18,.910E-19,
     6 .790E-19,.112E-18,.212E-18,.155E-18,.191E-18,.581E-18,.364E-18,
     7 .141E-18,.117E-18,.120E-18,.104E-18,.900E-19,.830E-19,.800E-19,
     8 .960E-19,.146E-18,.168E-18,.183E-18,.302E-18,.520E-18,.388E-18,
     9 .178E-18,.113E-18,.100E-18,.770E-19,.620E-19,.530E-19,.530E-19,
     a .500E-19,.580E-19,.800E-19,.960E-19,.113E-18,.159E-18,.210E-18,
     b .241E-18,.203E-18,.134E-18,.900E-19,.560E-19,.340E-19,.270E-19,
     c .200E-19,.150E-19,.110E-19,.600E-20,.100E-19,.400E-20,.000E+00/
       data y1no(1:nno)/0.410,0.410,0.411,0.421,0.432,0.443,0.454,
     1 0.464,0.475,0.486,0.496,0.507,0.518,0.529,
     2 0.539,0.550,0.561,0.571,0.582,0.593,0.604,
     3 0.614,0.625,0.636,0.646,0.657,0.668,0.679,
     4 0.689,0.700,0.711,0.721,0.732,0.743,0.754,
     5 0.764,0.775,0.786,0.796,0.807,0.818,0.829,
     6 0.839,0.850,0.861,0.871,0.882,0.893,0.904,
     7 0.914,0.925,0.936,0.946,0.957,0.968,0.979,
     8 0.989,1.000,1.000,1.000,1.000,1.000,1.000,
     9 1.000,1.000,1.000,1.000,1.000,1.000,1.000,
     a 1.000,1.000,1.000,1.000,1.000,1.000,1.000,
     b 1.000,1.000,1.000,1.000,1.000,1.000,1.000,
     c 1.000,1.000,1.000,1.000,1.000/
       data x1no2(1:nno2)/309.0,310.0,311.0,312.0,313.0,314.0,315.0,
     1 316.0,317.0,318.0,319.0,320.0,321.0,322.0,
     2 323.0,324.0,325.0,326.0,327.0,328.0,329.0,
     3 330.0,331.0,332.0,333.0,334.0,335.0,336.0,
     4 337.0,338.0,339.0,340.0,341.0,342.0,343.0,
     5 344.0,345.0,346.0,347.0,348.0,349.0,350.0,
     6 351.0,352.0,353.0,354.0,355.0,356.0,357.0,
     7 358.0,359.0,360.0,361.0,362.0,363.0,364.0,
     8 365.0,366.0/
       data s1no2(1:nno2)/.000E+00,.130E-19,.190E-19,.280E-19,.220E-19,
     1 .360E-19,.300E-19,.140E-19,.310E-19,.560E-19,.360E-19,.490E-19,
     2 .780E-19,.490E-19,.510E-19,.710E-19,.500E-19,.290E-19,.660E-19,
     3 .117E-18,.610E-19,.111E-18,.179E-18,.870E-19,.760E-19,.960E-19,
     4 .960E-19,.720E-19,.530E-19,.100E-18,.188E-18,.100E-18,.170E-18,
     5 .386E-18,.149E-18,.970E-19,.109E-18,.123E-18,.104E-18,.910E-19,
     6 .790E-19,.112E-18,.212E-18,.155E-18,.191E-18,.581E-18,.364E-18,
     7 .141E-18,.117E-18,.120E-18,.104E-18,.900E-19,.830E-19,.800E-19,
     8 .960E-19,.146E-18,.168E-18,.183E-18/
       data y1no2(1:nno2)/0.590,0.590,0.589,0.579,0.568,0.557,0.546,
     1 0.536,0.525,0.514,0.504,0.493,0.482,0.471,
     2 0.461,0.450,0.439,0.429,0.418,0.407,0.396,
     3 0.386,0.375,0.364,0.354,0.343,0.332,0.321,
     4 0.311,0.300,0.289,0.279,0.268,0.257,0.246,
     5 0.236,0.225,0.214,0.204,0.193,0.182,0.171,
     6 0.161,0.150,0.139,0.129,0.118,0.107,0.096,
     7 0.086,0.075,0.064,0.054,0.043,0.032,0.021,
     8 0.011,0.000/

* local

      REAL yg(kw),yg1(kw)
      REAL qy
      INTEGER i, iw, n
      INTEGER ierr

**************** HNO2-NO photodissociation
* cross section 
      j=j+1
      jlabel(j)='HNO2 + hv -> NO + OH'
      n = nno      
      x1(1:n)=x1no(1:n) 
      ywork(1:n)=s1no(1:n) 
      
      CALL addpnt(x1,ywork,kdata,n,x1(1)*(1.-deltax),0.)
      CALL addpnt(x1,ywork,kdata,n,               0.,0.)
      CALL addpnt(x1,ywork,kdata,n,x1(n)*(1.+deltax),0.)
      CALL addpnt(x1,ywork,kdata,n,           1.e+38,0.)
      CALL inter2(nw,wl,yg,n,x1,ywork,ierr)
      IF (ierr .NE. 0) THEN
         WRITE(*,*) ierr, jlabel(j)
         STOP
      ENDIF

* quantum yield 
      n=nno
      xwork(1:n)=x1no(1:n) 
      ywork(1:n)=y1no(1:n) 
       
      CALL addpnt(xwork,ywork,kdata,n,xwork(1)*(1.-deltax),0.)
      CALL addpnt(xwork,ywork,kdata,n,               0.,0.)
      CALL addpnt(xwork,ywork,kdata,n,xwork(n)*(1.+deltax),0.)
      CALL addpnt(xwork,ywork,kdata,n,           1.e+38,0.)
      CALL inter2(nw,wl,yg1,n,xwork,ywork,ierr)
      IF (ierr .NE. 0) THEN
         WRITE(*,*) ierr, jlabel(j)
         STOP
      ENDIF

      DO iw = 1, nw - 1
         DO i = 1, nz
            sq(j,i,iw) = yg(iw)*yg1(iw)
         ENDDO
      ENDDO
      
**************** HNO2-NO2 photodissociation
* cross section 
      j=j+1
      jlabel(j)='HNO2 + hv -> NO2 + HO2'
      n=nno2
      x1(1:n)=x1no2(1:n)
      ywork(1:n)=s1no2(1:n)
      
      CALL addpnt(x1,ywork,kdata,n,x1(1)*(1.-deltax),0.)
      CALL addpnt(x1,ywork,kdata,n,               0.,0.)
      CALL addpnt(x1,ywork,kdata,n,x1(n)*(1.+deltax),0.)
      CALL addpnt(x1,ywork,kdata,n,           1.e+38,0.)
      CALL inter2(nw,wl,yg,n,x1,ywork,ierr)
      IF (ierr .NE. 0) THEN
         WRITE(*,*) ierr, jlabel(j)
         STOP
      ENDIF

* quantum yield 
      n=nno2
      xwork(1:n)=x1no2(1:n)
      ywork(1:n)=y1no2(1:n)
      
      CALL addpnt(xwork,ywork,kdata,n,xwork(1)*(1.-deltax),0.)
      CALL addpnt(xwork,ywork,kdata,n,               0.,0.)
      CALL addpnt(xwork,ywork,kdata,n,xwork(n)*(1.+deltax),0.)
      CALL addpnt(xwork,ywork,kdata,n,           1.e+38,0.)
      CALL inter2(nw,wl,yg1,n,xwork,ywork,ierr)
      IF (ierr .NE. 0) THEN
         WRITE(*,*) ierr, jlabel(j)
         STOP
      ENDIF

      DO iw = 1, nw - 1
         DO i = 1, nz
            sq(j,i,iw) = yg(iw)*yg1(iw)
         ENDDO
      ENDDO
            
      END
      
      SUBROUTINE rn_HNO3(nw,wl,wc,nz,tlev,airlev,j,sq,jlabel)

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Provide product of (cross section) x (quantum yield) for HNO3 photolysis =*
*=        HNO3 + hv -> OH + NO2                                              =*
*=  Cross section: Burkholder et al., 1993                                   =*
*=  Quantum yield: Assumed to be unity                                       =*
*-----------------------------------------------------------------------------*
*=  PARAMETERS:                                                              =*
*=  NW     - INTEGER, number of specified intervals + 1 in working        (I)=*
*=           wavelength grid                                                 =*
*=  WL     - REAL, vector of lower limits of wavelength intervals in      (I)=*
*=           working wavelength grid                                         =*
*=  WC     - REAL, vector of center points of wavelength intervals in     (I)=*
*=           working wavelength grid                                         =*
*=  NZ     - INTEGER, number of altitude levels in working altitude grid  (I)=*
*=  TLEV   - REAL, temperature (K) at each specified altitude level       (I)=*
*=  AIRLEV - REAL, air density (molec/cc) at each altitude level          (I)=*
*=  J      - INTEGER, counter for number of weighting functions defined  (IO)=*
*=  SQ     - REAL, cross section x quantum yield (cm^2) for each          (O)=*
*=           photolysis reaction defined, at each defined wavelength and     =*
*=           at each defined altitude level                                  =*
*=  JLABEL - CHARACTER*40, string identifier for each photolysis reaction (O)=*
*=           defined                                                         =*
*-----------------------------------------------------------------------------*
*=  EDIT HISTORY:                                                            =*
*=  05/98  Original, adapted from former JSPEC1 subroutine                   =*
*-----------------------------------------------------------------------------*
*= This program is free software;  you can redistribute it and/or modify     =*
*= it under the terms of the GNU General Public License as published by the  =*
*= Free Software Foundation;  either version 2 of the license, or (at your   =*
*= option) any later version.                                                =*
*= The TUV package is distributed in the hope that it will be useful, but    =*
*= WITHOUT ANY WARRANTY;  without even the implied warranty of MERCHANTIBI-  =*
*= LITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public     =*
*= License for more details.                                                 =*
*= To obtain a copy of the GNU General Public License, write to:             =*
*= Free Software Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.   =*
*-----------------------------------------------------------------------------*

      IMPLICIT NONE
      INCLUDE 'tuv.params'

* input

      INTEGER nw
      REAL wl(kw), wc(kw)
      
      INTEGER nz

      REAL tlev(kz)
      REAL airlev(kz)

* weighting functions

      CHARACTER*40 jlabel(kj)
      REAL sq(kj,kz,kw)

* input/output:
      INTEGER j

* data arrays

      INTEGER kdata,ns
      PARAMETER(kdata=100,ns=83)

      INTEGER n1, n2
      REAL x1(kdata), x2(kdata)
      REAL y1(kdata), y2(kdata), ywork(kdata)
       
      data y1(1:ns)/1580.0,1480.0,1360.0,1225.0,1095.0, 940.0, 770.0,
     1  0.588E+03,0.447E+03, 0.328E+03,0.231E+03,0.156E+03,0.104E+03,
     2  0.675E+02,0.439E+02,0.292E+02,0.200E+02,0.149E+02,0.118E+02,
     3  0.961E+01,0.802E+01,0.682E+01,0.575E+01,0.487E+01,0.414E+01,
     4  0.336E+01,0.293E+01,0.258E+01,0.234E+01,0.216E+01,0.206E+01,
     5  0.200E+01,0.197E+01,0.196E+01,0.195E+01,0.195E+01,0.193E+01,
     6  0.191E+01,0.187E+01,0.183E+01,0.177E+01,0.170E+01,0.162E+01,
     7  0.153E+01,0.144E+01,0.133E+01,0.123E+01,0.112E+01,0.101E+01,
     8  0.909E+00,0.807E+00,0.709E+00,0.615E+00,0.532E+00,0.453E+00,
     9  0.381E+00,0.316E+00,0.263E+00,0.208E+00,0.167E+00,0.133E+00,
     a  0.105E+00,0.814E-01,0.628E-01,0.468E-01,0.362E-01,0.271E-01,
     b  0.197E-01,0.154E-01,0.108E-01,0.820E-02,0.613E-02,0.431E-02,
     c  0.319E-02,0.243E-02,0.196E-02,0.142E-02,0.103E-02,0.861E-03,
     d  0.694E-03,0.501E-03,0.415E-03,0.417E-03/
       data y2(1:ns)/1.70, 1.70, 1.70, 1.70, 1.70, 1.70, 1.65,
     1  1.66, 1.69, 1.74, 1.77, 1.85, 1.97, 2.08,
     2  2.17, 2.17, 2.21, 2.15, 2.06, 1.96, 1.84,
     3  1.78, 1.80, 1.86, 1.90, 1.97, 1.97, 1.97,
     4  1.88, 1.75, 1.61, 1.44, 1.34, 1.23, 1.18,
     5  1.14, 1.12, 1.14, 1.14, 1.18, 1.22, 1.25,
     6  1.45, 1.49, 1.56, 1.64, 1.69, 1.78, 1.87,
     7  1.94, 2.04, 2.15, 2.27, 2.38, 2.52, 2.70,
     8  2.92, 3.10, 3.24, 3.52, 3.77, 3.91, 4.23,
     9  4.70, 5.15, 5.25, 5.74, 6.45, 6.70, 7.16,
     a  7.55, 8.16, 9.75, 9.93, 9.60,10.50,10.80,
     b  11.80,11.80, 9.30,12.10,11.90, 9.30/
* local

      REAL yg1(kw), yg2(kw)
      INTEGER i, iw
      INTEGER ierr

**************** HNO3 photodissociation

       j = j + 1
       jlabel(j) = 'HNO3 -> OH + NO2'


* HNO3 cross section parameters from Burkholder et al. 1993

      n1 = ns
      n2 = n1
      DO i = 1, n1
         x1(i) = 184. + i*2.
         x2(i) = x1(i)
      END DO
      
      ywork(1:n1)=y1(1:n1) 
      CALL addpnt(x1,ywork,kdata,n1,x1(1)*(1.-deltax),0.)
      CALL addpnt(x1,ywork,kdata,n1,               0.,0.)
      CALL addpnt(x1,ywork,kdata,n1,x1(n1)*(1.+deltax),0.)
      CALL addpnt(x1,ywork,kdata,n1,            1.e+38,0.)
      CALL inter2(nw,wl,yg1,n1,x1,ywork,ierr)
      IF (ierr .NE. 0) THEN
         WRITE(*,*) ierr, jlabel(j)
         STOP
      ENDIF

      ywork(1:n2)=y2(1:n2) 
      CALL addpnt(x2,ywork,kdata,n2,x2(1)*(1.-deltax),ywork(1))
      CALL addpnt(x2,ywork,kdata,n2,               0.,ywork(1))
      CALL addpnt(x2,ywork,kdata,n2,x2(n2)*(1.+deltax),ywork(n2))
      CALL addpnt(x2,ywork,kdata,n2,            1.e+38,ywork(n2))
      CALL inter2(nw,wl,yg2,n2,x2,ywork,ierr)
      IF (ierr .NE. 0) THEN
         WRITE(*,*) ierr, jlabel(j)
         STOP
      ENDIF

* quantum yield = 1
* correct for temperature dependence

      DO iw = 1, nw - 1
         DO i = 1, nz
            sq(j,i,iw) = yg1(iw) * 1.E-20
     $           * exp( yg2(iw)/1.e3*(tlev(i)-298.) )
         ENDDO
      ENDDO

      END


      SUBROUTINE rn_HNO4(nw,wl,wc,nz,tlev,airlev,j,sq,jlabel)

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Provide product of (cross section) x (quantum yield) for HNO4 photolysis =*
*=       HNO4 + hv -> HO2 + NO2                                              =*
*=  Cross section:  from JPL97                                               =*
*=  Quantum yield:  Assumed to be unity                                      =*
*-----------------------------------------------------------------------------*
*=  PARAMETERS:                                                              =*
*=  NW     - INTEGER, number of specified intervals + 1 in working        (I)=*
*=           wavelength grid                                                 =*
*=  WL     - REAL, vector of lower limits of wavelength intervals in      (I)=*
*=           working wavelength grid                                         =*
*=  WC     - REAL, vector of center points of wavelength intervals in     (I)=*
*=           working wavelength grid                                         =*
*=  NZ     - INTEGER, number of altitude levels in working altitude grid  (I)=*
*=  TLEV   - REAL, temperature (K) at each specified altitude level       (I)=*
*=  AIRLEV - REAL, air density (molec/cc) at each altitude level          (I)=*
*=  J      - INTEGER, counter for number of weighting functions defined  (IO)=*
*=  SQ     - REAL, cross section x quantum yield (cm^2) for each          (O)=*
*=           photolysis reaction defined, at each defined wavelength and     =*
*=           at each defined altitude level                                  =*
*=  JLABEL - CHARACTER*40, string identifier for each photolysis reaction (O)=*
*=           defined                                                         =*
*-----------------------------------------------------------------------------*
*=  EDIT HISTORY:                                                            =*
*=  05/98  Original, adapted from former JSPEC1 subroutine                   =*
*-----------------------------------------------------------------------------*
*= This program is free software;  you can redistribute it and/or modify     =*
*= it under the terms of the GNU General Public License as published by the  =*
*= Free Software Foundation;  either version 2 of the license, or (at your   =*
*= option) any later version.                                                =*
*= The TUV package is distributed in the hope that it will be useful, but    =*
*= WITHOUT ANY WARRANTY;  without even the implied warranty of MERCHANTIBI-  =*
*= LITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public     =*
*= License for more details.                                                 =*
*= To obtain a copy of the GNU General Public License, write to:             =*
*= Free Software Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.   =*
*-----------------------------------------------------------------------------*

      IMPLICIT NONE
      INCLUDE 'tuv.params'

* input

      INTEGER nw
      REAL wl(kw), wc(kw)
      
      INTEGER nz

      REAL tlev(kz)
      REAL airlev(kz)

* weighting functions

      CHARACTER*40 jlabel(kj)
      REAL sq(kj,kz,kw)

* input/output:
      INTEGER j

* data arrays

      INTEGER kdata,ns
      PARAMETER(kdata=100,ns=30)

      REAL x1(kdata),y1(kdata),xwork(kdata),ywork(kdata)
       data x1(1:ns)/190.0,195.0,200.0,205.0,210.0,215.0,220.0,
     1  225.0,230.0,235.0,240.0,245.0,250.0,255.0,
     2  260.0,265.0,270.0,275.0,280.0,285.0,290.0,
     3  295.0,300.0,305.0,310.0,315.0,320.0,325.0,
     4  330.0,335.0/
       data y1(1:ns)/.101E-16,.816E-17,.563E-17,.367E-17,.239E-17,
     1 .161E-17,.118E-17,.932E-18,.788E-18,.680E-18,.579E-18,.497E-18,
     2 .411E-18,.349E-18,.284E-18,.229E-18,.180E-18,.133E-18,.930E-19,
     3 .620E-19,.390E-19,.240E-19,.140E-19,.850E-20,.530E-20,.390E-20,
     4 .240E-20,.150E-20,.900E-21,.000E+00/
 
C* local
 
      REAL yg(kw)
      REAL qy
      INTEGER i, iw, n
      INTEGER ierr

**************** HNO4 photodissociation

* cross section 

      j = j + 1
      jlabel(j) = 'HNO4 -> 0.61(HO2 + NO2) + 0.39(OH+NO3)'

      n = ns
      xwork(1:n)=x1(1:n)
      ywork(1:n)=y1(1:n)       
      CALL addpnt(xwork,ywork,kdata,n,xwork(1)*(1.-deltax),0.)
      CALL addpnt(xwork,ywork,kdata,n,               0.,0.)
      CALL addpnt(xwork,ywork,kdata,n,xwork(n)*(1.+deltax),0.)
      CALL addpnt(xwork,ywork,kdata,n,           1.e+38,0.)
      CALL inter2(nw,wl,yg,n,xwork,ywork,ierr)
      IF (ierr .NE. 0) THEN
         WRITE(*,*) ierr, jlabel(j)
         STOP
      ENDIF

* quantum yield = 1

      qy = 1.
      DO iw = 1, nw - 1
         DO i = 1, nz
            sq(j,i,iw) = yg(iw)*qy
         ENDDO
      ENDDO
      END
      
      
      SUBROUTINE rn_H2O2(nw,wl,wc,nz,tlev,airlev,j,sq,jlabel)

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Provide product of (cross section) x (quantum yield) for H2O2 photolysis =*
*=         H2O2 + hv -> 2 OH                                                 =*
*=  Cross section:  From JPL97, tabulated values @ 298K for <260nm, T-depend.=*
*=                  parameterization for 260-350nm                           =*
*=  Quantum yield:  Assumed to be unity                                      =*
*-----------------------------------------------------------------------------*
*=  PARAMETERS:                                                              =*
*=  NW     - INTEGER, number of specified intervals + 1 in working        (I)=*
*=           wavelength grid                                                 =*
*=  WL     - REAL, vector of lower limits of wavelength intervals in      (I)=*
*=           working wavelength grid                                         =*
*=  WC     - REAL, vector of center points of wavelength intervals in     (I)=*
*=           working wavelength grid                                         =*
*=  NZ     - INTEGER, number of altitude levels in working altitude grid  (I)=*
*=  TLEV   - REAL, temperature (K) at each specified altitude level       (I)=*
*=  AIRLEV - REAL, air density (molec/cc) at each altitude level          (I)=*
*=  J      - INTEGER, counter for number of weighting functions defined  (IO)=*
*=  SQ     - REAL, cross section x quantum yield (cm^2) for each          (O)=*
*=           photolysis reaction defined, at each defined wavelength and     =*
*=           at each defined altitude level                                  =*
*=  JLABEL - CHARACTER*40, string identifier for each photolysis reaction (O)=*
*=           defined                                                         =*
*-----------------------------------------------------------------------------*
*=  EDIT HISTORY:                                                            =*
*=  05/98  Original, adapted from former JSPEC1 subroutine                   =*
*-----------------------------------------------------------------------------*
*= This program is free software;  you can redistribute it and/or modify     =*
*= it under the terms of the GNU General Public License as published by the  =*
*= Free Software Foundation;  either version 2 of the license, or (at your   =*
*= option) any later version.                                                =*
*= The TUV package is distributed in the hope that it will be useful, but    =*
*= WITHOUT ANY WARRANTY;  without even the implied warranty of MERCHANTIBI-  =*
*= LITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public     =*
*= License for more details.                                                 =*
*= To obtain a copy of the GNU General Public License, write to:             =*
*= Free Software Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.   =*
*-----------------------------------------------------------------------------*

      IMPLICIT NONE
      INCLUDE 'tuv.params'

* input

      INTEGER nw
      REAL wl(kw), wc(kw)
      
      INTEGER nz

      REAL tlev(kz)
      REAL airlev(kz)

* weighting functions

      CHARACTER*40 jlabel(kj)
      REAL sq(kj,kz,kw)

* input/output:
      INTEGER j

* data arrays

      INTEGER kdata,ns
      PARAMETER(kdata=100,ns=33)

C     INTEGER n1, n2, n3, n4, n5
      REAL x1(kdata),y1(kdata),xwork(kdata),ywork(kdata)
      data x1(1:ns)/190.,195.,200.,205.,210.,215.,220.,
     1  225.,230.,235.,240.,245.,250.,255.,
     2  260.,265.,270.,275.,280.,285.,290.,
     3  295.,300.,305.,310.,315.,320.,325.,
     4  330.,335.,340.,345.,350./
       data y1(1:ns)/67.20,56.40,47.50,40.80,35.70,30.70,25.80,
     1   21.70,18.20,15.00,12.40,10.20, 8.30, 6.70,
     2   5.30, 4.20, 3.30, 2.60, 2.00, 1.50, 1.20,
     3   0.90, 0.68, 0.51, 0.39, 0.29, 0.22, 0.16,
     4   0.13, 0.10, 0.07, 0.05, 0.04/

* local

      REAL yg(kw)
      REAL qy
      REAL a0, a1, a2, a3, a4, a5, a6, a7
      REAL b0, b1, b2, b3, b4
      REAL xs
      REAL t
      INTEGER i, iw, n, idum
      INTEGER ierr
      REAL lambda
      REAL sumA, sumB, chi

**************** H2O2 photodissociation

* cross section from Lin et al. 1978

      j = j + 1
      jlabel(j) = 'H2O2 -> 2 OH'

* cross section from JPL94 (identical to JPL97)
* tabulated data up to 260 nm

      n=ns
      DO i = 1, n
         ywork(i) = y1(i) * 1.E-20
	 xwork(i) = x1(i)
      ENDDO

      CALL addpnt(xwork,ywork,kdata,n,xwork(1)*(1.-deltax),0.)
      CALL addpnt(xwork,ywork,kdata,n,               0.,0.)
      CALL addpnt(xwork,ywork,kdata,n,xwork(n)*(1.+deltax),0.)
      CALL addpnt(xwork,ywork,kdata,n,           1.e+38,0.)
      CALL inter2(nw,wl,yg,n,xwork,ywork,ierr)
      IF (ierr .NE. 0) THEN
         WRITE(*,*) ierr, jlabel(j)
         STOP
      ENDIF

      A0 = 6.4761E+04            
      A1 = -9.2170972E+02        
      A2 = 4.535649              
      A3 = -4.4589016E-03        
      A4 = -4.035101E-05         
      A5 = 1.6878206E-07
      A6 = -2.652014E-10
      A7 = 1.5534675E-13

      B0 = 6.8123E+03
      B1 = -5.1351E+01
      B2 = 1.1522E-01
      B3 = -3.0493E-05
      B4 = -1.0924E-07

* quantum yield = 1

      qy = 1.

      DO iw = 1, nw - 1

* Parameterization (JPL94)
* Range 260-350 nm; 200-400 K

         IF ((wl(iw) .GE. 260.) .AND. (wl(iw) .LT. 350.)) THEN

           lambda = wc(iw)
           sumA = ((((((A7*lambda + A6)*lambda + A5)*lambda + 
     >                  A4)*lambda +A3)*lambda + A2)*lambda + 
     >                  A1)*lambda + A0
           sumB = (((B4*lambda + B3)*lambda + B2)*lambda + 
     >               B1)*lambda + B0

           DO i = 1, nz
              t = MIN(MAX(tlev(i),200.),400.)            
              chi = 1./(1.+EXP(-1265./t))
              xs = (chi * sumA + (1.-chi)*sumB)*1E-21
              sq(j,i,iw) = xs*qy
           ENDDO
         ELSE
           DO i = 1, nz
              sq(j,i,iw) = yg(iw)*qy
           ENDDO
         ENDIF

      ENDDO

      END


      SUBROUTINE rn_HCHO(nw,wl,wc,nz,tlev,airlev,j,sq,jlabel)

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Provide product of (cross section) x (quantum yield) for CH2O photolysis =*
*=        (a) CH2O + hv -> 2HO2 + CO                                         =*
*=        (b) CH2O + hv -> CO                                                =* 
*=  Cross section: Choice between                                            =*
*=                 1) Bass et al., 1980 (resolution: 0.025 nm)               =*
*=                 2) Moortgat and Schneider (resolution: 1 nm)              =*
*=                 3) Cantrell et al. (orig res.) for > 301 nm,              =*
*=                    IUPAC 92, 97 elsewhere                                 =*
*=                 4) Cantrell et al. (2.5 nm res.) for > 301 nm,            =*
*=                    IUPAC 92, 97 elsewhere                                 =*
*=                 5) Rogers et al., 1990                                    =*
*=                 6) new NCAR recommendation, based on averages of          =*
*=                    Cantrell et al., Moortgat and Schneider, and Rogers    =*
*=                    et al.                                                 =*
*=  Quantum yield: Choice between                                            =*
*=                 1) Evaluation by Madronich 1991 (unpublished)             =*
*=                 2) IUPAC 89, 92, 97                                       =*
*=                 3) Madronich, based on 1), updated 1998.                  =*
*=                 tyh modify for SAPRC99   07/01                            =*
*-----------------------------------------------------------------------------*
*=  PARAMETERS:                                                              =*
*=  NW     - INTEGER, number of specified intervals + 1 in working        (I)=*
*=           wavelength grid                                                 =*
*=  WL     - REAL, vector of lower limits of wavelength intervals in      (I)=*
*=           working wavelength grid                                         =*
*=  WC     - REAL, vector of center points of wavelength intervals in     (I)=*
*=           working wavelength grid                                         =*
*=  NZ     - INTEGER, number of altitude levels in working altitude grid  (I)=*
*=  TLEV   - REAL, temperature (K) at each specified altitude level       (I)=*
*=  AIRLEV - REAL, air density (molec/cc) at each altitude level          (I)=*
*=  J      - INTEGER, counter for number of weighting functions defined  (IO)=*
*=  SQ     - REAL, cross section x quantum yield (cm^2) for each          (O)=*
*=           photolysis reaction defined, at each defined wavelength and     =*
*=           at each defined altitude level                                  =*
*=  JLABEL - CHARACTER*40, string identifier for each photolysis reaction (O)=*
*=           defined                                                         =*
*-----------------------------------------------------------------------------*
*=  EDIT HISTORY:                                                            =*
*=  05/98  Original, adapted from former JSPEC1 subroutine                   =*
*-----------------------------------------------------------------------------*
*= This program is free software;  you can redistribute it and/or modify     =*
*= it under the terms of the GNU General Public License as published by the  =*
*= Free Software Foundation;  either version 2 of the license, or (at your   =*
*= option) any later version.                                                =*
*= The TUV package is distributed in the hope that it will be useful, but    =*
*= WITHOUT ANY WARRANTY;  without even the implied warranty of MERCHANTIBI-  =*
*= LITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public     =*
*= License for more details.                                                 =*
*= To obtain a copy of the GNU General Public License, write to:             =*
*= Free Software Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.   =*
*-----------------------------------------------------------------------------*

      IMPLICIT NONE
      INCLUDE 'tuv.params'

      INTEGER kdata
      PARAMETER(kdata=160)

* input

      INTEGER nw
      REAL wl(kw), wc(kw)
      
      INTEGER nz

      REAL tlev(kz)
      REAL airlev(kz)

* weighting functions

      CHARACTER*40 jlabel(kj)
      REAL sq(kj,kz,kw)

* input/output:

      INTEGER j, iz, iw

* data arrays

      INTEGER n,nm,nr
      parameter(nr=101,nm=121)
      real x1r(kdata), x1m(kdata), s1r(kdata),s1m(kdata),xwork(kdata),
     1  y1r(kdata),y1m(kdata),x1(kdata),y1(kdata),ywork(kdata)
      INTEGER n1, n2, n3, n4, n5

       data x1r(1:nr)/240.0,241.0,242.0,243.0,244.0,245.0,246.0,
     1 247.0,248.0,249.0,250.0,251.0,252.0,253.0,
     2 254.0,255.0,256.0,257.0,258.0,259.0,260.0,
     3 261.0,262.0,263.0,264.0,265.0,266.0,267.0,
     4 268.0,269.0,270.0,271.0,272.0,273.0,274.0,
     5 275.0,276.0,277.0,278.0,279.0,280.0,281.0,
     6 282.0,283.0,284.0,285.0,286.0,287.0,288.0,
     7 289.0,290.0,291.0,292.0,293.0,294.0,295.0,
     8 296.0,297.0,298.0,299.0,300.0,301.0,302.0,
     9 303.0,304.0,305.0,306.0,307.0,308.0,309.0,
     a 310.0,311.0,312.0,313.0,314.0,315.0,316.0,
     b 317.0,318.0,319.0,320.0,321.0,322.0,323.0,
     c 324.0,325.0,326.0,327.0,328.0,329.0,330.0,
     d 331.0,332.0,333.0,334.0,335.0,336.0,337.0,
     e 338.0,339.0,340.0/
       data s1r(1:nr)/.640E-21,.560E-21,.105E-20,.115E-20,.820E-21,
     1 .103E-20,.980E-21,.135E-20,.191E-20,.282E-20,.205E-20,.170E-20,
     2 .288E-20,.255E-20,.255E-20,.360E-20,.509E-20,.339E-20,.226E-20,
     3 .504E-20,.505E-20,.549E-20,.520E-20,.933E-20,.823E-20,.430E-20,
     4 .495E-20,.124E-19,.111E-19,.878E-20,.936E-20,.179E-19,.123E-19,
     5 .645E-20,.656E-20,.223E-19,.242E-19,.140E-19,.105E-19,.255E-19,
     6 .208E-19,.148E-19,.881E-20,.107E-19,.449E-19,.359E-19,.196E-19,
     7 .130E-19,.336E-19,.284E-19,.130E-19,.175E-19,.832E-20,.373E-19,
     8 .654E-19,.395E-19,.233E-19,.151E-19,.404E-19,.287E-19,.871E-20,
     9 .172E-19,.106E-19,.320E-19,.690E-19,.491E-19,.463E-19,.210E-19,
     a .149E-19,.341E-19,.195E-19,.521E-20,.112E-19,.112E-19,.475E-19,
     b .525E-19,.290E-19,.537E-19,.298E-19,.918E-20,.126E-19,.153E-19,
     c .669E-20,.345E-20,.816E-20,.185E-19,.595E-19,.349E-19,.109E-19,
     d .335E-19,.332E-19,.107E-19,.289E-20,.215E-20,.171E-20,.143E-20,
     e .194E-20,.417E-20,.236E-19,.471E-19,.248E-19/
       data y1r(1:nr)/0.270,0.272,0.274,0.276,0.278,0.280,0.282,
     1 0.284,0.286,0.288,0.290,0.291,0.292,0.293,
     2 0.294,0.295,0.296,0.297,0.298,0.299,0.300,
     3 0.308,0.316,0.324,0.332,0.340,0.348,0.356,
     4 0.364,0.372,0.380,0.399,0.418,0.437,0.456,
     5 0.475,0.494,0.513,0.532,0.551,0.570,0.586,
     6 0.602,0.618,0.634,0.650,0.666,0.682,0.698,
     7 0.714,0.730,0.735,0.740,0.745,0.750,0.755,
     8 0.760,0.765,0.770,0.775,0.780,0.780,0.780,
     9 0.780,0.780,0.780,0.780,0.780,0.780,0.780,
     a 0.780,0.764,0.748,0.732,0.716,0.700,0.684,
     b 0.668,0.652,0.636,0.620,0.585,0.550,0.515,
     c 0.480,0.445,0.410,0.375,0.340,0.305,0.270,
     d 0.243,0.216,0.189,0.162,0.135,0.108,0.081,
     e 0.054,0.027,0.000/
       data x1m(1:nm)/240.0,241.0,242.0,243.0,244.0,245.0,246.0,
     1 247.0,248.0,249.0,250.0,251.0,252.0,253.0,
     2 254.0,255.0,256.0,257.0,258.0,259.0,260.0,
     3 261.0,262.0,263.0,264.0,265.0,266.0,267.0,
     4 268.0,269.0,270.0,271.0,272.0,273.0,274.0,
     5 275.0,276.0,277.0,278.0,279.0,280.0,281.0,
     6 282.0,283.0,284.0,285.0,286.0,287.0,288.0,
     7 289.0,290.0,291.0,292.0,293.0,294.0,295.0,
     8 296.0,297.0,298.0,299.0,300.0,301.0,302.0,
     9 303.0,304.0,305.0,306.0,307.0,308.0,309.0,
     a 310.0,311.0,312.0,313.0,314.0,315.0,316.0,
     b 317.0,318.0,319.0,320.0,321.0,322.0,323.0,
     c 324.0,325.0,326.0,327.0,328.0,329.0,330.0,
     d 331.0,332.0,333.0,334.0,335.0,336.0,337.0,
     e 338.0,339.0,340.0,341.0,342.0,343.0,344.0,
     f 345.0,346.0,347.0,348.0,349.0,350.0,351.0,
     g 352.0,353.0,354.0,355.0,356.0,357.0,358.0,
     h 359.0,360.0/
       data s1m(1:nm)/.640E-21,.560E-21,.105E-20,.115E-20,.820E-21,
     1 .103E-20,.980E-21,.135E-20,.191E-20,.282E-20,.205E-20,.170E-20,
     2 .288E-20,.255E-20,.255E-20,.360E-20,.509E-20,.339E-20,.226E-20,
     3 .504E-20,.505E-20,.549E-20,.520E-20,.933E-20,.823E-20,.430E-20,
     4 .495E-20,.124E-19,.111E-19,.878E-20,.936E-20,.179E-19,.123E-19,
     5 .645E-20,.656E-20,.223E-19,.242E-19,.140E-19,.105E-19,.255E-19,
     6 .208E-19,.148E-19,.881E-20,.107E-19,.449E-19,.359E-19,.196E-19,
     7 .130E-19,.336E-19,.284E-19,.130E-19,.175E-19,.832E-20,.373E-19,
     8 .654E-19,.395E-19,.233E-19,.151E-19,.404E-19,.287E-19,.871E-20,
     9 .172E-19,.106E-19,.320E-19,.690E-19,.491E-19,.463E-19,.210E-19,
     a .149E-19,.341E-19,.195E-19,.521E-20,.112E-19,.112E-19,.475E-19,
     b .525E-19,.290E-19,.537E-19,.298E-19,.918E-20,.126E-19,.153E-19,
     c .669E-20,.345E-20,.816E-20,.185E-19,.595E-19,.349E-19,.109E-19,
     d .335E-19,.332E-19,.107E-19,.289E-20,.215E-20,.171E-20,.143E-20,
     e .194E-20,.417E-20,.236E-19,.471E-19,.248E-19,.759E-20,.681E-20,
     f .195E-19,.114E-19,.323E-20,.113E-20,.660E-21,.122E-20,.320E-21,
     g .380E-21,.104E-20,.713E-20,.221E-19,.154E-19,.676E-20,.135E-20,
     h .360E-21,.570E-22,.580E-21,.820E-21/
       data y1m(1:nm)/0.490,0.490,0.490,0.490,0.490,0.490,0.490,
     1 0.490,0.490,0.490,0.490,0.490,0.490,0.490,
     2 0.490,0.490,0.490,0.490,0.490,0.490,0.490,
     3 0.484,0.478,0.472,0.466,0.460,0.454,0.448,
     4 0.442,0.436,0.430,0.419,0.408,0.397,0.386,
     5 0.375,0.364,0.353,0.342,0.331,0.320,0.312,
     6 0.304,0.296,0.288,0.280,0.272,0.264,0.256,
     7 0.248,0.240,0.237,0.234,0.231,0.228,0.225,
     8 0.222,0.219,0.216,0.213,0.210,0.211,0.212,
     9 0.213,0.214,0.215,0.216,0.217,0.218,0.219,
     a 0.220,0.236,0.252,0.268,0.284,0.300,0.316,
     b 0.332,0.348,0.364,0.380,0.408,0.436,0.464,
     c 0.492,0.520,0.548,0.576,0.604,0.632,0.660,
     d 0.650,0.640,0.630,0.620,0.610,0.600,0.590,
     e 0.580,0.570,0.560,0.525,0.490,0.455,0.420,
     f 0.385,0.350,0.315,0.280,0.245,0.210,0.192,
     g 0.174,0.156,0.138,0.120,0.102,0.084,0.066,
     h 0.048,0.000/

* local

      REAL yg(kw), yg1(kw), yg2(kw), yg3(kw), yg4(kw), yg5(kw)
      REAL a, b, c
      REAL a0, a1, a2, a3, a4, a5, a6, a7
      REAL b0, b1, b2, b3, b4
      REAL phi1, phi2, phi20, ak300, akt
      REAL qy, qy1, qy2, qy3

      REAL sigma, sig, slope
      REAL xs
      REAL t
      REAL dum
      INTEGER idum

      INTEGER i
      INTEGER irow, icol, irev
      INTEGER ierr

      INTEGER mopt1, mopt2


****************************************************************
**************** CH2O photodissociatation

      j = j+1
      jlabel(j) = 'HCHO -> 2HO2 + CO' 

* cross section 
      n = nr
      x1(1:n)=x1r(1:n) 
      ywork(1:n)=s1r(1:n)
      
      CALL addpnt(x1,ywork,kdata,n,x1(1)*(1.-deltax),0.)
      CALL addpnt(x1,ywork,kdata,n,               0.,0.)
      CALL addpnt(x1,ywork,kdata,n,x1(n)*(1.+deltax),0.)
      CALL addpnt(x1,ywork,kdata,n,           1.e+38,0.)
      CALL inter2(nw,wl,yg,n,x1,ywork,ierr)
      IF (ierr .NE. 0) THEN
         WRITE(*,*) ierr, jlabel(j)
         STOP
      ENDIF

* quantum yield 
      n=nr
      xwork(1:n)=x1r(1:n)
      ywork(1:n)=y1r(1:n)
       
      CALL addpnt(xwork,ywork,kdata,n,xwork(1)*(1.-deltax),0.)
      CALL addpnt(xwork,ywork,kdata,n,               0.,0.)
      CALL addpnt(xwork,ywork,kdata,n,xwork(n)*(1.+deltax),0.)
      CALL addpnt(xwork,ywork,kdata,n,           1.e+38,0.)
      CALL inter2(nw,wl,yg1,n,xwork,ywork,ierr)
      IF (ierr .NE. 0) THEN
         WRITE(*,*) ierr, jlabel(j)
         STOP
      ENDIF

      DO iw = 1, nw - 1
         DO i = 1, nz
            sq(j,i,iw) = yg(iw)*yg1(iw)
         ENDDO
      ENDDO

      j = j+1
      jlabel(j) = 'HCHO -> CO'

* cross section 
      n = nm
      x1(1:n)=x1m(1:n) 
      ywork(1:n)=s1m(1:n)
      
      CALL addpnt(x1,ywork,kdata,n,x1(1)*(1.-deltax),0.)
      CALL addpnt(x1,ywork,kdata,n,               0.,0.)
      CALL addpnt(x1,ywork,kdata,n,x1(n)*(1.+deltax),0.)
      CALL addpnt(x1,ywork,kdata,n,           1.e+38,0.)
      CALL inter2(nw,wl,yg,n,x1,ywork,ierr)
      IF (ierr .NE. 0) THEN
         WRITE(*,*) ierr, jlabel(j)
         STOP
      ENDIF

* quantum yield 
      n=nm
      xwork(1:n)=x1m(1:n)
      ywork(1:n)=y1m(1:n)
       
      CALL addpnt(xwork,ywork,kdata,n,xwork(1)*(1.-deltax),0.)
      CALL addpnt(xwork,ywork,kdata,n,               0.,0.)
      CALL addpnt(xwork,ywork,kdata,n,xwork(n)*(1.+deltax),0.)
      CALL addpnt(xwork,ywork,kdata,n,           1.e+38,0.)
      CALL inter2(nw,wl,yg1,n,xwork,ywork,ierr)
      IF (ierr .NE. 0) THEN
         WRITE(*,*) ierr, jlabel(j)
         STOP
      ENDIF

      DO iw = 1, nw - 1
         DO i = 1, nz
            sq(j,i,iw) = yg(iw)*yg1(iw)
         ENDDO
      ENDDO

      END


      SUBROUTINE rn_CCHO(nw,wl,wc,nz,tlev,airlev,j,sq,jlabel)

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Provide product (cross section) x (quantum yield) for CH3CHO photolysis: =*
*=           CH3CHO + hv -> CO + HO2 + C_O2                                  =*
*=  Cross section:  Choice between                                           =*
*=                   (1) IUPAC 97 data, from Martinez et al.                 =*
*=                   (2) Calvert and Pitts                                   =*
*=                   (3) Martinez et al., Table 1 scanned from paper         =*
*=                   (4) KFA tabulations                                     =*
*=  Quantum yields: Choice between                                           =*
*=                   (1) IUPAC 97, pressure correction using Horowith and    =*
*=                                 Calvert, 1982                             =*
*=                   (2) NCAR data file, from Moortgat, 1986                 =*
*-----------------------------------------------------------------------------*
*=  PARAMETERS:                                                              =*
*=  NW     - INTEGER, number of specified intervals + 1 in working        (I)=*
*=           wavelength grid                                                 =*
*=  WL     - REAL, vector of lower limits of wavelength intervals in      (I)=*
*=           working wavelength grid                                         =*
*=  WC     - REAL, vector of center points of wavelength intervals in     (I)=*
*=           working wavelength grid                                         =*
*=  NZ     - INTEGER, number of altitude levels in working altitude grid  (I)=*
*=  TLEV   - REAL, temperature (K) at each specified altitude level       (I)=*
*=  AIRLEV - REAL, air density (molec/cc) at each altitude level          (I)=*
*=  J      - INTEGER, counter for number of weighting functions defined  (IO)=*
*=  SQ     - REAL, cross section x quantum yield (cm^2) for each          (O)=*
*=           photolysis reaction defined, at each defined wavelength and     =*
*=           at each defined altitude level                                  =*
*=  JLABEL - CHARACTER*40, string identifier for each photolysis reaction (O)=*
*=           defined                                                         =*
*-----------------------------------------------------------------------------*
*=  EDIT HISTORY:                                                            =*
*=  05/98  Original, adapted from former JSPEC1 subroutine                   =*
*-----------------------------------------------------------------------------*
*= This program is free software;  you can redistribute it and/or modify     =*
*= it under the terms of the GNU General Public License as published by the  =*
*= Free Software Foundation;  either version 2 of the license, or (at your   =*
*= option) any later version.                                                =*
*= The TUV package is distributed in the hope that it will be useful, but    =*
*= WITHOUT ANY WARRANTY;  without even the implied warranty of MERCHANTIBI-  =*
*= LITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public     =*
*= License for more details.                                                 =*
*= To obtain a copy of the GNU General Public License, write to:             =*
*= Free Software Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.   =*
*-----------------------------------------------------------------------------*

      IMPLICIT NONE
      INCLUDE 'tuv.params'

* input

      INTEGER nw
      REAL wl(kw), wc(kw)
      
      INTEGER nz

      REAL tlev(kz)
      REAL airlev(kz)

* weighting functions

      CHARACTER*40 jlabel(kj)
      REAL sq(kj,kz,kw)

* input/output:

      INTEGER j

* data arrays

      INTEGER kdata,ns,ny,np
      PARAMETER(kdata=150,ns=106,ny=12,np=5)

      INTEGER i, n
      INTEGER n1, n2
      REAL xs1(kdata), ss1(kdata),xp(kdata),yp(kdata),xwork(kdata)
      REAL xy1(kdata), y1(kdata), y2(kdata),x2(kdata),ywork(kdata)
       data xs1(1:ns)/202.,206.,210.,214.,218.,222.,226.,
     1 230.,234.,238.,242.,246.,250.,254.,
     2 258.,262.,266.,270.,274.,278.,280.,
     3 281.,282.,283.,284.,285.,286.,287.,
     4 288.,289.,290.,291.,292.,293.,294.,
     5 295.,296.,297.,298.,299.,300.,301.,
     6 302.,303.,304.,305.,306.,307.,308.,
     7 309.,310.,311.,312.,313.,314.,315.,
     8 316.,317.,318.,319.,320.,321.,322.,
     9 323.,324.,325.,326.,327.,328.,329.,
     a 330.,331.,332.,333.,334.,335.,336.,
     b 337.,338.,339.,340.,341.,342.,343.,
     c 344.,345.,346.,347.,348.,349.,350.,
     d 351.,352.,353.,354.,355.,356.,357.,
     e 358.,359.,360.,361.,362.,363.,364.,
     f 365./
       data ss1(1:ns)/0.560E-21,0.530E-21,0.490E-21,0.480E-21,0.520E-21,
     1 0.650E-21,0.960E-21,0.151E-20,0.241E-20,0.375E-20,0.564E-20,
     2 0.818E-20,0.113E-19,0.152E-19,0.199E-19,0.244E-19,0.305E-19,
     3 0.342E-19,0.403E-19,0.419E-19,0.450E-19,0.469E-19,0.472E-19,
     4 0.475E-19,0.461E-19,0.449E-19,0.444E-19,0.459E-19,0.472E-19,
     5 0.477E-19,0.489E-19,0.478E-19,0.468E-19,0.453E-19,0.433E-19,
     6 0.427E-19,0.424E-19,0.438E-19,0.441E-19,0.426E-19,0.416E-19,
     7 0.399E-19,0.386E-19,0.372E-19,0.348E-19,0.342E-19,0.342E-19,
     8 0.336E-19,0.333E-19,0.314E-19,0.293E-19,0.276E-19,0.253E-19,
     9 0.247E-19,0.244E-19,0.220E-19,0.204E-19,0.207E-19,1.979E-20,
     a 1.874E-20,1.723E-20,1.484E-20,1.402E-20,1.244E-20,1.091E-20,
     b 1.136E-20,1.074E-20,0.858E-20,0.747E-20,0.707E-20,0.688E-20,
     c 0.588E-20,0.530E-20,0.398E-20,0.363E-20,0.350E-20,0.238E-20,
     d 0.222E-20,0.205E-20,0.219E-20,0.150E-20,0.740E-21,0.420E-21,
     e 0.310E-21,0.260E-21,0.210E-21,0.190E-21,0.150E-21,0.160E-21,
     f 0.100E-21,0.800E-22,0.700E-22,0.600E-22,0.500E-22,0.500E-22,
     g 0.400E-22,0.500E-22,0.300E-22,0.400E-22,0.200E-22,0.300E-22,
     h 0.200E-22,0.100E-22,0.000E+00,0.000E+00,0.000E+00/
       data xy1(1:ny)/260.,270.,280.,290.,295.,300.,305.,
     1  310.,315.,320.,325.,330./
       data y1(1:ny)/0.31,0.39,0.58,0.53,0.48,0.43,0.37,
     1  0.29,0.17,0.10,0.04,0.00/
       data y2(1:ny)/0.46,0.31,0.05,0.01,0.00,0.00,0.00,
     1  0.00,0.00,0.00,0.00,0.00/
       data xp(1:np)/290.,300.,313.,320.,331.2/
       data yp(1:np)/0.59,1.51,5.83,8.48,2.88/ 
       
* local

      REAL yg(kw), yg1(kw), yg2(kw), yg3(kw), yg4(kw)
      REAL qy1, qy2, qy3
      REAL sig
      REAL dum
      INTEGER ierr
      INTEGER  iz, iw

      INTEGER mabs, myld

****************************************************************
************************* CH3CHO photolysis

      j = j+1
      jlabel(j) = 'CH3CHO -> CO + HO2 + C_O2'

* options
* mabs for cross sections
* myld for quantum yields

* Absorption:
* 3:  Martinez et al., Table 1 scanned from paper


* Quantum yield
* 1:  DATAJ1/CH3CHO/CH3CHO_iup.yld
* pressure correction using Horowitz and Calvert 1982, based on slope/intercepth
* of Stern-Volmer plots

      n = ns      
      xwork(1:n)=xs1(1:n)
      ywork(1:n)=ss1(1:n)      

      CALL addpnt(xwork,ywork,kdata,n,xwork(1)*(1.-deltax),0.)
      CALL addpnt(xwork,ywork,kdata,n,               0.,0.)
      CALL addpnt(xwork,ywork,kdata,n,xwork(n)*(1.+deltax),0.)
      CALL addpnt(xwork,ywork,kdata,n,           1.e+38,0.)
      CALL inter2(nw,wl,yg,n,xwork,ywork,ierr)
      IF (ierr .NE. 0) THEN
        WRITE(*,*) ierr, jlabel(j)
         STOP
      ENDIF

* quantum yields


         n1 = ny
         n2 = ny
	 do i=1,ny
	  x2(i)=xy1(i)
	 enddo 

         xwork(1:n)=xy1(1:n)
	 ywork(1:n)=y1(1:n)
         CALL addpnt(xwork,ywork,kdata,n1,xwork(1)*(1.-deltax),0.)
         CALL addpnt(xwork,ywork,kdata,n1,               0.,0.)
         CALL addpnt(xwork,ywork,kdata,n1,xwork(n1)*(1.+deltax),0.)
         CALL addpnt(xwork,ywork,kdata,n1,           1.e+38,0.)
         CALL inter2(nw,wl,yg1,n1,xwork,ywork,ierr)
         IF (ierr .NE. 0) THEN
            WRITE(*,*) ierr, jlabel(j)
            STOP
         ENDIF

	 ywork(1:n)=y2(1:n)
         CALL addpnt(x2,ywork,kdata,n2,x2(1)*(1.-deltax),0.)
         CALL addpnt(x2,ywork,kdata,n2,               0.,0.)
         CALL addpnt(x2,ywork,kdata,n2,x2(n2)*(1.+deltax),0.)
         CALL addpnt(x2,ywork,kdata,n2,           1.e+38,0.)
         CALL inter2(nw,wl,yg2,n2,x2,ywork,ierr)
         IF (ierr .NE. 0) THEN
            WRITE(*,*) ierr, jlabel(j)
            STOP
         ENDIF

         DO iw = 1, nw-1
            yg3(iw) = 0.
         ENDDO

* pressure-dependence parameters
      
         n = np

         xwork(1:n)=xp(1:n)
	 ywork(1:n)=yp(1:n) 
         CALL addpnt(xwork,ywork,kdata,n,xwork(1)*(1.-deltax),0.)
         CALL addpnt(xwork,ywork,kdata,n,               0.,0.)
         CALL addpnt(xwork,ywork,kdata,n,xwork(n)*(1.+deltax),0.)
         CALL addpnt(xwork,ywork,kdata,n,           1.e+38,0.)
         CALL inter2(nw,wl,yg4,n,xwork,ywork,ierr)
         IF (ierr .NE. 0) THEN
            WRITE(*,*) ierr, jlabel(j)
            STOP
         ENDIF

* combine:

      DO iw = 1, nw - 1
         DO i = 1, nz

            sig = yg(iw)

* quantum yields:

            qy1 = yg1(iw)
c            qy2 = yg2(iw)
c            qy3 = yg3(iw)

* pressure correction for channel 1, CH3 + CHO
* based on Horowitz and Calvert 1982.

            qy1 = qy1 * (1. + yg4(iw))/(1. + yg4(iw)*airlev(i)/2.465E19)
            qy1 = MIN(1., qy1)
            qy1 = MAX(0., qy1)

c            sq(j-2,i,iw) = sig * qy1
c            sq(j-1,i,iw) = sig * qy2
c            sq(j  ,i,iw) = sig * qy3

            sq(j  ,i,iw) = sig * qy1
         ENDDO
      ENDDO

      END

      SUBROUTINE rn_RCHO(nw,wl,wc,nz,tlev,airlev,j,sq,jlabel)

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Provide product (cross section) x (quantum yield) for CH3CHO photolysis: =*
*=           CH3CHO + hv -> CCHO+RO2_R+CO+HO2                                =*
*=                    tyh modify for SAPRC99 07/01                           =*
*-----------------------------------------------------------------------------*
*=  PARAMETERS:                                                              =*
*=  NW     - INTEGER, number of specified intervals + 1 in working        (I)=*
*=           wavelength grid                                                 =*
*=  WL     - REAL, vector of lower limits of wavelength intervals in      (I)=*
*=           working wavelength grid                                         =*
*=  WC     - REAL, vector of center points of wavelength intervals in     (I)=*
*=           working wavelength grid                                         =*
*=  NZ     - INTEGER, number of altitude levels in working altitude grid  (I)=*
*=  TLEV   - REAL, temperature (K) at each specified altitude level       (I)=*
*=  AIRLEV - REAL, air density (molec/cc) at each altitude level          (I)=*
*=  J      - INTEGER, counter for number of weighting functions defined  (IO)=*
*=  SQ     - REAL, cross section x quantum yield (cm^2) for each          (O)=*
*=           photolysis reaction defined, at each defined wavelength and     =*
*=           at each defined altitude level                                  =*
*=  JLABEL - CHARACTER*40, string identifier for each photolysis reaction (O)=*
*=           defined                                                         =*
*-----------------------------------------------------------------------------*
*= This program is free software;  you can redistribute it and/or modify     =*
*= it under the terms of the GNU General Public License as published by the  =*
*= Free Software Foundation;  either version 2 of the license, or (at your   =*
*= option) any later version.                                                =*
*= The TUV package is distributed in the hope that it will be useful, but    =*
*= WITHOUT ANY WARRANTY;  without even the implied warranty of MERCHANTIBI-  =*
*= LITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public     =*
*= License for more details.                                                 =*
*= To obtain a copy of the GNU General Public License, write to:             =*
*= Free Software Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.   =*
*-----------------------------------------------------------------------------*

      IMPLICIT NONE
      INCLUDE 'tuv.params'

* input

      INTEGER nw
      REAL wl(kw), wc(kw)
      
      INTEGER nz

      REAL tlev(kz)
      REAL airlev(kz)

* weighting functions

      CHARACTER*40 jlabel(kj)
      REAL sq(kj,kz,kw)

* input/output:

      INTEGER j

* data arrays

      INTEGER kdata,ns
      PARAMETER(kdata=150,ns=50)

      INTEGER i, n
      REAL x1(kdata), s1(kdata), y1(kdata),xp(kdata), xwork(kdata), 
     1 ywork(kdata)

       data x1(1:ns)/294.0,295.0,296.0,297.0,298.0,299.0,300.0,
     1 301.0,302.0,303.0,304.0,305.0,306.0,307.0,
     2 308.0,309.0,310.0,311.0,312.0,313.0,314.0,
     3 315.0,316.0,317.0,318.0,319.0,320.0,321.0,
     4 322.0,323.0,324.0,325.0,326.0,327.0,328.0,
     5 329.0,330.0,331.0,332.0,333.0,334.0,335.0,
     6 336.0,337.0,338.0,339.0,340.0,341.0,342.0,
     7 343.0/
       data s1(1:ns)/.580E-19,.557E-19,.537E-19,.516E-19,.502E-19,
     1 .502E-19,.504E-19,.509E-19,.507E-19,.494E-19,.469E-19,.432E-19,
     2 .404E-19,.381E-19,.365E-19,.362E-19,.360E-19,.353E-19,.350E-19,
     3 .332E-19,.306E-19,.277E-19,.243E-19,.218E-19,.200E-19,.186E-19,
     4 .183E-19,.178E-19,.166E-19,.158E-19,.149E-19,.130E-19,.113E-19,
     5 .996E-20,.828E-20,.685E-20,.575E-20,.494E-20,.466E-20,.430E-20,
     6 .373E-20,.325E-20,.280E-20,.230E-20,.185E-20,.166E-20,.155E-20,
     7 .119E-20,.760E-21,.450E-21/
       data y1(1:ns)/0.890,0.885,0.880,0.875,0.870,0.865,0.860,
     1 0.855,0.850,0.818,0.786,0.755,0.723,0.691,
     2 0.659,0.627,0.596,0.564,0.532,0.500,0.480,
     3 0.460,0.440,0.420,0.400,0.380,0.360,0.340,
     4 0.320,0.300,0.280,0.260,0.248,0.236,0.223,
     5 0.211,0.199,0.187,0.174,0.162,0.150,0.133,
     6 0.117,0.100,0.083,0.067,0.050,0.033,0.017,
     7 0.000/

* local

      REAL yg(kw), yg1(kw), yg2(kw), yg3(kw), yg4(kw)
      REAL qy1, qy2, qy3
      REAL sig
      REAL dum
      INTEGER ierr
      INTEGER  iz, iw


************************* C2CHO photolysis

      j = j+1
      jlabel(j) = 'C2CHO -> CCHO + RO2_R + CO + HO2'

* cross section 
      n = ns
      xp(1:n)=x1(1:n) 
      ywork(1:n)=s1(1:n)
      
      CALL addpnt(xp,ywork,kdata,n,xp(1)*(1.-deltax),0.)
      CALL addpnt(xp,ywork,kdata,n,               0.,0.)
      CALL addpnt(xp,ywork,kdata,n,xp(n)*(1.+deltax),0.)
      CALL addpnt(xp,ywork,kdata,n,           1.e+38,0.)
      CALL inter2(nw,wl,yg,n,xp,ywork,ierr)
      IF (ierr .NE. 0) THEN
         WRITE(*,*) ierr, jlabel(j)
         STOP
      ENDIF

* quantum yield 
      n=ns
      xwork(1:n)=x1(1:n)
      ywork(1:n)=y1(1:n)
       
      CALL addpnt(xwork,ywork,kdata,n,xwork(1)*(1.-deltax),0.)
      CALL addpnt(xwork,ywork,kdata,n,               0.,0.)
      CALL addpnt(xwork,ywork,kdata,n,xwork(n)*(1.+deltax),0.)
      CALL addpnt(xwork,ywork,kdata,n,           1.e+38,0.)
      CALL inter2(nw,wl,yg1,n,xwork,ywork,ierr)
      IF (ierr .NE. 0) THEN
         WRITE(*,*) ierr, jlabel(j)
         STOP
      ENDIF

      DO iw = 1, nw - 1
         DO i = 1, nz
            sq(j,i,iw) = yg(iw)*yg1(iw)
         ENDDO
      ENDDO

      END

      SUBROUTINE rn_GLY(nw,wl,wc,nz,tlev,airlev,j,sq,jlabel)

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Provide the product (cross section) x (quantum yield) for CHOCHO         =*
*=  photolysis:                                                              =*
*=            ( CHOCHO + hv -> Products )                                    =*
*=              GLY + hv = CO + CO + HO2 + HO2                               =*
*=              GLY + hv = HCHO + CO                                         =*
*=    tyh created for SAPRC99 07/01                                          =*
*=                  (1) Plum et al., as tabulated by IUPAC 97                =*
*=                  (2) Plum et al., as tabulated by KFA.                    =*
*=  Quantum yield: IUPAC 97 recommendation                                   =*
*-----------------------------------------------------------------------------*
*=  PARAMETERS:                                                              =*
*=  NW     - INTEGER, number of specified intervals + 1 in working        (I)=*
*=           wavelength grid                                                 =*
*=  WL     - REAL, vector of lower limits of wavelength intervals in      (I)=*
*=           working wavelength grid                                         =*
*=  WC     - REAL, vector of center points of wavelength intervals in     (I)=*
*=           working wavelength grid                                         =*
*=  NZ     - INTEGER, number of altitude levels in working altitude grid  (I)=*
*=  TLEV   - REAL, temperature (K) at each specified altitude level       (I)=*
*=  AIRLEV - REAL, air density (molec/cc) at each altitude level          (I)=*
*=  J      - INTEGER, counter for number of weighting functions defined  (IO)=*
*=  SQ     - REAL, cross section x quantum yield (cm^2) for each          (O)=*
*=           photolysis reaction defined, at each defined wavelength and     =*
*=           at each defined altitude level                                  =*
*=  JLABEL - CHARACTER*40, string identifier for each photolysis reaction (O)=*
*=           defined                                                         =*
*-----------------------------------------------------------------------------*
*=  EDIT HISTORY:                                                            =*
*=  05/98  Original, adapted from former JSPEC1 subroutine                   =*
*-----------------------------------------------------------------------------*
*= This program is free software;  you can redistribute it and/or modify     =*
*= it under the terms of the GNU General Public License as published by the  =*
*= Free Software Foundation;  either version 2 of the license, or (at your   =*
*= option) any later version.                                                =*
*= The TUV package is distributed in the hope that it will be useful, but    =*
*= WITHOUT ANY WARRANTY;  without even the implied warranty of MERCHANTIBI-  =*
*= LITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public     =*
*= License for more details.                                                 =*
*= To obtain a copy of the GNU General Public License, write to:             =*
*= Free Software Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.   =*
*-----------------------------------------------------------------------------*

      IMPLICIT NONE
      INCLUDE 'tuv.params'

* input

      INTEGER nw
      REAL wl(kw), wc(kw)
      
      INTEGER nz

      REAL tlev(kz)
      REAL airlev(kz)

* weighting functions

      CHARACTER*40 jlabel(kj)
      REAL sq(kj,kz,kw)

* input/output:

      INTEGER j

* data arrays

      INTEGER kdata,nr,na
      PARAMETER(kdata=150,nr=70,na=122)

      INTEGER i, n
      REAL x1(kdata),x1r(kdata),x1a(kdata),s1r(kdata),s1a(kdata)
      REAL y1(kdata),xwork(kdata),ywork(kdata)

       data x1r(1:nr)/230.0,235.0,240.0,245.0,250.0,255.0,260.0,
     1 265.0,270.0,275.0,280.0,285.0,290.0,295.0,
     2 300.0,305.0,310.0,312.5,315.0,320.0,325.0,
     3 327.5,330.0,335.0,340.0,345.0,350.0,355.0,
     4 360.0,365.0,370.0,375.0,380.0,382.0,384.0,
     5 386.0,388.0,390.0,391.0,392.0,393.0,394.0,
     6 395.0,396.0,397.0,398.0,399.0,400.0,401.0,
     7 402.0,403.0,404.0,405.0,406.0,407.0,408.0,
     8 409.0,410.0,411.0,411.5,412.0,413.0,413.5,
     9 414.0,414.5,415.0,415.5,416.0,417.0,418.0/
       data s1r(1:nr)/.287E-20,.287E-20,.430E-20,.573E-20,.860E-20,
     1 .115E-19,.143E-19,.186E-19,.229E-19,.258E-19,.287E-19,.330E-19,
     2 .315E-19,.330E-19,.358E-19,.272E-19,.272E-19,.287E-19,.229E-19,
     3 .143E-19,.115E-19,.143E-19,.115E-19,.287E-20,.000E+00,.000E+00,
     4 .000E+00,.000E+00,.229E-20,.287E-20,.803E-20,.100E-19,.172E-19,
     5 .158E-19,.149E-19,.149E-19,.287E-19,.315E-19,.324E-19,.304E-19,
     6 .223E-19,.263E-19,.304E-19,.263E-19,.243E-19,.324E-19,.304E-19,
     7 .284E-19,.324E-19,.446E-19,.527E-19,.426E-19,.304E-19,.304E-19,
     8 .284E-19,.243E-19,.284E-19,.608E-19,.507E-19,.608E-19,.486E-19,
     9 .831E-19,.648E-19,.750E-19,.811E-19,.811E-19,.689E-19,.426E-19,
     a .486E-19,.588E-19/
       data x1a(1:na)/230.0,235.0,240.0,245.0,250.0,255.0,260.0,
     1 265.0,270.0,275.0,280.0,285.0,290.0,295.0,
     2 300.0,305.0,310.0,312.5,315.0,320.0,325.0,
     3 327.5,330.0,335.0,340.0,355.0,360.0,365.0,
     4 370.0,375.0,380.0,382.0,384.0,386.0,388.0,
     5 390.0,391.0,392.0,393.0,394.0,395.0,396.0,
     6 397.0,398.0,399.0,400.0,401.0,402.0,403.0,
     7 404.0,405.0,406.0,407.0,408.0,409.0,410.0,
     8 411.0,411.5,412.0,413.0,413.5,414.0,414.5,
     9 415.0,415.5,416.0,417.0,418.0,419.0,420.0,
     a 421.0,421.5,422.0,422.5,423.0,424.0,425.0,
     b 426.0,426.5,427.0,428.0,429.0,430.0,431.0,
     c 432.0,433.0,434.0,434.5,435.0,436.0,436.5,
     d 437.0,438.0,438.5,439.0,440.0,441.0,442.0,
     e 443.0,444.0,445.0,446.0,447.0,448.0,449.0,
     f 450.0,451.0,451.5,452.0,453.0,454.0,455.0,
     g 455.5,456.0,457.0,458.0,458.5,459.0,460.0,
     h 460.5,461.0,462.0/
       data s1a(1:na)/.287E-20,.287E-20,.430E-20,.573E-20,.860E-20,
     1 .115E-19,.143E-19,.186E-19,.229E-19,.258E-19,.287E-19,.330E-19,
     2 .315E-19,.330E-19,.358E-19,.272E-19,.272E-19,.287E-19,.229E-19,
     3 .143E-19,.115E-19,.143E-19,.115E-19,.287E-20,.000E+00,.000E+00,
     4 .229E-20,.287E-20,.803E-20,.100E-19,.172E-19,.158E-19,.149E-19,
     5 .149E-19,.287E-19,.315E-19,.324E-19,.304E-19,.223E-19,.263E-19,
     6 .304E-19,.263E-19,.243E-19,.324E-19,.304E-19,.284E-19,.324E-19,
     7 .446E-19,.527E-19,.426E-19,.304E-19,.304E-19,.284E-19,.243E-19,
     8 .284E-19,.608E-19,.507E-19,.608E-19,.486E-19,.831E-19,.648E-19,
     9 .750E-19,.811E-19,.811E-19,.689E-19,.426E-19,.486E-19,.588E-19,
     a .669E-19,.385E-19,.567E-19,.446E-19,.527E-19,.105E-18,.851E-19,
     b .608E-19,.729E-19,.118E-18,.130E-18,.107E-18,.166E-18,.405E-19,
     c .507E-19,.486E-19,.405E-19,.365E-19,.405E-19,.608E-19,.507E-19,
     d .811E-19,.113E-18,.527E-19,.101E-18,.138E-18,.770E-19,.247E-18,
     e .811E-19,.608E-19,.750E-19,.932E-19,.113E-18,.527E-19,.243E-19,
     f .284E-19,.385E-19,.608E-19,.109E-18,.932E-19,.122E-18,.239E-18,
     g .170E-18,.340E-18,.405E-18,.101E-18,.162E-19,.122E-19,.142E-19,
     h .405E-20,.405E-20,.608E-20,.203E-20,.000E+00/

* local

      REAL yg(kw),yg1(kw)
      REAL qy
      REAL sig
      INTEGER ierr
      INTEGER iw

      INTEGER mabs, myld

************************* CHOCHO photolysis
* 1:  CHOCHO

      j = j+1
      jlabel(j) = 'GLY ->  2CO + 2HO2'

* Absorption:
* cross section 
      n = nr
      x1(1:n)=x1r(1:n) 
      ywork(1:n)=s1r(1:n)
      
      CALL addpnt(x1,ywork,kdata,n,x1(1)*(1.-deltax),0.)
      CALL addpnt(x1,ywork,kdata,n,               0.,0.)
      CALL addpnt(x1,ywork,kdata,n,x1(n)*(1.+deltax),0.)
      CALL addpnt(x1,ywork,kdata,n,           1.e+38,0.)
      CALL inter2(nw,wl,yg,n,x1,ywork,ierr)
      IF (ierr .NE. 0) THEN
         WRITE(*,*) ierr, jlabel(j)
         STOP
      ENDIF
c      print*,'yg=',yg

* quantum yield 

      DO iw = 1, nw - 1
       if(wc(iw) .lt. 325.) then
          qy = 0.4
       else
          qy = 0.029
       endif
         DO i = 1, nz
            sq(j,i,iw) = yg(iw)*qy
         ENDDO
      ENDDO

      j = j+1
      jlabel(j) = 'GLY -> HCHO + CO'

* cross section 
      n = na
      x1(1:n)=x1a(1:n) 
      ywork(1:n)=s1a(1:n)
      
      CALL addpnt(x1,ywork,kdata,n,x1(1)*(1.-deltax),0.)
      CALL addpnt(x1,ywork,kdata,n,               0.,0.)
      CALL addpnt(x1,ywork,kdata,n,x1(n)*(1.+deltax),0.)
      CALL addpnt(x1,ywork,kdata,n,           1.e+38,0.)
      CALL inter2(nw,wl,yg,n,x1,ywork,ierr)
c      print*,'nw,wl=',nw,wl
c      print*,'x1=',x1
c      print*,'s1a=',s1a
      IF (ierr .NE. 0) THEN
         WRITE(*,*) ierr, jlabel(j)
         STOP
      ENDIF
c      print*,'yg=',yg
      
* quantum yield is 0.006 for all wavelength
      qy=0.006
      DO iw = 1, nw - 1
         DO i = 1, nz
            sq(j,i,iw) = yg(iw)*qy
         ENDDO
      ENDDO
c      open(111,file='tmp.111')
c      write(111,*)'sq=',(sq(j,1,iw),iw=1,nw-1)
       
      END

      SUBROUTINE rn_MGLY(nw,wl,wc,nz,tlev,airlev,j,sq,jlabel)

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Provide product of (cross section) x (quantum yield) for CH2O photolysis =*
*=            MGLY + hv -> HO2 + CO + CCO_O2                                 =*
*=                 tyh create for SAPRC99   07/01                            =*
*-----------------------------------------------------------------------------*
*=  PARAMETERS:                                                              =*
*=  NW     - INTEGER, number of specified intervals + 1 in working        (I)=*
*=           wavelength grid                                                 =*
*=  WL     - REAL, vector of lower limits of wavelength intervals in      (I)=*
*=           working wavelength grid                                         =*
*=  WC     - REAL, vector of center points of wavelength intervals in     (I)=*
*=           working wavelength grid                                         =*
*=  NZ     - INTEGER, number of altitude levels in working altitude grid  (I)=*
*=  TLEV   - REAL, temperature (K) at each specified altitude level       (I)=*
*=  AIRLEV - REAL, air density (molec/cc) at each altitude level          (I)=*
*=  J      - INTEGER, counter for number of weighting functions defined  (IO)=*
*=  SQ     - REAL, cross section x quantum yield (cm^2) for each          (O)=*
*=           photolysis reaction defined, at each defined wavelength and     =*
*=           at each defined altitude level                                  =*
*=  JLABEL - CHARACTER*40, string identifier for each photolysis reaction (O)=*
*=           defined                                                         =*
*-----------------------------------------------------------------------------*
*=  EDIT HISTORY:                                                            =*
*=  05/98  Original, adapted from former JSPEC1 subroutine                   =*
*-----------------------------------------------------------------------------*
*= This program is free software;  you can redistribute it and/or modify     =*
*= it under the terms of the GNU General Public License as published by the  =*
*= Free Software Foundation;  either version 2 of the license, or (at your   =*
*= option) any later version.                                                =*
*= The TUV package is distributed in the hope that it will be useful, but    =*
*= WITHOUT ANY WARRANTY;  without even the implied warranty of MERCHANTIBI-  =*
*= LITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public     =*
*= License for more details.                                                 =*
*= To obtain a copy of the GNU General Public License, write to:             =*
*= Free Software Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.   =*
*-----------------------------------------------------------------------------*

      IMPLICIT NONE
      INCLUDE 'tuv.params'

      INTEGER kdata
      PARAMETER(kdata=420)

* input

      INTEGER nw
      REAL wl(kw), wc(kw)
      
      INTEGER nz

      REAL tlev(kz)
      REAL airlev(kz)

* weighting functions

      CHARACTER*40 jlabel(kj)
      REAL sq(kj,kz,kw)

* input/output:

      INTEGER j, iz, iw

* data arrays

      INTEGER n,ns
      parameter(ns=405)
      real xs1(kdata), x1(kdata), s1(kdata),y1(kdata),ywork(kdata)

       data xs1(1:ns)/219.0,219.5,220.0,220.5,221.0,221.5,222.0,
     1  222.5,223.0,223.5,224.0,224.5,225.0,225.5,
     2  226.0,226.5,227.0,227.5,228.0,228.5,229.0,
     3  229.5,230.0,230.5,231.0,231.5,232.0,232.5,
     4  233.0,233.5,234.0,234.5,235.0,235.5,236.0,
     5  236.5,237.0,237.5,238.0,238.5,239.0,239.5,
     6  240.0,240.5,241.0,241.5,242.0,242.5,243.0,
     7  243.5,244.0,244.5,245.0,245.5,246.0,246.5,
     8  247.0,247.5,248.0,248.5,249.0,249.5,250.0,
     9  250.5,251.0,251.5,252.0,252.5,253.0,253.5,
     a  254.0,254.5,255.0,255.5,256.0,256.5,257.0,
     b  257.5,258.0,258.5,259.0,259.5,260.0,260.5,
     c  261.0,261.5,262.0,262.5,263.0,263.5,264.0,
     d  264.5,265.0,265.5,266.0,266.5,267.0,267.5,
     e  268.0,268.5,269.0,269.5,270.0,270.5,271.0,
     f  271.5,272.0,272.5,273.0,273.5,274.0,274.5,
     g  275.0,275.5,276.0,276.5,277.0,277.5,278.0,
     h  278.5,279.0,279.5,280.0,280.5,281.0,281.5,
     i  282.0,282.5,283.0,283.5,284.0,284.5,285.0,
     j  285.5,286.0,286.5,287.0,287.5,288.0,288.5,
     k  289.0,289.5,290.0,290.5,291.0,291.5,292.0,
     l  292.5,293.0,293.5,294.0,294.5,295.0,295.5,
     m  296.0,296.5,297.0,297.5,298.0,298.5,299.0,
     n  299.5,300.0,300.5,301.0,301.5,302.0,302.5,
     o  303.0,303.5,304.0,304.5,305.0,305.5,306.0,
     p  306.5,307.0,307.5,308.0,308.5,309.0,309.5,
     q  310.0,310.5,311.0,311.5,312.0,312.5,313.0,
     r  313.5,314.0,314.5,315.0,315.5,316.0,316.5,
     s  317.0,317.5,318.0,318.5,319.0,319.5,320.0,
     t  320.5,321.0,321.5,322.0,322.5,323.0,323.5,
     u  324.0,324.5,325.0,325.5,326.0,326.5,327.0,
     v  327.5,328.0,328.5,329.0,329.5,330.0,330.5,
     x 331.0,331.5,332.0,332.5,333.0,333.5,334.0,
     y 334.5,335.0,335.5,336.0,336.5,337.0,337.5,
     z 338.0,338.5,339.0,339.5,340.0,340.5,341.0,
     A 341.5,342.0,342.5,343.0,343.5,344.0,344.5,
     B 345.0,345.5,346.0,346.5,347.0,347.5,348.0,
     C 348.5,349.0,349.5,350.0,350.5,351.0,351.5,
     D 352.0,352.5,353.0,353.5,354.0,354.5,355.0,
     E 355.5,356.0,356.5,357.0,357.5,358.0,358.5,
     F 359.0,359.5,360.0,360.5,361.0,361.5,362.0,
     G 362.5,363.0,363.5,364.0,364.5,365.0,365.5,
     H 366.0,366.5,367.0,367.5,368.0,368.5,369.0,
     I 369.5,370.0,370.5,371.0,371.5,372.0,372.5,
     J 373.0,373.5,374.0,374.5,375.0,375.5,376.0,
     K 376.5,377.0,377.5,378.0,378.5,379.0,379.5,
     L 380.0,380.5,381.0,381.5,382.0,382.5,383.0,
     M 383.5,384.0,384.5,385.0,385.5,386.0,386.5,
     N 387.0,387.5,388.0,388.5,389.0,389.5,390.0,
     O 390.5,391.0,391.5,392.0,392.5,393.0,393.5,
     P 394.0,394.5,395.0,395.5,396.0,396.5,397.0,
     Q 397.5,398.0,398.5,399.0,399.5,400.0,400.5,
     R 401.0,401.5,402.0,402.5,403.0,403.5,404.0,
     S 404.5,405.0,405.5,406.0,406.5,407.0,407.5,
     T 408.0,408.5,409.0,409.5,410.0,410.5,411.0,
     U 411.5,412.0,412.5,413.0,413.5,414.0,414.5,
     V 415.0,415.5,416.0,416.5,417.0,417.5,418.0,
     W 418.5,419.0,419.5,420.0,420.5,421.0/
       data s1(1:ns)/.984E-20,.104E-19,.106E-19,.111E-19,.115E-19,
     1 .118E-19,.122E-19,.124E-19,.126E-19,.126E-19,.125E-19,.124E-19,
     2 .125E-19,.127E-19,.127E-19,.129E-19,.131E-19,.132E-19,.135E-19,
     3 .137E-19,.140E-19,.142E-19,.148E-19,.153E-19,.157E-19,.159E-19,
     4 .161E-19,.162E-19,.161E-19,.168E-19,.174E-19,.180E-19,.184E-19,
     5 .187E-19,.189E-19,.191E-19,.193E-19,.194E-19,.196E-19,.196E-19,
     6 .201E-19,.204E-19,.208E-19,.210E-19,.214E-19,.216E-19,.219E-19,
     7 .220E-19,.223E-19,.226E-19,.228E-19,.229E-19,.230E-19,.232E-19,
     8 .233E-19,.235E-19,.238E-19,.241E-19,.246E-19,.251E-19,.257E-19,
     9 .261E-19,.265E-19,.267E-19,.269E-19,.269E-19,.271E-19,.272E-19,
     a .273E-19,.274E-19,.276E-19,.278E-19,.282E-19,.287E-19,.293E-19,
     b .298E-19,.307E-19,.312E-19,.317E-19,.321E-19,.326E-19,.328E-19,
     c .329E-19,.331E-19,.333E-19,.334E-19,.336E-19,.338E-19,.342E-19,
     d .344E-19,.348E-19,.354E-19,.359E-19,.365E-19,.373E-19,.380E-19,
     e .387E-19,.395E-19,.402E-19,.408E-19,.413E-19,.417E-19,.420E-19,
     f .422E-19,.422E-19,.422E-19,.423E-19,.424E-19,.427E-19,.429E-19,
     g .431E-19,.433E-19,.437E-19,.442E-19,.448E-19,.456E-19,.464E-19,
     h .471E-19,.478E-19,.483E-19,.487E-19,.490E-19,.492E-19,.493E-19,
     i .494E-19,.492E-19,.490E-19,.486E-19,.483E-19,.479E-19,.476E-19,
     j .472E-19,.470E-19,.468E-19,.466E-19,.465E-19,.465E-19,.468E-19,
     k .473E-19,.478E-19,.484E-19,.489E-19,.492E-19,.492E-19,.490E-19,
     l .486E-19,.481E-19,.475E-19,.470E-19,.465E-19,.458E-19,.448E-19,
     m .438E-19,.427E-19,.417E-19,.407E-19,.399E-19,.394E-19,.388E-19,
     n .382E-19,.376E-19,.372E-19,.369E-19,.368E-19,.370E-19,.372E-19,
     o .374E-19,.374E-19,.375E-19,.371E-19,.362E-19,.351E-19,.338E-19,
     p .325E-19,.315E-19,.304E-19,.292E-19,.280E-19,.271E-19,.263E-19,
     q .252E-19,.243E-19,.234E-19,.225E-19,.219E-19,.212E-19,.206E-19,
     r .202E-19,.196E-19,.192E-19,.191E-19,.188E-19,.186E-19,.185E-19,
     s .186E-19,.187E-19,.187E-19,.187E-19,.183E-19,.175E-19,.169E-19,
     t .160E-19,.150E-19,.141E-19,.134E-19,.127E-19,.121E-19,.118E-19,
     u .114E-19,.108E-19,.101E-19,.962E-20,.928E-20,.875E-20,.849E-20,
     v .821E-20,.771E-20,.738E-20,.718E-20,.686E-20,.671E-20,.663E-20,
     x .646E-20,.629E-20,.621E-20,.618E-20,.620E-20,.549E-20,.521E-20,
     y .538E-20,.535E-20,.504E-20,.494E-20,.490E-20,.452E-20,.426E-20,
     z .411E-20,.376E-20,.361E-20,.358E-20,.347E-20,.332E-20,.322E-20,
     A .310E-20,.300E-20,.294E-20,.289E-20,.286E-20,.288E-20,.288E-20,
     B .289E-20,.291E-20,.295E-20,.300E-20,.308E-20,.318E-20,.325E-20,
     C .330E-20,.339E-20,.351E-20,.363E-20,.373E-20,.385E-20,.399E-20,
     D .427E-20,.447E-20,.463E-20,.478E-20,.492E-20,.507E-20,.523E-20,
     E .539E-20,.556E-20,.577E-20,.597E-20,.615E-20,.635E-20,.656E-20,
     F .676E-20,.695E-20,.720E-20,.744E-20,.764E-20,.789E-20,.815E-20,
     G .843E-20,.871E-20,.902E-20,.933E-20,.965E-20,.100E-19,.104E-19,
     H .108E-19,.111E-19,.115E-19,.119E-19,.123E-19,.127E-19,.131E-19,
     I .135E-19,.140E-19,.144E-19,.147E-19,.151E-19,.155E-19,.159E-19,
     J .164E-19,.170E-19,.173E-19,.177E-19,.181E-19,.186E-19,.190E-19,
     K .196E-19,.202E-19,.206E-19,.210E-19,.214E-19,.218E-19,.224E-19,
     L .230E-19,.237E-19,.242E-19,.247E-19,.254E-19,.262E-19,.269E-19,
     M .279E-19,.288E-19,.296E-19,.302E-19,.310E-19,.320E-19,.329E-19,
     N .339E-19,.351E-19,.362E-19,.369E-19,.370E-19,.377E-19,.388E-19,
     O .397E-19,.403E-19,.412E-19,.422E-19,.429E-19,.430E-19,.438E-19,
     P .447E-19,.455E-19,.456E-19,.459E-19,.467E-19,.480E-19,.487E-19,
     Q .496E-19,.508E-19,.519E-19,.523E-19,.539E-19,.546E-19,.554E-19,
     R .559E-19,.577E-19,.591E-19,.599E-19,.606E-19,.620E-19,.635E-19,
     S .652E-19,.654E-19,.664E-19,.693E-19,.715E-19,.719E-19,.732E-19,
     T .758E-19,.788E-19,.797E-19,.791E-19,.811E-19,.841E-19,.853E-19,
     U .859E-19,.860E-19,.880E-19,.904E-19,.945E-19,.934E-19,.937E-19,
     V .963E-19,.971E-19,.970E-19,.965E-19,.969E-19,.989E-19,.100E-18,
     W .102E-18,.100E-18,.102E-18,.101E-18,.101E-18,.103E-18,.101E-18,
     X .104E-18/
       data y1(1:ns)/1.000,1.000,1.000,1.000,1.000,1.000,1.000,
     1 1.000,1.000,1.000,1.000,1.000,1.000,1.000,
     2 1.000,1.000,1.000,1.000,1.000,1.000,1.000,
     3 1.000,1.000,1.000,1.000,1.000,1.000,1.000,
     4 1.000,1.000,1.000,1.000,1.000,1.000,1.000,
     5 1.000,1.000,1.000,1.000,1.000,1.000,1.000,
     6 1.000,1.000,1.000,1.000,1.000,1.000,1.000,
     7 1.000,1.000,1.000,1.000,1.000,1.000,1.000,
     8 1.000,1.000,1.000,1.000,1.000,1.000,1.000,
     9 1.000,1.000,1.000,1.000,1.000,1.000,1.000,
     a 1.000,1.000,1.000,1.000,1.000,1.000,1.000,
     b 1.000,1.000,1.000,1.000,1.000,1.000,1.000,
     c 1.000,1.000,1.000,1.000,1.000,1.000,1.000,
     d 1.000,1.000,1.000,1.000,1.000,1.000,1.000,
     e 1.000,1.000,1.000,1.000,1.000,1.000,1.000,
     f 1.000,1.000,1.000,1.000,1.000,1.000,1.000,
     g 1.000,1.000,1.000,1.000,1.000,1.000,1.000,
     h 1.000,1.000,1.000,1.000,1.000,1.000,1.000,
     i 1.000,1.000,1.000,1.000,1.000,1.000,1.000,
     j 1.000,1.000,1.000,1.000,1.000,1.000,1.000,
     k 1.000,1.000,1.000,1.000,1.000,1.000,1.000,
     l 1.000,1.000,1.000,1.000,1.000,1.000,1.000,
     m 1.000,1.000,1.000,1.000,1.000,1.000,1.000,
     n 1.000,1.000,1.000,1.000,1.000,1.000,1.000,
     o 1.000,1.000,1.000,1.000,1.000,1.000,1.000,
     p 1.000,1.000,1.000,1.000,1.000,1.000,1.000,
     q 1.000,1.000,1.000,1.000,1.000,1.000,1.000,
     r 1.000,1.000,1.000,1.000,1.000,1.000,1.000,
     s 1.000,1.000,1.000,1.000,1.000,1.000,1.000,
     t 1.000,1.000,1.000,1.000,1.000,1.000,1.000,
     u 1.000,1.000,1.000,1.000,1.000,1.000,1.000,
     v 1.000,1.000,1.000,1.000,1.000,1.000,1.000,
     x 1.000,1.000,1.000,1.000,1.000,1.000,1.000,
     y 1.000,1.000,1.000,1.000,1.000,1.000,1.000,
     z 1.000,1.000,1.000,1.000,1.000,1.000,1.000,
     A 1.000,1.000,1.000,1.000,1.000,0.992,0.984,
     B 0.976,0.968,0.960,0.953,0.945,0.937,0.929,
     C 0.921,0.913,0.905,0.897,0.889,0.881,0.873,
     D 0.865,0.858,0.850,0.842,0.834,0.826,0.818,
     E 0.810,0.802,0.794,0.786,0.778,0.770,0.763,
     F 0.755,0.747,0.739,0.731,0.723,0.715,0.707,
     G 0.699,0.691,0.683,0.675,0.668,0.660,0.652,
     H 0.644,0.636,0.628,0.620,0.612,0.604,0.596,
     I 0.588,0.580,0.573,0.565,0.557,0.549,0.541,
     J 0.533,0.525,0.517,0.509,0.501,0.493,0.486,
     K 0.478,0.470,0.462,0.454,0.446,0.438,0.430,
     L 0.422,0.414,0.406,0.398,0.391,0.383,0.375,
     M 0.367,0.359,0.351,0.343,0.335,0.327,0.319,
     N 0.311,0.303,0.296,0.288,0.280,0.272,0.264,
     O 0.256,0.248,0.240,0.232,0.224,0.216,0.208,
     P 0.201,0.193,0.185,0.177,0.169,0.161,0.153,
     Q 0.145,0.137,0.129,0.121,0.113,0.106,0.098,
     R 0.090,0.082,0.074,0.066,0.058,0.050,0.042,
     S 0.034,0.026,0.018,0.011,0.003,0.000,0.000,
     T 0.000,0.000,0.000,0.000,0.000,0.000,0.000,
     U 0.000,0.000,0.000,0.000,0.000,0.000,0.000,
     V 0.000,0.000,0.000,0.000,0.000,0.000,0.000,
     W 0.000,0.000,0.000,0.000,0.000,0.000/

* local

      REAL yg(kw), yg1(kw)
      REAL a, b, c
      REAL a0, a1, a2, a3, a4, a5, a6, a7
      REAL b0, b1, b2, b3, b4
      REAL phi1, phi2, phi20, ak300, akt
      REAL qy, qy1, qy2, qy3

      REAL sigma, sig, slope
      REAL xs
      REAL t
      REAL dum
      INTEGER idum

      INTEGER i
      INTEGER irow, icol, irev
      INTEGER ierr

      INTEGER mopt1, mopt2


****************************************************************
**************** CH2O photodissociatation

      j = j+1
      jlabel(j) = 'MGLY -> HO2 + CO + CCO_O2' 

* cross section 
      n = ns
      x1(1:n)=xs1(1:n) 
      ywork(1:n)=s1(1:n)
      
      CALL addpnt(x1,ywork,kdata,n,x1(1)*(1.-deltax),0.)
      CALL addpnt(x1,ywork,kdata,n,               0.,0.)
      CALL addpnt(x1,ywork,kdata,n,x1(n)*(1.+deltax),0.)
      CALL addpnt(x1,ywork,kdata,n,           1.e+38,0.)
      CALL inter2(nw,wl,yg,n,x1,ywork,ierr)
      IF (ierr .NE. 0) THEN
         WRITE(*,*) ierr, jlabel(j)
         STOP
      ENDIF

* quantum yield 
      n=ns
      x1(1:n)=xs1(1:n) 
      ywork(1:n)=y1(1:n)
       
      CALL addpnt(x1,ywork,kdata,n,x1(1)*(1.-deltax),0.)
      CALL addpnt(x1,ywork,kdata,n,               0.,0.)
      CALL addpnt(x1,ywork,kdata,n,x1(n)*(1.+deltax),0.)
      CALL addpnt(x1,ywork,kdata,n,           1.e+38,0.)
      CALL inter2(nw,wl,yg1,n,x1,ywork,ierr)
      IF (ierr .NE. 0) THEN
         WRITE(*,*) ierr, jlabel(j)
         STOP
      ENDIF

      DO iw = 1, nw - 1
         DO i = 1, nz
            sq(j,i,iw) = yg(iw)*yg1(iw)
         ENDDO
      ENDDO
      end
      

      SUBROUTINE rn_acet(nw,wl,wc,nz,tlev,airlev,j,sq,jlabel)

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Provide product (cross section) x (quantum yield) for CH3COCH3 photolysis=*
*=          CH3COCH3 + hv -> Products                                        =*
*=                                                                           =*
*=  Cross section:  Choice between                                           =*
*=                   (1) Calvert and Pitts                                   =*
*=                   (2) Martinez et al., 1991, alson in IUPAC 97            =*
*=                   (3) NOAA, 1998, unpublished as of 01/98                 =*
*=  Quantum yield:  Choice between                                           =*
*=                   (1) Gardiner et al, 1984                                =*
*=                   (2) IUPAC 97                                            =*
*=                   (3) McKeen et al., 1997                                 =*
*-----------------------------------------------------------------------------*
*=  PARAMETERS:                                                              =*
*=  NW     - INTEGER, number of specified intervals + 1 in working        (I)=*
*=           wavelength grid                                                 =*
*=  WL     - REAL, vector of lower limits of wavelength intervals in      (I)=*
*=           working wavelength grid                                         =*
*=  WC     - REAL, vector of center points of wavelength intervals in     (I)=*
*=           working wavelength grid                                         =*
*=  NZ     - INTEGER, number of altitude levels in working altitude grid  (I)=*
*=  TLEV   - REAL, temperature (K) at each specified altitude level       (I)=*
*=  AIRLEV - REAL, air density (molec/cc) at each altitude level          (I)=*
*=  J      - INTEGER, counter for number of weighting functions defined  (IO)=*
*=  SQ     - REAL, cross section x quantum yield (cm^2) for each          (O)=*
*=           photolysis reaction defined, at each defined wavelength and     =*
*=           at each defined altitude level                                  =*
*=  JLABEL - CHARACTER*40, string identifier for each photolysis reaction (O)=*
*=           defined                                                         =*
*-----------------------------------------------------------------------------*
*=  EDIT HISTORY:                                                            =*
*=  05/98  Original, adapted from former JSPEC1 subroutine                   =*
*-----------------------------------------------------------------------------*
*= This program is free software;  you can redistribute it and/or modify     =*
*= it under the terms of the GNU General Public License as published by the  =*
*= Free Software Foundation;  either version 2 of the license, or (at your   =*
*= option) any later version.                                                =*
*= The TUV package is distributed in the hope that it will be useful, but    =*
*= WITHOUT ANY WARRANTY;  without even the implied warranty of MERCHANTIBI-  =*
*= LITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public     =*
*= License for more details.                                                 =*
*= To obtain a copy of the GNU General Public License, write to:             =*
*= Free Software Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.   =*
*-----------------------------------------------------------------------------*

      IMPLICIT NONE
      INCLUDE 'tuv.params'

* input

      INTEGER nw
      REAL wl(kw), wc(kw)
      
      INTEGER nz

      REAL tlev(kz)
      REAL airlev(kz)

* weighting functions

      CHARACTER*40 jlabel(kj)
      REAL sq(kj,kz,kw)

* input/output:

      INTEGER j

* data arrays

      INTEGER kdata,ns
      PARAMETER(kdata=150,ns=96)

      INTEGER i, n
      INTEGER n1, n2, n3
      REAL x1(kdata), x2(kdata), x3(kdata), xwork(kdata), ywork(kdata)
      REAL y1(kdata), y2(kdata), y3(kdata)
      
       data x1(1:ns)/202.,206.,210.,214.,218.,222.,226.,
     1 230.,234.,238.,242.,246.,250.,254.,
     2 258.,262.,266.,270.,274.,278.,280.,
     3 281.,282.,283.,284.,285.,286.,287.,
     4 288.,289.,290.,291.,292.,293.,294.,
     5 295.,296.,297.,298.,299.,300.,301.,
     6 302.,303.,304.,305.,306.,307.,308.,
     7 309.,310.,311.,312.,313.,314.,315.,
     8 316.,317.,318.,319.,320.,321.,322.,
     9 323.,324.,325.,326.,327.,328.,329.,
     a 330.,331.,332.,333.,334.,335.,336.,
     b 337.,338.,339.,340.,341.,342.,343.,
     c 344.,345.,346.,347.,348.,349.,350.,
     d 351.,352.,353.,354.,355./
       data y1(1:ns)/0.533,0.125,0.104,0.120,0.163,0.242,0.361,
     1 0.533,0.774,1.086,1.479,1.944,2.470,3.040,
     2 3.610,4.150,4.580,4.910,5.060,5.070,5.050,
     3 5.010,4.940,4.860,4.760,4.680,4.580,4.500,
     4 4.410,4.290,4.190,4.080,3.940,3.810,3.670,
     5 3.520,3.350,3.200,3.070,2.910,2.770,2.660,
     6 2.530,2.370,2.240,2.110,1.952,1.801,1.663,
     7 1.537,1.408,1.276,1.173,1.081,0.967,0.858,
     8 0.777,0.699,0.608,0.530,0.467,0.407,0.344,
     9 0.287,0.243,0.205,0.168,0.135,0.108,0.086,
     a 0.067,0.051,0.040,0.031,0.026,0.017,0.014,
     b 0.011,0.009,0.006,0.005,0.005,0.003,0.004,
     c 0.002,0.002,0.001,0.002,0.001,0.001,0.001,
     d 0.000,0.001,0.000,0.001,0.000/

* local

      REAL yg(kw), yg1(kw), yg2(kw), yg3(kw)
      REAL qy
      REAL sig
      INTEGER ierr
      INTEGER iw

      REAL a, b, t
      INTEGER mabs, myld
**************** CH3COCH3 photodissociation

      j = j + 1
      jlabel(j) = 'CH3COCH3 -> CCO_O2 + C_O2'

* options
* mabs for cross sections
* myld for quantum yields

* Absorption:
* 2:  Martinez et al. 1991, also in IUPAC'97


* Quantum yield
* 3:  McKeen, S. A., T. Gierczak, J. B. Burkholder, P. O. Wennberg, T. F. Hanisco,
*       E. R. Keim, R.-S. Gao, S. C. Liu, A. R. Ravishankara, and D. W. Fahey, 
*       The photochemistry of acetone in the upper troposphere:  a source of 
*       odd-hydrogen radicals, Geophys. Res. Lett., 24, 3177-3180, 1997.

      mabs = 2
      myld = 3

      n = ns
      DO i = 1, n
       ywork(i) = y1(i) * 1.e-20
       xwork(i) = x1(i)
      ENDDO      
         
      CALL addpnt(xwork,ywork,kdata,n,xwork(1)*(1.-deltax),0.)
      CALL addpnt(xwork,ywork,kdata,n,               0.,0.)
      CALL addpnt(xwork,ywork,kdata,n,xwork(n)*(1.+deltax),0.)
      CALL addpnt(xwork,ywork,kdata,n,           1.e+38,0.)
      CALL inter2(nw,wl,yg,n,xwork,ywork,ierr)
      IF (ierr .NE. 0) THEN
        WRITE(*,*) ierr, jlabel(j)
        STOP
      ENDIF


      DO iw = 1, nw - 1

         DO i = 1, nz

            sig = yg(iw)

            IF (wc(iw) .LE. 292.) THEN
               qy = 1.
            ELSEIF (wc(iw) .GE. 292.  .AND. wc(iw) .LT. 308. ) THEN
               a = -15.696 + 0.05707*wc(iw)
               b = EXP(-88.81+0.15161*wc(iw))
               qy = 1./(a + b*airlev(i))
            ELSEIF (wc(iw) .GE. 308.  .AND. wc(iw) .LT. 337. ) THEN
               a = -130.2 + 0.42884*wc(iw)
               b = EXP(-55.947+0.044913*wc(iw))
               qy = 1./(a + b*airlev(i))
            ELSEIF (wc(iw) .GE. 337.) THEN
               qy = 0.
            ENDIF

               qy = max(0., qy)
               qy = min(1., qy)
            

            sq(j,i,iw) = sig*qy

         ENDDO
      ENDDO

      END

      SUBROUTINE rn_MEK(nw,wl,wc,nz,tlev,airlev,j,sq,jlabel)

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Provide product (cross section) x (quantum yield) for CH3CHO photolysis: =*
*=           Ketone + hv -> CCO_O2+CCHO+RO2_R                                =*
*=                    tyh modify for SAPRC99 07/01                           =*
*-----------------------------------------------------------------------------*
*=  PARAMETERS:                                                              =*
*=  NW     - INTEGER, number of specified intervals + 1 in working        (I)=*
*=           wavelength grid                                                 =*
*=  WL     - REAL, vector of lower limits of wavelength intervals in      (I)=*
*=           working wavelength grid                                         =*
*=  WC     - REAL, vector of center points of wavelength intervals in     (I)=*
*=           working wavelength grid                                         =*
*=  NZ     - INTEGER, number of altitude levels in working altitude grid  (I)=*
*=  TLEV   - REAL, temperature (K) at each specified altitude level       (I)=*
*=  AIRLEV - REAL, air density (molec/cc) at each altitude level          (I)=*
*=  J      - INTEGER, counter for number of weighting functions defined  (IO)=*
*=  SQ     - REAL, cross section x quantum yield (cm^2) for each          (O)=*
*=           photolysis reaction defined, at each defined wavelength and     =*
*=           at each defined altitude level                                  =*
*=  JLABEL - CHARACTER*40, string identifier for each photolysis reaction (O)=*
*=           defined                                                         =*
*-----------------------------------------------------------------------------*
*= This program is free software;  you can redistribute it and/or modify     =*
*= it under the terms of the GNU General Public License as published by the  =*
*= Free Software Foundation;  either version 2 of the license, or (at your   =*
*= option) any later version.                                                =*
*= The TUV package is distributed in the hope that it will be useful, but    =*
*= WITHOUT ANY WARRANTY;  without even the implied warranty of MERCHANTIBI-  =*
*= LITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public     =*
*= License for more details.                                                 =*
*= To obtain a copy of the GNU General Public License, write to:             =*
*= Free Software Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.   =*
*-----------------------------------------------------------------------------*

      IMPLICIT NONE
      INCLUDE 'tuv.params'

* input

      INTEGER nw
      REAL wl(kw), wc(kw)
      
      INTEGER nz

      REAL tlev(kz)
      REAL airlev(kz)

* weighting functions

      CHARACTER*40 jlabel(kj)
      REAL sq(kj,kz,kw)

* input/output:

      INTEGER j

* data arrays

      INTEGER kdata,ns
      PARAMETER(kdata=300,ns=290)

      INTEGER i, n
      REAL x1(kdata), s1(kdata), y1(kdata), xwork(kdata), ywork(kdata)
       data x1(1:ns)/198.5,199.0,199.5,200.0,200.5,201.0,201.5,
     1 202.0,202.5,203.0,203.5,204.0,204.5,205.0,
     2 205.5,206.0,206.5,207.0,207.5,208.0,208.5,
     3 209.0,209.5,210.0,210.5,211.0,211.5,212.0,
     4 212.5,213.0,213.5,214.0,214.5,215.0,215.5,
     5 216.0,216.5,217.0,217.5,218.0,218.5,219.0,
     6 219.5,220.0,220.5,221.0,221.5,222.0,222.5,
     7 223.0,223.5,224.0,224.5,225.0,225.5,226.0,
     8 226.5,227.0,227.5,228.0,228.5,229.0,229.5,
     9 230.0,230.5,231.0,231.5,232.0,232.5,233.0,
     a 233.5,234.0,234.5,235.0,235.5,236.0,236.5,
     b 237.0,237.5,238.0,238.5,239.0,239.5,240.0,
     c 240.5,241.0,241.5,242.0,242.5,243.0,243.5,
     d 244.0,244.5,245.0,245.5,246.0,246.5,247.0,
     e 247.5,248.0,248.5,249.0,249.5,250.0,250.5,
     f 251.0,251.5,252.0,252.5,253.0,253.5,254.0,
     g 254.5,255.0,255.5,256.0,256.5,257.0,257.5,
     h 258.0,258.5,259.0,259.5,260.0,260.5,261.0,
     i 261.5,262.0,262.5,263.0,263.5,264.0,264.5,
     j 265.0,265.5,266.0,266.5,267.0,267.5,268.0,
     k 268.5,269.0,269.5,270.0,270.5,271.0,271.5,
     l 272.0,272.5,273.0,273.5,274.0,274.5,275.0,
     m 275.5,276.0,276.5,277.0,277.5,278.0,278.5,
     n 279.0,279.5,280.0,280.5,281.0,281.5,282.0,
     o 282.5,283.0,283.5,284.0,284.5,285.0,285.5,
     p 286.0,286.5,287.0,287.5,288.0,288.5,289.0,
     q 289.5,290.0,290.5,291.0,291.5,292.0,292.5,
     r 293.0,293.5,294.0,294.5,295.0,295.5,296.0,
     s 296.5,297.0,297.5,298.0,298.5,299.0,299.5,
     t 300.0,300.5,301.0,301.5,302.0,302.5,303.0,
     u 303.5,304.0,304.5,305.0,305.5,306.0,306.5,
     v 307.0,307.5,308.0,308.5,309.0,309.5,310.0,
     x 310.5,311.0,311.5,312.0,312.5,313.0,313.5,
     y 314.0,314.5,315.0,315.5,316.0,316.5,317.0,
     z 317.5,318.0,318.5,319.0,319.5,320.0,320.5,
     A 321.0,321.5,322.0,322.5,323.0,323.5,324.0,
     B 324.5,325.0,325.5,326.0,326.5,327.0,327.5,
     C 328.0,328.5,329.0,329.5,330.0,330.5,331.0,
     D 331.5,332.0,332.5,333.0,333.5,334.0,334.5,
     E 335.0,335.5,336.0,336.5,337.0,337.5,338.0,
     F 338.5,339.0,339.5,340.0,340.5,341.0,341.5,
     G 342.0,342.5,343.0/
       data s1(1:ns)/.395E-18,.161E-18,.775E-19,.376E-19,.251E-19,
     1 .183E-19,.136E-19,.116E-19,.897E-20,.462E-20,.318E-20,.242E-20,
     2 .201E-20,.177E-20,.164E-20,.154E-20,.152E-20,.154E-20,.162E-20,
     3 .164E-20,.160E-20,.157E-20,.149E-20,.147E-20,.152E-20,.150E-20,
     4 .162E-20,.181E-20,.210E-20,.223E-20,.206E-20,.169E-20,.149E-20,
     5 .142E-20,.142E-20,.142E-20,.148E-20,.148E-20,.153E-20,.156E-20,
     6 .167E-20,.168E-20,.178E-20,.185E-20,.192E-20,.201E-20,.211E-20,
     7 .223E-20,.233E-20,.248E-20,.260E-20,.274E-20,.285E-20,.304E-20,
     8 .315E-20,.333E-20,.355E-20,.373E-20,.393E-20,.411E-20,.434E-20,
     9 .456E-20,.475E-20,.501E-20,.527E-20,.553E-20,.583E-20,.615E-20,
     a .645E-20,.673E-20,.702E-20,.742E-20,.783E-20,.811E-20,.845E-20,
     b .882E-20,.921E-20,.965E-20,.100E-19,.105E-19,.110E-19,.115E-19,
     c .120E-19,.123E-19,.128E-19,.132E-19,.138E-19,.144E-19,.150E-19,
     d .157E-19,.163E-19,.168E-19,.175E-19,.181E-19,.188E-19,.196E-19,
     e .203E-19,.211E-19,.219E-19,.225E-19,.233E-19,.240E-19,.248E-19,
     f .256E-19,.264E-19,.273E-19,.281E-19,.288E-19,.298E-19,.307E-19,
     g .316E-19,.325E-19,.334E-19,.343E-19,.351E-19,.359E-19,.367E-19,
     h .375E-19,.384E-19,.394E-19,.403E-19,.413E-19,.422E-19,.428E-19,
     i .433E-19,.441E-19,.449E-19,.457E-19,.465E-19,.472E-19,.478E-19,
     j .485E-19,.492E-19,.499E-19,.504E-19,.512E-19,.522E-19,.528E-19,
     k .534E-19,.541E-19,.546E-19,.551E-19,.555E-19,.559E-19,.563E-19,
     l .566E-19,.570E-19,.574E-19,.578E-19,.581E-19,.586E-19,.590E-19,
     m .593E-19,.596E-19,.597E-19,.598E-19,.598E-19,.599E-19,.599E-19,
     n .598E-19,.596E-19,.596E-19,.595E-19,.594E-19,.592E-19,.590E-19,
     o .588E-19,.586E-19,.583E-19,.579E-19,.575E-19,.571E-19,.567E-19,
     p .561E-19,.556E-19,.551E-19,.545E-19,.541E-19,.537E-19,.533E-19,
     q .527E-19,.521E-19,.515E-19,.508E-19,.499E-19,.489E-19,.482E-19,
     r .473E-19,.462E-19,.453E-19,.441E-19,.432E-19,.423E-19,.415E-19,
     s .411E-19,.401E-19,.394E-19,.388E-19,.377E-19,.369E-19,.363E-19,
     t .354E-19,.346E-19,.336E-19,.324E-19,.316E-19,.306E-19,.295E-19,
     u .282E-19,.270E-19,.259E-19,.249E-19,.242E-19,.234E-19,.228E-19,
     v .219E-19,.211E-19,.204E-19,.193E-19,.188E-19,.180E-19,.173E-19,
     x .166E-19,.158E-19,.148E-19,.142E-19,.134E-19,.126E-19,.117E-19,
     y .113E-19,.108E-19,.104E-19,.969E-20,.891E-20,.861E-20,.788E-20,
     z .725E-20,.692E-20,.643E-20,.607E-20,.564E-20,.519E-20,.466E-20,
     A .436E-20,.395E-20,.364E-20,.338E-20,.317E-20,.280E-20,.262E-20,
     B .229E-20,.213E-20,.193E-20,.170E-20,.158E-20,.148E-20,.124E-20,
     C .120E-20,.104E-20,.951E-21,.844E-21,.726E-21,.670E-21,.608E-21,
     D .515E-21,.456E-21,.413E-21,.356E-21,.330E-21,.297E-21,.267E-21,
     E .246E-21,.221E-21,.193E-21,.156E-21,.147E-21,.137E-21,.127E-21,
     F .119E-21,.109E-21,.101E-21,.909E-22,.822E-22,.766E-22,.743E-22,
     G .683E-22,.672E-22,.604E-22,.478E-22,.000E+00/

* local

      REAL yg(kw), yg1, yg2(kw), yg3(kw), yg4(kw)
      REAL qy1, qy2, qy3
      REAL sig
      REAL dum
      INTEGER ierr
      INTEGER  iz, iw


************************* MEK photolysis

      j = j+1
      jlabel(j) = 'MEK -> CCO_O2 + CCHO + RO2_R'

* cross section 
      n = ns
      xwork(1:n)=x1(1:n)
      ywork(1:n)=s1(1:n)
      
      CALL addpnt(xwork,ywork,kdata,n,xwork(1)*(1.-deltax),0.)
      CALL addpnt(xwork,ywork,kdata,n,               0.,0.)
      CALL addpnt(xwork,ywork,kdata,n,xwork(n)*(1.+deltax),0.)
      CALL addpnt(xwork,ywork,kdata,n,           1.e+38,0.)
      CALL inter2(nw,wl,yg,n,xwork,ywork,ierr)
      IF (ierr .NE. 0) THEN
         WRITE(*,*) ierr, jlabel(j)
         STOP
      ENDIF

* all quantum yield equal 0.15

      yg1=0.15 
      DO iw = 1, nw - 1             
         DO i = 1, nz	  
            sq(j,i,iw) = yg(iw)*yg1
         ENDDO
      ENDDO


      j = j+1
      jlabel(j) = 'PROD2->RO2_R+RO2_N+R2O2+CCO_O2+RCO_O2+HCHO+CCHO+RCHO'

* all quantum yield equal 0.02

      yg1=0.02 
      DO iw = 1, nw - 1             
         DO i = 1, nz	  
            sq(j,i,iw) = yg(iw)*yg1
         ENDDO
      ENDDO

      END

      
      SUBROUTINE rn_COOH(nw,wl,wc,nz,tlev,airlev,j,sq,jlabel)

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Provide product (cross section) x (quantum yield) for CH3OOH photolysis: =*
*=         CH3OOH + hv -> CH3O + OH                                          =*
*=                                                                           =*
*=  Cross section: Choice between                                            =*
*=                  (1) JPL 97 recommendation (based on Vaghjiana and        =*
*=                      Ravishankara, 1989), 10 nm resolution                =*
*=                  (2) IUPAC 97 (from Vaghjiana and Ravishankara, 1989),    =*
*=                      5 nm resolution                                      =*
*=                  (3) Cox and Tyndall, 1978; only for wavelengths < 280 nm =*
*=                  (4) Molina and Arguello, 1979;  might be 40% too high    =*
*=  Quantum yield: Assumed to be unity                                       =*
*-----------------------------------------------------------------------------*
*=  PARAMETERS:                                                              =*
*=  NW     - INTEGER, number of specified intervals + 1 in working        (I)=*
*=           wavelength grid                                                 =*
*=  WL     - REAL, vector of lower limits of wavelength intervals in      (I)=*
*=           working wavelength grid                                         =*
*=  WC     - REAL, vector of center points of wavelength intervals in     (I)=*
*=           working wavelength grid                                         =*
*=  NZ     - INTEGER, number of altitude levels in working altitude grid  (I)=*
*=  TLEV   - REAL, temperature (K) at each specified altitude level       (I)=*
*=  AIRLEV - REAL, air density (molec/cc) at each altitude level          (I)=*
*=  J      - INTEGER, counter for number of weighting functions defined  (IO)=*
*=  SQ     - REAL, cross section x quantum yield (cm^2) for each          (O)=*
*=           photolysis reaction defined, at each defined wavelength and     =*
*=           at each defined altitude level                                  =*
*=  JLABEL - CHARACTER*40, string identifier for each photolysis reaction (O)=*
*=           defined                                                         =*
*-----------------------------------------------------------------------------*
*=  EDIT HISTORY:                                                            =*
*=  05/98  Original, adapted from former JSPEC1 subroutine                   =*
*-----------------------------------------------------------------------------*
*= This program is free software;  you can redistribute it and/or modify     =*
*= it under the terms of the GNU General Public License as published by the  =*
*= Free Software Foundation;  either version 2 of the license, or (at your   =*
*= option) any later version.                                                =*
*= The TUV package is distributed in the hope that it will be useful, but    =*
*= WITHOUT ANY WARRANTY;  without even the implied warranty of MERCHANTIBI-  =*
*= LITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public     =*
*= License for more details.                                                 =*
*= To obtain a copy of the GNU General Public License, write to:             =*
*= Free Software Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.   =*
*-----------------------------------------------------------------------------*

      IMPLICIT NONE
      INCLUDE 'tuv.params'

* input

      INTEGER nw
      REAL wl(kw), wc(kw)
      
      INTEGER nz

      REAL tlev(kz)
      REAL airlev(kz)

* weighting functions

      CHARACTER*40 jlabel(kj)
      REAL sq(kj,kz,kw)

* input/output:

      INTEGER j

* data arrays

      INTEGER kdata,ns
      PARAMETER(kdata=100,ns=32)

      INTEGER i, n
      REAL x1(kdata),y1(kdata),xwork(kdata),ywork(kdata)

       data x1(1:ns)/210.,215.,220.,225.,230.,235.,240.,
     1  245.,250.,255.,260.,265.,270.,275.,
     2  280.,285.,290.,295.,300.,305.,310.,
     3  315.,320.,325.,330.,335.,340.,345.,
     4  350.,355.,360.,365./
       data y1(1:ns)/31.200,20.900,15.400,12.200, 9.620, 7.610, 6.050,
     1  4.880, 3.980, 3.230, 2.560, 2.110, 1.700, 1.390,
     2  1.090, 0.863, 0.691, 0.551, 0.413, 0.313, 0.239,
     3  0.182, 0.137, 0.105, 0.079, 0.061, 0.047, 0.035,
     4  0.027, 0.021, 0.016, 0.012/

* local

      REAL yg(kw)
      REAL qy
      INTEGER ierr
      INTEGER idum
      INTEGER iw

      INTEGER mabs


**************** CH3OOH photodissociation, ROOH is same
         j = j + 1
         jlabel(j) = 'COOH -> HCHO + HO2 + OH'
         j = j + 1
         jlabel(j) = 'ROOH -> RCHO + HO2 + OH'

* mabs: Absorption cross section options:
* 2:  IUPAC97 (from  Vaghjiani and Ravishankara (1989) at 5 nm resolution).


         n = ns
         DO i = 1, n
            ywork(i) = y1(i) * 1.E-20
	    xwork(i) = x1(i)
         ENDDO

         CALL addpnt(xwork,ywork,kdata,n,xwork(1)*(1.-deltax),0.)
         CALL addpnt(xwork,ywork,kdata,n,               0.,0.)
         CALL addpnt(xwork,ywork,kdata,n,xwork(n)*(1.+deltax),0.)
         CALL addpnt(xwork,ywork,kdata,n,           1.e+38,0.)
         CALL inter2(nw,wl,yg,n,xwork,ywork,ierr)
         IF (ierr .NE. 0) THEN
            WRITE(*,*) ierr, jlabel(j)
            STOP
         ENDIF

* quantum yield = 1

      qy = 1.
      DO iw = 1, nw - 1
         DO i = 1, nz
	    sq(j-1,i,iw) = yg(iw)*qy
            sq(j,i,iw) = yg(iw)*qy
         ENDDO
      ENDDO

      END

      SUBROUTINE rn_BACL (nw,wl,wc,nz,tlev,airlev,j,sq,jlabel)

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Provide product of (cross section) x (quantum yield) for CH2O photolysis =*
*=        (a) BACL + hv = 2 CCO_O2                                           =*
*=        (b) BALD + hv = products                                           =* 
*=                tyh created for SAPRC99   07/01                            =*
*-----------------------------------------------------------------------------*
*=  PARAMETERS:                                                              =*
*=  NW     - INTEGER, number of specified intervals + 1 in working        (I)=*
*=           wavelength grid                                                 =*
*=  WL     - REAL, vector of lower limits of wavelength intervals in      (I)=*
*=           working wavelength grid                                         =*
*=  WC     - REAL, vector of center points of wavelength intervals in     (I)=*
*=           working wavelength grid                                         =*
*=  NZ     - INTEGER, number of altitude levels in working altitude grid  (I)=*
*=  TLEV   - REAL, temperature (K) at each specified altitude level       (I)=*
*=  AIRLEV - REAL, air density (molec/cc) at each altitude level          (I)=*
*=  J      - INTEGER, counter for number of weighting functions defined  (IO)=*
*=  SQ     - REAL, cross section x quantum yield (cm^2) for each          (O)=*
*=           photolysis reaction defined, at each defined wavelength and     =*
*=           at each defined altitude level                                  =*
*=  JLABEL - CHARACTER*40, string identifier for each photolysis reaction (O)=*
*=           defined                                                         =*
*-----------------------------------------------------------------------------*
*=  EDIT HISTORY:                                                            =*
*=  05/98  Original, adapted from former JSPEC1 subroutine                   =*
*-----------------------------------------------------------------------------*
*= This program is free software;  you can redistribute it and/or modify     =*
*= it under the terms of the GNU General Public License as published by the  =*
*= Free Software Foundation;  either version 2 of the license, or (at your   =*
*= option) any later version.                                                =*
*= The TUV package is distributed in the hope that it will be useful, but    =*
*= WITHOUT ANY WARRANTY;  without even the implied warranty of MERCHANTIBI-  =*
*= LITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public     =*
*= License for more details.                                                 =*
*= To obtain a copy of the GNU General Public License, write to:             =*
*= Free Software Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.   =*
*-----------------------------------------------------------------------------*

      IMPLICIT NONE
      INCLUDE 'tuv.params'

      INTEGER kdata
      PARAMETER(kdata=160)

* input

      INTEGER nw
      REAL wl(kw), wc(kw)
      
      INTEGER nz

      REAL tlev(kz)
      REAL airlev(kz)

* weighting functions

      CHARACTER*40 jlabel(kj)
      REAL sq(kj,kz,kw)

* input/output:

      INTEGER j, iz, iw

* data arrays

      INTEGER n,nm,nr
      parameter(nr=98,nm=26)
      real x1r(kdata), x1m(kdata), s1r(kdata),s1m(kdata),
     1  y1r(kdata) ,x1(kdata),y1(kdata)
      INTEGER n1, n2, n3, n4, n5

       data x1r(1:nr)/230.0,232.5,235.0,237.5,240.0,242.5,245.0,
     1 247.5,250.0,252.5,255.0,257.5,260.0,262.5,
     2 265.0,267.5,270.0,272.5,275.0,277.5,280.0,
     3 282.5,285.0,287.5,290.0,292.5,295.0,297.5,
     4 300.0,302.5,305.0,307.5,310.0,312.5,315.0,
     5 317.5,320.0,322.5,325.0,327.5,330.0,332.5,
     6 335.0,337.5,340.0,342.5,345.0,347.5,350.0,
     7 352.5,355.0,357.5,360.0,362.5,365.0,367.5,
     8 370.0,372.5,375.0,377.5,380.0,382.5,385.0,
     9 387.5,390.0,392.5,395.0,397.5,400.0,402.5,
     a 405.0,407.5,410.0,412.5,415.0,417.5,420.0,
     b 422.5,425.0,427.5,430.0,432.5,435.0,437.5,
     c 440.0,442.5,445.0,447.5,450.0,452.5,455.0,
     d 457.5,460.0,462.5,465.0,467.5,470.0,472.5/
       data s1r(1:nr)/.130E-19,.146E-19,.168E-19,.184E-19,.216E-19,
     1 .249E-19,.265E-19,.271E-19,.303E-19,.346E-19,.346E-19,.357E-19,
     2 .395E-19,.417E-19,.417E-19,.422E-19,.460E-19,.454E-19,.433E-19,
     3 .422E-19,.444E-19,.433E-19,.390E-19,.357E-19,.325E-19,.292E-19,
     4 .260E-19,.216E-19,.179E-19,.173E-19,.146E-19,.108E-19,.920E-20,
     5 .703E-20,.649E-20,.541E-20,.541E-20,.541E-20,.433E-20,.325E-20,
     6 .379E-20,.379E-20,.433E-20,.487E-20,.541E-20,.595E-20,.649E-20,
     7 .703E-20,.812E-20,.757E-20,.920E-20,.974E-20,.108E-19,.119E-19,
     8 .141E-19,.151E-19,.179E-19,.200E-19,.211E-19,.233E-19,.260E-19,
     9 .281E-19,.314E-19,.346E-19,.390E-19,.411E-19,.433E-19,.438E-19,
     a .465E-19,.481E-19,.519E-19,.584E-19,.606E-19,.649E-19,.692E-19,
     b .687E-19,.682E-19,.671E-19,.649E-19,.595E-19,.573E-19,.628E-19,
     c .601E-19,.584E-19,.595E-19,.649E-19,.595E-19,.498E-19,.379E-19,
     d .281E-19,.173E-19,.108E-19,.541E-20,.379E-20,.216E-20,.108E-20,
     e .108E-20,.000E+00/
       data y1r(1:nr)/1.000,1.000,1.000,1.000,1.000,1.000,1.000,
     1 1.000,1.000,1.000,1.000,1.000,1.000,1.000,
     2 1.000,1.000,1.000,1.000,1.000,1.000,1.000,
     3 1.000,1.000,1.000,1.000,1.000,1.000,1.000,
     4 1.000,1.000,1.000,1.000,1.000,1.000,1.000,
     5 1.000,1.000,1.000,1.000,1.000,1.000,1.000,
     6 1.000,1.000,1.000,1.000,1.000,1.000,0.995,
     7 0.960,0.925,0.890,0.855,0.820,0.785,0.750,
     8 0.715,0.680,0.645,0.610,0.575,0.540,0.505,
     9 0.470,0.435,0.399,0.364,0.329,0.294,0.259,
     a 0.224,0.189,0.154,0.119,0.084,0.049,0.014,
     b 0.000,0.000,0.000,0.000,0.000,0.000,0.000,
     c 0.000,0.000,0.000,0.000,0.000,0.000,0.000,
     d 0.000,0.000,0.000,0.000,0.000,0.000,0.000/
       data x1m(1:nm)/299.0,304.0,306.0,309.0,313.0,314.0,318.0,
     1 325.0,332.0,338.0,342.0,346.0,349.0,354.0,
     2 355.0,364.0,368.0,369.0,370.0,372.0,374.0,
     3 376.0,377.0,380.0,382.0,386.0/
       data s1m(1:nm)/.178E-18,.740E-19,.691E-19,.641E-19,.691E-19,
     1 .691E-19,.641E-19,.839E-19,.765E-19,.888E-19,.888E-19,.789E-19,
     2 .789E-19,.913E-19,.814E-19,.567E-19,.666E-19,.839E-19,.839E-19,
     3 .345E-19,.321E-19,.247E-19,.247E-19,.358E-19,.990E-20,.000E+00/

* local

      REAL yg(kw), yg1(kw), yg2(kw), yg3(kw), yg4(kw), yg5(kw)
      REAL a, b, c
      REAL a0, a1, a2, a3, a4, a5, a6, a7
      REAL b0, b1, b2, b3, b4
      REAL phi1, phi2, phi20, ak300, akt
      REAL qy, qy1, qy2, qy3

      REAL sigma, sig, slope
      REAL xs
      REAL t
      REAL dum
      INTEGER idum

      INTEGER i
      INTEGER irow, icol, irev
      INTEGER ierr

      INTEGER mopt1, mopt2


****************************************************************

      j = j+1
      jlabel(j) = 'BACL -> 2 CCO_O2' 

* cross section 
      n = nr
      x1(1:n)=x1r(1:n) 
      y1(1:n)=s1r(1:n)
      
      CALL addpnt(x1,y1,kdata,n,x1(1)*(1.-deltax),0.)
      CALL addpnt(x1,y1,kdata,n,               0.,0.)
      CALL addpnt(x1,y1,kdata,n,x1(n)*(1.+deltax),0.)
      CALL addpnt(x1,y1,kdata,n,           1.e+38,0.)
      CALL inter2(nw,wl,yg,n,x1,y1,ierr)
      IF (ierr .NE. 0) THEN
         WRITE(*,*) ierr, jlabel(j)
         STOP
      ENDIF

* quantum yield 
      n=nr
      x1(1:n)=x1r(1:n) 
      y1(1:n)=y1r(1:n)
      
      CALL addpnt(x1,y1,kdata,n,x1(1)*(1.-deltax),0.)
      CALL addpnt(x1,y1,kdata,n,               0.,0.)
      CALL addpnt(x1,y1,kdata,n,x1(n)*(1.+deltax),0.)
      CALL addpnt(x1,y1,kdata,n,           1.e+38,0.)
      CALL inter2(nw,wl,yg1,n,x1,y1,ierr)
      IF (ierr .NE. 0) THEN
         WRITE(*,*) ierr, jlabel(j)
         STOP
      ENDIF

      DO iw = 1, nw - 1
         DO i = 1, nz
            sq(j,i,iw) = yg(iw)*yg1(iw)
         ENDDO
      ENDDO

      j = j+1
      jlabel(j) = 'BALD -> Products'

* cross section 
      n = nm
      x1(1:n)=x1m(1:n) 
      y1(1:n)=s1m(1:n)
      
      CALL addpnt(x1,y1,kdata,n,x1(1)*(1.-deltax),0.)
      CALL addpnt(x1,y1,kdata,n,               0.,0.)
      CALL addpnt(x1,y1,kdata,n,x1(n)*(1.+deltax),0.)
      CALL addpnt(x1,y1,kdata,n,           1.e+38,0.)
      CALL inter2(nw,wl,yg,n,x1,y1,ierr)
      IF (ierr .NE. 0) THEN
         WRITE(*,*) ierr, jlabel(j)
         STOP
      ENDIF

* quantum yield is 0.05 for all wavelength
      qy=0.05
       
      DO iw = 1, nw - 1
         DO i = 1, nz
            sq(j,i,iw) = yg(iw)*qy
         ENDDO
      ENDDO

      END

      SUBROUTINE rn_ACROLEIN(nw,wl,wc,nz,tlev,airlev,j,sq,jlabel)

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Provide product (cross section) x (quantum yield) for CH3OOH photolysis: =*
*=         METHACRO+hv=HO2+RO2_R+OH+CCO_O2+CO+HCHO+MA_RCO3                   =*
*=         MVK+hv=C_O2+CO+PROD2+MA_RCO3                                      =*
*=         ISOPROD+hv=HO2+CCO_O2+RCO_O2+CO+HCHO+CCHO+MEK                     =*
*=         DCB3+hv=RO2_R+CCO_O2+HO2+CO+R2O2+GLY+MGLY                         =*
*=                                                                           =*
*=  tyh created for SAPRC99                                                  =*
*-----------------------------------------------------------------------------*
*=  PARAMETERS:                                                              =*
*=  NW     - INTEGER, number of specified intervals + 1 in working        (I)=*
*=           wavelength grid                                                 =*
*=  WL     - REAL, vector of lower limits of wavelength intervals in      (I)=*
*=           working wavelength grid                                         =*
*=  WC     - REAL, vector of center points of wavelength intervals in     (I)=*
*=           working wavelength grid                                         =*
*=  NZ     - INTEGER, number of altitude levels in working altitude grid  (I)=*
*=  TLEV   - REAL, temperature (K) at each specified altitude level       (I)=*
*=  AIRLEV - REAL, air density (molec/cc) at each altitude level          (I)=*
*=  J      - INTEGER, counter for number of weighting functions defined  (IO)=*
*=  SQ     - REAL, cross section x quantum yield (cm^2) for each          (O)=*
*=           photolysis reaction defined, at each defined wavelength and     =*
*=           at each defined altitude level                                  =*
*=  JLABEL - CHARACTER*40, string identifier for each photolysis reaction (O)=*
*=           defined                                                         =*
*-----------------------------------------------------------------------------*
*=  EDIT HISTORY:                                                            =*
*=  05/98  Original, adapted from former JSPEC1 subroutine                   =*
*-----------------------------------------------------------------------------*
*= This program is free software;  you can redistribute it and/or modify     =*
*= it under the terms of the GNU General Public License as published by the  =*
*= Free Software Foundation;  either version 2 of the license, or (at your   =*
*= option) any later version.                                                =*
*= The TUV package is distributed in the hope that it will be useful, but    =*
*= WITHOUT ANY WARRANTY;  without even the implied warranty of MERCHANTIBI-  =*
*= LITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public     =*
*= License for more details.                                                 =*
*= To obtain a copy of the GNU General Public License, write to:             =*
*= Free Software Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.   =*
*-----------------------------------------------------------------------------*

      IMPLICIT NONE
      INCLUDE 'tuv.params'

* input

      INTEGER nw
      REAL wl(kw), wc(kw)
      
      INTEGER nz

      REAL tlev(kz)
      REAL airlev(kz)

* weighting functions

      CHARACTER*40 jlabel(kj)
      REAL sq(kj,kz,kw)

* input/output:

      INTEGER j

* data arrays

      INTEGER kdata,ns
      PARAMETER(kdata=150,ns=131)

      INTEGER i, n
      REAL x1(kdata), y1(kdata), xwork(kdata), ywork(kdata)

       data x1(1:ns)/250.0,252.0,253.0,254.0,255.0,256.0,257.0,
     1 258.0,259.0,260.0,261.0,262.0,263.0,264.0,
     2 265.0,266.0,267.0,268.0,269.0,270.0,271.0,
     3 272.0,273.0,274.0,275.0,276.0,277.0,278.0,
     4 279.0,280.0,281.0,282.0,283.0,284.0,285.0,
     5 286.0,287.0,288.0,289.0,290.0,291.0,292.0,
     6 293.0,294.0,295.0,296.0,297.0,298.0,299.0,
     7 300.0,301.0,302.0,303.0,304.0,305.0,306.0,
     8 307.0,308.0,309.0,310.0,311.0,312.0,313.0,
     9 314.0,315.0,316.0,317.0,318.0,319.0,320.0,
     a 321.0,322.0,323.0,324.0,325.0,326.0,327.0,
     b 328.0,329.0,330.0,331.0,332.0,333.0,334.0,
     c 335.0,336.0,337.0,338.0,339.0,340.0,341.0,
     d 342.0,343.0,344.0,345.0,346.0,347.0,348.0,
     e 349.0,350.0,351.0,352.0,353.0,354.0,355.0,
     f 356.0,357.0,358.0,359.0,360.0,361.0,362.0,
     g 363.0,364.0,365.0,366.0,367.0,368.0,369.0,
     h 370.0,371.0,372.0,373.0,374.0,375.0,376.0,
     i 377.0,378.0,379.0,380.0,381.0/
       data y1(1:ns)/.180E-20,.205E-20,.220E-20,.232E-20,.245E-20,
     1 .256E-20,.265E-20,.274E-20,.283E-20,.298E-20,.324E-20,.347E-20,
     2 .358E-20,.393E-20,.467E-20,.510E-20,.538E-20,.573E-20,.613E-20,
     3 .664E-20,.720E-20,.777E-20,.837E-20,.894E-20,.955E-20,.104E-19,
     4 .112E-19,.119E-19,.127E-19,.127E-19,.126E-19,.126E-19,.128E-19,
     5 .133E-19,.138E-19,.144E-19,.150E-19,.157E-19,.163E-19,.171E-19,
     6 .178E-19,.186E-19,.195E-19,.205E-19,.215E-19,.226E-19,.237E-19,
     7 .248E-19,.260E-19,.273E-19,.285E-19,.299E-19,.313E-19,.327E-19,
     8 .339E-19,.351E-19,.363E-19,.377E-19,.391E-19,.407E-19,.425E-19,
     9 .439E-19,.444E-19,.450E-19,.459E-19,.475E-19,.490E-19,.505E-19,
     a .519E-19,.531E-19,.543E-19,.552E-19,.560E-19,.567E-19,.567E-19,
     b .562E-19,.563E-19,.571E-19,.576E-19,.580E-19,.595E-19,.623E-19,
     c .639E-19,.638E-19,.624E-19,.601E-19,.579E-19,.563E-19,.556E-19,
     d .552E-19,.554E-19,.553E-19,.547E-19,.541E-19,.540E-19,.548E-19,
     e .590E-19,.608E-19,.600E-19,.553E-19,.503E-19,.450E-19,.403E-19,
     f .375E-19,.355E-19,.345E-19,.346E-19,.349E-19,.341E-19,.323E-19,
     g .295E-19,.281E-19,.291E-19,.325E-19,.354E-19,.330E-19,.278E-19,
     h .215E-19,.159E-19,.119E-19,.899E-20,.722E-20,.586E-20,.469E-20,
     i .372E-20,.357E-20,.355E-20,.283E-20,.169E-20,.829E-23,.000E+00/

* local

      REAL yg(kw)
      REAL qy
      INTEGER ierr
      INTEGER idum
      INTEGER iw

      INTEGER mabs


* Absorption cross section 

         n = ns
	 xwork(1:n)=x1(1:n)
	 ywork(1:n)=y1(1:n)

         CALL addpnt(xwork,ywork,kdata,n,xwork(1)*(1.-deltax),0.)
         CALL addpnt(xwork,ywork,kdata,n,               0.,0.)
         CALL addpnt(xwork,ywork,kdata,n,xwork(n)*(1.+deltax),0.)
         CALL addpnt(xwork,ywork,kdata,n,           1.e+38,0.)
         CALL inter2(nw,wl,yg,n,xwork,ywork,ierr)
         IF (ierr .NE. 0) THEN
            WRITE(*,*) ierr, jlabel(j)
            STOP
         ENDIF

      j=j+1
      jlabel(j)='METHACRO->HO2+RO2_R+OH+CCO_O2+CO+HCHO+MA_RCO3'

* quantum yield 0.0041

      qy = 0.0041
      DO iw = 1, nw - 1
         DO i = 1, nz
            sq(j,i,iw) = yg(iw)*qy
         ENDDO
      ENDDO
      

      j=j+1
      jlabel(j)='MVK->C_O2+CO+PROD2+MA_RCO3'

* quantum yield 0.0021

      qy = 0.0021
      DO iw = 1, nw - 1
         DO i = 1, nz
            sq(j,i,iw) = yg(iw)*qy
         ENDDO
      ENDDO       


      j=j+1
      jlabel(j)='ISOPROD->HO2+CCO_O2+RCO_O2+CO+HCHO+CCHO+MEK'

* quantum yield 0.0041

      qy = 0.0041
      DO iw = 1, nw - 1
         DO i = 1, nz
            sq(j,i,iw) = yg(iw)*qy
         ENDDO
      ENDDO       

      j=j+1
      jlabel(j)='DCB3->RO2_R+CCO_O2+HO2+CO+R2O2+GLY+MGLY'

* quantum yield 7.3

      qy = 7.3
      DO iw = 1, nw - 1
         DO i = 1, nz
            sq(j,i,iw) = yg(iw)*qy
         ENDDO
      ENDDO       

      END

      SUBROUTINE rn_RNO3(nw,wl,wc,nz,tlev,airlev,j,sq,jlabel)

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Provide product (cross section) x (quantum yield) for CH3OOH photolysis: =*
*=         RNO3 + hv -> Products                                             =*
*=                                                                           =*
*=    tyh created for SAPRC99                                                =*
*-----------------------------------------------------------------------------*
*=  PARAMETERS:                                                              =*
*=  NW     - INTEGER, number of specified intervals + 1 in working        (I)=*
*=           wavelength grid                                                 =*
*=  WL     - REAL, vector of lower limits of wavelength intervals in      (I)=*
*=           working wavelength grid                                         =*
*=  WC     - REAL, vector of center points of wavelength intervals in     (I)=*
*=           working wavelength grid                                         =*
*=  NZ     - INTEGER, number of altitude levels in working altitude grid  (I)=*
*=  TLEV   - REAL, temperature (K) at each specified altitude level       (I)=*
*=  AIRLEV - REAL, air density (molec/cc) at each altitude level          (I)=*
*=  J      - INTEGER, counter for number of weighting functions defined  (IO)=*
*=  SQ     - REAL, cross section x quantum yield (cm^2) for each          (O)=*
*=           photolysis reaction defined, at each defined wavelength and     =*
*=           at each defined altitude level                                  =*
*=  JLABEL - CHARACTER*40, string identifier for each photolysis reaction (O)=*
*=           defined                                                         =*
*-----------------------------------------------------------------------------*
*=  EDIT HISTORY:                                                            =*
*=  05/98  Original, adapted from former JSPEC1 subroutine                   =*
*-----------------------------------------------------------------------------*
*= This program is free software;  you can redistribute it and/or modify     =*
*= it under the terms of the GNU General Public License as published by the  =*
*= Free Software Foundation;  either version 2 of the license, or (at your   =*
*= option) any later version.                                                =*
*= The TUV package is distributed in the hope that it will be useful, but    =*
*= WITHOUT ANY WARRANTY;  without even the implied warranty of MERCHANTIBI-  =*
*= LITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public     =*
*= License for more details.                                                 =*
*= To obtain a copy of the GNU General Public License, write to:             =*
*= Free Software Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.   =*
*-----------------------------------------------------------------------------*

      IMPLICIT NONE
      INCLUDE 'tuv.params'

* input

      INTEGER nw
      REAL wl(kw), wc(kw)
      
      INTEGER nz

      REAL tlev(kz)
      REAL airlev(kz)

* weighting functions

      CHARACTER*40 jlabel(kj)
      REAL sq(kj,kz,kw)

* input/output:

      INTEGER j

* data arrays

      INTEGER kdata,ns
      PARAMETER(kdata=100,ns=32)

      INTEGER i, n
      REAL x1(kdata),y1(kdata),xwork(kdata),ywork(kdata)
      
       data x1(1:ns)/185.0,188.0,190.0,195.0,200.0,205.0,210.0,
     1  215.0,220.0,225.0,230.0,235.0,240.0,245.0,
     2  250.0,255.0,260.0,265.0,270.0,275.0,280.0,
     3  285.0,290.0,295.0,300.0,305.0,310.0,315.0,
     4  320.0,325.0,330.0,335.0/
       data y1(1:ns)/.179E-16,.181E-16,.179E-16,.161E-16,.126E-16,
     1 .867E-17,.498E-17,.247E-17,.117E-17,.580E-18,.310E-18,.180E-18,
     2 .110E-18,.700E-19,.570E-19,.520E-19,.490E-19,.460E-19,.410E-19,
     3 .360E-19,.290E-19,.230E-19,.170E-19,.120E-19,.810E-20,.520E-20,
     4 .320E-20,.190E-20,.110E-20,.610E-21,.370E-21,.000E+00/

* local

      REAL yg(kw)
      REAL qy
      INTEGER ierr
      INTEGER idum
      INTEGER iw

      INTEGER mabs


**************** CH3OOH photodissociation, ROOH is same
         j = j + 1
         jlabel(j) = 'RNO3 -> NO2+HO2+RO2_R....'

* mabs: Absorption cross section options:
* 2:  IUPAC97 (from  Vaghjiani and Ravishankara (1989) at 5 nm resolution).


         n = ns
	 xwork(1:n)=x1(1:n)
	 ywork(1:n)=y1(1:n)

         CALL addpnt(xwork,ywork,kdata,n,xwork(1)*(1.-deltax),0.)
         CALL addpnt(xwork,ywork,kdata,n,               0.,0.)
         CALL addpnt(xwork,ywork,kdata,n,xwork(n)*(1.+deltax),0.)
         CALL addpnt(xwork,ywork,kdata,n,           1.e+38,0.)
         CALL inter2(nw,wl,yg,n,xwork,ywork,ierr)
         IF (ierr .NE. 0) THEN
            WRITE(*,*) ierr, jlabel(j)
            STOP
         ENDIF

* quantum yield = 1

      qy = 1.
      DO iw = 1, nw - 1
         DO i = 1, nz
            sq(j,i,iw) = yg(iw)*qy
         ENDDO
      ENDDO

      END


      SUBROUTINE rn_DCB2(nw,wl,wc,nz,tlev,airlev,j,sq,jlabel)

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Provide product (cross section) x (quantum yield) for CH3OOH photolysis: =*
*=         DCB2 + hv -> Products                                             =*
*=                                                                           =* 
*=    tyh created for SAPRC99                                                =*
*=  Quantum yield: Assumed to be unity                                       =*
*-----------------------------------------------------------------------------*
*=  PARAMETERS:                                                              =*
*=  NW     - INTEGER, number of specified intervals + 1 in working        (I)=*
*=           wavelength grid                                                 =*
*=  WL     - REAL, vector of lower limits of wavelength intervals in      (I)=*
*=           working wavelength grid                                         =*
*=  WC     - REAL, vector of center points of wavelength intervals in     (I)=*
*=           working wavelength grid                                         =*
*=  NZ     - INTEGER, number of altitude levels in working altitude grid  (I)=*
*=  TLEV   - REAL, temperature (K) at each specified altitude level       (I)=*
*=  AIRLEV - REAL, air density (molec/cc) at each altitude level          (I)=*
*=  J      - INTEGER, counter for number of weighting functions defined  (IO)=*
*=  SQ     - REAL, cross section x quantum yield (cm^2) for each          (O)=*
*=           photolysis reaction defined, at each defined wavelength and     =*
*=           at each defined altitude level                                  =*
*=  JLABEL - CHARACTER*40, string identifier for each photolysis reaction (O)=*
*=           defined                                                         =*
*-----------------------------------------------------------------------------*
*=  EDIT HISTORY:                                                            =*
*=  05/98  Original, adapted from former JSPEC1 subroutine                   =*
*-----------------------------------------------------------------------------*
*= This program is free software;  you can redistribute it and/or modify     =*
*= it under the terms of the GNU General Public License as published by the  =*
*= Free Software Foundation;  either version 2 of the license, or (at your   =*
*= option) any later version.                                                =*
*= The TUV package is distributed in the hope that it will be useful, but    =*
*= WITHOUT ANY WARRANTY;  without even the implied warranty of MERCHANTIBI-  =*
*= LITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public     =*
*= License for more details.                                                 =*
*= To obtain a copy of the GNU General Public License, write to:             =*
*= Free Software Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.   =*
*-----------------------------------------------------------------------------*

      IMPLICIT NONE
      INCLUDE 'tuv.params'

* input

      INTEGER nw
      REAL wl(kw), wc(kw)
      
      INTEGER nz

      REAL tlev(kz)
      REAL airlev(kz)

* weighting functions

      CHARACTER*40 jlabel(kj)
      REAL sq(kj,kz,kw)

* input/output:

      INTEGER j

* data arrays

      INTEGER kdata,ns
      PARAMETER(kdata=560,ns=551)

      INTEGER i, n
      REAL x1(kdata),y1(kdata),xwork(kdata),ywork(kdata)
      
       data x1(1:ns)/219.0,219.5,220.0,220.5,221.0,221.5,222.0,
     1 222.5,223.0,223.5,224.0,224.5,225.0,225.5,
     2 226.0,226.5,227.0,227.5,228.0,228.5,229.0,
     3 229.5,230.0,230.5,231.0,231.5,232.0,232.5,
     4 233.0,233.5,234.0,234.5,235.0,235.5,236.0,
     5 236.5,237.0,237.5,238.0,238.5,239.0,239.5,
     6 240.0,240.5,241.0,241.5,242.0,242.5,243.0,
     7 243.5,244.0,244.5,245.0,245.5,246.0,246.5,
     8 247.0,247.5,248.0,248.5,249.0,249.5,250.0,
     9 250.5,251.0,251.5,252.0,252.5,253.0,253.5,
     a 254.0,254.5,255.0,255.5,256.0,256.5,257.0,
     b 257.5,258.0,258.5,259.0,259.5,260.0,260.5,
     c 261.0,261.5,262.0,262.5,263.0,263.5,264.0,
     d 264.5,265.0,265.5,266.0,266.5,267.0,267.5,
     e 268.0,268.5,269.0,269.5,270.0,270.5,271.0,
     f 271.5,272.0,272.5,273.0,273.5,274.0,274.5,
     g 275.0,275.5,276.0,276.5,277.0,277.5,278.0,
     h 278.5,279.0,279.5,280.0,280.5,281.0,281.5,
     i 282.0,282.5,283.0,283.5,284.0,284.5,285.0,
     j 285.5,286.0,286.5,287.0,287.5,288.0,288.5,
     k 289.0,289.5,290.0,290.5,291.0,291.5,292.0,
     l 292.5,293.0,293.5,294.0,294.5,295.0,295.5,
     m 296.0,296.5,297.0,297.5,298.0,298.5,299.0,
     n 299.5,300.0,300.5,301.0,301.5,302.0,302.5,
     o 303.0,303.5,304.0,304.5,305.0,305.5,306.0,
     p 306.5,307.0,307.5,308.0,308.5,309.0,309.5,
     q 310.0,310.5,311.0,311.5,312.0,312.5,313.0,
     r 313.5,314.0,314.5,315.0,315.5,316.0,316.5,
     s 317.0,317.5,318.0,318.5,319.0,319.5,320.0,
     t 320.5,321.0,321.5,322.0,322.5,323.0,323.5,
     u 324.0,324.5,325.0,325.5,326.0,326.5,327.0,
     v 327.5,328.0,328.5,329.0,329.5,330.0,330.5,
     x 331.0,331.5,332.0,332.5,333.0,333.5,334.0,
     y 334.5,335.0,335.5,336.0,336.5,337.0,337.5,
     z 338.0,338.5,339.0,339.5,340.0,340.5,341.0,
     A 341.5,342.0,342.5,343.0,343.5,344.0,344.5,
     B 345.0,345.5,346.0,346.5,347.0,347.5,348.0,
     C 348.5,349.0,349.5,350.0,350.5,351.0,351.5,
     D 352.0,352.5,353.0,353.5,354.0,354.5,355.0,
     E 355.5,356.0,356.5,357.0,357.5,358.0,358.5,
     F 359.0,359.5,360.0,360.5,361.0,361.5,362.0,
     G 362.5,363.0,363.5,364.0,364.5,365.0,365.5,
     H 366.0,366.5,367.0,367.5,368.0,368.5,369.0,
     I 369.5,370.0,370.5,371.0,371.5,372.0,372.5,
     J 373.0,373.5,374.0,374.5,375.0,375.5,376.0,
     K 376.5,377.0,377.5,378.0,378.5,379.0,379.5,
     L 380.0,380.5,381.0,381.5,382.0,382.5,383.0,
     M 383.5,384.0,384.5,385.0,385.5,386.0,386.5,
     N 387.0,387.5,388.0,388.5,389.0,389.5,390.0,
     O 390.5,391.0,391.5,392.0,392.5,393.0,393.5,
     P 394.0,394.5,395.0,395.5,396.0,396.5,397.0,
     Q 397.5,398.0,398.5,399.0,399.5,400.0,400.5,
     L 401.0,401.5,402.0,402.5,403.0,403.5,404.0,
     S 404.5,405.0,405.5,406.0,406.5,407.0,407.5,
     T 408.0,408.5,409.0,409.5,410.0,410.5,411.0,
     U 411.5,412.0,412.5,413.0,413.5,414.0,414.5,
     V 415.0,415.5,416.0,416.5,417.0,417.5,418.0,
     W 418.5,419.0,419.5,420.0,420.5,421.0,421.5,
     X 422.0,422.5,423.0,423.5,424.0,424.5,425.0,
     Y 425.5,426.0,426.5,427.0,427.5,428.0,428.5,
     Z 429.0,429.5,430.0,430.5,431.0,431.5,432.0,
     1 432.5,433.0,433.5,434.0,434.5,435.0,435.5,
     2 436.0,436.5,437.0,437.5,438.0,438.5,439.0,
     3 439.5,440.0,440.5,441.0,441.5,442.0,442.5,
     4 443.0,443.5,444.0,444.5,445.0,445.5,446.0,
     5 446.5,447.0,447.5,448.0,448.5,449.0,449.5,
     6 450.0,450.5,451.0,451.5,452.0,452.5,453.0,
     7 453.5,454.0,454.5,455.0,455.5,456.0,456.5,
     8 457.0,457.5,458.0,458.5,459.0,459.5,460.0,
     9 460.5,461.0,461.5,462.0,462.5,463.0,463.5,
     a 464.0,464.5,465.0,465.5,466.0,466.5,467.0,
     b 467.5,468.0,468.5,469.0,469.5,470.0,470.5,
     c 471.0,471.5,472.0,472.5,473.0,473.5,474.0,
     d 474.5,475.0,475.5,476.0,476.5,477.0,477.5,
     e 478.0,478.5,479.0,479.5,480.0,480.5,481.0,
     f 481.5,482.0,482.5,483.0,483.5,484.0,484.5,
     g 485.0,485.5,486.0,486.5,487.0,487.5,488.0,
     h 488.5,489.0,489.5,490.0,490.5,491.0,491.5,
     i 492.0,492.5,493.0,493.5,494.0/
       data y1(1:ns)/.984E-20,.104E-19,.106E-19,.111E-19,.115E-19,
     1 .118E-19,.122E-19,.124E-19,.126E-19,.126E-19,.125E-19,.124E-19,
     2 .125E-19,.127E-19,.127E-19,.129E-19,.131E-19,.132E-19,.135E-19,
     3 .137E-19,.140E-19,.142E-19,.148E-19,.153E-19,.157E-19,.159E-19,
     4 .161E-19,.162E-19,.161E-19,.168E-19,.174E-19,.180E-19,.184E-19,
     5 .187E-19,.189E-19,.191E-19,.193E-19,.194E-19,.196E-19,.196E-19,
     6 .201E-19,.204E-19,.208E-19,.210E-19,.214E-19,.216E-19,.219E-19,
     7 .220E-19,.223E-19,.226E-19,.228E-19,.229E-19,.230E-19,.232E-19,
     8 .233E-19,.235E-19,.238E-19,.241E-19,.246E-19,.251E-19,.257E-19,
     9 .261E-19,.265E-19,.267E-19,.269E-19,.269E-19,.271E-19,.272E-19,
     a .273E-19,.274E-19,.276E-19,.278E-19,.282E-19,.287E-19,.293E-19,
     b .298E-19,.307E-19,.312E-19,.317E-19,.321E-19,.326E-19,.328E-19,
     c .329E-19,.331E-19,.333E-19,.334E-19,.336E-19,.338E-19,.342E-19,
     d .344E-19,.348E-19,.354E-19,.359E-19,.365E-19,.373E-19,.380E-19,
     e .387E-19,.395E-19,.402E-19,.408E-19,.413E-19,.417E-19,.420E-19,
     f .422E-19,.422E-19,.422E-19,.423E-19,.424E-19,.427E-19,.429E-19,
     g .431E-19,.433E-19,.437E-19,.442E-19,.448E-19,.456E-19,.464E-19,
     h .471E-19,.478E-19,.483E-19,.487E-19,.490E-19,.492E-19,.493E-19,
     i .494E-19,.492E-19,.490E-19,.486E-19,.483E-19,.479E-19,.476E-19,
     j .472E-19,.470E-19,.468E-19,.466E-19,.465E-19,.465E-19,.468E-19,
     k .473E-19,.478E-19,.484E-19,.489E-19,.492E-19,.492E-19,.490E-19,
     l .486E-19,.481E-19,.475E-19,.470E-19,.465E-19,.458E-19,.448E-19,
     m .438E-19,.427E-19,.417E-19,.407E-19,.399E-19,.394E-19,.388E-19,
     n .382E-19,.376E-19,.372E-19,.369E-19,.368E-19,.370E-19,.372E-19,
     o .374E-19,.374E-19,.375E-19,.371E-19,.362E-19,.351E-19,.338E-19,
     p .325E-19,.315E-19,.304E-19,.292E-19,.280E-19,.271E-19,.263E-19,
     q .252E-19,.243E-19,.234E-19,.225E-19,.219E-19,.212E-19,.206E-19,
     r .202E-19,.196E-19,.192E-19,.191E-19,.188E-19,.186E-19,.185E-19,
     s .186E-19,.187E-19,.187E-19,.187E-19,.183E-19,.175E-19,.169E-19,
     t .160E-19,.150E-19,.141E-19,.134E-19,.127E-19,.121E-19,.118E-19,
     u .114E-19,.108E-19,.101E-19,.962E-20,.928E-20,.875E-20,.849E-20,
     v .821E-20,.771E-20,.738E-20,.718E-20,.686E-20,.671E-20,.663E-20,
     x .646E-20,.629E-20,.621E-20,.618E-20,.620E-20,.549E-20,.521E-20,
     y .538E-20,.535E-20,.504E-20,.494E-20,.490E-20,.452E-20,.426E-20,
     z .411E-20,.376E-20,.361E-20,.358E-20,.347E-20,.332E-20,.322E-20,
     A .310E-20,.300E-20,.294E-20,.289E-20,.286E-20,.288E-20,.288E-20,
     B .289E-20,.291E-20,.295E-20,.300E-20,.308E-20,.318E-20,.325E-20,
     C .330E-20,.339E-20,.351E-20,.363E-20,.373E-20,.385E-20,.399E-20,
     D .427E-20,.447E-20,.463E-20,.478E-20,.492E-20,.507E-20,.523E-20,
     E .539E-20,.556E-20,.577E-20,.597E-20,.615E-20,.635E-20,.656E-20,
     F .676E-20,.695E-20,.720E-20,.744E-20,.764E-20,.789E-20,.815E-20,
     G .843E-20,.871E-20,.902E-20,.933E-20,.965E-20,.100E-19,.104E-19,
     H .108E-19,.111E-19,.115E-19,.119E-19,.123E-19,.127E-19,.131E-19,
     I .135E-19,.140E-19,.144E-19,.147E-19,.151E-19,.155E-19,.159E-19,
     J .164E-19,.170E-19,.173E-19,.177E-19,.181E-19,.186E-19,.190E-19,
     K .196E-19,.202E-19,.206E-19,.210E-19,.214E-19,.218E-19,.224E-19,
     L .230E-19,.237E-19,.242E-19,.247E-19,.254E-19,.262E-19,.269E-19,
     M .279E-19,.288E-19,.296E-19,.302E-19,.310E-19,.320E-19,.329E-19,
     N .339E-19,.351E-19,.362E-19,.369E-19,.370E-19,.377E-19,.388E-19,
     O .397E-19,.403E-19,.412E-19,.422E-19,.429E-19,.430E-19,.438E-19,
     P .447E-19,.455E-19,.456E-19,.459E-19,.467E-19,.480E-19,.487E-19,
     Q .496E-19,.508E-19,.519E-19,.523E-19,.539E-19,.546E-19,.554E-19,
     L .559E-19,.577E-19,.591E-19,.599E-19,.606E-19,.620E-19,.635E-19,
     S .652E-19,.654E-19,.664E-19,.693E-19,.715E-19,.719E-19,.732E-19,
     T .758E-19,.788E-19,.797E-19,.791E-19,.811E-19,.841E-19,.853E-19,
     U .859E-19,.860E-19,.880E-19,.904E-19,.945E-19,.934E-19,.937E-19,
     V .963E-19,.971E-19,.970E-19,.965E-19,.969E-19,.989E-19,.100E-18,
     W .102E-18,.100E-18,.102E-18,.101E-18,.101E-18,.103E-18,.101E-18,
     X .104E-18,.105E-18,.106E-18,.104E-18,.105E-18,.105E-18,.101E-18,
     Y .101E-18,.105E-18,.103E-18,.102E-18,.101E-18,.977E-19,.981E-19,
     Z .100E-18,.102E-18,.989E-19,.985E-19,.104E-18,.108E-18,.105E-18,
     1 .102E-18,.964E-19,.101E-18,.106E-18,.109E-18,.104E-18,.103E-18,
     2 .107E-18,.116E-18,.109E-18,.111E-18,.981E-19,.971E-19,.106E-18,
     3 .116E-18,.108E-18,.105E-18,.970E-19,.101E-18,.104E-18,.107E-18,
     4 .102E-18,.968E-19,.100E-18,.114E-18,.113E-18,.103E-18,.974E-19,
     5 .846E-19,.870E-19,.997E-19,.101E-18,.915E-19,.941E-19,.899E-19,
     6 .110E-18,.912E-19,.856E-19,.828E-19,.615E-19,.556E-19,.647E-19,
     7 .727E-19,.575E-19,.508E-19,.438E-19,.381E-19,.361E-19,.361E-19,
     8 .313E-19,.272E-19,.244E-19,.222E-19,.182E-19,.143E-19,.132E-19,
     9 .105E-19,.895E-20,.890E-20,.794E-20,.704E-20,.646E-20,.563E-20,
     a .478E-20,.394E-20,.326E-20,.297E-20,.265E-20,.246E-20,.227E-20,
     b .208E-20,.186E-20,.176E-20,.160E-20,.144E-20,.134E-20,.120E-20,
     c .107E-20,.102E-20,.992E-21,.997E-21,.887E-21,.827E-21,.776E-21,
     d .715E-21,.671E-21,.667E-21,.610E-21,.617E-21,.554E-21,.522E-21,
     e .510E-21,.517E-21,.480E-21,.471E-21,.460E-21,.435E-21,.390E-21,
     f .371E-21,.362E-21,.352E-21,.305E-21,.305E-21,.286E-21,.253E-21,
     g .275E-21,.259E-21,.247E-21,.236E-21,.212E-21,.189E-21,.193E-21,
     h .186E-21,.182E-21,.175E-21,.174E-21,.172E-21,.166E-21,.175E-21,
     i .154E-21,.174E-21,.163E-21,.153E-21,.152E-21,.585E-22,.000E+00/

* local

      REAL yg(kw)
      REAL qy
      INTEGER ierr
      INTEGER idum
      INTEGER iw

      INTEGER mabs


**************** CH3OOH photodissociation, ROOH is same
         j = j + 1
         jlabel(j) = 'DCB2 -> RO2_R+CCO_O2+HO2+CO+R2O2+GLY+MGLY'

* mabs: Absorption cross section options:
* 2:  IUPAC97 (from  Vaghjiani and Ravishankara (1989) at 5 nm resolution).


         n = ns
	 xwork(1:n)=x1(1:n)
	 ywork(1:n)=y1(1:n)

         CALL addpnt(xwork,ywork,kdata,n,xwork(1)*(1.-deltax),0.)
         CALL addpnt(xwork,ywork,kdata,n,               0.,0.)
         CALL addpnt(xwork,ywork,kdata,n,xwork(n)*(1.+deltax),0.)
         CALL addpnt(xwork,ywork,kdata,n,           1.e+38,0.)
         CALL inter2(nw,wl,yg,n,xwork,ywork,ierr)
         IF (ierr .NE. 0) THEN
            WRITE(*,*) ierr, jlabel(j)
            STOP
         ENDIF

* quantum yield = 0.37

      qy = 0.37
      DO iw = 1, nw - 1
         DO i = 1, nz
            sq(j,i,iw) = yg(iw)*qy
         ENDDO
      ENDDO

      END

      SUBROUTINE rtlink(nz,z,
     $     iw, ag, zen, 
     $     dsdh, nid,
     $     dtrl, 
     $     dto3, 
     $     dto2, 
     $     dtso2, 
     $     dtno2, 
     $     dtcld, omcld, gcld,
     $     dtaer,omaer,gaer,
     $     edir, edn, eup, fdir, fdn, fup)
*_______________________________________________________________________

      IMPLICIT NONE

      INCLUDE 'tuv.params'

* input:

      INTEGER nz, iw
      REAL z(kz)
      REAL ag
      REAL zen
      REAL dtrl(kz,kw)
      REAL dto3(kz,kw), dto2(kz,kw), dtso2(kz,kw), dtno2(kz,kw)
      REAL dtcld(kz,kw), omcld(kz,kw), gcld(kz,kw)
      REAL dtaer(kz,kw), omaer(kz,kw), gaer(kz,kw)
      REAL dsdh(0:kz,kz)
      INTEGER nid(0:kz)


* output
      REAL edir(kz), edn(kz), eup(kz)
      REAL fdir(kz), fdn(kz), fup(kz)

* more program constants:
      REAL dr
      PARAMETER (dr = pi/180.)

* local:
      REAL dt(kz), om(kz), g(kz)
      REAL ediri(kz), edni(kz), eupi(kz)
      REAL fdiri(kz), fdni(kz), fupi(kz)
      REAL daaer, dtsct, dtabs, dsaer, dscld, dacld
      INTEGER i, ii
      LOGICAL delta

      DATA delta /.true./
*_______________________________________________________________________

* initialize:

      DO 5 i = 1, nz
         fdir(i) = 0.
         fup(i) = 0.
         fdn(i) = 0.
         edir(i) = 0.
         eup(i) = 0.
         edn(i) = 0.
 5    CONTINUE

*  set here any coefficients specific to rt scheme, 
* ----

      DO 10, i = 1, nz - 1

         dscld = dtcld(i,iw)*omcld(i,iw)
         dacld = dtcld(i,iw)*(1.-omcld(i,iw))

         dsaer = dtaer(i,iw)*omaer(i,iw)
         daaer = dtaer(i,iw)*(1.-omaer(i,iw))

         dtsct = dtrl(i,iw) + dscld + dsaer
         dtabs = dto3(i,iw) + dto2(i,iw) + dtso2(i,iw) + 
     >           dtno2(i,iw) + dacld + daaer

 	 dtabs = AMAX1(dtabs,1./largest)
 	 dtsct = AMAX1(dtsct,1./largest)

* invert z-coordinate:

         ii = nz - i
         dt(ii) = dtsct + dtabs
         om(ii) = dtsct/(dtsct + dtabs)
           IF(dtsct .EQ. 1./largest) om(ii) = 1./largest
         g(ii) = (gcld(i,iw)*dscld + gaer(i,iw)*dsaer)/dtsct
         if(g(ii).gt.1) then
	  print*,'g(ii) too high in rtlink, ii, i, nz,g(ii)=',
     1 	   ii,i,nz,g(ii)
          print*,'gcld(i,iw),dscld,gaer(i,iw),dsaer, dtsct,dtrl(i,iw)=',
     1	   gcld(i,iw),dscld,gaer(i,iw),dsaer, dtsct,dtrl(i,iw)
          print*,'dtcld(i,iw),omcld(i,iw),dtaer(i,iw),omaer(i,iw)=',
     1	     dtcld(i,iw),omcld(i,iw),dtaer(i,iw),omaer(i,iw)
          print*,'dtabs, dto3(i,iw), dto2(i,iw), dtso2(i,iw), 
     1	   dtno2(i,iw), dacld, daaer=',dtabs, dto3(i,iw), dto2(i,iw),
     2     dtso2(i,iw), dtno2(i,iw), dacld, daaer
          stop
	 endif 
   10 CONTINUE

*  call rt routine:
c      if(iw.eq.59) then
c      write(kout,*)'in rtlink g=',g
c      write(kout,*)'in rtlink om=',om
c      om(11)=0.25
c      endif
       
      CALL ps2str(nz,zen,ag,dt,om,g,
     $         dsdh, nid, delta,
     $         fdiri, fupi, fdni, ediri, eupi, edni)
      

* put on upright z-coordinate

      DO 20, i = 1, nz
         ii = nz - i + 1
         fdir(i) = fdiri(ii)
         fup(i) = fupi(ii)
         fdn(i) = fdni(ii)
c	 if(fdir(i).gt.100.or.fup(i).gt.100.or.fdn(i).gt.100) then
c	  print*,'fdir, fup, fdn =', fdir(i),fup(i),fdn(i)
c	  print*,' iz, iw, ag, nz, zen =', i, iw,ag, nz, zen
c	  print*,'nid =',nid
c	  print*,'z  =',z
c	  print*,'g =', g
c	  print*,'om=',om
c	  print*,'dt=',dt
c	  stop
c	 endif 
         edir(i) = ediri(ii)
         eup(i) = eupi(ii)
         edn(i) = edni(ii)
 20   CONTINUE
*_______________________________________________________________________

      RETURN
      END   
      SUBROUTINE schu(nz,o2col,secchi,iw,dto2,xscho2)

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Calculate the equivalent absorption cross section of O2 in the SR bands. =*
*=  The algorithm is based on:  G.Kockarts, Penetration of solar radiation   =*
*=  in the Schumann-Runge bands of molecular oxygen:  a robust approximation,=*
*=  Annales Geophysicae, v12, n12, pp. 1207ff, Dec 1994.  Calculation is     =*
*=  done on the wavelength grid used by Kockarts (1994).  Final values do    =*
*=  include effects from the Herzberg continuum.                             =*
*-----------------------------------------------------------------------------*
*=  PARAMETERS:                                                              =*
*=  NZ      - INTEGER, number of specified altitude levels in the working (I)=*
*=            grid                                                           =*
*=  O2COL   - REAL, slant overhead O2 column (molec/cc) at each specified (I)=*
*=            altitude                                                       =*
*=  IW      - INTEGER, index of current wavelength bin (Kockarts' grid)   (I)=*
*=  DTO2    - REAL, optical depth due to O2 absorption at each specified  (O)=*
*=            vertical layer at each specified wavelength                    =*
*=  XSCHO2  - REAL, molecular absorption cross section in SR bands at     (O)=*
*=            each specified wavelength.  Includes Herzberg continuum        =*
*-----------------------------------------------------------------------------*
*=  EDIT HISTORY:                                                            =*
*=  03/97  fix problem with last loop, define xscho2 at top level seperately =*
*=  10/96  converted to "true" double precision in all calculations, modified=*
*=         criterion to decide whether cross section should be zero          =*
*=  07/96  Force calculation on internal grid independent of user-defined    =*
*=         working wavelength grid                                           =*
*-----------------------------------------------------------------------------*
*= This program is free software;  you can redistribute it and/or modify     =*
*= it under the terms of the GNU General Public License as published by the  =*
*= Free Software Foundation;  either version 2 of the license, or (at your   =*
*= option) any later version.                                                =*
*= The TUV package is distributed in the hope that it will be useful, but    =*
*= WITHOUT ANY WARRANTY;  without even the implied warranty of MERCHANTIBI-  =*
*= LITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public     =*
*= License for more details.                                                 =*
*= To obtain a copy of the GNU General Public License, write to:             =*
*= Free Software Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.   =*
*-----------------------------------------------------------------------------*
*= To contact the authors, please mail to:                                   =*
*= Sasha Madronich, NCAR/ACD, P.O.Box 3000, Boulder, CO, 80307-3000, USA  or =*
*= send email to:  sasha@ucar.edu                                            =*
*-----------------------------------------------------------------------------*
*= Copyright (C) 1994,95,96  University Corporation for Atmospheric Research =*
*-----------------------------------------------------------------------------*

      IMPLICIT NONE
      INCLUDE 'tuv.params'

      REAL o2col(kz),dto2(kz,16), xscho2(kz,16)
      REAL secchi(kz)
      DOUBLE PRECISION a(16,12),b(16,12),rjm(kz),rjo2(kz)
      INTEGER iw, nz, i, j, lev


c  a(16,12)             coefficients for Rj(M) (Table 1 in Kockarts 1994)
c  b(16,12)                              Rj(O2)(Table 2 in Kockarts 1994)
c  rjm                  attenuation coefficients Rj(M)
c  rjo2                 Rj(O2)
      
      data((a(i,j),j=1,12),i=1,16)/
ca 57000-56500.5 cm-1
     l 1.13402D-01,1.00088D-20,3.48747D-01,2.76282D-20,3.47322D-01
     l,1.01267D-19
     l,1.67351D-01,5.63588D-19,2.31433D-02,1.68267D-18,0.00000D+00
     l,0.00000D+00
ca 56500-56000.5 cm-1
     l,2.55268D-03,1.64489D-21,1.85483D-01,2.03591D-21,2.60603D-01
     l,4.62276D-21
     l,2.50337D-01,1.45106D-20,1.92340D-01,7.57381D-20,1.06363D-01
     l,7.89634D-19
ca 56000-55500.5 cm-1
     l,4.21594D-03,8.46639D-22,8.91886D-02,1.12935D-21,2.21334D-01
     l,1.67868D-21
     l,2.84446D-01,3.94782D-21,2.33442D-01,1.91554D-20,1.63433D-01
     l,2.25346D-19
ca 55500-55000.5 cm-1
     l,3.93529D-03,6.79660D-22,4.46906D-02,9.00358D-22,1.33060D-01
     l,1.55952D-21
     l,3.25506D-01,3.43763D-21,2.79405D-01,1.62086D-20,2.10316D-01
     l,1.53883D-19
ca 55000-54500.5 cm-1
     l,2.60939D-03,2.33791D-22,2.08101D-02,3.21734D-22,1.67186D-01
     l,5.77191D-22
     l,2.80694D-01,1.33362D-21,3.26867D-01,6.10533D-21,1.96539D-01
     l,7.83142D-20
ca 54500-54000.5 cm-1
     l,9.33711D-03,1.32897D-22,3.63980D-02,1.78786D-22,1.46182D-01 
     l,3.38285D-22
     l,3.81762D-01,8.93773D-22,2.58549D-01,4.28115D-21,1.64773D-01
     l,4.67537D-20
ca 54000-53500.5 cm-1
     l,9.51799D-03,1.00252D-22,3.26320D-02,1.33766D-22,1.45962D-01
     l,2.64831D-22
     l,4.49823D-01,6.42879D-22,2.14207D-01,3.19594D-21,1.45616D-01
     l,2.77182D-20
ca 53500-53000.5 cm-1
     l,7.87331D-03,3.38291D-23,6.91451D-02,4.77708D-23,1.29786D-01
     l,8.30805D-23
     l,3.05103D-01,2.36167D-22,3.35007D-01,8.59109D-22,1.49766D-01
     l,9.63516D-21
ca 53000-52500.5 cm-1
     l,6.92175D-02,1.56323D-23,1.44403D-01,3.03795D-23,2.94489D-01 
     l,1.13219D-22
     l,3.34773D-01,3.48121D-22,9.73632D-02,2.10693D-21,5.94308D-02 
     l,1.26195D-20
ca 52500-52000.5 cm-1
     l,1.47873D-01,8.62033D-24,3.15881D-01,3.51859D-23,4.08077D-01 
     l,1.90524D-22
     l,8.08029D-02,9.93062D-22,3.90399D-02,6.38738D-21,8.13330D-03
     l,9.93644D-22
ca 52000-51500.5 cm-1
     l,1.50269D-01,1.02621D-23,2.39823D-01,3.48120D-23,3.56408D-01
     l,1.69494D-22
     l,1.61277D-01,6.59294D-22,8.89713D-02,2.94571D-21,3.25063D-03
     l,1.25548D-20
ca 51500-51000.5 cm-1
     l,2.55746D-01,8.49877D-24,2.94733D-01,2.06878D-23,2.86382D-01 
     l,9.30992D-23
     l,1.21011D-01,3.66239D-22,4.21105D-02,1.75700D-21,0.00000D+00
     l,0.00000D+00
ca 51000-50500.5 cm-1
     l,5.40111D-01,7.36085D-24,2.93263D-01,2.46742D-23,1.63417D-01
     l,1.37832D-22
     l,3.23781D-03,2.15052D-21,0.00000D+00,0.00000D+00,0.00000D+00
     l,0.00000D+00
ca 50500-50000.5 cm-1
     l,8.18514D-01,7.17937D-24,1.82262D-01,4.17496D-23,0.00000D+00 
     l,0.00000D+00
     l,0.00000D+00,0.00000D+00,0.00000D+00,0.00000D+00,0.00000D+00
     l,0.00000D+00
ca 50000-49500.5 cm-1
     l,8.73680D-01,7.13444D-24,1.25583D-01,2.77819D-23,0.00000D+00 
     l,0.00000D+00
     l,0.00000D+00,0.00000D+00,0.00000D+00,0.00000D+00,0.00000D+00 
     l,0.00000D+00
ca 49500-49000.5 cm-1
     l,3.32476D-04,7.00362D-24,9.89000D-01,6.99600D-24,0.00000D+00
     l,0.00000D+00
     l,0.00000D+00,0.00000D+00,0.00000D+00,0.00000D+00,0.00000D+00
     l,0.00000D+00/


      data((b(i,j),j=1,12),i=1,16)/
c  57000-56500.5 cm-1
     l 1.07382D-21,9.95029D-21,7.19430D-21,2.48960D-20,2.53735D-20
     l,7.54467D-20
     l,4.48987D-20,2.79981D-19,9.72535D-20,9.29745D-19,2.30892D-20
     l,4.08009D-17
c  56500-56000.5 cm-1
     l,3.16903D-22,1.98251D-21,5.87326D-22,3.44057D-21,2.53094D-21
     l,8.81484D-21
     l,8.82299D-21,4.17179D-20,2.64703D-20,2.43792D-19,8.73831D-20
     l,1.46371D-18
c  56000-55500.5 cm-1
     l,1.64421D-23,9.26011D-22,2.73137D-22,1.33640D-21,9.79188D-22
     l,2.99706D-21
     l,3.37768D-21,1.39438D-20,1.47898D-20,1.04322D-19,4.08014D-20
     l,6.31023D-19
c  55500-55000.5 cm-1
     l,8.68729D-24,7.31056D-22,8.78313D-23,1.07173D-21,8.28170D-22
     l,2.54986D-21
     l,2.57643D-21,9.42698D-21,9.92377D-21,5.21402D-20,3.34301D-20
     l,2.91785D-19
c  55000-54500.5 cm-1
     l,1.20679D-24,2.44092D-22,2.64326D-23,4.03998D-22,2.53514D-22
     l,8.53166D-22
     l,1.29834D-21,3.74482D-21,5.12103D-21,2.65798D-20,2.10948D-20
     l,2.35315D-19
c  54500-54000.5 cm-1
     l,2.79656D-24,1.40820D-22,3.60824D-23,2.69510D-22,4.02850D-22
     l,8.83735D-22
     l,1.77198D-21,6.60221D-21,9.60992D-21,8.13558D-20,4.95591D-21
     l,1.22858D-17
c  54000-53500.5 cm-1
     l,2.36959D-24,1.07535D-22,2.83333D-23,2.16789D-22,3.35242D-22
     l,6.42753D-22
     l,1.26395D-21,5.43183D-21,4.88083D-21,5.42670D-20,3.27481D-21
     l,1.58264D-17
c  53500-53000.5 cm-1
     l,8.65018D-25,3.70310D-23,1.04351D-23,6.43574D-23,1.17431D-22
     l,2.70904D-22
     l,4.88705D-22,1.65505D-21,2.19776D-21,2.71172D-20,2.65257D-21
     l,2.13945D-17
c  53000-52500.5 cm-1
     l,9.63263D-25,1.54249D-23,4.78065D-24,2.97642D-23,6.40637D-23
     l,1.46464D-22
     l,1.82634D-22,7.12786D-22,1.64805D-21,2.37376D-17,9.33059D-22
     l,1.13741D-20
c  52500-52000.5 cm-1
     l,1.08414D-24,8.37560D-24,9.15550D-24,2.99295D-23,9.38405D-23
     l,1.95845D-22
     l,2.84356D-22,3.39699D-21,1.94524D-22,2.72227D-19,1.18924D-21
     l,3.20246D-17
c  52000-51500.5 cm-1
     l,1.52817D-24,1.01885D-23,1.22946D-23,4.16517D-23,9.01287D-23 
     l,2.34869D-22
     l,1.93510D-22,1.44956D-21,1.81051D-22,5.17773D-21,9.82059D-22
     l,6.22768D-17
c  51500-51000.5 cm-1
     l,2.12813D-24,8.48035D-24,5.23338D-24,1.93052D-23,1.99464D-23 
     l,7.48997D-23
     l,4.96642D-22,6.15691D-17,4.47504D-23,2.76004D-22,8.26788D-23
     l,1.65278D-21
c  51000-50500.5 cm-1
     l,3.81336D-24,7.32307D-24,5.60549D-24,2.04651D-23,3.36883D-22
     l,6.15708D-17
     l,2.09877D-23,1.07474D-22,9.13562D-24,8.41252D-22,0.00000D+00
     l,0.00000D+00
c  50500-50000.5 cm-1
     l,5.75373D-24,7.15986D-24,5.90031D-24,3.05375D-23,2.97196D-22
     l,8.92000D-17
     l,8.55920D-24,1.66709D-17,0.00000D+00,0.00000D+00,0.00000D+00
     l,0.00000D+00
c  50000-49500.5 cm-1
     l,6.21281D-24,7.13108D-24,3.30780D-24,2.61196D-23,1.30783D-22 
     l,9.42550D-17
     l,2.69241D-24,1.46500D-17,0.00000D+00,0.00000D+00,0.00000D+00 
     l,0.00000D+00
c  49500-49000.5 cm-1
     l,6.81118D-24,6.98767D-24,7.55667D-25,2.75124D-23,1.94044D-22 
     l,1.45019D-16
     l,1.92236D-24,3.73223D-17,0.00000D+00,0.00000D+00,0.00000D+00 
     l,0.00000D+00/

c initialize R(M)

      DO lev = 1, nz
        rjm(lev) = 0.D+00
        rjo2(lev) = 0.D+00
      END DO

c calculate sum of exponentials (eqs 7 and 8 of Kockarts 1994)

      DO j = 1, 11, 2
        DO lev = 1, nz
           rjm(lev) = rjm(lev) + a(iw,j)*DEXP(-a(iw,j+1)*
     >                           DBLE(o2col(lev)))
           rjo2(lev) = rjo2(lev) + b(iw,j)*DEXP(-b(iw,j+1)*
     >                           DBLE(o2col(lev)))
        END DO
      END DO

      DO lev = 1, nz-1

         IF (rjm(lev) .GT. 1.D-100) THEN
            IF (rjo2(lev) .GT. 1.D-100) THEN
               xscho2(lev,iw) = rjo2(lev)/rjm(lev)
            ELSE
               xscho2(lev,iw) = 0.
            ENDIF
            IF (rjm(lev+1) .GT. 0.) THEN
               dto2(lev,iw) = LOG(rjm(lev+1)) / secchi(lev+1)
     $              - LOG(rjm(lev)) * secchi(lev)
            ELSE
               dto2(lev,iw) = 1000.
            ENDIF
         ELSE
            xscho2(lev,iw) = 0.
            dto2(lev,iw) = 1000.
         ENDIF
      END DO

      dto2(nz,iw) = 0.
      IF(rjm(nz) .GT. 1.D-100) THEN
         xscho2(nz,iw) = rjo2(nz)/rjm(nz)
      ELSE
         xscho2(nz,iw) = 0.
      ENDIF

      RETURN
      END



      SUBROUTINE setaer(nz,z,nw,wl,dtaer,omaer,gaer,mss,rhu,nzstem,daod)

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Set up an altitude profile of aerosols, and corresponding absorption     =*
*=  optical depths, single scattering albedo, and asymmetry factor.          =*
*=  Single scattering albedo and asymmetry factor can be selected for each   =*
*=  input aerosol layer (do not have to correspond to working altitude       =*
*=  grid).  See loop 27.                                                     =*
*-----------------------------------------------------------------------------*
*=  PARAMETERS:                                                              =*
*=  NZ      - INTEGER, number of specified altitude levels in the working (I)=*
*=            grid                                                           =*
*=  Z       - REAL, specified altitude working grid (km)                  (I)=*
*=  NW     - INTEGER, number of specified intervals + 1 in working        (I)=*
*=           wavelength grid                                                 =*
*=  WL     - REAL, vector of lower limits of wavelength intervals in      (I)=*
*=           working wavelength grid                                         =*
*=  DTAER   - REAL, optical depth due to absorption by aerosols at each   (O)=*
*=            altitude and wavelength                                        =*
*=  OMAER   - REAL, single scattering albedo due to aerosols at each      (O)=*
*=            defined altitude and wavelength                                =*
*=  GAER    - REAL, aerosol asymmetry factor at each defined altitude and (O)=*
*=            wavelength                                                     =*
*=  nzstem  - input data layers                                              =*
*=  mss     - inputted aerosol concentrations in molecular/cm3               =*   
*=  rhu     - inputted relative humidity for each layer                      =* 
*-----------------------------------------------------------------------------*
*=  EDIT HISTORY:                                                            =*
*=  Original                                                                 =*
*-----------------------------------------------------------------------------*
*= This program is free software;  you can redistribute it and/or modify     =*
*= it under the terms of the GNU General Public License as published by the  =*
*= Free Software Foundation;  either version 2 of the license, or (at your   =*
*= option) any later version.                                                =*
*= The TUV package is distributed in the hope that it will be useful, but    =*
*= WITHOUT ANY WARRANTY;  without even the implied warranty of MERCHANTIBI-  =*
*= LITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public     =*
*= License for more details.                                                 =*
*= To obtain a copy of the GNU General Public License, write to:             =*
*= Free Software Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.   =*
*-----------------------------------------------------------------------------*
*= To contact the authors, please mail to:                                   =*
*= Sasha Madronich, NCAR/ACD, P.O.Box 3000, Boulder, CO, 80307-3000, USA  or =*
*= send email to:  sasha@ucar.edu                                            =*
*-----------------------------------------------------------------------------*
*= Copyright (C) 1994,95,96  University Corporation for Atmospheric Research =*
*-----------------------------------------------------------------------------*

      IMPLICIT NONE
      INCLUDE 'tuv.params'

      INTEGER kdata,nzstem,naetp
      PARAMETER(kdata=51,naetp=5)
                    ! number of aerosol types 1=dust,2=watersoluble(50%sulfate),
     	            !  3=black carbon, 4=sea salt, 5=organic carbon
* input: (grids)
      REAL wl(kw)
      REAL z(kz),mss(nzstem,naetp),rhu(nzstem)
      INTEGER nz
      INTEGER nw

* output: (on converted grid)
      REAL dtaer(kz,kw), omaer(kz,kw), gaer(kz,kw), daod(nzstem,naetp+1)

* local:
      REAL zd(kdata), aer(kdata), x1(kdata), y1(kdata)
      real dtaerwl(kz,kw,naetp+1),omaerwl(kz,kw,naetp+1),absorp(kz,kw),
     1   gaerwl(kz,kw,naetp+1)       ! aerosol optical depth et al for each wavelength, each species
      REAL cd(kdata), omd(kdata), gd(kdata)
      REAL womd(kdata), wgd(kdata)

      REAL cz(kz)
      REAL omz(kz)
      REAL gz(kz)

      REAL wc, wscale
      INTEGER i, iw, nd

      LOGICAL aerosl

      REAL fsum
      EXTERNAL fsum
      
      INTEGER  ierr,n 
      INTEGER  IATP, ILAY, IRH  !loop control


      REAL     RH00(8),  RHUM
      REAL     WAVEI(7)
      
      REAL     Naer(nzstem)
      REAL     NDS(nzstem,NAETP), CONV(NAETP)
      
      REAL     SBEA(6,8,NAETP), SBAA(6,8,NAETP), SGA(6,8,NAETP)
      REAL     ABEAR(7,NAETP), AGAR(7,NAETP), AOMR(7,NAETP)
      REAL     extaer(7), albaer(7), asmaer(7)
      REAL     BEA(KW), OMA(KW), GA(KW)

*_______________________________________________________________________

* Aerosol data from Elterman (1968)
* These are vertical optical depths per km, in 1 km
* intervals from 0 km to 50 km, at 340 nm.
* This is one option.  User can specify different data set.

      DATA aer/
     1     2.40E-01,1.06E-01,4.56E-02,1.91E-02,1.01E-02,7.63E-03,
     2     5.38E-03,5.00E-03,5.15E-03,4.94E-03,4.82E-03,4.51E-03,
     3     4.74E-03,4.37E-03,4.28E-03,4.03E-03,3.83E-03,3.78E-03,
     4     3.88E-03,3.08E-03,2.26E-03,1.64E-03,1.23E-03,9.45E-04,
     5     7.49E-04,6.30E-04,5.50E-04,4.21E-04,3.22E-04,2.48E-04,
     6     1.90E-04,1.45E-04,1.11E-04,8.51E-05,6.52E-05,5.00E-05,
     7     3.83E-05,2.93E-05,2.25E-05,1.72E-05,1.32E-05,1.01E-05,
     8     7.72E-06,5.91E-06,4.53E-06,3.46E-06,2.66E-06,2.04E-06,
     9     1.56E-06,1.19E-06,9.14E-07/
      
      DATA WAVEI /.185E+3,.25E+3,.3E+3,
     >               .4E+3,.55E+3,.7E+3,1.5E+3/
      DATA RH00  /.0,.50,.70,.80,.90,.95,.98,.99/

C---------------------------------------------------------------------------
C** CONV convert ug/m3 to #/cm3 for dust 1, sulfate 2, soot 3, sea salt 4 **
C** aerosol size distributions follow OPAC 
C** organic carbon 5, optical from Liousse et al., 1996 & Cooke et al., 1999
C** r=0.0212micron+-2.24, density=1.8g/cm3, optcial reliable between
C** 300~1060nm wavelength
C** dust    300.142 particles/cm3 / 221.8ug/m3 = 1.3532    (#/cm3)/(ug/m3)
C** sulfate 28000   particles/cm3 / 56.0 ug/m3 = 500.0     (#/cm3)/(ug/m3)
C** soot    130,000 particles/cm3 / 7.8  ug/m3 = 16666.667 (#/cm3)/(ug/m3)
C** seasalt 20.0032 particles/cm3 / 39.5 ug/m3 = 0.50641   (#/cm3)/(ug/m3)
C** organic carbon  no reliable source, use WaterSoluble data in OPAC
C**                 1.34E-3(micro-g/m3/part/cm3)=746.27    (#/cm3)/(ug/m3)

c      DATA CONV  / 1.3532E+00, 5.0000E+02, 1.6667E+04,        ! for unit ug/m3
c     >             5.0641E-01, 7.4627E+02/
      DATA CONV  /2.6961E-11, 7.9695E-08, 3.3206E-07,       ! for unit molecular/cm3
     >             4.914E-11, 1.4868E-08/                    ! assuming dust has molecular weight 12         

C---------------------------------------------------------------------------
C** SBEA extinction coefficient (1/km) normalized to 1 particle/cm3
      DATA SBEA  /
C** mineral dust                                 **
     >  4.153E-04,4.177E-04,4.190E-04,4.182E-04,4.180E-04,4.032E-04,
     >  4.153E-04,4.177E-04,4.190E-04,4.182E-04,4.180E-04,4.032E-04,
     >  4.153E-04,4.177E-04,4.190E-04,4.182E-04,4.180E-04,4.032E-04,
     >  4.153E-04,4.177E-04,4.190E-04,4.182E-04,4.180E-04,4.032E-04,
     >  4.153E-04,4.177E-04,4.190E-04,4.182E-04,4.180E-04,4.032E-04,
     >  4.153E-04,4.177E-04,4.190E-04,4.182E-04,4.180E-04,4.032E-04,
     >  4.153E-04,4.177E-04,4.190E-04,4.182E-04,4.180E-04,4.032E-04,
     >  4.153E-04,4.177E-04,4.190E-04,4.182E-04,4.180E-04,4.032E-04,
C** watersoluable particles  (50% sulfate)       **
     >  9.520E-06,8.182E-06,6.011E-06,3.905E-06,2.649E-06,5.937E-07,
     >  1.511E-05,1.301E-05,9.654E-06,6.366E-06,4.382E-06,9.935E-07,
     >  1.835E-05,1.587E-05,1.188E-05,7.913E-06,5.493E-06,1.264E-06,
     >  2.169E-05,1.885E-05,1.423E-05,9.584E-06,6.711E-06,1.571E-06,
     >  2.945E-05,2.589E-05,1.994E-05,1.374E-05,9.794E-06,2.392E-06,
     >  4.076E-05,3.637E-05,2.874E-05,2.038E-05,1.485E-05,3.846E-06,
     >  6.157E-05,5.616E-05,4.605E-05,3.404E-05,2.560E-05,7.265E-06,
     >  7.993E-05,7.410E-05,6.230E-05,4.740E-05,3.644E-05,1.105E-05,
C** black carbon                                 **
     >  2.2035E-06,1.97E-06,13.995E-07,9.0E-07,6.375E-07,2.491E-07,
     >  2.2035E-06,1.97E-06,13.995E-07,9.0E-07,6.375E-07,2.491E-07,
     >  2.2035E-06,1.97E-06,13.995E-07,9.0E-07,6.375E-07,2.491E-07,
     >  2.2035E-06,1.97E-06,13.995E-07,9.0E-07,6.375E-07,2.491E-07,
     >  2.2035E-06,1.97E-06,13.995E-07,9.0E-07,6.375E-07,2.491E-07,
     >  2.2035E-06,1.97E-06,13.995E-07,9.0E-07,6.375E-07,2.491E-07,
     >  2.2035E-06,1.97E-06,13.995E-07,9.0E-07,6.375E-07,2.491E-07,
     >  2.2035E-06,1.97E-06,13.995E-07,9.0E-07,6.375E-07,2.491E-07,
C** sea salt
     >  9.132E-04,9.455E-04,9.934E-04,1.037E-03,1.040E-03,7.527E-04,
     >  2.255E-03,2.304E-03,2.412E-03,2.536E-03,2.602E-03,2.213E-03,
     >  2.814E-03,2.880E-03,3.005E-03,3.164E-03,3.261E-03,2.919E-03,
     >  3.378E-03,3.445E-03,3.593E-03,3.777E-03,3.908E-03,3.643E-03,
     >  4.757E-03,4.825E-03,5.006E-03,5.252E-03,5.459E-03,5.441E-03,
     >  6.929E-03,7.044E-03,7.249E-03,7.560E-03,7.840E-03,8.297E-03,
     >  1.199E-02,1.214E-02,1.238E-02,1.280E-02,1.318E-02,1.462E-02,
     >  1.832E-02,1.849E-02,1.878E-02,1.930E-02,1.975E-02,2.214E-02,
C** organic carbon from Liousse et al., 1996, RH dependency follow
C** watersoluable particles.
     >  1.405E-05,11.72E-06,8.969E-06,6.700E-06,5.259E-06,2.541E-06,                    
     >  2.228E-05,1.861E-05,1.424E-05,10.63E-06,8.346E-06,4.033E-06,                    
     >  2.707E-05,2.260E-05,1.728E-05,1.291E-05,10.14E-06,4.897E-06,                    
     >  3.200E-05,2.138E-05,2.044E-05,1.526E-05,11.98E-06,5.790E-06,                    
     >  4.345E-05,3.627E-05,2.775E-05,2.072E-05,1.626E-05,7.861E-06,                    
     >  6.014E-05,5.020E-05,3.840E-05,2.869E-05,2.251E-05,10.88E-06,                    
     >  9.083E-05,7.582E-05,5.801E-05,4.333E-05,3.401E-05,1.643E-05,                    
     >  11.79E-05,9.844E-05,7.530E-05,5.625E-05,4.415E-05,2.133E-05/                    
     
      DATA SBAA  /
C** mineral dust                                 **
     >  6.234E-01,6.539E-01,7.549E-01,8.732E-01,9.125E-01,9.313E-01,
     >  6.234E-01,6.539E-01,7.549E-01,8.732E-01,9.125E-01,9.313E-01,
     >  6.234E-01,6.539E-01,7.549E-01,8.732E-01,9.125E-01,9.313E-01,
     >  6.234E-01,6.539E-01,7.549E-01,8.732E-01,9.125E-01,9.313E-01,
     >  6.234E-01,6.539E-01,7.549E-01,8.732E-01,9.125E-01,9.313E-01,
     >  6.234E-01,6.539E-01,7.549E-01,8.732E-01,9.125E-01,9.313E-01,
     >  6.234E-01,6.539E-01,7.549E-01,8.732E-01,9.125E-01,9.313E-01,
     >  6.234E-01,6.539E-01,7.549E-01,8.732E-01,9.125E-01,9.313E-01,
C** watersoluable particles  (50% sulfate)       **
     >  8.360E-01,9.493E-01,9.685E-01,9.615E-01,9.522E-01,7.631E-01,
     >  8.918E-01,9.687E-01,9.808E-01,9.765E-01,9.708E-01,8.479E-01,
     >  9.094E-01,9.744E-01,9.844E-01,9.811E-01,9.765E-01,8.773E-01,
     >  9.223E-01,9.785E-01,9.871E-01,9.843E-01,9.806E-01,8.990E-01,
     >  9.415E-01,9.842E-01,9.907E-01,9.890E-01,9.865E-01,9.308E-01,
     >  9.568E-01,9.888E-01,9.935E-01,9.925E-01,9.910E-01,9.548E-01,
     >  9.707E-01,9.926E-01,9.959E-01,9.954E-01,9.946E-01,9.741E-01,
     >  9.773E-01,9.944E-01,9.970E-01,9.967E-01,9.962E-01,9.819E-01,
C** black carbon                                 **
     >  3.081E-01,3.128E-01,2.672E-01,2.088E-01,1.616E-01,4.407E-02,
     >  3.081E-01,3.128E-01,2.672E-01,2.088E-01,1.616E-01,4.407E-02,
     >  3.081E-01,3.128E-01,2.672E-01,2.088E-01,1.616E-01,4.407E-02,
     >  3.081E-01,3.128E-01,2.672E-01,2.088E-01,1.616E-01,4.407E-02,
     >  3.081E-01,3.128E-01,2.672E-01,2.088E-01,1.616E-01,4.407E-02,
     >  3.081E-01,3.128E-01,2.672E-01,2.088E-01,1.616E-01,4.407E-02,
     >  3.081E-01,3.128E-01,2.672E-01,2.088E-01,1.616E-01,4.407E-02,
     >  3.081E-01,3.128E-01,2.672E-01,2.088E-01,1.616E-01,4.407E-02,
C ** sea salt
     >  9.998E-01,9.999E-01,1.000E+00,1.000E+00,1.000E+00,9.956E-01,
     >  9.999E-01,1.000E+00,1.000E+00,1.000E+00,1.000E+00,9.970E-01,
     >  9.999E-01,1.000E+00,1.000E+00,1.000E+00,1.000E+00,9.971E-01,
     >  1.000E+00,1.000E+00,1.000E+00,1.000E+00,1.000E+00,9.971E-01,
     >  1.000E+00,1.000E+00,1.000E+00,1.000E+00,1.000E+00,9.969E-01,
     >  1.000E+00,1.000E+00,1.000E+00,1.000E+00,1.000E+00,9.966E-01,
     >  1.000E+00,1.000E+00,1.000E+00,1.000E+00,1.000E+00,9.957E-01,
     >  1.000E+00,1.000E+00,1.000E+00,1.000E+00,1.000E+00,9.947E-01,
C ** organic carbon from Cooke et al., 1999 only 550nm is available,
C ** at other wavelength and RH, is the same as 550nm, be careful.
     >  9.800E-01,9.800E-01,9.800E-01,9.800E-01,9.800E-01,9.800E-01,
     >  9.800E-01,9.800E-01,9.800E-01,9.800E-01,9.800E-01,9.800E-01,
     >  9.800E-01,9.800E-01,9.800E-01,9.800E-01,9.800E-01,9.800E-01,
     >  9.800E-01,9.800E-01,9.800E-01,9.800E-01,9.800E-01,9.800E-01,
     >  9.800E-01,9.800E-01,9.800E-01,9.800E-01,9.800E-01,9.800E-01,
     >  9.800E-01,9.800E-01,9.800E-01,9.800E-01,9.800E-01,9.800E-01,
     >  9.800E-01,9.800E-01,9.800E-01,9.800E-01,9.800E-01,9.800E-01,
     >  9.800E-01,9.800E-01,9.800E-01,9.800E-01,9.800E-01,9.800E-01/
     
     
      DATA SGA   /
C** mineral dust                                 **
     >  8.572E-01,8.363E-01,7.849E-01,7.329E-01,7.112E-01,6.933E-01,
     >  8.572E-01,8.363E-01,7.849E-01,7.329E-01,7.112E-01,6.933E-01,
     >  8.572E-01,8.363E-01,7.849E-01,7.329E-01,7.112E-01,6.933E-01,
     >  8.572E-01,8.363E-01,7.849E-01,7.329E-01,7.112E-01,6.933E-01,
     >  8.572E-01,8.363E-01,7.849E-01,7.329E-01,7.112E-01,6.933E-01,
     >  8.572E-01,8.363E-01,7.849E-01,7.329E-01,7.112E-01,6.933E-01,
     >  8.572E-01,8.363E-01,7.849E-01,7.329E-01,7.112E-01,6.933E-01,
     >  8.572E-01,8.363E-01,7.849E-01,7.329E-01,7.112E-01,6.933E-01,
C** watersoluable particles  (50% sulfate)       **
     >  6.900E-01,6.590E-01,6.390E-01,6.140E-01,5.890E-01,4.850E-01,
     >  7.290E-01,7.110E-01,6.950E-01,6.720E-01,6.470E-01,5.380E-01,
     >  7.400E-01,7.260E-01,7.120E-01,6.900E-01,6.660E-01,5.590E-01,
     >  7.480E-01,7.360E-01,7.240E-01,7.040E-01,6.810E-01,5.750E-01,
     >  7.590E-01,7.520E-01,7.430E-01,7.250E-01,7.040E-01,6.030E-01,
     >  7.670E-01,7.640E-01,7.580E-01,7.430E-01,7.240E-01,6.310E-01,
     >  7.730E-01,7.740E-01,7.720E-01,7.610E-01,7.450E-01,6.630E-01,
     >  7.760E-01,7.780E-01,7.780E-01,7.700E-01,7.560E-01,6.820E-01,
C** black carbon                                 **
     >  5.020E-01,4.530E-01,3.960E-01,3.360E-01,2.910E-01,1.580E-01,
     >  5.020E-01,4.530E-01,3.960E-01,3.360E-01,2.910E-01,1.580E-01,
     >  5.020E-01,4.530E-01,3.960E-01,3.360E-01,2.910E-01,1.580E-01,
     >  5.020E-01,4.530E-01,3.960E-01,3.360E-01,2.910E-01,1.580E-01,
     >  5.020E-01,4.530E-01,3.960E-01,3.360E-01,2.910E-01,1.580E-01,
     >  5.020E-01,4.530E-01,3.960E-01,3.360E-01,2.910E-01,1.580E-01,
     >  5.020E-01,4.530E-01,3.960E-01,3.360E-01,2.910E-01,1.580E-01,
     >  5.020E-01,4.530E-01,3.960E-01,3.360E-01,2.910E-01,1.580E-01,
C** sea salt
     >  7.219E-01,7.070E-01,6.990E-01,6.919E-01,6.968E-01,7.058E-01,
     >  7.856E-01,7.826E-01,7.757E-01,7.717E-01,7.736E-01,7.844E-01,
     >  7.985E-01,7.906E-01,7.837E-01,7.787E-01,7.806E-01,7.924E-01,
     >  8.034E-01,7.996E-01,7.896E-01,7.847E-01,7.856E-01,7.974E-01,
     >  8.163E-01,8.124E-01,8.006E-01,7.936E-01,7.906E-01,8.025E-01,
     >  8.272E-01,8.223E-01,8.115E-01,8.016E-01,7.976E-01,8.055E-01,
     >  8.380E-01,8.361E-01,8.283E-01,8.146E-01,8.066E-01,8.056E-01,
     >  8.447E-01,8.439E-01,8.382E-01,8.255E-01,8.156E-01,8.056E-01,
C** organic carbon no assymmetry factor is unavailable, use g value
C** for continental clean at 550nm from OPAC, then independent of
C** wavlength and RH, can not be used for radiative transfer calculation.
     >  7.090E-01,7.090E-01,7.090E-01,7.090E-01,7.090E-01,7.090E-01,
     >  7.090E-01,7.090E-01,7.090E-01,7.090E-01,7.090E-01,7.090E-01,
     >  7.090E-01,7.090E-01,7.090E-01,7.090E-01,7.090E-01,7.090E-01,
     >  7.090E-01,7.090E-01,7.090E-01,7.090E-01,7.090E-01,7.090E-01,
     >  7.090E-01,7.090E-01,7.090E-01,7.090E-01,7.090E-01,7.090E-01,
     >  7.090E-01,7.090E-01,7.090E-01,7.090E-01,7.090E-01,7.090E-01,
     >  7.090E-01,7.090E-01,7.090E-01,7.090E-01,7.090E-01,7.090E-01,
     >  7.090E-01,7.090E-01,7.090E-01,7.090E-01,7.090E-01,7.090E-01/
      
C-------------------------------------------------------------------------C
C   SBEA(6,8,5) extinction coefficient                                    C
C   SBAA(6,8,5) single scattering albedo                                  C
C   SGA (6,8,5) asymmetry factor                                          C
C        6 wavelength, .25, .3, .4, .55, .7, 1.5 micron                   C
C        8 relative humidity, .0, .5, .7, .8, .9, .95, .98, .99           C
C        5 aerosol type, 1-dust, 2-sulfate, 3-black carbon, 4-sea salt    C
C                        5-organice carbon                                C
C-------------------------------------------------------------------------C
C   optical parameters are calculated using OPAC                          **
C ** 1) mineral 300.142 particles/cm3, 221.8 ug/m3                        **
C **    number density: nuc 269.5, acc 30.5, coa 0.142 /cm3               **
C **    mass %: nuc 7.5 (3.4%), acc 168.7 (76.1%), coa 45.6 (20.5%) ug/m3 **
C ** 2) water soluble  28000 particles/cm3, 56.0 ug/m3                    **
C ** 3) black carbon 130,000 particles/cm3,  7.8 ug/m3                    **
C ** 4) sea salt 20.0032 part/cm3, 39.5ug/m3                              **
C **    number density: acc 20.(99.98%), coa 3.2E-3(0.02%)                **
C **    mass : acc 38.6 (97.72%), coa 0.9 (2.28%)                         ** 
C ** 5) organic carbon as Watersoluable in OPAC,1.34ug/m3/part/cm3        **
C-------------------------------------------------------------------------C
     
*_______________________________________________________________________

* initialize

      DO 15, iw = 1, nw - 1
         DO 10, i = 1, nz - 1
            dtaer(i,iw) = 0.
            omaer(i,iw) = 1.
            gaer(i,iw) = 0
   10    CONTINUE
   15 CONTINUE

* if dont want any aerosols, set AEROSL = .FALSE.

C      aerosl = .FALSE.
      aerosl = .TRUE.
      IF (.NOT. aerosl) THEN
c         WRITE(kout,*) 'no aerosols'
         RETURN
      ENDIF

* Altitudes corresponding to Elterman profile, from bottom to top:

c      WRITE(kout,*)'aerosols:  Elterman (1968)'
      nd = 51
      DO 22, i = 1, nd
         zd(i) = FLOAT(i-1)
   22 CONTINUE

* assume these are point values (at each level), so find column
* increments

      DO 27, i = 1, nd - 1
         cd(i) = (aer(i+1) + aer(i)) / 2.
         omd(i) = .99
         gd(i) = 0.6
   27 CONTINUE

*********** end data input.

* Compute integrals and averages over grid layers:
* for g and omega, use averages weigthed by optical depth

      DO 29, i = 1, nd-1
       womd(i) = omd(i) * cd(i)
       wgd(i) = gd(i)  * cd(i)
   29 CONTINUE
      CALL inter3(nz,z,cz, nd,zd,cd, 1)
      CALL inter3(nz,z,omz, nd, zd,womd, 1)
      CALL inter3(nz,z,gz , nd, zd,wgd, 1)
      DO 30, i = 1, nz-1
         IF (cz(i) .GT. 0.) THEN
            omz(i) = omz(i)/cz(i)
            gz(i)  = gz(i) /cz(i)
         ELSE
            omz(i) = 1.
            gz(i) = 0.
         ENDIF
   30 CONTINUE

* assign at all wavelengths
* (can move wavelength loop outside if want to vary with wavelength)

      DO 50, iw = 1, nw - 1
         wc = (wl(iw)+wl(iw+1))/2.

* Elterman's data are for 340 nm, so assume optical depth scales 
* inversely with first power of wavelength.

         wscale = 340./wc

* optical depths:

         DO 40, i = 1, nz - 1
            dtaer(i,iw) = cz(i)  * wscale
            omaer(i,iw) = omz(i)
            gaer(i,iw) = gz(i)
   40    CONTINUE

   50 CONTINUE


*________end of initialization

C--- transfer aerosol concentration ---------------------------------------C
      DO ILAY= 1, nzstem                   !vertical layer loop

       Naer(ILAY)      = 0.               !aerosol total number density
      
       DO IATP=1, NAETP                   !aerosol  loop
C        MSS(ILAY,IATP) = 0.               !mass concentration
        NDS(ILAY,IATP) = 0.               !number density
       END DO                             !aerosol  loop END

C--- aerosol mass concentration input here in unit of micro-g/m3 ---------C
C    IATP=1, dust;          IATP=2, watersoluble (50% sulfate);           C
C    IATP=3, black carbon;  IATP=4, sea salt;                             C
C    For example:                                                         C
C    MSS(ILAY,1) = total dust concentration at layer ILAY                 C
C--- aerosol mass concentration input END --------------------------------C
        
C--- convert aerosol mass ug/m3 to #/cm3 ---------------------------------C
       DO IATP=1, NAETP                     !aerosol type loop       
        IF(MSS(ILAY,IATP).LT.(1.E-3)) MSS(ILAY,IATP)=0.0
        NDS (ILAY,IATP)=MSS (ILAY,IATP)*CONV(IATP)
        Naer(ILAY)=Naer(ILAY)+NDS(ILAY,IATP)
       END DO                               !aerosol type loop END
       
      END DO                               !vertical layer loop END
C-------------------------------------------------------------------------C

C--- optical properties calculation at each level ------------------------C
      DO ILAY= 1,nzstem                   !vertical layer loop

C--- get relative humidity value -----------------------------------------C
C    RH value rhu(COL,ROW,ILAY) is required                               C
C      rhu(ILAY)= 0.
      
      RHUM = rhu(ILAY)/100.
      
      IF(RHUM.GT.0.99) RHUM=0.99
      IF(RHUM.LT.0.00) RHUM=0.00
       
      DO IATP=1, NAETP               !aerosol type loop
         ABEAR(1,IATP)=0.
         AOMR (1,IATP)=0.
         AGAR (1,IATP)=0.

       DO IW =2, 7                   !wavelength loop
         ABEAR(IW,IATP)=0.
         AOMR (IW,IATP)=0.
         AGAR (IW,IATP)=0.

       IF(IATP.EQ.2.OR.IATP.EQ.4) THEN        !RH sensitive
        DO IRH=1,7                            !relative humidity loop
        IF((RHUM.GE.RH00(IRH)).AND.(RHUM.LE.RH00(IRH+1))) THEN
         ABEAR(IW,IATP)=SBEA(IW-1,IRH,IATP)+ ((RHUM-RH00(IRH))/
     >            (RH00(IRH+1)-RH00(IRH)))*
     >            (SBEA(IW-1,IRH+1,IATP)-SBEA(IW-1,IRH,IATP))
         AOMR (IW,IATP)=SBAA(IW-1,IRH,IATP)+ ((RHUM-RH00(IRH))/
     >            (RH00(IRH+1)-RH00(IRH)))*
     >            (SBAA(IW-1,IRH+1,IATP)-SBAA(IW-1,IRH,IATP))
         AGAR (IW,IATP)=SGA (IW-1,IRH,IATP)+ ((RHUM-RH00(IRH))/
     >            (RH00(IRH+1)-RH00(IRH)))*
     >            (SGA (IW-1,IRH+1,IATP)-SGA (IW-1,IRH,IATP))
         END IF
        END DO                                !relative humidity loop END
       ELSE                                   !RH independent
         ABEAR(IW,IATP)=SBEA(IW-1,1,IATP)
         AOMR (IW,IATP)=SBAA(IW-1,1,IATP)
         AGAR (IW,IATP)=SGA (IW-1,1,IATP)
       END IF                                 !RH ENDIF
        
       END DO                        !wavelength loop END

       ABEAR(1,IATP)= ABEAR(2,IATP)+(ABEAR(2,IATP)-ABEAR(3,IATP))*
     >          (wavei(1)-wavei(2))/(wavei(2)-wavei(3))
       AOMR (1,IATP)= AOMR (2,IATP)+(AOMR (2,IATP) -AOMR(3,IATP))*
     >          (wavei(1)-wavei(2))/(wavei(2)-wavei(3))
       AGAR (1,IATP)= AGAR (2,IATP)+(AGAR (2,IATP) -AGAR(3,IATP))*
     >          (wavei(1)-wavei(2))/(wavei(2)-wavei(3))
     
       DO IW=1,7                     !wavelength loop
        extaer(IW)=ABEAR(IW,IATP)
        albaer(IW)=AOMR (IW,IATP)
        asmaer(IW)=AGAR (IW,IATP)
       END DO                        !wavelength loop END

      n=7
      x1(1:n)=wavei(1:n)
      y1(1:n)=extaer(1:n)
      CALL addpnt(x1,y1,kdata,n,x1(1)*(1.-deltax),y1(1))
      CALL addpnt(x1,y1,kdata,n,          0.,y1(1))
      CALL addpnt(x1,y1,kdata,n,x1(n)*(1.+deltax),   0.)
      CALL addpnt(x1,y1,kdata,n,      1.e+38,   0.)     
      CALL inter2(NW, wl, BEA, n, x1, y1, ierr)

      n=7
      x1(1:n)=wavei(1:n)
      y1(1:n)=albaer(1:n)
      CALL addpnt(x1,y1,kdata,n,x1(1)*(1.-deltax),y1(1))
      CALL addpnt(x1,y1,kdata,n,          0.,y1(1))
      CALL addpnt(x1,y1,kdata,n,x1(n)*(1.+deltax),   0.)
      CALL addpnt(x1,y1,kdata,n,      1.e+38,   0.)           
      CALL inter2(NW, wl, OMA, n, x1, y1, ierr)
      
      n=7
      x1(1:n)=wavei(1:n)
      y1(1:n)=asmaer(1:n)
      CALL addpnt(x1,y1,kdata,n,x1(1)*(1.-deltax),y1(1))
      CALL addpnt(x1,y1,kdata,n,          0.,y1(1))
      CALL addpnt(x1,y1,kdata,n,x1(n)*(1.+deltax),   0.)
      CALL addpnt(x1,y1,kdata,n,      1.e+38,   0.)
      CALL inter2(NW, wl, GA,  n, x1, y1, ierr)
       
      DO IW=1,NW-1                    !wavelength loop
       
       IF(BEA(IW).LT.0.OR.OMA(IW).LT.0.OR.GA(IW).LT.0) THEN
        print *, "Optical < 0 for aerosol ", IATP, " at layer ",ILAY
        STOP
       ENDIF

       if(ga(iw).gt.1) then
        print*,'ga(iw)>1 at iw,ilay=',iw,ilay, ga(iw)
	stop
       endif	

C--- calculating aerosol optical depth -----------------------------------C
C    Z(I) = altitude in km for each level, BEA = extin in 1/km per 1/cm3  C
C    watch out for the unit                                               C

       IF(ILAY.LE.nzstem-1) then
        dtaerwl(ILAY,IW,IATP)=BEA(IW)*
     1	((NDS(ILAY,IATP)+NDS(ILAY+1,IATP))/2.)*(Z(ILAY+1)-Z(ILAY))
        if(dtaerwl(ILAY,IW,IATP).lt.0) then
	 print*,'negative dtaerwl in setaer ',dtaerwl(ILAY,IW,IATP),
     1	  bea(IW),NDS(ILAY,IATP),NDS(ILAY+1,IATP),Z(ILAY+1),Z(ILAY),
     2    iLAY, IATP, IW
         print*,'nz,nztem,z=',nz,nzstem,z(1:nzstem)
         stop
	endif
       endif	 
       
       omaerwl (ILAY,IW,IATP) = OMA(IW)
       gaerwl  (ILAY,IW,IATP) = GA (IW)
       
       END DO                        !wavelength loop END
      END DO                         !aerosol type loop END
      
      IF(Naer(ILAY).GE.1.E-6) THEN   !IF Naer, avoid zero
C--- aerosol external mixture --------------------------------------------C
       DO IW=1,NW                    !wavelength loop
        dtaerwl(ILAY,IW,NAETP+1)=0.
	omaerwl(ILAY,IW,NAETP+1)=0.
	gaerwl (ILAY,IW,NAETP+1)=0.
	absorp (ILAY,IW)=0.
	
        DO IATP=1, NAETP             !aerosol type loop
         dtaerwl(ILAY,IW,NAETP+1)= dtaerwl(ILAY,IW,NAETP+1)
     >                         + dtaerwl(ILAY,IW,IATP)                 ! total aerosol optical depth  
         omaerwl(ILAY,IW,NAETP+1)= omaerwl(ILAY,IW,NAETP+1)
     >                         + omaerwl(ILAY,IW,IATP)*NDS(ILAY,IATP)
         absorp (ILAY,IW)= absorp(ILAY,IW)+dtaerwl(ILAY,IW,IATP)        ! absorption
     >                         * (1-omaerwl(ILAY,IW,IATP))
         gaerwl (ILAY,IW,NAETP+1)= gaerwl (ILAY,IW,NAETP+1)
     >                         + gaerwl(ILAY,IW,IATP)*NDS(ILAY,IATP)
        END DO                       !aerosol type loop END
             
         omaerwl(ILAY,IW,NAETP+1)= omaerwl (ILAY,IW,NAETP+1)/Naer(ILAY)   ! average back-scattering albedo
         gaerwl (ILAY,IW,NAETP+1)= gaerwl (ILAY,IW,NAETP+1)/Naer(ILAY)    ! average asymmetry factor
	 
         dtaer(ilay,iw)= amax1(dtaerwl(ILAY,IW,NAETP+1),dtaer(ilay,iw))   ! avoid dtaer from zero 
	 if(dtaer(ilay,iw).lt.1e-25) then
	  omaer(ilay,iw)= 1.
	 else 
 	  omaer(ilay,iw)= 1-absorp(ILAY,IW)/dtaer(ilay,iw)
	 endif 
cccc	   omaer(ilay,iw)= omaerwl(ILAY,IW,NAETP+1)
	 gaer(ilay,iw)=gaerwl(ILAY,IW,NAETP+1) 
         if(gaer(ilay,iw).gt.1.or.gaer(ilay,iw).lt.0.) then
          print*,'gaer wrong at setaer, iw,ilay=',iw,ilay,gaer(ilay,iw)
	  print*,'gaerwl(ilay,iw,1:5)=',gaerwl(ilay,iw,1:5)
	  print*,'nds(ilay,1:5),naer(ilay)=',nds(ilay,1:5),naer(ilay)
	  stop
         endif	
	 
       END DO                        !wavelength loop END

      END IF                         !Naer ENDIF
      END DO                         !vertical layer loop END 
  

C---- DAOD is at 550nm wavelength (IW=102)?-------------------------------C
      DO iw=1,nw
       if(wl(iw).ge.550) goto 23
      enddo
 23   if(abs(wl(iw-1)-550).lt.abs(wl(iw)-550)) iw=iw-1
   
      DO ILAY=1, NZSTEM-1
       DO IATP=1, NAETP+1      
        DAOD(ILAY,IATP)=dtaerwl(ilay,iw,iatp)/(Z(ILAY+1)-Z(ILAY))  ! extincting coefficient for each species
       END DO                                                      ! and total
c       DAOD(ILAY,NAETP+2)= omaer(ILAY,IW)   ! the 7th variable is Single Scaterring Albedo
      END DO
      
c      DO ILAY=NZSTEM-1, 1, -1
c      DAOD(ILAY,IATP+1)=DAOD(ILAY+1,IATP+1)           ! aerosol optical depth from top layer
c     >                       +dtaerwl(ILAY,IW,IATP+1)
c      END DO

*_______________________________________________________________________

      RETURN
      END



      SUBROUTINE setair(pmbnew,
     $     nz,z,nw,wl,
     $     airlev,dtrl,cz,airstem,nzstem)

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Set up an altitude profile of air molecules.  Subroutine includes a      =*
*=  shape-conserving scaling method that allows scaling of the entire        =*
*=  profile to a given sea-level pressure.                                   =*
*-----------------------------------------------------------------------------*
*=  PARAMETERS:                                                              =*
*=  PMBNEW  - REAL, sea-level pressure (mb) to which profile should be    (I)=*
*=            scaled.  If PMBNEW < 0, no scaling is done                     =*
*=  NZ      - INTEGER, number of specified altitude levels in the working (I)=*
*=            grid                                                           =*
*=  Z       - REAL, specified altitude working grid (km)                  (I)=*
*=  NW      - INTEGER, number of specified intervals + 1 in working       (I)=*
*=            wavelength grid                                                =*
*=  WL      - REAL, vector of lower limits of wavelength intervals in     (I)=*
*=            working wavelength grid                                        =*
*=  AIRLEV  - REAL, air density (molec/cc) at each specified altitude     (O)=* 
*=  DTRL    - REAL, Rayleigh optical depth at each specified altitude     (O)=*
*=            and each specified wavelength                                  =*
*=  CZ      - REAL, number of air molecules per cm^2 at each specified    (O)=*
*=            altitude layer                                                 =*
*-----------------------------------------------------------------------------*
*=  EDIT HISTORY:                                                            =*
*=  02/97  Read in profile from an input file                                =*
*-----------------------------------------------------------------------------*
*= This program is free software;  you can redistribute it and/or modify     =*
*= it under the terms of the GNU General Public License as published by the  =*
*= Free Software Foundation;  either version 2 of the license, or (at your   =*
*= option) any later version.                                                =*
*= The TUV package is distributed in the hope that it will be useful, but    =*
*= WITHOUT ANY WARRANTY;  without even the implied warranty of MERCHANTIBI-  =*
*= LITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public     =*
*= License for more details.                                                 =*
*= To obtain a copy of the GNU General Public License, write to:             =*
*= Free Software Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.   =*
*-----------------------------------------------------------------------------*
*= To contact the authors, please mail to:                                   =*
*= Sasha Madronich, NCAR/ACD, P.O.Box 3000, Boulder, CO, 80307-3000, USA  or =*
*= send email to:  sasha@ucar.edu                                            =*
*-----------------------------------------------------------------------------*
*= Copyright (C) 1994,95,96  University Corporation for Atmospheric Research =*
*-----------------------------------------------------------------------------*

      IMPLICIT NONE
      INCLUDE 'tuv.params'

      INTEGER kdata
      PARAMETER(kdata=150)

* input: (grids)
      REAL wl(kw)
      REAL z(kz)
      INTEGER nw
      INTEGER nz,nzstem
      REAL  pmbnew,airstem(nzstem)

* output:
* air density (molec cm-3) at each grid level
* Rayleigh optical depths

      REAL airlev(kz)
      REAL dtrl(kz,kw)

* local:
      REAL scale
      real airnew(kdata)
      REAL colold, colnew, pmbold
      REAL pconv
      PARAMETER(pconv = 980.665 * 1.E-3 * 28.9644 / 6.022169E23)
* specified data:
      REAL zd(kdata), air(kdata)
      REAL hscale
      REAL cd(kdata)

* other:
      REAL cz(kz)
      REAL srayl(kw)
      REAL deltaz
      REAL colz, pressz
      REAL wc, wmicrn, xx 
      INTEGER i, iw, nd
      parameter(nd=121)
      
       data zd(1:nd)/0.,  1.,  2.,  3.,  4.,  5.,  6.,  7.,  8.,
     1   9., 10., 11., 12., 13., 14., 15., 16., 17.,
     2  18., 19., 20., 21., 22., 23., 24., 25., 26.,
     3  27., 28., 29., 30., 31., 32., 33., 34., 35.,
     4  36., 37., 38., 39., 40., 41., 42., 43., 44.,
     5  45., 46., 47., 48., 49., 50., 51., 52., 53.,
     6  54., 55., 56., 57., 58., 59., 60., 61., 62.,
     7  63., 64., 65., 66., 67., 68., 69., 70., 71.,
     8  72., 73., 74., 75., 76., 77., 78., 79., 80.,
     9  81., 82., 83., 84., 85., 86., 87., 88., 89.,
     a  90., 91., 92., 93., 94., 95., 96., 97., 98.,
     b  99.,100.,101.,102.,103.,104.,105.,106.,107.,
     c  108.,109.,110.,111.,112.,113.,114.,115.,116.,
     d  117.,118.,119.,120./
       data air(1:nd)/.255E+20,.231E+20,.209E+20,.189E+20,
     1 .170E+20,.153E+20,.137E+20,.123E+20,.109E+20,.971E+19,.860E+19,
     2 .759E+19,.649E+19,.554E+19,.474E+19,.405E+19,.346E+19,.296E+19,
     3 .253E+19,.216E+19,.185E+19,.157E+19,.134E+19,.114E+19,.976E+18,
     4 .833E+18,.712E+18,.609E+18,.521E+18,.447E+18,.383E+18,.328E+18,
     5 .282E+18,.241E+18,.206E+18,.176E+18,.151E+18,.130E+18,.112E+18,
     6 .962E+17,.831E+17,.719E+17,.623E+17,.540E+17,.470E+17,.409E+17,
     7 .356E+17,.311E+17,.274E+17,.242E+17,.214E+17,.189E+17,.168E+17,
     8 .149E+17,.133E+17,.118E+17,.105E+17,.930E+16,.824E+16,.729E+16,
     9 .644E+16,.568E+16,.500E+16,.440E+16,.387E+16,.339E+16,.297E+16,
     a .260E+16,.227E+16,.198E+16,.172E+16,.150E+16,.130E+16,.112E+16,
     b .964E+15,.830E+15,.713E+15,.612E+15,.525E+15,.449E+15,.384E+15,
     c .327E+15,.279E+15,.237E+15,.202E+15,.171E+15,.144E+15,.121E+15,
     d .101E+15,.850E+14,.712E+14,.596E+14,.499E+14,.418E+14,.349E+14,
     e .292E+14,.244E+14,.204E+14,.170E+14,.142E+14,.119E+14,.999E+13,
     f .840E+13,.707E+13,.596E+13,.502E+13,.424E+13,.358E+13,.302E+13,
     g .255E+13,.214E+13,.180E+13,.152E+13,.130E+13,.111E+13,.968E+12,
     h .843E+12,.738E+12,.650E+12,.575E+12,.510E+12/

* External functions:
      REAL fsum
      EXTERNAL fsum

*_______________________________________________________________________

* read in air density profile

* compute column increments (logarithmic integrals)

      DO 6, i = 1, nd - 1
         deltaz = 1.E5 * (zd(i+1)-zd(i)) 
         cd(i) =  (air(i+1)-air(i)) /ALOG(air(i+1)/air(i)) * deltaz
C         cd(i) = (air(i+1)+air(i)) * deltaz / 2. 
    6 CONTINUE

* Include exponential tail integral from infinity to 50 km,
* fold tail integral into top layer
* specify scale height near top of data.

      hscale = 8.05e5
      cd(nd-1) = cd(nd-1) + hscale * air(nd)

* alternative input air density data could include, e.g., a read file here:

* If want, can rescale to any total pressure:

      colold = fsum(nd-1,cd)
      pmbold = colold * pconv 
c      WRITE(kout,100) colold, pmbold
  100 FORMAT(5x,'old sea level air column = ', 1pe11.4,1x,'# cm-2  = ',
     $     0pf8.2,' mbar')

* assign new sea level pressure

      if (pmbnew .lt. 0.) then
         scale = 1.
      else
         scale = pmbnew/pmbold
      endif

      DO i = 1, nd-1
         cd(i) = cd(i) * scale
         airnew(i) = air(i) * scale
      ENDDO
      airnew(nd) = air(nd) * scale
      
      colnew = fsum(nd-1,cd)
c      WRITE(kout,105) colnew, colnew * pconv
  105 FORMAT(5x,'new sea level air column = ', 1pe11.4,1x,'# cm-2  = ',
     $     0pf8.2,' mbar')

********************** end data input.

* Compute air density at each level

      CALL inter1(nz,z,airlev,nd,zd,airnew)
      airlev(1:nzstem)=airstem(1:nzstem)   ! air density from STEM 
      
* Compute column increments on standard z-grid.  

      CALL inter3(nz,z,cz, nd,zd,cd, 1) ! perform weight interpolation, the overlapped
                                        ! length is taken as weight. 
      DO i = 1, nzstem- 1
        deltaz = 1.E5 * (z(i+1)-z(i))
	if(abs(airstem(i+1)-airstem(i)).lt.1.and.i.le.(nzstem-2)) 
     1	  then
	 airstem(i+1)=airstem(i)+(airstem(i+2)-airstem(i))/(
     1	  z(i+2)-z(i))*(z(i+1)-z(i))        ! avoid zero difference of airstem
        endif
        cz(i) =  (airstem(i+1)-airstem(i))/         ! loading data from STEM
     1	   ALOG(airstem(i+1)/airstem(i)) * deltaz
         
	if(cz(i).lt.0.or.(.not.abs(cz(i)).lt.1e25)) then
	 print*,'cz wrong in setair ',i,cz(i),deltaz
	 print*,'z=',z(1:nzstem)
	 print*,'airstem=',airstem(1:nzstem)
	 stop
	endif 
      enddo  
      
      colz = fsum(nz-1,cz)
      pressz =  colz * pconv
c      write(kout,110) colz, pressz
 110  FORMAT(5x,'surface air column = ', 1pe11.4,1x,'# cm-2  = ',
     $     0pf8.2,' mbar')


* compute Rayleigh cross sections and depths:

      DO 30, iw = 1, nw - 1
         wc = (wl(iw) + wl(iw+1))/2.

* Rayleigh scattering cross section from WMO 1985 (originally from
* Nicolet, M., On the molecular scattering in the terrestrial atmosphere:
* An empirical formula for its calculation in the homoshpere, Planet.
* Space Sci., 32, 1467-1468, 1984.

         wmicrn =  wc/1.E3
         IF( wmicrn .LE. 0.55) THEN
            xx = 3.6772 + 0.389*wmicrn + 0.09426/wmicrn
         ELSE
            xx = 4. + 0.04
         ENDIF
         srayl(iw) = 4.02e-28/(wmicrn)**xx

* alternate (older) expression from
* Frohlich and Shaw, Appl.Opt. v.11, p.1773 (1980).
C     xx = 3.916 + 0.074*wmicrn + 0.050/wmicrn
C     srayl(iw) = 3.90e-28/(wmicrn)**xx

         DO 40, i = 1, nz - 1
           dtrl(i,iw) = cz(i)*srayl(iw)
c	   if(.not.(dtrl(i,iw).gt.-1.e36)) dtrl(i,iw)=0.   ! avoid overflow
	   
	   if(dtrl(i,iw).lt.0.or.(.not.abs(dtrl(i,iw)).lt.1e20)) then
	    print*,'wrong dtrl ',dtrl(i,iw),cz(i),srayl(iw)
	    print*,'airlev=',airlev(1:nz) 
	    print*,'z=',z(1:nzstem)
	    print*,'airstem=',airstem(1:nzstem)
	    stop
	   endif  
   40    CONTINUE

   30 CONTINUE
*_______________________________________________________________________

      RETURN
      END
      SUBROUTINE setalb(nw,wl,albedo)

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Set the albedo of the surface.  The albedo is assumed to be Lambertian,  =*
*=  i.e., the reflected light is isotropic, and idependt of the direction    =*
*=  of incidence of light.  Albedo can be chosen to be wavelength dependent. =*
*-----------------------------------------------------------------------------*
*=  PARAMETERS:                                                              =*
*=  NW      - INTEGER, number of specified intervals + 1 in working       (I)=*
*=            wavelength grid                                                =*
*=  WL     - REAL, vector of lower limits of wavelength intervals in      (I)=*
*=           working wavelength grid                                         =*
*=  ALBEDO  - REAL, surface albedo at each specified wavelength           (O)=*
*-----------------------------------------------------------------------------*
*=  EDIT HISTORY:                                                            =*
*=  Original                                                                 =*
*-----------------------------------------------------------------------------*
*= This program is free software;  you can redistribute it and/or modify     =*
*= it under the terms of the GNU General Public License as published by the  =*
*= Free Software Foundation;  either version 2 of the license, or (at your   =*
*= option) any later version.                                                =*
*= The TUV package is distributed in the hope that it will be useful, but    =*
*= WITHOUT ANY WARRANTY;  without even the implied warranty of MERCHANTIBI-  =*
*= LITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public     =*
*= License for more details.                                                 =*
*= To obtain a copy of the GNU General Public License, write to:             =*
*= Free Software Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.   =*
*-----------------------------------------------------------------------------*
*= To contact the authors, please mail to:                                   =*
*= Sasha Madronich, NCAR/ACD, P.O.Box 3000, Boulder, CO, 80307-3000, USA  or =*
*= send email to:  sasha@ucar.edu                                            =*
*-----------------------------------------------------------------------------*
*= Copyright (C) 1994,95,96  University Corporation for Atmospheric Research =*
*-----------------------------------------------------------------------------*

      IMPLICIT NONE
      INCLUDE 'tuv.params'

* input: (wavelength working grid data)
      INTEGER nw
      REAL wl(kw)

* output:
      REAL albedo(kw)

* local:
      INTEGER iw
      REAL alb
*_______________________________________________________________________

* set

      alb = 0.10
c      WRITE(kout,*)'wavelength-independent albedo = ', alb
      DO 10, iw = 1, nw - 1

         albedo(iw) = alb

   10 CONTINUE
*_______________________________________________________________________

      RETURN
      END
      SUBROUTINE setcld(nz,z,nw,wl,dtcld,omcld,gcld)

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Set cloud properties for each specified altitude layer.  Properties      =*
*=  may be wavelength dependent.                                             =*
*-----------------------------------------------------------------------------*
*=  PARAMETERS:                                                              =*
*=  NZ      - INTEGER, number of specified altitude levels in the working (I)=*
*=            grid                                                           =*
*=  Z       - REAL, specified altitude working grid (km)                  (I)=*
*=  NW      - INTEGER, number of specified intervals + 1 in working       (I)=*
*=            wavelength grid                                                =*
*=  WL     - REAL, vector of lower limits of wavelength intervals in      (I)=*
*=           working wavelength grid                                         =*
*=  DTCLD   - REAL, optical depth due to absorption by clouds at each     (O)=*
*=            altitude and wavelength                                        =*
*=  OMCLD   - REAL, single scattering albedo due to clouds at each        (O)=*
*=            defined altitude and wavelength                                =*
*=  GCLD    - REAL, cloud asymmetry factor at each defined altitude and   (O)=*
*=            wavelength                                                     =*
*-----------------------------------------------------------------------------*
*=  EDIT HISTORY:                                                            =*
*=  12/94  Bug fix                                                           =*
*-----------------------------------------------------------------------------*
*= This program is free software;  you can redistribute it and/or modify     =*
*= it under the terms of the GNU General Public License as published by the  =*
*= Free Software Foundation;  either version 2 of the license, or (at your   =*
*= option) any later version.                                                =*
*= The TUV package is distributed in the hope that it will be useful, but    =*
*= WITHOUT ANY WARRANTY;  without even the implied warranty of MERCHANTIBI-  =*
*= LITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public     =*
*= License for more details.                                                 =*
*= To obtain a copy of the GNU General Public License, write to:             =*
*= Free Software Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.   =*
*-----------------------------------------------------------------------------*
*= To contact the authors, please mail to:                                   =*
*= Sasha Madronich, NCAR/ACD, P.O.Box 3000, Boulder, CO, 80307-3000, USA  or =*
*= send email to:  sasha@ucar.edu                                            =*
*-----------------------------------------------------------------------------*
*= Copyright (C) 1994,95,96  University Corporation for Atmospheric Research =*
*-----------------------------------------------------------------------------*

      IMPLICIT NONE
      INCLUDE 'tuv.params'

      INTEGER kdata
      PARAMETER(kdata=51)

* input: (grids)
      REAL wl(kw)
      REAL z(kz)
      INTEGER nz
      INTEGER nw

* Output: 
      REAL dtcld(kz,kw), omcld(kz,kw), gcld(kz,kw)

* local:

* specified data:
      REAL zd(kdata), cd(kdata), omd(kdata), gd(kdata)
      REAL womd(kdata), wgd(kdata)

* other:
      REAL cz(kz)
      REAL omz(kz)
      REAL gz(kz)
      INTEGER i, iw, n

* External functions:
      REAL fsum
      EXTERNAL fsum
*_______________________________________________________________________


* cloud properties are set for each layer (not each level)

* Set as many clouds as want here:
* First choose a cloud grid, zd(n), in km above sea level
* Can allow altitude variation of omega, g:

      n = 4
      
      zd(1) = 5.
      cd(1) = 0.          ! cloud optical depth
      omd(1) = .9999
      gd(1) = .85

      zd(2) = 7.
      cd(2) = 0.
      omd(2) = .5
      gd(2) = .5

      zd(3) = 9.
      cd(3) = 0.
      omd(3) = .9999
      gd(3) = .85

      zd(4) = 11.

******************

* compute integrals and averages over grid layers:
* for g and omega, use averages weigthed by optical depth

C     DO 11, i = 1, n    !***** CHANGED!!See header!!*****
      DO 11, i = 1, n-1
         womd(i) = omd(i) * cd(i)
         wgd(i) = gd(i) * cd(i)
   11 CONTINUE
      CALL inter3(nz,z,cz,  n, zd,cd, 0)
      CALL inter3(nz,z,omz, n, zd,womd, 0)
      CALL inter3(nz,z,gz , n, zd,wgd, 0)

      DO 15, i = 1, nz-1
         IF (cz(i) .GT. 0.) THEN
            omz(i) = omz(i)/cz(i)
            gz(i)  = gz(i) /cz(i)
         ELSE
            omz(i) = 1.
            gz(i) = 0.
         ENDIF
   15 CONTINUE

c      WRITE(kout,*) 'Cloud: ', n, 'levels, tot opt. dep. = ', 
c     $     fsum(nz-1,cz)

* assign at all wavelengths
* (can move wavelength loop outside if want to vary with wavelength)

      DO 17, iw = 1, nw-1
         DO 16, i = 1, nz-1
            dtcld(i,iw) = cz(i)
            omcld(i,iw) = omz(i)
            gcld (i,iw) = gz(i)
   16    CONTINUE
   17 CONTINUE
*_______________________________________________________________________

      RETURN
      END
      SUBROUTINE setno2(no2new,
     $     nz,z,nw,wl,
     $     xsno2, tlay,
     $     dtno2,no2stem,nzstem)

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Set up an altitude profile of NO2 molecules, and corresponding absorption=*
*=  optical depths.  Subroutine includes a shape-conserving scaling method   =*
*=  that allows scaling of the entire profile to a given overhead NO2        =*
*=  column amount.                                                           =*
*-----------------------------------------------------------------------------*
*=  PARAMETERS:                                                              =*
*=  NO2NEW - REAL, overhead NO2 column amount (molec/cm^2) to which       (I)=*
*=           profile should be scaled.  If NO2NEW < 0, no scaling is done    =*
*=  NZ     - INTEGER, number of specified altitude levels in the working  (I)=*
*=           grid                                                            =*
*=  Z      - REAL, specified altitude working grid (km)                   (I)=*
*=  NW     - INTEGER, number of specified intervals + 1 in working        (I)=*
*=           wavelength grid                                                 =*
*=  WL     - REAL, vector of lower limits of wavelength intervals in      (I)=*
*=           working wavelength grid                                         =*
*=  XSNO2  - REAL, molecular absoprtion cross section (cm^2) of O2 at     (I)=*
*=           each specified wavelength                                       =*
*=  TLAY   - REAL, temperature (K) at each specified altitude layer       (I)=*
*=  DTNO2  - REAL, optical depth due to NO2 absorption at each            (O)=*
*=           specified altitude at each specified wavelength                 =*
*-----------------------------------------------------------------------------*
*=  EDIT HISTORY:                                                            =*
*=  03/97  fix DO-10 loop                                                    =*
*-----------------------------------------------------------------------------*
*= This program is free software;  you can redistribute it and/or modify     =*
*= it under the terms of the GNU General Public License as published by the  =*
*= Free Software Foundation;  either version 2 of the license, or (at your   =*
*= option) any later version.                                                =*
*= The TUV package is distributed in the hope that it will be useful, but    =*
*= WITHOUT ANY WARRANTY;  without even the implied warranty of MERCHANTIBI-  =*
*= LITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public     =*
*= License for more details.                                                 =*
*= To obtain a copy of the GNU General Public License, write to:             =*
*= Free Software Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.   =*
*-----------------------------------------------------------------------------*
*= To contact the authors, please mail to:                                   =*
*= Sasha Madronich, NCAR/ACD, P.O.Box 3000, Boulder, CO, 80307-3000, USA  or =*
*= send email to:  sasha@ucar.edu                                            =*
*-----------------------------------------------------------------------------*
*= Copyright (C) 1994,95,96  University Corporation for Atmospheric Research =*
*-----------------------------------------------------------------------------*

      IMPLICIT NONE
      INCLUDE 'tuv.params'

      INTEGER kdata
      PARAMETER(kdata=51)

********
* input:
********

* grids:

      REAL wl(kw)
      REAL z(kz)
      INTEGER nw, nz, nzstem
      REAL no2new, no2stem(nzstem)

* mid-layer temperature:

      REAL tlay(kz) 

********
* output:
********

      REAL dtno2(kz,kw)

********
* local:
********

* absorption cross sections 

      REAL xsno2(kw)
      REAL cz(kz)

* nitrogen dioxide profile data:

      REAL zd(kdata), no2(kdata)
      REAL cd(kdata)
      REAL hscale
      REAL colold, colnew
      REAL scale
      REAL sno2

* other:

      INTEGER i, l, nd

********
* External functions:
********
      REAL fsum
      EXTERNAL fsum

*_______________________________________________________________________
* Data input:

* Example:  set to 1 ppb in lowest 1 km, set to zero above that.
* - do by specifying concentration at 3 altitudes.

c      write(kout,*) 'NO2:  1 ppb in lowest 1 km, 0 above'

      nd = 3
      zd(1) = 0.
      no2(1) = 1. * 2.69e10

      zd(2) = 1.
      no2(2) = 1. * 2.69e10

      zd(3) = zd(2)* 1.000001
      no2(3) = 0.

C     zd(4) = zd(3)*1.1 
C     no2(4) = 0.

* compute column increments (alternatively, can specify these directly)

      DO 11, i = 1, nd - 1
         cd(i) = (no2(i+1)+no2(i)) * 1.E5 * (zd(i+1)-zd(i)) / 2. 
   11 CONTINUE

* Include exponential tail integral from top level to infinity.
* fold tail integral into top layer
* specify scale height near top of data (use ozone value)

      hscale = 4.50e5
      cd(nd-1) = cd(nd-1) + hscale * no2(nd)

***********
*********** end data input.

* Compute column increments on standard z-grid.  

      CALL inter3(nz,z,cz, nd,zd,cd, 1)
* scale values of cz(i) 

      colold = fsum(nz-1,cz)
c      WRITE(kout,100) colold, colold/2.687E16
  100 FORMAT(5x,'old NO2 Column = ', 1pe11.4,1x,'# cm-2  = ',
     $     0pf8.2, '  Dobson Units ')

      if ( (no2new .lt. 0.)  .or.  (colold .le. 0.) ) then
         scale = 1.
      else
         scale =  2.687e16*no2new/colold
      endif

      do i = 1, nz-1
         cz(i) = cz(i) * scale
      enddo
      DO i = 1, nzstem - 1
         cz(i) = (no2stem(i+1)+no2stem(i)) * 1.E5 *  ! load NO2 from STEM
     1	 (z(i+1)-z(i)) / 2. 
      enddo

      colnew = fsum(nz-1,cz)
c      WRITE(kout,105) colnew, colnew/2.687E16
  105 format(5x,'new NO2 Column = ', 1pe11.4,1x,'# cm-2  = ',
     $     0pf8.2, '  Dobson Units ')

************************************
* calculate optical depth for each layer, with temperature 
* correction.  Output, dtno2(kz,kw)

      DO 20, l = 1, nw-1
         sno2 = xsno2(l)
         DO 10, i = 1, nz-1
            dtno2(i,l) = cz(i)*sno2
   10    CONTINUE
   20 CONTINUE
*_______________________________________________________________________

      RETURN
      END

      SUBROUTINE seto2(nz,z,nw,wl,cz,vcol,scol,dto2,xso2)

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Compute equivalent optical depths for O2 absorption, parameterized in    =*
*=  the SR bands and the Lyman-alpha line.                                   =*
*-----------------------------------------------------------------------------*
*=  PARAMETERS:                                                              =*
*=  NZ      - INTEGER, number of specified altitude levels in the working (I)=*
*=            grid                                                           =*
*=  Z       - REAL, specified altitude working grid (km)                  (I)=*
*=  NW      - INTEGER, number of specified intervals + 1 in working       (I)=*
*=            wavelength grid                                                =*
*=  WL      - REAL, vector of lower limits of wavelength intervals in     (I)=*
*=            working wavelength grid                                        =*
*=  CZ      - REAL, number of air molecules per cm^2 at each specified    (I)=*
*=            altitude layer                                                 =*
*=  ZEN     - REAL, solar zenith angle                                    (I)=*
*=  DTO2    - REAL, optical depth due to O2 absorption at each specified  (O)=*
*=            vertical layer at each specified wavelength                    =*
*=  XSO2    - REAL, molecular absorption cross section in SR bands at     (O)=*
*=            each specified altitude and wavelength.  Includes Herzberg     =*
*=            continuum.                                                     =*
*-----------------------------------------------------------------------------*
*=  EDIT HISTORY:                                                            =*
*=  02/98  Included Lyman-alpha parameterization                             =*
*=  03/97  Fix dto2 problem at top level (nz)                                =*
*=  02/97  Changed offset for grid-end interpolation to relative number      =*
*=         (x * (1 +- deltax))                                               =*
*=  08/96  Modified for early exit, no redundant read of data and smaller    =*
*=         internal grid if possible;  internal grid uses user grid points   =*
*=         whenever possible                                                 =*
*=  07/96  Modified to work on internal grid and interpolate final values    =*
*=         onto the user-defined grid                                        =*
*-----------------------------------------------------------------------------*
*= This program is free software;  you can redistribute it and/or modify     =*
*= it under the terms of the GNU General Public License as published by the  =*
*= Free Software Foundation;  either version 2 of the license, or (at your   =*
*= option) any later version.                                                =*
*= The TUV package is distributed in the hope that it will be useful, but    =*
*= WITHOUT ANY WARRANTY;  without even the implied warranty of MERCHANTIBI-  =*
*= LITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public     =*
*= License for more details.                                                 =*
*= To obtain a copy of the GNU General Public License, write to:             =*
*= Free Software Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.   =*
*-----------------------------------------------------------------------------*
*= To contact the authors, please mail to:                                   =*
*= Sasha Madronich, NCAR/ACD, P.O.Box 3000, Boulder, CO, 80307-3000, USA  or =*
*= send email to:  sasha@ucar.edu                                            =*
*-----------------------------------------------------------------------------*
*= Copyright (C) 1994,95,96  University Corporation for Atmospheric Research =*
*-----------------------------------------------------------------------------*

      IMPLICIT NONE
      INCLUDE 'tuv.params'

      REAL wl(kw)
      REAL z(kz), cz(kz)
      INTEGER nz, nw

      REAL vcol(kz), scol(kz)
      REAL dto2(kz,kw), xso2(kz,kw)
      REAL secchi(kz)

* grid on which Kockarts' parameterization is defined
      INTEGER ngast
      PARAMETER (ngast = 17)
      REAL wlgast(ngast)
      SAVE wlgast
 
* O2 optical depth and equivalent cross section on Kockarts' grid
      REAL dto2k(kz,ngast-1), xso2k(kz,ngast-1)

* Lyman-alpha variables
      INTEGER nla
      PARAMETER (nla = 2)
      REAL wlla(nla)

* O2 optical depth and equivalent cross section in the Lyman-alpha region
      REAL dto2la(kz,nla-1), xso2la(kz,nla-1)

* internal grid and O2 cross section on internal grid
      INTEGER kdata
      PARAMETER (kdata = 200)
      REAL wlint(kdata), xso2int(kdata)
      SAVE xso2int
      INTEGER nwint
      SAVE nwint

* temporary one-dimensional storage for optical depth and cross section values
* XXtmp  - on internal grid
* XXuser - on user defined grid
      REAL dttmp(2*kw), xstmp(2*kw)
      REAL dtuser(kw), xsuser(kw)

      REAL o2col(kz)

* cross section data for use outside the SR-Bands (combined from
* Brasseur and Solomon and the JPL 1994 recommendation)
      INTEGER nosr,io2
      PARAMETER (nosr = 105,io2=101)
      REAL x1(nosr), y1(nosr)
       data x1(1:io2)/116.65,117.30,117.95,118.65,119.40,
     A 120.15,120.85,121.59,121.60,122.35,123.10,123.85,124.60,125.40,
     A 126.20,127.00,127.80,128.60,129.45,130.30,131.15,132.00,132.85,
     A 133.75,134.65,135.55,136.50,137.45,138.40,139.85,141.80,143.85,
     A 145.95,148.10,150.35,152.65,155.00,157.45,160.00,162.60,165.30,
     A 168.10,170.95,173.15,174.65,176.20,177.80,179.40,181.00,182.65,
     A 184.35,186.05,187.80,189.60,191.40,193.25,195.15,197.05,199.00,
     A 201.00,203.05,204.00,205.00,206.00,207.00,208.00,209.00,210.00,
     A 211.00,212.00,213.00,214.00,215.00,216.00,217.00,218.00,219.00,
     A 220.00,221.00,222.00,223.00,224.00,225.00,226.00,227.00,228.00,
     A 229.00,230.00,231.00,232.00,233.00,234.00,235.00,236.00,237.00,
     A 238.00,239.00,240.00,241.00,242.00,243.00/
       data y1(1:io2)/.2000E-19,.1250E-17,.2550E-18,.3000E-19,
     A .3750E-18,.4450E-17,.8350E-17,.1000E-19,.6000E-18,.2350E-18,
     A .4500E-18,.3350E-18,.1750E-16,.8950E-18,.4300E-18,.1100E-18,
     A .2050E-18,.4430E-18,.5550E-18,.4200E-18,.6850E-18,.1450E-17,
     A .2250E-17,.2300E-17,.4550E-17,.7230E-17,.9500E-17,.1230E-16,
     A .1320E-16,.1360E-16,.1400E-16,.1480E-16,.1410E-16,.1290E-16,
     A .1150E-16,.9910E-17,.8240E-17,.6580E-17,.4970E-17,.3450E-17,
     A .2080E-17,.1230E-17,.7220E-18,.4580E-18,.2740E-18,.7200E-20,
     A .3860E-20,.2160E-20,.1200E-20,.6600E-21,.3540E-21,.2000E-21,
     A .1200E-21,.5850E-22,.4600E-22,.3100E-22,.2500E-22,.1900E-22,
     A .1610E-22,.1320E-22,.1030E-22,.7500E-23,.7350E-23,.7130E-23,
     A .7050E-23,.6860E-23,.6680E-23,.6510E-23,.6240E-23,.6050E-23,
     A .5890E-23,.5720E-23,.5590E-23,.5350E-23,.5130E-23,.4880E-23,
     A .4640E-23,.4460E-23,.4260E-23,.4090E-23,.3890E-23,.3670E-23,
     A .3450E-23,.3210E-23,.2980E-23,.2770E-23,.2980E-23,.2770E-23,
     A .2630E-23,.2430E-23,.2250E-23,.2100E-23,.1940E-23,.1780E-23,
     A .1630E-23,.1480E-23,.1340E-23,.1220E-23,.1100E-23,.1010E-23,
     A .9000E-24/
       data wlgast/175.439,176.991,178.571,180.180,181.818,183.486,
     1  185.185,186.916,188.679,190.476,192.308,194.175,196.078,198.020,
     2  200.000,202.020,204.082/

* auxiliaries
      REAL x, y
      REAL dr
      PARAMETER (dr = pi/180.)
      REAL delO2
      INTEGER i, iw, igast, ierr, icount
      INTEGER iz
      INTEGER ifirst, n

      LOGICAL call1
      SAVE call1, icount
      DATA call1/.TRUE./

*-------------------------------------------------------------------------------



* check, whether user grid is in the O2 absorption band at all...
* if not, set cross section and optical depth values to zero and return

      CALL zero2(dto2,kz,kw)
      CALL zero2(xso2,kz,kw)
      IF (wl(1) .GT. 243.) RETURN

* sec Xhi or Chapman calculation
* for zen > 95 degrees, use zen = 95.  (this is only to compute effective O2
* cross sections. Still, better than setting dto2 = 0. as was done up to 
* version 4.0) sm 1/2000
* In future could replace with mu2(iz) (but mu2 is also wavelength-depenedent)
* or imporved chapman function 

* slant O2 column 

      DO i = 1, nz
         o2col(i) = 0.2095 * scol(i)
      ENDDO

* effective secant of solar zenith angle.  Use 2.0 if no direct sun. 
* For nz, use value at nz-1

      do i = 1, nz - 1
         secchi(i) = scol(i)/vcol(i)
         if(secchi(i) .eq. 0.) secchi(i) = 2.
      enddo
      secchi(nz) = secchi(nz-1)


* read O2 cross section data outside SR-bands only in the very first call
      IF (call1) THEN
************* O2 absorption cross sections:
* from 116 nm to 245 nm, including Schumann-Runge continumm
* from Brasseur and Solomon 1986.

* overwrite from 204 to 241 nm (Herzberg continuum)

        icount=io2
	 
* set values to zero outside the wavelength range defined by the data files

        CALL addpnt(x1,y1,nosr,icount,     x1(1)-deltax,0.)
        CALL addpnt(x1,y1,nosr,icount,               0.,0.)
        CALL addpnt(x1,y1,nosr,icount,x1(icount)+deltax,0.)
        CALL addpnt(x1,y1,nosr,icount,            1.e38,0.)

* set up the internal grid, use full resolution of the cross section data
* outside the SR bands, use Kockarts' grid inside the SR bands
* define Kockarts' grid points

* define the Lyman-Alpha grid
        wlla(1) = 121.4
        wlla(2) = 121.9

* put together the internal grid by "pasting" the Lyman-Alpha grid and 
* Kockarts' grid into the combination of Brasseur/Solomon and JPL grid
        nwint = 0
        DO iw = 1, 9
           nwint = nwint+1
           wlint(nwint) = x1(iw)
        ENDDO
        DO iw = 1, 2
           nwint = nwint+1
           wlint(nwint) = wlla(iw)
        ENDDO
        DO iw = 12, 47
           nwint = nwint + 1
           wlint(nwint) = x1(iw)
        ENDDO
        DO iw = 1, ngast
           nwint = nwint+1
           wlint(nwint) = wlgast(iw)
        ENDDO
        DO iw = 65, 105
           nwint = nwint+1
           wlint(nwint) = x1(iw)
        ENDDO


* interpolate Brasseur/Solomon and JPL data onto internal grid
        CALL inter2(nwint,wlint,xso2int, icount,x1,y1, ierr)

        IF (call1) call1 = .FALSE.

      ENDIF

* if necessary:
* do Kockarts' parameterization of the SR bands, output values of O2
* optical depth and O2 equivalent cross section are on his grid
      IF ((wl(1) .LT. wlgast(ngast)) .AND. 
     >    (wl(nw) .GT. wlgast(1))) THEN
        DO iw = 1, ngast-1
           CALL schu(nz,o2col,secchi,iw,dto2k,xso2k)
        ENDDO
      ENDIF

* do Lyman-Alpha parameterization, output values of O2 opticaldepth
* and O2 effective (equivalent) cross section

      IF ((wl(1) .LE. wlla(nla)) .AND. (wl(nw) .GE. wlla(1))) THEN
         CALL lymana(nz,o2col,secchi,dto2la,xso2la)
      ENDIF

* loop through the altitude levels 
      DO iz = 1, nz

         igast = 0
         delO2 = 0.2095 * cz(iz)    ! vertical O2 column

* loop through the internal wavelength grid
         DO iw = 1, nwint-1

* if outside Kockarts' grid and outside Lyman-Alpha, use the 
* JPL/Brasseur+Solomon data, if inside
* Kockarts' grid, use the parameterized values from the call to SCHU,
* if inside Lyman-Alpha, use the paraemterized values from call to LYMANA
           IF ((wlint(iw+1) .LE. wlgast(1)) .OR.
     >         (wlint(iw) .GE. wlgast(ngast))) THEN
             IF ((wlint(iw+1) .LE. wlla(1)) .OR.
     >           (wlint(iw) .GE. wlla(nla))) THEN
                IF (iz .EQ. nz) THEN
                  dttmp(iw) = 0.
                ELSE
                  dttmp(iw) = xso2int(iw) * delO2
                ENDIF
                xstmp(iw) = xso2int(iw)
             ELSE
                dttmp(iw) = dto2la(iz,1)
                xstmp(iw) = xso2la(iz,1)
             ENDIF
           ELSE
              igast = igast+1
              dttmp(iw) = dto2k(iz,igast)
              xstmp(iw) = xso2k(iz,igast)
           ENDIF

* compute the area in each bin (for correct interpolation purposes only!)
           dttmp(iw) = dttmp(iw) * (wlint(iw+1)-wlint(iw))
           xstmp(iw) = xstmp(iw) * (wlint(iw+1)-wlint(iw))

           if(dttmp(iw).lt.0) then
	    print*,'negative dttmp in seto2 ',iz,iw,dttmp(iw),
     1	      wlint(iw+1),wlint(iw)
	    print*,'dto2k=',igast,dto2k(1:nz,igast)
	    print*,'xso2int,delO2,dto2la=',xso2int(iw),delO2,
     1	      dto2la(iz,1)
            print*,'o2col=',o2col(1:nz)
	    print*,'scol=',scol(1:nz)
	    print*,'z=',z(1:nz)
	    stop
	   endif 
         ENDDO

* interpolate O2 optical depth from the internal grid onto the user grid
         CALL inter3(nw,wl,dtuser, nwint,wlint,dttmp, 0)
         DO iw = 1, nw-1
          dto2(iz,iw) = dtuser(iw)/(wl(iw+1)-wl(iw))
c	    if(.not.(dto2(iz,iw).gt.-1.e36)) dto2(iz,iw)=0.   ! avoid overflow
	    
	  if(dto2(iz,iw).lt.0.or.(.not.abs(dto2(iz,iw)).lt.1e20)) then
	   print*,'wrong dto2 ',iz,iw,dto2(iz,iw)
	   print*,'dtuser=',dtuser(iw)
	   print*,'wl(iw+1),wl(iw)=',wl(iw+1),wl(iw)
	   stop
	  endif  
	    
c	    if(dto2(iz,iw).lt.0) then
c	     print*,'negative dto2 ',iz,iw,dto2(iz,iw),dtuser(iw),
c     1	      wl(iw+1),wl(iw)
c             print*,'cz=',cz(1:nz)
c	     print*,'z=',z(1:nz)
c	     stop
c	    endif 
         ENDDO
      
* interpolate O2 cross section from the internal grid onto the user grid
         CALL inter3(nw,wl,xsuser, nwint,wlint,xstmp, 0)

         DO iw = 1, nw-1
            xso2(iz,iw) = xsuser(iw)/(wl(iw+1)-wl(iw))


         ENDDO

      ENDDO

      RETURN
      END
      
      SUBROUTINE setozo(dobnew,
     $     nz,z,nw,wl,
     $     xso3,s226,s263,s298,tlay,
     $     dto3,o3stem,nzstem)

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Set up an altitude profile of ozone, and corresponding absorption        =*
*=  optical depths.  Subroutine includes a shape-conserving scaling method   =*
*=  that allows scaling of the entire profile to a given overhead ozone      =*
*=  column amount.                                                           =*
*-----------------------------------------------------------------------------*
*=  PARAMETERS:                                                              =*
*=  DOBNEW - REAL, overhead ozone column amount (DU) to which profile     (I)=*
*=           should be scaled.  If DOBNEW < 0, no scaling is done            =*
*=  NZ     - INTEGER, number of specified altitude levels in the working  (I)=*
*=           grid                                                            =*
*=  Z      - REAL, specified altitude working grid (km)                   (I)=*
*=  NW     - INTEGER, number of specified intervals + 1 in working        (I)=*
*=           wavelength grid                                                 =*
*=  WL     - REAL, vector of lower limits of wavelength intervals in      (I)=*
*=           working wavelength grid                                         =*
*=  XSO3   - REAL, molecular absoprtion cross section (cm^2) of O3 at     (I)=*
*=           each specified wavelength (WMO value at 273)                    =*
*=  S226   - REAL, molecular absoprtion cross section (cm^2) of O3 at     (I)=*
*=           each specified wavelength (value from Molina and Molina at 226K)=*
*=  S263   - REAL, molecular absoprtion cross section (cm^2) of O3 at     (I)=*
*=           each specified wavelength (value from Molina and Molina at 263K)=*
*=  S298   - REAL, molecular absoprtion cross section (cm^2) of O3 at     (I)=*
*=           each specified wavelength (value from Molina and Molina at 298K)=*
*=  TLAY   - REAL, temperature (K) at each specified altitude layer       (I)=*
*=  DTO3   - REAL, optical depth due to ozone absorption at each          (O)=*
*=           specified altitude at each specified wavelength                 =*
*-----------------------------------------------------------------------------*
*=  EDIT HISTORY:                                                            =*
*=  02/97  Read in profile from an input file                                =*
*-----------------------------------------------------------------------------*
*= This program is free software;  you can redistribute it and/or modify     =*
*= it under the terms of the GNU General Public License as published by the  =*
*= Free Software Foundation;  either version 2 of the license, or (at your   =*
*= option) any later version.                                                =*
*= The TUV package is distributed in the hope that it will be useful, but    =*
*= WITHOUT ANY WARRANTY;  without even the implied warranty of MERCHANTIBI-  =*
*= LITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public     =*
*= License for more details.                                                 =*
*= To obtain a copy of the GNU General Public License, write to:             =*
*= Free Software Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.   =*
*-----------------------------------------------------------------------------*
*= To contact the authors, please mail to:                                   =*
*= Sasha Madronich, NCAR/ACD, P.O.Box 3000, Boulder, CO, 80307-3000, USA  or =*
*= send email to:  sasha@ucar.edu                                            =*
*-----------------------------------------------------------------------------*
*= Copyright (C) 1994,95,96  University Corporation for Atmospheric Research =*
*-----------------------------------------------------------------------------*

      IMPLICIT NONE
      INCLUDE 'tuv.params'

      INTEGER kdata
      PARAMETER(kdata=150)

********
* input:
********

* grids:
      INTEGER nw,nz,nzstem
      REAL wl(kw),o3stem(nzstem)
      REAL z(kz)
* ozone absorption cross sections at 226, 263, and 298 K:
      REAL xso3(kw), s226(kw),s263(kw),s298(kw)
      REAL dobnew
* mid-layer temperature:
      REAL tlay(kz) 

********
* output:
********
      REAL dto3(kz,kw)


********
* local:
********

      REAL cz(kz)

* ozone profile data:
      integer nd
      parameter(nd=39)
      REAL zd(kdata), o3(kdata)
      data zd(1:nd)/0., 1., 2., 4., 6., 8.,10.,12.,14.,
     1  16.,18.,20.,22.,24.,26.,28.,30.,32.,
     2  34.,36.,38.,40.,42.,44.,46.,48.,50.,
     3  52.,54.,56.,58.,60.,62.,64.,66.,68.,70.,72.,74./
     
       data o3(1:nd)/.102E+13,.920E+12,.680E+12,.580E+12,.570E+12,
     1  .650E+12,.113E+13,.202E+13,.235E+13,.295E+13,.404E+13,.477E+13,
     2  .486E+13,.454E+13,.403E+13,.324E+13,.252E+13,.203E+13,
     3  .158E+13,.122E+13,.873E+12,.607E+12,.398E+12,.274E+12,
     4  .169E+12,.103E+12,.664E+11,.384E+11,.255E+11,.161E+11,
     5 .112E+11,.733E+10,.481E+10,.317E+10,.172E+10,.750E+09,
     6 .540E+09,.220E+09,.170E+09/
     
      REAL cd(kdata)
      REAL hscale
      REAL dobold, scale, dobstem, doboldstem
      REAL colold, colnew
      REAL so3

* other:
      INTEGER i, iw

********
* External functions:
********
      REAL fsum
      EXTERNAL fsum
*_______________________________________________________________________


* read in ozone profile

c      WRITE(kout,*) 'ozone profile: USSA, 1976'

* compute column increments

      DO 11, i = 1, nd - 1
         cd(i) = (o3(i+1)+o3(i)) * 1.E5 * (zd(i+1)-zd(i)) / 2. 
   11 CONTINUE

* Include exponential tail integral from infinity to 50 km,
* fold tail integral into top layer
* specify scale height near top of data.

      hscale = 4.50e5
      cd(nd-1) = cd(nd-1) + hscale * o3(nd)

* alternative input ozone concentration data could include, e.g., 
* a read file here:

***********
*********** end data input.

* Compute column increments on standard z-grid.  

      CALL inter3(nz,z,cz, nd,zd,cd, 1)
       
* scale values of cz(i) by any dobson unit

      colold = fsum(nz-1,cz)
      dobold = colold/2.687e16
      doboldstem = fsum(nzstem-1,cz)/2.687e16    ! dobson depth below STEM height
c      WRITE(kout,100) colold, dobold
  100 FORMAT(5x,'old O3 Column = ', 1pe11.4,1x,'# cm-2  = ',
     $     0pf8.2, '  Dobson Units ')

      DO i = 1, nzstem - 1
         cz(i) = (o3stem(i+1)+o3stem(i)) * 1.E5 *  ! load ozone from STEM
     1	 (z(i+1)-z(i)) / 2. 
      enddo
      dobstem=fsum(nzstem-1,cz)/ 2.687e16        ! dobson depth of STEM ozone
            
c      if (dobnew .lt. 0.) then
c         scale = 1.
c      else
c         scale = dobnew/dobold
c      endif

      if((dobnew-dobstem).le.0) then
        scale = 1.
      else
        scale=(dobnew-dobstem)/(dobold-doboldstem)
      endif
        	
      do i = nzstem, nz-1
         cz(i) = cz(i) * scale
      enddo

      colnew = fsum(nz-1,cz)
c      WRITE(kout,105) colnew, colnew/2.687E16
  105 format(5x,'new O3 Column = ', 1pe11.4,1x,'# cm-2  = ',
     $     0pf8.2, '  Dobson Units ')

************************************
* calculate ozone optical depth for each layer, with temperature 
* correction.  Output, dto3(kz,kw)

      DO 20, iw = 1, nw-1
         so3 = xso3(iw)
         DO 10, i = 1, nz - 1

            IF ( wl(iw) .GT. 240.5  .AND. wl(iw+1) .LT. 350. ) THEN
               IF (tlay(i) .LT. 263.) THEN
                  so3 = s226(iw) + (s263(iw)-s226(iw)) / (263.-226.) *
     $                 (tlay(i)-226.)
               ELSE
                  so3 = s263(iw) + (s298(iw)-s263(iw)) / (298.-263.) *
     $              (tlay(i)-263.)
               ENDIF
            ENDIF

            dto3(i,iw) = cz(i)*so3

   10    CONTINUE
   20 CONTINUE
*_______________________________________________________________________

      RETURN
      END
      SUBROUTINE setso2(so2new,
     $     nz,z,nw,wl,
     $     xsso2, tlay,
     $     dtso2,so2stem,nzstem)

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Set up an altitude profile of SO2 molecules, and corresponding absorption=*
*=  optical depths.  Subroutine includes a shape-conserving scaling method   =*
*=  that allows scaling of the entire profile to a given overhead SO2        =*
*=  column amount.                                                           =*
*-----------------------------------------------------------------------------*
*=  PARAMETERS:                                                              =*
*=  SO2NEW - REAL, overhead SO2 column amount (molec/cm^2) to which       (I)=*
*=           profile should be scaled.  If SO2NEW < 0, no scaling is done    =*
*=  NZ     - INTEGER, number of specified altitude levels in the working  (I)=*
*=           grid                                                            =*
*=  Z      - REAL, specified altitude working grid (km)                   (I)=*
*=  NW     - INTEGER, number of specified intervals + 1 in working        (I)=*
*=           wavelength grid                                                 =*
*=  WL     - REAL, vector of lower limits of wavelength intervals in      (I)=*
*=           working wavelength grid                                         =*
*=  XSSO2  - REAL, molecular absoprtion cross section (cm^2) of O2 at     (I)=*
*=           each specified wavelength                                       =*
*=  TLAY   - REAL, temperature (K) at each specified altitude layer       (I)=*
*=  DTSO2  - REAL, optical depth due to SO2 absorption at each            (O)=*
*=           specified altitude at each specified wavelength                 =*
*-----------------------------------------------------------------------------*
*=  EDIT HISTORY:                                                            =*
*=  Original                                                                 =*
*-----------------------------------------------------------------------------*
*= This program is free software;  you can redistribute it and/or modify     =*
*= it under the terms of the GNU General Public License as published by the  =*
*= Free Software Foundation;  either version 2 of the license, or (at your   =*
*= option) any later version.                                                =*
*= The TUV package is distributed in the hope that it will be useful, but    =*
*= WITHOUT ANY WARRANTY;  without even the implied warranty of MERCHANTIBI-  =*
*= LITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public     =*
*= License for more details.                                                 =*
*= To obtain a copy of the GNU General Public License, write to:             =*
*= Free Software Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.   =*
*-----------------------------------------------------------------------------*
*= To contact the authors, please mail to:                                   =*
*= Sasha Madronich, NCAR/ACD, P.O.Box 3000, Boulder, CO, 80307-3000, USA  or =*
*= send email to:  sasha@ucar.edu                                            =*
*-----------------------------------------------------------------------------*
*= Copyright (C) 1994,95,96  University Corporation for Atmospheric Research =*
*-----------------------------------------------------------------------------*

      IMPLICIT NONE
      INCLUDE 'tuv.params'

      INTEGER kdata
      PARAMETER(kdata=51)

********
* input:
********

* grids:

      REAL wl(kw)
      REAL z(kz)
      INTEGER nw, nz, nzstem
      REAL so2new,so2stem(nzstem)

* mid-layer temperature:

      REAL tlay(kz) 

********
* output:
********

      REAL dtso2(kz,kw)

********
* local:
********

* absorption cross sections 

      REAL xsso2(kw)
      REAL cz(kz)

* sulfur dioxide profile data:

      REAL zd(kdata), so2(kdata)
      REAL cd(kdata)
      REAL hscale
      REAL colold, colnew
      REAL scale
      REAL sso2

* other:

      INTEGER i, l, nd

********
* External functions:
********
      REAL fsum
      EXTERNAL fsum

*_______________________________________________________________________
* Data input:

* Example:  set to 1 ppb in lowest 1 km, set to zero above that.
* - do by specifying concentration at 3 altitudes.

c      write(kout,*) 'SO2:  1 ppb in lowest 1 km, 0 above'

      nd = 3
      zd(1) = 0.
      so2(1) = 1. * 2.69e10

      zd(2) = 1.
      so2(2) = 1. * 2.69e10

      zd(3) = zd(2)* 1.000001
      so2(3) = 0.

C     zd(4) = zd(3)*1.1 
C     so2(4) = 0.

* compute column increments (alternatively, can specify these directly)

      DO 11, i = 1, nd - 1
         cd(i) = (so2(i+1)+so2(i)) * 1.E5 * (zd(i+1)-zd(i)) / 2. 
   11 CONTINUE

* Include exponential tail integral from top level to infinity.
* fold tail integral into top layer
* specify scale height near top of data (use ozone value)

      hscale = 4.50e5
      cd(nd-1) = cd(nd-1) + hscale * so2(nd)

***********
*********** end data input.

* Compute column increments on standard z-grid.  

      CALL inter3(nz,z,cz, nd,zd,cd, 1)

* scale values of cz(i) 

      colold = fsum(nz-1,cz)
c      WRITE(kout,100) colold, colold/2.687E16
  100 FORMAT(5x,'old SO2 Column = ', 1pe11.4,1x,'# cm-2  = ',
     $     0pf8.2, '  Dobson Units ')

      if ( (so2new .lt. 0.)  .or.  (colold .le. 0.) ) then
         scale = 1.
      else
         scale =  2.687e16*so2new/colold
      endif

      do i = 1, nz-1
         cz(i) = cz(i) * scale
      enddo
      DO i = 1, nzstem - 1
         cz(i) = (so2stem(i+1)+so2stem(i)) * 1.E5 *  ! load SO2 from STEM
     1	 (z(i+1)-z(i)) / 2. 
      enddo

      colnew = fsum(nz-1,cz)
c      WRITE(kout,105) colnew, colnew/2.687E16
  105 format(5x,'new SO2 Column = ', 1pe11.4,1x,'# cm-2  = ',
     $     0pf8.2, '  Dobson Units ')

************************************
* calculate sulfur optical depth for each layer, with temperature 
* correction.  Output, dtso2(kz,kw)

      DO 20, l = 1, nw-1
         sso2 = xsso2(l)
         DO 10, i = 1, nz - 1

c Leaving this part in in case i want to interpolate between 
c the 221K and 298K data.
c
c            IF ( wl(l) .GT. 240.5  .AND. wl(l+1) .LT. 350. ) THEN
c               IF (tlay(i) .LT. 263.) THEN
c                  sso2 = s221(l) + (s263(l)-s226(l)) / (263.-226.) *
c     $                 (tlay(i)-226.)
c               ELSE
c                  sso2 = s263(l) + (s298(l)-s263(l)) / (298.-263.) *
c     $              (tlay(i)-263.)
c               ENDIF
c            ENDIF

            dtso2(i,l) = cz(i)*sso2

   10    CONTINUE
   20 CONTINUE
*_______________________________________________________________________

      RETURN
      END

      SUBROUTINE settmp(nz,z,tlev,tlay,tempstem,nzstem)

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Set up an altitude profile of temperatures.  Temperature values are      =*
*=  needed to compute some cross sections and quantum yields.  Distinguish   =*
*=  between temperature at levels and layers.                                =*
*-----------------------------------------------------------------------------*
*=  PARAMETERS:                                                              =*
*=  NZ      - INTEGER, number of specified altitude levels in the working (I)=*
*=            grid                                                           =*
*=  Z       - REAL, specified altitude working grid (km)                  (I)=*
*=  TLEV    - REAL, temperature (K) at each specified altitude level      (O)=*
*=  TLAY    - REAL, temperature (K) at each specified altitude layer      (O)=*
*-----------------------------------------------------------------------------*
*=  EDIT HISTORY:                                                            =*
*=  02/97  Read in profile from an input file                                =*
*-----------------------------------------------------------------------------*
*= This program is free software;  you can redistribute it and/or modify     =*
*= it under the terms of the GNU General Public License as published by the  =*
*= Free Software Foundation;  either version 2 of the license, or (at your   =*
*= option) any later version.                                                =*
*= The TUV package is distributed in the hope that it will be useful, but    =*
*= WITHOUT ANY WARRANTY;  without even the implied warranty of MERCHANTIBI-  =*
*= LITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public     =*
*= License for more details.                                                 =*
*= To obtain a copy of the GNU General Public License, write to:             =*
*= Free Software Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.   =*
*-----------------------------------------------------------------------------*
*= To contact the authors, please mail to:                                   =*
*= Sasha Madronich, NCAR/ACD, P.O.Box 3000, Boulder, CO, 80307-3000, USA  or =*
*= send email to:  sasha@ucar.edu                                            =*
*-----------------------------------------------------------------------------*
*= Copyright (C) 1994,95,96  University Corporation for Atmospheric Research =*
*-----------------------------------------------------------------------------*

      IMPLICIT NONE
      INCLUDE 'tuv.params'

      INTEGER kdata
      PARAMETER(kdata=150)

* input: (altitude working grid)
      INTEGER nz,nzstem
      REAL z(kz),tempstem(nzstem)


* output:
      REAL tlev(kz), tlay(kz)

* local:
      integer  i,nd
      parameter(nd=121)
      REAL zd(kdata), td(kdata)
      data zd(1:nd)/  0.,  1.,  2.,  3.,  4.,  5.,  6.,  7.,  8.,
     1  9., 10., 11., 12., 13., 14., 15., 16., 17.,
     2  18., 19., 20., 21., 22., 23., 24., 25., 26.,
     3  27., 28., 29., 30., 31., 32., 33., 34., 35.,
     4  36., 37., 38., 39., 40., 41., 42., 43., 44.,
     5  45., 46., 47., 48., 49., 50., 51., 52., 53.,
     6  54., 55., 56., 57., 58., 59., 60., 61., 62.,
     7  63., 64., 65., 66., 67., 68., 69., 70., 71.,
     8  72., 73., 74., 75., 76., 77., 78., 79., 80.,
     9  81., 82., 83., 84., 85., 86., 87., 88., 89.,
     a  90., 91., 92., 93., 94., 95., 96., 97., 98.,
     b  99.,100.,101.,102.,103.,104.,105.,106.,107.,
     c 108.,109.,110.,111.,112.,113.,114.,115.,116.,
     d 117.,118.,119.,120./
       data td(1:nd)/288.150,281.651,275.154,268.659,262.166,255.676,
     1 249.187,242.700,236.215,229.733,223.252,216.774,216.650,216.650,
     2 216.650,216.650,216.650,216.650,216.650,216.650,216.650,
     3 217.581,218.574,219.567,220.560,221.552,222.544,223.536,
     4 224.527,225.518,226.509,227.500,228.490,230.973,233.743,
     5 236.513,239.282,242.050,244.818,247.584,250.350,253.114,
     6 255.878,258.641,261.403,264.164,266.925,269.684,270.650,
     7 270.650,270.650,270.650,269.031,266.277,263.524,260.771,
     8 258.019,255.268,252.518,249.769,247.021,244.274,241.524,
     9 238.781,236.036,233.292,230.549,227.807,225.065,222.325,
     a 219.585,216.846,214.263,212.308,210.353,208.399,206.446,
     b 204.493,202.541,200.590,198.639,196.688,194.739,192.790,
     c 190.841,188.893,186.870,186.870,186.870,186.870,186.870,
     d 186.870,186.960,187.250,187.740,188.420,189.310,190.400,
     e 191.720,193.280,195.080,197.160,199.530,202.230,205.310,
     f 208.840,212.890,217.630,223.290,230.330,240.000,252.000,
     g 264.000,276.000,288.000,300.000,312.000,324.000,336.000,
     h 348.000,360.000/

      
*_______________________________________________________________________


* read in temperature profile


* use constant temperature to infinity:  

      zd(nd) = 1.E10

* alternative input temperature data could include, e.g., a read file here:

***********
*********** end data input.

* interpolate onto z-grid

      CALL inter1(nz,z,tlev,nd,zd,td)

      tlev(1:nzstem)=tempstem(1:nzstem)  ! replace the lower layer with STEM temperature
* compute layer-averages

      DO 20, i = 1, nz - 1
         tlay(i) = (tlev(i+1) + tlev(i))/2.
 20   CONTINUE
*_______________________________________________________________________
      
      RETURN
      END

       SUBROUTINE sjo2(nz,nw,xso2,nj,sq)

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Update the weighting function (cross section x quantum yield) for O2     =*
*=  photolysis.  Effective O2 cross section depends upon solar zenith angle  =*
*=  (Schuman-Runge bands), so needs to be update at each time step.          =*
*-----------------------------------------------------------------------------*
*=  PARAMETERS:                                                              =*
*=  NZ     - INTEGER, number of altitude levels in working altitude grid  (I)=*
*=  NW     - INTEGER, number of specified intervals + 1 in working        (I)=*
*=           wavelength grid                                                 =*
*=  XSO2   - REAL, molecular absorption cross section in SR bands at      (I)=*
*=           each specified altitude and wavelength.  Includes Herzberg      =*
*=            continuum.                                                     =*
*=  NJ     - INTEGER, index of O2 photolysis in array SQ                  (I)=*
*=  SQ     - REAL, cross section x quantum yield (cm^2) for each          (O)=*
*=           photolysis reaction defined, at each defined wavelength and     =*
*=           at each defined altitude level                                  =*
*-----------------------------------------------------------------------------*
*=  EDIT HISTORY:                                                            =*
*=  08/96  Original                                                          =*
*-----------------------------------------------------------------------------*
*= This program is free software;  you can redistribute it and/or modify     =*
*= it under the terms of the GNU General Public License as published by the  =*
*= Free Software Foundation;  either version 2 of the license, or (at your   =*
*= option) any later version.                                                =*
*= The TUV package is distributed in the hope that it will be useful, but    =*
*= WITHOUT ANY WARRANTY;  without even the implied warranty of MERCHANTIBI-  =*
*= LITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public     =*
*= License for more details.                                                 =*
*= To obtain a copy of the GNU General Public License, write to:             =*
*= Free Software Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.   =*
*-----------------------------------------------------------------------------*
*= To contact the authors, please mail to:                                   =*
*= Sasha Madronich, NCAR/ACD, P.O.Box 3000, Boulder, CO, 80307-3000, USA  or =*
*= send email to:  sasha@ucar.edu                                            =*
*-----------------------------------------------------------------------------*
*= Copyright (C) 1994,95,96  University Corporation for Atmospheric Research =*
*-----------------------------------------------------------------------------*

      IMPLICIT NONE
      INCLUDE 'tuv.params'

* calling parameters

      INTEGER nz, nw, nj
      REAL xso2(kz,kw)
      REAL sq(kj,kz,kw)

* local

      INTEGER iw, iz
*______________________________________________________________________________

* O2 + hv -> O + O
* quantum yield assumed to be unity
* assign cross section values at all wavelengths and at all altitudes
*      qy = 1.

      DO iw = 1, nw-1
        DO iz = 1, nz
          sq(nj,iz,iw) = xso2(iz,iw)
        ENDDO
      ENDDO
*______________________________________________________________________________


      RETURN
      END
      SUBROUTINE sphers(nz, z, zen, dsdh, nid)

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Calculate slant path over vertical depth ds/dh in spherical geometry.    =*
*=  Calculation is based on:  A.Dahlback, and K.Stamnes, A new spheric model =*
*=  for computing the radiation field available for photolysis and heating   =*
*=  at twilight, Planet.Space Sci., v39, n5, pp. 671-683, 1991 (Appendix B)  =*
*-----------------------------------------------------------------------------*
*=  PARAMETERS:                                                              =*
*=  NZ      - INTEGER, number of specified altitude levels in the working (I)=*
*=            grid                                                           =*
*=  Z       - REAL, specified altitude working grid (km)                  (I)=*
*=  ZEN     - REAL, solar zenith angle (degrees)                          (I)=*
*=  DSDH    - REAL, slant path of direct beam through each layer crossed  (O)=*
*=            when travelling from the top of the atmosphere to layer i;     =*
*=            DSDH(i,j), i = 0..NZ-1, j = 1..NZ-1                            =*
*=  NID     - INTEGER, number of layers crossed by the direct beam when   (O)=*
*=            travelling from the top of the atmosphere to layer i;          =*
*=            NID(i), i = 0..NZ-1                                            =*
*-----------------------------------------------------------------------------*
*=  EDIT HISTORY:                                                            =*
*=  Original                                                                 =*
*-----------------------------------------------------------------------------*
*= This program is free software;  you can redistribute it and/or modify     =*
*= it under the terms of the GNU General Public License as published by the  =*
*= Free Software Foundation;  either version 2 of the license, or (at your   =*
*= option) any later version.                                                =*
*= The TUV package is distributed in the hope that it will be useful, but    =*
*= WITHOUT ANY WARRANTY;  without even the implied warranty of MERCHANTIBI-  =*
*= LITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public     =*
*= License for more details.                                                 =*
*= To obtain a copy of the GNU General Public License, write to:             =*
*= Free Software Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.   =*
*-----------------------------------------------------------------------------*
*= To contact the authors, please mail to:                                   =*
*= Sasha Madronich, NCAR/ACD, P.O.Box 3000, Boulder, CO, 80307-3000, USA  or =*
*= send email to:  sasha@ucar.edu                                            =*
*-----------------------------------------------------------------------------*
*= Copyright (C) 1994,95,96  University Corporation for Atmospheric Research =*
*-----------------------------------------------------------------------------*

      IMPLICIT NONE
      INCLUDE 'tuv.params'

* input
      INTEGER nz
      REAL zen, z(kz)

* output
      INTEGER nid(0:kz)
      REAL dsdh(0:kz,kz)

* more program constants
      REAL re, ze(kz)
      REAL  dr
      PARAMETER ( dr = pi/180.)

* local 

      REAL zenrad, rpsinz, rj, rjp1, dsj, dhj, ga, gb, sm
      INTEGER i, j, k
      INTEGER id

      INTEGER nlayer
      REAL zd(0:kz-1)

*-----------------------------------------------------------------------------

      zenrad = zen*dr

* number of layers:
      nlayer = nz - 1

* include the elevation above sea level to the radius of the earth:
      re = radius + z(1)
* correspondingly z changed to the elevation above earth surface:
      DO k = 1, nz
         ze(k) = z(k) - z(1)
      END DO

* inverse coordinate of z
      zd(0) = ze(nz)
      DO k = 1, nlayer
        zd(k) = ze(nz - k)
      END DO

* initialize dsdh(i,j), nid(i)
      DO i = 0, kz
       nid(i) = 0
       DO j = 1, kz
        dsdh(i,j) = 0.
       END DO
      END DO

* calculate ds/dh of every layer
      DO 100 i = 0, nlayer

        rpsinz = (re + zd(i)) * SIN(zenrad)
 
        IF ( (zen .GT. 90.0) .AND. (rpsinz .LT. re) ) THEN
           nid(i) = -1
        ELSE

*
* Find index of layer in which the screening height lies
*
           id = i 
           IF( zen .GT. 90.0 ) THEN
              DO 10 j = 1, nlayer
                 IF( (rpsinz .LT. ( zd(j-1) + re ) ) .AND.
     $               (rpsinz .GE. ( zd(j) + re )) ) id = j
 10           CONTINUE
           END IF
 
           DO 20 j = 1, id

             sm = 1.0
             IF(j .EQ. id .AND. id .EQ. i .AND. zen .GT. 90.0)
     $          sm = -1.0
 
             rj = re + zd(j-1)
             rjp1 = re + zd(j)
 
             dhj = zd(j-1) - zd(j)
 
             ga = rj*rj - rpsinz*rpsinz
             gb = rjp1*rjp1 - rpsinz*rpsinz
             IF (ga .LT. 0.0) ga = 0.0
             IF (gb .LT. 0.0) gb = 0.0
 
             IF(id.GT.i .AND. j.EQ.id) THEN
                dsj = SQRT( ga )
             ELSE
                dsj = SQRT( ga ) - sm*SQRT( gb )
             END IF
             dsdh(i,j) = dsj / dhj

 20        CONTINUE
 
           nid(i) = id
 
        END IF

 100  CONTINUE

*-----------------------------------------------------------------------------

      RETURN
      END
      SUBROUTINE sundis(idate,esrm2)

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Calculate Earth-Sun distance variation for a given date.  Based on       =*
*=  Fourier coefficients originally from:  Spencer, J.W., 1971, Fourier      =*
*=  series representation of the position of the sun, Search, 2:172          =*
*-----------------------------------------------------------------------------*
*=  PARAMETERS:                                                              =*
*=  IDATE  - INTEGER, specification of the date, from YYMMDD              (I)=*
*=  ESRM2  - REAL, variation of the Earth-sun distance                    (O)=*
*=           ESRM2 = (average e/s dist)^2 / (e/s dist on day IDATE)^2        =*
*-----------------------------------------------------------------------------*
*=  EDIT HISTORY:                                                            =*
*=  01/95  Changed computation of trig function values                       =*
*-----------------------------------------------------------------------------*
*= This program is free software;  you can redistribute it and/or modify     =*
*= it under the terms of the GNU General Public License as published by the  =*
*= Free Software Foundation;  either version 2 of the license, or (at your   =*
*= option) any later version.                                                =*
*= The TUV package is distributed in the hope that it will be useful, but    =*
*= WITHOUT ANY WARRANTY;  without even the implied warranty of MERCHANTIBI-  =*
*= LITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public     =*
*= License for more details.                                                 =*
*= To obtain a copy of the GNU General Public License, write to:             =*
*= Free Software Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.   =*
*-----------------------------------------------------------------------------*
*= To contact the authors, please mail to:                                   =*
*= Sasha Madronich, NCAR/ACD, P.O.Box 3000, Boulder, CO, 80307-3000, USA  or =*
*= send email to:  sasha@ucar.edu                                            =*
*-----------------------------------------------------------------------------*
*= Copyright (C) 1994,95,96  University Corporation for Atmospheric Research =*
*-----------------------------------------------------------------------------*

      IMPLICIT NONE

* input:
      INTEGER idate

* output:
      REAL esrm2

* internal:
      INTEGER iyear, imonth, iday, mday, month, jday
      REAL dayn, thet0
      REAL sinth, costh, sin2th, cos2th
      INTEGER imn(12)

      REAL pi
      PARAMETER(pi=3.1415926535898)
*_______________________________________________________________________

      DATA imn/31,28,31,30,31,30,31,31,30,31,30,31/             
*_______________________________________________________________________

* parse date to find day number (Julian day)
cgy
      if (idate.lt.1.or.idate.gt.991231) then
         write(*,*) 'Date must be between 000001 and 991231'
         write(*,*) 'date = ', idate
         stop
      endif
cgy

      iyear = int(idate/10000)
      imonth = int( (idate-10000*iyear)/100 )
      iday = idate - (10000*iyear + 100*imonth)

cgy
      if (imonth.gt.12) then
         write(*,*) 'Month in date exceeds 12'
         write(*,*) 'date = ', idate
         write(*,*) 'month = ', imonth
         stop
      endif
cgy

      IF ( MOD(iyear,4) .EQ. 0) THEN
         imn(2) = 29
      ELSE
         imn(2) = 28
      ENDIF

cgy
      if (iday.gt.imn(imonth)) then
         write(*,*) 'Day in date exceeds days in month'
         write(*,*) 'date = ', idate
         write(*,*) 'day = ', iday
         stop
      endif
cgy

      mday = 0
      DO 12, month = 1, imonth-1
         mday = mday + imn(month)	  	   
   12 CONTINUE
      jday = mday + iday
      dayn = FLOAT(jday - 1) + 0.5

* define angular day number and compute esrm2:

      thet0 = 2.*pi*dayn/365.

* calculate SIN(2*thet0), COS(2*thet0) from
* addition theoremes for trig functions for better
* performance;  the computation of sin2th, cos2th
* is about 5-6 times faster than the evaluation
* of the intrinsic functions SIN and COS
*
      sinth = SIN(thet0)
      costh = COS(thet0)
      sin2th = 2.*sinth*costh
      cos2th = costh*costh - sinth*sinth
      esrm2  = 1.000110 + 
     $         0.034221*costh  +  0.001280*sinth + 
     $         0.000719*cos2th +  0.000077*sin2th
*_______________________________________________________________________

      RETURN
      END

      SUBROUTINE zenith(lat,long,idate,ut,azim,zen)

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Calculate solar zenith angle and azimuth for a given time and location.  =*
*=  Calculation is based on equations given in:  Paltridge and Platt, Radia- =*
*=  tive Processes in Meteorology and Climatology, Elsevier, pp. 62,63, 1976.=*
*=  Fourier coefficients originally from:  Spencer, J.W., 1971, Fourier      =*
*=  series representation of the position of the sun, Search, 2:172.         =*
*=  Note:  This approximate program does not account fro changes from year   =*
*=  to year.                                                                 =*
*-----------------------------------------------------------------------------*
*=  PARAMETERS:                                                              =*
*=  LAT   - REAL, latitude of location (degrees)                          (I)=*
*=  LONG  - REAL, longitude of location (degrees)                         (I)=*
*=  IDATE - INTEGER, date in the form YYMMDD                              (I)=*
*=  UT    - REAL, local time in decimal UT (e.g., 16.25 means 15 minutes  (I)=*
*=          after 4 pm)                                                      =*
*=  AZIM  - REAL, azimuth (degrees)                                       (O)=*
*=  ZEN   - REAL, solar zenith angle (degrees)                            (O)=*
*-----------------------------------------------------------------------------*
*=  EDIT HISTORY:                                                            =*
*=  Original                                                                 =*
*-----------------------------------------------------------------------------*
*= This program is free software;  you can redistribute it and/or modify     =*
*= it under the terms of the GNU General Public License as published by the  =*
*= Free Software Foundation;  either version 2 of the license, or (at your   =*
*= option) any later version.                                                =*
*= The TUV package is distributed in the hope that it will be useful, but    =*
*= WITHOUT ANY WARRANTY;  without even the implied warranty of MERCHANTIBI-  =*
*= LITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public     =*
*= License for more details.                                                 =*
*= To obtain a copy of the GNU General Public License, write to:             =*
*= Free Software Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.   =*
*-----------------------------------------------------------------------------*
*= To contact the authors, please mail to:                                   =*
*= Sasha Madronich, NCAR/ACD, P.O.Box 3000, Boulder, CO, 80307-3000, USA  or =*
*= send email to:  sasha@ucar.edu                                            =*
*-----------------------------------------------------------------------------*
*= Copyright (C) 1994,95,96  University Corporation for Atmospheric Research =*
*-----------------------------------------------------------------------------*

      IMPLICIT NONE

* input:
      REAL lat,long
      REAL ut
      INTEGER idate

* output:
      REAL azim, zen

* local:
      REAL lbut,lzut
      REAL rlt
      REAL d, tz, rdecl, eqr, eqh, zpt
      REAL csz, zr, caz, raz 
      REAL sintz, costz, sin2tz, cos2tz, sin3tz, cos3tz

      INTEGER iiyear, imth, iday, ijd
      INTEGER imn(12)

      INTEGER i

* program constants:

      REAL pi, dr
      PARAMETER(pi=3.1415926535898)
      PARAMETER (dr=pi/180.D0)
*_______________________________________________________________________

      DATA imn/31,28,31,30,31,30,31,31,30,31,30,31/             
*_______________________________________________________________________

* convert to radians

      rlt = lat*dr

* parse date

      iiyear = idate/10000
      imth = (idate - iiyear*10000)/100
      iday = idate - iiyear*10000 - imth*100

* identify and correct leap years

      IF (MOD(iiyear,4) .EQ. 0) THEN
         imn(2) = 29
      ELSE
         imn(2) = 28
      ENDIF

* compute current (Julian) day of year IJD = 1 to 365

      ijd = 0
      DO 30, i = 1, imth - 1
         ijd = ijd + imn(i)
   30 CONTINUE
      ijd = ijd + iday

* calculate decimal Julian day from start of year:

      d = FLOAT(ijd-1) + ut/24.

* Equation 3.8 for "day-angle"

      tz = 2.*pi*d/365.

* Calculate sine and cosine from addition theoremes for 
* better performance;  the computation of sin2tz,
* sin3tz, cos2tz and cos3tz is about 5-6 times faster
* than the evaluation of the intrinsic functions 
*
* It is SIN(x+y) = SIN(x)*COS(y)+COS(x)*SIN(y)
* and   COS(x+y) = COS(x)*COS(y)-SIN(x)*SIN(y)
*
* sintz  = SIN(tz)      costz  = COS(tz)
* sin2tz = SIN(2.*tz)   cos2tz = SIN(2.*tz)
* sin3tz = SIN(3.*tz)   cos3tz = COS(3.*tz)
*
      sintz = SIN(tz)
      costz = COS(tz)
      sin2tz = 2.*sintz*costz
      cos2tz = costz*costz-sintz*sintz
      sin3tz = sintz*cos2tz + costz*sin2tz
      cos3tz = costz*cos2tz - sintz*sin2tz

* Equation 3.7 for declination in radians

      rdecl = 0.006918 - 0.399912*costz  + 0.070257*sintz 
     $                 - 0.006758*cos2tz + 0.000907*sin2tz    
     $                 - 0.002697*cos3tz + 0.001480*sin3tz

* Equation 3.11 for Equation of time  in radians

      eqr   = 0.000075 + 0.001868*costz  - 0.032077*sintz
     $		       - 0.014615*cos2tz - 0.040849*sin2tz

* convert equation of time to hours:

      eqh = eqr*24./(2.*pi) 

* calculate local hour angle (hours):

      lbut = 12. - eqh - long*24./360 

* convert to angle from UT

      lzut = 15.*(ut - lbut)
      zpt = lzut*dr

* Equation 2.4 for cosine of zenith angle 

      csz = SIN(rlt)*SIN(rdecl) + COS(rlt)*COS(rdecl)*COS(zpt)
      zr = ACOS(csz)
      zen = zr/dr

*   calc local solar azimuth

      caz = (SIN(rdecl) - SIN(rlt)*COS(zr))/(COS(rlt)*SIN(zr))
      raz = ACOS(caz)
      azim = raz/dr
*_______________________________________________________________________

      RETURN
      END
      subroutine zero1(x,m)

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Initialize all elements of a floating point vector with zero.            =*
*-----------------------------------------------------------------------------*
*=  PARAMETERS:                                                              =*
*=  X  - REAL, vector to be initialized                                   (O)=*
*=  M  - INTEGER, number of elements in X                                 (I)=*
*-----------------------------------------------------------------------------*
*=  EDIT HISTORY:                                                            =*
*=  Original                                                                 =*
*-----------------------------------------------------------------------------*
*= This program is free software;  you can redistribute it and/or modify     =*
*= it under the terms of the GNU General Public License as published by the  =*
*= Free Software Foundation;  either version 2 of the license, or (at your   =*
*= option) any later version.                                                =*
*= The TUV package is distributed in the hope that it will be useful, but    =*
*= WITHOUT ANY WARRANTY;  without even the implied warranty of MERCHANTIBI-  =*
*= LITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public     =*
*= License for more details.                                                 =*
*= To obtain a copy of the GNU General Public License, write to:             =*
*= Free Software Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.   =*
*-----------------------------------------------------------------------------*
*= To contact the authors, please mail to:                                   =*
*= Sasha Madronich, NCAR/ACD, P.O.Box 3000, Boulder, CO, 80307-3000, USA  or =*
*= send email to:  sasha@ucar.edu                                            =*
*-----------------------------------------------------------------------------*
*= Copyright (C) 1994,95,96  University Corporation for Atmospheric Research =*
*-----------------------------------------------------------------------------*

      implicit none
      integer i, m
      real x(m)
      do 1 i = 1, m
         x(i) = 0.
 1    continue
      return
      end

      subroutine zero2(x,m,n)

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Initialize all elements of a 2D floating point array with zero.          =*
*-----------------------------------------------------------------------------*
*=  PARAMETERS:                                                              =*
*=  X  - REAL, array to be initialized                                    (O)=*
*=  M  - INTEGER, number of elements along the first dimension of X,      (I)=*
*=       exactly as specified in the calling program                         =*
*=  N  - INTEGER, number of elements along the second dimension of X,     (I)=*
*=       exactly as specified in the calling program                         =*
*-----------------------------------------------------------------------------*
*=  EDIT HISTORY:                                                            =*
*=  Original                                                                 =*
*-----------------------------------------------------------------------------*
*= This program is free software;  you can redistribute it and/or modify     =*
*= it under the terms of the GNU General Public License as published by the  =*
*= Free Software Foundation;  either version 2 of the license, or (at your   =*
*= option) any later version.                                                =*
*= The TUV package is distributed in the hope that it will be useful, but    =*
*= WITHOUT ANY WARRANTY;  without even the implied warranty of MERCHANTIBI-  =*
*= LITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public     =*
*= License for more details.                                                 =*
*= To obtain a copy of the GNU General Public License, write to:             =*
*= Free Software Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.   =*
*-----------------------------------------------------------------------------*
*= To contact the authors, please mail to:                                   =*
*= Sasha Madronich, NCAR/ACD, P.O.Box 3000, Boulder, CO, 80307-3000, USA  or =*
*= send email to:  sasha@ucar.edu                                            =*
*-----------------------------------------------------------------------------*
*= Copyright (C) 1994,95,96  University Corporation for Atmospheric Research =*
*-----------------------------------------------------------------------------*

      implicit none
* m,n : dimensions of x, exactly as specified in the calling program
      integer i, j, m, n
      real x(m,n)
      do 1 j = 1, n
         do 2 i = 1, m
            x(i,j) = 0.
 2       continue
 1    continue
      return
      end
