56d55
< 			    ! 7 is the Single Scatter Albedo
146d144
<       save first
166c164
<        do while(z(nz).le.78)  ! set top to 82km
---
>        do while(z(nz).le.80)  ! set top to 82km
398d395
< 	       print*,'nz,z= ',nz,z(1:nz)
636,640d632
< c	       if(cz(nz-j).lt.0.or.dsdh(id,j).lt.0) then
< c	        print*,'wrong cz or dsdh =',id,nz-j,cz(nz-j),dsdh(id,j)
< c     1		  ,sum
< c		stop
< c	       endif	
651,657c643
< 	 	 
< c         if(id.gt.0.and.scol(nz-id+1).gt.sum) then
< c	  print*,'scol in wrong sequence ',nz-id, sum, scol(nz-id+1),
< c     1	    id, dsdh(id,j-1), dsdh(id,j),cz(nz-j),cz(nz-j+1)
< c          print*,'nid=',nid
< c	  stop
< c	 endif   
---
> 
659,663c645
< ctyh      
<        DO id = nz-2, nz-1
<         if(scol(id).le.scol(id+1)) 
<      1   scol(id+1)=scol(id)-0.5*(scol(id-1)-scol(id))  ! correct the unknown error
<        enddo	
---
> 
18506c18488
<       REAL dtaer(kz,kw), omaer(kz,kw), gaer(kz,kw), daod(nzstem,naetp+1)
---
>       REAL dtaer(kz,kw), omaer(kz,kw), gaer(kz,kw), daod(nzstem,naetp)
18510c18492
<       real dtaerwl(kz,kw,naetp+1),omaerwl(kz,kw,naetp+1),absorp(kz,kw),
---
>       real dtaerwl(kz,kw,naetp+1),omaerwl(kz,kw,naetp+1), 
18962,18972c18944,18946
<        IF(ILAY.LE.nzstem-1) then
<         dtaerwl(ILAY,IW,IATP)=BEA(IW)*
<      1	((NDS(ILAY,IATP)+NDS(ILAY+1,IATP))/2.)*(Z(ILAY+1)-Z(ILAY))
<         if(dtaerwl(ILAY,IW,IATP).lt.0) then
< 	 print*,'negative dtaerwl in setaer ',dtaerwl(ILAY,IW,IATP),
<      1	  bea(IW),NDS(ILAY,IATP),NDS(ILAY+1,IATP),Z(ILAY+1),Z(ILAY),
<      2    iLAY, IATP, IW
<          print*,'nz,nztem,z=',nz,nzstem,z(1:nzstem)
<          stop
< 	endif
<        endif	 
---
>        IF(ILAY.LE.nzstem-1)  dtaerwl(ILAY,IW,IATP)=BEA(IW)*
>      >               ((NDS(ILAY,IATP)+NDS(ILAY+1,IATP))/2.)*
>      >                          (Z(ILAY+1)-Z(ILAY))
18986d18959
< 	absorp (ILAY,IW)=0.
18993,18994d18965
<          absorp (ILAY,IW)= absorp(ILAY,IW)+dtaerwl(ILAY,IW,IATP)        ! absorption
<      >                         * (1-omaerwl(ILAY,IW,IATP))
19003,19008c18974
< 	 if(dtaer(ilay,iw).lt.1e-25) then
< 	  omaer(ilay,iw)= 1.
< 	 else 
<  	  omaer(ilay,iw)= 1-absorp(ILAY,IW)/dtaer(ilay,iw)
< 	 endif 
< cccc	   omaer(ilay,iw)= omaerwl(ILAY,IW,NAETP+1)
---
> 	 omaer(ilay,iw)= omaerwl(ILAY,IW,NAETP+1)
19029,19030c18995,18996
<       DO ILAY=1, NZSTEM-1
<        DO IATP=1, NAETP+1      
---
>       DO IATP=1, NAETP+1
>        DO ILAY=1, NZSTEM-1
19033d18998
< c       DAOD(ILAY,NAETP+2)= omaer(ILAY,IW)   ! the 7th variable is Single Scaterring Albedo
19232,19246c19197,19199
<         deltaz = 1.E5 * (z(i+1)-z(i))
< 	if(abs(airstem(i+1)-airstem(i)).lt.1.and.i.le.(nzstem-2)) 
<      1	  then
< 	 airstem(i+1)=airstem(i)+(airstem(i+2)-airstem(i))/(
<      1	  z(i+2)-z(i))*(z(i+1)-z(i))        ! avoid zero difference of airstem
<         endif
<         cz(i) =  (airstem(i+1)-airstem(i))/         ! loading data from STEM
<      1	   ALOG(airstem(i+1)/airstem(i)) * deltaz
<          
< 	if(cz(i).lt.0.or.(.not.abs(cz(i)).lt.1e25)) then
< 	 print*,'cz wrong in setair ',i,cz(i),deltaz
< 	 print*,'z=',z(1:nzstem)
< 	 print*,'airstem=',airstem(1:nzstem)
< 	 stop
< 	endif 
---
>          deltaz = 1.E5 * (z(i+1)-z(i)) 
>          cz(i) =  (airstem(i+1)-airstem(i))/         ! loading data from STEM
>      1	   ALOG(airstem(i+1)/airstem(i)) * deltaz      
19280,19289c19233
<            dtrl(i,iw) = cz(i)*srayl(iw)
< c	   if(.not.(dtrl(i,iw).gt.-1.e36)) dtrl(i,iw)=0.   ! avoid overflow
< 	   
< 	   if(dtrl(i,iw).lt.0.or.(.not.abs(dtrl(i,iw)).lt.1e20)) then
< 	    print*,'wrong dtrl ',dtrl(i,iw),cz(i),srayl(iw)
< 	    print*,'airlev=',airlev(1:nz) 
< 	    print*,'z=',z(1:nzstem)
< 	    print*,'airstem=',airstem(1:nzstem)
< 	    stop
< 	   endif  
---
>             dtrl(i,iw) = cz(i)*srayl(iw)
19982,19992d19925
<            if(dttmp(iw).lt.0) then
< 	    print*,'negative dttmp in seto2 ',iz,iw,dttmp(iw),
<      1	      wlint(iw+1),wlint(iw)
< 	    print*,'dto2k=',igast,dto2k(1:nz,igast)
< 	    print*,'xso2int,delO2,dto2la=',xso2int(iw),delO2,
<      1	      dto2la(iz,1)
<             print*,'o2col=',o2col(1:nz)
< 	    print*,'scol=',scol(1:nz)
< 	    print*,'z=',z(1:nz)
< 	    stop
< 	   endif 
19998,20014c19931
<           dto2(iz,iw) = dtuser(iw)/(wl(iw+1)-wl(iw))
< c	    if(.not.(dto2(iz,iw).gt.-1.e36)) dto2(iz,iw)=0.   ! avoid overflow
< 	    
< 	  if(dto2(iz,iw).lt.0.or.(.not.abs(dto2(iz,iw)).lt.1e20)) then
< 	   print*,'wrong dto2 ',iz,iw,dto2(iz,iw)
< 	   print*,'dtuser=',dtuser(iw)
< 	   print*,'wl(iw+1),wl(iw)=',wl(iw+1),wl(iw)
< 	   stop
< 	  endif  
< 	    
< c	    if(dto2(iz,iw).lt.0) then
< c	     print*,'negative dto2 ',iz,iw,dto2(iz,iw),dtuser(iw),
< c     1	      wl(iw+1),wl(iw)
< c             print*,'cz=',cz(1:nz)
< c	     print*,'z=',z(1:nz)
< c	     stop
< c	    endif 
---
>             dto2(iz,iw) = dtuser(iw)/(wl(iw+1)-wl(iw))
