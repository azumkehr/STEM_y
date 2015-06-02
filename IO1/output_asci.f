C********************************************************************
      subroutine prtini
c********************************************************************
      parameter (mxspg=250)
      include 'aqcon1.cmm'
      include 'aqsymb.cmm'
      do i=1,numl(1,3)
      write(71,'(2x,i3,2x,a20)') i,sname(i,1)
      enddo
      return
      end
c**********************************************************************
      subroutine prtout(idate,ut,ix,iy,iz,sg1,sl1,iout)
c*********************************************************************
      parameter (mxspg=250)
      include 'aqcon1.cmm'
      include 'aqsymb.cmm'
      dimension sg1(ix,iy,iz,1),sl1(ix,iy,iz,1),iout(1),idate(3)
      write(70,*) 'idate=',idate,'  ut=',ut
      write(70,100) 
     1 ((((sg1(i,j,k,l),i=1,ix),j=1,iy),k=1,iz),l=1,numl(1,2))
100   format(2x,6e10.3)
      return
      end
c************************************************************
      subroutine close_ioap
c************************************************************
      return
      end
