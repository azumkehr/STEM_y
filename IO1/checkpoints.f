! Name a checkpoint file
      subroutine fname_chkp2(MyId,fname)
      integer :: MyId
      character(len=16) :: fname
      character(len=4) :: buffer
      
      write(buffer,"(I4)") 1000+MyId
      
      fname = "Chkp2_"//buffer(2:4)//".chk"
      
      end subroutine fname_chkp2


! Open a checkpoint file
      subroutine open_chkp2_conc(MyId,iunit,fname,len_array)
      implicit none
      integer :: MyId, iunit, len_record, len_array, ios
      character(len=32) :: fname
      character(len=4) :: buffer
      character(len=48) :: filename
      
      write(buffer,"(I4)") 1000+MyId
      
      
      filename=trim(fname)//'_'//buffer(2:4)//'.chk'
      
      len_record = len_array*4
           
      ! print*,MyId,' Open Checkpoint file ',trim(fname),'&',filename
      open( UNIT=iunit, FILE=trim(filename), 
     &       ACCESS="DIRECT", 
     &       RECL=len_record, IOSTAT=ios )
      if (ios/=0) then
         print*,'Error opening Conc Checkpoint file ',filename
         print*,'      (iostat = ',ios,')'
	 stop
      end if
      
      end subroutine open_chkp2_conc

! Delete a checkpoint file
      subroutine rm_chkp2_conc(MyId,iunit,fname)
      implicit none
      integer :: MyId, iunit, len_record, ios
      character(len=32) :: fname
      character(len=4) :: buffer
      character(len=48) :: filename

      write(buffer,"(I4)") 1000+MyId


      filename=trim(fname)//'_'//buffer(2:4)//'.chk'

      close( iunit )

      call system( 'rm '//trim(filename) )

      end subroutine rm_chkp2_conc

      
      
! Write all concentrations in a checkpoint file, in the record numrec
      subroutine write_4D_chkp2(Nx,Ny,Nz,Ns,Conc,numrec,iunit)
      implicit none
      integer :: Nx, Ny, Nz, Ns
      integer :: numrec, iunit, ios
      real :: Conc(Nx,Ny,Nz,Ns)
      
      write(UNIT=iunit, REC=numrec, IOSTAT=ios) 
     &             Conc(1:Nx,1:Ny,1:Nz,1:Ns)
      if (ios/=0) then
         print*,'Error during write_4D_chkp2 unit ',iunit
	 print*,'      (iostat = ',ios,')'
         stop 'write_4D_chkp2 failed'
      end if
      
      end subroutine write_4D_chkp2
       
      
! read all concentrations in a checkpoint file, in the record numrec
      subroutine read_4D_chkp2(Nx,Ny,Nz,Ns,Conc,numrec,iunit)
      implicit none
      integer :: Nx, Ny, Nz, Ns
      integer :: numrec, iunit, ios
      real :: Conc(Nx,Ny,Nz,Ns)
      
      read(UNIT=iunit, REC=numrec, IOSTAT=ios) 
     &              Conc(1:Nx,1:Ny,1:Nz,1:Ns)
      if (ios/=0) then
         print*,'Error during read_4D_chkp2 unit ',iunit
	 print*,'      (iostat = ',ios,')'
         stop
      end if
           	 
      end subroutine read_4D_chkp2 
      
!      Write all concentrations in a checkpoint file, in the record numrec
      subroutine write_meteo_chkp2(Nx,Ny,Nz,Ns,u,v,w,kh,kv,
     &          vg,q,em,sx,sy,sz,wc,wr,sprc,rvel,cldod,kctop,
     &          ccover,dobson,numrec,iunit)
      implicit none
      integer :: Nx, Ny, Nz, Ns
      integer :: numrec, iunit, ios
      real ::     u(1:Nx,1:Ny,1:Nz),
     &             v(1:Nx,1:Ny,1:Nz),
     &             w(1:Nx,1:Ny,1:Nz),
     &             kh(1:Nx,1:Ny,1:Nz),
     &             kv(1:Nx,1:Ny,1:Nz),
     &             vg(1:Nx,1:Ny,1:Ns),
     &             q(1:Nx,1:Ny,1:Ns),
     &             em(1:Nx,1:Ny,1:Nz,1:Ns),
     &             sx(1:Ny,1:Nz,1:2,1:Ns),
     &             sy(1:Nx,1:Nz,1:2,1:Ns),
     &             sz(1:Nx,1:Ny,1:Ns),
     &             wc(1:Nx,1:Ny,1:Nz),
     &             wr(1:Nx,1:Ny,1:Nz),
     &             sprc(1:Nx,1:Ny),
     &             rvel(1:Nx,1:Ny,1:Nz),
     &             cldod(1:Nx,1:Ny,1:Nz),
     &             kctop(1:Nx,1:Ny),
     &             ccover(1:Nx,1:Ny,1:Nz),
     &             dobson(1:Nx,1:Ny)   
Conc(Nx,Ny,Nz,Ns)
      
      write(UNIT=iunit, REC=numrec, IOSTAT=ios) 
     &             u(1:Nx,1:Ny,1:Nz),
     &             v(1:Nx,1:Ny,1:Nz),
     &             w(1:Nx,1:Ny,1:Nz),
     &             kh(1:Nx,1:Ny,1:Nz),
     &             kv(1:Nx,1:Ny,1:Nz),
     &             vg(1:Nx,1:Ny,1:Ns),
     &             q(1:Nx,1:Ny,1:Ns),
     &             em(1:Nx,1:Ny,1:Nz,1:Ns),
     &             sx(1:Ny,1:Nz,1:2,1:Ns),
     &             sy(1:Nx,1:Nz,1:2,1:Ns),
     &             sz(1:Nx,1:Ny,1:Ns),
     &             wc(1:Nx,1:Ny,1:Nz),
     &             wr(1:Nx,1:Ny,1:Nz),
     &             sprc(1:Nx,1:Ny),
     &             rvel(1:Nx,1:Ny,1:Nz),
     &             cldod(1:Nx,1:Ny,1:Nz),
     &             kctop(1:Nx,1:Ny),
     &             ccover(1:Nx,1:Ny,1:Nz),
     &             dobson(1:Nx,1:Ny)   
     
      if (ios/=0) then
         print*,'Error during write_meteo_chkp2 ',iunit
	 stop
      end if
      
      end subroutine write_meteo_chkp2
       
      
! read all concentrations in a checkpoint file, in the record numrec
      subroutine read_meteo_chkp2(Nx,Ny,Nz,Ns,u,v,w,kh,kv,
     &          vg,q,em,sx,sy,sz,wc,wr,sprc,rvel,cldod,kctop,
     &          ccover,dobson,numrec,iunit)
      implicit none
      integer :: Nx, Ny, Nz, Ns
      integer :: numrec, iunit, ios
      real ::     u(1:Nx,1:Ny,1:Nz),
     &             v(1:Nx,1:Ny,1:Nz),
     &             w(1:Nx,1:Ny,1:Nz),
     &             kh(1:Nx,1:Ny,1:Nz),
     &             kv(1:Nx,1:Ny,1:Nz),
     &             vg(1:Nx,1:Ny,1:Ns),
     &             q(1:Nx,1:Ny,1:Ns),
     &             em(1:Nx,1:Ny,1:Nz,1:Ns),
     &             sx(1:Ny,1:Nz,1:2,1:Ns),
     &             sy(1:Nx,1:Nz,1:2,1:Ns),
     &             sz(1:Nx,1:Ny,1:Ns),
     &             wc(1:Nx,1:Ny,1:Nz),
     &             wr(1:Nx,1:Ny,1:Nz),
     &             sprc(1:Nx,1:Ny),
     &             rvel(1:Nx,1:Ny,1:Nz),
     &             cldod(1:Nx,1:Ny,1:Nz),
     &             kctop(1:Nx,1:Ny),
     &             ccover(1:Nx,1:Ny,1:Nz),
     &             dobson(1:Nx,1:Ny)   
Conc(Nx,Ny,Nz,Ns)
      
      read(UNIT=iunit, REC=numrec, IOSTAT=ios) 
     &             u(1:Nx,1:Ny,1:Nz),
     &             v(1:Nx,1:Ny,1:Nz),
     &             w(1:Nx,1:Ny,1:Nz),
     &             kh(1:Nx,1:Ny,1:Nz),
     &             kv(1:Nx,1:Ny,1:Nz),
     &             vg(1:Nx,1:Ny,1:Ns),
     &             q(1:Nx,1:Ny,1:Ns),
     &             em(1:Nx,1:Ny,1:Nz,1:Ns),
     &             sx(1:Ny,1:Nz,1:2,1:Ns),
     &             sy(1:Nx,1:Nz,1:2,1:Ns),
     &             sz(1:Nx,1:Ny,1:Ns),
     &             wc(1:Nx,1:Ny,1:Nz),
     &             wr(1:Nx,1:Ny,1:Nz),
     &             sprc(1:Nx,1:Ny),
     &             rvel(1:Nx,1:Ny,1:Nz),
     &             cldod(1:Nx,1:Ny,1:Nz),
     &             kctop(1:Nx,1:Ny),
     &             ccover(1:Nx,1:Ny,1:Nz),
     &             dobson(1:Nx,1:Ny)   
     
      if (ios/=0) then
         print*,'Error during read_meteo_chkp2 ',iunit
	 stop
      end if
           	 
      end subroutine read_meteo_chkp2
