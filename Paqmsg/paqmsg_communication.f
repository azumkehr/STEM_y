      module XYParallelCommunication

      use XYParallelDataMap
      use XYCommunicationLibrary
      !use XYParallelMemAlloc

c	versions of the communication library functions called
        integer :: dist_x_4D = 1, dist_y_4D = 1
	integer :: dist_x_3D = 1, dist_y_3D = 1
	integer :: dist_x_BDz = 1, dist_y_BDz = 1
	integer :: dist_x_2D = 1, dist_y_2D = 1
	integer :: dist_x_BD = 1, dist_y_BD = 1
	integer :: gath_x_4D = 1, gath_y_4D = 1
	integer :: shuf_x2y_4D = 1, shuf_y2x_4D = 1

      contains

c---------------------------------------------------------------------------------
      subroutine init_xy_versions()
      implicit none
      include 'mpif.h'
      
      integer :: ios,ierr,i
      character(20) :: version_name
      integer :: buf(14)

      if (myid .eq. 0) then
        OPEN (UNIT=11, FILE="XY_versions",STATUS="OLD",IOSTAT = ios)
        if (ios .ne. 0) then
          print*, 'Unable to open file XY_VERSIONS ',
     &		' default settings will be used'
	  buf(:) = 1
	else
  	  do i=1,14
            READ (11,*) version_name, buf(i)
	  end do
          CLOSE (11)
        end if

      end if
      
      call MPI_BCAST(buf, 14, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

      dist_x_4D = buf(1)
      dist_y_4D = buf(2)
      gath_x_4D = buf(3)
      gath_y_4D = buf(4)
      shuf_x2y_4D = buf(5)
      shuf_y2x_4D = buf(6)
      dist_x_2D = buf(7)
      dist_y_2D = buf(8)
      dist_x_BDz = buf(9)
      dist_y_BDz = buf(10)
      dist_x_3D = buf(11)
      dist_y_3D = buf(12)
      dist_x_BD = buf(13)
      dist_y_BD = buf(14)

      end subroutine init_xy_versions
      
      subroutine init_xy(ix,iy,iz,ixloc,iyloc,Ns)

      use XYCommDataTypes
      implicit none
      integer :: ix,iy,iz,ixloc,iyloc,Ns
      integer :: i, ierr
      
      call CreateMap(ix, iy, ixloc, iyloc)
      call init_xy_versions()
      call CreateCommDataTypes(ix,iy,iz,ixloc,iyloc,Ns)
      
      end subroutine init_xy
      
c---------------------------------------------------------------------------------
c     Call the right version of the different routines
c---------------------------------------------------------------------------------

      subroutine distrib_x_4D(ix,iy,iz,iyloc,Ns,s,s_x)
      implicit none

      integer :: ix, iy, iz, iyloc, Ns
      real :: s(ix,iy,iz,Ns)
      real :: s_x(ix,iyloc,iz,Ns)
      
      SELECT CASE (dist_x_4D)
      CASE(1) 
        call distrib_x_4D_v1(ix,iy,iz,iyloc,Ns,s,s_x)
      CASE(2) 
        call distrib_x_4D_v2(ix,iy,iz,iyloc,Ns,s,s_x)
      CASE(3) 
        call distrib_x_4D_v3(ix,iy,iz,iyloc,Ns,s,s_x)
      CASE DEFAULT 
        call distrib_x_4D_v1(ix,iy,iz,iyloc,Ns,s,s_x)
      END SELECT
      
      end subroutine distrib_x_4D
      

      subroutine distrib_y_4D(ix,iy,iz,ixloc,Ns,s,s_y)
      implicit none

      integer :: ix, iy, iz, ixloc, Ns
      real :: s(ix,iy,iz,Ns)
      real :: s_y(ixloc,iy,iz,Ns)
      
      SELECT CASE (dist_y_4D)
      CASE(1) 
        call distrib_y_4D_v1(ix,iy,iz,ixloc,Ns,s,s_y)
      CASE(2) 
        call distrib_y_4D_v2(ix,iy,iz,ixloc,Ns,s,s_y)
      CASE(3) 
        call distrib_y_4D_v3(ix,iy,iz,ixloc,Ns,s,s_y)
      CASE DEFAULT
        call distrib_y_4D_v1(ix,iy,iz,ixloc,Ns,s,s_y)
      END SELECT
      
      end subroutine distrib_y_4D


      subroutine distrib_x_3D(ix,iy,iz,iyloc,s,s_x)
      implicit none

      integer :: ix, iy, iz, iyloc
      real :: s(ix,iy,iz)
      real :: s_x(ix,iyloc,iz)
      
      SELECT CASE (dist_x_3D)
      CASE(1) 
        call distrib_x_3D_v1(ix,iy,iz,iyloc,s,s_x)
      CASE(2) 
        call distrib_x_3D_v2(ix,iy,iz,iyloc,s,s_x)
      CASE(3) 
        call distrib_x_3D_v3(ix,iy,iz,iyloc,s,s_x)
      CASE DEFAULT 
        call distrib_x_3D_v1(ix,iy,iz,iyloc,s,s_x)
      END SELECT
      
      end subroutine distrib_x_3D
      

      subroutine distrib_y_3D(ix,iy,iz,ixloc,s,s_y)
      implicit none

      integer :: ix, iy, iz, ixloc
      real :: s(ix,iy,iz)
      real :: s_y(ixloc,iy,iz)
      
      SELECT CASE (dist_y_3D)
      CASE(1) 
        call distrib_y_3D_v1(ix,iy,iz,ixloc,s,s_y)
      CASE(2) 
        call distrib_y_3D_v2(ix,iy,iz,ixloc,s,s_y)
      CASE(3) 
        call distrib_y_3D_v3(ix,iy,iz,ixloc,s,s_y)
      CASE DEFAULT
        call distrib_y_3D_v1(ix,iy,iz,ixloc,s,s_y)
      END SELECT
      
      end subroutine distrib_y_3D


      subroutine distrib_x_BDz(ix,iy,iyloc,Ns,s,s_x)
      implicit none

      integer :: ix, iy, iyloc, Ns
      real :: s(ix,iy,Ns)
      real :: s_x(ix,iyloc,Ns)
      
      SELECT CASE (dist_x_BDz)
      CASE(1) 
        call distrib_x_BDz_v1(ix,iy,iyloc,Ns,s,s_x)
      CASE(2) 
        call distrib_x_BDz_v2(ix,iy,iyloc,Ns,s,s_x)
      CASE(3) 
        call distrib_x_BDz_v3(ix,iy,iyloc,Ns,s,s_x)
      CASE DEFAULT 
        call distrib_x_BDz_v1(ix,iy,iyloc,Ns,s,s_x)
      END SELECT
      
      end subroutine distrib_x_BDz
      

      subroutine distrib_y_BDz(ix,iy,ixloc,Ns,s,s_y)
      implicit none

      integer :: ix, iy, ixloc, Ns
      real :: s(ix,iy,Ns)
      real :: s_y(ixloc,iy,Ns)
      
      SELECT CASE (dist_y_BDz)
      CASE(1) 
        call distrib_y_BDz_v1(ix,iy,ixloc,Ns,s,s_y)
      CASE(2) 
        call distrib_y_BDz_v2(ix,iy,ixloc,Ns,s,s_y)
      CASE(3) 
        call distrib_y_BDz_v3(ix,iy,ixloc,Ns,s,s_y)
      CASE DEFAULT
        call distrib_y_BDz_v1(ix,iy,ixloc,Ns,s,s_y)
      END SELECT
      
      end subroutine distrib_y_BDz


      subroutine distrib_x_2D(ix,iy,iyloc,s,s_x)
      implicit none

      integer :: ix, iy, iyloc
      real :: s(ix,iy)
      real :: s_x(ix,iyloc)
      
      SELECT CASE (dist_x_2D)
      CASE(1) 
        call distrib_x_2D_v1(ix,iy,iyloc,s,s_x)
      CASE(2) 
        call distrib_x_2D_v2(ix,iy,iyloc,s,s_x)
      CASE(3) 
        call distrib_x_2D_v3(ix,iy,iyloc,s,s_x)
      CASE DEFAULT 
        call distrib_x_2D_v1(ix,iy,iyloc,s,s_x)
      END SELECT
      
      end subroutine distrib_x_2D
      

      subroutine distrib_y_2D(ix,iy,ixloc,s,s_y)
      implicit none

      integer :: ix, iy, ixloc
      real :: s(ix,iy)
      real :: s_y(ixloc,iy)
      
      SELECT CASE (dist_y_2D)
      CASE(1) 
        call distrib_y_2D_v1(ix,iy,ixloc,s,s_y)
      CASE(2) 
        call distrib_y_2D_v2(ix,iy,ixloc,s,s_y)
      CASE(3) 
        call distrib_y_2D_v3(ix,iy,ixloc,s,s_y)
      CASE DEFAULT
        call distrib_y_2D_v1(ix,iy,ixloc,s,s_y)
      END SELECT
      
      end subroutine distrib_y_2D


      subroutine distrib_x_BD(iy,iz,iyloc,Ns,s,s_x)
      implicit none

      integer :: iy, iz, iyloc, Ns
      real :: s(iy,iz,2,Ns)
      real :: s_x(iyloc,iz,2,Ns)
      
      SELECT CASE (dist_x_BD)
      CASE(1) 
        call distrib_x_BD_v1(iy,iz,iyloc,Ns,s,s_x)
      CASE(2) 
        call distrib_x_BD_v2(iy,iz,iyloc,Ns,s,s_x)
      CASE(3) 
        call distrib_x_BD_v3(iy,iz,iyloc,Ns,s,s_x)
      CASE DEFAULT 
        call distrib_x_BD_v1(iy,iz,iyloc,Ns,s,s_x)
      END SELECT
      
      end subroutine distrib_x_BD
      

      subroutine distrib_y_BD(ix,iz,ixloc,Ns,s,s_y)
      implicit none

      integer :: ix, iz, ixloc, Ns
      real :: s(ix,iz,2,Ns)
      real :: s_y(ixloc,iz,2,Ns)
      
      SELECT CASE (dist_y_BD)
      CASE(1) 
        call distrib_y_BD_v1(ix,iz,ixloc,Ns,s,s_y)
      CASE(2) 
        call distrib_y_BD_v2(ix,iz,ixloc,Ns,s,s_y)
      CASE(3) 
        call distrib_y_BD_v3(ix,iz,ixloc,Ns,s,s_y)
      CASE DEFAULT
        call distrib_y_BD_v1(ix,iz,ixloc,Ns,s,s_y)
      END SELECT
      
      end subroutine distrib_y_BD


      subroutine gather_x_4D(ix,iy,iz,iyloc,Ns,s,s_x)
      implicit none

      integer :: ix, iy, iz, iyloc, Ns
      real :: s(ix,iy,iz,Ns)
      real :: s_x(ix,iyloc,iz,Ns)
      
      SELECT CASE (gath_x_4D)
      CASE(1) 
        call gather_x_4D_v1(ix,iy,iz,iyloc,Ns,s,s_x)
      CASE(2) 
        call gather_x_4D_v2(ix,iy,iz,iyloc,Ns,s,s_x)
      CASE(3) 
        call gather_x_4D_v3(ix,iy,iz,iyloc,Ns,s,s_x)
      CASE DEFAULT 
        call gather_x_4D_v1(ix,iy,iz,iyloc,Ns,s,s_x)
      END SELECT
      
      end subroutine gather_x_4D
      

      subroutine gather_y_4D(ix,iy,iz,ixloc,Ns,s,s_y)
      implicit none

      integer :: ix, iy, iz, ixloc, Ns
      real :: s(ix,iy,iz,Ns)
      real :: s_y(ixloc,iy,iz,Ns)
      
      SELECT CASE (gath_y_4D)
      CASE(1) 
        call gather_y_4D_v1(ix,iy,iz,ixloc,Ns,s,s_y)
      CASE(2) 
        call gather_y_4D_v2(ix,iy,iz,ixloc,Ns,s,s_y)
      CASE(3) 
        call gather_y_4D_v3(ix,iy,iz,ixloc,Ns,s,s_y)
      CASE DEFAULT
        call gather_y_4D_v1(ix,iy,iz,ixloc,Ns,s,s_y)
      END SELECT
      
      end subroutine gather_y_4D


      subroutine shuffle_x2y_4D(ix,iy,iz,ixloc,iyloc,Ns,s_x,s_y)
      
      implicit none
      integer :: ix, iy, iz, ixloc, iyloc, Ns
      real :: s_x(ix, iyloc, iz, Ns)
      real :: s_y(ixloc, iy, iz, Ns)
      
      SELECT CASE (shuf_x2y_4D)
      CASE(1)
        call shuffle_x2y_4D_v1(ix,iy,iz,ixloc,iyloc,Ns,s_x,s_y)
      CASE(2)
        call shuffle_x2y_4D_v2(ix,iy,iz,ixloc,iyloc,Ns,s_x,s_y)
      CASE(3)
        call shuffle_x2y_4D_v3(ix,iy,iz,ixloc,iyloc,Ns,s_x,s_y)
      CASE DEFAULT
        call shuffle_x2y_4D_v1(ix,iy,iz,ixloc,iyloc,Ns,s_x,s_y)
      END SELECT

      end subroutine shuffle_x2y_4D


      subroutine shuffle_y2x_4D(ix,iy,iz,ixloc,iyloc,Ns,s_x,s_y)
      
      implicit none
      integer :: ix, iy, iz, ixloc, iyloc, Ns
      real :: s_x(ix, iyloc, iz, Ns)
      real :: s_y(ixloc, iy, iz, Ns)
      
      SELECT CASE (shuf_y2x_4D)
      CASE(1)
        call shuffle_y2x_4D_v1(ix,iy,iz,ixloc,iyloc,Ns,s_x,s_y)
      CASE(2)
        call shuffle_y2x_4D_v2(ix,iy,iz,ixloc,iyloc,Ns,s_x,s_y)
      CASE(3)
        call shuffle_y2x_4D_v3(ix,iy,iz,ixloc,iyloc,Ns,s_x,s_y)
      CASE DEFAULT
        call shuffle_y2x_4D_v1(ix,iy,iz,ixloc,iyloc,Ns,s_x,s_y)
      END SELECT
      
      end subroutine shuffle_y2x_4D

      
      end module XYParallelCommunication

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      module HVParallelCommunication

      use HVParallelDataMap
      use HVCommunicationLibrary
      !use HVParallelMemAlloc

c	versions of the communication library functions called
        integer :: dist_h_4D = 1, dist_v_4D = 1
	integer :: dist_h_3D = 1, dist_v_3D = 1
	integer :: dist_h_BDz = 1, dist_v_BDz = 1
	integer :: dist_h_2D = 1, dist_v_2D = 1
	integer :: dist_xh_BD = 1, dist_yh_BD = 1
	integer :: gath_h_4D = 1, gath_v_4D = 1
	integer :: shuf_h2v_4D = 1, shuf_v2h_4D = 1
      
      contains


c---------------------------------------------------------------------------------
      subroutine init_hv_versions()
      implicit none
      include 'mpif.h'
      
      integer :: ios,ierr,i
      character(20) :: version_name
      integer :: buf(13)

      if (MyId .eq. 0) then
        OPEN (UNIT=11, FILE="HV_versions",STATUS="OLD",IOSTAT = ios)
        if (ios .ne. 0) then
          print*, 'Unable to open file HV_VERSIONS ',
     &		' default settings will be used'
	  buf(:) = 1	
        else
	  do i=1,13
            READ (11,*) version_name, buf(i)
	  end do
          CLOSE (11)
	end if


      end if
      
      call MPI_BCAST(buf, 13, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

      dist_h_4D = buf(1)
      dist_v_4D = buf(2)
      gath_h_4D = buf(3)
      gath_v_4D = buf(4)
      shuf_h2v_4D = buf(5)
      shuf_v2h_4D = buf(6)
      dist_v_2D = buf(7)
      dist_h_BDz = buf(8)
      dist_v_BDz = buf(9)
      dist_h_3D = buf(10)
      dist_v_3D = buf(11)
      dist_xh_BD = buf(12)
      dist_yh_BD = buf(13)

      end subroutine init_hv_versions

      subroutine init_hv(ix,iy,iz,izloc,icloc,Ns)

      use HVCommDataTypes
      implicit none
      integer :: ix,iy,iz,izloc,icloc,Ns
      integer :: i, ierr
      
      call CreateMap(ix, iy, iz, izloc, icloc)
      call init_hv_versions()
      call CreateCommDataTypes(ix,iy,iz,izloc,icloc,Ns)
      
      end subroutine init_hv

c---------------------------------------------------------------------------------
c     Call the right version of the different routines
c---------------------------------------------------------------------------------

      subroutine distrib_h_4D(ix,iy,iz,izloc,Ns,s,s_h)
      implicit none

      integer :: ix, iy, iz, izloc, Ns
      real :: s(ix,iy,iz,Ns)
      real :: s_h(ix,iy,izloc,Ns)
      
      SELECT CASE (dist_h_4D)
      CASE(1) 
        call distrib_h_4D_v1(ix,iy,iz,izloc,Ns,s,s_h)
      CASE(2) 
        call distrib_h_4D_v2(ix,iy,iz,izloc,Ns,s,s_h)
      CASE(3) 
        call distrib_h_4D_v3(ix,iy,iz,izloc,Ns,s,s_h)
      CASE DEFAULT 
        call distrib_h_4D_v1(ix,iy,iz,izloc,Ns,s,s_h)
      END SELECT
      
      end subroutine distrib_h_4D
      

      subroutine distrib_v_4D(ix,iy,iz,icloc,Ns,s,s_v)
      implicit none

      integer :: ix, iy, iz, icloc, Ns
      real :: s(ix,iy,iz,Ns)
      real :: s_v(1,icloc,iz,Ns)
      
      SELECT CASE (dist_v_4D)
      CASE(1) 
        call distrib_v_4D_v1(ix,iy,iz,icloc,Ns,s,s_v)
      CASE(2) 
        call distrib_v_4D_v2(ix,iy,iz,icloc,Ns,s,s_v)
      CASE(3) 
        call distrib_v_4D_v3(ix,iy,iz,icloc,Ns,s,s_v)
      CASE DEFAULT
        call distrib_v_4D_v1(ix,iy,iz,icloc,Ns,s,s_v)
      END SELECT
      
      end subroutine distrib_v_4D


      subroutine distrib_h_3D(ix,iy,iz,izloc,s,s_h)
      implicit none

      integer :: ix, iy, iz, izloc
      real :: s(ix,iy,iz)
      real :: s_h(ix,iy,izloc)
      
      SELECT CASE (dist_h_3D)
      CASE(1) 
        call distrib_h_3D_v1(ix,iy,iz,izloc,s,s_h)
      CASE(2) 
        call distrib_h_3D_v2(ix,iy,iz,izloc,s,s_h)
      CASE(3) 
        call distrib_h_3D_v3(ix,iy,iz,izloc,s,s_h)
      CASE DEFAULT 
        call distrib_h_3D_v1(ix,iy,iz,izloc,s,s_h)
      END SELECT
      
      end subroutine distrib_h_3D
      

      subroutine distrib_v_3D(ix,iy,iz,icloc,s,s_v)
      implicit none

      integer :: ix, iy, iz, icloc
      real :: s(ix,iy,iz)
      real :: s_v(1,icloc,iz)
      
      SELECT CASE (dist_v_3D)
      CASE(1) 
        call distrib_v_3D_v1(ix,iy,iz,icloc,s,s_v)
      CASE(2) 
        call distrib_v_3D_v2(ix,iy,iz,icloc,s,s_v)
      CASE(3) 
        call distrib_v_3D_v3(ix,iy,iz,icloc,s,s_v)
      CASE DEFAULT
        call distrib_v_3D_v1(ix,iy,iz,icloc,s,s_v)
      END SELECT
      
      end subroutine distrib_v_3D


      subroutine distrib_v_BDz(ix,iy,icloc,Ns,s,s_v)
      implicit none

      integer :: ix, iy, icloc, Ns
      real :: s(ix,iy,Ns)
      real :: s_v(1,icloc,Ns)
      
      SELECT CASE (dist_v_BDz)
      CASE(1) 
        call distrib_v_BDz_v1(ix,iy,icloc,Ns,s,s_v)
      CASE(2) 
        call distrib_v_BDz_v2(ix,iy,icloc,Ns,s,s_v)
      CASE(3) 
        call distrib_v_BDz_v3(ix,iy,icloc,Ns,s,s_v)
      CASE DEFAULT
        call distrib_v_BDz_v1(ix,iy,icloc,Ns,s,s_v)
      END SELECT
      
      end subroutine distrib_v_BDz


      subroutine distrib_h_2D(ix,iy,s)
      implicit none

      integer :: ix, iy
      real :: s(ix,iy)
      
        call distrib_h_2D_v1(ix,iy,s)
      
      end subroutine distrib_h_2D
      

      subroutine distrib_v_2D(ix,iy,icloc,s,s_v)
      implicit none

      integer :: ix, iy, icloc
      real :: s(ix,iy)
      real :: s_v(1,icloc)
      
      SELECT CASE (dist_v_2D)
      CASE(1) 
        call distrib_v_2D_v1(ix,iy,icloc,s,s_v)
      CASE(2) 
        call distrib_v_2D_v2(ix,iy,icloc,s,s_v)
      CASE(3) 
        call distrib_v_2D_v3(ix,iy,icloc,s,s_v)
      CASE DEFAULT
        call distrib_v_2D_v1(ix,iy,icloc,s,s_v)
      END SELECT
      
      end subroutine distrib_v_2D


      subroutine distrib_xh_BD(ix,iz,izloc,Ns,s,s_h)
      implicit none

      integer :: ix, iz, izloc, Ns
      real :: s(ix,iz,2,Ns)
      real :: s_h(ix,izloc,2,Ns)
      
      SELECT CASE (dist_xh_BD)
      CASE(1) 
        call distrib_xh_BD_v1(ix,iz,izloc,Ns,s,s_h)
      CASE(2) 
        call distrib_xh_BD_v2(ix,iz,izloc,Ns,s,s_h)
      CASE(3) 
        call distrib_xh_BD_v3(ix,iz,izloc,Ns,s,s_h)
      CASE DEFAULT
        call distrib_xh_BD_v1(ix,iz,izloc,Ns,s,s_h)
      END SELECT
      
      end subroutine distrib_xh_BD


      subroutine distrib_yh_BD(iy,iz,izloc,Ns,s,s_h)
      implicit none

      integer :: iy, iz, izloc, Ns
      real :: s(iy,iz,2,Ns)
      real :: s_h(iy,izloc,2,Ns)
      
      SELECT CASE (dist_yh_BD)
      CASE(1) 
        call distrib_yh_BD_v1(iy,iz,izloc,Ns,s,s_h)
      CASE(2) 
        call distrib_yh_BD_v2(iy,iz,izloc,Ns,s,s_h)
      CASE(3) 
        call distrib_yh_BD_v3(iy,iz,izloc,Ns,s,s_h)
      CASE DEFAULT 
        call distrib_yh_BD_v1(iy,iz,izloc,Ns,s,s_h)
      END SELECT
      
      end subroutine distrib_yh_BD


      subroutine gather_h_4D(ix,iy,iz,izloc,Ns,s,s_h)
      implicit none

      integer :: ix, iy, iz, izloc, Ns
      real :: s(ix,iy,iz,Ns)
      real :: s_h(ix,iy,izloc,Ns)
      
      SELECT CASE (gath_h_4D)
      CASE(1) 
        call gather_h_4D_v1(ix,iy,iz,izloc,Ns,s,s_h)
      CASE(2) 
        call gather_h_4D_v2(ix,iy,iz,izloc,Ns,s,s_h)
      CASE(3) 
        call gather_h_4D_v3(ix,iy,iz,izloc,Ns,s,s_h)
      CASE DEFAULT 
        call gather_h_4D_v1(ix,iy,iz,izloc,Ns,s,s_h)
      END SELECT
      
      end subroutine gather_h_4D
      

      subroutine gather_v_4D(ix,iy,iz,icloc,Ns,s,s_v)
      implicit none

      integer :: ix, iy, iz, icloc, Ns
      real :: s(ix,iy,iz,Ns)
      real :: s_v(1,icloc,iz,Ns)
      
      SELECT CASE (gath_v_4D)
      CASE(1) 
        call gather_v_4D_v1(ix,iy,iz,icloc,Ns,s,s_v)
      CASE(2) 
        call gather_v_4D_v2(ix,iy,iz,icloc,Ns,s,s_v)
      CASE(3) 
        call gather_v_4D_v3(ix,iy,iz,icloc,Ns,s,s_v)
      CASE DEFAULT
        call gather_v_4D_v1(ix,iy,iz,icloc,Ns,s,s_v)
      END SELECT
      
      end subroutine gather_v_4D


      subroutine shuffle_h2v_4D(ix,iy,iz,izloc,icloc,Ns,s_h,s_v)
      
      implicit none
      integer :: ix, iy, iz, izloc, icloc, Ns
      real :: s_h(ix, iy, izloc, Ns)
      real :: s_v(1,icloc, iz, Ns)
      
      SELECT CASE (shuf_h2v_4D)
      CASE(1)
        call shuffle_h2v_4D_v1(ix,iy,iz,izloc,icloc,Ns,s_h,s_v)
      CASE(2)
        call shuffle_h2v_4D_v2(ix,iy,iz,izloc,icloc,Ns,s_h,s_v)
      CASE(3)
        call shuffle_h2v_4D_v3(ix,iy,iz,izloc,icloc,Ns,s_h,s_v)
      CASE DEFAULT
        call shuffle_h2v_4D_v1(ix,iy,iz,izloc,icloc,Ns,s_h,s_v)
      END SELECT

      end subroutine shuffle_h2v_4D


      subroutine shuffle_v2h_4D(ix,iy,iz,izloc,icloc,Ns,s_h,s_v)
      
      implicit none
      integer :: ix, iy, iz, izloc, icloc, Ns
      real :: s_h(ix, iy, izloc, Ns)
      real :: s_v(1, icloc, iz, Ns)
      
      SELECT CASE (shuf_v2h_4D)
      CASE(1)
        call shuffle_v2h_4D_v1(ix,iy,iz,izloc,icloc,Ns,s_h,s_v)
      CASE(2)
        call shuffle_v2h_4D_v2(ix,iy,iz,izloc,icloc,Ns,s_h,s_v)
      CASE(3)
        call shuffle_v2h_4D_v3(ix,iy,iz,izloc,icloc,Ns,s_h,s_v)
      CASE DEFAULT
        call shuffle_v2h_4D_v1(ix,iy,iz,izloc,icloc,Ns,s_h,s_v)
      END SELECT
      
      end subroutine shuffle_v2h_4D

c
      end module HVParallelCommunication
