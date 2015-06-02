c---------------------------------------------------------------------------------------
c
c Three versions of the following functions are defined in this file:
c
c	distrib_x_4D  	Distribution of the concentration array in
c			x-slices from the master to the workers
c	distrib_y_4D 	Distribution of the concentration array in
c			y-slices from the master to the workers
c	gather_x_4D	Gathering of the concentration array from
c			x-slices from the workers at the master
c	gather_y_4D	Gathering of the concentration array from
c			y-slices from the workers at the master
c	shuffle_x2y_4D	Shuffling of the concentration array
c			from x to y-slices at each worker
c	shuffle_y2x_4D	Shuffling of the concentration array
c			from y to x-slices at each worker
c	distrib_x_2D	Distribution of the 2D (x,y) arrays in
c			x-slices from the master to the workers
c	distrib_y_2D	Distribution of the 2D (x,y) arrays in
c			y-slices from the master to the workers
c	distrib_x_BDz	Distribution of the BDz (x,y,Ns) arrays in
c			x-slices from the master to the workers
c	distrib_y_BDz	Distribution of the BDz (x,y,Ns) arrays in
c			y-slices from the master to the workers
c	distrib_x_3D	Distribution of the 3D (x,y,z) arrays in
c			x-slices from the master to the workers
c	distrib_y_3D	Distribution of the 3D (x,y,z) arrays in
c			y-slices from the master to the workers
c	distrib_x_BDx	Distribution of the BDx (y,z,Ns) boundary data 
c			arrays in x-slices from the master to the workers
c	distrib_y_BDy	Distribution of the BDy (x,z,Ns) boundary data 
c			arrays in y-slices from the master to the workers
c	distrib_x_BD	Distribution of the boundary data arrays in
c			x-slices from the master to the workers - Stem-III 		c			dependent: BDx2 (y,z,2,Ns)
c	distrib_y_BD	Distribution of the boundary data arrays in
c			y-slices from the master to the workers - Stem-III 		c			dependent: BDy2 (x,z,2,Ns)
c
c---------------------------------------------------------------------------------------

      module XYCommunicationLibrary

      use XYParallelDataMap
      use XYCommDataTypes
      include 'mpif.h'

      contains
      
c-----------------------------------------------------------------
c  Master distributes a 4D array of data in x-slice format to all xworkers
c  each x-slice seperately, standard version
c-----------------------------------------------------------------
      subroutine distrib_x_4D_v1(ix,iy,iz,iyloc,Ns,s,s_x)
c      
      implicit none
c
      integer :: ix, iy, iz, iyloc, Ns
      real :: s(ix,iy,iz,Ns)
      real :: s_x(ix,iyloc,iz,Ns)
      integer :: i, ierr 
      integer :: status(MPI_STATUS_SIZE, iyloc)
      integer :: request(iyloc)
c      
      if (Master) then
        do i=1,iy
          call MPI_SEND(s(1,i,1,1), 1, GLOBAL_4D_XSLICE, 
     & 		owner_of_xslice(i), i,  MPI_COMM_WORLD, ierr)
	end do
      else if (XWorker) then
        do i=1, no_of_xslices(myid)
          call MPI_IRECV(s_x(1,i,1,1), 1, LOCAL_4D_XSLICE, 0, 
     &		global_xslice_id(myid,i), MPI_COMM_WORLD, 
     &		request(i), ierr)
	end do
	call MPI_WAITALL(no_of_xslices(myid), request, status, ierr)
      end if
c
      end subroutine  distrib_x_4D_v1

c-----------------------------------------------------------------
c  Master distributes a 4D array of data in y-slice format to all yworkers
c  y-slice by y-slice
c-----------------------------------------------------------------
      subroutine distrib_y_4D_v1(ix,iy,iz,ixloc,Ns,s,s_y)
c      
      implicit none
c
      integer :: ix, iy, iz, ixloc, Ns
      real :: s(ix,iy,iz,Ns)
      real :: s_y(ixloc,iy,iz,Ns)
      integer :: i, ierr
      integer :: status(MPI_STATUS_SIZE, ixloc)
      integer :: request(ixloc)
c      
      if (Master) then
        do i=1,ix
          call MPI_SEND(s(i,1,1,1), 1, GLOBAL_4D_YSLICE, 
     &        owner_of_yslice(i), i,  MPI_COMM_WORLD, ierr)
	end do
      else if (YWorker) then
        do i=1, no_of_yslices(myid)
          call MPI_IRECV(s_y(i,1,1,1), 1, LOCAL_4D_YSLICE, 0,  
     &		global_yslice_id(myid,i), MPI_COMM_WORLD, 
     &		request(i), ierr)
	end do
	call MPI_WAITALL(no_of_yslices(myid), request, status, ierr)
      end if
c
      end subroutine  distrib_y_4D_v1
      
c-----------------------------------------------------------------
c  Master gathers a 4D array of data in x-slice format from all xworkers
c   xslice by xslice, standard version
c-----------------------------------------------------------------
      subroutine gather_x_4D_v1(ix,iy,iz,iyloc,Ns,s,s_x)
c      
      implicit none
c
      integer :: ix, iy, iz, iyloc, Ns
      real :: s(ix,iy,iz,Ns)
      real :: s_x(ix,iyloc,iz,Ns)
      integer :: i, ierr
      integer,pointer,dimension(:,:),save :: status
      integer,pointer,dimension(:),save :: request
      
c      
      if (Master) then
        if (.not. allocated(status)) then
	  allocate( status(MPI_STATUS_SIZE,iy),STAT=ierr )
	end if
	if (.not. allocated(request)) then
          allocate( request(iy),STAT=ierr )      
	end if
        do i=1,iy
        call MPI_IRECV(s(1,i,1,1), 1, GLOBAL_4D_XSLICE, 
     &      owner_of_xslice(i), i,  MPI_COMM_WORLD, request(i), ierr)
	end do
	call MPI_WAITALL(iy, request, status, ierr)
      else if (XWorker) then
        do i=1, no_of_xslices(myid)
        call MPI_SEND(s_x(1,i,1,1), 1, LOCAL_4D_XSLICE, 0, 
     &      global_xslice_id(myid,i), MPI_COMM_WORLD, ierr)
	end do
      end if
c
      end subroutine  gather_x_4D_v1
        
c-----------------------------------------------------------------
c  Master gathers a 4D array of data in x-slice format from all yworkers
c   yslice by yslice, standard version
c-----------------------------------------------------------------
      subroutine gather_y_4D_v1(ix,iy,iz,ixloc,Ns,s,s_y)
c      
      implicit none

      integer :: ix, iy, iz, ixloc, Ns
      real :: s(ix,iy,iz,Ns)
      real :: s_y(ixloc,iy,iz,Ns)
      integer :: i, ierr
      integer,pointer,dimension(:,:),save :: status
      integer,pointer,dimension(:),save :: request

      if (Master) then
        if (.not. allocated(status)) then
          allocate( status(MPI_STATUS_SIZE,ix),STAT=ierr )
	end if
	if (.not. allocated(request)) then
          allocate( request(ix),STAT=ierr )      
	end if
        do i=1,ix
          call MPI_IRECV(s(i,1,1,1), 1, GLOBAL_4D_YSLICE, 
     &        owner_of_yslice(i), i,  
     &	      MPI_COMM_WORLD, request(i), ierr)
	end do
	call MPI_WAITALL(ix, request, status, ierr)
      else if (YWorker) then
        do i=1, no_of_yslices(myid)
          call MPI_SEND(s_y(i,1,1,1), 1, LOCAL_4D_YSLICE, 0, 
     &        global_yslice_id(myid,i), MPI_COMM_WORLD, ierr)
	end do
      end if

      end subroutine  gather_y_4D_v1

c---------------------------------------------------------
c      Slaves exchange data:  - old format is x-slices (for sending)
c                             - new format is y-slices (for receiving)
c---------------------------------------------------------
      subroutine shuffle_x2y_4D_v1(ix,iy,iz,ixloc,iyloc,Ns,s_x,s_y)  
c      
      implicit none
c
      integer :: ix, iy, iz, ixloc, iyloc, Ns
      real   :: s(ix,iy,iz,Ns)
      real   :: s_x(ix,iyloc,iz,Ns)
      real   :: s_y(ixloc,iy,iz,Ns)
      integer :: i, j, Ierr
      integer,pointer,dimension(:,:),save :: status
      integer,pointer,dimension(:),save :: request
c      
      
      if (YWorker) then
        if (.not.allocated(status)) then
	  allocate( status(MPI_STATUS_SIZE,iy*no_of_yslices(myid)),
     &		STAT=ierr )
        end if
	if (.not.allocated(request)) then
          allocate( request(iy*no_of_yslices(myid)),STAT=ierr )
	end if

	do j=1,iy   
          do i=1,no_of_yslices(myid)   
             call MPI_IRECV(s_y(i,j,1,1), 1, LOCAL_4D_YCOLUMN,
     &           owner_of_xslice(j)-1, global_yslice_id(myid,i)+ix*j, 
     &           MPI_COMM_WORKERS, request(j+iy*(i-1)), ierr)
          end do
        end do
      end if

      if (XWorker) then
c      
        do j=1,no_of_xslices(myid)   
          do i=1,ix   
            call MPI_SEND(s_x(i,j,1,1), 1, LOCAL_4D_XCOLUMN, 
     &          owner_of_yslice(i)-1, i+ix*global_xslice_id(myid,j),  
     &          MPI_COMM_WORKERS, ierr)
          end do
        end do
      end if
      
      if (YWorker) then
	call MPI_WAITALL(iy*no_of_yslices(myid), request, status, ierr)
      end if
c	
      end subroutine shuffle_x2y_4D_v1
	
c---------------------------------------------------------
c      Slaves exchange data:  - old format is y-slices (for sending)
c                             - new format is x-slices (for receiving)
c---------------------------------------------------------
      subroutine shuffle_y2x_4D_v1(ix,iy,iz,ixloc,iyloc,Ns,s_x,s_y)  
c      
      implicit none
c
      integer :: ix, iy, iz, ixloc, iyloc, Ns
      real   :: s(ix,iy,iz,Ns)
      real   :: s_x(ix,iyloc,iz,Ns)
      real   :: s_y(ixloc,iy,iz,Ns)
      integer :: i, j, ierr
      integer,pointer,dimension(:,:),save :: status
      integer,pointer,dimension(:),save :: request
c      
c      
      if (XWorker) then
        if (.not.allocated(status)) then
  	  allocate( status(MPI_STATUS_SIZE,ix*no_of_xslices(myid)),
     &		STAT=ierr )
        end if
	if (.not.allocated(request)) then
          allocate( request(ix*no_of_xslices(myid)),STAT=ierr )      
	end if

        do j=1,no_of_xslices(myid)   
	  do i=1,ix   
            call MPI_IRECV(s_x(i,j,1,1), 1, LOCAL_4D_XCOLUMN, 
     &             owner_of_yslice(i)-1, i+ix*global_xslice_id(myid,j),
     &             MPI_COMM_WORKERS, request(i+ix*(j-1)), ierr)
          end do
        end do
	
      end if

      if (YWorker) then  ! Master has no role in this communication
	do j=1,iy   
          do i=1,no_of_yslices(myid)   
            call MPI_SEND(s_y(i,j,1,1), 1, LOCAL_4D_YCOLUMN, 
     &          owner_of_xslice(j)-1, global_yslice_id(myid,i)+j*ix,  
     &          MPI_COMM_WORKERS, ierr)
          end do
        end do
      end if

      if (XWorker) then
	call MPI_WAITALL(ix*no_of_xslices(myid), request, status, ierr)
      end if

      end subroutine shuffle_y2x_4D_v1
      
c--------------------------------------------------------------------------
c  Master distributes a 2D (ix.iy) array of data in x-slice format 
c  to all workers
c--------------------------------------------------------------------------
c      
      subroutine distrib_x_2D_v1(ix,iy,iyloc,s,s_x)

      implicit none
c
      integer :: ix, iy, iyloc
      real :: s(ix,iy)
      real :: s_x(ix,iyloc)
      integer :: j, ierr
      integer :: request(iyloc)
      integer :: status(MPI_STATUS_SIZE, iyloc)
c      
      if (Master) then
        do j=1,iy
        call MPI_SEND(s(1,j), ix, MPI_REAL, 
     &         owner_of_xslice(j), j,  MPI_COMM_WORLD, ierr)
	end do
      else if (XWorker) then
        do j=1, no_of_xslices(myid)
        call MPI_IRECV(s_x(1,j), ix, MPI_REAL, 0, 
     &			global_xslice_id(myid,j), 
     &                  MPI_COMM_WORLD, request(j), ierr)
	end do
	call MPI_WAITALL(no_of_xslices(myid), request, status, ierr)
      end if
c
      end subroutine  distrib_x_2D_v1
      
        
c---------------------------------------------------------------------------
c  Master distributes a 2D  (ix.iy) array of data in y-slice format 
c  to all workers
c---------------------------------------------------------------------------
c      
      subroutine distrib_y_2D_v1(ix,iy,ixloc,s,s_y)

      implicit none
c
      integer :: ix, iy, ixloc
      real :: s(ix,iy)
      real :: s_y(ixloc,iy)
      integer :: i, ierr
      integer :: request(ixloc)
      integer :: status(MPI_STATUS_SIZE, ixloc)
c      
      if (Master) then
        do i=1,ix
          call MPI_SEND(s(i,1), 1, GLOBAL_2D_YSLICE, 
     &         owner_of_yslice(i), i,  MPI_COMM_WORLD, ierr)
	end do
      else if (YWorker) then
        do i=1, no_of_yslices(myid)
          call MPI_IRECV(s_y(i,1), 1, LOCAL_2D_YSLICE, 0, 
     &			global_yslice_id(myid,i), 
     &                 MPI_COMM_WORLD, request(i), ierr)
	end do
	call MPI_WAITALL(no_of_yslices(myid), request, status, ierr)
      end if
c
      end subroutine  distrib_y_2D_v1
        
c      
c-------------------------------------------------------------------------
c  Master distributes a BDz (ix.iy.is) array of data in x-slice format 
c   to all workers
c-------------------------------------------------------------------------
c      
      subroutine distrib_x_BDz_v1(ix,iy,iyloc,Ns,s,s_x)

      implicit none
c
      integer :: ix, iy, Ns, iyloc
      real   :: s(ix,iy,Ns)
      real   :: s_x(ix,iyloc,Ns)
      integer :: j, ierr
      integer :: request(iyloc)
      integer :: status(MPI_STATUS_SIZE, iyloc)
c      
      if (Master) then
        do j=1,iy
        call MPI_SEND(s(1,j,1), 1, GLOBAL_BDz_XSLICE, 
     &         owner_of_xslice(j), j,  MPI_COMM_WORLD, ierr)
	end do
      else if (XWorker) then
        do j=1, no_of_xslices(myid)
        call MPI_IRECV(s_x(1,j,1), 1, LOCAL_BDz_XSLICE, 0, 
     &			global_xslice_id(myid,j), 
     &          	MPI_COMM_WORLD, request(j), ierr)
	end do
	call MPI_WAITALL(no_of_xslices(myid), request, status, ierr)
      end if
c
      end subroutine  distrib_x_BDz_v1
      
c-------------------------------------------------------------------------
c  Master distributes a BDz (ix.iy.is) array of data in y-slice format 
c   to all workers
c-------------------------------------------------------------------------
c      
      subroutine distrib_y_BDz_v1(ix,iy,ixloc,Ns,s,s_y)

      implicit none
c
      integer :: ix, iy, Ns, ixloc
      real :: s(ix,iy,Ns)
      real :: s_y(ixloc,iy,Ns)
      integer :: i, ierr
      integer :: request(ixloc)
      integer :: status(MPI_STATUS_SIZE, ixloc)
c      
      if (Master) then
        do i=1,ix
        call MPI_SEND(s(i,1,1), 1, GLOBAL_BDz_YSLICE, 
     &         owner_of_yslice(i), i,  MPI_COMM_WORLD, ierr)
	end do
      else if (YWorker) then
        do i=1, no_of_yslices(myid)
        call MPI_IRECV(s_y(i,1,1), 1, LOCAL_BDz_YSLICE, 0, 
     &			global_yslice_id(myid,i), 
     &                 MPI_COMM_WORLD, request(i), ierr)
	end do
	call MPI_WAITALL(no_of_yslices(myid), request, status, ierr)
      end if
c
      end subroutine  distrib_y_BDz_v1



              
c-------------------------------------------------------------------------
c  Master distributes a 3D (ix.iy.iz) array of data in x-slice format 
c   to all workers
c-------------------------------------------------------------------------
c      
      subroutine distrib_x_3D_v1(ix,iy,iz,iyloc,s,s_x)

      implicit none
c
      integer :: ix, iy, iz, iyloc
      real :: s(ix,iy,iz)
      real :: s_x(ix,iyloc,iz)
      integer :: i, ierr
      integer :: request(iyloc)
      integer :: status(MPI_STATUS_SIZE, iyloc)
c      
      if (Master) then
        do i=1,iy
          call MPI_SEND(s(1,i,1), 1, GLOBAL_3D_XSLICE, 
     &         owner_of_xslice(i), i,  MPI_COMM_WORLD, ierr)
          if (ierr.ne.MPI_SUCCESS) then
	    print*, 'Error: distrib_x_3D: send failed'
	    stop
	  end if
	end do
      else if (XWorker) then
        do i=1, no_of_xslices(myid)
          call MPI_IRECV(s_x(1,i,1), 1, LOCAL_3D_XSLICE, 0, 
     &			global_xslice_id(myid,i), 
     &                 MPI_COMM_WORLD, request(i), ierr)
          if (ierr.ne.MPI_SUCCESS) then
	    print*, 'Error: distrib_x_3D: receive failed'
	    stop
	  end if
	end do
	call MPI_WAITALL(no_of_xslices(myid), request, status, ierr)
      end if
c
      end subroutine  distrib_x_3D_v1

      
c--------------------------------------------------------------------------
c  Master distributes a 3D (ix.iy.iz) array of data in y-slice format 
c   to all workers
c--------------------------------------------------------------------------
c
      subroutine distrib_y_3D_v1(ix,iy,iz,ixloc,sg1,sgy)

      implicit none
c      
      integer :: ix, iy, iz, ixloc
      real   :: sg1(ix,iy,iz)
      real   :: sgy(ixloc,iy,iz)
      integer :: i, ierr
      integer :: request(ixloc)
      integer :: status(MPI_STATUS_SIZE, ixloc)
c      
      if (Master) then
        do i=1,ix
           call MPI_SEND(sg1(i,1,1), 1, GLOBAL_3D_YSLICE,
     &         owner_of_yslice(i), i,  MPI_COMM_WORLD, ierr)
        end do
      else if (YWorker) then
        do i=1, no_of_yslices(myid)
           call MPI_IRECV(sgy(i,1,1), 1, LOCAL_3D_YSLICE, 0, 
     &			global_yslice_id(myid,i), 
     &                 MPI_COMM_WORLD, request(i), ierr)
        end do
	call MPI_WAITALL(no_of_yslices(myid), request, status, ierr)
      endif
c	
      return
      end subroutine distrib_y_3D_v1

c---------------------------------------------------------------------------
c  Master distributes the x-boundary array  BDx (only for x-partitioning)
c---------------------------------------------------------------------------
c      
      subroutine distrib_x_BDx_v1(iy,iz,iyloc,Ns,s,sloc)

      implicit none
c
      integer :: iy, iz, Ns, iyloc
      real :: s(iy,iz,Ns)
      real :: sloc(iyloc,iz,Ns)
      integer :: j, ierr
      integer :: request(iyloc)
      integer :: status(MPI_STATUS_SIZE, iyloc)
c      
      if (Master) then
        do j=1,iy
        call MPI_SEND(s(j,1,1), 1, GLOBAL_BDx_XSLICE, 
     &         owner_of_xslice(j), j,  MPI_COMM_WORLD, ierr)
	end do
      else if (XWorker) then
        do j=1, no_of_xslices(myid)
        call MPI_IRECV(sloc(j,1,1), 1, LOCAL_BDx_XSLICE, 0, 
     &			global_xslice_id(myid,j), 
     &                 MPI_COMM_WORLD, request(j), ierr)
	end do
	call MPI_WAITALL(no_of_xslices(myid), request, status, ierr)
      end if
c
      end subroutine  distrib_x_BDx_v1

c---------------------------------------------------------------------------
c  Master distributes the x-boundary array BDx2
c---------------------------------------------------------------------------
c      
      subroutine distrib_x_BD_v1(iy,iz,iyloc,Ns,s,sloc)

      implicit none
c
      integer :: iy, iz, Ns, iyloc
      real :: s(iy,iz,2,Ns)
      real :: sloc(iyloc,iz,2,Ns)
      integer :: j, ierr
      integer :: request(iyloc)
      integer :: status(MPI_STATUS_SIZE, iyloc)
c      
      if (Master) then
        do j=1,iy
        call MPI_SEND(s(j,1,1,1), 1, GLOBAL_BD_XSLICE, 
     &         owner_of_xslice(j), j,  MPI_COMM_WORLD, ierr)
	end do
      else if (XWorker) then
        do j=1, no_of_xslices(myid)
        call MPI_IRECV(sloc(j,1,1,1), 1, LOCAL_BD_XSLICE, 0, 
     &			global_xslice_id(myid,j), 
     &                 MPI_COMM_WORLD, request(j), ierr)
	end do
	call MPI_WAITALL(no_of_xslices(myid), request, status, ierr)
      end if
c
      end subroutine  distrib_x_BD_v1
      
c---------------------------------------------------------------------------
c  Master distributes the y-boundary array  BDy (only for y-partitioning)
c---------------------------------------------------------------------------
c      
      subroutine distrib_y_BDy_v1(ix,iz,ixloc,Ns,s,sloc)

      implicit none
c
      integer :: ix, iz, Ns, ixloc
      real :: s(ix,iz,Ns)
      real :: sloc(ixloc,iz,Ns)
      integer :: i, ierr
      integer :: request(ixloc)
      integer :: status(MPI_STATUS_SIZE, ixloc)
c      
      if (Master) then
        do i=1,ix
        call MPI_SEND(s(i,1,1), 1, GLOBAL_BDy_YSLICE, 
     &         owner_of_yslice(i), i,  MPI_COMM_WORLD, ierr)
	end do
      else if (YWorker) then
        do i=1, no_of_yslices(myid)
        call MPI_IRECV(sloc(i,1,1), 1, LOCAL_BDy_YSLICE, 0,
     &			global_yslice_id(myid,i), 
     &                 MPI_COMM_WORLD, request(i), ierr)
	end do
	call MPI_WAITALL(no_of_yslices(myid), request, status, ierr)
      end if
c
      end subroutine  distrib_y_BDy_v1

c---------------------------------------------------------------------------
c  Master distributes the y-boundary array BDy2 (Stem-III dependent)
c---------------------------------------------------------------------------
c      
      subroutine distrib_y_BD_v1(ix,iz,ixloc,Ns,s,sloc)

      implicit none
c
      integer :: ix, iz, Ns, ixloc
      real :: s(ix,iz,2,Ns)
      real :: sloc(ixloc,iz,2,Ns)
      integer :: i, ierr
      integer :: request(ixloc)
      integer :: status(MPI_STATUS_SIZE, ixloc)
c      
      if (Master) then
        do i=1,ix
        call MPI_SEND(s(i,1,1,1), 1, GLOBAL_BD_YSLICE, 
     &         owner_of_yslice(i), i,  MPI_COMM_WORLD, ierr)
	end do
      else if (YWorker) then
        do i=1, no_of_yslices(myid)
        call MPI_IRECV(sloc(i,1,1,1), 1, LOCAL_BD_YSLICE, 0,
     &			global_yslice_id(myid,i), 
     &                 MPI_COMM_WORLD, request(i), ierr)
	end do
	call MPI_WAITALL(no_of_yslices(myid), request, status, ierr)
      end if
c
      end subroutine  distrib_y_BD_v1
      
c *******************************************************************
c ***************** single communication ****************************
c *******************************************************************

c-----------------------------------------------------------------
c  Master distributes a 4D array of data in x-slice format to all workers
c  all at once in one message
c-----------------------------------------------------------------
      subroutine distrib_x_4D_v2(ix,iy,iz,iyloc,Ns,s,s_x)
c      
      implicit none
c
      integer :: ix, iy, iz, iyloc, Ns
      real :: s(ix,iy,iz,Ns)
      real :: s_x(ix,iyloc,iz,Ns)
      integer :: i, ierr, status(MPI_STATUS_SIZE)

c The master sends the x-slice arrays to the workers using one
c  message for each worker, depending on the array size
      if (Master) then
        do i=1,NXworkers
	  if (no_of_xslices(i) .eq. iyloc) then
            call MPI_SEND(s(1,i,1,1), 1, GLOBAL_LARGE_4D_XSLICES, 
     &         owner_of_xslice(i), i,  MPI_COMM_WORLD, ierr)
	  else
            call MPI_SEND(s(1,i,1,1), 1, GLOBAL_SMALL_4D_XSLICES, 
     &         owner_of_xslice(i), i,  MPI_COMM_WORLD, ierr)
	  end if
	end do
c The workers receive the x-slice arrays from the master in one
c  message depending on the array size
      else if (XWorker) then
	if (no_of_xslices(myid) .eq. iyloc) then
          call MPI_RECV(s_x(1,1,1,1), 1, LOCAL_LARGE_4D_XSLICES, 0,
     &		myid, MPI_COMM_WORLD, status, ierr)
	else
          call MPI_RECV(s_x(1,1,1,1), 1, LOCAL_SMALL_4D_XSLICES, 0,
     &		myid, MPI_COMM_WORLD, status, ierr)
        end if
      end if
c
      end subroutine  distrib_x_4D_v2

c-----------------------------------------------------------------
c  Master distributes a 4D array of data in y-slice format to all workers
c  all in one message
c-----------------------------------------------------------------
      subroutine distrib_y_4D_v2(ix,iy,iz,ixloc,Ns,s,s_y)
c      
      implicit none
c
      integer :: ix, iy, iz, ixloc, Ns
      real :: s(ix,iy,iz,Ns)
      real :: s_y(ixloc,iy,iz,Ns)
      integer :: i, ierr, status(MPI_STATUS_SIZE)
c The master sends one y-slice array to each of the workers using
c  one message, depending on the size.
      if (Master) then
        do i=1,NYworkers
	  if (no_of_yslices(i) .eq. ixloc) then
	    call MPI_SEND(s(i,1,1,1), 1, GLOBAL_LARGE_4D_YSLICES, 
     &         owner_of_yslice(i), i,  MPI_COMM_WORLD, ierr)
	  else
	    call MPI_SEND(s(i,1,1,1), 1, GLOBAL_SMALL_4D_YSLICES, 
     &         owner_of_yslice(i), i,  MPI_COMM_WORLD, ierr)
          end if
	end do
c The workers receive their y-slice array at once depending on
c  the size
      else if (YWorker) then
	if (no_of_yslices(myid) .eq. ixloc) then
          call MPI_RECV(s_y(1,1,1,1), 1, LOCAL_LARGE_4D_YSLICES, 0, 
     &		myid, MPI_COMM_WORLD, status, ierr)
	else
          call MPI_RECV(s_y(1,1,1,1), 1, LOCAL_SMALL_4D_YSLICES, 0, 
     &		myid, MPI_COMM_WORLD, status, ierr)
        end if
      end if
c
      end subroutine  distrib_y_4D_v2
      
c-----------------------------------------------------------------
c  Master gathers a 4D array of data in x-slice format from all workers
c-----------------------------------------------------------------
      subroutine gather_x_4D_v2(ix,iy,iz,iyloc,Ns,s,s_x)
c      
      implicit none
c
      integer :: ix, iy, iz, iyloc, Ns
      real   :: s(ix,iy,iz,Ns)
      real   :: s_x(ix,iyloc,iz,Ns)
      integer :: i, ierr
      integer,pointer,dimension(:,:),save :: status
      integer,pointer,dimension(:),save :: request
c      
c The master receives a the x-slice array from each of the workers
c  and places it right in the global array, depending on the size
      if (Master) then
	if (.not. allocated(status)) then
          allocate( status(MPI_STATUS_SIZE,NXworkers),STAT=ierr )
	end if
	if (.not. allocated(request)) then
          allocate( request(NXworkers),STAT=ierr )      
	end if
        do i=1,NXworkers
	  if ( no_of_xslices(i) .eq. iyloc) then
	    call MPI_IRECV(s(1,i,1,1), 1, GLOBAL_LARGE_4D_XSLICES, 
     &         owner_of_xslice(i), i,  MPI_COMM_WORLD, 
     &         request(i), ierr)
     	  else
	    call MPI_IRECV(s(1,i,1,1), 1, GLOBAL_SMALL_4D_XSLICES, 
     &         owner_of_xslice(i), i,  MPI_COMM_WORLD, 
     &         request(i), ierr)
     	  end if
	end do
	call MPI_WAITALL(NXworkers, request, status, ierr)
	
c The xworkers send their x-slice array to the master in one message
c  depending on the size
      else if (XWorker) then
      	if (no_of_xslices(myid) .eq. iyloc) then
	  call MPI_SEND(s_x(1,1,1,1), 1, LOCAL_LARGE_4D_XSLICES, 
     &                 0, myid, 
     &                 MPI_COMM_WORLD, ierr)
        else
	  call MPI_SEND(s_x(1,1,1,1), 1, LOCAL_SMALL_4D_XSLICES, 
     &                 0, myid, 
     &                 MPI_COMM_WORLD, ierr)
     	end if
      end if
c
      end subroutine  gather_x_4D_v2  

c-----------------------------------------------------------------
c  Master gathers a 4D array of data in y-slice format from all yworkers
c-----------------------------------------------------------------
      subroutine gather_y_4D_v2(ix,iy,iz,ixloc,Ns,s,s_y)
c      
      implicit none
c
      integer :: ix, iy, iz, ixloc, Ns
      real   :: s(ix,iy,iz,Ns)
      real   :: s_y(ixloc,iy,iz,Ns)
      integer :: i, ierr
      integer,pointer,dimension(:,:),save :: status
      integer,pointer,dimension(:),save :: request
c      
c The master receives the y-slice array from each yworker and places
c  it in the right position, depending on the size
      if (Master) then
        if (.not. allocated(status)) then
          allocate( status(MPI_STATUS_SIZE,NYworkers),STAT=ierr )
	end if
	if (.not. allocated(request)) then
          allocate( request(NYworkers),STAT=ierr )      
	end if
        do i=1,NYworkers
	  if ( no_of_yslices(i) .eq. ixloc) then
	    call MPI_IRECV(s(i,1,1,1), 1, GLOBAL_LARGE_4D_YSLICES, 
     &         owner_of_yslice(i), i,  MPI_COMM_WORLD, 
     &         request(i), ierr)
     	  else
	    call MPI_IRECV(s(i,1,1,1), 1, GLOBAL_SMALL_4D_YSLICES, 
     &         owner_of_yslice(i), i,  MPI_COMM_WORLD, 
     &         request(i), ierr)
     	  end if
	end do

	call MPI_WAITALL(NYworkers, request, status, ierr)

c Each yworker sends the y-slice array to the master in one message
c  depending on the size
      else if (YWorker) then
      	if ( no_of_yslices(myid) .eq. ixloc ) then
          call MPI_SEND(s_y(1,1,1,1), 1, LOCAL_LARGE_4D_YSLICES, 
     &                 0, myid, 
     &                 MPI_COMM_WORLD, ierr)
        else
          call MPI_SEND(s_y(1,1,1,1), 1, LOCAL_SMALL_4D_YSLICES, 
     &                 0, myid, 
     &                 MPI_COMM_WORLD, ierr)
	end if
      end if
c
      end subroutine  gather_y_4D_v2  
        
c---------------------------------------------------------
c      Slaves exchange data:  - old format is x-slices (for sending)
c                             - new format is y-slices (for receiving)
c---------------------------------------------------------
      subroutine shuffle_x2y_4D_v2(ix,iy,iz,ixloc,iyloc,Ns,s_x,s_y)  
c      
      implicit none
c
      integer :: ix, iy, iz, ixloc, iyloc, Ns
      real   :: s_x(ix,iyloc,iz,Ns)
      real   :: s_y(ixloc,iy,iz,Ns)
      integer :: i, j, ierr
      integer,pointer,dimension(:,:),save :: status
      integer,pointer,dimension(:),save :: request

c      
c Each yworker receives one message from all the xworkers combining
c  all the data needed in one message, depending on x-slice
c  size at the sender and y-slice size at the receiver      
      if (YWorker) then
	if (.not. allocated(status)) then
          allocate( status(MPI_STATUS_SIZE,NXworkers),STAT=ierr )
	end if
	if (.not. allocated(request)) then
          allocate( request(NXworkers),STAT=ierr )      
	end if

	if (no_of_yslices(myid) .eq. ixloc) then
          do j=1, NXworkers ! ID's in worker group only
	  if (no_of_xslices(j) .eq. iyloc) then
c	print*,'proc',myid,'receives LXLY from',j
     	    call MPI_IRECV(s_y(1,j,1,1),1, LOCAL_4D_LXLY_YCOLS,
     &		j-1, myid, MPI_COMM_WORKERS, request(j), ierr)
	  else
c	print*,'proc',myid,'receives SXLY from',j
     	    call MPI_IRECV(s_y(1,j,1,1),1,LOCAL_4D_SXLY_YCOLS,
     &		j-1, myid, MPI_COMM_WORKERS, request(j), ierr)
	  end if
	  end do
	else
          do j=1, NXworkers ! ID's in worker group only
	  if (no_of_xslices(j) .eq. iyloc) then	
c	print*,'proc',myid,'receives LXSY from',j
     	    call MPI_IRECV(s_y(1,j,1,1),1,LOCAL_4D_LXSY_YCOLS,
     &		j-1, myid, MPI_COMM_WORKERS, request(j), ierr)
          else
c	print*,'proc',myid,'receives SXSY from',j
      	    call MPI_IRECV(s_y(1,j,1,1),1,LOCAL_4D_SXSY_YCOLS,
     &		j-1, myid, MPI_COMM_WORKERS, request(j), ierr)
	  end if
	  end do
	end if
      end if

      if (XWorker) then
c Each xworker sends one message to each
c  of the yworkers combining all the data that need
c  to be sent there in one message, depending on x-slice
c  size at the sender and y-slice size at the receiver      
	if (no_of_xslices(myid) .eq. iyloc) then
          do i=1, NYworkers ! ID's in worker group only
	  if (no_of_yslices(i) .eq. ixloc) then
c	print*,'proc',myid,'sends LXLY to',i
       	    call MPI_SEND(s_x(i,1,1,1),1,LOCAL_4D_LXLY_XCOLS,
     &		i-1, i, MPI_COMM_WORKERS, ierr)	  	
	  else
c	print*,'proc',myid,'sends LXSY to',i
       	    call MPI_SEND(s_x(i,1,1,1),1,LOCAL_4D_LXSY_XCOLS,
     &		i-1, i, MPI_COMM_WORKERS, ierr)
	  end if
	  end do
	else
          do i=1, NYworkers ! ID's in worker group only
	  if (no_of_yslices(i) .eq. ixloc) then	
c	print*,'proc',myid,'sends SXLY to',i
       	    call MPI_SEND(s_x(i,1,1,1),1,LOCAL_4D_SXLY_XCOLS,
     &		i-1, i, MPI_COMM_WORKERS, ierr)
          else
c	print*,'proc',myid,'sends SXSY to',i
       	    call MPI_SEND(s_x(i,1,1,1),1,LOCAL_4D_SXSY_XCOLS,
     &		i-1, i, MPI_COMM_WORKERS, ierr)
 	  end if
          end do
	end if
      end if 

      if (YWorker) call MPI_WAITALL(NXworkers, request, status, ierr)
	
      end subroutine shuffle_x2y_4D_v2

c---------------------------------------------------------
c      Slaves exchange data:  - old format is y-slices (for sending)
c                             - new format is x-slices (for receiving)
c---------------------------------------------------------
      subroutine shuffle_y2x_4D_v2(ix,iy,iz,ixloc,iyloc,Ns,s_x,s_y)  
c      
      implicit none
c
      integer :: ix, iy, iz, ixloc, iyloc, Ns
      real   :: s_x(ix,iyloc,iz,Ns)
      real   :: s_y(ixloc,iy,iz,Ns)
      integer :: i, j, ierr
      integer,pointer,dimension(:,:),save :: status
      integer,pointer,dimension(:),save :: request
c
    
      if (XWorker) then
c Each xworker receives one message from all the yworkers 
c  combining all the data needed in one message, depending 
c  on x-slice size at the sender and y-slice size at the 
c  receiver      
      if(.not. allocated(status)) then
        allocate( status(MPI_STATUS_SIZE,NYworkers),STAT=ierr )
      end if
      if (.not. allocated(request)) then
        allocate( request(NYworkers),STAT=ierr )      
      end if

      if (no_of_xslices(myid) .eq. iyloc) then
        do j=1, NYworkers
	  if (no_of_yslices(j) .eq. ixloc) then
c	print*, 'P', myid, ' receives LXLY from ', j
     	    call MPI_IRECV(s_x(j,1,1,1),1, LOCAL_4D_LXLY_XCOLS,
     &		j-1, myid, MPI_COMM_WORKERS, request(j), ierr)
     	  else
c	print*, 'P', myid, ' receives LXSY from ', j
     	    call MPI_IRECV(s_x(j,1,1,1),1,LOCAL_4D_LXSY_XCOLS,
     &		j-1, myid, MPI_COMM_WORKERS, request(j), ierr)
     	  end if
	end do
      else
	do j=1,NYworkers
	  if (no_of_yslices(j) .eq. ixloc) then	
c	print*, 'P', myid, ' receives SXLY from ', j
     	    call MPI_IRECV(s_x(j,1,1,1),1,LOCAL_4D_SXLY_XCOLS,
     &		j-1, myid, MPI_COMM_WORKERS, request(j), ierr)
     	  else
c	print*, 'P', myid, ' receives SXSY from ', j
     	    call MPI_IRECV(s_x(j,1,1,1),1,LOCAL_4D_SXSY_XCOLS,
     &		j-1, myid, MPI_COMM_WORKERS, request(j), ierr)
     	  end if
	end do
      end if

      end if
      
      if (YWorker) then
c Each yworker sends one message to each
c  of the xworkers combining all the data that needs
c  to be sent there in one message, depending on y-slice
c  size at the sender and x-slice size at the receiver      

      if (no_of_yslices(myid) .eq. ixloc) then
        do i=1, NXworkers
	  if (no_of_xslices(i) .eq. iyloc) then
c	print*, 'P', myid, ' sends LXLY to ', i
       	    call MPI_SEND(s_y(1,i,1,1),1,LOCAL_4D_LXLY_YCOLS,
     &		i-1, i, MPI_COMM_WORKERS, ierr)	  	
     	  else
c	print*, 'P', myid, ' sends SXLY to ', i
       	    call MPI_SEND(s_y(1,i,1,1),1,LOCAL_4D_SXLY_YCOLS,
     &		i-1, i, MPI_COMM_WORKERS, ierr)
          end if
     	end do
      else
	do i=1,NXworkers
	  if (no_of_xslices(i) .eq. iyloc) then	
c	print*, 'P', myid, ' sends LXSY to ', i
       	    call MPI_SEND(s_y(1,i,1,1),1,LOCAL_4D_LXSY_YCOLS,
     &		i-1, i, MPI_COMM_WORKERS, ierr)
     	  else
c	print*, 'P', myid, ' sends SXSY to ', i
       	    call MPI_SEND(s_y(1,i,1,1),1,LOCAL_4D_SXSY_YCOLS,
     &		i-1, i, MPI_COMM_WORKERS, ierr)
     	  end if
	end do
      end if
      
      end if

      if (XWorker) call MPI_WAITALL(NYworkers, request, status, ierr)
      
      end subroutine shuffle_y2x_4D_v2
	
c--------------------------------------------------------------------------
c  Master distributes a 2D (ix.iy) array of data in x-slice format 
c  to all xworkers
c--------------------------------------------------------------------------
c      
      subroutine distrib_x_2D_v2(ix,iy,iyloc,s,s_x)

      implicit none
c
      integer :: ix, iy, iyloc
      real :: s(ix,iy)
      real :: s_x(ix,iyloc)
      integer :: j, ierr, status(MPI_STATUS_SIZE)
c      
      if (Master) then
c The master sends the x-slice arrays to the xworkers using one
c  message for each xworker, depending on the array size
        do j=1,NXworkers
	  if (no_of_xslices(j) .eq. iyloc) then
            call MPI_SEND(s(1,j), 1, GLOBAL_LARGE_2D_XSLICES, 
     &         owner_of_xslice(j), j,  MPI_COMM_WORLD, ierr)
          else
            call MPI_SEND(s(1,j), 1, GLOBAL_SMALL_2D_XSLICES, 
     &         owner_of_xslice(j), j,  MPI_COMM_WORLD, ierr)
          end if
	end do
      else if (XWorker) then
c The xworkers receive the x-slice arrays from the master in one
c  message depending on the array size
        if (no_of_xslices(myid) .eq. iyloc) then
          call MPI_RECV(s_x(1,1), 1, LOCAL_LARGE_2D_XSLICES, 0, 
     &			myid, MPI_COMM_WORLD, status, ierr)
        else
          call MPI_RECV(s_x(1,1), 1, LOCAL_SMALL_2D_XSLICES, 0, 
     &			myid, MPI_COMM_WORLD, status, ierr)
	end if	  
      end if
c
      end subroutine  distrib_x_2D_v2
      
        
c---------------------------------------------------------------------------
c  Master distributes a 2D  (ix.iy) array of data in y-slice format 
c  to all yworkers
c---------------------------------------------------------------------------
c      
      subroutine distrib_y_2D_v2(ix,iy,ixloc,s,s_y)

      implicit none
c
      integer :: ix, iy, ixloc
      real :: s(ix,iy)
      real :: s_y(ixloc,iy)
      integer :: i, ierr, status(MPI_STATUS_SIZE)
c      
      if (Master) then
c The master sends one y-slice array to each of the yworkers using
c  one message, depending on the size.
        do i=1,NYworkers
	  if (no_of_yslices(i) .eq. ixloc) then
            call MPI_SEND(s(i,1), 1, GLOBAL_LARGE_2D_YSLICES, 
     &         owner_of_yslice(i), i,  MPI_COMM_WORLD, ierr)
          else
            call MPI_SEND(s(i,1), 1, GLOBAL_SMALL_2D_YSLICES, 
     &         owner_of_yslice(i), i,  MPI_COMM_WORLD, ierr)
          end if
	end do
      else if (YWorker) then
c The yworkers receive their y-slice array at once depending on
c  the size
        if (no_of_yslices(myid) .eq. ixloc) then
          call MPI_RECV(s_y(1,1), 1, LOCAL_LARGE_2D_YSLICES, 0, 
     &			myid, MPI_COMM_WORLD, status, ierr)
	else
          call MPI_RECV(s_y(1,1), 1, LOCAL_SMALL_2D_YSLICES, 0, 
     &			myid, MPI_COMM_WORLD, status, ierr)
	end if
      end if
c
      end subroutine  distrib_y_2D_v2
        
c      
c-------------------------------------------------------------------------
c  Master distributes a BDz (ix.iy.is) array of data in x-slice format 
c   to all xworkers
c-------------------------------------------------------------------------
c      
      subroutine distrib_x_BDz_v2(ix,iy,iyloc,Ns,s,s_x)

      implicit none
c
      integer :: ix, iy, Ns, iyloc
      real   :: s(ix,iy,Ns)
      real   :: s_x(ix,iyloc,Ns)
      integer :: j, ierr, status(MPI_STATUS_SIZE)
c      
      if (Master) then
c The master sends the x-slice arrays to the xworkers using one
c  message for each xworker, depending on the array size
        do j=1,NXworkers
	  if (no_of_xslices(j) .eq. iyloc) then
            call MPI_SEND(s(1,j,1), 1, GLOBAL_LARGE_BDz_XSLICES, 
     &         owner_of_xslice(j), j,  MPI_COMM_WORLD, ierr)
          else
            call MPI_SEND(s(1,j,1), 1, GLOBAL_SMALL_BDz_XSLICES, 
     &         owner_of_xslice(j), j,  MPI_COMM_WORLD, ierr)
          end if
	end do
      else if(XWorker) then
c The xworkers receive the x-slice arrays from the master in one
c  message depending on the array size
        if (no_of_xslices(myid) .eq. iyloc) then
          call MPI_RECV(s_x(1,1,1), 1, LOCAL_LARGE_BDz_XSLICES, 0, 
     &			myid, MPI_COMM_WORLD, status, ierr)
        else
          call MPI_RECV(s_x(1,1,1), 1, LOCAL_SMALL_BDz_XSLICES, 0, 
     &			myid, MPI_COMM_WORLD, status, ierr)
	end if	  
      end if
c
      end subroutine  distrib_x_BDz_v2
      
c-------------------------------------------------------------------------
c  Master distributes a BDz (ix.iy.is) array of data in y-slice format 
c   to all yworkers
c-------------------------------------------------------------------------
c      
      subroutine distrib_y_BDz_v2(ix,iy,ixloc,Ns,s,s_y)

      implicit none
c
      integer :: ix, iy, Ns, ixloc
      real :: s(ix,iy,Ns)
      real :: s_y(ixloc,iy,Ns)
      integer :: i, ierr, status(MPI_STATUS_SIZE)
c      
      if (Master) then
c The master sends one y-slice array to each of the yworkers using
c  one message, depending on the size.
        do i=1,NYworkers
	  if (no_of_yslices(i) .eq. ixloc) then
            call MPI_SEND(s(i,1,1), 1, GLOBAL_LARGE_BDz_YSLICES, 
     &         owner_of_yslice(i), i,  MPI_COMM_WORLD, ierr)
          else
            call MPI_SEND(s(i,1,1), 1, GLOBAL_SMALL_BDz_YSLICES, 
     &         owner_of_yslice(i), i,  MPI_COMM_WORLD, ierr)
          end if
	end do
      else if (YWorker) then
c The yworkers receive their y-slice array at once depending on
c  the size
        if (no_of_yslices(myid) .eq. ixloc) then
          call MPI_RECV(s_y(1,1,1), 1, LOCAL_LARGE_BDz_YSLICES, 0, 
     &			myid, MPI_COMM_WORLD, status, ierr)
	else
          call MPI_RECV(s_y(1,1,1), 1, LOCAL_SMALL_BDz_YSLICES, 0, 
     &			myid, MPI_COMM_WORLD, status, ierr)
	end if
      end if
c
      end subroutine  distrib_y_BDz_v2



              
c-------------------------------------------------------------------------
c  Master distributes a 3D (ix.iy.iz) array of data in x-slice format 
c   to all xworkers
c-------------------------------------------------------------------------
c      
      subroutine distrib_x_3D_v2(ix,iy,iz,iyloc,s,s_x)

      implicit none
c
      integer :: ix, iy, iz, iyloc
      real :: s(ix,iy,iz)
      real :: s_x(ix,iyloc,iz)
      integer :: j, ierr, status(MPI_STATUS_SIZE)
c      
      if (Master) then
c The master sends the x-slice arrays to the xworkers using one
c  message for each xworker, depending on the array size
        do j=1,NXworkers
	  if (no_of_xslices(j) .eq. iyloc) then
            call MPI_SEND(s(1,j,1), 1, GLOBAL_LARGE_3D_XSLICES, 
     &         owner_of_xslice(j), j,  MPI_COMM_WORLD, ierr)
          else
            call MPI_SEND(s(1,j,1), 1, GLOBAL_SMALL_3D_XSLICES, 
     &         owner_of_xslice(j), j,  MPI_COMM_WORLD, ierr)
          end if
	end do
      else if (XWorker) then
c The xworkers receive the x-slice arrays from the master in one
c  message depending on the array size
        if (no_of_xslices(myid) .eq. iyloc) then
          call MPI_RECV(s_x(1,1,1), 1, LOCAL_LARGE_3D_XSLICES, 0, 
     &			myid, MPI_COMM_WORLD, status, ierr)
        else
          call MPI_RECV(s_x(1,1,1), 1, LOCAL_SMALL_3D_XSLICES, 0, 
     &			myid, MPI_COMM_WORLD, status, ierr)
	end if	  
      end if
c
      end subroutine  distrib_x_3D_v2    

      
c--------------------------------------------------------------------------
c  Master distributes a 3D (ix.iy.iz) array of data in y-slice format 
c   to all yworkers
c--------------------------------------------------------------------------
c
      subroutine distrib_y_3D_v2(ix,iy,iz,ixloc,s,s_y)

      implicit none
c      
      integer :: ix, iy, iz, ixloc
      real   :: s(ix,iy,iz)
      real   :: s_y(ixloc,iy,iz)
      integer :: i, ierr, status(MPI_STATUS_SIZE)
c      
      if (Master) then
c The master sends one y-slice array to each of the yworkers using
c  one message, depending on the size.
        do i=1,NYworkers
	  if (no_of_yslices(i) .eq. ixloc) then
            call MPI_SEND(s(i,1,1), 1, GLOBAL_LARGE_3D_YSLICES, 
     &         owner_of_yslice(i), i,  MPI_COMM_WORLD, ierr)
          else
            call MPI_SEND(s(i,1,1), 1, GLOBAL_SMALL_3D_YSLICES, 
     &         owner_of_yslice(i), i,  MPI_COMM_WORLD, ierr)
          end if
	end do
      else if (YWorker) then
c The yworkers receive their y-slice array at once depending on
c  the size
        if (no_of_yslices(myid) .eq. ixloc) then
          call MPI_RECV(s_y(1,1,1), 1, LOCAL_LARGE_3D_YSLICES, 0, 
     &			myid, MPI_COMM_WORLD, status, ierr)
	else
          call MPI_RECV(s_y(1,1,1), 1, LOCAL_SMALL_3D_YSLICES, 0, 
     &			myid, MPI_COMM_WORLD, status, ierr)
	end if
      end if
c	
      end subroutine distrib_y_3D_v2 

c---------------------------------------------------------------------------
c  Master distributes the x-boundary array BDx
c---------------------------------------------------------------------------
c      
      subroutine distrib_x_BDx_v2(iy,iz,iyloc,Ns,s,sloc)

      implicit none
c
      integer :: iy, iz, Ns, iyloc
      real :: s(iy,iz,Ns)
      real :: sloc(iyloc,iz,Ns)
      integer :: j, ierr, status(MPI_STATUS_SIZE)
c      
      if (Master) then
        do j=1,NXworkers
c The master sends the x-slice arrays to the xworkers using one
c  message for each xworker, depending on the array size
	  if (no_of_xslices(j) .eq. iyloc) then
            call MPI_SEND(s(j,1,1), 1, GLOBAL_LARGE_BDx_XSLICES, 
     &         owner_of_xslice(j), j,  MPI_COMM_WORLD, ierr)
          else
            call MPI_SEND(s(j,1,1), 1, GLOBAL_SMALL_BDx_XSLICES, 
     &         owner_of_xslice(j), j,  MPI_COMM_WORLD, ierr)
          end if
	end do
      else if (XWorker) then
c The xworkers receive the x-slice arrays from the master in one
c  message depending on the array size
        if (no_of_xslices(myid) .eq. iyloc) then
          call MPI_RECV(sloc(1,1,1), 1, LOCAL_LARGE_BDx_XSLICES, 0, 
     &			myid, MPI_COMM_WORLD, status, ierr)
        else
          call MPI_RECV(sloc(1,1,1), 1, LOCAL_SMALL_BDx_XSLICES, 0, 
     &			myid, MPI_COMM_WORLD, status, ierr)
	end if	  
      end if
c
      end subroutine  distrib_x_BDx_v2
      
c---------------------------------------------------------------------------
c  Master distributes the x-boundary array BDx2
c---------------------------------------------------------------------------
c      
      subroutine distrib_x_BD_v2(iy,iz,iyloc,Ns,s,sloc)

      implicit none
c
      integer :: iy, iz, Ns, iyloc
      real :: s(iy,iz,2,Ns)
      real :: sloc(iyloc,iz,2,Ns)
      integer :: j, ierr, status(MPI_STATUS_SIZE)
c      
      if (Master) then
        do j=1,NXworkers
c The master sends the x-slice arrays to the xworkers using one
c  message for each xworker, depending on the array size
	  if (no_of_xslices(j) .eq. iyloc) then
            call MPI_SEND(s(j,1,1,1), 1, GLOBAL_LARGE_BD_XSLICES, 
     &         owner_of_xslice(j), j,  MPI_COMM_WORLD, ierr)
          else
            call MPI_SEND(s(j,1,1,1), 1, GLOBAL_SMALL_BD_XSLICES, 
     &         owner_of_xslice(j), j,  MPI_COMM_WORLD, ierr)
          end if
	end do
      else if (XWorker) then
c The xworkers receive the x-slice arrays from the master in one
c  message depending on the array size
        if (no_of_xslices(myid) .eq. iyloc) then
          call MPI_RECV(sloc(1,1,1,1), 1, LOCAL_LARGE_BD_XSLICES, 0, 
     &			myid, MPI_COMM_WORLD, status, ierr)
        else
          call MPI_RECV(sloc(1,1,1,1), 1, LOCAL_SMALL_BD_XSLICES, 0, 
     &			myid, MPI_COMM_WORLD, status, ierr)
	end if	  
      end if
c
      end subroutine  distrib_x_BD_v2
      
c---------------------------------------------------------------------------
c  Master distributes the y-boundary array BDy
c---------------------------------------------------------------------------
c      
      subroutine distrib_y_BDy_v2(ix,iz,ixloc,Ns,s,sloc)

      implicit none
c
      integer :: ix, iz, Ns, ixloc
      real :: s(ix,iz,Ns)
      real :: sloc(ixloc,iz,Ns)
      integer :: i, ierr, status(MPI_STATUS_SIZE)
c      
      if (Master) then
c The master sends one y-slice array to each of the yworkers using
c  one message, depending on the size.
        do i=1,NYworkers
	  if (no_of_yslices(i) .eq. ixloc) then
            call MPI_SEND(s(i,1,1), 1, GLOBAL_LARGE_BDy_YSLICES, 
     &         owner_of_yslice(i), i,  MPI_COMM_WORLD, ierr)
          else
            call MPI_SEND(s(i,1,1), 1, GLOBAL_SMALL_BDy_YSLICES, 
     &         owner_of_yslice(i), i,  MPI_COMM_WORLD, ierr)
          end if
	end do
      else if (YWorker) then
c The yworkers receive their y-slice array at once depending on
c  the size
        if (no_of_yslices(myid) .eq. ixloc) then
          call MPI_RECV(sloc(1,1,1), 1, LOCAL_LARGE_BDy_YSLICES, 0, 
     &			myid, MPI_COMM_WORLD, status, ierr)
	else
          call MPI_RECV(sloc(1,1,1), 1, LOCAL_SMALL_BDy_YSLICES, 0, 
     &			myid, MPI_COMM_WORLD, status, ierr)
	end if
      end if

      end subroutine  distrib_y_BDy_v2

c---------------------------------------------------------------------------
c  Master distributes the y-boundary array BDy2
c---------------------------------------------------------------------------
c      
      subroutine distrib_y_BD_v2(ix,iz,ixloc,Ns,s,sloc)

      implicit none
c
      integer :: ix, iz, Ns, ixloc
      real :: s(ix,iz,2,Ns)
      real :: sloc(ixloc,iz,2,Ns)
      integer :: i, ierr, status(MPI_STATUS_SIZE)
c      
      if (Master) then
c The master sends one y-slice array to each of the yworkers using
c  one message, depending on the size.
        do i=1,NYworkers
	  if (no_of_yslices(i) .eq. ixloc) then
            call MPI_SEND(s(i,1,1,1), 1, GLOBAL_LARGE_BD_YSLICES, 
     &         owner_of_yslice(i), i,  MPI_COMM_WORLD, ierr)
          else
            call MPI_SEND(s(i,1,1,1), 1, GLOBAL_SMALL_BD_YSLICES, 
     &         owner_of_yslice(i), i,  MPI_COMM_WORLD, ierr)
          end if
	end do
      else if (YWorker) then
c The yworkers receive their y-slice array at once depending on
c  the size
        if (no_of_yslices(myid) .eq. ixloc) then
          call MPI_RECV(sloc(1,1,1,1), 1, LOCAL_LARGE_BD_YSLICES, 0, 
     &			myid, MPI_COMM_WORLD, status, ierr)
	else
          call MPI_RECV(sloc(1,1,1,1), 1, LOCAL_SMALL_BD_YSLICES, 0, 
     &			myid, MPI_COMM_WORLD, status, ierr)
	end if
      end if

      end subroutine  distrib_y_BD_v2

c *******************************************************************
c ***************** collective communication ************************
c *******************************************************************
      
c-----------------------------------------------------------------
c  Master distributes a 4D array of data in x-slice format to all workers
c  using scatterv
c-----------------------------------------------------------------
      subroutine distrib_x_4D_v3(ix,iy,iz,iyloc,Ns,s,s_x)
c      
      implicit none
c
      integer :: ix, iy, iz, iyloc, Ns
      real :: s(ix,iy,iz,Ns)
      real :: s_x(ix,iyloc,iz,Ns)
      integer :: i, ierr
      integer :: sendcounts(Nprocs), displs(Nprocs)
      integer :: recvcount

      sendcounts(1) = 0
      sendcounts(2:NXworkers+1) = 1
      sendcounts(NXworkers+2:Nprocs) = 0
      displs(1:2) = 0
      do i=3, NXworkers+1
	displs(i) = displs(i-1) + ix
      end do
      if (.not. XWorker) then
        recvcount = 0
      else
        recvcount = 1
      end if
c The master scatters the x-slices in one call if all have the
c same amount of slices
      if (no_of_xslices(NXworkers) .eq. iyloc) then
        call MPI_SCATTERV(s(1,1,1,1), sendcounts, displs, 
     &		GLOB_LARGE_4D_XSLICES, s_x(1,1,1,1),
     &		recvcount, LOCAL_LARGE_4D_XSLICES, 0,
     &		MPI_COMM_WORLD, ierr)
      else
c Otherwise we need 2 calls.
c First everybody gets the minimal size everybody has
        call MPI_SCATTERV(s(1,1,1,1), sendcounts, displs, 
     &		GLOB_SMALL_4D_XSLICES, s_x(1,1,1,1),
     &		recvcount, LOCAL_SMALL_4D_XSLICES, 0,
     &		MPI_COMM_WORLD, ierr)
	sendcounts((mod(iy,NXworkers)+2) : Nprocs) = 0
	displs((mod(iy,NXworkers)+2) : Nprocs) = 0
        if (MyId .gt. mod(iy,NXworkers)) then
          recvcount = 0
        end if
c Then the rest processors get their further slices (max. one more)
        call MPI_SCATTERV(s(1,NXworkers*(iyloc-1)+1,1,1), sendcounts,  
     &		displs, GLOB_4D_XSLICE, s_x(1,iyloc,1,1),
     &		recvcount, LOCAL_4D_XSLICE, 0,
     &		MPI_COMM_WORLD, ierr)
      end if

      end subroutine  distrib_x_4D_v3

      
c-----------------------------------------------------------------
c  Master distributes a 4D array of data in y-slice format to all workers
c  using scatterv
c-----------------------------------------------------------------
      subroutine distrib_y_4D_v3(ix,iy,iz,ixloc,Ns,s,s_y)
c      
      implicit none
c
      integer :: ix, iy, iz, ixloc, Ns
      real :: s(ix,iy,iz,Ns)
      real :: s_y(ixloc,iy,iz,Ns)
      integer :: i, ierr
      integer :: sendcounts(Nprocs), displs(Nprocs)
      integer :: recvcount

      sendcounts(1) = 0
      sendcounts(2:NYworkers+1) = 1
      sendcounts(NYworkers+2:Nprocs) = 0
      displs(1:2) = 0
      do i=3, NYworkers+1
	displs(i) = displs(i-1) + 1
      end do
      if (.not. YWorker) then
        recvcount = 0
      else
        recvcount = 1
      end if
c If everyone has the same number of slices, one scatter is 
c  needed only
      if (no_of_yslices(NYworkers) .eq. ixloc) then
        call MPI_SCATTERV(s(1,1,1,1), sendcounts, displs, 
     &		GLOB_LARGE_4D_YSLICES, s_y(1,1,1,1),
     &		recvcount, LOCAL_LARGE_4D_YSLICES, 0,
     &		MPI_COMM_WORLD, ierr)
      else
c The matrix has to be scattered in 2 calls
c First everybody gets the minimal size everybody has
        call MPI_SCATTERV(s(1,1,1,1), sendcounts, displs, 
     &		GLOB_SMALL_4D_YSLICES, s_y(1,1,1,1),
     &		recvcount, LOCAL_SMALL_4D_YSLICES, 0,
     &		MPI_COMM_WORLD, ierr)
	sendcounts((mod(ix,NYworkers)+2) : Nprocs) = 0
	displs((mod(ix,NYworkers)+2) : Nprocs) = 0
        if (MyId .gt. mod(ix,NYworkers)+1) then
          recvcount = 0
        end if
c then the rest processors get their further slices (max. one more)        
	call MPI_SCATTERV(s(NYworkers*(ixloc-1)+1,1,1,1), sendcounts,  
     &		displs, GLOB_4D_YSLICE, s_y(ixloc,1,1,1),
     &		recvcount, LOCAL_4D_YSLICE, 0,
     &		MPI_COMM_WORLD, ierr)
      end if

      end subroutine  distrib_y_4D_v3

c-----------------------------------------------------------------
c  Master gathers a 4D array of data in x-slice format from all workers
c  using gatherv
c-----------------------------------------------------------------
      subroutine gather_x_4D_v3(ix,iy,iz,iyloc,Ns,s,s_x)
c      
      implicit none
c
      integer :: ix, iy, iz, iyloc, Ns
      real :: s(ix,iy,iz,Ns)
      real :: s_x(ix,iyloc,iz,Ns)
      integer :: i, ierr
      integer :: recvcounts(Nprocs), displs(Nprocs)
      integer :: sendcount

      recvcounts(1) = 0
      recvcounts(2:NXworkers+1) = 1
      recvcounts(NXworkers+2:Nprocs) = 0
      displs(1:2) = 0
      do i=3, NXworkers+1
	displs(i) = displs(i-1) + ix
      end do
      displs(NXworkers+2:Nprocs) = 0
      if (.not. XWorker) then
        sendcount = 0
      else
        sendcount = 1
      end if
c If everyone has the same number of slices as the last processor, 
c  one gather is needed only
      if (no_of_xslices(NXworkers) .eq. iyloc) then
        call MPI_GATHERV(s_x(1,1,1,1), sendcount,
     &		LOCAL_LARGE_4D_XSLICES, s(1,1,1,1),
     &		recvcounts, displs, GLOB_LARGE_4D_XSLICES, 0,
     &		MPI_COMM_WORLD, ierr)
      else
c The matrix has to be gathered in 2 calls
c First everybody sends the minimal size everybody has
        call MPI_GATHERV(s_x(1,1,1,1), sendcount,
     &		LOCAL_SMALL_4D_XSLICES, s(1,1,1,1),
     &		recvcounts, displs, GLOB_SMALL_4D_XSLICES, 0,
     &		MPI_COMM_WORLD, ierr)
	recvcounts((mod(iy,NXworkers)+2) : Nprocs) = 0
	displs((mod(iy,NXworkers)+2) : Nprocs) = 0
        if (MyId .gt. mod(iy,NXworkers)) then
          sendcount = 0
        end if
c then the rest processors send their further slices (max. one more)        
        call MPI_GATHERV(s_x(1,iyloc,1,1), sendcount,
     &		LOCAL_4D_XSLICE, s(1,NXworkers*(iyloc-1)+1,1,1),
     &		recvcounts, displs, GLOB_4D_XSLICE, 0,
     &		MPI_COMM_WORLD, ierr)
      end if

      end subroutine  gather_x_4D_v3

c-----------------------------------------------------------------
c  Master gathers a 4D array of data in y-slice format from all workers
c  using gatherv
c-----------------------------------------------------------------
      subroutine gather_y_4D_v3(ix,iy,iz,ixloc,Ns,s,s_y)
c      
      implicit none
c
      integer :: ix, iy, iz, ixloc, Ns
      real :: s(ix,iy,iz,Ns)
      real :: s_y(ixloc,iy,iz,Ns)
      integer :: i, ierr
      integer :: recvcounts(Nprocs), displs(Nprocs)
      integer :: sendcount

      recvcounts(1) = 0
      recvcounts(2:NYworkers+1) = 1
      recvcounts(NYworkers+2:Nprocs) = 0
      displs(1:2) = 0
      do i=3, NYworkers+1
	displs(i) = displs(i-1) + 1
      end do
      displs(NYworkers+2:Nprocs) = 0
      if (.not. YWorker) then
        sendcount = 0
      else
        sendcount = 1
      end if
c If everyone has the same number of slices, one gather is 
c  needed only
      if (no_of_yslices(NYworkers) .eq. ixloc) then
        call MPI_GATHERV(s_y(1,1,1,1), sendcount,
     &		LOCAL_LARGE_4D_YSLICES, s(1,1,1,1),
     &		recvcounts, displs, GLOB_LARGE_4D_YSLICES, 0,
     &		MPI_COMM_WORLD, ierr)
      else
c The matrix has to be gathered in 2 calls
c First everybody sends the minimal size everybody has
        call MPI_GATHERV(s_y(1,1,1,1), sendcount,
     &		LOCAL_SMALL_4D_YSLICES, s(1,1,1,1),
     &		recvcounts, displs, GLOB_SMALL_4D_YSLICES, 0,
     &		MPI_COMM_WORLD, ierr)
	recvcounts((mod(ix,NYworkers)+2) : Nprocs) = 0
	displs((mod(ix,NYworkers)+2) : Nprocs) = 0
        if (MyId .gt. mod(ix,NYworkers)) then
          sendcount = 0
        end if
c then the rest processors send their further slices (max. one more)        
        call MPI_GATHERV(s_y(ixloc,1,1,1), sendcount,
     &		LOCAL_4D_YSLICE, s(NYworkers*(ixloc-1)+1,1,1,1),
     &		recvcounts, displs, GLOB_4D_YSLICE, 0,
     &		MPI_COMM_WORLD, ierr)
      end if

      end subroutine  gather_y_4D_v3

c---------------------------------------------------------
c      Slaves exchange data:  - old format is x-slices (for sending)
c                             - new format is y-slices (for receiving)
c	using alltoallv and copying - easier, just one call!
c---------------------------------------------------------
      subroutine shuffle_x2y_4D_v3(ix,iy,iz,ixloc,iyloc,Ns,s_x,s_y)  
c      
      implicit none
c
      integer :: ix, iy, iz, ixloc, iyloc, Ns
      real   :: s_x(ix,iyloc,iz,Ns)
      real   :: s_y(ixloc,iy,iz,Ns)
      real   :: sbuf(ixloc,iyloc,iz,Ns,max(NXworkers,NYworkers))
      real   :: rbuf(ixloc,iyloc,iz,Ns,max(NXworkers,NYworkers))
      integer :: i, j, ierr, tag
      integer :: status(MPI_STATUS_SIZE)
      integer :: sendcounts(max(NXworkers,NYworkers)) 
      integer :: sdispls(max(NXworkers,NYworkers))
      integer :: recvcounts(max(NXworkers,NYworkers))
      integer :: rdispls(max(NXworkers,NYworkers))
      integer :: participants

c      
      if (XWorker .or. YWorker) then
        participants = max(NXworkers,NYworkers)
        sdispls(1) = 0
        rdispls(1) = 0
        do i=2,participants
          sdispls(i) = (i-1)*ixloc*iyloc*iz*Ns
	  rdispls(i) = sdispls(i)
        end do

	if (XWorker) then
	  do i=1,NYworkers
	    sbuf(1:no_of_yslices(i),1:iyloc,1:iz,1:Ns,i) = 
     &		s_x(owned_yslices(i,1:no_of_yslices(i)),
     &				1:iyloc,1:iz,1:Ns)
          end do
	end if
	sendcounts(1:participants) = ixloc*iyloc*iz*Ns
        recvcounts(1:participants) = ixloc*iyloc*iz*Ns

	call MPI_ALLTOALLV(sbuf, sendcounts, sdispls, 
     &		MPI_REAL, rbuf, recvcounts, rdispls,
     &		MPI_REAL, MPI_COMM_WORKERS, ierr)

	if (YWorker) then
	  do i=1,NXworkers
	    s_y(1:ixloc,
     &	   	owned_xslices(i,1:no_of_xslices(i)),1:iz,1:Ns) = 
     &	   	rbuf(1:ixloc,1:no_of_xslices(i),1:iz,1:Ns,i)
          end do
	end if
            
      end if
c	
      end subroutine shuffle_x2y_4D_v3

c---------------------------------------------------------
c      Slaves exchange data:  - old format is y-slices (for sending)
c                             - new format is x-slices (for receiving)
c	using alltoallv and copying - easier, just one call!
c	garbage sent
c---------------------------------------------------------
      subroutine shuffle_y2x_4D_v3(ix,iy,iz,ixloc,iyloc,Ns,s_x,s_y)  
c      
      implicit none
c
      integer :: ix, iy, iz, ixloc, iyloc, Ns
      real   :: s_x(ix,iyloc,iz,Ns)
      real   :: s_y(ixloc,iy,iz,Ns)
      real  :: sbuf(ixloc,iyloc,iz,Ns,max(NXworkers,NYworkers))
      real   :: rbuf(ixloc,iyloc,iz,Ns,max(NXworkers,NYworkers))
      integer :: i, j, ierr, tag
      integer :: status(MPI_STATUS_SIZE)
      integer :: sendcounts(max(NXworkers,NYworkers)) 
      integer :: sdispls(max(NXworkers,NYworkers))
      integer :: recvcounts(max(NXworkers,NYworkers))
      integer :: rdispls(max(NXworkers,NYworkers))
      integer :: participants

c      
      if (XWorker .or. YWorker) then
        participants = max(NXworkers,NYworkers)
        sdispls(1) = 0
        rdispls(1) = 0
        do i=2,participants
          sdispls(i) = (i-1)*ixloc*iyloc*iz*Ns
	  rdispls(i) = sdispls(i)
        end do
	
	if (YWorker) then
	  do i=1,NXworkers
	    sbuf(1:ixloc,1:no_of_xslices(i),1:iz,1:Ns,i) = 
     &		s_y(1:ixloc,owned_xslices(i,1:no_of_xslices(i)),
     &				1:iz,1:Ns)
          end do
	end if
	sendcounts(1:participants) = ixloc*iyloc*iz*Ns
        recvcounts(1:participants) = ixloc*iyloc*iz*Ns
        
	call MPI_ALLTOALLV(sbuf, sendcounts, sdispls, 
     &		MPI_REAL, rbuf, recvcounts, rdispls,
     &		MPI_REAL, MPI_COMM_WORKERS, ierr)

	if (XWorker) then
	  do i=1,NYworkers
	    s_x(owned_yslices(i,1:no_of_yslices(i)),
     &	   	1:iyloc,1:iz,1:Ns) = 
     &	   	rbuf(1:no_of_yslices(i),1:iyloc,1:iz,1:Ns,i)
          end do
	end if
            
      end if
c	
        end subroutine shuffle_y2x_4D_v3

c-----------------------------------------------------------------
c  Master distributes a 2D array of data in x-slice format to all workers
c  using scatterv
c-----------------------------------------------------------------
      subroutine distrib_x_2D_v3(ix,iy,iyloc,s,s_x)
c      
      implicit none
c
      integer :: ix, iy, iyloc
      real :: s(ix,iy)
      real :: s_x(ix,iyloc)
      integer :: i, ierr
      integer :: sendcounts(Nprocs), displs(Nprocs)
      integer :: recvcount

      sendcounts(1) = 0
      sendcounts(2:NXworkers+1) = 1
      sendcounts(NXworkers+2:Nprocs) = 0
      displs(1:2) = 0
      do i=3, NXworkers+1
	displs(i) = displs(i-1) + ix
      end do
      if (.not. XWorker) then
        recvcount = 0
      else
        recvcount = 1
      end if
c The master scatters the x-slices in one call if all have the
c same amount of slices
      if (no_of_xslices(NXworkers) .eq. iyloc) then
        call MPI_SCATTERV(s(1,1), sendcounts, displs, 
     &		GLOB_LARGE_2D_XSLICES, s_x(1,1),
     &		recvcount, LOCAL_LARGE_2D_XSLICES, 0,
     &		MPI_COMM_WORLD, ierr)
      else
c Otherwise we need 2 calls.
c First everybody gets the minimal size everybody has
        call MPI_SCATTERV(s(1,1), sendcounts, displs, 
     &		GLOB_SMALL_2D_XSLICES, s_x(1,1),
     &		recvcount, LOCAL_SMALL_2D_XSLICES, 0,
     &		MPI_COMM_WORLD, ierr)
	sendcounts(2:mod(iy,NXworkers)+1) = ix
	sendcounts((mod(iy,NXworkers)+2) : Nprocs) = 0
	displs((mod(iy,NXworkers)+2) : Nprocs) = 0
        if (MyId .gt. mod(iy,NXworkers)) then
          recvcount = 0
	else if (.not. Master) then
	  recvcount = ix
        end if
c Then the rest processors get their further slices (max. one more)
        call MPI_SCATTERV(s(1,NXworkers*(iyloc-1)+1), sendcounts,  
     &		displs, MPI_REAL, s_x(1,iyloc),
     &		recvcount, MPI_REAL, 0,
     &		MPI_COMM_WORLD, ierr)
      end if

      end subroutine  distrib_x_2D_v3

      
c-----------------------------------------------------------------
c  Master distributes a 2D array of data in y-slice format to all workers
c  using scatterv
c-----------------------------------------------------------------
      subroutine distrib_y_2D_v3(ix,iy,ixloc,s,s_y)
c      
      implicit none
c
      integer :: ix, iy, ixloc
      real :: s(ix,iy)
      real :: s_y(ixloc,iy)
      integer :: i, ierr
      integer :: sendcounts(Nprocs), displs(Nprocs)
      integer :: recvcount

      sendcounts(1) = 0
      sendcounts(2:NYworkers+1) = 1
      sendcounts(NYworkers+2:Nprocs) = 0
      displs(1:2) = 0
      do i=3, NYworkers+1
	displs(i) = displs(i-1) + 1
      end do
      if (.not. YWorker) then
        recvcount = 0
      else
        recvcount = 1
      end if
c If everyone has the same number of slices, one scatter is 
c  needed only
      if (no_of_yslices(NYworkers) .eq. ixloc) then
        call MPI_SCATTERV(s(1,1), sendcounts, displs, 
     &		GLOB_LARGE_2D_YSLICES, s_y(1,1),
     &		recvcount, LOCAL_LARGE_2D_YSLICES, 0,
     &		MPI_COMM_WORLD, ierr)
      else
c The matrix has to be scattered in 2 calls
c First everybody gets the minimal size everybody has
        call MPI_SCATTERV(s(1,1), sendcounts, displs, 
     &		GLOB_SMALL_2D_YSLICES, s_y(1,1),
     &		recvcount, LOCAL_SMALL_2D_YSLICES, 0,
     &		MPI_COMM_WORLD, ierr)
	sendcounts((mod(ix,NYworkers)+2) : Nprocs) = 0
	displs((mod(ix,NYworkers)+2) : Nprocs) = 0
        if (MyId .gt. mod(ix,NYworkers)) then
          recvcount = 0
        end if
c then the rest processors get their further slices (max. one more)        
	call MPI_SCATTERV(s(NYworkers*(ixloc-1)+1,1), sendcounts,  
     &		displs, GLOB_2D_YSLICE, s_y(ixloc,1),
     &		recvcount, LOCAL_2D_YSLICE, 0,
     &		MPI_COMM_WORLD, ierr)
      end if

      end subroutine  distrib_y_2D_v3

c-----------------------------------------------------------------
c  Master distributes a BDz array of data in x-slice format to all workers
c  using scatterv
c-----------------------------------------------------------------
      subroutine distrib_x_BDz_v3(ix,iy,iyloc,Ns,s,s_x)
c      
      implicit none
c
      integer :: ix, iy, iyloc, Ns
      real :: s(ix,iy,Ns)
      real :: s_x(ix,iyloc,Ns)
      integer :: i, ierr
      integer :: sendcounts(Nprocs), displs(Nprocs)
      integer :: recvcount

      sendcounts(1) = 0
      sendcounts(2:NXworkers+1) = 1
      sendcounts(NXworkers+2:Nprocs) = 0
      displs(1:2) = 0
      do i=3, NXworkers+1
	displs(i) = displs(i-1) + ix
      end do
      if (.not. XWorker) then
        recvcount = 0
      else
        recvcount = 1
      end if
c The master scatters the x-slices in one call if all have the
c same amount of slices
      if (no_of_xslices(NXworkers) .eq. iyloc) then
        call MPI_SCATTERV(s(1,1,1), sendcounts, displs, 
     &		GLOB_LARGE_BDz_XSLICES, s_x(1,1,1),
     &		recvcount, LOCAL_LARGE_BDz_XSLICES, 0,
     &		MPI_COMM_WORLD, ierr)
      else
c Otherwise we need 2 calls.
c First everybody gets the minimal size everybody has
        call MPI_SCATTERV(s(1,1,1), sendcounts, displs, 
     &		GLOB_SMALL_BDz_XSLICES, s_x(1,1,1),
     &		recvcount, LOCAL_SMALL_BDz_XSLICES, 0,
     &		MPI_COMM_WORLD, ierr)
	sendcounts((mod(iy,NXworkers)+2) : Nprocs) = 0
	displs((mod(iy,NXworkers)+2) : Nprocs) = 0
        if (MyId .gt. mod(iy,NXworkers)) then
          recvcount = 0
        end if
c Then the rest processors get their further slices (max. one more)
        call MPI_SCATTERV(s(1,NXworkers*(iyloc-1)+1,1), sendcounts,  
     &		displs, GLOB_BDz_XSLICE, s_x(1,iyloc,1),
     &		recvcount, LOCAL_BDz_XSLICE, 0,
     &		MPI_COMM_WORLD, ierr)
      end if

      end subroutine  distrib_x_BDz_v3

      
c-----------------------------------------------------------------
c  Master distributes a BDz array of data in y-slice format to all workers
c  using scatterv
c-----------------------------------------------------------------
      subroutine distrib_y_BDz_v3(ix,iy,ixloc,Ns,s,s_y)
c      
      implicit none
c
      integer :: ix, iy, ixloc, Ns
      real :: s(ix,iy,Ns)
      real :: s_y(ixloc,iy,Ns)
      integer :: i, ierr
      integer :: sendcounts(Nprocs), displs(Nprocs)
      integer :: recvcount

      sendcounts(1) = 0
      sendcounts(2:NYworkers+1) = 1
      sendcounts(NYworkers+2:Nprocs) = 0
      displs(1:2) = 0
      do i=3, NYworkers+1
	displs(i) = displs(i-1) + 1
      end do
      if (.not. YWorker) then
        recvcount = 0
      else
        recvcount = 1
      end if
c If everyone has the same number of slices, one scatter is 
c  needed only
      if (no_of_yslices(NYworkers) .eq. ixloc) then
        call MPI_SCATTERV(s(1,1,1), sendcounts, displs, 
     &		GLOB_LARGE_BDz_YSLICES, s_y(1,1,1),
     &		recvcount, LOCAL_LARGE_BDz_YSLICES, 0,
     &		MPI_COMM_WORLD, ierr)
      else
c The matrix has to be scattered in 2 calls
c First everybody gets the minimal size everybody has
        call MPI_SCATTERV(s(1,1,1), sendcounts, displs, 
     &		GLOB_SMALL_BDz_YSLICES, s_y(1,1,1),
     &		recvcount, LOCAL_SMALL_BDz_YSLICES, 0,
     &		MPI_COMM_WORLD, ierr)
	sendcounts((mod(ix,NYworkers)+2) : Nprocs) = 0
	displs((mod(ix,NYworkers)+2) : Nprocs) = 0
        if (MyId .gt. mod(ix,NYworkers)) then
          recvcount = 0
        end if
c then the rest processors get their further slices (max. one more)        
	call MPI_SCATTERV(s(NYworkers*(ixloc-1)+1,1,1), sendcounts,  
     &		displs, GLOB_BDz_YSLICE, s_y(ixloc,1,1),
     &		recvcount, LOCAL_BDz_YSLICE, 0,
     &		MPI_COMM_WORLD, ierr)
      end if

      end subroutine  distrib_y_BDz_v3

c-----------------------------------------------------------------
c  Master distributes a 3D array of data in x-slice format to all workers
c  using scatterv
c-----------------------------------------------------------------
      subroutine distrib_x_3D_v3(ix,iy,iz,iyloc,s,s_x)
c      
      implicit none
c
      integer :: ix, iy, iz, iyloc
      real :: s(ix,iy,iz)
      real :: s_x(ix,iyloc,iz)
      integer :: i, ierr
      integer :: sendcounts(Nprocs), displs(Nprocs)
      integer :: recvcount

      sendcounts(1) = 0
      sendcounts(2:NXworkers+1) = 1
      sendcounts(NXworkers+2:Nprocs) = 0
      displs(1:2) = 0
      do i=3, NXworkers+1
	displs(i) = displs(i-1) + ix
      end do
      if (.not. XWorker) then
        recvcount = 0
      else
        recvcount = 1
      end if
c The master scatters the x-slices in one call if all have the
c same amount of slices
      if (no_of_xslices(NXworkers) .eq. iyloc) then
        call MPI_SCATTERV(s(1,1,1), sendcounts, displs, 
     &		GLOB_LARGE_3D_XSLICES, s_x(1,1,1),
     &		recvcount, LOCAL_LARGE_3D_XSLICES, 0,
     &		MPI_COMM_WORLD, ierr)
      else
c Otherwise we need 2 calls.
c First everybody gets the minimal size everybody has
        call MPI_SCATTERV(s(1,1,1), sendcounts, displs, 
     &		GLOB_SMALL_3D_XSLICES, s_x(1,1,1),
     &		recvcount, LOCAL_SMALL_3D_XSLICES, 0,
     &		MPI_COMM_WORLD, ierr)
	sendcounts((mod(iy,NXworkers)+2) : Nprocs) = 0
	displs((mod(iy,NXworkers)+2) : Nprocs) = 0
        if (MyId .gt. mod(iy,NXworkers)) then
          recvcount = 0
        end if
c Then the rest processors get their further slices (max. one more)
        call MPI_SCATTERV(s(1,NXworkers*(iyloc-1)+1,1), sendcounts,  
     &		displs, GLOB_3D_XSLICE, s_x(1,iyloc,1),
     &		recvcount, LOCAL_3D_XSLICE, 0,
     &		MPI_COMM_WORLD, ierr)
      end if

      end subroutine  distrib_x_3D_v3

      
c-----------------------------------------------------------------
c  Master distributes a 3D array of data in y-slice format to all workers
c  using scatterv
c-----------------------------------------------------------------
      subroutine distrib_y_3D_v3(ix,iy,iz,ixloc,s,s_y)
c      
      implicit none
c
      integer :: ix, iy, iz, ixloc
      real :: s(ix,iy,iz)
      real :: s_y(ixloc,iy,iz)
      integer :: i, ierr
      integer :: sendcounts(Nprocs), displs(Nprocs)
      integer :: recvcount

      sendcounts(1) = 0
      sendcounts(2:NYworkers+1) = 1
      sendcounts(NYworkers+2:Nprocs) = 0
      displs(1:2) = 0
      do i=3, NYworkers+1
	displs(i) = displs(i-1) + 1
      end do
      if (.not. YWorker) then
        recvcount = 0
      else
        recvcount = 1
      end if
c If everyone has the same number of slices, one scatter is 
c  needed only
      if (no_of_yslices(NYworkers) .eq. ixloc) then
        call MPI_SCATTERV(s(1,1,1), sendcounts, displs, 
     &		GLOB_LARGE_3D_YSLICES, s_y(1,1,1),
     &		recvcount, LOCAL_LARGE_3D_YSLICES, 0,
     &		MPI_COMM_WORLD, ierr)
      else
c The matrix has to be scattered in 2 calls
c First everybody gets the minimal size everybody has
        call MPI_SCATTERV(s(1,1,1), sendcounts, displs, 
     &		GLOB_SMALL_3D_YSLICES, s_y(1,1,1),
     &		recvcount, LOCAL_SMALL_3D_YSLICES, 0,
     &		MPI_COMM_WORLD, ierr)
	sendcounts((mod(ix,NYworkers)+2) : Nprocs) = 0
	displs((mod(ix,NYworkers)+2) : Nprocs) = 0
        if (MyId .gt. mod(ix,NYworkers)) then
          recvcount = 0
        end if
c then the rest processors get their further slices (max. one more)        
	call MPI_SCATTERV(s(NYworkers*(ixloc-1)+1,1,1), sendcounts,  
     &		displs, GLOB_3D_YSLICE, s_y(ixloc,1,1),
     &		recvcount, LOCAL_3D_YSLICE, 0,
     &		MPI_COMM_WORLD, ierr)
      end if

      end subroutine  distrib_y_3D_v3

c-----------------------------------------------------------------
c  Master distributes a BDx array of data in x-slice format to all workers
c  using scatterv
c-----------------------------------------------------------------
      subroutine distrib_x_BDx_v3(iy,iz,iyloc,Ns,s,s_x)
c      
      implicit none
c
      integer :: iy, iz, iyloc, Ns
      real :: s(iy,iz,Ns)
      real :: s_x(iyloc,iz,Ns)
      integer :: i, ierr
      integer :: sendcounts(Nprocs), displs(Nprocs)
      integer :: recvcount

      sendcounts(1) = 0
      sendcounts(2:NXworkers+1) = 1
      sendcounts(NXworkers+2:Nprocs) = 0
      displs(1:2) = 0
      do i=3, NXworkers+1
	displs(i) = displs(i-1) + 1
      end do
      if (.not. XWorker) then
        recvcount = 0
      else
        recvcount = 1
      end if
c The master scatters the x-slices in one call if all have the
c same amount of slices
      if (no_of_xslices(NXworkers) .eq. iyloc) then
        call MPI_SCATTERV(s(1,1,1), sendcounts, displs, 
     &		GLOB_LARGE_BDx_XSLICES, s_x(1,1,1),
     &		recvcount, LOCAL_LARGE_BDx_XSLICES, 0,
     &		MPI_COMM_WORLD, ierr)
      else
c Otherwise we need 2 calls.
c First everybody gets the minimal size everybody has
        call MPI_SCATTERV(s(1,1,1), sendcounts, displs, 
     &		GLOB_SMALL_BDx_XSLICES, s_x(1,1,1),
     &		recvcount, LOCAL_SMALL_BDx_XSLICES, 0,
     &		MPI_COMM_WORLD, ierr)
	sendcounts((mod(iy,NXworkers)+2) : Nprocs) = 0
	displs((mod(iy,NXworkers)+2) : Nprocs) = 0
        if (MyId .gt. mod(iy,NXworkers)) then
          recvcount = 0
        end if
c Then the rest processors get their further slices (max. one more)
        call MPI_SCATTERV(s(NXworkers*(iyloc-1)+1,1,1), sendcounts,  
     &		displs, GLOB_BDx_XSLICE, s_x(iyloc,1,1),
     &		recvcount, LOCAL_BDx_XSLICE, 0,
     &		MPI_COMM_WORLD, ierr)
      end if

      end subroutine  distrib_x_BDx_v3

      
c-----------------------------------------------------------------
c  Master distributes a BDx2 array of data in x-slice format to all workers
c  using scatterv
c-----------------------------------------------------------------
      subroutine distrib_x_BD_v3(iy,iz,iyloc,Ns,s,s_x)
c      
      implicit none
c
      integer :: iy, iz, iyloc, Ns
      real :: s(iy,iz,2,Ns)
      real :: s_x(iyloc,iz,2,Ns)
      integer :: i, ierr
      integer :: sendcounts(Nprocs), displs(Nprocs)
      integer :: recvcount

      sendcounts(1) = 0
      sendcounts(2:NXworkers+1) = 1
      sendcounts(NXworkers+2:Nprocs) = 0
      displs(1:2) = 0
      do i=3, NXworkers+1
	displs(i) = displs(i-1) + 1
      end do
      if (.not. XWorker) then
        recvcount = 0
      else
        recvcount = 1
      end if
c The master scatters the x-slices in one call if all have the
c same amount of slices
      if (no_of_xslices(NXworkers) .eq. iyloc) then
        call MPI_SCATTERV(s(1,1,1,1), sendcounts, displs, 
     &		GLOB_LARGE_BD_XSLICES, s_x(1,1,1,1),
     &		recvcount, LOCAL_LARGE_BD_XSLICES, 0,
     &		MPI_COMM_WORLD, ierr)
      else
c Otherwise we need 2 calls.
c First everybody gets the minimal size everybody has
        call MPI_SCATTERV(s(1,1,1,1), sendcounts, displs, 
     &		GLOB_SMALL_BD_XSLICES, s_x(1,1,1,1),
     &		recvcount, LOCAL_SMALL_BD_XSLICES, 0,
     &		MPI_COMM_WORLD, ierr)
	sendcounts((mod(iy,NXworkers)+2) : Nprocs) = 0
	displs((mod(iy,NXworkers)+2) : Nprocs) = 0
        if (MyId .gt. mod(iy,NXworkers)) then
          recvcount = 0
        end if
c Then the rest processors get their further slices (max. one more)
        call MPI_SCATTERV(s(NXworkers*(iyloc-1)+1,1,1,1), sendcounts,  
     &		displs, GLOB_BD_XSLICE, s_x(iyloc,1,1,1),
     &		recvcount, LOCAL_BD_XSLICE, 0,
     &		MPI_COMM_WORLD, ierr)
      end if

      end subroutine  distrib_x_BD_v3

      
c-----------------------------------------------------------------
c  Master distributes a BDy array of data in y-slice format to all workers
c  using scatterv
c-----------------------------------------------------------------
      subroutine distrib_y_BDy_v3(ix,iz,ixloc,Ns,s,s_y)
c      
      implicit none
c
      integer :: ix, iz, ixloc, Ns
      real :: s(ix,iz,Ns)
      real :: s_y(ixloc,iz,Ns)
      integer :: i, ierr
      integer :: sendcounts(Nprocs), displs(Nprocs)
      integer :: recvcount

      sendcounts(1) = 0
      sendcounts(2:NYworkers+1) = 1
      sendcounts(NYworkers+2:Nprocs) = 0
      displs(1:2) = 0
      do i=3, NYworkers+1
	displs(i) = displs(i-1) + 1
      end do
      if (.not. YWorker) then
        recvcount = 0
      else
        recvcount = 1
      end if
c If everyone has the same number of slices, one scatter is 
c  needed only
      if (no_of_yslices(NYworkers) .eq. ixloc) then
        call MPI_SCATTERV(s(1,1,1), sendcounts, displs, 
     &		GLOB_LARGE_BDy_YSLICES, s_y(1,1,1),
     &		recvcount, LOCAL_LARGE_BDy_YSLICES, 0,
     &		MPI_COMM_WORLD, ierr)
      else
c The matrix has to be scattered in 2 calls
c First everybody gets the minimal size everybody has
        call MPI_SCATTERV(s(1,1,1), sendcounts, displs, 
     &		GLOB_SMALL_BDy_YSLICES, s_y(1,1,1),
     &		recvcount, LOCAL_SMALL_BDy_YSLICES, 0,
     &		MPI_COMM_WORLD, ierr)
	sendcounts((mod(ix,NYworkers)+2) : Nprocs) = 0
	displs((mod(ix,NYworkers)+2) : Nprocs) = 0
        if (MyId .gt. mod(ix,NYworkers)) then
          recvcount = 0
        end if
c then the rest processors get their further slices (max. one more)        
	call MPI_SCATTERV(s(NYworkers*(ixloc-1)+1,1,1), sendcounts,  
     &		displs, GLOB_BDy_YSLICE, s_y(ixloc,1,1),
     &		recvcount, LOCAL_BDy_YSLICE, 0,
     &		MPI_COMM_WORLD, ierr)
      end if

      end subroutine  distrib_y_BDy_v3

c-----------------------------------------------------------------
c  Master distributes a BDy2 array of data in y-slice format to all workers
c  using scatterv
c-----------------------------------------------------------------
      subroutine distrib_y_BD_v3(ix,iz,ixloc,Ns,s,s_y)
c      
      implicit none
c
      integer :: ix, iz, ixloc, Ns
      real :: s(ix,iz,2,Ns)
      real :: s_y(ixloc,iz,2,Ns)
      integer :: i, ierr
      integer :: sendcounts(Nprocs), displs(Nprocs)
      integer :: recvcount

      sendcounts(1) = 0
      sendcounts(2:NYworkers+1) = 1
      sendcounts(NYworkers+2:Nprocs) = 0
      displs(1:2) = 0
      do i=3, NYworkers+1
	displs(i) = displs(i-1) + 1
      end do
      if (.not. YWorker) then
        recvcount = 0
      else
        recvcount = 1
      end if
c If everyone has the same number of slices, one scatter is 
c  needed only
      if (no_of_yslices(NYworkers) .eq. ixloc) then
        call MPI_SCATTERV(s(1,1,1,1), sendcounts, displs, 
     &		GLOB_LARGE_BD_YSLICES, s_y(1,1,1,1),
     &		recvcount, LOCAL_LARGE_BD_YSLICES, 0,
     &		MPI_COMM_WORLD, ierr)
      else
c The matrix has to be scattered in 2 calls
c First everybody gets the minimal size everybody has
        call MPI_SCATTERV(s(1,1,1,1), sendcounts, displs, 
     &		GLOB_SMALL_BD_YSLICES, s_y(1,1,1,1),
     &		recvcount, LOCAL_SMALL_BD_YSLICES, 0,
     &		MPI_COMM_WORLD, ierr)
	sendcounts((mod(ix,NYworkers)+2) : Nprocs) = 0
	displs((mod(ix,NYworkers)+2) : Nprocs) = 0
        if (MyId .gt. mod(ix,NYworkers)) then
          recvcount = 0
        end if
c then the rest processors get their further slices (max. one more)        
	call MPI_SCATTERV(s(NYworkers*(ixloc-1)+1,1,1,1), sendcounts,  
     &		displs, GLOB_BD_YSLICE, s_y(ixloc,1,1,1),
     &		recvcount, LOCAL_BD_YSLICE, 0,
     &		MPI_COMM_WORLD, ierr)
      end if

      end subroutine  distrib_y_BD_v3

c *******************************************************************
c ***************** a few copy versions *****************************
c *******************************************************************

c-----------------------------------------------------------------
c  Master distributes a 4D array of data in x-slice format to all workers
c  using scatter
c using copying
c-----------------------------------------------------------------
      subroutine distrib_x_4D_v4(ix,iy,iz,iyloc,Ns,s,s_x)
c      
      implicit none
c
      integer :: ix, iy, iz, iyloc, Ns
      real :: s(ix,iy,iz,Ns)
      real :: s_x(ix,iyloc,iz,Ns)
      integer :: i, ierr
      integer :: sendcount, recvcount
      real,pointer,dimension(:,:,:,:,:) :: buf

      sendcount = ix*iyloc*iz*Ns
	print*, 'start scat4'
      if (MyId .eq. 0) then
	allocate( buf(ix,iyloc, iz, Ns, NXworkers+1), STAT=ierr)
	print*, 'init buffer'
        do i = 1,iy
	  buf(1:ix, local_xslice_id(i), 1:iz, 1:Ns, owner_of_xslice(i)+1) =
     &		s(1:ix, i, 1:iz, 1:Ns)
        end do
      else
	allocate( buf(ix,iyloc, iz, Ns, 1), STAT=ierr)
        recvcount = ix*iyloc*iz*Ns
      end if
      print*, 'call scatter'
      call MPI_SCATTER(buf(1,1,1,1,1), sendcount, MPI_REAL, buf, 
     &		recvcount, MPI_REAL, 0, MPI_COMM_WORLD, ierr)
      if (myid .ne. 0) then
        s_x(:,:,:,:) = buf(:,:,:,:,1)
      end if
      end subroutine  distrib_x_4D_v4

c---------------------------------------------------------
c      Slaves exchange data:  - old format is x-slices (for sending)
c                             - new format is y-slices (for receiving)
c	using alltoallw
c    commented out because not implemented on the SUN ENTERPRISE...
c---------------------------------------------------------
      subroutine shuffle_x2y_4D_v4(ix,iy,iz,ixloc,iyloc,Ns,s_x,s_y)  
c      
      implicit none
c
      integer :: ix, iy, iz, ixloc, iyloc, Ns
      real   :: s_x(ix,iyloc,iz,Ns)
      real   :: s_y(ixloc,iy,iz,Ns)
      integer :: i, j, ierr, tag
      integer :: status(MPI_STATUS_SIZE)
      integer :: sendcounts(max(NXworkers,NYworkers))
      integer :: sdispls(max(NXworkers,NYworkers))
      integer :: recvcounts(max(NXworkers, NYworkers))
      integer :: rdispls(max(NXworkers,NYworkers))
      integer :: sendtypes(max(NXworkers, NYworkers))
      integer :: recvtypes(max(NXworkers, NYworkers))
      integer :: lb, sizeofreal

      if (XWorker .or. YWorker) then
c        call MPI_TYPE_GET_EXTENT(MPI_REAL, lb,sizeofreal,ierr)
        sendcounts(1:NXworkers) = 1
        recvcounts(1:NYworkers) = 1
        sdispls(1) = 0
        rdispls(1) = 0
	sendtypes(1) = MPI_REAL
	recvtypes(1) = MPI_REAL
        do i=2,max(NXworkers,NYworkers)
          sdispls(i) = sdispls(i-1) + sizeofreal
	  rdispls(i) = rdispls(i-1) + ixloc*sizeofreal
	  if (no_of_xslices(myid) .eq. iyloc) then
	    if (no_of_yslices(i-1) .eq. ixloc) then
	      sendtypes(i) = LOCAL_4D_LXLY_XCOLS
	      recvtypes(i) = LOCAL_4D_LXLY_XCOLS
	    else
	      sendtypes(i) = LOCAL_4D_LXSY_XCOLS
	      recvtypes(i) = LOCAL_4D_LXSY_XCOLS
	    end if
	  else
	    if (no_of_yslices(i-1) .eq. ixloc) then
	      sendtypes(i) = LOCAL_4D_SXLY_XCOLS
	      recvtypes(i) = LOCAL_4D_SXLY_XCOLS
	    else
	      sendtypes(i) = LOCAL_4D_SXSY_XCOLS
	      recvtypes(i) = LOCAL_4D_SXSY_XCOLS
	    end if
	  end if
        end do
	
c          call MPI_ALLTOALLW(s_x, sendcounts, sdispls, 
c     &		sendtypes, s_y, recvcounts, rdispls,
c     &		recvtypes, MPI_COMM_WORKERS, ierr)
            
      end if
c	
      end subroutine shuffle_x2y_4D_v4      
	
	
      end module XYCommunicationLibrary
      

c---------------------------------------------------------------------------------------
c
c Three versions of the following functions are defined in this file:
c
c	distrib_h_4D  	Distribution of the concentration array in
c			h-slices from the master to the workers
c	distrib_v_4D 	Distribution of the concentration array in
c			v-columns from the master to the workers
c	gather_h_4D	Gathering of the concentration array from
c			h-slices from the workers at the master
c	gather_v_4D	Gathering of the concentration array from
c			v-columns from the workers at the master
c	shuffle_h2v_4D	Shuffling of the concentration array
c			from h-slices to v-columns at each worker
c	shuffle_v2h_4D	Shuffling of the concentration array
c			from v-columns to h-slices at each worker
c	distrib_h_2D	Distribution of the 2D (x,y) arrays in
c			h-slices from the master to the workers
c	distrib_v_2D	Distribution of the 2D (x,y) arrays in
c			v-columns from the master to the workers
c	distrib_v_BDz	Distribution of the BDz (x,y,Ns) arrays in
c			v-columns from the master to the workers
c	distrib_h_3D	Distribution of the 3D (x,y,z) arrays in
c			h-slices from the master to the workers
c	distrib_v_3D	Distribution of the 3D (x,y,z) arrays in
c			v-columns from the master to the workers
c	distrib_h_BDx	Distribution of the x-boundary data arrays in
c			h-slices from the master to the workers
c	distrib_h_BDy	Distribution of the y-boundary data arrays in
c			h-slices from the master to the workers
c
c	distrib_xh_BD	Distribution of the y-boundary data arrays in
c			h-slices from the master to the workers
c			Stem-III specific
c	distrib_yh_BD	Distribution of the x-boundary data arrays in
c			h-slices from the master to the workers
c			Stem-III specific
c
c---------------------------------------------------------------------------------------

      module HVCommunicationLibrary

      use HVParallelDataMap
      use HVCommDataTypes
      include 'mpif.h'

      contains

c-----------------------------------------------------------------
c  Master distributes a 4D array of data in h-slice 
c     format to all workers, standard version
c-----------------------------------------------------------------
      subroutine distrib_h_4D_v1(ix,iy,iz,izloc,Ns,s,s_h)
c
      implicit none
c
      integer :: ix, iy, iz, izloc, Ns
      real :: s(ix,iy,iz,Ns)
      real :: s_h(ix,iy,izloc,Ns)
      integer :: i, ierr
      integer :: status(MPI_STATUS_SIZE, izloc)
      integer :: request(izloc)
c
      if (Master) then
        do i=1,iz
          call MPI_SEND(s(1,1,i,1), 1, GLOBAL_4D_HSLICE,
     &         owner_of_hslice(i), i,  MPI_COMM_WORLD, ierr)
        end do
      else if (HWorker) then
        do i=1,no_of_hslices(MyId)
         call MPI_IRECV(s_h(1,1,i,1), 1, 
     &                 LOCAL_4D_HSLICE, 0, global_hslice_id(myid,i),
     &                 MPI_COMM_WORLD, request(i), ierr)
        end do
	call MPI_WAITALL(no_of_hslices(MyId), request, status, ierr)
      end if
c
c
      end subroutine  distrib_h_4D_v1



c-----------------------------------------------------------------
c  Master distributes vcolumns to all workers,
c     vcolumn by vcolumn, standard version
c-----------------------------------------------------------------
      subroutine distrib_v_4D_v1(ix,iy,iz,icloc,Ns,s,s_v)
c
      implicit none
c
      integer :: ix, iy, iz, icloc, Ns
      real :: s(ix,iy,iz,Ns)
      real :: s_v(1,icloc,iz,Ns)
      integer :: i, ierr
      integer :: status(MPI_STATUS_SIZE, icloc)
      integer :: request(icloc)

      if (Master) then
        do i=1,ix*iy
           call MPI_SEND(s(planar_vcol_id(i,1),planar_vcol_id(i,2),1,1), 
     & 		1, GLOBAL_4D_VCOL, owner_of_vcol(i), local_vcol_id(i), 
     &		MPI_COMM_WORLD, ierr)
        end do
      else if ( VWorker ) then
        do i=1,no_of_vcols(myid)
 	  call MPI_IRECV(s_v(1,i,1,1), 1, LOCAL_4D_VCOL, 
     &		0, i, MPI_COMM_WORLD, request(i), ierr)
	end do
	call MPI_WAITALL(no_of_vcols(MyId), request, status, ierr)
      end if
      
      end subroutine  distrib_v_4D_v1
      
c-----------------------------------------------------------------
c  Master gathers a 4D array of data in h-slice format from all workers
c  standard version
c-----------------------------------------------------------------
      subroutine gather_h_4D_v1(ix,iy,iz,izloc,Ns,s,s_h)
c      
      implicit none
c
      integer :: ix, iy, iz, izloc, Ns
      real   :: s(ix,iy,iz,Ns)
      real   :: s_h(ix,iy,izloc,Ns)
      integer :: i, ierr
      integer,pointer,dimension(:,:),save :: status
      integer,pointer,dimension(:),save :: request
c     
      if (Master) then
        if (.not. allocated(status)) then
          allocate( status(MPI_STATUS_SIZE,iz),STAT=ierr )
	end if
	if (.not. allocated(request)) then
          allocate( request(iz),STAT=ierr )      
	end if
        do i=1,iz
        call MPI_IRECV(s(1,1,i,1), 1, GLOBAL_4D_HSLICE, 
     &		owner_of_hslice(i), i,  MPI_COMM_WORLD, 
     &         request(i), ierr)
	end do
	call MPI_WAITALL(iz, request, status, ierr)
      else if (HWorker) then
        do i=1, no_of_hslices(myid)
          call MPI_SEND(s_h(1,1,i,1), 1, LOCAL_4D_HSLICE, 
     &                 0, global_hslice_id(MyId,i), 
     &                 MPI_COMM_WORLD, ierr)
	end do
      end if
      
      end subroutine  gather_h_4D_v1
        
c-----------------------------------------------------------------
c  Master gathers a 4D array of data in v-column format from all workers
c  standard version
c-----------------------------------------------------------------
      subroutine gather_v_4D_v1(ix,iy,iz,icloc,Ns,s,s_v)
c      
      implicit none
c
      integer :: ix, iy, iz, icloc, Ns
      real   :: s(ix,iy,iz,Ns)
      real   :: s_v(1,icloc,iz,Ns)
      integer :: i,ierr
      integer,pointer,dimension(:,:),save :: status
      integer,pointer,dimension(:),save :: request
c     
      if (Master) then
        if (.not. allocated(status)) then
          allocate( status(MPI_STATUS_SIZE,ix*iy),STAT=ierr )
	end if
	if (.not. allocated(request)) then
          allocate( request(ix*iy),STAT=ierr )      
	end if
        do i=1,ix*iy
          call MPI_IRECV(s(planar_vcol_id(i,1),planar_vcol_id(i,2),1,1),
     &		1, GLOBAL_4D_VCOL, owner_of_vcol(i), local_vcol_id(i), 
     &		MPI_COMM_WORLD, request(i), ierr)
	end do
	call MPI_WAITALL(ix*iy, request, status, ierr)
      else if (VWorker) then
        do i=1,no_of_vcols(myid)
 	  call MPI_SEND(s_v(1,i,1,1), 1, LOCAL_4D_VCOL, 
     &		0, i, MPI_COMM_WORLD, ierr)
	end do
      end if

      end subroutine  gather_v_4D_v1

c---------------------------------------------------------
c      Slaves exchange data:  - old format is h-slices (for sending)
c                             - new format is v-columns (for receiving)
c	using copying, sending blocks (single message to each proc)
c---------------------------------------------------------
      subroutine shuffle_h2v_4D_v1(ix,iy,iz,izloc,icloc,Ns,s_h,s_v)  
c      
      implicit none
c
      integer,intent(in) :: ix, iy, iz, izloc, icloc, Ns
      real,intent(in)  :: s_h(ix,iy,izloc,Ns)
      real,intent(out) :: s_v(1,icloc,iz,Ns)
      integer :: i, ierr, gid, p
      integer :: status(MPI_STATUS_SIZE)
      integer :: sstatus(MPI_STATUS_SIZE,NVworkers),srequest(NVworkers)
      real, pointer,dimension(:,:,:,:),save    :: sbuf
      real, pointer,dimension(:,:,:),save :: rbuf

      if (HWorker) then
        if (.not.allocated(sbuf)) then
	  allocate(sbuf(icloc, izloc, Ns,NVworkers),STAT=ierr)
	end if
        do p=1,NVworkers
          do i=1,no_of_vcols(p)
	    gid = global_vcol_id(p,i)
	    sbuf(i,1:no_of_hslices(MyId),1:Ns,p) = 
     &           s_h(planar_vcol_id(gid,1),planar_vcol_id(gid,2),
     &		1:no_of_hslices(MyId),1:Ns)
          end do
          call MPI_ISEND(sbuf(1,1,1,p), icloc*izloc*Ns, MPI_REAL, 
     &          p, p, MPI_COMM_WORLD, srequest(p), ierr)
        end do
      end if
      
      if (VWorker) then	
        if (.not.allocated(rbuf)) then
	  allocate(rbuf(icloc, izloc, Ns),STAT=ierr)
	end if
        do p=1,NHworkers
          call MPI_RECV(rbuf, icloc*izloc*Ns, MPI_REAL, 
     &             p, myid, MPI_COMM_WORLD, status, ierr)
	  s_v(1,1:no_of_vcols(MyId),
     &              owned_hslices(p,1:no_of_hslices(p)),1:Ns)=
     &	            rbuf(1:no_of_vcols(MyId),1:no_of_hslices(p),1:Ns)
        end do
      end if
      
      if (HWorker) call MPI_WAITALL(NVworkers,srequest,sstatus,ierr)
      
      end subroutine shuffle_h2v_4D_v1
      
c---------------------------------------------------------
c      Slaves exchange data:  - old format is v-columns (for sending)
c                             - new format is h-slices (for receiving)
c	using copying, sending blocks (single message to each proc)
c---------------------------------------------------------
      subroutine shuffle_v2h_4D_v1(ix,iy,iz,izloc,icloc,Ns,s_h,s_v)  
c      
      implicit none
c
      integer,intent(in) :: ix, iy, iz, izloc, icloc, Ns
      real,intent(out) :: s_h(ix,iy,izloc,Ns)
      real,intent(in)  :: s_v(1,icloc,iz,Ns)
      integer :: i, ierr, p
      integer :: status(MPI_STATUS_SIZE)
      integer :: sstatus(MPI_STATUS_SIZE,NHworkers),srequest(NHworkers)
      real, pointer, dimension(:,:,:,:),save    :: sbuf
      real, pointer, dimension(:,:,:),save    :: rbuf

      if (VWorker) then
        if (.not.allocated(sbuf)) then
          allocate( sbuf(icloc,izloc,Ns,NHworkers),STAT=ierr )
        endif
c      
        do p=1,NHworkers
	  sbuf(1:icloc,1:no_of_hslices(p),1:Ns,p)=
     &	      s_v(1,1:icloc,owned_hslices(p,1:no_of_hslices(p)),1:Ns)
          call MPI_ISEND(sbuf(1,1,1,p), icloc*izloc*Ns, MPI_REAL, 
     &          p, p,  
     &          MPI_COMM_WORLD, srequest(p), ierr)
        end do
      end if
      
      if (HWorker) then      
        if (.not.allocated(rbuf)) then
          allocate( rbuf(icloc,izloc,Ns),STAT=ierr )
        endif
        do p=1,NVworkers
          call MPI_RECV(rbuf, icloc*izloc*Ns, MPI_REAL, 
     &           p, myid, MPI_COMM_WORLD, status, ierr)
          do i = 1,no_of_vcols(p)
            s_h(planar_vcol_id(global_vcol_id(p,i),1),
     &		planar_vcol_id(global_vcol_id(p,i),2),
     &       	1:no_of_hslices(MyId),1:Ns)=
     &	        rbuf(i,1:no_of_hslices(MyId),1:Ns)
          enddo
        end do

      end if

      if (VWorker) call MPI_WAITALL(NHworkers,srequest,sstatus,ierr)

      end subroutine shuffle_v2h_4D_v1

c--------------------------------------------------------------------------
c  Master distributes a 2D (ix.iy) array of data in h-slice format 
c  to all workers - not necessary for STEM-III
c--------------------------------------------------------------------------
c      
      subroutine distrib_h_2D_v1(ix,iy,s)

      implicit none
c
      integer :: ix, iy
      real :: s(ix,iy)
      integer :: ierr
c      
      call MPI_BCAST(s, ix*iy, MPI_REAL, 0, MPI_COMM_WORLD, ierr)
      call MPI_BARRIER(MPI_COMM_WORLD, ierr)
c
      end subroutine  distrib_h_2D_v1
      
c--------------------------------------------------------------------------
c  Master distributes a 2D (ix.iy) array of data in v-column format 
c  to all workers 
c--------------------------------------------------------------------------
c      
      subroutine distrib_v_2D_v1(ix,iy,icloc,s,s_v)

      implicit none
c
      integer :: ix, iy, icloc
      real :: s(ix,iy)
      real :: s_v(1,icloc)
      integer :: i, ierr
      integer :: status(MPI_STATUS_SIZE, icloc)
      integer :: request(icloc)

      if (Master) then
        do i=1,ix*iy
          call MPI_SEND(s(planar_vcol_id(i,1),planar_vcol_id(i,2)), 
     &		1, MPI_REAL, 
     &         owner_of_vcol(i), i,  MPI_COMM_WORLD, ierr)
          if (ierr.ne.MPI_SUCCESS) then
	    print*, 'Error: distrib_v_2D: send failed'
	    stop
	  end if
	end do
      else if (VWorker) then
        do i=1, no_of_vcols(myid)
          call MPI_IRECV(s_v(1,i), 1, MPI_REAL, 0, 
     &			global_vcol_id(myid,i), 
     &                 MPI_COMM_WORLD, request(i), ierr)
          if (ierr.ne.MPI_SUCCESS) then
	    print*, 'Error: distrib_v_2D: receive failed'
	    stop
	  end if
	end do
	call MPI_WAITALL(no_of_vcols(MyId), request, status, ierr)
      end if
c
      end subroutine  distrib_v_2D_v1
      
c-------------------------------------------------------------------------
c  Master distributes a BDz (ix.iy.is) array of data in v-column format 
c   to all workers - only needed for v-transport
c-------------------------------------------------------------------------
c      
      subroutine distrib_v_BDz_v1(ix,iy,icloc,Ns,s,s_v)

      implicit none
c
      integer :: ix, iy, Ns, icloc
      real   :: s(ix,iy,Ns)
      real   :: s_v(1,icloc,Ns)
      integer :: j, ierr
      integer :: status(MPI_STATUS_SIZE, icloc)
      integer :: request(icloc)

      if (Master) then
        do j=1,ix*iy
        call MPI_SEND(s(planar_vcol_id(j,1),planar_vcol_id(j,2),1), 
     &		1, GLOBAL_BDz_VCOL, 
     &         owner_of_vcol(j), j,  MPI_COMM_WORLD, ierr)
	end do
      else if (VWorker) then
        do j=1, no_of_vcols(myid)
        call MPI_IRECV(s_v(1,j,1), 1, LOCAL_BDz_VCOL, 0, 
     &			global_vcol_id(myid,j), 
     &          	MPI_COMM_WORLD, request(j), ierr)
	end do
	call MPI_WAITALL(no_of_vcols(MyId), request, status, ierr)
      end if
c
      end subroutine  distrib_v_BDz_v1
      
c-------------------------------------------------------------------------
c  Master distributes a 3D (ix.iy.iz) array of data in h-slice format 
c   to all workers
c-------------------------------------------------------------------------
c      
      subroutine distrib_h_3D_v1(ix,iy,iz,izloc,s,s_h)

      implicit none
c
      integer :: ix, iy, iz, izloc
      real :: s(ix,iy,iz)
      real :: s_h(ix,iy,izloc)
      integer :: i, ierr
      integer :: status(MPI_STATUS_SIZE, izloc)
      integer :: request(izloc)

      if (Master) then
        do i=1,iz
          call MPI_SEND(s(1,1,i), ix*iy, MPI_REAL, 
     &         owner_of_hslice(i), i,  MPI_COMM_WORLD, ierr)
          if (ierr.ne.MPI_SUCCESS) then
	    print*, 'Error: distrib_h_3D: send failed'
	    stop
	  end if
	end do
      else if (HWorker) then
        do i=1, no_of_hslices(myid)
          call MPI_IRECV(s_h(1,1,i), ix*iy, MPI_REAL, 0, 
     &			global_hslice_id(myid,i), 
     &                 MPI_COMM_WORLD, request(i), ierr)
          if (ierr.ne.MPI_SUCCESS) then
	    print*, 'Error: distrib_h_3D: receive failed'
	    stop
	  end if
	end do
	call MPI_WAITALL(no_of_hslices(MyId), request, status, ierr)
      end if
c
      end subroutine  distrib_h_3D_v1    

      
c--------------------------------------------------------------------------
c  Master distributes a 3D (ix.iy.iz) array of data in v-column format 
c   to all workers
c--------------------------------------------------------------------------
c
      subroutine distrib_v_3D_v1(ix,iy,iz,icloc,s,s_v)

      implicit none
c      
      integer :: ix, iy, iz, icloc
      real   :: s(ix,iy,iz)
      real   :: s_v(1,icloc,iz)
      integer :: i, ierr
      integer :: status(MPI_STATUS_SIZE, icloc)
      integer :: request(icloc)

      if (Master) then
        do i=1,ix*iy
           call MPI_SEND(s(planar_vcol_id(i,1),planar_vcol_id(i,2),1), 
     &		1, GLOBAL_3D_VCOL,
     &         owner_of_vcol(i), i,  MPI_COMM_WORLD, ierr)
        end do
      else if (VWorker) then
        do i=1, no_of_vcols(myid)
           call MPI_IRECV(s_v(1,i,1), 1, LOCAL_3D_VCOL, 0, 
     &			global_vcol_id(myid,i), 
     &                 MPI_COMM_WORLD, request(i), ierr)
        end do
	call MPI_WAITALL(no_of_vcols(MyId), request, status, ierr)
      endif
c	
      return
      end subroutine distrib_v_3D_v1 

c---------------------------------------------------------------------------
c  Master distributes the y-boundary array BDy
c---------------------------------------------------------------------------
c      
      subroutine distrib_h_BDy_v1(ix,iz,izloc,Ns,s,sloc)

      implicit none
c
      integer :: ix, iz, Ns, izloc
      real :: s(ix,iz,Ns)
      real :: sloc(ix,izloc,Ns)
      integer :: j, ierr
      integer :: status(MPI_STATUS_SIZE, izloc)
      integer :: request(izloc)

      if (Master) then
        do j=1,iz
        call MPI_SEND(s(1,j,1), 1, GLOBAL_BDy_HSLICE, 
     &         owner_of_hslice(j), j,  MPI_COMM_WORLD, ierr)
	end do
      else if (HWorker) then
        do j=1, no_of_hslices(myid)
        call MPI_IRECV(sloc(1,j,1), 1, LOCAL_BDy_HSLICE, 0, 
     &			global_hslice_id(myid,j), 
     &                 MPI_COMM_WORLD, request(j), ierr)
	end do
	call MPI_WAITALL(no_of_hslices(MyId), request, status, ierr)
      end if
c
      end subroutine  distrib_h_BDy_v1
      
c---------------------------------------------------------------------------
c  Master distributes the y-boundary array BDy2 
c  - careful, X is used from
c  former different naming BDx...
c---------------------------------------------------------------------------
c      
      subroutine distrib_xh_BD_v1(ix,iz,izloc,Ns,s,sloc)

      implicit none
c
      integer :: ix, iz, Ns, izloc
      real :: s(ix,iz,2,Ns)
      real :: sloc(ix,izloc,2,Ns)
      integer :: j, ierr
      integer :: status(MPI_STATUS_SIZE, izloc)
      integer :: request(izloc)

      if (Master) then
        do j=1,iz
        call MPI_SEND(s(1,j,1,1), 1, GLOBAL_BD_XHSLICE, 
     &         owner_of_hslice(j), j,  MPI_COMM_WORLD, ierr)
	end do
      else if (HWorker) then
        do j=1, no_of_hslices(myid)
        call MPI_IRECV(sloc(1,j,1,1), 1, LOCAL_BD_XHSLICE, 0, 
     &			global_hslice_id(myid,j), 
     &                 MPI_COMM_WORLD, request(j), ierr)
	end do
	call MPI_WAITALL(no_of_hslices(MyId), request, status, ierr)
      end if
c
      end subroutine  distrib_xh_BD_v1
      
c---------------------------------------------------------------------------
c  Master distributes the x-boundary array BDx
c---------------------------------------------------------------------------
c      
      subroutine distrib_h_BDx_v1(iy,iz,izloc,Ns,s,sloc)

      implicit none
c
      integer :: iy, iz, Ns, izloc
      real :: s(iy,iz,Ns)
      real :: sloc(iy,izloc,Ns)
      integer :: i, ierr
      integer :: status(MPI_STATUS_SIZE, izloc)
      integer :: request(izloc)

      if (Master) then
        do i=1,iz
        call MPI_SEND(s(1,i,1), 1, GLOBAL_BDx_HSLICE, 
     &         owner_of_hslice(i), i,  MPI_COMM_WORLD, ierr)
	end do
      else if (HWorker) then
        do i=1, no_of_hslices(myid)
        call MPI_IRECV(sloc(1,i,1), 1, LOCAL_BDx_HSLICE, 0,
     &			global_hslice_id(myid,i), 
     &                 MPI_COMM_WORLD, request(i), ierr)
	end do
	call MPI_WAITALL(no_of_hslices(MyId), request, status, ierr)
      end if
c
      end subroutine  distrib_h_BDx_v1

c---------------------------------------------------------------------------
c  Master distributes the x-boundary array BDx2
c  - careful, Y is used from
c  former different naming BDy...
c---------------------------------------------------------------------------
c      
      subroutine distrib_yh_BD_v1(iy,iz,izloc,Ns,s,sloc)

      implicit none
c
      integer :: iy, iz, Ns, izloc
      real :: s(iy,iz,2,Ns)
      real :: sloc(iy,izloc,2,Ns)
      integer :: i, ierr
      integer :: status(MPI_STATUS_SIZE, izloc)
      integer :: request(izloc)

      if (Master) then
        do i=1,iz
        call MPI_SEND(s(1,i,1,1), 1, GLOBAL_BD_YHSLICE, 
     &         owner_of_hslice(i), i,  MPI_COMM_WORLD, ierr)
	end do
      else if (HWorker) then
        do i=1, no_of_hslices(myid)
        call MPI_IRECV(sloc(1,i,1,1), 1, LOCAL_BD_YHSLICE, 0,
     &			global_hslice_id(myid,i), 
     &                 MPI_COMM_WORLD, request(i), ierr)
	end do
	call MPI_WAITALL(no_of_hslices(MyId), request, status, ierr)
      end if
c
      end subroutine  distrib_yh_BD_v1

c *******************************************************************
c ***************** single message communication ********************
c *******************************************************************


c-----------------------------------------------------------------
c  Master distributes a 4D array of data in h-slice 
c     format to all workers
c Less communication than other versions - one message to each proc
c-----------------------------------------------------------------
      subroutine distrib_h_4D_v2(ix,iy,iz,izloc,Ns,s,s_h)
c
      implicit none
c
      integer :: ix, iy, iz, izloc, Ns
      real :: s(ix,iy,iz,Ns)
      real :: s_h(ix,iy,izloc,Ns)
      integer :: i, ierr, status(MPI_STATUS_SIZE)

      if (Master) then
        do i=1,NHworkers
 	  if (no_of_hslices(i) .eq. izloc) then
            call MPI_SEND(s(1,1,i,1), 1, GLOBAL_LARGE_4D_HSLICES, 
     &         owner_of_hslice(i), i,  MPI_COMM_WORLD, ierr)
	  else
            call MPI_SEND(s(1,1,i,1), 1, GLOBAL_SMALL_4D_HSLICES, 
     &         owner_of_hslice(i), i,  MPI_COMM_WORLD, ierr)
	  end if
	end do
      else if (Hworker) then
	if (no_of_hslices(myid) .eq. izloc) then
          call MPI_RECV(s_h(1,1,1,1), 1, LOCAL_LARGE_4D_HSLICES, 0,
     &		myid, MPI_COMM_WORLD, status, ierr)
	else
          call MPI_RECV(s_h(1,1,1,1), 1, LOCAL_SMALL_4D_HSLICES, 0,
     &		myid, MPI_COMM_WORLD, status, ierr)
        end if
      end if
      
      end subroutine  distrib_h_4D_v2

c-----------------------------------------------------------------
c  Master distributes a 4D array of data in v-column format to all workers
c  all in one message
c-----------------------------------------------------------------
      subroutine distrib_v_4D_v2(ix,iy,iz,icloc,Ns,s,s_v)
c      
      implicit none
c
      integer :: ix, iy, iz, icloc, Ns
      real :: s(ix,iy,iz,Ns)
      real :: s_v(1,icloc,iz,Ns)
      integer :: i, ierr, status(MPI_STATUS_SIZE)

      if (Master) then
        do i=1,NVworkers
	  if (no_of_vcols(i) .eq. icloc) then
c	print*,'send to',owner_of_vcol(i),'no',i,'large'
	    call MPI_SEND(s(planar_vcol_id(i,1),
     &		planar_vcol_id(i,2),1,1), 1, GLOBAL_LARGE_4D_VCOLS, 
     &         owner_of_vcol(i), i,  MPI_COMM_WORLD, ierr)
	  else
c	print*,'send to',owner_of_vcol(i),'no',i,'small'
	    call MPI_SEND(s(planar_vcol_id(i,1),
     &		planar_vcol_id(i,2),1,1), 1, GLOBAL_SMALL_4D_VCOLS, 
     &         owner_of_vcol(i), i,  MPI_COMM_WORLD, ierr)
          end if
	end do
      else if (VWorker) then
	if (no_of_vcols(myid) .eq. icloc) then
          call MPI_RECV(s_v(1,1,1,1), 1, LOCAL_LARGE_4D_VCOLS, 0, 
     &		myid, MPI_COMM_WORLD, status, ierr)
	else
          call MPI_RECV(s_v(1,1,1,1), 1, LOCAL_SMALL_4D_VCOLS, 0, 
     &		myid, MPI_COMM_WORLD, status, ierr)
        end if
      end if

      end subroutine  distrib_v_4D_v2

c-----------------------------------------------------------------
c  Master gathers a 4D array of data in h-slice format from all workers
c-----------------------------------------------------------------
      subroutine gather_h_4D_v2(ix,iy,iz,izloc,Ns,s,s_h)
c      
      implicit none
c
      integer :: ix, iy, iz, izloc, Ns
      real   :: s(ix,iy,iz,Ns)
      real   :: s_h(ix,iy,izloc,Ns)
      integer :: i, ierr
      integer,pointer,dimension(:,:),save :: status
      integer,pointer,dimension(:),save :: request
      
      if (Master) then
        if (.not. allocated(status)) then
  	  allocate( status(MPI_STATUS_SIZE,NHworkers),STAT=ierr )
	end if
	if (.not. allocated(request)) then
          allocate( request(NHworkers),STAT=ierr )      
	end if
        do i=1,NHworkers
	  if ( no_of_hslices(i) .eq. izloc) then
	    call MPI_IRECV(s(1,1,i,1), 1, GLOBAL_LARGE_4D_HSLICES, 
     &         owner_of_hslice(i), i,  MPI_COMM_WORLD, 
     &         request(i), ierr)
     	  else
	    call MPI_IRECV(s(1,1,i,1), 1, GLOBAL_SMALL_4D_HSLICES, 
     &         owner_of_hslice(i), i,  MPI_COMM_WORLD, 
     &         request(i), ierr)
     	  end if
	end do
	call MPI_WAITALL(NHworkers, request, status, ierr)

      else if (HWorker) then
      	if (no_of_hslices(myid) .eq. izloc) then
	  call MPI_SEND(s_h(1,1,1,1), 1, LOCAL_LARGE_4D_HSLICES, 
     &                 0, myid, 
     &                 MPI_COMM_WORLD, ierr)
        else
	  call MPI_SEND(s_h(1,1,1,1), 1, LOCAL_SMALL_4D_HSLICES, 
     &                 0, myid, 
     &                 MPI_COMM_WORLD, ierr)
     	end if
      end if
      
      end subroutine  gather_h_4D_v2

c-----------------------------------------------------------------
c  Master gathers a 4D array of data in v-column format from all workers
c-----------------------------------------------------------------
      subroutine gather_v_4D_v2(ix,iy,iz,icloc,Ns,s,s_v)
c      
      implicit none
c
      integer :: ix, iy, iz, icloc, Ns
      real   :: s(ix,iy,iz,Ns)
      real   :: s_v(1,icloc,iz,Ns)
      integer :: i, ierr
      integer,pointer,dimension(:,:),save :: status
      integer,pointer,dimension(:),save :: request
      
      if (Master) then
        if (.not. allocated(status)) then
          allocate( status(MPI_STATUS_SIZE,NVworkers),STAT=ierr )
	end if
	if (.not. allocated(request)) then
          allocate( request(NVworkers),STAT=ierr )      
	end if
        do i=1,NVworkers
	  if ( no_of_vcols(i) .eq. icloc) then
	    call MPI_IRECV(s(i,1,1,1), 1, GLOBAL_LARGE_4D_VCOLS, 
     &         owner_of_vcol(i), i,  MPI_COMM_WORLD, 
     &         request(i), ierr)
     	  else
	    call MPI_IRECV(s(i,1,1,1), 1, GLOBAL_SMALL_4D_VCOLS, 
     &         owner_of_vcol(i), i,  MPI_COMM_WORLD, 
     &         request(i), ierr)
     	  end if
	end do
	call MPI_WAITALL(NVworkers, request, status, ierr)

      else if (VWorker) then
      	if ( no_of_vcols(myid) .eq. icloc ) then
          call MPI_SEND(s_v(1,1,1,1), 1, LOCAL_LARGE_4D_VCOLS, 
     &                 0, myid, 
     &                 MPI_COMM_WORLD, ierr)
        else
          call MPI_SEND(s_v(1,1,1,1), 1, LOCAL_SMALL_4D_VCOLS, 
     &                 0, myid, 
     &                 MPI_COMM_WORLD, ierr)
	end if
      end if
      
      end subroutine  gather_v_4D_v2  

c---------------------------------------------------------
c      Slaves exchange data:  - old format is h-slices (for sending)
c                             - new format is v-columns (for receiving)
c---------------------------------------------------------
      subroutine shuffle_h2v_4D_v2(ix,iy,iz,izloc,icloc,Ns,s_h,s_v)  
c      
      implicit none
c
      integer :: ix, iy, iz, izloc, icloc, Ns
      real   :: s_h(ix,iy,izloc,Ns)
      real   :: s_v(1,icloc,iz,Ns)
      integer :: i, j, ierr
      integer :: status(MPI_STATUS_SIZE,NHworkers)
      integer :: request(NHworkers)

      if (VWorker) then      
      if (no_of_vcols(myid) .eq. icloc) then
        do j=1, NHworkers
	  if (no_of_hslices(j) .eq. izloc) then
c	print*,'proc',myid,'receives LHLV from',j
    	    call MPI_IRECV(s_v(1,1,j,1),1, LOCAL_4D_LHLV_VCOLS,
     &		j, j, MPI_COMM_WORLD, request(j), ierr)
     	  else
c	print*,'proc',myid,'receives SHLV from',j
     	    call MPI_IRECV(s_v(1,1,j,1),1,LOCAL_4D_SHLV_VCOLS,
     &		j, j, MPI_COMM_WORLD, request(j), ierr)
          end if
	end do
      else
        do j=1, NHworkers
	  if (no_of_hslices(j) .eq. izloc) then	
c 	print*,'proc',myid,'receives LHSV from',j
    	    call MPI_IRECV(s_v(1,1,j,1),1,LOCAL_4D_LHSV_VCOLS,
     &		j, j, MPI_COMM_WORLD, request(j), ierr)
          else
c	print*,'proc',myid,'receives SHSV from',j
     	    call MPI_IRECV(s_v(1,1,j,1),1,LOCAL_4D_SHSV_VCOLS,
     &		j, j, MPI_COMM_WORLD, request(j), ierr)
     	  end if
	end do
      end if
      end if
	
      if (HWorker) then
      if (no_of_hslices(myid) .eq. izloc) then
        do i=1, NVworkers
	  if (no_of_vcols(i) .eq. icloc) then
c	print*,'proc',myid,'sends LHLV to',i
       	    call MPI_SEND(s_h(i,1,1,1),1,LOCAL_4D_LHLV_HCOLS,
     &		i, myid, MPI_COMM_WORLD, ierr)	  	
     	  else
c	print*,'proc',myid,'sends LHSV to',i
       	    call MPI_SEND(s_h(i,1,1,1),1,LOCAL_4D_LHSV_HCOLS,
     &		i, myid, MPI_COMM_WORLD, ierr)
          end if
	end do
      else
        do i=1, NVworkers
c	print*,'proc',myid,'sends SHLV to',i
	  if (no_of_vcols(i) .eq. icloc) then	
       	    call MPI_SEND(s_h(i,1,1,1),1,LOCAL_4D_SHLV_HCOLS,
     &		i, myid, MPI_COMM_WORLD, ierr)
     	  else
c	print*,'proc',myid,'sends SHSV to',i
       	    call MPI_SEND(s_h(i,1,1,1),1,LOCAL_4D_SHSV_HCOLS,
     &		i, myid, MPI_COMM_WORLD, ierr)
          end if
	end do
      end if
      end if

      if (VWorker) call MPI_WAITALL(NHworkers, request, status, ierr)

      end subroutine shuffle_h2v_4D_v2

c---------------------------------------------------------
c      Slaves exchange data:  - old format is v-columns (for sending)
c                             - new format is h-slices (for receiving)
c---------------------------------------------------------
      subroutine shuffle_v2h_4D_v2(ix,iy,iz,izloc,icloc,Ns,s_h,s_v)  
c      
      implicit none
c
      integer :: ix, iy, iz, izloc, icloc, Ns
      real   :: s_h(ix,iy,izloc,Ns)
      real   :: s_v(1,icloc,iz,Ns)
      integer :: i, j, Ierr
      integer :: status(MPI_STATUS_SIZE,NVworkers)
      integer :: request(NVworkers)

      if (HWorker) then
      if (no_of_hslices(myid) .eq. izloc) then
        do i=1, NVworkers
	  if (no_of_vcols(i) .eq. icloc) then
c	print*,'proc',myid,'receives LHLV from',i
       	    call MPI_IRECV(s_h(i,1,1,1),1,LOCAL_4D_LHLV_HCOLS,
     &		i, i, MPI_COMM_WORLD, request(i), ierr)	  	
     	  else
c	print*,'proc',myid,'receives LHSV from',i
       	    call MPI_IRECV(s_h(i,1,1,1),1,LOCAL_4D_LHSV_HCOLS,
     &		i, i, MPI_COMM_WORLD, request(i), ierr)
          end if
	end do
      else
        do i=1, NVworkers
c	print*,'proc',myid,'receives SHLV from',i
	  if (no_of_vcols(i) .eq. icloc) then	
       	    call MPI_IRECV(s_h(i,1,1,1),1,LOCAL_4D_SHLV_HCOLS,
     &		i, i, MPI_COMM_WORLD, request(i), ierr)
     	  else
c	print*,'proc',myid,'receives SHSV from',i
       	    call MPI_IRECV(s_h(i,1,1,1),1,LOCAL_4D_SHSV_HCOLS,
     &		i, i, MPI_COMM_WORLD, request(i), ierr)
          end if
	end do
      end if
      end if
      
      if (VWorker) then
      if (no_of_vcols(myid) .eq. icloc) then
        do j=1, NHworkers
	  if (no_of_hslices(j) .eq. izloc) then
c	print*,'proc',myid,' sends LHLV to',j
    	    call MPI_SEND(s_v(1,1,j,1),1, LOCAL_4D_LHLV_VCOLS,
     &		j, myid, MPI_COMM_WORLD, ierr)
     	  else
c	print*,'proc',myid,'sends SHLV to',j
     	    call MPI_SEND(s_v(1,1,j,1),1,LOCAL_4D_SHLV_VCOLS,
     &		j, myid, MPI_COMM_WORLD, ierr)
          end if
	end do
      else
        do j=1, NHworkers
	  if (no_of_hslices(j) .eq. izloc) then	
c 	print*,'proc',myid,'sends LHSV to',j
    	    call MPI_SEND(s_v(1,1,j,1),1,LOCAL_4D_LHSV_VCOLS,
     &		j, myid, MPI_COMM_WORLD, ierr)
          else
c	print*,'proc',myid,'sends SHSV to',j
     	    call MPI_SEND(s_v(1,1,j,1),1,LOCAL_4D_SHSV_VCOLS,
     &		j, myid, MPI_COMM_WORLD, ierr)
     	  end if
	end do
      end if
      end if

      if (HWorker) call MPI_WAITALL(NVworkers, request, status, ierr)

      end subroutine shuffle_v2h_4D_v2

c--------------------------------------------------------------------------
c  Master distributes a 2D (ix.iy) array of data in vcolumn format 
c  to all workers 
c--------------------------------------------------------------------------
c      
      subroutine distrib_v_2D_v2(ix,iy,icloc,s,s_v)

      implicit none
c
      integer :: ix, iy, icloc
      real :: s(ix,iy)
      real :: s_v(1,icloc)
      integer :: i, ierr, status(MPI_STATUS_SIZE)
c      
      if (Master) then
        do i=1,NVworkers
	  if (no_of_vcols(i) .eq. icloc) then
            call MPI_SEND(s(i,1), 1, GLOBAL_LARGE_2D_VCOLS, 
     &         owner_of_vcol(i), i,  MPI_COMM_WORLD, ierr)
          else
            call MPI_SEND(s(i,1), 1, GLOBAL_SMALL_2D_VCOLS, 
     &         owner_of_vcol(i), i,  MPI_COMM_WORLD, ierr)
          end if
	end do
      else if (VWorker) then
        if (no_of_vcols(myid) .eq. icloc) then
          call MPI_RECV(s_v(1,1), 1, LOCAL_LARGE_2D_VCOLS, 0, 
     &			myid, MPI_COMM_WORLD, status, ierr)
        else
          call MPI_RECV(s_v(1,1), 1, LOCAL_SMALL_2D_VCOLS, 0, 
     &			myid, MPI_COMM_WORLD, status, ierr)
	end if
      end if
c
      end subroutine  distrib_v_2D_v2
      
c-------------------------------------------------------------------------
c  Master distributes a BDz (ix.iy.is) array of data in vcolumn format 
c   to all workers - only needed for v-transport
c-------------------------------------------------------------------------
c      
      subroutine distrib_v_BDz_v2(ix,iy,icloc,Ns,s,s_v)

      implicit none
c
      integer :: ix, iy, Ns, icloc
      real   :: s(ix,iy,Ns)
      real   :: s_v(1,icloc,Ns)
      integer :: j, ierr, status(MPI_STATUS_SIZE)
c      
      if (Master) then
        do j=1,NVworkers
	  if (no_of_vcols(j) .eq. icloc) then
            call MPI_SEND(s(j,1,1), 1, GLOBAL_LARGE_BDz_VCOLS, 
     &         owner_of_vcol(j), j,  MPI_COMM_WORLD, ierr)
          else
            call MPI_SEND(s(j,1,1), 1, GLOBAL_SMALL_BDz_VCOLS, 
     &         owner_of_vcol(j), j,  MPI_COMM_WORLD, ierr)
	  end if
	end do
      else if (VWorker) then
        if (no_of_vcols(myid) .eq. icloc) then
          call MPI_RECV(s_v(1,1,1), 1, LOCAL_LARGE_BDz_VCOLS, 0, 
     &			myid, MPI_COMM_WORLD, status, ierr)
        else
          call MPI_RECV(s_v(1,1,1), 1, LOCAL_SMALL_BDz_VCOLS, 0, 
     &			myid, MPI_COMM_WORLD, status, ierr)
	end if
      end if
c
      end subroutine  distrib_v_BDz_v2
      
c-------------------------------------------------------------------------
c  Master distributes a 3D (ix.iy.iz) array of data in h-slice format 
c   to all workers
c-------------------------------------------------------------------------
c      
      subroutine distrib_h_3D_v2(ix,iy,iz,izloc,s,s_h)

      implicit none
c
      integer :: ix, iy, iz, izloc
      real :: s(ix,iy,iz)
      real :: s_h(ix,iy,izloc)
      integer :: i, ierr, status(MPI_STATUS_SIZE)
c      
      if (Master) then
        do i=1,NHworkers
	  if (no_of_hslices(i) .eq. izloc) then
            call MPI_SEND(s(1,1,i), 1, GLOBAL_LARGE_3D_HSLICES, 
     &         owner_of_hslice(i), i,  MPI_COMM_WORLD, ierr)
          else
            call MPI_SEND(s(1,1,i), 1, GLOBAL_SMALL_3D_HSLICES, 
     &         owner_of_hslice(i), i,  MPI_COMM_WORLD, ierr)
          end if
	end do
      else if (HWorker) then
	if (no_of_hslices(myid) .eq. izloc) then
          call MPI_RECV(s_h(1,1,1), 1, LOCAL_LARGE_3D_HSLICES, 0, 
     &			myid, MPI_COMM_WORLD, status, ierr)
        else
          call MPI_RECV(s_h(1,1,1), 1, LOCAL_SMALL_3D_HSLICES, 0, 
     &			myid, MPI_COMM_WORLD, status, ierr)
	end if
      end if
c
      end subroutine  distrib_h_3D_v2   

      
c--------------------------------------------------------------------------
c  Master distributes a 3D (ix.iy.iz) array of data in v-column format 
c   to all workers
c--------------------------------------------------------------------------
c
      subroutine distrib_v_3D_v2(ix,iy,iz,icloc,s,s_v)

      implicit none
c      
      integer :: ix, iy, iz, icloc
      real   :: s(ix,iy,iz)
      real   :: s_v(1,icloc,iz)
      integer :: i, ierr, status(MPI_STATUS_SIZE)
c      
      if (Master) then
        do i=1,NVworkers
 	  if (no_of_vcols(i) .eq. icloc) then
            call MPI_SEND(s(i,1,1), 1, GLOBAL_LARGE_3D_VCOLS,
     &         owner_of_vcol(i), i,  MPI_COMM_WORLD, ierr)
          else
            call MPI_SEND(s(i,1,1), 1, GLOBAL_SMALL_3D_VCOLS,
     &         owner_of_vcol(i), i,  MPI_COMM_WORLD, ierr)
          end if
        end do
      else if (VWorker) then
        if (no_of_vcols(myid) .eq. icloc) then
           call MPI_RECV(s_v(1,1,1), 1, LOCAL_LARGE_3D_VCOLS, 0, 
     &			myid, MPI_COMM_WORLD, status, ierr)
        else
           call MPI_RECV(s_v(1,1,1), 1, LOCAL_SMALL_3D_VCOLS, 0, 
     &			myid, MPI_COMM_WORLD, status, ierr)
        end if
      endif
c	
      return
      end subroutine distrib_v_3D_v2

c---------------------------------------------------------------------------
c  Master distributes the y-boundary array BDy
c---------------------------------------------------------------------------
c      
      subroutine distrib_h_BDy_v2(ix,iz,izloc,Ns,s,sloc)

      implicit none
c
      integer :: ix, iz, Ns, izloc
      real :: s(ix,iz,Ns)
      real :: sloc(ix,izloc,Ns)
      integer :: j, ierr, status(MPI_STATUS_SIZE)
c      
      if (Master) then
        do j=1,NHworkers
	  if (no_of_hslices(j) .eq. izloc) then
            call MPI_SEND(s(1,j,1), 1, GLOBAL_LARGE_BDy_HSLICES, 
     &         owner_of_hslice(j), j,  MPI_COMM_WORLD, ierr)
          else
            call MPI_SEND(s(1,j,1), 1, GLOBAL_SMALL_BDy_HSLICES, 
     &         owner_of_hslice(j), j,  MPI_COMM_WORLD, ierr)
          end if
	end do
      else if (HWorker) then
	if (no_of_hslices(myid) .eq. izloc) then
          call MPI_RECV(sloc(1,1,1), 1, LOCAL_LARGE_BDy_HSLICES, 0, 
     &			myid, MPI_COMM_WORLD, status, ierr)
        else
          call MPI_RECV(sloc(1,1,1), 1, LOCAL_SMALL_BDy_HSLICES, 0, 
     &			myid, MPI_COMM_WORLD, status, ierr)
	end if
      end if
c
      end subroutine  distrib_h_BDy_v2
 
c---------------------------------------------------------------------------
c  Master distributes the y-boundary array BDy2
c  - careful, X is used from
c  former different naming BDx...
c---------------------------------------------------------------------------
c      
      subroutine distrib_xh_BD_v2(ix,iz,izloc,Ns,s,sloc)

      implicit none
c
      integer :: ix, iz, Ns, izloc
      real :: s(ix,iz,2,Ns)
      real :: sloc(ix,izloc,2,Ns)
      integer :: j, ierr, status(MPI_STATUS_SIZE)
c      
      if (Master) then
        do j=1,NHworkers
	  if (no_of_hslices(j) .eq. izloc) then
            call MPI_SEND(s(1,j,1,1), 1, GLOBAL_LARGE_BD_XHSLICES, 
     &         owner_of_hslice(j), j,  MPI_COMM_WORLD, ierr)
          else
            call MPI_SEND(s(1,j,1,1), 1, GLOBAL_SMALL_BD_XHSLICES, 
     &         owner_of_hslice(j), j,  MPI_COMM_WORLD, ierr)
          end if
	end do
      else if (HWorker) then
	if (no_of_hslices(myid) .eq. izloc) then
          call MPI_RECV(sloc(1,1,1,1), 1, LOCAL_LARGE_BD_XHSLICES, 0, 
     &			myid, MPI_COMM_WORLD, status, ierr)
        else
          call MPI_RECV(sloc(1,1,1,1), 1, LOCAL_SMALL_BD_XHSLICES, 0, 
     &			myid, MPI_COMM_WORLD, status, ierr)
	end if
      end if
c
      end subroutine  distrib_xh_BD_v2
      
c---------------------------------------------------------------------------
c  Master distributes the x-boundary array BDx
c---------------------------------------------------------------------------
c      
      subroutine distrib_h_BDx_v2(iy,iz,izloc,Ns,s,sloc)

      implicit none
c
      integer :: iy, iz, Ns, izloc
      real :: s(iy,iz,Ns)
      real :: sloc(iy,izloc,Ns)
      integer :: i, ierr, status(MPI_STATUS_SIZE)
c      
      if (Master) then
        do i=1,NHworkers
	  if (no_of_hslices(i) .eq. izloc) then
            call MPI_SEND(s(1,i,1), 1, GLOBAL_LARGE_BDx_HSLICES, 
     &         owner_of_hslice(i), i,  MPI_COMM_WORLD, ierr)
          else
            call MPI_SEND(s(1,i,1), 1, GLOBAL_SMALL_BDx_HSLICES, 
     &         owner_of_hslice(i), i,  MPI_COMM_WORLD, ierr)
          end if
     	end do
      else if (HWorker) then
	if (no_of_hslices(myid) .eq. izloc) then
          call MPI_RECV(sloc(1,1,1), 1, LOCAL_LARGE_BDx_HSLICES, 0,
     &			myid, MPI_COMM_WORLD, status, ierr)
        else
          call MPI_RECV(sloc(1,1,1), 1, LOCAL_SMALL_BDx_HSLICES, 0,
     &			myid, MPI_COMM_WORLD, status, ierr)
	end if
      end if
c
      end subroutine  distrib_h_BDx_v2

c---------------------------------------------------------------------------
c  Master distributes the x-boundary array 
c  - careful, Y is used from
c  former different naming BDy...
c---------------------------------------------------------------------------
c      
      subroutine distrib_yh_BD_v2(iy,iz,izloc,Ns,s,sloc)

      implicit none
c
      integer :: iy, iz, Ns, izloc
      real :: s(iy,iz,2,Ns)
      real :: sloc(iy,izloc,2,Ns)
      integer :: i, ierr, status(MPI_STATUS_SIZE)
c      
      if (Master) then
        do i=1,NHworkers
	  if (no_of_hslices(i) .eq. izloc) then
            call MPI_SEND(s(1,i,1,1), 1, GLOBAL_LARGE_BD_YHSLICES, 
     &         owner_of_hslice(i), i,  MPI_COMM_WORLD, ierr)
          else
            call MPI_SEND(s(1,i,1,1), 1, GLOBAL_SMALL_BD_YHSLICES, 
     &         owner_of_hslice(i), i,  MPI_COMM_WORLD, ierr)
          end if
     	end do
      else if (HWorker) then
	if (no_of_hslices(myid) .eq. izloc) then
          call MPI_RECV(sloc(1,1,1,1), 1, LOCAL_LARGE_BD_YHSLICES, 0,
     &			myid, MPI_COMM_WORLD, status, ierr)
        else
          call MPI_RECV(sloc(1,1,1,1), 1, LOCAL_SMALL_BD_YHSLICES, 0,
     &			myid, MPI_COMM_WORLD, status, ierr)
	end if
      end if
c
      end subroutine  distrib_yh_BD_v2

c *******************************************************************
c ***************** collective communication ************************
c *******************************************************************

c-----------------------------------------------------------------
c  Master distributes a 4D array of data in h-slice 
c     format to all workers
c   using scatterv
c-----------------------------------------------------------------
      subroutine distrib_h_4D_v3(ix,iy,iz,izloc,Ns,s,s_h)
c
      implicit none
c
      integer :: ix, iy, iz, izloc, Ns
      real :: s(ix,iy,iz,Ns)
      real :: s_h(ix,iy,izloc,Ns)
      integer :: i, ierr
      integer :: sendcounts(Nprocs), displs(Nprocs)
      integer :: recvcount

c The master scatters the h-slices in two calls
      sendcounts(1) = 0
      sendcounts(2:NHworkers+1) = 1
      sendcounts(NHworkers+2:Nprocs) = 0
      displs(1:2) = 0
      do i=3, NHworkers+1
	displs(i) = displs(i-1) + ix*iy
      end do
      if (HWorker) then
        recvcount = 1
      else
        recvcount = 0
      end if
c First everybody gets the minimal size everybody has, then
c the rest processors get their further slices (max. one more)
      if (no_of_hslices(NHworkers) .eq. izloc) then
        call MPI_SCATTERV(s(1,1,1,1), sendcounts, displs, 
     &		GLOB_LARGE_4D_HSLICES, s_h(1,1,1,1),
     &		recvcount, LOCAL_LARGE_4D_HSLICES, 0,
     &		MPI_COMM_WORLD, ierr)
      else
        call MPI_SCATTERV(s(1,1,1,1), sendcounts, displs, 
     &		GLOB_SMALL_4D_HSLICES, s_h(1,1,1,1),
     &		recvcount, LOCAL_SMALL_4D_HSLICES, 0,
     &		MPI_COMM_WORLD, ierr)
	sendcounts((mod(iz,NHworkers)+2) : Nprocs) = 0
	displs((mod(iz,NHworkers)+2) : Nprocs) = 0
        if (MyId .gt. mod(iz,NHworkers)) then
          recvcount = 0
        end if
        call MPI_SCATTERV(s(1,1,NHworkers*(izloc-1)+1,1), sendcounts,  
     &		displs, GLOB_4D_HSLICE, s_h(1,1,izloc,1),
     &		recvcount, LOCAL_4D_HSLICE, 0,
     &		MPI_COMM_WORLD, ierr)
      end if

      end subroutine  distrib_h_4D_v3

c-----------------------------------------------------------------
c  Master distributes a 4D array of data in v-column format to all workers
c  using scatterv
c-----------------------------------------------------------------
      subroutine distrib_v_4D_v3(ix,iy,iz,icloc,Ns,s,s_v)
c      
      implicit none
c
      integer :: ix, iy, iz, icloc, Ns
      real :: s(ix,iy,iz,Ns)
      real :: s_v(1,icloc,iz,Ns)
      integer :: i, ierr
      integer :: sendcounts(Nprocs), displs(Nprocs)
      integer :: recvcount

c The master scatters the v-columns in two calls
      sendcounts(1) = 0
      sendcounts(2:NVworkers+1) = 1
      sendcounts(NVworkers+2:Nprocs) = 0
      displs(1:2) = 0
      do i=3, NVworkers+1
	displs(i) = displs(i-1) + 1
      end do
      if (VWorker) then
        recvcount = 1
      else
        recvcount = 0
      end if
c First everybody gets the minimal size everybody has, then
c the rest processors get their further slices (max. one more)
      if (no_of_vcols(NVworkers) .eq. icloc) then
        call MPI_SCATTERV(s(1,1,1,1), sendcounts, displs, 
     &		GLOB_LARGE_4D_VCOLS, s_v(1,1,1,1),
     &		recvcount, LOCAL_LARGE_4D_VCOLS, 0,
     &		MPI_COMM_WORLD, ierr)
      else
        call MPI_SCATTERV(s(1,1,1,1), sendcounts, displs, 
     &		GLOB_SMALL_4D_VCOLS, s_v(1,1,1,1),
     &		recvcount, LOCAL_SMALL_4D_VCOLS, 0,
     &		MPI_COMM_WORLD, ierr)
	sendcounts((mod(ix*iy,NVworkers)+2) : Nprocs) = 0
	displs((mod(ix*iy,NVworkers)+2) : Nprocs) = 0
        if (MyId .gt. mod(ix*iy,NVworkers)) then
          recvcount = 0
        end if
        call MPI_SCATTERV(s(NVworkers*(icloc-1)+1,1,1,1), sendcounts,  
     &		displs, GLOB_4D_VCOL, s_v(1,icloc,1,1),
     &		recvcount, LOCAL_4D_VCOL, 0,
     &		MPI_COMM_WORLD, ierr)
      end if

      end subroutine  distrib_v_4D_v3

c-----------------------------------------------------------------
c  Master gathers a 4D array of data in h-slice format from all workers
c-----------------------------------------------------------------
      subroutine gather_h_4D_v3(ix,iy,iz,izloc,Ns,s,s_h)
c      
      implicit none
c
      integer :: ix, iy, iz, izloc, Ns
      real   :: s(ix,iy,iz,Ns)
      real   :: s_h(ix,iy,izloc,Ns)
      integer :: i, ierr
      integer :: recvcounts(Nprocs), displs(Nprocs)
      integer :: sendcount

c The master scatters the h-slices in two calls
      recvcounts(1) = 0
      recvcounts(2:NHworkers+1) = 1
      recvcounts(NHworkers+2:Nprocs) = 0
      displs(1:2) = 0
      do i=3, NHworkers+1
	displs(i) = displs(i-1) + ix*iy
      end do
      if (HWorker) then
        sendcount = 1
      else 
        sendcount = 0
      end if
c First everybody gets the minimal size everybody has, then
c the rest processors get their further slices (max. one more)
      if (no_of_hslices(NHworkers) .eq. izloc) then
        call MPI_GATHERV(s_h(1,1,1,1), sendcount,
     &		LOCAL_LARGE_4D_HSLICES, s(1,1,1,1),
     &		recvcounts, displs, GLOB_LARGE_4D_HSLICES, 0,
     &		MPI_COMM_WORLD, ierr)
      else
        call MPI_GATHERV(s_h(1,1,1,1), sendcount,
     &		LOCAL_SMALL_4D_HSLICES, s(1,1,1,1),
     &		recvcounts, displs, GLOB_SMALL_4D_HSLICES, 0,
     &		MPI_COMM_WORLD, ierr)
	recvcounts((mod(iz,NHworkers)+2) : Nprocs) = 0
	displs((mod(iz,NHworkers)+2) : Nprocs) = 0
        if (MyId .gt. mod(iz,NHworkers)) then
          sendcount = 0
        end if
        call MPI_GATHERV(s_h(1,1,izloc,1), sendcount,
     &		LOCAL_4D_HSLICE, s(1,1,NHworkers*(izloc-1)+1,1),
     &		recvcounts, displs, GLOB_4D_HSLICE, 0,
     &		MPI_COMM_WORLD, ierr)
      end if
c
      end subroutine  gather_h_4D_v3  

c-----------------------------------------------------------------
c  Master gathers a 4D array of data in vcolumn format from all workers
c-----------------------------------------------------------------
      subroutine gather_v_4D_v3(ix,iy,iz,icloc,Ns,s,s_v)
c      
      implicit none
c
      integer :: ix, iy, iz, icloc, Ns
      real   :: s(ix,iy,iz,Ns)
      real   :: s_v(1,icloc,iz,Ns)
      integer :: i, ierr
      integer :: recvcounts(Nprocs), displs(Nprocs)
      integer :: sendcount

c The master scatters the v-columns in two calls
      recvcounts(1) = 0
      recvcounts(2:NVworkers+1) = 1
      recvcounts(NVworkers+2:Nprocs) = 0
      displs(1:2) = 0
      do i=3, NVworkers+1
	displs(i) = displs(i-1) + 1
      end do
      if (VWorker) then
        sendcount = 1
      else
        sendcount = 0
      end if
c First everybody gets the minimal size everybody has, then
c the rest processors get their further slices (max. one more)
      if (no_of_vcols(NVworkers) .eq. icloc) then
        call MPI_GATHERV(s_v(1,1,1,1), sendcount,
     &		LOCAL_LARGE_4D_VCOLS, s(1,1,1,1),
     &		recvcounts, displs, GLOB_LARGE_4D_VCOLS, 0,
     &		MPI_COMM_WORLD, ierr)
      else
        call MPI_GATHERV(s_v(1,1,1,1), sendcount,
     &		LOCAL_SMALL_4D_VCOLS, s(1,1,1,1),
     &		recvcounts, displs, GLOB_SMALL_4D_VCOLS, 0,
     &		MPI_COMM_WORLD, ierr)
	recvcounts((mod(ix*iy,NVworkers)+2) : Nprocs) = 0
	displs((mod(ix*iy,NVworkers)+2) : Nprocs) = 0
        if (MyId .gt. mod(ix*iy,NVworkers)) then
          sendcount = 0
        end if
        call MPI_GATHERV(s_v(1,icloc,1,1), sendcount,
     &		LOCAL_4D_VCOL, s(NVworkers*(icloc-1)+1,1,1,1),
     &		recvcounts, displs, GLOB_4D_VCOL, 0,
     &		MPI_COMM_WORLD, ierr)
      end if
c
      end subroutine  gather_v_4D_v3  


c---------------------------------------------------------
c      Slaves exchange data:  - old format is h-slices (for sending)
c                             - new format is v-columns (for receiving)
c 	using MPI_ALLTOALL and copying
c---------------------------------------------------------
      subroutine shuffle_h2v_4D_v3(ix,iy,iz,izloc,icloc,Ns,s_h,s_v)  
c      
      implicit none
c
      integer,intent(in) :: ix, iy, iz, izloc, icloc, Ns
      real,intent(in)    :: s_h(ix,iy,izloc,Ns)
      real,intent(out)    :: s_v(1,icloc,iz,Ns)
      integer :: i,ierr, p, lk
      integer :: status(MPI_STATUS_SIZE)
      real,pointer, dimension(:,:,:,:), save :: sbuf,rbuf

      if (HWorker .or. VWorker) then
        if (.not.allocated(sbuf)) then
          allocate( sbuf(icloc,izloc,Ns,max(NHworkers,NVworkers)), 
     &		STAT=ierr)
        endif
        if (.not.allocated(rbuf)) then
          allocate( rbuf(icloc,izloc,Ns,max(NHworkers,NVworkers)), 
     &		STAT=ierr )
        endif

      if (HWorker) then

        do i=1,ix*iy
	  lk = local_vcol_id(i)
	  sbuf(lk,1:no_of_hslices(MyId),1:Ns,owner_of_vcol(i)) = 
     &		s_h(planar_vcol_id(i,1),
     &		planar_vcol_id(i,2), 1:no_of_hslices(MyId),1:Ns)
        end do
      end if

      call MPI_ALLTOALL(sbuf, icloc*izloc*Ns, MPI_REAL, 
     &          rbuf, icloc*izloc*Ns, MPI_REAL, 
     &          MPI_COMM_WORKERS, ierr)
     
      if (VWorker) then
        do p=1,NHworkers
  	  s_v(1,1:icloc,owned_hslices(p,1:no_of_hslices(p)),1:Ns)=
     &	            rbuf(1:icloc,1:no_of_hslices(p),1:Ns,p)
        enddo
      end if	
      end if
      
      end subroutine shuffle_h2v_4D_v3

c---------------------------------------------------------
c      Slaves exchange data:  - old format is v-columns (for sending)
c                             - new format is h-slices (for receiving)
c	using MPI_ALLTOALLV and copying
c---------------------------------------------------------
      subroutine shuffle_v2h_4D_v3(ix,iy,iz,izloc,icloc,Ns,s_h,s_v)  
c      
      implicit none
c
      integer,intent(in)  :: ix, iy, iz, izloc, icloc, Ns
      real,intent(out)  :: s_h(ix,iy,izloc,Ns)
      real,intent(in)   :: s_v(1,icloc,iz,Ns)
      integer :: i, ierr, p, lk
      integer :: status(MPI_STATUS_SIZE)
      real,pointer, dimension(:,:,:,:),save    :: sbuf,rbuf

      if (HWorker .or. VWorker) then
        if (.not.allocated(sbuf)) then
          allocate( sbuf(icloc,izloc,Ns,max(NHworkers,NVworkers)),
     &		STAT=ierr )
        endif
        if (.not.allocated(rbuf)) then
          allocate( rbuf(icloc,izloc,Ns,max(NVworkers,NHworkers)),
     &		STAT=ierr )
        endif

      if (VWorker) then
        do p=1,NHworkers
	  sbuf(1:icloc,1:no_of_hslices(p),1:Ns,p)=
     &	      s_v(1,1:icloc,owned_hslices(p,1:no_of_hslices(p)),1:Ns)
        end do
      end if
	
      call MPI_ALLTOALL(sbuf, icloc*izloc*Ns, MPI_REAL, 
     &          rbuf, icloc*izloc*Ns, MPI_REAL, 
     &          MPI_COMM_WORKERS, ierr)
     
      if (HWorker) then
        do p=1,NVworkers
          do i = 1,no_of_vcols(p)
             s_h(planar_vcol_id(global_vcol_id(p,i),1),
     &		planar_vcol_id(global_vcol_id(p,i),2),
     &	     	1:no_of_hslices(MyId),1:Ns)=
     &	        rbuf(i,1:no_of_hslices(MyId),1:Ns,p)
          enddo
        end do
	
      end if
      end if

      end subroutine shuffle_v2h_4D_v3

c-----------------------------------------------------------------
c  Master distributes a 2D array of data in v-column format to all workers
c  using scatterv
c-----------------------------------------------------------------
      subroutine distrib_v_2D_v3(ix,iy,icloc,s,s_v)
c      
      implicit none
c
      integer :: ix, iy, icloc
      real :: s(ix,iy)
      real :: s_v(1,icloc)
      integer :: i, ierr
      integer :: sendcounts(Nprocs), displs(Nprocs)
      integer :: recvcount

c The master scatters the v-columns in two calls
      sendcounts(1) = 0
      sendcounts(2:NVworkers+1) = 1
      sendcounts(NVworkers+2:Nprocs) = 0
      displs(1:2) = 0
      do i=3, NVworkers+1
	displs(i) = displs(i-1) + 1
      end do
      if (VWorker) then
        recvcount = 1
      else
        recvcount = 0
      end if
c First everybody gets the minimal size everybody has, then
c the rest processors get their further slices (max. one more)
      if (no_of_vcols(NVworkers) .eq. icloc) then
        call MPI_SCATTERV(s(1,1), sendcounts, displs, 
     &		GLOB_LARGE_2D_VCOLS, s_v(1,1),
     &		recvcount, LOCAL_LARGE_2D_VCOLS, 0,
     &		MPI_COMM_WORLD, ierr)
      else
        call MPI_SCATTERV(s(1,1), sendcounts, displs, 
     &		GLOB_SMALL_2D_VCOLS, s_v(1,1),
     &		recvcount, LOCAL_SMALL_2D_VCOLS, 0,
     &		MPI_COMM_WORLD, ierr)
	sendcounts((mod(ix*iy,NVworkers)+2) : Nprocs) = 0
	displs((mod(ix*iy,NVworkers)+2) : Nprocs) = 0
        if (MyId .gt. mod(ix*iy,NVworkers)) then
          recvcount = 0
        end if
        call MPI_SCATTERV(s(NVworkers*(icloc-1)+1,1), sendcounts,  
     &		displs, MPI_REAL, s_v(1,icloc),
     &		recvcount, MPI_REAL, 0,
     &		MPI_COMM_WORLD, ierr)
      end if

      end subroutine  distrib_v_2D_v3

c-----------------------------------------------------------------
c  Master distributes a BDz array of data in v-column format to all workers
c  using scatterv
c-----------------------------------------------------------------
      subroutine distrib_v_BDz_v3(ix,iy,icloc,Ns,s,s_v)
c      
      implicit none
c
      integer :: ix, iy, icloc, Ns
      real :: s(ix,iy,Ns)
      real :: s_v(1,icloc,Ns)
      integer :: i, ierr
      integer :: sendcounts(Nprocs), displs(Nprocs)
      integer :: recvcount

c The master scatters the v-columns in two calls
      sendcounts(1) = 0
      sendcounts(2:NVworkers+1) = 1
      sendcounts(NVworkers+2:Nprocs) = 0
      displs(1:2) = 0
      do i=3, NVworkers+1
	displs(i) = displs(i-1) + 1
      end do
      if (Vworker) then
        recvcount = 1
      else
        recvcount = 0
      end if
c First everybody gets the minimal size everybody has, then
c the rest processors get their further slices (max. one more)
      if (no_of_vcols(NVworkers) .eq. icloc) then
        call MPI_SCATTERV(s(1,1,1), sendcounts, displs, 
     &		GLOB_LARGE_BDz_VCOLS, s_v(1,1,1),
     &		recvcount, LOCAL_LARGE_BDz_VCOLS, 0,
     &		MPI_COMM_WORLD, ierr)
      else
        call MPI_SCATTERV(s(1,1,1), sendcounts, displs, 
     &		GLOB_SMALL_BDz_VCOLS, s_v(1,1,1),
     &		recvcount, LOCAL_SMALL_BDz_VCOLS, 0,
     &		MPI_COMM_WORLD, ierr)
	sendcounts((mod(ix*iy,NVworkers)+2) : Nprocs) = 0
	displs((mod(ix*iy,NVworkers)+2) : Nprocs) = 0
        if (MyId .gt. mod(ix*iy,NVworkers)) then
          recvcount = 0
        end if
        call MPI_SCATTERV(s(NVworkers*(icloc-1)+1,1,1), sendcounts,  
     &		displs, GLOB_BDz_VCOL, s_v(1,icloc,1),
     &		recvcount, LOCAL_BDz_VCOL, 0,
     &		MPI_COMM_WORLD, ierr)
      end if

      end subroutine  distrib_v_BDz_v3

c-----------------------------------------------------------------
c  Master distributes a 3D array of data in h-slice 
c     format to all workers
c   using scatterv
c-----------------------------------------------------------------
      subroutine distrib_h_3D_v3(ix,iy,iz,izloc,s,s_h)
c
      implicit none
c
      integer :: ix, iy, iz, izloc
      real :: s(ix,iy,iz)
      real :: s_h(ix,iy,izloc)
      integer :: i, ierr
      integer :: sendcounts(Nprocs), displs(Nprocs)
      integer :: recvcount

c The master scatters the h-slices in two calls
      sendcounts(1) = 0
      sendcounts(2:NHworkers+1) = 1
      sendcounts(NHworkers+2:Nprocs) = 0
      displs(1:2) = 0
      do i=3, NHworkers+1
	displs(i) = displs(i-1) + ix*iy
      end do
      if (HWorker) then
        recvcount = 1
      else
        recvcount = 0
      end if
c First everybody gets the minimal size everybody has, then
c the rest processors get their further slices (max. one more)
      if (no_of_hslices(NHworkers) .eq. izloc) then
        call MPI_SCATTERV(s(1,1,1), sendcounts, displs, 
     &		GLOB_LARGE_3D_HSLICES, s_h(1,1,1),
     &		recvcount, LOCAL_LARGE_3D_HSLICES, 0,
     &		MPI_COMM_WORLD, ierr)
      else
        call MPI_SCATTERV(s(1,1,1), sendcounts, displs, 
     &		GLOB_SMALL_3D_HSLICES, s_h(1,1,1),
     &		recvcount, LOCAL_SMALL_3D_HSLICES, 0,
     &		MPI_COMM_WORLD, ierr)
	sendcounts(2:(mod(iz,NHworkers)+1)) = ix*iy
	sendcounts((mod(iz,NHworkers)+2) : Nprocs) = 0
	displs((mod(iz,NHworkers)+2) : Nprocs) = 0
	if (.not.Master) recvcount = ix*iy
        if (MyId .gt. mod(iz,NHworkers)) then
          recvcount = 0
        end if
        call MPI_SCATTERV(s(1,1,NHworkers*(izloc-1)+1), sendcounts,  
     &		displs, MPI_REAL, s_h(1,1,izloc),
     &		recvcount, MPI_REAL, 0,
     &		MPI_COMM_WORLD, ierr)
      end if

      end subroutine  distrib_h_3D_v3

c-----------------------------------------------------------------
c  Master distributes a 3D array of data in v-column format to all workers
c  using scatterv
c-----------------------------------------------------------------
      subroutine distrib_v_3D_v3(ix,iy,iz,icloc,s,s_v)
c      
      implicit none
c
      integer :: ix, iy, iz, icloc
      real :: s(ix,iy,iz)
      real :: s_v(1,icloc,iz)
      integer :: i, ierr
      integer :: sendcounts(Nprocs), displs(Nprocs)
      integer :: recvcount

c The master scatters the v-columns in two calls
      sendcounts(1) = 0
      sendcounts(2:NVworkers+1) = 1
      sendcounts(NVworkers+2:Nprocs) = 0
      displs(1:2) = 0
      do i=3, NVworkers+1
	displs(i) = displs(i-1) + 1
      end do
      if (VWorker) then
        recvcount = 1
      else
        recvcount = 0
      end if
c First everybody gets the minimal size everybody has, then
c the rest processors get their further slices (max. one more)
      if (no_of_vcols(NVworkers) .eq. icloc) then
        call MPI_SCATTERV(s(1,1,1), sendcounts, displs, 
     &		GLOB_LARGE_3D_VCOLS, s_v(1,1,1),
     &		recvcount, LOCAL_LARGE_3D_VCOLS, 0,
     &		MPI_COMM_WORLD, ierr)
      else
        call MPI_SCATTERV(s(1,1,1), sendcounts, displs, 
     &		GLOB_SMALL_3D_VCOLS, s_v(1,1,1),
     &		recvcount, LOCAL_SMALL_3D_VCOLS, 0,
     &		MPI_COMM_WORLD, ierr)
	sendcounts((mod(ix*iy,NVworkers)+2) : Nprocs) = 0
	displs((mod(ix*iy,NVworkers)+2) : Nprocs) = 0
        if (MyId .gt. mod(ix*iy,NVworkers)) then
          recvcount = 0
        end if
        call MPI_SCATTERV(s(NVworkers*(icloc-1)+1,1,1), sendcounts,  
     &		displs, GLOB_3D_VCOL, s_v(1,icloc,1),
     &		recvcount, LOCAL_3D_VCOL, 0,
     &		MPI_COMM_WORLD, ierr)
      end if

      end subroutine  distrib_v_3D_v3

c-----------------------------------------------------------------
c  Master distributes a BDy array of data in h-slice 
c     format to all workers
c   using scatterv
c-----------------------------------------------------------------
      subroutine distrib_h_BDy_v3(ix,iz,izloc,Ns,s,s_h)
c
      implicit none
c
      integer :: ix, iy, iz, izloc, Ns
      real :: s(ix,iz,Ns)
      real :: s_h(ix,izloc,Ns)
      integer :: i, ierr
      integer :: sendcounts(Nprocs), displs(Nprocs)
      integer :: recvcount

c The master scatters the h-slices in two calls
      sendcounts(1) = 0
      sendcounts(2:NHworkers+1) = 1
      sendcounts(NHworkers+2:Nprocs) = 0
      displs(1:2) = 0
      do i=3, NHworkers+1
	displs(i) = displs(i-1) + ix
      end do
      if (Hworker) then
        recvcount = 1
      else
        recvcount = 0
      end if
c First everybody gets the minimal size everybody has, then
c the rest processors get their further slices (max. one more)
      if (no_of_hslices(NHworkers) .eq. izloc) then
        call MPI_SCATTERV(s(1,1,1), sendcounts, displs, 
     &		GLOB_LARGE_BDy_HSLICES, s_h(1,1,1),
     &		recvcount, LOCAL_LARGE_BDy_HSLICES, 0,
     &		MPI_COMM_WORLD, ierr)
      else
        call MPI_SCATTERV(s(1,1,1), sendcounts, displs, 
     &		GLOB_SMALL_BDy_HSLICES, s_h(1,1,1),
     &		recvcount, LOCAL_SMALL_BDy_HSLICES, 0,
     &		MPI_COMM_WORLD, ierr)
	sendcounts((mod(iz,NHworkers)+2) : Nprocs) = 0
	displs((mod(iz,NHworkers)+2) : Nprocs) = 0
        if (MyId .gt. mod(iz,NHworkers)) then
          recvcount = 0
        end if
        call MPI_SCATTERV(s(1,NHworkers*(izloc-1)+1,1), sendcounts,  
     &		displs, GLOB_BDy_HSLICE, s_h(1,izloc,1),
     &		recvcount, LOCAL_BDy_HSLICE, 0,
     &		MPI_COMM_WORLD, ierr)
      end if

      end subroutine  distrib_h_BDy_v3

c-----------------------------------------------------------------
c  Master distributes a BDy array of data in h-slice 
c     format to all workers
c   using scatterv
c  - careful, X is used from
c  former different naming BDx...
c-----------------------------------------------------------------
      subroutine distrib_xh_BD_v3(ix,iz,izloc,Ns,s,s_h)
c
      implicit none
c
      integer :: ix, iy, iz, izloc, Ns
      real :: s(ix,iz,2,Ns)
      real :: s_h(ix,izloc,2,Ns)
      integer :: i, ierr
      integer :: sendcounts(Nprocs), displs(Nprocs)
      integer :: recvcount

c The master scatters the h-slices in two calls
      sendcounts(1) = 0
      sendcounts(2:NHworkers+1) = 1
      sendcounts(NHworkers+2:Nprocs) = 0
      displs(1:2) = 0
      do i=3, NHworkers+1
	displs(i) = displs(i-1) + ix
      end do
      if (Hworker) then
        recvcount = 1
      else
        recvcount = 0
      end if
c First everybody gets the minimal size everybody has, then
c the rest processors get their further slices (max. one more)
      if (no_of_hslices(NHworkers) .eq. izloc) then
        call MPI_SCATTERV(s(1,1,1,1), sendcounts, displs, 
     &		GLOB_LARGE_BD_XHSLICES, s_h(1,1,1,1),
     &		recvcount, LOCAL_LARGE_BD_XHSLICES, 0,
     &		MPI_COMM_WORLD, ierr)
      else
        call MPI_SCATTERV(s(1,1,1,1), sendcounts, displs, 
     &		GLOB_SMALL_BD_XHSLICES, s_h(1,1,1,1),
     &		recvcount, LOCAL_SMALL_BD_XHSLICES, 0,
     &		MPI_COMM_WORLD, ierr)
	sendcounts((mod(iz,NHworkers)+2) : Nprocs) = 0
	displs((mod(iz,NHworkers)+2) : Nprocs) = 0
        if (MyId .gt. mod(iz,NHworkers)) then
          recvcount = 0
        end if
        call MPI_SCATTERV(s(1,NHworkers*(izloc-1)+1,1,1), sendcounts,  
     &		displs, GLOB_BD_XHSLICE, s_h(1,izloc,1,1),
     &		recvcount, LOCAL_BD_XHSLICE, 0,
     &		MPI_COMM_WORLD, ierr)
      end if

      end subroutine  distrib_xh_BD_v3

c-----------------------------------------------------------------
c  Master distributes a BDx array of data in h-slice 
c     format to all workers
c   using scatterv
c-----------------------------------------------------------------
      subroutine distrib_h_BDx_v3(iy,iz,izloc,Ns,s,s_h)
c
      implicit none
c
      integer :: ix, iy, iz, izloc, Ns
      real :: s(iy,iz,Ns)
      real :: s_h(iy,izloc,Ns)
      integer :: i, ierr
      integer :: sendcounts(Nprocs), displs(Nprocs)
      integer :: recvcount

c The master scatters the h-slices in two calls
      sendcounts(1) = 0
      sendcounts(2:NHworkers+1) = 1
      sendcounts(NHworkers+2:Nprocs) = 0
      displs(1:2) = 0
      do i=3, NHworkers+1
	displs(i) = displs(i-1) + iy
      end do
      if (HWorker) then
        recvcount = 1
      else
        recvcount = 0
      end if
c First everybody gets the minimal size everybody has, then
c the rest processors get their further slices (max. one more)
      if (no_of_hslices(NHworkers) .eq. izloc) then
        call MPI_SCATTERV(s(1,1,1), sendcounts, displs, 
     &		GLOB_LARGE_BDx_HSLICES, s_h(1,1,1),
     &		recvcount, LOCAL_LARGE_BDx_HSLICES, 0,
     &		MPI_COMM_WORLD, ierr)
      else
        call MPI_SCATTERV(s(1,1,1), sendcounts, displs, 
     &		GLOB_SMALL_BDx_HSLICES, s_h(1,1,1),
     &		recvcount, LOCAL_SMALL_BDx_HSLICES, 0,
     &		MPI_COMM_WORLD, ierr)
	sendcounts((mod(iz,NHworkers)+2) : Nprocs) = 0
	displs((mod(iz,NHworkers)+2) : Nprocs) = 0
        if (MyId .gt. mod(iz,NHworkers)) then
          recvcount = 0
        end if
        call MPI_SCATTERV(s(1,NHworkers*(izloc-1)+1,1), sendcounts,  
     &		displs, GLOB_BDx_HSLICE, s_h(1,izloc,1),
     &		recvcount, LOCAL_BDx_HSLICE, 0,
     &		MPI_COMM_WORLD, ierr)
      end if

      end subroutine  distrib_h_BDx_v3

c-----------------------------------------------------------------
c  Master distributes a BDx array of data in h-slice 
c     format to all workers
c   using scatterv
c  - careful, Y is used from
c  former different naming BDy...
c-----------------------------------------------------------------
      subroutine distrib_yh_BD_v3(iy,iz,izloc,Ns,s,s_h)
c
      implicit none
c
      integer :: ix, iy, iz, izloc, Ns
      real :: s(iy,iz,2,Ns)
      real :: s_h(iy,izloc,2,Ns)
      integer :: i, ierr
      integer :: sendcounts(Nprocs), displs(Nprocs)
      integer :: recvcount

c The master scatters the h-slices in two calls
      sendcounts(1) = 0
      sendcounts(2:NHworkers+1) = 1
      sendcounts(NHworkers+2:Nprocs) = 0
      displs(1:2) = 0
      do i=3, NHworkers+1
	displs(i) = displs(i-1) + iy
      end do
      if (HWorker) then
        recvcount = 1
      else
        recvcount = 0
      end if
c First everybody gets the minimal size everybody has, then
c the rest processors get their further slices (max. one more)
      if (no_of_hslices(NHworkers) .eq. izloc) then
        call MPI_SCATTERV(s(1,1,1,1), sendcounts, displs, 
     &		GLOB_LARGE_BD_YHSLICES, s_h(1,1,1,1),
     &		recvcount, LOCAL_LARGE_BD_YHSLICES, 0,
     &		MPI_COMM_WORLD, ierr)
      else
        call MPI_SCATTERV(s(1,1,1,1), sendcounts, displs, 
     &		GLOB_SMALL_BD_YHSLICES, s_h(1,1,1,1),
     &		recvcount, LOCAL_SMALL_BD_YHSLICES, 0,
     &		MPI_COMM_WORLD, ierr)
	sendcounts((mod(iz,NHworkers)+2) : Nprocs) = 0
	displs((mod(iz,NHworkers)+2) : Nprocs) = 0
        if (MyId .gt. mod(iz,NHworkers)) then
          recvcount = 0
        end if
        call MPI_SCATTERV(s(1,NHworkers*(izloc-1)+1,1,1), sendcounts,  
     &		displs, GLOB_BD_YHSLICE, s_h(1,izloc,1,1),
     &		recvcount, LOCAL_BD_YHSLICE, 0,
     &		MPI_COMM_WORLD, ierr)
      end if

      end subroutine  distrib_yh_BD_v3


      end module HVCommunicationLibrary

