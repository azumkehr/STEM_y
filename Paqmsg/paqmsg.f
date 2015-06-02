c---------------------------------------------------------------------------------------
c
c  The following function defines the mapping of the global
c  domain to local workers
c
c      owner_of_xslice(1..iy) = the process which gets that x-slice
c      owner_of_yslice(1..ix) = the process which gets that y-slice
c
c      no_of_xslices(p),  1<=p<=NXworkers = the number of xslices	 
c                         owned by process p
c      no_of_yslices(p),  1<=p<=NYworkers = the number of yslices	 
c                         owned by process p
c
c      local_xslice_id(1..iy) = local index of global xslice 1..iy
c                           (index seen by the process which owns it)
c      local_yslice_id(1..ix) = local index of global yslice 1..ix
c                           (index seen by the process which owns it)
c
c      global_xslice_id(p,i) 1<=p<=NXworkers, 1<=i<=iyloc
c               xslice of local index i on process p is the global slice ...
c      global_yslice_id(1..NYworkers,1..ixloc)
c               yslice of local index i on process p is the global slice ...
c
c      owned_xslices(p,i) = 1..iy, 1<=p<=NXworkers, 1<=i<=iyloc 
c                which global x-slices i are owned by slave p   
c      owned_yslices(p,i)= 1..ix, 1<=p<=NYworkers, 1<=i<=ixloc 
c                which y-slices i are owned by slave p   
c
c---------------------------------------------------------------------------------------
 
      module XYParallelDataMap
c       
c  No. of processes
        integer :: Nprocs       
c  No. of slave processes for x-slices and y-slices, respectively
        integer :: NXworkers       
        integer :: NYworkers       
c  Current process id
        integer :: MyId
c       
        logical :: Master, XWorker, YWorker

        integer, pointer, dimension(:)   :: no_of_xslices
	integer, pointer, dimension(:)   :: no_of_yslices 
        integer, pointer, dimension(:,:) :: global_xslice_id
        integer, pointer, dimension(:,:) :: global_yslice_id
        integer, pointer, dimension(:)   :: local_xslice_id
        integer, pointer, dimension(:)   :: local_yslice_id
        integer, pointer, dimension(:)   :: owner_of_xslice
	integer, pointer, dimension(:)   :: owner_of_yslice
	integer, pointer, dimension(:,:) :: owned_xslices
	integer, pointer, dimension(:,:) :: owned_yslices
c             
      contains 
c       
      subroutine CreateMap(ix, iy, ixloc, iyloc) 

      integer :: ix, iy, ixloc, iyloc
      integer :: i, id, p
      integer :: modx, mody

      Master   = MyId.eq.0
c     There cannot be more Xworkers than x-slices
      XWorker  = MyId.ne.0 .and. MyId.le.iy
c     There cannot be more Yworkers than y-slices
      YWorker  = MyId.ne.0 .and. MyId.le.ix

      if (Nprocs-1 .le. iy) then
        NXworkers = Nprocs-1
      else
        NXworkers = iy
      end if
      if (Nprocs-1 .le. ix) then
        NYworkers = Nprocs-1
      else
        NYworkers = ix
      end if

c The dimensions of the local data sets
      modx = mod(ix, NYworkers)
      if ( modx .eq. 0 ) then
         ixloc = ix/NYworkers
      else
         ixloc = (ix-modx)/NYworkers + 1
      endif	 	 
c
      mody = mod(iy, NXworkers)
      if ( mody .eq. 0 ) then
         iyloc = iy/NXworkers
      else
         iyloc = (iy-mody)/NXworkers + 1
      endif	

c Allocate stuff
      allocate( no_of_xslices(NXworkers),STAT=ierr )
      allocate( no_of_yslices(NYworkers),STAT=ierr )      
      allocate( global_xslice_id(NXworkers,iyloc),STAT=ierr ) 	 
      allocate( global_yslice_id(NYworkers,ixloc),STAT=ierr ) 	 
      allocate( local_xslice_id(iy),STAT=ierr )
      allocate( local_yslice_id(ix),STAT=ierr )
      allocate( owner_of_xslice(iy),STAT=ierr )
      allocate( owner_of_yslice(ix),STAT=ierr )
      allocate( owned_xslices(NXworkers,iyloc),STAT=ierr )
      allocate( owned_yslices(NYworkers,ixloc),STAT=ierr )
      
c Distribution of the x and y slices is set here - cyclic distribution
      do i=1,iy
         owner_of_xslice(i)=mod(i-1,NXworkers)+1
      enddo
      do i=1,ix
         owner_of_yslice(i)=mod(i-1,NYworkers)+1
      enddo

c Depending on distribution we set the global and local position      
c first for y
      do p=1,NYworkers
      id = 0
      do i=1,ix
         if ( owner_of_yslice(i) .eq. p) then
	   id = id+1
	    global_yslice_id(p,id) = i
	    local_yslice_id(i)     = id
	 endif
      enddo
      if (id .gt. ixloc) then
        print*,' Process ',p,' is given ',id,
     &          ' yslices, but ixloc = ',ixloc
        stop
      else
         no_of_yslices(p) = id
      endif	 
      enddo	
c then for x        
      do p=1,NXworkers
      id = 0
      do i=1,iy
         if ( owner_of_xslice(i) .eq. p) then
	   id = id+1
	    global_xslice_id(p,id) = i
	    local_xslice_id(i)     = id
	 endif
      enddo
      if (id .gt. iyloc) then
        print*,' Process ',p,' is given ',id,
     &          ' xslices, but iyloc = ',iyloc
        stop
      else
         no_of_xslices(p) = id
      endif	 
      enddo

c owned_array depends on the distribution of the slices
      do p=1,NXworkers
        do i=1,no_of_xslices(p)
          owned_xslices(p,i) = p+(i-1)*NXworkers
        enddo
        do i=no_of_xslices(p)+1,iyloc
          owned_xslices(p,i) = 0
        enddo
      enddo
      do p=1,NYworkers
        do i=1,no_of_yslices(p)
          owned_yslices(p,i) = p+(i-1)*NYworkers
        enddo
        do i=no_of_yslices(p)+1,ixloc
          owned_yslices(p,i) = 0
        enddo
      enddo
c
      end subroutine CreateMap
c      
      end module XYParallelDataMap

 
c-------------------------------------------------------------------------------
C Handles for special communication data types
c-------------------------------------------------------------------------------
      module XYCommDataTypes
      
        integer :: GROUP_WORKERS, MPI_COMM_WORKERS
c	for standard versions of X and Y distribution and gathering
c	for each type of X and Y-slice: 2D, 2DN, 3D, BD, 4D, local, global
        integer :: GLOBAL_4D_XSLICE, LOCAL_4D_XSLICE
        integer :: GLOBAL_4D_YSLICE, LOCAL_4D_YSLICE
        integer :: GLOBAL_4D_COLUMN
        integer :: LOCAL_4D_XCOLUMN, LOCAL_4D_YCOLUMN
	integer :: GLOBAL_2D_YSLICE, LOCAL_2D_YSLICE
	integer :: GLOBAL_2DN_XSLICE, LOCAL_2DN_XSLICE
	integer :: GLOBAL_2DN_YSLICE, LOCAL_2DN_YSLICE
	integer :: GLOBAL_3D_XSLICE, LOCAL_3D_XSLICE
	integer :: GLOBAL_3D_YSLICE, LOCAL_3D_YSLICE
	integer :: GLOBAL_BD_XSLICE, LOCAL_BD_XSLICE
	integer :: GLOBAL_BD_YSLICE, LOCAL_BD_YSLICE
c	for single versions of X and Y distribution and gathering:
c	- large means the local x/y-slice array is completely filled
c	- small means there is one x/y-slice less in the array than
c	    in others
	integer :: GLOBAL_LARGE_4D_XSLICES, LOCAL_LARGE_4D_XSLICES
	integer :: GLOBAL_SMALL_4D_XSLICES, LOCAL_SMALL_4D_XSLICES
	integer :: GLOBAL_LARGE_4D_YSLICES, LOCAL_LARGE_4D_YSLICES
	integer :: GLOBAL_SMALL_4D_YSLICES, LOCAL_SMALL_4D_YSLICES
	integer :: GLOBAL_LARGE_2D_XSLICES, LOCAL_LARGE_2D_XSLICES
	integer :: GLOBAL_SMALL_2D_XSLICES, LOCAL_SMALL_2D_XSLICES
	integer :: GLOBAL_LARGE_2D_YSLICES, LOCAL_LARGE_2D_YSLICES
	integer :: GLOBAL_SMALL_2D_YSLICES, LOCAL_SMALL_2D_YSLICES
	integer :: GLOBAL_LARGE_2DN_XSLICES, LOCAL_LARGE_2DN_XSLICES
	integer :: GLOBAL_SMALL_2DN_XSLICES, LOCAL_SMALL_2DN_XSLICES
	integer :: GLOBAL_LARGE_2DN_YSLICES, LOCAL_LARGE_2DN_YSLICES
	integer :: GLOBAL_SMALL_2DN_YSLICES, LOCAL_SMALL_2DN_YSLICES
	integer :: GLOBAL_LARGE_3D_XSLICES, LOCAL_LARGE_3D_XSLICES
	integer :: GLOBAL_SMALL_3D_XSLICES, LOCAL_SMALL_3D_XSLICES
	integer :: GLOBAL_LARGE_3D_YSLICES, LOCAL_LARGE_3D_YSLICES
	integer :: GLOBAL_SMALL_3D_YSLICES, LOCAL_SMALL_3D_YSLICES
	integer :: GLOBAL_LARGE_BD_XSLICES, LOCAL_LARGE_BD_XSLICES
	integer :: GLOBAL_SMALL_BD_XSLICES, LOCAL_SMALL_BD_XSLICES
	integer :: GLOBAL_LARGE_BD_YSLICES, LOCAL_LARGE_BD_YSLICES
	integer :: GLOBAL_SMALL_BD_YSLICES, LOCAL_SMALL_BD_YSLICES
c	for single versions of x2y and y2x shuffeling:
c	- LX: large in the x-slice array part (full)
c	- SX: small in the x-slice array part (one xslice less)
c	- LY: large in the y-slice array part (full)
c	- SY: small in the y-slice array part (one yslice less)
c	- all possibilities for sending and receiving are combined
	integer :: LOCAL_4D_LXLY_XCOLS, LOCAL_4D_LXSY_XCOLS
	integer :: LOCAL_4D_SXLY_XCOLS, LOCAL_4D_SXSY_XCOLS
	integer :: LOCAL_4D_LXLY_YCOLS, LOCAL_4D_LXSY_YCOLS
	integer :: LOCAL_4D_SXLY_YCOLS, LOCAL_4D_SXSY_YCOLS
c	for collective communication, types used for standard and
c	single communication with reset extents
	integer :: GLOB_LARGE_4D_XSLICES, GLOB_SMALL_4D_XSLICES
	integer :: GLOB_4D_XSLICE, GLOB_4D_YSLICE
	integer :: GLOB_LARGE_4D_YSLICES, GLOB_SMALL_4D_YSLICES
	integer :: GLOB_LARGE_2D_XSLICES, GLOB_SMALL_2D_XSLICES
	integer :: GLOB_2D_YSLICE
	integer :: GLOB_LARGE_2D_YSLICES, GLOB_SMALL_2D_YSLICES
	integer :: GLOB_LARGE_2DN_XSLICES, GLOB_SMALL_2DN_XSLICES
	integer :: GLOB_2DN_XSLICE, GLOB_2DN_YSLICE
	integer :: GLOB_LARGE_2DN_YSLICES, GLOB_SMALL_2DN_YSLICES
	integer :: GLOB_LARGE_3D_XSLICES, GLOB_SMALL_3D_XSLICES
	integer :: GLOB_3D_XSLICE, GLOB_3D_YSLICE
	integer :: GLOB_LARGE_3D_YSLICES, GLOB_SMALL_3D_YSLICES
	integer :: GLOB_LARGE_BD_XSLICES, GLOB_SMALL_BD_XSLICES
	integer :: GLOB_BD_XSLICE, GLOB_BD_YSLICE
	integer :: GLOB_LARGE_BD_YSLICES, GLOB_SMALL_BD_YSLICES

      contains
      
C Creates special communication data types
      subroutine CreateCommDataTypes( ix,iy,iz,ixloc,iyloc,is)
      use XYParallelDataMap
      implicit none 
      include 'mpif.h'
c      Type Vectors to help create the global data structures
c       for x/y-slice arrays locally, globally, large and small
c     - D: dimensions of the array already included (1,2,3 of 4)
c     - G: global structure, going through ix and iy respectively
c     - L: local structure, going through ixloc and iyloc only
c     - X: xslice structure, using iyloc
c     - Y: yslice structure, using ixloc
c     - l: large datastructure (ixloc or iyloc)
c     - s: small datastructure (ixloc-1 or iyloc-1)
c     - combinations of all the above needed can be found
      integer :: TWO_DGXl, TWO_DLXl, THREE_DGXl, THREE_DLXl
      integer :: TWO_DGXs, TWO_DLXs, THREE_DGXs, THREE_DLXs
      integer :: ONE_DGYl, ONE_DLYl, TWO_DGYl, TWO_DLYl 
      integer :: THREE_DGYl, THREE_DLYl
      integer :: ONE_DGYs, ONE_DLYs, TWO_DGYs, TWO_DLYs
      integer :: THREE_DGYs, THREE_DLYs
c     Type Vectors to help create local structures for x2y
c      and y2x shuffling, only local because master doesn't
c      take part in these communications
c     - LX/LY: large x/y-slice array, completely used
c     - SX/SY: small x/y-slice array, one column garbage
c     - 1D/../3D: consisting of 1 to 3 dimensions of the 4D array
c     - combinations of the above
      integer :: LOCAL_3D_LY_XCOLUMN, LOCAL_3D_SY_XCOLUMN
      integer :: LOCAL_3D_LY_YCOLUMN, LOCAL_3D_SY_YCOLUMN
      integer :: LOCAL_1D_LY_PART_XCOL, LOCAL_1D_SY_PART_XCOL
      integer :: LOCAL_2D_LXLY_PART_XCOL, LOCAL_2D_LXSY_PART_XCOL
      integer :: LOCAL_2D_SXLY_PART_XCOL, LOCAL_2D_SXSY_PART_XCOL
      integer :: LOCAL_3D_LXLY_PART_XCOL, LOCAL_3D_LXSY_PART_XCOL
      integer :: LOCAL_3D_SXLY_PART_XCOL, LOCAL_3D_SXSY_PART_XCOL
      integer :: LOCAL_2D_LXLY_PART_YCOL, LOCAL_2D_LXSY_PART_YCOL
      integer :: LOCAL_2D_SXLY_PART_YCOL, LOCAL_2D_SXSY_PART_YCOL
      integer :: LOCAL_3D_LXLY_PART_YCOL, LOCAL_3D_LXSY_PART_YCOL
      integer :: LOCAL_3D_SXLY_PART_YCOL, LOCAL_3D_SXSY_PART_YCOL

      integer :: ix,iy,iz,ixloc,iyloc,is
      integer :: i, ierr
      integer :: sizeofreal,lb
      integer :: idle, Nw, wid, GROUP_WORLD
      integer, pointer, dimension (:) :: ranks
      integer :: blength(2), types(2), tdispls(2)

c Set sizeofreal to the byte-size of a real number
      call MPI_TYPE_EXTENT(MPI_REAL, sizeofreal,ierr)
c      print*, 'lb = ',lb,' sizeofreal =', sizeofreal      
c Set the struct parameters in order to do the resizing of data types
      types(2) = MPI_UB
      tdispls(1) = 0
      tdispls(2) = sizeofreal
      blength(1:2) = 1
      
c-------------------------------------------------------------------------------
c    Create workers group and workers communicator
c-------------------------------------------------------------------------------
      call MPI_COMM_GROUP(MPI_COMM_WORLD, GROUP_WORLD, ierr)
      idle = Nprocs - max(NXworkers,NYworkers)
      allocate(ranks(idle), STAT=ierr)
      ranks(1) = 0
      if (idle .gt. 1) then
        do i=2, idle
	  ranks(i) = i - 1 + max(NXworkers,NYworkers)
	  if (Master) print*, 'Idle processor:',ranks(i)
	end do
      end if
      call MPI_GROUP_EXCL(GROUP_WORLD,idle,ranks,GROUP_WORKERS, ierr)
      call MPI_COMM_CREATE(MPI_COMM_WORLD,
     &              GROUP_WORKERS,MPI_COMM_WORKERS,ierr) 

      if (XWorker .or. YWorker) then
        call MPI_COMM_SIZE( MPI_COMM_WORKERS, Nw, ierr )
        call MPI_COMM_RANK( MPI_COMM_WORKERS, wid, ierr )
c        print *, 'Process ', myid, ' of ', Nprocs, ' is worker ',
c     &          wid,' of ',Nw
      endif

c-------------------------------------------------------------------------------
c    Define 4D column data types
c-------------------------------------------------------------------------------
      call MPI_TYPE_VECTOR(iz*is, 1, ix*iy, MPI_REAL, 
     &                      GLOBAL_4D_COLUMN, ierr)
      call MPI_TYPE_COMMIT(GLOBAL_4D_COLUMN, ierr)
      call MPI_TYPE_VECTOR(iz*is, 1, ix*iyloc, MPI_REAL, 
     &                      LOCAL_4D_XCOLUMN, ierr)
      call MPI_TYPE_COMMIT(LOCAL_4D_XCOLUMN, ierr)
      call MPI_TYPE_VECTOR(iz*is, 1, ixloc*iy, MPI_REAL, 
     &                      LOCAL_4D_YCOLUMN, ierr)
      call MPI_TYPE_COMMIT(LOCAL_4D_YCOLUMN, ierr)
c
c-------------------------------------------------------------------------------
c    Define 4D xslice data types
c-------------------------------------------------------------------------------
      call MPI_TYPE_VECTOR(iz*is, ix, ix*iy, MPI_REAL, 
     &                      GLOBAL_4D_XSLICE, ierr)
      call MPI_TYPE_COMMIT(GLOBAL_4D_XSLICE, ierr)
      
      call MPI_TYPE_VECTOR(iz*is, ix, ix*iyloc, MPI_REAL, 
     &                      LOCAL_4D_XSLICE, ierr)
      call MPI_TYPE_COMMIT(LOCAL_4D_XSLICE, ierr)

      types(1) = GLOBAL_4D_XSLICE
      call MPI_TYPE_STRUCT(2,blength,tdispls,
     &		types, GLOB_4D_XSLICE, ierr)
      call MPI_TYPE_COMMIT(GLOB_4D_XSLICE, ierr)
      
c
c-------------------------------------------------------------------------------
c    Define 4D yslice data types
c-------------------------------------------------------------------------------
      call MPI_TYPE_VECTOR(iy*iz*is, 1, ix, MPI_REAL, 
     &                      GLOBAL_4D_YSLICE, ierr)
      call MPI_TYPE_COMMIT(GLOBAL_4D_YSLICE, ierr)
      call MPI_TYPE_VECTOR(iy*iz*is, 1, ixloc, MPI_REAL, 
     &                      LOCAL_4D_YSLICE, ierr)
      call MPI_TYPE_COMMIT(LOCAL_4D_YSLICE, ierr)

      types(1) = GLOBAL_4D_YSLICE
      call MPI_TYPE_STRUCT(2,blength,tdispls,types,
     &		GLOB_4D_YSLICE, ierr)
      call MPI_TYPE_COMMIT(GLOB_4D_YSLICE, ierr)

c-------------------------------------------------------------------------------
c    Define 4D yslice array data types 
c     (large and small, differing in one yslice element)
c-------------------------------------------------------------------------------

c This version creates the array of yslices using the given
c  yslice structures. It takes more time during communication.
	if (.false.) then
      call MPI_TYPE_HVECTOR(ixloc, 1,
     &	NYworkers*sizeofreal, 
     &	GLOBAL_4D_YSLICE, GLOBAL_LARGE_4D_YSLICES, ierr)
      call MPI_TYPE_COMMIT(GLOBAL_LARGE_4D_YSLICES, ierr)
      call MPI_TYPE_HVECTOR(ixloc, 1,
     &	sizeofreal, 
     &	LOCAL_4D_YSLICE, LOCAL_LARGE_4D_YSLICES, ierr)
      call MPI_TYPE_COMMIT(LOCAL_LARGE_4D_YSLICES, ierr)

      call MPI_TYPE_HVECTOR(ixloc-1, 1,
     &	NYworkers*sizeofreal, 
     &	GLOBAL_4D_YSLICE, GLOBAL_SMALL_4D_YSLICES, ierr)
      call MPI_TYPE_COMMIT(GLOBAL_SMALL_4D_YSLICES, ierr)
      call MPI_TYPE_HVECTOR(ixloc-1, 1,
     &	sizeofreal, 
     &	LOCAL_4D_YSLICE, LOCAL_SMALL_4D_YSLICES, ierr)
      call MPI_TYPE_COMMIT(LOCAL_SMALL_4D_YSLICES, ierr)
	end if

c Below version creates the array of yslices completely new
c  going through each dimension one after the other.
c  Therefore it is faster than above method.

c      if (.false.) then      

c Large global y-slice array is created going through each dimension
      call MPI_TYPE_VECTOR(ixloc, 1, NYworkers, MPI_REAL, 
     &	ONE_DGYl, ierr)
      call MPI_TYPE_HVECTOR(iy, 1, 
     &	ix*sizeofreal, ONE_DGYl, TWO_DGYl, ierr)
      call MPI_TYPE_HVECTOR(iz, 1, 
     &	ix*iy*sizeofreal, TWO_DGYl, THREE_DGYl, ierr)
      call MPI_TYPE_HVECTOR(is, 1,
     &	ix*iy*iz*sizeofreal, THREE_DGYl,
     &	GLOBAL_LARGE_4D_YSLICES, ierr)
      call MPI_TYPE_COMMIT(GLOBAL_LARGE_4D_YSLICES, ierr)
c Large local y-slice array is created going through each dimension
      call MPI_TYPE_VECTOR(ixloc, 1, 1, MPI_REAL, 
     &	ONE_DLYl, ierr)
      call MPI_TYPE_HVECTOR(iy, 1, 
     &	ixloc*sizeofreal, ONE_DLYl, TWO_DLYl, ierr)
      call MPI_TYPE_HVECTOR(iz, 1, 
     &	ixloc*iy*sizeofreal, TWO_DLYl, THREE_DLYl, ierr)
      call MPI_TYPE_HVECTOR(is, 1,
     &	ixloc*iy*iz*sizeofreal, THREE_DLYl,
     &	LOCAL_LARGE_4D_YSLICES, ierr)
      call MPI_TYPE_COMMIT(LOCAL_LARGE_4D_YSLICES, ierr)
     
c Small global y-slice array is created going through each dimension     
      call MPI_TYPE_VECTOR(ixloc-1, 1, NYworkers, MPI_REAL, 
     &	ONE_DGYs, ierr)
      call MPI_TYPE_HVECTOR(iy, 1, 
     &	ix*sizeofreal, ONE_DGYs, TWO_DGYs, ierr)
      call MPI_TYPE_HVECTOR(iz, 1, 
     &	ix*iy*sizeofreal, TWO_DGYs, THREE_DGYs, ierr)
      call MPI_TYPE_HVECTOR(is, 1,
     &	ix*iy*iz*sizeofreal, THREE_DGYs,
     &	GLOBAL_SMALL_4D_YSLICES, ierr)
      call MPI_TYPE_COMMIT(GLOBAL_SMALL_4D_YSLICES, ierr)
c Small local y-slice array is created going through each dimension
      call MPI_TYPE_VECTOR(ixloc-1, 1, 1, MPI_REAL, 
     &	ONE_DLYs, ierr)
      call MPI_TYPE_HVECTOR(iy, 1, 
     &	ixloc*sizeofreal, ONE_DLYs, TWO_DLYs, ierr)
      call MPI_TYPE_HVECTOR(iz, 1, 
     &	ixloc*iy*sizeofreal, TWO_DLYs, THREE_DLYs, ierr)
      call MPI_TYPE_HVECTOR(is, 1,
     &	ixloc*iy*iz*sizeofreal, THREE_DLYs,
     &	LOCAL_SMALL_4D_YSLICES, ierr)
      call MPI_TYPE_COMMIT(LOCAL_SMALL_4D_YSLICES, ierr)      
c      end if

      types(1) = GLOBAL_LARGE_4D_YSLICES
      call MPI_TYPE_STRUCT(2,blength,tdispls,types,
     &		GLOB_LARGE_4D_YSLICES, ierr)
      call MPI_TYPE_COMMIT(GLOB_LARGE_4D_YSLICES, ierr)      
      types(1) = GLOBAL_SMALL_4D_YSLICES
      call MPI_TYPE_STRUCT(2,blength,tdispls,types,
     &		GLOB_SMALL_4D_YSLICES, ierr)
      call MPI_TYPE_COMMIT(GLOB_SMALL_4D_YSLICES, ierr)      

c-------------------------------------------------------------------------------
c    Define 4D xslice set data types (large and small, differing in one element)
c-------------------------------------------------------------------------------

c This version creates the array of xslices using the given
c  xslice structures. It takes more time during communication.
      if (.false.) then
      call MPI_TYPE_HVECTOR(iyloc, 1, 
     &	NXworkers*sizeofreal*ix, 
     &	GLOBAL_4D_XSLICE, GLOBAL_LARGE_4D_XSLICES, ierr)
      call MPI_TYPE_COMMIT(GLOBAL_LARGE_4D_XSLICES, ierr)
      call MPI_TYPE_HVECTOR(iyloc, 1, 
     &	sizeofreal*ix, 
     &	LOCAL_4D_XSLICE, LOCAL_LARGE_4D_XSLICES, ierr)
      call MPI_TYPE_COMMIT(LOCAL_LARGE_4D_XSLICES, ierr)

      call MPI_TYPE_HVECTOR(iyloc-1, 1,
     &	NXworkers*sizeofreal*ix, 
     &	GLOBAL_4D_XSLICE, GLOBAL_SMALL_4D_XSLICES, ierr)
      call MPI_TYPE_COMMIT(GLOBAL_SMALL_4D_XSLICES, ierr)
      call MPI_TYPE_HVECTOR(iyloc-1, 1, 
     &	sizeofreal*ix, 
     &	LOCAL_4D_XSLICE, LOCAL_SMALL_4D_XSLICES, ierr)
      call MPI_TYPE_COMMIT(LOCAL_SMALL_4D_XSLICES, ierr)

      end if

c Below version creates the array of xslices completely new
c  going through each dimension one after the other.
c  Therefore it is faster than above method.

c      if (.false.) then

c Large global x-slice array is created going through each dimension
      call MPI_TYPE_VECTOR(iyloc, ix, NXworkers*ix, MPI_REAL, 
     &	TWO_DGXl, ierr)
      call MPI_TYPE_HVECTOR(iz, 1, 
     &	ix*iy*sizeofreal, TWO_DGXl, THREE_DGXl, ierr)
      call MPI_TYPE_HVECTOR(is, 1, 
     & ix*iy*iz*sizeofreal, THREE_DGXl,
     & 	GLOBAL_LARGE_4D_XSLICES, ierr)
      call MPI_TYPE_COMMIT(GLOBAL_LARGE_4D_XSLICES, ierr)
c Large local x-slice array is created going through each dimension      
      call MPI_TYPE_VECTOR(iyloc, ix, ix, MPI_REAL, 
     &	TWO_DLXl, ierr)
      call MPI_TYPE_HVECTOR(iz, 1, 
     &	ix*iyloc*sizeofreal, TWO_DLXl, THREE_DLXl, ierr)
      call MPI_TYPE_HVECTOR(is, 1, 
     & ix*iyloc*iz*sizeofreal, THREE_DLXl,
     & 	LOCAL_LARGE_4D_XSLICES, ierr)
      call MPI_TYPE_COMMIT(LOCAL_LARGE_4D_XSLICES, ierr)
      
c Small global x-slice array is created going through each dimension      
      call MPI_TYPE_VECTOR(iyloc-1, ix, NXworkers*ix, MPI_REAL, 
     &	TWO_DGXs, ierr)
      call MPI_TYPE_HVECTOR(iz, 1, 
     &	ix*iy*sizeofreal, TWO_DGXs, THREE_DGXs, ierr)
      call MPI_TYPE_HVECTOR(is, 1, 
     & ix*iy*iz*sizeofreal, THREE_DGXs,
     & 	GLOBAL_SMALL_4D_XSLICES, ierr)
      call MPI_TYPE_COMMIT(GLOBAL_SMALL_4D_XSLICES, ierr)
c Small local x-slice array is created going through each dimension      
      call MPI_TYPE_VECTOR(iyloc-1, ix, ix, MPI_REAL, 
     &	TWO_DLXs, ierr)
      call MPI_TYPE_HVECTOR(iz, 1, 
     &	ix*iyloc*sizeofreal, TWO_DLXs, THREE_DLXs, ierr)
      call MPI_TYPE_HVECTOR(is, 1, 
     & ix*iyloc*iz*sizeofreal, THREE_DLXs,
     & 	LOCAL_SMALL_4D_XSLICES, ierr)
      call MPI_TYPE_COMMIT(LOCAL_SMALL_4D_XSLICES, ierr)
c      end if

      types(1) = GLOBAL_LARGE_4D_XSLICES
      call MPI_TYPE_STRUCT(2,blength,tdispls,types,
     &		GLOB_LARGE_4D_XSLICES, ierr)
      call MPI_TYPE_COMMIT(GLOB_LARGE_4D_XSLICES, ierr)
      types(1) = GLOBAL_SMALL_4D_XSLICES
      call MPI_TYPE_STRUCT(2,blength,tdispls,
     &		types, GLOB_SMALL_4D_XSLICES, ierr)
      call MPI_TYPE_COMMIT(GLOB_SMALL_4D_XSLICES, ierr)
      
c-------------------------------------------------------------------------------
c    Define column data types for shuffeling
c-------------------------------------------------------------------------------

c This version creates the array of xslices-parts for shuffeling 
c  using the given x/y-column structures. It takes more time 
c  during communication.
      if (.false.) then
      call MPI_TYPE_HVECTOR(ixloc, 1, 
     &	NYworkers*sizeofreal, LOCAL_4D_XCOLUMN,
     &	LOCAL_3D_LY_XCOLUMN, ierr)
      call MPI_TYPE_HVECTOR(iyloc, 1, ix*sizeofreal,
     &	LOCAL_3D_LY_XCOLUMN, LOCAL_4D_LXLY_XCOLS, ierr)
      call MPI_TYPE_COMMIT(LOCAL_4D_LXLY_XCOLS, ierr)
      call MPI_TYPE_HVECTOR(iyloc-1, 1, ix*sizeofreal,
     &	LOCAL_3D_LY_XCOLUMN, LOCAL_4D_SXLY_XCOLS, ierr)
      call MPI_TYPE_COMMIT(LOCAL_4D_SXLY_XCOLS, ierr)

      call MPI_TYPE_HVECTOR(ixloc-1, 1, 
     &	NYworkers*sizeofreal, LOCAL_4D_XCOLUMN,
     &	LOCAL_3D_SY_XCOLUMN, ierr)
      call MPI_TYPE_HVECTOR(iyloc, 1, ix*sizeofreal,
     &	LOCAL_3D_SY_XCOLUMN, LOCAL_4D_LXSY_XCOLS, ierr)
      call MPI_TYPE_COMMIT(LOCAL_4D_LXSY_XCOLS, ierr)
      call MPI_TYPE_HVECTOR(iyloc-1, 1, ix*sizeofreal,
     &	LOCAL_3D_SY_XCOLUMN, LOCAL_4D_SXSY_XCOLS, ierr)
      call MPI_TYPE_COMMIT(LOCAL_4D_SXSY_XCOLS, ierr)
      

      call MPI_TYPE_HVECTOR(ixloc, 1, sizeofreal,
     &	LOCAL_4D_YCOLUMN, LOCAL_3D_LY_YCOLUMN, ierr)
      call MPI_TYPE_HVECTOR(iyloc, 1,
     &	NXworkers*ixloc*sizeofreal, LOCAL_3D_LY_YCOLUMN,
     &	LOCAL_4D_LXLY_YCOLS, ierr)
      call MPI_TYPE_COMMIT(LOCAL_4D_LXLY_YCOLS, ierr)
      call MPI_TYPE_HVECTOR(iyloc-1, 1,
     &	NXworkers*ixloc*sizeofreal, LOCAL_3D_LY_YCOLUMN,
     &	LOCAL_4D_SXLY_YCOLS, ierr)
      call MPI_TYPE_COMMIT(LOCAL_4D_SXLY_YCOLS, ierr)

      call MPI_TYPE_HVECTOR(ixloc-1, 1, sizeofreal,
     &	LOCAL_4D_YCOLUMN, LOCAL_3D_SY_YCOLUMN, ierr)
      call MPI_TYPE_HVECTOR(iyloc, 1,
     &	NXworkers*ixloc*sizeofreal, LOCAL_3D_SY_YCOLUMN,
     &	LOCAL_4D_LXSY_YCOLS, ierr)
      call MPI_TYPE_COMMIT(LOCAL_4D_LXSY_YCOLS, ierr)
      call MPI_TYPE_HVECTOR(iyloc-1, 1,
     &	NXworkers*ixloc*sizeofreal, LOCAL_3D_SY_YCOLUMN,
     &	LOCAL_4D_SXSY_YCOLS, ierr)
      call MPI_TYPE_COMMIT(LOCAL_4D_SXSY_YCOLS, ierr)
      end if

c Below version creates the array of x/y-slice-parts completely new
c  going through each dimension one after the other.
c  Therefore it is faster than above method.
c      if (.false.) then
      call MPI_TYPE_VECTOR(ixloc, 1, NYworkers, MPI_REAL, 
     &                      LOCAL_1D_LY_PART_XCOL, ierr)
      call MPI_TYPE_HVECTOR(iyloc, 1, ix*sizeofreal, 
     &	LOCAL_1D_LY_PART_XCOL, LOCAL_2D_LXLY_PART_XCOL, ierr)
      call MPI_TYPE_HVECTOR(iz, 1, ix*iyloc*sizeofreal, 
     &	LOCAL_2D_LXLY_PART_XCOL, LOCAL_3D_LXLY_PART_XCOL, ierr)
      call MPI_TYPE_HVECTOR(is, 1, ix*iyloc*iz*sizeofreal, 
     &	LOCAL_3D_LXLY_PART_XCOL, LOCAL_4D_LXLY_XCOLS, ierr)
      call MPI_TYPE_COMMIT(LOCAL_4D_LXLY_XCOLS, ierr)
      
      call MPI_TYPE_HVECTOR(iyloc-1, 1, ix*sizeofreal, 
     &	LOCAL_1D_LY_PART_XCOL, LOCAL_2D_SXLY_PART_XCOL, ierr)
      call MPI_TYPE_HVECTOR(iz, 1, ix*iyloc*sizeofreal, 
     &	LOCAL_2D_SXLY_PART_XCOL, LOCAL_3D_SXLY_PART_XCOL, ierr)
      call MPI_TYPE_HVECTOR(is, 1, ix*iyloc*iz*sizeofreal, 
     &	LOCAL_3D_SXLY_PART_XCOL, LOCAL_4D_SXLY_XCOLS, ierr)
      call MPI_TYPE_COMMIT(LOCAL_4D_SXLY_XCOLS, ierr)

      call MPI_TYPE_VECTOR(ixloc-1, 1, NYworkers, MPI_REAL, 
     &                      LOCAL_1D_SY_PART_XCOL, ierr)
      call MPI_TYPE_HVECTOR(iyloc, 1, ix*sizeofreal, 
     &	LOCAL_1D_SY_PART_XCOL, LOCAL_2D_LXSY_PART_XCOL, ierr)
      call MPI_TYPE_HVECTOR(iz, 1, ix*iyloc*sizeofreal, 
     &	LOCAL_2D_LXSY_PART_XCOL, LOCAL_3D_LXSY_PART_XCOL, ierr)
      call MPI_TYPE_HVECTOR(is, 1, ix*iyloc*iz*sizeofreal, 
     &	LOCAL_3D_LXSY_PART_XCOL, LOCAL_4D_LXSY_XCOLS, ierr)
      call MPI_TYPE_COMMIT(LOCAL_4D_LXSY_XCOLS, ierr)
      
      call MPI_TYPE_HVECTOR(iyloc-1, 1, ix*sizeofreal, 
     &	LOCAL_1D_SY_PART_XCOL, LOCAL_2D_SXSY_PART_XCOL, ierr)
      call MPI_TYPE_HVECTOR(iz, 1, ix*iyloc*sizeofreal, 
     &	LOCAL_2D_SXSY_PART_XCOL, LOCAL_3D_SXSY_PART_XCOL, ierr)
      call MPI_TYPE_HVECTOR(is, 1, ix*iyloc*iz*sizeofreal, 
     &	LOCAL_3D_SXSY_PART_XCOL, LOCAL_4D_SXSY_XCOLS, ierr)
      call MPI_TYPE_COMMIT(LOCAL_4D_SXSY_XCOLS, ierr)


      call MPI_TYPE_VECTOR(iyloc, ixloc, NXworkers*ixloc, MPI_REAL, 
     &			LOCAL_2D_LXLY_PART_YCOL, ierr)
      call MPI_TYPE_HVECTOR(iz, 1, ixloc*iy*sizeofreal, 
     &	LOCAL_2D_LXLY_PART_YCOL, LOCAL_3D_LXLY_PART_YCOL, ierr)
      call MPI_TYPE_HVECTOR(is, 1, ixloc*iy*iz*sizeofreal, 
     &	LOCAL_3D_LXLY_PART_YCOL, LOCAL_4D_LXLY_YCOLS, ierr)
      call MPI_TYPE_COMMIT(LOCAL_4D_LXLY_YCOLS, ierr)
      
      call MPI_TYPE_VECTOR(iyloc-1, ixloc, NXworkers*ixloc, MPI_REAL, 
     &			LOCAL_2D_SXLY_PART_YCOL, ierr)
      call MPI_TYPE_HVECTOR(iz, 1, ixloc*iy*sizeofreal, 
     &	LOCAL_2D_SXLY_PART_YCOL, LOCAL_3D_SXLY_PART_YCOL, ierr)
      call MPI_TYPE_HVECTOR(is, 1, ixloc*iy*iz*sizeofreal, 
     &	LOCAL_3D_SXLY_PART_YCOL, LOCAL_4D_SXLY_YCOLS, ierr)
      call MPI_TYPE_COMMIT(LOCAL_4D_SXLY_YCOLS, ierr)

      call MPI_TYPE_VECTOR(iyloc, ixloc-1, NXworkers*ixloc, MPI_REAL, 
     &			LOCAL_2D_LXSY_PART_YCOL, ierr)
      call MPI_TYPE_HVECTOR(iz, 1, ixloc*iy*sizeofreal, 
     &	LOCAL_2D_LXSY_PART_YCOL, LOCAL_3D_LXSY_PART_YCOL, ierr)
      call MPI_TYPE_HVECTOR(is, 1, ixloc*iy*iz*sizeofreal, 
     &	LOCAL_3D_LXSY_PART_YCOL, LOCAL_4D_LXSY_YCOLS, ierr)
      call MPI_TYPE_COMMIT(LOCAL_4D_LXSY_YCOLS, ierr)
      
      call MPI_TYPE_VECTOR(iyloc-1, ixloc-1, NXworkers*ixloc, MPI_REAL, 
     &			LOCAL_2D_SXSY_PART_YCOL, ierr)
      call MPI_TYPE_HVECTOR(iz, 1, ixloc*iy*sizeofreal, 
     &	LOCAL_2D_SXSY_PART_YCOL, LOCAL_3D_SXSY_PART_YCOL, ierr)
      call MPI_TYPE_HVECTOR(is, 1, ixloc*iy*iz*sizeofreal, 
     &	LOCAL_3D_SXSY_PART_YCOL, LOCAL_4D_SXSY_YCOLS, ierr)
      call MPI_TYPE_COMMIT(LOCAL_4D_SXSY_YCOLS, ierr)
      
c      end if

c-------------------------------------------------------------------------------
c    Define 2D data types
c-------------------------------------------------------------------------------
      call MPI_TYPE_VECTOR(iy, 1, ix, MPI_REAL,
     &				GLOBAL_2D_YSLICE, ierr)
      call MPI_TYPE_COMMIT(GLOBAL_2D_YSLICE, ierr)
      call MPI_TYPE_VECTOR(iy, 1, ixloc, MPI_REAL,
     &				LOCAL_2D_YSLICE, ierr)
      call MPI_TYPE_COMMIT(LOCAL_2D_YSLICE, ierr)

      types(1) = GLOBAL_2D_YSLICE
      call MPI_TYPE_STRUCT(2,blength,tdispls,
     &		types, GLOB_2D_YSLICE, ierr)
      call MPI_TYPE_COMMIT(GLOB_2D_YSLICE, ierr)

c This version creates the array of x/y-slices for 2D distribution
c  using the given x/y-slice structures. 
      call MPI_TYPE_HVECTOR(iyloc, ix, 
     &		NXworkers*sizeofreal*ix, MPI_REAL,
     &		GLOBAL_LARGE_2D_XSLICES, ierr)
      call MPI_TYPE_COMMIT(GLOBAL_LARGE_2D_XSLICES, ierr)
      call MPI_TYPE_HVECTOR(iyloc-1, ix, 
     &		NXworkers*sizeofreal*ix, MPI_REAL,
     &		GLOBAL_SMALL_2D_XSLICES, ierr)
      call MPI_TYPE_COMMIT(GLOBAL_SMALL_2D_XSLICES, ierr)

      call MPI_TYPE_HVECTOR(iyloc, ix, 
     &		ix*sizeofreal, MPI_REAL,
     &		LOCAL_LARGE_2D_XSLICES, ierr)
      call MPI_TYPE_COMMIT(LOCAL_LARGE_2D_XSLICES, ierr)
      call MPI_TYPE_HVECTOR(iyloc-1, ix, 
     &		ix*sizeofreal, MPI_REAL,
     &		LOCAL_SMALL_2D_XSLICES, ierr)
      call MPI_TYPE_COMMIT(LOCAL_SMALL_2D_XSLICES, ierr)

      call MPI_TYPE_HVECTOR(ixloc, 1, 
     &		NYworkers*sizeofreal, GLOBAL_2D_YSLICE,
     &		GLOBAL_LARGE_2D_YSLICES, ierr)
      call MPI_TYPE_COMMIT(GLOBAL_LARGE_2D_YSLICES, ierr)
      call MPI_TYPE_HVECTOR(ixloc-1, 1, 
     &		NYworkers*sizeofreal, GLOBAL_2D_YSLICE,
     &		GLOBAL_SMALL_2D_YSLICES, ierr)
      call MPI_TYPE_COMMIT(GLOBAL_SMALL_2D_YSLICES, ierr)

      call MPI_TYPE_HVECTOR(ixloc, 1, 
     &		sizeofreal, LOCAL_2D_YSLICE,
     &		LOCAL_LARGE_2D_YSLICES, ierr)
      call MPI_TYPE_COMMIT(LOCAL_LARGE_2D_YSLICES, ierr)
      call MPI_TYPE_HVECTOR(ixloc-1, 1, 
     &		sizeofreal, LOCAL_2D_YSLICE,
     &		LOCAL_SMALL_2D_YSLICES, ierr)
      call MPI_TYPE_COMMIT(LOCAL_SMALL_2D_YSLICES, ierr)

      types(1) = GLOBAL_LARGE_2D_XSLICES
      call MPI_TYPE_STRUCT(2,blength,tdispls,
     &		types, GLOB_LARGE_2D_XSLICES, ierr)
      call MPI_TYPE_COMMIT(GLOB_LARGE_2D_XSLICES, ierr)
      types(1) = GLOBAL_SMALL_2D_XSLICES
      call MPI_TYPE_STRUCT(2,blength,tdispls,
     &		types, GLOB_SMALL_2D_XSLICES, ierr)
      call MPI_TYPE_COMMIT(GLOB_SMALL_2D_XSLICES, ierr)
      types(1) = GLOBAL_LARGE_2D_YSLICES
      call MPI_TYPE_STRUCT(2,blength,tdispls,
     &		types, GLOB_LARGE_2D_YSLICES, ierr)
      call MPI_TYPE_COMMIT(GLOB_LARGE_2D_YSLICES, ierr)      
      types(1) = GLOBAL_SMALL_2D_YSLICES
      call MPI_TYPE_STRUCT(2,blength,tdispls,
     &		types, GLOB_SMALL_2D_YSLICES, ierr)
      call MPI_TYPE_COMMIT(GLOB_SMALL_2D_YSLICES, ierr)      
	
c-------------------------------------------------------------------------------
c    Define 2DN data types
c-------------------------------------------------------------------------------
      call MPI_TYPE_VECTOR(is, ix, ix*iy, MPI_REAL,
     &				GLOBAL_2DN_XSLICE, ierr)
      call MPI_TYPE_COMMIT(GLOBAL_2DN_XSLICE, ierr)
      call MPI_TYPE_VECTOR(is, ix, ix*iyloc, MPI_REAL,
     &				LOCAL_2DN_XSLICE, ierr)
      call MPI_TYPE_COMMIT(LOCAL_2DN_XSLICE, ierr)

      call MPI_TYPE_VECTOR(iy*is, 1, ix, MPI_REAL,
     &				GLOBAL_2DN_YSLICE, ierr)
      call MPI_TYPE_COMMIT(GLOBAL_2DN_YSLICE, ierr)
      call MPI_TYPE_VECTOR(iy*is, 1, ixloc, MPI_REAL,
     &				LOCAL_2DN_YSLICE, ierr)
      call MPI_TYPE_COMMIT(LOCAL_2DN_YSLICE, ierr)

      types(1) = GLOBAL_2DN_XSLICE
      call MPI_TYPE_STRUCT(2,blength,tdispls,
     &		types, GLOB_2DN_XSLICE, ierr)
      call MPI_TYPE_COMMIT(GLOB_2DN_XSLICE, ierr)
      types(1) = GLOBAL_2DN_YSLICE
      call MPI_TYPE_STRUCT(2,blength,tdispls,
     &		types, GLOB_2DN_YSLICE, ierr)
      call MPI_TYPE_COMMIT(GLOB_2DN_YSLICE, ierr)

c This version creates the array of x/y-slices for 2DN distribution
c  using the given x/y-slice structures. 
      call MPI_TYPE_HVECTOR(iyloc, 1, 
     &		NXworkers*sizeofreal*ix, GLOBAL_2DN_XSLICE,
     &		GLOBAL_LARGE_2DN_XSLICES, ierr)
      call MPI_TYPE_COMMIT(GLOBAL_LARGE_2DN_XSLICES, ierr)
      call MPI_TYPE_HVECTOR(iyloc-1, 1, 
     &		NXworkers*sizeofreal*ix, GLOBAL_2DN_XSLICE,
     &		GLOBAL_SMALL_2DN_XSLICES, ierr)
      call MPI_TYPE_COMMIT(GLOBAL_SMALL_2DN_XSLICES, ierr)

      call MPI_TYPE_HVECTOR(iyloc, 1, 
     &		sizeofreal*ix, LOCAL_2DN_XSLICE,
     &		LOCAL_LARGE_2DN_XSLICES, ierr)
      call MPI_TYPE_COMMIT(LOCAL_LARGE_2DN_XSLICES, ierr)
      call MPI_TYPE_HVECTOR(iyloc-1, 1, 
     &		sizeofreal*ix, LOCAL_2DN_XSLICE,
     &		LOCAL_SMALL_2DN_XSLICES, ierr)
      call MPI_TYPE_COMMIT(LOCAL_SMALL_2DN_XSLICES, ierr)

      call MPI_TYPE_HVECTOR(ixloc, 1, 
     &		NYworkers*sizeofreal, GLOBAL_2DN_YSLICE,
     &		GLOBAL_LARGE_2DN_YSLICES, ierr)
      call MPI_TYPE_COMMIT(GLOBAL_LARGE_2DN_YSLICES, ierr)
      call MPI_TYPE_HVECTOR(ixloc-1, 1, 
     &		NYworkers*sizeofreal, GLOBAL_2DN_YSLICE,
     &		GLOBAL_SMALL_2DN_YSLICES, ierr)
      call MPI_TYPE_COMMIT(GLOBAL_SMALL_2DN_YSLICES, ierr)

      call MPI_TYPE_HVECTOR(ixloc, 1, 
     &		sizeofreal, LOCAL_2DN_YSLICE,
     &		LOCAL_LARGE_2DN_YSLICES, ierr)
      call MPI_TYPE_COMMIT(LOCAL_LARGE_2DN_YSLICES, ierr)
      call MPI_TYPE_HVECTOR(ixloc-1, 1, 
     &		sizeofreal, LOCAL_2DN_YSLICE,
     &		LOCAL_SMALL_2DN_YSLICES, ierr)
      call MPI_TYPE_COMMIT(LOCAL_SMALL_2DN_YSLICES, ierr)

      types(1) = GLOBAL_LARGE_2DN_XSLICES
      call MPI_TYPE_STRUCT(2,blength,tdispls,
     &		types, GLOB_LARGE_2DN_XSLICES, ierr)
      call MPI_TYPE_COMMIT(GLOB_LARGE_2DN_XSLICES, ierr)
      types(1) = GLOBAL_SMALL_2DN_XSLICES
      call MPI_TYPE_STRUCT(2,blength,tdispls,
     &		types, GLOB_SMALL_2DN_XSLICES, ierr)
      call MPI_TYPE_COMMIT(GLOB_SMALL_2DN_XSLICES, ierr)
      types(1) = GLOBAL_LARGE_2DN_YSLICES
      call MPI_TYPE_STRUCT(2,blength,tdispls,
     &		types, GLOB_LARGE_2DN_YSLICES, ierr)
      call MPI_TYPE_COMMIT(GLOB_LARGE_2DN_YSLICES, ierr)      
      types(1) = GLOBAL_SMALL_2DN_YSLICES
      call MPI_TYPE_STRUCT(2,blength,tdispls,
     &		types, GLOB_SMALL_2DN_YSLICES, ierr)
      call MPI_TYPE_COMMIT(GLOB_SMALL_2DN_YSLICES, ierr)      

c-------------------------------------------------------------------------------
c    Define 3D data types
c-------------------------------------------------------------------------------
      call MPI_TYPE_VECTOR(iz, ix, ix*iy, MPI_REAL,
     &				GLOBAL_3D_XSLICE, ierr)
      call MPI_TYPE_COMMIT(GLOBAL_3D_XSLICE, ierr)
      call MPI_TYPE_VECTOR(iz, ix, ix*iyloc, MPI_REAL,
     &				LOCAL_3D_XSLICE, ierr)
      call MPI_TYPE_COMMIT(LOCAL_3D_XSLICE, ierr)

      call MPI_TYPE_VECTOR(iy*iz, 1, ix, MPI_REAL,
     &				GLOBAL_3D_YSLICE, ierr)
      call MPI_TYPE_COMMIT(GLOBAL_3D_YSLICE, ierr)
      call MPI_TYPE_VECTOR(iy*iz, 1, ixloc, MPI_REAL,
     &				LOCAL_3D_YSLICE, ierr)
      call MPI_TYPE_COMMIT(LOCAL_3D_YSLICE, ierr)

      types(1) = GLOBAL_3D_XSLICE
      call MPI_TYPE_STRUCT(2,blength,tdispls,
     &		types, GLOB_3D_XSLICE, ierr)
      call MPI_TYPE_COMMIT(GLOB_3D_XSLICE, ierr)
      types(1) = GLOBAL_3D_YSLICE
      call MPI_TYPE_STRUCT(2,blength,tdispls,
     &		types, GLOB_3D_YSLICE, ierr)
      call MPI_TYPE_COMMIT(GLOB_3D_YSLICE, ierr)

c This version creates the array of x/y-slices for 3D distribution
c  using the given x/y-slice structures.
      call MPI_TYPE_HVECTOR(iyloc, 1, 
     &		NXworkers*sizeofreal*ix, GLOBAL_3D_XSLICE,
     &		GLOBAL_LARGE_3D_XSLICES, ierr)
      call MPI_TYPE_COMMIT(GLOBAL_LARGE_3D_XSLICES, ierr)
      call MPI_TYPE_HVECTOR(iyloc-1, 1, 
     &		NXworkers*sizeofreal*ix, GLOBAL_3D_XSLICE,
     &		GLOBAL_SMALL_3D_XSLICES, ierr)
      call MPI_TYPE_COMMIT(GLOBAL_SMALL_3D_XSLICES, ierr)

      call MPI_TYPE_HVECTOR(iyloc, 1, 
     &		sizeofreal*ix, LOCAL_3D_XSLICE,
     &		LOCAL_LARGE_3D_XSLICES, ierr)
      call MPI_TYPE_COMMIT(LOCAL_LARGE_3D_XSLICES, ierr)
      call MPI_TYPE_HVECTOR(iyloc-1, 1, 
     &		sizeofreal*ix, LOCAL_3D_XSLICE,
     &		LOCAL_SMALL_3D_XSLICES, ierr)
      call MPI_TYPE_COMMIT(LOCAL_SMALL_3D_XSLICES, ierr)

      call MPI_TYPE_HVECTOR(ixloc, 1, 
     &		NYworkers*sizeofreal, GLOBAL_3D_YSLICE,
     &		GLOBAL_LARGE_3D_YSLICES, ierr)
      call MPI_TYPE_COMMIT(GLOBAL_LARGE_3D_YSLICES, ierr)
      call MPI_TYPE_HVECTOR(ixloc-1, 1, 
     &		NYworkers*sizeofreal, GLOBAL_3D_YSLICE,
     &		GLOBAL_SMALL_3D_YSLICES, ierr)
      call MPI_TYPE_COMMIT(GLOBAL_SMALL_3D_YSLICES, ierr)

      call MPI_TYPE_HVECTOR(ixloc, 1, 
     &		sizeofreal, LOCAL_3D_YSLICE,
     &		LOCAL_LARGE_3D_YSLICES, ierr)
      call MPI_TYPE_COMMIT(LOCAL_LARGE_3D_YSLICES, ierr)
      call MPI_TYPE_HVECTOR(ixloc-1, 1, 
     &		sizeofreal, LOCAL_3D_YSLICE,
     &		LOCAL_SMALL_3D_YSLICES, ierr)
      call MPI_TYPE_COMMIT(LOCAL_SMALL_3D_YSLICES, ierr)

      types(1) = GLOBAL_LARGE_3D_XSLICES
      call MPI_TYPE_STRUCT(2,blength,tdispls,
     &		types, GLOB_LARGE_3D_XSLICES, ierr)
      call MPI_TYPE_COMMIT(GLOB_LARGE_3D_XSLICES, ierr)
      types(1) = GLOBAL_SMALL_3D_XSLICES
      call MPI_TYPE_STRUCT(2,blength,tdispls,
     &		types, GLOB_SMALL_3D_XSLICES, ierr)
      call MPI_TYPE_COMMIT(GLOB_SMALL_3D_XSLICES, ierr)
      types(1) = GLOBAL_LARGE_3D_YSLICES
      call MPI_TYPE_STRUCT(2,blength,tdispls,
     &		types, GLOB_LARGE_3D_YSLICES, ierr)
      call MPI_TYPE_COMMIT(GLOB_LARGE_3D_YSLICES, ierr)      
      types(1) = GLOBAL_SMALL_3D_YSLICES
      call MPI_TYPE_STRUCT(2,blength,tdispls,
     &		types, GLOB_SMALL_3D_YSLICES, ierr)
      call MPI_TYPE_COMMIT(GLOB_SMALL_3D_YSLICES, ierr)      

c-------------------------------------------------------------------------------
c    Define BD data types
c-------------------------------------------------------------------------------
      call MPI_TYPE_VECTOR(iz*2*is, 1, iy, MPI_REAL,
     &				GLOBAL_BD_XSLICE, ierr)
      call MPI_TYPE_COMMIT(GLOBAL_BD_XSLICE, ierr)
      call MPI_TYPE_VECTOR(iz*2*is, 1, iyloc, MPI_REAL,
     &				LOCAL_BD_XSLICE, ierr)
      call MPI_TYPE_COMMIT(LOCAL_BD_XSLICE, ierr)

      call MPI_TYPE_VECTOR(iz*2*is, 1, ix, MPI_REAL,
     &				GLOBAL_BD_YSLICE, ierr)
      call MPI_TYPE_COMMIT(GLOBAL_BD_YSLICE, ierr)
      call MPI_TYPE_VECTOR(iz*2*is, 1, ixloc, MPI_REAL,
     &				LOCAL_BD_YSLICE, ierr)
      call MPI_TYPE_COMMIT(LOCAL_BD_YSLICE, ierr)

      types(1) = GLOBAL_BD_XSLICE
      call MPI_TYPE_STRUCT(2,blength,tdispls,
     &		types, GLOB_BD_XSLICE, ierr)
      call MPI_TYPE_COMMIT(GLOB_BD_XSLICE, ierr)
      types(1) = GLOBAL_BD_YSLICE
      call MPI_TYPE_STRUCT(2,blength,tdispls,
     &		types, GLOB_BD_YSLICE, ierr)
      call MPI_TYPE_COMMIT(GLOB_BD_YSLICE, ierr)

c This version creates the array of x/y-slices for BD distribution
c  using the given x/y-slice structures.
      call MPI_TYPE_HVECTOR(iyloc, 1, 
     &		NXworkers*sizeofreal, GLOBAL_BD_XSLICE,
     &		GLOBAL_LARGE_BD_XSLICES, ierr)
      call MPI_TYPE_COMMIT(GLOBAL_LARGE_BD_XSLICES, ierr)
      call MPI_TYPE_HVECTOR(iyloc-1, 1, 
     &		NXworkers*sizeofreal, GLOBAL_BD_XSLICE,
     &		GLOBAL_SMALL_BD_XSLICES, ierr)
      call MPI_TYPE_COMMIT(GLOBAL_SMALL_BD_XSLICES, ierr)

      call MPI_TYPE_HVECTOR(iyloc, 1, 
     &		sizeofreal, LOCAL_BD_XSLICE,
     &		LOCAL_LARGE_BD_XSLICES, ierr)
      call MPI_TYPE_COMMIT(LOCAL_LARGE_BD_XSLICES, ierr)
      call MPI_TYPE_HVECTOR(iyloc-1, 1, 
     &		sizeofreal, LOCAL_BD_XSLICE,
     &		LOCAL_SMALL_BD_XSLICES, ierr)
      call MPI_TYPE_COMMIT(LOCAL_SMALL_BD_XSLICES, ierr)

      call MPI_TYPE_HVECTOR(ixloc, 1, 
     &		NYworkers*sizeofreal, GLOBAL_BD_YSLICE,
     &		GLOBAL_LARGE_BD_YSLICES, ierr)
      call MPI_TYPE_COMMIT(GLOBAL_LARGE_BD_YSLICES, ierr)
      call MPI_TYPE_HVECTOR(ixloc-1, 1, 
     &		NYworkers*sizeofreal, GLOBAL_BD_YSLICE,
     &		GLOBAL_SMALL_BD_YSLICES, ierr)
      call MPI_TYPE_COMMIT(GLOBAL_SMALL_BD_YSLICES, ierr)

      call MPI_TYPE_HVECTOR(ixloc, 1, 
     &		sizeofreal, LOCAL_BD_YSLICE,
     &		LOCAL_LARGE_BD_YSLICES, ierr)
      call MPI_TYPE_COMMIT(LOCAL_LARGE_BD_YSLICES, ierr)
      call MPI_TYPE_HVECTOR(ixloc-1, 1, 
     &		sizeofreal, LOCAL_BD_YSLICE,
     &		LOCAL_SMALL_BD_YSLICES, ierr)
      call MPI_TYPE_COMMIT(LOCAL_SMALL_BD_YSLICES, ierr)

      types(1) = GLOBAL_LARGE_BD_XSLICES
      call MPI_TYPE_STRUCT(2,blength,tdispls,
     &		types, GLOB_LARGE_BD_XSLICES, ierr)
      call MPI_TYPE_COMMIT(GLOB_LARGE_BD_XSLICES, ierr)
      types(1) = GLOBAL_SMALL_BD_XSLICES
      call MPI_TYPE_STRUCT(2,blength,tdispls,
     &		types, GLOB_SMALL_BD_XSLICES, ierr)
      call MPI_TYPE_COMMIT(GLOB_SMALL_BD_XSLICES, ierr)
      types(1) = GLOBAL_LARGE_BD_YSLICES
      call MPI_TYPE_STRUCT(2,blength,tdispls,
     &		types, GLOB_LARGE_BD_YSLICES, ierr)
      call MPI_TYPE_COMMIT(GLOB_LARGE_BD_YSLICES, ierr)      
      types(1) = GLOBAL_SMALL_BD_YSLICES
      call MPI_TYPE_STRUCT(2,blength,tdispls,
     &		types, GLOB_SMALL_BD_YSLICES, ierr)
      call MPI_TYPE_COMMIT(GLOB_SMALL_BD_YSLICES, ierr)      

      end subroutine CreateCommDataTypes 

      end module XYCommDataTypes
c---------------------------------------------------------------------------------------
c
c  The following function defines the mapping of the global
c  domain to local workers
c
c      owner_of_hslice(1..iz) = the process which gets that h-slice
c      owner_of_vcol(1..ix*iy) = the process which gets that v-column
c
c      no_of_hslices(p),  1<=p<=NHworkers = the number of h-slices	 
c                         owned by process p
c      no_of_vcols(p),  1<=p<=NVworkers = the number of v-columns	 
c                         owned by process p
c
c      local_hslice_id(1..iz) = local index of global hslice 1..iz
c                           (index seen by the process which owns it)
c      local_vcol_id(1..ix*iy) = local index of global vcolumn 1..ix*iy
c                           (index seen by the process which owns it)
c
c      global_hslice_id(p,i) 1<=p<=NHworkers, 1<=i<=izloc
c               hslice of local index i on process p is the global slice ...
c      global_yslice_id(1..NVworkers,1..icloc)
c               vcolumn of local index i on process p is the global slice ...
c
c      owned_hslices(p,i) = 1..iz, 1<=p<=NHworkers, 1<=i<=izloc 
c                which global h-slices i are owned by slave p   
c      owned_vcols(p,i)= 1..ix*iy, 1<=p<=Nvworkers, 1<=i<=icloc 
c                which v-columns i are owned by slave p   
c
c      planar_vcol_id(1..ix*iy,2) = x and y index of global vcolumn id
c
c      linear_vcol_id(ix,iy) = global vcolumn id for array index
c
c---------------------------------------------------------------------------------------

      module HVParallelDataMap
c       
c  No. of processes
        integer :: Nprocs       
c  No. of slave processes
        integer :: NHworkers, NVworkers
c  Current process id
        integer :: MyId
c Status as Master or Worker	
        logical :: Master, HWorker, VWorker
c       
        integer, pointer, dimension(:)   :: no_of_hslices
	integer, pointer, dimension(:)   :: no_of_vcols
        integer, pointer, dimension(:,:) :: global_hslice_id
        integer, pointer, dimension(:,:) :: global_vcol_id
        integer, pointer, dimension(:)   :: local_hslice_id
        integer, pointer, dimension(:)   :: local_vcol_id
        integer, pointer, dimension(:)   :: owner_of_hslice
	integer, pointer, dimension(:)   :: owner_of_vcol
	integer, pointer, dimension(:,:) :: owned_hslices
	integer, pointer, dimension(:,:) :: owned_vcols
        integer, pointer, dimension(:,:) :: planar_vcol_id, 
     &              linear_vcol_id


      contains 
c       
      subroutine CreateMap(ix,iy,iz,izloc,icloc) 
c    
      integer :: ix, iy, iz, izloc, icloc
      integer :: i, id, p
      integer :: modz, modc

      Master   = MyId.eq.0
      HWorker  = (MyId.ne.0) .and. (MyId.le.iz)
      VWorker  = (MyId.ne.0) .and. (MyId.le.(ix*iy))
      if (Nprocs<=1) return

      if (Nprocs-1 .le. iz) then
        NHworkers = Nprocs-1
      else
        NHworkers = iz
      end if
      if (Nprocs-1 .le. ix*iy) then
        NVworkers = Nprocs-1
      else
        NVworkers = ix*iy
      end if
	 

c The dimensions of the local data sets
      modz = mod(iz, NHworkers)
      if ( modz .eq. 0 ) then
         izloc = iz/NHworkers
      else
         izloc = (iz-modz)/NHworkers + 1
      endif	
c
      modc = mod(ix*iy, NVworkers)
      if ( modc .eq. 0 ) then
         icloc = ix*iy/NVworkers
      else
         icloc = (ix*iy-modc)/NVworkers + 1
      endif
c-----------------------------------------------------------------
c Allocate stuff
      allocate( no_of_hslices(NHworkers),STAT=ierr)
      allocate( no_of_vcols(NVworkers),STAT=ierr)      
      allocate( global_hslice_id(NHworkers,izloc),STAT=ierr) 	 
      allocate( global_vcol_id(NVworkers,icloc),STAT=ierr) 	 
      allocate( local_hslice_id(iz),STAT=ierr)
      allocate( local_vcol_id(ix*iy),STAT=ierr)
      allocate( owner_of_hslice(iz),STAT=ierr)
      allocate( owner_of_vcol(ix*iy),STAT=ierr)
      allocate( owned_hslices(NHworkers,izloc),STAT=ierr )
      allocate( owned_vcols(NVworkers,icloc),  STAT=ierr )
      allocate( planar_vcol_id(ix*iy,2), STAT=ierr )
      allocate( linear_vcol_id(ix,iy),   STAT=ierr )
c-------------------------------------------------      
      
      do i=1,iz
         owner_of_hslice(i)=mod(i-1,NHworkers)+1
      enddo
c The structure of the vslices is set here:
c  first through x, then through y
      k = 0
      do j=1,iy
        do i=1,ix
	  k = k+1
          owner_of_vcol(k)=mod(k-1,NVworkers)+1
          linear_vcol_id(i,j)=k
          planar_vcol_id(k,1:2)=(/i,j/)
        enddo
      enddo
c      do i=1,ix*iy
c        owner_of_vcol(i)=mod(i-1,NVworkers)+1
c        linear_vcol_id(mod(i-1,ix)+1,(i-1)/ix+1)=i
c        planar_vcol_id(i,1:2)=(/mod(i-1,ix)+1,(i-1)/ix+1/)
c      enddo

      do p=1,NVworkers
      id = 0
      do i=1,ix*iy
        if ( owner_of_vcol(i) .eq. p) then
          id = id+1
          global_vcol_id(p,id) = i
          local_vcol_id(i) = id
        endif
      enddo
      if (id .gt. icloc) then
        print*,' Process ',p,' is given ',id,
     &          ' vcols, but icloc = ',icloc
        stop
      else
        no_of_vcols(p) = id
      end if
      enddo
       
      do p=1,NHworkers
      id = 0
      do i=1,iz
         if ( owner_of_hslice(i) .eq. p) then
	   id = id+1
	    global_hslice_id(p,id) = i
	    local_hslice_id(i)     = id
	 endif
      enddo
      if (id .gt. izloc) then
        print*,' Process ',p,' is given ',id,
     &          ' hslices, but izloc = ',izloc
        stop
      else
         no_of_hslices(p) = id
      endif	 
      enddo

      do p=1,NHworkers
c
        do i=1,no_of_hslices(p)
          owned_hslices(p,i) = p+(i-1)*NHworkers
        enddo
        do i=no_of_hslices(p)+1,izloc
          owned_hslices(p,i) = 0
        enddo
      end do
      
      do p=1,NVworkers
        do i=1,no_of_vcols(p)
          owned_vcols(p,i) =  p+(i-1)*NVworkers
        enddo
        do i=no_of_vcols(p)+1,icloc
          owned_vcols(p,i) = 0
        enddo
c
      enddo

      end subroutine CreateMap
c      
      end module HVParallelDataMap

c-------------------------------------------------------------------------------
C Handles for special communication data types
c-------------------------------------------------------------------------------
      module HVCommDataTypes
        integer :: GROUP_WORKERS, MPI_COMM_WORKERS
c	for standard versions of H and V distribution and gathering
c	for necessary types of H-slice and V-column: 
c		2D, 2DN, 3D, BD, 4D, local, global
        integer :: GLOBAL_4D_VCOL, LOCAL_4D_VCOL
        integer :: GLOBAL_4D_HSLICE, LOCAL_4D_HSLICE
	integer :: GLOBAL_2DN_VCOL, LOCAL_2DN_VCOL
	integer :: GLOBAL_3D_VCOL, LOCAL_3D_VCOL
	integer :: GLOBAL_BD_XHSLICE, LOCAL_BD_XHSLICE
	integer :: GLOBAL_BD_YHSLICE, LOCAL_BD_YHSLICE
c	for single versions of H and V distribution and gathering:
c	- large means the local h-slice/v-column array is completely filled
c	- small means there is one h-slice/v-column less in the array than
c	    in others
	integer :: GLOBAL_LARGE_4D_HSLICES, LOCAL_LARGE_4D_HSLICES
	integer :: GLOBAL_LARGE_4D_VCOLS, LOCAL_LARGE_4D_VCOLS
	integer :: GLOBAL_SMALL_4D_HSLICES, LOCAL_SMALL_4D_HSLICES
	integer :: GLOBAL_SMALL_4D_VCOLS, LOCAL_SMALL_4D_VCOLS
	integer :: GLOBAL_LARGE_2D_VCOLS, LOCAL_LARGE_2D_VCOLS
	integer :: GLOBAL_SMALL_2D_VCOLS, LOCAL_SMALL_2D_VCOLS
	integer :: GLOBAL_LARGE_2DN_VCOLS, LOCAL_LARGE_2DN_VCOLS
	integer :: GLOBAL_SMALL_2DN_VCOLS, LOCAL_SMALL_2DN_VCOLS
	integer :: GLOBAL_LARGE_3D_HSLICES, LOCAL_LARGE_3D_HSLICES
	integer :: GLOBAL_SMALL_3D_HSLICES, LOCAL_SMALL_3D_HSLICES
	integer :: GLOBAL_LARGE_3D_VCOLS, LOCAL_LARGE_3D_VCOLS
	integer :: GLOBAL_SMALL_3D_VCOLS, LOCAL_SMALL_3D_VCOLS
	integer :: GLOBAL_LARGE_BD_XHSLICES, LOCAL_LARGE_BD_XHSLICES
	integer :: GLOBAL_SMALL_BD_XHSLICES, LOCAL_SMALL_BD_XHSLICES
	integer :: GLOBAL_LARGE_BD_YHSLICES, LOCAL_LARGE_BD_YHSLICES
	integer :: GLOBAL_SMALL_BD_YHSLICES, LOCAL_SMALL_BD_YHSLICES
c	for H2V and V2H shuffeling
c	- LH: large in the h-slice array part (full)
c	- SH: small in the h-slice array part (one hslice less)
c	- LV: large in the v-column array part (full)
c	- SV: small in the v-column array part (one vcolumn less)
c	- all possibilities for sending and receiving are combined
	integer :: LOCAL_4D_LHLV_HCOLS, LOCAL_4D_LHSV_HCOLS
	integer :: LOCAL_4D_SHLV_HCOLS, LOCAL_4D_SHSV_HCOLS
	integer :: LOCAL_4D_LHLV_VCOLS, LOCAL_4D_LHSV_VCOLS
	integer :: LOCAL_4D_SHLV_VCOLS, LOCAL_4D_SHSV_VCOLS
c	for collective communication, types used for standard and
c	single communication with reset extents
	integer :: GLOB_LARGE_4D_HSLICES, GLOB_SMALL_4D_HSLICES
	integer :: GLOB_4D_HSLICE, GLOB_4D_VCOL
	integer :: GLOB_LARGE_4D_VCOLS, GLOB_SMALL_4D_VCOLS
	integer :: GLOB_LARGE_2D_VCOLS, GLOB_SMALL_2D_VCOLS
	integer :: GLOB_2DN_VCOL
	integer :: GLOB_LARGE_2DN_VCOLS, GLOB_SMALL_2DN_VCOLS
	integer :: GLOB_LARGE_3D_HSLICES, GLOB_SMALL_3D_HSLICES
	integer :: GLOB_3D_VCOL
	integer :: GLOB_LARGE_3D_VCOLS, GLOB_SMALL_3D_VCOLS
	integer :: GLOB_LARGE_BD_XHSLICES, GLOB_SMALL_BD_XHSLICES
	integer :: GLOB_BD_XHSLICE, GLOB_BD_YHSLICE
	integer :: GLOB_LARGE_BD_YHSLICES, GLOB_SMALL_BD_YHSLICES
 
      contains

c 	Creates special communication data types
      subroutine CreateCommDataTypes( ix,iy,iz,izloc,icloc,is)
      use HVParallelDataMap
      include 'mpif.h'
      integer :: THREE_DGHl, THREE_DLHl
      integer :: THREE_DGHs, THREE_DLHs
      integer :: THREE_DGVl, THREE_DLVl
      integer :: THREE_DGVs, THREE_DLVs
      integer :: TWO_DGVl, TWO_DLVl
      integer :: TWO_DGVs, TWO_DLVs
      
      integer :: LOCAL_2D_LV_PART_HCOL, LOCAL_2D_SV_PART_HCOL
      integer :: LOCAL_3D_LHLV_PART_HCOL, LOCAL_3D_LHSV_PART_HCOL
      integer :: LOCAL_3D_SHLV_PART_HCOL, LOCAL_3D_SHSV_PART_HCOL
      integer :: LOCAL_2D_LV_PART_VCOL, LOCAL_2D_SV_PART_VCOL
      integer :: LOCAL_3D_LHLV_PART_VCOL, LOCAL_3D_LHSV_PART_VCOL
      integer :: LOCAL_3D_SHLV_PART_VCOL, LOCAL_3D_SHSV_PART_VCOL
      integer :: GLOBAL_3D_HSLICE
      
      integer :: ix,iy,iz,izloc,icloc,is
      integer :: i, ierr, sizeofreal, lb
      integer :: idle, Nw, wid, GROUP_WORLD
      integer, pointer, dimension (:) :: ranks
      integer :: blength(2), types(2), tdispls(2)

      call MPI_TYPE_EXTENT(MPI_REAL, sizeofreal,ierr)
c      print*, 'lb = ',lb,' sizeofreal =', sizeofreal      
c Set the struct parameters in order to do the resizing of data types
      types(2) = MPI_UB
      tdispls(1) = 0
      tdispls(2) = sizeofreal
      blength(1:2) = 1

c-------------------------------------------------------------------------------
c    Create workers group and workers communicator
c-------------------------------------------------------------------------------
      call MPI_COMM_GROUP(MPI_COMM_WORLD, GROUP_WORLD, ierr)
      idle = Nprocs - max(NHworkers,NVworkers)
      allocate(ranks(idle), STAT=ierr)
      ranks(1) = 0
      if (idle .gt. 1) then
        do i=2, idle
	  ranks(i) = i - 1 + max(NHworkers,NVworkers)
	  if (Master) print*, 'Idle processor:',ranks(i)
	end do
      end if
      call MPI_GROUP_EXCL(GROUP_WORLD,idle,ranks,GROUP_WORKERS, ierr)
      call MPI_COMM_CREATE(MPI_COMM_WORLD,
     &              GROUP_WORKERS,MPI_COMM_WORKERS,ierr) 

      if (HWorker .or. VWorker) then
        call MPI_COMM_SIZE( MPI_COMM_WORKERS, Nw, ierr )
        call MPI_COMM_RANK( MPI_COMM_WORKERS, wid, ierr )
c        print *, 'Process ', myid, ' of ', Nprocs, ' is worker ',
c     &          wid,' of ',Nw
      endif

c-------------------------------------------------------------------------------
c    Define 4D hslice data types
c-------------------------------------------------------------------------------
      call MPI_TYPE_VECTOR(is, ix*iy, ix*iy*iz, MPI_REAL, 
     &                      GLOBAL_4D_HSLICE, ierr)
      call MPI_TYPE_COMMIT(GLOBAL_4D_HSLICE, ierr)
      call MPI_TYPE_VECTOR(is, ix*iy, ix*iy*izloc, MPI_REAL, 
     &                      LOCAL_4D_HSLICE, ierr)
      call MPI_TYPE_COMMIT(LOCAL_4D_HSLICE, ierr)

      types(1) = GLOBAL_4D_HSLICE
      call MPI_TYPE_STRUCT(2,blength,tdispls,
     &		types, GLOB_4D_HSLICE, ierr)
      call MPI_TYPE_COMMIT(GLOB_4D_HSLICE, ierr)
      
c-------------------------------------------------------------------------------
c    Define 4D vcolumn data types
c-------------------------------------------------------------------------------
      call MPI_TYPE_VECTOR(iz*is, 1, ix*iy, MPI_REAL, 
     &                      GLOBAL_4D_VCOL, ierr)
      call MPI_TYPE_COMMIT(GLOBAL_4D_VCOL, ierr)
      call MPI_TYPE_VECTOR(iz*is, 1, icloc, MPI_REAL, 
     &                      LOCAL_4D_VCOL, ierr)
      call MPI_TYPE_COMMIT(LOCAL_4D_VCOL, ierr)

      types(1) = GLOBAL_4D_VCOL
      call MPI_TYPE_STRUCT(2,blength,tdispls,
     &		types, GLOB_4D_VCOL, ierr)
      call MPI_TYPE_COMMIT(GLOB_4D_VCOL, ierr)

c-------------------------------------------------------------------------------
c    Define 4D hslice set data types (large and small, differing in one element)
c-------------------------------------------------------------------------------

c	easy version to put slices together
c	faster version...
c	if (.false.) then
      call MPI_TYPE_HVECTOR(izloc, 1, 
     &	NHworkers*sizeofreal*ix*iy, 
     &	GLOBAL_4D_HSLICE, GLOBAL_LARGE_4D_HSLICES, ierr)
      call MPI_TYPE_COMMIT(GLOBAL_LARGE_4D_HSLICES, ierr)
      call MPI_TYPE_HVECTOR(izloc, 1, 
     &	sizeofreal*ix*iy, 
     &	LOCAL_4D_HSLICE, LOCAL_LARGE_4D_HSLICES, ierr)
      call MPI_TYPE_COMMIT(LOCAL_LARGE_4D_HSLICES, ierr)

      call MPI_TYPE_HVECTOR(izloc-1, 1,
     &	NHworkers*sizeofreal*ix*iy, 
     &	GLOBAL_4D_HSLICE, GLOBAL_SMALL_4D_HSLICES, ierr)
      call MPI_TYPE_COMMIT(GLOBAL_SMALL_4D_HSLICES, ierr)
      call MPI_TYPE_HVECTOR(izloc-1, 1, 
     &	sizeofreal*ix*iy, 
     &	LOCAL_4D_HSLICE, LOCAL_SMALL_4D_HSLICES, ierr)
      call MPI_TYPE_COMMIT(LOCAL_SMALL_4D_HSLICES, ierr)

c	end if

c	array dimension to put slices together...
c	more difficult... and slower here...

      if (.false.) then

      call MPI_TYPE_VECTOR(izloc, ix*iy, NHworkers*ix*iy, MPI_REAL, 
     &	THREE_DGHl, ierr)
      call MPI_TYPE_HVECTOR(is, 1, 
     &	ix*iy*iz*sizeofreal, THREE_DGHl, 
     &	GLOBAL_LARGE_4D_HSLICES, ierr)
      call MPI_TYPE_COMMIT(GLOBAL_LARGE_4D_HSLICES, ierr)
      
      call MPI_TYPE_VECTOR(izloc, ix*iy, ix*iy, MPI_REAL, 
     &	THREE_DLHl, ierr)
      call MPI_TYPE_HVECTOR(is, 1, 
     & ix*iy*izloc*sizeofreal, THREE_DLHl,
     & 	LOCAL_LARGE_4D_HSLICES, ierr)
      call MPI_TYPE_COMMIT(LOCAL_LARGE_4D_HSLICES, ierr)
      
      
      call MPI_TYPE_VECTOR(izloc-1, ix*iy, NHworkers*ix*iy, MPI_REAL, 
     &	THREE_DGHs, ierr)
      call MPI_TYPE_HVECTOR(is, 1, 
     & ix*iy*iz*sizeofreal, THREE_DGHs,
     & 	GLOBAL_SMALL_4D_HSLICES, ierr)
      call MPI_TYPE_COMMIT(GLOBAL_SMALL_4D_HSLICES, ierr)
      
      call MPI_TYPE_VECTOR(izloc-1, ix*iy, ix*iy, MPI_REAL, 
     &	THREE_DLHs, ierr)
      call MPI_TYPE_HVECTOR(is, 1, 
     & ix*iy*izloc*sizeofreal, THREE_DLHs,
     & 	LOCAL_SMALL_4D_HSLICES, ierr)
      call MPI_TYPE_COMMIT(LOCAL_SMALL_4D_HSLICES, ierr)
      
      end if

      types(1) = GLOBAL_LARGE_4D_HSLICES
      call MPI_TYPE_STRUCT(2,blength,tdispls,
     &		types, GLOB_LARGE_4D_HSLICES, ierr)
      call MPI_TYPE_COMMIT(GLOB_LARGE_4D_HSLICES, ierr)
      types(1) = GLOBAL_SMALL_4D_HSLICES
      call MPI_TYPE_STRUCT(2,blength,tdispls,
     &		types, GLOB_SMALL_4D_HSLICES, ierr)
      call MPI_TYPE_COMMIT(GLOB_SMALL_4D_HSLICES, ierr)

c-------------------------------------------------------------------------------
c    Define 4D V-column set data types (large and small, differing in one element)
c-------------------------------------------------------------------------------

c	easy version to create slices together
c	if (.false.) then
      call MPI_TYPE_HVECTOR(icloc, 1,
     &	NVworkers*sizeofreal, 
     &	GLOBAL_4D_VCOL, GLOBAL_LARGE_4D_VCOLS, ierr)
      call MPI_TYPE_COMMIT(GLOBAL_LARGE_4D_VCOLS, ierr)
      call MPI_TYPE_HVECTOR(icloc, 1,
     &	sizeofreal, 
     &	LOCAL_4D_VCOL, LOCAL_LARGE_4D_VCOLS, ierr)
      call MPI_TYPE_COMMIT(LOCAL_LARGE_4D_VCOLS, ierr)

      call MPI_TYPE_HVECTOR(icloc-1, 1,
     &	NVworkers*sizeofreal, 
     &	GLOBAL_4D_VCOL, GLOBAL_SMALL_4D_VCOLS, ierr)
      call MPI_TYPE_COMMIT(GLOBAL_SMALL_4D_VCOLS, ierr)
      call MPI_TYPE_HVECTOR(icloc-1, 1,
     &	sizeofreal, 
     &	LOCAL_4D_VCOL, LOCAL_SMALL_4D_VCOLS, ierr)
      call MPI_TYPE_COMMIT(LOCAL_SMALL_4D_VCOLS, ierr)
c	end if
c	array version, going through the dimensions...
      if (.false.) then
      
      call MPI_TYPE_VECTOR(icloc, 1, NVworkers, MPI_REAL, 
     &	TWO_DGVl, ierr)
      call MPI_TYPE_HVECTOR(iz, 1, 
     &	ix*iy*sizeofreal, TWO_DGVl, THREE_DGVl, ierr)
      call MPI_TYPE_HVECTOR(is, 1,
     &	ix*iy*iz*sizeofreal, THREE_DGVl,
     &	GLOBAL_LARGE_4D_VCOLS, ierr)
      call MPI_TYPE_COMMIT(GLOBAL_LARGE_4D_VCOLS, ierr)

      call MPI_TYPE_VECTOR(icloc, 1, 1, MPI_REAL, 
     &	TWO_DLVl, ierr)
      call MPI_TYPE_HVECTOR(iz, 1, 
     &	icloc*sizeofreal, TWO_DLVl, THREE_DLVl, ierr)
      call MPI_TYPE_HVECTOR(is, 1,
     &	icloc*iz*sizeofreal, THREE_DLVl,
     &	LOCAL_LARGE_4D_VCOLS, ierr)
      call MPI_TYPE_COMMIT(LOCAL_LARGE_4D_VCOLS, ierr)
     
     
      call MPI_TYPE_VECTOR(icloc-1, 1, NVworkers, MPI_REAL, 
     &	TWO_DGVs, ierr)
      call MPI_TYPE_HVECTOR(iz, 1, 
     &	ix*iy*sizeofreal, TWO_DGVs, THREE_DGVs, ierr)
      call MPI_TYPE_HVECTOR(is, 1,
     &	ix*iy*iz*sizeofreal, THREE_DGVs,
     &	GLOBAL_SMALL_4D_VCOLS, ierr)
      call MPI_TYPE_COMMIT(GLOBAL_SMALL_4D_VCOLS, ierr)

      call MPI_TYPE_VECTOR(icloc-1, 1, 1, MPI_REAL, 
     &	TWO_DLVs, ierr)
      call MPI_TYPE_HVECTOR(iz, 1, 
     &	icloc*sizeofreal, TWO_DLVs, THREE_DLVs, ierr)
      call MPI_TYPE_HVECTOR(is, 1,
     &	icloc*iz*sizeofreal, THREE_DLVs,
     &	LOCAL_SMALL_4D_VCOLS, ierr)
      call MPI_TYPE_COMMIT(LOCAL_SMALL_4D_VCOLS, ierr)
      
      end if

      types(1) = GLOBAL_LARGE_4D_VCOLS
      call MPI_TYPE_STRUCT(2,blength,tdispls,
     &		types, GLOB_LARGE_4D_VCOLS, ierr)
      call MPI_TYPE_COMMIT(GLOB_LARGE_4D_VCOLS, ierr)
      types(1) = GLOBAL_SMALL_4D_VCOLS
      call MPI_TYPE_STRUCT(2,blength,tdispls,
     &		types, GLOB_SMALL_4D_VCOLS, ierr)
      call MPI_TYPE_COMMIT(GLOB_SMALL_4D_VCOLS, ierr)
    
c-------------------------------------------------------------------------------
c    Define data types for shuffeling
c-------------------------------------------------------------------------------


c Below version creates the array of h-slice/v-column-parts 
c  going through each dimension one after the other.
      call MPI_TYPE_VECTOR(icloc, 1, NVworkers, MPI_REAL, 
     &                      LOCAL_2D_LV_PART_HCOL, ierr)
      call MPI_TYPE_HVECTOR(izloc, 1, ix*iy*sizeofreal, 
     &	LOCAL_2D_LV_PART_HCOL, LOCAL_3D_LHLV_PART_HCOL, ierr)
      call MPI_TYPE_HVECTOR(is, 1, ix*iy*izloc*sizeofreal, 
     &	LOCAL_3D_LHLV_PART_HCOL, LOCAL_4D_LHLV_HCOLS, ierr)
      call MPI_TYPE_COMMIT(LOCAL_4D_LHLV_HCOLS, ierr)
      
      call MPI_TYPE_HVECTOR(izloc-1, 1, ix*iy*sizeofreal, 
     &	LOCAL_2D_LV_PART_HCOL, LOCAL_3D_SHLV_PART_HCOL, ierr)
      call MPI_TYPE_HVECTOR(is, 1, ix*iy*izloc*sizeofreal, 
     &	LOCAL_3D_SHLV_PART_HCOL, LOCAL_4D_SHLV_HCOLS, ierr)
      call MPI_TYPE_COMMIT(LOCAL_4D_SHLV_HCOLS, ierr)

      call MPI_TYPE_VECTOR(icloc-1, 1, NVworkers, MPI_REAL, 
     &                      LOCAL_2D_SV_PART_HCOL, ierr)
      call MPI_TYPE_HVECTOR(izloc, 1, ix*iy*sizeofreal, 
     &	LOCAL_2D_SV_PART_HCOL, LOCAL_3D_LHSV_PART_HCOL, ierr)
      call MPI_TYPE_HVECTOR(is, 1, ix*iy*izloc*sizeofreal, 
     &	LOCAL_3D_LHSV_PART_HCOL, LOCAL_4D_LHSV_HCOLS, ierr)
      call MPI_TYPE_COMMIT(LOCAL_4D_LHSV_HCOLS, ierr)

      call MPI_TYPE_HVECTOR(izloc-1, 1, ix*iy*sizeofreal, 
     &	LOCAL_2D_SV_PART_HCOL, LOCAL_3D_SHSV_PART_HCOL, ierr)
      call MPI_TYPE_HVECTOR(is, 1, ix*iy*izloc*sizeofreal, 
     &	LOCAL_3D_SHSV_PART_HCOL, LOCAL_4D_SHSV_HCOLS, ierr)
      call MPI_TYPE_COMMIT(LOCAL_4D_SHSV_HCOLS, ierr)
      
      call MPI_TYPE_VECTOR(icloc, 1, 1, MPI_REAL, 
     &			LOCAL_2D_LV_PART_VCOL, ierr)
      call MPI_TYPE_HVECTOR(izloc, 1, 
     &	NHworkers*icloc*sizeofreal, 
     &	LOCAL_2D_LV_PART_VCOL, LOCAL_3D_LHLV_PART_VCOL, ierr)
      call MPI_TYPE_HVECTOR(is, 1, icloc*iz*sizeofreal, 
     &	LOCAL_3D_LHLV_PART_VCOL, LOCAL_4D_LHLV_VCOLS, ierr)
      call MPI_TYPE_COMMIT(LOCAL_4D_LHLV_VCOLS, ierr)
      
      call MPI_TYPE_HVECTOR(izloc-1, 1, 
     &	NHworkers*icloc*sizeofreal, 
     &	LOCAL_2D_LV_PART_VCOL, LOCAL_3D_SHLV_PART_VCOL, ierr)
      call MPI_TYPE_HVECTOR(is, 1, icloc*iz*sizeofreal, 
     &	LOCAL_3D_SHLV_PART_VCOL, LOCAL_4D_SHLV_VCOLS, ierr)
      call MPI_TYPE_COMMIT(LOCAL_4D_SHLV_VCOLS, ierr)

      call MPI_TYPE_VECTOR(icloc-1, 1, 1, MPI_REAL, 
     &			LOCAL_2D_SV_PART_VCOL, ierr)
      call MPI_TYPE_HVECTOR(izloc, 1, 
     &	NHworkers*icloc*sizeofreal, 
     &	LOCAL_2D_SV_PART_VCOL, LOCAL_3D_LHSV_PART_VCOL, ierr)
      call MPI_TYPE_HVECTOR(is, 1, icloc*iz*sizeofreal, 
     &	LOCAL_3D_LHSV_PART_VCOL, LOCAL_4D_LHSV_VCOLS, ierr)
      call MPI_TYPE_COMMIT(LOCAL_4D_LHSV_VCOLS, ierr)
      
      call MPI_TYPE_HVECTOR(izloc-1, 1, 
     &	NHworkers*icloc*sizeofreal, 
     &	LOCAL_2D_SV_PART_VCOL, LOCAL_3D_SHSV_PART_VCOL, ierr)
      call MPI_TYPE_HVECTOR(is, 1, icloc*iz*sizeofreal, 
     &	LOCAL_3D_SHSV_PART_VCOL, LOCAL_4D_SHSV_VCOLS, ierr)
      call MPI_TYPE_COMMIT(LOCAL_4D_SHSV_VCOLS, ierr)

c-------------------------------------------------------------------------------
c    Define 2D data types
c-------------------------------------------------------------------------------
      call MPI_TYPE_VECTOR(icloc, 1, NVworkers, MPI_REAL, 
     &		GLOBAL_LARGE_2D_VCOLS, ierr)
      call MPI_TYPE_COMMIT(GLOBAL_LARGE_2D_VCOLS, ierr)

      call MPI_TYPE_VECTOR(icloc-1, 1, NVworkers, MPI_REAL,
     &		GLOBAL_SMALL_2D_VCOLS, ierr)
      call MPI_TYPE_COMMIT(GLOBAL_SMALL_2D_VCOLS, ierr)

      call MPI_TYPE_VECTOR(icloc, 1, 1, MPI_REAL,
     &		LOCAL_LARGE_2D_VCOLS, ierr)
      call MPI_TYPE_COMMIT(LOCAL_LARGE_2D_VCOLS, ierr)
      call MPI_TYPE_VECTOR(icloc-1, 1, 1, MPI_REAL,
     &		LOCAL_SMALL_2D_VCOLS, ierr)
      call MPI_TYPE_COMMIT(LOCAL_SMALL_2D_VCOLS, ierr)

      types(1) = GLOBAL_LARGE_2D_VCOLS
      call MPI_TYPE_STRUCT(2,blength,tdispls,
     &		types, GLOB_LARGE_2D_VCOLS, ierr)
      call MPI_TYPE_COMMIT(GLOB_LARGE_2D_VCOLS, ierr)
      types(1) = GLOBAL_SMALL_2D_VCOLS
      call MPI_TYPE_STRUCT(2,blength,tdispls,
     &		types, GLOB_SMALL_2D_VCOLS, ierr)
      call MPI_TYPE_COMMIT(GLOB_SMALL_2D_VCOLS, ierr)

c-------------------------------------------------------------------------------
c    Define 2DN (column) data types
c-------------------------------------------------------------------------------
      call MPI_TYPE_VECTOR(is, 1, ix*iy, MPI_REAL,
     &				GLOBAL_2DN_VCOL, ierr)
      call MPI_TYPE_COMMIT(GLOBAL_2DN_VCOL, ierr)
      call MPI_TYPE_VECTOR(is, 1, icloc, MPI_REAL,
     &				LOCAL_2DN_VCOL, ierr)
      call MPI_TYPE_COMMIT(LOCAL_2DN_VCOL, ierr)

      types(1) = GLOBAL_2DN_VCOL
      call MPI_TYPE_STRUCT(2,blength,tdispls,
     &		types, GLOB_2DN_VCOL, ierr)
      call MPI_TYPE_COMMIT(GLOB_2DN_VCOL, ierr)

      call MPI_TYPE_HVECTOR(icloc, 1, 
     &		NVworkers*sizeofreal, GLOBAL_2DN_VCOL,
     &		GLOBAL_LARGE_2DN_VCOLS, ierr)
      call MPI_TYPE_COMMIT(GLOBAL_LARGE_2DN_VCOLS, ierr)
      call MPI_TYPE_HVECTOR(icloc-1, 1, 
     &		NVworkers*sizeofreal, GLOBAL_2DN_VCOL,
     &		GLOBAL_SMALL_2DN_VCOLS, ierr)
      call MPI_TYPE_COMMIT(GLOBAL_SMALL_2DN_VCOLS, ierr)

      call MPI_TYPE_HVECTOR(icloc, 1, 
     &		sizeofreal, LOCAL_2DN_VCOL,
     &		LOCAL_LARGE_2DN_VCOLS, ierr)
      call MPI_TYPE_COMMIT(LOCAL_LARGE_2DN_VCOLS, ierr)
      call MPI_TYPE_HVECTOR(icloc-1, 1, 
     &		sizeofreal, LOCAL_2DN_VCOL,
     &		LOCAL_SMALL_2DN_VCOLS, ierr)
      call MPI_TYPE_COMMIT(LOCAL_SMALL_2DN_VCOLS, ierr)

      types(1) = GLOBAL_LARGE_2DN_VCOLS
      call MPI_TYPE_STRUCT(2,blength,tdispls,
     &		types, GLOB_LARGE_2DN_VCOLS, ierr)
      call MPI_TYPE_COMMIT(GLOB_LARGE_2DN_VCOLS, ierr)
      types(1) = GLOBAL_SMALL_2DN_VCOLS
      call MPI_TYPE_STRUCT(2,blength,tdispls,
     &		types, GLOB_SMALL_2DN_VCOLS, ierr)
      call MPI_TYPE_COMMIT(GLOB_SMALL_2DN_VCOLS, ierr)

c-------------------------------------------------------------------------------
c    Define 3D (column) data types
c-------------------------------------------------------------------------------
      call MPI_TYPE_VECTOR(iz, 1, ix*iy, MPI_REAL,
     &				GLOBAL_3D_VCOL, ierr)
      call MPI_TYPE_COMMIT(GLOBAL_3D_VCOL, ierr)
      call MPI_TYPE_VECTOR(iz, 1, icloc, MPI_REAL,
     &				LOCAL_3D_VCOL, ierr)
      call MPI_TYPE_COMMIT(LOCAL_3D_VCOL, ierr)

      types(1) = GLOBAL_3D_VCOL
      call MPI_TYPE_STRUCT(2,blength,tdispls,
     &		types, GLOB_3D_VCOL, ierr)
      call MPI_TYPE_COMMIT(GLOB_3D_VCOL, ierr)

      call MPI_TYPE_VECTOR(izloc, ix*iy, NHworkers*ix*iy, 
     &			MPI_REAL, GLOBAL_LARGE_3D_HSLICES, ierr)
      call MPI_TYPE_COMMIT(GLOBAL_LARGE_3D_HSLICES, ierr)
      call MPI_TYPE_VECTOR(izloc, ix*iy, ix*iy, MPI_REAL,
     &				LOCAL_LARGE_3D_HSLICES, ierr)
      call MPI_TYPE_COMMIT(LOCAL_LARGE_3D_HSLICES, ierr)

      call MPI_TYPE_VECTOR(izloc-1, ix*iy, NHworkers*ix*iy, 
     &			MPI_REAL, GLOBAL_SMALL_3D_HSLICES, ierr)
      call MPI_TYPE_COMMIT(GLOBAL_LARGE_3D_HSLICES, ierr)
      call MPI_TYPE_VECTOR(izloc-1, ix*iy, ix*iy, MPI_REAL,
     &				LOCAL_SMALL_3D_HSLICES, ierr)
      call MPI_TYPE_COMMIT(LOCAL_LARGE_3D_HSLICES, ierr)

      call MPI_TYPE_HVECTOR(icloc, 1, 
     &		NVworkers*sizeofreal, GLOBAL_3D_VCOL,
     &		GLOBAL_LARGE_3D_VCOLS, ierr)
      call MPI_TYPE_COMMIT(GLOBAL_LARGE_3D_VCOLS, ierr)
      call MPI_TYPE_HVECTOR(icloc-1, 1, 
     &		NVworkers*sizeofreal, GLOBAL_3D_VCOL,
     &		GLOBAL_SMALL_3D_VCOLS, ierr)
      call MPI_TYPE_COMMIT(GLOBAL_SMALL_3D_VCOLS, ierr)

      call MPI_TYPE_HVECTOR(icloc, 1, 
     &		sizeofreal, LOCAL_3D_VCOL,
     &		LOCAL_LARGE_3D_VCOLS, ierr)
      call MPI_TYPE_COMMIT(LOCAL_LARGE_3D_VCOLS, ierr)
      call MPI_TYPE_HVECTOR(icloc-1, 1, 
     &		sizeofreal, LOCAL_3D_VCOL,
     &		LOCAL_SMALL_3D_VCOLS, ierr)
      call MPI_TYPE_COMMIT(LOCAL_SMALL_3D_VCOLS, ierr)

      types(1) = GLOBAL_LARGE_3D_HSLICES
      call MPI_TYPE_STRUCT(2,blength,tdispls,
     &		types, GLOB_LARGE_3D_HSLICES, ierr)
      call MPI_TYPE_COMMIT(GLOB_LARGE_3D_HSLICES, ierr)
      types(1) = GLOBAL_SMALL_3D_HSLICES
      call MPI_TYPE_STRUCT(2,blength,tdispls,
     &		types, GLOB_SMALL_3D_HSLICES, ierr)
      call MPI_TYPE_COMMIT(GLOB_SMALL_3D_HSLICES, ierr)
      types(1) = GLOBAL_LARGE_3D_VCOLS
      call MPI_TYPE_STRUCT(2,blength,tdispls,
     &		types, GLOB_LARGE_3D_VCOLS, ierr)
      call MPI_TYPE_COMMIT(GLOB_LARGE_3D_VCOLS, ierr)
      types(1) = GLOBAL_SMALL_3D_VCOLS
      call MPI_TYPE_STRUCT(2,blength,tdispls,
     &		types, GLOB_SMALL_3D_VCOLS, ierr)
      call MPI_TYPE_COMMIT(GLOB_SMALL_3D_VCOLS, ierr)

c-------------------------------------------------------------------------------
c    Define BD hslice data types (X and Y)
c-------------------------------------------------------------------------------
      call MPI_TYPE_VECTOR(2*is, ix, ix*iz, MPI_REAL, 
     &                      GLOBAL_BD_XHSLICE, ierr)
      call MPI_TYPE_COMMIT(GLOBAL_BD_XHSLICE, ierr)
      call MPI_TYPE_VECTOR(2*is, ix, ix*izloc, MPI_REAL, 
     &                      LOCAL_BD_XHSLICE, ierr)
      call MPI_TYPE_COMMIT(LOCAL_BD_XHSLICE, ierr)
      
      types(1) = GLOBAL_BD_XHSLICE
      call MPI_TYPE_STRUCT(2,blength,tdispls,
     &		types, GLOB_BD_XHSLICE, ierr)
      call MPI_TYPE_COMMIT(GLOB_BD_XHSLICE, ierr)

      call MPI_TYPE_VECTOR(2*is, iy, iy*iz, MPI_REAL, 
     &                      GLOBAL_BD_YHSLICE, ierr)
      call MPI_TYPE_COMMIT(GLOBAL_BD_YHSLICE, ierr)
      call MPI_TYPE_VECTOR(2*is, iy, iy*izloc, MPI_REAL, 
     &                      LOCAL_BD_YHSLICE, ierr)
      call MPI_TYPE_COMMIT(LOCAL_BD_YHSLICE, ierr)

      types(1) = GLOBAL_BD_YHSLICE
      call MPI_TYPE_STRUCT(2,blength,tdispls,
     &		types, GLOB_BD_YHSLICE, ierr)
      call MPI_TYPE_COMMIT(GLOB_BD_YHSLICE, ierr)

      call MPI_TYPE_HVECTOR(izloc, 1, 
     &		NHworkers*sizeofreal*ix, GLOBAL_BD_XHSLICE,
     &		GLOBAL_LARGE_BD_XHSLICES, ierr)
      call MPI_TYPE_COMMIT(GLOBAL_LARGE_BD_XHSLICES, ierr)
      call MPI_TYPE_HVECTOR(izloc-1, 1, 
     &		NHworkers*sizeofreal*ix, GLOBAL_BD_XHSLICE,
     &		GLOBAL_SMALL_BD_XHSLICES, ierr)
      call MPI_TYPE_COMMIT(GLOBAL_SMALL_BD_XHSLICES, ierr)

      call MPI_TYPE_HVECTOR(izloc, 1, 
     &		sizeofreal*ix, LOCAL_BD_XHSLICE,
     &		LOCAL_LARGE_BD_XHSLICES, ierr)
      call MPI_TYPE_COMMIT(LOCAL_LARGE_BD_XHSLICES, ierr)
      call MPI_TYPE_HVECTOR(izloc-1, 1, 
     &		sizeofreal*ix, LOCAL_BD_XHSLICE,
     &		LOCAL_SMALL_BD_XHSLICES, ierr)
      call MPI_TYPE_COMMIT(LOCAL_SMALL_BD_XHSLICES, ierr)

      call MPI_TYPE_HVECTOR(izloc, 1, 
     &		NHworkers*sizeofreal*iy, GLOBAL_BD_YHSLICE,
     &		GLOBAL_LARGE_BD_YHSLICES, ierr)
      call MPI_TYPE_COMMIT(GLOBAL_LARGE_BD_YHSLICES, ierr)
      call MPI_TYPE_HVECTOR(izloc-1, 1, 
     &		NHworkers*sizeofreal*iy, GLOBAL_BD_YHSLICE,
     &		GLOBAL_SMALL_BD_YHSLICES, ierr)
      call MPI_TYPE_COMMIT(GLOBAL_SMALL_BD_YHSLICES, ierr)

      call MPI_TYPE_HVECTOR(izloc, 1, 
     &		sizeofreal*iy, LOCAL_BD_YHSLICE,
     &		LOCAL_LARGE_BD_YHSLICES, ierr)
      call MPI_TYPE_COMMIT(LOCAL_LARGE_BD_YHSLICES, ierr)
      call MPI_TYPE_HVECTOR(izloc-1, 1, 
     &		sizeofreal*iy, LOCAL_BD_YHSLICE,
     &		LOCAL_SMALL_BD_YHSLICES, ierr)
      call MPI_TYPE_COMMIT(LOCAL_SMALL_BD_YHSLICES, ierr)

      types(1) = GLOBAL_LARGE_BD_XHSLICES
      call MPI_TYPE_STRUCT(2,blength,tdispls,
     &		types, GLOB_LARGE_BD_XHSLICES, ierr)
      call MPI_TYPE_COMMIT(GLOB_LARGE_BD_XHSLICES, ierr)
      types(1) = GLOBAL_SMALL_BD_XHSLICES
      call MPI_TYPE_STRUCT(2,blength,tdispls,
     &		types, GLOB_SMALL_BD_XHSLICES, ierr)
      call MPI_TYPE_COMMIT(GLOB_SMALL_BD_XHSLICES, ierr)
      types(1) = GLOBAL_LARGE_BD_YHSLICES
      call MPI_TYPE_STRUCT(2,blength,tdispls,
     &		types, GLOB_LARGE_BD_YHSLICES, ierr)
      call MPI_TYPE_COMMIT(GLOB_LARGE_BD_YHSLICES, ierr)
      types(1) = GLOBAL_SMALL_BD_YHSLICES
      call MPI_TYPE_STRUCT(2,blength,tdispls,
     &		types, GLOB_SMALL_BD_YHSLICES, ierr)
      call MPI_TYPE_COMMIT(GLOB_SMALL_BD_YHSLICES, ierr)

      end subroutine CreateCommDataTypes 

      end module HVCommDataTypes
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
c	distrib_x_2DN	Distribution of the 3D (x,y,Ns) arrays in
c			x-slices from the master to the workers
c	distrib_y_2DN	Distribution of the 3D (x,y,Ns) arrays in
c			y-slices from the master to the workers
c	distrib_x_3D	Distribution of the 3D (x,y,z) arrays in
c			x-slices from the master to the workers
c	distrib_y_3D	Distribution of the 3D (x,y,z) arrays in
c			y-slices from the master to the workers
c	distrib_x_BD	Distribution of the boundary data arrays in
c			x-slices from the master to the workers
c	distrib_y_BD	Distribution of the boundary data arrays in
c			y-slices from the master to the workers
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
      integer,pointer,dimension(:,:) :: status
      integer,pointer,dimension(:) :: request
      
c      
      if (Master) then
        allocate( status(MPI_STATUS_SIZE,iy),STAT=ierr )
        allocate( request(iy),STAT=ierr )      
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
      integer,pointer,dimension(:,:) :: status
      integer,pointer,dimension(:) :: request

      if (Master) then
        allocate( status(MPI_STATUS_SIZE,ix),STAT=ierr )
        allocate( request(ix),STAT=ierr )      
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
c  Master distributes a 3D (ix.iy.is) array of data in x-slice format 
c   to all workers
c-------------------------------------------------------------------------
c      
      subroutine distrib_x_2DN_v1(ix,iy,iyloc,Ns,s,s_x)

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
        call MPI_SEND(s(1,j,1), 1, GLOBAL_2DN_XSLICE, 
     &         owner_of_xslice(j), j,  MPI_COMM_WORLD, ierr)
	end do
      else if (XWorker) then
        do j=1, no_of_xslices(myid)
        call MPI_IRECV(s_x(1,j,1), 1, LOCAL_2DN_XSLICE, 0, 
     &			global_xslice_id(myid,j), 
     &          	MPI_COMM_WORLD, request(j), ierr)
	end do
	call MPI_WAITALL(no_of_xslices(myid), request, status, ierr)
      end if
c
      end subroutine  distrib_x_2DN_v1
      
c-------------------------------------------------------------------------
c  Master distributes a 3D (ix.iy.is) array of data in y-slice format 
c   to all workers
c-------------------------------------------------------------------------
c      
      subroutine distrib_y_2DN_v1(ix,iy,ixloc,Ns,s,s_y)

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
        call MPI_SEND(s(i,1,1), 1, GLOBAL_2DN_YSLICE, 
     &         owner_of_yslice(i), i,  MPI_COMM_WORLD, ierr)
	end do
      else if (YWorker) then
        do i=1, no_of_yslices(myid)
        call MPI_IRECV(s_y(i,1,1), 1, LOCAL_2DN_YSLICE, 0, 
     &			global_yslice_id(myid,i), 
     &                 MPI_COMM_WORLD, request(i), ierr)
	end do
	call MPI_WAITALL(no_of_yslices(myid), request, status, ierr)
      end if
c
      end subroutine  distrib_y_2DN_v1



              
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
c  Master distributes the x-boundary array 
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
c  Master distributes the y-boundary array 
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
      integer,pointer,dimension(:,:) :: status
      integer,pointer,dimension(:) :: request
c      
c The master receives a the x-slice array from each of the workers
c  and places it right in the global array, depending on the size
      if (Master) then
        allocate( status(MPI_STATUS_SIZE,NXworkers),STAT=ierr )
        allocate( request(NXworkers),STAT=ierr )      
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
      integer,pointer,dimension(:,:) :: status
      integer,pointer,dimension(:) :: request
c      
c The master receives the y-slice array from each yworker and places
c  it in the right position, depending on the size
      if (Master) then
        allocate( status(MPI_STATUS_SIZE,NYworkers),STAT=ierr )
        allocate( request(NYworkers),STAT=ierr )      
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
c  Master distributes a 3D (ix.iy.is) array of data in x-slice format 
c   to all xworkers
c-------------------------------------------------------------------------
c      
      subroutine distrib_x_2DN_v2(ix,iy,iyloc,Ns,s,s_x)

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
            call MPI_SEND(s(1,j,1), 1, GLOBAL_LARGE_2DN_XSLICES, 
     &         owner_of_xslice(j), j,  MPI_COMM_WORLD, ierr)
          else
            call MPI_SEND(s(1,j,1), 1, GLOBAL_SMALL_2DN_XSLICES, 
     &         owner_of_xslice(j), j,  MPI_COMM_WORLD, ierr)
          end if
	end do
      else if(XWorker) then
c The xworkers receive the x-slice arrays from the master in one
c  message depending on the array size
        if (no_of_xslices(myid) .eq. iyloc) then
          call MPI_RECV(s_x(1,1,1), 1, LOCAL_LARGE_2DN_XSLICES, 0, 
     &			myid, MPI_COMM_WORLD, status, ierr)
        else
          call MPI_RECV(s_x(1,1,1), 1, LOCAL_SMALL_2DN_XSLICES, 0, 
     &			myid, MPI_COMM_WORLD, status, ierr)
	end if	  
      end if
c
      end subroutine  distrib_x_2DN_v2
      
c-------------------------------------------------------------------------
c  Master distributes a 3D (ix.iy.is) array of data in y-slice format 
c   to all yworkers
c-------------------------------------------------------------------------
c      
      subroutine distrib_y_2DN_v2(ix,iy,ixloc,Ns,s,s_y)

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
            call MPI_SEND(s(i,1,1), 1, GLOBAL_LARGE_2DN_YSLICES, 
     &         owner_of_yslice(i), i,  MPI_COMM_WORLD, ierr)
          else
            call MPI_SEND(s(i,1,1), 1, GLOBAL_SMALL_2DN_YSLICES, 
     &         owner_of_yslice(i), i,  MPI_COMM_WORLD, ierr)
          end if
	end do
      else if (YWorker) then
c The yworkers receive their y-slice array at once depending on
c  the size
        if (no_of_yslices(myid) .eq. ixloc) then
          call MPI_RECV(s_y(1,1,1), 1, LOCAL_LARGE_2DN_YSLICES, 0, 
     &			myid, MPI_COMM_WORLD, status, ierr)
	else
          call MPI_RECV(s_y(1,1,1), 1, LOCAL_SMALL_2DN_YSLICES, 0, 
     &			myid, MPI_COMM_WORLD, status, ierr)
	end if
      end if
c
      end subroutine  distrib_y_2DN_v2



              
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
c  Master distributes the x-boundary array 
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
c  Master distributes the y-boundary array 
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
c  Master distributes a 2DN array of data in x-slice format to all workers
c  using scatterv
c-----------------------------------------------------------------
      subroutine distrib_x_2DN_v3(ix,iy,iyloc,Ns,s,s_x)
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
     &		GLOB_LARGE_2DN_XSLICES, s_x(1,1,1),
     &		recvcount, LOCAL_LARGE_2DN_XSLICES, 0,
     &		MPI_COMM_WORLD, ierr)
      else
c Otherwise we need 2 calls.
c First everybody gets the minimal size everybody has
        call MPI_SCATTERV(s(1,1,1), sendcounts, displs, 
     &		GLOB_SMALL_2DN_XSLICES, s_x(1,1,1),
     &		recvcount, LOCAL_SMALL_2DN_XSLICES, 0,
     &		MPI_COMM_WORLD, ierr)
	sendcounts((mod(iy,NXworkers)+2) : Nprocs) = 0
	displs((mod(iy,NXworkers)+2) : Nprocs) = 0
        if (MyId .gt. mod(iy,NXworkers)) then
          recvcount = 0
        end if
c Then the rest processors get their further slices (max. one more)
        call MPI_SCATTERV(s(1,NXworkers*(iyloc-1)+1,1), sendcounts,  
     &		displs, GLOB_2DN_XSLICE, s_x(1,iyloc,1),
     &		recvcount, LOCAL_2DN_XSLICE, 0,
     &		MPI_COMM_WORLD, ierr)
      end if

      end subroutine  distrib_x_2DN_v3

      
c-----------------------------------------------------------------
c  Master distributes a 4D array of data in y-slice format to all workers
c  using scatterv
c-----------------------------------------------------------------
      subroutine distrib_y_2DN_v3(ix,iy,ixloc,Ns,s,s_y)
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
     &		GLOB_LARGE_2DN_YSLICES, s_y(1,1,1),
     &		recvcount, LOCAL_LARGE_2DN_YSLICES, 0,
     &		MPI_COMM_WORLD, ierr)
      else
c The matrix has to be scattered in 2 calls
c First everybody gets the minimal size everybody has
        call MPI_SCATTERV(s(1,1,1), sendcounts, displs, 
     &		GLOB_SMALL_2DN_YSLICES, s_y(1,1,1),
     &		recvcount, LOCAL_SMALL_2DN_YSLICES, 0,
     &		MPI_COMM_WORLD, ierr)
	sendcounts((mod(ix,NYworkers)+2) : Nprocs) = 0
	displs((mod(ix,NYworkers)+2) : Nprocs) = 0
        if (MyId .gt. mod(ix,NYworkers)) then
          recvcount = 0
        end if
c then the rest processors get their further slices (max. one more)        
	call MPI_SCATTERV(s(NYworkers*(ixloc-1)+1,1,1), sendcounts,  
     &		displs, GLOB_2DN_YSLICE, s_y(ixloc,1,1),
     &		recvcount, LOCAL_2DN_YSLICE, 0,
     &		MPI_COMM_WORLD, ierr)
      end if

      end subroutine  distrib_y_2DN_v3

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
c  Master distributes a BD array of data in x-slice format to all workers
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
c  Master distributes a BD array of data in y-slice format to all workers
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
c	distrib_v_2DN	Distribution of the 3D (x,y,Ns) arrays in
c			v-columns from the master to the workers
c	distrib_h_3D	Distribution of the 3D (x,y,z) arrays in
c			h-slices from the master to the workers
c	distrib_v_3D	Distribution of the 3D (x,y,z) arrays in
c			v-columns from the master to the workers
c	distrib_hx_BD	Distribution of the x-boundary data arrays in
c			h-slices from the master to the workers
c	distrib_hy_BD	Distribution of the y-boundary data arrays in
c			h-slices from the master to the workers
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
      integer,pointer,dimension(:,:) :: status
      integer,pointer,dimension(:) :: request
c     
      if (Master) then
        allocate( status(MPI_STATUS_SIZE,iz),STAT=ierr )
        allocate( request(iz),STAT=ierr )      
        do i=1,iz
        call MPI_IRECV(s(1,1,i,1), 1, GLOBAL_4D_HSLICE, 
     &		owner_of_hslice(i), i,  MPI_COMM_WORLD, 
     &         request(i), ierr)
	end do
	call MPI_WAITALL(iz, request, status, ierr)
	deallocate(status, STAT=ierr)
	deallocate(request, STAT=ierr)
      else if (HWorker) then
        do i=1, no_of_hslices(myid)
          call MPI_SEND(s_h(1,1,i,1), 1, LOCAL_4D_HSLICE, 
     &                 0, global_hslice_id(MyId,i), 
     &                 MPI_COMM_WORLD, ierr)
	end do
      end if
      
      end subroutine  gather_h_4D_v1
        
c-----------------------------------------------------------------
c  Master gathers a 4D array of data in v-slice format from all workers
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
      integer,pointer,dimension(:,:) :: status
      integer,pointer,dimension(:) :: request
c     
      if (Master) then
        allocate( status(MPI_STATUS_SIZE,ix*iy),STAT=ierr )
        allocate( request(ix*iy),STAT=ierr )      
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
      real, pointer,dimension(:,:,:,:)    :: sbuf
      real, pointer,dimension(:,:,:) :: rbuf

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
      real, pointer, dimension(:,:,:,:)    :: sbuf
      real, pointer, dimension(:,:,:)    :: rbuf

      if (VWorker) then
        if (.not.allocated(sbuf)) then
          allocate( sbuf(icloc,izloc,Ns,NHworkers), STAT=ierr )
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
c  to all workers - only needed for xy-transport
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
c  Master distributes a 2D (ix.iy) array of data in h-slice format 
c  to all workers - only needed for xy-transport
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
c  Master distributes a 3D (ix.iy.is) array of data in hslice format 
c   to all workers - only needed for v-transport
c-------------------------------------------------------------------------
c      
      subroutine distrib_v_2DN_v1(ix,iy,icloc,Ns,s,s_v)

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
     &		1, GLOBAL_2DN_VCOL, 
     &         owner_of_vcol(j), j,  MPI_COMM_WORLD, ierr)
	end do
      else if (VWorker) then
        do j=1, no_of_vcols(myid)
        call MPI_IRECV(s_v(1,j,1), 1, LOCAL_2DN_VCOL, 0, 
     &			global_vcol_id(myid,j), 
     &          	MPI_COMM_WORLD, request(j), ierr)
	end do
	call MPI_WAITALL(no_of_vcols(MyId), request, status, ierr)
      end if
c
      end subroutine  distrib_v_2DN_v1
      
c-------------------------------------------------------------------------
c  Master distributes a 3D (ix.iy.iz) array of data in x-slice format 
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
c  Master distributes a 3D (ix.iy.iz) array of data in y-slice format 
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
c  Master distributes the x-boundary array 
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
c  Master distributes the y-boundary array 
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
      integer,pointer,dimension(:,:) :: status
      integer,pointer,dimension(:) :: request
      
      if (Master) then
        allocate( status(MPI_STATUS_SIZE,NHworkers),STAT=ierr )
        allocate( request(NHworkers),STAT=ierr )      
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
      integer,pointer,dimension(:,:) :: status
      integer,pointer,dimension(:) :: request
      
      if (Master) then
        allocate( status(MPI_STATUS_SIZE,NVworkers),STAT=ierr )
        allocate( request(NVworkers),STAT=ierr )      
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
c  Master distributes a 3D (ix.iy.is) array of data in vcolumn format 
c   to all workers - only needed for v-transport
c-------------------------------------------------------------------------
c      
      subroutine distrib_v_2DN_v2(ix,iy,icloc,Ns,s,s_v)

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
            call MPI_SEND(s(j,1,1), 1, GLOBAL_LARGE_2DN_VCOLS, 
     &         owner_of_vcol(j), j,  MPI_COMM_WORLD, ierr)
          else
            call MPI_SEND(s(j,1,1), 1, GLOBAL_SMALL_2DN_VCOLS, 
     &         owner_of_vcol(j), j,  MPI_COMM_WORLD, ierr)
	  end if
	end do
      else if (VWorker) then
        if (no_of_vcols(myid) .eq. icloc) then
          call MPI_RECV(s_v(1,1,1), 1, LOCAL_LARGE_2DN_VCOLS, 0, 
     &			myid, MPI_COMM_WORLD, status, ierr)
        else
          call MPI_RECV(s_v(1,1,1), 1, LOCAL_SMALL_2DN_VCOLS, 0, 
     &			myid, MPI_COMM_WORLD, status, ierr)
	end if
      end if
c
      end subroutine  distrib_v_2DN_v2
      
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
c  Master distributes the x-boundary array 
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
c  Master distributes the y-boundary array 
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
      real, dimension(:,:,:,:), allocatable :: sbuf,rbuf

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
      real, dimension(:,:,:,:), allocatable    :: sbuf,rbuf

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
c  Master distributes a 2DN array of data in v-column format to all workers
c  using scatterv
c-----------------------------------------------------------------
      subroutine distrib_v_2DN_v3(ix,iy,icloc,Ns,s,s_v)
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
     &		GLOB_LARGE_2DN_VCOLS, s_v(1,1,1),
     &		recvcount, LOCAL_LARGE_2DN_VCOLS, 0,
     &		MPI_COMM_WORLD, ierr)
      else
        call MPI_SCATTERV(s(1,1,1), sendcounts, displs, 
     &		GLOB_SMALL_2DN_VCOLS, s_v(1,1,1),
     &		recvcount, LOCAL_SMALL_2DN_VCOLS, 0,
     &		MPI_COMM_WORLD, ierr)
	sendcounts((mod(ix*iy,NVworkers)+2) : Nprocs) = 0
	displs((mod(ix*iy,NVworkers)+2) : Nprocs) = 0
        if (MyId .gt. mod(ix*iy,NVworkers)) then
          recvcount = 0
        end if
        call MPI_SCATTERV(s(NVworkers*(icloc-1)+1,1,1), sendcounts,  
     &		displs, GLOB_2DN_VCOL, s_v(1,icloc,1),
     &		recvcount, LOCAL_2DN_VCOL, 0,
     &		MPI_COMM_WORLD, ierr)
      end if

      end subroutine  distrib_v_2DN_v3

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
c  Master distributes a 4D array of data in v-column format to all workers
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
c  Master distributes a 4D array of data in h-slice 
c     format to all workers
c   using scatterv
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
c  Master distributes a 4D array of data in h-slice 
c     format to all workers
c   using scatterv
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

      module XYParallelCommunication

      use XYParallelDataMap
      use XYCommunicationLibrary
      use XYParallelMemAlloc

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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      module HVParallelCommunication

      use HVParallelDataMap
      use HVCommunicationLibrary
      use HVParallelMemAlloc

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
