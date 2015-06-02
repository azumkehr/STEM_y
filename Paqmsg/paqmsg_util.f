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
c	for each type of X and Y-slice: 2D, BDz, 3D, BD, 4D, local, global
        integer :: GLOBAL_4D_XSLICE, LOCAL_4D_XSLICE
        integer :: GLOBAL_4D_YSLICE, LOCAL_4D_YSLICE
        integer :: GLOBAL_4D_COLUMN
        integer :: LOCAL_4D_XCOLUMN, LOCAL_4D_YCOLUMN
	integer :: GLOBAL_2D_YSLICE, LOCAL_2D_YSLICE
	integer :: GLOBAL_BDz_XSLICE, LOCAL_BDz_XSLICE
	integer :: GLOBAL_BDz_YSLICE, LOCAL_BDz_YSLICE
	integer :: GLOBAL_3D_XSLICE, LOCAL_3D_XSLICE
	integer :: GLOBAL_3D_YSLICE, LOCAL_3D_YSLICE
	integer :: GLOBAL_BD_XSLICE, LOCAL_BD_XSLICE
	integer :: GLOBAL_BD_YSLICE, LOCAL_BD_YSLICE
	integer :: GLOBAL_BDx_XSLICE, LOCAL_BDx_XSLICE
	integer :: GLOBAL_BDy_YSLICE, LOCAL_BDy_YSLICE
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
	integer :: GLOBAL_LARGE_BDz_XSLICES, LOCAL_LARGE_BDz_XSLICES
	integer :: GLOBAL_SMALL_BDz_XSLICES, LOCAL_SMALL_BDz_XSLICES
	integer :: GLOBAL_LARGE_BDz_YSLICES, LOCAL_LARGE_BDz_YSLICES
	integer :: GLOBAL_SMALL_BDz_YSLICES, LOCAL_SMALL_BDz_YSLICES
	integer :: GLOBAL_LARGE_3D_XSLICES, LOCAL_LARGE_3D_XSLICES
	integer :: GLOBAL_SMALL_3D_XSLICES, LOCAL_SMALL_3D_XSLICES
	integer :: GLOBAL_LARGE_3D_YSLICES, LOCAL_LARGE_3D_YSLICES
	integer :: GLOBAL_SMALL_3D_YSLICES, LOCAL_SMALL_3D_YSLICES
	integer :: GLOBAL_LARGE_BD_XSLICES, LOCAL_LARGE_BD_XSLICES
	integer :: GLOBAL_SMALL_BD_XSLICES, LOCAL_SMALL_BD_XSLICES
	integer :: GLOBAL_LARGE_BD_YSLICES, LOCAL_LARGE_BD_YSLICES
	integer :: GLOBAL_SMALL_BD_YSLICES, LOCAL_SMALL_BD_YSLICES
	integer :: GLOBAL_LARGE_BDx_XSLICES, LOCAL_LARGE_BDx_XSLICES
	integer :: GLOBAL_SMALL_BDx_XSLICES, LOCAL_SMALL_BDx_XSLICES
	integer :: GLOBAL_LARGE_BDy_YSLICES, LOCAL_LARGE_BDy_YSLICES
	integer :: GLOBAL_SMALL_BDy_YSLICES, LOCAL_SMALL_BDy_YSLICES
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
	integer :: GLOB_LARGE_BDz_XSLICES, GLOB_SMALL_BDz_XSLICES
	integer :: GLOB_BDz_XSLICE, GLOB_BDz_YSLICE
	integer :: GLOB_LARGE_BDz_YSLICES, GLOB_SMALL_BDz_YSLICES
	integer :: GLOB_LARGE_3D_XSLICES, GLOB_SMALL_3D_XSLICES
	integer :: GLOB_3D_XSLICE, GLOB_3D_YSLICE
	integer :: GLOB_LARGE_3D_YSLICES, GLOB_SMALL_3D_YSLICES
	integer :: GLOB_LARGE_BD_XSLICES, GLOB_SMALL_BD_XSLICES
	integer :: GLOB_BD_XSLICE, GLOB_BD_YSLICE
	integer :: GLOB_LARGE_BD_YSLICES, GLOB_SMALL_BD_YSLICES
	integer :: GLOB_LARGE_BDx_XSLICES, GLOB_SMALL_BDx_XSLICES
	integer :: GLOB_BDx_XSLICE, GLOB_BDy_YSLICE
	integer :: GLOB_LARGE_BDy_YSLICES, GLOB_SMALL_BDy_YSLICES

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
c    Define BDz data types
c-------------------------------------------------------------------------------
      call MPI_TYPE_VECTOR(is, ix, ix*iy, MPI_REAL,
     &				GLOBAL_BDz_XSLICE, ierr)
      call MPI_TYPE_COMMIT(GLOBAL_BDz_XSLICE, ierr)
      call MPI_TYPE_VECTOR(is, ix, ix*iyloc, MPI_REAL,
     &				LOCAL_BDz_XSLICE, ierr)
      call MPI_TYPE_COMMIT(LOCAL_BDz_XSLICE, ierr)

      call MPI_TYPE_VECTOR(iy*is, 1, ix, MPI_REAL,
     &				GLOBAL_BDz_YSLICE, ierr)
      call MPI_TYPE_COMMIT(GLOBAL_BDz_YSLICE, ierr)
      call MPI_TYPE_VECTOR(iy*is, 1, ixloc, MPI_REAL,
     &				LOCAL_BDz_YSLICE, ierr)
      call MPI_TYPE_COMMIT(LOCAL_BDz_YSLICE, ierr)

      types(1) = GLOBAL_BDz_XSLICE
      call MPI_TYPE_STRUCT(2,blength,tdispls,
     &		types, GLOB_BDz_XSLICE, ierr)
      call MPI_TYPE_COMMIT(GLOB_BDz_XSLICE, ierr)
      types(1) = GLOBAL_BDz_YSLICE
      call MPI_TYPE_STRUCT(2,blength,tdispls,
     &		types, GLOB_BDz_YSLICE, ierr)
      call MPI_TYPE_COMMIT(GLOB_BDz_YSLICE, ierr)

c This version creates the array of x/y-slices for BDz distribution
c  using the given x/y-slice structures. 
      call MPI_TYPE_HVECTOR(iyloc, 1, 
     &		NXworkers*sizeofreal*ix, GLOBAL_BDz_XSLICE,
     &		GLOBAL_LARGE_BDz_XSLICES, ierr)
      call MPI_TYPE_COMMIT(GLOBAL_LARGE_BDz_XSLICES, ierr)
      call MPI_TYPE_HVECTOR(iyloc-1, 1, 
     &		NXworkers*sizeofreal*ix, GLOBAL_BDz_XSLICE,
     &		GLOBAL_SMALL_BDz_XSLICES, ierr)
      call MPI_TYPE_COMMIT(GLOBAL_SMALL_BDz_XSLICES, ierr)

      call MPI_TYPE_HVECTOR(iyloc, 1, 
     &		sizeofreal*ix, LOCAL_BDz_XSLICE,
     &		LOCAL_LARGE_BDz_XSLICES, ierr)
      call MPI_TYPE_COMMIT(LOCAL_LARGE_BDz_XSLICES, ierr)
      call MPI_TYPE_HVECTOR(iyloc-1, 1, 
     &		sizeofreal*ix, LOCAL_BDz_XSLICE,
     &		LOCAL_SMALL_BDz_XSLICES, ierr)
      call MPI_TYPE_COMMIT(LOCAL_SMALL_BDz_XSLICES, ierr)

      call MPI_TYPE_HVECTOR(ixloc, 1, 
     &		NYworkers*sizeofreal, GLOBAL_BDz_YSLICE,
     &		GLOBAL_LARGE_BDz_YSLICES, ierr)
      call MPI_TYPE_COMMIT(GLOBAL_LARGE_BDz_YSLICES, ierr)
      call MPI_TYPE_HVECTOR(ixloc-1, 1, 
     &		NYworkers*sizeofreal, GLOBAL_BDz_YSLICE,
     &		GLOBAL_SMALL_BDz_YSLICES, ierr)
      call MPI_TYPE_COMMIT(GLOBAL_SMALL_BDz_YSLICES, ierr)

      call MPI_TYPE_HVECTOR(ixloc, 1, 
     &		sizeofreal, LOCAL_BDz_YSLICE,
     &		LOCAL_LARGE_BDz_YSLICES, ierr)
      call MPI_TYPE_COMMIT(LOCAL_LARGE_BDz_YSLICES, ierr)
      call MPI_TYPE_HVECTOR(ixloc-1, 1, 
     &		sizeofreal, LOCAL_BDz_YSLICE,
     &		LOCAL_SMALL_BDz_YSLICES, ierr)
      call MPI_TYPE_COMMIT(LOCAL_SMALL_BDz_YSLICES, ierr)

      types(1) = GLOBAL_LARGE_BDz_XSLICES
      call MPI_TYPE_STRUCT(2,blength,tdispls,
     &		types, GLOB_LARGE_BDz_XSLICES, ierr)
      call MPI_TYPE_COMMIT(GLOB_LARGE_BDz_XSLICES, ierr)
      types(1) = GLOBAL_SMALL_BDz_XSLICES
      call MPI_TYPE_STRUCT(2,blength,tdispls,
     &		types, GLOB_SMALL_BDz_XSLICES, ierr)
      call MPI_TYPE_COMMIT(GLOB_SMALL_BDz_XSLICES, ierr)
      types(1) = GLOBAL_LARGE_BDz_YSLICES
      call MPI_TYPE_STRUCT(2,blength,tdispls,
     &		types, GLOB_LARGE_BDz_YSLICES, ierr)
      call MPI_TYPE_COMMIT(GLOB_LARGE_BDz_YSLICES, ierr)      
      types(1) = GLOBAL_SMALL_BDz_YSLICES
      call MPI_TYPE_STRUCT(2,blength,tdispls,
     &		types, GLOB_SMALL_BDz_YSLICES, ierr)
      call MPI_TYPE_COMMIT(GLOB_SMALL_BDz_YSLICES, ierr)      

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
c    Define BDx,BDy data types
c-------------------------------------------------------------------------------
      call MPI_TYPE_VECTOR(iz*is, 1, iy, MPI_REAL,
     &				GLOBAL_BDx_XSLICE, ierr)
      call MPI_TYPE_COMMIT(GLOBAL_BDx_XSLICE, ierr)
      call MPI_TYPE_VECTOR(iz*is, 1, iyloc, MPI_REAL,
     &				LOCAL_BDx_XSLICE, ierr)
      call MPI_TYPE_COMMIT(LOCAL_BDx_XSLICE, ierr)

      call MPI_TYPE_VECTOR(iz*is, 1, ix, MPI_REAL,
     &				GLOBAL_BDy_YSLICE, ierr)
      call MPI_TYPE_COMMIT(GLOBAL_BDy_YSLICE, ierr)
      call MPI_TYPE_VECTOR(iz*is, 1, ixloc, MPI_REAL,
     &				LOCAL_BDy_YSLICE, ierr)
      call MPI_TYPE_COMMIT(LOCAL_BDy_YSLICE, ierr)

      types(1) = GLOBAL_BDx_XSLICE
      call MPI_TYPE_STRUCT(2,blength,tdispls,
     &		types, GLOB_BDx_XSLICE, ierr)
      call MPI_TYPE_COMMIT(GLOB_BDx_XSLICE, ierr)
      types(1) = GLOBAL_BDy_YSLICE
      call MPI_TYPE_STRUCT(2,blength,tdispls,
     &		types, GLOB_BDy_YSLICE, ierr)
      call MPI_TYPE_COMMIT(GLOB_BDy_YSLICE, ierr)

c This version creates the array of x/y-slices for BD distribution
c  using the given x/y-slice structures.
      call MPI_TYPE_HVECTOR(iyloc, 1, 
     &		NXworkers*sizeofreal, GLOBAL_BDx_XSLICE,
     &		GLOBAL_LARGE_BDx_XSLICES, ierr)
      call MPI_TYPE_COMMIT(GLOBAL_LARGE_BDx_XSLICES, ierr)
      call MPI_TYPE_HVECTOR(iyloc-1, 1, 
     &		NXworkers*sizeofreal, GLOBAL_BDx_XSLICE,
     &		GLOBAL_SMALL_BDx_XSLICES, ierr)
      call MPI_TYPE_COMMIT(GLOBAL_SMALL_BDx_XSLICES, ierr)

      call MPI_TYPE_HVECTOR(iyloc, 1, 
     &		sizeofreal, LOCAL_BDx_XSLICE,
     &		LOCAL_LARGE_BDx_XSLICES, ierr)
      call MPI_TYPE_COMMIT(LOCAL_LARGE_BDx_XSLICES, ierr)
      call MPI_TYPE_HVECTOR(iyloc-1, 1, 
     &		sizeofreal, LOCAL_BDx_XSLICE,
     &		LOCAL_SMALL_BDx_XSLICES, ierr)
      call MPI_TYPE_COMMIT(LOCAL_SMALL_BDx_XSLICES, ierr)

      call MPI_TYPE_HVECTOR(ixloc, 1, 
     &		NYworkers*sizeofreal, GLOBAL_BDy_YSLICE,
     &		GLOBAL_LARGE_BDy_YSLICES, ierr)
      call MPI_TYPE_COMMIT(GLOBAL_LARGE_BDy_YSLICES, ierr)
      call MPI_TYPE_HVECTOR(ixloc-1, 1, 
     &		NYworkers*sizeofreal, GLOBAL_BDy_YSLICE,
     &		GLOBAL_SMALL_BDy_YSLICES, ierr)
      call MPI_TYPE_COMMIT(GLOBAL_SMALL_BDy_YSLICES, ierr)

      call MPI_TYPE_HVECTOR(ixloc, 1, 
     &		sizeofreal, LOCAL_BDy_YSLICE,
     &		LOCAL_LARGE_BDy_YSLICES, ierr)
      call MPI_TYPE_COMMIT(LOCAL_LARGE_BDy_YSLICES, ierr)
      call MPI_TYPE_HVECTOR(ixloc-1, 1, 
     &		sizeofreal, LOCAL_BDy_YSLICE,
     &		LOCAL_SMALL_BDy_YSLICES, ierr)
      call MPI_TYPE_COMMIT(LOCAL_SMALL_BDy_YSLICES, ierr)

      types(1) = GLOBAL_LARGE_BDx_XSLICES
      call MPI_TYPE_STRUCT(2,blength,tdispls,
     &		types, GLOB_LARGE_BDx_XSLICES, ierr)
      call MPI_TYPE_COMMIT(GLOB_LARGE_BDx_XSLICES, ierr)
      types(1) = GLOBAL_SMALL_BDx_XSLICES
      call MPI_TYPE_STRUCT(2,blength,tdispls,
     &		types, GLOB_SMALL_BDx_XSLICES, ierr)
      call MPI_TYPE_COMMIT(GLOB_SMALL_BDx_XSLICES, ierr)
      types(1) = GLOBAL_LARGE_BDy_YSLICES
      call MPI_TYPE_STRUCT(2,blength,tdispls,
     &		types, GLOB_LARGE_BDy_YSLICES, ierr)
      call MPI_TYPE_COMMIT(GLOB_LARGE_BDy_YSLICES, ierr)      
      types(1) = GLOBAL_SMALL_BDy_YSLICES
      call MPI_TYPE_STRUCT(2,blength,tdispls,
     &		types, GLOB_SMALL_BDy_YSLICES, ierr)
      call MPI_TYPE_COMMIT(GLOB_SMALL_BDy_YSLICES, ierr)      

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
c		2D, BDz, 3D, BD, 4D, local, global
        integer :: GLOBAL_4D_VCOL, LOCAL_4D_VCOL
        integer :: GLOBAL_4D_HSLICE, LOCAL_4D_HSLICE
	integer :: GLOBAL_BDz_VCOL, LOCAL_BDz_VCOL
	integer :: GLOBAL_3D_VCOL, LOCAL_3D_VCOL
	integer :: GLOBAL_BD_XHSLICE, LOCAL_BD_XHSLICE
	integer :: GLOBAL_BD_YHSLICE, LOCAL_BD_YHSLICE
	integer :: GLOBAL_BDx_HSLICE, LOCAL_BDx_HSLICE
	integer :: GLOBAL_BDy_HSLICE, LOCAL_BDy_HSLICE
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
	integer :: GLOBAL_LARGE_BDz_VCOLS, LOCAL_LARGE_BDz_VCOLS
	integer :: GLOBAL_SMALL_BDz_VCOLS, LOCAL_SMALL_BDz_VCOLS
	integer :: GLOBAL_LARGE_3D_HSLICES, LOCAL_LARGE_3D_HSLICES
	integer :: GLOBAL_SMALL_3D_HSLICES, LOCAL_SMALL_3D_HSLICES
	integer :: GLOBAL_LARGE_3D_VCOLS, LOCAL_LARGE_3D_VCOLS
	integer :: GLOBAL_SMALL_3D_VCOLS, LOCAL_SMALL_3D_VCOLS
	integer :: GLOBAL_LARGE_BD_XHSLICES, LOCAL_LARGE_BD_XHSLICES
	integer :: GLOBAL_SMALL_BD_XHSLICES, LOCAL_SMALL_BD_XHSLICES
	integer :: GLOBAL_LARGE_BD_YHSLICES, LOCAL_LARGE_BD_YHSLICES
	integer :: GLOBAL_SMALL_BD_YHSLICES, LOCAL_SMALL_BD_YHSLICES
	integer :: GLOBAL_LARGE_BDx_HSLICES, LOCAL_LARGE_BDx_HSLICES
	integer :: GLOBAL_SMALL_BDx_HSLICES, LOCAL_SMALL_BDx_HSLICES
	integer :: GLOBAL_LARGE_BDy_HSLICES, LOCAL_LARGE_BDy_HSLICES
	integer :: GLOBAL_SMALL_BDy_HSLICES, LOCAL_SMALL_BDy_HSLICES
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
	integer :: GLOB_BDz_VCOL
	integer :: GLOB_LARGE_BDz_VCOLS, GLOB_SMALL_BDz_VCOLS
	integer :: GLOB_LARGE_3D_HSLICES, GLOB_SMALL_3D_HSLICES
	integer :: GLOB_3D_VCOL
	integer :: GLOB_LARGE_3D_VCOLS, GLOB_SMALL_3D_VCOLS
	integer :: GLOB_LARGE_BD_XHSLICES, GLOB_SMALL_BD_XHSLICES
	integer :: GLOB_BD_XHSLICE, GLOB_BD_YHSLICE
	integer :: GLOB_LARGE_BD_YHSLICES, GLOB_SMALL_BD_YHSLICES
	integer :: GLOB_LARGE_BDx_HSLICES, GLOB_SMALL_BDx_HSLICES
	integer :: GLOB_BDx_HSLICE, GLOB_BDy_HSLICE
	integer :: GLOB_LARGE_BDy_HSLICES, GLOB_SMALL_BDy_HSLICES
 
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
c    Define BDz (column) data types
c-------------------------------------------------------------------------------
      call MPI_TYPE_VECTOR(is, 1, ix*iy, MPI_REAL,
     &				GLOBAL_BDz_VCOL, ierr)
      call MPI_TYPE_COMMIT(GLOBAL_BDz_VCOL, ierr)
      call MPI_TYPE_VECTOR(is, 1, icloc, MPI_REAL,
     &				LOCAL_BDz_VCOL, ierr)
      call MPI_TYPE_COMMIT(LOCAL_BDz_VCOL, ierr)

      types(1) = GLOBAL_BDz_VCOL
      call MPI_TYPE_STRUCT(2,blength,tdispls,
     &		types, GLOB_BDz_VCOL, ierr)
      call MPI_TYPE_COMMIT(GLOB_BDz_VCOL, ierr)

      call MPI_TYPE_HVECTOR(icloc, 1, 
     &		NVworkers*sizeofreal, GLOBAL_BDz_VCOL,
     &		GLOBAL_LARGE_BDz_VCOLS, ierr)
      call MPI_TYPE_COMMIT(GLOBAL_LARGE_BDz_VCOLS, ierr)
      call MPI_TYPE_HVECTOR(icloc-1, 1, 
     &		NVworkers*sizeofreal, GLOBAL_BDz_VCOL,
     &		GLOBAL_SMALL_BDz_VCOLS, ierr)
      call MPI_TYPE_COMMIT(GLOBAL_SMALL_BDz_VCOLS, ierr)

      call MPI_TYPE_HVECTOR(icloc, 1, 
     &		sizeofreal, LOCAL_BDz_VCOL,
     &		LOCAL_LARGE_BDz_VCOLS, ierr)
      call MPI_TYPE_COMMIT(LOCAL_LARGE_BDz_VCOLS, ierr)
      call MPI_TYPE_HVECTOR(icloc-1, 1, 
     &		sizeofreal, LOCAL_BDz_VCOL,
     &		LOCAL_SMALL_BDz_VCOLS, ierr)
      call MPI_TYPE_COMMIT(LOCAL_SMALL_BDz_VCOLS, ierr)

      types(1) = GLOBAL_LARGE_BDz_VCOLS
      call MPI_TYPE_STRUCT(2,blength,tdispls,
     &		types, GLOB_LARGE_BDz_VCOLS, ierr)
      call MPI_TYPE_COMMIT(GLOB_LARGE_BDz_VCOLS, ierr)
      types(1) = GLOBAL_SMALL_BDz_VCOLS
      call MPI_TYPE_STRUCT(2,blength,tdispls,
     &		types, GLOB_SMALL_BDz_VCOLS, ierr)
      call MPI_TYPE_COMMIT(GLOB_SMALL_BDz_VCOLS, ierr)

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
c    Define BDx,BDy hslice data types (X and Y)
c-------------------------------------------------------------------------------
      call MPI_TYPE_VECTOR(is, ix, ix*iz, MPI_REAL, 
     &                      GLOBAL_BDy_HSLICE, ierr)
      call MPI_TYPE_COMMIT(GLOBAL_BDy_HSLICE, ierr)
      call MPI_TYPE_VECTOR(is, ix, ix*izloc, MPI_REAL, 
     &                      LOCAL_BDy_HSLICE, ierr)
      call MPI_TYPE_COMMIT(LOCAL_BDy_HSLICE, ierr)
      
      types(1) = GLOBAL_BDy_HSLICE
      call MPI_TYPE_STRUCT(2,blength,tdispls,
     &		types, GLOB_BDy_HSLICE, ierr)
      call MPI_TYPE_COMMIT(GLOB_BDy_HSLICE, ierr)

      call MPI_TYPE_VECTOR(is, iy, iy*iz, MPI_REAL, 
     &                      GLOBAL_BDx_HSLICE, ierr)
      call MPI_TYPE_COMMIT(GLOBAL_BDx_HSLICE, ierr)
      call MPI_TYPE_VECTOR(is, iy, iy*izloc, MPI_REAL, 
     &                      LOCAL_BDx_HSLICE, ierr)
      call MPI_TYPE_COMMIT(LOCAL_BDx_HSLICE, ierr)

      types(1) = GLOBAL_BDx_HSLICE
      call MPI_TYPE_STRUCT(2,blength,tdispls,
     &		types, GLOB_BDx_HSLICE, ierr)
      call MPI_TYPE_COMMIT(GLOB_BDx_HSLICE, ierr)

      call MPI_TYPE_HVECTOR(izloc, 1, 
     &		NHworkers*sizeofreal*ix, GLOBAL_BDy_HSLICE,
     &		GLOBAL_LARGE_BDy_HSLICES, ierr)
      call MPI_TYPE_COMMIT(GLOBAL_LARGE_BDy_HSLICES, ierr)
      call MPI_TYPE_HVECTOR(izloc-1, 1, 
     &		NHworkers*sizeofreal*ix, GLOBAL_BDy_HSLICE,
     &		GLOBAL_SMALL_BDy_HSLICES, ierr)
      call MPI_TYPE_COMMIT(GLOBAL_SMALL_BDy_HSLICES, ierr)

      call MPI_TYPE_HVECTOR(izloc, 1, 
     &		sizeofreal*ix, LOCAL_BDy_HSLICE,
     &		LOCAL_LARGE_BDy_HSLICES, ierr)
      call MPI_TYPE_COMMIT(LOCAL_LARGE_BDy_HSLICES, ierr)
      call MPI_TYPE_HVECTOR(izloc-1, 1, 
     &		sizeofreal*ix, LOCAL_BDy_HSLICE,
     &		LOCAL_SMALL_BDy_HSLICES, ierr)
      call MPI_TYPE_COMMIT(LOCAL_SMALL_BDy_HSLICES, ierr)

      call MPI_TYPE_HVECTOR(izloc, 1, 
     &		NHworkers*sizeofreal*iy, GLOBAL_BDx_HSLICE,
     &		GLOBAL_LARGE_BDx_HSLICES, ierr)
      call MPI_TYPE_COMMIT(GLOBAL_LARGE_BDx_HSLICES, ierr)
      call MPI_TYPE_HVECTOR(izloc-1, 1, 
     &		NHworkers*sizeofreal*iy, GLOBAL_BDx_HSLICE,
     &		GLOBAL_SMALL_BDx_HSLICES, ierr)
      call MPI_TYPE_COMMIT(GLOBAL_SMALL_BDx_HSLICES, ierr)

      call MPI_TYPE_HVECTOR(izloc, 1, 
     &		sizeofreal*iy, LOCAL_BDx_HSLICE,
     &		LOCAL_LARGE_BDx_HSLICES, ierr)
      call MPI_TYPE_COMMIT(LOCAL_LARGE_BDx_HSLICES, ierr)
      call MPI_TYPE_HVECTOR(izloc-1, 1, 
     &		sizeofreal*iy, LOCAL_BDx_HSLICE,
     &		LOCAL_SMALL_BDx_HSLICES, ierr)
      call MPI_TYPE_COMMIT(LOCAL_SMALL_BDx_HSLICES, ierr)

      types(1) = GLOBAL_LARGE_BDy_HSLICES
      call MPI_TYPE_STRUCT(2,blength,tdispls,
     &		types, GLOB_LARGE_BDy_HSLICES, ierr)
      call MPI_TYPE_COMMIT(GLOB_LARGE_BDy_HSLICES, ierr)
      types(1) = GLOBAL_SMALL_BDy_HSLICES
      call MPI_TYPE_STRUCT(2,blength,tdispls,
     &		types, GLOB_SMALL_BDy_HSLICES, ierr)
      call MPI_TYPE_COMMIT(GLOB_SMALL_BDy_HSLICES, ierr)
      types(1) = GLOBAL_LARGE_BDx_HSLICES
      call MPI_TYPE_STRUCT(2,blength,tdispls,
     &		types, GLOB_LARGE_BDx_HSLICES, ierr)
      call MPI_TYPE_COMMIT(GLOB_LARGE_BDx_HSLICES, ierr)
      types(1) = GLOBAL_SMALL_BDx_HSLICES
      call MPI_TYPE_STRUCT(2,blength,tdispls,
     &		types, GLOB_SMALL_BDx_HSLICES, ierr)
      call MPI_TYPE_COMMIT(GLOB_SMALL_BDx_HSLICES, ierr)

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
