      program mpi_communication_setup
      use XYParallelDataMap
      use XYCommDataTypes
      use XYCommunicationLibrary
      
      implicit none
      character(20) :: version_name
      real, pointer, dimension(:,:,:,:) :: sg1
      real, pointer, dimension(:,:) :: g2d
      real, pointer, dimension(:,:,:) :: gBDz
      real, pointer, dimension(:,:,:) :: g3d
      real, pointer, dimension(:,:,:) :: gBDx
      real, pointer, dimension(:,:,:) :: gBDy
      real, pointer, dimension(:,:,:,:) :: gBDx2
      real, pointer, dimension(:,:,:,:) :: gBDy2

      real, pointer, dimension(:,:,:,:) ::  sg2
      real, pointer, dimension(:,:,:,:) :: sgx
      real, pointer, dimension(:,:,:,:) :: sgy
      real, pointer, dimension(:,:) :: s2dx, s2dy
      real, pointer, dimension(:,:,:) :: sBDzx, sBDzy
      real, pointer, dimension(:,:,:) :: s3dx, s3dy
      real, pointer, dimension(:,:,:) :: sBDx, sBDy
      real, pointer, dimension(:,:,:,:) :: sBDx2, sBDy2
      integer ix, iy, iz, is
      integer  i, j, k, l, rc, ierr
      integer ixloc, iyloc
      integer,parameter :: Ntests=3
      real error
      integer :: decision, dcount(3)
      integer, pointer, dimension(:) :: mdecision

c-------------------------------------------------------------------------------
c   MPI-realated declarations
c-------------------------------------------------------------------------------
      integer status(MPI_STATUS_SIZE)
      double precision :: Tstart, Tend
      double precision :: testTimes(Ntests,3), best

c-------------------------------------------------------------------------------
c   Initialize MPI Bull
c-------------------------------------------------------------------------------
      call MPI_INIT( ierr )
      call MPI_COMM_RANK( MPI_COMM_WORLD, myid, ierr )
      call MPI_COMM_SIZE( MPI_COMM_WORLD, Nprocs, ierr )

c     Master not defined yet - therefore myid = 0
      if (MyId .eq. 0) then
        print*, 'Please input the size of the 4D',
     &	' concentration array in the form (x,y,z,Ns) :'
        read*, ix, iy, iz, is
        print*, 'Initializing 4D array (',ix,iy,iz,is,
     &	') and other arrays'
      end if
      call MPI_BARRIER(MPI_COMM_WORLD, ierr)
      call MPI_Bcast(ix, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
      call MPI_Bcast(iy, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
      call MPI_Bcast(iz, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
      call MPI_Bcast(is, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

      

c      print *, 'Process ', myid, ' of ', Nprocs, ' is alive'

c-------------------------------------------------------------------------------
c    The data mapping scheme
c-------------------------------------------------
      call CreateMap(ix, iy, ixloc, iyloc)
c      if (XWorker .or. YWorker) 
c     &	print*,'P[',myid,'] has ixloc=',ixloc,'iyloc=',iyloc
c      print*, 'P[',myid,'] has xslices=',owned_xslices(myid,:)
c      print*, 'P[',myid,'] has yslices=',owned_yslices(myid,:)

c-------------------------------------------------------------------------------
c    Create Special Communication Patterns
c-------------------------------------------------
      call CreateCommDataTypes(ix,iy,iz,ixloc,iyloc,is)
       
c-------------------------------------------------------------------------------
c    Allocate the arrays
c-------------------------------------------------------------------------------
      allocate(sg1(ix,iy,iz,is),STAT=ierr)
      allocate(g2d(ix,iy),STAT=ierr)
      allocate(gBDz(ix,iy,is),STAT=ierr)
      allocate(g3d(ix,iy,iz),STAT=ierr)
      allocate(gBDx(iy,iz,is),STAT=ierr)
      allocate(gBDy(ix,iz,is),STAT=ierr)
      allocate(gBDx2(iy,iz,2,is),STAT=ierr)
      allocate(gBDy2(ix,iz,2,is),STAT=ierr)

      if (myid .eq. 0) then
        allocate(sg2(ix,iy,iz,is),STAT=ierr)

        allocate(sgx(1,1,1,1),STAT=ierr)
        allocate(sgy(1,1,1,1),STAT=ierr)
        allocate(s2dx(1,1),STAT=ierr)
        allocate(s2dy(1,1),STAT=ierr)
        allocate(sBDzx(1,1,1),STAT=ierr)
        allocate(sBDzy(1,1,1),STAT=ierr)
        allocate(s3dx(1,1,1),STAT=ierr)
        allocate(s3dy(1,1,1),STAT=ierr)
        allocate(sBDx(1,1,1),STAT=ierr)
        allocate(sBDy(1,1,1),STAT=ierr)
        allocate(sBDx2(1,1,2,1),STAT=ierr)
        allocate(sBDy2(1,1,2,1),STAT=ierr)
      else	
        allocate(sg2(1,1,1,1),STAT=ierr)
	
        allocate(sgx(ix,iyloc,iz,is),STAT=ierr)
        allocate(sgy(ixloc,iy,iz,is),STAT=ierr)

        allocate(s2dx(ix,iyloc),STAT=ierr)
        allocate(s2dy(ixloc,iy),STAT=ierr)
        allocate(sBDzx(ix,iyloc,is),STAT=ierr)
        allocate(sBDzy(ixloc,iy,is),STAT=ierr)
        allocate(s3dx(ix,iyloc,iz),STAT=ierr)
        allocate(s3dy(ixloc,iy,iz),STAT=ierr)
        allocate(sBDx(iyloc,iz,is),STAT=ierr)
        allocate(sBDy(ixloc,iz,is),STAT=ierr)
        allocate(sBDx2(iyloc,iz,2,is),STAT=ierr)
        allocate(sBDy2(ixloc,iz,2,is),STAT=ierr)
      endif

c-------------------------------------------------------------------------------
c    Master Initializes the concentration array - all initialize it here...
c-------------------------------------------------------------------------------
        do i=1,ix
	  do j=1,iy
	    do k=1,iz
	      do l=1,is
	        sg1(i,j,k,l) = l-1 + is*(k-1) + is*iz*(j-1) + is*iz*iy*(i-1)
     &			+ 1.0/3.0
		g2d(i,j) = j-1 + iy*(i-1)
     &			+ 1.0/3.0
		gBDz(i,j,l) = l-1 + is*(j-1) + is*iy*(i-1)
     &			+ 1.0/3.0
		g3d(i,j,k) = k-1 + iz*(j-1) + iz*iy*(i-1)
     &			+ 1.0/3.0
		gBDx(j,k,l) = l-1 + is*(k-1) + is*iz*(j-1)
     &			+ 1.0/3.0
		gBDx2(j,k,1,l) = l-1 + is*(k-1) + is*iz*(j-1)
     &			+ 1.0/3.0
		gBDx2(j,k,2,l) = is*iz*iy + l-1 + is*(k-1) + is*iz*(j-1)
     &			+ 1.0/3.0
		gBDy(i,k,l) = l-1 + is*(k-1) + is*iz*(i-1)
     &			+ 1.0/3.0
		gBDy2(i,k,1,l) = l-1 + is*(k-1) + is*iz*(i-1)
     &			+ 1.0/3.0
		gBDy2(i,k,2,l) = is*iz*ix + l-1 + is*(k-1) + is*iz*(i-1)
     &			+ 1.0/3.0
              end do		
            end do		
          end do		
        end do	

c-------------------------------------------------------------------------------
c     Testing communication methods starts here     
c-------------------------------------------------------------------------------

      if (Master) then
        allocate(mdecision(Nprocs), STAT=ierr)
      else
        allocate(mdecision(1), STAT=ierr)
      end if
c-------------------------------------------------------------------------------
c     Test 4D x-distribution     
c-------------------------------------------------------------------------------

      if (Master) print*, 'Standard 4D x-distribution test starts'
      do i=1, Ntests
        call MPI_BARRIER(MPI_COMM_WORLD, ierr)
	Tstart = MPI_WTIME()
	call distrib_x_4D_v1(ix,iy,iz,iyloc,is, sg1,sgx)
	if (XWorker) Tend = MPI_WTIME()
	call MPI_BARRIER(MPI_COMM_WORLD, ierr)
	if (Master) Tend = MPI_WTIME()
	testTimes(i,1) = Tend - Tstart
	if (XWorker) then
	  error = maxval(abs(
     &		sg1(:,owned_xslices(myid,1:no_of_xslices(myid)),:,:)
     &   	-sgx(:,1:no_of_xslices(myid),:,:)))
          if (error .gt. 0.0) then
	    print*, 'Proc',myid,'After standard 4D x-distribution',
     &		' max error = ', error, ' check implementation '
          end if
        end if
	sgx(:,:,:,:) = 0
      end do
      

      if (Master) print*, 'Single msg 4D x-distribution test starts'
      do i=1, Ntests
        call MPI_BARRIER(MPI_COMM_WORLD, ierr)
	Tstart = MPI_WTIME()
	call distrib_x_4D_v2(ix,iy,iz,iyloc,is, sg1,sgx)
	if (XWorker) Tend = MPI_WTIME()
	call MPI_BARRIER(MPI_COMM_WORLD, ierr)
	if (Master) Tend = MPI_WTIME()
	testTimes(i,2) = Tend - Tstart
	if (XWorker) then
	  error = maxval(abs(
     &		sg1(:,owned_xslices(myid,1:no_of_xslices(myid)),:,:)
     &   	-sgx(:,1:no_of_xslices(myid),:,:)))
          if (error .gt. 0.0) then
	    print*, 'Proc',myid,'After single msg 4D x-distribution',
     &		' max error = ', error, ' check implementation '
          end if
        end if
	sgx(:,:,:,:) = 0
      end do

      if (Master) print*, 'Collective 4D x-distribution test starts'
      do i=1, Ntests
        call MPI_BARRIER(MPI_COMM_WORLD, ierr)
	Tstart = MPI_WTIME()
	call distrib_x_4D_v3(ix,iy,iz,iyloc,is, sg1,sgx)
	if (XWorker) Tend = MPI_WTIME()
	call MPI_BARRIER(MPI_COMM_WORLD, ierr)
	if (Master) Tend = MPI_WTIME()
	testTimes(i,3) = Tend - Tstart
	if (XWorker) then
	  error = maxval(abs(
     &		sg1(:,owned_xslices(myid,1:no_of_xslices(myid)),:,:)
     &   	-sgx(:,1:no_of_xslices(myid),:,:)))
          if (error .gt. 0.0) then
	    print*, 'Proc',myid,'After collective 4D x-distribution',
     &		' max error = ', error, ' check implementation '
          end if
        end if
	sgx(:,:,:,:) = 0
      end do
      
      if (Master) then
        print*, 'Global standard 4D x-distribution time:',
     &		minval(testTimes(1:Ntests,1))
        print*, 'Global single msg 4D x-distribution time:',
     &		minval(testTimes(1:Ntests,2))
        print*, 'Global collective 4D x-distribution time:',
     &		minval(testTimes(1:Ntests,3))
      end if

      best = minval(testTimes(:,:))
      if (best .eq. minval(testTimes(1:Ntests,1))) then
	decision = 1
      else if (best .eq. minval(testTimes(1:Ntests,2))) then
        decision = 2
      else
	decision = 3
      end if
      
      call MPI_GATHER(decision, 1, MPI_INTEGER, mdecision, 1, 
     &		MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
      if (Master) then
        dcount(:) = 0
        do i=1,NXworkers+1
          if (mdecision(i) .eq. 1) then
	    dcount(1) = dcount(1) + 1
	  else if (mdecision(i) .eq. 2) then
	    dcount(2) = dcount(2) + 1
          else 
	    dcount(3) = dcount(3) + 1
	  end if
        end do

c	write to file XY_VERSIONS
        OPEN (UNIT=10, FILE='XY_versions')
	
        best = maxval(dcount(:))
	if (best .eq. dcount(1)) then
	  print*, 'Decision: standard 4D x-distribution.'
	  WRITE (10,*) 'distr_x_4D= ',1
	else if (best .eq. dcount(2)) then
	  print*, 'Decision: single msg 4D x-distribution.'
	  WRITE (10,*) 'distr_x_4D= ',2
	else 
	  print*, 'Decision: collective 4D x-distribution.'
	  WRITE (10,*) 'distr_x_4D= ',3
	end if
      end if
      
c-------------------------------------------------------------------------------
c     Test 4D y-distribution     
c-------------------------------------------------------------------------------
      if (Master) print*, 'Standard 4D y-distribution test starts'
      do i=1, Ntests
        call MPI_BARRIER(MPI_COMM_WORLD, ierr)
	Tstart = MPI_WTIME()
	call distrib_y_4D_v1(ix,iy,iz,ixloc,is, sg1,sgy)
	if (YWorker) Tend = MPI_WTIME()
	call MPI_BARRIER(MPI_COMM_WORLD, ierr)
	if (Master) Tend = MPI_WTIME()
	testTimes(i,1) = Tend - Tstart
	if (YWorker) then
	  error = maxval(abs(
     &		sg1(owned_yslices(myid,1:no_of_yslices(myid)),:,:,:)
     &   	-sgy(1:no_of_yslices(myid),:,:,:)))
          if (error .gt. 0.0) then
	    print*, 'Proc',myid,'After standard 4D y-distribution',
     &		' max error = ', error, ' check implementation '
          end if
        end if
	sgy(:,:,:,:) = 0
      end do
      

      if (Master) print*, 'Single msg 4D y-distribution test starts'
      do i=1, Ntests
        call MPI_BARRIER(MPI_COMM_WORLD, ierr)
	Tstart = MPI_WTIME()
	call distrib_y_4D_v2(ix,iy,iz,ixloc,is, sg1,sgy)
	if (YWorker) Tend = MPI_WTIME()
	call MPI_BARRIER(MPI_COMM_WORLD, ierr)
	if (Master) Tend = MPI_WTIME()
	testTimes(i,2) = Tend - Tstart
	if (YWorker) then
	  error = maxval(abs(
     &		sg1(owned_yslices(myid,1:no_of_yslices(myid)),:,:,:)
     &   	-sgy(1:no_of_yslices(myid),:,:,:)))
          if (error .gt. 0.0) then
	    print*, 'Proc',myid,'After single msg 4D y-distribution',
     &		' max error = ', error, ' check implementation '
          end if
        end if
	sgy(:,:,:,:) = 0
      end do

      if (Master) print*, 'Collective 4D y-distribution test starts'
      do i=1, Ntests
        call MPI_BARRIER(MPI_COMM_WORLD, ierr)
	Tstart = MPI_WTIME()
	call distrib_y_4D_v3(ix,iy,iz,ixloc,is, sg1,sgy)
	if (YWorker) Tend = MPI_WTIME()
	call MPI_BARRIER(MPI_COMM_WORLD, ierr)
	if (Master) Tend = MPI_WTIME()
	testTimes(i,3) = Tend - Tstart
	if (YWorker) then
	  error = maxval(abs(
     &		sg1(owned_yslices(myid,1:no_of_yslices(myid)),:,:,:)
     &   	-sgy(1:no_of_yslices(myid),:,:,:)))
          if (error .gt. 0.0) then
	    print*, 'Proc',myid,'After collective 4D y-distribution',
     &		' max error = ', error, ' check implementation '
          end if
        end if
	sgy(:,:,:,:) = 0
      end do
      
      if (Master) then
        print*, 'Global standard 4D y-distribution time:',
     &		minval(testTimes(1:Ntests,1))
        print*, 'Global single msg 4D y-distribution time:',
     &		minval(testTimes(1:Ntests,2))
        print*, 'Global collective 4D y-distribution time:',
     &		minval(testTimes(1:Ntests,3))
      end if

      best = minval(testTimes(:,:))
      if (best .eq. minval(testTimes(1:Ntests,1))) then
	decision = 1
      else if (best .eq. minval(testTimes(1:Ntests,2))) then
        decision = 2
      else
	decision = 3
      end if
      
      call MPI_GATHER(decision, 1, MPI_INTEGER, mdecision, 1, 
     &		MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
      if (Master) then
        dcount(:) = 0
        do i=1,NYworkers+1
          if (mdecision(i) .eq. 1) then
	    dcount(1) = dcount(1) + 1
	  else if (mdecision(i) .eq. 2) then
	    dcount(2) = dcount(2) + 1
          else 
	    dcount(3) = dcount(3) + 1
	  end if
        end do
      
        best = maxval(dcount(:))
	if (best .eq. dcount(1)) then
	  print*, 'Decision: standard 4D y-distribution.'
	  WRITE (10,*) 'distr_y_4D= ',1
	else if (best .eq. dcount(2)) then
	  print*, 'Decision: single msg 4D y-distribution.'
	  WRITE (10,*) 'distr_y_4D= ',2
	else 
	  print*, 'Decision: collective 4D y-distribution.'
	  WRITE (10,*) 'distr_y_4D= ',3
	end if
      end if

c-------------------------------------------------------------------------------
c     Test 4D x-gather     
c-------------------------------------------------------------------------------

      if (Master) print*, 'Standard 4D x-gather test starts'
c     in order to check the error we need to call distribution first
      call distrib_x_4D_v1(ix,iy,iz,iyloc,is, sg1,sgx)
      
      do i=1, Ntests
        call MPI_BARRIER(MPI_COMM_WORLD, ierr)
	Tstart = MPI_WTIME()
	call gather_x_4D_v1(ix,iy,iz,iyloc,is, sg2,sgx)
	Tend = MPI_WTIME()
	testTimes(i,1) = Tend - Tstart
	if (Master) then
	  error = maxval(abs(sg2-sg1))
          if (error .gt. 0.0) then
	    print*, 'Proc',myid,'After standard 4D x-gather',
     &		' max error = ', error, ' check implementation '
          end if
        end if
	sg2(:,:,:,:) = 0
      end do
      

      if (Master) print*, 'Single msg 4D x-gather test starts'
      do i=1, Ntests
        call MPI_BARRIER(MPI_COMM_WORLD, ierr)
	Tstart = MPI_WTIME()
	call gather_x_4D_v2(ix,iy,iz,iyloc,is, sg2,sgx)
	Tend = MPI_WTIME()
	testTimes(i,2) = Tend - Tstart
	if (Master) then
	  error = maxval(abs(sg2-sg1))
          if (error .gt. 0.0) then
	    print*, 'Proc',myid,'After single msg 4D x-gather',
     &		' max error = ', error, ' check implementation '
          end if
        end if
	sg2(:,:,:,:) = 0
      end do

      if (Master) print*, 'Collective 4D x-gather test starts'
      do i=1, Ntests
        call MPI_BARRIER(MPI_COMM_WORLD, ierr)
	Tstart = MPI_WTIME()
	call gather_x_4D_v3(ix,iy,iz,iyloc,is, sg2,sgx)
	Tend = MPI_WTIME()
	testTimes(i,3) = Tend - Tstart
	if (Master) then
	  error = maxval(abs(sg2-sg1))
          if (error .gt. 0.0) then
	    print*, 'Proc',myid,'After collective 4D x-gather',
     &		'max error = ', error, ' check implementation '
          end if
        end if
	sg2(:,:,:,:) = 0
      end do
      
      if (Master) then
        print*, 'Global standard 4D x-gather. time:',
     &		minval(testTimes(1:Ntests,1))
        print*, 'Global single msg 4D x-gather time:',
     &		minval(testTimes(1:Ntests,2))
        print*, 'Global collective 4D x-gather time:',
     &		minval(testTimes(1:Ntests,3))
      end if

      best = minval(testTimes(:,:))
      if (best .eq. minval(testTimes(1:Ntests,1))) then
	decision = 1
      else if (best .eq. minval(testTimes(1:Ntests,2))) then
        decision = 2
      else
	decision = 3
      end if
      
      call MPI_GATHER(decision, 1, MPI_INTEGER, mdecision, 1, 
     &		MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
      if (Master) then
        dcount(:) = 0
        do i=1,NXworkers+1
          if (mdecision(i) .eq. 1) then
	    dcount(1) = dcount(1) + 1
  	  else if (mdecision(i) .eq. 2) then
	    dcount(2) = dcount(2) + 1
          else 
	    dcount(3) = dcount(3) + 1
	  end if
        end do
      
        best = maxval(dcount(:))
	if (best .eq. dcount(1)) then
	  print*, 'Decision: standard 4D x-gather.'
	  WRITE (10,*) 'gather_x_4D= ',1
	else if (best .eq. dcount(2)) then
	  print*, 'Decision: single msg 4D x-gather.'
	  WRITE (10,*) 'gather_x_4D= ',2
	else 
	  print*, 'Decision: collective 4D x-gather.'
	  WRITE (10,*) 'gather_x_4D= ',3
	end if
      end if

c-------------------------------------------------------------------------------
c     Test 4D y-gather     
c-------------------------------------------------------------------------------
      if (Master) print*, 'Standard 4D y-gather test starts'
c     in order to check the error we need to call distribution first
      call distrib_y_4D_v1(ix,iy,iz,ixloc,is, sg1,sgy)

      do i=1, Ntests
        call MPI_BARRIER(MPI_COMM_WORLD, ierr)
	Tstart = MPI_WTIME()
	call gather_y_4D_v1(ix,iy,iz,ixloc,is, sg2,sgy)
	Tend = MPI_WTIME()
	testTimes(i,1) = Tend - Tstart
	if (Master) then
	  error = maxval(abs(sg2-sg1))
          if (error .gt. 0.0) then
	    print*, 'Proc',myid,'After standard 4D y-gather',
     &		' max error = ', error, ' check implementation '
          end if
        end if
	sg2(:,:,:,:) = 0
      end do
      

      if (Master) print*, 'Single msg 4D y-gather test starts'
      do i=1, Ntests
        call MPI_BARRIER(MPI_COMM_WORLD, ierr)
	Tstart = MPI_WTIME()
	call gather_y_4D_v2(ix,iy,iz,ixloc,is, sg2,sgy)
	Tend = MPI_WTIME()
	testTimes(i,2) = Tend - Tstart
	if (Master) then
	  error = maxval(abs(sg2-sg1))
          if (error .gt. 0.0) then
	    print*, 'Proc',myid,'After single msg 4D y-gather',
     &		' max error = ', error, ' check implementation '
          end if
        end if
	sg2(:,:,:,:) = 0
      end do

      if (Master) print*, 'Collective 4D y-gather test starts'
      do i=1, Ntests
        call MPI_BARRIER(MPI_COMM_WORLD, ierr)
	Tstart = MPI_WTIME()
	call gather_y_4D_v3(ix,iy,iz,ixloc,is, sg2,sgy)
	Tend = MPI_WTIME()
	testTimes(i,3) = Tend - Tstart
	if (Master) then
	  error = maxval(abs(sg2-sg1))
          if (error .gt. 0.0) then
	    print*, 'Proc',myid,'After collective 4D y-gather',
     &		' max error = ', error, ' check implementation '
          end if
        end if
	sg2(:,:,:,:) = 0
      end do
      
      if (Master) then
        print*, 'Global standard 4D y-gather time:',
     &		minval(testTimes(1:Ntests,1))
        print*, 'Global single msg 4D y-gather. time:',
     &		minval(testTimes(1:Ntests,2))
        print*, 'Global collective 4D y-gather time:',
     &		minval(testTimes(1:Ntests,3))
      end if

      best = minval(testTimes(:,:))
      if (best .eq. minval(testTimes(1:Ntests,1))) then
	decision = 1
      else if (best .eq. minval(testTimes(1:Ntests,2))) then
        decision = 2
      else
	decision = 3
      end if
      
      call MPI_GATHER(decision, 1, MPI_INTEGER, mdecision, 1, 
     &		MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
      if (Master) then
        dcount(:) = 0
        do i=1,NYworkers+1
          if (mdecision(i) .eq. 1) then
	    dcount(1) = dcount(1) + 1
	  else if (mdecision(i) .eq. 2) then
	    dcount(2) = dcount(2) + 1
          else 
	    dcount(3) = dcount(3) + 1
 	  end if
        end do
      
        best = maxval(dcount(:))
	if (best .eq. dcount(1)) then
	  print*, 'Decision: standard 4D y-gather.'
	  WRITE (10,*) 'gather_y_4D= ',1
	else if (best .eq. dcount(2)) then
	  print*, 'Decision: single msg 4D y-gather.'
	  WRITE (10,*) 'gather_y_4D= ',2
	else 
	  print*, 'Decision: collective 4D y-gather.'
	  WRITE (10,*) 'gather_y_4D= ',3
	end if
      end if

	
c-------------------------------------------------------------------------------
c     Test x2y shuffle    
c-------------------------------------------------------------------------------
      if (Master) print*, 'Standard 4D x2y shuffle test starts'
c     in order to check the error we need to call distribution first
      call distrib_x_4D_v1(ix,iy,iz,iyloc,is, sg1,sgx)
      sgy(:,:,:,:) = 0

      do i=1, Ntests
        call MPI_BARRIER(MPI_COMM_WORLD, ierr)
	Tstart = MPI_WTIME()
        call shuffle_x2y_4D_v1(ix,iy,iz,ixloc,iyloc,is,sgx,sgy)
	Tend = MPI_WTIME()
        call MPI_BARRIER(MPI_COMM_WORLD, ierr)
	if (Master) Tend = MPI_WTIME()
	testTimes(i,1) = Tend - Tstart
	if (YWorker) then
	  error = maxval(abs(
     &		sg1(owned_yslices(myid,1:no_of_yslices(myid)),:,:,:)
     &   -sgy(1:no_of_yslices(myid),:,:,:)))
          if (error .gt. 0.0) then
	    print*, 'Proc',myid,'After standard 4D x2y shuffle',
     &		' max error = ', error, ' check implementation '
          end if
        end if
	sgy(:,:,:,:) = 0
      end do
      

      if (Master) print*, 'Single msg 4D x2y shuffle test starts'
      do i=1, Ntests
        call MPI_BARRIER(MPI_COMM_WORLD, ierr)
	Tstart = MPI_WTIME()
        call shuffle_x2y_4D_v2(ix,iy,iz,ixloc,iyloc,is,sgx,sgy)
	Tend = MPI_WTIME()
        call MPI_BARRIER(MPI_COMM_WORLD, ierr)
	if (Master) Tend = MPI_WTIME()
	testTimes(i,2) = Tend - Tstart
	if (YWorker) then
	  error = maxval(abs(
     &		sg1(owned_yslices(myid,1:no_of_yslices(myid)),:,:,:)
     &   -sgy(1:no_of_yslices(myid),:,:,:)))
          if (error .gt. 0.0) then
	    print*, 'Proc',myid,'After single msg 4D x2y shuffle',
     &		' max error = ', error, ' check implementation '
          end if
        end if
	sgy(:,:,:,:) = 0
      end do

      if (Master) print*, 'Collective 4D x2y shuffle test starts'
      do i=1, Ntests
        call MPI_BARRIER(MPI_COMM_WORLD, ierr)
	Tstart = MPI_WTIME()
        call shuffle_x2y_4D_v3(ix,iy,iz,ixloc,iyloc,is,sgx,sgy)
	Tend = MPI_WTIME()
        call MPI_BARRIER(MPI_COMM_WORLD, ierr)
	if (Master) Tend = MPI_WTIME()
	testTimes(i,3) = Tend - Tstart
	if (YWorker) then
	  error = maxval(abs(
     &		sg1(owned_yslices(myid,1:no_of_yslices(myid)),:,:,:)
     &   -sgy(1:no_of_yslices(myid),:,:,:)))
          if (error .gt. 0.0) then
	    print*, 'Proc',myid,'After collective 4D x2y shuffle',
     &		' max error = ', error, ' check implementation '
          end if
        end if
	sgy(:,:,:,:) = 0
      end do
      
	      
      if (Master) then
        print*, 'Global standard 4D x2y shuffle time:',
     &		minval(testTimes(1:Ntests,1))
        print*, 'Global single msg 4D x2y shuffle time:',
     &		minval(testTimes(1:Ntests,2))
        print*, 'Global collective 4D x2y shuffle time:',
     &		minval(testTimes(1:Ntests,3))
      end if
      if (MyId .eq. 1) then
        print*, 'Proc 1 standard 4D x2y shuffle time:',
     &		minval(testTimes(1:Ntests,1))
        print*, 'Proc 1 single msg 4D x2y shuffle time:',
     &		minval(testTimes(1:Ntests,2))
        print*, 'Proc 1 collective 4D x2y shuffle time:',
     &		minval(testTimes(1:Ntests,3))
      end if

      best = minval(testTimes(:,:))
      if (best .eq. minval(testTimes(1:Ntests,1))) then
	decision = 1
      else if (best .eq. minval(testTimes(1:Ntests,2))) then
        decision = 2
      else
	decision = 3
      end if
      
      call MPI_GATHER(decision, 1, MPI_INTEGER, mdecision, 1, 
     &		MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
      if (Master) then
        dcount(:) = 0
        do i=2,max(NXworkers,NYworkers)+1
          if (mdecision(i) .eq. 1) then
	    dcount(1) = dcount(1) + 1
	  else if (mdecision(i) .eq. 2) then
	    dcount(2) = dcount(2) + 1
          else 
  	    dcount(3) = dcount(3) + 1
	  end if
        end do
      
        best = maxval(dcount(:))
	if (best .eq. dcount(1)) then
	  print*, 'Decision: standard 4D x2y shuffle.'
	  WRITE (10,*) 'shuffle_x2y_4D= ',1
	else if (best .eq. dcount(2)) then
	  print*, 'Decision: single msg 4D x2y shuffle.'
	  WRITE (10,*) 'shuffle_x2y_4D= ',2
	else 
	  print*, 'Decision: collective 4D x2y shuffle.'
	  WRITE (10,*) 'shuffle_x2y_4D= ',3
	end if
      end if

c-------------------------------------------------------------------------------
c     Test y2x shuffle    
c-------------------------------------------------------------------------------
      if (Master) print*, 'Standard 4D y2x shuffle test starts'
c     in order to check the error we need to call distribution first
      call distrib_y_4D_v1(ix,iy,iz,ixloc,is, sg1,sgy)
      sgx(:,:,:,:) = 0

      do i=1, Ntests
        call MPI_BARRIER(MPI_COMM_WORLD, ierr)
	Tstart = MPI_WTIME()
        call shuffle_y2x_4D_v1(ix,iy,iz,ixloc,iyloc,is,sgx,sgy)
	Tend = MPI_WTIME()
        call MPI_BARRIER(MPI_COMM_WORLD, ierr)
	if (Master) Tend = MPI_WTIME()
	testTimes(i,1) = Tend - Tstart
	if (XWorker) then
	  error = maxval(abs(
     &		sg1(:,owned_xslices(myid,1:no_of_xslices(myid)),:,:)
     &   -sgx(:,1:no_of_xslices(myid),:,:)))
          if (error .gt. 0.0) then
	    print*, 'Proc',myid,'After standard 4D y2x shuffle',
     &		' max error = ', error, ' check implementation '
          end if
        end if
	sgx(:,:,:,:) = 0
      end do
      

      if (Master) print*, 'Single msg 4D y2x shuffle test starts'
      do i=1, Ntests
        call MPI_BARRIER(MPI_COMM_WORLD, ierr)
	Tstart = MPI_WTIME()
        call shuffle_y2x_4D_v2(ix,iy,iz,ixloc,iyloc,is,sgx,sgy)
	Tend = MPI_WTIME()
        call MPI_BARRIER(MPI_COMM_WORLD, ierr)
	if (Master) Tend = MPI_WTIME()
	testTimes(i,2) = Tend - Tstart
	if (XWorker) then
	  error = maxval(abs(
     &		sg1(:,owned_xslices(myid,1:no_of_xslices(myid)),:,:)
     &   -sgx(:,1:no_of_xslices(myid),:,:)))
          if (error .gt. 0.0) then
	    print*, 'Proc',myid,'After single msg 4D y2x shuffle',
     &		' max error = ', error, ' check implementation '
          end if
        end if
	sgx(:,:,:,:) = 0
      end do

      if (Master) print*, 'Collective 4D y2x shuffle test starts'
      do i=1, Ntests
        call MPI_BARRIER(MPI_COMM_WORLD, ierr)
	Tstart = MPI_WTIME()
        call shuffle_y2x_4D_v3(ix,iy,iz,ixloc,iyloc,is,sgx,sgy)
	Tend = MPI_WTIME()
        call MPI_BARRIER(MPI_COMM_WORLD, ierr)
	if (Master) Tend = MPI_WTIME()
	testTimes(i,3) = Tend - Tstart
	if (XWorker) then
	  error = maxval(abs(
     &		sg1(:,owned_xslices(myid,1:no_of_xslices(myid)),:,:)
     &   -sgx(:,1:no_of_xslices(myid),:,:)))
          if (error .gt. 0.0) then
	    print*, 'Proc',myid,'After collective 4D y2x shuffle',
     &		' max error = ', error, ' check implementation '
          end if
        end if
	sgx(:,:,:,:) = 0
      end do
      
	      
      if (Master) then
        print*, 'Global standard 4D y2x shuffle time:',
     &		minval(testTimes(1:Ntests,1))
        print*, 'Global single msg 4D y2x shuffle time:',
     &		minval(testTimes(1:Ntests,2))
        print*, 'Global collective 4D y2x shuffle time:',
     &		minval(testTimes(1:Ntests,3))
      end if
      if (MyId .eq. 1) then
        print*, 'Proc 1 standard 4D y2x shuffle time:',
     &		minval(testTimes(1:Ntests,1))
        print*, 'Proc 1 single msg 4D y2x shuffle time:',
     &		minval(testTimes(1:Ntests,2))
        print*, 'Proc 1 collective 4D y2x shuffle time:',
     &		minval(testTimes(1:Ntests,3))
      end if

      best = minval(testTimes(:,:))
      if (best .eq. minval(testTimes(1:Ntests,1))) then
	decision = 1
      else if (best .eq. minval(testTimes(1:Ntests,2))) then
        decision = 2
      else
	decision = 3
      end if
      
      call MPI_GATHER(decision, 1, MPI_INTEGER, mdecision, 1, 
     &		MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
      if (Master) then
        dcount(:) = 0
        do i=2,max(NXworkers,NYworkers)+1
          if (mdecision(i) .eq. 1) then
	    dcount(1) = dcount(1) + 1
	  else if (mdecision(i) .eq. 2) then
	    dcount(2) = dcount(2) + 1
          else 
	    dcount(3) = dcount(3) + 1
  	  end if
        end do
      
        best = maxval(dcount(:))
	if (best .eq. dcount(1)) then
	  print*, 'Decision: standard 4D y2x shuffle.'
	  WRITE (10,*) 'shuffle_y2x_4D= ',1
	else if (best .eq. dcount(2)) then
	  print*, 'Decision: single msg 4D y2x shuffle.'
	  WRITE (10,*) 'shuffle_y2x_4D= ',2
	else 
	  print*, 'Decision: collective 4D y2x shuffle.'
	  WRITE (10,*) 'shuffle_y2x_4D= ',3
	end if
      end if

c-------------------------------------------------------------------------------
c     Test 2D x-distribution     
c-------------------------------------------------------------------------------
	
	
      if (Master) print*, 'Standard 2D x-distribution test starts'
      do i=1, Ntests
        call MPI_BARRIER(MPI_COMM_WORLD, ierr)
	Tstart = MPI_WTIME()
	call distrib_x_2D_v1(ix,iy,iyloc, g2d,s2dx)
	if (XWorker) Tend = MPI_WTIME()
	call MPI_BARRIER(MPI_COMM_WORLD, ierr)
	if (Master) Tend = MPI_WTIME()
	testTimes(i,1) = Tend - Tstart
	if (XWorker) then
	  error = maxval(abs(
     &		g2d(:,owned_xslices(myid,1:no_of_xslices(myid)))
     &   	-s2dx(:,1:no_of_xslices(myid))))
          if (error .gt. 0.0) then
	    print*, 'Proc',myid,'After standard 2D x-distribution',
     &		' max error = ', error, ' check implementation '
          end if
        end if
	s2dx(:,:) = 0
      end do
      

      if (Master) print*, 'Single msg 2D x-distribution test starts'
      do i=1, Ntests
        call MPI_BARRIER(MPI_COMM_WORLD, ierr)
	Tstart = MPI_WTIME()
	call distrib_x_2D_v2(ix,iy,iyloc, g2d,s2dx)
	if (XWorker) Tend = MPI_WTIME()
	call MPI_BARRIER(MPI_COMM_WORLD, ierr)
	if (Master) Tend = MPI_WTIME()
	testTimes(i,2) = Tend - Tstart
	if (XWorker) then
	  error = maxval(abs(
     &		g2d(:,owned_xslices(myid,1:no_of_xslices(myid)))
     &   	-s2dx(:,1:no_of_xslices(myid))))
          if (error .gt. 0.0) then
	    print*, 'Proc',myid,'After single msg 2D x-distribution',
     &		' max error = ', error, ' check implementation '
          end if
        end if
	s2dx(:,:) = 0
      end do

      if (Master) print*, 'Collective 2D x-distribution test starts'
      do i=1, Ntests
        call MPI_BARRIER(MPI_COMM_WORLD, ierr)
	Tstart = MPI_WTIME()
	call distrib_x_2D_v3(ix,iy,iyloc, g2d,s2dx)
	if (XWorker) Tend = MPI_WTIME()
	call MPI_BARRIER(MPI_COMM_WORLD, ierr)
	if (Master) Tend = MPI_WTIME()
	testTimes(i,3) = Tend - Tstart
	if (XWorker) then
	  error = maxval(abs(
     &		g2d(:,owned_xslices(myid,1:no_of_xslices(myid)))
     &   	-s2dx(:,1:no_of_xslices(myid))))
          if (error .gt. 0.0) then
	    print*, 'Proc',myid,'After collective 2D x-distribution',
     &		' max error = ', error, ' check implementation '
          end if
        end if
	s2dx(:,:) = 0
      end do
      if (Master) then
        print*, 'Global standard 2D x-distribution time:',
     &		minval(testTimes(1:Ntests,1))
        print*, 'Global single msg 2D x-distribution time:',
     &		minval(testTimes(1:Ntests,2))
        print*, 'Global collective 2D x-distribution time:',
     &		minval(testTimes(1:Ntests,3))
      end if

      best = minval(testTimes(:,:))
      if (best .eq. minval(testTimes(1:Ntests,1))) then
	decision = 1
      else if (best .eq. minval(testTimes(1:Ntests,2))) then
        decision = 2
      else
	decision = 3
      end if
      
      call MPI_GATHER(decision, 1, MPI_INTEGER, mdecision, 1, 
     &		MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
      if (Master) then
        dcount(:) = 0
        do i=1,NXworkers+1
          if (mdecision(i) .eq. 1) then
	    dcount(1) = dcount(1) + 1
	  else if (mdecision(i) .eq. 2) then
	    dcount(2) = dcount(2) + 1
          else 
	    dcount(3) = dcount(3) + 1
	  end if
        end do
        best = maxval(dcount(1:3))
	if (best .eq. dcount(1)) then
	  print*, 'Decision: standard 2D x-distribution.'
	  WRITE (10,*) 'distr_x_2D= ',1
	else if (best .eq. dcount(2)) then
	  print*, 'Decision: single msg 2D x-distribution.'
	  WRITE (10,*) 'distr_x_2D= ',2
	else 
	  print*, 'Decision: collective 2D x-distribution.'
	  WRITE (10,*) 'distr_x_2D= ',3
	end if
      end if


c-------------------------------------------------------------------------------
c     Test 2D y-distribution     
c-------------------------------------------------------------------------------
      if (Master) print*, 'Standard 2D y-distribution test starts'
      do i=1, Ntests
        call MPI_BARRIER(MPI_COMM_WORLD, ierr)
	Tstart = MPI_WTIME()
	call distrib_y_2D_v1(ix,iy,ixloc, g2d,s2dy)
	if (YWorker) Tend = MPI_WTIME()
	call MPI_BARRIER(MPI_COMM_WORLD, ierr)
	if (Master) Tend = MPI_WTIME()
	testTimes(i,1) = Tend - Tstart
	if (YWorker) then
	  error = maxval(abs(
     &		g2d(owned_yslices(myid,1:no_of_yslices(myid)),:)
     &   	-s2dy(1:no_of_yslices(myid),:)))
          if (error .gt. 0.0) then
	    print*, 'Proc',myid,'After standard 2D y-distribution',
     &		' max error = ', error, ' check implementation '
          end if
        end if
	s2dy(:,:) = 0
      end do
      

      if (Master) print*, 'Single msg 2D y-distribution test starts'
      do i=1, Ntests
        call MPI_BARRIER(MPI_COMM_WORLD, ierr)
	Tstart = MPI_WTIME()
	call distrib_y_2D_v2(ix,iy,ixloc, g2d,s2dy)
	if (YWorker) Tend = MPI_WTIME()
	call MPI_BARRIER(MPI_COMM_WORLD, ierr)
	if (Master) Tend = MPI_WTIME()
	testTimes(i,2) = Tend - Tstart
	if (YWorker) then
	  error = maxval(abs(
     &		g2d(owned_yslices(myid,1:no_of_yslices(myid)),:)
     &   	-s2dy(1:no_of_yslices(myid),:)))
          if (error .gt. 0.0) then
	    print*, 'Proc',myid,'After single msg 2D y-distribution',
     &		' max error = ', error, ' check implementation '
          end if
        end if
	s2dy(:,:) = 0
      end do

      if (Master) print*, 'Collective 2D y-distribution test starts'
      do i=1, Ntests
        call MPI_BARRIER(MPI_COMM_WORLD, ierr)
	Tstart = MPI_WTIME()
	call distrib_y_2D_v3(ix,iy,ixloc, g2d,s2dy)
	if (YWorker) Tend = MPI_WTIME()
	call MPI_BARRIER(MPI_COMM_WORLD, ierr)
	if (Master) Tend = MPI_WTIME()
	testTimes(i,3) = Tend - Tstart
	if (YWorker) then
	  error = maxval(abs(
     &		g2d(owned_yslices(myid,1:no_of_yslices(myid)),:)
     &   	-s2dy(1:no_of_yslices(myid),:)))
          if (error .gt. 0.0) then
	    print*, 'Proc',myid,'After collective 2D y-distribution',
     &		' max error = ', error, ' check implementation '
          end if
        end if
	s2dy(:,:) = 0
      end do
      
      if (Master) then
        print*, 'Global standard 2D y-distribution time:',
     &		minval(testTimes(1:Ntests,1))
        print*, 'Global single msg 2D y-distribution time:',
     &		minval(testTimes(1:Ntests,2))
        print*, 'Global collective 2D y-distribution time:',
     &		minval(testTimes(1:Ntests,3))
      end if

      best = minval(testTimes(:,:))
      if (best .eq. minval(testTimes(1:Ntests,1))) then
	decision = 1
      else if (best .eq. minval(testTimes(1:Ntests,2))) then
        decision = 2
      else
	decision = 3
      end if
      
      call MPI_GATHER(decision, 1, MPI_INTEGER, mdecision, 1, 
     &		MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
      if (Master) then
        dcount(:) = 0
        do i=1,NYworkers+1
          if (mdecision(i) .eq. 1) then
	    dcount(1) = dcount(1) + 1
	  else if (mdecision(i) .eq. 2) then
	    dcount(2) = dcount(2) + 1
          else 
	    dcount(3) = dcount(3) + 1
	  end if
        end do
      
        best = maxval(dcount(:))
	if (best .eq. dcount(1)) then
	  print*, 'Decision: standard 2D y-distribution.'
	  WRITE (10,*) 'distr_y_2D= ',1
	else if (best .eq. dcount(2)) then
	  print*, 'Decision: single msg 2D y-distribution.'
	  WRITE (10,*) 'distr_y_2D= ',2
	else 
	  print*, 'Decision: collective 2D y-distribution.'
	  WRITE (10,*) 'distr_y_2D= ',3
	end if
      end if

c-------------------------------------------------------------------------------
c     Test BDz x-distribution     
c-------------------------------------------------------------------------------

      if (Master) print*, 'Standard BDz x-distribution test starts'
      do i=1, Ntests
        call MPI_BARRIER(MPI_COMM_WORLD, ierr)
	Tstart = MPI_WTIME()
        call distrib_x_BDz_v1(ix,iy,iyloc,is, gBDz,sBDzx)
	if (XWorker) Tend = MPI_WTIME()
	call MPI_BARRIER(MPI_COMM_WORLD, ierr)
	if (Master) Tend = MPI_WTIME()
	testTimes(i,1) = Tend - Tstart
	if (XWorker) then
	  error = maxval(abs(
     &		gBDz(:,owned_xslices(myid,1:no_of_xslices(myid)),:)
     &   	-sBDzx(:,1:no_of_xslices(myid),:)))
          if (error .gt. 0.0) then
	    print*, 'Proc',myid,'After standard BDz x-distribution',
     &		' max error = ', error, ' check implementation '
          end if
        end if
	sBDzx(:,:,:) = 0
      end do
      

      if (Master) print*, 'Single msg BDz x-distribution test starts'
      do i=1, Ntests
        call MPI_BARRIER(MPI_COMM_WORLD, ierr)
	Tstart = MPI_WTIME()
        call distrib_x_BDz_v2(ix,iy,iyloc,is, gBDz,sBDzx)
	if (XWorker) Tend = MPI_WTIME()
	call MPI_BARRIER(MPI_COMM_WORLD, ierr)
	if (Master) Tend = MPI_WTIME()
	testTimes(i,2) = Tend - Tstart
	if (XWorker) then
	  error = maxval(abs(
     &		gBDz(:,owned_xslices(myid,1:no_of_xslices(myid)),:)
     &   	-sBDzx(:,1:no_of_xslices(myid),:)))
          if (error .gt. 0.0) then
	    print*, 'Proc',myid,'After single msg BDz x-distribution',
     &		' max error = ', error, ' check implementation '
      	    print*,'sBDzx',sBDzx
          end if
        end if
	sBDzx(:,:,:) = 0
      end do

      if (Master) print*, 'Collective BDz x-distribution test starts'
      do i=1, Ntests
        call MPI_BARRIER(MPI_COMM_WORLD, ierr)
	Tstart = MPI_WTIME()
        call distrib_x_BDz_v3(ix,iy,iyloc,is, gBDz,sBDzx)
	if (XWorker) Tend = MPI_WTIME()
	call MPI_BARRIER(MPI_COMM_WORLD, ierr)
	if (Master) Tend = MPI_WTIME()
	testTimes(i,3) = Tend - Tstart
	if (XWorker) then
	  error = maxval(abs(
     &		gBDz(:,owned_xslices(myid,1:no_of_xslices(myid)),:)
     &   	-sBDzx(:,1:no_of_xslices(myid),:)))
          if (error .gt. 0.0) then
	    print*, 'Proc',myid,'After collective BDz x-distribution',
     &		' max error = ', error, ' check implementation '
          end if
        end if
	sBDzx(:,:,:) = 0
      end do
      
      if (Master) then
        print*, 'Global standard BDz x-distribution time:',
     &		minval(testTimes(1:Ntests,1))
        print*, 'Global single msg BDz x-distribution time:',
     &		minval(testTimes(1:Ntests,2))
        print*, 'Global collective BDz x-distribution time:',
     &		minval(testTimes(1:Ntests,3))
      end if

      best = minval(testTimes(:,:))
      if (best .eq. minval(testTimes(1:Ntests,1))) then
	decision = 1
      else if (best .eq. minval(testTimes(1:Ntests,2))) then
        decision = 2
      else
	decision = 3
      end if
      
      call MPI_GATHER(decision, 1, MPI_INTEGER, mdecision, 1, 
     &		MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
      if (Master) then
        dcount(:) = 0
        do i=1,NXworkers+1
          if (mdecision(i) .eq. 1) then
	    dcount(1) = dcount(1) + 1
	  else if (mdecision(i) .eq. 2) then
	    dcount(2) = dcount(2) + 1
          else 
	    dcount(3) = dcount(3) + 1
	  end if
        end do
      
        best = maxval(dcount(:))
	if (best .eq. dcount(1)) then
	  print*, 'Decision: standard BDz x-distribution.'
	  WRITE (10,*) 'distr_x_BDz= ',1
	else if (best .eq. dcount(2)) then
	  print*, 'Decision: single msg BDz x-distribution.'
	  WRITE (10,*) 'distr_x_BDz= ',2
	else 
	  print*, 'Decision: collective BDz x-distribution.'
	  WRITE (10,*) 'distr_x_BDz= ',3
	end if
      end if
      
c-------------------------------------------------------------------------------
c     Test BDz y-distribution     
c-------------------------------------------------------------------------------
      if (Master) print*, 'Standard BDz y-distribution test starts'
      do i=1, Ntests
        call MPI_BARRIER(MPI_COMM_WORLD, ierr)
	Tstart = MPI_WTIME()
        call distrib_y_BDz_v1(ix,iy,ixloc,is, gBDz,sBDzy)
	if (YWorker) Tend = MPI_WTIME()
	call MPI_BARRIER(MPI_COMM_WORLD, ierr)
	if (Master) Tend = MPI_WTIME()
	testTimes(i,1) = Tend - Tstart
	if (YWorker) then
	  error = maxval(abs(
     &		gBDz(owned_yslices(myid,1:no_of_yslices(myid)),:,:)
     &   	-sBDzy(1:no_of_yslices(myid),:,:)))
          if (error .gt. 0.0) then
	    print*, 'Proc',myid,'After standard BDz y-distribution',
     &		' max error = ', error, ' check implementation '
          end if
        end if
	sBDzy(:,:,:) = 0
      end do
      

      if (Master) print*, 'Single msg BDz y-distribution test starts'
      do i=1, Ntests
        call MPI_BARRIER(MPI_COMM_WORLD, ierr)
	Tstart = MPI_WTIME()
	
        call distrib_y_BDz_v2(ix,iy,ixloc,is, gBDz,sBDzy)
	if (YWorker) Tend = MPI_WTIME()
	call MPI_BARRIER(MPI_COMM_WORLD, ierr)
	if (Master) Tend = MPI_WTIME()
	testTimes(i,2) = Tend - Tstart
	if (YWorker) then
	  error = maxval(abs(
     &		gBDz(owned_yslices(myid,1:no_of_yslices(myid)),:,:)
     &   	-sBDzy(1:no_of_yslices(myid),:,:)))
          if (error .gt. 0.0) then
	    print*, 'Proc',myid,'After single msg BDz y-distribution',
     &		' max error = ', error, ' check implementation '
     	    print*,'sBDzy',sBDzy
          end if
        end if
	sBDzy(:,:,:) = 0
      end do

      if (Master) print*, 'Collective BDz y-distribution test starts'
      do i=1, Ntests
        call MPI_BARRIER(MPI_COMM_WORLD, ierr)
	Tstart = MPI_WTIME()
        call distrib_y_BDz_v3(ix,iy,ixloc,is, gBDz,sBDzy)
	if (YWorker) Tend = MPI_WTIME()
	call MPI_BARRIER(MPI_COMM_WORLD, ierr)
	if (Master) Tend = MPI_WTIME()
	testTimes(i,3) = Tend - Tstart
	if (YWorker) then
	  error = maxval(abs(
     &		gBDz(owned_yslices(myid,1:no_of_yslices(myid)),:,:)
     &   	-sBDzy(1:no_of_yslices(myid),:,:)))
          if (error .gt. 0.0) then
	    print*, 'Proc',myid,'After collective BDz y-distribution',
     &		' max error = ', error, ' check implementation '
          end if
        end if
	sBDzy(:,:,:) = 0
      end do
      
      if (Master) then
        print*, 'Global standard BDz y-distribution time:',
     &		minval(testTimes(1:Ntests,1))
        print*, 'Global single msg BDz y-distribution time:',
     &		minval(testTimes(1:Ntests,2))
        print*, 'Global collective BDz y-distribution time:',
     &		minval(testTimes(1:Ntests,3))
      end if

      best = minval(testTimes(:,:))
      if (best .eq. minval(testTimes(1:Ntests,1))) then
	decision = 1
      else if (best .eq. minval(testTimes(1:Ntests,2))) then
        decision = 2
      else
	decision = 3
      end if
      
      call MPI_GATHER(decision, 1, MPI_INTEGER, mdecision, 1, 
     &		MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
      if (Master) then
        dcount(:) = 0
        do i=1,NYworkers+1
          if (mdecision(i) .eq. 1) then
	    dcount(1) = dcount(1) + 1
	  else if (mdecision(i) .eq. 2) then
	    dcount(2) = dcount(2) + 1
          else 
	    dcount(3) = dcount(3) + 1
	  end if
        end do
      
        best = maxval(dcount(:))
	if (best .eq. dcount(1)) then
	  print*, 'Decision: standard BDz y-distribution.'
	  WRITE (10,*) 'distr_y_BDz= ',1
	else if (best .eq. dcount(2)) then
	  print*, 'Decision: single msg BDz y-distribution.'
	  WRITE (10,*) 'distr_y_BDz= ',2
	else 
	  print*, 'Decision: collective BDz y-distribution.'
	  WRITE (10,*) 'distr_y_BDz= ',3
	end if
      end if

c-------------------------------------------------------------------------------
c     Test 3D x-distribution     
c-------------------------------------------------------------------------------

      if (Master) print*, 'Standard 3D x-distribution test starts'
      do i=1, Ntests
        call MPI_BARRIER(MPI_COMM_WORLD, ierr)
	Tstart = MPI_WTIME()
        call distrib_x_3D_v1(ix,iy,iz,iyloc,g3d,s3dx)
	if (XWorker) Tend = MPI_WTIME()
	call MPI_BARRIER(MPI_COMM_WORLD, ierr)
	if (Master) Tend = MPI_WTIME()
	testTimes(i,1) = Tend - Tstart
	if (XWorker) then
	  error = maxval(abs(
     &		g3d(:,owned_xslices(myid,1:no_of_xslices(myid)),:)
     &   	-s3dx(:,1:no_of_xslices(myid),:)))
          if (error .gt. 0.0) then
	    print*, 'Proc',myid,'After standard 3D x-distribution',
     &		' max error = ', error, ' check implementation '
          end if
        end if
	s3dx(:,:,:) = 0
      end do
      

      if (Master) print*, 'Single msg 3D x-distribution test starts'
      do i=1, Ntests
        call MPI_BARRIER(MPI_COMM_WORLD, ierr)
	Tstart = MPI_WTIME()
        call distrib_x_3D_v2(ix,iy,iz,iyloc,g3d,s3dx)
	if (XWorker) Tend = MPI_WTIME()
	call MPI_BARRIER(MPI_COMM_WORLD, ierr)
	if (Master) Tend = MPI_WTIME()
	testTimes(i,2) = Tend - Tstart
	if (XWorker) then
	  error = maxval(abs(
     &		g3d(:,owned_xslices(myid,1:no_of_xslices(myid)),:)
     &   	-s3dx(:,1:no_of_xslices(myid),:)))
          if (error .gt. 0.0) then
	    print*, 'Proc',myid,'After single msg 3D x-distribution',
     &		' max error = ', error, ' check implementation '
          end if
        end if
	s3dx(:,:,:) = 0
      end do

      if (Master) print*, 'Collective 3D x-distribution test starts'
      do i=1, Ntests
        call MPI_BARRIER(MPI_COMM_WORLD, ierr)
	Tstart = MPI_WTIME()
        call distrib_x_3D_v3(ix,iy,iz,iyloc,g3d,s3dx)
	if (XWorker) Tend = MPI_WTIME()
	call MPI_BARRIER(MPI_COMM_WORLD, ierr)
	if (Master) Tend = MPI_WTIME()
	testTimes(i,3) = Tend - Tstart
	if (XWorker) then
	  error = maxval(abs(
     &		g3d(:,owned_xslices(myid,1:no_of_xslices(myid)),:)
     &   	-s3dx(:,1:no_of_xslices(myid),:)))
          if (error .gt. 0.0) then
	    print*, 'Proc',myid,'After collective 3D x-distribution',
     &		' max error = ', error, ' check implementation '
          end if
        end if
	s3dx(:,:,:) = 0
      end do
      
      if (Master) then
        print*, 'Global standard 3D x-distribution time:',
     &		minval(testTimes(1:Ntests,1))
        print*, 'Global single msg 3D x-distribution time:',
     &		minval(testTimes(1:Ntests,2))
        print*, 'Global collective 3D x-distribution time:',
     &		minval(testTimes(1:Ntests,3))
      end if

      best = minval(testTimes(:,:))
      if (best .eq. minval(testTimes(1:Ntests,1))) then
	decision = 1
      else if (best .eq. minval(testTimes(1:Ntests,2))) then
        decision = 2
      else
	decision = 3
      end if
      
      call MPI_GATHER(decision, 1, MPI_INTEGER, mdecision, 1, 
     &		MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
      if (Master) then
        dcount(:) = 0
        do i=1,NXworkers+1
          if (mdecision(i) .eq. 1) then
	    dcount(1) = dcount(1) + 1
	  else if (mdecision(i) .eq. 2) then
	    dcount(2) = dcount(2) + 1
          else 
	    dcount(3) = dcount(3) + 1
	  end if
        end do
      
        best = maxval(dcount(:))
	if (best .eq. dcount(1)) then
	  print*, 'Decision: standard 3D x-distribution.'
	  WRITE (10,*) 'distr_x_3D= ',1
	else if (best .eq. dcount(2)) then
	  print*, 'Decision: single msg 3D x-distribution.'
	  WRITE (10,*) 'distr_x_3D= ',2
	else 
	  print*, 'Decision: collective 3D x-distribution.'
	  WRITE (10,*) 'distr_x_3D= ',3
	end if
      end if
      
c-------------------------------------------------------------------------------
c     Test 3D y-distribution     
c-------------------------------------------------------------------------------
      if (Master) print*, 'Standard 3D y-distribution test starts'
      do i=1, Ntests
        call MPI_BARRIER(MPI_COMM_WORLD, ierr)
	Tstart = MPI_WTIME()
        call distrib_y_3D_v1(ix,iy,iz,ixloc,g3d,s3dy)
	if (YWorker) Tend = MPI_WTIME()
	call MPI_BARRIER(MPI_COMM_WORLD, ierr)
	if (Master) Tend = MPI_WTIME()
	testTimes(i,1) = Tend - Tstart
	if (YWorker) then
	  error = maxval(abs(
     &		g3d(owned_yslices(myid,1:no_of_yslices(myid)),:,:)
     &   	-s3dy(1:no_of_yslices(myid),:,:)))
          if (error .gt. 0.0) then
	    print*, 'Proc',myid,'After standard 3D y-distribution',
     &		' max error = ', error, ' check implementation '
          end if
        end if
	s3dy(:,:,:) = 0
      end do
      

      if (Master) print*, 'Single msg 3D y-distribution test starts'
      do i=1, Ntests
        call MPI_BARRIER(MPI_COMM_WORLD, ierr)
	Tstart = MPI_WTIME()
        call distrib_y_3D_v2(ix,iy,iz,ixloc,g3d,s3dy)
	if (YWorker) Tend = MPI_WTIME()
	call MPI_BARRIER(MPI_COMM_WORLD, ierr)
	if (Master) Tend = MPI_WTIME()
	testTimes(i,2) = Tend - Tstart
	if (YWorker) then
	  error = maxval(abs(
     &		g3d(owned_yslices(myid,1:no_of_yslices(myid)),:,:)
     &   	-s3dy(1:no_of_yslices(myid),:,:)))
          if (error .gt. 0.0) then
	    print*, 'Proc',myid,'After single msg 3D y-distribution',
     &		' max error = ', error, ' check implementation '
          end if
        end if
	s3dy(:,:,:) = 0
      end do

      if (Master) print*, 'Collective 3D y-distribution test starts'
      do i=1, Ntests
        call MPI_BARRIER(MPI_COMM_WORLD, ierr)
	Tstart = MPI_WTIME()
        call distrib_y_3D_v3(ix,iy,iz,ixloc,g3d,s3dy)
	if (YWorker) Tend = MPI_WTIME()
	call MPI_BARRIER(MPI_COMM_WORLD, ierr)
	if (Master) Tend = MPI_WTIME()
	testTimes(i,3) = Tend - Tstart
	if (YWorker) then
	  error = maxval(abs(
     &		g3d(owned_yslices(myid,1:no_of_yslices(myid)),:,:)
     &   	-s3dy(1:no_of_yslices(myid),:,:)))
          if (error .gt. 0.0) then
	    print*, 'Proc',myid,'After collective 3D y-distribution',
     &		' max error = ', error, ' check implementation '
          end if
        end if
	s3dy(:,:,:) = 0
      end do
      
      if (Master) then
        print*, 'Global standard 3D y-distribution time:',
     &		minval(testTimes(1:Ntests,1))
        print*, 'Global single msg 3D y-distribution time:',
     &		minval(testTimes(1:Ntests,2))
        print*, 'Global collective 3D y-distribution time:',
     &		minval(testTimes(1:Ntests,3))
      end if

      best = minval(testTimes(:,:))
      if (best .eq. minval(testTimes(1:Ntests,1))) then
	decision = 1
      else if (best .eq. minval(testTimes(1:Ntests,2))) then
        decision = 2
      else
	decision = 3
      end if
      
      call MPI_GATHER(decision, 1, MPI_INTEGER, mdecision, 1, 
     &		MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
      if (Master) then
        dcount(:) = 0
        do i=1,NYworkers+1
          if (mdecision(i) .eq. 1) then
	    dcount(1) = dcount(1) + 1
	  else if (mdecision(i) .eq. 2) then
	    dcount(2) = dcount(2) + 1
          else 
	    dcount(3) = dcount(3) + 1
	  end if
        end do
      
        best = maxval(dcount(:))
	if (best .eq. dcount(1)) then
	  print*, 'Decision: standard 3D y-distribution.'
	  WRITE (10,*) 'distr_y_3D= ',1
	else if (best .eq. dcount(2)) then
	  print*, 'Decision: single msg 3D y-distribution.'
	  WRITE (10,*) 'distr_y_3D= ',2
	else 
	  print*, 'Decision: collective 3D y-distribution.'
	  WRITE (10,*) 'distr_y_3D= ',3
	end if
      end if

c-------------------------------------------------------------------------------
c     Test BDy x-distribution     
c-------------------------------------------------------------------------------

      if (Master) print*, 'Standard BDx x-distribution test starts'
      do i=1, Ntests
        call MPI_BARRIER(MPI_COMM_WORLD, ierr)
	Tstart = MPI_WTIME()
        call distrib_x_BDx_v1(iy,iz,iyloc,is, gBDx,sBDx)
	if (XWorker) Tend = MPI_WTIME()
	call MPI_BARRIER(MPI_COMM_WORLD, ierr)
	if (Master) Tend = MPI_WTIME()
	testTimes(i,1) = Tend - Tstart
	if (XWorker) then
	  error = maxval(abs(
     &		gBDx(owned_xslices(myid,1:no_of_xslices(myid)),:,:)
     &   	-sBDx(1:no_of_xslices(myid),:,:)))
          if (error .gt. 0.0) then
	    print*, 'Proc',myid,'After standard BDx x-distribution',
     &		' max error = ', error, ' check implementation '
          end if
        end if
	sBDx(:,:,:) = 0
      end do
      

      if (Master) print*, 'Single msg BDx x-distribution test starts'
      do i=1, Ntests
        call MPI_BARRIER(MPI_COMM_WORLD, ierr)
	Tstart = MPI_WTIME()
        call distrib_x_BDx_v2(iy,iz,iyloc,is, gBDx,sBDx)
	if (XWorker) Tend = MPI_WTIME()
	call MPI_BARRIER(MPI_COMM_WORLD, ierr)
	if (Master) Tend = MPI_WTIME()
	testTimes(i,2) = Tend - Tstart
	if (XWorker) then
	  error = maxval(abs(
     &		gBDx(owned_xslices(myid,1:no_of_xslices(myid)),:,:)
     &   	-sBDx(1:no_of_xslices(myid),:,:)))
          if (error .gt. 0.0) then
	    print*, 'Proc',myid,'After single msg BDx x-distribution',
     &		' max error = ', error, ' check implementation '
          end if
        end if
	sBDx(:,:,:) = 0
      end do

      if (Master) print*, 'Collective BDx x-distribution test starts'
      do i=1, Ntests
        call MPI_BARRIER(MPI_COMM_WORLD, ierr)
	Tstart = MPI_WTIME()
        call distrib_x_BDx_v3(iy,iz,iyloc,is, gBDx,sBDx)
	if (XWorker) Tend = MPI_WTIME()
	call MPI_BARRIER(MPI_COMM_WORLD, ierr)
	if (Master) Tend = MPI_WTIME()
	testTimes(i,3) = Tend - Tstart
	if (XWorker) then
	  error = maxval(abs(
     &		gBDx(owned_xslices(myid,1:no_of_xslices(myid)),:,:)
     &   	-sBDx(1:no_of_xslices(myid),:,:)))
          if (error .gt. 0.0) then
	    print*, 'Proc',myid,'After collective BDx x-distribution',
     &		' max error = ', error, ' check implementation '
          end if
        end if
	sBDx(:,:,:) = 0
      end do
      
      if (Master) then
        print*, 'Global standard BDx x-distribution time:',
     &		minval(testTimes(1:Ntests,1))
        print*, 'Global single msg BDx x-distribution time:',
     &		minval(testTimes(1:Ntests,2))
        print*, 'Global collective BDx x-distribution time:',
     &		minval(testTimes(1:Ntests,3))
      end if

      best = minval(testTimes(:,:))
      if (best .eq. minval(testTimes(1:Ntests,1))) then
	decision = 1
      else if (best .eq. minval(testTimes(1:Ntests,2))) then
        decision = 2
      else
	decision = 3
      end if
      
      call MPI_GATHER(decision, 1, MPI_INTEGER, mdecision, 1, 
     &		MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
      if (Master) then
        dcount(:) = 0
        do i=1,NXworkers+1
          if (mdecision(i) .eq. 1) then
	    dcount(1) = dcount(1) + 1
	  else if (mdecision(i) .eq. 2) then
	    dcount(2) = dcount(2) + 1
          else 
	    dcount(3) = dcount(3) + 1
	  end if
        end do
      
        best = maxval(dcount(:))
	if (best .eq. dcount(1)) then
	  print*, 'Decision: standard BDx x-distribution.'
c	  WRITE (10,*) 'distr_x_BDx= ',1
	else if (best .eq. dcount(2)) then
	  print*, 'Decision: single msg BDx x-distribution.'
c	  WRITE (10,*) 'distr_x_BDx= ',2
	else 
	  print*, 'Decision: collective BDx x-distribution.'
c	  WRITE (10,*) 'distr_x_BDx= ',3
	end if
      end if
      
c-------------------------------------------------------------------------------
c     Test BD y-distribution     
c-------------------------------------------------------------------------------
      if (Master) print*, 'Standard BDy y-distribution test starts'
      do i=1, Ntests
        call MPI_BARRIER(MPI_COMM_WORLD, ierr)
	Tstart = MPI_WTIME()
        call distrib_y_BDy_v1(ix,iz,ixloc,is,gBDy,sBDy)
	if (YWorker) Tend = MPI_WTIME()
	call MPI_BARRIER(MPI_COMM_WORLD, ierr)
	if (Master) Tend = MPI_WTIME()
	testTimes(i,1) = Tend - Tstart
	if (YWorker) then
	  error = maxval(abs(
     &		gBDy(owned_yslices(myid,1:no_of_yslices(myid)),:,:)
     &   	-sBDy(1:no_of_yslices(myid),:,:)))
          if (error .gt. 0.0) then
	    print*, 'Proc',myid,'After standard BDy y-distribution',
     &		' max error = ', error, ' check implementation '
          end if
        end if
	sBDy(:,:,:) = 0
      end do
      

      if (Master) print*, 'Single msg BDy y-distribution test starts'
      do i=1, Ntests
        call MPI_BARRIER(MPI_COMM_WORLD, ierr)
	Tstart = MPI_WTIME()
        call distrib_y_BDy_v2(ix,iz,ixloc,is,gBDy,sBDy)
	if (YWorker) Tend = MPI_WTIME()
	call MPI_BARRIER(MPI_COMM_WORLD, ierr)
	if (Master) Tend = MPI_WTIME()
	testTimes(i,2) = Tend - Tstart
	if (YWorker) then
	  error = maxval(abs(
     &		gBDy(owned_yslices(myid,1:no_of_yslices(myid)),:,:)
     &   	-sBDy(1:no_of_yslices(myid),:,:)))
          if (error .gt. 0.0) then
	    print*, 'Proc',myid,'After single msg BDy y-distribution',
     &		' max error = ', error, ' check implementation '
          end if
        end if
	sBDy(:,:,:) = 0
      end do

      if (Master) print*, 'Collective BDy y-distribution test starts'
      do i=1, Ntests
        call MPI_BARRIER(MPI_COMM_WORLD, ierr)
	Tstart = MPI_WTIME()
        call distrib_y_BDy_v3(ix,iz,ixloc,is,gBDy,sBDy)
	if (YWorker) Tend = MPI_WTIME()
	call MPI_BARRIER(MPI_COMM_WORLD, ierr)
	if (Master) Tend = MPI_WTIME()
	testTimes(i,3) = Tend - Tstart
	if (YWorker) then
	  error = maxval(abs(
     &		gBDy(owned_yslices(myid,1:no_of_yslices(myid)),:,:)
     &  	-sBDy(1:no_of_yslices(myid),:,:)))
          if (error .gt. 0.0) then
	    print*, 'Proc',myid,'After collective BDy y-distribution',
     &		' max error = ', error, ' check implementation '
          end if
        end if
	sBDy(:,:,:) = 0
      end do
      
      if (Master) then
        print*, 'Global standard BDy y-distribution time:',
     &		minval(testTimes(1:Ntests,1))
        print*, 'Global single msg BDy y-distribution time:',
     &		minval(testTimes(1:Ntests,2))
        print*, 'Global collective BDy y-distribution time:',
     &		minval(testTimes(1:Ntests,3))
      end if

      best = minval(testTimes(:,:))
      if (best .eq. minval(testTimes(1:Ntests,1))) then
	decision = 1
      else if (best .eq. minval(testTimes(1:Ntests,2))) then
        decision = 2
      else
	decision = 3
      end if
      
      call MPI_GATHER(decision, 1, MPI_INTEGER, mdecision, 1, 
     &		MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
      if (Master) then
        dcount(:) = 0
        do i=1,NYworkers+1
          if (mdecision(i) .eq. 1) then
	    dcount(1) = dcount(1) + 1
	  else if (mdecision(i) .eq. 2) then
	    dcount(2) = dcount(2) + 1
          else 
	    dcount(3) = dcount(3) + 1
	  end if
        end do
      
        best = maxval(dcount(:))
	if (best .eq. dcount(1)) then
	  print*, 'Decision: standard BDy y-distribution.'
c	  WRITE (10,*) 'distr_y_BDy= ',1
	else if (best .eq. dcount(2)) then
	  print*, 'Decision: single msg BDy y-distribution.'
c	  WRITE (10,*) 'distr_y_BDy= ',2
	else 
	  print*, 'Decision: collective BDy y-distribution.'
c	  WRITE (10,*) 'distr_y_BDy= ',3
	end if
	end if
c-------------------------------------------------------------------------------
c     Test BD x-distribution     
c-------------------------------------------------------------------------------

      if (Master) print*, 'Standard BD x-distribution test starts'
      do i=1, Ntests
        call MPI_BARRIER(MPI_COMM_WORLD, ierr)
	Tstart = MPI_WTIME()
        call distrib_x_BD_v1(iy,iz,iyloc,is, gBDx2,sBDx2)
	if (XWorker) Tend = MPI_WTIME()
	call MPI_BARRIER(MPI_COMM_WORLD, ierr)
	if (Master) Tend = MPI_WTIME()
	testTimes(i,1) = Tend - Tstart
	if (XWorker) then
	  error = maxval(abs(
     &		gBDx2(owned_xslices(myid,1:no_of_xslices(myid)),:,:,:)
     &   	-sBDx2(1:no_of_xslices(myid),:,:,:)))
          if (error .gt. 0.0) then
	    print*, 'Proc',myid,'After standard BD x-distribution',
     &		' max error = ', error, ' check implementation '
          end if
        end if
	sBDx2(:,:,:,:) = 0
      end do
      

      if (Master) print*, 'Single msg BD x-distribution test starts'
      do i=1, Ntests
        call MPI_BARRIER(MPI_COMM_WORLD, ierr)
	Tstart = MPI_WTIME()
        call distrib_x_BD_v2(iy,iz,iyloc,is, gBDx2,sBDx2)
	if (XWorker) Tend = MPI_WTIME()
	call MPI_BARRIER(MPI_COMM_WORLD, ierr)
	if (Master) Tend = MPI_WTIME()
	testTimes(i,2) = Tend - Tstart
	if (XWorker) then
	  error = maxval(abs(
     &		gBDx2(owned_xslices(myid,1:no_of_xslices(myid)),:,:,:)
     &   	-sBDx2(1:no_of_xslices(myid),:,:,:)))
          if (error .gt. 0.0) then
	    print*, 'Proc',myid,'After single msg BD x-distribution',
     &		' max error = ', error, ' check implementation '
          end if
        end if
	sBDx2(:,:,:,:) = 0
      end do

      if (Master) print*, 'Collective BD x-distribution test starts'
      do i=1, Ntests
        call MPI_BARRIER(MPI_COMM_WORLD, ierr)
	Tstart = MPI_WTIME()
        call distrib_x_BD_v3(iy,iz,iyloc,is, gBDx2,sBDx2)
	if (XWorker) Tend = MPI_WTIME()
	call MPI_BARRIER(MPI_COMM_WORLD, ierr)
	if (Master) Tend = MPI_WTIME()
	testTimes(i,3) = Tend - Tstart
	if (XWorker) then
	  error = maxval(abs(
     &		gBDx2(owned_xslices(myid,1:no_of_xslices(myid)),:,:,:)
     &   	-sBDx2(1:no_of_xslices(myid),:,:,:)))
          if (error .gt. 0.0) then
	    print*, 'Proc',myid,'After collective BD x-distribution',
     &		' max error = ', error, ' check implementation '
          end if
        end if
	sBDx2(:,:,:,:) = 0
      end do
      
      if (Master) then
        print*, 'Global standard BD x-distribution time:',
     &		minval(testTimes(1:Ntests,1))
        print*, 'Global single msg BD x-distribution time:',
     &		minval(testTimes(1:Ntests,2))
        print*, 'Global collective BD x-distribution time:',
     &		minval(testTimes(1:Ntests,3))
      end if

      best = minval(testTimes(:,:))
      if (best .eq. minval(testTimes(1:Ntests,1))) then
	decision = 1
      else if (best .eq. minval(testTimes(1:Ntests,2))) then
        decision = 2
      else
	decision = 3
      end if
      
      call MPI_GATHER(decision, 1, MPI_INTEGER, mdecision, 1, 
     &		MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
      if (Master) then
        dcount(:) = 0
        do i=1,NXworkers+1
          if (mdecision(i) .eq. 1) then
	    dcount(1) = dcount(1) + 1
	  else if (mdecision(i) .eq. 2) then
	    dcount(2) = dcount(2) + 1
          else 
	    dcount(3) = dcount(3) + 1
	  end if
        end do
      
        best = maxval(dcount(:))
	if (best .eq. dcount(1)) then
	  print*, 'Decision: standard BD x-distribution.'
	  WRITE (10,*) 'distr_x_BD= ',1
	else if (best .eq. dcount(2)) then
	  print*, 'Decision: single msg BD x-distribution.'
	  WRITE (10,*) 'distr_x_BD= ',2
	else 
	  print*, 'Decision: collective BD x-distribution.'
	  WRITE (10,*) 'distr_x_BD= ',3
	end if
      end if
      
c-------------------------------------------------------------------------------
c     Test BD y-distribution     
c-------------------------------------------------------------------------------
      if (Master) print*, 'Standard BD y-distribution test starts'
      do i=1, Ntests
        call MPI_BARRIER(MPI_COMM_WORLD, ierr)
	Tstart = MPI_WTIME()
        call distrib_y_BD_v1(ix,iz,ixloc,is,gBDy2,sBDy2)
	if (YWorker) Tend = MPI_WTIME()
	call MPI_BARRIER(MPI_COMM_WORLD, ierr)
	if (Master) Tend = MPI_WTIME()
	testTimes(i,1) = Tend - Tstart
	if (YWorker) then
	  error = maxval(abs(
     &		gBDy2(owned_yslices(myid,1:no_of_yslices(myid)),:,:,:)
     &   	-sBDy2(1:no_of_yslices(myid),:,:,:)))
          if (error .gt. 0.0) then
	    print*, 'Proc',myid,'After standard BD y-distribution',
     &		' max error = ', error, ' check implementation '
          end if
        end if
	sBDy2(:,:,:,:) = 0
      end do
      

      if (Master) print*, 'Single msg BD y-distribution test starts'
      do i=1, Ntests
        call MPI_BARRIER(MPI_COMM_WORLD, ierr)
	Tstart = MPI_WTIME()
        call distrib_y_BD_v2(ix,iz,ixloc,is,gBDy2,sBDy2)
	if (YWorker) Tend = MPI_WTIME()
	call MPI_BARRIER(MPI_COMM_WORLD, ierr)
	if (Master) Tend = MPI_WTIME()
	testTimes(i,2) = Tend - Tstart
	if (YWorker) then
	  error = maxval(abs(
     &		gBDy2(owned_yslices(myid,1:no_of_yslices(myid)),:,:,:)
     &   	-sBDy2(1:no_of_yslices(myid),:,:,:)))
          if (error .gt. 0.0) then
	    print*, 'Proc',myid,'After single msg BD y-distribution',
     &		' max error = ', error, ' check implementation '
          end if
        end if
	sBDy2(:,:,:,:) = 0
      end do

      if (Master) print*, 'Collective BD y-distribution test starts'
      do i=1, Ntests
        call MPI_BARRIER(MPI_COMM_WORLD, ierr)
	Tstart = MPI_WTIME()
        call distrib_y_BD_v3(ix,iz,ixloc,is,gBDy2,sBDy2)
	if (YWorker) Tend = MPI_WTIME()
	call MPI_BARRIER(MPI_COMM_WORLD, ierr)
	if (Master) Tend = MPI_WTIME()
	testTimes(i,3) = Tend - Tstart
	if (YWorker) then
	  error = maxval(abs(
     &		gBDy2(owned_yslices(myid,1:no_of_yslices(myid)),:,:,:)
     &  	-sBDy2(1:no_of_yslices(myid),:,:,:)))
          if (error .gt. 0.0) then
	    print*, 'Proc',myid,'After collective BD y-distribution',
     &		' max error = ', error, ' check implementation '
          end if
        end if
	sBDy2(:,:,:,:) = 0
      end do
      
      if (Master) then
        print*, 'Global standard BD y-distribution time:',
     &		minval(testTimes(1:Ntests,1))
        print*, 'Global single msg BD y-distribution time:',
     &		minval(testTimes(1:Ntests,2))
        print*, 'Global collective BD y-distribution time:',
     &		minval(testTimes(1:Ntests,3))
      end if

      best = minval(testTimes(:,:))
      if (best .eq. minval(testTimes(1:Ntests,1))) then
	decision = 1
      else if (best .eq. minval(testTimes(1:Ntests,2))) then
        decision = 2
      else
	decision = 3
      end if
      
      call MPI_GATHER(decision, 1, MPI_INTEGER, mdecision, 1, 
     &		MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
      if (Master) then
        dcount(:) = 0
        do i=1,NYworkers+1
          if (mdecision(i) .eq. 1) then
	    dcount(1) = dcount(1) + 1
	  else if (mdecision(i) .eq. 2) then
	    dcount(2) = dcount(2) + 1
          else 
	    dcount(3) = dcount(3) + 1
	  end if
        end do
      
        best = maxval(dcount(:))
	if (best .eq. dcount(1)) then
	  print*, 'Decision: standard BD y-distribution.'
	  WRITE (10,*) 'distr_y_BD= ',1
	else if (best .eq. dcount(2)) then
	  print*, 'Decision: single msg BD y-distribution.'
	  WRITE (10,*) 'distr_y_BD= ',2
	else 
	  print*, 'Decision: collective BD y-distribution.'
	  WRITE (10,*) 'distr_y_BD= ',3
	end if
	
	CLOSE (10)
      end if

           
      call MPI_FINALIZE(rc)
      stop
      end
