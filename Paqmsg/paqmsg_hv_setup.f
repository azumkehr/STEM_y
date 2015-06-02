      program mpi_communication_setup
      use HVParallelDataMap
      use HVCommDataTypes
      use HVCommunicationLibrary
      
      implicit none
      real, pointer, dimension(:,:,:,:) :: sg1
      real, pointer, dimension(:,:) :: g2d
      real, pointer, dimension(:,:,:) :: gBDz
      real, pointer, dimension(:,:,:) :: g3d
      real, pointer, dimension(:,:,:) :: gBDhx
      real, pointer, dimension(:,:,:) :: gBDhy
      real, pointer, dimension(:,:,:,:) :: gBDhx2
      real, pointer, dimension(:,:,:,:) :: gBDhy2

      real, pointer, dimension(:,:,:,:) ::  sg2
      real, pointer, dimension(:,:,:,:) :: sgh
      real, pointer, dimension(:,:,:,:) :: sgv
      real, pointer, dimension(:,:) :: s2dh, s2dv
      real, pointer, dimension(:,:,:) :: sBDzv
      real, pointer, dimension(:,:,:) :: s3dh, s3dv
      real, pointer, dimension(:,:,:) :: sBDhx, sBDhy
      real, pointer, dimension(:,:,:,:) :: sBDhx2, sBDhy2
      integer ix, iy, iz, is
      integer  i, i1, j, k, l, rc, ierr
      integer izloc, icloc
      real, dimension(:), pointer :: errc
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
      call CreateMap(ix,iy,iz,izloc, icloc)
c      if (HWorker .or. VWorker) 
c      print*, 'izloc=',izloc,'  icloc=',icloc
c      print*, 'P[',myid,'] has hslices=',owned_hslices(myid,:)
c      print*, 'P[',myid,'] has vcols=',owned_vcols(myid,:)

c-------------------------------------------------------------------------------
c    Create Special Communication Patterns
c-------------------------------------------------
      call CreateCommDataTypes(ix,iy,iz,izloc,icloc,is)
       
c-------------------------------------------------------------------------------
c    Allocate the arrays
c-------------------------------------------------------------------------------
      allocate(sg1(ix,iy,iz,is),STAT=ierr)
      allocate(g2d(ix,iy),STAT=ierr)
      allocate(gBDz(ix,iy,is),STAT=ierr)
      allocate(g3d(ix,iy,iz),STAT=ierr)
      allocate(gBDhx(iy,iz,is),STAT=ierr)
      allocate(gBDhy(ix,iz,is),STAT=ierr)
      allocate(gBDhx2(ix,iz,2,is),STAT=ierr)
      allocate(gBDhy2(iy,iz,2,is),STAT=ierr)

      if (Master) then
        allocate(sg2(ix,iy,iz,is),STAT=ierr)
        allocate(s2dh(ix,iy),STAT=ierr)

        allocate(sgh(1,1,1,1),STAT=ierr)
        allocate(sgv(1,1,1,1),STAT=ierr)
        allocate(s2dv(1,1),STAT=ierr)
        allocate(sBDzv(1,1,1),STAT=ierr)
        allocate(s3dh(1,1,1),STAT=ierr)
        allocate(s3dv(1,1,1),STAT=ierr)
        allocate(sBDhx2(1,1,2,1),STAT=ierr)
        allocate(sBDhy2(1,1,2,1),STAT=ierr)
        allocate(sBDhx(1,1,1),STAT=ierr)
        allocate(sBDhy(1,1,1),STAT=ierr)
      else	
        allocate(sg2(1,1,1,1),STAT=ierr)
	
        allocate(sgh(ix,iy,izloc,is),STAT=ierr)
        allocate(sgv(1,icloc,iz,is),STAT=ierr)

        allocate(s2dh(ix,iy),STAT=ierr)
        allocate(s2dv(1,icloc),STAT=ierr)
        allocate(sBDzv(1,icloc,is),STAT=ierr)
        allocate(s3dh(ix,iy,izloc),STAT=ierr)
        allocate(s3dv(1,icloc,iz),STAT=ierr)
        allocate(sBDhx(iy,izloc,is),STAT=ierr)
        allocate(sBDhy(ix,izloc,is),STAT=ierr)
        allocate(sBDhx2(ix,izloc,2,is),STAT=ierr)
        allocate(sBDhy2(iy,izloc,2,is),STAT=ierr)
      endif

      if (VWorker) then
        allocate(errc(no_of_vcols(myid)),STAT=ierr)
      else 
        allocate(errc(1),STAT=ierr)
      end if
      
c-------------------------------------------------------------------------------
c    Master Initializes the concentration array - all initialize it here...
c-------------------------------------------------------------------------------
        do i=1,ix
	  do j=1,iy
	    do k=1,iz
	      do l=1,is
	        sg1(i,j,k,l) = l-1 + is*(k-1) + is*iz*(j-1) + is*iz*iy*(i-1)
		g2d(i,j) = j-1 + iy*(i-1)
		gBDz(i,j,l) = l-1 + is*(j-1) + is*iy*(i-1)
		g3d(i,j,k) = k-1 + iz*(j-1) + iz*iy*(i-1)
		gBDhx(i,k,l) = l-1 + is*(k-1) + is*iz*(j-1)
		gBDhx2(j,k,1,l) = l-1 + is*(k-1) + is*iz*(j-1)
		gBDhx2(j,k,2,l) = is*iz*iy + l-1 + is*(k-1) + is*iz*(j-1)
		gBDhy(j,k,l) = l-1 + is*(k-1) + is*iz*(i-1)
		gBDhy2(i,k,1,l) = l-1 + is*(k-1) + is*iz*(i-1)
		gBDhy2(i,k,2,l) = is*iz*ix + l-1 + is*(k-1) + is*iz*(i-1)
              end do		
            end do		
          end do		
        end do	
	if (Master) then ! for broadcasting the array...
	  s2dh(:,:) = g2d(:,:)
	end if
c-------------------------------------------------------------------------------
c     Testing communication methods starts here     
c-------------------------------------------------------------------------------

      if (Master) then
        allocate(mdecision(Nprocs), STAT=ierr)
      else
        allocate(mdecision(1), STAT=ierr)
      end if

c-------------------------------------------------------------------------------
c     Test 4D h-distribution     
c-------------------------------------------------------------------------------
      if (Master) print*, 'Standard 4D h-distribution test starts'
      do i=1, Ntests
        call MPI_BARRIER(MPI_COMM_WORLD, ierr)
	Tstart = MPI_WTIME()
	call distrib_h_4D_v1(ix,iy,iz,izloc,is, sg1,sgh)
	if (HWorker) Tend = MPI_WTIME()
	call MPI_BARRIER(MPI_COMM_WORLD, ierr)
	if (Master) Tend = MPI_WTIME()
	testTimes(i,1) = Tend - Tstart
	if (HWorker) then
	  error = maxval(abs(
     &		sg1(:,:,owned_hslices(myid,1:no_of_hslices(myid)),:)
     &   	-sgh(:,:,1:no_of_hslices(myid),:)))
          if (error .gt. 0.0) then
	    print*,'P',myid, 'After standard 4D h-distribution',
     &		' max error = ', error, ' check implementation '
          end if
        end if
	sgh(:,:,:,:) = 0
      end do
      
      if (Master) print*, 'Single msg 4D h-distribution test starts'
      do i=1, Ntests
        call MPI_BARRIER(MPI_COMM_WORLD, ierr)
	Tstart = MPI_WTIME()
	call distrib_h_4D_v2(ix,iy,iz,izloc,is, sg1,sgh)
	if (HWorker) Tend = MPI_WTIME()
	call MPI_BARRIER(MPI_COMM_WORLD, ierr)
	if (Master) Tend = MPI_WTIME()
	testTimes(i,2) = Tend - Tstart
	if (HWorker) then
	  error = maxval(abs(
     &		sg1(:,:,owned_hslices(myid,1:no_of_hslices(myid)),:)
     &   	-sgh(:,:,1:no_of_hslices(myid),:)))
          if (error .gt. 0.0) then
	    print*,'P',myid, 'After single msg 4D h-distribution',
     &		' max error = ', error, ' check implementation '
          end if
        end if
	sgh(:,:,:,:) = 0
      end do

      if (Master) print*, 'Collective 4D h-distribution test starts'
      do i=1, Ntests
        call MPI_BARRIER(MPI_COMM_WORLD, ierr)
	Tstart = MPI_WTIME()
	call distrib_h_4D_v3(ix,iy,iz,izloc,is, sg1,sgh)
	if (HWorker) Tend = MPI_WTIME()
	call MPI_BARRIER(MPI_COMM_WORLD, ierr)
	if (Master) Tend = MPI_WTIME()
	testTimes(i,3) = Tend - Tstart
	if (HWorker) then
	  error = maxval(abs(
     &		sg1(:,:,owned_hslices(myid,1:no_of_hslices(myid)),:)
     &   	-sgh(:,:,1:no_of_hslices(myid),:)))
          if (error .gt. 0.0) then
	    print*,'P',myid, 'After collective 4D h-distribution',
     &		' max error = ', error, ' check implementation '
          end if
        end if
	sgh(:,:,:,:) = 0
      end do
      
      if (Master) then
        print*, 'Global standard 4D h-distribution time:',
     &		minval(testTimes(1:Ntests,1))
        print*, 'Global single msg 4D h-distribution time:',
     &		minval(testTimes(1:Ntests,2))
        print*, 'Global collective 4D h-distribution time:',
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
        do i=1,NHworkers+1
          if (mdecision(i) .eq. 1) then
	    dcount(1) = dcount(1) + 1
	  else if (mdecision(i) .eq. 2) then
	    dcount(2) = dcount(2) + 1
          else 
	    dcount(3) = dcount(3) + 1
	  end if
        end do

c	write to file HV_VERSIONS
        OPEN (UNIT=10, FILE='HV_versions')
      
        best = maxval(dcount(:))
	if (best .eq. dcount(1)) then
	  print*, 'Decision: standard 4D h-distribution.'
	  WRITE (10,*) 'distr_h_4D= ',1
	else if (best .eq. dcount(2)) then
	  print*, 'Decision: single msg 4D h-distribution.'
	  WRITE (10,*) 'distr_h_4D= ',2
	else 
	  print*, 'Decision: collective 4D h-distribution.'
	  WRITE (10,*) 'distr_h_4D= ',3
	end if
      end if
      
c-------------------------------------------------------------------------------
c     Test 4D v-distribution     
c-------------------------------------------------------------------------------

      if (Master) print*, 'Standard 4D v-distribution test starts'
      do i1=1, Ntests
        call MPI_BARRIER(MPI_COMM_WORLD, ierr)
	Tstart = MPI_WTIME()
	call distrib_v_4D_v1(ix,iy,iz,icloc,is, sg1,sgv)
	if (VWorker) Tend = MPI_WTIME()
	call MPI_BARRIER(MPI_COMM_WORLD, ierr)
	if (Master) Tend = MPI_WTIME()
	testTimes(i1,1) = Tend - Tstart
	if (VWorker) then
  	  do k=1,no_of_vcols(MyId)
            i = planar_vcol_id(global_vcol_id(myid,k),1)
    	    j =planar_vcol_id(global_vcol_id(myid,k),2)
            errc(k)=maxval(abs(sg1(i,j,:,:)-sgv(1,k,:,:)))
	  enddo
	  error = maxval(errc)
          if (error .gt. 0.0) then
	    print*,'P',myid, 'After standard 4D v-distribution',
     &		' max error = ', error, ' check implementation '
          end if
        end if
	sgv(1,:,:,:) = 0
      end do
      

      if (Master) print*, 'Single msg 4D v-distribution test starts'
      do i1=1, Ntests
        call MPI_BARRIER(MPI_COMM_WORLD, ierr)
	Tstart = MPI_WTIME()
	call distrib_v_4D_v2(ix,iy,iz,icloc,is, sg1,sgv)
	if (VWorker) Tend = MPI_WTIME()
	call MPI_BARRIER(MPI_COMM_WORLD, ierr)
	if (Master) Tend = MPI_WTIME()
	testTimes(i1,2) = Tend - Tstart
	if (VWorker) then
  	  do k=1,no_of_vcols(MyId)
            i = planar_vcol_id(global_vcol_id(myid,k),1)
    	    j =planar_vcol_id(global_vcol_id(myid,k),2)
            errc(k)=maxval(abs(sg1(i,j,:,:)-sgv(1,k,:,:)))
	  enddo
	  error = maxval(errc)
          if (error .gt. 0.0) then
	    print*,'P',myid, 'After single msg 4D v-distribution',
     &		' max error = ', error, ' check implementation '
          end if
        end if
	sgv(:,:,:,:) = 0
      end do

      if (Master) print*, 'Collective 4D v-distribution test starts'
      do i1=1, Ntests
        call MPI_BARRIER(MPI_COMM_WORLD, ierr)
	Tstart = MPI_WTIME()
	call distrib_v_4D_v3(ix,iy,iz,icloc,is, sg1,sgv)
	if (VWorker) Tend = MPI_WTIME()
	call MPI_BARRIER(MPI_COMM_WORLD, ierr)
	if (Master) Tend = MPI_WTIME()
	testTimes(i1,3) = Tend - Tstart
	if (VWorker) then
  	  do k=1,no_of_vcols(MyId)
            i = planar_vcol_id(global_vcol_id(myid,k),1)
    	    j =planar_vcol_id(global_vcol_id(myid,k),2)
            errc(k)=maxval(abs(sg1(i,j,:,:)-sgv(1,k,:,:)))
	  enddo
	  error = maxval(errc)
          if (error .gt. 0.0) then
	    print*,'P',myid, 'After collective 4D v-distribution',
     &		' max error = ', error, ' check implementation '
          end if
        end if
	sgv(:,:,:,:) = 0
      end do
      
      if (Master) then
        print*, 'Global standard 4D v-distribution time:',
     &		minval(testTimes(1:Ntests,1))
        print*, 'Global single msg 4D v-distribution time:',
     &		minval(testTimes(1:Ntests,2))
        print*, 'Global collective 4D v-distribution time:',
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
        do i=1,NVworkers+1
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
	  print*, 'Decision: standard 4D v-distribution.'
	  WRITE (10,*) 'distr_v_4D= ',1
	else if (best .eq. dcount(2)) then
	  print*, 'Decision: single msg 4D v-distribution.'
	  WRITE (10,*) 'distr_v_4D= ',2
	else 
	  print*, 'Decision: collective 4D v-distribution.'
	  WRITE (10,*) 'distr_v_4D= ',3
	end if
      end if

c-------------------------------------------------------------------------------
c     Test 4D h-gather     
c-------------------------------------------------------------------------------


      if (Master) print*, 'Standard 4D h-gather test starts'
c     in order to check the error we need to call distribution first
      call distrib_h_4D_v1(ix,iy,iz,izloc,is, sg1,sgh)
      
      do i=1, Ntests
        call MPI_BARRIER(MPI_COMM_WORLD, ierr)
	Tstart = MPI_WTIME()
	call gather_h_4D_v1(ix,iy,iz,izloc,is, sg2,sgh)
	if (HWorker) Tend = MPI_WTIME()
	call MPI_BARRIER(MPI_COMM_WORLD, ierr)
	if (Master) Tend = MPI_WTIME()
	testTimes(i,1) = Tend - Tstart
	if (Master) then
	  error = maxval(abs(sg2-sg1))
          if (error .gt. 0.0) then
	    print*,'P',myid, 'After standard 4D h-gather',
     &		' max error = ', error, ' check implementation '
          end if
        end if
	sg2(:,:,:,:) = 0
      end do
      

      if (Master) print*, 'Single msg 4D h-gather test starts'
      do i=1, Ntests
        call MPI_BARRIER(MPI_COMM_WORLD, ierr)
	Tstart = MPI_WTIME()
	call gather_h_4D_v2(ix,iy,iz,izloc,is, sg2,sgh)
	if (HWorker) Tend = MPI_WTIME()
	call MPI_BARRIER(MPI_COMM_WORLD, ierr)
	if (Master) Tend = MPI_WTIME()
	testTimes(i,2) = Tend - Tstart
	if (Master) then
	  error = maxval(abs(sg2-sg1))
          if (error .gt. 0.0) then
	    print*,'P',myid, 'After single msg 4D h-gather',
     &		' max error = ', error, ' check implementation '
          end if
        end if
	sg2(:,:,:,:) = 0
      end do

      if (Master) print*, 'Collective 4D h-gather test starts'
      do i=1, Ntests
        call MPI_BARRIER(MPI_COMM_WORLD, ierr)
	Tstart = MPI_WTIME()
	call gather_h_4D_v3(ix,iy,iz,izloc,is, sg2,sgh)
	if (HWorker) Tend = MPI_WTIME()
	call MPI_BARRIER(MPI_COMM_WORLD, ierr)
	if (Master) Tend = MPI_WTIME()
	testTimes(i,3) = Tend - Tstart
	if (Master) then
	  error = maxval(abs(sg2-sg1))
          if (error .gt. 0.0) then
	    print*,'P',myid, 'After collective 4D h-gather',
     &		'max error = ', error, ' check implementation '
          end if
        end if
	sg2(:,:,:,:) = 0
      end do
      
      if (Master) then
        print*, 'Global standard 4D h-gather. time:',
     &		minval(testTimes(1:Ntests,1))
        print*, 'Global single msg 4D h-gather time:',
     &		minval(testTimes(1:Ntests,2))
        print*, 'Global collective 4D h-gather time:',
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
        do i=1,NHworkers+1
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
	  print*, 'Decision: standard 4D h-gather.'
	  WRITE (10,*) 'gather_h_4D= ',1
	else if (best .eq. dcount(2)) then
	  print*, 'Decision: single msg 4D h-gather.'
	  WRITE (10,*) 'gather_h_4D= ',2
	else 
	  print*, 'Decision: collective 4D h-gather.'
	  WRITE (10,*) 'gather_h_4D= ',3
	end if
      end if

c-------------------------------------------------------------------------------
c     Test 4D v-gather     
c-------------------------------------------------------------------------------
      if (Master) print*, 'Standard 4D v-gather test starts'
c     in order to check the error we need to call distribution first
      call distrib_v_4D_v1(ix,iy,iz,icloc,is, sg1,sgv)

      do i=1, Ntests
        call MPI_BARRIER(MPI_COMM_WORLD, ierr)
	Tstart = MPI_WTIME()
	call gather_v_4D_v1(ix,iy,iz,icloc,is, sg2,sgv)
	if (VWorker) Tend = MPI_WTIME()
	call MPI_BARRIER(MPI_COMM_WORLD, ierr)
	if (Master) Tend = MPI_WTIME()
	testTimes(i,1) = Tend - Tstart
	if (Master) then
	  error = maxval(abs(sg2-sg1))
          if (error .gt. 0.0) then
	    print*,'P',myid, 'After standard 4D v-gather',
     &		' max error = ', error, ' check implementation '
          end if
        end if
	sg2(:,:,:,:) = 0
      end do
      

      if (Master) print*, 'Single msg 4D v-gather test starts'
      do i=1, Ntests
        call MPI_BARRIER(MPI_COMM_WORLD, ierr)
	Tstart = MPI_WTIME()
	call gather_v_4D_v2(ix,iy,iz,icloc,is, sg2,sgv)
	if (VWorker) Tend = MPI_WTIME()
	call MPI_BARRIER(MPI_COMM_WORLD, ierr)
	if (Master) Tend = MPI_WTIME()
	testTimes(i,2) = Tend - Tstart
	if (Master) then
	  error = maxval(abs(sg2-sg1))
          if (error .gt. 0.0) then
	    print*,'P',myid, 'After single msg 4D v-gather',
     &		' max error = ', error, ' check implementation '
          end if
        end if
	sg2(:,:,:,:) = 0
      end do

      if (Master) print*, 'Collective 4D v-gather test starts'
      do i=1, Ntests
        call MPI_BARRIER(MPI_COMM_WORLD, ierr)
	Tstart = MPI_WTIME()
	call gather_v_4D_v3(ix,iy,iz,icloc,is, sg2,sgv)
	if (VWorker) Tend = MPI_WTIME()
	call MPI_BARRIER(MPI_COMM_WORLD, ierr)
	if (Master) Tend = MPI_WTIME()
	testTimes(i,3) = Tend - Tstart
	if (Master) then
	  error = maxval(abs(sg2-sg1))
          if (error .gt. 0.0) then
	    print*,'P',myid, 'After collective 4D v-gather',
     &		' max error = ', error, ' check implementation '
          end if
        end if
	sg2(:,:,:,:) = 0
      end do
      
      if (Master) then
        print*, 'Global standard 4D v-gather time:',
     &		minval(testTimes(1:Ntests,1))
        print*, 'Global single msg 4D v-gather. time:',
     &		minval(testTimes(1:Ntests,2))
        print*, 'Global collective 4D v-gather time:',
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
        do i=1,NVworkers+1
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
	  print*, 'Decision: standard 4D v-gather.'
	  WRITE (10,*) 'gather_v_4D= ',1
	else if (best .eq. dcount(2)) then
	  print*, 'Decision: single msg 4D v-gather.'
	  WRITE (10,*) 'gather_v_4D= ',2
	else 
	  print*, 'Decision: collective 4D v-gather.'
	  WRITE (10,*) 'gather_v_4D= ',3
	end if
      end if

	
c-------------------------------------------------------------------------------
c     Test h2v shuffle    
c-------------------------------------------------------------------------------

      if (Master) print*, 'Standard 4D h2v shuffle test starts'
c     in order to check the error we need to call distribution first
      call distrib_h_4D_v1(ix,iy,iz,izloc,is, sg1,sgh)
      sgv(:,:,:,:) = 0

      do i1=1, Ntests
        call MPI_BARRIER(MPI_COMM_WORLD, ierr)
	Tstart = MPI_WTIME()
        call shuffle_h2v_4D_v1(ix,iy,iz,izloc,icloc,is,sgh,sgv)
	Tend = MPI_WTIME()
        call MPI_BARRIER(MPI_COMM_WORLD, ierr)
	if (Master) Tend = MPI_WTIME()
	testTimes(i1,1) = Tend - Tstart
        if (VWorker) then
	  do k=1,no_of_vcols(MyId)
	    i=planar_vcol_id(global_vcol_id(myid,k),1)
	    j=planar_vcol_id(global_vcol_id(myid,k),2)
            errc(k)=maxval(abs(sg1(i,j,:,:)-sgv(1,k,:,:)))
	  enddo
	  error = maxval(errc)
          if (error .gt. 0.0) then
	    print*,'P',myid, 'After standard 4D h2v shuffle',
     &		' max error = ', error, ' check implementation '
          end if
        end if
	sgv(:,:,:,:) = 0
      end do
      

      if (Master) print*, 'Single msg 4D h2v shuffle test starts'
      do i1=1, Ntests
        call MPI_BARRIER(MPI_COMM_WORLD, ierr)
	Tstart = MPI_WTIME()
        call shuffle_h2v_4D_v2(ix,iy,iz,izloc,icloc,is,sgh,sgv)
	Tend = MPI_WTIME()
        call MPI_BARRIER(MPI_COMM_WORLD, ierr)
	if (Master) Tend = MPI_WTIME()
	testTimes(i1,2) = Tend - Tstart
        if (VWorker) then
	  do k=1,no_of_vcols(MyId)
	    i=planar_vcol_id(global_vcol_id(myid,k),1)
	    j=planar_vcol_id(global_vcol_id(myid,k),2)
            errc(k)=maxval(abs(sg1(i,j,:,:)-sgv(1,k,:,:)))
	  enddo
	  error = maxval(errc)
          if (error .gt. 0.0) then
	    print*,'P',myid, 'After single msg 4D h2v shuffle',
     &		' max error = ', error, ' check implementation '
          end if
        end if
	sgv(:,:,:,:) = 0
      end do


      if (Master) print*, 'Collective 4D h2v shuffle test starts'
      do i1=1, Ntests
        call MPI_BARRIER(MPI_COMM_WORLD, ierr)
	Tstart = MPI_WTIME()
        call shuffle_h2v_4D_v3(ix,iy,iz,izloc,icloc,is,sgh,sgv)
	Tend = MPI_WTIME()
        call MPI_BARRIER(MPI_COMM_WORLD, ierr)
	if (Master) Tend = MPI_WTIME()
	testTimes(i1,3) = Tend - Tstart
        if (VWorker) then
	  do k=1,no_of_vcols(MyId)
	    i=planar_vcol_id(global_vcol_id(myid,k),1)
	    j=planar_vcol_id(global_vcol_id(myid,k),2)
            errc(k)=maxval(abs(sg1(i,j,:,:)-sgv(1,k,:,:)))
	  enddo
	  error = maxval(errc)
          if (error .gt. 0.0) then
	    print*,'P',myid, 'After collective 4D h2v shuffle',
     &		' max error = ', error, ' check implementation '
          end if
        end if
	sgv(:,:,:,:) = 0
      end do
      
	      
      if (Master) then
        print*, 'Global standard 4D h2v shuffle time:',
     &		minval(testTimes(1:Ntests,1))
        print*, 'Global single msg 4D h2v shuffle time:',
     &		minval(testTimes(1:Ntests,2))
        print*, 'Global collective 4D h2v shuffle time:',
     &		minval(testTimes(1:Ntests,3))
      end if
      if (MyId .eq. 1) then
        print*, 'Proc 1 standard 4D h2v shuffle time:',
     &		minval(testTimes(1:Ntests,1))
        print*, 'Proc 1 single msg 4D h2v shuffle time:',
     &		minval(testTimes(1:Ntests,2))
        print*, 'Proc 1 collective 4D h2v shuffle time:',
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
        do i=2,max(NHworkers,NVworkers)+1
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
	  print*, 'Decision: standard 4D h2v shuffle.'
	  WRITE (10,*) 'shuffle_h2v_4D= ',1
	else if (best .eq. dcount(2)) then
	  print*, 'Decision: single msg 4D h2v shuffle.'
	  WRITE (10,*) 'shuffle_h2v_4D= ',2
	else 
	  print*, 'Decision: collective 4D h2v shuffle.'
	  WRITE (10,*) 'shuffle_h2v_4D= ',3
	end if
      end if

c-------------------------------------------------------------------------------
c     Test v2h shuffle    
c-------------------------------------------------------------------------------
      if (Master) print*, 'Standard 4D v2h shuffle test starts'
c     in order to check the error we need to call distribution first
      call distrib_v_4D_v1(ix,iy,iz,icloc,is, sg1,sgv)
      sgh(:,:,:,:) = 0

      do i=1, Ntests
        call MPI_BARRIER(MPI_COMM_WORLD, ierr)
	Tstart = MPI_WTIME()
        call shuffle_v2h_4D_v1(ix,iy,iz,izloc,icloc,is,sgh,sgv)
	Tend = MPI_WTIME()
        call MPI_BARRIER(MPI_COMM_WORLD, ierr)
	if (Master) Tend = MPI_WTIME()
	testTimes(i,1) = Tend - Tstart
	if (HWorker) then
	  error = maxval(abs(
     &		sg1(:,:,owned_hslices(myid,1:no_of_hslices(myid)),:)
     & 		-sgh(:,:,1:no_of_hslices(myid),:)))
          if (error .gt. 0.0) then
	    print*,'P',myid, 'After standard 4D v2h shuffle',
     &		' max error = ', error, ' check implementation '
          end if
        end if
	sgh(:,:,:,:) = 0
      end do
      

      if (Master) print*, 'Single msg 4D v2h shuffle test starts'
      do i=1, Ntests
        call MPI_BARRIER(MPI_COMM_WORLD, ierr)
	Tstart = MPI_WTIME()
        call shuffle_v2h_4D_v2(ix,iy,iz,izloc,icloc,is,sgh,sgv)
	Tend = MPI_WTIME()
        call MPI_BARRIER(MPI_COMM_WORLD, ierr)
	if (Master) Tend = MPI_WTIME()
	testTimes(i,2) = Tend - Tstart
	if (HWorker) then
	  error = maxval(abs(
     &		sg1(:,:,owned_hslices(myid,1:no_of_hslices(myid)),:)
     & 		-sgh(:,:,1:no_of_hslices(myid),:)))
          if (error .gt. 0.0) then
	    print*,'P',myid, 'After single msg 4D v2h shuffle',
     &		' max error = ', error, ' check implementation '
          end if
        end if
	sgh(:,:,:,:) = 0
      end do

      if (Master) print*, 'Collective 4D v2h shuffle test starts'
      do i=1, Ntests
        call MPI_BARRIER(MPI_COMM_WORLD, ierr)
	Tstart = MPI_WTIME()
        call shuffle_v2h_4D_v3(ix,iy,iz,izloc,icloc,is,sgh,sgv)
	Tend = MPI_WTIME()
        call MPI_BARRIER(MPI_COMM_WORLD, ierr)
	if (Master) Tend = MPI_WTIME()
	testTimes(i,3) = Tend - Tstart
	if (HWorker) then
	  error = maxval(abs(
     &		sg1(:,:,owned_hslices(myid,1:no_of_hslices(myid)),:)
     & 		-sgh(:,:,1:no_of_hslices(myid),:)))
          if (error .gt. 0.0) then
	    print*,'P',myid, 'After collective 4D v2h shuffle',
     &		' max error = ', error, ' check implementation '
          end if
        end if
	sgh(:,:,:,:) = 0
      end do
      
	      
      if (Master) then
        print*, 'Global standard 4D v2h shuffle time:',
     &		minval(testTimes(1:Ntests,1))
        print*, 'Global single msg 4D v2h shuffle time:',
     &		minval(testTimes(1:Ntests,2))
        print*, 'Global collective 4D v2h shuffle time:',
     &		minval(testTimes(1:Ntests,3))
      end if
      if (MyId .eq. 1) then
        print*, 'Proc 1 standard 4D v2h shuffle time:',
     &		minval(testTimes(1:Ntests,1))
        print*, 'Proc 1 single msg 4D v2h shuffle time:',
     &		minval(testTimes(1:Ntests,2))
        print*, 'Proc 1 collective 4D v2h shuffle time:',
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
        do i=2,max(NHworkers,NVworkers)+1
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
	  print*, 'Decision: standard 4D v2h shuffle.'
	  WRITE (10,*) 'shuffle_v2h_4D= ',1
	else if (best .eq. dcount(2)) then
	  print*, 'Decision: single msg 4D v2h shuffle.'
	  WRITE (10,*) 'shuffle_v2h_4D= ',2
	else 
	  print*, 'Decision: collective 4D v2h shuffle.'
	  WRITE (10,*) 'shuffle_v2h_4D= ',3
	end if
      end if

c-------------------------------------------------------------------------------
c     Test 2D h-distribution     
c-------------------------------------------------------------------------------

      if (Master) print*, 'Standard 2D h-distribution test starts'
      do i=1, Ntests
        call MPI_BARRIER(MPI_COMM_WORLD, ierr)
	Tstart = MPI_WTIME()
	call distrib_h_2D_v1(ix,iy,s2dh)
	Tend = MPI_WTIME()
	testTimes(i,1) = Tend - Tstart
	if (HWorker) then
	  error = maxval(abs(g2d(:,:)-s2dh(:,:)))
          if (error .gt. 0.0) then
	    print*,'P',myid, 'After standard 2D h-distribution',
     &		' max error = ', error, ' check implementation '
          end if
	  s2dh(:,:) = 0
        end if
      end do
      
      if (Master) then
        print*, 'Global standard 2D h-distribution time:',
     &		minval(testTimes(1:Ntests,1))
      end if

      if (Master) then
	  print*, 'Decision: use the only version: 
     &	standard 2D h-distribution.'
	  WRITE (10,*) 'distr_h_2D= ',1
     
      end if


c-------------------------------------------------------------------------------
c     Test 2D v-distribution     
c-------------------------------------------------------------------------------
      if (Master) print*, 'Standard 2D v-distribution test starts'
      do i1=1, Ntests
        call MPI_BARRIER(MPI_COMM_WORLD, ierr)
	Tstart = MPI_WTIME()
	call distrib_v_2D_v1(ix,iy,icloc, g2d,s2dv)
	if (VWorker) Tend = MPI_WTIME()
	call MPI_BARRIER(MPI_COMM_WORLD, ierr)
	if (Master) Tend = MPI_WTIME()
	testTimes(i1,1) = Tend - Tstart
	if (VWorker) then
  	  do k=1,no_of_vcols(MyId)
	    i=planar_vcol_id(global_vcol_id(myid,k),1)
	    j=planar_vcol_id(global_vcol_id(myid,k),2)
            errc(k)=(abs(g2d(i,j)-s2dv(1,k)))
	  enddo
	  error = maxval(errc)
          if (error .gt. 0.0) then
	    print*,'P',myid, 'After standard 2D v-distribution',
     &		' max error = ', error, ' check implementation '
          end if
        end if
	s2dv(:,:) = 0
      end do
      

      if (Master) print*, 'Single msg 2D v-distribution test starts'
      do i1=1, Ntests
        call MPI_BARRIER(MPI_COMM_WORLD, ierr)
	Tstart = MPI_WTIME()
	call distrib_v_2D_v2(ix,iy,icloc, g2d,s2dv)
	if (VWorker) Tend = MPI_WTIME()
	call MPI_BARRIER(MPI_COMM_WORLD, ierr)
	if (Master) Tend = MPI_WTIME()
	testTimes(i1,2) = Tend - Tstart
	if (VWorker) then
  	  do k=1,no_of_vcols(MyId)
	    i=planar_vcol_id(global_vcol_id(myid,k),1)
	    j=planar_vcol_id(global_vcol_id(myid,k),2)
            errc(k)=(abs(g2d(i,j)-s2dv(1,k)))
	  enddo
	  error = maxval(errc)
          if (error .gt. 0.0) then
	    print*,'P',myid, 'After single msg 2D v-distribution',
     &		' max error = ', error, ' check implementation '
          end if
        end if
	s2dv(:,:) = 0
      end do

      if (Master) print*, 'Collective 2D v-distribution test starts'
      do i1=1, Ntests
        call MPI_BARRIER(MPI_COMM_WORLD, ierr)
	Tstart = MPI_WTIME()
	call distrib_v_2D_v3(ix,iy,icloc, g2d,s2dv)
	if (VWorker) Tend = MPI_WTIME()
	call MPI_BARRIER(MPI_COMM_WORLD, ierr)
	if (Master) Tend = MPI_WTIME()
	testTimes(i1,3) = Tend - Tstart
	if (VWorker) then
  	  do k=1,no_of_vcols(MyId)
	    i=planar_vcol_id(global_vcol_id(myid,k),1)
	    j=planar_vcol_id(global_vcol_id(myid,k),2)
            errc(k)=(abs(g2d(i,j)-s2dv(1,k)))
	  enddo
	  error = maxval(errc)
          if (error .gt. 0.0) then
	    print*,'P',myid, 'After collective 2D v-distribution',
     &		' max error = ', error, ' check implementation '
          end if
        end if
	s2dv(:,:) = 0
      end do
      
      if (Master) then
        print*, 'Global standard 2D v-distribution time:',
     &		minval(testTimes(1:Ntests,1))
        print*, 'Global single msg 2D v-distribution time:',
     &		minval(testTimes(1:Ntests,2))
        print*, 'Global collective 2D v-distribution time:',
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
        do i=1,NVworkers+1
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
	  print*, 'Decision: standard 2D v-distribution.'
	  WRITE (10,*) 'distr_v_2D= ',1
	else if (best .eq. dcount(2)) then
	  print*, 'Decision: single msg 2D v-distribution.'
	  WRITE (10,*) 'distr_v_2D= ',2
	else 
	  print*, 'Decision: collective 2D v-distribution.'
	  WRITE (10,*) 'distr_v_2D= ',3
	end if
      end if

c-------------------------------------------------------------------------------
c     Test BDz v-distribution     
c-------------------------------------------------------------------------------
      if (Master) print*, 'Standard BDz v-distribution test starts'
      do i1=1, Ntests
        call MPI_BARRIER(MPI_COMM_WORLD, ierr)
	Tstart = MPI_WTIME()
        call distrib_v_BDz_v1(ix,iy,icloc,is, gBDz,sBDzv)
	if (VWorker) Tend = MPI_WTIME()
	call MPI_BARRIER(MPI_COMM_WORLD, ierr)
	if (Master) Tend = MPI_WTIME()
	testTimes(i1,1) = Tend - Tstart
	if (VWorker) then
	  do k=1,no_of_vcols(MyId)
 	    i=planar_vcol_id(global_vcol_id(myid,k),1)
	    j=planar_vcol_id(global_vcol_id(myid,k),2)
            errc(k)=maxval(abs(gBDz(i,j,:)-sBDzv(1,k,:)))
	  enddo
	  error = maxval(abs(errc))
          if (error .gt. 0.0) then
	    print*,'P',myid, 'After standard BDz v-distribution',
     &		' max error = ', error, ' check implementation '
          end if
        end if
	sBDzv(:,:,:) = 0
      end do
      

      if (Master) print*, 'Single msg BDz v-distribution test starts'
      do i1=1, Ntests
        call MPI_BARRIER(MPI_COMM_WORLD, ierr)
	Tstart = MPI_WTIME()
        call distrib_v_BDz_v2(ix,iy,icloc,is, gBDz,sBDzv)
	if (VWorker) Tend = MPI_WTIME()
	call MPI_BARRIER(MPI_COMM_WORLD, ierr)
	if (Master) Tend = MPI_WTIME()
	testTimes(i1,2) = Tend - Tstart
	if (VWorker) then
	  do k=1,no_of_vcols(MyId)
 	    i=planar_vcol_id(global_vcol_id(myid,k),1)
	    j=planar_vcol_id(global_vcol_id(myid,k),2)
            errc(k)=maxval(abs(gBDz(i,j,:)-sBDzv(1,k,:)))
	  enddo
	  error = maxval(abs(errc))
          if (error .gt. 0.0) then
	    print*,'P',myid, 'After single msg BDz v-distribution',
     &		' max error = ', error, ' check implementation '
          end if
        end if
	sBDzv(:,:,:) = 0
      end do

      if (Master) print*, 'Collective BDz v-distribution test starts'
      do i1=1, Ntests
        call MPI_BARRIER(MPI_COMM_WORLD, ierr)
	Tstart = MPI_WTIME()
        call distrib_v_BDz_v3(ix,iy,icloc,is, gBDz,sBDzv)
	if (VWorker) Tend = MPI_WTIME()
	call MPI_BARRIER(MPI_COMM_WORLD, ierr)
	if (Master) Tend = MPI_WTIME()
	testTimes(i1,3) = Tend - Tstart
	if (VWorker) then
	  do k=1,no_of_vcols(MyId)
 	    i=planar_vcol_id(global_vcol_id(myid,k),1)
	    j=planar_vcol_id(global_vcol_id(myid,k),2)
            errc(k)=maxval(abs(gBDz(i,j,:)-sBDzv(1,k,:)))
	  enddo
	  error = maxval(abs(errc))
          if (error .gt. 0.0) then
	    print*,'P',myid, 'After collective BDz v-distribution',
     &		' max error = ', error, ' check implementation '
          end if
        end if
	sBDzv(:,:,:) = 0
      end do
      
      if (Master) then
        print*, 'Global standard BDz v-distribution time:',
     &		minval(testTimes(1:Ntests,1))
        print*, 'Global single msg BDz v-distribution time:',
     &		minval(testTimes(1:Ntests,2))
        print*, 'Global collective BDz v-distribution time:',
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
        do i=1,NVworkers+1
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
	  print*, 'Decision: standard BDz v-distribution.'
	  WRITE (10,*) 'distr_v_BDz= ',1
	else if (best .eq. dcount(2)) then
	  print*, 'Decision: single msg BDz v-distribution.'
	  WRITE (10,*) 'distr_v_BDz= ',2
	else 
	  print*, 'Decision: collective BDz v-distribution.'
	  WRITE (10,*) 'distr_v_BDz= ',3
	end if
      end if

c-------------------------------------------------------------------------------
c     Test 3D h-distribution     
c-------------------------------------------------------------------------------

      if (Master) print*, 'Standard 3D h-distribution test starts'
      do i=1, Ntests
        call MPI_BARRIER(MPI_COMM_WORLD, ierr)
	Tstart = MPI_WTIME()
        call distrib_h_3D_v1(ix,iy,iz,izloc,g3d,s3dh)
	if (HWorker) Tend = MPI_WTIME()
	call MPI_BARRIER(MPI_COMM_WORLD, ierr)
	if (Master) Tend = MPI_WTIME()
	testTimes(i,1) = Tend - Tstart
	if (HWorker) then
	  error = maxval(abs(
     &		g3d(:,:,owned_hslices(myid,1:no_of_hslices(myid)))
     &   	-s3dh(:,:,1:no_of_hslices(myid))))
          if (error .gt. 0.0) then
	    print*,'P',myid, 'After standard 3D h-distribution',
     &		' max error = ', error, ' check implementation '
          end if
        end if
	s3dh(:,:,:) = 0
      end do
      

      if (Master) print*, 'Single msg 3D h-distribution test starts'
      do i=1, Ntests
        call MPI_BARRIER(MPI_COMM_WORLD, ierr)
	Tstart = MPI_WTIME()
        call distrib_h_3D_v2(ix,iy,iz,izloc,g3d,s3dh)
	if (HWorker) Tend = MPI_WTIME()
	call MPI_BARRIER(MPI_COMM_WORLD, ierr)
	if (Master) Tend = MPI_WTIME()
	testTimes(i,2) = Tend - Tstart
	if (HWorker) then
	  error = maxval(abs(
     &		g3d(:,:,owned_hslices(myid,1:no_of_hslices(myid)))
     &   	-s3dh(:,:,1:no_of_hslices(myid))))
          if (error .gt. 0.0) then
	    print*,'P',myid, 'After single msg 3D h-distribution',
     &		' max error = ', error, ' check implementation '
          end if
        end if
	s3dh(:,:,:) = 0
      end do

      if (Master) print*, 'Collective 3D h-distribution test starts'
      do i=1, Ntests
        call MPI_BARRIER(MPI_COMM_WORLD, ierr)
	Tstart = MPI_WTIME()
        call distrib_h_3D_v3(ix,iy,iz,izloc,g3d,s3dh)
	if (HWorker) Tend = MPI_WTIME()
	call MPI_BARRIER(MPI_COMM_WORLD, ierr)
	if (Master) Tend = MPI_WTIME()
	testTimes(i,3) = Tend - Tstart
	if (HWorker) then
	  error = maxval(abs(
     &		g3d(:,:,owned_hslices(myid,1:no_of_hslices(myid)))
     &   	-s3dh(:,:,1:no_of_hslices(myid))))
          if (error .gt. 0.0) then
	    print*,'P',myid, 'After collective 3D h-distribution',
     &		' max error = ', error, ' check implementation '
          end if
        end if
	s3dh(:,:,:) = 0
      end do
      
      if (Master) then
        print*, 'Global standard 3D h-distribution time:',
     &		minval(testTimes(1:Ntests,1))
        print*, 'Global single msg 3D h-distribution time:',
     &		minval(testTimes(1:Ntests,2))
        print*, 'Global collective 3D h-distribution time:',
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
        do i=1,NHworkers+1
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
	  print*, 'Decision: standard 3D h-distribution.'
	  WRITE (10,*) 'distr_h_3D= ',1
	else if (best .eq. dcount(2)) then
	  print*, 'Decision: single msg 3D h-distribution.'
	  WRITE (10,*) 'distr_h_3D= ',2
	else 
	  print*, 'Decision: collective 3D h-distribution.'
	  WRITE (10,*) 'distr_h_3D= ',3
	end if
      end if
      
c-------------------------------------------------------------------------------
c     Test 3D v-distribution     
c-------------------------------------------------------------------------------
      if (Master) print*, 'Standard 3D v-distribution test starts'
      do i1=1, Ntests
        call MPI_BARRIER(MPI_COMM_WORLD, ierr)
	Tstart = MPI_WTIME()
        call distrib_v_3D_v1(ix,iy,iz,icloc,g3d,s3dv)
	if (VWorker) Tend = MPI_WTIME()
	call MPI_BARRIER(MPI_COMM_WORLD, ierr)
	if (Master) Tend = MPI_WTIME()
	testTimes(i1,1) = Tend - Tstart
        if (VWorker) then
	  do k=1,no_of_vcols(MyId)
	    i=planar_vcol_id(global_vcol_id(myid,k),1)
	    j=planar_vcol_id(global_vcol_id(myid,k),2)
            errc(k)=maxval(abs(g3d(i,j,:)-s3dv(1,k,:)))
	  enddo
	  error = maxval(abs(errc))
          if (error .gt. 0.0) then
	    print*,'P',myid, 'After standard 3D v-distribution',
     &		' max error = ', error, ' check implementation '
          end if
        end if
	s3dv(:,:,:) = 0
      end do
      

      if (Master) print*, 'Single msg 3D v-distribution test starts'
      do i1=1, Ntests
        call MPI_BARRIER(MPI_COMM_WORLD, ierr)
	Tstart = MPI_WTIME()
        call distrib_v_3D_v2(ix,iy,iz,icloc,g3d,s3dv)
	if (VWorker) Tend = MPI_WTIME()
	call MPI_BARRIER(MPI_COMM_WORLD, ierr)
	if (Master) Tend = MPI_WTIME()
	testTimes(i1,2) = Tend - Tstart
        if (VWorker) then
	  do k=1,no_of_vcols(MyId)
	    i=planar_vcol_id(global_vcol_id(myid,k),1)
	    j=planar_vcol_id(global_vcol_id(myid,k),2)
            errc(k)=maxval(abs(g3d(i,j,:)-s3dv(1,k,:)))
	  enddo
	  error = maxval(abs(errc))
          if (error .gt. 0.0) then
	    print*,'P',myid, 'After single msg 3D v-distribution',
     &		' max error = ', error, ' check implementation '
          end if
        end if
	s3dv(:,:,:) = 0
      end do

      if (Master) print*, 'Collective 3D v-distribution test starts'
      do i1=1, Ntests
        call MPI_BARRIER(MPI_COMM_WORLD, ierr)
	Tstart = MPI_WTIME()
        call distrib_v_3D_v3(ix,iy,iz,icloc,g3d,s3dv)
	if (VWorker) Tend = MPI_WTIME()
	call MPI_BARRIER(MPI_COMM_WORLD, ierr)
	if (Master) Tend = MPI_WTIME()
	testTimes(i1,3) = Tend - Tstart
        if (VWorker) then
	  do k=1,no_of_vcols(MyId)
	    i=planar_vcol_id(global_vcol_id(myid,k),1)
	    j=planar_vcol_id(global_vcol_id(myid,k),2)
            errc(k)=maxval(abs(g3d(i,j,:)-s3dv(1,k,:)))
	  enddo
	  error = maxval(abs(errc))
          if (error .gt. 0.0) then
	    print*,'P',myid, 'After collective 3D v-distribution',
     &		' max error = ', error, ' check implementation '
          end if
        end if
	s3dv(:,:,:) = 0
      end do
      
      if (Master) then
        print*, 'Global standard 3D v-distribution time:',
     &		minval(testTimes(1:Ntests,1))
        print*, 'Global single msg 3D v-distribution time:',
     &		minval(testTimes(1:Ntests,2))
        print*, 'Global collective 3D v-distribution time:',
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
        do i=1,NVworkers+1
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
	  print*, 'Decision: standard 3D v-distribution.'
	  WRITE (10,*) 'distr_v_3D= ',1
	else if (best .eq. dcount(2)) then
	  print*, 'Decision: single msg 3D v-distribution.'
	  WRITE (10,*) 'distr_v_3D= ',2
	else 
	  print*, 'Decision: collective 3D v-distribution.'
	  WRITE (10,*) 'distr_v_3D= ',3
	end if
      end if

c-------------------------------------------------------------------------------
c     Test BDx h-distribution     
c-------------------------------------------------------------------------------

      if (Master) print*, 'Standard BDx h-distribution test starts'
      do i=1, Ntests
        call MPI_BARRIER(MPI_COMM_WORLD, ierr)
	Tstart = MPI_WTIME()
        call distrib_h_BDx_v1(iy,iz,izloc,is, gBDhx,sBDhx)
	if (HWorker) Tend = MPI_WTIME()
	call MPI_BARRIER(MPI_COMM_WORLD, ierr)
	if (Master) Tend = MPI_WTIME()
	testTimes(i,1) = Tend - Tstart
	if (HWorker) then
	  error = maxval(abs(
     &		gBDhx(:,owned_hslices(myid,1:no_of_hslices(myid)),:)
     &   	-sBDhx(:,1:no_of_hslices(myid),:)))
          if (error .gt. 0.0) then
	    print*,'P',myid, 'After standard BDx h-distribution',
     &		' max error = ', error, ' check implementation '
          end if
        end if
	sBDhx(:,:,:) = 0
      end do
      

      if (Master) print*, 'Single msg BDx h-distribution test starts'
      do i=1, Ntests
        call MPI_BARRIER(MPI_COMM_WORLD, ierr)
	Tstart = MPI_WTIME()
        call distrib_h_BDx_v2(iy,iz,izloc,is, gBDhx,sBDhx)
	if (HWorker) Tend = MPI_WTIME()
	call MPI_BARRIER(MPI_COMM_WORLD, ierr)
	if (Master) Tend = MPI_WTIME()
	testTimes(i,2) = Tend - Tstart
	if (HWorker) then
	  error = maxval(abs(
     &		gBDhx(:,owned_hslices(myid,1:no_of_hslices(myid)),:)
     &   	-sBDhx(:,1:no_of_hslices(myid),:)))
          if (error .gt. 0.0) then
	    print*,'P',myid, 'After single msg BDx h-distribution',
     &		' max error = ', error, ' check implementation '
          end if
        end if
	sBDhx(:,:,:) = 0
      end do

      if (Master) print*, 'Collective BDx h-distribution test starts'
      do i=1, Ntests
        call MPI_BARRIER(MPI_COMM_WORLD, ierr)
	Tstart = MPI_WTIME()
        call distrib_h_BDx_v3(iy,iz,izloc,is, gBDhx,sBDhx)
	if (HWorker) Tend = MPI_WTIME()
	call MPI_BARRIER(MPI_COMM_WORLD, ierr)
	if (Master) Tend = MPI_WTIME()
	testTimes(i,3) = Tend - Tstart
	if (HWorker) then
	  error = maxval(abs(
     &		gBDhx(:,owned_hslices(myid,1:no_of_hslices(myid)),:)
     &   	-sBDhx(:,1:no_of_hslices(myid),:)))
          if (error .gt. 0.0) then
	    print*,'P',myid, 'After collective BDx h-distribution',
     &		' max error = ', error, ' check implementation '
          end if
        end if
	sBDhx(:,:,:) = 0
      end do
      
      if (Master) then
        print*, 'Global standard BDx h-distribution time:',
     &		minval(testTimes(1:Ntests,1))
        print*, 'Global single msg BDx h-distribution time:',
     &		minval(testTimes(1:Ntests,2))
        print*, 'Global collective BDx h-distribution time:',
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
        do i=1,NHworkers+1
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
	  print*, 'Decision: standard BDx h-distribution.'
c	  WRITE (10,*) 'distr_h_BDx= ',1
	else if (best .eq. dcount(2)) then
	  print*, 'Decision: single msg BDx h-distribution.'
c	  WRITE (10,*) 'distr_h_BDx= ',2
	else 
	  print*, 'Decision: collective BDx h-distribution.'
c	  WRITE (10,*) 'distr_h_BDx= ',3
	end if
      end if
      
c-------------------------------------------------------------------------------
c     Test BD hy-distribution     
c-------------------------------------------------------------------------------
      if (Master) print*, 'Standard BDy h-distribution test starts'
      do i=1, Ntests
        call MPI_BARRIER(MPI_COMM_WORLD, ierr)
	Tstart = MPI_WTIME()
        call distrib_h_BDy_v1(ix,iz,izloc,is,gBDhy,sBDhy)
	if (HWorker) Tend = MPI_WTIME()
	call MPI_BARRIER(MPI_COMM_WORLD, ierr)
	if (Master) Tend = MPI_WTIME()
	testTimes(i,1) = Tend - Tstart
	if (HWorker) then
	  error = maxval(abs(
     &		gBDhy(:,owned_hslices(myid,1:no_of_hslices(myid)),:)
     &   	-sBDhy(:,1:no_of_hslices(myid),:)))
          if (error .gt. 0.0) then
	    print*, 'P',myid,'After standard BDy h-distribution',
     &		' max error = ', error, ' check implementation '
          end if
        end if
	sBDhy(:,:,:) = 0
      end do
      

      if (Master) print*, 'Single msg BDy h-distribution test starts'
      do i=1, Ntests
        call MPI_BARRIER(MPI_COMM_WORLD, ierr)
	Tstart = MPI_WTIME()
        call distrib_h_BDy_v2(ix,iz,izloc,is,gBDhy,sBDhy)
	if (HWorker) Tend = MPI_WTIME()
	call MPI_BARRIER(MPI_COMM_WORLD, ierr)
	if (Master) Tend = MPI_WTIME()
	testTimes(i,2) = Tend - Tstart
	if (HWorker) then
	  error = maxval(abs(
     &		gBDhy(:,owned_hslices(myid,1:no_of_hslices(myid)),:)
     &   	-sBDhy(:,1:no_of_hslices(myid),:)))
          if (error .gt. 0.0) then
	    print*,'P',myid, 'After single msg BDy h-distribution',
     &		' max error = ', error, ' check implementation '
          end if
        end if
	sBDhy(:,:,:) = 0
      end do

      if (Master) print*, 'Collective BDy h-distribution test starts'
      do i=1, Ntests
        call MPI_BARRIER(MPI_COMM_WORLD, ierr)
	Tstart = MPI_WTIME()
        call distrib_h_BDy_v3(ix,iz,izloc,is,gBDhy,sBDhy)
	if (HWorker) Tend = MPI_WTIME()
	call MPI_BARRIER(MPI_COMM_WORLD, ierr)
	if (Master) Tend = MPI_WTIME()
	testTimes(i,3) = Tend - Tstart
	if (HWorker) then
	  error = maxval(abs(
     &		gBDhy(:,owned_hslices(myid,1:no_of_hslices(myid)),:)
     &   	-sBDhy(:,1:no_of_hslices(myid),:)))
          if (error .gt. 0.0) then
	    print*,'P',myid, 'After collective BDy h-distribution',
     &		' max error = ', error, ' check implementation '
          end if
        end if
	sBDhy(:,:,:) = 0
      end do
      
      if (Master) then
        print*, 'Global standard BDy h-distribution time:',
     &		minval(testTimes(1:Ntests,1))
        print*, 'Global single msg BDy h-distribution time:',
     &		minval(testTimes(1:Ntests,2))
        print*, 'Global collective BDy h-distribution time:',
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
        do i=1,NHworkers+1
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
	  print*, 'Decision: standard BDy h-distribution.'
c	  WRITE (10,*) 'distr_h_BDy= ',1
	else if (best .eq. dcount(2)) then
	  print*, 'Decision: single msg BDy h-distribution.'
c	  WRITE (10,*) 'distr_h_BDy= ',2
	else 
	  print*, 'Decision: collective BD h-distribution.'
c	  WRITE (10,*) 'distr_h_BDy= ',3
	end if
      end if
	
c-------------------------------------------------------------------------------
c     Test BD hx-distribution     
c-------------------------------------------------------------------------------

      if (Master) print*, 'Standard BD hx-distribution test starts'
      do i=1, Ntests
        call MPI_BARRIER(MPI_COMM_WORLD, ierr)
	Tstart = MPI_WTIME()
        call distrib_xh_BD_v1(ix,iz,izloc,is, gBDhx2,sBDhx2)
	if (HWorker) Tend = MPI_WTIME()
	call MPI_BARRIER(MPI_COMM_WORLD, ierr)
	if (Master) Tend = MPI_WTIME()
	testTimes(i,1) = Tend - Tstart
	if (HWorker) then
	  error = maxval(abs(
     &		gBDhx2(:,owned_hslices(myid,1:no_of_hslices(myid)),:,:)
     &   	-sBDhx2(:,1:no_of_hslices(myid),:,:)))
          if (error .gt. 0.0) then
	    print*,'P',myid, 'After standard BD hx-distribution',
     &		' max error = ', error, ' check implementation '
          end if
        end if
	sBDhx2(:,:,:,:) = 0
      end do
      

      if (Master) print*, 'Single msg BD hx-distribution test starts'
      do i=1, Ntests
        call MPI_BARRIER(MPI_COMM_WORLD, ierr)
	Tstart = MPI_WTIME()
        call distrib_xh_BD_v2(ix,iz,izloc,is, gBDhx2,sBDhx2)
	if (HWorker) Tend = MPI_WTIME()
	call MPI_BARRIER(MPI_COMM_WORLD, ierr)
	if (Master) Tend = MPI_WTIME()
	testTimes(i,2) = Tend - Tstart
	if (HWorker) then
	  error = maxval(abs(
     &		gBDhx2(:,owned_hslices(myid,1:no_of_hslices(myid)),:,:)
     &   	-sBDhx2(:,1:no_of_hslices(myid),:,:)))
          if (error .gt. 0.0) then
	    print*,'P',myid, 'After single msg BD hx-distribution',
     &		' max error = ', error, ' check implementation '
          end if
        end if
	sBDhx2(:,:,:,:) = 0
      end do

      if (Master) print*, 'Collective BD hx-distribution test starts'
      do i=1, Ntests
        call MPI_BARRIER(MPI_COMM_WORLD, ierr)
	Tstart = MPI_WTIME()
        call distrib_xh_BD_v3(ix,iz,izloc,is, gBDhx2,sBDhx2)
	if (HWorker) Tend = MPI_WTIME()
	call MPI_BARRIER(MPI_COMM_WORLD, ierr)
	if (Master) Tend = MPI_WTIME()
	testTimes(i,3) = Tend - Tstart
	if (HWorker) then
	  error = maxval(abs(
     &		gBDhx2(:,owned_hslices(myid,1:no_of_hslices(myid)),:,:)
     &   	-sBDhx2(:,1:no_of_hslices(myid),:,:)))
          if (error .gt. 0.0) then
	    print*,'P',myid, 'After collective BD hx-distribution',
     &		' max error = ', error, ' check implementation '
          end if
        end if
	sBDhx2(:,:,:,:) = 0
      end do
      
      if (Master) then
        print*, 'Global standard BD hx-distribution time:',
     &		minval(testTimes(1:Ntests,1))
        print*, 'Global single msg BD hx-distribution time:',
     &		minval(testTimes(1:Ntests,2))
        print*, 'Global collective BD hx-distribution time:',
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
        do i=1,NHworkers+1
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
	  print*, 'Decision: standard BD hx-distribution.'
	  WRITE (10,*) 'distr_xh_BD= ',1
	else if (best .eq. dcount(2)) then
	  print*, 'Decision: single msg BD hx-distribution.'
	  WRITE (10,*) 'distr_xh_BD= ',2
	else 
	  print*, 'Decision: collective BD hx-distribution.'
	  WRITE (10,*) 'distr_xh_BD= ',3
	end if
      end if
      
c-------------------------------------------------------------------------------
c     Test BD hy-distribution     
c-------------------------------------------------------------------------------
      if (Master) print*, 'Standard BD hy-distribution test starts'
      do i=1, Ntests
        call MPI_BARRIER(MPI_COMM_WORLD, ierr)
	Tstart = MPI_WTIME()
        call distrib_yh_BD_v1(iy,iz,izloc,is,gBDhy2,sBDhy2)
	if (HWorker) Tend = MPI_WTIME()
	call MPI_BARRIER(MPI_COMM_WORLD, ierr)
	if (Master) Tend = MPI_WTIME()
	testTimes(i,1) = Tend - Tstart
	if (HWorker) then
	  error = maxval(abs(
     &		gBDhy2(:,owned_hslices(myid,1:no_of_hslices(myid)),:,:)
     &   	-sBDhy2(:,1:no_of_hslices(myid),:,:)))
          if (error .gt. 0.0) then
	    print*, 'P',myid,'After standard BD hy-distribution',
     &		' max error = ', error, ' check implementation '
          end if
        end if
	sBDhy2(:,:,:,:) = 0
      end do
      

      if (Master) print*, 'Single msg BD hy-distribution test starts'
      do i=1, Ntests
        call MPI_BARRIER(MPI_COMM_WORLD, ierr)
	Tstart = MPI_WTIME()
        call distrib_yh_BD_v2(iy,iz,izloc,is,gBDhy2,sBDhy2)
	if (HWorker) Tend = MPI_WTIME()
	call MPI_BARRIER(MPI_COMM_WORLD, ierr)
	if (Master) Tend = MPI_WTIME()
	testTimes(i,2) = Tend - Tstart
	if (HWorker) then
	  error = maxval(abs(
     &		gBDhy2(:,owned_hslices(myid,1:no_of_hslices(myid)),:,:)
     &   	-sBDhy2(:,1:no_of_hslices(myid),:,:)))
          if (error .gt. 0.0) then
	    print*,'P',myid, 'After single msg BD hy-distribution',
     &		' max error = ', error, ' check implementation '
          end if
        end if
	sBDhy2(:,:,:,:) = 0
      end do

      if (Master) print*, 'Collective BD hy-distribution test starts'
      do i=1, Ntests
        call MPI_BARRIER(MPI_COMM_WORLD, ierr)
	Tstart = MPI_WTIME()
        call distrib_yh_BD_v3(iy,iz,izloc,is,gBDhy2,sBDhy2)
	if (HWorker) Tend = MPI_WTIME()
	call MPI_BARRIER(MPI_COMM_WORLD, ierr)
	if (Master) Tend = MPI_WTIME()
	testTimes(i,3) = Tend - Tstart
	if (HWorker) then
	  error = maxval(abs(
     &		gBDhy2(:,owned_hslices(myid,1:no_of_hslices(myid)),:,:)
     &   	-sBDhy2(:,1:no_of_hslices(myid),:,:)))
          if (error .gt. 0.0) then
	    print*,'P',myid, 'After collective BD hy-distribution',
     &		' max error = ', error, ' check implementation '
          end if
        end if
	sBDhy2(:,:,:,:) = 0
      end do
      
      if (Master) then
        print*, 'Global standard BD hy-distribution time:',
     &		minval(testTimes(1:Ntests,1))
        print*, 'Global single msg BD hy-distribution time:',
     &		minval(testTimes(1:Ntests,2))
        print*, 'Global collective BD hy-distribution time:',
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
        do i=1,NHworkers+1
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
	  print*, 'Decision: standard BD hy-distribution.'
	  WRITE (10,*) 'distr_yh_BD= ',1
	else if (best .eq. dcount(2)) then
	  print*, 'Decision: single msg BD hy-distribution.'
	  WRITE (10,*) 'distr_yh_BD= ',2
	else 
	  print*, 'Decision: collective BD hy-distribution.'
	  WRITE (10,*) 'distr_yh_BD= ',3
	end if
	
	CLOSE (10)
	
      end if

	           
      call MPI_FINALIZE(rc)
      stop
      end
