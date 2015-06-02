      ! Name of the executable model
      module execname
        character(len=*), parameter :: 
     &         exec_name = "./icartt.run "
      end module execname

      program aq_driver
      
c--------------------------------------------------------------------------
      include 'aqmax.param'
      include 'aqms.param'
      
      interface
      subroutine setulb(n, m, x, l, u, nbd, f, g, factr, pgtol, wa, iwa,
     +                 task, iprint, csave, lsave, isave, dsave)
 
      character*60     task, csave
      logical          lsave(4)
      integer          n, m, iprint, 
     +                 nbd(n), iwa(3*n), isave(44)
      double precision f, factr, pgtol, x(n), l(n), u(n), g(n),
     +                 wa(2*m*n+4*n+11*m*m+8*m), dsave(29)
      end subroutine setulb
      end interface

c--------------------------------------------------------------------------
      integer :: ix=ixm, iy=iym, iz=izm  ! Grid Dimensions
      integer :: N_gas=iLm               ! No. of gas species
c--------------------------------------------------------------------------
c The control Variables
      integer, parameter :: NUM_CTRL=1 !Number of species (at t=0) to adjust
      integer, parameter :: ICTRL(NUM_CTRL)=(/1/)
      real, parameter :: small=1e-6
      !The choice of ctrls are based upon the gradient magnitudes 4/20/04
c--------------------------------------------------------------------------
      real, dimension(:,:,:,:), pointer :: emi_fac,emi_grd,emi_unc
      character(len=32) :: fname_conc_ini, fname_lambda
      integer :: unit_emi = 78, unit_emi_grd = 79
      integer :: unit_howto=80, unit_bckg=81, unit_opt=82
      integer :: record_conc_ini = 0
      integer :: urpt = 188, uemiunc=189 
c-------------------------------------------------------------------------------
c  optimization variables and parameters
c-------------------------------------------------------------------------------
      double precision :: costfct,costfct_ini,bfgs_magn
      double precision :: cost_mismatch 
      double precision, dimension(:), pointer :: xoptim, gradient,
     &             lbfgs_l,  lbfgs_u, lbfgs_wa, xoptim_bck,xuncert
      integer, dimension(:), pointer :: lbfgs_nbd, lbfgs_iwa
      integer, dimension(:) :: i_map,j_map
      integer :: lbfgs_nmax, lbfgs_mmax, lbfgs_lenwa
      parameter (lbfgs_nmax=5000000,  lbfgs_mmax=17)
      parameter (lbfgs_lenwa=2*lbfgs_mmax*lbfgs_nmax +4*lbfgs_nmax 
     &           +11*lbfgs_mmax*lbfgs_mmax+8*lbfgs_mmax)
      character lbfgs_task*60,lbfgs_csave*60
      integer :: opt_iteration, model_runs,model_pre=0,model_max=25
      logical :: lbfgs_lsave(4),restart
      integer :: lbfgs_n,lbfgs_m,lbfgs_iprint,lbfgs_isave(44), Num_grid
      double precision :: lbfgs_factor,lbfgs_pgtol,lbfgs_dsave(29)

C ---- Open report file
       open(unit=urpt,file='Report.opt')
       !open(unit=uemiunc,file='Hg90x60_emi_unc.dat')
C ----- Allocate Memory -------------
       print*,"OPT: Allocate Memory"
       allocate( emi_fac(ix,iy,2,1),STAT=ierr )
       allocate( emi_grd(ix,iy,2,1),STAT=ierr )
       allocate( emi_unc(ix,iy,2,1),STAT=ierr )
C ---- Give the upper limit, to take advantage of lbfgs-b ------ 
       emi_fac=1.0
      
       allocate(i_map(ix*iy),STAT=ierr )
       allocate(j_map(ix*iy),STAT=ierr ) 
       lbfgs_n=0 
       do j=1,iy
       do i=1,ix
          !read(uemiunc,*) itemp,jtemp,unc_temp
          if(.true.) then
            itemp=i
	    jtemp=j
	    unc_temp=0.5
          endif 
          if((i-itemp)*(j-jtemp)>small) then
            write(*,*) 'reading error for uncertainty'
            stop
          endif
          emi_unc(i,j,1:2,1) = unc_temp 
          if(emi_unc(i,j)>small) then
            lbfgs_n=lbfgs_n+1 
            i_map(lbfgs_n)=i
            j_map(lbfgs_n)=j
          endif 
       enddo
       enddo
       lbfgs_n=lbfgs_n*2  ! Two sets of scaling factors, 1 for surface, one for upper
        
       do i=1,lbfgs_n/2
          print *, i, i_map(i), j_map(i), emi_unc(i_map(i),j_map(i))
       enddo 
c----------------------------------------------------------------------c
c                  Master Performs Optimization                        c           
c----------------------------------------------------------------------c
      print*,"OPT: Start Optimization Loop"

         bfgs_magn=1d0     ! 1d18 factor is used before adding uncertainty term in function       
         lbfgs_iprint=1
         lbfgs_factor=1.0d12
         lbfgs_pgtol=1.0d-16
         !lbfgs_n=ix*iy*1          ! dimension of the state vector
         lbfgs_m = 5               ! number of pairs used by the limited memory l-bfgs        
	 allocate( xoptim(lbfgs_n),  STAT=ierr)
         allocate( xuncert(lbfgs_n),  STAT=ierr)
         allocate( xoptim_bck(lbfgs_n),  STAT=ierr)
	 allocate( lbfgs_l(lbfgs_n),  STAT=ierr)
	 allocate( lbfgs_u(lbfgs_n),  STAT=ierr)
	 allocate( lbfgs_nbd(lbfgs_n),  STAT=ierr)
         lbfgs_nbd(1:lbfgs_n) = 2    !  Both lower and upper bounds
	 allocate( gradient(lbfgs_n), STAT=ierr)
	 allocate( lbfgs_wa((2*lbfgs_m + 4)*lbfgs_n + 11*lbfgs_m**2 + 
     &             8*lbfgs_m), STAT=ierr)
	 allocate( lbfgs_iwa(3*lbfgs_n), STAT=ierr)

       lbfgs_u=10.0
       lbfgs_l=0.1	 
       Num_grid=lbfgs_n

C----  Background /initial guess  =1
       do i=1,lbfgs_n/2
          xoptim(i)=emi_fac(i_map(i),j_map(i),1,1)     
          xuncert(i)=emi_unc(i_map(i),j_map(i),1,1)  
          xoptim(lbfgs_n/2+i)=emi_fac(i_map(i),j_map(i),2,1)
          xuncert(lbfgs_n/2+i)=emi_unc(i_map(i),j_map(i),2,1)
       enddo

c +++++  Comment this section out when doing emission inversion ++++++
cc   Just do the forward
c
       do i=1,lbfgs_n/2
          emi_fac(i_map(i),j_map(i),1,1)=xoptim(i)
          emi_fac(i_map(i),j_map(i),2,1)=xoptim(i+lbfgs_n/2)
       enddo
c
        call simulation('fwd',ix,iy,iz,N_gas,emi_fac,emi_grd,costfct) ! this can be fbw or fwd
c
        print*, 'Force stop in aq_driver_opt.f, just want forward'
        stop ! stop here if all want is forward
c ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


c +++++  Comment this section out when doing emission inversion ++++++
cc   Compute the forward model and cost function for getting emi_grd
c
c       do i=1,lbfgs_n/2
c          emi_fac(i_map(i),j_map(i),1,1)=xoptim(i)
c          emi_fac(i_map(i),j_map(i),2,1)=xoptim(i+lbfgs_n/2)
c       enddo
c
c         call simulation('fbw',ix,iy,iz,N_gas,emi_fac,emi_grd,costfct) ! this can be fbw or fwd
c       print*,'DEBUG emi_grd ',maxval(emi_grd)
c       print*,'DEBUG emi_grd ',maxval(emi_fac)
c        open(72,file='emi_grd.dat')
c        do isen = 1, ixm
c          do jsen = 1,iym
c                write(72,*) isen,jsen,emi_grd(isen,jsen,1,1)!isen,jsen,emi_grd(isen,jsen,2,1)
c          enddo
c        enddo
c        close(72)
c
c        print*, 'Force stop in aq_driver_opt.f, just want sensitivity'
c        stop ! stop here if all want is gradient for sensitivity
c ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


      print*,'__________START OPTIMIZATION__________'   
         model_runs=0
         lbfgs_task = 'START'
c     ------- The beginning of the loop ----------
  111    continue

         if(model_runs .ge. model_max  ) goto 777
 
       print*,"OPT: Run L-BFGS-B iteration no ", model_runs
       
       costfct=costfct*bfgs_magn
       gradient=gradient*bfgs_magn 
         call setulb(lbfgs_n,lbfgs_m,xoptim,lbfgs_l,lbfgs_u,
     &	 lbfgs_nbd, costfct,gradient,lbfgs_factor,lbfgs_pgtol,
     &   lbfgs_wa,lbfgs_iwa,
     &   lbfgs_task,lbfgs_iprint,lbfgs_csave,lbfgs_lsave,
     &   lbfgs_isave,lbfgs_dsave)
       costfct=costfct/bfgs_magn 
       gradient=gradient/bfgs_magn      
   
       print*,"     LBFGS-TASK = ",lbfgs_task
       write(urpt,*) "     LBFGS-TASK = ",lbfgs_task
  
      if (lbfgs_task(1:2) .eq. 'FG') then
 
c        The minimization routine has returned to request the
c        function f and gradient g values at the current x.

c   Compute the forward model and cost function
       
       do i=1,lbfgs_n/2
          emi_fac(i_map(i),j_map(i),1,1)=xoptim(i)
          emi_fac(i_map(i),j_map(i),2,1)=xoptim(i+lbfgs_n/2)     
       enddo   


         call simulation('fbw',ix,iy,iz,N_gas,emi_fac,emi_grd,costfct) 
         cost_mismatch=costfct
         model_runs=model_runs+1
         
	 if(model_runs == 1) then 
            costfct_ini = costfct 
         endif
          
         print *,'Mismatch part of the cost, costfct=', costfct  !molefraction^2
         costfct=costfct+SUM(((xoptim-1)/xuncert)**2)/2.0 
         print *,'After adding background error, costfct=', costfct
         !============================================================== 
         do i=1,lbfgs_n/2
          gradient(i)=emi_grd(i_map(i),j_map(i),1,1) 
     &      +(emi_fac(i_map(i),j_map(i),1,1)-1.)/xuncert(i)**2
          gradient(lbfgs_n/2+i)=emi_grd(i_map(i),j_map(i),2,1) 
     &      +(emi_fac(i_map(i),j_map(i),2,1)-1.)/xuncert(lbfgs_n/2+i)**2 
         enddo

         call print_report(NUM_CTRL,lbfgs_n,ICTRL,urpt,opt_iteration,
     &           model_runs, lbfgs_task,costfct,cost_mismatch,gradient)

         if (costfct .le. 0.001*costfct_ini) then
           write(urpt,*) 'Exit on fcn condition'
           print*, 'Exit on fcn condition'
           goto 777
         end if
         goto 111
   
        elseif (lbfgs_task(1:5) .eq. 'NEW_X') then
          do i=1,lbfgs_n/2
            emi_fac(i_map(i),j_map(i),1,1)=xoptim(i)
            emi_fac(i_map(i),j_map(i),2,1)=xoptim(i+lbfgs_n/2)
          enddo
         call print_report(NUM_CTRL,lbfgs_n,ICTRL,urpt,opt_iteration,
     &           model_runs, lbfgs_task,costfct,cost_mismatch,gradient)

         goto 111
        
        else

         call print_report(NUM_CTRL,lbfgs_n, ICTRL, urpt,opt_iteration,
     &    model_runs, lbfgs_task,costfct,cost_mismatch,gradient)

         if (lbfgs_iprint.le.-1 .and. lbfgs_task(1:4).ne.'STOP') then
	    print*, lbfgs_task
	 end if   
        endif
c     ---------- The end of the loop -------------      
777   continue
      open(124,file='q_fac.plt')
      write(124,*) 'ZONE F=point, I=',ix,' J=',iy
            do iit=1,1
            do j=1,iy
            do i=1,ix
                write(124,*) i,j,emi_fac(i,j,1,iit)
            enddo
            enddo
            enddo
      close(124)
      open(124,file='emi_fac.plt')
      write(124,*) 'ZONE F=point, I=',ix,' J=',iy
            do iit=1,1
            do j=1,iy
            do i=1,ix
                write(124,*) i,j,emi_fac(i,j,2,iit)
            enddo
            enddo
            enddo
      close(124)

       print*,'________END OPTIMIZATION____________'
       close(urpt)
	
       end program aq_driver


c--------------------------------------------------------------------------
      subroutine print_report(NUM_CTRL, Num_grid, ICTRL, urpt,
     &                        opt_iteration, model_runs, lbfgs_task,
     &                        costfct, cost_mismatch, gradient)
      integer :: NUM_CTRL, Num_grid, ICTRL(NUM_CTRL), urpt, 
     &           opt_iteration, model_runs,i_ctrl
      character*60 :: lbfgs_task
      double precision :: costfct,cost_mismatch, 
     &                    gradient(NUM_CTRL*Num_grid)
	 
      write(urpt,*) '***** ',lbfgs_task,' FOLLOWS ******'
      write(urpt,*) 'It=',opt_iteration,'Num of model runs', model_runs
      write(urpt,*) 'Cost=',costfct
      write(urpt,*) 'Misfit & bckg:',cost_mismatch,costfct-cost_mismatch
      write(urpt,*) 'Gradient=',
     &                sqrt(SUM(gradient**2)/dble(NUM_CTRL*Num_grid)) 
      do i_ctrl=1, NUM_CTRL
         write(urpt,*) 'CONTROL index:', ICTRL(i_ctrl)
         write(urpt,*) 'MAX abs GRAD:', 
     &    MAXVAL(abs(gradient((i_ctrl-1)*Num_grid+1:i_ctrl*Num_grid)))
         write(urpt,*) 'RMS of GRAD',
     &    sqrt(SUM(gradient((i_ctrl-1)*Num_grid+1:i_ctrl*Num_grid)**2)/
     &         dble(Num_grid))      
      enddo
      write(urpt,*) 
     
      end subroutine print_report
c--------------------------------------------------------------------------

c--------------------------------------------------------------------------
c  This subroutine calls STEM with:
c           - the initial concentration field sg1
c           - the final adjoint field Lambda
c  Upon return:
c           - sg1 holds the final concentration field
c           - Lambda holds the adjoint vars at initial time
c
c  Note: if mode == 'ini' then the model only reads the initial
c           values of sg1 and returns
c
c  mode   = 'fwd' : forward only run, to calculate cost function
c           'fbw' : both forward and backward runs (recommended)
c           'obs' : forward run to generate the observation file
c           'ini' : for initialization only
c
c  Purpose : it is a clean interface to the simulation code
c  Note:     it communicates with stem via 2 direct access files,
c            one for the conc variables and one for the adjoint   
c--------------------------------------------------------------------------
        subroutine simulation(mode,ix,iy,iz,N_gas,
     &                        emi_fac,emi_grd,costfct)
        use execname

C ----- Declare Variables -------------
        character(len=3) :: mode
	integer :: ix, iy, iz, N_gas
	real :: emi_fac(ix,iy,2,1), emi_grd(ix,iy,2,1)

C ----- Declare The file for initial concentrations -------------
        character(len=32) :: fname_emi_fac, fname_emi_grd
        integer :: unit_emi_fac = 78, unit_emi_grd = 79
	integer :: unit_mode = 80,unit_cost=81
        integer :: record_emi_fac = 1, record_emi_grd= 1
        double precision :: costfct
 
C --- Check if the value of mode is appropriate
        if ( (mode(1:3).ne.'fwd').and.(mode(1:3).ne.'fbw')
     &        .and.(mode(1:3).ne.'ini').and.(mode(1:3).ne.'obs') ) then
             print*,'Error in simulation: mode = <',mode,'>'
	     print*,'Accepted values are <fwd>, <fbw>, <ini> or <obs>'
	end if     		 

C ---- Open File for Mode -----
        open( unit=unit_mode, file='TmpMode' )
        write( unit_mode, fmt="(A3)" ) mode(1:3)
        close( unit_mode )

C ---- Direct Access Files -----
        fname_emi_fac = 'TmpEmiFac'
        fname_emi_grd = 'TmpEmiGrd'
        
C ---- Write Initial Conc. 
	if( mode=='fwd' .or. mode=='fbw'.or.mode=='obs') then
           open(unit_emi_fac,file=fname_emi_fac, access='direct',
     &         recl=4*ix*iy*2*1)
            write(unit_emi_fac) emi_fac
           close( unit_emi_fac)
        end if

 
C ---- Call simulation by executing Adjoint stem in a new shell -------------
        call system( exec_name )
        print *, 'after exec_name call'

       costfct=-1d0 ! initialize with a negative number
       if ( mode == 'fbw'.or.mode == 'fwd') then
        open(unit=unit_cost, file='costfct')
        read(unit_cost,*) costfct 
        close(unit_cost)
        call system('rm costfct') 
       endif 
        
C ---- Read The results -----

        if ( mode == 'fbw' ) then
           open(unit_emi_grd,file=fname_emi_grd, access='direct',
     &	       recl=4*ix*iy*2*1)
           read(unit_emi_grd) emi_grd 
           close(unit_emi_grd)
        end if
 
C ---- Clean the Temporary Files ----
        call system('rm TmpMode')
        call system('rm caca')
          
        end subroutine simulation

