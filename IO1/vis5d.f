!------------------------------------------------
! Module Name: vis5d 
! Description: Describes vis5d variables 
!
! Developed on: 10/16/02
!
! Usage Instructions: use aqmax
!
! Author: Adrian Sandu
!         C Belwal
!-----------------------------------------------
	module vis5d

	  integer::TOTAL_VARIABLES
	  integer::imax,kmax,jmax
	  real, dimension(:,:,:),pointer :: gg

	   integer begdate, begtime, tstep, numtimes         ! beginning date in YYYYDDD


	  !NOTE: These values have been taking from the v5df.h file
	  !      The file does not get included here 
	  parameter(MAXV=100,MAXT=400,MAXL=100,IMISS=-987654,MISS=1.0E35)
	  parameter(ipoff=10,jpoff=25,ipmax=100,jpmax=90,lj=11)   ! nesting ratio
	  integer::lmet,lchem

	  real, dimension(:,:), pointer::xlon,xlat,topo
          real, dimension(:,:,:), pointer::zheight
	  !---- common vars
	  real::paiv
	  real::ddx,dllat,dllon,xxstart,yystart,plat,plon
	  real::wlon,elon,slat,nlat,nllat,nllon
	  real::vglvs(20) !Change to Nz
	  integer::kout


	  integer nr, nc, nl(MAXV)
          integer numvars
          character*10 varname(MAXV)
          integer dates(MAXT)
          integer times(MAXT)
          integer compressmode
          integer projection
          real proj_args(100)
          integer vertical
          real vert_args(MAXL)


      data nr,nc / IMISS, IMISS /
      data (nl(i),i=1,MAXV) / MAXV*IMISS /
      data numtimes,numvars / IMISS, IMISS /
      data (varname(i),i=1,MAXV) / MAXV*"      " /         
      data (dates(i),i=1,MAXT) / MAXT*IMISS /
      data (times(i),i=1,MAXT) / MAXT*IMISS /
      data projection / IMISS /
      data (proj_args(i),i=1,100) / 100*MISS /
      data vertical / IMISS /
      data (vert_args(i),i=1,MAXL) / MAXL*MISS /


	end module
