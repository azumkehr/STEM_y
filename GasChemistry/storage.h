      integer nsmax,npmax
c     nsmax = maximum nr. of steps inside a time period dt
c     npmax = maximum number of time periods of length dt 
      parameter (nsmax = 500, npmax = 1000)
      integer iperiod,istep,istore
c     iperiod = counter for number of time periods of length dt
c     istep = counter for nr of steps inside one time period dt      
c     istore = flag for trajectory storage inside one
c              time period : istore = 1 => store 

      real*8 fwdtrajp(npmax,nvar)
c     fwdtrajp stores forward trajectory after each time period dt      
      real*8 stepvect(nsmax)
c     stepvect stores the number of steps inside one time period dt      
      real*8 fwdtraji(nsmax,nvar)
c     fwdtraji stores forward trajectory inside one time period dt          
      common /storage/ stepvect,fwdtrajp,fwdtraji,
     &                 iperiod,istep,istore
