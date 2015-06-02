53c53
<       !s1(1:ix,j,k,nstart:nend) = Conc(1:ix,1:num)
---
> 
104d103
<       !s1(i,1:iy,k,nstart:nend) =Conc(1:iy,1:num)
153a153
>       
164c164
< 
---
>              
166a167
> 			 			          
169d169
<       !s1(i,j,1:iz,nstart:nend) =Conc(1:iz,1:num)
232d231
<       !s1(i,j,1:iz,nstart:nend) = Conc(1:iz,1:num)
281c280
<       Sb(1:2,1:num) = 0!sx(j,k,1:2,nstart:nend)
---
>       Sb(1:2,1:num) = sx(j,k,1:2,nstart:nend)
335c334
<       Sb(1:2,1:num) = 0!sy(i,k,1:2,nstart:nend)
---
>       Sb(1:2,1:num) = sy(i,k,1:2,nstart:nend)
978c977
<   call ADVDIFF_JAC_FDZ_MF(N,X,U,K,Air,Jac)  
---
>   call ADVDIFF_JAC_FDZ_MF2(N,X,U,K,Air,Jac)  
989,990c988,989
<           call DGBTRS('T', N, kl, ku, 1, A, lda, IPIV, &
<           	       Lam(1,ispec), N, INFO )
---
>           !call DGBTRS('T', N, kl, ku, 1, A, lda, IPIV, &
>           !	       Lam(1,ispec), N, INFO )
994a994,995
>           call DGBTRS('T', N, kl, ku, 1, A, lda, IPIV, &
>                        Lam(1,ispec), N, INFO )
1100c1101
<     Jac(row(1,1),1)  =  & !-W(1)/DZ(1)       &
---
>     Jac(row(1,1),1)  =  -W(1)/DZ(1)       &
1136c1137
< 	     )/( Air(i)*2*DZ(N-1) )
---
> 	     )/( Air(N)*2*DZ(N-1) )
1141c1142
< 	     )/( Air(i)*2*DZ(N-1) )
---
> 	     )/( Air(N)*2*DZ(N-1) )
1146c1147
< 	 )/( Air(i)*2*DZ(N) )  
---
> 	 )/( Air(N)*2*DZ(N) )  
1150c1151
< 	 )/( Air(i)*2*DZ(N) )
---
> 	 )/( Air(N)*2*DZ(N) )
1167a1169,1254
> subroutine ADVDIFF_JAC_FDZ_MF2(N,Z,W,K,Air,Jac)
> !
> implicit none
> !
> integer, intent(in) :: N
> double precision, intent(in)  :: Z(N), W(N), K(N), Air(N)
> ! Jacobian for time derivative of the concentration
> integer, parameter :: kl=1, ku=1
> double precision, intent(out) :: Jac(kl+ku+1,N)
> 
> ! difflux/advflux = diffusive/advective fluxes through i-1/2
> integer :: i
> double precision :: difflux
> double precision :: AK(N), DZ(N)
> 
> AK = Air*K
> DZ(1:N-1) = Z(2:N)-Z(1:N-1); DZ(N)=DZ(N-1)
> 
> Jac = 0.d0
> 
> if ( W(1)>0.d0 ) then
>     Jac(row(1,1),1)  =  -W(1)/DZ(1)       &
>        -(AK(2)+AK(1))/2/DZ(1)**2/Air(1) 
>     Jac(row(1,2),2)  =        &
>        + (AK(2)+AK(1))/2/DZ(1)**2/Air(1) 
> else
>     Jac(row(1,1),1)  =  W(1)/DZ(1)  &
>        - (AK(2)+AK(1))/2/DZ(1)**2/Air(1)     
>     Jac(row(1,2),2)  = -W(1)/DZ(1)          &
>        + (AK(2)+AK(1))/2/DZ(1)**2/Air(1)     
> end if  
> 
> ! Intermediate Boundaries
> do i=2,N-1
>   if (W(i)>=0) then
>     Jac(row(i,i),i)     = -W(i)/DZ(i-1)
>     Jac(row(i,i-1),i-1) =  W(i)/DZ(i-1)
>   else
>     Jac(row(i,i),i)     =  W(i)/DZ(i)
>     Jac(row(i,i+1),i+1) = -W(i)/DZ(i)
>   end if  
>   Jac(row(i,i-1),i-1) = Jac(row(i,i-1),i-1) +    &
>          ( (AK(i)+AK(i-1))/DZ(i-1)               &
> 	 )/( Air(i)*(Z(i+1)-Z(i-1)) )
>   Jac(row(i,i),i)     = Jac(row(i,i),i) -        &
>          ( (AK(i)+AK(i-1))/DZ(i-1)               &
>           +(AK(i+1)+AK(i))/DZ(i)                 &
> 	 )/( Air(i)*(Z(i+1)-Z(i-1)) )
>   Jac(row(i,i+1),i+1) = Jac(row(i,i+1),i+1) +    &
>          ( (AK(i+1)+AK(i))/DZ(i)                 &
> 	 )/( Air(i)*(Z(i+1)-Z(i-1)) )
> end do
> 
> ! Top of the domain
> if ( W(N)<0 ) then ! inflow
>   Jac(row(N,N-1),N-1)    = Jac(row(N,N-1),N-1) + &
>              ( (AK(N)+AK(N-1))/DZ(N-1)           &
> 	     )/( Air(N)*2*DZ(N-1) )
>   Jac(row(N,N),N)     = Jac(row(N,N),N) +        &
>                W(N)/DZ(N-1) +                    &
>              ( -2*AK(N)/DZ(N-1)                  &
> 	      -(AK(N)+AK(N-1))/DZ(N-1)           &
> 	     )/( Air(N)*2*DZ(N-1) )
> else ! outflow
>   Jac(row(N,N-1),N-1)    = Jac(row(N,N-1),N-1) + &
>            W(N)/DZ(N) +                          &
>          ( (AK(N)+AK(N-1))/DZ(N)                 &
> 	 )/( Air(N)*2*DZ(N) )  
>   Jac(row(N,N),N)     = Jac(row(N,N),N) +  &
>          -W(N)/(Z(N)-Z(N-1)) +             &
>          ( -(AK(N)+AK(N-1))/DZ(N)          &
> 	 )/( Air(N)*2*DZ(N) )
> end if      
> 
> contains
>   
>   integer function row(i,j)
>   ! gives the row of the Blas banded format for pentadiagonal Jacobian
>   integer :: i, j
>   integer, parameter :: kl=1, ku=1
>   if ( (i<=0) .or. (j<=0) ) then
>      print*,'Error in ADVDIFF_JAC_FDZ_MF. i,j=',i,j
>      stop
>   end if
>   row = ku + 1 + i - j
>   end function row
1168a1256
> end subroutine ADVDIFF_JAC_FDZ_MF2
1195c1283
<      advflux = 0.d0 !W(1)*Bdry(1)
---
>      advflux = W(1)*Bdry(1)
