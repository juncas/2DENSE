subroutine boundary_rt_general
use flow_general
use var_index
implicit none
double precision :: rho,u,v,p,c
integer :: index_i,index_j,index_var,index_eq
integer :: bfsize_i,bfsize_j,ni,nj

ni = block_grid(num_block)%ni
nj = block_grid(num_block)%nj
bfsize_i=block_grid(num_block)%bfsize_i
bfsize_j=block_grid(num_block)%bfsize_j

! Bottom 
do index_j=1,bfsize_j
    do index_i=bfsize_i+1,ni+bfsize_i
        
        rho=2.d0
        u=0.d0
        v=0.d0
        p=1.d0

        block_grid(num_block)%cons(index_i,index_j,1) = rho
        block_grid(num_block)%cons(index_i,index_j,2) = rho*u
        block_grid(num_block)%cons(index_i,index_j,3) = rho*v
        block_grid(num_block)%cons(index_i,index_j,4) = p/(gamma-1.d0)+0.5d0*rho*(u**2+v**2)
        
        block_grid(num_block)%prim(index_i,index_j,1) = rho
        block_grid(num_block)%prim(index_i,index_j,2) = u
        block_grid(num_block)%prim(index_i,index_j,3) = v
        block_grid(num_block)%prim(index_i,index_j,4) = p
        block_grid(num_block)%prim(index_i,index_j,5) = sqrt(gamma*p/rho)
        
    enddo
enddo

! Top 
do index_j=nj+bfsize_j+1,nj+2*bfsize_j
    do index_i=bfsize_i+1,ni+bfsize_i
        
        rho=1.d0
        p=2.5d0
        u=0.d0
        v=0.d0
        block_grid(num_block)%cons(index_i,index_j,1) = rho
        block_grid(num_block)%cons(index_i,index_j,2) = rho*u
        block_grid(num_block)%cons(index_i,index_j,3) = rho*v
        block_grid(num_block)%cons(index_i,index_j,4) = p/(gamma-1.d0)+0.5d0*rho*(u**2+v**2)
        
        block_grid(num_block)%prim(index_i,index_j,1) = rho
        block_grid(num_block)%prim(index_i,index_j,2) = u
        block_grid(num_block)%prim(index_i,index_j,3) = v
        block_grid(num_block)%prim(index_i,index_j,4) = p
        block_grid(num_block)%prim(index_i,index_j,5) = sqrt(gamma*p/rho)
        
    enddo
enddo

! !left
! do index_j=bfsize_j+1,nj+bfsize_j
!     do index_i=1,bfsize_i
!     !Cons
!       block_grid(num_block)%cons(index_i,index_j,1) &
!     = block_grid(num_block)%cons(2*bfsize_i+index_i-1,index_j,1)

!       block_grid(num_block)%cons(index_i,index_j,2) &
!     =-block_grid(num_block)%cons(2*bfsize_i+index_i-1,index_j,2)

!       block_grid(num_block)%cons(index_i,index_j,3) &
!     = block_grid(num_block)%cons(2*bfsize_i+index_i-1,index_j,3)

!       block_grid(num_block)%cons(index_i,index_j,4) &
!     = block_grid(num_block)%cons(2*bfsize_i+index_i-1,index_j,4)

!     !Prims
!       block_grid(num_block)%prim(index_i,index_j,1) &
!     = block_grid(num_block)%prim(2*bfsize_i+index_i-1,index_j,1)

!       block_grid(num_block)%prim(index_i,index_j,2) &
!     =-block_grid(num_block)%prim(2*bfsize_i+index_i-1,index_j,2)

!       block_grid(num_block)%prim(index_i,index_j,3) &
!     = block_grid(num_block)%prim(2*bfsize_i+index_i-1,index_j,3)

!       block_grid(num_block)%prim(index_i,index_j,4) &
!     = block_grid(num_block)%prim(2*bfsize_i+index_i-1,index_j,4)

!       block_grid(num_block)%prim(index_i,index_j,5) &
!     = block_grid(num_block)%prim(2*bfsize_i+index_i-1,index_j,5)

!     enddo
! enddo

! !right
! do index_j=bfsize_j+1,nj+bfsize_j
!     do index_i=ni+bfsize_i+1,ni+2*bfsize_i
!     !Cons
!       block_grid(num_block)%cons(index_i,index_j,1) &
!     = block_grid(num_block)%cons(2*ni+2*bfsize_i+1-index_i,index_j,1)

!       block_grid(num_block)%cons(index_i,index_j,2) &
!     =-block_grid(num_block)%cons(2*ni+2*bfsize_i+1-index_i,index_j,2)

!       block_grid(num_block)%cons(index_i,index_j,3) &
!     = block_grid(num_block)%cons(2*ni+2*bfsize_i+1-index_i,index_j,3)

!       block_grid(num_block)%cons(index_i,index_j,4) &
!     = block_grid(num_block)%cons(2*ni+2*bfsize_i+1-index_i,index_j,4)

!     !Prims
!       block_grid(num_block)%prim(index_i,index_j,1) &
!     = block_grid(num_block)%prim(2*ni+2*bfsize_i+1-index_i,index_j,1)

!       block_grid(num_block)%prim(index_i,index_j,2) &
!     =-block_grid(num_block)%prim(2*ni+2*bfsize_i+1-index_i,index_j,2)

!       block_grid(num_block)%prim(index_i,index_j,3) &
!     = block_grid(num_block)%prim(2*ni+2*bfsize_i+1-index_i,index_j,3)

!       block_grid(num_block)%prim(index_i,index_j,4) &
!     = block_grid(num_block)%prim(2*ni+2*bfsize_i+1-index_i,index_j,4)

!       block_grid(num_block)%prim(index_i,index_j,5) &
!     = block_grid(num_block)%prim(2*ni+2*bfsize_i+1-index_i,index_j,5)
     
!      enddo
! enddo

do index_eq=1,neq
    do index_j=bfsize_j+1,nj+bfsize_j
        do index_i=1,bfsize_i
              block_grid(num_block)%cons(index_i,index_j,index_eq) &
            = block_grid(num_block)%cons(block_grid(num_block)%ni+index_i,index_j,index_eq)
              
              block_grid(num_block)%cons(block_grid(num_block)%ni+block_grid(num_block)%bfsize_i+index_i,index_j,index_eq) &
            = block_grid(num_block)%cons(block_grid(num_block)%bfsize_i+index_i,index_j,index_eq)
        enddo
    enddo
enddo

do index_var=1,prim_lastvar_index
    do index_j=bfsize_j+1,nj+bfsize_j
        do index_i=1,bfsize_i
              block_grid(num_block)%prim(index_i,index_j,index_var) &
            = block_grid(num_block)%prim(block_grid(num_block)%ni+index_i,index_j,index_var)
              
              block_grid(num_block)%prim(block_grid(num_block)%ni+block_grid(num_block)%bfsize_i+index_i,index_j,index_var) &
            = block_grid(num_block)%prim(block_grid(num_block)%bfsize_i+index_i,index_j,index_var)
        enddo
    enddo
enddo
 

block_grid(num_block)%bodyforce(:,:,3)=block_grid(num_block)%cons(:,:,1)
block_grid(num_block)%bodyforce(:,:,4)=block_grid(num_block)%cons(:,:,3)    
end subroutine
