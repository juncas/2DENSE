subroutine boundary_shockvortex_general
use flow_general
use var_index
implicit none
integer :: index_i,index_j,index_var,index_eq
integer :: bfsize_i,bfsize_j,ni,nj
double precision :: rho,u,v,p
double precision :: x,y

ni = block_grid(num_block)%ni
nj = block_grid(num_block)%nj
bfsize_i=block_grid(num_block)%bfsize_i
bfsize_j=block_grid(num_block)%bfsize_j

! Bottom reflective boundary
do index_j=1,bfsize_j
    do index_i=bfsize_i+1,ni+bfsize_i
              block_grid(num_block)%cons(index_i,index_j,1) &
            = block_grid(num_block)%cons(index_i,2*bfsize_j-index_j+1,1) 
              
              block_grid(num_block)%cons(index_i,index_j,2) &
            = block_grid(num_block)%cons(index_i,2*bfsize_j-index_j+1,2) 
              
              block_grid(num_block)%cons(index_i,index_j,3) &
            =-block_grid(num_block)%cons(index_i,2*bfsize_j-index_j+1,3) 
              
              block_grid(num_block)%cons(index_i,index_j,4) &
            = block_grid(num_block)%cons(index_i,2*bfsize_j-index_j+1,4) 
              
              block_grid(num_block)%prim(index_i,index_j,1) &
            = block_grid(num_block)%prim(index_i,2*bfsize_j-index_j+1,1) 
              
              block_grid(num_block)%prim(index_i,index_j,2) &
            = block_grid(num_block)%prim(index_i,2*bfsize_j-index_j+1,2) 
              
              block_grid(num_block)%prim(index_i,index_j,3) &
            =-block_grid(num_block)%prim(index_i,2*bfsize_j-index_j+1,3) 
              
              block_grid(num_block)%prim(index_i,index_j,4) &
            = block_grid(num_block)%prim(index_i,2*bfsize_j-index_j+1,4) 

              block_grid(num_block)%prim(index_i,index_j,5) &
            = block_grid(num_block)%prim(index_i,2*bfsize_j-index_j+1,5) 
    enddo
enddo

! Top reflective boundary
do index_j=nj+bfsize_j+1,nj+2*bfsize_j
    do index_i=bfsize_i+1,ni+bfsize_i

        block_grid(num_block)%cons(index_i,index_j,1) &
        = block_grid(num_block)%cons(index_i,2*nj+2*bfsize_j-index_j+1,1) 

        block_grid(num_block)%cons(index_i,index_j,2) &
        = block_grid(num_block)%cons(index_i,2*nj+2*bfsize_j-index_j+1,2) 

        block_grid(num_block)%cons(index_i,index_j,3) &
        =-block_grid(num_block)%cons(index_i,2*nj+2*bfsize_j-index_j+1,3) 

        block_grid(num_block)%cons(index_i,index_j,4) &
        = block_grid(num_block)%cons(index_i,2*nj+2*bfsize_j-index_j+1,4) 

        block_grid(num_block)%prim(index_i,index_j,1) &
        = block_grid(num_block)%prim(index_i,2*nj+2*bfsize_j-index_j+1,1) 

        block_grid(num_block)%prim(index_i,index_j,2) &
        = block_grid(num_block)%prim(index_i,2*nj+2*bfsize_j-index_j+1,2) 

        block_grid(num_block)%prim(index_i,index_j,3) &
        =-block_grid(num_block)%prim(index_i,2*nj+2*bfsize_j-index_j+1,3) 

        block_grid(num_block)%prim(index_i,index_j,4) &
        = block_grid(num_block)%prim(index_i,2*nj+2*bfsize_j-index_j+1,4) 

        block_grid(num_block)%prim(index_i,index_j,5) &
        = block_grid(num_block)%prim(index_i,2*nj+2*bfsize_j-index_j+1,5) 
    
    enddo
enddo
              
! Left outlet boundary
do index_eq = 1,neq
    do index_j=bfsize_j+1,nj+bfsize_j
        do index_i=bfsize_i,1,-1
              block_grid(num_block)%cons(index_i,index_j,index_eq) &
            = block_grid(num_block)%cons(index_i+1,index_j,index_eq)*4.0/3.0 &
             -block_grid(num_block)%cons(index_i+2,index_j,index_eq)/3.0
        enddo
    enddo
enddo

do index_eq = 1,prim_lastvar_index
    do index_j=bfsize_j+1,nj+bfsize_j
        do index_i=bfsize_i,1,-1
              block_grid(num_block)%prim(index_i,index_j,index_eq) &
            = block_grid(num_block)%prim(index_i+1,index_j,index_eq)*4.0/3.0 &
             -block_grid(num_block)%prim(index_i+2,index_j,index_eq)/3.0
        enddo
    enddo
enddo


! Right outlet boundary
do index_eq = 1,neq
    do index_j=bfsize_j+1,nj+bfsize_j
        do index_i=ni+bfsize_i+1,ni+2*bfsize_i
              block_grid(num_block)%cons(index_i,index_j,index_eq) &
            = block_grid(num_block)%cons(index_i-1,index_j,index_eq)*4.0/3.0 &
             -block_grid(num_block)%cons(index_i-2,index_j,index_eq)/3.0
        enddo
    enddo
enddo

do index_eq = 1,prim_lastvar_index
    do index_j=bfsize_j+1,nj+bfsize_j
        do index_i=ni+bfsize_i+1,ni+2*bfsize_i
              block_grid(num_block)%prim(index_i,index_j,index_eq) &
            = block_grid(num_block)%prim(index_i-1,index_j,index_eq)*4.0/3.0 &
             -block_grid(num_block)%prim(index_i-2,index_j,index_eq)/3.0
        enddo
    enddo
enddo
              

end subroutine
