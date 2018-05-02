subroutine boundary_ivc_general
use flow_general
use var_index
implicit none
integer :: index_i,index_j,index_eq,index_var
integer :: ni,nj,bfsize_i,bfsize_j

ni = block_grid(num_block)%ni
nj = block_grid(num_block)%nj
bfsize_i=block_grid(num_block)%bfsize_i
bfsize_j=block_grid(num_block)%bfsize_j

do index_eq=1,neq
	do index_j=1,block_grid(num_block)%bfsize_j
		do index_i=block_grid(num_block)%bfsize_i+1,block_grid(num_block)%ni+block_grid(num_block)%bfsize_i
			  block_grid(num_block)%cons(index_i,index_j,index_eq) &
			= block_grid(num_block)%cons(index_i,nj+index_j,index_eq) 
			  
			  block_grid(num_block)%cons(index_i,nj+bfsize_j+index_j,index_eq) &
			= block_grid(num_block)%cons(index_i,bfsize_j+index_j,index_eq) 
		enddo
	enddo
enddo

do index_var=1,prim_lastvar_index
	do index_j=1,bfsize_j
		do index_i=bfsize_i+1,ni+bfsize_i
			  block_grid(num_block)%prim(index_i,index_j,index_var) &
			= block_grid(num_block)%prim(index_i,nj+index_j,index_var) 
			  
			  block_grid(num_block)%prim(index_i,nj+bfsize_j+index_j,index_var) &
			= block_grid(num_block)%prim(index_i,bfsize_j+index_j,index_var) 
		enddo
	enddo
enddo

! do index_eq=1,neq
! 	do index_j=block_grid(num_block)%bfsize_j+1,block_grid(num_block)%nj+block_grid(num_block)%bfsize_j
! 		do index_i=1,block_grid(num_block)%bfsize_i
! 			  block_grid(num_block)%cons(index_i,index_j,index_eq) &
! 			= block_grid(num_block)%cons(block_grid(num_block)%ni+index_i,index_j,index_eq)
			  
! 			  block_grid(num_block)%cons(block_grid(num_block)%ni+block_grid(num_block)%bfsize_i+index_i,index_j,index_eq) &
! 			= block_grid(num_block)%cons(block_grid(num_block)%bfsize_i+index_i,index_j,index_eq)
! 		enddo
! 	enddo
! enddo

! do index_var=1,prim_lastvar_index
! 	do index_j=block_grid(num_block)%bfsize_j+1,block_grid(num_block)%nj+block_grid(num_block)%bfsize_j
! 		do index_i=1,block_grid(num_block)%bfsize_i
! 			  block_grid(num_block)%prim(index_i,index_j,index_var) &
! 			= block_grid(num_block)%prim(block_grid(num_block)%ni+index_i,index_j,index_var)
			  
! 			  block_grid(num_block)%prim(block_grid(num_block)%ni+block_grid(num_block)%bfsize_i+index_i,index_j,index_var) &
! 			= block_grid(num_block)%prim(block_grid(num_block)%bfsize_i+index_i,index_j,index_var)
! 		enddo
! 	enddo
! enddo
    
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
