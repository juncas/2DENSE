subroutine read_mesh_general
use flow_general
implicit none
integer :: index_i,index_j,index_block

open(21,file='grid.dat',form='formatted')

read(21) num_block

call allocate_block_grid

do index_block=1,num_block
    read(10) block_grid(index_block)%ni,       &
             block_grid(index_block)%nj,       &
             block_grid(index_block)%bfsize_i, &
             block_grid(index_block)%bfsize_j, &
             block_grid(index_block)%block_topo
enddo

do index_block=1,num_block
    
    call allocate_block(block_grid(index_block),          &
                        block_grid(index_block)%ni,       &
                        block_grid(index_block)%nj,       &
                        block_grid(index_block)%bfsize_i, &
                        block_grid(index_block)%bfsize_j, &
                        neq)
    read(10) block_grid(index_block)%x
    read(10) block_grid(index_block)%y
enddo
close(21)

if(if_curve_grid) then
    call build_jacobians
endif

end subroutine

subroutine build_jacobians
use flow_general
implicit none
integer :: index_i,index_j,index_block

read(10) block_grid(index_block)%xi_x
read(10) block_grid(index_block)%eta_x
read(10) block_grid(index_block)%xi_y
read(10) block_grid(index_block)%eta_y
read(10) block_grid(index_block)%jacob

end subroutine