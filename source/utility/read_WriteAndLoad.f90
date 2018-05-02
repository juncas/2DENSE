subroutine write_date_binary
use flow_general
implicit none
character(len=10) :: temp 
integer :: index_block
write(temp,'(i8)') current_step
open(10,file="2d_"//trim(adjustl(temp))//".dat2d",form='unformatted')     

write(10) current_step,current_time       
write(10) num_block 
do index_block=1,num_block
    write(10) block_grid(index_block)%index
    write(10) block_grid(index_block)%ni,block_grid(index_block)%nj
    write(10) block_grid(index_block)%bfsize_i,block_grid(index_block)%bfsize_j
	write(10) block_grid(index_block)%cons
	write(10) block_grid(index_block)%x
	write(10) block_grid(index_block)%y
	write(10) block_grid(index_block)%dx,block_grid(index_block)%dy
	if(if_curve_grid) then
		write(10) block_grid(index_block)%xi_x
		write(10) block_grid(index_block)%eta_x
		write(10) block_grid(index_block)%xi_y
		write(10) block_grid(index_block)%eta_y
		write(10) block_grid(index_block)%jacob
	endif
enddo
close(10)
end subroutine

subroutine load_date_binary
use flow_general
implicit none
character(len=10) :: temp  
integer :: index_block

print*,"Loading data at step=",start_step
write(temp,'(i8)') start_step
open(10,file="2d_"//trim(adjustl(temp))//".dat2d",form='unformatted')        
read(10) current_step,current_time        
read(10) num_block
call allocate_block_grid

print*,"Flow time=",current_time 
do index_block=1,num_block
    read(10) block_grid(index_block)%index
    read(10) block_grid(index_block)%ni,block_grid(index_block)%nj
    read(10) block_grid(index_block)%bfsize_i,block_grid(index_block)%bfsize_j

	call allocate_block(block_grid(index_block),          &
						block_grid(index_block)%ni,       &
						block_grid(index_block)%nj,       &
						block_grid(index_block)%bfsize_i, &
						block_grid(index_block)%bfsize_j, &
						neq)

	read(10) block_grid(index_block)%cons
	read(10) block_grid(index_block)%x
	read(10) block_grid(index_block)%y
	read(10) block_grid(index_block)%dx,block_grid(index_block)%dy
	if(if_curve_grid) then
		read(10) block_grid(index_block)%xi_x
		read(10) block_grid(index_block)%eta_x
		read(10) block_grid(index_block)%xi_y
		read(10) block_grid(index_block)%eta_y
		read(10) block_grid(index_block)%jacob
	endif
enddo
close(10)
print*,'Done!'
end subroutine
