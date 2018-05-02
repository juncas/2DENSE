subroutine write_tec_at_time
use flow_general
implicit none

if(if_curve_grid) then
	call write_tec_at_time_curve
else
	call write_tec_at_time_cartisian
endif
if(if_save_theta) then
	call write_theta_at_time_cartisian
endif
end subroutine


subroutine write_tec_at_time_cartisian
use flow_general
use var_index
implicit none
integer :: index_block,index_var
integer :: index_i,index_j
character(len=10) :: temp 

write(temp,'(i8)') current_step
open(10,file="2d_"//trim(adjustl(temp))//".dat",form='formatted')   

write(10,*)'title=','"flow time ',current_time,'"'
write(10,*)'variables="x","y","rho","u","v","p","alpha"'
do index_block=1,num_block
	write(10,*) 'ZONE I=',block_grid(num_block)%ni+1, &
				    ' J=',block_grid(num_block)%nj+1, &
				    'DATAPACKING=BLOCK, VARLOCATION=([3-7]=CELLCENTERED)'

	do index_j=1,block_grid(num_block)%nj+1
		do index_i=1,block_grid(num_block)%ni+1
			write(10,*) block_grid(num_block)%vet_x(index_i,index_j)
		end do
	end do


	do index_j=1,block_grid(num_block)%nj+1
		do index_i=1,block_grid(num_block)%ni+1
			write(10,*) block_grid(num_block)%vet_y(index_i,index_j)
		end do
	end do

	do index_var=1,prim_lastvar_index
		do index_j=block_grid(num_block)%bfsize_j+1,block_grid(num_block)%nj+block_grid(num_block)%bfsize_j
			do index_i=block_grid(num_block)%bfsize_i+1,block_grid(num_block)%ni+block_grid(num_block)%bfsize_i
				write(10,*) block_grid(num_block)%prim(index_i,index_j,index_var)
			end do
		end do
	enddo

enddo
close(10)

end subroutine


subroutine write_theta_at_time_cartisian
use flow_general
use var_index
implicit none
integer :: index_block,index_var
integer :: index_i,index_j
character(len=10) :: temp 

write(temp,'(i8)') current_step
open(10,file="2d_theta_"//trim(adjustl(temp))//".dat",form='formatted')   

write(10,*)'title=','"flow time ',current_time,'"'
write(10,*)'variables="x","y","xp","xn","ep","en"'
do index_block=1,num_block
	write(10,*) 'ZONE I=',block_grid(num_block)%ni+1, &
				    ' J=',block_grid(num_block)%nj+1, &
				    'DATAPACKING=BLOCK, VARLOCATION=([3-6]=CELLCENTERED)'

	do index_j=1,block_grid(num_block)%nj+1
		do index_i=1,block_grid(num_block)%ni+1
			write(10,*) block_grid(num_block)%vet_x(index_i,index_j)
		end do
	end do


	do index_j=1,block_grid(num_block)%nj+1
		do index_i=1,block_grid(num_block)%ni+1
			write(10,*) block_grid(num_block)%vet_y(index_i,index_j)
		end do
	end do

		do index_j=1,block_grid(num_block)%nj
			do index_i=1,block_grid(num_block)%ni
				write(10,*) max(block_grid(num_block)%theta_xip(1,index_i,index_j),block_grid(num_block)%theta_xip(1,index_i+1,index_j))
			end do
		end do

		do index_j=1,block_grid(num_block)%nj
			do index_i=1,block_grid(num_block)%ni
				write(10,*) max(block_grid(num_block)%theta_xin(1,index_i,index_j),block_grid(num_block)%theta_xin(1,index_i+1,index_j))
			end do
		end do

		do index_j=1,block_grid(num_block)%nj
			do index_i=1,block_grid(num_block)%ni
				write(10,*) max(block_grid(num_block)%theta_etap(1,index_j,index_i),block_grid(num_block)%theta_etap(1,index_j+1,index_i))
			end do
		end do

		do index_j=1,block_grid(num_block)%nj
			do index_i=1,block_grid(num_block)%ni
				write(10,*) max(block_grid(num_block)%theta_etan(1,index_j,index_i),block_grid(num_block)%theta_etan(1,index_j+1,index_i))
			end do
		end do

enddo
close(10)

end subroutine write_theta_at_time_cartisian

subroutine write_tec_at_time_curve
! use flow_general
! double precision :: rho,u,v,p,t
! character(len=10) :: temp 

! write(temp,'(i8)') current_step
! print*,temp
! write(10,*)'title=','"flow time ',current_time,'"'
! write(10,*)'variables="x","y","rho","u","v","p","t"'

! write(10,*)'zone i=',ni,'j=',jl
! do j=bfsize+1,jl+bfsize
! 	do i=bfsize+1,il+bfsize
! 		rho=block_grid(num_block)%cons(index_i,index_j,1)
! 		u=block_grid(num_block)%cons(index_i,index_j,2)/rho
! 		v=block_grid(num_block)%cons(index_i,index_j,3)/rho
! 		p=(gamma-1.)*(block_grid(num_block)%cons(index_i,index_j,4)-0.5*rho*(u**2+v**2))
! 		t=p/rho
! 		write(10,*) x(i,j),y(i,j),rho,u,v,p,t 
! 	end do
! end do
! close(10)

end subroutine

