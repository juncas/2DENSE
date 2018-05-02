subroutine init_rt_general
use flow_general
implicit none
double precision :: rho,u,v,p,c
double precision :: b,r
integer :: index_i,index_j

num_block = 1
if_bodyforce = .true.
gamma=5.d0/3.d0
call allocate_block_grid
call allocate_block(block_grid(num_block),single_ni,single_nj,4,4,neq)

block_grid(num_block)%xl=0.25d0
block_grid(num_block)%yl=1.d0
block_grid(num_block)%dx=block_grid(num_block)%xl/(block_grid(num_block)%ni-1.d0)
block_grid(num_block)%dy=block_grid(num_block)%yl/(block_grid(num_block)%nj-1.d0)

do index_j=1,block_grid(num_block)%nj
	do index_i=1,block_grid(num_block)%ni
		block_grid(num_block)%x(index_i,index_j)=(index_i-1.d0)*block_grid(num_block)%dx
		block_grid(num_block)%y(index_i,index_j)=(index_j-1.d0)*block_grid(num_block)%dy
	enddo
enddo


do index_j=1,block_grid(num_block)%nj+1
	do index_i=1,block_grid(num_block)%ni+1
		block_grid(num_block)%vet_x(index_i,index_j)=(index_i-1.5)*block_grid(num_block)%dx
		block_grid(num_block)%vet_y(index_i,index_j)=(index_j-1.5)*block_grid(num_block)%dy
	enddo
enddo

do index_j=1,block_grid(num_block)%nj
	do index_i=1,block_grid(num_block)%ni
		if(block_grid(num_block)%y(index_i,index_j).le.0.5) then
			rho=2.d0
			p=2.d0*(block_grid(num_block)%y(index_i,index_j))+1.d0
			u=0.d0
			c=sqrt(gamma*p/rho)
			v=-0.025d0*c*cos(8.d0*pie*block_grid(num_block)%x(index_i,index_j))
		else
			rho=1.d0
			p=block_grid(num_block)%y(index_i,index_j)+1.5d0
			u=0.d0
			c=sqrt(gamma*p/rho)
			v=-0.025*c*cos(8.d0*pie*block_grid(num_block)%x(index_i,index_j))
		endif

		block_grid(num_block)%cons(index_i+block_grid(num_block)%bfsize_i,index_j+block_grid(num_block)%bfsize_j,1)=rho
		block_grid(num_block)%cons(index_i+block_grid(num_block)%bfsize_i,index_j+block_grid(num_block)%bfsize_j,2)=rho*u
		block_grid(num_block)%cons(index_i+block_grid(num_block)%bfsize_i,index_j+block_grid(num_block)%bfsize_j,3)=rho*v
		block_grid(num_block)%cons(index_i+block_grid(num_block)%bfsize_i,index_j+block_grid(num_block)%bfsize_j,4)=p/(gamma-1.d0)+0.5*rho*(u**2+v**2)    
		
		block_grid(num_block)%prim(index_i+block_grid(num_block)%bfsize_i,index_j+block_grid(num_block)%bfsize_j,1)=rho
		block_grid(num_block)%prim(index_i+block_grid(num_block)%bfsize_i,index_j+block_grid(num_block)%bfsize_j,2)=u
		block_grid(num_block)%prim(index_i+block_grid(num_block)%bfsize_i,index_j+block_grid(num_block)%bfsize_j,3)=v
		block_grid(num_block)%prim(index_i+block_grid(num_block)%bfsize_i,index_j+block_grid(num_block)%bfsize_j,4)=p
		block_grid(num_block)%prim(index_i+block_grid(num_block)%bfsize_i,index_j+block_grid(num_block)%bfsize_j,5)=sqrt(gamma*p/rho)
	
		block_grid(num_block)%bodyforce(index_i+block_grid(num_block)%bfsize_i,index_j+block_grid(num_block)%bfsize_j,1)=0.d0
		block_grid(num_block)%bodyforce(index_i+block_grid(num_block)%bfsize_i,index_j+block_grid(num_block)%bfsize_j,2)=0.d0
		block_grid(num_block)%bodyforce(index_i+block_grid(num_block)%bfsize_i,index_j+block_grid(num_block)%bfsize_j,3)=rho
		block_grid(num_block)%bodyforce(index_i+block_grid(num_block)%bfsize_i,index_j+block_grid(num_block)%bfsize_j,4)=rho*v
	enddo
enddo   
end subroutine
