subroutine init_ivc_general
use flow_general
implicit none
double precision :: rho,u,v,p,t
double precision :: b,r
integer :: index_i,index_j

num_block = 1
call allocate_block_grid
call allocate_block(block_grid(num_block),single_ni,single_nj,4,4,neq)

b=0.5d0
block_grid(num_block)%xl=10.d0
block_grid(num_block)%yl=10.d0
block_grid(num_block)%dx=block_grid(num_block)%xl/(block_grid(num_block)%ni)
block_grid(num_block)%dy=block_grid(num_block)%yl/(block_grid(num_block)%nj)

do index_j=1,block_grid(num_block)%nj
	do index_i=1,block_grid(num_block)%ni
		block_grid(num_block)%x(index_i,index_j)=(index_i-0.5)*block_grid(num_block)%dx
		block_grid(num_block)%y(index_i,index_j)=(index_j-0.5)*block_grid(num_block)%dy
	enddo
enddo


do index_j=1,block_grid(num_block)%nj+1
	do index_i=1,block_grid(num_block)%ni+1
		block_grid(num_block)%vet_x(index_i,index_j)=(index_i-1.0)*block_grid(num_block)%dx
		block_grid(num_block)%vet_y(index_i,index_j)=(index_j-1.0)*block_grid(num_block)%dy
	enddo
enddo

do index_j=1,block_grid(num_block)%nj
	do index_i=1,block_grid(num_block)%ni
		r=sqrt((block_grid(num_block)%x(index_i,index_j)-5.d0)**2+(block_grid(num_block)%y(index_i,index_j)-5.d0)**2)
		
		rho=(1.d0-(gamma-1.d0)*b**2*exp(1-r**2)/(8.d0*gamma*pie**2))**2.5d0
		
		u=0.0-b/(2*pie)*exp((1-r**2)/2)*(block_grid(num_block)%y(index_i,index_j)-5.d0)
		
		v=b/(2*pie)*exp((1-r**2)/2)*(block_grid(num_block)%x(index_i,index_j)-5.d0)
		
		p=rho**(gamma)

		block_grid(num_block)%cons(index_i+block_grid(num_block)%bfsize_i,index_j+block_grid(num_block)%bfsize_j,1)=rho
		block_grid(num_block)%cons(index_i+block_grid(num_block)%bfsize_i,index_j+block_grid(num_block)%bfsize_j,2)=rho*u
		block_grid(num_block)%cons(index_i+block_grid(num_block)%bfsize_i,index_j+block_grid(num_block)%bfsize_j,3)=rho*v
		block_grid(num_block)%cons(index_i+block_grid(num_block)%bfsize_i,index_j+block_grid(num_block)%bfsize_j,4)=p/(gamma-1.d0)+0.5*rho*(u**2+v**2)    
		
		block_grid(num_block)%prim(index_i+block_grid(num_block)%bfsize_i,index_j+block_grid(num_block)%bfsize_j,1)=rho
		block_grid(num_block)%prim(index_i+block_grid(num_block)%bfsize_i,index_j+block_grid(num_block)%bfsize_j,2)=u
		block_grid(num_block)%prim(index_i+block_grid(num_block)%bfsize_i,index_j+block_grid(num_block)%bfsize_j,3)=v
		block_grid(num_block)%prim(index_i+block_grid(num_block)%bfsize_i,index_j+block_grid(num_block)%bfsize_j,4)=p
		block_grid(num_block)%prim(index_i+block_grid(num_block)%bfsize_i,index_j+block_grid(num_block)%bfsize_j,5)=sqrt(gamma*p/rho)
	
	enddo
enddo   
end subroutine
