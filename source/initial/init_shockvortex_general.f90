subroutine init_shockvortex_general
use flow_general
implicit none
double precision :: rho,u,v,p,t,x,y
double precision :: x_c=0.25,y_c=0.5
double precision :: r,tau
double precision :: epsl=0.3,r_c=0.05,alpha=0.204
integer :: index_i,index_j

num_block = 1
call allocate_block_grid
call allocate_block(block_grid(num_block),single_ni,single_nj,4,4,neq)

block_grid(num_block)%xl=2.d0
block_grid(num_block)%yl=1.d0
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
		x   =block_grid(num_block)%x(index_i,index_j)
		y   =block_grid(num_block)%y(index_i,index_j)
		r   =sqrt((x-x_c)**2+(y-y_c)**2)
		tau = r/r_c

		if(x<0.5) then
			rho = (1.d0-(gamma-1.d0)*epsl**2*exp(2*alpha*(1-tau**2))/(4*alpha*gamma))**(1.0/(gamma-1))
			u   = 1.1*sqrt(gamma)+epsl*tau*exp(alpha*(1.0-tau**2))*(y-y_c)/r
			v   =-epsl*tau*exp(alpha*(1.0-tau**2))*(x-x_c)/r
			p   = rho**gamma
		else
			rho = 1.169082d0
			u   = 1.113299d0
			v   = 0.0d0
			p   = 1.245000d0
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
	
	enddo
enddo   

end subroutine
