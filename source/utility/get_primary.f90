subroutine get_primary
use flow_general
implicit none
integer :: block_index,ni,nj,bfsize_i,bfsize_j

max_eig_i=0.0
max_eig_j=0.0
do block_index = 1,num_block
	ni 		 = block_grid(block_index)%ni
	nj 		 = block_grid(block_index)%nj
	bfsize_i = block_grid(block_index)%bfsize_i
	bfsize_j = block_grid(block_index)%bfsize_j
	call get_primary_block(block_grid(block_index)%cons, &
		ni, &
		nj, &
		neq, &
		bfsize_i, &
		bfsize_j, &
		block_grid(block_index)%prim, &
		max_eig_i, &
		max_eig_j)
enddo
end subroutine

subroutine get_primary_block(cons,ni,nj,neq,bfsize_i,bfsize_j,prim,max_eig_i,max_eig_j)
use flow_general,only: gamma
use var_index
implicit none
! inputs
integer,intent(in) :: ni,nj ! number of interior grids
integer,intent(in) :: bfsize_i,bfsize_j ! number of ghost grid
integer,intent(in) :: neq ! number equations
double precision,intent(in),dimension(ni+2*bfsize_i,nj+2*bfsize_j,neq) :: cons

!outputs
double precision,intent(out),dimension(ni+2*bfsize_i,nj+2*bfsize_j,neq) :: prim
double precision,intent(out) :: max_eig_i,max_eig_j

!local variables
double precision :: rho,u,v,p,alpha
double precision :: max_eig_i_local=0.0, &
					max_eig_j_local=0.0
integer :: index_i,index_j,index_eq


!$OMP PARALLEL DO default(shared) &
!$omp private(rho,u,v,p,alpha) &
!$omp private(max_eig_i_local,max_eig_j_local) &
!$omp private(index_i,index_j)
		do index_j=bfsize_j+1,nj+bfsize_j
			do index_i=bfsize_i+1,ni+bfsize_i

				rho = cons(index_i,index_j,1)
				u   = cons(index_i,index_j,2)/rho
				v   = cons(index_i,index_j,3)/rho
				p   = (gamma-1.d0)*(cons(index_i,index_j,4)-0.5d0*rho*(u**2+v**2))
				if(p.lt.0.0) then
					p = 1.0e-6
				endif
				alpha=sqrt(gamma*p/rho)

				prim(index_i,index_j,prim_rho_index)     = rho 
				prim(index_i,index_j,prim_u_index)       = u
				prim(index_i,index_j,prim_v_index)       = v
				prim(index_i,index_j,prim_p_index)       = p
				prim(index_i,index_j,prim_alpha_index)   = alpha
				prim(index_i,index_j,prim_entropy_index) = log(p)-log(rho)

				!$omp atomic
				max_eig_i=MAX(max_eig_i,abs(u)+alpha)

				!$omp atomic
				max_eig_j=MAX(max_eig_j,abs(v)+alpha)
			enddo
		enddo
!$OMP end PARALLEL DO

end subroutine get_primary_block
