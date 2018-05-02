subroutine time_rk3_general(irk)
use flow_general
implicit none
integer,intent(in) :: irk
integer :: index_i,index_j,index_eq,index_block
integer :: bfsize_i,bfsize_j,ni,nj
double precision :: dt_xi,dt_eta,dxi,deta
double precision :: dh_dxi,dh_deta,q,q0,q1
select case (irk)
case(1)
	do index_block=1,num_block
	    ni = block_grid(index_block)%ni
	    nj = block_grid(index_block)%nj
	    bfsize_i=block_grid(index_block)%bfsize_i
	    bfsize_j=block_grid(index_block)%bfsize_j

	    dxi = block_grid(index_block)%dx
	    deta = block_grid(index_block)%dy

		if(.not.if_const_dt) then 
			dt_xi=block_grid(index_block)%dx/max_eig_i
			dt_eta=block_grid(index_block)%dy/max_eig_j
			dt=cfl*dt_xi*dt_eta/(dt_xi+dt_eta)
		endif

		do index_eq=1,neq
			call copy_field(block_grid(index_block)%cons_temp(:,:,index_eq), &
							block_grid(index_block)%cons(:,:,index_eq),      &
							1,                                 &
							ni+2*bfsize_i,                     &
							1,                                 &
							nj+2*bfsize_j)
		enddo

		do index_eq=1,neq
			!$omp parallel do default(shared) &
			!$omp private(index_i,index_j)	  &
			!$omp private(q,q0) 			  &
			!$omp private(dh_dxi,dh_deta)
			do index_j=bfsize_j+1,nj+bfsize_j
				do index_i=bfsize_i+1,ni+bfsize_i
					q0 = block_grid(index_block)%cons(index_i,index_j,index_eq)
					
					dh_dxi = block_grid(index_block)%h_xi(index_i-bfsize_i+1,index_j-bfsize_j,index_eq) &
							-block_grid(index_block)%h_xi(index_i-bfsize_i,index_j-bfsize_j,index_eq)
					dh_deta = block_grid(index_block)%h_eta(index_j-bfsize_j+1,index_i-bfsize_i,index_eq) &
							-block_grid(index_block)%h_eta(index_j-bfsize_j,index_i-bfsize_i,index_eq)
					
					q = q0-dt/dxi*dh_dxi-dt/deta*dh_deta
					block_grid(index_block)%cons(index_i,index_j,index_eq)= q
				end do
			end do
			!$omp end parallel do
		end do

		if(if_source) then
			do index_eq=1,neq
			!$omp parallel do default(shared) &
			!$omp private(index_i,index_j)
			do index_j=block_grid(index_block)%bfsize_j+1,block_grid(index_block)%nj+block_grid(index_block)%bfsize_j
				do index_i=block_grid(index_block)%bfsize_i+1,block_grid(index_block)%ni+block_grid(index_block)%bfsize_i
					block_grid(index_block)%cons(index_i,index_j,index_eq)= &
					block_grid(index_block)%cons(index_i,index_j,index_eq) &
						+dt*block_grid(index_block)%source(index_i,index_j,index_eq)
				end do
			end do
			!$omp end parallel do
		end do
		endif

		if(if_bodyforce) then
			do index_eq=1,neq
			!$omp parallel do default(shared) &
			!$omp private(index_i,index_j)
			do index_j=block_grid(index_block)%bfsize_j+1,block_grid(index_block)%nj+block_grid(index_block)%bfsize_j
				do index_i=block_grid(index_block)%bfsize_i+1,block_grid(index_block)%ni+block_grid(index_block)%bfsize_i
					block_grid(index_block)%cons(index_i,index_j,index_eq)= &
					block_grid(index_block)%cons(index_i,index_j,index_eq) &
						+dt*block_grid(index_block)%bodyforce(index_i,index_j,index_eq)
				end do
			end do
			!$omp end parallel do
		end do
		endif

		if(if_viscous) then
			do index_eq=1,neq
			!$omp parallel do default(shared) &
			!$omp private(index_i,index_j)
			do index_j=block_grid(index_block)%bfsize_j+1,block_grid(index_block)%nj+block_grid(index_block)%bfsize_j
				do index_i=block_grid(index_block)%bfsize_i+1,block_grid(index_block)%ni+block_grid(index_block)%bfsize_i
					block_grid(index_block)%cons(index_i,index_j,index_eq)= &
					block_grid(index_block)%cons(index_i,index_j,index_eq) &
						-dt/block_grid(index_block)%dx*(block_grid(index_block)%vis_xi(index_i,index_j,index_eq) &
							    -block_grid(index_block)%vis_xi(index_i-1,index_j,index_eq)) &
						-dt/block_grid(index_block)%dy*(block_grid(index_block)%vis_eta(index_j,index_i,index_eq) &
							     -block_grid(index_block)%vis_eta(index_j-1,index_i,index_eq))
				end do
			end do
			!$omp end parallel do
		end do
		endif
	enddo
case(2)
	do index_block=1,num_block
	    ni = block_grid(index_block)%ni
	    nj = block_grid(index_block)%nj
	    bfsize_i=block_grid(index_block)%bfsize_i
	    bfsize_j=block_grid(index_block)%bfsize_j

	    dxi = block_grid(index_block)%dx
	    deta = block_grid(index_block)%dy
		do index_eq=1,neq
			!$omp parallel do default(shared) &
			!$omp private(index_i,index_j)	  &
			!$omp private(q,q0,q1) 			  &
			!$omp private(dh_dxi,dh_deta)
			do index_j=bfsize_j+1,nj+bfsize_j
				do index_i=bfsize_i+1,ni+bfsize_i
					q0 = block_grid(index_block)%cons_temp(index_i,index_j,index_eq)
					q1 = block_grid(index_block)%cons(index_i,index_j,index_eq)
					
					dh_dxi = block_grid(index_block)%h_xi(index_i-bfsize_i+1,index_j-bfsize_j,index_eq) &
							-block_grid(index_block)%h_xi(index_i-bfsize_i,index_j-bfsize_j,index_eq)
					dh_deta = block_grid(index_block)%h_eta(index_j-bfsize_j+1,index_i-bfsize_i,index_eq) &
							-block_grid(index_block)%h_eta(index_j-bfsize_j,index_i-bfsize_i,index_eq)
					
					q = 3.d0/4.d0*q0 &
					   +1.d0/4.d0*q1 &
					   -1.d0/4.d0*(dt/dxi*dh_dxi+dt/deta*dh_deta)
					block_grid(index_block)%cons(index_i,index_j,index_eq)= q
				end do
			end do
			!$omp end parallel do
		end do

		if(if_source) then
			do index_eq=1,neq
			!$omp parallel do default(shared) &
			!$omp private(index_i,index_j)
			do index_j=block_grid(index_block)%bfsize_j+1,block_grid(index_block)%nj+block_grid(index_block)%bfsize_j
				do index_i=block_grid(index_block)%bfsize_i+1,block_grid(index_block)%ni+block_grid(index_block)%bfsize_i
					block_grid(index_block)%cons(index_i,index_j,index_eq)= &
					block_grid(index_block)%cons(index_i,index_j,index_eq) &
						+1.d0/4.d0*dt*block_grid(index_block)%source(index_i,index_j,index_eq)
				end do
			end do
			!$omp end parallel do
		end do
		endif

		if(if_bodyforce) then
			do index_eq=1,neq
			!$omp parallel do default(shared) &
			!$omp private(index_i,index_j)
			do index_j=block_grid(index_block)%bfsize_j+1,block_grid(index_block)%nj+block_grid(index_block)%bfsize_j
				do index_i=block_grid(index_block)%bfsize_i+1,block_grid(index_block)%ni+block_grid(index_block)%bfsize_i
					block_grid(index_block)%cons(index_i,index_j,index_eq)= &
					block_grid(index_block)%cons(index_i,index_j,index_eq) &
						+1.d0/4.d0*dt*block_grid(index_block)%bodyforce(index_i,index_j,index_eq)
				end do
			end do
			!$omp end parallel do
		end do
		endif

		if(if_viscous) then
			do index_eq=1,neq
			!$omp parallel do default(shared) &
			!$omp private(index_i,index_j)
			do index_j=block_grid(index_block)%bfsize_j+1,block_grid(index_block)%nj+block_grid(index_block)%bfsize_j
				do index_i=block_grid(index_block)%bfsize_i+1,block_grid(index_block)%ni+block_grid(index_block)%bfsize_i
					block_grid(index_block)%cons(index_i,index_j,index_eq)= &
					block_grid(index_block)%cons(index_i,index_j,index_eq) &
						+1.d0/4.d0*dt/block_grid(index_block)%dx*(block_grid(index_block)%vis_xi(index_i,index_j,index_eq) &
							    -block_grid(index_block)%vis_xi(index_i-1,index_j,index_eq)) &
						+1.d0/4.d0*dt/block_grid(index_block)%dy*(block_grid(index_block)%vis_eta(index_j,index_i,index_eq) &
							     -block_grid(index_block)%vis_eta(index_j-1,index_i,index_eq))
				end do
			end do
			!$omp end parallel do
		end do
		endif
	enddo

case(3)
	do index_block=1,num_block
	    ni = block_grid(index_block)%ni
	    nj = block_grid(index_block)%nj
	    bfsize_i=block_grid(index_block)%bfsize_i
	    bfsize_j=block_grid(index_block)%bfsize_j

	    dxi = block_grid(index_block)%dx
	    deta = block_grid(index_block)%dy
		do index_eq=1,neq
			!$omp parallel do default(shared) &
			!$omp private(index_i,index_j)    &
			!$omp private(q,q0,q1) 			  &
			!$omp private(dh_dxi,dh_deta)
			do index_j=bfsize_j+1,nj+bfsize_j
				do index_i=bfsize_i+1,ni+bfsize_i
					q0 = block_grid(index_block)%cons_temp(index_i,index_j,index_eq)
					q1 = block_grid(index_block)%cons(index_i,index_j,index_eq)
					
					dh_dxi = block_grid(index_block)%h_xi(index_i-bfsize_i+1,index_j-bfsize_j,index_eq) &
							-block_grid(index_block)%h_xi(index_i-bfsize_i,index_j-bfsize_j,index_eq)
					dh_deta = block_grid(index_block)%h_eta(index_j-bfsize_j+1,index_i-bfsize_i,index_eq) &
							-block_grid(index_block)%h_eta(index_j-bfsize_j,index_i-bfsize_i,index_eq)
					
					q = 1.d0/3.d0*q0 &
					   +2.d0/3.d0*q1 &
					   -2.d0/3.d0*(dt/dxi*dh_dxi+dt/deta*dh_deta)
					block_grid(index_block)%cons(index_i,index_j,index_eq)= q
				end do
			end do
			!$omp end parallel do
		end do

		if(if_source) then
			do index_eq=1,neq
			!$omp parallel do default(shared) &
			!$omp private(index_i,index_j)
			do index_j=block_grid(index_block)%bfsize_j+1,block_grid(index_block)%nj+block_grid(index_block)%bfsize_j
				do index_i=block_grid(index_block)%bfsize_i+1,block_grid(index_block)%ni+block_grid(index_block)%bfsize_i
					block_grid(index_block)%cons(index_i,index_j,index_eq)= &
					block_grid(index_block)%cons(index_i,index_j,index_eq) &
						+2.0/3.0*dt*block_grid(index_block)%source(index_i,index_j,index_eq)
				end do
			end do
			!$omp end parallel do
		end do
		endif

		if(if_bodyforce) then
			do index_eq=1,neq
			!$omp parallel do default(shared) &
			!$omp private(index_i,index_j)
			do index_j=block_grid(index_block)%bfsize_j+1,block_grid(index_block)%nj+block_grid(index_block)%bfsize_j
				do index_i=block_grid(index_block)%bfsize_i+1,block_grid(index_block)%ni+block_grid(index_block)%bfsize_i
					block_grid(index_block)%cons(index_i,index_j,index_eq)= &
					block_grid(index_block)%cons(index_i,index_j,index_eq) &
						+2.0/3.0*dt*block_grid(index_block)%bodyforce(index_i,index_j,index_eq)
				end do
			end do
			!$omp end parallel do
		end do
		endif

		if(if_viscous) then
			do index_eq=1,neq
			!$omp parallel do default(shared) &
			!$omp private(index_i,index_j)
			do index_j=block_grid(index_block)%bfsize_j+1,block_grid(index_block)%nj+block_grid(index_block)%bfsize_j
				do index_i=block_grid(index_block)%bfsize_i+1,block_grid(index_block)%ni+block_grid(index_block)%bfsize_i
					block_grid(index_block)%cons(index_i,index_j,index_eq)= &
					block_grid(index_block)%cons(index_i,index_j,index_eq) &
						+2.0/3.0*dt/block_grid(index_block)%dx*(block_grid(index_block)%vis_xi(index_i,index_j,index_eq) &
							    -block_grid(index_block)%vis_xi(index_i-1,index_j,index_eq)) &
						+2.0/3.0*dt/block_grid(index_block)%dy*(block_grid(index_block)%vis_eta(index_j,index_i,index_eq) &
							     -block_grid(index_block)%vis_eta(index_j-1,index_i,index_eq))
				end do
			end do
			!$omp end parallel do
		end do
		endif
	enddo
end select

end subroutine

