subroutine scheme_upwind5_general
use flow_general
implicit none
integer :: index_block,index_eq

do index_block=1,num_block
	call scheme_upwind5_xip(block_grid(index_block)%flux_xip, &
						   block_grid(index_block)%ni,       &
						   block_grid(index_block)%nj,       &
						   block_grid(index_block)%bfsize_i, &
						   block_grid(index_block)%bfsize_j, &
						   neq,                              &
						   block_grid(index_block)%h_xip)

	call scheme_upwind5_xin(block_grid(index_block)%flux_xin, &
						   block_grid(index_block)%ni,       &
						   block_grid(index_block)%nj,       &
						   block_grid(index_block)%bfsize_i, &
						   block_grid(index_block)%bfsize_j, &
						   neq,                              &
						   block_grid(index_block)%h_xin)

	call scheme_upwind5_etap(block_grid(index_block)%flux_etap, &
						   block_grid(index_block)%ni,       &
						   block_grid(index_block)%nj,       &
						   block_grid(index_block)%bfsize_i, &
						   block_grid(index_block)%bfsize_j, &
						   neq,                              &
						   block_grid(index_block)%h_etap)

	call scheme_upwind5_etan(block_grid(index_block)%flux_etan, &
						   block_grid(index_block)%ni,       &
						   block_grid(index_block)%nj,       &
						   block_grid(index_block)%bfsize_i, &
						   block_grid(index_block)%bfsize_j, &
						   neq,                              &
						   block_grid(index_block)%h_etan)

	do index_eq=1,neq
		call sum_field(block_grid(index_block)%h_xi(:,:,index_eq),   &
					   block_grid(index_block)%h_xip(:,:,index_eq),  &
					   block_grid(index_block)%h_xin(:,:,index_eq),  &
						1,                                           &
						block_grid(index_block)%ni+1,                                        &
						1,                                           &
						block_grid(index_block)%nj)
	end do
	
	do index_eq=1,neq
		call sum_field(block_grid(index_block)%h_eta(:,:,index_eq),  &
					   block_grid(index_block)%h_etap(:,:,index_eq), &
					   block_grid(index_block)%h_etan(:,:,index_eq), &
						1,                                           &
						block_grid(index_block)%nj+1,                                        &
						1,                                           &
						block_grid(index_block)%ni)
	end do

	block_grid(index_block)%h_xi = block_grid(index_block)%h_xip+block_grid(index_block)%h_xin
	block_grid(index_block)%h_eta = block_grid(index_block)%h_etap+block_grid(index_block)%h_etan
enddo
end subroutine scheme_upwind5_general

subroutine scheme_upwind5_xip(flux_xip,ni,nj,bfsize_i,bfsize_j,neq,h_xip)
implicit none

integer,intent(in) :: ni,nj,bfsize_i,bfsize_j,neq
double precision,intent(in),dimension(ni+2*bfsize_i,nj,neq) :: flux_xip
double precision,intent(out),dimension(ni+1,nj,neq) :: h_xip

integer :: index_i,index_j,index_eq,m
double precision :: lf(-2:2)

do index_eq=1,neq
!$OMP PARALLEL DO default(shared) &
!$omp private(lf) 		      	  &
!$omp private(index_i,index_j,m)  


do index_j=1,nj
	do index_i=1,ni+1
        do m=-2,2
            lf(m)=flux_xip(index_i+m+bfsize_i-1,index_j,index_eq)
        enddo

		h_xip(index_i,index_j,index_eq) = (2.d0*lf(-2) &
										 -13.d0*lf(-1) &
										 +47.d0*lf(0) &
										 +27.d0*lf(1) &
										 - 3.d0*lf(2))/60.d0 

	enddo
enddo
!$OMP end PARALLEL DO
enddo
end subroutine

subroutine scheme_upwind5_xin(flux_xin,ni,nj,bfsize_i,bfsize_j,neq,h_xin)
implicit none
integer,intent(in) :: ni,nj,bfsize_i,bfsize_j,neq
double precision,intent(in),dimension(ni+2*bfsize_i,nj,neq) :: flux_xin
double precision,intent(out),dimension(ni+1,nj,neq) :: h_xin
integer :: index_i,index_j,index_eq,m
double precision :: lf(-2:2)

do index_eq=1,neq
!$OMP PARALLEL DO default(shared) &
!$omp private(lf) 		      	  &
!$omp private(index_i,index_j,m) 

do index_j=1,nj
	do index_i=1,ni+1
        do m=-2,2
            lf(m)=flux_xin(index_i-m+bfsize_i,index_j,index_eq)
        enddo

		h_xin(index_i,index_j,index_eq)=(2.d0*lf(-2) &
										 -13.d0*lf(-1) &
										 +47.d0*lf(0) &
										 +27.d0*lf(1) &
										 - 3.d0*lf(2))/60.d0 
	enddo
enddo
!$OMP end PARALLEL DO
enddo
end subroutine

subroutine scheme_upwind5_etap(flux_etap,ni,nj,bfsize_i,bfsize_j,neq,h_etap)
implicit none
integer,intent(in) :: ni,nj,bfsize_i,bfsize_j,neq
double precision,intent(in),dimension(nj+2*bfsize_i,ni,neq) :: flux_etap
double precision,intent(out),dimension(nj+1,ni,neq) :: h_etap
integer :: index_i,index_j,index_eq,m
double precision :: lf(-2:2)

do index_eq=1,neq
!$OMP PARALLEL DO default(shared) &
!$omp private(lf) 		      	  &
!$omp private(index_i,index_j,m) 

do index_i=1,ni
	do index_j=1,nj+1
        do m=-2,2
            lf(m)=flux_etap(index_j+m+bfsize_j-1,index_i,index_eq)
        enddo

		h_etap(index_j,index_i,index_eq)= (2.d0*lf(-2) &
										 -13.d0*lf(-1) &
										 +47.d0*lf(0) &
										 +27.d0*lf(1) &
										 - 3.d0*lf(2))/60.d0 
	enddo
enddo
!$OMP end PARALLEL DO
enddo
end subroutine

subroutine scheme_upwind5_etan(flux_etan,ni,nj,bfsize_i,bfsize_j,neq,h_etan)
implicit none
integer,intent(in) :: ni,nj,bfsize_i,bfsize_j,neq
double precision,intent(in),dimension(nj+2*bfsize_i,ni,neq) :: flux_etan
double precision,intent(out),dimension(nj+1,ni,neq) :: h_etan
integer :: index_i,index_j,index_eq,m
double precision :: lf(-2:2)

do index_eq=1,neq
!$OMP PARALLEL DO default(shared) &
!$omp private(lf) 		      	  &
!$omp private(index_i,index_j,m)

do index_i=1,ni
	do index_j=1,nj+1
        do m=-2,2
            lf(m)=flux_etan(index_j-m+bfsize_j,index_i,index_eq)
        enddo

		h_etan(index_j,index_i,index_eq)=(2.d0*lf(-2) &
										 -13.d0*lf(-1) &
										 +47.d0*lf(0) &
										 +27.d0*lf(1) &
										 - 3.d0*lf(2))/60.d0 
	enddo
enddo
!$OMP end PARALLEL DO
enddo
end subroutine
