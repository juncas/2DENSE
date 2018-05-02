subroutine scheme_weno5js_general
use flow_general
implicit none
integer :: index_block,index_eq

do index_block=1,num_block
	call scheme_weno5js_xip(block_grid(index_block)%flux_xip, &
						   block_grid(index_block)%ni,       &
						   block_grid(index_block)%nj,       &
						   block_grid(index_block)%bfsize_i, &
						   block_grid(index_block)%bfsize_j, &
						   neq,                              &
						   block_grid(index_block)%h_xip)

	call scheme_weno5js_xin(block_grid(index_block)%flux_xin, &
						   block_grid(index_block)%ni,       &
						   block_grid(index_block)%nj,       &
						   block_grid(index_block)%bfsize_i, &
						   block_grid(index_block)%bfsize_j, &
						   neq,                              &
						   block_grid(index_block)%h_xin)

	call scheme_weno5js_etap(block_grid(index_block)%flux_etap, &
						   block_grid(index_block)%ni,       &
						   block_grid(index_block)%nj,       &
						   block_grid(index_block)%bfsize_i, &
						   block_grid(index_block)%bfsize_j, &
						   neq,                              &
						   block_grid(index_block)%h_etap)

	call scheme_weno5js_etan(block_grid(index_block)%flux_etan, &
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
enddo
end subroutine scheme_weno5js_general

subroutine scheme_weno5js_xip(flux_xip,ni,nj,bfsize_i,bfsize_j,neq,h_xip)
use weno
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
		call weno5js_stencil(lf,h_xip(index_i,index_j,index_eq))
	enddo
enddo
!$OMP end PARALLEL DO
enddo
end subroutine

subroutine scheme_weno5js_xin(flux_xin,ni,nj,bfsize_i,bfsize_j,neq,h_xin)
use weno
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
		call weno5js_stencil(lf,h_xin(index_i,index_j,index_eq))
	enddo
enddo
!$OMP end PARALLEL DO
enddo
end subroutine


subroutine scheme_weno5js_etap(flux_etap,ni,nj,bfsize_i,bfsize_j,neq,h_etap)
use weno
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
		call weno5js_stencil(lf,h_etap(index_j,index_i,index_eq))
	enddo
enddo
!$OMP end PARALLEL DO
enddo
end subroutine

subroutine scheme_weno5js_etan(flux_etan,ni,nj,bfsize_i,bfsize_j,neq,h_etan)
use weno
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

		call weno5js_stencil(lf,h_etan(index_j,index_i,index_eq))
	enddo
enddo
!$OMP end PARALLEL DO
enddo
end subroutine

