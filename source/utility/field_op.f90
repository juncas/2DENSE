subroutine sum_field(field_sum,field0,field1,ni0,ni1,nj0,nj1)
implicit none
integer,intent(in) :: ni0,ni1 
integer,intent(in) :: nj0,nj1 
double precision,intent(in),dimension(ni0:ni1,nj0:nj1) :: field0,field1 
double precision,intent(out),dimension(ni0:ni1,nj0:nj1) :: field_sum
integer :: index_i,index_j


!$omp parallel do default(shared) &
!$omp private(index_i,index_j)	  
do index_j=nj0,nj1
	do index_i=ni0,ni1
		field_sum(index_i,index_j) &
	  = field0(index_i,index_j)    &
	   +field1(index_i,index_j)
	end do
end do
!$omp end parallel do

end subroutine sum_field

subroutine copy_field(dest,src,ni0,ni1,nj0,nj1)
implicit none
integer,intent(in) :: ni0,ni1 
integer,intent(in) :: nj0,nj1 
double precision,intent(in),dimension(ni0:ni1,nj0:nj1) :: src 
double precision,intent(out),dimension(ni0:ni1,nj0:nj1) :: dest
integer :: index_i,index_j


!$omp parallel do default(shared) &
!$omp private(index_i,index_j)	  
do index_j=nj0,nj1
	do index_i=ni0,ni1
		dest(index_i,index_j) &
	  = src(index_i,index_j)
	end do
end do
!$omp end parallel do

end subroutine copy_field
