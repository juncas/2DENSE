subroutine flux_splitting_LLF
use flow_general
implicit none

if(if_curve_grid) then
    call flux_splitting_LLF_curve
else 
    call flux_splitting_LLF_cartisian
endif

end subroutine

subroutine flux_splitting_LLF_cartisian
use flow_general
use var_index
implicit none
double precision :: rho,u,v,p,alpha
double precision :: max_eig_i_local,max_eig_j_local
integer :: index_i,index_j,index_block
integer :: bfsize_i,bfsize_j,ni,nj

do index_block=1,num_block
    ni = block_grid(index_block)%ni
    nj = block_grid(index_block)%nj
    bfsize_i=block_grid(index_block)%bfsize_i
    bfsize_j=block_grid(index_block)%bfsize_j
    !$OMP PARALLEL DO default(shared) &
    !$omp private(u,p,alpha,max_eig_i_local) &
    !$omp private(index_i,index_j)

    do index_j=1,nj
        do index_i=1,ni+2*bfsize_i
            u = block_grid(index_block)%prim(index_i,index_j+bfsize_j,prim_u_index)
            p = block_grid(index_block)%prim(index_i,index_j+bfsize_j,prim_p_index)
            alpha = block_grid(index_block)%prim(index_i,index_j+bfsize_j,prim_alpha_index)
           
            max_eig_i_local = abs(u)+alpha
           
            block_grid(index_block)%flux_xip(index_i,index_j,1)=0.5*((u+max_eig_i_local)*block_grid(index_block)%cons(index_i,index_j+bfsize_j,1))
            block_grid(index_block)%flux_xip(index_i,index_j,2)=0.5*((u+max_eig_i_local)*block_grid(index_block)%cons(index_i,index_j+bfsize_j,2)+p)
            block_grid(index_block)%flux_xip(index_i,index_j,3)=0.5*((u+max_eig_i_local)*block_grid(index_block)%cons(index_i,index_j+bfsize_j,3))
            block_grid(index_block)%flux_xip(index_i,index_j,4)=0.5*((u+max_eig_i_local)*block_grid(index_block)%cons(index_i,index_j+bfsize_j,4)+u*p)
            
            block_grid(index_block)%flux_xin(index_i,index_j,1)=0.5*((u-max_eig_i_local)*block_grid(index_block)%cons(index_i,index_j+bfsize_j,1))
            block_grid(index_block)%flux_xin(index_i,index_j,2)=0.5*((u-max_eig_i_local)*block_grid(index_block)%cons(index_i,index_j+bfsize_j,2)+p)
            block_grid(index_block)%flux_xin(index_i,index_j,3)=0.5*((u-max_eig_i_local)*block_grid(index_block)%cons(index_i,index_j+bfsize_j,3))
            block_grid(index_block)%flux_xin(index_i,index_j,4)=0.5*((u-max_eig_i_local)*block_grid(index_block)%cons(index_i,index_j+bfsize_j,4)+u*p)
        enddo
    enddo
    !$OMP END PARALLEL DO

    !$OMP PARALLEL DO default(shared) &
    !$omp private(v,p,alpha,max_eig_j_local) &
    !$omp private(index_i,index_j)

    do index_i=1,ni
        do index_j=1,nj+2*bfsize_j
            v = block_grid(index_block)%prim(index_i+bfsize_i,index_j,prim_v_index)
            p = block_grid(index_block)%prim(index_i+bfsize_i,index_j,prim_p_index)
            alpha = block_grid(index_block)%prim(index_i+bfsize_i,index_j,prim_alpha_index)
            
            max_eig_j_local = abs(v)+alpha
            
            block_grid(index_block)%flux_etap(index_j,index_i,1)=0.5*((v+max_eig_j_local)*block_grid(index_block)%cons(index_i+bfsize_i,index_j,1))
            block_grid(index_block)%flux_etap(index_j,index_i,2)=0.5*((v+max_eig_j_local)*block_grid(index_block)%cons(index_i+bfsize_i,index_j,2))
            block_grid(index_block)%flux_etap(index_j,index_i,3)=0.5*((v+max_eig_j_local)*block_grid(index_block)%cons(index_i+bfsize_i,index_j,3)+p)
            block_grid(index_block)%flux_etap(index_j,index_i,4)=0.5*((v+max_eig_j_local)*block_grid(index_block)%cons(index_i+bfsize_i,index_j,4)+v*p)
            
            block_grid(index_block)%flux_etan(index_j,index_i,1)=0.5*((v-max_eig_j_local)*block_grid(index_block)%cons(index_i+bfsize_i,index_j,1))
            block_grid(index_block)%flux_etan(index_j,index_i,2)=0.5*((v-max_eig_j_local)*block_grid(index_block)%cons(index_i+bfsize_i,index_j,2))
            block_grid(index_block)%flux_etan(index_j,index_i,3)=0.5*((v-max_eig_j_local)*block_grid(index_block)%cons(index_i+bfsize_i,index_j,3)+p)
            block_grid(index_block)%flux_etan(index_j,index_i,4)=0.5*((v-max_eig_j_local)*block_grid(index_block)%cons(index_i+bfsize_i,index_j,4)+v*p)
        enddo
    enddo
    !$OMP END PARALLEL DO
enddo
end subroutine


subroutine flux_splitting_LLF_curve
use flow_general
use var_index
implicit none
double precision :: rho,u,v,p,alpha
double precision :: max_eig_i_local,max_eig_j_local
integer :: index_i,index_j,index_block
integer :: bfsize_i,bfsize_j,ni,nj

end subroutine


