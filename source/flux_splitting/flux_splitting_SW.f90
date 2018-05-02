subroutine flux_splitting_SW
use flow_general
implicit none

if(if_curve_grid) then
    call flux_splitting_SW_curve
else 
    call flux_splitting_SW_cartisian
endif


end subroutine

subroutine flux_splitting_SW_cartisian
use flow_general
use var_index
implicit none
double precision :: rho,u,v,p,alpha
double precision :: lambda(3),lambda_pm(3),h,uu2,wii
double precision,parameter :: epsl = 0.125d0
integer :: index_i,index_j,index_block
integer :: bfsize_i,bfsize_j,ni,nj

do index_block=1,num_block
    ni = block_grid(index_block)%ni
    nj = block_grid(index_block)%nj
    bfsize_i=block_grid(index_block)%bfsize_i
    bfsize_j=block_grid(index_block)%bfsize_j

    !$OMP PARALLEL DO default(shared) &
    !$omp private(rho,u,v,p,alpha,uu2,wii) &
    !$omp private(lambda,lambda_pm) &
    !$omp private(index_i,index_j)
    do index_j=1,nj
        do index_i=1,ni+2*bfsize_i            
        	rho = block_grid(index_block)%prim(index_i,index_j+bfsize_j,prim_rho_index)
            u = block_grid(index_block)%prim(index_i,index_j+bfsize_j,prim_u_index)
            v = block_grid(index_block)%prim(index_i,index_j+bfsize_j,prim_v_index)
            p = block_grid(index_block)%prim(index_i,index_j+bfsize_j,prim_p_index)
            alpha = block_grid(index_block)%prim(index_i,index_j+bfsize_j,prim_alpha_index)

            uu2 = u**2+v**2
            lambda(1) = u
            lambda(2) = u-alpha
            lambda(3) = u+alpha

            lambda_pm(1) = (lambda(1)+sqrt(lambda(1)**2+epsl**2))/2.0
            lambda_pm(2) = (lambda(2)+sqrt(lambda(2)**2+epsl**2))/2.0
            lambda_pm(3) = (lambda(3)+sqrt(lambda(3)**2+epsl**2))/2.0
			
			wii=(3.d0-gamma)*(lambda_pm(2)+lambda_pm(3))*alpha**2/(2.d0*(gamma-1.d0))

			block_grid(index_block)%flux_xip(index_i,index_j,1) =(lambda_pm(1)*2.d0*(gamma-1.d0) &
			                                                     +lambda_pm(2) &
			                                                     +lambda_pm(3))*rho/(2.d0*gamma)

			block_grid(index_block)%flux_xip(index_i,index_j,2) =(lambda_pm(1)*2.d0*(gamma-1.d0)*u &
												                 +lambda_pm(2)*(u-alpha) &
												                 +lambda_pm(3)*(u+alpha))*rho/(2.d0*gamma)

			block_grid(index_block)%flux_xip(index_i,index_j,3) =(lambda_pm(1)*2.d0*(gamma-1.d0)*v &
			                                                     +lambda_pm(2)*v &
												                 +lambda_pm(3)*v)*rho/(2.d0*gamma)

			block_grid(index_block)%flux_xip(index_i,index_j,4) =(lambda_pm(1)*(gamma-1.d0)*uu2 &
												                 +lambda_pm(2)*0.5d0*((u-alpha)**2+v**2) &
												                 +lambda_pm(3)*0.5d0*((u+alpha)**2+v**2) &
												                 +wii)*rho/(2.d0*gamma)
           
            lambda_pm(1) = (lambda(1)-sqrt(lambda(1)**2+epsl**2))/2.0
            lambda_pm(2) = (lambda(2)-sqrt(lambda(2)**2+epsl**2))/2.0
            lambda_pm(3) = (lambda(3)-sqrt(lambda(3)**2+epsl**2))/2.0

			wii=(3.d0-gamma)*(lambda_pm(2)+lambda_pm(3))*alpha**2/(2.d0*(gamma-1.d0))

			block_grid(index_block)%flux_xin(index_i,index_j,1) =(lambda_pm(1)*2.d0*(gamma-1.d0) &
			                                                     +lambda_pm(2) &
			                                                     +lambda_pm(3))*rho/(2.d0*gamma)

			block_grid(index_block)%flux_xin(index_i,index_j,2) =(lambda_pm(1)*2.d0*(gamma-1.d0)*u &
												                 +lambda_pm(2)*(u-alpha) &
												                 +lambda_pm(3)*(u+alpha))*rho/(2.d0*gamma)

			block_grid(index_block)%flux_xin(index_i,index_j,3) =(lambda_pm(1)*2.d0*(gamma-1.d0)*v &
			                                                     +lambda_pm(2)*v &
												                 +lambda_pm(3)*v)*rho/(2.d0*gamma)

			block_grid(index_block)%flux_xin(index_i,index_j,4) =(lambda_pm(1)*(gamma-1.d0)*uu2 &
												                 +lambda_pm(2)*0.5d0*((u-alpha)**2+v**2) &
												                 +lambda_pm(3)*0.5d0*((u+alpha)**2+v**2) &
												                 +wii)*rho/(2.d0*gamma)
           
 		enddo
    enddo
    !$OMP END PARALLEL DO

    !$OMP PARALLEL DO default(shared) &
    !$omp private(rho,u,v,p,alpha,uu2,wii) &
    !$omp private(lambda,lambda_pm) &
    !$omp private(index_i,index_j)

    do index_i=1,ni
        do index_j=1,nj+2*bfsize_j
            rho = block_grid(index_block)%prim(index_i+bfsize_i,index_j,prim_rho_index)
            u = block_grid(index_block)%prim(index_i+bfsize_i,index_j,prim_u_index)
            v = block_grid(index_block)%prim(index_i+bfsize_i,index_j,prim_v_index)
            p = block_grid(index_block)%prim(index_i+bfsize_i,index_j,prim_p_index)
            alpha = block_grid(index_block)%prim(index_i+bfsize_i,index_j,prim_alpha_index)

            uu2 = u**2+v**2
            lambda(1) = v
            lambda(2) = v-alpha
            lambda(3) = v+alpha

            lambda_pm(1) = (lambda(1)+sqrt(lambda(1)**2+epsl**2))/2.0
            lambda_pm(2) = (lambda(2)+sqrt(lambda(2)**2+epsl**2))/2.0
            lambda_pm(3) = (lambda(3)+sqrt(lambda(3)**2+epsl**2))/2.0
			
			wii=(3.d0-gamma)*(lambda_pm(2)+lambda_pm(3))*alpha**2/(2.d0*(gamma-1.d0))

			block_grid(index_block)%flux_etap(index_j,index_i,1) =(lambda_pm(1)*2.d0*(gamma-1.d0) &
			                                                      +lambda_pm(2) &
			                                                      +lambda_pm(3))*rho/(2.d0*gamma)

			block_grid(index_block)%flux_etap(index_j,index_i,2) =(lambda_pm(1)*2.d0*(gamma-1.d0)*u &
												                  +lambda_pm(2)*u &
												                  +lambda_pm(3)*u)*rho/(2.d0*gamma)

			block_grid(index_block)%flux_etap(index_j,index_i,3) =(lambda_pm(1)*2.d0*(gamma-1.d0)*v &
			                                                      +lambda_pm(2)*(v-alpha) &
												                  +lambda_pm(3)*(v+alpha))*rho/(2.d0*gamma)

			block_grid(index_block)%flux_etap(index_j,index_i,4) =(lambda_pm(1)*(gamma-1.d0)*uu2 &
												                  +lambda_pm(2)*0.5d0*((v-alpha)**2+u**2) &
												                  +lambda_pm(3)*0.5d0*((v+alpha)**2+u**2) &
												                  +wii)*rho/(2.d0*gamma)
           

            lambda_pm(1) = (lambda(1)-sqrt(lambda(1)**2+epsl**2))/2.0
            lambda_pm(2) = (lambda(2)-sqrt(lambda(2)**2+epsl**2))/2.0
            lambda_pm(3) = (lambda(3)-sqrt(lambda(3)**2+epsl**2))/2.0
			
			wii=(3.d0-gamma)*(lambda_pm(2)+lambda_pm(3))*alpha**2/(2.d0*(gamma-1.d0))

			block_grid(index_block)%flux_etan(index_j,index_i,1) =(lambda_pm(1)*2.d0*(gamma-1.d0) &
			                                                      +lambda_pm(2) &
			                                                      +lambda_pm(3))*rho/(2.d0*gamma)

			block_grid(index_block)%flux_etan(index_j,index_i,2) =(lambda_pm(1)*2.d0*(gamma-1.d0)*u &
												                  +lambda_pm(2)*u &
												                  +lambda_pm(3)*u)*rho/(2.d0*gamma)

			block_grid(index_block)%flux_etan(index_j,index_i,3) =(lambda_pm(1)*2.d0*(gamma-1.d0)*v &
			                                                      +lambda_pm(2)*(v-alpha) &
												                  +lambda_pm(3)*(v+alpha))*rho/(2.d0*gamma)

			block_grid(index_block)%flux_etan(index_j,index_i,4) =(lambda_pm(1)*(gamma-1.d0)*uu2 &
												                  +lambda_pm(2)*0.5d0*((v-alpha)**2+u**2) &
												                  +lambda_pm(3)*0.5d0*((v+alpha)**2+u**2) &
												                  +wii)*rho/(2.d0*gamma)
           

		enddo
    enddo
    !$OMP END PARALLEL DO
enddo

end subroutine flux_splitting_SW_cartisian

subroutine flux_splitting_SW_curve
use flow_general
implicit none

end subroutine flux_splitting_SW_curve
