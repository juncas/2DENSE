subroutine weno_character_5th_xi_p(prim,flux_xip,ni,nj,bfsize_i,bfsize_j,neq,h_xip)
use flow_general,only:gamma,delta
use var_index
use weno
implicit none
integer,intent(in) :: ni,nj,bfsize_i,bfsize_j,neq
double precision,intent(in),dimension(ni+2*bfsize_i,nj+2*bfsize_j,neq) :: prim
double precision,intent(in),dimension(ni+2*bfsize_i,nj,neq) :: flux_xip

double precision,intent(inout),dimension(ni+1,nj,neq) :: h_xip

double precision :: w0,w1,w2
double precision :: q30,q31,q32

double precision :: rhor,rhol
double precision :: ur,ul
double precision :: vr,vl
double precision :: pr,pl
double precision :: hr,hl
double precision :: h0,u0,v0,alpha0

double precision :: char_var(neq)
double precision,dimension(-2:2) :: lf
double precision ::r(neq,neq),l(neq,neq)
double precision :: h_xip_temp
integer :: index_i,index_j,index_eq,m,k


!$OMP PARALLEL DO default(shared) &
!$omp private(w0,w1,w2) &
!$omp private(q30,q31,q32) &
!$omp private(lf,r,l) &
!$omp private(rhor,rhol) &
!$omp private(ur,ul) &
!$omp private(vr,vl) &
!$omp private(pr,pl) &
!$omp private(hr,hl) &
!$omp private(h0,u0,v0,alpha0) &
!$omp private(char_var) &
!$omp private(h_xip_temp) &
!$omp private(index_i,index_j,index_eq,m,k)
do index_j=1,nj
	do index_i=1,ni+1
		rhol = prim(index_i+bfsize_i-1,index_j+bfsize_j,prim_rho_index)
		rhor = prim(index_i+bfsize_i  ,index_j+bfsize_j,prim_rho_index)

		ul   = prim(index_i+bfsize_i-1,index_j+bfsize_j,prim_u_index)
		ur   = prim(index_i+bfsize_i  ,index_j+bfsize_j,prim_u_index)

		vl   = prim(index_i+bfsize_i-1,index_j+bfsize_j,prim_v_index)
		vr   = prim(index_i+bfsize_i  ,index_j+bfsize_j,prim_v_index)

		pl   = prim(index_i+bfsize_i-1,index_j+bfsize_j,prim_p_index)   
		pr   = prim(index_i+bfsize_i  ,index_j+bfsize_j,prim_p_index)  

		hl   = pl/(gamma-1)/rhol+0.5*(ul*ul+vl*vl)+pl/rhol
		hr   = pr/(gamma-1)/rhor+0.5*(ur*ur+vr*vr)+pr/rhor

		u0=(sqrt(rhor)*ur+sqrt(rhol)*ul)/(sqrt(rhor)+sqrt(rhol))
		v0=(sqrt(rhor)*vr+sqrt(rhol)*vl)/(sqrt(rhor)+sqrt(rhol))
		h0=(sqrt(rhor)*hr+sqrt(rhol)*hl)/(sqrt(rhor)+sqrt(rhol))
		alpha0=sqrt((gamma-1.0)*(h0-0.5d0*(u0**2+v0**2)))

		call get_l_r_vec_xi(u0,v0,h0,alpha0,l,r,neq)

		do index_eq=1,neq
			do m=-2,2       
					lf(m)= l(1,index_eq)*flux_xip(index_i+m+bfsize_i-1,index_j,1) &          
						  +l(2,index_eq)*flux_xip(index_i+m+bfsize_i-1,index_j,2) &          
				          +l(3,index_eq)*flux_xip(index_i+m+bfsize_i-1,index_j,3) &          
				          +l(4,index_eq)*flux_xip(index_i+m+bfsize_i-1,index_j,4)           
			enddo
			call weno5z_stencil(lf,char_var(index_eq))
		enddo

		do index_eq=1,neq
			h_xip(index_i,index_j,index_eq) = r(1,index_eq)*char_var(1) &
			                                 +r(2,index_eq)*char_var(2) &
			                                 +r(3,index_eq)*char_var(3) &
			                                 +r(4,index_eq)*char_var(4)
		enddo
	enddo
enddo
!$OMP END PARALLEL DO

end subroutine

subroutine weno_character_5th_xi_n(prim,flux_xin,ni,nj,bfsize_i,bfsize_j,neq,h_xin)
use flow_general,only:gamma,delta
use var_index
use weno
implicit none
integer,intent(in) :: ni,nj,bfsize_i,bfsize_j,neq
double precision,intent(in),dimension(ni+2*bfsize_i,nj+2*bfsize_j,neq) :: prim
double precision,intent(in),dimension(ni+2*bfsize_i,nj,neq) :: flux_xin

double precision,intent(inout),dimension(ni+1,nj,neq) :: h_xin

double precision :: w0,w1,w2
double precision :: q30,q31,q32

double precision :: rhor,rhol
double precision :: ur,ul
double precision :: vr,vl
double precision :: pr,pl
double precision :: hr,hl
double precision :: h0,u0,v0,alpha0

double precision :: char_var(neq)
double precision,dimension(-2:2) :: lf
double precision ::r(neq,neq),l(neq,neq)
double precision :: h_xin_temp
integer :: index_i,index_j,index_eq,m,k


!$OMP PARALLEL DO default(shared) &
!$omp private(w0,w1,w2) &
!$omp private(q30,q31,q32) &
!$omp private(lf,r,l) &
!$omp private(rhor,rhol) &
!$omp private(ur,ul) &
!$omp private(vr,vl) &
!$omp private(pr,pl) &
!$omp private(hr,hl) &
!$omp private(h0,u0,v0,alpha0) &
!$omp private(h_xin_temp) &
!$omp private(char_var) &
!$omp private(index_i,index_j,index_eq,m,k)
do index_j=1,nj
	do index_i=1,ni+1
		rhol = prim(index_i+bfsize_i-1,index_j+bfsize_j,prim_rho_index)
		rhor = prim(index_i+bfsize_i  ,index_j+bfsize_j,prim_rho_index)

		ul   = prim(index_i+bfsize_i-1,index_j+bfsize_j,prim_u_index)
		ur   = prim(index_i+bfsize_i  ,index_j+bfsize_j,prim_u_index)

		vl   = prim(index_i+bfsize_i-1,index_j+bfsize_j,prim_v_index)
		vr   = prim(index_i+bfsize_i  ,index_j+bfsize_j,prim_v_index)

		pl   = prim(index_i+bfsize_i-1,index_j+bfsize_j,prim_p_index)   
		pr   = prim(index_i+bfsize_i  ,index_j+bfsize_j,prim_p_index)  

		hl   = pl/(gamma-1)/rhol+0.5*(ul*ul+vl*vl)+pl/rhol
		hr   = pr/(gamma-1)/rhor+0.5*(ur*ur+vr*vr)+pr/rhor

		u0=(sqrt(rhor)*ur+sqrt(rhol)*ul)/(sqrt(rhor)+sqrt(rhol))
		v0=(sqrt(rhor)*vr+sqrt(rhol)*vl)/(sqrt(rhor)+sqrt(rhol))
		h0=(sqrt(rhor)*hr+sqrt(rhol)*hl)/(sqrt(rhor)+sqrt(rhol))
		alpha0=sqrt((gamma-1)*(h0-0.5d0*(u0**2+v0**2)))

		call get_l_r_vec_xi(u0,v0,h0,alpha0,l,r,neq)

		do index_eq=1,neq
			do m=-2,2     
					lf(m)= l(1,index_eq)*flux_xin(index_i+bfsize_i-m,index_j,1) &
					      +l(2,index_eq)*flux_xin(index_i+bfsize_i-m,index_j,2) &
						  +l(3,index_eq)*flux_xin(index_i+bfsize_i-m,index_j,3) &
						  +l(4,index_eq)*flux_xin(index_i+bfsize_i-m,index_j,4)           
				
			enddo
			call weno5z_stencil(lf,char_var(index_eq))
		enddo

		do index_eq=1,neq
			h_xin(index_i,index_j,index_eq) = r(1,index_eq)*char_var(1) &
			                                 +r(2,index_eq)*char_var(2) &
			                                 +r(3,index_eq)*char_var(3) &
			                                 +r(4,index_eq)*char_var(4)
		enddo
	enddo
enddo
!$OMP END PARALLEL DO

end subroutine

subroutine weno_character_5th_eta_p(prim,flux_etap,ni,nj,bfsize_i,bfsize_j,neq,h_etap)
use flow_general,only:gamma,delta
use var_index
use weno
implicit none
integer,intent(in) :: ni,nj,bfsize_i,bfsize_j,neq
double precision,intent(in),dimension(ni+2*bfsize_i,nj+2*bfsize_j,neq) :: prim
double precision,intent(in),dimension(nj+2*bfsize_j,ni,neq) :: flux_etap

double precision,intent(inout),dimension(nj+1,ni,neq) :: h_etap

double precision :: w0,w1,w2
double precision :: q30,q31,q32

double precision :: rhor,rhol
double precision :: ur,ul
double precision :: vr,vl
double precision :: pr,pl
double precision :: hr,hl
double precision :: h0,u0,v0,alpha0

double precision :: char_var(neq)
double precision,dimension(-2:2) :: lf
double precision ::r(neq,neq),l(neq,neq)
double precision :: h_etap_temp
integer :: index_i,index_j,index_eq,m,k


!$OMP PARALLEL DO default(shared) &
!$omp private(w0,w1,w2) &
!$omp private(q30,q31,q32) &
!$omp private(lf,r,l) &
!$omp private(rhor,rhol) &
!$omp private(ur,ul) &
!$omp private(vr,vl) &
!$omp private(pr,pl) &
!$omp private(hr,hl) &
!$omp private(h0,u0,v0,alpha0) &
!$omp private(h_etap_temp) &
!$omp private(index_i,index_j,index_eq,m,k)
do index_i=1,ni
	do index_j=1,nj+1     
		rhol = prim(index_i+bfsize_i,index_j+bfsize_j-1,prim_rho_index)
		rhor = prim(index_i+bfsize_i,index_j+bfsize_j  ,prim_rho_index)

		ul   = prim(index_i+bfsize_i,index_j+bfsize_j-1,prim_u_index)
		ur   = prim(index_i+bfsize_i,index_j+bfsize_j  ,prim_u_index)

		vl   = prim(index_i+bfsize_i,index_j+bfsize_j-1,prim_v_index)
		vr   = prim(index_i+bfsize_i,index_j+bfsize_j  ,prim_v_index)

		pl   = prim(index_i+bfsize_i,index_j+bfsize_j-1,prim_p_index)   
		pr   = prim(index_i+bfsize_i,index_j+bfsize_j  ,prim_p_index)  

		hl=pl/(gamma-1)/rhol+0.5*(ul*ul+vl*vl)+pl/rhol
		hr=pr/(gamma-1)/rhor+0.5*(ur*ur+vr*vr)+pr/rhor

		u0=(sqrt(rhor)*ur+sqrt(rhol)*ul)/(sqrt(rhor)+sqrt(rhol))
		v0=(sqrt(rhor)*vr+sqrt(rhol)*vl)/(sqrt(rhor)+sqrt(rhol))
		h0=(sqrt(rhor)*hr+sqrt(rhol)*hl)/(sqrt(rhor)+sqrt(rhol))
		alpha0=sqrt((gamma-1)*(h0-0.5d0*(u0**2+v0**2)))

		call get_l_r_vec_eta(u0,v0,h0,alpha0,l,r,neq)

		do index_eq=1,neq
			do m=-2,2        
				lf(m)= l(1,index_eq)*flux_etap(index_j+bfsize_j-1+m,index_i,1) &           
					  +l(2,index_eq)*flux_etap(index_j+bfsize_j-1+m,index_i,2) &         
			          +l(3,index_eq)*flux_etap(index_j+bfsize_j-1+m,index_i,3) &         
			          +l(4,index_eq)*flux_etap(index_j+bfsize_j-1+m,index_i,4)          
			enddo
			call weno5z_stencil(lf,char_var(index_eq))
		enddo

		do index_eq=1,neq
		h_etap(index_j,index_i,index_eq) = r(1,index_eq)*char_var(1) &
			                              +r(2,index_eq)*char_var(2) &
			                              +r(3,index_eq)*char_var(3) &
			                              +r(4,index_eq)*char_var(4)
		enddo
	enddo
enddo
!$OMP END PARALLEL DO

end subroutine

subroutine weno_character_5th_eta_n(prim,flux_etan,ni,nj,bfsize_i,bfsize_j,neq,h_etan)
use flow_general,only:gamma,delta
use var_index
use weno
implicit none
integer,intent(in) :: ni,nj,bfsize_i,bfsize_j,neq
double precision,intent(in),dimension(ni+2*bfsize_i,nj+2*bfsize_j,neq) :: prim
double precision,intent(in),dimension(nj+2*bfsize_j,ni,neq) :: flux_etan

double precision,intent(inout),dimension(nj+1,ni,neq) :: h_etan

double precision :: w0,w1,w2
double precision :: q30,q31,q32

double precision :: rhor,rhol
double precision :: ur,ul
double precision :: vr,vl
double precision :: pr,pl
double precision :: hr,hl
double precision :: h0,u0,v0,alpha0

double precision :: char_var(neq)
double precision,dimension(-2:2) :: lf
double precision ::r(neq,neq),l(neq,neq)
double precision :: h_etan_temp
integer :: index_i,index_j,index_eq,m,k


!$OMP PARALLEL DO default(shared) &
!$omp private(w0,w1,w2) &
!$omp private(q30,q31,q32) &
!$omp private(lf,r,l) &
!$omp private(rhor,rhol) &
!$omp private(ur,ul) &
!$omp private(vr,vl) &
!$omp private(pr,pl) &
!$omp private(hr,hl) &
!$omp private(h0,u0,v0,alpha0) &
!$omp private(h_etan_temp) &
!$omp private(index_i,index_j,index_eq,m,k)
do index_i=1,ni
	do index_j=1,nj+1
		rhol = prim(index_i+bfsize_i,index_j+bfsize_j-1,prim_rho_index)
		rhor = prim(index_i+bfsize_i,index_j+bfsize_j  ,prim_rho_index)

		ul   = prim(index_i+bfsize_i,index_j+bfsize_j-1,prim_u_index)
		ur   = prim(index_i+bfsize_i,index_j+bfsize_j  ,prim_u_index)

		vl   = prim(index_i+bfsize_i,index_j+bfsize_j-1,prim_v_index)
		vr   = prim(index_i+bfsize_i,index_j+bfsize_j  ,prim_v_index)

		pl   = prim(index_i+bfsize_i,index_j+bfsize_j-1,prim_p_index)   
		pr   = prim(index_i+bfsize_i,index_j+bfsize_j  ,prim_p_index)  

		hl=pl/(gamma-1)/rhol+0.5*(ul*ul+vl*vl)+pl/rhol
		hr=pr/(gamma-1)/rhor+0.5*(ur*ur+vr*vr)+pr/rhor

		u0=(sqrt(rhor)*ur+sqrt(rhol)*ul)/(sqrt(rhor)+sqrt(rhol))
		v0=(sqrt(rhor)*vr+sqrt(rhol)*vl)/(sqrt(rhor)+sqrt(rhol))
		h0=(sqrt(rhor)*hr+sqrt(rhol)*hl)/(sqrt(rhor)+sqrt(rhol))
		alpha0=sqrt((gamma-1)*(h0-0.5d0*(u0**2+v0**2)))

		call get_l_r_vec_eta(u0,v0,h0,alpha0,l,r,neq)

		do index_eq=1,neq
			do m=-2,2        
				lf(m)= l(1,index_eq)*flux_etan(index_j+bfsize_j-m,index_i,1) &        
			          +l(2,index_eq)*flux_etan(index_j+bfsize_j-m,index_i,2) &   
			          +l(3,index_eq)*flux_etan(index_j+bfsize_j-m,index_i,3) &     
			          +l(4,index_eq)*flux_etan(index_j+bfsize_j-m,index_i,4)          
			enddo
			
			call weno5z_stencil(lf,char_var(index_eq))
		enddo

		do index_eq=1,neq
		h_etan(index_j,index_i,index_eq)= r(1,index_eq)*char_var(1) &
			                             +r(2,index_eq)*char_var(2) &
			                             +r(3,index_eq)*char_var(3) &
			                             +r(4,index_eq)*char_var(4)
		enddo
	enddo
enddo
!$OMP END PARALLEL DO

end subroutine
