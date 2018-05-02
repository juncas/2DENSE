subroutine get_theta_xi_p(theta_xi_p,var,ni,nj,bfsize_i,bfsize_j)
use weno
implicit none
!inputs
integer,intent(in) :: ni,nj,bfsize_i,bfsize_j
double precision,intent(in),dimension(ni+2*bfsize_i,nj) :: var

!outputs
double precision,intent(out),dimension(4,ni+1,nj) :: theta_xi_p

!local variables
double precision :: lf(-2:2)
integer :: index_i,index_j,m

!$OMP PARALLEL DO default(shared) &
!$omp private(lf) &
!$omp private(index_i,index_j,m)
do index_j=1,nj
	do index_i=1,ni+1
        do m=-2,2
            lf(m)=var(index_i+m+bfsize_i-1,index_j)
        enddo

        call weno5z_stencil_omega(lf,                            &
        						  theta_xi_p(2,index_i,index_j), &
        						  theta_xi_p(3,index_i,index_j), &
        						  theta_xi_p(4,index_i,index_j), &
        						  theta_xi_p(1,index_i,index_j))
	enddo
enddo
!$OMP END PARALLEL DO
end subroutine

subroutine get_theta_xi_n(theta_xi_n,var,ni,nj,bfsize_i,bfsize_j)
use weno
implicit none
!inputs
integer,intent(in) :: ni,nj,bfsize_i,bfsize_j
double precision,intent(in),dimension(ni+2*bfsize_i,nj) :: var

!outputs
double precision,intent(out),dimension(4,ni+1,nj) :: theta_xi_n

!local variables
double precision :: lf(-2:2)
integer :: index_i,index_j,m

!$OMP PARALLEL DO default(shared) &
!$omp private(lf) &
!$omp private(index_i,index_j,m)
do index_j=1,nj
	do index_i=1,ni+1
        do m=-2,2
            lf(m)=var(index_i-m+bfsize_i,index_j)
        enddo
        call weno5z_stencil_omega(lf,                            &
        						  theta_xi_n(2,index_i,index_j), &
        						  theta_xi_n(3,index_i,index_j), &
        						  theta_xi_n(4,index_i,index_j), &
        						  theta_xi_n(1,index_i,index_j))
	enddo
enddo
!$OMP END PARALLEL DO
end subroutine

subroutine weno_semi_chart_5th_xi_p(prim,flux_xip,theta_xi_p,ni,nj,bfsize_i,bfsize_j,neq,h_xip)
use flow_general,only:gamma,delta
use var_index
use weno
implicit none
integer,intent(in) :: ni,nj,bfsize_i,bfsize_j,neq
double precision,intent(in),dimension(ni+2*bfsize_i,nj+2*bfsize_j,neq) :: prim
double precision,intent(in),dimension(ni+2*bfsize_i,nj,neq) :: flux_xip
double precision,intent(in),dimension(4,ni+1,nj) :: theta_xi_p

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
		if(theta_xi_p(1,index_i,index_j)>delta) then
	        w0=theta_xi_p(2,index_i,index_j)
	        w1=theta_xi_p(3,index_i,index_j)
	        w2=theta_xi_p(4,index_i,index_j)

			do index_eq=1,neq
		        do m=-2,2
		            lf(m)=flux_xip(index_i+m+bfsize_i-1,index_j,index_eq)
		        enddo
		        q30=a300*lf(-2)+a301*lf(-1)+a302*lf(0)
		        q31=a310*lf(-1)+a311*lf(0) +a312*lf(1)
		        q32=a320*lf(0) +a321*lf(1) +a322*lf(2)

		        h_xip(index_i,index_j,index_eq)=w0*q30+w1*q31+w2*q32
		    enddo
		else
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
		endif
	enddo
enddo
!$OMP END PARALLEL DO

end subroutine

subroutine weno_semi_chart_5th_xi_n(prim,flux_xin,theta_xi_n,ni,nj,bfsize_i,bfsize_j,neq,h_xin)
use flow_general,only:gamma,delta
use var_index
use weno
implicit none
integer,intent(in) :: ni,nj,bfsize_i,bfsize_j,neq
double precision,intent(in),dimension(ni+2*bfsize_i,nj+2*bfsize_j,neq) :: prim
double precision,intent(in),dimension(ni+2*bfsize_i,nj,neq) :: flux_xin
double precision,intent(in),dimension(4,ni+1,nj) :: theta_xi_n

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
		if(theta_xi_n(1,index_i,index_j)>delta) then
	        w0=theta_xi_n(2,index_i,index_j)
	        w1=theta_xi_n(3,index_i,index_j)
	        w2=theta_xi_n(4,index_i,index_j)

			do index_eq=1,neq
		        do m=-2,2
		            lf(m)=flux_xin(index_i+bfsize_i-m,index_j,index_eq)
		        enddo 
		        q30=a300*lf(-2)+a301*lf(-1)+a302*lf(0)
		        q31=a310*lf(-1)+a311*lf(0) +a312*lf(1)
		        q32=a320*lf(0) +a321*lf(1) +a322*lf(2)

		        h_xin(index_i,index_j,index_eq)=w0*q30+w1*q31+w2*q32
		    enddo
		else
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
		endif
	enddo
enddo
!$OMP END PARALLEL DO

end subroutine

subroutine get_l_r_vec_xi(u0,v0,h0,alpha0,l,r,neq)
use flow_general,only:gamma
implicit none
integer,intent(in) :: neq
double precision,intent(in) :: h0,u0,v0,alpha0
double precision,intent(out) ::r(neq,neq),l(neq,neq)
double precision ::uu0,hh0

uu0 = 0.5*(u0**2+v0**2)
hh0 = alpha0**2/(gamma-1.d0)

	r(1,1)=1.d0
	r(2,1)=1.d0
	r(3,1)=0.d0
	r(4,1)=1.d0

	r(1,2)=u0-alpha0
	r(2,2)=u0
	r(3,2)=0.d0
	r(4,2)=u0+alpha0

	r(1,3)=v0
	r(2,3)=v0
	r(3,3)=1.0
	r(4,3)=v0

	r(1,4)=h0-alpha0*u0
	r(2,4)=uu0
	r(3,4)=v0
	r(4,4)=h0+alpha0*u0

	l(1,1)= 0.5*uu0/hh0+0.5*u0/alpha0
	l(2,1)=-0.5*u0/hh0-0.5/alpha0
	l(3,1)=-0.5*v0/hh0
	l(4,1)= 0.5/hh0

	l(1,2)=-uu0/hh0+1.0
	l(2,2)= u0/hh0
	l(3,2)= v0/hh0
	l(4,2)=-1.0/hh0

	l(1,3)=-v0
	l(2,3)=0.d0
	l(3,3)=1.d0
	l(4,3)=0.d0

	l(1,4)= 0.5*uu0/hh0-0.5*u0/alpha0
	l(2,4)=-0.5*u0/hh0+0.5/alpha0
	l(3,4)=-0.5*v0/hh0
	l(4,4)= 0.5/hh0

end subroutine

subroutine get_theta_eta_p(theta_eta_p,var,ni,nj,bfsize_i,bfsize_j)
use weno
implicit none
!inputs
integer,intent(in) :: ni,nj,bfsize_i,bfsize_j
double precision,intent(in),dimension(nj+2*bfsize_j,ni) :: var

!outputs
double precision,intent(out),dimension(4,nj+1,ni) :: theta_eta_p

!local variables
double precision :: lf(-2:2)
integer :: index_i,index_j,m

!$OMP PARALLEL DO default(shared) &
!$omp private(lf) 				  &
!$omp private(index_i,index_j,m)
do index_i=1,ni
	do index_j=1,nj+1
        do m=-2,2
            lf(m)=var(index_j+m+bfsize_j-1,index_i)
        enddo
        call weno5z_stencil_omega(lf,                             &
        						  theta_eta_p(2,index_j,index_i), &
        						  theta_eta_p(3,index_j,index_i), &
        						  theta_eta_p(4,index_j,index_i), &
        						  theta_eta_p(1,index_j,index_i))
	enddo
enddo
!$OMP END PARALLEL DO	

end subroutine

subroutine get_theta_eta_n(theta_eta_n,var,ni,nj,bfsize_i,bfsize_j)
use weno
implicit none
!inputs
integer,intent(in) :: ni,nj,bfsize_i,bfsize_j
double precision,intent(in),dimension(nj+2*bfsize_j,ni) :: var

!outputs
double precision,intent(out),dimension(4,nj+1,ni) :: theta_eta_n

!local variables
double precision :: lf(-2:2)
integer :: index_i,index_j,m

!$OMP PARALLEL DO default(shared) &
!$omp private(lf) &
!$omp private(index_i,index_j,m)
do index_i=1,ni
	do index_j=1,nj+1
        do m=-2,2
            lf(m)=var(index_j-m+bfsize_j,index_i)
        enddo
        call weno5z_stencil_omega(lf,                             &
        						  theta_eta_n(2,index_j,index_i), &
        						  theta_eta_n(3,index_j,index_i), &
        						  theta_eta_n(4,index_j,index_i), &
        						  theta_eta_n(1,index_j,index_i))
	enddo
enddo
!$OMP END PARALLEL DO
end subroutine

subroutine weno_semi_chart_5th_eta_p(prim,flux_etap,theta_eta_p,ni,nj,bfsize_i,bfsize_j,neq,h_etap)
use flow_general,only:gamma,delta
use var_index
use weno
implicit none
integer,intent(in) :: ni,nj,bfsize_i,bfsize_j,neq
double precision,intent(in),dimension(ni+2*bfsize_i,nj+2*bfsize_j,neq) :: prim
double precision,intent(in),dimension(nj+2*bfsize_j,ni,neq) :: flux_etap
double precision,intent(in),dimension(4,nj+1,ni) :: theta_eta_p

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

		if(theta_eta_p(1,index_j,index_i)>delta) then

	        w0=theta_eta_p(2,index_j,index_i)
	        w1=theta_eta_p(3,index_j,index_i)
	        w2=theta_eta_p(4,index_j,index_i)

			do index_eq=1,neq
		        do m=-2,2
		            lf(m)=flux_etap(index_j+bfsize_j-1+m,index_i,index_eq)
		        enddo
		        q30=a300*lf(-2)+a301*lf(-1)+a302*lf(0)
		        q31=a310*lf(-1)+a311*lf(0) +a312*lf(1)
		        q32=a320*lf(0) +a321*lf(1) +a322*lf(2)

		        h_etap(index_j,index_i,index_eq)=w0*q30+w1*q31+w2*q32
		    enddo
		else           
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
		endif

	enddo
enddo
!$OMP END PARALLEL DO

end subroutine

subroutine weno_semi_chart_5th_eta_n(prim,flux_etan,theta_eta_n,ni,nj,bfsize_i,bfsize_j,neq,h_etan)
use flow_general,only:gamma,delta
use var_index
use weno
implicit none
integer,intent(in) :: ni,nj,bfsize_i,bfsize_j,neq
double precision,intent(in),dimension(ni+2*bfsize_i,nj+2*bfsize_j,neq) :: prim
double precision,intent(in),dimension(nj+2*bfsize_j,ni,neq) :: flux_etan
double precision,intent(in),dimension(4,nj+1,ni) :: theta_eta_n

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
		if(theta_eta_n(1,index_j,index_i)>delta) then
	        w0=theta_eta_n(2,index_j,index_i)
	        w1=theta_eta_n(3,index_j,index_i)
	        w2=theta_eta_n(4,index_j,index_i)

			do index_eq=1,neq
		        do m=-2,2
		            lf(m)=flux_etan(index_j+bfsize_j-m,index_i,index_eq)
		        enddo
		        q30=a300*lf(-2)+a301*lf(-1)+a302*lf(0)
		        q31=a310*lf(-1)+a311*lf(0) +a312*lf(1)
		        q32=a320*lf(0) +a321*lf(1) +a322*lf(2)

		        h_etan(index_j,index_i,index_eq)=w0*q30+w1*q31+w2*q32
		    enddo
		else
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
		endif
	enddo
enddo
!$OMP END PARALLEL DO

end subroutine

subroutine get_l_r_vec_eta(u0,v0,h0,alpha0,l,r,neq)
use flow_general,only:gamma
implicit none
integer,intent(in) :: neq
double precision,intent(in) :: h0,u0,v0,alpha0
double precision,intent(out) ::r(neq,neq),l(neq,neq)
double precision ::uu0,hh0

uu0 = 0.5*(u0**2+v0**2)
hh0 = alpha0**2/(gamma-1.d0)

	r(1,1)=1.d0
	r(2,1)=1.d0
	r(3,1)=0.d0
	r(4,1)=1.d0

	r(1,2)=u0
	r(2,2)=u0
	r(3,2)=1.0
	r(4,2)=u0

	r(1,3)=v0-alpha0
	r(2,3)=v0
	r(3,3)=0.d0
	r(4,3)=v0+alpha0

	r(1,4)=h0-alpha0*v0
	r(2,4)=uu0
	r(3,4)=u0
	r(4,4)=h0+alpha0*v0

	l(1,1)= 0.5*uu0/hh0+0.5*v0/alpha0
	l(2,1)=-0.5*u0/hh0
	l(3,1)=-0.5*v0/hh0-0.5/alpha0
	l(4,1)= 0.5/hh0

	l(1,2)=-uu0/hh0+1.0
	l(2,2)= u0/hh0
	l(3,2)= v0/hh0
	l(4,2)=-1.0/hh0

	l(1,3)=-u0
	l(2,3)= 1.d0
	l(3,3)= 0.d0
	l(4,3)= 0.d0

	l(1,4)= 0.5*uu0/hh0-0.5*v0/alpha0
	l(2,4)=-0.5*u0/hh0
	l(3,4)=-0.5*v0/hh0+0.5/alpha0
	l(4,4)= 0.5/hh0

end subroutine