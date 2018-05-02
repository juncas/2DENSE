module weno
implicit none
double precision,parameter :: a300 = 1.d0/3.d0
double precision,parameter :: a301 =-7.d0/6.d0
double precision,parameter :: a302 = 11.d0/6.d0
double precision,parameter :: a310 =-1.d0/6.d0
double precision,parameter :: a311 = 5.d0/6.d0
double precision,parameter :: a312 = 1.d0/3.d0
double precision,parameter :: a320 = 1.d0/3.d0
double precision,parameter :: a321 = 5.d0/6.d0
double precision,parameter :: a322 =-1.d0/6.d0
double precision,parameter :: c30  = 0.1d0
double precision,parameter :: c31  = 0.6d0
double precision,parameter :: c32  = 0.3d0

double precision,parameter :: a400=-1.d0/4.d0
double precision,parameter :: a401=13.d0/12.d0
double precision,parameter :: a402=-23.d0/12.d0
double precision,parameter :: a403=25.d0/12.d0
double precision,parameter :: a410=1.d0/12.d0
double precision,parameter :: a411=-5.d0/12.d0
double precision,parameter :: a412=13.d0/12.d0
double precision,parameter :: a413=3.d0/12.d0
double precision,parameter :: a420=-1.d0/12.d0
double precision,parameter :: a421=7.d0/12.d0
double precision,parameter :: a422=7.d0/12.d0
double precision,parameter :: a423=-1.d0/12.d0
double precision,parameter :: a430=1.d0/4.d0
double precision,parameter :: a431=13.d0/12.d0
double precision,parameter :: a432=-5.d0/12.d0
double precision,parameter :: a433=1.d0/12.d0
double precision,parameter :: c40=1.d0/35.d0
double precision,parameter :: c41=12.d0/35.d0
double precision,parameter :: c42=18.d0/35.d0
double precision,parameter :: c43=4.d0/35.d0

double precision,parameter :: ss   = 1.0e-6

interface weno5js_stencil_omega
	module procedure weno5js_stencil_omega_only
	module procedure weno5js_stencil_omega_theta
end interface weno5js_stencil_omega

interface weno5z_stencil_omega
	module procedure weno5z_stencil_omega_only
	module procedure weno5z_stencil_omega_theta
end interface weno5z_stencil_omega
contains

	subroutine weno5z_stencil(lf,h)
	implicit none
	double precision,intent(in) :: lf(-2:2)
	double precision,intent(out) :: h
	double precision :: q30,q31,q32
	double precision :: w0,w1,w2

	q30=a300*lf(-2)+a301*lf(-1)+a302*lf(0)
	q31=a310*lf(-1)+a311*lf(0) +a312*lf(1)
	q32=a320*lf(0) +a321*lf(1) +a322*lf(2)

	call weno5z_stencil_omega(lf,w0,w1,w2)

	h=w0*q30+w1*q31+w2*q32
	end subroutine


	subroutine weno5z_stencil_omega_only(lf,w0,w1,w2)
	implicit none
	double precision,intent(in) :: lf(-2:2)
	double precision,intent(out) :: w0,w1,w2
	double precision :: aa0,aa1,aa2,aa_sum
	
	call weno5z_stencil_alpha(lf,aa0,aa1,aa2)
	aa_sum = aa0+aa1+aa2

	w0=aa0/aa_sum
	w1=aa1/aa_sum
	w2=aa2/aa_sum 
	end subroutine


	subroutine weno5z_stencil_omega_theta(lf,w0,w1,w2,theta)
	implicit none
	double precision,intent(in) :: lf(-2:2)
	double precision,intent(out) :: w0,w1,w2,theta
	double precision :: aa0,aa1,aa2,aa_sum
	
	call weno5z_stencil_alpha(lf,aa0,aa1,aa2)
	aa_sum = aa0+aa1+aa2

	w0=aa0/aa_sum
	w1=aa1/aa_sum
	w2=aa2/aa_sum 
	theta = 1.0/(1.0+(aa_sum-1.0)**2)
	end subroutine


	subroutine weno5js_stencil(lf,h)
	implicit none
	double precision,intent(in) :: lf(-2:2)
	double precision,intent(out) :: h
	double precision :: q30,q31,q32
	double precision :: w0,w1,w2

	q30=a300*lf(-2)+a301*lf(-1)+a302*lf(0)
	q31=            a310*lf(-1)+a311*lf(0) +a312*lf(1)
	q32=                        a320*lf(0) +a321*lf(1) +a322*lf(2)

	call weno5js_stencil_omega(lf,w0,w1,w2)

	h=w0*q30+w1*q31+w2*q32
	end subroutine


	subroutine weno5js_stencil_omega_only(lf,w0,w1,w2)
	implicit none
	double precision,intent(in) :: lf(-2:2)
	double precision,intent(out) :: w0,w1,w2
	double precision :: aa0,aa1,aa2,aa_sum
	
	call weno5js_stencil_alpha(lf,aa0,aa1,aa2)
	aa_sum = aa0+aa1+aa2

	w0=aa0/aa_sum
	w1=aa1/aa_sum
	w2=aa2/aa_sum 
	end subroutine


	subroutine weno5js_stencil_omega_theta(lf,w0,w1,w2,theta)
	implicit none
	double precision,intent(in) :: lf(-2:2)
	double precision,intent(out) :: w0,w1,w2,theta
	double precision :: aa0,aa1,aa2,aa_sum
	
	call weno5js_stencil_alpha(lf,aa0,aa1,aa2)
	aa_sum = aa0+aa1+aa2

	w0=aa0/aa_sum
	w1=aa1/aa_sum
	w2=aa2/aa_sum 
	theta = 1.0/(1.0+(aa_sum-1.0)**2)
	end subroutine


	subroutine weno5js_stencil_alpha(lf,aa0,aa1,aa2)
	implicit none
	double precision,intent(in) :: lf(-2:2)
	double precision,intent(out) :: aa0,aa1,aa2
	double precision :: is0,is1,is2

	is0= 13.d0*(lf(-2)-2.d0*lf(-1)+lf(0))**2/12.d0 &
	    +(lf(-2)-4.d0*lf(-1)+3.d0*lf(0))**2/4.d0

	is1= 13.d0*(lf(-1)-2.d0*lf(0)+lf(1))**2/12.d0 &
		+(lf(-1)-lf(1))**2/4.d0

	is2= 13.d0*(lf(0)-2.d0*lf(1)+lf(2))**2/12.d0 &
		+(3.d0*lf(0)-4.d0*lf(1)+lf(2))**2/4.d0


	aa0=c30/(ss+is0)**2
	aa1=c31/(ss+is1)**2
	aa2=c32/(ss+is2)**2
	end subroutine


	subroutine weno5z_stencil_alpha(lf,aa0,aa1,aa2)
	implicit none
	double precision,intent(in) :: lf(-2:2)
	double precision,intent(out) :: aa0,aa1,aa2
	double precision :: is0,is1,is2
	double precision :: tau5

	is0= 13.d0*(lf(-2)-2.d0*lf(-1)+lf(0))**2/12.d0 &
	    +(lf(-2)-4.d0*lf(-1)+3.d0*lf(0))**2/4.d0

	is1= 13.d0*(lf(-1)-2.d0*lf(0)+lf(1))**2/12.d0 &
		+(lf(-1)-lf(1))**2/4.d0

	is2= 13.d0*(lf(0)-2.d0*lf(1)+lf(2))**2/12.d0 &
		+(3.d0*lf(0)-4.d0*lf(1)+lf(2))**2/4.d0

	tau5 = abs(is2-is0)

    aa0=c30*(1.0+(tau5/(ss+is0))**2)
    aa1=c31*(1.0+(tau5/(ss+is1))**2)
    aa2=c32*(1.0+(tau5/(ss+is2))**2)
	end subroutine
end module