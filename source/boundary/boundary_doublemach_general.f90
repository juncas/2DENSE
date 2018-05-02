subroutine boundary_doublemach_general
use flow_general
use var_index
implicit none
integer :: index_i,index_j,index_var,index_eq
integer :: bfsize_i,bfsize_j,ni,nj
double precision :: rho,u,v,p
double precision :: x,y

ni = block_grid(num_block)%ni
nj = block_grid(num_block)%nj
bfsize_i=block_grid(num_block)%bfsize_i
bfsize_j=block_grid(num_block)%bfsize_j

! Bottom reflective boundary
do index_j=1,bfsize_j
    do index_i=bfsize_i+1,ni+bfsize_i
        x = block_grid(num_block)%x(index_i-bfsize_i,1)
        if(x.ge.1.d0/6.d0) then
              block_grid(num_block)%cons(index_i,index_j,1) &
            = block_grid(num_block)%cons(index_i,2*bfsize_j-index_j+1,1) 
              
              block_grid(num_block)%cons(index_i,index_j,2) &
            = block_grid(num_block)%cons(index_i,2*bfsize_j-index_j+1,2) 
              
              block_grid(num_block)%cons(index_i,index_j,3) &
            =-block_grid(num_block)%cons(index_i,2*bfsize_j-index_j+1,3) 
              
              block_grid(num_block)%cons(index_i,index_j,4) &
            = block_grid(num_block)%cons(index_i,2*bfsize_j-index_j+1,4) 
              
              block_grid(num_block)%prim(index_i,index_j,1) &
            = block_grid(num_block)%prim(index_i,2*bfsize_j-index_j+1,1) 
              
              block_grid(num_block)%prim(index_i,index_j,2) &
            = block_grid(num_block)%prim(index_i,2*bfsize_j-index_j+1,2) 
              
              block_grid(num_block)%prim(index_i,index_j,3) &
            =-block_grid(num_block)%prim(index_i,2*bfsize_j-index_j+1,3) 
              
              block_grid(num_block)%prim(index_i,index_j,4) &
            = block_grid(num_block)%prim(index_i,2*bfsize_j-index_j+1,4) 

              block_grid(num_block)%prim(index_i,index_j,5) &
            = block_grid(num_block)%prim(index_i,2*bfsize_j-index_j+1,5) 
        else
            rho=8.d0
            u=8.25*sqrt(3.d0)*0.5
            v=-8.25*0.5d0
            p=116.5d0
              
            block_grid(num_block)%cons(index_i,index_j,1)=rho
            block_grid(num_block)%cons(index_i,index_j,2)=rho*u
            block_grid(num_block)%cons(index_i,index_j,3)=rho*v
            block_grid(num_block)%cons(index_i,index_j,4)=p/(gamma-1.d0)+0.5*rho*(u**2+v**2)    
            
            block_grid(num_block)%prim(index_i,index_j,1)=rho
            block_grid(num_block)%prim(index_i,index_j,2)=u
            block_grid(num_block)%prim(index_i,index_j,3)=v
            block_grid(num_block)%prim(index_i,index_j,4)=p
            block_grid(num_block)%prim(index_i,index_j,5)=sqrt(gamma*p/rho)
        endif
    enddo
enddo

! Top boundary
do index_j=nj+bfsize_j+1,nj+2*bfsize_j
    do index_i=bfsize_i+1,ni+bfsize_i
        x = block_grid(num_block)%x(index_i-bfsize_i,nj)

        if(x.lt.(1.d0/6.d0+(1.d0+20.d0*current_time)/sqrt(3.d0))) then
            rho=8.d0
            u=8.25*sqrt(3.d0)*0.5
            v=-8.25*0.5d0
            p=116.5d0
        else
            rho=1.4d0
            u=0.d0
            v=0.d0
            p=1.d0
        endif       
    
        block_grid(num_block)%cons(index_i,index_j,1)=rho
        block_grid(num_block)%cons(index_i,index_j,2)=rho*u
        block_grid(num_block)%cons(index_i,index_j,3)=rho*v
        block_grid(num_block)%cons(index_i,index_j,4)=p/(gamma-1.d0)+0.5*rho*(u**2+v**2)    
        
        block_grid(num_block)%prim(index_i,index_j,1)=rho
        block_grid(num_block)%prim(index_i,index_j,2)=u
        block_grid(num_block)%prim(index_i,index_j,3)=v
        block_grid(num_block)%prim(index_i,index_j,4)=p
        block_grid(num_block)%prim(index_i,index_j,5)=sqrt(gamma*p/rho)
    
    enddo
enddo
              

! left inlet boundary
do index_j=bfsize_j+1,nj+bfsize_j
    do index_i=1,bfsize_i
        rho=8.d0
        u=8.25*sqrt(3.d0)*0.5
        v=-8.25*0.5d0
        p=116.5d0
          
        block_grid(num_block)%cons(index_i,index_j,1)=rho
        block_grid(num_block)%cons(index_i,index_j,2)=rho*u
        block_grid(num_block)%cons(index_i,index_j,3)=rho*v
        block_grid(num_block)%cons(index_i,index_j,4)=p/(gamma-1.d0)+0.5*rho*(u**2+v**2)    
        
        block_grid(num_block)%prim(index_i,index_j,1)=rho
        block_grid(num_block)%prim(index_i,index_j,2)=u
        block_grid(num_block)%prim(index_i,index_j,3)=v
        block_grid(num_block)%prim(index_i,index_j,4)=p
        block_grid(num_block)%prim(index_i,index_j,5)=sqrt(gamma*p/rho)
    
    enddo
enddo

! Right outlet boundary
do index_eq = 1,neq
    do index_j=bfsize_j+1,nj+bfsize_j
        do index_i=ni+bfsize_i+1,ni+2*bfsize_i
              block_grid(num_block)%cons(index_i,index_j,index_eq) &
            = block_grid(num_block)%cons(index_i-1,index_j,index_eq)*4.0/3.0 &
             -block_grid(num_block)%cons(index_i-2,index_j,index_eq)/3.0
        enddo
    enddo
enddo

do index_eq = 1,prim_lastvar_index
    do index_j=bfsize_j+1,nj+bfsize_j
        do index_i=ni+bfsize_i+1,ni+2*bfsize_i
              block_grid(num_block)%prim(index_i,index_j,index_eq) &
            = block_grid(num_block)%prim(index_i-1,index_j,index_eq)*4.0/3.0 &
             -block_grid(num_block)%prim(index_i-2,index_j,index_eq)/3.0
        enddo
    enddo
enddo
              

end subroutine
