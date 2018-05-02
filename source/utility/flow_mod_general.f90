module data_types
    TYPE :: block
        integer :: index
        integer :: ni,nj
        integer :: bfsize_i,bfsize_j
        double precision :: x_corner(2,2)
        double precision :: y_corner(2,2)
        double precision :: xl,yl,dx,dy
        double precision :: max_eig_i,max_eig_j
        
        ! Geometry data 
        double precision,pointer :: x(:,:),y(:,:)
        double precision,pointer :: vet_x(:,:),vet_y(:,:)

        double precision,pointer :: xi_x(:,:),eta_x(:,:)
        double precision,pointer :: xi_y(:,:),eta_y(:,:)
        double precision,pointer :: jacob(:,:)

        integer :: block_topo(4) ! first index: lower,right,upper,left

        ! Flow field data 
        double precision,pointer :: cons(:,:,:)
        double precision,pointer :: prim(:,:,:)

        double precision,pointer :: cons_temp(:,:,:)

        double precision,pointer :: flux_xi(:,:,:),flux_eta(:,:,:)
        double precision,pointer :: flux_xip(:,:,:),flux_xin(:,:,:), &
                                    flux_etap(:,:,:),flux_etan(:,:,:)

        double precision,pointer :: h_xip(:,:,:),h_xin(:,:,:), &
                                    h_etan(:,:,:),h_etap(:,:,:)
        double precision,pointer :: h_xi(:,:,:),h_eta(:,:,:)
        double precision,pointer :: theta_xip(:,:,:),theta_etap(:,:,:)
        double precision,pointer :: theta_xin(:,:,:),theta_etan(:,:,:)

        double precision,pointer :: vis_xi(:,:,:),vis_eta(:,:,:)

        double precision,pointer :: bodyforce(:,:,:)
        double precision,pointer :: source(:,:,:)

        double precision,pointer :: rhs(:,:,:)

    end TYPE block
end module data_types

module flow_general
use data_types
implicit none
    integer :: omp_thread_num = 4
    NAMELIST/omp_pars/ omp_thread_num

    integer :: neq = 4 ! number of equations
    double precision :: mach_inf = 1.d0,reynolds = 500.d0
    logical :: if_source =.false.
    logical :: if_bodyforce =.false.
    logical :: if_viscous = .false.
    integer :: EOS_index = 1
    integer :: num_material = 1
    integer :: benchmark_index = 0
    ! Benchmark:
    !   0 - User defined initial condition
    !   1 - Isentropic vortex convection
    !   2 - Sedov problem
    !   3 - Rayleigh-Taylor instability
    !   4 - Richtmyer-Meshkov instability
    !   5 - Double mach reflection
    !   6 - Shock/shear-layer interaction
    !   7 - Forward step
    !   8 - Backward step

    NAMELIST/model_pars/ neq,if_source,if_viscous,if_bodyforce,mach_inf,reynolds, &
                         benchmark_index

    type(block),dimension(:),allocatable :: block_grid
    integer :: num_block = 1
    integer :: single_ni=64,single_nj=64
    logical :: if_curve_grid=.false.
    logical :: if_read_from_file=.false.
    NAMELIST/grid_pars/ if_curve_grid,if_read_from_file,single_ni,single_nj

    integer :: stages = 3
    double precision :: total_time=0.0,current_time=0.0, &
                        elasp_time = 0.0
    double precision :: dt = 0.01,cfl = 0.5
    logical :: if_const_dt = .false.
    logical :: if_save_tec = .true.
    logical :: if_save_binary = .true.


    integer :: total_step = 1, start_step = 0, &
               current_step=0,save_step=100 
    double precision :: dt_save=0.01
    NAMELIST/time_step_pars/ stages,total_time,dt_save,dt,cfl, &
                             total_step,save_step,if_const_dt, &
                             if_save_binary,if_save_tec

    ! Flux splitting methods:
    !   #Vector splitting methods:
    !       101 - local LF
    !       102 - global LF
    !       103 - Steger-warmming
    !       104 - AUSM-FVS
    !       105 - SemiCh
    !
    !   #Approximate Riemann solver methods:
    !       111 - ROE
    !       112 - AUSM
    !       113 - HLLC
    !
    integer :: scheme_splitting_index = 101

    ! Reconstruction Methods:
    !   #First order schemes   :
    !       210 - first order upwind
    !
    !   #Second order schemes  :
    !       220 - second order center 
    !       221 - MUSCL
    !
    !   #Third order schemes   :  
    !       230 - third order compact
    !
    !       231 - third order upwind
    !       233 - third order WENO-JS
    !
    !   #Fourth order schemes  :
    !       240 - Fourth order center
    !
    !       242 - Fourth order compact
    !       
    !   #Fifth order schemes   :
    !       250 - fifth order compact
    !       252 - fifth order CRWENO
    !
    !       251 - fifth order upwind
    !       253 - fifth order WENO-JS
    !       255 - fifth order WENO-Z
    !       257 - fifth order WENO-SC
    !   
    !   #Sixth order schemes   :
    !       260 - sixth order compact
    !
    !       261 - sixth order center
    !
    !   #Seventh order schemes :
    !       271 - seventh order compact
    !       273 - seventh order HCCS
    !
    !       270 - seventh order upwind
    !       272 - seventh order WENO-JS
    !       274 - seventh order WENO-Z
    !
    integer :: scheme_inv_index = 201
    logical :: if_save_theta= .false.
    double precision :: delta = 0.5d0

    ! Methods for second order terms:
    !   #Second order schemes  :
    !       320 - second order center 
    !
    !   #Fourth order schemes  :
    !       340 - Fourth order center 
    !
    !   #Sixth order schemes  :
    !       360 - Sixth order center 
    !
    !   #Eigth order schemes  :
    !       380 - Eigth order center 
    !
    integer :: scheme_vis_index = 301
    NAMELIST/scheme_pars/ scheme_splitting_index, &
                          scheme_inv_index, &
                          scheme_vis_index, &
                          if_save_theta, &
                          delta

    double precision :: gamma=1.4d0
    NAMELIST/gas_pars/ gamma

    double precision,parameter :: pie = atan(1.d0)*4.d0
    double precision :: max_eig_i,max_eig_j

contains
    subroutine allocate_block_grid
    implicit none
    allocate(block_grid(num_block))
    end subroutine allocate_block_grid

    subroutine allocate_block(myblock,ni,nj,bfsize_i,bfsize_j,neq)    
    use var_index
    implicit none
    integer,intent(in) :: ni,nj
    integer,intent(in) :: bfsize_i,bfsize_j
    integer,intent(in) :: neq

    type(block),intent(inout) :: myblock

    myblock%ni=ni
    myblock%nj=nj
    myblock%bfsize_i=bfsize_i
    myblock%bfsize_j=bfsize_j


        allocate(myblock%cons(ni+2*bfsize_i,nj+2*bfsize_j,cons_lastvar_index))
        allocate(myblock%prim(ni+2*bfsize_i,nj+2*bfsize_j,prim_lastvar_index))
        myblock%cons = 0.d0
        myblock%prim = 0.d0
        
        allocate(myblock%cons_temp(ni+2*bfsize_i,nj+2*bfsize_j,neq))
        myblock%cons_temp= 0.d0
        
        allocate(myblock%flux_xi(ni+2*bfsize_i,nj,neq), &
                 myblock%flux_eta(nj+2*bfsize_j,ni,neq))
        myblock%flux_xi= 0.d0
        myblock%flux_eta= 0.d0
        
        allocate(myblock%flux_xip(ni+2*bfsize_i,nj,neq), &
                 myblock%flux_xin(ni+2*bfsize_i,nj,neq), &
                 myblock%flux_etap(nj+2*bfsize_j,ni,neq), &
                 myblock%flux_etan(nj+2*bfsize_j,ni,neq))
        myblock%flux_xip = 0.d0
        myblock%flux_xin = 0.d0
        myblock%flux_etap = 0.d0
        myblock%flux_etan = 0.d0

        allocate(myblock%h_xip(ni+1,nj,neq), &
                 myblock%h_xin(ni+1,nj,neq), &
                 myblock%h_etap(nj+1,ni,neq),&
                 myblock%h_etan(nj+1,ni,neq))
        myblock%h_xip = 0.d0
        myblock%h_xin = 0.d0
        myblock%h_etan = 0.d0
        myblock%h_etap = 0.d0

        allocate(myblock%h_xi(ni+1,nj,neq), &
                 myblock%h_eta(nj+1,ni,neq))
        myblock%h_xi = 0.d0
        myblock%h_eta = 0.d0

        allocate(myblock%theta_xip(4,ni+1,nj), &
                 myblock%theta_etap(4,nj+1,ni))
        myblock%theta_xip = 0.d0
        myblock%theta_etap = 0.d0
        
        allocate(myblock%theta_xin(4,ni+1,nj), &
                 myblock%theta_etan(4,nj+1,ni))
        myblock%theta_xin = 0.d0
        myblock%theta_etan = 0.d0

        allocate(myblock%x(ni,nj), &
                 myblock%y(ni,nj))
        myblock%x = 0.d0
        myblock%y = 0.d0

        allocate(myblock%vet_x(ni+1,nj+1), &
                 myblock%vet_y(ni+1,nj+1))
        myblock%vet_x = 0.d0
        myblock%vet_y = 0.d0

        if(if_curve_grid) then
            allocate(myblock%xi_x(ni+2*bfsize_i,nj+2*bfsize_j), &
                     myblock%eta_x(ni+2*bfsize_i,nj+2*bfsize_j))
            allocate(myblock%xi_y(ni+2*bfsize_i,nj+2*bfsize_j), &
                     myblock%eta_y(ni+2*bfsize_i,nj+2*bfsize_j))
            allocate(myblock%jacob(ni+2*bfsize_i,nj+2*bfsize_j))

            myblock%xi_x = 0.d0
            myblock%eta_x =0.d0
            myblock%xi_y = 0.d0
            myblock%eta_y =0.d0
            myblock%jacob = 0.d0
        endif

        if(if_viscous) then
            allocate(myblock%vis_xi(ni+2*bfsize_i,nj+2*bfsize_j,neq), &
                     myblock%vis_eta(nj+2*bfsize_j,ni+2*bfsize_i,neq))
            myblock%vis_xi = 0.d0
            myblock%vis_eta = 0.d0
        endif

        if(if_bodyforce) then 
            allocate(myblock%bodyforce(ni+2*bfsize_i,nj+2*bfsize_j,neq))
            myblock%bodyforce = 0.d0
        endif
        
        if(if_source) then 
            allocate(myblock%source(ni+2*bfsize_i,nj+2*bfsize_j,neq))
            myblock%source = 0.d0
        endif


            allocate(myblock%rhs(ni,nj,neq))
            myblock%rhs = 0.0d0
    end subroutine allocate_block

    subroutine read_par_namelist
    implicit none
    logical :: if_exist = .false.
    inquire(FILE ='input.in',EXIST = if_exist )
    if(.not.if_exist) then
        print*,'Input file input.in does not exist!'
        print*,'Write your own input file by modifying input.in.example.'
        open(10,FILE='input.in.example')
        write(10,NML=omp_pars)
        write(10,NML=model_pars)
        write(10,NML=grid_pars)
        write(10,NML=time_step_pars)
        write(10,NML=scheme_pars)
        write(10,NML=gas_pars)    
        close(10)
    else
        open(10,FILE='input.in')
        read(10,NML=omp_pars)
        !print*,"omp_pars has been read..."
        !write(*,NML=omp_pars)
        read(10,NML=model_pars)
        !print*,"model_pars has been read..."
        !write(*,NML=model_pars)
        read(10,NML=grid_pars)
        !print*,"grid_pars has been read..."
        !write(*,NML=grid_pars)
        read(10,NML=time_step_pars)
        !print*,"time_step_pars has been read..."
        !write(*,NML=time_step_pars)
        read(10,NML=scheme_pars)
        !print*,"scheme_pars has been read..."
        !write(*,NML=scheme_pars)
        read(10,NML=gas_pars)    
        !print*,"gas_pars has been read..."
        !write(*,NML=gas_pars)
        close(10)
    endif
    end subroutine read_par_namelist

    subroutine write_par_namelist
    implicit none
    logical :: if_exist = .false.
    inquire(FILE ='input.in',EXIST = if_exist )
    if(.not.if_exist) then
        print*,'Input file input.in does not exist!'
        print*,'Write your own input file by modifying input.in.example.'
        open(10,FILE='input.in.example')
        write(10,NML=omp_pars)
        write(10,NML=model_pars)
        write(10,NML=grid_pars)
        write(10,NML=time_step_pars)
        write(10,NML=scheme_pars)
        write(10,NML=gas_pars)    
        close(10)
    endif
    end subroutine write_par_namelist

end module flow_general
