program two_dimension_parallel
	use flow_general
	use omp_lib
	implicit none
	integer :: irk,i,j,nthread
	double precision :: cpu_t,cpu_t0,cpu_t1

	call read_par_namelist

    if(start_step.eq.0) then
        print*,'Start with initial conditon...'
	    call init_general
    else
        print*,'Start with data from file...'
        call load_date_binary
    endif

	call omp_set_num_threads(omp_thread_num)
	write(*,*) 'Number of Threads=',omp_thread_num
	open(20,file="CpuTime.dat",ACCESS='append')

	cpu_t = 0.0
	do while((current_time.lt.total_time).and. &
			 (current_step.lt.total_step))
		cpu_t0=OMP_get_wtime()
		do irk=1,stages
			call get_primary
			call boundary_general  
			call flux_splitting_general
			call time_rk3_general(irk)
		enddo
		cpu_t1=OMP_get_wtime()
		cpu_t=cpu_t+cpu_t1-cpu_t0
		write(*,*)  '=========','CPU time:',cpu_t1-cpu_t0,'Seconds','=========='
		current_step=current_step + 1
		current_time=current_time + dt
		elasp_time = elasp_time + dt
		write(*,*) current_step,'Flow time= ',current_time,'dt=', dt 

		if((MOD(current_step,save_step).eq.0).or.(elasp_time.ge.dt_save)) then
			print*,"Saving binary file at time step = ",current_step
            
            if(if_save_binary) then
            	call write_date_binary
            endif
            if(if_save_tec) then
				print*,"Saving TEC data file at time = ",current_time
				call get_primary
				call write_tec_at_time
				elasp_time = 0.0
			endif
			print*,"Saving finished!"
			call read_par_namelist
		endif
	enddo
	write(20,*) 'Total time and steps:',cpu_t,current_step
	write(20,*) 'averaged time:',cpu_t/current_step
	close(20)
	print*,'Writing final data... '
	call get_primary
	call write_tec_at_time
	print*,'Calculation is completed! '
end
