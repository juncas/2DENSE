subroutine init_general
use flow_general
!   1 - isentropic vortex convection
!   2 - Sedov problem
!   3 - Rayleigh-Taylor instability
!   4 - Richtmyer-Meshkov instability
!   5 - Double Mach refection
!   6 - Shock-shear layer interaction
!   7 - Forward step
!   8 - Backward Step
!   9 - Shock-vortex interaction
!   0 - user defined
select case(benchmark_index)
case(1)
    call init_ivc_general
case(2)
    call init_sedov_general
case(3)
    call init_rt_general
case(4)
    call init_rm_general
case(5)
    call init_doublemach_general
case(6)
    call init_ssl_general
case(7)
    !call init_forwardstep_general
case(8)
    !call init_backwardstep_general
case(9)
    call init_shockvortex_general
case(0)
    !call init_user_general
case default 
        print*,"No such Benchmark index!"
        stop
end select
end subroutine