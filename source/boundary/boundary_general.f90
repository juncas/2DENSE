subroutine boundary_general
use flow_general
implicit none

select case(benchmark_index)
case(1)
    call boundary_ivc_general
case(2)
    call boundary_sedov_general
case(3)
    call boundary_rt_general
case(4)
    call boundary_rm_general
case(5)
    call boundary_doublemach_general
case(6)
    call boundary_ssl_general
case(7)
    !call boundary_forwardstep_general
case(8)
    !call boundary_backwardstep_general
case(9)
    call boundary_shockvortex_general
case(0)
    !call boundary_user_general
	case default 
        print*,"No such Benchmark index!"
        stop
end select

call set_buffer

end subroutine

subroutine set_buffer
use flow_general
implicit none


end subroutine set_buffer