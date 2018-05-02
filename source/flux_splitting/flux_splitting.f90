subroutine flux_splitting_general
use flow_general    
implicit none
select case (scheme_splitting_index)
	case(101)
		  call flux_splitting_LLF   
		  call scheme_invis_general
	case(102)
		  call flux_splitting_GLF   
		  call scheme_invis_general
	case(103)
		  call flux_splitting_SW  
		  call scheme_invis_general
	case(104)
		  call flux_splitting_CH  
	case(105)
		  call flux_splitting_SC  
	case(106)
		! call flux_ROE   
		! call scheme_invis_general
	case(107)
		! call flux_ROE   
		! call scheme_invis_general
	case(108)
		! call flux_AUSM  
		! call scheme_invis_general 
	case(109)
		! call flux_HLLC      
		! call scheme_invis_general 
    case default 
        print*,"No such Selection"
        stop
end select 
end subroutine
