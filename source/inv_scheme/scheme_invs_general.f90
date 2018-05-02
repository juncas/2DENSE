subroutine scheme_invis_general
use flow_general
implicit none
select case (scheme_inv_index)
	case(251)
		call scheme_upwind5_general
	case(253)
		call scheme_weno5js_general
	case(255)
		call scheme_weno5z_general
	case(2)
		!call scheme_crweno_general
	case(3)
		!call scheme_mscw_general
	case(4)
		!call scheme_hccs57_general
	case(5)
		!call scheme_hccs77_general
	case(6)
		!call scheme_mscw_m_general
	case(7)
		!call scheme_weno7_general
	case(8)
		!call scheme_weno3_general
	case(9)
		!call scheme_hccs_switch_general
	end select
end subroutine scheme_invis_general
