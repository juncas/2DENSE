module var_index
implicit none
! this module should be modified if new primative variables are added.

integer,parameter :: prim_rho_index=1
integer,parameter :: prim_u_index=2
integer,parameter :: prim_v_index=3
integer,parameter :: prim_p_index=4
integer,parameter :: prim_alpha_index=5
integer,parameter :: prim_entropy_index=6
integer,parameter :: prim_lastvar_index=6

integer,parameter :: cons_rho_index=1
integer,parameter :: cons_rhou_index=2
integer,parameter :: cons_rhov_index=3
integer,parameter :: cons_rhoe_index=4
integer,parameter :: cons_lastvar_index=4

end module