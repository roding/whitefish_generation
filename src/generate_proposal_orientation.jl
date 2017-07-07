function generate_proposal_orientation(	q0::Float64, 
									q1::Float64, 
									q2::Float64,  
									q3::Float64,
									sigma_rotation::Float64)
	
	ux_rot::Float64 = randn()
	uy_rot::Float64 = randn()
	uz_rot::Float64 = randn()
	u_norm::Float64 = sqrt( ux_rot^2 + uy_rot^2 + uz_rot^2 )
	ux_rot /= u_norm
	uy_rot /= u_norm
	uz_rot /= u_norm
	
	theta_rot::Float64 = sigma_rotation * randn()
	
	q0_rot::Float64 = cos(0.5 * theta_rot)
	q1_rot::Float64 = ux_rot * sin(0.5 * theta_rot)
	q2_rot::Float64 = uy_rot * sin(0.5 * theta_rot)
	q3_rot::Float64 = uz_rot * sin(0.5 * theta_rot)
	
	(q0_star::Float64, q1_star::Float64, q2_star::Float64, q3_star::Float64) = quaternion_mult(q0_rot, q1_rot, q2_rot, q3_rot, q0, q1, q2, q3)
	
	# Normalize proposal quaternion to avoid numerical problems.
	q_norm::Float64 = sqrt( q0_star^2 + q1_star^2 + q2_star^2 + q3_star^2 )
	q0_star /= q_norm
	q1_star /= q_norm
	q2_star /= q_norm
	q3_star /= q_norm
	
	return (q0_star, q1_star, q2_star, q3_star)
end
		
		
	