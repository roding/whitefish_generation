function generate_random_unit_quaternion()
	
	q0::Float64 = randn()
	q1::Float64 = randn()
	q2::Float64 = randn()
	q3::Float64 = randn()
	
	q_norm::Float64 = sqrt( q0^2 + q1^2 + q2^2 + q3^2 )

	q0 /= q_norm
	q1 /= q_norm
	q2 /= q_norm
	q3 /= q_norm
	
	return (q0, q1, q2, q3)
end