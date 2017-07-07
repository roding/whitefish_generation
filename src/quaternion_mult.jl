function quaternion_mult(	p0::Float64, 
						p1::Float64, 
						p2::Float64, 
						p3::Float64,
						q0::Float64, 
						q1::Float64, 
						q2::Float64, 
						q3::Float64)
						
	r0::Float64 = p0 * q0 - p1 * q1 - p2 * q2 - p3 * q3
	r1::Float64 = p0 * q1 + p1 * q0 - p2 * q3 + p3 * q2
	r2::Float64 = p0 * q2 + p2 * q0 + p1 * q3 - p3 * q1
	r3::Float64 = p0 * q3 - p1 * q2 + p2 * q1 + p3 * q0
	
	return (r0, r1, r2, r3)
end
