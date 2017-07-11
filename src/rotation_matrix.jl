function rotation_matrix(q0::Float64, q1::Float64, q2::Float64, q3::Float64)

	a11::Float64 = q0^2 + q1^2 - q2^2 - q3^2
	a12::Float64 = 2 * q1 * q2 + 2 * q0 * q3
	a13::Float64 = 2 * q1 * q3 - 2 * q0 * q2
	a21::Float64 = 2 * q1 * q2 - 2 * q0 * q3
	a22::Float64 = q0^2 - q1^2 + q2^2 - q3^2
	a23::Float64 = 2 * q2 * q3 + 2 * q0 * q1
	a31::Float64 = 2 * q1 * q3 + 2 * q0 * q2
	a32::Float64 = 2 * q2 * q3 - 2 * q0 * q1
	a33::Float64 = q0^2 - q1^2 - q2^2 + q3^2
	
	return (a11, a12, a13, a21, a22, a23, a31, a32, a33)

end