function overlap_cuboid(	xAB::Float64,
						yAB::Float64,
						zAB::Float64,
						a11::Float64,
						a12::Float64,
						a13::Float64,
						a21::Float64,
						a22::Float64,
						a23::Float64,
						a31::Float64,
						a32::Float64,
						a33::Float64,
						b11::Float64,
						b12::Float64,
						b13::Float64,
						b21::Float64,
						b22::Float64,
						b23::Float64,
						b31::Float64,
						b32::Float64,
						b33::Float64,
						r1a::Float64,
						r2a::Float64,
						r3a::Float64,
						r1b::Float64,
						r2b::Float64,
						r3b::Float64)

	overlapfun::Float64 = 0.0
	overlapfun_binary::Float64 = 0.0

	overlapfun_binary = overlap_cuboid_binary(xAB, yAB, zAB, a11, a12, a13, a21, a22, a23, a31, a32, a33, b11, b12, b13, b21, b22, b23, b31, b32, b33, r1a, r2a, r3a, r1b, r2b, r3b)

	return overlapfun_binary

	if overlapfun_binary == 0.0
		return 0.0
	else
		gamma_min::Float64 = 0.0
		gamma_max::Float64 = 1.0
		gamma::Float64 = 0.5
		number_of_iterations::Int64 = 10

		for current_iteration = 1:number_of_iterations
			if overlap_cuboid_binary(xAB, yAB, zAB, a11, a12, a13, a21, a22, a23, a31, a32, a33, b11, b12, b13, b21, b22, b23, b31, b32, b33, gamma * r1a, gamma * r2a, gamma * r3a, gamma * r1b, gamma * r2b, gamma * r3b) == 1.0
				gamma_max = gamma
			else
				gamma_min = gamma
			end
			gamma = 0.5 * (gamma_min + gamma_max)
		end

		overlapfun = 1.0 - gamma
	end
end

# xAB = 0.01
# yAB = xAB
# zAB = xAB
# a11 = 1.0
# a12 = 0.0
# a13 = 0.0
# a21 = 0.0
# a22 = 1.0
# a23 = 0.0
# a31 = 0.0
# a32 = 0.0
# a33 = 1.0
# b11 = 1.0
# b12 = 0.0
# b13 = 0.0
# b21 = 0.0
# b22 = 1.0
# b23 = 0.0
# b31 = 0.0
# b32 = 0.0
# b33 = 1.0
# r1a = 1.0
# r2a = 0.75
# r3a = 1.0
# r1b = 1.0
# r2b = 1.5
# r3b = 2.0

# overlap_cuboid(	xAB,
				# yAB,
				# zAB,
				# a11,
				# a12,
				# a13,
				# a21,
				# a22,
				# a23,
				# a31,
				# a32,
				# a33,
				# b11,
				# b12,
				# b13,
				# b21,
				# b22,
				# b23,
				# b31,
				# b32,
				# b33,
				# r1a,
				# r2a,
				# r3a,
				# r1b,
				# r2b,
				# r3b)
