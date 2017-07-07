function signed_distance_mod(x1::Float64, x2::Float64, L::Float64)
	x12::Float64 = x2 - x1
	if x12 < -0.5 * L
		x12 += L
	elseif x12 > 0.5 * L
		x12 -= L
	end
	
	return x12
end