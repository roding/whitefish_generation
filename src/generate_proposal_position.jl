function generate_proposal_position(	x::Float64, 
									y::Float64, 
									z::Float64,  
									Lx::Float64, 
									Ly::Float64, 
									Lz::Float64, 
									sigma_translation::Float64)
	
	x_star::Float64 = x + sigma_translation * randn()
	y_star::Float64 = y + sigma_translation * randn()
	z_star::Float64 = z + sigma_translation * randn()
	
	if x_star < 0.0
		x_star += Lx
	elseif x_star > Lx
		x_star -= Lx
	end
	
	if y_star < 0.0
		y_star += Ly
	elseif y_star > Ly
		y_star -= Ly
	end
	
	if z_star < 0.0
		z_star += Lz
	elseif z_star > Lz
		z_star -= Lz
	end
	
	return (x_star, y_star, z_star)
end
		
		
	