function generate_proposal_positions(	X::Array{Float64, 1}, 
									Y::Array{Float64, 1}, 
									Z::Array{Float64, 1},  
									Lx::Float64, 
									Ly::Float64, 
									Lz::Float64, 
									sigma_translation::Float64)
	
	number_of_particles::Int64 = length(X)
	
	X_star::Array{Float64, 1} = X + sigma_translation * randn(number_of_particles)
	Y_star::Array{Float64, 1} = Y + sigma_translation * randn(number_of_particles)
	Z_star::Array{Float64, 1} = Z + sigma_translation * randn(number_of_particles)
	
	for current_particle = 1:number_of_particles
		if X_star[current_particle] < 0.0
			X_star[current_particle] += Lx
		elseif X_star[current_particle] > Lx
			X_star[current_particle] -= Lx
		end
		
		if Y_star[current_particle] < 0.0
			Y_star[current_particle] += Ly
		elseif Y_star[current_particle] > Ly
			Y_star[current_particle] -= Ly
		end
		
		if Z_star[current_particle] < 0.0
			Z_star[current_particle] += Lz
		elseif Z_star[current_particle] > Lz
			Z_star[current_particle] -= Lz
		end
	end
	
	return (X_star, Y_star, Z_star)
end
		
		
	