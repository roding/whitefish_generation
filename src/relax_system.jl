function relax_system(		Lx::Float64, 
						Ly::Float64, 
						Lz::Float64, 
						R1::Array{Float64, 1}, 
						R2::Array{Float64, 1}, 
						R3::Array{Float64, 1}, 
						X::Array{Float64, 1}, 
						Y::Array{Float64, 1}, 
						Z::Array{Float64, 1}, 
						Q0::Array{Float64, 1}, 
						Q1::Array{Float64, 1}, 
						Q2::Array{Float64, 1}, 
						Q3::Array{Float64, 1}, 
						A11::Array{Float64, 1},
						A12::Array{Float64, 1},
						A13::Array{Float64, 1},
						A21::Array{Float64, 1},
						A22::Array{Float64, 1},
						A23::Array{Float64, 1},
						A31::Array{Float64, 1},
						A32::Array{Float64, 1},
						A33::Array{Float64, 1},
						sigma_translation::Float64,
						sigma_rotation::Float64)

	number_of_particles::Int64 = length(R1)

	# Pre-computed max radii.
	RMAX::Array{Float64, 1} = zeros(number_of_particles)
	for current_particle = 1:number_of_particles
		RMAX[current_particle] = maximum( (R1[current_particle], R2[current_particle], R3[current_particle]) )
	end
	
	# Preallocation.
	#X_star::Array{Float64, 1} = zeros(number_of_particles)
	#Y_star::Array{Float64, 1} = zeros(number_of_particles)
	#Z_star::Array{Float64, 1} = zeros(number_of_particles)
	Q0_star::Array{Float64, 1} = zeros(number_of_particles)
	Q1_star::Array{Float64, 1} = zeros(number_of_particles)
	Q2_star::Array{Float64, 1} = zeros(number_of_particles)
	Q3_star::Array{Float64, 1} = zeros(number_of_particles)
	A11_star::Array{Float64, 1} = zeros(number_of_particles)
	A12_star::Array{Float64, 1} = zeros(number_of_particles)
	A13_star::Array{Float64, 1} = zeros(number_of_particles)
	A21_star::Array{Float64, 1} = zeros(number_of_particles)
	A22_star::Array{Float64, 1} = zeros(number_of_particles)
	A23_star::Array{Float64, 1} = zeros(number_of_particles)
	A31_star::Array{Float64, 1} = zeros(number_of_particles)
	A32_star::Array{Float64, 1} = zeros(number_of_particles)
	A33_star::Array{Float64, 1} = zeros(number_of_particles)
	
	x_star::Float64 = 0.0
	y_star::Float64 = 0.0
	z_star::Float64 = 0.0
	
	energy_system::Float64 = 0.0
	energy_particle::Float64 = 0.0
	energy_particle_star::Float64 = 0.0
	current_sweep::Int64 = 0
	acceptance_probability_translation::Float64 = 0.0
	acceptance_probability_rotation::Float64 = 0.0
	xAB::Float64 = 0.0
	yAB::Float64 = 0.0
	zAB::Float64 = 0.0
	overlapfun::Float64 = 0.0
	
	energy_system = 1.0
	while energy_system > 0.0
		current_sweep += 1
		println(join(["   Sweep ", string(current_sweep)]))
		
		acceptance_probability_translation = 0.0
		acceptance_probability_rotation = 0.0
		
		energy_system = 0.0
		
		for currentA = 1:number_of_particles
			# Compute current local energy.
			energy_particle = 0.0
			for currentB = [1:currentA-1 ; currentA+1:number_of_particles]
				xAB = signed_distance_mod(X[currentA], X[currentB], Lx)
				yAB = signed_distance_mod(Y[currentA], Y[currentB], Ly)
				zAB = signed_distance_mod(Z[currentA], Z[currentB], Lz)

				if xAB^2 + yAB^2 + zAB^2 < (RMAX[currentA] + RMAX[currentB])^2
					overlapfun = overlap_function(xAB, yAB, zAB, A11[currentA], A12[currentA], A13[currentA], A21[currentA], A22[currentA], A23[currentA], A31[currentA], A32[currentA], A33[currentA], A11[currentB], A12[currentB], A13[currentB], A21[currentB], A22[currentB], A23[currentB], A31[currentB], A32[currentB], A33[currentB], R1[currentA]^2 * R2[currentA]^2 * R3[currentA]^2)
					
					if overlapfun < 1.0
						energy_particle += (1.0-overlapfun)^2
					end
				end
			end

			# Generate random proposal position and compute new local energy with translation.
			(x_star, y_star, z_star) = generate_proposal_position(X[currentA], Y[currentA], Z[currentA], Lx, Ly, Lz, sigma_translation)
			energy_particle_star = 0.0
			for currentB = [1:currentA-1;currentA+1:number_of_particles]
				xAB = signed_distance_mod(x_star, X[currentB], Lx)
				yAB = signed_distance_mod(y_star, Y[currentB], Ly)
				zAB = signed_distance_mod(z_star, Z[currentB], Lz)

				if xAB^2 + yAB^2 + zAB^2 < (RMAX[currentA] + RMAX[currentB])^2
					overlapfun = overlap_function(xAB, yAB, zAB, A11[currentA], A12[currentA], A13[currentA], A21[currentA], A22[currentA], A23[currentA], A31[currentA], A32[currentA], A33[currentA], A11[currentB], A12[currentB], A13[currentB], A21[currentB], A22[currentB], A23[currentB], A31[currentB], A32[currentB], A33[currentB], R1[currentA]^2 * R2[currentA]^2 * R3[currentA]^2)
					
					if overlapfun < 1.0
						energy_particle_star += (1.0-overlapfun)^2
					end
				end
			end
			
			if energy_particle_star <= energy_particle
				X[currentA] = x_star
				Y[currentA] = y_star
				Z[currentA] = z_star
				acceptance_probability_translation += 1.0
				energy_particle = energy_particle_star
			end
			
			# Generate random proposal orientation and compute new local energy with rotation.
			(x_star, y_star, z_star) = generate_proposal_orientation(X[currentA], Y[currentA], Z[currentA], Lx, Ly, Lz, sigma_translation)
			
			
			
			
			
			energy_system += energy_particle
		
		end
		
		
		
		
				
				
				
		acceptance_probability_translation /= number_of_particles
		acceptance_probability_rotation /= number_of_particles
				
		################# REDO WITH MAX LIMITS!!!!
		if acceptance_probability_translation <= acceptance_probability_target
			sigma_translation *= 0.95
		else
			sigma_translation *= 1.05
			sigma_translation = min(sigma_translation, 10.0)
			#println(sigma_translation)
		end
		
		if acceptance_probability_rotation <= acceptance_probability_target
			sigma_rotation *= 0.95
		else
			sigma_rotation *= 1.05
			sigma_rotation = min(sigma_rotation, 1.0)
		end	
	end
	
	return (X, Y, Z, Q0, Q1, Q2, Q3, A11, A12, A13, A21, A22, A23, A31, A32, A33, sigma_translation, sigma_rotation)
end
