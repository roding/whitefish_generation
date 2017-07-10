function compress_system(	Lx::Float64, 
						Ly::Float64, 
						Lz::Float64,
						particle_type::String,
						R::Array{Float64, 2},
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
						sigma_translation_ub::Float64,
						sigma_rotation::Float64,
						sigma_rotation_ub::Float64)

	number_of_particles::Int64 = size(R, 1)
	
	acceptance_probability_target::Float64 = 0.25

	# Pre-computed max radii (i.e. radius of bounding sphere).
	RMAX::Array{Float64, 1} = zeros(number_of_particles)
	if particle_type == "sphere"
		for current_particle = 1:number_of_particles
			RMAX[current_particle] = R[current_particle, 1]
		end
	elseif particle_type == "ellipse"
		for current_particle = 1:number_of_particles
			RMAX[current_particle] = maximum( R[current_particle, 1:2] )
		end
	elseif particle_type == "ellipsoid"
		for current_particle = 1:number_of_particles
			RMAX[current_particle] = maximum( R[current_particle, 1:3] )
		end
	end
	
	# Preallocation.
	x_star::Float64 = 0.0
	y_star::Float64 = 0.0
	z_star::Float64 = 0.0
	q0_star::Float64 = 0.0
	q1_star::Float64 = 0.0
	q2_star::Float64 = 0.0
	q3_star::Float64 = 0.0
	a11_star::Float64 = 0.0
	a12_star::Float64 = 0.0
	a13_star::Float64 = 0.0
	a21_star::Float64 = 0.0
	a22_star::Float64 = 0.0
	a23_star::Float64 = 0.0
	a31_star::Float64 = 0.0
	a32_star::Float64 = 0.0
	a33_star::Float64 = 0.0
	
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
					if particle_type == "sphere"
						overlapfun = (RMAX[currentA] + RMAX[currentB])^2 - (xAB^2 + yAB^2 + zAB^2)
						energy_particle += overlapfun
					elseif particle_type == "ellipse"
						overlapfun = overlap_ellipse(xAB, yAB, zAB, A11[currentA], A12[currentA], A13[currentA], A21[currentA], A22[currentA], A23[currentA], A31[currentA], A32[currentA], A33[currentA], A11[currentB], A12[currentB], A13[currentB], A21[currentB], A22[currentB], A23[currentB], A31[currentB], A32[currentB], A33[currentB])
					
						if overlapfun < 1.0
							energy_particle += (1.0 - overlapfun)^2
						end
					elseif particle_type == "ellipsoid"
						overlapfun = overlap_ellipsoid(xAB, yAB, zAB, A11[currentA], A12[currentA], A13[currentA], A21[currentA], A22[currentA], A23[currentA], A31[currentA], A32[currentA], A33[currentA], A11[currentB], A12[currentB], A13[currentB], A21[currentB], A22[currentB], A23[currentB], A31[currentB], A32[currentB], A33[currentB], R[currentA, 1]^2 * R[currentA, 2]^2 * R[currentA, 3]^2)
					
						if overlapfun < 1.0
							energy_particle += (1.0 - overlapfun)^2
						end
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
					if particle_type == "sphere"
						overlapfun = (RMAX[currentA] + RMAX[currentB])^2 - (xAB^2 + yAB^2 + zAB^2)
						energy_particle_star += overlapfun
					elseif particle_type == "ellipse"
						overlapfun = overlap_ellipse(xAB, yAB, zAB, A11[currentA], A12[currentA], A13[currentA], A21[currentA], A22[currentA], A23[currentA], A31[currentA], A32[currentA], A33[currentA], A11[currentB], A12[currentB], A13[currentB], A21[currentB], A22[currentB], A23[currentB], A31[currentB], A32[currentB], A33[currentB])
					
						if overlapfun < 1.0
							energy_particle_star += (1.0 - overlapfun)^2
						end
					elseif particle_type == "ellipsoid"
						overlapfun = overlap_ellipsoid(xAB, yAB, zAB, A11[currentA], A12[currentA], A13[currentA], A21[currentA], A22[currentA], A23[currentA], A31[currentA], A32[currentA], A33[currentA], A11[currentB], A12[currentB], A13[currentB], A21[currentB], A22[currentB], A23[currentB], A31[currentB], A32[currentB], A33[currentB], R[currentA, 1]^2 * R[currentA, 2]^2 * R[currentA, 3]^2)
					
						if overlapfun < 1.0
							energy_particle_star += (1.0 - overlapfun)^2
						end
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
			if particle_type != "sphere"
				(q0_star, q1_star, q2_star, q3_star) = generate_proposal_orientation(Q0[currentA], Q1[currentA], Q2[currentA], Q3[currentA], sigma_rotation)
				
				if particle_type == "ellipse"
					(a11_star, a12_star, a13_star, a21_star, a22_star, a23_star, a31_star, a32_star, a33_star) = characteristic_matrix_ellipse(q0_star, q1_star, q2_star, q3_star, R[currentA, 1], R[currentA, 2])
				elseif particle_type == "ellipsoid"
					(a11_star, a12_star, a13_star, a21_star, a22_star, a23_star, a31_star, a32_star, a33_star) = characteristic_matrix_ellipsoid(q0_star, q1_star, q2_star, q3_star, R[currentA, 1], R[currentA, 2], R[currentA, 3])
				end
				
				energy_particle_star = 0.0
				for currentB = [1:currentA-1;currentA+1:number_of_particles]
					xAB = signed_distance_mod(X[currentA], X[currentB], Lx)
					yAB = signed_distance_mod(Y[currentA], Y[currentB], Ly)
					zAB = signed_distance_mod(Z[currentA], Z[currentB], Lz)

					if xAB^2 + yAB^2 + zAB^2 < (RMAX[currentA] + RMAX[currentB])^2
						if particle_type == "ellipse"
							overlapfun = overlap_ellipse(xAB, yAB, zAB, a11_star, a12_star, a13_star, a21_star, a22_star, a23_star, a31_star, a32_star, a33_star, A11[currentB], A12[currentB], A13[currentB], A21[currentB], A22[currentB], A23[currentB], A31[currentB], A32[currentB], A33[currentB])
						
							if overlapfun < 1.0
								energy_particle_star += (1.0 - overlapfun)^2
							end
						elseif particle_type == "ellipsoid"
							overlapfun = overlap_ellipsoid(xAB, yAB, zAB, a11_star, a12_star, a13_star, a21_star, a22_star, a23_star, a31_star, a32_star, a33_star, A11[currentB], A12[currentB], A13[currentB], A21[currentB], A22[currentB], A23[currentB], A31[currentB], A32[currentB], A33[currentB], R[currentA, 1]^2 * R[currentA, 2]^2 * R[currentA, 3]^2)
						
							if overlapfun < 1.0
								energy_particle_star += (1.0 - overlapfun)^2
							end
						end
						
					end
				end
				
				if energy_particle_star <= energy_particle
					Q0[currentA] = q0_star
					Q1[currentA] = q1_star
					Q2[currentA] = q2_star
					Q3[currentA] = q3_star
					
					A11[currentA] = a11_star
					A12[currentA] = a12_star
					A13[currentA] = a13_star
					A21[currentA] = a21_star
					A22[currentA] = a22_star
					A23[currentA] = a23_star
					A31[currentA] = a31_star
					A32[currentA] = a32_star
					A33[currentA] = a33_star
					
					acceptance_probability_rotation += 1.0
					energy_particle = energy_particle_star
				end
			end
			
			energy_system += energy_particle
		
		end
			
		# Update sigma_translation and sigma_rotation based on acceptance probabilities.
		acceptance_probability_translation /= number_of_particles		
		if acceptance_probability_translation <= acceptance_probability_target
			sigma_translation *= 0.95
		else
			sigma_translation = min(1.05 * sigma_translation, sigma_translation_ub)
		end
		
		if particle_type != "sphere"
			acceptance_probability_rotation /= number_of_particles
			if acceptance_probability_rotation <= acceptance_probability_target
				sigma_rotation *= 0.95
			else
				sigma_rotation = min(1.05 * sigma_rotation, sigma_rotation_ub)
			end
		end
	end
	
	return (X, Y, Z, Q0, Q1, Q2, Q3, A11, A12, A13, A21, A22, A23, A31, A32, A33, sigma_translation, sigma_rotation)
end
