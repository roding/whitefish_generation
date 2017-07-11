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
						sigma_rotation_ub::Float64,
						delta_phi::Float64,
						phi_target::Float64,
						number_of_sweeps_ub::Int64)
			
	number_of_particles::Int64 = size(R, 1)
	
	acceptance_probability_target::Float64 = 0.25

	# Pre-computed max radii (i.e. radius of bounding sphere).
	RMAX::Array{Float64, 1} = zeros(number_of_particles)
	if particle_type == "sphere"
		for current_particle = 1:number_of_particles
			RMAX[current_particle] = R[current_particle, 1]
		end
	elseif particle_type == "ellipsoid"
		for current_particle = 1:number_of_particles
			RMAX[current_particle] = maximum( R[current_particle, 1:3] )
		end
	elseif particle_type == "cuboid"
		for current_particle = 1:number_of_particles
			RMAX[current_particle] = sqrt(R[current_particle, 1]^2 + R[current_particle, 2]^2 + R[current_particle, 3]^2)
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
	
	phi::Float64 = 0.0
	if particle_type == "sphere"
		phi = sum(4.0 * pi / 3.0 * R[:, 1].^3) / (Lx * Ly * Lz)
	elseif particle_type == "ellipsoid"
		phi = sum(4.0 * pi / 3.0 * R[:, 1] .* R[:, 2] .* R[:, 3]) / (Lx * Ly * Lz)
	elseif particle_type == "cuboid"
		phi = sum(8.0 * R[:, 1] .* R[:, 2] .* R[:, 3]) / (Lx * Ly * Lz)
	end
	
	is_converged::Bool = false
	phi_prim::Float64 = 0.0
	Lx_prim::Float64 = 0.0
	Ly_prim::Float64 = 0.0
	Lz_prim::Float64 = 0.0
	
	X_prim::Array{Float64, 1} = zeros(number_of_particles)
	Y_prim::Array{Float64, 1} = zeros(number_of_particles)
	Z_prim::Array{Float64, 1} = zeros(number_of_particles)
	Q0_prim::Array{Float64, 1} = zeros(number_of_particles)
	Q1_prim::Array{Float64, 1} = zeros(number_of_particles)
	Q2_prim::Array{Float64, 1} = zeros(number_of_particles)
	Q3_prim::Array{Float64, 1} = zeros(number_of_particles)
	A11_prim::Array{Float64, 1} = zeros(number_of_particles)
	A12_prim::Array{Float64, 1} = zeros(number_of_particles)
	A13_prim::Array{Float64, 1} = zeros(number_of_particles)
	A21_prim::Array{Float64, 1} = zeros(number_of_particles)
	A22_prim::Array{Float64, 1} = zeros(number_of_particles)
	A23_prim::Array{Float64, 1} = zeros(number_of_particles)
	A31_prim::Array{Float64, 1} = zeros(number_of_particles)
	A32_prim::Array{Float64, 1} = zeros(number_of_particles)
	A33_prim::Array{Float64, 1} = zeros(number_of_particles)
	
	while !is_converged && phi < phi_target
		phi_prim = phi + delta_phi
		println(phi_prim)

		Lx_prim = (phi / phi_prim)^(1/3) * Lx
		Ly_prim = (phi / phi_prim)^(1/3) * Ly
		Lz_prim = (phi / phi_prim)^(1/3) * Lz
		
		X_prim = Lx_prim / Lx * X
		Y_prim = Ly_prim / Ly * Y
		Z_prim = Lz_prim / Lz * Z
		
		Q0_prim = Q0
		Q1_prim = Q1
		Q2_prim = Q2
		Q3_prim = Q3
		A11_prim = A11
		A12_prim = A12
		A13_prim = A13
		A21_prim = A21
		A22_prim = A22
		A23_prim = A23
		A31_prim = A31
		A32_prim = A32
		A33_prim = A33
	
		energy_system = 1.0
		current_sweep = 0
		while energy_system > 0.0 && current_sweep < number_of_sweeps_ub
			current_sweep += 1
			#println(join(["   Sweep ", string(current_sweep)]))
			
			acceptance_probability_translation = 0.0
			acceptance_probability_rotation = 0.0
			
			energy_system = 0.0
			
			for currentA = 1:number_of_particles
				# Compute current local energy.
				energy_particle = 0.0
				for currentB = [1:currentA-1 ; currentA+1:number_of_particles]
					xAB = signed_distance_mod(X_prim[currentA], X_prim[currentB], Lx_prim)
					yAB = signed_distance_mod(Y_prim[currentA], Y_prim[currentB], Ly_prim)
					zAB = signed_distance_mod(Z_prim[currentA], Z_prim[currentB], Lz_prim)

					if xAB^2 + yAB^2 + zAB^2 < (RMAX[currentA] + RMAX[currentB])^2
						if particle_type == "sphere"
							overlapfun = (RMAX[currentA] + RMAX[currentB])^2 - (xAB^2 + yAB^2 + zAB^2)
							energy_particle += overlapfun
						elseif particle_type == "ellipse"
							overlapfun = overlap_ellipse(xAB, yAB, zAB, A11_prim[currentA], A12_prim[currentA], A13_prim[currentA], A21_prim[currentA], A22_prim[currentA], A23_prim[currentA], A31_prim[currentA], A32_prim[currentA], A33_prim[currentA], A11_prim[currentB], A12_prim[currentB], A13_prim[currentB], A21_prim[currentB], A22_prim[currentB], A23_prim[currentB], A31_prim[currentB], A32_prim[currentB], A33_prim[currentB])
						
							if overlapfun < 1.0
								energy_particle += (1.0 - overlapfun)^2
							end
						elseif particle_type == "ellipsoid"
							overlapfun = overlap_ellipsoid(xAB, yAB, zAB, A11_prim[currentA], A12_prim[currentA], A13_prim[currentA], A21_prim[currentA], A22_prim[currentA], A23_prim[currentA], A31_prim[currentA], A32_prim[currentA], A33_prim[currentA], A11_prim[currentB], A12_prim[currentB], A13_prim[currentB], A21_prim[currentB], A22_prim[currentB], A23_prim[currentB], A31_prim[currentB], A32_prim[currentB], A33_prim[currentB], R[currentA, 1]^2 * R[currentA, 2]^2 * R[currentA, 3]^2)
						
							if overlapfun < 1.0
								energy_particle += (1.0 - overlapfun)^2
							end
						elseif particle_type == "cuboid"
							overlapfun = overlap_cuboid(xAB, yAB, zAB, A11_prim[currentA], A12_prim[currentA], A13_prim[currentA], A21_prim[currentA], A22_prim[currentA], A23_prim[currentA], A31_prim[currentA], A32_prim[currentA], A33_prim[currentA], A11_prim[currentB], A12_prim[currentB], A13_prim[currentB], A21_prim[currentB], A22_prim[currentB], A23_prim[currentB], A31_prim[currentB], A32_prim[currentB], A33_prim[currentB], R[currentA, 1], R[currentA, 2], R[currentA, 3], R[currentB, 1], R[currentB, 2], R[currentB, 3])

							energy_particle += overlapfun
						end
						
					end
				end
				#println(energy_particle)

				# Generate random proposal position and compute new local energy with translation.
				(x_star, y_star, z_star) = generate_proposal_position(X_prim[currentA], Y_prim[currentA], Z_prim[currentA], Lx_prim, Ly_prim, Lz_prim, sigma_translation)
				energy_particle_star = 0.0
				for currentB = [1:currentA-1;currentA+1:number_of_particles]
					xAB = signed_distance_mod(x_star, X_prim[currentB], Lx_prim)
					yAB = signed_distance_mod(y_star, Y_prim[currentB], Ly_prim)
					zAB = signed_distance_mod(z_star, Z_prim[currentB], Lz_prim)

					if xAB^2 + yAB^2 + zAB^2 < (RMAX[currentA] + RMAX[currentB])^2
						if particle_type == "sphere"
							overlapfun = (RMAX[currentA] + RMAX[currentB])^2 - (xAB^2 + yAB^2 + zAB^2)
							energy_particle_star += overlapfun
						elseif particle_type == "ellipse"
							overlapfun = overlap_ellipse(xAB, yAB, zAB, A11_prim[currentA], A12_prim[currentA], A13_prim[currentA], A21_prim[currentA], A22_prim[currentA], A23_prim[currentA], A31_prim[currentA], A32_prim[currentA], A33_prim[currentA], A11_prim[currentB], A12_prim[currentB], A13_prim[currentB], A21_prim[currentB], A22_prim[currentB], A23_prim[currentB], A31_prim[currentB], A32_prim[currentB], A33_prim[currentB])
						
							if overlapfun < 1.0
								energy_particle_star += (1.0 - overlapfun)^2
							end
						elseif particle_type == "ellipsoid"
							overlapfun = overlap_ellipsoid(xAB, yAB, zAB, A11_prim[currentA], A12_prim[currentA], A13_prim[currentA], A21_prim[currentA], A22_prim[currentA], A23_prim[currentA], A31_prim[currentA], A32_prim[currentA], A33_prim[currentA], A11_prim[currentB], A12_prim[currentB], A13_prim[currentB], A21_prim[currentB], A22_prim[currentB], A23_prim[currentB], A31_prim[currentB], A32_prim[currentB], A33_prim[currentB], R[currentA, 1]^2 * R[currentA, 2]^2 * R[currentA, 3]^2)
						
							if overlapfun < 1.0
								energy_particle_star += (1.0 - overlapfun)^2
							end
						elseif particle_type == "cuboid"
							overlapfun = overlap_cuboid(xAB, yAB, zAB, A11_prim[currentA], A12_prim[currentA], A13_prim[currentA], A21_prim[currentA], A22_prim[currentA], A23_prim[currentA], A31_prim[currentA], A32_prim[currentA], A33_prim[currentA], A11_prim[currentB], A12_prim[currentB], A13_prim[currentB], A21_prim[currentB], A22_prim[currentB], A23_prim[currentB], A31_prim[currentB], A32_prim[currentB], A33_prim[currentB], R[currentA, 1], R[currentA, 2], R[currentA, 3], R[currentB, 1], R[currentB, 2], R[currentB, 3])
							
							energy_particle_star += overlapfun
						end
						
					end
				end
				#println(energy_particle_star)
				
				if energy_particle_star <= energy_particle
					X_prim[currentA] = x_star
					Y_prim[currentA] = y_star
					Z_prim[currentA] = z_star
					
					acceptance_probability_translation += 1.0
					energy_particle = energy_particle_star
				end
				
				# Generate random proposal orientation and compute new local energy with rotation.
				if particle_type != "sphere"
					(q0_star, q1_star, q2_star, q3_star) = generate_proposal_orientation(Q0_prim[currentA], Q1_prim[currentA], Q2_prim[currentA], Q3_prim[currentA], sigma_rotation)
					
					if particle_type == "ellipse"
						(a11_star, a12_star, a13_star, a21_star, a22_star, a23_star, a31_star, a32_star, a33_star) = characteristic_matrix_ellipse(q0_star, q1_star, q2_star, q3_star, R[currentA, 1], R[currentA, 2])
					elseif particle_type == "ellipsoid"
						(a11_star, a12_star, a13_star, a21_star, a22_star, a23_star, a31_star, a32_star, a33_star) = characteristic_matrix_ellipsoid(q0_star, q1_star, q2_star, q3_star, R[currentA, 1], R[currentA, 2], R[currentA, 3])
					elseif particle_type == "cuboid"
						(a11_star, a12_star, a13_star, a21_star, a22_star, a23_star, a31_star, a32_star, a33_star) = rotation_matrix(q0_star, q1_star, q2_star, q3_star)
					end
					
					energy_particle_star = 0.0
					for currentB = [1:currentA-1;currentA+1:number_of_particles]
						xAB = signed_distance_mod(X_prim[currentA], X_prim[currentB], Lx_prim)
						yAB = signed_distance_mod(Y_prim[currentA], Y_prim[currentB], Ly_prim)
						zAB = signed_distance_mod(Z_prim[currentA], Z_prim[currentB], Lz_prim)

						if xAB^2 + yAB^2 + zAB^2 < (RMAX[currentA] + RMAX[currentB])^2
							if particle_type == "ellipse"
								overlapfun = overlap_ellipse(xAB, yAB, zAB, a11_star, a12_star, a13_star, a21_star, a22_star, a23_star, a31_star, a32_star, a33_star, A11_prim[currentB], A12_prim[currentB], A13_prim[currentB], A21_prim[currentB], A22_prim[currentB], A23_prim[currentB], A31_prim[currentB], A32_prim[currentB], A33_prim[currentB])
							
								if overlapfun < 1.0
									energy_particle_star += (1.0 - overlapfun)^2
								end
							elseif particle_type == "ellipsoid"
								overlapfun = overlap_ellipsoid(xAB, yAB, zAB, a11_star, a12_star, a13_star, a21_star, a22_star, a23_star, a31_star, a32_star, a33_star, A11_prim[currentB], A12_prim[currentB], A13_prim[currentB], A21_prim[currentB], A22_prim[currentB], A23_prim[currentB], A31_prim[currentB], A32_prim[currentB], A33_prim[currentB], R[currentA, 1]^2 * R[currentA, 2]^2 * R[currentA, 3]^2)
							
								if overlapfun < 1.0
									energy_particle_star += (1.0 - overlapfun)^2
								end
							elseif particle_type == "cuboid"
								overlapfun = overlap_cuboid(xAB, yAB, zAB, a11_star, a12_star, a13_star, a21_star, a22_star, a23_star, a31_star, a32_star, a33_star, A11[currentB], A12_prim[currentB], A13_prim[currentB], A21_prim[currentB], A22_prim[currentB], A23_prim[currentB], A31_prim[currentB], A32_prim[currentB], A33_prim[currentB], R[currentA, 1], R[currentA, 2], R[currentA, 3], R[currentB, 1], R[currentB, 2], R[currentB, 3])
							
								energy_particle_star += overlapfun
							end
							
						end
					end
					
					if energy_particle_star <= energy_particle
						Q0_prim[currentA] = q0_star
						Q1_prim[currentA] = q1_star
						Q2_prim[currentA] = q2_star
						Q3_prim[currentA] = q3_star
						
						A11_prim[currentA] = a11_star
						A12_prim[currentA] = a12_star
						A13_prim[currentA] = a13_star
						A21_prim[currentA] = a21_star
						A22_prim[currentA] = a22_star
						A23_prim[currentA] = a23_star
						A31_prim[currentA] = a31_star
						A32_prim[currentA] = a32_star
						A33_prim[currentA] = a33_star
						
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
		
		if energy_system == 0.0
			phi = phi_prim
			Lx = Lx_prim
			Ly = Ly_prim
			Lz = Lz_prim
			
			X = X_prim
			Y = Y_prim
			Z = Z_prim
			
			Q0 = Q0_prim
			Q1 = Q1_prim
			Q2 = Q2_prim
			Q3 = Q3_prim
			A11 = A11_prim
			A12 = A12_prim
			A13 = A13_prim
			A21 = A21_prim
			A22 = A22_prim
			A23 = A23_prim
			A31 = A31_prim
			A32 = A32_prim
			A33 = A33_prim
		else
			is_converged = true
		end
		println(join(["   Current volume fraction: ", string(phi)]))
	end
	
	return (Lx, Ly, Lz, X, Y, Z, Q0, Q1, Q2, Q3, A11, A12, A13, A21, A22, A23, A31, A32, A33, sigma_translation, sigma_rotation)
end
