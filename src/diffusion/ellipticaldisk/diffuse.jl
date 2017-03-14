include("cell_lists.jl")

@inbounds function diffuse(	X::Array{Float64,1},
							Y::Array{Float64,1},
							Z::Array{Float64,1},
							THETA1::Array{Float64,1},
							THETA2::Array{Float64,1},
							THETA3::Array{Float64,1},
							R1::Array{Float64,1}, 
							R2::Array{Float64,1},
							Lx::Float64,
							Ly::Float64,
							Lz::Float64,
							D0::Float64,
							deltat_coarse::Float64,
							number_of_time_points_coarse::Int64,
							number_of_time_points_fine_per_coarse::Int64,
							number_of_diffusers::Int64,
							number_of_cells_x::Int64,
							number_of_cells_y::Int64,
							number_of_cells_z::Int64,
							silent_mode::Bool)
		
	if !silent_mode
		println("Preparing simulation...")
	end
	# Elliptical disk parameters.
	const number_of_particles::Int64 = length(X)
	const RMAX = Array(Float64,number_of_particles)
	for current_particle = 1:number_of_particles
		RMAX[current_particle] = maximum((R1[current_particle], R2[current_particle]))
	end
	
	const c1::Array{Float64,1} = cos(THETA1)
	const c2::Array{Float64,1} = cos(THETA2)
	const c3::Array{Float64,1} = cos(THETA3)
	const s1::Array{Float64,1} = sin(THETA1)
	const s2::Array{Float64,1} = sin(THETA2)
	const s3::Array{Float64,1} = sin(THETA3)
	
	const q11::Array{Float64,1} = c2.*c3
	const q12::Array{Float64,1} = -c2.*s3
	const q13::Array{Float64,1} = s2
	const q21::Array{Float64,1} = c1.*s3 + c3.*s1.*s2
	const q22::Array{Float64,1} = c1.*c3 - s1.*s2.*s3
	const q23::Array{Float64,1} = -c2.*s1
	const q31::Array{Float64,1} = s1.*s3 - c1.*c3.*s2
	const q32::Array{Float64,1} = c3.*s1 + c1.*s2.*s3
	const q33::Array{Float64,1} = c1.*c2

	# Standard deviation of Gaussian jumps.
	const deltat_fine::Float64 = deltat_coarse / convert(Float64, number_of_time_points_fine_per_coarse)
	const sigma::Float64 = sqrt(2 * D0 * deltat_fine)
	
	# Create cell lists.
	if !silent_mode
		println("Creating cell lists...")
	end
	cell_overlap::Float64 = 6 * sigma
	lists = cell_lists(	X,
						Y,
						Z,
						THETA1,
						THETA2,
						THETA3,
						R1, 
						R2,
						Lx,
						Ly,
						Lz,
						number_of_cells_x,
						number_of_cells_y,
						number_of_cells_z,
						cell_overlap)

	# Simulate diffusion.
	if !silent_mode
		println("Starting simulation...")
	end
	current_particle::Int64 = 0
	is_ok::Bool = true
	
	x::Float64 = 0.0
	y::Float64 = 0.0
	z::Float64 = 0.0
	x_abs::Float64 = 0.0
	y_abs::Float64 = 0.0
	z_abs::Float64 = 0.0
	x_star::Float64 = 0.0
	y_star::Float64 = 0.0
	z_star::Float64 = 0.0
	deltax::Float64 = 0.0
	deltay::Float64 = 0.0
	deltaz::Float64 = 0.0
	vx::Float64 = 0.0
	vy::Float64 = 0.0
	vz::Float64 = 0.0
	vx_star::Float64 = 0.0
	vy_star::Float64 = 0.0
	vz_star::Float64 = 0.0
	w1::Float64 = 0.0
	w1_star::Float64 = 0.0
	w2::Float64 = 0.0
	w2_star::Float64 = 0.0
	w3::Float64 = 0.0
	w3_star::Float64 = 0.0
	w1_intersection::Float64 = 0.0
	w2_intersection::Float64 = 0.0
	alpha::Float64 = 0.0
	trajectory_x::Array{Float64} = zeros(number_of_time_points_coarse)
	trajectory_y::Array{Float64} = zeros(number_of_time_points_coarse)
	trajectory_z::Array{Float64} = zeros(number_of_time_points_coarse)
	msd_x::Array{Float64} = zeros(number_of_time_points_coarse)
	msd_y::Array{Float64} = zeros(number_of_time_points_coarse)
	msd_z::Array{Float64} = zeros(number_of_time_points_coarse)
	D0_empirical::Float64 = 0.0
	
	t_start_diffusion::Float64 = convert(Float64, time_ns()) / 1e9
	t_elapsed_diffusion::Float64 = 0.0
	#fraction_done::Float64 = 0.0
	#percent_done::Float64 = 0.0
	#t_remaining_diffusion::Float64 = 0.0
	chunk::Int64 = 0
	println("   Elapsed time (hh:mm:ss)   Done (%)   Est. time remaining (hh:mm:ss)")
	for current_diffuser = 1:number_of_diffusers
		#println(current_diffuser)
		
		# By definition, a random position end up in void almost surely (w.p. 1),
		# so as initial position we just pick one randomly.
		x = Lx * rand()
		y = Ly * rand()
		z = Lz * rand()
		
		x_abs = x
		y_abs = y
		z_abs = z
		
		trajectory_x[1] = x_abs
		trajectory_y[1] = y_abs
		trajectory_z[1] = z_abs
		
		# Starting diffusion.
		for current_time_coarse = 2:number_of_time_points_coarse
			for current_time_fine = 1:number_of_time_points_fine_per_coarse
				# Determine current cell.
				current_cell_x = convert(Int64, ceil(x / Lx * convert(Float64, number_of_cells_x)))
				current_cell_y = convert(Int64, ceil(y / Ly * convert(Float64, number_of_cells_y)))
				current_cell_z = convert(Int64, ceil(z / Lz * convert(Float64, number_of_cells_z)))
				#println((current_cell_x,current_cell_y,current_cell_z))
				#println(size(lists))
				number_of_particles_current_cell = length(lists[current_cell_x, current_cell_y, current_cell_z])
				#println(number_of_particles_current_cell)
				#return :TEST	
				list = lists[current_cell_x, current_cell_y, current_cell_z]
				# Random proposal displacement.
				deltax = sigma * randn()
				deltay = sigma * randn()
				deltaz = sigma * randn()
				
				# Calculate proposed new position.
				x_star = x + deltax
				y_star = y + deltay
				z_star = z + deltaz
				
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
				
				# Check for diffuser-particle intersections.
				current_particle = 0
				is_ok = true
				while is_ok & ( current_particle < number_of_particles_current_cell )
					current_particle += 1
						
					# Coordinates of current diffuser position relative to particle.
					vx = x - X[list[current_particle]]
					if vx < -0.5*Lx
						vx = vx + Lx
					elseif vx > 0.5*Lx
						vx = vx - Lx
					end
					
					vy = y - Y[list[current_particle]]
					if vy < -0.5*Ly
						vy = vy + Ly
					elseif vy > 0.5*Ly
						vy = vy - Ly
					end
					
					vz = z - Z[list[current_particle]]
					if vz < -0.5*Lz
						vz = vz + Lz
					elseif vz > 0.5*Lz
						vz = vz - Lz
					end
						
					# Coordinates of candidate diffuser position relative to particle.
					vx_star = vx + deltax
					vy_star = vy + deltay
					vz_star = vz + deltaz
					
					#if ( vx^2 + vy^2 + vz^2 < RMAX[list[current_particle]]^2 ) || ( vx_star^2 + #vy_star^2 + vz_star^2 < RMAX[list[current_particle]]^2 )
						#println("Within bounding sphere of ellipse...")
						
						# Coordinate of current and candidate diffuser position relative to particle 
						# in the intrinsic coordinates of the particle, orthogonal to ellipse plane.
						w3 = q31[list[current_particle]] * vx + q32[list[current_particle]] * vy + q33[list[current_particle]] * vz
						w3_star = q31[list[current_particle]] * vx_star + q32[list[current_particle]] * vy_star + q33[list[current_particle]] * vz_star
		
						if sign(w3_star) != sign(w3)
							#println("   On different sides of ellipse...")
							#println((z, z_star, w3, w3_star))
							# Coordinate of current and candidate diffuser position relative to particle 
							# in the intrinsic coordinates of the particle, parallel to ellipse plane.
							w1 = q11[list[current_particle]] * vx + q12[list[current_particle]] * vy + q13[list[current_particle]] * vz
							w1_star = q11[list[current_particle]] * vx_star + q12[list[current_particle]] * vy_star + q13[list[current_particle]] * vz_star
							w2 = q21[list[current_particle]] * vx + q22[list[current_particle]] * vy + q23[list[current_particle]] * vz
							w2_star = q21[list[current_particle]] * vx_star + q22[list[current_particle]] * vy_star + q23[list[current_particle]] * vz_star
							
							# Compute intersection point between diffuser trajectory and ellipse plane.
							alpha = - w3_star / (w3 - w3_star)
							w1_intersection = alpha * w1 + (1 - alpha) * w1_star
							w2_intersection = alpha * w2 + (1 - alpha) * w2_star
							
							# If intrinsic elliptical planar coordinates lie satisfy ellipse equation we have hit the ellipse.
							if (w1_intersection/R1[list[current_particle]])^2 + (w2_intersection/R2[list[current_particle]])^2 < 1.0
								#println("      Intersection point within ellipse...")
								#println(join(("      ", string((current_time, current_particle)))))
								is_ok = false
							end
						end
					#end
				end
			
				if is_ok
					x = x_star
					y = y_star
					z = z_star
					
					x_abs = x_abs + deltax
					y_abs = y_abs + deltay
					z_abs = z_abs + deltaz
					
					D0_empirical = D0_empirical + deltax^2 + deltay^2 + deltaz^2
				end
			end
			
			trajectory_x[current_time_coarse] = x_abs
			trajectory_y[current_time_coarse] = y_abs
			trajectory_z[current_time_coarse] = z_abs
			
			msd_x[current_time_coarse] += (trajectory_x[current_time_coarse] - trajectory_x[1])^2
			msd_y[current_time_coarse] += (trajectory_y[current_time_coarse] - trajectory_y[1])^2
			msd_z[current_time_coarse] += (trajectory_z[current_time_coarse] - trajectory_z[1])^2
		end
		
		t_elapsed_diffusion = convert(Float64, time_ns()) / 1e9 - t_start_diffusion
		if !silent_mode && convert(Int64, floor(t_elapsed_diffusion / 10.0)) > chunk
			print_progress_diffusion(t_elapsed_diffusion, current_diffuser, number_of_diffusers)
		end
		chunk = convert(Int64, floor(t_elapsed_diffusion / 10.0))
	end
	print_progress_diffusion(t_elapsed_diffusion, number_of_diffusers, number_of_diffusers)
	
	msd_x = msd_x ./ convert(Float64, number_of_diffusers)
	msd_y = msd_y ./ convert(Float64, number_of_diffusers)
	msd_z = msd_z ./ convert(Float64, number_of_diffusers)
	
	D0_empirical = D0_empirical / (3.0 * convert(Float64, number_of_diffusers * (number_of_time_points_coarse-1) * number_of_time_points_fine_per_coarse) * 2.0 * deltat_fine)
		
	return (msd_x, msd_y, msd_z, D0_empirical)
end
