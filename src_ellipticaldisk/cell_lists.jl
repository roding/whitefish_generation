function cell_lists(	X::Array{Float64,1},
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
					number_of_cells_x::Int64,
					number_of_cells_y::Int64,
					number_of_cells_z::Int64,
					cell_overlap::Float64)
		
	# Elliptical disk parameters.
	const number_of_particles::Int64 = length(X)
	
	const RMAX = Array(Float64, number_of_particles)
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
	
	# Compute bounding box for all elliptical disks.
	lbx_bounding_box = Array(Float64, number_of_particles)
	ubx_bounding_box = Array(Float64, number_of_particles)
	lby_bounding_box = Array(Float64, number_of_particles)
	uby_bounding_box = Array(Float64, number_of_particles)
	lbz_bounding_box = Array(Float64, number_of_particles)
	ubz_bounding_box = Array(Float64, number_of_particles)
	
	for current_particle = 1:number_of_particles
		(	lbx_bounding_box[current_particle], 
			ubx_bounding_box[current_particle],
			lby_bounding_box[current_particle],
			uby_bounding_box[current_particle],
			lbz_bounding_box[current_particle],
			ubz_bounding_box[current_particle]) = bounding_box(	X[current_particle],
															Y[current_particle],
															Z[current_particle],
															THETA1[current_particle],
															THETA2[current_particle],
															THETA3[current_particle],
															R1[current_particle], 
															R2[current_particle])
	end

	# Create cell dimension data structures.
	const cell_bounds_x::Array{Float64, 1} = linspace(0.0, Lx, number_of_cells_x + 1)
	#const cell_centers_x::Array{Float64, 1} = 0.5 * cell_bounds_x[1:end-1] + 0.5 * cell_bounds_x[2:end]
	const lbx_cell::Array{Float64, 1} = cell_bounds_x[1:end-1]
	const ubx_cell::Array{Float64, 1} = cell_bounds_x[2:end]
	
	const cell_bounds_y::Array{Float64, 1} = linspace(0.0, Ly, number_of_cells_y + 1)
	#const cell_centers_y::Array{Float64, 1} = 0.5 * cell_bounds_y[1:end-1] + 0.5 * cell_bounds_y[2:end]
	const lby_cell::Array{Float64, 1} = cell_bounds_y[1:end-1]
	const uby_cell::Array{Float64, 1} = cell_bounds_y[2:end]
	
	const cell_bounds_z::Array{Float64, 1} = linspace(0.0, Lz, number_of_cells_z + 1)
	#const cell_centers_z::Array{Float64, 1} = 0.5 * cell_bounds_z[1:end-1] + 0.5 * cell_bounds_z[2:end]
	const lbz_cell::Array{Float64, 1} = cell_bounds_z[1:end-1]
	const ubz_cell::Array{Float64, 1} = cell_bounds_z[2:end]
	
	# Initialize cell list data structure.
	cell_lists = Array{Any}(number_of_cells_x, number_of_cells_y, number_of_cells_z)
	for current_cell_x = 1:number_of_cells_x
		for current_cell_y = 1:number_of_cells_y
			for current_cell_z = 1:number_of_cells_z
				cell_lists[current_cell_x, current_cell_y, current_cell_z] = Array{Int64}(0)
			end
		end
	end
	
	# Compute cell lists.
	#x_cell_vertex::Array{Float64, 1} = zeros(8)
	#y_cell_vertex::Array{Float64, 1} = zeros(8)
	#z_cell_vertex::Array{Float64, 1} = zeros(8)
	
	#vx::Array{Float64, 1} = zeros(8)
	#vy::Array{Float64, 1} = zeros(8)
	#vz::Array{Float64, 1} = zeros(8)
	#w3::Array{Float64, 1} = zeros(8)
	
	for current_cell_x = 1:number_of_cells_x
		for current_cell_y = 1:number_of_cells_y
			for current_cell_z = 1:number_of_cells_z
				# Compute vertices (corners) of current cell.
			#	vertex_count = 0
			#	for bx = current_cell_x:current_cell_x+1
			#		for by = current_cell_y:current_cell_y+1
			#			for bz = current_cell_z:current_cell_z+1
			#				vertex_count += 1
			#				if bx == current_cell_x
			#					x_cell_vertex[vertex_count] = cell_bounds_x[bx] - cell_overlap
			#				else
			#					x_cell_vertex[vertex_count] = cell_bounds_x[bx] + cell_overlap
			#				end
			#				
			#				y_cell_vertex[vertex_count] = cell_bounds_y[by]
			#				if by == current_cell_y
			#					y_cell_vertex[vertex_count] = cell_bounds_y[by] - cell_overlap
			#				else
			#					y_cell_vertex[vertex_count] = cell_bounds_y[by] + cell_overlap
			#				end
			#				
			#				z_cell_vertex[vertex_count] = cell_bounds_z[bz]
			#				if bz == current_cell_z
			#					z_cell_vertex[vertex_count] = cell_bounds_z[bz] - cell_overlap
			#				else
			#					z_cell_vertex[vertex_count] = cell_bounds_z[bz] + cell_overlap
			#				end
			#			end
			#		end
			#	end
				
				# Check particle.
				for current_particle = 1:number_of_particles				




					if box_box_intersect(	lbx_bounding_box[current_particle],
										ubx_bounding_box[current_particle],
										lby_bounding_box[current_particle],
										uby_bounding_box[current_particle],
										lbz_bounding_box[current_particle],
										ubz_bounding_box[current_particle],
										lbx_cell[current_cell_x],
										ubx_cell[current_cell_x],
										lby_cell[current_cell_y],
										uby_cell[current_cell_y],
										lbz_cell[current_cell_z],
										ubz_cell[current_cell_z],
										
										
					# Coordinates of vertex positions relative to particle.
					vx = x_cell_vertex - X[current_particle]
					if mean(vx) < -0.5*Lx
						vx = x_cell_vertex - X[current_particle] + Lx
					elseif mean(vx) > 0.5*Lx
						vx = x_cell_vertex - X[current_particle] - Lx
					end
					
					vy = y_cell_vertex - Y[current_particle]
					if mean(vy) < -0.5*Ly
						vy = y_cell_vertex - Y[current_particle] + Ly
					elseif mean(vy) > 0.5*Ly
						vy = y_cell_vertex - Y[current_particle] - Ly
					end
					
					vz = z_cell_vertex - Z[current_particle]
					if mean(vz) < -0.5*Lz
						vz = z_cell_vertex - Z[current_particle] + Lz
					elseif mean(vz) > 0.5*Lz
						vz = z_cell_vertex - Z[current_particle] - Lz
					end
						
					w3 = q31[current_particle] * vx + q32[current_particle] * vy + q33[current_particle] * vz	
						
					if sign(maximum(w3)) - sign(minimum(w3)) >= 1.0 # Different sides of plane or one point intersecting plane
						if minimum(vx.^2 + vy.^2 + vz.^2) <= RMAX[current_particle]^2
							push!(cell_lists[current_cell_x, current_cell_y, current_cell_z], current_particle)
						end
					end
				end
			end
		end
	end
	return (cell_lists, cell_bounds_x, cell_bounds_y, cell_bounds_z)
end