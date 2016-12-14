include("bounding_box.jl")
include("box_box_intersect.jl")
include("plane_box_intersect.jl")

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
	const lbx_cell::Array{Float64, 1} = cell_bounds_x[1:end-1] - cell_overlap
	const ubx_cell::Array{Float64, 1} = cell_bounds_x[2:end] + cell_overlap
	
	const cell_bounds_y::Array{Float64, 1} = linspace(0.0, Ly, number_of_cells_y + 1)
	const lby_cell::Array{Float64, 1} = cell_bounds_y[1:end-1] - cell_overlap
	const uby_cell::Array{Float64, 1} = cell_bounds_y[2:end] + cell_overlap
	
	const cell_bounds_z::Array{Float64, 1} = linspace(0.0, Lz, number_of_cells_z + 1)
	const lbz_cell::Array{Float64, 1} = cell_bounds_z[1:end-1] - cell_overlap
	const ubz_cell::Array{Float64, 1} = cell_bounds_z[2:end] + cell_overlap
	
	# Initialize cell list data structure.
	lists = Array{Array{Int64, 1}}(number_of_cells_x, number_of_cells_y, number_of_cells_z)
	for current_cell_x = 1:number_of_cells_x
		for current_cell_y = 1:number_of_cells_y
			for current_cell_z = 1:number_of_cells_z
				lists[current_cell_x, current_cell_y, current_cell_z] = Array{Int64}(0)
			end
		end
	end
	
	# Compute cell lists.
	for current_cell_x = 1:number_of_cells_x
		for current_cell_y = 1:number_of_cells_y
			for current_cell_z = 1:number_of_cells_z
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
										Lx,
										Ly,
										Lz)
						if plane_box_intersect(	X[current_particle],
												Y[current_particle],
												Z[current_particle],
												q31[current_particle],
												q32[current_particle],
												q33[current_particle],
												lbx_cell[current_cell_x],
												ubx_cell[current_cell_x], 
												lby_cell[current_cell_y],
												uby_cell[current_cell_y], 
												lbz_cell[current_cell_z],
												ubz_cell[current_cell_z],
												Lx,
												Ly,
												Lz)
							push!(lists[current_cell_x, current_cell_y, current_cell_z], current_particle)
						end
					end
				end
			end
		end
	end
	#println(lists)
	return lists
end