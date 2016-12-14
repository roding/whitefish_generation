include("create_cell_lists.jl")

function test_create_cell_lists()
	Lx::Float64 = 100.0
	Ly::Float64 = 100.0
	Lz::Float64 = 200.0
	
	number_of_particles::Int64 = 500
	X::Array{Float64,1} = Lx * rand(number_of_particles)
	Y::Array{Float64,1} = Ly * rand(number_of_particles)
	Z::Array{Float64,1} = Lz * rand(number_of_particles)
	THETA1::Array{Float64,1} = 2.0 * pi * rand(number_of_particles)
	THETA2::Array{Float64,1} = 2.0 * pi * rand(number_of_particles)
	THETA3::Array{Float64,1} = 2.0 * pi * rand(number_of_particles)
	R1::Array{Float64,1} = 10.0 * ones(number_of_particles)
	R2::Array{Float64,1} = 10.0 * ones(number_of_particles)
	
	number_of_cells_x::Int64 = 10
	number_of_cells_y::Int64 = 10
	number_of_cells_z::Int64 = 20
	cell_overlap::Float64 = 0.0 # same unit as simulation domain dimensions (µm).


	(cell_lists, cell_bounds_x, cell_bounds_y, cell_bounds_z) = create_cell_lists(	X,
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
																			
	#println(cell_lists)
	#println(cell_bounds_x)
	#println(cell_bounds_y)
	#println(cell_bounds_z)
	
	count = 0.0
	for i = 1:(number_of_cells_x*number_of_cells_y*number_of_cells_z)
		count += length(cell_lists[i])
	end
	println(count/(number_of_cells_x*number_of_cells_y*number_of_cells_z))
	
	nothing
end

test_create_cell_lists()