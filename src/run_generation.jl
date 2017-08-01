include("file_io/read_input.jl")
include("file_io/read_key.jl")

include("initialize_system.jl")

include("generate_random_unit_quaternion.jl")
include("characteristic_matrix_ellipsoid.jl")
include("characteristic_matrix_ellipse.jl")
include("rotation_matrix.jl")
include("relax_system.jl")
include("equilibrate_system.jl")
include("compress_system.jl")
include("generate_proposal_position.jl")
include("generate_proposal_orientation.jl")
include("overlap_ellipsoid.jl")
include("overlap_ellipse.jl")
include("overlap_cuboid.jl")
include("signed_distance_mod.jl")
include("quaternion_mult.jl")

function run_generation()
	# Inititalization of random number generation device.
	random_seed::Int64 = convert(Int64, time_ns())
	srand(random_seed)
	
	# Change current folder to the folder where this script lies.
	(program_file_dir::String, program_file_name::String) = splitdir(PROGRAM_FILE)
	program_file_dir = abspath(program_file_dir)
	cd(program_file_dir)

	# Assert that input is file and store path.
	input_file_path::String = ""
	if isfile(ARGS[1])
		input_file_path = ARGS[1]
	else
		println("No input file specified or specified input file does not exist. Aborting.")
		return nothing
	end
	
	# Read input from file.
	println(join(("Reading input from file ", input_file_path, "...")))
	(	particle_type::String, 
		number_of_particles::Int64, 
		R::Array{Float64, 2},
		Lx::Float64, 
		Ly::Float64, 
		Lz::Float64, 
		phi_initial::Float64, 
		phi_target::Float64,
		position_constraint_axis::String, 
		position_constraint_lower::Float64, 
		position_constraint_upper::Float64, 
		orientation_constraint_axis::Array{Float64, 1},
		orientation_constraint_lower::Float64,
		orientation_constraint_upper::Float64,
		number_of_equilibration_sweeps::Int64, 
		output_file_path::String) = read_input(input_file_path)
	
	println((Lx,Ly,Lz))
	println(particle_type)
	println(number_of_particles)
	println(R)
	println(position_constraint_axis)
	println(position_constraint_lower)
	println(position_constraint_upper)
	println(orientation_constraint_axis)
	println(orientation_constraint_lower)
	println(orientation_constraint_upper)
	println(number_of_equilibration_sweeps)
	println(output_file_path)
	
	
	
		
	return 0
	
	nothing
end

run_generation()