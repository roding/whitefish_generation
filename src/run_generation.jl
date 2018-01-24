include("file_io/read_input.jl")
include("file_io/read_key.jl")
include("file_io/write_output.jl")
include("file_io/write_key.jl")

include("initialize_system.jl")
include("evolve_system.jl")
include("compress_system.jl")

include("generate_random_unit_quaternion.jl")
include("characteristic_matrix_ellipsoid.jl")
include("characteristic_matrix_ellipse.jl")
include("rotation_matrix.jl")

include("generate_proposal_position.jl")
include("generate_proposal_orientation.jl")
include("overlap_ellipsoid.jl")
include("overlap_ellipse.jl")

include("overlap_cuboid.jl")
include("overlap_cuboid_binary.jl")

include("signed_distance_mod.jl")
include("quaternion_mult.jl")

include("vtk/voxel_structure.jl")
include("vtk/write_vtk.jl")

function run_generation()
	# Inititalization of random number generation device.
	random_seed::Int64 = convert(Int64, time_ns())
	srand(random_seed)

	# Start time.
	t_start_ns::Int64 = convert(Int64, time_ns())

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
		R::Array{Float64, 2},
		Lx::Float64,
		Ly::Float64,
		Lz::Float64,
		phi::Float64,
		phi_target::Float64,
		position_constraint_axis::String,
		position_constraint_lower::Float64,
		position_constraint_upper::Float64,
		orientation_axis::Array{Float64, 1},
		orientation_constraint_axis::Array{Float64, 1},
		orientation_constraint_lower::Float64,
		orientation_constraint_upper::Float64,
		sigma_translation_max::Float64,
		sigma_rotation_max::Float64,
		sigma_ratio::Float64,
		number_of_equilibration_sweeps::Int64,
		delta_phi::Float64,
		output_file_path::String) = read_input(input_file_path)

	# Initialize system.
	(	Lx,
		Ly,
		Lz,
		phi,
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
		A33::Array{Float64, 1}) = initialize_system(	particle_type,
														R,
														Lx,
														Ly,
														Lz,
														phi,
														position_constraint_axis,
														position_constraint_lower,
														position_constraint_upper,
														orientation_axis,
														orientation_constraint_axis,
														orientation_constraint_lower,
														orientation_constraint_upper)

	# Translation and rotation speeds.
	sigma_translation::Float64 = sigma_translation_max
	sigma_rotation::Float64 = 0.0
	if sigma_ratio >= 0.0
		sigma_rotation = sigma_ratio * sigma_translation
	else
		sigma_rotation = sigma_rotation_max
	end

	# Relax and equlibrate system.
	number_of_relaxation_sweeps_max::Int64 = typemax(Int64)
	(	X,
		Y,
		Z,
		Q0,
		Q1,
		Q2,
		Q3,
		A11,
		A12,
		A13,
		A21,
		A22,
		A23,
		A31,
		A32,
		A33,
		sigma_translation,
		sigma_rotation,
		is_relaxed::Bool) = evolve_system(	particle_type,
										R,
										Lx,
										Ly,
										Lz,
										X,
										Y,
										Z,
										Q0,
										Q1,
										Q2,
										Q3,
										A11,
										A12,
										A13,
										A21,
										A22,
										A23,
										A31,
										A32,
										A33,
										position_constraint_axis,
										position_constraint_lower,
										position_constraint_upper,
										orientation_axis,
										orientation_constraint_axis,
										orientation_constraint_lower,
										orientation_constraint_upper,
										sigma_translation,
										sigma_translation_max,
										sigma_rotation,
										sigma_rotation_max,
										sigma_ratio,
										number_of_relaxation_sweeps_max,
										number_of_equilibration_sweeps)

	# Compress system.
	number_of_relaxation_sweeps_max = 1000
	if particle_type != "ellipse"
		if phi_target > phi
			(	Lx,
				Ly,
				Lz,
				phi,
				X,
				Y,
				Z,
				Q0,
				Q1,
				Q2,
				Q3,
				A11,
				A12,
				A13,
				A21,
				A22,
				A23,
				A31,
				A32,
				A33,
				sigma_translation,
				sigma_rotation) = compress_system(	particle_type,
												R,
												Lx,
												Ly,
												Lz,
												phi_target,
												X,
												Y,
												Z,
												Q0,
												Q1,
												Q2,
												Q3,
												A11,
												A12,
												A13,
												A21,
												A22,
												A23,
												A31,
												A32,
												A33,
												position_constraint_axis,
												position_constraint_lower,
												position_constraint_upper,
												orientation_axis,
												orientation_constraint_axis,
												orientation_constraint_lower,
												orientation_constraint_upper,
												sigma_translation,
												sigma_translation_max,
												sigma_rotation,
												sigma_rotation_max,
												sigma_ratio,
												delta_phi,
												number_of_relaxation_sweeps_max)
		else
			println("Target volume fraction not larger than initial volume fraction. Not compressing system.")
		end
	else
		if phi_target > phi
			println("System with particle type 'ellipse' not eligible for compression. Not compressing system.")
		end
	end

	# Write output.
	t_exec_ns::Int64 = convert(Int64, time_ns()) - t_start_ns
	t_exec_s::Float64 = convert(Float64, t_exec_ns) / 1e9
	write_output(	output_file_path,
				particle_type,
				R,
				Lx,
				Ly,
				Lz,
				phi,
				X,
				Y,
				Z,
				Q0,
				Q1,
				Q2,
				Q3,
				t_exec_s)

	nothing
end

run_generation()
