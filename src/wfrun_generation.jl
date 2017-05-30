include("../io/read_xml_key.jl")
include("../io/write_xml_key.jl")

include("../io/read_xml_input_generation.jl")
include("../io/write_xml_output_generation.jl")

include("../io/print_header.jl")
include("../io/print_simulation_stats_generation.jl")

include("ellipticaldisk/generate.jl")

function wfrun_generation()
	# Inititalization of random number generation device.
	const random_seed::Int64 = convert(Int64, time_ns())
	srand(random_seed)
	
	# Change current folder to the folder where this script lies.
	(program_file_dir::String, program_file_name::String) = splitdir(PROGRAM_FILE)
	program_file_dir = abspath(program_file_dir)
	cd(program_file_dir)
	
	# Process input arguments.
	silent_mode::Bool = false
	
	input_generation_path::String = ""
	input_generation_specified::Bool = false
	
	number_of_arguments::Int64 = length(ARGS)
	
	for current_argument = 1:number_of_arguments
		if ARGS[current_argument] == "-silent"
			silent_mode = true
		elseif isfile(ARGS[current_argument])
			input_generation_path = ARGS[current_argument]
			input_generation_specified = true
		end
	end
	
	# Print text header.
	if !silent_mode
		print_header()
	end
	
	# Error checking for input file specification.
	if !input_generation_specified
		println("No input file specified or specified input does not exist. Aborting.")
		return nothing
	end

	# Process input file.
	if !silent_mode
		println(join(("Reading input from file ", input_generation_path, "...")))
	end
	(Lx::Float64, Ly::Float64, Lz::Float64, particle_type::String, number_of_particles::Int64, lbz::Float64, ubz::Float64, ubangle::Float64, R1::Array{Float64,1}, R2::Array{Float64,1}, number_of_equilibration_sweeps::Int64, output_generation_path::String) = read_xml_input_generation(input_generation_path)
	
	# Print simulation stats.
	if !silent_mode
		print_simulation_stats_generation(particle_type, number_of_particles, Lx, Ly, Lz)
	end
	
	# Generate system.
	acceptance_probability_target::Float64 = 0.25
	number_of_iterations_overlap_criterion::Int64 = 20
	t_start_ns::Int64 = convert(Int64, time_ns())
	(X::Array{Float64,1}, Y::Array{Float64,1}, Z::Array{Float64,1}, THETA1::Array{Float64,1}, THETA2::Array{Float64,1}, THETA3::Array{Float64,1}) = generate(R1, R2, Lx, Ly, Lz, lbz, ubz, ubangle, number_of_equilibration_sweeps, acceptance_probability_target, number_of_iterations_overlap_criterion, silent_mode)
	t_finish_ns::Int64 = convert(Int64, time_ns())
	t_exec::Float64 = convert(Float64, t_finish_ns - t_start_ns) / 1e9
	
	# Write output.
	write_xml_output_generation(output_generation_path, Lx, Ly, Lz, number_of_particles, X, Y, Z, THETA1, THETA2, THETA3, R1, R2, t_exec)
	
	# Print output information.
	if !silent_mode
		println(join(("Output written to ", output_generation_path, ".")))
		println("Finished.")
	end
	
	nothing
end

wfrun_generation()
