include("io/get_version.jl")
include("io/print_header.jl")
include("io/read_xml_input_generation.jl")
include("io/write_xml_output_generation.jl")

#include("io/read_xml_system.jl")
#include("io/write_xml_output.jl")
#include("io/print_simulation_stats.jl")
#include("io/print_progress.jl")

include("generation/ellipticaldisk/generate.jl")

function wfrun_generation()
	# Change current folder to the folder where this script lies.
	(program_dir::String, program_file_name::String) = splitdir(PROGRAM_FILE)
	if program_dir != ""
		cd(program_dir)
	end
	
	# Process input arguments.
	silent_mode::Bool = false
	
	input_file_path::String = ""
	input_file_specified::Bool = false
	
	number_of_arguments::Int64 = length(ARGS)
	
	for current_argument = 1:number_of_arguments
		if ARGS[current_argument] == "-silent"
			silent_mode = true
		elseif isfile(ARGS[current_argument])
			input_file_path = ARGS[current_argument]
			input_file_specified = true
		end
	end
	
	# Print text header.
	if !silent_mode
		print_header()
	end
	
	# Error checking for input file specification.
	if !input_file_specified
		println("No input file specified or specified input does not exist. Aborting.")
	end

	# Process input file.
	if !silent_mode
		println(join(("Reading input from file ", input_file_path, "...")))
	end
	(Lx::Float64, Ly::Float64, Lz::Float64, particle_type::String, number_of_particles::Int64, lbz::Float64, ubz::Float64, ubangle::Float64, number_of_equilibration_sweeps::Int64, output_file_path::String) = read_xml_input_generation(input_file_path)
	
	# Print simulation stats.
	if !silent_mode
		print_simulation_stats_generation(particle_type, number_of_particles, Lx, Ly, Lz)
	end
	
	# Generate system
	
	# Write output.
	write_xml_output_generation(output_file_path, Lx, Ly, Lz, number_of_particles, X, Y, Z, THETA1, THETA2, THETA3, R1, R2)
	
	# Print output information.
	if !silent_mode
		println(join(("Output written to ", output_file_path, ".")))
		println("Finished.")
	end
	
	nothing
end

wfrun_generation()