include("io/get_version.jl")
include("io/print_header.jl")
include("io/read_xml_input_generation.jl")
include("io/read_xml_key.jl")



#include("io/read_xml_system.jl")
#include("io/write_xml_output.jl")
#include("io/print_simulation_stats.jl")
#include("io/print_progress.jl")

#include("generation/ellipticaldisk/generate.jl")

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
	(Lx::Float64, Ly::Float64, Lz::Float64, number_of_particles::Int64, lbz::Float64, ubz::Float64, ubangle::Float64, number_of_equilibration_sweeps::Int64, output_file_path::String) = read_xml_input_generation(input_file_path)
	
	# Print simulation stats.
#	if !silent_mode
#		print_simulation_stats(particles, number_of_particles, number_of_diffusers, Lx, Ly, Lz, number_of_cells_x, n#umber_of_cells_y, number_of_cells_z)
#	end
	
	# Simulate diffusion.
#	(msd_x::Array{Float64, 1}, msd_y::Array{Float64, 1}, msd_z::Array{Float64, 1}, D0_empirical::Float64) = s#imulate_diffusion(X, Y, Z, THETA1, THETA2, THETA3, R1, R2, Lx, Ly, Lz, D0, deltat_coarse, number_of_time_points_coarse, number_of_time_points_fine_per_coarse, number_of_diffusers, number_of_cells_x, number_of_cells_y, number_of_cells_z, silent_mode)	
	
	# Write output.
	write_xml_output_generation(output_path, D0, D0_empirical, deltat_coarse, number_of_time_points_coarse, msd_x, msd_y, msd_z)
	
	
	# Print output information.
	if !silent_mode
		println(join(("Output written to ", output_file_path, ".")))
		println("Finished.")
	end
	
	nothing
end

wfrun_generation()