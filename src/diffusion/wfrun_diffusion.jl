include("../io/get_version.jl")
include("../io/print_header.jl")
include("../io/read_xml_system.jl")
include("../io/read_xml_input.jl")
include("../io/write_xml_output.jl")
include("../io/print_simulation_stats.jl")
include("../io/print_progress.jl")

include("ellipticaldisk/simulate_diffusion.jl")

function wfrun_diffusion()
	# Change current folder to the folder where this script lies.
	(program_dir::String, script_file_name::String) = splitdir(PROGRAM_FILE)
	cd(program_dir)
	
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
	(particle_system_path::String, output_path::String, D0::Float64, deltat_coarse::Float64, number_of_time_points_coarse::Int64, number_of_time_points_fine_per_coarse::Int64, number_of_diffusers::Int64, number_of_cells_x::Int64, number_of_cells_y::Int64, number_of_cells_z::Int64) = read_xml_input(input_file_path)
	
	# Process particle system file.
	if !silent_mode
		println(join(("Reading particle system configuration from file ", particle_system_path, "...")))
	end
	(particles::String, X::Array{Float64, 1}, Y::Array{Float64, 1}, Z::Array{Float64, 1}, THETA1::Array{Float64, 1}, THETA2::Array{Float64, 1}, THETA3::Array{Float64, 1}, R1::Array{Float64, 1}, R2::Array{Float64, 1}, Lx::Float64, Ly::Float64, Lz::Float64) = read_xml_system(particle_system_path)
	
	number_of_particles::Int64	 = length(X)
		
	# Print simulation stats.
	if !silent_mode
		print_simulation_stats(particles, number_of_particles, number_of_diffusers, Lx, Ly, Lz, number_of_cells_x, number_of_cells_y, number_of_cells_z)
	end
	
	# Simulate diffusion.
	(msd_x::Array{Float64, 1}, msd_y::Array{Float64, 1}, msd_z::Array{Float64, 1}, D0_empirical::Float64) = simulate_diffusion(X, Y, Z, THETA1, THETA2, THETA3, R1, R2, Lx, Ly, Lz, D0, deltat_coarse, number_of_time_points_coarse, number_of_time_points_fine_per_coarse, number_of_diffusers, number_of_cells_x, number_of_cells_y, number_of_cells_z, silent_mode)	
	
	# Write output.
	write_xml_output(output_path, D0, D0_empirical, deltat_coarse, number_of_time_points_coarse, msd_x, msd_y, msd_z)
	
	# Print output information.
	if !silent_mode
		println(join(("Output written to ", output_path, ".")))
		println("Finished.")
	end
	
	nothing
end

wfrun_diffusion()