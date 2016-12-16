#workspace()

include("io/get_version.jl")
include("io/print_header.jl")
include("io/read_xml_system.jl")
include("io/read_xml_input.jl")
include("io/print_simulation_stats.jl")

function run()
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
	
	# Error checking for input file specification.
	if !input_file_specified
		println("No input file specified or specified input does not exist. Aborting.")
	end
	
	# Print text header.
	if !silent_mode
		print_header()
	end
	
	# Process input file.
	(particle_system_path::String, output_path::String, D0::Float64, deltat_coarse::Float64, number_of_time_points_coarse::Int64, number_of_time_points_fine_per_coarse::Int64, number_of_diffusers::Int64, number_of_cells_x::Int64, number_of_cells_y::Int64, number_of_cells_z::Int64) = read_xml_input(input_file_path)
	
	# Print simulation stats.
	if !silent_mode
		print_simulation_stats()
	end
	
	
	
	
	
	
	
	
	
	
	# Print output information.
	if !silent_mode
		println(join(("Output written to ", output_path)))
	end
	
	nothing
end

run()