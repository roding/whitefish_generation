#workspace()

include("io/get_version.jl")
include("io/print_header.jl")
include("io/read_xml_system.jl")
include("io/read_xml_input.jl")

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
	
	
	
	
	
	
	
	
	
	
	
	nothing
end

run()