workspace()

include("io/get_version.jl")
include("io/print_header.jl")
include("io/read_xml_system.jl")

#include("ellipticaldisk/*.jl")

function run()
	silent_mode::Bool = false
	input_file_path::String = ""
	input_file_ok::Bool = false
	
	number_of_arguments::Int64 = length(ARGS)
	
	for current_argument = 1:number_of_arguments
		if ARGS[current_argument] == "silent"
			silent_mode = true
		elseif isfile(ARGS[current_argument])
			input_file_path = ARGS[current_argument]
			input_file_ok = true
		end
	end
	
	if !input_file_ok
		error("No input file specified.")
	else
		println("Input file specified.")
	end
	
	if !silent_mode
		print_header()
	end
	
	
	
	
	nothing
end

run()