workspace()

include("io/get_version.jl")
include("io/print_header.jl")
include("io/read_xml_system.jl")

#include("ellipticaldisk/*.jl")

function run()
	silent_mode::Bool = false
		
	if !silent_mode
		print_header()
	end
	
	println(ARGS)
	
	
	nothing
end

run()