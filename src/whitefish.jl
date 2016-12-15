include("get_version.jl")
include("print_header.jl")
#include("ellipticaldisk/*.jl")

function whitefish()
	silent_mode::Bool = false
	
	
	if !silent_mode
		print_header()
	end
	
	
	
	
	nothing
end

whitefish()