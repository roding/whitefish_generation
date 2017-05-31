include("get_version.jl")

function print_header()
	version_string::String = get_version()
	println(join(("This is Whitefish generation module version ", version_string, ".")))
	println(join(("Executed on Julia version ", string(Sys.VERSION), " (", Sys.MACHINE, ").")))
	
	nothing	
end
