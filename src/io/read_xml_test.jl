function read_xml_test()
	file_name::String = "particle_system.xml"
	file_stream::IOStream = open(file_name, "r")
	file_lines::Array{String, 1} = readlines(file_stream)
	close(file_stream)
	#stream_ = IOBuffer(string)
	
	
	
	println(file_lines)






	nothing
end

read_xml_test()