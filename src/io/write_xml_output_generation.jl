include("write_xml_key.jl")

function write_xml_output_generation(output_file_path::String)
	file_stream::IOStream = open(output_file_path, "w")
	
	@printf(file_stream, "%s", "<output_generation>\n")

	write_xml_key(file_stream, "key", 3.4)
	
	@printf(file_stream, "%s", "</output_generation>")
	
	close(file_stream)

	nothing
end

#write_xml_output_generation("test_output.xml")