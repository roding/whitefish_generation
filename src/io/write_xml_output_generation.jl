function write_xml_output_generation(output_file_path::String)
	file_stream::IOStream = open(output_file_path, "w")
	
	@printf(file_stream, "%s", "<output_generation>\n")

	
	
	@printf(file_stream, "%s", "<output_generation>\n")
	
	close(file_stream)

	nothing
end