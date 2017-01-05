function write_xml_key(file_stream::IOStream, key_name::String, key_value)
	
	@printf(file_stream, "%s", "\u0009")
	@printf(file_stream, "%s", join(("<", key_name, ">")))
	@printf(file_stream, "%0.3f", key_value)
	@printf(file_stream, "%s", join(("</", key_name, ">")))
	@printf(file_stream, "%s", "\n")
	
	nothing
end