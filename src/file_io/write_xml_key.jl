function write_xml_key(file_stream::IOStream, key_name::String, key_value)
	
	@printf(file_stream, "%s", "\u0009")
	@printf(file_stream, "%s", join(("<", key_name, ">")))
	
	if typeof(key_value) == Float64
		@printf(file_stream, "%0.6f", key_value)
	elseif typeof(key_value) == Array{Float64, 1}
		for i = 1:length(key_value)-1
			@printf(file_stream, "%0.6f", key_value[i])
			@printf(file_stream, "%s", ",")
		end
		@printf(file_stream, "%0.6f", key_value[end])
	elseif typeof(key_value) == Int64
		@printf(file_stream, "%d", key_value)
	elseif typeof(key_value) == String
		@printf(file_stream, "%s", key_value)
	else
		error("Incompatible data type.")
	end
	
	@printf(file_stream, "%s", join(("</", key_name, ">")))
	@printf(file_stream, "%s", "\n")
	
	nothing
end