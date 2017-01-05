function get_xml_key(file_string::String, key_name::String, fmt::DataType)
	
	ind_before::UnitRange{Int64} = search(file_string, join(("<", key_name, ">")))
	ind_after::UnitRange{Int64} = search(file_string, join(("</", key_name, ">")))
	ind::UnitRange{Int64} = ind_before[end]+1:ind_after[1]-1
	
	key_string::String = file_string[ind]
	
	if fmt == String
		return key_string
	else
		return parse(fmt, key_string)
	end
end