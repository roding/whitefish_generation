function read_key(file_string::String, key_name::String, fmt::DataType)
	
	ind_before::UnitRange{Int64} = search(file_string, join(("<", key_name, ">")))
	ind_after::UnitRange{Int64} = search(file_string, join(("</", key_name, ">")))
	ind::UnitRange{Int64} = ind_before[end]+1:ind_after[1]-1
	
	key_string::String = file_string[ind]
	
	if fmt == String
		return key_string
	elseif fmt == Int64 || fmt == Float64
		return parse(fmt, key_string)
	elseif fmt == Array{Float64, 1}
		key_string_array = split(key_string, ",")
		number_of_values = length(key_string_array)
		data = Array{Float64}(number_of_values)
		for current_value = 1:number_of_values
			data[current_value] = parse(Float64, key_string_array[current_value])
		end
		return data
	else
		error("Incompatible data type.")
		nothing
	end
end