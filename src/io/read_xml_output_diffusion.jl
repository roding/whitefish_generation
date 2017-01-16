function read_xml_output_diffusion(file_path::String)
	file_stream::IOStream = open(file_path, "r")
	file_string::String = readstring(file_stream)
	close(file_stream)
			
	diagnostic_diffusion_coefficient_ratio::Float64 = read_xml_key(file_string, "diagnostic_diffusion_coefficient_ratio", Float64)
	t::Array{Float64, 1} = read_xml_key(file_string, "time", Array{Float64, 1})
	msd_x::Array{Float64, 1} = read_xml_key(file_string, "mean_square_displacement_x", Array{Float64, 1})
	msd_y::Array{Float64, 1} = read_xml_key(file_string, "mean_square_displacement_y", Array{Float64, 1})
	msd_z::Array{Float64, 1} = read_xml_key(file_string, "mean_square_displacement_z", Array{Float64, 1})
	t_exec::Float64 = read_xml_key(file_string, "execution_time", Float64)
	
	return (diagnostic_diffusion_coefficient_ratio, t, msd_x, msd_y, msd_z, t_exec)
end