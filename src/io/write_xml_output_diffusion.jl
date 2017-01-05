function write_xml_output(file_name::String, D0::Float64, D0_empirical::Float64, deltat_coarse::Float64, number_of_time_points_coarse::Int64, msd_x::Array{Float64, 1}, msd_y::Array{Float64, 1}, msd_z::Array{Float64, 1})
	
	file_stream::IOStream = open(file_name, "w")
		
	t::Array{Float64, 1} = 0.0:deltat_coarse:(convert(Float64, number_of_time_points_coarse-1) * deltat_coarse)
	
	@printf(file_stream, "%s", "<output>\n")
	
	@printf(file_stream, "%s", "\u0009")
	@printf(file_stream, "%s", "<diagnostic_diffusion_coefficient_ratio>")
	@printf(file_stream, "%0.3f", D0_empirical/D0) 
	@printf(file_stream, "%s", "</diagnostic_diffusion_coefficient_ratio>\n")
	
	@printf(file_stream, "%s", "\u0009")
	@printf(file_stream, "%s", "<time>")
	for i = 1:length(t)-1
		@printf(file_stream, "%0.3f", t[i])
		@printf(file_stream, "%s", ",")
	end
	@printf(file_stream, "%0.3f", t[end])
	@printf(file_stream, "%s", "</time>\n")
	
	@printf(file_stream, "%s", "\u0009")
	@printf(file_stream, "%s", "<mean_square_displacement_x>")
	for i = 1:length(t)-1
		@printf(file_stream, "%0.3f", msd_x[i])
		@printf(file_stream, "%s", ",")
	end
	@printf(file_stream, "%0.3f", msd_x[end])
	@printf(file_stream, "%s", "</mean_square_displacement_x>\n")
	
	@printf(file_stream, "%s", "\u0009")
	@printf(file_stream, "%s", "<mean_square_displacement_y>")
	for i = 1:length(t)-1
		@printf(file_stream, "%0.3f", msd_y[i])
		@printf(file_stream, "%s", ",")
	end
	@printf(file_stream, "%0.3f", msd_y[end])
	@printf(file_stream, "%s", "</mean_square_displacement_y>\n")
	
	@printf(file_stream, "%s", "\u0009")
	@printf(file_stream, "%s", "<mean_square_displacement_z>")
	for i = 1:length(t)-1
		@printf(file_stream, "%0.3f", msd_z[i])
		@printf(file_stream, "%s", ",")
	end
	@printf(file_stream, "%0.3f", msd_z[end])
	@printf(file_stream, "%s", "</mean_square_displacement_z>\n")
	
	@printf(file_stream, "%s", "</output>\n")
	
	close(file_stream)
	
	nothing
end

#read_xml_system("particle_system.xml")