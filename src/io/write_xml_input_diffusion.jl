function write_xml_input_diffusion(file_path::String, output_generation_path::String, inherent_diffusion_coefficient::Float64, deltat_coarse::Float64, number_of_time_points_coarse::Int64, number_of_time_points_fine_per_coarse::Int64, number_of_diffusers::Int64, number_of_cells_x::Int64, number_of_cells_y::Int64, number_of_cells_z::Int64, output_diffusion_path::String)
	file_stream::IOStream = open(file_path, "w")
	
	@printf(file_stream, "%s", "<input_diffusion>\n")

	write_xml_key(file_stream, "output_generation_path", output_generation_path)
	write_xml_key(file_stream, "inherent_diffusion_coefficient", inherent_diffusion_coefficient)
	write_xml_key(file_stream, "deltat_coarse", deltat_coarse)
	write_xml_key(file_stream, "number_of_time_points_coarse", number_of_time_points_coarse)
	write_xml_key(file_stream, "number_of_time_points_fine_per_coarse", number_of_time_points_fine_per_coarse)
	write_xml_key(file_stream, "number_of_diffusers", number_of_diffusers)
	write_xml_key(file_stream, "number_of_cells_x", number_of_cells_x)
	write_xml_key(file_stream, "number_of_cells_y", number_of_cells_y)
	write_xml_key(file_stream, "number_of_cells_z", number_of_cells_z)
	write_xml_key(file_stream, "output_diffusion_path", output_diffusion_path)
	
	@printf(file_stream, "%s", "</input_diffusion>")
	
	close(file_stream)

	nothing
end