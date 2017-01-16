function read_xml_input_diffusion(file_path::String)
	file_stream::IOStream = open(file_path, "r")
	file_string::String = readstring(file_stream)
	close(file_stream)
	
	output_generation_path::String = read_xml_key(file_string, "output_generation_path", String)
	inherent_diffusion_coefficient::Float64 = read_xml_key(file_string, "inherent_diffusion_coefficient", Float64)
	deltat_coarse::Float64 = read_xml_key(file_string, "deltat_coarse", Float64)
	number_of_time_points_coarse::Int64 = read_xml_key(file_string, "number_of_time_points_coarse", Int64)
	number_of_time_points_fine_per_coarse::Int64 = read_xml_key(file_string, "number_of_time_points_fine_per_coarse", Int64)
	number_of_diffusers::Int64 = read_xml_key(file_string, "number_of_diffusers", Int64)
	number_of_cells_x::Int64 = read_xml_key(file_string, "number_of_cells_x", Int64)
	number_of_cells_y::Int64 = read_xml_key(file_string, "number_of_cells_y", Int64)
	number_of_cells_z::Int64 = read_xml_key(file_string, "number_of_cells_z", Int64)
	output_diffusion_path::String = read_xml_key(file_string, "output_diffusion_path", String)
	
	return (output_generation_path, inherent_diffusion_coefficient, deltat_coarse, number_of_time_points_coarse, number_of_time_points_fine_per_coarse, number_of_diffusers, number_of_cells_x, number_of_cells_y, number_of_cells_z, output_diffusion_path)
end

#read_xml_input_diffusion("../io_test_files/input_diffusion.xml")