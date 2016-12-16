function read_xml_input(input_file_path::String)
	file_stream::IOStream = open(input_file_path, "r")
	file_string::String = readstring(file_stream)
	close(file_stream)
	
	# Read path to particle system configuration.
	ind_before = search(file_string, "<particle_system_path>")
	ind_after = search(file_string, "</particle_system_path>")
	ind = ind_before[end]+1:ind_after[1]-1
	particle_system_path::String = file_string[ind]
	
	# Read path to output file.
	ind_before = search(file_string, "<output_path>")
	ind_after = search(file_string, "</output_path>")
	ind = ind_before[end]+1:ind_after[1]-1
	output_path::String = file_string[ind]
	
	# Read inherent diffusion coefficient.
	ind_before = search(file_string, "<inherent_diffusion_coefficient>")
	ind_after = search(file_string, "</inherent_diffusion_coefficient>")
	ind = ind_before[end]+1:ind_after[1]-1
	D0::Float64 = parse(Float64, file_string[ind])
	
	# Read deltat_coarse.
	ind_before = search(file_string, "<deltat_coarse>")
	ind_after = search(file_string, "</deltat_coarse>")
	ind = ind_before[end]+1:ind_after[1]-1
	deltat_coarse::Float64 = parse(Float64, file_string[ind])
	
	# Read number_of_time_points_coarse.
	ind_before = search(file_string, "<number_of_time_points_coarse>")
	ind_after = search(file_string, "</number_of_time_points_coarse>")
	ind = ind_before[end]+1:ind_after[1]-1
	number_of_time_points_coarse::Int64 = parse(Int64, file_string[ind])
	
	# Read number_of_time_points_fine_per_coarse.
	ind_before = search(file_string, "<number_of_time_points_fine_per_coarse>")
	ind_after = search(file_string, "</number_of_time_points_fine_per_coarse>")
	ind = ind_before[end]+1:ind_after[1]-1
	number_of_time_points_fine_per_coarse::Int64 = parse(Int64, file_string[ind])
	
	# Read number_of_diffusers.
	ind_before = search(file_string, "<number_of_diffusers>")
	ind_after = search(file_string, "</number_of_diffusers>")
	ind = ind_before[end]+1:ind_after[1]-1
	number_of_diffusers::Int64 = parse(Int64, file_string[ind])
	
	# Read number_of_cells_x.
	ind_before = search(file_string, "<number_of_cells_x>")
	ind_after = search(file_string, "</number_of_cells_x>")
	ind = ind_before[end]+1:ind_after[1]-1
	number_of_cells_x::Int64 = parse(Int64, file_string[ind])
	
	# Read number_of_cells_y.
	ind_before = search(file_string, "<number_of_cells_y>")
	ind_after = search(file_string, "</number_of_cells_y>")
	ind = ind_before[end]+1:ind_after[1]-1
	number_of_cells_y::Int64 = parse(Int64, file_string[ind])
	
	# Read number_of_cells_z.
	ind_before = search(file_string, "<number_of_cells_z>")
	ind_after = search(file_string, "</number_of_cells_z>")
	ind = ind_before[end]+1:ind_after[1]-1
	number_of_cells_z::Int64 = parse(Int64, file_string[ind])
	
	return (particle_system_path, output_path, D0, deltat_coarse, number_of_time_points_coarse, number_of_time_points_fine_per_coarse, number_of_diffusers, number_of_cells_x, number_of_cells_y, number_of_cells_z)
end

#read_xml_input("../input.xml")