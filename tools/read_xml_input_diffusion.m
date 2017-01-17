function [output_generation_path, inherent_diffusion_coefficient, deltat_coarse, number_of_time_points_coarse, number_of_time_points_fine_per_coarse, number_of_diffusers, number_of_cells_x, number_of_cells_y, number_of_cells_z, output_diffusion_path] = read_xml_input_diffusion(file_path)

file_string = fileread(file_path);

output_generation_path = read_xml_key(file_string, 'output_generation_path', 'string');
inherent_diffusion_coefficient = read_xml_key(file_string, 'inherent_diffusion_coefficient', 'scalar');
deltat_coarse = read_xml_key(file_string, 'deltat_coarse', 'scalar');
number_of_time_points_coarse = read_xml_key(file_string, 'number_of_time_points_coarse', 'scalar');
number_of_time_points_fine_per_coarse = read_xml_key(file_string, 'number_of_time_points_fine_per_coarse', 'scalar');
number_of_diffusers = read_xml_key(file_string, 'number_of_diffusers', 'scalar');
number_of_cells_x = read_xml_key(file_string, 'number_of_cells_x', 'scalar');
number_of_cells_y = read_xml_key(file_string, 'number_of_cells_y', 'scalar');
number_of_cells_z = read_xml_key(file_string, 'number_of_cells_z', 'scalar');
output_diffusion_path = read_xml_key(file_string, 'output_diffusion_path', 'string');

end

