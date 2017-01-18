function [diagnostic_diffusion_coefficient_ratio, t, msd_x, msd_y, msd_z, t_exec_diffusion] = read_xml_output_diffusion(file_path)

file_string = fileread(file_path);

diagnostic_diffusion_coefficient_ratio = read_xml_key(file_string, 'diagnostic_diffusion_coefficient_ratio', 'scalar');
t = read_xml_key(file_string, 'time', 'array');
msd_x = read_xml_key(file_string, 'mean_square_displacement_x', 'array');
msd_y = read_xml_key(file_string, 'mean_square_displacement_y', 'array');
msd_z = read_xml_key(file_string, 'mean_square_displacement_z', 'array');
t_exec_diffusion = read_xml_key(file_string, 'execution_time', 'scalar');

end

