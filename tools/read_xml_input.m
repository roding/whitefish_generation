function [Lx, Ly, Lz, particle_type, number_of_particles, lbz, ubz, ubangle, R1, R2, number_of_equilibration_sweeps, output_generation_path] = read_xml_input_generation(file_path)

file_string = fileread(file_path);

Lx = read_xml_key(file_string, 'domain_size_x', 'scalar');
Ly = read_xml_key(file_string, 'domain_size_y', 'scalar');
Lz = read_xml_key(file_string, 'domain_size_z', 'scalar');
particle_type = read_xml_key(file_string, 'particle_type', 'string');
number_of_particles = read_xml_key(file_string, 'number_of_particles', 'scalar');
lbz = read_xml_key(file_string, 'lower_bound_z', 'scalar');
ubz = read_xml_key(file_string, 'upper_bound_z', 'scalar');
ubangle = read_xml_key(file_string, 'upper_bound_angle_to_z_axis', 'scalar');
R1 = read_xml_key(file_string, 'R1', 'array');
R2 = read_xml_key(file_string, 'R2', 'array');
number_of_equilibration_sweeps = read_xml_key(file_string, 'number_of_equilibration_sweeps', 'scalar');
output_generation_path = read_xml_key(file_string, 'output_file_path', 'string');

end

