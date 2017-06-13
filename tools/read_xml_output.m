function [Lx, Ly, Lz, particle_type, number_of_particles, X, Y, Z, THETA1, THETA2, THETA3, R1, R2, t_exec_generation] = read_xml_output_generation(file_path)

file_string = fileread(file_path);

Lx = read_xml_key(file_string, 'domain_size_x', 'scalar');
Ly = read_xml_key(file_string, 'domain_size_y', 'scalar');
Lz = read_xml_key(file_string, 'domain_size_z', 'scalar');
particle_type = read_xml_key(file_string, 'particle_type', 'string');
number_of_particles = read_xml_key(file_string, 'number_of_particles', 'scalar');
X = read_xml_key(file_string, 'X', 'array');
Y = read_xml_key(file_string, 'Y', 'array');
Z = read_xml_key(file_string, 'Z', 'array');
THETA1 = read_xml_key(file_string, 'THETA1', 'array');
THETA2 = read_xml_key(file_string, 'THETA2', 'array');
THETA3 = read_xml_key(file_string, 'THETA3', 'array');
R1 = read_xml_key(file_string, 'R1', 'array');
R2 = read_xml_key(file_string, 'R2', 'array');
t_exec_generation = read_xml_key(file_string, 'execution_time', 'scalar');

end

