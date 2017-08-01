function [  particle_type, ...
            R, ...
            Lx, ...
            Ly, ...
            Lz, ...
            phi, ...
            X, ...
            Y, ...
            Z, ...
            Q0, ...
            Q1, ...
            Q2, ...
            Q3, ...
            execution_time] = read_output(file_path)

file_string = fileread(file_path);

particle_type = read_key(file_string, 'particle_type', 'string');

R = read_key(file_string, 'R', 'array');
switch particle_type
    case "sphere"
        number_of_properties = 1;
    case "ellipse"
        number_of_properties = 2;
    case "ellipsoid"
        number_of_properties = 3;
    case "cuboid"
        number_of_properties = 3;
end
number_of_particles = numel(R) / number_of_properties;
R = reshape(R, [number_of_particles, number_of_properties]);

Lx = read_key(file_string, 'domain_size_x', 'scalar');
Ly = read_key(file_string, 'domain_size_y', 'scalar');
Lz = read_key(file_string, 'domain_size_z', 'scalar');
phi = read_key(file_string, 'phi', 'scalar');
X = read_key(file_string, 'X', 'array');
Y = read_key(file_string, 'Y', 'array');
Z = read_key(file_string, 'Z', 'array');
Q0 = read_key(file_string, 'Q0', 'array');
Q1 = read_key(file_string, 'Q1', 'array');
Q2 = read_key(file_string, 'Q2', 'array');
Q3 = read_key(file_string, 'Q3', 'array');
execution_time = read_key(file_string, 'execution_time', 'scalar');

end

