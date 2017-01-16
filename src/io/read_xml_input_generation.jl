include("read_xml_key.jl")

function read_xml_input_generation(file_path::String)
	file_stream::IOStream = open(file_path, "r")
	file_string::String = readstring(file_stream)
	close(file_stream)
			
	Lx::Float64 = read_xml_key(file_string, "domain_size_x", Float64)
	Ly::Float64 = read_xml_key(file_string, "domain_size_y", Float64)
	Lz::Float64 = read_xml_key(file_string, "domain_size_z", Float64)
	particle_type::String = read_xml_key(file_string, "particle_type", String)
	number_of_particles::Int64 = read_xml_key(file_string, "number_of_particles", Int64)
	lbz::Float64 = read_xml_key(file_string, "lower_bound_z", Float64)
	ubz::Float64 = read_xml_key(file_string, "upper_bound_z", Float64)
	ubangle::Float64 = read_xml_key(file_string, "upper_bound_angle_to_z_axis", Float64)
	R1::Array{Float64,1} = read_xml_key(file_string, "R1", Array{Float64, 1})
	R2::Array{Float64,1} = read_xml_key(file_string, "R2", Array{Float64, 1})
	number_of_equilibration_sweeps::Int64 = read_xml_key(file_string, "number_of_equilibration_sweeps", Int64)
	output_file_path::String = read_xml_key(file_string, "output_file_path", String)
	
	return (Lx, Ly, Lz, particle_type, number_of_particles, lbz, ubz, ubangle, R1, R2, number_of_equilibration_sweeps, output_file_path)
end