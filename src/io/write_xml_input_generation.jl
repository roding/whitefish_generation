include("write_xml_key.jl")

function write_xml_input_generation(output_file_path::String, Lx::Float64, Ly::Float64, Lz::Float64, number_of_particles::Int64, lbz::Float64, ubz::Float64, ubangle::Float64, R1::Array{Float64, 1}, R2::Array{Float64, 1})
	file_stream::IOStream = open(output_file_path, "w")
	
	@printf(file_stream, "%s", "<input_generation>\n")

	write_xml_key(file_stream, "domain_size_x", Lx)
	write_xml_key(file_stream, "domain_size_y", Ly)
	write_xml_key(file_stream, "domain_size_z", Lz)
	write_xml_key(file_stream, "particle_type", "ellipticaldisk")
	write_xml_key(file_stream, "number_of_particles", number_of_particles)
	write_xml_key(file_stream, "lower_bound_z", lbz)
	write_xml_key(file_stream, "upper_bound_z", ubz)
	write_xml_key(file_stream, "upper_bound_angle_to_z_axis", ubangle)
	write_xml_key(file_stream, "R1", R1)
	write_xml_key(file_stream, "R2", R2)
	write_xml_key(file_stream, "number_of_equilibration_sweeps", number_of_equilibration_sweeps)
	write_xml_key(file_stream, "output_file_path", output_file_path)
	
	@printf(file_stream, "%s", "</input_generation>")
	
	close(file_stream)

	nothing
end