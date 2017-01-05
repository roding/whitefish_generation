include("write_xml_key.jl")

function write_xml_output_generation(output_file_path::String, Lx::Float64, Ly::Float64, Lz::Float64, number_of_particles::Int64, X::Array{Float64, 1}, Y::Array{Float64, 1}, Z::Array{Float64, 1}, THETA1::Array{Float64, 1}, THETA2::Array{Float64, 1}, THETA3::Array{Float64, 1}, R1::Array{Float64, 1}, R2::Array{Float64, 1})
	file_stream::IOStream = open(output_file_path, "w")
	
	@printf(file_stream, "%s", "<output_generation>\n")

	write_xml_key(file_stream, "domain_size_x", Lx)
	write_xml_key(file_stream, "domain_size_y", Ly)
	write_xml_key(file_stream, "domain_size_z", Lz)
	write_xml_key(file_stream, "particle_type", "ellipticaldisk")
	write_xml_key(file_stream, "number_of_particles", number_of_particles)
	write_xml_key(file_stream, "X", X)
	write_xml_key(file_stream, "Y", Y)
	write_xml_key(file_stream, "Z", Z)
	write_xml_key(file_stream, "THETA1", THETA1)
	write_xml_key(file_stream, "THETA2", THETA2)
	write_xml_key(file_stream, "THETA3", THETA3)
	write_xml_key(file_stream, "R1", R1)
	write_xml_key(file_stream, "R1", R2)
	
	@printf(file_stream, "%s", "</output_generation>")
	
	close(file_stream)

	nothing
end