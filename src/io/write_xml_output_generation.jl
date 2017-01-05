include("write_xml_key.jl")

function write_xml_output_generation(output_file_path::String)
	file_stream::IOStream = open(output_file_path, "w")
	
	@printf(file_stream, "%s", "<output_generation>\n")

	write_xml_key(file_stream, "domain_size_x", Lx)
	write_xml_key(file_stream, "domain_size_y", Ly)
	write_xml_key(file_stream, "domain_size_z", Lz)
	write_xml_key(file_stream, "particle_type", "ellipticaldisk")
	write_xml_key(file_stream, "number_of_particles", number_of_particles)
	write_xml_key(file_stream, "X", X)
	
	@printf(file_stream, "%s", "</output_generation>")
	
	close(file_stream)

	nothing
end

#write_xml_output_generation("test_output.xml")

<particle_system>
	<domain_size_x>100.0</domain_size_x>
	<domain_size_y>100.0</domain_size_y>
	<domain_size_z>100.0</domain_size_z>
	<type>ellipticaldisk</type>
	<number_of_particles>1</number_of_particles>
	<X>10.0</X>
	<Y>10.0</Y>
	<Z>10.0</Z>
	<THETA1>0.0</THETA1>
	<THETA2>1.0</THETA2>
	<THETA3>2.0</THETA3>
	<R1>10.0</R1>
	<R2>10.0</R2>
</particle_system>