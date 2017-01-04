function write_xml_system(file_name::String, Lx::Float64, Ly::Float64, Lz::Float64, X::Array{Float64, 1}, Y::Array{Float64, 1}, Z::Array{Float64, 1}, THETA1::Array{Float64, 1}, THETA2::Array{Float64, 1}, THETA3::Array{Float64, 1}, R1::Array{Float64, 1}, R2::Array{Float64, 1})
	file_stream::IOStream = open(file_name, "w")
	
	@printf(file_stream, "%s", "<particle_system>\n")
	
	@printf(file_stream, "%s", "\u0009")
	@printf(file_stream, "%s", "<domain_size_x>")
	@printf(file_stream, "%0.6f", Lx) 
	@printf(file_stream, "%s", "</domain_size_x>\n")
	
	@printf(file_stream, "%s", "\u0009")
	@printf(file_stream, "%s", "<domain_size_y>")
	@printf(file_stream, "%0.6f", Ly) 
	@printf(file_stream, "%s", "</domain_size_y>\n")
	
	@printf(file_stream, "%s", "\u0009")
	@printf(file_stream, "%s", "<domain_size_z>")
	@printf(file_stream, "%0.6f", Lz) 
	@printf(file_stream, "%s", "</domain_size_z>\n")
	
	@printf(file_stream, "%s", "\u0009")
	@printf(file_stream, "%s", "<type>ellipticaldisk</type>\n")
	
	number_of_particles::Int64 = length(X)
	@printf(file_stream, "%s", "\u0009")
	@printf(file_stream, "%s", "<number_of_particles>")
	@printf(file_stream, "%d", number_of_particles)
	@printf(file_stream, "%s", "</number_of_particles>\n")
	
	@printf(file_stream, "%s", "\u0009")
	@printf(file_stream, "%s", "<X>")
	for i = 1:number_of_particles-1
		@printf(file_stream, "%0.6f", X[i])
		@printf(file_stream, "%s", ",")
	end
	@printf(file_stream, "%0.6f", X[end])
	@printf(file_stream, "%s", "</X>\n")
	
	@printf(file_stream, "%s", "\u0009")
	@printf(file_stream, "%s", "<Y>")
	for i = 1:number_of_particles-1
		@printf(file_stream, "%0.6f", Y[i])
		@printf(file_stream, "%s", ",")
	end
	@printf(file_stream, "%0.6f", Y[end])
	@printf(file_stream, "%s", "</Y>\n")
	
	@printf(file_stream, "%s", "\u0009")
	@printf(file_stream, "%s", "<Z>")
	for i = 1:number_of_particles-1
		@printf(file_stream, "%0.6f", Z[i])
		@printf(file_stream, "%s", ",")
	end
	@printf(file_stream, "%0.6f", Z[end])
	@printf(file_stream, "%s", "</Z>\n")
	
	@printf(file_stream, "%s", "\u0009")
	@printf(file_stream, "%s", "<THETA1>")
	for i = 1:number_of_particles-1
		@printf(file_stream, "%0.6f", THETA1[i])
		@printf(file_stream, "%s", ",")
	end
	@printf(file_stream, "%0.6f", THETA1[end])
	@printf(file_stream, "%s", "</THETA1>\n")
	
	@printf(file_stream, "%s", "\u0009")
	@printf(file_stream, "%s", "<THETA2>")
	for i = 1:number_of_particles-1
		@printf(file_stream, "%0.6f", THETA2[i])
		@printf(file_stream, "%s", ",")
	end
	@printf(file_stream, "%0.6f", THETA2[end])
	@printf(file_stream, "%s", "</THETA2>\n")
	
	@printf(file_stream, "%s", "\u0009")
	@printf(file_stream, "%s", "<THETA3>")
	for i = 1:number_of_particles-1
		@printf(file_stream, "%0.6f", THETA3[i])
		@printf(file_stream, "%s", ",")
	end
	@printf(file_stream, "%0.6f", THETA3[end])
	@printf(file_stream, "%s", "</THETA3>\n")
	
	@printf(file_stream, "%s", "\u0009")
	@printf(file_stream, "%s", "<R1>")
	for i = 1:number_of_particles-1
		@printf(file_stream, "%0.6f", R1[i])
		@printf(file_stream, "%s", ",")
	end
	@printf(file_stream, "%0.6f", R1[end])
	@printf(file_stream, "%s", "</R1>\n")
	
	@printf(file_stream, "%s", "\u0009")
	@printf(file_stream, "%s", "<R2>")
	for i = 1:number_of_particles-1
		@printf(file_stream, "%0.6f", R2[i])
		@printf(file_stream, "%s", ",")
	end
	@printf(file_stream, "%0.6f", R2[end])
	@printf(file_stream, "%s", "</R2>\n")
	
	@printf(file_stream, "%s", "</particle_system>\n")
	
	close(file_stream)
	
	nothing
end

#write_xml_system("test.xml", 10.0, 10.0, 10.0, [10.0,20.0], [10.0,20.0], [10.0,20.0], [10.0,20.0], [10.0,20.0], [10.0,20.0], [10.0,20.0], [10.0,20.0])