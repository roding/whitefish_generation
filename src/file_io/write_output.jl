function write_output(	file_path::String, 
					particle_type::String,
					R::Array{Float64, 2},
					Lx::Float64, 
					Ly::Float64, 
					Lz::Float64,
					phi::Float64,
					X::Array{Float64, 1}, 
					Y::Array{Float64, 1}, 
					Z::Array{Float64, 1}, 
					Q0::Array{Float64, 1}, 
					Q1::Array{Float64, 1}, 
					Q2::Array{Float64, 1}, 
					Q3::Array{Float64, 1}, 
					t_exec::Float64)

	file_stream::IOStream = open(file_path, "w")
	
	@printf(file_stream, "%s", "<output_generation>\n")

	write_key(file_stream, "particle_type", particle_type)
	write_key(file_stream, "R", R[:])
	write_key(file_stream, "domain_size_x", Lx)
	write_key(file_stream, "domain_size_y", Ly)
	write_key(file_stream, "domain_size_z", Lz)
	write_key(file_stream, "phi", phi)
	write_key(file_stream, "X", X)
	write_key(file_stream, "Y", Y)
	write_key(file_stream, "Z", Z)
	write_key(file_stream, "Q0", Q0)
	write_key(file_stream, "Q1", Q1)
	write_key(file_stream, "Q2", Q2)
	write_key(file_stream, "Q3", Q3)
	write_key(file_stream, "execution_time", t_exec)
	
	@printf(file_stream, "%s", "</output_generation>")
	
	close(file_stream)

	nothing
end