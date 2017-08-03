function read_output(file_path::String)
	file_stream::IOStream = open(file_path, "r")
	file_string::String = readstring(file_stream)
	close(file_stream)

	# Particle data.
	particle_type::String = read_key(file_string, "particle_type", String)
	R_temp::Array{Float64, 1} = read_key(file_string, "R", Array{Float64, 1})
	number_of_properties::Int64 = 0
	if particle_type == "sphere"
		number_of_properties = 1
	elseif particle_type == "ellipse"
		number_of_properties = 2
	elseif particle_type == "ellipsoid"
		number_of_properties = 3
	elseif particle_type == "cuboid"
		number_of_properties = 3
	end
	number_of_particles::Int64 = length(R_temp) / number_of_properties
	R::Array{Float64, 2} = reshape(R_temp, number_of_particles, number_of_properties)

	# System size.
	Lx::Float64 = read_key(file_string, "domain_size_x", Float64)
	Ly::Float64 = read_key(file_string, "domain_size_y", Float64)
	Lz::Float64 = read_key(file_string, "domain_size_z", Float64)
	phi::Float64 = read_key(file_string, "phi", Float64)

	# Positions.
	X::Array{Float64, 1} = read_key(file_string, "X", Array{Float64, 1})
	Y::Array{Float64, 1} = read_key(file_string, "Y", Array{Float64, 1})
	Z::Array{Float64, 1} = read_key(file_string, "Z", Array{Float64, 1})

	# Orientations.
	Q0::Array{Float64, 1} = read_key(file_string, "Q0", Array{Float64, 1})
	Q1::Array{Float64, 1} = read_key(file_string, "Q1", Array{Float64, 1})
	Q2::Array{Float64, 1} = read_key(file_string, "Q2", Array{Float64, 1})
	Q3::Array{Float64, 1} = read_key(file_string, "Q3", Array{Float64, 1})

	# Other.
	execution_time::Float64 = read_key(file_string, "execution_time", Float64)

	return (
		particle_type,
		R,
		Lx,
		Ly,
		Lz,
		phi,
		X,
		Y,
		Z,
		Q0,
		Q1,
		Q2,
		Q3,
		execution_time)
end
