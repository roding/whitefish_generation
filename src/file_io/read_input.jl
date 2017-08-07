function read_input(file_path::String)
	file_stream::IOStream = open(file_path, "r")
	file_string::String = readstring(file_stream)
	close(file_stream)

	# Read particle data.
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

	# Read input on initial system size.
	Lx::Float64 = read_key(file_string, "domain_size_x", Float64)
	Ly::Float64 = read_key(file_string, "domain_size_y", Float64)
	Lz::Float64 = read_key(file_string, "domain_size_z", Float64)
	phi::Float64 = read_key(file_string, "phi", Float64)

	# Read input on target system size.
	phi_target::Float64 = read_key(file_string, "phi_target", Float64)

	# Position constraint.
	position_constraint_axis::String = read_key(file_string, "position_constraint_axis", String)
	position_constraint_lower::Float64 = read_key(file_string, "position_constraint_lower", Float64)
	position_constraint_upper::Float64 = read_key(file_string, "position_constraint_upper", Float64)

	# Orientation constraint.
	orientation_axis::Array{Float64, 1} = read_key(file_string, "orientation_axis", Array{Float64, 1})
	orientation_constraint_axis::Array{Float64, 1} = read_key(file_string, "orientation_constraint_axis", Array{Float64, 1})
	orientation_constraint_lower::Float64 = read_key(file_string, "orientation_constraint_lower", Float64)
	orientation_constraint_upper::Float64 = read_key(file_string, "orientation_constraint_upper", Float64)

	# Other input.
	sigma_translation_max::Float64 = read_key(file_string, "sigma_translation_max", Float64)
	sigma_rotation_max::Float64 = read_key(file_string, "sigma_rotation_max", Float64)
	sigma_ratio::Float64 = read_key(file_string, "sigma_ratio", Float64)
	number_of_equilibration_sweeps::Int64 = read_key(file_string, "number_of_equilibration_sweeps", Int64)
	delta_phi::Float64 = read_key(file_string, "delta_phi", Float64)

	output_file_path::String = read_key(file_string, "output_file_path", String)

	return (
		particle_type,
		R,
		Lx,
		Ly,
		Lz,
		phi,
		phi_target,
		position_constraint_axis,
		position_constraint_lower,
		position_constraint_upper,
		orientation_axis,
		orientation_constraint_axis,
		orientation_constraint_lower,
		orientation_constraint_upper,
		sigma_translation_max,
		sigma_rotation_max,
		sigma_ratio,
		number_of_equilibration_sweeps,
		delta_phi,
		output_file_path)
end
