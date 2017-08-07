function write_input(	file_path::String,
					particle_type::String,
					R::Array{Float64, 2},
					Lx::Float64,
					Ly::Float64,
					Lz::Float64,
					phi::Float64,
					phi_target::Float64,
					position_constraint_axis::String,
					position_constraint_lower::Float64,
					position_constraint_upper::Float64,
					orientation_axis::Array{Float64, 1},
					orientation_constraint_axis::Array{Float64, 1},
					orientation_constraint_lower::Float64,
					orientation_constraint_upper::Float64,
					sigma_translation_max::Float64,
					sigma_rotation_max::Float64,
					sigma_ratio::Float64,
					number_of_equilibration_sweeps::Int64,
					delta_phi::Float64,
					output_file_path::String)

	file_stream::IOStream = open(file_path, "w")

	@printf(file_stream, "%s", "<input_generation>\n")

	write_key(file_stream, "particle_type", particle_type)
	write_key(file_stream, "R", R[:])
	write_key(file_stream, "domain_size_x", Lx)
	write_key(file_stream, "domain_size_y", Ly)
	write_key(file_stream, "domain_size_z", Lz)
	write_key(file_stream, "phi", phi)
	write_key(file_stream, "phi_target", phi_target)
	write_key(file_stream, "position_constraint_axis", position_constraint_axis)
	write_key(file_stream, "position_constraint_lower", position_constraint_lower)
	write_key(file_stream, "position_constraint_upper", position_constraint_upper)
	write_key(file_stream, "orientation_axis", orientation_axis)
	write_key(file_stream, "orientation_constraint_axis", orientation_constraint_axis)
	write_key(file_stream, "orientation_constraint_lower", orientation_constraint_lower)
	write_key(file_stream, "orientation_constraint_upper", orientation_constraint_upper)
	write_key(file_stream, "sigma_translation_max", sigma_translation_max)
	write_key(file_stream, "sigma_rotation_max", sigma_rotation_max)
	write_key(file_stream, "sigma_ratio", sigma_ratio)
	write_key(file_stream, "number_of_equilibration_sweeps", number_of_equilibration_sweeps)
	write_key(file_stream, "delta_phi", delta_phi)
	write_key(file_stream, "output_file_path", output_file_path)

	@printf(file_stream, "%s", "</input_generation>")

	close(file_stream)

	nothing
end
