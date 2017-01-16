function print_simulation_stats_diffusion(particle_type::String, number_of_particles::Int64, number_of_diffusers::Int64, Lx::Float64, Ly::Float64, Lz::Float64, numbers_of_cells_x::Int64, numbers_of_cells_y::Int64, numbers_of_cells_z::Int64)
	println("Simulation stats:")
	println(join(("   Particle type:               ", particle_type)))
	println(join(("   Simulation domain size:   x: ", string(Lx))))
	println(join(("                             y: ", string(Ly))))
	println(join(("                             z: ", string(Lz))))
	println(join(("   Number of particles:         ", string(number_of_particles))))
	println(join(("   Number of diffusers:         ", string(number_of_diffusers))))
	println(join(("   Number of cells:             ", string(numbers_of_cells_x), " x ", string(numbers_of_cells_y), " x ", string(numbers_of_cells_z))))
	
	nothing
end
