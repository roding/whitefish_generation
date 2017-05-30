function print_simulation_stats(particle_type::String, number_of_particles::Int64, Lx::Float64, Ly::Float64, Lz::Float64)
	println("Simulation stats:")
	println(join(("   Particle type:               ", particle_type)))
	println(join(("   Simulation domain size:   x: ", string(Lx))))
	println(join(("                             y: ", string(Ly))))
	println(join(("                             z: ", string(Lz))))
	println(join(("   Number of particles:         ", string(number_of_particles))))
	
	nothing
end
