function print_simulation_stats(number_of_particles::Int64, number_of_diffusers::Int64, Lx::Float64, Ly::Float64, Lz::Float64)
	println("Simulation stats:")
	println(join(("   Simulation domain size:   x: ", string(Lx))))
	println(join(("                             y: ", string(Ly))))
	println(join(("                             z: ", string(Lz))))
	println(join(("   Number of particles:         ", string(number_of_particles))))
	println(join(("   Number of diffusers:         ", string(number_of_diffusers))))
	
	nothing
end
