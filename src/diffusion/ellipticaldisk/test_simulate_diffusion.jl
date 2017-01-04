include("simulate_diffusion.jl")

function test_simulate_diffusion()

	# Initialization.
	const t_start::Int64 = time_ns()
	const seed::Int64 = t_start
	srand(seed)
	println(seed)
	
	# Simulation domain parameters.
	Lx::Float64 = 100.0
	Ly::Float64 = 100.0
	Lz::Float64 = 200.0
	
	# Particle parameters (lengths are in µm).
	number_of_particles::Int64 = 500
	X::Array{Float64,1} = Lx * rand(number_of_particles)
	Y::Array{Float64,1} = Ly * rand(number_of_particles)
	Z::Array{Float64,1} = 100.0 * zeros(number_of_particles)#Lz * rand(number_of_particles)
	THETA1::Array{Float64,1} = zeros(number_of_particles)#randn(number_of_particles)
	THETA2::Array{Float64,1} = zeros(number_of_particles)#randn(number_of_particles)
	THETA3::Array{Float64,1} = zeros(number_of_particles)#randn(number_of_particles)
	R1::Array{Float64,1} = 10.0 * ones(number_of_particles)
	R2::Array{Float64,1} = 10.0 * ones(number_of_particles)
	
		# Simulation parameters.
	const D0::Float64 = 1.0 # Intrinsic diffusion coefficient.	
	const deltat_coarse::Float64 = 10.0 # Time step.
	const number_of_time_points_coarse::Int64 = 2500
	const number_of_time_points_fine_per_coarse::Int64 = 10
	const number_of_diffusers::Int64 = 10000 # Number of diffusing particles.
	
	number_of_cells_x::Int64 = 10
	number_of_cells_y::Int64 = 10
	number_of_cells_z::Int64 = 20
	
	silent_mode::Bool = false

	# Run simulation.
	println("   Simulating diffusion...")
	t_start_diffusion::Int64 = convert(Int64, time_ns())
	(msd_x, msd_y, msd_z, D0_empirical) = simulate_diffusion(	X,
															Y,
															Z,
															THETA1,
															THETA2,
															THETA3,
															R1,
															R2,
															Lx,
															Ly,
															Lz,
															D0,
															deltat_coarse,
															number_of_time_points_coarse,
															number_of_time_points_fine_per_coarse,
															number_of_diffusers,
															number_of_cells_x,
															number_of_cells_y,
															number_of_cells_z,
															silent_mode)
																				
	const t::Array{Float64, 1} = 0.0:deltat_coarse:(convert(Float64, number_of_time_points_coarse-1) * deltat_coarse)
	#const deltat_fine::Float64 = deltat_coarse / convert(Float64, number_of_time_points_fine_per_coarse)
	#const sigma::Float64 = sqrt(2 * D0 * deltat_fine)
	#println(sigma_empirical/sigma)
	
	t_exec_diffusion::Float64 = convert(Float64, time_ns() - t_start_diffusion) / 1e9
	file_name_diffusion = "diffusion.dat"
	file_stream_diffusion = open(file_name_diffusion,"w")
	writedlm(file_stream_diffusion, t_exec_diffusion, ",")
	writedlm(file_stream_diffusion, D0_empirical/D0, ",")
	writedlm(file_stream_diffusion, [t msd_x msd_y msd_z], ",")
	close(file_stream_diffusion)
	
	println(join(("      Finished in ", string(t_exec_diffusion), " seconds.")))
	
	nothing
	
end

test_simulate_diffusion()
