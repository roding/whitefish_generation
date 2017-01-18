include("../src/io/write_xml_key.jl")
include("../src/io/write_xml_input_generation.jl")
include("../src/io/write_xml_input_diffusion.jl")

function run_generation_and_diffusion()
	# Inititalization of random number generation device.
	const random_seed::Int64 = convert(Int64, time_ns())
	srand(random_seed)
	
	# Change current folder to the folder where this script lies.
	(program_file_dir::String, program_file_name::String) = splitdir(PROGRAM_FILE)
	program_file_dir = abspath(program_file_dir)
	cd(program_file_dir)
	
	# Main output directory.
	main_output_dir::String = abspath(joinpath(program_file_dir, "../../output_layer"))
	if !isdir(main_output_dir)
        mkdir(main_output_dir)
    end
	
	# Output directory.
	output_dir::String = abspath(joinpath(main_output_dir, join(("output_", string(random_seed)))))
	if !isdir(output_dir)
		mkdir(output_dir)
	end
	
	# Generation parameters.
	Lx::Float64 = 100.0 # µm.
	Ly::Float64 = 100.0 # µm.
	Lz::Float64 = 200.0 # µm.
	particle_type::String = "ellipticaldisk"
	number_of_particles::Int64 = 100
	thickness::Float64 = 50.0 + (Lz - 50.0) * rand() # µm.
	lbz::Float64 = 0.5 * (Lz - thickness) # µm.
	ubz::Float64 = 0.5 * (Lz + thickness) # µm.
	ubangle::Float64 = pi # Radians.
	R1::Array{Float64,1} = 10.0 * ones(number_of_particles) # µm.
	R2::Array{Float64,1} = 10.0 * ones(number_of_particles) # µm.
	number_of_equilibration_sweeps::Int64 = 1000
	output_generation_path::String = join((output_dir, "/", "output_generation.xml"))
	
	# Write input file for generation.
	input_generation_path::String = join((output_dir, "/", "input_generation.xml"))
	write_xml_input_generation(input_generation_path, Lx, Ly, Lz, particle_type, number_of_particles, lbz, ubz, ubangle, R1, R2,number_of_equilibration_sweeps, output_generation_path)

	# Run generation.
	program_generation_path::String = "../src/generation/wfrun_generation.jl"
	cmd_generation::Cmd = `julia $program_generation_path $input_generation_path`
	run(cmd_generation)
	
	# Diffusion parameters.
	D0::Float64 = 1.0
	deltat_coarse::Float64 = 5.0
	number_of_time_points_coarse::Int64 = 5000
	number_of_time_points_fine_per_coarse::Int64 = 100
	number_of_diffusers::Int64 = 500000
	number_of_cells_x::Int64 = 10
	number_of_cells_y::Int64 = 10
	number_of_cells_z::Int64 = 20
	output_diffusion_path::String = join((output_dir, "/", "output_diffusion.xml"))
	
	# Write input file for diffusion.
	input_diffusion_path::String = join((output_dir, "/", "input_diffusion.xml"))
	write_xml_input_diffusion(input_diffusion_path, output_generation_path, D0, deltat_coarse, number_of_time_points_coarse, number_of_time_points_fine_per_coarse, number_of_diffusers, number_of_cells_x, number_of_cells_y, number_of_cells_z, output_diffusion_path)
	
	# Run diffusion.
	program_diffusion_path::String = "../src/diffusion/wfrun_diffusion.jl"
	cmd_diffusion::Cmd = `julia $program_diffusion_path $input_diffusion_path`
	run(cmd_diffusion)
	
	# Exit.
	nothing
	
end

while true
	run_generation_and_diffusion()
end
