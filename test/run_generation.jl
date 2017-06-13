include("../src/file_io/write_xml_key.jl")
include("../src/file_io/write_xml_input.jl")

function run_generation()
	# Inititalization of random number generation device.
	random_seed::Int64 = convert(Int64, time_ns())
	srand(random_seed)

	# Change current folder to the folder where this script lies.
	(program_file_dir::String, program_file_name::String) = splitdir(PROGRAM_FILE)
	program_file_dir = abspath(program_file_dir)
	cd(program_file_dir)

	# Output directory.
	output_dir::String = abspath(joinpath(program_file_dir, "../../output"))
	if !isdir(output_dir)
        mkdir(output_dir)
    end

	# Generation parameters.
	Lx::Float64 = 100.0 # µm.
	Ly::Float64 = 100.0 # µm.
	Lz::Float64 = 100.0 # µm.
	particle_type::String = "ellipticaldisk"
	number_of_particles::Int64 = 50
	thickness::Float64 = 25.0#Lz#25.0 + (Lz - 25.0) * rand() # µm.
	lbz::Float64 = 0.5 * (Lz - thickness) # µm.
	ubz::Float64 = 0.5 * (Lz + thickness) # µm.
	ubangle::Float64 = pi # Radians.
	R1::Array{Float64,1} = 7.5 * ones(number_of_particles) # µm.
	R2::Array{Float64,1} = 7.5 * ones(number_of_particles) # µm.
	number_of_equilibration_sweeps::Int64 = 1#000

	# Write input file for generation.
	output_file_path::String = abspath(joinpath(output_dir, "output_generation.xml"))
	input_file_path::String = abspath(joinpath(output_dir, "input_generation.xml"))
	write_xml_input(input_file_path, Lx, Ly, Lz, particle_type, number_of_particles, lbz, ubz, ubangle, R1, R2, number_of_equilibration_sweeps, output_file_path)

	# Run generation.
	program_path::String = abspath("../src/wfrun_generation.jl")
	cmd::Cmd = `julia $program_path $input_file_path`
	run(cmd)

	# Exit.
	nothing

end

run_generation()
