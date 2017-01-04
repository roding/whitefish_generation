include("generate_elliptical_disk_system.jl")
include("write_xml_system.jl")

function run()

	# Initialization.
	const t_start::Int64 = time_ns()
	const random_seed::Int64 = t_start
	srand(random_seed)
	
	# Main output directory.
	main_output_dir = join((pwd(), "/", "output"))
	if !isdir(main_output_dir)
        mkdir(main_output_dir)
    end
	
	# Output directory.
	output_dir = join((main_output_dir, "/", "output_", string(random_seed)))
	if !isdir(output_dir)
		mkdir(output_dir)
	end
	
	# Simulation domain parameters (lengths are in µm).
	Lx::Float64 = 100.0
	Ly::Float64 = 100.0
	Lz::Float64 = 200.0
	
	# Particle parameters (lengths are in µm).
	number_of_particles::Int64 = 500#500
	R1::Array{Float64,1} = 10.0 * ones(number_of_particles)
	R2::Array{Float64,1} = 10.0 * ones(number_of_particles)
	
	# Bounds for particle positions.
	thickness::Float64 = Lz#10.0#50.0 + (Lz - 50.0) * rand()
	lbz::Float64 = 0.5 * (Lz - thickness)
	ubz::Float64 = 0.5 * (Lz + thickness)

	# Bound for angular deviation from z-axis.
	ubangle::Float64 = 0.0#pi#0.0 * pi
	
	# Compute weight fraction (physical).
	#w_bulk::Float64 = (Lx * Ly * Lz) * 1e-12 # µm^3 times 1e-12 g/µm^3 (the latter taken as a rough average of polyvinylalcohol and polyethene)
	#average_number_of_layers::Float64 = 10.0
	#w_graphene::Float64 = pi * sum(R1 .* R2) * 7.7e-16 * average_number_of_layers # 7.7e-16 g/µm^2 is the weight of graphene
	
	#psi::Float64 = w_graphene / w_bulk
	#println(psi)
	
	number_of_equilibration_sweeps::Int64 = 10#1000
	acceptance_probability_target::Float64 = 0.25
	number_of_iterations_overlap_criterion::Int64 = 20
	silent_mode::Bool = false
	
	println("   Generating structure...")
	(X, Y, Z, THETA1, THETA2, THETA3) = generate_elliptical_disk_system(R1, R2, Lx, Ly, Lz, lbz, ubz, ubangle, number_of_equilibration_sweeps, acceptance_probability_target, number_of_iterations_overlap_criterion, silent_mode)
	
	# Save structure in XML format.
	file_path_structure::String = join((output_dir, "/", "structure.xml"))
	write_xml_system(file_path_structure, Lx, Ly, Lz, X, Y, Z, THETA1, THETA2, THETA3, R1, R2)
	
	# Simulate diffusion.
	file_name_template = "input_template.xml"
	file_stream_template = open(file_name_template, "r")
	data_template = readstring(file_stream_template)
	data_template = replace(data_template, "PARTICLE_SYSTEM_PATH", file_path_structure)
	output_file_path = join((output_dir, "/", "diffusion.xml"))
	data_template = replace(data_template, "OUTPUT_PATH", output_file_path)
	
	file_name_input = join((output_dir, "/", "input.xml"))
	file_stream_input = open(file_name_input, "w")
	write(file_stream_input, data_template)
	close(file_stream_input)

	std_out_file_path = join((output_dir, "/", "wf_output.txt"))
	
	julia_script_path = "V:/whitefish/dev/whitefish/src/wfrun.jl"
	cmd = `julia $julia_script_path $file_name_input`
	#cmd = `mpiexec -n $Sys.CPU_CORES gesualdo $file_name_input`
	
	#proc = spawn(pipeline(cmd, stdout = std_out_file_path))
	proc = spawn(cmd)
	wait(proc)
	
	nothing
	
end

run()
