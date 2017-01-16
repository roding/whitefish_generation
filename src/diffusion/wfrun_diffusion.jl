
include("../io/read_xml_key.jl")
include("../io/write_xml_key.jl")

include("../io/read_xml_input_diffusion.jl")
include("../io/read_xml_output_generation.jl")
include("../io/write_xml_output_diffusion.jl")

include("../io/print_header.jl")
include("../io/print_simulation_stats_diffusion.jl")
include("../io/print_progress_diffusion.jl")

include("ellipticaldisk/diffuse.jl")

function wfrun_diffusion()
	# Inititalization of random number generation device.
	const random_seed::Int64 = convert(Int64, time_ns())
	srand(random_seed)
	
	# Change current folder to the folder where this script lies.
	(program_file_dir::String, program_file_name::String) = splitdir(PROGRAM_FILE)
	program_file_dir = abspath(program_file_dir)
	cd(program_file_dir)
	
	# Process input arguments.
	silent_mode::Bool = false
	
	input_diffusion_path::String = ""
	input_diffusion_specified::Bool = false
	
	number_of_arguments::Int64 = length(ARGS)
	
	for current_argument = 1:number_of_arguments
		if ARGS[current_argument] == "-silent"
			silent_mode = true
		elseif isfile(ARGS[current_argument])
			input_diffusion_path = ARGS[current_argument]
			input_diffusion_specified = true
		end
	end
	
	# Print text header.
	if !silent_mode
		print_header()
	end
	
	# Error checking for input file specification.
	if !input_diffusion_specified
		println("No input file specified or specified input does not exist. Aborting.")
		return nothing
	end

	# Process input file.
	if !silent_mode
		println(join(("Reading input from file ", input_diffusion_path, "...")))
	end
	(output_generation_path::String, D0::Float64, deltat_coarse::Float64, number_of_time_points_coarse::Int64, number_of_time_points_fine_per_coarse::Int64, number_of_diffusers::Int64, number_of_cells_x::Int64, number_of_cells_y::Int64, number_of_cells_z::Int64, output_diffusion_path::String) = read_xml_input_diffusion(input_diffusion_path)
	
	# Process particle system file.
	if !silent_mode
		println(join(("Reading particle system configuration from file ", output_generation_path, "...")))
	end
	(Lx::Float64, Ly::Float64, Lz::Float64, particle_type::String, number_of_particles::Int64, X::Array{Float64, 1}, Y::Array{Float64, 1}, Z::Array{Float64, 1}, THETA1::Array{Float64, 1}, THETA2::Array{Float64, 1}, THETA3::Array{Float64, 1}, R1::Array{Float64, 1}, R2::Array{Float64, 1}) = read_xml_output_generation(output_generation_path)
		
	# Print simulation stats.
	if !silent_mode
		print_simulation_stats_diffusion(particle_type, number_of_particles, number_of_diffusers, Lx, Ly, Lz, number_of_cells_x, number_of_cells_y, number_of_cells_z)
	end
	
	# Simulate diffusion.
	t_start_ns::Int64 = convert(Int64, time_ns())
	(msd_x::Array{Float64, 1}, msd_y::Array{Float64, 1}, msd_z::Array{Float64, 1}, D0_empirical::Float64) = diffuse(X, Y, Z, THETA1, THETA2, THETA3, R1, R2, Lx, Ly, Lz, D0, deltat_coarse, number_of_time_points_coarse, number_of_time_points_fine_per_coarse, number_of_diffusers, number_of_cells_x, number_of_cells_y, number_of_cells_z, silent_mode)	
	t_finish_ns::Int64 = convert(Int64, time_ns())
	t_exec::Float64 = convert(Float64, t_finish_ns - t_start_ns) / 1e9
	
	# Write output.
	write_xml_output_diffusion(output_diffusion_path, D0, D0_empirical, deltat_coarse, number_of_time_points_coarse, msd_x, msd_y, msd_z, t_exec)
	
	# Print output information.
	if !silent_mode
		println(join(("Output written to ", output_diffusion_path, ".")))
		println("Finished.")
	end
	
	nothing
end

wfrun_diffusion()
