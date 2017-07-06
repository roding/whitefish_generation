include("characteristic_matrix_ellipsoid.jl")

function run_ellipsoids()

	Lx::Float64 = 100.0
	Ly::Float64 = 100.0
	Lz::Float64 = 100.0
	
	number_of_particles::Int64 = 10
	
	# Positions.
	X::Array{Float64, 1} = Lx * rand(number_of_particles)
	Y::Array{Float64, 1} = Ly * rand(number_of_particles)
	Z::Array{Float64, 1} = Lz * rand(number_of_particles)
	
	# Radii.
	R1::Array{Float64, 1} = ones(number_of_particles)
	R2::Array{Float64, 1} = ones(number_of_particles)
	R3::Array{Float64, 1} = ones(number_of_particles)
	
	# Pre-computed max radii.
	RMAX = Array(Float64, number_of_particles)
	for current_particle = 1:number_of_particles
		RMAX[current_particle] = maximum( (R1[current_particle], R2[current_particle], R3[current_particle]) )
	end
	
	# Quaternions.
	Q0::Array{Float64, 1} = zeros(number_of_particles)
	Q1::Array{Float64, 1} = zeros(number_of_particles)
	Q2::Array{Float64, 1} = zeros(number_of_particles)
	Q3::Array{Float64, 1} = zeros(number_of_particles)
	
	q_norm::Float64 = 0.0
	
	# Generate random orientations.
	for current_particle = 1:number_of_particles
		Q0[current_particle] = randn()
		Q1[current_particle] = randn()
		Q2[current_particle] = randn()
		Q3[current_particle] = randn()
		q_norm = sqrt( Q0[current_particle]^2 + Q1[current_particle]^2 + Q2[current_particle]^2 + Q3[current_particle]^2 )
		Q0[current_particle] /= q_norm
		Q1[current_particle] /= q_norm
		Q2[current_particle] /= q_norm
		Q3[current_particle] /= q_norm
	end
	
	# Simulation parameters.
	sigma_translation_max::Float64 = 0.05
	sigma_rotation_max::Float64 = 0.01
	sigma_translation::Float64 = sigma_translation_max
	sigma_rotation::Float64 = sigma_rotation_max
	
	# Ellipsoid characteristic matrix entries.
	A33::Array{Float64, 1} = zeros(number_of_particles)
	A12::Array{Float64, 1} = zeros(number_of_particles)
	A13::Array{Float64, 1} = zeros(number_of_particles)
	A21::Array{Float64, 1} = zeros(number_of_particles)
	A22::Array{Float64, 1} = zeros(number_of_particles)
	A23::Array{Float64, 1} = zeros(number_of_particles)
	A31::Array{Float64, 1} = zeros(number_of_particles)
	A32::Array{Float64, 1} = zeros(number_of_particles)
	A33::Array{Float64, 1} = zeros(number_of_particles)
	
	for current_particle = 1:number_of_particles
		(a11, a12, a13, a21, a22, a23, a31, a32, a33) = characteristic_matrix_ellipsoid(Q0[current_particle], Q1[current_particle], Q2[current_particle], Q3[current_particle], R1[current_particle], R2[current_particle], R3[current_particle])
		A11[current_particle] = a11
		A12[current_particle] = a12
		A13[current_particle] = a13
		A21[current_particle] = a21
		A22[current_particle] = a22
		A23[current_particle] = a23
		A31[current_particle] = a31
		A32[current_particle] = a32
		A33[current_particle] = a33
	end
	
	
	
	
end

run_ellipsoids()