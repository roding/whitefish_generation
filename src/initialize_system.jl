function initialize_system(	particle_type::String,
							R::Array{Float64, 2},
							Lx::Float64, 
							Ly::Float64, 
							Lz::Float64,
							phi_initial::Float64,
							position_constraint_axis::String, 
							position_constraint_lower::Float64, 
							position_constraint_upper::Float64, 
							orientation_constraint_axis::Array{Float64, 1},
							orientation_constraint_lower::Float64,
							orientation_constraint_upper::Float64)

	number_of_particles::Int64 = size(R, 1)
	
	# Determine initial system size.
	if phi_initial > 0.0 # If phi_initial specified, compute (Lx, Ly, Lz), assuming cubic domain.
		if particle_type == "sphere"
			Lx = ( sum(4.0 * pi / 3.0 * R[:, 1].^3) / phi_initial )^(1/3)
		elseif particle_type == "ellipse"
			warn("Initial volume fraction specified but not used for particle type 'ellipse'.")
		elseif particle_type == "ellipsoid"
			Lx = ( sum(4.0 * pi / 3.0 * R[:, 1] .* R[:, 2] .* R[:, 3]) / phi_initial )^(1/3)
		elseif particle_type == "cuboid"
			Lx = ( sum(8.0 * R[:, 1] .* R[:, 2] .* R[:, 3]) / phi_initial )^(1/3)
		end
		Ly = Lx
		Lz = Lx
	else # If (Lx, Ly, Lz) specified, compute phi_initial.
		if particle_type == "sphere"
			phi_initial = sum(4.0 * pi / 3.0 * R[:, 1].^3 ) / (Lx * Ly * Lz)
		elseif particle_type == "ellipse"
			warn("Not computing initial volume fraction for particle type 'ellipse'.")
		elseif particle_type == "ellipsoid"
			phi_initial = sum(4.0 * pi / 3.0 * R[:, 1] .* R[:, 2] .* R[:, 3]) / (Lx * Ly * Lz)
		elseif particle_type == "cuboid"
			phi_initial = sum(8.0 * R[:, 1] .* R[:, 2] .* R[:, 3]) / (Lx * Ly * Lz)
		end
	end

	# Generate random positions.
	warn("Ignoring position constraints.")
	X::Array{Float64, 1} = Lx * rand(number_of_particles)
	Y::Array{Float64, 1} = Ly * rand(number_of_particles)
	Z::Array{Float64, 1} = Lz * rand(number_of_particles)
	
	# Generate random orientations.
	warn("Ignoring orientation constraints.")
	Q0::Array{Float64, 1} = zeros(number_of_particles)
	Q1::Array{Float64, 1} = zeros(number_of_particles)
	Q2::Array{Float64, 1} = zeros(number_of_particles)
	Q3::Array{Float64, 1} = zeros(number_of_particles)
	q0::Float64 = 0.0
	q1::Float64 = 0.0
	q2::Float64 = 0.0
	q3::Float64 = 0.0
	if particle_type != "sphere"
		for current_particle = 1:number_of_particles
			(q0, q1, q2, q3) = generate_random_unit_quaternion()
			Q0[current_particle] = q0
			Q1[current_particle] = q1
			Q2[current_particle] = q2
			Q3[current_particle] = q3
		end
	end
	
	# Characteristic/rotation matrix entries.
	A11::Array{Float64, 1} = zeros(number_of_particles)
	A12::Array{Float64, 1} = zeros(number_of_particles)
	A13::Array{Float64, 1} = zeros(number_of_particles)
	A21::Array{Float64, 1} = zeros(number_of_particles)
	A22::Array{Float64, 1} = zeros(number_of_particles)
	A23::Array{Float64, 1} = zeros(number_of_particles)
	A31::Array{Float64, 1} = zeros(number_of_particles)
	A32::Array{Float64, 1} = zeros(number_of_particles)
	A33::Array{Float64, 1} = zeros(number_of_particles)
	
	a11::Float64 = 0.0
	a12::Float64 = 0.0
	a13::Float64 = 0.0
	a21::Float64 = 0.0
	a22::Float64 = 0.0
	a23::Float64 = 0.0
	a31::Float64 = 0.0
	a32::Float64 = 0.0
	a33::Float64 = 0.0
	if particle_type == "ellipse"
		for current_particle = 1:number_of_particles
			(a11, a12, a13, a21, a22, a23, a31, a32, a33) = characteristic_matrix_ellipse(Q0[current_particle], Q1[current_particle], Q2[current_particle], Q3[current_particle], R[current_particle, 1], R[current_particle, 2])
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
	elseif particle_type == "ellipsoid"
		for current_particle = 1:number_of_particles
			(a11, a12, a13, a21, a22, a23, a31, a32, a33) = characteristic_matrix_ellipsoid(Q0[current_particle], Q1[current_particle], Q2[current_particle], Q3[current_particle], R[current_particle, 1], R[current_particle, 2], R[current_particle, 3])
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
	elseif particle_type == "cuboid"
		for current_particle = 1:number_of_particles
			(a11, a12, a13, a21, a22, a23, a31, a32, a33) = rotation_matrix(Q0[current_particle], Q1[current_particle], Q2[current_particle], Q3[current_particle])
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
	
	return (	Lx, 
			Ly, 
			Lz, 
			phi_initial, 
			X, 
			Y, 
			Z, 
			Q0, 
			Q1, 
			Q2, 
			Q3, 
			A11, 
			A12, 
			A13, 
			A21, 
			A22, 
			A23, 
			A31, 
			A32, 
			A33)
end
