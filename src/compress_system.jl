function compress_system(	particle_type::String,
						R::Array{Float64, 2},
						Lx::Float64,
						Ly::Float64,
						Lz::Float64,
						phi_initial::Float64,
						phi_target::Float64,
						X::Array{Float64, 1},
						Y::Array{Float64, 1},
						Z::Array{Float64, 1},
						Q0::Array{Float64, 1},
						Q1::Array{Float64, 1},
						Q2::Array{Float64, 1},
						Q3::Array{Float64, 1},
						A11::Array{Float64, 1},
						A12::Array{Float64, 1},
						A13::Array{Float64, 1},
						A21::Array{Float64, 1},
						A22::Array{Float64, 1},
						A23::Array{Float64, 1},
						A31::Array{Float64, 1},
						A32::Array{Float64, 1},
						A33::Array{Float64, 1},
						position_constraint_axis::String,
						position_constraint_lower::Float64,
						position_constraint_upper::Float64,
						orientation_axis::Array{Float64, 1},
						orientation_constraint_axis::Array{Float64, 1},
						orientation_constraint_lower::Float64,
						orientation_constraint_upper::Float64,
						sigma_translation::Float64,
						sigma_translation_max::Float64,
						sigma_rotation::Float64,
						sigma_rotation_max::Float64,
						sigma_ratio::Float64,
						delta_phi::Float64,
						number_of_relaxation_sweeps_max::Int64)

	number_of_particles::Int64 = size(R, 1)

	phi::Float64 = 0.0
	if particle_type == "sphere"
		phi = sum(4.0 * pi / 3.0 * R[:, 1].^3) / (Lx * Ly * Lz)
	elseif particle_type == "ellipsoid"
		phi = sum(4.0 * pi / 3.0 * R[:, 1] .* R[:, 2] .* R[:, 3]) / (Lx * Ly * Lz)
	elseif particle_type == "cuboid"
		phi = sum(8.0 * R[:, 1] .* R[:, 2] .* R[:, 3]) / (Lx * Ly * Lz)
	end

	is_converged::Bool = false
	is_relaxed::Bool = false
	phi_prim::Float64 = 0.0
	Lx_prim::Float64 = 0.0
	Ly_prim::Float64 = 0.0
	Lz_prim::Float64 = 0.0

	X_prim::Array{Float64, 1} = zeros(number_of_particles)
	Y_prim::Array{Float64, 1} = zeros(number_of_particles)
	Z_prim::Array{Float64, 1} = zeros(number_of_particles)
	Q0_prim::Array{Float64, 1} = zeros(number_of_particles)
	Q1_prim::Array{Float64, 1} = zeros(number_of_particles)
	Q2_prim::Array{Float64, 1} = zeros(number_of_particles)
	Q3_prim::Array{Float64, 1} = zeros(number_of_particles)
	A11_prim::Array{Float64, 1} = zeros(number_of_particles)
	A12_prim::Array{Float64, 1} = zeros(number_of_particles)
	A13_prim::Array{Float64, 1} = zeros(number_of_particles)
	A21_prim::Array{Float64, 1} = zeros(number_of_particles)
	A22_prim::Array{Float64, 1} = zeros(number_of_particles)
	A23_prim::Array{Float64, 1} = zeros(number_of_particles)
	A31_prim::Array{Float64, 1} = zeros(number_of_particles)
	A32_prim::Array{Float64, 1} = zeros(number_of_particles)
	A33_prim::Array{Float64, 1} = zeros(number_of_particles)

	while !is_converged && phi < phi_target
		phi_prim = phi + delta_phi

		Lx_prim = (phi / phi_prim)^(1/3) * Lx
		Ly_prim = (phi / phi_prim)^(1/3) * Ly
		Lz_prim = (phi / phi_prim)^(1/3) * Lz

		for current_particle = 1:number_of_particles
			X_prim[current_particle] = Lx_prim / Lx * X[current_particle]
			Y_prim[current_particle] = Ly_prim / Ly * Y[current_particle]
			Z_prim[current_particle] = Lz_prim / Lz * Z[current_particle]

			Q0_prim[current_particle] = Q0[current_particle]
			Q1_prim[current_particle] = Q1[current_particle]
			Q2_prim[current_particle] = Q2[current_particle]
			Q3_prim[current_particle] = Q3[current_particle]
			A11_prim[current_particle] = A11[current_particle]
			A12_prim[current_particle] = A12[current_particle]
			A13_prim[current_particle] = A13[current_particle]
			A21_prim[current_particle] = A21[current_particle]
			A22_prim[current_particle] = A22[current_particle]
			A23_prim[current_particle] = A23[current_particle]
			A31_prim[current_particle] = A31[current_particle]
			A32_prim[current_particle] = A32[current_particle]
			A33_prim[current_particle] = A33[current_particle]
		end

		(	X_prim,
			Y_prim,
			Z_prim,
			Q0_prim,
			Q1_prim,
			Q2_prim,
			Q3_prim,
			A11_prim,
			A12_prim,
			A13_prim,
			A21_prim,
			A22_prim,
			A23_prim,
			A31_prim,
			A32_prim,
			A33_prim,
			sigma_translation,
			sigma_rotation,
			is_relaxed) = evolve_system(particle_type,
										R,
										Lx_prim,
										Ly_prim,
										Lz_prim,
										X_prim,
										Y_prim,
										Z_prim,
										Q0_prim,
										Q1_prim,
										Q2_prim,
										Q3_prim,
										A11_prim,
										A12_prim,
										A13_prim,
										A21_prim,
										A22_prim,
										A23_prim,
										A31_prim,
										A32_prim,
										A33_prim,
										position_constraint_axis,
										position_constraint_lower,
										position_constraint_upper,
										orientation_axis,
										orientation_constraint_axis,
										orientation_constraint_lower,
										orientation_constraint_upper,
										sigma_translation,
										sigma_translation_max,
										sigma_rotation,
										sigma_rotation_max,
										sigma_ratio,
										number_of_relaxation_sweeps_max,
										0)

		println(is_relaxed)
		if is_relaxed # If system has been successfully relaxed.
			phi = phi_prim
			Lx = Lx_prim
			Ly = Ly_prim
			Lz = Lz_prim

			for current_particle = 1:number_of_particles
				X[current_particle] = X_prim[current_particle]
				Y[current_particle] = Y_prim[current_particle]
				Z[current_particle] = Z_prim[current_particle]

				Q0[current_particle] = Q0_prim[current_particle]
				Q1[current_particle] = Q1_prim[current_particle]
				Q2[current_particle] = Q2_prim[current_particle]
				Q3[current_particle] = Q3_prim[current_particle]
				A11[current_particle] = A11_prim[current_particle]
				A12[current_particle] = A12_prim[current_particle]
				A13[current_particle] = A13_prim[current_particle]
				A21[current_particle] = A21_prim[current_particle]
				A22[current_particle] = A22_prim[current_particle]
				A23[current_particle] = A23_prim[current_particle]
				A31[current_particle] = A31_prim[current_particle]
				A32[current_particle] = A32_prim[current_particle]
				A33[current_particle] = A33_prim[current_particle]
			end

			println(join(["   Current volume fraction: ", string(phi)]))
		else
			is_converged = true
		end
	end

	return (Lx, Ly, Lz, phi, X, Y, Z, Q0, Q1, Q2, Q3, A11, A12, A13, A21, A22, A23, A31, A32, A33, sigma_translation, sigma_rotation)
end
