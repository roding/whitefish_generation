function voxel_structure(	particle_type::String,
						R::Array{Float64, 2},
						Lx::Float64,
						Ly::Float64,
						Lz::Float64,
						X::Array{Float64, 1},
						Y::Array{Float64, 1},
						Z::Array{Float64, 1},
						Q0::Array{Float64, 1},
						Q1::Array{Float64, 1},
						Q2::Array{Float64, 1},
						Q3::Array{Float64, 1},
						voxel_size::Float64)

	Nx::Int64 = convert(Int64, round(Lx/voxel_size))
	Ny::Int64 = convert(Int64, round(Ly/voxel_size))
	Nz::Int64 = convert(Int64, round(Lz/voxel_size))

	M::BitArray{3} = falses(Nx, Ny, Nz)

    X /= voxel_size
	Y /= voxel_size
	Z /= voxel_size
    R /= voxel_size
	Lx = convert(Float64, Nx)
	Ly = convert(Float64, Ny)
	Lz = convert(Float64, Nz)

	number_of_particles::Int64 = length(X)

	# Pre-computed max radii (i.e. radius of bounding sphere).
	RMAX::Array{Float64, 1} = zeros(number_of_particles)
	if particle_type == "sphere"
		for current_particle = 1:number_of_particles
			RMAX[current_particle] = R[current_particle, 1]
		end
	elseif particle_type == "ellipsoid"
		for current_particle = 1:number_of_particles
			RMAX[current_particle] = maximum( R[current_particle, 1:3] )
		end
	elseif particle_type == "cuboid"
		for current_particle = 1:number_of_particles
			RMAX[current_particle] = sqrt(R[current_particle, 1]^2 + R[current_particle, 2]^2 + R[current_particle, 3]^2)
		end
	end

	x::Float64 = 0.0
	y::Float64 = 0.0
	z::Float64 = 0.0
	x_min::Int64 = 0
	x_max::Int64 = 0
	y_min::Int64 = 0
	y_max::Int64 = 0
	z_min::Int64 = 0
	z_max::Int64 = 0
	a11::Float64 = 0.0
	a12::Float64 = 0.0
	a13::Float64 = 0.0
	a21::Float64 = 0.0
	a22::Float64 = 0.0
	a23::Float64 = 0.0
	a31::Float64 = 0.0
	a32::Float64 = 0.0
	a33::Float64 = 0.0
	f::Float64 = 0.0

	x_prim::Float64 = 0.0
	y_prim::Float64 = 0.0
	z_prim::Float64 = 0.0
	X_prim::Array{Float64, 1} = zeros(3)

	for current_particle = 1:number_of_particles
		for i = -1:1
			for j = -1:1
				for k = -1:1
					x = X[current_particle] + convert(Float64, i) * Lx
					y = Y[current_particle] + convert(Float64, j) * Ly
					z = Z[current_particle] + convert(Float64, k) * Lz
					if all([x, y, z] .< [Lx, Ly, Lz] + RMAX[current_particle]) && all([x, y, z] .> - RMAX[current_particle])
						x_min = max(1, convert(Int64, floor(x - RMAX[current_particle])))
						x_max = min(Nx, convert(Int64, ceil(x + RMAX[current_particle])))
						y_min = max(1, convert(Int64, floor(y - RMAX[current_particle])))
						y_max = min(Ny, convert(Int64, ceil(y + RMAX[current_particle])))
						z_min = max(1, convert(Int64, floor(z - RMAX[current_particle])))
						z_max = min(Nz, convert(Int64, ceil(z + RMAX[current_particle])))

						for xx = x_min:x_max
							for yy =  y_min:y_max
								for zz =  z_min:z_max
									if particle_type == "sphere"
										if (x - convert(Float64, xx))^2 + (y - convert(Float64, yy))^2 + (z - convert(Float64, zz))^2 <= R[current_particle]^2
											M[xx, yy, zz] = true
										end
									elseif particle_type == "ellipsoid"
										(a11, a12, a13, a21, a22, a23, a31, a32, a33) = characteristic_matrix_ellipsoid(Q0[current_particle], Q1[current_particle], Q2[current_particle], Q3[current_particle], R[current_particle, 1], R[current_particle, 2], R[current_particle, 3])
										f = [x - convert(Float64, xx), y - convert(Float64, yy), z - convert(Float64, zz)]' * ([a11 a12 a13; a21 a22 a23; a31 a32 a33] \ [x - convert(Float64, xx), y - convert(Float64, yy), z - convert(Float64, zz)])
										if f <= 1.0
											M[xx, yy, zz] = true
										end
									elseif particle_type == "cuboid"
										(a11, a12, a13, a21, a22, a23, a31, a32, a33) = rotation_matrix(Q0[current_particle], Q1[current_particle], Q2[current_particle], Q3[current_particle])
										X_prim = [a11 a12 a13; a21 a22 a23; a31 a32 a33] \ [x - convert(Float64, xx), y - convert(Float64, yy), z - convert(Float64, zz)]

										x_prim = X_prim[1]#a11 * (x - convert(Float64, xx)) + a12 * (y - convert(Float64, yy)) + a13 * (z - convert(Float64, zz))
										y_prim = X_prim[2]#a21 * (x - convert(Float64, xx)) + a22 * (y - convert(Float64, yy)) + a23 * (z - convert(Float64, zz))
										z_prim = X_prim[3]#a31 * (x - convert(Float64, xx)) + a32 * (y - convert(Float64, yy)) + a33 * (z - convert(Float64, zz))
										if abs(x_prim) <= R[current_particle, 1] && abs(y_prim) <= R[current_particle, 2] && abs(z_prim) <= R[current_particle, 3]
											M[xx, yy, zz] = true
										end
									end
								end
							end
						end
					end
				end
			end
		end
	end

	return M
end
