@inbounds function generate(R1::Array{Float64,1}, R2::Array{Float64,1}, Lx::Float64, Ly::Float64, Lz::Float64, lbz::Float64, ubz::Float64, ubangle::Float64, number_of_equilibration_sweeps::Int64, acceptance_probability_target::Float64, number_of_iterations_overlap_criterion::Int64, silent::Bool)

	# Particle system parameters.
	const number_of_particles::Int64 = length(R1)
	if !silent
		println(join(["Initializing random configuration of ", string(number_of_particles), " elliptical disks..."]))
	end
	X::Array{Float64,1} = Lx * rand(number_of_particles)
	Y::Array{Float64,1} = Ly * rand(number_of_particles)
	Z::Array{Float64,1} = lbz + (ubz - lbz) * rand(number_of_particles)

	THETA1::Array{Float64,1} = zeros(number_of_particles)
	THETA2::Array{Float64,1} = zeros(number_of_particles)
	THETA3::Array{Float64,1} = zeros(number_of_particles)
	angle_is_ok::Bool = false
	if ubangle < pi # Angular constraint.
		for current_particle = 1:number_of_particles
			angle_is_ok = false
			while !angle_is_ok
					THETA1[current_particle] = 2.0 * pi * rand()
					THETA2[current_particle] = 2.0 * pi * rand()
					THETA3[current_particle] = 2.0 * pi * rand()
					angle_to_z_axis = acos((cos(THETA1[current_particle])*cos(THETA2[current_particle]))/(abs(sin(THETA2[current_particle]))^2 + abs(cos(THETA1[current_particle])*cos(THETA2[current_particle]))^2 + abs(cos(THETA2[current_particle])*sin(THETA1[current_particle]))^2)^(1/2))
					if angle_to_z_axis <= ubangle
						angle_is_ok = true
					end
				end
			end
	else # No angular constraint, the angles can be what they want.
		THETA1 = 2.0 * pi * rand(number_of_particles)
		THETA2 = 2.0 * pi * rand(number_of_particles)
		THETA3 = 2.0 * pi * rand(number_of_particles)
	end


	# Simulation parameters.
	sigma_translation_max::Float64 = 0.05
	sigma_rotation_max::Float64 = 0.01

	sigma_translation::Float64 = sigma_translation_max
	sigma_rotation::Float64 = sigma_rotation_max


	# Stored coefficients.
	const RMAX = Array(Float64,number_of_particles)
	for current_particle = 1:number_of_particles
		RMAX[current_particle] = maximum((R1[current_particle], R2[current_particle]))
	end

	cA1::Array{Float64,1} = cos(THETA1)
	cA2::Array{Float64,1} = cos(THETA2)
	sA1::Array{Float64,1} = sin(THETA1)
	sA2::Array{Float64,1} = sin(THETA2)
	cA1d::Array{Float64,1} = cos(2.0*THETA1)
	cA2d::Array{Float64,1} = cos(2.0*THETA2)
	cA3d::Array{Float64,1} = cos(2.0*THETA3)
	sA1d::Array{Float64,1} = sin(2.0*THETA1)
	sA2d::Array{Float64,1} = sin(2.0*THETA2)
	sA3d::Array{Float64,1} = sin(2.0*THETA3)
	a11::Array{Float64,1} = (cA2d*0.25 + cA3d*0.25 + (cA2d.*cA3d)*0.25 + 0.25).*R1.^2 + (cA2d*0.25 - cA3d*0.25 - (cA2d.*cA3d)*0.25 + 0.25).*R2.^2
	a12::Array{Float64,1} = ((sA1.*sA2d)*0.25 + (cA1.*cA2.*sA3d)*0.5 + (cA3d.*sA1.*sA2d)*0.25).*R1.^2 + ((sA1.*sA2d)*0.25 - (cA1.*cA2.*sA3d)*0.5 - (cA3d.*sA1.*sA2d)*0.25).*R2.^2
	a13::Array{Float64,1} = ((cA2.*sA1.*sA3d)*0.5 - (cA1.*cA3d.*sA2d)*0.25 - (cA1.*sA2d)*0.25).*R1.^2 + ((cA1.*cA3d.*sA2d)*0.25 - (cA1.*sA2d)*0.25 - (cA2.*sA1.*sA3d)*0.5).*R2.^2
	a21::Array{Float64,1} = a12
	a22::Array{Float64,1} = ((cA2.^2.*cA1d)*0.25 - (cA1d.*cA3d)*0.5 - (cA2.^2.*cA3d)*0.25 - cA2.^2.0*0.25 + (sA2.*sA1d.*sA3d)*0.5 + (cA2.^2.*cA1d.*cA3d)*0.25 + 0.5).*R1.^2 + ((cA1d.*cA3d)*0.5 + (cA2.^2.*cA1d)*0.25 + (cA2.^2.*cA3d)*0.25 - cA2.^2.0*0.25 - (sA2.*sA1d.*sA3d)*0.5 - (cA2.^2.*cA1d.*cA3d)*0.25 + 0.5).*R2.^2
	a23::Array{Float64,1} = ((cA2.^2.*sA1d)*0.25 - (cA3d.*sA1d)*0.5 - (cA1d.*sA2.*sA3d)*0.5 + (cA2.^2.*cA3d.*sA1d)*0.25).*R1.^2 + ((cA3d.*sA1d)*0.5 + (cA2.^2.*sA1d)*0.25 + (cA1d.*sA2.*sA3d)*0.5 - (cA2.^2.*cA3d.*sA1d)*0.25).*R2.^2
	a31::Array{Float64,1} = a13
	a32::Array{Float64,1} = a23
	a33::Array{Float64,1} = ((cA1d.*cA3d)*0.5 - (cA2.^2.*cA1d)*0.25 - (cA2.^2.*cA3d)*0.25 - cA2.^2.0*0.25 - (sA2.*sA1d.*sA3d)*0.5 - (cA2.^2.*cA1d.*cA3d)*0.25 + 0.5).*R1.^2 + ((cA2.^2.*cA3d)*0.25 - (cA2.^2.*cA1d)*0.25 - (cA1d.*cA3d)*0.5 - cA2.^2.0*0.25 + (sA2.*sA1d.*sA3d)*0.5 + (cA2.^2.*cA1d.*cA3d)*0.25 + 0.5).*R2.^2

	# Preallocation.
	X_star = Array(Float64,number_of_particles)
	Y_star = Array(Float64,number_of_particles)
	Z_star = Array(Float64,number_of_particles)
	THETA1_star = Array(Float64,number_of_particles)
	THETA2_star = Array(Float64,number_of_particles)
	THETA3_star = Array(Float64,number_of_particles)
	cA1_star = Array(Float64,number_of_particles)
	ca2_star = Array(Float64,number_of_particles)
	sA1_star = Array(Float64,number_of_particles)
	sA2_star = Array(Float64,number_of_particles)
	cA1d_star = Array(Float64,number_of_particles)
	cA2d_star = Array(Float64,number_of_particles)
	cA3d_star = Array(Float64,number_of_particles)
	sA1d_star = Array(Float64,number_of_particles)
	sA2d_star = Array(Float64,number_of_particles)
	sA3d_star = Array(Float64,number_of_particles)
	a11_star = Array(Float64,number_of_particles)
	a12_star = Array(Float64,number_of_particles)
	a13_star = Array(Float64,number_of_particles)
	a21_star = Array(Float64,number_of_particles)
	a22_star = Array(Float64,number_of_particles)
	a23_star = Array(Float64,number_of_particles)
	a31_star = Array(Float64,number_of_particles)
	a32_star = Array(Float64,number_of_particles)
	a33_star = Array(Float64,number_of_particles)
	energy_system::Float64 = 0.0
	energy_particle::Float64 = 0.0
	energy_particle_star::Float64 = 0.0
	current_sweep::Int64 = 0
	acceptance_probability_translation::Float64 = 0.0
	acceptance_probability_rotation::Float64 = 0.0
	xAB::Float64 = 0.0
	yAB::Float64 = 0.0
	zAB::Float64 = 0.0
	const1t::Float64 = 0.0
	const2t::Float64 = 0.0
	const3t::Float64 = 0.0
	const4t::Float64 = 0.0
	const1n::Float64 = 0.0
	const2n::Float64 = 0.0
	const3n::Float64 = 0.0
	a::Float64 = 0.0
	b::Float64 = 0.0
	c::Float64 = 0.0
	x::Float64 = 0.0
	fb::Float64 = 0.0
	fx::Float64 = 0.0
	overlapfun::Float64 = 0.0
	iter::Int64 = 0

	# Relax configuration.
	if !silent
		println("Relaxing configuration (removing all overlaps)...")
	end
	energy_system = 1.0
	while energy_system > 0.0
		current_sweep += 1
		if !silent
			println(join(["   Sweep ", string(current_sweep)]))
		end
		randn!(X_star)
		randn!(Y_star)
		randn!(Z_star)
		randn!(THETA1_star)
		randn!(THETA2_star)
		randn!(THETA3_star)
		X_star = X + sigma_translation * X_star
		Y_star = Y + sigma_translation * Y_star
		Z_star = Z + sigma_translation * Z_star
		for current_particle = 1:number_of_particles
			if X_star[current_particle] < 0.0
				X_star[current_particle] += Lx
			elseif X_star[current_particle] > Lx
				X_star[current_particle] -= Lx
			end
			if Y_star[current_particle] < 0.0
				Y_star[current_particle] += Ly
			elseif Y_star[current_particle] > Ly
				Y_star[current_particle] -= Ly
			end
			if Z_star[current_particle] < 0.0
				Z_star[current_particle] += Lz
			elseif Z_star[current_particle] > Lz
				Z_star[current_particle] -= Lz
			end
		end
		THETA1_star = THETA1 + sigma_rotation * THETA1_star
		THETA2_star = THETA2 + sigma_rotation * THETA2_star
		THETA3_star = THETA3 + sigma_rotation * THETA3_star

		cA1_star = cos(THETA1_star)
		cA2_star = cos(THETA2_star)
		sA1_star = sin(THETA1_star)
		sA2_star = sin(THETA2_star)
		cA1d_star = cos(2.0*THETA1_star)
		cA2d_star = cos(2.0*THETA2_star)
		cA3d_star = cos(2.0*THETA3_star)
		sA1d_star = sin(2.0*THETA1_star)
		sA2d_star = sin(2.0*THETA2_star)
		sA3d_star = sin(2.0*THETA3_star)
		a11_star = (cA2d_star*0.25 + cA3d_star*0.25 + (cA2d_star.*cA3d_star)*0.25 + 0.25).*R1.^2 + (cA2d_star*0.25 - cA3d_star*0.25 - (cA2d_star.*cA3d_star)*0.25 + 0.25).*R2.^2
		a12_star = ((sA1_star.*sA2d_star)*0.25 + (cA1_star.*cA2_star.*sA3d_star)*0.5 + (cA3d_star.*sA1_star.*sA2d_star)*0.25).*R1.^2 + ((sA1_star.*sA2d_star)*0.25 - (cA1_star.*cA2_star.*sA3d_star)*0.5 - (cA3d_star.*sA1_star.*sA2d_star)*0.25).*R2.^2
		a13_star = ((cA2_star.*sA1_star.*sA3d_star)*0.5 - (cA1_star.*cA3d_star.*sA2d_star)*0.25 - (cA1_star.*sA2d_star)*0.25).*R1.^2 + ((cA1_star.*cA3d_star.*sA2d_star)*0.25 - (cA1_star.*sA2d_star)*0.25 - (cA2_star.*sA1_star.*sA3d_star)*0.5).*R2.^2
		a21_star = a12_star
		a22_star = ((cA2_star.^2.*cA1d_star)*0.25 - (cA1d_star.*cA3d_star)*0.5 - (cA2_star.^2.*cA3d_star)*0.25 - cA2_star.^2.0*0.25 + (sA2_star.*sA1d_star.*sA3d_star)*0.5 + (cA2_star.^2.*cA1d_star.*cA3d_star)*0.25 + 0.5).*R1.^2 + ((cA1d_star.*cA3d_star)*0.5 + (cA2_star.^2.*cA1d_star)*0.25 + (cA2_star.^2.*cA3d_star)*0.25 - cA2_star.^2.0*0.25 - (sA2_star.*sA1d_star.*sA3d_star)*0.5 - (cA2_star.^2.*cA1d_star.*cA3d_star)*0.25 + 0.5).*R2.^2
		a23_star = ((cA2_star.^2.*sA1d_star)*0.25 - (cA3d_star.*sA1d_star)*0.5 - (cA1d_star.*sA2_star.*sA3d_star)*0.5 + (cA2_star.^2.*cA3d_star.*sA1d_star)*0.25).*R1.^2 + ((cA3d_star.*sA1d_star)*0.5 + (cA2_star.^2.*sA1d_star)*0.25 + (cA1d_star.*sA2_star.*sA3d_star)*0.5 - (cA2_star.^2.*cA3d_star.*sA1d_star)*0.25).*R2.^2
		a31_star = a13_star
		a32_star = a23_star
		a33_star = ((cA1d_star.*cA3d_star)*0.5 - (cA2_star.^2.*cA1d_star)*0.25 - (cA2_star.^2.*cA3d_star)*0.25 - cA2_star.^2.0*0.25 - (sA2_star.*sA1d_star.*sA3d_star)*0.5 - (cA2_star.^2.*cA1d_star.*cA3d_star)*0.25 + 0.5).*R1.^2 + ((cA2_star.^2.*cA3d_star)*0.25 - (cA2_star.^2.*cA1d_star)*0.25 - (cA1d_star.*cA3d_star)*0.5 - cA2_star.^2.0*0.25 + (sA2_star.*sA1d_star.*sA3d_star)*0.5 + (cA2_star.^2.*cA1d_star.*cA3d_star)*0.25 + 0.5).*R2.^2

		energy_system = 0.0
		acceptance_probability_translation = 0.0
		acceptance_probability_rotation = 0.0
		for currentA = 1:number_of_particles
			# Compute current local energy.
			energy_particle = 0.0
			for currentB = [1:currentA-1;currentA+1:number_of_particles]
				xAB = X[currentB] - X[currentA]
				if xAB < -0.5*Lx
					xAB += Lx
				elseif xAB > 0.5*Lx
					xAB -= Lx
				end

				yAB = Y[currentB] - Y[currentA]
				if yAB < -0.5*Ly
					yAB += Ly
				elseif yAB > 0.5*Ly
					yAB -= Ly
				end

				zAB = Z[currentB] - Z[currentA]
				if zAB < -0.5*Lz
					zAB += Lz
				elseif zAB > 0.5*Lz
					zAB -= Lz
				end

				if xAB^2 + yAB^2 + zAB^2 < (RMAX[currentA] + RMAX[currentB])^2
					const4t = (a12[currentA]*a21[currentA] - a11[currentA]*a22[currentA] + a11[currentA]*a22[currentB] - a12[currentA]*a21[currentB] - a21[currentA]*a12[currentB] + a22[currentA]*a11[currentB] - a11[currentB]*a22[currentB] + a12[currentB]*a21[currentB])*zAB^2 + (yAB*(a11[currentA]*a23[currentA] - a13[currentA]*a21[currentA] + a11[currentA]*a32[currentA] - a12[currentA]*a31[currentA] - a11[currentA]*a23[currentB] + a13[currentA]*a21[currentB] + a21[currentA]*a13[currentB] - a23[currentA]*a11[currentB] - a11[currentA]*a32[currentB] + a12[currentA]*a31[currentB] + a31[currentA]*a12[currentB] - a32[currentA]*a11[currentB] + a11[currentB]*a23[currentB] - a13[currentB]*a21[currentB] + a11[currentB]*a32[currentB] - a12[currentB]*a31[currentB]) - xAB*(a12[currentA]*a23[currentA] - a13[currentA]*a22[currentA] + a21[currentA]*a32[currentA] - a22[currentA]*a31[currentA] - a12[currentA]*a23[currentB] + a13[currentA]*a22[currentB] + a22[currentA]*a13[currentB] - a23[currentA]*a12[currentB] - a21[currentA]*a32[currentB] + a22[currentA]*a31[currentB] + a31[currentA]*a22[currentB] - a32[currentA]*a21[currentB] + a12[currentB]*a23[currentB] - a13[currentB]*a22[currentB] + a21[currentB]*a32[currentB] - a22[currentB]*a31[currentB]))*zAB + (a23[currentA]*a32[currentA] - a22[currentA]*a33[currentA] + a22[currentA]*a33[currentB] - a23[currentA]*a32[currentB] - a32[currentA]*a23[currentB] + a33[currentA]*a22[currentB] - a22[currentB]*a33[currentB] + a23[currentB]*a32[currentB])*xAB^2 + (a12[currentA]*a33[currentA] - a13[currentA]*a32[currentA] + a21[currentA]*a33[currentA] - a23[currentA]*a31[currentA] - a12[currentA]*a33[currentB] + a13[currentA]*a32[currentB] + a32[currentA]*a13[currentB] - a33[currentA]*a12[currentB] - a21[currentA]*a33[currentB] + a23[currentA]*a31[currentB] + a31[currentA]*a23[currentB] - a33[currentA]*a21[currentB] + a12[currentB]*a33[currentB] - a13[currentB]*a32[currentB] + a21[currentB]*a33[currentB] - a23[currentB]*a31[currentB])*xAB*yAB + (a13[currentA]*a31[currentA] - a11[currentA]*a33[currentA] + a11[currentA]*a33[currentB] - a13[currentA]*a31[currentB] - a31[currentA]*a13[currentB] + a33[currentA]*a11[currentB] - a11[currentB]*a33[currentB] + a13[currentB]*a31[currentB])*yAB^2
					const3t = (3.0*a11[currentA]*a22[currentA] - 3.0*a12[currentA]*a21[currentA] - 2.0*a11[currentA]*a22[currentB] + 2.0*a12[currentA]*a21[currentB] + 2.0*a21[currentA]*a12[currentB] - 2.0*a22[currentA]*a11[currentB] + a11[currentB]*a22[currentB] - a12[currentB]*a21[currentB])*zAB^2 + (xAB*(3.0*a12[currentA]*a23[currentA] - 3.0*a13[currentA]*a22[currentA] + 3.0*a21[currentA]*a32[currentA] - 3.0*a22[currentA]*a31[currentA] - 2.0*a12[currentA]*a23[currentB] + 2.0*a13[currentA]*a22[currentB] + 2.0*a22[currentA]*a13[currentB] - 2.0*a23[currentA]*a12[currentB] - 2.0*a21[currentA]*a32[currentB] + 2.0*a22[currentA]*a31[currentB] + 2.0*a31[currentA]*a22[currentB] - 2.0*a32[currentA]*a21[currentB] + a12[currentB]*a23[currentB] - a13[currentB]*a22[currentB] + a21[currentB]*a32[currentB] - a22[currentB]*a31[currentB]) - yAB*(3.0*a11[currentA]*a23[currentA] - 3.0*a13[currentA]*a21[currentA] + 3.0*a11[currentA]*a32[currentA] - 3.0*a12[currentA]*a31[currentA] - 2.0*a11[currentA]*a23[currentB] + 2.0*a13[currentA]*a21[currentB] + 2.0*a21[currentA]*a13[currentB] - 2.0*a23[currentA]*a11[currentB] - 2.0*a11[currentA]*a32[currentB] + 2.0*a12[currentA]*a31[currentB] + 2.0*a31[currentA]*a12[currentB] - 2.0*a32[currentA]*a11[currentB] + a11[currentB]*a23[currentB] - a13[currentB]*a21[currentB] + a11[currentB]*a32[currentB] - a12[currentB]*a31[currentB]))*zAB + (3.0*a22[currentA]*a33[currentA] - 3.0*a23[currentA]*a32[currentA] - 2.0*a22[currentA]*a33[currentB] + 2.0*a23[currentA]*a32[currentB] + 2.0*a32[currentA]*a23[currentB] - 2.0*a33[currentA]*a22[currentB] + a22[currentB]*a33[currentB] - a23[currentB]*a32[currentB])*xAB^2 + (3.0*a13[currentA]*a32[currentA] - 3.0*a12[currentA]*a33[currentA] - 3.0*a21[currentA]*a33[currentA] + 3.0*a23[currentA]*a31[currentA] + 2.0*a12[currentA]*a33[currentB] - 2.0*a13[currentA]*a32[currentB] - 2.0*a32[currentA]*a13[currentB] + 2.0*a33[currentA]*a12[currentB] + 2.0*a21[currentA]*a33[currentB] - 2.0*a23[currentA]*a31[currentB] - 2.0*a31[currentA]*a23[currentB] + 2.0*a33[currentA]*a21[currentB] - a12[currentB]*a33[currentB] + a13[currentB]*a32[currentB] - a21[currentB]*a33[currentB] + a23[currentB]*a31[currentB])*xAB*yAB + (3.0*a11[currentA]*a33[currentA] - 3.0*a13[currentA]*a31[currentA] - 2.0*a11[currentA]*a33[currentB] + 2.0*a13[currentA]*a31[currentB] + 2.0*a31[currentA]*a13[currentB] - 2.0*a33[currentA]*a11[currentB] + a11[currentB]*a33[currentB] - a13[currentB]*a31[currentB])*yAB^2
					const2t = (3.0*a12[currentA]*a21[currentA] - 3.0*a11[currentA]*a22[currentA] + a11[currentA]*a22[currentB] - a12[currentA]*a21[currentB] - a21[currentA]*a12[currentB] + a22[currentA]*a11[currentB])*zAB^2 + (yAB*(3.0*a11[currentA]*a23[currentA] - 3.0*a13[currentA]*a21[currentA] + 3.0*a11[currentA]*a32[currentA] - 3.0*a12[currentA]*a31[currentA] - a11[currentA]*a23[currentB] + a13[currentA]*a21[currentB] + a21[currentA]*a13[currentB] - a23[currentA]*a11[currentB] - a11[currentA]*a32[currentB] + a12[currentA]*a31[currentB] + a31[currentA]*a12[currentB] - a32[currentA]*a11[currentB]) - xAB*(3.0*a12[currentA]*a23[currentA] - 3.0*a13[currentA]*a22[currentA] + 3.0*a21[currentA]*a32[currentA] - 3.0*a22[currentA]*a31[currentA] - a12[currentA]*a23[currentB] + a13[currentA]*a22[currentB] + a22[currentA]*a13[currentB] - a23[currentA]*a12[currentB] - a21[currentA]*a32[currentB] + a22[currentA]*a31[currentB] + a31[currentA]*a22[currentB] - a32[currentA]*a21[currentB]))*zAB + (3.0*a23[currentA]*a32[currentA] - 3.0*a22[currentA]*a33[currentA] + a22[currentA]*a33[currentB] - a23[currentA]*a32[currentB] - a32[currentA]*a23[currentB] + a33[currentA]*a22[currentB])*xAB^2 + (3.0*a12[currentA]*a33[currentA] - 3.0*a13[currentA]*a32[currentA] + 3.0*a21[currentA]*a33[currentA] - 3.0*a23[currentA]*a31[currentA] - a12[currentA]*a33[currentB] + a13[currentA]*a32[currentB] + a32[currentA]*a13[currentB] - a33[currentA]*a12[currentB] - a21[currentA]*a33[currentB] + a23[currentA]*a31[currentB] + a31[currentA]*a23[currentB] - a33[currentA]*a21[currentB])*xAB*yAB + (3.0*a13[currentA]*a31[currentA] - 3.0*a11[currentA]*a33[currentA] + a11[currentA]*a33[currentB] - a13[currentA]*a31[currentB] - a31[currentA]*a13[currentB] + a33[currentA]*a11[currentB])*yAB^2
					const1t = (a11[currentA]*a22[currentA] - a12[currentA]*a21[currentA])*zAB^2 + (xAB*(a12[currentA]*a23[currentA] - a13[currentA]*a22[currentA] + a21[currentA]*a32[currentA] - a22[currentA]*a31[currentA]) - yAB*(a11[currentA]*a23[currentA] - a13[currentA]*a21[currentA] + a11[currentA]*a32[currentA] - a12[currentA]*a31[currentA]))*zAB + (a22[currentA]*a33[currentA] - a23[currentA]*a32[currentA])*xAB^2 + (a13[currentA]*a32[currentA] - a12[currentA]*a33[currentA] - a21[currentA]*a33[currentA] + a23[currentA]*a31[currentA])*xAB*yAB + (a11[currentA]*a33[currentA] - a13[currentA]*a31[currentA])*yAB^2
					const3n = a11[currentA]*a23[currentA]*a32[currentA] - a11[currentA]*a22[currentA]*a33[currentA] + a12[currentA]*a21[currentA]*a33[currentA] - a12[currentA]*a23[currentA]*a31[currentA] - a13[currentA]*a21[currentA]*a32[currentA] + a13[currentA]*a22[currentA]*a31[currentA] + a11[currentA]*a22[currentA]*a33[currentB] - a11[currentA]*a23[currentA]*a32[currentB] - a11[currentA]*a32[currentA]*a23[currentB] + a11[currentA]*a33[currentA]*a22[currentB] - a12[currentA]*a21[currentA]*a33[currentB] + a12[currentA]*a23[currentA]*a31[currentB] + a12[currentA]*a31[currentA]*a23[currentB] - a12[currentA]*a33[currentA]*a21[currentB] + a13[currentA]*a21[currentA]*a32[currentB] - a13[currentA]*a22[currentA]*a31[currentB] - a13[currentA]*a31[currentA]*a22[currentB] + a13[currentA]*a32[currentA]*a21[currentB] + a21[currentA]*a32[currentA]*a13[currentB] - a21[currentA]*a33[currentA]*a12[currentB] - a22[currentA]*a31[currentA]*a13[currentB] + a22[currentA]*a33[currentA]*a11[currentB] + a23[currentA]*a31[currentA]*a12[currentB] - a23[currentA]*a32[currentA]*a11[currentB] - a11[currentA]*a22[currentB]*a33[currentB] + a11[currentA]*a23[currentB]*a32[currentB] + a12[currentA]*a21[currentB]*a33[currentB] - a12[currentA]*a23[currentB]*a31[currentB] - a13[currentA]*a21[currentB]*a32[currentB] + a13[currentA]*a22[currentB]*a31[currentB] + a21[currentA]*a12[currentB]*a33[currentB] - a21[currentA]*a13[currentB]*a32[currentB] - a22[currentA]*a11[currentB]*a33[currentB] + a22[currentA]*a13[currentB]*a31[currentB] + a23[currentA]*a11[currentB]*a32[currentB] - a23[currentA]*a12[currentB]*a31[currentB] - a31[currentA]*a12[currentB]*a23[currentB] + a31[currentA]*a13[currentB]*a22[currentB] + a32[currentA]*a11[currentB]*a23[currentB] - a32[currentA]*a13[currentB]*a21[currentB] - a33[currentA]*a11[currentB]*a22[currentB] + a33[currentA]*a12[currentB]*a21[currentB] + a11[currentB]*a22[currentB]*a33[currentB] - a11[currentB]*a23[currentB]*a32[currentB] - a12[currentB]*a21[currentB]*a33[currentB] + a12[currentB]*a23[currentB]*a31[currentB] + a13[currentB]*a21[currentB]*a32[currentB] - a13[currentB]*a22[currentB]*a31[currentB]
					const2n = 3.0*a11[currentA]*a22[currentA]*a33[currentA] - 3.0*a11[currentA]*a23[currentA]*a32[currentA] - 3.0*a12[currentA]*a21[currentA]*a33[currentA] + 3.0*a12[currentA]*a23[currentA]*a31[currentA] + 3.0*a13[currentA]*a21[currentA]*a32[currentA] - 3.0*a13[currentA]*a22[currentA]*a31[currentA] - 2.0*a11[currentA]*a22[currentA]*a33[currentB] + 2.0*a11[currentA]*a23[currentA]*a32[currentB] + 2.0*a11[currentA]*a32[currentA]*a23[currentB] - 2.0*a11[currentA]*a33[currentA]*a22[currentB] + 2.0*a12[currentA]*a21[currentA]*a33[currentB] - 2.0*a12[currentA]*a23[currentA]*a31[currentB] - 2.0*a12[currentA]*a31[currentA]*a23[currentB] + 2.0*a12[currentA]*a33[currentA]*a21[currentB] - 2.0*a13[currentA]*a21[currentA]*a32[currentB] + 2.0*a13[currentA]*a22[currentA]*a31[currentB] + 2.0*a13[currentA]*a31[currentA]*a22[currentB] - 2.0*a13[currentA]*a32[currentA]*a21[currentB] - 2.0*a21[currentA]*a32[currentA]*a13[currentB] + 2.0*a21[currentA]*a33[currentA]*a12[currentB] + 2.0*a22[currentA]*a31[currentA]*a13[currentB] - 2.0*a22[currentA]*a33[currentA]*a11[currentB] - 2.0*a23[currentA]*a31[currentA]*a12[currentB] + 2.0*a23[currentA]*a32[currentA]*a11[currentB] + a11[currentA]*a22[currentB]*a33[currentB] - a11[currentA]*a23[currentB]*a32[currentB] - a12[currentA]*a21[currentB]*a33[currentB] + a12[currentA]*a23[currentB]*a31[currentB] + a13[currentA]*a21[currentB]*a32[currentB] - a13[currentA]*a22[currentB]*a31[currentB] - a21[currentA]*a12[currentB]*a33[currentB] + a21[currentA]*a13[currentB]*a32[currentB] + a22[currentA]*a11[currentB]*a33[currentB] - a22[currentA]*a13[currentB]*a31[currentB] - a23[currentA]*a11[currentB]*a32[currentB] + a23[currentA]*a12[currentB]*a31[currentB] + a31[currentA]*a12[currentB]*a23[currentB] - a31[currentA]*a13[currentB]*a22[currentB] - a32[currentA]*a11[currentB]*a23[currentB] + a32[currentA]*a13[currentB]*a21[currentB] + a33[currentA]*a11[currentB]*a22[currentB] - a33[currentA]*a12[currentB]*a21[currentB]
					const1n = 3.0*a11[currentA]*a23[currentA]*a32[currentA] - 3.0*a11[currentA]*a22[currentA]*a33[currentA] + 3.0*a12[currentA]*a21[currentA]*a33[currentA] - 3.0*a12[currentA]*a23[currentA]*a31[currentA] - 3.0*a13[currentA]*a21[currentA]*a32[currentA] + 3.0*a13[currentA]*a22[currentA]*a31[currentA] + a11[currentA]*a22[currentA]*a33[currentB] - a11[currentA]*a23[currentA]*a32[currentB] - a11[currentA]*a32[currentA]*a23[currentB] + a11[currentA]*a33[currentA]*a22[currentB] - a12[currentA]*a21[currentA]*a33[currentB] + a12[currentA]*a23[currentA]*a31[currentB] + a12[currentA]*a31[currentA]*a23[currentB] - a12[currentA]*a33[currentA]*a21[currentB] + a13[currentA]*a21[currentA]*a32[currentB] - a13[currentA]*a22[currentA]*a31[currentB] - a13[currentA]*a31[currentA]*a22[currentB] + a13[currentA]*a32[currentA]*a21[currentB] + a21[currentA]*a32[currentA]*a13[currentB] - a21[currentA]*a33[currentA]*a12[currentB] - a22[currentA]*a31[currentA]*a13[currentB] + a22[currentA]*a33[currentA]*a11[currentB] + a23[currentA]*a31[currentA]*a12[currentB] - a23[currentA]*a32[currentA]*a11[currentB]

					b = 0.5

					fb = (((const4t*b+const3t)*b+const2t)*b+const1t)/((const3n*b+const2n)*b+const1n)
					if fb >= 1.0
						overlapfun = 1.0
					else
						iter = 1
						a = 0.0
						c = 1.0
						while fb < 1 && iter < number_of_iterations_overlap_criterion
							# Pick point in upper half.
							x = 0.5*(b+c)
							fx = (((const4t*x+const3t)*x+const2t)*x+const1t)/((const3n*x+const2n)*x+const1n)
							if fx < fb
								c = x
							else
								b = x
								fb = fx
							end
							# Pick point in lower half.
							x = 0.5*(a+b)
							fx = (((const4t*x+const3t)*x+const2t)*x+const1t)/((const3n*x+const2n)*x+const1n)

							if fx < fb
								a = x
							else
								b = x
								fb = fx
							end

							iter = iter + 1
						end
						overlapfun = fb
					end
					if overlapfun < 1.0
						energy_particle += (1.0-overlapfun)^2
					end
				end
			end

			# Compute new local energy with translation.
			energy_particle_star = 0.0
			for currentB = [1:currentA-1;currentA+1:number_of_particles]
				xAB = X[currentB] - X_star[currentA]
				if xAB < -0.5*Lx
					xAB += Lx
				elseif xAB > 0.5*Lx
					xAB -= Lx
				end

				yAB = Y[currentB] - Y_star[currentA]
				if yAB < -0.5*Ly
					yAB += Ly
				elseif yAB > 0.5*Ly
					yAB -= Ly
				end

				zAB = Z[currentB] - Z_star[currentA]
				if zAB < -0.5*Lz
					zAB += Lz
				elseif zAB > 0.5*Lz
					zAB -= Lz
				end

				if xAB^2 + yAB^2 + zAB^2 < (RMAX[currentA] + RMAX[currentB])^2
					const4t = (a12[currentA]*a21[currentA] - a11[currentA]*a22[currentA] + a11[currentA]*a22[currentB] - a12[currentA]*a21[currentB] - a21[currentA]*a12[currentB] + a22[currentA]*a11[currentB] - a11[currentB]*a22[currentB] + a12[currentB]*a21[currentB])*zAB^2 + (yAB*(a11[currentA]*a23[currentA] - a13[currentA]*a21[currentA] + a11[currentA]*a32[currentA] - a12[currentA]*a31[currentA] - a11[currentA]*a23[currentB] + a13[currentA]*a21[currentB] + a21[currentA]*a13[currentB] - a23[currentA]*a11[currentB] - a11[currentA]*a32[currentB] + a12[currentA]*a31[currentB] + a31[currentA]*a12[currentB] - a32[currentA]*a11[currentB] + a11[currentB]*a23[currentB] - a13[currentB]*a21[currentB] + a11[currentB]*a32[currentB] - a12[currentB]*a31[currentB]) - xAB*(a12[currentA]*a23[currentA] - a13[currentA]*a22[currentA] + a21[currentA]*a32[currentA] - a22[currentA]*a31[currentA] - a12[currentA]*a23[currentB] + a13[currentA]*a22[currentB] + a22[currentA]*a13[currentB] - a23[currentA]*a12[currentB] - a21[currentA]*a32[currentB] + a22[currentA]*a31[currentB] + a31[currentA]*a22[currentB] - a32[currentA]*a21[currentB] + a12[currentB]*a23[currentB] - a13[currentB]*a22[currentB] + a21[currentB]*a32[currentB] - a22[currentB]*a31[currentB]))*zAB + (a23[currentA]*a32[currentA] - a22[currentA]*a33[currentA] + a22[currentA]*a33[currentB] - a23[currentA]*a32[currentB] - a32[currentA]*a23[currentB] + a33[currentA]*a22[currentB] - a22[currentB]*a33[currentB] + a23[currentB]*a32[currentB])*xAB^2 + (a12[currentA]*a33[currentA] - a13[currentA]*a32[currentA] + a21[currentA]*a33[currentA] - a23[currentA]*a31[currentA] - a12[currentA]*a33[currentB] + a13[currentA]*a32[currentB] + a32[currentA]*a13[currentB] - a33[currentA]*a12[currentB] - a21[currentA]*a33[currentB] + a23[currentA]*a31[currentB] + a31[currentA]*a23[currentB] - a33[currentA]*a21[currentB] + a12[currentB]*a33[currentB] - a13[currentB]*a32[currentB] + a21[currentB]*a33[currentB] - a23[currentB]*a31[currentB])*xAB*yAB + (a13[currentA]*a31[currentA] - a11[currentA]*a33[currentA] + a11[currentA]*a33[currentB] - a13[currentA]*a31[currentB] - a31[currentA]*a13[currentB] + a33[currentA]*a11[currentB] - a11[currentB]*a33[currentB] + a13[currentB]*a31[currentB])*yAB^2
					const3t = (3.0*a11[currentA]*a22[currentA] - 3.0*a12[currentA]*a21[currentA] - 2.0*a11[currentA]*a22[currentB] + 2.0*a12[currentA]*a21[currentB] + 2.0*a21[currentA]*a12[currentB] - 2.0*a22[currentA]*a11[currentB] + a11[currentB]*a22[currentB] - a12[currentB]*a21[currentB])*zAB^2 + (xAB*(3.0*a12[currentA]*a23[currentA] - 3.0*a13[currentA]*a22[currentA] + 3.0*a21[currentA]*a32[currentA] - 3.0*a22[currentA]*a31[currentA] - 2.0*a12[currentA]*a23[currentB] + 2.0*a13[currentA]*a22[currentB] + 2.0*a22[currentA]*a13[currentB] - 2.0*a23[currentA]*a12[currentB] - 2.0*a21[currentA]*a32[currentB] + 2.0*a22[currentA]*a31[currentB] + 2.0*a31[currentA]*a22[currentB] - 2.0*a32[currentA]*a21[currentB] + a12[currentB]*a23[currentB] - a13[currentB]*a22[currentB] + a21[currentB]*a32[currentB] - a22[currentB]*a31[currentB]) - yAB*(3.0*a11[currentA]*a23[currentA] - 3.0*a13[currentA]*a21[currentA] + 3.0*a11[currentA]*a32[currentA] - 3.0*a12[currentA]*a31[currentA] - 2.0*a11[currentA]*a23[currentB] + 2.0*a13[currentA]*a21[currentB] + 2.0*a21[currentA]*a13[currentB] - 2.0*a23[currentA]*a11[currentB] - 2.0*a11[currentA]*a32[currentB] + 2.0*a12[currentA]*a31[currentB] + 2.0*a31[currentA]*a12[currentB] - 2.0*a32[currentA]*a11[currentB] + a11[currentB]*a23[currentB] - a13[currentB]*a21[currentB] + a11[currentB]*a32[currentB] - a12[currentB]*a31[currentB]))*zAB + (3.0*a22[currentA]*a33[currentA] - 3.0*a23[currentA]*a32[currentA] - 2.0*a22[currentA]*a33[currentB] + 2.0*a23[currentA]*a32[currentB] + 2.0*a32[currentA]*a23[currentB] - 2.0*a33[currentA]*a22[currentB] + a22[currentB]*a33[currentB] - a23[currentB]*a32[currentB])*xAB^2 + (3.0*a13[currentA]*a32[currentA] - 3.0*a12[currentA]*a33[currentA] - 3.0*a21[currentA]*a33[currentA] + 3.0*a23[currentA]*a31[currentA] + 2.0*a12[currentA]*a33[currentB] - 2.0*a13[currentA]*a32[currentB] - 2.0*a32[currentA]*a13[currentB] + 2.0*a33[currentA]*a12[currentB] + 2.0*a21[currentA]*a33[currentB] - 2.0*a23[currentA]*a31[currentB] - 2.0*a31[currentA]*a23[currentB] + 2.0*a33[currentA]*a21[currentB] - a12[currentB]*a33[currentB] + a13[currentB]*a32[currentB] - a21[currentB]*a33[currentB] + a23[currentB]*a31[currentB])*xAB*yAB + (3.0*a11[currentA]*a33[currentA] - 3.0*a13[currentA]*a31[currentA] - 2.0*a11[currentA]*a33[currentB] + 2.0*a13[currentA]*a31[currentB] + 2.0*a31[currentA]*a13[currentB] - 2.0*a33[currentA]*a11[currentB] + a11[currentB]*a33[currentB] - a13[currentB]*a31[currentB])*yAB^2
					const2t = (3.0*a12[currentA]*a21[currentA] - 3.0*a11[currentA]*a22[currentA] + a11[currentA]*a22[currentB] - a12[currentA]*a21[currentB] - a21[currentA]*a12[currentB] + a22[currentA]*a11[currentB])*zAB^2 + (yAB*(3.0*a11[currentA]*a23[currentA] - 3.0*a13[currentA]*a21[currentA] + 3.0*a11[currentA]*a32[currentA] - 3.0*a12[currentA]*a31[currentA] - a11[currentA]*a23[currentB] + a13[currentA]*a21[currentB] + a21[currentA]*a13[currentB] - a23[currentA]*a11[currentB] - a11[currentA]*a32[currentB] + a12[currentA]*a31[currentB] + a31[currentA]*a12[currentB] - a32[currentA]*a11[currentB]) - xAB*(3.0*a12[currentA]*a23[currentA] - 3.0*a13[currentA]*a22[currentA] + 3.0*a21[currentA]*a32[currentA] - 3.0*a22[currentA]*a31[currentA] - a12[currentA]*a23[currentB] + a13[currentA]*a22[currentB] + a22[currentA]*a13[currentB] - a23[currentA]*a12[currentB] - a21[currentA]*a32[currentB] + a22[currentA]*a31[currentB] + a31[currentA]*a22[currentB] - a32[currentA]*a21[currentB]))*zAB + (3.0*a23[currentA]*a32[currentA] - 3.0*a22[currentA]*a33[currentA] + a22[currentA]*a33[currentB] - a23[currentA]*a32[currentB] - a32[currentA]*a23[currentB] + a33[currentA]*a22[currentB])*xAB^2 + (3.0*a12[currentA]*a33[currentA] - 3.0*a13[currentA]*a32[currentA] + 3.0*a21[currentA]*a33[currentA] - 3.0*a23[currentA]*a31[currentA] - a12[currentA]*a33[currentB] + a13[currentA]*a32[currentB] + a32[currentA]*a13[currentB] - a33[currentA]*a12[currentB] - a21[currentA]*a33[currentB] + a23[currentA]*a31[currentB] + a31[currentA]*a23[currentB] - a33[currentA]*a21[currentB])*xAB*yAB + (3.0*a13[currentA]*a31[currentA] - 3.0*a11[currentA]*a33[currentA] + a11[currentA]*a33[currentB] - a13[currentA]*a31[currentB] - a31[currentA]*a13[currentB] + a33[currentA]*a11[currentB])*yAB^2
					const1t = (a11[currentA]*a22[currentA] - a12[currentA]*a21[currentA])*zAB^2 + (xAB*(a12[currentA]*a23[currentA] - a13[currentA]*a22[currentA] + a21[currentA]*a32[currentA] - a22[currentA]*a31[currentA]) - yAB*(a11[currentA]*a23[currentA] - a13[currentA]*a21[currentA] + a11[currentA]*a32[currentA] - a12[currentA]*a31[currentA]))*zAB + (a22[currentA]*a33[currentA] - a23[currentA]*a32[currentA])*xAB^2 + (a13[currentA]*a32[currentA] - a12[currentA]*a33[currentA] - a21[currentA]*a33[currentA] + a23[currentA]*a31[currentA])*xAB*yAB + (a11[currentA]*a33[currentA] - a13[currentA]*a31[currentA])*yAB^2
					const3n = a11[currentA]*a23[currentA]*a32[currentA] - a11[currentA]*a22[currentA]*a33[currentA] + a12[currentA]*a21[currentA]*a33[currentA] - a12[currentA]*a23[currentA]*a31[currentA] - a13[currentA]*a21[currentA]*a32[currentA] + a13[currentA]*a22[currentA]*a31[currentA] + a11[currentA]*a22[currentA]*a33[currentB] - a11[currentA]*a23[currentA]*a32[currentB] - a11[currentA]*a32[currentA]*a23[currentB] + a11[currentA]*a33[currentA]*a22[currentB] - a12[currentA]*a21[currentA]*a33[currentB] + a12[currentA]*a23[currentA]*a31[currentB] + a12[currentA]*a31[currentA]*a23[currentB] - a12[currentA]*a33[currentA]*a21[currentB] + a13[currentA]*a21[currentA]*a32[currentB] - a13[currentA]*a22[currentA]*a31[currentB] - a13[currentA]*a31[currentA]*a22[currentB] + a13[currentA]*a32[currentA]*a21[currentB] + a21[currentA]*a32[currentA]*a13[currentB] - a21[currentA]*a33[currentA]*a12[currentB] - a22[currentA]*a31[currentA]*a13[currentB] + a22[currentA]*a33[currentA]*a11[currentB] + a23[currentA]*a31[currentA]*a12[currentB] - a23[currentA]*a32[currentA]*a11[currentB] - a11[currentA]*a22[currentB]*a33[currentB] + a11[currentA]*a23[currentB]*a32[currentB] + a12[currentA]*a21[currentB]*a33[currentB] - a12[currentA]*a23[currentB]*a31[currentB] - a13[currentA]*a21[currentB]*a32[currentB] + a13[currentA]*a22[currentB]*a31[currentB] + a21[currentA]*a12[currentB]*a33[currentB] - a21[currentA]*a13[currentB]*a32[currentB] - a22[currentA]*a11[currentB]*a33[currentB] + a22[currentA]*a13[currentB]*a31[currentB] + a23[currentA]*a11[currentB]*a32[currentB] - a23[currentA]*a12[currentB]*a31[currentB] - a31[currentA]*a12[currentB]*a23[currentB] + a31[currentA]*a13[currentB]*a22[currentB] + a32[currentA]*a11[currentB]*a23[currentB] - a32[currentA]*a13[currentB]*a21[currentB] - a33[currentA]*a11[currentB]*a22[currentB] + a33[currentA]*a12[currentB]*a21[currentB] + a11[currentB]*a22[currentB]*a33[currentB] - a11[currentB]*a23[currentB]*a32[currentB] - a12[currentB]*a21[currentB]*a33[currentB] + a12[currentB]*a23[currentB]*a31[currentB] + a13[currentB]*a21[currentB]*a32[currentB] - a13[currentB]*a22[currentB]*a31[currentB]
					const2n = 3.0*a11[currentA]*a22[currentA]*a33[currentA] - 3.0*a11[currentA]*a23[currentA]*a32[currentA] - 3.0*a12[currentA]*a21[currentA]*a33[currentA] + 3.0*a12[currentA]*a23[currentA]*a31[currentA] + 3.0*a13[currentA]*a21[currentA]*a32[currentA] - 3.0*a13[currentA]*a22[currentA]*a31[currentA] - 2.0*a11[currentA]*a22[currentA]*a33[currentB] + 2.0*a11[currentA]*a23[currentA]*a32[currentB] + 2.0*a11[currentA]*a32[currentA]*a23[currentB] - 2.0*a11[currentA]*a33[currentA]*a22[currentB] + 2.0*a12[currentA]*a21[currentA]*a33[currentB] - 2.0*a12[currentA]*a23[currentA]*a31[currentB] - 2.0*a12[currentA]*a31[currentA]*a23[currentB] + 2.0*a12[currentA]*a33[currentA]*a21[currentB] - 2.0*a13[currentA]*a21[currentA]*a32[currentB] + 2.0*a13[currentA]*a22[currentA]*a31[currentB] + 2.0*a13[currentA]*a31[currentA]*a22[currentB] - 2.0*a13[currentA]*a32[currentA]*a21[currentB] - 2.0*a21[currentA]*a32[currentA]*a13[currentB] + 2.0*a21[currentA]*a33[currentA]*a12[currentB] + 2.0*a22[currentA]*a31[currentA]*a13[currentB] - 2.0*a22[currentA]*a33[currentA]*a11[currentB] - 2.0*a23[currentA]*a31[currentA]*a12[currentB] + 2.0*a23[currentA]*a32[currentA]*a11[currentB] + a11[currentA]*a22[currentB]*a33[currentB] - a11[currentA]*a23[currentB]*a32[currentB] - a12[currentA]*a21[currentB]*a33[currentB] + a12[currentA]*a23[currentB]*a31[currentB] + a13[currentA]*a21[currentB]*a32[currentB] - a13[currentA]*a22[currentB]*a31[currentB] - a21[currentA]*a12[currentB]*a33[currentB] + a21[currentA]*a13[currentB]*a32[currentB] + a22[currentA]*a11[currentB]*a33[currentB] - a22[currentA]*a13[currentB]*a31[currentB] - a23[currentA]*a11[currentB]*a32[currentB] + a23[currentA]*a12[currentB]*a31[currentB] + a31[currentA]*a12[currentB]*a23[currentB] - a31[currentA]*a13[currentB]*a22[currentB] - a32[currentA]*a11[currentB]*a23[currentB] + a32[currentA]*a13[currentB]*a21[currentB] + a33[currentA]*a11[currentB]*a22[currentB] - a33[currentA]*a12[currentB]*a21[currentB]
					const1n = 3.0*a11[currentA]*a23[currentA]*a32[currentA] - 3.0*a11[currentA]*a22[currentA]*a33[currentA] + 3.0*a12[currentA]*a21[currentA]*a33[currentA] - 3.0*a12[currentA]*a23[currentA]*a31[currentA] - 3.0*a13[currentA]*a21[currentA]*a32[currentA] + 3.0*a13[currentA]*a22[currentA]*a31[currentA] + a11[currentA]*a22[currentA]*a33[currentB] - a11[currentA]*a23[currentA]*a32[currentB] - a11[currentA]*a32[currentA]*a23[currentB] + a11[currentA]*a33[currentA]*a22[currentB] - a12[currentA]*a21[currentA]*a33[currentB] + a12[currentA]*a23[currentA]*a31[currentB] + a12[currentA]*a31[currentA]*a23[currentB] - a12[currentA]*a33[currentA]*a21[currentB] + a13[currentA]*a21[currentA]*a32[currentB] - a13[currentA]*a22[currentA]*a31[currentB] - a13[currentA]*a31[currentA]*a22[currentB] + a13[currentA]*a32[currentA]*a21[currentB] + a21[currentA]*a32[currentA]*a13[currentB] - a21[currentA]*a33[currentA]*a12[currentB] - a22[currentA]*a31[currentA]*a13[currentB] + a22[currentA]*a33[currentA]*a11[currentB] + a23[currentA]*a31[currentA]*a12[currentB] - a23[currentA]*a32[currentA]*a11[currentB]

					b = 0.5

					fb = (((const4t*b+const3t)*b+const2t)*b+const1t)/((const3n*b+const2n)*b+const1n)
					if fb >= 1.0
						overlapfun = 1.0
					else
						iter = 1
						a = 0.0
						c = 1.0
						while fb < 1 && iter < number_of_iterations_overlap_criterion
							# Pick point in upper half.
							x = 0.5*(b+c)
							fx = (((const4t*x+const3t)*x+const2t)*x+const1t)/((const3n*x+const2n)*x+const1n)
							if fx < fb
								c = x
							else
								b = x
								fb = fx
							end
							# Pick point in lower half.
							x = 0.5*(a+b)
							fx = (((const4t*x+const3t)*x+const2t)*x+const1t)/((const3n*x+const2n)*x+const1n)

							if fx < fb
								a = x
							else
								b = x
								fb = fx
							end

							iter = iter + 1
						end
						overlapfun = fb
					end
					if overlapfun < 1.0
						energy_particle_star += (1.0-overlapfun)^2
					end
				end
			end

			if Z_star[currentA] < lbz || Z_star[currentA] > ubz
				energy_particle_star = energy_particle + 1.0
			end

			if energy_particle_star <= energy_particle
				X[currentA] = X_star[currentA]
				Y[currentA] = Y_star[currentA]
				Z[currentA] = Z_star[currentA]
				acceptance_probability_translation += 1.0
				energy_particle = energy_particle_star
			end

			# Compute new local energy with rotation.
			energy_particle_star = 0.0
			for currentB = [1:currentA-1;currentA+1:number_of_particles]
				xAB = X[currentB] - X[currentA]
				if xAB < -0.5*Lx
					xAB += Lx
				elseif xAB > 0.5*Lx
					xAB -= Lx
				end

				yAB = Y[currentB] - Y[currentA]
				if yAB < -0.5*Ly
					yAB += Ly
				elseif yAB > 0.5*Ly
					yAB -= Ly
				end

				zAB = Z[currentB] - Z[currentA]
				if zAB < -0.5*Lz
					zAB += Lz
				elseif zAB > 0.5*Lz
					zAB -= Lz
				end

				if xAB^2 + yAB^2 + zAB^2 < (RMAX[currentA] + RMAX[currentB])^2
					const4t = (a12_star[currentA]*a21_star[currentA] - a11_star[currentA]*a22_star[currentA] + a11_star[currentA]*a22[currentB] - a12_star[currentA]*a21[currentB] - a21_star[currentA]*a12[currentB] + a22_star[currentA]*a11[currentB] - a11[currentB]*a22[currentB] + a12[currentB]*a21[currentB])*zAB^2 + (yAB*(a11_star[currentA]*a23_star[currentA] - a13_star[currentA]*a21_star[currentA] + a11_star[currentA]*a32_star[currentA] - a12_star[currentA]*a31_star[currentA] - a11_star[currentA]*a23[currentB] + a13_star[currentA]*a21[currentB] + a21_star[currentA]*a13[currentB] - a23_star[currentA]*a11[currentB] - a11_star[currentA]*a32[currentB] + a12_star[currentA]*a31[currentB] + a31_star[currentA]*a12[currentB] - a32_star[currentA]*a11[currentB] + a11[currentB]*a23[currentB] - a13[currentB]*a21[currentB] + a11[currentB]*a32[currentB] - a12[currentB]*a31[currentB]) - xAB*(a12_star[currentA]*a23_star[currentA] - a13_star[currentA]*a22_star[currentA] + a21_star[currentA]*a32_star[currentA] - a22_star[currentA]*a31_star[currentA] - a12_star[currentA]*a23[currentB] + a13_star[currentA]*a22[currentB] + a22_star[currentA]*a13[currentB] - a23_star[currentA]*a12[currentB] - a21_star[currentA]*a32[currentB] + a22_star[currentA]*a31[currentB] + a31_star[currentA]*a22[currentB] - a32_star[currentA]*a21[currentB] + a12[currentB]*a23[currentB] - a13[currentB]*a22[currentB] + a21[currentB]*a32[currentB] - a22[currentB]*a31[currentB]))*zAB + (a23_star[currentA]*a32_star[currentA] - a22_star[currentA]*a33_star[currentA] + a22_star[currentA]*a33[currentB] - a23_star[currentA]*a32[currentB] - a32_star[currentA]*a23[currentB] + a33_star[currentA]*a22[currentB] - a22[currentB]*a33[currentB] + a23[currentB]*a32[currentB])*xAB^2 + (a12_star[currentA]*a33_star[currentA] - a13_star[currentA]*a32_star[currentA] + a21_star[currentA]*a33_star[currentA] - a23_star[currentA]*a31_star[currentA] - a12_star[currentA]*a33[currentB] + a13_star[currentA]*a32[currentB] + a32_star[currentA]*a13[currentB] - a33_star[currentA]*a12[currentB] - a21_star[currentA]*a33[currentB] + a23_star[currentA]*a31[currentB] + a31_star[currentA]*a23[currentB] - a33_star[currentA]*a21[currentB] + a12[currentB]*a33[currentB] - a13[currentB]*a32[currentB] + a21[currentB]*a33[currentB] - a23[currentB]*a31[currentB])*xAB*yAB + (a13_star[currentA]*a31_star[currentA] - a11_star[currentA]*a33_star[currentA] + a11_star[currentA]*a33[currentB] - a13_star[currentA]*a31[currentB] - a31_star[currentA]*a13[currentB] + a33_star[currentA]*a11[currentB] - a11[currentB]*a33[currentB] + a13[currentB]*a31[currentB])*yAB^2
					const3t = (3.0*a11_star[currentA]*a22_star[currentA] - 3.0*a12_star[currentA]*a21_star[currentA] - 2.0*a11_star[currentA]*a22[currentB] + 2.0*a12_star[currentA]*a21[currentB] + 2.0*a21_star[currentA]*a12[currentB] - 2.0*a22_star[currentA]*a11[currentB] + a11[currentB]*a22[currentB] - a12[currentB]*a21[currentB])*zAB^2 + (xAB*(3.0*a12_star[currentA]*a23_star[currentA] - 3.0*a13_star[currentA]*a22_star[currentA] + 3.0*a21_star[currentA]*a32_star[currentA] - 3.0*a22_star[currentA]*a31_star[currentA] - 2.0*a12_star[currentA]*a23[currentB] + 2.0*a13_star[currentA]*a22[currentB] + 2.0*a22_star[currentA]*a13[currentB] - 2.0*a23_star[currentA]*a12[currentB] - 2.0*a21_star[currentA]*a32[currentB] + 2.0*a22_star[currentA]*a31[currentB] + 2.0*a31_star[currentA]*a22[currentB] - 2.0*a32_star[currentA]*a21[currentB] + a12[currentB]*a23[currentB] - a13[currentB]*a22[currentB] + a21[currentB]*a32[currentB] - a22[currentB]*a31[currentB]) - yAB*(3.0*a11_star[currentA]*a23_star[currentA] - 3.0*a13_star[currentA]*a21_star[currentA] + 3.0*a11_star[currentA]*a32_star[currentA] - 3.0*a12_star[currentA]*a31_star[currentA] - 2.0*a11_star[currentA]*a23[currentB] + 2.0*a13_star[currentA]*a21[currentB] + 2.0*a21_star[currentA]*a13[currentB] - 2.0*a23_star[currentA]*a11[currentB] - 2.0*a11_star[currentA]*a32[currentB] + 2.0*a12_star[currentA]*a31[currentB] + 2.0*a31_star[currentA]*a12[currentB] - 2.0*a32_star[currentA]*a11[currentB] + a11[currentB]*a23[currentB] - a13[currentB]*a21[currentB] + a11[currentB]*a32[currentB] - a12[currentB]*a31[currentB]))*zAB + (3.0*a22_star[currentA]*a33_star[currentA] - 3.0*a23_star[currentA]*a32_star[currentA] - 2.0*a22_star[currentA]*a33[currentB] + 2.0*a23_star[currentA]*a32[currentB] + 2.0*a32_star[currentA]*a23[currentB] - 2.0*a33_star[currentA]*a22[currentB] + a22[currentB]*a33[currentB] - a23[currentB]*a32[currentB])*xAB^2 + (3.0*a13_star[currentA]*a32_star[currentA] - 3.0*a12_star[currentA]*a33_star[currentA] - 3.0*a21_star[currentA]*a33_star[currentA] + 3.0*a23_star[currentA]*a31_star[currentA] + 2.0*a12_star[currentA]*a33[currentB] - 2.0*a13_star[currentA]*a32[currentB] - 2.0*a32_star[currentA]*a13[currentB] + 2.0*a33_star[currentA]*a12[currentB] + 2.0*a21_star[currentA]*a33[currentB] - 2.0*a23_star[currentA]*a31[currentB] - 2.0*a31_star[currentA]*a23[currentB] + 2.0*a33_star[currentA]*a21[currentB] - a12[currentB]*a33[currentB] + a13[currentB]*a32[currentB] - a21[currentB]*a33[currentB] + a23[currentB]*a31[currentB])*xAB*yAB + (3.0*a11_star[currentA]*a33_star[currentA] - 3.0*a13_star[currentA]*a31_star[currentA] - 2.0*a11_star[currentA]*a33[currentB] + 2.0*a13_star[currentA]*a31[currentB] + 2.0*a31_star[currentA]*a13[currentB] - 2.0*a33_star[currentA]*a11[currentB] + a11[currentB]*a33[currentB] - a13[currentB]*a31[currentB])*yAB^2
					const2t = (3.0*a12_star[currentA]*a21_star[currentA] - 3.0*a11_star[currentA]*a22_star[currentA] + a11_star[currentA]*a22[currentB] - a12_star[currentA]*a21[currentB] - a21_star[currentA]*a12[currentB] + a22_star[currentA]*a11[currentB])*zAB^2 + (yAB*(3.0*a11_star[currentA]*a23_star[currentA] - 3.0*a13_star[currentA]*a21_star[currentA] + 3.0*a11_star[currentA]*a32_star[currentA] - 3.0*a12_star[currentA]*a31_star[currentA] - a11_star[currentA]*a23[currentB] + a13_star[currentA]*a21[currentB] + a21_star[currentA]*a13[currentB] - a23_star[currentA]*a11[currentB] - a11_star[currentA]*a32[currentB] + a12_star[currentA]*a31[currentB] + a31_star[currentA]*a12[currentB] - a32_star[currentA]*a11[currentB]) - xAB*(3.0*a12_star[currentA]*a23_star[currentA] - 3.0*a13_star[currentA]*a22_star[currentA] + 3.0*a21_star[currentA]*a32_star[currentA] - 3.0*a22_star[currentA]*a31_star[currentA] - a12_star[currentA]*a23[currentB] + a13_star[currentA]*a22[currentB] + a22_star[currentA]*a13[currentB] - a23_star[currentA]*a12[currentB] - a21_star[currentA]*a32[currentB] + a22_star[currentA]*a31[currentB] + a31_star[currentA]*a22[currentB] - a32_star[currentA]*a21[currentB]))*zAB + (3.0*a23_star[currentA]*a32_star[currentA] - 3.0*a22_star[currentA]*a33_star[currentA] + a22_star[currentA]*a33[currentB] - a23_star[currentA]*a32[currentB] - a32_star[currentA]*a23[currentB] + a33_star[currentA]*a22[currentB])*xAB^2 + (3.0*a12_star[currentA]*a33_star[currentA] - 3.0*a13_star[currentA]*a32_star[currentA] + 3.0*a21_star[currentA]*a33_star[currentA] - 3.0*a23_star[currentA]*a31_star[currentA] - a12_star[currentA]*a33[currentB] + a13_star[currentA]*a32[currentB] + a32_star[currentA]*a13[currentB] - a33_star[currentA]*a12[currentB] - a21_star[currentA]*a33[currentB] + a23_star[currentA]*a31[currentB] + a31_star[currentA]*a23[currentB] - a33_star[currentA]*a21[currentB])*xAB*yAB + (3.0*a13_star[currentA]*a31_star[currentA] - 3.0*a11_star[currentA]*a33_star[currentA] + a11_star[currentA]*a33[currentB] - a13_star[currentA]*a31[currentB] - a31_star[currentA]*a13[currentB] + a33_star[currentA]*a11[currentB])*yAB^2
					const1t = (a11_star[currentA]*a22_star[currentA] - a12_star[currentA]*a21_star[currentA])*zAB^2 + (xAB*(a12_star[currentA]*a23_star[currentA] - a13_star[currentA]*a22_star[currentA] + a21_star[currentA]*a32_star[currentA] - a22_star[currentA]*a31_star[currentA]) - yAB*(a11_star[currentA]*a23_star[currentA] - a13_star[currentA]*a21_star[currentA] + a11_star[currentA]*a32_star[currentA] - a12_star[currentA]*a31_star[currentA]))*zAB + (a22_star[currentA]*a33_star[currentA] - a23_star[currentA]*a32_star[currentA])*xAB^2 + (a13_star[currentA]*a32_star[currentA] - a12_star[currentA]*a33_star[currentA] - a21_star[currentA]*a33_star[currentA] + a23_star[currentA]*a31_star[currentA])*xAB*yAB + (a11_star[currentA]*a33_star[currentA] - a13_star[currentA]*a31_star[currentA])*yAB^2
					const3n = a11_star[currentA]*a23_star[currentA]*a32_star[currentA] - a11_star[currentA]*a22_star[currentA]*a33_star[currentA] + a12_star[currentA]*a21_star[currentA]*a33_star[currentA] - a12_star[currentA]*a23_star[currentA]*a31_star[currentA] - a13_star[currentA]*a21_star[currentA]*a32_star[currentA] + a13_star[currentA]*a22_star[currentA]*a31_star[currentA] + a11_star[currentA]*a22_star[currentA]*a33[currentB] - a11_star[currentA]*a23_star[currentA]*a32[currentB] - a11_star[currentA]*a32_star[currentA]*a23[currentB] + a11_star[currentA]*a33_star[currentA]*a22[currentB] - a12_star[currentA]*a21_star[currentA]*a33[currentB] + a12_star[currentA]*a23_star[currentA]*a31[currentB] + a12_star[currentA]*a31_star[currentA]*a23[currentB] - a12_star[currentA]*a33_star[currentA]*a21[currentB] + a13_star[currentA]*a21_star[currentA]*a32[currentB] - a13_star[currentA]*a22_star[currentA]*a31[currentB] - a13_star[currentA]*a31_star[currentA]*a22[currentB] + a13_star[currentA]*a32_star[currentA]*a21[currentB] + a21_star[currentA]*a32_star[currentA]*a13[currentB] - a21_star[currentA]*a33_star[currentA]*a12[currentB] - a22_star[currentA]*a31_star[currentA]*a13[currentB] + a22_star[currentA]*a33_star[currentA]*a11[currentB] + a23_star[currentA]*a31_star[currentA]*a12[currentB] - a23_star[currentA]*a32_star[currentA]*a11[currentB] - a11_star[currentA]*a22[currentB]*a33[currentB] + a11_star[currentA]*a23[currentB]*a32[currentB] + a12_star[currentA]*a21[currentB]*a33[currentB] - a12_star[currentA]*a23[currentB]*a31[currentB] - a13_star[currentA]*a21[currentB]*a32[currentB] + a13_star[currentA]*a22[currentB]*a31[currentB] + a21_star[currentA]*a12[currentB]*a33[currentB] - a21_star[currentA]*a13[currentB]*a32[currentB] - a22_star[currentA]*a11[currentB]*a33[currentB] + a22_star[currentA]*a13[currentB]*a31[currentB] + a23_star[currentA]*a11[currentB]*a32[currentB] - a23_star[currentA]*a12[currentB]*a31[currentB] - a31_star[currentA]*a12[currentB]*a23[currentB] + a31_star[currentA]*a13[currentB]*a22[currentB] + a32_star[currentA]*a11[currentB]*a23[currentB] - a32_star[currentA]*a13[currentB]*a21[currentB] - a33_star[currentA]*a11[currentB]*a22[currentB] + a33_star[currentA]*a12[currentB]*a21[currentB] + a11[currentB]*a22[currentB]*a33[currentB] - a11[currentB]*a23[currentB]*a32[currentB] - a12[currentB]*a21[currentB]*a33[currentB] + a12[currentB]*a23[currentB]*a31[currentB] + a13[currentB]*a21[currentB]*a32[currentB] - a13[currentB]*a22[currentB]*a31[currentB]
					const2n = 3.0*a11_star[currentA]*a22_star[currentA]*a33_star[currentA] - 3.0*a11_star[currentA]*a23_star[currentA]*a32_star[currentA] - 3.0*a12_star[currentA]*a21_star[currentA]*a33_star[currentA] + 3.0*a12_star[currentA]*a23_star[currentA]*a31_star[currentA] + 3.0*a13_star[currentA]*a21_star[currentA]*a32_star[currentA] - 3.0*a13_star[currentA]*a22_star[currentA]*a31_star[currentA] - 2.0*a11_star[currentA]*a22_star[currentA]*a33[currentB] + 2.0*a11_star[currentA]*a23_star[currentA]*a32[currentB] + 2.0*a11_star[currentA]*a32_star[currentA]*a23[currentB] - 2.0*a11_star[currentA]*a33_star[currentA]*a22[currentB] + 2.0*a12_star[currentA]*a21_star[currentA]*a33[currentB] - 2.0*a12_star[currentA]*a23_star[currentA]*a31[currentB] - 2.0*a12_star[currentA]*a31_star[currentA]*a23[currentB] + 2.0*a12_star[currentA]*a33_star[currentA]*a21[currentB] - 2.0*a13_star[currentA]*a21_star[currentA]*a32[currentB] + 2.0*a13_star[currentA]*a22_star[currentA]*a31[currentB] + 2.0*a13_star[currentA]*a31_star[currentA]*a22[currentB] - 2.0*a13_star[currentA]*a32_star[currentA]*a21[currentB] - 2.0*a21_star[currentA]*a32_star[currentA]*a13[currentB] + 2.0*a21_star[currentA]*a33_star[currentA]*a12[currentB] + 2.0*a22_star[currentA]*a31_star[currentA]*a13[currentB] - 2.0*a22_star[currentA]*a33_star[currentA]*a11[currentB] - 2.0*a23_star[currentA]*a31_star[currentA]*a12[currentB] + 2.0*a23_star[currentA]*a32_star[currentA]*a11[currentB] + a11_star[currentA]*a22[currentB]*a33[currentB] - a11_star[currentA]*a23[currentB]*a32[currentB] - a12_star[currentA]*a21[currentB]*a33[currentB] + a12_star[currentA]*a23[currentB]*a31[currentB] + a13_star[currentA]*a21[currentB]*a32[currentB] - a13_star[currentA]*a22[currentB]*a31[currentB] - a21_star[currentA]*a12[currentB]*a33[currentB] + a21_star[currentA]*a13[currentB]*a32[currentB] + a22_star[currentA]*a11[currentB]*a33[currentB] - a22_star[currentA]*a13[currentB]*a31[currentB] - a23_star[currentA]*a11[currentB]*a32[currentB] + a23_star[currentA]*a12[currentB]*a31[currentB] + a31_star[currentA]*a12[currentB]*a23[currentB] - a31_star[currentA]*a13[currentB]*a22[currentB] - a32_star[currentA]*a11[currentB]*a23[currentB] + a32_star[currentA]*a13[currentB]*a21[currentB] + a33_star[currentA]*a11[currentB]*a22[currentB] - a33_star[currentA]*a12[currentB]*a21[currentB]
					const1n = 3.0*a11_star[currentA]*a23_star[currentA]*a32_star[currentA] - 3.0*a11_star[currentA]*a22_star[currentA]*a33_star[currentA] + 3.0*a12_star[currentA]*a21_star[currentA]*a33_star[currentA] - 3.0*a12_star[currentA]*a23_star[currentA]*a31_star[currentA] - 3.0*a13_star[currentA]*a21_star[currentA]*a32_star[currentA] + 3.0*a13_star[currentA]*a22_star[currentA]*a31_star[currentA] + a11_star[currentA]*a22_star[currentA]*a33[currentB] - a11_star[currentA]*a23_star[currentA]*a32[currentB] - a11_star[currentA]*a32_star[currentA]*a23[currentB] + a11_star[currentA]*a33_star[currentA]*a22[currentB] - a12_star[currentA]*a21_star[currentA]*a33[currentB] + a12_star[currentA]*a23_star[currentA]*a31[currentB] + a12_star[currentA]*a31_star[currentA]*a23[currentB] - a12_star[currentA]*a33_star[currentA]*a21[currentB] + a13_star[currentA]*a21_star[currentA]*a32[currentB] - a13_star[currentA]*a22_star[currentA]*a31[currentB] - a13_star[currentA]*a31_star[currentA]*a22[currentB] + a13_star[currentA]*a32_star[currentA]*a21[currentB] + a21_star[currentA]*a32_star[currentA]*a13[currentB] - a21_star[currentA]*a33_star[currentA]*a12[currentB] - a22_star[currentA]*a31_star[currentA]*a13[currentB] + a22_star[currentA]*a33_star[currentA]*a11[currentB] + a23_star[currentA]*a31_star[currentA]*a12[currentB] - a23_star[currentA]*a32_star[currentA]*a11[currentB]

					b = 0.5

					fb = (((const4t*b+const3t)*b+const2t)*b+const1t)/((const3n*b+const2n)*b+const1n)
					if fb >= 1.0
						overlapfun = 1.0
					else
						iter = 1
						a = 0.0
						c = 1.0
						while fb < 1 && iter < number_of_iterations_overlap_criterion
							# Pick point in upper half.
							x = 0.5*(b+c)
							fx = (((const4t*x+const3t)*x+const2t)*x+const1t)/((const3n*x+const2n)*x+const1n)
							if fx < fb
								c = x
							else
								b = x
								fb = fx
							end
							# Pick point in lower half.
							x = 0.5*(a+b)
							fx = (((const4t*x+const3t)*x+const2t)*x+const1t)/((const3n*x+const2n)*x+const1n)

							if fx < fb
								a = x
							else
								b = x
								fb = fx
							end

							iter = iter + 1
						end
						overlapfun = fb
					end
					if overlapfun < 1.0
						energy_particle_star += (1.0-overlapfun)^2
					end
				end
			end

			angle_to_z_axis = acos((cos(THETA1_star[currentA])*cos(THETA2_star[currentA]))/(abs(sin(THETA2_star[currentA]))^2 + abs(cos(THETA1_star[currentA])*cos(THETA2_star[currentA]))^2 + abs(cos(THETA2_star[currentA])*sin(THETA1_star[currentA]))^2)^(1/2))
			#if angle_to_z_axis > ubangle
			#	energy_particle_star = energy_particle + 1.0 # Just a value that prevents the rotation.
			#end

			if energy_particle_star <= energy_particle && angle_to_z_axis <= ubangle
				THETA1[currentA] = THETA1_star[currentA]
				THETA2[currentA] = THETA2_star[currentA]
				THETA3[currentA] = THETA3_star[currentA]

				cA1[currentA] = cA1_star[currentA]
				cA2[currentA] = cA2_star[currentA]
				sA1[currentA] = sA1_star[currentA]
				sA2[currentA] = sA2_star[currentA]
				cA1d[currentA] = cA1d_star[currentA]
				cA2d[currentA] = cA2d_star[currentA]
				cA3d[currentA] = cA3d_star[currentA]
				sA1d[currentA] = sA1d_star[currentA]
				sA2d[currentA] = sA2d_star[currentA]
				sA3d[currentA] = sA3d_star[currentA]
				a11[currentA] = a11_star[currentA]
				a12[currentA] = a12_star[currentA]
				a13[currentA] = a13_star[currentA]
				a21[currentA] = a21_star[currentA]
				a22[currentA] = a22_star[currentA]
				a23[currentA] = a23_star[currentA]
				a31[currentA] = a31_star[currentA]
				a32[currentA] = a32_star[currentA]
				a33[currentA] = a33_star[currentA]

				acceptance_probability_rotation += 1.0
				energy_particle = energy_particle_star
			end

			energy_system += energy_particle
		end
		acceptance_probability_translation /= number_of_particles
		acceptance_probability_rotation /= number_of_particles

		if acceptance_probability_translation <= acceptance_probability_target
			sigma_translation *= 0.95
		else
			sigma_translation *= 1.05
			sigma_translation = min(sigma_translation, sigma_translation_max)
		end

		if acceptance_probability_rotation <= acceptance_probability_target
			sigma_rotation *= 0.95
		else
			sigma_rotation *= 1.05
			sigma_rotation = min(sigma_rotation, sigma_rotation_max)
		end

		m1x = 0.0
		m2x = 0.0
		m1y = 0.0
		m2y = 0.0
		m1z = 0.0
		m2z = 0.0
		for currentA = 1:number_of_particles
			angle_to_x_axis = acos((cos(THETA2[currentA])*cos(THETA3[currentA]))/(abs(sin(THETA3[currentA]))^2 + abs(cos(THETA2[currentA])*cos(THETA3[currentA]))^2 + abs(cos(THETA3[currentA])*sin(THETA2[currentA]))^2)^(1/2))
			angle_to_y_axis = acos((cos(THETA3[currentA])*cos(THETA1[currentA]))/(abs(sin(THETA1[currentA]))^2 + abs(cos(THETA3[currentA])*cos(THETA1[currentA]))^2 + abs(cos(THETA1[currentA])*sin(THETA3[currentA]))^2)^(1/2))
			angle_to_z_axis = acos((cos(THETA1[currentA])*cos(THETA2[currentA]))/(abs(sin(THETA2[currentA]))^2 + abs(cos(THETA1[currentA])*cos(THETA2[currentA]))^2 + abs(cos(THETA2[currentA])*sin(THETA1[currentA]))^2)^(1/2))

			m1x += angle_to_x_axis
			m2x += angle_to_x_axis^2
			m1y += angle_to_y_axis
			m2y += angle_to_y_axis^2
			m1z += angle_to_z_axis
			m2z += angle_to_z_axis^2

		end
		m1x /= number_of_particles
		m2x /= number_of_particles
		m1y /= number_of_particles
		m2y /= number_of_particles
		m1z /= number_of_particles
		m2z /= number_of_particles
		sx = sqrt(m2x - m1x^2)
		sy = sqrt(m2y - m1y^2)
		sz = sqrt(m2z - m1z^2)
		println((sx, sy, sz))
		#max_val = 0.0
		#for currentA = 1:number_of_particles
		#	angle_to_z_axis = acos((cos(THETA1[currentA])*cos(THETA2[currentA]))/(abs(sin(THETA2[currentA]))^2 + abs(cos(THETA1[currentA])*cos(THETA2[currentA]))^2 + abs(cos(THETA2[currentA])*sin(THETA1[currentA]))^2)^(1/2))
		#	max_val = max(max_val, angle_to_z_axis)
		#end
		#println(max_val)
	end

	# Equilibrate configuration.
	if !silent
		println("Equilibrating configuration (performing random moves)...")
	end
	for current_sweep = 1:number_of_equilibration_sweeps
		if !silent
			println(join(["   Sweep ", string(current_sweep), " of ", string(number_of_equilibration_sweeps)]))
		end
		randn!(X_star)
		randn!(Y_star)
		randn!(Z_star)
		randn!(THETA1_star)
		randn!(THETA2_star)
		randn!(THETA3_star)
		X_star = X + sigma_translation * X_star
		Y_star = Y + sigma_translation * Y_star
		Z_star = Z + sigma_translation * Z_star
		for current_particle = 1:number_of_particles
			if X_star[current_particle] < 0.0
				X_star[current_particle] += Lx
			elseif X_star[current_particle] > Lx
				X_star[current_particle] -= Lx
			end
			if Y_star[current_particle] < 0.0
				Y_star[current_particle] += Ly
			elseif Y_star[current_particle] > Ly
				Y_star[current_particle] -= Ly
			end
			if Z_star[current_particle] < 0.0
				Z_star[current_particle] += Lz
			elseif Z_star[current_particle] > Lz
				Z_star[current_particle] -= Lz
			end
		end
		THETA1_star = THETA1 + sigma_rotation * THETA1_star
		THETA2_star = THETA2 + sigma_rotation * THETA2_star
		THETA3_star = THETA3 + sigma_rotation * THETA3_star

		cA1_star = cos(THETA1_star)
		cA2_star = cos(THETA2_star)
		sA1_star = sin(THETA1_star)
		sA2_star = sin(THETA2_star)
		cA1d_star = cos(2.0*THETA1_star)
		cA2d_star = cos(2.0*THETA2_star)
		cA3d_star = cos(2.0*THETA3_star)
		sA1d_star = sin(2.0*THETA1_star)
		sA2d_star = sin(2.0*THETA2_star)
		sA3d_star = sin(2.0*THETA3_star)
		a11_star = (cA2d_star*0.25 + cA3d_star*0.25 + (cA2d_star.*cA3d_star)*0.25 + 0.25).*R1.^2 + (cA2d_star*0.25 - cA3d_star*0.25 - (cA2d_star.*cA3d_star)*0.25 + 0.25).*R2.^2
		a12_star = ((sA1_star.*sA2d_star)*0.25 + (cA1_star.*cA2_star.*sA3d_star)*0.5 + (cA3d_star.*sA1_star.*sA2d_star)*0.25).*R1.^2 + ((sA1_star.*sA2d_star)*0.25 - (cA1_star.*cA2_star.*sA3d_star)*0.5 - (cA3d_star.*sA1_star.*sA2d_star)*0.25).*R2.^2
		a13_star = ((cA2_star.*sA1_star.*sA3d_star)*0.5 - (cA1_star.*cA3d_star.*sA2d_star)*0.25 - (cA1_star.*sA2d_star)*0.25).*R1.^2 + ((cA1_star.*cA3d_star.*sA2d_star)*0.25 - (cA1_star.*sA2d_star)*0.25 - (cA2_star.*sA1_star.*sA3d_star)*0.5).*R2.^2
		a21_star = a12_star
		a22_star = ((cA2_star.^2.*cA1d_star)*0.25 - (cA1d_star.*cA3d_star)*0.5 - (cA2_star.^2.*cA3d_star)*0.25 - cA2_star.^2.0*0.25 + (sA2_star.*sA1d_star.*sA3d_star)*0.5 + (cA2_star.^2.*cA1d_star.*cA3d_star)*0.25 + 0.5).*R1.^2 + ((cA1d_star.*cA3d_star)*0.5 + (cA2_star.^2.*cA1d_star)*0.25 + (cA2_star.^2.*cA3d_star)*0.25 - cA2_star.^2.0*0.25 - (sA2_star.*sA1d_star.*sA3d_star)*0.5 - (cA2_star.^2.*cA1d_star.*cA3d_star)*0.25 + 0.5).*R2.^2
		a23_star = ((cA2_star.^2.*sA1d_star)*0.25 - (cA3d_star.*sA1d_star)*0.5 - (cA1d_star.*sA2_star.*sA3d_star)*0.5 + (cA2_star.^2.*cA3d_star.*sA1d_star)*0.25).*R1.^2 + ((cA3d_star.*sA1d_star)*0.5 + (cA2_star.^2.*sA1d_star)*0.25 + (cA1d_star.*sA2_star.*sA3d_star)*0.5 - (cA2_star.^2.*cA3d_star.*sA1d_star)*0.25).*R2.^2
		a31_star = a13_star
		a32_star = a23_star
		a33_star = ((cA1d_star.*cA3d_star)*0.5 - (cA2_star.^2.*cA1d_star)*0.25 - (cA2_star.^2.*cA3d_star)*0.25 - cA2_star.^2.0*0.25 - (sA2_star.*sA1d_star.*sA3d_star)*0.5 - (cA2_star.^2.*cA1d_star.*cA3d_star)*0.25 + 0.5).*R1.^2 + ((cA2_star.^2.*cA3d_star)*0.25 - (cA2_star.^2.*cA1d_star)*0.25 - (cA1d_star.*cA3d_star)*0.5 - cA2_star.^2.0*0.25 + (sA2_star.*sA1d_star.*sA3d_star)*0.5 + (cA2_star.^2.*cA1d_star.*cA3d_star)*0.25 + 0.5).*R2.^2

		acceptance_probability_translation = 0.0
		acceptance_probability_rotation = 0.0
		for currentA = 1:number_of_particles
			# Compute current local energy.
			energy_particle = 0.0
			for currentB = [1:currentA-1;currentA+1:number_of_particles]
				xAB = X[currentB] - X[currentA]
				if xAB < -0.5*Lx
					xAB += Lx
				elseif xAB > 0.5*Lx
					xAB -= Lx
				end

				yAB = Y[currentB] - Y[currentA]
				if yAB < -0.5*Ly
					yAB += Ly
				elseif yAB > 0.5*Ly
					yAB -= Ly
				end

				zAB = Z[currentB] - Z[currentA]
				if zAB < -0.5*Lz
					zAB += Lz
				elseif zAB > 0.5*Lz
					zAB -= Lz
				end

				if xAB^2 + yAB^2 + zAB^2 < (RMAX[currentA] + RMAX[currentB])^2
					const4t = (a12[currentA]*a21[currentA] - a11[currentA]*a22[currentA] + a11[currentA]*a22[currentB] - a12[currentA]*a21[currentB] - a21[currentA]*a12[currentB] + a22[currentA]*a11[currentB] - a11[currentB]*a22[currentB] + a12[currentB]*a21[currentB])*zAB^2 + (yAB*(a11[currentA]*a23[currentA] - a13[currentA]*a21[currentA] + a11[currentA]*a32[currentA] - a12[currentA]*a31[currentA] - a11[currentA]*a23[currentB] + a13[currentA]*a21[currentB] + a21[currentA]*a13[currentB] - a23[currentA]*a11[currentB] - a11[currentA]*a32[currentB] + a12[currentA]*a31[currentB] + a31[currentA]*a12[currentB] - a32[currentA]*a11[currentB] + a11[currentB]*a23[currentB] - a13[currentB]*a21[currentB] + a11[currentB]*a32[currentB] - a12[currentB]*a31[currentB]) - xAB*(a12[currentA]*a23[currentA] - a13[currentA]*a22[currentA] + a21[currentA]*a32[currentA] - a22[currentA]*a31[currentA] - a12[currentA]*a23[currentB] + a13[currentA]*a22[currentB] + a22[currentA]*a13[currentB] - a23[currentA]*a12[currentB] - a21[currentA]*a32[currentB] + a22[currentA]*a31[currentB] + a31[currentA]*a22[currentB] - a32[currentA]*a21[currentB] + a12[currentB]*a23[currentB] - a13[currentB]*a22[currentB] + a21[currentB]*a32[currentB] - a22[currentB]*a31[currentB]))*zAB + (a23[currentA]*a32[currentA] - a22[currentA]*a33[currentA] + a22[currentA]*a33[currentB] - a23[currentA]*a32[currentB] - a32[currentA]*a23[currentB] + a33[currentA]*a22[currentB] - a22[currentB]*a33[currentB] + a23[currentB]*a32[currentB])*xAB^2 + (a12[currentA]*a33[currentA] - a13[currentA]*a32[currentA] + a21[currentA]*a33[currentA] - a23[currentA]*a31[currentA] - a12[currentA]*a33[currentB] + a13[currentA]*a32[currentB] + a32[currentA]*a13[currentB] - a33[currentA]*a12[currentB] - a21[currentA]*a33[currentB] + a23[currentA]*a31[currentB] + a31[currentA]*a23[currentB] - a33[currentA]*a21[currentB] + a12[currentB]*a33[currentB] - a13[currentB]*a32[currentB] + a21[currentB]*a33[currentB] - a23[currentB]*a31[currentB])*xAB*yAB + (a13[currentA]*a31[currentA] - a11[currentA]*a33[currentA] + a11[currentA]*a33[currentB] - a13[currentA]*a31[currentB] - a31[currentA]*a13[currentB] + a33[currentA]*a11[currentB] - a11[currentB]*a33[currentB] + a13[currentB]*a31[currentB])*yAB^2
					const3t = (3.0*a11[currentA]*a22[currentA] - 3.0*a12[currentA]*a21[currentA] - 2.0*a11[currentA]*a22[currentB] + 2.0*a12[currentA]*a21[currentB] + 2.0*a21[currentA]*a12[currentB] - 2.0*a22[currentA]*a11[currentB] + a11[currentB]*a22[currentB] - a12[currentB]*a21[currentB])*zAB^2 + (xAB*(3.0*a12[currentA]*a23[currentA] - 3.0*a13[currentA]*a22[currentA] + 3.0*a21[currentA]*a32[currentA] - 3.0*a22[currentA]*a31[currentA] - 2.0*a12[currentA]*a23[currentB] + 2.0*a13[currentA]*a22[currentB] + 2.0*a22[currentA]*a13[currentB] - 2.0*a23[currentA]*a12[currentB] - 2.0*a21[currentA]*a32[currentB] + 2.0*a22[currentA]*a31[currentB] + 2.0*a31[currentA]*a22[currentB] - 2.0*a32[currentA]*a21[currentB] + a12[currentB]*a23[currentB] - a13[currentB]*a22[currentB] + a21[currentB]*a32[currentB] - a22[currentB]*a31[currentB]) - yAB*(3.0*a11[currentA]*a23[currentA] - 3.0*a13[currentA]*a21[currentA] + 3.0*a11[currentA]*a32[currentA] - 3.0*a12[currentA]*a31[currentA] - 2.0*a11[currentA]*a23[currentB] + 2.0*a13[currentA]*a21[currentB] + 2.0*a21[currentA]*a13[currentB] - 2.0*a23[currentA]*a11[currentB] - 2.0*a11[currentA]*a32[currentB] + 2.0*a12[currentA]*a31[currentB] + 2.0*a31[currentA]*a12[currentB] - 2.0*a32[currentA]*a11[currentB] + a11[currentB]*a23[currentB] - a13[currentB]*a21[currentB] + a11[currentB]*a32[currentB] - a12[currentB]*a31[currentB]))*zAB + (3.0*a22[currentA]*a33[currentA] - 3.0*a23[currentA]*a32[currentA] - 2.0*a22[currentA]*a33[currentB] + 2.0*a23[currentA]*a32[currentB] + 2.0*a32[currentA]*a23[currentB] - 2.0*a33[currentA]*a22[currentB] + a22[currentB]*a33[currentB] - a23[currentB]*a32[currentB])*xAB^2 + (3.0*a13[currentA]*a32[currentA] - 3.0*a12[currentA]*a33[currentA] - 3.0*a21[currentA]*a33[currentA] + 3.0*a23[currentA]*a31[currentA] + 2.0*a12[currentA]*a33[currentB] - 2.0*a13[currentA]*a32[currentB] - 2.0*a32[currentA]*a13[currentB] + 2.0*a33[currentA]*a12[currentB] + 2.0*a21[currentA]*a33[currentB] - 2.0*a23[currentA]*a31[currentB] - 2.0*a31[currentA]*a23[currentB] + 2.0*a33[currentA]*a21[currentB] - a12[currentB]*a33[currentB] + a13[currentB]*a32[currentB] - a21[currentB]*a33[currentB] + a23[currentB]*a31[currentB])*xAB*yAB + (3.0*a11[currentA]*a33[currentA] - 3.0*a13[currentA]*a31[currentA] - 2.0*a11[currentA]*a33[currentB] + 2.0*a13[currentA]*a31[currentB] + 2.0*a31[currentA]*a13[currentB] - 2.0*a33[currentA]*a11[currentB] + a11[currentB]*a33[currentB] - a13[currentB]*a31[currentB])*yAB^2
					const2t = (3.0*a12[currentA]*a21[currentA] - 3.0*a11[currentA]*a22[currentA] + a11[currentA]*a22[currentB] - a12[currentA]*a21[currentB] - a21[currentA]*a12[currentB] + a22[currentA]*a11[currentB])*zAB^2 + (yAB*(3.0*a11[currentA]*a23[currentA] - 3.0*a13[currentA]*a21[currentA] + 3.0*a11[currentA]*a32[currentA] - 3.0*a12[currentA]*a31[currentA] - a11[currentA]*a23[currentB] + a13[currentA]*a21[currentB] + a21[currentA]*a13[currentB] - a23[currentA]*a11[currentB] - a11[currentA]*a32[currentB] + a12[currentA]*a31[currentB] + a31[currentA]*a12[currentB] - a32[currentA]*a11[currentB]) - xAB*(3.0*a12[currentA]*a23[currentA] - 3.0*a13[currentA]*a22[currentA] + 3.0*a21[currentA]*a32[currentA] - 3.0*a22[currentA]*a31[currentA] - a12[currentA]*a23[currentB] + a13[currentA]*a22[currentB] + a22[currentA]*a13[currentB] - a23[currentA]*a12[currentB] - a21[currentA]*a32[currentB] + a22[currentA]*a31[currentB] + a31[currentA]*a22[currentB] - a32[currentA]*a21[currentB]))*zAB + (3.0*a23[currentA]*a32[currentA] - 3.0*a22[currentA]*a33[currentA] + a22[currentA]*a33[currentB] - a23[currentA]*a32[currentB] - a32[currentA]*a23[currentB] + a33[currentA]*a22[currentB])*xAB^2 + (3.0*a12[currentA]*a33[currentA] - 3.0*a13[currentA]*a32[currentA] + 3.0*a21[currentA]*a33[currentA] - 3.0*a23[currentA]*a31[currentA] - a12[currentA]*a33[currentB] + a13[currentA]*a32[currentB] + a32[currentA]*a13[currentB] - a33[currentA]*a12[currentB] - a21[currentA]*a33[currentB] + a23[currentA]*a31[currentB] + a31[currentA]*a23[currentB] - a33[currentA]*a21[currentB])*xAB*yAB + (3.0*a13[currentA]*a31[currentA] - 3.0*a11[currentA]*a33[currentA] + a11[currentA]*a33[currentB] - a13[currentA]*a31[currentB] - a31[currentA]*a13[currentB] + a33[currentA]*a11[currentB])*yAB^2
					const1t = (a11[currentA]*a22[currentA] - a12[currentA]*a21[currentA])*zAB^2 + (xAB*(a12[currentA]*a23[currentA] - a13[currentA]*a22[currentA] + a21[currentA]*a32[currentA] - a22[currentA]*a31[currentA]) - yAB*(a11[currentA]*a23[currentA] - a13[currentA]*a21[currentA] + a11[currentA]*a32[currentA] - a12[currentA]*a31[currentA]))*zAB + (a22[currentA]*a33[currentA] - a23[currentA]*a32[currentA])*xAB^2 + (a13[currentA]*a32[currentA] - a12[currentA]*a33[currentA] - a21[currentA]*a33[currentA] + a23[currentA]*a31[currentA])*xAB*yAB + (a11[currentA]*a33[currentA] - a13[currentA]*a31[currentA])*yAB^2
					const3n = a11[currentA]*a23[currentA]*a32[currentA] - a11[currentA]*a22[currentA]*a33[currentA] + a12[currentA]*a21[currentA]*a33[currentA] - a12[currentA]*a23[currentA]*a31[currentA] - a13[currentA]*a21[currentA]*a32[currentA] + a13[currentA]*a22[currentA]*a31[currentA] + a11[currentA]*a22[currentA]*a33[currentB] - a11[currentA]*a23[currentA]*a32[currentB] - a11[currentA]*a32[currentA]*a23[currentB] + a11[currentA]*a33[currentA]*a22[currentB] - a12[currentA]*a21[currentA]*a33[currentB] + a12[currentA]*a23[currentA]*a31[currentB] + a12[currentA]*a31[currentA]*a23[currentB] - a12[currentA]*a33[currentA]*a21[currentB] + a13[currentA]*a21[currentA]*a32[currentB] - a13[currentA]*a22[currentA]*a31[currentB] - a13[currentA]*a31[currentA]*a22[currentB] + a13[currentA]*a32[currentA]*a21[currentB] + a21[currentA]*a32[currentA]*a13[currentB] - a21[currentA]*a33[currentA]*a12[currentB] - a22[currentA]*a31[currentA]*a13[currentB] + a22[currentA]*a33[currentA]*a11[currentB] + a23[currentA]*a31[currentA]*a12[currentB] - a23[currentA]*a32[currentA]*a11[currentB] - a11[currentA]*a22[currentB]*a33[currentB] + a11[currentA]*a23[currentB]*a32[currentB] + a12[currentA]*a21[currentB]*a33[currentB] - a12[currentA]*a23[currentB]*a31[currentB] - a13[currentA]*a21[currentB]*a32[currentB] + a13[currentA]*a22[currentB]*a31[currentB] + a21[currentA]*a12[currentB]*a33[currentB] - a21[currentA]*a13[currentB]*a32[currentB] - a22[currentA]*a11[currentB]*a33[currentB] + a22[currentA]*a13[currentB]*a31[currentB] + a23[currentA]*a11[currentB]*a32[currentB] - a23[currentA]*a12[currentB]*a31[currentB] - a31[currentA]*a12[currentB]*a23[currentB] + a31[currentA]*a13[currentB]*a22[currentB] + a32[currentA]*a11[currentB]*a23[currentB] - a32[currentA]*a13[currentB]*a21[currentB] - a33[currentA]*a11[currentB]*a22[currentB] + a33[currentA]*a12[currentB]*a21[currentB] + a11[currentB]*a22[currentB]*a33[currentB] - a11[currentB]*a23[currentB]*a32[currentB] - a12[currentB]*a21[currentB]*a33[currentB] + a12[currentB]*a23[currentB]*a31[currentB] + a13[currentB]*a21[currentB]*a32[currentB] - a13[currentB]*a22[currentB]*a31[currentB]
					const2n = 3.0*a11[currentA]*a22[currentA]*a33[currentA] - 3.0*a11[currentA]*a23[currentA]*a32[currentA] - 3.0*a12[currentA]*a21[currentA]*a33[currentA] + 3.0*a12[currentA]*a23[currentA]*a31[currentA] + 3.0*a13[currentA]*a21[currentA]*a32[currentA] - 3.0*a13[currentA]*a22[currentA]*a31[currentA] - 2.0*a11[currentA]*a22[currentA]*a33[currentB] + 2.0*a11[currentA]*a23[currentA]*a32[currentB] + 2.0*a11[currentA]*a32[currentA]*a23[currentB] - 2.0*a11[currentA]*a33[currentA]*a22[currentB] + 2.0*a12[currentA]*a21[currentA]*a33[currentB] - 2.0*a12[currentA]*a23[currentA]*a31[currentB] - 2.0*a12[currentA]*a31[currentA]*a23[currentB] + 2.0*a12[currentA]*a33[currentA]*a21[currentB] - 2.0*a13[currentA]*a21[currentA]*a32[currentB] + 2.0*a13[currentA]*a22[currentA]*a31[currentB] + 2.0*a13[currentA]*a31[currentA]*a22[currentB] - 2.0*a13[currentA]*a32[currentA]*a21[currentB] - 2.0*a21[currentA]*a32[currentA]*a13[currentB] + 2.0*a21[currentA]*a33[currentA]*a12[currentB] + 2.0*a22[currentA]*a31[currentA]*a13[currentB] - 2.0*a22[currentA]*a33[currentA]*a11[currentB] - 2.0*a23[currentA]*a31[currentA]*a12[currentB] + 2.0*a23[currentA]*a32[currentA]*a11[currentB] + a11[currentA]*a22[currentB]*a33[currentB] - a11[currentA]*a23[currentB]*a32[currentB] - a12[currentA]*a21[currentB]*a33[currentB] + a12[currentA]*a23[currentB]*a31[currentB] + a13[currentA]*a21[currentB]*a32[currentB] - a13[currentA]*a22[currentB]*a31[currentB] - a21[currentA]*a12[currentB]*a33[currentB] + a21[currentA]*a13[currentB]*a32[currentB] + a22[currentA]*a11[currentB]*a33[currentB] - a22[currentA]*a13[currentB]*a31[currentB] - a23[currentA]*a11[currentB]*a32[currentB] + a23[currentA]*a12[currentB]*a31[currentB] + a31[currentA]*a12[currentB]*a23[currentB] - a31[currentA]*a13[currentB]*a22[currentB] - a32[currentA]*a11[currentB]*a23[currentB] + a32[currentA]*a13[currentB]*a21[currentB] + a33[currentA]*a11[currentB]*a22[currentB] - a33[currentA]*a12[currentB]*a21[currentB]
					const1n = 3.0*a11[currentA]*a23[currentA]*a32[currentA] - 3.0*a11[currentA]*a22[currentA]*a33[currentA] + 3.0*a12[currentA]*a21[currentA]*a33[currentA] - 3.0*a12[currentA]*a23[currentA]*a31[currentA] - 3.0*a13[currentA]*a21[currentA]*a32[currentA] + 3.0*a13[currentA]*a22[currentA]*a31[currentA] + a11[currentA]*a22[currentA]*a33[currentB] - a11[currentA]*a23[currentA]*a32[currentB] - a11[currentA]*a32[currentA]*a23[currentB] + a11[currentA]*a33[currentA]*a22[currentB] - a12[currentA]*a21[currentA]*a33[currentB] + a12[currentA]*a23[currentA]*a31[currentB] + a12[currentA]*a31[currentA]*a23[currentB] - a12[currentA]*a33[currentA]*a21[currentB] + a13[currentA]*a21[currentA]*a32[currentB] - a13[currentA]*a22[currentA]*a31[currentB] - a13[currentA]*a31[currentA]*a22[currentB] + a13[currentA]*a32[currentA]*a21[currentB] + a21[currentA]*a32[currentA]*a13[currentB] - a21[currentA]*a33[currentA]*a12[currentB] - a22[currentA]*a31[currentA]*a13[currentB] + a22[currentA]*a33[currentA]*a11[currentB] + a23[currentA]*a31[currentA]*a12[currentB] - a23[currentA]*a32[currentA]*a11[currentB]

					b = 0.5

					fb = (((const4t*b+const3t)*b+const2t)*b+const1t)/((const3n*b+const2n)*b+const1n)
					if fb >= 1.0
						overlapfun = 1.0
					else
						iter = 1
						a = 0.0
						c = 1.0
						while fb < 1 && iter < number_of_iterations_overlap_criterion
							# Pick point in upper half.
							x = 0.5*(b+c)
							fx = (((const4t*x+const3t)*x+const2t)*x+const1t)/((const3n*x+const2n)*x+const1n)
							if fx < fb
								c = x
							else
								b = x
								fb = fx
							end
							# Pick point in lower half.
							x = 0.5*(a+b)
							fx = (((const4t*x+const3t)*x+const2t)*x+const1t)/((const3n*x+const2n)*x+const1n)

							if fx < fb
								a = x
							else
								b = x
								fb = fx
							end

							iter = iter + 1
						end
						overlapfun = fb
					end
					if overlapfun < 1.0
						energy_particle += (1.0-overlapfun)^2
					end
				end
			end

			# Compute new local energy with translation.
			energy_particle_star = 0.0
			for currentB = [1:currentA-1;currentA+1:number_of_particles]
				xAB = X[currentB] - X_star[currentA]
				if xAB < -0.5*Lx
					xAB += Lx
				elseif xAB > 0.5*Lx
					xAB -= Lx
				end

				yAB = Y[currentB] - Y_star[currentA]
				if yAB < -0.5*Ly
					yAB += Ly
				elseif yAB > 0.5*Ly
					yAB -= Ly
				end

				zAB = Z[currentB] - Z_star[currentA]
				if zAB < -0.5*Lz
					zAB += Lz
				elseif zAB > 0.5*Lz
					zAB -= Lz
				end

				if xAB^2 + yAB^2 + zAB^2 < (RMAX[currentA] + RMAX[currentB])^2
					const4t = (a12[currentA]*a21[currentA] - a11[currentA]*a22[currentA] + a11[currentA]*a22[currentB] - a12[currentA]*a21[currentB] - a21[currentA]*a12[currentB] + a22[currentA]*a11[currentB] - a11[currentB]*a22[currentB] + a12[currentB]*a21[currentB])*zAB^2 + (yAB*(a11[currentA]*a23[currentA] - a13[currentA]*a21[currentA] + a11[currentA]*a32[currentA] - a12[currentA]*a31[currentA] - a11[currentA]*a23[currentB] + a13[currentA]*a21[currentB] + a21[currentA]*a13[currentB] - a23[currentA]*a11[currentB] - a11[currentA]*a32[currentB] + a12[currentA]*a31[currentB] + a31[currentA]*a12[currentB] - a32[currentA]*a11[currentB] + a11[currentB]*a23[currentB] - a13[currentB]*a21[currentB] + a11[currentB]*a32[currentB] - a12[currentB]*a31[currentB]) - xAB*(a12[currentA]*a23[currentA] - a13[currentA]*a22[currentA] + a21[currentA]*a32[currentA] - a22[currentA]*a31[currentA] - a12[currentA]*a23[currentB] + a13[currentA]*a22[currentB] + a22[currentA]*a13[currentB] - a23[currentA]*a12[currentB] - a21[currentA]*a32[currentB] + a22[currentA]*a31[currentB] + a31[currentA]*a22[currentB] - a32[currentA]*a21[currentB] + a12[currentB]*a23[currentB] - a13[currentB]*a22[currentB] + a21[currentB]*a32[currentB] - a22[currentB]*a31[currentB]))*zAB + (a23[currentA]*a32[currentA] - a22[currentA]*a33[currentA] + a22[currentA]*a33[currentB] - a23[currentA]*a32[currentB] - a32[currentA]*a23[currentB] + a33[currentA]*a22[currentB] - a22[currentB]*a33[currentB] + a23[currentB]*a32[currentB])*xAB^2 + (a12[currentA]*a33[currentA] - a13[currentA]*a32[currentA] + a21[currentA]*a33[currentA] - a23[currentA]*a31[currentA] - a12[currentA]*a33[currentB] + a13[currentA]*a32[currentB] + a32[currentA]*a13[currentB] - a33[currentA]*a12[currentB] - a21[currentA]*a33[currentB] + a23[currentA]*a31[currentB] + a31[currentA]*a23[currentB] - a33[currentA]*a21[currentB] + a12[currentB]*a33[currentB] - a13[currentB]*a32[currentB] + a21[currentB]*a33[currentB] - a23[currentB]*a31[currentB])*xAB*yAB + (a13[currentA]*a31[currentA] - a11[currentA]*a33[currentA] + a11[currentA]*a33[currentB] - a13[currentA]*a31[currentB] - a31[currentA]*a13[currentB] + a33[currentA]*a11[currentB] - a11[currentB]*a33[currentB] + a13[currentB]*a31[currentB])*yAB^2
					const3t = (3.0*a11[currentA]*a22[currentA] - 3.0*a12[currentA]*a21[currentA] - 2.0*a11[currentA]*a22[currentB] + 2.0*a12[currentA]*a21[currentB] + 2.0*a21[currentA]*a12[currentB] - 2.0*a22[currentA]*a11[currentB] + a11[currentB]*a22[currentB] - a12[currentB]*a21[currentB])*zAB^2 + (xAB*(3.0*a12[currentA]*a23[currentA] - 3.0*a13[currentA]*a22[currentA] + 3.0*a21[currentA]*a32[currentA] - 3.0*a22[currentA]*a31[currentA] - 2.0*a12[currentA]*a23[currentB] + 2.0*a13[currentA]*a22[currentB] + 2.0*a22[currentA]*a13[currentB] - 2.0*a23[currentA]*a12[currentB] - 2.0*a21[currentA]*a32[currentB] + 2.0*a22[currentA]*a31[currentB] + 2.0*a31[currentA]*a22[currentB] - 2.0*a32[currentA]*a21[currentB] + a12[currentB]*a23[currentB] - a13[currentB]*a22[currentB] + a21[currentB]*a32[currentB] - a22[currentB]*a31[currentB]) - yAB*(3.0*a11[currentA]*a23[currentA] - 3.0*a13[currentA]*a21[currentA] + 3.0*a11[currentA]*a32[currentA] - 3.0*a12[currentA]*a31[currentA] - 2.0*a11[currentA]*a23[currentB] + 2.0*a13[currentA]*a21[currentB] + 2.0*a21[currentA]*a13[currentB] - 2.0*a23[currentA]*a11[currentB] - 2.0*a11[currentA]*a32[currentB] + 2.0*a12[currentA]*a31[currentB] + 2.0*a31[currentA]*a12[currentB] - 2.0*a32[currentA]*a11[currentB] + a11[currentB]*a23[currentB] - a13[currentB]*a21[currentB] + a11[currentB]*a32[currentB] - a12[currentB]*a31[currentB]))*zAB + (3.0*a22[currentA]*a33[currentA] - 3.0*a23[currentA]*a32[currentA] - 2.0*a22[currentA]*a33[currentB] + 2.0*a23[currentA]*a32[currentB] + 2.0*a32[currentA]*a23[currentB] - 2.0*a33[currentA]*a22[currentB] + a22[currentB]*a33[currentB] - a23[currentB]*a32[currentB])*xAB^2 + (3.0*a13[currentA]*a32[currentA] - 3.0*a12[currentA]*a33[currentA] - 3.0*a21[currentA]*a33[currentA] + 3.0*a23[currentA]*a31[currentA] + 2.0*a12[currentA]*a33[currentB] - 2.0*a13[currentA]*a32[currentB] - 2.0*a32[currentA]*a13[currentB] + 2.0*a33[currentA]*a12[currentB] + 2.0*a21[currentA]*a33[currentB] - 2.0*a23[currentA]*a31[currentB] - 2.0*a31[currentA]*a23[currentB] + 2.0*a33[currentA]*a21[currentB] - a12[currentB]*a33[currentB] + a13[currentB]*a32[currentB] - a21[currentB]*a33[currentB] + a23[currentB]*a31[currentB])*xAB*yAB + (3.0*a11[currentA]*a33[currentA] - 3.0*a13[currentA]*a31[currentA] - 2.0*a11[currentA]*a33[currentB] + 2.0*a13[currentA]*a31[currentB] + 2.0*a31[currentA]*a13[currentB] - 2.0*a33[currentA]*a11[currentB] + a11[currentB]*a33[currentB] - a13[currentB]*a31[currentB])*yAB^2
					const2t = (3.0*a12[currentA]*a21[currentA] - 3.0*a11[currentA]*a22[currentA] + a11[currentA]*a22[currentB] - a12[currentA]*a21[currentB] - a21[currentA]*a12[currentB] + a22[currentA]*a11[currentB])*zAB^2 + (yAB*(3.0*a11[currentA]*a23[currentA] - 3.0*a13[currentA]*a21[currentA] + 3.0*a11[currentA]*a32[currentA] - 3.0*a12[currentA]*a31[currentA] - a11[currentA]*a23[currentB] + a13[currentA]*a21[currentB] + a21[currentA]*a13[currentB] - a23[currentA]*a11[currentB] - a11[currentA]*a32[currentB] + a12[currentA]*a31[currentB] + a31[currentA]*a12[currentB] - a32[currentA]*a11[currentB]) - xAB*(3.0*a12[currentA]*a23[currentA] - 3.0*a13[currentA]*a22[currentA] + 3.0*a21[currentA]*a32[currentA] - 3.0*a22[currentA]*a31[currentA] - a12[currentA]*a23[currentB] + a13[currentA]*a22[currentB] + a22[currentA]*a13[currentB] - a23[currentA]*a12[currentB] - a21[currentA]*a32[currentB] + a22[currentA]*a31[currentB] + a31[currentA]*a22[currentB] - a32[currentA]*a21[currentB]))*zAB + (3.0*a23[currentA]*a32[currentA] - 3.0*a22[currentA]*a33[currentA] + a22[currentA]*a33[currentB] - a23[currentA]*a32[currentB] - a32[currentA]*a23[currentB] + a33[currentA]*a22[currentB])*xAB^2 + (3.0*a12[currentA]*a33[currentA] - 3.0*a13[currentA]*a32[currentA] + 3.0*a21[currentA]*a33[currentA] - 3.0*a23[currentA]*a31[currentA] - a12[currentA]*a33[currentB] + a13[currentA]*a32[currentB] + a32[currentA]*a13[currentB] - a33[currentA]*a12[currentB] - a21[currentA]*a33[currentB] + a23[currentA]*a31[currentB] + a31[currentA]*a23[currentB] - a33[currentA]*a21[currentB])*xAB*yAB + (3.0*a13[currentA]*a31[currentA] - 3.0*a11[currentA]*a33[currentA] + a11[currentA]*a33[currentB] - a13[currentA]*a31[currentB] - a31[currentA]*a13[currentB] + a33[currentA]*a11[currentB])*yAB^2
					const1t = (a11[currentA]*a22[currentA] - a12[currentA]*a21[currentA])*zAB^2 + (xAB*(a12[currentA]*a23[currentA] - a13[currentA]*a22[currentA] + a21[currentA]*a32[currentA] - a22[currentA]*a31[currentA]) - yAB*(a11[currentA]*a23[currentA] - a13[currentA]*a21[currentA] + a11[currentA]*a32[currentA] - a12[currentA]*a31[currentA]))*zAB + (a22[currentA]*a33[currentA] - a23[currentA]*a32[currentA])*xAB^2 + (a13[currentA]*a32[currentA] - a12[currentA]*a33[currentA] - a21[currentA]*a33[currentA] + a23[currentA]*a31[currentA])*xAB*yAB + (a11[currentA]*a33[currentA] - a13[currentA]*a31[currentA])*yAB^2
					const3n = a11[currentA]*a23[currentA]*a32[currentA] - a11[currentA]*a22[currentA]*a33[currentA] + a12[currentA]*a21[currentA]*a33[currentA] - a12[currentA]*a23[currentA]*a31[currentA] - a13[currentA]*a21[currentA]*a32[currentA] + a13[currentA]*a22[currentA]*a31[currentA] + a11[currentA]*a22[currentA]*a33[currentB] - a11[currentA]*a23[currentA]*a32[currentB] - a11[currentA]*a32[currentA]*a23[currentB] + a11[currentA]*a33[currentA]*a22[currentB] - a12[currentA]*a21[currentA]*a33[currentB] + a12[currentA]*a23[currentA]*a31[currentB] + a12[currentA]*a31[currentA]*a23[currentB] - a12[currentA]*a33[currentA]*a21[currentB] + a13[currentA]*a21[currentA]*a32[currentB] - a13[currentA]*a22[currentA]*a31[currentB] - a13[currentA]*a31[currentA]*a22[currentB] + a13[currentA]*a32[currentA]*a21[currentB] + a21[currentA]*a32[currentA]*a13[currentB] - a21[currentA]*a33[currentA]*a12[currentB] - a22[currentA]*a31[currentA]*a13[currentB] + a22[currentA]*a33[currentA]*a11[currentB] + a23[currentA]*a31[currentA]*a12[currentB] - a23[currentA]*a32[currentA]*a11[currentB] - a11[currentA]*a22[currentB]*a33[currentB] + a11[currentA]*a23[currentB]*a32[currentB] + a12[currentA]*a21[currentB]*a33[currentB] - a12[currentA]*a23[currentB]*a31[currentB] - a13[currentA]*a21[currentB]*a32[currentB] + a13[currentA]*a22[currentB]*a31[currentB] + a21[currentA]*a12[currentB]*a33[currentB] - a21[currentA]*a13[currentB]*a32[currentB] - a22[currentA]*a11[currentB]*a33[currentB] + a22[currentA]*a13[currentB]*a31[currentB] + a23[currentA]*a11[currentB]*a32[currentB] - a23[currentA]*a12[currentB]*a31[currentB] - a31[currentA]*a12[currentB]*a23[currentB] + a31[currentA]*a13[currentB]*a22[currentB] + a32[currentA]*a11[currentB]*a23[currentB] - a32[currentA]*a13[currentB]*a21[currentB] - a33[currentA]*a11[currentB]*a22[currentB] + a33[currentA]*a12[currentB]*a21[currentB] + a11[currentB]*a22[currentB]*a33[currentB] - a11[currentB]*a23[currentB]*a32[currentB] - a12[currentB]*a21[currentB]*a33[currentB] + a12[currentB]*a23[currentB]*a31[currentB] + a13[currentB]*a21[currentB]*a32[currentB] - a13[currentB]*a22[currentB]*a31[currentB]
					const2n = 3.0*a11[currentA]*a22[currentA]*a33[currentA] - 3.0*a11[currentA]*a23[currentA]*a32[currentA] - 3.0*a12[currentA]*a21[currentA]*a33[currentA] + 3.0*a12[currentA]*a23[currentA]*a31[currentA] + 3.0*a13[currentA]*a21[currentA]*a32[currentA] - 3.0*a13[currentA]*a22[currentA]*a31[currentA] - 2.0*a11[currentA]*a22[currentA]*a33[currentB] + 2.0*a11[currentA]*a23[currentA]*a32[currentB] + 2.0*a11[currentA]*a32[currentA]*a23[currentB] - 2.0*a11[currentA]*a33[currentA]*a22[currentB] + 2.0*a12[currentA]*a21[currentA]*a33[currentB] - 2.0*a12[currentA]*a23[currentA]*a31[currentB] - 2.0*a12[currentA]*a31[currentA]*a23[currentB] + 2.0*a12[currentA]*a33[currentA]*a21[currentB] - 2.0*a13[currentA]*a21[currentA]*a32[currentB] + 2.0*a13[currentA]*a22[currentA]*a31[currentB] + 2.0*a13[currentA]*a31[currentA]*a22[currentB] - 2.0*a13[currentA]*a32[currentA]*a21[currentB] - 2.0*a21[currentA]*a32[currentA]*a13[currentB] + 2.0*a21[currentA]*a33[currentA]*a12[currentB] + 2.0*a22[currentA]*a31[currentA]*a13[currentB] - 2.0*a22[currentA]*a33[currentA]*a11[currentB] - 2.0*a23[currentA]*a31[currentA]*a12[currentB] + 2.0*a23[currentA]*a32[currentA]*a11[currentB] + a11[currentA]*a22[currentB]*a33[currentB] - a11[currentA]*a23[currentB]*a32[currentB] - a12[currentA]*a21[currentB]*a33[currentB] + a12[currentA]*a23[currentB]*a31[currentB] + a13[currentA]*a21[currentB]*a32[currentB] - a13[currentA]*a22[currentB]*a31[currentB] - a21[currentA]*a12[currentB]*a33[currentB] + a21[currentA]*a13[currentB]*a32[currentB] + a22[currentA]*a11[currentB]*a33[currentB] - a22[currentA]*a13[currentB]*a31[currentB] - a23[currentA]*a11[currentB]*a32[currentB] + a23[currentA]*a12[currentB]*a31[currentB] + a31[currentA]*a12[currentB]*a23[currentB] - a31[currentA]*a13[currentB]*a22[currentB] - a32[currentA]*a11[currentB]*a23[currentB] + a32[currentA]*a13[currentB]*a21[currentB] + a33[currentA]*a11[currentB]*a22[currentB] - a33[currentA]*a12[currentB]*a21[currentB]
					const1n = 3.0*a11[currentA]*a23[currentA]*a32[currentA] - 3.0*a11[currentA]*a22[currentA]*a33[currentA] + 3.0*a12[currentA]*a21[currentA]*a33[currentA] - 3.0*a12[currentA]*a23[currentA]*a31[currentA] - 3.0*a13[currentA]*a21[currentA]*a32[currentA] + 3.0*a13[currentA]*a22[currentA]*a31[currentA] + a11[currentA]*a22[currentA]*a33[currentB] - a11[currentA]*a23[currentA]*a32[currentB] - a11[currentA]*a32[currentA]*a23[currentB] + a11[currentA]*a33[currentA]*a22[currentB] - a12[currentA]*a21[currentA]*a33[currentB] + a12[currentA]*a23[currentA]*a31[currentB] + a12[currentA]*a31[currentA]*a23[currentB] - a12[currentA]*a33[currentA]*a21[currentB] + a13[currentA]*a21[currentA]*a32[currentB] - a13[currentA]*a22[currentA]*a31[currentB] - a13[currentA]*a31[currentA]*a22[currentB] + a13[currentA]*a32[currentA]*a21[currentB] + a21[currentA]*a32[currentA]*a13[currentB] - a21[currentA]*a33[currentA]*a12[currentB] - a22[currentA]*a31[currentA]*a13[currentB] + a22[currentA]*a33[currentA]*a11[currentB] + a23[currentA]*a31[currentA]*a12[currentB] - a23[currentA]*a32[currentA]*a11[currentB]

					b = 0.5

					fb = (((const4t*b+const3t)*b+const2t)*b+const1t)/((const3n*b+const2n)*b+const1n)
					if fb >= 1.0
						overlapfun = 1.0
					else
						iter = 1
						a = 0.0
						c = 1.0
						while fb < 1 && iter < number_of_iterations_overlap_criterion
							# Pick point in upper half.
							x = 0.5*(b+c)
							fx = (((const4t*x+const3t)*x+const2t)*x+const1t)/((const3n*x+const2n)*x+const1n)
							if fx < fb
								c = x
							else
								b = x
								fb = fx
							end
							# Pick point in lower half.
							x = 0.5*(a+b)
							fx = (((const4t*x+const3t)*x+const2t)*x+const1t)/((const3n*x+const2n)*x+const1n)

							if fx < fb
								a = x
							else
								b = x
								fb = fx
							end

							iter = iter + 1
						end
						overlapfun = fb
					end
					if overlapfun < 1.0
						energy_particle_star += (1.0-overlapfun)^2
					end
				end
			end

			if Z_star[currentA] < lbz || Z_star[currentA] > ubz
				energy_particle_star = energy_particle + 1.0
			end

			if energy_particle_star <= energy_particle
				X[currentA] = X_star[currentA]
				Y[currentA] = Y_star[currentA]
				Z[currentA] = Z_star[currentA]
				acceptance_probability_translation += 1.0
				energy_particle = energy_particle_star
			end

			# Compute new local energy with rotation.
			energy_particle_star = 0.0
			for currentB = [1:currentA-1;currentA+1:number_of_particles]
				xAB = X[currentB] - X[currentA]
				if xAB < -0.5*Lx
					xAB += Lx
				elseif xAB > 0.5*Lx
					xAB -= Lx
				end

				yAB = Y[currentB] - Y[currentA]
				if yAB < -0.5*Ly
					yAB += Ly
				elseif yAB > 0.5*Ly
					yAB -= Ly
				end

				zAB = Z[currentB] - Z[currentA]
				if zAB < -0.5*Lz
					zAB += Lz
				elseif zAB > 0.5*Lz
					zAB -= Lz
				end

				if xAB^2 + yAB^2 + zAB^2 < (RMAX[currentA] + RMAX[currentB])^2
					const4t = (a12_star[currentA]*a21_star[currentA] - a11_star[currentA]*a22_star[currentA] + a11_star[currentA]*a22[currentB] - a12_star[currentA]*a21[currentB] - a21_star[currentA]*a12[currentB] + a22_star[currentA]*a11[currentB] - a11[currentB]*a22[currentB] + a12[currentB]*a21[currentB])*zAB^2 + (yAB*(a11_star[currentA]*a23_star[currentA] - a13_star[currentA]*a21_star[currentA] + a11_star[currentA]*a32_star[currentA] - a12_star[currentA]*a31_star[currentA] - a11_star[currentA]*a23[currentB] + a13_star[currentA]*a21[currentB] + a21_star[currentA]*a13[currentB] - a23_star[currentA]*a11[currentB] - a11_star[currentA]*a32[currentB] + a12_star[currentA]*a31[currentB] + a31_star[currentA]*a12[currentB] - a32_star[currentA]*a11[currentB] + a11[currentB]*a23[currentB] - a13[currentB]*a21[currentB] + a11[currentB]*a32[currentB] - a12[currentB]*a31[currentB]) - xAB*(a12_star[currentA]*a23_star[currentA] - a13_star[currentA]*a22_star[currentA] + a21_star[currentA]*a32_star[currentA] - a22_star[currentA]*a31_star[currentA] - a12_star[currentA]*a23[currentB] + a13_star[currentA]*a22[currentB] + a22_star[currentA]*a13[currentB] - a23_star[currentA]*a12[currentB] - a21_star[currentA]*a32[currentB] + a22_star[currentA]*a31[currentB] + a31_star[currentA]*a22[currentB] - a32_star[currentA]*a21[currentB] + a12[currentB]*a23[currentB] - a13[currentB]*a22[currentB] + a21[currentB]*a32[currentB] - a22[currentB]*a31[currentB]))*zAB + (a23_star[currentA]*a32_star[currentA] - a22_star[currentA]*a33_star[currentA] + a22_star[currentA]*a33[currentB] - a23_star[currentA]*a32[currentB] - a32_star[currentA]*a23[currentB] + a33_star[currentA]*a22[currentB] - a22[currentB]*a33[currentB] + a23[currentB]*a32[currentB])*xAB^2 + (a12_star[currentA]*a33_star[currentA] - a13_star[currentA]*a32_star[currentA] + a21_star[currentA]*a33_star[currentA] - a23_star[currentA]*a31_star[currentA] - a12_star[currentA]*a33[currentB] + a13_star[currentA]*a32[currentB] + a32_star[currentA]*a13[currentB] - a33_star[currentA]*a12[currentB] - a21_star[currentA]*a33[currentB] + a23_star[currentA]*a31[currentB] + a31_star[currentA]*a23[currentB] - a33_star[currentA]*a21[currentB] + a12[currentB]*a33[currentB] - a13[currentB]*a32[currentB] + a21[currentB]*a33[currentB] - a23[currentB]*a31[currentB])*xAB*yAB + (a13_star[currentA]*a31_star[currentA] - a11_star[currentA]*a33_star[currentA] + a11_star[currentA]*a33[currentB] - a13_star[currentA]*a31[currentB] - a31_star[currentA]*a13[currentB] + a33_star[currentA]*a11[currentB] - a11[currentB]*a33[currentB] + a13[currentB]*a31[currentB])*yAB^2
					const3t = (3.0*a11_star[currentA]*a22_star[currentA] - 3.0*a12_star[currentA]*a21_star[currentA] - 2.0*a11_star[currentA]*a22[currentB] + 2.0*a12_star[currentA]*a21[currentB] + 2.0*a21_star[currentA]*a12[currentB] - 2.0*a22_star[currentA]*a11[currentB] + a11[currentB]*a22[currentB] - a12[currentB]*a21[currentB])*zAB^2 + (xAB*(3.0*a12_star[currentA]*a23_star[currentA] - 3.0*a13_star[currentA]*a22_star[currentA] + 3.0*a21_star[currentA]*a32_star[currentA] - 3.0*a22_star[currentA]*a31_star[currentA] - 2.0*a12_star[currentA]*a23[currentB] + 2.0*a13_star[currentA]*a22[currentB] + 2.0*a22_star[currentA]*a13[currentB] - 2.0*a23_star[currentA]*a12[currentB] - 2.0*a21_star[currentA]*a32[currentB] + 2.0*a22_star[currentA]*a31[currentB] + 2.0*a31_star[currentA]*a22[currentB] - 2.0*a32_star[currentA]*a21[currentB] + a12[currentB]*a23[currentB] - a13[currentB]*a22[currentB] + a21[currentB]*a32[currentB] - a22[currentB]*a31[currentB]) - yAB*(3.0*a11_star[currentA]*a23_star[currentA] - 3.0*a13_star[currentA]*a21_star[currentA] + 3.0*a11_star[currentA]*a32_star[currentA] - 3.0*a12_star[currentA]*a31_star[currentA] - 2.0*a11_star[currentA]*a23[currentB] + 2.0*a13_star[currentA]*a21[currentB] + 2.0*a21_star[currentA]*a13[currentB] - 2.0*a23_star[currentA]*a11[currentB] - 2.0*a11_star[currentA]*a32[currentB] + 2.0*a12_star[currentA]*a31[currentB] + 2.0*a31_star[currentA]*a12[currentB] - 2.0*a32_star[currentA]*a11[currentB] + a11[currentB]*a23[currentB] - a13[currentB]*a21[currentB] + a11[currentB]*a32[currentB] - a12[currentB]*a31[currentB]))*zAB + (3.0*a22_star[currentA]*a33_star[currentA] - 3.0*a23_star[currentA]*a32_star[currentA] - 2.0*a22_star[currentA]*a33[currentB] + 2.0*a23_star[currentA]*a32[currentB] + 2.0*a32_star[currentA]*a23[currentB] - 2.0*a33_star[currentA]*a22[currentB] + a22[currentB]*a33[currentB] - a23[currentB]*a32[currentB])*xAB^2 + (3.0*a13_star[currentA]*a32_star[currentA] - 3.0*a12_star[currentA]*a33_star[currentA] - 3.0*a21_star[currentA]*a33_star[currentA] + 3.0*a23_star[currentA]*a31_star[currentA] + 2.0*a12_star[currentA]*a33[currentB] - 2.0*a13_star[currentA]*a32[currentB] - 2.0*a32_star[currentA]*a13[currentB] + 2.0*a33_star[currentA]*a12[currentB] + 2.0*a21_star[currentA]*a33[currentB] - 2.0*a23_star[currentA]*a31[currentB] - 2.0*a31_star[currentA]*a23[currentB] + 2.0*a33_star[currentA]*a21[currentB] - a12[currentB]*a33[currentB] + a13[currentB]*a32[currentB] - a21[currentB]*a33[currentB] + a23[currentB]*a31[currentB])*xAB*yAB + (3.0*a11_star[currentA]*a33_star[currentA] - 3.0*a13_star[currentA]*a31_star[currentA] - 2.0*a11_star[currentA]*a33[currentB] + 2.0*a13_star[currentA]*a31[currentB] + 2.0*a31_star[currentA]*a13[currentB] - 2.0*a33_star[currentA]*a11[currentB] + a11[currentB]*a33[currentB] - a13[currentB]*a31[currentB])*yAB^2
					const2t = (3.0*a12_star[currentA]*a21_star[currentA] - 3.0*a11_star[currentA]*a22_star[currentA] + a11_star[currentA]*a22[currentB] - a12_star[currentA]*a21[currentB] - a21_star[currentA]*a12[currentB] + a22_star[currentA]*a11[currentB])*zAB^2 + (yAB*(3.0*a11_star[currentA]*a23_star[currentA] - 3.0*a13_star[currentA]*a21_star[currentA] + 3.0*a11_star[currentA]*a32_star[currentA] - 3.0*a12_star[currentA]*a31_star[currentA] - a11_star[currentA]*a23[currentB] + a13_star[currentA]*a21[currentB] + a21_star[currentA]*a13[currentB] - a23_star[currentA]*a11[currentB] - a11_star[currentA]*a32[currentB] + a12_star[currentA]*a31[currentB] + a31_star[currentA]*a12[currentB] - a32_star[currentA]*a11[currentB]) - xAB*(3.0*a12_star[currentA]*a23_star[currentA] - 3.0*a13_star[currentA]*a22_star[currentA] + 3.0*a21_star[currentA]*a32_star[currentA] - 3.0*a22_star[currentA]*a31_star[currentA] - a12_star[currentA]*a23[currentB] + a13_star[currentA]*a22[currentB] + a22_star[currentA]*a13[currentB] - a23_star[currentA]*a12[currentB] - a21_star[currentA]*a32[currentB] + a22_star[currentA]*a31[currentB] + a31_star[currentA]*a22[currentB] - a32_star[currentA]*a21[currentB]))*zAB + (3.0*a23_star[currentA]*a32_star[currentA] - 3.0*a22_star[currentA]*a33_star[currentA] + a22_star[currentA]*a33[currentB] - a23_star[currentA]*a32[currentB] - a32_star[currentA]*a23[currentB] + a33_star[currentA]*a22[currentB])*xAB^2 + (3.0*a12_star[currentA]*a33_star[currentA] - 3.0*a13_star[currentA]*a32_star[currentA] + 3.0*a21_star[currentA]*a33_star[currentA] - 3.0*a23_star[currentA]*a31_star[currentA] - a12_star[currentA]*a33[currentB] + a13_star[currentA]*a32[currentB] + a32_star[currentA]*a13[currentB] - a33_star[currentA]*a12[currentB] - a21_star[currentA]*a33[currentB] + a23_star[currentA]*a31[currentB] + a31_star[currentA]*a23[currentB] - a33_star[currentA]*a21[currentB])*xAB*yAB + (3.0*a13_star[currentA]*a31_star[currentA] - 3.0*a11_star[currentA]*a33_star[currentA] + a11_star[currentA]*a33[currentB] - a13_star[currentA]*a31[currentB] - a31_star[currentA]*a13[currentB] + a33_star[currentA]*a11[currentB])*yAB^2
					const1t = (a11_star[currentA]*a22_star[currentA] - a12_star[currentA]*a21_star[currentA])*zAB^2 + (xAB*(a12_star[currentA]*a23_star[currentA] - a13_star[currentA]*a22_star[currentA] + a21_star[currentA]*a32_star[currentA] - a22_star[currentA]*a31_star[currentA]) - yAB*(a11_star[currentA]*a23_star[currentA] - a13_star[currentA]*a21_star[currentA] + a11_star[currentA]*a32_star[currentA] - a12_star[currentA]*a31_star[currentA]))*zAB + (a22_star[currentA]*a33_star[currentA] - a23_star[currentA]*a32_star[currentA])*xAB^2 + (a13_star[currentA]*a32_star[currentA] - a12_star[currentA]*a33_star[currentA] - a21_star[currentA]*a33_star[currentA] + a23_star[currentA]*a31_star[currentA])*xAB*yAB + (a11_star[currentA]*a33_star[currentA] - a13_star[currentA]*a31_star[currentA])*yAB^2
					const3n = a11_star[currentA]*a23_star[currentA]*a32_star[currentA] - a11_star[currentA]*a22_star[currentA]*a33_star[currentA] + a12_star[currentA]*a21_star[currentA]*a33_star[currentA] - a12_star[currentA]*a23_star[currentA]*a31_star[currentA] - a13_star[currentA]*a21_star[currentA]*a32_star[currentA] + a13_star[currentA]*a22_star[currentA]*a31_star[currentA] + a11_star[currentA]*a22_star[currentA]*a33[currentB] - a11_star[currentA]*a23_star[currentA]*a32[currentB] - a11_star[currentA]*a32_star[currentA]*a23[currentB] + a11_star[currentA]*a33_star[currentA]*a22[currentB] - a12_star[currentA]*a21_star[currentA]*a33[currentB] + a12_star[currentA]*a23_star[currentA]*a31[currentB] + a12_star[currentA]*a31_star[currentA]*a23[currentB] - a12_star[currentA]*a33_star[currentA]*a21[currentB] + a13_star[currentA]*a21_star[currentA]*a32[currentB] - a13_star[currentA]*a22_star[currentA]*a31[currentB] - a13_star[currentA]*a31_star[currentA]*a22[currentB] + a13_star[currentA]*a32_star[currentA]*a21[currentB] + a21_star[currentA]*a32_star[currentA]*a13[currentB] - a21_star[currentA]*a33_star[currentA]*a12[currentB] - a22_star[currentA]*a31_star[currentA]*a13[currentB] + a22_star[currentA]*a33_star[currentA]*a11[currentB] + a23_star[currentA]*a31_star[currentA]*a12[currentB] - a23_star[currentA]*a32_star[currentA]*a11[currentB] - a11_star[currentA]*a22[currentB]*a33[currentB] + a11_star[currentA]*a23[currentB]*a32[currentB] + a12_star[currentA]*a21[currentB]*a33[currentB] - a12_star[currentA]*a23[currentB]*a31[currentB] - a13_star[currentA]*a21[currentB]*a32[currentB] + a13_star[currentA]*a22[currentB]*a31[currentB] + a21_star[currentA]*a12[currentB]*a33[currentB] - a21_star[currentA]*a13[currentB]*a32[currentB] - a22_star[currentA]*a11[currentB]*a33[currentB] + a22_star[currentA]*a13[currentB]*a31[currentB] + a23_star[currentA]*a11[currentB]*a32[currentB] - a23_star[currentA]*a12[currentB]*a31[currentB] - a31_star[currentA]*a12[currentB]*a23[currentB] + a31_star[currentA]*a13[currentB]*a22[currentB] + a32_star[currentA]*a11[currentB]*a23[currentB] - a32_star[currentA]*a13[currentB]*a21[currentB] - a33_star[currentA]*a11[currentB]*a22[currentB] + a33_star[currentA]*a12[currentB]*a21[currentB] + a11[currentB]*a22[currentB]*a33[currentB] - a11[currentB]*a23[currentB]*a32[currentB] - a12[currentB]*a21[currentB]*a33[currentB] + a12[currentB]*a23[currentB]*a31[currentB] + a13[currentB]*a21[currentB]*a32[currentB] - a13[currentB]*a22[currentB]*a31[currentB]
					const2n = 3.0*a11_star[currentA]*a22_star[currentA]*a33_star[currentA] - 3.0*a11_star[currentA]*a23_star[currentA]*a32_star[currentA] - 3.0*a12_star[currentA]*a21_star[currentA]*a33_star[currentA] + 3.0*a12_star[currentA]*a23_star[currentA]*a31_star[currentA] + 3.0*a13_star[currentA]*a21_star[currentA]*a32_star[currentA] - 3.0*a13_star[currentA]*a22_star[currentA]*a31_star[currentA] - 2.0*a11_star[currentA]*a22_star[currentA]*a33[currentB] + 2.0*a11_star[currentA]*a23_star[currentA]*a32[currentB] + 2.0*a11_star[currentA]*a32_star[currentA]*a23[currentB] - 2.0*a11_star[currentA]*a33_star[currentA]*a22[currentB] + 2.0*a12_star[currentA]*a21_star[currentA]*a33[currentB] - 2.0*a12_star[currentA]*a23_star[currentA]*a31[currentB] - 2.0*a12_star[currentA]*a31_star[currentA]*a23[currentB] + 2.0*a12_star[currentA]*a33_star[currentA]*a21[currentB] - 2.0*a13_star[currentA]*a21_star[currentA]*a32[currentB] + 2.0*a13_star[currentA]*a22_star[currentA]*a31[currentB] + 2.0*a13_star[currentA]*a31_star[currentA]*a22[currentB] - 2.0*a13_star[currentA]*a32_star[currentA]*a21[currentB] - 2.0*a21_star[currentA]*a32_star[currentA]*a13[currentB] + 2.0*a21_star[currentA]*a33_star[currentA]*a12[currentB] + 2.0*a22_star[currentA]*a31_star[currentA]*a13[currentB] - 2.0*a22_star[currentA]*a33_star[currentA]*a11[currentB] - 2.0*a23_star[currentA]*a31_star[currentA]*a12[currentB] + 2.0*a23_star[currentA]*a32_star[currentA]*a11[currentB] + a11_star[currentA]*a22[currentB]*a33[currentB] - a11_star[currentA]*a23[currentB]*a32[currentB] - a12_star[currentA]*a21[currentB]*a33[currentB] + a12_star[currentA]*a23[currentB]*a31[currentB] + a13_star[currentA]*a21[currentB]*a32[currentB] - a13_star[currentA]*a22[currentB]*a31[currentB] - a21_star[currentA]*a12[currentB]*a33[currentB] + a21_star[currentA]*a13[currentB]*a32[currentB] + a22_star[currentA]*a11[currentB]*a33[currentB] - a22_star[currentA]*a13[currentB]*a31[currentB] - a23_star[currentA]*a11[currentB]*a32[currentB] + a23_star[currentA]*a12[currentB]*a31[currentB] + a31_star[currentA]*a12[currentB]*a23[currentB] - a31_star[currentA]*a13[currentB]*a22[currentB] - a32_star[currentA]*a11[currentB]*a23[currentB] + a32_star[currentA]*a13[currentB]*a21[currentB] + a33_star[currentA]*a11[currentB]*a22[currentB] - a33_star[currentA]*a12[currentB]*a21[currentB]
					const1n = 3.0*a11_star[currentA]*a23_star[currentA]*a32_star[currentA] - 3.0*a11_star[currentA]*a22_star[currentA]*a33_star[currentA] + 3.0*a12_star[currentA]*a21_star[currentA]*a33_star[currentA] - 3.0*a12_star[currentA]*a23_star[currentA]*a31_star[currentA] - 3.0*a13_star[currentA]*a21_star[currentA]*a32_star[currentA] + 3.0*a13_star[currentA]*a22_star[currentA]*a31_star[currentA] + a11_star[currentA]*a22_star[currentA]*a33[currentB] - a11_star[currentA]*a23_star[currentA]*a32[currentB] - a11_star[currentA]*a32_star[currentA]*a23[currentB] + a11_star[currentA]*a33_star[currentA]*a22[currentB] - a12_star[currentA]*a21_star[currentA]*a33[currentB] + a12_star[currentA]*a23_star[currentA]*a31[currentB] + a12_star[currentA]*a31_star[currentA]*a23[currentB] - a12_star[currentA]*a33_star[currentA]*a21[currentB] + a13_star[currentA]*a21_star[currentA]*a32[currentB] - a13_star[currentA]*a22_star[currentA]*a31[currentB] - a13_star[currentA]*a31_star[currentA]*a22[currentB] + a13_star[currentA]*a32_star[currentA]*a21[currentB] + a21_star[currentA]*a32_star[currentA]*a13[currentB] - a21_star[currentA]*a33_star[currentA]*a12[currentB] - a22_star[currentA]*a31_star[currentA]*a13[currentB] + a22_star[currentA]*a33_star[currentA]*a11[currentB] + a23_star[currentA]*a31_star[currentA]*a12[currentB] - a23_star[currentA]*a32_star[currentA]*a11[currentB]

					b = 0.5

					fb = (((const4t*b+const3t)*b+const2t)*b+const1t)/((const3n*b+const2n)*b+const1n)
					if fb >= 1.0
						overlapfun = 1.0
					else
						iter = 1
						a = 0.0
						c = 1.0
						while fb < 1 && iter < number_of_iterations_overlap_criterion
							# Pick point in upper half.
							x = 0.5*(b+c)
							fx = (((const4t*x+const3t)*x+const2t)*x+const1t)/((const3n*x+const2n)*x+const1n)
							if fx < fb
								c = x
							else
								b = x
								fb = fx
							end
							# Pick point in lower half.
							x = 0.5*(a+b)
							fx = (((const4t*x+const3t)*x+const2t)*x+const1t)/((const3n*x+const2n)*x+const1n)

							if fx < fb
								a = x
							else
								b = x
								fb = fx
							end

							iter = iter + 1
						end
						overlapfun = fb
					end
					if overlapfun < 1.0
						energy_particle_star += (1.0-overlapfun)^2
					end
				end
			end

			angle_to_z_axis = acos((cos(THETA1_star[currentA])*cos(THETA2_star[currentA]))/(abs(sin(THETA2_star[currentA]))^2 + abs(cos(THETA1_star[currentA])*cos(THETA2_star[currentA]))^2 + abs(cos(THETA2_star[currentA])*sin(THETA1_star[currentA]))^2)^(1/2))
			#if angle_to_z_axis > ubangle
			#	energy_particle_star = energy_particle + 1.0 # Just a  value that prevents the rotation.
			#end

			if energy_particle_star <= energy_particle && angle_to_z_axis <= ubangle
				THETA1[currentA] = THETA1_star[currentA]
				THETA2[currentA] = THETA2_star[currentA]
				THETA3[currentA] = THETA3_star[currentA]

				cA1[currentA] = cA1_star[currentA]
				cA2[currentA] = cA2_star[currentA]
				sA1[currentA] = sA1_star[currentA]
				sA2[currentA] = sA2_star[currentA]
				cA1d[currentA] = cA1d_star[currentA]
				cA2d[currentA] = cA2d_star[currentA]
				cA3d[currentA] = cA3d_star[currentA]
				sA1d[currentA] = sA1d_star[currentA]
				sA2d[currentA] = sA2d_star[currentA]
				sA3d[currentA] = sA3d_star[currentA]
				a11[currentA] = a11_star[currentA]
				a12[currentA] = a12_star[currentA]
				a13[currentA] = a13_star[currentA]
				a21[currentA] = a21_star[currentA]
				a22[currentA] = a22_star[currentA]
				a23[currentA] = a23_star[currentA]
				a31[currentA] = a31_star[currentA]
				a32[currentA] = a32_star[currentA]
				a33[currentA] = a33_star[currentA]

				acceptance_probability_rotation += 1.0
				energy_particle = energy_particle_star
			end

		end
		acceptance_probability_translation /= number_of_particles
		acceptance_probability_rotation /= number_of_particles

		if acceptance_probability_translation <= acceptance_probability_target
			sigma_translation *= 0.95
		else
			sigma_translation *= 1.05
			sigma_translation = min(sigma_translation, sigma_translation_max)
		end

		if acceptance_probability_rotation <= acceptance_probability_target
			sigma_rotation *= 0.95
		else
			sigma_rotation *= 1.05
			sigma_rotation = min(sigma_rotation, sigma_rotation_max)
		end

		m1x = 0.0
		m2x = 0.0
		m1y = 0.0
		m2y = 0.0
		m1z = 0.0
		m2z = 0.0
		for currentA = 1:number_of_particles
			angle_to_x_axis = acos((cos(THETA2[currentA])*cos(THETA3[currentA]))/(abs(sin(THETA3[currentA]))^2 + abs(cos(THETA2[currentA])*cos(THETA3[currentA]))^2 + abs(cos(THETA3[currentA])*sin(THETA2[currentA]))^2)^(1/2))
			angle_to_y_axis = acos((cos(THETA3[currentA])*cos(THETA1[currentA]))/(abs(sin(THETA1[currentA]))^2 + abs(cos(THETA3[currentA])*cos(THETA1[currentA]))^2 + abs(cos(THETA1[currentA])*sin(THETA3[currentA]))^2)^(1/2))
			angle_to_z_axis = acos((cos(THETA1[currentA])*cos(THETA2[currentA]))/(abs(sin(THETA2[currentA]))^2 + abs(cos(THETA1[currentA])*cos(THETA2[currentA]))^2 + abs(cos(THETA2[currentA])*sin(THETA1[currentA]))^2)^(1/2))

			m1x += angle_to_x_axis
			m2x += angle_to_x_axis^2
			m1y += angle_to_y_axis
			m2y += angle_to_y_axis^2
			m1z += angle_to_z_axis
			m2z += angle_to_z_axis^2

		end
		m1x /= number_of_particles
		m2x /= number_of_particles
		m1y /= number_of_particles
		m2y /= number_of_particles
		m1z /= number_of_particles
		m2z /= number_of_particles
		sx = sqrt(m2x - m1x^2)
		sy = sqrt(m2y - m1y^2)
		sz = sqrt(m2z - m1z^2)
		println((sx, sy, sz))
		#max_val = 0.0
		#for currentA = 1:number_of_particles
		#	angle_to_z_axis = acos((cos(THETA1[currentA])*cos(THETA2[currentA]))/(abs(sin(THETA2[currentA]))^2 + #abs(cos(THETA1[currentA])*cos(THETA2[currentA]))^2 + abs(cos(THETA2[currentA])*sin(THETA1[currentA]))^2)^(1/2))
		#	max_val = max(max_val, angle_to_z_axis)
		#end
		#println(max_val)
  end

	# Translate angles to 0 <= theta < 2*pi.
	THETA1 = mod(THETA1, 2 * pi)
	THETA2 = mod(THETA2, 2 * pi)
	THETA3 = mod(THETA3, 2 * pi)

	return (X, Y, Z, THETA1, THETA2, THETA3)
end
