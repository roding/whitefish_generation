function bounding_box(		X::Array{Float64,1},
							Y::Array{Float64,1},
							Z::Array{Float64,1},
							THETA1::Array{Float64,1},
							THETA2::Array{Float64,1},
							THETA3::Array{Float64,1},
							R1::Array{Float64,1}, 
							R2::Array{Float64,1})
		
	# Ellipse parameters.
	const number_of_particles::Int64 = length(X)
	
	const RMAX = Array(Float64, number_of_particles)
	for current_particle = 1:number_of_particles
		RMAX[current_particle] = maximum((R1[current_particle], R2[current_particle]))
	end
	
	const c1::Array{Float64,1} = cos(THETA1)
	const c2::Array{Float64,1} = cos(THETA2)
	const c3::Array{Float64,1} = cos(THETA3)
	const s1::Array{Float64,1} = sin(THETA1)
	const s2::Array{Float64,1} = sin(THETA2)
	const s3::Array{Float64,1} = sin(THETA3)
	
	const q11::Array{Float64,1} = c2.*c3
	const q12::Array{Float64,1} = -c2.*s3
	const q13::Array{Float64,1} = s2
	const q21::Array{Float64,1} = c1.*s3 + c3.*s1.*s2
	const q22::Array{Float64,1} = c1.*c3 - s1.*s2.*s3
	const q23::Array{Float64,1} = -c2.*s1
	const q31::Array{Float64,1} = s1.*s3 - c1.*c3.*s2
	const q32::Array{Float64,1} = c3.*s1 + c1.*s2.*s3
	const q33::Array{Float64,1} = c1.*c2

end