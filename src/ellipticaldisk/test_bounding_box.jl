include("bounding_box.jl")

function test_bounding_box()
	
	X::Float64 = randn()
	Y::Float64 = randn()
	Z::Float64 = randn()
	THETA1::Float64 = 0.0#randn()
	THETA2::Float64 = 0.0#randn()
	THETA3::Float64 = 0.0#randn()
	R1::Float64 = 1.0
	R2::Float64 = 2.0
	
	(	lbx::Float64, 
		ubx::Float64, 
		lby::Float64, 
		uby::Float64, 
		lbz::Float64, 
		ubz::Float64) = bounding_box(	X,
									Y,
									Z,
									THETA1,
									THETA2,
									THETA3,
									R1, 
									R2)
	
	t::Array{Float64, 1} = linspace(0, 2*pi, 100000)
	
	x::Array{Float64, 1} = X + R1*cos(t)*cos(THETA2)*cos(THETA3) - R2*cos(THETA2)*sin(t)*sin(THETA3)
	y::Array{Float64, 1} = Y + R1*cos(t)*(cos(THETA1)*sin(THETA3) + cos(THETA3)*sin(THETA1)*sin(THETA2)) + R2*sin(t)*(cos(THETA1)*cos(THETA3) - sin(THETA1)*sin(THETA2)*sin(THETA3))
	z::Array{Float64, 1} = Z + R1*cos(t)*(sin(THETA1)*sin(THETA3) - cos(THETA1)*cos(THETA3)*sin(THETA2)) + R2*sin(t)*(cos(THETA3)*sin(THETA1) + cos(THETA1)*sin(THETA2)*sin(THETA3))
	
	error::Float64 = norm([minimum(x) - lbx, maximum(x) - ubx, minimum(y) - lby, maximum(y) - uby, minimum(z) - lbz, maximum(z) - ubz])
	println((THETA1,THETA2,THETA3))
	println(error)

	nothing
end

test_bounding_box()
