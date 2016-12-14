include("bounding_box.jl")

function test_bounding_box()
	
	X::Float64 = randn()
	Y::Float64 = randn()
	Z::Float64 = randn()
	THETA1::Float64 = randn()
	THETA2::Float64 = randn()
	THETA3::Float64 = randn()
	R1::Array{Float64,1} = 1.0
	R2::Array{Float64,1} = 2.0
	
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
	x = xc + r1*cos(t)*cos(theta2)*cos(theta3) - r2*cos(theta2)*sin(t)*sin(theta3);
y = yc + r1*cos(t)*(cos(theta1)*sin(theta3) + cos(theta3)*sin(theta1)*sin(theta2)) + r2*sin(t)*(cos(theta1)*cos(theta3) - sin(theta1)*sin(theta2)*sin(theta3));
z = zc + r1*cos(t)*(sin(theta1)*sin(theta3) - cos(theta1)*cos(theta3)*sin(theta2)) + r2*sin(t)*(cos(theta3)*sin(theta1) + cos(theta1)*sin(theta2)*sin(theta3));


	nothing
end

test_create_cell_lists()