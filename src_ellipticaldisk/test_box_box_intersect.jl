include("box_box_intersect.jl")

function test_box_box_intersect()
	
	lbx1::Float64 = 0.0
	ubx1::Float64 = 1.0
	lby1::Float64 = 0.0
	uby1::Float64 = 1.0
	lbz1::Float64 = 0.0
	ubz1::Float64 = 1.0
	
	lbx2::Float64 = -1.0
	ubx2::Float64 = -0.01
	lby2::Float64 = -1.0
	uby2::Float64 = 0.0
	lbz2::Float64 = -1.0
	ubz2::Float64 = 0.0
	
	is_intersecting = box_box_intersect(	lbx1,
										ubx1,
										lby1,
										uby1,
										lbz1,
										ubz1,
										lbx2,
										ubx2,
										lby2,
										uby2,
										lbz2,
										ubz2)
										
	println(is_intersecting)

end

test_box_box_intersect()