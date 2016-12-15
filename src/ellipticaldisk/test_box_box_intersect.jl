include("box_box_intersect.jl")

function test_box_box_intersect()
	
	lbx1::Float64 = 0.0
	ubx1::Float64 = 1.0
	lby1::Float64 = 0.0
	uby1::Float64 = 1.0
	lbz1::Float64 = 0.0
	ubz1::Float64 = 1.0
	
	lbx2::Float64 = 9.0
	ubx2::Float64 = 10.0
	lby2::Float64 = 9.0
	uby2::Float64 = 10.0
	lbz2::Float64 = 9.0
	ubz2::Float64 = 10.5
	
	Lx::Float64 = 10.0
	Ly::Float64 = 10.0
	Lz::Float64 = 10.0
	
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
										ubz2,
										Lx,
										Ly,
										Lx)
										
	println(is_intersecting)

end

test_box_box_intersect()