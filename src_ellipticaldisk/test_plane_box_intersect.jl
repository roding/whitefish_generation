function test_plane_box_intersect()
	
	px::Float64 = 0.0
	py::Float64 = 0.0
	pz::Float64 = 0.0
	nx::Float64 = 1.0
	ny::Float64 = 0.0
	nz::Float64 = 0.0
	lbx::Float64 = 0.0
	ubx::Float64 = 1.0
	lby::Float64 = 0.0
	uby::Float64 = 1.0
	lbz::Float64 = 0.0
	ubz::Float64 = 1.0
	
	is_intersecting = plane_box_intersect(	px
										py
										pz
										nx
										ny
										nz
										lbx
										ubx 
										lby
										uby 
										lbz
										ubz)
										
	println(is_intersecting)

end

is_intersecting()