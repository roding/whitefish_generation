###############################################################################
#
# PLANE_BOX_INTERSECT
#
# Returns whether a plane, specified by a point (px, py, pz) and a normal
# vector (nx, ny, nz), intersects with a box, specified by its lower and upper 
# bounds, (lbx, ubx, lby, uby, lbz, ubz).
#
###############################################################################

function plane_box_intersect(	px::Float64,
							py::Float64,
							pz::Float64,
							nx::Float64,
							ny::Float64,
							nz::Float64,
							lbx::Float64,
							ubx::Float64, 
							lby::Float64,
							uby::Float64, 
							lbz::Float64,
							ubz::Float64)
	
	# Compute absolute coordinates of the 8 box vertices.
	x_vertex::Array{Float64, 1} = [lbx, lbx, lbx, lbx, ubx, ubx, ubx, ubx]
	y_vertex::Array{Float64, 1} = [lby, lby, uby, uby, lby, lby, uby, uby]
	z_vertex::Array{Float64, 1} = [lbz, ubz, lbz, ubz, lbz, ubz, lbz, ubz]
	
	# Compute vertex coordinates relative to point in plane.
	x_vertex = x_vertex - px
	y_vertex = y_vertex - py
	z_vertex = z_vertex - pz
	
	# Compute coordinates along axis defined by normal vector and relative to point.
	w::Array{Float64, 1} = nx * x_vertex + ny * y_vertex + nz * z_vertex	
		
	# Check if vertices are on same or opposite sides of plane. If the maximum sign 
	# difference is 0.0, they are on the same side (or both are in the plane, but 
	# this case is of no relevance), if it is 1.0, one point is outside the plane and
	# one inside, if it is 2.0, points are  strictly on different sides.
	if sign(maximum(w)) - sign(minimum(w)) >= 1.0
		return true
	else
		return false
	end

end