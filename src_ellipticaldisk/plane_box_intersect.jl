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
	
	
	
	
	w3 = q31[current_particle] * vx + q32[current_particle] * vy + q33[current_particle] * vz	
						
					if sign(maximum(w3)) - sign(minimum(w3)) >= 1.0 # Different sides of plane or one point intersecting plane
						if minimum(vx.^2 + vy.^2 + vz.^2) <= RMAX[current_particle]^2
							push!(cell_lists[current_cell_x, current_cell_y, current_cell_z], current_particle)
						end
					end
	

end