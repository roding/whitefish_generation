###############################################################################
#
# BOUNDING_BOX
#
# Computes lower and upper bounds for the bounding box of an elliptical disk.
#
###############################################################################

function bounding_box(X::Float64,
					Y::Float64,
					Z::Float64,
					THETA1::Float64,
					THETA2::Float64,
					THETA3::Float64,
					R1::Float64, 
					R2::Float64)
		
	t_extrema = Array(Float64, 2)
	val_extrema = Array(Float64, 2)
	
	# Bounds in x direction.
	t_extrema[1] = atan(-R2/R1 * tan(THETA3))
	t_extrema[2] = t_extrema[1] + pi
	
	val_extrema = X + R1*cos(t_extrema)*cos(THETA2)*cos(THETA3) - R2*sin(t_extrema)*cos(THETA2)*sin(THETA3)

	lbx::Float64 = minimum(val_extrema)
	ubx::Float64 = maximum(val_extrema)

	# Bounds in y direction.
	t_extrema[1] = atan(R2/R1 * (cos(THETA1)*cos(THETA3) - sin(THETA1)*sin(THETA2)*sin(THETA3)) / (cos(THETA1)*sin(THETA3) + cos(THETA3)*sin(THETA1)*sin(THETA2)))
	t_extrema[2] = t_extrema[1] + pi
	
	val_extrema = Y + R1*cos(t_extrema)*(cos(THETA1)*sin(THETA3) + cos(THETA3)*sin(THETA1)*sin(THETA2)) + R2*sin(t_extrema)*(cos(THETA1)*cos(THETA3) - sin(THETA1)*sin(THETA2)*sin(THETA3));

	lby::Float64 = minimum(val_extrema)
	uby::Float64 = maximum(val_extrema)

	# Bounds in z direction.
	if isnan( (cos(THETA3)*sin(THETA1) + cos(THETA1)*sin(THETA2)*sin(THETA3)) / (sin(THETA1)*sin(THETA3) - cos(THETA1)*cos(THETA3)*sin(THETA2)) )
		t_extrema[1] = atan(R2/R1 * 1/tan(THETA3))
	else
		t_extrema[1] = atan(R2/R1 * (cos(THETA3)*sin(THETA1) + cos(THETA1)*sin(THETA2)*sin(THETA3)) / (sin(THETA1)*sin(THETA3) - cos(THETA1)*cos(THETA3)*sin(THETA2)))
	end
	t_extrema[2] = t_extrema[1] + pi

	val_extrema = Z + R1*cos(t_extrema)*(sin(THETA1)*sin(THETA3) - cos(THETA1)*cos(THETA3)*sin(THETA2)) + R2*sin(t_extrema)*(cos(THETA3)*sin(THETA1) + cos(THETA1)*sin(THETA2)*sin(THETA3))
	
	lbz::Float64 = minimum(val_extrema)
	ubz::Float64 = maximum(val_extrema)

	return (lbx, ubx, lby, uby, lbz, ubz)
end