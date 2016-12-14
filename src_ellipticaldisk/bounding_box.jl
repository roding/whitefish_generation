###############################################################################
#
# BOUNDING_BOX
#
# Returns lower and upper bounds for the bounding box of an elliptical disk.
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
	lbx::Float64 = 0.0
	ubx::Float64 = 0.0
	lby::Float64 = 0.0
	uby::Float64 = 0.0
	lbz::Float64 = 0.0
	ubz::Float64 = 0.0
	
	# Bounds in x direction.
	t_extrema[1] = atan(-R2/R1 * tan(THETA3))
	t_extrema[2] = t_extrema[1] + pi
	
	val_extrema = X + R1 * cos(t_extrema)*cos(THETA2)*cos(THETA3) - r2*cos(THETA2)*sin(t_extrema)*sin(THETA3);
x_min = min(x_extrema);
x_max = max(x_extrema);

disp([min(x) x_min max(x) x_max])

t_extrema(1) = atan(r2/r1 * (cos(THETA1)*cos(THETA3) - sin(THETA1)*sin(THETA2)*sin(THETA3)) / (cos(THETA1)*sin(THETA3) + cos(THETA3)*sin(THETA1)*sin(THETA2)));
t_extrema(2) = t_extrema(1) + pi;
y_extrema = yc + r1*cos(t_extrema)*(cos(THETA1)*sin(THETA3) + cos(THETA3)*sin(THETA1)*sin(THETA2)) + r2*sin(t_extrema)*(cos(THETA1)*cos(THETA3) - sin(THETA1)*sin(THETA2)*sin(THETA3));
y_min = min(y_extrema);
y_max = max(y_extrema);

disp([min(y) y_min max(y) y_max])

if isnan((cos(THETA3)*sin(THETA1) + cos(THETA1)*sin(THETA2)*sin(THETA3)) / (sin(THETA1)*sin(THETA3) - cos(THETA1)*cos(THETA3)*sin(THETA2)))
    t_extrema(1) = atan(r2/r1 * 1/tan(THETA3));
else
    t_extrema(1) = atan(r2/r1 * (cos(THETA3)*sin(THETA1) + cos(THETA1)*sin(THETA2)*sin(THETA3)) / (sin(THETA1)*sin(THETA3) - cos(THETA1)*cos(THETA3)*sin(THETA2)));
end
t_extrema(2) = t_extrema(1) + pi;
z_extrema = zc + r1*cos(t_extrema)*(sin(THETA1)*sin(THETA3) - cos(THETA1)*cos(THETA3)*sin(THETA2)) + r2*sin(t_extrema)*(cos(THETA3)*sin(THETA1) + cos(THETA1)*sin(THETA2)*sin(THETA3));
z_min = min(z_extrema);
z_max = max(z_extrema);
	

end