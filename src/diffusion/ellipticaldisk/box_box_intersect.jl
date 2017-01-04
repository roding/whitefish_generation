###############################################################################
#
# BOX_BOX_INTERSECT
#
# Returns whether two boxes, specified by their lower and upper bounds
# (lbx1, ubx1, lby1, uby1, lbz1, ubz1) and (lbx2, ubx2, lby2, uby2, lbz2, ubz2)
# are intersecting. The calculation is performed "modulo the simulation domain"
# i.e. using the closest box image. 
#
###############################################################################

function box_box_intersect(	lbx1::Float64,
							ubx1::Float64, 
							lby1::Float64,
							uby1::Float64, 
							lbz1::Float64,
							ubz1::Float64,
							lbx2::Float64,
							ubx2::Float64, 
							lby2::Float64,
							uby2::Float64, 
							lbz2::Float64,
							ubz2::Float64,
							Lx::Float64,
							Ly::Float64,
							Lz::Float64)
		
	# Move box 1 to its image closest to box 2.
	if 0.5 * (lbx1 + ubx1) - 0.5 * (lbx2 + ubx2) < -0.5*Lx
		lbx1 = lbx1 + Lx
		ubx1 = ubx1 + Lx
	elseif 0.5 * (lbx1 + ubx1) - 0.5 * (lbx2 + ubx2) > 0.5*Lx
		lbx1 = lbx1 - Lx
		ubx1 = ubx1 - Lx
	end
	
	if 0.5 * (lby1 + uby1) - 0.5 * (lby2 + uby2) < -0.5*Ly
		lby1 = lby1 + Ly
		uby1 = uby1 + Ly
	elseif 0.5 * (lby1 + uby1) - 0.5 * (lby2 + uby2) > 0.5*Ly
		lby1 = lby1 - Ly
		uby1 = uby1 - Ly
	end
	
	if 0.5 * (lbz1 + ubz1) - 0.5 * (lbz2 + ubz2) < -0.5*Lz
		lbz1 = lbz1 + Lz
		ubz1 = ubz1 + Lz
	elseif 0.5 * (lbz1 + ubz1) - 0.5 * (lbz2 + ubz2) > 0.5*Lz
		lbz1 = lbz1 - Lz
		ubz1 = ubz1 - Lz
	end
	
	#println((lbx1,ubx1,lby1,uby1,lbz1,ubz1))

	return 	(abs( (lbx1 + ubx1) - (lbx2 + ubx2)) <= ubx1 - lbx1 + ubx2 - lbx2 ) &&
			(abs( (lby1 + uby1) - (lby2 + uby2)) <= uby1 - lby1 + uby2 - lby2 ) &&
			(abs( (lbz1 + ubz1) - (lbz2 + ubz2)) <= ubz1 - lbz1 + ubz2 - lbz2 )
end