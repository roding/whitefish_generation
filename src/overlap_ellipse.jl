function overlap_ellipse(	xAB::Float64,
						yAB::Float64,
						zAB::Float64,
						a11::Float64,
						a12::Float64,
						a13::Float64,
						a21::Float64,
						a22::Float64,
						a23::Float64,
						a31::Float64,
						a32::Float64,
						a33::Float64,
						b11::Float64,
						b12::Float64,
						b13::Float64,
						b21::Float64,
						b22::Float64,
						b23::Float64,
						b31::Float64,
						b32::Float64,
						b33::Float64)

	number_of_iterations_overlap_criterion::Int64 = 20
	
	a::Float64 = 0.0
	b::Float64 = 0.0
	c::Float64 = 0.0
	x::Float64 = 0.0
	fb::Float64 = 0.0
	fx::Float64 = 0.0
	overlapfun::Float64 = 0.0
	iter::Int64 = 0
	
	const4t::Float64 = (a12 * a21 - a11 * a22 + a11 * b22 - a12 * b21 - a21 * b12 + a22 * b11 - b11 * b22 + b12 * b21) * zAB^2 + (yAB * (a11 * a23 - a13 * a21 + a11 * a32 - a12 * a31 - a11 * b23 + a13 * b21 + a21 * b13 - a23 * b11 - a11 * b32 + a12 * b31 + a31 * b12 - a32 * b11 + b11 * b23 - b13 * b21 + b11 * b32 - b12 * b31) - xAB * (a12 * a23 - a13 * a22 + a21 * a32 - a22 * a31 - a12 * b23 + a13 * b22 + a22 * b13 - a23 * b12 - a21 * b32 + a22 * b31 + a31 * b22 - a32 * b21 + b12 * b23 - b13 * b22 + b21 * b32 - b22 * b31)) * zAB + (a23 * a32 - a22 * a33 + a22 * b33 - a23 * b32 - a32 * b23 + a33 * b22 - b22 * b33 + b23 * b32) * xAB^2 + (a12 * a33 - a13 * a32 + a21 * a33 - a23 * a31 - a12 * b33 + a13 * b32 + a32 * b13 - a33 * b12 - a21 * b33 + a23 * b31 + a31 * b23 - a33 * b21 + b12 * b33 - b13 * b32 + b21 * b33 - b23 * b31) * xAB * yAB + (a13 * a31 - a11 * a33 + a11 * b33 - a13 * b31 - a31 * b13 + a33 * b11 - b11 * b33 + b13 * b31) * yAB^2
	const3t::Float64 = (3.0 * a11 * a22 - 3.0 * a12 * a21 - 2.0 * a11 * b22 + 2.0 * a12 * b21 + 2.0 * a21 * b12 - 2.0 * a22 * b11 + b11 * b22 - b12 * b21) * zAB^2 + (xAB * (3.0 * a12 * a23 - 3.0 * a13 * a22 + 3.0 * a21 * a32 - 3.0 * a22 * a31 - 2.0 * a12 * b23 + 2.0 * a13 * b22 + 2.0 * a22 * b13 - 2.0 * a23 * b12 - 2.0 * a21 * b32 + 2.0 * a22 * b31 + 2.0 * a31 * b22 - 2.0 * a32 * b21 + b12 * b23 - b13 * b22 + b21 * b32 - b22 * b31) - yAB * (3.0 * a11 * a23 - 3.0 * a13 * a21 + 3.0 * a11 * a32 - 3.0 * a12 * a31 - 2.0 * a11 * b23 + 2.0 * a13 * b21 + 2.0 * a21 * b13 - 2.0 * a23 * b11 - 2.0 * a11 * b32 + 2.0 * a12 * b31 + 2.0 * a31 * b12 - 2.0 * a32 * b11 + b11 * b23 - b13 * b21 + b11 * b32 - b12 * b31)) * zAB + (3.0 * a22 * a33 - 3.0 * a23 * a32 - 2.0 * a22 * b33 + 2.0 * a23 * b32 + 2.0 * a32 * b23 - 2.0 * a33 * b22 + b22 * b33 - b23 * b32) * xAB^2 + (3.0 * a13 * a32 - 3.0 * a12 * a33 - 3.0 * a21 * a33 + 3.0 * a23 * a31 + 2.0 * a12 * b33 - 2.0 * a13 * b32 - 2.0 * a32 * b13 + 2.0 * a33 * b12 + 2.0 * a21 * b33 - 2.0 * a23 * b31 - 2.0 * a31 * b23 + 2.0 * a33 * b21 - b12 * b33 + b13 * b32 - b21 * b33 + b23 * b31) * xAB * yAB + (3.0 * a11 * a33 - 3.0 * a13 * a31 - 2.0 * a11 * b33 + 2.0 * a13 * b31 + 2.0 * a31 * b13 - 2.0 * a33 * b11 + b11 * b33 - b13 * b31) * yAB^2
	const2t::Float64 = (3.0 * a12 * a21 - 3.0 * a11 * a22 + a11 * b22 - a12 * b21 - a21 * b12 + a22 * b11) * zAB^2 + (yAB * (3.0 * a11 * a23 - 3.0 * a13 * a21 + 3.0 * a11 * a32 - 3.0 * a12 * a31 - a11 * b23 + a13 * b21 + a21 * b13 - a23 * b11 - a11 * b32 + a12 * b31 + a31 * b12 - a32 * b11) - xAB * (3.0 * a12 * a23 - 3.0 * a13 * a22 + 3.0 * a21 * a32 - 3.0 * a22 * a31 - a12 * b23 + a13 * b22 + a22 * b13 - a23 * b12 - a21 * b32 + a22 * b31 + a31 * b22 - a32 * b21)) * zAB + (3.0 * a23 * a32 - 3.0 * a22 * a33 + a22 * b33 - a23 * b32 - a32 * b23 + a33 * b22) * xAB^2 + (3.0 * a12 * a33 - 3.0 * a13 * a32 + 3.0 * a21 * a33 - 3.0 * a23 * a31 - a12 * b33 + a13 * b32 + a32 * b13 - a33 * b12 - a21 * b33 + a23 * b31 + a31 * b23 - a33 * b21) * xAB * yAB + (3.0 * a13 * a31 - 3.0 * a11 * a33 + a11 * b33 - a13 * b31 - a31 * b13 + a33 * b11) * yAB^2
	const1t::Float64 = (a11 * a22 - a12 * a21) * zAB^2 + (xAB * (a12 * a23 - a13 * a22 + a21 * a32 - a22 * a31) - yAB * (a11 * a23 - a13 * a21 + a11 * a32 - a12 * a31)) * zAB + (a22 * a33 - a23 * a32) * xAB^2 + (a13 * a32 - a12 * a33 - a21 * a33 + a23 * a31) * xAB * yAB + (a11 * a33 - a13 * a31) * yAB^2
	const3n::Float64 = a11 * a23 * a32 - a11 * a22 * a33 + a12 * a21 * a33 - a12 * a23 * a31 - a13 * a21 * a32 + a13 * a22 * a31 + a11 * a22 * b33 - a11 * a23 * b32 - a11 * a32 * b23 + a11 * a33 * b22 - a12 * a21 * b33 + a12 * a23 * b31 + a12 * a31 * b23 - a12 * a33 * b21 + a13 * a21 * b32 - a13 * a22 * b31 - a13 * a31 * b22 + a13 * a32 * b21 + a21 * a32 * b13 - a21 * a33 * b12 - a22 * a31 * b13 + a22 * a33 * b11 + a23 * a31 * b12 - a23 * a32 * b11 - a11 * b22 * b33 + a11 * b23 * b32 + a12 * b21 * b33 - a12 * b23 * b31 - a13 * b21 * b32 + a13 * b22 * b31 + a21 * b12 * b33 - a21 * b13 * b32 - a22 * b11 * b33 + a22 * b13 * b31 + a23 * b11 * b32 - a23 * b12 * b31 - a31 * b12 * b23 + a31 * b13 * b22 + a32 * b11 * b23 - a32 * b13 * b21 - a33 * b11 * b22 + a33 * b12 * b21 + b11 * b22 * b33 - b11 * b23 * b32 - b12 * b21 * b33 + b12 * b23 * b31 + b13 * b21 * b32 - b13 * b22 * b31
	const2n::Float64 = 3.0 * a11 * a22 * a33 - 3.0 * a11 * a23 * a32 - 3.0 * a12 * a21 * a33 + 3.0 * a12 * a23 * a31 + 3.0 * a13 * a21 * a32 - 3.0 * a13 * a22 * a31 - 2.0 * a11 * a22 * b33 + 2.0 * a11 * a23 * b32 + 2.0 * a11 * a32 * b23 - 2.0 * a11 * a33 * b22 + 2.0 * a12 * a21 * b33 - 2.0 * a12 * a23 * b31 - 2.0 * a12 * a31 * b23 + 2.0 * a12 * a33 * b21 - 2.0 * a13 * a21 * b32 + 2.0 * a13 * a22 * b31 + 2.0 * a13 * a31 * b22 - 2.0 * a13 * a32 * b21 - 2.0 * a21 * a32 * b13 + 2.0 * a21 * a33 * b12 + 2.0 * a22 * a31 * b13 - 2.0 * a22 * a33 * b11 - 2.0 * a23 * a31 * b12 + 2.0 * a23 * a32 * b11 + a11 * b22 * b33 - a11 * b23 * b32 - a12 * b21 * b33 + a12 * b23 * b31 + a13 * b21 * b32 - a13 * b22 * b31 - a21 * b12 * b33 + a21 * b13 * b32 + a22 * b11 * b33 - a22 * b13 * b31 - a23 * b11 * b32 + a23 * b12 * b31 + a31 * b12 * b23 - a31 * b13 * b22 - a32 * b11 * b23 + a32 * b13 * b21 + a33 * b11 * b22 - a33 * b12 * b21
	const1n::Float64 = 3.0 * a11 * a23 * a32 - 3.0 * a11 * a22 * a33 + 3.0 * a12 * a21 * a33 - 3.0 * a12 * a23 * a31 - 3.0 * a13 * a21 * a32 + 3.0 * a13 * a22 * a31 + a11 * a22 * b33 - a11 * a23 * b32 - a11 * a32 * b23 + a11 * a33 * b22 - a12 * a21 * b33 + a12 * a23 * b31 + a12 * a31 * b23 - a12 * a33 * b21 + a13 * a21 * b32 - a13 * a22 * b31 - a13 * a31 * b22 + a13 * a32 * b21 + a21 * a32 * b13 - a21 * a33 * b12 - a22 * a31 * b13 + a22 * a33 * b11 + a23 * a31 * b12 - a23 * a32 * b11

	b = 0.5

	fb = (((const4t * b + const3t) * b + const2t) * b + const1t) * b/(((const3n * b + const2n) * b + const1n) * b)
	if fb >= 1.0
		overlapfun = 1.0
	else
		iter = 1
		a = 0.0
		c = 1.0
		while fb < 1 && iter < number_of_iterations_overlap_criterion
			# Pick point in upper half.
			x = 0.5 * (b + c)
			fx = (((const4t * x + const3t) * x + const2t) * x + const1t) * x/(((const3n * x + const2n) * x + const1n) * x)
			if fx < fb
				c = x
			else
				b = x
				fb = fx
			end
			# Pick point in lower half.
			x = 0.5 * (a + b)
			fx = (((const4t * x + const3t) * x + const2t) * x + const1t) * x/(((const3n * x + const2n) * x + const1n) * x)
		
			if fx < fb
				a = x
			else
				b = x
				fb = fx
			end
			
			iter = iter + 1
		end
		overlapfun = fb
	end
	
	return overlapfun
end