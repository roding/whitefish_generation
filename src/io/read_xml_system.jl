function read_xml_system(file_name::String)
	file_stream::IOStream = open(file_name, "r")
	file_string::String = readstring(file_stream)
	close(file_stream)
		
	# Read simulation domain size in x direction.
	ind_before = search(file_string, "<domain_size_x>")
	ind_after = search(file_string, "</domain_size_x>")
	ind = ind_before[end]+1:ind_after[1]-1
	Lx::Float64 = parse(Float64, file_string[ind])
	
	# Read simulation domain size in y direction.
	ind_before = search(file_string, "<domain_size_y>")
	ind_after = search(file_string, "</domain_size_y>")
	ind = ind_before[end]+1:ind_after[1]-1
	Ly::Float64 = parse(Float64, file_string[ind])
	
	# Read simulation domain size in z direction.
	ind_before = search(file_string, "<domain_size_z>")
	ind_after = search(file_string, "</domain_size_z>")
	ind = ind_before[end]+1:ind_after[1]-1
	Lz::Float64 = parse(Float64, file_string[ind])
	
	# Read particle type.
	ind_before = search(file_string, "<type>")
	ind_after = search(file_string, "</type>")
	ind = ind_before[end]+1:ind_after[1]-1
	particle_category::Symbol = Symbol(file_string[ind])
	
	# Read number of particles.
	ind_before = search(file_string, "<number_of_particles>")
	ind_after = search(file_string, "</number_of_particles>")
	ind = ind_before[end]+1:ind_after[1]-1
	number_of_particles::Int64 = parse(Int64, file_string[ind])
	
	# FOR NOW WE ASSUME PARTICLES ARE ELLIPICAL DISKS.
	
	ind_before = search(file_string, "<X>")
	ind_after = search(file_string, "</X>")
	ind = ind_before[end]+1:ind_after[1]-1
	string_array = split(file_string[ind], ",")
	X = Array(Float64, number_of_particles)
	for current_particle = 1:number_of_particles
		X[current_particle] = parse(Float64, string_array[current_particle])
	end
	
	ind_before = search(file_string, "<Y>")
	ind_after = search(file_string, "</Y>")
	ind = ind_before[end]+1:ind_after[1]-1
	string_array = split(file_string[ind], ",")
	Y = Array(Float64, number_of_particles)
	for current_particle = 1:number_of_particles
		Y[current_particle] = parse(Float64, string_array[current_particle])
	end
	
	ind_before = search(file_string, "<Z>")
	ind_after = search(file_string, "</Z>")
	ind = ind_before[end]+1:ind_after[1]-1
	string_array = split(file_string[ind], ",")
	Z = Array(Float64, number_of_particles)
	for current_particle = 1:number_of_particles
		Z[current_particle] = parse(Float64, string_array[current_particle])
	end
	
	ind_before = search(file_string, "<THETA1>")
	ind_after = search(file_string, "</THETA1>")
	ind = ind_before[end]+1:ind_after[1]-1
	string_array = split(file_string[ind], ",")
	THETA1 = Array(Float64, number_of_particles)
	for current_particle = 1:number_of_particles
		THETA1[current_particle] = parse(Float64, string_array[current_particle])
	end
	
	ind_before = search(file_string, "<THETA2>")
	ind_after = search(file_string, "</THETA2>")
	ind = ind_before[end]+1:ind_after[1]-1
	string_array = split(file_string[ind], ",")
	THETA2 = Array(Float64, number_of_particles)
	for current_particle = 1:number_of_particles
		THETA2[current_particle] = parse(Float64, string_array[current_particle])
	end
	
	ind_before = search(file_string, "<THETA3>")
	ind_after = search(file_string, "</THETA3>")
	ind = ind_before[end]+1:ind_after[1]-1
	string_array = split(file_string[ind], ",")
	THETA3 = Array(Float64, number_of_particles)
	for current_particle = 1:number_of_particles
		THETA3[current_particle] = parse(Float64, string_array[current_particle])
	end
	
	ind_before = search(file_string, "<R1>")
	ind_after = search(file_string, "</R1>")
	ind = ind_before[end]+1:ind_after[1]-1
	string_array = split(file_string[ind], ",")
	R1 = Array(Float64, number_of_particles)
	for current_particle = 1:number_of_particles
		R1[current_particle] = parse(Float64, string_array[current_particle])
	end
	
	ind_before = search(file_string, "<R2>")
	ind_after = search(file_string, "</R2>")
	ind = ind_before[end]+1:ind_after[1]-1
	string_array = split(file_string[ind], ",")
	R2 = Array(Float64, number_of_particles)
	for current_particle = 1:number_of_particles
		R2[current_particle] = parse(Float64, string_array[current_particle])
	end
	
	return (X, Y, Z, THETA1, THETA2, THETA3, R1, R2, Lx, Ly, Lz)
end

#read_xml_system("particle_system.xml")