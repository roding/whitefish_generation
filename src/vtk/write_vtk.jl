function write_vtk(M, file_path::String)

	file_stream::IOStream = open(file_path, "w")
	@printf(file_stream, "%s", "# vtk DataFile Version 3.0\n")
	@printf(file_stream, "%s", "Generated from within Julia\n")
	@printf(file_stream, "%s", "ASCII\n")
	@printf(file_stream, "%s", "DATASET STRUCTURED_POINTS\n")
	@printf(file_stream, "%s", join(("DIMENSIONS ", string(size(M, 1)), " ", string(size(M, 2)), " ", string(size(M, 3)), "\n")))
	@printf(file_stream, "%s", "ORIGIN 0 0 0\n")
	@printf(file_stream, "%s", "SPACING 1 1 1\n")
	@printf(file_stream, "%s", join(("POINT_DATA ", string(prod(size(M))), "\n")))
	@printf(file_stream, "%s", "SCALARS values float 1\n")
	@printf(file_stream, "%s", "LOOKUP_TABLE default\n")
	for i = 1:length(M)
		@printf(file_stream, "%s\n", M[i])
	end
	close(file_stream)

	nothing
end
