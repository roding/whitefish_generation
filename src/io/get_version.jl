function get_version()
	file_name_version = "../version"
	file_stream_version = open(file_name_version, "r")
	version::Array{Int64, 1} = read(file_stream_version, Int64, 3)
	close(file_stream_version)
		
	version_major::Int64 = version[1]
	version_minor::Int64 = version[2]
	version_patch::Int64 = version[3]
	
	version_string::String = join((string(version_major), ".", string(version_minor), ".", string(version_patch)))
	
	return version_string
end