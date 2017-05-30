function set_version()
	version_major::Int64 = 0
	version_minor::Int64 = 0
	version_patch::Int64 = 0

	file_name_version = "../version"
	file_stream_version = open(file_name_version, "w")
	write(file_stream_version, version_major)
	write(file_stream_version, version_minor)
	write(file_stream_version, version_patch)
	close(file_stream_version)
end

set_version()