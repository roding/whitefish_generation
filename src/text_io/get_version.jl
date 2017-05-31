function get_version()
	cmd::Cmd = `git describe --always`
	version_string::String = readstring(cmd)
	
	if version_string[1] == 'v'
		version_string = version_string[2:end]
	end
	
	if version_string[end] == '\n'
		version_string = version_string[1:end-1]
	end
	
	return version_string
end