function get_version()
	cmd::Cmd = `git describe --always`
	
	version_string::String = readstring(cmd)
	version_string = version_string[1:end-2]
	
	return version_string
end