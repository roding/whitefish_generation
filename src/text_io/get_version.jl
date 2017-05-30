function get_version()
	cmd::Cmd = `git describe`
	run(cmd)
	
	version_string::String = ""
	
	return version_string
end