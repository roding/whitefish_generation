function print_progress_diffusion(t_elapsed_diffusion::Float64, current_diffuser::Int64, number_of_diffusers::Int64)
	time_elapsed_seconds::Int64 = convert(Int64, floor(t_elapsed_diffusion))
	
	time_elapsed_hours::Int64 = fld(time_elapsed_seconds, 3600)
	time_elapsed_seconds = time_elapsed_seconds - time_elapsed_hours * 3600
	
	time_elapsed_minutes::Int64 = fld(time_elapsed_seconds, 60)
	time_elapsed_seconds = time_elapsed_seconds - time_elapsed_minutes * 60
	
	time_elapsed_hours_string::String = string(time_elapsed_hours)
	if length(time_elapsed_hours_string) == 1
		time_elapsed_hours_string = join(("0", time_elapsed_hours_string))
	end
	
	time_elapsed_minutes_string::String = string(time_elapsed_minutes)
	if length(time_elapsed_minutes_string) == 1
		time_elapsed_minutes_string = join(("0", time_elapsed_minutes_string))
	end
	
	time_elapsed_seconds_string::String = string(time_elapsed_seconds)
	if length(time_elapsed_seconds_string) == 1
		time_elapsed_seconds_string = join(("0", time_elapsed_seconds_string))
	end
	
	time_elapsed_string::String = join((time_elapsed_hours_string, ":", time_elapsed_minutes_string, ":", time_elapsed_seconds_string))
	
	fraction_done::Float64 = convert(Float64, current_diffuser) / convert(Float64, number_of_diffusers)
	
	time_remaining_seconds::Int64 = convert(Int64, round((1 - fraction_done) / fraction_done * t_elapsed_diffusion))
	
	time_remaining_hours::Int64 = fld(time_remaining_seconds, 3600)
	time_remaining_seconds = time_remaining_seconds - time_remaining_hours * 3600
	
	time_remaining_minutes::Int64 = fld(time_remaining_seconds, 60)
	time_remaining_seconds = time_remaining_seconds - time_remaining_minutes * 60
	
	time_remaining_hours_string::String = string(time_remaining_hours)
	if length(time_remaining_hours_string) == 1
		time_remaining_hours_string = join(("0", time_remaining_hours_string))
	end
	
	time_remaining_minutes_string::String = string(time_remaining_minutes)
	if length(time_remaining_minutes_string) == 1
		time_remaining_minutes_string = join(("0", time_remaining_minutes_string))
	end
	
	time_remaining_seconds_string::String = string(time_remaining_seconds)
	if length(time_remaining_seconds_string) == 1
		time_remaining_seconds_string = join(("0", time_remaining_seconds_string))
	end
	
	time_remaining_string::String = join((time_remaining_hours_string, ":", time_remaining_minutes_string, ":", time_remaining_seconds_string))
	
	percent_done_string::String = @sprintf("%2.2f", 100.0 * fraction_done)
	if length(percent_done_string) == 4
		percent_done_string = join(("00", percent_done_string))
	elseif length(percent_done_string) == 5
		percent_done_string = join(("0", percent_done_string))
	end
	
	output_str::String = join(("   ", time_elapsed_string, "                  ", percent_done_string, "     ", time_remaining_string))
	println(output_str)
	
	nothing
	#percent_done::Float64 = 100.0 * fraction_done
	
			
end
			