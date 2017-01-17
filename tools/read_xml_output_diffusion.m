clear
close all hidden
clc

folders = dir('output');
folders = folders(3:end); % Remove . and ..
number_of_folders = numel(folders);

for current_folder = 1:number_of_folders
    str = fileread(['output' '/' folders(current_folder).name '/' 'diffusion.xml']);

    pattern_before = '<diagnostic_diffusion_coefficient_ratio>';
    pattern_after = [pattern_before(1) '/' pattern_before(2:end)];
    ind_start = strfind(str, pattern_before) + numel(pattern_before);
    ind_end = strfind(str, pattern_after) - 1;
    diagnostic_diffusion_coefficient_ratio{current_folder} = str2double(str(ind_start:ind_end));

    pattern_before = '<time>';
    pattern_after = [pattern_before(1) '/' pattern_before(2:end)];
    ind_start = strfind(str, pattern_before) + numel(pattern_before);
    ind_end = strfind(str, pattern_after) - 1;
    substr = str(ind_start:ind_end);
    ind_separator = strfind(substr, ',');
    t{current_folder} = zeros(numel(ind_separator)+1, 1);
    ind_separator = [0 ind_separator numel(substr)+1];
    for i = 1:numel(ind_separator)-1
        t{current_folder}(i) = str2double(substr(ind_separator(i)+1:ind_separator(i+1)-1));
    end

    pattern_before = '<mean_square_displacement_x>';
    pattern_after = [pattern_before(1) '/' pattern_before(2:end)];
    ind_start = strfind(str, pattern_before) + numel(pattern_before);
    ind_end = strfind(str, pattern_after) - 1;
    substr = str(ind_start:ind_end);
    ind_separator = strfind(substr, ',');
    msd_x{current_folder} = zeros(numel(ind_separator)+1, 1);
    ind_separator = [0 ind_separator numel(substr)+1];
    for i = 1:numel(ind_separator)-1
        msd_x{current_folder}(i) = str2double(substr(ind_separator(i)+1:ind_separator(i+1)-1));
    end

    pattern_before = '<mean_square_displacement_y>';
    pattern_after = [pattern_before(1) '/' pattern_before(2:end)];
    ind_start = strfind(str, pattern_before) + numel(pattern_before);
    ind_end = strfind(str, pattern_after) - 1;
    substr = str(ind_start:ind_end);
    ind_separator = strfind(substr, ',');
    msd_y{current_folder} = zeros(numel(ind_separator)+1, 1);
    ind_separator = [0 ind_separator numel(substr)+1];
    for i = 1:numel(ind_separator)-1
        msd_y{current_folder}(i) = str2double(substr(ind_separator(i)+1:ind_separator(i+1)-1));
    end

    pattern_before = '<mean_square_displacement_z>';
    pattern_after = [pattern_before(1) '/' pattern_before(2:end)];
    ind_start = strfind(str, pattern_before) + numel(pattern_before);
    ind_end = strfind(str, pattern_after) - 1;
    substr = str(ind_start:ind_end);
    ind_separator = strfind(substr, ',');
    msd_z{current_folder} = zeros(numel(ind_separator)+1, 1);
    ind_separator = [0 ind_separator numel(substr)+1];
    for i = 1:numel(ind_separator)-1
        msd_z{current_folder}(i) = str2double(substr(ind_separator(i)+1:ind_separator(i+1)-1));
    end
    
    figure
    Deff = msd_z{current_folder}./(2.*t{current_folder}); 
    Deff(1) = 1; 
    plot(t{current_folder}, Deff)
end
    