clear
close all hidden
clc

main_folder = '/home/sms/output_layer';
folders = dir(main_folder);
folders = folders(3:end); % Remove . and ..
number_of_folders = numel(folders);

is_failed = false(number_of_folders, 1);

for current_folder = 1:number_of_folders
    disp(current_folder)
    try
        file_path = [main_folder '/' folders(current_folder).name '/' 'input_generation.xml'];
        [Lx_, Ly_, Lz_, particle_type_, number_of_particles_, lbz_, ubz_, ubangle_, R1_, R2_, number_of_equilibration_sweeps_, output_generation_path_] = read_xml_input_generation(file_path);
        file_path = [main_folder '/' folders(current_folder).name '/' 'output_generation.xml'];
        [Lx_, Ly_, Lz_, particle_type_, number_of_particles_, X_, Y_, Z_, THETA1_, THETA2_, THETA3_, R1_, R2_, t_exec_generation_] = read_xml_output_generation(file_path);
        file_path = [main_folder '/' folders(current_folder).name '/' 'input_diffusion.xml'];
        [output_generation_path_, inherent_diffusion_coefficient_, deltat_coarse_, number_of_time_points_coarse_, number_of_time_points_fine_per_coarse_, number_of_diffusers_, number_of_cells_x_, number_of_cells_y_, number_of_cells_z_, output_diffusion_path_] = read_xml_input_diffusion(file_path);
        file_path = [main_folder '/' folders(current_folder).name '/' 'output_diffusion.xml'];
        [diagnostic_diffusion_coefficient_ratio_, t_, msd_x_, msd_y_, msd_z_, t_exec_diffusion_] = read_xml_output_diffusion(file_path);

        lbz{current_folder} = lbz_;
        ubz{current_folder} = ubz_;
        ubangle{current_folder} = ubangle_;

        diagnostic_diffusion_coefficient_ratio{current_folder} = diagnostic_diffusion_coefficient_ratio_;
        t{current_folder} = t_;
        msd_x{current_folder} = msd_x_;
        msd_y{current_folder} = msd_y_;
        msd_z{current_folder} = msd_z_;

        t_exec_generation{current_folder} = t_exec_generation_;
        t_exec_diffusion{current_folder} = t_exec_diffusion_;
    catch
        warning('Reading failure');
        is_failed(current_folder) = true;
    end
end

lbz = lbz(~is_failed);
ubz = ubz(~is_failed);
ubangle = ubangle(~is_failed);
diagnostic_diffusion_coefficient_ratio = diagnostic_diffusion_coefficient_ratio(~is_failed);
t = t(~is_failed);
msd_x = msd_x(~is_failed);
msd_y = msd_y(~is_failed);
msd_z = msd_z(~is_failed);
t_exec_generation = t_exec_generation(~is_failed);
t_exec_diffusion = t_exec_diffusion(~is_failed);

clear *_ main_folder folders current_folder number_of_folders file_path is_failed

save('consolidated_data.mat')

    