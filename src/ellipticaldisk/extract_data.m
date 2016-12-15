clear
clc
close all hidden

filename = 'diffusion.dat';
delimiter = ',';
formatSpec = '%f%f%f%f%[^\n\r]';
fileID = fopen(filename,'r');
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'EmptyValue' ,NaN, 'ReturnOnError', false);
fclose(fileID);
diffusion = [dataArray{1:end-1}];
clearvars filename delimiter formatSpec fileID dataArray ans;

t = diffusion(3:end, 1);
msd_x = diffusion(3:end, 2);
msd_y = diffusion(3:end, 3);
msd_z = diffusion(3:end, 4);
D0 = 1;
D_x = ones(size(t));
D_x(2:end) = msd_x(2:end) ./ (2 * t(2:end));
D_y = ones(size(t));
D_y(2:end) = msd_y(2:end) ./ (2 * t(2:end));
D_z = ones(size(t));
D_z(2:end) = msd_z(2:end) ./ (2 * t(2:end));

figure, hold on
plot(t, msd_x)
plot(t, msd_y)
plot(t, msd_z)
title('MSD')
legend('x','y','z')

figure, hold on
plot(t, D_x)
plot(t, D_y)
plot(t, D_z)
title('D/D0')
legend('x','y','z')



