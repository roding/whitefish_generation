clear
clc
close all hidden

figure, hold on

alpha = linspace(0,pi,10000);
    angle_to_z_axis = acos(Rq(:, 1) * );
P = 1/2*(1 - cos(alpha));
P_lim = 1/2*(1 - cos(pi));
plot(alpha, P / P_lim)
[P_emp, x_emp] = ecdf(angle_to_z_axis);
plot(x_emp, P_emp)


% angle_to_x_axis = acos((cos(THETA2).*cos(THETA3))./(abs(sin(THETA3)).^2 + abs(cos(THETA2).*cos(THETA3)).^2 + abs(cos(THETA3).*sin(THETA2)).^2).^(1/2));
% angle_to_y_axis = acos((cos(THETA3).*cos(THETA1))./(abs(sin(THETA1)).^2 + abs(cos(THETA3).*cos(THETA1)).^2 + abs(cos(THETA1).*sin(THETA3)).^2).^(1/2));
% figure, hist(angle_to_x_axis,100)
% figure, hist(angle_to_y_axis,100)
% figure, hist(angle_to_z_axis,100)