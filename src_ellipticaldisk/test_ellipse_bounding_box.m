clear
clc
close all hidden

seed = sum(1e6*clock())%1;
s = RandStream('mt19937ar', 'Seed', seed);
RandStream.setGlobalStream(s);

xc = rand();
yc = rand();
zc = rand();

r1 = 1; 
r2 = 2;

% theta1 = 0%1e-9*randn();
% theta2 = 0;
% theta3 = randn();

% theta1 = randsample(0:pi/4:2*pi, 1)
% theta2 = randsample(0:pi/4:2*pi, 1)
% theta3 = randsample(0:pi/4:2*pi, 1)

theta1 = 2*pi*rand();
theta2 = 2*pi*rand();
theta3 = 2*pi*rand(); 

t = linspace(0,2*pi,100000);

x = xc + r1*cos(t)*cos(theta2)*cos(theta3) - r2*cos(theta2)*sin(t)*sin(theta3);
y = yc + r1*cos(t)*(cos(theta1)*sin(theta3) + cos(theta3)*sin(theta1)*sin(theta2)) + r2*sin(t)*(cos(theta1)*cos(theta3) - sin(theta1)*sin(theta2)*sin(theta3));
z = zc + r1*cos(t)*(sin(theta1)*sin(theta3) - cos(theta1)*cos(theta3)*sin(theta2)) + r2*sin(t)*(cos(theta3)*sin(theta1) + cos(theta1)*sin(theta2)*sin(theta3));

% figure 
% hold on
% plot3(x, y, z)
% xlabel('x')
% ylabel('y')
% zlabel('z')

t_extrema(1) = atan(-r2/r1*tan(theta3));
t_extrema(2) = t_extrema(1) + pi;
x_extrema = xc + r1*cos(t_extrema)*cos(theta2)*cos(theta3) - r2*cos(theta2)*sin(t_extrema)*sin(theta3);
x_min = min(x_extrema);
x_max = max(x_extrema);

disp([min(x) x_min max(x) x_max])

t_extrema(1) = atan(r2/r1 * (cos(theta1)*cos(theta3) - sin(theta1)*sin(theta2)*sin(theta3)) / (cos(theta1)*sin(theta3) + cos(theta3)*sin(theta1)*sin(theta2)));
t_extrema(2) = t_extrema(1) + pi;
y_extrema = yc + r1*cos(t_extrema)*(cos(theta1)*sin(theta3) + cos(theta3)*sin(theta1)*sin(theta2)) + r2*sin(t_extrema)*(cos(theta1)*cos(theta3) - sin(theta1)*sin(theta2)*sin(theta3));
y_min = min(y_extrema);
y_max = max(y_extrema);

disp([min(y) y_min max(y) y_max])

if isnan((cos(theta3)*sin(theta1) + cos(theta1)*sin(theta2)*sin(theta3)) / (sin(theta1)*sin(theta3) - cos(theta1)*cos(theta3)*sin(theta2)))
    t_extrema(1) = atan(r2/r1 * 1/tan(theta3));
else
    t_extrema(1) = atan(r2/r1 * (cos(theta3)*sin(theta1) + cos(theta1)*sin(theta2)*sin(theta3)) / (sin(theta1)*sin(theta3) - cos(theta1)*cos(theta3)*sin(theta2)));
end
t_extrema(2) = t_extrema(1) + pi;
z_extrema = zc + r1*cos(t_extrema)*(sin(theta1)*sin(theta3) - cos(theta1)*cos(theta3)*sin(theta2)) + r2*sin(t_extrema)*(cos(theta3)*sin(theta1) + cos(theta1)*sin(theta2)*sin(theta3));
z_min = min(z_extrema);
z_max = max(z_extrema);

disp([min(z) z_min max(z) z_max])

norm([min(x)-x_min, max(x)-x_max, min(y)-y_min, max(y)-y_max, min(z)-z_min, max(z)-z_max])


