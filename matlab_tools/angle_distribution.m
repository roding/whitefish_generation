clear
clc
close all hidden

file_name = '../src/output.dat';
file_data = dlmread(file_name, ',');
Lx = file_data(1, 1);
Ly = file_data(1, 2);
Lz = file_data(1, 3);

R1 = file_data(2:end, 1);
R2 = file_data(2:end, 2);
R3 = file_data(2:end, 3);
X = file_data(2:end, 4);
Y = file_data(2:end, 5);
Z = file_data(2:end, 6);
Q0 = file_data(2:end, 7);
Q1 = file_data(2:end, 8);
Q2 = file_data(2:end, 9);
Q3 = file_data(2:end, 10);

number_of_particles = numel(R1);

angle = zeros(number_of_particles, 1);

v = randn(3, 1);
v = v / norm(v);

for current_particle = 1:number_of_particles
    Rq = rotation_matrix(Q0(current_particle), Q1(current_particle), Q2(current_particle), Q3(current_particle));
    angle(current_particle) = acos( Rq(:, 1)' * v );
end

figure, hold on
alpha = linspace(0, pi, 10000);
P = 1/2 * (1 - cos(alpha));
plot(alpha, P)
[P_emp, x_emp] = ecdf(angle);
plot(x_emp, P_emp)