clear
clc
close all hidden

file_path = '../../test_files/output_generation.xml';
[   particle_type, ...
    R, ...
    Lx, ...
    Ly, ...
    Lz, ...
    phi, ...
    X, ...
    Y, ...
    Z, ...
    Q0, ...
    Q1, ...
    Q2, ...
    Q3, ...
    execution_time] = read_output(file_path);

number_of_particles = numel(X);

fig = figure();

h = axes();
hold on

n = 40;
[xsi,eta,zeta] = sphere(n);

count = 0;
for current_particle = 1:number_of_particles
    x = X(current_particle);
    y = Y(current_particle);
    z = Z(current_particle);
    r = R(current_particle);
    
    for i = -1:1
        for j = -1:1
            for k = -1:1
                if (i*Lx + x >= - r) && (i*Lx + x <= Lx + r) && ...
                   (j*Ly + y >= - r) && (j*Ly + y <= Ly + r) && ...
                   (k*Lz + z >= - r) && (k*Lz + z <= Lz + r)
                    hs = surf(i * Lx + x + r * xsi, j * Ly + y + r * eta, k * Lz + z + r * zeta);
                    hs.LineStyle = 'none';
                    hs.SpecularExponent = .7;
                    hs.SpecularStrength = .5;
                    
                    count = count + 1;
                    disp(count)
                end
            end
        end
    end
end
                    
map = repmat([.2 .2 .8],[64 1]);
colormap(map)

h.XLim = [0 Lx];
h.YLim = [0 Ly];
h.ZLim = [0 Lz];

h.XTick = [];
h.YTick = [];
h.ZTick = [];

h.Box = 'on';
h.BoxStyle = 'full';

h.Projection = 'perspective';
h.View = [60, 20];

axis vis3d tight
camlight left; 
lighting flat

axis 'equal'
axis([0 Lx 0 Ly 0 Lz])