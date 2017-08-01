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
    q0 = Q0(current_particle);
    q1 = Q1(current_particle);
    q2 = Q2(current_particle);
    q3 = Q3(current_particle);
    r1 = R(current_particle, 1);
    r2 = R(current_particle, 2);
    
    Rq = rotation_matrix(q0, q1, q2, q3);
    XSI = Rq * [r1*xsi(:)' ; r2*eta(:)' ; 0.01*zeta(:)'];
    
    xsi_sc_rot = reshape(XSI(1,:), [n+1 n+1]);
    eta_sc_rot = reshape(XSI(2,:), [n+1 n+1]);
    zeta_sc_rot = reshape(XSI(3,:), [n+1 n+1]);
    
    rmax = max(r1, r2);
    for i = -1:1
        for j = -1:1
            for k = -1:1
                if (i*Lx + x >= - rmax) && (i*Lx + x <= Lx + rmax) && ...
                   (j*Ly + y >= - rmax) && (j*Ly + y <= Ly + rmax) && ...
                   (k*Lz + z >= - rmax) && (k*Lz + z <= Lz + rmax)
                    hs = surf(i*Lx + x + xsi_sc_rot, j*Ly + y + eta_sc_rot, k*Lz + z + zeta_sc_rot);
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