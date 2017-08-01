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

% disp('mirrored')
% X = Lx - X;
% Y = Ly - Y;
% Z = Lz - Z;

number_of_particles = numel(X);

xsi = 2 * [0 1 1 0 0 0 ; 1 1 0 0 1 1 ; 1 1 0 0 1 1 ; 0 1 1 0 0 0] - 1;
eta = 2 * [0 0 1 1 0 0 ; 0 1 1 0 0 0 ; 0 1 1 0 1 1 ; 0 0 1 1 1 1] - 1;
zeta = 2 * [0 0 0 0 0 1 ; 0 0 0 0 0 1 ; 1 1 1 1 0 1 ; 1 1 1 1 0 1] - 1;

fig = figure();

h = axes();
hold on

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
    r3 = R(current_particle, 3);
    
    Rq = rotation_matrix(q0, q1, q2, q3);
    XSI = Rq * [r1*xsi(:)' ; r2*eta(:)' ; r3*zeta(:)'];
    
    xsi_sc_rot = reshape(XSI(1,:), [4, 6]);
    eta_sc_rot = reshape(XSI(2,:), [4, 6]);
    zeta_sc_rot = reshape(XSI(3,:), [4, 6]);
    
    rmax = max([r1, r2, r3]);
    for i = -1:1
        for j = -1:1
            for k = -1:1
                if (i*Lx + x >= - rmax) && (i*Lx + x <= Lx + rmax) && ...
                   (j*Ly + y >= - rmax) && (j*Ly + y <= Ly + rmax) && ...
                   (k*Lz + z >= - rmax) && (k*Lz + z <= Lz + rmax)
                    
                   for current_facet = 1:6
                       hp = patch(i * Lx + x + xsi_sc_rot(:, current_facet), j * Ly + y + eta_sc_rot(:, current_facet), k * Lz + z + zeta_sc_rot(:, current_facet), [.5 .5 .5]);
                   end

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
% h.View = [60, 20];

h.View = [-34, 1];
% h.View = [-70, 15];


axis vis3d tight
camlight left; 
lighting flat

axis 'equal'
axis([0 Lx 0 Ly 0 Lz])
% 
check_overlap_cuboid()