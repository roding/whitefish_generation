
for i = 1:number_of_particles
    for j = [1:i-1 i+1:number_of_particles]
%         disp([i j])
        RA = rotation_matrix(Q0(i), Q1(i), Q2(i), Q3(i));
        PA = 2 * [0 0 0 0 1 1 1 1 ; 0 0 1 1 0 0 1 1 ; 0 1 0 1 0 1 0 1] - 1;
        PA(1, :) = R(i, 1) * PA(1, :);
        PA(2, :) = R(i, 2) * PA(2, :);
        PA(3, :) = R(i, 3) * PA(3, :);
        PA = RA * PA;
        
        xAB = X(j) - X(i);
        if xAB < -0.5 * Lx
            xAB = xAB + Lx;
        elseif xAB > 0.5 * Lx
            xAB = xAB - Lx;
        end
        yAB = Y(j) - Y(i);
        if yAB < -0.5 * Ly
            yAB = yAB + Ly;
        elseif yAB > 0.5 * Ly
            yAB = yAB - Ly;
        end
        zAB = Z(j) - Z(i);
        if zAB < -0.5 * Lz
            zAB = zAB + Lz;
        elseif zAB > 0.5 * Lz
            zAB = zAB - Lz;
        end
%         disp([xAB, yAB, zAB])
        RB = rotation_matrix(Q0(j), Q1(j), Q2(j), Q3(j));
        PB = 2 * [0 0 0 0 1 1 1 1 ; 0 0 1 1 0 0 1 1 ; 0 1 0 1 0 1 0 1] - 1;
        PB(1, :) = R(j, 1) * PB(1, :);
        PB(2, :) = R(j, 2) * PB(2, :);
        PB(3, :) = R(j, 3) * PB(3, :);
        PB = RB * PB;
        PB(1, :) = xAB + PB(1, :);
        PB(2, :) = yAB + PB(2, :);
        PB(3, :) = zAB + PB(3, :);

        % First axis.
        is_intersecting = true;
        for k = 1:3
            cA = RA(:, k)' * PA;
            cB = RA(:, k)' * PB;
            if min(cA) >= max(cB) || max(cA) <= min(cB)
        %         disp('yeay')
                is_intersecting = false;
            else
%                 min(cA) - max(cB)
%                 max(cA) - min(cB)
            end
        end

        for k = 1:3
            cA = RB(:, k)' * PA;
            cB = RB(:, k)' * PB;
            if min(cA) >= max(cB) || max(cA) <= min(cB)
        %         disp('yeay')
                is_intersecting = false;
            else
%                 min(cA) - max(cB)
%                 max(cA) - min(cB)
            end
        end

        for k = 1:3
            for l = 1:3
                cross_prod = cross(RA(:, k), RB(:, l));
                cross_prod = cross_prod / sqrt(cross_prod(1)^2 + cross_prod(2)^2 + cross_prod(3)^2);
                cA = cross_prod' * PA;
                cB = cross_prod' * PB;
                if min(cA) >= max(cB) || max(cA) <= min(cB)
        %             disp('yeay')
                    is_intersecting = false;
                else
%                     min(cA) - max(cB)
%                     max(cA) - min(cB)
                end
            end
        end

        if is_intersecting
            disp([i j])
        end
    end
end