
for i = 1:number_of_particles
    for j = i+1:number_of_particles
%         disp([i j])
        RA = rotation_matrix(Q0(i), Q1(i), Q2(i), Q3(i));
        PA = 2 * [0 0 0 0 1 1 1 1 ; 0 0 1 1 0 0 1 1 ; 0 1 0 1 0 1 0 1] - 1;
        PA(1, :) = R1(i) * PA(1, :);
        PA(2, :) = R2(i) * PA(2, :);
        PA(3, :) = R3(i) * PA(3, :);
        PA = RA * PA;
        PA(1, :) = X(i) + PA(1, :);
        PA(2, :) = Y(i) + PA(2, :);
        PA(3, :) = Z(i) + PA(3, :);

        RB = rotation_matrix(Q0(j), Q1(j), Q2(j), Q3(j));
        PB = 2 * [0 0 0 0 1 1 1 1 ; 0 0 1 1 0 0 1 1 ; 0 1 0 1 0 1 0 1] - 1;
        PB(1, :) = R1(j) * PB(1, :);
        PB(2, :) = R2(j) * PB(2, :);
        PB(3, :) = R3(j) * PB(3, :);
        PB = RB * PB;
        PB(1, :) = X(j) + PB(1, :);
        PB(2, :) = Y(j) + PB(2, :);
        PB(3, :) = Z(j) + PB(3, :);

        % First axis.
        is_intersecting = true;
        for k = 1:3
            cA = RA(:, k)' * PA;
            cB = RA(:, k)' * PB;
            if min(cA) > max(cB) || max(cA) < min(cB)
        %         disp('yeay')
                is_intersecting = false;
            else
        %         disp('nooo')
            end
        end

        for k = 1:3
            cA = RB(:, k)' * PA;
            cB = RB(:, k)' * PB;
            if min(cA) > max(cB) || max(cA) < min(cB)
        %         disp('yeay')
                is_intersecting = false;
            else
        %         disp('nooo')
            end
        end

        for k = 1:3
            for l = 1:3
                cross_prod = cross(RA(:, k), RB(:, l));
                cross_prod = cross_prod / sqrt(cross_prod(1)^2 + cross_prod(2)^2 + cross_prod(3)^2);
                cA = cross_prod' * PA;
                cB = cross_prod' * PB;
                if min(cA) > max(cB) || max(cA) < min(cB)
        %             disp('yeay')
                    is_intersecting = false;
                else
        %             disp('nooo')
                end
            end
        end

        if is_intersecting
            disp([i j])
        end
    end
end