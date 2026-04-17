function x = mylsm(a, b)
    at = mytranspose(a);
    U = mymult(at, a);
    delta = 0.01; % regularization parameter
    U = U + delta * eye(size(U));
    V = mymult(at, b);
    x = V;
    
    [rows, ~] = size(U);
    for eq = 1:size(V, 2)
        v = V(:, eq);
        u = U;
        for var = 1:rows
            good = true;
            if (u(var, var) == 0)
                good = false;
                for row = 1:rows
                    if (u(row, var) ~= 0)
                        tmp = u(var, :);
                        u(var, :) = u(row, :);
                        u(row, :) = tmp;
    
                        tmp = v(var);
                        v(var) = u(row);
                        v(row) = tmp;
    
                        good = true;
                        break;
                    end
                end
            end
    
            if ~good
               error('No solution');
            end
    
            coef = u(var, var);
            u(var, :) = u(var, :) ./ coef;
            v(var) = v(var) / coef;
    
            for i = 1:rows
                if (i == var)
                    continue;
                end
    
                coef = u(i, var);
                u(i, :) = u(i, :) - u(var, :) * coef;
                v(i) = v(i) - v(var) * coef;
            end
        end
        x(:, eq) = v;
    end
end