load 'real_data/sprottb.mat'

[N, vc] = size(sprottb);
N = N / 2;
h = h * 2;
x = sprottb((1:N) * 2, :);

t = (1:N) * h;

deg = 2;
sigma = deglexord(deg, vc);
mc = size(sigma, 1);

% y = diff4(x, t);
% B = EvalPoly(eye(mc), x, sigma);
% delta = 1e-1;
% H0 = (B'*B + delta*eye(mc))\B'*y;

H0 = orthpoly_sindy(t, x, sigma, 2);

a0 = 1 * ones(1, vc);
p0 = [a0; H0];

[x3, err2] = fminunc(@(params) sprottb_objective(params, t, x, sigma), p0);

plot3(x3(:, 1), x3(:, 2), x3(:, 3));


%% Функции

function err = sprottb_objective(params, t, x, sigma)
    [tspan, y] = sprottb_simulate(params, t, x(1, :), sigma);
    
    err = sqrt(sum(vecnorm(y - x, 2, 2).^2)) / length(tspan);
end

function [t, y] = sprottb_simulate(params, t, x0, sigma)
    h = t(2) - t(1);
    
    alpha = params(1, :);
    H = params(2:end, :);
    
    N = length(t);
    vc = length(x0);
    y = x0 .* ones(N, vc);
    
    for j = 2:N
        Tuntil = t(j);
        for tcur = h:h:Tuntil
            % Euler
            y(j, :) = y(j, :) + h * (Tuntil - tcur).^(alpha - 1)./gamma(alpha).*EvalPoly(H, y(j, :), sigma);
        end
    end
end