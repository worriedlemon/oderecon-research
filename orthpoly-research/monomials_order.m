rng_i default;
warning off;
close all

% Used system
%system = @Rossler
system = @Lorenz

[Href, deg, vc] = load_href(func2str(system)); % Used coefficients
Ttrans = 100; % Transient time
Tmax = 5; % Time of experimental series for Lorenz
%Tmax = 10; % Time of experimental series for Rossler
h = 1e-2; % Step
start_point = [4 -2 0]; % Initial point

% %first, go transient
[~, x] = ode78(system, 0:h:Ttrans, start_point);
start_point = x(end, :);
%then, go experiment
[t, x] = ode78(system, 0:h:Tmax, start_point);


y = transpose(system(0, x'));
N = size(x, 1);
eqv = vc;

sigma = deglexord(deg, vc);
mc = size(sigma, 1);
revidx = mc:-1:1;
randidx = randperm(mc);
sigma_rev = sigma(revidx, :);
sigma_rand = sigma(randidx, :);

Href_rev = Href(revidx, :);
Href_rand = Href(randidx, :);

noises = logspace(-5, 1, 70);
erro = zeros(length(noises), 3);
tol = 1e-2;
for k = 1:length(noises)
    noise_amp = noises(k);
    rx = x + noise_amp * randn(N, vc);
    ry = diff4(rx, t); % fourth order
    
    F_dir = orthpoly_t(sigma, t, rx);
    F_rev = orthpoly_t(sigma_rev, t, rx);
    F_rand = orthpoly_t(sigma_rand, t, rx);
    
    Ho_dir = zeros(mc, eqc);
    Ho_rev = Ho_dir;
    Ho_rand = Ho_dir;

    for j = 1:eqc
        [~, ~, ~, Ho_dir(:, j)] = delMinorTerms_dy(t, x(:, j), x, y(:, j), F_dir, sigma, tol, 0);
        [~, ~, ~, Ho_rev(:, j)] = delMinorTerms_dy(t, x(:, j), x, y(:, j), F_rev, sigma_rev, tol, 0);
        [~, ~, ~, Ho_rand(:, j)] = delMinorTerms_dy(t, x(:, j), x, y(:, j), F_rand, sigma_rand, tol, 0);
    end
    
    % E_dir = EvalPoly(F_dir', rx, sigma);
    % E_rev = EvalPoly(F_rev', rx, sigma_rev);
    % E_rand = EvalPoly(F_rand', rx, sigma_rand);
    % for j = 1:mc
    %     for i = 1:eqc
    %         Ho_dir(j, i) = intdiff4(rx(:, i), E_dir(:, j));
    %         Ho_rev(j, i) = intdiff4(rx(:, i), E_rev(:, j));
    %         Ho_rand(j, i) = intdiff4(rx(:, i), E_rand(:, j));
    %     end
    % end
    % 
    % Ho_dir =  F_dir' * Ho_dir;
    % Ho_rev =  F_rev' * Ho_rev;
    % Ho_rand =  F_rand' * Ho_rand;

    erro(k, 1) = norm(Ho_dir - Href);
    erro(k, 2) = norm(Ho_rev - Href_rev);
    erro(k, 3) = norm(Ho_rand - Href_rand);
end

figure(1);
loglog(noises, erro(:, 1), 'b-', Marker='.', DisplayName='Direct order');
hold on; grid on;
loglog(noises, erro(:, 2), 'r--', Marker='o', DisplayName='Reverse order');
loglog(noises, erro(:, 3), 'g-.', Marker='p', DisplayName='Random order');
legend show;
xlabel('Noise magnitude $\sigma$', Interpreter='latex')
ylabel('Coefficients error $\zeta$', Interpreter='latex')
set(gca, 'TickLabelInterpreter', 'latex')
title('Monomials order sensitivity')