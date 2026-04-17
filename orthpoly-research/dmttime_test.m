close all;
warning off;

% System Simulation
sys = @Lorenz

start_point = [4 -2 0]; % Initial point
Tmax = 100; % Time end
h = 1e-2; % Step

mcs = [];
el_time_lsm = [];
el_time_orth = [];

[t, x] = ode45(sys, 0:h:Tmax, start_point);
y = diff4(x,t);
tol = 1e-3;

for deg = 2:6 % Degree
    disp(deg);
    vc = size(x, 2); % Variables count
    eqc = vc; % Equations count
    sigma = deglexord(deg, vc);
    mc = size(sigma, 1); % Monomials count
    mcs = [mcs mc];
    Ho = cell(1, vc);
    Ht = cell(1, vc);
    
    tic
    F = orthpoly_F(sigma, t, x, eye(mc), 1); % Getting orthogonal polynomials matrix
    
    for j = 1:vc
        [~, ~, ~, Ho{1, j}] = delMinorTerms_dy(t, x(:, j), x, y(:, j), F, sigma, tol);
    end
    t_end = toc;
    el_time_orth = [el_time_orth t_end];
    disp('Orthpoly calculation ended')

    tic
    for j = 1:vc
        [Ht{1, j}, ~] = delMinorTerms(x, y(:, j), sigma, tol);
    end
    t_end = toc;
    el_time_lsm = [el_time_lsm t_end];
    disp('LSM calculation ended')
end

figure(1);
plot(mcs, el_time_lsm, 'DisplayName', 'LSM'); hold on; grid on;
plot(mcs, el_time_orth, 'DisplayName', 'Orthogonal polynomials');
xtickformat('$%g$'); ytickformat('$%g$');
set(gca,'TickLabelInterpreter','latex');
xticks(mcs);
xlabel('Monomials count $L$','Interpreter','latex');
ylabel('Time, s','Interpreter','latex');
legend show
