close all;
warning off;

% System Simulation
sysname = 'SprottB';
sys_f = 'P7_1.txt';

addpath('..\int_diff');

lsm_on = 0;     % LSM comparison
method = "dmt"; % "raw" | "dmt" - raw method or using delMinorTerms
tol = 1;
phase_on = 1;   % Phase portrait
sync_on = 0;    % Synchronization

disp(['Reconstructing ', sysname, ' system, Variant ', sys_f]);
x_r = readmatrix(['research\IGN_', sysname, '/', sys_f]);
h_pred = 6e-6;

h = h_pred; 

t = h * x_r(:, 1);
x_r = x_r(:, 2 * (1:3)) ./ [1 (7.07 * 10000 / 120) 1];

start_point = x_r(1, :); % Initial point
%[~, x_orig] = ode45(@SprottB, t, start_point);

deg = 2; % Degree of reconstructed function
vc = size(x_r, 2); % Variables count

eqc = vc; % Equations count

sigma = deglexord(deg, vc);
mc = size(sigma, 1); % Monomials count
T = mat2cell(repmat(sigma, 1, eqc), mc, repmat(vc, 1, eqc));

x = sgolayfilt(x_r, deg, 2 * deg + 1);

F = orthpoly_F(sigma, t, x, eye(mc), 1);
%F = orthpoly_t4(sigma, t, x, 1);

if method == "raw"
    coefs = zeros(mc, eqc);
    E = EvalPoly(F', x, sigma);
    for i = 1:eqc
        %coefs(:, i) = trapz(x(:, i), E));
        for j = 1:mc
            coefs(j, i) = intdiff4(x(:, i), E(:, j));
        end
    end

    Ho_c = mat2cell(F' * coefs, mc, ones(1, eqc));
elseif method == "dmt"
    y = diff4(x);
    Ho_c = cell(1, vc);
    for j = 1:vc
        [~, ~, ~, Ho_c{1, j}] = delMinorTerms_dy(t, x(:, j), x, y(:, j), F, sigma, tol);
    end
end

Ho_c{1, :}

[~, x1] = ode78(@(t, x)oderecon(Ho_c, T, t, x), t, start_point);

if (lsm_on)
    y = diff4(x, t);
    B = EvalPoly(eye(mc), x, sigma);
    delta = 0;
    
    if method == "raw"
        Ht = (B'*B + delta*eye(mc))\B'*y;
        Ht_c = mat2cell(Ht, mc, ones(1, eqc));
    elseif method == "dmt"
        Ht_c = cell(1, vc);
        for j = 1:vc
            [Ht_c{1, j}, ~] = delMinorTerms(x, y(:, j), sigma, tol);
        end
    end
    [~, x2] = ode78(@(t, x)oderecon(Ht_c, T, t, x), t, start_point);
end

%x = sgolayfilt(x, 1, 7);

if (phase_on)
    figure(1);
    plot3(x(:, 1), x(:, 2), x(:, 3), 'b', 'DisplayName', 'Smoothed original data');
    hold on; grid on;
    plot3(x1(:, 1), x1(:, 2), x1(:, 3), 'r', 'DisplayName', 'Orthogonal Polynomials reconstruction');
    if (lsm_on)
        %plot3(x2(:, 1), x2(:, 2), x2(:, 3), 'g', 'DisplayName', 'LSM');
        plot3(x2(:, 1), x2(:, 2), x2(:, 3), 'r', 'DisplayName', 'LSM');
    end
    xtickformat('$%g$'); ytickformat('$%g$'); ztickformat('$%g$')
    set(gca,'TickLabelInterpreter','latex');
    xlabel('$x$','Interpreter','latex');
    ylabel('$y$','Interpreter','latex');
    zlabel('$z$','Interpreter','latex');
    legend show
end

if (sync_on)
    sync_coef = 10e4;
    %xo_slave = RK4SyncFOH(x, Ho_c, T, t, sync_coef);
    xo_slave = RK4sync(x, Ho_c, T, t, sync_coef);
    
    figure(2);
    err_orth = vecnorm(x - xo_slave, 2, 2);
    semilogy(t, err_orth, 'Color', [0 0 1 0.2], 'LineWidth', 0.75);
    hold on; grid on;
    mean_orth = mean(err_orth)
    semilogy(t, (t * 0) + mean_orth, 'Color', [0 0 1 1], 'LineWidth', 1);
    legend('', 'Orthpoly (Mean)');
    
    if (lsm_on)
        %xt_slave = RK4SyncFOH(x, Ht_c, T, t, sync_coef);
        xt_slave = RK4sync(x, Ht_c, T, t, sync_coef);

        figure(2);
        err_lsm = vecnorm(x - xt_slave, 2, 2);
        semilogy(t, err_lsm, 'Color', [1 0 0 0.2], 'LineWidth', 0.75);
        mean_lsm = mean(err_lsm)
        semilogy(t, (t * 0) + mean_lsm, 'Color', [1 0 0 1], 'LineWidth', 1);
        legend('', 'Orthpoly (Mean)', '', 'LSM (Mean)');
    end

    figure(2);
    legend show;
    xtickformat('$%g$'); ytickformat('$%g$'); ztickformat('$%g$');
    set(gca, 'TickLabelInterpreter', 'latex');
    %xlabel('Time $t$, s', 'Interpreter', 'latex');
    ylabel('Synchronization error $\overline{\zeta}$', 'Interpreter', 'latex');
    xlim([t(end - 1000), t(end)]);
end