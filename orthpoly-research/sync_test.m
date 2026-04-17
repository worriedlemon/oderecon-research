rng_i default;
warning off;

% Initial coefficients and points
Plorenz = [0.67963 1.29782 7.83735];
Prossler = [4.57073 -3.0225 0.13198];

% Used system
system = @Rossler

sysname = func2str(system);
[Href, deg, vc] = load_href(sysname); % Used coefficients
Pref = eval(['P', lower(sysname)]);

Tmax = 250; % Time end
h = 1e-2; % Step

[t, x] = ode45(system, 0:h:Tmax, Pref);

N = size(x, 1);
eqc = vc;

sigma = deglexord(deg, vc);
mc = size(sigma, 1);

rx = x;
ry = [diff(rx) / h; (rx(end, :) - rx(end - 1, :)) / h]; % first order
%ry = diff2(rx) / h; % second order
%ry = diff4(rx, t); % fourth order

delta = 0.01; % regularization parameter

B = EvalPoly(eye(mc), rx, sigma);
Ht = (B'*B + delta*eye(mc))\B'*ry;

F = orthpoly_t(sigma, t, rx);

Ho = zeros(mc, eqc);
E = EvalPoly(F', rx, sigma);
for i = 1:eqc
    Ho(:, i) = trapz(rx(:, i), E);
end
Ho = F' * Ho;

Ho_c = mat2cell(Ho, mc, ones(1, eqc));
Ht_c = mat2cell(Ht, mc, ones(1, eqc));
T = mat2cell(repmat(sigma, 1, eqc), mc, repmat(vc, 1, eqc));

t = 0:h:Tmax;

[~, x] = ode45(system, t, Pref);
sync_coef = 1.5;

% xt_slave = [Pref; zeros(N - 1, vc)];
% xo_slave = xt_slave;
% N = length(t);
% for i = 1:N - 1
%     % Runge-Kutta (2-nd order)
%     syncv = sync_coef * [1 1 1];
%     k1 = h / 2 * (EvalPoly(Ht, xt_slave(i, :), sigma) + syncv .* (x(i, :) - xt_slave(i, :)));
%     k2 = h * (EvalPoly(Ht, xt_slave(i, :) + k1, sigma) + syncv .* (x(i, :) - (xt_slave(i, :) + k1)));
%     xt_slave(i + 1, :) = xt_slave(i, :) + k2;
% 
%     k1 = h / 2 * (EvalPoly(Ho, xo_slave(i, :), sigma) + syncv .* (x(i, :) - xo_slave(i, :)));
%     k2 = h * (EvalPoly(Ho, xo_slave(i, :) + k1, sigma) + syncv .* (x(i, :) - (xo_slave(i, :) + k1)));
%     xo_slave(i + 1, :) = xo_slave(i, :) + k2;
% end

xo_slave = RK4SyncFOH(x, Ho_c, T, t, sync_coef);
xt_slave = RK4SyncFOH(x, Ht_c, T, t, sync_coef);

figure(1)
semilogy(t, vecnorm(x - xt_slave, 2, 2), 'Color', [1 0 0 0.2], 'LineWidth', 0.75);
hold on;
semilogy(t, vecnorm(x - xo_slave, 2, 2), 'Color', [0 0 1 0.2], 'LineWidth', 0.75);
semilogy(t, (t * 0) + mean(vecnorm(x - xt_slave, 2, 2)), 'Color', [1 0 0 1], 'LineWidth', 1);
semilogy(t, (t * 0) + mean(vecnorm(x - xo_slave, 2, 2)), 'Color', [0 0 1 1], 'LineWidth', 1);
grid on;
%title(['Synchronization error (', sysname, ')']);
legend('', '', 'LSM (mean)', 'Orthogonal polynomials (mean)');
xtickformat('$%g$'); ytickformat('$%g$'); ztickformat('$%g$');
set(gca, 'TickLabelInterpreter', 'latex');
%xlabel('Time $t$, s', 'Interpreter', 'latex');
ylabel('Synchronization error $\overline{\zeta}$', 'Interpreter', 'latex');
xlim([t(end - 25 / h), t(end)]);