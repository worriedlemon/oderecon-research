rng_i default;
close all;
warning off;

sys = @Rossler %#ok

sysname = func2str(sys);
[Href, deg, vc] = load_href(sysname); % Used coefficients

start_point = [4 -2 0]; % Initial point
Tmax = 100; % Time end
h = 1e-2; % Step

[t, x] = ode45(sys, 0:h:Tmax, start_point);

eqc = vc; % Equations count
sigma = deglexord(deg, vc);

F = orthpoly_t(sigma, t, x) %#ok Getting orthogonal polynomials matrix
mc = size(F, 1); % Monomials count

coefs = zeros(mc, eqc);
coefs_reg = zeros(mc, eqc);

tols = logspace(-7, -1, 56);
mcs = zeros(1, length(tols));
errs = mcs;
y = diff4(x,t);
for i = 1:length(tols)
    tol = tols(i);
    for j = 1:vc
        [coefs(:, j), ~, ~, coefs_reg(:, j)] = delMinorTerms_dy(t, x(:, j), x, y(:, j), F, sigma, tol, 0);
        mcs(i) = sum(coefs ~= 0, "all");
    end

    errs(i) = norm(coefs_reg - Href);
end

act = mcs * 0 + sum(Href ~= 0, "all");

figure(1);

subplot(2, 1, 1);
semilogx(tols, mcs, "b-", Marker='.'); hold on;
semilogx(tols, act, "m");
legend('', 'Actual monomials count', Interpreter='latex');
grid on;
xtickformat('$%g$'); ytickformat('$%g$');
set(gca, TickLabelInterpreter='latex', xdir='reverse');
xlabel('Tolerance $\eta$','Interpreter','latex');
ylabel('Monomials count $L$','Interpreter','latex');
title(['Monomials count depending on tolerance (', sysname, ')']);
xlim([tols(1), tols(end)]); ylim([0, mc * vc]);

subplot(2, 1, 2);
loglog(tols, errs, "r-", Marker='.');
grid on;
xtickformat('$%g$'); ytickformat('$%g$');
set(gca, TickLabelInterpreter='latex', xdir='reverse');
xlabel('Tolerance $\eta$','Interpreter','latex');
ylabel('Coefficients error $\zeta$','Interpreter','latex');
xlim([tols(1), tols(end)]);
