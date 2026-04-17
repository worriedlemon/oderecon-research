warning off;
%close all

%sys = @Lorenz
sys = @Rossler
[Href, deg, vc] = load_href(func2str(sys));

Tmax = 5; % Rossler
hs = logspace(-4,-1,50);
% Tmax = 30; % Lorenz
% hs = logspace(-4,-1,50);

eqc = vc;
sigma = deglexord(deg, vc);
mc = size(sigma, 1);

N = length(hs);

delta = 1e-7; % regularization parameter

start_point = [4 -2 0]; % Initial point

% first, go transient
Ttrans = 101;
[~, x] = ode78(sys, 0:0.001:Ttrans, start_point);
start_point = x(end,:);

errt = zeros(1, N);
errt1 = errt;
erro = errt;
for k = 1:N
    disp(k);
    h = hs(k);
    t = 0:h:Tmax;
    [~, x] = ode78(sys, t, start_point);
    y1 = [diff(x) / h; (x(end, :) - x(end - 1, :)) / h];
    y = diff2(x)/h;
    F = orthpoly_t(sigma, t, x);
    E = EvalPoly(F', x, sigma);

    Ho = zeros(mc, eqc);
    for i = 1:eqc
        Ho(:, i) = trapz(x(:, i), E);
    end

    Ho = F' * Ho;

    E = EvalPoly(eye(mc), x, sigma);
    Ht = (E'*E + delta*eye(mc))\E'*y;
    Ht1 = (E'*E + delta*eye(mc))\E'*y1;

    errt(k) = norm(Ht - Href);
    errt1(k) = norm(Ht1 - Href);
    erro(k) = norm(Ho - Href);
end

figure(1)
loglog(hs, errt1 + 1e-14, 'c', hs, errt + 1e-14, 'r', hs, erro + 1e-14, 'b');
grid on;
hold on
title(func2str(sys));
xlabel('Simulation time step \it{h}'); ylabel('Coefficients error \zeta');

errt = zeros(1, N);
erro = errt;
for k = 1:N
    disp(k);
    t = 0:hs(k):Tmax;
    [~, x] = ode78(sys, t, start_point);
    y = diff4(x,t);
    F = orthpoly_t(sigma, t, x);
    E = EvalPoly(F', x, sigma);

    Ho = zeros(mc, eqc);
    for i = 1:eqc
        for j = 1:mc
            Ho(j, i) = intdiff4(x(:, i),E(:, j));
        end
    end

    Ho = F' * Ho;

    E = EvalPoly(eye(mc), x, sigma);
    Ht = (E'*E + delta*eye(mc))\E'*y;

    errt(k) = norm(Ht - Href);
    erro(k) = norm(Ho - Href);
end

figure(1)
loglog(hs, errt + 1e-14, 'm-', hs, erro + 1e-14, 'g-');
grid on;
hold on
%title([sysname, ' - 4 order diff + int']);
legend('LSM 1 order','LSM 2 order', 'Orthogonal Polynomials 2 order','LSM 4 order', 'Orthogonal Polynomials 4 order');
xlabel('Simulation time step \it{h}'); ylabel('Coefficients error \zeta');

set(gcf,'position',[181  118  500  250]);
xtickformat('$%g$'); ytickformat('$%g$');
set(gca, 'TickLabelInterpreter', 'latex');
