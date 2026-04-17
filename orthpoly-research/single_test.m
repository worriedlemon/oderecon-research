warning off;

sys = @Rossler;

[Href, deg, vc] = load_href(func2str(sys));

Tmax = 50;

hs = logspace(-4, -1, 20);
N = length(hs);

eqc = vc;
sigma = deglexord(deg, vc);
mc = size(sigma, 1);

delta = 0.01; % regularization parameter

errt = zeros(1, N);
erro = errt;
for k = 1:N
    disp(k);
    
    t = single(0:hs(k):Tmax);
    [~, x] = ode45(sys, t, [4 -2 0]);
    y = diff4(x, t);
    %y = [diff(x) / h; (x(end, :) - x(end - 1, :)) / h];

    F = single(orthpoly_t(sigma, t, x));

    E = single(EvalPoly(F', x, sigma));
    Ho = zeros(mc, eqc);
    for i = 1:eqc
        Ho(:, i) = trapz(x(:, i), E);
    end

    Ho = F' * Ho;

    E = EvalPoly(eye(mc), x, sigma);
    Ht = (E'*E + delta*eye(mc))\E'*y;

    errt(k) = norm(Ht - Href);
    erro(k) = norm(Ho - Href);
end

figure(1)
loglog(hs, errt + 1e-14, 'r', hs, erro + 1e-14, 'b');
grid on;
title(['Reconstruction with reduced precision (', sysname, ')']);
legend('LSM', 'Orthogonal Polynomials');
xlabel('Simulation time step \it{h}'); ylabel('Coefficients error \zeta');