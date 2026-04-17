rng_i default;
warning off;

%system = @Rossler
system = @Lorenz

[Href, deg, vc] = load_href(func2str(system)); % Used coefficients

Tmax = 10;
Tmaxs = 10:10:100; % Time end
hs = 10.^(-3:0.01:-1);
h = 1e-3; % Step
start_point = [4 -2 0]; % Initial point

delta = 1e-5; % regularization parameter

eqc = vc;

sigma = deglexord(deg, vc);
mc = size(sigma, 1);

nrm = zeros(3, length(hs));
for i = 1:length(hs)
    disp(i)
    [t, x] = ode78(system, 0:hs(i):Tmax, start_point);
    y = [diff(x); (x(end, :) - x(end - 1, :))] / h; % first order
    %y = diff4(x, t); % fourth order
    
    B = EvalPoly(eye(mc), x, sigma);
    Ht = (B'*B + delta*eye(mc))\B'*y;

    nrm(1, i) = norm(Ht - Href);
    
    F = orthpoly_t(sigma, t, x);

    Ho_t = zeros(mc, eqc);
    Ho_x = Ho_t;
    for k = 1:eqc
        for j = 1:mc
            %Ho_x(j, k) = trapz(t, EvalPoly(F(j, :)', x, sigma) .* y(:, k));
            Ho_x(j, k) = trapz(x(:, k), EvalPoly(F(j, :)', x, sigma));
            %Ho_x(j, k) = integrate_simpvar(x(:, k), EvalPoly(F(j, :)', x, sigma));
            %Ho_x(j, k) = integrate_simp(t, EvalPoly(F(j, :)', x, sigma) .* y(:, k));
        end
    end

    Ho_t = F' * Ho_t;
    Ho_x = F' * Ho_x;

    %nrm(2, i) = norm(Ho_t - Href);
    nrm(3, i) = norm(Ho_x - Href);
end

figure(1);
%plot(Tmaxs, nrm(1, :), 'b--', Tmaxs, nrm(2, :), 'r--', Tmaxs, nrm(3, :), 'g');
loglog(hs, nrm(1, :), 'r', hs, nrm(3, :), 'b');
grid on;
%title(['Reconstruction error (', sysname, ')']);
xtickformat('$%g$'); ytickformat('$%g$'); ztickformat('$%g$')
set(gca,'TickLabelInterpreter','latex');
%xlabel('Time end $T_\text{max}$, s');
%xlabel('Time step $h$, s','Interpreter','latex');
ylabel('Error $\zeta$','Interpreter','latex');
legend('LSM', 'Orthogonal polynomials')

hold on


% figure(2);
% %plot(Tmaxs, nrm(1, :), 'b--', Tmaxs, nrm(2, :), 'r--');
% semilogx(hs, nrm(1, :), 'b--', hs, nrm(2, :), 'r--');
% grid on;
% title(['Reconstruction accuracy (', sysname, ')']);
% %xlabel('Time end \it{T}_{max}, s');
% xlabel('Time step \it\Delta{t}, s');
% ylabel('Accuracy \it{\xi}');
% legend('LSM', 'Orthogonal polynomials (time)');