rng_i default;
warning off;
close all

% Used system
%system = @Rossler
system = @Lorenz

[Href, deg, vc] = load_href(func2str(system)); % Used coefficients
Tmax = 10; % Time end
h = 1e-3; % Step
start_point = [4 -2 0]; % Initial point

%first, go transient
Ttrans = 100;
[~, x] = ode78(system, 0:h:Ttrans, start_point);
start_point = x(end,:);

%then, get full dataset
[t, x] = ode78(system, 0:h:Tmax, start_point);
y = transpose(system(0, x'));

N = size(x, 1);
eqc = vc;

sigma = deglexord(deg, vc);
mc = size(sigma, 1);

delta = 1e-3; % regularization parameter

noises = logspace(-4, 1, 50);
noise_M = pinknoise(N, vc);
noise_M_mean = mean(abs(noise_M))
noise_M = noise_M ./ noise_M_mean;

errt = []; erro = [];
for noise_amp = noises
    rx = x + noise_amp * noise_M;
    %ry = [diff(rx) / h; (rx(end, :) - rx(end - 1, :)) / h]; % first order
    %ry = diff2(rx) / h; % second order
    ry = diff4(rx, t); % fourth order
    
    F = orthpoly_t(sigma, t, rx);
    
    Ho = zeros(mc, eqc);
    E = EvalPoly(F', rx, sigma);
    for i = 1:eqc
        for j = 1:mc
            %Ho(j, i) = trapz(t, E(:, j) .* ry(:, i));
            %Ho(j, i) = trapz(rx(:, i), E(:, j));
            
            %val = integrate_simp(E(:, j) .* ry(:, i),0,h);
            %Ho(j, i) = val(end);
            %Ho(j, i) = integrate_simpvar(t, E(:, j) .* ry(:, i));

            %Ho(j, i) = integrate_simpvar(rx(:, i),E(:, j));
            Ho(j, i) = intdiff4(rx(:, i),E(:, j));
            %Ho(j, i) = simps(rx(:, i),E(:, j));
        end
    end

    B = EvalPoly(eye(mc), rx, sigma);
    Ht = (B'*B + delta*eye(mc))\B'*ry;
    errt = [errt norm(Ht - Href)];

    Ht1 = F' * Ho;
    erro = [erro norm(Ht1 - Href)];
end

figure(1);
loglog(noises, errt, 'r', noises, erro, 'b');
hold on; grid on;
%title(['Noise resistance (', func2str(system), ')']);
legend('LSM', 'Orthogonal polynomials');
xtickformat('$%g$'); ytickformat('$%g$'); ztickformat('$%g$');
xlabel('Average pink noise scaling factor $K$', 'Interpreter', 'latex');
ylabel('Coefficients error $\zeta$', 'Interpreter', 'latex');
set(gca, 'TickLabelInterpreter', 'latex');
xlim([noises(1), noises(end)]);