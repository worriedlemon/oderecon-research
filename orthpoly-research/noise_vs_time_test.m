rng_i default;
warning off;
close all

% System coefficients
Hlorenz = [0 -10 10 0 0 0 0 0 0 0; 0 28 -1 0 0 0 -1 0 0 0; 0 0 0 -8/3 0 1 0 0 0 0]';
Hrossler = [0 0 -1 -1 0 0 0 0 0 0; 0 1 0.3 0 0 0 0 0 0 0; 0.3 0 0 -5.7 0 0 1 0 0 0]';

% Used system
%system = @Rossler
system = @Lorenz

Href = eval(['H', lower(func2str(system))]); % Used coefficients
Ttrans = 100; % Transient time

h = 1e-2; % Step
start_point = [4 -2 0]; % Initial point
noise_amp = 1e-2;

% %first, go transient
[~, x] = ode78(system, 0:h:Ttrans, start_point);
start_point = x(end,:);

eqc = size(y, 2);
deg = 2;

sigma = deglexord(deg, vc);
mc = size(sigma, 1);

delta = 1e-7; % regularization parameter

%Tmax = 5; % Time of experimental series for Lorenz
%Tmax = 10; % Time of experimental series for Rossler

times = logspace(0, 2, 200);

errt = []; erro = [];
hw = waitbar(0,'Please wait...','Name','Calculate error dependence on data...');
for Tmax = times
    %then, go experiment
    [t, x] = ode78(system, 0:h:Tmax, start_point);
    y = transpose(system(0, x'));

    waitbar(Tmax/times(end),hw,['Processing, T = ',num2str(Tmax)]);

    % figure(99);
    % plot(x(:,1),x(:,2));
    % xlabel('x'); ylabel('y');

    [N, vc] = size(x);

    rx = x + noise_amp * randn(N, vc);
    %ry = [diff(rx) / h; (rx(end, :) - rx(end - 1, :)) / h]; % first order
    %ry = diff2(rx) / h; % second order
    ry = diff4(rx, t); % fourth order
    
    F = orthpoly_t4(sigma, t, rx);
    
    Ho = zeros(mc, eqc);
    E = EvalPoly(F', rx, sigma);
    for i = 1:eqc
        Ho(:, i) = trapz(rx(:, i), E);
    end

    B = EvalPoly(eye(mc), rx, sigma);
    Ht = (B'*B + delta*eye(mc))\B'*ry;
    errt = [errt norm(Ht - Href)];

    Ht1 = F' * Ho;
    erro = [erro norm(Ht1 - Href)];
end
close(hw);

figure(1);
loglog(times, errt, 'r', times, erro, 'b');

xlim([times(1), times(end)])
hold on; grid on;
legend('LSM', 'OrthPoly');
xlabel('$T$', 'Interpreter', 'latex');
ylabel('Coefficients error $\zeta$', 'Interpreter', 'latex');
set(gca, 'TickLabelInterpreter', 'latex');

set(gcf,'position',[450  230  372  240]);
xtickformat('$%g$'); ytickformat('$%g$');
set(gca, 'TickLabelInterpreter', 'latex');