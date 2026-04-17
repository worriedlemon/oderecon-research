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
start_point = x(end,:);
%then, go experiment
[t, x] = ode78(system, 0:h:Tmax, start_point);

y = transpose(system(0, x'));

N = size(x, 1);
eqc = vc;

sigma = deglexord(deg, vc);
mc = size(sigma, 1);

delta = 1e-7;     % regularization parameter
tol = 5e-2;       % DMT tolerance
lambda = 1e-1;    % SINDy sparsification parameter
lambda_orth = 1;  % Orthpoly + SINDy sparsification parameter

noises = logspace(-5, 1, 200);

% Homoscedasticity test
errt = []; erro = []; errs = []; erro2 = []; erros = [];
Fs = cell(1, vc);
for noise_amp = noises
    rx = x + noise_amp * randn(N, vc);
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
            %Ho(j, i) = E(end, j)*rx(end, i) - E(1, j)*rx(1, i) -  integrate_simpvar(E(:, j), rx(:, i));


            %Ho(j, i) = integrate_simpvar(t, E(:, j) .* ry(:, i));
            Ho(j, i) = intdiff4(rx(:, i),E(:, j));
           
            
            %Ho(j, i) = integrate_simpvar(rx(:, i), E(:, j));
            % tmp = integrate_bool(E(:, j) .* ry(:, i),0,h);
            % Ho(j, i) = tmp(end);
        end
    end

    B = EvalPoly(eye(mc), rx, sigma);
    Ht = (B'*B + delta*eye(mc))\B'*ry;
    errt = [errt norm(Ht - Href)];

    Ht1 = F' * Ho;
    erro = [erro norm(Ht1 - Href)];

    Ht1dmt = Ht1;
    tol = 5e-2;
    for j = 1:vc
        [~, ~, ~, Ht1dmt(:, j)] = delMinorTerms_dy(t, x(:, j), x, y(:, j), F, sigma, tol, 0);
    end
    erro2 = [erro2 norm(Ht1dmt - Href)];
    
    Ht2 = sindy(rx, ry, sigma, lambda);
    errs = [errs norm(Ht2 - Href)];

    Ht3 = Ho;
    while (1)
        smallinds = (abs(Ht3) < lambda_orth);
        if (all(Ht3(smallinds) == 0, "all")) % all small values are already zeros, stop
            break
        end
    
        Ht3(smallinds) = 0;
        for ind = 1:vc
            biginds = ~smallinds(:, ind);
            sigma_temp = sigma(biginds, :);
            F_temp = zeros(size(F));
            F_temp(biginds, biginds) = orthpoly_t(sigma_temp, t, rx);
            E = EvalPoly(F_temp', rx, sigma);
            for j = 1:mc
                Ht3(j, ind) = trapz(rx(:, ind), E(:, j));
            end
            Fs{1, ind} = F_temp;
        end
    end

    for i = 1:vc
        Ht3(:, i) = Fs{1, i}' * Ht3(:, i);
    end
    erros = [erros norm(Ht3 - Href)];
end

figure(1);
loglog(noises, errt, 'r', DisplayName='LSM (Homoscedastic)');
hold on; grid on;
loglog(noises, errs, 'b', DisplayName='SINDy (Homoscedastic)');
loglog(noises, erro, 'g', DisplayName='OrthPoly (Homoscedastic)');
loglog(noises, erro2, 'm', DisplayName='OrthPoly DMT (Homoscedastic)');
loglog(noises, erros, 'c', DisplayName='OrthPoly-Sindy (Homoscedastic)');
%title(['Noise resistance (', func2str(system), ')']);
legend show;
xtickformat('$%g$'); ytickformat('$%g$'); ztickformat('$%g$');
xlabel('$\sigma$', 'Interpreter', 'latex');
ylabel('Coefficients error $\zeta$', 'Interpreter', 'latex');
set(gca, 'TickLabelInterpreter', 'latex');
xlim([noises(1), noises(end)])

% uncomment to skip heteroscedasticity calculation
% return

% Heteroscedasticity test
errt = []; erro = []; errs = []; erro2 = []; erros = [];
for noise_amp = noises
    dmx = x - mean(x);
    rx = x + noise_amp * randn(N, vc) .* dmx;
    %ry = [diff(rx) / h; (rx(end, :) - rx(end - 1, :)) / h];
    %ry = diff2(rx) / h; % second order
    ry = diff4(rx, t); % fourth order
    F = orthpoly_t(sigma, t, rx);
    
    Ho = zeros(mc, eqc);
    E = EvalPoly(F', rx, sigma);
    for i = 1:eqc 
        for j = 1:mc
            %Ho(j, i) = trapz(rx(:, i), E(:, j));
            %Ho(j, i) = integrate_simpvar(t, E(:, j) .* ry(:, i));
            %Ho(j, i) = intdiff2(rx(:, i),E(:, j));
            %Ho(j, i) = trapz(t, E(:, j) .* ry(:, i));
            Ho(j, i) = intdiff4(rx(:, i), E(:, j));
        end
    end

    B = EvalPoly(eye(mc), rx, sigma);
    Ht = (B'*B)\B'*ry;
    errt = [errt norm(Ht - Href)];

    Ht1 = F' * Ho;
    erro = [erro norm(Ht1 - Href)];

    Ht1dmt = Ht1;
    tol = 5e-2;
    for j = 1:vc
        [~, ~, ~, Ht1dmt(:, j)] = delMinorTerms_dy(t, x(:, j), x, y(:, j), F, sigma, tol, 0);
    end
    erro2 = [erro2 norm(Ht1dmt - Href)];
    
    Ht2 = sindy(rx, ry, sigma, lambda);
    errs = [errs norm(Ht2 - Href)];

    Ht3 = Ho;
    while (1)
        smallinds = (abs(Ht3) < lambda_orth);
        if (all(Ht3(smallinds) == 0, "all")) % all small values are already zeros, stop
            break
        end
    
        Ht3(smallinds) = 0;
        for ind = 1:vc
            biginds = ~smallinds(:, ind);
            sigma_temp = sigma(biginds, :);
            F_temp = zeros(size(F));
            F_temp(biginds, biginds) = orthpoly_t(sigma_temp, t, rx);
            E = EvalPoly(F_temp', rx, sigma);
            Ht3(:, ind) = trapz(rx(:, ind), E);
            Fs{1, ind} = F_temp;
        end
    end

    for i = 1:vc
        Ht3(:, i) = Fs{1, i}' * Ht3(:, i);
    end
    erros = [erros norm(Ht3 - Href)];
end

figure(2);
loglog(noises, errt, 'r', DisplayName='LSM (Heteroscedastic)');
hold on; grid on;
loglog(noises, errs, 'b', DisplayName='SINDy (Heteroscedastic)');
loglog(noises, erro, 'g', DisplayName='OrthPoly (Heteroscedastic)');
loglog(noises, erro2, 'm', DisplayName='OrthPoly DMT (Heteroscedastic)');
loglog(noises, erros, 'c', DisplayName='OrthPoly-Sindy (Heteroscedastic)');
%title(['Noise resistance (', func2str(system), ')']);
legend show;
xtickformat('$%g$'); ytickformat('$%g$'); ztickformat('$%g$');
xlabel('$\sigma$', 'Interpreter', 'latex');
ylabel('Coefficients error $\zeta$', 'Interpreter', 'latex');
set(gca, 'TickLabelInterpreter', 'latex');
xlim([noises(1), noises(end)])