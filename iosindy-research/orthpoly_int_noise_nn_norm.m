rng default;
warning off;
close all

% Used system
system = @Rossler

% Used system coefficients
[Href, deg, vc] = load_href(func2str(system));

Ttrans = 100; % Transient time
%Tmax = 5; % Time of experimental series for Lorenz
Tmax = 15; % Time of experimental series for Rossler
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
T = mat2cell(repmat(sigma, 1, eqc), mc, repmat(vc, 1, eqc));

delta = 1e-7;     % regularization parameter
%tol = 5e-2;       % DMT tolerance
tol = 4e-3;
eta =4e-3; %LSM + DMT
lambda = 1e-2;    % SINDy sparsification parameter
lambda_orth = 1e-2;  % Orthpoly + SINDy sparsification parameter

noises = logspace(-5, 1, 50);

% Homoscedasticity test
errt = []; errs = []; erro = []; erros = []; errois = [];
Fs = cell(1, vc);
for noise_amp = noises
    rx = x + noise_amp * randn(N, vc);
    ry = diff4(rx, t);
    
    F = orthpoly_t_integration(sigma,t,rx);
    
    %% LSM
    B = EvalPoly(eye(mc), rx, sigma);
    Ht = (B'*B)\B'*ry;
    Ht = mat2cell(Ht, mc, ones(1, eqc));
    [~, x_new] = ode78(@(t, x)oderecon(Ht, T, t, x), t, start_point);
    errt = [errt nn_norm(t, x_new, x)];
    
    %% SINDy
    Hs = sindy(rx, ry, sigma, lambda);
    Hs = mat2cell(Hs, mc, ones(1, eqc));
    [~, x_new] = ode78(@(t, x)oderecon(Hs, T, t, x), t, start_point);
    errs = [errs nn_norm(t, x_new, x)];

    %% OrthPoly
    Ho = zeros(mc, eqc);
    E = EvalPoly(F', rx, sigma);
    for i = 1:eqc
        for j = 1:mc
            Ho(j, i) = h_ji_byint(E(:, j), t, rx(:, i));
        end
    end
    Ho = F' * Ho;

    Ho = mat2cell(Ho, mc, ones(1, eqc));
    [~, x_new] = ode78(@(t, x)oderecon(Ho, T, t, x), t, start_point);
    erro = [erro nn_norm(t, x_new, x)];
    
    %% OrthPoly-SINDy
    % Hos = orthpoly_sindy(t,rx,sigma,lambda_orth);
    % Hos = mat2cell(Hos, mc, ones(1, eqc));
    % [~, x_new] = ode78(@(t, x)oderecon(Hos, T, t, x), t, start_point);
    % erros = [erros nn_norm(t, x_new, x)];

    %% OrthPoly-Int-SINDy
    Hois = orthpoly_int_sindy(t,rx,sigma,lambda_orth);
    Hois = mat2cell(Hois, mc, ones(1, eqc));
    [~, x_new] = ode78(@(t, x)oderecon(Hois, T, t, x), t, start_point);
    errois = [errois nn_norm(t, x_new, x)];
end

figure(1)
loglog(noises, errt, 'r', DisplayName='LSM (Homoscedastic)');
hold on; grid on;
loglog(noises, errs, 'b', DisplayName='SINDy (Homoscedastic)');
loglog(noises, erro, 'g', DisplayName='OrthPoly (Homoscedastic)');
loglog(noises, erros, 'c--', DisplayName='OrthPoly-Sindy (Homoscedastic)');
loglog(noises, errois, '--', 'Color',[0 0.7 0.7],DisplayName='OrthPoly-Int-Sindy (Heteroscedastic)');
%title(['Noise resistance (', func2str(system), ')']);
legend show;
xtickformat('$%g$'); ytickformat('$%g$'); ztickformat('$%g$');
xlabel('$\sigma$', 'Interpreter', 'latex');
ylabel('Nearest Neighbor norm $\Xi$', 'Interpreter', 'latex');
set(gca, 'TickLabelInterpreter', 'latex');
xlim([noises(1), noises(end)])

