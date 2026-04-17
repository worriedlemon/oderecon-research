rng_i default;
warning off;
try
    clf(1);
catch Me
    disp('Figures are already closed')
end

% Used system
sys = @Lorenz

sysname = func2str(sys);

% Used system coefficients
[Href, deg, vc, mc] = load_href(sysname);

Tmax = 30; % Time of experimental series
h = 1e-2; % Fixed step
start_point = load_pref(sysname); % Initial point

% go experiment
[t, x] = ode78(sys, 0:h:Tmax, start_point);

y = transpose(sys(0, x'));
N = size(x, 1);
eqc = vc;   

sigma = deglexord(deg, vc);
T = mat2cell(repmat(sigma, 1, eqc), mc, repmat(vc, 1, eqc));

lambda = min(abs(Href(Href ~= 0))) * 0.8; % SINDy sparsification parameter
nk = 141;
snrs = linspace(90, 20, nk);

noises = snr2amp(x, snrs);

errs = zeros(2, nk);
errw = zeros(2, nk);
erro = zeros(2, nk);
erros = zeros(2, nk);
errois = zeros(2, nk);

other = 1;
Fs = cell(1, vc);
local_attempts = 3;
for k = 1:nk
    disp(['Test #', num2str(k)]);
    for test_ = 1:local_attempts
        rx = x + noises(k) * randn(N, vc);
        ry = diff4(rx, t);
        if (other)
            %% SINDy
            Hs = sindy(rx, ry, sigma, lambda);
            errs(1, k) = norm(Href - Hs);
            Hs = mat2cell(Hs, mc, ones(1, eqc));
            [~, x_new] = ode45(@(t, x)oderecon(Hs, T, t, x), t, start_point);
            errs(2, k) = nn_norm(t, x_new, x);
        
            %% WSINDy
            Hw = wsindy(t, rx, sigma, lambda);
            errw(1, k) = norm(Href - Hw);
            Hw = mat2cell(Hw, mc, ones(1, eqc));
            [~, x_new] = ode45(@(t, x)oderecon(Hw, T, t, x), t, start_point);
            errw(2, k) = nn_norm(t, x_new, x);

            %% OrthPoly
            Ho = zeros(mc, eqc);
            F = orthpoly_t(sigma, t, x);
            E = EvalPoly(F', rx, sigma);
            for i = 1:eqc
                for j = 1:mc
                    Ho(j, i) = h_ji_byint(E(:, j), t, rx(:, i));
                end
            end
            Ho = F' * Ho;
            erro(1, k) = norm(Href - Ho);
        
            Ho = mat2cell(Ho, mc, ones(1, eqc));
            [~, x_new] = ode45(@(t, x)oderecon(Ho, T, t, x), t, start_point);
            erro(2, k) = nn_norm(t, x_new, x);
            
            %% OrthPoly-SINDy
            Hos = orthpoly_sindy(t,rx,sigma,lambda);
            erros(1, k) = norm(Href - Hos);
            Hos = mat2cell(Hos, mc, ones(1, eqc));
            [~, x_new] = ode45(@(t, x)oderecon(Hos, T, t, x), t, start_point);
            erros(2, k) = nn_norm(t, x_new, x);
        end
    
        %% OrthPoly-Int-SINDy
        Hois = orthpoly_int_sindy(t,rx,sigma,lambda);
        errois(1, k) = errois(1, k) + norm(Href - Hois);
        Hois = mat2cell(Hois, mc, ones(1, eqc));
        [~, x_new] = ode45(@(t, x)oderecon(Hois, T, t, x), t, start_point);
        errois(2, k) = errois(2, k) + nn_norm(t, x_new, x);
    end

    errs(:, k) = errs(:, k) / local_attempts;
    erro(:, k) = erro(:, k) / local_attempts;
    erros(:, k) = erros(:, k) / local_attempts;
    errois(:, k) = errois(:, k) / local_attempts;
end

figure(1)
subplot(2, 1, 1);
semilogy(snrs, errois(1, :), DisplayName='OrthPoly-Int-Sindy');
hold on; grid on;
if (other)
    semilogy(snrs, errs(1, :), DisplayName='SINDy');
    semilogy(snrs, errw(1, :), DisplayName='WSINDy');
    semilogy(snrs, erro(1, :), DisplayName='OrthPoly');
    semilogy(snrs, erros(1, :), DisplayName='OrthPoly-Sindy');
    legend show;
end
xtickformat('$%g$'); ytickformat('$%g$'); ztickformat('$%g$');
xlabel('SNR $\sigma$, db', Interpreter='latex');
ylabel('Coefficients error $\xi$', Interpreter='latex');
set(gca, TickLabelInterpreter='latex');
xlim([min(snrs), max(snrs)])
set(gca, Xdir='reverse');

subplot(2, 1, 2);
semilogy(snrs, errois(2, :), DisplayName='OrthPoly-Int-Sindy');
hold on; grid on;
if (other)
    semilogy(snrs, errs(2, :), DisplayName='SINDy');
    semilogy(snrs, errw(2, :), DisplayName='WSINDy');
    semilogy(snrs, erro(2, :), DisplayName='OrthPoly');
    semilogy(snrs, erros(2, :), DisplayName='OrthPoly-Sindy');
    legend show;
end
xtickformat('$%g$'); ytickformat('$%g$'); ztickformat('$%g$');
xlabel('SNR $\sigma$, db', Interpreter='latex');
ylabel('NNRMS error $\mathcal{{R}}$', Interpreter='latex');
set(gca, TickLabelInterpreter='latex');
xlim([min(snrs), max(snrs)])
set(gca, Xdir='reverse');