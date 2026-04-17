rng_i shuffle % maybe default is better
warning off
close all

% Used system
sys = @Lorenz

sysname = func2str(sys);

[Href, deg, vc, mc] = load_href(sysname); % Used coefficients

T_range = [10; 30];    % relatively large T
h_range = [-3; -1];  % relatively small h
x_range = [-5; 5];
snr_range = [30; 70];

sigma = deglexord(deg, vc);

addpath('research2/sindy_interface/')

% Methods:
% 1. SINDY
% 2. WSINDy
% 3. Orthpoly-SINDY
% 4. Ohrtpoly-Int-SINDY
methods = get_methods;
mtdn = size(methods, 2);

% Lambda parameters
default_lambda = min(abs(Href(Href ~= 0))) * 0.8; % taking lambda as 0.8 of the smallest existing coefficient

tests = 300;

slices = 10000;
h = logspace(h_range(1), h_range(2), slices);
snr = linspace(snr_range(1), snr_range(2), slices);
Tmax = linspace(T_range(1), T_range(2), slices);

data = zeros(tests, 7);
for p = 1:tests
    ind = randi(slices, 3, 1);
    h_cur = h(ind(1));
    snr_cur = snr(ind(2));
    Tmax_cur = Tmax(ind(3));

    data(p, 1) = h_cur;
    data(p, 2) = snr_cur;
    data(p, 3) = Tmax_cur;

    start_point = affine_transform(rand(1, vc), x_range, [0; 1]);
    fprintf("Test #%d: Tmax = %f, h = %f, SNR = %fdb\n", p, Tmax_cur, h_cur, snr_cur);

    [t, x_ref] = RK7(sys, 0:h_cur:Tmax_cur, start_point');
    snr_amp = snr2amp(x_ref, snr_cur);
    x_noised = x_ref + snr_amp * randn(size(x_ref));
    
    data_temp = zeros(mtdn, 2);
    for k = 1:mtdn
        [xi_temp, r_temp] = get_norms(Href, t, x_noised, sigma, default_lambda, methods{1, k});
        data_temp(k, :) = [xi_temp, r_temp];
    end

    [data(p, 6:7), data(p, 4:5)] = min(data_temp, [], 1);
end

indexes = nchoosek(1:3, 2);
log_plot = logical([1, 0, 0]);

method_colors = hsv(mtdn);
opacity = 0.5;

labels = { 'Step size $h$', 'SNR, db', 'Simulation time $T_{max}$, s' };

for fig = 1:2
    figure(fig);
    if (fig == 1)
        sgtitle('Coefficients norm $\Xi$', Interpreter='latex');
    else
        sgtitle('NNRMS norm $\mathcal{R}$', Interpreter='latex');
    end
    for d = 1:mtdn
        draw_data = data(data(:, fig + 3) == d, :);

        for i = 1:3
            subplot(1, 3, i);

            F = plot(draw_data(:, indexes(i, 1)), draw_data(:, indexes(i, 2)), MarkerEdgeColor='none', MarkerFaceColor=method_colors(d, :), LineStyle='none', Marker='o');
            if (isempty(F))
                continue
            end

            if log_plot(indexes(i, 1))
                F.Parent.XScale = 'log';
                F.Parent.XScale = 'log';
            end
            if log_plot(indexes(i, 2))
                F.Parent.YScale = 'log';
            end

            xlabel(labels{1, indexes(i, 1)}, Interpreter='latex');
            ylabel(labels{1, indexes(i, 2)}, Interpreter='latex');
            hold on;
        end
    end
    
    leg = cell(1, mtdn);
    for d = 1:mtdn
        leg{1, d} = method2str(methods{1, d});
    end
    legend(leg(1, :), Interpreter='none');
end

%rmpath('research2/sindy_interface/')