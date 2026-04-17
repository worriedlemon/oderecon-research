%% PREPARATION

rng_i default;
warning off;
%close all
figure(1);
clf

% Used system
sys = @Rossler
sysname = func2str(sys);

% Used system coefficients and initial point
[Href, deg, vc] = load_href(sysname);
start_point = load_pref(sysname);

Tmax = 25; % Time of experimental series
h = 1e-2; % Step
[t, x_ref] = ode78(sys, 0:h:Tmax, start_point);

N = size(x_ref, 1);
eqc = vc;

sigma = deglexord(deg, vc);
mc = size(sigma, 1);

T = mat2cell(repmat(sigma, 1, eqc), mc, repmat(vc, 1, eqc));

%% NOISING
% signal to noise ratio in decibels
snrdb = 30;

% get noise amplification factor
noise_amp = snr2amp(x_ref, snrdb);
x_noised = x_ref + noise_amp * randn(size(x_ref)); 

%% CALCULATION
lambda = 1e-1;
Hoi = orthpoly_int_sindy(t, x_noised, sigma, lambda)

H = mat2cell(Hoi, mc, ones(1, eqc));

prettyABM(H, T, 1);

% if (~input('To continue write 1, to stop write 0: '))
%    return
% end

[~, x] = ode78(@(t,x)oderecon(H, T, t, x), t, start_point);

%% PLOTS
enable_title = 0;
noised_opacity = 0.55;
fig = figure(1);
set(fig, 'Color', 'white');
used_letters = 'xy';

if vc == 2 % 2-dimensional graphics
    plot(x_ref(:, 1), x_ref(:, 2), Color=[0 0 1], DisplayName='Original System');
    hold on; grid on;
    plot(x_noised(:, 1), x_noised(:, 2), Color=[0 1 0.5 noised_opacity], DisplayName='Noised Data');
    plot(x(:, 1), x(:, 2), 'r', DisplayName='Reconstructed System');
else % >3-dimensional graphics
    axis_to_plot = [1, 2, 3];
    axis_letters = 'xyzuv';
    used_letters = axis_letters(axis_to_plot);
    
    plot3(x_ref(:, axis_to_plot(1)), x_ref(:, axis_to_plot(2)), x_ref(:, axis_to_plot(3)), Color=[0 0 1], DisplayName='Original System');
    hold on; grid on;
    plot3(x_noised(:, axis_to_plot(1)), x_noised(:, axis_to_plot(2)), x_noised(:, axis_to_plot(3)), Color=[0 1 0.5 noised_opacity], DisplayName='Noised Data');
    plot3(x(:, axis_to_plot(1)), x(:, axis_to_plot(2)), x(:, axis_to_plot(3)), 'r', DisplayName='Reconstructed System');
    ztickformat('$%g$');
    zlabel(['$', used_letters(3), '$'], Interpreter='latex');
end

legend show;
enable_title && title(['Phase portrait for the ', func2str(sys), ' system']) == ""; % hack for disabling title
set(gca, 'TickLabelInterpreter', 'latex');
xtickformat('$%g$'); ytickformat('$%g$');
xlabel(['$', used_letters(1), '$'], Interpreter='latex');
ylabel(['$', used_letters(2), '$'], Interpreter='latex');