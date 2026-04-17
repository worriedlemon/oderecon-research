rng_i shuffle % maybe default is better
warning off
close all

% Used system
sys = @Zarei

sysname = func2str(sys);

[Href, deg, vc, mc] = load_href(sysname); % Used coefficients

T_range = [5; 15];    % relatively large T
h_range = [1e-4; 1e-2];  % relatively small h
x_range = [-5; 5];
snr_range = [70; 90];
c01 = [0; 1];

sigma = deglexord(deg, vc);

addpath('research2/sindy_interface/')

% Methods:
% 1. SINDY
% 2. Orthpoly-SINDY
% 3. Ohrtpoly-Int-SINDY
methods = get_methods;
mtdn = size(methods, 2);

% Lambda parameters
default_lambda = min(abs(Href(Href ~= 0))) * 0.8; % taking lambda as 0.8 of the smallest existing coefficient

% Pairwise compared methods
methods_pairs = nchoosek(1:mtdn, 2);
methods_pairs_size = size(methods_pairs, 1);

tests = 50;
pair_resultsXi = -ones(mtdn);
pair_resultsR = pair_resultsXi;

data = cell(tests * methods_pairs_size, 2);
for p = 1:methods_pairs_size
    [first_method, second_method] = methods{1, methods_pairs(p, :)};
    disp(['Test method pair ', num2str(p), ': ', method2str(first_method), ' vs ', method2str(second_method)])
    
    winsXi = 0;
    winsR = 0;
    for k = 1:tests        
        % Take random
        Tmax = affine_transform(rand(1), T_range, c01);
        start_point = affine_transform(rand(1, vc), x_range, c01);
        h = affine_transform(rand(1), h_range, c01);
        snr = affine_transform(rand(1), snr_range, c01);

        [t, x_ref] = RK4(sys, 0:h:Tmax, start_point');
        snr_amp = snr2amp(x_ref, snr);
        x_noised = x_ref + snr_amp * randn(size(x_ref));

        fprintf("Test #%d: Tmax = %f, h = %f, SNR = %fdb\n", k, Tmax, h, snr);
    
        [Xi1, Rscript1] = get_norms(Href, t, x_noised, sigma, default_lambda, first_method);
        [Xi2, Rscript2] = get_norms(Href, t, x_noised, sigma, default_lambda, second_method);
    
        winsXi = winsXi + (Xi1 < Xi2);
        winsR = winsR + (Rscript1 < Rscript2);
    end

    pair_resultsXi(methods_pairs(p, 1), methods_pairs(p, 2)) = winsXi;
    pair_resultsXi(methods_pairs(p, 2), methods_pairs(p, 1)) = tests - winsXi;
    pair_resultsR(methods_pairs(p, 1), methods_pairs(p, 2)) = winsR;
    pair_resultsR(methods_pairs(p, 2), methods_pairs(p, 1)) = tests - winsR;
end

disp('Results matrixes:');
Xi = cat(2, pair_resultsXi, sum(pair_resultsXi, 2) + 1)
Rho = cat(2, pair_resultsR, sum(pair_resultsR, 2) + 1)

rmpath('research2/sindy_interface/')