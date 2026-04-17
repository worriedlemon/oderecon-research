rng_i default;
close all;
warning off;

addpath("research\mycalc")

% Dynamic system Simulation
start_point = [4 -2 0]; % Initial point
h = 0.01; % Step
Tspan = [10, 100:500:10000] * h * 10;

deg = 2; % Degree of reconstructed function
vc = 3; % Variables count
eqc = 3; % Equations count
mc = nchoosek(deg + vc, vc); % Monomials count

sigma = deglexord(deg, vc);

Ns = zeros(length(Tspan), 1);
Es = [Ns, Ns];
Orths = Ns;

for i = 1:length(Tspan) % Time max
    disp(i)
    Tmax = Tspan(i);
    [t, x] = ode45(@Rossler, 0:h:Tmax, start_point);
    y = transpose(Rossler(0, x'));
    Ns(i) = length(t);

    tic; % Timer start
    F = mytranspose(orthpoly_t(sigma, t, x));
    Orths(i) = toc;

    tic; % Timer start
    E = EvalPoly(F, x, sigma);
    
    coefs = zeros(mc, eqc);
    for eq = 1:eqc
        coefs(:, eq) = trapz(x(:, eq), E);
    end

    coefs = mymult(F, coefs);
    Es(i, 1) = toc;
    
    tic;
    E = EvalPoly(eye(mc), x, sigma);
    coefs = mylsm(E, y);
    Es(i, 2) = toc; % timer end
end

semilogy(Ns, Es(:, 1), 'b', Ns, Orths + Es(:, 1), 'g', Ns, Es(:, 2), 'r');
grid on;
xlabel('Points count \it{N}'); ylabel('Elapsed time \it{t}, s');
title('Speed comparison between LSM and Orth');
legend('Orthogonal polynomials (algorithm only)', 'Orthogonal polynomials', 'LSM');

rmpath("research\mycalc")