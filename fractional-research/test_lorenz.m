rng_i default

sys = @Lorenz;
Tmax = 30;
h = 1e-3;
start_point = [4 -2 0];

alpha = 1.16;
%alpha = [1.085 0.982 1.024]; % exotic alphas
[t, x] = fractode(sys, 0:h:Tmax, start_point, alpha);

plot3(x(:, 1), x(:, 2), x(:, 3));
hold on; grid on;

% suppose, x, alpha and t are known
y = fractdiff(x, t, alpha, 'gl');

% Guessed max degree
deg = 2;
vc = length(start_point);
sigma = deglexord(deg, vc);
mc = size(sigma, 1);

% SINDy reconstruction
Hcoef = sindy(x, y, sigma, 0.5);

disp('Coefficients are:');
disp(Hcoef);

[H, T] = HTmat2cell(Hcoef, sigma);

[~, xr] = fractode(@(t,x)oderecon(H, T, t, x), t, start_point, alpha);

% Our data
plot3(xr(:, 1), xr(:, 2), xr(:, 3));
title(sprintf('Reconstructed %s system with \\alpha = [%s]', func2str(sys), join(string(alpha), ', ')));
legend('Original data', 'Reconstructed data');