load 'real_data/sprottb.mat'

[N, vc] = size(sprottb);
N = N / 2;
h = h * 2;
x = sprottb((1:N) * 2, :);

[Href, deg, eqc, mc] = load_href('SprottB');

t = (1:N) * h;

sigma = deglexord(deg, vc);

% LSM poly basis
B = EvalPoly(eye(mc), x, sigma);
delta = 1e-1;
B = (B'*B + delta*eye(mc))\B';

alpha_n = 50;
alpha_candidates = linspace(0.2, 0.9, alpha_n);

min_err = +Inf;
for alpha = alpha_candidates
    y = fractdiff(x, t, alpha);
    Htemp = sindy(x, y, sigma, 1);

    [Hc, Tc] = HTmat2cell(Htemp, sigma);
    [~, xr] = fractode(@(t,x)oderecon(Hc, Tc, t, x), t, x(1, :), alpha);
    if (any(isnan(xr) | isinf(xr)))
        fprintf(1, 'No good solution for alpha = %.3f\n', alpha);
    else
        err = nn_norm(t, x, xr);
        if (err < min_err)
            min_err = err;
            H = Htemp;
            xmin = xr;
        end
        fprintf(1, 'Error for alpha = %.3f: %.3f\n', alpha, err);
    end
end

plot3(x(:, 1), x(:, 2), x(:, 3));
hold on; grid on;
plot3(xmin(:, 1), xmin(:, 2), xmin(:, 3));