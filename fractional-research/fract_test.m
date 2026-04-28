func = @(x) x.^3;
x = linspace(0.1, 5, 500)';
fx = func(x);

% method
alphas = [0, 0.33, 0.66, 1, 1.5, 2];
for alpha = alphas
    da_x = fractdiff(fx, x, alpha);

    A = [ ones(size(x)), log(abs(x)) ];
    B = log(abs(da_x)); 
    idx = sum(isnan(A) | isinf(A) | isnan(B) | isinf(B), 2) == 0;
    A = A(idx, :);
    B = B(idx, :);

    X = (A'*A)\A'*B;
    k = exp(X(1));
    p = X(2);

    plot(x, da_x, DisplayName=sprintf('$\\alpha = %.2f, f(x) = %.2f \\cdot x^{%.2f}$', alpha, k, p));
    lgd = legend;
    lgd.Interpreter = 'latex';
    lgd.Visible = 'on';
    hold on; grid on;
end
