function [Xi, R] = get_norms(Href, t, x_ref, sigma, lambda_orth, obj)
    start_point = x_ref(1, :);
    [mc, eqc] = size(sigma);
    vc = eqc;
    T = mat2cell(repmat(sigma, 1, eqc), mc, repmat(vc, 1, eqc));

    H = obj(t, x_ref, sigma, lambda_orth);
    Xi = norm(H - Href);

    H = mat2cell(H, mc, ones(1, eqc));
    [~, x_new] = RK7(@(t, x)oderecon(H, T, t, x), t, start_point');
    R = nn_norm(t, x_new, x_ref);
end

