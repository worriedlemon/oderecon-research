function amp = snr2amp(x, snr)
    N = size(x, 1);
    amp = sqrt(sum(vecnorm(x, 2, 2) .^ 2) / N * 10.^(-snr ./ 10));
end

