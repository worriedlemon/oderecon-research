function H = sindy_real(t, x, sigma, lambda)
%SINDY_REAL convenience function for SINDY()
    H = sindy(x, diff4(x, t), sigma, lambda);
end

