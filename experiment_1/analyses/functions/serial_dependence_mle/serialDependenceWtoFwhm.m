function fwhm = serialDependenceWtoFwhm(w)
% serialDependenceWtoFwhm  Map S&S DoG width w (1/deg) to DoG-lobe FWHM in deg.
%
% After reparameterization to Sheehan & Serences (2022) Eq. 1
% (y = x*A*w*c*exp(-(w*x)^2), c = sqrt(2e), w in 1/deg), S&S report
% "the resulting FWHM estimated numerically." We follow that wording by
% solving for the width of the positive DoG lobe at half peak height:
%   sqrt(2e) * u * exp(-u^2) = 0.5, where u = w*x.

    w = w(:);
    fwhm = nan(size(w));
    ok = isfinite(w) & w > 0;
    fwhm(ok) = localDoGLobeFwhmCoefficient() ./ w(ok);

end

function coeff = localDoGLobeFwhmCoefficient()
    persistent cached_coeff
    if isempty(cached_coeff)
        half_height = @(u) sqrt(2 * exp(1)) .* u .* exp(-(u .^ 2)) - 0.5;
        peak_u = 1 / sqrt(2);
        try
            u1 = fzero(half_height, [eps, peak_u]);
            u2 = fzero(half_height, [peak_u, 3]);
            cached_coeff = u2 - u1;
        catch
            % Stable numerical solution for the S&S normalized DoG lobe.
            cached_coeff = 1.133150790039097;
        end
    end
    coeff = cached_coeff;
end
