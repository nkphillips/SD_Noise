function fwhm = unbinnedWtoFwhm(w)
% unbinnedWtoFwhm  Map S&S DoG width w (1/deg) to Gaussian-envelope FWHM in deg.
%
% After reparameterization to Sheehan & Serences (2022) Eq. 1
% (y = x*A*w*c*exp(-(w*x)^2), c = sqrt(2e), w in 1/deg), the Gaussian-envelope
% FWHM is 2*sqrt(log(2))/w. This matches the convention used in the legacy
% bootstrap helper bootstrapSuperSubjectSerialDependenceBySubject.m and the
% ~42.9 deg behavioral FWHM reported by S&S 2022 (peak A at 1/(w*sqrt(2))).

    w = w(:);
    fwhm = nan(size(w));
    ok = isfinite(w) & w > 0;
    fwhm(ok) = (2 * sqrt(log(2))) ./ w(ok);

end
