function y = dogIsolated(dtheta, A, w)
% dogIsolated  Sheehan & Serences (2022) DoG, no baseline.
%
%   y = dtheta * A * w * c * exp(-(w*dtheta)^2),  c = sqrt(2*e)
%
% Matches Eq. 1 of Sheehan & Serences 2022 PLoS Biology and the legacy
% calcDoG.m parameterization (w in 1/deg). Peak amplitude = A at
% dtheta = 1/(w*sqrt(2)); reported FWHM is the numerically estimated
% DoG-lobe width at half peak height.
%
% dtheta : column vector (deg). A, w : scalars (or conformable for broadcasting).

    dtheta = dtheta(:);
    c = sqrt(2 * exp(1));
    y = (dtheta .* A .* w .* c) .* exp(-(w .* dtheta).^2);

end
