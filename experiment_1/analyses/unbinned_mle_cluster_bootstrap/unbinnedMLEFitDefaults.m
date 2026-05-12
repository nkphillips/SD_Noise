function d = unbinnedMLEFitDefaults()
% unbinnedMLEFitDefaults  Central defaults for unbinned DoG + psychometric MLE.
%
% Parameters are ordered as [A; w; sigma; beta]. FWHM bounds are expressed
% in degrees under the numerically estimated S&S DoG-lobe FWHM convention.

    d.fwhm_min_deg = 8;
    d.fwhm_max_deg = 90;

    d.A_ub = 30;
    d.A_lb = -d.A_ub;

    d.sigma_lb = 1;
    d.sigma_ub = 90;

    d.beta_ub = 10;
    d.beta_lb = -d.beta_ub;

    d.w_lb = unbinnedFwhmToW(d.fwhm_max_deg);
    d.w_ub = unbinnedFwhmToW(d.fwhm_min_deg);

    d.lb = [d.A_lb; d.w_lb; d.sigma_lb; d.beta_lb];
    d.ub = [d.A_ub; d.w_ub; d.sigma_ub; d.beta_ub];
    d.x0 = [1; 0.1; 5; 0];

    d.fit_opts = struct('lb', d.lb, 'ub', d.ub, 'x0', d.x0);
end
