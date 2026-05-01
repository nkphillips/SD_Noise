function plotWidthSummary(analysis_date, n_back_list, opts)
% plotWidthSummary
%
% Poster Figure 4 - DoG FWHM vs n-back across the full 3x3 design.
% Same structure as plotAmplitudeSummary but converts the width parameter
% w (1/deg) to Gaussian FWHM (deg): FWHM = 2*sqrt(ln 2) / w.
%
% Because FWHM is a monotone-decreasing function of w, the CI bounds are
% inverted: lo_FWHM = c / hi_w, hi_FWHM = c / lo_w.
%
% Usage:
%   plotWidthSummary('11.24.2025', [1 2 3]);

    if nargin < 3, opts = struct(); end
    if ~isfield(opts, 'objective'),        opts.objective        = 'sse'; end
    if ~isfield(opts, 'save'),             opts.save             = true;  end

    plt_opts = posterPlotOpts();
    p        = defaultLabels();

    this_dir = fileparts(mfilename('fullpath'));
    analyses_dir = fullfile(this_dir, '..');
    if ~isfield(opts, 'fig_dir') || isempty(opts.fig_dir)
        opts.fig_dir = fullfile(analyses_dir, 'figures', 'poster', analysis_date);
    end
    if opts.save && ~exist(opts.fig_dir, 'dir'); mkdir(opts.fig_dir); end

    if isfield(opts, 'estimates') && ~isempty(opts.estimates)
        estimates = opts.estimates;
    else
        estimates = loadEstimates(analysis_date, n_back_list, opts);
    end

    c = 2 * sqrt(log(2));
    w_to_fwhm    = @(w) c ./ max(w, eps);
    fwhm_ci_flip = @(lo_w, hi_w) deal(c ./ max(hi_w, eps), c ./ max(lo_w, eps));

    plotSdSummary(estimates, n_back_list, plt_opts, p, opts, ...
        2, 'FWHM (°)', 'Poster Fig4 FWHM Summary', ...
        w_to_fwhm, fwhm_ci_flip);
end

function p = defaultLabels()
    p.cond_names = {'Contrast','Precision'};
    p.contrast   = {'90%','50%','25%'};
    p.precision  = {'2°','40°','80°'};
end
