function plotAmplitudeSummary(analysis_date, n_back_list, opts)
% plotAmplitudeSummary
%
% Poster Figure 3 - DoG amplitude vs n-back across the full 3x3 design.
% One figure is written per condition. Rows are previous level, columns are
% current level, and each cell plots amplitude across n-back with CI bars.
%
% Error bars come from sd_ci.lo / sd_ci.hi.
%
% Usage:
%   plotAmplitudeSummary('11.24.2025', [1 2 3]);

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

    estimates = loadEstimates(analysis_date, n_back_list, opts);

    id_transform = @(v) v;
    id_ci       = @(lo, hi) deal(lo, hi);

    plotSdSummary(estimates, n_back_list, plt_opts, p, opts, ...
        1, 'Amplitude (°)', 'Poster Fig3 Amplitude Summary', ...
        id_transform, id_ci);
end

function p = defaultLabels()
    p.cond_names = {'Contrast','Precision'};
    p.contrast   = {'90%','50%','25%'};
    p.precision  = {'2°','40°','80°'};
end
