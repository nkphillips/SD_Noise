%%% rerender_requested_figures
% Rerender figures affected by late plotting-only edits:
% - Poster Fig2/3/4 legends on both condition subplots
% - Super Subj SD Baseline Scatter error bars

clear; clc; close all;

try
    addpath('functions');
    addpath('plotting');

    analysis_date = '04.24.2026';
    n_back_list = [1 2 3];

    disp('Rerendering Poster Fig2 Sigma Summary...');
    plotSigmaSummary(analysis_date, n_back_list);

    disp('Rerendering Poster Fig3 Amplitude Summary...');
    plotAmplitudeSummary(analysis_date, n_back_list);

    disp('Rerendering Poster Fig4 FWHM Summary...');
    plotWidthSummary(analysis_date, n_back_list);

    disp('Loading estimates for baseline serial-dependence rerender...');
    estimates = loadEstimates(analysis_date, n_back_list);

    p = struct();
    p.cond_names = {'Contrast','Precision'};
    p.contrast = {'90%','50%','25%'};
    p.precision = {'2°','40°','80°'};
    p.num_chunks = 8;
    p.fmincon_options = optimoptions('fmincon','Display','off');
    p.rb_bounds = [20, 90; -20, 0.1];
    p.guess_rate = 0.25;
    p.sd_init_params = [1, 0.1, 0];
    fwhm_min_deg = 8;
    fwhm_max_deg = 120;
    w_lb = 1.6651 / fwhm_max_deg;
    w_ub = 1.6651 / fwhm_min_deg;
    p.sd_bounds = [p.rb_bounds(1,1), w_ub, 5; 1, w_lb, -5];
    p.sd_objective = 'sse';

    num = struct();
    num.levels = 3;
    num.conds = 2;

    bootstrap = struct();
    bootstrap.B = 10;
    bootstrap.ci = [2.5, 97.5];

    toggles = struct();
    toggles.disp_on = 1;
    toggles.parallelization = 0;

    plt_opts = posterPlotOpts();
    plt_opts.save_sup_figures = 1;
    plt_opts.fg_type = 'pdf';

    for n_back = n_back_list
        key = sprintf('n%d', n_back);
        est = estimates.(key);
        delta_theta_centers = est.delta_theta_centers;
        num.delta_theta_windows = numel(delta_theta_centers);

        if isempty(est.sd_ci)
            disp(sprintf('Recomputing sd_ci for %d-back...', n_back)); %#ok<DSPS>
            delta_theta_windows = makeFiniteWindowMask(est.rb, num);
            sd_ci = bootstrapSerialDependence(delta_theta_windows, ...
                delta_theta_centers, num, p, est.rb, bootstrap, toggles);
            save(est.file, 'sd_ci', '-append');
            est.sd_ci = sd_ci;
        end

        plt_opts.sup_figure_path = fullfile('figures', 'super', ...
            analysis_date, sprintf('%d_back', n_back));
        if ~exist(plt_opts.sup_figure_path, 'dir')
            mkdir(plt_opts.sup_figure_path);
        end

        disp(sprintf('Rerendering baseline scatter for %d-back...', n_back)); %#ok<DSPS>
        fg = figure('Visible','off','Color','w');
        plotSerialDependence(est.sd.all.params_est, 3, 'Baseline', p, ...
            plt_opts, fg, est.sd_ci.lo, est.sd_ci.hi);
        close(fg);
    end

    disp('Done.');
catch ME
    disp(getReport(ME, 'extended', 'hyperlinks', 'off'));
    rethrow(ME);
end

function delta_theta_windows = makeFiniteWindowMask(rb, num)
% Create the minimal delta_theta_windows.all.delta_thetas cell array needed
% for SSE-mode serial-dependence bootstrapping when reloading estimates from
% disk. The bootstrap only needs to know which windows had valid mu fits.

    delta_theta_windows = struct();
    delta_theta_windows.all = struct();
    delta_theta_windows.all.delta_thetas = cell(num.levels, num.levels, ...
        num.conds, num.delta_theta_windows);

    for cond = 1:num.conds
        for prev_lvl = 1:num.levels
            for curr_lvl = 1:num.levels
                mu = squeeze(rb.all.params_est(prev_lvl, curr_lvl, cond, :, 1));
                for iw = 1:num.delta_theta_windows
                    if isfinite(mu(iw))
                        delta_theta_windows.all.delta_thetas{prev_lvl, curr_lvl, cond, iw} = 1;
                    end
                end
            end
        end
    end
end
