function plt_opts = posterPlotOpts()
% posterPlotOpts
%
% Shared plotting defaults used by the Phase 2 poster-oriented scripts
% (plotResponseBiasSigma, plotSigmaSummary, plotAmplitudeSummary,
%  plotWidthSummary, plotAmplitudeWidthScatter, plotBiasDogDiagonal).
%
% Mirrors the `ps` struct returned by functions/plotSettings.m (SD_Noise_Analyses_And_Figures).
% so that figures produced post-hoc match the aesthetic of the super-subject
% figures written during the main run.

    plt_opts = struct();

    plt_opts.axis_square   = 1;
    plt_opts.tick_length   = 0.020;
    plt_opts.line_width    = 1;
    plt_opts.marker_size   = 0;
    plt_opts.marker_size_scatter = 50;
    plt_opts.alpha_lvl     = 0.75;
    plt_opts.fg_type       = 'pdf';

    plt_opts.colors.blue   = [38 71 237]/255;
    plt_opts.colors.red    = [204 0 0]/255;
    plt_opts.colors.green  = [0 153 0]/255;
    plt_opts.colors.black  = [0 0 0];
    plt_opts.colors.white  = [1 1 1];
    plt_opts.colors.purple = [102 51 204]/255;
    plt_opts.colors.orange = [255 128 0]/255;
    plt_opts.colors.yellow = [242 214 53]/255;
    plt_opts.colors.gray   = [128 128 128]/255;

    plt_opts.figure_color        = plt_opts.colors.white;
    plt_opts.rb_subtract_baseline = 0;
    plt_opts.sd_colorbar_global  = 0;
end
