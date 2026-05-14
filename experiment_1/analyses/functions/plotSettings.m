function ps = plotSettings()
% plotSettings  Central plot options for experiment_1 analyses (SD_Noise_Analyses_And_Figures).
% Figure directory fields (ind_figure_path, grp_figure_path, sup_figure_path) are filled by
% init_paths after analysis_date is known.

ps.plot_ind_figures   = 0;
ps.plot_grp_figures   = 0;
ps.plot_sup_figures   = 1;
ps.plot_ci_diagnostics = 1; % saves bootstrap CI diagnostics under the super-subject CI folder

ps.save_ind_figures   = 0;
ps.save_grp_figures   = 0;
ps.save_sup_figures   = 1;
ps.save_ci_diagnostics = 1;

ps.axis_square = 1;

ps.tick_length       = 0.020;
ps.line_width        = 1;
ps.marker_size       = 0;
ps.marker_size_scatter = 50;
ps.marker_size_bin   = 5;
ps.marker_size_polarplot = 350;

ps.colors.blue   = [38 71 237]/255;
ps.colors.red    = [204 0 0]/255;
ps.colors.green  = [0 153 0]/255;
ps.colors.black  = [0 0 0];
ps.colors.white  = [1 1 1];
ps.colors.purple = [102 51 204]/255;
ps.colors.orange = [255 128 0]/255;
ps.colors.yellow = [242 214 53]/255;
ps.colors.gray   = [128 128 128]/255;

ps.font_type = 'Helvetica';
ps.axes_label_font_size = 14;
ps.axes_tick_font_size = 13;

ps.figure_color = ps.colors.white;
ps.alpha_lvl    = 0.75;
ps.fg_type      = 'pdf';

ps.sd_colorbar_global = 0; % SD grids: 1=global min/max, 0=per-subplot
ps.rb_subtract_baseline = 0; % response bias: subtract DoG baseline (b) if provided

end
