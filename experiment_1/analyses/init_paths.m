%%% init_paths

%% Set directories

script_dir = pwd;

data_dir = '../data';
func_dir = 'functions';
estimates_path = 'estimates/';
reports_path = 'reports/';

figures_folder = 'figures/';
plt_opts.ind_figure_path = [figures_folder 'subjects/' analysis_date];
plt_opts.grp_figure_path = [figures_folder 'group/' analysis_date];
plt_opts.sup_figure_path = [figures_folder 'super/' analysis_date];
figure_file_type = '.pdf';

%% Check/Add paths

addpath(func_dir);

if plt_opts.save_ind_figures
    if ~exist(plt_opts.ind_figure_path,'dir'), mkdir(plt_opts.ind_figure_path); end
    addpath(plt_opts.ind_figure_path);
end

if plt_opts.save_grp_figures
    if ~exist(plt_opts.grp_figure_path,'dir'), mkdir(plt_opts.grp_figure_path); end
    addpath(plt_opts.grp_figure_path);
end

if plt_opts.save_sup_figures
    if ~exist(plt_opts.sup_figure_path,'dir'), mkdir(plt_opts.sup_figure_path); end
    addpath(plt_opts.sup_figure_path);
end

if ~exist('init/','dir'), mkdir('init/'); end
addpath('init/');

if ~exist('loading/','dir'), mkdir('loading/'); end
addpath([script_dir '/loading']);

if ~exist('plotting/','dir'), mkdir('plotting/'); end
addpath('plotting/');

if ~exist('tests/','dir'), mkdir('tests/'); end
addpath([script_dir '/tests']);

% Ensure estimates directory exists
if ~exist(estimates_path, 'dir')
    mkdir(estimates_path);
end

% Ensure reports directory exists
if ~exist(reports_path, 'dir')
    mkdir(reports_path);
end

