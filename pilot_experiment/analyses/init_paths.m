%%% init_paths

%% Set directories

script_dir = pwd;

data_dir = '../data';
func_dir = '../../functions';

figures_folder = 'figures/';
ind_figure_path = [figures_folder 'subjects/' analysis_date];
grp_figure_path = [figures_folder 'group/' analysis_date];
figure_file_type = '.pdf';

%% Check/Add paths

addpath(func_dir);

if ~exist(ind_figure_path,'dir'), mkdir(ind_figure_path); end
addpath(ind_figure_path);

if ~exist(grp_figure_path,'dir'), mkdir(grp_figure_path); end
addpath(grp_figure_path);

if ~exist('init/','dir'), mkdir('init/'); end
addpath('init/');

if ~exist('loading/','dir'), mkdir('loading/'); end
addpath([script_dir '/loading']);

if ~exist('plotting/','dir'), mkdir('plotting/'); end
addpath('plotting/');

if ~exist('tests/','dir'), mkdir('tests/'); end
addpath([script_dir '/tests']);

if ~exist('exploratory/','dir'), mkdir('exploratory/'); end
addpath([script_dir '/exploratory']);
