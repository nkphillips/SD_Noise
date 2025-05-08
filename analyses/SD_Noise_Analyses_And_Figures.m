%%% SD_SF_Adapt_Figures

%% Prepare workspace

clc; close all; clear all;

%% Toggles

plot_ind = 0;
plot_grp = 0;
plot_spr = 1;

save_ind_figures = 0;
save_grp_figures = 0;
save_spr_figures = 1;

estimate_response_bias = 1;
estimate_sd_bias = 1;

toggles.disp_on = 1;

%% Hard coded variables

which_setup = '3329C_ASUS';
analysis_date = '05.06.2025';

subj_IDs = {'001', '002'};  
cond_names = {'Contrast' 'Filter'};

num.subjs = length(subj_IDs);
num.conds = length(cond_names);
num.blocks = 6;
num.levels = 3;
num.blocks_per_cond = num.blocks/num.conds;
num.cond_combos = num.conds * num.levels;

% Model fitting options
options = optimoptions('fmincon');
lower_bounds = [-15, 1]; % lower bounds for [mu, sigma]
upper_bounds = [15, 20];

%% Initialize paths

init_paths

%% Load experiment runs

load_matfiles

unique_probe_offsets = unique(all_runs{1}(1).p.trial_events(:,2) - all_runs{1}(1).p.trial_events(:,1));

%% Plot settings

axis_square = 1;

tick_length = 0.020;
line_width = 1;
marker_size = 0;
marker_size_scatter = 50;
marker_size_bin = 5;
marker_size_polarplot = 350;

blue = [0 76 152]/255;
red = [204 0 0]/255;
green = [0 153 0]/255;
black = [0 0 0];
white = [1 1 1];
purple = [102 51 204]/255;
orange = [255 128 0]/255;

figure_color = white;
cond_colors = [green; blue; red];
alpha_lvl = 0.75;

% Calculate subplot dimensions 
[subplotX, subplotY] = get_subplot_dimensions(num.subjs);

%% Analyze behavioral data

analyzeData

%% Stats

% test_performance

%% plot info


%% Subject Performance

% if plot_ind

%     fg  = figure('Color',figure_color);
%     set(0, 'CurrentFigure', fg)

%     for subj = 1:num.subjs

%         plot_performance_ind

%     end

%     close(fg)

% end

%% Group Performance

% if plot_grp

%     fg  = figure('Color',figure_color);
%     set(0, 'CurrentFigure', fg)

%     plot_performance_grp

%     close(fg)

% end