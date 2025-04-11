%%% SD_SF_Adapt_Figures

%% Prepare workspace

clc; close all; clear all;

which_setup = 'simulated';
analysis_date = '11.18.2024';

init_paths

subj_IDs = {'900' '901' '902' '903' '904' '905' '906' '907' '908' '909'};  
cond_names = {'contrast' 'filter'};

num.subjs = length(subj_IDs);
num.conds = length(cond_names);

%% Toggles

plot_ind = 1;
plot_grp = 1;

save_ind_figures = 1;
save_grp_figures = 1;

estimate_response_bias = 1;
estimate_sd_bias = 1;

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

%% Hard coded variables


%% Load experiment runs

load_matfiles

%% Analyze behavioral data

init_behav_perf

%% Stats

test_performance

%% plot info


%% Subject Performance

if plot_ind

    fg  = figure('Color',figure_color);
    set(0, 'CurrentFigure', fg)

    for subj = 1:num.subjs

        plot_performance_ind

    end

    close(fg)

end

%% Group Performance

if plot_grp

    fg  = figure('Color',figure_color);
    set(0, 'CurrentFigure', fg)

    plot_performance_grp

    close(fg)

end