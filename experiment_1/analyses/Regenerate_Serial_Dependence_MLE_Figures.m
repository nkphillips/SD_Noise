%%% Regenerate_Serial_Dependence_MLE_Figures
%
% Re-render serial-dependence MLE figures from a saved sd_noise analysis object.
% This script does not rerun trial-table construction, MLE, bootstrap, or
% jackknife. It only loads saved results and calls plotting functions.
%
% Usage: cd to experiment_1/analyses, edit the configuration below if needed,
% then run this script.

close all;
clc;

%% --- User configuration ---

% Leave empty to use the newest saved sd_noise file.
analysis_datetime = '';
sd_noise_file = '';
fit_method = 'hierarchical_map';   % 'pooled' | 'hierarchical_map'

% false: write into figures/<analysis_datetime>/regen_<timestamp>/
% true : overwrite/regenerate figures under figures/<analysis_datetime>/
overwrite_original_figures = false;

regenerate_super_subject = true;
regenerate_subject = true;
regenerate_folded_delta_theta = true;
show_middle_level_endpoint = false;  % faded, unconnected L2 reference dot in endpoint-effect figures
ci_figure_methods = {'percentile', 'bca'};  % writes percentile_ci/ and bca_ci/ figure subfolders

%% Paths

addpath('functions');
addpath(fullfile('functions', 'serial_dependence_mle'));
addpath(fullfile('functions', 'serial_dependence_mle', 'lapse_sensitivity'));
addpath(fullfile('plotting', 'serial_dependence_mle'));
addpath(fullfile('plotting', 'serial_dependence_mle', 'lapse_sensitivity'));
addpath('plotting');

ps = plotSettings();

%% Locate and load sd_noise

if isempty(sd_noise_file)
    estimates_root = fullfile('estimates', char(fit_method));
    if isempty(analysis_datetime)
        files = dir(fullfile(estimates_root, '*', 'SD_Noise_SerialDependence_sd_noise_*.mat'));
    else
        files = dir(fullfile(estimates_root, analysis_datetime, 'SD_Noise_SerialDependence_sd_noise_*.mat'));
    end
    if isempty(files)
        error('Regenerate_Serial_Dependence_MLE_Figures:noSavedAnalysis', ...
            'No saved sd_noise files found under %s.', estimates_root);
    end
    [~, ord] = sort([files.datenum], 'descend');
    files = files(ord);
    sd_noise_file = fullfile(files(1).folder, files(1).name);
end

S = load(sd_noise_file, 'sd_noise');
if ~isfield(S, 'sd_noise')
    error('Regenerate_Serial_Dependence_MLE_Figures:missingSdNoise', ...
        'File does not contain variable sd_noise: %s', sd_noise_file);
end
sd_noise = S.sd_noise;

fprintf('Loaded sd_noise: %s\n', sd_noise_file);

%% Output location

if overwrite_original_figures
    output_root = sd_noise.paths.figure_root;
else
    regen_datetime = datestr(now, 'yyyymmdd_HHMMSS');
    output_root = fullfile(sd_noise.paths.figure_root, ['regen_' regen_datetime]);
end

if ~exist(output_root, 'dir')
    mkdir(output_root);
end

regen_opts = struct( ...
    'output_root', output_root, ...
    'ps', ps, ...
    'n_back_list', sd_noise.config.n_back_list, ...
    'regenerate_folded_delta_theta', regenerate_folded_delta_theta, ...
    'show_middle_level_endpoint', show_middle_level_endpoint, ...
    'ci_figure_methods', {ci_figure_methods});

%% Regenerate figures

if regenerate_super_subject
    regenSuperSubjectFigs(sd_noise, regen_opts);
end

if regenerate_subject
    regenSubjectFigs(sd_noise, regen_opts);
end

fprintf('Regenerated figures under: %s\n', output_root);
