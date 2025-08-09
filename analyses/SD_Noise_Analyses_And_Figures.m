%%% SD_Noise_Analyses_And_Figures_v2

%% Prepare workspace

close all; 
clear all; 
clc; 

%% Toggles

toggles.parallelization = 1;
toggles.sd_objective = 'sse'; % minimize 'nll' or 'sse' for serial dependence estimation
toggles.disp_on = 1;
toggles.save_estimates = 1;

%% Plot settings

plt_settings.plot_ind_figures = 1;
plt_settings.plot_grp_figures = 1;
plt_settings.plot_sup_figures = 1;
plt_settings.save_ind_figures = 1;
plt_settings.save_grp_figures = 1;
plt_settings.save_sup_figures = 1;

plt_settings.axis_square = 1;

plt_settings.tick_length = 0.020;
plt_settings.line_width = 1;
plt_settings.marker_size = 0;
plt_settings.marker_size_scatter = 50;
plt_settings.marker_size_bin = 5;
plt_settings.marker_size_polarplot = 350;

plt_settings.colors.blue = [0 76 152]/255;
plt_settings.colors.red = [204 0 0]/255;
plt_settings.colors.green = [0 153 0]/255;
plt_settings.colors.black = [0 0 0];
plt_settings.colors.white = [1 1 1];
plt_settings.colors.purple = [102 51 204]/255;
plt_settings.colors.orange = [255 128 0]/255;
plt_settings.colors.yellow = [242 214 53]/255;
plt_settings.colors.gray = [128 128 128]/255;

plt_settings.figure_color = plt_settings.colors.white;
plt_settings.alpha_lvl = 0.75;
plt_settings.fg_type = 'pdf';

%% Open figure handle

if plt_settings.plot_ind_figures || plt_settings.plot_sup_figures
    fg = figure('Visible','on','Color','w');
    set(0,'CurrentFigure',fg); 
end

%% Hard-coded vars

which_setup = '3329C_ASUS';
analysis_date = datestr(now, 'mm.dd.yyyy'); % automatically pull current date from system

p.subj_IDs = {'001', '002' '003', '004' '006', '007', '008'};  
p.cond_names = {'Contrast' 'Precision'};

% Define contrast and precision values for axis labels
p.contrast = {'90%', '50%', '25%'}; % high, medium, low contrast
p.precision = {'2°', '40°', '80°'}; % low, medium, high precision

num.subjs = length(p.subj_IDs);
num.conds = length(p.cond_names);
num.blocks = 6;
num.levels = 3;
num.blocks_per_cond = num.blocks/num.conds;
num.cond_combos = num.conds * num.levels;

%% Initialize paths and load experiment runs

init_paths
load_matfiles

%% Count trials across all subjects

if toggles.disp_on == 1
   
    disp(' ')
    disp('Counting trials across all subjects...')

end

% Initialize trial count arrays
trial_counts = nan(num.levels, num.levels, num.conds, num.subjs); % [prev_lvl, curr_lvl, cond, subj]
total_trials_per_cond = nan(num.conds, num.subjs); % [cond, subj]

% Count trials for each subject
for subj = 1:num.subjs
    subj_total_trials = 0;
    
    for n_run = 1:num.runs(subj)
        subj_p = all_runs{subj}(n_run).p;
        
        for cond = 1:num.conds
            curr_cond_blocks = find(subj_p.cond_order == cond);
            
            for n_block = 1:length(curr_cond_blocks)
                curr_lvls = subj_p.trial_events(:,3,curr_cond_blocks(n_block));
                
                for prev_lvl = 1:num.levels
                    for curr_lvl = 1:num.levels
                        curr_lvl_pair_indx = curr_lvls(1:end-1) == prev_lvl & curr_lvls(2:end) == curr_lvl;
                        trial_counts(prev_lvl, curr_lvl, cond, subj) = sum(trial_counts(prev_lvl, curr_lvl, cond, subj), 'omitnan') + sum(curr_lvl_pair_indx);
                    end
                end
                
                % Count total trials for this condition
                total_trials_per_cond(cond, subj) = sum(total_trials_per_cond(cond, subj), 'omitnan') + length(curr_lvls) - 1; % -1 because first trial has no previous trial
            end
        end
        
        subj_total_trials = subj_total_trials + size(subj_p.trial_events, 1) - 1; % -1 because first trial has no previous trial
    end
end

if toggles.disp_on == 1
    disp(' ');

    % Display trial count summary
    disp('=== TRIAL COUNT SUMMARY ===');
    disp(['Total subjects: ' num2str(num.subjs)]);
    disp(['Total runs: ' num2str(sum(num.runs))]);

    % Display per-subject totals
    disp(' ');
    disp('Trials per subject:');
    for subj = 1:num.subjs
        subj_total = sum(total_trials_per_cond(:, subj), 'omitnan');
        disp(['  Subject ' p.subj_IDs{subj} ': ' num2str(subj_total) ' trials (' num2str(num.runs(subj)) ' runs)']);
    end

    % Display per-condition totals
    disp(' ');
    disp('Trials per condition (across all subjects):');
    for cond = 1:num.conds
        cond_total = sum(total_trials_per_cond(cond, :), 'omitnan');
        cond_name = p.cond_names{cond};
        disp(['  ' cond_name ': ' num2str(cond_total) ' trials']);
    end

    % Display level combination totals
    disp(' ');
    disp('Trials per level combination (across all subjects):');
    for cond = 1:num.conds
        disp(['  ' p.cond_names{cond} ':']);
        % Columns represent current level; rows represent previous level
        for curr_lvl = 1:num.levels
            for prev_lvl = 1:num.levels
                combo_total = sum(trial_counts(prev_lvl, curr_lvl, cond, :), 'omitnan');
                if cond == 1
                    prev_label = p.contrast{prev_lvl};
                    curr_label = p.contrast{curr_lvl};
                else
                    prev_label = p.precision{prev_lvl};
                    curr_label = p.precision{curr_lvl};
                end
                disp(['    ' prev_label ' -> ' curr_label ': ' num2str(combo_total) ' trials']);
            end
        end
    end

    disp('================================');
end

%% Sort data into delta theta windows

delta_theta_width = 32;
delta_theta_steps = 1;
delta_theta_centers = -90:delta_theta_steps:90;

num.delta_theta_windows = length(delta_theta_centers);

delta_theta_windows = makeDeltaThetaWindows(delta_theta_centers, delta_theta_width, all_runs, num, p, plt_settings);

%% Define model bounds and parameters

p.fmincon_options = optimoptions('fmincon','Display','off');

%%% Response bias model parameters %%%
p.rb_init_params = [0, 1]; % [mu, sigma]
p.rb_bounds = [20, 30; -20, 0.1]; % [upper; lower]
p.guess_rate = 0.25; % Assuming constant guess rate (Sheehan & Serences 2022 PLOS Biology)

%%% Serial dependence model parameters %%%
p.sd_init_params = [1, 0.1, 0, 1]; % [amplitude, width (1/deg), baseline, sigma]

% Set width bounds via FWHM bounds (more interpretable), then convert to w = 1.6651 / FWHM
fwhm_min_deg = 8; 
fwhm_max_deg = 120; 
w_lb = 1.6651 / fwhm_max_deg; 
w_ub = 1.6651 / fwhm_min_deg;

p.sd_bounds = [p.rb_bounds(1,1), w_ub, 5, p.rb_bounds(1,2); 4, w_lb, -5, p.rb_bounds(2,2)]; % [upper; lower]
p.sd_objective = toggles.sd_objective; % 'nll' or 'sse' for serial dependence estimation

if strcmp(p.sd_objective, 'sse')
    p.sd_init_params = p.sd_init_params(1:3);
    p.sd_bounds = p.sd_bounds(:, 1:3);
end

%% Setup parallelization

p.num_workers = 9;
p.num_chunks = p.num_workers - 1;

if toggles.parallelization
    maxNumCompThreads(p.num_workers);
    if p.num_workers > 1
        parpool('local', p.num_chunks);
    end
end

%% Estimate response bias

if toggles.disp_on, disp(' '); disp('Estimating response bias for super subject...'); end

overall_start_time = tic;
rb_start_time = tic;

%%% Initialize structs %%%

rb = makeResponseBias(num);

%%% Create linear task list from 4D cell array %%%

if toggles.disp_on, disp('Creating task list from delta theta windows...'); end

task_list = {};
task_indices = [];

task_count = 0;
for prev_lvl = 1:num.levels
    for curr_lvl = 1:num.levels
        for cond = 1:num.conds
            for delta_theta_window_indx = 1:num.delta_theta_windows
                probe_offsets = delta_theta_windows.all.probe_offsets{prev_lvl, curr_lvl, cond, delta_theta_window_indx};
                responses = delta_theta_windows.all.responses{prev_lvl, curr_lvl, cond, delta_theta_window_indx};
                
                % Only include tasks with data
                if ~isempty(probe_offsets) && ~isempty(responses)
                    task_list{end+1} = {probe_offsets, responses};
                    task_indices = [task_indices; prev_lvl, curr_lvl, cond, delta_theta_window_indx];
                    task_count = task_count + 1;
                end 

            end
        end
    end
end

num_tasks = length(task_list);

if toggles.disp_on 
    disp(['✓ Task list created: ' num2str(num_tasks) ' tasks with data']); 
    disp(['  - Filtered out ' num2str(num.levels * num.levels * num.conds * num.delta_theta_windows - num_tasks) ' empty windows']); 
    disp(' '); 
end

%%% Process tasks using unified chunking function %%%

if toggles.disp_on
    disp('Starting response bias estimation...');
    if toggles.parallelization
        disp(['  - Using parallel processing with ' num2str(p.num_chunks) ' chunks']);
    else
        disp('  - Using sequential processing');
    end
end

[all_results, ~] = processTasksInChunks(task_list, p.num_chunks, toggles.parallelization, @processResponseBiasTask, toggles, p);

if toggles.disp_on, disp('✓ Response bias estimation completed'); end

%%% Extract results from the unified output %%%

if toggles.disp_on, disp('Extracting and compiling results...'); end

all_start_params = nan(num_tasks, 2);
all_start_nll = nan(num_tasks, 1);
all_params_est = nan(num_tasks, 2);
all_nll = nan(num_tasks, 1);
all_null_nll = nan(num_tasks, 1);
all_better_than_null = nan(num_tasks, 1);
all_exit_flag = nan(num_tasks, 1);

for i_task = 1:num_tasks
    all_start_params(i_task,:) = all_results{i_task}.start_params;
    all_start_nll(i_task) = all_results{i_task}.start_nll;
    all_params_est(i_task,:) = all_results{i_task}.params_est;
    all_nll(i_task) = all_results{i_task}.nll;
    % null_nll and better_than_null are computed later per window; keep as NaN here
    all_exit_flag(i_task) = all_results{i_task}.exit_flag;
end

if toggles.disp_on, disp('✓ Results extracted'); end

%%% Reconstruct original 4D structure %%%

if toggles.disp_on, disp('Reconstructing 4D response bias structure...'); end

for i_task = 1:num_tasks
    prev_lvl = task_indices(i_task, 1);
    curr_lvl = task_indices(i_task, 2);
    cond = task_indices(i_task, 3);
    delta_theta_window_indx = task_indices(i_task, 4);    

    % Get the original data for null NLL calculation
    probe_offsets = delta_theta_windows.all.probe_offsets{prev_lvl, curr_lvl, cond, delta_theta_window_indx};
    
    % Compile response bias
    rb.all.start_params(prev_lvl, curr_lvl, cond, delta_theta_window_indx,:) = all_start_params(i_task,:);
    rb.all.start_nll(prev_lvl, curr_lvl, cond, delta_theta_window_indx) = all_start_nll(i_task);
    rb.all.params_est(prev_lvl, curr_lvl, cond, delta_theta_window_indx,:) = all_params_est(i_task,:);
    rb.all.nll(prev_lvl, curr_lvl, cond, delta_theta_window_indx) = all_nll(i_task);
    rb.all.null_nll(prev_lvl, curr_lvl, cond, delta_theta_window_indx) = - size(probe_offsets,1) * log(0.5); % - N * log(0.5) for 2AFC random guessing
    rb.all.better_than_null(prev_lvl, curr_lvl, cond, delta_theta_window_indx) = all_nll(i_task) < rb.all.null_nll(prev_lvl, curr_lvl, cond, delta_theta_window_indx);
    rb.all.exit_flag(prev_lvl, curr_lvl, cond, delta_theta_window_indx) = all_exit_flag(i_task);
end

rb_duration = toc(rb_start_time);

if toggles.disp_on, disp('✓ Response bias structure reconstructed'); end

%%% Display summary %%%

if toggles.disp_on
    
    disp(' ');
    disp('=== RESPONSE BIAS ANALYSIS SUMMARY ===');
    disp(['Time: ~' num2str(round(rb_duration/60, 1)) ' minutes']);
    disp(['Tasks processed: ' num2str(num_tasks) ' estimations']);
    disp(['Data windows: ' num2str(num.delta_theta_windows) ' delta theta windows']);
    disp(['Conditions: ' num2str(num.conds) ' (Contrast, Precision)']);
    disp(['Level combinations: ' num2str(num.levels * num.levels) ' (3×3)']);

    if toggles.parallelization
        disp('Processing mode: Parallel');
        disp(['Chunks used: ' num2str(p.num_chunks)]);
    else
        disp('Processing mode: Sequential');
    end
    disp('=====================================');

end

%% Estimate serial dependence

if toggles.disp_on
    disp(' ');
    disp('Estimating serial dependence for super subject...');
end

%%% Initialize structs %%%

sd = makeSerialDependence(num);

if toggles.disp_on
    disp(' ');
    disp('Starting serial dependence estimation...');
end
    
sd_start_time = tic;

%%% Create linear task list from 4D cell array %%%
sd_task_list = {};
sd_task_indices = [];
sd_count = 0;

for cond = 1:num.conds
    for prev_lvl = 1:num.levels
        for curr_lvl = 1:num.levels
            sd_count = sd_count + 1;

            if strcmp(p.sd_objective, 'sse')
                % SSE mode: use observed mu per window from response-bias fits
                mu_per_window = squeeze(rb.all.params_est(prev_lvl, curr_lvl, cond, :, 1));
                % Determine which windows had data
                has_data = false(num.delta_theta_windows,1);
                for iw = 1:num.delta_theta_windows
                    has_data(iw) = ~isempty(delta_theta_windows.all.delta_thetas{prev_lvl, curr_lvl, cond, iw});
                end
                valid_windows = has_data & ~isnan(mu_per_window);
                if any(valid_windows)
                    x_dt = delta_theta_centers(valid_windows)';
                    y_mu = mu_per_window(valid_windows)';
                    
                    task_data.probe_offsets = zeros(size(x_dt)); % unused in SSE
                    task_data.responses = y_mu;                   % observed mu
                    task_data.delta_thetas = x_dt;                % window centers
                    task_data.condition_info = [prev_lvl, curr_lvl, cond];

                    sd_task_list{end+1} = task_data;
                    sd_task_indices = [sd_task_indices; prev_lvl, curr_lvl, cond];
                end
            else
                % NLL mode: use trial-level binary responses
                curr_probe_offsets = delta_theta_windows.all.probe_offsets(prev_lvl, curr_lvl, cond, :);
                curr_responses = delta_theta_windows.all.responses(prev_lvl, curr_lvl, cond, :);
                curr_delta_thetas = delta_theta_windows.all.delta_thetas(prev_lvl, curr_lvl, cond, :);

                probe_offsets = vertcat(curr_probe_offsets{:});
                responses = vertcat(curr_responses{:});
                delta_thetas = vertcat(curr_delta_thetas{:});

                if ~isempty(probe_offsets) && ~isempty(responses) && ~isempty(delta_thetas)
                    task_data.probe_offsets = probe_offsets;
                    task_data.responses = responses;
                    task_data.delta_thetas = delta_thetas;
                    task_data.condition_info = [prev_lvl, curr_lvl, cond];

                    sd_task_list{end+1} = task_data;
                    sd_task_indices = [sd_task_indices; prev_lvl, curr_lvl, cond];
                end
            end
        end
    end
end

num_sd_tasks = length(sd_task_list);
if toggles.disp_on
    disp(['  - Created ' num2str(num_sd_tasks) ' serial dependence tasks']);
end

if toggles.parallelization && num_sd_tasks > 1
    disp(['  - Using parallel processing with ' num2str(min(p.num_chunks, num_sd_tasks)) ' chunks']);
    [sd_results, ~] = processTasksInChunks(sd_task_list, min(p.num_chunks, num_sd_tasks), true, @processSerialDependenceTask, toggles, p);
else
    disp('  - Using sequential processing');
    [sd_results, ~] = processTasksInChunks(sd_task_list, 1, false, @processSerialDependenceTask, toggles, p);
end

%%% Compile results %%%

if toggles.disp_on, disp('  - Compiling serial dependence results...'); end

for i_task = 1:num_sd_tasks
    prev_lvl = sd_task_indices(i_task, 1);
    curr_lvl = sd_task_indices(i_task, 2);
    cond = sd_task_indices(i_task, 3);
    
    start_params = sd_results{i_task}.start_params;
    params_est = sd_results{i_task}.params_est;
    
    if strcmp(p.sd_objective, 'sse')
        sd.all.start_params(prev_lvl, curr_lvl, cond, 1:3) = start_params;
    else
        sd.all.start_params(prev_lvl, curr_lvl, cond, :) = start_params;
    end
    sd.all.start_nll(prev_lvl, curr_lvl, cond) = sd_results{i_task}.start_metric;
    if strcmp(p.sd_objective, 'sse')
        sd.all.params_est(prev_lvl, curr_lvl, cond, 1:3) = params_est;
    else
        sd.all.params_est(prev_lvl, curr_lvl, cond, :) = params_est;
    end
    sd.all.nll(prev_lvl, curr_lvl, cond) = sd_results{i_task}.final_metric;
    sd.all.exit_flag(prev_lvl, curr_lvl, cond) = sd_results{i_task}.exit_flag;

    % Calculate R² using appropriate observable
    if ~isempty(sd_task_list{i_task})
        if strcmp(p.sd_objective, 'nll')
            % Binary responses vs predicted p(CW)
            y_observed = sd_task_list{i_task}.responses; % binary 0/1
            dt = sd_task_list{i_task}.delta_thetas;
            po = sd_task_list{i_task}.probe_offsets;
            amp = sd_results{i_task}.params_est(1);
            w = sd_results{i_task}.params_est(2);
            base = sd_results{i_task}.params_est(3);
            sigma = sd_results{i_task}.params_est(4);
            mu = calcDoG(dt, [amp, w, base]);
            pCW = (1 - p.guess_rate) * normcdf(po, mu, sigma) + 0.5 * p.guess_rate;
            sd.all.r2(prev_lvl, curr_lvl, cond) = calcR2(y_observed, pCW);
        else
            % Window-level observed mu vs DoG-predicted mu
            y_observed = sd_task_list{i_task}.responses(:); % observed mu per window
            dt = sd_task_list{i_task}.delta_thetas(:);
            dog_params = sd_results{i_task}.params_est(1:3);
            y_fitted = calcDoG(dt, dog_params);
            sd.all.r2(prev_lvl, curr_lvl, cond) = calcR2(y_observed, y_fitted);
        end
    end

end

sd_duration = toc(sd_start_time);
if toggles.disp_on, disp('✓ Serial dependence estimation completed'); end

%%% Display summary %%%

if toggles.disp_on
    disp(' ');
    disp('=== SERIAL DEPENDENCE ANALYSIS SUMMARY ===');
    disp(['Time: ~' num2str(round(sd_duration/60, 1)) ' minutes']);
    disp(['Tasks processed: ' num2str(num_sd_tasks) ' conditions']);
    disp(['Conditions: ' num2str(num.conds) ' (Contrast, Precision)']);
    disp(['Level combinations: ' num2str(num.levels * num.levels) ' (3×3)']);
    if strcmp(p.sd_objective, 'sse')
        disp('Optimization metric: SSE (Sum of Squared Errors)');
    else
        disp('Optimization metric: NLL (Negative Log-Likelihood)');
    end
    if toggles.parallelization && num_sd_tasks > 1
        disp('Processing mode: Parallel');
    else
        disp('Processing mode: Sequential');
    end
    if toggles.parallelization && num_sd_tasks > 1
        disp(['Chunks used: ' num2str(min(p.num_chunks, num_sd_tasks))]);
    end

    % Display R² summary
    valid_r2 = sd.all.r2(~isnan(sd.all.r2));
    if ~isempty(valid_r2)
        disp(['R² summary: Mean = ' num2str(mean(valid_r2), '%.3f') ', Range = [' num2str(min(valid_r2), '%.3f') ', ' num2str(max(valid_r2), '%.3f') ']']);
    end
    disp('==========================================');
end

%% Plot results

if plt_settings.plot_sup_figures
    
    %% Performance (percent correct and percent CCW)

    plotPerformance(delta_theta_centers, delta_theta_windows.all.responses, delta_theta_windows.all.probe_offsets, p, plt_settings, 'Super Subj');
    clf(gcf);

    %% Response bias

    for cond = 1:num.conds

        fg_name = ['Super Subj Bias ' p.cond_names{cond}];

        % Columns represent current level; rows represent previous level
        for prev_lvl = 1:num.levels
            for curr_lvl = 1:num.levels

                % Create title with actual values
                if cond == 1
                    % Contrast condition
                    fg_title = [p.contrast{prev_lvl} ' -> ' p.contrast{curr_lvl}];
                else
                    % Precision condition
                    fg_title = [p.precision{prev_lvl} ' -> ' p.precision{curr_lvl}];
                end

                subplot(num.levels, num.levels, curr_lvl + (prev_lvl-1)*num.levels);
                set(0, 'CurrentFigure', fg);

                % Pass condition info to plotResponseBias for color selection
                plotResponseBias(delta_theta_centers, squeeze(rb.all.params_est(prev_lvl, curr_lvl, cond, :,1)), plt_settings, cond);
                
                % Format figure
                axis square;
                title(fg_title);
                
                % Only show x-axis labels on bottom row 
                if prev_lvl == num.levels
                    xlabel('\Delta\theta (°)');
                    xticks(-90:45:90);
                else
                    xticks(-90:45:90);
                    set(gca, 'XTickLabel', []);
                end
                
                % Only show y-axis label on left subplot
                if curr_lvl == 1 
                    ylabel('Bias (°)');
                end
                
                set(gca, 'TickDir', 'out', 'TickLength', [plt_settings.tick_length, plt_settings.tick_length]);
                xlim([-90 90]);
                line([min(xlim), max(xlim)], [0, 0], 'LineWidth', 1, 'Color', 'k');
                line([0, 0], [p.rb_bounds(2,1), p.rb_bounds(1,1)], 'LineWidth', 1, 'Color', 'k');
                ylim([p.rb_bounds(2,1) p.rb_bounds(1,1)]);
                yticks(p.rb_bounds(2,1):5:p.rb_bounds(1,1));
                box off;

                hold on;

            end
        end

        % Save figure
        if plt_settings.save_sup_figures
            saveas(gcf, fullfile(plt_settings.sup_figure_path, [fg_name '.' plt_settings.fg_type]));
        end

        clf(gcf);
    end

    %% Serial dependence

    % Create separate figures for each parameter
    if strcmp(p.sd_objective, 'sse')
        param_names = {'Amplitude', 'Width', 'Baseline'};
        num_params_to_plot = 3;
    else
        param_names = {'Amplitude', 'Width', 'Baseline', 'Sigma'};
        num_params_to_plot = 4;
    end
    
    for param_idx = 1:num_params_to_plot

        % Plot the parameter
        plotSerialDependence(sd.all.params_est, param_idx, param_names{param_idx}, p, plt_settings);
        
        % Save figure
        if plt_settings.save_sup_figures
            fg_name = ['Super Subj SD ' param_names{param_idx}];
            saveas(gcf, fullfile(plt_settings.sup_figure_path, [fg_name '.' plt_settings.fg_type]));
        end

        clf(gcf);
    end

    %% Response bias with best-fit DoG curves

    for cond = 1:num.conds

        fg_name = ['Super Subj Bias with DoG ' p.cond_names{cond}];

        % Columns represent current level; rows represent previous level
        for prev_lvl = 1:num.levels
            for curr_lvl = 1:num.levels

                % Create title with actual values
                if cond == 1
                    % Contrast condition
                    fg_title = [p.contrast{prev_lvl} ' -> ' p.contrast{curr_lvl}];
                else
                    % Precision condition
                    fg_title = [p.precision{prev_lvl} ' -> ' p.precision{curr_lvl}];
                end

                subplot(num.levels, num.levels, curr_lvl + (prev_lvl-1)*num.levels);
                set(0, 'CurrentFigure', fg);

                % Get DoG parameters from serial dependence analysis
                sd_params = squeeze(sd.all.params_est(prev_lvl, curr_lvl, cond, 1:3));
                
                % Plot response bias with DoG curve
                plotResponseBias(delta_theta_centers, squeeze(rb.all.params_est(prev_lvl, curr_lvl, cond, :, 1)), plt_settings, cond);
                
                % Plot DoG curve from existing estimates
                if ~isempty(sd_params)
                    
                    % Extract DoG parameters [amplitude, width, baseline]
                    dog_params = sd_params(1:3);
                    
                    % Generate smooth curve for plotting
                    delta_smooth = linspace(-90, 90, 100);
                    dog_fit = calcDoG(delta_smooth, dog_params);
                    
                    % Plot DoG fit with dashed line
                    if cond == 1
                        % Contrast condition - dashed blue
                        plot(delta_smooth, dog_fit, '--', 'LineWidth', plt_settings.line_width, 'Color', plt_settings.colors.black);
                    elseif cond == 2
                        % Precision condition - dashed green
                        plot(delta_smooth, dog_fit, '--', 'LineWidth', plt_settings.line_width, 'Color', plt_settings.colors.black);
                    end

                    % Annotate amplitude and FWHM (bottom-right)
                    amp_est = dog_params(1);
                    w_est = max(dog_params(2), eps);
                    fwhm_est = (2 * sqrt(log(2))) / w_est; % degrees
                    text(0.95, 0.10, sprintf('A = %.2f°\nFWHM = %.1f°\nb = %.2f°', amp_est, fwhm_est, dog_params(3)), ...
                        'Units', 'normalized', 'HorizontalAlignment', 'right', ...
                        'VerticalAlignment', 'bottom', 'FontWeight', 'normal', 'FontSize', 8, ...
                        'Color', plt_settings.colors.black);
                end

                % Format figure
                axis square;
                title(fg_title);
                
                % Only show x-axis labels on bottom row (prev level = last)
                if prev_lvl == num.levels
                    xlabel('\Delta\theta (°)');
                    xticks(-90:45:90);
                else
                    xticks(-90:45:90);
                    set(gca, 'XTickLabel', []);
                end
                
                % Only show y-axis label on bottom-left subplot (prev=last row, curr=first col)
                if curr_lvl == 1
                    ylabel('Bias (°)');
                end
                
                set(gca, 'TickDir', 'out', 'TickLength', [plt_settings.tick_length, plt_settings.tick_length]);
                xlim([-90 90]);
                line([min(xlim), max(xlim)], [0, 0], 'LineWidth', 1, 'Color', 'k');
                line([0, 0], [p.rb_bounds(2,1), p.rb_bounds(1,1)], 'LineWidth', 1, 'Color', 'k');
                ylim([p.rb_bounds(2,1) p.rb_bounds(1,1)]);
                yticks(p.rb_bounds(2,1):5:p.rb_bounds(1,1));
                
                % Annotate R^2 in top-left quadrant (smaller, regular weight)
                curr_r2 = sd.all.r2(prev_lvl, curr_lvl, cond);
                if ~isnan(curr_r2)
                    text(0.05, 0.90, sprintf('R^2 = %.2f', curr_r2), ...
                        'Units', 'normalized', 'HorizontalAlignment', 'left', ...
                        'VerticalAlignment', 'top', 'FontWeight', 'normal', 'FontSize', 8, ...
                        'Color', plt_settings.colors.black);
                end
                box off;
                hold on;

            end
        end

        % Save figure
        if plt_settings.save_sup_figures
            saveas(gcf, fullfile(plt_settings.sup_figure_path, [fg_name '.' plt_settings.fg_type]));
        end

        clf(fg);
    end
    
end

%% Clean up
 
if toggles.parallelization
    poolobj = gcp('nocreate');
    if ~isempty(poolobj)
        delete(poolobj);
    end
end

%% Overall analysis summary

overall_duration = toc(overall_start_time);
if toggles.disp_on
    disp(' ');
    disp('=== OVERALL ANALYSIS SUMMARY ===');
    disp(['Total time: ~' num2str(round(overall_duration/60, 1)) ' minutes']);
    disp(['Response bias: ' num2str(num_tasks) ' tasks completed (~' num2str(round(rb_duration/60, 1)) ' min)']);
    disp(['Serial dependence: ' num2str(num_sd_tasks) ' tasks completed (~' num2str(round(sd_duration/60, 1)) ' min)']);
    disp('================================');
end

%% Save workspace
% Store rb and sd structs into estimates folder 
% Check if file already exists, if so, present a GUI to ask if you want to overwrite

if toggles.save_estimates 
    
    % Get current time in HHMMSS format
    current_time = datestr(now, 'HHMMSS');
    
    filename = ['SD_Noise_Estimates_' analysis_date '_' current_time '_' toggles.sd_objective '.mat'];
    if exist(fullfile(estimates_path, filename), 'file')
        overwrite = questdlg('File already exists. Overwrite?', 'Overwrite?', 'Yes', 'No', 'Yes');
        if strcmp(overwrite, 'Yes')
            save(fullfile(estimates_path, filename), 'rb', 'sd');
            disp(['✓ Saved estimates to: ' filename]);
        end
    else
        save(fullfile(estimates_path, filename), 'rb', 'sd');
        disp(['✓ Saved estimates to: ' filename]);
    end

end