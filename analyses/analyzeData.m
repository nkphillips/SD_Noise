% analyzeData

%% Initialize analyzed_data structure

analyzed_data = struct('ind',[],'grp',[],'super',[]);

tmp_array = cell(1, num.subjs);
analyzed_data.ind = struct('counts',tmp_array,'performance_within_run',tmp_array,'performance_within_run_sem',tmp_array,'pCW_within_run',tmp_array,'pCW_within_run_sem',tmp_array);

tmp_array = cell(1);
analyzed_data.grp = struct('performance_within_run',tmp_array,'performance_within_run_sem',tmp_array,'pCW_within_run',tmp_array,'pCW_within_run_sem',tmp_array);

tmp_array = cell(1);
analyzed_data.super = struct('performance_within_run',tmp_array,'performance_within_run_sem',tmp_array,'pCW_within_run',tmp_array,'pCW_within_run_sem',tmp_array);

struct_fieldnames = fieldnames(analyzed_data.ind);

%% Initialize delta theta bins 

delta_theta_bin_width = 32;
delta_theta_bin_steps = 4;

delta_theta_bin_centers = (-180+delta_theta_bin_width/2):delta_theta_bin_steps:(180-delta_theta_bin_width/2);
num.delta_theta_bins = length(delta_theta_bin_centers);
delta_theta_bin_edges = [delta_theta_bin_centers - delta_theta_bin_width/2; delta_theta_bin_centers + delta_theta_bin_width/2];

delta_theta_bins = makeDeltaThetaBins(delta_theta_bin_edges, all_runs, unique_probe_offsets, num);

%% Initialize arrays for the super subject serial dependence

analyzed_data.super.serial_dependence_params_est = nan(num.levels, num.levels, num.conds, 2);
analyzed_data.super.serial_dependence_r2 = nan(num.levels, num.levels, num.conds);
analyzed_data.super.serial_dependence_estimated_bias = nan(num.levels, num.levels, num.conds, length(-180:180));

% Initialize cell arrays for storing data
analyzed_data.super.resp_per_offset_per_delta_theta_bin = cell(num.levels, num.levels, num.conds, length(unique_probe_offsets), num.delta_theta_bins);
analyzed_data.super.perf_per_delta_theta_bin = cell(num.levels, num.levels, num.conds, num.delta_theta_bins);
analyzed_data.super.pCW_per_delta_theta_bin = cell(num.levels, num.levels, num.conds, num.delta_theta_bins);
analyzed_data.super.mu_per_delta_theta_bin = cell(num.levels, num.levels, num.conds, num.delta_theta_bins);
analyzed_data.super.sigma_per_delta_theta_bin = cell(num.levels, num.levels, num.conds, num.delta_theta_bins);
analyzed_data.super.r2_per_delta_theta_bin = cell(num.levels, num.levels, num.conds, num.delta_theta_bins);

% Initialize cell array for storing all delta_thetas for the super subject
analyzed_data.super.all_delta_thetas = cell(num.levels, num.levels, num.conds);

% Initialize storage for super subject data accumulation
super_perf_per_delta_theta = cell(num.levels, num.levels, num.conds, num.delta_theta_bins);
super_pCW_per_delta_theta = cell(num.levels, num.levels, num.conds, num.delta_theta_bins);
super_responses_per_offset_delta_theta = cell(num.levels, num.levels, num.conds, length(unique_probe_offsets), num.delta_theta_bins);
all_delta_thetas = cell(num.levels, num.levels, num.conds);

%% Open figure handle

if plot_ind || plot_grp || plot_spr
    fg = figure('Visible','on','Color','w');
    set(0,'CurrentFigure',fg); 
end

%% Subject

for subj = 1:num.subjs

    if save_ind_figures
        figure_path = [ind_figure_path '/S' subj_IDs{subj}];
        if exist(figure_path, 'dir') == 0, mkdir(figure_path); end
    end

    % Initialize storage for performance and pCW within runs
    performance_mean_within_run = nan(num.levels, num.levels, num.conds, num.runs(subj));
    performance_sem_within_run = performance_mean_within_run;
    pCW_mean_within_run = nan(num.levels, num.levels, num.conds, num.runs(subj));
    pCW_sem_within_run = pCW_mean_within_run;

    % Initialize storage for pCW binned per offset and delta theta across runs
    pCW_per_probe_offset = cell(num.levels, num.levels, num.conds, length(unique_probe_offsets));
    pCW_per_probe_offset_mean = nan(num.levels, num.levels, num.conds, length(unique_probe_offsets));
    pCW_per_probe_offset_sem = pCW_per_probe_offset_mean;
    
    % Initialize storage for performance per delta theta bin
    perf_per_delta_theta_bin = cell(num.levels, num.levels, num.conds, num.delta_theta_bins);
    
    % Initialize storage for pCW binned per offset and delta theta across runs, adding a dimension for runs
    current_subj_pCW_binned_runs = nan(num.levels, num.levels, num.conds, length(unique_probe_offsets), num.delta_theta_bins, num.runs(subj));
    
    % Loop through runs
    for n_run = 1:num.runs(subj)

        % if toggles.disp_on, disp(['Run ' num2str(n_run)]); end

        subj_p = all_runs{subj}(n_run).p;
        behav_data = all_runs{subj}(n_run).behav_data;

        counts(n_run) = getConditionDistributions(subj_p, behav_data, unique_probe_offsets);

        delta_thetas = squeeze(subj_p.trial_events(1:end-1,1,:) - subj_p.trial_events(2:end,1,:));
        probe_offsets = squeeze(subj_p.trial_events(:,2,:) - subj_p.trial_events(:,1,:));

        for cond = 1:num.conds
            
            % if toggles.disp_on, disp(['Condition ' num2str(cond)]); end
            curr_cond_blocks = find(subj_p.cond_order == cond);

            for n_block = 1:length(curr_cond_blocks)

                curr_correct = behav_data.correct(:,curr_cond_blocks(n_block));
                curr_response = behav_data.response(:,curr_cond_blocks(n_block)) == 2;
                curr_delta_theta = delta_thetas(:,curr_cond_blocks(n_block));
                curr_probe_offsets = probe_offsets(:,curr_cond_blocks(n_block));
                curr_lvls = subj_p.trial_events(:,3,curr_cond_blocks(n_block));

                for prev_lvl = 1:num.levels
                    for curr_lvl = 1:num.levels

                        % if toggles.disp_on, disp(['Level Pair ' num2str(prev_lvl) '-' num2str(curr_lvl)]); end

                        curr_lvl_pair_indices = find( curr_lvls(1:end-1) == prev_lvl & curr_lvls(2:end) == curr_lvl );

  
                        curr_lvl_pair_correct = curr_correct(curr_lvl_pair_indices);
                        performance_mean_within_run(prev_lvl, curr_lvl, cond, n_run) = mean(curr_lvl_pair_correct);
                        performance_sem_within_run(prev_lvl, curr_lvl, cond, n_run) = std(curr_lvl_pair_correct)/sqrt(length(curr_lvl_pair_correct));

                        curr_lvl_pair_response_CW = curr_response(curr_lvl_pair_indices);
                        pCW_mean_within_run(prev_lvl, curr_lvl, cond, n_run) = mean(curr_lvl_pair_response_CW);
                        pCW_sem_within_run(prev_lvl, curr_lvl, cond, n_run) = std(curr_lvl_pair_response_CW)/sqrt(length(curr_lvl_pair_response_CW));

                        curr_lvl_pair_delta_theta = curr_delta_theta(curr_lvl_pair_indices);
                        curr_lvl_pair_probe_offset = curr_probe_offsets(curr_lvl_pair_indices);
                        
                        % Check if any data exists for this level pair/condition
                        if isempty(curr_lvl_pair_delta_theta)
                            continue; % Skip if no trials
                        end
    
                        % Accumulate delta_thetas for the super subject
                        all_delta_thetas{prev_lvl, curr_lvl, cond} = [all_delta_thetas{prev_lvl, curr_lvl, cond}; curr_lvl_pair_delta_theta];

                        for delta_theta_bin_indx = 1:num.delta_theta_bins

                            % if toggles.disp_on, disp(['Delta Theta bin ' num2str(delta_theta_bin_indx)]); end

                            curr_delta_theta_bin = [delta_theta_bin_edges(1, delta_theta_bin_indx), delta_theta_bin_edges(2, delta_theta_bin_indx)];
    
                            % Find trials within the current delta theta bin
                            delta_theta_mask = (curr_lvl_pair_delta_theta >= curr_delta_theta_bin(1)) & (curr_lvl_pair_delta_theta <= curr_delta_theta_bin(2));
    
                            % Accumulate for super subject
                            super_perf_per_delta_theta{prev_lvl, curr_lvl, cond, delta_theta_bin_indx} = [super_perf_per_delta_theta{prev_lvl, curr_lvl, cond, delta_theta_bin_indx}; ...
                                curr_lvl_pair_correct(delta_theta_mask)];
                            super_pCW_per_delta_theta{prev_lvl, curr_lvl, cond, delta_theta_bin_indx} = [super_pCW_per_delta_theta{prev_lvl, curr_lvl, cond, delta_theta_bin_indx}; ...
                                curr_lvl_pair_response_CW(delta_theta_mask)];
                            
                            % Also store for individual subjects
                            perf_per_delta_theta_bin{prev_lvl, curr_lvl, cond, delta_theta_bin_indx} = [perf_per_delta_theta_bin{prev_lvl, curr_lvl, cond, delta_theta_bin_indx}; ...
                                curr_lvl_pair_correct(delta_theta_mask)];
                            
                            for offset_indx = 1:length(unique_probe_offsets)

                                % if toggles.disp_on, disp(['Offset ' num2str(offset_indx)]); end

                                curr_offset = unique_probe_offsets(offset_indx);
    
                                % Find trials matching the current probe offset *within* the delta theta bin
                                offset_mask = (curr_lvl_pair_probe_offset == curr_offset);
                                combined_mask = delta_theta_mask & offset_mask;
    
                                % Get responses for the selected trials
                                responses_in_bin = curr_lvl_pair_response_CW(combined_mask);

                                % Accumulate for super subject
                                super_responses_per_offset_delta_theta{prev_lvl, curr_lvl, cond, offset_indx, delta_theta_bin_indx} = [super_responses_per_offset_delta_theta{prev_lvl, curr_lvl, cond, offset_indx, delta_theta_bin_indx}; ...
                                    responses_in_bin];
    
                                % Calculate pCW for this specific bin
                                if ~isempty(responses_in_bin)
                                    pCW_in_bin = mean(responses_in_bin, 'omitnan');
                                else
                                    pCW_in_bin = NaN; % Assign NaN if no trials in this bin
                                end
    
                                % Store the calculated pCW for this run in the numeric array
                                current_subj_pCW_binned_runs(prev_lvl, curr_lvl, cond, offset_indx, delta_theta_bin_indx, n_run) = pCW_in_bin;
    
                                if delta_theta_bin_indx == 1 % This seems unrelated to the delta_theta binning logic, might need review? Assumed correct for now.
                                pCW_per_probe_offset{prev_lvl, curr_lvl,cond,offset_indx} = [pCW_per_probe_offset{prev_lvl, curr_lvl,cond,offset_indx}, ...
                                    counts(n_run).pCW_per_probe_offset{prev_lvl, curr_lvl,cond, offset_indx}];
                                end
    
                            end
                        end

                    end
                end
            
            end

        end

    end

    % Store data in the analyzed_data structure
    analyzed_data.ind(subj).counts = counts;

    analyzed_data.ind(subj).performance_within_run = performance_mean_within_run;
    analyzed_data.ind(subj).performance_within_run_sem = performance_sem_within_run;
    analyzed_data.ind(subj).performance_across_runs = mean(analyzed_data.ind(subj).performance_within_run, 4);
    analyzed_data.ind(subj).performance_across_runs_sem = std(analyzed_data.ind(subj).performance_within_run, [], 4)/sqrt(num.runs(subj));

    analyzed_data.ind(subj).pCW_within_run = pCW_mean_within_run;
    analyzed_data.ind(subj).pCW_within_run_sem = pCW_sem_within_run;
    analyzed_data.ind(subj).pCW_across_runs = mean(analyzed_data.ind(subj).pCW_within_run, 4);
    analyzed_data.ind(subj).pCW_across_runs_sem = std(analyzed_data.ind(subj).pCW_within_run, [], 4)/sqrt(num.runs(subj));

    analyzed_data.ind(subj).pCW_per_offset_delta_theta_all_runs = current_subj_pCW_binned_runs; % Store raw per-run data
    analyzed_data.ind(subj).pCW_per_probe_offset = pCW_per_probe_offset;
    analyzed_data.ind(subj).perf_per_delta_theta_bin = perf_per_delta_theta_bin;

    % Transfer super subject data 
    analyzed_data.super.all_delta_thetas = all_delta_thetas;

    for prev_lvl = 1:num.levels
        for curr_lvl = 1:num.levels
            for cond = 1:num.conds
                
                
                for delta_theta_bin_indx = 1:num.delta_theta_bins
                    analyzed_data.super.perf_per_delta_theta_bin{prev_lvl, curr_lvl, cond, delta_theta_bin_indx} = [analyzed_data.super.perf_per_delta_theta_bin{prev_lvl, curr_lvl, cond, delta_theta_bin_indx}; ...
                        super_perf_per_delta_theta{prev_lvl, curr_lvl, cond, delta_theta_bin_indx}];
                    analyzed_data.super.pCW_per_delta_theta_bin{prev_lvl, curr_lvl, cond, delta_theta_bin_indx} = [analyzed_data.super.pCW_per_delta_theta_bin{prev_lvl, curr_lvl, cond, delta_theta_bin_indx}; ...
                        super_pCW_per_delta_theta{prev_lvl, curr_lvl, cond, delta_theta_bin_indx}];
                    for offset_indx = 1:length(unique_probe_offsets)
                        analyzed_data.super.resp_per_offset_per_delta_theta_bin{prev_lvl, curr_lvl, cond, offset_indx, delta_theta_bin_indx} = [analyzed_data.super.resp_per_offset_per_delta_theta_bin{prev_lvl, curr_lvl, cond, offset_indx, delta_theta_bin_indx}; ...
                            super_responses_per_offset_delta_theta{prev_lvl, curr_lvl, cond, offset_indx, delta_theta_bin_indx}];
                    end
                
                end

            end
        end
    end

    % Calculate mean and SEM for pCW binned per offset and delta theta across runs
    tmp_pCW_binned_across_runs = current_subj_pCW_binned_runs; % Use the numeric array containing data from all runs

    % Calculate mean across runs (dimension 6), ignoring NaNs
    analyzed_data.ind(subj).pCW_per_offset_delta_theta_mean = mean(tmp_pCW_binned_across_runs, 6, 'omitnan');

    % Calculate the number of non-NaN runs for each bin
    n_runs_with_data = sum(~isnan(tmp_pCW_binned_across_runs), 6); 

    % Calculate SEM across runs (dimension 6), ignoring NaNs
    analyzed_data.ind(subj).pCW_per_offset_delta_theta_sem = std(tmp_pCW_binned_across_runs, 0, 6, 'omitnan') ./ sqrt(n_runs_with_data);
    analyzed_data.ind(subj).pCW_per_offset_delta_theta_sem(n_runs_with_data < 2) = NaN; % Set SEM to NaN if less than 2 data points

    %% Fit response bias model to each delta theta bin
    
    if toggles.disp_on, disp('Fitting response bias model to each delta theta bin...'); end

    % Initialize storage for fitted parameters per delta theta bin
    mu_per_delta_theta = nan(num.levels, num.levels, num.conds, num.delta_theta_bins);
    sigma_per_delta_theta = nan(num.levels, num.levels, num.conds, num.delta_theta_bins);
    r2_per_delta_theta = nan(num.levels, num.levels, num.conds, num.delta_theta_bins);
    pCW_per_probe_offset_mean = nan(num.levels, num.levels, num.conds, length(unique_probe_offsets));
    pCW_per_probe_offset_sem = pCW_per_probe_offset_mean;

    for prev_lvl = 1:num.levels
        if toggles.disp_on, disp(['Previous LVL ' num2str(prev_lvl)]); end
        for curr_lvl = 1:num.levels
            if toggles.disp_on, disp(['Current LVL ' num2str(curr_lvl)]); end
            for cond = 1:num.conds
                if toggles.disp_on, disp([cond_names{cond} ' Condition']); end
                
                for offset_indx = 1:length(unique_probe_offsets)
                    pCW_per_probe_offset_mean(prev_lvl, curr_lvl, cond, offset_indx) = mean(pCW_per_probe_offset{prev_lvl, curr_lvl, cond, offset_indx},'omitnan');
                    pCW_per_probe_offset_sem(prev_lvl, curr_lvl, cond, offset_indx) = std(pCW_per_probe_offset{prev_lvl, curr_lvl, cond, offset_indx},'omitnan')/sqrt(num.runs(subj));
                end

                %% Fit model to pCW per probe offset for each delta theta bin
                for delta_theta_bin_indx = 1:num.delta_theta_bins
                    
                    curr_probe_offsets = unique_probe_offsets;

                    % Extract pCW for the current delta theta bin across all probe offsets
                    p_CW = squeeze(analyzed_data.ind(subj).pCW_per_offset_delta_theta_mean(prev_lvl, curr_lvl, cond, :, delta_theta_bin_indx));

                    rmv_rows = isnan(p_CW);
                    curr_probe_offsets(rmv_rows) = [];
                    p_CW(rmv_rows) = [];

                    % Skip fitting if there are NaN values (insufficient data for this bin)
                    if isempty(p_CW) || length(p_CW) < 3
                        if toggles.disp_on, disp(['Not enough data for delta theta bin ' num2str(delta_theta_bin_edges(1, delta_theta_bin_indx)) ' to ' num2str(delta_theta_bin_edges(2, delta_theta_bin_indx)) '°']); end
                        continue;
                    end

                    fixed_params = cell(1,2);
                    fixed_params{1} = [curr_probe_offsets, p_CW]; % Offsets and corresponding pCW
                    fixed_params{2} = p.guess_rate;

                    response_bias = estimateResponseBias(p.response_bias_init_params, fixed_params, p);

                    mu_per_delta_theta(prev_lvl, curr_lvl, cond, delta_theta_bin_indx) = response_bias.params_est(1);
                    sigma_per_delta_theta(prev_lvl, curr_lvl, cond, delta_theta_bin_indx) = response_bias.params_est(2);
                    r2_per_delta_theta(prev_lvl, curr_lvl, cond, delta_theta_bin_indx) = response_bias.r2;

                    clear fixed_params
                end 
            end 
        end 
    end 

    % Store all response bias parameters in the structure
    analyzed_data.ind(subj).mu_per_delta_theta = mu_per_delta_theta;
    analyzed_data.ind(subj).sigma_per_delta_theta = sigma_per_delta_theta;
    analyzed_data.ind(subj).r2_per_delta_theta = r2_per_delta_theta;
    analyzed_data.ind(subj).pCW_per_probe_offset_mean = pCW_per_probe_offset_mean;
    analyzed_data.ind(subj).pCW_per_probe_offset_sem = pCW_per_probe_offset_sem;
    
    % Estimate serial dependence for individual subject
    serial_dependence_params_est = nan(num.levels, num.levels, num.conds, 2);
    serial_dependence_r2 = nan(num.levels, num.levels, num.conds);
    serial_dependence_estimated_bias = nan(num.levels, num.levels, num.conds, length(-180:180));
    
    % Calculate performance per delta theta mean and SEM
    perf_per_delta_theta_bin_mean = nan(num.levels, num.levels, num.conds, num.delta_theta_bins);
    perf_per_delta_theta_bin_sem = nan(num.levels, num.levels, num.conds, num.delta_theta_bins);
    
    % Calculate pCW per delta theta bin mean and SEM (averaged across probe offsets)
    pCW_per_delta_theta_bin_mean = nan(num.levels, num.levels, num.conds, num.delta_theta_bins);
    pCW_per_delta_theta_bin_sem = nan(num.levels, num.levels, num.conds, num.delta_theta_bins);
    
    for prev_lvl = 1:num.levels
        for curr_lvl = 1:num.levels
            for cond = 1:num.conds
                % Performance calculations
                for win_idx = 1:num.delta_theta_bins
                    if ~isempty(perf_per_delta_theta_bin{prev_lvl, curr_lvl, cond, win_idx})
                        perf_per_delta_theta_bin_mean(prev_lvl, curr_lvl, cond, win_idx) = mean(perf_per_delta_theta_bin{prev_lvl, curr_lvl, cond, win_idx}, 'omitnan');
                        if sum(~isnan(perf_per_delta_theta_bin{prev_lvl, curr_lvl, cond, win_idx})) >= 2
                            perf_per_delta_theta_bin_sem(prev_lvl, curr_lvl, cond, win_idx) = std(perf_per_delta_theta_bin{prev_lvl, curr_lvl, cond, win_idx}, 'omitnan') / ...
                                sqrt(sum(~isnan(perf_per_delta_theta_bin{prev_lvl, curr_lvl, cond, win_idx})));
                        else
                            perf_per_delta_theta_bin_sem(prev_lvl, curr_lvl, cond, win_idx) = NaN;
                        end
                        
                    end
                    
                    % pCW calculations - average across probe offsets for each delta theta bin
                    if isfield(analyzed_data.ind(subj), 'pCW_per_offset_delta_theta_mean')
                        % Extract pCW for all probe offsets in this delta theta bin
                        curr_pCW = squeeze(analyzed_data.ind(subj).pCW_per_offset_delta_theta_mean(prev_lvl, curr_lvl, cond, :, win_idx));
                        
                        % Calculate mean pCW for this delta theta bin (across all probe offsets)
                        if any(~isnan(curr_pCW))
                            pCW_per_delta_theta_bin_mean(prev_lvl, curr_lvl, cond, win_idx) = mean(curr_pCW, 'omitnan');
                            
                            % Calculate SEM using the probe offset SEM values
                            curr_sem = squeeze(analyzed_data.ind(subj).pCW_per_offset_delta_theta_sem(prev_lvl, curr_lvl, cond, :, win_idx));
                            num_valid = sum(~isnan(curr_pCW));
                            if num_valid > 0
                                pCW_per_delta_theta_bin_sem(prev_lvl, curr_lvl, cond, win_idx) = sqrt(sum(curr_sem.^2, 'omitnan')) / num_valid;
                            end
                        end
                    end
                end
            end
        end
    end
    
    % Store statistics in the structure
    analyzed_data.ind(subj).perf_per_delta_theta_bin_mean = perf_per_delta_theta_bin_mean;
    analyzed_data.ind(subj).perf_per_delta_theta_bin_sem = perf_per_delta_theta_bin_sem;
    analyzed_data.ind(subj).pCW_per_delta_theta_bin_mean = pCW_per_delta_theta_bin_mean;
    analyzed_data.ind(subj).pCW_per_delta_theta_bin_sem = pCW_per_delta_theta_bin_sem;

    % Loop through levels and conditions to estimate serial dependence
    for cond = 1:num.conds
        for prev_lvl = 1:num.levels
            for curr_lvl = 1:num.levels
                % Get mu values across delta theta bins
                curr_delta_theta = delta_theta_bin_centers';
                curr_mu = squeeze(mu_per_delta_theta(prev_lvl, curr_lvl, cond, :));
                
                % Remove NaN values
                rmv_rows = isnan(curr_mu);
                curr_delta_theta(rmv_rows) = [];
                curr_mu(rmv_rows) = [];
                
                % Only skip if absolutely no data is available
                if isempty(curr_mu)
                    if toggles.disp_on
                        disp(['Subject ' subj_IDs{subj} ': No data for serial dependence model (Prev Level: ' num2str(prev_lvl) ', Curr Level: ' num2str(curr_lvl) ', Condition: ' cond_names{cond} ')']);
                    end
                    continue;
                end
                
                % Estimate serial dependence - passing data directly, not in a cell
                fixed_params = [curr_delta_theta, curr_mu];
                serial_dependence = estimateSerialDependence(p.serial_dependence_init_params, fixed_params, p);
                
                serial_dependence_params_est(prev_lvl, curr_lvl, cond, :) = serial_dependence.params_est;
                serial_dependence_r2(prev_lvl, curr_lvl, cond) = serial_dependence.r2;
                serial_dependence_estimated_bias(prev_lvl, curr_lvl, cond, :) = gaussianPrime(serial_dependence.params_est, -180:180);
            end
        end
    end
    
    % Store serial dependence parameters in the structure
    analyzed_data.ind(subj).serial_dependence_params_est = serial_dependence_params_est;
    analyzed_data.ind(subj).serial_dependence_r2 = serial_dependence_r2;
    analyzed_data.ind(subj).serial_dependence_estimated_bias = serial_dependence_estimated_bias;
end

%% Group


%% Super

% Calculate mean pCW and SEM for the super subject from pooled raw responses
pCW_per_offset_delta_theta_bin_mean = nan(num.levels, num.levels, num.conds, length(unique_probe_offsets), num.delta_theta_bins);
pCW_per_offset_delta_theta_bin_sem = nan(num.levels, num.levels, num.conds, length(unique_probe_offsets), num.delta_theta_bins);
resp_per_offset_per_delta_theta_bin_mean = nan(num.levels, num.levels, num.conds, length(unique_probe_offsets), num.delta_theta_bins);
perf_per_delta_theta_bin_mean = nan(num.levels, num.levels, num.conds, num.delta_theta_bins);
perf_per_delta_theta_bin_sem = nan(num.levels, num.levels, num.conds, num.delta_theta_bins);
pCW_per_delta_theta_bin_mean = nan(num.levels, num.levels, num.conds, num.delta_theta_bins);
pCW_per_delta_theta_bin_sem = nan(num.levels, num.levels, num.conds, num.delta_theta_bins);

for cond = 1:num.conds
    for delta_theta_bin_indx = 1:num.delta_theta_bins
        for prev_lvl = 1:num.levels
            for curr_lvl = 1:num.levels

                for offset_indx = 1:length(unique_probe_offsets)

                    pooled_responses = analyzed_data.super.resp_per_offset_per_delta_theta_bin{prev_lvl, curr_lvl, cond, offset_indx, delta_theta_bin_indx};
                    
                    if ~isempty(pooled_responses)
                        pCW_per_offset_delta_theta_bin_mean(prev_lvl, curr_lvl, cond, offset_indx, delta_theta_bin_indx) = mean(pooled_responses, 'omitnan');
                        resp_per_offset_per_delta_theta_bin_mean(prev_lvl, curr_lvl, cond, offset_indx, delta_theta_bin_indx) = mean(pooled_responses, 'omitnan');
                        if length(pooled_responses) >= 2
                            pCW_per_offset_delta_theta_bin_sem(prev_lvl, curr_lvl, cond, offset_indx, delta_theta_bin_indx) = std(pooled_responses, 'omitnan') / sqrt(length(pooled_responses));
                        else
                            pCW_per_offset_delta_theta_bin_sem(prev_lvl, curr_lvl, cond, offset_indx, delta_theta_bin_indx) = NaN;
                        end
                    else
                        pCW_per_offset_delta_theta_bin_mean(prev_lvl, curr_lvl, cond, offset_indx, delta_theta_bin_indx) = NaN;
                        pCW_per_offset_delta_theta_bin_sem(prev_lvl, curr_lvl, cond, offset_indx, delta_theta_bin_indx) = NaN;
                        resp_per_offset_per_delta_theta_bin_mean(prev_lvl, curr_lvl, cond, offset_indx, delta_theta_bin_indx) = NaN;
                    end
                end

                curr_responses = analyzed_data.super.perf_per_delta_theta_bin{prev_lvl, curr_lvl, cond, delta_theta_bin_indx};
                perf_per_delta_theta_bin_mean(prev_lvl, curr_lvl, cond, delta_theta_bin_indx) = mean(curr_responses, 'omitnan');
                perf_per_delta_theta_bin_sem(prev_lvl, curr_lvl, cond, delta_theta_bin_indx) = std(curr_responses, 'omitnan') ./ sqrt(sum(~isnan(curr_responses)));

                curr_pCW = analyzed_data.super.pCW_per_delta_theta_bin{prev_lvl, curr_lvl, cond, delta_theta_bin_indx};
                pCW_per_delta_theta_bin_mean(prev_lvl, curr_lvl, cond, delta_theta_bin_indx) = mean(curr_pCW, 'omitnan');
                pCW_per_delta_theta_bin_sem(prev_lvl, curr_lvl, cond, delta_theta_bin_indx) = std(curr_pCW, 'omitnan') ./ sqrt(sum(~isnan(curr_pCW)));

            end
        end
    end
end

% Store the results in the analyzed_data structure
analyzed_data.super.pCW_per_offset_delta_theta_bin_mean = pCW_per_offset_delta_theta_bin_mean;
analyzed_data.super.pCW_per_offset_delta_theta_bin_sem = pCW_per_offset_delta_theta_bin_sem;
analyzed_data.super.resp_per_offset_per_delta_theta_bin_mean = resp_per_offset_per_delta_theta_bin_mean;
analyzed_data.super.perf_per_delta_theta_bin_mean = perf_per_delta_theta_bin_mean;
analyzed_data.super.perf_per_delta_theta_bin_sem = perf_per_delta_theta_bin_sem;
analyzed_data.super.pCW_per_delta_theta_bin_mean = pCW_per_delta_theta_bin_mean;
analyzed_data.super.pCW_per_delta_theta_bin_sem = pCW_per_delta_theta_bin_sem;

% Initialize response bias arrays for super subject
mu_per_delta_theta_bin = nan(num.levels, num.levels, num.conds, num.delta_theta_bins);
sigma_per_delta_theta_bin = nan(num.levels, num.levels, num.conds, num.delta_theta_bins);
r2_per_delta_theta_bin = nan(num.levels, num.levels, num.conds, num.delta_theta_bins);

serial_dependence_params_est = nan(num.levels, num.levels, num.conds, 2);
serial_dependence_r2 = nan(num.levels, num.levels, num.conds);
serial_dependence_estimated_bias = nan(num.levels, num.levels, num.conds, length(-180:180));

for cond = 1:num.conds
    % disp(['Super Subject: ' cond_names{cond} ' Condition']);
    
    for prev_lvl = 1:num.levels
        % disp(['Super Subject: Previous LVL ' num2str(prev_lvl)]);
        
        for curr_lvl = 1:num.levels
            % disp(['Super Subject: Current LVL ' num2str(curr_lvl)]);
            
            % Fit response bias model for each delta theta bin
            for delta_theta_bin_indx = 1:num.delta_theta_bins % 1:round(num.delta_theta_bins/3):num.delta_theta_bins

                curr_probe_offsets = unique_probe_offsets;
                
                % Extract pCW for the current delta theta bin across all probe offsets for super subject
                p_CW = squeeze(pCW_per_offset_delta_theta_bin_mean(prev_lvl, curr_lvl, cond, :, delta_theta_bin_indx));
                
                rmv_rows = isnan(p_CW);
                curr_probe_offsets(rmv_rows) = [];
                p_CW(rmv_rows) = [];

                % Skip fitting if there are NaN values (insufficient data for this bin)
                if isempty(p_CW)
                    if toggles.disp_on, disp(['Super Subject: Not enough data for delta theta bin ' num2str(delta_theta_bin_edges(1, delta_theta_bin_indx)) ' to ' num2str(delta_theta_bin_edges(2, delta_theta_bin_indx)) '°']); end
                    continue;
                end

                % Estimate response bias
                try
                    fixed_params = cell(1,2);
                    fixed_params{1} = [curr_probe_offsets, p_CW]; % Offsets and corresponding pCW
                    fixed_params{2} = p.guess_rate;
                    
                    response_bias = estimateResponseBias(p.response_bias_init_params, fixed_params, p);
                    
                    mu_per_delta_theta_bin(prev_lvl, curr_lvl, cond, delta_theta_bin_indx) = response_bias.params_est(1);
                    sigma_per_delta_theta_bin(prev_lvl, curr_lvl, cond, delta_theta_bin_indx) = response_bias.params_est(2);
                    r2_per_delta_theta_bin(prev_lvl, curr_lvl, cond, delta_theta_bin_indx) = response_bias.r2;
                catch
                    if toggles.disp_on
                        disp(['Failed to estimate response bias for Super Subject, Prev Level: ' num2str(prev_lvl) ', Curr Level: ' num2str(curr_lvl) ', Condition: ' cond_names{cond} ', Delta Theta Bin: ' num2str(delta_theta_bin_indx)]);
                    end
                end
                
                clear fixed_params
            end

            % Fit serial dependence model
            curr_delta_theta = delta_theta_bin_centers';
            curr_mu = squeeze(mu_per_delta_theta_bin(prev_lvl, curr_lvl, cond, :));

            rmv_rows = isnan(curr_mu);
            curr_delta_theta(rmv_rows) = [];
            curr_mu(rmv_rows) = [];

            % Only skip if absolutely no data is available
            if isempty(curr_mu) 
                if toggles.disp_on, disp('Super Subject: No data for serial dependence model'); end
                continue;
            end

            % Estimate serial dependence - passing data directly, not in a cell
            fixed_params = [curr_delta_theta, curr_mu];
            serial_dependence = estimateSerialDependence(p.serial_dependence_init_params, fixed_params, p);

            serial_dependence_params_est(prev_lvl, curr_lvl, cond,:) = serial_dependence.params_est;
            serial_dependence_r2(prev_lvl, curr_lvl, cond) = serial_dependence.r2;
            serial_dependence_estimated_bias(prev_lvl, curr_lvl, cond,:) = gaussianPrime(serial_dependence.params_est, -180:180);
        end
    end
end

% Store the results in the analyzed_data structure
analyzed_data.super.mu_per_delta_theta_bin = mu_per_delta_theta_bin;
analyzed_data.super.sigma_per_delta_theta_bin = sigma_per_delta_theta_bin;
analyzed_data.super.r2_per_delta_theta_bin = r2_per_delta_theta_bin;
analyzed_data.super.serial_dependence_params_est = serial_dependence_params_est;
analyzed_data.super.serial_dependence_r2 = serial_dependence_r2;
analyzed_data.super.serial_dependence_estimated_bias = serial_dependence_estimated_bias;

%% Plot performance as a function of delta theta for each level pair

if plot_spr 
    % Set figure size for performance plots
    if exist('fg', 'var')
        set(fg, 'Position', [100 100 900 700]);
    end
    plotPerformancePerDeltaTheta(delta_theta_bin_centers, analyzed_data.super, 'super', num, cond_names, spr_figure_path, save_spr_figures, fg);
end

%% Plot pCW as a function of delta theta for each level pair     

if plot_spr
    % Set figure size for pCCW plots
    if exist('fg', 'var')
        set(fg, 'Position', [100 100 900 700]);
    end
    plotPCCWPerDeltaTheta(delta_theta_bin_centers, analyzed_data.super, 'super', num, cond_names, spr_figure_path, save_spr_figures, fg);
end

%% Plot response bias as a function of delta theta for each level pair
if plot_spr
    % Set figure size for response bias plots
    if exist('fg', 'var')
        set(fg, 'Position', [100 100 900 700]);
    end
    plotResponseBiasPerDeltaTheta(delta_theta_bin_centers, analyzed_data.super, 'super', num, cond_names, p, green, spr_figure_path, save_spr_figures, fg);
end

%% Plot Behavioral variance as a function of delta theta for each level pair
if plot_spr
    % Set figure size for behavioral variance plots
    if exist('fg', 'var')
        set(fg, 'Position', [100 100 900 700]);
    end
    plotBehavioralVariancePerDeltaTheta(delta_theta_bin_centers, analyzed_data.super, 'super', num, cond_names, p, spr_figure_path, save_spr_figures, fg);
end

%% Plot amplitude of serial dependence for each level pair
if plot_spr
    % Set figure size for serial dependence amplitude plots
    if exist('fg', 'var')
        set(fg, 'Position', [100 100 1000 500]);
    end
    plotSerialDependenceAmplitude(analyzed_data.super, 'super', num, cond_names, p, spr_figure_path, save_spr_figures, fg);
end

%% Plot width of serial dependence for each level pair
if plot_spr
    % Set figure size for serial dependence width plots
    if exist('fg', 'var')
        set(fg, 'Position', [100 100 1000 500]);
    end
    plotSerialDependenceWidth(analyzed_data.super, 'super', num, cond_names, spr_figure_path, save_spr_figures, fg);
end

%% Plot Super Subject Delta Theta Counts
if plot_spr
    % Set figure size for delta theta count plots
    if exist('fg', 'var')
        set(fg, 'Position', [100 100 1200 800]);
    end
    plotSuperDeltaThetaCount(num, analyzed_data, cond_names, spr_figure_path, save_spr_figures, fg);
end

%% Plot Individual Subject Data
if plot_ind
    for subj = 1:num.subjs
        figure_path = [ind_figure_path '/S' subj_IDs{subj}];
        if exist(figure_path, 'dir') == 0, mkdir(figure_path); end
        
        % Plot performance per delta theta
        if ~isempty(analyzed_data.ind(subj).pCW_per_offset_delta_theta_mean)
            % Set figure size for performance plots
            if exist('fg', 'var')
                set(fg, 'Position', [100 100 900 700]);
            end
            plotPerformancePerDeltaTheta(delta_theta_bin_centers, analyzed_data.ind(subj), ['S' subj_IDs{subj}], num, cond_names, figure_path, save_ind_figures, fg);
        end
        
        % Plot pCCW per delta theta
        if ~isempty(analyzed_data.ind(subj).pCW_per_offset_delta_theta_mean)
            % Set figure size for pCCW plots
            if exist('fg', 'var')
                set(fg, 'Position', [100 100 900 700]);
            end
            plotPCCWPerDeltaTheta(delta_theta_bin_centers, analyzed_data.ind(subj), ['S' subj_IDs{subj}], num, cond_names, figure_path, save_ind_figures, fg);
        end
        
        % Plot response bias per delta theta
        if any(~isnan(analyzed_data.ind(subj).mu_per_delta_theta(:)))
            % Set figure size for response bias plots
            if exist('fg', 'var')
                set(fg, 'Position', [100 100 900 700]);
            end
            plotResponseBiasPerDeltaTheta(delta_theta_bin_centers, analyzed_data.ind(subj), ['S' subj_IDs{subj}], num, cond_names, p, green, figure_path, save_ind_figures, fg);
        end
        
        % Plot behavioral variance per delta theta
        if any(~isnan(analyzed_data.ind(subj).sigma_per_delta_theta(:)))
            % Set figure size for behavioral variance plots
            if exist('fg', 'var')
                set(fg, 'Position', [100 100 900 700]);
            end
            plotBehavioralVariancePerDeltaTheta(delta_theta_bin_centers, analyzed_data.ind(subj), ['S' subj_IDs{subj}], num, cond_names, p, figure_path, save_ind_figures, fg);
        end
        
        % Plot serial dependence amplitude and width if available and not all NaN
        if any(~isnan(analyzed_data.ind(subj).serial_dependence_params_est(:)))
            % Set figure size for serial dependence amplitude plots
            if exist('fg', 'var')
                set(fg, 'Position', [100 100 1000 500]);
            end
            plotSerialDependenceAmplitude(analyzed_data.ind(subj), ['S' subj_IDs{subj}], num, cond_names, p, figure_path, save_ind_figures, fg);
            
            % Set figure size for serial dependence width plots
            if exist('fg', 'var')
                set(fg, 'Position', [100 100 1000 500]);
            end
            plotSerialDependenceWidth(analyzed_data.ind(subj), ['S' subj_IDs{subj}], num, cond_names, figure_path, save_ind_figures, fg);
        end
        
        % Plot delta theta counts for individual subjects
        % Set figure size for delta theta count plots
        if exist('fg', 'var')
            set(fg, 'Position', [100 100 1200 800]);
        end
        plotSubjectDeltaThetaCount(num, analyzed_data, subj, subj_IDs{subj}, cond_names, figure_path, save_ind_figures, fg);
    end
end

%%

if plot_spr || plot_ind || plot_grp && exist('fg', 'var'), close(fg); end