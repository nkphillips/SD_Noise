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

%% Define response bias model parameters

p.response_bias_bounds = [20, 10; -20, 0]; % upper, lower bounds for mu and sigma
p.guess_rate = 0.25; % Assuming guess rate is constant
p.fmincon_options = optimoptions('fmincon','Display','off');

%% Initialize delta theta windows

delta_theta_window_width = 32;
delta_theta_window_centers = (-180+delta_theta_window_width/2):(180-delta_theta_window_width/2);
num.delta_theta_windows = length(delta_theta_window_centers);

delta_theta_window_edges = [delta_theta_window_centers - delta_theta_window_width/2; delta_theta_window_centers + delta_theta_window_width/2];

analyzed_data.super.responses_by_offset_delta_theta_window = cell(num.levels, num.levels, num.conds, length(unique_probe_offsets), num.delta_theta_windows);
analyzed_data.super.perf_by_delta_theta_window = cell(num.levels, num.levels, num.conds, num.delta_theta_windows);
analyzed_data.super.pCW_by_delta_theta_window = cell(num.levels, num.levels, num.conds, num.delta_theta_windows);
analyzed_data.super.mu_by_delta_theta_window = cell(num.levels, num.levels, num.conds, num.delta_theta_windows);
analyzed_data.super.sigma_by_delta_theta_window = cell(num.levels, num.levels, num.conds, num.delta_theta_windows);
analyzed_data.super.r2_by_delta_theta_window = cell(num.levels, num.levels, num.conds, num.delta_theta_windows);


%% 

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
    lvl_pair_perf_mean = nan(num.levels, num.levels, num.conds, num.runs(subj));
    lvl_pair_perf_sem = lvl_pair_perf_mean;
    lvl_pair_pCW_mean = nan(num.levels, num.levels, num.conds, num.runs(subj));
    lvl_pair_pCW_sem = lvl_pair_pCW_mean;

    % Initialize storage for pCW binned by offset and delta theta across runs, adding a dimension for runs
    lvl_pair_pCW_by_probe_offset = cell(num.levels, num.levels, num.conds, length(unique_probe_offsets));
    lvl_pair_pCW_by_probe_offset_mean = nan(num.levels, num.levels, num.conds, length(unique_probe_offsets));
    lvl_pair_pCW_by_probe_offset_sem = lvl_pair_pCW_by_probe_offset_mean;
    % Initialize storage for pCW binned by offset and delta theta across runs, adding a dimension for runs
    current_subj_pCW_binned_runs = nan(num.levels, num.levels, num.conds, length(unique_probe_offsets), num.delta_theta_windows, num.runs(subj)); 

    % Loop through runs
    for n_run = 1:num.runs(subj)

        % if toggles.disp_on, disp(['Run ' num2str(n_run)]); end

        subj_p = all_runs{subj}(n_run).p;
        behav_data = all_runs{subj}(n_run).behav_data;

        counts(n_run) = getConditionDistributions(subj_p, behav_data, unique_probe_offsets, ind_figure_path);

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

                for prev_lvl = 1:num.levels
                    for curr_lvl = 1:num.levels

                       % if toggles.disp_on, disp(['Level Pair ' num2str(prev_lvl) '-' num2str(curr_lvl)]); end

                        curr_lvl_pair_indices = counts(n_run).lvl_pair_indices{prev_lvl, curr_lvl, n_block, cond};

                        curr_lvl_pair_correct = curr_correct(curr_lvl_pair_indices);
                        lvl_pair_perf_mean(prev_lvl, curr_lvl, cond, n_run) = mean(curr_lvl_pair_correct);
                        lvl_pair_perf_sem(prev_lvl, curr_lvl, cond, n_run) = std(curr_lvl_pair_correct)/sqrt(length(curr_lvl_pair_correct));

                        curr_lvl_pair_response_CW = curr_response(curr_lvl_pair_indices);
                        lvl_pair_pCW_mean(prev_lvl, curr_lvl, cond, n_run) = mean(curr_lvl_pair_response_CW);
                        lvl_pair_pCW_sem(prev_lvl, curr_lvl, cond, n_run) = std(curr_lvl_pair_response_CW)/sqrt(length(curr_lvl_pair_response_CW));

                        curr_lvl_pair_delta_theta = curr_delta_theta(curr_lvl_pair_indices-1);
                        curr_lvl_pair_probe_offset = curr_probe_offsets(curr_lvl_pair_indices);
                        
                        % Check if any data exists for this level pair/condition
                        if isempty(curr_lvl_pair_delta_theta)
                            continue; % Skip if no trials
                        end
    
                        for delta_theta_window_indx = 1:num.delta_theta_windows

                            % if toggles.disp_on, disp(['Delta Theta Window ' num2str(delta_theta_window_indx)]); end

                            curr_delta_theta_window = [delta_theta_window_edges(1, delta_theta_window_indx), delta_theta_window_edges(2, delta_theta_window_indx)];
    
                            % Find trials within the current delta theta window
                            delta_theta_mask = (curr_lvl_pair_delta_theta >= curr_delta_theta_window(1)) & (curr_lvl_pair_delta_theta <= curr_delta_theta_window(2));
    
                            analyzed_data.super.perf_by_delta_theta_window{prev_lvl, curr_lvl, cond, delta_theta_window_indx} = [analyzed_data.super.perf_by_delta_theta_window{prev_lvl, curr_lvl, cond, delta_theta_window_indx}; ...
                                curr_lvl_pair_correct(delta_theta_mask)];
                            analyzed_data.super.pCW_by_delta_theta_window{prev_lvl, curr_lvl, cond, delta_theta_window_indx} = [analyzed_data.super.pCW_by_delta_theta_window{prev_lvl, curr_lvl, cond, delta_theta_window_indx}; ...
                                curr_lvl_pair_response_CW(delta_theta_mask)];
                            
                            for offset_indx = 1:length(unique_probe_offsets)

                                % if toggles.disp_on, disp(['Offset ' num2str(offset_indx)]); end

                                curr_offset = unique_probe_offsets(offset_indx);
    
                                % Find trials matching the current probe offset *within* the delta theta window
                                offset_mask = (curr_lvl_pair_probe_offset == curr_offset);
                                combined_mask = delta_theta_mask & offset_mask;
    
                                % Get responses for the selected trials
                                responses_in_bin = curr_lvl_pair_response_CW(combined_mask);

                                analyzed_data.super.responses_by_offset_delta_theta_window{prev_lvl, curr_lvl, cond, offset_indx, delta_theta_window_indx} = [analyzed_data.super.responses_by_offset_delta_theta_window{prev_lvl, curr_lvl, cond, offset_indx, delta_theta_window_indx}; ...
                                    responses_in_bin];
    
                                % Calculate pCW for this specific bin
                                if ~isempty(responses_in_bin)
                                    pCW_in_bin = mean(responses_in_bin, 'omitnan');
                                else
                                    pCW_in_bin = NaN; % Assign NaN if no trials in this bin
                                end
    
                                % Store the calculated pCW for this run in the numeric array
                                current_subj_pCW_binned_runs(prev_lvl, curr_lvl, cond, offset_indx, delta_theta_window_indx, n_run) = pCW_in_bin;
    
                                if delta_theta_window_indx == 1 % This seems unrelated to the delta_theta binning logic, might need review? Assumed correct for now.
                                lvl_pair_pCW_by_probe_offset{prev_lvl, curr_lvl,cond,offset_indx} = [lvl_pair_pCW_by_probe_offset{prev_lvl, curr_lvl,cond,offset_indx}, ...
                                    counts(n_run).pCW_by_probe_offset{prev_lvl, curr_lvl,cond, offset_indx}];
                                end
    
                            end
                        end

                    end
                end
            
            end

        end

    end

    analyzed_data.ind(subj).counts = counts;

    analyzed_data.ind(subj).performance_within_run = lvl_pair_perf_mean;
    analyzed_data.ind(subj).performance_within_run_sem = lvl_pair_perf_sem;
    analyzed_data.ind(subj).performance_across_runs = mean(analyzed_data.ind(subj).performance_within_run, 4);
    analyzed_data.ind(subj).performance_across_runs_sem = std(analyzed_data.ind(subj).performance_within_run, [], 4)/sqrt(num.runs(subj));

    analyzed_data.ind(subj).pCW_within_run = lvl_pair_pCW_mean;
    analyzed_data.ind(subj).pCW_within_run_sem = lvl_pair_pCW_sem;
    analyzed_data.ind(subj).pCW_across_runs = mean(analyzed_data.ind(subj).pCW_within_run, 4);
    analyzed_data.ind(subj).pCW_across_runs_sem = std(analyzed_data.ind(subj).pCW_within_run, [], 4)/sqrt(num.runs(subj));

    analyzed_data.ind(subj).pCW_by_offset_delta_theta_all_runs = current_subj_pCW_binned_runs; % Store raw per-run data
    analyzed_data.ind(subj).pCW_by_probe_offset = lvl_pair_pCW_by_probe_offset;

    % --- Calculate mean and SEM for pCW binned by offset and delta theta across runs ---
    tmp_pCW_binned_across_runs = current_subj_pCW_binned_runs; % Use the numeric array containing data from all runs

    % Calculate mean across runs (dimension 6), ignoring NaNs
    analyzed_data.ind(subj).pCW_by_offset_delta_theta_mean = mean(tmp_pCW_binned_across_runs, 6, 'omitnan');

    % Calculate the number of non-NaN runs for each bin
    n_runs_with_data = sum(~isnan(tmp_pCW_binned_across_runs), 6); 

    % Calculate SEM across runs (dimension 6), ignoring NaNs
    analyzed_data.ind(subj).pCW_by_offset_delta_theta_sem = std(tmp_pCW_binned_across_runs, 0, 6, 'omitnan') ./ sqrt(n_runs_with_data);
    analyzed_data.ind(subj).pCW_by_offset_delta_theta_sem(n_runs_with_data < 2) = NaN; % Set SEM to NaN if less than 2 data points
    % ------------------------------------------------------------------------------------

    %% Fit model to pCW by probe offset for each delta theta bin
    
    % Initialize storage for fitted parameters per delta theta bin
    analyzed_data.ind(subj).mu_by_delta_theta = nan(num.levels, num.levels, num.conds, num.delta_theta_windows);
    analyzed_data.ind(subj).sigma_by_delta_theta = nan(num.levels, num.levels, num.conds, num.delta_theta_windows);
    analyzed_data.ind(subj).r2_by_delta_theta = nan(num.levels, num.levels, num.conds, num.delta_theta_windows);

    for prev_lvl = 1:num.levels
        disp(['Previous LVL ' num2str(prev_lvl)]);
        for curr_lvl = 1:num.levels
            disp(['Current LVL ' num2str(curr_lvl)]);
            for cond = 1:num.conds
                disp([cond_names{cond} ' Condition']);
                
                for offset_indx = 1:length(unique_probe_offsets)
                    analyzed_data.ind(subj).pCW_by_probe_offset_mean(prev_lvl, curr_lvl, cond, offset_indx) = mean(analyzed_data.ind(subj).pCW_by_probe_offset{prev_lvl, curr_lvl, cond, offset_indx},'omitnan');
                    analyzed_data.ind(subj).pCW_by_probe_offset_sem(prev_lvl, curr_lvl, cond, offset_indx) = std(analyzed_data.ind(subj).pCW_by_probe_offset{prev_lvl, curr_lvl, cond, offset_indx},'omitnan')/sqrt(num.runs(subj));
                end

                %% Fit model to pCW by probe offset for each delta theta bin

                for delta_theta_window_indx = 1:num.delta_theta_windows
                    
                    curr_probe_offsets = unique_probe_offsets;

                    % Extract pCW for the current delta theta window across all probe offsets
                    p_CW = squeeze(analyzed_data.ind(subj).pCW_by_offset_delta_theta_mean(prev_lvl, curr_lvl, cond, :, delta_theta_window_indx));

                    rmv_rows = isnan(p_CW);
                    curr_probe_offsets(rmv_rows) = [];
                    p_CW(rmv_rows) = [];

                    % Skip fitting if there are NaN values (insufficient data for this bin)
                    if isempty(p_CW) % sum(~isnan(p_CW)) < 2 || any(isnan(p_CW))
                        disp(['Not enough data for delta theta bin ' num2str(delta_theta_window_edges(1, delta_theta_window_indx)) ' to ' num2str(delta_theta_window_edges(2, delta_theta_window_indx)) '°']);
                        continue;
                    end

                    fixed_params{1} = [curr_probe_offsets, p_CW]; % Offsets and corresponding pCW
                    fixed_params{2} = p.guess_rate;

                    response_bias = estimateResponseBias(fixed_params, p, toggles);

                    mu = response_bias.params_est(1);
                    sigma = response_bias.params_est(2);
                    r2 = response_bias.r2;

                    analyzed_data.ind(subj).mu_by_delta_theta(prev_lvl, curr_lvl, cond, delta_theta_window_indx) = mu;
                    analyzed_data.ind(subj).sigma_by_delta_theta(prev_lvl, curr_lvl, cond, delta_theta_window_indx) = sigma;
                    analyzed_data.ind(subj).r2_by_delta_theta(prev_lvl, curr_lvl, cond, delta_theta_window_indx) = r2;

                    %% Plot response bias 
                    %{

                    if plot_ind
                        x = -20:0.01:20;
                        cdf_response = calc_pCW(x, mu, sigma, p.guess_rate);
                        pdf_response = (1 - p.guess_rate) * normpdf(x, mu, sigma);
                        pdf_response = pdf_response / max(pdf_response); % normalize to peak at 1

                        figure_name = ['S' subj_IDs{subj} '_response_bias_L' num2str(prev_lvl) '-L' num2str(curr_lvl) '_' cond_names{cond} '_DeltaThetaBin_' num2str(delta_theta_window_edges(1, delta_theta_window_indx)) '–' num2str(delta_theta_window_edges(2, delta_theta_window_indx)) '°'];
                    
                        plot(x, cdf_response', 'LineWidth', 1,'Color', 'b');
                        hold on;
                        plot(x, pdf_response', 'LineWidth', 1,'LineStyle', '--','Color', 'b');
                        scatter(curr_probe_offsets(curr_probe_offsets <= 0), p_CW(curr_probe_offsets <= 0), 30, 'filled','MarkerFaceColor', 'r');
                        scatter(curr_probe_offsets(curr_probe_offsets > 0), p_CW(curr_probe_offsets > 0), 30, 'filled','MarkerFaceColor', 'g');

                        line([0,0], [0,1], 'LineWidth', 1, 'Color', 'k');
                        line([-20, 20], [0.5, 0.5], 'LineWidth', 1, 'Color', 'k');
                        title(['Delta Theta Bin ' num2str(delta_theta_window_indx) ': \mu = ' num2str(round(mu, 2)) ', \sigma = ' num2str(round(sigma, 2)) ' (R^2 = ' num2str(round(r2, 2)) ')']);
                        xlabel('Probe Offset');
                        ylabel('P(Resp|CW)');
                        ylim([0, 1]);
                        xlim([-20, 20]);
                        box off
                        set(gca, 'TickDir', 'out');

                        if save_ind_figures
                            saveas(gcf, [figure_path '/' figure_name '.png']); clf;
                            disp(['Saved ' figure_name '.png']);
                        end

                    end

                    %}
                end 

            end 
        end 
    end 

end 

if plot_ind, close(fg); end

%% Group


%% Super

% Calculate mean pCW and SEM for the super subject from pooled raw responses
analyzed_data.super.pCW_by_offset_delta_theta_window_mean = nan(num.levels, num.levels, num.conds, length(unique_probe_offsets), num.delta_theta_windows);
analyzed_data.super.pCW_by_offset_delta_theta_window_sem = nan(num.levels, num.levels, num.conds, length(unique_probe_offsets), num.delta_theta_windows);
analyzed_data.super.responses_by_offset_delta_theta_window_mean = nan(num.levels, num.levels, num.conds, length(unique_probe_offsets), num.delta_theta_windows);

for cond = 1:num.conds
    for delta_theta_window_indx = 1:num.delta_theta_windows
        for prev_lvl = 1:num.levels
            for curr_lvl = 1:num.levels

                for offset_indx = 1:length(unique_probe_offsets)

                    pooled_responses = analyzed_data.super.responses_by_offset_delta_theta_window{prev_lvl, curr_lvl, cond, offset_indx, delta_theta_window_indx};
                    
                    if ~isempty(pooled_responses)
                        analyzed_data.super.pCW_by_offset_delta_theta_window_mean(prev_lvl, curr_lvl, cond, offset_indx, delta_theta_window_indx) = mean(pooled_responses, 'omitnan');
                        analyzed_data.super.responses_by_offset_delta_theta_window_mean(prev_lvl, curr_lvl, cond, offset_indx, delta_theta_window_indx) = mean(pooled_responses, 'omitnan');
                        if length(pooled_responses) >= 2
                            analyzed_data.super.pCW_by_offset_delta_theta_window_sem(prev_lvl, curr_lvl, cond, offset_indx, delta_theta_window_indx) = std(pooled_responses, 'omitnan') / sqrt(length(pooled_responses));
                        else
                            analyzed_data.super.pCW_by_offset_delta_theta_sem(prev_lvl, curr_lvl, cond, offset_indx, delta_theta_window_indx) = NaN;
                        end
                    else
                        analyzed_data.super.pCW_by_offset_delta_theta_window_mean(prev_lvl, curr_lvl, cond, offset_indx, delta_theta_window_indx) = NaN;
                        analyzed_data.super.pCW_by_offset_delta_theta_window_sem(prev_lvl, curr_lvl, cond, offset_indx, delta_theta_window_indx) = NaN;
                        analyzed_data.super.responses_by_offset_delta_theta_window_mean(prev_lvl, curr_lvl, cond, offset_indx, delta_theta_window_indx) = NaN;
                    end
                end

                curr_responses = analyzed_data.super.perf_by_delta_theta_window{prev_lvl, curr_lvl, cond, delta_theta_window_indx};
                analyzed_data.super.perf_by_delta_theta_window_mean(prev_lvl, curr_lvl, cond, delta_theta_window_indx) = mean(curr_responses, 'omitnan');
                analyzed_data.super.perf_by_delta_theta_window_sem(prev_lvl, curr_lvl, cond, delta_theta_window_indx) = std(curr_responses, 'omitnan') ./ sqrt(sum(~isnan(curr_responses)));

                curr_pCW = analyzed_data.super.pCW_by_delta_theta_window{prev_lvl, curr_lvl, cond, delta_theta_window_indx};
                analyzed_data.super.pCW_by_delta_theta_window_mean(prev_lvl, curr_lvl, cond, delta_theta_window_indx) = mean(curr_pCW, 'omitnan');
                analyzed_data.super.pCW_by_delta_theta_window_sem(prev_lvl, curr_lvl, cond, delta_theta_window_indx) = std(curr_pCW, 'omitnan') ./ sqrt(sum(~isnan(curr_pCW)));

            end
        end
    end
end

% Fit model to pCW by probe offset for each delta theta bin for the super subject
analyzed_data.super.mu_by_delta_theta_window = nan(num.levels, num.levels, num.conds, num.delta_theta_windows);
analyzed_data.super.sigma_by_delta_theta_window = nan(num.levels, num.levels, num.conds, num.delta_theta_windows);
analyzed_data.super.r2_by_delta_theta_window = nan(num.levels, num.levels, num.conds, num.delta_theta_windows);

for cond = 1:num.conds
    % disp(['Super Subject: ' cond_names{cond} ' Condition']);
    
    for delta_theta_window_indx = 1:num.delta_theta_windows % 1:round(num.delta_theta_windows/3):num.delta_theta_windows
        
        for prev_lvl = 1:num.levels
            % disp(['Super Subject: Previous LVL ' num2str(prev_lvl)]);
            
            for curr_lvl = 1:num.levels
                % disp(['Super Subject: Current LVL ' num2str(curr_lvl)]);
        
                curr_probe_offsets = unique_probe_offsets;
                
                % Extract pCW for the current delta theta window across all probe offsets for super subject
                p_CW = squeeze(analyzed_data.super.pCW_by_offset_delta_theta_window_mean(prev_lvl, curr_lvl, cond, :, delta_theta_window_indx));
                
                rmv_rows = isnan(p_CW);
                curr_probe_offsets(rmv_rows) = [];
                p_CW(rmv_rows) = [];

                % Skip fitting if there are NaN values (insufficient data for this bin)
                if isempty(p_CW) % sum(~isnan(p_CW)) < 2 || any(isnan(p_CW))
                    disp(['Super Subject: Not enough data for delta theta bin ' num2str(delta_theta_window_edges(1, delta_theta_window_indx)) ' to ' num2str(delta_theta_window_edges(2, delta_theta_window_indx)) '°']);
                    continue;
                end

                fixed_params{1} = [curr_probe_offsets, p_CW]; % Offsets and corresponding pCW
                fixed_params{2} = p.guess_rate;
                
                response_bias = estimateResponseBias(fixed_params, p, toggles);
                
                mu = response_bias.params_est(1);
                sigma = response_bias.params_est(2);
                r2 = response_bias.r2;
                
                analyzed_data.super.mu_by_delta_theta_window(prev_lvl, curr_lvl, cond, delta_theta_window_indx) = mu;
                analyzed_data.super.sigma_by_delta_theta_window(prev_lvl, curr_lvl, cond, delta_theta_window_indx) = sigma;
                analyzed_data.super.r2_by_delta_theta_window(prev_lvl, curr_lvl, cond, delta_theta_window_indx) = r2;
                

                %% Plot response bias per delta theta window
                %{
                
                if plot_spr

                    x = -20:0.01:20;
                    cdf_response = calc_pCW(x, mu, sigma, p.guess_rate);
                    pdf_response = (1 - p.guess_rate) * normpdf(x, mu, sigma);
                    pdf_response = pdf_response / max(pdf_response); % normalize to peak at 1

                    figure_name = ['Super Response Bias ' cond_names{cond} ' L' num2str(prev_lvl) '-L' num2str(curr_lvl) ' DeltaTheta Bin ' num2str(delta_theta_window_edges(1, delta_theta_window_indx)) '–' num2str(delta_theta_window_edges(2, delta_theta_window_indx)) '°'];
                
                    plot(x, cdf_response', 'LineWidth', 1,'Color', 'b');
                    hold on;
                    plot(x, pdf_response', 'LineWidth', 1,'LineStyle', '--','Color', 'b');
                    scatter(curr_probe_offsets(curr_probe_offsets <= 0), p_CW(curr_probe_offsets <= 0), 30, 'filled','MarkerFaceColor', 'r');
                    scatter(curr_probe_offsets(curr_probe_offsets > 0), p_CW(curr_probe_offsets > 0), 30, 'filled','MarkerFaceColor', 'g');

                    line([0,0], [0,1], 'LineWidth', 1, 'Color', 'k');
                    line([-20, 20], [0.5, 0.5], 'LineWidth', 1, 'Color', 'k');
                    title(['\mu = ' num2str(round(mu, 2)) ', \sigma = ' num2str(round(sigma, 2)) ', R^2 = ' num2str(round(r2, 2))]);
                    xlabel('Probe Offset');
                    ylabel('P(Resp|CW)');
                    ylim([0, 1]);
                    xlim([-20, 20]);
                    box off
                    set(gca, 'TickDir', 'out');
                
                    if save_spr_figures
                         saveas(gcf, [spr_figure_path '/' figure_name '.png']); clf;
                         disp(['Saved ' figure_name '.png']);
                    end

                end
                %}

            end

        end

    end

    %% Plot performance as a function of delta theta for each level pair

    %
    if plot_spr 
        for curr_lvl = 1:num.levels
            for prev_lvl = 1:num.levels
                
                figure_name = ['Super Performance ' cond_names{cond} ' by Delta Level'];

                subplot(num.levels, num.levels, prev_lvl + (curr_lvl-1)*num.levels);
                shadedErrorBar(delta_theta_window_centers, squeeze(analyzed_data.super.perf_by_delta_theta_window_mean(prev_lvl, curr_lvl, cond, :)), squeeze(analyzed_data.super.perf_by_delta_theta_window_sem(prev_lvl, curr_lvl, cond, :)),'lineprops', '-k');

                title([num2str(prev_lvl) ' -> ' num2str(curr_lvl)]);
                if curr_lvl == 3, xlabel('\Delta \theta (°)'); end
                if prev_lvl == 1, ylabel('Correct'); end
                ylim([0, 1]);
                xlim([-100 100]);
                xticks(min(xlim):50:max(xlim));
                line([min(xlim), max(xlim)], [0.5, 0.5], 'LineWidth', 1, 'Color', 'k');
                box off
                set(gca, 'TickDir', 'out');

            end
        end
    end

    if save_spr_figures
        saveas(gcf, [spr_figure_path '/' figure_name '.png']); clf;
        disp(['Saved ' figure_name '.png']);
    end
    %}

    %% Plot pCW as a function of delta theta for each level pair    

    %
    if plot_spr

        for prev_lvl = 1:num.levels
            for curr_lvl = 1:num.levels

                figure_name = ['Super pCW ' cond_names{cond} ' by Delta Level'];

                subplot(num.levels, num.levels, prev_lvl + (curr_lvl-1)*num.levels);

                shadedErrorBar(delta_theta_window_centers, squeeze(analyzed_data.super.pCW_by_delta_theta_window_mean(prev_lvl, curr_lvl, cond, :)), squeeze(analyzed_data.super.pCW_by_delta_theta_window_sem(prev_lvl, curr_lvl, cond, :)),'lineprops', '-b');
                hold on;
                
                title([num2str(prev_lvl) ' -> ' num2str(curr_lvl)]);
                if curr_lvl == 3, xlabel('\Delta \theta (°)'); end
                if prev_lvl == 1, ylabel('p(Resp|CW)'); end
                ylim([0, 1]);
                xlim([-100 100]);
                xticks(min(xlim):50:max(xlim));
                line([min(xlim), max(xlim)], [0.5, 0.5], 'LineWidth', 1, 'Color', 'k');

                box off
                set(gca, 'TickDir', 'out');

            end
        end

    end

    if save_spr_figures
        saveas(gcf, [spr_figure_path '/' figure_name '.png']); clf;
        disp(['Saved ' figure_name '.png']);
    end
    %}

    %% Plot response bias as a function of delta theta for each level pair
    %
    if plot_spr

        for prev_lvl = 1:num.levels
            for curr_lvl = 1:num.levels

                figure_name = ['Super Response Bias ' cond_names{cond} ' by Delta Level'];

                subplot(num.levels, num.levels, prev_lvl + (curr_lvl-1)*num.levels);

                plot(delta_theta_window_centers, squeeze(analyzed_data.super.mu_by_delta_theta_window(prev_lvl, curr_lvl, cond, :)), 'LineWidth', 1,'Color', 'g');

                title([num2str(prev_lvl) ' -> ' num2str(curr_lvl)]);
                if curr_lvl == 3, xlabel('\Delta \theta (°)'); end
                if prev_lvl == 1, ylabel('\mu'); end
                ylim([-20 20]);
                xlim([-100 100]);
                xticks(min(xlim):50:max(xlim));
                line([min(xlim), max(xlim)], [0, 0], 'LineWidth', 1, 'Color', 'k');

                box off
                set(gca, 'TickDir', 'out');

            end
        end

    end

    if save_spr_figures
        saveas(gcf, [spr_figure_path '/' figure_name '.png']); clf;
        disp(['Saved ' figure_name '.png']);
    end
    %}

end

if plot_spr, close(fg); end