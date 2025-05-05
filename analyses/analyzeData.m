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

%% Initialize delta theta windows

delta_theta_window_width = 32;
delta_theta_window_centers = (-180+delta_theta_window_width/2):delta_theta_window_width:(180-delta_theta_window_width/2);
num.delta_theta_windows = length(delta_theta_window_centers);

delta_theta_window_edges = [delta_theta_window_centers - delta_theta_window_width/2; delta_theta_window_centers + delta_theta_window_width/2];

%%

analyzed_data.super.responses_by_offset_delta_theta = cell(num.levels, num.levels, num.conds, length(unique_probe_offsets), num.delta_theta_windows);

%% Subject

if plot_ind 
    fg = figure('Visible','on','Color','w');
    set(0,'CurrentFigure',fg); 
end

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
    lvl_pair_pCW_by_offset_delta_theta_all_runs = nan(num.levels, num.levels, num.conds, length(unique_probe_offsets), num.delta_theta_windows, num.runs(subj)); 

    % Loop through runs
    for n_run = 1:num.runs(subj)

        if toggles.disp_on, disp(['Run ' num2str(n_run)]); end

        subj_p = all_runs{subj}(n_run).p;
        behav_data = all_runs{subj}(n_run).behav_data;

        counts(n_run) = getCondDist(subj_p, behav_data, unique_probe_offsets, ind_figure_path);

        delta_thetas = squeeze(subj_p.trial_events(1:end-1,1,:) - subj_p.trial_events(2:end,1,:));
        probe_offsets = squeeze(subj_p.trial_events(:,2,:) - subj_p.trial_events(:,1,:));

        for cond = 1:num.conds
            
            if toggles.disp_on, disp(['Condition ' num2str(cond)]); end
            curr_cond_blocks = find(subj_p.cond_order == cond);

            for n_block = 1:length(curr_cond_blocks)

                curr_correct = behav_data.correct(:,curr_cond_blocks(n_block));
                curr_response = behav_data.response(:,curr_cond_blocks(n_block));
                curr_delta_theta = delta_thetas(:,curr_cond_blocks(n_block));
                curr_probe_offsets = probe_offsets(:,curr_cond_blocks(n_block));

                for prev_lvl = 1:num.levels
                    for curr_lvl = 1:num.levels

                        if toggles.disp_on, disp(['Level Pair ' num2str(prev_lvl) '-' num2str(curr_lvl)]); end

                        curr_lvl_pair_indices = counts(n_run).lvl_pair_indices{prev_lvl, curr_lvl, n_block, cond};

                        curr_lvl_pair_correct = curr_correct(curr_lvl_pair_indices);
                        lvl_pair_perf_mean(prev_lvl, curr_lvl, cond, n_run) = mean(curr_lvl_pair_correct);
                        lvl_pair_perf_sem(prev_lvl, curr_lvl, cond, n_run) = std(curr_lvl_pair_correct)/sqrt(length(curr_lvl_pair_correct));

                        curr_lvl_pair_response_CW = curr_response(curr_lvl_pair_indices) == 2;
                        lvl_pair_pCW_mean(prev_lvl, curr_lvl, cond, n_run) = mean(curr_lvl_pair_response_CW);
                        lvl_pair_pCW_sem(prev_lvl, curr_lvl, cond, n_run) = std(curr_lvl_pair_response_CW)/sqrt(length(curr_lvl_pair_response_CW));

                        curr_lvl_pair_delta_theta = curr_delta_theta(curr_lvl_pair_indices-1);
                        curr_lvl_pair_probe_offset = curr_probe_offsets(curr_lvl_pair_indices);
                        
                        % Check if any data exists for this level pair/condition
                        if isempty(curr_lvl_pair_delta_theta)
                            continue; % Skip if no trials
                        end
    
                        for delta_theta_window_indx = 1:num.delta_theta_windows

                            if toggles.disp_on, disp(['Delta Theta Window ' num2str(delta_theta_window_indx)]); end

                            curr_delta_theta_window = [delta_theta_window_edges(1, delta_theta_window_indx), delta_theta_window_edges(2, delta_theta_window_indx)];
    
                            % Find trials within the current delta theta window
                            delta_theta_mask = (curr_lvl_pair_delta_theta >= curr_delta_theta_window(1)) & (curr_lvl_pair_delta_theta <= curr_delta_theta_window(2));
    
                            for offset_indx = 1:length(unique_probe_offsets)

                                if toggles.disp_on, disp(['Offset ' num2str(offset_indx)]); end

                                curr_offset = unique_probe_offsets(offset_indx);
    
                                % Find trials matching the current probe offset *within* the delta theta window
                                offset_mask = (curr_lvl_pair_probe_offset == curr_offset);
                                combined_mask = delta_theta_mask & offset_mask;
    
                                % Get responses for the selected trials
                                responses_in_bin = curr_lvl_pair_response_CW(combined_mask);

                                analyzed_data.super.responses_by_offset_delta_theta{prev_lvl, curr_lvl, cond, offset_indx, delta_theta_window_indx} = [analyzed_data.super.responses_by_offset_delta_theta{prev_lvl, curr_lvl, cond, offset_indx, delta_theta_window_indx}; ...
                                    responses_in_bin];
    
                                % Calculate pCW for this specific bin (mean of response == 2)
                                if ~isempty(responses_in_bin)
                                    pCW_in_bin = mean(responses_in_bin, 'omitnan');
                                else
                                    pCW_in_bin = NaN; % Assign NaN if no trials in this bin
                                end
    
                                % Store the calculated pCW for this run in the numeric array
                                lvl_pair_pCW_by_offset_delta_theta_all_runs(prev_lvl, curr_lvl, cond, offset_indx, delta_theta_window_indx, n_run) = pCW_in_bin;
    
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

    analyzed_data.ind(subj).pCW_by_offset_delta_theta_all_runs = lvl_pair_pCW_by_offset_delta_theta_all_runs; % Store raw per-run data
    analyzed_data.ind(subj).pCW_by_probe_offset = lvl_pair_pCW_by_probe_offset;

    % --- Calculate mean and SEM for pCW binned by offset and delta theta across runs ---
    tmp_pCW_binned_across_runs = lvl_pair_pCW_by_offset_delta_theta_all_runs; % Use the numeric array containing data from all runs

    % Calculate mean across runs (dimension 6), ignoring NaNs
    analyzed_data.ind(subj).pCW_by_offset_delta_theta_mean = mean(tmp_pCW_binned_across_runs, 6, 'omitnan');

    % Calculate the number of non-NaN runs for each bin
    n_runs_with_data = sum(~isnan(tmp_pCW_binned_across_runs), 6); 

    % Calculate SEM across runs (dimension 6), ignoring NaNs
    analyzed_data.ind(subj).pCW_by_offset_delta_theta_sem = std(tmp_pCW_binned_across_runs, 0, 6, 'omitnan') ./ sqrt(n_runs_with_data);
    analyzed_data.ind(subj).pCW_by_offset_delta_theta_sem(n_runs_with_data < 2) = NaN; % Set SEM to NaN if less than 2 data points
    % ------------------------------------------------------------------------------------

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

                guess_rate = 0.25; % Assuming guess rate is constant

                for delta_theta_window_indx = 1:num.delta_theta_windows

                    % Extract pCW for the current delta theta window across all probe offsets
                    p_CW = squeeze(analyzed_data.ind(subj).pCW_by_offset_delta_theta_mean(prev_lvl, curr_lvl, cond, :, delta_theta_window_indx))';

                    % Skip fitting if there are NaN values (insufficient data for this bin)
                    if any(isnan(p_CW)) || isempty(p_CW)
                        disp(['Not enough data for delta theta bin ' num2str(delta_theta_window_edges(1, delta_theta_window_indx)) ' to ' num2str(delta_theta_window_edges(2, delta_theta_window_indx)) '°']);
                        continue;
                    end

                    fixed_params{1} = [unique_probe_offsets', p_CW']; % Offsets and corresponding pCW
                    fixed_params{2} = guess_rate;

                    % Find best starting parameters (assuming grid_search works for this)
                    free_params = grid_search(unique_probe_offsets, p_CW, guess_rate);

                    % Define the model function handle
                    response_model = @(params) calc_fit(params, fixed_params);

                    % Define lower and upper bounds for mu and sigma if not already defined globally
                    % Example bounds (adjust as needed):
                    lower_bounds = [-20, 0]; % mu can be anything, sigma must be positive
                    upper_bounds = [20, 10]; 
                    options = optimoptions('fmincon','Display','off'); % Suppress fmincon output

                    [params_est, sse, exit_flag] = fmincon(response_model, free_params, [], [], [], [], lower_bounds, upper_bounds, [], options);

                    % Store estimated parameters if optimization was successful
                    if exit_flag > 0
                        mu = params_est(1);
                        sigma = params_est(2);

                        % Calculate R^2
                        x_exp = (1 - guess_rate) * normcdf(unique_probe_offsets, mu, sigma) + 0.5 * guess_rate;
                        ss_res = sum((p_CW - x_exp).^2);
                        ss_tot = sum((p_CW - mean(p_CW)).^2);
                        if ss_tot == 0 % Handle case where p_CW is constant
                           if ss_res < 1e-10 % Perfect fit to constant
                               r2 = 1;
                           else 
                               r2 = 0; % No variance explained
                           end
                        else
                            r2 = 1 - (ss_res / ss_tot);
                        end
                        

                        analyzed_data.ind(subj).mu_by_delta_theta(prev_lvl, curr_lvl, cond, delta_theta_window_indx) = mu;
                        analyzed_data.ind(subj).sigma_by_delta_theta(prev_lvl, curr_lvl, cond, delta_theta_window_indx) = sigma;
                        analyzed_data.ind(subj).r2_by_delta_theta(prev_lvl, curr_lvl, cond, delta_theta_window_indx) = r2;
                    else
                         % Store NaN if fitting failed
                        analyzed_data.ind(subj).mu_by_delta_theta(prev_lvl, curr_lvl, cond, delta_theta_window_indx) = NaN;
                        analyzed_data.ind(subj).sigma_by_delta_theta(prev_lvl, curr_lvl, cond, delta_theta_window_indx) = NaN;
                        analyzed_data.ind(subj).r2_by_delta_theta(prev_lvl, curr_lvl, cond, delta_theta_window_indx) = NaN;                       
                    end


                    % --- Plotting (Commented out to avoid excessive figures) ---
                    x = -20:0.01:20;
                    cdf_response = (1 - guess_rate) * normcdf(x, mu, sigma) + 0.5 * guess_rate;
                    pdf_response = (1 - guess_rate) * normpdf(x, mu, sigma);
                    pdf_response = pdf_response / max(pdf_response); % normalize to peak at 1

                    figure_name = ['S' subj_IDs{subj} '_response_bias_L' num2str(prev_lvl) '-L' num2str(curr_lvl) '_' cond_names{cond} '_DeltaThetaBin_' num2str(delta_theta_window_edges(1, delta_theta_window_indx)) '–' num2str(delta_theta_window_edges(2, delta_theta_window_indx)) '°'];
                    figure('Name', figure_name, 'Visible', 'off'); % Create new figure for each plot if needed
                    plot(x, cdf_response', 'LineWidth', 1,'Color', 'b');
                    hold on;
                    plot(x, pdf_response', 'LineWidth', 1,'LineStyle', '--','Color', 'b');
                    scatter(unique_probe_offsets(1:end/2), p_CW(1:end/2), 30, 'filled','MarkerFaceColor', 'r');
                    scatter(unique_probe_offsets(end/2+1:end), p_CW(end/2+1:end), 30, 'filled','MarkerFaceColor', 'g');

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
                        saveas(gcf, [figure_path '/' figure_name '.png']);
                        disp(['Saved ' figure_name '.png']);
                    end
                    close(gcf); % Close figure after saving
                    %}
                    % --- End Plotting ---

                end % End of delta_theta_window loop

            end % End of cond loop
        end % End of curr_lvl loop
    end % End of prev_lvl loop

end % End of subject loop

if plot_ind, close(fg); end

%% Group


%% Super

