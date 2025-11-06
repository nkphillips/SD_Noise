% makeDeltaThetaWindows

function [delta_theta_windows, all_delta_thetas] = makeDeltaThetaWindows(delta_theta_centers, delta_theta_width, all_runs, num, p, plt_settings, n_back)

%% Initialize structure

% Calculate number of delta theta windows
num.delta_theta_windows = length(delta_theta_centers);

delta_theta_windows = struct('ind',[],'grp',[],'all',[]);

% Subject level arrays
struct_size = cell(1,num.subjs);
delta_theta_windows.ind = struct('delta_thetas', struct_size, 'probe_offsets', struct_size, ...
    'responses', struct_size, 'pCW', struct_size, 'pCW_sem', struct_size, ...
    'correct_trials', struct_size, 'performance', struct_size, 'performance_sem', struct_size);

% Group level arrays (averages across subject-level averages)
struct_size = cell(1);
delta_theta_windows.grp = struct('delta_thetas', struct_size, 'probe_offsets', struct_size, ...
    'responses', struct_size, 'pCW', struct_size, 'pCW_sem', struct_size, ...
    'correct_trials', struct_size, 'performance', struct_size, 'performance_sem', struct_size);

% Super subject level arrays (similar to subject level but accumulating data across subjects to treat as one subject)
delta_theta_windows.all = struct('delta_thetas', struct_size, 'probe_offsets', struct_size, ...
    'responses', struct_size, 'pCW', struct_size, 'pCW_sem', struct_size, ...
    'correct_trials', struct_size, 'performance', struct_size, 'performance_sem', struct_size);

%% Pre-allocate delta_theta_windows.ind, delta_theta_windows.grp, and delta_theta_windows.all fields
% Performance and sem are arrays, as they present summary statistics
% Delta thetas, probe offsets, responses, and correct trials are cells, as they present individual trials, which vary in length for each condition

fieldname_list = fieldnames(delta_theta_windows.ind);

for subj = 1:num.subjs
    for i_field = 1:numel(fieldname_list)

        if contains(fieldname_list{i_field},{'performance','sem','pCW'})

            delta_theta_windows.ind(subj).(fieldname_list{i_field}) = nan(num.levels, num.levels, num.conds, num.delta_theta_windows);
            if subj == 1
                delta_theta_windows.all.(fieldname_list{i_field}) = nan(num.levels, num.levels, num.conds, num.delta_theta_windows);
                delta_theta_windows.grp.(fieldname_list{i_field}) = nan(num.levels, num.levels, num.conds, num.delta_theta_windows);
            end

        elseif contains(fieldname_list{i_field}, {'delta_thetas','probe_offsets','responses','correct_trials'})

            delta_theta_windows.ind(subj).(fieldname_list{i_field}) = cell(num.levels, num.levels, num.conds, num.delta_theta_windows);
            if subj == 1
                delta_theta_windows.all.(fieldname_list{i_field}) = cell(num.levels, num.levels, num.conds, num.delta_theta_windows);
                delta_theta_windows.grp.(fieldname_list{i_field}) = cell(num.levels, num.levels, num.conds, num.delta_theta_windows);
            end

        end

    end

end


%% Subject level and 'all' level

all_delta_thetas = cell(num.levels, num.levels, num.conds); % Store all delta thetas for each condition for histogram

for subj = 1:num.subjs

    for n_run = 1:num.runs(subj)

        subj_p = all_runs{subj}(n_run).p;
        behav_data = all_runs{subj}(n_run).behav_data;

        % Calculate orientation difference (range: -90 to 90)
        delta_thetas = calcOrientationDiff(squeeze(subj_p.trial_events(1:end-n_back,1,:)), squeeze(subj_p.trial_events(n_back+1:end,1,:)));
        probe_offsets = calcOrientationDiff(subj_p.trial_events(:,2,:), subj_p.trial_events(:,1,:));

        for cond = 1:num.conds

            curr_cond_blocks = find(subj_p.cond_order == cond);

            for n_block = 1:length(curr_cond_blocks)

                curr_delta_thetas = delta_thetas(:,curr_cond_blocks(n_block));
                curr_lvls = subj_p.trial_events(:,3,curr_cond_blocks(n_block));
                curr_probe_offsets = probe_offsets(n_back+1:end,curr_cond_blocks(n_block));

                curr_CW_response = behav_data.response(n_back+1:end,curr_cond_blocks(n_block)) == 2; % originally, 1 = CCW; 2 = CW, which becomes 0 = CCW and 1 = CW

                for prev_lvl = 1:num.levels
                    for curr_lvl = 1:num.levels

                        curr_lvl_pair_indx = curr_lvls(1:end-n_back) == prev_lvl & curr_lvls(n_back+1:end) == curr_lvl;
                        % Accumulate per level pair once for overall histogram (avoid double-counting across windows)
                        pair_delta_thetas = curr_delta_thetas(curr_lvl_pair_indx);
                        all_delta_thetas{prev_lvl, curr_lvl, cond} = [all_delta_thetas{prev_lvl, curr_lvl, cond}; pair_delta_thetas];

                        for i_window = 1:num.delta_theta_windows

                            left_edge = delta_theta_centers(i_window) - delta_theta_width/2;
                            right_edge = delta_theta_centers(i_window) + delta_theta_width/2;

                            % Get current trials that satisfy the current lvl pair and window
                            if left_edge < -90
                                % Wrap around from positive to negative
                                curr_window_indx = curr_delta_thetas >= (left_edge + 180) | curr_delta_thetas <= right_edge;
                            elseif right_edge > 90
                                % Wrap around from negative to positive
                                curr_window_indx = curr_delta_thetas <= (right_edge - 180) | curr_delta_thetas >= left_edge;
                            else
                                % No wrapping
                                curr_window_indx = curr_delta_thetas >= left_edge & curr_delta_thetas <= right_edge;
                            end

                            curr_trials_indx = curr_window_indx & curr_lvl_pair_indx; % current trials that satisfy the current lvl pair and window

                            % Delta thetas
                            delta_theta_windows.ind(subj).delta_thetas{prev_lvl, curr_lvl, cond, i_window} = curr_delta_thetas(curr_trials_indx);
                            delta_theta_windows.all.delta_thetas{prev_lvl, curr_lvl, cond, i_window} = [delta_theta_windows.all.delta_thetas{prev_lvl, curr_lvl, cond, i_window}; ...
                                curr_delta_thetas(curr_trials_indx)];
                            % Group level: aggregate data across subjects
                            delta_theta_windows.grp.delta_thetas{prev_lvl, curr_lvl, cond, i_window} = [delta_theta_windows.grp.delta_thetas{prev_lvl, curr_lvl, cond, i_window}; ...
                                curr_delta_thetas(curr_trials_indx)];

                            % Probe offsets
                            delta_theta_windows.ind(subj).probe_offsets{prev_lvl, curr_lvl, cond, i_window} = curr_probe_offsets(curr_trials_indx);
                            delta_theta_windows.all.probe_offsets{prev_lvl, curr_lvl, cond, i_window} = [delta_theta_windows.all.probe_offsets{prev_lvl, curr_lvl, cond, i_window}; ...
                                curr_probe_offsets(curr_trials_indx)];
                            % Group level: aggregate data across subjects
                            delta_theta_windows.grp.probe_offsets{prev_lvl, curr_lvl, cond, i_window} = [delta_theta_windows.grp.probe_offsets{prev_lvl, curr_lvl, cond, i_window}; ...
                                curr_probe_offsets(curr_trials_indx)];

                            % Responses
                            CW_response = curr_CW_response(curr_trials_indx);
                            delta_theta_windows.ind(subj).responses{prev_lvl, curr_lvl, cond, i_window} = CW_response;
                            delta_theta_windows.all.responses{prev_lvl, curr_lvl, cond, i_window} = [delta_theta_windows.all.responses{prev_lvl, curr_lvl, cond, i_window}; ...
                                CW_response];
                            % Group level: aggregate data across subjects
                            delta_theta_windows.grp.responses{prev_lvl, curr_lvl, cond, i_window} = [delta_theta_windows.grp.responses{prev_lvl, curr_lvl, cond, i_window}; ...
                                CW_response];

                        end
                    end
                end

            end

        end

    end

end

%% Calculate performance metrics for all levels

% Individual subject and super-subject performance
for prev_lvl = 1:num.levels
    for curr_lvl = 1:num.levels
        for cond = 1:num.conds
            for i_window = 1:num.delta_theta_windows

                % Individual subject performance
                for subj = 1:num.subjs
                    subj_resp = delta_theta_windows.ind(subj).responses{prev_lvl, curr_lvl, cond, i_window};
                    subj_probe = delta_theta_windows.ind(subj).probe_offsets{prev_lvl, curr_lvl, cond, i_window};

                    if ~isempty(subj_resp) && ~isempty(subj_probe)
                        % Match lengths
                        n = min(length(subj_resp), length(subj_probe));
                        subj_resp = subj_resp(1:n);
                        subj_probe = subj_probe(1:n);

                        % Calculate percent correct
                        is_correct = (subj_resp == (subj_probe > 0));
                        delta_theta_windows.ind(subj).performance(prev_lvl, curr_lvl, cond, i_window) = 100 * mean(is_correct, 'omitnan');

                        % Calculate percent CCW (pCW)
                        delta_theta_windows.ind(subj).pCW(prev_lvl, curr_lvl, cond, i_window) = 100 * mean(1 - subj_resp, 'omitnan');

                        % Calculate SEM
                        n_trials = length(subj_resp);
                        if n_trials > 0
                            p_correct = mean(is_correct, 'omitnan');
                            delta_theta_windows.ind(subj).performance_sem(prev_lvl, curr_lvl, cond, i_window) = 100 * sqrt(p_correct * (1 - p_correct) / n_trials);

                            p_ccw = mean(1 - subj_resp, 'omitnan');
                            delta_theta_windows.ind(subj).pCW_sem(prev_lvl, curr_lvl, cond, i_window) = 100 * sqrt(p_ccw * (1 - p_ccw) / n_trials);
                        end
                    end
                end

                % Super-subject performance (from aggregated data)
                all_resp = delta_theta_windows.all.responses{prev_lvl, curr_lvl, cond, i_window};
                all_probe = delta_theta_windows.all.probe_offsets{prev_lvl, curr_lvl, cond, i_window};

                if ~isempty(all_resp) && ~isempty(all_probe)
                    % Match lengths
                    n = min(length(all_resp), length(all_probe));
                    all_resp = all_resp(1:n);
                    all_probe = all_probe(1:n);

                    % Calculate percent correct
                    is_correct = (all_resp == (all_probe > 0));
                    delta_theta_windows.all.performance(prev_lvl, curr_lvl, cond, i_window) = 100 * mean(is_correct, 'omitnan');

                    % Calculate percent CCW (pCW)
                    delta_theta_windows.all.pCW(prev_lvl, curr_lvl, cond, i_window) = 100 * mean(1 - all_resp, 'omitnan');

                    % Calculate SEM
                    n_trials = length(all_resp);
                    if n_trials > 0
                        p_correct = mean(is_correct, 'omitnan');
                        delta_theta_windows.all.performance_sem(prev_lvl, curr_lvl, cond, i_window) = 100 * sqrt(p_correct * (1 - p_correct) / n_trials);

                        p_ccw = mean(1 - all_resp, 'omitnan');
                        delta_theta_windows.all.pCW_sem(prev_lvl, curr_lvl, cond, i_window) = 100 * sqrt(p_ccw * (1 - p_ccw) / n_trials);
                    end
                end

                % Group-level performance (average across subject-level averages)
                subj_performance = nan(num.subjs, 1);
                subj_pCW = nan(num.subjs, 1);
                subj_performance_sem = nan(num.subjs, 1);
                subj_pCW_sem = nan(num.subjs, 1);

                for subj = 1:num.subjs
                    subj_performance(subj) = delta_theta_windows.ind(subj).performance(prev_lvl, curr_lvl, cond, i_window);
                    subj_pCW(subj) = delta_theta_windows.ind(subj).pCW(prev_lvl, curr_lvl, cond, i_window);
                    subj_performance_sem(subj) = delta_theta_windows.ind(subj).performance_sem(prev_lvl, curr_lvl, cond, i_window);
                    subj_pCW_sem(subj) = delta_theta_windows.ind(subj).pCW_sem(prev_lvl, curr_lvl, cond, i_window);
                end

                % Average across subjects (group-level average)
                delta_theta_windows.grp.performance(prev_lvl, curr_lvl, cond, i_window) = mean(subj_performance, 'omitnan');
                delta_theta_windows.grp.pCW(prev_lvl, curr_lvl, cond, i_window) = mean(subj_pCW, 'omitnan');

                % Group-level SEM: standard error of the mean across subjects
                valid_subjs = ~isnan(subj_performance);
                if sum(valid_subjs) > 1
                    delta_theta_windows.grp.performance_sem(prev_lvl, curr_lvl, cond, i_window) = std(subj_performance(valid_subjs), 'omitnan') / sqrt(sum(valid_subjs));
                    delta_theta_windows.grp.pCW_sem(prev_lvl, curr_lvl, cond, i_window) = std(subj_pCW(valid_subjs), 'omitnan') / sqrt(sum(valid_subjs));
                end

            end
        end
    end
end

end