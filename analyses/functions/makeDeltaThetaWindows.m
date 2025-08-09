% makeDeltaThetaWindows

function delta_theta_windows = makeDeltaThetaWindows(delta_theta_centers, delta_theta_width, all_runs, num, p, plt_settings)

    %% Initialize structure

    delta_theta_windows = struct('ind',[],'grp',[],'all',[]);
    
    % Subject level arrays
    struct_size = cell(1,num.subjs);
    delta_theta_windows.ind = struct('delta_thetas', struct_size, 'probe_offsets', struct_size, ...
        'responses', struct_size, 'pCW', struct_size, 'pCW_sem', struct_size, ...
        'correct_trials', struct_size, 'performance', struct_size, 'performance_sem', struct_size);

    % Group level arrays (mostly averages)
    struct_size = cell(1);
    delta_theta_windows.grp = struct('performance', struct_size, 'pCW', struct_size, 'pCW_sem', struct_size);
    
    % Super subject level arrays (similar to subject level but accumulating data across subjects to treat as one subject)
    delta_theta_windows.all = struct('delta_thetas', struct_size, 'probe_offsets', struct_size, ...
    'responses', struct_size, 'pCW', struct_size, 'pCW_sem', struct_size, ...
    'correct_trials', struct_size, 'performance', struct_size, 'performance_sem', struct_size);

    %% Pre-allocate delta_theta_windows.ind and delta_theta_windows.all fields
    % Performance and sem are arrays, as they present summary statistics
    % Delta thetas, probe offsets, responses, and correct trials are cells, as they present individual trials, which vary in length for each condition

    fieldname_list = fieldnames(delta_theta_windows.ind);

    for subj = 1:num.subjs
        for i_field = 1:numel(fieldname_list)

            if contains(fieldname_list{i_field},{'performance','sem'})
               
                delta_theta_windows.ind(subj).(fieldname_list{i_field}) = nan(num.levels, num.levels, num.conds, num.delta_theta_windows);
                if subj == 1
                    delta_theta_windows.all.(fieldname_list{i_field}) = nan(num.levels, num.levels, num.conds, num.delta_theta_windows);
                end
            
            elseif contains(fieldname_list{i_field}, {'delta_thetas','probe_offsets','responses','correct_trials'})
             
                delta_theta_windows.ind(subj).(fieldname_list{i_field}) = cell(num.levels, num.levels, num.conds, num.delta_theta_windows);
                if subj == 1
                    delta_theta_windows.all.(fieldname_list{i_field}) = cell(num.levels, num.levels, num.conds, num.delta_theta_windows);
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
            delta_thetas = calcOrientationDiff(squeeze(subj_p.trial_events(1:end-1,1,:)), squeeze(subj_p.trial_events(2:end,1,:)));
            probe_offsets = calcOrientationDiff(subj_p.trial_events(:,2,:), subj_p.trial_events(:,1,:));

            for cond = 1:num.conds

                curr_cond_blocks = find(subj_p.cond_order == cond);

                for n_block = 1:length(curr_cond_blocks)
                                        
                    curr_delta_thetas = delta_thetas(:,curr_cond_blocks(n_block));
                    curr_lvls = subj_p.trial_events(:,3,curr_cond_blocks(n_block));
                    curr_probe_offsets = probe_offsets(2:end,curr_cond_blocks(n_block));

                    curr_CW_response = behav_data.response(2:end,curr_cond_blocks(n_block)) == 2; % originally, 1 = CCW; 2 = CW, which becomes 0 = CCW and 1 = CW

                    for prev_lvl = 1:num.levels
                        for curr_lvl = 1:num.levels

                            curr_lvl_pair_indx = curr_lvls(1:end-1) == prev_lvl & curr_lvls(2:end) == curr_lvl;
                            % Accumulate per level pair once for overall histogram (avoid double-counting across windows)
                            pair_delta_thetas = curr_delta_thetas(curr_lvl_pair_indx);
                            all_delta_thetas{prev_lvl, curr_lvl, cond} = [all_delta_thetas{prev_lvl, curr_lvl, cond}; pair_delta_thetas];

                            for i_window = 1:length(delta_theta_centers)

                                % Get current trials that satisfy the current lvl pair and window
                                curr_window_indx = (delta_theta_centers(i_window) - delta_theta_width/2 <= curr_delta_thetas & delta_theta_centers(i_window) + delta_theta_width/2 >= curr_delta_thetas);
                                curr_trials_indx = curr_window_indx & curr_lvl_pair_indx; % current trials that satisfy the current lvl pair and window

                                % Delta thetas
                                delta_theta_windows.ind(subj).delta_thetas{prev_lvl, curr_lvl, cond, i_window} = curr_delta_thetas(curr_trials_indx);
                                delta_theta_windows.all.delta_thetas{prev_lvl, curr_lvl, cond, i_window} = [delta_theta_windows.all.delta_thetas{prev_lvl, curr_lvl, cond, i_window}; ...
                                    curr_delta_thetas(curr_trials_indx)];

                                % Probe offsets
                                delta_theta_windows.ind(subj).probe_offsets{prev_lvl, curr_lvl, cond, i_window} = curr_probe_offsets(curr_trials_indx);
                                delta_theta_windows.all.probe_offsets{prev_lvl, curr_lvl, cond, i_window} = [delta_theta_windows.all.probe_offsets{prev_lvl, curr_lvl, cond, i_window}; ...
                                    curr_probe_offsets(curr_trials_indx)];

                                % Responses
                                CW_response = curr_CW_response(curr_trials_indx);
                                delta_theta_windows.ind(subj).responses{prev_lvl, curr_lvl, cond, i_window} = CW_response;
                                delta_theta_windows.all.responses{prev_lvl, curr_lvl, cond, i_window} = [delta_theta_windows.all.responses{prev_lvl, curr_lvl, cond, i_window}; ...
                                    CW_response];

                            end
                        end
                    end

                end

            end

        end
    
    end

    %% Generate histogram of delta thetas

    %%% Super subject level %%%

    % Define consistent bin edges and find max count across subplots (ignoring windows)
    edges = -90:5:90; % 5-degree bins across [-90, 90]
    max_y = 0;
    for cond = 1:num.conds
        for prev_lvl = 1:num.levels 
            for curr_lvl = 1:num.levels
                counts = histcounts(all_delta_thetas{prev_lvl, curr_lvl, cond}, edges);
                if ~isempty(counts)
                    max_y = max(max_y, max(counts));
                end
            end
        end
    end

    % Plot
    if plt_settings.plot_sup_figures
        for cond = 1:num.conds

            fg_name = ['Super Subj Delta Thetas ' p.cond_names{cond}];

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

                    subplot(num.levels, num.levels, prev_lvl + (curr_lvl-1)*num.levels);

                    % Plot histogram
                    histogram(all_delta_thetas{prev_lvl, curr_lvl, cond}, 'BinEdges', edges, 'Normalization', 'count');

                    % Format figure
                    axis square;
                    title(fg_title);
                    
                    % Only show x-axis labels on bottom row
                    if curr_lvl == num.levels
                        xlabel('\Delta\theta (°)');
                        xticks(-90:45:90);
                    else
                        xticks(-90:45:90);
                        set(gca, 'XTickLabel', []);
                    end
                    
                    % Only show y-axis label on bottom left subplot
                    if prev_lvl == 1 && curr_lvl == num.levels
                        ylabel('Count');
                    end
                    
                    set(gca, 'TickDir', 'out', 'TickLength', [plt_settings.tick_length, plt_settings.tick_length]);
                    xlim([-90 90]);
                    ylim([0, max_y]);
                    line([0, 0], [0, max_y], 'LineWidth', plt_settings.line_width, 'Color', plt_settings.colors.black);
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
    end

end