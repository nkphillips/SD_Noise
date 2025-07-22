% makeDeltaThetaBins


function delta_theta_bins = makeDeltaThetaBins(bin_edges, all_runs, unique_probe_offsets, num)

    %% Initialize structure

    delta_theta_bins = struct('ind',[],'grp',[],'all',[]);
    
    % Subject level arrays
    struct_size = cell(1,num.subjs);
    delta_theta_bins.ind = struct('delta_thetas', struct_size, 'probe_offsets', struct_size, ...
        'performance', struct_size, 'performance_mean', struct_size, 'performance_sem', struct_size,...
    'pCW', struct_size, 'pCW_mean', struct_size, 'pCW_sem', struct_size, ...
    'pCCW', struct_size, 'pCCW_mean', struct_size, 'pCCW_sem', struct_size,...
    'response_bias', struct_size);

    % Group level arrays (mostly averages)
    struct_size = cell(1);
    delta_theta_bins.grp = struct('performance', struct_size, 'pCW', struct_size, 'pCCW', struct_size, ...
    'response_bias', struct_size);
    
    % Super subject level arrays (similar to subject level but accumulating data across subjects to treat as one subject)
    delta_theta_bins.all = struct('delta_thetas', struct_size, 'probe_offsets', struct_size, 'performance', struct_size, ...
    'pCW', struct_size, 'pCCW', struct_size, 'response_bias', struct_size);

    %% Pre-allocate delta_theta_bins.ind fields

    fieldname_list = fieldnames(delta_theta_bins.ind);
    for subj = 1:num.subjs
        for i_field = 1:numel(fieldname_list)
            if contains(fieldname_list{i_field},{'mean', 'sem'})
                delta_theta_bins.ind(subj).(fieldname_list{i_field}) = nan(num.levels, num.levels, num.conds, num.delta_theta_bins);
            elseif strcmp(fieldname_list{i_field}, {'performance','delta_theta'})
                delta_theta_bins.ind(subj).(fieldname_list{i_field}) = cell(num.levels, num.levels, num.conds, num.delta_theta_bins);
            end
        end
    end
    
    %% Subject level

    for subj = 1:num.subjs

        for n_run = 1:num.runs(subj)

            subj_p = all_runs{subj}(n_run).p;
            behav_data = all_runs{subj}(n_run).behav_data;

            delta_thetas = squeeze(subj_p.trial_events(1:end-1,1,:) - subj_p.trial_events(2:end,1,:));
            probe_offsets = squeeze(subj_p.trial_events(:,2,:) - subj_p.trial_events(:,1,:));

            for cond = 1:num.conds

                curr_cond_blocks = find(subj_p.cond_order == cond);

                for n_block = 1:length(curr_cond_blocks)
                                        
                    curr_correct = behav_data.correct(2:end,curr_cond_blocks(n_block));
                    curr_response = behav_data.response(2:end,curr_cond_blocks(n_block)) == 2; % 0 = CCW; 1 = CW
                    curr_probe_offsets = probe_offsets(2:end,curr_cond_blocks(n_block));
                    curr_delta_theta = delta_thetas(:,curr_cond_blocks(n_block));
                    curr_lvls = subj_p.trial_events(:,3,curr_cond_blocks(n_block));

                    for prev_lvl = 1:num.levels
                        for curr_lvl = 1:num.levels

                            curr_lvl_pair_indx = curr_lvls(1:end-1) == prev_lvl & curr_lvls(2:end) == curr_lvl;

                            for i_bin = 1:num.delta_theta_bins

                                curr_bin_indx = (bin_edges(1,i_bin) <= curr_delta_theta & bin_edges(2,i_bin) >= curr_delta_theta);
                                curr_trials_indx = curr_bin_indx & curr_lvl_pair_indx; % current trials that satisfy the current lvl pair and bin

                                % delta thetas
                                delta_theta_bins.ind(subj).delta_thetas{prev_lvl, curr_lvl, cond, i_bin} = curr_delta_theta(curr_trials_indx);
    
                                % probe offsets
                                delta_theta_bins.ind(subj).probe_offsets{prev_lvl, curr_lvl, cond, i_bin} = curr_probe_offsets(curr_trials_indx);
    
                                % performance 
                                curr_trials_correct = curr_correct(curr_trials_indx);
                                delta_theta_bins.ind(subj).performance{prev_lvl, curr_lvl, cond, i_bin} = curr_trials_correct;
                                delta_theta_bins.ind(subj).performance_mean(prev_lvl, curr_lvl, cond, i_bin) = mean(curr_trials_correct,'omitnan');
                                delta_theta_bins.ind(subj).performance_sem(prev_lvl, curr_lvl, cond, i_bin) = calcSEM(curr_trials_correct);
    
                                % pCW
                                curr_trials_CW = curr_response(curr_trials_indx);
                                delta_theta_bins.ind(subj).pCW{prev_lvl, curr_lvl, cond, i_bin} = curr_trials_CW;
                                delta_theta_bins.ind(subj).pCW_mean(prev_lvl, curr_lvl, cond, i_bin) = mean(curr_trials_CW);
                                delta_theta_bins.ind(subj).pCW_sem(prev_lvl, curr_lvl, cond, i_bin) = calcSEM(curr_trials_CW);

                                % pCCW
                                curr_trials_CCW = ~curr_response(curr_trials_indx);
                                delta_theta_bins.ind(subj).pCCW{prev_lvl, curr_lvl, cond, i_bin} = curr_trials_CCW;
                                delta_theta_bins.ind(subj).pCCW_mean(prev_lvl, curr_lvl, cond, i_bin) = mean(curr_trials_CCW);
                                delta_theta_bins.ind(subj).pCCW_sem(prev_lvl, curr_lvl, cond, i_bin) = calcSEM(curr_trials_CCW);


                                for i_probe_offset =1:length(unique_probe_offsets)
                                    curr_offsets_indx = curr_trials_indx & (curr_probe_offsets == unique_probe_offsets(i_probe_offset));
                                   
                                end
                            end
                        end
                    end

                end

            end

        end
    
    end

    %% Group level



    %% Super subject ("all") level


end