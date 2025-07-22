function counts = getConditionDistributions(p, behav_data, unique_probe_offsets)

    counts.lvl_pair_count = zeros(p.num_levels, p.num_levels, p.num_conds);
    counts.lvl_pair_indices = cell(p.num_levels, p.num_levels, p.num_blocks_per_cond, p.num_conds);
    counts.delta_theta = cell(p.num_levels, p.num_levels, p.num_conds);
    counts.probe_offset = cell(p.num_levels, p.num_levels, p.num_conds);
    counts.pCW_per_probe_offset = cell(p.num_levels, p.num_levels, p.num_conds, length(unique_probe_offsets));

    for cond = 1:p.num_conds
        
        curr_blocks = find(p.cond_order == cond);
        curr_trials = p.trial_events(:,:,curr_blocks);
        
        for prev_lvl = 1:p.num_levels 
            for curr_lvl = 1:p.num_levels 

                for n_block = 1:length(curr_blocks)
                    for n_trial = 2:size(curr_trials,1)
                        
                        prev_trial_lvl = curr_trials(n_trial-1,3,n_block);
                        curr_trial_lvl = curr_trials(n_trial,3,n_block);
                        
                        if prev_trial_lvl == prev_lvl && curr_trial_lvl == curr_lvl
                            
                            delta_theta = curr_trials(n_trial-1,1,n_block) - curr_trials(n_trial,1,n_block);
                            counts.delta_theta{prev_lvl, curr_lvl, cond} = [counts.delta_theta{prev_lvl, curr_lvl, cond}, delta_theta];

                            probe_offset = curr_trials(n_trial,2,n_block) - curr_trials(n_trial,1,n_block);
                            counts.probe_offset{prev_lvl, curr_lvl, cond} = [counts.probe_offset{prev_lvl, curr_lvl, cond}, probe_offset];

                            offset_indx = find(unique_probe_offsets == probe_offset);
                            counts.pCW_per_probe_offset{prev_lvl, curr_lvl, cond, offset_indx} = [counts.pCW_per_probe_offset{prev_lvl, curr_lvl, cond, offset_indx}, behav_data.response(n_trial, curr_blocks(n_block)) == 2];

                            counts.lvl_pair_count(prev_lvl, curr_lvl, cond) = counts.lvl_pair_count(prev_lvl, curr_lvl, cond) + 1;
                            counts.lvl_pair_indices{prev_lvl, curr_lvl, n_block, cond} = [counts.lvl_pair_indices{prev_lvl, curr_lvl, n_block, cond}, n_trial];
    
                        end
                    end
                
                end

            end
        end



    end


end