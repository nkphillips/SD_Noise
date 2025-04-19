function [counts] = getTrialDistribution(p)


    counts.lvl_pair_count = zeros(p.num_levels, p.num_levels, p.num_conds);
    counts.lvl_pair_indices = cell(p.num_levels, p.num_levels, p.num_blocks_per_cond, p.num_conds);
    counts.delta_theta = cell(p.num_levels, p.num_levels, p.num_conds);
    
    save_fig = 1;
    if save_fig
        fg = figure('Color', 'w');
        set(0, 'CurrentFigure', fg);
    end

    for cond = 1:p.num_conds
        
        curr_blocks = find(p.block_order == cond);
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

                            counts.lvl_pair_count(prev_lvl, curr_lvl, cond) = counts.lvl_pair_count(prev_lvl, curr_lvl, cond) + 1;
                            counts.lvl_pair_indices{prev_lvl, curr_lvl, n_block, cond} = [counts.lvl_pair_indices{prev_lvl, curr_lvl, n_block, cond}, n_trial];
    
                        end
                    end

                end
                
                if save_fig && cond == 1

                    figure_name = ['Delta theta distribution (' num2str(p.num_trials_per_cond) ' trials per cond)'];

                    subplot(p.num_levels, p.num_levels, prev_lvl + (curr_lvl-1)*p.num_levels);
                    histogram(counts.delta_theta{prev_lvl, curr_lvl, cond},'BinWidth',2)

                    % Format figure
                    title([num2str(prev_lvl) ' -> ' num2str(curr_lvl)]);
                    if curr_lvl == 3, xlabel('\Delta \theta (°)'); end
                    if prev_lvl == 1, ylabel('Count'); end
                    box off;
                    set(gca,'TickDir','out');
                    xlim([-90 90]);
                    xticks([-90:45:90]);

                end

            end
        end


        if save_fig
            if ~exist('figures', 'dir'), mkdir('figures'); end
            saveas(gcf, ['figures/' figure_name '.png']); 
        end

    end

    if save_fig
        close(fg);
    end

end