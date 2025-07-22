
if plot_ind

    fg = figure('Color', 'w');
    set(0, 'CurrentFigure', fg);
    figure_path = [ind_figure_path '/S' p.subj_ID];

    for cond = 1:p.num_conds
    
        figure_name = ['Delta theta distribution (' num2str(p.num_trials_per_cond) ' trials per cond) ' cond_names{cond}];

        for prev_lvl = 1:p.num_levels
            for curr_lvl = 1:p.num_levels
                
                curr_counts = counts.delta_theta{prev_lvl, curr_lvl, cond};

                subplot(p.num_levels, p.num_levels, prev_lvl + (curr_lvl-1)*p.num_levels);
                histogram(curr_counts,'BinWidth',3)

                % Format figure
                title([num2str(prev_lvl) ' -> ' num2str(curr_lvl)]);
                if curr_lvl == 3, xlabel('\Delta \theta (°)'); end
                if prev_lvl == 1, ylabel('Count'); end
                box off;
                set(gca,'TickDir','out');
                xlim([-90 90]);
                xticks(-90:45:90);

            end
        end
    end

    if save_ind_figures
        if ~exist(figure_path, 'dir'), mkdir(figure_path); end
        saveas(gcf, [figure_path '/' figure_name '.png']); 
        disp(['Saved ' figure_name '.png']);
    end

end
