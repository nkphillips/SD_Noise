%%% plot_performance_ind

%% Set figure path

figure_path = [ind_figure_path '/' subj_IDs{subj} '/performance/main'];

if ~exist(figure_path,'dir')
    mkdir(figure_path)
    disp([figure_path ' created.'])
end

%% Condition vs Performance
%
figure_name = ['S' subj_IDs{subj} ' Performance Per Condition'];

x = 1:(length(behav_perf.ind(subj).mean_across_runs) + 1);
y = [behav_perf.ind(subj).mean_across_runs behav_perf.ind(subj).mean_across_conds] * 100;
err = [behav_perf.ind(subj).sem_across_runs behav_perf.ind(subj).sem_across_conds] * 100;

h = bar(x,y);
hold on;
errorbar(x,y,err,'LineStyle','none','CapSize',0,'LineWidth',line_width,'Color',black)

% Format figure
text(max(xlim),max(ylim),['n = ' num2str(num.runs(subj))],'FontSize',12);
h.FaceColor = 'flat';
h.CData = [cond_colors; black]; h.EdgeColor = white;
set(gca,'TickDir','out'); box off

ylim([0 100])
ylabel('% Test SF reported as higher')
xlabel('Condition')

xticks(x)
xticklabels([cond_names 'Total'])

if save_ind_figures, saveas(gcf,[figure_path '/' figure_name figure_file_type]); end
clf
%}
%% Test SF vs Performance
%
figure_name = ['S' subj_IDs{subj} ' Performance Per Test SF'];

x = sf_ratios;
y = behav_perf.ind(subj).mean_within_test_sf_across_runs * 100;
err = behav_perf.ind(subj).sem_within_test_sf_across_runs * 100;

% h = errorbar(x,y,err,'o','MarkerSize',marker_size,'CapSize',0,'LineStyle','-');
h = errorbar(x,y,err,'CapSize',0,'LineStyle','-');
hold on;

% Format figure
for i_h = 1:numel(h)
    h(i_h).Color = cond_colors(i_h,:);
    h(i_h).MarkerFaceColor = cond_colors(i_h,:);
    h(i_h).MarkerEdgeColor = white;
end
set(gca,'TickDir','out'); box off

ylim([0 100])
ylabel('% Test SF reported as higher')
xlabel({'SF ratio (octaves)','log_2[SF_T/SF_R]'})

xticks(x)
xticklabels(num2cell(round(x,2)));

% line([0 0],[min(ylim) max(ylim)],'LineStyle','-','Color', black)
line([min(xlim) max(xlim)], [50 50], 'LineStyle', '--', 'Color', black)

text(max(xlim),min(ylim)+max(ylim)/10,['n = ' num2str(num.runs(subj))],'FontSize',12);

if save_ind_figures, saveas(gcf,[figure_path '/' figure_name figure_file_type]); end
clf
%}
%% Delta SF vs Performance
%{
figure_name = ['S' subj_IDs{subj} ' Delta SF vs Performance'];

[subplot_x, subplot_y] = get_subplot_dimensions(num.test_sfs);

for n_sf = 1:num.test_sfs

    subplot(subplot_x, subplot_y, n_sf);

    x = delta_sfs(:,n_sf)';
    y = behav_perf.ind(subj).mean_within_delta_sf(:,:,n_sf) * 100;
    err = behav_perf.ind(subj).sem_within_delta_sf(:,:,n_sf) * 100;

    % h = errorbar(x,y,err,'o','MarkerSize',marker_size,'CapSize',0,'LineStyle','-');
    h = errorbar(x,y,err,'CapSize',0,'LineStyle','-');

    % Format figure
    title(num2str(round(sf_ratios(n_sf),2)));
    for i_h = 1:numel(h)
        h(i_h).Color = cond_colors(i_h,:);
        h(i_h).MarkerFaceColor = cond_colors(i_h,:);
        h(i_h).MarkerEdgeColor = white;
    end
    set(gca,'TickDir','out'); box off
    
    xlim([min_delta_sf max_delta_sf])
    ylim([0 100])

    if n_sf == 1, xlabel({'\Delta SF (octaves)','log_2[SF_{N-1}/SF_N]'}); ylabel('% Test SF reported as higher'); end
       
    xticks(min(xlim):0.5:max(xlim))
    % xticklabels(num2cell(round(delta_sfs(:,n_sf),2))');
    
    % line([0 0],[min(ylim) max(ylim)],'LineStyle','-','Color', black)
    line([min(xlim) max(xlim)], [50 50], 'LineStyle', '--', 'Color', black)

end

if save_ind_figures, saveas(gcf,[figure_path '/' figure_name figure_file_type]); end
clf
%}

%% Delta SF vs Performance (Per condition)
%{

for cond = 1:num.conds

    figure_name = ['S' subj_IDs{subj} ' Delta SF vs Performance (' cond_names{cond} ')'];

    [subplot_x, subplot_y] = get_subplot_dimensions(num.test_sfs);

    for n_sf = 1:num.test_sfs

        subplot(subplot_x, subplot_y, n_sf);

        x = delta_sfs(:,n_sf)';
        y = behav_perf.ind(subj).mean_within_delta_sf(:,cond,n_sf) * 100;
        err = behav_perf.ind(subj).sem_within_delta_sf(:,cond,n_sf) * 100;

        % h = errorbar(x,y,err,'o','MarkerSize',marker_size,'CapSize',0,'LineStyle','-');
        h = errorbar(x,y,err,'CapSize',0,'LineStyle','-');

        % Format figure
        title(num2str(round(sf_ratios(n_sf),2)));
        h.Color = cond_colors(cond,:);
        h.MarkerFaceColor = cond_colors(cond,:);
        h.MarkerEdgeColor = white;
        set(gca,'TickDir','out'); box off

        xlim([min_delta_sf max_delta_sf])
        ylim([0 100])

        if n_sf == 1, xlabel({'\Delta SF (octaves)','log_2[SF_{N-1}/SF_N]'}); ylabel('% Test SF reported as higher'); end

        xticks(min(xlim):0.5:max(xlim))
        % xticklabels(num2cell(round(delta_sfs(:,n_sf),2))');

        % line([0 0],[min(ylim) max(ylim)],'LineStyle','-','Color', black)
        line([min(xlim) max(xlim)], [50 50], 'LineStyle', '--', 'Color', black)

    end

    if save_ind_figures, saveas(gcf,[figure_path '/' figure_name figure_file_type]); end
    clf

end
%}


%% Trials per Delta SF

figure_name = ['S' subj_IDs{subj} ' Number of Trials per Delta SF'];

for cond = 1:num.conds

    subplot(1,num.conds,cond)

    x = behav_perf.ind(subj).trials_per_delta_sf(:,:,cond);

    imagesc(x);
    colormap('jet');  
    colorbar;

    % Format figure
    xticks(1:num.test_sfs)
    xtick_labels = num2cell(round(delta_sfs(:,round(num.test_sfs/2)),2));
    xticklabels(xtick_labels)
    xlabel('Test SF')

    yticks(1:num.test_sfs)
    yticklabels(xtick_labels)
    ylabel('Test SF')

    set(gca,'TickDir','out'); axis square; box off;

end

if save_ind_figures, saveas(gcf,[figure_path '/' figure_name figure_file_type]); end
clf