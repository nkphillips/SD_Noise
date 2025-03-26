%%% test_performance

%% Is performance significantly different from chance?

% Performance as a function of test SF
data = cat(3,behav_perf.ind(:).mean_within_test_sf_across_runs);

h = nan(1, num.conds);
p = nan(1, num.conds);
ci = nan(2, num.conds);

for sf = 1:num.test_sfs
    for cond = 1:num.conds

        x = squeeze(data(sf,cond,:));
        y = 0.5;

        [h(sf,cond), p(sf,cond), ci(:,cond,sf)] = ttest(x,y);

    end
end


% Performance as a function of delta SF and test SF
data = cat(4,behav_perf.ind.mean_within_delta_sf);

h = nan(num.test_sfs,num.conds,num.test_sfs);
p = nan(num.test_sfs,num.conds,num.test_sfs);

for sf_t = 1:num.test_sfs
    for cond = 1:num.conds
        for delta = 1:num.test_sfs

            x = squeeze(data(delta,cond,sf_t,:));
            y = 0.5;

            [h(delta,cond,sf_t), p(delta,cond,sf_t)] = ttest(x,y);

        end
    end
end


%% Is performance significantly different between adaptation conditions?

% Overall performance
data = cat(1,behav_perf.ind(:).mean_across_runs);

% 1-way ANOVA
anova_p = anova1(data,cond_names,'off');

% Paired-sample t-test
h = nan(num.conds);
p = nan(num.conds);

for cond_a = 1:num.conds
    for cond_b = 1:num.conds
        if cond_a ~= cond_b

            x = data(:,cond_a);
            y = data(:,cond_b);

            [h(cond_a,cond_b), p(cond_a,cond_b)] = ttest(x,y);

        end

    end
end


% Performance as a function of test SF
data = cat(3,behav_perf.ind(:).mean_within_test_sf_across_runs);

% 1-way ANOVA
anova_p = nan(1,num.test_sfs);
for sf = 1:num.test_sfs

    x = squeeze(data(sf,:,:))';
    anova_p(sf) = anova1(x,cond_names,'off');

end
anova_h = anova_p < 0.05;

% Paired-sample t-test
h = nan(num.conds,num.conds,num.test_sfs);
p = nan(num.conds,num.conds,num.test_sfs);
for cond_a = 1:num.conds
    for cond_b = 1:num.conds
        if cond_a ~= cond_b
            
            for sf = 1:num.test_sfs

                x = squeeze(data(sf,cond_a,:));
                y = squeeze(data(sf,cond_b,:));

                [h(cond_a,cond_b,sf), p(cond_a,cond_b,sf)] = ttest(x,y);

            end

        end
    end
end



% Performance as a function of delta SF
data = cat(4,behav_perf.ind(:).mean_within_delta_sf);

% 1-way ANOVA
anova_p = nan(num.test_sfs);

for sf_i = 1:num.test_sfs
    for sf_j = 1:num.test_sfs

        x = squeeze(data(sf_j,:,sf_i,:))';
        anova_p(sf_j,sf_i) = anova1(x,cond_names,"off");

    end
end
anova_h = anova_p < 0.05;

figure('Color',white,'Name','1-way ANOVA Across Conditions'); 

imagesc(anova_h)

colormap("default");
colorbar;

% Format figure
xticks(1:num.test_sfs)
xtick_labels = num2cell(round(delta_sfs(:,round(num.test_sfs/2)),2));
xticklabels(xtick_labels)
xlabel('SF_{N}')

yticks(1:num.test_sfs)
yticklabels(xtick_labels)
ylabel('SF_{N-1}')

set(gca,'TickDir','out'); axis square; box off;