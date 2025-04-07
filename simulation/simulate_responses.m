function responses = simulate_responses(p)

    % Calculate orientation difference
    orient_diff = p.trial_events(:, 2) - p.trial_events(:, 1);

    % Calculate difference between previous and current trial orientation
    delta_theta = p.trial_events(1:end-1,1) - p.trial_events(2:end,1);

    % Get response bias from serial dependence
    bias = nan(length(delta_theta),1);
    for n = 1:length(delta_theta)
        curr_level = p.trial_events(n,end);
        bias(n) = gaussian_prime(p.DoG_params(curr_level,:), delta_theta(n));
    end
    bias = [0; bias];

    % Get internal orientation difference (actual diff + bias + noise)
    internal_orient_diff = orient_diff + bias + p.DoG_params(3) .* randn(length(orient_diff), 1);

    % Calculate probability of correct discrimination using cumulative normal
    mu = p.psychometric_params(1); % ie threshold
    sigma = p.psychometric_params(2); % ie slope
    guess_rate = p.psychometric_params(3);
    p_correct = (1 - guess_rate) * normcdf(internal_orient_diff, mu, sigma) + 0.5 * guess_rate;
    
    % Simulate response based on discrimination probability
    responses = rand() < p_correct;

end



% === Compute probe offset and condition labels ===




trial_events = experiments(n_exp).p(subj).trial_events;

probe_offsets = trial_events(:,2) - trial_events(:,1);  % Test - Ref

contrast_levels = stimuli.contrast(trial_events(:,3));

filter_widths = stimuli.bp_filter_width(trial_events(:,3));

responses = experiments(n_exp).subj_data(subj).responses;

% === Build table of responses ===
resp_table = table(probe_offsets, contrast_levels, filter_widths, responses, ...
    'VariableNames', {'Offset', 'Contrast', 'Filter', 'Response'});

% === Summarize p(CW) ===
summary = groupsummary(resp_table, {'Offset', 'Contrast', 'Filter'}, 'mean', 'Response');

summary.Properties.VariableNames{'mean_Response'} = 'respCW';

% === Store for this subject ===
experiments(n_exp).subj_data(subj).summary = summary;



%% Gridsearch

function [best_params, best_sse] = grid_search_mu_sigma(offsets, pCW)
    mu_vals = -10:0.5:10;
    sigma_vals = 0.5:0.5:10;
    best_sse = inf;

    for mu = mu_vals
        for sigma = sigma_vals
            preds = normcdf(offsets, mu, sigma);
            sse = sum((pCW - preds).^2);
            if sse < best_sse
                best_sse = sse;
                best_params = [mu, sigma];
            end
        end
    end
end

%% Fit Subject Data

summary = experiments(n_exp).subj_data(subj).summary;
offsets = summary.Offset;
pCW = summary.respCW;

% Initial fit
[init_params, ~] = grid_search_mu_sigma(offsets, pCW);

% Refine using nonlinear optimization
fitfun = @(params, x) normcdf(x, params(1), params(2));
errfun = @(params) sum((pCW - fitfun(params, offsets)).^2);

opt_params = fminsearch(errfun, init_params);
experiments(n_exp).subj_data(subj).fit_params = opt_params;

%% Evaluate Param Estimates
predicted = normcdf(offsets, opt_params(1), opt_params(2));
SSE = sum((pCW - predicted).^2);
R2 = 1 - SSE / sum((pCW - mean(pCW)).^2);

experiments(n_exp).subj_data(subj).fit_stats = struct('SSE', SSE, 'R2', R2);

%% CDF PDF Plots
x = linspace(min(offsets), max(offsets), 100);
cdf_fit = normcdf(x, opt_params(1), opt_params(2));
pdf_fit = normpdf(x, opt_params(1), opt_params(2));

figure;
subplot(1,2,1); hold on;
plot(x, cdf_fit, 'LineWidth', 2);
scatter(offsets, pCW, 'filled');
xlabel('Probe Offset'); ylabel('p(CW)');
title('Psychometric CDF Fit');

subplot(1,2,2);
plot(x, pdf_fit, 'LineWidth', 2);
xlabel('Probe Offset'); ylabel('Density');
title('PDF of Fit');

