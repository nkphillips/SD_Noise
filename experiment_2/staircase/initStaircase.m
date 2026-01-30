function sc = initStaircase(feature, rule, n_trials, stim_bounds)

%% 

sc.rule = rule;
sc.stim_boudnds = stim_bounds;
sc.responses = nan(n_trials,1);
sc.stim_vals = nan(n_trials,1);
sc.step_down = nan(n_trials,1);
sc.step_up = nan(n_trials,1);

%% step size

switch feature
    case 'contrast'
        
        sc.step_down = sc.step_down * sc.step_ratio; 

%%

sc.stim_vals(1) = stim_bounds(2);


end

