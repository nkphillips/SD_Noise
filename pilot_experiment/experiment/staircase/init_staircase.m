%%% init_staircase

% Initialize struct
staircases = struct();
staircases.num_staircases_per_cond = 2;

% Probe offset range
staircases.max_offset = 15;
staircases.min_offset = 1;
staircases.init_probe_offset = [staircases.min_offset staircases.max_offset];

% Step size
staircases.min_step_size = 1;
staircases.max_step_size = 5;
staircases.init_step_size = 5 * ones(1, staircases.num_staircases_per_cond); % both staircases start at with same step size

% Reversals
staircases.max_reversals = 8;
staircases.num_reversals_to_consider = round(staircases.max_reversals * 0.7);

% Trials
staircases.max_trials_per_sc = 50; % default = 50
staircases.num_trials_per_sf = staircases.max_trials_per_sc * staircases.num_staircases_per_cond;

% Arrays
staircases.probe_offsets = nan(staircases.num_staircases_per_cond , staircases.max_trials_per_sc, p.num_levels);
staircases.step_size = nan(staircases.num_staircases_per_cond, staircases.max_trials_per_sc, p.num_levels);
staircases.direction = nan(staircases.num_staircases_per_cond, 2, p.num_levels);

staircases.responses = nan(staircases.num_staircases_per_cond , staircases.max_trials_per_sc, p.num_levels);
staircases.correct_count = zeros(staircases.num_staircases_per_cond, p.num_levels);

staircases.trial_order = nan(staircases.max_trials_per_sc * staircases.num_staircases_per_cond, p.num_levels);
staircases.trial_indices = nan(staircases.num_staircases_per_cond, staircases.max_trials_per_sc, p.num_levels);

staircases.num_reversals = zeros(staircases.num_staircases_per_cond, p.num_levels);
staircases.reversal_indices = nan(staircases.num_staircases_per_cond, staircases.max_trials_per_sc, p.num_levels);

% Initialize staircases
for n_sf = 1:p.num_levels

    % Initialize trial order
    trial_order = repmat(1:staircases.num_staircases_per_cond, 1 , staircases.max_trials_per_sc);
    staircases.trial_order(:, n_sf) = trial_order(randperm(length(trial_order)));

    % Initialize direction
    staircases.direction(:, :, n_sf) = [1 1; -1 -1];

    for n_sc = 1:staircases.num_staircases_per_cond

        staircases.trial_indices(n_sc, :, n_sf) = find(staircases.trial_order(:, n_sf) == n_sc);
        staircases.probe_offsets(n_sc, 1, n_sf) = staircases.init_probe_offset(n_sc);
        staircases.step_size(n_sc, 1, n_sf) = staircases.init_step_size(n_sc);

    end

end
