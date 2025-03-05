%%% init_staircase

p.num_conds_names = {'contrast', 'filter'};
p.num_conds = numel(p.num_conds_names);

% Initialize struct
staircases = struct();
staircases.num_staircases_per_cond = 2;

% Probe offset range
staircases.max_probe_offset = 15;
staircases.min_probe_offset = 1;
staircases.init_probe_offset = [staircases.min_probe_offset staircases.max_probe_offset];

% Step size
staircases.min_step_size = 1;
staircases.max_step_size = 5;
staircases.init_step_size = 5 * ones(1, staircases.num_staircases_per_cond); % both staircases start at with same step size

% Reversals
staircases.max_reversals = 8;
staircases.num_reversals_to_consider = round(staircases.max_reversals * 0.7);

% Trials
staircases.max_trials_per_sc = 50; % default = 50

% Arrays
staircases.probe_offsets = nan(staircases.num_staircases_per_cond , staircases.max_trials_per_sc, p.num_levels, p.num_conds);
staircases.step_size = nan(staircases.num_staircases_per_cond, staircases.max_trials_per_sc, p.num_levels, p.num_conds);
staircases.direction = nan(staircases.num_staircases_per_cond, 2, p.num_levels, p.num_conds);

staircases.test_orientation = nan(staircases.max_trials_per_sc * staircases.num_staircases_per_cond, p.num_levels, p.num_conds);

staircases.responses = nan(staircases.num_staircases_per_cond , staircases.max_trials_per_sc, p.num_levels, p.num_conds);
staircases.correct_count = zeros(staircases.num_staircases_per_cond, p.num_levels, p.num_conds);

staircases.trial_order = nan(staircases.max_trials_per_sc * staircases.num_staircases_per_cond, p.num_levels, p.num_conds);
staircases.trial_indices = nan(staircases.num_staircases_per_cond, staircases.max_trials_per_sc, p.num_levels, p.num_conds);

staircases.num_reversals = zeros(staircases.num_staircases_per_cond, p.num_levels, p.num_conds);
staircases.reversal_indices = nan(staircases.num_staircases_per_cond, staircases.max_trials_per_sc, p.num_levels, p.num_conds);

% Initialize staircases
p.num_trials = numel(staircases.trial_order);
p.num_trials_per_block = staircases.max_trials_per_sc*2;
p.num_blocks = p.num_conds * p.num_levels;

for cond = 1:p.num_conds
    for lvl = 1:p.num_levels

        % Initialize trial order
        trial_order = repmat(1:staircases.num_staircases_per_cond, 1 , staircases.max_trials_per_sc);
        staircases.trial_order(:, lvl, cond) = trial_order(randperm(length(trial_order)));

        % Initialize direction
        staircases.direction(:, :, lvl, cond) = [1 1; -1 -1];

        for n_sc = 1:staircases.num_staircases_per_cond

            staircases.trial_indices(n_sc, :, lvl, cond) = find(staircases.trial_order(:, lvl, cond) == n_sc);

            staircases.probe_offsets(n_sc, 1, lvl, cond) = staircases.init_probe_offset(n_sc); % set the initial probe offset
            staircases.step_size(n_sc, 1, lvl, cond) = staircases.init_step_size(n_sc); % set the initial step size

        end

        % Generate test orientations
        % use the datasample() function to generate a sequence of test orientations from a min and max range of orientation values,
        % equal in size to the staircase structure (eg., staircases.probe_offsets)
        data = stimuli.orientation_min:stimuli.orientation_max;
        staircases.test_orientation(:, lvl, cond) = datasample(data, staircases.max_trials_per_sc * staircases.num_staircases_per_cond); 

    end
end


%

