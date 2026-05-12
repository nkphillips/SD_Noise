function plan_row = makeBootstrapPlan(bootstrap_unit, n_subj, n_trials, seed)
% makeBootstrapPlan  Build one bootstrap resampling plan row.

    if nargin < 4 || isempty(seed)
        stream = [];
    else
        stream = RandStream('mt19937ar', 'Seed', seed);
    end

    bootstrap_unit = lower(char(bootstrap_unit));
    switch bootstrap_unit
        case 'subject'
            plan_row = local_randi(stream, n_subj, [1, n_subj]);
        case 'trial'
            plan_row = local_randi(stream, n_trials, [1, n_trials]);
        otherwise
            error('makeBootstrapPlan:badUnit', ...
                'bootstrap_unit must be ''subject'' or ''trial''.');
    end
end

function x = local_randi(stream, imax, sz)
    if isempty(stream)
        x = randi(imax, sz(1), sz(2));
    else
        x = randi(stream, imax, sz(1), sz(2));
    end
end
