function packed = generateBootstrapSample(tbl_trials, bootstrap_plan_row, bootstrap_unit)
% generateBootstrapSample  Dispatch subject-cluster or trial-level resampling.

    if nargin < 3 || isempty(bootstrap_unit)
        bootstrap_unit = 'subject';
    end
    bootstrap_unit = lower(char(bootstrap_unit));

    switch bootstrap_unit
        case 'subject'
            packed = generateClusterSample(tbl_trials, bootstrap_plan_row);
        case 'trial'
            packed = generateTrialSample(tbl_trials, bootstrap_plan_row);
        otherwise
            error('generateBootstrapSample:badUnit', ...
                'bootstrap_unit must be ''subject'' or ''trial''.');
    end
end
