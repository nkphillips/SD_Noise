function [pfit, exitf] = fitConditionByMethod(tbl_cond, fit_opts, fit_method)
% fitConditionByMethod  Dispatch one condition fit by requested estimator.

    if nargin < 2 || isempty(fit_opts)
        fit_opts = struct();
    end
    if nargin < 3 || isempty(fit_method)
        fit_method = 'pooled';
    end

    fit_method = lower(char(fit_method));
    switch fit_method
        case 'pooled'
            [pfit, exitf] = fitConditionMLE( ...
                tbl_cond.delta_theta, tbl_cond.x_probe, tbl_cond.response, fit_opts);

        case 'hierarchical_map'
            [pfit, subject_params] = fitConditionHierarchical(tbl_cond, fit_opts);
            if all(isfinite(pfit)) && any(all(isfinite(subject_params), 2))
                exitf = 1;
            else
                exitf = -99;
            end

        otherwise
            error('fitConditionByMethod:badFitMethod', ...
                'fit_method must be ''pooled'' or ''hierarchical_map''.');
    end
end
