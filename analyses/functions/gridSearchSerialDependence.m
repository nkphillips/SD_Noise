% gridSearchSerialDependence
%
% Grid search for serial dependence parameters
%
% Inputs:
%   init_params: initial parameters
%   fixed_params: [probe_offsets, responses, delta_thetas]
%   which_search: 'coarse' or 'fine'
%   p: parameter structure
%
% Outputs:
%   best_params: best parameters
%   best_nll: best negative log-likelihood

function [best_params, best_nll] = gridSearchSerialDependence(init_params, fixed_params, which_search, p)

    %% Define adaptive search steps based on data size
        
    switch which_search
        case 'coarse'
            % Fixed resolution per request
            search_steps = 10;
            amplitude_vals = linspace(p.sd_bounds(2,1), p.sd_bounds(1,1), search_steps);
            
            % Sample width uniformly in FWHM (i.e., uniformly in 1/w)
            fwhm_vals = linspace(1.6651/p.sd_bounds(1,2), 1.6651/p.sd_bounds(2,2), search_steps);
            w_vals = 1.6651 ./ fwhm_vals;
            baseline_vals = linspace(p.sd_bounds(2,3), p.sd_bounds(1,3), search_steps);
            
            % Define sigma_vals only if not using SSE
            if ~isfield(p, 'sd_objective') || ~strcmp(p.sd_objective, 'sse')
                sigma_vals = linspace(p.sd_bounds(2,4), p.sd_bounds(1,4), search_steps);
            else
                sigma_vals = []; % Not used for SSE, but needs to be defined
            end
        
        case 'fine'
            if init_params(1)/2 < p.sd_bounds(2,1)
                amplitude_lb = p.sd_bounds(2,1);
            else
                amplitude_lb = init_params(1)/2;
            end

            if init_params(1)*1.5 > p.sd_bounds(1,1)
                amplitude_ub = p.sd_bounds(1,1);
            else
                amplitude_ub = init_params(1)*1.5;
            end

            % Expand width range in FWHM space around the current init
            init_fwhm = 1.6651 / init_params(2);
            fwhm_lb = max(1.6651/p.sd_bounds(1,2), init_fwhm/2);
            fwhm_ub = min(1.6651/p.sd_bounds(2,2), init_fwhm*1.5);
            w_lb = 1.6651 / fwhm_ub;
            w_ub = 1.6651 / fwhm_lb;

            % Fixed resolution per request
            search_steps = 50;
            amplitude_vals = linspace(amplitude_lb, amplitude_ub, search_steps);
            fwhm_vals = linspace(fwhm_lb, fwhm_ub, search_steps);
            w_vals = 1.6651 ./ fwhm_vals;
            baseline_vals = linspace(p.sd_bounds(2,3), p.sd_bounds(1,3), search_steps);
            
            % Define sigma_vals only if not using SSE
            if ~isfield(p, 'sd_objective') || ~strcmp(p.sd_objective, 'sse')
                sigma_vals = linspace(p.sd_bounds(2,4), p.sd_bounds(1,4), search_steps);
            else
                sigma_vals = []; % Not used for SSE, but needs to be defined
            end
    end

    %% Vectorized evaluation
    
    % Create parameter grid based on objective
    if isfield(p, 'sd_objective') && strcmp(p.sd_objective, 'sse')
        % For SSE, we only need [amplitude, width, baseline] - no sigma
        [amplitude_grid, w_grid, baseline_grid] = ndgrid(amplitude_vals, w_vals, baseline_vals);
        params_matrix = [amplitude_grid(:), w_grid(:), baseline_grid(:)];
    else
        % For NLL, we need [amplitude, width, baseline, sigma]
        [amplitude_grid, w_grid, baseline_grid, sigma_grid] = ndgrid(amplitude_vals, w_vals, baseline_vals, sigma_vals);
        params_matrix = [amplitude_grid(:), w_grid(:), baseline_grid(:), sigma_grid(:)];
    end
    
    % Vectorized evaluation
    nll_vec = zeros(size(params_matrix, 1), 1);
    
    for i = 1:size(params_matrix, 1)
        amplitude = params_matrix(i, 1);
        w = params_matrix(i, 2);
        baseline = params_matrix(i, 3);
        
        % Calculate metric based on objective using unified entry point
        if isfield(p, 'sd_objective') && strcmp(p.sd_objective, 'sse')
            sd_params = [amplitude, w, baseline];
        else
            sigma = params_matrix(i, 4);
            sd_params = [amplitude, w, baseline, sigma];
        end
        current_metric = calcSerialDependenceFit(sd_params, fixed_params, p);
        nll_vec(i) = current_metric;
    end
    
    % Reshape back to grid
    nll_grid = reshape(nll_vec, size(amplitude_grid));

    %% Find best combination of parameters

    [best_nll, min_idx] = min(nll_grid(:));
    
    if isfield(p, 'sd_objective') && strcmp(p.sd_objective, 'sse')
        % For SSE, we only have [amplitude, width, baseline]
        [best_amplitude_idx, best_w_idx, best_baseline_idx] = ind2sub(size(nll_grid), min_idx);
        best_params = [amplitude_vals(best_amplitude_idx), w_vals(best_w_idx), baseline_vals(best_baseline_idx)];
    else
        % For NLL, we have [amplitude, width, baseline, sigma]
        [best_amplitude_idx, best_w_idx, best_baseline_idx, best_sigma_idx] = ind2sub(size(nll_grid), min_idx);
        best_params = [amplitude_vals(best_amplitude_idx), w_vals(best_w_idx), baseline_vals(best_baseline_idx), sigma_vals(best_sigma_idx)];
    end

end 