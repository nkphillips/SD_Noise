function specs = buildStandardContrasts(varargin)
% buildStandardContrasts  Build a curated set of within-manipulation contrasts
% on the 3 (prev) x 3 (curr) design for both manipulations. By default returns
% contrasts on amplitude A only; pass 'params', {'A','w','sigma','beta'} (any
% subset) to compute on more than A.
%
% Cell indexing convention (conditionIndexUnbinned): c = (m-1)*9 + (prev-1)*3 + curr.
% For each manipulation m in {1=contrast, 2=precision}:
%   prev=1 corresponds to the highest-information level (90% contrast or 2 deg precision).
%   prev=3 corresponds to the lowest-information level (25% contrast or 80 deg precision).
%
% Returned contrasts (defaults, on amplitude A):
%   <manip>.prev_main      -- A_(prev=high) - A_(prev=low), marginalized over curr
%   <manip>.curr_main      -- A_(curr=high) - A_(curr=low), marginalized over prev
%   <manip>.diagonal_HH_LL -- A_(prev=high,curr=high) - A_(prev=low,curr=low)
%   <manip>.diagonal_HL_LH -- A_(prev=high,curr=low) - A_(prev=low,curr=high)
%   <manip>.cell_grand     -- mean of A across all 9 cells in this manipulation
%                             (used as a "is the manipulation-marginal A different from 0" check)
%
% Optional name-value pair:
%   'params' -- cell array of parameter names from {'A','w','sigma','beta'}; default {'A'}.

    ip = inputParser;
    addParameter(ip, 'params', {'A'}, @iscell);
    parse(ip, varargin{:});
    pnames_to_use = ip.Results.params;

    pname_to_idx = struct('A', 1, 'w', 2, 'sigma', 3, 'beta', 4);

    cidx = @(m, prev, curr) (m - 1) * 9 + (prev - 1) * 3 + curr;

    specs = struct('name', {}, 'weights', {}, 'param_idx', {});

    manip_labels = {'contrast', 'precision'};

    for ip_loop = 1:numel(pnames_to_use)
        pname = pnames_to_use{ip_loop};
        if ~isfield(pname_to_idx, pname)
            error('buildStandardContrasts:badParam', ...
                'Unknown parameter name: %s. Use one of A, w, sigma, beta.', pname);
        end
        k = pname_to_idx.(pname);

        for m = 1:2
            mname = manip_labels{m};

            % --- Main effect of prev: prev=1 (high) - prev=3 (low), marginalized over curr ---
            w = zeros(1, 18);
            w(cidx(m, 1, 1):cidx(m, 1, 3)) = 1/3;
            w(cidx(m, 3, 1):cidx(m, 3, 3)) = -1/3;
            specs(end+1) = struct('name', sprintf('%s.%s.prev_main(L1-L3)', pname, mname), ...
                                  'weights', w, 'param_idx', k); %#ok<AGROW>

            % --- Main effect of curr: curr=1 (high) - curr=3 (low), marginalized over prev ---
            w = zeros(1, 18);
            w([cidx(m,1,1), cidx(m,2,1), cidx(m,3,1)]) = 1/3;
            w([cidx(m,1,3), cidx(m,2,3), cidx(m,3,3)]) = -1/3;
            specs(end+1) = struct('name', sprintf('%s.%s.curr_main(L1-L3)', pname, mname), ...
                                  'weights', w, 'param_idx', k); %#ok<AGROW>

            % --- Diagonal: high->high vs low->low ---
            w = zeros(1, 18);
            w(cidx(m, 1, 1)) = 1;
            w(cidx(m, 3, 3)) = -1;
            specs(end+1) = struct('name', sprintf('%s.%s.diag(L1L1-L3L3)', pname, mname), ...
                                  'weights', w, 'param_idx', k); %#ok<AGROW>

            % --- Off-diagonal: high->low vs low->high ---
            w = zeros(1, 18);
            w(cidx(m, 1, 3)) = 1;
            w(cidx(m, 3, 1)) = -1;
            specs(end+1) = struct('name', sprintf('%s.%s.offdiag(L1L3-L3L1)', pname, mname), ...
                                  'weights', w, 'param_idx', k); %#ok<AGROW>

            % --- Cell grand-mean: average of all 9 cells in this manipulation ---
            w = zeros(1, 18);
            for prev = 1:3
                for curr = 1:3
                    w(cidx(m, prev, curr)) = 1/9;
                end
            end
            specs(end+1) = struct('name', sprintf('%s.%s.cell_grand', pname, mname), ...
                                  'weights', w, 'param_idx', k); %#ok<AGROW>
        end
    end
end
