function specs = buildStandardContrasts(varargin)
% buildStandardContrasts  Within-manipulation contrasts on the 3 (prev) x 3 (curr)
% factorial design for both manipulations. Reports ALL level-pair comparisons:
% three pairwise main effects per axis, three pairwise diagonal comparisons,
% three pairwise off-diagonal comparisons, and the manipulation-grand mean.
% Defaults to amplitude only; pass 'params', {'A','w','sigma','beta'} for more.
%
% Cell indexing convention (conditionIndexUnbinned): c = (m-1)*9 + (prev-1)*3 + curr.
% Levels are ordered so that L1 = highest-information level (90% contrast or 2 deg precision)
% and L3 = lowest-information level (25% contrast or 80 deg precision).
%
% Returned contrasts per parameter, per manipulation (13 total):
%   <param>.<manip>.prev_main(L1-L2), prev_main(L1-L3), prev_main(L2-L3)
%   <param>.<manip>.curr_main(L1-L2), curr_main(L1-L3), curr_main(L2-L3)
%   <param>.<manip>.diag(L1L1-L2L2), diag(L1L1-L3L3), diag(L2L2-L3L3)
%   <param>.<manip>.offdiag(L1L2-L2L1), offdiag(L1L3-L3L1), offdiag(L2L3-L3L2)
%   <param>.<manip>.cell_grand
%
% Total per parameter: 13 contrasts/manip x 2 manips = 26.
% With all 4 parameters: 104 contrasts.
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

    % All pairwise level comparisons for a 3-level design
    level_pairs = {[1 2], [1 3], [2 3]};
    pair_labels = {'L1-L2', 'L1-L3', 'L2-L3'};

    for ip_loop = 1:numel(pnames_to_use)
        pname = pnames_to_use{ip_loop};
        if ~isfield(pname_to_idx, pname)
            error('buildStandardContrasts:badParam', ...
                'Unknown parameter name: %s. Use one of A, w, sigma, beta.', pname);
        end
        k = pname_to_idx.(pname);

        for m = 1:2
            mname = manip_labels{m};

            % --- Main effect of prev (all three pairwise contrasts, marginalized over curr) ---
            for ipair = 1:numel(level_pairs)
                la = level_pairs{ipair}(1);
                lb = level_pairs{ipair}(2);
                w = zeros(1, 18);
                for c = 1:3
                    w(cidx(m, la, c)) = 1/3;
                    w(cidx(m, lb, c)) = -1/3;
                end
                specs(end+1) = struct('name', sprintf('%s.%s.prev_main(%s)', pname, mname, pair_labels{ipair}), ...
                                      'weights', w, 'param_idx', k); %#ok<AGROW>
            end

            % --- Main effect of curr (all three pairwise contrasts, marginalized over prev) ---
            for ipair = 1:numel(level_pairs)
                la = level_pairs{ipair}(1);
                lb = level_pairs{ipair}(2);
                w = zeros(1, 18);
                for prev = 1:3
                    w(cidx(m, prev, la)) = 1/3;
                    w(cidx(m, prev, lb)) = -1/3;
                end
                specs(end+1) = struct('name', sprintf('%s.%s.curr_main(%s)', pname, mname, pair_labels{ipair}), ...
                                      'weights', w, 'param_idx', k); %#ok<AGROW>
            end

            % --- Diagonal comparisons: cell(L_a, L_a) - cell(L_b, L_b) for each pair ---
            for ipair = 1:numel(level_pairs)
                la = level_pairs{ipair}(1);
                lb = level_pairs{ipair}(2);
                w = zeros(1, 18);
                w(cidx(m, la, la)) =  1;
                w(cidx(m, lb, lb)) = -1;
                specs(end+1) = struct('name', sprintf('%s.%s.diag(L%dL%d-L%dL%d)', pname, mname, la, la, lb, lb), ...
                                      'weights', w, 'param_idx', k); %#ok<AGROW>
            end

            % --- Off-diagonal comparisons: cell(L_a, L_b) - cell(L_b, L_a) for each pair ---
            for ipair = 1:numel(level_pairs)
                la = level_pairs{ipair}(1);
                lb = level_pairs{ipair}(2);
                w = zeros(1, 18);
                w(cidx(m, la, lb)) =  1;
                w(cidx(m, lb, la)) = -1;
                specs(end+1) = struct('name', sprintf('%s.%s.offdiag(L%dL%d-L%dL%d)', pname, mname, la, lb, lb, la), ...
                                      'weights', w, 'param_idx', k); %#ok<AGROW>
            end

            % --- Cell-grand mean: average across all 9 cells in this manipulation ---
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
