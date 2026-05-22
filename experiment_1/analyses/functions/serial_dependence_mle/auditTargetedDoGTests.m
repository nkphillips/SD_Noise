function audit = auditTargetedDoGTests(targeted_tests, varargin)
% auditTargetedDoGTests  Report impossible or fragile targeted DoG test rows.
%
% Usage:
%   audit = auditTargetedDoGTests('targeted_dog_hypothesis_tests_all_nback.csv');
%   audit = auditTargetedDoGTests(tests_table, 'summary_root', output_root);
%
% This helper is read-only. It does not modify CSVs or figures.

    ip = inputParser;
    addParameter(ip, 'summary_root', '', @(x) ischar(x) || isstring(x));
    addParameter(ip, 'max_abs_z0', 2, @(x) isnumeric(x) && isscalar(x));
    addParameter(ip, 'min_admit_fraction', 0.20, @(x) isnumeric(x) && isscalar(x));
    parse(ip, varargin{:});

    tests = local_readTable(targeted_tests);
    audit = struct();
    audit.tests = local_auditTests(tests, ip.Results.max_abs_z0, ip.Results.min_admit_fraction);
    audit.summary = local_auditSummaries(char(ip.Results.summary_root));

    local_printAudit(audit);
end

function tbl = local_readTable(src)
    if istable(src)
        tbl = src;
    elseif ischar(src) || isstring(src)
        tbl = readtable(char(src));
    else
        error('auditTargetedDoGTests:badInput', ...
            'targeted_tests must be a table or path to a targeted tests CSV.');
    end
end

function out = local_auditTests(tests, max_abs_z0, min_admit_fraction)
    n = height(tests);
    estimate_outside_pc = false(n, 1);
    estimate_outside_bca = false(n, 1);
    bca_pc_conflict = false(n, 1);
    low_admit_fraction = false(n, 1);
    extreme_z0 = false(n, 1);

    for i = 1:n
        estimate_outside_pc(i) = ~local_inInterval(tests.estimate(i), tests.pc_lo(i), tests.pc_hi(i));
        estimate_outside_bca(i) = ~local_inInterval(tests.estimate(i), tests.bca_lo(i), tests.bca_hi(i));
        bca_pc_conflict(i) = local_ciConflict(tests.pc_lo(i), tests.pc_hi(i), tests.bca_lo(i), tests.bca_hi(i));
        if ismember('admit_fraction', tests.Properties.VariableNames)
            low_admit_fraction(i) = tests.admit_fraction(i) < min_admit_fraction;
        elseif ismember('n_admit', tests.Properties.VariableNames) && ismember('B', tests.Properties.VariableNames)
            low_admit_fraction(i) = tests.n_admit(i) ./ tests.B(i) < min_admit_fraction;
        end
        if ismember('z0', tests.Properties.VariableNames)
            extreme_z0(i) = isfinite(tests.z0(i)) && abs(tests.z0(i)) > max_abs_z0;
        end
    end

    flagged = estimate_outside_pc | estimate_outside_bca | bca_pc_conflict | ...
        low_admit_fraction | extreme_z0;
    out = tests(flagged, :);
    if ~isempty(out) && height(out) > 0
        out.audit_estimate_outside_pc = estimate_outside_pc(flagged);
        out.audit_estimate_outside_bca = estimate_outside_bca(flagged);
        out.audit_bca_pc_conflict = bca_pc_conflict(flagged);
        out.audit_low_admit_fraction = low_admit_fraction(flagged);
        out.audit_extreme_z0 = extreme_z0(flagged);
    end
end

function out = local_auditSummaries(summary_root)
    out = table();
    if isempty(summary_root) || ~exist(summary_root, 'dir')
        return
    end

    files = dir(fullfile(summary_root, '**', 'summary_table.csv'));
    file_col = strings(0, 1);
    row_col = zeros(0, 1);
    variable_col = strings(0, 1);
    lo_col = zeros(0, 1);
    hi_col = zeros(0, 1);

    pairs = {'A_ci_lo', 'A_ci_hi', 'A'; ...
        'w_ci_lo', 'w_ci_hi', 'w'; ...
        'sigma_ci_lo', 'sigma_ci_hi', 'sigma'; ...
        'beta_ci_lo', 'beta_ci_hi', 'beta'; ...
        'fwhm_ci_lo', 'fwhm_ci_hi', 'FWHM'};

    for i_file = 1:numel(files)
        path_i = fullfile(files(i_file).folder, files(i_file).name);
        tbl = readtable(path_i);
        for i_pair = 1:size(pairs, 1)
            lo_name = pairs{i_pair, 1};
            hi_name = pairs{i_pair, 2};
            if ~ismember(lo_name, tbl.Properties.VariableNames) || ~ismember(hi_name, tbl.Properties.VariableNames)
                continue
            end
            bad = isfinite(tbl.(lo_name)) & isfinite(tbl.(hi_name)) & tbl.(lo_name) > tbl.(hi_name);
            bad_idx = find(bad);
            for j = 1:numel(bad_idx)
                r = bad_idx(j);
                file_col(end+1, 1) = string(path_i); %#ok<AGROW>
                row_col(end+1, 1) = r; %#ok<AGROW>
                variable_col(end+1, 1) = string(pairs{i_pair, 3}); %#ok<AGROW>
                lo_col(end+1, 1) = tbl.(lo_name)(r); %#ok<AGROW>
                hi_col(end+1, 1) = tbl.(hi_name)(r); %#ok<AGROW>
            end
        end
    end

    out = table(file_col, row_col, variable_col, lo_col, hi_col, ...
        'VariableNames', {'file', 'row', 'variable', 'lo', 'hi'});
end

function local_printAudit(audit)
    fprintf('Targeted DoG test audit: %d flagged rows.\n', height(audit.tests));
    if ~isempty(audit.tests) && height(audit.tests) > 0
        show_cols = {'n_back', 'name', 'estimate', 'bca_lo', 'bca_hi', 'pc_lo', 'pc_hi'};
        show_cols = show_cols(ismember(show_cols, audit.tests.Properties.VariableNames));
        disp(audit.tests(:, show_cols));
    end

    fprintf('Summary-table CI audit: %d reversed intervals.\n', height(audit.summary));
    if ~isempty(audit.summary) && height(audit.summary) > 0
        disp(audit.summary);
    end
end

function tf = local_inInterval(x, lo, hi)
    if ~isfinite(x) || ~isfinite(lo) || ~isfinite(hi)
        tf = false;
        return
    end
    tf = x >= min(lo, hi) && x <= max(lo, hi);
end

function tf = local_excludesZero(lo, hi)
    if ~isfinite(lo) || ~isfinite(hi)
        tf = false;
        return
    end
    tf = min(lo, hi) > 0 || max(lo, hi) < 0;
end

function tf = local_ciConflict(pc_lo, pc_hi, bca_lo, bca_hi)
    pc_sig = local_excludesZero(pc_lo, pc_hi);
    bca_sig = local_excludesZero(bca_lo, bca_hi);
    if pc_sig ~= bca_sig
        tf = true;
    elseif pc_sig && bca_sig
        tf = sign(mean([pc_lo, pc_hi])) ~= sign(mean([bca_lo, bca_hi]));
    else
        tf = false;
    end
end
