function estimates = loadEstimates(analysis_date, n_back_list, opts)
% loadEstimates
%
% Shared loader for the per-n-back estimates produced by
% SD_Noise_Analyses_And_Figures.m. Returns a struct keyed by n-back
% containing rb, sd, rb_ci, sd_ci, sd_ci_cluster, perf_ci, meta and a convenience
% delta_theta_centers field pulled from meta when available.
%
% Usage:
%   est = loadEstimates('11.24.2025', [1 2 3]);
%   est.n1.rb.all.params_est(...)
%
% Optional opts:
%   opts.objective   - 'sse' (default) or 'nll'; filters the file suffix.
%   opts.which       - 'latest' (default) picks the newest matching file
%                      per n-back; 'all' returns a struct array of all
%                      matches per n-back.
%   opts.root        - override root estimates dir; defaults to
%                      ../estimates relative to this file.

    if nargin < 3, opts = struct(); end
    if ~isfield(opts, 'objective') || isempty(opts.objective)
        opts.objective = 'sse';
    end
    if ~isfield(opts, 'which') || isempty(opts.which)
        opts.which = 'latest';
    end

    this_dir = fileparts(mfilename('fullpath'));
    if ~isfield(opts, 'root') || isempty(opts.root)
        opts.root = fullfile(this_dir, '..', 'estimates');
    end

    estimates = struct();
    estimates.analysis_date = analysis_date;
    estimates.objective     = opts.objective;
    estimates.root          = opts.root;

    for i = 1:numel(n_back_list)
        n_back = n_back_list(i);
        key = sprintf('n%d', n_back);
        folder = fullfile(opts.root, analysis_date, sprintf('%d_back', n_back));

        if ~exist(folder, 'dir')
            warning('loadEstimates:missingFolder', ...
                'Estimates folder missing: %s', folder);
            estimates.(key) = [];
            continue
        end

        pattern = fullfile(folder, ...
            sprintf('SD_Noise_Estimates_*_%s.mat', opts.objective));
        files = dir(pattern);

        if isempty(files)
            warning('loadEstimates:noFiles', ...
                'No %s estimates found under %s', opts.objective, folder);
            estimates.(key) = [];
            continue
        end

        % Newest first
        [~, ord] = sort([files.datenum], 'descend');
        files = files(ord);

        if strcmp(opts.which, 'latest')
            estimates.(key) = load_one(fullfile(files(1).folder, files(1).name));
        else
            arr = repmat(load_one(fullfile(files(1).folder, files(1).name)), 1, numel(files));
            for k = 2:numel(files)
                arr(k) = load_one(fullfile(files(k).folder, files(k).name));
            end
            estimates.(key) = arr;
        end
    end
end

function out = load_one(full_path)
    S = load(full_path);

    out.file = full_path;

    fields = {'rb','sd','rb_ci','sd_ci','sd_ci_cluster','sd_boot_cluster','perf_ci'};
    for f = 1:numel(fields)
        if isfield(S, fields{f})
            out.(fields{f}) = S.(fields{f});
        else
            out.(fields{f}) = [];
        end
    end

    if isfield(S, 'meta')
        out.meta = S.meta;
    else
        out.meta = struct();
    end

    % Convenience: surface delta_theta_centers at the top level when present.
    if isfield(out.meta, 'delta_theta_centers')
        out.delta_theta_centers = out.meta.delta_theta_centers;
    else
        out.delta_theta_centers = [];
    end
end
