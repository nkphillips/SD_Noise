function printBootstrapAdmissionSummary(sd_boot_cluster, num)
% printBootstrapAdmissionSummary
%
% Prints a one-block summary of how many subject-cluster bootstrap replicates
% were admitted vs discarded (and the breakdown by reason) after CI
% computation. Surfaces the admission rate so it's visible at runtime without
% having to inspect the saved estimates.

    if nargin < 2 || isempty(num)
        num = struct();
    end

    if isempty(sd_boot_cluster) || ~isfield(sd_boot_cluster, 'discard_summary') || ...
            isempty(sd_boot_cluster.discard_summary)
        return
    end

    summary = sd_boot_cluster.discard_summary;
    B = summary.B_total;

    if isfield(num, 'levels') && isfield(num, 'conds')
        cells = num.levels * num.levels * num.conds;
    else
        cells = numel(summary.admitted_per_cell);
    end

    total_fits = B * cells;
    admitted   = summary.global_admitted;
    discarded  = summary.global_discarded;
    if total_fits > 0
        pct_adm = 100 * admitted / total_fits;
        pct_dis = 100 * discarded / total_fits;
    else
        pct_adm = NaN; pct_dis = NaN;
    end

    by_reason = summary.by_reason;
    disp(' ');
    disp('=== BOOTSTRAP CI ADMISSION SUMMARY ===');
    disp(['B_total: ' num2str(B) ' x ' num2str(cells) ' cells = ' num2str(total_fits) ' fits']);
    disp(['Admitted:  ' num2str(admitted)  ' (' num2str(pct_adm, '%.1f') '%)']);
    disp(['Discarded: ' num2str(discarded) ' (' num2str(pct_dis, '%.1f') '%) - ' ...
        'min_windows: '  num2str(by_reason.min_windows)  ', ' ...
        'missing_lobe: ' num2str(by_reason.missing_lobe) ', ' ...
        'failed_fit: '   num2str(by_reason.failed_fit)   ', ' ...
        'at_bound: '     num2str(by_reason.at_bound)]);
    disp('======================================');
end
