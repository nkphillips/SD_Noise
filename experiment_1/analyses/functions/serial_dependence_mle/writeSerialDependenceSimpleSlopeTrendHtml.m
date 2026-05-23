function html_path = writeSerialDependenceSimpleSlopeTrendHtml(fig_dir, simple_tests, ci_method)
% writeSerialDependenceSimpleSlopeTrendHtml  Write browser-readable simple-slope report.

    if nargin < 3 || isempty(ci_method)
        ci_method = local_tableCIMethod(simple_tests);
    end
    if ~exist(fig_dir, 'dir')
        mkdir(fig_dir);
    end
    html_path = fullfile(fig_dir, 'simple_slope_trend_results.html');

    fid = fopen(html_path, 'w');
    if fid < 0
        error('writeSerialDependenceSimpleSlopeTrendHtml:openFailed', ...
            'Could not open HTML report for writing: %s', html_path);
    end
    cleaner = onCleanup(@() fclose(fid));

    fprintf(fid, '<!doctype html><html><head><meta charset="utf-8">\n');
    fprintf(fid, '<title>Simple-slope trend results</title>\n');
    fprintf(fid, '<style>\n');
    fprintf(fid, 'body{font-family:Helvetica,Arial,sans-serif;margin:28px;color:#222;line-height:1.35} ');
    fprintf(fid, 'h1{font-size:24px;margin-bottom:4px} h2{font-size:18px;margin-top:28px} ');
    fprintf(fid, '.note{color:#555;max-width:980px}.links a{display:inline-block;margin:0 12px 8px 0} ');
    fprintf(fid, 'table{border-collapse:collapse;width:100%%;font-size:12px;margin-top:12px} ');
    fprintf(fid, 'th,td{border-bottom:1px solid #ddd;padding:5px 7px;text-align:left;white-space:nowrap} ');
    fprintf(fid, 'th{background:#f3f3f3;position:sticky;top:0}.supported{background:#ecf7ec}.fragile{color:#777}.num{text-align:right} ');
    fprintf(fid, '</style></head><body>\n');
    fprintf(fid, '<h1>Simple-slope trend results</h1>\n');
    fprintf(fid, '<p class="note">Each row tests a three-point linear trend within one subplot line. Levels are coded L1=-1, L2=0, L3=+1, so each estimate is a slope per level step. CIs are 95%% %s intervals. Supported rows require the slope CI to exclude zero and the point estimates to be monotonic across L1, L2, and L3.</p>\n', local_escape(local_ciLabel(ci_method)));
    fprintf(fid, '<div class="links">\n');
    local_link(fid, 'simple_slope_trends_amplitude_by_past_level.pdf', 'Amplitude by past level');
    local_link(fid, 'simple_slope_trends_amplitude_by_current_level.pdf', 'Amplitude by current level');
    local_link(fid, 'simple_slope_trends_fwhm_by_past_level.pdf', 'FWHM by past level');
    local_link(fid, 'simple_slope_trends_fwhm_by_current_level.pdf', 'FWHM by current level');
    local_link(fid, 'simple_slope_trend_tests.csv', 'CSV');
    fprintf(fid, '</div>\n');

    if isempty(simple_tests) || height(simple_tests) == 0
        fprintf(fid, '<p>No simple-slope rows available.</p></body></html>\n');
        return
    end

    supported = simple_tests(logical(simple_tests.supports_effect), :);
    fprintf(fid, '<h2>Supported simple slopes</h2>\n');
    if isempty(supported) || height(supported) == 0
        fprintf(fid, '<p class="note">No rows passed the support flags. Inspect the full table for descriptive slopes and fragile rows.</p>\n');
    else
        local_table(fid, supported);
    end

    fprintf(fid, '<h2>All simple slopes</h2>\n');
    local_table(fid, simple_tests);
    fprintf(fid, '</body></html>\n');
end

function local_link(fid, href, label)
    fprintf(fid, '<a href="%s">%s</a>\n', local_escape(href), local_escape(label));
end

function local_table(fid, T)
    cols = {'n_back', 'parameter', 'manipulation', 'fixed_axis', 'fixed_level_label', ...
        'slope_axis', 'estimate', 'bca_lo', 'bca_hi', 'middle_deviation', ...
        'point_monotonic', 'p_bca_label', ...
        'p_fdr_bh_label', 'admit_fraction', 'valid_for_inference'};
    cols = cols(ismember(cols, T.Properties.VariableNames));
    fprintf(fid, '<table><thead><tr>');
    for i = 1:numel(cols)
        fprintf(fid, '<th>%s</th>', local_escape(cols{i}));
    end
    fprintf(fid, '</tr></thead><tbody>\n');

    for r = 1:height(T)
        cls = '';
        if ismember('supports_effect', T.Properties.VariableNames) && logical(T.supports_effect(r))
            cls = ' class="supported"';
        elseif ismember('valid_for_inference', T.Properties.VariableNames) && ~logical(T.valid_for_inference(r))
            cls = ' class="fragile"';
        end
        fprintf(fid, '<tr%s>', cls);
        for i = 1:numel(cols)
            val = T.(cols{i})(r);
            if isnumeric(val) || islogical(val)
                if isfinite(double(val))
                    if contains(cols{i}, 'fraction')
                        txt = sprintf('%.3f', double(val));
                    elseif contains(cols{i}, 'n_back') || contains(cols{i}, 'fixed_level')
                        txt = sprintf('%.0f', double(val));
                    else
                        txt = sprintf('%.3f', double(val));
                    end
                else
                    txt = 'NA';
                end
                fprintf(fid, '<td class="num">%s</td>', txt);
            else
                fprintf(fid, '<td>%s</td>', local_escape(char(string(val))));
            end
        end
        fprintf(fid, '</tr>\n');
    end
    fprintf(fid, '</tbody></table>\n');
end

function ci_method = local_tableCIMethod(T)
    if ~isempty(T) && height(T) > 0 && ismember('ci_method', T.Properties.VariableNames)
        ci_method = char(T.ci_method(1));
    else
        ci_method = 'bootstrap';
    end
end

function s = local_ciLabel(ci_method)
    ci_method = lower(char(ci_method));
    switch ci_method
        case 'bca'
            s = 'BCa';
        case 'percentile'
            s = 'percentile';
        otherwise
            s = char(ci_method);
    end
end

function s = local_escape(s)
    s = char(string(s));
    s = strrep(s, '&', '&amp;');
    s = strrep(s, '<', '&lt;');
    s = strrep(s, '>', '&gt;');
    s = strrep(s, '"', '&quot;');
    s = strrep(s, '''', '&#39;');
end
