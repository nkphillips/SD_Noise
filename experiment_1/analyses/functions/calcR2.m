function R2 = calcR2(y_observed, y_fitted)
% calcR2 - Robust R-squared that tolerates row/column shapes and NaNs

% Ensure column vectors and matched lengths
y_observed = y_observed(:);
y_fitted = y_fitted(:);
n = min(numel(y_observed), numel(y_fitted));
y_observed = y_observed(1:n);
y_fitted = y_fitted(1:n);

% Remove NaNs pairwise
valid_idx = ~(isnan(y_observed) | isnan(y_fitted));
y_observed = y_observed(valid_idx);
y_fitted = y_fitted(valid_idx);

if isempty(y_observed)
    R2 = NaN;
    return;
end

% Calculate the mean of the observed data
y_mean = mean(y_observed);

% Compute the total sum of squares (proportional to variance)
SS_tot = sum((y_observed - y_mean).^2, 'omitnan');

% Compute the residual sum of squares
SS_res = sum((y_observed - y_fitted).^2, 'omitnan');

% Compute R-squared
R2 = 1 - (SS_res / SS_tot);

end
