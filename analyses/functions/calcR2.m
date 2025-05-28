function R2 = calcR2(y_observed, y_fitted)
   
% Calculate the mean of the observed data
y_mean = mean(y_observed);

% Compute the total sum of squares (proportional to variance)
SS_tot = sum((y_observed - y_mean).^2);

% Compute the residual sum of squares
SS_res = sum((y_observed - y_fitted).^2);

% Compute R-squared
R2 = 1 - (SS_res / SS_tot);

end
