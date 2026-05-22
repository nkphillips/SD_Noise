function y = medianColsOmitNan(A)
% medianColsOmitNan  Row vector of column medians, ignoring NaNs (compat with older MATLAB).

    [~, nCol] = size(A);
    y = nan(1, nCol);
    for j = 1:nCol
        v = A(:, j);
        v = v(~isnan(v));
        if ~isempty(v)
            y(j) = median(v);
        end
    end

end
