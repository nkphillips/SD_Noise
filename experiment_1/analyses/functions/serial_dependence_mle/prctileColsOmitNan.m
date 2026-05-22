function y = prctileColsOmitNan(A, p)
% prctileColsOmitNan  Row vector of column percentiles, ignoring NaNs (compat with older MATLAB).

    [~, nCol] = size(A);
    y = nan(1, nCol);
    for j = 1:nCol
        v = A(:, j);
        v = v(~isnan(v));
        if ~isempty(v)
            y(j) = prctile(v, p);
        end
    end

end
