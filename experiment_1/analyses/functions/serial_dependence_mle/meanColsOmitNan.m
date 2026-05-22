function y = meanColsOmitNan(A)
% meanColsOmitNan  Row vector of column means, ignoring NaNs (compat with older MATLAB).

    [~, nCol] = size(A);
    y = nan(1, nCol);
    for j = 1:nCol
        v = A(:, j);
        v = v(~isnan(v));
        if ~isempty(v)
            y(j) = mean(v);
        end
    end

end
