function y = calcDoG(x, params)
    % calcDoG
    %
    % Computes the first derivative of gaussian function with baseline term
    % x = difference between previous and current orientation (delta theta)
    % params = [amplitude, width, baseline] or matrix of [amplitude, width, baseline] pairs
    %
    % Outputs:
    %   y = predicted bias for each trial (vector or matrix)
    
    c = sqrt(2)/exp(-0.5);
    
    if size(params, 1) == 1
        % Single parameter case (current usage)
        A = params(1);
        w = params(2);
        b = params(3); % baseline term
        y = (x * A * w * c) .* exp(-(w * x).^2) + b;
    else
        % Multiple parameter case (for vectorization)
        y = zeros(length(x), size(params, 1));
        for i = 1:size(params, 1)
            A = params(i, 1);
            w = params(i, 2);
            b = params(i, 3); % baseline term
            y(:, i) = (x * A * w * c) .* exp(-(w * x).^2) + b;
        end
    end
end 