% Computes the first derivative of gaussian function
% x = difference between previous and current stimulus
% amplitude = amplitude of DoG curve
% c scales the amplitude and is a constant
% width_param = the width of the curve (e.g., standard deviation of the underlying Gaussian)

function y = gaussianPrime(params, x)

    c = sqrt(2)/exp(-0.5);

    amplitude = params(1);
    width_param = params(2); % Direct width parameter
    
    if width_param == 0 % Avoid division by zero; treat as infinitely narrow
        w = 1 / 0.01;; 
    else
        w = 1 / width_param; % Inverse of the curve width
    end

    y = x .* amplitude .* w .* c .* exp(-((w .* x).^2));
 
end


