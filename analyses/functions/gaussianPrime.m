% Computes the first derivative of gaussian function
% x = difference between previous and current stimulus
% amplitude = amplitude of DoG curve
% c scales the amplitude and is a constant
% w = the inverse of the curve width

function y = gaussianPrime(params, x)

    c = sqrt(2)/exp(-0.5);

    amplitude = params(1);
    w = params(2);

    y = x * amplitude * w * c .* exp(-((w * x).^2));
 
end


