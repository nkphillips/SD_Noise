% Computes the first derivative of gaussian function
% delta = difference between previous and current stimulus
% amplitude = amplitude of DoG curve
% c scales the amplitude and is a constant
% w = the inverse of the curve width

function y = gaussian_prime_split(params, delta)

    c = sqrt(2)/exp(-0.5);

    amplitude_neg = params(1);
    amplitude_pos = params(2);
    
    w_neg = params(3);
    w_pos = params(4);

    y = zeros(size(delta)); 

    y(delta < 0) = delta(delta < 0) .* amplitude_neg .* w_neg * c .* exp(-((w_neg .* delta(delta < 0)).^2));
    y(delta >= 0) = delta(delta >= 0) .* amplitude_pos .* w_pos * c .* exp(-((w_pos .* delta(delta >= 0)).^2));

end


