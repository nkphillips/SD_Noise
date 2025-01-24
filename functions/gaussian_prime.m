% Computes the first derivative of gaussian function
% delta = difference between previous and current stimulus
% amplitude = amplitude of DoG curve
% c scales the amplitude and is a constant
% w = the inverse of the curve width

function y = gaussian_prime(params,delta)

    c = sqrt(2)/exp(-0.5);

    amplitude = params(1);
    w = params(2);

    y = delta .* amplitude .* w * c .* exp(-((w.*delta).^2));

end


