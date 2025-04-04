function sse = calc_fit(free_params, fixed_params)

    x = fixed_params{1}(:,1);
    p_CW = fixed_params{1}(:,2);
    guess_rate = fixed_params{2};

    mu = free_params(1);
    sigma = free_params(2);

    p_CW_est = calc_pCW(x, mu, sigma, guess_rate);

    sse = sum((p_CW_est - p_CW).^2); % sum of squared errors

end