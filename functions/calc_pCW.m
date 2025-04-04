function p_CW = calc_pCW(x, mu, sigma, guess_rate)

    p_CW = (1 - guess_rate) * normcdf(x, mu, sigma) + (0.5 * guess_rate);

end