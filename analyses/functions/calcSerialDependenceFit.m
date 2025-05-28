function sse = calcSerialDependenceFit(free_params, fixed_params)

    delta_theta = fixed_params(:,1);
    measured_bias = fixed_params(:,2);

    amplitude = free_params(1);
    w = free_params(2);

    estimated_bias = gaussianPrime([amplitude, w], delta_theta);

    sse = calcSSE(measured_bias, estimated_bias);

end