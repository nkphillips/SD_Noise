function serial_dependence = estimateSerialDependence(init_params, fixed_params, p)

    %% Grid search for starting parameters

    start_params = gridSearchSerialDependence(init_params, fixed_params, 'coarse', p.serial_dependence_bounds);
    [start_params, start_sse] = gridSearchSerialDependence(start_params, fixed_params, 'fine', p.serial_dependence_bounds);
    free_params = start_params;

    %% Define the model function handle

    serial_dependence_model = @(free_params) calcSerialDependenceFit(free_params, fixed_params);
    
    %% Fit the model

    [params_est, sse, exit_flag] = fmincon(serial_dependence_model, free_params, [], [], [], [], p.serial_dependence_bounds(2,:), p.serial_dependence_bounds(1,:), [], p.fmincon_options);

    %% Evaluate the model

    estimated_bias = gaussianPrime(params_est, fixed_params(:,1));
    r2 = calcR2(fixed_params(:,2), estimated_bias);

    %% Compile serial dependence

    serial_dependence.start_params = start_params;
    serial_dependence.start_sse = start_sse;
    serial_dependence.params_est = params_est;
    serial_dependence.sse = sse;
    serial_dependence.r2 = r2;
    serial_dependence.exit_flag = exit_flag;

end