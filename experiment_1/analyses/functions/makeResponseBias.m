function rb = makeResponseBias(num)

    %% Initialize structure

    rb = struct('ind',[],'grp',[],'all',[]);

    %% Subject

    struct_size = cell(1,num.subjs);
    rb.ind = struct('start_params', struct_size, 'start_nll', struct_size, 'params_est', struct_size, 'nll', struct_size, 'null_nll', struct_size, 'better_than_null', struct_size, 'exit_flag', struct_size);

    fieldname_list = fieldnames(rb.ind);

    for subj = 1:num.subjs

        for i_field = 1:numel(fieldname_list)

            if contains(fieldname_list{i_field}, {'start_params','params_est'})
                rb.ind(subj).(fieldname_list{i_field}) = nan(num.levels, num.levels, num.conds, num.delta_theta_windows, 2);
            elseif contains(fieldname_list{i_field}, {'start_nll','nll','exit_flag','null_nll','better_than_null'})
                rb.ind(subj).(fieldname_list{i_field}) = nan(num.levels, num.levels, num.conds, num.delta_theta_windows);
            end

        end

    end 

    %% Super subject

    rb.all.start_params = nan(num.levels, num.levels, num.conds, num.delta_theta_windows, 2);
    rb.all.start_nll = nan(num.levels, num.levels, num.conds, num.delta_theta_windows);
    rb.all.params_est = nan(num.levels, num.levels, num.conds, num.delta_theta_windows, 2);
    rb.all.nll = nan(num.levels, num.levels, num.conds, num.delta_theta_windows);
    rb.all.null_nll = nan(num.levels, num.levels, num.conds, num.delta_theta_windows);
    rb.all.better_than_null = nan(num.levels, num.levels, num.conds, num.delta_theta_windows);
    rb.all.exit_flag = nan(num.levels, num.levels, num.conds, num.delta_theta_windows);

end