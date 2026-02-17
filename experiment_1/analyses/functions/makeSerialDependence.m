function sd = makeSerialDependence(num)

    %% Initialize structure

    sd = struct('ind',[],'grp',[],'all',[]);

    struct_size = cell(1,num.subjs);
    sd.ind = struct('start_params', struct_size, 'start_nll', struct_size, 'params_est', struct_size, 'nll', struct_size, 'exit_flag', struct_size);

    fieldname_list = fieldnames(sd.ind);

    %% Subject

    for subj = 1:num.subjs

        for i_field = 1:numel(fieldname_list)

            if contains(fieldname_list{i_field}, {'start_params','params_est'})
                sd.ind(subj).(fieldname_list{i_field}) = nan(num.levels, num.levels, num.conds, 4);
            elseif contains(fieldname_list{i_field}, {'start_nll','nll','exit_flag','r2'})
                sd.ind(subj).(fieldname_list{i_field}) = nan(num.levels, num.levels, num.conds);
            end

        end

    end

    %% Super subject

    sd.all.start_params = nan(num.levels, num.levels, num.conds, 4);
    sd.all.start_nll = nan(num.levels, num.levels, num.conds);
    sd.all.params_est = nan(num.levels, num.levels, num.conds, 4);
    sd.all.nll = nan(num.levels, num.levels, num.conds);
    sd.all.exit_flag = nan(num.levels, num.levels, num.conds);
    sd.all.r2 = nan(num.levels, num.levels, num.conds);

end