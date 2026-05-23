function sd = deriveSdStructFromSerialDependenceResult(result)
% deriveSdStructFromSerialDependenceResult  Build legacy-style sd struct for plot helpers.

    sd = struct();
    sd.all = struct();

    if isempty(result) || ~isfield(result, 'summary_table') || isempty(result.summary_table)
        sd.all.params_est = nan(3, 3, 2, 4);
        return
    end

    sd.all.params_est = packSummaryTableToSdParamsEst(result.summary_table);

end
