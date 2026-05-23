function c = conditionIndexSerialDependence(manipulation, prev_lvl, curr_lvl)
% Map manipulation (1=contrast, 2=precision) and 3x3 levels to 1..18

    c = (manipulation - 1) * 9 + (prev_lvl - 1) * 3 + curr_lvl;

end
