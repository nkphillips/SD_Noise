function [new_probe_offsets] = update_probe_offset(correct_responses, lvl_order, current_probe_offsets, step_size, min_probe_offset, max_probe_offset)
% UPDATE_PROBE_OFFSET Adjusts probe offset based on block performance using 2-down 1-up rule
%   [new_probe_offsets] = update_probe_offset(correct_responses, lvl_order, current_probe_offsets, step_size, min_probe_offset, max_probe_offset)
%
%   Inputs:
%       correct_responses: Vector of 1s (correct) and 0s (incorrect)
%       lvl_order: Vector of level order
%       current_probe_offsets: Vector of current probe offset values for each level
%       step_size: Current step size for adjustments
%       min_probe_offset: Minimum allowed probe offset
%       max_probe_offset: Maximum allowed probe offset
%
%   Output:
%       new_probe_offsets: Updated probe offset values for each level

new_probe_offsets = current_probe_offsets;

unique_levels = unique(lvl_order);

for lvl = 1:length(unique_levels)
    
    lvl_responses = correct_responses(lvl_order == lvl);
    correct_count = sum(lvl_responses);
    performance = correct_count / length(lvl_responses);
    
    if performance >= 0.72
        new_probe_offsets(lvl) = round(max(min_probe_offset, current_probe_offsets(lvl) - step_size(lvl)));
    elseif performance < 0.68
        new_probe_offsets(lvl) = round(min(max_probe_offset, current_probe_offsets(lvl) + step_size(lvl)));
    end
end

end 