% what if responses is a logical that says whether the two trials are
% correct.
function contrast = func_name(curr_contrast, step, responses)

if len(responses) == 1
    if responses == 1
        contrast = curr_contrast - step;
    else 
        contrast = curr_contrast + step;
    end
else
    if responses(2) == 0
        contrast = curr_contrast + step;
    elseif responses(1) == 1 && responses(2) == 1
        contrast = curr_contrast - step;
    end
end

end