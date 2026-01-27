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

% staircase rule (2-up-1-down)
if t == 1
    if resp(t) == 1
        contrast(t+1) = contrast(t) - step;  % harder
    else
        contrast(t+1) = contrast(t) + step;  % easier
    end
elseif resp(t-1) == 1 && resp(t) == 1
    contrast(t+1) = contrast(t) - step;  % harder
else
    contrast(t+1) = contrast(t) + step;  % easier
end