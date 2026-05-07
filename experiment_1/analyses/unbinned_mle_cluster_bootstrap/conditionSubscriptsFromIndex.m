function [manipulation, prev_lvl, curr_lvl] = conditionSubscriptsFromIndex(c)
% Linear index c in 1..18 -> manipulation 1|2, prev 1..3, curr 1..3

    manipulation = floor((c - 1) / 9) + 1;
    rem = c - (manipulation - 1) * 9;
    prev_lvl = floor((rem - 1) / 3) + 1;
    curr_lvl = rem - (prev_lvl - 1) * 3;

end
