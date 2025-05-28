function [subplotX, subplotY] = get_subplot_dimensions(num_subplots)

if num_subplots > 3
    tmp = 1:num_subplots;
    subplot_dims = tmp(~(rem(numel(tmp), tmp))); %
    if mod(length(subplot_dims),2) ~= 0 % if number subplot dims. is odd, grab the middle value (m x m subplot)
        subplotX = median(subplot_dims);
        subplotY = subplotX;
    else % if even, take the middle left and right (m x n)
        subplotX = subplot_dims(length(subplot_dims)/2);
        subplotY = subplot_dims(length(subplot_dims)/2+1);
    end
else
    subplotX = 1;
    subplotY = num_subplots;
end

end