
function normalized_array = normalize_array(array, method)
    switch method
    
        case 'z-score'
            
            normalized_array = (array - mean2(array)) ./ std2(array);

        case 'min-max'

            normalized_array = (array - min(array(:))) ./ (max(array(:)) - min(array(:)));
    
    end
end