function img = centerTextureContrast(img, contrast, gray)

    img = img - mean(img(:));
    max_abs_val = max(abs(img(:)));
    
    if max_abs_val > 0
        img = img / (2 * max_abs_val);
    end 

    img = gray * (1 + img * 2 * contrast);

end