function img_contrast = calc_contrast(img,contrast_type)

switch contrast_type

    case 'Michelson'

        img_contrast = ( max(img(:)) - min(img(:)) ) / ( max(img(:)) + min(img(:)) );

    case 'RMS'

        % Image is assumed to have pixel values between 0 and 1
        img_contrast = std2(img);

end

end