function maskedmean = mean_forTTM(img, mask)
    mask = mask(:) > .5;
    img  = img(mask);
    maskedmean = mean(img(:));
    %try
    %maskedmean = sum(img.*mask)/sum(mask);
    %catch
    %    keyboard
    %end
return