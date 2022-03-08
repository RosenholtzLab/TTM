function maskedskew = skew_forTTM(img, mask, mu, sigSq)
    
    
    if (~exist('mu','var') || isempty(mu))
        mu = mean_forTTM(img, mask);
    end
    if (~exist('sigSq','var') || isempty(sigSq))
        sigSq = var_forTTM(img, mask, mu);
    end
    mask = mask(:) > .5;
    img = img(mask);
    if (isreal(img))
        img = (img(:)-mu).^3;        
        maskedskew = sum((img)/sigSq^(3/2)) / sum(mask);
    else
        imgr = real(img(:)-mu).^3;
        imgi = imag(img(:)-mu).^3;
        maskedskewr = sum((imgr)/real(sigSq)^(3/2)) / sum(mask);
        maskedskewi = sum((imgi)/imag(sigSq)^(3/2)) / sum(mask);
        maskedskew  = maskedskewr + 1i*maskedskewi;
    end
return