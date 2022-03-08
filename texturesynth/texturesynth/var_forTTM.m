function outvar = var_forTTM(img)
% outvar = var_forTTM(img)
% 2-D variance with added fudge factor so > 0.
    outvar = var(img(:))+1e-10;
return