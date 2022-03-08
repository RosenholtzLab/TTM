function recons = integrateImage(recons, componentInfo, H, W)
% recons = integrateImage(input, componentInfo, H, W)
% 
% Used by color texture synthesis routines to convert ica-colorspace image back to rgb 
    c1 = vector(recons(:,:,1))';
    c2 = vector(recons(:,:,2))';
    c3 = vector(recons(:,:,3))';
    recons = integrateComponents(c1,c2,c3,componentInfo);
    try
        recons = reshape(permute(recons, [2 1]), [H, W, 3]);
    catch
        keyboard
    end
return

function [data] = integrateComponents(c1,c2,c3,componentInfo)

    switch(lower(componentInfo.method))
        case 'ica'
            data = cat(1,c1,c2,c3);
            data = componentInfo.mix * data;
        case 'rgb'
            data = cat(1,c1,c2,c3);
        otherwise
            error('unknown method');
    end
    
return