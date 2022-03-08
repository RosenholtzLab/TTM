function [img componentInfo] = separateImage(img, colorMethod)
% [out componentInfo] = separateImage(img, colorMethod)
% 
% Convert input image to another colorspace, e.g. ica space
% img: input image
% colorMethod: 'ica' or 'rgb'; or previously computed independent
%       components
% out: converted image in the new colorspace
% componentInfo: indicates independent components used for conversion

    [H W D] = size(img);
    vectorImg = reshape(permute(img,[3 1 2]), [D, H*W]);    
    [c1 c2 c3 componentInfo] = separateComponents(vectorImg, colorMethod);
    img(:,:,1) = reshape(c1, [H W]);
    img(:,:,2) = reshape(c2, [H W]);
    img(:,:,3) = reshape(c3, [H W]);
return

function [c1,c2,c3,componentInfo] = separateComponents(data, method)

    if (ischar(method))
        switch(lower(method))
            case 'ica'
                [c123, mix, sep] = fastica(data,'numOfIC',3,'approach','symm','g','tanh');
                componentInfo.method = method;
                componentInfo.mix = mix;
                componentInfo.sep = sep;
            case 'rgb'
                c123 = data;
                componentInfo.method = method;
                componentInfo.mix = eye(3);
                componentInfo.sep = eye(3);
            otherwise
                error('unknown method');
        end
    elseif (isstruct(method))
        c123 = method.sep * data;
        componentInfo = method;
    else
        error('Unknown method. Please give a string, or a structure');
    end
       
    c1 = c123(1,:);    
    if size(c123,1)>1
        c2 = c123(2,:);
        if size(c123,1)>2
            c3 = c123(3,:);
        else
            c3 = zeros(1,size(c123,2));
        end
    else
        c2 = zeros(1,size(c123,2));
        c3 = zeros(1,size(c123,2));
    end
return
