function corr = computeWeightedAutocorrelation(img, Na, mask, normalize)
% assert(mod(Na,2) == 1);
% half = floor(Na/2);
% corr = zeros(Na,Na);
% 
% img = img - mean(img(:));
% 
% [ht wt] = size(img);
% 
% [ngrid mgrid] = meshgrid(-half:half, -half:half);
% 
% 
% numPixels = ((ht-half) - (half+1) + 1) * ((wt-half) - (half+1) + 1);
% K = ceil(Na*Na/2);
% for direction = 1 : K
%     n = ngrid(direction);
%     m = mgrid(direction);
%     for ii = (half+1) : (ht-half)
%     for jj = (half+1) : (wt-half)
%         d = img(ii,jj)*img(ii+n,jj+m);
%         corr(n+half+1,m+half+1) = corr(n+half+1,m+half+1) + d;
%     end 
%     end
%     % be symmetric
%     corr(-n+half+1,-m+half+1) = corr(n+half+1,m+half+1);
% end
% corr = corr/numPixels;


% assert(mod(Na,2) == 1);
half = floor(Na/2);
corr = zeros(Na,Na);

%img = (img - mean_forTTM(img, mask)).*mask;
numPixels = sum(mask(:));
wImg = img.*mask;
wMu  = sum(wImg(:))/numPixels;
%img  = wImg-wMu;

[ht wt] = size(img);

[ngrid mgrid] = meshgrid(-half:half, -half:half);

K = ceil(Na*Na/2);
for direction = 1 : K
    n = ngrid(direction);
    m = mgrid(direction);
    corr(n+half+1, m+half+1) = sum(vector(wImg.*circshift(wImg, [n m])));
    
    if (normalize)
        nfact = sum(vector(mask.*circshift(mask, [n m])));
        corr(n+half+1, m+half+1) = corr(n+half+1, m+half+1) / numPixels - wMu^2;
    else
        corr(n+half+1, m+half+1) = corr(n+half+1, m+half+1) - wMu^2*numPixels;
    end
    
    % be symmetric
    corr(-n+half+1,-m+half+1) = corr(n+half+1,m+half+1);
end