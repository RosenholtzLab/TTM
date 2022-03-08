function [paddedimg, fx, fy, poolingRegions, padding] = ...
    padToFixation(img, fixPt, foveaSize, poolingRate, radialOverlap, numAngular, padval, out_dir)
% [paddedimg, fx, fy, poolingRegions] = padToFixation(img, fixPt, foveaSize, 
%       poolingRate, radialOverlap, numAngular, padval, out_dir)
%
% Computes pooling regions and pads the input image such that the outermost 
% regions in each radial direction from fixation lie completely outside the 
% original image.
% (In early experiments, this seemed to give the best mongrel results.)
% 
% INPUTS
% img: image to be mongrelized
% fixPt: [fx fy] fixation in original image
% foveaSize: radius, in pixels, of the "fovea"
% poolingRate: poolingRate * eccentricity gives poolingRegion size
% radialOverlap: amount that regions overlap in a radial direction
% numAngular: number of pooling regions at a given eccentricity
% padval: what color to pad with
% out_dir: the directory where the figure showing pooling regions is saved
%
% OUTPUTS
% paddedimg: image after padding
% fx, fy: the same fixation, in padded image coordinates
% poolingRegions: a list of pooling regions, each line of the form:
%       [p_x, p_y, a, b, buffer]
%       where (p_x, p_y) is the center, a the lower bound, b the upper bound, 
%       and buffer a small region used for blending neighboring regions
% padding: [padLeft, padRight, padTop, padBot], amount of padding on each
%       side
%
% Copyright (C) 2019  Ruth Rosenholtz, Alvin Raj, & Lavanya Sharan
% See COPYING.txt for license details

if nargin == 7
    out_dir = [];
end
fx = fixPt(1); fy = fixPt(2);
[h,w,d] = size(img);

%% 1. DEFINE POOLING REGION CENTERS AND SIZES

% geometric ratio is a different way of thinking about the increase in size
% of the pooling regions
% note that poolingRate specifies the rate at which the *width* of the 
% pooling regions grow
geoRatio = (2+poolingRate-radialOverlap*poolingRate*2)/(2-poolingRate);

% pooling regions radiate outwards from the "fovea" at a fixed number of angles
% here we create a list of those angles
aList = [0:2*pi/numAngular:2*pi]; aList = aList(1:end-1);
% notdone is a flag for each angle
% starting value is 1, meaning not finished with this radial, yet
% set to 0 (done) when reach the first pooling region at that
% angle that falls completely outside the original image
notdone = ones(length(aList)); 
poolingRegions = [];
itr = 0;    % index of radial rings, starting with 0
ang_idx = []; % will hold angle indices
while sum(notdone)>0,   % if some angles still not done
    % log size and placement of regions
    r = foveaSize*geoRatio^itr;
    for ita=1:length(aList)
        if notdone(ita) % still working on this angle
            % offset of poolingRegion center from fixation point
            p_x = round(r*cos(aList(ita)));
            p_y = round(r*sin(aList(ita)));
            s = round(r*poolingRate); 
            a = round(s/2); b = s-a; 
            buff = max(3,ceil(a/4)); 
            
            % save this pooling region in our list
            % note these are pooling regions relative to fixation
            % we'll convert them to pixels once we figure out the padding, 
            % and the coordinates of the fixation after padding
            poolingRegions = [poolingRegions;p_x,p_y,a,b,buff];
            
            % check if this pooling region is the first at this angle to
            % lie entirely outside the original image, i.e. if 
            % (min x > w OR max x < 1 OR min y > h OR max y < 1)
            if (fx+p_x-a-buff+1 > w || fx+p_x+b+buff < 1 || ...
                    fy+p_y-a-buff+1 > h || fy+p_y+b+buff < 1)
                notdone(ita) = 0; % notdone is false, we're done with this angle
            end
            ang_idx = [ang_idx;ita];
        end
    end % do another angle
    itr = itr + 1; % increment the radius to the next ring
end

%% 2. CHOOSE PADDING TO ENSURE THAT THE PADDED IMAGE JUST BARELY INCLUDES ALL POOLING REGIONS

% to find the right padding, need to find the minimum and maximum values of
% all pooling regions, in both x and y
% convert this to a pad amount, and pad

% unpacking the poolingRegions array just so the code is clearer
clear p_x p_y a b buff
p_x = poolingRegions(:,1);
p_y = poolingRegions(:,2);
a = poolingRegions(:,3);
b = poolingRegions(:,4);
buff = poolingRegions(:,5);
start_x = fx + p_x - a - buff + 1;
end_x = fx + p_x + b + buff;
start_y = fy + p_y - a - buff + 1;
end_y = fy + p_y + b + buff;

min_x = min(start_x);
max_x = max(end_x);
min_y = min(start_y);
max_y = max(end_y);

% pad[Left/Right/Bot/Top] is the amount to add in a particular direction
% for min_x and min_y, add a pixel; the minimum will be 
% something like -2, which means we need to account for -2, -1, and 0 
% before we get to the image at 1, ...
padLeft = abs(min_x)+1; 
padTop = abs(min_y)+1;
% for max_x and max_y, just pad by the amount that the max is greater than
% the current dimension. E.G., if max is 5 and w is 3, then pad by 2
padRight = max_x-w;
padBot = max_y-h;

% actually pad here:
for i=1:d,
    tmp{i} = padarray(img(:,:,i),[0 padLeft],padval(i),'pre') ;
    tmp{i} = padarray(tmp{i},[0 padRight],padval(i),'post') ;
    tmp{i} = padarray(tmp{i},[padTop 0],padval(i),'pre') ;
    tmp{i} = padarray(tmp{i},[padBot 0],padval(i),'post') ;
    paddedimg(:,:,i) = tmp{i};
end
padding = [padLeft padRight padTop padBot];

% now that the image is padded, need to add the padding to fx, fy, and 
% the pooling regions, to convert them to their new pixel values
fx = fx + padLeft;
fy = fy + padTop;
poolingRegions(:,1) = poolingRegions(:,1) + fx;
poolingRegions(:,2) = poolingRegions(:,2) + fy;

% mark pooling region boundaries with buffer zone added and plot on screen
[H, W, ~] = size(paddedimg);
figure; image(paddedimg); axis('image');
% length(poolingRegions);
for itp1=1:length(poolingRegions)
    buffPixels = poolingRegions(itp1,5);            
        p_x = poolingRegions(itp1, 1);
        p_y = poolingRegions(itp1, 2);
        a = poolingRegions(itp1, 3);
        b = poolingRegions(itp1, 4);
        if (p_x-a-buffPixels+1 > 0 & p_x+b+buffPixels < W & p_y-a-buffPixels+1 > 0 & p_y+b+buffPixels < H & ang_idx(itp1) == 1)
            buffPixels = 0;
            X = [p_x-a-buffPixels+1 p_x-a-buffPixels+1;
                p_x+b+buffPixels, p_x+b+buffPixels;
                p_x-a-buffPixels+1 p_x+b+buffPixels;
                p_x-a-buffPixels+1 p_x+b+buffPixels]';
            Y = [p_y-a-buffPixels+1 p_y+b+buffPixels;
                p_y-a-buffPixels+1 p_y+b+buffPixels;
                p_y-a-buffPixels+1 p_y-a-buffPixels+1;
                p_y+b+buffPixels p_y+b+buffPixels]';
            line(X, Y,'LineWidth',1)
        end      
end
% save a picture of the pooling regions 
save_name = strcat(out_dir,'/poolingRegion.png');
saveas(gcf,save_name);

return;