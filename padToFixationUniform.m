function [paddedimg, fx, fy, poolingRegions, padding] = ...
    padToFixationUniform(img, fixPt, foveaSize, poolingRate, radialOverlap, numAngular, padval, out_dir, latticeType)
% [paddedimg, fx, fy, poolingRegions] = padToFixation(img, fixPt, foveaSize, 
%       poolingRate, radialOverlap, numAngular, padval, out_dir, latticeType)
%
% Computes pooling regions and pads the input image such that the outermost 
% regions in each radial direction from fixation lie completely outside the 
% original image.
% (In early experiments, this seemed to give the best mongrel results.)
% 
% INPUTS
% img: image to be mongrelized
% ecc: [fx fy] fixation in original image
% foveaSize: radius, in pixels, of the "fovea"
% poolingRate: poolingRate * eccentricity gives poolingRegion size
% radialOverlap: amount that regions overlap in a radial direction
% numAngular: number of pooling regions at a given eccentricity
% padval: what color to pad with
% out_dir: the directory where the figure showing pooling regions is saved
% latticeType: how to align the pooling regions (square, rhombic, or hexagonal)
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
r = fixPt(1);
[h,w,d] = size(img);

%if not(isempty(fx) | isempty(fx))
%    error('Must Call PadToFixation not PadToFixation Uniform when fixation points fx and fy are not empty lists.');
%end

%% 1. DEFINE POOLING REGION CENTERS AND SIZES

% geometric ratio is a different way of thinking about the increase in size
% of the pooling regions
% note that poolingRate specifies the rate at which the *width* of the 
% pooling regions grow
%geoRatio = (2+poolingRate-radialOverlap*poolingRate*2)/(2-poolingRate);

% pooling regions radiate outwards from the "fovea" at a fixed number of angles
% here we create a list of those angles
%aList = [0:2*pi/numAngular:2*pi]; aList = aList(1:end-1);
% notdone is a flag for each angle
% starting value is 1, meaning not finished with this radial, yet
% set to 0 (done) when reach the first pooling region at that
% angle that falls completely outside the original image

% set pooling region size parameters
s = round(r*poolingRate); % diameter
a = round(s/2); b = s-a; % a is radius of pooling region. b is also radius plus spill-over
buff = max(3,ceil(a/4)); % buffer for edge of pooling region


%% Calulate Lattice:
latticeType
radialOverlap
dist = (1-radialOverlap)*s;
if strcmp(latticeType,'hexagonal') | strcmp(latticeType,'grid')
    % get centers of each pooling region: center regions at img center
    img_center = [round(w/2),round(h/2)];
    px_left = round(img_center(1):-(1-radialOverlap)*s:1-a);
    px_right = round(img_center(1):(1-radialOverlap)*s:w+b);
    px = cat(2,flip(px_left),px_right(2:end));
    py_left = round(img_center(2):-(1-radialOverlap)*s:1-a);
    py_right = round(img_center(2):(1-radialOverlap)*s:h+b);
    py = cat(2,flip(py_left),py_right(2:end));

    %old code to center on upper left edge
    %px = round(1:(1-radialOverlap)*s:w+b);
    %py = round(1:(1-radialOverlap)*s:h+b);
elseif strcmp(latticeType,'rhombic')
    img_center = [round(w/2),round(h/2)];
    % make the left side longer by one pooling region (2*a instead of a)
    px_left = round(img_center(1):-2*dist/sqrt(2):1-2*a);
    px_right = round(img_center(1):2*dist/sqrt(2):w+b);
    px = cat(2,flip(px_left),px_right(2:end));
    py_left = round(img_center(2):-dist/sqrt(2):1-a);
    py_right = round(img_center(2):dist/sqrt(2):h+b);
    py = cat(2,flip(py_left),py_right(2:end));
else
    error('Lattice Type must be grid, hexagonal or rhombic')
end

%create meshgrid of x and y coordinates, and calcaulate size
[p_x,p_y] = meshgrid(px,py)
num_pools_x = length(px);
num_pools_y = length(py);

px,py

%if we have a rhombic lattice, shift every other x coordiante down by half index
if strcmp(latticeType,'rhombic')
    p_x(1:2:end,1:end) = p_x(1:2:end,1:end)+round(dist/sqrt(2));
    % shifting causes an asymmetry on the left side
    % we need to remove every other pooling reigion on the left edge
    % set reigons to NaN and remove later
    p_x(2:2:end,1) = NaN;
    p_y(2:2:end,1) = NaN;

elseif strcmp(latticeType,'hexagonal')
    p_x(1:2:end,1:end) = p_x(1:2:end,1:end)+round(dist/2);
    %p_y(1:end,1:2:end) = p_y(1:end,1:2:end)+round(dist/2);
    %p_x(1:2:end,end) = p_x(1:2:end,1)-dist; %wrap last onto top
end

%% 2. CHOOSE PADDING TO ENSURE THAT THE PADDED IMAGE JUST BARELY INCLUDES ALL POOLING REGIONS

% to find the right padding, need to find the minimum and maximum values of
% all pooling regions, in both x and y
% convert this to a pad amount, and pad

% unpacking the poolingRegions array just so the code is clearer
%clear p_x p_y a b buff

p_x = p_x(:);
p_y = p_y(:);

% remove extra pooling regions
p_x = p_x(~isnan(p_x));
p_y = p_y(~isnan(p_y));


start_x = p_x - a - buff + 1;
end_x = p_x + b + buff;
start_y = p_y - a - buff + 1;
end_y = p_y + b + buff;

min_x = min(start_x);
max_x = max(end_x);
min_y = min(start_y);
max_y = max(end_y);

p_x,p_y

%plotting center pool region dots
figure(1);
imshow(img);
hold on;
plot(p_x,p_y,'.');
hold off;



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

numPoolingRegions = length(p_x);

poolingRegions = cat(2,p_x,p_y,repmat(a,numPoolingRegions,1),repmat(b,numPoolingRegions,1),repmat(buff,numPoolingRegions,1));

poolingRegions(:,1) = poolingRegions(:,1) + padLeft;
poolingRegions(:,2) = poolingRegions(:,2) + padTop;

% mark pooling region boundaries with buffer zone added and plot on screen
[H, W, ~] = size(paddedimg);
figure; image(paddedimg); 
hold on;
for itp1=1:numPoolingRegions
    buffPixels = poolingRegions(itp1,5);            
    p_x = poolingRegions(itp1, 1);
    p_y = poolingRegions(itp1, 2);
    a = poolingRegions(itp1, 3);
    b = poolingRegions(itp1, 4);
    buffPixelsViz = 0;
    X = [p_x-a-buffPixelsViz+1 p_x-a-buffPixelsViz+1;
        p_x+b+buffPixelsViz, p_x+b+buffPixelsViz;
        p_x-a-buffPixelsViz+1 p_x+b+buffPixelsViz;
        p_x-a-buffPixelsViz+1 p_x+b+buffPixelsViz]';
    Y = [p_y-a-buffPixelsViz+1 p_y+b+buffPixelsViz;
        p_y-a-buffPixelsViz+1 p_y+b+buffPixelsViz;
        p_y-a-buffPixelsViz+1 p_y-a-buffPixelsViz+1;
        p_y+b+buffPixelsViz p_y+b+buffPixelsViz]';
    line(X, Y,'LineWidth',1)  
end
axis('image');
hold off;
% save a picture of the pooling regions 
save_name = strcat(out_dir,'/poolingRegion.png');

[X, map] = frame2im(getframe(gcf));                                              

imwrite(X,save_name); 
return;