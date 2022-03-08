function [fix_center,fovea_radius] = convertImageForMongrels(im_name,varargin)
% [fix_center,fovea_radius] = convertImageForMongrels(im_name,varargin)
%
% This function takes in an image and allows you to select the fixation and
% fovea size by clicking on the image. It generates a text file for which
% each row specifies the parameters to run TTM in the following format:
% directory(if not in current folder)/filename [tab] fix_x [tab] fix_y [tab] radius-of-fovea [tab] mongrel_index(from 1 to num_mongrels)
    
% Required parameter:
% im_name: name of the image file, including the path and extension

% Optional parameters are below. Set as you would plotting properties, e.g. 
% convert_image_for_mongrels('foo.img','num_mongrels',10,'do_zoom',1)

% num_mongrels: pass the number of mongrels to generate from the image. If
% left unset, defaults to 1.

% manual_fix: pass -1 for center fixation, or a 1 by 2 vector [X Y].
% Do not set if you want to select the center by hand.

% fovea_radius: pass the width of the fovea in pixels, used by TTM as
% pixels-per-degree. Do not set if you want to select the radius by hand.
% This is the RADIUS, not the radius, of the fovea.

% do_zoom: pass 1 if you want to first click near the fixation and
% zoom in before clicking the actual fixation (useful for big images).

% destination_folder: pass the directory to write the text file to. If 
% left unset it defaults to the input image's directory. 

% Note that this will OVERWRITE an existing text file WITHOUT prompting

% License stuff:
% convertImageForMongrels: Select fixation and fovea size for modeling 
%      human peripheral vision
% Copyright (C) 2019  Ruth Rosenholtz, Alvin Raj, & Lavanya Sharan.
%
% This program is free software; you can redistribute it and/or
% modify it under the terms of the GNU General Public License
% as published by the Free Software Foundation; either version 2
% of the License, or (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.

fprintf('Texture Tiling Model release 1.0, Copyright (C) 2019\n');
fprintf('Ruth Rosenholtz, Alvin Raj, and Lavanya Sharan\n');
fprintf('The Texture Tiling Model comes with ABSOLUTELY NO WARRANTY.\n');
fprintf('This is free software, and you are welcome to redistribute it\n')
fprintf('under certain conditions. See COPYING.txt for details.\n');

assert(mod(length(varargin),2)==0,'Optional arguments must have names and values')

possible_inputs = {'num_mongrels','do_zoom','manual_fix',...
    'fovea_radius','destination_folder'};
    
[filepath,name,ext] = fileparts(im_name);

%convert_to_grayscale = 0;
num_mongrels = 1;
do_zoom = 0;
manual_fix = [];
fovea_radius = [];
destination_folder = filepath;

% read in optional parameters
for ii = 1:2:length(varargin)

    index = find(strcmp(varargin{ii},possible_inputs), 1);
    assert(~isempty(index),['Could not parse parameter: ' varargin{ii}])
    if ischar(varargin{ii+1})
        destination_folder = varargin{ii+1};
    else
        eval([varargin{ii} '=[' num2str(varargin{ii+1}) '];']);
    end
end

% read in image
[im,~] = imread(im_name);

% set manual fixation parameters if given
if ~isempty(manual_fix)
    if manual_fix == -1
        fix_center = round([size(im,2) size(im,1)]/2);
    else
        fix_center = manual_fix;
    end
    
    % set the zoom to center on the fixation center
    zoom_center = fix_center;
    
else % input fixation parameters using ginput
    
    % if requested zooming for ease of use
    if do_zoom
        imshow(im); title('CLICK TO ZOOM');
        set(gcf,'units','normalized','outerposition',[0 0 1 1])
        zoom_center = round(ginput(1));
    end
    
    % show the image, allow user to click fixation center
    imshow(im); title('CLICK FIXATION CENTER');
    set(gcf,'units','normalized','outerposition',[0 0 1 1])
    if do_zoom % if zooming
        xlim([zoom_center(1)-70 zoom_center(1)+70])
        ylim([zoom_center(2)-70 zoom_center(2)+70])
    end
    
    fix_center = round(ginput(1));
    
end

if isempty(fovea_radius)
    
    if do_zoom % if zooming
        xlim([zoom_center(1)-70 zoom_center(1)+70])
        ylim([zoom_center(2)-70 zoom_center(2)+70])
    end

    % show the image, allow user to click fixation boundaries
    imshow(im); title('CLICK WHEN RADIUS IS CORRECT');
    hold on
    plot(fix_center(1),fix_center(2),'r+') % show fixation center
    
    set(gcf,'units','normalized','outerposition',[0 0 1 1])
    if do_zoom
        xlim([zoom_center(1)-70 zoom_center(1)+70])
        ylim([zoom_center(2)-70 zoom_center(2)+70])
    end
    
    % get the fovea
    set (gcf, 'WindowButtonMotionFcn', @mouseMove);
    set (gcf, 'name', num2str(fix_center));
    k = waitforbuttonpress;
    C = get (gca, 'CurrentPoint');
    
    fovea_radius = round(sqrt((fix_center(1)-C(1,1))^2 + (fix_center(2)-C(1,2))^2));
end

close

% create data file to pass to generate_multiple_mongrels
fprintf('Saving file: %s\n',fullfile(destination_folder,[name ext '_mongrel_list.txt']));
fid = fopen(fullfile(destination_folder,[name ext '_mongrel_list.txt']),'w+');

% add rows
for ii = 1:num_mongrels
    fprintf(fid, '%s\t%d\t%d\t%d\t%d',im_name,fix_center(1),fix_center(2),fovea_radius,ii);
    fprintf(fid,'\n');
end
fclose(fid);

end

% for plotting fovea
function mouseMove (object, eventdata)
C = get (gca, 'CurrentPoint'); % get the mouse position
fix_center = str2num(get(gcf, 'name')); % get the fixation center
r = sqrt((fix_center(1)-C(1,1))^2 + (fix_center(2)-C(1,2))^2); % compute the radius

h = findobj('type','line');
delete(h(1)) % find and delete the old line

t = linspace(0,2*pi);
plot(r*cos(t)+fix_center(1),r*sin(t)+fix_center(2),'r'); % plot the circle in red

drawnow
end