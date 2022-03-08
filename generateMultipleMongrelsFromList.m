    function generateMultipleMongrelsFromList(list_file,job_file)
%  generateMultipleMongrelsFromList(list_file,job_file)
%
% Generates one mongrel per line in a text file (columns separated by tabs
% or spaces)
% 
% list_file: filename of file with one line per image to be synthesized.
%       Each line specifies the original image file, the x and y 
%       coordinates of the fixation, the fovea radius in pixels, and an 
%       index number (used to number the mongrels when one is 
%       making multiple from the same fixation point and input image. When 
%       specifying the original imagefile, if the original image is not in 
%       the same folder with this script, please combine image directory 
%       with image name, eg. "imageMaterial/image1.png".
% job_file: filename of file that specifies parameters for synthesis. These
%       are parameters that are not image-specific, like pooling rate
%       (pooling region width as a function of eccentricity), etc. If this
%       field is not provided, then the progaram automatically runs the
%       default jobfile.
% im_dir and job_file are simply passed on to synthesize_mongrel

% License stuff:
% generateMultipleMongrelsFromList: Run Texture Tiling Model to predict
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

if nargin == 1
    job_file = 'default.job';
end

%% 1. START PARALLEL PROCESSING
% p = gcp('nocreate');
% if isempty(p)
%     c = parcluster('local'); % build the 'local' cluster object
%     nw = c.NumWorkers; % get the number of workers
%     parpool(nw-1); % open parallel pool, leaving one worker for CPU managment
% end

%% 2. LOAD THE LIST OF IMAGES AND PARAMETERS
fid = fopen(list_file);
% mongrel list is a cell structure with one element per column
mongrel_list = textscan(fid, '%s%d%d%d%d', 'CommentStyle', '%', 'MultipleDelimsAsOne', 1);
fclose(fid);
num_mongrels = size(mongrel_list{1}, 1);
im_name = mongrel_list{1};
im_fixation_x = double(mongrel_list{2});
im_fixation_y = double(mongrel_list{3});
fovea_size = double(mongrel_list{4});
mongrel_index = mongrel_list{5};
% mongrel_index is used to number output filenames if generating a bunch of
% mongrels from the same input image + fixation point

%% 3. GENERATE THE MONGRELS
for ii = 1:num_mongrels 
    synthesizeMongrel(im_name{ii}, ...
        im_fixation_x(ii), im_fixation_y(ii), fovea_size(ii), ...
        mongrel_index(ii),0,job_file);
end

