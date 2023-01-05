function synthesizeMongrel(img_file,fx,fy,foveaSize,mongrel_idx, verbose,job_file, output_folder)   
% synthesizeMongrel(img_file,fx,fy,foveaSize,mongrel_idx, verbose,job_file,output_folder)   
%
% Synthesizes a single mongrel
% 
% img_file: input image
% fx, fy: fixation coordinates 
% foveaSize: radius of the fovea in pixels 
% mongrel_idx: when generating multiple mongrels with the same parameters, 
%       including fixation, this is used to number the output files
% verbose: specifies how verbose to be in messages to the terminal
% job_file: filename of file that specifies parameters for synthesis
% 
% Presumes that the code needed sits within the current directory and
% sub-directories

% License stuff:
% synthesizeMongrel: Creates a single 'mongrel', a visualization of the
%      information available in peripheral vision, given a fixation
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

no_save = 1;

addpath(genpath('./'));

if nargin < 5
    mongrel_idx = 1;
end
if nargin < 6
    verbose = 0;
end
if nargin < 7
    job_file = 'default.job';
end 
prm = loadStructureFile(job_file);


%% 1. LOAD IMAGE, SPECIFY FIXATION
% load original image
cmap = [];
im = im2double(imread(img_file));
%check if this is a greyscale image with three color channels repeated.
disp('loaded image equal mannual tolerance large tol')
if size(im,3)==3
   channel_12 = abs(im(:,:,1)-im(:,:,2)) < 1e-1;
   channel_13 = abs(im(:,:,1)-im(:,:,3)) < 1e-1;
   sum_12 = sum(channel_12,'all')
   sum_13 = sum(channel_13,'all')
   % disp(norm(im(:,:,1)-im(:,:,2),2))
   % disp(norm(im(:,:,1)-im(:,:,3),2))
   % disp(norm(im(:,:,3)-im(:,:,2),2))
   % if norm(im(:,:,1)-im(:,:,2),2)<=25 & norm(im(:,:,1)-im(:,:,3),2)<=25 & norm(im(:,:,3)-im(:,:,2),2) <=25
   if sum(channel_12,'all') == numel(channel_12) | sum(channel_13,'all') == numel(channel_13)
   % if abs(sum(channel_12,'all') - numel(channel_12)) < 100 & abs(sum(channel_13,'all'), numel(channel_13))<100
      im = mean(im,3);
      prm.colorSynth = 0;
      disp('gray')
   end
end

x = round(fx);
y = round( fy);

image_size = size(im)
image_name = img_file


%% 2. CHECK IF IMAGE IS COLOR/GRAYSCALE
% convert indexed-color images to rgb
if ~isempty(cmap)
    im = (ind2rgb(im,cmap));
end
% perform colorSynth based on whether the image is actually color/grayscale
if size(im,3) == 1
    prm.colorSynth = 0;
    componentInfo = [];
    msgbox('Performing grayscale synthesis.');
    fprintf('Performing grayscale synthesis');
% elseif (sum(sum(abs(im(:,:,1)-im(:,:,2))))<=5 && sum(sum(abs(im(:,:,2)-im(:,:,3))))<=5)
elseif (sum(sum(abs(im(:,:,1)-im(:,:,2))))==0 && sum(sum(abs(im(:,:,2)-im(:,:,3))))==0)
%elseif ((sum(sum(ismembertol(im(:,:,1), im(:,:,2))))==0) && ((sum(sum(ismembertol(im(:,:,2),im(:,:,3)))))==0))
    im = mean(im,3);
    prm.colorSynth = 0;
    componentInfo = [];
    msgbox('Performing grayscale synthesis.');
    fprintf('Performing grayscale sysnthesis.');
else
    prm.colorSynth = 1;
    fprintf('Perforring color synthesis');
    msgbox('Performing color synthesis.');
    % get axes of color space and convert to that space
    if prm.colorSynth
        [im_ica componentInfo] = separateImage(im, 'ica');
    end
end


%% 3. DIRECTORY HANDLING AND OTHER ADMINISTRATION
% output to a folder (within the same folder of the code) named
% image_fixation_date_mongrelIndex
[~, imname] = fileparts(img_file);
%out_dir = strcat(imname,'_X',num2str(round(fx)),'_Y',num2str(round(fy)),'_',datestr(now,29),'_',char(mongrel_idx));
out_dir = strcat(output_folder,'/ecc_',num2str(round(fx)));
% create folder as needed

curr_dir = dir;
tmp_i = 1;

while sum(strcmp({curr_dir.name},out_dir))
    namesplit = strsplit(out_dir,'_');
    if ~strcmp(namesplit{end},'1') && tmp_i == 1
    out_dir = strcat(out_dir,'_',num2str(tmp_i));
    else
        tmp_i = tmp_i+1;
    out_dir = strcat(out_dir(1:end-1),num2str(tmp_i));
    end
end
mkdir(out_dir);
im_pad_path = strcat(out_dir, '/pad_',imname,'_',string(mongrel_idx),'.png');
img_mask_path = strcat(out_dir, '/imgmask',imname,'_',string(mongrel_idx),'.png');
output_parameter_path = strcat(out_dir, '/synthesis_parameter_',imname,'_',string(mongrel_idx),'.txt'); % for storing parameters.
% out_path =  strcat(out_dir, '/mongrelized_',imname,'_',string(mongrel_idx),'.jpg'); % for storing synthesized mongrels.
out_path =  strcat(out_dir, '/mongrel_',imname,'_ecc_',num2str(round(fx)),'_',string(mongrel_idx),'.jpg');
parameter_output_path =  strcat(out_dir, '/parameters_',imname,'_',string(mongrel_idx),'.txt'); % for storing parameters used in the synthesis.
parameter_output_path =  strcat(out_dir, '/parameters_',imname,'_',string(mongrel_idx),'.txt'); % for storing 
debug_path = strcat(out_dir, '/debug/');
if prm.colorSynth
    addpath(genpath('FastICA/'));
end

% set seed, parameters necessary for generating a random mongrel
% 'shuffle' mode sets the seed based on the current time
% If you don't do something like this then silly MATLAB will generate the
% same mongrel every time it starts up, e.g. you will start 10 processes to
% generate 10 mongrels and get the same mongrel 10 times.
if strcmp(prm.random_mode,'shuffle')
    rng('shuffle');
else
    % sometimes you don't want to set the seed randomly, because you are
    % trying to regenerate a mongrel starting with different parameters, or
    % you are generating a mongrel for frame n+1 of a video and want to use
    % the same seed as frame n, etc. 
    % in these cases, read in a pre-defined seed
    try
        rng(prm.random_seed,prm.random_type);
    catch
        sprintf('Random seed set incorrectly. Use random_mode "shuffle" or define random_type and random_seed in jobfile.');
    end
end
curr_random_seed = rng; 
curr_random_seed_Type = curr_random_seed.Type;
curr_random_seed_Seed = curr_random_seed.Seed;
prm.savename = sprintf('intermediateResults_%s',imname); %folder name for saving intermediate results

%% 4. PAD ORIGINAL IMAGE AND BLUR TO MIMIC ACUITY LOSS
% pad the image, and save it (used as input to TTM code)
% somewhat arbitrarily make the padding the same color as the top left pixel
%   (this is a good idea for typical psychophysical displays with a blank 
%   background, but it is ok for the user to change this)
padding_color = im(1,1,:);

if isnan(y)
    latticeType = 'rhombic'; %grid rhombic hexagonal
    disp(sprintf('Second Fixation coordinate is nan, running uniform pooling at %d pixels eccentricity.',x));
    disp(sprintf('Using %s Lattice',latticeType));
    [im_pad,fx,fy,poolingRegions,padding] = padToFixationUniform(im,[x,y], foveaSize,prm.poolingRate, prm.radialOverlap, prm.numAngular, padding_color,out_dir,latticeType);
    
    fixPt = [fx,fy]; % after padding, the fixation point has new coordinates

    % add peripheral blur to original image
    im_pad = simulatePeripheralBlurUniform(im_pad, fixPt, foveaSize);
    
else
    disp(sprintf('Fixation coordinate is %d,%d, running foveated pooling.\n',x,y));
    [im_pad,fx,fy,poolingRegions,padding] = padToFixation(im,[x,y], foveaSize,prm.poolingRate, prm.radialOverlap, prm.numAngular, padding_color,out_dir);
    
    fixPt = [fx,fy]; % after padding, the fixation point has new coordinates

    % add peripheral blur to original image
    im_pad = simulatePeripheralBlur(im_pad, fixPt, foveaSize);
    
end


% since some old routines grab the padding color from the corner, reinstate
% that here, after the blurring operation which might have changed it
im_pad(1,1,:) = padding_color; 
% save the padded, acuity-adjusted input image
if no_save == 0
    imwrite(im_pad, im_pad_path);
end
% *** remove saving ***
% imwrite(im_pad, im_pad_path);
% *** remove saving ***

%% 5. CREATE A MASK SHOWING ORIGINAL IMAGE VS. PADDING
% if useImgMask is set (typical use), then on some synthesis iterations we
% reimpose the constraint that we know where the original (unpadded) image
% was located. In informal experiments this speeded convergence without
% greatly affecting the results
useImgMask = 1;
if ~useImgMask
    maskimg = [];
else    
    % create maskimg, which indicates location of the original image
    % padded region indicated in black (0), image region is white (1)
    maskimg = padarray(ones(size(im,1),size(im,2)), [0 padding(1)], 0, 'pre') ;
    maskimg = padarray(maskimg,[0 padding(2)], 0, 'post') ;
	maskimg = padarray(maskimg,[padding(3) 0], 0, 'pre') ;
	maskimg = padarray(maskimg,[padding(4) 0], 0, 'post') ;
    if no_save == 0
        imwrite(maskimg, img_mask_path);
    end
    % *** remove saving ***
    % imwrite(maskimg, img_mask_path);
    % *** remove saving ***
end

%% 6. CREATE OUTPUT PARAMETERS FILE
prm.foveaSize = foveaSize; 
prm.image = im; % input image
prm.image_name = imname;
prm.paddingColor = padding_color;
prm.input_path = img_file;
prm.out_dir = out_dir;
prm.im_pad_path = im_pad_path;
prm.image_mask_path = img_mask_path;
prm.parameter_output_path = output_parameter_path;
prm.out_path = out_path;
prm.debugPath = debug_path;
prm.intermediate_result_folder = prm.savename;
prm.source_img_fixation_x = x;
prm.source_img_fixation_y = y;
prm.curr_random_seed = curr_random_seed;
prm.componentInfo = componentInfo;
prm.poolingRegions = poolingRegions;
prm.useImgMask = useImgMask;
try
    if no_save == 0
        save(strcat(out_dir,'/parameter.mat'),'prm');
    end
    % *** remove saving *** 
    %save(strcat(out_dir,'/parameter.mat'),'prm');
    % *** remove saving ***
catch
    keyboard
end

%check if file is open and wait until not open to use
filename = strcat(output_folder,"/imgname_ecc_randomseeds_",num2str(prm.source_img_fixation_x),".txt");
fid = fopen(filename,'a');
fwrite(fid,sprintf('%s\t%s\t%d\n', prm.image_name,num2str(round(fx)),prm.curr_random_seed.Seed));
%fwrite(fid,sprintf('%s\t',num2str(round(fx))));
%fwrite(fid,sprintf('%d\n', prm.curr_random_seed.Seed));
fclose(fid);

% print parameters to output file
% fid = fopen(parameter_output_path,'wt');
% fwrite(fid,sprintf('useImgMask ''%d''\n', prm.useImgMask)); 
% fwrite(fid,sprintf('random_mode ''%s''\n', prm.random_mode)); 
% fwrite(fid,sprintf('random_type ''%s''\n', prm.curr_random_seed.Type)); 
% fwrite(fid,sprintf('random_seed ''%d''\n', prm.curr_random_seed.Seed)); 
% fwrite(fid,sprintf('saveIntermediate ''%d''\n', prm.saveIntermediate)); 
% fwrite(fid,sprintf('poolingRate ''%d''\n', prm.poolingRate)); 
% fprintf(fid, sprintf('scalesToRun [%d, %d, %d, %d]\n', ...
%     prm.scalesToRun(1), prm.scalesToRun(2), prm.scalesToRun(3), prm.scalesToRun(4)));
% fwrite(fid,sprintf('Nor ''%d''\n', prm.Nor)); 
% fwrite(fid,sprintf('Na ''%d''\n', prm.Na)); 
% fprintf(fid, sprintf('nIters [%d, %d, %d, %d]\n', ...
%     prm.nIters(1), prm.nIters(2), prm.nIters(3), prm.nIters(4)));
% fwrite(fid,sprintf('prIters ''%d''\n', prm.prIters)); 
% fwrite(fid,sprintf('radialOverlap ''%d''\n', prm.radialOverlap)); 
% fwrite(fid,sprintf('numAngular ''%d''\n', prm.numAngular)); 
% fwrite(fid,sprintf('colorSynth ''%d''\n', prm.colorSynth)); 
% fwrite(fid,sprintf('foveaSize ''%d''\n', prm.foveaSize)); 
% fwrite(fid,sprintf('fix_x ''%d''\n', prm.source_img_fixation_x)); 
% fwrite(fid,sprintf('fix_y ''%d''\n', prm.source_img_fixation_y)); 
% fwrite(fid,sprintf('image_name ''%s''\n', prm.image_name)); 
% fwrite(fid,sprintf('out_directory ''%s''\n', prm.out_dir)); 
% fclose(fid);
'sizes'
size(im_pad,3)
size(im,3)
%% 7. RUN TTM
% run the job
% "stitch" together multiple pooling regions to create a continuous mongrel 
% consistent with all measured statistics (approximately, of course) 
reconsb = [];
for ii = 1 : length(prm.scalesToRun)
    maxscale = prm.scalesToRun(ii); % lower numbers are coarser scales
    nIterStart = 1;
    [reconsb] = stitchBareBones(im_pad, componentInfo, [maxscale, max(prm.scalesToRun)], ...
        fixPt, reconsb, maskimg, prm.Nor, prm.Na, [prm.nIters(ii),nIterStart], ...
        prm.foveaSize, prm.prIters, ...
        prm.poolingRegions, prm.savename, ...
        prm.saveIntermediate, prm.debugPath, verbose,prm.colorSynth);
end

%% 8. SAVE RESULT
outname = prm.out_path;
if ~isempty(reconsb)
    % remove padding
    res = reconsb(padding(3)+1:padding(3)+size(im,1), padding(1)+1:padding(1)+size(im,2), :);
    imwrite(res, outname, 'jpg', 'Quality', 100);

end
end
