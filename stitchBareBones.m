function [synth] = stitchBareBones(img, componentInfo, ...
    maxscale, fixPt, reconSeed, imgMask, Nor, Na, ...
    nItersVec, foveaSize, prIters, poolingRegions, savename, ...
    saveIntermediateFlag, debugPath, verbose, colorSynth)
% [synth] = stitchBareBones(img, componentInfo, ...
%     maxscale, fixPt, reconSeed, seedMask, Nor, Na, ...
%     nItersVec, foveaSize, prIters, poolingRegions, savename, ...
%     saveIntermediateFlag, debugPath, verbose)
%
% Synthesize ("stitch") a mongrel with a bare bones model that uses square 
% pooling regions (not biologically plausible nor in line with psychophysics, 
% but fast-ish and low on artifacts. It treats Portilla & Simoncelli's 
% texture analysis/synthesis as mostly a black box to be run in each 
% pooling region.
%
% This version applies a mask specifying the boundaries of the original
% stimulus on each global iteration. In our experience this speeds up
% convergence but does not appreciably change the results.
%
% Intended only to be called by synthesize_mongrel
% 
%   img: (Pre-padded) input image
%   componentInfo: Specifies conversion from RGB to ICA color space
%   maxscale: Has to do what what scales we are currently
%           synthesizing. 2-vector.
%   fixPt: [fx, fy] image coordinates of fixation
%   reconSeed: % Result so far, or [] if the first time through
%   imgMask: Shows where the original (pre-padding) image lay
%   Nor: P&S (Portilla & Simoncelli) number of orientations
%   Na: P&S number of taps of autocorrelation to use
%   nItersVec: How many iterations to do at this set of scales
%   foveaSize: Radius of the fovea, in pixels 
%   prIters: Number of local iterations for each pooling region
%   poolingRegions: List of precomputed pooling regions
%   savename: Name of file used to save intermediate results
%   saveIntermediateFlag: Save intermediate results if == 1
%   debugPath: Directory in which to save intermediate results
%
%   Returns synthesis for the current set of scales
%
% Copyright (C) 2019  Ruth Rosenholtz, Alvin Raj, & Lavanya Sharan
% See COPYING.txt for license details

%% 1. INITIAL SETUP
[H, W, D] = size(img)
stitch = "stitch"
%if D==3
%    colorSynth = 1;
%else
%    colorSynth = 0;
%end
[x, y] = meshgrid(1:W,1:H); 

%check if we are using a fovea
if isnan(fixPt(2))
    uniform_pooling = true;
else
    uniform_pooling = false;
end;

%create our seed image
foveaImage = zeros(size(img));

%fill foveal region of seed image with original image if this is a foveaed mongrel
if uniform_pooling==false
    % save fovea information
    % we presume that the "fovea" -- a region about fixation specified by
    % foveaSize, is preserved perfectly. On each global iteration we reapply
    % this known "fovea" information.
    % (somewhat wasteful, could just read in fovea image after the first time)
    foveaMask = (x-fixPt(1)).^2+(y-fixPt(2)).^2<foveaSize^2; 
    for i = 1:D
        foveaImage(:,:,i) = img(:,:,i).*foveaMask;
    end
    if saveIntermediateFlag>0 && maxscale(1)==1 % 1st time, save fovea image
        if ~exist(sprintf('%s%s',debugPath,savename),'dir')
            mkdir(sprintf('%s%s',debugPath,savename));
        end
        if maxscale(1)==1 && nItersVec(2)==1
            imwrite(foveaImage,sprintf('%s%s/fovea.png',debugPath,savename));
        end
    end
end

% set up seed image, i.e. what image the synthesis starts with 
nItersStop = nItersVec(1);
nItersStart = nItersVec(2);
if maxscale(1)==1 && nItersStart==1 % first pooling region at first scale
    if ~isempty(reconSeed)
        synth = reconSeed;
    else
        synth = foveaImage;
    end
else
    % if this is not your first time through this code, reconSeed will be 
    % set based on previous syntheses; use that
    synth = reconSeed;
end

% this reapplies the mask specifying the boundary of the original image
% used because useImgMask = 1 (currently hard-coded)
% seed mask is reapplied on first iteration and after every reapplyMask iterations
reapplyMask = 1; 

% initialize color seed
if colorSynth
    % convert img and foveaImage to ICA space
    % (somewhat wasteful, could do this once)
    img = separateImage(img, componentInfo); 
    synth = separateImage(synth, componentInfo);       
    foveaImage = separateImage(foveaImage, componentInfo); 
    if ~isempty(reconSeed)
        reconSeed = separateImage(reconSeed, componentInfo);  
    end
end
paddingColor = img(1,1,:); % paddingColor in ica space, if colorSynth. 
% image is already padded, so this is the actual padding color, not the desired color.

%% 2. GET TEXTURE DESCRIPTORS FOR EACH POOLING REGION
% (legacy code, wasteful, recomputing stats every time)
textureDescriptor = cell(size(poolingRegions,1),1);
globalMaxScale = maxscale(2); currMaxScale = maxscale(1);
largestShrink = globalMaxScale - 1;
for itp=1:length(poolingRegions)
    % pull out variables just to make the code easier to read
    buffPixels = poolingRegions(itp,5);            
    p_x = poolingRegions(itp,1);
    p_y = poolingRegions(itp,2);
    a = poolingRegions(itp,3);
    b = poolingRegions(itp,4);
    try
        p_y-a-buffPixels+1;
        p_y+b+buffPixels;
        p_x-a-buffPixels+1;
        p_x+b+buffPixels;
        size(img);
        impatch = img(p_y-a-buffPixels+1:p_y+b+buffPixels, p_x-a-buffPixels+1:p_x+b+buffPixels, :);
    catch
        keyboard
    end

    % resize to make the patch size a power of 2
    nearestpow2 = ceil(log2(size(impatch,1)));
    nearestpow2 = max([nearestpow2, largestShrink+3]);
        
    % resize to ensure coarse-to-fine processing works correctly
    impatch_sc = imresize(impatch,(2^nearestpow2)*[1,1],'bicubic');
    shrinkFactor = 2^(globalMaxScale-currMaxScale);
    if shrinkFactor>1
        impatch_sc = imresize(impatch_sc,1/shrinkFactor,'bicubic');
    end
    % get texture descriptor using P-S code
    for i=1:D,
        textureDescriptor{itp}.p(i) = textureAnalysis_forTTM(impatch_sc(:,:,i),currMaxScale,Nor,Na);
    end 
    if colorSynth % Also save stats about correlations between bands
        c123 = [reshape(impatch_sc(:,:,1),[],1) reshape(impatch_sc(:,:,2),[],1) reshape(impatch_sc(:,:,3),[],1)];
        textureDescriptor{itp}.cmu = mean(c123);
        textureDescriptor{itp}.cmat = cov(c123);
    end
    
end

% saving textureDescriptor (the feature vectors for each pooling region)
% only saving at 4th scale this is the only time that the call to textureAnalysis_for_TTM
% gets 4 scales as the scale parameter (which is all the scales to run)
if currMaxScale==4
    % get the output directory (remove the debug folder)
    path_parts = regexp(debugPath,'/','split');
    output_dir = sprintf('%s',path_parts{1});
    % if ~exist(sprintf('%s',debugPath),'dir')
    %     mkdir(sprintf('%s',debugPath));
    % end
    % save(sprintf('%s/poolingDescriptor.mat',output_dir),'textureDescriptor', 'poolingRegions');
    % save(sprintf('%s/parameter.mat',output_dir),'poolingRegions', '-append');
    % disp(sprintf('saved pooling region information'))
end

if saveIntermediateFlag>0
    if ~exist(sprintf('%s%s/intermediate',debugPath,savename),'dir')
        mkdir(sprintf('%s%s/intermediate',debugPath,savename));
    end
end

%% 3. SYNTHESIZE POOLING REGIONS 
for its=nItersStart:nItersStop
    disp(sprintf('Scale = %02d, Sweep = %02d',currMaxScale,its));
    listPR = [1:length(poolingRegions)];
    
    % check if this round is spiraling out or in
    % if mod(its,2)
    %     listPR = [1:length(poolingRegions)];
    % else
    %     listPR = [length(poolingRegions):-1:1];
    % end
    
    %shuffle randomly
    listPR = randperm(length(poolingRegions));

    % synthesize pooling regions one by one
    for itp=listPR
        buffPixels = poolingRegions(itp,5);  
        p_x = poolingRegions(itp,1);
        p_y = poolingRegions(itp,2);
        a = poolingRegions(itp,3);
        b = poolingRegions(itp,4);

        % get seed for synthesis of this patch, from previous
        % iterations over this or other pooling regions, if available
        if its==1 && currMaxScale==1 
            if ~isempty(reconSeed)  % If you're passed a seed the first time, use it
                seedpatch = reconSeed(p_y-a-buffPixels+1:p_y+b+buffPixels, p_x-a-buffPixels+1:p_x+b+buffPixels, :);
            else % if it's your first time through and you don't have a seed, use a random one
                % set random seed to fixed value
                seedpatch = randn(a+b+2*buffPixels, a+b+2*buffPixels, D);
            end
        else % if this is not your first time through, use what you know from previous syntheses
            seedpatch = synth(p_y-a-buffPixels+1:p_y+b+buffPixels, p_x-a-buffPixels+1:p_x+b+buffPixels, :);
        end
    
        % apply image mask, i.e. reimpose the boundaries of the original image
        skipPR = 0;
        if (mod(its-1,reapplyMask)==0)&&(~isempty(imgMask))
            mask = imgMask(p_y-a-buffPixels+1:p_y+b+buffPixels, p_x-a-buffPixels+1:p_x+b+buffPixels);
            for i=1:D,
                seedpatch(:,:,i) = (seedpatch(:,:,i).*mask)+(paddingColor(i)*(1-mask));
            end
            if sum(mask(:))==0,
                skipPR = 1;
                res = seedpatch;
            end
        end
    
        if ~skipPR,
            % resize to make the seed patch size a power of 2
            nearestpow2 = ceil(log2(size(seedpatch,1)));
            nearestpow2 = max([nearestpow2, largestShrink+3]);
            
            % resize to ensure coarse-to-fine processing works correctly
            seedpatch_sc = imresize(seedpatch,(2^nearestpow2)*[1,1],'bicubic');
            shrinkFactor = 2^(globalMaxScale-currMaxScale);
            if shrinkFactor>1
                seedpatch_sc = imresize(seedpatch_sc,1/shrinkFactor,'bicubic');
            end
            res_sc = [];
            for i=1:D,
%                 disp([its itp i]); % RRR DEBUG
                res_sc(:,:,i) = textureSynthesis_forTTM(textureDescriptor{itp}.p(i), seedpatch_sc(:,:,i), prIters, verbose);
            end
            
            [rh, rw, d] = size(res_sc);
            
            if colorSynth
                data = [reshape(res_sc(:,:,1),[],1) reshape(res_sc(:,:,2),[],1) reshape(res_sc(:,:,3),[],1)];
                mudata = mean(data);
                % zero mean it first
                data = data - repmat(mudata, [size(data,1), 1]);
                % fit the covariance next
                data = fitCovariance(data, textureDescriptor{itp}.cmat);
                % add the desired mean next
                data = data + repmat(textureDescriptor{itp}.cmu, [size(data,1), 1]);
                for i=1:D,
                    res_sc(:,:,i) = reshape(data(:,i),[rh rw]);
                end
            end
            
            % undo resize so synthesized region fits back into image
            res = imresize(res_sc,(a+b+2*buffPixels)*[1,1],'bicubic');
        end

        % put synthesized region back into image
        % make a feathered mask
        % scale-dependent smoothing
        fEdge = max(3,ceil((a+b)/(2^(1+currMaxScale))));
        prMask = zeros(a+b+2*buffPixels);  
        prMask(buffPixels+1:end-buffPixels,buffPixels+1:end-buffPixels) = 1;
        filt = fspecial('gaussian', [fEdge, fEdge], fEdge/3);
        prMask = imfilter(prMask, filt);

        % paste synthesis back using feathered mask
        for i=1:D,
            synth(p_y-a-buffPixels+1:p_y+b+buffPixels, p_x-a-buffPixels+1:p_x+b+buffPixels, i)...
                = (1-prMask).*synth(p_y-a-buffPixels+1:p_y+b+buffPixels, p_x-a-buffPixels+1:p_x+b+buffPixels, i)... % Outside the mask, fill in what was there before
                + prMask.*res(:,:,i); % The new stuff
        end
                
        if verbose
            disp(sprintf('Scale = %02d, Sweep = %02d, Pooling Region = %04d',currMaxScale,its,itp));
        end        
    end % move on to the next pooling region
    
    if uniform_pooling==false
        % apply fovea to result if this is a foveated image
        for i=1:D,
            synth(:,:,i) = synth(:,:,i).*(1-foveaMask) + foveaImage(:,:,i).*foveaMask;
        end
    end
    
    % save intermediate results
    if saveIntermediateFlag>0
        if colorSynth
            % convert back to rgb
            synthRGB = integrateImage(synth, componentInfo, H, W);        
            imwrite(synthRGB,sprintf('%s%s/intermediate/result_scale_%02d_sweep%02d.png',debugPath,savename,currMaxScale,its));   
        else
            imwrite(synth,sprintf('%s%s/intermediate/result_scale_%02d_sweep%02d.png',debugPath,savename,currMaxScale,its));   
        end
    end
end % do another iteration

if colorSynth
    % convert back to rgb
    synth = integrateImage(synth, componentInfo, H, W);  
end
