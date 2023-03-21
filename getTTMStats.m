function stats = getTTMStats(img_file,fx,fy,foveaSize,job_file,output_folder)

    %read in job file stuff

    addpath(genpath('./'));
    if nargin < 7
        job_file = 'default_fulliter_60olap.job';
        output_folder = './'
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
    % get rid of strcat because we are reading it from the textfile now
    % out_dir = strcat(output_folder,'/ecc_',num2str(round(fx)));
    out_dir = output_folder;
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
    im_pad_path = strcat(out_dir, '/pad_',imname,'.png');
    img_mask_path = strcat(out_dir, '/imgmask',imname,'.png');
    output_parameter_path = strcat(out_dir, '/synthesis_parameter_',imname,'.txt'); % for storing parameters.

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
            disp('Random seed set incorrectly. Use random_mode "shuffle" or define random_type and random_seed in jobfile.');
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

    end

    disp('finished peripheral blur')

    % since some old routines grab the padding color from the corner, reinstate
    % that here, after the blurring operation which might have changed it
    im_pad(1,1,:) = padding_color; 
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
    prm.source_img_fixation_x = x;
    prm.source_img_fixation_y = y;
    prm.curr_random_seed = curr_random_seed;
    prm.componentInfo = componentInfo;
    prm.poolingRegions = poolingRegions;
    prm.useImgMask = useImgMask;

    disp('saving random seeds')
    %check if file is open and wait until not open to use
    filename = strcat(output_folder,"/imgname_ecc_randomseeds_",num2str(prm.source_img_fixation_x),".txt");
    fid = fopen(filename,'a');
    fwrite(fid,sprintf('%s\t%s\t%d\n', prm.image_name,num2str(round(fx)),prm.curr_random_seed.Seed));
    fclose(fid);

    %% 7. RUN TTM
    % run the job
    % "stitch" together multiple pooling regions to create a continuous mongrel 
    % consistent with all measured statistics (approximately, of course) 
    reconsb = [];
    %maxscale = prm.scalesToRun(length(prm.scalesToRun)); % lower numbers are coarser scales
    %nIterStart = 1;
    
    %convert to stitch bare bones variable names
    img = im_pad;
    maxscale = [1, max(prm.scalesToRun)];
    reconSeed = reconsb;
    imgMask = maskimg;
    Nor = prm.Nor;
    Na = prm.Na;
    poolingRegions = prm.poolingRegions;
    savename = prm.savename;
    colorSynth = prm.colorSynth;
    out_path = output_folder
    output_dir = out_path

    %%%%%%% Call to StitchBareBones

    disp('calling stitchbarebones')
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
    synth = foveaImage;
    
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
        %itp
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


    %save pooling descriptor and stats
    save(sprintf('%s/poolingDescriptor.mat',output_dir),'textureDescriptor', 'poolingRegions');
    save(sprintf('%s/poolingDescriptor.mat',output_dir),'poolingRegions', '-append');
    disp(sprintf('saved pooling region information'))
    disp('all done')

end

