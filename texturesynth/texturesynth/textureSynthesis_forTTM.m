function [im,snrP,imS] = textureSynthesis_forTTM(params, im0, Niter, verbose, cmask, displayFlag)
% [res,snrP,imS] = textureSynthesis_forTTM(params, initialIm, Niter, verbose, cmask, displayFlag)
%
% Synthesize texture using a slightly modified Portilla-Simoncelli model/algorithm.
%
% params: structure containing texture parameters (as returned by textureAnalysis).
%
% im0: initial image, OR a vector (Ydim, Xdim, [SEED]) containing
% dimensions of desired image and an optional seed for the random
% number generator.  If dimensions are passed, initial image is
% Gaussian white noise.
%
% Niter (optional): Number of iterations.  Default = 50.
%
% cmask (optional): binary column vector (4x1) indicating which sets of
% constraints we want to apply in the synthesis. The four sets are:
%               1) Marginal statistics (mean, var, skew, kurt, range)
%               2) Correlation of subbands (space, orientation, scale)
%               3) Correlation of magnitude responses (sp, or, sc)
%               4) Relative local phase
%
% snrP (optional):	Set of adjustment values (in dB) of the parameters.
% imS (optional):	Sequence of synthetic images, from niter = 1 to 2^n, being
% 			n = floor(log2(Niter)).

% Originally by Javier Portilla and Eero Simoncelli.
% Work described in:
%  "A Parametric Texture Model based on Joint Statistics of Complex Wavelet Coefficients".
%  J Portilla and E P Simoncelli. Int'l Journal of Computer Vision,
%  vol.40(1), pp. 49-71, Dec 2000.
%
% Please refer to this publication if you use the program for research or
% for technical applications. Thank you.
%
% Copyright, Center for Neural Science, New York University, January 2001.
% All rights reserved.
%
% LS: Added displayFlag to allow visualization of intermediate results
% Various other changes made by the Rosenholtz group, most non-substantive,
% bug fixes, etc. Nonetheless, please use this version for TTM, not the
% original. 


Warn = 0;  % Set to 1 if you want to see warning messages

if (params.isConstant)
    im = zeros(size(im0)) + params.constval;
    imS = im;
    return;
end

% Check required args are passed:
if (nargin < 2)
    error('Function called with too few input arguments');
end

if ( ~exist('Niter','var') | isempty(Niter) )
    Niter = 50;
end

if (exist('cmask','var') & ~isempty(cmask) )
    cmask = (cmask > 0.5);  % indices of ones in mask
else
    cmask = ones(4,1);
end

if ~exist('displayFlag','var') | isempty(displayFlag)
    displayFlag = 0;
end

% Extract parameters
statg0 = params.pixelStats;
mean0 = statg0(1); var0 =  statg0(2);
skew0 =  statg0(3); kurt0 =  statg0(4);
mn0 =  statg0(5);  mx0 = statg0(6);
statsLPim = params.pixelLPStats;
skew0p = statsLPim(:,1);
kurt0p = statsLPim(:,2);
vHPR0 = params.varianceHPR;
acr0 = params.autoCorrReal;
ace0 = params.autoCorrMag;
magMeans0 = params.magMeans;
C0 = params.cousinMagCorr;
Cx0 = params.parentMagCorr;
Crx0 = params.parentRealCorr;

% Extract {Nsc, Nor, Na} from params
tmp = size(params.autoCorrMag);
Na = tmp(1); Nsc = tmp(3);
Nor = tmp(length(tmp))*(length(tmp)==4) + (length(tmp)<4);
la = (Na-1)/2;

% If im0 is a vector of length 2, create Gaussian white noise image of this
% size, with desired pixel mean and variance.  If vector length is
% 3,  use the 3rd element to seed the random number generator.
if ( length(im0) <= 3 )
    if ( length(im0) == 3)
        randn('state', im0(3)); % Reset Seed
        im0 = im0(1:2);
    end
    im = mean0 + sqrt(var0)*randn(im0);
else
    im = im0;
end

% RRR DEBUG
% foo = isnan(im);
% if sum(foo(:))>0 
%     disp('warning, already have NaNs at start of texture synth of this patch');
% end

% If the spatial neighborhood Na is too big for the lower scales,
% "modacor22_forTTM.m" will make it as big as the spatial support at
% each scale:
[Ny,Nx] = size(im);
nth = log2(min(Ny,Nx)/Na);
if nth<Nsc+1 & Warn,
    fprintf(1,'Warning: Na will be cut off for levels above #%d !\n',floor(nth));
end

if displayFlag
    imf = max(1,gcf-1); snrf = imf+1;
    figure(imf);  clf
    subplot(1,2,1); grayRange = showIm(im,'auto',1); title('Starting image');
    drawnow
end

prev_im=im;

if displayFlag
    imf = max(1,gcf-1);
    figure(imf);
    clf;showIm(im,'auto',1); title(sprintf('iteration 0'));
end

nq = 0;
Nq = floor(log2(Niter));
imS = zeros(Ny,Nx,Nq);
snr1 = zeros(Niter,Nsc*Nor+2);
snr2 = zeros(Niter,Nsc+1);
snr3 = zeros(Niter,Nsc);
snr4 = zeros(Niter,Nsc+1);
snr4r = zeros(Niter,Nsc);
snr6 = zeros(Niter);
snr7 = zeros(Niter,2*(Nsc+1)+4);

% MAIN LOOP
for niter = 1:Niter
%     disp(niter); % RRR DEBUG
    p = 1;
    
    % Build the steerable pyramid
    [pyr,pind] = buildSCFpyr(im,Nsc,Nor-1);
    
    if ( any(vector(mod(pind,4))) )
        error('Algorithm will fail: band dimensions are not all multiples of 4!');
    end
    
    % Subtract mean of lowBand:
    nband = size(pind,1);
    pyr(pyrBandIndices(pind,nband)) = ...
        pyrBand(pyr,pind,nband) - mean2(pyrBand(pyr,pind,nband));
    
    apyr = abs(pyr);
    
    % Adjust autoCorr of lowBand
    
    nband = size(pind,1);
    ch = pyrBand(pyr,pind,nband);
    Sch = min(size(ch)/2);
    nz = sum(sum(~isnan(acr0(:,:,Nsc+1))));
    lz = (sqrt(nz)-1)/2;
    le = min(Sch/2-1,lz);
    im = real(ch);  %Reconstructed image: initialize to lowband
    [mpyr,mpind] = buildSFpyr(im,0,0);
    im = pyrBand(mpyr,mpind,2);
    % RRR DEBUG
%     foo = isnan(im);
%     if sum(foo(:))>0 
%         disp('warning, NaNs after autoCorr of lowBand');
%     end
    
    vari = var2(im);
    
    if cmask(2),
        
        if vari/var0 > 1e-4 && vari ~= 0,
            acr_nonan = getBiggestNonNanSquare(acr0(:,:,Nsc+1), min(size(im))); 
%             acr_nonan = acr0(:,:,Nsc+1);
            if (sum(isnan(im(:)))>0 || sum(isnan(acr_nonan(:)))>0) && verbose,
                disp(sprintf('First call to modacorr22_forTTM: NaNs in im (1st argument) = %d',sum(isnan(im(:)))));
                disp(sprintf('First call to modacorr22_forTTM: NaNs in acr_nonan (1st argument) = %d',sum(isnan(acr_nonan(:)))));
            end
            if verbose,
                disp(sprintf('First call to modacorr22, vari=%.2f, var0=%.2f',vari, var0));
            end
            [im, snr2(niter,Nsc+1)] = ...
                modacor22_forTTM(im, acr_nonan);
            % RRR DEBUG
%             foo = isnan(im);
%             if sum(foo(:))>0 
%                 disp('warning, NaNs after modacor22');
%             end

        end
        im = real(im);
    end % cmask(2)

    % RRR CHANGE
%     if cmask(1),
%         if vari/var0 > 1e-4 && vari ~= 0 && verbose 
%             if verbose,
%                 disp(sprintf('Imposing histogram, vari=%.2f, var0=%.2f',vari, var0));
%             end
%             [im,hflag,snr_sk,snr_kt] = imposeHistogram(params.pixelLPHist{Nsc+1}, im, verbose, 0);
%             % RRR debugging
%             foo = isnan(im);
%             if sum(foo(:))>0 
%                 disp('warning, NaNs after imposeHistogram');
%             end
% 
%             snr7(niter,2*(Nsc+1)-1) = snr_sk; snr7(niter,2*(Nsc+1)) = snr_kt;
%         end
%     end	% cmask(1)
    
    % Subtract mean of magnitude
    if cmask(3),
        magMeans = zeros(size(pind,1), 1);
        for nband = 1:size(pind,1)
            indices = pyrBandIndices(pind,nband);
            magMeans(nband) = mean2(apyr(indices));
            apyr(indices) = apyr(indices) - magMeans(nband);
        end
    end	% cmask(3)
    
    % Coarse-to-fine loop:
    for nsc = Nsc:-1:1
        
        firstBnum = (nsc-1)*Nor+2;
        cousinSz = prod(pind(firstBnum,:));
        ind = pyrBandIndices(pind,firstBnum);
        cousinInd = ind(1) + [0:Nor*cousinSz-1];
        
        % Interpolate parents
        parents = zeros(cousinSz,Nor);
        rparents = zeros(cousinSz,Nor*2);
        
        if (cmask(3) | cmask(4)),
            if (nsc<Nsc)
                for nor = 1:Nor
                    nband = (nsc+1-1)*Nor+nor+1;
                    
                    tmp = expand(pyrBand(pyr, pind, nband),2)/4;
                    rtmp = real(tmp); itmp = imag(tmp);
                    tmp = sqrt(rtmp.^2 + itmp.^2) .* exp(2 * sqrt(-1) * atan2(rtmp,itmp));
                    rparents(:,nor) = vector(real(tmp));
                    rparents(:,Nor+nor) = vector(imag(tmp));
                    
                    tmp = abs(tmp);
                    parents(:,nor) = vector(tmp - mean2(tmp));
                end
            else
                rparents = [];
                parents = [];
            end
        end % if (cmask(3) | cmask(4))
        
        if cmask(3),
            % Adjust cross-correlation with MAGNITUDES at other orientations/scales:
            cousins = reshape(apyr(cousinInd), [cousinSz Nor]);
            nc = size(cousins,2);   np = size(parents,2);
            if (np == 0)
                [cousins, snr3(niter,nsc)] = adjustCorr1s_forTTM(cousins,C0(1:nc,1:nc,nsc),verbose,2);
                % RRR DEBUG
%                 foo = isnan(cousins);
%                 if sum(foo(:))>0 
%                     disp('warning, NaNs after adjustCorr1s, coarse-to-fine');
%                 end

            else
                [cousins, snr3(niter,nsc), snr4(niter,nsc)] = ...
                    adjustCorr2s_forTTM(cousins, C0(1:nc,1:nc,nsc), parents, Cx0(1:nc,1:np,nsc), 3);
                % RRR DEBUG
%                 foo = isnan(cousins);
%                 if sum(foo(:))>0 
%                     disp('warning, NaNs after adjustCorr2s, coarse-to-fine');
%                 end

            end
            cousins = real(cousins);
            ind = cousinInd;
            apyr(ind) = vector(cousins);
            
            % Adjust autoCorr of mag responses
            nband = (nsc-1)*Nor+2;
            Sch = min(pind(nband,:)/2);
            nz = sum(sum(~isnan(ace0(:,:,nsc,1))));
            lz = (sqrt(nz)-1)/2;
            le = min(Sch/2-1,lz);
            for nor = 1:Nor,
                nband = (nsc-1)*Nor+nor+1;
                ch = pyrBand(apyr,pind,nband);
                acr_nonan = ace0(:,:,nsc,nor);
                if (sum(isnan(im(:)))>0 || sum(isnan(acr_nonan(:)))>0) && verbose,
                    disp(sprintf('Second call to modacor22_forTTM: NaNs in im (1st argument) = %d, scale = %d',sum(isnan(im(:))), nsc));
                    disp(sprintf('Second call to modacor22_forTTM: NaNs in acr_nonan (1st argument) = %d',sum(isnan(acr_nonan(:)))));
                end
                
                [ch, snr1(niter,nband-1)] = ...
                    modacor22_forTTM(ch, acr_nonan);
                % RRR DEBUG
%                 foo = isnan(ch);
%                 if sum(foo(:))>0 
%                     disp('warning, NaNs after 2nd call to modacor22');
%                 end

                ch = real(ch);
                ind = pyrBandIndices(pind,nband);
                apyr(ind) = ch;
                % Impose magnitude:
                mag = apyr(ind) + magMeans0(nband);
                mag = mag .* (mag>0);
                pyr(ind) = pyr(ind) .* (mag./(abs(pyr(ind))+(abs(pyr(ind))<eps)));
            end
        end   % cmask(3)
        
        % Adjust cross-correlation of REAL PARTS at other orientations/scales:
        cousins = reshape(real(pyr(cousinInd)), [cousinSz Nor]);
        Nrc = size(cousins,2);  Nrp = size(rparents,2);
        if  cmask(4) & (Nrp ~= 0)
            a3 = 0; a4 = 0;
            for nrc = 1:Nrc,
                cou = cousins(:,nrc);
                [cou, s3, s4] = ...
                    adjustCorr2s_forTTM(cou,mean(cou.^2),rparents,Crx0(nrc,1:Nrp,nsc), 3);
                % RRR DEBUG
%                 foo = isnan(cou);
%                 if sum(foo(:))>0 
%                     disp('warning, NaNs after 2nd call to adjustCorr2s');
%                 end

                a3 = s3 + a3;
                a4 = s4 + a4;
                cousins(:,nrc) = cou;
            end
            snr4r(niter,nsc) = a4/Nrc;
        end
        % always update the pyramid.
        pyr(cousinInd) = vector(cousins(1:Nor*cousinSz));
                
        % Re-create analytic subbands
        dims = pind(firstBnum,:);
        ctr = ceil((dims+0.5)/2);
        ang = mkAngle(dims, 0, ctr);
        ang(ctr(1),ctr(2)) = -pi/2;
        for nor = 1:Nor,
            nband = (nsc-1)*Nor+nor+1;
            ind = pyrBandIndices(pind,nband);
            ch = pyrBand(pyr, pind, nband);
            ang0 = pi*(nor-1)/Nor;
            xang = mod(ang-ang0+pi, 2*pi) - pi;
            amask = 2*(abs(xang) < pi/2) + (abs(xang) == pi/2);
            amask(ctr(1),ctr(2)) = 1;
            amask(:,1) = 1;
            amask(1,:) = 1;
            amask = fftshift(amask);
            ch = ifft2(amask.*fft2(ch));	% "Analytic" version
            pyr(ind) = ch;
        end
        
        % Combine ori bands
        bandNums = [1:Nor] + (nsc-1)*Nor+1;  %ori bands only
        ind1 = pyrBandIndices(pind, bandNums(1));
        indN = pyrBandIndices(pind, bandNums(Nor));
        bandInds = [ind1(1):indN(length(indN))];
        % Make fake pyramid, containing dummy hi, ori, lo
        fakePind = pind([bandNums(1), bandNums, bandNums(Nor)+1],:);
        fakePyr = [zeros(prod(fakePind(1,:)),1);...
            real(pyr(bandInds)); zeros(prod(fakePind(size(fakePind,1),:)),1)];
        ch = reconSFpyr(fakePyr, fakePind, [1]);     % recon ori bands only
        im = real(expand(im,2))/4;
        im = im + ch;
        vari =  acr0(la+1:la+1,la+1:la+1,nsc);
        if cmask(2),
            if vari/var0 > 1e-4,
                acr_nonan = acr0(:,:,nsc);
                if (sum(isnan(im(:)))>0 || sum(isnan(acr_nonan(:)))>0) && verbose,
                    disp(sprintf('Third call to modacor22_forTTM: NaNs in im (1st argument) = %d, scale = %d',sum(isnan(im(:))),nsc));
                    disp(sprintf('Third call to modacor22_forTTM: NaNs in acr_nonan (1st argument) = %d',sum(isnan(acr_nonan(:)))));
                end
                
                [im, snr2(niter,nsc)] = ...
                    modacor22_forTTM(im, acr_nonan);
                % RRR DEBUG
%                 foo = isnan(im);
%                 if sum(foo(:))>0 
%                     disp('warning, NaNs after 3rd call to modacor22');
%                 end

            end
        end	% cmask(2)
        im = real(im);
        
        % RRR CHANGE
%         if cmask(1),
%             % Fix marginal stats
%             if vari/var0 > 1e-4 && vari ~= 0 && verbose,
%                 % [im,snr7(niter,2*nsc-1)] = modskew(im,skew0p(nsc),p);       % Adjusts skewness
%                 % [im,snr7(niter,2*nsc)] = modkurt(im,kurt0p(nsc),p);         % Adjusts kurtosis
%                 [im,hflag,snr_sk,snr_kt] = imposeHistogram(params.pixelLPHist{nsc}, im, verbose, 0);
%                 % RRR DEBUG
% %                 foo = isnan(im);
% %                 if sum(foo(:))>0 
% %                     disp('warning, NaNs after 2nd call to imposeHistogram');
% %                 end
% 
%                 snr7(niter,2*nsc-1) = snr_sk; snr7(niter,2*nsc) = snr_kt;
%             end
%         end	% cmask(1)
        
    end  %END Coarse-to-fine loop
    
    % Adjust variance in HP, if higher than desired
    if (cmask(2)|cmask(3)|cmask(4)),
        ind = pyrBandIndices(pind,1);
        ch = pyr(ind);
        vHPR = mean2(ch.^2);
        if vHPR > vHPR0,
            ch = ch * sqrt(vHPR0/vHPR);
            pyr(ind) = ch;
        end
    end % cmask
    im = im + reconSFpyr(real(pyr), pind, [0]);  %recon hi only
    
    % Pixel statistics
    means = mean2(im);
    vars = var_forTTM(im);
    snr7(niter,2*(Nsc+1)+1) = snr(var0,var0-vars);
    im = im-means;  			% Adjusts mean and variance
    [mns mxs] = range2(im + mean0);
    snr7(niter,2*(Nsc+1)+2) = snr(mx0-mn0,sqrt((mx0-mxs)^2+(mn0-mns)^2));
    if cmask(1),
        im = im*sqrt(((1-p)*vars + p*var0)/vars);
    end	% cmask(1)
    im = im+mean0;
    if cmask(1),
        %         [im, snr7(niter,2*(Nsc+1)+3)] = modskew(im,skew0,p);	% Adjusts skewness (keep mean and variance)
        %         [im, snr7(niter,2*(Nsc+1)+4)] = modkurt(im,kurt0,p);	% Adjusts kurtosis (keep mean and variance,
        [im,hflag,snr_sk,snr_kt] = imposeHistogram(params.pixelHist, im, verbose, 0);
        snr7(niter,2*(Nsc+1)+3) = snr_sk; snr7(niter,2*(Nsc+1)+4) = snr_kt;
        im = max(min(im,(1-p)*max(max(im))+p*mx0),...
            (1-p)*min(min(im))+p*mn0);		% Adjusts range (affects everything)
    else
        snr7(niter,2*(Nsc+1)+3) = snr(skew0,skew0-skew2(im));
        snr7(niter,2*(Nsc+1)+4) = snr(kurt0,kurt0-kurt2(im));
    end	% cmask(1)
    % RRR DEBUG
%     foo = isnan(cousins);
%     if sum(foo(:))>0 
%         disp('warning, NaNs after adjusting pixel stats');
%     end
    
    snr6(niter,1) = snr(im-mean0,im-prev_im);
    
    if floor(log2(niter))==log2(niter),
        nq = nq + 1;
        imS(:,:,nq) = im;
    end
    
    im = real(im);
    
    tmp = prev_im;
    prev_im=im;
    
    if displayFlag
        imf = max(1,gcf-1);
        figure(imf);
        clf;showIm(im,'auto',1); title(sprintf('iteration %02d',niter));
    end
    
    
    % accelerator
    alpha = 0.8;
    im = im + alpha*(im - tmp);
    
end


im = prev_im;
if (any(imag(im) ~= 0))
    keyboard;
end

snrP = struct();
snrP.snr1 = snr1;
snrP.snr2 = snr2;
snrP.snr3 = snr3;
snrP.snr4 = snr4;
snrP.snr4r = snr4r;
snrP.snr6 = snr6;
snrP.snr7 = snr7;

return;
