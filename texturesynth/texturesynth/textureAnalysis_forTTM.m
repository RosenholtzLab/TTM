function [params] = textureAnalysis_forTTM(im0, Nsc, Nor, Na)

% Analyze texture for application of Portilla-Simoncelli model/algorithm.
%
% [params] = textureAnalysis(im0, Nsc, Nor, Na);
% 	im0: 	original image
% 	Nsc: 	number of scales
% 	Nor: 	number of orientations
% 	Na:	spatial neighborhood considered (Na x Na)
%
% Example: Nsc=4; Nor=4; Na=7;
%
% See also textureSynthesis.

% Javier Portilla and Eero Simoncelli.
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
% Slightly adapted for TTM, to minimize crashes on psychophysical stimuli
% that P&S was not built for.


Warn = 0;  % Set to 1 if you want to see warning messages

% Check required args are passed
if (nargin < 4)
    error('Function called with too few input arguments');
end

% 1D interpolation filter, for scale cross-correlations:
interp = [-1/16 0 9/16 1 9/16 0 -1/16]/sqrt(2);

if ( mod(Na,2) == 0 )
    error('Na is not an odd integer');
end

% If the spatial neighborhood Na is too big for the lower scales,
% "modacor22.m" will make it as big as the spatial support at
% each scale:

[Ny,Nx] = size(im0);
nth = log2(min(Ny,Nx)/Na);
if nth<Nsc & Warn,
    fprintf(1,'Warning: Na will be cut off for levels above #%d !\n', floor(nth+1));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

la = floor((Na-1)/2);

bins = 256;

% Pixel statistics
[mn0 mx0] = range2(im0);
diff = mx0-mn0;
if (diff < 1e-7)
    params.isConstant = 1;
    params.constval = mn0;
    return;
else
    params.isConstant = 0;
    params.constval = NaN;
end

[mean0 var0 skew0 kurt0] = marginalStats_forTTM(im0);
hist0 = computeHistogram_forTTM(im0, bins);
statg0 = [mean0 var0 skew0 kurt0 mn0 mx0];


% Build the steerable pyramid
[pyr0,pind0] = buildSCFpyr(im0,Nsc,Nor-1);

if ( any(vector(mod(pind0,2))) )
    error('Algorithm will fail: Some bands have odd dimensions!');
end

% Subtract mean of lowBand:
nband = size(pind0,1);
% remove mean from whole band
pyr0(pyrBandIndices(pind0,nband)) = ...
    real(pyrBand(pyr0,pind0,nband)) - mean2(real(pyrBand(pyr0,pind0,nband)));

rpyr0 = real(pyr0);
apyr0 = abs(pyr0);

% Subtract mean of magnitude:
magMeans0 = zeros(size(pind0,1), 1);
for nband = 1:size(pind0,1)
    indices = pyrBandIndices(pind0,nband);
    magMeans0(nband) = mean2(apyr0(indices));
    apyr0(indices) = apyr0(indices) - magMeans0(nband);
end

% Compute central autoCorr of lowband
acr = NaN * ones(Na,Na,Nsc+1);
nband = size(pind0,1);
ch = pyrBand(pyr0,pind0,nband);

[mpyr,mpind] = buildSFpyr(real(ch),0,0);
im = pyrBand(mpyr,mpind,2);

acr(:,:,Nsc+1) = computeWeightedAutocorrelation_mex(im, Na, ones(size(im)), 1);

mui = mean2(im);
vari = var_forTTM(im);

skew0p = zeros(Nsc+1,1);
kurt0p = zeros(Nsc+1,1);
hist0p = cell(Nsc+1,1);
if vari/var0 > 1e-6,
    skew0p(Nsc+1) = skew2_forTTM(im, mui, vari);
    kurt0p(Nsc+1) = kurt2_forTTM(im, mui, vari);
    hist0p{Nsc+1} = computeHistogram_forTTM(im, bins);
else
    skew0p(Nsc+1) = 0;
    kurt0p(Nsc+1) = 3;
end

% Compute  central autoCorr of each Mag band, and the autoCorr of the
% combined (non-oriented) band.
ace = NaN * ones(Na,Na,Nsc,Nor);
for nsc = Nsc:-1:1,
    for nor = 1:Nor,
        
        nband = (nsc-1)*Nor+nor+1;
        ch = pyrBand(apyr0,pind0,nband);
        rmsk = ones(size(ch));
        ace(:,:,nsc, nor) = computeWeightedAutocorrelation_mex(ch, Na, rmsk,1);
    end
    
    % Combine ori bands
        
    ch = combineOrientationBands(nsc, Nor, rpyr0, pind0);
    
    % add it to the low pass band
    im = real(expand(im,2))/4;
    im = im + ch;
    
    rmsk = ones(size(im));
    acr(:,:,nsc) = computeWeightedAutocorrelation_mex(im, Na, rmsk,1);
        
    mui = mean2(im);
    vari = var_forTTM(im);
    
    if vari/var0 > 1e-6,
        skew0p(nsc) = skew2_forTTM(im, mui, vari);
        kurt0p(nsc) = kurt2_forTTM(im, mui, vari);
        hist0p{nsc} = computeHistogram_forTTM(im, bins);
    else
        skew0p(nsc) = 0;
        kurt0p(nsc) = 3;
        hist0p{nsc} = [];    
    end
    
end

% Compute the cross-correlation matrices of the coefficient magnitudes
% pyramid at the different levels and orientations

C0 = zeros(Nor,Nor,Nsc+1);
Cx0 = zeros(Nor,Nor,Nsc);

Cr0 = zeros(2*Nor,2*Nor,Nsc+1);
Crx0 = zeros(2*Nor,2*Nor,Nsc);

for nsc = 1:Nsc,
    firstBnum = (nsc-1)*Nor+2;
    cousinSz = prod(pind0(firstBnum,:));
    ind = pyrBandIndices(pind0,firstBnum);
    cousinInd = ind(1) + (0:Nor*cousinSz-1);
    
    if (nsc<Nsc)
        parents = zeros(cousinSz,Nor);
        rparents = zeros(cousinSz,Nor*2);
        for nor=1:Nor,
            nband = (nsc-1+1)*Nor+nor+1;
            
            tmp = expand(pyrBand(pyr0, pind0, nband),2)/4;
            rtmp = real(tmp); itmp = imag(tmp);
            % Double phase:
            tmp = sqrt(rtmp.^2 + itmp.^2) .* exp(2 * sqrt(-1) * atan2(rtmp,itmp));
            rparents(:,nor) = vector(real(tmp));
            rparents(:,Nor+nor) = vector(imag(tmp));
            
            tmp = abs(tmp);
            parents(:,nor) = vector(tmp - mean2(tmp));
        end
    else
        tmp = real(expand(pyrLow(rpyr0,pind0),2))/4;
        rparents = [vector(tmp),...
            vector(shift(tmp,[0 1])), vector(shift(tmp,[0 -1])), ...
            vector(shift(tmp,[1 0])), vector(shift(tmp,[-1 0]))];
        parents = [];
    end
    cousins = reshape(apyr0(cousinInd), [cousinSz Nor]);
    
    cousins = reshape(apyr0(cousinInd), [cousinSz Nor]);

    nc = size(cousins,2);   np = size(parents,2);
    
    C0(1:nc,1:nc,nsc) = innerProd(cousins)/cousinSz;
    if (np > 0)
        Cx0(1:nc,1:np,nsc) = (cousins'*parents)/cousinSz;
        if (nsc==Nsc)
            C0(1:np,1:np,Nsc+1) = innerProd(parents)/(cousinSz/4);
        end
    end
    
    cousins = reshape(real(pyr0(cousinInd)), [cousinSz Nor]);
    nrc = size(cousins,2);   nrp = size(rparents,2);
    Cr0(1:nrc,1:nrc,nsc) = innerProd(cousins)/cousinSz;
    if (nrp > 0)
        Crx0(1:nrc,1:nrp,nsc) = (cousins'*rparents)/cousinSz;
        if (nsc==Nsc)
            Cr0(1:nrp,1:nrp,Nsc+1) = innerProd(rparents)/(cousinSz/4);
        end
    end
end

% Calculate the mean, range and variance of the LF and HF residuals' energy.

channel = pyr0(pyrBandIndices(pind0,1));
vHPR0 = var_forTTM(channel);

statsLPim = [skew0p kurt0p];

params.pixelStats = statg0;
params.pixelLPStats = statsLPim;
params.autoCorrReal = acr;
params.autoCorrMag = ace;
params.magMeans = magMeans0;
params.cousinMagCorr = C0;
params.parentMagCorr = Cx0;
params.cousinRealCorr = Cr0;
params.parentRealCorr = Crx0;
params.varianceHPR = vHPR0;

params.pixelHist = hist0;
params.pixelLPHist = hist0p;
