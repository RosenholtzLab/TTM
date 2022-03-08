function [Y,snrV,Chf]=modacor22_forTTM(X,Cy)
% [Y,snr,Chf]=modacor22(X,Cy,p);
% 
% Imposes the desired autocorrelation in the given (central) samples (Cy) to
% an image X, convolving it with an even filter of size(Cy), in such a way
% that the image containts change as less as possible, in a LSE sense.
%	Chf:	Fourier transform of the filter that forces the autocorrelation
%	p [OPTIONAL]:	mixing proportion between Cx and Cy
%			it imposes (1-p)*Cx + p*Cy,
%			being Cx the actual autocorrelation.
%			DEFAULT: p = 1; 

% JPM, 10/97, working with EPS, NYU
% RR group may have modified slightly to fix crashes. Please use this
% version for TTM syntheses.

Warn = 0;  % Set to 1 if you want to see warning messages
p = 1;

C0 = Cy;

% Compute the autocorrelation function of the original image

[Ny,Nx]=size(X);
Nc=size(Cy,1); 	% Normally Nc<<Nx, only the low indices of the autocorrelation	
if (2*Nc-1 > Nx) & Warn
  warning('Autocorrelation neighborhood too large for image: reducing');
  Nc = 2*floor(Nx/4)-1;
  first = (size(Cy,1)-Nc)/2;
  Cy = Cy(first+1:first+Nc, first+1:first+Nc);
end

mask = ones(size(X));
[Cx, NormX]  = computeWeightedAutocorrelation_mex(X, 2*Nc-1, mask, 0);
if sum(isnan(Cx(:)))>0, 
    disp(sprintf('NaNs in Cx (after computeweightedautocorr call) = %d',sum(isnan(Cx(:)))));
end
    
Lc=(Nc-1)/2;
try
Cx_sy = Cx(Nc + (-Lc:Lc), Nc + (-Lc:Lc));
catch
    keyboard
end
% Unnormalize the previously normalized correlation
Cy = Cy .* NormX(Nc + (-Lc:Lc), Nc + (-Lc:Lc));

cy=Ny/2+1;
cx=Nx/2+1;

Cy0 = Cy;
Cy = p*Cy + (1-p)*Cx_sy;

% Compare the actual correlation with the desired one
snrV=10*log10(sum(sum(Cy0.^2))/sum(sum((Cy0-Cx_sy).^2)));

% Build the matrix that performs the convolution Cy1=Tcx*Ch1

Ncx=4*Lc+1;
M=(Nc^2+1)/2;
Tcx=zeros(M);

for i=Lc+1:2*Lc,
	for j=Lc+1:3*Lc+1,
		nm=(i-Lc-1)*(2*Lc+1)+j-Lc;
		ccx=Cx(i-Lc:i+Lc,j-Lc:j+Lc);
		ccxi=ccx(2*Lc+1:-1:1,2*Lc+1:-1:1);
		ccx=ccx+ccxi;
		ccx(Lc+1,Lc+1)=ccx(Lc+1,Lc+1)/2;
		ccx=vector(ccx');
		Tcx(nm,:)=ccx(1:M)';
	end
end
i=2*Lc+1;
for j=Lc+1:2*Lc+1,
	nm=(i-Lc-1)*(2*Lc+1)+j-Lc;
	ccx=Cx(i-Lc:i+Lc,j-Lc:j+Lc);
	ccxi=ccx(2*Lc+1:-1:1,2*Lc+1:-1:1);
	ccx=ccx+ccxi;
	ccx(Lc+1,Lc+1)=ccx(Lc+1,Lc+1)/2;
	ccx=vector(ccx');
	Tcx(nm,:)=ccx(1:M)';
end
if sum(isnan(Tcx(:)))>0,
    disp(sprintf('NaNs in Tcx = %d',sum(isnan(Tcx(:)))));
end

% Rearrange Cy indices and solve the equation

Cy1=vector(Cy');
Cy1=Cy1(1:M);

Ch1=pinv(Tcx)*Cy1;

% Rearrange Ch1

Ch1=[Ch1;Ch1(length(Cy1)-1:-1:1)];
Ch=reshape(Ch1,Nc,Nc)';

% Compute Y as conv(X,H) in the Fourier domain
aux=zeros(Ny,Nx);
aux(cy-Lc:cy+Lc,cx-Lc:cx+Lc)=Ch;
Ch=fftshift(aux);
Chf=real(fft2(Ch));
Xf=fft2(X);
Xf2=abs(Xf).^2;

Yf=Xf.*sqrt(abs(Chf));
Y=ifft2(Yf);


return