function [out,histFlag,snr_skew,snr_kurt,snr_skew_within,snr_kurt_within] = imposeHistogram(ht, seed, verbose, alwaysHistogram)
% [out, histFlag, snr, snr_skew, snr_kurt, snr_skew_within, snr_kurt_within] = 
%       imposeHistogram(ht, seed, verbose, alwaysHistogram)
%
% Standard Portilla & Simoncelli runs into trouble on big blank regions
% because parameters like skew and kurtosis can be undefined. Instead of
% matching skew and kurtosis in those conditions, this instead does a
% histogram match.
%  
%   ht: desired statistics
%   seed: starting point -- matrix that will be adjusted to match stats
%   verbose: display status updates or no
%   alwaysHistogram: if = 1, don't even try to match skew and kurtosis;
%       default = 1
%   out: output once adjusted statistics
%   histFlag: 1 if matched histogram instead of skew and kurtosis
%   snr, snr_skew, etc: legacy from P&S code

histFlag = 0;

if (~exist('alwaysHistogram','var') || isempty(alwaysHistogram))
    alwaysHistogram = 1;
end

% disp(sprintf('NaNs in histogram N to be imposed = %d',sum(isnan(ht.n(:)))));
% disp(sprintf('NaNs in histogram X to be imposed = %d',sum(isnan(ht.x(:)))));
% disp(sprintf('NaNs in image whose histogram is to be modified = %d',sum(isnan(seed(:)))));

X = seed;

if (alwaysHistogram)
    X = histoMatch(X, ht.n, ht.x);    
    %was X = imposeMomentStatistics(ht.moments, X, ones(size(X)), 1);
    out = X;
else

    Xcopy = X;
    try
        [X,snr_skew,snr_skew_within,problemFlagSkew] = modskew(X, ht.moments(3), 1);
        [X,snr_kurt,snr_kurt_within,problemFlagKurt] = modkurt(X, ht.moments(4), 1); 
        if problemFlagSkew + problemFlagKurt > 0
            disp('NaNs detected: Matching full histogram instead of only skewness and kurtosis');
            X = Xcopy;

            snr_skew_within(1) = snr(ht.moments(3),ht.moments(3)-skew2_forTTM(Xcopy));
            snr_kurt_within(1) = snr(ht.moments(4),ht.moments(4)-kurt2_forTTM(Xcopy));

            if (~isempty(ht))
                X = histoMatch(X, ht.n, ht.x);
            end
            histFlag = 1;

            snr_skew = snr(ht.moments(3),ht.moments(3)-skew2_forTTM(X)); snr_skew_within(2) = snr_skew;
            snr_kurt = snr(ht.moments(4),ht.moments(4)-kurt2_forTTM(X)); snr_kurt_within(2) = snr_kurt;    
        end
    catch
        % suppress this print statement
        % there was a verbose flag for it previously
        % disp('Catch block: Matching full histogram instead of only skewness and kurtosis');
        X = Xcopy;
    
        snr_skew_within(1) = snr(ht.moments(3),ht.moments(3)-skew2_forTTM(Xcopy));
        snr_kurt_within(1) = snr(ht.moments(4),ht.moments(4)-kurt2_forTTM(Xcopy));
    
        if (~isempty(ht) && var(X(:))>0)
            X = histoMatch(X, ht.n, ht.x);
        end
        histFlag = 1;
    
        snr_skew = snr(ht.moments(3),ht.moments(3)-skew2_forTTM(X)); snr_skew_within(2) = snr_skew;
        snr_kurt = snr(ht.moments(4),ht.moments(4)-kurt2_forTTM(X)); snr_kurt_within(2) = snr_kurt; 
    end

    out = X;
end
