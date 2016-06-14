function [W2, T] = cvm_statistic(meas, real_cdf, two_sample, watson_test) % ok
% Compute the Cramer von Mises (CvM) statistic (or its modifcation - Watson 
% statistic) in 1d between one dataset and the known model cdf or between 
% two datasets
% - MJ

% input:
% meas, the first set of datapoints (in increasing order!)
% real_cdf, if two_sample == 0 then these are the values of the model cdf 
% computed at ordered datapoints meas, if two_sample == 1 these are the 
% second set of datapoints (also in increasing order!) 
% two_sample, 1 if two sample version of the test, 0 if one sample version 
% i.e. the model cdf is known and is given as input
% watson_test, 1 if Watson statistic is computed instead of CvM
%
% output
% W2, T, both are the usual values of the statistic, T is scaled ~ sqrt(n)
%
% Note:
% 1) if two_sample == 0 then meas is not needed and real_cdf is only needed
% to compute the statistic
% 2) if two_sample == 1 then real_cdf are the second set of datapoints and not the 
% cdf values. Datasets may not be of the same size.

if nargin < 4
    watson_test = 0;
end
if nargin < 3
    two_sample = 0;
end

if ~two_sample
    n = length(real_cdf);
    T = 1/(12*n) + sum(((2*(1:n)-1)/(2*n) - real_cdf).^2);
    W2 = T/n;
    if watson_test % use the Watson test statistic instead of CvM
        Fbar = sum(real_cdf)/n;
        U2 = W2 - n*(Fbar - 1/2)^2;
        T = U2;
        W2 = U2/n;
    end
    
else % two sample case
    n1 = length(meas);
    n2 = length(real_cdf);
    N = n1 + n2;
    
    % pool the samples
    all = [meas, real_cdf];
    
%     % compute the ranks i.e. ordinal numbers of the original samples in the 
%     % pooled sample
%     % (assume samples are unique!)
%     [~, iall] = sort(all);
%     a = 1:N;
%     r = a(iall <= n1);
%     s = a(iall > n1);
%     %r,s
%     
%     % now compute U
%     U = n1*sum((r-(1:n1)).^2) + n2*sum((s-(1:n2)).^2);
%     D = U;
%     
%     % and finally compute T
%     T = U/(n1*n2*N) - (4*n1*n2 - 1)/(6*N);
    
    % simpler to compute the same thing (without uniqueness assumption or 
    % dealing with the ranks)
    all = sort(all);
    ecdf1 = emp_cdf(meas,all);
    ecdf2 = emp_cdf(real_cdf,all);
    W2 = 1/N * sum((ecdf1 - ecdf2).^2);
    T = n1*n2/N * W2;
end

end




