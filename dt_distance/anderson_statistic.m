function [W2, A2] = anderson_statistic(meas, real_cdf, two_sample)
% Computes Anderson-Darling statistic in 1d between one dataset and the
% known model cdf or between two datasets
% - MJ

% input:
% meas, the first set of datapoints (in increasing order!)
% real_cdf, if two_sample == 0 then these are the values of the model cdf 
% computed at ordered datapoints meas, if two_sample == 1 these are the 
% second set of datapoints (also in increasing order!) 
% two_sample, 1 if two sample version of the test, 0 if one sample version 
% i.e. the model cdf is known and is given as input
% 
% output:
% W2, A2, both are the usual values of the statistic, A2 is scaled ~ sqrt(n)
%
% Note:
% 1) if two_sample == 0 then meas is not needed and real_cdf is only needed
% to compute the statistic
% 2) if two_sample == 1 then real_cdf are the second set of datapoints and not the 
% cdf values. Datasets may not be of the same size.
% 3) some inaccuracies may occur if samples are not unique in two sample
% case
% 4) if the cdf attains value that equals or is close to 0 or 1, numerical
% issues can emerge! This happens because of the weighting and since the
% probability of generating out-of-support data value is zero if the model
% is correctly specified and this may cause difficulties with optimization.

if nargin < 3
    two_sample = 0;
end

if ~two_sample
    n = length(real_cdf);
    S = sum((2*(1:n)-1)/n .* (log(real_cdf(1:n)) + ...
        log(1 - real_cdf(n+1-(1:n))))); 
    % check the value of S, see note 4
    if(~isfinite(S))
        %disp('Inf value in Anderson-Darling statistic!');
        S = -1e10;
    end
    A2 = -n - S;
    W2 = A2/n; 
    
else % two sample case
    n1 = length(meas);
    n2 = length(real_cdf);
    N = n1 + n2;
    
    % pooled sample 
    all = [meas, real_cdf];
    [~, iall] = sort(all);
    
    % compute mi's i.e. the #samples in the first set <= i:th smallest in 
    % the pooled samples 
    mi = zeros(1,N-1);
    mi(iall <= n1) = 1;
    mi = cumsum(mi);
    
    % compute the statistic
    S = sum((mi(1:N-1)*N - n1*(1:N-1)).^2./((1:N-1).*(N-(1:N-1))));
    A2 = 1/(n1*n2) * S;
    W2 = N/(n1*n2)^2 * S; 
end

end


%     S = 0;
%     for i = 1:n
%        S = S + (2*i-1)/n * (log(real_cdf(i)) + log(1 - real_cdf(n+1-i))); 
%     end

%     for i = 1:N-1
%         if iall(i) <= n1 % current point belongs to the first set
%             mi(i) = mi(i) + 1;
%         end
%     end

%     S = 0;
%     for i = 1:N-1
%        S = S + (mi(i)*N - n1*i)^2/(i*(N-i)); 
%     end





