function C = weighted_covariance(A,w)
% Treats the columns of A as vectors i.e. each row is s sample vector, 
% returning a sample cov matrix. The weights (not necessarily normalised) 
% are give in vector w.
% The weights are not the number of repetitions but tell the importance of
% the corresponding sample
% - MJ

% input:
% A, samples
% w, weights for the samples
%
% output:
% C, weighted sample covariance matrix

if size(A,1) == 1 || nargin < 2 || isempty(w)
    C = cov(A);
    return;
end

% length(w) must be the same as size(A,1)
if length(w) ~= size(A,1)
    disp('Warning! Size of weights do not match with the samples!');
    disp('Computed unweighted covariance matrix instead!');
    C = cov(A);
    return;
end
    
% normalise weights just in case
w = w./sum(w);
w = w(:);

% use the standard formula to compute the weighted covariance
M = weighted_mean(A,w);
nconst = sum(w)/(sum(w)^2 - sum(w.^2));
C = 0;
% todo: remove this loop, it can be slow
for i = 1:size(A,1)
    C = C + w(i)*(A(i,:) - M)' * (A(i,:) - M);
end
C = nconst * C;

end




