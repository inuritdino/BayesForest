function M = weighted_mean(A,w)
% treats the columns of A as vectors i.e. each row is s sample vector,
% returning a row vector of mean valueswith the weights (not necessarily 
% normalised) in vector w
% The weights are not the number of repetitions but tell the importance of
% the corresponding sample
% - MJ

% input:
% A, samples
% w, weights for the samples
%
% output:
% M, weighted mean

if size(A,1) == 1 || nargin < 2 || isempty(w)
    M = mean(A);
    return;
end

% length(w) must be the same as size(A,1)
if length(w) ~= size(A,1)
    disp('Warning! Size of weights do not match with the samples!');
    disp('Computed unweighted mean instead!');
    M = mean(A);
    return;
end

w = w./sum(w);
w = w(:);
M = sum(A.*repmat(w,1,size(A,2)),1);

end




