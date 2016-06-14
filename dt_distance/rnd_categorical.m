function r = rnd_categorical(n, p)
% Return n iid samples from categorical distribution that has its 
% weights in p, dimension equals length(p)
% - MJ

% input:
% n, number of samples to be generated
% p, weights (probabilities) of the categories
%
% output:
% r, requested samples

p = abs(p);
p = p/sum(p);

% p = cumsum(p);
% [~, i] = histc(rand(1,n),[0 p]);
% v = 1:length(p);
% r = v(i);

r = datasample(1:length(p),n,'weights',p);

end




