function sa = rand_ab(n,dim,a,b)
% Generate uniform random samples from a cube that has its limits in 
% vectors a and b
% - MJ

% n, #samples to be generated
% dim, dimension of the rv's, must equal (or be smaller) the length(a) and length(b)
% a & b, limits for the uniform distribution, must hold that a < b (componentwise)

sa = NaN(dim,n);
nr = rand(dim,n);
for i = 1:dim
   sa(i,:) = a(i) + (b(i) - a(i))*nr(i,:); 
end

end




