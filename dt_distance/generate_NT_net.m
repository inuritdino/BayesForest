function NT = generate_NT_net(k, s, method)
% Generates gp set of k numbers on unit-cube C_s
% (pp. 26-27 in Number-theoretic methods in statistics book)
% - MJ

% input:
% s: dimension
% k: # samples
% method: the algorithm for generating the quasi-random numbers in C_s:
% 1, gp-1 set
% 2, gp-2 set
% 3, cyclotomic field method
% 4, Halton sequence
% 5, Sobol sequence 
% 
% output:
% NT: the requested quasi-random numbers
%
% note:
% All the algorithms 1...5 should work for dimension <= 100, some
% algorithms will work in any dimension. Uniformity of the algorithms
% depend on dimension and some fixed settings, the settings used in this
% function may not be optimal for all cases.


% --- settings ---
default_method = 2; % this should be ok choice
% --- end of settings ---

if nargin < 3
    method = default_method; 
end

NT = NaN(k,s);
if method == 1
    % method (a) i.e. gp-1
    primes = [ ...
      2      3      5      7     11     13     17     19     23     29 ... 
     31     37     41     43     47     53     59     61     67     71 ... 
     73     79     83     89     97    101    103    107    109    113 ... 
    127    131    137    139    149    151    157    163    167    173 ... 
    179    181    191    193    197    199    211    223    227    229 ... 
    233    239    241    251    257    263    269    271    277    281 ... 
    283    293    307    311    313    317    331    337    347    349 ... 
    353    359    367    373    379    383    389    397    401    409 ... 
    419    421    431    433    439    443    449    457    461    463 ... 
    467    479    487    491    499    503    509    521    523    541];
    if s > length(primes)
        disp('Error in generating NT-net, one needs to set more primes!');
        return;
    end
    primes = primes(1:s);
    lambda = sqrt(primes);
    
elseif method == 2
    % method (b) i.e. gp-2
    p = 2; % choose some prime (what prime is the best? - unknown)
    q = p^(1/(s+1));
    lambda = q.^(1:s);
    
elseif method == 3
    % method (c) i.e. cyclotomic field method
    % choose some prime p
    % use smallest prime by default that is large enough
    primes = [ ...
      2      3      5      7     11     13     17     19     23     29 ... 
     31     37     41     43     47     53     59     61     67     71 ... 
     73     79     83     89     97    101    103    107    109    113 ... 
    127    131    137    139    149    151    157    163    167    173 ... 
    179    181    191    193    197    199    211    223    227    229 ... 
    233    239    241    251    257    263    269    271    277    281 ... 
    283    293    307    311    313    317    331    337    347    349 ... 
    353    359    367    373    379    383    389    397    401    409 ... 
    419    421    431    433    439    443    449    457    461    463 ... 
    467    479    487    491    499    503    509    521    523    541];

    ind = find(primes > 2*s + 3,1,'first'); 
    if isempty(ind)
    	disp('Error in generating NT-net, prime p must be set larger!');
        return;
    end
    p = primes(ind);
    lambda = rem(abs(2*cos(2*pi*(1:s)./p)),1);
    
end

if method == 1 || method == 2 || method == 3
    for i = 1:s
        NT(:,i) = rem(lambda(i)*(1:k),1);
    end
    
% from matlab:
elseif method == 4
    hs = haltonset(s); 
    NT = net(hs,k);
    
elseif method == 5
    ss = sobolset(s); 
    NT = net(ss,k);
    
else
    disp('Error! Incorrect algorithm for quasi-MC method.');
    
end
end




