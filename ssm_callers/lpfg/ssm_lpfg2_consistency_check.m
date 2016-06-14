function Ok = ssm_lpfg2_consistency_check( N, X, ssm_lpfg2_h )
%SSM_LPFG2_CONSISTENCY_CHECK Check the consistency of the output from the
%ssm_lpfg2() function. Namely, it checks whether the output is the same for
%N simulations given the same parameter set. This is crucial for the
%optimization procedure in the Bayes Forest.
% USAGE:
%     OK = SSM_LPFG2_CONSISTENCY_CHECK( N, X, SSM_LPFG2_H )
% N: number of simulations to run
% X: arguments to the
% SSM_LPFG2_H: function handle to the ssm_lpfg2() function specifying
% arguments.
% OK: the output is consistent (true) or not (false)

Ok = true;% assume consistency
in = ssm_lpfg2_h(X);
for ii = 2:N
    out = ssm_lpfg2_h(X);
    n = numel(out);
    if(n == numel(in))
        for jj = 1:n
            if(any(size(out{jj}) ~= size(in{jj})))
                Ok = false;
                break;
            end
        end
    end
    if(Ok == false)
        break;
    end
end
end

