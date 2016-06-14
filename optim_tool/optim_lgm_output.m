function [state, options,optchanged] = optim_lgm_output(options,state,flag)
% A simple function to remember best solution in a generation.

optchanged = false;

fname = 'gaOut.dat';
switch flag
    case 'init'
        fid = fopen('gaOut.dat','w');
        fprintf(fid,'%s\t%s\t%s\n','Generations','Best_score','Population');
        fclose(fid);
    case {'iter','interrupt'}
        fid = fopen(fname,'a');
        ii = state.Generation;
        best = state.Best(ii);
        [~,I] = min(abs(state.Score - best));
        fprintf(fid,'%d\t%g\t%s\n',ii,best,num2str(state.Population(I,:)));
        fclose(fid);
    case 'done'
end

end