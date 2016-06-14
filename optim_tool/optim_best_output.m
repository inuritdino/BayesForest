function [state, options,optchanged] = optim_best_output(options,state,flag)
% A simple function to remember best solution in each generation.

optchanged = false;

fname = 'gaOut.dat';
switch flag
    case 'init'
        fid = fopen('gaOut.dat','w');
        fprintf(fid,'%s\t%s\t%s\t%s\n','Generation','Best_score','Mean_score','Population');
        fclose(fid);
    case {'iter','interrupt'}
        fid = fopen(fname,'a');
        ii = state.Generation;
        best = state.Best(ii);
        avg = mean(state.Score);
        [~,I] = min(abs(state.Score - best));
        fprintf(fid,'%d\t%g\t%g\t%s\n',ii,best,avg,num2str(state.Population(I,:)));
        fclose(fid);
    case 'done'
end

end