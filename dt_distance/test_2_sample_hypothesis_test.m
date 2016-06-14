function [] = test_2_sample_hypothesis_test()
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% a quick code to test hypothesis test procedures
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% - MJ

%% test hypothesis testing for computing several times
if 1
    close all;
    
    n = 100;
    dim = 2;
    repeats = 10;
    settings.stat1d = 3;
    settings.norm_method = 1;
    settings.directions = 100;
    
    H = NaN(repeats,1);
    p_value = NaN(repeats,1);
    tic;
    for i = 1:repeats
        i
        [H(i), p_value(i), p_dist] = tree_significance_test(0.1+randn(dim,n), randn(dim,n), settings);
    end
    toc
    
    disp('results:');
    p_value'
    sum(H)
end


%% test how much time the hypothesis test will take
if 0
    close all;
    
    n = 10000;
    dim = 6;
    settings.stat1d = 1;
    settings.norm_method = 1;
    settings.directions = 100;
    
    tic;
    [H, p_value, p_dist] = tree_significance_test(randn(dim,n), randn(dim,n), settings)
    toc
    
    disp('results:');
    H
    
    % example times:
    % 10000 points, dim = 6, sector_lines = 4 -> MC draws 116, time: 123s
    % 20000 points, dim = 6, sector_lines = 4 -> MC draws 268, time: 515s
end


%% test the case where we have multiple datasets
if 0
    close all;
    
    n1 = 1000;
    n2 = 100;
    dim1 = 2;
    dim2 = 3;
    settings.stat1d = 1;
    settings.norm_method = 1;
    settings.directions = 100;
    
    % set data
    data1{1} = randn(dim1,n1);
    data2{1} = randn(dim1,n1);
    data1{2} = rand(dim2,n2);
    data2{2} = 4+rand(dim2,n2);
    
    tic;
    [H, p_value, p_dist] = tree_significance_test(data1, data2, settings)
    toc
    
    disp('results:');
    H
end

end







