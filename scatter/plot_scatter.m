function plot_scatter(data,model)
% Plot the DATA and MODEL scatters over each other for visual comparison.
% USAGE:
%       PLOT_SCATTER(DATA,MODEL)
%
% 
% Function plots only several scatters on top of each other.
% NOTE: DATA/MODEL must be matrix of size NxM, where N is the number of
% variates, M is the number of samples (basically wide matrices ncol >
% nrow). Several scatters can be combined into a cell array then every
% DATA/MODEL{I} contains such a matrix.
%
% For every scatter/cell array I the function plots a row of
% figures/subplot(s). In each row the corresponding variates (2D max) for
% the given scatter I are plotted against each other for DATA and MODEL
% such that number of columns of the figure is equal N-1 (the following
% pairs of variates are plotted: 1-2, 2-3, 3-4,... (N-1)-N, there are N-1 such
% pairs in total).
%
% Restrictions:
% I <= 5 (not more than 5 rows of the figure)
% N-1 <= 5 (not more than 5 columns of the figure)
%


if(isa(data,'cell'))
    N = length(data);
else
    N = 1;
end
% Dimension check
if(isa(model,'cell'))
    if(length(model) ~= N)
        error('Error: different dimensions');
    end
else
    if(N ~= 1)
       error('Error: different dimensions');
    end
end
% if(N == 1)
%     model = {model};
%     data = {data};
% end
for ii = 1:N
    if(size(model{ii},1) > size(model{ii},2))% invert if columns=vars
        model{ii} = model{ii}';
    end
    if(size(data{ii},1) > size(data{ii},2))% invert if columns=vars
        data{ii} = data{ii}';
    end
end
% Restrict number of datasets to 5
if(N > 5)
    N = 5;
end
% % Further find out multi-dimensional data
% if(isa(data,'cell'))
%     
% end
% New figure
h = figure;
set(h,'Position',[0 0 1000 250*N]);
%
for ii = 1:N
    n = size(data{ii},1);
    if(n > 1)
        M = n-1;
        if(M > 5), M = 5; end % restrict to 5 columns
        for jj = 1:M
            subplot(N,M,(ii-1)*M+jj)
            plot(data{ii}(jj,:),data{ii}(jj+1,:),'ob','MarkerSize',5,'LineWidth',1);
            hold on;
            plot(model{ii}(jj,:),model{ii}(jj+1,:),'or','MarkerSize',5,'LineWidth',1);
            hold off
            xlabel(['Dataset ' num2str(ii) ', var ' num2str(jj)]);
            ylabel(['Dataset ' num2str(ii) ', var ' num2str(jj+1)]);
            if(jj == 1)
                legend('DATA','MODEL');
                title(['$n_{\mathrm{data}}$ = ' num2str(size(data{ii},2)) ...
                    ', $n_{\mathrm{model}}$ = ' num2str(size(model{ii},2))],'interpreter','latex');
            end
        end
    elseif(n == 1)
        subplot(N,M,ii);
        low = min(min(data{ii}),min(model{ii}));
        high = max(max(data{ii}),max(model{ii}));
        step = (high-low)/10;
        x = low:step:high;
        n1 = histc(data{ii},x);
        n2 = histc(model{ii},x);
        b = bar(x,[n1;n2]');
        b(1).FaceColor = 'b'; b(1).EdgeColor = 'b';
        b(2).FaceColor = 'r'; b(2).EdgeColor = 'r';
        xlabel(['Dataset ' num2str(ii)]);
    end
end


end
