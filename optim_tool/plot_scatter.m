function plot_scatter(data,model)
% Plot the DATA and MODEL scatters over each other for visual comparison.
%
%
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
% Restrict number of plots to 5
if(N > 5)
    N = 5;
end
% % Further find out multi-dimensional data
% if(isa(data,'cell'))
%     
% end
% New figure
h = figure;
set(h,'Position',[0 0 1000 350]);
%
for ii = 1:N
    subplot(1,N,ii);
    if(size(data{ii},1) > 1)
        plot(data{ii}(1,:),data{ii}(2,:),'ob','MarkerSize',8,'LineWidth',2);
        hold on;
        plot(model{ii}(1,:),model{ii}(2,:),'or','MarkerSize',8,'LineWidth',2);
        hold off
        legend('DATA','MODEL');
        xlabel(['Dataset ' num2str(ii) ', var 1']);
        ylabel(['Dataset ' num2str(ii) ', var 2']);
    elseif(size(data{ii},1) == 1)
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