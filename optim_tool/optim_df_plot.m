function optim_df_plot(data,model,hin)
% Simple plotting utility for the DATA and MODEL scatters fitting.
% USAGE:
%       OPTIM_DF_PLOT(DATA,MODEL)
%
% 

%% Preliminaries
if(nargin < 3)
    hin = [];
end
%% Find the number of scatters=DF's
if(isa(data,'cell'))
    N = length(data);
else
    N = 1;
end
if(isa(model,'cell'))
    if(length(model) ~= N)
        error('Error: dimensions of DATA and MODEL do not match.');
    end
else
    if(N ~= 1)
        error('Error: dimensions of DATA and MODEL do not match.');
    end
end

%% Arrange plots
if(isempty(hin))
    figure;
else
    figure(hin);
end
for ii = 1:N
    subplot(1,N,ii);
    if(N == 1 && size(data,1) == 1)% Histogram
    elseif(N ~= 1 && size(data{ii},1) == 1)% Histogram
    elseif(N == 1 && size(data,1) > 1)% Scatter
    elseif(N ~= 1 && size(data,1) > 1)% Scatter
    end
end
%     if(size(data,1) == 1)
%         subplot(121);
%         hist(data);
%         subplot(122);
%         hist(model_data);
%         title(['Dist = ' num2str(z)]);
%         label = [];
%         for ii=1:length(X)
%             label = [label 'X(' num2str(ii) ')=' num2str(X(ii)) 10];
%         end
%         xl = xlim;
%         yl = ylim;
%         text(0.6*xl(2),0.7*yl(2),label);
%     else% for higher dimensions
%         scatter(data(1,:),data(2,:));
%         hold on;
%         scatter(model_data(1,:),model_data(2,:),'r');
%         hold off;
%         title(['Dist = ' num2str(z)]);
%         label = [];
%         for ii=1:length(X)
%             label = [label 'X(' num2str(ii) ')=' num2str(X(ii)) 10];
%         end
%         xl = xlim;
%         yl = ylim;
%         text(0.6*xl(2),0.7*yl(2),label);
%     end


end