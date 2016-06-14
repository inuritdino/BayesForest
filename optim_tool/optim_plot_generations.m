function mov = optim_plot_generations(gaInput,model_fun,data_tr)
% Plot the solution trees by generations
%
%

%% Initial
nline = 0;
pop = [];% population matrix
best = [];% best score array
gens = [];% generations
averg = [];% average score array

%% Get the input
fid = fopen(gaInput,'r');
line = fgetl(fid);
while(ischar(line))
    if(~isempty(line))
        if(nline == 0)
            C = textscan(line,'%s');
            ncols = length(C{1});
            if(ncols > 4)
                error('ERROR: #columns exceeds 4.');
            end
            %colnames = C{1};
            nline = nline + 1;
        else
            C = textscan(line,'%f','delimiter','\t');
            gens = cat(2,gens,C{1}(1));
            best = cat(2,best,C{1}(2));
            if(ncols == 4)
                averg = cat(2,averg,C{1}(3));
            end
            pop = cat(1,pop,C{1}(ncols:end)');
            nline = nline + 1;
        end
    end
    line = fgetl(fid);
end

%% Plotting
% Positioning in the plot
DataSideViewPos = {1, 3, 1};% data tree side view
ModelSideViewPos = {1, 3, 2};% model tree side view
DataTopViewPos = {3, 3, 3};% data top view
ModelTopViewPos = {3, 3, 6};% model top view
ScorePlotPos = {3, 3, 9};% score plot

h = figure('position',[100 100 1000 500]);
mov(1:length(gens)) = struct('cdata',[],'colormap',[]);

% Score plot labels
%subplot(231);grid on;
subplot(ScorePlotPos{:});grid on;
xlim([1 max(gens)]);
mn = min(best);
mx = max(best);
% if(~isempty(averg))
%     mx = max(averg);
% else
%     mx = max(best);
% end
dy = mx - mn;
% if(ncols ~= 4)
%     ylabel('Best score');
% else
%     ylabel('Best(blue)/Average(green) score');
% end
ylabel('Best score');
xlabel('Generations');
ylim([mn-0.5*dy mx+0.5*dy]);

% Plot the data tree, both views
data_tr = data_tr.move_tree([0 0 0]);
%subplot(232);
subplot(DataSideViewPos{:});
data_tr.draw; view(0,0);
%subplot(235);
subplot(DataTopViewPos{:});
data_tr.draw;
view(0,90);
% Determine the data spans for correct scaling
[xmin0,xmax0,ymin0,ymax0,zmax0] = span(data_tr);

for ii = 1:length(best);
    subplot(ScorePlotPos{:});hold on;% Best score plot
    plot(gens(1:ii-1),best(1:ii-1),'b.','MarkerSize',10);
    plot(gens(ii),best(ii),'r.','MarkerSize',10);
    title(['Best = ' num2str(best(ii))]);
%     title_str = ['Best = ' num2str(best(ii))];
%     if(ncols == 4)% average score is present
%         plot(gens(1:ii-1),averg(1:ii-1),'g.','MarkerSize',10);
%         plot(gens(ii),averg(ii),'k.','MarkerSize',10);
%         title_str_app = ['; avg = ' num2str(averg(ii))];
%     end
%     title([title_str title_str_app]);
    % GET THE MODEL TREE
    if(ii > 1 && all(abs(pop(ii,:)-pop(ii-1,:)) < 1e-04) && all(abs(best(ii)-best(ii-1)) < 1e-04))
        fprintf('Ignore generation %d.\n',gens(ii));
        mov(ii) = getframe(h);
        continue;
    end
    model_fun(pop(ii,:));
    tr = read_mtg('out.mtg');
    [xmin1,xmax1,ymin1,ymax1,zmax1] = span(tr);% spans
    %%%%%%%%%%%%%%%
    % Plot at view (0,0)
    subplot(ModelSideViewPos{:});
    tr.draw;
    view(0,0);
    xlim([min(xmin0,xmin1) max(xmax0,xmax1)]);
    zlim([0 max(zmax0,zmax1)]);
    subplot(DataSideViewPos{:});% change the corresponding view of the data tree
    xlim([min(xmin0,xmin1) max(xmax0,xmax1)]);
    zlim([0 max(zmax0,zmax1)]);
    % Plot at view (0,90)
    subplot(ModelTopViewPos{:});
    tr.draw;
    view(0,90);
    xlim([min(xmin0,xmin1) max(xmax0,xmax1)]);
    ylim([min(ymin0,ymin1) max(ymax0,ymax1)]);
    subplot(DataTopViewPos{:});% change the corresponding view of the data tree
    xlim([min(xmin0,xmin1) max(xmax0,xmax1)]);
    ylim([min(ymin0,ymin1) max(ymax0,ymax1)]);
    % Get the frame
    drawnow;% update the frame immediately
    mov(ii) = getframe(h);
end


%% Close the file
fclose(fid);

end

function [xmin,xmax,ymin,ymax,zmax] = span(tr)
% Find the span in directions
ADDUP = 0.5;
if(min(tr.end_point(:,1)) < min(tr.start_point(:,1)))
    [~,xmin] = min(tr.end_point(:,1));
    xmin = tr.end_point(xmin,1) - ADDUP;
else
    [~,xmin] = min(tr.start_point(:,1));
    xmin = tr.start_point(xmin,1) - ADDUP;
end
if(max(tr.end_point(:,1)) > max(tr.start_point(:,1)))
    [~,xmax] = max(tr.end_point(:,1));
    xmax = tr.end_point(xmax,1) + ADDUP;
else
    [~,xmax] = max(tr.start_point(:,1));
    xmax = tr.start_point(xmax,1) + ADDUP;
end

if(min(tr.end_point(:,2)) < min(tr.start_point(:,2)))
    [~,ymin] = min(tr.end_point(:,2));
    ymin = tr.end_point(ymin,2) - ADDUP;
else
    [~,ymin] = min(tr.start_point(:,2));
    ymin = tr.start_point(ymin,2) - ADDUP;
end
if(max(tr.end_point(:,2)) > max(tr.start_point(:,2)))
    [~,ymax] = max(tr.end_point(:,2));
    ymax = tr.end_point(ymax,2) + ADDUP;
else
    [~,ymax] = max(tr.start_point(:,2));
    ymax = tr.start_point(ymax,2) + ADDUP;
end

if(max(tr.end_point(:,3)) > max(tr.start_point(:,3)))
    [~,zmax] = max(tr.end_point(:,3));
    zmax = tr.end_point(zmax,3) + ADDUP;
else
    [~,zmax] = max(tr.start_point(:,3));
    zmax = tr.start_point(zmax,3) + ADDUP;
end
end