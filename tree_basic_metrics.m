function [h,g,c,d] = tree_basic_metrics(tree,varargin)
% Tree overall information extraction system for an object TREE
% USE:
%	[H,G,C] = TREE_METRICS(TREE,VARARGIN)
%
% TREE is a tree object (see elsewhere)
% H: HEIGHT of the tree, that is the highest point over Z-axis.
% G: GIRTH of the tree, that is diameter at certain height of the trunk
% (1.3 m for trees with trunk, at the base, 0m, for non-trunk trees etc.)
% C: CROWN spread of the tree, that is projection area of the crown.
% Methods vary.
%

% Error check
if(~isa(tree,'tree'))
    error('Error: Object class must be tree.');
end

% Init
girth_height = 1.3;
dang = 1.0;

% Process input arguments
tf = strcmp('base',varargin);
if(find(tf))
    girth_height = 0;
end
tf = strcmpi('ang',varargin);
if(find(tf))
    dang = varargin{find(tf)+1};
end

h = height(tree);
g = girth(tree,girth_height);
c = crown(tree,dang);

end

function c = crown(tree,dang)
dn = dang;
deg = -180:dn:180;
pnts = zeros(length(deg)-1,2);
spokes = zeros(length(deg)-1);
m = length(deg);
angles = (180/pi)*atan2(tree.end_point(:,2),tree.end_point(:,1));
center = [tree.start_point(1,1) tree.start_point(1,2)];
dists = sqrt((tree.end_point(:,1)-center(1)).^2 + (tree.end_point(:,2)-center(2)).^2);
for jj = 2:m
    ii = find(angles < deg(jj) & angles > deg(jj-1));
    if(isempty(ii)), continue; end
    %pnts(jj-1,:) = [tree.end_point(ii(idx),1) tree.end_point(ii(idx),2)];
    spokes(jj-1) = max(dists(ii));
end
%pnts = pnts(pnts(:,1) ~= 0,:);% remove not found angular slices
spokes = spokes(spokes ~= 0);
% plot(pnts(:,1),pnts(:,2),'ob');
% hold on;
% lineX = [];
% lineY = [];
% for ii = 1:length(pnts(:,1))
%     lineX = cat(2,lineX,[0 pnts(ii,1)]);
%     lineY = cat(2,lineY,[0 pnts(ii,2)]);
% end
% line(lineX,lineY,'Color','r');
% plot(0,0,'og','LineWidth',2);
% hold off;
c = 2*mean(spokes);
end

function g = girth(tree,girth_height)
ii = 1;
g = 0.0;
if( (abs(tree.start_point(ii,3)-tree.start_point(1,3)) <= girth_height) && ...
        (abs(tree.end_point(ii,3)-tree.start_point(1,3)) >= girth_height) )
    g = 2 * tree.radius(ii);
    return;
end
while( tree.extension(ii) )
    ii = tree.extension(ii);
    if( (abs(tree.start_point(ii,3)-tree.start_point(1,3)) <= girth_height) && ...
            (abs(tree.end_point(ii,3)-tree.start_point(1,3)) >= girth_height) )
        g = 2 * tree.radius(ii);
        return;
    end
end
end

function h = height(tree)
h = max(max(tree.start_point(:,3)),max(tree.end_point(:,3)));
h = h - tree.start_point(1,3);% assuming the 1st cyl is a ground cyl
end
