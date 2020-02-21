function [A,nMatches,Error] = PointCloudRegisterSplitAnchor(o,y0, x0, nTiles, nSplitChannels,D0)     %MADE A THE SAME FOR ALL TILES
% o = o.PointCloudRegister(y, x, Options)
% 
% Perform point cloud registration to map points x onto points y by
% iterative closest point: repeatedly finding the best y for each x, 
% and doing linear regression to find the A that maps best maps x to y
%
% inputs:
% y0 is a cell containing the non centered YXZ location of spots in the 
% full anchor channel for all tiles
%
% x0 is a cell containing the non centered YXZ location of spots in the 
% split anchor channels for all tiles
%
% nSplitChannels is the number of channels to use in the anchor round
%
% D0(t,:,b) is the initial shift between the full anchor round and split
% anchor channel b for tile t
%%
nD = 3;

%centre anchor channel spots
y = cell(nTiles,1);
x = cell(nTiles,nSplitChannels);

for t=1:nTiles
    if o.EmptyTiles(t); continue; end
    %Center and scale z direction. Also append array of ones for translation
    y(t) = {y0{t}.*[1,1,o.Zpixelsize/o.XYpixelsize]};
    for b=1:nSplitChannels
        x(t,b) = {[x0{t,b}.*[1,1,o.Zpixelsize/o.XYpixelsize],ones(size(x0{t,b},1),1)]};
    end
end

if isempty(o.PcDist)
    o.PcDist = inf;
end

%Initialize variables
A = zeros(4,3,nTiles,nSplitChannels);
for t=1:nTiles    
    for b=1:nSplitChannels
        A(1:3,:,t,b) = eye(3);
        A(4,:,t,b) = D0(t,:,b);
    end
end


fprintf('\nPCR - Finding well isolated points');
% find well isolated points as those whose second neighbor is far
for t=1:nTiles
    if o.EmptyTiles(t); continue; end
    % make kd tree - default options!
    k0 = KDTreeSearcher(y{t});
    [~, d2] = k0.knnsearch(y{t}, 'k', 2);
    if isfinite(o.PcDist) && size(y{t},1) > 1
        y(t) = {y{t}(d2(:,2)>o.PcDist*2,:)};
    end
    
end

fprintf('\nPCR - Making kd trees');
%Make kd trees out of these well isolated points
k = cell(nTiles,1);
for t=1:nTiles
    if o.EmptyTiles(t); continue; end
    k(t) = {KDTreeSearcher(y{t})};
end


%%
UseMe = cell(nTiles,nSplitChannels);           %nP DIFFERENT FOR DIFFERENT TILES!!!
Neighbor = cell(nTiles,nSplitChannels);
MyNeighb = cell(nTiles,nSplitChannels);
xTransformedAnchor = cell(nTiles,nSplitChannels);
nMatches = zeros(nTiles,nSplitChannels);
Error = zeros(nTiles,nSplitChannels);

if isempty(o.ToPlot)
    fprintf('\nPCR - Iteration   ');
end
for i=1:o.PcIter
    if isempty(o.ToPlot)
        if i<10
            fprintf('\b%d', i);
        else
            fprintf('\b\b%d', i);
        end
        if i ==o.PcIter
            fprintf('\nPCR - Max number of iterations reached');
        end
    end
    
    LastNeighbor = Neighbor;
    
    for t=1:nTiles
        if o.EmptyTiles(t); continue; end
        for b=1:nSplitChannels
            xTransformedAnchor(t,b) = {(x{t,b}*A(:,:,t,b))'};
        end
    end
    
    %This part finds new neighbours and new estimates for A
    for t=1:nTiles
        xR = [];
        yR = [];
        if o.EmptyTiles(t); continue; end
        for b=1:nSplitChannels
            Neighbor(t,b) = {k{t}.knnsearch(xTransformedAnchor{t,b}')};
            [~,Dist] = k{t}.knnsearch(xTransformedAnchor{t,b}');
            UseMe(t,b) = {Dist<o.PcDist};
            MyNeighb(t,b) = {Neighbor{t,b}(UseMe{t,b}>0)};
            nMatches(t,b) = sum(UseMe{t,b});
            Error(t,b) = sqrt(mean(Dist(UseMe{t,b}>0).^2));
            
            A(:,:,t,b) = x{t,b}(UseMe{t,b}>0,:)\y{t}(MyNeighb{t,b},:);
        end
    end
    
    
    if min(min(min(cellfun(@isequal, Neighbor, LastNeighbor)))) == 1; break; end
    
end
fprintf('\n')

