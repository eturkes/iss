function o = SplitAnchorFindSpots(o)
% o=iss_SplitAnchorFindSpots(o)
%
% Detects spots in split anchor channels.
% Register the split anchor channels to the full anchor channel.
% Transform the split anchor spot coordiantes to the full anchor channel.
% Replace the full anchor spots with these transformed split anchor spots.
%
% 
% Kenneth D. Harris, 29/3/17
% GPL 3.0 https://www.gnu.org/licenses/gpl-3.0.en.html


%% Detect the spots
nSplitChannels = size(o.SplitAnchorChannels,2);

rr = o.ReferenceRound;
NonemptyTiles = find(~o.EmptyTiles)';

[nY, nX] = size(o.EmptyTiles);
nTiles = nY*nX;

RawLocalFullAnchorYXZ = o.RawLocalYXZ;  % cell array, giving spots in local coordinates
RawFullAnchorIsolated = o.RawIsolated;

RawLocalSplitAnchorYXZ = cell(nTiles,3);
RawSplitAnchorIsolated = cell(nTiles,3);

fprintf('\nDetecting Split Anchor spots in tile   ');
for t=NonemptyTiles
    if t<10
        fprintf('\b%d',t);
    else
        fprintf('\b\b%d',t);
    end 
    
    [y,x] = ind2sub([nY nX], t);
    for b = o.SplitAnchorChannels
        AnchorIm = o.load_3D(rr,y,x,b)-o.TilePixelValueShift;
        [RawLocalSplitAnchorYXZ{t,find(o.SplitAnchorChannels==b)},...
            RawSplitAnchorIsolated{t,find(o.SplitAnchorChannels==b)}] = o.detect_spots(AnchorIm,t,b,rr);
    end
end
fprintf('\n');

save(fullfile(o.OutputDirectory, 'SplitAnchorSpots.mat'), 'o', 'RawLocalFullAnchorYXZ','RawFullAnchorIsolated',...
    'RawLocalSplitAnchorYXZ','RawSplitAnchorIsolated');

%% Find initial shift between FullAnchor and SplitAnchor

D0 = zeros(nTiles,3,nSplitChannels);
Scores = zeros(nTiles,nSplitChannels);
ChangedSearch = zeros(nSplitChannels,1);
OutlierShifts = zeros(nTiles,3,nSplitChannels);

SearchRange.Y = -50:o.FindSpotsStep(1):50;
SearchRange.X = -50:o.FindSpotsStep(2):50;
SearchRange.Z = 0;

for t=1:nTiles
    if o.EmptyTiles(t); continue; end
    for b = 1:nSplitChannels
        [D0(t,:,b), Scores(t,b),tChangedSearch] = o.get_initial_shift2(RawLocalFullAnchorYXZ{t},...
            RawLocalSplitAnchorYXZ{t,b}, SearchRange,'Register');
        ChangedSearch(b) = ChangedSearch(b)+tChangedSearch;
        
        fprintf('Tile %d, shift from full anchor round to Split Anchor %d: [%d %d %d], score %f\n', t, b, D0(t,:,b),...
            Scores(t,b));        
    end
end

save(fullfile(o.OutputDirectory, 'SplitAnchorSpots.mat'), 'o', 'RawLocalFullAnchorYXZ','RawFullAnchorIsolated',...
    'RawLocalSplitAnchorYXZ','RawSplitAnchorIsolated','D0','Scores');

%% Do PCR

[A,nMatches,Error] = PointCloudRegisterSplitAnchor(o,RawLocalFullAnchorYXZ, RawLocalSplitAnchorYXZ, nTiles, nSplitChannels,D0);

save(fullfile(o.OutputDirectory, 'SplitAnchorSpots.mat'), 'o', 'RawLocalFullAnchorYXZ','RawFullAnchorIsolated',...
    'RawLocalSplitAnchorYXZ','RawSplitAnchorIsolated','D0','Scores','A','nMatches','Error');

%% Transform Split Anchor coordinates
%Scale Spots in Z direction then transform and descale
ScaledRawLocalSplitAnchorYXZ = cell(nTiles,3);
ScaledTransformedLocalSplitAnchorYXZ = cell(nTiles,3);
TransformedLocalSplitAnchorYXZ = cell(nTiles,3);
for t=1:nTiles
    if o.EmptyTiles(t); continue; end
    for b=1:nSplitChannels
        ScaledRawLocalSplitAnchorYXZ(t,b) = {[RawLocalSplitAnchorYXZ{t,b}.*[1,1,o.Zpixelsize/o.XYpixelsize],...
            ones(size(RawLocalSplitAnchorYXZ{t,b},1),1)]};
        ScaledTransformedLocalSplitAnchorYXZ(t,b) = {(ScaledRawLocalSplitAnchorYXZ{t,b}*A(:,:,t,b))};  %Transform
        TransformedLocalSplitAnchorYXZ(t,b) = {round(ScaledTransformedLocalSplitAnchorYXZ{t,b}.*[1,1,o.XYpixelsize/o.Zpixelsize])};
    end
end

save(fullfile(o.OutputDirectory, 'SplitAnchorSpots.mat'), 'o', 'RawLocalFullAnchorYXZ','RawFullAnchorIsolated',...
    'RawLocalSplitAnchorYXZ','RawSplitAnchorIsolated','D0','Scores','A','nMatches','Error','ScaledRawLocalSplitAnchorYXZ',...
    'ScaledTransformedLocalSplitAnchorYXZ','TransformedLocalSplitAnchorYXZ');

%% Set o.RawLocalYXZ and o.RawIsolated as the values from the SplitAnchor
o.RawLocalYXZ = cell(nTiles,1);
o.RawIsolated = cell(nTiles,1);
o.RawChannel = cell(nTiles,1);

for t=1:nTiles
    if o.EmptyTiles(t); continue; end
    o.RawLocalYXZ(t,1) = {[TransformedLocalSplitAnchorYXZ{t,1};TransformedLocalSplitAnchorYXZ{t,2};TransformedLocalSplitAnchorYXZ{t,3}]};
    o.RawIsolated(t,1) = {[RawSplitAnchorIsolated{t,1};RawSplitAnchorIsolated{t,2};RawSplitAnchorIsolated{t,3}]};
    RawChannel = ones(size(TransformedLocalSplitAnchorYXZ{t,1},1),1)*o.SplitAnchorChannels(1);
    for c = 2:nSplitChannels
        RawChannel = [RawChannel;ones(size(TransformedLocalSplitAnchorYXZ{t,c},1),1)*o.SplitAnchorChannels(c)];
    end
    o.RawChannel(t,1) = {RawChannel};
end

end
