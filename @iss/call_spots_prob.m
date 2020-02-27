function o = call_spots_prob(o)
%% o = o.call_spots
% calls spots to codes for in-situ sequencing. Run this after find_spots
% 
% produces SpotGene{Spot}: name of gene for each spot
% SpotCode{Spot}: text representation of code for each spot 
% SpotScore(Spot): score saying how well the code fits (0...1)
% SpotIntensity(Spot): RMS intensity of the spot
% 
% Using o.UseChannels and o.UseRounds, you can do spot calling
% without using certain rounds and colour channels.
%
% Kenneth D. Harris, 29/3/17
% GPL 3.0 https://www.gnu.org/licenses/gpl-3.0.en.html

%% Make Bleed Matrices
%Only using channels and rounds given by o.UseChannels and o.UseRounds
if isempty(o.UseChannels)
    o.UseChannels = 1:o.nBP;
end
    
if isempty(o.UseRounds)
    o.UseRounds = 1:o.nRounds;
end

nChans = size(o.UseChannels,2);
nRounds = size(o.UseRounds,2);

%Need to normalise rounds and channels to get bleed matrix
p = prctile(o.cSpotColors, o.SpotNormPrctile);
SpotColors = bsxfun(@rdivide, o.cSpotColors, p);

% now we cluster the intensity vectors to estimate the Bleed Matrix
NormBleedMatrix = zeros(nChans,nChans,nRounds); % (Measured, Real, Round)
if strcmpi(o.BleedMatrixType,'Separate')
    for r=o.UseRounds
        m = squeeze(SpotColors(o.cSpotIsolated,o.UseChannels,r)); % data: nCodes by nBases
    
        [Cluster, v, s2] = ScaledKMeans(m, eye(nChans));
        for i=1:nChans
            NormBleedMatrix(:,i,find(o.UseRounds==r)) = v(i,:) * sqrt(s2(i));
        end
    end
    
elseif strcmpi(o.BleedMatrixType,'Single')
    m = permute(squeeze(squeeze(SpotColors(o.cSpotIsolated,o.UseChannels,o.UseRounds))),[1 3 2]);
    m = squeeze(reshape(m,[],size(m,1)*nRounds,nChans));
    [Cluster, v, s2] = ScaledKMeans(m, eye(nChans));
    for i=1:nChans
        NormBleedMatrix(:,i,1) = v(i,:) * sqrt(s2(i));
    end
    for r=2:nRounds
        NormBleedMatrix(:,:,r) = NormBleedMatrix(:,:,1);
    end
    
else
    warning('Wrong o.BleedMatrixType entry, should be either Separate or Single')
end

%Unnormalise Bleed matrices by multiplying rows by percentiles
BleedMatrix = zeros(nChans,nChans,nRounds);
for r=1:o.nRounds
    for b=1:o.nBP
        BleedMatrix(b,:,r) = p(:,b,r)*NormBleedMatrix(b,:,r);
    end
end

if o.Graphics
    figure(98043764); clf
    for i=1:nRounds
        subplot(ceil(nRounds/3),3,i); 
        imagesc(BleedMatrix(:,:,i)); 
        %caxis([0 1]); 
        title(sprintf('Cycle %d', o.UseRounds(i))); 
        set(gca, 'xtick', 1:nChans);
        set(gca, 'XTickLabel', o.bpLabels(o.UseChannels));
        set(gca, 'ytick', 1:nChans);
        set(gca, 'yTickLabel', o.bpLabels(o.UseChannels));
        if i==4
            xlabel('Actual')
            ylabel('Measured');
        end
    end
%     subplot(2,3,6);
%     caxis([0 1]); 
%     axis off
%     colormap hot
% %     colorbar
end

o.BleedMatrix = NormBleedMatrix;
o.pBleedMatrix = BleedMatrix;

%Load in anchor channel for each gene
GeneAnchorChannelInfo = load(o.GeneAnchorChannelFile);
o.GeneAnchorChannel =  GeneAnchorChannelInfo.T.Channel+1;

%% Get BledCodes for each gene
% now load in the code book and apply bleeds to it
%codebook_raw = importdata(o.CodeFile);
%CharCode = codebook_raw.textdata(2:end,5);
%GeneName = codebook_raw.textdata(2:end,3);
GeneName = {};
CharCode = {};
fp = fopen(o.CodeFile, 'r');
tmp = textscan(fp, '%s %s', inf);
GeneName=tmp{1};
CharCode=tmp{2};
fclose(fp);

% bit of a hack to get rid of Sst and Npy (assume always in the end)
nCodes = size(CharCode,1) - nnz(cellfun(@(v) strcmp(v(1:2),'SW'), CharCode));

% put them into object o but without the extras
o.CharCodes=CharCode(1:nCodes);
o.GeneNames=GeneName(1:nCodes);

% create numerical code (e.g. 33244 for CCGAA)
NumericalCode = zeros(nCodes, o.nRounds);
for r=1:o.nRounds
    if r<=o.nRounds-o.nRedundantRounds
        for c=1:nCodes
            [~, NumericalCode(c,r)] = ismember(CharCode{c}(r), o.bpLabels);
        end
    else
        % redundant round - compute codes automatically
        % find pseudobases for this code
        for c=1:nCodes
            PseudoCode = repmat('0',1,o.nRounds-o.nRedundantRounds);
            for p = 1:length(o.RedundantPseudobases)
                PseudoCode(1,ismember(CharCode{c}, o.RedundantPseudobases{p}))=('0'+p);
            end
            % now match them to the redundant codes
            for cc=1:nChans
                rrn = r-o.nRounds+o.nRedundantRounds;
                if ~isempty(regexp(PseudoCode, o.RedundantCodes{rrn,cc}, 'once'))
                    NumericalCode(c,r)=cc;
                end
            end
        end
    end
end

BledCodes = zeros(nCodes, o.nBP*o.nRounds);
UnbledCodes = zeros(nCodes, o.nBP*o.nRounds);
% make starting point using bleed vectors (means for each base on each day)
for i=1:nCodes
    for r=1:nRounds
        if any(o.UseChannels == NumericalCode(i,o.UseRounds(r))) == 0 continue; end
        BledCodes(i,o.UseChannels+o.nBP*(r-1)) = BleedMatrix(:, find(o.UseChannels == NumericalCode(i,o.UseRounds(r))), r);
        UnbledCodes(i,o.UseChannels(find(o.UseChannels == NumericalCode(i,o.UseRounds(r))))+o.nBP*(r-1)) = 1;
    end
end

o.pBledCodes = BledCodes;
o.UnbledCodes = UnbledCodes;

%% Find probability of each spot to each gene
%Comes from P(f) = P_lambda(lambda)P_hist(f-lambda*g) as f = lambda*g+background

%Load histogram data - background prob distribution
nBins = length(o.HistValues);
HistCounts = o.HistCounts;
nPixels = sum(HistCounts(:,1,1));
HistProbs = HistCounts/nPixels;    
o.HistProbs = (HistProbs+o.alpha)./(1+nBins*o.alpha);

%Get Lambda probability distribution for all genes
x = min(o.cSpotColors(:))-1:max(o.cSpotColors(:))-1;    %subsitiution x=lambda*g, -1 due to matlab indexing
o.ZeroIndex = find(x==0);     %need when looking up conv values
o.LambdaDist = zeros(length(x),nCodes,o.nBP,o.nRounds);

for GeneNo = 1:nCodes
    BledCode = reshape(o.pBledCodes(GeneNo,:),[o.nBP,o.nRounds]);
    numCharCode = str2double(regexp(cell2mat(o.CharCodes(GeneNo)),'\d','match'))+1;

    for b=1:o.nBP
        for r=1:o.nRounds
            g = BledCode(b,r);
            if numCharCode(r)==b
                %for b/r in CharCodes, expect non zero lambda.
                %g always >0 in this case
                o.LambdaDist(x>0,GeneNo,b,r) = raylpdf(x(x>0)/g,o.RaylConst)/g;
            else
                %for b/r not in CharCodes, expect approx zero lambda.
                o.LambdaDist(:,GeneNo,b,r) = (o.ExpConst/2)*exp(-o.ExpConst*abs(x/g))/abs(g);
            end
        end
    end
end

%Store convolution results as look up table
LookupTable = zeros(length(x),nCodes,o.nBP,o.nRounds);
fprintf('\nDoing convolutions for Channel  ');
for b=1:o.nBP
    fprintf('\b%d',b);
    for r=1:o.nRounds
        LookupTable(:,:,b,r)=log(conv2(o.LambdaDist(:,:,b,r),o.HistProbs(:,b,r),'same'));
    end
end
fprintf('\n');

%Get log probs for each spot 
nSpots = size(o.cSpotColors,1);
LogProb = zeros(nSpots,nCodes);

gChannelIndex = repmat(1:o.nBP,1,o.nRounds);
gRoundIndex = repelem(1:o.nRounds,1,o.nBP);
ChannelIndex = repmat(gChannelIndex,1,nCodes);
RoundIndex = repmat(gRoundIndex,1,nCodes);
GeneIndex = repelem(1:nCodes,1,o.nRounds*o.nBP);

for s=1:nSpots
    SpotIndex = repmat(o.ZeroIndex-1+o.cSpotColors(s,:),1,nCodes); %-1 due to matlab indexing I think
    Indices = sub2ind(size(LookupTable),SpotIndex,GeneIndex,ChannelIndex,RoundIndex);
    LogProb(s,:)=sum(reshape(LookupTable(Indices),[o.nRounds*o.nBP,nCodes]));
end
[LogProb,SpotCodeNo] = sort(LogProb,2,'descend');

o.pLogProb = LogProb(:,1);
o.pSpotCodeNo = SpotCodeNo(:,1);
o.pSpotScore = LogProb(:,1)-LogProb(:,2);
%Store deviation in spot scores - can rule out matches based on a low
%deviation.
o.pSpotScoreDev = std(LogProb,[],2);
o.pSpotIntensity = o.get_spot_intensity(o.pSpotCodeNo);
end