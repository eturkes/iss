function QualOK = quality_threshold_prob(o)
% QualOK = o.quality_threshold
% quick function that returns a binary saying which spots are above quality
% threshold

QualOK = (o.pSpotScore>o.pScoreThresh | ...
o.pSpotIntensity>o.pIntensityThresh & o.pLogProb>o.pLogProbThresh & o.pSpotScore+o.pSpotScoreDev>o.pDevThresh);

%Ensure spot was detected on correct channel
if ~isempty(o.cSpotAnchorChannel)
    CorrectSpotAnchorChannel = o.GeneAnchorChannel(o.SpotCodeNo);
    QualOK = (QualOK & o.cSpotAnchorChannel == CorrectSpotAnchorChannel);
end

% % HACK ALERT
% QualOK = QualOK & o.cSpotIsolated;

nCombiCodes = sum(~strcmp(o.CharCodes, 'EXTRA'));

% now extras - they have their own thresholds, set manually for each type
for i=1:size(o.ExtraCodes,1)
    MySpots = (o.SpotCodeNo == nCombiCodes+i);
    QualOK(MySpots) = o.SpotIntensity(MySpots)>o.ExtraCodes{i,4};
end