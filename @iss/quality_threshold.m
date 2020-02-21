function QualOK = quality_threshold(o)
% QualOK = o.quality_threshold
% quick function that returns a binary saying which spots are above quality
% threshold.
% If are using split channels, also ensures all spots detected on correct
% anchor channel.

QualOK = (o.SpotCombi & o.SpotScore>o.CombiQualThresh & o.SpotIntensity>o.CombiIntensityThresh & o.SpotScoreDev>o.CombiDevThresh);

if ~isempty(o.cSpotAnchorChannel)
    CorrectSpotAnchorChannel = o.GeneAnchorChannel(o.SpotCodeNo);
    %QualOK = (QualOK & o.cSpotAnchorChannel == CorrectSpotAnchorChannel);
end

% % HACK ALERT
% QualOK = QualOK & o.cSpotIsolated;

nCombiCodes = sum(~strcmp(o.CharCodes, 'EXTRA'));

% now extras - they have their own thresholds, set manually for each type
for i=1:size(o.ExtraCodes,1)
    MySpots = (o.SpotCodeNo == nCombiCodes+i);
    QualOK(MySpots) = o.SpotIntensity(MySpots)>o.ExtraCodes{i,4};
end