function SpotNo = iss_view_codes(o, FigNo, Norm, SpotNum)
    %Now have option for different normalisations
    %Norm = 1: Raw Colors
    %Norm = 2: Normalised by o.SpotNormPrctile in each colour channel and round,
    %then if o.CallSpotsCodeNorm='WholeCode', normalise so whole code has L2 norm
    %of 1 but if o.CallSpotsCodeNorm='Round', normalise so each round has L2 norm of 1.
    %Norm = 3: Normalised by percentile for each color channel across all
    %rounds

    
    if nargin>=4
        SpotNo = SpotNum;
    else
        if nargin>=2
            figure(FigNo);
        end
        CrossHairColor = [1,1,1];   %Make white as black background        
        xy = ginput_modified(1,CrossHairColor);        
        S = evalin('base', 'issPlot3DObject');  
        InRoi = all(int64(round(S.SpotYXZ))>=S.Roi([3 1 5]) & round(S.SpotYXZ)<=S.Roi([4 2 6]),2);
        PlotSpots = find(InRoi & S.QualOK);         %Only consider spots that can be seen in current plot
        [~,SpotIdx] = min(sum(abs(o.SpotGlobalYXZ(PlotSpots,1:2)-[xy(2),xy(1)]),2));
        SpotNo = PlotSpots(SpotIdx);
    end
    CodeNo = o.SpotCodeNo(SpotNo);
    
    if nargin<3 || isempty(Norm)
        Norm = 1;
    end
    
    %Different Normalisations
    if isempty(Norm) || Norm == 1
        cSpotColors = o.cSpotColors;
        cBledCodes = o.pBledCodes;
        if isempty(o.pBledCodes)
            cBledCodes = o.BledCodes;
        end
    elseif Norm == 2
        cSpotColors = o.cNormSpotColors;
        cBledCodes = o.NormBledCodes;
    elseif Norm == 3
        cSpotColors = o.cSpotColors;
        NewBleedMatrix = o.pBleedMatrix;
        for b = 1:o.nBP
            bSpotColors = o.cSpotColors(:,b,:);
            p = prctile(bSpotColors(:), o.SpotNormPrctile);
            cSpotColors(:,b,:) = cSpotColors(:,b,:)/p;
            NewBleedMatrix(b,:,:) = o.pBleedMatrix(b,:,:)/p;                        
        end
        cBledCodes = change_bled_codes(o,NewBleedMatrix);
    end
            
    
    MeasuredCode = squeeze(cSpotColors(SpotNo,:,:));
    CodeShape = size(MeasuredCode);
    
    figure(930476530)
    subplot(2,1,1);
    imagesc(MeasuredCode); colorbar
    caxis([0 max(MeasuredCode(:))]);
    title(sprintf('Measured code: match %.3f to %s', o.SpotScore(SpotNo), o.GeneNames{CodeNo}));
    
    set(gca, 'ytick', 1:o.nBP);
    set(gca, 'YTickLabel', o.bpLabels);
    ylabel('Color Channel');
    
    subplot(2,1,2)
    BledCode = cBledCodes(CodeNo,:);
    imagesc(reshape(BledCode, CodeShape)); colorbar
    caxis([0 max(BledCode(:))]);

    title(sprintf('Predicted Code for %s, code #%d', o.GeneNames{CodeNo}, CodeNo));
    
    
    set(gca, 'ytick', 1:o.nBP);
    set(gca, 'YTickLabel', o.bpLabels);
    ylabel('Color Channel');
    xlabel('Round');

    fprintf('Spot %d at yxz=(%d,%d,%d): code %d, %s\n', ...
        SpotNo, o.SpotGlobalYXZ(SpotNo,1),o.SpotGlobalYXZ(SpotNo,2),round(o.SpotGlobalYXZ(SpotNo,3)), CodeNo, o.GeneNames{CodeNo});

    
end
    