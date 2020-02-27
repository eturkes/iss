function SpotNo = iss_view_prob(o, FigNo, Norm, SpotNum)
    %This function lets you view the spot code, gene code,
    %spotcode-lambda*bledcode and ln(Prob) for the best match for a chosen
    %spot.
    %Now have option for different normalisations
    %Norm = 1: Raw Colors
    %Norm = 2: Normalised by percentile for each round and color channel
    %and then so each round has unit norm
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
        try
            zPlane = evalin('base', 'issPlot3DZPlane');     
            [~,SpotNo] = min(sum(abs(o.SpotGlobalYXZ-[xy(2),xy(1),zPlane]),2));
        catch
            [~,SpotNo] = min(sum(abs(o.SpotGlobalYXZ(:,1:2)-[xy(2),xy(1)]),2));
        end
    end
    CodeNo = o.pSpotCodeNo(SpotNo);
    
    %Different Normalisations
    if isempty(Norm) || Norm == 1
        cSpotColors = o.cSpotColors;
        cBledCodes = o.pBledCodes;
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
    BledCode = cBledCodes(CodeNo,:);
    ProbMatrix = get_prob_matrix(o,squeeze(o.cSpotColors(SpotNo,:,:)),CodeNo);
    
    try
        clf(430476533)
        figure(430476533)
    catch
        figure(430476533)
    end
    subplot(3,1,1);
    imagesc(MeasuredCode); colorbar
    caxis([0 max(MeasuredCode(:))]);
    title(sprintf('Spot %.0f Code',SpotNo));
    set(gca, 'ytick', 1:o.nBP);
    set(gca, 'YTickLabel', o.bpLabels);
    ylabel('Color Channel');
    
    subplot(3,1,2)
    imagesc(reshape(BledCode, CodeShape)); colorbar
    %caxis([0 max(cBledCode(:))]);
    title(sprintf('Predicted Code for %s, code #%d', o.GeneNames{CodeNo}, CodeNo));
    set(gca, 'ytick', 1:o.nBP);
    set(gca, 'YTickLabel', o.bpLabels);
    ylabel('Color Channel');
    xlabel('Round');
    
    
    subplot(3,1,3);
    imagesc(ProbMatrix); colorbar
    %caxis([min(ProbMatrix(:)) max(ProbMatrix(:))]);
    set(gca, 'ytick', 1:o.nBP);
    set(gca, 'YTickLabel', o.bpLabels);
    ylabel('Color Channel');
    xlabel('Round');
    title(sprintf('Log Probability'));
    
    %Color different parameters depending if over threshold
    if o.pSpotScore(SpotNo)>o.pScoreThresh
        c1 = [0,0.7,0]; else; c1 = [0,0,0];end
    if o.pLogProb(SpotNo)<o.pLogProbThresh
        c2 = [1,0,0]; else; c2 = [0,0,0];end
    if o.pSpotScore(SpotNo)+o.pSpotScoreDev(SpotNo)<o.pDevThresh
        c3 = [1,0,0]; else; c3 = [0,0,0];end
    if o.pSpotIntensity(SpotNo)<o.pIntensityThresh
        c4 = [1,0,0]; else; c4 = [0,0,0];end

    
    set(gcf,'Position',[350 100 1000 850])
    figtitle = sgtitle('', 'interpreter', 'tex');   %'tex' required for colors
    if isempty(o.cSpotAnchorChannel)
        figtitle.String = sprintf('%s{%f %f %f}Score = %.1f, %s{%f %f %f}LogProb = %.0f, %s{%f %f %f}Score Deviation = %.1f, %s{%f %f %f}Intensity = %.0f',...
            '\color[rgb]',c1,o.pSpotScore(SpotNo),'\color[rgb]',c2, o.pLogProb(SpotNo),...
            '\color[rgb]',c3,o.pSpotScoreDev(SpotNo),'\color[rgb]',c4,o.pSpotIntensity(SpotNo));
    else
        if o.cSpotAnchorChannel(SpotNo) ~= o.GeneAnchorChannel(CodeNo)
            c5 = [1,0,0]; else; c5 = [0,0,0];end
        figtitle.String = sprintf('%s{%f %f %f}Score = %.1f, %s{%f %f %f}LogProb = %.0f, %s{%f %f %f}Score Deviation = %.1f, %s{%f %f %f}Intensity = %.0f, %s{%f %f %f}Spot SplitAnchor Channel = %.0f, %s SplitAnchor Channel = %.0f',...
            '\color[rgb]',c1,o.pSpotScore(SpotNo),'\color[rgb]',c2, o.pLogProb(SpotNo),...
            '\color[rgb]',c3,o.pSpotScoreDev(SpotNo),'\color[rgb]',c4,o.pSpotIntensity(SpotNo),...
            c5,o.cSpotAnchorChannel(SpotNo),o.GeneNames{CodeNo},o.GeneAnchorChannel(CodeNo));
    end

    drawnow
    
    fprintf('Spot %d at yxz=(%d,%d,%d): code %d, %s\n', ...
        SpotNo, o.SpotGlobalYXZ(SpotNo,1),o.SpotGlobalYXZ(SpotNo,2),...
        round(o.SpotGlobalYXZ(SpotNo,3)), CodeNo, o.GeneNames{CodeNo});
end
