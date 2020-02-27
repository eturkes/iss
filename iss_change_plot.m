function iss_change_plot(o,Method)
%Given issPlot3DObject, this function lets you change the details
%of the plot without closing the figure e.g. you can change
%o.CombiQualThresh or issPlot3DObject.ZThick to change the threshold value and the number
%of Z planes you see. 
%If Method == 'Prob', this changes the gene assignments to those given by
%the probability method rather than the dot product. In this case
%o.pScoreThresh is the threshold.

S = evalin('base', 'issPlot3DObject');
figure(S.FigNo);
h = findobj('type','line'); %KEY LINES: DELETE EXISTING SCATTER PLOTS SO CHANGE_SYMBOLS WORKS
delete(h);
issPlot3DZPlane = evalin('base', 'issPlot3DZPlane');  %Curent Zplane


if strcmpi('Prob',Method)
    S.SpotGeneName = o.GeneNames(o.pSpotCodeNo);
    S.uGenes = unique(S.SpotGeneName);
    % which ones pass quality threshold (combi first)
    S.QualOK = o.quality_threshold_prob;
else
    S.SpotGeneName = o.GeneNames(o.SpotCodeNo);
    S.uGenes = unique(S.SpotGeneName);
    % which ones pass quality threshold (combi first)
    S.QualOK = o.quality_threshold;
end
%S.SpotYXZ = o.SpotGlobalYXZ;
%S.Roi is the Roi for the current Z plane
S.Roi(5:6) = [issPlot3DZPlane-S.ZThick,issPlot3DZPlane+S.ZThick];
InRoi = all(int64(round(S.SpotYXZ))>=S.Roi([3 1 5]) & round(S.SpotYXZ)<=S.Roi([4 2 6]),2);
PlotSpots = find(InRoi & S.QualOK);
[~, S.GeneNo] = ismember(S.SpotGeneName(PlotSpots), S.uGenes);
S.h = zeros(size(S.uGenes));

%hold on; GET RID OF HOLD AND JUST DELETE PLOTS USING DELETE MEANS THAT THE
%ZOOM IS KEPT BETWEEN Z PLANES
for i=1:length(S.uGenes)
    MySpots = PlotSpots(S.GeneNo==i);
    if any(MySpots)
        S.h(i) = plot(S.SpotYXZ(MySpots,2), S.SpotYXZ(MySpots,1), '.');
    end
end 
%hold off

legend(S.h(S.h~=0), S.uGenes(S.h~=0));
legend off;

set(gca, 'Clipping', 'off');

if ~isempty(PlotSpots)
    change_gene_symbols(0);
else
    set(gcf, 'color', 'k');
    set(gcf, 'InvertHardcopy', 'off');    
end

assignin('base','issPlot3DZPlane',S.Roi(5)+S.ZThick)
assignin('base','issPlot3DObject',S)