function iss_change_plot_contour(o,Method,GenesToShow,UseSpots)
%% iss_change_plot(o,Method,GenesToShow,UseSpots)
%
% Given issPlot3DObject, this function lets you change the details
% of the plot without closing the figure e.g. you can change
% o.CombiQualThresh, o.pIntensityThresh, o.pScoreThresh, o.pScoreThresh2,
% o.pLogProbThresh
%
% o: iss object
% Method: 'DotProduct','Prob' or 'Pixel' to consider gene assignments given
% by o.SpotCodeNo, o.pSpotCodeNo and o.pxSpotCodeNo respectively.
% GenesToShow: cell of gene names that you want to see e.g.
% [{'Npy'},{'Pvalb'}]. It is case sensitive.
% UseSpots: if you want to use your own thresholding, not
% o.quality_threshold. Logical array e.g. o.pSpotScore>0 or
% get_gene_clusters(o)

%%
S = evalin('base', 'issPlot2DObject');
figure(S.FigNo);
h = findobj('type','line'); %KEY LINES: DELETE EXISTING SCATTER PLOTS SO CHANGE_SYMBOLS WORKS
delete(h);


if nargin<3 || isempty(GenesToShow)
    GenesToShow = o.GeneNames;
    if ~isfield(S,'GeneNoToShow')
        %Only change if not previosuly given GenesToShow
        S.GeneNoToShow = find(ismember(o.GeneNames,GenesToShow));
    end
else
    S.GeneNoToShow = find(ismember(o.GeneNames,GenesToShow));
end


if nargin<2 || isempty(Method)  
    if strcmpi(S.CallMethod,'DotProduct')
        if nargin>=4 && length(UseSpots)==length(o.SpotScore) && islogical(UseSpots)
            S.QualOK = UseSpots & ismember(o.SpotCodeNo,S.GeneNoToShow);
        else
            if nargin>=4; warning('UseSpots not valid, using o.quality_threshold');end
            S.QualOK = o.quality_threshold(S.CallMethod) & ismember(o.SpotCodeNo,S.GeneNoToShow);
        end
    elseif strcmpi(S.CallMethod,'Prob')
        if nargin>=4 && length(UseSpots)==length(o.pSpotScore) && islogical(UseSpots)
            S.QualOK = UseSpots & ismember(o.pSpotCodeNo,S.GeneNoToShow);    
        else
            if nargin>=4; warning('UseSpots not valid, using o.quality_threshold');end
            S.QualOK = o.quality_threshold(S.CallMethod) & ismember(o.pSpotCodeNo,S.GeneNoToShow);
        end
    elseif strcmpi(S.CallMethod,'Pixel')
        if nargin>=4 && length(UseSpots)==length(o.pxSpotScore) && islogical(UseSpots)
            S.QualOK = UseSpots & ismember(o.pxSpotCodeNo,S.GeneNoToShow);    
        else
            if nargin>=4; warning('UseSpots not valid, using o.quality_threshold');end
            S.QualOK = o.quality_threshold(S.CallMethod) & ismember(o.pxSpotCodeNo,S.GeneNoToShow);
        end
    end
else
    if strcmpi('Prob',Method)
        S.CallMethod = 'Prob';
        S.SpotGeneName = o.GeneNames(o.pSpotCodeNo);
        S.uGenes = unique(S.SpotGeneName);
        if nargin>=4 && length(UseSpots)==length(o.pSpotScore) && max(UseSpots)==1
            S.QualOK = UseSpots & ismember(o.pSpotCodeNo,S.GeneNoToShow);  
        else
            % which ones pass quality threshold (combi first)
            if nargin>=4; warning('UseSpots not valid, using o.quality_threshold_prob');end
            S.QualOK = o.quality_threshold(S.CallMethod) & ismember(o.pSpotCodeNo,S.GeneNoToShow);
        end
    elseif strcmpi('Pixel',Method)
        S.CallMethod = 'Pixel';
        S.SpotGeneName = o.GeneNames(o.pxSpotCodeNo);
        S.uGenes = unique(S.SpotGeneName);
        if nargin>=4 && length(UseSpots)==length(o.pxSpotScore) && max(UseSpots)==1
            S.QualOK = UseSpots & ismember(o.pxSpotCodeNo,S.GeneNoToShow);  
        else
            % which ones pass quality threshold (combi first)
            if nargin>=4; warning('UseSpots not valid, using o.quality_threshold_prob');end
            S.QualOK = o.quality_threshold(S.CallMethod) & ismember(o.pxSpotCodeNo,S.GeneNoToShow);
        end
    else
        S.CallMethod = 'DotProduct';
        S.SpotGeneName = o.GeneNames(o.SpotCodeNo);
        S.uGenes = unique(S.SpotGeneName);
        if nargin>=4 && length(UseSpots)==length(o.SpotScore) && max(UseSpots)==1
            S.QualOK = UseSpots & ismember(o.SpotCodeNo,S.GeneNoToShow);    
        else
            % which ones pass quality threshold (combi first)
            if nargin>=4; warning('UseSpots not valid, using o.quality_threshold');end
            S.QualOK = o.quality_threshold(S.CallMethod) & ismember(o.SpotCodeNo,S.GeneNoToShow);
        end
    end
end

if strcmpi('Pixel',Method)
    S.SpotYX = o.pxSpotGlobalYX;
else
    S.SpotYX = o.SpotGlobalYX;
end
        
InRoi = all(int64(round(S.SpotYX))>=S.Roi([3 1]) & round(S.SpotYX)<=S.Roi([4 2]),2);
PlotSpots = find(InRoi & S.QualOK);
[~, S.GeneNo] = ismember(S.SpotGeneName(PlotSpots), S.uGenes);
S.h = zeros(size(S.uGenes));

%hold on; GET RID OF HOLD AND JUST DELETE PLOTS USING DELETE MEANS THAT THE
%ZOOM IS KEPT BETWEEN Z PLANES
for i=1:length(S.uGenes)
    MySpots = PlotSpots(S.GeneNo==i);
    if any(MySpots)
        S.h(i) = plot(S.SpotYX(MySpots,2), S.SpotYX(MySpots,1), '.');
    end
end 
%hold off

legend(S.h(S.h~=0), S.uGenes(S.h~=0));
legend off;

set(gca, 'Clipping', 'off');

if ~isempty(PlotSpots)
    change_gene_symbols_individual(2, 6);
else
    set(gcf, 'color', 'k');
    set(gcf, 'InvertHardcopy', 'off');    
end

spot_no = find(ismember(o.GeneNames, GenesToShow));
thresh = o.quality_threshold('Pixel') & ismember(o.pxSpotCodeNo, spot_no);

spots = o.pxSpotGlobalYX(thresh,:);
codes = o.pxSpotCodeNo(thresh,:);
subset = size(spots);
rng(1);
subset = randsample(subset(1), subset(1));
sub_spots = spots(subset,:);

% The following section is adapted from content on Stack Overflow.
% https://stackoverflow.com/questions/9134014/contour-plot-coloured-by-clustering-of-points-matlab
% Question asked by: HCAI https://stackoverflow.com/users/1134241/hcai
% Answer given by: Vidar https://stackoverflow.com/users/346645/vidar
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Border = 0;
Sigma = 500;
stepSize = 100;

X=(sub_spots(:,2) / 1)';
Y=(sub_spots(:,1) / 1)';
D = [X' Y'];
N = length(X);

Xrange = [min(X)-Border max(X)+Border];
Yrange = [min(Y)-Border max(Y)+Border];

% Set up coordinate grid.
%%%%%%%%%%%%%%%%%%%%%%%%%
[XX, YY] = meshgrid(Xrange(1):stepSize:Xrange(2), Yrange(1):stepSize:Yrange(2));
YY = flipud(YY);
%%%%%%%%%%%%%%%%%%%%%%%%%

% Parzen parameters and function handle.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pf1 = @(C1,C2) (1/N)*(1/((2*pi)*Sigma^2)).* ...
         exp(-( (C1(1)-C2(1))^2+ (C1(2)-C2(2))^2)/(2*Sigma^2));

PPDF1 = zeros(size(XX));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Populate coordinate surface.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[R, C] = size(PPDF1);
NN = length(D);
for c = 1:C
   for r = 1:R
       for d = 1:N
            PPDF1(r,c) = PPDF1(r,c) + ...
                pf1([XX(1,c) YY(r,1)], [D(d,1) D(d,2)]);
       end
   end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Normalize data.
%%%%%%%%%%%%%%%%%
m1 = max(PPDF1(:));
PPDF1 = PPDF1 / m1;
%%%%%%%%%%%%%%%%%

% Add PDF estimates to figure.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
orig_ax = gca;
contour_ax = axes;
[~, contour_h] = contourf(XX, YY, PPDF1);
colormap(gca, 'hot')
set(gca, 'YDir', 'reverse');
set(gca, 'XDir', 'reverse');
set(gca, 'color', 'none');
set(gca,'XTick',[], 'YTick', []);

drawnow; % Required for the below to work.
fills_h = contour_h.FacePrims;
[fills_h.ColorType] = deal('truecoloralpha');
for i = 1:numel(fills_h)
   fills_h(i).ColorData(4) = 100;
end

set(gcf, 'currentaxes', orig_ax);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

assignin('base','issPlot2DObject',S)
