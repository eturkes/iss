function iss_view_codes(o, FigNo)
% iss_view_codes(o, FigNo)
% 
% run this after iss.plot, then left-click on a gene read to see its actual 
% code value and the bled template for this code. When you are done, press
% any key and it will quit.
%
% note it is NOT a member function, you need to provide o as an argument.

    if nargin>=2
        figure(FigNo);
    else
        FigNo = gcf;
    end
    
    % use normed SpotColors that are actually used to determine spot scores
    cSpotColors = o.cNormSpotColors;
    
    while 1
        figure(FigNo);
        set(gca, 'Color', [1 1 1]*.2);
        [x, y, button] = ginput(1);
        if button~=1
            break;
        end
        set(gca, 'color', 'k');

        %If specify z plane in workspace as should be from plot3D, find spots closest to (y,x,ZPlane)
        try
            zPlane = evalin('base', 'issPlot3DZPlane');     
            [~,SpotNo] = min(sum(abs(o.SpotGlobalYXZ-[y,x,zPlane]),2));
        catch
            [~,SpotNo] = min(sum(abs(o.SpotGlobalYXZ(:,1:2)-[y,x]),2));
        end
        
        CodeNo = o.SpotCodeNo(SpotNo);

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
        cBledCode = o.BledCodes(CodeNo,:);
        imagesc(reshape(cBledCode, CodeShape)); colorbar
        caxis([0 max(cBledCode(:))]);

        title(sprintf('Predicted Code for %s, code #%d (%s)', o.GeneNames{CodeNo}, CodeNo, o.CharCodes{CodeNo}));


        set(gca, 'ytick', 1:o.nBP);
        set(gca, 'YTickLabel', o.bpLabels);
        ylabel('Color Channel');
        xlabel('Round');

        fprintf('Spot %d at yxz=(%d,%d,%d): code %d, %s\n', ...
            SpotNo, o.SpotGlobalYXZ(SpotNo,1),o.SpotGlobalYXZ(SpotNo,2),round(o.SpotGlobalYXZ(SpotNo,3)), CodeNo, o.GeneNames{CodeNo});

    end
end
    