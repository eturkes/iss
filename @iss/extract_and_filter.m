function o = extract_and_filter(o)
% create tiff files for each tile that are top-hat filtered versions of
% original czi files

    o.TileFiles = cell(o.nRounds+o.nExtraRounds,1,1,1); % 1,1,1 because we don't yet know how many tiles
    
    %New filter
    h = -hanning(o.ExtractR2*2+1);
    h = -h/sum(h);
    h(o.ExtractR2+1-o.ExtractR1:o.ExtractR2+1+o.ExtractR1) = ...
        h(o.ExtractR2+1-o.ExtractR1:o.ExtractR2+1+o.ExtractR1)+hanning(o.ExtractR1*2+1)/sum(hanning(o.ExtractR1*2+1));
    SE = ftrans2(h');
%     h2D = ftrans2(h');
%     hzdirection = hanning(3);
%     hzdirection = reshape(hzdirection,[1,1,3]);
%     SE = h2D.*hzdirection;
    
    
    for r = 1:o.nRounds+o.nExtraRounds       
        imfile = fullfile(o.InputDirectory, [o.FileBase{r}, o.RawFileExtension]);

        % construct a Bio-Formats reader with the Memoizer wrapper
        bfreader = loci.formats.Memoizer(bfGetReader(), 0);
        % initiate reader
        bfreader.setId(imfile);

        % get some basic image metadata
        [nSeries, nSerieswPos, nChannels, o.nZ, xypos, o.XYpixelsize,o.Zpixelsize] = ...
            get_ome_tilepos(bfreader);
        if isempty(xypos) || size(xypos, 1)==1
            if r == 1
                warning('first round xypos empty - using values from initial manual input')
                assert(~isempty(o.TileInitialPosXY), 'xypos unavailable')
                xypos = o.TileInitialPosXY;
                xyposOld = xypos;
            else
                warning('xypos empty - using values from previous round')
                xypos = xyposOld;
            end
            nSerieswPos = size(xypos,1);
        else
            xyposOld = xypos;
        end
        
        scene = nSeries/nSerieswPos;

        bfreader.close();
        
        if r == 1
            if isempty(o.AutoThresh)
                o.AutoThresh = zeros(nSerieswPos,nChannels,o.nRounds+o.nExtraRounds);  
            end
            
            % find x and y grid spacing as median of distances that are about
            % right
            dx = xypos(:,1)-xypos(:,1)'; % all pairs of x distances
            xStep = median(dx(abs(1- dx(:)/o.MicroscopeStepSize)<.5));
            dy = xypos(:,1)-xypos(:,1)'; % all pairs of y distances
            yStep = median(dy(abs(1- dy(:)/o.MicroscopeStepSize)<.5));
        
        
            % find coordinates for each tile
            o.TileInitialPosYX = fliplr(1+round((xypos - min(xypos))./[xStep yStep]));
            TilePosYX = o.TileInitialPosYX;
            %Below is a safeguard incase wrong positions found - can do
            %this as we know what the answer should be.
            MaxY = max(TilePosYX(:,1));
            MaxX = max(TilePosYX(:,2));
            if MaxY*MaxX ~= nSeries
                warning('Number of tiles (%d) is not equal to maximum Y position (%d) multiplied by maximum X position (%d)'...
                    , nSeries, MaxY, MaxX)
                break
            else
                TilePosY = flip(repelem(1:MaxY,MaxX));
                TilePosYX(:,1) = TilePosY;
                TilePosX = repmat([flip(1:MaxX),1:MaxX],1,ceil(MaxY/2));
                TilePosYX(1:nSeries,2) = TilePosX(1:nSeries);
            end
        end
        
        o.TilePosYXC = zeros(nSerieswPos*nChannels,3);

        % set up filename grid for this round
        fName = cell(nSerieswPos*nChannels,1);
        
        Index = 1;
        %parfor t = 1:nSerieswPos  
        for t = 1:nSerieswPos  
                       
            % a new reader per worker
            bfreader = javaObject('loci.formats.Memoizer', bfGetReader(), 0);
            % use the memo file cached before
            bfreader.setId(imfile);

            bfreader.setSeries(scene*t-1);
            for c = 1:nChannels
                tic
                fName{Index} = fullfile(o.TileDirectory, ...
                    [o.FileBase{r}, '_t', num2str(t),'c', num2str(c), '.tif']);  
                
                if exist(fName{Index}, 'file')
                    fprintf('Round %d tile %d already done.\n', r, t);
                    o.TilePosYXC(Index,:) = [TilePosYX(t,:),c];          %Think first Z plane is the highest
                    o.TileFiles{r,o.TilePosYXC(Index,1), o.TilePosYXC(Index,2),o.TilePosYXC(Index,3)} = fName{Index};
                    if o.AutoThresh(t,c,r) == 0  
                        if c == o.DapiChannel && r == o.ReferenceRound; continue; end
                        IFS = o.load_3D(r,o.TilePosYXC(Index,1),o.TilePosYXC(Index,2),c)-o.TilePixelValueShift;
                        o.AutoThresh(t,c,r) = median(abs(IFS(:)))*o.AutoThreshMultiplier;
                    end
                    Index = Index+1;
                    continue;
                elseif (r == o.ReferenceRound && c ~= o.AnchorChannel) && (r == o.ReferenceRound && c ~= o.DapiChannel)
                    %Only need anchor and dapi tiles in reference round
                    continue;
                end
                                                                        
                %TopHat SE
%                 if c == o.DapiChannel && r == o.ReferenceRound    
%                         %SE = strel3D_2(20,10);       % I.e. set to 8 microns for DAPI
%                         SE = get_3DSE(o.DapiR1YX,o.DapiR1Z,o.DapiR2YX,o.DapiR2Z);       
%                 else
%                         %SE = strel3D_2(3,3);    %I.e. Set to 1 micron
%                         SE = get_3DSE(o.ExtractR1YX,o.ExtractR1Z,o.ExtractR2YX,o.ExtractR2Z);
%                 end

                I = zeros(o.TileSz,o.TileSz,o.nZ); 
                for z = 1:o.nZ
                    iPlane = bfreader.getIndex(z-1, c-1, 0)+1;
                    I(:,:,z) = bfGetPlane(bfreader, iPlane);
                end       
                
                %Make noise white first by divding amplitude of FT
                
                %FT = fftn(I);
                %Norm_FT = FT ./ abs(FT);
                %filter = fspecial3('gaussian',size(I),2);%DESCRIBE BETTER!!!!
                %Shiftfilter = fftshift(filter);     %Shift for FT so centered on 0
                %FT_filter = fftn(Shiftfilter);          
                %NormFT_filter = FT_filter ./ abs(FT);
                %Final_FT = Norm_FT .* NormFT_filter;
                %IFS = ifftn(Final_FT);
                
                %I = ifftn(Norm_FT);
                I = padarray(I,(size(SE)-1)/2,'replicate','both');
                IFS = convn(I,SE,'valid');    
                
                %Scaling so fills uint16 range.
                if c == o.DapiChannel && r == o.ReferenceRound  
                    if strcmpi(o.DapiScale, 'auto')
                        o.DapiScale = 10000/max(max(max(IFS)));
                    end
                    IFS = IFS*o.DapiScale;
                else
                    %Finds o.ExtractScale from first image and uses this
                    %value for the rest
                    if strcmpi(o.ExtractScale, 'auto')
                        o.ExtractScale = 10000/max(max(max(IFS)));
                    end
                    IFS = IFS*o.ExtractScale;
                    
                    %Determine auto thresholds
                    o.AutoThresh(t,c,r) = median(abs(IFS(:)))*o.AutoThreshMultiplier;
                end
                
                
                
                %Append each z plane to same tiff image
                %Add 2^16/2 so keep negative pixels for background analysis
                IFS = IFS + o.TilePixelValueShift;
                for z = 1:o.nZ
                    imwrite(uint16(IFS(:,:,z)),...  %Not sure if uint16 is correct, wasnt working without
                            fullfile(o.TileDirectory,...
                            [o.FileBase{r}, '_t', num2str(t),'c', num2str(c), '.tif']),...
                            'tiff', 'writemode', 'append');
                end

                o.TilePosYXC(Index,:) = [TilePosYX(t,:),c];          %Think first Z plane is the highest
                o.TileFiles{r,o.TilePosYXC(Index,1), o.TilePosYXC(Index,2),o.TilePosYXC(Index,3)} = fName{Index};
                fprintf('Round %d tile %d colour channel %d finished.\n', r, t, c);                                               
                Index = Index+1; 
                toc
            end
            bfreader.close();
            
        end
        
    
    o.EmptyTiles = cellfun(@isempty, squeeze(o.TileFiles(o.ReferenceRound,:,:,1)))*0;

    end
    
    %Plot boxplots showing distribution af AutoThresholds
    if o.Graphics
        Thresholds = [];
        group = [];
        index = 1;
        for c=1:nChannels
            for r=1:o.nRounds
                Thresholds = [Thresholds;o.AutoThresh(:,c,r)];
                group = [group;index*ones(size(o.AutoThresh(:,1,1)))];
                index = index+1;
            end
        end
        %Add anchor
        Thresholds = [Thresholds;o.AutoThresh(:,o.AnchorChannel,o.ReferenceRound)];
        group = [group;index*ones(size(o.AutoThresh(:,1,1)))];
        
        figure(43290);
        colors = colormap(lines(nChannels));
        Colors = repelem(colors,o.nRounds,1);
        Colors = [Colors;repelem([0,0,0],nChannels,1)];        
        boxplot(Thresholds,group,'Colors',Colors, 'plotstyle', 'compact','labels', [string(repmat(1:o.nRounds,1,nChannels)),'Anchor']);
        set(gca,'TickLength',[0 0])
        ylabel('AutoThreshold');
        xlabel('Round');
        hold on
        for c=1:nChannels
            plot(NaN,1,'color', colors(c,:), 'LineWidth', 4);
        end
        leg = legend(o.bpLabels);
        title(leg,'Color Channel');
        hold off
    end
end

function SE = get_3DSE(r1YX,r1Z,r2YX,r2Z)
    % structuring element for convlolution filtering
    % Positive inner circle radius r1 and negative outer annulus radius r2. Overall sums to
    % zero. 
    SE = zeros(r2YX*2+1,r2YX*2+1,r2Z*2+1);
    SE(r2YX+1-r1YX:r2YX+1+r1YX,r2YX+1-r1YX:r2YX+1+r1YX,r2Z+1-r1Z:r2Z+1+r1Z) = fspecial3('ellipsoid',[r1YX,r1YX,r1Z]);
    SE = SE - fspecial3('ellipsoid',[r2YX,r2YX,r2Z]);
end