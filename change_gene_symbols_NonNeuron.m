function New_symbols = change_gene_symbols_NonNeuron(MarkerSize, FontSize, MultiCol)
% ChangeGeneSymbols(MarkerSize, FontSize, nPerCol);
%
% changes gene symbols so in situ plots look nice. 
% MarkerSize defaults to 6 - if 0, won't change existing sizes
%
% FontSize is font size for legend
%
% nPerCol says how many legend entries per column (0 for one-column)
% 
% Kenneth D. Harris, 29/3/17
% GPL 3.0 https://www.gnu.org/licenses/gpl-3.0.en.html
 
if nargin<1 || isempty(MarkerSize)
    MarkerSize = 6;
end
MarkerSize = 6;
if nargin<2 || isempty(FontSize)
    FontSize = 5;
end

if nargin<3
    MultiCol = 49;
end

%Want 34 distinguishable RGB colors, that will show up on black[0,0,0]
%and that are not grey [0.0667,0.0667,0.0667] as these are neurons.
BackgroundColors = [0,0,0;0.0667,0.0667,0.0667;1,1,1];
Colors = distinguishable_colors(33,BackgroundColors);

% Colors based on Class Assignments

% blood vessels etc
VECA = hsv2rgb([0.02, 1, 1]); %vascular endothelial cells arterial
VECV = hsv2rgb([0.98, 1, 1]); %vascular endothelial cells venous
VLMC = hsv2rgb([.95, .5, .6]); %vascular leptomengial cells
ABC = hsv2rgb([.95, .5, .6]); %more  vascular leptomengial cells
PER = hsv2rgb([0, .5, 1]); % pericytes
VSMCA = hsv2rgb([.0, .8, .6]); % vascular smooth muscle

% microglia
MGL = hsv2rgb([.1 1 .6]);
PVM = hsv2rgb([.1 1 .6]); %Colors(10,:);

% oligodendrocytes
COP = hsv2rgb([.6, .2, .5]); %committed oligo precursors
NFOL = hsv2rgb([.6, .2, .5]); %newly formed oligos (same color as different symbols)
MOL = hsv2rgb([.65, .1, .75]); % mature oligos
MOL2 = hsv2rgb([.7, 0, 1]); % more mature oligos
MFOL = hsv2rgb([.7, 0, 1]); %myelin forming oligos (same color as different symbols)

% astrocytes
ACMB = hsv2rgb(.5,1,1); % all same color, different symbols
ACBG = hsv2rgb(.5,1,1); % 
ACTE = hsv2rgb(.5,1,1); % 
ACOB = hsv2rgb(.5,1,1);

% choroid
CHOR = hsv2rgb([ .13 .6 .8]);
% ependymal
EPEN = hsv2rgb([ .19 1 1]);
EPSC = hsv2rgb([ .19 1 1]);




% neurons 1
CBPC = hsv2rgb([1/3 1 1]);
MSN = hsv2rgb([1/3 1 1]);
HBGLU = hsv2rgb([1/3 1 1]);
HBSER =hsv2rgb([1/3 1 1]);
MBDOP = Colors(30,:);

% neurons 2
HBINH = hsv2rgb([.27 1 .7]);
HBCHO = hsv2rgb([.27 1 .7]);
DECHO = hsv2rgb([.27 1 .7]);

% neuroblasts
DGNBL = hsv2rgb([.4 .5 .5]);
OBNBL = hsv2rgb([.4 .5 .5]);
SEPNBL = hsv2rgb([.4 .5 .5]);
SZNBL = hsv2rgb([.4 .5 .5]);
DETPH = hsv2rgb([.4 .5 .5]);

% sattelite glia and ensheathing cells
SATG = hsv2rgb([2/3 1 1]);
OEC = hsv2rgb([2/3 1 1]);


% All symbols: +o*.xsd^v<>ph - have them in that order to not get confused
New_symbols = {...
	%vascular
    'Rgs5', VECA, 'x'; ...
    'Flt1', VECA, 'p'; ...
    'Igf2', VECA, '+'; ...
    'Epas1', VECA, 'o'; ...
    'Bsg', VECV, 'd'; ...
    'Itm2a', VECV, '<'; ...
    'Ly6c1', VECV, '>'; ...
    'Car4', VECV, '^'; ...
    'Pltp', VECV, 'v'; ...
	'Vim', VECA, 'h'; ...

	
    'Crip1', VSMCA, 'p'; ...
    'Acta2', VSMCA, 'h'; ...
    'Myl9', VSMCA, 's'; ...
    'Tagln', VSMCA, 'o'; ...   
	
	'Ptgds', VLMC, '>'; ...
    'Igfbp2', VLMC, '^'; ...
    'Dcn', VLMC, '<'; ...
    'Nupr1', VLMC, 'v'; ...
	'Mgp', ABC, 'd'; ...
    'Rbp1', ABC, '.'; ...
    'Fn1', ABC, 's'; ...
	'Igfbp6', ABC, '*'; ...

	'Cldn5', PER, '+'; ...
    'Igfbp7', PER, 'o'; ...
    'Vtn', PER, '*'; ...
    'Sparc', PER, 'x'; ...
    'Ctla2a', PER, '.'; ...
    'Higd1b', PER, 's'; ...
    'Pglyrp1', PER, 'd'; ...
    'Myl12a', PER, '^'; ...
    'Gng11', PER, 'v'; ...
    'Serpinh1', PER, '>'; ...
    'Rhoc', PER, 'h'; ...
    'Col4a1', PER, '<'; ...

	% microglia
	'Hexb', MGL, 'p'; ...
    'Ctss', MGL, 'o'; ...
    'C1qa', MGL, 'h'; ...
    'P2ry12', MGL, '*'; ...
    'C1qb', MGL, 'x'; ...
    'Ccl4', MGL, 's'; ...
	% perivascular macrophages
	'Lyz2', PVM, '<'; ...
    'Pf4', PVM, '>'; ...
    'C1qc', PVM, '^'; ...
    'Fcer1g', PVM, 'v'; ...
    'Tyrobp', PVM, 'd'; ...

	% oligos
	'Cd9', COP, 'o'; ...
    'Frmd4a', COP, '+'; ...
    'Gpr17', COP, '*'; ...
    'Nfasc', NFOL, '.'; ...
    'Sirt2', NFOL, 'x'; ...
    'Marcks', NFOL, 's'; ...
    'Rras2', NFOL, 'd'; ...
    'Arpc1b', NFOL, '^'; ...
    'Tsc22d4', NFOL, 'v'; ...
    'Qk', NFOL, '<'; ...
    'Gng12', NFOL, '>'; ...
    'Gm15440', NFOL, 'p'; ...
    'Trf', MOL, '+'; ...
    'Mog', MOL, 'o'; ...
    'Apod', MOL, '*'; ...
    'Plp1', MOL, 'x'; ...
    'S100b', MOL, '.'; ...
    'Gltp', MOL, 's'; ...
    'Ppp1r14a', MOL, 'd'; ...
    'Car2', MOL, '^'; ...
    'Fa2h', MOL, 'v'; ...
    'Klk6', MOL, '>'; ...
    'Gjb1', MOL, '<'; ...
    'Sepp1', MOL, 'h'; ...
    'Gsn', MOL, 'p'; ...
	'Mal', MOL2, '+'; ...
	'Pllp', MOL2, 'o'; ...
    'Cdkn1c', MOL2, '*'; ...
    'Mag', MFOL, '.'; ...
    'Mobp', MFOL, 'x'; ...
    'Opalin', MFOL, 's'; ...
    'Tspan2', MFOL, 'd'; ...
    'Grb14', MFOL, 'v'; ...
    'Olig1', MFOL, '^'; ...
    'Tmem141', MFOL, '<'; ...
    'Pdlim2', MFOL, '>'; ...
    'Cd81', MFOL, 'p'; ...
	
	% astros

    'Slc1a3', ACBG, '+'; ...
    'Gfap', ACMB, '*'; ...
    'Aqp4', ACMB, 'o'; ...
    'Slc6a11', ACOB, 's'; ...
    'Slc1a2', ACTE, 'p'; ...
    'Aldoc', ACTE, 'h'; ...
	
	% choroid
    'Ttr', CHOR, '+'; ...
    'Enpp2', CHOR, 'o'; ...
    '1500015O10Rik', CHOR, '*'; ...
    'Fxyd1', CHOR, '.'; ...
    'Sostdc1', CHOR, 'x'; ...
	
	% ependymal
    '1110017D15Rik', EPEN, 's'; ...
    'Ccdc153', EPEN, 'd'; ...
    'Tmem212', EPEN, '<'; ...
    'Foxj1', EPEN, '>'; ...
    'Rarres2', EPEN, 'v'; ...
    'Dbi', EPSC, '^'; ...
    'Fos', EPSC, 'p'; ...
    'Ifitm3', EPSC, 'h'; ...
    'Rsph1', EPSC, '*'; ...
    'Jun', EPSC, '+'; ...
    'Mia', EPSC, 'o'; ...

	
	% neurons
	
    'Car8', CBPC, '*'; ...
	'Ppp1r1b', MSN, '+'; ...
    'Arpp21', MSN, 'o'; ...
    'Pde10a', MSN, '.'; ...
    'Efr3a', HBGLU, 'x'; ...
    'Slc17a6', HBGLU, 's'; ...
    'Snap25', HBGLU, 'd'; ...
    'Cabp7', HBGLU, '<'; ...
    'Gls', HBGLU, '>'; ...
    'Peg3', HBSER, '^'; ...
    'Cntn1', HBSER, 'h'; ...
    'Hs6st2', HBSER, 'p'; ...
    'En1', MBDOP, 'v'; ...
	
    'Irs4', HBINH, 'v'; ...
    'Lhx8', DECHO, 'p'; ...
    'Syt4', HBCHO, 'x'; ...
    'Vat1l', HBCHO, 's'; ...
	
	% neuroblasts
	'Pde6g', DETPH, '*'; ...
    'Sag', DETPH, '.'; ...
    'Gnb3', DETPH, 'x'; ...
    'Igfbpl1', DGNBL, 's'; ...
   'Doc2g', OBNBL, 'd'; ...
    'Cdhr1', OBNBL, '<'; ...
    'Nmb', OBNBL, '>'; ...
    'A230065H16Rik', SEPNBL, 'h'; ...
    'Hmgb2', SZNBL, 'p'; ...


    'Fabp7', OEC, 's'; ...
    'Rgcc', SATG, '+'; ...
    'Sfrp5', SATG, 'o'; ...
    'Fbln5', SATG, '*'; ...
    'Fbln2', SATG, '>'; ...
    'Ube2c', SATG, 'x'; ...


};

if nargout<1
    % delete any existing legend
    fc = get(gcf, 'Children');
    for i=1:length(fc)
        if strcmp(get(fc(i), 'UserData'), 'key')
            delete(fc(i));
        end
    end
    
    MainAxes = gca;
    
    n =  size(New_symbols,1);
    
    gc = get(MainAxes, 'children');
    MyChildren = [];
    for i=1:length(gc)
        if (strcmp(gc(i).Type, 'line') || strcmp(gc(i).Type, 'scatter')) ...
                && ~isempty(gc(i).DisplayName)
            MyChildren = [MyChildren; i];
        end
    end
    DisplayNames = {gc(MyChildren).DisplayName};
    % get first word of display name as gene
    GeneNames = cell(size(DisplayNames));
    for i=1:length(DisplayNames)
        GeneNames{i} = strtok(DisplayNames{i});
    end
    
    clear h s;
    j=1;
    Present = [];
    for i=1:n
        MyGeneName = New_symbols{i,1};
        l = find(strcmp(MyGeneName, GeneNames));
        if ~isempty(l)
            h(j) = gc(MyChildren(l));
            if strcmp(h(j).Type, 'line')
                set(h(j), 'Color', New_symbols{i,2});
            elseif strcmp(h(j).Type, 'scatter')
                set(h(j), 'CData', New_symbols{i,2});
            end
            set(h(j), 'Marker', New_symbols{i,3});
            
            if MarkerSize>0
                if strcmp(gc(l).Type, 'line')
                    set(h(j), 'MarkerSize', MarkerSize);
                elseif strcmp(gc(l).Type, 'scatter')
                    set(h(j), 'SizeData', MarkerSize);
                end
            end
            Present(j) = i;
            j=j+1;
        end
    end
    
    other_h = setdiff(gc(MyChildren), h);
    other_symbols = {other_h.DisplayName};
    
    all_h = [h(:); other_h(:)];
    all_sym = {New_symbols{Present,1}, other_symbols{:}};
    
    % lh = legend([h(:); other_h(:)], ...
    %     {New_symbols{s,1}, other_symbols{:}}, ...
    %     'color', 'k', 'textcolor', 'w', 'fontsize', FontSize);
    % set(lh, 'color', 'k');
    %
    % return;
    
    if MultiCol==0
        lh = legend(all_h, all_sym, 'color', 'k', 'textcolor', 'w', 'fontsize', FontSize);
        set(lh, 'color', 'k');
    else
        ah = axes('Position', [.925 .13 .05 .8]);
        set(ah, 'color', 'k'); cla; hold on; box off
        set(ah, 'UserData', 'key');
        for j=1:length(Present)
            i = Present(j);
            plot(ceil(j/MultiCol)+.1, mod(j-1,MultiCol), New_symbols{i,3}, 'Color', New_symbols{i,2});
            text(ceil(j/MultiCol)+.3, mod(j-1,MultiCol), New_symbols{i,1}, 'color', 'w', 'fontsize', FontSize);
        end
        if ~isempty(setdiff(GeneNames,New_symbols(:,1)))
            j=j+1;
            plot(ceil(j/MultiCol)+.1, mod(j-1,MultiCol), '.', 'Color', hsv2rgb([0,0,0.5]));
            text(ceil(j/MultiCol)+.3, mod(j-1,MultiCol), 'Neuron', 'color', 'w', 'fontsize', FontSize);
        end
        ylim([-1 MultiCol]);
        set(ah, 'xtick', []);
        set(ah, 'ytick', []);
        set(ah, 'ydir', 'reverse');
    end
    %     for c=1:nCols
    %         rr=((c-1)*50 + 1):min(c*50, length(all_h));
    %         if c==1
    %             ah(c) = gca;
    %             lh(c) = legend(all_h(rr), all_sym(rr), 'color', 'k', 'textcolor', 'w', 'fontsize', FontSize, 'location', 'east');
    %             set(lh(c), 'color', 'k');
    %             pos(c,:) = get(lh(c), 'position');
    %         else
    %             ah(c) = axes('position',get(gca,'position'), 'visible','off');
    %             lh(c) = legend(ah(c), all_h(rr), all_sym(rr), 'color', 'k', 'textcolor', 'w', 'fontsize', FontSize, 'location', 'east');
    %             set(lh(c), 'position', pos(c-1,:) + [1.1 0 0 0]*pos(c-1,3));
    %             uistack(lh(c), 'top');
    %         end
    %     end
    %     axes(ah(1));
    
    %    error('multicolumn not done yet!');
    % end
    %     for i=1:nCols
    %         first = 1+(i-1)*nCols;
    %         last = min(i*nCols,length(all_h));
    %         lh = legend(all_h(first:last), ...
    %             all_sym{first:last});%, ...
    %             %'color', 'k', 'textcolor', 'w', 'fontsize', FontSize);
    %         set(lh, 'color', 'k');
    %     end
    set(gcf, 'color', 'k');
    set(gcf, 'InvertHardcopy', 'off');
    
    axes(MainAxes)
    uistack(ah, 'top');
end

end