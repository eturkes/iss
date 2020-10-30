function change_gene_symbols_individual_MG150(MarkerSize, FontSize, MultiCol)
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

if nargin<2 || isempty(FontSize)
    FontSize = 5;
end

if nargin<3
    MultiCol = 49;
end


% colors
c1a =    hsv2rgb([ 1 1 1]);
c1b =    hsv2rgb([1 .6 1]);
c1c =    hsv2rgb([1 .3 .9]);
c2a =    hsv2rgb([ .13 1 1]);
c3a =    hsv2rgb([1/3 1 1]);
c3b =    hsv2rgb([.27 1 .7]);
c4a =     hsv2rgb([.1 1 .6]);
c4b =     hsv2rgb([.1 .4 .6]);
c5a =    hsv2rgb([.55 1 1]);
c5b =    hsv2rgb([2/3 1 1]);
c5c =    hsv2rgb([.7 .8 1]);
c5d =    hsv2rgb([.75 1 .7]);
cNAa =   hsv2rgb([0 0 1]);

% All symbols: +o*.xsd^v<>ph - have them in that order to not get confused
New_symbols = {...
    
    'Tmem178',     c1a, '.'; ...
    'Fcrl1',    c1a, '.'; ...
    'C2',     c1a, '.'; ...
    'C1qb',    c1a, '.'; ...
    'C1qc',     c1a, '.'; ...
    'Dcstamp',   c1a, '.'; ...
    'Cd300lf',   c1a, '.'; ...
    'C1qa',   c1a, '.'; ...
    'Clec4a2',   c1a, '.'; ...
    'Rsad2',   c1a, '.'; ...
    'Upk1b',      c1a, '.'; ...
    'Cd52',   c1a, '.'; ...
    'Gpr160',    c1a, '.'; ...

    'Cd79b',    c1b, '.'; ...
    'Cd48',      c1b, '.'; ...
    'Tnfrsf17',     c1b, '.'; ...
    'Ccrl2',     c1b, '.'; ...
    'Cxcl13',     c1b, '.'; ...
    'Cxcl16',     c1b, '.'; ...
    'Gbp2',     c1b, '.'; ...
    'Bst2',     c1b, '.'; ...
    'Ifitm6',     c1b, '.'; ...
    'Cd300lb',       c1b, '.'; ...
    'Ifit3',    c1b, '.'; ...
    'Siglece',      c1b, '.'; ...
    'Colec12',      c1b, '.'; ...

    'Cd22',    c1c, '.'; ...
    'Irf7',    c1c, '.'; ...
    'Xaf1',  c1c, '.'; ...
    'Oasl',    c1c, '.'; ...
    'Psmb8',     c1c, '.'; ...
    'Usp18',   c1c, '.'; ...
    

    'Atp6v0d1',     c2a, '.'; ...
    'Atp6v1d',      c2a, '.'; ...
    'Rab39',   c2a, '.'; ...
    'Slc1a5',     c2a, '.'; ...
    'Ykt6',     c2a, '.'; ...
    'Atp6v0d2',      c2a, '.'; ...
    'Slc16a3', c2a, '.'; ...
    'Kcnk13',   c2a, '.'; ...
    'Tmem37', c2a, '.'; ...
    'Atp6ap1',   c2a, '.'; ...
    'Atp6v0c',       c2a, '.'; ...
    'Atp6ap2',  c2a, '.'; ...


    'Apoc2', c3a, '.'; ...
    'Ctsd',     c3a, '.'; ...
    'Pld4',     c3a, '.'; ...
    'Prune2',     c3a, '.'; ...
    'Nceh1',     c3a, '.'; ...
    'Gusb',     c3a, '.'; ...
    'Pld3',     c3a, '.'; ...
    'Cmas',    c3a, '.'; ...
    'Mymk',   c3a, '.'; ...
    'Dkk2',   c3a, '.'; ...
    'Ecscr',      c3a, '.'; ...
    'Fam20c',   c3a, '.'; ...
    'Tmem119',   c3a, '.'; ...

    'Ctsb',   c3b, '.'; ...
    'Cst7',    c3b, '.'; ...
    'Ctss',     c3b, '.'; ...
    'Acaca',    c3b, '.'; ...
    'Dna2',    c3b, '.'; ...
    'Ctsz',  c3b, '.'; ...
    'Lgmn',    c3b, '.'; ...
    'Ctsl',      c3b, '.'; ...
    'Rps6ka1',     c3b, '.'; ...
    'Nans',     c3b, '.'; ...
    'Mgat1',     c3b, '.'; ...
    'Mbnl1',    c3b, '.';...
    'Tpi1',     c3b, '.'; ...


    'Capg',    c4a, '.';...
    'Cotl1',    c4a, '.'; ...
    'Abi1',   c4a, '.'; ...
    'Actr3b',     c4a, '.'; ...
    'Baiap2l2',    c4a, '.'; ...
    'Arpc2',    c4a, '.'; ...
    'Fscn1',     c4a, '.'; ...
    'Cfl1',      c4a, '.'; ...
    'Pfn1',    c4a, '.'; ...
    'Cyfip1',   c4a, '.'; ...
    'Nckap1l',   c4a, '.'; ...
    'Bin1',    c4a, '.'; ...
    'Bin2',     c4a, '.'; ...

    'Plek',      c4b, '.'; ...


    'Abi3',     c5a, '.'; ...
    'Cd63',    c5a, '.'; ...
    'Abca7',    c5a, '.'; ...
    'Sparc',     c5a, '.';...
    'Gsn',    c5a, '.'; ...
    'Rraga',    c5a, '.'; ...
    'Egln3',    c5a, '.'; ...
    'Hcar2',    c5a, '.'; ...
    'Alox5ap',      c5a, '.'; ...
    'Pilra',   c5a, '.'; ...
    'Laptm5',     c5a, '.'; ...
    'Fcgr2b',  c5a, '.'; ...
    'Trem2',    c5a, '.'; ...

    'Ltf',     c5b, '.'; ...
    'Tnfrsf1b',    c5b, '.'; ...
    'Tgfbr1',   c5b, '.'; ...
    'Csf1r',    c5b, '.'; ...
    'Ptk2b',    c5b, '.'; ...
    'Nr1h3',   c5b, '.'; ...
    'Plcg2',   c5b, '.'; ...
    'Nr1d1',      c5b, '.'; ...
    'Nr1h4',   c5b, '.'; ...
    'Grn',     c5b, '.'; ...
    'Gnai2',  c5b, '.'; ...
    'P2ry12',    c5b, '.'; ...
    'Ccr5',     c5b, '.'; ...

    'Cx3cr1',    c5c, '.'; ...
    'Gpnmb',    c5c, '.'; ...
    'Il6ra',   c5c, '.'; ...
    'Coro1a',    c5c, '.'; ...
    'Itgb2',    c5c, '.'; ...
    'Cd47',   c5c, '.'; ...
    'Cd33',      c5c, '.'; ...
    'Cd9',   c5c, '.'; ...
    'Hpse',     c5c, '.'; ...
    'Itgax',  c5c, '.'; ...
    'Clec7a',   c5c, '.'; ...
    'Fcgr4',     c5c, '.'; ...
    'H2-Ob',     c5c, '.'; ...

    'H2-DMb2',  c5d, '.'; ...
    'H2-DMa',   c5d, '.'; ...
    'Ncf4', c5d, '.'; ...
    'Itgam',    c5d, '.'; ...
    'Clec4a3',   c5d, '.'; ...
    'Slamf9',    c5d, '.'; ...

    'Lgi2',    cNAa, '.'; ...
    'Nr1c1',    cNAa, '.'; ...
    'Nr1c2',    cNAa, '.'; ...
    'Nr1c3',   cNAa, '.'; ...
    'Nr2b1',    cNAa, '.'; ...
    'Nr2b2',    cNAa, '.'; ...
    'Nr2b3',   cNAa, '.'; ...
    'Tm2d3',   cNAa, '.'; ...
    'Tmem163',    cNAa, '.'; ...
    'Olfml3',   cNAa, '.'; ...

    };

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
    ah = axes('Position', [.01 .13 .06 .8]);
    set(ah, 'color', 'k'); cla; hold on; box off
    set(ah, 'UserData', 'key');
    for j=1:length(Present)
        i = Present(j);
        plot(ceil(j/MultiCol)+.1, mod(j-1,MultiCol), New_symbols{i,3}, 'Color', New_symbols{i,2}, 'MarkerSize', MarkerSize);
        text(ceil(j/MultiCol)+.3, mod(j-1,MultiCol), New_symbols{i,1}, 'color', 'w', 'fontsize', FontSize);
    end
    ylim([-1 MultiCol]);
    set(ah, 'Visible', 'off')
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
