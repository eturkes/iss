function change_gene_symbols_DAM_MG150(MarkerSize, FontSize, MultiCol)
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
c4a =    hsv2rgb([.1 1 .6]);
c4b =    hsv2rgb([.1 .4 .6]);
c5a =    hsv2rgb([.55 1 1]);
c5b =    hsv2rgb([2/3 1 1]);
c5c =    hsv2rgb([.7 .8 1]);
c5d =    hsv2rgb([.75 1 .7]);
cNAa =   hsv2rgb([0 0 1]);

% All symbols: +o*.xsd^v<>ph - have them in that order to not get confused
New_symbols = {...

    'Cd52',   c1a, '+'; ...

    'Cxcl16',     c1a, 'o'; ...

    'Cd22',    c1a, '*'; ...

    'Atp6v0d2',      c2a, '+'; ...

    'Ctsd',     c3a, '+'; ...
    'Prune2',     c3a, 'o'; ...
    'Pld3',     c3a, '*'; ...

    'Ctsb',   c3a, '.'; ...
    'Cst7',    c3a, 'x'; ...
    'Ctsz',  c3a, 's'; ...
    'Ctsl',      c3a, 'd'; ...
    'Tpi1',     c3a, '^'; ...


    'Capg',    c4a, '+';...
    'Baiap2l2',    c4a, 'o'; ...

    'Cd63',    c5a, '+'; ...

    'Gpnmb',    c5a, 'o'; ...
    'Cd9',   c5a, '*'; ...
    'Clec7a',   c5a, '.'; ...
    'Fcgr4',     c5a, 'x'; ...

    'H2-DMa',   c5a, 's'; ...
    'Slamf9',    c5a, 'd'; ...

    'Lgi2',    cNAa, '+'; ...

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
