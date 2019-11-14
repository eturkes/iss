function [BestShift,BestScore,ChangedSearchRange] = get_initial_shift2(o, y, x0,Search,section)
% Finds the initial shifts to give to PCR algorithm between rounds for each
% tile. Does this by finding the colour channel with the most spots for
% each round. Then uses ImRegFft2 to find the shift between this base binary
% image and the anchor channel binary image
%
% inputs:
% y is a cell containing the centered YXZ location of all spots in round r 
% , colour channel c for all tiles in units of XY pixels in find_spots. In
% register, y is equivalent to x0 but for another tile.
%
% x0 is a cell containing the non centered YXZ location of spots in the 
% anchor channel for the current tile. Z units are Z pixels.
%
% Search.Y,Search.X and Search.Z are the ranges in XY and Z pixel size 
% respectively of shifts to search.
%
% section specifies which part of the pipeline we are on: Register or
% FindSpots
%
% Now searches a wide range with gaps then does a finer search with
% no gaps about the maxima found.

%% 
if strcmpi(section, 'FindSpots')
    %centre anchor channel spots
    x = (x0 - o.CentreCorrection).*[1,1,o.Zpixelsize/o.XYpixelsize];
    MinScore = o.FindSpotsMinScore;
    Step = o.FindSpotsStep;
    WidenSearch = o.FindSpotsWidenSearch;
    RefineSearch = o.FindSpotsRefinedSearch;
    RefineStep = o.FindSpotsRefinedStep;
elseif strcmpi(section, 'Register')
    x = x0.*[1,1,o.Zpixelsize/o.XYpixelsize];
    y = y.*[1,1,o.Zpixelsize/o.XYpixelsize];
    MinScore = o.RegMinScore;
    Step = o.RegStep;
    WidenSearch = o.RegWidenSearch;
    RefineSearch = o.RegRefinedSearch;
    RefineStep = o.RegisterRefinedStep;
else
    warning('Have not specified which part of pipeline is being run')
end

% make kd tree - default options!
k = KDTreeSearcher(y);

%do initial large search
[shifts,Score] = get_score(Search,x,k,o);
BestScore = max(Score);
BestShift = shifts(Score == BestScore,:);

ChangedSearchRange = 0;

if strcmpi('auto',MinScore)
    MinScore = median(Score)+5*iqr(Score);
end

%if maxima below MinScore, widen search
if BestScore < MinScore
    warning('BestScore %d is less than the minimum score required, %d. Now widening search range.',BestScore,MinScore);
    Search.Y = min(Search.Y)-WidenSearch(1)+1:Step(1):max(Search.Y)+WidenSearch(1)+1;
    Search.X = min(Search.X)-WidenSearch(2)+1:Step(2):max(Search.X)+WidenSearch(2)+1;
    Search.Z = min(Search.Z)-WidenSearch(3)+1:Step(3):max(Search.Z)+WidenSearch(3)+1;
    [NewShifts,NewScore] = get_score(Search,x,k,o);
    shifts = [shifts;NewShifts];
    Score = [Score;NewScore];
    BestScore = max(Score);
    BestShift = shifts(Score == BestScore,:);
    ChangedSearchRange = 1;
end

if BestScore < MinScore
    warning('BestScore %d is still less than the minimum score required, %d, even after widening search range. SHIFT PROBABLY WRONG',BestScore,MinScore);
    ChangedSearchRange = 0;
end

%Refine Search around maxima
%Refined Search with smaller step
Search.Y = BestShift(1)-RefineSearch(1):RefineStep(1):BestShift(1)+RefineSearch(1);
Search.X = BestShift(2)-RefineSearch(2):RefineStep(2):BestShift(2)+RefineSearch(2);
ZShift = BestShift(3)*o.XYpixelsize/o.Zpixelsize;
Search.Z = ZShift-RefineSearch(3):RefineStep(3):ZShift+RefineSearch(3);
[NewShifts,NewScore] = get_score(Search,x,k,o);
shifts = [shifts;NewShifts];
Score = [Score;NewScore];
BestScore = max(Score);
BestShift = shifts(Score == BestScore,:);

if size(BestShift,1)>1
    BestShift = BestShift(1,:);
end

%Final search with no step in single Z plane
Search.Y = BestShift(1)-RefineStep(1):BestShift(1)+RefineStep(1);
Search.X = BestShift(2)-RefineStep(2):BestShift(2)+RefineStep(2);
ZShift = BestShift(3)*o.XYpixelsize/o.Zpixelsize;
Search.Z = ZShift-RefineStep(3):ZShift+RefineStep(3);
[NewShifts,NewScore] = get_score(Search,x,k,o);
shifts = [shifts;NewShifts];
Score = [Score;NewScore];
BestScore = max(Score);
BestShift = shifts(Score == BestScore,:);


if size(BestShift,1)>1
    BestShift = BestShift(1,:);
end


if o.Graphics == 2
    plotShiftSearch(int64(shifts.*[1,1,o.XYpixelsize/o.Zpixelsize]),Score,...
        strcat('Scores for different shifts'));
end


fprintf('\n');
return
end

%% 
function [shifts, Score] = get_score(Search,x,tree,o)
Z = Search.Z.*o.Zpixelsize/o.XYpixelsize;
%Find all permutations of shifts to try
[A,B,C] = meshgrid(Search.Y,Search.X,Z);
c=cat(3,A,B,C);
shifts=reshape(c,[],3);

Score = zeros(size(shifts,1),1);
for i = 1:size(shifts,1)
    xShifted = x+shifts(i,:);
    [~,Dist] = tree.knnsearch(xShifted);
    Score(i) = sum(exp(-Dist.^2/(2*o.ShiftScoreThresh^2)));
end
return
end

