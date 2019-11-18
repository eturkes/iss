function [NewShift,Outliers] = AmendShifts(o,shifts,score,section)
%For a particular set of shifts, this considers the score for each shift.
%Then if a particular score is too low and the shift is an outlier, the
%shift is changed to the average of all the other shifts.
%For regiter.m, a set of shifts corresponds to a particular direction
%(South or East).
%For find_spots.m, a set of shifts corresponds to a particular round.
%In z direction, shifts are often very similar hence outlier should be
%found if difference to median more than an absolute value. This absolute
%value is taken to be o.OutlierThresh.
if strcmpi(section, 'FindSpots')
    AbsoluteMinScore = o.FindSpotsAbsoluteMinScore;
elseif strcmpi(section, 'Register')
    AbsoluteMinScore = o.RegAbsoluteMinScore;
end

Outliers = zeros(size(shifts));
AnomalousScores = score<o.OutlierMinScore;
if max(AnomalousScores)>0
    warning('Looking at anomalous shifts');
    AnomalousShiftArray = cat(3,isoutlier(shifts(:,1),'median','ThresholdFactor',o.OutlierThresh)...
        ,isoutlier(shifts(:,2),'median','ThresholdFactor',o.OutlierThresh)...
        ,abs(shifts(:,3)-median(shifts(:,3)))>=o.OutlierThresh);
    AnomalousShift = max(AnomalousShiftArray,[],3);
    AwfulScore = score < AbsoluteMinScore;
    for i=1:size(score,1)
        if min([AnomalousScores(i),AnomalousShift(i)+AwfulScore(i)])>0
            Outliers(i,:) = shifts(i,:);
            shifts(i,:) = round(mean(shifts(AnomalousScores==0,:)));    
            warning('shift(%d) changed',i);
        end
    end
end

NewShift = shifts;
end