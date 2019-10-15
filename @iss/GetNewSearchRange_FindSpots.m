function o = GetNewSearchRange_FindSpots(o,t,r)

%For a particular round r, this considers all shifts found so far up to
%tile t and adjusts the search range based on them. The Y and X search
%ranges are always changed when this is called, but the Z range is only
%changed if the previous range was more severely wrong.

%Find extremes of shifts found so far
MinY = min(o.D0(1:t,1,r));
MaxY = max(o.D0(1:t,1,r));
MinX = min(o.D0(1:t,2,r));
MaxX = max(o.D0(1:t,2,r));
MinZ = min(o.D0(1:t,3,r))*o.XYpixelsize/o.Zpixelsize;
MaxZ = max(o.D0(1:t,3,r))*o.XYpixelsize/o.Zpixelsize;

%Adjust search ranges
o.FindSpotsSearch{r}.Y = MinY-2*o.FindSpotsStep(1):o.FindSpotsStep(1):MaxY+2*o.FindSpotsStep(1);
o.FindSpotsSearch{r}.X = MinX-2*o.FindSpotsStep(2):o.FindSpotsStep(2):MaxX+2*o.FindSpotsStep(2);

if MaxZ - MinZ > 2*o.FindSpotsStep(3) && o.FindSpotsStep(3)<3
    o.FindSpotsStep(3) = o.FindSpotsStep(3)+1;
    o.FindSpotsSearch{r}.Z = MinZ:o.FindSpotsStep(3):MaxZ;
elseif MinZ <= min(o.FindSpotsSearch{r}.Z)-o.FindSpotsStep(3) || ...
        MaxZ >= max(o.FindSpotsSearch{r}.Z)+o.FindSpotsStep(3)
    o.FindSpotsSearch{r}.Z = MinZ:o.FindSpotsStep(3):MaxZ;
end