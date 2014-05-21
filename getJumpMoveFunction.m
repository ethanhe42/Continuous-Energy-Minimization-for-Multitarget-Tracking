function moveFct=getJumpMoveFunction(m)
% returns the name of the function
% corresponding to jump move m
% 


switch(m)
    case 1
        moveFct='makeMergeMove';
    case 2
        moveFct='makeSplitMove';
    case 3
        moveFct='makeGrowMove';
    case 4
        moveFct='makeShrinkMove';
    case 5
        moveFct='makeRemoveMove';
    case 6
        moveFct='makeAddMove';
end
end