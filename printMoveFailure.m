function printMoveFailure(m)
% 

failstr='unknown move';
    switch(m)
        case 1
            failstr='merging does no good!';
        case 2
            failstr='split does no good!';
        case 3
            failstr='growing does no good!';
        case 4
            failstr='shrinking does no good!';
        case 5
            failstr='purge does no good!';
        case 6
            failstr='adding does no good!';
    end
    printMessage(2,'%s\n',failstr);
end