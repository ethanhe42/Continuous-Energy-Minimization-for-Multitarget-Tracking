function [EdetValue EdynValue EexcValue EappValue EperValue EregValue,EoriValue]= ...
    printEnergies(opt, i, stateVec, stateInfo)        

global cemStartTime globiter 

if opt.verbosity>=3
    [fx, ~, EdetValue EdynValue EexcValue EappValue EperValue EregValue EoriValue] = E(stateVec, stateInfo);
    
    % weight energy values
    EdetValue=EdetValue*opt.wtEdet;
    EdynValue=EdynValue*opt.wtEdyn;
    EexcValue=EexcValue*opt.wtEexc;
    EappValue=EappValue*opt.wtEapp;
    EperValue=EperValue*opt.wtEper;
    EregValue=EregValue*opt.wtEreg;
    EoriValue=EoriValue*opt.wtEori;
    
    printMessage(3,'%4i|%3i|%6.1f|%8.1f|%8.1f|%6.1f|%6.1f|%6.1f|%6.1f|%6.1f|%6.1f|||', ...
        globiter,i, toc(cemStartTime)/60,fx,EdetValue,EdynValue,EexcValue,EappValue,EperValue,EregValue,EoriValue);
%     LOG_allens(globiter,:)=[EdetValue EdynValue EexcValue EappValue EperValue EregValue,EoriValue];
end
        
end