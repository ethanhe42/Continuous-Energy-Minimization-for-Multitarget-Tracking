function compileAndQuit(fileToCompile)
% Try to compile (mcc) and then quit matlab

fprintf('compiling %s...\n',fileToCompile);
[pt fl ex]=fileparts(fileToCompile);

compstr=sprintf('mcc -a /gris/gris-f/home/aandriye/software/lightspeed/xrepmat.m -a E.m -a makeAddMove.m -a makeGrowMove.m -a makeMergeMove.m -a makeRemoveMove.m -a makeShrinkMove.m -a makeSplitMove.m -I ~/software/lightspeed -I ./utils -I ./utils/splinefit -I ./utils/camera -I ./utils/mex/bin -I ./mex/bin -R -singleCompThread -R -nodisplay -C %s%s -m %s;',fl,ex,fl);
compstr

while true
    try
	fprintf('trying...\n');
	eval(compstr);
	quit;
    catch err
        fprintf('not this time\n');
        pause(15);
    end
end
