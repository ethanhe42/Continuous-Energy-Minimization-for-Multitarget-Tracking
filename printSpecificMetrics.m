%%
metricsInfo.names.short = {'Rcll','Prcn','FAR', ...
    'GT','MT','PT','ML','FP', 'FN', 'IDs', 'FM','MOTA','MOTP', 'MOTAL'};
metricsInfo.names.long= {'Rcll','Prcn','FAR', ...
    'GT','MT','PT','ML','FP', 'FN', 'IDs', 'FM','MOTA','MOTP', 'MOTAL'};
metricsInfo.widths.long = [6 9 16 9 14 17 11 15 15 11 14 5 5 8];
metricsInfo.widths.short = [5 5 5 3 3 3 3 4 4 3 3 5 5 5];
metricsInfo.format.long = {'.1f','.1f','.2f', ...
    'i','i','i','i','i','i','i','i','.1f','.1f','.1f'};

metricsInfo.format.short=metricsInfo.format.long;

useexp=1;
load([resdir sprintf('res_%03d.mat',useexp-1)]);

allscen=[103 113]; % joint
allscen=[104 114]; % tracklets new

fprintf('\n');
scen=allscen(1);
printMetrics(mets3d(scen,:,6),metricsInfo,1,[1 2 12 13 5 6 7 8 9 11 10]);fprintf('\n');

global gtInfo sceneInfo opt scenario
scenario=scen;
opt=getConOptions;
sceneInfo=getSceneInfo(scen);  
evopt.eval3d=1;
stateInfo=infos(scen,6).stateInfo;


% t1=350; t2=370; idtr=3;
% stateInfo.X(t1:t2,idtr)=0;stateInfo.Y(t1:t2,idtr)=0;
% [stateInfo.X stateInfo.Y stateInfo]=cleanState(stateInfo.X, stateInfo.Y,stateInfo);
% % pause
% [metr,metrInfo,addInf]=CLEAR_MOT(gtInfo,stateInfo,evopt);
% printMetrics(metr,metrInfo,0,[1 2 12 13 5 6 7 8 9 11 10]);    fprintf('\n');
% 
% t1=362; t2=369; idtr=102;
% stateInfo.X(t1:t2,idtr)=0;stateInfo.Y(t1:t2,idtr)=0;
% [stateInfo.X stateInfo.Y stateInfo]=cleanState(stateInfo.X, stateInfo.Y,stateInfo);
% % pause
% [metr,metrInfo,addInf]=CLEAR_MOT(gtInfo,stateInfo,evopt);
% printMetrics(metr,metrInfo,0,[1 2 12 13 5 6 7 8 9 11 10]);    fprintf('\n');


scen=allscen(2);
printMetrics(mets3d(scen,:,6),metricsInfo,0,[1 2 12 13 5 6 7 8 9 11 10]);fprintf('\n');