%%

fpath='/home/aanton/visinf/theses/ongoing/andriyenko-phd/Chapters/Cont/numbers/';
fid=fopen([fpath 'start-gt.tex'],'w');


allscen=[23 25 71];
allm={'S2L1 (GT)','     (EKF)','S2L2 (GT)','     (EKF)','S1L1 (GT)','     (EKF)'};
scncnt=0;
for scen=allscen
    scncnt=scncnt+1;
    col=getColorFromID(scncnt);
    load(sprintf('%ss%d-from-gt.mat',exppath,scen));    

    sallens=sum(LOG_allens,2);
    inds=find(sallens);
    ie=inds(end);
    gte=sallens(end);
    nensvalues=length(inds);
    outstr=sprintf('%s & %.1f & %.1f & %.1f & %.1f & %d & %d & %d\\\\\\\\ \n', ...
        char(allm(scncnt*2-1)),sallens(1),sallens(ie),allmets3d(1,12),allmets3d(ie,12), ...
        allmets3d(ie,8),allmets3d(ie,9),allmets3d(ie,10));

    fprintf(fid,sprintf('%s',outstr));
    
    load(sprintf('%ss%d-from-ekf.mat',exppath,scen));    

    sallens=sum(LOG_allens,2);inds=find(sallens);    ie=inds(end);

    ekfe=sallens(end);
    nensvalues=length(inds);
    outstr=sprintf('%s & %.1f & %.1f & %.1f & %.1f & %d & %d & %d\\\\\\\\ \n', ...
        char(allm(scncnt*2)),sallens(1),sallens(ie),allmets3d(1,12),allmets3d(ie,12), ...
        allmets3d(ie,8),allmets3d(ie,9),allmets3d(ie,10));
    fprintf(fid,sprintf('%s',outstr));
    fprintf(fid,sprintf('\\\\hline \n'));

end
fclose(fid);