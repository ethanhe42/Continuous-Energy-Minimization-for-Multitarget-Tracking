function runGridSearchCluster(jobname,jobid)

settingsDir=strsplit(jobname,'-');
settingsDir=char(settingsDir{1});
confdir=sprintf('config/%s',settingsDir);

jobid=str2double(jobid);

resdir=sprintf('results/%s',settingsDir);
if ~exist(resdir,'dir'), mkdir(resdir); end

conffile=fullfile(confdir,'inis',sprintf('%04d.ini',jobid));


allscen=dlmread(fullfile(confdir,'doscens.txt'));

gtI.Xi=[];gtI.Yi=[];gtI.W=[];gtI.H=[];
stI.Xi=[];stI.Yi=[];stI.W=[];stI.H=[];
% to combine all sequences into one file, evaluate and delete individual ones
maxrandruns=1;
maxexpruns=1;
allmets2d=zeros(max(allscen),14,maxrandruns);
allmets3d=zeros(max(allscen),14,maxrandruns);
allens=zeros(max(allscen),maxrandruns);
allInfo=[];



global gtInfo opt;
for scenario=allscen
    [metrics2d, metrics3d, allens, stateInfo]=cemTracker(scenario,conffile);
    
    resultsfile=fullfile(resdir,sprintf('res_%03d-scen%04d.mat',jobid,scenario));
    save(resultsfile,'metrics2d','metrics3d','allens','stateInfo');
    
    % create one state Info for all sequences
    [F,N]=size(stateInfo.Xi);    [Fo,No]=size(stI.Xi);
    newf=Fo+1:Fo+F;    newi=No+1:No+N;
    stI.Xi(newf,newi)=stateInfo.Xi;    stI.Yi(newf,newi)=stateInfo.Yi;
    stI.W(newf,newi)=stateInfo.W;    stI.H(newf,newi)=stateInfo.H;
        
    % create one groundtruth for all sequences
    [F,N]=size(gtInfo.Xi);    [Fo,No]=size(gtI.Xi);
    newf=Fo+1:Fo+F;    newi=No+1:No+N;
    gtI.Xi(newf,newi)=gtInfo.Xi;    gtI.Yi(newf,newi)=gtInfo.Yi;
    gtI.W(newf,newi)=gtInfo.W;    gtI.H(newf,newi)=gtInfo.H;
    
    randrun=opt.startsol;
    allmets2d(scenario,:,randrun)=metrics2d;
    allmets3d(scenario,:,randrun)=metrics3d;
    allens(scenario,randrun)=sum(allens);
    allInfo(scenario,randrun).stateInfo=stateInfo;    
end

gtI.frameNums=1:size(gtI.Xi,1);stI.frameNums=1:size(stI.Xi,1);
gtI.X=gtI.Xi;gtI.Y=gtI.Yi; stI.X=stI.Xi;stI.Y=stI.Yi;
stI.Xgp=stI.X;stI.Ygp=stI.Y;gtI.Xgp=gtI.X;gtI.Ygp=gtI.Y;

    
mets2d=CLEAR_MOT(gtI,stI);
mets3d=zeros(size(mets2d));
if opt.track3d, mets3d=CLEAR_MOT(gtI,stI,struct('eval3d','1')); end

resultsfile=fullfile(resdir,sprintf('res_%03d.mat',jobid));
save(resultsfile,'allmets*','allens','allInfo','mets*');

%% delete intermediate
for scenario=allscen
    resultsfile=fullfile(resdir,sprintf('res_%03d-scen%04d.mat',jobid,scenario));
    if exist(resultsfile,'file'),        delete(resultsfile);    end
end

end