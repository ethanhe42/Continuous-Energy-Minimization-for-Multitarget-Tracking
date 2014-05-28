global opt
allscen=700:720;

settingsDir='KITTI_PED';
jobid=73;

confdir=sprintf('config/%s',settingsDir);

resdir=sprintf('results/%s',settingsDir);
if ~exist(resdir,'dir'), mkdir(resdir); end
conffile=fullfile(confdir,'inis',sprintf('%04d.ini',jobid));


maxrandruns=1;
maxexpruns=1;
allmets2d=zeros(max(allscen),14,maxrandruns);
allmets3d=zeros(max(allscen),14,maxrandruns);
allens=zeros(max(allscen),maxrandruns);
allInfo=[];

jobid=73;
tic;
for scenario=allscen
    [metrics2d, metrics3d, allens, stateInfo]=cemTracker(scenario,conffile);
    
    resultsfile=fullfile(resdir,sprintf('res_%03d-scen%04d.mat',jobid,scenario));
    save(resultsfile,'metrics2d','metrics3d','allens','stateInfo');
    randrun=opt.startsol;
    allmets2d(scenario,:,randrun)=metrics2d;
    allmets3d(scenario,:,randrun)=metrics3d;
    allens(scenario,randrun)=sum(allens);
    allInfo(scenario,randrun).stateInfo=stateInfo;    
    
end
toc

mets=[];
resultsfile=fullfile(resdir,sprintf('res_%03d.mat',jobid));
save(resultsfile,'allmets*','allens','allInfo','mets*');

%% delete intermediate
for scenario=allscen
    resultsfile=fullfile(resdir,sprintf('res_%03d-scen%04d.mat',jobid,scenario));
    if exist(resultsfile,'file'),        delete(resultsfile);    end
end



%%
resultsfile=fullfile(resdir,sprintf('res_%03d.mat',jobid));

thiswd=pwd;
load(resultsfile);

pathToKittiDevkit='/home/amilan/storage/databases/KITTI/tracking/devkit_tracking/python';


% allscen=750:778;
for scenario=allscen
    stateInfo=allInfo(scenario,6).stateInfo;
    cd(thiswd)
    tracklets=convertToKITTI(stateInfo,'Pedestrian');
    cd(pathToKittiDevkit)
    
    writeLabels(tracklets,'results/a/data',scenario-700);
end

