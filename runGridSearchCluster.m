function runGridSearchCluster(jobname,jobid)

settingsDir=strsplit(jobname,'-');
settingsDir=char(settingsDir{1});
confdir=sprintf('config/%s',settingsDir);

jobid=str2double(jobid);

resdir=sprintf('results/%s',settingsDir);
if ~exist(resdir,'dir'), mkdir(resdir); end

conffile=fullfile(confdir,sprintf('%04d.ini',jobid));


allscen=dlmread(fullfile(confdir,'doscens.txt'));
for scenario=allscen
    [metrics2d, metrics3d, allens, stateInfo]=cemTracker(scenario,conffile);
    
    resultsfile=fullfile(resdir,sprintf('res_%03d-scen%04d.mat',jobid,scenario));
    save(resultsfile,'metrics2d','metrics3d','allens','stateInfo');
end


end