function trainModel(jobname,jobid,maxexper)

%% determine paths for config, logs, etc...
addpath(genpath('../motutils/'));
format compact

[~,hname]=system('hostname')
settingsDir=strsplit(jobname,'-');
runname=char(settingsDir{1})
learniter=char(settingsDir{2})
jid=char(settingsDir{3}); % unused

settingsDir=[runname '-' learniter];

confdir=sprintf('config/%s',settingsDir);

jobid=str2double(jobid);
confdir

resdir=sprintf('results/%s',settingsDir);
if ~exist(resdir,'dir'), mkdir(resdir); end
resdir

trainStartTime=tic;
resultsfile=sprintf('%s/res_%03d.mat',resdir,jobid);

% if computed alread, just load it
if exist(resultsfile,'file')
  load(resultsfile);
else

  conffile=fullfile(confdir,sprintf('%04d.ini',jobid));
  conffile

  inifile=fullfile(confdir,'0001.ini');
  inifile
  if ~exist(inifile,'file')
	  error('You must provide initial options file 0001.ini');
  end

  opt=readConOptions(inifile);

  %% zip up logs (previous)
  if jobid==1
    prevSetting=sprintf('%s-%d',runname,str2double(learniter)-1)
    zipstr=sprintf('!sh ziplogs.sh %s',prevSetting)
    eval(zipstr);
  end


  % take care of parameters
  jobid

  rng(jobid);
  % if jobid==1, leave as is
  if jobid==1
  % otherwise, randomize and write out new ini file
  else
	  ini=IniConfig();
	  ini.ReadFile(inifile);

	  params=[];
	  % we are only interested in [Parameters]
	  sec='Parameters';
	  keys = ini.GetKeys(sec);

	  for k=1:length(keys)
	      key=char(keys{k});
	      %params = setfield(params,key,ini.GetValues(sec,key));
          val=ini.GetValues(sec,key);
          
          % e-notation is read in as a string
          if ischar(val), val=str2double(val); end
	      params = [params val];
	  end

	  rnmeans = params; % mean values are the starting point
	  rmvars = rnmeans ./ 10; % variance is one tenth

	  if jobid <= maxexper/2
		  params = 2*rand(1,length(params)) .* rnmeans; % uniform [0, 2*max]
	  else
		  params = abs(rnmeans + rmvars .* randn(1,length(rnmeans))); % normal sampling
	  end


	  for k=1:length(keys)
	      key=char(keys{k});
	      %params = setfield(params,key,ini.GetValues(sec,key));
	      opt = setfield(opt, key, params(k));
	  end

	  % write out new opt file
	  status = writeConOptions(opt,conffile);
	  

	  
  end
  rng(1);

  allscensfile=fullfile(confdir,'doscens.txt');
  if ~exist(allscensfile)
	  warning('Doing the standard PETS TUD combo...');
	  dlmwrite(fullfile(confdir,'doscens.txt'),[23 25 27 71 72 42]);
  end
  allscen=dlmread(fullfile(confdir,'doscens.txt'));
  allscen

  learniter=str2double(learniter)

  mets2d=zeros(max(allscen),14);
  mets3d=zeros(max(allscen),14);  
  ens=zeros(max(allscen),5);

  for scenario=allscen
	  fprintf('jobid: %d,   learn iteration %d\n',jobid,learniter);
	  scensolfile=sprintf('%s/prt_res_%03d-scen%02d.mat',resdir,jobid,scenario)

	  try
	    load(scensolfile);
	  catch err
	    fprintf('Could not load result: %s\n',err.message);
	    [metrics2d, metrics3d, energies, stateInfo]=cemTracker(scenario,conffile);
	    save(scensolfile,'stateInfo','metrics2d','metrics3d','energies');
	  end	  

	  mets2d(scenario,:)=metrics2d;
	  mets3d(scenario,:)=metrics3d;
	  ens(scenario,:)=double(energies);
	  infos(scenario).stateInfo=stateInfo;

  end

  % KITTI

%    if ~isempty(intersect(scenario,[500:899, 1500:1899]))
    if isfield(opt,'KITTI')
      a=infos(allscen);
      kittiDir=sprintf('%s-%d',settingsDir,jobid)
      metricsKITTI=evalKITTI(a,kittiDir,opt.KITTI);
      metricsKITTI(1:12)=100*metricsKITTI(1:12); % percentages
      KITTIToMineMapping=[6 7 9 19 10 11 12 14 15 16 17 1 2 3]
      myMets=metricsKITTI(KITTIToMineMapping);
      myMets(4:11)=round(myMets(4:11)); % rounding, carefull, also GT,MT, etc...
    
      for scenario=allscen
  	mets2d(scenario,:)=myMets;
      end
      
    end  

    
  save(resultsfile,'opt','mets2d','mets3d','ens','infos','allscen');
  
  
  % remove temp scene files
  for scenario=allscen
    scensolfile=sprintf('%s/prt_res_%03d-scen%02d.mat',resdir,jobid,scenario)
    if exist(scensolfile,'file')
      delete(scensolfile);
    end
  end
end

printMessage(1,'Job done (%.2f min = %.2fh = %.2f sec per sequence)\n', ...
    toc(trainStartTime)/60,toc(trainStartTime)/3600,toc(trainStartTime)/numel(allscen));

% evaluate what we have so far
if isempty(intersect(allscen,1001))
	[bestexper,bestmota]=combineResultsRemote(settingsDir);
else
	[bestexper,bestmota]=combineResultsBenchmark(settingsDir,jobid,maxexper);
end

querystring=sprintf('qstat -t | grep %s | wc -l',settingsDir);
[rs,rjobs] = system(querystring); rjobs=str2double(rjobs)-1; % subtract currently running

fprintf('%d other jobs still running\n',rjobs);

% save to bis allres file
if rjobs==0
    try
        fid=fopen('alltraining.txt','a');
        fprintf(fid,'%s\t%s\t%.1f (%d)\n',datestr(now),settingsDir,bestmota,bestexper);
        fclose(fid);
    catch
    end
end

% if last one, resubmit
if bestexper==1
	fprintf('Training Done!');
else
  if rjobs<=0
	fprintf('resubmitting ... \n');
	newSetting=sprintf('%s-%d',runname,learniter+1)	
	newConfdir=sprintf('config/%s',newSetting)
	cpstr=sprintf('!cp -R %s %s',confdir,newConfdir)
	fprintf('copy config dir\n');
	eval(cpstr);
	
	% copy relevant config into first one
	conffile=sprintf('%s/%04d.ini',confdir,bestexper)
	cpstr=sprintf('!cp %s %s/0001.ini',conffile,newConfdir)
	fprintf('copy best config file');
	eval(cpstr);
	
	% 
	submitstr=sprintf('!ssh -o UserKnownHostsFile=/dev/null -o StrictHostKeyChecking=no moby  \"cd research/projects/contracking; sh submitTrain.sh %s\"',newSetting)
	fprintf('submit: %s\n',newSetting)
  	eval(submitstr);	
  else
    fprintf('waiting for other jobs to finish\n');
  end
end

