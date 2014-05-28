function [metrics2d, metrics3d, allens, stateInfo]=cemTracker(scen, options)
% Multi-Target Tracking by Continuous Energy Minimization
%
%
%
% STATE  VECTOR
% There are two ways in which the state is represented.
% (1) compact
% One is a d-x-1 vector, containing x,y coordinates in
% the order
%
% Target 1        Target 2        Target N
% x1,y1,x2,y2,...,x1,y1,x2,y2,...,xN,yN
%
% (2) visual
% The other representation consists of two matrices
% for x and y locations respectively
% row i column j corresponds to target j at frame i
%
%
%
% This code contains modifications compared to the 
% one that was used to produce results for our 
% CVPR 2011, VS 2011 and PAMI 2014 papers
%
%
% (C) Anton Milan, 2010-2014
%
% The code may be used free of charge for non-commercial and
% educational purposes, the only requirement is that this text is
% preserved within the derivative work. For any other purpose you
% must contact the authors for permission. This code may not be
% redistributed without written permission from the authors.



% clear all; clear global
% set start time
% clear all
% clear global

global cemStartTime
cemStartTime=tic;


if ~isdeployed,    addpath(genpath('.')); end
addPaths;
homefolder=getHomeFolder();


% addpath('../dctracking')
% addpath('../dctracking/splinefit')
%% declare global variables
global detMatrices detections sceneInfo opt globiter gtInfo scenario
globiter=0;

global LOG_allens LOG_allmets2d LOG_allmets3d %for debug output
global stateInfoTMP % debugging
global allmets3d allmets2d
allmets3d=[]; allmets2d=[]; LOG_allens=[];

%% setup options and scene
% what dataset/sequence?
scenario=80;
if nargin, scenario=scen; end

% fill options struct

% fill options struct with default if not given as parameter
% opt=getDCOptions;
if nargin<2, options='config/default2d.ini'; end
if isstruct(options)
    opt=options;
elseif ischar(options)
    opt=readConOptions(options);
%         opt
else
    error('options parameters not recognized')
end


startsol=opt.startsol;

% what dataset/sequence?
scenario=80;
if nargin, scenario=scen; end
% fill scene info

sceneInfo=getSceneInfo(scenario);

frames=1:length(sceneInfo.frameNums);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% frames=10:20; % do a part of the whole sequence
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if exist('./doframes.txt','file'), frl=load('./doframes.txt'); frames=frl(1):frl(2); end
if exist('/gris/gris-f','dir') && ~exist('frl','var'), frames=1:length(sceneInfo.frameNums); end % if remote and no file

if isfield(opt,'frames'), frames=opt.frames; end
if length(sceneInfo.frameNums)<frames(end), frames=1:length(sceneInfo.frameNums); end % if section outside range

sceneInfo.frameNums=sceneInfo.frameNums(frames);

% detections=detections(frames);
% for t=1:length(detections), detections(t).sc=1; end

%% remove unnecessary frames from GT
gtInfo=cropFramesFromGT(sceneInfo,gtInfo,frames,opt);

% if scenario==97
%     frames=1:3:750;
%     gtInfo.frameNums=1:3:750;
% end
%% load detections
[detections, nPoints]=parseDetections(sceneInfo,frames,opt.detThreshold);

%%%%%%%%%%%%%
% gt detections
% fprintf('!!! WARNING: USING GT DETECTIONS! For debugging only !!!\n');
% [Fgt Ngt]=size(gtInfo.X);
% detections=[];
% for t=1:Fgt
% 	exgt=find(gtInfo.X(t,:));
% 	
%     detections(t).xw=gtInfo.Xgp(t,exgt);	detections(t).yw=gtInfo.Ygp(t,exgt);
%     
%     detections(t).xi=gtInfo.Xi(t,exgt);	detections(t).yi=gtInfo.Yi(t,exgt);
%     detections(t).wd=gtInfo.W(t,exgt);	detections(t).ht=gtInfo.H(t,exgt);
% 	detections(t).sc=ones(1,length(exgt));
%     detections(t).xp=gtInfo.Xi(t,exgt);	detections(t).yp=gtInfo.Yi(t,exgt);
%     if opt.track3d
%         detections(t).xp=gtInfo.Xgp(t,exgt);	detections(t).yp=gtInfo.Ygp(t,exgt);
%     end
% end
% nPoints=numel(find(gtInfo.X(:)));
%%%%%%%%%%%%%%
% 

if scenario>190 || scenario==41
    for t=1:length(detections)
%         detections(t).sc(:)=1;
    end
end
% if scenario==97
%     [detections nPoints]=parseDetections(sceneInfo,1:3:750,opt.detThreshold);
% end
for t=1:length(detections)
%     detections(t).sc=detections(t).sc - opt.detThreshold;
%     detections(t).sc=detections(t).sc ./ max([detections(:).sc]);
end

% for t=1:length(detections)
%     detections(t).sc=detections(t).sc/2;
%     detections(t).sc=detections(t).sc +0.5;
% end


% for t=1:length(detections)
%      detections(t).dirx=detections(t).dirgtxw;
%      detections(t).diry=detections(t).dirgtyw;
%      detections(t).dirxw=detections(t).dirgtxw;
%      detections(t).diryw=detections(t).dirgtyw;
% %     detections(t).sc=detections(t).sc ./ max([detections(:).sc]);
% end
[detections, nPoints]=cutDetections(detections,nPoints,sceneInfo,opt);
% detMatrices=getDetectionMatrices(detections);

% load('/home/aanton/diss/others/yangbo/TUD/TUD_Stadtmitte.avi.detection.mat')
% load('/home/aanton/diss/others/yangbo/PETS09/PETS09_View001_S2_L1_000to794.avi.detection.mat');
F=size(detections,2);
stateInfo.F=F;                          % number of frames

if isempty([detections(:).xi])
    fprintf('no detections present\n');
    [metrics2d metrics3d m2i m3i addInfo2d addInfo3d]=getMetricsForEmptySolution();
    
stateInfo.Xi=zeros(F,0);stateInfo.Yi=zeros(F,0);stateInfo.Xgp=zeros(F,0);stateInfo.Ygp=zeros(F,0);stateInfo.X=zeros(F,0);stateInfo.Y=zeros(F,0);
    stateInfo.X=zeros(F,0);stateInfo.Y=zeros(F,0); stateInfo.W=zeros(F,0); stateInfo.H=zeros(F,0);
    stateInfo.stateVec=[];
    stateInfo.sceneInfo=sceneInfo;    stateInfo.opt=opt;
    stateInfo.splines=[];    stateInfo.outlierLabel=0;    stateInfo.labeling=[];
    
    allens=zeros(1,7);
    return;
end

% experimental, remove certain detections
% Field = fieldnames(detections);
% nDets=0;
% for t=1:F
%     tokeep=find(detections(t).sc>=0.6);
%     nDets=nDets+length(tokeep);
%
%     for iField = 1:length(Field)
%         fcontent=detections(t).(char(Field(iField)));
%         fcontent=fcontent(tokeep);
%         detections(t).(char(Field(iField)))=fcontent;
%     end
% end

detMatrices=getDetectionMatrices(detections);

%% top image limit
sceneInfo.imTopLimit=min([detections(:).yi]);
if scenario>300 && scenario<400
    sceneInfo.imTopLimit = sceneInfo.trackingArea(3);
end
sceneInfo=computeImBordersOnGroundPlane(opt,sceneInfo,detections);
% if ~opt.track3d, sceneInfo.trackingArea(2)
% sceneInfo.trackingArea(3:4)=[320 380];

%% compute histograms if appearance is used
global allhist
allhist=precomputeHistograms(scenario,opt,sceneInfo);


%% print header
printHeader(sceneInfo,scenario,startsol);
printSceneInfo;
printParams(opt);

% evaluate detections
% evaluateDetections(detMatrices,gtInfo);

%% init solution
% ...
% initsolfile=fullfile(getHomeFolder,'diss','ekftracking','output','s0080','e00012d.mat');
% initsolfile=fullfile(getHomeFolder,'diss','ekftracking','output','s0080','e0003.mat');
% initsolfile=fullfile(getHomeFolder,'diss','ekftracking','output','s0023','e0003.mat');
% initsolfile=fullfile(getHomeFolder,'diss','ekftracking','output','s0042','e0003.mat');

X=[]; Y=[];
frtokeep=frames;

X0=X; Y0=Y;
if startsol<6
    initsolfile=fullfile(opt.EKFDir,sprintf('s%04d',scenario),sprintf('e%04d.mat',startsol));
    if exist(initsolfile,'file')
        load(initsolfile);
        printMessage(2,'EKF file loaded...\n');
        EKFlength=size(X,1);
        X0=X; Y0=Y;
        frtokeep=intersect(1:EKFlength,frames);
        X0=X(frtokeep,:); Y0=Y(frtokeep,:);
    end
    startPT=[];
    startPT.frameNums=sceneInfo.frameNums;startPT.F=length(startPT.frameNums);

%     startPT=cropFramesFromGT(sceneInfo,startPT,frames,opt);
    
    startPT.X=X0;startPT.Y=Y0;
    
    
    [startPT.X, startPT.Y, startPT]=cleanState(X0, Y0, startPT);
    [startPT.X, startPT.Y]=checkInitSolution(startPT.X,startPT.Y,startPT.F);
        
    if opt.track3d
        startPT.Xgp=startPT.X; startPT.Ygp=startPT.Y;
        [startPT.Xi startPT.Yi]=projectToImage(startPT.Xgp,startPT.Ygp,sceneInfo);
    else
        startPT.Xi=startPT.X; startPT.Yi=startPT.Y;
    end
    
   


    if opt.track3d && opt.cutToTA,        startPT=cutStateToTrackingArea(startPT);    end
    [startPT.X, startPT.Y, startPT]=cleanState(startPT.X, startPT.X, startPT);
    
    startPT=getBBoxesFromState(startPT);
   
    fprintf('EKF Result: \n');

    [metrics2d, metrics3d, addInfo2d, addInfo3d]=printFinalEvaluation(startPT, gtInfo, sceneInfo, opt);
    

end


% for startsol=1
%     initsolfile=fullfile(getHomeFolder,'diss','ekftracking','output',sprintf('s%04d',scenario),'e0001.mat');
%     if exist(initsolfile,'file')
%         load(initsolfile);
%         X0=[X0 X]; Y0=[Y0 Y];
%     end
% end

% initsolfile=fullfile(getHomeFolder,'diss','others','okuma','BPF','myresult.mat');
% initsolfile=fullfile(getHomeFolder,'diss','others','okuma','BPF',sprintf('myresult-s%04d.mat',scenario));
% initsolfile='tmp.mat';
if startsol==6
    %% Pirsiavash
%     initsolfile=sprintf('%s/startPT-pir-s%04d.mat',opt.DPDir,scenario);
        
    % just compute it on the fly
    pOpt=getPirOptions;
    [metrics2d, metrics3d, allene, startPT]=runDP(scenario,pOpt,opt);
    
%     if exist(initsolfile,'file')
%         load(initsolfile);
        
        if opt.track3d
            [startPT.X,startPT.Y]=projectToGroundPlane(startPT.Xi,startPT.Yi,sceneInfo);
            startPT.Xgp=startPT.X;startPT.Ygp=startPT.Y;
        end
        startPT.X=startPT.X;startPT.Y=startPT.Y;
        
        if opt.track3d && opt.cutToTA,        startPT=cutStateToTrackingArea(startPT,sceneInfo, opt);    end
        startPT=cropFramesFromGT(sceneInfo,startPT,frames,opt);
        
        fprintf('Pirsiavash Result: \n');
        [metrics2d, metrics3d, addInfo2d, addInfo3d]=printFinalEvaluation(startPT, gtInfo, sceneInfo, opt);
        
        X=startPT.Xi;Y=startPT.Yi;
        if opt.track3d
            [X Y]=projectToGroundPlane(startPT.Xi,startPT.Yi,sceneInfo);
        end
        X0=[X0 X]; Y0=[Y0 Y];
        
        
        % X=X(fr,:);Y=Y(fr,:);H=H(fr,:); % !!!!!! REMOVE
%     else
%         warning('DP result file not found: %s\n',initsolfile);
%     end
    
end

%% ground truth
if startsol==7
    fprintf('!!! WARNING: STARTING FROM GROUND TRUTH! For debugging only !!!\n');
    X0=gtInfo.X;Y0=gtInfo.Y;
    startPT.Xi=gtInfo.Xi;startPT.Yi=gtInfo.Yi;startPT.H=gtInfo.H;startPT.W=gtInfo.W;
    if opt.track3d
        X0=gtInfo.Xgp;Y0=gtInfo.Ygp;
        startPT.Xgp=X0; startPT.Ygp=Y0;
    end
    startPT.X=X0;startPT.Y=Y0;startPT.frameNums=gtInfo.frameNums;
    
    
    if opt.track3d && opt.cutToTA,        startPT=cutStateToTrackingArea(startPT);    end
    [metrics2d metrics3d addInfo2d addInfo3d]=printFinalEvaluation(startPT, gtInfo, sceneInfo, opt);

end

X=X0; Y=Y0;
if ~isempty(X)
    [X, Y, stateInfo]=cleanState(X, Y, stateInfo);
    [X Y]=checkInitSolution(X,Y,stateInfo.F);
end


% sampling
% alldpoints=createAllDetPoints;
% nInit=round(nDets/F);
% mhs=getSplineProposals(alldpoints,nInit,F);
% for id=1:length(mhs)
%     tt=mhs(id).start:mhs(id).end;
%     allvals=ppval(mhs(id),tt);
%     X(tt',id)=allvals(1,:)';Y(tt',id)=allvals(2,:)';
% end
% [X, Y, stateInfo]=cleanState(X, Y, stateInfo);
% [X Y]=checkInitSolution(X,Y,stateInfo.F);






% global gtInfo
% X=gtInfo.X; Y=gtInfo.Y;

if isempty(X), X=[]; end
stateInfo.N=size(X,2);
stateInfo.targetsExist=getTracksLifeSpans(X);
stateInfo.frameNums=sceneInfo.frameNums;
stateInfo=matricesToVector(X,Y,stateInfo);
[stateVec N F targetsExist stateInfo.X stateInfo.Y]=getStateInfo(stateInfo);

assert(isempty(X) || size(X,1)==length(detections),'Initial Solution must be full length');


% stateInfo.X
%% initial gradient descent (if initial solution available)
if ~isempty(stateInfo.stateVec)
    [stateInfo.stateVec Evalue nIterations]=minimize(stateInfo.stateVec,'E',opt.maxIterCGD,stateInfo);
end

%% now do main optimization
global allshr allgrow % HACK FOR NOW
converged=false;
epoch=0;
while ~converged && epoch<opt.maxEpochs
    globiter=globiter+1;
    if opt.verbosity>=3 && stateInfo.N>0
        printProgressHeader(scenario,sceneInfo,opt);
        [EdetValue EdynValue EexcValue EappValue EperValue EregValue,EoriValue]= ...
            printEnergies(opt, 0, stateInfo.stateVec, stateInfo);
        LOG_allens(globiter,:)=[EdetValue EdynValue EexcValue EappValue EperValue EregValue,EoriValue];
        printProgressMetrics(opt,sceneInfo,stateInfo.stateVec, stateInfo, 0);
        fprintf('\n');
    end
        
    epoch=epoch+1;
    printMessage(1,'---- JUMP  MOVES  ROUND %i ----\n',epoch);
    jumpExecuted=false(1,6);
    
    for jumpMove=opt.jumpsOrder
        stateInfoOld=stateInfo;
        eval(sprintf('stateInfo=%s(stateInfo);',getJumpMoveFunction(jumpMove)))
        
        % did we jump?
        if isequal(stateInfoOld,stateInfo) % no
            printMoveFailure(jumpMove)
            
            %  perform conjugate gradient descent if move was successful
        else
            jumpExecuted(jumpMove)=1;
            [stateInfo.stateVec Evalue nIterations]=minimize(stateInfo.stateVec,'E',opt.maxIterCGD,stateInfo);
        end
        stateInfoTMP=stateInfo;
        globiter=globiter+1;
    end
    
    % HACK FOR NOW
    if isequal(size(allshr),size(allgrow))
        if any(allshr(:)) && isequal(allshr,allgrow)
            converged=true;
            printMessage(1,'stuck in a loop after %i epochs.\n',epoch);            
        end
    end
%     if opt.verbosity>=3
%         printProgressHeader(scenario,sceneInfo,opt);
%         [EdetValue EdynValue EexcValue EappValue EperValue EregValue,EoriValue]= ...
%             printEnergies(opt, 0, stateInfo.stateVec, stateInfo);
%         LOG_allens(globiter,:)=[EdetValue EdynValue EexcValue EappValue EperValue EregValue,EoriValue];
%         printProgressMetrics(opt,sceneInfo,stateInfo.stateVec, stateInfo, 0);
%         fprintf('\n');
%     end
        
        
    % if no jumps were perfomed, we're done
    if all(~jumpExecuted)
        converged=true;
        printMessage(1,'No jumps were executed. Optimization has converged after %i epochs.\n',epoch);
    end
    
    % also, if we have ground truth and MOTA is negative or
    % much lower than when we started, abort
    % this clearly means that the energy is wrong
    if sceneInfo.gtAvailable
%         if (opt.track3d && allmets3d(end,12)<max(allmets3d(:,12))/2)
%             converged=true;
%             printMessage(1, 'Current result is very poor. Abort optimization.\n');
%         end
%         if (~opt.track3d && allmets2d(find(allmets2d(:,12),1,'last'),12)<max(allmets2d(:,12))/2)
%             converged=true;
%             printMessage(1, 'Current result is very poor. Abort optimization.\n');
%         end
    end
    
    
    % if last epoch
    if epoch>=opt.maxEpochs
        printMessage(1,'Max number of rounds reached.\n');
    end
end % converged
% basically we are done
printMessage(1,'All done (%.2f min = %.2fh = %.2f sec per frame)\n',toc(cemStartTime)/60,toc(cemStartTime)/3600,toc(cemStartTime)/F);


[stateVec N F targetsExist stateInfo.X stateInfo.Y]=getStateInfo(stateInfo);
% stateInfo.X
% stateInfo
if ~isempty(intersect(scenario,[195:199]))
    [a b c]=intersect(stateInfo.frameNums,gtInfo.frameNums);
    
    frkeep=b';
    stateInfo.frameNums=stateInfo.frameNums(frkeep);
    stateInfo.X=stateInfo.X(frkeep,:);stateInfo.Y=stateInfo.Y(frkeep,:);
    [stateInfo.F stateInfo.N]=size(stateInfo.X);
    stateInfo.targetsExist=getTracksLifeSpans(stateInfo.X);    
%     stateInfo.X
%     stateInfo.Y
    
    stateInfo=matricesToVector(stateInfo.X,stateInfo.Y,stateInfo);
    
            sceneInfo.frameNums=stateInfo.frameNums;
            [detections nPoints]=parseDetections(sceneInfo,sceneInfo.frameNums,opt.detThreshold);
    
    
%     stateInfo
%     stateInfo.W=stateInfo.W(frkeep,:);stateInfo.H=stateInfo.H(frkeep,:);
%     if options.eval3d
%         stateInfo.Xgp=stateInfo.Xgp(frkeep,:);stateInfo.Ygp=stateInfo.Ygp(frkeep,:);
%     end
%     stateInfo.Xi=stateInfo.Xi(frkeep,:);stateInfo.Yi=stateInfo.Yi(frkeep,:);


end

%% post processing
% get X Y matrices
[stateVec N F targetsExist stateInfo.X stateInfo.Y]=getStateInfo(stateInfo);
stateInfo=postProcessState(stateInfo);
% stateInfo
% pause

if ~isempty(intersect(scenario,[195:199]))
    [a b c]=intersect(stateInfo.frameNums,gtInfo.frameNums);
        if opt.visOptim

            
            prepFigure();
            plotDetections(1);
            plotTrajectories(stateInfo.stateVec,stateInfo,1); pause(.001);
        end
end

%% if we have ground truth, evaluate results
[metrics2d, metrics3d]=printFinalEvaluation(stateInfo, gtInfo, sceneInfo, opt);

% TODO: LOGGING!
% %% clean up (remove zero rows from logs)
% itinfo=find(sum(LOG_allens,2));
% LOG_allmets2d=LOG_allmets2d(~~sum(LOG_allmets2d,2),:);
% LOG_allmets3d=LOG_allmets3d(~~sum(LOG_allmets3d,2),:);
% LOG_allens=LOG_allens(~~sum(LOG_allens,2),:);

if isempty(LOG_allens), LOG_allens=zeros(1,7); end
allens=LOG_allens;
allens=allens(end,:);

stateInfo.sceneInfo=sceneInfo;
stateInfo.opt=opt;

% you can display the results with
% displayTrackingResult(sceneInfo,stateInfo)

end
