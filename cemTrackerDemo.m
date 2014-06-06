% function cemTracker(detfile, opt)
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



clear all; clear global
% set start time


% download demo sequence
% if ~exist('demo','dir'), mkdir('demo'); end
if ~exist('demo.zip','file')
    fprintf('Downloading demo sequence (ca. 6MB) ...\n');
    try urlwrite('http://research.milanton.net/files/demo-tracking.zip','demo.zip')
    catch err, error('Download failed! %s',err.identifier);
    end
end

unzip('demo.zip');



global cemStartTime
cemStartTime=tic;

addpath(genpath('.'))
%% declare global variables
global detMatrices detections sceneInfo opt globiter gtInfo scenario
globiter=0; scenario=80;

global LOG_allens LOG_allmets2d LOG_allmets3d %for debug output

%% setup options and scene
% fill options struct
opt=readConOptions('config/default3d.ini');

% fill scene info
sceneInfo=getSceneInfoConDemo;
frames=1:length(sceneInfo.frameNums);
sceneInfo.scenario=scenario;

%% cut ground truth
gtInfo=cutGTToTrackingArea(gtInfo,sceneInfo);

%% load detections
[detections nDets]=parseDetections(sceneInfo); fr=1:length(detections);
stateInfo.F=size(detections,2);

% cut detections to tracking area if needed
[detections nDets]=cutDetections(detections,nDets,sceneInfo, opt);
detMatrices=getDetectionMatrices(detections);


%% init solution
X=[]; Y=[];
if opt.startsol==6
    %% Pirsiavash
    %     initsolfile=sprintf('%s/startPT-pir-s%04d.mat',opt.DPDir,scenario);
    
    % just compute it on the fly
    pOpt=getPirOptions;
    [metrics2d, metrics3d, allene, startPT]=runDP(scenario,pOpt,opt);
    
    if opt.track3d
        [startPT.X,startPT.Y]=projectToGroundPlane(startPT.Xi,startPT.Yi,sceneInfo);
        startPT.Xgp=startPT.X;startPT.Ygp=startPT.Y;
    end
%     startPT.X=startPT.Xi;startPT.Y=startPT.Yi;
    
    if opt.track3d && opt.cutToTA,        startPT=cutStateToTrackingArea(startPT,sceneInfo, opt);    end
    startPT=cropFramesFromGT(sceneInfo,startPT,frames,opt);
    
    fprintf('Pirsiavash Result: \n');
    [metrics2d, metrics3d, addInfo2d, addInfo3d]=printFinalEvaluation(startPT, gtInfo, sceneInfo, opt);
    
    X=startPT.Xi;Y=startPT.Yi;
    if opt.track3d
        [X Y]=projectToGroundPlane(startPT.Xi,startPT.Yi,sceneInfo);
    end    
end


if ~isempty(X)
    [X, Y, stateInfo]=cleanState(X, Y, stateInfo);
    [X Y]=checkInitSolution(X,Y,stateInfo.F);
end

stateInfo.N=size(X,2);
stateInfo.targetsExist=getTracksLifeSpans(X);
stateInfo.frameNums=sceneInfo.frameNums;
stateInfo=matricesToVector(X,Y,stateInfo);
[stateVec N F targetsExist stateInfo.X stateInfo.Y]=getStateInfo(stateInfo);

assert(size(X,1)==length(detections),'Initial Solution must be full length');
printSceneInfo;
%% initial gradient descent (if initial solution available)
if ~isempty(stateInfo.stateVec)
    [stateInfo.stateVec Evalue nIterations]=minimize(stateInfo.stateVec,'E',opt.maxIterCGD,stateInfo);
end

%% now do main optimization
converged=false;
epoch=0;
while ~converged && epoch<opt.maxEpochs
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
    end
    
    % if no jumps were perfomed, we're done
    if all(~jumpExecuted)
        converged=true;
        printMessage(1,'No jumps were executed. Optimization has converged after %i epochs.\n',epoch);
    end
    
    % if last epoch
    if epoch>=opt.maxEpochs
        printMessage(1,'Max number of rounds reached.\n');
    end
end % converged
% basically we are done
printMessage(1,'All done (%.2f min = %.2fh = %.2f sec per frame)\n',toc(cemStartTime)/60,toc(cemStartTime)/3600,toc(cemStartTime)/F);




%% post processing
% get X Y matrices
[stateVec N F targetsExist stateInfo.X stateInfo.Y]=getStateInfo(stateInfo);
stateInfo=postProcessState(stateInfo);



%% if we have ground truth, evaluate results
printFinalEvaluation(stateInfo, gtInfo, sceneInfo, opt);

%% clean up (remove zero rows from logs)
itinfo=find(sum(LOG_allens,2));
LOG_allmets2d=LOG_allmets2d(~~sum(LOG_allmets2d,2),:);
LOG_allmets3d=LOG_allmets3d(~~sum(LOG_allmets3d,2),:);
LOG_allens=LOG_allens(~~sum(LOG_allens,2),:);

stateInfo.sceneInfo=sceneInfo;
stateInfo.opt=opt;

% you can display the results with
displayTrackingResult(sceneInfo,stateInfo)

% end