function [metrics2d metrics3d]=cemTrackerSearch(scenario,options)
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



% clear all
global cemStartTime
cemStartTime=tic;

addpath(genpath('.'))

%% declare global variables
global detMatrices detections sceneInfo opt globiter;
globiter=0;

% fill options struct
opt=getConOptions(options);

% fill scene info
% scenario=23;
sceneInfo=getSceneInfo(scenario);


%% load detections
detections=parseDetections(sceneInfo); fr=1:length(detections);
% fr=1:50; detections=detections(fr);  % !!!!!!! REMOVE
% for t=1:length(detections), detections(t).sc=1; end

% load('/home/aanton/diss/others/yangbo/TUD/TUD_Stadtmitte.avi.detection.mat')
% load('/home/aanton/diss/others/yangbo/PETS09/PETS09_View001_S2_L1_000to794.avi.detection.mat');
F=size(detections,2);
stateInfo.F=F;                          % number of frames

% do 2d
% for t=1:F
%     detections(t).xp=detections(t).xi;
%     detections(t).yp=detections(t).yi;
% end

% for t=1:F
%     detections(t).sc(:)=1;
% end

detMatrices=getDetectionMatrices(detections);

%% init solution
% ...
% initsolfile=fullfile(getHomeFolder,'diss','ekftracking','output','s0080','e00012d.mat');
% initsolfile=fullfile(getHomeFolder,'diss','ekftracking','output','s0080','e0003.mat');
% initsolfile=fullfile(getHomeFolder,'diss','ekftracking','output','s0023','e0003.mat');
% initsolfile=fullfile(getHomeFolder,'diss','ekftracking','output','s0042','e0003.mat');

X=[]; Y=[];

initsolfile=fullfile(getHomeFolder,'diss','ekftracking','output',sprintf('s%04d',scenario),'e0003.mat');
% initsolfile=fullfile(getHomeFolder,'diss','others','okuma','BPF','myresult.mat');
% initsolfile=fullfile(getHomeFolder,'diss','others','okuma','BPF',sprintf('myresult-s%04d.mat',scenario));
% initsolfile='tmp.mat';
if 1%exist(initsolfile,'file')
load(initsolfile);

% X=X(fr,:);Y=Y(fr,:);H=H(fr,:); % !!!!!! REMOVE
[X, Y, stateInfo]=cleanState(X, Y, stateInfo);
[X Y]=checkInitSolution(X,Y,stateInfo.F);
end

% global gtInfo
% X=gtInfo.X; Y=gtInfo.Y;

stateInfo.N=size(X,2);
stateInfo.targetsExist=getTracksLifeSpans(X);
stateInfo.frameNums=sceneInfo.frameNums;

stateInfo=matricesToVector(X,Y,stateInfo);
[stateVec N F targetsExist stateInfo.X stateInfo.Y]=getStateInfo(stateInfo);

% pause

%% do we have ground truth?
% ...
% global gtInfo
% gtInfo=parseGT('/storage/databases/PETS2009/Crowd_PETS09/S3/Multiple_Flow/Time_12-43/View_001/GT2d_full.xml');
% sceneInfo.gtAvailable=1;





% printSceneInfo;

%% initial gradient descent
stateVec=stateInfo.stateVec;
[stateVec Evalue nIterations]=minimize(stateVec,'E',opt.maxIterCGD,stateInfo);
stateInfo.stateVec=stateVec;

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
%         if opt.track3d, 
%             [~,~,~,~,stateInfo.X stateInfo.Y]=getStateInfo(stateInfo);
%             stateInfo=cutStateToTrackingArea(stateInfo); 
%         end
         
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

if opt.track3d
    stateInfo=cutStateToTrackingArea(stateInfo);
end

% if we tracked on image, Xi = X
if ~opt.track3d
    stateInfo.Xi=stateInfo.X; stateInfo.Yi=stateInfo.Y;
% otherwise project back
else
    stateInfo.Xgp=stateInfo.X; stateInfo.Ygp=stateInfo.Y;
    [stateInfo.Xi stateInfo.Yi]=projectToImage(stateInfo.X,stateInfo.Y,sceneInfo);
end

%% get bounding boxes from corresponding detections
stateInfo=getBBoxesFromState(stateInfo);

% if we have ground truth, evaluate results
if sceneInfo.gtAvailable
    global gtInfo
    printMessage(1,'\nEvaluation 2D:\n');
    [metrics2d metricsInfo2d]=CLEAR_MOT(gtInfo,stateInfo);
    printMetrics(metrics2d,metricsInfo2d,1);
    
    if opt.track3d
        printMessage(1,'\nEvaluation 3D:\n');
        evopt.eval3d=1;
        [metrics3d metricsInfo3d]=CLEAR_MOT(gtInfo,stateInfo,evopt);
        printMetrics(metrics3d,metricsInfo3d,1);
        
    end
end


% you can display the results with
% displayTrackingResult(sceneInfo,stateInfo)

% end