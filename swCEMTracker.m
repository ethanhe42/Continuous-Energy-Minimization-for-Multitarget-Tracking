function [metrics2d, metrics3d, allens, stateInfo]=swCEMTracker(scen,options,swfile)

% scen=71;
global scenario
scenario=scen;

% default
windowSize=50;minWindowSize=15;overlapSize=5;
windowSize=20;minWindowSize=5;overlapSize=2;

if nargin<3, swfile='swparams.txt'; end
if exist(swfile,'file')
    swparams=load(swfile);
    windowSize=swparams(1);minWindowSize=swparams(2);overlapSize=swparams(3);
end

allstInfo=[];

if ~isdeployed
    addpath(genpath('./utils/'));  
end

global opt;
opt=getConOptions;
opti=opt;
if nargin==2
    opti=options;
end
sceneInfo=getSceneInfo(scenario);
allframeNums=sceneInfo.frameNums;
F=length(allframeNums);
fromframe=1; toframe=windowSize;

wincnt=0;
allwins=[];
allallens=zeros(0,7);
while toframe<=F
    
    wincnt=wincnt+1;
    fprintf('from %4d to %4d = %4d\n',fromframe,toframe,length(fromframe:toframe));
    
    opti.frames=fromframe:toframe;
    [metrics2d, metrics3d, allens, stateInfo]=cemTracker(scen,opti);
    allallens(wincnt,:)=allens;
%     allstInfo
%     stateInfo
    allstInfo=[allstInfo stateInfo];
    allwins(wincnt,:)=[fromframe toframe];

    % now adjust new time frame
    % if already at end, break
    if toframe==F,        break;    end
    
    % otherwise slide window and make bigger if needed
    % at the end
    fromframe=fromframe+windowSize-overlapSize;    
    newend=toframe+windowSize-overlapSize;
    if newend > F-minWindowSize
        toframe=F;
    else
        toframe=newend;
    end
end

%% now stich em all together (just stack)

X=zeros(F,0);Y=zeros(F,0);Xgp=zeros(F,0);Ygp=zeros(F,0);
Xi=zeros(F,0);Yi=zeros(F,0);W=zeros(F,0);H=zeros(F,0);

curN=1;
for w=1:wincnt
    winframes=(allwins(w,1):allwins(w,2))';    
    N=size(allstInfo(w).X,2);
    newIDs=curN:curN+N-1;
    X(winframes,newIDs)=0;   
    
    X(winframes,newIDs)=allstInfo(w).X;Y(winframes,newIDs)=allstInfo(w).Y;
    Xi(winframes,newIDs)=allstInfo(w).Xi;Yi(winframes,newIDs)=allstInfo(w).Yi;
    if isfield(allstInfo,'Xgp')
        Xgp(winframes,newIDs)=allstInfo(w).Xgp;Ygp(winframes,newIDs)=allstInfo(w).Ygp;
    end
    W(winframes,newIDs)=allstInfo(w).W;H(winframes,newIDs)=allstInfo(w).H;
    
    curN=curN+N;
end

%% finish up
global scenario
global gtInfo
scenario=scen;
sceneInfo=getSceneInfo(scenario);

stInfo.sceneInfo=sceneInfo;
opti.frames=allframeNums;
stInfo.opt=opti;
stInfo.frameNums=allframeNums;
stInfo.F=length(stInfo.frameNums);
stInfo.X=X;stInfo.Y=Y;
stInfo.Xi=Xi;stInfo.Yi=Yi;stInfo.W=W;stInfo.H=H;
if isfield(allstInfo,'Xgp') ,stInfo.Xgp=Xgp;stInfo.Ygp=Ygp; end

stateInfo=stitchTemporalWindows(allstInfo,stInfo,allwins,overlapSize);

opt.frames=1:length(stateInfo.frameNums);
[metrics2d, metrics3d]=printFinalEvaluation(stateInfo, gtInfo, sceneInfo, stInfo.opt);

allens=mean(allallens);
%%
% clf
% hold on
% wx=-13000;
% for w=1:wincnt
%     plot(allwins(w,1):allwins(w,2),allstInfo(w).X,'.');
%     line([allwins(w,1)-.5 allwins(w,1)-.5],[wx wx+10000])
%     line([allwins(w,2)+.5 allwins(w,2)+.5],[wx wx+10000])
%     line([allwins(w,1)-.5 allwins(w,2)+.5],[wx wx])
%     line([allwins(w,1)-.5 allwins(w,2)+.5],[wx+10000 wx+10000])
%     wx=wx+100;
% end

end
