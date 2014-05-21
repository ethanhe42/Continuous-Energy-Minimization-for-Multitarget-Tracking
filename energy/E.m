function [fx dfx EdetValue EdynValue EexcValue EappValue EperValue EregValue EoriValue] = E(stateVec, stateInfo)
% the objective function E
%
% inputs:
% state vector, state info and scene info and options and detections
%
% outputs
% Energy value and its derivative vector
% and optionally all individual (not weighted) energy components
%



global sceneInfo opt;

F=stateInfo.F; N=stateInfo.N;
targetsExist=stateInfo.targetsExist;
tiToInd=stateInfo.tiToInd;
% stateVec=stateInfo.stateVec;

% if no targets present, E=0 and return
if ~N
    fx=0;    dfx=0;
    EdetValue=0; EdynValue =0;  EexcValue=0;  EappValue=0;   EperValue=0; EregValue=0; EoriValue=0;
    return;
end


targetSize=sceneInfo.targetSize;
areaLimits=sceneInfo.trackingArea;

% params=getOptimParameters;
% global params;
% global gridStep itToInd Xd Yd Sd;
% global areaLimits;
% [X,Y]=vectorToMatrices(stateVec,stateInfo);
[X, Y]=vectorToMatrices_mex(stateVec,stateInfo.tiToInd,stateInfo.F,stateInfo.N);
% global Xgt Ygt;


% detections
global detMatrices;



% %%%%%%%%%%%%%%
% % Detections %
% %%%%%%%%%%%%%%
EdetValue=0;
dEdet=zeros(length(stateVec),1);
% EdetValue=Edet(x);
% if params.alpha>0
%     if params.det3d
if opt.mex
    % vis needed for Eapp!!!
    if nargout>1 || opt.wtEapp
        %%%%% PRML !!!! OFF
        if sceneInfo.scenario<=400
            [EdetValue, dEdet, visv, vis, visx, visy, ddvix, ddviy]= ...
                Edet_mex(X,Y,detMatrices.Xd,detMatrices.Yd,detMatrices.Sd, ...
                sceneInfo.targetSize,opt.lambda,stateInfo.targetsExist,length(stateVec),opt.occ);
        else
            %                     save('tmpvars.mat','X','Y','detMatrices','sceneInfo','opt','stateInfo','stateVec');
            [EdetValue, dEdet]= ...
                Edet_sig_mex(X,Y,detMatrices.Xd,detMatrices.Yd,detMatrices.Sd, ...
                sceneInfo.targetSize,opt.lambda,stateInfo.targetsExist,length(stateVec),opt.occ);
            %                 pause
        end
        
    else
        if sceneInfo.scenario<=400
            EdetValue=Edet_mex(X,Y,detMatrices.Xd,detMatrices.Yd,detMatrices.Sd, ...
                sceneInfo.targetSize,opt.lambda,stateInfo.targetsExist,length(stateVec),opt.occ);
        else
            [EdetValue, dEdet]= ...
                Edet_sig_mex(X,Y,detMatrices.Xd,detMatrices.Yd,detMatrices.Sd, ...
                sceneInfo.targetSize,opt.lambda,stateInfo.targetsExist,length(stateVec),opt.occ);
        end
    end
else
    if nargout>1, %[EdetValue dEdet ds VIS]=Edet(x);
        [EdetValue dEdet]=Edet(stateVec,stateInfo);
    else EdetValue=Edet(stateVec,stateInfo);
    end
end

% [EdetValuec, dEdetc]= ...
%     Edet_mex(X,Y,detMatrices.Xd,detMatrices.Yd,detMatrices.Sd, ...
%     sceneInfo.targetSize,opt.lambda,stateInfo.targetsExist,length(stateVec),opt.occ);
% [EdetValues, dEdets]= ...
%     Edet_sig_mex(X,Y,detMatrices.Xd,detMatrices.Yd,detMatrices.Sd, ...
%     sceneInfo.targetSize,opt.lambda,stateInfo.targetsExist,length(stateVec),opt.occ);
% [EdetValuec EdetValues EdetValuec-EdetValues]
% pause


%     else
%         if nargout>1, %[EdetValue dEdet ds VIS]=Edet(x);
%             [EdetValue dEdet visv visx visy ddvix ddviy]=Edet2D(x);
%         else EdetValue=Edet2D(x);
%         end
%     end
% end
% [EdetValue dEdet visv visx visy ddvix ddviy ds]=Edet(x);
%
% % pause
% stateInfo.X=X; stateInfo.Y=Y;
% stateInfo.Xgp=stateInfo.X; stateInfo.Ygp=stateInfo.Y;
% [stateInfo.Xi stateInfo.Yi]=projectToImage(stateInfo.X,stateInfo.Y,sceneInfo);
% a=getBBoxesFromPrior(stateInfo);
%
% %%%%%%%%%%%%%%
% % appearance %
% %%%%%%%%%%%%%%
EappValue = 0; dEapp=zeros(length(stateVec),1);
if opt.wtEapp>0
    if nargout>1, [EappValue dEapp]=Eapp(stateVec, stateInfo, vis, ddvix, ddviy);
    else
        EappValue=Eapp(stateVec,stateInfo,vis);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Edyn - constant velocity %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
EdynValue=0; dEdyn=zeros(length(stateVec),1);
if opt.wtEdyn>0
    if opt.mex
        [EdynValue dEdyn]=Edyn_mex(X,Y,targetSize,stateInfo.targetsExist,length(stateVec));
    else
        if nargout>1, [EdynValue dEdyn]=Edyn(stateVec,stateInfo);
        else EdynValue=Edyn(stateVec,stateInfo);
        end
    end
    
end

% if ~isempty(find(isnan(stateVec)))
%     X
%     stateVec'
%     pause
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Eori -  orientation    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
EoriValue=0; dEori=zeros(length(stateVec),1);
if opt.wtEori>0
    if opt.mex
        [EoriValue dEori]=Eori_mex(X,Y,detMatrices.Xd,detMatrices.Yd,detMatrices.Dx,detMatrices.Dy,targetSize,stateInfo.targetsExist,length(stateVec));
    else
        if nargout>1, [EoriValue dEori]=Eori(stateVec,stateInfo);
        else EoriValue=Eori(stateVec,stateInfo);
        end
    end
    
end
% EoriValue

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Eexc - distance between objects %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
EexcValue = 0; dEexc=zeros(length(stateVec),1);

% [size(X) size(itToInd)]
% itToInd
if N>1 && opt.wtEexc>0
    if opt.mex
        if nargout>1
            [EexcValue dEexc]=Eexc_mex(X,Y,targetSize,tiToInd,length(stateVec));
        else
            EexcValue=Eexc_mex(X,Y,targetSize,tiToInd,length(stateVec));
        end
    else
        if nargout>1, [EexcValue dEexc]=Eexc(stateVec,stateInfo);
        else EexcValue=Eexc(stateVec, stateInfo);
        end
    end
    
end
%
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Eper - persistent tracks %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%
EperValue = 0; dEper=zeros(length(stateVec),1);

if opt.wtEper>0
    if opt.mex
        [EperValue dEper]=Eper_mex(X,Y,areaLimits,targetSize,targetsExist,length(stateVec));
    else
        if nargout>1, [EperValue dEper]=Eper(stateVec,stateInfo);
        else EperValue=Eper(stateVec, stateInfo);
        end
    end
end


%%%%%%%%%%%%%%%%%%
% regularization %
%%%%%%%%%%%%%%%%%%
Eregularization1=N;
% Eregularization1=(N-4)^2;
% Eregularization2=-sum(sqrt(diff(targetsExist,[],2)+1));
Eregularization2=sum(1./(diff(targetsExist,[],2)+1));
EregValue=Eregularization1+1*Eregularization2;


%%%%%%%%%%%%%%%%%%
% Exper with Ngt %
%%%%%%%%%%%%%%%%%%
% global gtInfo;
% Xgt=gtInfo.X;
% Ndev=sum(abs((sum(~~X,2)-sum(~~Xgt,2))));
% [sum(~~X,2)';sum(~~Xgt,2)']
% Ndev
% pause
% EregValue=EregValue+opt.wtEcnt*Ndev;


% final value is a linear combination of all terms
fx= opt.wtEdet*EdetValue + ...
    opt.wtEdyn*EdynValue + ...
    opt.wtEexc*EexcValue + ...
    opt.wtEper*EperValue + ...
    opt.wtEreg*EregValue + ...
    opt.wtEapp*EappValue + ...
    opt.wtEori*EoriValue;

% and the gradient
if nargout>1
    dfx = ...
        opt.wtEdet*dEdet + ...
        opt.wtEdyn*dEdyn + ...
        opt.wtEexc*dEexc + ...
        opt.wtEper*dEper + ...
        opt.wtEapp*dEapp + ...
        opt.wtEori*dEori;
    
end