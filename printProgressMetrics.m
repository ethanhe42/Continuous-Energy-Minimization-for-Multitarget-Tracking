function printProgressMetrics(opt, sceneInfo, X, stateInfo, i)
global allmets3d allmets2d globiter gtInfo

if sceneInfo.gtAvailable
    
    
    stateInfo.stateVec=X;
    [~,~,~,~, stateInfo.X stateInfo.Y]=getStateInfo(stateInfo);
    
    % set state struct
    if opt.track3d
        stateInfo=cutStateToTrackingArea(stateInfo);
    end
    %             [Xm Ym]=vectorToMatrices(X,stateInfo);
    %             stateInfo.X=Xm; stateInfo.Y=Ym;
    if opt.track3d
        stateInfo.Xgp=stateInfo.X; stateInfo.Ygp=stateInfo.Y;
        [stateInfo.Xi stateInfo.Yi]=projectToImage(stateInfo.Xgp,stateInfo.Ygp,sceneInfo);
    else
        stateInfo.Xi=stateInfo.X; stateInfo.Yi=stateInfo.Y;
    end
    
    
    % evaluate 3D CLEAR MOT
    if opt.track3d
        evopt.eval3d=1;
        [metrics3d metricsInfo3d]=CLEAR_MOT(gtInfo,stateInfo,evopt);
        printMetrics(metrics3d,metricsInfo3d,0,[12 13 4 5 7:11 1 2]);
        allmets3d(globiter,:)=metrics3d;
    end
    
%     printMessage(3,'|||');
    
     if ~mod(globiter,20) 
        % evaluate 2D CLEAR MOT
        % get bounding boxes from corresponding detections
        stateInfo=getBBoxesFromState(stateInfo);
        
        evopt.eval3d=0;
        [metrics2d metricsInfo2d]=CLEAR_MOT(gtInfo,stateInfo,evopt);
        printMetrics(metrics2d,metricsInfo2d,0,[12 13 4 5 7:11 1 2]);
        allmets2d(globiter,:)=metrics2d;
    end
    
end
end