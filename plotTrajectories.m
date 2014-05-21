function plotTrajectories(stateVec,stateInfo,plot3d)
% plot trajectories
% 


hold on

% plot3d=1;
if ~exist('plot3d','var')
    plot3d=0;
end

lw=2;

N=stateInfo.N;
targetsExist=stateInfo.targetsExist;
[X Y]=vectorToMatrices(stateVec,stateInfo);

if plot3d
     for id=1:N
        tarFrames=(targetsExist(id,1):targetsExist(id,2))';
%         tarFrames=tarFrames(tarFrames<=50);
        plot3(X(tarFrames,id),Y(tarFrames,id),tarFrames, ...
            'color',getColorFromID(id),'linewidth',lw);
    end   
    
else
    for id=1:N
        tarFrames=(targetsExist(id,1):targetsExist(id,2))';
        plot(X(tarFrames,id),Y(tarFrames,id), ...
            'color',getColorFromID(id),'linewidth',lw);
    end
end
drawnow
end