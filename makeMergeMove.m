function [stateInfo allmergeenergy]=makeMergeMove(stateInfo)
% merge two trajectories if possible
% 


% tic;

% how do we interpolate?
% interpmeth='spline';
interpmeth='linear';

global sceneInfo;

stateVecOld=stateInfo.stateVec;
N=stateInfo.N;F=stateInfo.F;
targetsExist=stateInfo.targetsExist;
gridStep=sceneInfo.targetSize;


[X Y]=vectorToMatrices(stateVecOld,stateInfo);

x=stateVecOld;
origenergy=E(x,stateInfo);


allmergeenergy=zeros(N);


printMessage(2,'merging...');
itcnt=0;
while itcnt<1000 && N
%     fprintf('.');
    itcnt=itcnt+1;
    te=targetsExist;
    stateInfoOld=stateInfo;

    mergeCost=Inf*ones(N);
    
    %% look for potential connections first (speedlimit)
    potentialConnections=false(N,N);
    for id=1:N
        otheri=setdiff(1:N,id);
        for id2=otheri
            tend=targetsExist(id,2);
            tsta=targetsExist(id2,1);
            if tend<tsta
                framesbetween=tsta-tend;
                distance=norm([X(tend,id) Y(tend,id)]-[X(tsta,id2) Y(tsta,id2)]);
                speed=distance/framesbetween;
                if speed<3*gridStep && distance<5*gridStep && framesbetween<15 % avoid unplausible connections
%                     fprintf('%id and %id can be connected, distance: %f, speed p/f: %f\n',id,id2,distance,speed);
                    potentialConnections(id,id2)=1;
                end
            end
        end
    end
    %% now look which connection is best
    potidx=find(potentialConnections);
%     fprintf('%id c:',numel(potidx));
    if ~isempty(potidx)

        for pi=potidx'
            [id id2]=ind2sub([N N],pi);

                Xt=X; Yt=Y;
                fromframe1=targetsExist(id,1);
                toframe1=targetsExist(id,2);
                fromframe2=targetsExist(id2,1);
                toframe2=targetsExist(id2,2);
                
                newframes=[fromframe1:toframe1 fromframe2:toframe2];
                
                onecolX=[Xt(fromframe1:toframe1,id);Xt(fromframe2:toframe2,id2)];
                onecolY=[Yt(fromframe1:toframe1,id);Yt(fromframe2:toframe2,id2)];
                
                onecolX=interp1(newframes,onecolX,newframes(1):newframes(end),interpmeth);
                
                onecolY=interp1(newframes,onecolY,newframes(1):newframes(end),interpmeth);
                onecolX=onecolX'; onecolY=onecolY';
                
                Xt(fromframe1:toframe2,id)=onecolX;Yt(fromframe1:toframe2,id)=onecolY;
                
                % try to smooth
                sw=[0.2 0.3 0.3 0.2]; % smoothing weights
%                 sw=[0 1 0 0]; % smoothing weights
                Xt(toframe1+1,id)=sw(1)*Xt(toframe1-1,id)+sw(2)*Xt(toframe1,id)+sw(3)*Xt(toframe1+2,id)+sw(4)*Xt(toframe1+3,id);
                Yt(toframe1+1,id)=sw(1)*Yt(toframe1-1,id)+sw(2)*Yt(toframe1,id)+sw(3)*Yt(toframe1+2,id)+sw(4)*Yt(toframe1+3,id);
                Xt(fromframe2-1,id)=sw(1)*Xt(fromframe2-3,id)+sw(2)*Xt(fromframe2-2,id)+sw(3)*Xt(fromframe2,id)+sw(4)*Xt(fromframe2+1,id);
                Yt(fromframe2-1,id)=sw(1)*Yt(fromframe2-3,id)+sw(2)*Yt(fromframe2-2,id)+sw(3)*Yt(fromframe2,id)+sw(4)*Yt(fromframe2+1,id);
                
 

%                 Yt(toframe1+1,id)=0.5*(Yt(toframe1,id)+Yt(toframe1+2,id));
%                 Xt(fromframe2-1,id)=0.5*(Xt(fromframe2-2,id)+Xt(fromframe2,id));
%                 Yt(fromframe2-1,id)=0.5*(Yt(fromframe2-2,id)+Yt(fromframe2,id));

                otheri=setdiff(1:N,id2);
                Xt=Xt(:,otheri);Yt=Yt(:,otheri);
                
                targetsExist(id,1)=fromframe1;targetsExist(id,2)=toframe2;
                targetsExist=targetsExist(otheri,:);
                
                N=N-1;
                
                stateInfo.N=N; stateInfo.targetsExist=targetsExist;
                stateInfo=matricesToVector(Xt,Yt,stateInfo);
                
                xt=stateInfo.stateVec;
                newenergy=E(xt,stateInfo);
                
                mergeCost(id,id2)=newenergy-origenergy;
                
                N=N+1;
                targetsExist=te;
                
                stateInfo.N=N; stateInfo.targetsExist=targetsExist;
        end
    end
    [c mi]=min(mergeCost(:));
    if itcnt==1
        allmergeenergy=mergeCost;
    end
%     mergeCost

    if c>=0
%         fprintf('best %id and %id. Energy gain: %f\n',id,id2,c);
        printMessage(2,'...done\n');
        break;
    else
%         c
        
        [id id2]=ind2sub([N N],mi);
        if id==id2
            printMessage(2,'...done\n');
            break;
        end
        
        fromframe1=targetsExist(id,1);
        toframe1=targetsExist(id,2);
        fromframe2=targetsExist(id2,1);
        toframe2=targetsExist(id2,2);
        
        
        newframes=[fromframe1:toframe1 fromframe2:toframe2];
        printMessage(3,' %i+%i (%i fr)',id,id2,fromframe2-toframe1-1);
%         fprintf('merging %id and %id over %id frames. Energy gain: %f\n',id,id2,fromframe2-toframe1-1,c);
        onecolX=[X(fromframe1:toframe1,id);X(fromframe2:toframe2,id2)];
        onecolY=[Y(fromframe1:toframe1,id);Y(fromframe2:toframe2,id2)];
        
        onecolX=interp1(newframes,onecolX,newframes(1):newframes(end),interpmeth);
        onecolY=interp1(newframes,onecolY,newframes(1):newframes(end),interpmeth);
        onecolX=onecolX'; onecolY=onecolY';
        
        
        
%         X(fromframe1:toframe2,id)=onecolX;Y(fromframe1:toframe2,id)=onecolY;
%         otheri=setdiff(1:N,id2);
%         X=X(:,otheri);Y=Y(:,otheri);
%         
%         targetsExist(id,1)=fromframe1;targetsExist(id,2)=toframe2;
%         targetsExist=targetsExist(otheri,:);
%         
%         N=N-1;
        
        %%% SMOOTH HERE AS WELL

                 X(fromframe1:toframe2,id)=onecolX;Y(fromframe1:toframe2,id)=onecolY;

                 % try to smooth
                sw=[0.2 0.3 0.3 0.2]; % smoothing weights
%                 sw=[0 1 0 0]; % smoothing weights
                X(toframe1+1,id)=sw(1)*X(toframe1-1,id)+sw(2)*X(toframe1,id)+sw(3)*X(toframe1+2,id)+sw(4)*X(toframe1+3,id);
                Y(toframe1+1,id)=sw(1)*Y(toframe1-1,id)+sw(2)*Y(toframe1,id)+sw(3)*Y(toframe1+2,id)+sw(4)*Y(toframe1+3,id);
                X(fromframe2-1,id)=sw(1)*X(fromframe2-3,id)+sw(2)*X(fromframe2-2,id)+sw(3)*X(fromframe2,id)+sw(4)*X(fromframe2+1,id);
                Y(fromframe2-1,id)=sw(1)*Y(fromframe2-3,id)+sw(2)*Y(fromframe2-2,id)+sw(3)*Y(fromframe2,id)+sw(4)*Y(fromframe2+1,id);
                
                otheri=setdiff(1:N,id2);
                X=X(:,otheri);Y=Y(:,otheri);
                
                targetsExist(id,1)=fromframe1;targetsExist(id,2)=toframe2;
                targetsExist=targetsExist(otheri,:);
                
                N=N-1;
                stateInfo.N=N; stateInfo.targetsExist=targetsExist;
  
        
    end
end
stateInfo=matricesToVector(X,Y,stateInfo);

% toc;
end
