function [stateInfo allpurgeenergy]=makeRemoveMove(stateInfo, doTargetN)
%% removal move
%



% get state info
[~, N, ~, targetsExist, X, Y]=getStateInfo(stateInfo);

% original state
Nold=N;
stateInfoOld=stateInfo;

allpurgeenergy=zeros(1,N);

x=stateInfo.stateVec;
origenergy=E(x,stateInfo); % current energy with removed frames

% in which order should we process?
% tarorder=randperm(Nold);
tarorder=1:Nold;

if nargin>1
     tarorder=doTargetN;
end

printMessage(2,'removing targets...');
for i=tarorder
    if i>N, break; end;
        
    printMessage(2,'.');
    te=targetsExist;
    Xt=X;Yt=Y;
    
    %what would the world look like without target i
    otheris=setdiff(1:N,i);
    
    Xt=Xt(:,otheris);Yt=Yt(:,otheris);
    targetsExist=targetsExist(otheris,:);
    
    N=N-1;
    
    stateInfo.N=N;
    stateInfo.targetsExist=targetsExist;
    stateInfo=matricesToVector(Xt,Yt,stateInfo);
    xt=stateInfo.stateVec;
    newenergy=E(xt,stateInfo); % what is the new energy
        
    allpurgeenergy(i)=newenergy-origenergy;
    if newenergy<origenergy
        printMessage(2,', %i',i);
%         fprintf('remove target %i. Energy gain: %f\n',i,newenergy-origenergy);
        X=Xt; Y=Yt;
        te=targetsExist;
        origenergy=newenergy;
        stateInfoOld=stateInfo;
        
    else
%         fprintf('removing target %i does no good\n',i);
        targetsExist=te;        
        N=N+1;
        stateInfo=stateInfoOld;
    end
    
end

stateInfo=matricesToVector(X,Y,stateInfo);
printMessage(2,'...done\n');
% [Xold X]
% pause
end