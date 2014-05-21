function [stateInfo allsplitenergy]=makeSplitMove(stateInfo, doTargetN)
% checks for each location, whether a split of
% the trajectory yields lower energy
% 


% get state info
[~, N F targetsExist X Y]=getStateInfo(stateInfo);

% original state
Xold=X; Yold=Y;Nold=N;
stateInfoOld=stateInfo;

global opt % we need Eper weight

allsplitenergy=zeros(1,N);

Nold=N; % remember original number of targets

%tarorder=randperm(Nold);
tarorder=1:Nold;

if nargin>1
     tarorder=doTargetN;
end

printMessage(2,'splitting...');

% for each target...
for i=tarorder
    printMessage(2,'.');
    
    % copy current configuration
    te=targetsExist;
    X=Xold; Y=Yold;
    
    % where in time does this target exist?
    existframes=targetsExist(i,1):targetsExist(i,2);
    fromframe=existframes(1);
    toframe=existframes(end);
    
    % if length < 6 then no split possible
    if length(existframes)>=6
        
        N=1;                        % global N, needed for evaluation
        X=X(:,i); Y=Y(:,i);
        targetsExist=targetsExist(i,:);
        
        % what is the energy of this target?
        stateInfo.N=N; stateInfo.targetsExist=targetsExist;
        stateInfo=matricesToVector(X,Y,stateInfo);
        x=stateInfo.stateVec;
        origenergy=E(x,stateInfo);
        
        
        minenergy=origenergy;
        bestmove=0;
        
        
        N=2;                        % global N, needed for evaluation
        
        % now check for all frames, where a split makes most sense
        for ft=fromframe+3:toframe-2 % ft = new targets' start frame
            
            % shift part of first column of X,Y to the second
            Xt=X; Yt=Y;
            Xt(ft:toframe,2)=X(ft:toframe,1);
            Yt(ft:toframe,2)=Y(ft:toframe,1);
            Xt(ft:end,1)=0;Yt(ft:end,1)=0;
            
            
            targetsExist(1,2)=ft-1;
            targetsExist(2,1:2)=[ft toframe];
            
            % what is the new energy?
            stateInfo.targetsExist=targetsExist;
            stateInfo=matricesToVector(Xt,Yt,stateInfo);
            xt=stateInfo.stateVec;
            newenergy=E(xt,stateInfo); % what is the new energy
            
            % if best, save
            if newenergy<minenergy
                minenergy=newenergy;
                bestmove=ft;           % new target starts here
            end
        end
        
        % restore old
        X=Xold; Y=Yold;
        N=Nold;
        targetsExist=te;
        stateInfo=stateInfoOld;
        
        allsplitenergy(i)=minenergy-origenergy;
        
        % if split makes sense... do so
        if minenergy-origenergy<0
            N=N+1;                  % add new target
            printMessage(3,', %i (%i)',i,bestmove);
            %             printMessage(2,'split target %i at frame %i. Energy gain: %f \n',i,bestmove,minenergy-origenergy);
            X(bestmove:toframe,N)=X(bestmove:toframe,i);
            Y(bestmove:toframe,N)=Y(bestmove:toframe,i);
            X(bestmove:end,i)=0;Y(bestmove:end,i)=0;
            
            targetsExist(i,2)=bestmove-1;
            targetsExist(N,1:2)=[bestmove toframe];
            
            Xold=X; Yold=Y; % replace X Y
            Nold=N;
            stateInfo.N=N;
            stateInfo.targetsExist=targetsExist;
            stateInfo=matricesToVector(X,Y,stateInfo);
            stateInfoOld=stateInfo;
            
            xt=stateInfo.stateVec;
            
        end
    end
    
end

% generate new state vector
stateInfo=matricesToVector(X,Y,stateInfo);
printMessage(2,'...done\n');

end
