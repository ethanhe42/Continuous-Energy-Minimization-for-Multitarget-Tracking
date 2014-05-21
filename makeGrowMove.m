function [stateInfo allgrowenergy]=makeGrowMove(stateInfo, doTargetN, doTargetN2)
%% extension move
% 


% get state info
[~, N F targetsExist X Y]=getStateInfo(stateInfo);

% original state
Xold=X; Yold=Y;Nold=N;
stateInfoOld=stateInfo;

global opt % we need Eper weight
global allgrow

% how do we extrapolate?
% extrapmeth='spline';
extrapmeth='linear';

% how many frames to extend max?
MAXEXT=50;

allgrowenergy=zeros(2,Nold);
allgrow=zeros(2,Nold);

% params=getOptimParameters();

% in which order should we process?
% tarorderbck=randperm(Nold);
% tarorderfrw=randperm(Nold);
tarorderbck=1:Nold;
tarorderfrw=1:Nold;

if nargin==2
    tarorderbck=doTargetN;
    tarorderfrw=doTargetN;
elseif nargin==3
    tarorderbck=doTargetN;
    tarorderfrw=doTargetN2;
end

% save info about growth of each trajectory
targrownbck=zeros(1,Nold);
targrownfrw=zeros(1,Nold);

if nargin==1
    printMessage(2,'Target:  ')
    printMessage(2,'%4i',tarorderbck);
    printMessage(2,' Tars Frms');
    printMessage(2,'\ngrow bck:');
end

for id=tarorderbck
    %     printMessage(2,'.');
    te=targetsExist;
    X=Xold; Y=Yold;
    
    
    existframes=targetsExist(id,1):targetsExist(id,2);
    if existframes(1) > 1
        %         plot3(X(targetsExist(id,1):targetsExist(id,2),id),Y(targetsExist(id,1):targetsExist(id,2),id),targetsExist(id,1):targetsExist(id,2),'color',getColorFromID(id));
        %         pause
        % Remove targets with no influence
        notimportanti=find(targetsExist(:,1)>=existframes(1));        % all that start after current begins
        notimportanti=union(notimportanti,find(targetsExist(:,2)<existframes(1)-MAXEXT));       % all that end MAXEXT before current begins
        
        newid=id-length(find(notimportanti<id)); % new index
        
        notimportanti=setdiff(notimportanti,id); % make sure to keep current
        importanti=setdiff(1:Nold,notimportanti); % all important ones
        
        targetsExist=targetsExist(importanti,:); % keep important
        N=size(targetsExist,1); % new number of targets
        targetsExist(:,1)=max(targetsExist(:,1),existframes(1)-MAXEXT); % shorten from past
        targetsExist(:,2)=min(targetsExist(:,2),existframes(1)+2); % shorten future
        
        
        targetsExist(:,1)=targetsExist(:,2)-max(2,diff(targetsExist,[],2)); % take care of special case
        
        X=X(:,importanti);Y=Y(:,importanti); % adjust X and Y
        
        for ii=1:N % fill rest of frames with zeros
            X(1:targetsExist(ii,1)-1,ii)=0; Y(1:targetsExist(ii,1)-1,ii)=0;
            X(targetsExist(ii,2)+1:end,ii)=0; Y(targetsExist(ii,2)+1:end,ii)=0;
        end
        
        % what is the old energy
        stateInfo.N=N; stateInfo.targetsExist=targetsExist;
        stateInfo=matricesToVector(X,Y,stateInfo);
        x=stateInfo.stateVec;
        origenergy=E(x,stateInfo); % current energy with removed frames

        %% consider cut off ends for Ereg2
        origenergy=origenergy ...
            -opt.wtEper*sum(1./(diff(targetsExist,[],2)+1)) ...
            +opt.wtEper*sum(1./(diff(te,[],2)+1));
        
        Xt=X;Yt=Y;
        
        Xi=Xt(:,newid); Yi=Yt(:,newid);
        fromframe=existframes(1);
        
        minenergy=origenergy;
        bestmove=[];
        
        possible_moves=fromframe-1:-1:max(1,fromframe-MAXEXT); % store individual energy gains for debugging
        moveenergy=zeros(1,length(possible_moves)); % store individual energy gains for debugging
        
        for ft=fromframe-1:-1:max(1,fromframe-MAXEXT) % find best extension for current target
            %                             printMessage(2,'extrapolate %id frames backwards\n',fromframe-ft);
%             Xi(ft:fromframe-1)=interp1(existframes,Xi(existframes),ft:fromframe-1,extrapmeth,'extrap');
            Xi(ft:fromframe-1)=myLinExtrap(Xi(existframes(1:2)),ft-fromframe)';
            
%             Yi(ft:fromframe-1)=interp1(existframes,Yi(existframes),ft:fromframe-1,extrapmeth,'extrap');
            Yi(ft:fromframe-1)=myLinExtrap(Yi(existframes(1:2)),ft-fromframe)';
            
            Xt(:,newid)=Xi; Yt(:,newid)=Yi;
            
            targetsExist(newid,1)=targetsExist(newid,1)-1; % adjust
            te2=te; te2(id,1)=te2(id,1)-1;

            stateInfo.targetsExist=targetsExist;
            stateInfo=matricesToVector(Xt,Yt,stateInfo);
            xt=stateInfo.stateVec;
            newenergy=E(xt,stateInfo); % what is the new energy
            
            %% consider cut off ends for Ereg2 (the unimportant ones)
            newenergy=newenergy ...
                - opt.wtEper*sum(1./(diff(targetsExist,[],2)+1)) ...
                + opt.wtEper*sum(1./(diff(te2,[],2)+1));
            moveenergy(fromframe-ft) = newenergy;
            
            if newenergy<minenergy
                minenergy=newenergy;
                bestmove=ft:fromframe-1;
            end
            %             plot3(Xt(targetsExist(id,1):targetsExist(id,2),id),Yt(targetsExist(id,1):targetsExist(id,2),id),targetsExist(id,1):targetsExist(id,2),'color',getColorFromID(id));
            %             pause
        end
        
        % restore old
        X=Xold; Y=Yold;
        N=Nold;
        targetsExist=te;
        stateInfo=stateInfoOld;
        
        
        allgrowenergy(1,id)=minenergy-origenergy;
        
        if ~isempty(bestmove) %if energy gain, extend
            %             printMessage(2,', +%if<-%id',bestmove(end)+1-bestmove(1),id);
            
            targrownbck(id)=bestmove(end)+1-bestmove(1);
            printMessage(2,'%4i',targrownbck(id))
            
            X(bestmove,id)=interp1(existframes,Xi(existframes),bestmove,extrapmeth,'extrap');
            Y(bestmove,id)=interp1(existframes,Yi(existframes),bestmove,extrapmeth,'extrap');
            targetsExist(id,1)=bestmove(1);
            targetsExist(id,2)=te(id,2);
            
            Xold=X; Yold=Y; % replace X Y
            
            stateInfo.targetsExist=targetsExist;
            stateInfo=matricesToVector(X,Y,stateInfo);
            stateInfoOld=stateInfo;
            
            xt=stateInfo.stateVec;
           

        else
            printMessage(2,'   -');
        end
    else
        printMessage(2,'   -');
    end
end
printMessage(2,'%5i%5i',numel(find(targrownbck)),sum(targrownbck));

if nargin==1 % if no order given, use standard
    if ~all(size(tarorderbck)==size(tarorderfrw)) || (all(size(tarorderbck)==size(tarorderfrw)) && ~all(tarorderbck==tarorderfrw))
        printMessage(2,'\n\nTarget:  ')
        printMessage(2,'%4i',tarorderfrw)
        printMessage(2,' Tars Frms');   
    end
    printMessage(2,'\ngrow frw:');
end

% % forwards
for id=tarorderfrw
    %     printMessage(2,'.');
    te=targetsExist;
    X=Xold; Y=Yold;
%     
    existframes=targetsExist(id,1):targetsExist(id,2);
    if existframes(end) < F
        
        % Remove targets with no influence
        notimportanti=find(targetsExist(:,2)<=existframes(end));        % all that end before current ends
        notimportanti=union(notimportanti,find(targetsExist(:,1)>existframes(end)+MAXEXT));       % all that end MAXEXT before current begins
        
        newid=id-length(find(notimportanti<id)); % new index
        
        notimportanti=setdiff(notimportanti,id); % make sure to keep current
        importanti=setdiff(1:Nold,notimportanti); % all important ones
        
        targetsExist=targetsExist(importanti,:); % keep important
        N=size(targetsExist,1); % new number of targets
        
        targetsExist(:,1)=max(targetsExist(:,1),existframes(end)-2); % shorten from past
        targetsExist(:,2)=min(targetsExist(:,2),existframes(end)+MAXEXT); % shorten future
        
        
        targetsExist(:,2)=targetsExist(:,1)+max(2,diff(targetsExist,[],2)); % take care of special case
        
        X=X(:,importanti);Y=Y(:,importanti); % adjust X and Y
        
        
        
        for ii=1:N % fill rest of frames with zeros
            X(1:targetsExist(ii,1)-1,ii)=0; Y(1:targetsExist(ii,1)-1,ii)=0;
            X(targetsExist(ii,2)+1:end,ii)=0; Y(targetsExist(ii,2)+1:end,ii)=0;
        end
        
        % what is the old energy
        stateInfo.N=N; stateInfo.targetsExist=targetsExist;
        stateInfo=matricesToVector(X,Y,stateInfo);
        x=stateInfo.stateVec;
        origenergy=E(x,stateInfo); % current energy with removed frames
        
        %% consider cut off ends for Ereg2
        origenergy=origenergy ...
            -opt.wtEper*sum(1./(diff(targetsExist,[],2)+1)) ...
            +opt.wtEper*sum(1./(diff(te,[],2)+1));
        
        Xt=X;Yt=Y;
        
        
        
        Xi=Xt(:,newid); Yi=Yt(:,newid);
        toframe=existframes(end);
        
        minenergy=origenergy;
        bestmove=[];
        
        possible_moves=toframe+1:min(F,toframe+MAXEXT); % store individual energy gains for debugging
        moveenergy=zeros(1,length(possible_moves)); % store individual energy gains for debugging
        
        for ft=toframe+1:min(F,toframe+MAXEXT) % find best extension for current target
            %             printMessage(2,'extrapolate %id frames forward\n',ft-toframe);
%             Xi(toframe+1:ft)=interp1(existframes,Xi(existframes),toframe+1:ft,extrapmeth,'extrap');
            Xi(toframe+1:ft)=myLinExtrap(Xi(existframes(end-1:end)),ft-toframe)';
            
%             Yi(toframe+1:ft)=interp1(existframes,Yi(existframes),toframe+1:ft,extrapmeth,'extrap');
            Yi(toframe+1:ft)=myLinExtrap(Yi(existframes(end-1:end)),ft-toframe)';
            Xt(:,newid)=Xi; Yt(:,newid)=Yi;
            
            targetsExist(newid,2)=targetsExist(newid,2)+1;
            te2=te; te2(id,2)=te2(id,2)+1;
            
            stateInfo.targetsExist=targetsExist;
            stateInfo=matricesToVector(Xt,Yt,stateInfo);
            xt=stateInfo.stateVec;
            newenergy=E(xt,stateInfo); % what is the new energy
            
            %% consider cut off ends for Ereg2 (the unimportant ones)
            newenergy=newenergy ...
                - opt.wtEper*sum(1./(diff(targetsExist,[],2)+1)) ...
                + opt.wtEper*sum(1./(diff(te2,[],2)+1));
            moveenergy(ft-toframe) = newenergy;
                       
            if newenergy<minenergy
                minenergy=newenergy;
                bestmove=toframe+1:ft;
            end
        end
        
        % restore old
        X=Xold; Y=Yold;
        N=Nold;
        targetsExist=te;
        stateInfo=stateInfoOld;
        
        allgrowenergy(2,id)=minenergy-origenergy;

        
%         if newenergy<-1000000
%             printMessage(2,'something wrong. bg: %f, best: %i %i %i\n',bestgain, best);

        if ~isempty(bestmove) %if energy gain, extend
            %             printMessage(2,'extrapolate best for target %id: %id frames forwards, Energy gain: %f\n',id,bestmove(end)-bestmove(1)+1,minenergy-origenergy);
            %             printMessage(2,', %id->+%if',id,bestmove(end)-bestmove(1)+1);
            
            targrownfrw(id)=bestmove(end)-bestmove(1)+1;
            printMessage(2,'%4i',targrownfrw(id))
            
            X(bestmove,id)=interp1(existframes,Xi(existframes),bestmove,extrapmeth,'extrap');
            Y(bestmove,id)=interp1(existframes,Yi(existframes),bestmove,extrapmeth,'extrap');
            targetsExist(id,2)=bestmove(end);
            targetsExist(id,1)=te(id,1);
            Xold=X; Yold=Y;
            stateInfo.targetsExist=targetsExist;
            stateInfo=matricesToVector(X,Y,stateInfo);
            stateInfoOld=stateInfo;
            
            xt=stateInfo.stateVec;
            
        else
            printMessage(2,'   -');
        end
    else
        printMessage(2,'   -');
    end
end
printMessage(2,'%5i%5i',numel(find(targrownfrw)),sum(targrownfrw));
allgrow(1,:)=targrownbck;
allgrow(2,:)=targrownfrw;



stateInfo=matricesToVector(X,Y,stateInfo);

printMessage(2,'...done\n');
% printMessage(2,'Trajectory length');
% printMessage(2,'%4i',diff(targetsExist,[],2)')
% [Xold X]
% pause
% toc;
end