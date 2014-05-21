%% shrink move
function [stateInfo allshrinkenergy]=makeShrinkMove(stateInfo,doTargetN,doTargetN2)
% 


% get state info
[~, N F targetsExist X Y]=getStateInfo(stateInfo);

% original state
Xold=X; Yold=Y;Nold=N;
stateInfoOld=stateInfo;

global opt % we need Eper weight
global allshr

MAXSHR=50;
Nold=N;

allshrinkenergy=zeros(2,Nold);
allshr=zeros(2,Nold);

% in which order should we process?
tarorderbck=randperm(Nold);
tarorderfrw=randperm(Nold);
tarorderbck=1:Nold;
tarorderfrw=1:Nold;

if nargin==2
    tarorderbck=doTargetN;
    tarorderfrw=doTargetN;
elseif nargin==3
    tarorderbck=doTargetN;
    tarorderfrw=doTargetN2;
end

tarshrnkbck=zeros(1,Nold);
tarshrnkfrw=zeros(1,Nold);

if nargin==1
printMessage(2,'Target:  ')
printMessage(2,'%4i',tarorderbck);
printMessage(2,' Tars Frms');
printMessage(2,'\nshrk pst:');
end
% from past
for id=tarorderbck
%     printMessage(2,'.');
    te=targetsExist;
    X=Xold; Y=Yold;
    
    existframes=targetsExist(id,1):targetsExist(id,2);
    
    % Remove targets with no influence
    notimportanti=find(targetsExist(:,1)>existframes(1)+MAXSHR);        % all that start MAXSHR after current begins
    notimportanti=union(notimportanti,find(targetsExist(:,2)<existframes(1)));       % all that end before current begins
    
    newid=id-length(find(notimportanti<id)); % new index
    %notimportanti=setdiff(notimportanti,id); % make sure to keep current
    importanti=setdiff(1:Nold,notimportanti); % all important ones
    
    targetsExist=targetsExist(importanti,:); % keep important
    N=size(targetsExist,1); % new number of targets
    targetsExist(:,1)=max(targetsExist(:,1),existframes(1)); % shorten from past
    targetsExist(:,1)=targetsExist(:,2)-max(2,diff(targetsExist,[],2)); % take care of special case
    
    targetsExist(:,2)=min(targetsExist(:,2),existframes(1)+MAXSHR+2); % shorten future
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
    existframes=targetsExist(newid,1):targetsExist(newid,2);
    fromframe=existframes(1);
    toframe=existframes(end);
    
    
    minenergy=origenergy;
    bestmove=[];
    
    possible_moves=fromframe:min(fromframe+MAXSHR-1,toframe-3); % store individual energy gains for debugging
    moveenergy=zeros(1,length(possible_moves)); % store individual energy gains for debugging
    
    for ft=fromframe:min(fromframe+MAXSHR-1,toframe-3) % dont shrink to less than 3 frames
        %             printMessage(2,'remove %id frames from past\n',ft-fromframe+1);
        Xi(fromframe:ft)=0;
        Yi(fromframe:ft)=0;
        Xt(:,newid)=Xi; Yt(:,newid)=Yi;
        
        targetsExist(newid,1)=targetsExist(newid,1)+1;
            
        stateInfo.targetsExist=targetsExist;
        stateInfo=matricesToVector(Xt,Yt,stateInfo);
        xt=stateInfo.stateVec;
        newenergy=E(xt,stateInfo); % what is the new energy
            
        te2=te; te2(id,2)=te2(id,1)+1;
        
        newenergy=newenergy ...
            - opt.wtEper*sum(1./(diff(targetsExist,[],2)+1)) ...
            + opt.wtEper*sum(1./(diff(te2,[],2)+1));
        
        moveenergy(ft-fromframe+1) = newenergy;

        if newenergy<minenergy
            minenergy=newenergy;
            bestmove=fromframe:ft;
        end
    end
    
    % restore old
    X=Xold; Y=Yold;
    N=Nold;
    targetsExist=te;
    stateInfo=stateInfoOld;
    
    allshrinkenergy(1,id)=minenergy-origenergy;

    
    if ~isempty(bestmove)
%         printMessage(2,'remove %id frames from past from target %id\n',length(bestmove),id);
%         printMessage(2,', -%if<-%id',length(bestmove),id);
        printMessage(2,'%4i',length(bestmove))
        tarshrnkbck(id)=length(bestmove);
        allshr(1,id)=length(bestmove);

        X(bestmove,id)=0;
        Y(bestmove,id)=0;
        targetsExist(id,1)=bestmove(end)+1;
        Xold=X; Yold=Y; % replace X Y
        stateInfo.targetsExist=targetsExist;
        stateInfo=matricesToVector(X,Y,stateInfo);
        stateInfoOld=stateInfo;
        
        xt=stateInfo.stateVec;
        
    else
            printMessage(2,'   -');
    end
    
end
if nargin==1
printMessage(2,'%5i%5i',numel(find(tarshrnkbck)),sum(tarshrnkbck));

if ~all(size(tarorderbck)==size(tarorderfrw)) || (all(size(tarorderbck)==size(tarorderfrw)) && ~all(tarorderbck==tarorderfrw))
    printMessage(2,'\n\nTarget:  ')
    printMessage(2,'%4i',tarorderfrw)
    printMessage(2,' Tars Frms');
end

printMessage(2,'\nshrk ftr:');
end
for id=tarorderfrw
%     printMessage(2,'.');
    te=targetsExist;
    X=Xold; Y=Yold;
    
    existframes=targetsExist(id,1):targetsExist(id,2);
    
    % Remove targets with no influence
    notimportanti=find(targetsExist(:,1)>existframes(end));        % all that start after current ends
    notimportanti=union(notimportanti,find(targetsExist(:,2)<=existframes(end)-MAXSHR));       % all that end MAXSHR before current begins
    
    newid=id-length(find(notimportanti<id)); % new index
    notimportanti=setdiff(notimportanti,id); % make sure to keep current
    importanti=setdiff(1:Nold,notimportanti); % all important ones
    
    targetsExist=targetsExist(importanti,:); % keep important
    N=size(targetsExist,1); % new number of targets
    targetsExist(:,1)=max(targetsExist(:,1),existframes(end)-MAXSHR-2); % shorten from past
    targetsExist(:,1)=targetsExist(:,2)-max(2,diff(targetsExist,[],2)); % take care of special case
    
    targetsExist(:,2)=min(targetsExist(:,2),existframes(end)); % shorten future
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
    existframes=targetsExist(newid,1):targetsExist(newid,2);
    fromframe=existframes(1);
    toframe=existframes(end);
    
    
    minenergy=origenergy;
    bestmove=[];
    
    for ft=toframe:-1:max(toframe-MAXSHR+1,fromframe+3) % dont shrink to less than 3 frames
        %             printMessage(2,'remove %id frames from past\n',ft-fromframe+1);
        Xi(ft:toframe)=0;
        Yi(ft:toframe)=0;
        Xt(:,newid)=Xi; Yt(:,newid)=Yi;
        
        targetsExist(newid,2)=targetsExist(newid,2)-1;
        
        stateInfo.targetsExist=targetsExist;
        stateInfo=matricesToVector(Xt,Yt,stateInfo);
        xt=stateInfo.stateVec;
        newenergy=E(xt,stateInfo); % what is the new energy

        te2=te; te2(id,2)=te2(id,2)-1;
        %% consider cut off ends for Ereg2 (the unimportant ones)
        newenergy=newenergy ...
            - opt.wtEper*sum(1./(diff(targetsExist,[],2)+1)) ...
            + opt.wtEper*sum(1./(diff(te2,[],2)+1));
        moveenergy(ft-fromframe+1) = newenergy;
        
        if newenergy<minenergy
            minenergy=newenergy;
            bestmove=ft:toframe;
        end
    end
    % restore old
    X=Xold; Y=Yold;
    N=Nold;
    targetsExist=te;
    stateInfo=stateInfoOld;
    
    allshrinkenergy(2,id)=minenergy-origenergy;
    
    if ~isempty(bestmove)
%         printMessage(2,'remove %id frames from future from target %id\n',length(bestmove),id);
        printMessage(2,'%4i',length(bestmove))
        tarshrnkfrw(id)=length(bestmove);
        allshr(2,id)=length(bestmove);

%         printMessage(2,', %id->-%if',length(bestmove),id);
        X(bestmove,id)=0;
        Y(bestmove,id)=0;
        targetsExist(id,2)=bestmove(1)-1;
        Xold=X; Yold=Y; % replace X Y
        stateInfo.targetsExist=targetsExist;
        stateInfo=matricesToVector(X,Y,stateInfo);
        stateInfoOld=stateInfo;
        
        xt=stateInfo.stateVec;
    else
            printMessage(2,'   -');
    end
end
printMessage(2,'%5i%5i',numel(find(tarshrnkfrw)),sum(tarshrnkfrw));

stateInfo=matricesToVector(X,Y,stateInfo);

printMessage(2,'...done\n');
% toc
% [Xold X]
% pause
end