function [stateInfo alladdenergy]=makeAddMove(stateInfo,createN)
%% add move
% add new targets at unused detections
% 


global detMatrices;
global sceneInfo;

% get state info
[~, N F targetsExist X Y]=getStateInfo(stateInfo);

% global scene info and detections
gridStep=sceneInfo.targetSize;
Xd=detMatrices.Xd; Yd=detMatrices.Yd; Sd=detMatrices.Sd;

% original state
Xold=X; Yold=Y;Nold=N;

addnoise=0*gridStep/10;


MAXNEW=8;   % maximum number of targets to be added
printMessage(2,'creating targets...');

if exist('createN','var')
    MAXNEW=createN;
end
alladdenergy=zeros(1,MAXNEW);

for j=1:MAXNEW
    printMessage(2,'.');        
    
    best=[];
    cntbest=0;
    bestgain=Inf;
    
    % in each frame...
    for t=1:F
        
        % memorize current configuration
        X=Xold; Y=Yold;
        stateInfoOld=stateInfo;
        te=targetsExist;
        
        % first of all, find out, which detections are still 'free'
        exdet=find(Xd(t,:));
        potdet=~~Xd(t,:);
        % get rid of the ones occupied
        extar=find(X(t,:));
        for d=exdet
            for i=extar
                if norm([X(t,i) Y(t,i)]-[Xd(t,d) Yd(t,d)])<gridStep*2
%                     exdet=setdiff(exdet,d);
                    potdet(d)=0;
                end
                
            end
        end
        exdet=find(potdet);
        
        
        % sort 
        [s ix]=sort(Sd(t,exdet),'descend');
        exdet=exdet(ix);
%         potdet
%         Xd(t,exdet)
%         isequal(exdet,potdet)
%         exdet
%         potdet
%         pause
        
                
        
        if ~isempty(exdet)
            
            % new target extends over 3 frames: (t-1, t, t+1)
            newframes=t-1:t+1;
            if t==1, newframes=t:t+2; end
            if t==F, newframes=t-2:t; end
            cf=newframes(2);                                           % center frame
            
            % for speed up, remove all targets before t-1 and after t+1
            % because they are irrelevant for energy computation
            notimportanti=find(targetsExist(:,2)<newframes(1));        % all that end before newframes
            notimportanti=union(notimportanti,find(targetsExist(:,1)>newframes(3)));       % all that start after newframes
%             importanti=setdiff(1:Nold,notimportanti);                  % all important ones
            allOld=true(1,Nold);
            allOld(notimportanti)=0;
            importanti=find(allOld);
%             importanti
%             importanti2
%             pause
            
            targetsExist=targetsExist(importanti,:);                    % keep important
            targetsExist(:,1)=max(targetsExist(:,1),newframes(1));      % shorten from past
            
            targetsExist(:,1)=targetsExist(:,2)-max(2,diff(targetsExist,[],2)); % take care of special case
            
            targetsExist(:,2)=min(targetsExist(:,2),newframes(2));      % shorten future
            targetsExist(:,2)=targetsExist(:,1)+max(2,diff(targetsExist,[],2)); % take care of special case
            
            Xt=X(:,importanti);Yt=Y(:,importanti);
            
            N=size(targetsExist,1);
            %             printMessage(2,'oldc %i\n',size(targetsExist,1));
            
            for ii=1:N % fill rest of frames with zeros
                Xt(1:targetsExist(ii,1)-1,ii)=0; Yt(1:targetsExist(ii,1)-1,ii)=0;
                Xt(targetsExist(ii,2)+1:end,ii)=0; Yt(targetsExist(ii,2)+1:end,ii)=0;
            end
            
            % what is the old energy
            stateInfo.targetsExist=targetsExist;
            stateInfo.N=N;

            if N
                stateInfo=matricesToVector(Xt,Yt,stateInfo);
                x=stateInfo.stateVec;
                origenergy=E(x,stateInfo);
            else
                origenergy=0;
            end
            
            minbest=origenergy;
%             [Edet(x) Eexc(x) Edyn(x) Econ(x)]

            %try to insert target on current detection
            N=N+1;
            targetsExist(N,:)=[newframes(1) newframes(end)];
            
            for d=exdet
                for trytrick=[1]
                    Xt(newframes,N)=Xd(t,d);Yt(newframes,N)=Yd(t,d);
                    %% try some trick here
                    %find closest forward and backward
                    usedtrick=0;
                    if trytrick
                        if t==1
                            exxdn=find(Xd(t+1,:));
                            exxdnn=find(Xd(t+2,:));
                            mindist=Inf;
                            closestn=[];
                            closestnn=[];
                            for n=exxdn
                                thisdist=norm([Xd(t,d) Yd(t,d)]-[Xd(t+1,n) Yd(t+1,n)]);
                                if thisdist<mindist;
                                    mindist=thisdist; closestn=n;
                                end
                            end
                            if ~isempty(closestn)
                                
                                mindist=Inf;
                                for nn=exxdnn
                                    thisdist=norm([Xd(t+1,n) Yd(t+1,n)]-[Xd(t+2,nn) Yd(t+2,nn)]);
                                    if thisdist<mindist;
                                        mindist=thisdist; closestnn=nn;
                                    end
                                end
                                if ~isempty(closestnn)
                                    p=[Xd(t,d) Yd(t,d)]-addnoise;
                                    p1=[Xd(t+1,closestn) Yd(t+1,closestn)]+addnoise;
                                    p2=[Xd(t+2,closestnn) Yd(t+2,closestnn)];
                                    v1=p1-p; v2=(p2-p)/2;
                                    v=v1+v2;
                                    p1=p+v; p2=p1+v;
                                    
                                    % nasty bug here... If detections are
                                    % on a grid (e.g. pixels), the new
                                    % locations might actually be 0. In
                                    % this case, discard.
                                    
                                    if all([p1 p2])
                                        %                             [norm(v1) norm(v2) norm(v)]
                                        if norm(v)<2*gridStep
                                            Xt(newframes(2),N)=p1(1); Yt(newframes(2),N)=p1(2);
                                            Xt(newframes(3),N)=p2(1); Yt(newframes(3),N)=p2(2);
                                            usedtrick=1;
                                        end
                                    end
                                end
                            end
                        elseif t==F
                            exxdn=find(Xd(t-1,:));
                            exxdnn=find(Xd(t-2,:));
                            mindist=Inf;
                            closestn=[];
                            closestnn=[];
                            for n=exxdn
                                thisdist=norm([Xd(t,d) Yd(t,d)]-[Xd(t-1,n) Yd(t-1,n)]);
                                if thisdist<mindist;
                                    mindist=thisdist; closestn=n;
                                end
                            end
                            if ~isempty(closestn)
                                
                                mindist=Inf;
                                for nn=exxdnn
                                    thisdist=norm([Xd(t-1,n) Yd(t-1,n)]-[Xd(t-2,nn) Yd(t-2,nn)]);
                                    if thisdist<mindist;
                                        mindist=thisdist; closestnn=nn;
                                    end
                                end
                                if ~isempty(closestnn)
                                    p=[Xd(t,d) Yd(t,d)]-addnoise;
                                    p1=[Xd(t-1,closestn) Yd(t-1,closestn)]+addnoise;
                                    p2=[Xd(t-2,closestnn) Yd(t-2,closestnn)];
                                    v1=p1-p; v2=(p2-p)/2;
                                    v=v1+v2;
                                    p1=p+v; p2=p1+v;
                                    
                                    %                             [norm(v1) norm(v2) norm(v)]
                                    
                                    % nasty bug here... If detections are
                                    % on a grid (e.g. pixels), the new
                                    % locations might actually be 0. In
                                    % this case, discard.
                                    
                                    if all([p1 p2])
                                        if norm(v)<2*gridStep
                                            Xt(newframes(2),N)=p1(1); Yt(newframes(2),N)=p1(2);
                                            Xt(newframes(1),N)=p2(1); Yt(newframes(1),N)=p2(2);
                                            usedtrick=1;
                                        end
                                    end
                                    
                                end
                            end
                        else
                            exxdn=find(Xd(t-1,:));
                            exxdnn=find(Xd(t+1,:));
                            mindist=Inf;
                            closestn=[];
                            closestnn=[];
                            for n=exxdn
                                thisdist=norm([Xd(t,d) Yd(t,d)]-[Xd(t-1,n) Yd(t-1,n)]);
                                if thisdist<mindist;
                                    mindist=thisdist; closestn=n;
                                end
                            end
                            if ~isempty(closestn)
                                
                                mindist=Inf;
                                for nn=exxdnn
                                    thisdist=norm([Xd(t,n) Yd(t,n)]-[Xd(t+1,nn) Yd(t+1,nn)]);
                                    if thisdist<mindist;
                                        mindist=thisdist; closestnn=nn;
                                    end
                                end
                                if ~isempty(closestnn)
                                    p=[Xd(t,d) Yd(t,d)]-addnoise;
                                    p1=[Xd(t-1,closestn) Yd(t-1,closestn)]+addnoise;
                                    p2=[Xd(t+1,closestnn) Yd(t+1,closestnn)];
                                    v1=p-p1; v2=p2-p;
                                    v=(v1+v2)/2;
                                    p1=p-v; p2=p+v;
                                    
%                                     [norm(v1) norm(v2) norm(v)]
                                    % nasty bug here... If detections are
                                    % on a grid (e.g. pixels), the new
                                    % locations might actually be 0. In
                                    % this case, discard.
                                    
                                    if all([p1 p2])
                                    if norm(v)<2*gridStep
                                        Xt(newframes(1),N)=p1(1); Yt(newframes(1),N)=p1(2);
                                        Xt(newframes(3),N)=p2(1); Yt(newframes(3),N)=p2(2);
                                        usedtrick=1;
                                    end
                                    end
                                    
                                end
                            end
                        end
                    end
                    %                 printMessage(2,'new  %i\n',N);                 
                    stateInfo.N=N; stateInfo.targetsExist=targetsExist;
                    stateInfo=matricesToVector(Xt,Yt,stateInfo);
                    xt=stateInfo.stateVec;
%                     xt=generateXfromXY(Xt,Yt);
                    newenergy=E(xt,stateInfo);            
                    
%                     [Edet(xt) Eexc(xt) Edyn(xt) Econ(xt)]
%                     [newenergy origenergy]
                    
                    
                    if newenergy-origenergy<bestgain && newenergy>-1e5 % numerical issue
                        minbest=newenergy;
                        best=[t d cf];
                        bestX=Xt(newframes,N);
                        bestY=Yt(newframes,N);
                        bestnewframes=newframes;
                        cntbest=cntbest+1;
                        bestgain=newenergy-origenergy;
                        bestnewgain=newenergy;
                        bestoriggain=origenergy;
                        bestimportanti=importanti;
                        woohoo=0;
                        if usedtrick
                            woohoo=1;
                        end
                    end
                    %                 printMessage(2,'new2 %i\n',N);
                end
            end
            %             printMessage(2,'old? %i\n',N);
            targetsExist=te;
            stateInfo=stateInfoOld;
        end
        
        if cntbest>5
            break;
            
            
            %             N=size(te,1);
            %             printMessage(2,'break old %i\n',N);
        end
    end
    N=Nold;
    targetsExist=te;
    stateInfo=stateInfoOld;
    
    % there is a bug, sometimes...
    alladdenergy(j)=bestgain;

    if bestgain<-1000000
        printMessage(2,'something wrong. bg: %f, best: %i %i %i\n',bestgain, best);
        
    elseif bestgain<0
        
        N=N+1;
        stateInfo.N=N;
        
        d=best(2);
        cf=best(3);
        t=best(1);
%         printMessage(2,'add target %i around frame %i. Energy gain: %f (%f %f). Trick: %i\n',N,cf,bestgain,bestnewgain,bestoriggain,woohoo);
%         printMessage(2,', %i (%i,%i, %.2f|%.2f)',N,cf,woohoo,bestX(2),bestY(2));
        printMessage(3,', %i',N);
        newframes=cf-1:cf+1;
        X(newframes,N)=bestX;Y(newframes,N)=bestY;
        targetsExist(N,:)=[newframes(1) newframes(end)];
        stateInfo.targetsExist=targetsExist;
        
        Xold=X; Yold=Y; % replace X Y
        Nold=N;
        stateInfoOld=stateInfo;
    else
%                 printMessage(2,'do not add target %i around frame %i. Energy gain: %f \n',N+1,cf,bestgain);
        break;
    end
    
end

stateInfo=matricesToVector(X,Y,stateInfo);
% toc
printMessage(2,'...done\n');
% [Xold X]
% pause
end