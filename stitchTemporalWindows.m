function stNew=stitchTemporalWindows(allstInfo,stInfo,allwins,overlapSize)

% for w=1:wincnt

stNew=stInfo;
global sceneInfo

Noffset=0;
for w=1:size(allwins,1)-1
    Xtmp=zeros(overlapSize,0);    Ytmp=zeros(overlapSize,0);
    Wtmp=zeros(overlapSize,0);    Htmp=zeros(overlapSize,0);
    Xitmp=zeros(overlapSize,0);    Yitmp=zeros(overlapSize,0);
    Xgptmp=zeros(overlapSize,0);    Ygptmp=zeros(overlapSize,0);
    
    
    winframes=(allwins(w,1):allwins(w,2))';
    winframes2=(allwins(w,2)+1:allwins(w+1,2))';
    olframes=(allwins(w+1,1):allwins(w,2))';
    
    curN=1+Noffset;
    N=size(allstInfo(w).X,2);newIDs=1:N+Noffset;
    N1=N;
    newIDs1=newIDs;
    
%     Xtmp(:,newIDs)=allstInfo(w).X(end-overlapSize+1:end,:);
%     Ytmp(:,newIDs)=allstInfo(w).Y(end-overlapSize+1:end,:);
    Xtmp(:,newIDs)=stNew.X(olframes,newIDs);Ytmp(:,newIDs)=stNew.Y(olframes,newIDs);
    Wtmp(:,newIDs)=stNew.W(olframes,newIDs);Htmp(:,newIDs)=stNew.H(olframes,newIDs);
    Xitmp(:,newIDs)=stNew.Xi(olframes,newIDs);    Yitmp(:,newIDs)=stNew.Yi(olframes,newIDs);
    if isfield(stNew,'Xgp')
        Xgptmp(:,newIDs)=stNew.Xgp(olframes,newIDs);    Ygptmp(:,newIDs)=stNew.Ygp(olframes,newIDs);
    end
%     Xtmp
    
    curN=curN+N1+Noffset;
    N=size(allstInfo(w+1).X,2);newIDs=N1+1+Noffset:N1+N+Noffset;
    Xtmp(:,newIDs)=stNew.X(olframes,newIDs);    Ytmp(:,newIDs)=stNew.Y(olframes,newIDs);
    Wtmp(:,newIDs)=stNew.W(olframes,newIDs);    Htmp(:,newIDs)=stNew.H(olframes,newIDs);
    Xitmp(:,newIDs)=stNew.Xi(olframes,newIDs);    Yitmp(:,newIDs)=stNew.Yi(olframes,newIDs);
    if isfield(stNew,'Xgp')
        Xgptmp(:,newIDs)=stNew.Xgp(olframes,newIDs);    Ygptmp(:,newIDs)=stNew.Ygp(olframes,newIDs);
    end
%     Xtmp
    
    stTmp.X=Xtmp;stTmp.Y=Ytmp;
    mh=zeros(1,size(Xtmp,2));
    [prox, proxt, proxcost]=getSplineProximity(stTmp);
    
    tomerge=[];
    
    for id=newIDs1
        [m, minc]=min(prox(id,:));
        if m<sceneInfo.targetSize/2
            tomerge=[tomerge; id minc];
        end
    end
%     tomerge
    %%
    for ol=1:size(tomerge,1)
        id1=tomerge(ol,1);id2=tomerge(ol,2);
        
        noT1=~Xtmp(:,id1);noT2=~Xtmp(:,id2);    Xtmp(noT1,id1)=Xtmp(noT1,id2);Xtmp(noT2,id2)=Xtmp(noT2,id1);
        noT1=~Ytmp(:,id1);noT2=~Ytmp(:,id2);    Ytmp(noT1,id1)=Ytmp(noT1,id2);Ytmp(noT2,id2)=Ytmp(noT2,id1);
        noT1=~Xitmp(:,id1);noT2=~Xitmp(:,id2);  Xitmp(noT1,id1)=Xitmp(noT1,id2);Xitmp(noT2,id2)=Xitmp(noT2,id1);
        noT1=~Yitmp(:,id1);noT2=~Yitmp(:,id2);  Yitmp(noT1,id1)=Yitmp(noT1,id2);Yitmp(noT2,id2)=Yitmp(noT2,id1);
        if isfield(stNew,'Xgp')
            noT1=~Xgptmp(:,id1);noT2=~Xgptmp(:,id2);Xgptmp(noT1,id1)=Xgptmp(noT1,id2);Xgptmp(noT2,id2)=Xgptmp(noT2,id1);
            noT1=~Ygptmp(:,id1);noT2=~Ygptmp(:,id2);Ygptmp(noT1,id1)=Ygptmp(noT1,id2);Ygptmp(noT2,id2)=Ygptmp(noT2,id1);
        end
        noT1=~Wtmp(:,id1);noT2=~Wtmp(:,id2);    Wtmp(noT1,id1)=Wtmp(noT1,id2);Wtmp(noT2,id2)=Wtmp(noT2,id1);
        noT1=~Htmp(:,id1);noT2=~Htmp(:,id2);    Htmp(noT1,id1)=Htmp(noT1,id2);Htmp(noT2,id2)=Htmp(noT2,id1);

        stNew.X(olframes,id1)=mean([Xtmp(:,id1) Xtmp(:,id2)],2);
        stNew.X(winframes2,id1)=stNew.X(winframes2,id2);
        stNew.X(olframes,id2)=0;stNew.X(winframes2,id2)=0;

        stNew.Y(olframes,id1)=mean([Ytmp(:,id1) Ytmp(:,id2)],2);
        stNew.Y(winframes2,id1)=stNew.Y(winframes2,id2);
        stNew.Y(olframes,id2)=0;stNew.Y(winframes2,id2)=0;
        
        stNew.W(olframes,id1)=mean([Wtmp(:,id1) Wtmp(:,id2)],2);
        stNew.W(winframes2,id1)=stNew.W(winframes2,id2);
        stNew.W(olframes,id2)=0;stNew.W(winframes2,id2)=0;

        stNew.H(olframes,id1)=mean([Htmp(:,id1) Htmp(:,id2)],2);
        stNew.H(winframes2,id1)=stNew.H(winframes2,id2);
        stNew.H(olframes,id2)=0;stNew.H(winframes2,id2)=0;

        stNew.Xi(olframes,id1)=mean([Xitmp(:,id1) Xitmp(:,id2)],2);
        stNew.Xi(winframes2,id1)=stNew.Xi(winframes2,id2);
        stNew.Xi(olframes,id2)=0;stNew.Xi(winframes2,id2)=0;

        stNew.Yi(olframes,id1)=mean([Yitmp(:,id1) Yitmp(:,id2)],2);
        stNew.Yi(winframes2,id1)=stNew.Yi(winframes2,id2);
        stNew.Yi(olframes,id2)=0;stNew.Yi(winframes2,id2)=0;

        if isfield(stNew,'Xgp')
            stNew.Xgp(olframes,id1)=mean([Xgptmp(:,id1) Xgptmp(:,id2)],2);
            stNew.Xgp(winframes2,id1)=stNew.Xgp(winframes2,id2);
            stNew.Xgp(olframes,id2)=0;stNew.Xgp(winframes2,id2)=0;

            stNew.Ygp(olframes,id1)=mean([Ygptmp(:,id1) Ygptmp(:,id2)],2);
            stNew.Ygp(winframes2,id1)=stNew.Ygp(winframes2,id2);
            stNew.Ygp(olframes,id2)=0;stNew.Ygp(winframes2,id2)=0;
        end
        %         allstInfo(w+1).X(end-overlapSize+1:end,:)
        
        
%         allstInfo(w).X(end-overlapSize+1:end,:)=stNew.X(olframes,1:N1);
%         allstInfo(w).Y(end-overlapSize+1:end,:)=stNew.Y(olframes,1:N1);
        
%         allstInfo(w+1).X(1:overlapSize,:)=stNew.X(olframes,:);
%         allstInfo(w+1).Y(1:overlapSize,:)=stNew.Y(olframes,:);
    end
    Noffset=Noffset+N1;
end

[X Y stNew]=cleanState(stNew.X, stNew.Y, stNew);

end