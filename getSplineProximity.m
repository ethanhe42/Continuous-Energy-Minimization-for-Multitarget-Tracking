function [prox, proxt, proxcost, merpairs]=getSplineProximity(stInfo)
% compute pairwise distances along all splines
%
%

[T,F]=size(stInfo.X);

% global opt sceneInfo
proxcost=zeros(F);
prox=Inf*ones(F);
proxt=zeros(F);
merpairs=[];


% pairwise distance between traj
X=stInfo.X; Y=stInfo.Y;

Xo=X; Yo=Y;

%
% isactive=false(1,F);
% for m=1:F
%     if any(labeling==m), isactive(m)=1;end
% end

% for m1=1:F
%     rounded(m1).breaks=round(mh(m1).breaks);
% end

% for m1=1:F
%     s1=mh(m1).start; e1=mh(m1).end;
%     s1=find(X(:,m1),1,'first');e1=find(X(:,m1),1,'last');
%     for m2=m1+1:F
%         s2=mh(m2).start; e2=mh(m2).end;
%         s2=find(X(:,m2),1,'first');e2=find(X(:,m2),1,'last');
%
%         mus=max(s1,s2);        mue=min(e1,e2);
%         if mus<=mue
%
%             evalPerFr=1;
%             tt=mus:mue;
%             %%%% CAREFUL! only for integer tt!!!
%             splinexy1=[X(tt,m1)';Y(tt,m1)'];
%             splinexy2=[X(tt,m2)';Y(tt,m2)'];
%
%             d=splinexy1-splinexy2;
%             squares=(d.^2);
%
%             sqrts=sqrt(sum(squares));
%             distalongtraj=sqrts/1;
%             [mindist atframe]=min(distalongtraj);
%
% %             [m1 m2]
% %             distalongtraj
%             prox(m1,m2)=mindist;
%             proxt(m1,m2)=mus+atframe/evalPerFr-1;
%         end
%     end
%
% end

% vectorized version
X(~X)=Inf; Y(~Y)=Inf;

prox2=Inf*ones(F,F,T);
for t=1:T
    tmp=repmat(X(t,:),F,1);
    ax=tmp(:)';
    tmp=tmp';
    bx=tmp(:)';
    
    tmp=repmat(Y(t,:),F,1);
    ay=tmp(:)';
    tmp=tmp';
    by=tmp(:)';
    
    dx=(ax-bx).^2;
    dy=(ay-by).^2;
    d=sqrt(dx+dy);
    
    % inf-inf = NaN
%     if ~isempty(find(isnan(d)))
%         d
%         ax
%         ay
%     end
    d(isnan(d))=Inf;
    
%     [size(prox2) F]
% if t==14
%     size(d)
%     t
%     fprintf('%.15f\n',ax);
%     fprintf('%.15f\n',bx);
%     fprintf('%.15f\n',ay);
%     fprintf('%.15f\n',by);
%     fprintf('%.15f\n',dx);
%     fprintf('%.15f\n',dy);
%     fprintf('%.15f\n',d);
% ax
% ay
% bx
% by
%     dx
%     dy
%      d
% end
%     pause
    prox2(:,:,t)=reshape(d,F,F);
end
% prox2
proxfromloop=prox2;
[prox2, proxt2]=min(prox2,[],3);


prox2=triu(prox2); prox2(isinf(prox2))=Inf;
prox2(~~(tril(ones(F))))=Inf;
% prox
% prox2
prox=prox2;

end