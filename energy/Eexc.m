function [fx dfx ds]=Eexc(x,stateInfo)

% Exclusion (collision avoidance) term.
% pairwise distances between all targets in all frames
% are penalized according to 1/distance^2 and scaled with 
% target size (usually 350)
% 


global sceneInfo;


F=stateInfo.F;
tiToInd=stateInfo.tiToInd;
gridStep=sceneInfo.targetSize;
% x=stateInfo.stateVec;

% grab state vector
[X Y]=vectorToMatrices(x,stateInfo);

% initialize with 0
fx=0;
dfx=zeros(length(x),1);
ds=zeros(size(X));
gs_squared=gridStep*gridStep;

% if ~all(size(X)==size(tiToInd)), pause; end

cnt=0;
% in each frame
for t=1:F
    
    % for each existing target
    targets=find(X(t,:));
    for i1=1:numel(targets)
        i=targets(i1);
        xind=tiToInd(t,i);
        yind=xind+1;
        
        cnt=cnt+1;
        tmpds=0;
        % for all other targets j ~= i
        otheri=[targets(1:i1-1) targets(i1+1:end)];
        for j=otheri

            % target i = (a,b) 
            % target j = (c,d)
            a=X(t,i); b=Y(t,i);
            c=X(t,j); d=Y(t,j);

            % exclusion term
            term = gs_squared / (a^2 + b^2 + c^2 + d^2 -2*a*c - 2*b*d); %  gs^2 / |Xi-Xj|^2           

            % add to energy
            fx=fx+term;
            tmpds=tmpds+term;
            
            % compute gradient if needed
            if nargout > 1
                denom=(a^2-2*a*c+b^2-2*b*d+c^2+d^2)^2;
                tmp=gs_squared/denom;
                dda=-(2*a-2*c)*tmp;
                ddb=-(2*b-2*d)*tmp;
%                 [length(x) xind]
%                 pause
                dfx(xind)=dfx(xind)+dda+dda;
                dfx(yind)=dfx(yind)+ddb+ddb;
            end

        end
        ds(t,i)=tmpds;
        
    end
end

end