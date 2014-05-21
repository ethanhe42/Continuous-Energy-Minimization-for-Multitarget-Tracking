function [fx dfx]=Edet(x,stateInfo)

% The Detection Energy Term
%
% The energy should be minimized when the targets pass through
% detections [Xd Yd]


global detections sceneInfo opt;

N=stateInfo.N;
F=stateInfo.F;
targetsExist=stateInfo.targetsExist;
gridStep=sceneInfo.targetSize;
% x=stateInfo.stateVec;

% convert state vector to matrix representation
[X Y]=vectorToMatrices(x,stateInfo);

% initialize return values with 0
fx=0;
dfx=zeros(length(x),1);


% if no targets are present, return 0
if isempty(x) || ~N
    return;
end

% otherwise compute the actual Edet-value
dfxind=0;


csig=gridStep*gridStep;
lambda=opt.lambda;

occ=ones(F,N);
occx=zeros(F,N);
occy=zeros(F,N);
GG=zeros(F,N);
GGx=zeros(F,N);
GGy=zeros(F,N);
if opt.occ && nargout > 1
    [occ vis occx occy ddvix ddviy]=computeOcclusions2(X,Y);
elseif opt.occ
    [occ vis]=computeOcclusions2(X,Y);
end

for id=1:N
    for t=targetsExist(id,1):targetsExist(id,2)
        
        % current x,y position
        x=X(t,id); y=Y(t,id);
        
        % v * lambda
        d=occ(t,id)*lambda;
        
        % matrix implementation
        XXd=detections(t).xp-x;
        YYd=detections(t).yp-y;
        SSd=detections(t).sc;
        
        XXdsq=XXd.^2; YYdsq=YYd.^2;
        
        GG(t,id)=-sum((csig*SSd) ./ ((XXdsq + YYdsq) + csig));
        
        % add current value to total energy
        fx=fx+d+GG(t,id);
        
        
        % if derivative requested
        if nargout>1
            ddfx=occx(t,id)*lambda;
            ddfy=occy(t,id)*lambda;
            
            GGx(t,id)=-sum(2*csig*SSd.*XXd ./ ((XXdsq + YYdsq) + csig).^2);
            GGy(t,id)=-sum(2*csig*SSd.*YYd ./ ((XXdsq + YYdsq) + csig).^2);
            
            ddfx=ddfx+GGx(t,id);
            ddfy=ddfy+GGy(t,id);
            
            
            dfxind=dfxind+1;
            dfx(dfxind)=dfx(dfxind)+ddfx;
            dfxind=dfxind+1;
            dfx(dfxind)=dfx(dfxind)+ddfy;         
            
        end        
        
    end
end

end




