function [fx dfx ds]=Edyn(x,stateInfo)

% The dynamic model
% This is a simple constant velocity model
% that is minimized when succesive velocity vectors
% of the same target are the same
% 


global sceneInfo;

N=stateInfo.N;
targetsExist=stateInfo.targetsExist;
gridStep=sceneInfo.targetSize;
% x=stateInfo.stateVec;

% convert state vector to matrix representation
[X Y]=vectorToMatrices(x, stateInfo);


fx=0;
dfx=zeros(length(x),1);
ds=zeros(size(X));
cnt=0;
xind=1;
for i=1:N
    tlength=diff(targetsExist(i,:))+1;

    xind=xind+2;
    yind=xind+1;

    
	% Edyn is 0 for first and last frame
    for t=targetsExist(i,1)+1:targetsExist(i,2)-1
		% (a,b) = past frame position
		% (c,d) = current frame position
		% (e,f) = next frame position
        a=X(t-1,i);
        b=Y(t-1,i);
        c=X(t,i);
        d=Y(t,i);
        e=X(t+1,i);
        f=Y(t+1,i);
        
		% squared norm of the two succesive vectors
        diffterm= a^2 + a*(2*e - 4*c) + b^2 + b*(2*f - 4*d) + 4*c^2  - 4*c*e + 4*d^2  - 4*d*f + e^2  + f^2;
        cnt=cnt+1;ds(t,i)=diffterm;
        fx=fx+diffterm;
        
		% derivative
        if nargout>1
            dfx(xind-2)=dfx(xind-2) +   (2*a-4*c+2*e);
            dfx(yind-2)=dfx(yind-2) +   (2*b-4*d+2*f);
            dfx(xind)=dfx(xind)     +   (-4*a+8*c-4*e);
            dfx(yind)=dfx(yind)     +   (-4*b+8*d-4*f);
            dfx(xind+2)=dfx(xind+2) +   (2*a-4*c+2*e);
            dfx(yind+2)=dfx(yind+2) +   (2*b-4*d+2*f);
        end
               
        xind=xind+2;
        yind=xind+1;
        
    end
    if tlength>1
        xind=xind+2;
    end
end

fx=fx/gridStep;		% normalize
dfx=dfx/gridStep;	% normalize
ds=ds/gridStep;
end