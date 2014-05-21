%%%% WRONG GRADIENT  ???? %%%%

function [fx dfx ds]=Eper(x,stateInfo)

% The persistence term
% it calculates the distance of the
% start and end point of each trajectory
% to the closest border (image or tracking
% area) and penalizes with a sigmoid
% 


global sceneInfo;

% get state info
[~, N F targetsExist X Y]=getStateInfo(stateInfo);
gridStep=sceneInfo.targetSize;


fx=0;
dfx=zeros(length(x),1);
dfxind=0;

itToInd=stateInfo.tiToInd;
areaLimits=sceneInfo.trackingArea;

minX=areaLimits(1);
maxX=areaLimits(2);
minY=areaLimits(3);
maxY=areaLimits(4);

A=1*gridStep; %% offset
B=1/gridStep; %% scale
AB=A*B;



ds=zeros(size(X));
cnt=0;
for i=1:N
   
   st=targetsExist(i,1);
   en=targetsExist(i,2);
   
   if st>1
       
       x=X(st,i);
       y=Y(st,i);
       
       % dist left
       
       dl=abs(minX-x);
       % dist right
       dr=abs(maxX-x);
       % dist up
       du=abs(minY-y);
       % dist down
       dd=abs(maxY-y);
       distances=[dl dr du dd];
       
       [dm ci]=min(distances);
       d=1/(1+exp(AB-B*dm));
       
	   ds(st,i)=d;
      
       fx=fx+d;
       if nargout>1
           if ci<=2 % x
               dfx(itToInd(st,i))=B*exp(AB-B*dm)/((exp(AB-B*dm)+1)^2);
           else % y
               dfx(itToInd(st,i)+1)=B*exp(AB-B*dm)/((exp(AB-B*dm)+1)^2);
           end
       end
       
   end
   
   
   if en<F
       
       x=X(en,i);
       y=Y(en,i);
       
       % dist left
       dl=abs(minX-x);
       % dist right
       dr=abs(maxX-x);
       % dist up
       du=abs(minY-y);
       % dist down
       dd=abs(maxY-y);
       
       distances=[dl dr du dd];
       
       [dm ci]=min(distances);
       d=1/(1+exp(AB-B*dm));
       ds(en,i)=d;
       fx=fx+d;
       if nargout>1
           if ci<=2 % x
               dfx(itToInd(en,i))=B*exp(AB-B*dm)/((exp(AB-B*dm)+1)^2);
           else % y
               dfx(itToInd(en,i)+1)=B*exp(AB-B*dm)/((exp(AB-B*dm)+1)^2);
           end
       end
   end
   
   
   cnt=cnt+1;
   %         ds(cnt)=d;
end


