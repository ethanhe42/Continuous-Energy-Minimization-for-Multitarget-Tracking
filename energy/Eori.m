function [fx dfx ds]=Eori(x,stateInfo)
%
% Orientation
%


global detections detMatrices sceneInfo;

N=stateInfo.N;
targetsExist=stateInfo.targetsExist;
gridStep=sceneInfo.targetSize;
csig=gridStep;
% x=stateInfo.stateVec;

% convert state vector to matrix representation
[X Y]=vectorToMatrices(x, stateInfo);


fx=0;
dfx=zeros(length(x),1);
ds=zeros(size(X));
cnt=0;
xind=1;

d1=[1 0]; d1=d1./norm(d1);
d2=[1 0]; d2=d2./norm(d2);
so=2;
epsilon=1e-1;

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
        %         d1x=d1(1); d1y=d1(2);
        %         d2x=d2(1); d2y=d2(2);
        
        
        
        v1=([c; d] - [a; b]);
        v2=([e; f] - [c; d]);
        
        %         norm(v1)
        m1=sqrt(v1(1)^2 + v1(2)^2+epsilon);
        m2=sqrt(v2(1)^2 + v2(2)^2+epsilon);
        %         o1=v1./m1;
        %         o2=v2./m2;
        %         if m1<1e-5 || m2<1e-5
        %             fprintf('%.15f\n',[a b c d e f]);
        %             [m1 m2]
        %             pause
        %
        %         end
        
        ndet2=numel(find(detMatrices.Xd(t,:)));
        dets2=[detMatrices.Xd(t,1:ndet2);detMatrices.Yd(t,1:ndet2)];
        reppt=repmat([c;d],1,ndet2);
        %         dets2
        %         reppt
        ddists2=sqrt(sum((dets2-reppt).^2));
        %         ddists2
        ddists2=find(ddists2<csig*10);
        %         ddists2
        %         pause
        
        for det2=ddists2
            d2x=detMatrices.Dx(t,det2);d2y=detMatrices.Dy(t,det2);d2=[d2x d2y];
            det2x=dets2(1,det2);det2y=dets2(2,det2);
            dist2=sqrt((c-det2x)^2 + (d-det2y)^2);
            
            f2=1/(1+exp(-m2+so));
            f2=f2*csig/(dist2^2+csig);
            f2acos=f2*acos((v2(1)*d2x + v2(2)*d2y) / m2 );
            
            ndet1=numel(find(detMatrices.Xd(t-1,:)));
            
            dets1=[detMatrices.Xd(t-1,1:ndet1);detMatrices.Yd(t-1,1:ndet1)];
            dirs1=[detMatrices.Dx(t-1,1:ndet1);detMatrices.Dy(t-1,1:ndet1)];
            for det1=1:ndet1
                d1x=dirs1(1,det1);d1y=dirs1(2,det1);d1=[d1x d1y];
                det1x=dets1(1,det1);det1y=dets1(2,det1);
                dist1=sqrt((a-det1x)^2 + (b-det1y)^2);
                
                if dist1<csig*10
                    f1=1/(1+exp(-m1+so));
                    
                    
                    
                    f1=f1*csig/(dist1^2+csig);
                    
                    f1acos=f1*acos((v1(1)*d1x + v1(2)*d1y) / m1);
                    
                    %                     if m1<1e-5, f1acos=0; end
                    %                     if m2<1e-5, f2acos=0; end
                    
                    
                    diffterm=f1acos+f2acos;
                    
                    
                    cnt=cnt+1;ds(t,i)=diffterm;
                    fx=fx+diffterm;
                    
                    % derivative
                    if nargout>1
                        % diff1 =(acos(- (d1x*(a - c))/((a - c)^2 + (b - d)^2)^(1/2) - (d1y*(b - d))/((a - c)^2 + (b - d)^2)^(1/2))*exp(so - ((a - c)^2 + (b - d)^2)^(1/2))*(2*a - 2*c))/(2*((a - c)^2 + (b - d)^2)^(1/2)*(exp(so - ((a - c)^2 + (b - d)^2)^(1/2)) + 1)^2) - ((d1x*(a - c)*(2*a - 2*c))/(2*((a - c)^2 + (b - d)^2)^(3/2)) - d1x/((a - c)^2 + (b - d)^2)^(1/2) + (d1y*(b - d)*(2*a - 2*c))/(2*((a - c)^2 + (b - d)^2)^(3/2)))/((exp(so - ((a - c)^2 + (b - d)^2)^(1/2)) + 1)*(1 - ((d1x*(a - c))/((a - c)^2 + (b - d)^2)^(1/2) + (d1y*(b - d))/((a - c)^2 + (b - d)^2)^(1/2))^2)^(1/2));
                        % diff2 =(acos(- (d1x*(a - c))/((a - c)^2 + (b - d)^2)^(1/2) - (d1y*(b - d))/((a - c)^2 + (b - d)^2)^(1/2))*exp(so - ((a - c)^2 + (b - d)^2)^(1/2))*(2*b - 2*d))/(2*((a - c)^2 + (b - d)^2)^(1/2)*(exp(so - ((a - c)^2 + (b - d)^2)^(1/2)) + 1)^2) - ((d1x*(a - c)*(2*b - 2*d))/(2*((a - c)^2 + (b - d)^2)^(3/2)) - d1y/((a - c)^2 + (b - d)^2)^(1/2) + (d1y*(b - d)*(2*b - 2*d))/(2*((a - c)^2 + (b - d)^2)^(3/2)))/((exp(so - ((a - c)^2 + (b - d)^2)^(1/2)) + 1)*(1 - ((d1x*(a - c))/((a - c)^2 + (b - d)^2)^(1/2) + (d1y*(b - d))/((a - c)^2 + (b - d)^2)^(1/2))^2)^(1/2));
                        % diff3 =((d1x*(a - c)*(2*a - 2*c))/(2*((a - c)^2 + (b - d)^2)^(3/2)) - d1x/((a - c)^2 + (b - d)^2)^(1/2) + (d1y*(b - d)*(2*a - 2*c))/(2*((a - c)^2 + (b - d)^2)^(3/2)))/((exp(so - ((a - c)^2 + (b - d)^2)^(1/2)) + 1)*(1 - ((d1x*(a - c))/((a - c)^2 + (b - d)^2)^(1/2) + (d1y*(b - d))/((a - c)^2 + (b - d)^2)^(1/2))^2)^(1/2)) - ((d2x*(c - e)*(2*c - 2*e))/(2*((c - e)^2 + (d - f)^2)^(3/2)) - d2x/((c - e)^2 + (d - f)^2)^(1/2) + (d2y*(d - f)*(2*c - 2*e))/(2*((c - e)^2 + (d - f)^2)^(3/2)))/((exp(so - ((c - e)^2 + (d - f)^2)^(1/2)) + 1)*(1 - ((d2x*(c - e))/((c - e)^2 + (d - f)^2)^(1/2) + (d2y*(d - f))/((c - e)^2 + (d - f)^2)^(1/2))^2)^(1/2)) - (acos(- (d1x*(a - c))/((a - c)^2 + (b - d)^2)^(1/2) - (d1y*(b - d))/((a - c)^2 + (b - d)^2)^(1/2))*exp(so - ((a - c)^2 + (b - d)^2)^(1/2))*(2*a - 2*c))/(2*((a - c)^2 + (b - d)^2)^(1/2)*(exp(so - ((a - c)^2 + (b - d)^2)^(1/2)) + 1)^2) + (exp(so - ((c - e)^2 + (d - f)^2)^(1/2))*acos(- (d2x*(c - e))/((c - e)^2 + (d - f)^2)^(1/2) - (d2y*(d - f))/((c - e)^2 + (d - f)^2)^(1/2))*(2*c - 2*e))/(2*((c - e)^2 + (d - f)^2)^(1/2)*(exp(so - ((c - e)^2 + (d - f)^2)^(1/2)) + 1)^2);
                        % diff4 =((d1x*(a - c)*(2*b - 2*d))/(2*((a - c)^2 + (b - d)^2)^(3/2)) - d1y/((a - c)^2 + (b - d)^2)^(1/2) + (d1y*(b - d)*(2*b - 2*d))/(2*((a - c)^2 + (b - d)^2)^(3/2)))/((exp(so - ((a - c)^2 + (b - d)^2)^(1/2)) + 1)*(1 - ((d1x*(a - c))/((a - c)^2 + (b - d)^2)^(1/2) + (d1y*(b - d))/((a - c)^2 + (b - d)^2)^(1/2))^2)^(1/2)) - ((d2x*(c - e)*(2*d - 2*f))/(2*((c - e)^2 + (d - f)^2)^(3/2)) - d2y/((c - e)^2 + (d - f)^2)^(1/2) + (d2y*(d - f)*(2*d - 2*f))/(2*((c - e)^2 + (d - f)^2)^(3/2)))/((exp(so - ((c - e)^2 + (d - f)^2)^(1/2)) + 1)*(1 - ((d2x*(c - e))/((c - e)^2 + (d - f)^2)^(1/2) + (d2y*(d - f))/((c - e)^2 + (d - f)^2)^(1/2))^2)^(1/2)) - (acos(- (d1x*(a - c))/((a - c)^2 + (b - d)^2)^(1/2) - (d1y*(b - d))/((a - c)^2 + (b - d)^2)^(1/2))*exp(so - ((a - c)^2 + (b - d)^2)^(1/2))*(2*b - 2*d))/(2*((a - c)^2 + (b - d)^2)^(1/2)*(exp(so - ((a - c)^2 + (b - d)^2)^(1/2)) + 1)^2) + (exp(so - ((c - e)^2 + (d - f)^2)^(1/2))*acos(- (d2x*(c - e))/((c - e)^2 + (d - f)^2)^(1/2) - (d2y*(d - f))/((c - e)^2 + (d - f)^2)^(1/2))*(2*d - 2*f))/(2*((c - e)^2 + (d - f)^2)^(1/2)*(exp(so - ((c - e)^2 + (d - f)^2)^(1/2)) + 1)^2);
                        % diff5 =((d2x*(c - e)*(2*c - 2*e))/(2*((c - e)^2 + (d - f)^2)^(3/2)) - d2x/((c - e)^2 + (d - f)^2)^(1/2) + (d2y*(d - f)*(2*c - 2*e))/(2*((c - e)^2 + (d - f)^2)^(3/2)))/((exp(so - ((c - e)^2 + (d - f)^2)^(1/2)) + 1)*(1 - ((d2x*(c - e))/((c - e)^2 + (d - f)^2)^(1/2) + (d2y*(d - f))/((c - e)^2 + (d - f)^2)^(1/2))^2)^(1/2)) - (exp(so - ((c - e)^2 + (d - f)^2)^(1/2))*acos(- (d2x*(c - e))/((c - e)^2 + (d - f)^2)^(1/2) - (d2y*(d - f))/((c - e)^2 + (d - f)^2)^(1/2))*(2*c - 2*e))/(2*((c - e)^2 + (d - f)^2)^(1/2)*(exp(so - ((c - e)^2 + (d - f)^2)^(1/2)) + 1)^2);
                        % diff6 =((d2x*(c - e)*(2*d - 2*f))/(2*((c - e)^2 + (d - f)^2)^(3/2)) - d2y/((c - e)^2 + (d - f)^2)^(1/2) + (d2y*(d - f)*(2*d - 2*f))/(2*((c - e)^2 + (d - f)^2)^(3/2)))/((exp(so - ((c - e)^2 + (d - f)^2)^(1/2)) + 1)*(1 - ((d2x*(c - e))/((c - e)^2 + (d - f)^2)^(1/2) + (d2y*(d - f))/((c - e)^2 + (d - f)^2)^(1/2))^2)^(1/2)) - (exp(so - ((c - e)^2 + (d - f)^2)^(1/2))*acos(- (d2x*(c - e))/((c - e)^2 + (d - f)^2)^(1/2) - (d2y*(d - f))/((c - e)^2 + (d - f)^2)^(1/2))*(2*d - 2*f))/(2*((c - e)^2 + (d - f)^2)^(1/2)*(exp(so - ((c - e)^2 + (d - f)^2)^(1/2)) + 1)^2);
                        
%                         t1=(epsilon + ((a - c)^2 + (b - d)^2)^(1/2));
%                         t2=((2*a - 2*c)*(d1x*(a - c) + d1y*(b - d)));
%                         t3=(2*((a - c)^2 + (b - d)^2)^(1/2)*t1^2);
                        diff1 = (csig*(d1x/(epsilon + (a - c)^2 + (b - d)^2)^(1/2) - ((2*a - 2*c)*(d1x*(a - c) + d1y*(b - d)))/(2*(epsilon + (a - c)^2 + (b - d)^2)^(3/2))))/((exp(so - (epsilon + (a - c)^2 + (b - d)^2)^(1/2)) + 1)*(1 - (d1x*(a - c) + d1y*(b - d))^2/(epsilon + (a - c)^2 + (b - d)^2))^(1/2)*(csig + (a - det1x)^2 + (b - det1y)^2)) - (csig*(pi - acos((d1x*(a - c) + d1y*(b - d))/(epsilon + (a - c)^2 + (b - d)^2)^(1/2)))*(2*a - 2*det1x))/((exp(so - (epsilon + (a - c)^2 + (b - d)^2)^(1/2)) + 1)*(csig + (a - det1x)^2 + (b - det1y)^2)^2) + (csig*exp(so - (epsilon + (a - c)^2 + (b - d)^2)^(1/2))*(pi - acos((d1x*(a - c) + d1y*(b - d))/(epsilon + (a - c)^2 + (b - d)^2)^(1/2)))*(2*a - 2*c))/(2*(exp(so - (epsilon + (a - c)^2 + (b - d)^2)^(1/2)) + 1)^2*(csig + (a - det1x)^2 + (b - det1y)^2)*(epsilon + (a - c)^2 + (b - d)^2)^(1/2));
                        diff2 = (csig*(d1y/(epsilon + (a - c)^2 + (b - d)^2)^(1/2) - ((2*b - 2*d)*(d1x*(a - c) + d1y*(b - d)))/(2*(epsilon + (a - c)^2 + (b - d)^2)^(3/2))))/((exp(so - (epsilon + (a - c)^2 + (b - d)^2)^(1/2)) + 1)*(1 - (d1x*(a - c) + d1y*(b - d))^2/(epsilon + (a - c)^2 + (b - d)^2))^(1/2)*(csig + (a - det1x)^2 + (b - det1y)^2)) - (csig*(pi - acos((d1x*(a - c) + d1y*(b - d))/(epsilon + (a - c)^2 + (b - d)^2)^(1/2)))*(2*b - 2*det1y))/((exp(so - (epsilon + (a - c)^2 + (b - d)^2)^(1/2)) + 1)*(csig + (a - det1x)^2 + (b - det1y)^2)^2) + (csig*exp(so - (epsilon + (a - c)^2 + (b - d)^2)^(1/2))*(pi - acos((d1x*(a - c) + d1y*(b - d))/(epsilon + (a - c)^2 + (b - d)^2)^(1/2)))*(2*b - 2*d))/(2*(exp(so - (epsilon + (a - c)^2 + (b - d)^2)^(1/2)) + 1)^2*(csig + (a - det1x)^2 + (b - det1y)^2)*(epsilon + (a - c)^2 + (b - d)^2)^(1/2));
                        diff3 = (csig*(d2x/(epsilon + (c - e)^2 + (d - f)^2)^(1/2) - ((2*c - 2*e)*(d2x*(c - e) + d2y*(d - f)))/(2*(epsilon + (c - e)^2 + (d - f)^2)^(3/2))))/((1 - (d2x*(c - e) + d2y*(d - f))^2/(epsilon + (c - e)^2 + (d - f)^2))^(1/2)*(exp(so - (epsilon + (c - e)^2 + (d - f)^2)^(1/2)) + 1)*(csig + (c - det2x)^2 + (d - det2y)^2)) - (csig*(pi - acos((d2x*(c - e) + d2y*(d - f))/(epsilon + (c - e)^2 + (d - f)^2)^(1/2)))*(2*c - 2*det2x))/((exp(so - (epsilon + (c - e)^2 + (d - f)^2)^(1/2)) + 1)*(csig + (c - det2x)^2 + (d - det2y)^2)^2) - (csig*(d1x/(epsilon + (a - c)^2 + (b - d)^2)^(1/2) - ((2*a - 2*c)*(d1x*(a - c) + d1y*(b - d)))/(2*(epsilon + (a - c)^2 + (b - d)^2)^(3/2))))/((exp(so - (epsilon + (a - c)^2 + (b - d)^2)^(1/2)) + 1)*(1 - (d1x*(a - c) + d1y*(b - d))^2/(epsilon + (a - c)^2 + (b - d)^2))^(1/2)*(csig + (a - det1x)^2 + (b - det1y)^2)) + (csig*exp(so - (epsilon + (c - e)^2 + (d - f)^2)^(1/2))*(pi - acos((d2x*(c - e) + d2y*(d - f))/(epsilon + (c - e)^2 + (d - f)^2)^(1/2)))*(2*c - 2*e))/(2*(exp(so - (epsilon + (c - e)^2 + (d - f)^2)^(1/2)) + 1)^2*(csig + (c - det2x)^2 + (d - det2y)^2)*(epsilon + (c - e)^2 + (d - f)^2)^(1/2)) - (csig*exp(so - (epsilon + (a - c)^2 + (b - d)^2)^(1/2))*(pi - acos((d1x*(a - c) + d1y*(b - d))/(epsilon + (a - c)^2 + (b - d)^2)^(1/2)))*(2*a - 2*c))/(2*(exp(so - (epsilon + (a - c)^2 + (b - d)^2)^(1/2)) + 1)^2*(csig + (a - det1x)^2 + (b - det1y)^2)*(epsilon + (a - c)^2 + (b - d)^2)^(1/2));
                        diff4 = (csig*(d2y/(epsilon + (c - e)^2 + (d - f)^2)^(1/2) - ((2*d - 2*f)*(d2x*(c - e) + d2y*(d - f)))/(2*(epsilon + (c - e)^2 + (d - f)^2)^(3/2))))/((1 - (d2x*(c - e) + d2y*(d - f))^2/(epsilon + (c - e)^2 + (d - f)^2))^(1/2)*(exp(so - (epsilon + (c - e)^2 + (d - f)^2)^(1/2)) + 1)*(csig + (c - det2x)^2 + (d - det2y)^2)) - (csig*(pi - acos((d2x*(c - e) + d2y*(d - f))/(epsilon + (c - e)^2 + (d - f)^2)^(1/2)))*(2*d - 2*det2y))/((exp(so - (epsilon + (c - e)^2 + (d - f)^2)^(1/2)) + 1)*(csig + (c - det2x)^2 + (d - det2y)^2)^2) - (csig*(d1y/(epsilon + (a - c)^2 + (b - d)^2)^(1/2) - ((2*b - 2*d)*(d1x*(a - c) + d1y*(b - d)))/(2*(epsilon + (a - c)^2 + (b - d)^2)^(3/2))))/((exp(so - (epsilon + (a - c)^2 + (b - d)^2)^(1/2)) + 1)*(1 - (d1x*(a - c) + d1y*(b - d))^2/(epsilon + (a - c)^2 + (b - d)^2))^(1/2)*(csig + (a - det1x)^2 + (b - det1y)^2)) + (csig*exp(so - (epsilon + (c - e)^2 + (d - f)^2)^(1/2))*(pi - acos((d2x*(c - e) + d2y*(d - f))/(epsilon + (c - e)^2 + (d - f)^2)^(1/2)))*(2*d - 2*f))/(2*(exp(so - (epsilon + (c - e)^2 + (d - f)^2)^(1/2)) + 1)^2*(csig + (c - det2x)^2 + (d - det2y)^2)*(epsilon + (c - e)^2 + (d - f)^2)^(1/2)) - (csig*exp(so - (epsilon + (a - c)^2 + (b - d)^2)^(1/2))*(pi - acos((d1x*(a - c) + d1y*(b - d))/(epsilon + (a - c)^2 + (b - d)^2)^(1/2)))*(2*b - 2*d))/(2*(exp(so - (epsilon + (a - c)^2 + (b - d)^2)^(1/2)) + 1)^2*(csig + (a - det1x)^2 + (b - det1y)^2)*(epsilon + (a - c)^2 + (b - d)^2)^(1/2));
                        diff5 = - (csig*(d2x/(epsilon + (c - e)^2 + (d - f)^2)^(1/2) - ((2*c - 2*e)*(d2x*(c - e) + d2y*(d - f)))/(2*(epsilon + (c - e)^2 + (d - f)^2)^(3/2))))/((1 - (d2x*(c - e) + d2y*(d - f))^2/(epsilon + (c - e)^2 + (d - f)^2))^(1/2)*(exp(so - (epsilon + (c - e)^2 + (d - f)^2)^(1/2)) + 1)*(csig + (c - det2x)^2 + (d - det2y)^2)) - (csig*exp(so - (epsilon + (c - e)^2 + (d - f)^2)^(1/2))*(pi - acos((d2x*(c - e) + d2y*(d - f))/(epsilon + (c - e)^2 + (d - f)^2)^(1/2)))*(2*c - 2*e))/(2*(exp(so - (epsilon + (c - e)^2 + (d - f)^2)^(1/2)) + 1)^2*(csig + (c - det2x)^2 + (d - det2y)^2)*(epsilon + (c - e)^2 + (d - f)^2)^(1/2));
                        diff6 = - (csig*(d2y/(epsilon + (c - e)^2 + (d - f)^2)^(1/2) - ((2*d - 2*f)*(d2x*(c - e) + d2y*(d - f)))/(2*(epsilon + (c - e)^2 + (d - f)^2)^(3/2))))/((1 - (d2x*(c - e) + d2y*(d - f))^2/(epsilon + (c - e)^2 + (d - f)^2))^(1/2)*(exp(so - (epsilon + (c - e)^2 + (d - f)^2)^(1/2)) + 1)*(csig + (c - det2x)^2 + (d - det2y)^2)) - (csig*exp(so - (epsilon + (c - e)^2 + (d - f)^2)^(1/2))*(pi - acos((d2x*(c - e) + d2y*(d - f))/(epsilon + (c - e)^2 + (d - f)^2)^(1/2)))*(2*d - 2*f))/(2*(exp(so - (epsilon + (c - e)^2 + (d - f)^2)^(1/2)) + 1)^2*(csig + (c - det2x)^2 + (d - det2y)^2)*(epsilon + (c - e)^2 + (d - f)^2)^(1/2));
                        
                        dfx(xind-2)=dfx(xind-2) +   diff1;
                        dfx(yind-2)=dfx(yind-2) +   diff2;
                        dfx(xind)=dfx(xind)     +   diff3;
                        dfx(yind)=dfx(yind)     +   diff4;
                        dfx(xind+2)=dfx(xind+2) +   diff5;
                        dfx(yind+2)=dfx(yind+2) +   diff6;
                        
                    end
%                                         dfx(isnan(dfx))=1e-5*rand;
%                                         dfx(isinf(dfx))=1e-5*rand;
                    if ~isempty(find(~isfinite(dfx)))
                        [i t]
                        (v1(1)*d1x + v1(2)*d1y) / m1
                        (v2(1)*d2x + v2(2)*d2y) / m2
                        f1
                        f2
                        f1acos
                        f2acos
                        d1
                        d2
                        v1
                        v2
                        diffterm
                        [diff1 diff2 diff3 diff4 diff5 diff6]
                        [a b c d e f]
                        t1
                        t2
                        t3

                        
%                         pause
                    end
                end
            end
        end
        
        xind=xind+2;
        yind=xind+1;
        
    end
    if tlength>1
        xind=xind+2;
    end
end
% dfx'
% fx=fx/gridStep;		% normalize
% dfx=dfx/gridStep;	% normalize
% ds=ds/gridStep;
end