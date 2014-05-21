%%
function [f df ds]=Eapp(xv,stateInfo,vis,ddvix,ddviy)

f=0;

global scenario opt sceneInfo
global allhist;


nbins=opt.app.nbins;
xstrade=opt.app.xstrade;
tstrade=opt.app.tstrade;


camPar=sceneInfo.camPar;
N=stateInfo.N;
F=stateInfo.F;
targetsExist=stateInfo.targetsExist;
tiToInd=stateInfo.tiToInd;

% global xtil
% cam=1;
% camPar=getCameraParameters(getCameraCalibrationFile(scenario,cam));
[mR mT]=getRotTrans(camPar);
imSizeY =sceneInfo.imgHeight;
imSizeX =sceneInfo.imgWidth;
scaleX=camPar.mGeo.mImgWidth/imSizeX;
scaleY=camPar.mGeo.mImgHeight/imSizeY;


df=zeros(F,N,3);
% convert state vector to matrix representation
[X Y]=vectorToMatrices(xv,stateInfo);

% convert all to image
zw=850;
% zw=1275;
% zw=0;
Z=zw*ones(size(X));

[Xi Yi allCT]=allWorldToImageWithImDer_mex(X,Y,Z, ...
    camPar.mGeo.mDpx, camPar.mGeo.mDpy, ...
    camPar.mInt.mSx, camPar.mInt.mCx, camPar.mInt.mCy, camPar.mInt.mFocal, camPar.mInt.mKappa1, ...
    mR,mT);

[Xiz Yiz]=allWorldToImage_mex(X,Y,zeros(size(X)), ...
    camPar.mGeo.mDpx, camPar.mGeo.mDpy, ...
    camPar.mInt.mSx, camPar.mInt.mCx, camPar.mInt.mCy, camPar.mInt.mFocal, camPar.mInt.mKappa1, ...
    mR,mT);
[Xit Yit]=allWorldToImage_mex(X,Y,1700*ones(size(X)), ...
    camPar.mGeo.mDpx, camPar.mGeo.mDpy, ...
    camPar.mInt.mSx, camPar.mInt.mCx, camPar.mInt.mCy, camPar.mInt.mFocal, camPar.mInt.mKappa1, ...
    mR,mT);

% allCT(:,:,end,end)
%%% !!!!!!!!!! TAKE CARE OF IMAGE SCALE!!!
Xi=Xi/scaleX; Yi=Yi/scaleY;
Xiz=Xiz/scaleX; Yiz=Yiz/scaleY;
Xit=Xit/scaleX; Yit=Yit/scaleY;
allboxheights=70*ones(1,N);
for i=1:N
    boxHeight=mean(Yiz(targetsExist(i,1):targetsExist(i,2),i)-Yit(targetsExist(i,1):targetsExist(i,2),i));
    allboxheights(i)=boxHeight;
end
% Xi(4,1)
% vis

rho=1;
truevis=1+log(vis)/rho;
truevis(truevis<0)=0;

allvismeans=ones(F,N);
if opt.occ
    for i=1:N
        for t=targetsExist(i,1):tstrade:targetsExist(i,2)-tstrade
            allvismeans(t,i)=vismean(truevis(t,i),truevis(t+tstrade,i));
        end
    end
end


% bin center points
m=1:2:2*nbins;
m=m/(2*nbins);
minval=-10; maxval=10;
minval=0; maxval=1;
valran=maxval-minval;
m=m*valran+minval;
d=1/nbins;
SIG=150*eye(2); SIG(1)=750;
%%
SIG=50*eye(2); SIG(1)=500;
% SIG=10*eye(2); SIG(1)=100;
if scenario==42, SIG=4*SIG; end
% if cexperiment==74, SIG=2*SIG; end
% SIG=200*eye(2); SIG(1)=2000;
% SIG=40*eye(2); SIG(1)=160;

sa=SIG(1,1);
sb=SIG(2,1);
sc=SIG(1,2);
sd=SIG(2,2);
detSIG=sa*sd - sb*sc;
oneOverDetSIG=1/detSIG;
SIGINV=oneOverDetSIG * [SIG(4) -SIG(3); -SIG(2) SIG(1)];
sia=SIGINV(1); sid=SIGINV(4);
twodetSIG=(2*detSIG);
twoPiSqrtDet=(2*pi*sqrt(detSIG));
oneOverTwoPiSqrtDet=1/twoPiSqrtDet;
normfac=1/(2*pi*sqrt(detSIG));

sigshift=5.564;  sigscale=40.73;
sigshift=7.235;  sigscale=33.71;

%%
%%% CAREFUL! THIS IS ALL IN MATRIX COORDINATES:
%    y  --->
%  x
% ||
% ||
% \/
%%%%%%%%%%%%%%%%%%%%%%%%%
% winradX=60;
% winradY=30;
%                 k=rand(1,nbins);
%                 k=k/sum(k);

channels=1:3;
if length(channels)<3; fprintf('channels = %i \n',channels); end

nch=length(channels);
% channels=3;
fchan=zeros(1,max(channels));
dfchan=zeros(length(xv),max(channels));
ds=zeros(F,N,max(channels));
% ds=zeros(size(X));

global ds1 ds2 ds3 dsm allvm;
ds1=zeros(length(xv)/2,max(channels));
ds2=zeros(F,N,max(channels));
ds3=zeros(F,N,max(channels));
dsm=zeros(F,N,max(channels));
allvm=zeros(F,N);


% global  winradX winradY  ;

dhx=zeros(1,nbins);
dkx=zeros(1,nbins);
dhy=zeros(1,nbins);
dky=zeros(1,nbins);


allBC=zeros(F-1,N,nch);
xstrade=xstrade; % integer!
ystrade=xstrade; % integer!
norm2=xstrade*ystrade;

% tstrade=2;
tempmatch=0;
showmasked=0;

%     dstep=0.0001;
%     minlimit=-20;
%     suppoints=minlimit:dstep:-0.001;
%     shiftright=numel(find(suppoints<=0));
%     shiftright=numel(suppoints);
%     expLUT=exp(suppoints);
    
for chan=channels
    cnt=0;
    for i=1:N
        
        if diff(targetsExist(i,:))>=tstrade
            h=zeros(1,nbins);
            hb=zeros(1,nbins);
            dhx=zeros(1,nbins);
            dhy=zeros(1,nbins);
            
            boxHeight=allboxheights(i);
            SIG=[(boxHeight/2)^2 0; 0 (boxHeight/2)^2/9]/8;
            sa=SIG(1,1); sb=SIG(2,1); sc=SIG(1,2); sd=SIG(2,2);
            detSIG=sa*sd - sb*sc;
            oneOverDetSIG=1/detSIG;
            SIGINV=oneOverDetSIG * [SIG(4) -SIG(3); -SIG(2) SIG(1)];

            twoPiSqrtDet=(2*pi*sqrt(detSIG));
            oneOverTwoPiSqrtDet=1/twoPiSqrtDet;
            normfac=1/(2*pi*sqrt(detSIG));
            sia=SIGINV(1); sid=SIGINV(4);
            winradX=round(boxHeight/4*3);
            winradY=round(winradX/3);
            

            for t=targetsExist(i,1):tstrade:targetsExist(i,2)-tstrade
                if t==targetsExist(i,1)
                    
                    x=Yi(t,i);
                    y=Xi(t,i);
                    %                 CT=imDerivatives(X(t,i),Y(t,i),zw,mR,mT,camPar);
                    CT=allCT(:,:,t,i);
                    
                    x1 = max(1,round(x)-winradX):xstrade:min(imSizeY,round(x)+winradX);
                    x2 = max(1,round(y)-winradY):ystrade:min(imSizeX,round(y)+winradY);
                    if showmasked
                                    II=double(imread(getImageFile(scenario,1,frameNums(t))))/255;
                                    II=II(x1,x2,:);
                    end
                    
                    xx = x1(:).'; % Make sure x is a full row vector.
                    yy = x2(:);   % Make sure y is a full column vector.
                    nx = length(xx); ny = length(yy);
                    X1 = xx(ones(ny, 1),:);
                    X2 = yy(:,ones(1, nx));
                    X1=X1';
                    X2=X2';
                    
                    %% speed up mvnpdf. only valid for matrices with offdiagonal = 0 !!!
                    xlessm=X1(:)-x;
                    ylessm=X2(:)-y;
                    
                    xlmsq=sia*xlessm.^2;
                    ylmsq=sid*ylessm.^2;
                    expon=-0.5*(xlmsq+ylmsq);
                    
                    expterm=exp(expon);

                    
                    expterm=reshape(expterm,length(x1),length(x2));
                    w=normfac * expterm * norm2;

                    
                    if nargout>1
                        ylessyy=-ylessm;
                        xlessxx=-xlessm;
                        
                        %                     dwx=-(exp((ylessyy.*((sc*xlessxx)/detSIG - (sa*ylessyy)/detSIG))/2 - (xlessxx.*((sd*xlessxx)/detSIG - (sb*ylessyy)/detSIG))/2).*((sb*ylessyy)/twodetSIG - (sd*xlessxx)/detSIG + (sc*ylessyy)/twodetSIG))/(2*pi*detSIG^(1/2));
                        %                     dwy=-(exp((ylessyy.*((sc*xlessxx)/detSIG - (sa*ylessyy)/detSIG))/2 - (xlessxx.*((sd*xlessxx)/detSIG - (sb*ylessyy)/detSIG))/2).*((sb*xlessxx)/twodetSIG + (sc*xlessxx)/twodetSIG - (sa*ylessyy)/detSIG))/(2*pi*detSIG^(1/2));
                        
                        dxlx=sd*xlessxx;
                        ayly=sa*ylessyy;
                        tmp2=oneOverDetSIG*ayly;
                        tmp3=0.5*(ylessyy.*(-tmp2)) - 0.5*(xlessxx.*(oneOverDetSIG*dxlx));
                        tmp5=exp(tmp3);
                        
                        
                        dwx=-oneOverTwoPiSqrtDet*(tmp5.*(- oneOverDetSIG*dxlx));
                        dwy=-oneOverTwoPiSqrtDet*(tmp5.*(- oneOverDetSIG*ayly));
                        
                        dwx=reshape(norm2*dwx,length(x1),length(x2));
                        dwy=reshape(norm2*dwy,length(x1),length(x2));                        
                        
                        Dwy=(CT(1)*dwy) + (CT(3)*dwx);
                        Dwx=(CT(2)*dwy) + (CT(4)*dwx);
                        %                     Dwx=dwx; Dwy=dwy;
                        
                        
                        
                    end
                    
                    
                    hh=allhist(x1,x2,chan,t);
                    
                    
                    h=gausshist(w,hh,nbins);
                    %                 hb=h;
                    hb=gausshist(1/numel(w)*ones(size(w)),hh,nbins);

                    if nargout>1
                        [dhx dhy]=gausshistdx(Dwx,Dwy,hh,nbins);
                    end
                else
                    k2=h;
                    kb2=hb;
                    
                    if ~tempmatch
                        h=k;
                        hb=kb;
                        dhx=dkx;
                        dhy=dky;
                    end
                                    if showmasked, II=IIn; end
                    
                    
                end

                
                
                x=Yi(t+tstrade,i);
                y=Xi(t+tstrade,i);
                CT=allCT(:,:,t+tstrade,i);
                
                x1 = max(1,round(x)-winradX):xstrade:min(imSizeY,round(x)+winradX);
                x2 = max(1,round(y)-winradY):ystrade:min(imSizeX,round(y)+winradY);
                if showmasked
                                IIn=double(imread(getImageFile(scenario,1,frameNums(t+tstrade))))/255;
                                IIn=IIn(x1,x2,:);
                end
                
                
                
                xx = x1(:).'; % Make sure x is a full row vector.
                yy = x2(:);   % Make sure y is a full column vector.
                nx = length(xx); ny = length(yy);
                X1 = xx(ones(ny, 1),:);
                X2 = yy(:,ones(1, nx));
                
                
                X1=X1';
                X2=X2';

                %
                %% speed up mvnpdf. only valid for matrices with offdiagonal = 0 !!!
                xlessm=X1(:)-x;
                ylessm=X2(:)-y;
                
                xlmsq=sia*xlessm.^2;
                ylmsq=sid*ylessm.^2;
                expon=-0.5*(xlmsq+ylmsq);
                expterm=exp(expon);
%                 exptermind=floor(expon/dstep+shiftright);
%                 exptermLUT=expLUT(exptermind)';
%                 norm(expterm-exptermLUT)
%                 pause
                
                wp=w;
                expterm=reshape(expterm,length(x1),length(x2));
                w=normfac * expterm * norm2;
                %             w=w/sum(sum(w));
                %             sum(sum(w))
%                             w=ones(size(w))/numel(w);
                

                if nargout>1
                    ylessyy=-ylessm;
                    xlessxx=-xlessm;
                    
                    dxlx=sd*xlessxx;
                    ayly=sa*ylessyy;                    
                    
                    tmp2=oneOverDetSIG*ayly;
                    
                    tmp3=0.5*(ylessyy.*(-tmp2)) - 0.5*(xlessxx.*(oneOverDetSIG*dxlx));
                    
                    tmp5=exp(tmp3);
                    
                    
                    dwx=-oneOverTwoPiSqrtDet*(tmp5.*(- oneOverDetSIG*dxlx));
                    dwy=-oneOverTwoPiSqrtDet*(tmp5.*(- oneOverDetSIG*ayly));
                    
                    
                    dwx=reshape(norm2*dwx,length(x1),length(x2));
                    dwy=reshape(norm2*dwy,length(x1),length(x2));
                    %                 dwx=reshape(norm2*dwxaut,length(x1),length(x2));
                    %                 dwy=reshape(norm2*dwyaut,length(x1),length(x2));
                    
                    
                    
                    Dw2y=(CT(1)*dwy) + (CT(3)*dwx);
                    Dw2x=(CT(2)*dwy) + (CT(4)*dwx);
                end
                
                
                hh=allhist(x1,x2,chan,t+tstrade);
                k=gausshist(w,hh,nbins);
                kb=gausshist(1/numel(w)*ones(size(w)),hh,nbins);
                if nargout>1
                    [dkx dky]=gausshistdx(Dw2x,Dw2y,hh,nbins);
                end
                
                
                
                
                
                curbc=1-sum(sqrt(h.*k));
                %             curbc=0.5-sum(sqrt(h.*k));
                
                %             curbc=-sum(sqrt(vt*vtt*h.*k));
                %             curbc=-1;
                allBC(t,i,chan)=curbc;
                mv=allvismeans(t,i);
                %             mv=1;
                %             mv
                
                %             title(sprintf('%f %f',curbc,mv*curbc));
                %             pause
                mvc=mv*curbc;
                fac=mvc;
%                 fac=fac*fac;
                fac=1/(1+exp(sigshift-sigscale*mvc));
%                 if ~isempty(intersect(cexperiment,[78,86,79,87]))
%                     fac=1/(1+exp(sigshift-sigscale*mvc));
%                 end
                    %             fac=1/(1+exp(10-200*fac));
                fchan(chan)=fchan(chan)+fac;
                
                if showmasked
                    figure(1);
                    clf;
                    subplot(1,4,1); 
                    imshow(II);title(sprintf('BC %.2f',curbc));
                    subplot(1,4,2); 
                    w3p=repmat(w/max(max(wp)),[1 1 3]);
                    w3=repmat(w/max(max(w)),[1 1 3]);
                    IImasked=II.*w3p;
                    imshow(IImasked);title(sprintf('mvc %.2f',mvc));
                    
                    subplot(1,4,3);
                    imshow(IIn);title(sprintf('sigm %.2f',fac));
                    subplot(1,4,4);
                    
                    
                    IImasked=IIn.*w3;
                    imshow(IImasked);
                    pause
                end
                
                %             disp([chan t i mv curbc curbc*mv])
                cnt=cnt+1;
                ds(t,i,chan)=fac;
                ds1(cnt,chan)=sqrt(1-sum(sqrt(h.*k)));
                ds2(t,i,chan)=curbc;
                %             ds2(t,i,chan)=curbc;
                %             ds2(t,i,chan)=mv*curbc;
                %             ds3(cnt,chan)=sum(((h-k).^2)./(h+k));
                ds3(t,i,chan)=mvc;
                %             dsm(t,i,chan)=abs(mean(mean(II(:,:,chan)))-mean(mean(IIn(:,:,chan))));
                allvm(t,i)=mv;
                %         sqrt(1-sum(sqrt(h.*k)))
                
                if ~tempmatch
                    %%%%%%%%%%%%%
                    % dh
                    h(h==0)=1e-5; %%% !!!!
                    k(k==0)=1e-5; %%% !!!!
                    term1x=-sum((sqrt(k).*dhx)./(2*sqrt(h)));
                    term1y=-sum((sqrt(k).*dhy)./(2*sqrt(h)));
                    term2=1/(2*sqrt(1-sum(sqrt(h.*k))));
                    term2=1;
                    dfx=term1y*term2; % horizontal
                    dfy=term1x*term2; % vert
                    
                    %             dfx=-sum((vt*dhx.*sqrt(vtt*k))./(2*sqrt(vt*h)));
                    
                    %             dmvx=sum(ddvix(t,:,i));
                    dfx=mv*dfx;
                    dfy=mv*dfy;
                    
%                     dfx=dfx*2*mvc;
%                     dfy=dfy*2*mvc;
                    
                    dfx=dfx*(sigscale*exp(sigshift - sigscale*mvc))/(exp(sigshift - sigscale*mvc) + 1)^2;                 
                    dfy=dfy*(sigscale*exp(sigshift - sigscale*mvc))/(exp(sigshift - sigscale*mvc) + 1)^2;                 
                    
                    dfchan(tiToInd(t,i),chan)=dfchan(tiToInd(t,i),chan)+dfx;
                    dfchan(tiToInd(t,i)+1,chan)=dfchan(tiToInd(t,i)+1,chan)+dfy;
                    
                    %% derivative past frame!
                    if t>targetsExist(i,1)
                        % dh
                        
                        k2(k2==0)=1e-5; %%% !!!!
                        term1x=-sum((sqrt(k2).*dhx)./(2*sqrt(h)));
                        term1y=-sum((sqrt(k2).*dhy)./(2*sqrt(h)));
                        term2=1/(2*sqrt(1-sum(sqrt(h.*k2))));
                        term2=1;
                        dfx=term1y*term2; % horizontal
                        dfy=term1x*term2; % vert
                        
                        mv=allvismeans(t-tstrade,i);
                        mvc=mv*allBC(t-tstrade,i,chan);
                        
                        dfx=mv*dfx;
                        dfy=mv*dfy;
                        
%                         dfx=dfx*2*mvc;
%                         dfy=dfy*2*mvc;
                        dfx=dfx*(sigscale*exp(sigshift - sigscale*mvc))/(exp(sigshift - sigscale*mvc) + 1)^2;
                        dfy=dfy*(sigscale*exp(sigshift - sigscale*mvc))/(exp(sigshift - sigscale*mvc) + 1)^2;
                        
                        %             df(tiToInd(t-1,i))=dfx;
                        %             df(tiToInd(t-1,i)+1)=dfy;
                        %                 [dfx dfy mv]
                        dfchan(tiToInd(t,i),chan)=dfchan(tiToInd(t,i),chan)+dfx;
                        dfchan(tiToInd(t,i)+1,chan)=dfchan(tiToInd(t,i)+1,chan)+dfy;
                    end
                    
                else
                    h(h==0)=1e-5; %%% !!!!
                    k(k==0)=1e-5; %%% !!!!
                    term1x=-sum((sqrt(k).*dhx)./(2*sqrt(h)));
                    term1y=-sum((sqrt(k).*dhy)./(2*sqrt(h)));
                    term2=1/(2*sqrt(1-sum(sqrt(h.*k))));
                    term2=1;
                    dfx=term1y*term2; % horizontal
                    dfy=term1x*term2; % vert
                    
                    dfx=mv*dfx;    dfy=mv*dfy;
                    
                    dfchan(tiToInd(targetsExist(i,1),i),chan)=dfchan(tiToInd(targetsExist(i,1),i),chan)+dfx;
                    dfchan(tiToInd(targetsExist(i,1),i)+1,chan)=dfchan(tiToInd(targetsExist(i,1),i)+1,chan)+dfy;
                    
                    term1x=-sum((sqrt(h).*dkx)./(2*sqrt(k)));
                    term1y=-sum((sqrt(h).*dky)./(2*sqrt(k)));
                    term2=1/(2*sqrt(1-sum(sqrt(k.*h))));
                    term2=1;
                    dfx=term1y*term2; % horizontal
                    dfy=term1x*term2; % vert
                    
                    mv=allvismeans(t,i);
                    
                    dfx=mv*dfx;           dfy=mv*dfy;
                    
                    dfchan(tiToInd(t+tstrade,i),chan)=dfchan(tiToInd(t+tstrade,i),chan)+dfx;
                    dfchan(tiToInd(t+tstrade,i)+1,chan)=dfchan(tiToInd(t+tstrade,i)+1,chan)+dfy;
                    
                end
                %             if t==targetsExist(i,2)-tstrade
                %                 cnt=cnt+1; ds(cnt,chan)=0;
                %             end
            end
            
            %% last frame
            lf=targetsExist(i,1):tstrade:targetsExist(i,2);
            lf=lf(end);
            kt=h;
            h=k;
            dhx=dkx;
            dhy=dky;
            k=kt;
            
            % dh
            h(h==0)=1e-5; %%% !!!!
            k(k==0)=1e-5; %%% !!!!
            term1x=-sum((sqrt(k).*dhx)./(2*sqrt(h)));
            term1y=-sum((sqrt(k).*dhy)./(2*sqrt(h)));
            term2=1/(2*sqrt(1-sum(sqrt(h.*k))));
            term2=1;
            dfx=term1y*term2; % horizontal
            dfy=term1x*term2; % vert
            
            mv=allvismeans(lf-tstrade,i);
            mvc=mv*allBC(lf-tstrade,i,chan);
            %         mv=1;
            
            dfx=mv*dfx;
            dfy=mv*dfy;
            
%             dfx=dfx*2*mvc;
%             dfy=dfy*2*mvc;
            dfx=dfx*(sigscale*exp(sigshift - sigscale*mvc))/(exp(sigshift - sigscale*mvc) + 1)^2;
            dfy=dfy*(sigscale*exp(sigshift - sigscale*mvc))/(exp(sigshift - sigscale*mvc) + 1)^2;
            
            dfchan(tiToInd(lf,i),chan)=dfx;
            dfchan(tiToInd(lf,i)+1,chan)=dfy;
        end
    end
    
    %     if nargout>1
    %         for k=1:N
    %             for p=targetsExist(k,1):targetsExist(k,2)-1
    %                 exobj=find(X(p,:));
    %                 dx=0;
    %                 dy=0;
    %                 for i=exobj
    %                     dx=dx+ddvix(p,i,k)*allBC(p,i,chan);
    %                     dy=dy+ddviy(p,i,k)*allBC(p,i,chan);
    %                 end
    %
    %                 dfchan(tiToInd(p,k),chan)=dfchan(tiToInd(p,k),chan)+dx;
    %                 dfchan(tiToInd(p,k)+1,chan)=dfchan(tiToInd(p,k)+1,chan)+dy;
    %             end
    %         end
    %     end
    
end



f=mean(fchan(channels));
df=mean(dfchan(:,channels),2);


% df(:)=0;
% size(df)
% size(ds)
% ds
ds=mean(ds(:,:,channels),3);
ds1=mean(ds1(:,channels),2);
ds2=mean(ds2(:,:,channels),3);
ds3=mean(ds3(:,:,channels),3);
dsm=mean(dsm(:,:,channels),3);

end

function ret=H(I,x,y,mk,d)
ret=abs(I(x,y)-mk)<=d/2;
end

function ret=HH(I,mk,d)
ret=abs(I-mk)<=d/2;
end

function ret=vismean(vt,vtt)
%     ret=0.5*(vt+vtt); % arithmetic
ret=sqrt(vt*vtt); % geometric
%     ret=1;
end

% function ret=H(I,x,mk,d)
% ret=0;
% if abs(I(x)-mk)<d/2
%     ret=1;
% end

