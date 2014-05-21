function [v tvis occx occy ddvix ddviy]=computeOcclusions2d(X,Y)

% NOT YET IMPLEMENTED !!!
%
% compute all mutual occlusions between all targets (on image)
%
% input
% X,Y is the current state vector in matrix form
%
% output
% v     - a FxN matrix, containing the rho-adjusted visibility portion of each target in each frame
% tvis  - same as v, but not adjusted with rho. Could potentially be > 1
% occx  - x-derivative of each v
% occy  - y-derivative of each v
% ddvix - x-derivative of visibility between each pair of targets
% ddviy - y-derivative of visibility between each pair of targets
% 
% 

% global N F;
[F N]=size(X);
global sceneInfo opt;

htobj=sceneInfo.htobj;

rho=1;

% convert all to image
Xi=X; Yi=Y;

extar=find(Xi);
heights=feval(htobj,Xi(extar),Yi(extar));

Z=zeros(size(X));
Z(extar)=heights);
Xi
Z
pause
% Z=0*Z;

% allCT=zeros(2,2,F,N);
% for t=1:F
%     existobj=find(X(t,:));
%     for i=existobj
%         xw=X(t,i); yw=Y(t,i);
%                 [Xi2(t,i) Yi2(t,i)]=worldToImage(xw,yw,zw,mR,mT,camPar.mInt,camPar.mGeo);
% %         [Xi(t,i) Yi(t,i)]=worldToImageNRD(xw,yw,zw,mR,mT,camPar.mInt,camPar.mGeo);
%     end
% end

% [Xi Yi ]=worldToImageNRD_mex(X,Y,Z, ...
%     camPar.mGeo.mDpx, camPar.mGeo.mDpy, ...
%     camPar.mInt.mSx, camPar.mInt.mCx, camPar.mInt.mCy, camPar.mInt.mFocal, camPar.mInt.mKappa1, ...
%     mR,mT);

[Xi Yi]=allWorldToImage_mex(X,Y,Z, ...
    camPar.mGeo.mDpx, camPar.mGeo.mDpy, ...
    camPar.mInt.mSx, camPar.mInt.mCx, camPar.mInt.mCy, camPar.mInt.mFocal, camPar.mInt.mKappa1,...
    mR,mT);

%%% !!!!!!!!!! TAKE CARE OF IMAGE SCALE!!!
Xi=Xi/scaleX; Yi=Yi/scaleY;
% Xit=Xit/scaleX; Yit=Yit/scaleY;


% Xi=X;
% Yi=Y;
% sum(sum(abs([Xi;Yi]-[Xi2;Yi2])))
ol=zeros(N);
v=zeros(F,N);
tvis=zeros(F,N);
occx=zeros(F,N);
occy=zeros(F,N);
ddvix=zeros(F,N,N);
ddviy=zeros(F,N,N);
sf=zeros(N,N);
% size(X)
% p.mex=0;
if opt.mex
%     if cexperiment==11 || cexperiment==12 || cexperiment == 13
%         [v tvis occx occy ddvix ddviy]=computeOcclusions2new_mex(X,Y,Xi,Yi,Dpx,Dpy,Sx,focal,mR11,mR12,mR13,mR21,mR22,mR23,mR31,mR32,mR33,tx,ty,tz,txi,tyi,tzi,rho);
%     else
%         [v tvis occx occy ddvix ddviy]=computeOcclusions2_mex(X,Y,Xi,Yi,Dpx,Dpy,Sx,focal,mR11,mR12,mR13,mR21,mR22,mR23,mR31,mR32,mR33,tx,ty,tz,rho);
%     end
    [v tvis occx occy ddvix ddviy]=computeOcclusions2new_mex(X,Y,Xi,Yi,Dpx,Dpy,Sx,focal,mR11,mR12,mR13,mR21,mR22,mR23,mR31,mR32,mR33,tx,ty,tz,txi,tyi,tzi,rho);
else
    for t=1:F
        exobj=find(X(t,:));
        
        icnt=0;
        occmat=zeros(N);
        occmatdX=zeros(N);
        occmatdY=zeros(N);
        occmatdX2=zeros(N);
        occmatdY2=zeros(N);
        ol=zeros(N);
        
        
        for i=exobj
            
            icnt=icnt+1;
            
            a=[Xi(t,i) Yi(t,i)];
            Xwa=X(t,i); Ywa=Y(t,i);
            Zwa=900; Zwb=900;
            
            
%             if cexperiment==11 || cexperiment==12 || cexperiment == 13
%                 Za=sqrt((txi-Xwa)^2 + (tyi-Ywa)^2 +(tzi-Zwa)^2);
%             else
%                 Za=sqrt((tx-Xwa)^2 + (ty-Ywa)^2 +(tz-Zwa)^2);
%             end
            Za=sqrt((txi-Xwa)^2 + (tyi-Ywa)^2 +(tzi-Zwa)^2);
            sizeOnImagea=1800*focal/Za/Dpy;
%             figure(2)
%             clf
%             global frameNums;
%             [frameNums(t) t i]
%             im=imread(getImageFile(scenario,1,frameNums(t)));
%             aa=round(a);
%             x1=max(1,aa(1)-30);
%             x2=min(768,aa(1)+30);
%             x2=min(640,aa(1)+30);
%             
%             a
%             sizeOnImagea            
%             y2=aa(2); y1=y2-round(sizeOnImagea);
%             [x1 x2 y1 y2]
%             im=im(y1:y2,x1:x2,:);
%             imshow(im);
% 
%             pause
            A=[(sizeOnImagea/2)^2/4 0; 0 (sizeOnImagea/2)^2];
            
            for j=exobj
                if j>i
                    Xwb=X(t,j); Ywb=Y(t,j);
                    b=[Xi(t,j) Yi(t,j)];
                    
%                     if cexperiment==11 || cexperiment==12
%                         Zb=sqrt((txi-Xwb)^2 + (tyi-Ywb)^2 +(tzi-Zwb)^2);
%                     else
%                         Zb=sqrt((tx-Xwb)^2 + (ty-Ywb)^2 +(tz-Zwb)^2);
%                     end
                    Zb=sqrt((txi-Xwb)^2 + (tyi-Ywb)^2 +(tzi-Zwb)^2);
                    sizeOnImageb=1800*focal/Zb/Dpy;
                    
                    B=[(sizeOnImageb/2)^2/4 0; 0 (sizeOnImageb/2)^2];
                    
                    C=A+B;
                    alessb=a-b;
                    ol(i,j)=exp(-.5*(alessb * (C \ alessb')));
                end
            end
            
            
            for j=exobj
                if j<i
                    ol(i,j)=ol(j,i);
                end
                
                if i~=j
                    Xwb=X(t,j); Ywb=Y(t,j);
                    
                    sigfac=1/(1+exp(Yi(t,i)-Yi(t,j)));
                    sf(i,j)=sigfac;
                    
                    occmat(i,j)=ol(i,j)*sigfac;
                    
                    if ol(i,j)>0.0001 && nargout > 1 %% source for inexact gradient !!!
%                                                             occmatdX(i,j)=dfdXwa(Dpx,Dpy,Sx,Xwa,Xwb,Ywa,Ywb,Zwa,Zwb,focal,kappa,mR11,mR12,mR13,mR21,mR22,mR23,mR31,mR32,mR33,tx,ty,tz);
%                                                             occmatdY(i,j)=dfdYwa(Dpx,Dpy,Sx,Xwa,Xwb,Ywa,Ywb,Zwa,Zwb,focal,kappa,mR11,mR12,mR13,mR21,mR22,mR23,mR31,mR32,mR33,tx,ty,tz);
%                                                             occmatdX2(i,j)=dfdXwb(Dpx,Dpy,Sx,Xwa,Xwb,Ywa,Ywb,Zwa,Zwb,focal,kappa,mR11,mR12,mR13,mR21,mR22,mR23,mR31,mR32,mR33,tx,ty,tz);
%                                                             occmatdY2(i,j)=dfdYwb(Dpx,Dpy,Sx,Xwa,Xwb,Ywa,Ywb,Zwa,Zwb,focal,kappa,mR11,mR12,mR13,mR21,mR22,mR23,mR31,mR32,mR33,tx,ty,tz);
%                                             occmatdX(i,j)=dfnddXwa(Dpx,Dpy,Sx,Xwa,Xwb,Ywa,Ywb,Zwa,Zwb,focal,mR11,mR12,mR13,mR21,mR22,mR23,mR31,mR32,mR33,tx,ty,tz);
%                                             occmatdY(i,j)=dfnddYwa(Dpx,Dpy,Sx,Xwa,Xwb,Ywa,Ywb,Zwa,Zwb,focal,mR11,mR12,mR13,mR21,mR22,mR23,mR31,mR32,mR33,tx,ty,tz);
%                                             occmatdX2(i,j)=dfnddXwb(Dpx,Dpy,Sx,Xwa,Xwb,Ywa,Ywb,Zwa,Zwb,focal,mR11,mR12,mR13,mR21,mR22,mR23,mR31,mR32,mR33,tx,ty,tz);
%                                             occmatdY2(i,j)=dfnddYwb(Dpx,Dpy,Sx,Xwa,Xwb,Ywa,Ywb,Zwa,Zwb,focal,mR11,mR12,mR13,mR21,mR22,mR23,mR31,mR32,mR33,tx,ty,tz);
%                                             occmatdX(i,j)=dfnddXwa_new(Dpx,Dpy,Sx,Xwa,Xwb,Ywa,Ywb,Zwa,Zwb,focal,mR11,mR12,mR13,mR21,mR22,mR23,mR31,mR32,mR33,tx,txi,ty,tyi,tz,tzi);
%                                             occmatdY(i,j)=dfnddYwa_new(Dpx,Dpy,Sx,Xwa,Xwb,Ywa,Ywb,Zwa,Zwb,focal,mR11,mR12,mR13,mR21,mR22,mR23,mR31,mR32,mR33,tx,txi,ty,tyi,tz,tzi);
%                                             occmatdX2(i,j)=dfnddXwb_new(Dpx,Dpy,Sx,Xwa,Xwb,Ywa,Ywb,Zwa,Zwb,focal,mR11,mR12,mR13,mR21,mR22,mR23,mR31,mR32,mR33,tx,txi,ty,tyi,tz,tzi);
%                                             occmatdY2(i,j)=dfnddYwb_new(Dpx,Dpy,Sx,Xwa,Xwb,Ywa,Ywb,Zwa,Zwb,focal,mR11,mR12,mR13,mR21,mR22,mR23,mR31,mR32,mR33,tx,txi,ty,tyi,tz,tzi);
                        %                     [dfndXa dfndYa dfndXb dfndYb]=dfnddALL(Dpx,Dpy,Sx,Xwa,Xwb,Ywa,Ywb,Zwa,Zwb,focal,mR11,mR12,mR13,mR21,mR22,mR23,mR31,mR32,mR33,tx,ty,tz);
                        
                        [dfndXa dfndYa dfndXb dfndYb]=dfnddALL_mex(Dpx,Dpy,Sx,Xwa,Xwb,Ywa,Ywb,Zwa,Zwb,focal,mR11,mR12,mR13,mR21,mR22,mR23,mR31,mR32,mR33,tx,ty,tz);
%                         [dfndXa dfndYa dfndXb dfndYb]=dfnddALLnew_mex(Dpx,Dpy,Sx,Xwa,Xwb,Ywa,Ywb,Zwa,Zwb,focal,mR11,mR12,mR13,mR21,mR22,mR23,mR31,mR32,mR33,tx,ty,tz,txi,tyi,tzi);
                        occmatdX(i,j)=dfndXa;
                        occmatdY(i,j)=dfndYa;
                        occmatdX2(i,j)=dfndXb;
                        occmatdY2(i,j)=dfndYb;
                        
                    end
                    
                end
            end
            
            
        end
        
        if nargout > 2
            for i=exobj
                sumdx=0;
                sumdy=0;
                for j=exobj
                    
                    if j~=i
                        sumdx=sumdx-rho*exp(-rho*sum(occmat(j,:)))*(occmatdX2(j,i));
                        ddvix(t,i,j) = -rho*exp(-rho*sum(occmat(j,:)))*(occmatdX2(j,i));
                        
                        sumdy=sumdy-rho*exp(-rho*sum(occmat(j,:)))*(occmatdY2(j,i));
                        ddviy(t,i,j) = -rho*exp(-rho*sum(occmat(j,:)))*(occmatdY2(j,i));
                    else
                        sumdx=sumdx-rho*exp(-rho*sum(occmat(i,:)))*(sum(occmatdX(i,:)));
                        ddvix(t,i,j) = -rho*exp(-rho*sum(occmat(i,:)))*(sum(occmatdX(i,:)));
                        
                        sumdy=sumdy-rho*exp(-rho*sum(occmat(i,:)))*(sum(occmatdY(i,:)));
                        ddviy(t,i,j) = -rho*exp(-rho*sum(occmat(i,:)))*(sum(occmatdY(i,:)));
                    end
                end
                occx(t,i)=sumdx;
                occy(t,i)=sumdy;
                %             if t==15
                %                 [t i]
                %                 [sumdx sumdy]
                %             end
            end
        end
        v(t,:)=exp(-rho*sum(occmat,2)');
        tvis(t,:)=exp(-sum(occmat,2)');
        %     t
        %     [occmat 1-sum(occmat,2) exp(-rho*sum(occmat,2))]
        %     pause
    end
    
    
end
% a=sum(sum(abs(v-v2)));
% b=sum(sum(abs(tvis-tvis2)));
% c=sum(sum(abs(occx-occx2)));
% % occx
% % occx2
% d=sum(sum(abs(occy-occy2)));
% e=sum(sum(abs(sum(ddvix-ddvix2))));
% f=sum(sum(abs(sum(ddviy-ddviy2))));
% % [a b c d e f]
% if sum([a b c d e f]) > 1e-5
%     [a b c d e f]
%     pause
% end
% sf
end
