%%
gtInfom=gtInfo;
load('/storage/databases/PNNL/ParkingLot/GT_PNNL_ParkingLot.mat')
a=1:3:stateInfo.frameNums(end);
a=find(a==stateInfo.frameNums(1)):find(a==stateInfo.frameNums(end)); a=a';
gtInfo.X=gtInfo.X(a,:);gtInfo.Y=gtInfo.Y(a,:);gtInfo.W=gtInfo.W(a,:);gtInfo.H=gtInfo.H(a,:);
gtInfo.Xi=gtInfo.Xi(a,:);gtInfo.Yi=gtInfo.Yi(a,:);
gtInfo.frameNums=stateInfo.frameNums;
gtInfo=cleanGT(gtInfo);
printFinalEvaluation(stateInfo,gtInfo,sceneInfo,opt)
[metrics2d_ metrics3d_]=printFinalEvaluation(stateInfo,gtInfo,sceneInfo,opt);
% gtInfo=gtInfom;