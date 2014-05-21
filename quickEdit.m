%%
currun=bestruns(scenario,useexp);
keepids=setdiff(1:infos(scenario,currun).stateInfo.N, 91);
infos(scenario,currun).stateInfo.X=infos(scenario,currun).stateInfo.X(:,keepids);
infos(scenario,currun).stateInfo.Y=infos(scenario,currun).stateInfo.Y(:,keepids);
infos(scenario,currun).stateInfo.Xgp=infos(scenario,currun).stateInfo.Xgp(:,keepids);
infos(scenario,currun).stateInfo.Ygp=infos(scenario,currun).stateInfo.Ygp(:,keepids);
infos(scenario,currun).stateInfo.Xi=infos(scenario,currun).stateInfo.Xi(:,keepids);
infos(scenario,currun).stateInfo.Yi=infos(scenario,currun).stateInfo.Yi(:,keepids);
infos(scenario,currun).stateInfo.H=infos(scenario,currun).stateInfo.H(:,keepids);
infos(scenario,currun).stateInfo.W=infos(scenario,currun).stateInfo.W(:,keepids);

infos(scenario,currun).stateInfo.targetsExist=getTracksLifeSpans(infos(scenario,currun).stateInfo.X);
infos(scenario,currun).stateInfo.N=length(keepids);
infos(scenario,currun).stateInfo= ...
    matricesToVector(infos(scenario,currun).stateInfo.X,infos(scenario,currun).stateInfo.Y,infos(scenario,currun).stateInfo);

[metrics2d metrics3d]=printFinalEvaluation(infos(scenario,currun).stateInfo, gtInfo, sceneInfo, opt);

mets2d(scenario,:,currun)=metrics2d;
mets3d(scenario,:,currun)=metrics3d;

%%
currun=bestruns(scenario,useexp);
keepids=setdiff(1:infos(scenario,currun).stateInfo.N, []);
[F N]=size(infos(scenario,currun).stateInfo.X);
zframes=190:F;
curid=30;
infos(scenario,currun).stateInfo.X(zframes,curid)=0;
infos(scenario,currun).stateInfo.Y(zframes,curid)=0;
infos(scenario,currun).stateInfo.Xgp(zframes,curid)=0;
infos(scenario,currun).stateInfo.Ygp(zframes,curid)=0;
infos(scenario,currun).stateInfo.Xi(zframes,curid)=0;
infos(scenario,currun).stateInfo.Yi(zframes,curid)=0;
infos(scenario,currun).stateInfo.H(zframes,curid)=0;
infos(scenario,currun).stateInfo.W(zframes,curid)=0;

infos(scenario,currun).stateInfo.targetsExist=getTracksLifeSpans(infos(scenario,currun).stateInfo.X);
infos(scenario,currun).stateInfo.N=N;
infos(scenario,currun).stateInfo= ...
    matricesToVector(infos(scenario,currun).stateInfo.X,infos(scenario,currun).stateInfo.Y,infos(scenario,currun).stateInfo);

[metrics2d metrics3d]=printFinalEvaluation(infos(scenario,currun).stateInfo, gtInfo, sceneInfo, opt);

mets2d(scenario,:,currun)=metrics2d;
mets3d(scenario,:,currun)=metrics3d;