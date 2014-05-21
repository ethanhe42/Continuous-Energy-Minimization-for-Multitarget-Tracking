function [stateVec N F targetsExist X Y]=getStateInfo(stateInfo)
% unmake stateInfo
% get number of targets, number of frames
% state vector, mapping and matrices
% 


stateVec=stateInfo.stateVec;
N=stateInfo.N;F=stateInfo.F;
targetsExist=stateInfo.targetsExist;

[X Y]=vectorToMatrices(stateVec,stateInfo);