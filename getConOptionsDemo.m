function opt=getConOptionsDemo(w)
% 



% general
opt.track3d=1;                  % set to 1 for track estimation on ground plane
opt.verbosity=2;                % 0=silent, 1=short info, 2=long info, 3=all
opt.mex=1;                      % use mex
opt.visOptim=1;                 % visualize optimization
opt.occ=0;                      % compute occlusions [Andriyenko et al. ICCV VS Workshop 2011]
                                % only works for 3d tracking for now!
opt.cutToTA=1;                  % cut detections, ground truth and result to tracking area


% optimization
opt.jumpsOrder=[1 3 4 2 6 5];   % standard: merge grow shrink split add remove
opt.maxEpochs=15;                % max global iterations (rounds)
opt.maxIterCGD=30;              % max iterations for each gradient descent

% energy weights (default 2d)
opt.wtEdet=1;               % should be kept at 1
opt.wtEdyn=1;
opt.wtEexc=1;
opt.wtEper=.5;
opt.wtEreg=1;
opt.wtEapp=0;

% energy weights (default 3d)
if opt.track3d
    opt.wtEdet=1;               % should be kept at 1
    opt.wtEdyn=.02;
    opt.wtEexc=.5;
    opt.wtEper=.5;
    opt.wtEreg=1;

end

% other parameters
opt.lambda=0.125;

if nargin==1
    opt.wtEdet=w(1);               % should be kept at 1
    opt.wtEdyn=w(2);
    opt.wtEexc=w(3);
    opt.wtEper=w(4);
    opt.wtEreg=w(5);
    opt.wtEreg=w(6);
    
end


% opt.wtEdet=1;
% opt.wtEdyn=1;
% opt.wtEexc=1;
% opt.wtEper=1;
% opt.wtEreg=1;




end