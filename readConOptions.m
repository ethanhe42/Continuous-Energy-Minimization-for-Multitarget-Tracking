function opt=readConOptions(inifile)
% parse configuration for Continuous Tracking

if ~nargin, inifile=[]; end;

if ~exist(inifile,'file')  || isempty(inifile)
    fprintf('WARNING! Config file %s does not exist! Using default setting...\n',inifile);
    inifile='config/default2d.ini';
end


ini=IniConfig();
ini.ReadFile(inifile);

opt=[];


% Main Parameters
opt=fillInOpt(opt,ini,'Parameters');


% General
opt=fillInOpt(opt,ini,'General');
% take care of frames vector
if isfield(opt,'ff') && isfield(opt,'lf')
    opt.frames=opt.ff:opt.lf;
    opt=rmfield(opt,'ff');opt=rmfield(opt,'lf');
end

% Initialization (EKF, DP, ...)
opt=fillInOpt(opt,ini,'Initialization');

% Detections
opt=fillInOpt(opt,ini,'Detections');
if isfield(opt,'sigA') && isfield(opt,'sigB')
    opt.detScale.sigA=opt.sigA;
    opt.detScale.sigB=opt.sigB;
    opt=rmfield(opt,'sigA');opt=rmfield(opt,'sigB');
end

% Appearance
app=fillInOpt(opt,ini,'Appearance');
opt.app=app;

% Misc
opt=fillInOpt(opt,ini,'Miscellaneous');
end


function opt = fillInOpt(opt, ini, sec)
% loop through all keys in section and
% append to struct

keys = ini.GetKeys(sec);
for k=1:length(keys)
    key=char(keys{k});
    opt = setfield(opt,key,ini.GetValues(sec,key));
end

end