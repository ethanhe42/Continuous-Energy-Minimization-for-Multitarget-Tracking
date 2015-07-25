%% generate parameter files for grid search
% function generateOptGrid
opt=getConOptions;
searchpar.wtEdyn=1:2;
searchpar.wtEexc=.5:.5:1.5;
searchpar.wtEper=.5:.5:1;
searchpar.wtEreg=.25:.25:.75;
searchpar.lambda=0.05:0.05:0.15;

outdir='config/KITTI';


%%%%%%%
opt.track3d=1; opt.cutToTA=1;
searchpar.wtEdyn=[.01, .03, .05];
searchpar.wtEexc=.5:.5:1.5;
searchpar.wtEper=.5:.5:1;
searchpar.wtEreg=.25:.25:.75;
searchpar.lambda=0.05:0.05:0.15;
outdir='config/PETSTUD/inis';


fnames=fieldnames(searchpar);
fn=length(fnames);
parHCube=zeros(1,fn);

% get dimensions
for f=1:fn
    parHCube(f) = length(getfield(searchpar,char(fnames(f))));
end

npar=prod(parHCube);
fprintf('We will generate %d parameter files...\n',npar);


HCubeSize=parHCube;
parHCube=zeros(HCubeSize);

ini=IniConfig();


% now fill par cube
for p=1:npar
    outpar=sprintf('p%d,',1:fn); outpar=outpar(1:end-1);
    outpar=['[' outpar ']'];
    evstr=[outpar '=ind2sub(HCubeSize,p);']
    eval(evstr);
    eval(sprintf('ind=%s;',outpar));
    for f=1:fn
        fname=char(fnames(f));
        pstr=sprintf('opt.%s=searchpar.%s(%d);',fname,fname,ind(f))
        eval(pstr);
%         pause
%         eval(sprintf('p%d',f));
    end
    wstatus = writeConOptions(opt,sprintf('%s/%04d.ini',outdir,p));
%     eval()
%     outpar
%     [a,b,c,d,e]=ind2sub(HCubeSize,p)
end
