
%% compile mex files

srcFiles={'Edet_mex','Edyn_mex','Eexc_mex','Eper_mex','Eori_mex','CLEAR_MOT_mex','vectorToMatrices_mex', ...
   'computeOcclusions2new_mex','gausshist','gausshistdx'};
srcExt='c';
srcdir=fullfile('mex','src');
outdir=fullfile('mex','bin');

if ~exist(outdir,'dir'), mkdir(outdir); end

opts='-silent'; if ispc, opts=''; end

for k=1:length(srcFiles)
    eval(sprintf('mex %s -outdir %s %s.%s',opts,outdir,fullfile(srcdir,char(srcFiles(k))),srcExt));
end

%% compile utils
utildir=fullfile('external','motutils');
srcFiles={'allWorldToImage_mex','allWorldToImageWithImDer_mex'};
srcExt='c';
srcdir=fullfile(utildir,'mex','src');
outdir=fullfile(utildir,'mex','bin');

if ~exist(outdir,'dir'), mkdir(outdir); end
for k=1:length(srcFiles)
    eval(sprintf('mex %s -outdir %s %s.%s',opts,outdir,fullfile(srcdir,char(srcFiles(k))),srcExt));
end