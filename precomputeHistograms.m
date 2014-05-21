function allhist=precomputeHistograms(scenario,opt,sceneInfo)


allhist=[];
if ~opt.wtEapp
    return;
end


filtersigma=opt.app.filtersigma;
filtersize=opt.app.filtersize;
ycb=opt.app.ycb;
nbins=opt.app.nbins;
frameNums=sceneInfo.frameNums;
F=length(frameNums);

fH = fspecial('gaussian', filtersize,filtersigma);
[isX isY isc]=size(imread([sceneInfo.imgFolder sprintf(sceneInfo.imgFileFormat,frameNums(1))]));

fprintf('frameNums: %i:%i:%i\n',frameNums(1),diff(frameNums(1:2)),frameNums(end));
fprintf('Precomputing histograms...');


outfolder=sprintf('input/s%04d',scenario);
if ~exist(outfolder,'dir'), mkdir(outfolder); end

if ycb
    allhistfile=sprintf('%s/allhist-%04d-%i-%i-y-%i-%i-%i.mat',outfolder,scenario,nbins,filtersigma,frameNums(1),frameNums(2)-frameNums(1),frameNums(end));
else
    allhistfile=sprintf('%s/allhist-%04d-%i-%i-r-%i-%i-%i.mat',outfolder,scenario,nbins,filtersigma,frameNums(1),frameNums(2)-frameNums(1),frameNums(end));
end

if exist(allhistfile,'file')
    fprintf('loading...');
    load(allhistfile);
else
    % bin center points
    m=1:2:2*nbins;
    m=m/(2*nbins);
    minval=0; maxval=1;
    if ycb, minval=16/255; maxval=240/255; end % YCbCr
    valran=maxval-minval;
    m=m*valran+minval;
    d=1/nbins;
    allhist=zeros(isX,isY,isc,F,'uint8');
    for t=1:F
        if ~mod(t,25), fprintf('.'); end
        
        im=double(imread([sceneInfo.imgFolder sprintf(sceneInfo.imgFileFormat,frameNums(1))]))/255;
        if ycb, im=rgb2ycbcr(im); end
        im = imfilter(im, fH);
        for chan=1:isc
            hh=zeros(size(im,1),size(im,2));
            for bin=1:nbins
                hi=abs(im(:,:,chan)-m(bin))<=d/2;
                if (numel(hi))
                    hh(hi)=bin;
                end
            end
            allhist(:,:,chan,t)=hh;
        end
    end
    save(allhistfile,'allhist');
    fprintf('results saved to %s\n',allhistfile);
end



fprintf('done!\n');

end
