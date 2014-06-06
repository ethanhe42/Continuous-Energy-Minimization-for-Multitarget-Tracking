%%
allmets=[];

for j=1:108
    j
%     load(sprintf('results/KITTI/res_%03d.mat',j));
    resultsfile=sprintf('results/KITTI/res_%03d.mat',j);
    metricsKITTI=evalKITTI(resultsfile);
    allmets=[allmets; metricsKITTI];
%     allmets=[allmets; mets2d];
end