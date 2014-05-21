%%
function [bestruns useexp bestmets allm2d expdone reshmets]=paramtest(filename,forceRandrun)
% filename='testrun3d.mat';
if nargin<1
    filename='testrun.mat';
end
if nargin<2
    forceRandrun=0;
end

if ~isdeployed,    addpath(genpath('../contracking')); end
metricsInfo.names.long = {'Recall','Precision','False Alarm Rate', ...
    'GT Tracks','Mostly Tracked','Partially Tracked','Mostly Lost', ...
    'False Positives', 'False Negatives', 'ID Switches', 'Fragmentations', ...
    'MOTA','MOTP', 'MOTA Log'};

metricsInfo.names.short = {'Rcll','Prcn','FAR', ...
    'GT','MT','PT','ML', ...
    'FP', 'FN', 'IDs', 'FM', ...
    'MOTA','MOTP', 'MOTAL'};

metricsInfo.widths.long = [6 9 16 9 14 17 11 15 15 11 14 5 5 8];
metricsInfo.widths.short = [5 5 5 3 3 3 3 4 4 3 3 5 5 5];

metricsInfo.format.long = {'.1f','.1f','.2f', ...
    'i','i','i','i', ...
    'i','i','i','i', ...
    '.1f','.1f','.1f'};

metricsInfo.format.short=metricsInfo.format.long;
homefolder='~'; if ispc, homefolder='D:'; end
% load('/home/aanton/grishome/visinf/projects/ongoing/dctracking/testrun.mat')

filetoload=[homefolder '/visinf/projects/ongoing/contracking/' filename];
if ~exist(filetoload,'file'),    filetoload=filename; end
load(filetoload,'allmets*','allens','allscen','allopts')
% allscen=[23 25 27 71 72];
% allscen=[40 41 42 23];
expdone=find(allmets2d(allscen(1),4,:,size(allmets2d,4)))';
% expdone=find(allmets3d(allscen(1),4,:,1))';
% expdone=1:3;
if usejava('desktop'),    figure(4); end
if ~forceRandrun
    clf
    grid
end
hold on

% get lowest ens
bestruns=zeros(size(allmets2d,1),max(expdone));
allm2d=zeros(size(allmets2d,1),size(allmets2d,2),size(allmets2d,3));
% allscen=setdiff(allscen,allscen(1));
allscen
expdone
for scen=allscen
    le=ones(1,length(expdone));
    for exper=expdone
        [m le]=min(sum(allens(scen,:,exper,:),2));       
%         [m le]=max(allmets3d(scen,12,exper,:)); % TODO % WARNING CHEATING
%         le=6;
% le
        if forceRandrun
            le=forceRandrun;
        end
        
        bestruns(scen,exper)=le;
        allm2d(scen,:,exper)=allmets2d(scen,:,exper,le);
        if howToTrack(scen)==1
            allm2d(scen,:,exper)=allmets3d(scen,:,exper,le);
        end
        
        
    end
    

end
allmets2d=allm2d;
%     pause

%% fill legend

legtext=cell(0);
cnt=0;
for scen=allscen
    cnt=cnt+1;
    legtext(cnt)=cellstr(getSequenceFromScenario(scen));
end


linetype='-';
linestyleorder={'-','.-','-x','-o','-s','v-','^-'};

% set(gca,'LineStyleOrder', '-*|-x|o')
addlw=0;
if ~forceRandrun, addlw=2; end

for met=[12]
    if met==10, linetype='--';
    else linetype='-';
    end
    reshmets=reshape(allm2d(allscen,met,expdone),length(allscen),length(expdone));    
%     reshmets=reshape(allm2d(allscen,met,expdone),length(allscen),length(expdone));    
    

    if met==12
        plot(expdone,reshmets','linewidth',1+addlw);   
    end
% 	expdone
% 	mean(reshmets,1)
	%sum(sum(diff(reshmets')))
    
    if met==10, meanidsw=mean(reshmets,1);
    elseif met==12, meanmota=mean(reshmets,1);
    end

    reshmets
    mean(reshmets,1)
    plot(expdone,mean(reshmets,1),[linetype 'k'],'linewidth',2+addlw);
    if length(allscen)==7
    plot(expdone,mean(reshmets(1:3,:),1),['--' 'k'],'linewidth',2+addlw-1);
    plot(expdone,mean(reshmets(4:7,:),1),[':' 'k'],'linewidth',2+addlw-1);
    end
    legtext(end+1)=cellstr('mean');
%     legend('Campus','Crossing','Stadtmitte','PETS');
    legend(legtext);
    cnt=0;
    for scen=allscen
        cnt=cnt+1;
        ex=0;
        for exper=expdone
            ex=ex+1;
%             text(exper,reshmets(cnt,ex),sprintf('%g',bestruns(scen,exper)));
        end
    end
    
    [maxv bestmax]=max(mean(reshmets,1));
    [minv bestmin]=min(mean(reshmets,1));
    if numel(intersect(met,[1 2 5 6 12 13 14]))
        bestset=bestmax; bestval=maxv;
    else
        bestset=bestmin; bestval=minv;
    end
    bestset=expdone(bestset);
    
    
    for scen=allscen
        prhead=0; if scen==allscen(1), prhead=1; end
        printMetrics(allm2d(scen,:,bestset),metricsInfo,prhead,[1 2 12 13 5 6 7 8 9 11 10]);fprintf('(%d,%d) \n',scen,bestruns(scen,bestset));
    end
    
    fprintf('best: %d %f\n',bestset,bestval);
    useexp=bestset;
    
    bestmets=mean(allmets2d(allscen,:,bestset),1);
    bestmets(4:11)=round(bestmets(4:11));
    printMetrics(bestmets,metricsInfo,1,[1 2 12 13 5 6 7 8 9 11 10]);
    fprintf('\n');
end

if length(allscen)==7
    bestmets=mean(allmets2d(allscen(1:3),:,bestset),1);bestmets(4:11)=round(bestmets(4:11));
    printMetrics(bestmets,metricsInfo,1,[1 2 12 13 5 6 7 8 9 11 10]);
    fprintf('\n');
    bestmets=mean(allmets2d(allscen(4:7),:,bestset),1);bestmets(4:11)=round(bestmets(4:11));
    printMetrics(bestmets,metricsInfo,1,[1 2 12 13 5 6 7 8 9 11 10]);
    fprintf('\n');
end

if 0
for met=[1]
    if met==10, linetype='--';
    else linetype='-';
    end
    reshmets=reshape(allm2d(allscen,met,expdone),length(allscen),length(expdone));    
%     reshmets=reshape(allm2d(allscen,met,expdone),length(allscen),length(expdone));    
    

    if met==12
        plot(expdone,reshmets','linewidth',1+addlw);   
    end
% 	expdone
% 	mean(reshmets,1)
	%sum(sum(diff(reshmets')))
    
    if met==10, meanidsw=mean(reshmets,1);
    elseif met==12, meanmota=mean(reshmets,1);
    end

    reshmets
    mean(reshmets,1)
    plot(expdone,mean(reshmets,1),[linetype 'c'],'linewidth',2+addlw);
    legtext(end+1)=cellstr('mean rcl');
%     legend('Campus','Crossing','Stadtmitte','PETS');
    legend(legtext);
    cnt=0;
    for scen=allscen
        cnt=cnt+1;
        ex=0;
        for exper=expdone
            ex=ex+1;
%             text(exper,reshmets(cnt,ex),sprintf('%g',bestruns(scen,exper)));
        end
    end
    
    [maxv bestmax]=max(mean(reshmets,1));
    [minv bestmin]=min(mean(reshmets,1));
    if numel(intersect(met,[1 2 5 6 12 13 14]))
        bestset=bestmax; bestval=maxv;
    else
        bestset=bestmin; bestval=minv;
    end
    bestset=expdone(bestset);
    
    for scen=allscen
        prhead=0; if scen==allscen(1), prhead=1; end
        printMetrics(allm2d(scen,:,bestset),metricsInfo,prhead,[1 2 12 13 5 6 7 8 9 11 10]);fprintf('(%d,%d) \n',scen,bestruns(scen,bestset));
    end
    
    fprintf('best: %d %f\n',bestset,bestval);
    useexp=bestset;
    
    bestmets=mean(allmets2d(allscen,:,bestset),1);
    bestmets(4:11)=round(bestmets(4:11));
    printMetrics(bestmets,metricsInfo,1,[1 2 12 13 5 6 7 8 9 11 10]);
    fprintf('\n');
end

for met=[2]
    if met==10, linetype='--';
    else linetype='-';
    end
    reshmets=reshape(allm2d(allscen,met,expdone),length(allscen),length(expdone));    
%     reshmets=reshape(allm2d(allscen,met,expdone),length(allscen),length(expdone));    
    

    if met==12
        plot(expdone,reshmets','linewidth',1+addlw);   
    end
% 	expdone
% 	mean(reshmets,1)
	%sum(sum(diff(reshmets')))
    
    if met==10, meanidsw=mean(reshmets,1);
    elseif met==12, meanmota=mean(reshmets,1);
    end

    reshmets
    mean(reshmets,1)
    plot(expdone,mean(reshmets,1),[linetype 'm'],'linewidth',2+addlw);
    legtext(end+1)=cellstr('mean prc');
%     legend('Campus','Crossing','Stadtmitte','PETS');
    legend(legtext);
    
    [maxv bestmax]=max(mean(reshmets,1));    [minv bestmin]=min(mean(reshmets,1));
    if numel(intersect(met,[1 2 5 6 12 13 14]))
        bestset=bestmax; bestval=maxv;
    else
        bestset=bestmin; bestval=minv;
    end
    bestset=expdone(bestset);
    
    for scen=allscen
        prhead=0; if scen==allscen(1), prhead=1; end
        printMetrics(allm2d(scen,:,bestset),metricsInfo,prhead,[1 2 12 13 5 6 7 8 9 11 10]);fprintf('(%d,%d) \n',scen,bestruns(scen,bestset));
    end
    
    fprintf('best: %d %f\n',bestset,bestval);
    useexp=bestset;
    
    bestmets=mean(allmets2d(allscen,:,bestset),1);
    bestmets(4:11)=round(bestmets(4:11));
    printMetrics(bestmets,metricsInfo,1,[1 2 12 13 5 6 7 8 9 11 10]);
    fprintf('\n');
end
end

if exist('meanidsw','var')
    [~,bestdiff]=max(meanmota-meanidsw)
end


hdl = findobj(gca,'Type','line');
MarkerSet = 'so^d*+hvp><.xv';
for q=1:length(hdl)
    modded=mod(q,length(MarkerSet))+1;
    set(hdl(q), 'marker', MarkerSet(modded)); % Assign markers to each line handle
end

reshmets=reshape(allmets2d(allscen,12,expdone),length(allscen),length(expdone));
% ylim([max(20,min(reshmets(:))) 100]);
xlim([expdone(1)-1 expdone(end)+1]);
