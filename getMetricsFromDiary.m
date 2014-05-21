%%
fid = fopen('testingparams.txt');
lncnt=0;
tline = fgets(fid);
allmets=[];
while ischar(tline)     && lncnt<2447
    lncnt=lncnt+1;
    tline = fgets(fid);    
    if ~mod(lncnt+1,6)
%         disp(tline)
        rdline=sscanf(tline,' %f %f %f| %d  %d  %d  %d|  %d  %d  %d  %d|  %f  %f  %f');
        allmets=[allmets; rdline'];
    end
end

fclose(fid);