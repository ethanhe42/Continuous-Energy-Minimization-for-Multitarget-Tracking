%% standard
syms a b c d e f v1 v2 m1 m2 f1 f2 o1 det1x det1y det2x det2y d1x d1y d2x d2y so dist1 dist2 csig epsilon g1 g2 kappa wt intv distrconst
syms ff1b ff2b f1b f2b

v1=([c; d] - [a; b]);
v2=([e; f] - [c; d]);

m1=sqrt((c-a)^2 + (d-b)^2+epsilon);
m2=sqrt((e-c)^2 + (f-d)^2+epsilon);
o1=v1./m1;
o2=v2./m2;

ff1=1/(1+exp(-m1+so));
ff2=1/(1+exp(-m2+so));

dist1=sqrt((a-det1x)^2 + (b-det1y)^2);
dist2=sqrt((c-det2x)^2 + (d-det2y)^2);
f1=ff1*csig/(dist1^2+csig);
f2=ff2*csig/(dist2^2+csig);

ff1b=1-ff1;
ff2b=1-ff2;
f1b=ff1b * (csig - csig/(dist1^2+csig));
f2b=ff2b * (csig - csig/(dist2^2+csig));

d1=[d1x d1y]; d2=[d2x d2y];

% g1=acos((v1(1)*d1x + v1(2)*d1y) / (m1));
% g2=acos((v2(1)*d2x + v2(2)*d2y) / (m2));

% kappa=2;
% van misses
% wt=0.8;
% vmd1=exp(sum(kappa*xx.* muu));
% vmd2=exp(sum(kappa*xx.* -muu));

% intv1=sum(vmd1)*(theta(2)-theta(1));
% intv2=sum(vmd2)*(theta(2)-theta(1));

% vmd1=vmd1/intv1;vmd2=vmd2/intv2;
% g1 = (wt*exp(kappa*(o1(1)*d1x + o1(2)*d1y)) + (1-wt)*exp(kappa*(-o1(1)*d1x - o1(2)*d1y)));
% g2 = (wt*exp(kappa*(o2(1)*d2x + o2(2)*d2y)) + (1-wt)*exp(kappa*(-o2(1)*d2x - o2(2)*d2y)));
g1=exp(kappa*(o1(1)*d1x + o1(2)*d1y));
g2=exp(kappa*(o2(1)*d2x + o2(2)*d2y));

g1=-log(g1/intv);g2=-log(g2/intv);
% vmdm=(wt*vmd1 + (1-wt)*vmd2);

% Eori=f1 * (-(g1-pi/2).^2 + pi^2/4) + ...
%     f2* (-(g2-pi/2).^2 + pi^2/4);

Eori=f1 * g1 +  f2* g2 ...
    + f1b*distrconst + f2b*distrconst;

da=diff(Eori, a);
db=diff(Eori, b);
dc=diff(Eori, c);
dd=diff(Eori, d);
de=diff(Eori, e);
df=diff(Eori, f);

fprintf('diff1 = %s;\n',char(da));
fprintf('diff2 = %s;\n',char(db));
fprintf('diff3 = %s;\n',char(dc));
fprintf('diff4 = %s;\n',char(dd));
fprintf('diff5 = %s;\n',char(de));
fprintf('diff6 = %s;\n',char(df));

%%
cda=ccode(da);
cdb=ccode(db);
cdc=ccode(dc);
cdd=ccode(dd);
cde=ccode(de);
cdf=ccode(df);

charcda=char(cda); charcda=charcda(8:end);charcdb=char(cdb); charcdb=charcdb(8:end);
charcdc=char(cdc); charcdc=charcdc(8:end);charcdd=char(cdd); charcdd=charcdd(8:end);
charcde=char(cde); charcde=charcde(8:end);charcdf=char(cdf); charcdf=charcdf(8:end);



fprintf('\n');
fprintf('diff1 = %s;\n',charcda);
fprintf('diff2 = %s;\n',charcdb);
fprintf('diff3 = %s;\n',charcdc);
fprintf('diff4 = %s;\n',charcdd);
fprintf('diff5 = %s;\n',charcde);
fprintf('diff6 = %s;\n',charcdf);


%% do some tricks
oldstr= {'(a-c)','(b-d)','(c-e)','(d-f)'};
newstr= {'amc','bmd','cme','dmf'};

addoldstr0={'(a*2.0-c*2.0)','(b*2.0-d*2.0)','(c*2.0-e*2.0)','(d*2.0-f*2.0)'};
addnewstr0={'a2c2','b2d2','c2e2','d2f2'};
for addstr=1:length(addoldstr0)
    oldstr{end+1}= addoldstr0{addstr};
    newstr{end+1}= addnewstr0{addstr};
end

addoldstr1={'pow(a-c,2.0)','pow(b-d,2.0)','pow(c-e,2.0)','pow(d-f,2.0)'};
addnewstr1={'amcs','bmds','cmes','dmfs'};
for addstr=1:length(addoldstr1)
    oldstr{end+1}= addoldstr1{addstr};
    newstr{end+1}= addnewstr1{addstr};
end


addoldstr2={'epsilon+amcs+bmds','epsilon+cmes+dmfs','sqrt(exp1)','sqrt(exp2)','pow(exp1,3.0/2.0)','pow(exp2,3.0/2.0)'};
addnewstr2={'exp1','exp2','exp3','exp4','exp5','exp6'};
for addstr=1:length(addoldstr2)
    oldstr{end+1}= addoldstr2{addstr};
    newstr{end+1}= addnewstr2{addstr};
end

addoldstr3={'1.0/exp3','1.0/exp4','1.0/exp5','1.0/exp6'};
addnewstr3={'exp7','exp8','exp9','exp10'};
for addstr=1:length(addoldstr3)
    oldstr{end+1}= addoldstr3{addstr};
    newstr{end+1}= addnewstr3{addstr};
end
addoldstr4={'pow(a-det1x,2.0)','pow(b-det1y,2.0)','pow(c-det2x,2.0)','pow(d-det2y,2.0)'};
addnewstr4={'exp11','exp12','exp13','exp14'};
for addstr=1:length(addoldstr4)
    oldstr{end+1}= addoldstr4{addstr};
    newstr{end+1}= addnewstr4{addstr};
end

addoldstr5={'exp(so-exp3)','exp(so-exp4)','d1x*amc','d1y*bmd','d2x*cme','d2y*dmf'};
addnewstr5={'exp15','exp16','exp17','exp18','exp19','exp20'};
for addstr=1:length(addoldstr5)
    oldstr{end+1}= addoldstr5{addstr};
    newstr{end+1}= addnewstr5{addstr};
end

addoldstr6={'(exp17*exp7+exp18*exp7)','(exp19*exp8+exp20*exp8)','1.0/pow(csig+exp11+exp12,2.0)','1.0/pow(exp16+1.0,2.0)'};
addnewstr6={'exp21','exp22','exp23','exp24'};
for addstr=1:length(addoldstr6)
    oldstr{end+1}= addoldstr6{addstr};
    newstr{end+1}= addnewstr6{addstr};
end

addoldstr7={'exp(-kappa*exp21)','exp(kappa*exp21)','exp(-kappa*exp22)','exp(kappa*exp22)'};
addnewstr7={'exp25','exp26','exp27','exp28'};
for addstr=1:length(addoldstr7)
    oldstr{end+1}= addoldstr7{addstr};
    newstr{end+1}= addnewstr7{addstr};
end
addoldstr8={'(wt-1.0)','1.0/pow(exp15+1.0,2.0)'};
addnewstr8={'exp29','exp30'};
for addstr=1:length(addoldstr8)
    oldstr{end+1}= addoldstr8{addstr};
    newstr{end+1}= addnewstr8{addstr};
end


for repnr=1:length(oldstr)
    charcda=strrep(charcda,oldstr{repnr},newstr{repnr});
    charcdb=strrep(charcdb,oldstr{repnr},newstr{repnr});
    charcdc=strrep(charcdc,oldstr{repnr},newstr{repnr});
    charcdd=strrep(charcdd,oldstr{repnr},newstr{repnr});
    charcde=strrep(charcde,oldstr{repnr},newstr{repnr});
    charcdf=strrep(charcdf,oldstr{repnr},newstr{repnr});
end
    
fprintf('\n');
fprintf('diff1 = %s\n',charcda);
fprintf('diff2 = %s\n',charcdb);
fprintf('diff3 = %s\n',charcdc);
fprintf('diff4 = %s\n',charcdd);
fprintf('diff5 = %s\n',charcde);
fprintf('diff6 = %s\n',charcdf);

%%
for newexp=1:length(oldstr)
    fprintf('%s = %s;\n',char(newstr{newexp}),char(oldstr{newexp}));
end

%% first frame
fprintf('\n\n');

% Eori=acos(d1 * o1);
% diff(Eori, a)
% diff(Eori, b)

