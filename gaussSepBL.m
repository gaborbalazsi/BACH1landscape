% important code: defines BL invasion landscape for paper
% it is not only Gaussian, as the filename would sugest
% Produces SI modeling panels A and B
% It produces text files utilized by ouBLinLandscUPD.m

clearvars; close all;

doxRange=[    0    0.1    0.3    0.35    0.5    0.6    1    2    10];

%average fluorescence (BL cells)
GFPave=[...
    2.979016133
    3.398626132
    3.648960343
    3.788605752
    3.809999739
    3.823816377
    4.062188809
    4.202669121
    4.322541832];

%GFPave=GFPave/log(10);

GFPcv=[...   
    0.0670
    0.0690
    0.0757
    0.0785 %0.0637 - corrected
    0.0751
    0.0764
    0.0797
    0.0802
    0.0902];

GFPstd=GFPcv.*GFPave;

%invasiveness (Low mNF-BACH1)
bInv=10;
sInv=100;
InvFL=[...
    27.55 	25.66 	24.27 
    30.92 	30.29 	33.33 
    35.93 	31.15 	35.86 
    30.39 	26.63 	25.13 
    13.97 	11.42 	15.74 
    16.25 	17.58 	18.06 
    16.62 	21.42 	22.65 
    23.79 	29.70 	31.43 
    38.35 	36.83 	39.37];

mInv=mean(InvFL')/100;

%define a function for the seeded cells' CV versus the mean
% p1 =    0.008908;
% p2 =    -0.05025;
% p3 =      0.1378;
p1 =    0.007183;
p2 =     -0.0376;
p3 =      0.1152;
quadfit=@(x) p1*x.^2 + p2*x + p3;

%average GFP fluorescence of Invaded BL cells
invGFPave=[...
    3.37356879
    3.533857401
    3.609787441
    3.735868215
    3.785768841
    3.859890899
    4.394765724
    4.352679483
    4.502049469];

%CV of Invaded BL cells
invGFPcv=[...
    0.054
    0.0682
    0.0744
    0.0427
    0.1238
    0.1140
    0.0377
    0.0701
    0.0611];

GFPstd=GFPcv.*GFPave;

%define landscape
xx=linspace(1.5,9,4000);
% LandSc=pchip(GFPave,[bInv mInv sInv],xx);
% RateSc=-log(1-LandSc/100)/24;

%length of input data or function
numA=length(xx);
%xEdges=linspace(2.5,5,50);yEdges=xEdges;
xEdges=linspace(0.5,6,100);yEdges=xEdges;
xCenters=(xEdges(1:end-1)+xEdges(2:end))/2;yCenters=(yEdges(1:end-1)+yEdges(2:end))/2;

% fit parameters for each individual Dox concentration
fitparsBL=load('fitparsBL.txt');
slopSr=fitparsBL(:,1);  %slope
intcSr=fitparsBL(:,2);  %intercept
slopIr=fitparsBL(:,3);  %slope
intcIr=fitparsBL(:,4);  %intercept

%reserve memory for storing landscapes
LandscExp=zeros(length(GFPave),length(xCenters));
LandscEst=zeros(length(GFPave),length(xCenters));
LandscOri=LandscExp;
LandscWgt=LandscExp;

% for each induction level
for currind=1:length(GFPave)
    currind
    doxC=doxRange(currind);
    muNF=GFPave(currind);
    %sigmaNF=GFPcv(currind)*muNF;
    cvEst=quadfit(muNF);    % quadratic estimate of CV from the mean
    sigmaNF=GFPstd(currind);
    
    %generate histograms and statistics from indidividual seeded & invading cell fluorescence data
    sheetname=sprintf('%d',currind);
    
    %seeded
    ExpData=xlsread('./data231/MB231_1.1BLdataGB.xlsx',sheetname);
    ExpSeeded=ExpData(:,1:3);   %get seeded cell fluorescence values
    ExpSeeded1=ExpSeeded(~isnan(ExpSeeded(:,1)),1);
    ExpSeeded2=ExpSeeded(~isnan(ExpSeeded(:,2)),2);
    ExpSeeded3=ExpSeeded(~isnan(ExpSeeded(:,3)),3);
    ExpSeeded=ExpSeeded(~isnan(ExpSeeded(:)));
    [sBLexp1(currind,:),~]=histcounts(ExpSeeded1,xEdges,'Normalization','pdf');
    [sBLexp2(currind,:),~]=histcounts(ExpSeeded2,xEdges,'Normalization','pdf');
    [sBLexp3(currind,:),~]=histcounts(ExpSeeded3,xEdges,'Normalization','pdf');
    [sBLexp(currind,:),~]=histcounts(ExpSeeded,xEdges,'Normalization','pdf');
    sBLexpAVE=mean([sBLexp1(currind,:);sBLexp2(currind,:);sBLexp3(currind,:)]);
    sBLexpSTD=std([sBLexp1(currind,:);sBLexp2(currind,:);sBLexp3(currind,:)]);
    sBLexpSTA=smooth(sBLexpAVE).^slopSr(currind)/exp(-intcSr(currind));
    sBLexpCVR=sBLexpSTD./sBLexpAVE;
    sBLexpCVA=sBLexpSTA./smooth(sBLexpAVE);

    %invading
    ExpInvded=ExpData(:,4:6);
    ExpInvded1=ExpInvded(~isnan(ExpInvded(:,1)),1);
    ExpInvded2=ExpInvded(~isnan(ExpInvded(:,2)),2);
    ExpInvded3=ExpInvded(~isnan(ExpInvded(:,3)),3);
    ExpInvded=ExpInvded(~isnan(ExpInvded(:)));
    [iBLexp1(currind,:),~]=histcounts(ExpInvded1,xEdges,'Normalization','pdf');
    [iBLexp2(currind,:),~]=histcounts(ExpInvded2,xEdges,'Normalization','pdf');
    [iBLexp3(currind,:),~]=histcounts(ExpInvded3,xEdges,'Normalization','pdf');
    [iBLexp(currind,:),~]=histcounts(ExpInvded,xEdges,'Normalization','pdf');
    iBLexp(currind,:)=mInv(currind).*iBLexp(currind,:);
    iBLexp1(currind,:)=mInv(currind).*iBLexp1(currind,:);
    iBLexp2(currind,:)=mInv(currind).*iBLexp2(currind,:);
    iBLexp3(currind,:)=mInv(currind).*iBLexp3(currind,:);
    iBLexpAVE=mean([iBLexp1(currind,:);iBLexp2(currind,:);iBLexp3(currind,:)]);
    iBLexpSTD=std([iBLexp1(currind,:);iBLexp2(currind,:);iBLexp3(currind,:)]);
    iBLexpSTA=smooth(iBLexpAVE).^slopIr(currind)/exp(-intcIr(currind));
    iBLexpCVR=iBLexpSTD./iBLexpAVE;
    iBLexpCVA=iBLexpSTA./smooth(iBLexpAVE);
    
    %normalizing histograms for cutting the tails later
%     [sBLexpC(currind,:),~]=histcounts(mean(ExpSeeded)+(ExpSeeded-mean(ExpSeeded))/std(ExpSeeded),xEdges,'Normalization','pdf');
%     [iBLexpC(currind,:),~]=histcounts(mean(ExpInvded)+(ExpInvded-mean(ExpInvded))/std(ExpInvded),xEdges,'Normalization','pdf');
    
    %defining Gaussians for cutting the tails later
    sBLexpC(currind,:)=pdf('Normal',xCenters,mean(ExpSeeded),1);
    iBLexpC(currind,:)=pdf('Normal',xCenters,mean(ExpInvded),1);
    
    %compare Gaussian historgrams with estimated std to actual histograms 
    figure;hold on;
    stdEst=muNF*quadfit(muNF);xEst=pdf('Normal',xCenters,muNF,stdEst);plot(xCenters,xEst,'k--');
    xExp=pdf('Normal',xCenters,muNF,sigmaNF);plot(xCenters,xExp,'k+');
    yExp=mInv(currind).*pdf('Normal',xCenters,invGFPave(currind),invGFPave(currind).*invGFPcv(currind));plot(xCenters,yExp,'m+');
    plot(xCenters,sBLexp(currind,:),'k');
    plot(xCenters,iBLexp(currind,:),'m');
    errorbar(xCenters,sBLexpAVE,sBLexpSTD,'ko');
    errorbar(xCenters,iBLexpAVE,iBLexpSTD,'mo');
    
    [trapz(xCenters,xExp) trapz(xCenters,yExp) trapz(xCenters,sBLexp(currind,:)) trapz(xCenters,iBLexp(currind,:))]
    
    %normal pdf ratios
%     yxRatio=smooth(yExp,'sgolay')./smooth(xExp,'sgolay');
%     LandscEst(currind,:)=yxRatio;
    %if any(LandscEst(currind,:)>1)
    %LandscEst(currind,:)=LandscEst(currind,:)/max(LandscEst(currind,:)); end    %rescaling between 0 and 1
    %LandscEst(currind,:)=min(LandscEst(currind,:),ones(size(LandscEst(currind,:))));    %cutoff at 1
%     LandscEst(currind,:)=max(LandscEst(currind,:),zeros(size(LandscEst(currind,:))));
%     indDrop= (smooth(yExp)+smooth(xExp))<0.05;
%     LandscEst(currind,indDrop)=NaN;
    
    %invaded vs. seeded ratio
    siBLexp=smooth(xCenters,iBLexp(currind,:));ssBLexp=smooth(xCenters,sBLexp(currind,:));
    siBLexpC=smooth(xCenters,iBLexpC(currind,:));ssBLexpC=smooth(xCenters,sBLexpC(currind,:));
    isRatio=siBLexp./ssBLexp;isRatio=isRatio';
    %isRatio=iBLexpAVE./sBLexpAVE;  %simple ratio
    %LandscExp(currind,:)=isRatio;
    %LandscExp(currind,:)=isRatio.*(1+1./sBLexp');%corrected ratio, correction #1
    LandscExp(currind,:)=isRatio./(1+sBLexpCVA'.^2/3); %corrected ratio, correction #2
    
    %if any(isRatio>1) LandscExp(currind,:)=LandscExp(currind,:)/max(isRatio(isfinite(isRatio))); end
    
    %LandscExp(currind,:)=min(LandscExp(currind,:),ones(size(isRatio)));
    LandscExp(currind,:)=max(LandscExp(currind,:),zeros(size(isRatio)));
%     LandscExp(currind,:)=min(isRatio,ones(size(isRatio)));
%     LandscExp(currind,:)=max(isRatio,zeros(size(isRatio)));   %this line is problematic. It is undoing the previous line.
    

    %indDrop= (siBLexp+ssBLexp)<0.5;
    %indDrop= (siBLexp/max(siBLexp)+ssBLexp/max(ssBLexp))<1.0;
    %indDrop= (ssBLexp/max(ssBLexp))<0.05;
    %indDrop= (siBLexp/max(siBLexp))<0.18;
    
    %WeightFun=1./sqrt(mInv(currind)./siBLexp + 1./ssBLexp); %1/CV
    %WeightFun=(1+1./ssBLexp)./sqrt( mInv(currind)./siBLexp + 1./ssBLexp + 3*mInv(currind)./siBLexp./ssBLexp + 8*1./ssBLexp./ssBLexp ); %1/CV corrected
    %WeightFun=1./sqrt(iBLexpCVR.^2 + sBLexpCVR.^2); %1/CV
    %WeightFun=(1+sBLexpCVR.^2)./sqrt( iBLexpCVR.^2 + sBLexpCVR.^2 + 3*iBLexpCVR.^2.*sBLexpCVR.^2 + 8*sBLexpCVR.^4 ); %1/CV corrected
    WeightFun=(1+sBLexpCVA.^2)./sqrt( iBLexpCVA.^2 + sBLexpCVA.^2 + 3*iBLexpCVA.^2.*sBLexpCVA.^2 + 8*sBLexpCVA.^4 ); %1/CV corrected
    %WeightFun= siBLexp/mInv(currind);
    %WeightFun= siBLexp.*ssBLexp/mInv(currind);
    
    %cut off unreliable tails from histograms
    %CutFun= 1./abs(WeightFun);indDrop=find(CutFun>0.3);
    CutFun= siBLexpC.*ssBLexpC;indDrop=find(CutFun<0.145);
    %CutFun= siBLexpC;indDrop=find(CutFun<0.38);
    %CutFun= siBLexp.*ssBLexp/max(ssBLexp)/max(siBLexp)<0.01;
    
    LandscOri(currind,:)=LandscExp(currind,:);
    LandscExp(currind,indDrop)=NaN;
    WeightMat(currind,:)=WeightFun;
    WeightMat(currind,indDrop)=NaN;
    %LandscWgt(currind,:)=LandscExp(currind,:).*WeightMat(currind,:);
    
    legend('NFctr','NFest','NFinv','Location','NorthWest')
    xlabel ( 'BACH1 level (arb. units)', 'FontSize', 24)
    ylabel ( 'Probability', 'FontSize', 24);
    tstr=sprintf("Dox=%.2f",doxC);
    title ( tstr, 'FontSize', 24 )
    grid on; hold on;
    set(gca,'XLim',[2 6], 'FontSize', 24);
    %set(gca,'FontSize', 24);
    
    figure;hold on;
    plot(log(sBLexpAVE),log(sBLexpSTD),'ko');
    plot(log(iBLexpAVE),log(iBLexpSTD),'ro');
    plot(log(sBLexpAVE),log(sBLexpSTA),'k+');
    plot(log(iBLexpAVE),log(iBLexpSTA),'r+');
    legend('NFseed','NFinv','Location','southeast')
    xlabel ( 'log(<barheight>) (arb. units)', 'FontSize', 24)
    ylabel ( 'log(bar std)', 'FontSize', 24);
%     tstr=sprintf("Dox=%.2f",doxC);
    title ( tstr, 'FontSize', 24 )
    grid on;
%     set(gca,'XLim',[2 6], 'FontSize', 24);

    %linear fit to log(BLexpSTD) versus log(BLexpAVE) when neither is zero, ensured by indices (BLexpSTD) & (BLexpAVE)
    [cs,qs]=polyfit(log(sBLexpAVE((sBLexpSTD) & (sBLexpAVE))),log(sBLexpSTD((sBLexpSTD) & (sBLexpAVE))),1);   %linear fit to seeded
    [ci,qi]=polyfit(log(iBLexpAVE((iBLexpSTD) & (iBLexpAVE))),log(iBLexpSTD((iBLexpSTD) & (iBLexpAVE))),1);   %linear fit to invded
    slopS(currind)=cs(1);intcS(currind)=cs(2);  %obtain fit parameters for std versus mean of seeded cells
    slopI(currind)=ci(1);intcI(currind)=ci(2);  %obtain fit parameters for std versus mean of invading cells
    plot(-10:0.1:2,slopS(currind)*(-10:0.1:2)+intcS(currind),'k');
    plot(-10:0.1:2,slopI(currind)*(-10:0.1:2)+intcI(currind),'r');
    plot(-10:0.1:2,slopSr(currind)*(-10:0.1:2)+intcSr(currind),'k--');
    plot(-10:0.1:2,slopIr(currind)*(-10:0.1:2)+intcIr(currind),'r--');
end;
% LandscExp(5,:)=LandscExp(5,:)-min(LandscExp(5,:));
% LandscExp(6,:)=LandscExp(6,:)-min(LandscExp(6,:));

% LandscExp(2,:)=circshift(LandscExp(2,:),-7);
% WeightMat(2,:)=circshift(WeightMat(2,:),-7);

LandscWgt=LandscExp.*WeightMat;

% LandscWgt(1,1:9)=0;
% LandscWgt(2,:)=circshift(LandscWgt(2,:),-7);
% WeightMat(2,:)=circshift(WeightMat(2,:),-7);

% get correlations between successive landscape estimates, to check consistency
sMin=-10;sMax=25;sRange=sMin:sMax;
LowrCorr=NaN*zeros(size(LandscExp,1)-2,sMax-sMin+1);
UpprCorr=LowrCorr;
for r=2:size(LandscExp,1)-1
    for s=sMin:sMax
%         LowrTmp=corrcoef(LandscExp(r-1,:),circshift(LandscExp(r,:),-s),'rows','complete');
%         UpprTmp=corrcoef(LandscExp(r,:),circshift(LandscExp(r+1,:),-s),'rows','complete');
        
        LowrTmp=corrcoef(LandscWgt(r-1,:),circshift(LandscWgt(r,:),-s),'rows','complete');
        UpprTmp=corrcoef(LandscWgt(r,:),circshift(LandscWgt(r+1,:),-s),'rows','complete');
        LowrCorr(r-1,s-sMin+1)=LowrTmp(1,2);
        UpprCorr(r-1,s-sMin+1)=UpprTmp(1,2);
    end
end
figure;imagesc(sRange,2:size(LandscExp,1)-1,LowrCorr);colorbar;title('LowrCorr');
figure;imagesc(sRange,2:size(LandscExp,1)-1,UpprCorr);colorbar;title('UpprCorr');

%plot individual landscape estimates
figure;plot(xCenters,LandscEst);
legCell=cellstr(num2str(doxRange', 'Dox=%.2f'));legend(legCell,'Location','EastOutside')
title('Gaussian estimate landscape');

%plot individual landscape estimates
figure;plot(xCenters,LandscExp,'LineWidth',2);
legCell=cellstr(num2str(doxRange', 'Dox=%.2f'));legend(legCell,'Location','EastOutside')
title('Local inv. landscape');
set(gca,'FontSize',24,'YLim',[0 1])

%plot weighted landscapes
figure;plot(xCenters,LandscWgt);
legCell=cellstr(num2str(doxRange', 'Dox=%.2f'));legend(legCell,'Location','EastOutside')
title('Weighted landscape');

%plot histogram product used for cutoff
figure;plot(xCenters,sBLexpC.*iBLexpC);
legCell=cellstr(num2str(doxRange', 'Dox=%.2f'));legend(legCell,'Location','EastOutside')
title('Cutting histogram');

% plot experimental landscape, estimated as corrected invaded vs. seeded histogram ratio
figure;h=imagesc(xCenters,GFPave,LandscExp);colorbar;set(h, 'AlphaData', ~isnan(LandscExp))
text(2.5,2.85,'Dox=');
for ctr=1:length(doxRange)
    dox=doxRange(ctr);
    doxstr=sprintf('%.2f',dox);
    text(2.55,min(GFPave)+(ctr-1)*((max(GFPave)-min(GFPave))/(length(GFPave)-1)),doxstr);
end
title('Local inv. landscape');

% plot experimental landscape, estimated as weighted invaded vs. seeded histogram ratio
figure;h=imagesc(xCenters,GFPave,LandscWgt);colorbar;set(h, 'AlphaData', ~isnan(LandscWgt))
text(2.5,2.85,'Dox=');
for ctr=1:length(doxRange)
    dox=doxRange(ctr);
    doxstr=sprintf('%.2f',dox);
    text(2.55,min(GFPave)+(ctr-1)*((max(GFPave)-min(GFPave))/(length(GFPave)-1)),doxstr);
end
title('LandscWgt');

%aveLandscExp=nanmean(LandscExp);           %arithmetic mean
%aveLandscExp=exp(nanmean(log(LandscExp)));  %geometric mean
aveLandscExp=nansum(LandscWgt)./nansum(WeightMat);           %arithmetic mean
aveLandscExp(1:find(~isnan(aveLandscExp), 1 )-1)=0;
fillmissing(aveLandscExp,'linear');
aveLandscExp(find(~isnan(aveLandscExp), 1, 'last' )+1:end)=1;   % here you set the high end
%aveLandscExp(1:48)=aveLandscExp(1:48)+xCenters(1:48)*(aveLandscExp(49)-0.18)/xCenters(49)+0.18; % a possible linear correction on low end

%aveLandscExp=exp(nanmean(log(LandscExp)));  %geometric mean
figure;plot(xCenters,aveLandscExp)
xlabel('log(BACH1 expression)');ylabel('Invasiveness');

% save average consolidated landscape, x values, and individual landscapes
save -ASCII -TABS 'aveLandscBL.txt' aveLandscExp;
save -ASCII -TABS 'xCentersBL.txt' xCenters;
save -ASCII -TABS 'LandscIndBL.txt' LandscOri;

%smLandscExp=smooth(aveLandscExp,'sgolay');
smLandscExp=smooth(aveLandscExp);
figure;plot(xCenters,smLandscExp)
xlabel('log(BACH1 expression)');ylabel('Invasiveness');

% save smoothed average consolidated landscape
save -ASCII -TABS 'smLandscBL.txt' smLandscExp;


%back-estimating invaded histogram (bHSTest) and average invasiveness (bINVest)
for currind=1:length(GFPave)
    currind
    doxC=doxRange(currind);
    muNF=GFPave(currind);
    %sigmaNF=GFPcv(currind)*muNF;
    cvEst=quadfit(muNF);
    sigmaNF=GFPstd(currind);
    
    figure;hold on;
    siBLexp=smooth(iBLexp(currind,:));ssBLexp=smooth(sBLexp(currind,:));    %smoothed hisograms
    bHSTest(currind,:)=ssBLexp.*LandscOri(currind,:)';  %landscape estimation based on individual landscapes
    %bHSTest(currind,:)=ssBLexp.*aveLandscExp'; %landscape estimation based on average consolidated landscape
    bINVest(currind)=trapz(xCenters(~isnan(bHSTest(currind,:))),bHSTest(currind,~isnan(bHSTest(currind,:))));
    plot(xCenters,ssBLexp,'k','LineWidth',2);
    plot(xCenters,siBLexp,'m','LineWidth',2);
    plot(xCenters,bHSTest(currind,:),'mo');
    title(sprintf('back est. Dox=%.2f',doxC));
    legend('expS','expI','bckI');
    %legend('Seeded','Inv.','Location','NorthWest');
    set(gca,'FontSize',24);
end

%plot xperimental versus back-estimated invasion landscapes
figure;hold on;
plot(GFPave,mInv,'m+--');
plot(GFPave,bINVest,'ro')
legend('expI','bckI','Location','NorthWest');
title('measured versus back-estimated invasiveness')

% this is the statistics of log(bar heights) in flow cytometry histograms
% after plotting each bar's log(standard deviation) vesrus the bar's log(mean height)
% we observe linear dependence and perform linear fits
% we plot the fit parameters for seeded and invading cells

% plot slopes of bars' log(std) versus log(mean) values for seeded and invaded
figure;hold on;
plot(slopS,'ko');hold on;plot(slopI,'ro');
legend('seeded','invading');
xlabel('induction level');ylabel('slopes, log(mean height) vs. log(std height)')

% plot intercepts of bars' log(std) versus log(mean) values for seeded and invaded
figure;hold on;
plot(intcS,'ko');hold on;plot(intcI,'ro');
legend('seeded','invading');
xlabel('induction level');ylabel('intcepts, log(mean height) vs. log(std height)')

[mean(slopS) std(slopS) mean(intcS) std(intcS)]
[mean(slopI) std(slopI) mean(intcI) std(intcI)]

% fitparsBL=[slopS' intcS' slopI' intcI']
% save -ASCII -TABS 'fitparsBL.txt' fitparsBL
