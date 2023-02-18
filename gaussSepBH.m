% important code: defines BH invasion landscape for paper
% it is not only Gaussian, as the filename would sugest
% Produces SI figure comparing BL and BH invasion landscapes
% It produces text files utilized by ouBHinLandscUPD.m

clearvars; close all;

%doxRange=[    0    0.1    0.3    0.35    0.5    0.6    1    2    10];
doxRange=[    0    0.1    0.3    0.35    0.5    0.6    1    10];

%average fluorescence (BH cells)
GFPave=[...
    2.993826666
    3.488107435
    3.725761448
    3.75081275
    3.759890103
    3.828228019
    3.967407157
%    4.287067948
    4.196857051];

%GFPave=GFPave/log(10);

GFPcv=[...   
    0.0776
    0.0848
    0.0855
    0.0888
    0.0838
    0.0815
    0.1082
%    0.0867
    0.1100];

GFPstd=GFPcv.*GFPave;

%invasiveness (Low mNF-BACH1)
bInv=10;
sInv=100;
InvFL=[...
    27.80 	23.33 	29.09 
    28.17 	31.08 	24.74 
    25.29 	30.24 	30.41 
    24.31 	23.17 	19.35 
    18.30 	21.59 	18.10 
    21.14 	22.34 	20.52 
    26.48 	20.08 	25.00 
%    39.59 	29.32 	29.22 
    43.35 	55.02 	46.87 ];

mInv=mean(InvFL')/100;

%define a function for the seeded cells' CV versus the mean
% p1 =    0.008908;
% p2 =    -0.05025;
% p3 =      0.1378;
p1 =    0.007183;
p2 =     -0.0376;
p3 =      0.1152;
quadfit=@(x) p1*x.^2 + p2*x + p3;

%average GFP fluorescence of Invaded BH cells
invGFPave=[...
    3.313116624
    3.55380731
    3.678800632
    3.660852415
    3.749328014
    3.806324972
    4.306124068
%    4.428730165
    4.525960097];

%CV of Invaded BH cells
invGFPcv=[...
    0.0699
    0.0862
    0.0855
    0.0847
    0.1264
    0.1257
    0.0439
%    0.0733
    0.0796];

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

fitparsBL=load('fitparsBH.txt');
slopSr=fitparsBL(:,1);
intcSr=fitparsBL(:,2);
slopIr=fitparsBL(:,3);
intcIr=fitparsBL(:,4);

LandscExp=zeros(length(GFPave),length(xCenters));
LandscEst=zeros(length(GFPave),length(xCenters));
LandscOri=LandscExp;
LandscWgt=LandscExp;
for currind=1:length(GFPave)
    currind
    doxC=doxRange(currind);
    muNF=GFPave(currind);
    %sigmaNF=GFPcv(currind)*muNF;
    cvEst=quadfit(muNF);
    sigmaNF=GFPstd(currind);
    
    sheetname=sprintf('%d',currind);
    ExpData=xlsread('./data231/MB231_1.1BHdataGB.xlsx',sheetname);
    ExpSeeded=ExpData(:,1:3);
    ExpSeeded1=ExpSeeded(~isnan(ExpSeeded(:,1)),1);
    ExpSeeded2=ExpSeeded(~isnan(ExpSeeded(:,2)),2);
    ExpSeeded3=ExpSeeded(~isnan(ExpSeeded(:,3)),3);
    ExpSeeded=ExpSeeded(~isnan(ExpSeeded(:)));
    [sBHexp(currind,:),~]=histcounts(ExpSeeded,xEdges,'Normalization','pdf');
    [sBHexp1(currind,:),~]=histcounts(ExpSeeded1,xEdges,'Normalization','pdf');
    [sBHexp2(currind,:),~]=histcounts(ExpSeeded2,xEdges,'Normalization','pdf');
    [sBHexp3(currind,:),~]=histcounts(ExpSeeded3,xEdges,'Normalization','pdf');
    ExpInvded=ExpData(:,4:6);
    ExpInvded1=ExpInvded(~isnan(ExpInvded(:,1)),1);
    ExpInvded2=ExpInvded(~isnan(ExpInvded(:,2)),2);
    ExpInvded3=ExpInvded(~isnan(ExpInvded(:,3)),3);
    ExpInvded=ExpInvded(~isnan(ExpInvded(:)));
    [iBHexp(currind,:),~]=histcounts(ExpInvded,xEdges,'Normalization','pdf');
    [iBHexp1(currind,:),~]=histcounts(ExpInvded1,xEdges,'Normalization','pdf');
    [iBHexp2(currind,:),~]=histcounts(ExpInvded2,xEdges,'Normalization','pdf');
    [iBHexp3(currind,:),~]=histcounts(ExpInvded3,xEdges,'Normalization','pdf');
    iBHexp(currind,:)=mInv(currind).*iBHexp(currind,:);
    iBHexp1(currind,:)=mInv(currind).*iBHexp1(currind,:);
    iBHexp2(currind,:)=mInv(currind).*iBHexp2(currind,:);
    iBHexp3(currind,:)=mInv(currind).*iBHexp3(currind,:);
    sBHexpAVE=mean([sBHexp1(currind,:);sBHexp2(currind,:);sBHexp3(currind,:)]);
    sBHexpSTD=std([sBHexp1(currind,:);sBHexp2(currind,:);sBHexp3(currind,:)]);
    sBHexpSTA=smooth(sBHexpAVE).^slopSr(currind)/exp(-intcSr(currind));
    sBHexpCVR=sBHexpSTD./sBHexpAVE;
    sBHexpCVA=sBHexpSTA./smooth(sBHexpAVE);
    iBHexpAVE=mean([iBHexp1(currind,:);iBHexp2(currind,:);iBHexp3(currind,:)]);
    iBHexpSTD=std([iBHexp1(currind,:);iBHexp2(currind,:);iBHexp3(currind,:)]);
    iBHexpSTA=smooth(iBHexpAVE).^slopIr(currind)/exp(-intcIr(currind));
    iBHexpCVR=iBHexpSTD./iBHexpAVE;
    iBHexpCVA=iBHexpSTA./smooth(iBHexpAVE);
    
     %normalizing histograms for cutting the tails later
%     [sBHexpC(currind,:),~]=histcounts(mean(ExpSeeded)+(ExpSeeded-mean(ExpSeeded))/std(ExpSeeded),xEdges,'Normalization','pdf');
%     [iBHexpC(currind,:),~]=histcounts(mean(ExpInvded)+(ExpInvded-mean(ExpInvded))/std(ExpInvded),xEdges,'Normalization','pdf');
    
    %defining Gaussians for cutting the tails later
    sBHexpC(currind,:)=pdf('Normal',xCenters,mean(ExpSeeded),1);
    iBHexpC(currind,:)=pdf('Normal',xCenters,mean(ExpInvded),1);
    
    figure;hold on;
    stdEst=muNF*quadfit(muNF);xEst=pdf('Normal',xCenters,muNF,stdEst);plot(xCenters,xEst,'k--');
    xExp=pdf('Normal',xCenters,muNF,sigmaNF);plot(xCenters,xExp,'k+');
    yExp=mInv(currind).*pdf('Normal',xCenters,invGFPave(currind),invGFPave(currind).*invGFPcv(currind));plot(xCenters,yExp,'m+');
    plot(xCenters,sBHexp(currind,:),'k');
    plot(xCenters,iBHexp(currind,:),'m');
    errorbar(xCenters,sBHexpAVE,sBHexpSTD,'ko');
    errorbar(xCenters,iBHexpAVE,iBHexpSTD,'mo');
    
    [trapz(xCenters,xExp) trapz(xCenters,yExp) trapz(xCenters,sBHexp(currind,:)) trapz(xCenters,iBHexp(currind,:))]
    
    %normal pdf ratios
%     yxRatio=smooth(yExp,'sgolay')./smooth(xExp,'sgolay');
%     LandscEst(currind,:)=yxRatio;
    %if any(LandscEst(currind,:)>1)
    %LandscEst(currind,:)=LandscEst(currind,:)/max(LandscEst(currind,:)); end    %rescaling between ö and 1
    %LandscEst(currind,:)=min(LandscEst(currind,:),ones(size(LandscEst(currind,:))));    %cutoff at 1
%     LandscEst(currind,:)=max(LandscEst(currind,:),zeros(size(LandscEst(currind,:))));
%     indDrop= (smooth(yExp)+smooth(xExp))<0.05;
%     LandscEst(currind,indDrop)=NaN;
    
    %invaded vs. seeded ratio
    siBHexp=smooth(xCenters,iBHexp(currind,:));ssBHexp=smooth(xCenters,sBHexp(currind,:));
    siBHexpC=smooth(xCenters,iBHexpC(currind,:));ssBHexpC=smooth(xCenters,sBHexpC(currind,:));
    isRatio=siBHexp./ssBHexp;isRatio=isRatio';
    %isRatio=iBHexpAVE./sBHexpAVE;
    %LandscExp(currind,:)=isRatio;
    %LandscExp(currind,:)=isRatio.*(1+1./sBHexp');
    %LandscExp(currind,:)=isRatio./(1+sBHexpCVR.^2/3);
    LandscExp(currind,:)=isRatio./(1+sBHexpCVA'.^2/3);
    
    %if any(isRatio>1) LandscExp(currind,:)=LandscExp(currind,:)/max(isRatio(isfinite(isRatio))); end
    
    %LandscExp(currind,:)=min(LandscExp(currind,:),ones(size(isRatio)));
    LandscExp(currind,:)=max(LandscExp(currind,:),zeros(size(isRatio)));
%     LandscExp(currind,:)=min(isRatio,ones(size(isRatio)));
%     LandscExp(currind,:)=max(isRatio,zeros(size(isRatio)));   %this line is problematic. It is undoing the previous line.
    

    %indDrop= (siBHexp+ssBHexp)<0.5;
    %indDrop= (siBHexp/max(siBHexp)+ssBHexp/max(ssBHexp))<1.0;
    %indDrop= (ssBHexp/max(ssBHexp))<0.05;
    %indDrop= (siBHexp/max(siBHexp))<0.18;
    
    %WeightFun=1./sqrt(mInv(currind)./siBHexp + 1./ssBHexp); %1/CV
    %WeightFun=(1+1./ssBHexp)./sqrt( mInv(currind)./siBHexp + 1./ssBHexp + 3*mInv(currind)./siBHexp./ssBHexp + 8*1./ssBHexp./ssBHexp ); %1/CV corrected
    %WeightFun=1./sqrt(iBHexpCVR.^2 + sBHexpCVR.^2); %1/CV
    %WeightFun=(1+sBHexpCVR.^2)./sqrt( iBHexpCVR.^2 + sBHexpCVR.^2 + 3*iBHexpCVR.^2.*sBHexpCVR.^2 + 8*sBHexpCVR.^4 ); %1/CV corrected
    WeightFun=(1+sBHexpCVA.^2)./sqrt( iBHexpCVA.^2 + sBHexpCVA.^2 + 3*iBHexpCVA.^2.*sBHexpCVA.^2 + 8*sBHexpCVA.^4 ); %1/CV corrected
    %WeightFun= siBHexp/mInv(currind);
    %WeightFun= siBHexp.*ssBHexp/mInv(currind);
    
    %CutFun= 1./abs(WeightFun);indDrop=find(CutFun>0.3);
    CutFun= siBHexpC.*ssBHexpC;indDrop=find(CutFun<0.145);
    %CutFun= siBHexpC;indDrop=find(CutFun<0.38);
    %CutFun= siBHexp.*ssBHexp/max(ssBHexp)/max(siBHexp)<0.01;
    
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
    plot(log(sBHexpAVE),log(sBHexpSTD),'ko');
    plot(log(iBHexpAVE),log(iBHexpSTD),'ro');
    plot(log(sBHexpAVE),log(sBHexpSTA),'k+');
    plot(log(iBHexpAVE),log(iBHexpSTA),'r+');
    legend('NFseed','NFinv','Location','southeast')
    xlabel ( 'log(<barheight>) (arb. units)', 'FontSize', 24)
    ylabel ( 'log(bar std)', 'FontSize', 24);
%     tstr=sprintf("Dox=%.2f",doxC);
    title ( tstr, 'FontSize', 24 )
    grid on;
%     set(gca,'XLim',[2 6], 'FontSize', 24);
    [cs,qs]=polyfit(log(sBHexpAVE((sBHexpSTD) & (sBHexpAVE))),log(sBHexpSTD((sBHexpSTD) & (sBHexpAVE))),1);   %linear fit to seeded
    [ci,qi]=polyfit(log(iBHexpAVE((iBHexpSTD) & (iBHexpAVE))),log(iBHexpSTD((iBHexpSTD) & (iBHexpAVE))),1);   %linear fit to invded
    slopS(currind)=cs(1);intcS(currind)=cs(2);
    slopI(currind)=ci(1);intcI(currind)=ci(2);
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

figure;plot(xCenters,LandscEst);
legCell=cellstr(num2str(doxRange', 'Dox=%.2f'));legend(legCell,'Location','EastOutside')
title('Gaussian estimate landscape');

figure;plot(xCenters,LandscExp);
legCell=cellstr(num2str(doxRange', 'Dox=%.2f'));legend(legCell,'Location','EastOutside')
title('Experimental landscape');

figure;plot(xCenters,LandscWgt);
legCell=cellstr(num2str(doxRange', 'Dox=%.2f'));legend(legCell,'Location','EastOutside')
title('Weighted landscape');

figure;plot(xCenters,sBHexpC.*iBHexpC);
legCell=cellstr(num2str(doxRange', 'Dox=%.2f'));legend(legCell,'Location','EastOutside')
title('Cutting histogram');

figure;h=imagesc(xCenters,GFPave,LandscExp);colorbar;set(h, 'AlphaData', ~isnan(LandscExp))
text(2.5,2.85,'Dox=');
for ctr=1:length(doxRange)
    dox=doxRange(ctr);
    doxstr=sprintf('%.2f',dox);
    text(2.55,min(GFPave)+(ctr-1)*((max(GFPave)-min(GFPave))/(length(GFPave)-1)),doxstr);
end
title('LandscExp');

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
aveLandscExp(find(~isnan(aveLandscExp), 1, 'last' )+1:end)=1;
%aveLandscExp=exp(nanmean(log(LandscExp)));  %geometric mean
figure;plot(xCenters,aveLandscExp)
xlabel('log(BACH1 expression)');ylabel('Invasiveness');
save -ASCII -TABS 'aveLandscBH.txt' aveLandscExp;
save -ASCII -TABS 'xCentersBH.txt' xCenters;
save -ASCII -TABS 'LandscIndBH.txt' LandscOri;

figure;hold on;
%smLandscExp=smooth(aveLandscExp,'sgolay');
BLcons=load('smLandscBL.txt');
xCentBL=load('xCentersBL.txt');
plot(xCentBL,BLcons,'LineWidth',2);
smLandscExp=smooth(aveLandscExp);
plot(xCenters,smLandscExp,'Color',[0.1 0.5 1],'LineWidth',2)
legend('BL','BH','Location','NorthWest');
xlabel('log(BACH1 expression)', 'FontSize', 24);ylabel('Invasiveness', 'FontSize', 24);
save -ASCII -TABS 'smLandscBH.txt' smLandscExp;
set(gca,'FontSize', 24,'XLim',[2.5 5]);
title('Cons. inv. landscapes')

%back-estimating invaded histogram (bHSTest) and average invasiveness (bINVest)
for currind=1:length(GFPave)
    currind
    doxC=doxRange(currind);
    muNF=GFPave(currind);
    %sigmaNF=GFPcv(currind)*muNF;
    cvEst=quadfit(muNF);
    sigmaNF=GFPstd(currind);
    
    figure;hold on;
    siBHexp=smooth(iBHexp(currind,:));ssBHexp=smooth(sBHexp(currind,:));
    bHSTest(currind,:)=ssBHexp.*LandscOri(currind,:)';
    %bHSTest(currind,:)=ssBHexp.*aveLandscExp';
    bINVest(currind)=trapz(xCenters(~isnan(bHSTest(currind,:))),bHSTest(currind,~isnan(bHSTest(currind,:))));
    plot(xCenters,ssBHexp,'k');
    plot(xCenters,siBHexp,'m');
    plot(xCenters,bHSTest(currind,:),'mo');
    title(sprintf('%.2f',doxC));
    legend('expS','expI','bckI');
end;

figure;hold on;
plot(GFPave,mInv,'m+--');
plot(GFPave,bINVest,'ro')
legend('expI','bckI','Location','NorthWest');

figure;hold on;
plot(slopS,'ko');hold on;plot(slopI,'ro');
figure;hold on;
plot(intcS,'ko');hold on;plot(intcI,'ro');
[mean(slopS) std(slopS) mean(intcS) std(intcS)]
[mean(slopI) std(slopI) mean(intcI) std(intcI)]

fitparsBH=[slopS' intcS' slopI' intcI']
save -ASCII -TABS 'fitparsBH.txt' fitparsBH
