clearvars; close all;

% reads in text files creaed by gausSepBH.m
%same as ouBLindLandsc, but different x bins.
%USED in the paper, SI modeling panel D, E

%*****************************************************************************80
%
%% OU_EXACT applies the exact method to the Ornstein-Uhlenbeck SDE.
%
%  Original code by John Burkardt was modified.   
%
%  Reference:
%
%    Daniel T. Gillespie,
%    Exact numerical simulation of the Ornstein-Uhlenbeck process 
%    and its integral,
%    Physical review. E, 54(2):2084-2091 · September 1996.
%
%  Parameters:
%
%    Input, real TMAX, NTSTEPS, NCELLS, THRESHOLD, SEED, the value of problem parameters.
%
%
%    Input, real TMAX, the final time.
%
%    Input, integer NTSTEPS, the number of time steps.
%
%    Input, integer NCELLS, the maximum number of cells.
%
%    Input, real THRESHOLD, the threshold for survival.
%
%    Input, integer SEED, a seed for the random number generator.
%
% for example: ou_second_NFPFstats(150, 10000, 1000, 0.24, 'shuffle');

tmax=50;
ntsteps=tmax*10;
ncells=5000;
AmpFac=1;   %amplification factor for landscape
CVfac=1.0;
% thetaNF = 10/(tmax/ntsteps);
thetaNF = 100.0;   %actual memory is: thetaNF*tmax/ntsteps
ProbAmp=0.25;
seed='shuffle';

% quadratic fit estimates experimental CV from experiemntal mean
p1 =    0.007183;
p2 =     -0.0376;
p3 =      0.1152;
quadfit=@(x) p1*x.^2 + p2*x + p3;

facL=20;    %for theory

doxRange=[    0    0.1    0.3    0.35    0.5    0.6    1    10];

%binned ladscape estimates
BHdata=load('LandscIndBH.txt'); %Landscape data for each individual Dox
%BHcons=load('aveLandscProdCut13BH.txt');
BHcons=load('smLandscBH.txt');
xCentExp=load('xCentersBH.txt');
figure;
plot(xCentExp,BHdata);
legCell=cellstr(num2str(doxRange', 'Dox=%.2f'));legend(legCell,'Location','EastOutside')

%inducer concentrations (padded at both ends)
dox=[...
    -1
    0
    0.1
    0.3
    0.35
    0.5
    0.6
    1
%    2
    10
    20];

%average fluorescence of Seeded BH cells (padded at both ends)
GFPave=[...
    2.0
    2.993826666
    3.488107435
    3.725761448
    3.75081275
    3.759890103
    3.828228019
    3.967407157
%    4.287067948
    4.196857051
    8];

%CV of Seeded BH cells (padded at both ends)
GFPcv=[...
    0.05    
    0.0776
    0.0848
    0.0855
    0.0888
    0.0838
    0.0815
    0.1082
%    0.0867
    0.1100
    0.1];

GFPstd=GFPcv.*GFPave/CVfac;

%invasiveness (Low mNF-BACH1)  (padded at both ends with bInv, sInv)
bInv=0;
sInv=150;
InvFL=[...
    bInv    bInv    bInv
    27.80 	23.33 	29.09 
    28.17 	31.08 	24.74 
    25.29 	30.24 	30.41 
    24.31 	23.17 	19.35 
    18.30 	21.59 	18.10 
    21.14 	22.34 	20.52 
    26.48 	20.08 	25.00 
%    39.59 	29.32 	29.22 
    43.35 	55.02 	46.877 
    sInv    sInv    sInv];

avInv=mean(InvFL(2:end-1,:)');   % mean of 3 invasion replicates, unpadded
stInv=std(InvFL(2:end-1,:)');  % std of 3 invasion replicates, unpadded
minInv=min(avInv);  % overall landscape miminmum, unpadded
nInv=AmpFac*(avInv-minInv);  %normalized landscape with minimum reset to = 0
%mInv=nInv;  %measured Invsiveness, unpadded
mInv=avInv;


%average GFP fluorescence of Invaded BH cells (not padded)
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

%CV of Invaded BH cells (not padded)
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

GFPstd=GFPcv.*GFPave/CVfac;


%these will be the simulated values
percInv=zeros(8,10);
aveInv=zeros(810);
stdInv=zeros(8,10);
cvrInv=zeros(8,10);

xx=linspace(0,11,1000);

ctr=1;
for doxC=doxRange;
%current Dox concentration
    indD=find(dox==doxC);   %where on the landscape?
    currind=indD-1;
    gfpC=GFPave(indD);      %what mean expression level?
    sigmaNF=GFPstd(indD);      %what std expression level?
    muNF = gfpC;
    varNF=sigmaNF*sigmaNF;
    varChi2NF=2*varNF^4*(1+2*muNF^2);
    
    ExpLandsc=BHdata(currind,:);    %individual landscape
    %ExpLandsc=BHcons;              %consolidated landscape
   
    xRange=xx(xx >= muNF-sigmaNF/4 & xx < muNF+sigmaNF/4);
    LandSc=pchip(GFPave,[bInv mInv sInv],xRange);
    [c1,q1]=polyfit(xRange,LandSc,1);   %linear fit to current landsc. region
    [c2,q2]=polyfit(xRange,LandSc,2);   %quadratic fit to current region

    a1=c1(1);b1=c1(2);  %slope and intercept
    a2=c2(1);b2=c2(2);  %quad coefficients
    
    Dmean=a1*varNF/facL;    %expected change in mean (from theory)
    Dvar=-(a1*varNF/facL )^2 + a2*varChi2NF/facL;   %expect. change in variance
%    Dvar=a2*varChi2NF/facL;
    
    Nmean(ctr)=muNF+Dmean;  %new mean (from theory) for invaded cells
    Nvar(ctr)=varNF+Dvar;   %new variance (from theory)
    
    %read in actual data
    xEdges=linspace(2.5,5,50);yEdges=xEdges;
    sheetname=sprintf('%d',currind);
    ExpData=xlsread('./data231/MB231_1.1BHdataGB.xlsx',sheetname);
    ExpSeeded=ExpData(:,1:3);ExpSeeded=ExpSeeded(~isnan(ExpSeeded(:)));
    [sBHexp(currind,:),~]=histcounts(ExpSeeded,xEdges,'Normalization','pdf');
    ExpInvded=ExpData(:,4:6);ExpInvded=ExpInvded(~isnan(ExpInvded(:)));
    [iBHexp(currind,:),~]=histcounts(ExpInvded,xEdges,'Normalization','pdf');
    iBHexp(currind,:)=mInv(currind).*iBHexp(currind,:);

    NMAX=ncells*2;   %maximum number of cells
    growthRate=0.0;

    % simulate NF
    for rep=1:10
    %
    %  Initialize the random number generator.
    %  The following way to initialize the random number generator 
    %  may not be available in older versions of MATLAB.
    %
        rng ( seed )
    %
    %  Set the discrete time stepsize.
    %
        dt = tmax / ntsteps;
    %
    %  Compute the Brownian increments.
    %
        dw = randn ( NMAX, ntsteps );
    %
    %  Carry out the exact simulation for NF.
    %
        t = linspace ( 0, tmax, ntsteps + 1 );
        xNF = NaN*zeros ( NMAX, ntsteps + 1 );
        NcurrentNF = zeros ( 1 , ntsteps + 1 );
        NinvNF = zeros ( 1 , ntsteps + 1 );

        xNF(1:ncells,1) = muNF + sigmaNF * randn(ncells,1); %prepopulate gene expression 
        yNF=[];
        NcurrentNF(1)=ncells;
        for j = 1 : ntsteps
            xNF(:,j+1) = muNF + (xNF(:,j)-muNF)*exp(-dt/thetaNF) + sqrt( (1-exp(-2*dt/thetaNF))*sigmaNF^2 ) * dw(:,j);
            yNF = muNF + (yNF-muNF)*exp(-dt/thetaNF) + sqrt( (1-exp(-2*dt/thetaNF))*sigmaNF^2 ) * randn ( size(yNF) );
            prbInv=ProbAmp * 0.01*pchip(xCentExp,ExpLandsc,xNF(:,j+1));
            %prbInv=-ProbAmp *log(1-pchip(xCentExp,ExpLandsc,xNF(:,j+1)));
            rndInv=rand(size(prbInv));
            indInv=find(prbInv>rndInv);
            yNF=[yNF; xNF(indInv,j+1)]; % add next set of invading cells to yNF
            xNF(indInv,j+1)=NaN;        % remove next set of invading cells from xNF
            NcurrentNF(j+1)=sum(isfinite(xNF(:,j+1)));
            NinvNF(j+1)=sum(isfinite(yNF));
        end
%         plot ( t, NcurrentNF, 'b-', 'LineWidth', 2 );
%         plot ( t, NinvNF, 'c-', 'LineWidth', 2 );
        percInv(ctr,rep)=NinvNF(end)/ncells;
        aveInv(ctr,rep)=mean(yNF);
        stdInv(ctr,rep)=std(yNF);
        cvrInv(ctr,rep)=std(yNF)/mean(yNF);
        xNFp=reshape(xNF,numel(xNF),1);xNFp=xNFp(~isnan(xNFp));
        [xNFh(rep,:),xEdges]=histcounts(xNFp,linspace(2.5,5,100),'Normalization','Probability');
        [yNFh(rep,:),yEdges]=histcounts(yNF,linspace(2.5,5,100),'Normalization','Probability');
    end
    
    %plot control and invasion histograms
    figure; hold on;
    xCenters=(xEdges(1:end-1)+xEdges(2:end))/2;yCenters=(yEdges(1:end-1)+yEdges(2:end))/2;
    plot(xCenters,mean(xNFh),'LineWidth',2,'Color',[0.5 0.5 0.5]);
    plot(xCenters,mean(yNFh),'LineWidth',2,'Color',[0.53 0 1]);
    stdEst=muNF*quadfit(muNF);xEst=pdf('Normal',xCenters,muNF,stdEst);plot(xCenters,xEst*max(mean(xNFh))/max(xEst),'r+');
    xExp=pdf('Normal',xCenters,muNF,sigmaNF);plot(xCenters,xExp*max(mean(xNFh))/max(xExp),'k+');
    yExp=pdf('Normal',xCenters,invGFPave(ctr),invGFPave(ctr).*invGFPcv(ctr));plot(xCenters,yExp*max(mean(yNFh))/max(yExp),'m+');
    
    %Experimental histograms
%     plot(xCentExp,sBHexp(currind,:)*max(mean(xNFh))/max(sBHexp(currind,:)),'k--');
%     plot(xCentExp,iBHexp(currind,:)*max(mean(yNFh))/max(iBHexp(currind,:)),'m--');
    
    legend('simSEED','simINV','gftSEED','gxpSEED','gxpINV','expSEED','expINV','Location','NorthWest')
    xlabel ( 'BACH1 level (arb. units)', 'FontSize', 24)
    ylabel ( 'Probability', 'FontSize', 24);
    tstr=sprintf("Dox=%.2f",doxC);
    title ( tstr, 'FontSize', 24 )
    grid on; hold on;
    set(gca,'XLim',[2.5 5],'YLim',[0 0.1], 'FontSize', 24)
    
    aveCtr(ctr)=muNF;
    stdCtr(ctr)=sigmaNF;
    cvrCtr(ctr)=sigmaNF/muNF;
    doxC
    ctr=ctr+1;
end

%prepare a normal distribution to see where it is on landscape
xx=linspace(0,7,1000);

%plot invasion data and interpolated landscape
figure;hold on;
errorbar(GFPave(2:end-1),mInv,stInv,'bo','LineWidth',2);
%plot(xx,pchip(GFPave,[bInv mInv sInv],xx),'k--','LineWidth',2);
%plot(xx,pchip(BHdata(:,1),100*BHdata(:,2),xx),'k--','LineWidth',2);
plot(GFPave(2:end-1),100*mean(percInv'),'c*','LineWidth',1);
plot(xCentExp,100*BHcons,'b--','LineWidth',2);
xlabel('log_{10}(GFP MFI) (arb. units)');ylabel('Invasiveness(%)');
legend('expt.','sim.','inf.','Location','SouthEast');
box on;yMax=70;
set(gca,'FontSize',24,'XLim',[2.5 5.0],'YLim',[0 yMax]);
% text(2.5,yMax+5,'Dox=','FontSize',16,'Color','b');
% %text(GFPave(2)-0.1,52,'0.0','FontSize',16,'Color','b');
% text(GFPave(2)-0.08,yMax+5,'0.0','FontSize',16,'Color','b');
% % text(GFPave(4)-0.05,105,'0.3','FontSize',20);
% text(GFPave(6)-0.08,yMax+5,'0.5','FontSize',16,'Color','b');
% % text(GFPave(7)-0.05,105,'1.0','FontSize',20);
% text(GFPave(8)-0.06,yMax+5,'10','FontSize',16,'Color','b');

%plot invasion data and CHIP interpolated landscape
figure;hold on;errorbar(10.^(GFPave(2:end-1)),mInv,stInv,'bo','LineWidth',2);
plot(10.^(xx),pchip(GFPave,[25 mInv sInv],xx),'b--','LineWidth',2);
%plot(xx,pchip(BHdata(:,1),100*BHdata(:,2),xx),'k--','LineWidth',2);
%plot(GFPave(2:end-1),100*mean(percInv'),'c*','LineWidth',2);
% plot(xCentExp,100*BLcons,'b--','LineWidth',2);
xlabel('GFP MFI (arb. u.)');ylabel('Invasiveness(%)');
%legend('expt.','sim.','inf.','Location','SouthEast');
% box on;yMax=50;
set(gca,'FontSize',24,'XLim',[200 2.5e4],'YLim',[0 yMax]);
% text(2.5,yMax+2,'Dox=','FontSize',16,'Color','b');

%inducer post-invasion
% doxPL=[...
%     0
%     0.1
%     0.3
%     0.35
%     0.5
%     0.6
%     1
%     2
%     10];

%average fluorescence post-invasion
% GFPavePL=[...
%     3.37356879
%     3.533857401
%     3.609787441
%     3.735868215
%     3.785768841
%     3.859890899
%     4.394765724
%     4.352679483
%     4.502049469];
% 
% figure;hold on;
% plot([doxRange(1:end-1) 2],aveCtr,'ko','LineWidth',2);
% errorbar([doxRange(1:end-1) 2],mean(aveInv'),std(aveInv'),'ms','LineWidth',2);
% plot([doxRange(1:end-1) 2],Nmean,'rs','MarkerSize',12);
% plot([doxPL(1:end-1)' 2],GFPavePL,'o','MarkerSize',12,'Color',[0.53 0 1]);
% xlabel ( 'Dox (ng/ml)', 'FontSize', 24)
% ylabel ( 'mean of lg(GFP)', 'FontSize', 24);
% grid on;
% set(gca,'FontSize', 24,'XLim',[0 2.5],'YLim',[3.3 4.65],'XTick',[0 1 2],'XTickLabel',[0 1 10]);
% L1=legend('ctr.','sim.','th.','expt.','Location','EastOutside');
% line([1.2 1.4],[3.3 3.5],'LineWidth',2,'Color','k')
% line([1.3 1.5],[3.3 3.5],'LineWidth',2,'Color','k')
% line([1.2 1.4],[4.45 4.65],'LineWidth',2,'Color','k')
% line([1.3 1.5],[4.45 4.65],'LineWidth',2,'Color','k')
% set(L1,'String',{'ctr.' 'sim.' 'th.' 'expt.' })
% title('mNF, low noise');

%fluorescence CV post-invasion
% GFPcvPL=[...
%     0.0540
%     0.0682
%     0.0744
%     0.0427
%     0.1238
%     0.1140
%     0.0377
%     0.0701
%     0.0611];
% 
% GFPcvPL=GFPcvPL/CVfac;
% 
% figure;hold on;
% plot([doxRange(1:end-1) 2],cvrCtr,'ko','LineWidth',2);
% errorbar([doxRange(1:end-1) 2],mean(cvrInv'),std(cvrInv'),'ms','LineWidth',2);
% plot([doxRange(1:end-1) 2],sqrt(Nvar)./Nmean,'rs','MarkerSize',12);
% plot([doxPL(1:end-1)' 2],GFPcvPL,'o','MarkerSize',12,'Color',[0.53 0 1]);
% xlabel ( 'Dox (ng/ml)', 'FontSize', 24)
% ylabel ( 'CV of lg(GFP)', 'FontSize', 24);
% grid on;
% set(gca,'FontSize', 24,'XLim',[0 2.5],'XTick',[0 1 2],'XTickLabel',[0 1 10])
% L2=legend('ctr.','sim.','th.','expt.','Location','EastOutside');
% line([1.2 1.4],[0.0 0.02],'LineWidth',2,'Color','k')
% line([1.3 1.5],[0.0 0.02],'LineWidth',2,'Color','k')
% line([1.2 1.4],[0.08 0.1],'LineWidth',2,'Color','k')
% line([1.3 1.5],[0.08 0.1],'LineWidth',2,'Color','k')
% set(L2,'String',{'ctr.' 'sim.' 'th.' 'expt.' })
% title('mNF, low noise');

% figure;hold on;
% plot ( t, NcurrentNF, 'b-', 'LineWidth', 2 );
% plot ( t, NinvNF, 'c-', 'LineWidth', 2 );
% xlabel ( 't (hours)', 'FontSize', 20 )
% ylabel ( 'N(t)', 'FontSize', 20, 'HorizontalAlignment', 'right' )
% legend('NF','NFinv');
% set(gca,'FontSize',20);
% 
% figure;
% plot ( t, xNF(1:100,:), 'k-', 'LineWidth', 1 )
% xlabel ( 't (hours)', 'FontSize', 16 )
% ylabel ( 'X(t)', 'FontSize', 16, 'Rotation', 0, 'HorizontalAlignment', 'right' )
% title ( 'Exact simulation of O-U SDE', 'FontSize', 16 )
% grid on;set(gca,'FontSize',20);

% plot O-U simulation time courses
figure;
plot(t,xNF(1:5,:),'LineWidth',2);
xlabel('time (h)');ylabel('X(t)');
title('Exact simulation of O-U SDE');
set(gca,'FontSize',20,'XLim',[ 0 50]);
legend('Cell 1','Cell 2','Cell 3','Cell 4','Cell 5','Location','bestoutside')

