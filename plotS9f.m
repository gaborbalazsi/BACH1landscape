% Produces Ext. figure panels S9f 

clearvars; close all;

%same as ouBLindLandsc, but different x bins.

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
thetaNF = 1000.0;
ProbAmp=0.25;
seed='shuffle';

doxRange=[    0    0.1    0.3    0.35    0.5    0.6    1    2    10];

%binned ladscape estimates
BLdata=load('LandscIndBL.txt'); %Landscape data for each individual Dox
%BLcons=load('aveLandscProdCut13BL.txt');
%BLcons=load('aveLandscBL.txt');
BLcons=load('smLandscBL.txt');
xCentExp=load('xCentersBL.txt');

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
    2
    10
    20];

%average fluorescence of Seeded BL cells (padded at both ends)
GFPave=[...
    2.0
    2.979016133
    3.398626132
    3.648960343
    3.788605752
    3.809999739
    3.823816377
    4.062188809
    4.202669121
    4.322541832
    10];

%CV of Seeded BL cells (padded at both ends)
GFPcv=[...
    0.05    
    0.0670
    0.0690
    0.0757
    0.0637
    0.0751
    0.0764
    0.0797
    0.0802
    0.0902
    0.1];

GFPstd=GFPcv.*GFPave/CVfac;

%invasiveness (Low mNF-BACH1)  (padded at both ends with bInv, sInv)
bInv=0;
sInv=100;
InvFL=[...
    bInv    bInv    bInv
    27.55 	25.66 	24.27 
    30.92 	30.29 	33.33 
    35.93 	31.15 	35.86 
    30.39 	26.63 	25.13 
    13.97 	11.42 	15.74 
    16.25 	17.58 	18.06 
    16.62 	21.42 	22.65 
    23.79 	29.70 	31.43 
    38.35 	36.83 	39.37 
    sInv    sInv    sInv];

avInv=mean(InvFL(2:end-1,:)');   % mean of 3 invasion replicates, unpadded
stInv=std(InvFL(2:end-1,:)');  % std of 3 invasion replicates, unpadded
minInv=min(avInv);  % overall landscape miminmum, unpadded
nInv=AmpFac*(avInv-minInv);  %normalized landscape with minimum reset to = 0
%mInv=nInv;  %measured Invsiveness, unpadded
mInv=avInv;


%average GFP fluorescence of Invaded BL cells (not padded)
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

%CV of Invaded BL cells (not padded)
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

GFPstd=GFPcv.*GFPave/CVfac;


%these will be the simulated values
percInv=zeros(9,10);
aveInv=zeros(9,10);
stdInv=zeros(9,10);
cvrInv=zeros(9,10);

xx=linspace(0,11,1000);

ctr=1;
for doxC=doxRange(end)
%current Dox concentration
    indD=find(dox==doxC);   %where on the landscape?
    currind=indD-1;
    gfpC=GFPave(indD);      %what mean expression level?
    sigmaNF=GFPstd(indD);      %what std expression level?
    muNF = gfpC;
    varNF=sigmaNF*sigmaNF;
    varChi2NF=2*varNF*(varNF+2*muNF^2);
    
    ExpLandsc=BLdata(currind,:);    %individual landscape
    %ExpLandsc=BLcons;              %consolidated landscape
    
    %read in actual data
    xEdges=linspace(2.5,5,50);yEdges=xEdges;
    sheetname=sprintf('%d',currind);
    ExpData=xlsread('./data231/MB231_1.1BLdataGB.xlsx',sheetname);
    ExpSeeded=ExpData(:,1:3);ExpSeeded=ExpSeeded(~isnan(ExpSeeded(:)));
    [sBLexp(currind,:),~]=histcounts(ExpSeeded,xEdges,'Normalization','pdf');
    ExpInvded=ExpData(:,4:6);ExpInvded=ExpInvded(~isnan(ExpInvded(:)));
    [iBLexp(currind,:),~]=histcounts(ExpInvded,xEdges,'Normalization','pdf');
    iBLexp(currind,:)=mInv(currind).*iBLexp(currind,:);

    ExpSeeded=muNF+sigmaNF*randn(21000,1);

    figure;hold on;

    NMAX=ncells*2;   %maximum number of cells
    growthRate=0.0;

    % simulate NF
    for rep=1:1
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
            plot ( t, NcurrentNF, 'b-', 'LineWidth', 2 );
            plot ( t, NinvNF, 'c-', 'LineWidth', 2 );
        percInv(ctr,rep)=NinvNF(end)/ncells;
        aveInv(ctr,rep)=mean(yNF);
        stdInv(ctr,rep)=std(yNF);
        cvrInv(ctr,rep)=std(yNF)/mean(yNF);
        xNFp=reshape(xNF,numel(xNF),1);xNFp=xNFp(~isnan(xNFp));
        [xNFh(rep,:),xEdges]=histcounts(xNFp,linspace(2.5,5,100),'Normalization','Probability');
        [yNFh(rep,:),yEdges]=histcounts(yNF,linspace(2.5,5,100),'Normalization','Probability');
    end
      
    %experimental pre-invasion values
    aveCtr(ctr)=muNF;
    stdCtr(ctr)=sigmaNF;
    cvrCtr(ctr)=sigmaNF/muNF;
    doxC
    ctr=ctr+1;
end
xlabel('t (hours)');
ylabel('N(t)')
legend('NF','NFinv')
set(gca,'FontSize',20,'XLim',[0 60],'YLim',[0 5000]);
