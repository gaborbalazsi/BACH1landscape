%produces 2D invasion landscape versus the CV and the mean

clearvars; close all;

doxRange=[0.1 0.3 0.5 1.0 10.0];

%inducer concentrations
dox=[...
    -1
    0.1
    0.3
    0.5
    1
    10
    20];

%average fluorescence
GFPaveH=[...
    2.4
    3.488107435
    3.725761448
    3.759890103
    3.967407157
    4.196857051
    5.0];

%GFPave=GFPave/log(10);

GFPcvH=[...
    0.05    
    0.084847147
    0.085538987
    0.08376257
    0.108231423
    0.109985837
    0.16];

%invasiveness (High mNF-BACH1)
InvFH=[...
    25       25     25
    28.17	31.08	24.74
    25.29	30.24	30.41
    18.3	21.59	18.1
    26.48	20.08	25
    43.35	55.02	46.87
    100     100     100];

mInvFH=mean(InvFH')';


%average fluorescence
GFPaveL=[...
    2.4
    3.398626132
    3.648960343
    3.809999739
    4.062188809
    4.322541832
    5.0];

%GFPave=GFPave/log(10);

GFPcvL=[...
    0.05    
    0.068988301
    0.075706523
    0.075070879
    0.079746859
    0.090183528
    0.1];

%invasiveness (Low mNF-BACH1)
InvFL=[...
    25      25      25
    30.92	30.29	33.33
    35.93	31.15	35.86
    13.97	11.42	15.74
    16.62	21.42	22.65
    38.35	36.83	39.37
    100     100     100];

mInvFL=mean(InvFL')';

mInvF=[mInvFL mInvFH];

figure;hold on;
plot3(GFPaveL,GFPcvL,mInvFL,'o','LineWidth',2,'MarkerSize',12);
plot3(GFPaveH,GFPcvH,mInvFH,'ro','LineWidth',2,'MarkerSize',12);
% x=[GFPaveL;GFPaveH;5];
% y=[GFPcvL;GFPcvH;0.06];
% z=[mInvFL;mInvFH;100];
x=[GFPaveL;GFPaveH;4.4;4.4;3.2;3.2];
y=[GFPcvL;GFPcvH;0.09;0.11;0.06;0.1];
z=[mInvFL;mInvFH;60;70;10;10];
XI=3:0.01:5;
YI=0.03:0.001:0.16;
YI=YI';
ZI = griddata(x,y,z,XI,YI,'v4');
surf(XI,YI,ZI);alpha 0.5
shading interp;colormap gray;
xlabel('mean, lg(BACH1)');
ylabel('CV, lg(BACH1)');
zlabel('Invasiveness (%)')
L1=legend('mNF-BACH1-L','mNF-BACH1-H');
set(gca,'FontSize',24,'XLim',[3.22 4.4],'YLim',[0.03 0.16],'ZLim',[0 100],'View',[-15.5604 32.8190]);
%set(gca,'FontSize',24,'View',[-13.0637 41.4925]);
set(L1,'String',{'mNF-BACH1-L','mNF-BACH1-H'})

