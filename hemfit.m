clearvars; close all;
h=[0 10 20 30 50 80];

% bgr=mean([664, 672, 594])-190;
bgr=mean([410, 390, 434]);


mbLd00=[...
    1013	1020	952
    711	783	739
    729	458	503
    455	454	457
    468	462	463
    470	470	469];

mbLd10=[...
    17336	17448	17491
    6072	6072	6132
    4007	3922	3938
    3389	3356	3408
    3350	3267	3336
    3038	3114	3063];

mbHd00=[...
    825	832	827
    547	538	545
    498	506	508
    482	485	482
    471	473	474
    473	481	484];

mbHd10=[...
    18610	18518	18243
    6019	5871	5879
    4142	4093	4064
    3714	3498	3582
    3336	3249	3405
    3222	3230	3149];

figure;hold on;
hh=0:80;

mbLd00a=mean(mbLd00,2)-bgr;
mbLd00s=std(mbLd00,0,2);
errorbar(h,mbLd00a,mbLd00s,'bv','MarkerSize',12,'LineWidth',2);

mbLd10a=mean(mbLd10,2)-bgr;
mbLd10s=std(mbLd10,0,2);
errorbar(h,mbLd10a,mbLd10s,'cv','MarkerSize',12,'LineWidth',2);

mbHd00a=mean(mbHd00,2)-bgr;
mbHd00s=std(mbHd00,0,2);
errorbar(h,mbHd00a,mbHd00s,'r^','MarkerSize',12,'LineWidth',2);

mbHd10a=mean(mbHd10,2)-bgr;
mbHd10s=std(mbHd10,0,2);
errorbar(h,mbHd10a,mbHd10s,'m^','MarkerSize',12,'LineWidth',2);

a =          40;
b =       9.286;
% a =       32.48;
% b =       2.871;
plot(hh,mbLd00a(1)*(a+hh)./(a+b*hh),'b--','LineWidth',2);

a =       27.56;
b =       9.046;
% a =       23.16;
% b =       7.602;
plot(hh,mbLd10a(1)*(a+hh)./(a+b*hh),'c--','LineWidth',2);

a =       26.13;
b =       9.052;
% a =       5.225;
% b =        1.81;
plot(hh,mbHd00a(1)*(a+hh)./(a+b*hh),'r--','LineWidth',2);

a =       21.82;
b =       8.562;
% a =       18.67;
% b =       7.327;
plot(hh,mbHd10a(1)*(a+hh)./(a+b*hh),'m--','LineWidth',2);

xlabel('[hemin] (\muM)');ylabel('BACH1 MFI (arb. u.)');
set(gca,'FontSize',24,'YScale','Log','YLim',[1 1e10]);
legend('BL Dox=0','BL Dox=10','BH Dox=0','BH Dox=10','Location','NorthEast')
title('MB231-mNF-BACH1')
