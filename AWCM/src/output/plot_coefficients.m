clc, clf, clear all

data0=load('coeff0.dat')
data1=load('coeff1.dat')
data2=load('coeff2.dat')
data3=load('coeff3.dat')
data4=load('coeff4.dat')
data5=load('coeff5.dat')
% data6=load('coeff6.dat')


h=figure(1)
axis([-1 1 0 7]); hold on
h=plot(data0(:,1),data0(:,2),'ko'), hold on
set(h(1),'MarkerEdgeColor','k','MarkerFaceColor','k')
h=plot(data1(:,1),data1(:,2),'ko'), hold on
set(h(1),'MarkerEdgeColor','k','MarkerFaceColor','k')
h=plot(data2(:,1),data2(:,2),'ko'), hold on
set(h(1),'MarkerEdgeColor','k','MarkerFaceColor','k')
h=plot(data3(:,1),data3(:,2),'ko'), hold on
set(h(1),'MarkerEdgeColor','k','MarkerFaceColor','k')
h=plot(data4(:,1),data4(:,2),'ko'), hold on
set(h(1),'MarkerEdgeColor','k','MarkerFaceColor','k')
h=plot(data5(:,1),data5(:,2),'ko'), hold on
set(h(1),'MarkerEdgeColor','k','MarkerFaceColor','k')
% h=plot(data6(:,1),data6(:,2),'ko'), hold on
% set(h(1),'MarkerEdgeColor','k','MarkerFaceColor','k')
set(gca,'ygrid','on')

xlabel('$d_{k}^{j}$','Interpreter','LaTex')