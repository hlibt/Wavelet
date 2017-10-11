clc, clf, clear all

data0=load('coeff0.dat')
data1=load('coeff1.dat')
data2=load('coeff2.dat')
data3=load('coeff3.dat')
data4=load('coeff4.dat')
data5=load('coeff5.dat')
% data6=load('coeff6.dat')
% data7=load('coeff7.dat')
% data8=load('coeff8.dat')

h=figure(1);
axis([-1 1 0 9]); hold on

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
% h=plot(data7(:,1),data7(:,2),'ko'), hold on
% set(h(1),'MarkerEdgeColor','k','MarkerFaceColor','k')
% h=plot(data8(:,1),data8(:,2),'ko'), hold on
% set(h(1),'MarkerEdgeColor','k','MarkerFaceColor','k')
set(gca,'ygrid','on')
set(gca,'Ytick',0:1:9);

ylabel('$j$','Interpreter','LaTex')
xlabel('$x_{k}^{j}$','Interpreter','Latex')

  