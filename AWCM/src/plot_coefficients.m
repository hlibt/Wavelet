clc, clf, clear all

data0=load('solution0.dat')
data1=load('solution1.dat')
data2=load('solution2.dat')
data3=load('solution3.dat')
data4=load('solution4.dat')
data5=load('solution5.dat')
data6=load('solution6.dat')


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
h=plot(data6(:,1),data6(:,2),'ko'), hold on
set(h(1),'MarkerEdgeColor','k','MarkerFaceColor','k')
set(gca,'ygrid','on')

xlabel('$d_{k}^{j}$','Interpreter','LaTex')

  