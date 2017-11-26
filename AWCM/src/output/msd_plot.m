clc, clear all, clf
data = load('src/msd_output.dat');
fs = 19;
lw = 3;
figure(1)
set(gca,'Fontsize',fs);
set(0,'defaulttextinterpreter','latex')
plot( data(:,1) * 0.001, data(:,2), 'k', 'Linewidth',lw), grid on
xlabel('$t$ (fs)','Fontsize',fs+5), ylabel('$\langle r(t)^2 \rangle$ ($\sigma^2$)','Fontsize',fs+5)
matlab2tikz()