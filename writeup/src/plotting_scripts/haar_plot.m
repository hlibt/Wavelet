clc, clear all, clf

% define domain
x = linspace(-.2,1.2,200);

% compute haar scaling function
for i=1:length(x)
    phi(i) = ( x(i) >= 0 & x(i) <= 1 );
    psi(i) = ( x(i) >= 0 & x(i) <=.5 ) - ( x(i) >.5 & x(i) <= 1 );
end

A = [x; phi];
B = [x; psi];

% output data into file to be fed into gnuplot
fileID = fopen('haar_scaling.dat','w');
fprintf(fileID,'%12.8f %12.8f\n',A);
fclose(fileID);

fileID = fopen('haar_wavelet.dat','w');
fprintf(fileID,'%12.8f %12.8f\n',B);
fclose(fileID);

%{
% plot the wavelet function
subplot(122)
h = plot(x,psi,'k'); 
set(h,'Linewidth',4);
axis( [-.2 1.2 -1.2 1.2] );
xlabel('$x$','fontsize',30);
ylabel('$\psi(x)$','fontsize',30);
h.XTicks = [0 1];
h.YTicks = [-1 1];
set(gca,'fontsize',30);
box off %}
