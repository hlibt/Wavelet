clc, clear all
set(0,'defaulttextinterpreter','latex')
fs=16;
set(gca,'fontsize',fs-3)

% number of steps
n = 1000;

% input data and make movie
figure(1)
for i=n:n
    filename = sprintf('_soln_files/u%d.dat',i-1);
    U = load(filename);
    plot( U(:,1) , U(:,2) );
    axis( [ -1 1 0 1 ] );
    pause(0.000001);
end
