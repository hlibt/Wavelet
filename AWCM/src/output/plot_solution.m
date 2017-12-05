clc, clear all
set(0,'defaulttextinterpreter','latex')
fs=16;
set(gca,'fontsize',fs-3)

% number of steps
n = 400;

% input data and make movie
figure(1)
for i=1:n
    filename = sprintf('_soln_files/u%d.dat',i-1);
    U = load(filename);
    plot( U(:,1) , U(:,2) );
end
