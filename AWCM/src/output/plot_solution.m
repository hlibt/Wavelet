clc, clear all
set(0,'defaulttextinterpreter','latex')
fs=16;
set(gca,'fontsize',fs-3)

% number of steps
n = 9000;

% input data and make movie
figure(1)
for i=0:n
    if mod(i,20)==0
        filename = sprintf('_soln_files/u%d.dat',i);
        U = load(filename);
        plot( U(:,1) , U(:,2) ); grid on; hold on;
        plot( U(:,1) , U(:,2) , '.r' );
        hold off;
        axis( [ -1 1 0 1 ] );
        pause(0.001);
    end
end
