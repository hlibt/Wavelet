clc, clear all, clf
set(0,'defaulttextinterpreter','latex')
fs=16;
set(gca,'fontsize',fs-3)

% number of steps
n = 6000;

vxn = 2400;
vx = linspace(-1,1,vxn);
vtn = n;
vt = linspace(0,1,vtn);
vu = burgers_solution(0.01,vxn,vx,vtn,vt);

% input data and make movie
figure(1)
for i=0:n-1
    if ( mod(i,5)==0 )
        filename = sprintf('_soln_files/u%d.dat',i);
        U = load(filename);
        plot( U(:,1) , U(:,2), 'b' ); grid on; hold on;
        plot( vx, vu(:,i+1) ,'r'); hold on;
        plot( U(:,1) , 0 , '.b' );
        hold off;
        axis( [ -1 1 -1 1 ] );
        pause(0.001);
    end
end
