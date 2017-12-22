clc, clear all, clf
set(0,'defaulttextinterpreter','latex')
fs=16;
set(gca,'fontsize',fs-3)

% number of steps
n = 8000;

%vxn = 1000;
%vx = linspace(-1,1,vxn);
%vtn = n;
%vt = linspace(0,1,vtn);
%vu = burgers_solution(0.01,vxn,vx,vtn,vt);

% input data and make movie
figure(1)
for i=0:n-1
    if ( mod(i,2)==0 )
        filename = sprintf('_soln_files/u%d.dat',i);
        U = load(filename);
        [x,I] = sort( U(:,1) );
        u = U(:,2);
        plot( x , u(I), 'b' ); grid on; hold on;
%        plot( vx, vu(:,i+1) ,'r'); hold on;
        plot( x , 0, '.b' );
        hold off;
        axis( [ -1 1 -1 1 ] );
        pause(0.0000001);
    end
end
