clc, clear all, clf
set(0,'defaulttextinterpreter','latex')
fs=16;
set(gca,'fontsize',fs-3)

% input data
u_tru=load('u_tru.dat');
d0=load('u_eps_2.dat');
d1=load('u_eps_3.dat');
d2=load('u_eps_4.dat');
d3=load('u_eps_6.dat');
d4=load('u_eps_8.dat');
d5=load('u_eps_9.dat');
d6=load('u_eps_10.dat');
d7=load('u_eps_12.dat');

% organize data
x=u_tru(:,1);
u=u_tru(:,2);

% true solution
x0=1/3;
v=10^(-2);
u_vtrue=-tanh((x+x0)/(2*v)) + exp(-64^2*(x-x0).^2);


% compute norms
Norm=[  norm(u-d0(:,2),'inf')
        norm(u-d1(:,2),'inf')
        norm(u-d2(:,2),'inf')
        norm(u-d3(:,2),'inf')
        norm(u-d4(:,2),'inf')
        norm(u-d5(:,2),'inf')
        norm(u-d6(:,2),'inf')
        norm(u-d7(:,2),'inf') ]';

activPnt=[ 32 52 89 259 760 1290 2332 4874 ];

figure(1)
loglog(activPnt,Norm), grid on
xlabel('Number of points kept, $N$','fontsize',fs)
ylabel('$||f^J(x) - f_{\geq}^{J}||_{\infty}$','fontsize',fs)
