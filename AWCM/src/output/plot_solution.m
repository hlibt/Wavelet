clc, clear all
set(0,'defaulttextinterpreter','latex')
fs=16;
set(gca,'fontsize',fs-3)

% input data
U=load('u_t1.dat');
Ux=load('ux_t1.dat');
Uxx=load('uxx_t1.dat');

% organize data
x=U(:,1);
u=U(:,2);
ux=Ux(:,2);
uxx=Uxx(:,2);

% true solution
u_true=sin(pi*x);
ux_true=pi*cos(pi*x);
uxx_true=-pi^2*sin(pi*x);

% compute norms
Norm = norm(u_true-u,'inf') / norm(u_true)

% plots
figure(1)
% plot(x,u,'bo-',x,u_true,'r')
plot(x,u,'b-')
% figure(2)
% plot(x,ux,'bo-',x,ux_true,'r')
% figure(3)
% plot(x,uxx,'bo-',x,uxx_true,'r')
