clc, clear all
set(0,'defaulttextinterpreter','latex')
fs=16;
set(gca,'fontsize',fs-3)

% input data
U0=load('_soln_files/u.dat');
U1=load('_soln_files/u1.dat');

% organize data
x=U0(:,1);
u=U0(:,2);
u1=U1(:,2);

% plots
figure(1)
plot(x,u,'bo-',x,u1,'r-');
