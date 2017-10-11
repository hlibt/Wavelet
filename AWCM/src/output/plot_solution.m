clc, clear all, clf

% input data:
U=load('solution.dat');
u=U(:,2);
D=load('derivative.dat');
x=D(:,1);
ux=D(:,2);
figure(1)
plot(x,u,'bo-',x,ux,'r+-')