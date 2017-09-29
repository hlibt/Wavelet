clc, clear all, clf

% input data:
D=load('solution.dat');
x=D(:,1);
u=D(:,2);
figure(1)
% utrue=-pi*sin(pi*x);
plot(x,u,'b')
% axis([-1 1 -1.2 1.2])