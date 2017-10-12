clc, clear all, clf

% input data:
U=load('solution.dat');
u=U(:,2);
D=load('derivative.dat');
x=D(:,1);
ux=D(:,2);
xtru=linspace(-1,1,101);
figure(1)
plot(x,u,'bo-',x,ux,'ro-',xtru,cos(pi*xtru),'k',xtru,-pi*sin(pi*xtru),'g')