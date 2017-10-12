clc, clear all, clf

% input data:
U=load('solution.dat');
u=U(:,2);
D=load('derivative.dat');
x=D(:,1);
ux=D(:,2);
xtru=linspace(-1,1,4000);
utru=cos(80*pi*xtru).*exp(-64*xtru.^2);
figure(1)
% plot(x,u,'bo-',x,ux,'r+-')
plot(x,u,'b-')
% plot(xtru,utru,'r-')