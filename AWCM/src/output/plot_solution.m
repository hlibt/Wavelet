clc, clear all, clf

% input data:
U=load('solution.dat');
u=U(:,2);
D=load('derivative.dat');
x=D(:,1);
ux=D(:,2);
xtru=linspace(-1,1,400);
utru=-pi^2*cos(pi*xtru);
figure(1)
plot(x,ux,'b+-',xtru,utru,'r'), hold on
plot(x,u,'bo-'), hold on