clc, clear all, clf

% input data:
U=load('solution.dat');
u=U(:,2);
D=load('derivative.dat');
x=D(:,1);
ux=D(:,2);
<<<<<<< HEAD
xtru=linspace(-1,1,101);
figure(1)
plot(x,u,'bo-',x,ux,'ro-',xtru,cos(pi*xtru),'k',xtru,-pi*sin(pi*xtru),'g')
=======
xtru=linspace(-1,1,4000);
utru=cos(80*pi*xtru).*exp(-64*xtru.^2);
figure(1)
% plot(x,u,'bo-',x,ux,'r+-')
plot(x,u,'b-')
% plot(xtru,utru,'r-')
>>>>>>> fa413aa2a4c1cb7b44ef9ae9e9e0975556d75844
