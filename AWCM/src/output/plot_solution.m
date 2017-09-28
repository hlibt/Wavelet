clc, clear all, clf

% input data:
D=load('solution.dat');
j=1;
for i=1:length(D)
    if (mod(i,2)==0)
       ux(j)=D(i,2);
       x(j)=D(i,1);
       j=j+1;
    end
end
soln=-pi*sin(pi*x);
figure(1)
plot(x,ux,'b',x,soln,'r')