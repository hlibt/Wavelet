clc, clear all, clf

% input data:
D=load('derivative.dat');
x=D(:,1);
for i=2:length(x)
    h=x(i)-x(i-1)
end
u=D(:,2);
figure(1)
for i=1:length(x)
    plot(x(i),u(i),'bo'), hold on
end