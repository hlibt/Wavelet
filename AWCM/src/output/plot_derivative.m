clc, clear all, clf

data0=load('derivative0.dat');
data1=load('derivative1.dat');
data2=load('derivative2.dat');
data3=load('derivative3.dat');
data4=load('derivative4.dat');
data5=load('derivative5.dat');
data6=load('derivative6.dat');
data7=load('derivative7.dat');

figure(1)
for i=1:length(data0)
    plot(data0(i,1),data0(i,2)), hold on
end
for i=1:length(data1)
    plot(data1(i,1),data1(i,2),'ms'), hold on
end
for i=1:length(data2)
    plot(data2(i,1),data2(i,2),'bo'), hold on
end
for i=1:length(data3)
    plot(data3(i,1),data3(i,2),'r+'), hold on
end
for i=1:length(data4)
    plot(data4(i,1),data4(i,2)), hold on
end
for i=1:length(data5)
    plot(data5(i,1),data5(i,2),'r'), hold on
end
for i=1:length(data6)
    plot(data6(i,1),data6(i,2),'r'), hold on
end
for i=1:length(data7)
    plot(data7(i,1),data7(i,2),'r'), hold on
end