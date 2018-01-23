clc, clear all

% load data output from c++ program
data = load('error.dat');

% plot the data
loglog(data(:,1),data(:,2),'ks--',data(:,1),data(:,3),'k^--',data(:,1),data(:,4),'r');
