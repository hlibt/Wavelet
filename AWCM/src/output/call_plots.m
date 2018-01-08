clc, clear all

% number of timesteps
n = 14000;

% plot or do not plot the solution
psoln = true;

% plot or do not plot the coefficients
pcoeffs = false;

% plotting frequency
fqz = 20;

% call function to produce plots
post_process(n,psoln,pcoeffs,fqz)
