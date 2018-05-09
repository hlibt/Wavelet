clc, clear all

% define maximum level
Jmax = 12;

% define stepsizes (include some shift so h_0 = 0.25)
for j=0:Jmax
    stepsize(j+1) = single( 2^(-(j+2)) );
end

% define some test point
x = single( 1.9 );

% define some test function f(x) and its derivative df(x)
f = @(x) single( cos( pi * x ) );
df = @(x) single( -pi * sin( pi * x ) );

% perform the 2-point centered finite difference ...
% ... scheme for all levels j
for j=0:Jmax
    h = single( stepsize(j+1) );
    cntrd2pnt(j+1) = single( ( f(x+h) - f(x-h) ) / (2*h) );
end

% perform the 4-point centered finite difference ...
% ... scheme for all levels j
for j=0:Jmax
    h = single( stepsize(j+1) );
    cntrd4pnt(j+1) = single( ( f(x-2*h) - 8*f(x-h) ...
    + 8*f(x+h) - f(x+2*h) ) / ( 12 * h ) );
end

% compute richardson extrapolation not using different ...
% ... levels but only a 4-point stencil on j-1 ...
% ... and compute the finite difference using a 2-pnt FD
t = 2; % refer to wiki page for richardson extrapolation
for j=0:Jmax
    h = single( stepsize((j+1<=3)*(j+1)+(j+1>3)*3) );
    extrap2pnt(j+1) = ( t^2*single( ( f(x+h/t) - f(x-h/t) ) / (2*h/t) ) ...
                      - single( ( f(x+h) - f(x-h) ) / (2*h) ) ) / (t^2 - 1);
end

% compute the error
for j=0:Jmax
    resid2pnt(j+1) = abs( df(x) - cntrd2pnt(j+1) );
    resid4pnt(j+1) = abs( df(x) - cntrd4pnt(j+1) );
    residRE(j+1) = abs( df(x) - extrap2pnt(j+1) );
end


% plot error
figure(1)
semilogy(0:Jmax,resid2pnt,'k+--',0:Jmax,resid4pnt,'k.--',...
            0:Jmax,residRE,'ro-');
set(0,'defaulttextinterpreter','latex');
set(gca,'fontsize',16);
xlabel('Simulated Wavelet level $j$')
ylabel('$|f_x(x) - \overline{f}_x(x)|$');
legend('2-pnt','4-pnt','2-pnt RE');
