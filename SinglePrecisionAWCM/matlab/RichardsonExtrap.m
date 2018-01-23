clc, clear all

% define maximum level
Jmax = 12;

% define stepsizes (include some shift so h_0 = 0.25)
for j=0:Jmax
    stepsize(j+1) = single( 2^(-(j+1)) );
end

% define some test point
xvec = single( 0:0.0001:1 );

% define some test function f(x) and its derivative df(x)
f = @(x) single( exp(4*x) + cos(x) );
df = @(x) single( 4*exp(4*x) - sin(x) );

for i = 1:length(xvec)

    % current point of evaluation
    x = xvec(i);

    % perform the 2-point and 4-point centered finite difference ...
    % ... schemes for all levels j
    for j=0:Jmax
        h = single( stepsize(j+1) );
        cntrd2pnt(j+1) = single( ( f(x+h) - f(x-h) ) / (2*h) );
        cntrd4pnt(j+1) = single( ( f(x-2*h) - 8*f(x-h) ...
        + 8*f(x+h) - f(x+2*h) ) / ( 12 * h ) );
    end

    % compute richardson extrapolation estimation of first derivative ...
    % ... based on a 2-point formula
    for j = 0:Jmax
        h = single( stepsize(j+1) );
        a = single( RE(0,h,f,x) );
        b = single( RE(0,2*h,f,x) );
        extrap2pnt(j+1) = double( 4*double(a) - double(b) ) / 3;
    end

    % compute the error
    for j=0:Jmax
        resid2pnt(j+1,i) = abs( df(x) - cntrd2pnt(j+1) );
        resid4pnt(j+1,i) = abs( df(x) - cntrd4pnt(j+1) );
        residRE(j+1,i) = abs( df(x) - extrap2pnt(j+1) );
    end

end

% compute norms
for j=0:Jmax
    norm2pnt(j+1) = mean( resid2pnt(j+1,:) );
    norm4pnt(j+1) = mean( resid4pnt(j+1,:) );
    normRE(j+1) = mean( residRE(j+1,:) );
end

% plot error
figure(1)
semilogy(0:Jmax,norm2pnt,'k+--',0:Jmax,norm4pnt,'k.--',...
            0:Jmax,normRE,'ro-');
set(0,'defaulttextinterpreter','latex');
set(gca,'fontsize',16);
xlabel('Simulated Wavelet level $j$')
ylabel('$|f_x(x) - \overline{f}_x(x)|$');
legend('2-pnt','4-pnt','2-pnt RE');
