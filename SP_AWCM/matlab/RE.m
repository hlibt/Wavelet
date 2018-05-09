function approx = RE(j,h,f,x)

    % define factor k
    k = 2 + 2*(j-1);

    % if j==0 compute 2-point centered finite difference
    if j==0
        approx = ( f(x+h) - f(x-h) ) / ( 2 * h );
    else
        approx =  RE(j-1,h/2,f,x) + ( RE(j-1,h/2,f,x) - ...
                  RE(j-1,h,f,x) ) / ( 2^k - 1 );
    end

end
