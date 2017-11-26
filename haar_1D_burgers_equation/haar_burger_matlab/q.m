function c = q(x,Xi1,Xi2,Xi3,m)
    if (x>=Xi1 && x<Xi2 )
        c=.5*(x-Xi1)^2;
    elseif (x>=Xi2 && x<Xi3)
        c=1/(4*m^2)-.5*(Xi3-x)^2;
    elseif (x>=Xi3 && x<1)
        c=1/(4*m^2);
    else
        c=0;
    end