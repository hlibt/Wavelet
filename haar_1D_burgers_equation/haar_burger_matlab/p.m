function b = p(x,Xi1,Xi2,Xi3)
    if (x>=Xi1 && x<Xi2 )
        b=x-Xi1;
    elseif (x>=Xi2 && x<Xi3)
        b=Xi3-x;
    else
        b=0;
    end
end