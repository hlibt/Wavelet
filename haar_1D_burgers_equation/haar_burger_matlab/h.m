function a = h(x,Xi1,Xi2,Xi3)
    if (x>=Xi1 && x<Xi2 )
        a=1;
    elseif (x>=Xi2 && x<Xi3)
        a=-1;
    else
        a=0;
    end
end