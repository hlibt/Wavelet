% Exact solution:
t2=0;
for s=1:N
    t2=t2+dt;
    x2=0;
    for j=1:101
        sum1=0;
        sum2=0;
        for i=1:100
            A0=@(y) exp(-1/(3*v).*(3*y.^2-2*y.^3));
            An=@(y) exp(-1/(3*v).*(3*y.^2-2*y.^3)).*cos(i*pi*y);
            sum1=sum1+integral(An,0,1)*exp(-i^2*pi^2*v*t2)*i*sin(i*pi*x2);
            sum2=sum2+integral(An,0,1)*exp(-i^2*pi^2*v*t2)*cos(i*pi*x2);
        end
        soln(j,s)=2*pi*v*sum1/(integral(A0,0,1)+sum2);
        x2=x2+0.01;
    end
end
xtru=0:0.01:1;