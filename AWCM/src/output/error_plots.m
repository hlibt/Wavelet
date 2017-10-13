clc, clear all

J2=load('J2_N3_ux.dat');
J3=load('J3_N3_ux.dat');
J4=load('J4_N3_ux.dat');
J5=load('J5_N3_ux.dat');
J6=load('J6_N3_ux.dat');

uxtru2=-pi*sin(pi*J2(:,1));
uxtru3=-pi*sin(pi*J3(:,1));
uxtru4=-pi*sin(pi*J4(:,1));
uxtru5=-pi*sin(pi*J5(:,1));
uxtru6=-pi*sin(pi*J6(:,1));

E2=norm(uxtru2-J2(:,2),2)/norm(uxtru2,2);
E3=norm(uxtru3-J3(:,2),2)/norm(uxtru3,2);
E4=norm(uxtru4-J4(:,2),2)/norm(uxtru4,2);
E5=norm(uxtru5-J5(:,2),2)/norm(uxtru5,2);
E6=norm(uxtru6-J6(:,2),2)/norm(uxtru6,2);
Err=[E2 E3 E4 E5 E6];

for i=1:5 
    dx(i)=.125/2^(i-1);
end

figure(1)
loglog(dx,Err)