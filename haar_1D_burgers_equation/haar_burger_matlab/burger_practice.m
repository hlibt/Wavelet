%-------------------------------------------------------------------%
%                                                                   %
% SOLUTION TO 1-D BURGERS EQUATION USING THE HAAR WAVELET BASIS     %
%                                                                   %
% AUTHOR : BRANDON GUSTO                                            %
% CREATED: JULY 05, 2017                                            %
%                                                                   %
%-------------------------------------------------------------------%

% Clear data:
clc, clear all

% Initial parameters:
v=0.1;            % Coefficient of kinematic viscosity
M=2;             % Number (half?) of collocation points
J=log2(M);         % Maximum scale value

% Temporal discretization:
N=10;             % Number of timesteps
t_i=0;             % Initial time
t_f=1;             % Final time
dt=(t_f-t_i)/N;    % Step size
t=t_i:dt:t_f;

% Initialize matrices:
A=zeros(2*M,2*M);
H=zeros(2*M,2*M);
P=zeros(2*M,2*M);
Q=zeros(2*M,2*M);
Q_tilda=zeros(2*M,2*M);
b=zeros(2*M,1);
% c=zeros(2*M,1);
U=zeros(N,2*M);
Ux=zeros(N,2*M);
Uxx=zeros(N,2*M);

% Spatial discretization:
% % % % dx1=1/M;
% % % % dx2=1/(3*M+.5);
% % % % for n=1:M/2
% % % %     x(n)=dx1*n;
% % % %     xstar=x(n);
% % % % end
% % % % for n=1:3*M/2
% % % %     x(n+M/2)=xstar+dx2*n;
% % % % end
% % n=1:2*M;
% % x=(n-0.5)/(2*M);  
j=1;
for n=0:2^(j+1)-1
    x(n+1,j+1)=n*2^-(j+1);
end
x

% Boundary and initial conditions:
% f=@(x) sin(2*pi*x);                     % u(x,0)
% f_prime=@(x) 2*pi*cos(2*pi*x);          % u'(x,0)
% f_dblprime=@(x) -4*pi^2*sin(2*pi*x);    % u''(x,0)
f=@(x) 4*x*(1-x);                     % u(x,0)
fx=@(x) 4-8*x;          % u'(x,0)
fxx=@(x) -8;    % u''(x,0)
f1=@(t) 0;                              % u(0,t)
f1t=@(t) 0;    
f2=@(t) 0;                              % u(1,t)
f2t=@(t) 0;
for l=1:2*M
    U(1,l)=f(x(l));
    Ux(1,l)=fx(x(l));
    Uxx(1,l)=fxx(x(l));
end

% Advance simulation in time:
for s=1:N-1
    dt*(s-1)
    
    % Populate 'A' with scaling function at each collocation point (i=1):
    for l=1:2^(1)-1
        x(l,1)=l*2^(-(j+1));
        A(l,1)=q1(x(l))-x(l)*q_tilda(1,0)+dt*(-v*h1(x(l))+...
        Ux(s,l)*(q1(x(l))-x(l)*q_tilda(1,0))+...
        U(s,l)*(p1(x(l))-q_tilda(1,0))); 
        H(l,1)=h1(x(l));
        P(l,1)=p1(x(l));
        Q(l,1)=q1(x(l));
        Q_tilda(l,1)=q_tilda(1,0);
    end
    
    % Solve for c1 here?
    
    % Populate 'A' with first detail wavelets:
       m=2^0;
       for k=0:m-1
           Xi1=k/m;
           Xi2=(k+.5)/m3=;
           Xi3=(k+1)/m;
           i=1+m+k;
           for l=1:2^(j+1)-1
               x(l,1)=l*2^(-(j+1));
               A(l,i)=q(x(l),Xi1,Xi2,Xi3,m)-x(l)*q_tilda(i,m)+...
                   dt*(-v*h(x(l),Xi1,Xi2,Xi3)+Ux(s,l)*(q(x(l),Xi1,Xi2,Xi3,m)-x(l)*q_tilda(i,m))+...
                   U(s,l)*(p(x(l),Xi1,Xi2,Xi3)-q_tilda(i,m)));  
               H(l,i)=h(x(l),Xi1,Xi2,Xi3);
               P(l,i)=p(x(l),Xi1,Xi2,Xi3);
               Q(l,i)=q(x(l),Xi1,Xi2,Xi3,m);
               Q_tilda(l,i)=q_tilda(i,m);
           end
       end
    end
    
% % % % % % % % %     % Populate 'A' with remaining functions:
% % % % % % % % %     for j=0:0
% % % % % % % % %        m=2^j;
% % % % % % % % %        for n=1:2^(j+1)-1 
% % % % % % % % %            x(n,j+1)=n*2^-(j+1);
% % % % % % % % %        end
% % % % % % % % %        for k=0:m-1
% % % % % % % % %            Xi1=k/m;
% % % % % % % % %            Xi2=(k+.5)/m;
% % % % % % % % %            Xi3=(k+1)/m;
% % % % % % % % %            i=1+m+k;
% % % % % % % % %            for l=1:2^(j+1)-1
% % % % % % % % %                x(l,1)=l*2^(-(j+1));
% % % % % % % % %                A(l,i)=q(x(l),Xi1,Xi2,Xi3,m)-x(l)*q_tilda(i,m)+...
% % % % % % % % %                    dt*(-v*h(x(l),Xi1,Xi2,Xi3)+U_prime(s,l)*(q(x(l),Xi1,Xi2,Xi3,m)-x(l)*q_tilda(i,m))+...
% % % % % % % % %                    U(s,l)*(p(x(l),Xi1,Xi2,Xi3)-q_tilda(i,m)));  
% % % % % % % % %                H(l,i)=h(x(l),Xi1,Xi2,Xi3);
% % % % % % % % %                P(l,i)=p(x(l),Xi1,Xi2,Xi3);
% % % % % % % % %                Q(l,i)=q(x(l),Xi1,Xi2,Xi3,m);
% % % % % % % % %                Q_tilda(l,i)=q_tilda(i,m);
% % % % % % % % %            end
% % % % % % % % %        end
% % % % % % % % %     end
    
    % Construct rhs of system with the vector 'b':
    for l=1:2*M
        b(l)=-f1t(t(s+1))-x(l)*(f2t(t(s+1))-f1t(t(s+1)))...
            +v*Uxx(s,l)-Ux(s,l)*(f1(t(s))-f1(t(s+1)))...
            +x(l)*(-f2(t(s+1))+f1(t(s+1))-f2(t(s))+f1(t(s)))*Ux(s,l)...
            +U(s,l)*(f1(t(s+1))+f1(t(s))-f2(t(s))-Ux(s,l));
    end
     
% % %     for i=1:size(A)
% % %         D(i)=det(A(1:i,1:i))
% % %     end

    % Solve for wavelet coefficients in vector 'c': 
    c=mldivide(A,b)';
    
    % Advance solution in timei
    for l=1:2*M
        Uxx(s+1,l)=dt*c*H(l,:)'+Uxx(s,l);
        Ux(s+1,l)=dt*c*(P(l,:)-Q_tilda(l,:))'+f2(t(s+1))-f1(t(s+1))...
            +f1(t(s))-f2(t(s))+Ux(s,l);
        U(s+1,l)=dt*c*(Q(l,:)-x(l)*Q_tilda(l,:))'+x(l)*(f2(t(s+1))...
            -f1(t(s+1))-f2(t(s))+f1(t(s)))+f1(t(s+1))+U(s,l)-f1(t(s));
    end
end

% % % % Fix ends:
% % % xnew=zeros(1,2*M+2);
% % % for l=1:2*M
% % %     xnew(l+1)=x(l);
% % % end
% % % xnew(1)=0;xnew(2*M+2)=1;
% % % Unew=zeros(N,2*M+2);
% % % for s=1:N
% % %     for l=1:2*M;
% % %     Unew(s,l+1)=U(s,l);
% % %     end
% % %     Unew(s,1)=0;
% % %     Unew(s,2*M+2)=0;
% % % end
% % % 
% % % % Plot movie:
% % % v = VideoWriter('movie.avi');
% % % open(v);
% % % for i = 1:N
% % %     % Exact solution and solver solution:
% % %     plot(xnew,Unew(i,:));
% % %     axis([-1 1.5 -1 1]);
% % %     
% % %     % Store the frame
% % %     frame=getframe(gcf); 
% % %     writeVideo(v,frame);
% % % end
% % % close(v);