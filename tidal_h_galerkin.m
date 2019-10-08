clearvars;
% gravity flag
flag = 0; % =0: no gravity; =1: with gravity

back_fd =0; 
forward_fd =0; 
central_fd = 1; 

% universal gravity constant
G = 6.67408e-11; % m3/kg/s2

% planet radius
a = 2000e3; % m
rho = 2840; % kg/m3
Mplanet = 4/3*pi*a^3*rho; % planet mass
gplanet = G*Mplanet/a^2; % planet surface gravity

vs = 1200; % m/s
vp = 3000; % m/s

mu = rho*vs*vs; % shear modulus
lambda = rho*vp*vp - 2 *mu; % Lame constant

B = lambda /mu;

% moon orbit radius
rp = 10*a;
% moon mass
moon = 1e16; % kg
fac = sqrt(5/(4*pi));

A = G*moon/rp^3/fac;

%% define number of parameters
N=3;
Q=zeros(2*N+2,2*N);
b=zeros(2*N+2,1);
%%
Nr=600;
r=linspace(0,a,Nr)/a;
dr=r(2)-r(1);
%%
for k=1:N
    for n=1:N
        c(k,n)=sum(r.^k.*r.^n)*dr;
    end

end
%% construct matrix
krow=0;
for k=1:N
    krow=krow+1;
    for n=1:N
        Q(krow,n)=Q(krow,n)+mu*(-10-2*B+4*n+2*B*n+2*(n-1)*n+B*n*(n-1))*c(k,n);
        Q(krow,n+N)=Q(krow,n+N)+mu*(18+6*B-6*n-6*B*n)*c(k,n);
    end
    %b(krow)=-sum(rho*2*A*r.^3.*r.^k)*dr*a^3;
    b(krow)=-2*A*a^3*rho*c(k,3);
end

for k=1:N
    krow=krow+1;
    for n=1:N
        Q(krow,n)=Q(krow,n)+mu*(4+2*B+n+B*n)*c(k,n);
        Q(krow,n+N)=Q(krow,n+N)+mu*(-12-6*B+2*n+(n-1)*n)*c(k,n);
    end
    %b(krow)=-sum(rho*A*r.^3.*r.^k)*dr*a^3;
    b(krow)=-A*a^3*rho*c(k,3);
end
%BC1
krow=krow+1;
for n=1:N
    Q(krow,n)=2*B+n*(2+B);
    Q(krow,n+N)=-6*B;
end
% BC2
krow=krow+1;
for n=1:N
    Q(krow,n)=1;
    Q(krow,n+N)=-1+n;
end
%%
x=Q\b;
U1=x(1)
U3=x(3)
%%
Ua=sum(x(1:N));
Vnorthpole = G * moon / rp^3 * a^2;
hhydro = Vnorthpole/gplanet;
h2 = Ua / hhydro *fac 

