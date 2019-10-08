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
vp = 3000*10; % m/s

mu = rho*vs*vs; % shear modulus
lambda = rho*vp*vp - 2 *mu; % Lame constant

B = lambda /mu;

% moon orbit radius
rp = 10*a;
% moon mass
moon = 1e16; % kg
fac = sqrt(5/(4*pi));

A = G*moon/rp^3/fac;
gra_flag=1;
%% define number of parameters
N=6;
Q=sparse(zeros(3*N+3,3*N));
b=zeros(3*N+3,1);
%%
Nr=3000;
r=linspace(0,a,Nr)/a;
dr=r(2)-r(1);
%% 
c=zeros(N+2);
for k=1:N+2
    for n=1:N+2
        c(k,n)=sum(r.^k.*r.^n)*dr;
        %c(k,n)=1/(n+k);
    end

end
%% construct matrix
krow=0;
%eq1
for k=1:N
    krow=krow+1;
    for n=1:N
        Q(krow,n)=Q(krow,n)+mu*(-10-2*B+4*n+2*B*n+2*(n-1)*n+B*n*(n-1))*c(k,n);
        Q(krow,n+N)=Q(krow,n+N)+mu*(18+6*B-6*n-6*B*n)*c(k,n);
        if gra_flag==1
            Q(krow,n)=Q(krow,n)+4/3*G*pi*a^2*rho^2*c(k,n+2);
            Q(krow,n+N)=Q(krow,n+N)-8*G*pi*a^2*rho^2*c(k,n+2);
            Q(krow,n+2*N)=Q(krow,n+2*N)+n*a*rho*c(k,n+1);
        end
    end
    %b(krow)=-sum(rho*2*A*r.^3.*r.^k)*dr*a^3;
    b(krow)=-rho*A*2*a^3*c(k,3);
end
%eq2
for k=1:N
    krow=krow+1;
    for n=1:N
        Q(krow,n)=Q(krow,n)+mu*(4+2*B+n+B*n)*c(k,n);
        Q(krow,n+N)=Q(krow,n+N)+mu*(-12-6*B+2*n+(n-1)*n)*c(k,n);
        if gra_flag==1
            Q(krow,n)=Q(krow,n)-4/3*G*pi*a^2*rho^2*c(k,n+2);
            Q(krow,n+2*N)=Q(krow,n+2*N)+a*rho*c(k,n+1);
        end
    end
    b(krow)=-rho*A*a^3*c(k,3);
end
%eq3
if gra_flag==1
    for k=1:N
        krow=krow+1;
        for n=1:N
            Q(krow,n+2*N)=Q(krow,n+2*N)+(-6+2*n+n*(n-1))*c(k,n);
            Q(krow,n)=Q(krow,n)+(-8*G*pi*a*rho-4*G*pi*a*n*rho)*c(k,n+1);
            Q(krow,n+N)=Q(krow,n+N)+24*G*pi*a*rho*c(k,n+1);
        end
        
    end
end
% BC1
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
%BC3

if gra_flag==1
    krow=krow+1;
    for n=1:N
        Q( krow ,n+2*N)=Q( krow ,n+2*N)+3;
        Q( krow ,n)=Q( krow ,n)+(-4*G*pi*rho*a);
        Q( krow ,n+2*N)=Q( krow ,n+2*N)+n;
    end
%     krow=krow+1;
%     Q( krow ,1+2*N)=1;
end
%%
x=Q\b;
%%
Ua=sum(x(1:N));
figure
bar(x(1:N));
Vnorthpole = G * moon / rp^3 * a^2;
hhydro = Vnorthpole/gplanet;
h2 = Ua / hhydro *fac 
%%
figure
Ur=0;
subplot(3,1,1)
for i=1:N
    Ur=Ur+x(i)*r.^i;
end
plot(r,Ur)
subplot(3,1,3)
Kr=0;
for i=1:N
    Kr=Kr+x(i+2*N)*r.^i;
end
plot(r,Kr)
subplot(3,1,2)
Vr=0;
for i=1:N
    Vr=Vr+x(i+N)*r.^i;
end
plot(r,Vr)
%%