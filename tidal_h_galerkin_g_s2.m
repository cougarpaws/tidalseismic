clearvars;
set(0,'defaultaxesfontsize',14);
set(0,'defaulttextfontsize',14);
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
gra_flag=1;
w=0;
%% devide to M segment
M=2;
N=20;
Nr=5000/M;
dr=1/M;
bdw=sqrt(mu);
%% coupling coefficients of base fnction
c=zeros(N+2,N+2,M);
for m=1:M
    rl=(m-1)*dr;
    rh=m*dr;
    r=linspace(rl,rh,Nr);
    for k=1:N+2
        for n=1:N+2
            %c(k,n,m)=sum(r.^k.*r.^n)*dr;
            %c(k,n,m)=1/(n-1+k)*(rh^(k+n)-rl^(k+n-1));
            if k==1 && n==1
                c(k,n,m)=rh-rl;
            else
                k0=k-1;
                n0=n-1;
                c(k,n,m)=1/(k0+n0)*(rh^(k0+n0+1)-rl^(k0+n0));
            end
        end
        
    end
end
%% construct matrix
Q=zeros((3*N+6)*M+3,3*N*M);
b=zeros((3*N+6)*M+3,1);
krow=0;
for m=1:M
    %eq1
    for k=1:N
        krow=krow+1;
        for n=1:N
            n0=n-1;
            Q(krow,n+(m-1)*N*3)=Q(krow,n+(m-1)*N*3)+mu*(-10-2*B+4*n0+2*B*n0...
                +2*(n0-1)*n0+B*n0*(n0-1))*c(k,n,m)+rho*w^2*c(k,n+2,m);
            Q(krow,n+N+(m-1)*N*3)=Q(krow,n+N+(m-1)*N*3)+mu*(18+6*B-6*n0...
                -6*B*n0)*c(k,n,m);
            if gra_flag==1
                Q(krow,n+(m-1)*N*3)=Q(krow,n+(m-1)*N*3)+4/3*G*pi*a^2*rho^2*c(k,n+2,m);
                Q(krow,n+N+(m-1)*N*3)=Q(krow,n+N+(m-1)*N*3)-8*G*pi*a^2*rho^2*c(k,n+2,m);
                Q(krow,n+2*N+(m-1)*N*3)=Q(krow,n+2*N+(m-1)*N*3)+n0*a*rho*c(k,n+1,m);
            end
        end
        %b(krow)=-sum(rho*2*A*r.^3.*r.^k)*dr*a^3;
        b(krow)=-rho*A*2*a^3*c(k,4,m);
    end
    %eq2
    for k=1:N
        krow=krow+1;
        for n=1:N
            n0=n-1;
            Q(krow,n+(m-1)*N*3)=Q(krow,n+(m-1)*N*3)+mu*(4+2*B+n0+B*n0)*c(k,n,m);
            Q(krow,n+N+(m-1)*N*3)=Q(krow,n+N+(m-1)*N*3)+mu*(-12-6*B+2*n0+...
                (n0-1)*n0)*c(k,n,m)+rho*w^2*c(k,n+2,m);
            if gra_flag==1
                Q(krow,n+(m-1)*N*3)=Q(krow,n+(m-1)*N*3)-4/3*G*pi*a^2*rho^2*c(k,n+2,m);
                Q(krow,n+2*N+(m-1)*N*3)=Q(krow,n+2*N+(m-1)*N*3)+a*rho*c(k,n+1,m);
            end
        end
        b(krow)=-rho*A*a^3*c(k,4,m);
    end
    %eq3
    if gra_flag==1
        for k=1:N
            krow=krow+1;
            for n=1:N
                n0=n-1;
                Q(krow,n+2*N+(m-1)*N*3)=Q(krow,n+2*N+(m-1)*N*3)+(-6+2*n0+n0*(n0-1))*c(k,n,m);
                Q(krow,n+(m-1)*N*3)=Q(krow,n+(m-1)*N*3)+(-8*G*pi*a*rho-4*G*pi*a*n0*rho)*c(k,n+1,m);
                Q(krow,n+N+(m-1)*N*3)=Q(krow,n+N+(m-1)*N*3)+24*G*pi*a*rho*c(k,n+1,m);
            end
            
        end
    end
    
    %BC
    
    if m<M 
        % BC for U
        krow=krow+1;
        for n=1:N
            n0=n-1;
            Q(krow,n+(m-1)*N*3)=Q(krow,n+(m-1)*N*3)+(m/M)^n0*bdw;
            Q(krow,n+(m)*N*3)=Q(krow,n+(m)*N*3)+(-(m/M)^n0)*bdw;
        end
        krow=krow+1;
        for n=1:N
            n0=n-1;
            Q(krow,n+(m-1)*N*3)=Q(krow,n+(m-1)*N*3)+n0*(m/M)^(n0-1)*bdw;
            Q(krow,n+(m)*N*3)=Q(krow,n+(m)*N*3)+(-n0*(m/M)^(n0-1))*bdw;
        end
        krow=krow+1;
        % BC for V
        for n=1:N
            n0=n-1;
            Q(krow,n+N+(m-1)*N*3)=Q(krow,n+N+(m-1)*N*3)+(m/M)^n0*bdw;
            Q(krow,n+N+(m)*N*3)=Q(krow,n+N+(m)*N*3)+(-(m/M)^n0)*bdw;
        end
        krow=krow+1;
        for n=1:N
            n0=n-1;
            Q(krow,n+N+(m-1)*N*3)=Q(krow,n+N+(m-1)*N*3)+n0*(m/M)^(n0-1)*bdw;
            Q(krow,n+N+(m)*N*3)=Q(krow,n+N+(m)*N*3)+(-n0*(m/M)^(n0-1))*bdw;
        end
        krow=krow+1;
        % BC for K
        if gra_flag==1
            for n=1:N
                n0=n-1;
                Q(krow,n+2*N+(m-1)*N*3)=Q(krow,n+2*N+(m-1)*N*3)+(m/M)^n0*bdw;
                Q(krow,n+2*N+(m)*N*3)=Q(krow,n+2*N+(m)*N*3)+(-(m/M)^n0)*bdw;
            end
        end
        krow=krow+1;
        if gra_flag==1
            for n=1:N
                n0=n-1;
                Q(krow,n+2*N+(m-1)*N*3)=Q(krow,n+2*N+(m-1)*N*3)+n0*(m/M)^(n0-1)*bdw;
                Q(krow,n+2*N+(m)*N*3)=Q(krow,n+2*N+(m)*N*3)+(-n0*(m/M)^(n0-1))*bdw;
            end
        end
    end
end
% BC1
krow=krow+1;
for n=1:N
    n0=n-1;
    Q(krow,n+(M-1)*N*3)=Q(krow,n+(M-1)*N*3)+(2*B+n0*(2+B))*bdw;
    Q(krow,n+N+(M-1)*N*3)=Q(krow,n+N+(M-1)*N*3)-6*B*bdw;
end
% BC2
krow=krow+1;
for n=1:N
    n0=n-1;
    Q(krow,n+(M-1)*N*3)=Q(krow,n+(M-1)*N*3)+1*bdw;
    Q(krow,n+N+(M-1)*N*3)=Q(krow,n+N+(M-1)*N*3)+(-1+n0)*bdw;
end
%BC3

if gra_flag==1
    krow=krow+1;
    for n=1:N
        n0=n-1;
        Q( krow ,n+2*N+(M-1)*N*3)=Q( krow ,n+2*N+(M-1)*N*3)+3*bdw;
        Q( krow ,n+(M-1)*N*3)=Q( krow ,n+(M-1)*N*3)+(-4*G*pi*rho*a)*bdw;
        Q( krow ,n+2*N+(M-1)*N*3)=Q( krow ,n+2*N+(M-1)*N*3)+n0*bdw;
    end
     krow=krow+1;
     Q( krow ,1+2*N)=1*bdw;
end
krow=krow+1;
Q( krow ,1)=1*bdw;
krow=krow+1;
Q( krow ,1+N)=1*bdw;
%%
x=Q\b;
%%
Ua=sum(x(1+(M-1)*N*3:N+(M-1)*N*3));
Vnorthpole = G * moon / rp^3 * a^2;
hhydro = Vnorthpole/gplanet;
h2 = Ua / hhydro *fac 
%%
figure

subplot(3,1,1)
Ur=[];
rall=[];
for m=1:M
    rl=(m-1)*dr;
    rh=m*dr;
    r=linspace(rl,rh,Nr);
    Ur0=0;
    for n=1:N
       Ur0=Ur0+x(n+(m-1)*N*3)*r.^n;
    end
    Ur=[Ur,Ur0];
    rall=[rall,r];
end
plot(rall,Ur)
title('U(r)')
subplot(3,1,2)
Vr=[];

for m=1:M
    rl=(m-1)*dr;
    rh=m*dr;
    r=linspace(rl,rh,Nr);
    Vr0=0;
    for n=1:N
       Vr0=Vr0+x(n+N+(m-1)*N*3)*r.^n;
    end
    Vr=[Vr,Vr0];
end
plot(rall,Vr)
title('V(r)')
subplot(3,1,3)
Kr=[];

for m=1:M
    rl=(m-1)*dr;
    rh=m*dr;
    r=linspace(rl,rh,Nr);
    Kr0=0;
    for n=1:N
       Kr0=Kr0+x(n+2*N+(m-1)*N*3)*r.^n;
    end
    Kr=[Kr,Kr0];
end
plot(rall,Kr)
title('K(r)')