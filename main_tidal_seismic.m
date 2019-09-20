
%% this is our first github collaborative project 

clear;
close all;
clc;
set(0,'defaultaxesfontsize',14);
set(0,'defaulttextfontsize',14);

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

%% define number of radial points 
% N = 100; 
N = 300; 
dr = a/(N-1); % discretization interval 

% init [U, V, Kr] along the radial direction 
U(1:N) = 0; 
V(1:N) = 0; 
Kr(1:N) =0 ; 
% stiffness matrix Q(3*N, 3*N) 
Q = zeros (3*N); 
b = zeros(3*N,1) ;

%% assemble matrices 
% (N-2) equations for each ODE 

% equation 11 (yuan tian job) 
krow = 0; 
ku = 0; 
kv = N; 
kk = 2*N; 
for j =2: N-1 
    rj = (j-1)*dr; 
    krow = krow + 1; % row counter number 
    
    b(krow) = - 2 * A * rj^3 * rho; 
    Q(krow, ku+j) = Q(krow, ku+j) + 4/3 * G * pi * rj^2 * rho^2  - 10 * mu - 2 * B * mu;
    Q(krow, kv+j) = Q(krow, kv+j) - 8 * G * pi * rj^2 * rho^2 + 18 * mu + 6 *B * mu; 
    Q(krow, kk+j-1:kk+j) = Q(krow, kk+j-1:kk+j)+ rj^2 * rho *[-1 1]/dr;    
    Q(krow, ku+j-1:ku+j) =  Q(krow, ku+j-1:ku+j) + ( 4 * rj * mu + 2 * B * rj * mu ) * [-1 1]/dr ;
    Q(krow, kv+j-1:kv+j) =  Q(krow, kv+j-1:kv+j) + ( -6 * rj * rho - 6 * B * rj * mu ) * [-1 1]/dr ;
    Q(krow, ku+j-1:ku+j+1) = Q (krow, ku+j-1:ku+j+1)+ ( 2*rj^2*mu + B * rj^2 * mu ) *  [1 -2 1]/dr^2;
end 

% equation 22 (Hao Hu job); 
krow =(N-2); 
ku = 0; 
kv = N; 
kk = 2*N; 
for j =2: N-1 
    rj = (j-1)*dr; 
    krow = krow + 1; % row counter number 
    
    b(krow) = - A * rj^3 * rho; 
    Q(krow, kk+j) = Q(krow, kk+j) + rj * rho;
    Q(krow, ku+j) = Q(krow, ku+j) - 4/3 * G * pi * rj^2 * rho^2  + 4 * mu + 2 * B * mu; 
    Q(krow, kv+j) = Q(krow, kv+j) - 12 * mu - 6 * B * mu ; 
    Q(krow, ku+j-1:ku+j) =  Q(krow, ku+j-1:ku+j) + ( rj * mu + B * rj * mu ) * [-1 1]/dr ;
    Q(krow, kv+j-1:kv+j) =  Q(krow, kv+j-1:kv+j) + ( 2 * rj * mu ) * [-1 1]/dr ;
    Q(krow, kv+j-1:kv+j+1) = Q(krow, kv+j-1:kv+j+1)+ rj^2*mu*[1 -2 1]/dr^2;    
end 

% equation 33 (Zheng job)  
krow =2*(N-2); 
ku=0; 
kv = N; 
kk = 2*N; 
for j =2: N-1 
    rj = (j-1)*dr; 
    krow = krow + 1; % row number 
    Q (krow, kk+j) =Q (krow, kk+j) -6; 
    Q (krow, ku+j) = Q (krow, ku+j) -8*G*pi*rj*rho; 
    Q (krow, kv+j) = Q (krow, kv+j)+ 24*G*pi*rj*rho; 
    Q (krow, kk+j-1:kk+j) =Q (krow, kk+j-1:kk+j)+ 2*rj*[-1 1]/dr; 
    Q (krow, ku+j-1:ku+j) =  Q (krow, ku+j-1:ku+j)-4*G*pi*rj^2*rho*[-1 1]/dr; 
    Q (krow, ku+j-1:ku+j+1) = Q (krow, ku+j-1:ku+j+1)+ rj^2*[1 -2 1]/dr^2; 
end 


%% Boundary condtions 
% u(0) = 0, v(0) = 0
krow = krow + 1; 
Q(krow,1) = 1; 
krow = krow + 1; 
Q(krow,N+1) = 1; 

% dk/dr(r=0) = 0
krow = krow + 1; 
% Q(krow,kk+1:kk+2) = Q(krow,kk+1:kk+2) + [-1 1]/dr; 
Q(krow,kk+1) = Q(krow,kk+1) + 1;  % k(0)=0

% dk/dr(r=a) - 4Gpi rho U + 3 K/a 
krow = krow + 1; 
Q(krow,3*N) = Q(krow,3*N) + 3/a; % k 
Q(krow,N) = Q(krow,N) - 4 * G * pi * rho ; 
Q(krow,3*N-1:3*N) = Q(krow,3*N-1:3*N) + [-1 1]/dr; % k 

% trr = 0, 
krow = krow + 1; 
Q(krow,N) = Q(krow,N) + 2; 
Q(krow,2*N) = Q(krow,2*N) - 6; 
Q(krow,N-1:N) = Q(krow,N-1:N) + 3*a*[-1 1]/dr; 

%trtheta = 0;
krow = krow + 1; 
Q(krow,N) = Q(krow,N) + 1; 
Q(krow,2*N) = Q(krow,2*N) - 1; 
Q(krow,2*N-1:2*N) = Q(krow,2*N-1:2*N) + a*[-1 1]/dr; 


%% solve Q*X = b 
x = inv(Q) * b;
U = x(1:N); 
V = x(N+1:2*N); 
Kr = x(2*N+1:3*N); 

rr = dr * N; 
figure;
subplot(3,1,1)
plot(U);
subplot(3,1,2)
plot(V);
subplot(3,1,3)
plot(Kr);


% calculate the h2 
Vnorthpole = G * moon / rp^3 * a^2
hhydro = Vnorthpole/gplanet
h2 = U(end) / hhydro * fac
















 