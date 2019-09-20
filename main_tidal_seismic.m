
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
N = 100; 
dr = a/(N-1); % discretization interval 

% init [U, V, Kr] along the radial direction 
U(1:N) = 0; 
V(1:N) = 0; 
Kr(1:N) =0 ; 
% stiffness matrix Q(3*N, 3*N) 
Q = zeros (3*N); 

<<<<<<< HEAD
<<<<<<< HEAD

%% assemble matrices 
% (N-2) equations for each ODE 

% equation 11 (yuan tian job) 

% equation 22 (Hao Hu job); 


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


 