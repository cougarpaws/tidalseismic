
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
vs = 1200; % m/s 
vp = 3000; % m/s 

mu = rho*vs*vs; % shear modulus 
lambda = rho*vp*vp - 2 *mu; % Lame constant  


