% Flow in z-dir
% Magnetic field in x-dir

clear all; close all;
output_precision(9);
set(0, "defaultAxesFontSize", 24)
set(0, "defaultTextFontSize", 24)
set(0, "defaultLineLineWidth", 2)
N = 20;


Re = 2.1e3;
Ha = 1e1;

rho = 1;
mu = 1/Re;

% By = Ha/sqrt(sigma/mu)/L
By = 1;
sigma = Ha^2/Re;
sigma_w = 0.0455531;

W = 0.1315;

alpha = 0.8; % Wave number in z

[Ah,Bh,Ch,Dh,z,w] = semhat(N);
% [U_base,Phi_base,f] = solve_steady(N, mu, sigma, By);
% [uvec,vvec,wvec,pvec,phivec,gamma] = solve_lin(N, rho, mu, sigma, U_base,Phi_base, By, alpha);

[U_base,Phi_base,f] = solve_steady_wall(N, mu, sigma, sigma_w, W, By);
[uvec,vvec,wvec,pvec,phivec] = solve_lin_wall(N, rho, mu, sigma,sigma_w, W, U_base, Phi_base, By, alpha);
f
c = gamma*1i/alpha


save("results.mat")
