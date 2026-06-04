% Flow in z-dir
% Magnetic field in x-dir

clear all; close all;
output_precision(9);
set(0, "defaultAxesFontSize", 24)
set(0, "defaultTextFontSize", 24)
set(0, "defaultLineLineWidth", 2)
N = 20;


Re = 6e3;
Ha = 15;
c = 1;
W = 0.1;

rho = 1;
mu = 1/Re;

% By = Ha/sqrt(sigma/mu)/L
By = 1;
sigma = Ha^2/Re;

sigma_w = c*sigma/W;


alpha = 1.0; % Wave number in z

[Ah,Bh,Ch,Dh,z,w] = semhat(N);
% [U_base,Phi_base,f] = solve_steady(N, mu, sigma, By);
% [uvec,vvec,wvec,pvec,phivec,gamma] = solve_lin(N, rho, mu, sigma, U_base,Phi_base, By, alpha);

[U_base,Phi_base,f] = solve_steady_wall(N, mu, sigma, sigma_w, W, By);
[uvec,vvec,wvec,pvec,phivec,gamma] = solve_lin_wall(N, rho, mu, sigma,sigma_w, W, U_base, Phi_base, By, alpha);
f
c = gamma*1i/alpha


save("results.mat")
