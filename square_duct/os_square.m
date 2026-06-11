% Flow in z-dir
% Magnetic field in x-dir

% clear all; close all;
% output_precision(9);
set(0, "defaultAxesFontSize", 24)
set(0, "defaultTextFontSize", 24)
set(0, "defaultLineLineWidth", 2)

N = 30;

Re = 1e4;
Ha = 15;
c = 1.0;
W = 0.1;

alpha = 1.0; % Wave number in z

rho = 1;
mu = 1/Re;
sigma = Ha^2/Re
sigma_w = c*sigma/W
By = 1;


% [U_base,Phi_base,f] = solve_steady(N, mu, sigma, By);
% [uvec,vvec,wvec,pvec,phivec,gamma] = solve_lin(N, rho, mu, sigma, U_base,Phi_base, By, alpha);

[U_base,Phi_base,f] = solve_steady_wall(N, mu, sigma, sigma_w, W, By);
f
[uvec,vvec,wvec,pvec,phivec,gamma] = solve_lin_wall(N, rho, mu, sigma,sigma_w, W, U_base, Phi_base, By, alpha);
gamma


save("results_wall.mat")
