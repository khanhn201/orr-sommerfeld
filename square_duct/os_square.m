% Flow in z-dir
% Magnetic field in x-dir

% clear all; close all;
output_precision(9);
set(0, "defaultAxesFontSize", 24)
set(0, "defaultTextFontSize", 24)
set(0, "defaultLineLineWidth", 2)

N = 20;

Re = 3e4;
Ha = 100;
c = 1.0;
W = 0.1;


% Re = 1e4;
% Ha = 15;
% c = 1.0;
% W = 0.1;

alpha = 1.0; % Wave number in z


rho = 1;
mu = 1/Re;
sigma = Ha^2/Re
sigma_w = c*sigma/W
By = 1;

% meshVel  = [-1.0, -0.95, -0.9, -0.8, -0.6, -0.4, -0.2, 0.0, ...
%              0.2, 0.4, 0.6, 0.8, 0.9, 0.95, 1.0];
% meshWall = [0.0, 0.25, 0.5, 1.0]*W + 1.0;

% meshVel  = [-1.0, -0.8, -0.4, 0.0,...
%              0.4, 0.8, 1.0];
meshVel  = [-1.0, -0.5, 0.0, 0.5, 1.0];
% meshVel  = [-1.0, 1.0];
meshWall = [0.0, 1.0]*W + 1.0;


% [U_base,Phi_base,f] = solve_steady(N, mu, sigma, By);
% [uvec,vvec,wvec,pvec,phivec,gamma] = solve_lin(N, rho, mu, sigma, U_base,Phi_base, By, alpha);

% [U_base1,Phi_base,f] = solve_steady_wall(N, mu, sigma, sigma_w, W, By);
% [uvec1,vvec,wvec,pvec,phivec,gamma] = solve_lin_wall(N, rho, mu, sigma,sigma_w, W, U_base1, Phi_base, By, alpha);

[U_base,Phi_base,f] = solve_steady_wall_multi(N, meshVel, meshWall, mu, sigma, sigma_w, W, By);
% [uvec,vvec,wvec,pvec,phivec,gamma] = solve_lin_wall_multi(N, meshVel, meshWall, rho, mu, sigma,sigma_w, W, U_base, Phi_base, By, alpha);

% f
% gamma
plot_vmesh_multi(N, meshVel, meshWall, U_base)

save("results_wall.mat")


