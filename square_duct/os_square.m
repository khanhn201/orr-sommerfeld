% Flow in z-dir
% Magnetic field in x-dir

clear all; close all;
output_precision(9);
set(0, "defaultAxesFontSize", 24)
set(0, "defaultTextFontSize", 24)
set(0, "defaultLineLineWidth", 2)
N = 24;


Re = 2e3;
Ha = 1e1;

rho = 1;
mu = 1/Re;
% sigma = 0.0277924;
% sigma_w = 0.0455531;
W = 0.1315;

% By = Ha/sqrt(sigma/mu)/L
By = 1;
sigma = Ha^2/Re;


alpha = 1.0; % Wave number in z
alpha = 0.8; % Wave number in z

[Ah,Bh,Ch,Dh,z,w] = semhat(N);
[U,Phi,f] = solve_steady(N, mu, sigma, By);
[uvec,vvec,wvec,pvec,phivec,gamma] = solve_lin(N, rho, mu, sigma, U,Phi, By, alpha);
c = gamma*1i/alpha

% [U,Phi] = solve_steady_wall(N, mu, sigma, sigma_w, W, By, f);
% [uvec,vvec,wvec,pvec,phivec] = solve_lin_wall(N, rho, mu, sigma,sigma_w, W, U,Phi, By, f, alpha);

[zd,wd]=zwgl(N-1);


% Plot steady
Phi = reshape(Phi,N+1,N+1);
[y,x] = meshgrid(z,z);
figure;
surf(x,y,abs(Phi));
title('phi');

U = reshape(U,N+1,N+1);
[y,x] = meshgrid(z,z);
figure;
surf(x,y,abs(U));
title('u');
% plot_phi_3b3(N, Phi, W, z);



% Plot perturb mode
U = reshape(uvec,N+1,N+1);
[y,x] = meshgrid(z,z);
figure;
surf(x,y,abs(U));
title('u');
U = reshape(vvec,N+1,N+1);
[y,x] = meshgrid(z,z);
figure;
surf(x,y,abs(U));
title('v');
% P = reshape(pvec,N-1,N-1);
% [y,x] = meshgrid(zd,zd);
% figure;
% surf(x,y,abs(P));
% shading interp;
% title('p');

% plot_phi_3b3(N, phivec, W, z);
Phi = reshape(phivec,N+1,N+1);
[y,x] = meshgrid(z,z);
figure;
surf(x,y,abs(Phi));
title('phi');
