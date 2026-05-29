% Flow in z-dir
% Magnetic field in x-dir

clear all; close all;
output_precision(9);
set(0, "defaultAxesFontSize", 24)
set(0, "defaultTextFontSize", 24)
set(0, "defaultLineLineWidth", 2)
N = 20;


Re = 1e5;
Ha = 10;

rho = 1;
mu = 1/Re;
sigma = 1.0;
sigma_w = 0.53308;
L = 1.0;
W = 0.1315;

By = 1.0;
% By = Ha/sqrt(sigma/mu)/L;
f = 2.0/Re; % Acceleration or pressure gradient
f = 1.0;


alpha = 1.0; % Wave number in z

[Ah,Bh,Ch,Dh,z,w] = semhat(N);
% [uvec,phivec] = solve_steady(N, mu, sigma, By, f);
[U,Phi] = solve_steady_wall(N, mu, sigma, sigma_w, W, By, f);
[uvec,vvec,wvec,pvec,phivec] = solve_lin(N, rho, mu, sigma,sigma_w, W, U,Phi, By, f, alpha);

[zd,wd]=zwgl(N-1);


% Plot steady
U = reshape(U,N+1,N+1);
[y,x] = meshgrid(z,z);
figure;
surf(x,y,abs(U));
shading interp;
title('u');

Nx = 3;
Ny = 3;
n  = N + 1;
n2 = n^2;

x_edges = [-1-W, -1, 1, 1+W];
y_edges = x_edges;

figure; hold on;

for e = 1:9
    idx = (e-1)*n2 + (1:n2);
    Phie = reshape(Phi(idx), N+1, N+1);

    ex = mod(e-1,3);
    ey = floor((e-1)/3);

    xa = x_edges(ex+1);
    xb = x_edges(ex+2);
    ya = y_edges(ey+1);
    yb = y_edges(ey+2);

    [Zy, Zx] = meshgrid(z,z);
    X = (xb-xa)/2 * (Zx + 1) + xa;
    Y = (yb-ya)/2 * (Zy + 1) + ya;

    surf(X, Y, Phie, 'EdgeColor','none');
end

shading interp;
title('\phi');
colorbar;
view(3);




% Plot perturb mode
U = reshape(uvec,N+1,N+1);
[y,x] = meshgrid(z,z);
figure;
surf(x,y,abs(U));
shading interp;
title('u');
U = reshape(vvec,N+1,N+1);
[y,x] = meshgrid(z,z);
figure;
surf(x,y,abs(U));
shading interp;
title('v');

P = reshape(pvec,N-1,N-1);
[y,x] = meshgrid(zd,zd);
figure;
surf(x,y,abs(P));
shading interp;
title('p');

Nx = 3;
Ny = 3;
n  = N + 1;
n2 = n^2;

x_edges = [-1-W, -1, 1, 1+W];
y_edges = x_edges;

figure; hold on;

for e = 1:9
    idx = (e-1)*n2 + (1:n2);
    Phie = reshape(phivec(idx), N+1, N+1);

    ex = mod(e-1,3);
    ey = floor((e-1)/3);

    xa = x_edges(ex+1);
    xb = x_edges(ex+2);
    ya = y_edges(ey+1);
    yb = y_edges(ey+2);

    [Zy, Zx] = meshgrid(z,z);
    X = (xb-xa)/2 * (Zx + 1) + xa;
    Y = (yb-ya)/2 * (Zy + 1) + ya;

    surf(X, Y, abs(Phie), 'EdgeColor','none');
end

shading interp;
title('\phi');
colorbar;
view(3);
