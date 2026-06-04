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
% sigma = 0.0277924;
% sigma_w = 0.0455531;
W = 0.1315;

% By = Ha/sqrt(sigma/mu)/L
By = 1;
sigma = Ha^2/Re;


alpha = 1.0; % Wave number in z
alpha = 0.8; % Wave number in z

[Ah,Bh,Ch,Dh,z,w] = semhat(N);
[U_base,Phi_base,f] = solve_steady(N, mu, sigma, By);
f
[uvec,vvec,wvec,pvec,phivec,gamma] = solve_lin(N, rho, mu, sigma, U_base,Phi_base, By, alpha);
c = gamma*1i/alpha

% [U,Phi] = solve_steady_wall(N, mu, sigma, sigma_w, W, By, f);
% [uvec,vvec,wvec,pvec,phivec] = solve_lin_wall(N, rho, mu, sigma,sigma_w, W, U,Phi, By, f, alpha);

[zd,wd]=zwgl(N-1);


% Plot steady
% U = reshape(U_base,N+1,N+1);
% [y,x] = meshgrid(z,z);
% figure;
% surf(x,y,abs(U));
% title('U');
% Phi = reshape(Phi_base,N+1,N+1);
% [y,x] = meshgrid(z,z);
% figure;
% surf(x,y,abs(Phi));
% title('Phi');

% plot_phi_3b3(N, Phi, W, z);



% Plot perturb mode
% U = reshape(uvec,N+1,N+1);
% [y,x] = meshgrid(z,z);
% figure;
% surf(x,y,abs(U));
% title('u');
% U = reshape(vvec,N+1,N+1);
% [y,x] = meshgrid(z,z);
% figure;
% surf(x,y,abs(U));
% title('v');
% U = reshape(wvec,N+1,N+1);
% [y,x] = meshgrid(z,z);
% figure;
% surf(x,y,abs(U));
% title('w');
% P = reshape(pvec,N-1,N-1);
% [y,x] = meshgrid(zd,zd);
% figure;
% surf(x,y,abs(P));
% title('p');
% Phi = reshape(phivec,N+1,N+1);
% [y,x] = meshgrid(z,z);
% figure;
% surf(x,y,abs(Phi));
% title('phi');

% plot_phi_3b3(N, phivec, W, z);


% Interp and export to Nek mesh
Nf = 6; % lx1
el_pos = [-1.0, -0.95, -0.9, -0.5, -0.25, 0.0, 0.25, 0.5, 0.9, 0.95, 1.0];

nely = length(el_pos) - 1;
[zff,wf] = zwgll(Nf-1);
zf = zeros(nely*Nf,1);
delta_el = diff(el_pos);
el_center = (el_pos(2:end) + el_pos(1:end-1))/2.0;
for i=1:nely
    zf((i-1)*Nf+1:i*Nf) = el_center(i) + zff*delta_el(i)/2.0;
end

[Ah,Bhh,Ch,Dhh,z,w] = semhat(N);
J = interp_mat(zf,z);
J2D = kron(J,J);

uf   = J2D*uvec;
vf   = J2D*vvec;
wf   = J2D*wvec;
phif = J2D*phivec;
data = [real(uf), imag(uf), real(vf), imag(vf), real(wf), imag(wf), real(phif), imag(phif)];


fid = fopen("u2.txt", "w");
fprintf(fid, "%d  ! nelxy\n", nely);
fprintf(fid, "% .16e  % .16e  % .16e  % .16e  ! rho, mu, sigma\n", rho, mu, sigma);
fprintf(fid, "% .16e  % .16e  ! f, By\n", f, By);
fprintf(fid, "% .16e  ! alpha\n", alpha);
fprintf(fid, "% .16e  % .16e  ! gamma\n", real(gamma), imag(gamma));
fprintf(fid, "% .16e  % .16e  % .16e  % .16e  % .16e  % .16e  % .16e  % .16e\n", data.');
fclose(fid);


wBasef = J2D*U_base;
phiBasef = J2D*Phi_base;
data = [wBasef, phiBasef];
fid = fopen("ubase.txt", "w");
fprintf(fid, "%d  ! nelxy\n", nely);
fprintf(fid, "% .16e  % .16e  % .16e  % .16e  ! rho, mu, sigma\n", rho, mu, sigma);
fprintf(fid, "% .16e  % .16e  ! f, By\n", f, By);
fprintf(fid, "% .16e  ! alpha\n", alpha);
fprintf(fid, "% .16e  % .16e  ! gamma\n", real(gamma), imag(gamma));
fprintf(fid, "% .16e  % .16e\n", data.');
fclose(fid);

