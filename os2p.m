% Based on Giannakis (2009)
% https://doi.org/10.1016/j.jcp.2008.10.016
% but extends to two phases flow

clear all; close all;
output_precision(9);
set(0, "defaultAxesFontSize", 24)
set(0, "defaultTextFontSize", 24)
set(0, "defaultLineLineWidth", 2)
N = 250;


Re = 3e4;
Hz = 14.0;
ratio = 1e8;
Ca = 0.07;
Ga = 8.3e7;

rhos = [1, 1/ratio];
mus = [1/Re, 1/Re/ratio];
% sigmas = [1/(mu_mag/Re*Pm), 1/(mu_mag/Re*Pm)/ratio];
sigmas = [1/1, 1/1/ratio];
st = mus(1)/Ca; % surface tension
g = Ga*mus(1)^2/rhos(1)^2;
alpha = 1;

h = 1; % Height of the air domain
mode = 1; % Which of the unstable modes to output
top_bc = 'W';
accel_type = 'S';
Bz = Hz;
Bx = 0;

[xs,umat,vmat,amat,gamma,f,r] = solve_os2p(alpha, N, rhos, mus, sigmas, st, g, h, top_bc, accel_type, Bx, Bz);
c = gamma*1i/alpha;
lambda = alpha*c;


% figure;
% plot(real(c), imag(c), 'xr', "markersize", 15); hold on;
% A = load("giannakis_rosner.dat");
% x = A(:,1);
% y = A(:,2);
% plot(x, y, "ok", "markersize", 15)
% legend("Density ratio = 10", "Giannakis 2009")
% xlabel('Im(gamma)');
% ylabel('Re(gamma)');
% xlim([-0.2, 1.5]);
% ylim([-0.00, 0.01]);
% grid on;
%
%
% asdfzcxvxcv

N

% unstable = find(imag(c) > 0.0);
unstable = find(real(gamma) > -1e-1);
c_unstable = c(unstable)
unstable = unstable(mode);
a = amat(unstable)
v = vmat(:, unstable);
u = umat(:, unstable);
v = reshape(v, [N+1, 2]);
u = reshape(u, [N+1, 2]);
% figure;
% plot(abs(v), xs', 'linewidth', 2)
% title('V')
% figure;
% plot(abs(u), xs', 'linewidth', 2);
% title('U')


% figure; hold on;
% for x=xs(1,:);
%   plot([0, 1], [x,x], '-k');
% end
% for x=xs(2,:);
%   plot([0, 1], [x,x], '-k');
% end
% plot([0, 1], [0,0], '-r', 'linewidth', 4);
% plot([0, 1], [1,1], '-b', 'linewidth', 4);
% plot([0, 1], [-1,-1], '-b', 'linewidth', 4);

% Interpolate to Nek mesh
nely = 20; % Nely in 1 phase; must be divisible by 2 for geometric
nelx = -30;
Nf = 8; % lx1
el_ratio = 16.0000;
zf = zeros(nely*Nf,1);
[zff,wf] = zwgll(Nf-1);
el_pos = linspace(-1, 1, nely + 1); % Linear element

geom_fac = (el_ratio)^(1/(nely/2));
el_pos_tmp = geom_fac.^(0:nely/2) - 1.0;
el_pos_tmp = el_pos_tmp/el_pos_tmp(end) - 1.0;
el_pos_tmp2 = flip(-el_pos_tmp);
el_pos(1:nely/2+1) = el_pos_tmp;
el_pos(nely/2+1:end) = el_pos_tmp2;

delta_el = diff(el_pos);
min_delta = min(delta_el)/2
max_delta = max(delta_el)/2
el_center = (el_pos(2:end) + el_pos(1:end-1))/2.0;
for i=1:nely
    zf((i-1)*Nf+1:i*Nf) = el_center(i) + zff*delta_el(i)/2.0;
end

[Ah,Bhh,Ch,Dhh,z,w] = semhat(N);
J = interp_mat(zf,z);
vf = J*v; vf = vf(:);
uf = J*u; uf = uf(:);
data = [real(uf), imag(uf), real(vf), imag(vf)];

% u2.txt
fid = fopen("u2.txt", "w");
fprintf(fid, "%d  %d  ! nelx, nely\n", -nelx, nely*2);
fprintf(fid, "% .16e  % .16e  % .16e  % .16e  ! rhol, rhog, mul, mug\n", rhos(1), rhos(2), mus(1), mus(2));
fprintf(fid, "% .16e  % .16e  % .16e  % .16e  ! f1, f2, g, st\n", f, f*r, g, st);
fprintf(fid, "% .16e  ! alpha\n", alpha);
fprintf(fid, "% .16e  % .16e  ! gamma\n", real(gamma(unstable)), imag(gamma(unstable)));
fprintf(fid, "% .16e  % .16e  ! a\n", real(a), imag(a));
fprintf(fid, "% .16e  % .16e  ! U1 = 1.0 + A*y + B*y^2\n", 1-f/2*rhos(1)/mus(1), -f/2*rhos(1)/mus(1)); % U profile in water
fprintf(fid, "% .16e  % .16e  ! U2 = 1.0 + A*y + B*y^2\n", mus(1)/mus(2)-f/2*rhos(1)/mus(2),  -f*r/2*rhos(2)/mus(2)); % U profile in air
fprintf(fid, "% .16e  % .16e  % .16e  % .16e\n", data.');
fclose(fid);


% Box file
el_pos_real = 1:(nely*2 + 1);
el_pos_real(1:nely+1) = (el_pos - 1.0)/2.0;
el_pos_real(nely+1:end) = (el_pos + 1.0)/2.0;
fid = fopen("fs_g.box", "w");
fprintf(fid, "-2\n");
fprintf(fid, "5\n");
fprintf(fid, "Box\n");
fprintf(fid, "%d %d\n", nelx, nely*2);
fprintf(fid, "0 6.0 1.\n");
fprintf(fid, "%f\n", el_pos_real);
fprintf(fid, "P  ,P  ,W  ,SYM\n");
fprintf(fid, "P  ,P  ,I  ,I  \n");
fprintf(fid, "P  ,P  ,I  ,I  \n");
fprintf(fid, "P  ,P  ,I  ,I  \n");
fprintf(fid, "P  ,P  ,I  ,I  \n");
fclose(fid);
