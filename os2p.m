% Based on Giannakis (2009)
% https://doi.org/10.1016/j.jcp.2008.10.016
% but extends to two phases flow

clear all; close all;
output_precision(9);
N = 200;

rhos = [1, 1e-3];
mus = [1/3e4, 1/3e4*1e-3];
sigmas = 0.06;
sigmas = mus(1)/0.07;
% sigmas = 0.0;
g = 9.81;
g = 8.3e7*mus(1)^2/rhos(1)^2;

alpha = 1;
h = 1; % Height of the air domain
mode = 1; % Which of the unstable modes to output

top_bc = 'W'
accel_type = 'S'

if top_bc == 'W'
    if accel_type == 'S' % same acceleration on both domains
        f = 2*(mus(2) + mus(1)*h)/(rhos(1)*h + rhos(2)*h^2);
        r = 1.0;
    elseif accel_type == 'F' % Prescribed accel on fluid; ratio is computed correspondingly
        f = mus(1)/rhos(1)*2;
        r = (2*(mus(2)+mus(1)*h)/f-rhos(1)*h)/(h^2*rhos(2));
    else % Prescribed f2/f1 ratio
        r = 1/h;
        f = 2*(mus(2)+mus(1)*h)/(r*(h^2*rhos(2)) + rhos(1)*h);
    end
else % Symmetry top bc
    if accel_type == 'S' % same acceleration on both domains
        f = mus(1)/(rhos(2)*h + rhos(1)/2)
        r = 1.0
    elseif accel_type == 'F' % Prescribed accel on fluid; ratio is computed correspondingly
        f = mus(1)/rhos(1)*2;
        r = 0.0;
    else % Prescribed f2/f1 ratio
        r = 1/h;
        f = mus(1)/mus(2)/(rhos(1)/mus(2)/2 + h*r*rhos(2)/mus(2));
    end
end



[Ah,Bhh,Ch,Dhh,z,w] = semhat(N);


Nelem = 2;
nh = N + 1;
Ng = Nelem * (N-1) + 2;  % total global nodes
Q = sparse([], [], [], Nelem*nh, Ng, 2*Ng);
for e = 1:Nelem
    gidx = (e-1)*(nh-2) + (1:nh);
    Q((e-1)*nh + (1:nh), gidx) = speye(nh);
end

xs = [];
Uxplot = [];
Uplot = [];
for e = 1:Nelem
    rho = rhos(e);
    mu = mus(e);

    if e == 1
        x = (z-1.0)/2; % x = [-1, 0]
        U = 1.0 + (1-f/2*rho/mu)*x - (f/2*rho/mu)*x.^2;
        Dh = 2*Dhh;
        Bh = 1/2*Bhh;
    else
        x = (z+1.0)/2*h; % x = [0, h]
        U = 1.0 + (mus(1)/mu - f/2*rhos(1)/mu)*x - (r*f/2*rho/mu)*x.^2;
        Dh = 2/h*Dhh;
        Bh = h/2*Bhh;
    end
    T = speye(nh);
    T(2, :) = Dh(1,:);
    T(end-1, :) = speye(nh)(end, :);
    T(end, :) = Dh(end,:);
    T = inv(T);

    xs = [xs; x'];
    Uxplot = [Uxplot;x];
    Uplot = [Uplot;U];
    DU = Dh * U;
    D2U = Dh * DU;

    if e == 1
        S = 0*(1:nh); S(end) = 1;
        Sh = diag(S);
        S = zeros(nh,1); S(end,1) = 1;
    else
        S = 0*(1:nh); S(1) = -1;
        Sh = diag(S);
        S = zeros(nh,1); S(1,1) = -1;
    end
    Kuu0 = -mu*Dh'*Dh'*Bh*Dh*Dh...
           -mu*2*alpha^2*Dh'*Bh*Dh...
           -mu*alpha^4*Bh;
    KuuU = -1i*alpha*rho*(alpha^2*Bh*diag(U)...
                    +Dh'*Bh*diag(U)*Dh...
                    -Dh'*Bh*diag(DU));
    KuuS = -alpha^2*mu*(Dh'*Sh + Sh*Dh);
    Kuu = Kuu0 + KuuU + KuuS;

    Kua = -alpha^2*S*rho*g...
          -0.5*abs(S)*alpha^4*sigmas... % 0.5 to counter multiplicity
          +1i*alpha*mu*Dh'*Sh*D2U;
    Kau = 0.5*abs(S');  % 0.5 to counter multiplicity

    Muu = rho*(Dh'*Bh*Dh + alpha^2*Bh);

    Kuu_block((e-1)*(N+1)+1:e*(N+1), (e-1)*(N+1)+1:e*(N+1)) = T'*Kuu*T;
    Muu_block((e-1)*(N+1)+1:e*(N+1), (e-1)*(N+1)+1:e*(N+1)) = T'*Muu*T;
    Kau_block(1,(e-1)*(N+1)+1:e*(N+1)) =  Kau*T;
    Kua_block((e-1)*(N+1)+1:e*(N+1),1) =  T'*Kua;
    T_block((e-1)*(N+1)+1:e*(N+1), (e-1)*(N+1)+1:e*(N+1)) = T;
    D_block((e-1)*(N+1)+1:e*(N+1), (e-1)*(N+1)+1:e*(N+1)) = Dh;
end
% figure;
% plot(Uplot, Uxplot, 'linewidth', 2);


Ih = speye(Ng); 
if top_bc == 'W'
    R=Ih(3:end-2,:);
else
    R=Ih(3:end-1,:);
    R(end, :) = Ih(end,:);
end


Kuu_global = R*Q'*Kuu_block*Q*R';
Kua_global = R*Q'*Kua_block;
Kau_global = Kau_block*Q*R';
Muu_global = R*Q'*Muu_block*Q*R';

K = [
  Kuu_global, Kua_global;
  Kau_global, -1i*alpha*1.0; % U(0) = 1.0
];
M = [
  Muu_global, R*Q'*zeros(Nelem*nh, 1);
  zeros(1, Nelem*nh)*Q*R', 1;
];

[vecs, gamma] = eig(K, M, 'vector');
c = gamma*1i/alpha;

% figure;
% plot(real(c), imag(c), 'ok', 'linewidth', 2);
% xlabel('Re(c)');
% ylabel('Im(c)');
% xlim([-0.2, 1.8]);
% ylim([-1, 0.2]);
% grid on;
% axis equal;

N
res = norm(K*vecs-M*vecs*diag(gamma))
condK = cond(K)
condM = cond(M)


unstable = find(imag(c) > 0.0);
c_unstable = c(unstable)
unstable = unstable(mode);
vecs(:, unstable) = vecs(:, unstable)/abs(vecs(end, unstable));
a = vecs(end, unstable)
v = T_block*Q*R'*vecs(1:end-1, unstable);
u = alpha*1i*D_block*v;
v = reshape(v, [nh, 2]);
u = reshape(u, [nh, 2]);
% figure;
% plot(abs(v), xs', 'linewidth', 2)
% title('V')
% figure;
% plot(abs(u), xs', 'linewidth', 2);
% title('U')


% Interpolate to Nek mesh
nely = 40; % Nely in 1 phase; must be divisible by 2 for geometric
nelx = -100;
Nf = 8; % lx1
el_ratio = 8.0000;
zf = zeros(nely*Nf,1);
[zff,wf] = zwgll(Nf-1);
el_pos = linspace(-1, 1, nely + 1); % Linear element
% geom_fac = (el_ratio)^(1/(nely/2));
%
% el_pos_tmp = geom_fac.^(0:nely/2) - 1.0;
% el_pos_tmp = el_pos_tmp/el_pos_tmp(end) - 1.0;
% el_pos_tmp2 = flip(-el_pos_tmp);
% el_pos(1:nely/2+1) = el_pos_tmp;
% el_pos(nely/2+1:end) = el_pos_tmp2;

delta_el = diff(el_pos);
min_delta = min(delta_el)/2
max_delta = max(delta_el)/2
el_center = (el_pos(2:end) + el_pos(1:end-1))/2.0;
for i=1:nely
    zf((i-1)*Nf+1:i*Nf) = el_center(i) + zff*delta_el(i)/2.0;
end

J = interp_mat(zf,z);
vf = J*v; vf = vf(:);
uf = J*u; uf = uf(:);
data = [real(uf), imag(uf), real(vf), imag(vf)];

% u2.txt
fid = fopen("u2.txt", "w");
fprintf(fid, "%d  %d  ! nelx, nely\n", -nelx, nely*2);
fprintf(fid, "% .16e  % .16e  % .16e  % .16e  ! rhol, rhog, mul, mug\n", rhos(1), rhos(2), mus(1), mus(2));
fprintf(fid, "% .16e  % .16e  % .16e  % .16e  ! f1, f2, g, sigma\n", f, f*r, g, sigmas);
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
