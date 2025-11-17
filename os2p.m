clear all; close all;
N = 200;

rhos = [1, 1e-1];
mus = [5e-5, 5e-6];
sigmas = 0.06;
sigmas = mus(1)/0.07;
sigmas = 0.0;
g = 9.81;
g = 8.3e7*mus(1)^2/rhos(1)^2;

alpha = 1;

top_bc = 'W'

if top_bc == 'W'
    f = 2*(mus(1)+mus(2))/(rhos(1)+rhos(2)); % If top is W BC 
else
    f = mus(1)/(rhos(1)/2 + rhos(2)); % If top is SYM BC
end

[Ah,Bh,Ch,Dh,z,w] = semhat(N);
x = (z-1.0)/2.0;
Lx2=1.0/2.0; Lxi=1./Lx2;

Dh = Lxi*Dh;
Bh = Lx2*Bh;

Nelem = 2;
nh = N + 1;
Ng = Nelem * (N-1) + 2;  % total global nodes
Q = sparse([], [], [], Nelem*nh, Ng, 2*Ng);
for e = 1:Nelem
    gidx = (e-1)*(nh-2) + (1:nh);
    Q((e-1)*nh + (1:nh), gidx) = speye(nh);
end

T = speye(nh);
T(2, :) = Dh(1,:);
T(end-1, :) = speye(nh)(end, :);
T(end, :) = Dh(end,:);
T = inv(T);

xs = [];
Uxplot = [];
Uplot = [];
for e = 1:Nelem
    rho = rhos(e);
    mu = mus(e);

    x = -1.0 + 2.0*e/Nelem + (z-1.0)/Nelem;
    xs = [xs; x'];
    if e == 1
        U = 1.0 + (1-f/2*rho/mu)*x - (f/2*rho/mu)*x.^2;
    else
        U = 1.0 + (mus(1)/mu - f/2*rhos(1)/mu)*x - (f/2*rho/mu)*x.^2;
    end
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
end
figure;
plot(Uplot, Uxplot, 'linewidth', 2);


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

figure;
plot(real(c), imag(c), 'ok', 'linewidth', 2);
xlabel('Re(c)');
ylabel('Im(c)');
xlim([-0.2, 1.8]);
ylim([-1, 0.2]);
grid on;
axis equal;


unstable = find(imag(c) > 0.0);
figure;
a = vecs(end, unstable)
v = T_block*Q*R'*vecs(1:end-1, unstable);
v = reshape(v, [nh, 2]);
u = alpha*1i*Dh*v;
plot(abs(v), xs', 'linewidth', 2)
title('V')
figure;
plot(abs(u), xs', 'linewidth', 2);
title('U')

c(unstable)
N
res = norm(K*vecs-M*vecs*diag(gamma))
condK = cond(K)
condM = cond(M)


% Interpolate to Nek mesh
nely = 20; % Nely in 1 phase
Nf = 7; % lx1-1
zf = zeros(nely*(Nf+1),1);
[zff,wf] = zwgll(Nf);
zff = zff/nely;
for i=1:nely
    zf((i-1)*(Nf+1)+1:i*(Nf+1)) = ((i-1)/(nely-1)*2 - 1)*(1-1/nely) + zff;
end
J = interp_mat(zf,z);
vf = J*v; vf = vf(:);
uf = J*u; uf = uf(:);
data = [real(uf), imag(uf), real(vf), imag(vf)];
fid = fopen("u2.txt", "w");
fprintf(fid, "% .16e  % .16e  % .16e  % .16e  ! rhol, rhog, mul, mug\n", rhos(1), rhos(2), mus(1), mus(2));
fprintf(fid, "% .16e  % .16e  % .16e  ! f, g, sigma\n", f, g, sigmas);
fprintf(fid, "% .16e  ! alpha\n", alpha);
fprintf(fid, "% .16e  % .16e  ! gamma\n", real(gamma(unstable)), imag(gamma(unstable)));
fprintf(fid, "% .16e  % .16e  ! a\n", real(a), imag(a));
fprintf(fid, "% .16e  % .16e  ! U1 = 1.0 + A*y + B*y^2\n", 1-f/2*rhos(1)/mus(1), -f/2*rhos(1)/mus(1)); % U profile in water
fprintf(fid, "% .16e  % .16e  ! U2 = 1.0 + A*y + B*y^2\n", mus(1)/mus(2)-f/2*rhos(1)/mus(2),  -f/2*rhos(2)/mus(2)); % U profile in air
fprintf(fid, "% .16e  % .16e  % .16e  % .16e\n", data.');
fclose(fid);
