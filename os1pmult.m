clear all; close all;
N = 10;
Nelem = 10;

Re = 3.e4;
alpha = 1.0;

Ca = 0.07;
Ga = 8.3e7;

[Ah,Bh,Ch,Dh,z,w] = semhat(N);

nh = N + 1;
Ng = Nelem * N + 1;  % total global nodes
Q = sparse([], [], [], Nelem*nh, Ng, 2*Ng);
for e = 1:Nelem
    gidx = (e-1)*N + (1:nh);
    Q((e-1)*nh + (1:nh), gidx) = speye(nh);
end

Lx2=1.0/Nelem/2.0; Lxi=1./Lx2;

Dh = Lxi*Dh;
Bh = Lx2*Bh;

nh = N+1;

xs = [];


for e = 1:Nelem
    x = -1.0 + e/Nelem + (z-1.0)/2.0/Nelem;
    xs = [xs, x'];
    U = 1.0 - x.^2;
    DU = Dh * U;
    D2U = Dh * DU;
    if e == Nelem
        S = 0*(1:nh); S(end) = 1;
        Sh = diag(S);
        S = zeros(nh,1); S(end,1) = 1;
    else
        S = 0*(1:nh);
        Sh = diag(S);
        S = zeros(nh,1);
    end

    Kuu0 = -Dh'*Dh'*Bh*Dh*Dh...
           -2*alpha^2*Dh'*Bh*Dh...
           -alpha^4*Bh;
    KuuU = -1i*alpha*Re*(alpha^2*Bh*diag(U)...
                    +Dh'*Bh*diag(U)*Dh...
                    -Dh'*Bh*diag(DU));
    KuuS = -alpha^2*(Dh'*Sh+Sh*Dh);
    Kuu = Kuu0 + KuuU + KuuS;
    Muu = Re*(Dh'*Bh*Dh + alpha^2*Bh);

    Kua = -alpha^2*S*(Ga/Re + alpha^2/Ca - 2i*alpha*S'*DU)...
          +1i*alpha*S'*D2U*Dh'*S;
    Kau = S';
    Kaa = -1i*alpha*S'*U;

    Kuu_block((e-1)*(N+1)+1:e*(N+1), (e-1)*(N+1)+1:e*(N+1)) = Kuu;
    Muu_block((e-1)*(N+1)+1:e*(N+1), (e-1)*(N+1)+1:e*(N+1)) = Muu;
    Kau_block(1,(e-1)*(N+1)+1:e*(N+1)) =  Kau;
    Kua_block((e-1)*(N+1)+1:e*(N+1),1) =  Kua;
end


Ih = speye(Ng); 
R=Ih(2:end,:);
Kuu_global = R*Q'*Kuu_block*Q*R';
Kua_global = R*Q'*Kua_block;
Kau_global = Kau_block*Q*R';
Muu_global = R*Q'*Muu_block*Q*R';

Maa = 1;

K = [
  Kuu_global, Kua_global;
  Kau_global, Kaa;
];
M = [
  Muu_global, R*Q'*zeros(Nelem*nh, 1);
  zeros(1, Nelem*nh)*Q*R', Maa;
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

figure;
unstable = find(imag(c) > -0.01);
plot(abs(R'*vecs(1:end-1, unstable)), Q'*xs', 'linewidth', 2)
labels = arrayfun(@(g) sprintf('gamma = %.2f + %.2fi, a =  %.2f + %.2fi', 
                                real(gamma(g)), 
                                imag(gamma(g)), 
                                real(vecs(end, g)),
                                imag(vecs(end, g)) ), unstable, 'UniformOutput', false);
xlabel('mag(V(y))');
ylabel('y');
legend(labels);
gamma(unstable)

N
res = norm(K*vecs-M*vecs*diag(gamma))
condK = cond(K)
condM = cond(M)
