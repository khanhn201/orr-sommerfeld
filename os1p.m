% clear all; close all;
N = 200;

Re = 3.e4;
alpha = 1.0;
Ca = 0.07;
Ga = 8.3e7;

% Pg = 1.1e-4;
% rho = 1.0;
% mu = 1./Re;
% g = (Re*Pg)^(-2);
% sigma = 0.011269332539972


[Ah,Bh,Ch,Dh,z,w] = semhat(N);
x = (z-1.0)/2.0;
Lx2=1.0/2.0; Lxi=1./Lx2;

U = 1.0 - x.^2;
Dh = Lxi*Dh;
Bh = Lx2*Bh;
DU = Dh * U;
D2U = Dh * DU;

nh = N+1; Ih = speye(nh); 
R=Ih(2:end,:);

S = 0*(1:nh); S(end) = 1;
Sh = diag(S);
S = zeros(nh,1); S(end,1) = 1;

% Kuu0 = -mu/rho*R*Dh'*Dh'*Bh*Dh*Dh*R'...
%        -mu/rho*2*alpha^2*R*Dh'*Bh*Dh*R'...
%        -mu/rho*alpha^4*R*Bh*R';
% KuuU = -1i*alpha*(alpha^2*R*Bh*diag(U)*R'...
%                 +R*Dh'*Bh*diag(U)*Dh*R'
%                 -R*Dh'*Bh*diag(DU)*R');
% KuuS = -mu/rho*alpha^2*(R*Dh'*Sh*R'+R*Sh*Dh*R');
% Kuu = Kuu0 + KuuU + KuuS;
% Kua = -alpha^2*(g+sigma/rho*alpha^2-2i*mu/rho*alpha*DU(end))*R*S...
%       +1i*alpha*mu/rho*D2U(end)*R*Dh'*S;
% Kau = S'*R';
% Kaa = -1i*alpha*DU(end);
% Muu = R*Dh'*Bh*Dh*R' + alpha^2*R*Bh*R';
% Maa = 1;

Kuu0 = -R*Dh'*Dh'*Bh*Dh*Dh*R'...
       -2*alpha^2*R*Dh'*Bh*Dh*R'...
       -alpha^4*R*Bh*R';
KuuU = -1i*alpha*Re*(alpha^2*R*Bh*diag(U)*R'...
                +R*Dh'*Bh*diag(U)*Dh*R'...
                -R*Dh'*Bh*diag(DU)*R');
KuuS = -alpha^2*(R*Dh'*Sh*R'+R*Sh*Dh*R');
Kuu = Kuu0 + KuuU + KuuS;
Kua = -alpha^2*R*S*(Ga/Re + alpha^2/Ca - 2i*alpha*S'*DU)...
      +1i*alpha*R*Dh'*Sh*D2U;
Kau = S'*R';
Kaa = -1i*alpha*S'*U;
Muu = Re*(R*Dh'*Bh*Dh*R' + alpha^2*R*Bh*R');
Maa = 1;

K = [
  Kuu, Kua;
  Kau, Kaa;
];
M = [
  Muu, zeros(N, 1);
  zeros(1, N), Maa;
];



save -mat -v7 "KM.mat" K M



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


unstable = find(imag(c) > -0.00);

figure;
plot(abs(R'*vecs(1:end-1, unstable)), x, 'linewidth', 2)
labels = arrayfun(@(g) sprintf('c = %.2f + %.2fi, a =  %.2f + %.2fi', 
                                real(c(g)), 
                                imag(c(g)), 
                                real(vecs(end, g)),
                                imag(vecs(end, g)) ), unstable, 'UniformOutput', false);
xlabel('mag(V(y))');
ylabel('y');
legend(labels);
c(unstable)

N
res = norm(K*vecs-M*vecs*diag(gamma))
condK = cond(K)
condM = cond(M)
