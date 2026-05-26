% Add pressure and phi restriction? or subtract mean out every step
% Solve eigenvalue with Schur compl?
function [uvec,vvec,wvec,pvec,phivec] = solve_lin(N, rho, mu, sigma, sigma_w, W, U,Phi, By, f, alpha)
  [Ah,Bh,Ch,Dh,z,w] = semhat(N);
  Nelx = 3;
  Nely = 3;
  Nelem = Nelx*Nely;

  n = N + 1;
  n2 = n^2;
  Ng = n2 * 13;  % total: 3 for v, 1 for p, 9 for phi

  K = zeros(Ng, Ng);
  M = zeros(Ng, 1);

  Ih = speye(N+1);
  R = Ih(2:end-1,:);
  R2D = kron(R, R);

  Ih = speye(N+1);
  Rp = kron(Ih, Ih);
  Rp = Rp(2:end,:);

  B2D = kron(Bh,Bh);
  A2D = kron(Bh,Ah) + kron(Ah,Bh);
  DX = kron(Ih,Dh);
  DY = kron(Dh,Ih);


  % Velocity part
  Kuu = 1i*alpha*rho*B2D*diag(U) + mu*A2D + sigma*By^2*B2D;
  Kuphi = -1i*alpha*sigma*By*B2D;
  Kup = B2D*DX;
  Mu = -rho*B2D;

  Kvv = 1i*alpha*rho*B2D*diag(U) + mu*A2D;
  Kvphi = 0*B2D;
  Kvp = B2D*DY;
  Mv = -rho*B2D;


  Kww = 1i*alpha*rho*B2D*diag(U) + mu*A2D + sigma*By^2*B2D;
  Kwu = rho*B2D*DX*diag(U);
  Kwv = rho*B2D*DY*diag(U);
  Kwphi = sigma*By*B2D*DX;
  Kwp = 1i*alpha*B2D;
  Mw = -rho*B2D;

  Kpu = B2D*DX;
  Kpv = B2D*DY;
  Kpw = 1i*alpha*B2D;

  Kphiu = 1i*alpha*sigma*By*B2D;
  Kphiv = 0*B2D;
  Kphiw = -sigma*By*B2D*DX;

  KVV = [
        R2D*Kuu*R2D', 0*R2D*Kuu*R2D', 0*R2D*Kuu*R2D';
      0*R2D*Kvv*R2D',   R2D*Kvv*R2D', 0*R2D*Kvv*R2D';
        R2D*Kwu*R2D',   R2D*Kwv*R2D',   R2D*Kww*R2D';
  ];
  KVPhi = [
      R2D*Kuphi;
      R2D*Kvphi;
      R2D*Kwphi;
  ];
  KVP = [
      R2D*Kup*Rp';  
      R2D*Kvp*Rp';  
      R2D*Kwp*Rp';  
  ];
  KPhiV = [Kphiu*R2D', Kphiv*R2D', Kphiw*R2D'];
  KPV = [Rp*Kpu*R2D', Rp*Kpv*R2D', Rp*Kpw*R2D'];
  MV = [
    R2D*Mu*R2D', 0*R2D*Mu*R2D', 0*R2D*Mu*R2D';
    0*R2D*Mu*R2D', R2D*Mv*R2D', 0*R2D*Mu*R2D';
    0*R2D*Mu*R2D', 0*R2D*Mu*R2D', R2D*Mw*R2D';
  ];

  % phi part
  for ey = 1:Nely
  for ex = 1:Nelx
      Ih = speye(N+1);
      R = Ih(2:end-1,:);

      sigmal = sigma;
      if ey != 2 || ex != 2
         sigmal = sigma_w;
      end
      Ax = kron(Bh,Ah);
      Ay = kron(Ah,Bh);
      if ex != 2
         Ax = 1/W*Ax;
      end
      if ey != 2
         Ay = 1/W*Ay;
      end

      e = (ey-1)*Nelx + ex;
      idx_start = (e-1)*n2+1;
      idx_end = e*n2;
      Kphiphi(idx_start:idx_end,idx_start:idx_end) = sigmal*(Ax+Ay);
  end
  end

  Kuphi_global = zeros(3*size(R2D,1), n2*9);
  Kphiu_global = zeros(n2*9, 3*size(R2D,1));
  Kuphi_global(:, n2*4+1:n2*5) = KVPhi;
  Kphiu_global(n2*4+1:n2*5, :) = KPhiV;

  [Q,glo_num]=set_tp_semq(Nelx,Nely,N);
  Rphi2D = speye(size(Q,2));
  Rphi2D = Rphi2D(2:end,:);

  K = [
    KVV, KVP, Kuphi_global*Q*Rphi2D';
    KPV, 0*Rp*B2D*Rp', zeros(size(KPV,1), size(Rphi2D',2));
    Rphi2D*Q'*Kphiu_global, zeros(size(Rphi2D,1),size(KVP,2)), Rphi2D*Q'*Kphiphi*Q*Rphi2D';
  ];
  npphi = size(K,1)-size(MV,1);
  M = [
    MV                     ,zeros(size(MV,1),npphi);
    zeros(npphi,size(MV,1)),zeros(npphi,npphi);
  ];
  pinvM = pinv(K);
  % [vecs, gamma] = eigs(pinvM*K, 5, "lr");
  % gamma = diag(gamma)
  % [vecs, gamma] = eig(pinvM*K, 'vector');
  
  [vecs, gamma] = eigs(pinvM*K, pinvM*M, 5, 'lr');
  % size(vecs)
  c = gamma*1i/alpha;
  [~, unstable] = max(imag(c));
  vu = vecs(:, unstable);
  vu = vu(:);
  
  uvec   =      R2D'*vu(0*(N-1)^2+1       :1*(N-1)^2        );
  vvec   =      R2D'*vu(1*(N-1)^2+1       :2*(N-1)^2        );
  wvec   =      R2D'*vu(2*(N-1)^2+1       :3*(N-1)^2        );
  pvec   =       Rp'*vu(3*(N-1)^2+1       :3*(N-1)^2+(N+1)^2-1);
  phivec = Q*Rphi2D'*vu(3*(N-1)^2+(N+1)^2 :end              );
end
