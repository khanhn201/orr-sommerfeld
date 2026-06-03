function [uvec,vvec,wvec,pvec,phivec] = solve_lin_wall(N, rho, mu, sigma, sigma_w, W, U,Phi, By, f, alpha)
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
  R2D = kron(R, R); % Dirichlet min/max x and min/max y for velocity


  B2D = kron(Bh,Bh);
  A2D = kron(Bh,Ah) + kron(Ah,Bh);
  DX = kron(Ih,Dh);
  DY = kron(Dh,Ih);

  [zd,wd]=zwgl(N-1); Bd=diag(wd); % For pressure
  % [zd,wd]=zwgll(N); Bd=diag(wd); % For pressure
  JM=interp_mat(zd,z);
  JMT=interp_mat(z,zd);
  JM2D = kron(JM,JM);
  JMT2D = kron(JMT,JMT);
  Bd2D = kron(Bd,Bd);
  Ihp = speye(size(zd,1));
  Rp = kron(Ihp, Ihp);
  % Rp = Rp(2:end,:);  % Dirichlet first node for pressure

  Dph = deriv_mat(zd);
  DpX = kron(Ihp,Dph);
  DpY = kron(Dph,Ihp);


  % Velocity part
  Kuu = 1i*alpha*rho*B2D*diag(U) + mu*A2D + mu*alpha^2*B2D + sigma*By^2*B2D;
  Kuphi = -1i*alpha*sigma*By*B2D;
  % Kup = B2D*JMT2D*DpX;
  Kup = -DX'*B2D*JMT2D;
  Mu = -rho*B2D;

  Kvv = 1i*alpha*rho*B2D*diag(U) + mu*A2D + mu*alpha^2*B2D;
  Kvphi = 0*B2D;
  % Kvp = B2D*JMT2D*DpY;
  Kvp = -DY'*B2D*JMT2D;
  Mv = -rho*B2D;


  Kww = 1i*alpha*rho*B2D*diag(U) + mu*A2D + mu*alpha^2*B2D + sigma*By^2*B2D;
  Kwu = rho*B2D*diag(DX*U);
  Kwv = rho*B2D*diag(DY*U);
  Kwphi = sigma*By*B2D*DX;
  Kwp = 1i*alpha*B2D*JMT2D;
  Mw = -rho*B2D;

  % Kpu = Bd2D*JM2D*DX;
  % Kpv = Bd2D*JM2D*DY;
  % Kpw = 1i*alpha*Bd2D*JM2D;
  Kpu = -Bd2D*JM2D*DX;
  Kpv = -Bd2D*JM2D*DY;
  Kpw = -1i*alpha*Bd2D*JM2D;

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
      Kphiphi(idx_start:idx_end,idx_start:idx_end) = sigmal*(Ax+Ay) + sigmal*alpha^2*B2D;
  end
  end

  Kuphi_global = zeros(3*size(R2D,1), n2*9);
  Kphiu_global = zeros(n2*9, 3*size(R2D,1));
  Kuphi_global(:, n2*4+1:n2*5) = KVPhi;
  Kphiu_global(n2*4+1:n2*5, :) = KPhiV;

  [Q,glo_num]=set_tp_semq(Nelx,Nely,N);
  Rphi2D = speye(size(Q,2));
  % Rphi2D = Rphi2D(2:end,:); % Dirichlet first node for phi

  KPhiU = Rphi2D*Q'*Kphiu_global;
  KUPhi = Kuphi_global*Q*Rphi2D';
  KPhiPhi = Rphi2D*Q'*Kphiphi*Q*Rphi2D';
  K = [
    KVV  , KVP                              , KUPhi;
    KPV  , 0*Rp*Bd2D*Rp'                     , zeros(size(KPV,1), size(Rphi2D',2));
    KPhiU, zeros(size(Rphi2D,1),size(KVP,2)), KPhiPhi;
  ];
  npphi = size(K,1)-size(MV,1);
  M = [
    MV                     ,zeros(size(MV,1),npphi);
    zeros(npphi,size(MV,1)),zeros(npphi,npphi);
  ];
  % pinvM = pinv(K);
  % [vecs, gamma] = eigs(pinvM*K, 5, "lr");
  % gamma = diag(gamma)
  % [vecs, gamma] = eig(pinvM*K, 'vector');
  
  % [vecs, gamma] = eigs(pinvM*K, pinvM*M, 5, 'lr');
  % size(vecs)

  KVV2 = KVV - KUPhi*(KPhiPhi \ KPhiU);

  M2 = (KPV*(KVV2\KVP)) \ (KPV*(KVV2\MV));

  M3 = MV - KVP*M2;
  S = KPV*(KVV2\KVP);
  Sinv = inv(S);
  tmp = eye(size(KUPhi,1)) - KVP*Sinv*(KPV*inv(KVV2));

  % [vecs, gamma] = eigs(KVV2, M3, 5, "lr");
  % gamma = diag(gamma)
  % M4 = KVV2\M3;
  [vecs, gamma] = eig(KVV2, M3, 'vector');
  % gamma = 1./gamma
  res = zeros(size(gamma));

  for i = 1:length(gamma)
      v = vecs(:,i);
      r = KVV2*v - gamma(i)*M3*v;
      res(i) = norm(r) / ...
          (norm(KVV2*v) + norm(M3*v) + eps);
  end

  good = res < 1e-8;

  gamma = gamma(good);
  vecs  = vecs(:,good);
  res   = res(good);
  % gamma

  c = gamma*1i/alpha;
  plot(real(c), imag(c), 'o');
  [~, unstable] = max(imag(c));
  c(unstable)
  vu = vecs(:, unstable);
  uvec   =      R2D'*vu(0*(N-1)^2+1       :1*(N-1)^2        );
  vvec   =      R2D'*vu(1*(N-1)^2+1       :2*(N-1)^2        );
  wvec   =      R2D'*vu(2*(N-1)^2+1       :3*(N-1)^2        );
  % pvec   =       Rp'*vu(3*(N-1)^2+1       :3*(N-1)^2+(N+1)^2-1);
  % phivec = Q*Rphi2D'*vu(3*(N-1)^2+(N+1)^2 :end              );
  pvec = Rp'*gamma(unstable)*M2*vu;
  phivec = -Q*Rphi2D'*(KPhiPhi \ KPhiU)*vu;
end
