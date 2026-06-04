function [uvec,vvec,wvec,pvec,phivec,gamma] = solve_lin(N, rho, mu, sigma, U,Phi, By, alpha)
  [Ah,Bh,Ch,Dh,z,w] = semhat(N);
  Nelx = 1;
  Nely = 1;
  Nelem = Nelx*Nely;

  n = N + 1;
  n2 = n^2;
  Ih = speye(N+1);
  R = Ih(2:end-1,:);
  R2D = kron(R, R); % Dirichlet min/max x and min/max y for velocity


  B2D = kron(Bh,Bh);
  A2D = kron(Bh,Ah) + kron(Ah,Bh);
  DX = kron(Ih,Dh);
  DY = kron(Dh,Ih);

  [zd,wd]=zwgl(N-1); Bd=diag(wd); % For pressure

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
  Kphiphi = sigma*A2D + sigma*alpha^2*B2D;

  Rphi = Ih(2:end-1,:); % Perfectly conducting
  Rphi2D = kron(Rphi, Ih);
  % Rphi2D = Rphi2D(1:end,:);


  KPhiU = Rphi2D*KPhiV;
  KUPhi = KVPhi*Rphi2D';
  KPhiPhi = Rphi2D*Kphiphi*Rphi2D';
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

  KVV2 = KVV - KUPhi*(KPhiPhi \ KPhiU);
  M2 = (KPV*(KVV2\KVP)) \ (KPV*(KVV2\MV));

  M3 = MV - KVP*M2;

  % [vecs, gamma] = eigs(KVV2, M3, 5, "lr");
  % gamma = diag(gamma)
  % M4 = KVV2\M3;
  % [i,j,v] = find(KVV2);
  % rv = real(v);
  % iv = imag(v);
  % save('A_i.dat','i','-ascii');
  % save('A_j.dat','j','-ascii');
  % save('A_v_real.dat','rv','-ascii');
  % save('A_v_imag.dat','iv','-ascii');
  % [i,j,v] = find(M3);
  % rv = real(v);
  % iv = imag(v);
  % save('B_i.dat','i','-ascii');
  % save('B_j.dat','j','-ascii');
  % save('B_v_real.dat','rv','-ascii');
  % save('B_v_imag.dat','iv','-ascii');


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

  [re_sorted, idx] = sort(real(gamma), 'descend');
  disp(gamma(idx(1:min(5,end))))
  [~, unstable] = max(real(gamma));
  gamma = gamma(unstable);
  vu = vecs(:, unstable);
  uvec   =      R2D'*vu(0*(N-1)^2+1       :1*(N-1)^2        );
  vvec   =      R2D'*vu(1*(N-1)^2+1       :2*(N-1)^2        );
  wvec   =      R2D'*vu(2*(N-1)^2+1       :3*(N-1)^2        );
  % pvec   =       Rp'*vu(3*(N-1)^2+1       :3*(N-1)^2+(N+1)^2-1);
  % phivec = Rphi2D'*vu(3*(N-1)^2+(N+1)^2 :end              );
  pvec = Rp'*gamma*M2*vu;
  phivec = -Rphi2D'*(KPhiPhi \ KPhiU)*vu;


  check_error(N, alpha, uvec,vvec,wvec,pvec,phivec, By);
end

function check_error(N, alpha, uvec, vvec, wvec, pvec, phivec, By)
  [Ah,Bh,Ch,Dh,z,w] = semhat(N);
  [zd,wd]=zwgl(N-1); Bd=diag(wd); % For pressure
  B2D = kron(Bh,Bh);
  A2D = kron(Bh,Ah) + kron(Ah,Bh);


  Bd2D = kron(Bd,Bd);
  JMT=interp_mat(z,zd);
  JMT2D = kron(JMT,JMT);
  JM=interp_mat(zd,z);
  JM2D = kron(JM,JM);

  Ih = speye(N+1);
  DX = kron(Ih,Dh);
  DY = kron(Dh,Ih);
  div = DX*uvec + DY*vvec + 1i*alpha*wvec;
  div = norm(Bd2D*JM2D*div)


  philap = - DX'*B2D*DX*phivec - DY'*B2D*DY*phivec -alpha^2*B2D*phivec + By*B2D*DX*wvec - By*1i*alpha*B2D*uvec;
  % philap = -A2D*phivec -alpha^2*B2D*phivec + By*B2D*DX*wvec - By*1i*alpha*B2D*uvec;
  philap = norm(philap)
end
