function [uvec,phivec] = solve_steady_wall(N, mu, sigma, sigma_w, W, By, f)
  [Ah,Bh,Ch,Dh,z,w] = semhat(N);
  Nelx = 3;
  Nely = 3;
  Nelem = Nelx*Nely;

  n = N + 1;
  n2 = n^2;
  Ng = n2 * 10;  % total: 1 for v, 9 for phi

  K = zeros(Ng, Ng);
  M = zeros(Ng, 1);

  Ih = speye(N+1);
  R = Ih(2:end-1,:);
  R2D = kron(R, R);

  B2D = kron(Bh,Bh);
  A2D = kron(Bh,Ah) + kron(Ah,Bh);
  DX = kron(Ih,Dh);


  % Velocity part
  Kuu = mu*A2D + sigma*By^2*B2D;
  Kuphi = sigma*By*B2D*DX;
  Kphiu = -sigma*By*B2D*DX;


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

  Kuphi_global = zeros(n2, n2*9);
  Kphiu_global = zeros(n2*9, n2);
  Kuphi_global(:, n2*4+1:n2*5) = Kuphi;
  Kphiu_global(n2*4+1:n2*5, :) = Kphiu;

  [Q,glo_num]=set_tp_semq(Nelx,Nely,N);
  Rphi2D = speye(size(Q,2));
  % Rphi2D = Rphi2D(2:end,:);
  % spy(Q)

  K = [
    R2D*Kuu*R2D', R2D*Kuphi_global*Q*Rphi2D';
    Rphi2D*Q'*Kphiu_global*R2D', Rphi2D*Q'*Kphiphi*Q*Rphi2D';
  ];
  rhs = [
    R2D*B2D*(f*ones(n2, 1));
    Rphi2D*Q'*zeros(n2*9,1)
  ];
  % sol = K\rhs;
  M = spdiags(diag(K),0,size(K,1),size(K,2));
  restart = [];
  tol     = 1e-8;
  maxit   = 2000;
  [sol,flag,relres,iter,resvec] = gmres( ...
      K, rhs, restart, tol, maxit, M);
  iter
  relres

  uvec = R2D'*sol(1:(N-1)^2);
  phivec = Q*Rphi2D'*sol((N-1)^2+1:end);
end
