function [uvec,phivec,f] = solve_steady(N, mu, sigma, By)
  f = 1;
  [Ah,Bh,Ch,Dh,z,w] = semhat(N);
  Ih = speye(N+1);
  R = Ih(2:end-1,:);
  
  n = size(R, 1);

  Rphi = Ih(2:end-1,:); % Perfectly conducting
  Rphi2D = kron(Rphi, Ih);
  % Rphi2D = Rphi2D(2:end,:);
  R2D = kron(R, R);
  B2D = kron(Bh,Bh);
  A2D = kron(Bh,Ah) + kron(Ah,Bh);
  DX = kron(Ih,Dh);

  K = [
    R2D*(mu*A2D + sigma*By^2*B2D)*R2D', sigma*By*R2D*B2D*DX*Rphi2D';
    -Rphi2D*By*B2D*DX*R2D', Rphi2D*A2D*Rphi2D';
  ];

  rhs = [
    R2D*B2D*(f*ones((N+1)^2, 1));
    Rphi2D*B2D*zeros((N+1)^2,1)
  ];

  % sol = K\rhs;
  M = spdiags(diag(K),0,size(K,1),size(K,2));
  restart = [];
  tol     = 1e-8;
  maxit   = min(2000, N^2);
  [sol,flag,relres,iter,resvec] = gmres( ...
      K, rhs, restart, tol, maxit, M);
  iter
  relres


  uvec = R2D'*sol(1:n*n);
  phivec = Rphi2D'*sol(n*n+1:end);
  umax = max(uvec);
  uvec = uvec/umax;
  phivec = phivec/umax;
  f = f/umax;
  rhs = [
    R2D*B2D*(f*ones((N+1)^2, 1));
    Rphi2D*B2D*zeros((N+1)^2,1)
  ];
  [sol,flag,relres,iter,resvec] = gmres( ...
      K, rhs, restart, tol, maxit, M);
  iter
  relres


  uvec = R2D'*sol(1:n*n);
  phivec = Rphi2D'*sol(n*n+1:end);
end
