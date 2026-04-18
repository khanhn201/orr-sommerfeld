function [xs,umat,vmat,amat,gamma,f,U] = solve_osls_mhd(alpha, N, eps, rhos, mus, sigmas, st, g, Bx, Bz)
  [Ah,Bhh,Ch,Dhh,z,w] = semhat(N);
  [xs,f,Uss] = solve_steady(alpha, N, eps, rhos, mus, sigmas, st, g, Bx, Bz);

  Nelem = 2;
  nh = N + 1;
  Ng = Nelem * (N-1) + 2;  % total global nodes
  Q = sparse([], [], [], Nelem*nh, Ng, 2*Ng);
  for e = 1:Nelem
      gidx = (e-1)*(nh-2) + (1:nh);
      Q((e-1)*nh + (1:nh), gidx) = speye(nh);
  end

  for e = 1:Nelem
      Dh = 2*Dhh;
      Bh = 1/2*Bhh;

      T = speye(nh);
      T(2, :) = Dh(1,:);
      T(end-1, :) = speye(nh)(end, :);
      T(end, :) = Dh(end,:);
      T = inv(T);

      if e == 1
          x = (z-1.0)/2; % x = [-1, 0]
      else
          x = (z+1.0)/2; % x = [0, 1]
      end
      U = Uss(:,e);
      phi = 0.5*(tanh(x/2/eps)+1);
      rho = (1-phi)*rhos(1) + phi*rhos(2);
      mu = (1-phi)*mus(1) + phi*mus(2);
      sigma = (1-phi)*sigmas(1) + phi*sigmas(2);
      delta = (1/4/eps)./(cosh(x/2/eps).^2);

      DU = Dh * U;
      D2U = Dh * DU;

      if e == 1
          S = zeros(nh,1); S(end,1) = 1;
          Sh = diag(S);
      else
          S = zeros(nh,1); S(1,1) = -1;
          Sh = diag(S);
      end
      Kuu0 = -Dh'*Dh'*Bh*diag(mu)*Dh*Dh...
             -4*alpha^2*Dh'*Bh*diag(mu)*Dh...
             -alpha^2*Bh*diag(mu)*Dh*Dh...
             -alpha^2*Dh'*Dh'*Bh*diag(mu)...
             -alpha^4*Bh*diag(mu);

      KuuU = 1i*alpha*Dh'*Bh*diag(rho)*diag(DU)...
             -1i*alpha*Dh'*Bh*diag(rho)*diag(U)*Dh...
             -1i*alpha^3*Bh*diag(rho)*diag(U);

      KuuB = -Bx^2*alpha^2*Bh*diag(sigma)... 
             + Bx*Bz*1i*alpha*Bh*diag(sigma)*Dh...
             - Bx*Bz*1i*alpha*Dh'*Bh*diag(sigma)...
             - Bz^2*Dh'*Bh*diag(sigma)*Dh;
      Kuu = Kuu0 + KuuU + KuuB;

      Kua = -alpha^4*st*Bh*delta...
            -1i*alpha*(rhos(1)-rhos(2))*f*Dh'*Bh*delta...
            -alpha^2*(rhos(1)-rhos(2))*g*Bh*delta...
            +1i*(mus(1)-mus(2))*(Dh'*Dh'*Bh*diag(delta)*DU...
                 + alpha^2*Bh*diag(delta)*DU)...
            +1i*alpha*(sigmas(1)-sigmas(2))*Bz^2*(Dh'*Bh*diag(delta)*U)...
            +alpha^2*(sigmas(1)-sigmas(2))*Bx*Bz*Bh*diag(delta)*U;

      Kau = 0.5*abs(S');  % 0.5 to counter multiplicity
      Muu = Dh'*Bh*diag(rho)*Dh + alpha^2*Bh*diag(rho);

      Kuu_block((e-1)*(N+1)+1:e*(N+1), (e-1)*(N+1)+1:e*(N+1)) = T'*Kuu*T;
      Muu_block((e-1)*(N+1)+1:e*(N+1), (e-1)*(N+1)+1:e*(N+1)) = T'*Muu*T;
      Kau_block(1,(e-1)*(N+1)+1:e*(N+1)) =  Kau*T;
      Kua_block((e-1)*(N+1)+1:e*(N+1),1) =  T'*Kua;
      T_block((e-1)*(N+1)+1:e*(N+1), (e-1)*(N+1)+1:e*(N+1)) = T;
      D_block((e-1)*(N+1)+1:e*(N+1), (e-1)*(N+1)+1:e*(N+1)) = Dh;
  end

  Ih = speye(Ng); 
  R=Ih(3:end-2,:);

  Kuu_global = R*Q'*Kuu_block*Q*R';
  Kua_global = R*Q'*Kua_block;
  Kau_global = Kau_block*Q*R';
  Muu_global = R*Q'*Muu_block*Q*R';

  K = [
    Kuu_global, Kua_global;
    Kau_global, -1i*alpha*Uss(end,1);
  ];
  M = [
    Muu_global, R*Q'*zeros(Nelem*nh, 1);
    zeros(1, Nelem*nh)*Q*R', 1;
  ];
  [vecs, gamma] = eig(K, M, 'vector');
  for k = 1:size(vecs,2)
    vecs(:,k) = vecs(:,k) / abs(vecs(end,k));
  end
  amat = vecs(end, :);
  vmat = T_block*Q*R'*vecs(1:end-1, :);
  umat = 1i/alpha*D_block*vmat;

  U1 = Uss(:,1);
  U2 = Uss(:,2);
  % test_error(alpha, N, rhos, mus, sigmas, st, g, umat, vmat, amat, gamma, Bx, Bz, U1, U2, Dh);
end

function [xs,f,Uss] = solve_steady(alpha, N, eps, rhos, mus, sigmas, st, g, Bx, Bz)
  [Ah,Bhh,Ch,Dhh,z,w] = semhat(N);

  k1 = sqrt(sigmas(1)/mus(1))*Bz;
  k2 = sqrt(sigmas(2)/mus(2))*Bz;
  rat = sqrt(sigmas(1)*mus(1)/sigmas(2)/mus(2));
  f = Bz^2/(...
      rhos(1)/sigmas(1)...
      -(rhos(2)/sigmas(2) + rat*rhos(1)/sigmas(1)/sinh(k1)*sinh(k2) + (rhos(1)/sigmas(1) - rhos(2)/sigmas(2))*cosh(k2))...
      /(cosh(k2) + rat*coth(k1)*sinh(k2))...
  );
  A1 = f/Bz^2*(-rhos(2)/sigmas(2) - rat*rhos(1)/sigmas(1)/sinh(k1)*sinh(k2) - (rhos(1)/sigmas(1)-rhos(2)/sigmas(2))*cosh(k2))...
        /(cosh(k2) + rat*coth(k1)*sinh(k2));
  A2 = (f*rhos(1)/sigmas(1)/Bz^2 + A1*cosh(k1))/sinh(k1);
  A3 = A1 + f/Bz^2*(rhos(1)/sigmas(1) - rhos(2)/sigmas(2));
  A4 = rat*A2;
  A = [A1, A2, A3, A4];

  Nelem = 2;
  nh = N + 1;
  Ng = Nelem * N + 1;  % total global nodes
  Q = sparse([], [], [], Nelem*nh, Ng, 2*Ng);
  for e = 1:Nelem
      gidx = (e-1)*(nh-1) + (1:nh);
      Q((e-1)*nh + (1:nh), gidx) = speye(nh);
  end
  xs = [];
  rho_arr = [];
  Uplot = [];
  for e = 1:Nelem
      Dh = 2*Dhh;
      Bh = 1/2*Bhh;

      if e == 1
          x = (z-1.0)/2; % x = [-1, 0]
          Usharp = A1*cosh(k1*x) + A2*sinh(k1*x) + f*rhos(1)/sigmas(1)/Bz^2;
      else
          x = (z+1.0)/2; % x = [0, 1]
          Usharp = A3*cosh(k2*x) + A4*sinh(k2*x) + f*rhos(2)/sigmas(2)/Bz^2;
      end
      xs = [xs; x'];
      Uplot = [Uplot;Usharp'];
      phi = 0.5*(tanh(x/2/eps)+1);
      rho = (1-phi)*rhos(1) + phi*rhos(2);
      mu = (1-phi)*mus(1) + phi*mus(2);
      sigma = (1-phi)*sigmas(1) + phi*sigmas(2);
      rho_arr = [rho_arr; rho];

      % K = -Dh'*Bh*diag(mu)*Dh - Bz^2*Bh*diag(sigma);
      K = Bh*diag(Dh*mu)*Dh + Bh*diag(mu)*Dh*Dh  - Bz^2*Bh*diag(sigma);
      K = -Dh'*Bh*diag(mu)*Dh  - Bz^2*Bh*diag(sigma);

      Kuu_block((e-1)*(N+1)+1:e*(N+1), (e-1)*(N+1)+1:e*(N+1)) = K;
      rhs((e-1)*(N+1)+1:e*(N+1)) = -Bh*f*rho;
  end

  Ih = speye(Ng); 
  R=Ih(2:end-1,:);
  Kuu_global = R*Q'*Kuu_block*Q*R';
  rhs = R*Q'*rhs';

  U = Q*R'*(Kuu_global \ rhs);
  Uss = reshape(U, [nh, 2]);
  % figure
  % plot(xs',Uss); hold on;
  % plot(xs',Uplot');hold off;
  % plot(xs',rho_arr);
end

function test_error(alpha, N, rhos, mus, sigmas, st, g, umat, vmat, amat, gamma_temp, Bx, Bz, U1, U2, Dh)
  [~, idx] = sort(real(gamma_temp), 'descend');
  unstable = idx(1:min(4, length(idx)));
  v1 = vmat(1:N+1,unstable);
  v2 = vmat(N+2:end,unstable);
  DU1 = Dh*U1;
  DU2 = Dh*U2;
  D2U1 = Dh*DU1;
  D2U2 = Dh*DU2;
  D2v1 = Dh*Dh*v1;
  D2v2 = Dh*Dh*v2;
  Dv1 = Dh*v1;
  Dv2 = Dh*v2;
  D3v1 = Dh*Dh*Dv1;
  D3v2 = Dh*Dh*Dv2;
  a = amat(unstable);
  gamma = gamma_temp(unstable).';

  shear_stress_error = (mus(1)*D2U1(end)-mus(2)*D2U2(1))*a...
                       + mus(1)*1i*alpha*v1(end, :) - mus(2)*1i*alpha*v2(1, :)...
                       - mus(1)/1i/alpha*D2v1(end, :) + mus(2)/1i/alpha*D2v2(1, :)
  normal_stress_error = (rhos(1)-rhos(2))*g*a...
                       + rhos(1)*(gamma/alpha^2 - U1(end)/1i/alpha).*Dv1(end, :)...
                       - rhos(2)*(gamma/alpha^2 - U2(1)/1i/alpha).*Dv2(1, :)...
                       - mus(1)/alpha^2*D3v1(end, :) + mus(2)/alpha^2*D3v2(1, :)...
                       + 3*mus(1)*Dv1(end, :) - 3*mus(2)*Dv2(1, :)...
                       + rhos(1)*DU1(end)/1i/alpha*v1(end,:)...
                       - rhos(2)*DU2(1)/1i/alpha*v2(1,:)...
                       + sigmas(1)/alpha^2*Bz^2*Dv1(end,:)...
                       - sigmas(2)/alpha^2*Bz^2*Dv2(1,:)...
                       - sigmas(1)/1i/alpha*Bz*Bx*v1(end,:)...
                       - sigmas(2)/1i/alpha*Bz*Bx*v2(1,:)...
                       + st*alpha^2*a
  kinematic_error = v1(end, :) - (gamma + 1i*alpha*U1(end)).*a

end
