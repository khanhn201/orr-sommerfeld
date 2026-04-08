function [xs,umat,vmat,amat,gamma,f] = solve_os1p_hartmann(alpha, N, rho, mu, sigma, st, g, Bx, Bz)
  [Ah,Bhh,Ch,Dhh,z,w] = semhat(N);

  k1 = sqrt(sigma/mu)*Bz;
  f = Bz^2/rho*sigma*cosh(k1)/(cosh(k1)-1);


  Nelem = 1;
  nh = N + 1;
  Ng = Nelem * (N-1) + 2;  % total global nodes
  xs = [];
  Uxplot = [];
  Uplot = [];
      x = (z-1.0)/2; % x = [-1, 0]
      U = -cosh(k1*x)/(cosh(k1)-1) + f*rho/sigma/Bz^2;
      U1 = U;
      Dh = 2*Dhh;
      Bh = 1/2*Bhh;
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

      S = zeros(nh,1); S(end,1) = 1;
      Sh = diag(S);
      Kuu0 = -mu*Dh'*Dh'*Bh*Dh*Dh...
             -mu*2*alpha^2*Dh'*Bh*Dh...
             -mu*alpha^4*Bh;
      KuuU = -1i*alpha*rho*(alpha^2*Bh*diag(U)...
                      +Dh'*Bh*diag(U)*Dh...
                      -Dh'*Bh*diag(DU));
      KuuS = -alpha^2*mu*(Dh'*Sh + Sh*Dh);
      KuuB = sigma*(-Bx^2*alpha^2*Bh + Bx*Bz*1i*alpha*(Bh*Dh - Dh'*Bh) - Bz^2*Dh'*Bh*Dh);
      Kuu = Kuu0 + KuuU + KuuS + KuuB;

      Kua = -alpha^2*S*rho*g...
            -abs(S)*alpha^4*st... % 0.5 to counter multiplicity
            +2*1i*alpha^3*mu*Sh*DU...
            +sigma*Bx*Bz*Sh*U...
            +1i*alpha*mu*Dh'*Sh*D2U;
      Kau = abs(S');  % 0.5 to counter multiplicity

      Muu = rho*(Dh'*Bh*Dh + alpha^2*Bh);

      Kuu = T'*Kuu*T;
      Muu = T'*Muu*T;
      Kau =  Kau*T;
      Kua =  T'*Kua;


  Ih = speye(Ng); 
  R=Ih(3:end,:);


  Kuu_global = R*Kuu*R';
  Kua_global = R*Kua;
  Kau_global = Kau*R';
  Muu_global = R*Muu*R';

  K = [
    Kuu_global, Kua_global;
    Kau_global, -1i*alpha*1.0; % U(0) = 1.0
  ];
  M = [
    Muu_global, R*zeros(Nelem*nh, 1);
    zeros(1, Nelem*nh)*R', 1;
  ];


  [vecs, gamma] = eig(K, M, 'vector');
  for k = 1:size(vecs,2)
    vecs(:,k) = vecs(:,k) / abs(vecs(end,k));
  end
  amat = vecs(end, :);
  vmat = T*R'*vecs(1:end-1, :);
  umat = 1i/alpha*Dh*vmat;
  % test_error(alpha, N, rhos, mus, sigmas, st, g, umat, vmat, amat, gamma, Bx, Bz, U1, U2, Dh);
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
