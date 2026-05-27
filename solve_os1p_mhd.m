function [xs,umat,vmat,amat,gamma,f, A] = solve_os1p_mhd(alpha, N, rhos, mus, sigmas, st, g, Bx, Bz)
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


  Nelem = 1;
  nh = N + 1;
  Ng = Nelem * (N-1) + 2;  % total global nodes

  xs = [];
  Uxplot = [];
  Uplot = [];
  for e = 1:Nelem
      rho = rhos(e);
      mu = mus(e);
      sigma = sigmas(e);

      if e == 1
          x = (z-1.0)/2; % x = [-1, 0]
          U = A1*cosh(k1*x) + A2*sinh(k1*x) + f*rho/sigma/Bz^2;
          U1 = U;
          Dh = 2*Dhh;
          Bh = 1/2*Bhh;
      else
          x = (z+1.0)/2; % x = [0, 1]
          U = A3*cosh(k2*x) + A4*sinh(k2*x) + f*rho/sigma/Bz^2;
          U2 = U;
          Dh = 2*Dhh;
          Bh = 1/2*Bhh;
      end
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

      if e == 1
          S = zeros(nh,1); S(end,1) = 1;
          Sh = diag(S);
      end
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
            -1.0*abs(S)*alpha^4*st... % 0.5 to counter multiplicity
            +1i*alpha*mu*Dh'*Sh*D2U;
      Kau = 1.0*abs(S');  % 0.5 to counter multiplicity

      Muu = rho*(Dh'*Bh*Dh + alpha^2*Bh);

      Kuu_block((e-1)*(N+1)+1:e*(N+1), (e-1)*(N+1)+1:e*(N+1)) = T'*Kuu*T;
      Muu_block((e-1)*(N+1)+1:e*(N+1), (e-1)*(N+1)+1:e*(N+1)) = T'*Muu*T;
      Kau_block(1,(e-1)*(N+1)+1:e*(N+1)) =  Kau*T;
      Kua_block((e-1)*(N+1)+1:e*(N+1),1) =  T'*Kua;
      T_block((e-1)*(N+1)+1:e*(N+1), (e-1)*(N+1)+1:e*(N+1)) = T;
      D_block((e-1)*(N+1)+1:e*(N+1), (e-1)*(N+1)+1:e*(N+1)) = Dh;
  end


  Ih = speye(Ng); 
  R=Ih(3:end,:);


  Kuu_global = R*Kuu_block*R';
  Kua_global = R*Kua_block;
  Kau_global = Kau_block*R';
  Muu_global = R*Muu_block*R';

  K = [
    Kuu_global, Kua_global;
    Kau_global, -1i*alpha*U1(end);
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
  vmat = T_block*R'*vecs(1:end-1, :);
  umat = 1i/alpha*D_block*vmat;
end
