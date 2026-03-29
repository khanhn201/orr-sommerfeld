function [xs,umat,vmat,amat,gamma,f,r] = solve_os2p(alpha, N, rhos, mus, st, g, h, top_bc, accel_type, Bx, Bz)
  if top_bc == 'W'
      if accel_type == 'S' % same acceleration on both domains
          f = 2*(mus(2) + mus(1)*h)/(rhos(1)*h + rhos(2)*h^2);
          r = 1.0;
      elseif accel_type == 'F' % Prescribed accel on fluid; ratio is computed correspondingly
          f = mus(1)/rhos(1)*2;
          r = (2*(mus(2)+mus(1)*h)/f-rhos(1)*h)/(h^2*rhos(2));
      else % Prescribed f2/f1 ratio
          r = 1/h;
          f = 2*(mus(2)+mus(1)*h)/(r*(h^2*rhos(2)) + rhos(1)*h);
      end
  else % Symmetry top bc
      if accel_type == 'S' % same acceleration on both domains
          f = mus(1)/(rhos(2)*h + rhos(1)/2)
          r = 1.0
      elseif accel_type == 'F' % Prescribed accel on fluid; ratio is computed correspondingly
          f = mus(1)/rhos(1)*2;
          r = 0.0;
      else % Prescribed f2/f1 ratio
          r = 1/h;
          f = mus(1)/mus(2)/(rhos(1)/mus(2)/2 + h*r*rhos(2)/mus(2));
      end
  end

  [Ah,Bhh,Ch,Dhh,z,w] = semhat(N);

  Nelem = 2;
  nh = N + 1;
  Ng = Nelem * (N-1) + 2;  % total global nodes
  Q = sparse([], [], [], Nelem*nh, Ng, 2*Ng);
  for e = 1:Nelem
      gidx = (e-1)*(nh-2) + (1:nh);
      Q((e-1)*nh + (1:nh), gidx) = speye(nh);
  end

  xs = [];
  Uxplot = [];
  Uplot = [];
  for e = 1:Nelem
      rho = rhos(e);
      mu = mus(e);

      if e == 1
          x = (z-1.0)/2; % x = [-1, 0]
          U = 1.0 + (1-f/2*rho/mu)*x - (f/2*rho/mu)*x.^2;
          U1 = U;
          Dh = 2*Dhh;
          Bh = 1/2*Bhh;
      else
          x = (z+1.0)/2*h; % x = [0, h]
          U = 1.0 + (mus(1)/mu - f/2*rhos(1)/mu)*x - (r*f/2*rho/mu)*x.^2;
          U2 = U;
          Dh = 2/h*Dhh;
          Bh = h/2*Bhh;
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
      else
          S = zeros(nh,1); S(1,1) = -1;
          Sh = diag(S);
      end
      Kuu0 = -mu*Dh'*Dh'*Bh*Dh*Dh...
             -mu*2*alpha^2*Dh'*Bh*Dh...
             -mu*alpha^4*Bh;
      KuuU = -1i*alpha*rho*(alpha^2*Bh*diag(U)...
                      +Dh'*Bh*diag(U)*Dh...
                      -Dh'*Bh*diag(DU));
      KuuS = -alpha^2*mu*(Dh'*Sh + Sh*Dh);
      Kuu = Kuu0 + KuuU + KuuS;

      Kua = -alpha^2*S*rho*g...
            -0.5*abs(S)*alpha^4*st... % 0.5 to counter multiplicity
            +1i*alpha*mu*Dh'*Sh*D2U;
      Kau = 0.5*abs(S');  % 0.5 to counter multiplicity

      Muu = rho*(Dh'*Bh*Dh + alpha^2*Bh);

      Kuu_block((e-1)*(N+1)+1:e*(N+1), (e-1)*(N+1)+1:e*(N+1)) = T'*Kuu*T;
      Muu_block((e-1)*(N+1)+1:e*(N+1), (e-1)*(N+1)+1:e*(N+1)) = T'*Muu*T;
      Kau_block(1,(e-1)*(N+1)+1:e*(N+1)) =  Kau*T;
      Kua_block((e-1)*(N+1)+1:e*(N+1),1) =  T'*Kua;
      T_block((e-1)*(N+1)+1:e*(N+1), (e-1)*(N+1)+1:e*(N+1)) = T;
      D_block((e-1)*(N+1)+1:e*(N+1), (e-1)*(N+1)+1:e*(N+1)) = Dh;
  end


  Ih = speye(Ng); 
  if top_bc == 'W'
      R=Ih(3:end-2,:);
  else
      R=Ih(3:end-1,:);
      R(end, :) = Ih(end,:);
  end


  Kuu_global = R*Q'*Kuu_block*Q*R';
  Kua_global = R*Q'*Kua_block;
  Kau_global = Kau_block*Q*R';
  Muu_global = R*Q'*Muu_block*Q*R';

  K = [
    Kuu_global, Kua_global;
    Kau_global, -1i*alpha*1.0; % U(0) = 1.0
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
  test_error(alpha, N, rhos, mus, st, g, umat, vmat, amat, gamma, U1, U2, Dh);
end

function test_error(alpha, N, rhos, mus, st, g, umat, vmat, amat, gamma_temp, U1, U2, Dh)
  unstable = find(real(gamma_temp) > 0.0);
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
                       + st*alpha^2*a

  kinematic_error = v1(end, :) - (gamma + 1i*alpha*U1(end)).*a
end
