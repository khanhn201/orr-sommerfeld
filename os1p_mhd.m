N = 250;


Re = 3e4;
Ca = 0.07;
Ga = 8.3e7;

rho = 1;
mu = 1/Re;
sigma = 1/Re;
st = mu/Ca; % surface tension
g = Ga*mu^2/rho^2;
alpha = 1.13;

Bz = 5.0;
Bx = 0.0;

[xs,umat,vmat,amat,gamma,f] = solve_os1p_hartmann(alpha, N, rho, mu, sigma, st, g, Bx, Bz);
