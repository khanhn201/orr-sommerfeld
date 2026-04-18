clear all; close all;
output_precision(9);
set(0, "defaultAxesFontSize", 24)
set(0, "defaultTextFontSize", 24)
set(0, "defaultLineLineWidth", 2)

N = 250;

Re = 3e4;
ratio = 1e2;
Ca = 0.07;
Ga = 8.3e7;

rhos = [1, 1/ratio];
mus = [1/Re, 1/Re/ratio];
sigmas = [1/1, 1/1/ratio]*1/Re;
st = mus(1)/Ca; % surface tension
g = Ga*mus(1)^2/rhos(1)^2;
alpha = 1.0;

Bz = 1.0;
Bx = 0.0;


mode = 1; % Which of the unstable modes to output


[xs,umat,vmat,amat,gamma,f,A] = solve_os2p_mhd(alpha, N, rhos, mus, sigmas, st, g, Bx, Bz);
c = gamma*1i/alpha;
[~, unstable] = max(imag(c));
c_unstable_2p = c(unstable)
a = amat(unstable);
v = vmat(:, unstable);
u = umat(:, unstable);
v = reshape(v, [N+1, 2]);
u = reshape(u, [N+1, 2]);

eps_list = [];
err_list = [];
for k = -16:0
  eps = 10^(k/2);
  eps_list(end+1) = eps;
  [xs,umat,vmat,amat,gamma,f,U] = solve_osls_mhd(alpha, N, eps, rhos, mus, sigmas, st, g, Bx, Bz);
  % figure;
  % plot(abs(v), xs', 'linewidth', 2)
  % title('V'); hold on;
  % plot(abs(u), xs', 'linewidth', 2);
  c = gamma*1i/alpha;
  lambda = alpha*c;
  [~, unstable] = max(imag(c));
  c_unstable = c(unstable)
  err = abs(imag(c_unstable) - imag(c_unstable_2p)) / abs(imag(c_unstable_2p));
  err_list(end+1) = err;
  % a = amat(unstable);
  % v = vmat(:, unstable);
  % u = umat(:, unstable);
  % v = reshape(v, [N+1, 2]);
  % u = reshape(u, [N+1, 2]);
  % plot(abs(v), xs', 'linewidth', 2)
  % plot(abs(u), xs', 'linewidth', 2);
end
plot(eps_list, err_list)
