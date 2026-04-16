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

eps = 0.01; % interface thickness
mode = 1; % Which of the unstable modes to output




[xs,umat,vmat,amat,gamma,f,A] = solve_os2p_hartmann(alpha, N, rhos, mus, sigmas, st, g, Bx, Bz);
c = gamma*1i/alpha;
lambda = alpha*c;
