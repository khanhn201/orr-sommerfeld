
[Ah,Bh,Ch,Dh,z,w] = semhat(N);
Nx = 3;
Ny = 3;
n  = N + 1;
n2 = n^2;
x_edges = [-1-W, -1, 1, 1+W];
y_edges = x_edges;
figure; hold on;
for e = 1:9
    idx = (e-1)*n2 + (1:n2);
    Phie = Phi_base(idx);

    ex = mod(e-1,3);
    ey = floor((e-1)/3);

    xa = x_edges(ex+1);
    xb = x_edges(ex+2);
    ya = y_edges(ey+1);
    yb = y_edges(ey+2);

    [Zy, Zx] = meshgrid(z,z);
    X = (xb-xa)/2 * (Zx + 1) + xa;
    Y = (yb-ya)/2 * (Zy + 1) + ya;

    sigmal = sigma;
    ex = ex+1;
    ey = ey+1;
    if (ey != 2) || (ex != 2)
       sigmal = sigma_w;
    end

    Dhx = Dh;
    if ex != 2
       Dhx = Dhx/W*2.0;
    end

    Ih = speye(N+1);
    Phie = sigmal*kron(Ih,Dhx)*Phie;
    Phie = reshape(Phie, N+1, N+1);
    surf(X, Y, real(Phie), 'EdgeColor', 'none');
end
title('\phi');
colorbar;
view(3);
