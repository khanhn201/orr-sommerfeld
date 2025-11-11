function [z, w] = legendre_quadrature(p)
    % Returns Gauss-Legendre nodes and weights for L_p(x)
    beta = 0.5 ./ sqrt(1 - (2*(1:p)).^(-2));
    T = diag(beta,1) + diag(beta,-1);
    [V, D] = eig(T);
    z = diag(D);
    w = 2 * (V(1,:).^2)';
end

