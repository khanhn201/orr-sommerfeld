load("modeII_wall.mat")
U = reshape(uvec,N+1,N+1);
[y,x] = meshgrid(z,z);
figure;
surf(x,y,abs(U));
title('u');
U = reshape(vvec,N+1,N+1);
[y,x] = meshgrid(z,z);
figure;
surf(x,y,abs(U));
title('v');
U = reshape(wvec,N+1,N+1);
[y,x] = meshgrid(z,z);
figure;
surf(x,y,abs(U));
title('w');
P = reshape(pvec,N-1,N-1);
[zd,wd]=zwgl(N-1);
[y,x] = meshgrid(zd,zd);
figure;
surf(x,y,abs(P));
title('p');
plot_phi_3b3(N, phivec, W)

Nf = 10; % lx1
el_pos = [-1.0, -0.95, -0.9, -0.5, -0.25, 0.0, 0.25, 0.5, 0.9, 0.95, 1.0];

nely = length(el_pos) - 1;
[zff,wf] = zwgll(Nf-1);
zf = zeros(nely*Nf,1);
delta_el = diff(el_pos);
el_center = (el_pos(2:end) + el_pos(1:end-1))/2.0;


uf = [];
vf = [];
wf = [];
phif = [];
% wBasef = [];
% phiBasef = [];
for j=1:nely
for i=1:nely
    zf((i-1)*Nf+1:i*Nf) = el_center(i) + zff*delta_el(i)/2.0;
    zfx =  el_center(i) + zff*delta_el(i)/2.0;
    zfy =  el_center(j) + zff*delta_el(j)/2.0;
    Jx = interp_mat(zfx,z);
    Jy = interp_mat(zfy,z);
    J2D = kron(Jy,Jx);
    uf   = [uf; J2D*uvec];
    vf   = [vf; J2D*vvec];
    wf   = [wf; J2D*wvec];
    % phif = [phif; J2D*phivec];
    wBasef = [wBasef; J2D*U_base];
    % phiBasef = [phiBasef; J2D*Phi_base];

end
end
data = [real(uf), imag(uf), real(vf), imag(vf), real(wf), imag(wf), real(phif), imag(phif)];


fid = fopen("u2.txt", "w");
fprintf(fid, "%d  ! nelxy\n", nely);
fprintf(fid, "% .16e  % .16e  % .16e  ! rho, mu, sigma\n", rho, mu, sigma);
fprintf(fid, "% .16e  % .16e  ! f, By\n", f, By);
fprintf(fid, "% .16e  ! alpha\n", alpha);
fprintf(fid, "% .16e  % .16e  ! gamma\n", real(gamma), imag(gamma));
fprintf(fid, "% .16e  % .16e  % .16e  % .16e  % .16e  % .16e  % .16e  % .16e\n", data.');
fclose(fid);


data = [wBasef, phiBasef];
fid = fopen("ubase.txt", "w");
fprintf(fid, "%d  ! nelxy\n", nely);
fprintf(fid, "% .16e  % .16e  % .16e  ! rho, mu, sigma\n", rho, mu, sigma);
fprintf(fid, "% .16e  % .16e  ! f, By\n", f, By);
fprintf(fid, "% .16e  ! alpha\n", alpha);
fprintf(fid, "% .16e  % .16e  ! gamma\n", real(gamma), imag(gamma));
fprintf(fid, "% .16e  % .16e\n", data.');
fclose(fid);

