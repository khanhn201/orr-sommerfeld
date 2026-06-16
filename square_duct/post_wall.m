output_precision(9);
load("ha100n40.mat")
[Ah,Bh,Ch,Dh,z,w] = semhat(N);

U = reshape(uvec,N+1,N+1);
[y,x] = meshgrid(z,z);
figure;
surf(x,y,real(U));
title('u');
% U = reshape(vvec,N+1,N+1);
% [y,x] = meshgrid(z,z);
% figure;
% surf(x,y,real(U));
% title('v');
% U = reshape(wvec,N+1,N+1);
% [y,x] = meshgrid(z,z);
% figure;
% surf(x,y,real(U));
% title('w');
% P = reshape(pvec,N-1,N-1);
% [zd,wd]=zwgl(N-1);
% [y,x] = meshgrid(zd,zd);
% figure;
% surf(x,y,real(P));
% title('p');
plot_phi_3b3(N, phivec, W)

Nf = 8; % lx1
el_pos = [-1.0, -0.95, -0.9, -0.5, -0.25, 0.0, 0.25, 0.5, 0.9, 0.95, 1.0];
mid = (el_pos(1:end-1) + el_pos(2:end))/2;
out = zeros(1, 2*numel(el_pos)-1);
out(1:2:end) = el_pos;
out(2:2:end) = mid;
el_pos = out;
mid = (el_pos(1:end-1) + el_pos(2:end))/2;
out = zeros(1, 2*numel(el_pos)-1);
out(1:2:end) = el_pos;
out(2:2:end) = mid;
el_pos = out;

nelxy = length(el_pos) - 1;
[zff,wf] = zwgll(Nf-1);
delta_el = diff(el_pos);
el_center = (el_pos(2:end) + el_pos(1:end-1))/2.0;


nelz = 1;
z_length = 10;
data_r = zeros(nelxy*nelxy*nelz, Nf*Nf*Nf, 6); % xyzuvw
data_i = zeros(nelxy*nelxy*nelz, Nf*Nf*Nf, 6); % xyzuvw
data_base = zeros(nelxy*nelxy*nelz, Nf*Nf*Nf, 6); % xyzuvw
for j=1:nelxy
for i=1:nelxy
    el = i + (j-1)*nelxy;
    zfx =  el_center(i) + zff*delta_el(i)/2.0;
    zfy =  el_center(j) + zff*delta_el(j)/2.0;
    zfz =  z_length/2 + zff*z_length/2.0; % [0, z_length]

    Jx = interp_mat(zfx,z);
    Jy = interp_mat(zfy,z);
    J2D = kron(Jy,Jx);
    uel = reshape(J2D*uvec, Nf,Nf);
    vel = reshape(J2D*vvec, Nf,Nf);
    wel = reshape(J2D*wvec, Nf,Nf);
    Uel = reshape(J2D*U_base, Nf,Nf); % (Nf,Nf)
    uel = repmat(uel, 1, 1, Nf);
    vel = repmat(vel, 1, 1, Nf);
    wel = repmat(wel, 1, 1, Nf);
    Uel = repmat(Uel, 1, 1, Nf); % (Nf,Nf,Nf)
    uel = uel(:);
    vel = vel(:);
    wel = wel(:);
    Uel = Uel(:); % (Nf^3, 1)



    [X,Y,Z] = ndgrid(zfx,zfy,zfz);

    data_r(el,:,1) = X(:);
    data_r(el,:,2) = Y(:);
    data_r(el,:,3) = Z(:);
    data_i(el,:,1) = X(:);
    data_i(el,:,2) = Y(:);
    data_i(el,:,3) = Z(:);
    data_base(el,:,1) = X(:);
    data_base(el,:,2) = Y(:);
    data_base(el,:,3) = Z(:);
    data_r(el,:,4) = real(uel);
    data_r(el,:,5) = real(vel);
    data_r(el,:,6) = real(wel);
    data_i(el,:,4) = imag(uel);
    data_i(el,:,5) = imag(vel);
    data_i(el,:,6) = imag(wel);
    data_base(el,:,6) = Uel;
end
end






E = nelxy*nelxy*nelz;
lr1 = [Nf,Nf,Nf];
elmap = (1:E);
istep = 0;
time= 0;
fields = ['XU'];
emode = 'le';
wdsz = 8; % precision
etag = 6.54321;

status = writenek('u_real0.f00001',
                  data_r,
                  lr1,
                  elmap,
                  time,
                  istep,
                  fields,
                  emode,
                  wdsz,
                  etag);
status = writenek('u_imag0.f00001',
                  data_i,
                  lr1,
                  elmap,
                  time,
                  istep,
                  fields,
                  emode,
                  wdsz,
                  etag);
status = writenek('u_base0.f00001',
                  data_base,
                  lr1,
                  elmap,
                  time,
                  istep,
                  fields,
                  emode,
                  wdsz,
                  etag);
