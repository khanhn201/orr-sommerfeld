gr = 7.5891076682694115e-03;

cases = {
  'R3',  'r3/r3.log';
  % 'R2',  'r2/r2.log';
  'H06', 'h06/h06.log';
  % 'H06R1', 'h06r1/h06r1.log';
  % 'H06N10', 'h06n10/h06n10.log';
  % 'H06N10R3', 'h06n10r3/h06n10r3.log';
  % 'H06N10I', 'h06n10i/h06n10i.log';
  'H05', 'h05/h05.log';
  % 'H05E2', 'h05e2/h05e2.log';
  % 'H05A2', 'h05a2/h05a2.log';
  % 'H05A2F', 'h05a2f/h05a2f.log';
  % 'H05A4F', 'h05a4f/h05a4f.log';
  % 'H05A4R1F', 'h05a4r1f/h05a4r1f.log';
  % 'H05A4R2F', 'h05a4r2f/h05a4r2f.log';
  % 'H05A4R2C5F', 'h05a4r2c5f/h05a4r2c5f.log';
  % 'H05A4R3C2F', 'h05a4r3c2f/h05a4r3c2f.log';
  % 'M2X10', 'm2x10/m2x10.log';
  'M2X10F', 'm2x10f/m2x10f.log';
  % 'M8X10', 'm8x10/m8x10.log';
  'M8X05', 'm8x05/m8x05.log';
};


function data = read_log(fname)
  fid = fopen(fname,'r');
  C = textscan(fid, '%d %f %f %f %f %f %s');
  fclose(fid);

  data.t = C{2}; % Time
  data.g = C{5}; % Step-wise growth rate
  data.e = C{6}; % Culmulative relative error in growth rate
end

nc = size(cases,1);
D  = struct([]);

for k = 1:nc
  D(k).name = cases{k,1};
  D(k).data = read_log(cases{k,2});
end

% Plot
figpos = [100 100 1200 700];

% figure('Position', figpos);
% hold on;
% for k = 1:nc
%   t = D(k).data.t;
%   g = D(k).data.g;
%   semilogy(t, min(abs(g - gr)/gr, 1.0), 'linewidth', 2);
% end
%
% set(gca, 'fontsize', 20, 'linewidth', 2);
% legend({D.name});
% xlabel('t');
% ylabel('Error');
% title('Step-wise error')


figure('Position', figpos);
hold on;
for k = 1:nc
  t = D(k).data.t;
  e = D(k).data.e;
  semilogy(t, max(min(e, 1.0),1e-7), 'linewidth', 2);
end

set(gca, 'fontsize', 20, 'linewidth', 2);
legend({D.name});
xlim([0 15]);
xlabel('t');
ylabel('Error');
title('Culmulative error')
