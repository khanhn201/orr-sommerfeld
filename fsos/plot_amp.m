gr = 7.5891076682694115e-03;

cases = {
  'R3',  'r3/r3.log';
  'R2',  'r2/r2.log';
  'H06', 'h06/h06.log';
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

figure;
hold on;
for k = 1:nc
  t = D(k).data.t;
  g = D(k).data.g;
  plot(t, abs(g - gr)/gr, 'linewidth', 2);
end

set(gca, 'fontsize', 20, 'linewidth', 2);
legend({D.name}, 'location', 'best');
xlabel('t');
ylabel('Error');
title('Step-wise error')


figure;
hold on;
for k = 1:nc
  t = D(k).data.t;
  e = D(k).data.e;
  plot(t, e, 'linewidth', 2);
end

set(gca, 'fontsize', 20, 'linewidth', 2);
legend({D.name}, 'location', 'best');
xlabel('t');
ylabel('Error');
title('Culmulative error')
