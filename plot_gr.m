fid = fopen('grF.log','r');
C = textscan(fid, '%f %f %f %f %f %f %s');
fclose(fid);

t   = C{2};
amp = C{6};

figure;
semilogy(t, amp, '-', 'linewidth', 2);
xlabel('Time');
ylabel('Relative error');
grid on;

