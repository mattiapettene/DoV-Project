load('kus_th.mat');
x_vec = linspace(0, 1, 1500000);
plot(x_vec, kus_th);
xlim([0 1]);
ylim([-0.006 0]);