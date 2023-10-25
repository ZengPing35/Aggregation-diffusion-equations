t_start=0;
t_end=1;
time_interval = 0.02;
free_energy_x = 0:time_interval:t_end;
load('free_energy_32.mat');
Free_energy = free_energy;
plot(free_energy_x, Free_energy, 'r');
xlabel('time');
ylabel('free energy');
title('Free energy');
savefig('Free_energy_32.fig')